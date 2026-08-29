"""Statistical helpers shared by campaign and sensitivity scripts."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats


def load_sensitivity_ranges(path: Path, keep: Iterable[str] | None = None) -> dict[str, tuple[float, float, str]]:
    """Load registered Latin-hypercube bounds from CSV (single source of truth)."""
    ranges: dict[str, tuple[float, float, str]] = {}
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            ranges[row["parameter"]] = (float(row["lower"]), float(row["upper"]), row["scale"])
    if keep is not None:
        keep_set = set(keep)
        ranges = {k: v for k, v in ranges.items() if k in keep_set}
    return ranges


def latin_hypercube(n: int, ranges: dict[str, tuple[float, float, str]], seed: int) -> list[dict[str, float]]:
    """Generate a reproducible Latin-hypercube design."""
    if n <= 0:
        raise ValueError("n must be positive")
    rng = np.random.default_rng(seed)
    names = list(ranges)
    dim = len(names)
    u = np.zeros((n, dim))
    for j in range(dim):
        perm = rng.permutation(n)
        u[:, j] = (perm + rng.random(n)) / n

    rows: list[dict[str, float]] = []
    for i in range(n):
        rec: dict[str, float] = {}
        for j, name in enumerate(names):
            lo, hi, scale = ranges[name]
            if scale == "log":
                if lo <= 0 or hi <= 0:
                    raise ValueError(f"log-scale bounds must be positive for {name}")
                rec[name] = float(10 ** (np.log10(lo) + u[i, j] * (np.log10(hi) - np.log10(lo))))
            elif scale == "linear":
                rec[name] = float(lo + u[i, j] * (hi - lo))
            else:
                raise ValueError(f"unknown scale {scale!r} for {name}")
        rows.append(rec)
    return rows


def partial_rank_correlation(
    df: pd.DataFrame,
    param_col: str,
    outcome_col: str,
    parameter_cols: Iterable[str],
) -> tuple[float, float]:
    """Compute PRCC and a t-test p-value.

    All specified parameters and the outcome are rank transformed. The focal
    parameter and outcome are each residualized against the *other sampled
    parameters only*, then their Pearson correlation is computed. The p-value
    uses df = n - k - 2, where k is the number of controlled parameters.
    """
    parameter_cols = list(parameter_cols)
    if param_col not in parameter_cols:
        raise ValueError(f"{param_col!r} is not in parameter_cols")
    required = parameter_cols + [outcome_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(f"missing columns for PRCC: {missing}")

    work = df[required].replace([np.inf, -np.inf], np.nan).dropna()
    n = len(work)
    controls = [c for c in parameter_cols if c != param_col]
    k = len(controls)
    if n <= k + 2:
        return float("nan"), float("nan")

    ranked = work.rank(method="average")
    x = ranked[param_col].to_numpy(dtype=float)
    y = ranked[outcome_col].to_numpy(dtype=float)

    if controls:
        z = ranked[controls].to_numpy(dtype=float)
        z = np.column_stack([np.ones(n), z])
        bx, *_ = np.linalg.lstsq(z, x, rcond=None)
        by, *_ = np.linalg.lstsq(z, y, rcond=None)
        rx = x - z @ bx
        ry = y - z @ by
    else:
        rx = x - np.mean(x)
        ry = y - np.mean(y)

    sx = float(np.linalg.norm(rx))
    sy = float(np.linalg.norm(ry))
    if sx <= 0.0 or sy <= 0.0:
        return float("nan"), float("nan")

    r = float(np.dot(rx, ry) / (sx * sy))
    r = float(np.clip(r, -1.0, 1.0))
    dof = n - k - 2
    if dof <= 0:
        return r, float("nan")
    if abs(r) >= 1.0:
        return r, 0.0
    t_stat = r * np.sqrt(dof / max(1.0 - r * r, np.finfo(float).tiny))
    p_value = float(2.0 * stats.t.sf(abs(t_stat), df=dof))
    return r, p_value
