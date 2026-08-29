"""Integrate Model v1.2.1 and write run manifests."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np
from scipy.integrate import solve_ivp

from model import MODEL_VERSION, STATE, Scenario, default_params, initial_state, rhs

ROOT = Path(__file__).resolve().parents[1]


def git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "NO_GIT"


def file_hash(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()[:16]


def child_seed(master: int, run_id: str) -> int:
    digest = hashlib.sha256(f"{master}:{run_id}".encode()).digest()
    return int.from_bytes(digest[:4], "little")


def load_registry_baselines(csv_path: Optional[Path] = None) -> Dict[str, float]:
    import csv

    path = csv_path or ROOT / "parameters" / "parameter_registry.csv"
    out: Dict[str, float] = {}
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            try:
                out[row["parameter"]] = float(row["baseline"])
            except (ValueError, KeyError):
                continue
    return out


def resolve_params(overrides: Optional[Dict[str, float]] = None) -> Dict[str, float]:
    base = load_registry_baselines()
    p = default_params(base)
    if overrides:
        unknown = sorted(set(overrides) - set(p))
        if unknown:
            raise KeyError(f"Unknown parameter override(s): {unknown}")
        p.update({k: float(v) for k, v in overrides.items()})
        p["lambda_T"] = p["d_T"] * p["T0"]
        from model import k_syn0_from_ss, receptor_infectability

        p["k_syn0_R5"] = k_syn0_from_ss(p["Rs5_0"], p["k_int0_R5"], p["k_rec0_R5"], p["k_deg_R5"])
        p["k_syn0_R4"] = k_syn0_from_ss(p["Rs4_0"], p["k_int0_R4"], p["k_rec0_R4"], p["k_deg_R4"])
        F5 = receptor_infectability(p["Rs5_0"], p["K_R5"], p["h_R5"])
        p["beta"] = p["R0_target"] * p["c"] * p["delta_I"] / (p["T0"] * p["p"] * max(F5, 1e-12))
        import math

        p["delta_A"] = math.log(2.0) / p["t_half_A"]
        p["delta_B"] = math.log(2.0) / p["t_half_B"]
        p["delta_C"] = math.log(2.0) / p["t_half_C"]
    return p


def _integrate_segment(
    t0: float,
    t1: float,
    y0: np.ndarray,
    params: Dict[str, float],
    scenario: Scenario,
    method: str,
) -> Tuple[np.ndarray, np.ndarray]:
    sol = solve_ivp(
        lambda t, y: rhs(t, y, params, scenario),
        (t0, t1),
        y0,
        method=method,
        rtol=params["rtol"],
        atol=params["atol"],
        dense_output=False,
        max_step=0.05,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    y_out = sol.y.T
    # Fail closed on meaningful negative trajectory states, while tolerating and
    # clipping tiny solver-scale excursions that can occur in stiff delay states.
    neg_tol = max(float(params["atol"]) * 100.0, 1e-8)
    min_val = float(np.min(y_out))
    if min_val < -neg_tol:
        i_t, i_s = np.unravel_index(int(np.argmin(y_out)), y_out.shape)
        raise ValueError(
            f"Negative state after integration at t={sol.t[i_t]:.6g}: "
            f"{STATE[i_s]}={y_out[i_t, i_s]:.6g} (tolerance {neg_tol:.3g})"
        )
    y_out = np.maximum(y_out, 0.0)
    return sol.t, y_out


def simulate(
    scenario: Scenario,
    params: Optional[Dict[str, float]] = None,
    method: str = "LSODA",
    n_out: int = 800,
) -> Dict[str, Any]:
    params = params or resolve_params()
    y0 = initial_state(params, scenario)
    t_end = scenario.t_end
    t_shift = scenario.secondary_x4_day if scenario.secondary_x4_day is not None else scenario.tropism_shift_day
    t_release = scenario.trigger_removal_day

    if t_release is not None and t_release < t_end:
        if t_shift is not None and t_shift < t_end:
            raise ValueError("trigger_removal_day and secondary_x4_day are not supported in the same scenario")
        t1, y1 = _integrate_segment(0.0, t_release, y0, params, scenario, method)
        y_mid = y1[-1].copy()
        from model import IDX

        if scenario.clear_infection_on_trigger_removal:
            # Remove the upstream infection source at the registered release
            # event, but preserve A/B/C, Q1-Q3, and receptor states. Any delayed
            # confirmation already present must therefore decay naturally.
            for state in ("E5", "E4", "I5", "I4", "V5", "V4"):
                y_mid[IDX[state]] = 0.0

        t2, y2 = _integrate_segment(t_release, t_end, y_mid, params, scenario, method)
        # The release is a deliberate state discontinuity. Retain the post-release
        # state at t_release so interpolation does not bridge across the reset.
        t = np.concatenate([t1[:-1], t2])
        y = np.concatenate([y1[:-1], y2], axis=0)
    elif t_shift is None or t_shift >= t_end:
        t, y = _integrate_segment(0.0, t_end, y0, params, scenario, method)
    else:
        t1, y1 = _integrate_segment(0.0, t_shift, y0, params, scenario, method)
        y_mid = y1[-1].copy()
        from model import IDX

        v4_seed = scenario.secondary_x4_V4 if scenario.secondary_x4_day is not None else scenario.tropism_shift_V4
        y_mid[IDX["V4"]] += v4_seed
        t2, y2 = _integrate_segment(t_shift, t_end, y_mid, params, scenario, method)
        t = np.concatenate([t1, t2[1:]])
        y = np.concatenate([y1, y2[1:]], axis=0)

    t_grid = np.linspace(0.0, t_end, n_out)
    y_grid = np.vstack([np.interp(t_grid, t, y[:, i]) for i in range(y.shape[1])]).T
    return {"t": t_grid, "y": y_grid, "t_native": t, "y_native": y, "params": params, "method": method}


def save_run(
    run_id: str,
    scenario: Scenario,
    result: Dict[str, Any],
    extras: Optional[Dict[str, Any]] = None,
) -> Path:
    raw = ROOT / "results" / "raw" / run_id
    raw.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(raw / "trajectory.npz", t=result["t"], y=result["y"], state=np.array(STATE))
    p = result["params"]
    manifest = {
        "run_id": run_id,
        "git_commit": git_commit(),
        "random_seed": child_seed(int(p["master_seed"]), run_id),
        "model_version": MODEL_VERSION,
        "solver": result["method"],
        "rtol": p["rtol"],
        "atol": p["atol"],
        "simulation_start": 0.0,
        "simulation_end": scenario.t_end,
        "parameter_file_hash": file_hash(ROOT / "parameters" / "parameter_registry.csv"),
        "scenario": scenario.name,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "circuit_on": scenario.circuit_on,
        "parameters": {k: float(v) if isinstance(v, (int, float, np.floating)) else v for k, v in p.items()},
    }
    if extras:
        manifest.update(extras)
    (raw / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return raw
