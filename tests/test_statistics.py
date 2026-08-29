"""Regression tests for sensitivity statistics."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from stats_utils import partial_rank_correlation


def test_prcc_controls_only_registered_parameters_and_recovers_signal():
    rng = np.random.default_rng(1234)
    n = 500
    z = rng.normal(size=n)
    x = rng.normal(size=n)
    y = 2.0 * x + 3.0 * z + rng.normal(scale=0.25, size=n)
    nuisance_output = y + rng.normal(scale=0.01, size=n)
    df = pd.DataFrame({"x": x, "z": z, "y": y, "nuisance_output": nuisance_output})

    rho, p = partial_rank_correlation(df, "x", "y", ["x", "z"])
    assert rho > 0.9
    assert p < 1e-20


def test_prcc_removes_confounding():
    rng = np.random.default_rng(5678)
    n = 600
    z = rng.normal(size=n)
    x = z + rng.normal(scale=0.4, size=n)
    y = z + rng.normal(scale=0.4, size=n)
    df = pd.DataFrame({"x": x, "z": z, "y": y})

    rho, _ = partial_rank_correlation(df, "x", "y", ["x", "z"])
    assert abs(rho) < 0.15
