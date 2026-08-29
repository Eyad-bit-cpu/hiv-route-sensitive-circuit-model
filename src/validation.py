"""Internal sanity checks (Level 1) plus component checks (Level 2)."""

from __future__ import annotations

import numpy as np

from model import IDX, Scenario, default_params, rhs
from simulate import resolve_params, simulate
from scenarios import const, ligand_pulse


def assert_nonnegative(y: np.ndarray, names) -> None:
    if np.any(y < -1e-8):
        bad = np.where(y.min(axis=0) < -1e-8)[0]
        raise AssertionError(f"negative states: {[names[i] for i in bad]}")


def backbone_r0(params: dict) -> float:
    from model import receptor_infectability

    F5 = receptor_infectability(params["Rs5_0"], params["K_R5"], params["h_R5"])
    return params["beta"] * F5 * params["T0"] * params["p"] / (params["c"] * params["delta_I"])
