"""Model v1.2.1: within-host HIV + temporally gated coreceptor circuit.

Physical outputs: cells/mL, virions/mL, receptor molecules/cell, time in days.
Synthetic signals A, B, C are dimensionless effective signaling variables; they are not
constrained a priori to the interval [0, 1]. Their downstream actions are bounded where
appropriate through Hill-response or saturating regulatory functions.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, Mapping, Optional, Sequence

import numpy as np

MODEL_VERSION = "1.2.1"

STATE = [
    "T_E",
    "T_U",
    "E5",
    "E4",
    "I5",
    "I4",
    "V5",
    "V4",
    "Rs5",
    "Ri5",
    "Rs4",
    "Ri4",
    "A5",
    "A4",
    "B",
    "C5",
    "C4",
    "Q1",
    "Q2",
    "Q3",
]

IDX = {name: i for i, name in enumerate(STATE)}


def hill(x: np.ndarray | float, k: float, n: float) -> np.ndarray | float:
    x = np.maximum(x, 0.0)
    kn = k**n
    return (x**n) / (kn + x**n + 1e-30)


def receptor_infectability(R: float, K: float, h: float) -> float:
    return float(hill(R, K, h))


def k_syn0_from_ss(Rs0: float, k_int0: float, k_rec0: float, k_deg: float) -> float:
    """Synthesis that yields surface density Rs0 at C=0, L=0."""
    return Rs0 * k_int0 * k_deg / max(k_rec0 + k_deg, 1e-12)


def ri0_from_ss(Rs0: float, k_int0: float, k_rec0: float, k_deg: float) -> float:
    """Internal pool consistent with the two-compartment unstimulated steady state."""
    return Rs0 * k_int0 / max(k_rec0 + k_deg, 1e-12)


@dataclass
class Scenario:
    name: str
    V5_0: float = 0.0
    V4_0: float = 0.0
    P5: Callable[[float], float] = field(default_factory=lambda: (lambda t: 0.02))
    P4: Callable[[float], float] = field(default_factory=lambda: (lambda t: 0.02))
    L5: Callable[[float], float] = field(default_factory=lambda: (lambda t: 0.0))
    L4: Callable[[float], float] = field(default_factory=lambda: (lambda t: 0.0))
    t_end: float = 40.0
    circuit_on: bool = True
    constitutive_R5: Optional[float] = None
    hiv_triggered_both: bool = False
    k_intC_scale: float = 1.0
    k_synC_scale: float = 1.0
    k_recC_scale: float = 1.0  # independently scales the recycling-suppression component,
                                 # allowing true fast-only and sustained-only Run-9 ablations.
    f_E: Optional[float] = None
    # Canonical terminology: secondary X4 emergence into an established R5 infection.
    secondary_x4_day: Optional[float] = None
    secondary_x4_V4: float = 10.0
    # Deprecated compatibility aliases retained for old manifests/scripts.
    tropism_shift_day: Optional[float] = None
    tropism_shift_V4: float = 10.0
    # Optional two-phase protocol for recovery experiments. At the registered
    # switch time, exogenous route forcing can already turn off via P5/P4, and
    # the simulator can remove the active infection source while preserving
    # receptor and circuit-memory states so their relaxation is measured cleanly.
    trigger_removal_day: Optional[float] = None
    clear_infection_on_trigger_removal: bool = False
    extra: Dict = field(default_factory=dict)


def default_params(registry_row: Mapping[str, float] | None = None) -> Dict[str, float]:
    """Numeric parameter vector. registry_row overrides baselines."""
    p = {
        "d_T": 0.01,
        "T0": 1.0e6,
        "delta_I": 0.7,
        "p": 4000.0,
        "c": 23.0,
        "k_E": 1.0,
        "R0_target": 8.0,
        "Rs5_0": 10000.0,
        "Rs4_0": 20000.0,
        "K_R5": 1500.0,
        "h_R5": 2.0,
        "K_R4": 3000.0,
        "h_R4": 2.0,
        "k_int0_R5": 2.218,
        "k_rec0_R5": 8.318,
        "k_deg_R5": 2.218,
        "k_int0_R4": 14.4,
        "k_rec0_R4": 11.09,
        "k_deg_R4": 2.0,
        "k_int_lig_R5": 66.36,
        "k_int_lig_R4": 185.2,
        "t_half_A": 0.08333,
        "T_A": 0.35,
        "n_A": 3.0,
        "tau_B": 0.25,
        "t_half_B": 0.5,
        "t_half_C": 0.25,
        "k_f": 20.0,
        "k_r": 5.0,
        "K_B": 1000.0,
        "n_B": 2.0,
        "alpha_B": 8.0,
        "eta": 5.0e-6,
        "k_intC": 80.0,
        "k_recC": 8.0,
        "k_synC": 8.0,
        "f_E": 1.0,
        "V_inoculum": 10.0,
        "P5": 0.02,
        "P4": 0.02,
        "rtol": 1.0e-6,
        "atol": 1.0e-8,
        "master_seed": 20260827.0,
    }
    if registry_row:
        p.update({k: float(v) for k, v in registry_row.items()})
    p["lambda_T"] = p["d_T"] * p["T0"]
    p["k_syn0_R5"] = k_syn0_from_ss(p["Rs5_0"], p["k_int0_R5"], p["k_rec0_R5"], p["k_deg_R5"])
    p["k_syn0_R4"] = k_syn0_from_ss(p["Rs4_0"], p["k_int0_R4"], p["k_rec0_R4"], p["k_deg_R4"])
    F5 = receptor_infectability(p["Rs5_0"], p["K_R5"], p["h_R5"])
    # beta from R0 at baseline surface CCR5 (R5 challenge). Same beta used for X4.
    p["beta"] = p["R0_target"] * p["c"] * p["delta_I"] / (p["T0"] * p["p"] * max(F5, 1e-12))
    p["delta_A"] = math.log(2.0) / p["t_half_A"]
    p["delta_B"] = math.log(2.0) / p["t_half_B"]
    p["delta_C"] = math.log(2.0) / p["t_half_C"]
    return p


def initial_state(params: Mapping[str, float], scenario: Scenario) -> np.ndarray:
    f_E = scenario.f_E if scenario.f_E is not None else params["f_E"]
    y = np.zeros(len(STATE))
    y[IDX["T_E"]] = f_E * params["T0"]
    y[IDX["T_U"]] = (1.0 - f_E) * params["T0"]
    y[IDX["V5"]] = scenario.V5_0
    y[IDX["V4"]] = scenario.V4_0
    y[IDX["Rs5"]] = params["Rs5_0"]
    y[IDX["Ri5"]] = ri0_from_ss(params["Rs5_0"], params["k_int0_R5"], params["k_rec0_R5"], params["k_deg_R5"])
    y[IDX["Rs4"]] = params["Rs4_0"]
    y[IDX["Ri4"]] = ri0_from_ss(params["Rs4_0"], params["k_int0_R4"], params["k_rec0_R4"], params["k_deg_R4"])
    if scenario.constitutive_R5 is not None:
        y[IDX["Rs5"]] = scenario.constitutive_R5 * params["Rs5_0"]
        y[IDX["Ri5"]] = ri0_from_ss(y[IDX["Rs5"]], params["k_int0_R5"], params["k_rec0_R5"], params["k_deg_R5"])
    return y


def rhs(
    t: float,
    y: np.ndarray,
    params: Mapping[str, float],
    scenario: Scenario,
) -> np.ndarray:
    # Adaptive ODE solvers can evaluate intermediate trial states a tiny distance
    # outside the nonnegative orthant. Use nonnegative values for rates here and
    # perform a fail-closed trajectory check after a successful integration.
    y = np.maximum(y, 0.0)
    
    p = params
    s = scenario

    T_E, T_U = y[IDX["T_E"]], y[IDX["T_U"]]
    E5, E4 = y[IDX["E5"]], y[IDX["E4"]]
    I5, I4 = y[IDX["I5"]], y[IDX["I4"]]
    V5, V4 = y[IDX["V5"]], y[IDX["V4"]]
    Rs5, Ri5 = y[IDX["Rs5"]], y[IDX["Ri5"]]
    Rs4, Ri4 = y[IDX["Rs4"]], y[IDX["Ri4"]]
    A5, A4 = y[IDX["A5"]], y[IDX["A4"]]
    B, C5, C4 = y[IDX["B"]], y[IDX["C5"]], y[IDX["C4"]]
    Q1, Q2, Q3 = y[IDX["Q1"]], y[IDX["Q2"]], y[IDX["Q3"]]

    if s.constitutive_R5 is not None:
        # The actual solver state is initialized at the constitutive clamp and
        # both receptor derivatives are held at zero below. Do not mutate y here:
        # solve_ivp passes transient arrays and such mutations are not persistent.
        Rs5 = s.constitutive_R5 * p["Rs5_0"]
        Ri5 = ri0_from_ss(Rs5, p["k_int0_R5"], p["k_rec0_R5"], p["k_deg_R5"])

    F5_E = receptor_infectability(Rs5, p["K_R5"], p["h_R5"])
    F4_E = receptor_infectability(Rs4, p["K_R4"], p["h_R4"])
    F5_U = receptor_infectability(p["Rs5_0"], p["K_R5"], p["h_R5"])
    F4_U = receptor_infectability(p["Rs4_0"], p["K_R4"], p["h_R4"])

    inf5_E = p["beta"] * F5_E * T_E * V5
    inf4_E = p["beta"] * F4_E * T_E * V4
    inf5_U = p["beta"] * F5_U * T_U * V5
    inf4_U = p["beta"] * F4_U * T_U * V4

    f_E = s.f_E if s.f_E is not None else p["f_E"]
    lam_E = f_E * p["lambda_T"]
    lam_U = (1.0 - f_E) * p["lambda_T"]

    dy = np.zeros_like(y)
    dy[IDX["T_E"]] = lam_E - p["d_T"] * T_E - inf5_E - inf4_E
    dy[IDX["T_U"]] = lam_U - p["d_T"] * T_U - inf5_U - inf4_U
    dy[IDX["E5"]] = inf5_E + inf5_U - p["k_E"] * E5
    dy[IDX["E4"]] = inf4_E + inf4_U - p["k_E"] * E4
    dy[IDX["I5"]] = p["k_E"] * E5 - p["delta_I"] * I5
    dy[IDX["I4"]] = p["k_E"] * E4 - p["delta_I"] * I4
    dy[IDX["V5"]] = p["p"] * I5 - p["c"] * V5
    dy[IDX["V4"]] = p["p"] * I4 - p["c"] * V4

    circuit = s.circuit_on
    C5_eff = C5 if circuit else 0.0
    C4_eff = C4 if circuit else 0.0
    if s.hiv_triggered_both:
        # Route-insensitive comparator: B drives both suppressors.
        C5_eff = B
        C4_eff = B

    k_int5 = p["k_int0_R5"] + s.k_intC_scale * p["k_intC"] * C5_eff + p["k_int_lig_R5"] * s.L5(t)
    k_int4 = p["k_int0_R4"] + s.k_intC_scale * p["k_intC"] * C4_eff + p["k_int_lig_R4"] * s.L4(t)
    k_rec5 = p["k_rec0_R5"] / (1.0 + s.k_recC_scale * p["k_recC"] * C5_eff)
    k_rec4 = p["k_rec0_R4"] / (1.0 + s.k_recC_scale * p["k_recC"] * C4_eff)
    k_syn5 = p["k_syn0_R5"] / (1.0 + s.k_synC_scale * p["k_synC"] * C5_eff)
    k_syn4 = p["k_syn0_R4"] / (1.0 + s.k_synC_scale * p["k_synC"] * C4_eff)

    if s.constitutive_R5 is None:
        dy[IDX["Rs5"]] = k_syn5 + k_rec5 * Ri5 - k_int5 * Rs5
        dy[IDX["Ri5"]] = k_int5 * Rs5 - k_rec5 * Ri5 - p["k_deg_R5"] * Ri5
    else:
        # Clamp Rs5 to constitutive value - no dynamics
        dy[IDX["Rs5"]] = 0.0
        dy[IDX["Ri5"]] = 0.0
    dy[IDX["Rs4"]] = k_syn4 + k_rec4 * Ri4 - k_int4 * Rs4
    dy[IDX["Ri4"]] = k_int4 * Rs4 - k_rec4 * Ri4 - p["k_deg_R4"] * Ri4

    P5 = s.P5(t)
    P4 = s.P4(t)
    Eact5 = p["eta"] * inf5_E + P5
    Eact4 = p["eta"] * inf4_E + P4
    dy[IDX["A5"]] = p["delta_A"] * Eact5 - p["delta_A"] * A5
    dy[IDX["A4"]] = p["delta_A"] * Eact4 - p["delta_A"] * A4
    # With this scaling, constant input u produces A* = u (dimensionless).

    I_tot = I5 + I4
    prod_B = p["alpha_B"] * hill(I_tot, p["K_B"], p["n_B"])
    tau = max(float(p["tau_B"]), 0.0)
    if tau <= 1e-9:
        delayed = prod_B
        dy[IDX["Q1"]] = 0.0
        dy[IDX["Q2"]] = 0.0
        dy[IDX["Q3"]] = 0.0
    else:
        kdel = 3.0 / tau
        dy[IDX["Q1"]] = kdel * (prod_B - Q1)
        dy[IDX["Q2"]] = kdel * (Q1 - Q2)
        dy[IDX["Q3"]] = kdel * (Q2 - Q3)
        delayed = Q3

    dy[IDX["B"]] = delayed - p["delta_B"] * B

    HA5 = hill(A5, p["T_A"], p["n_A"])
    HA4 = hill(A4, p["T_A"], p["n_A"])
    HB = hill(B, 0.3, 2.0)
    dy[IDX["C5"]] = p["k_f"] * HA5 * HB - (p["k_r"] + p["delta_C"]) * C5
    dy[IDX["C4"]] = p["k_f"] * HA4 * HB - (p["k_r"] + p["delta_C"]) * C4

    if not circuit:
        dy[IDX["A5"]] = -p["delta_A"] * A5
        dy[IDX["A4"]] = -p["delta_A"] * A4
        dy[IDX["B"]] = -p["delta_B"] * B
        dy[IDX["C5"]] = -(p["k_r"] + p["delta_C"]) * C5
        dy[IDX["C4"]] = -(p["k_r"] + p["delta_C"]) * C4
        dy[IDX["Q1"]] = -Q1
        dy[IDX["Q2"]] = -Q2
        dy[IDX["Q3"]] = -Q3

    return dy
