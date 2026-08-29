"""pytest tests — fail closed: no science until these pass."""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from model import IDX, STATE, Scenario, default_params, hill, receptor_infectability, rhs
from scenarios import const, ligand_pulse
from simulate import resolve_params, simulate
from validation import assert_nonnegative, backbone_r0


def test_r0_matches_ribeiro_median_mapping():
    p = resolve_params()
    r0 = backbone_r0(p)
    assert abs(r0 - 8.0) < 1e-6


def test_positivity_run0a():
    p = resolve_params()
    sc = Scenario("t", V5_0=10.0, circuit_on=False, P5=const(0.02), P4=const(0.02), t_end=20.0)
    out = simulate(sc, p)
    assert_nonnegative(out["y"], STATE)


def test_no_virus_no_infection():
    p = resolve_params()
    sc = Scenario("t", V5_0=0.0, V4_0=0.0, circuit_on=False, t_end=30.0, P5=const(0.0), P4=const(0.0))
    out = simulate(sc, p)
    I = out["y"][:, IDX["I5"]] + out["y"][:, IDX["I4"]]
    assert I[-1] < 1e-6


def test_disable_circuit_keeps_A_B_C_off():
    p = resolve_params()
    sc = Scenario("t", V5_0=10.0, circuit_on=False, t_end=15.0)
    out = simulate(sc, p)
    assert out["y"][:, IDX["C5"]].max() < 1e-8
    assert out["y"][:, IDX["A5"]].max() < 1e-8


def test_and_gate_rhs():
    p = default_params()
    y = np.zeros(len(STATE))
    y[IDX["T_E"]] = p["T0"]
    y[IDX["Rs5"]] = p["Rs5_0"]
    y[IDX["Rs4"]] = p["Rs4_0"]
    sc = Scenario("gate", circuit_on=True, P5=const(0.0), P4=const(0.0))

    y[IDX["A5"]] = 0.0
    y[IDX["B"]] = 0.8
    dy = rhs(0.0, y, p, sc)
    assert dy[IDX["C5"]] <= 1e-9

    y[IDX["A5"]] = 0.9
    y[IDX["B"]] = 0.0
    dy = rhs(0.0, y, p, sc)
    assert dy[IDX["C5"]] <= 1e-9

    y[IDX["A5"]] = 0.9
    y[IDX["B"]] = 0.8
    dy = rhs(0.0, y, p, sc)
    assert dy[IDX["C5"]] > 0.1


def test_A_decay_and_accumulation():
    p = resolve_params()
    sc = Scenario(
        "A",
        circuit_on=True,
        t_end=1.5,
        P5=lambda t: 1.0 if 0.05 <= t <= 0.08 else (1.0 if 0.10 <= t <= 0.13 else 0.0),
        P4=const(0.0),
    )
    out = simulate(sc, p, n_out=1500)
    A = out["y"][:, IDX["A5"]]
    t = out["t"]
    # after first pulse A rises then falls; second nearby pulse should reach higher than first peak neighborhood
    i1 = np.searchsorted(t, 0.09)
    i2 = np.searchsorted(t, 0.14)
    assert A[i1] > 0.05
    assert A[i2] >= A[i1] * 0.8


def test_cxcr4_fast_internalization_order_of_magnitude():
    """Signoret 1997: 50% surface CXCR4 in ~5 min under SDF-1."""
    p = resolve_params()
    sc = Scenario("c", circuit_on=False, t_end=0.02, L4=ligand_pulse(0.0, 0.02), P5=const(0), P4=const(0))
    out = simulate(sc, p, n_out=400)
    t = out["t"]
    Rs = out["y"][:, IDX["Rs4"]]
    i5 = np.searchsorted(t, 5.0 / 1440.0)
    frac = Rs[i5] / Rs[0]
    assert 0.25 < frac < 0.75


def test_ccr5_slower_than_cxcr4_ligand():
    p = resolve_params()
    sc5 = Scenario("5", circuit_on=False, t_end=0.02, L5=ligand_pulse(0.0, 0.02), P5=const(0), P4=const(0))
    sc4 = Scenario("4", circuit_on=False, t_end=0.02, L4=ligand_pulse(0.0, 0.02), P5=const(0), P4=const(0))
    r5 = simulate(sc5, p, n_out=400)
    r4 = simulate(sc4, p, n_out=400)
    i = np.searchsorted(r5["t"], 5.0 / 1440.0)
    loss5 = 1.0 - r5["y"][i, IDX["Rs5"]] / r5["y"][0, IDX["Rs5"]]
    loss4 = 1.0 - r4["y"][i, IDX["Rs4"]] / r4["y"][0, IDX["Rs4"]]
    assert loss4 > loss5


def test_hill_saturation_thirty_percent_cut():
    """Platt high-CD4 mapping: 10000 → 7000 can remain on the flat part of F."""
    p = resolve_params()
    f_hi = receptor_infectability(10000.0, p["K_R5"], p["h_R5"])
    f_cut = receptor_infectability(7000.0, p["K_R5"], p["h_R5"])
    assert f_hi > 0.95
    assert (f_hi - f_cut) / f_hi < 0.1


def test_no_B_implies_C_stays_small_without_I():
    p = resolve_params()
    sc = Scenario("n", V5_0=0.0, circuit_on=True, P5=const(1.0), P4=const(0.0), t_end=3.0)
    out = simulate(sc, p)
    assert out["y"][:, IDX["B"]].max() < 1e-6
    assert out["y"][:, IDX["C5"]].max() < 0.05


def test_solver_lsoda_vs_bdf():
    p = resolve_params()
    sc = Scenario("s", V5_0=10.0, circuit_on=False, t_end=10.0)
    a = simulate(sc, p, method="LSODA")
    b = simulate(sc, p, method="BDF")
    v1 = a["y"][:, IDX["V5"]]
    v2 = np.interp(a["t"], b["t"], b["y"][:, IDX["V5"]])
    rel = np.max(np.abs(v1 - v2) / (np.maximum(v1, 1.0)))
    assert rel < 0.05


def test_receptor_baseline_is_actual_steady_state():
    p = resolve_params()
    sc = Scenario("steady", V5_0=0.0, V4_0=0.0, circuit_on=False, P5=const(0.0), P4=const(0.0))
    from model import initial_state
    y0 = initial_state(p, sc)
    dy = rhs(0.0, y0, p, sc)
    for name in ["Rs5", "Ri5", "Rs4", "Ri4"]:
        assert abs(dy[IDX[name]]) < 1e-8

    out = simulate(sc, p, n_out=100)
    assert np.allclose(out["y"][:, IDX["Rs5"]], p["Rs5_0"], rtol=0, atol=1e-6)
    assert np.allclose(out["y"][:, IDX["Rs4"]], p["Rs4_0"], rtol=0, atol=1e-6)


def test_constitutive_r5_comparator_is_clamped_in_stored_state():
    from scenarios import scenario_by_name
    p = resolve_params()
    out = simulate(scenario_by_name("comparator_ccr5_ko", p), p, n_out=100)
    ratio = out["y"][:, IDX["Rs5"]] / p["Rs5_0"]
    assert np.allclose(ratio, 0.05, rtol=0, atol=1e-10)


def test_recovery_targets_suppressed_receptor_not_either_receptor():
    from metrics import recovery_time
    t = np.array([0.0, 1.0, 2.0, 3.0])
    y = np.zeros((len(t), len(STATE)))
    y[:, IDX["Rs5"]] = np.array([100.0, 20.0, 40.0, 95.0])
    y[:, IDX["Rs4"]] = 100.0
    assert recovery_time(t, y, 100.0, 100.0, t_off=1.0, receptor="R5") == pytest.approx(1.9090909)
    assert recovery_time(t, y, 100.0, 100.0, t_off=1.0, receptor="R4") == pytest.approx(0.0)


def test_unknown_parameter_override_fails_closed():
    with pytest.raises(KeyError):
        resolve_params({"typo_parameter": 1.0})


def test_recovery_protocol_removes_trigger_source_but_preserves_memory_states():
    """Recovery release removes active infection, not the stored circuit state."""
    from scenarios import scenario_by_name
    p = resolve_params()
    sc = scenario_by_name("run_recovery_test", p)
    out = simulate(sc, p, n_out=2000)
    t = out["t"]
    y = out["y"]
    i_after = np.searchsorted(t, 14.05)
    # Infection source is removed at the release and cannot regrow without virus.
    for state in ["E5", "E4", "I5", "I4", "V5", "V4"]:
        assert np.max(y[i_after:, IDX[state]]) < 1e-6
    # Stored confirmation/circuit states are not zeroed by hand; they relax naturally.
    i_rel = np.searchsorted(t, 14.0)
    assert y[i_rel, IDX["B"]] > 0.0 or y[i_rel, IDX["Q3"]] > 0.0 or y[i_rel, IDX["C5"]] > 0.0
    assert y[-1, IDX["B"]] < y[i_rel, IDX["B"]] + 1e-9


def test_run3_registered_forcing_is_extreme_relative_to_internal_route_input_and_sweep_is_graded():
    """Reviewer hardening: Run 3 stress forcing is internally scaled and false loss rises with forcing."""
    from dataclasses import replace
    from scenarios import scenario_by_name, const
    from receptors import receptor_infectability
    from metrics import summarize
    p = resolve_params()
    r5 = simulate(scenario_by_name("run1_r5_circuit", p), p, n_out=800)
    y = r5["y"]
    F5 = np.array([receptor_infectability(r, p["K_R5"], p["h_R5"]) for r in y[:, IDX["Rs5"]]])
    drive = p["eta"] * p["beta"] * F5 * y[:, IDX["T_E"]] * y[:, IDX["V5"]]
    ref = float(np.max(drive))
    assert 1.2 / ref > 5.0

    base = scenario_by_name("run3_false_route", p)
    losses = []
    for strength in [0.02, 0.5 * ref, ref, 1.2]:
        out = simulate(replace(base, P5=const(float(strength))), p, n_out=800)
        sm = summarize(out["t"], out["y"], p, p["Rs5_0"], p["Rs4_0"])
        losses.append(sm["mean_R5_loss"])
    assert losses[0] < 0.05
    assert losses[0] < losses[1] < losses[2] < losses[3]
