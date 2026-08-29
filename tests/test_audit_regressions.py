"""Regression protection for issues found during the post-v1.2.1 forensic audit."""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from metrics import efficacy, summarize, suppression_share_series
from model import IDX
from scenarios import scenario_by_name
from simulate import resolve_params, simulate


def _summary(out, p):
    return summarize(out["t"], out["y"], p, p["Rs5_0"], p["Rs4_0"])


def test_mixed_suppression_share_masks_undefined_early_times():
    loss5 = np.array([0.0, 1e-12, 0.20, 0.30])
    loss4 = np.array([0.0, 1e-12, 0.20, 0.10])
    share = suppression_share_series(loss5, loss4, min_total=1e-6)
    assert np.isnan(share[0]) and np.isnan(share[1])
    assert share[2] == 0.5
    assert np.isclose(share[3], 0.75)


def test_run6_is_secondary_x4_emergence_not_imposed_switch():
    p = resolve_params()
    sc = scenario_by_name("run6_secondary_x4_emergence", p)
    assert sc.V5_0 > 0 and sc.V4_0 == 0
    assert sc.secondary_x4_day == 10.0
    out = simulate(sc, p)
    t = out["t"]
    y = out["y"]
    i10 = int(np.argmin(np.abs(t - 10.0)))
    assert y[i10, IDX["V5"]] > 1e4 * max(y[i10, IDX["V4"]], 1.0)
    after = t >= 10.0
    assert np.any(after & (y[:, IDX["C4"]] > y[:, IDX["C5"]]))


def test_tauB_changes_false_suppression_onset_not_late_extent():
    onsets, late = [], []
    for tau in (0.0, 3.0):
        p = resolve_params({"tau_B": tau})
        out = simulate(scenario_by_name("run3_false_route", p), p)
        t = out["t"]
        loss = 1.0 - out["y"][:, IDX["Rs5"]] / p["Rs5_0"]
        hit = np.where(loss >= 0.80)[0]
        onsets.append(float(t[hit[0]]))
        late.append(float(np.mean(loss[t > 30.0])))
    assert onsets[1] - onsets[0] > 1.0
    assert abs(late[1] - late[0]) < 1e-3
    assert min(late) > 0.99


def test_run9_true_arm_isolation_scales():
    p = resolve_params()
    fast = scenario_by_name("run9_fast_only", p)
    sustained = scenario_by_name("run9_sustained_only", p)
    both = scenario_by_name("run9_both", p)
    assert fast.k_intC_scale == 1.0
    assert fast.k_synC_scale == 0.0 and fast.k_recC_scale == 0.0
    assert sustained.k_intC_scale == 0.0
    assert sustained.k_synC_scale == 1.0 and sustained.k_recC_scale == 1.0
    assert both.k_intC_scale == both.k_synC_scale == both.k_recC_scale == 1.0


def test_run3_false_route_failure_robust_across_representative_sweep_corners():
    vals = []
    effs = []
    for th, TA in ((5/1440, 0.1), (5/1440, 0.7), (2.0, 0.1), (2.0, 0.7)):
        p = resolve_params({"t_half_A": th, "T_A": TA})
        fout = simulate(scenario_by_name("run3_false_route", p), p)
        cout = simulate(scenario_by_name("run2_x4_control", p), p)
        fs, cs = _summary(fout, p), _summary(cout, p)
        vals.append(fs["mean_R5_loss"])
        effs.append(efficacy(fs["auc_V"], cs["auc_V"]))
    assert min(vals) > 0.79
    assert max(vals) - min(vals) < 0.01
    assert max(effs) - min(effs) > 0.2
