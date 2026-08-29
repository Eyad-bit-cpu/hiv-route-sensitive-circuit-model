#!/usr/bin/env python
"""Deterministic campaign: unit modules → challenges → sweeps → figures.

Does not retune Class A biology. Engineered parameters are swept, not secretly fitted.
"""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from model import IDX, STATE, Scenario, receptor_infectability
from metrics import classify, efficacy, summarize, recovery_time, suppression_share_series
from scenarios import const, scenario_by_name
from simulate import resolve_params, save_run, simulate
from stats_utils import latin_hypercube, load_sensitivity_ranges, partial_rank_correlation

FIG = ROOT / "results" / "figures"
SUM = ROOT / "results" / "summaries"
FIG.mkdir(parents=True, exist_ok=True)
SUM.mkdir(parents=True, exist_ok=True)


def criteria():
    with (ROOT / "parameters" / "success_criteria_v1.yaml").open(encoding="utf-8") as f:
        return yaml.safe_load(f)


def series(result, name):
    return result["t"], result["y"][:, IDX[name]]


def savefig(fig, name, dpi=160):
    fig.tight_layout()
    fig.savefig(FIG / name, dpi=dpi)
    plt.close(fig)


def run_named(run_id, scenario, params, extras=None):
    out = simulate(scenario, params)
    save_run(run_id, scenario, out, extras)
    sm = summarize(out["t"], out["y"], params, params["Rs5_0"], params["Rs4_0"])
    sm["run_id"] = run_id
    sm["scenario"] = scenario.name
    return out, sm


def plot_states(result, keys, title, fname, logy=False):
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    t = result["t"]
    for k in keys:
        y = result["y"][:, IDX[k]]
        ax.plot(t, y, label=k)
    ax.set_xlabel("time (days)")
    ax.set_title(title)
    if logy:
        ax.set_yscale("log")
        ax.set_ylim(bottom=1e-2)
    ax.legend(fontsize=8, ncol=2)
    savefig(fig, fname)



def main():
    params = resolve_params()
    crit = criteria()
    rows = []

    # --- Phase I/II: Run 0 ---
    jobs = [
        ("run0a_backbone", "run0a_backbone"),
        ("run0b_ccr5_trafficking", "run0b_ccr5_trafficking"),
        ("run0c_cxcr4_trafficking", "run0c_cxcr4_trafficking"),
        ("run0d_A_memory", "run0d_A_memory"),
        ("run0e_and_gate", "run0e_and_gate"),
    ]
    results = {}
    for rid, sname in jobs:
        sc = scenario_by_name(sname, params)
        out, sm = run_named(rid, sc, params)
        results[rid] = out
        rows.append(sm)

    plot_states(results["run0a_backbone"], ["T_E", "E5", "I5", "V5"], "HIV backbone dynamics (circuit OFF)", "fig_run0a_backbone.png", logy=False)
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    t, v = series(results["run0a_backbone"], "V5")
    ax.plot(t, v)
    ax.set_yscale("log")
    ax.set_xlabel("time (days)")
    ax.set_ylabel("V5 (virions/mL)")
    ax.set_title("Free R5 virus dynamics (circuit disabled)")
    savefig(fig, "fig_run0a_V5_log.png")

    # Normalized receptor plots for trafficking
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    t = results["run0b_ccr5_trafficking"]["t"]
    Rs5 = results["run0b_ccr5_trafficking"]["y"][:, IDX["Rs5"]]
    Ri5 = results["run0b_ccr5_trafficking"]["y"][:, IDX["Ri5"]]
    # Calculate baseline internal pool
    Ri5_0 = params["Rs5_0"] * params["k_int0_R5"] / (params["k_rec0_R5"] + params["k_deg_R5"])
    ax.plot(t, Rs5 / params["Rs5_0"], label="R_s5/R_s5,0")
    ax.plot(t, Ri5 / Ri5_0, label="R_i5/R_i5,0")
    ax.set_xlabel("time (days)")
    ax.set_ylabel("Normalized receptor abundance")
    ax.set_title("CCR5 trafficking response to ligand pulse")
    ax.legend()
    savefig(fig, "fig_run0b_ccr5.png")
    
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    t = results["run0c_cxcr4_trafficking"]["t"]
    Rs4 = results["run0c_cxcr4_trafficking"]["y"][:, IDX["Rs4"]]
    Ri4 = results["run0c_cxcr4_trafficking"]["y"][:, IDX["Ri4"]]
    # Calculate baseline internal pool
    Ri4_0 = params["Rs4_0"] * params["k_int0_R4"] / (params["k_rec0_R4"] + params["k_deg_R4"])
    ax.plot(t, Rs4 / params["Rs4_0"], label="R_s4/R_s4,0")
    ax.plot(t, Ri4 / Ri4_0, label="R_i4/R_i4,0")
    ax.set_xlabel("time (days)")
    ax.set_ylabel("Normalized receptor abundance")
    ax.set_title("CXCR4 trafficking response to ligand pulse")
    ax.legend()
    savefig(fig, "fig_run0c_cxcr4.png")
    plot_states(results["run0d_A_memory"], ["A5"], "Route-memory signal integration (A5 dynamics)", "fig_run0d_A.png")

    # infectability curve (literature mapping, not a fit)
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    R = np.linspace(0, 25000, 400)
    F = [receptor_infectability(r, params["K_R5"], params["h_R5"]) for r in R]
    ax.plot(R, F, color="k")
    ax.axvline(10000, ls="--", color="C0", label="baseline 10,000 (Reynes-scale)")
    ax.axvline(1500, ls=":", color="C2", label="K_R5=1500 (model-mapped midpoint)")
    ax.axvline(10000, ymin=0, ymax=0)  # keep
    ax.set_xlabel("surface CCR5 (molecules/cell)")
    ax.set_ylabel(r"$F_5(R_{s5})$")
    ax.set_title("Nonlinear dependence of modeled infectability on CCR5 surface abundance")
    ax.legend(fontsize=8)
    savefig(fig, "fig_F5_hill_platt_mapping.png")

    # --- Run 1–2 ---
    ctrl1, sm = run_named("run1_r5_control", scenario_by_name("run1_r5_control", params), params)
    rows.append(sm)
    cir1, sm = run_named("run1_r5_circuit", scenario_by_name("run1_r5_circuit", params), params)
    rows.append(sm)
    ev1 = efficacy(sm["auc_V"], rows[-2]["auc_V"])
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.0))
    for ax, key in zip(axes.ravel(), ["V5", "Rs5", "A5", "C5"]):
        if key == "Rs5":
            # Normalize receptor
            ax.plot(ctrl1["t"], ctrl1["y"][:, IDX[key]] / params["Rs5_0"], label="control")
            ax.plot(cir1["t"], cir1["y"][:, IDX[key]] / params["Rs5_0"], label="circuit")
            ax.set_ylabel("R_s5/R_s5,0")
        else:
            ax.plot(ctrl1["t"], ctrl1["y"][:, IDX[key]], label="control")
            ax.plot(cir1["t"], cir1["y"][:, IDX[key]], label="circuit")
            ax.set_title(key)
        ax.legend(fontsize=8)
    fig.suptitle("R5 HIV infection: control vs circuit")
    savefig(fig, "fig_run1_r5.png")

    ctrl2, smc = run_named("run2_x4_control", scenario_by_name("run2_x4_control", params), params)
    rows.append(smc)
    cir2, smx = run_named("run2_x4_circuit", scenario_by_name("run2_x4_circuit", params), params)
    rows.append(smx)
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4))
    for ax, key in zip(axes, ["V4", "Rs4", "Rs5"]):
        if key in ["Rs4", "Rs5"]:
            # Normalize receptors
            baseline = params["Rs4_0"] if key == "Rs4" else params["Rs5_0"]
            ax.plot(ctrl2["t"], ctrl2["y"][:, IDX[key]] / baseline, label="control")
            ax.plot(cir2["t"], cir2["y"][:, IDX[key]] / baseline, label="circuit")
            ax.set_ylabel(f"{key}/{key}_0")
        else:
            ax.plot(ctrl2["t"], ctrl2["y"][:, IDX[key]], label="control")
            ax.plot(cir2["t"], cir2["y"][:, IDX[key]], label="circuit")
            ax.set_title(key)
        ax.legend(fontsize=8)
    fig.suptitle("X4 HIV infection: control vs circuit")
    savefig(fig, "fig_run2_x4.png")

    # --- Run 3 heatmap ---
    t_halves = np.geomspace(5 / 1440, 2.0, 8)
    TAs = np.linspace(0.1, 0.7, 8)
    EV = np.zeros((len(TAs), len(t_halves)))
    FF = np.zeros_like(EV)
    for i, TA in enumerate(TAs):
        for j, th in enumerate(t_halves):
            p2 = resolve_params({"t_half_A": float(th), "T_A": float(TA)})
            sc_c = scenario_by_name("run2_x4_control", p2)
            sc_f = scenario_by_name("run3_false_route", p2)
            c_out = simulate(sc_c, p2)
            f_out = simulate(sc_f, p2)
            sc_sum = summarize(c_out["t"], c_out["y"], p2, p2["Rs5_0"], p2["Rs4_0"])
            fs = summarize(f_out["t"], f_out["y"], p2, p2["Rs5_0"], p2["Rs4_0"])
            # Use X4 control for X4 false-route scenario
            EV[i, j] = efficacy(fs["auc_V"], sc_sum["auc_V"])
            FF[i, j] = fs["mean_R5_loss"]  # false CCR5 suppression during X4 + high P5
            rows.append(
                {
                    "run_id": f"run3_{i}_{j}",
                    "scenario": "run3_false_route",
                    "t_half_A": th,
                    "T_A": TA,
                    "E_V": EV[i, j],
                    "false_R5_loss": FF[i, j],
                    **{k: fs[k] for k in ("auc_V", "peak_V", "mean_R4_loss")},
                }
            )
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    extent = [np.log10(t_halves[0]), np.log10(t_halves[-1]), TAs[0], TAs[-1]]
    im0 = axes[0].imshow(EV, origin="lower", aspect="auto", extent=extent, vmin=0, vmax=1, cmap="viridis")
    axes[0].set_title("Antiviral efficacy E_V (X4 + high P5)")
    axes[0].set_xlabel("log10 t_half_A (days)")
    axes[0].set_ylabel("T_A")
    plt.colorbar(im0, ax=axes[0], fraction=0.046)
    im1 = axes[1].imshow(FF, origin="lower", aspect="auto", extent=extent, vmin=0, vmax=1, cmap="magma")
    axes[1].set_title("False CCR5 suppression")
    axes[1].set_xlabel("log10 t_half_A (days)")
    plt.colorbar(im1, ax=axes[1], fraction=0.046)
    savefig(fig, "fig_run3_false_route_heatmap.png")

    # Efficacy variation under robust persistent false-route suppression
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    ax.scatter(FF.ravel(), EV.ravel(), c="k", s=22)
    ax.axvline(crit["false_suppression"]["F_max_success"], ls="--", color="C1")
    ax.axhline(crit["antiviral_efficacy"]["E_V_success"], ls="--", color="C2")
    ax.set_xlabel("false CCR5 suppression")
    ax.set_ylabel("E_V")
    ax.set_title("Efficacy variation under persistent false-route suppression")
    savefig(fig, "fig_run3_robust_false_route.png")

    # Graded wrong-route forcing sweep, scaled to the model's own infection-derived
    # R5 route-input magnitude. This provides an internal model scale for the
    # otherwise dimensionless P5 forcing without claiming physiological calibration.
    r5_ref = simulate(scenario_by_name("run1_r5_circuit", params), params)
    tr = r5_ref["t"]
    yr = r5_ref["y"]
    F5_ref = np.array([receptor_infectability(r, params["K_R5"], params["h_R5"]) for r in yr[:, IDX["Rs5"]]])
    inf5E_ref = params["beta"] * F5_ref * yr[:, IDX["T_E"]] * yr[:, IDX["V5"]]
    route_ref_peak = float(np.max(params["eta"] * inf5E_ref))
    route_ref_tpeak = float(tr[int(np.argmax(params["eta"] * inf5E_ref))])

    forcing_abs = np.array([0.0, 0.02, 0.05, 0.10, route_ref_peak, 0.20, 0.40, 0.80, 1.20], dtype=float)
    forcing_ratio = forcing_abs / max(route_ref_peak, 1e-30)
    forcing_false = []
    forcing_final = []
    for strength in forcing_abs:
        sc = replace(scenario_by_name("run3_false_route", params), P5=const(float(strength)))
        out = simulate(sc, params)
        sm = summarize(out["t"], out["y"], params, params["Rs5_0"], params["Rs4_0"])
        forcing_false.append(float(sm["mean_R5_loss"]))
        forcing_final.append(float(max(1.0 - out["y"][-1, IDX["Rs5"]] / params["Rs5_0"], 0.0)))

    forcing_df = pd.DataFrame({
        "P5": forcing_abs,
        "P5_over_peak_infection_drive": forcing_ratio,
        "mean_false_R5_loss": forcing_false,
        "final_false_R5_loss": forcing_final,
    })
    forcing_df.to_csv(SUM / "run3_wrong_route_forcing_sweep.csv", index=False)
    (SUM / "run3_wrong_route_forcing_reference.json").write_text(
        json.dumps({
            "reference": "peak eta*inf5_E in baseline run1_r5_circuit",
            "peak_infection_derived_R5_route_input": route_ref_peak,
            "time_of_peak_days": route_ref_tpeak,
            "registered_run3_P5": 1.2,
            "registered_run3_P5_over_reference_peak": 1.2 / route_ref_peak,
            "registered_run3_duration_days": float(scenario_by_name("run3_false_route", params).t_end),
            "baseline_background_P5": float(params["P5"]),
            "baseline_background_P5_over_reference_peak": float(params["P5"] / route_ref_peak),
            "interpretation": "Internal model-scale comparison only; P5 is dimensionless and is not calibrated to a measured physiological CCR5 signal.",
        }, indent=2), encoding="utf-8")

    fig, ax = plt.subplots(figsize=(6.2, 4.5))
    ax.plot(forcing_ratio, forcing_false, marker="o", label="40-day mean false R5 loss")
    ax.plot(forcing_ratio, forcing_final, marker="s", ls="--", label="final-time false R5 loss")
    ax.axhline(crit["false_suppression"]["F_max_success"], ls=":", color="C2", label="full-success ceiling")
    ax.axhline(crit["false_suppression"]["F_max_partial"], ls=":", color="C1", label="partial ceiling")
    ax.axvline(params["P5"] / route_ref_peak, ls="--", color="0.45", label="baseline background P5")
    ax.axvline(1.2 / route_ref_peak, ls="--", color="C3", label="registered stress P5=1.2")
    ax.set_xlabel(r"constant wrong-route forcing $P_5$ / peak baseline infection-derived $\eta\,\mathrm{inf}_{5,E}$")
    ax.set_ylabel("False CCR5 suppression")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("Graded persistent wrong-route forcing")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7, loc="lower right")
    savefig(fig, "fig_run3_wrong_route_forcing_sweep.png")

    # --- Run 4: external forcing patterns during ongoing infection ---
    for sname in ["run4_isolated_pulse", "run4_repeated_pulses", "run4_sustained"]:
        scenario4 = scenario_by_name(sname, params)
        out, sm = run_named(sname, scenario4, params)
        rows.append(sm)
        # Three-panel diagnostic exposing the two additive drivers of A5:
        # scripted P5(t) and infection-driven eta*inf5_E, plus downstream B/C5
        # and receptor phenotype. This prevents the external input waveform from
        # being mistaken for the total circuit drive once infection is established.
        fig, axes = plt.subplots(3, 1, figsize=(8.0, 8.4), sharex=True)
        t = out["t"]
        y = out["y"]
        A5 = y[:, IDX["A5"]]
        B = y[:, IDX["B"]]
        C5 = y[:, IDX["C5"]]
        Rs5 = y[:, IDX["Rs5"]]
        T_E = y[:, IDX["T_E"]]
        V5 = y[:, IDX["V5"]]
        F5 = np.array([receptor_infectability(r, params["K_R5"], params["h_R5"]) for r in Rs5])
        inf5_E = params["beta"] * F5 * T_E * V5
        scripted = np.array([scenario4.P5(float(tt)) for tt in t])
        infection_drive = params["eta"] * inf5_E

        axes[0].plot(t, scripted, label=r"scripted $P_5(t)$")
        axes[0].plot(t, infection_drive, label=r"infection-driven $\eta\,\mathrm{inf}_{5,E}$")
        axes[0].set_ylabel("Drive to $A_5$")
        friendly = {
            "run4_isolated_pulse": "isolated scripted pulse",
            "run4_repeated_pulses": "repeated scripted pulses",
            "run4_sustained": "sustained scripted input",
        }[sname]
        axes[0].set_title(f"External forcing + infection feedback: {friendly}")
        axes[0].legend(fontsize=8)
        axes[0].grid(True, alpha=0.25)

        axes[1].plot(t, A5, label=r"$A_5$")
        axes[1].plot(t, C5, label=r"$C_5$")
        axb = axes[1].twinx()
        axb.plot(t, B, color="C3", alpha=0.75, label=r"$B$")
        axes[1].set_ylabel(r"$A_5, C_5$")
        axb.set_ylabel(r"$B$")
        lines = axes[1].get_lines() + axb.get_lines()
        axes[1].legend(lines, [ln.get_label() for ln in lines], fontsize=8, loc="upper right")
        axes[1].grid(True, alpha=0.25)

        axes[2].plot(t, Rs5 / params["Rs5_0"], label=r"$R_{s5}/R_{s5,0}$")
        axes[2].set_xlabel("time (days)")
        axes[2].set_ylabel("Normalized CCR5")
        axes[2].set_ylim(-0.02, 1.05)
        axes[2].legend(fontsize=8)
        axes[2].grid(True, alpha=0.25)

        plt.tight_layout()
        savefig(fig, f"fig_{sname}.png")

    # --- Run 5 mixed ---
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for sname in ["run5_mix_90_10", "run5_mix_50_50", "run5_mix_10_90"]:
        out, sm = run_named(sname, scenario_by_name(sname, params), params)
        rows.append(sm)
        V5, V4 = out["y"][:, IDX["V5"]], out["y"][:, IDX["V4"]]
        loss5 = np.maximum(1.0 - out["y"][:, IDX["Rs5"]] / params["Rs5_0"], 0.0)
        loss4 = np.maximum(1.0 - out["y"][:, IDX["Rs4"]] / params["Rs4_0"], 0.0)
        S = suppression_share_series(loss5, loss4, min_total=1e-6)
        trop = V5 / (V5 + V4 + 1e-30)
        ax.plot(out["t"], S, label=f"S CCR5 {sname[-5:]}")
        ax.plot(out["t"], trop, ls="--", label=f"V5 fraction {sname[-5:]}")
    ax.set_xlabel("time (days)")
    ax.set_ylim(0, 1)
    ax.set_title("Suppression share vs viral tropism mix")
    ax.legend(fontsize=7, ncol=2)
    savefig(fig, "fig_run5_mixed.png")

    # --- Run 6: secondary X4 emergence into established R5 infection ---
    out, sm = run_named("run6_secondary_x4_emergence", scenario_by_name("run6_secondary_x4_emergence", params), params)
    rows.append(sm)
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(out["t"], out["y"][:, IDX["C5"]], label="C5")
    ax.plot(out["t"], out["y"][:, IDX["C4"]], label="C4")
    ax.axvline(10.0, ls="--", color="k", label="X4 introduced")
    ax.legend()
    ax.set_xlabel("time (days)")
    ax.set_title("Secondary X4 emergence: $C_4$ vs $C_5$ dynamics")
    savefig(fig, "fig_run6_secondary_x4_emergence.png")
    c5, c4, t = out["y"][:, IDX["C5"]], out["y"][:, IDX["C4"]], out["t"]
    after = t >= 10
    crossover_days = None
    if np.any(after & (c4 > c5)):
        crossover_days = float(t[np.where(after & (c4 > c5))[0][0]] - 10.0)
    (SUM / "run6_crossover_time.json").write_text(
        json.dumps({
            "secondary_x4_seed_day": 10.0,
            "c4_exceeds_c5_after_days": crossover_days,
            "interpretation": "secondary X4 emergence into established R5 infection; not an imposed route switch",
        }, indent=2), encoding="utf-8"
    )

    # --- Run 7 tau_B ---
    taus = [0.0, 3 / 24, 6 / 24, 12 / 24, 1.0, 2.0, 3.0]
    evs = []
    early_evs = []  # Early-window efficacy
    peak_vs = []    # Peak viral load
    t_peaks = []    # Time to peak
    i_day7s = []    # Infected cells at day 7
    ctrl = simulate(scenario_by_name("run1_r5_control", params), params)
    ctrl_sm = summarize(ctrl["t"], ctrl["y"], params, params["Rs5_0"], params["Rs4_0"])
    auc_c = ctrl_sm["auc_V"]
    auc_c_early = ctrl_sm["auc_V_early"]
    for tau in taus:
        p2 = resolve_params({"tau_B": float(tau)})
        out = simulate(scenario_by_name("run1_r5_circuit", p2), p2)
        sm = summarize(out["t"], out["y"], p2, p2["Rs5_0"], p2["Rs4_0"])
        e = efficacy(sm["auc_V"], auc_c)
        e_early = efficacy(sm["auc_V_early"], auc_c_early)
        evs.append(e)
        early_evs.append(e_early)
        peak_vs.append(sm["peak_V"])
        t_peaks.append(sm["t_peak"])
        i_day7s.append(sm["I_day7"])
        rows.append({"run_id": f"run7_tauB_{tau}", "scenario": "run7", "tau_B": tau, "E_V": e, "E_V_early": e_early, "peak_V": sm["peak_V"], "t_peak": sm["t_peak"], "I_day7": sm["I_day7"], **sm})
    fig, axes = plt.subplots(2, 3, figsize=(15.0, 8.0))
    # 40-day efficacy
    axes[0, 0].plot([x * 24 for x in taus], evs, marker="o")
    axes[0, 0].axhline(crit["antiviral_efficacy"]["E_V_partial"], ls="--", color="C1")
    axes[0, 0].set_xlabel("tau_B (hours)")
    axes[0, 0].set_ylabel("E_V (40-day AUC)")
    axes[0, 0].set_title("HIV recognition delay - 40-day efficacy")
    # Early-window efficacy
    axes[0, 1].plot([x * 24 for x in taus], early_evs, marker="s", color="C2")
    axes[0, 1].axhline(crit["antiviral_efficacy"]["E_V_partial"], ls="--", color="C1")
    axes[0, 1].set_xlabel("tau_B (hours)")
    axes[0, 1].set_ylabel("E_V (days 0-10 AUC)")
    axes[0, 1].set_title("HIV recognition delay - early-window efficacy")
    # Peak viral load
    axes[0, 2].plot([x * 24 for x in taus], peak_vs, marker="^", color="C3")
    axes[0, 2].set_xlabel("tau_B (hours)")
    axes[0, 2].set_ylabel("Peak V (virions/mL)")
    axes[0, 2].set_title("HIV recognition delay - peak viral load")
    # Time to peak
    axes[1, 0].plot([x * 24 for x in taus], t_peaks, marker="d", color="C4")
    axes[1, 0].set_xlabel("tau_B (hours)")
    axes[1, 0].set_ylabel("Time to peak (days)")
    axes[1, 0].set_title("HIV recognition delay - time to peak")

    # False-route onset timing and late-window suppression under run3_false_route,
    # swept over tau_B. This directly tests whether tau_B changes WHEN false suppression
    # begins vs whether it changes the EVENTUAL false-suppression level -- the previous
    # 4-panel version of this figure never computed the false-route scenario at all,
    # despite the caption claiming the panel showed "false-route risk."
    onset_times = []
    late_losses = []
    for tau in taus:
        p2 = resolve_params({"tau_B": float(tau)})
        s3 = scenario_by_name("run3_false_route", p2)
        out3 = simulate(s3, p2)
        t3 = out3["t"]
        Rs5_3 = out3["y"][:, IDX["Rs5"]]
        frac_loss = 1.0 - Rs5_3 / p2["Rs5_0"]
        above80 = np.where(frac_loss >= 0.80)[0]
        onset_times.append(float(t3[above80[0]]) if len(above80) else float("nan"))
        late_mask = t3 > 30.0
        late_losses.append(float(np.mean(frac_loss[late_mask])) if np.any(late_mask) else float("nan"))

    axes[1, 1].plot([x * 24 for x in taus], onset_times, marker="o", color="C5")
    axes[1, 1].set_xlabel("tau_B (hours)")
    axes[1, 1].set_ylabel("Time to 80% false R5 loss (days)")
    axes[1, 1].set_title("False-route suppression onset timing")

    axes[1, 2].plot([x * 24 for x in taus], late_losses, marker="o", color="C6")
    axes[1, 2].set_ylim(0, 1)
    axes[1, 2].set_xlabel("tau_B (hours)")
    axes[1, 2].set_ylabel("Mean R5 loss, days 30-40")
    axes[1, 2].set_title("False-route suppression, late-window (near-asymptotic)")

    plt.tight_layout()
    savefig(fig, "fig_run7_tauB.png")
    tau_crit = None
    for tau, e in zip(taus, evs):
        if e < crit["antiviral_efficacy"]["E_V_partial"]:
            tau_crit = tau
            break
    (SUM / "run7_tauB_critical.json").write_text(json.dumps({"tau_B_critical_days": tau_crit, "taus": taus, "E_V": evs, "E_V_early": early_evs, "peak_V": peak_vs, "t_peak": t_peaks, "I_day7": i_day7s, "false_route_onset_t80": onset_times, "false_route_late_window_loss": late_losses}, indent=2), encoding="utf-8")

    # --- Run 8 A half-life (memory test with noisy physiological bursts) ---
    ths = np.geomspace(5 / 1440, 2.0, 12)
    ev8, f8 = [], []
    for th in ths:
        p2 = resolve_params({"t_half_A": float(th)})
        c_out = simulate(scenario_by_name("run1_r5_circuit", p2), p2)
        # Use noisy bursts scenario instead of constant forcing
        f_out = simulate(scenario_by_name("run8_memory_test", p2), p2)
        smc = summarize(c_out["t"], c_out["y"], p2, p2["Rs5_0"], p2["Rs4_0"])
        smf = summarize(f_out["t"], f_out["y"], p2, p2["Rs5_0"], p2["Rs4_0"])
        ev8.append(efficacy(smc["auc_V"], auc_c))
        f8.append(smf["mean_R5_loss"])
        rows.append({"run_id": f"run8_{th}", "t_half_A": th, "E_V_r5": ev8[-1], "false_R5": f8[-1]})
    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    ax.plot(ths * 24, ev8, marker="o", label="E_V (R5)")
    ax.plot(ths * 24, f8, marker="s", label="false R5 loss (X4+noisy P5)")
    ax.set_xscale("log")
    ax.set_xlabel("t_half_A (hours)")
    ax.legend()
    ax.set_title("Dependence of antiviral efficacy and off-target suppression on route-memory half-life")
    savefig(fig, "fig_run8_A_half_life.png")

    # --- Run 9 ---
    labels = []
    ev9 = []
    for sname, lab in [("run9_fast_only", "fast (internalization) only"), ("run9_sustained_only", "sustained (synthesis+recycling) only"), ("run9_both", "both")]:
        out, sm = run_named(sname, scenario_by_name(sname, params), params)
        rows.append(sm)
        ev9.append(efficacy(sm["auc_V"], auc_c))
        labels.append(lab)
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.bar(labels, ev9, color=["C0", "C1", "C2"])
    ax.set_ylabel("E_V")
    ax.set_title("Contribution of fast and sustained receptor-restriction mechanisms")
    ax.tick_params(axis="x", labelrotation=15)
    plt.tight_layout()
    savefig(fig, "fig_run9_effectors.png")

    # --- Run 10 f_E ---
    fes = np.linspace(0.1, 1.0, 10)
    evf = []
    diagnostic_results = {}  # Store trajectories for diagnostic plots
    for fe in fes:
        p2 = resolve_params({"f_E": float(fe)})
        out = simulate(scenario_by_name("run1_r5_circuit", p2), p2)
        sm = summarize(out["t"], out["y"], p2, p2["Rs5_0"], p2["Rs4_0"])
        e = efficacy(sm["auc_V"], auc_c)
        evf.append(e)
        rows.append({"run_id": f"run10_fE_{fe:.2f}", "f_E": fe, "E_V": e, **sm})
        # Store diagnostic data for key f_E values
        if fe in [0.1, 0.3, 0.5, 0.8, 1.0]:
            diagnostic_results[fe] = out
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.plot(fes, evf, marker="o")
    ax.axhline(crit["antiviral_efficacy"]["E_V_partial"], ls="--")
    ax.axhline(0.0, ls="-", color="k", alpha=0.3)  # Zero efficacy line
    ax.set_xlabel("f_E")
    ax.set_ylabel("E_V")
    ax.set_title("Dependence of antiviral efficacy on engineered cell coverage")
    savefig(fig, "fig_run10_fE.png")
    
    # Diagnostic plots for key f_E values to explain negative efficacy
    fig, axes = plt.subplots(3, 1, figsize=(8.0, 10.0))
    key_fes = [0.1, 0.3, 0.5, 0.8, 1.0]
    colors = ["C0", "C1", "C2", "C3", "C4"]
    for fe, color in zip(key_fes, colors):
        if fe in diagnostic_results:
            out = diagnostic_results[fe]
            t = out["t"]
            T_total = out["y"][:, IDX["T_E"]] + out["y"][:, IDX["T_U"]]
            I = out["y"][:, IDX["I5"]] + out["y"][:, IDX["I4"]]
            V = out["y"][:, IDX["V5"]] + out["y"][:, IDX["V4"]]
            axes[0].plot(t, T_total, label=f"f_E={fe:.1f}", color=color)
            axes[1].plot(t, I, label=f"f_E={fe:.1f}", color=color)
            axes[2].plot(t, V, label=f"f_E={fe:.1f}", color=color)
    
    axes[0].set_ylabel("Total target cells (T_E + T_U)")
    axes[0].set_title("Target cell dynamics by engineered coverage")
    axes[0].legend(fontsize=7)
    axes[1].set_ylabel("Infected cells (I)")
    axes[1].set_title("Infected cell dynamics by engineered coverage")
    axes[1].legend(fontsize=7)
    axes[2].set_ylabel("Virus (V)")
    axes[2].set_xlabel("time (days)")
    axes[2].set_title("Viral dynamics by engineered coverage")
    axes[2].legend(fontsize=7)
    plt.tight_layout()
    savefig(fig, "fig_run10_fE_diagnostics.png")
    
    fcrit = None
    for fe, e in zip(fes, evf):
        if e >= crit["antiviral_efficacy"]["E_V_partial"]:
            fcrit = fe
            break
    (SUM / "run10_fE_critical.json").write_text(json.dumps({"f_E_critical_partial": fcrit, "f_E": list(map(float, fes)), "E_V": evf}, indent=2), encoding="utf-8")

    # comparators
    for sname in ["comparator_ccr5_ko", "comparator_hiv_triggered"]:
        out, sm = run_named(sname, scenario_by_name(sname, params), params)
        sm["E_V"] = efficacy(sm["auc_V"], auc_c)
        rows.append(sm)
    
    # --- Recovery test ---
    recovery_scenario = scenario_by_name("run_recovery_test", params)
    recovery_out, recovery_sm = run_named(
        "run_recovery_test",
        recovery_scenario,
        params,
        extras={
            "recovery_protocol": "two_phase_trigger_source_removal",
            "trigger_removal_day": 14.0,
            "cleared_states_at_release": ["E5", "E4", "I5", "I4", "V5", "V4"],
            "preserved_states_at_release": ["Rs5", "Ri5", "Rs4", "Ri4", "A5", "A4", "B", "C5", "C4", "Q1", "Q2", "Q3"],
        },
    )
    # Recovery is measured only after BOTH scripted route forcing and the active
    # infection source are removed at day 14. Stored circuit/receptor states are
    # preserved, so this is a release/relaxation test rather than a recovery-under-
    # ongoing-infection challenge.
    t_release = 14.0
    t_rec = recovery_time(
        recovery_out["t"], recovery_out["y"],
        params["Rs5_0"], params["Rs4_0"],
        t_off=t_release, receptor="R5",
    )
    recovery_sm["t_recovery_days"] = t_rec
    rows.append(recovery_sm)

    # Plot the clean two-phase protocol: stored synthetic memory, explicit loss
    # of the active infection source, and receptor recovery.
    fig, axes = plt.subplots(3, 1, figsize=(8.2, 10.2), sharex=True)
    t = recovery_out["t"]
    yrec = recovery_out["y"]
    Rs5 = yrec[:, IDX["Rs5"]]
    Rs4 = yrec[:, IDX["Rs4"]]
    A5 = yrec[:, IDX["A5"]]
    B = yrec[:, IDX["B"]]
    C5 = yrec[:, IDX["C5"]]
    I5 = yrec[:, IDX["I5"]]
    V5 = yrec[:, IDX["V5"]]

    axes[0].plot(t, A5, label="A5")
    axes[0].plot(t, B, label="B")
    axes[0].plot(t, C5, label="C5")
    axes[0].axvline(t_release, ls="--", color="k", alpha=0.6, label="Trigger sources removed")
    axes[0].set_ylabel("Synthetic state\n(dimensionless)")
    axes[0].set_title("Stored circuit states decay after complete trigger-source removal")
    axes[0].legend(fontsize=8, ncol=2)
    axes[0].grid(True, alpha=0.3)

    # log10(1+x) allows infected cells and virions to be shown on one finite axis
    # while making the imposed source-removal discontinuity explicit.
    axes[1].plot(t, np.log10(1.0 + I5), label=r"$\log_{10}(1+I_5)$")
    axes[1].plot(t, np.log10(1.0 + V5), label=r"$\log_{10}(1+V_5)$")
    axes[1].axvline(t_release, ls="--", color="k", alpha=0.6, label="Trigger sources removed")
    axes[1].set_ylabel("log-transformed\ninfection state")
    axes[1].set_title("Active infection source is removed at the release event")
    axes[1].legend(fontsize=8, ncol=2)
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(t, Rs5 / params["Rs5_0"], label=r"$R_{s5}/R_{s5,0}$")
    axes[2].plot(t, Rs4 / params["Rs4_0"], label=r"$R_{s4}/R_{s4,0}$")
    axes[2].axvline(t_release, ls="--", color="k", alpha=0.6, label="Trigger sources removed")
    axes[2].axhline(0.9, ls=":", alpha=0.8, label="90% recovery threshold")
    if t_rec is not None:
        axes[2].axvline(t_release + t_rec, ls="-.", alpha=0.8, label=f"R5 recovery: {t_rec:.2f} d")
    axes[2].set_xlabel("time (days)")
    axes[2].set_ylabel("Normalized receptor abundance")
    axes[2].set_title("CCR5 recovery after trigger-source removal")
    axes[2].legend(fontsize=8, ncol=2)
    axes[2].grid(True, alpha=0.3)

    savefig(fig, "fig_recovery_test.png", dpi=300)

    (SUM / "recovery_results.json").write_text(
        json.dumps({
            "model_version": "1.2.1",
            "protocol": "two_phase_trigger_source_removal",
            "trigger_removal_day": t_release,
            "scripted_route_input_after_release": 0.0,
            "infection_source_states_cleared_at_release": ["E5", "E4", "I5", "I4", "V5", "V4"],
            "circuit_and_receptor_states_preserved_at_release": True,
            "t_recovery_days": t_rec,
            "recovery_achieved": t_rec is not None,
            "recovery_threshold": 0.9,
        }, indent=2), encoding="utf-8"
    )

    # --- LHS: manuscript N=128 eight-parameter global screening analysis ---
    n_lhs = 128
    keep = ["t_half_A", "T_A", "tau_B", "k_intC", "k_synC", "eta", "f_E", "k_rec0_R4"]
    ranges = load_sensitivity_ranges(ROOT / "parameters" / "sensitivity_ranges_v1.csv", keep=keep)
    samples = latin_hypercube(n_lhs, ranges, seed=20260827)
    lhs_rows = []
    for i, rec in enumerate(samples):
        p2 = resolve_params(rec)
        out = simulate(scenario_by_name("run3_false_route", p2), p2)
        ctrl_lhs = simulate(scenario_by_name("run2_x4_control", p2), p2)
        sm = summarize(out["t"], out["y"], p2, p2["Rs5_0"], p2["Rs4_0"])
        ctrl_lhs_sm = summarize(ctrl_lhs["t"], ctrl_lhs["y"], p2, p2["Rs5_0"], p2["Rs4_0"])
        e = efficacy(sm["auc_V"], ctrl_lhs_sm["auc_V"])
        lhs_rows.append({**rec, "E_V": e, "false_R5": sm["mean_R5_loss"], "peak_V": sm["peak_V"], "auc_V": sm["auc_V"]})
    lhs_df = pd.DataFrame(lhs_rows)
    lhs_df.to_csv(SUM / "lhs_v1_n128.csv", index=False)

    prcc = []
    parameter_cols = list(ranges)
    for col in parameter_cols:
        rho, pval = partial_rank_correlation(lhs_df, col, "E_V", parameter_cols)
        rho2, p2 = partial_rank_correlation(lhs_df, col, "false_R5", parameter_cols)
        prcc.append({"parameter": col, "PRCC_E_V": rho, "p_E_V": pval, "PRCC_false_R5": rho2, "p_false_R5": p2})
    pd.DataFrame(prcc).to_csv(SUM / "prcc_v1.csv", index=False)
    # Publication Figure 8: two outcomes from the same N=128 global screen.
    # Keep the parameter order and x-axis scale identical so the panels are directly comparable.
    names = [r["parameter"] for r in prcc]
    display_names = {
        "t_half_A": r"$t_{1/2,A}$",
        "T_A": r"$T_A$",
        "tau_B": r"$\tau_B$",
        "k_intC": r"$k_{\mathrm{int},C}$",
        "k_synC": r"$k_{\mathrm{syn},C}$",
        "eta": r"$\eta$",
        "f_E": r"$f_E$",
        "k_rec0_R4": r"$k_{\mathrm{rec},0,4}$",
    }
    labels = [display_names.get(n, n) for n in names]
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6), sharey=True)

    axes[0].barh(labels, [r["PRCC_E_V"] for r in prcc])
    axes[0].axvline(0.0, linewidth=0.8)
    axes[0].set_xlim(-1.0, 1.0)
    axes[0].set_xlabel("PRCC vs $E_V$")
    axes[0].set_title("(a) Antiviral efficacy $E_V$")
    axes[0].grid(axis="x", alpha=0.25)

    axes[1].barh(labels, [r["PRCC_false_R5"] for r in prcc])
    axes[1].axvline(0.0, linewidth=0.8)
    axes[1].set_xlim(-1.0, 1.0)
    axes[1].set_xlabel("PRCC vs false-route R5 suppression")
    axes[1].set_title("(b) False-route R5 suppression")
    axes[1].grid(axis="x", alpha=0.25)

    fig.suptitle(f"N={n_lhs} Latin-hypercube / PRCC global sensitivity screen")
    savefig(fig, "fig_lhs_prcc.png", dpi=300)

    # classify exploratory baseline (not a claim that the circuit works)
    fr = simulate(scenario_by_name("run3_false_route", params), params)
    F_false = summarize(fr["t"], fr["y"], params, params["Rs5_0"], params["Rs4_0"])["mean_R5_loss"]
    # Use actual recovery time from recovery test
    verdict = classify(ev1, F_false, t_rec, crit)
    (SUM / "baseline_verdict.json").write_text(
        json.dumps(
            {
                "note": "Exploratory engineered midpoint — NOT a tuned win. Negative results are retained.",
                "E_V_run1": ev1,
                "false_R5_run3": F_false,
                "class": verdict,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    pd.DataFrame(rows).to_csv(SUM / "all_runs_v1.csv", index=False)
    print("Wrote figures to", FIG)
    print("Baseline exploratory class:", verdict, "E_V", ev1)


if __name__ == "__main__":
    main()
