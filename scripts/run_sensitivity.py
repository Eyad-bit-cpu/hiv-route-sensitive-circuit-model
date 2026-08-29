#!/usr/bin/env python
"""Latin hypercube + PRCC sensitivity analysis.

Default N=128 for a laptop check; use --n 10000 for the registered paper-scale run.
Each sampled circuit run is compared with a matched control using the same sampled
parameter vector, so efficacy does not mix treatment effects with baseline changes.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from metrics import efficacy, summarize
from scenarios import scenario_by_name
from simulate import resolve_params, simulate
from stats_utils import latin_hypercube, load_sensitivity_ranges, partial_rank_correlation


def control_name_for(scenario_name: str) -> str:
    """Return the matched untreated comparator for a challenge scenario."""
    if scenario_name in {"run2_x4_circuit", "run3_false_route", "run8_memory_test"}:
        return "run2_x4_control"
    return "run1_r5_control"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=128)
    ap.add_argument("--seed", type=int, default=20260827)
    ap.add_argument("--scenario", default="run3_false_route")
    args = ap.parse_args()

    ranges = load_sensitivity_ranges(ROOT / "parameters" / "sensitivity_ranges_v1.csv")
    if args.n < 1000:
        keep = ["t_half_A", "T_A", "tau_B", "k_intC", "k_synC", "eta", "f_E", "k_rec0_R4"]
        ranges = {k: ranges[k] for k in keep if k in ranges}

    samples = latin_hypercube(args.n, ranges, args.seed)
    recs = []
    ctrl_name = control_name_for(args.scenario)
    for rec in samples:
        p = resolve_params(rec)
        out = simulate(scenario_by_name(args.scenario, p), p)
        ctrl = simulate(scenario_by_name(ctrl_name, p), p)
        sm = summarize(out["t"], out["y"], p, p["Rs5_0"], p["Rs4_0"])
        ctrl_sm = summarize(ctrl["t"], ctrl["y"], p, p["Rs5_0"], p["Rs4_0"])
        recs.append({**rec, **sm, "E_V": efficacy(sm["auc_V"], ctrl_sm["auc_V"]), "false_R5": sm["mean_R5_loss"]})

    df = pd.DataFrame(recs)
    outp = ROOT / "results" / "summaries" / f"lhs_n{args.n}.csv"
    outp.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(outp, index=False)

    parameter_cols = list(ranges)
    prcc = []
    for col in parameter_cols:
        rho_ev, p_ev = partial_rank_correlation(df, col, "E_V", parameter_cols)
        rho_false, p_false = partial_rank_correlation(df, col, "false_R5", parameter_cols)
        prcc.append({
            "parameter": col,
            "PRCC_E_V": rho_ev,
            "p_E_V": p_ev,
            "PRCC_false_R5": rho_false,
            "p_false_R5": p_false,
        })
    prcc_df = pd.DataFrame(prcc)
    prcc_df.to_csv(ROOT / "results" / "summaries" / f"prcc_n{args.n}.csv", index=False)

    # For the registered false-route scenario, visualize both paper outcomes side by side.
    if args.scenario == "run3_false_route":
        names = prcc_df["parameter"].tolist()
        fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6), sharey=True)
        axes[0].barh(names, prcc_df["PRCC_E_V"].to_numpy())
        axes[0].axvline(0.0, linewidth=0.8)
        axes[0].set_xlim(-1.0, 1.0)
        axes[0].set_xlabel("PRCC vs $E_V$")
        axes[0].set_title("(a) Antiviral efficacy $E_V$")
        axes[0].grid(axis="x", alpha=0.25)

        axes[1].barh(names, prcc_df["PRCC_false_R5"].to_numpy())
        axes[1].axvline(0.0, linewidth=0.8)
        axes[1].set_xlim(-1.0, 1.0)
        axes[1].set_xlabel("PRCC vs false-route R5 suppression")
        axes[1].set_title("(b) False-route R5 suppression")
        axes[1].grid(axis="x", alpha=0.25)

        fig.suptitle(f"N={args.n} Latin-hypercube / PRCC global sensitivity screen")
        fig.tight_layout()
        fig.savefig(ROOT / "results" / "figures" / f"fig_prcc_n{args.n}.png", dpi=200)
        plt.close(fig)

    print(outp)


if __name__ == "__main__":
    main()
