#!/usr/bin/env python
"""One-parameter sweep. Example:
python scripts/run_sweep.py t_half_A --lo 0.003472 --hi 2 --n 12 --log
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from metrics import efficacy, summarize
from scenarios import scenario_by_name
from simulate import resolve_params, simulate


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("parameter")
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--scenario", default="run1_r5_circuit")
    ap.add_argument("--control", default="run1_r5_control")
    ap.add_argument("--log", action="store_true")
    ap.add_argument("--lo", type=float, required=True)
    ap.add_argument("--hi", type=float, required=True)
    args = ap.parse_args()
    # Controls are recomputed for each sampled parameter value so efficacy is
    # matched on the same underlying parameter vector.
    xs = np.geomspace(args.lo, args.hi, args.n) if args.log else np.linspace(args.lo, args.hi, args.n)
    recs = []
    for x in xs:
        p = resolve_params({args.parameter: float(x)})
        out = simulate(scenario_by_name(args.scenario, p), p)
        ctrl = simulate(scenario_by_name(args.control, p), p)
        sm = summarize(out["t"], out["y"], p, p["Rs5_0"], p["Rs4_0"])
        ctrl_sm = summarize(ctrl["t"], ctrl["y"], p, p["Rs5_0"], p["Rs4_0"])
        sm[args.parameter] = float(x)
        sm["E_V"] = efficacy(sm["auc_V"], ctrl_sm["auc_V"])
        recs.append(sm)
    outp = ROOT / "results" / "summaries" / f"sweep_{args.parameter}.csv"
    outp.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(recs).to_csv(outp, index=False)
    print(outp)


if __name__ == "__main__":
    main()
