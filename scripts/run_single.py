#!/usr/bin/env python
"""Run one named scenario. Example: python scripts/run_single.py run1_r5_circuit"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from metrics import summarize
from scenarios import scenario_by_name
from simulate import resolve_params, save_run, simulate


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scenario")
    ap.add_argument("--run-id", default=None)
    ap.add_argument("--method", default="LSODA", choices=["LSODA", "BDF", "Radau"])
    args = ap.parse_args()
    params = resolve_params()
    sc = scenario_by_name(args.scenario, params)
    rid = args.run_id or args.scenario
    out = simulate(sc, params, method=args.method)
    path = save_run(rid, sc, out)
    sm = summarize(out["t"], out["y"], params, params["Rs5_0"], params["Rs4_0"])
    print(json.dumps({"saved": str(path), **sm}, indent=2))


if __name__ == "__main__":
    main()
