# Figure 8 patch

This archive updates the N=128 LHS/PRCC visualization used as Figure 8 in the manuscript.

Changes:
- `scripts/run_campaign.py` now generates a two-panel PRCC figure with:
  - panel (a): PRCC vs antiviral efficacy `E_V`;
  - panel (b): PRCC vs false-route R5 suppression;
- both panels use the same parameter ordering and the same x-axis limits `[-1, 1]`;
- `scripts/run_sensitivity.py` now preserves and reports both PRCC outcomes for the false-route scenario and can generate the same two-panel layout;
- `results/figures/fig_lhs_prcc.png` was regenerated from the existing frozen N=128 PRCC table;
- no parameter values, simulation equations, N=128 samples, PRCC coefficients, or p-values were changed.

Verification: `pytest -q` -> 21 passed.
