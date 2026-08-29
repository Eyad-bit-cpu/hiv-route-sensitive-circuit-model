# Regenerated Results — Model v1.2.1

Generated from the repaired Model v1.2.1 codebase by running:

```bash
python scripts/run_campaign.py
```

## Regeneration status

- Campaign: completed successfully
- Fresh figures: 22
- Summary files: 8
- Raw run directories: 22
- Baseline registered class: **FAILED**

## Core regenerated outputs

| Quantity | v1.2.1 regenerated value |
|---|---:|
| Run 1 R5 control AUC | 2.67135e+08 |
| Run 1 R5 circuit AUC | 1.46171e+08 |
| Run 1 R5 efficacy E_V | 0.452821 |
| Run 1 R5 control peak V | 5.79139e+07 |
| Run 1 R5 circuit peak V | 5.26307e+06 |
| Run 1 mean R5 loss | 0.699530 |
| Run 1 mean off-route R4 loss | 0.007438 |
| Run 2 X4 control AUC | 2.67135e+08 |
| Run 2 X4 circuit AUC | 1.92061e+08 |
| Run 2 X4 efficacy E_V | 0.281034 |
| Run 2 mean R4 loss | 0.690143 |
| Run 2 mean off-route R5 loss | 0.017301 |
| Baseline Run-3 false R5 suppression | 0.801419 |
| Clean-release R5 recovery to 90% | 5.85 days |
| Recovery achieved within 14-day criterion? | True |

## Sensitivity snapshot

- Strongest absolute PRCC with E_V: `eta` = 0.5199 (p=9.88e-10)
- Strongest absolute PRCC with false-route R5 loss: `tau_B` = -0.9802 (p=1.43e-85)

## Figure inventory

- `fig_F5_hill_platt_mapping.png`
- `fig_lhs_prcc.png`
- `fig_recovery_test.png`
- `fig_run0a_V5_log.png`
- `fig_run0a_backbone.png`
- `fig_run0b_ccr5.png`
- `fig_run0c_cxcr4.png`
- `fig_run0d_A.png`
- `fig_run10_fE.png`
- `fig_run10_fE_diagnostics.png`
- `fig_run1_r5.png`
- `fig_run2_x4.png`
- `fig_run3_false_route_heatmap.png`
- `fig_run3_robust_false_route.png`
- `fig_run4_isolated_pulse.png`
- `fig_run4_repeated_pulses.png`
- `fig_run4_sustained.png`
- `fig_run5_mixed.png`
- `fig_run6_secondary_x4_emergence.png`
- `fig_run7_tauB.png`
- `fig_run8_A_half_life.png`
- `fig_run9_effectors.png`

## Important interpretation note

The v1.2.1 recovery protocol additionally removes both scripted route forcing and the active infection source at day 14 while preserving receptor, A/B/C, and delay-chain states. This isolates post-trigger relaxation and yields R5 recovery to 90% in approximately 5.85 days. All non-recovery headline values and the N=128 PRCC screen are unchanged. No parameters or thresholds were retuned to force a favorable outcome.
