# Simulation protocol — Model v1.2.1

## What is being simulated

A **mechanistic within-host effective-compartment model of early HIV infection** in a representative susceptible CD4⁺ target-cell compartment. Quantities are concentrations per mL of effective compartment. This is **not** an anatomically resolved human body.

The synthetic circuit is grafted onto a literature-backed HIV backbone (eclipse-phase cells, saturating coreceptor infectability, surface/internal receptor pools).

## Units

Master time unit: **day**. Example: 5 min = 0.003472 day. Solver: adaptive BDF/LSODA (`scipy.integrate.solve_ivp`, method `LSODA` by default). Benchmarks may be repeated with `BDF`.

Synthetic species `A5`, `A4`, `B`, `C5`, `C4` are **dimensionless effective signaling variables**; they are not constrained a priori to the interval `[0, 1]`. Their downstream actions are bounded where appropriate through Hill-response or saturating regulatory functions. Cells, virions, and receptor molecules/cell remain physical.

## Parameter classes

| Class | May we tune it so the circuit “wins”? | Role |
| --- | --- | --- |
| A — literature-anchored biology | **No.** Change only by issuing `baseline_v2.csv` with a written reason. | Priors/ranges. |
| B — engineered circuit | **No secret tuning.** Sweep and **report the required region**. | Paper predictions. |
| C — scenario / stress tests | Designed challenges, not hidden knobs. | Adversarial worlds. |

## Pre-registered success criteria (frozen before headline sweeps)

Defined in `parameters/success_criteria_v1.yaml`. Do not move goalposts after viewing sweep heatmaps.

- **Successful:** `E_V >= 0.5` AND false-route suppression `F <= 0.25` AND recovery time `t_rec < 14 d` in a dedicated two-phase release scenario. At the release event, both the scripted route input and the active infection-source states (`E5/E4/I5/I4/V5/V4`) are removed, while receptor, `A/B/C`, and delay-chain states are preserved and allowed to relax. This isolates post-trigger reversibility rather than recovery during continuing infection.
- **Partial:** `E_V >= 0.25` with `F <= 0.40`, or viral control without recovery/false-route compliance.
- **Failed:** otherwise.

`E_V = 1 - AUC_V,circuit / AUC_V,control` over the registered horizon.

## Negative-result policy

Parameter regions producing poor viral control, excessive physiological receptor suppression, unstable dynamics, failure to track tropism, or unrealistic kinetic requirements are **retained and reported**, not excluded.

## Run hierarchy

| Phase | Runs | Gate |
| --- | --- | --- |
| I Biological calibration | 0A HIV backbone; 0B CCR5 trafficking; 0C CXCR4 trafficking | If 0A is biologically insane, **stop**. |
| II Circuit unit tests | 0D A-memory; 0E AND gate | Modules must pass tests. |
| III Nominal challenges | 1 R5; 2 X4; 5 mixed | Descriptive, not optimized. |
| IV Adversarial | 3 false-route; 4 pulses; 6 secondary X4 emergence; 7 τ_B | Includes the central heatmap. |
| V One-at-a-time sweeps | 8 A half-life; 9 internalization vs transcription; 10 f_E | Discover constraints. |
| VI Global sensitivity | LHS (+ optional Sobol) | Interactions. |
| VII Comparators | constitutive CCR5 KO; HIV-triggered non-route; no circuit | Advantage test. |
| VIII Freeze + figures | `scripts/run_campaign.py` | Paper panels from saved `run_id`s only. |

## Reproducibility metadata (every run)

`run_id`, `git_commit` (or `NO_GIT`), `random_seed`, `model_version`, `solver`, `rtol`, `atol`, `t_start`, `t_end`, `parameter_file_hash`, `scenario`, plus the full resolved parameter vector.

Master seed: **20260827**. Child seeds are `sha256(master + run_id)` truncated to uint32.

## Solver cross-check

Run 0A and Run 1 control are repeated with `LSODA` and `BDF`. Trajectories must agree within a relative envelope documented in `tests/test_model.py::test_solver_lsoda_vs_bdf`.


## Run 3 forcing interpretation and graded robustness analysis (v1.2.1 reviewer hardening)

The registered persistent wrong-route stress condition uses constant \(P_5=1.2\) for the full 40-day X4-only simulation. \(P_5\) is a dimensionless effective route-associated forcing and is **not calibrated to a measured physiological CCR5 signal**. In the baseline route-matched R5 circuit simulation, the peak infection-derived route input \(\eta\,\mathrm{inf}_{5,E}\) is approximately 0.12549, so the registered stress forcing is approximately 9.56-fold that internal model reference peak. It is therefore interpreted as an intentionally extreme adversarial stress test rather than a prediction of ordinary physiological co-activation.

The campaign also runs a graded constant-forcing sweep and writes `results/summaries/run3_wrong_route_forcing_sweep.csv` plus `run3_wrong_route_forcing_reference.json`. This sweep is scaled to the model's own infection-derived reference and is used only to map the onset/severity of false suppression within the effective model; it does not establish a physiological conversion for \(P_5\).
