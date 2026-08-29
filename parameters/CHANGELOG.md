## Model v1.2.1 audit hardening (2026-08-28)

- Added complete mathematical documentation of the HIV-confirmation state B, including infected-cell input, saturating transform, three-stage delay equations, and separate delay/persistence parameters.
- Added a graded Run-3 wrong-route forcing sweep scaled to the peak infection-derived route-input in the baseline R5 simulation. The registered P5=1.2 stress condition is ~9.56x that internal model reference and is now explicitly classified as an extreme adversarial stress test rather than a physiological prediction.
- Relabeled P5/P4 registry descriptions as uncalibrated effective scenario forcings rather than physiological activity measurements.
- Renamed Run 3 figure language to emphasize robust persistent false-route suppression rather than a nonexistent efficacy-specificity trade-off.
- Expanded Run 4 figures to expose scripted $P_5(t)$ versus infection-driven $\eta\,\mathrm{inf}_{5,E}$, plus $B$, $C_5$, and CCR5 feedback.
- Canonicalized Run 6 as `run6_secondary_x4_emergence`; legacy `run6_tropism_shift` remains a compatibility alias. Output is now `run6_crossover_time.json` and the figure is `fig_run6_secondary_x4_emergence.png`.
- Figure 8 now uses manuscript mathematical notation for parameter labels.
- Removed remaining v1.2/"memory" terminology ghosts where $\tau_B$ is specifically an HIV-confirmation delay.
- Added regression tests covering mixed-route undefined-share masking, Run-6 emergence semantics, $\tau_B$ onset-vs-late suppression, Run-9 arm isolation, and Run-3 robust false-route failure.

## Model v1.2.1 (2026-08-27) - Recovery Protocol Correction

- Corrected the registered recovery test so that day-14 release removes **both** scripted route forcing and the active infection source. The previous scenario removed only the scripted `P5` forcing while live infection remained able to regenerate `B`, confounding a clean reversibility test.
- At release, infection-source states (`E5`, `E4`, `I5`, `I4`, `V5`, `V4`) are set to zero; receptor states, `A/B/C`, and `Q1-Q3` are preserved and decay/recover dynamically.
- Added a regression test verifying that infection cannot regrow after release and that stored circuit memory is not artificially zeroed.
- Recovery Figure now displays circuit-state decay, explicit infection-source removal, and receptor recovery.
- No biological baseline parameter, success threshold, or N=128 sensitivity sample was retuned.

# Parameter-set versions

`baseline_v1.csv` is frozen. If a literature-anchored biological baseline must change, create `baseline_v2.csv` and record the reason here. Never overwrite v1 because a figure looked bad.

## v1 (2026-08-27)

- First freeze of Model v1.0 Class A biology.
- Infection-rate constant `beta` is **derived** from Ribeiro et al. 2010 median R0 = 8 together with the frozen T0, p, c, delta_I pair — not tuned to make the circuit win.
- Engineered-circuit numbers in `baseline_v1.csv` are **exploratory midpoints**, not biological estimates. They are not cited as if measured.

## Model v1.1.0 (2026-08-27) - Major Technical Fixes

### Documentation Corrections
- Removed incorrect claim that A,B,C ∈ [0,1]. Updated to reflect that synthetic state variables are dimensionless effective signaling variables without a priori bounds.

### Statistical Analysis Fixes
- Implemented actual Partial Rank Correlation Coefficient (PRCC) instead of Spearman rank correlation.
- Previous implementation was simply Spearman correlation, not true partial correlation controlling for other parameters.

### Logic Bug Fixes
- Fixed recovery classification: missing recovery time (t_rec = None) now correctly prevents full success classification.
- Fixed constitutive CCR5 comparator: now clamps both surface and internal receptor pools to match infection dynamics.
- Fixed suppression-share metric: added numerator clipping to prevent negative values when receptors exceed baseline.

### Scenario Improvements
- Fixed Run 3 false-route efficacy to use X4 control AUC instead of R5 control AUC.
- Redesigned Run 8 with transient/noisy physiological activation to properly test temporal memory discrimination.
- Added early-window endpoints to τ_B analysis (days 0-10 AUC, viral peak, time to peak).
- Added diagnostic plots for f_E analysis to explain negative efficacy at intermediate coverage.
- Added dedicated recovery/OFF scenario with explicit recovery time measurement.

### Figure Improvements
- Updated all figure titles to publication style.
- Normalized receptor axes (R_s/R_s,0) throughout.
- Improved panel layouts for better temporal logic visualization.

### Code Quality
- Added explicit positivity checks with meaningful error reporting.
- Added cross-platform ZIP packaging script.
- Updated model version to 1.1.0.

## Model v1.2.0 (2026-08-27) - Correctness and Reproducibility Repair

### Core dynamics
- Corrected the two-compartment receptor steady-state algebra. Baseline internal pools now use `k_int*Rs/(k_rec+k_deg)`, and baseline synthesis uses `k_int*k_deg*Rs/(k_rec+k_deg)`.
- Constitutive CCR5 comparator is now initialized at the clamped surface/internal state, so the stored trajectory and the effective infection calculation agree.
- Positivity enforcement moved from intermediate RHS trial evaluations to a post-integration fail-closed trajectory check; tiny solver-scale excursions are clipped only after a successful integration.

### Metrics and statistics
- Recovery is now measured for an explicitly selected receptor (`R5`, `R4`, or `both`); an untouched receptor can no longer create a false immediate recovery.
- PRCC is centralized in `src/stats_utils.py`: rank-transform → residualize against other sampled parameters only → Pearson correlation of residuals, with partial-correlation degrees of freedom for p-values.
- Campaign and sensitivity scripts share the same PRCC and Latin-hypercube implementations.

### Controls and sensitivity
- Run-3 LHS now uses a matched X4 control for every sampled parameter vector.
- `run_sensitivity.py` now recomputes the appropriate matched control for every LHS sample, including full biological-parameter sensitivity runs.
- `run_sweep.py` now recomputes the selected control for each swept parameter value.
- Campaign LHS ranges are loaded from `parameters/sensitivity_ranges_v1.csv` instead of being duplicated in code.

### Reproducibility and regression protection
- Added steady-state, constitutive-clamp, recovery-target, unknown-override, PRCC, and config-alignment tests.
- Run-8 stochastic forcing now derives a stable per-scenario seed from the registered master seed.
- Git commit lookup no longer prints stderr noise outside a Git checkout.
- Documentation updated to Model v1.2, corrected state count, executable figure path, and actual solver-consistency test location.
