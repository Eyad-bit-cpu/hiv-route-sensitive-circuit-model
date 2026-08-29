# Model v1.2.1 forensic-audit hardening

This release implements the six remaining issues identified in the post-recovery forensic audit without retuning biological parameters or frozen success criteria.

1. **Run 3 / Figure 4**: removed efficacy-specificity "trade-off" framing. The figure now emphasizes that false-route suppression remains persistently high while efficacy varies across the registered 2D sweep.
2. **Run 4 / Figure 5**: retained the scientific question of external forcing during ongoing infection, but made the feedback explicit by plotting scripted `P5(t)` separately from infection-driven `eta*inf5_E`, together with `B`, `C5`, and normalized CCR5.
3. **Run 6**: canonicalized the experiment as secondary X4 emergence into established R5 infection. New canonical source/output names are `run6_secondary_x4_emergence`, `fig_run6_secondary_x4_emergence.png`, and `run6_crossover_time.json`; the former scenario name remains a backwards-compatible alias.
4. **Figure 8**: replaced code-style parameter names with manuscript mathematical notation.
5. **Terminology/manuscript cleanup**: removed residual `v1.2` state-count wording, removed `tau_B` memory terminology where it referred to recognition delay, and removed internal bug-history prose from the scientific manuscript.
6. **Regression hardening**: added tests for Run-5 undefined-share masking, Run-6 secondary-X4 semantics/crossover, Run-7 onset-vs-late suppression behavior, Run-9 true arm isolation, and Run-3 robust false-route failure.

The full campaign completes successfully after these changes.
