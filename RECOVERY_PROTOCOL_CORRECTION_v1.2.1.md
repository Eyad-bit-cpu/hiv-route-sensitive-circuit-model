# Recovery protocol correction - Model v1.2.1

The v1.2.0 recovery scenario was not a clean trigger-removal experiment. It ended scripted R5 forcing at day 14, but retained live infected-cell and virion states. Because infection itself drives HIV-confirmation signal B, receptor recovery could re-enable entry, amplify infection, and re-trigger suppression. The resulting endpoint therefore mixed receptor reversibility with an ongoing infection feedback loop.

Model v1.2.1 uses a two-phase release protocol. Phase 1 is unchanged through day 14. At day 14, scripted route forcing is set to zero and the active infection-source states E5/E4/I5/I4/V5/V4 are cleared. Receptor states, A5/A4, B, C5/C4, and Q1-Q3 are preserved exactly at their release values and then evolve under the unchanged ODEs. Thus any already-stored route/confirmation memory must decay naturally, and receptor recovery is measured from the release event to the first time Rs5 reaches 90% of baseline.

With the unchanged baseline parameters, R5 reaches 90% of baseline approximately 5.86 days after release, satisfying the registered 14-day recovery criterion. The overall baseline classification remains FAILED because the false-route suppression criterion is still violated.
