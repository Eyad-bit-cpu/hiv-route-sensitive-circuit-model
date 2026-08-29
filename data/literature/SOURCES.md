# Literature sources used to freeze Model v1.0 Class A parameters

This file records **what was actually read** for biological anchors. Engineered-circuit parameters have **no** experimental source and must not be given fake citations.

## HIV within-host backbone

| Claim | Value used | Source | Notes |
| --- | --- | --- | --- |
| Standard T–I–V ODE structure | — | Perelson & Ribeiro, *BMC Biol* 2013; 11:96. PMID 24020860. PMC3765939 | Review of the classical within-host model. |
| Uninfected target-cell death `d_T` | 0.01 day⁻¹ | Kadelka et al., *PLoS Comput Biol* 2024; 20(6):e1012129. PMID 38900725. PMC11189221 | Fixed `d = 0.01 /day` in acute-HIV fits; they cite prior viral-dynamics practice. |
| Effective susceptible density `T0` | 10⁶ cells/mL | Same paper | Fixed `T0 = 10^6` cells/mL, “assuming that 1 in 1000 cells are available as targets”. |
| Target supply | `lambda = d_T * T0` | Same paper | Pre-infection equilibrium. |
| Plasma virion clearance `c` | 23 day⁻¹ (range 9.1–36) | Ramratnam et al., *Lancet* 1999; 354:1782–1786. PMID 10577640 | Apheresis: `c` 9.1–36 day⁻¹, mean ~23 day⁻¹; half-life 28–110 min. |
| Earlier lower `c` | ~3 day⁻¹ | Perelson et al., *Science* 1996; 271:1582–1586. PMID 8599114 | Mean virion life-span 0.3 d under treatment-perturbation analysis. **Kept as uncertainty, not overwritten.** |
| Infected-cell loss `delta_I` (typical simulation value) | 0.7 day⁻¹ | Ikeda et al., *Theor Biol Med Model* 2014; 11:22. PMID 24886060. PMC4035760 | Paper uses typical `δ = 0.7 /day`, `p = 4000 /day`, `c = 23 /day` in methods; mouse fits averaged `δ ≈ 0.61 /day` (range ~0.3–0.76). |
| Infected-cell life-span (human, 1996) | 2.2 d → δ ≈ 0.45 day⁻¹ | Perelson et al. 1996, PMID 8599114 | Lower bound of our Class A range. |
| Viral production `p` | 4000 virions cell⁻¹ day⁻¹ | Ikeda et al. 2014 typical parameterization | **Not** an independent human burst-size measurement. Burst-size literature is larger (see Chen et al. below); `p` and infectious fraction are partially non-identifiable with `c` and `β`. |
| Eclipse / generation time | mean eclipse 1 d; min lifecycle 1.2 d; generation 2.6 d | Perelson et al. 1996; Ribeiro et al. 2010 used a 24 h eclipse | `k_E = 1.0 day⁻¹` is a mapping of a ~1-day intracellular delay. |
| Acute `R0` (used to **derive** `β`, not to fit the circuit) | median 8.0 (IQR 4.9–11) | Ribeiro, Qin, Chavez, Li, Self, Perelson, *J Virol* 2010; 84:6096–6102. PMID 20357090. PMC2876646 | 47 plasma donors. With 24 h eclipse, mixed-effects group `R0 ≈ 8.8`. |
| Acute peak VL (validation target, not a fitted cost) | median 5.8 log10 copies/mL | Ribeiro et al. 2010 | Peak ~14 days after virus became quantifiable (LOD 50 copies/mL), not days post-exposure. |
| Another acute `R0` estimate | mean 7.1 (no delay) / 19.3 (24 h delay) | Little et al., *J Exp Med* 1999; 190:841–850. PMID 10499922 | Shows `R0` depends on eclipse assumptions — carried as uncertainty. |

Chen HY, Di Mascio M, Perelson AS, Ho DD, Zhang L. Determination of virus burst size in vivo using a single-cycle SIV in rhesus macaques. *PNAS* 2007. PMID 18000036. Burst size on the order of 5×10⁴ virions/cell. We **do not** retune `p` to this number in v1 because infectious-to-RNA ratio is unknown; it is recorded as an uncertainty bound.

## Coreceptor abundance and infectability (not linear in R)

| Claim | Value used | Source | Notes |
| --- | --- | --- | --- |
| CCR5 molecules/cell on primary CD4⁺ T cells | thousands; inter-individual ~5-fold; examples ~4×10³–1.9×10⁴ | Reynes et al., *J Infect Dis* 2000; 181:927–932. PMID 10720514 | Quantitative flow; density correlated with plasma RNA (`r = 0.666`, `P = 0.009`). Examples in text: 19,147 and 9,132 (controls); 4,116 (lowest WT/WT control); 9,461/9,601 (patients). |
| Saturation / threshold vs CCR5 density | high-CD4 cells: 7×10²–2×10³ CCR5/cell already maximal; low-CD4 cells: ~1–2×10⁴ needed | Platt, Wehrly, Kuhmann, Chesebro, Kabat, *J Virol* 1998; 72:2855–2864. PMID 9525605 | HeLa-CD4/CCR5 clones, 7.0×10²–1.3×10⁵ CCR5/cell. **This forbids a linear `β R V T` infection term.** |
| Low-CD4 threshold language | ~10⁴ CCR5/cell | Reynes et al. 2000 discussing Platt | Used as an **alternative** Hill `K_R5` scenario, not the only truth. |
| Chemotaxis vs CCR5 density | not identical to infection Hill | Desmetz et al., *J Immunol* 2005 / PMC2265826 | Cited for the broader point that CCR5 density maps nonlinearly onto function. |

v1 baseline Hill `K_R5 = 1500` molecules/cell is a **mapping** of Platt’s high-CD4 saturation window (700–2000), because the target cells are CD4⁺ T cells. It is **not** a statistical fit of Platt’s foci data. `K_R5 = 10000` is provided as a documented alternative (Reynes’s reading of Platt’s low-CD4 threshold).

CXCR4 molecules/cell: no single primary-cell mean is frozen as gospel. v1 uses `Rs4_0 = 20000` molecules/cell as an order-of-magnitude value (CXCR4 is readily detected and often more abundant than CCR5 on naïve/resting T-cell subsets) with a wide range. `K_R4` is **not** given a fake “CXCR4-Platt” citation; it is an analogous saturating map with high uncertainty.

## Receptor trafficking timescales

| Claim | Value used | Source | Notes |
| --- | --- | --- | --- |
| CXCR4 constitutive internalization | 1.0% of surface pool per minute | Signoret et al., *J Cell Biol* 1997; 139:651–664. PMID 9348282. PMC2141706 | Human T-cell lines SupT1 and BC7. |
| SDF-1 down-modulation | 50% in ~5 min; ~20% remaining at 30 min | Same paper, Fig. 8A, 500 nM SDF-1 | Used to map ligand-stimulated `k_int` for Run 0C. |
| PMA | endocytosis rate increased >6-fold; 60–90% surface loss over 120 min; recycling after PMA removal | Same paper | Recycling exists; exact 1–3 h recovery is **not** a single number in this paper. |
| CCR5 constitutive protein turnover | t½ 6–9 h | Signoret et al., *J Cell Biol* 2000; 151:1281–1298. PMID 11121442. PMC2190598 | CHX chase in CCR5-CHO cells; similar to Mirzabekov et al. 1999. This is **total cellular CCR5 content**, not a plasma-membrane-only rate. |
| Ligand-induced CCR5 endocytosis + recycling after RANTES removal | qualitative: rapid internalization; recycle after RANTES, not AOP-RANTES | Same 2000 paper | Quantitative minute-scale `k_int` for CCR5 is less tightly pinned than CXCR4/SDF-1; Run 0B uses a documented mapping (see registry `k_int_lig_R5`) with high uncertainty. |
| Distinct CXCR4 vs CCR5 endocytic signals | — | Signoret et al., *J Cell Sci* 1998; 111:2819–2830. PMID 9730970 | CCR5 internalizes with RANTES, not PMA; CXCR4 has separable ligand vs PMA signals. |

## What simulation is allowed to discover (Class B)

Half-lives of `A`, thresholds `T_A`, HIV-confirmation delay `tau_B`, suppressor strength, engineered fraction `f_E`, and Hill coefficients of the **synthetic** coincidence gate are **not** imported from HIV or GPCR papers. They are sweep variables. Baseline values are log-range midpoints so the code runs, and they are labeled `class=engineered`.

## Calibration vs testing (Level 2)

- **Trafficking calibration (train):** Signoret 1997 CXCR4 constitutive %/min and SDF-1 5-min half-downmodulation; Signoret 2000 CCR5 6–9 h protein turnover.
- **Trafficking test (held out from fitting):** qualitative recycling after stimulus removal (Signoret 1997 PMA washout; Signoret 2000 RANTES washout). We do **not** numerically optimize recycling rates to both papers at once.
- **HIV backbone calibration:** `d_T`, `T0`, `c`, `p`, `delta_I` frozen from the table above; `beta` derived from Ribeiro median `R0`.
- **HIV backbone test (not used to retune):** order-of-magnitude acute peak (~10⁶–10⁷ RNA/mL) and peak timing (days–weeks), compared to Ribeiro 2010 / Stafford-type primary infection, with the explicit caveat that this effective-compartment model is **not** a plasma-only clinical simulator.

## HIV-associated information classes relevant to abstract confirmation state B

These papers are used only to support the statement that HIV-associated molecular information can be sensed intracellularly; they do not calibrate or specify an implementation of the synthetic B state.

- Berg RK et al. 2012. *Genomic HIV RNA induces innate immune responses through RIG-I-dependent sensing of secondary-structured RNA*. PLoS ONE 7:e29291. DOI 10.1371/journal.pone.0029291. PMID 22235281.
- Jakobsen MR et al. 2013. *IFI16 senses DNA forms of the lentiviral replication cycle and controls HIV-1 replication*. PNAS 110:E4571-E4580. DOI 10.1073/pnas.1311669110. PMID 24154727.
- Eschbach JE et al. 2024. *HIV-1 capsid stability and reverse transcription are finely balanced to minimize sensing of reverse transcription products via the cGAS-STING pathway*. mBio 15:e00348-24. DOI 10.1128/mbio.00348-24. PMID 38530034.
