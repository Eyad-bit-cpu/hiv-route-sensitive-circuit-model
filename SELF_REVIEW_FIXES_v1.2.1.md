# Self-review major-comment fixes - v1.2.1

1. **HIV-confirmation state B**
   - Main manuscript now defines the full infected-cell input, Hill transform, Q1-Q3 delay chain, tau_B=0 bypass, B decay, and baseline effective parameters.
   - tau_B is explicitly a confirmation-production delay; t_half_B is the B persistence timescale.
   - Biological discussion uses HIV RNA/DNA sensing literature only to demonstrate plausible *information classes*, not to claim a ready-made molecular implementation.

2. **Persistent wrong-route challenge**
   - Registered stress condition is explicitly P5=1.2 for 40 days.
   - Baseline route-matched R5 simulation gives peak eta*inf5_E = 0.125486654, so P5=1.2 is 9.56277x that internal model reference.
   - P5 is uncalibrated to physiology; 80.1% false suppression is therefore described as an extreme adversarial robustness failure.
   - New graded forcing sweep maps false suppression from background through the stress regime and is included in source results, main Figure 4, and supplementary data.

Verification: source test suite = 28 passing tests; main and supplementary PDFs compile with the documented XeLaTeX route.
