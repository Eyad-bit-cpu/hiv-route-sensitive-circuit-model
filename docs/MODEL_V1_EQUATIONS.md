# Model v1.2.1 equations

Time unit: day. `A, B, C` dimensionless. Receptors in molecules/cell. Cells and virions per mL.

## Infectability (not linear in R)

\[
F_x(R_{s,x}) = \frac{R_{s,x}^{h_x}}{K_{R,x}^{h_x} + R_{s,x}^{h_x}}
\]

Baseline `K_{R5} = 1500` maps Platt et al. 1998 high-CD4 saturation (700–2000 molecules/cell). Alternative `K_{R5} = 10000` is Reynes’s reading of the low-CD4 threshold and is **not** the v1 default.

## Infection (engineered vs unmodified)

\[
\mathrm{inf}_{5E} = \beta F_5(R_{s5}) T_E V_5,\quad
\mathrm{inf}_{5U} = \beta F_5(R_{s5,0}) T_U V_5
\]

and analogously for X4. Unmodified cells keep baseline surface density.

\[
\dot T_E = f_E\lambda - d_T T_E - \mathrm{inf}_{5E} - \mathrm{inf}_{4E}
\]

\[
\dot E_5 = \mathrm{inf}_{5E}+\mathrm{inf}_{5U} - k_E E_5,\quad
\dot I_5 = k_E E_5 - \delta_I I_5,\quad
\dot V_5 = p I_5 - c V_5
\]

\(\beta\) is derived so that \(R_0 = \beta F_5(R_{s5,0}) T_0 p / (c \delta_I)\) equals Ribeiro 2010 median 8.

## Receptor trafficking

\[
k_{\mathrm{int},x} = k_{\mathrm{int},0,x} + k_{\mathrm{int},C} C_x + k_{\mathrm{int,lig},x} L_x(t)
\]

\[
k_{\mathrm{rec},x} = \frac{k_{\mathrm{rec},0,x}}{1+k_{\mathrm{rec},C} C_x},\quad
k_{\mathrm{syn},x} = \frac{k_{\mathrm{syn},0,x}}{1+k_{\mathrm{syn},C} C_x}
\]

\[
\dot R_{s,x} = k_{\mathrm{syn},x} + k_{\mathrm{rec},x} R_{i,x} - k_{\mathrm{int},x} R_{s,x}
\]

\[
\dot R_{i,x} = k_{\mathrm{int},x} R_{s,x} - k_{\mathrm{rec},x} R_{i,x} - k_{\mathrm{deg},x} R_{i,x}
\]

At the unstimulated steady state,

\[
R_{i,0}=\frac{k_{\mathrm{int},0}R_{s,0}}{k_{\mathrm{rec},0}+k_{\mathrm{deg}}},
\qquad
k_{\mathrm{syn},0}=\frac{k_{\mathrm{int},0}k_{\mathrm{deg}}R_{s,0}}{k_{\mathrm{rec},0}+k_{\mathrm{deg}}}.
\]

These expressions are used directly to initialize the internal pools and derive baseline synthesis.

## Circuit

Qualifying activity \(E_x = \eta\,\mathrm{inf}_{xE} + P_x(t)\). Memory:

\[
\dot A_x = \delta_A (E_x - A_x)
\]

so a constant input \(u\) yields \(A^*=u\). Half-life \(t_{1/2,A} = \ln 2/\delta_A\).

HIV confirmation is driven by the total productively infected-cell state,\n\n\[\nI_{\mathrm{tot}} = I_5 + I_4,\qquad\n\pi_B = \alpha_B H(I_{\mathrm{tot}};K_B,n_B),\qquad\nH(z;K,n)=\frac{z^n}{K^n+z^n}.\n\]\n\nFor \(\tau_B>0\), the confirmation drive passes through a three-stage Erlang delay chain with rate \(k_{\mathrm{del}}=3/\tau_B\):\n\n\[\n\dot Q_1=k_{\mathrm{del}}(\pi_B-Q_1),\qquad\n\dot Q_2=k_{\mathrm{del}}(Q_1-Q_2),\qquad\n\dot Q_3=k_{\mathrm{del}}(Q_2-Q_3),\n\]\n\nand \(\pi_B^{\mathrm{delayed}}=Q_3\). For \(\tau_B=0\), the delay chain is bypassed and \(\pi_B^{\mathrm{delayed}}=\pi_B\). The confirmation state then obeys\n\n\[\n\dot B=\pi_B^{\mathrm{delayed}}-\delta_B B,\qquad\n\delta_B=\frac{\ln 2}{t_{1/2,B}}.\n\]\n\nThus \(\tau_B\) is the mean production/recognition delay, whereas \(t_{1/2,B}\) controls persistence of the accumulated confirmation state. Baseline values are \(\alpha_B=8\), \(K_B=1000\) cells mL\(^{-1}\), \(n_B=2\), \(\tau_B=0.25\) d, and \(t_{1/2,B}=0.5\) d. These are exploratory effective circuit quantities, not measured kinetics of a specific HIV sensor.

Coincidence:

\[
\dot C_x = k_f H(A_x; T_A, n_A)\, H(B; 0.3, 2) - (k_r+\delta_C) C_x
\]

Circuit-off comparator zeros synthetic production. Constitutive CCR5-KO comparator clamps surface CCR5. HIV-triggered comparator sets \(C_5=C_4=B\).
