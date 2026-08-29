"""Registered outcome metrics. Weights are not changed after seeing sweeps."""

from __future__ import annotations

from typing import Dict

import numpy as np
from numpy.typing import NDArray

from model import IDX, receptor_infectability


def trapz(t: NDArray, y: NDArray) -> float:
    fn = getattr(np, "trapezoid", None) or np.trapz
    return float(fn(y, t))


def summarize(t: NDArray, y: NDArray, params: dict, rs5_0: float, rs4_0: float) -> Dict[str, float]:
    V = y[:, IDX["V5"]] + y[:, IDX["V4"]]
    I = y[:, IDX["I5"]] + y[:, IDX["I4"]]
    Rs5 = y[:, IDX["Rs5"]]
    Rs4 = y[:, IDX["Rs4"]]
    nR5 = Rs5 / rs5_0
    nR4 = Rs4 / rs4_0
    auc_v = trapz(t, V)
    peak_v = float(np.max(V))
    t_peak = float(t[int(np.argmax(V))]) if len(t) else 0.0
    loss5 = trapz(t, np.maximum(1.0 - nR5, 0.0)) / max(t[-1] - t[0], 1e-12)
    loss4 = trapz(t, np.maximum(1.0 - nR4, 0.0)) / max(t[-1] - t[0], 1e-12)
    loss5_arr = np.maximum(1.0 - nR5, 0.0)
    loss4_arr = np.maximum(1.0 - nR4, 0.0)
    total_loss = loss5_arr + loss4_arr
    # Suppression share is mathematically undefined before any meaningful suppression has
    # occurred (both losses ~0, e.g. early timepoints before the circuit engages). The old
    # 1e-9 epsilon silently resolved this to S=0 rather than "undefined," which is an
    # artifact, not a result -- e.g. at total_loss~1e-11 (both receptors still ~baseline)
    # the epsilon forced S to a specific small number instead of NaN. Threshold chosen well
    # below any epsilon-scale noise floor but well above genuine early partial suppression.
    meaningful = total_loss > 1e-6
    S = np.where(meaningful, loss5_arr / np.where(meaningful, total_loss, 1.0), np.nan)
    Vsum = y[:, IDX["V5"]] + y[:, IDX["V4"]] + 1e-30
    trop_v = y[:, IDX["V5"]] / Vsum
    
    # Early-window metrics (days 0-10)
    early_mask = t <= 10.0
    if np.any(early_mask):
        t_early = t[early_mask]
        V_early = V[early_mask]
        I_early = I[early_mask]
        auc_v_early = trapz(t_early, V_early)
        peak_v_early = float(np.max(V_early))
        t_peak_early = float(t_early[int(np.argmax(V_early))]) if len(t_early) else 0.0
        auc_i_early = trapz(t_early, I_early)
        # Infected cells at day 7 and 10
        i_day7 = float(I[np.argmin(np.abs(t - 7.0))]) if len(t) else 0.0
        i_day10 = float(I[np.argmin(np.abs(t - 10.0))]) if len(t) else 0.0
    else:
        auc_v_early = 0.0
        peak_v_early = 0.0
        t_peak_early = 0.0
        auc_i_early = 0.0
        i_day7 = 0.0
        i_day10 = 0.0
    
    return {
        "auc_V": auc_v,
        "peak_V": peak_v,
        "t_peak": t_peak,
        "auc_I": trapz(t, I),
        "mean_R5_loss": float(loss5),
        "mean_R4_loss": float(loss4),
        "mean_S_ccr5": float(np.nanmean(S)) if np.any(~np.isnan(S)) else float("nan"),  # nanmean:
                                                # S is undefined (NaN) before meaningful suppression
                                                # begins; if a whole trajectory never reaches meaningful
                                                # suppression, mean_S_ccr5 is legitimately undefined too
                                                # (not a bug) -- explicit guard avoids the benign but
                                                # noisy "Mean of empty slice" RuntimeWarning.
        "mean_trop_V5": float(np.mean(trop_v)),
        "max_A5": float(np.max(y[:, IDX["A5"]])),
        "max_A4": float(np.max(y[:, IDX["A4"]])),
        "max_B": float(np.max(y[:, IDX["B"]])),
        "max_C5": float(np.max(y[:, IDX["C5"]])),
        "max_C4": float(np.max(y[:, IDX["C4"]])),
        "min_T": float(np.min(y[:, IDX["T_E"]] + y[:, IDX["T_U"]])),
        "F5_end": receptor_infectability(float(Rs5[-1]), params["K_R5"], params["h_R5"]),
        "F5_start": receptor_infectability(rs5_0, params["K_R5"], params["h_R5"]),
        # Early-window metrics
        "auc_V_early": auc_v_early,
        "peak_V_early": peak_v_early,
        "t_peak_early": t_peak_early,
        "auc_I_early": auc_i_early,
        "I_day7": i_day7,
        "I_day10": i_day10,
    }



def suppression_share_series(loss_target: NDArray, loss_other: NDArray, min_total: float = 1e-6) -> NDArray:
    """Return target suppression share, masking times when total suppression is negligible.

    A share is undefined when both losses are essentially zero. Returning NaN avoids
    epsilon-driven early-time artifacts in mixed-route plots and summaries.
    """
    lt = np.maximum(np.asarray(loss_target, dtype=float), 0.0)
    lo = np.maximum(np.asarray(loss_other, dtype=float), 0.0)
    total = lt + lo
    out = np.full_like(total, np.nan, dtype=float)
    mask = total > float(min_total)
    out[mask] = lt[mask] / total[mask]
    return out

def recovery_time(
    t: NDArray,
    y: NDArray,
    rs5_0: float,
    rs4_0: float,
    t_off: float,
    recovery_threshold: float = 0.9,
    receptor: str = "R5",
) -> float | None:
    """Time from trigger-source removal until the selected receptor recovers.

    Crossing time is linearly interpolated between output samples so the reported
    recovery endpoint is not tied to the plotting/output grid resolution.
    """
    if receptor not in {"R5", "R4", "both"}:
        raise ValueError("receptor must be 'R5', 'R4', or 'both'")

    nR5 = y[:, IDX["Rs5"]] / rs5_0
    nR4 = y[:, IDX["Rs4"]] / rs4_0

    def crossing(values: NDArray) -> float | None:
        if len(t) == 0 or t[-1] < t_off:
            return None
        i0 = int(np.searchsorted(t, t_off, side="left"))
        i0 = min(i0, len(t) - 1)
        # If already at/above threshold at release, recovery time is zero.
        if t[i0] == t_off and values[i0] >= recovery_threshold:
            return 0.0
        for i in range(i0, len(t)):
            if t[i] < t_off or values[i] < recovery_threshold:
                continue
            if i == 0:
                return max(0.0, float(t[i] - t_off))
            t1, t2 = float(t[i - 1]), float(t[i])
            v1, v2 = float(values[i - 1]), float(values[i])
            if v2 == v1:
                tcross = t2
            else:
                frac = (recovery_threshold - v1) / (v2 - v1)
                frac = min(max(frac, 0.0), 1.0)
                tcross = t1 + frac * (t2 - t1)
            return max(0.0, tcross - t_off)
        return None

    if receptor == "R5":
        return crossing(nR5)
    if receptor == "R4":
        return crossing(nR4)
    r5 = crossing(nR5)
    r4 = crossing(nR4)
    if r5 is None or r4 is None:
        return None
    return max(r5, r4)


def efficacy(auc_circuit: float, auc_control: float) -> float:
    if auc_control <= 0:
        return 0.0
    return 1.0 - auc_circuit / auc_control


def classify(E_V: float, F_false: float, t_rec: float | None, criteria: dict) -> str:
    ev_s = criteria["antiviral_efficacy"]["E_V_success"]
    ev_p = criteria["antiviral_efficacy"]["E_V_partial"]
    f_s = criteria["false_suppression"]["F_max_success"]
    f_p = criteria["false_suppression"]["F_max_partial"]
    t_max = criteria["recovery"]["t_recovery_max_days"]
    # If recovery was not measured, cannot classify as full success
    rec_ok = t_rec is not None and t_rec <= t_max
    if E_V >= ev_s and F_false <= f_s and rec_ok:
        return "successful"
    if E_V >= ev_p and F_false <= f_p:
        return "partial"
    return "failed"
