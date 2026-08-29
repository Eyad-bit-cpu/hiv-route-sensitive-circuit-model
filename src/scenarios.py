"""Physiological forcing and challenge scenarios (Class C)."""

from __future__ import annotations

import hashlib
import math
from typing import Callable

from model import Scenario



def scenario_seed(master_seed: int, scenario_name: str) -> int:
    """Stable per-scenario child seed derived from the registered master seed."""
    digest = hashlib.sha256(f"{int(master_seed)}:{scenario_name}".encode()).digest()
    return int.from_bytes(digest[:4], "little")

def const(value: float) -> Callable[[float], float]:
    return lambda t, v=value: v


def gaussian_pulse(t0: float, sigma: float, amp: float, baseline: float = 0.02) -> Callable[[float], float]:
    def f(t: float) -> float:
        return baseline + amp * math.exp(-((t - t0) ** 2) / (2.0 * sigma**2))

    return f


def repeated_pulses(t0s: list[float], sigma: float, amp: float, baseline: float = 0.02) -> Callable[[float], float]:
    pulses = [gaussian_pulse(t0, sigma, amp, 0.0) for t0 in t0s]

    def f(t: float) -> float:
        return baseline + sum(p(t) for p in pulses)

    return f


def step(t_on: float, t_off: float, amp: float, baseline: float = 0.02) -> Callable[[float], float]:
    def f(t: float) -> float:
        return baseline + (amp if t_on <= t <= t_off else 0.0)

    return f


def ligand_pulse(t_on: float, t_off: float) -> Callable[[float], float]:
    """Unit ligand occupancy during a window (Run 0B/0C)."""

    def f(t: float) -> float:
        return 1.0 if t_on <= t <= t_off else 0.0

    return f


def noisy_bursts(t_on: float, t_off: float, n_bursts: int, amp: float, baseline: float = 0.02, seed: int = 42) -> Callable[[float], float]:
    """Random effective route-associated bursts during a window to test memory discrimination."""
    import random
    rng = random.Random(seed)
    
    # Generate random burst times within the window
    burst_times = sorted([rng.uniform(t_on, t_off) for _ in range(n_bursts)])
    
    def f(t: float) -> float:
        if not (t_on <= t <= t_off):
            return baseline
        # Add contributions from nearby bursts
        total = baseline
        for bt in burst_times:
            # Each burst is a gaussian with width 0.1 days
            total += amp * math.exp(-((t - bt) ** 2) / (2.0 * 0.1**2))
        return total
    
    return f


def step_off(t_on: float, t_off: float, amp: float, baseline: float = 0.02) -> Callable[[float], float]:
    """Step function that turns OFF after t_off to test recovery."""
    def f(t: float) -> float:
        if t_on <= t <= t_off:
            return baseline + amp
        return baseline
    
    return f


def scenario_by_name(name: str, params: dict) -> Scenario:
    P5 = float(params.get("P5", 0.02))
    P4 = float(params.get("P4", 0.02))
    inoc = float(params.get("V_inoculum", 10.0))
    tend = float(params.get("t_end", 40.0)) if "t_end" in params else 40.0

    low = const(0.02)
    if name == "run0a_backbone":
        return Scenario(name, V5_0=inoc, V4_0=0.0, P5=low, P4=low, t_end=40.0, circuit_on=False)
    if name == "run0b_ccr5_trafficking":
        return Scenario(
            name,
            t_end=2.0,
            circuit_on=False,
            P5=low,
            P4=low,
            L5=ligand_pulse(0.05, 0.15),
        )
    if name == "run0c_cxcr4_trafficking":
        return Scenario(
            name,
            t_end=0.5,
            circuit_on=False,
            P5=low,
            P4=low,
            L4=ligand_pulse(0.02, 0.06),
        )
    if name == "run0d_A_memory":
        return Scenario(name, t_end=2.0, circuit_on=True, P5=repeated_pulses([0.1, 0.15, 0.8], 0.01, 1.0, 0.0), P4=const(0.0))
    if name == "run0e_and_gate":
        return Scenario(name, t_end=5.0, circuit_on=True, P5=const(0.0), P4=const(0.0))
    if name == "run1_r5_control":
        return Scenario(name, V5_0=inoc, circuit_on=False, P5=const(P5), P4=const(P4), t_end=tend)
    if name == "run1_r5_circuit":
        return Scenario(name, V5_0=inoc, circuit_on=True, P5=const(P5), P4=const(P4), t_end=tend)
    if name == "run2_x4_control":
        return Scenario(name, V4_0=inoc, circuit_on=False, P5=const(P5), P4=const(P4), t_end=tend)
    if name == "run2_x4_circuit":
        return Scenario(name, V4_0=inoc, circuit_on=True, P5=const(P5), P4=const(P4), t_end=tend)
    if name == "run3_false_route":
        return Scenario(
            name,
            V4_0=inoc,
            circuit_on=True,
            P5=const(1.2),
            P4=const(0.02),
            t_end=tend,
        )
    if name == "run4_isolated_pulse":
        return Scenario(name, V5_0=inoc, P5=gaussian_pulse(5.0, 0.3, 1.5, P5), P4=const(P4), t_end=tend)
    if name == "run4_repeated_pulses":
        return Scenario(
            name,
            V5_0=inoc,
            P5=repeated_pulses([4.0, 5.0, 6.0, 12.0], 0.15, 1.2, P5),
            P4=const(P4),
            t_end=tend,
        )
    if name == "run4_sustained":
        return Scenario(name, V5_0=inoc, P5=step(3.0, 20.0, 1.0, P5), P4=const(P4), t_end=tend)
    if name == "run5_mix_90_10":
        return Scenario(name, V5_0=0.9 * inoc, V4_0=0.1 * inoc, t_end=tend)
    if name == "run5_mix_50_50":
        return Scenario(name, V5_0=0.5 * inoc, V4_0=0.5 * inoc, t_end=tend)
    if name == "run5_mix_10_90":
        return Scenario(name, V5_0=0.1 * inoc, V4_0=0.9 * inoc, t_end=tend)
    if name == "run6_secondary_x4_emergence":
        return Scenario(name, V5_0=inoc, V4_0=0.0, secondary_x4_day=10.0, secondary_x4_V4=inoc, t_end=tend)
    # Backwards-compatible alias for pre-audit naming.
    if name == "run6_tropism_shift":
        return scenario_by_name("run6_secondary_x4_emergence", params)
    if name.startswith("run7_tauB"):
        return Scenario(name, V5_0=inoc, t_end=tend)
    if name == "run8_memory_test":
        # Use noisy physiological bursts to test memory discrimination
        return Scenario(
            name,
            V4_0=inoc,
            circuit_on=True,
            P5=noisy_bursts(
                0.0, 40.0, n_bursts=8, amp=1.2, baseline=0.02,
                seed=scenario_seed(int(params.get("master_seed", 20260827)), name),
            ),
            P4=const(0.02),
            t_end=tend,
        )
    if name == "run_recovery_test":
        # Clean release-phase recovery protocol. Phase 1 establishes suppression
        # under R5 challenge. At day 14 the scripted route input is removed AND
        # the active infection source is cleared by the simulator, while A/B/C,
        # delay-chain, and receptor states are retained and allowed to relax.
        # This tests trigger removal -> circuit decay -> receptor recovery, rather
        # than asking recovery to occur during continuing self-sustaining infection.
        return Scenario(
            name,
            V5_0=inoc,
            circuit_on=True,
            P5=step_off(0.0, 14.0, 1.0, 0.0),
            P4=const(0.0),
            t_end=35.0,
            trigger_removal_day=14.0,
            clear_infection_on_trigger_removal=True,
        )
    if name == "run9_fast_only":
        # Fast trafficking arm ONLY: internalization active; both sustained-arm
        # effects (synthesis suppression and recycling suppression) are zeroed.
        return Scenario(name, V5_0=inoc, k_synC_scale=0.0, k_recC_scale=0.0, t_end=tend)
    if name == "run9_sustained_only":
        # Sustained arm ONLY: synthesis suppression and recycling suppression are
        # active, while fast internalization is zeroed (Sec 2.4: the
        # sustained arm = synthesis + recycling suppression) -- so this scenario's
        # BEHAVIOR was already conceptually closer to correct than its former name/
        # its sibling scenario suggested; renamed for clarity and to make the true
        # fast-only ablation (above) exist as its proper counterpart.
        return Scenario(name, V5_0=inoc, k_intC_scale=0.0, t_end=tend)
    if name == "run9_both":
        return Scenario(name, V5_0=inoc, t_end=tend)
    # Back-compatibility aliases for the old (mislabeled) run9 scenario names --
    # kept so existing configs/scripts referencing them don't break, but they now
    # point at the SAME corrected definitions above rather than silently keeping
    # the old leaky-recycling-arm behavior under the old names.
    if name == "run9_int_only":
        return scenario_by_name("run9_fast_only", params)
    if name == "run9_syn_only":
        return scenario_by_name("run9_sustained_only", params)
    if name == "comparator_ccr5_ko":
        return Scenario(name, V5_0=inoc, circuit_on=False, constitutive_R5=0.05, t_end=tend)
    if name == "comparator_hiv_triggered":
        return Scenario(name, V5_0=inoc, hiv_triggered_both=True, t_end=tend)
    raise KeyError(name)
