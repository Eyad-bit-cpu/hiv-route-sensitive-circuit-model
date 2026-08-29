"""Ensure descriptive YAML intents do not silently drift from scenario definitions."""

import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from scenarios import scenario_by_name
from simulate import resolve_params


def _load(rel):
    return yaml.safe_load((ROOT / rel).read_text(encoding="utf-8"))


def test_run0a_config_matches_factory():
    p = resolve_params()
    cfg = _load("configs/run0_validation/run0a.yaml")
    sc = scenario_by_name("run0a_backbone", p)
    assert sc.V5_0 == pytest.approx(cfg["V5_0"])
    assert sc.V4_0 == pytest.approx(cfg["V4_0"])
    assert sc.t_end == pytest.approx(cfg["t_end"])
    assert sc.circuit_on is cfg["circuit_on"]


def test_run1_config_matches_factory():
    p = resolve_params()
    cfg = _load("configs/run1_r5/config.yaml")
    for name in ["run1_r5_control", "run1_r5_circuit"]:
        sc = scenario_by_name(name, p)
        assert sc.V5_0 == pytest.approx(cfg["V5_0"])
        assert sc.V4_0 == pytest.approx(cfg["V4_0"])
        assert sc.t_end == pytest.approx(cfg["t_end"])
        assert sc.P5(0.0) == pytest.approx(cfg["P5"])
        assert sc.P4(0.0) == pytest.approx(cfg["P4"])


def test_run2_config_matches_factory():
    p = resolve_params()
    cfg = _load("configs/run2_x4/config.yaml")
    for name in ["run2_x4_control", "run2_x4_circuit"]:
        sc = scenario_by_name(name, p)
        assert sc.V5_0 == pytest.approx(cfg["V5_0"])
        assert sc.V4_0 == pytest.approx(cfg["V4_0"])
        assert sc.t_end == pytest.approx(cfg["t_end"])


def test_run3_config_matches_factory():
    p = resolve_params()
    cfg = _load("configs/run3_false_route/config.yaml")
    sc = scenario_by_name("run3_false_route", p)
    assert sc.V4_0 == pytest.approx(cfg["V4_0"])
    assert sc.P5(0.0) == pytest.approx(cfg["P5"])
    assert sc.P4(0.0) == pytest.approx(cfg["P4"])
    assert sc.circuit_on is cfg["circuit_on"]
