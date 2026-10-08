"""Tests of the Python interface: run with  python3 -m pytest python/tests"""

import copy
import math

import numpy as np
import pytest

import phasetracer as pt
from phasetracer.models import OneDimModel, TwoDimModel, Z2ScalarSingletModel


def precise_config(stop_after=pt.Stage.TransitionFinder):
    config = pt.Config()
    config.phase_finder.seed = 1
    config.phase_finder.find_min_x_tol_rel = 1e-8
    config.phase_finder.find_min_x_tol_abs = 1e-8
    config.transition_finder.TC_tol_rel = 1e-16
    config.pipeline.stop_after = stop_after
    config.pipeline.to_print = False
    return config


class QuarticThermal(pt.Potential):
    """V = D (T^2 - T0^2) phi^2 - E T phi^3 + lambda/4 phi^4, with analytic TC."""

    def __init__(self, D=0.1, E=0.01, lam=0.1, T0=100.0):
        super().__init__()
        self.D, self.E, self.lam, self.T0 = D, E, lam, T0

    def V(self, phi, T):
        x = phi[0]
        return self.D * (T**2 - self.T0**2) * x**2 - self.E * T * x**3 + 0.25 * self.lam * x**4

    def get_n_scalars(self):
        return 1

    def dV_dx(self, phi, T):
        x = phi[0]
        return np.array([2 * self.D * (T**2 - self.T0**2) * x - 3 * self.E * T * x**2 + self.lam * x**3])

    def forbidden(self, phi):
        return phi[0] < -0.1

    def TC(self):
        return self.T0 / math.sqrt(1 - self.E**2 / (self.lam * self.D))


def test_one_dim_model_up_to_transition_finder():
    model = OneDimModel()
    runner = pt.Runner(model, precise_config())
    assert runner.run_id == 0

    status = runner.run()
    assert status
    assert status.code == pt.StatusCode.Success
    assert runner.run_id == 1

    assert len(runner.get_phases()) > 0
    transitions = runner.get_transitions()
    assert len(transitions) == 1
    assert transitions[0].TC == pytest.approx(model.get_TC_from_expression(), rel=1e-8)

    assert runner.has(pt.Stage.TransitionFinder)
    assert not runner.has(pt.Stage.ThermoFinder)
    with pytest.raises(RuntimeError):
        runner.thermo_finder()


def test_config_edits_and_copies():
    config = pt.Config()
    config.phase_finder.seed = 3
    config.phase_finder.upper_bounds = [500.0]
    assert config.phase_finder.seed == 3
    assert config.phase_finder.upper_bounds == [500.0]

    runner = pt.Runner(OneDimModel(), config)
    runner.config.phase_finder.t_high = 321.0
    assert runner.config.phase_finder.t_high == 321.0

    clone = copy.deepcopy(config)
    clone.phase_finder.seed = 7
    assert config.phase_finder.seed == 3

    assert config.validate()
    config.gravwave.min_frequency, config.gravwave.max_frequency = 1.0, 1e-3
    status = config.validate()
    assert not status
    assert status.code == pt.StatusCode.InvalidConfig
    assert "min_frequency" in status.message


def test_python_potential():
    model = QuarticThermal()
    config = precise_config()
    config.phase_finder.t_high = 200.0
    pt.set_num_threads(1)
    runner = pt.Runner(model, config)
    del model  # the Runner keeps the model alive
    status = runner.run()
    assert status, str(status)
    transitions = runner.get_transitions()
    assert len(transitions) == 1
    assert transitions[0].TC == pytest.approx(QuarticThermal().TC(), rel=1e-6)


def test_filter_rejecting_everything():
    config = precise_config(stop_after=pt.Stage.GravWave)
    config.thermo_finder.transition_filter = lambda transitions: []
    runner = pt.Runner(OneDimModel(), config)
    status = runner.run()
    assert status.code == pt.StatusCode.NoThermalParameters
    assert status.stage == pt.Stage.ThermoFinder
    assert runner.get_thermal_parameters() == []


def test_throw_on_error():
    config = precise_config()
    config.gravwave.min_frequency, config.gravwave.max_frequency = 1.0, 1e-3
    config.pipeline.throw_on_error = True
    runner = pt.Runner(OneDimModel(), config)
    with pytest.raises(pt.RunnerError) as error:
        runner.run()
    assert error.value.status.code == pt.StatusCode.InvalidConfig
    assert error.value.status.stage == pt.Stage.Config


def test_stale_handles():
    runner = pt.Runner(OneDimModel(), precise_config())
    runner.run()
    phase_finder = runner.phase_finder()
    phases = runner.get_phases()
    assert len(phase_finder.get_phases()) == len(phases)

    runner.run()
    with pytest.raises(RuntimeError, match="stale"):
        phase_finder.get_phases()
    assert len(phases) > 0  # copies stay valid
    assert len(runner.phase_finder().get_phases()) == len(phases)


def test_handle_keeps_runner_alive():
    runner = pt.Runner(OneDimModel(), precise_config())
    runner.run()
    transition_finder = runner.transition_finder()
    del runner
    assert len(transition_finder.get_transitions()) == 1


def test_full_pipeline_two_dim_model():
    config = pt.Config()
    config.phase_finder.seed = 1
    config.pipeline.to_print = False
    runner = pt.Runner(TwoDimModel(), config)
    status = runner.run()
    assert status, str(status)

    sets = runner.get_thermal_parameters()
    assert len(sets) >= 1
    percolation = sets[0].percolation
    assert percolation.status == pt.MilestoneStatus.YES
    assert 0 < percolation.temperature < sets[0].TC
    assert math.isfinite(sets[0].decay_rate().get_action(percolation.temperature))
    energy_false, energy_true = sets[0].equation_of_state().get_energy(percolation.temperature)
    assert energy_false > energy_true

    spectra = runner.get_spectra()
    assert len(spectra) >= 1
    assert len(spectra[0].frequency) == config.gravwave.num_frequency
    assert len(spectra[0].SNR) == 2

    snr = runner.gravwave_calculator().get_SNR(percolation)
    assert len(snr) == 2


def test_z2_model_reports_thermal_failure():
    """The Z2 high-temperature model has a transition, but too weak for the thermal stage."""
    model = Z2ScalarSingletModel()
    config = pt.Config()
    config.pipeline.to_print = False
    runner = pt.Runner(model, config)
    status = runner.run()
    assert runner.get_transitions()[0].TC == pytest.approx(model.get_TC_from_expression(), rel=1e-3)
    if not status:
        assert status.code == pt.StatusCode.NoThermalParameters
        assert status.warnings


def test_log_level():
    pt.set_log_level("warning")
    pt.set_log_level(pt.LogLevel.fatal)
    with pytest.raises(ValueError):
        pt.set_log_level("loud")


def test_python_prefactor_function():
    calls = []

    def prefactor(T, S_over_T, result):
        calls.append(T)
        return T**4

    config = pt.Config()
    config.phase_finder.seed = 1
    config.pipeline.to_print = False
    config.pipeline.stop_after = pt.Stage.ThermoFinder
    config.thermo_finder.prefactor_function = prefactor
    runner = pt.Runner(TwoDimModel(), config)
    assert runner.run()
    assert calls
    decay_rate = runner.get_thermal_parameters()[0].decay_rate()
    T = 0.5 * (decay_rate.t_min + decay_rate.t_max)
    assert decay_rate.get_prefactor(T) == pytest.approx(T**4, rel=1e-6)
