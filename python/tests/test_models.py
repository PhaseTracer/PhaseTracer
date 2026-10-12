"""Tests of phasetracer.models: the compiled models and the models written in Python."""

import math
import os
import subprocess
import sys

import pytest

import phasetracer as pt
import phasetracer.models as models

# example/TestThermoFinder parameters (lambda_hs as in model_params.json)
XSM_MSBAR = dict(lambda_hs=1.05, lambda_s=1.0, ms=125.0, Q=100.0, xi=1.0, use_covariant_gauge=False,
                 use_1L_EWSB_in_0L_mass=False, use_Goldstone_resum=True, use_tree_level_tadpole=False)
# example/xSM/run_xSM_HT.cpp and run_xSM_OSlike.cpp
XSM_SMALL = dict(lambda_hs=0.3, lambda_s=0.1)
# example/Benchmarks/benchmark_THDM.cpp
THDM_BENCHMARK = dict(tan_beta=21.3949, m12_sq=1573.171, lambda1=0.28894, lambda2=0.26237,
                      lambda3=5.80958300958301, lambda4=-2.17494, lambda5=-2.24170, yukawa_type=1)
RSS_POINT = dict(muh_sq=-7823.0, mus_sq=-1.0e4, lh=0.13, ls=0.1, c1=0.1, c2=0.1, b3=-10.0, vs=50.0, ms=100.0,
                 theta=0.1)


def quiet_config(stop_after=pt.Stage.TransitionFinder):
    config = pt.Config()
    config.phase_finder.seed = 1
    config.pipeline.stop_after = stop_after
    config.pipeline.to_print = False
    return config


def all_models():
    """One instance of every compiled model, with the field point to evaluate it at."""
    rss = dict(RSS_POINT, muh_sq_tree_EWSB=RSS_POINT["muh_sq"], mus_sq_tree_EWSB=RSS_POINT["mus_sq"])
    return {
        "OneDimModel": (models.OneDimModel(), [10.0]),
        "TwoDimModel": (models.TwoDimModel(), [100.0, 100.0]),
        "Z2ScalarSingletModel": (models.Z2ScalarSingletModel(), [100.0, 10.0]),
        "xSM_MSbar": (models.xSM_MSbar.from_tadpoles(**XSM_MSBAR), [246.0, 0.0]),
        "xSM_OSlike": (models.xSM_OSlike(ms=65.0, **XSM_SMALL), [246.0, 0.0]),
        "xSM_HT": (models.xSM_HT(m_s=65.0, **XSM_SMALL), [246.0, 0.0]),
        "RSS": (models.RSS(**rss), [246.0, 50.0]),
        "RSS_HT": (models.RSS_HT(**RSS_POINT), [246.0, 50.0]),
        "ConformalDarkU1": (models.ConformalDarkU1(g=0.72, y=0.1, v=100.0), [100.0]),
        "SuperCoolModel": (models.SuperCoolModel(kappa=-0.5, mu_sq=-100.0, lambda_=0.1, mu0_sq=1.0), [50.0]),
        "ScalarYukawa3D": (models.ScalarYukawa3D(), [0.5]),
        "THDM": (models.THDM(**THDM_BENCHMARK), [10.0, 200.0]),
    }


@pytest.mark.parametrize("name", list(all_models()))
def test_compiled_models_construct_and_evaluate(name):
    model, phi = all_models()[name]
    assert model.get_n_scalars() == len(phi)
    assert math.isfinite(model.V(phi, 50.0))
    assert name in models.__all__


def test_xsm_msbar_from_tadpoles_finds_singlet_to_higgs_transition():
    model = models.xSM_MSbar.from_tadpoles(**XSM_MSBAR)
    model.set_daisy_method(pt.DaisyMethod.Parwani)
    assert model.iteration_converged
    runner = pt.Runner(model, quiet_config())
    assert runner.run()

    def singlet_to_higgs(t):
        tv, fv = t.true_vacuum, t.false_vacuum
        return abs(tv[0]) > 5 and abs(tv[1]) < 1e-3 and abs(fv[1]) > 5 and abs(fv[0]) < 1e-3

    assert any(singlet_to_higgs(t) for t in runner.get_transitions())


def test_xsm_ht_analytic_tc():
    model = models.xSM_HT(m_s=65.0, **XSM_SMALL)
    runner = pt.Runner(model, quiet_config())
    assert runner.run()
    tcs = [t.TC for t in runner.get_transitions()]
    assert any(tc == pytest.approx(model.get_TC_from_expression(), rel=1e-3) for tc in tcs)


def test_supercool_constructor_argument_order():
    model = models.SuperCoolModel(kappa=-0.5, mu_sq=-100.0, lambda_=0.1, mu0_sq=1.0)
    assert model.kappa == -0.5
    assert model.mu_sq == -100.0
    assert model.lambda_ == 0.1
    assert model.mu0_sq == 1.0


def test_thdm_needs_parameters():
    with pytest.raises(TypeError):
        models.THDM()
    counterterms = models.THDM(**THDM_BENCHMARK).get_counterterms()
    assert counterterms and all(math.isfinite(value) for value in counterterms.values())


def test_rss_set_parameters():
    model = models.RSS()
    rss = dict(RSS_POINT, muh_sq_tree_EWSB=RSS_POINT["muh_sq"], mus_sq_tree_EWSB=RSS_POINT["mus_sq"])
    model.setParameters(**rss)
    assert math.isfinite(model.V([246.0, 50.0], 50.0))
    assert len(model.get_EW_VEV()) == 2


def test_python_models():
    quartic = models.QuarticThermal()
    assert isinstance(quartic, pt.Potential)
    assert quartic.TC() == pytest.approx(100.50378152592121)
    dark_higgs = models.DarkHiggs()
    assert isinstance(dark_higgs, pt.OneLoopPotential)
    assert math.isfinite(dark_higgs.V([200.0], 100.0))


def run_isolated(code, tl_root):
    """Run code in a fresh interpreter with TRANSITIONLISTENER_ROOT pointing at tl_root."""
    env = dict(os.environ, TRANSITIONLISTENER_ROOT=str(tl_root))
    return subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True)


def test_tl2hdm_is_optional(tmp_path):
    # importing phasetracer and its models never needs TransitionListener ...
    result = run_isolated("import phasetracer; from phasetracer.models import *; print(DarkHiggs, xSM_MSbar)", tmp_path)
    assert result.returncode == 0, result.stderr
    # ... and asking for TL2HDM without it gives a clear error
    result = run_isolated("from phasetracer.models import TL2HDM", tmp_path)
    assert result.returncode != 0
    assert "TRANSITIONLISTENER_ROOT" in result.stderr
