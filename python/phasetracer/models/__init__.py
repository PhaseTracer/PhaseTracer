"""Effective potentials shipped with PhaseTracer.

    from phasetracer.models import TwoDimModel, xSM_MSbar, DarkHiggs

Compiled (C++) models:
    test models     OneDimModel, TwoDimModel, Z2ScalarSingletModel
    xSM             xSM_MSbar, xSM_OSlike, xSM_HT (and their base xSM_base)
    RSS             RSS, RSS_HT
    other BSM       ConformalDarkU1, SuperCoolModel, ScalarYukawa3D
    2HDM            THDM

Models written in Python (templates for your own, see quartic_thermal.py and dark_higgs.py):
    QuarticThermal  polynomial thermal potential with an analytic TC
    DarkHiggs       one-loop dark Higgs with a dark gauge boson and a Yukawa fermion
    TL2HDM          TransitionListener's 2HDM; needs a TransitionListener checkout
                    (TRANSITIONLISTENER_ROOT) and is only imported when first used

Your own models do not need to live here: subclass phasetracer.Potential or
phasetracer.OneLoopPotential anywhere and pass an instance to phasetracer.Runner.
"""

from .._phasetracer import models as _models
from .dark_higgs import DarkHiggs
from .quartic_thermal import QuarticThermal

# compiled models
OneDimModel = _models.OneDimModel
TwoDimModel = _models.TwoDimModel
Z2ScalarSingletModel = _models.Z2ScalarSingletModel
xSM_base = _models.xSM_base
xSM_MSbar = _models.xSM_MSbar
xSM_OSlike = _models.xSM_OSlike
xSM_HT = _models.xSM_HT
RSS = _models.RSS
RSS_HT = _models.RSS_HT
ConformalDarkU1 = _models.ConformalDarkU1
SuperCoolModel = _models.SuperCoolModel
ScalarYukawa3D = _models.ScalarYukawa3D
THDM = _models.THDM

# TL2HDM is not listed, so that `from phasetracer.models import *` does not need TransitionListener
__all__ = [
    "OneDimModel", "TwoDimModel", "Z2ScalarSingletModel",
    "xSM_base", "xSM_MSbar", "xSM_OSlike", "xSM_HT",
    "RSS", "RSS_HT",
    "ConformalDarkU1", "SuperCoolModel", "ScalarYukawa3D",
    "THDM",
    "QuarticThermal", "DarkHiggs",
]


def __getattr__(name):
    # models with optional dependencies, imported on first use
    if name == "TL2HDM":
        from .tl2hdm import TL2HDM
        return TL2HDM
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
