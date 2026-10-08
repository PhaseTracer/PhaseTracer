"""Effective potentials shipped with PhaseTracer.

    from phasetracer.models import OneDimModel, TwoDimModel, Z2ScalarSingletModel

To write a model of your own, subclass phasetracer.Potential or phasetracer.OneLoopPotential.
"""

from ._phasetracer import models as _models

OneDimModel = _models.OneDimModel
TwoDimModel = _models.TwoDimModel
Z2ScalarSingletModel = _models.Z2ScalarSingletModel

__all__ = ["OneDimModel", "TwoDimModel", "Z2ScalarSingletModel"]
