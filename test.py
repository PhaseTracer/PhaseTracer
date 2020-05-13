import cppyy
import os

# Load PhaseTracer library and relevant headers
dir_ = os.path.dirname(__file__)
cppyy.load_library(os.path.join(dir_, './lib/libphasetracer.so'))
cppyy.add_include_path(os.path.join(dir_, './include/'))
cppyy.add_include_path(os.path.join(dir_, './EffectivePotential/include/effectivepotential/'))
cppyy.include(os.path.join(dir_, './include/phase_finder.hpp'))
cppyy.include(os.path.join(dir_, './include/transition_finder.hpp'))


# Load some example models from EffectivePotential
cppyy.load_library(os.path.join(dir_, './EffectivePotential/lib/libeffectivepotential.so'))
cppyy.include(os.path.join(dir_, './EffectivePotential/include/models/1D_test_model.hpp'))
cppyy.include(os.path.join(dir_, './EffectivePotential/include/models/2D_test_model.hpp'))

# Load examples from effective potential

from cppyy.gbl import PhaseTracer, EffectivePotential

# model = EffectivePotential.OneDimModel()
model = EffectivePotential.TwoDimModel()
pf = PhaseTracer.PhaseFinder(model)
pf.find_phases()  # This actually works!



