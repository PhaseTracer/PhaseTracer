import phasetracer as pt
from phasetracer.models import TwoDimModel

config = pt.Config()
config.pipeline.to_print = True

runner = pt.Runner(TwoDimModel(), config)
status = runner.run()

print(status)

phases = runner.get_phases()
tps = runner.get_thermal_parameters()