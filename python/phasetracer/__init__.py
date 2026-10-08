"""Python interface to PhaseTracer.

Run the full pipeline for a model with a Runner:

    import phasetracer as pt
    from phasetracer.models import Z2ScalarSingletModel

    config = pt.Config()
    config.pipeline.to_print = False
    runner = pt.Runner(Z2ScalarSingletModel(), config)
    status = runner.run()
    if status:
        for tps in runner.get_thermal_parameters():
            print(tps.TC, tps.percolation.temperature)

Define a model in Python by subclassing pt.Potential (or pt.OneLoopPotential).
"""

from . import _phasetracer
from ._phasetracer import (
    # pipeline
    Runner,
    Config,
    PhaseFinderConfig,
    TransitionFinderConfig,
    ActionCalculatorConfig,
    ThermoFinderConfig,
    GravWaveConfig,
    PipelineConfig,
    RunStatus,
    RunnerError,
    Stage,
    StatusCode,
    # potentials
    Potential,
    OneLoopPotential,
    DaisyMethod,
    # settings enums
    LogLevel,
    ActionMethod,
    NLoptAlgorithm,
    PrintSettings,
    ValidateMethod,
    MilestoneType,
    GravWaveMethod,
    # results
    Point,
    Phase,
    PhaseEnd,
    Transition,
    Message,
    Profile1D,
    ActionResult,
    TransitionMilestone,
    MilestoneStatus,
    NucleationType,
    ThermalProfiles,
    NucleationHistory,
    FluidProfile,
    GravWaveSpectrum,
    # stage objects (obtained from a Runner)
    PhaseFinder,
    ActionCalculator,
    TransitionFinder,
    ThermoFinder,
    ThermalParameterSet,
    FalseVacuumDecayRate,
    EquationOfState,
    GravWaveCalculator,
    # utilities
    set_log_level,
    set_num_threads,
    get_max_threads,
)
from . import models

__all__ = [name for name in dir() if not name.startswith("_")]
