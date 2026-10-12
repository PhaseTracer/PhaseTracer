"""Running models written in Python through PhaseTracer.

1. QuarticThermal: a polynomial thermal potential (subclass of pt.Potential) with an analytic TC.
2. DarkHiggs: a one-loop model (subclass of pt.OneLoopPotential) of a scalar charged under a dark
   gauge group with a Yukawa fermion. PhaseTracer adds the Coleman-Weinberg and thermal corrections
   to the tree-level potential, and the full pipeline gives a gravitational wave spectrum.

Both are defined in phasetracer/models/ (quartic_thermal.py, dark_higgs.py), which also serve as
templates for writing your own model. Your model can live in your own script or module: subclass
pt.Potential or pt.OneLoopPotential and pass an instance to pt.Runner.

Python potentials are evaluated many thousands of times, so they are slower than C++ models;
supplying analytic derivatives (dV_dx) helps.
"""

import phasetracer as pt
from phasetracer.models import DarkHiggs, QuarticThermal


def main():
    # a Python V holds the GIL, so extra OpenMP threads would only wait
    pt.set_num_threads(1)

    # 1. critical temperature of the polynomial model against the analytic result
    quartic = QuarticThermal()
    config = pt.Config()
    config.phase_finder.t_high = 200.0
    config.pipeline.stop_after = pt.Stage.TransitionFinder
    config.pipeline.to_print = False
    runner = pt.Runner(quartic, config)  # the Runner keeps the model alive
    print("QuarticThermal:", runner.run())
    print(f"  TC = {runner.get_transitions()[0].TC:.5f} (analytic {quartic.TC():.5f})")

    # 2. full pipeline for the one-loop model
    dark_higgs = DarkHiggs()
    config = pt.Config()
    config.phase_finder.seed = 1
    config.phase_finder.t_high = 500.0
    config.pipeline.to_print = False
    runner = pt.Runner(dark_higgs, config)
    status = runner.run()
    print("DarkHiggs:", status)
    if not status:
        return 1
    tps = runner.get_thermal_parameters()[0]
    p = tps.percolation
    spectrum = runner.get_spectra()[0]
    print(f"  TC = {tps.TC:.3f}, percolation at T = {p.temperature:.3f}, alpha = {p.alpha:.4g}, beta/H = {p.betaH:.4g}")
    print(f"  GW peak f = {spectrum.peak_frequency:.3g} Hz, h^2 Omega = {spectrum.peak_amplitude:.3g}, "
          f"SNR LISA = {spectrum.SNR[0]:.3g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
