"""Full PhaseTracer pipeline from Python: phases -> transitions -> thermal parameters -> GW spectrum.

    python3 run_runner.py           # print the results
    python3 run_runner.py --plot    # also save the spectrum to spectrum.png (needs matplotlib)
"""

import sys

import phasetracer as pt
from phasetracer.models import TwoDimModel


def main():
    model = TwoDimModel()

    config = pt.Config()
    config.phase_finder.seed = 1
    config.thermo_finder.percolation_target = 0.71
    config.gravwave.min_frequency = 1e-4
    config.gravwave.max_frequency = 1e0
    config.pipeline.to_print = False  # set True to print every stage as the C++ examples do

    runner = pt.Runner(model, config)
    status = runner.run()
    print(status)
    if not status:
        return 1

    for transition in runner.get_transitions():
        print(f"transition at TC = {transition.TC:.3f}: {transition.false_vacuum} -> {transition.true_vacuum}")

    for tps in runner.get_thermal_parameters():
        p = tps.percolation
        print(f"TC = {tps.TC:.3f}: percolation at T = {p.temperature:.3f}, "
              f"alpha = {p.alpha:.4g}, beta/H = {p.betaH:.4g}, vw = {p.vw:.3g}")
        print(f"  bounce action S/T at percolation = {tps.decay_rate().get_action(p.temperature) / p.temperature:.2f}")

    spectrum = runner.get_spectra()[0]
    print(f"GW peak: f = {spectrum.peak_frequency:.3g} Hz, h^2 Omega = {spectrum.peak_amplitude:.3g}, "
          f"SNR LISA = {spectrum.SNR[0]:.3g}, Taiji = {spectrum.SNR[1]:.3g}")

    if "--plot" in sys.argv:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plt.loglog(spectrum.frequency, spectrum.total_amplitude, label="total")
        plt.loglog(spectrum.frequency, spectrum.sound_wave, "--", label="sound waves")
        plt.loglog(spectrum.frequency, spectrum.lisa_noise, ":", label="LISA noise")
        plt.xlabel("f [Hz]")
        plt.ylabel(r"$h^2\Omega_\mathrm{GW}$")
        plt.legend()
        plt.savefig("spectrum.png", dpi=150)
        print("spectrum saved to spectrum.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
