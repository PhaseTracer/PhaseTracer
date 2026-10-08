"""Defining new models in Python and running them through PhaseTracer.

1. QuarticThermal: a polynomial thermal potential (subclass of pt.Potential) with an analytic TC.
2. DarkHiggs: a one-loop model (subclass of pt.OneLoopPotential) of a scalar charged under a dark
   gauge group with a Yukawa fermion. PhaseTracer adds the Coleman-Weinberg and thermal corrections
   to the tree-level potential, and the full pipeline gives a gravitational wave spectrum.

Python potentials are evaluated many thousands of times, so they are slower than C++ models;
supplying analytic derivatives (dV_dx) helps.
"""

import math

import numpy as np

import phasetracer as pt


class QuarticThermal(pt.Potential):
    """V = D (T^2 - T0^2) phi^2 - E T phi^3 + lambda/4 phi^4"""

    def __init__(self, D=0.1, E=0.01, lam=0.1, T0=100.0):
        super().__init__()  # required: initialises the C++ base class
        self.D, self.E, self.lam, self.T0 = D, E, lam, T0

    # required
    def V(self, phi, T):  # phi is a numpy array of length get_n_scalars()
        x = phi[0]
        return self.D * (T**2 - self.T0**2) * x**2 - self.E * T * x**3 + 0.25 * self.lam * x**4

    def get_n_scalars(self):
        return 1

    # optional: analytic gradient (the default is numerical) and a forbidden region
    def dV_dx(self, phi, T):
        x = phi[0]
        return np.array([2 * self.D * (T**2 - self.T0**2) * x - 3 * self.E * T * x**2 + self.lam * x**3])

    def forbidden(self, phi):
        return phi[0] < -0.1

    def TC(self):
        """Degenerate minima: E^2 T^2 = lambda D (T^2 - T0^2)"""
        return self.T0 / math.sqrt(1 - self.E**2 / (self.lam * self.D))


class DarkHiggs(pt.OneLoopPotential):
    """V0 = -mu^2/2 phi^2 + lambda/4 phi^4 with a vector of mass g phi and a fermion of mass y phi / sqrt(2)"""

    def __init__(self, lam=0.05, v=246.0, g=1.0, y=0.5):
        super().__init__()
        self.lam, self.v, self.g, self.y = lam, v, g, y
        self.mu_sq = lam * v**2  # tree-level vev at v
        self.set_renormalization_scale(v)
        self.set_daisy_method(pt.DaisyMethod.NoDaisy)

    def V0(self, phi):
        x = phi[0]
        return -0.5 * self.mu_sq * x**2 + 0.25 * self.lam * x**4

    def get_n_scalars(self):
        return 1

    def apply_symmetry(self, phi):  # phi -> -phi
        return [-phi]

    def get_scalar_masses_sq(self, phi, xi):
        return [-self.mu_sq + 3 * self.lam * phi[0] ** 2]

    def get_scalar_dofs(self):
        return [1.0]

    def get_vector_masses_sq(self, phi):
        return [self.g**2 * phi[0] ** 2]

    def get_vector_dofs(self):
        return [3.0]

    def get_fermion_masses_sq(self, phi):
        return [0.5 * self.y**2 * phi[0] ** 2]

    def get_fermion_dofs(self):  # positive; PhaseTracer applies the fermionic sign
        return [4.0]


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
