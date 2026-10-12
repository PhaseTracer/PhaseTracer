"""QuarticThermal: a polynomial thermal potential written in Python, with an analytic TC.

A template for a model subclassing phasetracer.Potential: only V and get_n_scalars are required.
"""

import math

import numpy as np

from .._phasetracer import Potential


class QuarticThermal(Potential):
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
