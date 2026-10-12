"""DarkHiggs: a one-loop model written in Python.

A template for a model subclassing phasetracer.OneLoopPotential: provide the tree-level potential V0
and the field-dependent masses with their degrees of freedom, and PhaseTracer adds the
Coleman-Weinberg and thermal corrections.
"""

from .._phasetracer import DaisyMethod, OneLoopPotential


class DarkHiggs(OneLoopPotential):
    """V0 = -mu^2/2 phi^2 + lambda/4 phi^4 with a vector of mass g phi and a fermion of mass y phi / sqrt(2)"""

    def __init__(self, lam=0.05, v=246.0, g=1.0, y=0.5):
        super().__init__()
        self.lam, self.v, self.g, self.y = lam, v, g, y
        self.mu_sq = lam * v**2  # tree-level vev at v
        self.set_renormalization_scale(v)
        self.set_daisy_method(DaisyMethod.NoDaisy)

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
