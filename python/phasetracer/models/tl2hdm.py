"""TL2HDM: TransitionListener's 2HDM as a PhaseTracer potential.

TL2HDM wraps the R2HDM class of TransitionListener v2 (models/TL_2HDM.py), so PhaseTracer runs on
exactly the potential that TransitionListener uses: V(phi, T) is TL's Vtot without the
field-independent radiation term, which PhaseTracer accounts for in the equation of state instead
(config.thermo_finder.eos_background_dof, see TL_BACKGROUND_DOF in examples/benchmark_THDM.py).

    from phasetracer.models import TL2HDM
    model = TL2HDM(tan_beta=21.3949, m12_sq=1573.171, lambda1=0.28894, lambda2=0.26237,
                   lambda3=5.8096, lambda4=-2.17494, lambda5=-2.24170, yukawa_type=1)
    runner = pt.Runner(model, config)

TransitionListener is found through the environment variable TRANSITIONLISTENER_ROOT, which
defaults to ~/PhysicsCodes/TransitionListener.

The gradient, Hessian and d(gradient)/dT use the same finite-difference stencils as PhaseTracer's
defaults (step `h`, second or fourth order from `set_h_4`), but evaluate all stencil points in a
single vectorised call to TL's Vtot, rather than one Python call per point.
"""

import contextlib
import importlib.util
import io
import os
import sys
from pathlib import Path

import numpy as np

from .._phasetracer import Potential

TL_ROOT = Path(os.environ.get("TRANSITIONLISTENER_ROOT", Path.home() / "PhysicsCodes" / "TransitionListener"))


def _load_r2hdm():
    """Import R2HDM from TL_ROOT/models/TL_2HDM.py, adding TL's src/ to the path if needed."""
    try:
        import transitionlistener  # noqa: F401
    except ImportError:
        sys.path.insert(0, str(TL_ROOT / "src"))
    path = TL_ROOT / "models" / "TL_2HDM.py"
    if not path.is_file():
        raise ImportError(f"TL2HDM needs TransitionListener, but {path} was not found; set the environment "
                          "variable TRANSITIONLISTENER_ROOT to the root of the TransitionListener repository")
    # loaded by path, so that TL's top-level `models` package does not shadow anything
    spec = importlib.util.spec_from_file_location("tl_models_2hdm", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.R2HDM


R2HDM = _load_r2hdm()

# PhaseTracer's finite-difference stencils (EffectivePotential/include/effectivepotential/potential.hpp):
# (steps in units of h, weights) for first derivatives / mixed second derivatives, and for the
# diagonal of the Hessian.
_XY = {False: ([-1., 1.], [-0.5, 0.5]),
       True: ([-2., -1., 1., 2.], [1 / 12, -8 / 12, 8 / 12, -1 / 12])}
_XX = {False: ([-1., 0., 1.], [1., -2., 1.]),
       True: ([-2., -1., 0., 1., 2.], [-1 / 12, 16 / 12, -30 / 12, 16 / 12, -1 / 12])}


def _gradient_stencil(n, h_4):
    """Offsets (P, n) in units of h and weights (n, P), such that grad = W @ V(phi + h offsets) / h."""
    steps, coeffs = _XY[h_4]
    offsets, weights = [], []
    for i in range(n):
        for s, c in zip(steps, coeffs):
            dx = np.zeros(n)
            dx[i] = s
            offsets.append(dx)
            row = np.zeros(n)
            row[i] = c
            weights.append(row)
    return np.array(offsets), np.array(weights).T


def _hessian_stencil(n, h_4):
    """Offsets (P, n) in units of h and weights (n*n, P), such that H = W @ V(phi + h offsets) / h^2."""
    xx_steps, xx_coeffs = _XX[h_4]
    xy_steps, xy_coeffs = _XY[h_4]
    offsets, weights = [], []

    def add(dx, i, j, c):
        offsets.append(dx)
        row = np.zeros((n, n))
        row[i, j] = row[j, i] = c
        weights.append(row.ravel())

    for i in range(n):
        for s, c in zip(xx_steps, xx_coeffs):
            dx = np.zeros(n)
            dx[i] = s
            add(dx, i, i, c)
    for i in range(n):
        for j in range(i):
            for si, ci in zip(xy_steps, xy_coeffs):
                for sj, cj in zip(xy_steps, xy_coeffs):
                    dx = np.zeros(n)
                    dx[i], dx[j] = si, sj
                    add(dx, i, j, ci * cj)
    return np.array(offsets), np.array(weights).T


def _gradient_dt_stencil(n, h_4):
    """Offsets (P, n + 1) in units of h, the last column in T, and weights (n, P)."""
    steps, coeffs = _XY[h_4]
    offsets, weights = [], []
    for i in range(n):
        for sx, cx in zip(steps, coeffs):
            for st, ct in zip(steps, coeffs):
                d = np.zeros(n + 1)
                d[i], d[n] = sx, st
                offsets.append(d)
                row = np.zeros(n)
                row[i] = cx * ct
                weights.append(row)
    return np.array(offsets), np.array(weights).T


class TL2HDM(Potential):
    """CP-conserving 2HDM of TransitionListener v2, fields (phi1, phi2) in GeV.

    The arguments are those of R2HDM, with m12_sq in GeV^2. The TransitionListener object is
    available as `self.tl`, e.g. `self.tl.derived_parameters` for the masses and counterterms.
    Set `quiet=False` to see R2HDM's output during initialisation.
    """

    def __init__(self, tan_beta, m12_sq, lambda1, lambda2, lambda3, lambda4, lambda5,
                 yukawa_type=2, v=246.21965079413735, quiet=True):
        super().__init__()  # required: initialises the C++ base class
        params = dict(lambda1=lambda1, lambda2=lambda2, lambda3=lambda3, lambda4=lambda4, lambda5=lambda5,
                      m12_sq_GeV2=m12_sq, tan_beta=tan_beta, yukawa_type=yukawa_type, v_GeV=v)
        with contextlib.redirect_stdout(io.StringIO()) if quiet else contextlib.nullcontext():
            self.tl = R2HDM(params)
        self._stencils = {}

    def _stencil(self, kind, h_4):
        key = (kind, h_4)
        if key not in self._stencils:
            build = {"grad": _gradient_stencil, "hess": _hessian_stencil, "grad_dt": _gradient_dt_stencil}[kind]
            self._stencils[key] = build(2, h_4)
        return self._stencils[key]

    def _Vtot(self, phi, T):
        # TL adds -pi^2/90 g_eff T^4 for the radiation bath unless include_radiation=False
        return self.tl.Vtot(phi, T, include_radiation=False)

    # required
    def V(self, phi, T):
        return float(self._Vtot(np.asarray(phi, dtype=float), T))

    def get_n_scalars(self):
        return 2

    # the potential is symmetric under (phi1, phi2) -> (-phi1, -phi2), as in THDM.hpp
    def apply_symmetry(self, phi):
        return [-np.asarray(phi)]

    # derivatives, each a single vectorised evaluation of TL's Vtot
    def dV_dx(self, phi, T):
        h = self.h
        offsets, weights = self._stencil("grad", self.get_h_4())
        return weights @ self._Vtot(np.asarray(phi) + h * offsets, T) / h

    def d2V_dx2(self, phi, T):
        h = self.h
        offsets, weights = self._stencil("hess", self.get_h_4())
        return (weights @ self._Vtot(np.asarray(phi) + h * offsets, T)).reshape(2, 2) / h**2

    def d2V_dxdt(self, phi, T):
        h = self.h
        offsets, weights = self._stencil("grad_dt", self.get_h_4())
        points = np.asarray(phi) + h * offsets[:, :2]
        return weights @ self._Vtot(points, T + h * offsets[:, 2]) / h**2
