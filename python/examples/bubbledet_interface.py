"""Thermal parameters of the conformal dark U(1) model with a BubbleDet prefactor.

The decay rate is Gamma = A exp(-S3/T). By default PhaseTracer uses the analytic
A = T^4 (S3/2piT)^{3/2}; here A = exp(-S1) is the one-loop functional determinant
computed by BubbleDet (https://github.com/bubbledet/bubbledet) on PhaseTracer's bounce.

BubbleDet is called in-process from config.thermo_finder.prefactor_function, which
PhaseTracer invokes once per action knot with (T, S3/T, ActionResult). The action knots
are still computed in parallel; only the BubbleDet calls are serialised by the GIL.

    python3 bubbledet_interface.py

Needs BubbleDet, numpy and scipy (e.g. ~/.venvs/phasetracer_venv).
"""

import math
import sys
import threading
import warnings

import numpy as np
from scipy.interpolate import CubicSpline

from BubbleDet import BubbleConfig, BubbleDeterminant, ParticleConfig

import phasetracer as pt
from phasetracer.models import ConformalDarkU1

# success code - indicates a successful run
SUCCESS_CODE = pt.StatusCode.Success

# valid codes - these are physical and should return some null likelihood
VALID_CODES = [
    pt.StatusCode.NoPhases,
    pt.StatusCode.NoTransitions,
    pt.StatusCode.NoThermalParameters,
    pt.StatusCode.NoSpectra]

# invalid codes - these are due to errors and should be distinguished
INVALID_CODES = [
    pt.StatusCode.InvalidConfig,
    pt.StatusCode.PhaseFinderFailed,
    pt.StatusCode.ActionCalculatorFailed,
    pt.StatusCode.TransitionFinderFailed,
    pt.StatusCode.ThermoFinderFailed,
    pt.StatusCode.GravWaveFailed,
]


def analytic_prefactor(T, S_over_T):
    """PhaseTracer's default prefactor, T^4 (S3/2piT)^{3/2}."""
    return T**4 * (S_over_T / (2. * math.pi))**1.5


class ConformalU1ForBubbleDet:
    """Presents ConformalDarkU1 and a PhaseTracer bounce to BubbleDet.

    BubbleDet's thermal (dim=3) determinant must be fed the 3d-normalised field
    phi3 = phi / sqrt(T) and potential V3 = V / T, so that its bounce action is S3/T --
    the zero-mode Jacobian depends on the field normalisation, and 4d-normalised input
    makes A too large by T^{3/2}. Masses squared are unchanged (V3'' = V'').

    The fluctuating fields are the light (zero Matsubara mode) 3d fields:
      - radial phi: W = V'', carries the translational and negative modes ("Higgs")
      - Goldstone: lambda phi^2 + Pi_s T^2
    Fermions have no zero Matsubara mode and so are absent from the 3d theory.

    The dark photon is left out by default: its loops generate the barrier of this
    conformal model, so they are already in PhaseTracer's V(phi, T), and its mass g phi
    is far above the bubble scale. Including it double counts them and gives
    |S1| >> S3/T. include_gauge=True adds it anyway (BubbleDet's spin 1 counts
    dim-1 transverse polarisations, i.e. Landau gauge):
      - spatial dark photon A'_i: g^2 phi^2
      - temporal dark photon A'_0: g^2 phi^2 + m_D^2
    """

    def __init__(self, model, include_gauge=False, n_grid=200, margin=0.05, massless_threshold=1e-6):
        self.model = model
        self.include_gauge = include_gauge
        self.n_grid = n_grid
        self.margin = margin
        self.massless_threshold = massless_threshold

    def _field(self, phi3, T):
        return np.array([math.sqrt(T) * phi3])

    def profile(self, T, result):
        """Bounce profile in 3d units, with the origin r=0 prepended (BubbleDet needs it)."""
        prof = result.bubble_profile
        R = np.asarray(prof.R, dtype=float)
        Phi = np.asarray(prof.Phi, dtype=float)
        dPhi = np.asarray(prof.dPhi, dtype=float)
        if R[0] > 0.:
            # regularity at the origin: dPhi(0) = 0, Phi(0) = Phi(rmin) to leading order
            R = np.concatenate(([0.], R))
            Phi = np.concatenate(([Phi[0]], Phi))
            dPhi = np.concatenate(([0.], dPhi))
        return R, Phi / math.sqrt(T), dPhi / math.sqrt(T)

    def _sample(self, Phi3, values_of_phi3):
        """Cubic spline of values_of_phi3 on a grid spanning the profile plus a margin."""
        lo, hi = Phi3.min(), Phi3.max()
        pad = self.margin * ((hi - lo) or 1.)
        grid = np.linspace(lo - pad, hi + pad, self.n_grid)
        cs = CubicSpline(grid, [values_of_phi3(x) for x in grid])

        def f(phi3):
            out = cs(phi3)
            return float(out) if np.ndim(out) == 0 else out

        return f

    def bubble_config(self, T, result):
        R, Phi3, dPhi3 = self.profile(T, result)
        phi_meta = Phi3[-1]
        V_meta = self.model.V(self._field(phi_meta, T), T)
        rt = math.sqrt(T)

        V3 = self._sample(Phi3, lambda x: (self.model.V(self._field(x, T), T) - V_meta) / T)
        dV3 = self._sample(Phi3, lambda x: self.model.dV_dx(self._field(x, T), T)[0] / rt)
        d2V3 = self._sample(Phi3, lambda x: self.model.d2V_dx2(self._field(x, T), T)[0, 0])

        return BubbleConfig(V3, phi_meta, R, Phi3, 3, dV=dV3, d2V=d2V3, dPhi=dPhi3)

    def particles(self, T, bubble):
        model = self.model
        xi = model.xi
        Phi3 = np.asarray(bubble.Phi)
        mass_sq = {"goldstone": lambda x: model.get_scalar_debye_sq(self._field(x, T), xi, T)[1]}
        spins = {"goldstone": 0, "dark_photon_spatial": 1, "dark_photon_temporal": 0}
        if self.include_gauge:
            mass_sq["dark_photon_spatial"] = lambda x: model.get_vector_masses_sq(self._field(x, T))[0]
            mass_sq["dark_photon_temporal"] = lambda x: model.get_vector_debye_sq(self._field(x, T), T)[1]

        particles = [ParticleConfig(W_Phi=bubble.d2V, spin=0, dof_internal=1, zero_modes="Higgs")]
        for name, m2 in mass_sq.items():
            W = self._sample(Phi3, m2)
            # a field that stays massless across the wall gives a degenerate determinant
            if np.max(np.abs(W(Phi3))) < self.massless_threshold:
                continue
            particles.append(ParticleConfig(W_Phi=W, spin=spins[name], dof_internal=1, zero_modes="None"))
        return particles

    @staticmethod
    def action(bubble):
        """S3/T from the 3d-normalised profile; must reproduce PhaseTracer's S/T."""
        R = np.asarray(bubble.R)
        integrand = R**2 * (0.5 * np.asarray(bubble.dPhi)**2 + bubble.V(np.asarray(bubble.Phi)))
        return 4. * math.pi * np.trapezoid(integrand, R)


class BubbleDetPrefactor:
    """Callable for config.thermo_finder.prefactor_function returning A = exp(-S1).

    Never raises: if BubbleDet fails at a temperature, the analytic prefactor is
    returned and the failure recorded. If exp(-S1) would overflow a double (near TC,
    where Gamma is negligible anyway) NaN is returned, so PhaseTracer drops the knot.
    Each call is logged in self.log.

    PhaseTracer calls this from several OpenMP threads; the lock keeps one BubbleDet
    evaluation at a time. The GIL alone does not, since numpy releases it, and BubbleDet
    temporarily turns warnings into errors through the process-wide warnings filters.
    """

    MAX_LOG_A = 700.  # exp(709) is the largest finite double

    def __init__(self, adapter, action_rtol=0.05, quiet=True):
        self.adapter = adapter
        self.action_rtol = action_rtol
        self.quiet = quiet
        self.log = []
        self._lock = threading.Lock()

    def __call__(self, T, S_over_T, result):
        with self._lock:
            return self._evaluate(T, S_over_T, result)

    def _evaluate(self, T, S_over_T, result):
        A_analytic = analytic_prefactor(T, S_over_T)
        entry = {"T": T, "S_over_T": S_over_T, "A_analytic": A_analytic, "fallback": False}
        try:
            with warnings.catch_warnings():
                if self.quiet:
                    warnings.simplefilter("ignore")
                bubble = self.adapter.bubble_config(T, result)
                entry["S_over_T_bubbledet"] = self.adapter.action(bubble)
                det = BubbleDeterminant(bubble, self.adapter.particles(T, bubble), thermal=True)
                S1, S1_err = det.findDeterminant()
            if not math.isfinite(S1):
                raise ValueError(f"non-finite S1 = {S1}")
            A = math.exp(-S1) if -S1 < self.MAX_LOG_A else math.nan
            entry.update(S1=S1, S1_err=S1_err, A=A)
        except Exception as exc:  # noqa: BLE001 -- fall back rather than lose the knot
            entry.update(A=A_analytic, fallback=True, error=f"{type(exc).__name__}: {exc}")
        self.log.append(entry)
        return entry["A"]

    def summary(self):
        lines = [f"{'T':>12} {'S/T':>10} {'S/T (BD)':>10} {'S1':>11} {'A/A_analytic':>13}"]
        for e in sorted(self.log, key=lambda e: e["T"]):
            if e["fallback"]:
                lines.append(f"{e['T']:12.6g} {e['S_over_T']:10.4g}  fallback: {e['error']}")
                continue
            if math.isnan(e["A"]):
                lines.append(f"{e['T']:12.6g} {e['S_over_T']:10.4g} {e['S_over_T_bubbledet']:10.4g} "
                             f"{e['S1']:11.4g}  exp(-S1) overflows: knot dropped")
                continue
            flag = ""
            if abs(e["S_over_T_bubbledet"] / e["S_over_T"] - 1.) > self.action_rtol:
                flag = "  <- action mismatch"
            lines.append(f"{e['T']:12.6g} {e['S_over_T']:10.4g} {e['S_over_T_bubbledet']:10.4g} "
                         f"{e['S1']:11.4g} {e['A'] / e['A_analytic']:13.4g}{flag}")
        n_fail = sum(e["fallback"] for e in self.log)
        lines.append(f"{len(self.log)} BubbleDet calls, {n_fail} fell back to the analytic prefactor")
        return "\n".join(lines)


def get_thermal_params(model, config):

    config.pipeline.stop_after = pt.Stage.ThermoFinder

    runner = pt.Runner(model, config)
    status = runner.run()

    output_code = status.code
    print(f"Output code: {output_code}, message: {status.message}")

    if output_code == SUCCESS_CODE:
        milestone = runner.get_thermal_parameters()[0].percolation
        alpha = milestone.alpha
        Tp = milestone.temperature
        Rs = milestone.Rs
        vw = milestone.vw
        g_eff = milestone.g_eff
        return [alpha, Tp, Rs, vw, g_eff]
    elif output_code in VALID_CODES:
        # A valid code, I'll leave you to judge how to handle it
        return None
    elif output_code in INVALID_CODES:
        # An invalid code, I'll leave you to judge how to handle it
        raise RuntimeError(f"Invalid status code: {output_code}")
    else:
        raise RuntimeError(f"Unknown status code: {output_code}")

    
def get_snr(model, config):

    config.pipeline.stop_after = pt.Stage.GravWave

    runner = pt.Runner(model, config)
    status = runner.run()

    output_code = status.code
    print(f"Output code: {output_code}, message: {status.message}")

    if output_code == SUCCESS_CODE:
        spectra = runner.get_spectra()[0]
        snr = spectra.peak_amplitude
        return [snr]
    elif output_code in VALID_CODES:
        # A valid code, I'll leave you to judge how to handle it
        return None
    elif output_code in INVALID_CODES:
        # An invalid code, I'll leave you to judge how to handle it
        raise RuntimeError(f"Invalid status code: {output_code}")
    else:
        raise RuntimeError(f"Unknown status code: {output_code}")

def make_config():

    config = pt.Config()

    config.phase_finder.seed = 1

    config.thermo_finder.percolation_target = 0.71

    # This corresponds to adding an extra term propto g_eff * T^4 to the EoS,
    # used to include any contributions from radiation not included in Jb and Jf.
    # Set accordingly.
    config.thermo_finder.eos_background_dof = 20

    config.action.PD_phitol = 1e-10
    config.action.PD_xtol = 1e-10

    config.gravwave.min_frequency = 1e-10
    config.gravwave.max_frequency = 1e-5
    config.gravwave.num_frequency = 500

    config.pipeline.scale = pt.Scale.MEV
    config.pipeline.to_print = False

    return config

def main(argc, argv):

    g = 0.72
    v = 100
    y = 0.1
    model = ConformalDarkU1(g=g, v=v, y=y)

    results = {}
    for label in ("analytic", "BubbleDet"):
        config = make_config()
        if label == "BubbleDet":
            prefactor = BubbleDetPrefactor(ConformalU1ForBubbleDet(model))
            config.thermo_finder.prefactor_function = prefactor
        print(f"--- {label} prefactor")
        # results[label] = get_thermal_params(model, config)
        results[label] = get_snr(model, config)
        print(results[label])

    # print(prefactor.summary())

    # for label, thermal_params in results.items():
    #     if thermal_params is not None:
    #         alpha, Tp, Rs, vw, g_eff = thermal_params
    #         print(f"{label:>9}: alpha: {alpha}, Tp: {Tp}, Rs: {Rs}, vw: {vw}, g_eff: {g_eff}")

    return 0


if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
