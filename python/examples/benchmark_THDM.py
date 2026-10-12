"""2HDM benchmark scan over lambda3 through the Python interface.

The Python counterpart of example/Benchmarks/benchmark_THDM.cpp: the same benchmark point, scan and
PhaseTracer settings (helper_functions.hpp). The potential is either
    --potential tl2 (default): TransitionListener v2's own R2HDM, wrapped for PhaseTracer by
                               phasetracer.models.TL2HDM (needs TransitionListener, see
                               phasetracer/models/tl2hdm.py);
    --potential cpp:           the compiled C++ potential of the C++ benchmark,
                               phasetracer.models.THDM (EffectivePotential/include/models/THDM.hpp).
The two agree to rounding, so `cpp` measures the overhead of the Python interface itself against
bin/run_THDM_benchmark, and `tl2` the cost of a potential evaluated in Python.
The phasetracer package must be importable, e.g. PYTHONPATH=build/python for a development build.

    python3 benchmark_THDM.py                                  # scan with the default settings
    python3 benchmark_THDM.py --potential cpp                  # the same with the C++ potential
    python3 benchmark_THDM.py --tol 1e-6 --fratio 0.01 --no-warm
    python3 benchmark_THDM.py --jobs 8                         # 8 points at a time
    python3 benchmark_THDM.py --lambda3 5.80958300958301       # one point, with debug output

TL's potential is evaluated in Python, so a point takes minutes rather than the seconds of the C++
benchmark; --jobs runs the independent points in separate processes.

The scan writes one row per point to example/Benchmarks/data/THDM_TL2_results.csv (tl2) or
THDM_cpp_python_results.csv (cpp):
    status,lambda3,Tc,Tn,Tp,Treh,alpha,betaH,snrLISA,runtime
with the same status flags as the C++ benchmark (0 success, -1 ... -7 the stage that failed) and
the runtime in ms.
"""

import argparse
import contextlib
import multiprocessing
import sys
import time
from pathlib import Path

import phasetracer as pt
from phasetracer.models import THDM

DATA_DIR = Path(__file__).resolve().parents[2] / "example" / "Benchmarks" / "data"

# Background dof that make PhaseTracer's equation of state match TL's radiation term
# -pi^2/90 g_eff T^4 (example/Benchmarks/tools/background_dof.py).
TL_BACKGROUND_DOF = 60.15297440

# benchmark point, lambda3 BP = 5.80958300958301
BENCHMARK = dict(tan_beta=21.3949, m12_sq=1573.171, lambda1=0.28894, lambda2=0.26237,
                 lambda4=-2.17494, lambda5=-2.24170, yukawa_type=1)

# status codes -> the flags of benchmark_THDM.cpp
FLAGS = {
    pt.StatusCode.Success: 0,
    pt.StatusCode.InvalidConfig: -1,
    pt.StatusCode.PhaseFinderFailed: -1,
    pt.StatusCode.NoPhases: -1,
    pt.StatusCode.TransitionFinderFailed: -2,
    pt.StatusCode.ActionCalculatorFailed: -3,
    pt.StatusCode.ThermoFinderFailed: -4,
    pt.StatusCode.NoTransitions: -5,
    pt.StatusCode.NoThermalParameters: -5,
    pt.StatusCode.GravWaveFailed: -6,
    pt.StatusCode.NoSpectra: -7,
}

# --potential -> the file name stem of the outputs
OUTPUT_STEMS = {"tl2": "THDM_TL2", "cpp": "THDM_cpp_python"}

COLUMNS = ["status", "lambda3", "Tc", "Tn", "Tp", "Treh", "alpha", "betaH", "snrLISA", "runtime"]


def transition_filter(transitions):
    """Keep transitions from the origin, with a false phase extending > 1 GeV below TC, and TC < 150 GeV."""
    def valid(t):
        fv_at_origin = abs(t.false_vacuum[0]) < 1e-3 and abs(t.false_vacuum[1]) < 1e-3
        valid_temp = t.TC - t.false_phase.T[0] > 1
        return fv_at_origin and valid_temp and t.TC < 150
    return [t for t in transitions if valid(t)]


def make_config(action_tol=1e-6, fratio_conv=0.01, deformation_npoints=150, action_steps=50, warm_start=True,
                to_print=False):
    """The settings of example/Benchmarks/helper_functions.hpp."""
    config = pt.Config()

    pf = config.phase_finder
    pf.seed = 0
    pf.check_hessian_singular = False
    pf.hessian_singular_rel_tol = 1e-4
    pf.check_vacuum_at_high = False
    pf.hessian_eig_max_rel_change = 0.
    pf.check_midpoint_hessian = False

    ac = config.action
    ac.method = pt.ActionMethod.PathDeformation
    ac.PD_xtol = action_tol
    ac.PD_phitol = action_tol
    ac.PD_fRatioConv = fratio_conv
    ac.PD_deformation_npoints = deformation_npoints

    tf = config.thermo_finder
    tf.transition_filter = transition_filter
    tf.warm_start_chunk_size = 0 if warm_start else 1
    tf.action_spline_evaluations = action_steps
    tf.update_percolation_temperature = True
    tf.action_smoothing_window = 1
    tf.action_smoothing_order = 3
    tf.action_laurent_tail = True
    tf.percolation_target = 1 - 0.28957  # TL convention
    tf.eos_background_dof = TL_BACKGROUND_DOF
    tf.percolation_print_setting = pt.PrintSettings.VERBOSE
    tf.nucleation_print_setting = pt.PrintSettings.MINIMAL

    gw = config.gravwave
    gw.gw_method = pt.GravWaveMethod.FitFormulae
    gw.default_milestone = pt.MilestoneType.PERCOLATION
    gw.vw = 1.0
    gw.SNR_f_max = 1e1
    gw.max_frequency = 1e1
    gw.SNR_f_min = 1e-5
    gw.min_frequency = 1e-5
    gw.num_frequency = 500

    config.pipeline.to_print = to_print
    return config


def make_model(lambda3, potential):
    if potential == "cpp":
        return THDM(lambda3=lambda3, **BENCHMARK)
    from phasetracer.models import TL2HDM  # imports TransitionListener, so only when it is used
    return TL2HDM(lambda3=lambda3, **BENCHMARK)


def write_debug(runner, potential, n_steps=200):
    """Write S(T), the prefactor and Gamma(T) of the first transition, and its GW spectrum."""
    decay_rate = runner.get_thermal_parameters()[0].decay_rate()
    t_min, t_max = decay_rate.t_min, decay_rate.t_max
    stem = OUTPUT_STEMS[potential]
    with open(DATA_DIR / f"{stem}_debug_action.csv", "w") as f:
        f.write("temperature,action,prefactor,decay_rate\n")
        for i in range(n_steps):
            t = t_min + i * (t_max - t_min) / (n_steps - 1)
            f.write(f"{t},{decay_rate.get_action(t)},{decay_rate.get_prefactor(t)},{decay_rate.get_gamma(t)}\n")
    # the Friedmann evolution written by the C++ benchmark is not exposed to Python
    if runner.get_spectra():
        runner.gravwave_calculator().write_spectrum_to_text(0, str(DATA_DIR / f"{stem}_debug_spectrum.csv"))


def get_results(lambda3, config, potential="tl2", debug=False):
    """Run the pipeline at one lambda3; returns the row [status, lambda3, Tc, ..., runtime]."""
    start = time.perf_counter()
    try:
        model = make_model(lambda3, potential)
    except Exception as e:  # e.g. TL rejects the parameters
        print(f"{potential} model initialisation failed at lambda3 = {lambda3}: {e}", file=sys.stderr)
        return [-1, lambda3] + [0.] * 8

    runner = pt.Runner(model, config)
    status = runner.run()
    if debug:
        print(status)
        if runner.get_thermal_parameters():
            write_debug(runner, potential)
    if not status:
        return [FLAGS.get(status.code, -1), lambda3] + [0.] * 8

    tp = runner.get_thermal_parameters()[0]
    n, p = tp.nucleation, tp.percolation
    snr_lisa = runner.get_spectra()[0].SNR[0]
    runtime_ms = round((time.perf_counter() - start) * 1e3)
    return [0, lambda3, tp.TC, n.temperature, p.temperature, p.reheating_temperature, p.alpha_bar, p.betaH_eff,
            snr_lisa, runtime_ms]


# Scan points run in worker processes; the Config holds a Python callable (the transition filter),
# so each worker builds its own.
_worker_config = None
_worker_potential = None


def _init_worker(settings, potential):
    global _worker_config, _worker_potential
    pt.set_num_threads(1)
    pt.set_log_level("fatal")
    _worker_config = make_config(*settings)
    _worker_potential = potential


def _run_point(lambda3):
    return get_results(lambda3, _worker_config, _worker_potential)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--potential", choices=OUTPUT_STEMS, default="tl2",
                        help="TransitionListener's Python potential (tl2) or the compiled C++ one (cpp)")
    parser.add_argument("--lambda3", type=float, help="run a single point with debug output instead of the scan")
    parser.add_argument("--tol", type=float, default=1e-8, help="path deformation xtol and phitol")
    parser.add_argument("--fratio", type=float, default=0.005, help="path deformation fRatioConv")
    parser.add_argument("--n-deformations", type=int, default=150, help="path deformation npoints")
    parser.add_argument("--action-steps", type=int, default=50, help="action spline evaluations")
    parser.add_argument("--no-warm", action="store_true", help="disable warm starts of the action")
    parser.add_argument("--n-points", type=int, default=200, help="number of lambda3 values")
    parser.add_argument("--lambda3-min", type=float, default=5.5)
    parser.add_argument("--lambda3-max", type=float, default=5.875)
    parser.add_argument("--output", type=Path, help="default: data/THDM_TL2_results.csv or THDM_cpp_python_results.csv")
    parser.add_argument("--jobs", type=int, default=1, help="number of points run in parallel processes")
    args = parser.parse_args()

    # A Python potential holds the GIL, so extra OpenMP threads would only wait. The C++ potential
    # keeps the default threads, like bin/run_THDM_benchmark, unless points already run in parallel.
    if args.potential == "tl2" or args.jobs > 1:
        pt.set_num_threads(1)
    output = args.output or DATA_DIR / f"{OUTPUT_STEMS[args.potential]}_results.csv"
    warm = not args.no_warm
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    if args.lambda3 is not None:
        pt.set_log_level("debug")
        config = make_config(args.tol, args.fratio, args.n_deformations, args.action_steps, warm, to_print=True)
        res = get_results(args.lambda3, config, args.potential, debug=True)
        print(f"Result for lambda3 = {args.lambda3}: " + ",".join(str(x) for x in res))
        return 0

    pt.set_log_level("fatal")
    settings = (f"tol={args.tol}, fRatioConv={args.fratio}, nDeformations={args.n_deformations}, "
                f"nPoints={args.action_steps}, warm={int(warm)}, potential={args.potential}")
    print(f"Running scan with settings: {settings}, jobs={args.jobs}")
    config_args = (args.tol, args.fratio, args.n_deformations, args.action_steps, warm)
    lambda3s = [args.lambda3_min + i * (args.lambda3_max - args.lambda3_min) / (args.n_points - 1)
                for i in range(args.n_points)]

    with contextlib.ExitStack() as stack:
        if args.jobs > 1:
            # spawn rather than fork: forking after OpenMP has started can deadlock the child
            pool = stack.enter_context(multiprocessing.get_context("spawn").Pool(
                args.jobs, initializer=_init_worker, initargs=(config_args, args.potential)))
            results = pool.imap(_run_point, lambda3s)  # in order, as they complete
        else:
            config = make_config(*config_args)
            results = (get_results(lambda3, config, args.potential) for lambda3 in lambda3s)

        out = stack.enter_context(open(output, "w"))
        out.write("# " + ",".join(COLUMNS) + "\n")
        out.write(f"# {settings}\n")
        for i, res in enumerate(results):
            out.write(",".join(f"{x:.15g}" for x in res) + "\n")
            out.flush()
            print(f"Ran point {i + 1} with lambda3 = {lambda3s[i]}")
            print("Result: " + ",".join(str(x) for x in res))
    return 0


if __name__ == "__main__":
    sys.exit(main())
