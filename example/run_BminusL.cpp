/**
        Generates and outputs the phase structure for minimal B-L (lps =) and its singlet extension (lps !=0).
        This illustrates how to change settings for models with varying scales.
*/

#include "models/BminusL.hpp"
#include <iostream>
#include <string.h>
#include <cmath>
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"
#include "transition_graph_util.hpp"
#include <Eigen/Eigenvalues>
#include "thermal_function.hpp"
#include "phasetracer.hpp"
#include "gravwave_calculator.hpp"

void printPaths(const std::vector<TransitionGraph::Path> &paths) {
  if (paths.size() == 0) {
    std::cout << "Found no paths!" << std::endl
              << std::endl;
    return;
  }

  std::cout << "Found " << paths.size() << " paths:" << std::endl;

  for (int i = 0; i < paths.size(); ++i) {
    std::cout << "Path " << i + 1 << ": " << paths[i] << std::endl;
  }

  std::cout << std::endl;
}

int main(int argc, char *argv[]) {

  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <input-file> <output-folder> [options]" << std::endl;
    return 1;
  }

  bool check_subcrit = false;
  bool allow_phase_oscillation = false;
  bool bDebug = false;
  bool bTrace = false;
  bool bPlot = false;
  bool bNoTransitionPathFinding = false;

  std::vector<std::string> args;

  if (argc > 1) {
    // See https://stackoverflow.com/questions/15344714/convert-command-line-argument-to-string
    // +1 for the starting pointer so we skip over the executable name argv[0].
    args.assign(argv, argv + argc);
  }

  std::string inputFileName = args[1];
  std::string outputFolderName = args[2];
  // Check for additional input configuration settings.
  for (int i = 3; i < argc; ++i) {
    // std::cout << i << ": " << argv[i] << " " << args[i] << std::endl;
    // std::cout << i << ": " << argv[i] << std::endl;
    if (!bDebug && !bTrace && strcmp(argv[i], "-debug") == 0) {
      LOGGER(debug);
      bDebug = true;
      continue;
    }

    if (!bTrace && strcmp(argv[i], "-trace") == 0) {
      LOGGER(trace);
      bTrace = true;
      continue;
    }

    if (!bPlot && strcmp(argv[i], "-plot") == 0) {
      bPlot = true;
      continue;
    }

    if (!check_subcrit && strcmp(argv[i], "-subcrit") == 0) {
      check_subcrit = true;
      continue;
    }

    if (!allow_phase_oscillation && strcmp(argv[i], "-osc") == 0) {
      allow_phase_oscillation = true;
      continue;
    }

    if (!bNoTransitionPathFinding && strcmp(argv[i], "-no_tpf") == 0) {
      bNoTransitionPathFinding = true;
      continue;
    }
  }

  // Set level of screen output
  if (!bDebug && !bTrace) {
    LOGGER(fatal);
  }

  // Construct model
  EffectivePotential::BminusL model(inputFileName);
  /// get mass scale to rescale settings
  const double mass_scale = model.get_vphi();
  const double small_scale = mass_scale * pow(10, -3);

  /// Rescale step size used in derivatives of the potential
  model.set_h(0.001 * small_scale);

  PhaseTracer::PhaseFinder pf(model);
  pf.set_t_high(1.5 * mass_scale); // Essential

  pf.set_upper_bounds({10 * mass_scale});  // Essential
  pf.set_lower_bounds({-10 * mass_scale}); // Essential

  pf.set_x_abs_identical(1 * small_scale); // Essential to scale
  pf.set_seed(10.23);
  pf.set_x_abs_jump(0.5 * small_scale);
  // passed to nlopt:opt
  pf.set_find_min_x_tol_abs(0.0001 * small_scale);
  // This doesn't affect highscale conformal BminusL as longs a its small enough  // For low scales well below EW scale you may need to rescale this also
  pf.set_find_min_min_step(1e-4);

  /// initial mimum step step size for tracing
  pf.set_find_min_trace_abs_step(1. * small_scale);
  /// initial minimum step step size for locating,
  pf.set_find_min_locate_abs_step(1. * small_scale); // Essential

  /// set to the smaller of time_scale*dt_max_rel and dt_max_abs
  //  so if this stays fixed as scale increases can be a problem,
  // dt steps get smaller and smaller relative to T
  // dt_max is set to minimum of dt_max_rel * time_scale and dt_max_abs
  pf.set_dt_max_abs(50 * small_scale); // Essential when mass scale is much higher than EW
  // dt_min is set to maximum of dt_min_rel * time_scale and dt_min_abs
  pf.set_dt_min_abs(1.e-10 * small_scale); // Essential when mass scale is much lower than EW

  /** Tolerance for checking whether Hessian was singular */
  // pf.hessian_singular_rel_tol(1.e-2);
  /** Tolerance for checking solutions to linear algebra */
  // pf.linear_algebra_rel_tol(1.e-3);

  pf.set_check_vacuum_at_high(true);
  pf.set_check_hessian_singular(true);
  pf.set_check_dx_min_dt(true);
  pf.find_phases();

  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(pf);
  tf.set_check_subcritical_transitions(check_subcrit);
  tf.set_assume_only_one_transition(!allow_phase_oscillation);
  tf.find_transitions();
  tf.find_transition_paths(false);

  std::cout << pf;
  std::cout << tf;

  if (!bNoTransitionPathFinding) {
    LOG(debug) << pf;
    LOG(debug) << "Finding transition paths...";
    tf.find_transition_paths(true);
  }

  if (bDebug) {
    std::cout << pf;
    std::cout << tf;

    if (!bNoTransitionPathFinding) {
      printPaths(tf.get_transition_paths());
    }
  }

  PhaseTracer::phase_plotter(tf, outputFolderName, "phase_structure", bPlot);

  return 0;
}
