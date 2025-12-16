#include "models/RSS.hpp"
#include "models/RSS_HT.hpp"
#include "transition_finder.hpp"
#include "phase_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"
#include "potential_plotter.hpp"
#include "potential_line_plotter.hpp"
#include "transition_graph_util.hpp"

#include <iostream>
#include <iomanip>
#include <eigen3/Eigen/Eigenvalues>
#include <sstream>
#include <vector>
#include <iterator>
#include <fstream>
#include <limits>
#include <assert.h>
#include "boost/filesystem.hpp"

// Based on https://stackoverflow.com/questions/236129/how-do-i-iterate-over-the-words-of-a-string
template <typename DataType>
void split(const std::string &inputString, char delimiter, DataType result, int keepCount) {
  std::istringstream inputStringStream(inputString);
  std::string substring;
  int i = 0;

  while (std::getline(inputStringStream, substring, delimiter) && i++ < keepCount) {
    *result++ = std::stod(substring);
  }
}

std::vector<double> split(const std::string &inputString, char delimiter, int keepCount) {
  std::vector<double> elements;

  split(inputString, delimiter, std::back_inserter(elements), keepCount);

  return elements;
}

void printPaths(std::vector<TransitionGraph::Path> paths) {
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

EffectivePotential::RSS &getModel(bool bHighTempExpansion, const std::vector<double> &data) {
  if (bHighTempExpansion) {
    // Need to do it this way so that we don't return the parent cast. This allows the use of the child classes.
    EffectivePotential::RSS_HT &model = *(new EffectivePotential::RSS_HT());
    model.setParameters(data);
    return model;
  } else {
    EffectivePotential::RSS &model = *(new EffectivePotential::RSS());
    model.setParameters(data);
    return model;
  }
}

void getCriticalTemperatureDataForPoint(std::string inputFileName, std::string outputFolderName, bool check_subcrit,
                                        bool allow_phase_oscillation, bool bPlot, bool bDebug, bool bNoTransitionPathFinding, bool check_Hessian,
                                        bool bNoGSResum, bool bUseBoltzmannSuppression, bool bHighTempExpansion, bool bModStep, bool bMergeGaps, double maxT) {
  std::ifstream inputFile(inputFileName);

  if (!inputFile) {
    std::cerr << "Cannot open the file: " << inputFileName << std::endl;
    std::cerr << "Absolute path: " << boost::filesystem::system_complete(inputFileName) << std::endl;
  }

  std::string line;
  std::vector<double> data;

  std::cout.precision(std::numeric_limits<double>::max_digits10);

  int lineIndex = 0;

  int dataValues = bHighTempExpansion ? 10 : 12;

  std::getline(inputFile, line);
  if (line.size() == 0) {
    std::cerr << "Parameter point line is empty!" << std::endl;
    return;
  }

  // Keep only the first 10 or 12 values depending on whether the high-temperature expansion has been employed.
  // These values are the Lagrangian and input parameters. Any following values are probably transition
  // analysis results.
  data = split(line, ' ', dataValues);

  inputFile.close();

  EffectivePotential::RSS &model = getModel(bHighTempExpansion, data);

  model.set_daisy_method(EffectivePotential::DaisyMethod::Parwani);
  model.set_xi(0);

  if (!bHighTempExpansion) {
    model.set_useGSResummation(!bNoGSResum);
  }

  model.set_bUseBoltzmannSuppression(bUseBoltzmannSuppression);

  Eigen::VectorXd origin(2);
  origin << 0.0, 0.0;
  Eigen::VectorXd vev = model.get_EW_VEV();

  auto prevPrecision = std::cout.precision(std::numeric_limits<double>::max_digits10);

  model.set_DEBUG(true);
  std::cout << "Printing values of the potential in PhaseTracer" << std::endl;
  std::cout << "V0(0 , 0)     : " << model.V0(origin) << std::endl;
  std::cout << "V0(vh, vs)    : " << model.V0(vev) << std::endl;
  std::cout << "V(0 , 0 , 0)  : " << model.V(origin, 0.) << std::endl;
  std::cout << "V(vh, vs, 0)  : " << model.V(vev, 0.) << std::endl;
  std::cout << "V(0 , 0 , 100): " << model.V(origin, 100.) << std::endl;
  std::cout << "V(vh, vs, 100): " << model.V(vev, 100.) << std::endl;
  model.set_DEBUG(false);

  // return;

  std::cout.precision(prevPrecision);

  std::vector<PhaseTracer::Transition> transitions;

  double maxTemp = maxT;

  if (maxT < 0) {
    double temperatureScale = model.get_temperature_scale();
    // Start just above the maximum temperature at which there are two phases so we don't have an erroneous high
    // temperature subcritical transition.
    maxTemp = model.getMaxTemp(temperatureScale, temperatureScale * 10, 500, temperatureScale * 0.01,
                               temperatureScale * 0.001) *
              1.01;
  }

  if (bDebug) {
    std::cout << "Max temperature is: " << maxTemp << std::endl;
  }

  PhaseTracer::PhaseFinder pf(model);
  pf.set_t_high(maxTemp);
  pf.set_check_vacuum_at_high(false);
  pf.set_seed(0);

  if (bModStep) {
    pf.set_x_abs_jump(10.);
    pf.set_x_rel_jump(3.e-2);
    pf.set_t_jump_rel(1.e-4);
    pf.set_dt_min_rel(1.e-9);
    pf.set_dt_min_abs(1.e-11);
  }

  // If check_Hessian is false:
  // Fixes cases where there are slight discontinuities in otherwise presumably second-order transitions that are
  // treated as two distinct phases with no transition between them. This causes issues where no path is found from
  // the high temperature phase to the EW VEV due to this discontinuity.
  pf.set_check_hessian_singular(check_Hessian);

  pf.set_check_merge_phase_gaps(bMergeGaps);

  PhaseTracer::TransitionFinder tf(pf);
  tf.set_check_subcritical_transitions(check_subcrit);
  tf.set_assume_only_one_transition(!allow_phase_oscillation);

  pf.find_phases();
  tf.find_transitions();

  if (!bNoTransitionPathFinding) {
    LOG(debug) << pf;
    LOG(debug) << "Finding transition paths...";
    tf.find_transition_paths(model, true);
  }

  if (bDebug) {
    std::cout << pf;
    std::cout << tf;

    if (!bNoTransitionPathFinding) {
      printPaths(tf.get_transition_paths());
    }
  }

  PhaseTracer::phase_plotter(tf, outputFolderName, "phase_structure", bPlot);
}

int main(int argc, char *argv[]) {
  LOGGER(fatal);

  bool check_subcrit = false;
  bool allow_phase_oscillation = false;
  bool bDebug = false;
  bool bTrace = false;
  bool bPlot = false;
  bool bNoTransitionPathFinding = false;
  bool bUseBoltzmannSuppression = false;
  bool check_Hessian = false;
  bool bNoGSResum = false;
  bool bHighTempExpansion = false;
  bool bModStep = false;
  bool bMergeGaps = false;
  double maxT = -1.;

  // Need to have a list of strings so we can use string operations. Char arrays are quite limited, and we
  // can't then handle maxT easily.
  std::vector<std::string> args;

  if (argc > 1) {
    // See https://stackoverflow.com/questions/15344714/convert-command-line-argument-to-string
    args.assign(argv, argv + argc);
  }

  std::string inputFileName = args[1];
  std::string outputFolderName = args[2];

  for (int i = 3; i < argc; ++i) {
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

    if (!bUseBoltzmannSuppression && strcmp(argv[i], "-boltz") == 0) {
      bUseBoltzmannSuppression = true;
      continue;
    }

    if (!check_Hessian && strcmp(argv[i], "-hessian") == 0) {
      check_Hessian = true;
      continue;
    }

    if (!check_Hessian && strcmp(argv[i], "-noGSResum") == 0) {
      bNoGSResum = true;
      continue;
    }

    if (!bHighTempExpansion && strcmp(argv[i], "-ht") == 0) {
      bHighTempExpansion = true;
      continue;
    }

    if (!bModStep && strcmp(argv[i], "-modstep") == 0) {
      bModStep = true;
      continue;
    }

    if (!bMergeGaps && strcmp(argv[i], "-merge") == 0) {
      bMergeGaps = true;
      continue;
    }

    if (maxT < 0 && args[i].compare(0, 6, "-maxT=") == 0) {
      maxT = std::stod(args[i].substr(6, args[i].size() - 6));
      continue;
    }
  }

  // Set level of screen output
  if (!bDebug && !bTrace) {
    LOGGER(fatal);
  }

  getCriticalTemperatureDataForPoint(inputFileName, outputFolderName, check_subcrit, allow_phase_oscillation, bPlot,
                                     bDebug || bTrace, bNoTransitionPathFinding, check_Hessian, bNoGSResum, bUseBoltzmannSuppression,
                                     bHighTempExpansion, bModStep, bMergeGaps, maxT);
}
