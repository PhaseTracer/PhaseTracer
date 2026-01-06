/**
  1D example program for PhaseTracer.
*/

#include <iostream>

#include "models/SubcriticalTransitionExample.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"
#include "transition_graph_util.hpp"

#include <eigen3/Eigen/Eigenvalues>

// Copied from transition_finder.cpp (called changed()), because it's a private method in TransitionFinder.
/*std::vector<bool> hasVEVChanged(const Eigen::VectorXd& true_vacuum, const Eigen::VectorXd&
        false_vacuum, const PhaseTracer::TransitionFinder& tf)
{
        auto delta = true_vacuum - false_vacuum;
        std::vector<bool> changed;

        for (unsigned int i = 0; i < delta.size(); ++i)
        {
                changed.push_back(std::abs(delta[i]) > tf.get_change_abs_tol()
                        + tf.get_change_rel_tol()*std::max(true_vacuum.norm(), false_vacuum.norm()));
        }

        return changed;
}*/

// Transitions isn't const because we modify it.
/*void appendSubcriticalTransitions(const std::vector<PhaseTracer::Phase>& phases, std::vector<PhaseTracer::Transition>&
        transitions, const PhaseTracer::PhaseFinder& pf, const PhaseTracer::TransitionFinder& tf)
{
        // Loop through phases and check if any phase has lower energy than any other phase when it first appears (i.e. at
        // its Tmax).
        for(int i = 0; i < phases.size()-1; ++i)
        {
                double Tmax = phases[i].T.back();
                double energyAtTmax = phases[i].V.back();

                for(int j = j+1; j < phases.size(); ++j)
                {
                        PhaseTracer::Point phasejAtTmax = pf.phase_at_T(phases[j], Tmax);

                        // If the phases do not overlap in temperature, there cannot be a transition between them.
                        if(abs(phasejAtTmax.t - Tmax) > 1)
                        {
                                continue;
                        }

                        // If this second phase has a higher energy, then we have a subcritical transition.
                        // We flag subcritical transitions with a negative critical temperature and gamma.
                        if(phasejAtTmax.potential > energyAtTmax)
                        {
                                Eigen::VectorXd trueVacuum = phases[i].X.back();
                                Eigen::VectorXd falseVacuum = phasejAtTmax.x;
                                double gamma = (trueVacuum - falseVacuum).norm() / Tmax;
                                std::vector<bool> changed = hasVEVChanged(trueVacuum, falseVacuum, tf);
                                double deltaPotential = energyAtTmax - phasejAtTmax.potential;

                                std::cout << "Found subcritical transition from phase " << j << " to phase " << i << " at T=" << Tmax << std::endl;

                                transitions.push_back({PhaseTracer::SUCCESS, -Tmax, phases[i], phases[j], trueVacuum, falseVacuum,
                                        -gamma, changed, deltaPotential, transitions.size()});
                        }
                }
        }
}*/

/*void findTransitionPaths(const PhaseTracer::PhaseFinder& pf, const PhaseTracer::TransitionFinder& tf, bool
        knownHighTPhase)
{
        std::vector<PhaseTracer::Phase> symmetrisedPhases;
        std::vector<PhaseTracer::Transition> symmetrisedTransitions;

        TransitionGraph::extractExplicitSymmetricPhasesAndTransitions(pf.get_phases(), tf.get_transitions(), {0},
                symmetrisedPhases, symmetrisedTransitions);

        Eigen::VectorXd EWVEV(1);
        EWVEV << 1.5;

        //TransitionGraph::PhaseStructureData phaseStructureData = TransitionGraph::extractPhaseStructureData(
        //	pf.get_phases(), tf.get_transitions(), EWVEV);
        TransitionGraph::PhaseStructureData phaseStructureData = TransitionGraph::extractPhaseStructureData(
                symmetrisedPhases, symmetrisedTransitions, EWVEV, pf.get_t_high(), knownHighTPhase);

        if(phaseStructureData.validAtZeroT)
        {
                std::cout << "Phase structure is valid at T=0." << std::endl;
                std::cout << "Phase " << phaseStructureData.EWVEVIndex << " is the EW VEV phase at T=0." << std::endl;
                std::cout << "Phases ";

                for(int i = 0; i < phaseStructureData.highTPhaseIndices.size(); ++i)
                {
                        std::cout << phaseStructureData.highTPhaseIndices[i] << " ";
                }

                if(phaseStructureData.highTPhaseIndices.size() == 1)
                {
                        std::cout << "is the high temperature phase at T=";
                }
                else
                {
                        std::cout << "are the high temperature phases at T=";
                }

                std::cout << pf.get_t_high() << "." << std::endl;

                std::cout << std::endl;
                std::cout << "Finding all transition paths from the high temperature phase to the EW VEV..." << std::endl;

                std::vector<std::vector<TransitionGraph::Path>> paths = TransitionGraph::getTransitionPaths(symmetrisedPhases,
                        symmetrisedTransitions, phaseStructureData);

                int numPaths = 0;

                for(int i = 0; i < paths.size(); ++i)
                {
                        for(int j = 0; j < paths[i].size(); ++j)
                        {
                                numPaths += 1;
                        }
                }

                std::cout << "Found " << numPaths << " paths:" << std::endl;

                int pathID = 1;

                for(int i = 0; i < paths.size(); ++i)
                {
                        for(int j = 0; j < paths[i].size(); ++j)
                        {
                                bool valid = paths[i][j].phases.back() == phaseStructureData.EWVEVIndex;
                                std::cout << (valid ? "[Valid]  " : "[Invalid]") << " Path " << ++pathID  << ": " << paths[i][j]
                                        << std::endl;
                                //std::cout << paths[i] << std::endl;
                        }
                }

                std::cout << std::endl;
        }
        else
        {
                std::cout << "Phase structure is invalid at T=0. It does not describe our Universe." << std::endl;
        }
}*/

/*void printPaths(std::vector<std::vector<TransitionGraph::Path>> paths)
{
        int numPaths = 0;

        for(int i = 0; i < paths.size(); ++i)
        {
                for(int j = 0; j < paths[i].size(); ++j)
                {
                        numPaths += 1;
                }
        }

        std::cout << "Found " << numPaths << " paths:" << std::endl;

        int pathID = 1;
        //bool constrainedLowTPhases = model.get_low_t_phases().size() > 0;

        for(int i = 0; i < paths.size(); ++i)
        {
                for(int j = 0; j < paths[i].size(); ++j)
                {
                        //std::cout << (valid ? "[Valid]  " : "[Invalid]") << " Path " << pathID++  << ": " << paths[i][j]
                        //	<< std::endl;
                        std::cout << "Path " << pathID++ << ": " << paths[i][j] << std::endl;
                        //std::cout << paths[i] << std::endl;
                }
        }

        std::cout << std::endl;
}*/

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

int main(int argc, char *argv[]) {
  bool check_subcrit = false;
  bool allow_phase_oscillation = false;
  bool bDebug = false;
  bool bTrace = false;
  bool bPlot = false;
  bool bNoTransitionPathFinding = false;
  bool bUseBoltzmannSuppression = false;
  double dx = -1;
  double dt = -1;

  std::vector<std::string> args;

  if (argc > 1) {
    // See https://stackoverflow.com/questions/15344714/convert-command-line-argument-to-string
    // +1 for the starting pointer so we skip over the executable name argv[0].
    args.assign(argv, argv + argc);
  }

  std::string outputFolderName = args[1];

  // Check for additional input configuration settings.
  for (int i = 2; i < argc; ++i) {
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

    if (!bUseBoltzmannSuppression && strcmp(argv[i], "-boltz") == 0) {
      bUseBoltzmannSuppression = true;
      continue;
    }

    if (dx < 0 && args[i].compare(0, 3, "dx=") == 0) {
      dx = std::stod(args[i].substr(3, args[i].size() - 3));
      continue;
    }

    if (dt < 0 && args[i].compare(0, 3, "dt=") == 0) {
      dt = std::stod(args[i].substr(3, args[i].size() - 3));
      continue;
    }
  }

  // Set level of screen output
  if (!bDebug && !bTrace) {
    LOGGER(fatal);
  }

  // Construct model
  EffectivePotential::SubcriticalTransitionExample model;

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
  pf.set_t_high(2);
  pf.set_check_vacuum_at_high(false);
  pf.set_check_hessian_singular(false);
  pf.set_x_abs_identical(0.01);
  pf.set_x_abs_jump(0.01);
  pf.set_dt_max_abs(0.01);
  //	pf.set_phase_min_length(0.01);
  pf.set_find_min_trace_abs_step(0.001);
  pf.set_find_min_locate_abs_step(0.001);
  pf.find_phases();

  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(pf);
  tf.set_check_subcritical_transitions(check_subcrit);
  tf.set_assume_only_one_transition(!allow_phase_oscillation);
  tf.find_transitions();
  tf.find_transition_paths(false);

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
