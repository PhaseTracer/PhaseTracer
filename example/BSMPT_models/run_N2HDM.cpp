/**
  RN2HDM in BSMPT
*/

#include <iostream>
#include <string>
#include <vector>

#include "Minimizer.h"
#include "IncludeAllModels.h"
#include "ClassPotentialRN2HDM.h"

#include "models/BSMPT.hpp"
#include "transition_finder.hpp"
#include "phase_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"
#include "potential_line_plotter.hpp"

int main() {

  LOGGER(debug);

  // Construct our model - for BSMPT models, we provide a string containing
  // parameter values
  std::string parameters = "1	0.300812	0.321809	-0.133425	4.11105	-3.84178	9.46329	-0.750455	0.743982	293.035	5.91129	4842.28";
  EffectivePotential::BSMPTPotential<Class_Potential_RN2HDM> model(parameters);

  // Make PhaseFinder object
  PhaseTracer::PhaseFinder pf(model);
  pf.set_t_high(200);
  // At T=200 GeV, singlet VEV is not at origin
  pf.set_check_vacuum_at_high(false);
  
  // Guesses for location of vacuum
  Eigen::VectorXd guess(5);
  guess << 0., 0., 41., 242., 293.;
  pf.set_guess_points({guess, Eigen::VectorXd::Zero(5)});
  pf.set_n_test_points(0);
  
  // Find the phases
  pf.find_phases();
  std::cout << pf;

  // Make TransitionFinder object and find the transitions
  // This takes the phase finder object as an argument. We must have already
  // populated the phases by e.g. the find_phases method, as above
  PhaseTracer::TransitionFinder tf(pf);
  tf.find_transitions();
  std::cout << tf;

  // Print the data in a particular format for plotting
  PhaseTracer::phase_plotter(tf, "N2HDM");
  PhaseTracer::potential_line_plotter(model, tf.get_transitions(), "N2HDM");
}
