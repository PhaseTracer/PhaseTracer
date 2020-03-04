/**
  C2HDM in BSMPT
*/

#include <iostream>
#include <string>
#include <vector>

#include "Minimizer.h"
#include "IncludeAllModels.h"
#include "ClassPotentialC2HDM.h"

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
  std::string parameters = "1	3.29771	0.274365	4.71019	-2.23056	-2.43487	0.124948	2706.86	4.64487";
  EffectivePotential::BSMPTPotential<Class_Potential_C2HDM> model(parameters);

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);

  // Model-specific settings for finding the phases
  pf.set_t_high(200.);
  
  // Guesses for location of vacuum
  Eigen::VectorXd guess(4);
  guess << 0., 52, 240, 0.;
  pf.set_guess_points({guess, Eigen::VectorXd::Zero(4)});
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
  PhaseTracer::phase_plotter(tf, "C2HDM");
  PhaseTracer::potential_line_plotter(model, tf.get_transitions(), "C2HDM");
}
