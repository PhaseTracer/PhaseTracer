/**
  R2HDM in BSMPT
*/

#include <iostream>
#include <string>
#include <vector>

#include "Minimizer.h"
#include "IncludeAllModels.h"
#include "ClassPotentialR2HDM.h"

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
  std::string parameters = "1	2.740594787	0.2423556498	5.534491052	-2.585467181	-2.225991025	7738.56	4.63286";
  EffectivePotential::BSMPTPotential<BSMPT::Models::Class_Potential_R2HDM> model(parameters);

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
  pf.set_t_high(200);
  
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
  PhaseTracer::phase_plotter(tf, "R2HDM");
  PhaseTracer::potential_line_plotter(model, tf.get_transitions(), "R2HDM");
}
