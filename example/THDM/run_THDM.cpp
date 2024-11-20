/**
  2HDM+singlet DM
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>

#include "THDM.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"
#include "thermal_function.hpp"


int main() {

  const bool debug_mode = true;
  
  // Set level of screen  output
  if (debug_mode) {
    LOGGER(debug);
  } else {
    LOGGER(fatal);
  }
  
  // Construct our model
  EffectivePotential::THDM model;
  model.init_params(5., 383.4, 0.258, 0.258, 11.411, -5.576, -5.576);

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
  pf.set_check_vacuum_at_high(false);
  pf.find_phases();
  std::cout << pf;
  
  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(pf);
  tf.find_transitions();
  std::cout << tf;
  
  return 0;

}
