/**
 The 2d example BubbleProfiler
 ./run_BP_2d
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "BP_2d.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"


int main(int argc, char* argv[]) {

  LOGGER(debug);
    
  // Construct our model
  EffectivePotential::BP_2d model;

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
//  pf.set_seed(0);
  pf.set_check_hessian_singular(false);
  pf.set_check_vacuum_at_high(false);
  pf.set_guess_points({Eigen::VectorXd::Zero(2)});
  
  pf.find_phases();
  std::cout << pf;

  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(pf);
//  tf.find_transitions();
//  std::cout << tf;
  
  const auto phases = pf.get_phases();
  auto phase1 = phases[1];
  auto phase2 = phases[0];
  auto s = tf.get_action(phase1, phase2, 0);
  auto vacua = tf.get_vacua_at_T(phase1, phase2, 0);
  std::cout << vacua[0][0] << " " << vacua[1][0] << " " << s << std::endl;
  
  return 0;
}
