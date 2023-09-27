/**
 The scale example BubbleProfiler
 ./run_BP_scale 1. 0.1 2.
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "BP_scale.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"


int main(int argc, char* argv[]) {

  double E, alpha, scale;
  
  if ( argc == 4 ) {
    E = atof(argv[1]);
    alpha = atof(argv[2]);
    scale = atof(argv[3]);
  } else {
    std::cout << "Use ./run_BP_scale E alpha scale" << std::endl;
    return 0;
  }
  
  LOGGER(debug);
    
  // Construct our model
  EffectivePotential::BP_scale model(E, alpha, scale);

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
//  pf.set_seed(0);
  pf.set_check_hessian_singular(false);
  pf.set_check_vacuum_at_high(false);
  pf.set_guess_points({Eigen::VectorXd::Zero(1)});
  
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
