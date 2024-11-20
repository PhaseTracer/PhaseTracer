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
#include "action_calculator.hpp"


int main(int argc, char* argv[]) {

  LOGGER(debug);
    
  // Construct our model
  EffectivePotential::BP_2d model;

  // Make PhaseFinder object and find the vacuums
  PhaseTracer::PhaseFinder pf(model);
  pf.set_guess_points({Eigen::VectorXd::Zero(2)});
  auto vacuums = pf.get_minima_at_t_low();
  
  // Make ActionCalculator object and calcualte the action
  PhaseTracer::ActionCalculator ac(model);
  ac.set_action_calculator(PhaseTracer::ActionMethod::BubbleProfiler);
  
  double action = ac.get_action(vacuums[0].x,Eigen::VectorXd::Zero(2),0); // 0 is the temperature, which is not used in this potential
  
  std::cout << "S = " << action << std::endl;
  return 0;
}
