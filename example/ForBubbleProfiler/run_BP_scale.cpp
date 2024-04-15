/**
 The scale example BubbleProfiler
 ./run_BP_scale 1. 0.6 200.
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "BP_scale.hpp"
#include "action_calculator.hpp"


int main(int argc, char* argv[]) {

  double E, alpha, scale;
  
  if ( argc == 4 ) {
    E = atof(argv[1]);
    alpha = atof(argv[2]);
    scale = atof(argv[3]);
  } else {
    std::cout << "Use ./run_BP_scale E alpha scale" << std::endl;
    std::cout << "E.g. ./run_BP_scale 1 0.6 200" << std::endl;
    return 0;
  }
  
  LOGGER(debug);
    
  // Construct our model
  EffectivePotential::BP_scale model(E, alpha, scale);

  // Make ActionCalculator object and calcualte the action
  PhaseTracer::ActionCalculator ac(model);
  
  Eigen::VectorXd true_vacuum(1);
  true_vacuum << 0.;
  Eigen::VectorXd false_vacuum(1);
  false_vacuum << scale;
  
  double action = ac.get_action(true_vacuum,false_vacuum,0); // 0 is the temperature, which is not used in this potential
  
  std::cout << "S = " << action << std::endl;
  
  return 0;
}
