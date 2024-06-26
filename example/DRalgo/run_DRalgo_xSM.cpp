/**
 The xSM in  DRalgo

*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "DRalgo_xSM.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"

std::string toString(std::vector<double> in, std::vector<double> out) {
  std::stringstream data_str;
  for (auto i : in    ) data_str << i << "\t";
  for (auto i : out   ) data_str << i << "\t";
  return data_str.str();;
}

int main(int argc, char* argv[]) {

  std::ofstream output_file;
  output_file.open("output.txt");

  bool debug_mode = false;
  double Ms, lambda_hs, lambda_s;
  
  if ( argc == 1 ) {
    debug_mode = true;
    lambda_hs = 1.3;
    lambda_s = 1.0;
    Ms = 160;
  } else if ( argc >= 4 ) {
    lambda_hs = atof(argv[1]);
    lambda_s = atof(argv[2]);
    Ms = atof(argv[3]);
  } else {
    std::cout << "Use ./run_DR_xSM lambda_hs lambda_s Ms" << std::endl;
    return 0;
  }
  
  if (debug_mode){
    LOGGER(debug);
    std::cout << "Ms = " << Ms << std::endl
              << "lambda_s = " << lambda_s << std::endl
              << "lambda_hs = " << lambda_hs << std::endl;
  } else {
    LOGGER(fatal);
  }
    
  // Construct our model
  EffectivePotential::DR_xSM model(lambda_hs, lambda_s, Ms);
  std::vector<double> in = {Ms, lambda_s, lambda_hs};

  
  if (debug_mode){
    Eigen::VectorXd phi(2);
    phi[0]=100;
    phi[1]=100;
    std::cout << "V(100,100;T=300)=" << model.V(phi,300/2.) << std::endl;
//    return 0;
  }

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
  
  pf.set_seed(0);
  pf.set_check_hessian_singular(false);
  pf.set_check_vacuum_at_high(false);
  pf.set_t_low(40);
  
  try {
    pf.find_phases();
  } catch (...) {
    std::cout << "Ms = " << Ms << ",\t"
              << "lambda_s = " << lambda_s << ",\t"
              << "lambda_hs = " << lambda_hs << "\t"
              << "encounters bug!" << std::endl;
    std::vector<double> out = {-1, 0, 0, 0, 0, 0};
    output_file << toString(in, out) << std::endl;
    return 0;
  }
  std::cout << pf;

  // Make PhaseFinder object
  PhaseTracer::ActionCalculator ac(model);
  ac.set_use_BubbleProfiler(false);
  
  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(pf,ac);
  tf.find_transitions();
  std::cout << tf;

  auto t = tf.get_transitions();
  if (t.size()==0){
    std::cout << "Ms = " << Ms << ",\t"
              << "lambda_s = " << lambda_s << ",\t"
              << "lambda_hs = " << lambda_hs << "\t"
              << "found 0 transition!" << std::endl;
    std::vector<double> out = {-2, 0, 0, 0, 0, 0};
    output_file << toString(in, out) << std::endl;
    return 0;
  }
  
  std::vector<double> out = {(float)t.size(), t[0].TC, t[0].true_vacuum[0], t[0].false_vacuum[0]};

  output_file << toString(in, out) << std::endl;
  output_file.close();
  return 0;
}
