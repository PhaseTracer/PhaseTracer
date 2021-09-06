/**
 Z2 real scalar singlet extension of
 the Standard Model
 
  High-temperature approximation

*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "models/xSM_HT.hpp" 
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
  double ms, lambda_s, lambda_hs;
  
  if ( argc == 1 ) {
    debug_mode = true;
    lambda_hs = 0.3;
    ms = 65;
    lambda_s =  0.1;
  } else if ( argc >= 4 ) {
    ms = atof(argv[1]);
    lambda_s = atof(argv[2]);
    lambda_hs = atof(argv[3]);
  } else {
    std::cout << "Use ./run_xSM_HT ms lambda_s lambda_hs" << std::endl;
    return 0;
  }
  
  if (debug_mode){
    LOGGER(debug);
    std::cout << "ms = " << ms << std::endl
              << "lambda_s = " << lambda_s << std::endl
              << "lambda_hs = " << lambda_hs << std::endl;
  } else {
    LOGGER(fatal);
  }
    
  // Construct our model
  EffectivePotential::xSM_HT model(lambda_hs, lambda_s, ms);
  std::vector<double> in = {ms, lambda_s, lambda_hs};

  
  if (debug_mode){
    std::cout << "Check = " << model.check() << std::endl;
  }
    
  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
      
  pf.set_seed(0);
  pf.set_check_hessian_singular(false);
  
  try {
    pf.find_phases();
  } catch (...) {
    std::cout << "ms = " << ms << ",\t"
              << "lambda_s = " << lambda_s << ",\t"
              << "lambda_hs = " << lambda_hs << "\t"
              << "encounters bug!" << std::endl;
    std::vector<double> out = {-1, 0, 0, 0, 0, 0};
    output_file << toString(in, out) << std::endl;
    return 0;
  }
  if (debug_mode) std::cout << pf;

  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(pf);
  tf.find_transitions();
  if (debug_mode) std::cout << tf;
      
  auto t = tf.get_transitions();
  if (t.size()==0){
    std::cout << "ms = " << ms << ",\t"
              << "lambda_s = " << lambda_s << ",\t"
              << "lambda_hs = " << lambda_hs << "\t"
              << "found 0 transition!" << std::endl;
    std::vector<double> out = {-2, 0, 0, 0, 0, 0};
    output_file << toString(in, out) << std::endl;
    return 0;
  }
  
  // Find the transition with largest gamma from (0,vs) -> (vh,0)
  int jj = -1;
  double gamme_max = 0.;
  for (int i=0; i<t.size(); i++) {
    double gamma = t[i].gamma;
    if (gamme_max < gamma and abs(t[i].false_vacuum[0])<1. and abs(t[i].true_vacuum[1])<1.){
      jj = i;
      gamme_max = gamma;
    }
  }
  
  if (jj<0) {
    std::vector<double> out = {-3, 0, 0, 0, 0, 0};
    output_file << toString(in, out) << std::endl;
    return 0;
  }
  
  std::vector<double> out = {(float)t.size(), t[jj].TC, t[jj].true_vacuum[0], t[jj].true_vacuum[1], t[jj].false_vacuum[0], t[jj].false_vacuum[1]};
  
  output_file << toString(in, out) << std::endl;
  output_file.close();
  return 0;
}
