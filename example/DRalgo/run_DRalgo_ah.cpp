/**
 The Abelian Higgs Model in  DRalgo

*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "DRalgo_ah.hpp" 
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"
#include "property_calculator.hpp"


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
  double M, gsq, lam;
  
  if ( argc == 1 ) {
    debug_mode = true;
    M = 100.;
    gsq = 0.42;
    lam =  0.005;
  } else if ( argc >= 4 ) {
    M = atof(argv[1]);
    gsq = atof(argv[2]);
    lam = atof(argv[3]);
  } else {
    std::cout << "Use ./run_DR_ah M gsq lam" << std::endl;
    return 0;
  }
  
  if (debug_mode){
    LOGGER(debug);
    std::cout << "M = " << M << std::endl
              << "gsq = " << gsq << std::endl
              << "lam = " << lam << std::endl;
  } else {
    LOGGER(fatal);
  }
    
  // Construct our model
  EffectivePotential::DR_ah model(M, gsq, lam);
  std::vector<double> in = {M, gsq, lam};

  
  if (debug_mode){
    Eigen::VectorXd phi(1);
    phi[0]=100;
    std::cout << "V(100,300)=" << model.V(phi,500) << std::endl;
    
  }

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
  
  pf.set_seed(0);
  pf.set_check_hessian_singular(false);
  pf.set_t_low(100);
  pf.set_lower_bounds({-1E5});
  pf.set_upper_bounds({1E5});
  
  try {
    pf.find_phases();
  } catch (...) {
    std::cout << "M = " << M << ",\t"
              << "gsq = " << gsq << ",\t"
              << "lam = " << lam << "\t"
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
    std::cout << "M = " << M << ",\t"
              << "gsq = " << gsq << ",\t"
              << "lam = " << lam << "\t"
              << "found 0 transition!" << std::endl;
    std::vector<double> out = {-2, 0, 0, 0, 0, 0};
    output_file << toString(in, out) << std::endl;
    return 0;
  }
  
  // Make PropertyCalculator object and calculate the transition properties
  PhaseTracer::PropertyCalculator pc(tf);
  
  auto profile = pc.findProfile();
  
  std::ofstream file("test_data.txt");
  for (int jj=0; jj< profile.R.size(); jj++){
    std::cout << "R, phi, dphi = " << profile.R[jj] << ", " << profile.Phi(jj) << ", "  << profile.dPhi(jj) << std::endl;
    file << profile.R[jj] << ", " << profile.Phi(jj) << ", "  << profile.dPhi(jj) << std::endl;
  }
  file.close();
  
  std::vector<double> out = {(float)t.size(), t[0].TC, t[0].true_vacuum[0], t[0].false_vacuum[0]};

  output_file << toString(in, out) << std::endl;
  output_file.close();
  return 0;
}
