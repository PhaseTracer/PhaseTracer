/**
  Z2 real scalar singlet extension of
  the Standard Model 
  
  OS-like
  
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>

#include "models/xSM_OSlike.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"
#include "thermal_function.hpp"

std::string toString(std::vector<double> in, std::vector<double> out, std::vector<double> flags) {
  std::stringstream data_str;
  for (auto i : in    ) data_str << i << "\t";
  for (auto i : out   ) data_str << i << "\t";
  for (auto i : flags ) data_str << i << "\t";
  return data_str.str();;
}

int main(int argc, char* argv[]) {

  std::ofstream output_file;  
  output_file.open("output.txt");

  bool debug_mode = false;
  double ms, lambda_s, lambda_hs;
  double daisy_flag;
  
  if ( argc == 1 ) {
    debug_mode = true;
    // Match choices in 1808.01098
//    lambda_hs = 0.24;
//    ms = 0.5 * SM::mh;
//    double lambda_s_min = 2. / square(SM::mh * SM::v) *
//                          square(square(ms) - 0.5 * lambda_hs * square(SM::v));
//    lambda_s =  lambda_s_min + 0.1;
//    daisy_flag = 1;
    
    ms = 65.;
    lambda_s =  0.1;
    lambda_hs = 0.3;
    daisy_flag = 2;
  } else if ( argc >= 9 ) {
    ms = atof(argv[1]);
    lambda_s = atof(argv[2]);
    lambda_hs = atof(argv[3]);
    daisy_flag = atoi(argv[6]);
  } else {
    std::cout << "Use ./run_xSM_OSlike ms lambda_s lambda_hs" << std::endl;
    return 0;
  }

  if (debug_mode){
    LOGGER(debug);
    std::cout << "ms = " << ms << std::endl
              << "lambda_s = " << lambda_s << std::endl
              << "lambda_hs = " << lambda_hs << std::endl
              << "daisy_term = " << ( daisy_flag == 0  ? "None" : ( daisy_flag == 1 ? "Parwani" : "ArnoldEspinosa")) << std::endl;

  } else {
    LOGGER(fatal);
  }
  
  // Construct our model
  EffectivePotential::xSM_OSlike model(lambda_hs, lambda_s, ms);
  
  std::vector<double> in ={ms, lambda_s, lambda_hs};
  std::vector<double> flags ={0, daisy_flag, 1, model.get_renormalization_scale()};
  
  // Choose Daisy method 
  if (daisy_flag == 0){
    model.set_daisy_method(EffectivePotential::DaisyMethod::None);
  } else if (daisy_flag == 1){
    model.set_daisy_method(EffectivePotential::DaisyMethod::Parwani);
  } else if (daisy_flag == 2){
    model.set_daisy_method(EffectivePotential::DaisyMethod::ArnoldEspinosa);
  } else {
    std::cout << "Wrong daisy flag" << std::endl;
  }
  
  if (debug_mode) {
  
    std::cout << std::setprecision(16);
    std::cout << std::endl;
    std::cout << "@ after applying the 1L EWSB condition" << std::endl;
    std::cout << "muH2       = "<< model.get_muh_sq() << std::endl;
    std::cout << "muS2       = "<< model.get_mus_sq() << std::endl;   
    std::cout << "lambda_h   = "<< model.get_lambda_h() << std::endl;      
    std::cout << "lambda_s   = "<< model.get_lambda_s() << std::endl;  
    std::cout << "lambda_hs  = "<< model.get_lambda_hs() << std::endl; 
  
      Eigen::VectorXd x(2);
      x <<  SM::v, 0;
      std::cout << "V0=" << model.V0(x) << std::endl;
//      std::cout << "V1=" << model.V1(x,0) << std::endl;
      std::cout << "V1T=" << model.V1T(x,0) << std::endl;
      std::cout << "V=" << model.V(x,0) << std::endl;
      std::cout << "V1=" << model.V(x,0) - model.V0(x) << std::endl;
      
      std::cout << "Numerically derivatives of the full potential at EW VEV:" << std::endl;
      auto d2Vdh2 = model.d2V_dx2(x,0);
      std::cout << "Sqrt[d^2V/dh^2] = "<< std::sqrt(abs(d2Vdh2(0,0))) << std::endl;
      std::cout << "Sqrt[d^2V/ds^2] = "<< std::sqrt(abs(d2Vdh2(1,1))) << std::endl;
  }  
  
  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
      
  pf.set_check_vacuum_at_high(false);
  pf.set_seed(0);
//  pf.set_check_hessian_singular(false);
    
  try {
    pf.find_phases();
  } catch (...) {
    std::cout << "ms = " << ms << ",\t"
              << "lambda_s = " << lambda_s << ",\t"
              << "lambda_hs = " << lambda_hs << "\t"
              << "encounters bug!" << std::endl;
    std::vector<double> out = {-1, 0, 0, 0, 0, 0};
    output_file << toString(in, out, flags) << std::endl;
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
    output_file << toString(in, out, flags) << std::endl;
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
    output_file << toString(in, out, flags) << std::endl;
    return 0;
  }
  
  std::vector<double> out = {(float)t.size(), t[jj].TC, t[jj].true_vacuum[0], t[jj].true_vacuum[1], t[jj].false_vacuum[0], t[jj].false_vacuum[1]};
  
  output_file << toString(in, out, flags) << std::endl;
  output_file.close();  
  return 0;
}
