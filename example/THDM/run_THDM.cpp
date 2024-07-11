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

  const bool debug_mode = false;
  
  // Set level of screen  output
  if (debug_mode) {
    LOGGER(debug);
  } else {
    LOGGER(fatal);
  }
  
  // Construct our model
  EffectivePotential::THDM model;
  model.set_daisy_method(EffectivePotential::DaisyMethod::ArnoldEspinosa);
  
  model.init_params(3., 383.4, 0.258, 0.258, 11.411, -5.576, -5.576);

//  Print component of eff potential
Eigen::VectorXd x(2);
x << 200, 100;
//auto SM = model.get_scalar_masses_sq(x, 0);
std::cout << std::setprecision(16);
//std::cout << SM[0] << std::endl;
//std::cout << SM[1] << std::endl;
//std::cout << SM[2] << std::endl;
//std::cout << SM[3] << std::endl;
//std::cout << SM[4] << std::endl;
//std::cout << SM[5] << std::endl;
//std::cout <<"------------"<<std::endl;
//auto VM = model.get_vector_masses_sq(x);
//std::cout << VM[0] << std::endl;
//std::cout << VM[1] << std::endl;
//std::cout << VM[2] << std::endl;
//std::cout << VM[3] << std::endl;
//std::cout << VM[4] << std::endl;
//std::cout << VM[5] << std::endl;
//std::cout <<"------------"<<std::endl;
//auto FM = model.get_fermion_masses_sq(x);
//std::cout << FM[0] << std::endl;
//std::cout << FM[1] << std::endl;

//auto CT_term = model.get_counter_term(246., 1e-1);
//std::cout << CT_term[0]<<std::endl;
//std::cout << CT_term[1]<<std::endl;
//std::cout << CT_term[2]<<std::endl;
//std::cout << CT_term[3]<<std::endl;
//std::cout << CT_term[4]<<std::endl;

std::cout << "V0     = "<< model.V0(x)<< std::endl;
std::cout << "V1    = "<< model.V1(x) << std::endl;
std::cout << "Vct = "<< model.counter_term(x,0) << std::endl;
std::cout << "Vdaisy = "<< model.daisy(x,200) << std::endl;
std::cout << "V1T   = "<< model.V1T(x,200) << std::endl;
std::cout << "V     = "<< model.V(x,200) << std::endl;

//PhaseTracer::PhaseFinder pf(model);
//pf.set
//pf.set_check_vacuum_at_high(false);
//pf.find_phases();
//std::cout << pf;

return 0;

}
