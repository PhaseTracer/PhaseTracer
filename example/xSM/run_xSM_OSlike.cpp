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

void help_info(){
  std::cout << "Wrong command! Please run 'run_xSM_OSlike X' "<< std::endl;
  std::cout << "  X=1:  Parwani method, lambda_hs = 0.31 "<< std::endl;
  std::cout << "  X=2:  Parwani method, lambda_hs = 0.2 ~ 0.4 "<< std::endl;
  std::cout << "  X=3:  ArnoldEspinosa method, lambda_hs = 0.31 "<< std::endl;
  std::cout << "  X=4:  ArnoldEspinosa method, lambda_hs = 0.2 ~ 0.4 "<< std::endl;
}

int main(int argc, char* argv[]) {

  if (argc < 2){
    help_info();
    return 0;
  }

  std::ofstream output_file;
  bool debug_mode=false;
  bool Parwani=true;
  std::string plotname = "OSlike_Parwani";
  switch (atoi(argv[1])) {
    case 1:
      debug_mode=true;
      break;
    case 2:
      output_file.open("OSlike_Parwani.txt");
      break;
    case 3:
      debug_mode=true;
      Parwani=false;
      plotname="OSlike_ArnoldEspinosa";
      break;
    case 4:
      Parwani=false;
      output_file.open("OSlike_ArnoldEspinosa.txt");
      break;
    default:
      help_info();
      return 0;
  }
  
  double bins_lambda_hs;
  double lambda_hs;
  if (debug_mode){
    LOGGER(debug);
    bins_lambda_hs = 1;
    lambda_hs = 0.28;
  }else {
    bins_lambda_hs = 100;
    LOGGER(fatal);
  }
  

  // Construct our model
  EffectivePotential::xSM_OSlike model;
  if (Parwani){
    model.set_daisy_method(EffectivePotential::DaisyMethod::Parwani);
  } else {
    model.set_daisy_method(EffectivePotential::DaisyMethod::ArnoldEspinosa);
  }
  model.set_m_s(SM::mh/2.);
  
  for (double ii = 0; ii < bins_lambda_hs; ii++) {
    if (not debug_mode){
      lambda_hs = 0.4 / bins_lambda_hs * ii+0.2;
    }
    model.set_lambda_hs(lambda_hs);
    model.solve_Q();
  
    if (debug_mode) {
      Eigen::VectorXd x(2);
      x <<  SM::v, 0;
      std::cout << "Numerically derivatives of the full potential at EW VEV:" << std::endl;
      auto d2Vdh2 = model.d2V_dx2(x,0);
      std::cout << "Sqrt[d^2V/dh^2] = "<< std::sqrt(abs(d2Vdh2(0,0))) << std::endl;
      std::cout << "Sqrt[d^2V/ds^2] = "<< std::sqrt(abs(d2Vdh2(1,1))) << std::endl;
    }
      
    // Make PhaseFinder object and find the phases
    PhaseTracer::PhaseFinder pf(model);
      
    pf.set_check_vacuum_at_high(false);
    pf.set_seed(0);
    pf.set_check_hessian_singular(false);
    
    try {
      pf.find_phases();
    } catch (...) {
      std::cout << "lambda_hs = " << lambda_hs << " encounters bug!" << std::endl;
      continue;
    }
    
    if (debug_mode) std::cout << pf;

    // Make TransitionFinder object and find the transitions
    PhaseTracer::TransitionFinder tf(pf);
    tf.find_transitions();
    if (debug_mode) std::cout << tf;

    if (not debug_mode) {
      auto t = tf.get_transitions();
      int jj=0;
      double gamme_max=0;
      for (int i=0; i<t.size(); i++) {
        double gamma = t[i].gamma;
        if (gamme_max < gamma){
            jj = i;
            gamme_max = gamma;
        }
      }
      output_file << lambda_hs << "\t" << t.size() << "\t";
      output_file << t[jj].TC << "\t" << t[jj].true_vacuum[0] << "\t"
                << t[jj].true_vacuum[1] << "\t" << t[jj].false_vacuum[0] << "\t"
                << t[jj].false_vacuum[1] << "\t";
      output_file << std::endl;
    }
    // Print the data in a particular format for plotting
    if (debug_mode) PhaseTracer::phase_plotter(tf, plotname);
  }
  output_file.close();
  return 0;
}
