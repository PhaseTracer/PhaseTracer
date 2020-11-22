/**
  Z2 real scalar singlet extension of
  the Standard Model 
  
  MSbar
  
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>

#include "models/xSM_MSbar.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"
#include "thermal_function.hpp"

void help_info(){
  std::cout << "Wrong command! Please run 'run_xSM_MSbar X Q' "<< std::endl;
  std::cout << "  X=1:  Parwani method, lambda_hs = 0.31 "<< std::endl;
  std::cout << "  X=2:  Parwani method, lambda_hs = 0.2 ~ 0.4 "<< std::endl;
  std::cout << "  X=3:  ArnoldEspinosa method, lambda_hs = 0.31 "<< std::endl;
  std::cout << "  X=4:  ArnoldEspinosa method, lambda_hs = 0.2 ~ 0.4 "<< std::endl;
  
  std::cout << "  Q=0:  renormalization scale = 0.5*m_top "<< std::endl;
  std::cout << "  Q=1:  renormalization scale = 1.0*m_top "<< std::endl;
  std::cout << "  Q=2:  renormalization scale = 2.0*m_top "<< std::endl;
  
}

int main(int argc, char* argv[]) {

  if (argc < 3){
    help_info();
    return 0;
  }

  std::ofstream output_file;
  bool debug_mode=false;
  bool Parwani=true;
  std::string out_name = "MSbar_";
  switch (atoi(argv[1])) {
    case 1:
      debug_mode=true;
      break;
    case 2:
      break;
    case 3:
      debug_mode=true;
      Parwani=false;
      break;
    case 4:
      Parwani=false;
      break;
    default:
      help_info();
      return 0;
  }
  
  if (Parwani) {
    out_name += "Parwani";
  } else {
    out_name += "ArnoldEspinosa";
  }
  
  double Q = SM::mtop;
  switch (atoi(argv[2])) {
    case 0:
      Q = 0.5*SM::mtop;
      out_name += "_0.5mt";
      break;
    case 1:
      Q = 1.0*SM::mtop;
      out_name += "_mt";
      break;
    case 2:
      Q = 2.0*SM::mtop;
      out_name += "_2mt";
      break;
    default:
      help_info();
      return 0;
  }
  
  output_file.open(out_name+".txt");
  
  double bins_lambda_hs;
  double lambda_hs;
  if (debug_mode){
    LOGGER(debug);
    bins_lambda_hs = 1;
    lambda_hs = 0.31;
  }else {
    bins_lambda_hs = 50;
    LOGGER(fatal);
  }

  const double xi = 0;
  const bool tree_level_tadpoles = false;
  const bool tree_ewsb = false;
  
  for (double ii = 0; ii < bins_lambda_hs; ii++) {
    if (not debug_mode){
      lambda_hs = 0.2 / bins_lambda_hs * ii+0.2;
    }
    // Construct our model
    auto model = EffectivePotential::xSM_MSbar::from_tadpoles(lambda_hs, Q, xi, tree_level_tadpoles, tree_ewsb);
    
    model.set_daisy_method(EffectivePotential::DaisyMethod::None);
//  if (Parwani){
//    model.set_daisy_method(EffectivePotential::DaisyMethod::Parwani);
//  } else {
//    model.set_daisy_method(EffectivePotential::DaisyMethod::ArnoldEspinosa);
//  }


    if (debug_mode) {
      Eigen::VectorXd x(2);
      x <<  SM::v, 0;
      std::cout << "Numerically derivatives of the full potential at EW VEV:" << std::endl;
      auto d2Vdh2 = model.d2V_dx2(x,0);
      std::cout << std::setprecision(16);
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
    if (debug_mode) PhaseTracer::phase_plotter(tf, out_name);
  }
  output_file.close();
  return 0;
}
