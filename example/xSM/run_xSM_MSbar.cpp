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
#include "potential_plotter.hpp"

void help_info(){
  std::cout << "Wrong command! Please run 'run_xSM_MSbar x1 x2 x3 x4 x5' "<< std::endl;
  
  
  std::cout << "  x1=0: lambda_hs = 0.31, lambda_s = lambda_s^{min}+0.1, m_s=m_h/2"<< std::endl;
  std::cout << "  x1=1: lambda_hs = 0.2 ~ 0.4, lambda_s = lambda_s^{min}+0.1, m_s=m_h/2"<< std::endl;
  std::cout << "  x1=2: lambda_hs = 0.2 ~ 0.4, lambda_s = 0.01~0.2, m_s=m_h/2"<< std::endl;
  std::cout << "  x1=3: lambda_hs = 0.2 ~ 0.4, lambda_s = 0.1, m_s= 10~100 GeV"<< std::endl;
  
  std::cout << "  x2=0: Parwani method"<< std::endl;
  std::cout << "  x2=1: ArnoldEspinosa method"<< std::endl;
  
  std::cout << "  x3=0:  renormalization scale = 0.5*m_top "<< std::endl;
  std::cout << "  x3=1:  renormalization scale = 1.0*m_top "<< std::endl;
  std::cout << "  x3=2:  renormalization scale = 2.0*m_top "<< std::endl;
  
  std::cout << "  x4=0:  tree_ewsb = false "<< std::endl;
  std::cout << "  x4=1:  tree_ewsb = true "<< std::endl;
  
  std::cout << "  x5:  xi "<< std::endl;
}

int main(int argc, char* argv[]) {

  double xi = 0;
  if (argc < 5){
    help_info();
    return 0;
  } else if (argc > 5 ){
    xi = atof(argv[5]);
  }
  
  std::ofstream output_file;
  std::string out_name = "MSbar_";
    
  bool debug_mode = atoi(argv[1]) == 0;
  bool scan_lambda_hs = atoi(argv[1]) == 1;
  bool scan_lambda_hs_s = atoi(argv[1]) == 2;
  bool scan_lambda_hs_ms = atoi(argv[1]) == 3;
  
  bool Parwani = atoi(argv[2]) == 0;
  if (Parwani) {
    out_name += "Parwani";
  } else {
    out_name += "ArnoldEspinosa";
  }
  double Q = SM::mtop;
  switch (atoi(argv[3])) {
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

  const bool tree_ewsb = atoi(argv[4]) == 1;
  if (tree_ewsb) {
    out_name += "_TreeEWSB";
  } else {
    out_name += "_NoTreeEWSB";
  }
  
  if (scan_lambda_hs_s) {
    out_name += "_lhs_ls";
  }
  if (scan_lambda_hs_ms) {
    out_name += "_lhs_ms";  
  }
  
  if (xi > 0){
    out_name += "_xi";
  }
  
  output_file.open(out_name+".txt");
  
  double n_bin_lambda_hs;
  double n_bin_y=0;
  double lambda_hs;
  double lambda_s;
  double ms;
  if (debug_mode){
    LOGGER(debug);
    n_bin_lambda_hs = 0;
    lambda_hs = 0.24;
    // Match choices in 1808.01098
    ms = 0.5 * SM::mh;
    double lambda_s_min = 2. / square(SM::mh * SM::v) *
    square(square(ms) - 0.5 * lambda_hs * square(SM::v));
    lambda_s =  lambda_s_min + 0.1;
    
  }else {
    n_bin_lambda_hs = 100;
    if (not scan_lambda_hs){
      n_bin_y = 50;
    }
    LOGGER(fatal);
  }

  std::cout << "lambda_hs  = " << (debug_mode ? std::to_string(lambda_hs) : "0.2 ~ 0.4") << std::endl;
  std::cout << "renormal Q = " << Q/SM::mtop << "*m_top" << std::endl;
  std::cout << "daisy_term = " << (Parwani  ? "Parwani" : "ArnoldEspinosa") << std::endl;
  std::cout << "tree_ewsb  = " << (tree_ewsb  ? "true" : "false") << std::endl;
  std::cout << "xi  = " << xi << std::endl;
  
  const bool tree_level_tadpoles = false;
  
  for (double ii = 0; ii <= n_bin_lambda_hs; ii++) {
    for (double jj = 0; jj <= n_bin_y; jj++) {
      if (not debug_mode){
        lambda_hs = 0.4 / n_bin_lambda_hs * ii+0.2;
        if (scan_lambda_hs_s) {
          lambda_s = 0.2 / n_bin_y * jj +0.01;
          ms = 0.5 * SM::mh;
        }
        if (scan_lambda_hs_ms) {
          lambda_s = 0.1;
          ms = 90 / n_bin_lambda_hs * jj + 10;
        }
        if (scan_lambda_hs){
          double lambda_s_min = 2. / square(SM::mh * SM::v) *
          square(square(ms) - 0.5 * lambda_hs * square(SM::v));
          lambda_s =  lambda_s_min + 0.1;
        }
      }
    
      std::cout << "Runing lambda_hs  = " << lambda_hs << ", lambda_s  = " << lambda_s << ", m_s  = " << ms << std::endl;
    
      // Construct our model
    
      auto model = EffectivePotential::xSM_MSbar::from_tadpoles(lambda_hs, lambda_s, ms, Q, xi, tree_level_tadpoles, tree_ewsb);
      std::cout << "iteration converged = " << model.iteration_converged << std::endl;

      
//    model.set_daisy_method(EffectivePotential::DaisyMethod::None);
      if (Parwani){
        model.set_daisy_method(EffectivePotential::DaisyMethod::Parwani);
      } else {
        model.set_daisy_method(EffectivePotential::DaisyMethod::ArnoldEspinosa);
      }


      if (debug_mode) {
        Eigen::VectorXd x(2);
        x <<  SM::v, 0;
        std::cout << "Numerically derivatives of the full potential at EW VEV:" << std::endl;
        auto d2Vdh2 = model.d2V_dx2(x,0);
        std::cout << std::setprecision(16);
        std::cout << "Sqrt[d^2V/dh^2] = "<< std::sqrt(abs(d2Vdh2(0,0))) << std::endl;
        std::cout << "Sqrt[d^2V/ds^2] = "<< std::sqrt(abs(d2Vdh2(1,1))) << std::endl;
      
//      PhaseTracer::potential_plotter(model, 254, "potential", -5., 5, 0.01, -5., 40., 0.1);
//      PhaseTracer::potential_plotter(model, 142.35, "potential", 0., 160, 0.2, -2., 160., 0.2);
//      return 0;
      }
      
      // Make PhaseFinder object and find the phases
      PhaseTracer::PhaseFinder pf(model);
      
      pf.set_check_vacuum_at_high(false);
      pf.set_seed(1);
//      pf.set_check_hessian_singular(false);
//      pf.set_hessian_singular_rel_tol(1.e-6);
    
      try {
        pf.find_phases();
      } catch (...) {
        std::cout << "lambda_hs = " << lambda_hs << " encounters bug!" << std::endl;
        output_file << lambda_hs << "\t" << 0 << "\t";
        output_file << -1 << "\t" << 0 << "\t"
                    << 0 << "\t" << 0 << "\t"
                    << 0 << "\t";
        output_file << lambda_s << "\t" << ms;
        output_file << std::endl;
        continue;
      }
    
      if (debug_mode) std::cout << pf;

      // Make TransitionFinder object and find the transitions
      PhaseTracer::TransitionFinder tf(pf);
      tf.find_transitions();
      if (debug_mode) std::cout << tf;
    
      if (not debug_mode) {
        auto t = tf.get_transitions();
        if (t.size()==0){
          std::cout << "lambda_hs = " << lambda_hs << " found 0 transition!" << std::endl;
          output_file << lambda_hs << "\t" << 0 << "\t";
          output_file << -1 << "\t" << 0 << "\t"
                      << 0 << "\t" << 0 << "\t"
                      << 0 << "\t";
          output_file << lambda_s << "\t" << ms;
          output_file << std::endl;
          continue;
        }
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
        output_file << lambda_s << "\t" << ms;
        output_file << std::endl;
      }
      // Print the data in a particular format for plotting
      if (debug_mode) PhaseTracer::phase_plotter(tf, out_name);
    }
  }
  output_file.close();
  return 0;
}
