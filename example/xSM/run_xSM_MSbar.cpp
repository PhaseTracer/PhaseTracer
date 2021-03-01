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
#include <random>

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
  std::cout << "  x1=2: lambda_hs = 0.1 ~ 0.4, lambda_s = 0.01~0.2, m_s=m_h/2"<< std::endl;
  std::cout << "  x1=3: lambda_hs = 0.1 ~ 0.4, lambda_s = 0.1, m_s= 10~100 GeV"<< std::endl;
  std::cout << "  x1=4: lambda_hs = 0.3, lambda_s = 0.01~0.2, m_s= 10~100 GeV"<< std::endl;
  std::cout << "  x1=5: lambda_hs = 0.1 ~ 0.4, lambda_s = 0.01~0.2, m_s= 10~100 GeV"<< std::endl;
      
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
  bool scan_lhs = atoi(argv[1]) == 1;
  bool scan_lhs_ls = atoi(argv[1]) == 2;
  bool scan_lhs_ms = atoi(argv[1]) == 3;
  bool scan_ls_ms = atoi(argv[1]) == 4;
  bool scan_lhs_ls_ms = atoi(argv[1]) == 5;
  
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
  
  if (scan_lhs_ls) {
    out_name += "_lhs_ls";
  } else if (scan_lhs_ms) {
    out_name += "_lhs_ms";  
  } else if (scan_ls_ms) {
    out_name += "_ls_ms";  
  } else if (scan_lhs_ls_ms) {
    out_name += "_lhs_ls_ms";  
  }
  
  if (xi > 0){
    out_name += "_xi";
  }
  
  output_file.open(out_name+".txt");
  
  double n_bin_x;
  double n_bin_y;
  if (debug_mode){
    LOGGER(debug);
    n_bin_x = 0;
    n_bin_y = 0;
  }else {
    n_bin_x = 50;
    if (scan_lhs){
      n_bin_y = 0;
    } else if (scan_lhs_ms or scan_lhs_ls or scan_ls_ms) {
      n_bin_y = 50;
    } else if (scan_lhs_ls_ms){
      n_bin_x = 5000;
      n_bin_y = 0;
    }
    LOGGER(fatal);
  }

  std::cout << "renormal Q = " << Q/SM::mtop << "*m_top" << std::endl;
  std::cout << "daisy_term = " << (Parwani  ? "Parwani" : "ArnoldEspinosa") << std::endl;
  std::cout << "tree_ewsb  = " << (tree_ewsb  ? "true" : "false") << std::endl;
  std::cout << "xi  = " << xi << std::endl;
  
  const bool tree_level_tadpoles = false;
  double lambda_hs;
  double lambda_s;
  double ms;
  double lhs_min = 0.1;
  double lhs_del = 0.4;
  double ls_min = 0.01;
  double ls_del = 0.2;
  double ms_min = 10;
  double ms_del = 110;
  
  std::default_random_engine random(1);
  std::uniform_real_distribution<double> rand_lhs(lhs_min, lhs_min+lhs_del);
  std::uniform_real_distribution<double> rand_ls(ls_min, ls_min+ls_del);
  std::uniform_real_distribution<double> rand_ms(ms_min, ms_min+ms_del);
  
  for (double ii = 0; ii <= n_bin_x; ii++) {
    for (double jj = 0; jj <= n_bin_y; jj++) {
      if (debug_mode){
//        lambda_hs = 0.24;
//        // Match choices in 1808.01098
//        ms = 0.8 * SM::mh;
//        double lambda_s_min = 2. / square(SM::mh * SM::v) *
//        square(square(ms) - 0.5 * lambda_hs * square(SM::v));
//        lambda_s =  lambda_s_min + 0.1;
        // BK point
//        lambda_hs = 0.4;
//        ms = 60;
//        lambda_s =  0.15;
        // test point for FS
        lambda_hs = 0.25;
        ms = 67.3656;
        lambda_s =  0.1;
        
      } else {
        if (scan_lhs){
          lambda_hs = 0.4 / n_bin_x * ii+0.2;
          // Match choices in 1808.01098
          ms = 0.5 * SM::mh;
          double lambda_s_min = 2. / square(SM::mh * SM::v) *
          square(square(ms) - 0.5 * lambda_hs * square(SM::v));
          lambda_s =  lambda_s_min + 0.1;
        } else if (scan_lhs_ls) {
          lambda_hs = lhs_del / n_bin_x * ii + lhs_min;
          lambda_s = ls_del / n_bin_y * jj + ls_min;
          ms = 0.5 * SM::mh;
        } else if (scan_lhs_ms) {
          lambda_hs = lhs_del / n_bin_x * ii + lhs_min;
          lambda_s = 0.1;
          ms = ms_del / n_bin_x * jj + ms_min;
        } else if (scan_ls_ms) {
          lambda_hs = 0.3;
          lambda_s = ls_del / n_bin_y * ii + ls_min;
          ms = ms_del / n_bin_x * jj + ms_min;
        } else if (scan_lhs_ls_ms){
          lambda_hs = rand_lhs(random);
          lambda_s = rand_ls(random);
          ms = rand_ms(random);
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
