/**
  Z2 real scalar singlet extension of
  the Standard Model 
  
  MSbar
  
*/

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
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
  double Q, xi, daisy_flag;
  bool use_1L_EWSB_in_0L_mass;
  bool use_Goldstone_resum = true;
  bool tree_level_tadpoles = false;
  
  if ( argc == 1 ) {
    debug_mode = true;
    // Compare with run_ScalarSingletZ2DMMhInput_withSingletVEVinPT
    ms = 66.50567996545465;
    lambda_s =  0.1;
    lambda_hs = 0.25;
    Q = 173.;
    xi = 1;
    daisy_flag = 1;
    use_1L_EWSB_in_0L_mass = false;  
    use_Goldstone_resum = false;
    
//    // Match choices in 1808.01098
//    lambda_hs = 0.24;
//    ms = 0.5 * SM::mh;
//    double lambda_s_min = 2. / square(SM::mh * SM::v) *
//                          square(square(ms) - 0.5 * lambda_hs * square(SM::v));
//    lambda_s =  lambda_s_min + 0.1;
//    
//    // BK point
//    lambda_hs = 0.4;
//    ms = 60;
//    lambda_s =  0.15;
    
  } else if ( argc == 8 ) {
    ms = atof(argv[1]);
    lambda_s = atof(argv[2]);
    lambda_hs = atof(argv[3]);
    Q = atof(argv[4]);
    xi = atof(argv[5]);

    daisy_flag = atoi(argv[6]);
    use_1L_EWSB_in_0L_mass = atoi(argv[7]);
  } else {
    std::cout << "Use ./run_xSM_MSbar ms lambda_s lambda_hs Q xi daisy_flag use_1L_EWSB_in_0L_mass " << std::endl;
    return 0;
  }

  std::vector<double> in ={ms, lambda_s, lambda_hs};
  std::vector<double> flags ={Q, xi, daisy_flag, (float)use_1L_EWSB_in_0L_mass};

  if (debug_mode){
    LOGGER(debug);
    std::cout << "ms = " << ms << std::endl
              << "lambda_s = " << lambda_s << std::endl
              << "lambda_hs = " << lambda_hs << std::endl
              << "Q = " << Q << std::endl
              << "xi = " << xi << std::endl
              << "daisy_term = " << ( daisy_flag == 0  ? "None" : ( daisy_flag == 1 ? "Parwani" : "ArnoldEspinosa")) << std::endl
              << "tree-level tadpoles = " << tree_level_tadpoles << std::endl
              << "use 1-level ewsb in tree-level masses = " << use_1L_EWSB_in_0L_mass << std::endl;

  } else {
    LOGGER(fatal);
  }
  
  // Construct our model
  auto model = EffectivePotential::xSM_MSbar::from_tadpoles(lambda_hs, lambda_s, ms, Q, xi, tree_level_tadpoles, use_1L_EWSB_in_0L_mass, use_Goldstone_resum);
  if (debug_mode) std::cout << "1-L EWSB iteration converged = " << model.iteration_converged << std::endl;

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
    std::cout << "m_s        = "<< ms << std::endl;
    std::cout << "muH2       = "<< model.get_muh_sq() << std::endl;
    std::cout << "muS2       = "<< model.get_mus_sq() << std::endl;   
    std::cout << "lambda_h   = "<< model.get_lambda_h() << std::endl;      
    std::cout << "lambda_s   = "<< model.get_lambda_s() << std::endl;  
    std::cout << "lambda_hs  = "<< model.get_lambda_hs() << std::endl; 
    
    std::cout << std::endl;
    std::cout << "@ EWSB VEV" << std::endl;
    Eigen::VectorXd test(2);
    test <<  SM::v, 0;
    double Ttest = 100;
    
    auto mh_check =  model.get_scalar_masses_sq(test,1);
    auto mV_check = model.get_vector_masses_sq(test);
    auto mf_check = model.get_fermion_masses_sq(test);
    std::cout << "mh1 = "<< std::sqrt(mh_check[0]) << std::endl;
    std::cout << "mh2 = "<< std::sqrt(mh_check[1]) << std::endl;
    std::cout << "mg0 = "<< std::sqrt(mh_check[2]) << std::endl;
    std::cout << "mg+ = "<< std::sqrt(mh_check[3]) << std::endl;
    std::cout << "MW = "<< std::sqrt(mV_check[0]) << std::endl;
    std::cout << "MZ = "<< std::sqrt(mV_check[1]) << std::endl;
    std::cout << "Mphoton = "<< std::sqrt(mV_check[2]) << std::endl;
    std::cout << "mt = " << std::sqrt(mf_check[0]) << std::endl;
    std::cout << "mb = " << std::sqrt(mf_check[1]) << std::endl;
    std::cout << "mtau = " << std::sqrt(mf_check[2]) << std::endl;
    
    double Vtree = model.V0(test);
    double VCW = model.V1(test);
    double V1T = model.V1T(test, Ttest);
    double Vtot = model.V(test, Ttest);
    std::cout << "Vtree      = "<< Vtree << std::endl;
    std::cout << "VCW        = "<< VCW << std::endl;   
    std::cout << "V1T(T=100) = "<< V1T << std::endl;      
    std::cout << "V(T=100)   = "<< Vtot << std::endl;  
    std::cout << std::endl;
    
    std::cout << "Numerically derivatives of the full potential at EW VEV:" << std::endl;
    auto d2Vdh2 = model.d2V_dx2(test,0);
    std::cout << std::setprecision(16);
    std::cout << "Sqrt[d^2V/dh^2] = "<< std::sqrt(abs(d2Vdh2(0,0))) << std::endl;
    std::cout << "Sqrt[d^2V/ds^2] = "<< std::sqrt(abs(d2Vdh2(1,1))) << std::endl;
    
//    PhaseTracer::potential_plotter(model, 254, "potential", -5., 5, 0.01, -5., 40., 0.1);
//    PhaseTracer::potential_plotter(model, 142.35, "potential", 0., 160, 0.2, -2., 160., 0.2);
//    return 0;
  }
      
  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);    
  pf.set_check_vacuum_at_high(false);
  pf.set_seed(1);
//  pf.set_check_hessian_singular(false);
//  pf.set_hessian_singular_rel_tol(1.e-6);
    
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
  
  int jj=0;
  double gamme_max=0;
  for (int i=0; i<t.size(); i++) {
    double gamma = t[i].gamma;
    if (gamme_max < gamma){
      jj = i;
      gamme_max = gamma;
    }
  }
  std::vector<double> out = {(float)t.size(), t[jj].TC, t[jj].true_vacuum[0], t[jj].true_vacuum[1], t[jj].false_vacuum[0], t[jj].false_vacuum[1]};
  
  output_file << toString(in, out, flags) << std::endl;
  output_file.close();  
  // Print the data in a particular format for plotting
  if (debug_mode) PhaseTracer::phase_plotter(tf, "xSM_MSbar");
  return 0;
}
