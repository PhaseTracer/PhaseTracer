/**
  Example program for PhaseTracer.

  We scan a parameter in a two field model and show the phases and phase transitions
  for each value of the parameter.
*/

#include <iostream>
#include <string>
#include <vector>

#include "ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT.hpp"
#include "transition_finder.hpp"
#include "h_bar_expansion.hpp"
#include "phase_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"
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
  bool use_Goldstone_resum = false;
  bool use_tree_level_tadpoles = false;
  bool use_tree_level_beta = false;
  const double Qin = 173;
  const double MhInput = 125.25;
  
  if ( argc == 1 ) {
    debug_mode = true;
    // Compare with xSM_MSbar
    ms = 65;
    lambda_s =  0.1;
    lambda_hs = 0.26;
    Q = 173.;
    xi = 1;
    daisy_flag = 0;
    use_1L_EWSB_in_0L_mass = false;  
    use_Goldstone_resum = true;
    use_tree_level_tadpoles = true;
    use_tree_level_beta = true;
  } else if ( argc >= 9 ) {
    ms = atof(argv[1]);
    lambda_s = atof(argv[2]);
    lambda_hs = atof(argv[3]);
    Q = atof(argv[4]);
    xi = atof(argv[5]);
    daisy_flag = atoi(argv[6]);
    use_1L_EWSB_in_0L_mass = atoi(argv[7]);
    use_Goldstone_resum = atoi(argv[8]);
    if ( argc > 9 ){
      if ( atoi(argv[9]) == 1 ) {
        use_tree_level_tadpoles = true;
      }
      if ( atoi(argv[9]) == 3 ) {
        use_tree_level_tadpoles = true;
        use_tree_level_beta = true;
      }
    }
    
  } else {
    std::cout << "Use ./run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT ms lambda_s lambda_hs Q xi daisy_flag use_1L_EWSB_in_0L_mass " << std::endl;
    return 0;
  }

    if (xi != 1) {
      std::cout << "xi is set to 1 in ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT." << std::endl;
    }
    
    
  std::vector<double> in ={ms, lambda_s, lambda_hs};
  std::vector<double> flags ={Q, xi, daisy_flag, (float)use_1L_EWSB_in_0L_mass, (float)use_Goldstone_resum};

  if (debug_mode){
    LOGGER(debug);
    std::cout << "ms = " << ms << std::endl
              << "lambda_s = " << lambda_s << std::endl
              << "lambda_hs = " << lambda_hs << std::endl
              << "Q = " << Q << std::endl
              << "xi = " << xi << std::endl
              << "daisy_term = " << ( daisy_flag == 0  ? "None" : ( daisy_flag == 1 ? "Parwani" : "ArnoldEspinosa")) << std::endl
              << "use 1-level ewsb in tree-level masses = " << use_1L_EWSB_in_0L_mass << std::endl
              << "use Goldstone resum = " << use_Goldstone_resum << std::endl;

  } else {
    LOGGER(fatal);
  }

  // Construct our model
  EffectivePotential::ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT model;
//  if (debug_mode) model.set_debug(true);

  model.set_use_tree_level_tadpoles(use_tree_level_tadpoles);
  model.set_use_tree_level_beta(use_tree_level_beta);
  model.set_use_1L_EWSB_in_0L_mass(use_1L_EWSB_in_0L_mass);
  model.set_use_Goldstone_resum(use_Goldstone_resum);
  model.set_input({Qin, Qin, MhInput, ms, lambda_s, lambda_hs});
  model.Run_pars_to(Q,1E-2);
  // Choose Daisy method
  if (daisy_flag == 0){
    model.set_daisy_method(EffectivePotential::DaisyMethod::None);
  } else {
    std::cout << "No daisy correction for PRM method!" << std::endl;
    return 0;
  }
  
  
  if (debug_mode) {
    std::cout << std::setprecision(16);

    std::cout << "EW VEV = " << model.get_EW_VEV() << std::endl;
    std::cout << "g = " << model.get_g() << std::endl;
    std::cout << "gp = " << model.get_gp() << std::endl;
    std::cout << "yt = " << model.get_yt() << std::endl;
    std::cout << "yb = " << model.get_yb() << std::endl;
    std::cout << "ytau = " << model.get_ytau() << std::endl;

    std::cout << std::endl;
    std::cout << "@ after applying the 1L EWSB condition" << std::endl;
    std::cout << "m_s        = "<< model.get_ms() << std::endl;
    std::cout << "muH2       = "<< model.get_muh_sq() << std::endl;
    std::cout << "muS2       = "<< model.get_mus_sq() << std::endl;   
    std::cout << "lambda_h   = "<< model.get_lambda_h() << std::endl;      
    std::cout << "lambda_s   = "<< model.get_lambda_s() << std::endl;  
    std::cout << "lambda_hs  = "<< model.get_lambda_hs() << std::endl; 
  
    std::cout << std::endl;
    std::cout << "@ EWSB VEV" << std::endl;
    Eigen::VectorXd test(2);
    test <<  model.get_EW_VEV(), 0;
    double Ttest = 100;
    
    auto mh_check =  model.get_scalar_masses_sq(test,1);
    auto mV_check = model.get_vector_masses_sq(test);
    auto mf_check = model.get_fermion_masses_sq(test);
    std::cout << "ms = "<< std::sqrt(mh_check[0]) << std::endl;
    std::cout << "mh = "<< std::sqrt(mh_check[4]) << std::endl;
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
//    return 0;
  }
    
  // Make PhaseFinder object and find the phases
  PhaseTracer::HbarExpansion hb(model);
  hb.set_seed(0);
  Eigen::ArrayXd pseudo(2);
  pseudo << 0., model.get_v_tree_s();
  hb.add_pseudo_phase(pseudo);
  hb.find_phases();
  
  
  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(hb);
  tf.set_TC_tol_rel(1e-12);
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
  
  // Find gamma using high-temperature expansion
  PhaseTracer::HTExpansion ht(model);
  ht.set_seed(0);
  const double TC = tf.get_transitions()[0].TC;
  const auto ht_minima = ht.find_minima_at_t(TC);


  // Use minima with greatest Higgs/scalar
  double vs, vh;
  for (const auto& m : ht_minima) {
    if (debug_mode) {
      std::cout << "HT minimum = " << m.x << std::endl;
    }
    vh = std::max(vh, std::abs(m.x(0)));
    vs = std::max(vs, std::abs(m.x(1)));
  }
  
  if (debug_mode) {
    std::cout << hb;
    std::cout << std::setprecision(15) << tf;
    std::cout << "gamma_HT = " << vh / TC << std::endl;
  }
  
  if (ht_minima.size()==1){
    std::vector<double> out = {0, TC, vh, vs, vh, vs};
    output_file << toString(in, out, flags) << std::endl;
  } else {
    std::vector<double> out = {1, TC, 0.0, vs, vh, 0.0};
    output_file << toString(in, out, flags) << std::endl;
  }
  output_file.close();
  
  return 0;
}
