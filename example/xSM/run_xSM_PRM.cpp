/**
  Example of h-bar expansion using xSM model.
*/

#include <iostream>

#include "models/xSM_MSbar.hpp"
#include "models/SM_parameters.hpp"
#include "transition_finder.hpp"
#include "h_bar_expansion.hpp"
#include "phase_plotter.hpp"
#include "potential_plotter.hpp"
#include "potential_line_plotter.hpp"
#include "logger.hpp"

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
  bool use_covariant_gauge = false;
  
  if ( argc == 1 ) {
    debug_mode = true;
    // Compare with run_ScalarSingletZ2DMMhInput_withSingletVEVinPT
    ms = 65.;
    lambda_s =  0.1;
    lambda_hs = 0.3;
    Q = 173.;
    xi = 1.;
    daisy_flag = 0;
    use_1L_EWSB_in_0L_mass = false;
    use_Goldstone_resum = true;
    tree_level_tadpoles = true;
    
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
      // default
      // tree_level_tadpoles = false
      // use_covariant_gauge = true
      if ( atoi(argv[9]) == 1 ){
        tree_level_tadpoles = true;
      } else if ( atoi(argv[9]) == 2 ){
        use_covariant_gauge = true;
      } else if ( atoi(argv[9]) == 3 ){
        use_covariant_gauge = true;
        tree_level_tadpoles = true;
      } 
    }
  } else {
    std::cout << "Use ./run_xSM_MSbar ms lambda_s lambda_hs Q xi daisy_flag use_1L_EWSB_in_0L_mass use_Goldstone_resum tree_level_tadpoles" << std::endl;
    return 0;
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
              << "tree-level tadpoles = " << tree_level_tadpoles << std::endl
              << "use 1-level ewsb in tree-level masses = " << use_1L_EWSB_in_0L_mass << std::endl
              << "use Goldstone resum = " << use_Goldstone_resum << std::endl;

  } else {
    LOGGER(fatal);
  }
    
  // Construct our model
  auto model = EffectivePotential::xSM_MSbar::from_tadpoles(lambda_hs, lambda_s, ms, Q, xi, use_covariant_gauge, use_1L_EWSB_in_0L_mass, use_Goldstone_resum, tree_level_tadpoles);

  // Choose Daisy method
  if (daisy_flag == 0){
    model.set_daisy_method(EffectivePotential::DaisyMethod::None);
  } else {
    std::cout << "No daisy correction for PRM method!" << std::endl;
    return 0;
  }
  
  if (debug_mode) {
      std::cout << "muh_sq = " << model.get_muh_sq() << std::endl
                << "mus_sq = " << model.get_mus_sq() << std::endl
                << "lambda_h = " << model.get_lambda_h() << std::endl;
  }
  
  // Make PhaseFinder object and find the phases
  PhaseTracer::HbarExpansion hb(model);
  hb.set_seed(0);
  Eigen::ArrayXd pseudo(2);
  pseudo << 0., model.get_v_tree_s();
  hb.add_pseudo_phase(pseudo);
  hb.find_phases();
  if (debug_mode) std::cout << hb;
  
  
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
  double vs=0, vh=0;
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
