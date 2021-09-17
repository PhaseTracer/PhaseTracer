/**
  Example of covariant gauge using xSM model.
*/

#include <fstream>
#include <iostream>
#include <iomanip>

#include "models/xSM_MSbar_covariant.hpp"
#include "models/SM_parameters.hpp"
#include "transition_finder.hpp"
#include "phase_finder.hpp"
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
  if ( argc == 1 ) {
    debug_mode = true;
    ms = 65.;
    lambda_s =  0.1;
    lambda_hs = 0.3;
    Q = 173.;
    xi = 1;
    daisy_flag = 1;
    use_1L_EWSB_in_0L_mass = false;
    if ( xi==0 and not use_1L_EWSB_in_0L_mass )
      use_Goldstone_resum = true;
    else 
      use_Goldstone_resum = false;
    
  } else if ( argc >= 9 ) {
    ms = atof(argv[1]);
    lambda_s = atof(argv[2]);
    lambda_hs = atof(argv[3]);
    Q = atof(argv[4]);
    xi = atof(argv[5]);

    daisy_flag = atoi(argv[6]);
    use_1L_EWSB_in_0L_mass = atoi(argv[7]);
    use_Goldstone_resum = atoi(argv[8]);

  } else {
    std::cout << "Use ./run_xSM_covaraint_gaugue ms lambda_s lambda_hs" << std::endl;
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
  auto model = EffectivePotential::xSM_MSbar_covariant::from_tadpoles(lambda_hs, lambda_s, ms, Q, xi, tree_level_tadpoles);
  model.set_use_1L_EWSB_in_0L_mass(use_1L_EWSB_in_0L_mass);

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

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
  pf.set_check_vacuum_at_high(false);
  pf.set_seed(1);
  
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
  // Print the data in a particular format for plotting
//  if (debug_mode) PhaseTracer::phase_plotter(tf, "xSM_MSbar");
  return 0;
  
}
