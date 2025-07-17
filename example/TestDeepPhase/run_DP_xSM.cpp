#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include <sys/stat.h>
#include <filesystem>
#include <sstream>
#include <chrono>
#include <omp.h>
#include <array>
#include <memory>
#include <unistd.h>

#include "dataanalysis.h"

#include "DP_xSM.hpp"
#include "phasetracer.hpp"
#include "DeepPhase/include/deepphase.hpp"

using json = nlohmann::json;

/*
  Main
*/
int main(int argc, char* argv[]) {

  // default output
  std::string output_filename = "example/Comparison/results_xsm.csv";
  bool file_has_data = file_exists_and_not_empty(output_filename);

  std::ofstream output_file(output_filename, std::ios::app);
  if (!file_has_data) {
    // output_file << "mS,lambda_s,lambda_hs,alpha,beta_H,Tref,vw,cs_true,cs_false\n";
    output_file << "mS,lambda_s,lambda_hs,vw,Tn,alpha_tn,beta_tn,cs_true_tn,cs_false_tn,Tp,alpha_tp,beta_tp,cs_true_tp,cs_false_tp\n";
  }

  double Ms, lambda_s, lambda_hs;
  bool running = true;
  int potential = 1;
  int matching = 1;

  LOGGER(fatal);

  std::string json_filename;
  if ( argc == 1 ) {
    json_filename = "example/TestDeepPhase/modelParams_xSM.json";
  } else {
    json_filename = argv[1];
    std::cout << json_filename << std::endl;
  }
  json modelParams = readFile(json_filename);

  try {
    Ms = modelParams["mS"].get<double>();
    lambda_s = modelParams["lambda_s"].get<double>();
    lambda_hs = modelParams["lambda_hs"].get<double>();
  } catch (...) {
    std::cerr << "Input file might not be formatted correctly! Quitting now." << std::endl;
    // check find_phases suceeded
    output_file << formatOutput({-1.,-1.,-1.}, -10);
    output_file.close();
    return 1;
  }

  std::vector<double> in = {Ms, lambda_s, lambda_hs};
  EffectivePotential::DR_xsm model(Ms, lambda_hs, lambda_s);
  model.set_matching_flag(matching);
  model.set_potential_flag(potential);
  model.set_running_flag(running);

  PhaseTracer::PhaseFinder pf(model);
  pf.set_seed(0);
  pf.set_check_hessian_singular(false);
  pf.set_check_vacuum_at_high(false);
  pf.set_t_high(1000);
  pf.set_t_low(15);
  pf.set_dt_max_rel(1.); // makes dt_max = rel * (t_high - t_low)
  pf.set_dt_max_abs(2.5); // get a finer resolution for splines

  try {
    pf.find_phases();
  } catch (...) {
    // check find_phases suceeded
    output_file << formatOutput(in, -1);
    output_file.close();
    return 1;
  }
  auto p1 = pf.get_phases();

  // check for phases
  if (p1.size() == 0) {
    output_file << formatOutput(in, -2);
    output_file.close();
    return 1;
  }

  PhaseTracer::ActionCalculator ac(model);
  PhaseTracer::TransitionFinder tf(pf, ac);
  tf.set_calculate_action(false); // This initialises with ac but stops it from calculating TN.
  tf.find_transitions();

  std::cout << pf;
  std::cout << tf;

  auto t = tf.get_transitions();
  // check for transitions
  if (t.size() == 0) {
    output_file << formatOutput(in, -3);
    output_file.close();
    return 1;
  }
  

  // extract the correct transition (probably unnecessary)
  int flag = -1;
  for ( auto trans : t) {
    if ( 
      (abs(trans.false_vacuum[0]) < 1e-4 && abs(trans.false_vacuum[1]) > 5.) &&
      (abs(trans.true_vacuum[0]) > 10. &&  abs(trans.true_vacuum[1]) < 1e-4 ) ) {
        flag = trans.key;
      }
  }
  // should place a check for flag = -1 here.
  if ( flag == -1 ) {
    output_file << formatOutput(in, -4);
    output_file.close();
    return 1;
  }

  // testing
  PhaseTracer::ThermalParameters tp(tf);
  tp.find_thermal_parameters();
  std::cout << tp;

  auto tps = tp.get_thermal_params();

  double vw = 0.623436;
  double alpha = tps[0].alpha_tp;
  double beta = tps[0].beta_tp;
  double cs_true_tp = 1/sqrt(3.);
  double cs_false_tp = 1/sqrt(3.);

  // interface with deepphase
  double dtau = 1/beta;
  double wN = 1.71;
  std::string PTmodel = "bag";
  std::string nuc_type = "exp";
  PhaseTransition::Universe un;

  PhaseTransition::PTParams paramsPT(vw, alpha, beta, dtau, wN, PTmodel, nuc_type, un);
  std::vector<double> momentumVec = logspace(1e-3, 1e+3, 100);
  Spectrum::PowerSpec OmegaGW = Spectrum::GWSpec(momentumVec, paramsPT);
  OmegaGW.write("gw_spec.csv");

  PhaseTracer::GravWaveCalculator gc(tf);
  gc.set_min_frequency(1e-3);
  gc.set_max_frequency(1e+3);
  gc.set_num_frequency(100);

  PhaseTracer::GravWaveSpectrum gw = gc.calc_spectrum(alpha, tps[0].betaH_tp, tps[0].TP);
  gc.write_spectrum_to_text(gw, "gw_spec_pt.csv");

  return 0;
}