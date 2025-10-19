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
#include <interpolation.h>

#include "DP_xSM.hpp"
#include "models/xSM_HT.hpp" 
#include "phasetracer.hpp"
#include "helperIncludes/gwTools.hpp"
#include "helperIncludes/wallGoWrapper.hpp"
#include "deep_phase.hpp"

using json = nlohmann::json;

void write_spectra(std::vector<double> fit_amp, std::vector<double> fit_freq, std::vector<double> ssm_amp, std::vector<double> ssm_freq, std::string path = "spec.csv") {
  std::ofstream output_file(path);
  if (!output_file.is_open()) {
    std::cerr << "Error: Could not open file " << path << " for writing" << std::endl;
    return;
  }

  size_t fit_length = fit_amp.size();
  size_t ssm_length = fit_amp.size();
  if ( fit_length != ssm_length ) {
    std::cerr << "Sizes incompatible" << std::endl;
  }
  output_file << "fit_freq,fit_amp,ssm_freq,ssm_amp,\n";
  for ( size_t i = 0; i < fit_length; i += 1) {
    output_file << fit_freq[i] << "," << fit_amp[i] << "," << ssm_freq[i] << "," << ssm_amp[i] << "\n";
  }
  output_file.close();
}

std::vector<PhaseTracer::Transition>
validate_transitions(std::vector<PhaseTracer::Transition> input)
{
  std::vector<PhaseTracer::Transition> output;
  for (auto t : input) 
  {
    // vac = [h, s]
    auto true_vac = t.true_vacuum;
    auto false_vac = t.false_vacuum;
    auto changed = t.changed;

    if(changed[0] && changed[1])
    {
      double h_true = true_vac[0];
      double h_false = false_vac[0];
      double s_true = true_vac[1];
      double s_false = false_vac[1];
      if( (h_true > 5. && abs(s_true) < 1e-3) && (s_false > 5. && abs(h_false) < 1e-3) )
      {
        output.push_back(t);
      }     
    }
  }
  return output;
}

std::string
format_thermal_params(std::vector<double> in, int flag, double vw, double vwLTE, PhaseTracer::ThermalParams tp)
{
  std::string output = "";

  for ( double val : in ) { output += std::to_string(val) + ","; }

  output += std::to_string(flag) + ",";
  double Tf = tp.completes == PhaseTracer::MilestoneStatus::YES ? tp.TF : -1;
  output += std::to_string(vw) + ",";
  output += std::to_string(vwLTE) + ",";
  output += std::to_string(Tf) + ",";

  if(tp.percolates == PhaseTracer::MilestoneStatus::YES) {
    output += std::to_string(tp.TP) + ","
    + std::to_string(tp.alpha_tp) + ","
    + std::to_string(tp.betaH_tp) + ","
    + std::to_string(tp.we_tp) + ","
    + std::to_string(tp.cs_true_tp) + ","
    + std::to_string(tp.cs_false_tp) + ",";
  } else {
    output += "-1,0,0,0,0,0,";
  }

  if(tp.nucleates == PhaseTracer::MilestoneStatus::YES) {
    output += std::to_string(tp.TN) + ","
    + std::to_string(tp.alpha_tn) + ","
    + std::to_string(tp.betaH_tn) + ","
    + std::to_string(tp.we_tn) + ","
    + std::to_string(tp.cs_true_tn) + ","
    + std::to_string(tp.cs_false_tn) + ",";
  } else {
    output += "-1,0,0,0,0,0,";
  }

  if ( !output.empty() ) { output.pop_back(); }
  output += "\n";

  return output;
}

std::string
format_fail_case(std::vector<double> in, int flag)
{
  std::string output = "";

  for ( double val : in ) { output += std::to_string(val) + ","; }

  output += std::to_string(flag) + ",";

  output += "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0\n";

  return output;
}

/*
  Main
*/
int main(int argc, char* argv[]) {

  // default output
  std::string output_filename = "z_eos_scan/ht_output.csv";
  bool file_has_data = file_exists_and_not_empty(output_filename);

  std::ofstream output_file(output_filename, std::ios::app);
  if (!file_has_data) {
    output_file << "mS,lambda_s,lambda_hs,";
    output_file << "success,vw,vwLTE,Tf,";
    output_file << "Tp,alpha_tp,betaH_tp,we_tp,cs_true_tp,cs_false_tp,";
    output_file << "Tn,alpha_tn,betaH_tn,we_tn,cs_true_tn,cs_false_tn\n";
  }

  double vw, vwLTE, dtauRs;

  double mS, lambda_s, lambda_hs;
  double Q, xi;
  double daisy_flag = 2;
  bool use_1L_EWSB_in_0L_mass = false;;
  bool use_Goldstone_resum = true;
  bool tree_level_tadpoles = false;
  bool use_covariant_gauge = false;
  std::vector<double> SM_parameters ={};

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
    vw = modelParams["vw"].get<double>();
    dtauRs = modelParams["dtauRs"].get<double>();
    mS = modelParams["mS"].get<double>();
    lambda_hs = modelParams["lambda_hs"].get<double>();
    lambda_s = modelParams["lambda_s"].get<double>();
  } catch (...) {
    std::cerr << "Input file might not be formatted correctly! Quitting now." << std::endl;
    // check find_phases suceeded
    output_file << format_fail_case({-1.,-1.,-1.}, -10);
    output_file.close();
    return 1;
  }

  std::vector<double> in = {mS, lambda_s, lambda_hs};
  EffectivePotential::xSM_HT model(lambda_hs, lambda_s, mS);

  PhaseTracer::PhaseFinder pf(model);
  pf.set_seed(0);
  pf.set_check_hessian_singular(true);
  pf.set_check_vacuum_at_high(false);
  pf.set_t_high(1000);
  pf.set_t_low(10);
  pf.set_dt_max_rel(1.);
  pf.set_dt_max_abs(2.5);

  try {
    pf.find_phases();
  } catch (const std::exception& e) {
    std::cerr << "Exception in pf.find_phases(): " << e.what() << std::endl;
    output_file << format_fail_case(in, -1);
    output_file.close();
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception in pf.find_phases()" << std::endl;
    output_file << format_fail_case(in, -1);
    output_file.close();
    return 1;
  }
  auto p1 = pf.get_phases();

  // check for phases
  if (p1.size() == 0) {
    output_file << format_fail_case(in, -2);
    output_file.close();
    return 1;
  }

  std::cout << pf;

  PhaseTracer::ActionCalculator ac(model);
  PhaseTracer::TransitionFinder tf(pf, ac);
  tf.set_calculate_action(false); // This initialises with ac but stops it from calculating TN.
  try {
    tf.find_transitions();
  } catch (const std::exception& e) {
    std::cerr << "Exception in tf.find_transitions(): " << e.what() << std::endl;
    output_file << format_fail_case(in, -3);
    output_file.close();
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception in tf.find_transitions()" << std::endl;
    output_file << format_fail_case(in, -3);
    output_file.close();
    return 1;
  }
  std::cout << tf;

  auto t = tf.get_transitions();
  // check for transitions
  if (t.size() == 0) {
    output_file << format_fail_case(in, -4);
    output_file.close();
    return 1;
  }

  auto valid_t = validate_transitions(t);
  if (valid_t.size() == 0) {
    std::cout << "No valid transitions for required phases." << std::endl;
    output_file << format_fail_case(in, -5);
    output_file.close();
    return 1;
  }

  PhaseTracer::ThermalParameters tp(tf);
  tp.set_compute_debug(false);
  tp.set_vw(vw);
  tp.set_dof(107.75);

  PhaseTracer::ThermalParams tps;
  try {
    tps = tp.process_transition(valid_t[0]);
  } catch (...) {
    std::cout << "Error calculating thermal params." << std::endl;
    output_file << format_fail_case(in, -6);
    output_file.close();
    return 1;
  }

  if(tps.percolates != PhaseTracer::MilestoneStatus::YES && tps.nucleates != PhaseTracer::MilestoneStatus::YES) {
    std::cout << "No vaid reference temperature." << std::endl;
    output_file << format_fail_case(in, -7);
    output_file.close();
    return 1;
  }

  // optional, use wallgo to calculate vw and replace
  double tperc = tps.TP;
  const auto params4d = model.get_4d_parameter_map(tperc);
  auto wallVelocities = std::vector<double>{};
  try {
    wallVelocities = wallGoWrapper::getWallVelocity(valid_t[0], params4d, tperc, "ht_xSM");
  } catch (...) { wallVelocities = wallVelocities;}
  if ( wallVelocities.size() == 2) {
    vw = wallVelocities[0];
    vwLTE = wallVelocities[1];
  } else {
    vw = vw;
    vwLTE = -1;
  }

  std::cout << "========================\n";
  std::cout << "Tref = " << tps.TP << "\n";
  std::cout << "alpha = " << tps.alpha_tp << "\n";
  std::cout << "betaH = " << tps.betaH_tp << "\n";
  std::cout << "H = " << tps.H_tp << "\n";
  std::cout << "w/e = " << tps.we_tp << "\n";
  std::cout << "vw = " << vw << "\n";
  std::cout << "cs_true (tp) = " << tps.cs_true_tp << "\n";
  std::cout << "cs_false (tp) = " << tps.cs_false_tp << "\n";

  // const auto eos_path = "z_eos_scan/eos/ht_l_hs_" + std::to_string(lambda_hs) + "_l_s_" + std::to_string(lambda_s) + "_m_s_" + std::to_string(mS) + ".csv";
  // const auto eos_path = "z_eos_scan/ht_bp1.csv";
  // tps.eos.write_EoS(eos_path);

  output_file << format_thermal_params(in, 1, vw, vwLTE, tps);
  output_file.close();

    // try deep phase fluid profiles

  auto eos_dp = PhaseTracer::pt_EoS_to_dp_EoS(tps.eos);
  auto universe = PhaseTransition::Universe(tps.TP, 107.75, tps.H_tp);

  double RsH = pow(8.*M_PI, 1./3.) * vw / tps.betaH_tp;
  double dtau = RsH * dtauRs / tps.H_tp;

  auto pt_params_bag = PhaseTransition::PTParams_Bag(vw, tps.alpha_tp, tps.TP, tps.betaH_tp * tps.H_tp, dtau, "exp", universe, tps.cs_true_tp, tps.cs_false_tp);

  const Hydrodynamics::FluidProfile profile_bag(pt_params_bag);
  profile_bag.write("ht_fp_bag.csv");

  auto pt_params_veff = PhaseTransition::PTParams_Veff(vw, tps.alpha_tp, tps.TP, tps.betaH_tp * tps.H_tp, dtau, "exp", universe, eos_dp);

  const Hydrodynamics::FluidProfile profile_veff(pt_params_veff);
  profile_veff.write("ht_fp_veff.csv");

  const auto kRs_vals = logspace(1e-3, 1e+3, 100);
  Spectrum::PowerSpec OmegaGW_bag = Spectrum::GWSpec2(kRs_vals, pt_params_bag);
  OmegaGW_bag.write("ht_spec_bag.csv");
  Spectrum::PowerSpec OmegaGW_veff = Spectrum::GWSpec2(kRs_vals, pt_params_veff);
  OmegaGW_veff.write("ht_spec_veff.csv");
  return 0;
}