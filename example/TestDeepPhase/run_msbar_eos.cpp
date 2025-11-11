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
#include "models/xSM_MSbar.hpp"
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
      if( (h_true > 10. && abs(s_true) < 1e-3) && (s_false > 10. && abs(h_false) < 1e-3) )
      {
        output.push_back(t);
      }     
    }
  }
  return output;
}

std::string
format_thermal_params(std::vector<double> in, int flag, double TC, PhaseTracer::TransitionMilestone milestone, wallGoWrapper::wallGoResults vwResults)
{
  std::string output = "";

  for ( double val : in ) { output += std::to_string(val) + ","; }

  output += std::to_string(flag) + ",";
  double Tc = TC;
  output += std::to_string(Tc) + ",";

  output += std::to_string(vwResults.vwLTE) + ","
  + std::to_string(vwResults.vJ) + ","
  + std::to_string(vwResults.vw) + ","
  + std::to_string(vwResults.vw_flag) + ","
  + std::to_string(vwResults.vw_det) + ","
  + std::to_string(vwResults.vw_det_flag) + ",";

  if(milestone.status == PhaseTracer::MilestoneStatus::YES) {
    output += std::to_string(milestone.temperature) + ","
    + std::to_string(milestone.alpha) + ","
    + std::to_string(milestone.betaH) + ","
    + std::to_string(milestone.Rs) + ","
    + std::to_string(milestone.we) + ","
    + std::to_string(milestone.cs_plus) + ","
    + std::to_string(milestone.cs_minus) + ",";
  } else {
    output += "0,0,0,0,0,0,0,";
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

  output += "0,0,0,0,0,0,0,0,0,0,0,0,0,0\n";

  return output;
}

std::string 
format_debug_fail_case(std::vector<double> in, int flag) 
{
  std::string output = "";

  for ( double val : in ) { output += std::to_string(val) + ","; }

  output += std::to_string(flag) + ",";

  output += "0,0,0,0,0,0,0,0,0,0\n";

  return output;
}

/*
  Main
*/
int main(int argc, char* argv[]) {

  // default output
  std::string output_filename = "comparisonData/parameter_scans/data/msbar.csv";
  bool file_has_data = file_exists_and_not_empty(output_filename);

  std::ofstream output_file(output_filename, std::ios::app);
  if (!file_has_data) {
    output_file << "mS,lambda_s,lambda_hs,";
    output_file << "success,Tc,";
    output_file << "vwLTE,vJ,vw,vwFlag,vwDet,vwDetFlag,";
    output_file << "T,alpha,betaH,RsH,we,cs_plus,cs_minus\n";
  }

  double vw, vwLTE, vwDet, vJ, dtauRs;

  double mS, lambda_s, lambda_hs;
  double Q, xi;
  double daisy_flag = 2;
  bool use_1L_EWSB_in_0L_mass = false;
  bool use_Goldstone_resum = true;
  bool tree_level_tadpoles = false;
  bool use_covariant_gauge = false;
  std::vector<double> SM_parameters ={};

  LOGGER(debug);

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
    Q = modelParams["Q"].get<double>();
    xi = modelParams["xi"].get<double>();
    daisy_flag = modelParams["daisy"].get<int>();
    tree_level_tadpoles = modelParams["tree_level_tadpoles"].get<bool>();
  } catch (...) {
    std::cerr << "Input file might not be formatted correctly! Quitting now." << std::endl;
    // check find_phases suceeded
    output_file << format_fail_case({-1.,-1.,-1.}, -10);
    output_file.close();
    return 1;
  }

  std::vector<double> in = {mS, lambda_s, lambda_hs};
  auto model = EffectivePotential::xSM_MSbar::from_tadpoles(lambda_hs, lambda_s, mS, Q, xi, use_covariant_gauge, use_1L_EWSB_in_0L_mass, use_Goldstone_resum, tree_level_tadpoles, SM_parameters);
  if (daisy_flag == 0){
    model.set_daisy_method(EffectivePotential::DaisyMethod::None);
  } else if (daisy_flag == 1){
    model.set_daisy_method(EffectivePotential::DaisyMethod::Parwani);
  } else if (daisy_flag == 2){
    model.set_daisy_method(EffectivePotential::DaisyMethod::ArnoldEspinosa);
  } else {
    std::cout << "Wrong daisy flag" << std::endl;
  }

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

  std::cout << valid_t[0];

  PhaseTracer::ThermoFinder tm(tf, ac);
  tm.set_vw(vw);
  tm.set_dof(107.75);
  tm.set_background_dof(107.75 - 38.5);
  tm.set_percolation_print_setting(PhaseTracer::PrintSettings::VERBOSE);
  tm.set_nucleation_print_setting(PhaseTracer::PrintSettings::VERBOSE);

  PhaseTracer::TransitionMilestone milestone;
  std::optional<PhaseTracer::EquationOfState> eos_opt;

  try {
    auto tps = tm.get_thermal_parameter_set(valid_t[0]);

    milestone = tps.percolation;
    eos_opt = std::move(tps.eos);
  } catch(...) {
    std::cout << "Error during calculation of thermal parameters." << std::endl;
    output_file << format_fail_case(in, -6);
    output_file.close();
    return 1;
  }
  
  if(milestone.status != PhaseTracer::MilestoneStatus::YES) {
    std::cout << "Nucleation milestone not achieved." << std::endl;
    output_file << format_fail_case(in, -7);
    output_file.close();
    return 1;
  }

  std::cout << "============================" << "\n";
  std::cout << milestone;

  // optional, use wallgo to calculate vw and replace
  double tnuc = milestone.temperature;
  const auto params4d = model.get_4d_parameter_map(tnuc);
  const auto wallGoOut = wallGoWrapper::getWallVelocity(valid_t[0], params4d, tnuc, "msbar_xSM");
  vwLTE = wallGoOut.vwLTE;
  vw = (wallGoOut.vw_flag == 1) ? wallGoOut.vw : vw;
  vJ = (wallGoOut.vw_flag == 1) ? wallGoOut.vJ : vw;
  vwDet = wallGoOut.vw_det;

  std::cout << "============================" << "\n";
  std::cout << "  vJ = " << vJ << "\n";
  std::cout << "  vw (LTE) = " << vwLTE << "\n";
  std::cout << "  vw (deflagration) = " << vw << "\n";
  std::cout << "  vw (detonation) = " << vwDet << "\n";
  std::cout << "============================" << "\n";

  output_file << format_thermal_params(in, 1, valid_t[0].TC, milestone, wallGoOut);
  output_file.close();
  return 0;

  PhaseTracer::GravWaveCalculator gc(tf);
  gc.set_min_frequency(1e-5);
  gc.set_max_frequency(1e0);
  gc.set_num_frequency(500);
  gc.set_vw(vw);

  auto gw = gc.calc_spectrum(milestone.alpha, milestone.betaH, milestone.temperature);
  gc.write_spectrum_to_text(gw, "comparisonData/spectra/msbar/naive.csv");

  namespace DPInt = PhaseTracer::DeepPhaseInterface;

  eos_opt->write("comparisonData/eos/msbar_eos.csv");

  if(vw > 0)
  {
    auto deflagaration_bag = DPInt::DeepPhaseResults(milestone, *eos_opt, 107.75, vw, dtauRs, DPInt::EoSModel::BAG);
    auto def_bag_profile = deflagaration_bag.get_profile_bag();
    def_bag_profile.write("comparisonData/profiles/msbar/def_bag.csv");
    auto def_bag_spectrum = deflagaration_bag.get_spectrum_bag();
    def_bag_spectrum.write("comparisonData/spectra/msbar/def_bag.csv");

    auto deflagaration_munu = DPInt::DeepPhaseResults(milestone, *eos_opt, 107.75, vw, dtauRs, DPInt::EoSModel::MUNU);
    auto def_munu_profile = deflagaration_munu.get_profile_munu();
    def_munu_profile.write("comparisonData/profiles/msbar/def_munu.csv");
    auto def_munu_spectrum = deflagaration_munu.get_spectrum_munu();
    def_munu_spectrum.write("comparisonData/spectra/msbar/def_munu.csv");

    auto deflagaration_veff = DPInt::DeepPhaseResults(milestone, *eos_opt, 107.75, vw, dtauRs, DPInt::EoSModel::VEFF);
    auto def_veff_profile = deflagaration_veff.get_profile_veff();
    def_veff_profile.write("comparisonData/profiles/msbar/def_veff.csv");
    auto def_veff_spectrum = deflagaration_veff.get_spectrum_veff();
    def_veff_spectrum.write("comparisonData/spectra/msbar/def_veff.csv");
  }
  
  if(vwDet > 0)
  {
    auto detonation_bag = DPInt::DeepPhaseResults(milestone, *eos_opt, 107.75, vwDet, dtauRs, DPInt::EoSModel::BAG);
    auto det_bag_profile = detonation_bag.get_profile_bag();
    det_bag_profile.write("comparisonData/profiles/msbar/det_bag.csv");
    auto det_bag_spectrum = detonation_bag.get_spectrum_bag();
    det_bag_spectrum.write("comparisonData/spectra/msbar/det_bag.csv");

    auto detonation_munu = DPInt::DeepPhaseResults(milestone, *eos_opt, 107.75, vwDet, dtauRs, DPInt::EoSModel::MUNU);
    auto det_munu_profile = detonation_munu.get_profile_munu();
    det_munu_profile.write("comparisonData/profiles/msbar/det_munu.csv");
    auto det_munu_spectrum = detonation_munu.get_spectrum_munu();
    det_munu_spectrum.write("comparisonData/spectra/msbar/det_munu.csv");

    auto detonation_veff = DPInt::DeepPhaseResults(milestone, *eos_opt, 107.75, vwDet, dtauRs, DPInt::EoSModel::VEFF);
    auto det_veff_profile = detonation_veff.get_profile_veff();
    det_veff_profile.write("comparisonData/profiles/msbar/det_veff.csv");
    auto det_veff_spectrum = detonation_veff.get_spectrum_veff();
    det_veff_spectrum.write("comparisonData/spectra/msbar/det_veff.csv");
  }

  return 0;
}