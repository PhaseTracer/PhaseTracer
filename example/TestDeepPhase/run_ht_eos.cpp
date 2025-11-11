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
#include <nlopt.hpp>

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
    + std::to_string(milestone.dt) + ","
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

  output += "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0\n";

  return output;
}

/*
  Main
*/
int main(int argc, char* argv[]) {

  // default output
  std::string output_filename = "comparisonData/parameter_scans/data/ht.csv";
  bool file_has_data = file_exists_and_not_empty(output_filename);

  std::ofstream output_file(output_filename, std::ios::app);
  if (!file_has_data) {
    output_file << "mS,lambda_s,lambda_hs,";
    output_file << "success,Tc,";
    output_file << "vwLTE,vJ,vw,vwFlag,vwDet,vwDetFlag,";
    output_file << "T,alpha,betaH,RsH,dtH,we,cs_plus,cs_minus\n";
  }

  double vw, vwDet, vJ, vwLTE, dtauRs;

  double mS, lambda_s, lambda_hs;
  double Q, xi;
  double daisy_flag = 2;
  bool use_1L_EWSB_in_0L_mass = false;;
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

  PhaseTracer::ActionCalculator ac(pf);

  PhaseTracer::TransitionFinder tf(pf);
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

  // std::cout << valid_t[0];

  // const auto true_phase = valid_t[0].true_phase;
  // const auto false_phase = valid_t[0].false_phase;

  // // ========== CONSISTENT TEMPERATURE GRID MINIMIZATION ==========
  // std::cout << "\n========== Finding Minima on Consistent Temperature Grid ==========\n";
  
  // // Define temperature range
  // double t_min = valid_t[0].false_phase.T.front();
  // double t_max = valid_t[0].TC;
  // int num_temps = 50;  // Number of temperature points
  
  // std::cout << "Temperature range: " << t_min << " to " << t_max << "\n";
  // std::cout << "Number of temperature points: " << num_temps << "\n\n";
  
  // // Wrapper structure to pass model and temperature to NLopt
  // struct PotentialData {
  //   const EffectivePotential::xSM_HT* model;
  //   double T;
  // };
  
  // // Objective function for NLopt (must match nlopt signature)
  // auto objective_func = [](const std::vector<double>& x, std::vector<double>& grad, void* data) -> double {
  //   PotentialData* pd = static_cast<PotentialData*>(data);
  //   Eigen::VectorXd phi(2);
  //   phi << x[0], x[1];
    
  //   // Optionally compute gradient numerically if needed
  //   if (!grad.empty()) {
  //     double eps = 1e-6;
  //     for (size_t i = 0; i < 2; ++i) {
  //       Eigen::VectorXd phi_plus = phi;
  //       Eigen::VectorXd phi_minus = phi;
  //       phi_plus[i] += eps;
  //       phi_minus[i] -= eps;
  //       grad[i] = (pd->model->V(phi_plus, pd->T) - pd->model->V(phi_minus, pd->T)) / (2.0 * eps);
  //     }
  //   }
    
  //   return pd->model->V(phi, pd->T);
  // };
  
  // // Helper function to find initial guess by interpolation
  // auto get_interpolated_guess = [](double T_target, const PhaseTracer::Phase& phase) -> Eigen::VectorXd {
  //   const auto& T_list = phase.T;
  //   const auto& X_list = phase.X;
    
  //   // Find the two nearest temperatures that bracket T_target
  //   size_t n = T_list.size();
    
  //   // Handle edge cases
  //   if (T_target <= T_list[0]) return X_list[0];
  //   if (T_target >= T_list[n-1]) return X_list[n-1];
    
  //   // Find bracketing indices
  //   size_t i_upper = 0;
  //   for (size_t i = 0; i < n; ++i) {
  //     if (T_list[i] >= T_target) {
  //       i_upper = i;
  //       break;
  //     }
  //   }
  //   size_t i_lower = (i_upper > 0) ? i_upper - 1 : 0;
    
  //   // Linear interpolation
  //   double T_lower = T_list[i_lower];
  //   double T_upper = T_list[i_upper];
  //   double alpha = (T_target - T_lower) / (T_upper - T_lower + 1e-15);  // avoid division by zero
    
  //   Eigen::VectorXd guess = (1.0 - alpha) * X_list[i_lower] + alpha * X_list[i_upper];
  //   return guess;
  // };
  
  // // Helper function to find minimum using NLopt
  // auto find_minimum = [&](double T, const Eigen::VectorXd& initial_guess) -> std::optional<Eigen::VectorXd> {
  //   nlopt::opt opt(nlopt::LN_COBYLA, 2);
    
  //   PotentialData pd{&model, T};
  //   opt.set_min_objective(objective_func, &pd);
  //   opt.set_xtol_rel(1e-8);
  //   opt.set_maxeval(1000);
    
  //   std::vector<double> x = {initial_guess[0], initial_guess[1]};
    
  //   try {
  //     double minf;
  //     nlopt::result result = opt.optimize(x, minf);
      
  //     if (result > 0) {  // Positive results indicate success
  //       Eigen::VectorXd minimum(2);
  //       minimum << x[0], x[1];
  //       return minimum;
  //     }
  //   } catch (...) {
  //     return std::nullopt;
  //   }
    
  //   return std::nullopt;
  // };
  
  // // Create consistent temperature grid
  // std::vector<double> T_grid;
  // for (int i = 0; i < num_temps; ++i) {
  //   double T = t_max - (t_max - t_min) * i / (num_temps - 1.0);
  //   T_grid.push_back(T);
  // }
  
  // // Storage for results
  // std::vector<double> T_consistent;
  // std::vector<Eigen::VectorXd> true_minima;
  // std::vector<Eigen::VectorXd> false_minima;
  
  // // Find minima for both phases at each temperature
  // std::cout << "Finding minima for both phases...\n";
  // int success_count = 0;
  // int failed_count = 0;
  
  // for (size_t i = 0; i < T_grid.size(); ++i) {
  //   double T = T_grid[i];
    
  //   // Get initial guesses via interpolation
  //   Eigen::VectorXd true_guess = get_interpolated_guess(T, true_phase);
  //   Eigen::VectorXd false_guess = get_interpolated_guess(T, false_phase);
    
  //   // Find minima
  //   auto true_min_opt = find_minimum(T, true_guess);
  //   auto false_min_opt = find_minimum(T, false_guess);
    
  //   if (true_min_opt && false_min_opt) {
  //     T_consistent.push_back(T);
  //     true_minima.push_back(*true_min_opt);
  //     false_minima.push_back(*false_min_opt);
  //     success_count++;
      
  //     if (i % 10 == 0 || i < 3 || i >= T_grid.size() - 3) {
  //       std::cout << "T = " << T << ": ✓\n";
  //       std::cout << "  True:  [" << (*true_min_opt)[0] << ", " << (*true_min_opt)[1] << "]\n";
  //       std::cout << "  False: [" << (*false_min_opt)[0] << ", " << (*false_min_opt)[1] << "]\n";
  //     }
  //   } else {
  //     failed_count++;
  //     std::cout << "T = " << T << ": ✗ Failed to find minimum\n";
  //   }
  // }
  
  // std::cout << "\n========== Consistent Grid Results ==========\n";
  // std::cout << "Requested points: " << num_temps << "\n";
  // std::cout << "Successful: " << success_count << "\n";
  // std::cout << "Failed: " << failed_count << "\n";
  // std::cout << "============================================\n\n";
  
  // // Optional: Print all results to a file
  // if (success_count > 0) {
  //   std::ofstream minima_file("consistent_minima.csv");
  //   minima_file << "T,true_phi0,true_phi1,false_phi0,false_phi1,V_true,V_false,delta_V\n";
  //   for (size_t i = 0; i < T_consistent.size(); ++i) {
  //     double T = T_consistent[i];
  //     double V_true = model.V(true_minima[i], T);
  //     double V_false = model.V(false_minima[i], T);
  //     minima_file << T << ","
  //                 << true_minima[i][0] << "," << true_minima[i][1] << ","
  //                 << false_minima[i][0] << "," << false_minima[i][1] << ","
  //                 << V_true << "," << V_false << "," << (V_false - V_true) << "\n";
  //   }
  //   minima_file.close();
  //   std::cout << "Results written to: consistent_minima.csv\n\n";
  // }

  // return 0;

  // // ========== OLD TEST CODE (keeping for reference) ==========
  // std::cout << "\n========== Testing NLopt Minimization ==========\n";
  // // Test ALL temperature points from the true_phase
  // size_t num_tests = true_phase.T.size();
  // size_t num_success = 0;
  // size_t num_failed = 0;
  // double max_diff_0 = 0.0;
  // double max_diff_1 = 0.0;
  // double max_diff_V = 0.0;
  
  // std::cout << "Testing " << num_tests << " temperature points...\n\n";
  
  // for (size_t idx = 0; idx < num_tests; ++idx) {
  //   double T = true_phase.T[idx];
  //   Eigen::VectorXd phi_known = true_phase.X[idx];
    
  //   // Print progress every 10 points
  //   if (idx % 10 == 0 || idx < 3 || idx >= num_tests - 3) {
  //     std::cout << "Test " << idx + 1 << "/" << num_tests 
  //               << " (T = " << T << "): ";
  //   }
    
  //   // Setup NLopt
  //   nlopt::opt opt(nlopt::LN_COBYLA, 2);  // 2D optimization, gradient-free COBYLA
    
  //   PotentialData pd{&model, T};
  //   opt.set_min_objective(objective_func, &pd);
  //   opt.set_xtol_rel(1e-8);
  //   opt.set_maxeval(1000);
    
  //   // Use known minimum as initial guess
  //   std::vector<double> x = {phi_known[0], phi_known[1]};
    
  //   try {
  //     double minf;
  //     nlopt::result result = opt.optimize(x, minf);
      
  //     // Compare with known minimum
  //     double diff_0 = std::abs(x[0] - phi_known[0]);
  //     double diff_1 = std::abs(x[1] - phi_known[1]);
  //     double diff_V = std::abs(minf - model.V(phi_known, T));
      
  //     // Track maximum differences
  //     max_diff_0 = std::max(max_diff_0, diff_0);
  //     max_diff_1 = std::max(max_diff_1, diff_1);
  //     max_diff_V = std::max(max_diff_V, diff_V);
      
  //     if (diff_0 < 1e-4 && diff_1 < 1e-4) {
  //       num_success++;
  //       if (idx % 10 == 0 || idx < 3 || idx >= num_tests - 3) {
  //         std::cout << "✓ (diffs: " << diff_0 << ", " << diff_1 << ")\n";
  //       }
  //     } else {
  //       num_failed++;
  //       std::cout << "\n✗ FAILED at T = " << T << "\n";
  //       std::cout << "  Known minimum: [" << phi_known[0] << ", " << phi_known[1] << "]\n";
  //       std::cout << "  Found minimum: [" << x[0] << ", " << x[1] << "]\n";
  //       std::cout << "  Difference in phi[0]: " << diff_0 << "\n";
  //       std::cout << "  Difference in phi[1]: " << diff_1 << "\n";
  //       std::cout << "  Difference in V: " << diff_V << "\n";
  //     }
      
  //   } catch (std::exception& e) {
  //     num_failed++;
  //     std::cerr << "\n✗ NLopt exception at T = " << T << ": " << e.what() << "\n";
  //   }
  // }
  
  // std::cout << "\n========== NLopt Test Summary ==========\n";
  // std::cout << "Total tests: " << num_tests << "\n";
  // std::cout << "Successful: " << num_success << " (" 
  //           << (100.0 * num_success / num_tests) << "%)\n";
  // std::cout << "Failed: " << num_failed << "\n";
  // std::cout << "Max difference in phi[0]: " << max_diff_0 << "\n";
  // std::cout << "Max difference in phi[1]: " << max_diff_1 << "\n";
  // std::cout << "Max difference in V: " << max_diff_V << "\n";
  // std::cout << "========================================\n\n";

  // return 0;

  PhaseTracer::ThermoFinder tm(ac);
  tm.set_vw(vw);
  tm.set_dof(107.75);
  tm.set_background_dof(107.75);
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
    std::cout << "Milestone not achieved." << std::endl;
    output_file << format_fail_case(in, -7);
    output_file.close();
    return 1;
  }

  std::cout << "============================" << "\n";
  std::cout << milestone;

  // optional, use wallgo to calculate vw and replace
  double tnuc = milestone.temperature;
  const auto params4d = model.get_4d_parameter_map(tnuc);
  const auto wallGoOut = wallGoWrapper::getWallVelocity(valid_t[0], params4d, tnuc, "ht_xSM");
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
  gc.write_spectrum_to_text(gw, "comparisonData/spectra/ht/naive.csv");

  namespace DPInt = PhaseTracer::DeepPhaseInterface;

  eos_opt->write("comparisonData/eos/ht_eos.csv");

  if(vw > 0)
  {
    auto deflagaration_bag = DPInt::DeepPhaseResults(milestone, *eos_opt, 107.75, vw, dtauRs, DPInt::EoSModel::BAG);
    auto def_bag_profile = deflagaration_bag.get_profile_bag();
    def_bag_profile.write("comparisonData/profiles/ht/def_bag.csv");
    auto def_bag_spectrum = deflagaration_bag.get_spectrum_bag();
    def_bag_spectrum.write("comparisonData/spectra/ht/def_bag.csv");

    auto deflagaration_munu = DPInt::DeepPhaseResults(milestone, *eos_opt, 107.75, vw, dtauRs, DPInt::EoSModel::MUNU);
    auto def_munu_profile = deflagaration_munu.get_profile_munu();
    def_munu_profile.write("comparisonData/profiles/ht/def_munu.csv");
    auto def_munu_spectrum = deflagaration_munu.get_spectrum_munu();
    def_munu_spectrum.write("comparisonData/spectra/ht/def_munu.csv");

    auto deflagaration_veff = DPInt::DeepPhaseResults(milestone, *eos_opt, 107.75, vw, dtauRs, DPInt::EoSModel::VEFF);
    auto def_veff_profile = deflagaration_veff.get_profile_veff();
    def_veff_profile.write("comparisonData/profiles/ht/def_veff.csv");
    auto def_veff_spectrum = deflagaration_veff.get_spectrum_veff();
    def_veff_spectrum.write("comparisonData/spectra/ht/def_veff.csv");
  }
  
  if(vwDet > 0)
  {
    auto detonation_bag = DPInt::DeepPhaseResults(milestone, *eos_opt, 107.75, vwDet, dtauRs, DPInt::EoSModel::BAG);
    auto det_bag_profile = detonation_bag.get_profile_bag();
    det_bag_profile.write("comparisonData/profiles/ht/det_bag.csv");
    auto det_bag_spectrum = detonation_bag.get_spectrum_bag();
    det_bag_spectrum.write("comparisonData/spectra/ht/det_bag.csv");

    auto detonation_munu = DPInt::DeepPhaseResults(milestone, *eos_opt, 107.75, vwDet, dtauRs, DPInt::EoSModel::MUNU);
    auto det_munu_profile = detonation_munu.get_profile_munu();
    det_munu_profile.write("comparisonData/profiles/ht/det_munu.csv");
    auto det_munu_spectrum = detonation_munu.get_spectrum_munu();
    det_munu_spectrum.write("comparisonData/spectra/ht/det_munu.csv");

    auto detonation_veff = DPInt::DeepPhaseResults(milestone, *eos_opt, 107.75, vwDet, dtauRs, DPInt::EoSModel::VEFF);
    auto det_veff_profile = detonation_veff.get_profile_veff();
    det_veff_profile.write("comparisonData/profiles/ht/det_veff.csv");
    auto det_veff_spectrum = detonation_veff.get_spectrum_veff();
    det_veff_spectrum.write("comparisonData/spectra/ht/det_veff.csv");
  }

  return 0;
}