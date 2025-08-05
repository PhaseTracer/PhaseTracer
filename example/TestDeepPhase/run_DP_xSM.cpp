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
#include "phasetracer.hpp"

using json = nlohmann::json;

double L2diff(PhaseTracer::GravWaveSpectrum &sp) {

  if( sp.total_amplitude.size() == 0  || sp.amplitude_ssm.size() == 0) { 
    throw std::runtime_error("Error in GravWaveSpectrum amplitude data.");
  }

  size_t n = sp.total_amplitude.size();
  alglib::real_1d_array integrand_array, log_freq_array;
  integrand_array.setlength(n);
  log_freq_array.setlength(n);
  for ( size_t i = 0; i < n; i++) {
    double fg = sp.amplitude_ssm[i] != 0.0 ? sp.total_amplitude[i]/sp.amplitude_ssm[i] : 1.0;
    double lg = log(fg);
    integrand_array[i] = lg*lg * sp.frequency[i];
    log_freq_array[i] = log(sp.frequency[i]);
  }

  alglib::spline1dinterpolant spline;
  alglib::spline1dbuildcubic(log_freq_array, integrand_array, spline);

  double min_log_freq = log_freq_array[0];
  double max_log_freq = log_freq_array[n-1];
  
  if (min_log_freq > max_log_freq) {
    std::swap(min_log_freq, max_log_freq);
  }
  
  double norm = alglib::spline1dintegrate(spline, max_log_freq) - alglib::spline1dintegrate(spline, min_log_freq);
  LOG(debug) << "Integration bounds: [" << min_log_freq << ", " << max_log_freq << "], norm = " << norm << std::endl;
  return - norm;
}

void write_tps_debug(const PhaseTracer::ThermalParams &tp, const std::string &filename) {

  size_t i = tp.debug_info.temp.size();
  std::ofstream output_file(filename);
  if (!output_file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << " for writing" << std::endl;
    return;
  }
  
  output_file << "T,pf,nt,gamma,H,t" << std::endl;
  for (size_t j = 0; j < i; j++) {
    output_file << tp.debug_info.temp[j] << "," << tp.debug_info.pf[j]
                << "," << tp.debug_info.nt[j] << "," << tp.debug_info.gamma[j]
                << "," << tp.debug_info.hubble[j] << "," << tp.debug_info.t[j]
                << "\n";
  }
  output_file.close();
}

void write_eos(const PhaseTracer::ThermalParams &tp, const std::string &filename) {
  std::ofstream output_file(filename);
  if (!output_file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << " for writing" << std::endl;
    return;
  }

  output_file << "T,p,e,w,s\n";
  for (size_t i = 0; i < tp.eos.temp.size(); i++) {
    output_file << tp.eos.temp[i] << ","
                << tp.eos.pressure[i] << ","
                << tp.eos.energy[i] << ","
                << tp.eos.enthalpy[i] << ","
                << tp.eos.entropy[i] << "\n";
  }
  output_file.close();
}

/*
  Main
*/
int main(int argc, char* argv[]) {

  // default output
  std::string output_filename = "zData/vw_raw.csv";
  bool file_has_data = file_exists_and_not_empty(output_filename);

  std::ofstream output_file(output_filename, std::ios::app);
  if (!file_has_data) {
    output_file << "mS,lambda_s,lambda_hs,vw,Tref,alpha,beta,norm\n";
  }

  double Ms, lambda_s, lambda_hs, vw;
  bool running = true;
  int potential = 1;
  int matching = 1;

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
    Ms = modelParams["mS"].get<double>();
    lambda_s = modelParams["lambda_s"].get<double>();
    lambda_hs = modelParams["lambda_hs"].get<double>();
    vw = modelParams["vw"].get<double>();
  } catch (...) {
    std::cerr << "Input file might not be formatted correctly! Quitting now." << std::endl;
    // check find_phases suceeded
    // output_file << formatOutputNN({-1.,-1.,-1.}, -10);
    output_file << formatOutputNN({-1.,-1.,-1.});
    output_file.close();
    return 1;
  }

  std::vector<double> in = {Ms, lambda_s, lambda_hs, vw};
  EffectivePotential::DR_xsm model(Ms, lambda_hs, lambda_s);
  model.set_matching_flag(matching);
  model.set_potential_flag(potential);
  model.set_running_flag(running);

  PhaseTracer::PhaseFinder pf(model);
  pf.set_seed(0);
  pf.set_check_hessian_singular(true);
  pf.set_check_vacuum_at_high(false);
  pf.set_t_high(1000);
  pf.set_t_low(15);
  pf.set_dt_max_rel(1.);
  pf.set_dt_max_abs(2.5);

  try {
    pf.find_phases();
  } catch (const std::exception& e) {
    std::cerr << "Exception in pf.find_phases(): " << e.what() << std::endl;
    // output_file << formatOutput(in, -1);
    output_file << formatOutputNN(in);
    output_file.close();
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception in pf.find_phases()" << std::endl;
    // output_file << formatOutput(in, -1);
    output_file << formatOutputNN(in);
    output_file.close();
    return 1;
  }
  auto p1 = pf.get_phases();

  // check for phases
  if (p1.size() == 0) {
    // output_file << formatOutput(in, -2);
    output_file << formatOutputNN(in);
    output_file.close();
    return 1;
  }
  // for ( size_t i = 0; i < p1.size(); i++) {
  //   std::cout << "Phase " << i << "\n";
  //   for ( size_t j = 0; j < p1[i].V.size(); j++) {
  //     std::cout << p1[i].V[j] << " ";
  //     if( j % 5 == 4) {
  //       std::cout << "\n";
  //     }
  //   }
  // }
  std::cout << pf;

  PhaseTracer::ActionCalculator ac(model);
  PhaseTracer::TransitionFinder tf(pf, ac);
  tf.set_calculate_action(false); // This initialises with ac but stops it from calculating TN.
  try {
    tf.find_transitions();
  } catch (const std::exception& e) {
    std::cerr << "Exception in tf.find_transitions(): " << e.what() << std::endl;
    // output_file << formatOutput(in, -3);
    output_file << formatOutputNN(in);
    output_file.close();
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception in tf.find_transitions()" << std::endl;
    // output_file << formatOutput(in, -3);
    output_file << formatOutputNN(in);
    output_file.close();
    return 1;
  }
  std::cout << tf;

  auto t = tf.get_transitions();
  // check for transitions
  if (t.size() == 0) {
    // output_file << formatOutput(in, -4);
    output_file << formatOutputNN(in);
    output_file.close();
    return 1;
  }

  // testing
  PhaseTracer::ThermalParameters tp(tf);
  tp.set_compute_debug(true);
  tp.set_vw(vw);
  try {
    tp.find_thermal_parameters();
  } catch (const std::exception& e) {
    std::cerr << "Exception in tp.find_thermal_parameters(): " << e.what() << std::endl;
    // output_file << formatOutput(in, -5);
    output_file << formatOutputNN(in);
    output_file.close();
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception in tp.find_thermal_parameters()" << std::endl;
    // output_file << formatOutput(in, -5);
    output_file << formatOutputNN(in);
    output_file.close();
    return 1;
  }
  std::cout << tp;

  auto tps = tp.get_thermal_params();
  if (tps.size() == 0) {
    // output_file << formatOutput(in, -6);
    output_file << formatOutputNN(in);
    output_file.close();
    return 1;
  }

  // debugging information for tps.
  write_tps_debug(tps[0], "zData/tps_debug.csv");
  write_eos(tps[0], "zData/eos.csv");

  PhaseTracer::GravWaveCalculator gc(tp);
  gc.set_vw(vw);
  gc.set_min_frequency(1e-4);
  gc.set_max_frequency(1e+1);
  gc.set_num_frequency(1000);

  try {
    gc.calc_spectrums();
  } catch (const std::exception& e) {
    std::cerr << "Exception in gc.calc_spectrums(): " << e.what() << std::endl;
    // output_file << formatOutput(in, -7);
    output_file << formatOutputNN(in);
    output_file.close();
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception in gc.calc_spectrums()" << std::endl;
    // output_file << formatOutput(in, -7);
    output_file << formatOutputNN(in);
    output_file.close();
    return 1;
  }

  std::cout << gc;

  auto gw = gc.get_spectrums();

  double norm;
  try {
    norm = L2diff(gw[0]);
  } catch (const std::exception& e) {
    std::cerr << "Exception in L2diff(): " << e.what() << std::endl;
    // output_file << formatOutput(in, -8);
    output_file << formatOutputNN(in);
    output_file.close();
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception in L2diff()" << std::endl;
    // output_file << formatOutput(in, -8);
    output_file << formatOutputNN(in);
    output_file.close();
    return 1;
  }


  gc.write_spectrum_to_text(gw[0], "gw_spec.csv");

  std::vector<double> out = {gw[0].Tref, gw[0].alpha, gw[0].beta_H, norm};
  output_file << formatOutput(in, out);
  output_file.close();
  
  return 0;
}