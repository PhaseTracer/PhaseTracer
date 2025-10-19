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

/*
  Main
*/
int main(int argc, char* argv[]) {

  // default output
  std::string output_filename = "scanData/output.csv";
  bool file_has_data = file_exists_and_not_empty(output_filename);

  std::ofstream output_file(output_filename, std::ios::app);
  if (!file_has_data) {
    output_file << "mS,lambda_s,lambda_hs,Tref,alpha,betaH,dtau,we,vw,snr_fit,snr_ssm\n";
  }

  double Ms, lambda_s, lambda_hs, vw, dtauRs;
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
    vw = modelParams["vw"].get<double>();
    dtauRs = modelParams["dtauRs"].get<double>();
  } catch (...) {
    std::cerr << "Input file might not be formatted correctly! Quitting now." << std::endl;
    // check find_phases suceeded
    output_file << formatOutput({-1.,-1.,-1.}, -10);
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
    output_file << formatOutput(in, -1);
    output_file.close();
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception in pf.find_phases()" << std::endl;
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

  // std::cout << pf;

  PhaseTracer::ActionCalculator ac(model);
  PhaseTracer::TransitionFinder tf(pf, ac);
  tf.set_calculate_action(false); // This initialises with ac but stops it from calculating TN.
  try {
    tf.find_transitions();
  } catch (const std::exception& e) {
    std::cerr << "Exception in tf.find_transitions(): " << e.what() << std::endl;
    output_file << formatOutput(in, -3);
    output_file.close();
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception in tf.find_transitions()" << std::endl;
    output_file << formatOutput(in, -3);
    output_file.close();
    return 1;
  }
  // std::cout << tf;

  auto t = tf.get_transitions();
  // check for transitions
  if (t.size() == 0) {
    output_file << formatOutput(in, -4);
    output_file.close();
    return 1;
  }

  // testing
  PhaseTracer::ThermalParameters tp(tf);
  tp.set_compute_debug(false);
  tp.set_vw(vw);
  tp.set_dof(107.75);
  try {
    tp.find_thermal_parameters();
  } catch (const std::exception& e) {
    std::cerr << "Exception in tp.find_thermal_parameters(): " << e.what() << std::endl;
    output_file << formatOutput(in, -5);
    output_file.close();
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception in tp.find_thermal_parameters()" << std::endl;
    output_file << formatOutput(in, -5);
    output_file.close();
    return 1;
  }
  // std::cout << tp;

  auto tps = tp.get_thermal_params();
  if (tps.size() == 0) {
    output_file << formatOutput(in, -6);
    output_file.close();
    return 1;
  }

  // optional, use wallgo to calculate vw and replace
  std::cout << "original vw = " << vw << "\n";
  double tperc = tps[0].TP;
  const auto params4d_map = model.get_4d_parameter_map(tperc);
  vw = wallGoWrapper::getWallVelocity(t[0], params4d_map, tperc);
  if ( vw == 1.0 ) { vw -= 1e-6; }
  std::cout << "new vw = " << vw << "\n";

  /*
    Here we manually calculate the sound spectrum using both PhaseTracer and DeepPhase.
    This can more easily by done by calling gc.find_spectrums()
  */
  PhaseTracer::GravWaveCalculator gc(tp);
  PhaseTracer::GW_DeepPhase gcDeepPhase;
  gc.set_vw(vw);
  gc.set_dof(107.75);
  gc.set_min_frequency(1e-4);
  gc.set_max_frequency(1e+1);
  gc.set_num_frequency(200);
  gc.set_num_frequency_ssm(200);

  // calculate gravitational wave spectrum using DeepPhase
  // note this is normally wrapped inside GravWaveCalculator
  std::cout << "Calculating spectra with vw = " << gc.get_vw() << "\n";
  double RsH = pow(8.*M_PI, 1./3.) * vw / tps[0].betaH_tp;
  double dtau = RsH * dtauRs;
  auto ssm_result = gcDeepPhase.get_ssm_amplitude(tps[0], gc.get_vw(), gc.get_dof(), dtau, gc.get_min_frequency(), gc.get_max_frequency(), gc.get_num_frequency());
  std::vector<double> ssm_freq = ssm_result.first;
  std::vector<double> ssm_amp = ssm_result.second;
  auto ssm_peaks = gcDeepPhase.obtain_peaks(ssm_freq, ssm_amp);
  double ssm_peak_freq = ssm_peaks.first;
  double ssm_peak_amp = ssm_peaks.second;
  double zp = gcDeepPhase.zp;
  double Gamma = tps[0].we_tp;

  std::vector<double> fit_freq, fit_amp;
  double logMin = std::log10(gc.get_min_frequency());
  double logMax = std::log10(gc.get_max_frequency());
  double logInterval = (logMax - logMin) / (gc.get_num_frequency()- 1);
  double fit_peak_amp = 0;
  double fit_peak_freq = 0;
  for (int i = 0; i < gc.get_num_frequency(); ++i) 
  {
    double fq = std::pow(10, logMin + i * logInterval);
    fit_freq.push_back(fq);
    double sound_wave = gc.GW_sound_wave(fq, tps[0].alpha_tp, tps[0].betaH_tp, tps[0].TP, 10.0, 4./3.);
    fit_amp.push_back(sound_wave);
    if (sound_wave > fit_peak_amp) {
      fit_peak_amp = sound_wave;
      fit_peak_freq = fq;
    }
  }

  std::cout << "========================\n";
  std::cout << "Tref = " << tps[0].TP << "\n";
  std::cout << "alpha = " << tps[0].alpha_tp << "\n";
  std::cout << "betaH = " << tps[0].betaH_tp << "\n";
  std::cout << "H = " << tps[0].H_tp << "\n";
  std::cout << "w/e = " << tps[0].we_tp << "\n";
  std::cout << "zp = " << zp << "\n";
  std::cout << "vw = " << vw << "\n";
  std::cout << "dtau = " << dtau << "\n";
  std::cout << "cs_true (tp) = " << tps[0].cs_true_tp << "\n";
  std::cout << "cs_false (tp) = " << tps[0].cs_false_tp << "\n";
  std::cout << "Peak amplitude (fit) = " << fit_peak_amp << "\n";
  std::cout << "Peak amplitude (ssm) = " << ssm_peak_amp << "\n";
  std::cout << "Peak frequency (fit) = " << fit_peak_freq << "\n";
  std::cout << "Peak frequency (ssm) = " << ssm_peak_freq << "\n";

  // write_spectra(fit_amp, fit_freq, ssm_amp, ssm_freq, "scanData/spec/spec_" + std::to_string(lambda_hs) + "_" + std::to_string(lambda_s) + ".csv");
  // tps[0].eos.write_EoS("scanData/eos/eos_" + std::to_string(lambda_hs) + "_" + std::to_string(lambda_s) + ".csv");

  write_spectra(fit_amp, fit_freq, ssm_amp, ssm_freq, "spec_" + std::to_string(lambda_hs) + "_" + std::to_string(lambda_s) + ".csv");
  tps[0].eos.write_EoS("eos_" + std::to_string(lambda_hs) + "_" + std::to_string(lambda_s) + ".csv");

  double snr_fit = gwTools::getSNR(fit_freq, fit_amp);
  double snr_ssm = gwTools::getSNR(ssm_freq, ssm_amp);

  std::vector<double> out = {tps[0].TP, tps[0].alpha_tp, tps[0].betaH_tp, tps[0].dtf_tp, tps[0].we_tp, vw, log10(snr_fit), log10(snr_ssm)};
  output_file << formatOutput(in, out);
  output_file.close();
  return 0;
}