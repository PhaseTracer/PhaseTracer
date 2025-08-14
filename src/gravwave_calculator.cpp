// ====================================================================
// This file is part of PhaseTracer

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ====================================================================

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <tuple>
#include <vector>

#include <boost/math/quadrature/trapezoidal.hpp>

#include "gravwave_calculator.hpp"

namespace PhaseTracer {

std::ostream &operator<<(std::ostream &o, const GravWaveCalculator &a) {
  if (a.spectrums.empty()) {
    o << "found no spectrums" << std::endl;
    return o;
  }

  o << "found " << a.spectrums.size() << " spectrum";
  if (a.spectrums.size() > 1) {
    o << "s";
  }

  o << std::endl
    << std::endl;

  for (const auto &t : a.spectrums) {
    o << t << std::endl;
  }

  o << "=== total gravitational wave spectrum ===" << std::endl
    << "peak frequency = " << a.total_spectrum.peak_frequency << std::endl
    << "peak amplitude = " << a.total_spectrum.peak_amplitude << std::endl;
  return o;
}

double GravWaveCalculator::rho_R(double T) const {
  return std::pow(T, 4) * M_PI * M_PI * dof / 30.;
}

double GravWaveCalculator::H(double T) const {
  return std::sqrt(8. * M_PI * G / 3. * rho_R(T));
}

double GravWaveCalculator::V(const Eigen::VectorXd &phi, const double T) const {
  return tf->pf.P.V(phi, T);
}
double GravWaveCalculator::dVdT(const Eigen::VectorXd &phi, const double T) const {
  return (-V(phi, T + 2 * h_dVdT) + 8 * V(phi, T + h_dVdT) - 8 * V(phi, T - h_dVdT) + V(phi, T - 2 * h_dVdT)) / (12.0 * h_dVdT);
}
double GravWaveCalculator::rho(const Eigen::VectorXd &phi, const double T) const {
  return V(phi, T) - 0.25 * T * dVdT(phi, T);
}
double GravWaveCalculator::get_alpha(const Eigen::VectorXd &vacuum_1, const Eigen::VectorXd &vacuum_2, double T) const {
  return (rho(vacuum_1, T) - rho(vacuum_2, T)) / rho_R(T);
}
double GravWaveCalculator::V(const Phase &phase, const double T) const {
  return tf->pf.phase_at_T(phase, T).potential;
}
double GravWaveCalculator::dVdT(const Phase &phase, const double T) const {
  return (-V(phase, T + 2 * h_dVdT) + 8 * V(phase, T + h_dVdT) - 8 * V(phase, T - h_dVdT) + V(phase, T - 2 * h_dVdT)) / (12.0 * h_dVdT);
}
double GravWaveCalculator::rho(const Phase &phase, const double T) const {
  return V(phase, T) - 0.25 * T * dVdT(phase, T);
}
double GravWaveCalculator::get_alpha(const Phase &phase1, const Phase &phase2, double T) const {
  return (rho(phase1, T) - rho(phase2, T)) / rho_R(T);
}

Eigen::Vector2d GravWaveCalculator::LinearRegression(std::vector<double> &x_, std::vector<double> &y_) const {
  Eigen::Map<Eigen::VectorXd> x(x_.data(), x_.size());
  Eigen::Map<Eigen::VectorXd> y(y_.data(), y_.size());
  int n = x.size();
  Eigen::MatrixXd A(n, 2);
  A.col(0) = Eigen::VectorXd::Ones(n);
  A.col(1) = x;
  // Calculate (A^T * A)^(-1) * A^T * y
  Eigen::VectorXd coeff = (A.transpose() * A).ldlt().solve(A.transpose() * y);
  return coeff;
}

double GravWaveCalculator::S3T(const Phase &phase1, const Phase &phase2, double T, size_t i_unique) const {
  return tf->get_action(phase1, phase2, T, i_unique) / T;
}

double GravWaveCalculator::dSdT(const Phase &phase1, const Phase &phase2, double T, size_t i_unique) const {
  std::vector<double> x, y;
  for (int ii = 0; ii <= np_dSdT; ii++) {
    double Ti = T - h_dSdT + ii * 2. * h_dSdT / np_dSdT;
    double S = S3T(phase1, phase2, Ti, i_unique);
    if ((!std::isnan(S)) && (!std::isinf(S))) {
      x.push_back(Ti);
      y.push_back(S);
    }
  }
  if (x.size() < 2) {
    throw std::runtime_error("Not enough valid action values in the calculation of dS/dT");
  }
  Eigen::Vector2d coeff = LinearRegression(x, y);
  return coeff[1];
}

double GravWaveCalculator::get_beta_H(const Phase &phase1, const Phase &phase2, double T, size_t i_unique) const {
  return T * dSdT(phase1, phase2, T, i_unique);
}

double GravWaveCalculator::GW_bubble_collision(double f, double alpha, double beta_H, double T_ref) const {
  double omega_env;
  double s_env;
  double kappa = 1 / (1 + 0.715 * alpha) * (0.715 * alpha + 4. / 27 * sqrt(3 * alpha / 2));
  double f_peak_beta = 0.35 / (1 + 0.069 * vw + 0.69 * pow(vw, 4.));
  double f_env = 1.65e-5 * f_peak_beta * beta_H * (T_ref / 100) * pow(dof / 100, 1. / 6);
  double delta = 0.48 * pow(vw, 3.) / (1 + 5.3 * pow(vw, 2.) + 5 * pow(vw, 4.));
  s_env = pow(0.064 * pow(f / f_env, -3.) + (1 - 0.064 - 0.48) * pow(f / f_env, -1.) + 0.48 * (f / f_env), -1);
  omega_env = 1.67e-5 * delta * pow(beta_H, -2.) * pow(kappa * alpha / (1 + alpha), 2.) * pow(100 / dof, 1 / 3.) * s_env;
  return omega_env;
}

double GravWaveCalculator::GW_sound_wave(double f, double alpha, double beta_H, double T_ref) const {
  double F_gw0 = 3.57e-5 * pow(100/dof, 1./3.);
  double zp = 10;
  double Gamma = 4./3.;
  double HRs = pow(8.*M_PI, 1./3.) * vw / beta_H;
  double f_peak_sw = 2.6e-5 * (zp/10) * (T_ref/100.) * pow(dof/100., 1./6.) / HRs;
  double S_sw = (f/f_peak_sw) * (f/f_peak_sw) * (f/f_peak_sw) * pow(7./(4. + 3.*(f/f_peak_sw)*(f/f_peak_sw)), 7./2.);
  double K_sw = Kappa_sound_wave(alpha) * alpha / (1 + alpha);
  double omega_sw_peak = 2.061 * 0.678*0.678 * F_gw0 * K_sw*K_sw * HRs * 0.012;
  double omega_sw = omega_sw_peak * S_sw;
  double H_tau = HRs / sqrt(K_sw / Gamma);
  omega_sw = omega_sw * std::min(1.0, H_tau);
  return omega_sw;
}

double GravWaveCalculator::GW_turbulence(double f, double alpha, double beta_H, double T_ref) const {
  double omega_turb;
  double hn = 1.65e-5 * (T_ref / 100) * pow(dof / 100, 1. / 6);
  double f_peak_turb = 2.7e-5 / vw * beta_H * (T_ref / 100) * pow(dof / 100, 1. / 6);
  double kappa_turb = Kappa_sound_wave(alpha) * epsilon;
  omega_turb = 3.35e-4 * pow(beta_H, -1.) * pow(kappa_turb * alpha / (1 + alpha), 3. / 2) * pow(100 / dof, 1. / 3) * vw * pow(f / f_peak_turb, 3) / (pow(1 + f / f_peak_turb, 11. / 3) * (1 + 8 * 3.1415926 * f / hn));
  return omega_turb;
}

double GravWaveCalculator::Kappa_sound_wave(double alpha) const {
  double cs = sqrt(1 / 3.);
  double v_cj = 1 / (1 + alpha) * (cs + sqrt(pow(alpha, 2.) + 2. / 3 * alpha));
  double kappa_a = pow(vw, 6. / 5) * 6.9 * alpha / (1.36 - 0.037 * sqrt(alpha) + alpha);
  double kappa_b = pow(alpha, 2. / 5) / (0.017 + pow(0.997 + alpha, 2. / 5));
  double kappa_c = sqrt(alpha) / (0.135 + sqrt(0.98 + alpha));
  double kappa_d = alpha / (0.73 + 0.083 * sqrt(alpha) + alpha);
  double delta_kappa = -0.9 * log(sqrt(alpha) / (1 + sqrt(alpha)));

  if (0. < vw && vw <= cs) {
    return pow(cs, 11. / 5) * kappa_a * kappa_b / ((pow(cs, 11. / 5) - pow(vw, 11. / 5)) * kappa_b + vw * pow(cs, 6. / 5) * kappa_a);
  } else if (cs < vw && vw < v_cj) {
    return kappa_b + (vw - cs) * delta_kappa + pow(vw - cs, 3.) / pow(v_cj - cs, 3.) * (kappa_c - kappa_b - (v_cj - cs) * delta_kappa);
  } else if (v_cj <= vw && vw <= 1.) {
    return pow(v_cj - 1, 3.) * pow(v_cj, 5. / 2) * pow(vw, -5. / 2) * kappa_c * kappa_d / ((pow(v_cj - 1, 3.) - pow(vw - 1, 3.)) * pow(v_cj, 5. / 2) * kappa_c + pow(vw - 1, 3.) * kappa_d);
  }
  throw std::runtime_error("Invalid bubble wall velocity (vw > 1)");
}

GravWaveSpectrum GravWaveCalculator::calc_spectrum(double alpha, double beta_H, double Tref) {
  if (num_frequency < 2) {
    throw std::runtime_error("Number of frequencies must be greater than 1");
  }
  if (max_frequency < min_frequency) {
    throw std::runtime_error("max_frequency < min_frequency");
  }

  GravWaveSpectrum sp;
  sp.Tref = Tref;
  sp.alpha = alpha;
  sp.beta_H = beta_H;
  double logMin = std::log10(max_frequency);
  double logMax = std::log10(min_frequency);
  double logInterval = (logMax - logMin) / (num_frequency - 1);
  double peak_frequency = 0;
  double peak_amplitude = 0;
  for (int i = 0; i < num_frequency; ++i) {
    double fq = std::pow(10, logMin + i * logInterval);
    sp.frequency.push_back(fq);
    double sound_wave = GW_sound_wave(fq, alpha, beta_H, Tref);
    sp.sound_wave.push_back(sound_wave);
    double turbulence = GW_turbulence(fq, alpha, beta_H, Tref);
    sp.turbulence.push_back(turbulence);
    double bubble_collision = Tref < T_threshold_bubble_collision ? GW_bubble_collision(fq, alpha, beta_H, Tref) : 0;
    sp.bubble_collision.push_back(bubble_collision);
    double total_amplitude = sound_wave + turbulence + bubble_collision;
    sp.total_amplitude.push_back(total_amplitude);
    if (total_amplitude > peak_amplitude) {
      peak_amplitude = total_amplitude;
      peak_frequency = fq;
    }
  }
  sp.peak_frequency = peak_frequency;
  sp.peak_amplitude = peak_amplitude;
  sp.SNR = get_SNR(SNR_f_min, SNR_f_max, run_time_LISA, run_time_Taiji, alpha, beta_H, Tref);
  return sp;
}

GravWaveSpectrum GravWaveCalculator::sum_spectrums(const std::vector<GravWaveSpectrum> &sps) const {

  GravWaveSpectrum summed_sp;

  if (sps.empty()) {
    throw std::runtime_error("No GW spectrums were given - cannot sum them");
  }

  double peak_frequency = 0;
  double peak_amplitude = 0;

  for (int ii = 0; ii < sps[0].frequency.size(); ii++) {
    summed_sp.frequency.push_back(sps[0].frequency[ii]);
    double sound_wave = 0;
    double turbulence = 0;
    double bubble_collision = 0;
    double total_amplitude = 0;
    for (int jj = 0; jj < sps.size(); jj++) {
      sound_wave += sps[jj].sound_wave[ii];
      turbulence += sps[jj].turbulence[ii];
      bubble_collision += sps[jj].bubble_collision[ii];
      sound_wave += sps[jj].sound_wave[ii];
      total_amplitude += sps[jj].total_amplitude[ii];
    }
    if (total_amplitude > peak_amplitude) {
      peak_amplitude = total_amplitude;
      peak_frequency = sps[0].frequency[ii];
    }
    summed_sp.sound_wave.push_back(sound_wave);
    summed_sp.turbulence.push_back(turbulence);
    summed_sp.bubble_collision.push_back(bubble_collision);
    summed_sp.total_amplitude.push_back(total_amplitude);
  }

  return summed_sp;
}

std::vector<GravWaveSpectrum> GravWaveCalculator::calc_spectrums() {
  if(mode == GWCalcMode::FromTransition) {
    if (!tf) { throw std::runtime_error("TransitionFinder is not set for GravWaveCalculator");}
    for (const auto &ti : trans) {
      double Tref = ti.TN;
      if (Tref < 1.5 * h_dSdT) {
        continue;
      }
      std::vector<Eigen::VectorXd> vacua = tf->get_vacua_at_T(ti.true_phase, ti.false_phase, Tref, ti.key);
      double alpha = get_alpha(vacua[0], vacua[1], Tref);
      double beta_H = get_beta_H(ti.true_phase, ti.false_phase, Tref, ti.key);
      GravWaveSpectrum spi = calc_spectrum(alpha, beta_H, Tref);
      spectrums.push_back(spi);
    }
  } else if (mode == GWCalcMode::FromThermalParameters) {
    if(!tp) {throw std::runtime_error("ThermalParams is not set for GravWaveCalculator");}
    #ifdef BUILD_WITH_DP
    GW_DeepPhase dp;
    #endif
    for (const auto &tps : thermal_params) {
      double Tref, alpha, beta_H;
      if (tps.nucleates == MilestoneStatus::YES) {
        Tref = tps.TN;
        alpha = tps.alpha_tn;
        beta_H = tps.betaH_tn;
      } else if (tps.percolates == MilestoneStatus::YES) {
        Tref = tps.TP;
        alpha = tps.alpha_tp;
        beta_H = tps.betaH_tp;
      } // no handling of failure in both
      GravWaveSpectrum spi = calc_spectrum(alpha, beta_H, Tref);
      #ifdef BUILD_WITH_DP
      std::pair<std::vector<double>, std::vector<double>> ssm_spectrum = dp.get_ssm_amplitude(tps, vw, dof, min_frequency, max_frequency, num_frequency_ssm);
      spi.freq_ssm = ssm_spectrum.first;
      spi.amplitude_ssm = ssm_spectrum.second;
      auto maxes = dp.get_ssm_max(spi.amplitude_ssm, spi.freq_ssm);
      LOG(debug) << "SSM peak amplitude = " << maxes.first << ", frequency = " << maxes.second;
      if (maxes.first > spi.peak_amplitude) {
        spi.peak_amplitude = maxes.first;
        spi.peak_frequency = maxes.second;
      }
      #endif
      spectrums.push_back(spi);
    }
  }
  total_spectrum = sum_spectrums(spectrums);
  return spectrums;
}

void GravWaveCalculator::write_spectrum_to_text(const GravWaveSpectrum &sp, const std::string &filename) const {
  std::ofstream file(filename);
  for (int ii = 0; ii < sp.frequency.size(); ii++) {
    file << sp.frequency[ii] << ", " << sp.total_amplitude[ii] << ", " << sp.sound_wave[ii] << ", " << sp.turbulence[ii] << ", " << sp.bubble_collision[ii];
    #ifdef BUILD_WITH_DP
    file << ", " <<  sp.freq_ssm[ii] << ", " << sp.amplitude_ssm[ii];
    #endif
    file << std::endl;
  }

  LOG(debug) << "GW spectrum has been written to " << filename;
}

void GravWaveCalculator::write_spectrum_to_text(int i, const std::string &filename) const {
  write_spectrum_to_text(spectrums[i], filename);
}

void GravWaveCalculator::write_spectrum_to_text(const std::string &filename) const {
  LOG(info) << "writing " << spectrums.size() << "spectrums to text";
  for (int ii = 0; ii < spectrums.size(); ii++) {
    write_spectrum_to_text(spectrums[ii], std::to_string(ii) + "_" + filename);
  }
}

double GravWaveCalculator::intergrand_SNR_LISA(double f, double alpha, double beta_H, double T_ref) const {
  double P_oms = 3.6e-41;
  double P_acc = 1.44e-48 / pow(2 * M_PI * f, 4) * (1 + pow(0.4e-3 / f, 2));
  double S_A = sqrt(2) * 20. / 3 * (P_oms + 4 * P_acc) * (1 + pow(f / (2.54e-2), 2));
  double H_0 = 67.4 / (3.086e19);
  double omegahsq_lisa = 4 * M_PI * M_PI / (3 * H_0 * H_0) * pow(f, 3) * S_A * 0.674 * 0.674;
  double sound_wave = GW_sound_wave(f, alpha, beta_H, T_ref);
  double turbulence = GW_turbulence(f, alpha, beta_H, T_ref);
  double bubble_collision = T_ref < T_threshold_bubble_collision ? GW_bubble_collision(f, alpha, beta_H, T_ref) : 0;
  double omegahsq = sound_wave + turbulence + bubble_collision;
  return omegahsq * omegahsq / (omegahsq_lisa * omegahsq_lisa);
}

double GravWaveCalculator::intergrand_SNR_Taiji(double f, double alpha, double beta_H, double T_ref) const {
  double P_oms = 64e-24 * (1 + pow(2e-3 / f, 4)) * (2 * M_PI * f / 3e8);
  double P_acc = 9e-30 * (1 + pow(0.4e-3 / f, 2)) * (1 + pow(f / 8e-3, 4)) * (1 / (2 * M_PI * f * 3e8));
  double f_star = 3e8 / (2 * M_PI * 3e9);
  double N_A = 8 * std::sin(f / f_star) * std::sin(f / f_star) * (4 * P_acc * (1 + std::cos(f / f_star) + pow(std::cos(f / f_star), 2)) + P_oms * (2 + std::cos(f / f_star)));
  double R_A = 9. / 20 * 1. / (1 + pow(3 * f / (4 * f_star), 2)) * (2 - 2 * std::cos(2 * f / f_star));
  /*inverse-noise weighted sensitivity, 1/sqrt(2) accounts for the S_E contribution */
  double S_A = N_A / R_A * 1 / sqrt(2);
  double H_0 = 67.4 / (3.086e19);
  double omegahsq_taiji = 4 * M_PI * M_PI / (3 * H_0 * H_0) * pow(f, 3) * S_A * 0.674 * 0.674;
  double sound_wave = GW_sound_wave(f, alpha, beta_H, T_ref);
  double turbulence = GW_turbulence(f, alpha, beta_H, T_ref);
  double bubble_collision = T_ref < T_threshold_bubble_collision ? GW_bubble_collision(f, alpha, beta_H, T_ref) : 0;
  double omegahsq = sound_wave + turbulence + bubble_collision;
  return omegahsq * omegahsq / (omegahsq_taiji * omegahsq_taiji);
}

std::vector<double> GravWaveCalculator::get_SNR(double f_min, double f_max, double run_time_LISA, double run_time_Taiji, double alpha, double beta_H, double T_ref) const {
  std::vector<double> SNR_vector;
  double T_obs_LISA = run_time_LISA * 365.25 * 86400;
  double T_obs_Taiji = run_time_Taiji * 365.25 * 86400;
  auto intergrand_SNR_LISA_fix_var = std::bind(&GravWaveCalculator::intergrand_SNR_LISA, this, std::placeholders::_1, alpha, beta_H, T_ref);
  auto intergrand_SNR_Taiji_fix_var = std::bind(&GravWaveCalculator::intergrand_SNR_Taiji, this, std::placeholders::_1, alpha, beta_H, T_ref);
  double SNR_LISA_sq = boost::math::quadrature::trapezoidal(intergrand_SNR_LISA_fix_var, f_min, f_max);
  double SNR_Taiji_sq = boost::math::quadrature::trapezoidal(intergrand_SNR_Taiji_fix_var, f_min, f_max);
  SNR_vector.push_back(sqrt(SNR_LISA_sq * T_obs_LISA));
  SNR_vector.push_back(sqrt(SNR_Taiji_sq * T_obs_Taiji));
  return SNR_vector;
}
} // namespace PhaseTracer
