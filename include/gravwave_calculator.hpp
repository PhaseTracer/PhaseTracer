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

#ifndef PHASETRACER_GRAVWAVECALCULATOR_HPP_
#define PHASETRACER_GRAVWAVECALCULATOR_HPP_

#include <cmath>
#include <fstream>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "transition_finder.hpp"
#include "thermal_parameters.hpp"
#ifdef BUILD_WITH_DP
#include "deep_phase.hpp"
#endif

namespace PhaseTracer {

struct GravWaveSpectrum {
  double Tref;
  double alpha;
  double beta_H;
  double peak_frequency;
  double peak_amplitude;
  std::vector<double> frequency;
  std::vector<double> sound_wave;
  std::vector<double> turbulence;
  std::vector<double> bubble_collision;
  std::vector<double> total_amplitude;
  std::vector<double> SNR;

  std::vector<double> freq_ssm;
  std::vector<double> amplitude_ssm;
  
  /** Pretty-printer for GravWaveSpectrum */
  friend std::ostream &operator<<(std::ostream &o, const GravWaveSpectrum &a) {
    o << "=== gravitational wave spectrum generated at T = " << a.Tref << " ===" << "\n"
      << "alpha = " << a.alpha << "\n"
      << "beta over H = " << a.beta_H << "\n"
      << "peak frequency = " << a.peak_frequency << "\n"
      << "peak amplitude = " << a.peak_amplitude << "\n"
      << "signal to noise ratio for LISA = " << a.SNR[0] << std::endl;
    return o;
  }
};

/** Whether gc is initialused with TransitionFinder or ThermalParameters */
enum class GWCalcMode { FromTransition, FromThermalParameters };

class GravWaveCalculator {

public:
  explicit GravWaveCalculator(const TransitionFinder &tf_) : tf(std::make_unique<TransitionFinder>(tf_)), tp(nullptr), mode(GWCalcMode::FromTransition) {
    for (const auto &t : tf->get_transitions()) {
      if (std::isnan(t.TN))
        LOG(debug) << "Nucleation temperature dose not exist. GW will not be calculated !";
      else
        trans.push_back(t);
    }
  }

  explicit GravWaveCalculator(const ThermalParameters &tp_) : tf(nullptr), tp(std::make_unique<ThermalParameters>(tp_)), mode(GWCalcMode::FromThermalParameters) {
    for (const auto &tps : tp->get_thermal_params()) {
      if (std::isnan(tps.TP))
        LOG(debug) << "Percolation temperature dose not exist. GW will not be calculated!";
      else
        thermal_params.push_back(tps);
    }
  }

  /** Pretty-printer for set of transitions in this object */
  friend std::ostream &operator<<(std::ostream &o, const GravWaveCalculator &a);

  /** Functions to calculate amplitude at a fixed frequency */
  double Kappa_sound_wave(double alpha) const;
  double GW_bubble_collision(double f, double alpha, double beta_H, double T_ref) const;
  double GW_sound_wave(double f, double alpha, double beta_H, double T_ref) const;
  double GW_turbulence(double f, double alpha, double beta_H, double T_ref) const;

  /** Calculate GW spectrum for one transition */
  GravWaveSpectrum calc_spectrum(double alpha, double beta_H, double Tref);
  /** Calculate GW spectrums for all the transitions */
  std::vector<GravWaveSpectrum> calc_spectrums();
  /** Return GW spectrums for all the transitions */
  std::vector<GravWaveSpectrum> get_spectrums() const { return spectrums; }

  /** Sum GW spectrums */
  GravWaveSpectrum sum_spectrums(const std::vector<GravWaveSpectrum> &spectrums) const;
  /** Return the summed GW spectrum */
  GravWaveSpectrum get_total_spectrum() const { return total_spectrum; }

  /** Write a GW spectrum to a text file */
  void write_spectrum_to_text(const GravWaveSpectrum &sp, const std::string &filename) const;
  void write_spectrum_to_text(int i, const std::string &filename) const;
  void write_spectrum_to_text(const std::string &filename) const;

  /** Calcualte alpha(phase transition strength) with fixed phi */
  double get_alpha(const Eigen::VectorXd &vacuum_1, const Eigen::VectorXd &vacuum_2, double T) const;
  /** Calcualte alpha(phase transition strength) along the phases */
  double get_alpha(const Phase &phase1, const Phase &phase2, double T) const;

  /** Calcualte beta/H (Inverse phase transition duration) */
  double get_beta_H(const Phase &phase1, const Phase &phase2, double T, size_t i_unique) const;
  /** Sensitivity for LISA */
  double intergrand_SNR_LISA(double f, double alpha, double beta_H, double T_ref) const;
  /** Sensitivity for Taiji */
  double intergrand_SNR_Taiji(double f, double alpha, double beta_H, double T_ref) const;
  std::vector<double> get_SNR(double f_min, double f_max, double T_obs_LISA, double T_obs_Taiji, double alpha, double beta_H, double T_ref) const;

private:
  std::unique_ptr<TransitionFinder> tf;
  std::unique_ptr<ThermalParameters> tp;
  GWCalcMode mode;

  /** Degree of freedom */
  PROPERTY(double, dof, 106.75);
  /** Velocity of the bubble wall */
  PROPERTY(double, vw, 0.3);
  /** Ratio of efficiency factor of turbulence to the one of sound wave */
  PROPERTY(double, epsilon, 0.1);
  /** Gravitational constant */
  const double G = 6.7088e-39;

  /**Effecive run time**/
  PROPERTY(double, run_time_LISA, 3);
  PROPERTY(double, run_time_Taiji, 3);
  /**Acquisition time for LISA, Taiji, units second */
  PROPERTY(double, T_obs_LISA, 3 * 365.25 * 86400);
  PROPERTY(double, T_obs_Taiji, 3 * 365.25 * 86400);

  /**Integrate bound for SNR calculation */
  PROPERTY(double, SNR_f_min, 1e-5);
  PROPERTY(double, SNR_f_max, 1e-1);

  /** GW spectrums for all the transitions */
  std::vector<GravWaveSpectrum> spectrums;
  /** Summed GW spectrum*/
  GravWaveSpectrum total_spectrum;
  /** All transitions with valid TN*/
  std::vector<Transition> trans;
  /** All thermal params with valid temp*/
  std::vector<ThermalParams> thermal_params;

  /** Lower bound on the frequency of the GW spectrum */
  PROPERTY(double, min_frequency, 1e-4);
  /** Upper bound on the frequency of the GW spectrum  */
  PROPERTY(double, max_frequency, 1e1);
  /** Number of points for  the frequency of the GW spectrum */
  PROPERTY(int, num_frequency, 500);
  /** Number of points for  the frequency of the SSM GW spectrum */
  PROPERTY(int, num_frequency_ssm, 100);
  /** Temperature threshold for using bubble collision */
  PROPERTY(double, T_threshold_bubble_collision, 10);

  /** The step-size in numerical derivative for dVdT */
  PROPERTY(double, h_dVdT, 1e-2);
  /** The step-size in numerical derivative for dSdT */
  PROPERTY(double, h_dSdT, 1e-1);
  /** The number of points used in numerical derivative of action */
  PROPERTY(double, np_dSdT, 5);

  /**  */
  double rho_R(double T) const;

  /** Hubble constant */
  double H(double T) const;

  /** Functions for calculating alpha with fixed phi */
  double V(const Eigen::VectorXd &phi, const double T) const;
  double dVdT(const Eigen::VectorXd &phi, const double T) const;
  double rho(const Eigen::VectorXd &phi, const double T) const;

  /** Functions for calcualting alpha along the phase */
  Eigen::VectorXd vacuum_of_phase_at_T(const Phase &phase, double T) const;
  double V(const Phase &phase, const double T) const;
  double dVdT(const Phase &phase, const double T) const;
  double rho(const Phase &phase, const double T) const;

  /** Using linear regression for solving the derivative of the action functional */
  Eigen::Vector2d LinearRegression(std::vector<double> &x_, std::vector<double> &y_) const;
  /** Get action over temperature */
  double S3T(const Phase &phase1, const Phase &phase2, double T, size_t i_unique) const;
  /** Calculate derivative of S3/T to T */
  double dSdT(const Phase &phase1, const Phase &phase2, double T, size_t i_unique) const;
};

} // namespace PhaseTracer

#endif // PHASETRACER_GRAVWAVECALCULATOR_HPP_
