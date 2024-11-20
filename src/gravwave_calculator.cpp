#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
#include <fstream>
#include <sstream>
#include <random>
#include "gravwave_calculator.hpp"
#include <boost/math/quadrature/trapezoidal.hpp>

namespace PhaseTracer {



std::ostream& operator << (std::ostream& o, const GravWaveCalculator& a){
  if (a.spectrums.empty()) {
    o << "found no spectrums" << std::endl;
    return o;
  }

  o << "found " << a.spectrums.size() << " spectrum";
  if (a.spectrums.size() > 1) {
    o << "s";
  }

  o << std::endl << std::endl;

  for (const auto &t : a.spectrums) {
    o << t << std::endl;
  }

  o << "=== total gravitational wave spectrum  ===" << std::endl;
  o << "peak frequency = " << a.total_spectrum.peak_frequency << std::endl
    << "peak amplitude = " << a.total_spectrum.peak_amplitude << std::endl;
  return o;
}

const double GravWaveCalculator::rho_R(double T){
  return std::pow(T, 4) * M_PI * M_PI * dof / 30.;
}

const double GravWaveCalculator::H(double T){
  return std::sqrt(8. * M_PI * G / 3. * rho_R(T));
}

const double GravWaveCalculator::V(const Eigen::VectorXd phi, const double T) {
  return tf.pf.P.V(phi,T);
}
const double GravWaveCalculator::dVdT(const Eigen::VectorXd phi, const double T) {
  return (-V(phi,T+2*h_dVdT) + 8*V(phi,T+h_dVdT) - 8*V(phi,T-h_dVdT) + V(phi,T-2*h_dVdT)) / (12.0 * h_dVdT);
}
const double GravWaveCalculator::rho(const Eigen::VectorXd phi, const double T){
  return V(phi,T) - 0.25*T*dVdT(phi,T);
}
const double GravWaveCalculator::get_alpha(const Eigen::VectorXd vacuum_1, const Eigen::VectorXd vacuum_2, double T){
  return ( rho(vacuum_1,T) - rho(vacuum_2,T) ) / rho_R(T);
}

Eigen::VectorXd GravWaveCalculator::vacuum_of_phase_at_T(Phase phase1, double T){
  return tf.pf.phase_at_T(phase1, T).x;
}
const double GravWaveCalculator::V(Phase phase1, const double T) {
  return tf.pf.P.V(vacuum_of_phase_at_T(phase1,T),T);
}
const double GravWaveCalculator::dVdT(Phase phase1, const double T) {
  return (-V(phase1,T+2*h_dVdT) + 8*V(phase1,T+h_dVdT) - 8*V(phase1,T-h_dVdT) + V(phase1,T-2*h_dVdT)) / (12.0 * h_dVdT);
}
const double GravWaveCalculator::rho(Phase phase1, const double T){
  return V(phase1,T) - 0.25*T*dVdT(phase1,T);
}
const double GravWaveCalculator::get_alpha(Phase phase1, Phase phase2, double T){
  return ( rho(phase1,T) - rho(phase1,T) ) / rho_R(T);
}

Eigen::Vector2d GravWaveCalculator::LinearRegression(std::vector<double> x_, std::vector<double> y_) {
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

const double GravWaveCalculator::S3T(Phase phase1, Phase phase2, double T, size_t i_unique){
  return tf.get_action(phase1, phase2, T, i_unique)/T;
}

const double GravWaveCalculator::dSdT(Phase phase1, Phase phase2, double T, size_t i_unique){
  std::vector<double> x,y;
  for (int ii=0; ii <= np_dSdT; ii++ ){
    double Ti = T - h_dSdT + ii * 2. * h_dSdT / np_dSdT;
    double S = S3T(phase1,phase2,Ti,i_unique);
    if ( (not std::isnan(S)) and (not std::isinf(S)) ) {
      x.push_back(Ti);
      y.push_back(S);
    }
  }
  if (x.size()<2){
    LOG(fatal) << "No enough valid S values in the calculation of dSdT.";
    exit(EXIT_FAILURE);
  }
  Eigen::Vector2d coeff = LinearRegression(x, y); // coeff[0] + coeff[1]*x
//      std::cout << "  T = np.array(" << x << ")" << std::endl;
//      std::cout << "  S = np.array(" << y << ")" << std::endl;
//      std::cout << "  coeff = np.array(" << coeff[0] << ", "<< coeff[1] << "])" << std::endl;
  return  coeff[1];
}

const double GravWaveCalculator::get_beta_H(Phase phase1, Phase phase2, double T, size_t i_unique){
  return T * dSdT(phase1, phase2, T, i_unique);
}

double GravWaveCalculator::GW_bubble_collision(double f, double alpha, double beta_H, double T_ref) {
  double omega_env;
  double s_env;
  double kappa = 1 / (1 + 0.715 * alpha) * (0.715 * alpha + 4. / 27 * sqrt(3 * alpha / 2));
  double f_peak_beta = 0.35 / (1 + 0.069 * vw + 0.69 * pow(vw, 4.));
  double f_env = 1.65e-5 * f_peak_beta * beta_H * (T_ref / 100) * pow(dof / 100, 1./ 6);
  double delta = 0.48 * pow(vw, 3.) / (1 + 5.3 * pow(vw, 2.) + 5 * pow(vw, 4.));
  s_env = pow(0.064 * pow(f/f_env, -3.) + (1 - 0.064 - 0.48) * pow(f / f_env, -1.) + 0.48 * (f / f_env),-1);
  omega_env = 1.67e-5 * delta * pow(beta_H, -2.) * pow(kappa * alpha / (1 + alpha), 2.) * pow(100 / dof, 1 / 3.) * s_env;
  return omega_env;
}

double GravWaveCalculator::GW_sound_wave(double f, double alpha, double beta_H, double T_ref) {
  double omega_sw;
  double omega_sw_peak;
  double Hstar_R = pow(beta_H, -1.) * pow(8 * 3.1415926, 1./ 3) * vw;
  double f_peak_sw = 2.6e-5 / Hstar_R * (T_ref / 100) * pow(dof / 100, 1./ 6);
  double K_sw = Kappa_sound_wave(alpha) * alpha / (1 + alpha);
  double H_tau = Hstar_R / sqrt(3 * K_sw /4);
  omega_sw_peak = 2.061 * 0.678 * 0.678 * 3.57e-5 * pow(100 / dof, 1./ 3) * 0.012 * Hstar_R * pow(K_sw, 2.);
  omega_sw = omega_sw_peak * pow(f / f_peak_sw, 3.) *
             pow(7. / (4. + 3. * pow(f / f_peak_sw, 2.)),7. / 2);
  omega_sw = omega_sw * std::min(1.0, H_tau);
  return omega_sw;
}

double GravWaveCalculator::GW_turbulence(double f, double alpha, double beta_H, double T_ref) {
  double omega_turb;
  double hn = 1.65e-5 * (T_ref / 100) * pow(dof/100, 1./ 6);
  double f_peak_turb = 2.7e-5 / vw * beta_H * (T_ref / 100) * pow(dof / 100, 1./ 6);
  double kappa_turb = Kappa_sound_wave(alpha) * epsilon;
  omega_turb = 3.35e-4 * pow(beta_H, -1.) * pow(kappa_turb * alpha / (1 + alpha), 3./2) * pow(100 / dof, 1./3)
            * vw * pow(f / f_peak_turb, 3) / (pow(1 + f / f_peak_turb, 11./3) * (1 + 8 * 3.1415926 * f / hn));
  return omega_turb;
}

double GravWaveCalculator::Kappa_sound_wave(double alpha) {
  double kappa_sw;
  double cs = sqrt(1 / 3.);
  double v_cj = 1 / (1 + alpha) * (cs + sqrt(pow(alpha,2.) + 2./3 * alpha));
  double kappa_a = pow(vw, 6./5) * 6.9 * alpha / (1.36 - 0.037 * sqrt(alpha) + alpha);
  double kappa_b = pow(alpha, 2./5) / (0.017 + pow(0.997 + alpha, 2./5));
  double kappa_c = sqrt(alpha) / (0.135 + sqrt(0.98 + alpha ));
  double kappa_d = alpha / (0.73 + 0.083 * sqrt(alpha) + alpha);
  double delta_kappa = -0.9 * log(sqrt(alpha) / (1 + sqrt(alpha)));

  if (0. < vw && vw <= cs){
    kappa_sw = pow(cs, 11./ 5) * kappa_a * kappa_b / ((pow(cs, 11./ 5) - pow(vw, 11./ 5)) * kappa_b + vw * pow(cs, 6./ 5) * kappa_a);
  }
  else if (cs < vw && vw < v_cj){
    kappa_sw = kappa_b + (vw - cs) * delta_kappa + pow(vw-cs, 3.) / pow(v_cj-cs, 3.) * (kappa_c - kappa_b - (v_cj -cs) * delta_kappa);
  }
  else if (v_cj <= vw && vw <= 1.){
    kappa_sw = pow(v_cj-1, 3.) * pow(v_cj, 5./ 2) * pow(vw, -5./ 2) * kappa_c * kappa_d / ((pow(v_cj-1, 3.) - pow(vw-1,3.)) * pow(v_cj, 5./ 2) * kappa_c + pow(vw-1,3.) * kappa_d);
  }
  else{
    LOG(fatal) << "Wrong bubble velocity.";
    exit(EXIT_FAILURE);
  }
  return kappa_sw;
}

GravWaveSpectrum GravWaveCalculator::calc_spectrum(double alpha, double beta_H, double Tref){
  if (num_frequency < 2){
    LOG(fatal) << "Number of frequencies must be greater than 1.";
    exit(EXIT_FAILURE);
  }
  if (max_frequency < min_frequency){
    LOG(fatal) << "Error: max_frequency < min_frequency";
    exit(EXIT_FAILURE);
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
    double total_amplitude = sound_wave+turbulence+bubble_collision;
    sp.total_amplitude.push_back(total_amplitude);
    if (total_amplitude > peak_amplitude){
        peak_amplitude = total_amplitude;
        peak_frequency = fq;
    }
  }
  sp.peak_frequency = peak_frequency;
  sp.peak_amplitude = peak_amplitude;
  sp.SNR = get_SNR(SNR_f_min, SNR_f_max, run_time_LISA, run_time_Taiji, alpha, beta_H, Tref);
  return sp;
}

GravWaveSpectrum GravWaveCalculator::sum_spectrums(std::vector<GravWaveSpectrum> sps){

  GravWaveSpectrum summed_sp;
  
  if (sps.empty()){
    LOG(fatal) << "No GW spectrums is generated. Can not sum spectrums";
    exit(EXIT_FAILURE);
  }
  
  double peak_frequency = 0;
  double peak_amplitude = 0;
  
  for (int ii=0; ii<sps[0].frequency.size(); ii++){
    summed_sp.frequency.push_back(sps[0].frequency[ii]);
    double sound_wave=0;
    double turbulence=0;
    double bubble_collision=0;
    double total_amplitude=0;
    for(int jj=0; jj<sps.size(); jj++){
      sound_wave += sps[jj].sound_wave[ii];
      turbulence += sps[jj].turbulence[ii];
      bubble_collision += sps[jj].bubble_collision[ii];
      sound_wave += sps[jj].sound_wave[ii];
      total_amplitude += sps[jj].total_amplitude[ii];
    }
    if (total_amplitude > peak_amplitude){
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

std::vector<GravWaveSpectrum> GravWaveCalculator::calc_spectrums(){
  for(const auto& ti : trans){
    double Tref = ti.TN;
    std::vector<Eigen::VectorXd>  vacua = tf.get_vacua_at_T(ti.true_phase,ti.false_phase,Tref,ti.key);
    double alpha = get_alpha(vacua[0],vacua[1],Tref);
    double beta_H = get_beta_H(ti.true_phase,ti.false_phase,Tref,ti.key);
    GravWaveSpectrum spi = calc_spectrum(alpha,beta_H,Tref);
    spectrums.push_back(spi);
  }
  total_spectrum = sum_spectrums(spectrums);
  return spectrums;
}

void GravWaveCalculator::write_spectrum_to_text(GravWaveSpectrum sp, const std::string &filename){
  std::ofstream file(filename);
  for (int ii=0; ii<sp.frequency.size(); ii++){
    file << sp.frequency[ii] << ", " << sp.total_amplitude[ii]<< ", "
         << sp.sound_wave[ii] << ", " << sp.turbulence[ii] << ", " << sp.bubble_collision[ii] << std::endl;
  }
  
  LOG(debug) << "GW spectrum has been written to " << filename;
}

void GravWaveCalculator::write_spectrum_to_text(int i, const std::string &filename){
  write_spectrum_to_text(spectrums[i], filename);
}

void GravWaveCalculator::write_spectrum_to_text(const std::string &filename){
  LOG(fatal) << spectrums.size();
  for (int ii=0; ii<spectrums.size(); ii++){
    write_spectrum_to_text(spectrums[ii], std::to_string(ii) + "_" + filename);
  }
}

double GravWaveCalculator::intergrand_SNR_LISA(double f, double alpha, double beta_H, double T_ref){
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

double GravWaveCalculator::intergrand_SNR_Taiji(double f, double alpha, double beta_H, double T_ref){
	double P_oms = 64e-24 * (1 + pow(2e-3 / f, 4)) * (2 * M_PI * f / 3e8);
	double P_acc = 9e-30 * (1 + pow(0.4e-3 / f, 2)) * (1 + pow(f/8e-3, 4)) * (1 / (2 * M_PI * f * 3e8));
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

std::vector<double> GravWaveCalculator::get_SNR(double f_min, double f_max, double run_time_LISA, double run_time_Taiji, double alpha, double beta_H, double T_ref){
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
}
