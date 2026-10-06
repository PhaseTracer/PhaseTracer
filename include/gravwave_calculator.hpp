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
#include <functional>
#include <limits>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "transition_finder.hpp"
#include "thermo_finder.hpp"

namespace PhaseTracer {

/**
 * @brief Enumeration of the backends used to compute a GW spectrum.
 *
 * FitFormulae uses the default fitting formulae originally shipped with PhaseTracer2.
 * SoundShell uses the full HydroGrav pipeline, including a calculation of the 
 * associated fluid profiles. Because HydroGrav uses the SSM, only the acoustic
 * contribution is included.
 */
enum class GravWaveMethod 
{
	FitFormulae,
	SoundShell
};

inline std::string to_string(GravWaveMethod m) 
{
  	return m == GravWaveMethod::SoundShell ? "sound shell model" : "fit formulae";
}

/**
 * @struct FluidProfile
 * @brief Self-similar fluid profile across the bubble wall.
 */
struct FluidProfile 
{
	/** @brief Self-similar coordinate xi = r/t. */
	std::vector<double> xi;

	/** @brief Fluid velocity v(xi). */
	std::vector<double> v;

	/** @brief Enthalpy density w(xi), normalised to its value at nucleation. */
	std::vector<double> w;

	/** @brief Lambda profile. */
	std::vector<double> lambda;

	/** @brief Temperature T(xi)/T_N. */
	std::vector<double> T;

	/** @brief Inner boundary of the integration domain: vw for deflagrations, c_- otherwise. */
	double xi_min = std::numeric_limits<double>::quiet_NaN();

	/** @brief Outer boundary of the integration domain, i.e. the shock position. */
	double xi_max = std::numeric_limits<double>::quiet_NaN();

	/** @brief Sound speed squared in the symmetric phase. */
	double cs_plus_sq = std::numeric_limits<double>::quiet_NaN();

	/** @brief Sound speed squared in the broken phase. */
	double cs_minus_sq = std::numeric_limits<double>::quiet_NaN();

	/** @brief Hydrodynamic mode: 0 = deflagration, 1 = hybrid, 2 = detonation, -1 = none. */
	int mode = -1;

	/** @brief Whether the shock front converged; if false a mu-nu fallback was used. */
	bool shock_converged = false;

	bool empty() const { return xi.empty(); }

	/** @brief Name of the hydrodynamic mode. */
	std::string mode_str() const;

	/** @brief Write the profile to a text file as xi, v, w, lambda, T */
	void write_profile_to_text(const std::string &filename) const;

	/** @brief Pretty-printer for FluidProfile */
	friend std::ostream &operator<<(std::ostream &o, const FluidProfile &a);
};

struct GravWaveSpectrum 
{
	/** @brief Thermal parameters used to calculate spectra. */
	double Tref = std::numeric_limits<double>::quiet_NaN();
	double Treh = std::numeric_limits<double>::quiet_NaN();
	double vw = std::numeric_limits<double>::quiet_NaN();
	double alpha = std::numeric_limits<double>::quiet_NaN();
	double alpha_fit = std::numeric_limits<double>::quiet_NaN();
	double g_eff = std::numeric_limits<double>::quiet_NaN();
	double h_eff = std::numeric_limits<double>::quiet_NaN();
	double cs = std::numeric_limits<double>::quiet_NaN();
	double beta_H = std::numeric_limits<double>::quiet_NaN();

	/** @brief Peak frequency of the gravitational wave spectrum. */
	double peak_frequency = std::numeric_limits<double>::quiet_NaN();

	/** @brief Peak amplitude of the gravitational wave spectrum. */
	double peak_amplitude = std::numeric_limits<double>::quiet_NaN();

	/** @brief Frequency grid for the spectrum. */
	std::vector<double> frequency;

	/** @brief Sound wave contribution to the spectrum. */
	std::vector<double> sound_wave;

	/** @brief Turbulence contribution to the spectrum. */
	std::vector<double> turbulence;

	/** @brief Bubble collision contribution to the spectrum. */
	std::vector<double> bubble_collision;

	/** @brief Total amplitude of the spectrum. */
	std::vector<double> total_amplitude;

	/** @brief LISA noise curve Omega h^2 on the frequency grid */
	std::vector<double> lisa_noise;

	/** @brief Taiji noise curve Omega h^2 on the frequency grid */
	std::vector<double> taiji_noise;

	/** @brief Signal-to-noise ratio for the spectrum. */
	std::vector<double> SNR;

	/** @brief Backend that produced the spectrum. */
	GravWaveMethod method = GravWaveMethod::FitFormulae;

	/** @brief Dimensionless momentum grid kRs. SoundShell only. */
	std::vector<double> kRs;

	/** @brief Fluid profile behind the spectrum. SoundShell only. */
	FluidProfile profile;

	/** @brief Sound wave lifetime. SoundShell only; NaN otherwise. */
	double dtau = std::numeric_limits<double>::quiet_NaN();

	/** @brief Pretty-printer for GravWaveSpectrum */
	friend std::ostream &operator<<(std::ostream &o, const GravWaveSpectrum &a) {
		o << "=== gravitational wave spectrum generated at T = " << a.Tref << " ===" << "\n"
		<< "method = " << to_string(a.method) << "\n"
		<< "vw = " << a.vw << "\n"
		<< "alpha = " << a.alpha << "\n";
		if (a.method == GravWaveMethod::SoundShell && !std::isnan(a.alpha_fit)) {
		o << "alpha for fit-formula contributions = " << a.alpha_fit << "\n";
		}
		o << "beta over H = " << a.beta_H << "\n"
		<< "peak frequency = " << a.peak_frequency << "\n"
		<< "peak amplitude = " << a.peak_amplitude << "\n";
		if (a.SNR.size() > 1) {
		o << "SNR for LISA = " << a.SNR[0] << "\n"
			<< "SNR for Taiji = " << a.SNR[1] << "\n";
		} else if (!a.SNR.empty()) {
		o << "SNR for LISA = " << a.SNR[0] << "\n";
		}
		if (!a.profile.empty()) {
		o << a.profile;
		}
		o << std::flush;
		return o;
	}
};

class GravWaveCalculator {

public:
	/* Backwards compatible version, will do nothing but throw an error. */
	explicit GravWaveCalculator(const TransitionFinder &tf_) : tf(tf_) {
		LOG(fatal) << "GravWaveCalculator constructed with TransitionFinder. "
		<< "This is no longer supported in PhaseTracer3. Please use the grav wave calculator class.";
		throw std::runtime_error("GravWaveCalculator constructed with TransitionFinder.");
	}

	explicit GravWaveCalculator(ThermoFinder &tm_) : tm(&tm_) {
		LOG(debug) << "GravWaveCalculator constructed with ThermoFinder";
	}

	/** Pretty-printer for set of transitions in this object */
	friend std::ostream &operator<<(std::ostream &o, const GravWaveCalculator &a);

	/** @brief Select the backend used to compute the spectrum. */
	void set_gw_method(GravWaveMethod m) 
	{
#ifndef BUILD_WITH_HG
		if (m == GravWaveMethod::SoundShell) 
		{
			LOG(fatal) << "Enable HydroGrav in CMake configuration before using it.";
			throw std::runtime_error("HydroGrav is not installed.");
		}
#endif
		gw_method = m;
	}

	/** 
	 * @brief Legacy functions to calculate amplitude at a fixed frequency.
	 * @param f The frequency at which to evaluate the gravitational wave contribution.
	 * @param alpha The phase transition strength parameter.
	 * @param beta_H The inverse phase transition duration parameter.
	 * @param T_ref The reference temperature for the phase transition.
	 * @return The gravitational wave contribution from the specified source at the given frequency.
	 * 
	 * @note These functions where included in PhaseTracer2, but have been replaced with
	 * the release of PhaseTracer3. They can be enabled using 'use_legacy_gw_methods'.
	 * They are kept for the purpose of backward replication of results.
	 */
	double GW_bubble_collision_legacy(double f, double alpha, double beta_H, double T_ref, double g_eff, double vw=0.3, double cs=std::sqrt(1/3)) const;
	double GW_sound_wave_legacy(double f, double alpha, double beta_H, double T_ref, double g_eff, double vw=0.3, double cs=std::sqrt(1/3)) const;
	double GW_turbulence_legacy(double f, double alpha, double beta_H, double T_ref, double g_eff, double vw=0.3, double cs=std::sqrt(1/3)) const;

	/** 
	 * @brief Functions to calculate the gravitational wave contributions from different sources.
	 * @param f The frequency at which to evaluate the gravitational wave contribution.
	 * @param alpha The phase transition strength parameter.
	 * @param beta_H The inverse phase transition duration parameter.
	 * @param T_ref The reference temperature for the phase transition.
	 * @return The gravitational wave contribution from the specified source at the given frequency.
	 */
	double GW_bubble_collision(double f, double alpha, double beta_H, double T_ref, double T_reh, double vw, double cs, double g_eff, double h_eff) const;
	double GW_sound_wave(double f, double alpha, double beta_H, double T_ref, double T_reh, double vw, double cs, double g_eff, double h_eff) const;
	double GW_turbulence(double f, double alpha, double beta_H, double T_ref, double T_reh, double vw, double cs, double g_eff, double h_eff) const;

	/** Calculate GW spectrum for one transition */
	GravWaveSpectrum calc_spectrum(const TransitionMilestone &milestone);

	/** Calculate GW spectrums for all the transitions */
	std::vector<GravWaveSpectrum> calc_spectrums();
	/** Return GW spectrums for all the transitions */
	std::vector<GravWaveSpectrum> get_spectrums() const { return spectrums; }

	/** Sum GW spectrums */
	GravWaveSpectrum sum_spectrums(const std::vector<GravWaveSpectrum> &spectrums) const;
	/** Return the summed GW spectrum */
	GravWaveSpectrum get_total_spectrum() const { return total_spectrum; }

	/**
	 * @brief Noise energy density Omega h^2 of the LISA sensitivity curve at frequency f.
	 * @param f Frequency at which to evaluate the LISA noise curve.
	 * @return Noise energy density Omega h^2 of the LISA sensitivity curve at frequency f.
	 */
	double noise_omega_LISA(double f) const;

	/**
	 * @brief Noise energy density Omega h^2 of the LISA sensitivity curve at frequency f.
	 * @param f Frequency at which to evaluate the legacy LISA noise curve.
	 * @return Noise energy density Omega h^2 of the legacy LISA sensitivity curve at frequency f.
	 * 
	 * @note This is the version used in PhaseTracer2.
	 */
	double noise_omega_LISA_legacy(double f) const;

	/**
	 * @brief Noise energy density Omega h^2 of the Taiji sensitivity curve at frequency f.
	 * @param f Frequency at which to evaluate the Taiji noise curve.
	 * @return Noise energy density Omega h^2 of the Taiji sensitivity curve at frequency f.
	 */
	double noise_omega_Taiji(double f) const;

	/** 
	 * @brief Total fit-formula amplitude at one frequency: sound wave + turbulence + collision
	 * @param f Frequency at which to evaluate the fit-formula amplitude.
	 * @param alpha Strength of the phase transition.
	 * @param beta_H Inverse duration of the phase transition normalized by the Hubble rate.
	 * @param T_ref Reference temperature of the phase transition.
	 * @param T_reh Reheating temperature after the phase transition.
	 * @param vw Bubble wall velocity.
	 * @param cs Sound speed in the plasma.
	 * @return Total fit-formula amplitude at the given frequency.
	 */
	double fit_omega(double f, double alpha, double beta_H, double T_ref, double T_reh, double vw, double cs, double g_eff, double h_eff) const;

	/**
	 * @brief Convenience overload of fit_omega that takes a TransitionMilestone.
	 * @param f Frequency at which to evaluate the fit-formula amplitude.
	 * @param m TransitionMilestone containing the phase transition parameters.
	 * @return Total fit-formula amplitude at the given frequency.
	 */
	double fit_omega(double f, TransitionMilestone m) const
	{
		return fit_omega(f, m.alpha, m.betaH_eff, m.temperature, m.reheating_temperature, m.vw, m.cs_plus, m.g_eff, m.h_eff);
	}

	/**
	 * @brief Signal-to-noise ratio (SNR) for a tabulated spectrum.
	 * @param frequency Vector of frequencies at which the spectrum is evaluated.
	 * @param omega Vector of energy densities corresponding to the frequencies.
	 * @return SNR for the tabulated spectrum as a vector {LISA, Taiji}.
	 */
	std::vector<double> get_SNR_tabulated(const std::vector<double> &frequency, const std::vector<double> &omega) const;

	/** 
	 * @brief Signal-to-noise ratio (SNR) for the fit-formula spectrum.
	 * @param alpha Strength of the phase transition.
	 * @param beta_H Inverse duration of the phase transition normalized by the Hubble rate.
	 * @param T_ref Reference temperature of the phase transition.
	 * @param T_reh Reheating temperature after the phase transition.
	 * @param vw Bubble wall velocity.
	 * @param cs Sound speed in the plasma.
	 * @return SNR for the fit-formula spectrum as a vector {LISA, Taiji}.
	 */
	std::vector<double> get_SNR(double alpha, double beta_H, double T_ref, double T_reh, double vw, double cs, double g_eff, double h_eff) const;

	/**
	 * @brief Convenience overload of get_SNR that takes a TransitionMilestone.
	 * @param m TransitionMilestone containing the phase transition parameters.
	 * @return SNR for the fit-formula spectrum as a vector {LISA, Taiji}.
	 */
	std::vector<double> get_SNR(TransitionMilestone m) const
	{
		return get_SNR(m.alpha, m.betaH_eff, m.temperature, m.reheating_temperature, m.vw, m.cs_plus, m.g_eff, m.h_eff);
	}

	/** 
	 * @brief Write a GW spectrum to a text file.
	 * @param sp The GW spectrum to write.
	 * @param filename The name of the text file to write to.
	 */
	void write_spectrum_to_text(const GravWaveSpectrum &sp, const std::string &filename) const;

	/**
	 * @brief Write the GW spectrum of the i-th transition to a text file.
	 * @param i Index of the transition.
	 * @param filename The name of the text file to write to.
	 */
	void write_spectrum_to_text(int i, const std::string &filename) const;

	/**
	 * @brief Write the total GW spectrum to a text file.
	 * @param filename The name of the text file to write to.
	 */
	void write_spectrum_to_text(const std::string &filename) const;

private:

	std::optional<TransitionFinder> tf;
	/** Non-owning; the ThermoFinder must outlive this calculator. */
	ThermoFinder *tm = nullptr;

	/** GW spectrums for all the transitions */
	std::vector<GravWaveSpectrum> spectrums;

	/** Summed GW spectrum*/
	GravWaveSpectrum total_spectrum;

	/** All transitions with valid TN */
	std::vector<Transition> trans;

	/** The milestone of a thermal parameter set selected by default_milestone.*/
	const TransitionMilestone *milestone_of(const ThermalParameterSet &tps) const;

	/** Turn the noise-weighted integrals into {LISA, Taiji} SNRs using the run times */
	std::vector<double> SNR_from_integrals(double snr_sq_LISA, double snr_sq_Taiji) const;

#ifdef BUILD_WITH_HG
	/** Calculate a GW spectrum for one transition with HydroGrav's sound shell model */
	GravWaveSpectrum calc_spectrum_ssm(const ThermalParameterSet &tps, const TransitionMilestone &milestone) const;
#endif

	/**
	 * Evaluate the turbulence and bubble collision fits on a spectrum's own
	 * frequency grid and store them in its turbulence and bubble_collision.
	 */
	void add_fit_contributions(GravWaveSpectrum &sp, double alpha_fit) const;

	/**
	 * Sum the three contributions into total_amplitude, then set the peak and SNR
	 * from it. Call after every contribution is in place.
	 */
	void finalise_spectrum(GravWaveSpectrum &sp) const;

	/** Fill lisa_noise and taiji_noise on the spectrum's frequency grid */
	void add_noise_curves(GravWaveSpectrum &sp) const;

	/**
	 * Omega h^2 of a Robson, Cornish & Liu (arXiv:1803.01944) sky-averaged sensitivity curve with
	 * arm length L [m], transfer frequency f_star [Hz], optical-metrology noise P_oms [m/sqrt(Hz)] and
	 * acceleration noise P_acc [m s^-2/sqrt(Hz)], plus the 4-yr galactic confusion-noise fit.
	 */
	double noise_omega_RCL(double f, double L, double f_star, double P_oms, double P_acc) const;

	/**
	 * @brief Helper to calculate RH from betaH.
	 * @param betaH The inverse time scale of the phase transition.
	 * @return The Hubble radius corresponding to the given betaH.
	 */
    double get_RH(const double& betaH) const;

	std::array<double, 1> get_collision_peaks(const double &H0_star, const double &RH) const;
	std::array<double, 2> get_sound_wave_peaks(const double& H0_star, const double& RH, const double& vw=std::sqrt(1/3), const double& cs=std::sqrt(1/3)) const;
	std::array<double, 3> get_turbulence_peaks(const double& H0_star, const double& RH, const double& K_turb) const;

	/**
	 * @brief Helper to calculate K from the given parameters.
	 * @param A A constant factor.
	 * @param kappa The efficiency factor.
	 * @param alpha The strength of the phase transition.
	 * @return The value of K based on the given parameters.
	 */
	double get_K(const double& A, const double& kappa, const double& alpha) const;

	/**
	 * @brief Helper to calculate the efficiency factor for sound waves.
	 * @param alpha The strength of the phase transition.
	 * @param vw The bubble wall velocity.
	 * @param cs The speed of sound in the plasma (default is 1/sqrt(3)).
	 * @return The efficiency factor for sound waves based on the given parameters.
	 */
	double get_kappa_sw(const double& alpha, const double& vw, const double& cs=std::sqrt(1./3.)) const;

	/**
	 * @brief Helper to calculate the efficiency factor for turbulence.
	 * @param alpha The strength of the phase transition.
	 * @param vw The bubble wall velocity.
	 * @param cs The speed of sound in the plasma (default is 1/sqrt(3)).
	 * @return The efficiency factor for turbulence based on the given parameters.
	 */
	double get_kappa_turb(const double& alpha, const double& vw, const double& cs=std::sqrt(1./3.)) const;

	/**
	 * @brief Helper to calculate the efficiency factor for bubble collisions.
	 * @param alpha The strength of the phase transition.
	 * @return The efficiency factor for bubble collisions based on the given parameters.
	 */
	double get_kappa_col(const double& alpha) const;

	/**
	 * @brief Helper to calculate the prefactor used in the gravitational wave spectrum calculations.
	 * @return The calculated prefactor based on the current cosmological parameters.
	 */
	double get_prefactor(const double& g_eff, const double& h_eff) const;

	/**
	 * @brief Helper to calculate the Hubble rate today based on the reference temperature.
	 * @param Tref The reference temperature.
	 * @return The Hubble rate today corresponding to the given reference temperature.
	 */
	// Can this be obtained by redshifting Hstar obtained from FE?
	double get_Hubble_rate_today(const double& T, const double& g_eff, const double& h_eff) const;

	/**
	 * @brief Helper to calculate the normalization factor for the sound wave contribution to the gravitational wave spectrum.
	 * @param peak_freqs An array containing the peak frequencies of the sound wave contribution.
	 * @return The calculated normalization factor for the sound wave contribution based on the given peak frequencies.
	 */
	double get_sound_wave_N(const std::array<double, 2>& peak_freqs) const;

	/**
	 * @brief Helper to calculate the Y parameter for sound waves.
	 * @param RH The Hubble radius at the time of the phase transition.
	 * @param K The efficiency factor for sound waves.
	 * @return The calculated Y parameter for sound waves based on the given parameters.
	 */
	double get_sound_wave_Y(const double& RH, const double& K) const;

	/** @brief Singly broken power law template function.
	 * @param f The frequency at which to evaluate the power law.
	 * @param f0 The first break frequency.
	 * @param f1 The second break frequency.
	 * @param n0 The power law index before the first break.
	 * @param n1 The power law index between the first and second breaks.
	 * @param a1 The parameter controlling the smoothness of the first break.
	 * @param b1 Constant appearing in first break.
	 * @return The value of the singly broken power law at the given frequency.
	 */
	double singly_broken_power_law(
		const double& f, 
		const double& f0, 
		const double& f1, 
		const double& n0, 
		const double& n1, 
		const double& a1, 
		const double& b1) const;

	/**
	 * @brief Doubly broken power law template function.
	 * @param f The frequency at which to evaluate the power law.
	 * @param f0 The first break frequency.
	 * @param f1 The second break frequency.
	 * @param f2 The third break frequency.
	 * @param n0 The power law index before the first break.
	 * @param n1 The power law index between the first and second breaks.
	 * @param n2 The power law index after the second break.
	 * @param a1 The parameter controlling the smoothness of the first break.
	 * @param a2 The parameter controlling the smoothness of the second break.
	 * @return The value of the doubly broken power law at the given frequency.
	 */
	double doubly_broken_power_law(
		const double& f, const double& f0, const double& f1, const double& f2,
		const double& n0, const double& n1, const double& n2, 
		const double& a1, const double& a2) const;

	/** Amplitude constants for different gravitational wave sources. */
	constexpr static double A_col = 0.0026;
	constexpr static double A_sw = 0.11;
	constexpr static double A_turb = 0.255;

	/** Whether to use the old PhaseTracer2 fitting formulas */
	PROPERTY(bool, use_legacy_gw_methods, false);

	/** Lower bound on the kRs values of the SSM spectrum */
	PROPERTY(double, min_kRs_value, 1e-3);
	/** Upper bound on the kRs values of the SSM spectrum */
	PROPERTY(double, max_kRs_value, 1e3);
	/** Number of points in the kRs grid */
	PROPERTY(int, n_kRs_value, 200);

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

	/** ThermoFinder Milestone */
	PROPERTY(PhaseTracer::MilestoneType, default_milestone, PhaseTracer::MilestoneType::PERCOLATION);

	/** Backend used by calc_spectrums; set through set_gw_method */
	PROPERTY_CUSTOM_SETTER(GravWaveMethod, gw_method, GravWaveMethod::FitFormulae);

	/** Include collisions and turbulence in SSM spectrum */
	PROPERTY(bool, include_col_and_turb_in_ssm, false);

	
	/** Relativistic degrees of freedom at present */
	PROPERTY(double, g_0, 2.0);
	/** Entropy degrees of freedom at present*/
	PROPERTY(double, h_0, 3.91);
	
	/**
	 * @brief Effective degrees of freedom.
	 * 
	 * The effective relativistic degrees of freedom. These overwrite
	 * the value contained in the transition milestone. If set to zero, the value
	 * in the milestone is used.
	 */ 
	PROPERTY(double, g_eff, 0.0);
	PROPERTY(double, h_eff, 0.0);
	
	/** Neutrino energy density factor */
	PROPERTY(double, omega_hsq_neutrino, 2.473e-5);

	/** Entropy injection factor */
	PROPERTY(double, D, 1.0);

	/** 
	 * @brief Bubble wall velocity.
	 * 
	 * The bubble wall velocity used to calculate the GW spectrum. This overwrites
	 * the value contained in the transition milestone. If set to zero, the value
	 * in the milestone is used.
	 */
	PROPERTY(double, vw, 0.0);

	/** Ratio of efficiency factor of turbulence to the one of sound wave */
	PROPERTY(double, epsilon, 0.1);

	/**Effective observation time in years for LISA, Taiji**/
	PROPERTY(double, run_time_LISA, 4);
	PROPERTY(double, run_time_Taiji, 3);
	/** Use the old LISA noise curve (noise_omega_LISA_legacy) instead of Robson, Cornish & Liu */
	PROPERTY(bool, use_legacy_LISA_noise, false);

	/**Integrate bound for SNR calculation */
	PROPERTY(double, SNR_f_min, 1e-5);
	PROPERTY(double, SNR_f_max, 1e-1);
	/** Number of ln f grid points per decade used by get_SNR to tabulate the fit-formula spectrum */
	PROPERTY(double, SNR_steps_per_decade, 200);

};

} // namespace PhaseTracer

#endif // PHASETRACER_GRAVWAVECALCULATOR_HPP_
