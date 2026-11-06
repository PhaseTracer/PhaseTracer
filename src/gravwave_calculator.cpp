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

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <tuple>
#include <vector>

#include "gravwave_calculator.hpp"
#include "scale.hpp"

namespace PhaseTracer {

namespace {

/* The redshift fits below are calibrated in GeV; temperatures arrive in units of PhaseTracer::scale */
double to_GeV(double T) { return T / scale(); }

/* Composite Simpson's rule for y(x) on an ascending, possibly non-uniform grid. An odd number
	of intervals is closed with the end-interval correction used by scipy.integrate.simpson */
double simpson(const std::vector<double> &x, const std::vector<double> &y)
{
	if (x.size() < 2) return 0.;
	const size_t n = x.size() - 1;
	if (n == 1) return 0.5 * (x[1] - x[0]) * (y[0] + y[1]);

	double result = 0.;
	const size_t n_even = n - n % 2;
	for (size_t ii = 0; ii + 2 <= n_even; ii += 2)
	{
		const double h0 = x[ii + 1] - x[ii];
		const double h1 = x[ii + 2] - x[ii + 1];
		const double hsum = h0 + h1;
		result += hsum / 6. * (y[ii] * (2. - h1 / h0) + y[ii + 1] * hsum * hsum / (h0 * h1) + y[ii + 2] * (2. - h0 / h1));
	}
	if (n % 2 == 1)
	{
		const double h0 = x[n - 1] - x[n - 2];
		const double h1 = x[n] - x[n - 1];
		const double alpha = (2. * h1 * h1 + 3. * h0 * h1) / (6. * (h0 + h1));
		const double beta = (h1 * h1 + 3. * h0 * h1) / (6. * h0);
		const double eta = h1 * h1 * h1 / (6. * h0 * (h0 + h1));
		result += alpha * y[n] + beta * y[n - 1] - eta * y[n - 2];
	}
	return result;
}

} // namespace

// GravWaveSpectrum Definitions

const SNRResult &
GravWaveSpectrum::snr_for(const std::string &name) const
{
	for (const auto &r : snr) 
	{
		if (r.detector == name) return r;
	}
	throw std::out_of_range("No SNR result for detector " + name);
}

// Fluid Profile Definitions

std::string 
FluidProfile::mode_str() const 
{
    switch (mode) {
    case 0: 	return "deflagration";
    case 1: 	return "hybrid";
    case 2: 	return "detonation";
    default: 	return "none";
    }
}

void 
FluidProfile::write_profile_to_text(const std::string &filename) const 
{
	std::ofstream file(filename);
	file << "xi,v,w,lambda,T" << std::endl;
	for (int ii = 0; ii < xi.size(); ii++) 
	{
		file << xi[ii] << "," << v[ii] << "," << w[ii] << "," << lambda[ii] << "," << T[ii];
		file << std::endl;
	}

  LOG(debug) << "Fluid profile has been written to " << filename;
}

std::ostream &operator<<(std::ostream &o, const FluidProfile &a) 
{
	if (a.empty()) 
	{
		o << "no fluid profile" << std::endl;
		return o;
	}

	o << "hydrodynamic mode = " << a.mode_str() << "\n"
		<< "xi range = [" << a.xi_min << ", " << a.xi_max << "]" << "\n"
		<< "shock converged = " << (a.shock_converged ? "yes" : "no") << std::endl;
  return o;
}

std::ostream &operator<<(std::ostream &o, const GravWaveCalculator &a) 
{
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

// GravWaveCalculator Definitions

double 
GravWaveCalculator::GW_bubble_collision_legacy(double f, double alpha, double beta_H, double T_ref, double g_eff, double vw, double cs) const 
{
	double omega_env;
	double s_env;
	double kappa = 1 / (1 + 0.715 * alpha) * (0.715 * alpha + 4. / 27 * sqrt(3 * alpha / 2));
	double f_peak_beta = 0.35 / (1 + 0.069 * vw + 0.69 * pow(vw, 4.));
	double f_env = 1.65e-5 * f_peak_beta * beta_H * (to_GeV(T_ref) / 100) * pow(g_eff / 100, 1. / 6);
	double delta = 0.48 * pow(vw, 3.) / (1 + 5.3 * pow(vw, 2.) + 5 * pow(vw, 4.));
	s_env = pow(0.064 * pow(f / f_env, -3.) + (1 - 0.064 - 0.48) * pow(f / f_env, -1.) + 0.48 * (f / f_env), -1);
	omega_env = 1.67e-5 * delta * pow(beta_H, -2.) * pow(kappa * alpha / (1 + alpha), 2.) * pow(100 / g_eff, 1 / 3.) * s_env;
	return omega_env;
}

double 
GravWaveCalculator::GW_sound_wave_legacy(double f, double alpha, double beta_H, double T_ref, double g_eff, double vw, double cs) const 
{
	double zp = 10.0;
	double Gamma = 4./3.;
	double F_gw0 = 3.57e-5 * pow(100/g_eff, 1./3.);
	double HRs = pow(8.*M_PI, 1./3.) * vw / beta_H;
	double f_peak_sw = 2.6e-5 * (zp/10) * (to_GeV(T_ref)/100.) * pow(g_eff/100., 1./6.) / HRs;
	double S_sw = (f/f_peak_sw) * (f/f_peak_sw) * (f/f_peak_sw) * pow(7./(4. + 3.*(f/f_peak_sw)*(f/f_peak_sw)), 7./2.);
	double K_sw = get_kappa_sw(alpha, vw, cs) * alpha / (1 + alpha);
	double omega_sw_peak = 2.061 * 0.678*0.678 * F_gw0 * K_sw*K_sw * HRs * 0.012;
	double omega_sw = omega_sw_peak * S_sw;
	double H_tau = HRs / sqrt(K_sw / Gamma);
	omega_sw = omega_sw * std::min(1.0, H_tau);
	return omega_sw;
}

double 
GravWaveCalculator::GW_turbulence_legacy(double f, double alpha, double beta_H, double T_ref, double g_eff, double vw, double cs) const 
{
	double legacy_vw = (vw==0) ? 0.3 : vw;
	double omega_turb;
	double hn = 1.65e-5 * (to_GeV(T_ref) / 100) * pow(g_eff / 100, 1. / 6);
	double f_peak_turb = 2.7e-5 / vw * beta_H * (to_GeV(T_ref) / 100) * pow(g_eff / 100, 1. / 6);
	double kappa_turb = get_kappa_turb(alpha, vw, cs);
	omega_turb = 3.35e-4 * pow(beta_H, -1.) * pow(kappa_turb * alpha / (1 + alpha), 3. / 2) * pow(100 / g_eff, 1. / 3) * vw * pow(f / f_peak_turb, 3) / (pow(1 + f / f_peak_turb, 11. / 3) * (1 + 8 * 3.1415926 * f / hn));
	return omega_turb;
}

double 
GravWaveCalculator::GW_bubble_collision(double f, double alpha, double beta_H, double T_ref, double T_reh, double vw, double cs, double g_eff, double h_eff) const
{
	const double RH = get_RH(beta_H);
	const double prefactor = get_prefactor(g_eff, h_eff);
	const double H0_star = get_Hubble_rate_today(T_reh, g_eff, h_eff);
	const double kappa_col = get_kappa_col(alpha);
	const double K_col = get_K(1., kappa_col, alpha);

	auto [f_col] = get_collision_peaks(H0_star, RH);

	const double S_col = singly_broken_power_law(f, f_col, f_col, 2.4, -4, 1.2, 0.5);

	return prefactor * A_col * K_col*K_col * RH*RH * S_col;
}

double 
GravWaveCalculator::GW_sound_wave(double f, double alpha, double beta_H, double T_ref, double T_reh, double vw, double cs, double g_eff, double h_eff) const
{
	const double RH = get_RH(beta_H);
	const double prefactor = get_prefactor(g_eff, h_eff);
	const double H0_star = get_Hubble_rate_today(T_reh, g_eff, h_eff);
	const double kappa_sw = get_kappa_sw(alpha, vw, cs);
	const double K_sw = get_K(0.6, kappa_sw, alpha);

	const auto sound_wave_peaks = get_sound_wave_peaks(H0_star, RH, vw, cs);
	const auto [f_sw_1, f_sw_2] = sound_wave_peaks;

	const auto N = get_sound_wave_N(sound_wave_peaks);
	const auto Y_sw = get_sound_wave_Y(RH, K_sw);

	const double S_sw = doubly_broken_power_law(f, f_sw_1, f_sw_1, f_sw_2, 3, -1, -1, 2, 4);

	return prefactor * A_sw * K_sw*K_sw * N * Y_sw * RH * S_sw;
}

double 
GravWaveCalculator::GW_turbulence(double f, double alpha, double beta_H, double T_ref, double T_reh, double vw, double cs, double g_eff, double h_eff) const
{
	const double RH = get_RH(beta_H);
	const double prefactor = get_prefactor(g_eff, h_eff);
	const double H0_star = get_Hubble_rate_today(T_reh, g_eff, h_eff);
	const double kappa_turb = get_kappa_turb(alpha, vw, cs);
	const double K_turb = get_K(0.6, kappa_turb, alpha);

	const auto turbulence_peaks = get_turbulence_peaks(H0_star, RH, K_turb);
	const auto [f_turb_1, f_turb_2, f_turb_3] = turbulence_peaks;

	auto log_factor = [&H0_star](double freq)
	{
	double temp = std::log(1 + H0_star/(2.*M_PI*freq));
	return temp*temp;
	};
	double mod = f <= f_turb_3 ? log_factor(f_turb_3) : log_factor(f);

	const double S_turb = singly_broken_power_law(f, f_turb_1, f_turb_2, 3, -7.9, 2.15, 1.);

	return prefactor * A_turb * K_turb*K_turb * RH*RH*RH * S_turb * mod;
}

GravWaveSpectrum 
GravWaveCalculator::calc_spectrum(const TransitionMilestone &milestone) 
{
	if (num_frequency < 2) 
	{
		throw std::runtime_error("Number of frequencies must be greater than 1");
	}
	if (max_frequency < min_frequency) 
	{
		throw std::runtime_error("max_frequency < min_frequency");
	}

	auto use_custom_check = [this](double input, std::string name) {
		bool use_custom = input > 0.0;
		if(use_custom) {LOG(debug) << "Creating GW spectrum with user defined " + name + " = " << input;}
		return use_custom;
	};

	bool use_custom_vw = use_custom_check(vw, "vw");
	bool use_custom_g_eff = use_custom_check(g_eff, "g_eff");
	bool use_custom_h_eff = use_custom_check(h_eff, "g_eff");

	GravWaveSpectrum sp;
	sp.Tref = milestone.temperature;
	sp.Treh = milestone.reheating_temperature;
	sp.alpha = milestone.alpha;
	sp.alpha_fit = milestone.alpha;
	sp.beta_H = milestone.betaH_eff;
	sp.cs = milestone.cs_plus;
	sp.vw = use_custom_vw ? vw : milestone.vw;
	sp.g_eff = use_custom_g_eff ? g_eff : milestone.g_eff;
	sp.h_eff = use_custom_h_eff ? h_eff : milestone.h_eff;
	LOG(debug) << "Calculating GW spectrum for milestone: temperature = " << milestone.temperature
			   << ", reheating_temperature = " << milestone.reheating_temperature
	           << ", alpha = " << milestone.alpha
	           << ", betaH_eff = " << milestone.betaH_eff
	           << ", cs_plus = " << milestone.cs_plus
	           << ", vw = " << (use_custom_vw ? vw : milestone.vw)
			   << ", g_eff = " << (use_custom_g_eff ? g_eff : milestone.g_eff)
			   << ", h_eff = " << (use_custom_h_eff ? h_eff : milestone.h_eff);

	// these can all be precomputed, but only if not using legacy
	double RH, prefactor, H0_star;
	double kappa_col, kappa_sw, kappa_turb;
	double K_col, K_sw, K_turb;
	double N_sw, Y_sw;
	double f_col;
	std::array<double, 2> sound_wave_peaks;
	std::array<double, 3> turb_peaks;

	auto log_factor = [&H0_star](double freq)
	{
		double temp = std::log(1 + H0_star/(2.*M_PI*freq));
		return temp*temp;
	};

	if(!use_legacy_gw_methods)
	{
		RH = get_RH(sp.beta_H);
		prefactor = get_prefactor(sp.g_eff, sp.h_eff);
		H0_star = get_Hubble_rate_today(sp.Treh, sp.g_eff, sp.h_eff);

		kappa_col = get_kappa_col(sp.alpha);
		K_col = get_K(1., kappa_col, sp.alpha);
		f_col = get_collision_peaks(H0_star, RH)[0];

		kappa_sw = get_kappa_sw(sp.alpha, sp.vw, sp.cs); 
		K_sw = get_K(0.6, kappa_sw, sp.alpha);
		sound_wave_peaks = get_sound_wave_peaks(H0_star, RH, sp.vw, sp.cs);

		N_sw = get_sound_wave_N(sound_wave_peaks);
		Y_sw = get_sound_wave_Y(RH, K_sw);

		kappa_turb = get_kappa_turb(sp.alpha, sp.vw, sp.cs); 
		K_turb = get_K(0.6, kappa_turb, sp.alpha);
		turb_peaks = get_turbulence_peaks(H0_star, RH, K_turb);
	}

	double logMin = std::log10(max_frequency);
	double logMax = std::log10(min_frequency);
	double logInterval = (logMax - logMin) / (num_frequency - 1);
	double peak_frequency = 0;
	double peak_amplitude = 0;

	for (int i = 0; i < num_frequency; ++i) 
	{
		double fq = std::pow(10, logMin + i * logInterval);
		sp.frequency.push_back(fq);

		double sound_wave, turbulence, bubble_collision;
		bool include_col = to_GeV(sp.Tref) < T_threshold_bubble_collision;
		if(use_legacy_gw_methods) 
		{
			sound_wave = GW_sound_wave_legacy(fq, sp.alpha, sp.beta_H, sp.Tref, sp.g_eff, sp.vw, sp.cs);
			turbulence = GW_turbulence_legacy(fq, sp.alpha, sp.beta_H, sp.Tref, sp.g_eff, sp.vw, sp.cs);
			bubble_collision = include_col ? GW_bubble_collision_legacy(fq, sp.alpha, sp.beta_H, sp.Tref, sp.g_eff, sp.vw, sp.cs) : 0;
		} else {
			const auto [f_sw_1, f_sw_2] = sound_wave_peaks;
			const double S_sw = doubly_broken_power_law_LISA(fq, f_sw_1, f_sw_2, 3, 1, -3, 2, 4);
			sound_wave = prefactor * A_sw * K_sw*K_sw * N_sw * Y_sw * RH * S_sw;

			const auto [f_turb_1, f_turb_2, f_turb_3] = turb_peaks;
			const double S_turb = singly_broken_power_law(fq, f_turb_1, f_turb_2, 3, -7.9, 2.15, 1.);
			const double mod = fq <= f_turb_3 ? log_factor(f_turb_3) : log_factor(fq);
			turbulence = prefactor * A_turb * K_turb*K_turb * RH*RH*RH * S_turb * mod;

			const double S_col = include_col ? singly_broken_power_law(fq, f_col, f_col, 2.4, -4, 1.2, 0.5) : 0;
			bubble_collision = prefactor * A_col * K_col*K_col * RH*RH * S_col;
		}
		
		sp.sound_wave.push_back(sound_wave);
		sp.turbulence.push_back(turbulence);
		sp.bubble_collision.push_back(bubble_collision);
		
		double total_amplitude = sound_wave + turbulence + bubble_collision;
		sp.total_amplitude.push_back(total_amplitude);
		if (total_amplitude > peak_amplitude) 
		{
			peak_amplitude = total_amplitude;
			peak_frequency = fq;
		}
	}
	sp.peak_frequency = peak_frequency;
	sp.peak_amplitude = peak_amplitude;
	/* Evaluate the fit formulae over each detector's whole band, which the frequency grid need not cover */
	sp.snr = get_SNR(sp.alpha, sp.beta_H, sp.Tref, sp.Treh, sp.vw, sp.cs, sp.g_eff, sp.h_eff);
	add_noise_curves(sp);
	LOG(debug) << "GW spectrum: peak f = " << peak_frequency << " Hz, peak h^2 Omega = " << peak_amplitude;
	for (const auto &r : sp.snr) 
	{
		LOG(debug) << r;
	}
	return sp;
}

std::vector<GravWaveSpectrum> GravWaveCalculator::calc_spectrums() 
{
	if (tf)
	{
		LOG(fatal) << "calc_spectrums called on an instance of GravWaveCalculator constructed using TransitionFinder.";
		throw std::runtime_error("GravWaveCalculator constructed with TransitionFinder.");
	}

	if (tm)
	{
		for (const auto &tps : tm->get_thermal_parameters()) 
		{
			const TransitionMilestone *milestone = milestone_of(tps);
			if (milestone == nullptr) 
			{
				continue;
			}
			if (gw_method == GravWaveMethod::SoundShell) {
#ifdef BUILD_WITH_HG
        		spectrums.push_back(calc_spectrum_ssm(tps, *milestone));
#else
        		throw std::runtime_error("HydroGrav is not installed.");
#endif
			} else {
				spectrums.push_back(calc_spectrum(*milestone));
			}
		}
		total_spectrum = sum_spectrums(spectrums);
		return spectrums;
	} else {
		throw std::runtime_error("No TransitionFinder or ThermoFinder provided to GravWaveCalculator");
	}
}

GravWaveSpectrum
GravWaveCalculator::sum_spectrums(const std::vector<GravWaveSpectrum> &sps) const 
{
	GravWaveSpectrum summed_sp;

	if (sps.empty()) 
	{
		throw std::runtime_error("No GW spectrums were given - cannot sum them");
	}

	for (const auto &sp : sps) 
	{
		if (sp.frequency != sps[0].frequency) 
		{
			throw std::runtime_error("GW spectrums are on different frequency grids - cannot sum them");
		}
	}

	summed_sp.method = sps[0].method;

	double peak_frequency = 0;
	double peak_amplitude = 0;

	for (int ii = 0; ii < sps[0].frequency.size(); ii++) 
	{
		summed_sp.frequency.push_back(sps[0].frequency[ii]);
		double sound_wave = 0;
		double turbulence = 0;
		double bubble_collision = 0;
		double total_amplitude = 0;

		for (int jj = 0; jj < sps.size(); jj++) 
		{
			sound_wave += sps[jj].sound_wave[ii];
			turbulence += sps[jj].turbulence[ii];
			bubble_collision += sps[jj].bubble_collision[ii];
			total_amplitude += sps[jj].total_amplitude[ii];
		}

		if (total_amplitude > peak_amplitude) 
		{
			peak_amplitude = total_amplitude;
			peak_frequency = sps[0].frequency[ii];
		}

		summed_sp.sound_wave.push_back(sound_wave);
		summed_sp.turbulence.push_back(turbulence);
		summed_sp.bubble_collision.push_back(bubble_collision);
		summed_sp.total_amplitude.push_back(total_amplitude);
	}

	summed_sp.peak_frequency = peak_frequency;
	summed_sp.peak_amplitude = peak_amplitude;
	const bool all_fit = std::all_of(sps.begin(), sps.end(),
		[](const GravWaveSpectrum &sp) { return sp.method == GravWaveMethod::FitFormulae; });
	if (sps.size() == 1)
	{
		summed_sp.snr = sps[0].snr;
	} else if (all_fit) {
		/* SNR is quadratic in Omega, so sum the spectra before integrating */
		summed_sp.snr = SNR_of([&](double f) {
			double o = 0.;
			for (const auto &sp : sps)
			{
				o += fit_omega(f, sp.alpha, sp.beta_H, sp.Tref, sp.Treh, sp.vw, sp.cs, sp.g_eff, sp.h_eff);
			}
			return o;
		});
	} else {
		summed_sp.snr = get_SNR_tabulated(summed_sp.frequency, summed_sp.total_amplitude);
	}
	add_noise_curves(summed_sp);

	return summed_sp;
}

double 
GravWaveCalculator::fit_omega(double f, double alpha, double beta_H, double T_ref, double T_reh, double vw, double cs, double g_eff, double h_eff) const 
{
	if(use_legacy_gw_methods)
	{
		double sound_wave = GW_sound_wave_legacy(f, alpha, beta_H, T_ref, g_eff, vw, cs);
		double turbulence = GW_turbulence_legacy(f, alpha, beta_H, T_ref, g_eff, vw, cs);
		double bubble_collision = to_GeV(T_ref) < T_threshold_bubble_collision ? GW_bubble_collision_legacy(f, alpha, beta_H, T_ref, g_eff, vw, cs) : 0;
		return sound_wave + turbulence + bubble_collision;
	} else {
		double sound_wave = GW_sound_wave(f, alpha, beta_H, T_ref, T_reh, vw, cs, g_eff, h_eff);
		double turbulence = GW_turbulence(f, alpha, beta_H, T_ref, T_reh, vw, cs, g_eff, h_eff);
		double bubble_collision = to_GeV(T_ref) < T_threshold_bubble_collision ? GW_bubble_collision(f, alpha, beta_H, T_ref, T_reh, vw, cs, g_eff, h_eff) : 0;
		return sound_wave + turbulence + bubble_collision;
	}
}

std::vector<SNRResult> 
GravWaveCalculator::get_SNR_tabulated(const std::vector<double> &frequency, const std::vector<double> &omega) const 
{
	if (frequency.size() != omega.size()) 
	{
		throw std::runtime_error("Frequency and amplitude grids differ in length - cannot compute SNR");
	}

	if (frequency.size() < 2) 
	{
		throw std::runtime_error("Need at least two frequency points to compute SNR");
	}

	std::vector<double> f_asc(frequency);
	std::vector<double> o_asc(omega);

	if (f_asc.front() > f_asc.back()) 
	{
		std::reverse(f_asc.begin(), f_asc.end());
		std::reverse(o_asc.begin(), o_asc.end());
	}

	std::vector<SNRResult> results;
	for (const auto &d : detector_list) 
	{
		results.push_back(SNR_in(*d, f_asc, o_asc));
	}
	return results;
}

SNRResult 
GravWaveCalculator::SNR_in(const Detector &d, const std::vector<double> &f_asc, const std::vector<double> &o_asc) const 
{
	SNRResult r;
	r.detector = d.name();
	r.threshold = d.threshold();
	r.T_obs = d.T_obs();
	r.f_min = d.f_min();
	r.f_max = d.f_max();
	r.band_covered = f_asc.front() <= d.f_min() * (1 + 1e-9) && f_asc.back() >= d.f_max() * (1 - 1e-9);

	if (!r.band_covered) 
	{
		LOG(warning) << "Spectrum covers [" << f_asc.front() << ", " << f_asc.back()
					<< "] Hz but the " << d.name() << " SNR band is [" << d.f_min() << ", " << d.f_max()
					<< "] Hz. The SNR is integrated over the overlap only, so it is a lower bound.";
	}

	/* Integrate in u = ln f, where the grid is (usually) uniform: int (Omega/N)^2 df = int f (Omega/N)^2 du */
	std::vector<double> u, g;
	for (size_t ii = 0; ii < f_asc.size(); ii++) 
	{
		const double f = f_asc[ii];
		if (f <= 0. || f < d.f_min() * (1 - 1e-9) || f > d.f_max() * (1 + 1e-9)) 
		{
			continue;
		}
		const double o = o_asc[ii] > 0. ? o_asc[ii] : 0.;
		const double n = d.omega_noise(f);
		u.push_back(std::log(f));
		g.push_back(f * o * o / (n * n));
	}

	if (u.size() < 2) 
	{
		LOG(warning) << "Fewer than two spectrum frequencies lie in the " << d.name() << " SNR band [" 
					<< d.f_min() << ", " << d.f_max() << "] Hz. Setting the SNR to zero.";
		r.snr = 0.;
	} else {
		const double T_obs_s = d.T_obs() * 365.25 * 86400;
		r.snr = std::sqrt(d.ndet() * simpson(u, g) * T_obs_s);
	}
	r.detected = r.snr >= r.threshold;
	return r;
}

std::vector<SNRResult>
GravWaveCalculator::get_SNR(double alpha, double beta_H, double T_ref, double T_reh, double vw, double cs, double g_eff, double h_eff) const
{
	return SNR_of([&](double f) { return fit_omega(f, alpha, beta_H, T_ref, T_reh, vw, cs, g_eff, h_eff); });
}

std::vector<SNRResult>
GravWaveCalculator::SNR_of(const std::function<double(double)> &omega_of) const
{
	/* Tabulate the spectrum separately on each detector's band, so every band gets the same resolution */
	std::vector<SNRResult> results;
	for (const auto &d : detector_list) 
	{
		const double log_lo = std::log10(d->f_min());
		const double log_hi = std::log10(d->f_max());
		const int n = std::max(3, static_cast<int>(std::ceil((log_hi - log_lo) * SNR_steps_per_decade)) + 1);

		std::vector<double> frequency(n), omega(n);
		for (int ii = 0; ii < n; ii++) 
		{
			frequency[ii] = std::pow(10., log_lo + ii * (log_hi - log_lo) / (n - 1));
			omega[ii] = omega_of(frequency[ii]);
		}
		results.push_back(SNR_in(*d, frequency, omega));
	}
	return results;
}

void 
GravWaveCalculator::write_spectrum_to_text(const GravWaveSpectrum &sp, const std::string &filename) const 
{
	std::ofstream file(filename);
	file << "frequency,total_amplitude,sound_wave,turbulence,bubble_collision";
	for (const auto &r : sp.snr) 
	{
		file << "," << r.detector << "_noise";
	}
	file << std::endl;

	for (int ii = 0; ii < sp.frequency.size(); ii++) 
	{
		file << sp.frequency[ii] << "," << sp.total_amplitude[ii] << "," 
		<< sp.sound_wave[ii] << "," << sp.turbulence[ii] << "," 
		<< sp.bubble_collision[ii];
		for (const auto &r : sp.snr) 
		{
			file << "," << r.noise[ii];
		}
		file << std::endl;
	}

  	LOG(debug) << "GW spectrum has been written to " << filename;
}

void 
GravWaveCalculator::write_spectrum_to_text(int i, const std::string &filename) const 
{
  	write_spectrum_to_text(spectrums[i], filename);
}

void 
GravWaveCalculator::write_spectrum_to_text(const std::string &filename) const 
{
  LOG(debug) << "writing " << spectrums.size() << "spectrums to text at " << filename;

  for (int ii = 0; ii < spectrums.size(); ii++) 
  {
    write_spectrum_to_text(spectrums[ii], std::to_string(ii) + "_" + filename);
  }
}

void 
GravWaveCalculator::check_detectors() const 
{
	for (size_t ii = 0; ii < detector_list.size(); ii++) 
	{
		check_detector(*detector_list[ii]);
		for (size_t jj = 0; jj < ii; jj++) 
		{
			if (detector_list[jj]->name() == detector_list[ii]->name()) 
			{
				throw std::invalid_argument("Two detectors are named " + detector_list[ii]->name());
			}
		}
	}
}

const TransitionMilestone 
*GravWaveCalculator::milestone_of(const ThermalParameterSet &tps) const 
{
	const TransitionMilestone *milestone = nullptr;
	switch (default_milestone) {
		case MilestoneType::ONSET:
			milestone = &tps.onset;
			break;
		case MilestoneType::PERCOLATION:
			milestone = &tps.percolation;
			break;
		case MilestoneType::COMPLETION:
			milestone = &tps.completion;
			break;
		case MilestoneType::NUCLEATION:
			milestone = &tps.nucleation;
			break;
		default:
			LOG(debug) << "Invalid milestone type. GW will not be calculated !";
			return nullptr;
	}

	if (milestone->status != MilestoneStatus::YES) 
	{
		LOG(debug) << "No " << static_cast<int>(milestone->type) << " milestone found for transition with TC = " << tps.TC;
		return nullptr;
	}

	LOG(debug) << "Found " << static_cast<int>(milestone->type) << " milestone with T = " << milestone->temperature;
	return milestone;
}

void 
GravWaveCalculator::add_fit_contributions(GravWaveSpectrum &sp, double alpha_fit) const 
{
	const bool use_collision = to_GeV(sp.Tref) < T_threshold_bubble_collision;

	sp.alpha_fit = alpha_fit;
	sp.turbulence.assign(sp.frequency.size(), 0.);
	sp.bubble_collision.assign(sp.frequency.size(), 0.);

	for (size_t ii = 0; ii < sp.frequency.size(); ii++) 
	{
		const double f = sp.frequency[ii];
		if (use_legacy_gw_methods)
		{
			sp.turbulence[ii] = GW_turbulence_legacy(f, alpha_fit, sp.beta_H, sp.Tref, sp.g_eff);
			if (use_collision) 
			{
				sp.bubble_collision[ii] = GW_bubble_collision_legacy(f, alpha_fit, sp.beta_H, sp.Tref, sp.g_eff);
			}
		} else {
			sp.turbulence[ii] = GW_turbulence(f, alpha_fit, sp.beta_H, sp.Tref, sp.Treh, sp.vw, sp.cs, sp.g_eff, sp.h_eff);
			if (use_collision) 
			{
				sp.bubble_collision[ii] = GW_bubble_collision(f, alpha_fit, sp.beta_H, sp.Tref, sp.Treh, sp.vw, sp.cs, sp.g_eff, sp.h_eff);
			}
		}
	}

	if (!use_collision) 
	{
		LOG(debug) << "Tref = " << sp.Tref << " is above T_threshold_bubble_collision = "
			<< T_threshold_bubble_collision << ", so only turbulence was added.";
	}
}

void 
GravWaveCalculator::finalise_spectrum(GravWaveSpectrum &sp) const 
{
	const size_t n = sp.frequency.size();
	if (sp.turbulence.size() != n) sp.turbulence.assign(n, 0.);
	if (sp.bubble_collision.size() != n) sp.bubble_collision.assign(n, 0.);

	sp.total_amplitude.assign(n, 0.);

	double peak_frequency = 0.;
	double peak_amplitude = 0.;

	for (size_t ii = 0; ii < n; ii++) 
	{
		sp.total_amplitude[ii] = sp.sound_wave[ii] + sp.turbulence[ii] + sp.bubble_collision[ii];
		if (sp.total_amplitude[ii] > peak_amplitude) 
		{
			peak_amplitude = sp.total_amplitude[ii];
			peak_frequency = sp.frequency[ii];
		}
	}

	sp.peak_frequency = peak_frequency;
	sp.peak_amplitude = peak_amplitude;
	sp.snr = get_SNR_tabulated(sp.frequency, sp.total_amplitude);
	add_noise_curves(sp);
}

void 
GravWaveCalculator::add_noise_curves(GravWaveSpectrum &sp) const 
{
	for (auto &r : sp.snr) 
	{
		const auto d = std::find_if(detector_list.begin(), detector_list.end(),
			[&r](const std::shared_ptr<const Detector> &det) { return det->name() == r.detector; });
		if (d == detector_list.end()) 
		{
			throw std::runtime_error("SNR result for unknown detector " + r.detector);
		}
		r.noise.resize(sp.frequency.size());
		for (size_t ii = 0; ii < sp.frequency.size(); ii++) 
		{
			r.noise[ii] = (*d)->omega_noise(sp.frequency[ii]);
		}
	}
}

double
GravWaveCalculator::get_RH(const double& betaH) const
{
	constexpr double f_perc = 0.28957;
	const double cs = std::sqrt(1./3.);
	return std::cbrt(8 * M_PI / f_perc) * std::max(vw, cs) / betaH;
}

std::array<double, 1> 
GravWaveCalculator::get_collision_peaks(const double &H0_star, const double &RH) const
{
	const double f_col = 0.49 * H0_star/RH;
	return {f_col};
}

std::array<double, 2> 
GravWaveCalculator::get_sound_wave_peaks(const double& H0_star, const double& RH, const double& vw, const double& cs) const
{
	const double Delta_w = std::abs(vw - cs)/std::max(vw, cs);
	const double f_sw_1 = 0.2 * H0_star / RH;
	const double f_sw_2 = 0.5/Delta_w * H0_star / RH;
	return {f_sw_1, f_sw_2};
}

std::array<double, 3> 
GravWaveCalculator::get_turbulence_peaks(const double& H0_star, const double& RH, const double& K_turb) const
{
	const double f_turb_1 = H0_star;
	const double f_turb_2 = 2.2 * H0_star/RH;
	const double f_turb_3 = std::sqrt(3.*K_turb)/4. * H0_star/RH;
	return {f_turb_1, f_turb_2, f_turb_3};
}

double
GravWaveCalculator::get_K(const double& A, const double& kappa, const double& alpha) const
{
	return A * kappa * alpha/(1+alpha);
}

double 
GravWaveCalculator::get_kappa_sw(const double& alpha, const double& vw, const double& cs) const
{
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

double 
GravWaveCalculator::get_kappa_turb(const double& alpha, const double& vw, const double& cs) const 
{
	return epsilon * get_kappa_sw(alpha, vw, cs);
}

double 
GravWaveCalculator::get_kappa_col(const double& alpha) const 
{
	return 1 / (1 + 0.715 * alpha) * (0.715 * alpha + 4. / 27 * sqrt(3 * alpha / 2));
}

double 
GravWaveCalculator::get_prefactor(const double& g_eff, const double& h_eff) const
{
	const double term1 = omega_hsq_neutrino * std::pow(D, -4./3.);
	const double term2 = std::pow(h_0/h_eff, 4./3.);
	const double term3 = g_eff/g_0;
	return term1*term2*term3;
}

double 
GravWaveCalculator::get_Hubble_rate_today(const double& T, const double& g_eff, const double& h_eff) const
{
	const double term1 = (11.2e-9) * std::pow(D, -1./3.);
	const double term2 = to_GeV(T) / (100 * 1e-3);
	const double term3 = std::pow(g_eff/10, 1./2.);
	const double term4 = std::pow(10/h_eff, 1./3.);

	return term1*term2*term3*term4;
}

double
GravWaveCalculator::get_sound_wave_N(const std::array<double, 2>& peak_freqs) const
{
	auto [f_sw_1, f_sw_2] = peak_freqs;
	const double N_num = f_sw_1 * (f_sw_1*f_sw_1*f_sw_1*f_sw_1 + f_sw_2*f_sw_2*f_sw_2*f_sw_2);
	const double N_den = f_sw_2*f_sw_2*f_sw_2 * (std::sqrt(2)*f_sw_1*f_sw_1 - 2.*f_sw_1*f_sw_2 + std::sqrt(2)*f_sw_2*f_sw_2);
	const double N = 4./M_PI * N_num/N_den;

	return N;
}

double
GravWaveCalculator::get_sound_wave_Y(const double& RH, const double& K) const
{
	return std::min(1., 2 * RH / std::sqrt(3 * K));
}

double 
GravWaveCalculator::singly_broken_power_law(
	const double& f, const double& f0, const double& f1, 
	const double& n0, const double& n1, 
	const double& a1, const double& b1) const
{
	const double term0 = std::pow(f/f0, n0);
	const double term1 = std::pow(b1 + b1*std::pow(f/f1, a1), n1);
	return term0 * term1;
}

double 
GravWaveCalculator::doubly_broken_power_law(
	const double& f, const double& f0, const double& f1, const double& f2,
	const double& n0, const double& n1, const double& n2, 
	const double& a1, const double& a2) const
{
	const double term0 = std::pow(f/f0, n0);
	const double term1 = std::pow(1 + std::pow(f/f1, a1), n1);
	const double term2 = std::pow(1 + std::pow(f/f2, a2), n2);
	return term0 * term1 * term2;
}

double 
GravWaveCalculator::doubly_broken_power_law_LISA(
	const double& f, const double& f1, const double& f2,
	const double& n1, const double& n2, const double& n3, 
	const double& a1, const double& a2) const
{
	const double term0 = std::pow(f/f1, n1);
	const double term1 = std::pow(1 + std::pow(f/f1, a1), (-n1+n2)/a1);
	const double term2 = std::pow(1 + std::pow(f/f2, a2), (-n2+n3)/a2);
	return term0 * term1 * term2;
}

} // namespace PhaseTracer
