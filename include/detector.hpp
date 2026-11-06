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

#ifndef PHASETRACER_DETECTOR_HPP_
#define PHASETRACER_DETECTOR_HPP_

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

namespace PhaseTracer {

/**
 * @brief A gravitational wave detector used in the SNR calculation.
 */
class Detector {
public:
    virtual ~Detector() = default;

	/** @brief Detector name in plain English. */
    virtual std::string name() const = 0;

	/** @brief SNR^2 = ndet * T_obs * int (Omega/Omega_n)^2 df: 2 for a cross-correlated network, 1 for a single detector. */
    virtual double ndet() const { return 1.0; }

	/** @brief Detector noise in units of omega h^2. */
    virtual double omega_noise(double f) const = 0;  

	/** @brief SNR threshold for detection. */
	virtual double threshold() const { return 10.; }

	/** @brief Observation time in years. */
	virtual double T_obs() const = 0;

	/** @brief Lower bound of frequency band. */
	virtual double f_min() const = 0;

	/** @brief Upper bound of frequency band. */
	virtual double f_max() const = 0;

	friend std::ostream &operator<<(std::ostream &o, const Detector &d);
};

/** @brief Robson, Cornish & Liu noise model, used for LISA and Taiji */
class RCL : public Detector {
public:

	RCL(double L, double f_star, double P_oms, double P_acc, double threshold, double t_obs, const std::string &name)
		: L_(L), f_star_(f_star), P_oms_(P_oms), P_acc_(P_acc), threshold_(threshold), t_obs_(t_obs), name_(name) {}

	RCL(double L, double P_oms, double P_acc, double threshold, double t_obs, const std::string &name)
		: L_(L), f_star_(2.99792458e8 / (2 * M_PI * L)), P_oms_(P_oms), P_acc_(P_acc), threshold_(threshold), t_obs_(t_obs), name_(name) {}

    std::string name() const override { return name_; }

	double omega_noise(double f) const override 
	{
		/* Robson, Cornish & Liu (arXiv:1803.01944) */
		const double S_oms = P_oms_ * P_oms_ * (1 + pow(2e-3 / f, 4));
		const double S_acc = P_acc_ * P_acc_ * (1 + pow(0.4e-3 / f, 2)) * (1 + pow(f / 8e-3, 4));
		const double S_n = 10. / (3 * L_ * L_) * (S_oms + 2 * (1 + pow(std::cos(f / f_star_), 2)) * S_acc / pow(2 * M_PI * f, 4)) * (1 + 0.6 * pow(f / f_star_, 2));

		double S_c = 0.;
		if (f < 1e-2) {
			S_c = 9e-45 * pow(f, -7. / 3) * std::exp(-pow(f, 0.138) - 221 * f * std::sin(521 * f)) * (1 + std::tanh(1680 * (1.13e-3 - f)));
		}

		const double H_100 = 100. / 3.0857e19;
		return 2 * M_PI * M_PI / (3 * H_100 * H_100) * pow(f, 3) * (S_n + S_c);
	}

	double threshold() const override { return threshold_; }

	double T_obs() const override { return t_obs_; }

	double f_min() const override { return 1e-5; }

	double f_max() const override { return 1e0; }

private:
	double L_; // arm length in metres
	double f_star_; // characteristic frequency
	double P_oms_; // optical metrology system noise
	double P_acc_; // acceleration noise
	double t_obs_; // observation time in years
	double threshold_; // snr threshold
	std::string name_; // name

	static constexpr double days_to_seconds = 24.*3600.;
	static constexpr double years_to_seconds = 365.25*24.*3600.;
	static constexpr double H0_conv = (3.241e-18); // 100 km/s/Mpc in 1/s
};

/** @brief Generic PTA class, used for SKA, NANAOGrav, and EPTA */
class PTA : public Detector {
public:

	PTA(int N_pulsars, double t_obs, double sampling_rate, double sigma, double threshold, const std::string &name)
		: N_pulsars_(N_pulsars), t_obs_(t_obs), sampling_rate_(sampling_rate), sigma_(sigma), threshold_(threshold), name_(name) {}

    std::string name() const override { return name_; }

	double omega_noise(double f) const override 
	{
		double pref_1 = std::sqrt(2. / (N_pulsars_ * (N_pulsars_-1)));
		double pref_2 = 64. * std::sqrt(3.) * M_PI*M_PI*M_PI*M_PI * sigma_*sigma_ * (sampling_rate_*days_to_seconds) / (H0_conv*H0_conv);
		return pref_1 * pref_2 * f*f*f*f*f;
	}

	/* Cross-correlated pulsar pairs */
	double ndet() const override { return 2.; }

	double threshold() const override { return threshold_; }

	double T_obs() const override { return t_obs_; }

	double f_min() const override { return 1.0 / t_obs_ / years_to_seconds; }

	double f_max() const override { return 1.0 / sampling_rate_ / days_to_seconds; }

private:
	int N_pulsars_; // number of observed pulsars
	double t_obs_; // observation time in years
	double sampling_rate_;  // observation interval in days
	double sigma_; // timing uncertainty in seconds
	double threshold_; // snr threshold
	std::string name_; // name

	static constexpr double days_to_seconds = 24.*3600.;
	static constexpr double years_to_seconds = 365.25*24.*3600.;
	static constexpr double H0_conv = 100. / 3.0857e19; // 100 km/s/Mpc in 1/s
};

/**
 * @brief Base for the space interferometers of ArXiv:1811.11175 (B-DECIGO, DECIGO, BBO).
 *
 * Converts a strain noise S_eff(f) in 1/Hz to h^2 Omega. The h^2 is absorbed by using H_100.
 */
class SpaceInterferometer : public Detector {
public:
	double omega_noise(double f) const override 
	{ 
		const double H_100 = 100. / 3.0857e19;
		return 2. * M_PI * M_PI / (3. * H_100 * H_100) * f*f*f * S_eff(f);
	}

	double T_obs() const override { return 4.; }

protected:
	/** Effective strain noise in 1/Hz */
	virtual double S_eff(double f) const = 0;

	static constexpr double years_to_seconds = 365.25*24.*3600.;
};

/**
 * @brief Base for DECIGO and BBO: instrument noise plus the white dwarf binary foregrounds of 1811.11175.
 *
 * S_eff = 5 [min(S_inst exp(kappa T^-1 dN/df), S_inst + S_gal) + S_exgal]
 */
class WDForegroundInterferometer : public SpaceInterferometer {
public:
	/* Cross-correlated detector network */
	double ndet() const override { return 2.; }

protected:
	/** Instrument noise in 1/Hz */
	virtual double S_inst(double f) const = 0;

	double S_eff(double f) const override
	{
		const double kappa = 4.5;
		const double dNdf = 2e-3 * std::pow(f, -11./3.);
		const double T = T_obs() * years_to_seconds;
		const double S_inst_ = S_inst(f);
		/* At low f the exponential overflows to inf, so the min correctly picks S_inst + S_gal */
		return 5. * (std::min(S_inst_ * std::exp(kappa * dNdf / T), S_inst_ + S_gal(f)) + S_exgal(f));
	}

private:
	static double cutoff(double f) { return std::exp(-2. * (f / 0.05) * (f / 0.05)); }

	/** Galactic white dwarf binary foreground */
	static double S_gal(double f) { return 2.1e-45 * std::pow(f, -7./3.) * cutoff(f); }

	/** Extragalactic white dwarf binary foreground */
	static double S_exgal(double f) { return 4.2e-47 * std::pow(f, -7./3.) * cutoff(f); }
};

/** @brief B_DECIGO, based on ArXiv:1811.11175. */
class B_DECIGO : public SpaceInterferometer {
public:
	std::string name() const override { return "B_DECIGO"; }

	double threshold() const override { return 8.; }

	double f_min() const override { return 1e-2; }

	double f_max() const override { return 1e2; }

protected:
	double S_eff(double f) const override
	{
		double A = 2.02e-45;
		double B = 1.584e-2;
		double C = 1.584e-3;

		double invf = 1. / f;
		double invf2 = invf * invf;

		return A * (1  + B * invf2*invf2 + C * f*f);
	}
};

/** @brief BBO, based on ArXiv:1811.11175. */
class BBO : public WDForegroundInterferometer {
public:
	std::string name() const override { return "BBO"; }

	double threshold() const override { return 10.; }

	double f_min() const override { return 1e-3; }

	double f_max() const override { return 1e2; }

protected:
	double S_inst(double f) const override
	{
		double A = 2.9e-49;
		double B = 9.2e-52;
		double C = 1.8e-49;

		double invf = 1. / f;
		double inv2f = invf * invf;
		double inv4f = inv2f * inv2f;

		return A + B * inv4f + C * f*f;
	}
};

/** @brief DECIGO, based on ArXiv:1811.11175. */
class DECIGO : public WDForegroundInterferometer {
public:
	std::string name() const override { return "DECIGO"; }

	double threshold() const override { return 10.; }

	double f_min() const override { return 1e-3; }

	double f_max() const override { return 1e2; }

protected:
	double S_inst(double f) const override
	{
		double fp = 7.36;
		double f_fp = f/fp;

		double A = 5.3e-48;
		double B = 2.3e-7;
		double C = 2.6e-8;

		double invf_fp = 1. / f_fp;
		double inv2f_fp  = invf_fp  * invf_fp;
		double inv4f_fp  = inv2f_fp  * inv2f_fp;

		return A * ( (1 + f_fp*f_fp) + B * inv4f_fp * 1./(1. + f_fp*f_fp) + C * inv4f_fp);
	}
};

/**
 * @brief Check that a detector is physically sensible.
 * @throws std::invalid_argument naming the detector and the problem.
 */
void check_detector(const Detector &d);

/**
 * @struct SNRResult
 * @brief Signal-to-noise ratio of one spectrum in one detector.
 */
struct SNRResult
{
	/** @brief Name of the detector. */
	std::string detector;

	/** @brief Signal-to-noise ratio. */
	double snr = std::numeric_limits<double>::quiet_NaN();

	/** @brief SNR above which the signal counts as detected. */
	double threshold = std::numeric_limits<double>::quiet_NaN();

	/** @brief Whether snr >= threshold. */
	bool detected = false;

	/** @brief Observation time in years. */
	double T_obs = std::numeric_limits<double>::quiet_NaN();

	/** @brief SNR integration band of the detector in Hz. */
	double f_min = std::numeric_limits<double>::quiet_NaN();
	double f_max = std::numeric_limits<double>::quiet_NaN();

	/** @brief Whether the spectrum spans [f_min, f_max]. If false, snr is a lower bound. */
	bool band_covered = false;

	/** @brief Noise curve Omega h^2 on the spectrum's frequency grid. Empty if there is no spectrum grid. */
	std::vector<double> noise;

	/** @brief Pretty-printer for SNRResult */
	friend std::ostream &operator<<(std::ostream &o, const SNRResult &r);
};

} // namespace PhaseTracer

#endif // PHASETRACER_DETECTOR_HPP_
