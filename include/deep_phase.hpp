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

#ifndef PHASETRACER_DP_HPP_
#define PHASETRACER_DP_HPP_

#include "DeepPhase/include/deepphase.hpp"
#include "DeepPhase/include/maths_ops.hpp"
#include "thermal_parameters.hpp"
#include "logger.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <streambuf>

namespace PhaseTracer {

// Custom stream buffer that redirects output to LOG(debug)
class LogStreamBuffer : public std::streambuf {
private:
	std::ostringstream buffer_;
    
protected:
	int overflow(int c) override {
		if (c != EOF) {
			buffer_ << static_cast<char>(c);
			// If we encounter a newline, flush the buffer to LOG(debug)
			if (c == '\n') {
				sync();
			}
		}
		return c;
	}
	
	int sync() override {
		std::string line = buffer_.str();
		if (!line.empty()) {
			// Remove trailing newline for cleaner logging
			if (line.back() == '\n') {
				line.pop_back();
			}
			if (!line.empty()) {
				LOG(debug) << line;
			}
			buffer_.str("");
			buffer_.clear();
		}
		return 0;
	}
};

class GW_DeepPhase {

public:

  GW_DeepPhase() = default;
  ~GW_DeepPhase() = default;

	std::vector<double> get_ssm_amplitude(ThermalParams tps, double vw, double dof, double min_frequency, double max_frequency, int num_frequency) {

		// Redirect stdout to LOG(debug) instead of suppressing it
		std::streambuf* orig_buf = std::cout.rdbuf();
		LogStreamBuffer log_buffer;
		
		try {
			std::cout.rdbuf(&log_buffer);
			
			auto params = get_pt_params(tps, vw, dof);
			params.print();

			double Hconf = params.un().Hs() /(std::pow(params.un().gs() / params.un().g0(), 1./3.) * params.un().Ts() / params.un().T0());
			double R0 = 1./(1.65e-5 / (params.Rs() * Hconf) * (params.un().Ts() / 100) * std::pow(params.un().gs() / 100, 1./6.));
			std::cout << "R0 = " << R0 << "\n";
			std::cout << "Hconf = " << Hconf << "\n";
			// R0 = 1; // TODO

			double kmin = 2*M_PI*min_frequency;
			double kmax = 2*M_PI*max_frequency;

			double Kmin = kmin*R0; std::cout << "Kmin = " << Kmin << "\n";
			double Kmax = kmax*R0; std::cout << "Kmax = " << Kmax << "\n";
			std::vector<double> momentumVec = logspace(Kmin, Kmax, static_cast<std::size_t>(num_frequency));

			Spectrum::PowerSpec OmegaGW = Spectrum::GWSpec2(momentumVec, params);

			// OmegaGW.write("gw_spec.csv");

			// Restore stdout
			std::cout.rdbuf(orig_buf);
			// Flush any remaining content in the log buffer
			log_buffer.pubsync();

			return OmegaGW.P();
			
		} catch (...) {
			// Ensure stdout is restored even if an exception occurs
			std::cout.rdbuf(orig_buf);
			log_buffer.pubsync();
			throw;
		}
	}

	std::pair<double, double> get_ssm_max(const std::vector<double> amp, const std::vector<double> freq) {
		double max_freq = 0;
		double max_amp = 0;
		for (size_t i = 0; i < amp.size(); i++) {
			if( amp[i] > max_amp ) {
				max_amp = amp[i];
				max_freq = freq[i];
			}
		}
		return {max_amp, max_freq};
	}

private :

	PhaseTransition::Universe get_universe(ThermalParams tps, double dof) {
		if ( tps.nucleates == MilestoneStatus::YES ) {
			// Validate inputs
			if (tps.TN <= 0 || tps.H_tn <= 0) {
				throw std::invalid_argument("Invalid thermal parameters for universe creation");
			}
			return PhaseTransition::Universe(tps.TN, tps.H_tn, dof);
		} else if (tps.percolates == MilestoneStatus::YES) {
			// Validate inputs
			if (tps.TP <= 0 || tps.H_tp <= 0) {
				throw std::invalid_argument("Invalid thermal parameters for universe creation");
			}
			return PhaseTransition::Universe(tps.TP, tps.H_tp, dof);
		}
		PhaseTransition::Universe un;
		return un;
	}

	PhaseTransition::PTParams get_pt_params(ThermalParams tps, double vw, double dof) {
		if (tps.nucleates == MilestoneStatus::YES) {
			// Validate inputs
			if (tps.alpha_tn <= 0 || tps.betaH_tn <= 0 || tps.beta_tn <= 0) {
				throw std::invalid_argument("Invalid thermal parameters for PTParams creation");
			}
			double alpha = tps.alpha_tn;
			double betaH = tps.betaH_tn;
			double Tref = tps.TN;
			double wN = tps.we_tn;
			PhaseTransition::Universe un = get_universe(tps, dof);
			return PhaseTransition::PTParams(vw, alpha, betaH, 0.001, wN, "exp", un);
		} else if (tps.percolates == MilestoneStatus::YES) {
			// Validate inputs
			if (tps.alpha_tp <= 0 || tps.betaH_tp <= 0 || tps.beta_tp <= 0) {
				throw std::invalid_argument("Invalid thermal parameters for PTParams creation");
			}
			double alpha = tps.alpha_tp;
			double betaH = tps.betaH_tp;
			double Tref = tps.TP;
			double wN = tps.we_tp;
			PhaseTransition::Universe un = get_universe(tps, dof);
			return PhaseTransition::PTParams(vw, alpha, betaH, 0.001, wN, "exp", un);
		}
		PhaseTransition::PTParams pt;
		return pt;
	}

};

} // namespace PhaseTracer

#endif // PHASETRACER_DP_HPP_