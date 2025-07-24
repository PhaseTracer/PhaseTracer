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

std::vector<double> get_ssm_amplitude(ThermalParams tps, double vw, double dof, double min_frequency, double max_frequency, int num_frequency) {

	// Validate inputs
	if (tps.beta_tp == 0.0) {
		throw std::invalid_argument("beta_tp cannot be zero");
	}
	if (min_frequency <= 0 || max_frequency <= 0 || min_frequency >= max_frequency) {
		throw std::invalid_argument("Invalid frequency range");
	}
	if (num_frequency < 2) {
		throw std::invalid_argument("num_frequency must be at least 2");
	}

	// Redirect stdout to LOG(debug) instead of suppressing it
	std::streambuf* orig_buf = std::cout.rdbuf();
	LogStreamBuffer log_buffer;
	
	try {
		std::cout.rdbuf(&log_buffer);
		
		// Create PT params object
		double alpha = tps.alpha_tp;
		double beta = tps.beta_tp;
		double Tref = tps.TP;
		double Hstar = tps.H_tp;
		double dtau = 1.0 / beta;
		double wN = 1.71;
		std::string PTmodel = "bag";
		std::string nuc_type = "exp";
		PhaseTransition::Universe un;
		PhaseTransition::PTParams paramsPT(vw, alpha, beta, dtau, wN, PTmodel, nuc_type, un);
		paramsPT.print();
		double Rs = paramsPT.Rs();

		double den = pow(8 * M_PI, 1./3.) * vw / tps.betaH_tp;
		std::cout << "den = den = " << den << "\n";

		double shift_factor = (1.65e-5) / (den) * (Tref/100) * pow(dof/100., 1./6.);
		// double shift_factor = (1.65e-5) / (Rs * Hstar) * (Tref/100) * pow(dof/100., 1./6.);

		double kmin = 2*M_PI*min_frequency;
		double kmax = 2*M_PI*max_frequency;

		double Kmin = kmin / shift_factor; std::cout << "Kmin = " << Kmin << "\n";
		double Kmax = kmax / shift_factor; std::cout << "Kmax = " << Kmax << "\n";
		std::vector<double> momentumVec = logspace(Kmin, Kmax, static_cast<std::size_t>(num_frequency));

		Spectrum::PowerSpec OmegaGW = Spectrum::GWSpec(momentumVec, paramsPT);

		OmegaGW.write("gw_spec.csv");

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

} // namespace PhaseTracer

#endif // PHASETRACER_DP_HPP_