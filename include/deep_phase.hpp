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
#include "thermal_parameters.hpp"

namespace PhaseTracer {

std::vector<double> get_ssm_amplitude(ThermalParams tps, double vw, double dof, double min_frequency, double max_frequency, double num_frequency) {

	// created pt params object
	double alpha = tps.alpha_tp;
	double beta = tps.beta_tp;
	double Tref = tps.TP;
	double cs_true_tp = 1/sqrt(3);
	double cs_false_tp = 1/sqrt(3);
	double dtau = 1/beta;
	double wN = 1.71;
	std::string PTmodel = "bag";
	std::string nuc_type = "exp";
	PhaseTransition::Universe un;
	PhaseTransition::PTParams paramsPT(vw, alpha, beta, dtau, wN, PTmodel, nuc_type, un);
	double Rs = paramsPT.Rs();

	// red shift factor
	double a_shift = pow(dof/3.014, 1./3.) * (Tref/2.35e-13); // remove magic numbers
	double fac = 2 * M_PI * Rs * a_shift;

	double Kmin = min_frequency * fac;
	double Kmax = max_frequency * fac;
	std::vector<double> momentumVec = logspace(Kmin, Kmax, num_frequency);

  Spectrum::PowerSpec OmegaGW = Spectrum::GWSpec(momentumVec, paramsPT);

	return OmegaGW.P();
}

} // namespace PhaseTracer

#endif // PHASETRACER_DP_HPP_