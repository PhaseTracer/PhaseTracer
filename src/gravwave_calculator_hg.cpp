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

// The HydroGrav backend of GravWaveCalculator. Kept apart from
// gravwave_calculator.cpp so that the fit formula code never includes HydroGrav.

#ifdef BUILD_WITH_HG

#include <cmath>
#include <stdexcept>

#include "gravwave_calculator.hpp"
#include "hydrograv_interface.hpp"

namespace PhaseTracer {

GravWaveSpectrum GravWaveCalculator::calc_spectrum_ssm(const ThermalParameterSet &tps,
                                                       const TransitionMilestone &milestone) const {
  if (n_kRs_value < 2) {
    throw std::runtime_error("Number of kRs points must be greater than 1");
  }
  if (max_kRs_value < min_kRs_value) {
    throw std::runtime_error("max_kRs_value < min_kRs_value");
  }

  const auto pt_params = HydroGravBridge::to_hydrograv_pt_params(tps, milestone, vw, dof);

  LOG(debug) << "Calculating sound shell spectrum at T = " << milestone.temperature
             << ", alpha_munu = " << milestone.alpha_munu << ", vw = " << vw;

  const auto kRs_values = logspace(std::log10(min_kRs_value), std::log10(max_kRs_value), n_kRs_value);
  const auto spectrum = HydroGrav::Spectrum::GWSpec(kRs_values, pt_params);

  GravWaveSpectrum sp;
  sp.Tref = milestone.temperature;
  sp.alpha = milestone.alpha_munu;
  sp.beta_H = milestone.betaH;
  HydroGravBridge::fill_spectrum(sp, spectrum);

  return sp;
}

} // namespace PhaseTracer

#endif // BUILD_WITH_HG
