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

#ifndef HYDROGRAV_INTERFACE_HPP_INCLUDED
#define HYDROGRAV_INTERFACE_HPP_INCLUDED

#include "HydroGrav/include/hydrograv.hpp"
#include "HydroGrav/include/maths.hpp"
#include "thermo_finder.hpp"
#include "equation_of_state.hpp"
#include "gravwave_calculator.hpp"
#include "transition_milestones.hpp"
#include "logger.hpp"
#include <string>

namespace HydroGrav {
namespace config = ::config;
namespace PhaseTransition = ::PhaseTransition;
namespace Hydrodynamics = ::Hydrodynamics;
namespace Spectrum = ::Spectrum;

enum class EoSModel {
    BAG,
	MUNU,
    VEFF,
	ALL
};

}  // namespace HydroGrav

namespace PhaseTracer {

/**
 * Translators between PhaseTracer and HydroGrav.
 *
 * @warning HydroGrav's PowerSpec and FluidProfile both keep a non-owning
 * pointer to the PTParams they were built from, and dereference it in write().
 * Whatever a caller does with a PowerSpec, the PTParams_Veff must outlive it.
 */
namespace HydroGravBridge {

/** @brief Translates equation of state. */
HydroGrav::PhaseTransition::EquationOfState to_hydrograv_eos(const PhaseTracer::EquationOfState &eos);

/** @brief Builds a universe from a TransitionMilestone. */
HydroGrav::PhaseTransition::Universe to_hydrograv_universe(const TransitionMilestone &milestone, double dof);

/** @brief Translates nucleation type. */
const char *to_hydrograv_nuc_type(NucleationType nucleation_type);

/** @brief Build PT_params. */
HydroGrav::PhaseTransition::PTParams_Veff to_hydrograv_pt_params(
	const ThermalParameterSet &tps, 
	const TransitionMilestone &milestone, 
	double vw, double dof);

/** @brief Translates fluid profiles. */
FluidProfile to_phasetracer_profile(const HydroGrav::Hydrodynamics::FluidProfile &fp);

/** @brief Translates the GW power spectrum. */
void fill_spectrum(GravWaveSpectrum &sp, const HydroGrav::Spectrum::PowerSpec &spec_in);

} // namespace HydroGravBridge

}  // namespace PhaseTracer


#endif // HYDROGRAV_INTERFACE_HPP_INCLUDED
