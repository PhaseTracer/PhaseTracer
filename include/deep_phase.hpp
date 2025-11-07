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
#include <variant>
#include <memory>

namespace PhaseTracer::DeepPhaseInterface {

enum class EoSModel {
    BAG,
	BAGext,
    VEFF,
	ALL
};

using PTParamsVariant = std::variant<PhaseTransition::PTParams_Bag, PhaseTransition::PTParams_Veff>;

inline std::pair<double, double>
obtain_peaks(const std::vector<double>& xVals, const std::vector<double>& yVals) 
{
	double xMax = 0;
	double yMax = 0;
	for (size_t i = 0; i < yVals.size(); ++i) 
	{
		if (yVals[i] > yMax) 
		{
			yMax = yVals[i];
			xMax = xVals[i];
		}
	}
	return {xMax, yMax};
}

inline
PhaseTransition::EquationOfState
phasetracer_EoS_to_deepphase_EoS(PhaseTracer::EoS eos)
{
	std::vector<double> t_vals = eos.temp;
	std::vector<double> ps_vals = eos.pressure_false;
	std::vector<double> pb_vals = eos.pressure_true;
	std::vector<double> es_vals = eos.energy_false;
	std::vector<double> eb_vals = eos.energy_true;

	PhaseTransition::EquationOfState output(t_vals, ps_vals, pb_vals, es_vals, eb_vals);
	return output;
}

inline PhaseTransition::Universe
get_universe_from_thermal_params(const ThermalParams& tps, const double& dof = 107.75) 
{
	if ( tps.nucleates == MilestoneStatus::YES ) 
	{
		if (tps.TN <= 0 || tps.H_tn <= 0) 
		{
			throw std::invalid_argument("Invalid thermal parameters for universe creation");
		}
		return PhaseTransition::Universe(tps.TN, dof, tps.H_tn);
	} else if (tps.percolates == MilestoneStatus::YES) 
	{
		if (tps.TP <= 0 || tps.H_tp <= 0) 
		{
			throw std::invalid_argument("Invalid thermal parameters for universe creation");
		}
		return PhaseTransition::Universe(tps.TP, dof, tps.H_tp);
	}
	PhaseTransition::Universe un;
	return un;
}

inline PTParamsVariant
get_pt_params_from_thermal_params(
	const ThermalParams& tps,
	const PhaseTransition::Universe& un,
	const double& vw,
	const double& dtauRs,
	EoSModel model,
	const PhaseTransition::EquationOfState* eos = nullptr)
{

	if (tps.nucleates != MilestoneStatus::YES && tps.percolates != MilestoneStatus::YES) {
		throw std::invalid_argument("No valid thermal parameters: neither nucleates nor percolates");
	}

	const bool use_nucleation = (tps.nucleates == MilestoneStatus::YES);

	const double alpha = use_nucleation ? tps.alpha_tn : tps.alpha_tp;
	const double betaH = use_nucleation ? tps.betaH_tn : tps.betaH_tp;
	const double H = use_nucleation ? tps.H_tn : tps.H_tp;
	const double Tref = use_nucleation ? tps.TN : tps.TP;
	const double cs_true = use_nucleation ? tps.cs_true_tn : tps.cs_true_tp;
	const double cs_false = use_nucleation ? tps.cs_false_tn : tps.cs_false_tp;

	if (alpha <= 0 || betaH <= 0) {
		throw std::invalid_argument("Invalid thermal parameters: alpha or betaH <= 0");
	}

	const double beta = betaH * H;
	const double Rs = std::pow(8. * M_PI, 1./3.) * vw / beta;
	const double dtau = Rs * dtauRs;

	if (model == EoSModel::BAG) {
		return PhaseTransition::PTParams_Bag(vw, alpha, Tref, beta, dtau, "exp", un, 1./3., 1./3.);
		// return PhaseTransition::PTParams_Bag(vw, alpha);
	} else if (model == EoSModel::BAGext) {
		return PhaseTransition::PTParams_Bag(vw, alpha, Tref, beta, dtau, "exp", un, cs_true*cs_true, cs_false*cs_false);
	} else {
		if (!eos) {
			throw std::invalid_argument("EquationOfState required for VEFF model but not provided");
		}
		return PhaseTransition::PTParams_Veff(vw, alpha, Tref, beta, dtau, "exp", un, *eos);
	}
}

inline PTParamsVariant
get_pt_params_from_thermal_params(
	const ThermalParams& tps,
	const PhaseTransition::Universe& un,
	const double& vw,
	const double& dtauRs,
	EoSModel model = EoSModel::BAG)
{
	return get_pt_params_from_thermal_params(tps, un, vw, dtauRs, model, nullptr);
}

inline PTParamsVariant
get_pt_params_from_thermal_params(
	const ThermalParams& tps,
	const double& vw,
	const double& dtauRs,
	EoSModel model = EoSModel::BAG,
	const double& dof = 107.75)
{
	const PhaseTransition::Universe un = get_universe_from_thermal_params(tps, dof);
	
	if (model == EoSModel::VEFF) {
		const PhaseTransition::EquationOfState eos = phasetracer_EoS_to_deepphase_EoS(tps.eos);
		return get_pt_params_from_thermal_params(tps, un, vw, dtauRs, model, &eos);
	}
	
	return get_pt_params_from_thermal_params(tps, un, vw, dtauRs, model, nullptr);
}

// helper methods for extracting and checking variant object
inline bool 
is_bag_model(const PTParamsVariant& params) 
{
	return std::holds_alternative<PhaseTransition::PTParams_Bag>(params);
}

inline bool 
is_veff_model(const PTParamsVariant& params) 
{
	return std::holds_alternative<PhaseTransition::PTParams_Veff>(params);
}

inline PhaseTransition::PTParams_Bag& 
get_bag(PTParamsVariant& params) 
{
	return std::get<PhaseTransition::PTParams_Bag>(params);
}

inline const PhaseTransition::PTParams_Bag& 
get_bag(const PTParamsVariant& params) 
{
	return std::get<PhaseTransition::PTParams_Bag>(params);
}

inline PhaseTransition::PTParams_Veff& 
get_veff(PTParamsVariant& params) 
{
	return std::get<PhaseTransition::PTParams_Veff>(params);
}

inline const PhaseTransition::PTParams_Veff& 
get_veff(const PTParamsVariant& params) 
{
	return std::get<PhaseTransition::PTParams_Veff>(params);
}

struct DeepPhaseResults
{
	EoSModel eos_model;

	std::unique_ptr<PhaseTransition::PTParams_Bag> ptparams_bag;
	std::unique_ptr<Hydrodynamics::FluidProfile> profile_bag;
	std::unique_ptr<Spectrum::PowerSpec> spectrum_bag;

	std::unique_ptr<PhaseTransition::PTParams_Bag> ptparams_bag_ext;
	std::unique_ptr<Hydrodynamics::FluidProfile> profile_bag_ext;
	std::unique_ptr<Spectrum::PowerSpec> spectrum_bag_ext;

	std::unique_ptr<PhaseTransition::PTParams_Veff> ptparams_veff;
	std::unique_ptr<Hydrodynamics::FluidProfile> profile_veff;
	std::unique_ptr<Spectrum::PowerSpec> spectrum_veff;

	DeepPhaseResults(const ThermalParams& tps_, const double& vw_, const double& dtauRs_, EoSModel eos_model_)
	: eos_model(eos_model_)
	{
		auto kRs_vals = logspace(1e-3, 1e+3, 100);
		
		if (eos_model == EoSModel::BAG || eos_model == EoSModel::ALL) 
		{
			auto ptparams_variant = get_pt_params_from_thermal_params(tps_, vw_, dtauRs_, EoSModel::BAG);
			ptparams_bag = std::make_unique<PhaseTransition::PTParams_Bag>(get_bag(ptparams_variant));
			profile_bag = std::make_unique<Hydrodynamics::FluidProfile>(*ptparams_bag);
			spectrum_bag = std::make_unique<Spectrum::PowerSpec>(Spectrum::GWSpec2(kRs_vals, *ptparams_bag));
		}

		if (eos_model == EoSModel::BAGext || eos_model == EoSModel::ALL) 
		{
			auto ptparams_variant = get_pt_params_from_thermal_params(tps_, vw_, dtauRs_, EoSModel::BAGext);
			ptparams_bag_ext = std::make_unique<PhaseTransition::PTParams_Bag>(get_bag(ptparams_variant));
			profile_bag_ext = std::make_unique<Hydrodynamics::FluidProfile>(*ptparams_bag_ext);
			spectrum_bag_ext = std::make_unique<Spectrum::PowerSpec>(Spectrum::GWSpec2(kRs_vals, *ptparams_bag_ext));
		}
		
		if (eos_model == EoSModel::VEFF || eos_model == EoSModel::ALL)
		{
			auto ptparams_variant = get_pt_params_from_thermal_params(tps_, vw_, dtauRs_, EoSModel::VEFF);
			ptparams_veff = std::make_unique<PhaseTransition::PTParams_Veff>(get_veff(ptparams_variant));
			profile_veff = std::make_unique<Hydrodynamics::FluidProfile>(*ptparams_veff);
			spectrum_veff = std::make_unique<Spectrum::PowerSpec>(Spectrum::GWSpec2(kRs_vals, *ptparams_veff));
		}
	}
	
	// Check if BAG model results are available
	bool has_bag() const { return ptparams_bag != nullptr; }

	// Check if BAGext model results are available
	bool has_bag_ext() const { return ptparams_bag_ext != nullptr; }

	// Check if VEFF model results are available
	bool has_veff() const { return ptparams_veff != nullptr; }
	
	// Accessor methods for BAG model (throws if not available)
	PhaseTransition::PTParams_Bag& get_ptparams_bag() { 
		if (!ptparams_bag) throw std::runtime_error("BAG model results not available");
		return *ptparams_bag; 
	}
	
	const PhaseTransition::PTParams_Bag& get_ptparams_bag() const { 
		if (!ptparams_bag) throw std::runtime_error("BAG model results not available");
		return *ptparams_bag; 
	}
	
	Hydrodynamics::FluidProfile& get_profile_bag() { 
		if (!profile_bag) throw std::runtime_error("BAG model results not available");
		return *profile_bag; 
	}
	
	const Hydrodynamics::FluidProfile& get_profile_bag() const { 
		if (!profile_bag) throw std::runtime_error("BAG model results not available");
		return *profile_bag; 
	}
	
	Spectrum::PowerSpec& get_spectrum_bag() { 
		if (!spectrum_bag) throw std::runtime_error("BAG model results not available");
		return *spectrum_bag; 
	}
	
	const Spectrum::PowerSpec& get_spectrum_bag() const { 
		if (!spectrum_bag) throw std::runtime_error("BAG model results not available");
		return *spectrum_bag; 
	}

	// Accessor methods for BAGext model (throws if not available)
	PhaseTransition::PTParams_Bag& get_ptparams_bag_ext() { 
		if (!ptparams_bag_ext) throw std::runtime_error("BAG model results not available");
		return *ptparams_bag_ext; 
	}
	
	const PhaseTransition::PTParams_Bag& get_ptparams_bag_ext() const { 
		if (!ptparams_bag_ext) throw std::runtime_error("BAGext model results not available");
		return *ptparams_bag_ext; 
	}
	
	Hydrodynamics::FluidProfile& get_profile_bag_ext() { 
		if (!profile_bag_ext) throw std::runtime_error("BAGext model results not available");
		return *profile_bag_ext; 
	}
	
	const Hydrodynamics::FluidProfile& get_profile_bag_ext() const { 
		if (!profile_bag_ext) throw std::runtime_error("BAGext model results not available");
		return *profile_bag_ext; 
	}
	
	Spectrum::PowerSpec& get_spectrum_bag_ext() { 
		if (!spectrum_bag_ext) throw std::runtime_error("BAGext model results not available");
		return *spectrum_bag_ext; 
	}
	
	const Spectrum::PowerSpec& get_spectrum_bag_ext() const { 
		if (!spectrum_bag_ext) throw std::runtime_error("BAGext model results not available");
		return *spectrum_bag_ext; 
	}
	
	// Accessor methods for VEFF model (throws if not available)
	PhaseTransition::PTParams_Veff& get_ptparams_veff() { 
		if (!ptparams_veff) throw std::runtime_error("VEFF model results not available");
		return *ptparams_veff; 
	}
	
	const PhaseTransition::PTParams_Veff& get_ptparams_veff() const { 
		if (!ptparams_veff) throw std::runtime_error("VEFF model results not available");
		return *ptparams_veff; 
	}
	
	Hydrodynamics::FluidProfile& get_profile_veff() { 
		if (!profile_veff) throw std::runtime_error("VEFF model results not available");
		return *profile_veff; 
	}
	
	const Hydrodynamics::FluidProfile& get_profile_veff() const { 
		if (!profile_veff) throw std::runtime_error("VEFF model results not available");
		return *profile_veff; 
	}
	
	Spectrum::PowerSpec& get_spectrum_veff() { 
		if (!spectrum_veff) throw std::runtime_error("VEFF model results not available");
		return *spectrum_veff; 
	}
	
	const Spectrum::PowerSpec& get_spectrum_veff() const { 
		if (!spectrum_veff) throw std::runtime_error("VEFF model results not available");
		return *spectrum_veff; 
	}
};

} // namespace PhaseTracer::DeepPhaseInterface

#endif // PHASETRACER_DP_HPP_