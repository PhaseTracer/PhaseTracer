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
#include "thermo_finder.hpp"
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
	MUNU,
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

inline double 
refine_lower_bound(const PhaseTracer::EquationOfState& eos, double step_size = 0.001, size_t max_iterations = 100)
{
    double t = eos.get_t_min();
    const double t_max = eos.get_t_max();
    
    if (t >= t_max) {
        throw std::invalid_argument("refine_lower_bound: t_min >= t_max");
    }

    auto [p_plus, p_minus] = eos.get_pressure(t);
    auto [e_plus, e_minus] = eos.get_energy(t);

    size_t iterations = 0;
    
    while ((p_plus < 0 || p_minus < 0 || e_plus < 0 || e_minus < 0) && t < t_max && iterations < max_iterations)
    {
        t += step_size;
        
        auto [p_plus_new, p_minus_new] = eos.get_pressure(t);
        auto [e_plus_new, e_minus_new] = eos.get_energy(t);
        
        p_plus = p_plus_new;
        p_minus = p_minus_new;
        e_plus = e_plus_new;
        e_minus = e_minus_new;
        
        ++iterations;
    }
    
    if (iterations >= max_iterations) {
        throw std::runtime_error("refine_lower_bound: Maximum iterations reached without finding valid thermodynamic region");
    }
    
    if (t >= t_max) {
        throw std::runtime_error("refine_lower_bound: No valid thermodynamic region found in [t_min, t_max]");
    }

    return t;
}

inline PhaseTransition::EquationOfState
phasetracer_EoS_to_deepphase_EoS(const PhaseTracer::EquationOfState& eos)
{
	double t_min = eos.get_t_min();
	double t_max = eos.get_t_max();
	int n_temp = eos.get_n_temp();
	double dt = (t_max - t_min) / (n_temp - 1);

	std::vector<double> t_vals, p_plus_vals, p_minus_vals, e_plus_vals, e_minus_vals;

	for (int i = 0; i < n_temp; ++i) 
	{
		double T = t_min + i * dt;
		auto [p_plus, p_minus] = eos.get_pressure(T);
		auto [e_plus, e_minus] = eos.get_energy(T);

		if(p_plus < 0 || p_minus < 0 || e_plus < 0 || e_minus < 0) { continue; }

		t_vals.push_back(T);
		p_plus_vals.push_back(p_plus);
		p_minus_vals.push_back(p_minus);
		e_plus_vals.push_back(e_plus);
		e_minus_vals.push_back(e_minus);
	}

	PhaseTransition::EquationOfState output(t_vals, p_plus_vals, p_minus_vals, e_plus_vals, e_minus_vals);
	return output;
}

inline PhaseTransition::Universe
get_universe_from_transition_milestone(const TransitionMilestone& milestone, const double& dof = 107.75) 
{
	if (milestone.status == MilestoneStatus::YES) 
	{
		if (milestone.temperature <= 0 || milestone.H <= 0) 
		{
			throw std::invalid_argument("Invalid thermal parameters for universe creation");
		}
		return PhaseTransition::Universe(milestone.temperature, dof, milestone.H);
	}

	PhaseTransition::Universe un;
	return un;
}

// general constructor
inline PTParamsVariant
get_pt_params_from_transition_milestone(
	const TransitionMilestone& milestone,
	const PhaseTransition::Universe& un,
	const double& vw,
	const double& dtauRs,
	EoSModel model,
	const PhaseTransition::EquationOfState* eos_ptr = nullptr)
{

	if (milestone.status != MilestoneStatus::YES) {
		throw std::invalid_argument("No valid thermal parameters: neither nucleates nor percolates");
	}

	const double alpha = milestone.alpha;
	const double betaH = milestone.betaH;
	const double RsH = milestone.Rs;
	const double H = milestone.H;
	const double Tref = milestone.temperature;
	const double cs_plus = milestone.cs_plus;
	const double cs_minus = milestone.cs_minus;

	char* nuc_type; // why isn't this a string in deepphase?
	if (milestone.type == MilestoneType::PERCOLATION)
	{
		if (milestone.nucleation_type == PhaseTracer::NucleationType::EXPONENTIAL)
		{
			nuc_type = "exp";
		} else 
		{
			nuc_type = "sim";
		}
	} else {
		nuc_type = "exp";
	}

	const double beta = betaH * H;
	const double Rs = RsH / H;
	const double dtau = Rs * dtauRs;

	if (alpha <= 0 || Rs <= 0) {
		throw std::invalid_argument("Invalid thermal parameters: alpha or betaH <= 0");
	}

	if (model == EoSModel::BAG) {
		return PhaseTransition::PTParams_Bag(vw, alpha, Tref, beta, Rs, dtau, nuc_type, un, 1./3., 1./3.);
	} else if (model == EoSModel::MUNU) {
		return PhaseTransition::PTParams_Bag(vw, alpha, Tref, beta, Rs, dtau, nuc_type, un, cs_plus*cs_plus, cs_minus*cs_minus);
	} else {
		if (!eos_ptr) {
			throw std::invalid_argument("EquationOfState required for VEFF model but not provided");
		}
		return PhaseTransition::PTParams_Veff(vw, alpha, Tref, beta, Rs, dtau, nuc_type, un, *eos_ptr);
	}
}

// takes universe as input
inline PTParamsVariant
get_pt_params_from_transition_milestone(
	const TransitionMilestone& milestone,
	const EquationOfState& eos,
	const PhaseTransition::Universe& un,
	const double& vw,
	const double& dtauRs,
	EoSModel model = EoSModel::BAG)
{
	return get_pt_params_from_transition_milestone(milestone, un, vw, dtauRs, model, nullptr);
}

// takes dof and constructs the universe object
inline PTParamsVariant
get_pt_params_from_transition_milestone(
	const TransitionMilestone& milestone,
	const EquationOfState& eos,
	const double& dof,
	const double& vw,
	const double& dtauRs,
	EoSModel model = EoSModel::BAG
	)
{
	const PhaseTransition::Universe un = get_universe_from_transition_milestone(milestone, dof);
	
	if (model == EoSModel::VEFF) {
		const PhaseTransition::EquationOfState eos_dp = phasetracer_EoS_to_deepphase_EoS(eos);
		return get_pt_params_from_transition_milestone(milestone, un, vw, dtauRs, model, &eos_dp);
	}
	
	return get_pt_params_from_transition_milestone(milestone, un, vw, dtauRs, model, nullptr);
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
	TransitionMilestone milestone;
	const EquationOfState& eos;
	double dof;
	double vw;
	double dtauRs;

	std::unique_ptr<PhaseTransition::PTParams_Bag> ptparams_bag;
	std::unique_ptr<Hydrodynamics::FluidProfile> profile_bag;
	std::unique_ptr<Spectrum::PowerSpec> spectrum_bag;

	std::unique_ptr<PhaseTransition::PTParams_Bag> ptparams_munu;
	std::unique_ptr<Hydrodynamics::FluidProfile> profile_munu;
	std::unique_ptr<Spectrum::PowerSpec> spectrum_munu;

	std::unique_ptr<PhaseTransition::PTParams_Veff> ptparams_veff;
	std::unique_ptr<Hydrodynamics::FluidProfile> profile_veff;
	std::unique_ptr<Spectrum::PowerSpec> spectrum_veff;

	DeepPhaseResults(const TransitionMilestone& milestone_, const EquationOfState& eos_, const double& dof_, const double& vw_, const double& dtauRs_, EoSModel eos_model_)
	: eos_model(eos_model_), milestone(milestone_), eos(eos_), dof(dof_), vw(vw_), dtauRs(dtauRs_)
	{
		auto kRs_vals = logspace(1e-2, 5e+3, 100);
		
		if (eos_model == EoSModel::BAG || eos_model == EoSModel::ALL) 
		{
			auto ptparams_variant = get_pt_params_from_transition_milestone(milestone, eos, dof, vw, dtauRs, EoSModel::BAG);
			ptparams_bag = std::make_unique<PhaseTransition::PTParams_Bag>(get_bag(ptparams_variant));
			profile_bag = std::make_unique<Hydrodynamics::FluidProfile>(*ptparams_bag);
			std::cout << profile_bag->xi_min() << "\n";
			std::cout << profile_bag->xi_max() << "\n";
			std::cout << profile_bag->mode_str() << "\n";
			spectrum_bag = std::make_unique<Spectrum::PowerSpec>(Spectrum::GWSpec(kRs_vals, *ptparams_bag));
		}

		if (eos_model == EoSModel::MUNU || eos_model == EoSModel::ALL) 
		{
			auto ptparams_variant = get_pt_params_from_transition_milestone(milestone, eos, dof, vw, dtauRs, EoSModel::MUNU);
			ptparams_munu = std::make_unique<PhaseTransition::PTParams_Bag>(get_bag(ptparams_variant));
			profile_munu = std::make_unique<Hydrodynamics::FluidProfile>(*ptparams_munu);
			spectrum_munu = std::make_unique<Spectrum::PowerSpec>(Spectrum::GWSpec(kRs_vals, *ptparams_munu));
		}
		
		if (eos_model == EoSModel::VEFF || eos_model == EoSModel::ALL)
		{
			auto ptparams_variant = get_pt_params_from_transition_milestone(milestone, eos, dof, vw, dtauRs, EoSModel::VEFF);
			ptparams_veff = std::make_unique<PhaseTransition::PTParams_Veff>(get_veff(ptparams_variant));
			profile_veff = std::make_unique<Hydrodynamics::FluidProfile>(*ptparams_veff);
			std::cout << profile_veff->xi_min() << "\n";
			std::cout << profile_veff->xi_max() << "\n";
			spectrum_veff = std::make_unique<Spectrum::PowerSpec>(Spectrum::GWSpec(kRs_vals, *ptparams_veff));
		}
	}
	
	// Check if BAG model results are available
	bool has_bag() const { return ptparams_bag != nullptr; }

	// Check if BAGext model results are available
	bool has_munu() const { return ptparams_munu != nullptr; }

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
	PhaseTransition::PTParams_Bag& get_ptparams_munu() { 
		if (!ptparams_munu) throw std::runtime_error("BAG model results not available");
		return *ptparams_munu; 
	}
	
	const PhaseTransition::PTParams_Bag& get_ptparams_munu() const { 
		if (!ptparams_munu) throw std::runtime_error("MUNU model results not available");
		return *ptparams_munu; 
	}
	
	Hydrodynamics::FluidProfile& get_profile_munu() { 
		if (!profile_munu) throw std::runtime_error("MUNU model results not available");
		return *profile_munu; 
	}
	
	const Hydrodynamics::FluidProfile& get_profile_munu() const { 
		if (!profile_munu) throw std::runtime_error("MUNU model results not available");
		return *profile_munu; 
	}
	
	Spectrum::PowerSpec& get_spectrum_munu() { 
		if (!spectrum_munu) throw std::runtime_error("MUNU model results not available");
		return *spectrum_munu; 
	}
	
	const Spectrum::PowerSpec& get_spectrum_munu() const { 
		if (!spectrum_munu) throw std::runtime_error("MUNU model results not available");
		return *spectrum_munu; 
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