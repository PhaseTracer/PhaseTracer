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

#ifndef PHASETRACER_FALSE_VACUUM_DECAY_RATE_HPP_
#define PHASETRACER_FALSE_VACUUM_DECAY_RATE_HPP_

#include <cmath>
#include <vector>
#include <optional>
#include <stdexcept>
#include <functional>
#include <interpolation.h>

#include "property.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"

namespace PhaseTracer {

class FalseVacuumDecayRate {

private:

public:
    
    FalseVacuumDecayRate(Transition t_in, ActionCalculator ac_in)
    : t(t_in), ac(ac_in), t_min(t_in.false_phase.T.front()), t_max(t_in.TC), spline_evaluations(50),
      prefactor_function(default_decay_rate_prefactor())
    {
        get_splines();
    }

    FalseVacuumDecayRate(Transition t_in, ActionCalculator ac_in, double t_min_in, double t_max_in, int spline_evaluations_in)
    : t(t_in), ac(ac_in), t_min(t_min_in), t_max(t_max_in), spline_evaluations(spline_evaluations_in),
      prefactor_function(default_decay_rate_prefactor())
    {
        get_splines();
    }

    FalseVacuumDecayRate(Transition t_in, ActionCalculator ac_in, 
                         std::function<double(double, double)> custom_prefactor)
    : t(t_in), ac(ac_in), t_min(t_in.false_phase.T.front()), t_max(t_in.TC), spline_evaluations(50),
      prefactor_function(custom_prefactor)
    {
        get_splines();
    }

    FalseVacuumDecayRate(Transition t_in, ActionCalculator ac_in, double t_min_in, double t_max_in, 
                         int spline_evaluations_in, std::function<double(double, double)> custom_prefactor)
    : t(t_in), ac(ac_in), t_min(t_min_in), t_max(t_max_in), spline_evaluations(spline_evaluations_in),
      prefactor_function(custom_prefactor)
    {
        get_splines();
    }

    /**
     * @brief Retrieves the minimum temperature for which the decay rate is computed.
     * @return The minimum temperature.
     */
    const double get_t_min() const {return t_min;}

    /**
     * @brief Retrieves the maximum temperature for which the decay rate is computed.
     * @return The maximum temperature.
     */
    const double get_t_max() const {return t_max;}

    /** 
     * @brief Computes the action at a given temperature using the precomputed spline.
     * @param temperature The temperature at which to evaluate the action.
     * @return The action at the specified temperature.
    */
    double get_action(const double& temperature) const;

    /** 
     * @brief Computes d(S/T)/dT at a given temperature using the precomputed spline.
     * @param temperature The temperature at which to evaluate the action.
     * @return The derivative of the action at the specified temperature.
    */
    double get_action_deriv(const double& temperature) const;

    /** 
     * @brief Computes the false vacuum decay rate at a given temperature using the precomputed spline.
     * @param temperature The temperature at which to evaluate the action.
     * @return The false vacuum decay rate at the specified temperature.
    */
    double get_gamma(const double& temperature) const;

    /** 
     * @brief Set a custom decay rate prefactor function.
     * @param custom_prefactor A function taking (temperature, action_on_T) and returning the prefactor.
    */
    void set_prefactor_function(std::function<double(double, double)> custom_prefactor) {
        prefactor_function = custom_prefactor;
    }

    /** 
     * @brief Get the current decay rate prefactor function.
     * @return The prefactor function.
    */
    const std::function<double(double, double)>& get_prefactor_function() const {
        return prefactor_function;
    }

    /** 
     * @brief Compute the decay rate prefactor using the current prefactor function.
     * @param temperature The temperature at which to evaluate.
     * @param action_on_T The action divided by temperature (S/T).
     * @return The decay rate prefactor.
    */
    double decay_rate_prefactor(double temperature, double action_on_T) const;

    /** 
     * @brief Default decay rate prefactor function following standard bounce action formula.
     * @return A function object that computes the standard prefactor.
    */
    static std::function<double(double, double)> default_decay_rate_prefactor();

private:

    /** Compute splines for action and log(gamma) */
    void get_splines();

    /** Local instance of ActionCalculator class */
    ActionCalculator ac;

    /** Transition for which the decay rate is computed */
    Transition t;

    /** Number of spline evaluations for building the splines */
    int spline_evaluations;

    /** Minimum and maximum temperatures for which the decay rate is computed */
    double t_min, t_max;

    /** Splines for action and log(gamma) */
    alglib::spline1dinterpolant action_spline, log_gamma_spline;
    
    /** Function for computing the decay rate prefactor */
    std::function<double(double, double)> prefactor_function;

}; // class FalseVacuumDecayRate

} // namespace PhaseTracer

#endif // PHASETRACER_FALSE_VACUUM_DECAY_RATE_HPP_