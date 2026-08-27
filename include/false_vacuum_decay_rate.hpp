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
#include "action_calculator.hpp"

namespace PhaseTracer {

class FalseVacuumDecayRate {

public:

    /**
     * Signature of a decay-rate prefactor A(T). It receives the temperature,
     * the action-on-temperature S/T, and the full ActionResult (bounce action,
     * profile and tunneling path) so that a prefactor needing the bounce
     * solution -- e.g. a one-loop functional determinant via BubbleDet -- can
     * access it. The default analytic prefactor ignores the ActionResult.
     */
    using PrefactorFunction = std::function<double(double temperature, double action_on_T, const ActionResult& bounce)>;

    
    // Delete copy constructor and copy assignment to prevent shallow copies of ALGLIB splines
    FalseVacuumDecayRate(const FalseVacuumDecayRate&) = delete;
    FalseVacuumDecayRate& operator=(const FalseVacuumDecayRate&) = delete;
    
    // Allow move semantics
    FalseVacuumDecayRate(FalseVacuumDecayRate&&) = default;
    FalseVacuumDecayRate& operator=(FalseVacuumDecayRate&&) = default;
    
    FalseVacuumDecayRate(Transition t_in, const ActionCalculator& ac_in)
    : ac(ac_in), t(t_in), t_min(t_in.false_phase.T.front()), t_max(t_in.TC) {}

    /**
     * @brief Solves the bounce action over [t_min, t_max] and fits the
     *        log(action), log(prefactor) and log(gamma) splines.
     *
     * This is the expensive part of the class. It must be called before any of
     * the get_action/get_prefactor/get_gamma accessors, which will otherwise 
     * throw a logic error.
     */
    void calculate();

    /** @brief Whether calculate() has completed successfully. */
    bool is_calculated() const { return calculated; }

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
     * @brief Computes d^2(S/T)/dT^2 at a given temperature using the precomputed spline.
     * @param temperature The temperature at which to evaluate the action.
     * @return The second derivative of the action at the specified temperature.
    */
    double get_action_double_deriv(const double& temperature) const;

    /** 
     * @brief Computes the false vacuum decay rate at a given temperature using the precomputed spline.
     * @param temperature The temperature at which to evaluate the action.
     * @return The false vacuum decay rate at the specified temperature.
    */
    double get_gamma(const double& temperature) const;

    /** 
     * @brief Evaluates the decay rate prefactor A(T) using the precomputed spline.
     * @param temperature The temperature at which to evaluate the prefactor.
     * @return The prefactor A(T).
    */
    double get_prefactor(const double& temperature) const;

    /** 
     * @brief Set a custom decay rate prefactor function.
     * @param custom_prefactor A function taking (temperature, action_on_T) and returning the prefactor.
    */
    void set_prefactor_function(PrefactorFunction custom_prefactor) {
        prefactor_function = custom_prefactor;
    }

    /**
     * @brief Get the current decay rate prefactor function.
     * @return The prefactor function.
    */
    const PrefactorFunction& get_prefactor_function() const {
        return prefactor_function;
    }

    /**
     * @brief Compute the decay rate prefactor using the current prefactor function.
     * @param temperature The temperature at which to evaluate.
     * @param action_on_T The action divided by temperature (S/T).
     * @param bounce The full bounce solution (action, profile, path) at this temperature.
     * @return The decay rate prefactor.
    */
    double decay_rate_prefactor(double temperature, double action_on_T, const ActionResult& bounce) const;

    /**
     * @brief Default decay rate prefactor function following standard bounce action formula.
     * @return A function object that computes the standard prefactor.
    */
    static PrefactorFunction default_decay_rate_prefactor();

    /** 
     * @brief Computes the bubble profile at a given temperature.
     *        Triggers an action evaluation at the specified temperature,
     *        which populates the bubble profile cache in ActionCalculator.
     * @param temperature The temperature at which to compute the bubble profile.
     * @return The bubble profile at the specified temperature.
    */
    Profile1D get_bubble_profile(const double& temperature) {
        if (temperature < t_min || temperature > t_max) {
            throw std::out_of_range("Temperature is outside the valid range [t_min, t_max].");
        }
        return ac.get_action_full(t.true_phase, t.false_phase, temperature).bubble_profile;
    }

    /**
     * @brief Write the decay rate data to a file.
     * @param filename The name of the file to write to.
     * @param n_steps The number of steps to use in the output.
    */
    void write(const std::string& filename, const int& n_steps=100);

private:

    /** Compute splines for action and log(gamma) */
    void get_splines();

    /** Throws std::logic_error if calculate() has not been called */
    void require_calculated(const char* caller) const;

    /** Reference to ActionCalculator class */
    const ActionCalculator& ac;

    /** Transition for which the decay rate is computed */
    Transition t;

    /** Splines for log(action), log(prefactor), and log(gamma) */
    alglib::spline1dinterpolant log_action_spline, log_prefactor_spline, log_gamma_spline;

    /** Function for computing the decay rate prefactor */
    PrefactorFunction prefactor_function = default_decay_rate_prefactor();

    /** Set by calculate() once the splines are built */
    bool calculated = false;

    /** Minimum and maximum temperatures for which the decay rate is computed. */
    PROPERTY(double, t_min, 0.0)

    PROPERTY(double, t_max, 0.0)

    /** Number of spline evaluations for building the splines */
    PROPERTY(int, spline_evaluations, 50)

}; // class FalseVacuumDecayRate

} // namespace PhaseTracer

#endif // PHASETRACER_FALSE_VACUUM_DECAY_RATE_HPP_