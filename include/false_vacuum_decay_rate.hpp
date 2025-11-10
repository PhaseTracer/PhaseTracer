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
#include <interpolation.h>

#include "property.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"

namespace PhaseTracer {

class FalseVacuumDecayRate {

private:
    TransitionFinder tf;
    Transition t;
    int spline_evaluations;

public:
    double minimum_temp, maximum_temp;
    alglib::spline1dinterpolant action_spline, log_gamma_spline;

    FalseVacuumDecayRate(Transition t_in, TransitionFinder tf_in)
    : t(t_in), tf(tf_in), minimum_temp(t_in.false_phase.T.front()), maximum_temp(t_in.TC), spline_evaluations(50) 
    {
        get_splines();
    }

    FalseVacuumDecayRate(Transition t_in, TransitionFinder tf_in, double minimum_temp_in, double maximum_temp_in, int spline_evaluations_in)
    : t(t_in), tf(tf_in), minimum_temp(minimum_temp_in), maximum_temp(maximum_temp_in), spline_evaluations(spline_evaluations_in)
    {
        get_splines();
    }

    /**
     * @brief Computes the bounce and log(gamma) splines for given transition.
     * @return void.
     */
    void get_splines();

    const double get_t_min() const {return minimum_temp;}

    const double get_t_max() const {return maximum_temp;}

private:
    /** Calculate action at given temperature */
    std::optional<double> 
    calculate_action(double temperature) const 
    {
        double action = tf.get_action(t.true_phase, t.false_phase, temperature) / temperature;
        if (std::isnan(action) || std::isinf(action) || action > 1e150) {
        return std::nullopt;
        }
        return action;
    }

    /** Calculate log(gamma) at given temperature and action */
    std::optional<double> 
    calculate_log_gamma(double temperature, double action) const 
    {
        double log_gamma = 4*log(temperature) + 1.5*log(action/(2.*M_PI)) - action;
        return log_gamma < -700 ? std::nullopt : std::optional<double>(log_gamma);
    }
};

} // namespace PhaseTracer

#endif // PHASETRACER_FALSE_VACUUM_DECAY_RATE_HPP_