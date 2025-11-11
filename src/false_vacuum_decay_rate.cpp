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

#include <cmath>
#include "logger.hpp"
#include "false_vacuum_decay_rate.hpp"

namespace PhaseTracer {

    void FalseVacuumDecayRate::get_splines() {

        if (spline_evaluations < 2) {
            throw std::runtime_error("Spline evaluations must be at least 2.");
        }

        std::vector<double> valid_temps, valid_actions, valid_log_gammas;
        double dt = (t_max - t_min) / (spline_evaluations - 1);
        const double log_gamma_min = -700;

        for (int i = 0; i < spline_evaluations; i++) 
        {
            double tt = t_min + i * dt;
            auto action_opt = calculate_action(tt);
            if (action_opt.has_value()) 
            {
                double action = action_opt.value();
                auto log_gamma_opt = calculate_log_gamma(tt, action);
                double log_gamma = log_gamma_min;
                if (log_gamma_opt.has_value()) {
                    log_gamma = std::max(log_gamma_opt.value(), log_gamma_min);
                }
                valid_temps.push_back(tt);
                valid_actions.push_back(action);
                valid_log_gammas.push_back(log_gamma);
            } else {
                LOG(debug) << "t = " << tt << ", action = INVALID" << std::endl;
            }
        }

        if (valid_temps.size() < 2) 
        {
            throw std::runtime_error("Not enough valid action points to build spline.");
        }

        alglib::real_1d_array temp_array, action_array, log_gamma_array;
        temp_array.setcontent(valid_temps.size(), valid_temps.data());
        action_array.setcontent(valid_actions.size(), valid_actions.data());
        log_gamma_array.setcontent(valid_log_gammas.size(), valid_log_gammas.data());

        LOG(debug) << "Building action and gamma arrays";
        alglib::spline1dbuildcubic(temp_array, action_array, this->action_spline);
        alglib::spline1dbuildcubic(temp_array, log_gamma_array, this->log_gamma_spline);
        LOG(debug) << "Action and log(gamma) arrays built.";
    }

    std::optional<double> 
    FalseVacuumDecayRate::calculate_action(double temperature) const
    {
        double action = tf.get_action(t.true_phase, t.false_phase, temperature) / temperature;
        if (std::isnan(action) || std::isinf(action) || action > 1e150) {
        return std::nullopt;
        }
        return action;
    }

    std::optional<double>
    FalseVacuumDecayRate::calculate_log_gamma(double temperature, double action) const
    {
        double prefactor = decay_rate_prefactor(temperature, action);
        double log_gamma = log(prefactor) - action;
        return log_gamma < -700 ? std::nullopt : std::optional<double>(log_gamma);
    }

    double
    FalseVacuumDecayRate::decay_rate_prefactor(double temperature, double action_on_T) const
    {
        return prefactor_function(temperature, action_on_T);
    }

    std::function<double(double, double)>
    FalseVacuumDecayRate::default_decay_rate_prefactor()
    {
        return [](double temperature, double action_on_T) -> double {
            double t4 = temperature * temperature * temperature * temperature;
            double ratio = std::pow(action_on_T / (2. * M_PI), 1.5);
            return t4 * ratio;
        };
    }

    double
    FalseVacuumDecayRate::get_action(const double& temperature) const
    {
        double action_on_T = alglib::spline1dcalc(action_spline, temperature);
        return action_on_T * temperature;
    }

    double 
    FalseVacuumDecayRate::get_action_deriv(const double& temperature) const
    {
        double y, dy, ddy;
        alglib::spline1ddiff(action_spline, temperature, y, dy, ddy);
        return dy;
    }

    double
    FalseVacuumDecayRate::get_gamma(const double& temperature) const
    {
        double log_gamma = alglib::spline1dcalc(log_gamma_spline, temperature);
        return exp(log_gamma);
    }

} // namespace PhaseTracer