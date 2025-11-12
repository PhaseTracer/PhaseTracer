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
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "logger.hpp"
#include "false_vacuum_decay_rate.hpp"

namespace PhaseTracer {

    void FalseVacuumDecayRate::get_splines() {

        if (spline_evaluations < 2) {
            throw std::runtime_error("Spline evaluations must be at least 2.");
        }

        const double log_gamma_min = -700;
        double dt = (t_max - t_min) / (spline_evaluations - 1);

        std::vector<double> temp_results(spline_evaluations);
        std::vector<double> action_results(spline_evaluations);
        std::vector<double> log_gamma_results(spline_evaluations);
        std::vector<bool> valid_flags(spline_evaluations, false);
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Note: OpenMP parallelization disabled because ActionCalculator has mutable state
        // and is not thread-safe (contains mutable bubble_profile, phi_for_profile, tunneling_path)
        // #ifdef _OPENMP
        // #pragma omp parallel for schedule(dynamic)
        // #endif

        for (int i = 0; i < spline_evaluations; i++) 
        {
            double tt = t_min + i * dt;
            
            double action;
            try {
                action = ac.get_action(t.true_phase, t.false_phase, tt) / tt;
            } catch ( ... )
            {
                action = 1.;
                valid_flags[i] = false;
                continue;
            }
            
            if (std::isnan(action) || std::isinf(action) || action > 1e150) 
            {
                valid_flags[i] = false;
                continue;
            }
            
            double prefactor = decay_rate_prefactor(tt, action);
            double log_gamma = log(prefactor) - action;
            
            if (log_gamma < -700) 
            {
                log_gamma = log_gamma_min;
            } else {
                log_gamma = std::max(log_gamma, log_gamma_min);
            }
            
            temp_results[i] = tt;
            action_results[i] = action;
            log_gamma_results[i] = log_gamma;
            valid_flags[i] = true;
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        LOG(info) << "Action calculation loop completed in " << duration.count() << " ms" << std::endl;
        
        std::vector<double> valid_temps, valid_actions, valid_log_gammas;
        for (int i = 0; i < spline_evaluations; i++) 
        {
            if (valid_flags[i]) 
            {
                valid_temps.push_back(temp_results[i]);
                valid_actions.push_back(action_results[i]);
                valid_log_gammas.push_back(log_gamma_results[i]);
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

        alglib::spline1dbuildcubic(temp_array, action_array, this->action_spline);
        alglib::spline1dbuildcubic(temp_array, log_gamma_array, this->log_gamma_spline);
        
        LOG(debug) << "Built splines with " << valid_temps.size() << " valid points" << std::endl;
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