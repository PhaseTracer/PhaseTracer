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

    void 
    FalseVacuumDecayRate::calculate() 
    {
        calculated = false;
        get_splines();
        calculated = true;
    }

    void 
    FalseVacuumDecayRate::require_calculated(const char* caller) const 
    {
        if (!calculated) 
        {
            throw std::logic_error(std::string("FalseVacuumDecayRate::") + caller + " called before calculate().");
        }
    }

    void 
    FalseVacuumDecayRate::get_splines() 
    {
        if (spline_evaluations < 2) {
            throw std::runtime_error("Spline evaluations must be at least 2.");
        }

        const double log_gamma_min = -700;
        double dt = (t_max - t_min) / (spline_evaluations - 1);

        std::vector<double> temp_results(spline_evaluations);
        std::vector<double> log_action_results(spline_evaluations);
        std::vector<double> log_prefactor_results(spline_evaluations);
        std::vector<double> log_gamma_results(spline_evaluations);
        std::vector<char> valid_flags(spline_evaluations, 0);

        auto start_time = std::chrono::high_resolution_clock::now();

        // The bounce is solved on contiguous blocks of temperatures rather than
        // one temperature at a time, so that the converged tunneling path at
        // one temperature can seed the next. The path varies smoothly with T,
        // so this removes most of the deformation iterations that would
        // otherwise be spent rediscovering it from a straight line.
        //
        // Blocks are kept short and scheduled dynamically: the work per
        // temperature is uneven, and on hybrid CPUs the cores are too, so long
        // static blocks would leave fast cores idle. A chunk of 1 disables
        // warm starting entirely.
        int n_threads = 1;
        #ifdef _OPENMP
        n_threads = omp_get_max_threads();
        #endif
        int chunk = warm_start_chunk;
        if (chunk <= 0)
        {
            chunk = std::max(1, static_cast<int>(
                std::ceil(static_cast<double>(spline_evaluations) / n_threads)));
        }
        const int n_chunks = (spline_evaluations + chunk - 1) / chunk;

        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (int c = 0; c < n_chunks; c++)
        {
            const int begin = c * chunk;
            const int end = std::min(spline_evaluations, begin + chunk);
            
            std::vector<Eigen::VectorXd> path_guess;
            for (int i = end - 1; i >= begin; i--)
            {
                double tt = t_min + i * dt;

                ActionResult bounce;
                double action;
                try {
                    bounce = ac.get_action_full(t.true_phase, t.false_phase, tt, 0, path_guess);
                    action = bounce.action / tt; // NB action is S_3(T)/T
                } catch (const std::exception& e)
                {
                    LOG(warning) << "Action evaluation failed at T = " << tt << ": " << e.what();
                    path_guess.clear();
                    continue;
                } catch (...)
                {
                    LOG(warning) << "Action evaluation failed at T = " << tt << ": unknown error";
                    path_guess.clear();
                    continue;
                }

                if (std::isnan(action) || std::isinf(action) || action > 1e150 || action < 0)
                {
                    LOG(debug) << "Rejected non-physical action " << action
                               << " at T = " << tt;
                    path_guess.clear();
                    continue;
                }

                // Only carry a path forward once it has produced a usable action.
                path_guess = bounce.tunneling_path;

                double prefactor = decay_rate_prefactor(tt, action, bounce);
                if (!std::isfinite(prefactor) || prefactor <= 0.0)
                {
                    LOG(debug) << "Rejected non-physical prefactor " << prefactor
                               << " at T = " << tt;
                    continue;
                }
                double log_prefactor = std::log(prefactor);
                double log_gamma = std::max(log_prefactor - action, log_gamma_min);

                temp_results[i] = tt;
                log_action_results[i] = std::log(action);
                log_prefactor_results[i] = log_prefactor;
                log_gamma_results[i] = log_gamma;
                valid_flags[i] = 1;
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        LOG(debug) << "Action calculation loop completed in " << duration.count() << " ms";
        
        std::vector<double> valid_temps, valid_log_actions, valid_log_prefactors, valid_log_gammas;
        for (int i = 0; i < spline_evaluations; i++) 
        {
            if (valid_flags[i]) 
            {
                valid_temps.push_back(temp_results[i]);
                valid_log_actions.push_back(log_action_results[i]);
                valid_log_prefactors.push_back(log_prefactor_results[i]);
                valid_log_gammas.push_back(log_gamma_results[i]);
            }
        }
        
        if (valid_temps.size() < 2) 
        {
            throw std::runtime_error("Not enough valid action points to build spline.");
        }

        alglib::real_1d_array temp_array, log_action_array, log_prefactor_array, log_gamma_array;
        temp_array.setcontent(valid_temps.size(), valid_temps.data());
        log_action_array.setcontent(valid_log_actions.size(), valid_log_actions.data());
        log_prefactor_array.setcontent(valid_log_prefactors.size(), valid_log_prefactors.data());
        log_gamma_array.setcontent(valid_log_gammas.size(), valid_log_gammas.data());

        alglib::spline1dbuildcubic(temp_array, log_action_array, this->log_action_spline);
        alglib::spline1dbuildcubic(temp_array, log_prefactor_array, this->log_prefactor_spline);
        alglib::spline1dbuildcubic(temp_array, log_gamma_array, this->log_gamma_spline);
        
        LOG(debug) << "Built action splines with " << valid_temps.size() << " valid points";
    }

    double
    FalseVacuumDecayRate::decay_rate_prefactor(double temperature, double action_on_T, const ActionResult& bounce) const
    {
        return prefactor_function(temperature, action_on_T, bounce);
    }

    FalseVacuumDecayRate::PrefactorFunction
    FalseVacuumDecayRate::default_decay_rate_prefactor()
    {
        return [](double temperature, double action_on_T, const ActionResult&) -> double {
            double t4 = temperature * temperature * temperature * temperature;
            double ratio = std::pow(action_on_T / (2. * M_PI), 1.5);
            return t4 * ratio;
        };
    }

    double
    FalseVacuumDecayRate::get_action(const double& temperature) const
    {
        require_calculated("get_action");
        double log_action_on_T = alglib::spline1dcalc(log_action_spline, temperature);
        return exp(log_action_on_T) * temperature;
    }

    double 
    FalseVacuumDecayRate::get_action_deriv(const double& temperature) const
    {
        require_calculated("get_action_deriv");
        double y, dy, ddy;
        alglib::spline1ddiff(log_action_spline, temperature, y, dy, ddy);
        return dy*exp(y);
    }

    double 
    FalseVacuumDecayRate::get_action_double_deriv(const double& temperature) const
    {
        require_calculated("get_action_double_deriv");
        double y, dy, ddy;
        alglib::spline1ddiff(log_action_spline, temperature, y, dy, ddy);
        return (dy*dy + ddy) * exp(y);
    }

    double
    FalseVacuumDecayRate::get_gamma(const double& temperature) const
    {
        require_calculated("get_gamma");
        double log_gamma = alglib::spline1dcalc(log_gamma_spline, temperature);
        return exp(log_gamma);
    }

    double
    FalseVacuumDecayRate::get_prefactor(const double& temperature) const
    {
        require_calculated("get_prefactor");
        double log_prefactor = alglib::spline1dcalc(log_prefactor_spline, temperature);
        return exp(log_prefactor);
    }

    void
    FalseVacuumDecayRate::write(const std::string& filename, const int& n_steps)
    {
        require_calculated("write");
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filename);
        }

        file << "# Temperature,Action,Prefactor,Gamma\n";
        double dt = (t_max - t_min) / (n_steps - 1);
        for (int i = 0; i < n_steps; ++i) {
            double temperature = t_min + i * dt;
            double action = get_action(temperature);
            double prefactor = get_prefactor(temperature);
            double gamma = get_gamma(temperature);
            file << temperature << "," << action << "," << prefactor << "," << gamma << "\n";
        }
        file.close();
    }

} // namespace PhaseTracer