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
#include "thermo_finder.hpp"

namespace PhaseTracer {

    ThermalParameterSet 
    ThermoFinder::get_thermal_parameter_set(Transition t) 
    {
        ThermalParameterSet output(
            t, 
            ac,
            n_temp_action,
            n_temp_eos,
            vw,
            background_dof,
            dof,
            use_pf_in_nt_integrand,
            use_bag_dtdT,
            percolation_target,
            completion_target,
            onset_target,
            nucleation_target,
            temperature_abs_tol
        );

        output.onset = output.transition_metrics.onset_milestone;
        output.onset.set_print_setting(onset_print_setting);
        add_thermal_parameter_values(output.onset, output.decay_rate, output.eos, output.transition_metrics);

        output.percolation = output.transition_metrics.percolation_milestone;
        output.percolation.set_print_setting(percolation_print_setting);
        add_thermal_parameter_values(output.percolation, output.decay_rate, output.eos, output.transition_metrics);

        output.completion = output.transition_metrics.completion_milestone;
        output.completion.set_print_setting(completion_print_setting);
        add_thermal_parameter_values(output.completion, output.decay_rate, output.eos, output.transition_metrics);

        output.nucleation = output.transition_metrics.nucleation_milestone;
        output.nucleation.set_print_setting(nucleation_print_setting);
        add_thermal_parameter_values(output.nucleation, output.decay_rate, output.eos, output.transition_metrics);

        if(compute_profiles)
        {
            ThermalProfiles profile_out;
            double t_min = output.decay_rate.get_t_min();
            double t_max = output.decay_rate.get_t_max();
            double dt = (t_max - t_min)/(n_temp_profiles-1);

            for(double tt = t_min; tt < t_max; tt += dt)
            {
                double dtdT, dt, H, action, gamma, vext, pf, nt, n, Rs, Rbar;

                try {
                    dtdT = output.transition_metrics.get_dtdT(tt);
                    dt = get_dt(tt, output.transition_metrics);
                    H = get_H(tt, output.transition_metrics);
                    action = output.decay_rate.get_action(tt)/tt;
                    gamma = output.decay_rate.get_gamma(tt);
                    vext = output.transition_metrics.get_extended_volume_from_spline(tt);
                    pf = output.transition_metrics.get_false_vacuum_fraction(tt);
                    nt =  output.transition_metrics.get_nucleation_rate(tt);
                    n = get_n(tt, output.transition_metrics);
                    Rs = std::pow(n, -1./3.) * H;
                    Rbar = get_Rbar_integral(tt, output.transition_metrics)/n * H; 
                } catch (const std::exception& e) {
                    LOG(debug) << "Error computing thermal profile values at T = " << tt << ": " << e.what();
                    continue;
                } catch (...) {
                    LOG(debug) << "Unknown error computing thermal profile values at T = " << tt;
                    continue;
                }

                profile_out.temperature.push_back(tt);
                profile_out.dtdT.push_back(dtdT);
                profile_out.time.push_back(dt);
                profile_out.hubble_rate.push_back(H);
                profile_out.bounce_action.push_back(action);
                profile_out.false_vacuum_decay_rate.push_back(gamma);
                profile_out.extended_volume.push_back(vext);
                profile_out.false_vacuum_fraction.push_back(pf);
                profile_out.nucleation_rate.push_back(nt);
                profile_out.mean_bubble_separation.push_back(Rs);
                profile_out.mean_bubble_radius.push_back(Rbar);
            }
            output.profiles = profile_out;
        }

        return output;
    }

    const void
    ThermoFinder::add_thermal_parameter_values(TransitionMilestone& milestone, const FalseVacuumDecayRate& decay_rate, const EquationOfState& eos, TransitionMetrics& tm)
    {
        if(milestone.status == MilestoneStatus::YES) 
        {
            // TODO try catch with default values
            const auto temp = milestone.temperature;
            double alpha = get_alpha(temp, eos);
            milestone.alpha = alpha;
            double betaH = get_betaH(temp, decay_rate);
            milestone.betaH = betaH;
            double H = get_H(temp, tm);
            milestone.H = H;
            double we = get_we(temp, eos);
            milestone.we = we;
            std::pair<double, double> cs = get_cs(temp, eos);
            milestone.cs_plus = cs.first;
            milestone.cs_minus = cs.second;
            double n = get_n(temp, tm);
            milestone.n = n;
            milestone.Rs = std::pow(n, -1./3.) * H;
            milestone.Rbar = get_Rbar_integral(temp, tm)/n * H;

            double dt;
            dt = get_dt(temp, tm);
            milestone.dt = dt;
        } else {
            return;
        }
    }

    const double 
    ThermoFinder::get_alpha(const double& temperature, const EquationOfState& eos)
    {
        const auto theta = eos.get_theta(temperature);
        const auto w = eos.get_enthalpy_plus(temperature);
        return abs(theta.first - theta.second)/w * 4./3.;
    }

    const double
    ThermoFinder::get_betaH(const double& temperature, const FalseVacuumDecayRate& decay_rate)
    {
        double dy = decay_rate.get_action_deriv(temperature);
        return temperature * dy;
    }

    const double
    ThermoFinder::get_H(const double& temperature, TransitionMetrics& tm)
    {
        return tm.get_hubble_rate(temperature);
    }

    const double 
    ThermoFinder::get_we(const double& temperature, const EquationOfState& eos)
    {
        const double w_p = eos.get_enthalpy_plus(temperature);
        const double e_p = eos.get_energy_plus(temperature);
        return w_p/e_p;
    }

    const std::pair<double, double> 
    ThermoFinder::get_cs(const double& temperature, const EquationOfState& eos)
    {
        return eos.get_sound_speed(temperature);
    }

    const double 
    ThermoFinder::get_n(const double& temperature, TransitionMetrics& tm)
    {
        return tm.get_bubble_density(temperature);
    }

    const double 
    ThermoFinder::get_Rbar_integral(const double& temperature, TransitionMetrics& tm)
    {
        return tm.get_bubble_radius_integral(temperature);
    }

    const double
    ThermoFinder::get_dt(const double& temperature, TransitionMetrics& tm)
    {
        return tm.get_duration(temperature);
    }

} // namespace PhaseTracer