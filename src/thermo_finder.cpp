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

    const std::vector<ThermalParameterSet>&
    ThermoFinder::get_thermal_parameters()
    {
        find_thermal_parameters();
        return thermal_parameters;
    }

    void
    ThermoFinder::find_thermal_parameters()
    {
        if(calculated_thermal_parameters) { return; }

        if(!tf)
        {
            throw std::logic_error("find_thermal_parameters requires a ThermoFinder constructed with a TransitionFinder");
        }

        calculated_thermal_parameters = true;
        auto transitions = tf->get_transitions();

        auto valid_transitions = transition_filter ? transition_filter(transitions)
                                                   : default_transition_filter(transitions);

        LOG(debug) << "Found " << valid_transitions.size() << " valid transitions out of " << transitions.size() << " total transitions.";

        thermal_parameters.reserve(valid_transitions.size());

        for(const auto& t : valid_transitions)
        {
            try {
                thermal_parameters.push_back(get_thermal_parameter_set(t));
            } catch (const std::exception& e) {
                LOG(debug) << "Error computing thermal parameters for transition @ TC = " << t.TC << ": " << e.what();
            } catch (...) {
                LOG(debug) << "Unknown error computing thermal parameters for transition @ TC = " << t.TC;
            }
        }

        LOG(debug) << "Ran find_thermal_parameters.";
    }

    std::vector<Transition>
    ThermoFinder::default_transition_filter(const std::vector<Transition>& transitions)
    {
        std::vector<Transition> valid_transitions;
        
        if(default_validation_method==PhaseTracer::ValidateMethod::TEMP)
        {
            for(const auto& t : transitions)
            {
                const double t_max = t.TC;
                const double t_min = t.false_phase.T.front();

                if(t_max - t_min > temperature_threshold)
                {
                    valid_transitions.push_back(t);
                }
            }
        }

        else if(default_validation_method==PhaseTracer::ValidateMethod::VEV)
        {
            for(const auto& t : transitions)
            {
                const auto vev = t.true_vacuum - t.false_vacuum;
                if(vev.norm() > vev_threshold)
                {
                    valid_transitions.push_back(t);
                }
            }
        }

        return valid_transitions;
    }

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
            temperature_abs_tol,
            prefactor_function
        );

        output.onset = output.friedmann_evolution->onset_milestone;
        output.onset.set_print_setting(onset_print_setting);
        add_thermal_parameter_values(output.onset, *output.decay_rate, *output.eos, *output.friedmann_evolution);

        output.percolation = output.friedmann_evolution->percolation_milestone;
        output.percolation.set_print_setting(percolation_print_setting);

        // update the percolation temperature 
        if(update_percolation_temperature)
        {
            try{
                revise_percolation_temperature(output.percolation, *output.eos, *output.friedmann_evolution);
            } catch (const std::exception& e) {
                LOG(debug) << "Error updating percolation temperature: " << e.what();
            } catch (...) {
                LOG(debug) << "Unknown error updating percolation temperature.";
            }
        }
        
        add_thermal_parameter_values(output.percolation, *output.decay_rate, *output.eos, *output.friedmann_evolution);

        output.completion = output.friedmann_evolution->completion_milestone;
        output.completion.set_print_setting(completion_print_setting);
        add_thermal_parameter_values(output.completion, *output.decay_rate, *output.eos, *output.friedmann_evolution);

        output.nucleation = output.friedmann_evolution->nucleation_milestone;
        output.nucleation.set_print_setting(nucleation_print_setting);
        add_thermal_parameter_values(output.nucleation, *output.decay_rate, *output.eos, *output.friedmann_evolution);

        output.nucleation_history = output.friedmann_evolution->nucleation_history;
        fill_nucleation_history(output.nucleation_history, output.percolation, output.nucleation, *output.decay_rate, *output.friedmann_evolution);

        if(compute_profiles)
        {
            ThermalProfiles profile_out;
            double t_min = output.friedmann_evolution->get_t_min();
            double t_max = output.friedmann_evolution->get_t_max();
            double dt = (t_max - t_min)/(n_temp_profiles-1);

            for(double tt = t_min; tt < t_max; tt += dt)
            {
                double dtdT, dt, H, action, gamma, vext, pf, d_pf, nt, n, Rs, Rbar;

                try {
                    dtdT = output.friedmann_evolution->get_time_temperature_false(tt);
                    dt = get_dt(tt, *output.friedmann_evolution);
                    H = get_H(tt, *output.friedmann_evolution);
                    action = output.decay_rate->get_action(tt)/tt;
                    gamma = output.decay_rate->get_gamma(tt);
                    pf = output.friedmann_evolution->get_false_vacuum_fraction(tt);
                    vext = -log(pf);
                    // d_pf = output.friedmann_evolution->get_d_false_vacuum_fraction_dT(tt);
                    nt =  output.friedmann_evolution->get_nucleation_rate(tt);
                    n = get_n(tt, *output.friedmann_evolution);
                    Rs = std::pow(n, -1./3.) * H;
                    Rbar = get_Rbar(tt, *output.friedmann_evolution) * H; 
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
                profile_out.d_false_vacuum_fraction.push_back(d_pf);
                profile_out.nucleation_rate.push_back(nt);
                profile_out.mean_bubble_separation.push_back(Rs);
                profile_out.mean_bubble_radius.push_back(Rbar);
            }
            output.profiles = profile_out;
        }

        return output;
    }

    const void
    ThermoFinder::add_thermal_parameter_values(TransitionMilestone& milestone, const FalseVacuumDecayRate& decay_rate, const EquationOfState& eos, FriedmannEvolution& tm)
    {
        if(milestone.status == MilestoneStatus::YES) 
        {
            // TODO try catch with default values
            const auto temp = milestone.temperature;
            double alpha = get_alpha(temp, eos);
            milestone.alpha = alpha;
            double alpha_munu = get_alpha(temp, eos, true);
            milestone.alpha_munu = alpha_munu;
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
            milestone.Rbar = get_Rbar(temp, tm) * H;

            double betaH_eff = get_betaH_eff(vw, milestone.Rs);
            milestone.betaH_eff = betaH_eff;

            double dt;
            dt = get_dt(temp, tm);
            milestone.dt = dt * H;
        }
    }

    void
    ThermoFinder::fill_nucleation_history(
        NucleationHistory& history,
        TransitionMilestone& percolation, 
        TransitionMilestone& nucleation, 
        const FalseVacuumDecayRate& decay_rate, 
        FriedmannEvolution& tm)
    {
        
        if(percolation.status == MilestoneStatus::YES)
        {
            const double Tref = percolation.temperature;
            const double T_m = history.T_m;
            const double betaH_1 = get_betaH_1(Tref, decay_rate, tm);
            const double betaH_2 = get_betaH_2(Tref, decay_rate, tm);

            LOG(debug) << "betaH_1 (Tp) = " << betaH_1 << ", betaH_1 (Tn) = " << get_betaH_1(nucleation.temperature, decay_rate, tm);
            LOG(debug) << "betaH_2 (Tp) = " << betaH_2 << ", betaH_2 (Tm) = " << get_betaH_2(T_m, decay_rate, tm);

            history.betaH_1 = betaH_1;
            history.betaH_2 = betaH_2;
        }
    }

    const double 
    ThermoFinder::get_alpha(const double& temperature, const EquationOfState& eos, bool use_munu)
    {
        const auto theta = eos.get_theta(temperature, use_munu);
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
    ThermoFinder::get_betaH_eff(const double& vw, const double& RsH)
    {
        return std::pow(8.*M_PI, 1./3.) * vw/RsH;
    }

    const double
    ThermoFinder::get_betaH_1(const double& temperature, const FalseVacuumDecayRate& decay_rate, FriedmannEvolution& tm)
    {
        const auto betas = tm.get_action_expansion(temperature);
        const double H = tm.get_hubble_rate(temperature);
        return betas.first/H;
    }

    const double
    ThermoFinder::get_betaH_2(const double& temperature, const FalseVacuumDecayRate& decay_rate, FriedmannEvolution& tm)
    {
        const auto betas = tm.get_action_expansion(temperature);
        const double H = tm.get_hubble_rate(temperature);
        return betas.second/H;
    }

    const double
    ThermoFinder::get_H(const double& temperature, FriedmannEvolution& tm)
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
    ThermoFinder::get_n(const double& temperature, FriedmannEvolution& tm)
    {
        return tm.get_bubble_density(temperature);
    }

    const double 
    ThermoFinder::get_Rbar(const double& temperature, FriedmannEvolution& tm)
    {
        return tm.get_mean_bubble_radius(temperature);
    }

    const double
    ThermoFinder::get_dt(const double& temperature, FriedmannEvolution& tm)
    {
        return tm.get_t(temperature);
    }

    const double
    ThermoFinder::get_percolation_temperature_wrapper(const double& vw, const double& percolation_target, const FriedmannEvolution& tm)
    {
        auto tm_copy = tm;
        tm_copy.set_vw(vw);
        tm_copy.set_percolation_target(percolation_target);

        const auto percolation_milestone = tm_copy.get_transition_milestone(MilestoneType::PERCOLATION);
        if (percolation_milestone.status == MilestoneStatus::YES || percolation_milestone.status == MilestoneStatus::FAST)
        {
            return percolation_milestone.temperature;
        }

        std::ostringstream message;
        message << "Failed to find percolation temperature for vw = " << vw
                << " and percolation_target = " << percolation_target;
        throw std::runtime_error(message.str());
    }

    const double
    ThermoFinder::get_vw_wrapper(const double& temperature, const FriedmannEvolution& tm, const EquationOfState& eos)
    {
        // 2303.10171
        
        const double cb = get_cs(temperature, eos).first; 
        const double alpha = get_alpha(temperature, eos, false);
        const double T_true = tm.get_T_true(temperature);
        const double wb = eos.get_enthalpy_minus(T_true);
        const double wp = eos.get_enthalpy_plus(temperature);
        const double Psi = std::min(wp/wb, 1.0);

        LOG(debug) << "In vw calculation: cb = " << cb << ", alpha = " << alpha << ", Psi = " << Psi;

        const double vJ = cb * (1 + sqrt(3 * alpha * (1 - cb*cb + 3*cb*cb*alpha)))/(1 + 3*cb*cb*alpha);
        
        const double v_low = sqrt((3*alpha + Psi - 1)/(2*(2 - 3*Psi + Psi*Psi*Psi)));

        const double a = 0.2233;
        const double b = 1.704;

        const double v_high = vJ * (1 - a * pow(1-Psi, b)/alpha);

        if(isnan(vJ) || isnan(v_low) || isnan(v_high))
        {
            std::ostringstream message;
            message << "Failed to compute vw: vJ = " << vJ << ", v_low = " << v_low << ", v_high = " << v_high << ", cb = " << cb << ", alpha = " << alpha << ", Psi = " << Psi;
            throw std::runtime_error(message.str());
        }

        LOG(debug) << "In vw calculation: vJ = " << vJ << ", v_low = " << v_low << ", v_high = " << v_high;

        const double p = -3.433;

        const double vw = pow(pow(abs(v_low), p) + pow(abs(v_high), p), 1/p);

        return vw;
    }

    const void
    ThermoFinder::revise_percolation_temperature(TransitionMilestone& percolation, const EquationOfState& eos, const FriedmannEvolution& tm)
    {
        const double vw_initial = vw;
        const double percolation_temp_initial = percolation.temperature;

        double vw_updated = vw_initial;
        double percolation_temp_updated = percolation_temp_initial;

        const int max_iter = 100;
        const double tol      = 1e-10;
        bool converged = false;

        LOG(debug) << "Updated Tp and vw using fixed-point iteration. Initial guess: [vw, Tp] = [" << vw_initial << ", " << percolation_temp_initial << "]";

        for (int i = 0; i < max_iter; i++) 
        {
            const double vw_new = get_vw_wrapper(percolation_temp_updated, tm, eos);
            const double Tp_new = get_percolation_temperature_wrapper(vw_new, percolation_target, tm);

            const double delta_vw = std::abs(vw_new - vw_updated) / (std::abs(vw_updated) + 1e-30);
            const double delta_Tp = std::abs(Tp_new - percolation_temp_updated) / (std::abs(percolation_temp_updated) + 1e-30);

            vw_updated = vw_new;
            percolation_temp_updated = Tp_new;

            LOG(debug) << std::setprecision(10) << "Iteration " << i << ": [vw, Tp] = [" << vw_updated << ", " << percolation_temp_updated << "], |delta_vw| = " << delta_vw << ", |delta_Tp| = " << delta_Tp;

            if (delta_vw < tol && delta_Tp < tol) 
            {
                converged = true;
                LOG(debug) << "Converged after " << i + 1 << " iterations.\n";
                break;
            }
        }

        if (!converged) 
        {
            LOG(debug) << "Warning: fixed-point iteration did not converge within " << max_iter << " iterations.\n";
        }

        LOG(debug) << "Final values after iteration: [vw, Tp] = [" << vw_updated << ", " << percolation_temp_updated << "]\n";

        percolation.temperature = percolation_temp_updated;
        set_vw(vw_updated);
    }

    std::ostream &operator<<(std::ostream &o, const ThermoFinder &a)
    {
        if(!a.calculated_thermal_parameters)
        {
            o << "no thermal parameters calculated yet\n"
            << "\n";
            return o; 
        }

        if(a.thermal_parameters.empty())
        {
            o << "found no thermal parameters"
            << "\n";
            return o; 
        }

        o << "found " << a.thermal_parameters.size() << " thermal parameter set";
        if (a.thermal_parameters.size() > 1) {
            o << "s";
        }
        o << "\n";

        for (const auto &t : a.thermal_parameters) 
        {
            o << t << "\n";
        }

        return o;
    }

} // namespace PhaseTracer