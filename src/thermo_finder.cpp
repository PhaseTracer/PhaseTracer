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

        // update the percolation temperature 
        if(update_percolation_temperature)
        {
            try{
                revise_percolation_temperature(output.percolation, output.eos, output.transition_metrics);
            } catch (const std::exception& e) {
                LOG(debug) << "Error updating percolation temperature: " << e.what();
            } catch (...) {
                LOG(debug) << "Unknown error updating percolation temperature.";
            }
        }
        add_thermal_parameter_values(output.percolation, output.decay_rate, output.eos, output.transition_metrics);

        output.completion = output.transition_metrics.completion_milestone;
        output.completion.set_print_setting(completion_print_setting);
        add_thermal_parameter_values(output.completion, output.decay_rate, output.eos, output.transition_metrics);

        output.nucleation = output.transition_metrics.nucleation_milestone;
        output.nucleation.set_print_setting(nucleation_print_setting);
        add_thermal_parameter_values(output.nucleation, output.decay_rate, output.eos, output.transition_metrics);

        output.nucleation_history = output.transition_metrics.nucleation_history;
        // output.nucleation_history.set_print_setting(nucleation_history_print_setting);
        fill_nucleation_history(output.nucleation_history, output.percolation, output.nucleation, output.decay_rate, output.transition_metrics);

        if(compute_profiles)
        {
            ThermalProfiles profile_out;
            double t_min = output.transition_metrics.get_t_min();
            double t_max = output.transition_metrics.get_t_max();
            double dt = (t_max - t_min)/(n_temp_profiles-1);

            for(double tt = t_min; tt < t_max; tt += dt)
            {
                double dtdT, dt, H, action, gamma, vext, pf, d_pf, nt, n, Rs, Rbar;

                try {
                    dtdT = output.transition_metrics.get_time_temperature_false(tt);
                    dt = get_dt(tt, output.transition_metrics);
                    H = get_H(tt, output.transition_metrics);
                    action = output.decay_rate.get_action(tt)/tt;
                    gamma = output.decay_rate.get_gamma(tt);
                    pf = output.transition_metrics.get_false_vacuum_fraction(tt);
                    vext = -log(pf);
                    // d_pf = output.transition_metrics.get_d_false_vacuum_fraction_dT(tt);
                    nt =  output.transition_metrics.get_nucleation_rate(tt);
                    n = get_n(tt, output.transition_metrics);
                    Rs = std::pow(n, -1./3.) * H;
                    // Rbar = get_Rbar_integral(tt, output.transition_metrics)/n * H; 
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
    ThermoFinder::add_thermal_parameter_values(TransitionMilestone& milestone, const FalseVacuumDecayRate& decay_rate, const EquationOfState& eos, TransitionMetrics& tm)
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
            // milestone.Rbar = get_Rbar_integral(temp, tm)/n * H;

            double betaH_eff = get_betaH_eff(vw, milestone.Rs);
            milestone.betaH_eff = betaH_eff;

            double dt;
            dt = get_dt(temp, tm);
            milestone.dt = dt * H;
        }
    }

    void
    ThermoFinder::add_history_lists(
        NucleationHistory& history,
        const FalseVacuumDecayRate& decay_rate, 
        TransitionMetrics& tm
    )
    {
        const int n = 100;
        const auto T_min = decay_rate.get_t_min();
        const auto T_max = decay_rate.get_t_max();
        const double dt = (T_max - T_min)/(n - 1);
        double conformal_time_val = 0.0;
        double prev_t = 0.0;
        double prev_a = 0.0;
        bool first_iter = true;
        for(double TT = T_max; TT > T_min; TT -= dt)
        {
            double t = get_dt(TT, tm);
            // double a = tm.get_atop_abottom(T_max, TT);
            double a = tm.get_scale_factor_ratio(T_max, TT);

            if(first_iter)
            {
                conformal_time_val = 0.0;
                first_iter = false;
            }
            else
            {
                conformal_time_val += 0.5 * (1.0/prev_a + 1.0/a) * (t - prev_t);
            }

            history.temperature.push_back(TT);
            history.time.push_back(t);
            history.scale_factor.push_back(a);
            history.conformal_time.push_back(conformal_time_val);

            prev_t = t;
            prev_a = a;
        }
    }

    void
    ThermoFinder::fill_nucleation_history(
        NucleationHistory& history,
        TransitionMilestone& percolation, 
        TransitionMilestone& nucleation, 
        const FalseVacuumDecayRate& decay_rate, 
        TransitionMetrics& tm)
    {
        add_history_lists(history, decay_rate, tm);
        

        if(percolation.status == MilestoneStatus::YES)
        {
            const double Tref = percolation.temperature;

            if(percolation.nucleation_type == NucleationType::SIMULTANEOUS)
            {
                const double T_m = history.T_m;
                const double T_p = (percolation.status == MilestoneStatus::YES) ? percolation.temperature : Tref;
                const double T_n = (nucleation.status == MilestoneStatus::YES) ? nucleation.temperature : Tref;

                const double betaH_1_perc = get_betaH_1(T_p, decay_rate, tm);
                const double betaH_1_nuc = get_betaH_1(T_n, decay_rate, tm);
                const double RsH_exp_perc = get_RsH_exp(betaH_1_perc);
                const double RsH_exp_nuc = get_RsH_exp(betaH_1_nuc);

                const double Gamma_0_perc = get_gamma_on_H4(T_p, decay_rate, tm);
                const double Gamma_0_nuc = get_gamma_on_H4(T_n, decay_rate, tm);
                const double Gamma_M = get_gamma_on_H4(T_m, decay_rate, tm);

                const double betaH_2_perc = get_betaH_2(T_p, decay_rate, tm);
                const double betaH_2_nuc = get_betaH_2(T_n, decay_rate, tm);
                const double betaH_2_m = get_betaH_2(T_m, decay_rate, tm);
                const double RsH_sim_perc = get_RsH_sim(Gamma_0_perc, betaH_2_perc);
                const double RsH_sim_nuc = get_RsH_sim(Gamma_0_nuc, betaH_2_perc);
                const double RsH_sim_m = get_RsH_sim(Gamma_M, betaH_2_m);

                history.betaH_1_perc = betaH_1_perc;
                history.betaH_1_nuc = betaH_1_nuc;
                history.RsH_exp_perc = RsH_exp_perc;
                history.RsH_exp_nuc = RsH_exp_nuc;

                history.betaH_2_perc = betaH_2_perc;
                history.betaH_2_nuc = betaH_2_nuc;
                history.RsH_sim_perc = RsH_sim_perc;
                history.RsH_sim_nuc = RsH_sim_nuc;

                history.betaH_2_m = betaH_2_m;
                history.RsH_sim_m = RsH_sim_m;

                history.Gamma_0_perc = Gamma_0_perc;
                history.Gamma_0_nuc = Gamma_0_nuc;
                history.Gamma_M = Gamma_M;


            } else 
            {
                const double T_p = (percolation.status == MilestoneStatus::YES) ? percolation.temperature : Tref;
                const double T_n = (nucleation.status == MilestoneStatus::YES) ? nucleation.temperature : Tref;

                const double betaH_1_perc = get_betaH_1(T_p, decay_rate, tm);
                const double betaH_1_nuc = get_betaH_1(T_n, decay_rate, tm);
                const double RsH_exp_perc = get_RsH_exp(betaH_1_perc);
                const double RsH_exp_nuc = get_RsH_exp(betaH_1_nuc);
                const double Gamma_0_perc = get_gamma_on_H4(T_p, decay_rate, tm);
                const double Gamma_0_nuc = get_gamma_on_H4(T_p, decay_rate, tm);

                history.betaH_1_perc = betaH_1_perc;
                history.betaH_1_nuc = betaH_1_nuc;
                history.RsH_exp_perc = RsH_exp_perc;
                history.RsH_exp_nuc = RsH_exp_nuc;
                history.Gamma_0_perc = Gamma_0_perc;
                history.Gamma_0_nuc = Gamma_0_nuc;
            }
        }
    }

    const double
    ThermoFinder::get_gamma_on_H4(const double& temperature, const FalseVacuumDecayRate& decay_rate, TransitionMetrics& tm)
    {
        const double gamma_m = decay_rate.get_gamma(temperature);
        const double h_m = tm.get_hubble_rate(temperature);
        return gamma_m / (h_m*h_m*h_m*h_m);
    }

    const double
    ThermoFinder::get_RsH_sim(const double& gammaH4, const double& betaH)
    {
        const double nb = std::sqrt(2.*M_PI) * gammaH4 / betaH;
        return std::pow(nb, -1./3.);
    }

    const double 
    ThermoFinder::get_RsH_exp(const double& betaH)
    {
        // const double nb = (betaH*betaH*betaH) / (8.*M_PI*vw*vw*vw);
        return std::pow(8.*M_PI, 1./3.) * vw/betaH;
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
    ThermoFinder::get_betaH_1(const double& temperature, const FalseVacuumDecayRate& decay_rate, TransitionMetrics& tm)
    {
        const double dy = decay_rate.get_action_deriv(temperature);
        const double dtdT = tm.get_time_temperature_false(temperature);
        const double H = tm.get_hubble_rate(temperature);
        return - dy/(dtdT*H);
    }

    const double
    ThermoFinder::get_betaH_2(const double& temperature, const FalseVacuumDecayRate& decay_rate, TransitionMetrics& tm)
    {
        const double dtdT = tm.get_time_temperature_false(temperature);
        const double ddSdTT2 = decay_rate.get_action_double_deriv(temperature);
        const double H = tm.get_hubble_rate(temperature);
        return std::sqrt(ddSdTT2/(dtdT*dtdT*H*H));
    }

    const double
    ThermoFinder::get_decay_rate_FWHM(const double& target_maximum, const double& target_temperature, const FalseVacuumDecayRate& decay_rate, TransitionMetrics& tm)
    {
        int bits = std::numeric_limits<double>::digits;
        boost::uintmax_t max_iter = 100;

        const double half_max = 0.5*target_maximum;

        auto half_max_func = [&decay_rate, &tm, half_max](double T) 
        {
            return decay_rate.get_gamma(T) / pow(tm.get_hubble_rate(T), 4) - half_max; 
        };

        auto upper_root_pair = boost::math::tools::toms748_solve(half_max_func, target_temperature, decay_rate.get_t_max(), [=](double l, double u){ return std::abs(u - l) < 1e-4; }, max_iter);
        double upper_root = (upper_root_pair.first + upper_root_pair.second) / 2.0;

        auto lower_root_pair = boost::math::tools::toms748_solve(half_max_func, decay_rate.get_t_min(), target_temperature, [=](double l, double u){ return std::abs(u - l) < 1e-4; }, max_iter);
        double lower_root = (lower_root_pair.first + lower_root_pair.second) / 2.0;

        auto integrand = [&tm](double T) 
        {
            return tm.get_time_temperature_false(T) * tm.get_hubble_rate(T);
        };
        const double FWHM = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, upper_root, lower_root, 5, 1e-5);

        return FWHM;
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

    // const double 
    // ThermoFinder::get_Rbar_integral(const double& temperature, TransitionMetrics& tm)
    // {
    //     return tm.get_bubble_radius_integral(temperature);
    // }

    const double
    ThermoFinder::get_dt(const double& temperature, TransitionMetrics& tm)
    {
        return tm.get_t(temperature);
    }

    const double
    ThermoFinder::get_percolation_temperature_wrapper(const double& vw, const double& percolation_target, const TransitionMetrics& tm)
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
    ThermoFinder::get_vw_wrapper(const double& temperature, const TransitionMetrics& tm, const EquationOfState& eos)
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
    ThermoFinder::revise_percolation_temperature(TransitionMilestone& percolation, const EquationOfState& eos, const TransitionMetrics& tm)
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

} // namespace PhaseTracer