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
#include <limits>
#include "logger.hpp"
#include "transition_metrics.hpp"

namespace PhaseTracer {

    const double
    TransitionMetrics::get_hubble_rate(const double& T) const
    {
        double false_vacuum_fraction = 1.0;
        if (log_Vext_spline_computed)
        {
            const double Vext = get_extended_volume_from_spline(T);
            false_vacuum_fraction = std::exp(-Vext);
        }
        double true_vacuum_fraction = 1.0 - false_vacuum_fraction;
        
        double e_false_vacuum = abs(eos.get_energy_plus(T));
        double e_true_vacuum = get_e_true(T);

        const double e_averaged = false_vacuum_fraction * e_false_vacuum + true_vacuum_fraction * e_true_vacuum;
        const double Hsq = 8. * M_PI * newtonG/3. * e_averaged;

        return sqrt(Hsq);
    }

    const double
    TransitionMetrics::get_hubble_rate(const double& true_vacuum_fraction, const double& e_false, const double& e_true) const
    {
        const double e_averaged = (1-true_vacuum_fraction)*e_false + true_vacuum_fraction*e_true;
        const double hubble_sq = 8. * M_PI * newtonG/3. * e_averaged;
        return std::sqrt(hubble_sq);
    }

    const double
    TransitionMetrics::get_time_temperature_false(const double& T) const
    {
        const double prefac = - 3. * T * get_hubble_rate(T);
        const double cs_false = eos.get_sound_speed_plus(T);
        const double cs_sq = cs_false*cs_false;
        const double dT_dt = prefac * cs_sq;
        return 1. / dT_dt;
    }

    double 
    TransitionMetrics::get_T_true(const double& e_true, double tol, boost::uintmax_t max_iter)
    {
        // LOG(debug) << "get_T_true called for e_true = " << e_true;
        std::pair<double, double> bracket = {t_min, t_max};

        const double e_true_max = eos.get_energy_minus(t_max);
        if (e_true > e_true_max)
        {
            return t_max;
        }

        auto target_function = [this, e_true](double T_true)
        {
            double e_true_eos = eos.get_energy_minus(T_true);
            return e_true_eos - e_true;
        };

        double f_low = target_function(bracket.first);
        double f_high = target_function(bracket.second);
        if (f_low * f_high > 0.0)        
        {
            return t_max;
        }

        auto root_pair = boost::math::tools::toms748_solve(
            target_function,
            bracket.first, bracket.second,
            [=](double l, double u){ return std::abs(u - l) < tol; },
            max_iter
        );

        double root = (root_pair.first + root_pair.second) / 2.0;
        return root;
    }

    double 
    TransitionMetrics::get_T_false(const double& e_false, double tol, boost::uintmax_t max_iter)
    {
        // LOG(debug) << "get_T_false called for e_false = " << e_false;
        std::pair<double, double> bracket = {t_min, t_max};

        const double e_false_min = eos.get_energy_plus(t_min);
        if (e_false < e_false_min)
        {
            return t_min;
        }

        auto target_function = [this, e_false](double T_false)
        {
            double e_false_eos = eos.get_energy_plus(T_false);
            return e_false_eos - e_false;
        };

        double f_low = target_function(bracket.first);
        double f_high = target_function(bracket.second);
        if (f_low * f_high > 0.0)        
        {
            return t_max;
        }

        auto root_pair = boost::math::tools::toms748_solve(
            target_function,
            bracket.first, bracket.second,
            [=](double l, double u){ return std::abs(u - l) < tol; },
            max_iter
        );

        double root = (root_pair.first + root_pair.second) / 2.0;
        return root;
    }

    void
    TransitionMetrics::evolve_friedmann()
    {
        const double e_false_min = eos.get_energy_plus(t_min);
        const double e_true_min = eos.get_energy_minus(t_min);

        namespace odeint = boost::numeric::odeint;
        using state_type = std::array<double, 7>;

        auto rhs = [&](const state_type& state, state_type& dstate, double tau)
        {
            const double e_false = state[0];
            const double e_true  = state[1];
            const double a       = state[2];
            const double I_0     = state[3];
            const double I_1     = state[4];
            const double I_2     = state[5];
            const double I_3     = state[6];

            const double T_true = get_T_true(e_true);
            const double T_false = get_T_false(e_false);

            const double p_false = eos.get_pressure_plus(T_false);
            const double p_true = eos.get_pressure_minus(T_true);
            const double latent_heat = e_false - eos.get_energy_minus(T_false);

            const double false_vacuum_fraction = get_false_vacuum_fraction_from_I3(I_3);
            const double true_vacuum_fraction = 1 - false_vacuum_fraction;

            const double hubble = get_hubble_rate(true_vacuum_fraction, e_false - e_true_min, e_true - e_true_min);
            const double gamma = decay_rate.get_gamma(T_false);

            const double time = std::exp(tau); // time is log-time: tau = ln(t)

            dstate[2] = time * a * hubble; // da/d(ln t) = a * H

            dstate[3] = time * (gamma       - 3.0 * hubble * I_0); // d(I_0)/d(ln t)
            dstate[4] = time * (      I_0   - 2.0 * hubble * I_1); // d(I_1)/d(ln t)
            dstate[5] = time * (2.0 * I_1   - 1.0 * hubble * I_2); // d(I_2)/d(ln t)
            dstate[6] = time * (3.0 * I_2);                        // d(I_3)/d(ln t)

            const double deriv_true_vacuum_fraction = 4.0/3.0*M_PI*vw*vw*vw * false_vacuum_fraction * dstate[6];
            const double reheating = (true_vacuum_fraction < 1e-30) ? 0.0 : deriv_true_vacuum_fraction/true_vacuum_fraction * latent_heat;

            dstate[0] = time * (- 3.0 * hubble * (e_false + p_false));           // d(e_false)/d(ln t)
            dstate[1] = time * (- 3.0 * hubble * (e_true + p_true)) + reheating; // d(e_true)/d(ln t)
        };

        auto observer = [&](const state_type& state, double tau)
        {
            const double e_false = state[0];
            const double e_true  = state[1];
            const double a       = state[2];
            const double I_0     = state[3];
            const double I_1     = state[4];
            const double I_2     = state[5];
            const double I_3     = state[6];

            const double T_false = get_T_false(e_false);
            const double T_true  = get_T_true(e_true);

            // Stop if either temperature has reached (or gone below) t_min
            if (T_false <= t_min || T_true <= t_min)
            {
                throw TransitionCompleteException{};
            }

            const double p_false = eos.get_pressure_plus(T_false);
            const double p_true  = eos.get_pressure_minus(T_true);
            const double s_false = eos.get_entropy_plus(T_false);
            const double s_true  = eos.get_entropy_minus(T_true);

            const double false_vacuum_fraction = get_false_vacuum_fraction_from_I3(I_3);
            const double true_vacuum_fraction = 1 - false_vacuum_fraction;

            const double hubble = get_hubble_rate(true_vacuum_fraction, e_false - e_true_min, e_true - e_true_min);
            const double gamma  = decay_rate.get_gamma(T_false);

            const double t      = std::exp(tau);

            system.log_time.push_back(tau);
            system.time.push_back(t);
            system.e_f.push_back(e_false);
            system.e_t.push_back(e_true);
            system.p_f.push_back(p_false);
            system.p_t.push_back(p_true);
            system.w_f.push_back(e_false + p_false);
            system.w_t.push_back(e_true  + p_true);
            system.s_f.push_back(s_false);
            system.s_t.push_back(s_true);
            system.T_f.push_back(T_false);
            system.T_t.push_back(T_true);
            system.hubble.push_back(hubble);
            system.a.push_back(a);
            system.gamma.push_back(gamma);
            system.I_0.push_back(I_0);
            system.I_1.push_back(I_1);
            system.I_2.push_back(I_2);
            system.I_3.push_back(I_3);
        };

        state_type initial_state = {eos.get_energy_plus(t_max), eos.get_energy_plus(t_max), 1.0, 0.0, 0.0, 0.0, 0.0};

        // set initial time to be just below critical temp.
        double d_temp = 1e-3*(t_max - t_min);
        double T_initial = t_max - d_temp;
        double initial_time = std::log(0 - d_temp*get_time_temperature_false(T_initial));

        // estimate final time to be at t_min, but add a buffer to ensure we capture the full transition
        double final_time = std::log(get_t(t_min)) + 1.0;

        double dt = (final_time - initial_time)/(250.0-1.0);

        LOG(debug) << "Starting Friedmann evolution from T_initial = " << T_initial << " GeV at time " << exp(initial_time) << " GeV^-1, with estimated final time " << exp(final_time) << " GeV^-1";

        try
        {
            odeint::integrate_adaptive
            (
                odeint::make_controlled<odeint::runge_kutta_dopri5<state_type>>(1e-6, 1e-6),
                rhs, 
                initial_state, 
                initial_time, 
                final_time, 
                dt, 
                observer
            );
        }
        catch (const TransitionCompleteException&) {
            LOG(debug) << "Transition complete: integration stopped.";
        }
    }

    const void
    TransitionMetrics::fit_friedmann_splines() const
    {
        if (system.time.empty())
        {
            friedmann_splines_computed = false;
            return;
        }

        alglib::real_1d_array log_time_array, T_false_array, T_true_array, scale_factor_array;
        log_time_array.setlength(system.time.size());
        T_false_array.setlength(system.time.size());
        T_true_array.setlength(system.time.size());
        scale_factor_array.setlength(system.time.size());

        for (std::size_t i = 0; i < system.time.size(); ++i)
        {
            log_time_array[i] = std::log(system.time[i]);
            T_false_array[i] = system.T_f[i];
            T_true_array[i] = system.T_t[i];
            scale_factor_array[i] = system.a[i];
        }

        alglib::spline1dbuildcubic(T_false_array, log_time_array, T_false_spline);
        alglib::spline1dbuildcubic(T_false_array, T_true_array, reheating_spline);
        alglib::spline1dbuildcubic(log_time_array, scale_factor_array, scale_factor_spline);

        friedmann_splines_computed = true;
    }

    const double
    TransitionMetrics::get_scale_factor(const double& T_false) const
    {
        double log_time = alglib::spline1dcalc(T_false_spline, T_false);
        double scale_factor = alglib::spline1dcalc(scale_factor_spline, log_time);
        return scale_factor;
    }

    const double
    TransitionMetrics::get_scale_factor_log_time(const double& log_t) const
    {
        double scale_factor = alglib::spline1dcalc(scale_factor_spline, log_t);
        return scale_factor;
    }

    const double 
    TransitionMetrics::get_scale_factor_ratio(const double& Ttop, const double& Tbottom) const
    {
        if (use_bag_dtdT) { return Tbottom/Ttop; }
        
        if (!friedmann_splines_computed)
        {
            LOG(warning) << "Friedmann splines not computed. Falling back to bag model approximation for scale factor ratio.";
            return Tbottom/Ttop;
        }

        double a_top = get_scale_factor(Ttop);
        double a_bottom = get_scale_factor(Tbottom);
        return a_bottom/a_top;
    }

    const double 
    TransitionMetrics::get_scale_factor_ratio_log_time(const double& log_t_top, const double& log_t_bottom) const
    {
        if (use_bag_dtdT) { return exp(log_t_bottom - log_t_top); }
        if (!friedmann_splines_computed)
        {
            LOG(warning) << "Friedmann splines not computed. Falling back to bag model approximation for scale factor ratio.";
            return exp(log_t_bottom - log_t_top);
        }
        double a_top = get_scale_factor_log_time(log_t_top);
        double a_bottom = get_scale_factor_log_time(log_t_bottom);
        return a_bottom/a_top;
    }

    void
    TransitionMetrics::calculate_distributions() // TODO 4th order
    {
        double nucleation_rate = 0.0;
        double bubble_density_int = 0.0;

        const size_t n = system.log_time.size();
        system.nucleation_rate.resize(n);
        system.number_density.resize(n);

        for(size_t i = 0; i < n; ++i)
        {
            double d_log_time = (i==0) ? system.log_time[0] : system.log_time[i] - system.log_time[i-1];
            double time = std::exp(system.log_time[i]);

            double false_vacuum_fraction = std::exp(-4.0/3.0 * M_PI * vw*vw*vw * system.I_3[i]);
            double gamma = system.gamma[i];
            double hubble = system.hubble[i];
            double scale_factor = system.a[i];

            double d_nucleation_rate = 4.0/3.0 * M_PI * time * gamma * false_vacuum_fraction / (hubble*hubble*hubble) * d_log_time;
            double d_number_density_int = time * gamma * false_vacuum_fraction * scale_factor*scale_factor*scale_factor * d_log_time;

            nucleation_rate  += d_nucleation_rate;
            bubble_density_int += d_number_density_int;

            system.nucleation_rate[i] = nucleation_rate;
            system.number_density[i]  = bubble_density_int / (scale_factor*scale_factor*scale_factor);
        }
    }

    std::vector<double>
    TransitionMetrics::calculate_lifetime_distribution(const double beta, const std::vector<double>& lifetime_grid) const
    {
        const size_t n_grid = system.log_time.size();
        const double n_asymptotic = system.number_density.back();
        const double Vext_prefac = 4.0/3.0 * M_PI * vw*vw*vw;

        // Build spline for I_3(log_t) so we can interpolate and differentiate
        alglib::real_1d_array log_time_arr, I3_arr;
        log_time_arr.setlength(n_grid);
        I3_arr.setlength(n_grid);
        for (size_t i = 0; i < n_grid; ++i)
        {
            log_time_arr[i] = system.log_time[i];
            I3_arr[i]       = system.I_3[i];
        }
        alglib::spline1dinterpolant I3_spline;
        alglib::spline1dbuildcubic(log_time_arr, I3_arr, I3_spline);

        const double log_t_min = system.log_time.front();
        const double log_t_max = system.log_time.back();

        // df/dt at arbitrary time t_eval (zero outside the grid range)
        auto df_dt = [&](double t_eval) -> double
        {
            const double log_t_eval = std::log(t_eval);
            if (log_t_eval < log_t_min || log_t_eval > log_t_max) return 0.0;

            double I3_val, dI3_dlnt, d2I3;
            alglib::spline1ddiff(I3_spline, log_t_eval, I3_val, dI3_dlnt, d2I3);

            const double h        = std::exp(-Vext_prefac * I3_val);   // false vacuum fraction
            const double dVdt     = Vext_prefac * dI3_dlnt / t_eval;   // dV_ext/dt
            return h * dVdt;
        };

        // For each lifetime value t, integrate over t' in the grid
        std::vector<double> result(lifetime_grid.size());
        for (size_t k = 0; k < lifetime_grid.size(); ++k)
        {
            const double tau = lifetime_grid[k];
            double integral = 0.0;

            for (size_t i = 0; i < n_grid; ++i)
            {
                const double t_prime  = system.time[i];
                const double d_log_t  = (i == 0) ? system.log_time[i] : system.log_time[i] - system.log_time[i-1];
                const double dt_prime = t_prime * d_log_t;

                integral += system.gamma[i] * df_dt(t_prime + tau) * dt_prime;
            }

            result[k] = integral / (beta * n_asymptotic);
        }

        return result;
    }

    const double
    TransitionMetrics::get_atop_abottom(const double& Ttop, const double& Tbottom) const
    {
        if(use_bag_dtdT) { return Tbottom/Ttop; }
        if(include_optimisations && scale_factor_spline_computed)
        {
            double A_Ttop = alglib::spline1dcalc(scale_factor_spline, Ttop);
            double A_Tbottom = alglib::spline1dcalc(scale_factor_spline, Tbottom);
            return exp(A_Tbottom - A_Ttop);
        }

        auto integrand = [this](double t)
        {
            auto potential = eos.eval_false_potential(t);
            return potential[2] / (3. * potential[1]);
        };

        double result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, Ttop, Tbottom, 5, 1e-5);
        return std::exp(result);
    }

    void
    TransitionMetrics::make_volume_term_integral_spline() const
    {
        if (volume_term_integration_steps < 2)
        {
            volume_term_integral_spline_computed = false;
            return;
        }

        auto integrand = [this](double t)
        {
            double log_a = use_bag_dtdT ? std::log(t) : alglib::spline1dcalc(scale_factor_spline, t);
            return get_time_temperature_false(t) * exp(-log_a);
        };

        integrate_and_fit_spline(volume_term_integral_spline, integrand, volume_term_integration_steps);
        volume_term_integral_spline_computed = true;
    }

    const double
    TransitionMetrics::get_volume_term(const double& T1, const double& T2) const
    {
        if (include_optimisations && volume_term_integral_spline_computed)
        {
            double log_a_T2 = alglib::spline1dcalc(scale_factor_spline, T2);
            double K_T1 = alglib::spline1dcalc(volume_term_integral_spline, T1);
            double K_T2 = alglib::spline1dcalc(volume_term_integral_spline, T2);
            return exp(log_a_T2) * (K_T2 - K_T1);
        }

        auto integrand = [this, T2](double Tdash) 
        {
            double dtdT = get_time_temperature_false(Tdash);
            double scale_factor_ratio = get_atop_abottom(T2, Tdash);
            return dtdT * scale_factor_ratio;
        };

        double result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, T1, T2, 5, 1e-5);
        return result;
    }

    const double 
    TransitionMetrics::extended_volume_integrand(const double& T1, const double& T2) const
    {
        double dtdT = get_time_temperature_false(T1);
        double gamma = decay_rate.get_gamma(T1);
        double aT1_on_aT2 = get_atop_abottom(T1, T2);
        double volume_term = get_volume_term(T1, T2);

        return dtdT * gamma * aT1_on_aT2*aT1_on_aT2*aT1_on_aT2 * volume_term*volume_term*volume_term;
    }

    const double 
    TransitionMetrics::get_extended_volume(const double& T) const
    {
        auto integrand = [this, T](double Tdash) 
        {
            return extended_volume_integrand(Tdash, T);
        };

        if(T >= t_max) { return 0.0; }

        double result = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(integrand, t_max, T, 12, 1e-6);
        return  result;
    }

    void
    TransitionMetrics::compute_log_extended_volume_spline()
    {
        if (include_optimisations)
        {
            make_volume_term_integral_spline();
        }
        else
        {
            volume_term_integral_spline_computed = false;
        }

        alglib::real_1d_array temp_array, Vext_array;
        temp_array.setlength(total_number_temp_steps);
        Vext_array.setlength(total_number_temp_steps);

        double dt = (t_max - t_min) / (total_number_temp_steps - 1);
        for (int i = 0; i < total_number_temp_steps; i++) 
        {
            double tt = t_min + i * dt;
            double log_Vext = log(get_extended_volume(tt));
            if(std::isnan(log_Vext) || std::isinf(log_Vext)) { log_Vext = -700; }
            temp_array[i] = tt;
            Vext_array[i] = log_Vext;
        }

        alglib::spline1dbuildcubic(temp_array, Vext_array, log_Vext_spline);
    }

    const double
    TransitionMetrics::get_extended_volume_from_spline(const double& T) const
    {
        return 4 * M_PI * vw*vw*vw / 3 * exp(alglib::spline1dcalc(log_Vext_spline, T));
    }

    const double
    TransitionMetrics::get_false_vacuum_fraction_from_I3(const double& I3) const
    {
        double Veff = 4.0 * M_PI * vw*vw*vw / 3.0 * I3;
        return exp(-Veff);
    }

    const double
    TransitionMetrics::get_false_vacuum_fraction(const double& T) const
    {
        double Vext = include_optimisations ? get_extended_volume_from_spline(T) : 4 * M_PI * vw*vw*vw / 3 * get_extended_volume(T);
        return exp(-Vext);
    }

    const double
    TransitionMetrics::get_d_false_vacuum_fraction_dT(const double& T) const
    {
        double Vext = get_extended_volume_from_spline(T);
        double y, dy, ddy;
        alglib::spline1ddiff(log_Vext_spline, T, y, dy, ddy);

        return -exp(-Vext) * Vext * dy;
    }

    void TransitionMetrics::refine_temperature_bounds()
    {
        const int N = 1000;
        const double dT = (t_max - t_min) / (N - 1);
        const double monotonicity_tol = 10.0 * temperature_abs_tol;

        double e_false_prev = eos.get_energy_plus(t_max);
        double e_true_prev  = eos.get_energy_minus(t_max);
        double p_false_prev = eos.get_pressure_plus(t_max);
        double p_true_prev  = eos.get_pressure_minus(t_max);

        for(int i = 1; i < N; ++i)
        {
            const double T = t_max - i * dT;

            const double e_false = eos.get_energy_plus(T);
            const double e_true  = eos.get_energy_minus(T);
            const double p_false = eos.get_pressure_plus(T);
            const double p_true  = eos.get_pressure_minus(T);

            const bool monotonic_decreasing_broken =
                e_false > e_false_prev + monotonicity_tol ||
                e_true  > e_true_prev  + monotonicity_tol ||
                p_false > p_false_prev + monotonicity_tol ||
                p_true  > p_true_prev  + monotonicity_tol;

            if(monotonic_decreasing_broken)
            {
                t_min = T + dT;
                LOG(debug) << "Refining temperature bounds: setting t_min to " << t_min << " GeV to keep e(T) and p(T) monotonically decreasing";
                break;
            }

            e_false_prev = e_false;
            e_true_prev  = e_true;
            p_false_prev = p_false;
            p_true_prev  = p_true;
        }
    }

    const double
    TransitionMetrics::get_e_true(const double& T_false) const
    {
        double T_true = T_true_spline_computed ? alglib::spline1dcalc(T_false_spline, T_false) : T_false;
        const double t_min_eos = eos.get_t_min();
        const double t_max_eos = eos.get_t_max();

        if (T_true < t_min_eos) { T_true = t_min_eos;}
        else if (T_true > t_max_eos) { T_true = t_max_eos; }

        double e_true = abs(eos.get_energy_minus(T_true));
        return e_true;
    }

    const double 
    TransitionMetrics::get_nucleation_rate(const double& T) // TODO
    {
        auto integrand = [this](double Tdash) 
        {
            double Pf = use_pf_in_nt_integrand ? get_false_vacuum_fraction(Tdash) : 1.0;
            double gamma = decay_rate.get_gamma(Tdash);
            double hubble = get_hubble_rate(Tdash);
            double dtdT = get_time_temperature_false(Tdash);
            return Pf * gamma * dtdT / (hubble*hubble*hubble);
        };

        if(T >= t_max) { return 0.0; }

        double result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, t_max, T, 5, 1e-5);
        return 4 * M_PI * vw / 3 * result;
    }

    const double
    TransitionMetrics::get_bubble_density(const double& T)
    {
        auto integrand = [this, T](double Tdash)
        {
            double dtdT = get_time_temperature_false(Tdash);
            double gamma = decay_rate.get_gamma(Tdash);
            double pf = get_false_vacuum_fraction(Tdash);
            double a_ratio = get_atop_abottom(Tdash, T);
            return dtdT * gamma * pf * a_ratio*a_ratio*a_ratio;
        };

        if(T >= t_max) { return 0.0; }

        double result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, t_max, T, 5, 1e-5);
        return result;
    }

    const double 
    TransitionMetrics::bubble_radius_integrand(const double& T1, const double& T2) const
    {
        double dtdT = get_time_temperature_false(T1);
        double gamma = decay_rate.get_gamma(T1);
        double pf = get_false_vacuum_fraction(T1);
        double aT1_on_aT2 = get_atop_abottom(T1, T2);
        double volume_term = get_volume_term(T1, T2);

        return dtdT * gamma *pf *  aT1_on_aT2*aT1_on_aT2*aT1_on_aT2 * volume_term;
    }

    const double 
    TransitionMetrics::get_bubble_radius_integral(const double& T)
    {
        auto integrand = [this, T](double Tdash) 
        {
            return bubble_radius_integrand(Tdash, T);
        };

        if(T >= t_max) { return 0.0; }

        double result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, t_max, T, 5, 1e-5);
        return vw * result;
    }

    const double
    TransitionMetrics::get_t(const double& T)
    {
        auto integrand = [this](double Tdash) 
        {
            double dtdT = get_time_temperature_false(Tdash);
            return dtdT;
        };

        if(T >= t_max) { return 0.0; }
        
        double result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, t_max, T, 5, 1e-5);
        return result;
    }

    const RadiiDistribution 
    TransitionMetrics::get_radii_distribution(const double& temperature)
    {
        int n = 200; // TODO
        std::vector<double> temp, radii, dndR, log_dndR;

        double H = get_hubble_rate(temperature);

        LOG(debug) << "Calculating radii distribution at T = " << temperature << " GeV";

        double delta_T = t_max - temperature;
        double log_min = std::log(1.0);
        double log_max = std::log(1.0 + delta_T);
        double dlog = (log_max - log_min) / (n - 1);
        
        for(int i = 0; i < n; ++i)
        {
            double log_val = log_min + i * dlog;
            double tt = temperature + (std::exp(log_val) - 1.0);
            
            double rad = get_volume_term(tt, temperature)*H;

            double gamma = decay_rate.get_gamma(tt);
            double pf = get_false_vacuum_fraction(tt);
            double a_ratio = get_atop_abottom(tt, temperature);

            double dndR_at_tt = (gamma * pf / vw * a_ratio*a_ratio*a_ratio*a_ratio)/(H*H*H*H);

            temp.push_back(tt);
            radii.push_back(rad);
            dndR.push_back(dndR_at_tt);
            log_dndR.push_back(std::log(dndR_at_tt));

            LOG(debug) << "Temp: " << tt << ", Radius: " << rad << ", dn/dR: " << dndR_at_tt;
        }

        RadiiDistribution output(temperature, temp, radii, dndR, log_dndR);

        return output;
    };

    const LifetimeDistribution
    TransitionMetrics::get_lifetime_distribution(const double& timescale, const double& lifetime_min_fraction)
    {
        if (timescale <= 0.0) {
            throw std::invalid_argument(
                "TransitionMetrics::get_lifetime_distribution requires a positive timescale.");
        }

        const int n = 200;

        // -------------------------------------------------------------------------
        // The maximum elapsed time is tau(t_min) — the full duration of the transition.
        // -------------------------------------------------------------------------
        const double max_time     = std::abs(get_t(t_min));
        const double max_lifetime = max_time * timescale;

        // -------------------------------------------------------------------------
        // Pre-compute log(gamma) on a temperature grid and build a spline.
        // Splining log(gamma) then exp()-ing inside the integrand avoids
        // overflow/underflow when gamma spans many orders of magnitude.
        // -------------------------------------------------------------------------
        const int    n_spline        = 400;
        const double LOG_GAMMA_FLOOR = -700.0;
        const double dT_spline       = (t_max - t_min) / (n_spline - 1);

        std::vector<double> T_grid(n_spline), log_gamma_grid(n_spline);
        for (int k = 0; k < n_spline; ++k)
        {
            const double T     = t_min + k * dT_spline;
            T_grid[k]          = T;
            const double g     = decay_rate.get_gamma(T);
            log_gamma_grid[k]  = (g > 0.0 && std::isfinite(g))
                ? std::max(std::log(g), LOG_GAMMA_FLOOR)
                : LOG_GAMMA_FLOOR;
        }

        alglib::real_1d_array alg_T, alg_log_gamma;
        alg_T.setcontent(n_spline, T_grid.data());
        alg_log_gamma.setcontent(n_spline, log_gamma_grid.data());
        alglib::spline1dinterpolant log_gamma_spline;
        alglib::spline1dbuildcubic(alg_T, alg_log_gamma, log_gamma_spline);

        auto gamma_at = [&](double T) -> double {
            return std::exp(alglib::spline1dcalc(log_gamma_spline, T));
        };

        // -------------------------------------------------------------------------
        // Build a duration spline: tau(T) = get_t(T), then invert it to
        // give T as a function of tau. This lets us cheaply find T_shifted from
        // tau(T_current) + Delta_tau without a root search per quadrature point.
        // -------------------------------------------------------------------------
        std::vector<double> duration_grid(n_spline);
        for (int k = 0; k < n_spline; ++k)
            duration_grid[k] = get_t(T_grid[k]);

        // duration_grid increases as k increases (lower T = more elapsed time).
        // Collect strictly increasing (tau, T) pairs for the inverse spline.
        std::vector<double> tau_inc, T_inc;
        tau_inc.reserve(n_spline);
        T_inc.reserve(n_spline);
        for (int k = n_spline - 1; k >= 0; --k)
        {
            const double tau = duration_grid[k];
            if (tau_inc.empty() || tau > tau_inc.back() + 1e-30)
            {
                tau_inc.push_back(tau);
                T_inc.push_back(T_grid[k]);
            }
        }

        if (tau_inc.size() < 2)
            throw std::runtime_error(
                "get_lifetime_distribution: duration is not monotone over [t_min, t_max].");

        alglib::real_1d_array alg_tau, alg_T_from_tau;
        alg_tau.setcontent(tau_inc.size(), tau_inc.data());
        alg_T_from_tau.setcontent(T_inc.size(), T_inc.data());
        alglib::spline1dinterpolant T_from_tau_spline;
        alglib::spline1dbuildcubic(alg_tau, alg_T_from_tau, T_from_tau_spline);

        const double tau_min = tau_inc.front();
        const double tau_max = tau_inc.back();

        // -------------------------------------------------------------------------
        // Build a log-spaced lifetime axis.
        // lifetime = 0 is excluded (log undefined); the smallest bin is one
        // linear-equivalent step above zero, i.e. max_lifetime / (n-1).
        // -------------------------------------------------------------------------
        const double lifetime_min = lifetime_min_fraction * max_lifetime;
        const double log_lmin     = std::log(lifetime_min);
        const double log_lmax     = std::log(max_lifetime);
        const double dlog_l       = (log_lmax - log_lmin) / (n - 1);

        std::vector<double> lifetime_axis(n);
        std::vector<double> lifetime_pdf(n, 0.0);
        for (int i = 0; i < n; ++i)
            lifetime_axis[i] = std::exp(log_lmin + i * dlog_l);

        // -------------------------------------------------------------------------
        // Sample the PDF at each lifetime bin.
        //
        // pdf(Delta_tau) = integral over T_current of:
        //   -dpf/dt(T_shifted) * gamma(T_current) * (a_current/a_shifted)^3 * dtdT(T_current)
        //
        // where T_shifted satisfies tau(T_shifted) = tau(T_current) + Delta_tau.
        //
        // get_atop_abottom(Ttop, Tbottom) = a(Tbottom)/a(Ttop).
        // T_current is cooler (later) than T_shifted (earlier/hotter), so:
        //   a_current/a_shifted = get_atop_abottom(T_current, T_shifted).
        //
        // dtdT(T_current) is the Jacobian converting dT -> dt; it is negative
        // (T decreases with time), which combined with -dpf/dt > 0 gives a
        // positive integrand without any manual sign flip.
        // -------------------------------------------------------------------------
        for (int i = 0; i < n; ++i)
        {
            const double Delta_tau = lifetime_axis[i] / timescale;

            const double tau_current_upper = tau_max - Delta_tau;
            if (tau_current_upper <= tau_min)
            {
                lifetime_pdf[i] = 0.0;
                continue;
            }

            const double T_current_lower = alglib::spline1dcalc(T_from_tau_spline, tau_current_upper);

            auto integrand = [&](double T_current) -> double
            {
                if (T_current < t_min || T_current > t_max) return 0.0;

                const double tau_current = get_t(T_current);
                const double tau_shifted = tau_current + Delta_tau;

                if (tau_shifted > tau_max || tau_shifted < tau_min) return 0.0;

                const double T_shifted = alglib::spline1dcalc(T_from_tau_spline, tau_shifted);
                if (T_shifted < t_min || T_shifted > t_max) return 0.0;

                const double dpf_dT      = get_d_false_vacuum_fraction_dT(T_shifted);
                const double dtdT_shift  = get_time_temperature_false(T_shifted);
                if (std::abs(dtdT_shift) < 1e-14) return 0.0;
                const double dpf_dt      = dpf_dT / dtdT_shift;

                const double gamma_c     = gamma_at(T_current);

                const double a_ratio     = get_atop_abottom(T_current, T_shifted);
                if (!std::isfinite(a_ratio) || a_ratio <= 0.0) return 0.0;

                const double dtdT_current = get_time_temperature_false(T_current);

                return -dpf_dt * gamma_c
                    * a_ratio * a_ratio * a_ratio
                    * dtdT_current;
            };

            const double result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(
                integrand, T_current_lower, t_max, 7, 1e-6);

            lifetime_pdf[i] = std::isfinite(result) ? result : 0.0;
        }

        // -------------------------------------------------------------------------
        // Normalise to unit area using the trapezoid rule.
        // The variable bin widths from the log-spaced axis are handled correctly
        // since we use lifetime_axis[i] - lifetime_axis[i-1] directly.
        // -------------------------------------------------------------------------
        double normalisation = 0.0;
        for (int i = 1; i < n; ++i)
            normalisation += 0.5 * (lifetime_pdf[i-1] + lifetime_pdf[i])
                                * (lifetime_axis[i] - lifetime_axis[i-1]);

        if (normalisation <= 0.0 || !std::isfinite(normalisation))
            throw std::runtime_error(
                "get_lifetime_distribution: normalisation failed ("
                + std::to_string(normalisation) + "). "
                "Check get_false_vacuum_fraction, get_gamma, and get_atop_abottom "
                "over [t_min=" + std::to_string(t_min)
                + ", t_max=" + std::to_string(t_max) + "].");

        for (double& v : lifetime_pdf)
            v /= normalisation;

        return LifetimeDistribution{lifetime_axis, lifetime_pdf};
    }

    const double 
    TransitionMetrics::find_temperature(std::function<double(double)> target_function, double tol, boost::uintmax_t max_iter)
    {
        std::pair<double, double> bracket = {t_min, t_max};

        auto root_pair = boost::math::tools::toms748_solve(
            target_function,
            bracket.first,
            bracket.second,
            [=](double l, double u){ return std::abs(u - l) < tol; },
            max_iter
        );

        double root = (root_pair.first + root_pair.second) / 2.0;
        return root;
    }

    std::function<double(double)> 
    TransitionMetrics::get_target_function(const MilestoneType type)
    {
        switch (type) 
        {
            case MilestoneType::ONSET:
                return [this](double T) {return get_false_vacuum_fraction(T) - onset_target;};
            case MilestoneType::PERCOLATION:
                return [this](double T) {return get_false_vacuum_fraction(T) - percolation_target;};
            case MilestoneType::COMPLETION:
                return [this](double T) {return get_false_vacuum_fraction(T) - completion_target;};
            case MilestoneType::NUCLEATION:
                // code assumes the target function is monotonically decreasing, so multiply by -1
                return [this](double T) {return -(get_nucleation_rate(T) - nucleation_target);};
            default:
                throw std::invalid_argument("Invalid MilestoneType provided.");
        }
    }

    const TransitionMilestone 
    TransitionMetrics::get_transition_milestone(const MilestoneType type)
    {
        auto target_function = get_target_function(type);

        TransitionMilestone output(type);

        const auto valid = valid_lower_bound(target_function);
        if(valid)
        {
            double t = find_temperature(target_function);
            if (t_max - t < 1e-8)
            {
                output.status = MilestoneStatus::FAST;
                output.temperature = t;
            } else
            {
                output.status = MilestoneStatus::YES;
                output.temperature = t;
            }
        } else 
        {
            output.status = MilestoneStatus::NO;
        }

        return output;
    }

    void
    TransitionMetrics::compute_nucleation_history(const double& t_min, const double& t_max)
    {
        if(percolation_milestone.status == PhaseTracer::MilestoneStatus::YES) // && nucleation_milestone.status == PhaseTracer::MilestoneStatus::YES
        {
            const double t_guess = percolation_milestone.temperature;

            // TODO this chops off edge effects, can be removed once a better action calc is ready.
            double lower_bound = t_min + 0.05 * (t_guess - t_min);
            double upper_bound = t_max - 0.05 * (t_max - t_guess);

            auto action_func = [this](double T){
                return this->decay_rate.get_action(T) / T;
            };

            int bits = std::numeric_limits<double>::digits;
            boost::uintmax_t max_iter = 100;

            NucleationType type_out;

            /*
                The algorithm below searches for a minimum in the action curve, between the 'safe' temperature
                bounds defined above (these just try to limit edge effects, but should be refined in future).

                If a minimum is found, but it is at the lower bound (indicating monotonic decrease), it 
                is assumed to be exponential. Otherwise, we assume it is simultaneous. This is a 
                simplification.
            */
            
            try 
            {
                LOG(debug) << "Finding minimum T_m between lower bound " << lower_bound << " and upper bound " << upper_bound;

                auto result = boost::math::tools::brent_find_minima(action_func, lower_bound, upper_bound, bits, max_iter);
                
                double T_m = result.first;
                LOG(info) << "Found minimum T_mm = " << T_m;

                if (abs(lower_bound - T_m) < 1e-4)
                {
                    LOG(debug) << "TM is at lower bound, indicating exponential nucleation. Computing beta.";
                    type_out = NucleationType::EXPONENTIAL;
                } else 
                {
                    LOG(debug) << "TM is not at lower bound, indicating simultaneous nucleation. Computing beta2.";
                    type_out = NucleationType::SIMULTANEOUS;
                    nucleation_history.T_m = T_m; // only store if successful.
                }
            } catch (const std::exception& e) 
            {
                LOG(error) << "Error during minima finding: " << e.what();

                /*
                    If the above fails, we simply check the gradient at t_min.
                    This is succeptible to outliers.
                */

                bool action_gradient_negative_at_t_min = decay_rate.get_action_deriv(t_min) > 0.0;
                if (action_gradient_negative_at_t_min)
                {
                    type_out = NucleationType::EXPONENTIAL;
                } else 
                {
                    type_out = NucleationType::SIMULTANEOUS;
                    nucleation_history.T_m = 0.0;
                }
            }

            nucleation_history.nucleation_type = type_out;
            percolation_milestone.nucleation_type = type_out;
            nucleation_milestone.nucleation_type = type_out;
        } else {
            nucleation_history.nucleation_type = NucleationType::EXPONENTIAL;
            percolation_milestone.nucleation_type = NucleationType::EXPONENTIAL;
            nucleation_milestone.nucleation_type = NucleationType::EXPONENTIAL;
        }
    }

    std::vector<double> 
    TransitionMetrics::cumulative_simpson(const std::function<double(double)>& integrand, const std::vector<double>& x, double F_initial) const
    {
        const int N = x.size();
        assert(N >= 2);

        auto simpson_step = [&integrand](double x_lo, double x_hi, double f_lo, double f_hi)
        {
            double h     = x_hi - x_lo;
            double f_mid = integrand(x_lo + h * 0.5);
            return (h / 6.0) * (f_lo + 4.0 * f_mid + f_hi);
        };

        std::vector<double> F(N);
        F[0] = F_initial;

        double f_prev = integrand(x[0]);
        for (int i = 1; i < N; ++i)
        {
            double f_hi = integrand(x[i]);
            F[i] = F[i-1] + simpson_step(x[i-1], x[i], f_prev, f_hi);
            f_prev = f_hi;
        }

        return F;
    }

    alglib::real_1d_array 
    TransitionMetrics::cumulative_simpson(const std::function<double(double)>& integrand, const alglib::real_1d_array& x, double F_initial) const
    {
        const int N = x.length();
        assert(N >= 2);

        auto simpson_step = [&integrand](double x_lo, double x_hi, double f_lo, double f_hi)
        {
            double h     = x_hi - x_lo;
            double f_mid = integrand(x_lo + h * 0.5);
            return (h / 6.0) * (f_lo + 4.0 * f_mid + f_hi);
        };

        alglib::real_1d_array F;
        F.setlength(N);
        F[0] = F_initial;

        double f_prev = integrand(x[0]);
        for (int i = 1; i < N; ++i)
        {
            double f_hi = integrand(x[i]);
            F[i] = F[i-1] + simpson_step(x[i-1], x[i], f_prev, f_hi);
            f_prev = f_hi;
        }

        return F;
    }

    void
    TransitionMetrics::integrate_and_fit_spline 
    (
        alglib::spline1dinterpolant& spline, 
        const std::function<double(double)>& integrand, 
        int steps,
        double F_initial
    ) const
    {
        alglib::real_1d_array temp_grid;
        temp_grid.setlength(steps);
        double dt = (t_max - t_min) / (volume_term_integration_steps - 1);
        for (int i = 0; i < volume_term_integration_steps; ++i)
            temp_grid[i] = t_min + i * dt;

        alglib::real_1d_array integral_array = cumulative_simpson(integrand, temp_grid);

        alglib::spline1dbuildlinear(temp_grid, integral_array, spline);
    }

} // namespace PhaseTracer