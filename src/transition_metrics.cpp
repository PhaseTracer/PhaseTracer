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
    TransitionMetrics::get_hubble_rate(const double& true_vacuum_fraction, const double& e_false, const double& e_true) const
    {
        const double e_averaged = (1-true_vacuum_fraction)*e_false + true_vacuum_fraction*e_true;
        const double hubble_sq = 8. * M_PI * newtonG/3. * e_averaged;
        return std::sqrt(hubble_sq);
    }

    const double
    TransitionMetrics::get_time_temperature_false(const double& T_false) const
    {
        const double prefac = - 3. * T_false * get_hubble_rate(T_false);
        const double cs_false = eos.get_sound_speed_plus(T_false);
        const double cs_sq = cs_false*cs_false;
        const double dT_dt = prefac * cs_sq;
        return 1. / dT_dt;
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
    TransitionMetrics::get_false_vacuum_fraction_from_I3(const double& I3) const
    {
        double Veff = 4.0 * M_PI * vw*vw*vw / 3.0 * I3;
        return exp(-Veff);
    }

    const double
    TransitionMetrics::get_d_false_vacuum_fraction_from_I3(const double& I3, const double& I3_dot) const
    {
        double pf = get_false_vacuum_fraction_from_I3(I3);
        return - 4.0 * M_PI * vw*vw*vw / 3.0 * pf * I3_dot;
    }

    double 
    TransitionMetrics::match_T_true(const double& e_true, double tol, boost::uintmax_t max_iter)
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
    TransitionMetrics::match_T_false(const double& e_false, double tol, boost::uintmax_t max_iter)
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
        using state_type = std::array<double, 10>;

        auto rhs = [&](const state_type& state, state_type& dstate, double tau)
        {
            const double e_false = state[0];
            const double e_true  = state[1];
            const double a       = state[2];
            const double I_0     = state[3];
            const double I_1     = state[4];
            const double I_2     = state[5];
            const double I_3     = state[6];
            const double nucleation_rate = state[7];
            const double number_density = state[8];
            const double J = std::max(0.0, state[9]); // J = number_density * mean_bubble_radius

            const double T_true = match_T_true(e_true);
            const double T_false = match_T_false(e_false);

            const double p_false = eos.get_pressure_plus(T_false);
            const double p_true = eos.get_pressure_minus(T_true);
            const double latent_heat = e_false - e_true;

            const double false_vacuum_fraction = get_false_vacuum_fraction_from_I3(I_3);
            const double true_vacuum_fraction = 1 - false_vacuum_fraction;

            const double e_average = false_vacuum_fraction * e_false + true_vacuum_fraction * e_true;
            const double x_f = false_vacuum_fraction * e_false / e_average;
            const double x_t = true_vacuum_fraction * e_true / e_average;
            // LOG(debug) << "Tau = " << tau << ", x_f = " << x_f << ", x_t = " << x_t;

            const double hubble = get_hubble_rate(true_vacuum_fraction, e_false, e_true);
            const double gamma = decay_rate.get_gamma(T_false);

            const double time = std::exp(tau);

            dstate[2] = time * a * hubble;

            dstate[3] = time * (gamma       - 3.0 * hubble * I_0); // d(I_0)/d(ln t)
            dstate[4] = time * (      I_0   - 2.0 * hubble * I_1); // d(I_1)/d(ln t)
            dstate[5] = time * (2.0 * I_1   - 1.0 * hubble * I_2); // d(I_2)/d(ln t)
            dstate[6] = time * (3.0 * I_2);                        // d(I_3)/d(ln t)

            const double deriv_true_vacuum_fraction = 4.0/3.0*M_PI*vw*vw*vw * false_vacuum_fraction * dstate[6];
            const double reheating = (true_vacuum_fraction < 1e-30) ? 0.0 : deriv_true_vacuum_fraction/true_vacuum_fraction * latent_heat;

            dstate[0] = time * (- 3.0 * hubble * (e_false + p_false));           // d(e_false)/d(ln t)
            dstate[1] = time * (- 3.0 * hubble * (e_true + p_true)) + reheating; // d(e_true)/d(ln t)

            dstate[7] = time * (4.0/3.0 * M_PI * gamma * false_vacuum_fraction / (hubble*hubble*hubble));
            dstate[8] = time * (- 3.0 * hubble * number_density + gamma * false_vacuum_fraction);
            dstate[9] = (number_density < 1e-100) ? 0.0 : time * (number_density - 2.0 * hubble * J);
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
            const double nucleation_rate = state[7];
            const double number_density = state[8];
            const double J = std::max(0.0, state[9]); // J = number_density * mean_bubble_radius
            const double mean_bubble_radius = (number_density > 1e-100) ? J / number_density : 0.0;

            const double T_false = match_T_false(e_false);
            const double T_true  = match_T_true(e_true);

            // Stop if either temperature has reached (or gone below) t_min
            if (T_false <= t_min || T_true <= t_min || e_false <= e_false_min || e_true <= e_true_min)
            {
                throw TransitionCompleteException{};
            }

            const double p_false = eos.get_pressure_plus(T_false);
            const double p_true  = eos.get_pressure_minus(T_true);
            const double s_false = eos.get_entropy_plus(T_false);
            const double s_true  = eos.get_entropy_minus(T_true);

            const double false_vacuum_fraction = get_false_vacuum_fraction_from_I3(I_3);
            const double true_vacuum_fraction = 1 - false_vacuum_fraction;

            const double hubble = get_hubble_rate(true_vacuum_fraction, e_false, e_true);
            const double gamma  = decay_rate.get_gamma(T_false);

            const double t = std::exp(tau);

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
            system.nucleation_rate.push_back(nucleation_rate);
            system.number_density.push_back(number_density);
            system.mean_bubble_radius.push_back(mean_bubble_radius);
        };

        state_type initial_state = {eos.get_energy_plus(t_max), eos.get_energy_plus(t_max), 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

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
        } catch (const std::domain_error& e) {
            LOG(debug) << "Boost rootfinder error: integration stopped.";
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

        alglib::real_1d_array T_false_array;
        alglib::real_1d_array T_true_array; 
        alglib::real_1d_array log_time_array;
        alglib::real_1d_array scale_factor_array;
        alglib::real_1d_array hubble_rate_array;
        alglib::real_1d_array log_I_3_array;
        alglib::real_1d_array log_nucleation_rate_array;
        alglib::real_1d_array log_bubble_number_density_array;
        alglib::real_1d_array log_mean_bubble_radius_array;

        T_false_array.setlength(system.time.size());
        T_true_array.setlength(system.time.size());
        log_time_array.setlength(system.time.size());
        scale_factor_array.setlength(system.time.size());
        hubble_rate_array.setlength(system.time.size());
        log_I_3_array.setlength(system.time.size());
        log_nucleation_rate_array.setlength(system.time.size());
        log_bubble_number_density_array.setlength(system.time.size());
        log_mean_bubble_radius_array.setlength(system.time.size());

        for (std::size_t i = 0; i < system.time.size(); ++i)
        {
            T_false_array[i] = system.T_f[i];
            T_true_array[i] = system.T_t[i];
            log_time_array[i] = system.log_time[i];
            scale_factor_array[i] = system.a[i];
            hubble_rate_array[i] = system.hubble[i];
            log_I_3_array[i] = (i==0) ? -700 : std::log(system.I_3[i]);
            log_nucleation_rate_array[i] = (i==0) ? -700 : std::log(system.nucleation_rate[i]);
            log_bubble_number_density_array[i] = (i==0) ? -700 : std::log(system.number_density[i]);
            log_mean_bubble_radius_array[i] = (system.mean_bubble_radius[i]>0) ? std::log(system.mean_bubble_radius[i]) : -700;

            // LOG(debug) << "Friedmann data point " << i 
            //     << ": T_false = " << T_false_array[i] 
            //     << ", T_true = " << T_true_array[i] 
            //     << ", log_time = " << log_time_array[i] 
            //     << ", scale_factor = " << scale_factor_array[i] 
            //     << ", hubble_rate = " << hubble_rate_array[i] 
            //     << ", log_I_3 = " << log_I_3_array[i] 
            //     << ", log_nucleation_rate = " << log_nucleation_rate_array[i] 
            //     << ", log_bubble_number_density = " << log_bubble_number_density_array[i] 
            //     << ", log_mean_bubble_radius = " << log_mean_bubble_radius_array[i];
        }

        alglib::spline1dbuildcubic(T_false_array, T_true_array, reheating_spline);
        alglib::spline1dbuildcubic(T_false_array, log_time_array, log_time_spline);
        alglib::spline1dbuildcubic(T_false_array, scale_factor_array, scale_factor_spline);
        alglib::spline1dbuildcubic(T_false_array, hubble_rate_array, hubble_rate_spline);
        alglib::spline1dbuildcubic(T_false_array, log_I_3_array, log_I_3_spline);
        alglib::spline1dbuildcubic(T_false_array, log_nucleation_rate_array, log_nucleation_rate_spline);
        alglib::spline1dbuildcubic(T_false_array, log_bubble_number_density_array, log_bubble_number_density_spline);
        alglib::spline1dbuildcubic(T_false_array, log_mean_bubble_radius_array, log_mean_bubble_radius_spline);

        friedmann_splines_computed = true;
    }

    const double
    TransitionMetrics::get_scale_factor(const double& T_false) const
    {
        if(!friedmann_splines_computed)
        {
            LOG(warning) << "Friedmann splines not computed. Cannot compute scale factor.";
            return 1.0;
        }
        double scale_factor = alglib::spline1dcalc(scale_factor_spline, T_false);
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
    TransitionMetrics::get_hubble_rate(const double& T_false) const
    {
        if(friedmann_splines_computed)
        {
            double H = alglib::spline1dcalc(hubble_rate_spline, T_false);
            return H;
        }
        const double e_false = abs(eos.get_energy_plus(T_false));
        const double H_sq = 8. * M_PI * newtonG/3. * e_false;
        return std::sqrt(H_sq);
    }

    const double
    TransitionMetrics::get_false_vacuum_fraction(const double& T_false) const
    {
        if(!friedmann_splines_computed)
        {
            LOG(warning) << "Friedmann splines not computed. Cannot compute false vacuum fraction.";
            return 1.0;
        }
        double log_I_3 = alglib::spline1dcalc(log_I_3_spline, T_false);
        double I_3 = std::exp(log_I_3);
        return get_false_vacuum_fraction_from_I3(I_3);
    }

    const double
    TransitionMetrics::get_nucleation_rate(const double& T_false) const
    {
        if(!friedmann_splines_computed)
        {
            LOG(warning) << "Friedmann splines not computed. Cannot compute nucleation rate.";
            return 0.0;
        }
        double log_N = alglib::spline1dcalc(log_nucleation_rate_spline, T_false);
        return std::exp(log_N);
    }

    const double
    TransitionMetrics::get_bubble_density(const double& T_false) const
    {
        if(!friedmann_splines_computed)
        {
            LOG(warning) << "Friedmann splines not computed. Cannot compute bubble number density.";
            return 0.0;
        }
        double log_n = alglib::spline1dcalc(log_bubble_number_density_spline, T_false);
        return std::exp(log_n);
    }

    const double
    TransitionMetrics::get_mean_bubble_radius(const double& T_false) const
    {
        if(!friedmann_splines_computed)
        {
            LOG(warning) << "Friedmann splines not computed. Cannot compute mean bubble radius.";
            return 0.0;
        }
        double log_Rbar = alglib::spline1dcalc(log_mean_bubble_radius_spline, T_false);
        return std::exp(log_Rbar);
    }

    const double
    TransitionMetrics::get_t(const double& T_false) const
    {
        if(friedmann_splines_computed)
        {
            double log_t = alglib::spline1dcalc(log_time_spline, T_false);
            return std::exp(log_t);
        }

        auto integrand = [this](double Tdash) 
        {
            double dtdT = get_time_temperature_false(Tdash);
            return dtdT;
        };

        if(T_false >= t_max) { return 0.0; }
        
        double result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, t_max, T_false, 5, 1e-5);
        return result;
    }

    const double
    TransitionMetrics::get_T_true(const double& T_false) const
    {
        if(!friedmann_splines_computed)
        {
            LOG(warning) << "Friedmann splines not computed. Cannot compute T_true.";
            return T_false;
        }
        double T_true = alglib::spline1dcalc(reheating_spline, T_false);
        return T_true;
    }

    const LifetimeDistribution
    TransitionMetrics::get_lifetime_distribution(const double& timescale, const double& lifetime_min_fraction)
    {
        LifetimeDistribution distribution_out;

        const size_t n_grid = system.log_time.size();
        const double n_asymptotic = system.number_density.back();

        LOG(debug) << "Calculating lifetime distribution with timescale = " << timescale << " and lifetime_min_fraction = " << lifetime_min_fraction;
        LOG(debug) << "Number of grid points: " << n_grid;
        LOG(debug) << "Asymptotic number density: " << n_asymptotic;

        alglib::real_1d_array h_array, log_gamma_array, log_time_array;
        h_array.setlength(n_grid);
        log_time_array.setlength(n_grid);
        log_gamma_array.setlength(n_grid);
        for (size_t i = 0; i < n_grid; ++i)
        {
            double I3 = system.I_3[i];
            log_time_array[i] = system.log_time[i];
            h_array[i] = get_false_vacuum_fraction_from_I3(I3);
            log_gamma_array[i] = (system.gamma[i] > 0) ? std::log(system.gamma[i]) : -700;
        }

        alglib::spline1dbuildcubic(log_time_array, h_array, distribution_out.h_spline);
        alglib::spline1dbuildcubic(log_time_array, log_gamma_array, distribution_out.log_gamma_spline);

        const double log_t_min = log_time_array[0];
        const double log_t_max = log_time_array[n_grid - 1];
        const double t_max_linear = std::exp(log_t_max);

        // build lifetime grid from lifetime_min_fraction * timescale to 10 * timescale, with 100 points spaced logarithmically
        const double lifetime_min = lifetime_min_fraction * timescale;
        const double lifetime_max = 100.0 * timescale;
        const size_t n_lifetime_grid = 500;
        std::vector<double> lifetime_grid(n_lifetime_grid);
        for (size_t i = 0; i < n_lifetime_grid; ++i)
        {
            lifetime_grid[i] = lifetime_min * std::pow(lifetime_max/lifetime_min, static_cast<double>(i)/static_cast<double>(n_lifetime_grid-1));
        }

        distribution_out.lifetime_values = lifetime_grid;
        distribution_out.distribution_values.resize(n_lifetime_grid);
        for (size_t k = 0; k < n_lifetime_grid; ++k)
        {
            const double tau = lifetime_grid[k];

            // Upper limit on t': nucleation must happen early enough that t' + tau <= t_max_linear
            const double t_prime_max = t_max_linear - tau;
            if (t_prime_max <= std::exp(log_t_min))
            {
                distribution_out.distribution_values[k] = 0.0;
                LOG(debug) << "Lifetime grid point " << k << ": tau/timescale = " << tau/timescale << ", out of range — setting to 0";
                continue;
            }
            const double log_t_prime_max = std::log(t_prime_max);

            auto integrand = [&](double s) -> double
            {
                const double t_prime = std::exp(s);
                const double t_double_prime = t_prime + tau;
                const double log_t_double_prime = std::log(t_double_prime);

                double h, dh_ds, d2h_ds2;
                alglib::spline1ddiff(distribution_out.h_spline, log_t_double_prime, h, dh_ds, d2h_ds2);
                double dh_dt = std::min(0.0, dh_ds / t_double_prime);

                double log_gamma = alglib::spline1dcalc(distribution_out.log_gamma_spline, s);
                double gamma = std::exp(log_gamma);

                return t_prime * gamma * dh_dt / (n_asymptotic) * timescale;
            };

            double integral = boost::math::quadrature::gauss_kronrod<double, 31>::integrate(integrand, log_t_min, log_t_prime_max, 5, 1e-5);
            distribution_out.distribution_values[k] = - integral ;
            distribution_out.lifetime_values[k] = tau;
            distribution_out.chi_values[k] = tau / timescale;

            LOG(debug) << "Lifetime grid point " << k 
                << ": tau = " << tau / timescale
                << ", distribution value = " << distribution_out.distribution_values[k];
        }

        return distribution_out;
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