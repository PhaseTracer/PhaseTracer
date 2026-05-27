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
        double false_vacuum_fraction = log_Vext_spline_computed ? get_false_vacuum_fraction(T) : 1.0;
        double true_vacuum_fraction = 1.0 - false_vacuum_fraction;
        
        double e_false_vacuum = abs(eos.get_energy_plus(T));
        double e_true_vacuum = get_e_true(T);

        const double e_averaged = false_vacuum_fraction * e_false_vacuum + true_vacuum_fraction * e_true_vacuum;
        const double Hsq = 8. * M_PI * newtonG/3. * e_averaged;

        return sqrt(Hsq);
    }

    const double
    TransitionMetrics::get_hubble_rate(const double& T, const double& e_true) const
    {
        const double false_vacuum_fraction = log_Vext_spline_computed ? get_false_vacuum_fraction(T) : 1.0;
        const double true_vacuum_fraction = 1.0 - false_vacuum_fraction;
        const double e_false = abs(eos.get_energy_minus(T));

        const double e_bar = false_vacuum_fraction * e_false + true_vacuum_fraction * e_true;
        const double Hsq = 8. * M_PI * newtonG / 3. * e_bar;
        return std::sqrt(std::max(Hsq, 0.0));
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

    void
    TransitionMetrics::make_scale_factor_ratio_integrand_spline()
    {
        alglib::real_1d_array temp_array, integrand_array;
        temp_array.setlength(total_number_temp_steps);
        integrand_array.setlength(total_number_temp_steps);

        double dt = (t_max - t_min) / (total_number_temp_steps - 1);
        for (int i = 0; i < total_number_temp_steps; i++) 
        {
            double tt = t_min + i * dt;
            auto potential = eos.eval_false_potential(tt);
            double integrand = potential[2] / (3. * potential[1]);
            temp_array[i] = tt;
            integrand_array[i] = integrand;
        }

        alglib::spline1dbuildcubic(temp_array, integrand_array, a2a1_integrand_spline);
    }

    const double
    TransitionMetrics::get_atop_abottom(const double& Ttop, const double& Tbottom) const
    {
        double integral = alglib::spline1dintegrate(a2a1_integrand_spline, Tbottom) - alglib::spline1dintegrate(a2a1_integrand_spline, Ttop);
        return use_bag_dtdT ? Tbottom/Ttop : exp(integral);
    }

    const double
    TransitionMetrics::get_volume_term(const double& T1, const double& T2)
    {
        auto integrand = [this, T2](double Tdash) 
        {
            double dtdT = get_time_temperature_false(Tdash);
            double aT2_on_aTdash = get_atop_abottom(T2, Tdash);
            return dtdT * aT2_on_aTdash;
        };
        double result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, T1, T2, 5, 1e-5);
        return result;
    }

    const double 
    TransitionMetrics::extended_volume_integrand(const double& T1, const double& T2)
    {
        double dtdT = get_time_temperature_false(T1);
        double gamma = decay_rate.get_gamma(T1);
        double aT1_on_aT2 = get_atop_abottom(T1, T2);
        double volume_term = get_volume_term(T1, T2);

        return dtdT * gamma * aT1_on_aT2*aT1_on_aT2*aT1_on_aT2 * volume_term*volume_term*volume_term;
    }

    const double 
    TransitionMetrics::get_extended_volume(const double& T) 
    {
        auto integrand = [this, T](double Tdash) 
        {
            return extended_volume_integrand(Tdash, T);
        };

        if(T >= t_max) { return 0.0; }

        double result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, t_max, T, 5, 1e-5);
        return  result;
    }

    void
    TransitionMetrics::compute_log_extended_volume_spline()
    {
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
    TransitionMetrics::get_false_vacuum_fraction(const double& T) const
    {
        double Vext = get_extended_volume_from_spline(T);
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

    void
    TransitionMetrics::make_T_true_spline()
    {
        const int N = 100;
        const double dT_false = (t_max - t_min) / (N - 1);

        std::vector<double> T_false_grid(N);
        std::vector<double> T_true_grid(N);

        T_false_grid[N - 1] = t_max;
        T_true_grid[N - 1]  = t_max;

        const double tol = 1e-8;
        boost::uintmax_t max_iter = 100;

        bool printed_ic = false;
        bool printed_reheat_start = false;
        bool printed_reheat_end = false;
        bool adiabatic_seen = false;
        bool reheating_seen = false;
        bool last_reheat_valid = false;
        bool last_adiabatic_valid = false;
        double last_reheat_Pt = 0.0;
        double last_reheat_Pf = 0.0;
        double last_reheat_T_false = 0.0;
        double last_reheat_T_true = 0.0;
        double last_adiabatic_Pt = 0.0;
        double last_adiabatic_Pf = 0.0;
        double last_adiabatic_T_false = 0.0;
        double last_adiabatic_T_true = 0.0;

        if (!printed_ic)
        {
            const double Pf_ic = get_false_vacuum_fraction(t_max);
            const double Pt_ic = 1.0 - Pf_ic;
            LOG(info)
                << "[IC] "
                << "Pt = " << Pt_ic
                << ", Pf = " << Pf_ic
                << ", T_false = " << t_max
                << ", T_true = " << t_max;
            printed_ic = true;
        }

        for(int i = N - 2; i >= 0; --i)
        {
            const double T_false_prev = T_false_grid[i + 1];
            const double T_true_prev = T_true_grid[i + 1];

            // iterate T_false down
            const double T_false = T_false_prev - dT_false;
            T_false_grid[i] = T_false;

            // to track the evolutin of the false vacuum fraction
            const double Pf_prev = get_false_vacuum_fraction(T_false_prev);
            const double Pf = get_false_vacuum_fraction(T_false);
            const double Pt_prev = 1.0 - Pf_prev;
            const double Pt = 1.0 - Pf;
            const double dPt = Pt - Pt_prev;

            if(Pt < 1e-16)
            {
                const double T_true = false ? get_T_true_matching(T_false, tol, max_iter) : T_false;
                T_true_grid[i] = T_true;
                continue;
            }

            const double reheating_pt_min = 1e-8;
            if (Pt < 1.0 - 1e-16)
            {
                const double e_false = std::abs(eos.get_energy_plus(T_false));
                // const double e_old = std::abs(eos.get_energy_minus(T_true_prev));
                const double e_old = get_e_true(T_true_prev);
                const double target = e_old * (Pt - dPt) + dPt * e_false;

                auto reheating_eq = [this, Pt, target](double T_trial)
                {
                    // const double e_new = std::abs(eos.get_energy_minus(T_trial));
                    const double e_new = get_e_true(T_trial);
                    return e_new * Pt - target;
                };

                const double T_low = eos.get_t_min();
                const double T_high = eos.get_t_max();
                const int n_scan = 60;
                double best_T_abs = T_true_prev;
                double best_abs = std::numeric_limits<double>::infinity();
                double best_T_root = T_true_prev;
                double best_root_dist = std::numeric_limits<double>::infinity();
                bool found_root = false;

                if (T_low < T_high)
                {
                    double prev_T = T_low;
                    double prev_f = reheating_eq(prev_T);
                    best_T_abs = prev_T;
                    best_abs = std::abs(prev_f);
                    if (prev_f == 0.0)
                    {
                        found_root = true;
                        best_T_root = prev_T;
                        best_root_dist = std::abs(prev_T - T_true_prev);
                    }

                    for (int j = 1; j < n_scan; ++j)
                    {
                        const double T = T_low + (T_high - T_low) * static_cast<double>(j) / (n_scan - 1);
                        const double f = reheating_eq(T);

                        const double abs_f = std::abs(f);
                        if (abs_f < best_abs)
                        {
                            best_abs = abs_f;
                            best_T_abs = T;
                        }

                        if (f == 0.0)
                        {
                            const double dist = std::abs(T - T_true_prev);
                            if (!found_root || dist < best_root_dist)
                            {
                                found_root = true;
                                best_T_root = T;
                                best_root_dist = dist;
                            }
                            prev_T = T;
                            prev_f = f;
                            continue;
                        }

                        if ((prev_f < 0.0 && f > 0.0) || (prev_f > 0.0 && f < 0.0))
                        {
                            auto root_pair = boost::math::tools::toms748_solve(
                                reheating_eq,
                                prev_T, T,
                                [=](double l, double u){ return std::abs(u - l) < tol; },
                                max_iter
                            );

                            const double root = (root_pair.first + root_pair.second) / 2.0;
                            const double dist = std::abs(root - T_true_prev);
                            if (!found_root || dist < best_root_dist)
                            {
                                found_root = true;
                                best_T_root = root;
                                best_root_dist = dist;
                            }
                        }

                        prev_T = T;
                        prev_f = f;
                    }
                }

                T_true_grid[i] = found_root ? best_T_root : best_T_abs;

                reheating_seen = true;
                last_reheat_valid = true;
                last_reheat_Pt = Pt;
                last_reheat_Pf = Pf;
                last_reheat_T_false = T_false;
                last_reheat_T_true = T_true_grid[i];

                if (!printed_reheat_start)
                {
                    LOG(info)
                        << "[REHEATING START] "
                        << "Pt = " << Pt
                        << ", Pf = " << Pf
                        << ", T_false = " << T_false
                        << ", T_true = " << T_true_grid[i];
                    printed_reheat_start = true;
                }
                continue;
            }

            if (reheating_seen && last_reheat_valid && !printed_reheat_end)
            {
                LOG(info)
                    << "[REHEATING END] "
                    << "Pt = " << last_reheat_Pt
                    << ", Pf = " << last_reheat_Pf
                    << ", T_false = " << last_reheat_T_false
                    << ", T_true = " << last_reheat_T_true;
                printed_reheat_end = true;
            }

            const double T_true_adiabatic = get_T_true_adiabatic(T_false, T_false_prev, T_true_prev, tol, max_iter);
            T_true_grid[i] = T_true_adiabatic;

            adiabatic_seen = true;
            last_adiabatic_valid = true;
            last_adiabatic_Pt = Pt;
            last_adiabatic_Pf = Pf;
            last_adiabatic_T_false = T_false;
            last_adiabatic_T_true = T_true_adiabatic;
        }

        if (adiabatic_seen && last_adiabatic_valid)
        {
            LOG(info)
                << "[ADIABATIC END] "
                << "Pt = " << last_adiabatic_Pt
                << ", Pf = " << last_adiabatic_Pf
                << ", T_false = " << last_adiabatic_T_false
                << ", T_true = " << last_adiabatic_T_true;
        }

        alglib::real_1d_array T_false_arr, T_true_arr;
        T_false_arr.setcontent(N, T_false_grid.data());
        T_true_arr.setcontent(N, T_true_grid.data());
        alglib::spline1dbuildcubic(T_false_arr, T_true_arr, T_true_spline);
    }

    const double 
    TransitionMetrics::get_T_true_matching(const double& T_false, double tol, boost::uintmax_t max_iter) const
    {
        std::pair<double, double> bracket = {T_false, t_max};

        if (std::abs(T_false - t_max) < tol) 
        {
            return T_false;
        }

        auto target_function = [this, T_false](double T_true)
        {
            double e_true_from_ode = abs(eos.get_energy_minus(T_true));;
            double e_true_from_eos = abs(eos.get_energy_plus(T_false));
            return e_true_from_eos - e_true_from_ode;
        };

        auto root_pair = boost::math::tools::toms748_solve(
            target_function,
            bracket.first, bracket.second,
            [=](double l, double u){ return std::abs(u - l) < tol; },
            max_iter
        );

        double root = (root_pair.first + root_pair.second) / 2.0;
        return root;
    }

    const double 
    TransitionMetrics::get_T_true_adiabatic(const double& T_false, const double& T_false_prev, const double& T_true_prev, double tol, boost::uintmax_t max_iter) const
    {
        auto adiabatic_eq = [this, T_false, T_false_prev, T_true_prev](double T_trial)
        {
            const double s_trial = std::abs(eos.get_entropy_minus(T_trial));
            const double a_ratio = get_atop_abottom(T_false_prev, T_false);
            const double s_true = std::abs(eos.get_entropy_minus(T_true_prev));
            return s_trial - s_true * a_ratio*a_ratio*a_ratio;
        };

        const double T_low = eos.get_t_min();
        const double T_high = eos.get_t_max();
        if (T_low >= T_high) { return T_true_prev; }

        const int n_scan = 60;
        double best_T = T_true_prev;
        double best_abs = std::numeric_limits<double>::infinity();

        double prev_T = T_low;
        double prev_f = adiabatic_eq(prev_T);
        best_T = prev_T;
        best_abs = std::abs(prev_f);

        for (int i = 1; i < n_scan; ++i)
        {
            const double T = T_low + (T_high - T_low) * static_cast<double>(i) / (n_scan - 1);
            const double f = adiabatic_eq(T);

            const double abs_f = std::abs(f);
            if (abs_f < best_abs)
            {
                best_abs = abs_f;
                best_T = T;
            }

            if (prev_f == 0.0) { return prev_T; }
            if (f == 0.0) { return T; }

            if ((prev_f < 0.0 && f > 0.0) || (prev_f > 0.0 && f < 0.0))
            {
                auto root_pair = boost::math::tools::toms748_solve(
                    adiabatic_eq,
                    prev_T, T,
                    [=](double l, double u){ return std::abs(u - l) < tol; },
                    max_iter
                );

                return (root_pair.first + root_pair.second) / 2.0;
            }

            prev_T = T;
            prev_f = f;
        }

        return best_T;
    }

    const double
    TransitionMetrics::get_e_true(const double& T_false) const
    {
        double T_true = T_true_spline_computed ? alglib::spline1dcalc(T_true_spline, T_false) : T_false;
        const double t_min_eos = eos.get_t_min();
        const double t_max_eos = eos.get_t_max();

        if (T_true < t_min_eos) { T_true = t_min_eos;}
        else if (T_true > t_max_eos) { T_true = t_max_eos; }

        double e_true = abs(eos.get_energy_minus(T_true));
        return e_true;
    }

    void
    TransitionMetrics::calculate_false_vacuum_fraction()
    {
        double prev_percolation_temperature = t_max;

        {
            auto t0 = std::chrono::high_resolution_clock::now();
            make_scale_factor_ratio_integrand_spline();
            auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0);
            LOG(info) << "  scale_factor_ratio_integrand spline: " << dt.count() << " ms";
        }

        for (int iter = 0; iter < max_extended_volume_refinements; ++iter)
        {
            LOG(info) << "Running false vacuum iteration " << iter;

            {
                auto t0 = std::chrono::high_resolution_clock::now();
                compute_log_extended_volume_spline();
                log_Vext_spline_computed = true;
                auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0);
                LOG(info) << "  time to build Vext spline: " << dt.count() << " ms"; 
            }

            if (iter > 0 && include_reheating)
            {
                auto t0 = std::chrono::high_resolution_clock::now();
                make_T_true_spline();
                T_true_spline_computed = true;
                auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0);
                LOG(info) << "  T_true spline: " << dt.count() << " ms";
            }

            if (!refine_extended_volume_spline) { break; }

            const auto perc = get_transition_milestone(MilestoneType::PERCOLATION);
            const double percolation_temperature = (perc.status == MilestoneStatus::YES || perc.status == MilestoneStatus::FAST) ? perc.temperature : t_min;
            const double d_perc_temp = std::abs(percolation_temperature - prev_percolation_temperature) / (std::abs(prev_percolation_temperature) + 1e-30);

            LOG(info) << "  Tp = " << percolation_temperature << ", d(Tp) = " << d_perc_temp;

            if (iter > 0 && d_perc_temp < extended_volume_t_perc_tolerance)
            {
                LOG(info) << "  false vacuum refinement converged after " << iter + 1 << " iterations.";
                break;
            }

            prev_percolation_temperature = percolation_temperature;
        }
    }

    const double 
    TransitionMetrics::get_nucleation_rate(const double& T)
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
    TransitionMetrics::bubble_radius_integrand(const double& T1, const double& T2)
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

} // namespace PhaseTracer