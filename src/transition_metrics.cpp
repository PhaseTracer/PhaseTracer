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
            // Avoid recursion when include_optimisations is false.
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
    TransitionMetrics::get_hubble_rate(const double& T, const double& e_true) const
    {
        double false_vacuum_fraction = 1.0;
        if (log_Vext_spline_computed)
        {
            const double Vext = get_extended_volume_from_spline(T);
            false_vacuum_fraction = std::exp(-Vext);
        }
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
    TransitionMetrics::make_scale_factor_ratio_spline() const
    {
        if (volume_term_integration_steps < 2)
        {
            scale_factor_spline_computed = false;
            return;
        }

        auto integrand = [this](double t)
        {
            auto potential = eos.eval_false_potential(t);
            return potential[2] / (3. * potential[1]);
        };

        integrate_and_fit_spline(scale_factor_spline, integrand, volume_term_integration_steps);
        scale_factor_spline_computed = true;

        // test, compare atop_abottom from spline to direct integration
        // evaluate at points not on the spline grid
        // for (int i = 0; i < 10; i++)
        // {
        //     double Ttop = t_min + (t_max - t_min) * (i + 0.5) / 10;
        //     double Tbottom = t_min + (t_max - t_min) * (i + 1.0) / 10;
        //     double ratio_spline = exp(alglib::spline1dcalc(scale_factor_spline, Tbottom) - alglib::spline1dcalc(scale_factor_spline, Ttop));
        //     double ratio_direct = get_atop_abottom(Ttop, Tbottom);
        //     std::cout << "Testing scale factor ratio spline at Ttop = " << Ttop << ", Tbottom = " << Tbottom << std::endl;
        //     std::cout << "Ratio from spline: " << ratio_spline << ", Ratio from direct integration: " << ratio_direct << std::endl;
        //     if (std::abs(ratio_spline - ratio_direct) > 1e-3)
        //     {
        //         LOG(warning) << "Scale factor ratio spline differs from direct integration by more than 0.1% at Ttop = " << Ttop << ", Tbottom = " << Tbottom;
        //         LOG(warning) << "Ratio from spline: " << ratio_spline << ", Ratio from direct integration: " << ratio_direct;
        //     }
        // }
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

        // compare volume term integral spline to direct integration
        // evaluate at points not on the spline grid
        // for (int i = 0; i < 10; i++)
        // {
        //     double T1 = t_min + (t_max - t_min) * (i + 0.5) / 10;
        //     double T2 = t_min + (t_max - t_min) * (i + 1.0) / 10;
        //     double K_spline_T1 = alglib::spline1dcalc(volume_term_integral_spline, T1);
        //     double K_spline_T2 = alglib::spline1dcalc(volume_term_integral_spline, T2);
        //     double K_direct_T1 = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(integrand, t_min, T1, 12, 1e-8);
        //     double K_direct_T2 = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(integrand, t_min, T2, 12, 1e-8);
        //     std::cout << "Testing volume term integral spline at T1 = " << T1 << ", T2 = " << T2 << std::endl;
        //     std::cout << "K from spline: " << K_spline_T1 << ", K from direct integration: " << K_direct_T1 << std::endl;
        //     std::cout << "K from spline: " << K_spline_T2 << ", K from direct integration: " << K_direct_T2 << std::endl;
        //     if (std::abs(K_spline_T1 - K_direct_T1) > 1e-3)
        //     {
        //         LOG(warning) << "Volume term integral spline differs from direct integration by more than 0.1% at T = " << T1;
        //         LOG(warning) << "K from spline: " << K_spline_T1 << ", K from direct integration: " << K_direct_T1;
        //     }
        //     if (std::abs(K_spline_T2 - K_direct_T2) > 1e-3)
        //     {
        //         LOG(warning) << "Volume term integral spline differs from direct integration by more than 0.1% at T = " << T2;
        //         LOG(warning) << "K from spline: " << K_spline_T2 << ", K from direct integration: " << K_direct_T2;
        //     }
        // }
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

    void
    TransitionMetrics::make_T_true_spline()
    {
        const int N = 500;
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
            LOG(debug)
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
                const double T_true = true ? get_T_true_matching(T_false, tol, max_iter) : T_false;
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
                    LOG(debug)
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
                LOG(debug)
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
            LOG(debug)
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

    void
    TransitionMetrics::solve_friedmann()
    {
        const int N = 250; // TODO

        const double dT_false = (t_min - t_max)/(N - 1);

        // cache values
        std::vector<double> T_false_grid(N);
        std::vector<double> T_true_grid(N);
        std::vector<double> t_grid(N);
        std::vector<double> false_vacuum_grid(N);
        std::vector<double> true_vacuum_grid(N);
        std::vector<double> e_false_grid(N);
        std::vector<double> e_true_grid(N);
        std::vector<double> p_true_grid(N);
        
        // initial conditions
        T_false_grid[N-1] = t_max;
        T_true_grid[N-1] = t_max;
        t_grid[N-1] = 0.0; // by definition
        false_vacuum_grid[N-1] = 1.0; // by definition
        true_vacuum_grid[N-1] = 0.0; // by definition
        e_false_grid[N-1] = std::abs(eos.get_energy_plus(t_max));
        e_true_grid[N-1] = std::abs(eos.get_energy_plus(t_max));
        p_true_grid[N-1] = std::abs(eos.get_pressure_minus(t_max));

        std::cout << "initial conditions: "
            << " T_false = " << T_false_grid[N-1]
            << ", T_true = " << T_true_grid[N-1]
            << ", Pf = " << false_vacuum_grid[N-1]
            << ", Pt = " << true_vacuum_grid[N-1]
            << ", e_false = " << e_false_grid[N-1]
            << ", e_true = " << e_true_grid[N-1]
            << ", p_true = " << p_true_grid[N-1]
            << std::endl;

        const double tol = 1e-8;
        boost::uintmax_t max_iter = 100;

        bool transition_complete = false;

        for(int i = N - 2; i >= 0; --i)
        {
            /*
                1. Step down in T_false by dT_false
            */
            const double T_false_prev = T_false_grid[i+1];
            const double T_false_current = T_false_prev + dT_false;
            T_false_grid[i] = T_false_current;


            /*
                2.  Evaluate the false vacuum fraction at new temp.

                    As a sanity check, we expect d(P_false)/d(T_false) to be positive, because as
                    we move towards T_f -> T_c, P_f increases towards 1. Likewise, d(P_true)/d(T_false)
                    should be negative, because as T_f -> T_c, P_t decreases towards 0.
            */
            const double false_vacuum_prev = false_vacuum_grid[i+1];
            const double false_vacuum_current = get_false_vacuum_fraction(T_false_current);
            const double true_vacuum_prev = 1.0 - false_vacuum_prev;
            const double true_vacuum_current = 1.0 - false_vacuum_current;
            const double d_false_vacuum = false_vacuum_current - false_vacuum_prev;
            const double d_true_vacuum = - d_false_vacuum;
            false_vacuum_grid[i] = false_vacuum_current;
            true_vacuum_grid[i] = true_vacuum_current;
            if ( d_false_vacuum > 0.0 )
            {
                LOG(warning) << "Warning: d(P_false) is positive at T_false = " << T_false_current;
            }
            if ( d_true_vacuum < 0.0 )
            {
                LOG(warning) << "Warning: d(P_true) is negative at T_false = " << T_false_current;
            }

            /*
                2.5 Check if the transition is complete. If so, we simply step down T_true.


            */
            if(!transition_complete && true_vacuum_current > 1.0 - 1e-6)
            {
                transition_complete = true;
                LOG(info) << "Transition complete at T_false = " << T_false_current;
            }

            if(transition_complete)
            {
                const double T_true_prev = T_true_grid[i+1];
                const double dT_true = dT_false; // same step size, same direction
                const double T_true_current = T_true_prev + dT_true;
                const double e_true_current = std::abs(eos.get_energy_minus(T_true_current));
                const double e_false_current = 0.0;
                T_true_grid[i] = T_true_current;
                e_true_grid[i] = e_true_current;
                e_false_grid[i] = e_false_current;
                p_true_grid[i] = std::abs(eos.get_pressure_minus(T_true_current));

                std::cout << "T_false = " << T_false_current
                    << ", T_true = " << T_true_current
                    << ", Pf = " << false_vacuum_current
                    << ", Pt = " << true_vacuum_current
                    << " [POST-TRANSITION]"
                    << ", e_true = " << e_true_current
                    << std::endl;
                continue;
            }

            /*
                3.  Evaluate the rate of change of d_true and the decay_rate dot(f)/f.
                    Set reheating flag based on the size of dot(f)/f. We manually make
                    this negative to agree with 2.
            */
            const double d_true_vacuum_dT_false = d_true_vacuum/dT_false;
            double transfer_rate = (true_vacuum_prev == 0.0) ? 0.0 : d_true_vacuum_dT_false/true_vacuum_prev; 
            const bool reheatingQ = std::abs(transfer_rate) > 1e-10;

            /*
                4.a If we are not reheating, solving the matching equation for T_true
                    this assumes instant thermalisation of the true vacuum.

                    As a sanity check, we expect e_true(T_false) < e_false(T_false). Hence, matching 
                    e(T_true) = e(T_false) forces T_true > T_false.
            */
            if(!reheatingQ)
            {
                const double T_true_current = get_T_true_matching(T_false_current, 1e-12, max_iter);
                const double e_false_current = std::abs(eos.get_energy_plus(T_false_current));
                const double e_true_current = std::abs(eos.get_energy_minus(T_true_current));
                const double p_true_current = std::abs(eos.get_pressure_minus(T_true_current));
                e_false_grid[i] = e_false_current;
                e_true_grid[i] = e_true_current;
                p_true_grid[i] = p_true_current;
                T_true_grid[i] = T_true_current;
                std::cout << "T_false = " << T_false_current
                    << ", T_true = " << T_true_current
                    << ", Pf = " << false_vacuum_current
                    << ", Pt = " << true_vacuum_current
                    << ", transfer_rate = " << transfer_rate
                    << (reheatingQ ? " [REHEATING]" : " [ADIABATIC]")
                    << ", e_false = " << e_false_current
                    << ", e_true = " << e_true_current
                    << std::endl;
                continue;
            }

            /*
                4.b If we are reheating, first evaluate false vacuum energy density.
            */
            const double e_false_current = std::abs(eos.get_energy_plus(T_false_current));
            e_false_grid[i] = e_false_current;

            /*
                5. Evaluate d(e_true)/dT_f. This has contributions from Hubble expansion
                   and the injected energy. 

                   The hubble term is - 3 H (e_t + p_t), which is negative such that as t
                   increases, e_t will lose energy due to redshift. But, changing 
                   d/dt = dT_f/dt d/dT_f = dot(T_f) d/dT_f introduces a factor 
                   of - 3 H (e_t + p_t) / dot(T_f). dot(T_f) is negative because 
                   such that the term is positive; as T_f increases, t decreases,
                   so the hubble term is positive.

                   The injected energy, dot(f)/f (e_f - e_t) is positive, provided
                   e_f > e_t (e_t is more stable). Hence, as time progresses, we gain 
                   energy due to reheating. But, we cancel the dT_f/dt term to get 
                   (df/dT_f)/f (e_f - e_t). As above, df/dT_f is negative, so as T_f
                   increases, t decreases, and we lose energy.
            */
            // double latent_heat = e_false_grid[i+1] - e_true_grid[i+1];
            double latent_heat = e_false_grid[i+1] - std::abs(eos.get_energy_minus(T_false_prev));
            if(std::abs(latent_heat) <= 1e-4) { latent_heat = 0.0; }
            const double injected_energy = transfer_rate * latent_heat;
            const double hubble = get_hubble_rate(T_false_prev);
            const double dt_dT_false = get_time_temperature_false(T_false_prev);
            const double redshifted_energy = -3.0 * dt_dT_false * hubble * (e_true_grid[i+1] + p_true_grid[i+1]);
            double d_e_true_dT_false = injected_energy + redshifted_energy;

            /*
                6. Update e_true by adding the above term.
            */
            const double d_e_true = d_e_true_dT_false * dT_false;
            const double e_true_current = e_true_grid[i+1] + d_e_true;
            e_true_grid[i] = e_true_current;

            /*
                7. Solve the matching condition for e_true, and update p_true
            */
            auto target_function = [this, e_true_current](double T_trial)
            {
                double e_true = abs(eos.get_energy_minus(T_trial));
                return e_true_current / e_true - 1.0;
            };

            auto root_pair = boost::math::tools::toms748_solve(
                target_function,
                eos.get_t_min(), eos.get_t_max(),
                [=](double l, double u){ return std::abs(u - l) < tol; },
                max_iter
            );

            const double T_true_current = (root_pair.first + root_pair.second) / 2.0;
            T_true_grid[i] = T_true_current;
            p_true_grid[i] = std::abs(eos.get_pressure_minus(T_true_current));

            std::cout << "T_false = " << T_false_current
                    << ", T_true = " << T_true_current
                    << ", e_false = " << e_false_current
                    << ", e_true = " << e_true_current
                    << ", Pt = " << true_vacuum_current
                    << (reheatingQ ? " [REHEATING]" : " [ADIABATIC]")
                    << ", latent_heat = " << latent_heat
                    << ", sign(injected) = " << (injected_energy > 0.0 ? "+" : "-")
                    << ", sign(redshifted) = " << (redshifted_energy > 0.0 ? "+" : "-")
                    << ", sign(d_e_true) = " << (d_e_true > 0.0 ? "+" : "-")
                    << std::endl;
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
            double e_true = abs(eos.get_energy_minus(T_true));
            double e_false = abs(eos.get_energy_plus(T_false));
            return e_false / e_true - 1.0;
        };

        // there should, in principle, always be a solution
        // this could be due to numerical error in the EOS
        double f_low = target_function(bracket.first);
        double f_high = target_function(bracket.second);
        if (f_low * f_high > 0.0)        
        {
            std::cout << "Target function does not change sign in the bracket. Returning T_false = " << T_false << std::endl;
            return T_false;
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
            make_scale_factor_ratio_spline();
            auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0);
            LOG(debug) << "  scale_factor_ratio_integrand spline: " << dt.count() << " ms";
        }
        

        for (int iter = 0; iter < max_extended_volume_refinements; ++iter)
        {
            LOG(debug) << "Running false vacuum iteration " << iter;

            {
                auto t0 = std::chrono::high_resolution_clock::now();
                compute_log_extended_volume_spline();
                log_Vext_spline_computed = true;
                auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0);
                LOG(debug) << "  time to build Vext spline: " << dt.count() << " ms"; 
            }

            if (iter > 0 && include_reheating)
            {
                auto t0 = std::chrono::high_resolution_clock::now();
                // make_T_true_spline();
                solve_friedmann();
                T_true_spline_computed = true;
                auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0);
                LOG(debug) << "  T_true spline: " << dt.count() << " ms";
            }

            if (!refine_extended_volume_spline) { break; }

            const auto perc = get_transition_milestone(MilestoneType::PERCOLATION);
            const double percolation_temperature = (perc.status == MilestoneStatus::YES || perc.status == MilestoneStatus::FAST) ? perc.temperature : t_min;
            const double d_perc_temp = std::abs(percolation_temperature - prev_percolation_temperature) / (std::abs(prev_percolation_temperature) + 1e-30);

            LOG(debug) << "  Tp = " << percolation_temperature << ", d(Tp) = " << d_perc_temp;

            if (iter > 0 && d_perc_temp < extended_volume_t_perc_tolerance)
            {
                LOG(debug) << "  false vacuum refinement converged after " << iter + 1 << " iterations.";
                break;
            }

            prev_percolation_temperature = percolation_temperature;
        }

        // write T_true as a function of T_false to file for debugging
        {
            std::ofstream file("example/TestThermalParameters/T_true_vs_T_false.csv");
            file << "# T_false,T_true,time,e_false,e_true,Pf,Pt\n";
            for (double T_false = t_min; T_false <= t_max; T_false += 0.1)
            {
                double T_true = T_true_spline_computed ? alglib::spline1dcalc(T_true_spline, T_false) : T_false;
                double e_false = std::abs(eos.get_energy_plus(T_false));
                double time = get_t(T_false);
                double e_true = get_e_true(T_false);
                double Pf = get_false_vacuum_fraction(T_false);
                double Pt = 1.0 - Pf;
                file << T_false << "," << T_true << "," << time << "," << e_false << "," << e_true << "," << Pf << "," << Pt << "\n";
            }
            file.close();
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