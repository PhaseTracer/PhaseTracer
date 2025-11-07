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
#include "transition_metrics.hpp"

namespace PhaseTracer {

    const double
    TransitionMetrics::get_hubble_rate(const double& T)
    {
        const double e_radiation = dof * M_PI * M_PI * T*T*T*T / 30.;
        const double e_vacuum = abs(eos.get_energy_plus(T) - eos.get_energy_minus(T)); // TODO

        const double Hsq = 8. * M_PI * newtonG/3. * (e_radiation + e_vacuum);

        return sqrt(Hsq);
    }

    const double
    TransitionMetrics::get_dtdT(const double& T)
    {
        const double prefac = -1./3. * 1/get_hubble_rate(T);

        const auto pressure_derivs = eos.get_pressure_derivs_plus(T);
        const double pressure_ratio = pressure_derivs[2]/pressure_derivs[1];

        return prefac * pressure_ratio;
    }

    void
    TransitionMetrics::make_a2a1_integrand_spline()
    {
        alglib::real_1d_array temp_array, integrand_array;
        temp_array.setlength(total_number_temp_steps);
        integrand_array.setlength(total_number_temp_steps);

        double dt = (t_max - t_min) / (total_number_temp_steps - 1);
        for (int i = 0; i < total_number_temp_steps; i++) {
            double tt = t_min + i * dt;
            double integrand = get_dtdT(tt) * get_hubble_rate(tt);
            temp_array[i] = tt;
            integrand_array[i] = integrand;
        }

        alglib::spline1dbuildcubic(temp_array, integrand_array, a2a1_integrand_spline);
    }

    const double
    TransitionMetrics::get_atop_abottom(const double& Ttop, const double& Tbottom)
    {
        double integral = alglib::spline1dintegrate(a2a1_integrand_spline, Tbottom) - alglib::spline1dintegrate(a2a1_integrand_spline, Ttop);
        return exp(integral);
    }

    const double
    TransitionMetrics::get_volume_term(const double& T1, const double& T2)
    {
        auto integrand = [this, T2](double Tdash) {
            double dtdT = get_dtdT(Tdash);
            double aT2_on_aTdash = get_atop_abottom(T2, Tdash);
            return dtdT * aT2_on_aTdash;
        };
        // double result = boost::math::quadrature::trapezoidal(integrand, T1, T2);
        double result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, T1, T2, 5, 1e-5);
        return result;
    }

    const double 
    TransitionMetrics::extended_volume_integrand(const double& T1, const double& T2)
    {
        double dtdT = get_dtdT(T1);
        double gamma = exp(alglib::spline1dcalc(decay_rate.log_gamma_spline, T1));
        double aT1_on_aT2 = get_atop_abottom(T1, T2);
        double volume_term = get_volume_term(T1, T2);

        return dtdT * gamma * aT1_on_aT2*aT1_on_aT2*aT1_on_aT2 * volume_term*volume_term*volume_term;
    }

    const double 
    TransitionMetrics::get_extended_volume(const double& T) 
    {
        auto integrand = [this, T](double Tdash) {
            return extended_volume_integrand(Tdash, T);
        };
        // double result = boost::math::quadrature::trapezoidal(integrand, t_max, T);
        double result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, t_max, T, 5, 1e-5);
        return 4 * M_PI * vw*vw*vw / 3 * result;
    }

    void
    TransitionMetrics::compute_Pf_spline()
    {
        alglib::real_1d_array temp_array, Pf_array, Nt_integrand_array;
        temp_array.setlength(total_number_temp_steps);
        Pf_array.setlength(total_number_temp_steps);
        Nt_integrand_array.setlength(total_number_temp_steps);

        double dt = (t_max - t_min) / (total_number_temp_steps - 1);
        for (int i = 0; i < total_number_temp_steps; i++) {
            double tt = t_min + i * dt;
            double Pf = exp(-get_extended_volume(tt));
            temp_array[i] = tt;
            Pf_array[i] = Pf;
        }

        alglib::spline1dbuildcubic(temp_array, Pf_array, Pf_spline);
    }

    const double 
    TransitionMetrics::get_nucleation_rate(const double& T)
    {
        auto integrand = [this](double Tdash) {
            double Pf = alglib::spline1dcalc(Pf_spline, Tdash);
            double gamma = exp(alglib::spline1dcalc(decay_rate.log_gamma_spline, Tdash));
            double hubble = get_hubble_rate(Tdash);
            double dtdT = get_dtdT(Tdash);
            return Pf * gamma * dtdT / (hubble*hubble*hubble);
        };

        double result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, t_max, T, 5, 1e-5);
        return 4 * M_PI * vw / 3 * result;
    }

    const double 
    TransitionMetrics::find_pf_temperature(const double& target) {

        auto func = [this, target] (double T) {
            return alglib::spline1dcalc(Pf_spline, T) - target;
        };

        std::pair<double, double> bracket = {t_min, t_max};

        boost::uintmax_t max_iter = 100;
        double tol = 1e-8;

        auto root_pair = boost::math::tools::toms748_solve(
            func,
            bracket.first,   // Lower bound
            bracket.second,
            [=](double l, double u){ return std::abs(u - l) < tol; },  // Termination condition
            max_iter         // Max iterations
        );

        double root = (root_pair.first + root_pair.second) / 2.0;  // Take midpoint of bracket

        return root;
    }

    const double 
    TransitionMetrics::find_nucleation_temperature() {

        auto func = [this] (double T) {
            return get_nucleation_rate(T) - 1.0;
        };

        std::pair<double, double> bracket = {t_min, t_max};

        boost::uintmax_t max_iter = 100;
        double tol = 1e-8;

        auto root_pair = boost::math::tools::toms748_solve(
            func,
            bracket.first,   // Lower bound
            bracket.second,
            [=](double l, double u){ return std::abs(u - l) < tol; },  // Termination condition
            max_iter         // Max iterations
        );

        double root = (root_pair.first + root_pair.second) / 2.0;  // Take midpoint of bracket

        return root;
    }

} // namespace PhaseTracer