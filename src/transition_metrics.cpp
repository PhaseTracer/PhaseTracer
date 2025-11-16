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
        const double pressure_ratio = use_bag_dtdT ? 3. / T : pressure_derivs[2]/pressure_derivs[1];

        if (pressure_ratio < 0) { return 3./T;} // TODO

        return prefac * pressure_ratio;
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
        return use_bag_dtdT ? Tbottom/Ttop : exp(integral);
    }

    const double
    TransitionMetrics::get_volume_term(const double& T1, const double& T2)
    {
        auto integrand = [this, T2](double Tdash) 
        {
            double dtdT = get_dtdT(Tdash);
            double aT2_on_aTdash = get_atop_abottom(T2, Tdash);
            return dtdT * aT2_on_aTdash;
        };
        double result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, T1, T2, 5, 1e-5);
        return result;
    }

    const double 
    TransitionMetrics::extended_volume_integrand(const double& T1, const double& T2)
    {
        double dtdT = get_dtdT(T1);
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
        return 4 * M_PI * vw*vw*vw / 3 * result;
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
    TransitionMetrics::get_extended_volume_from_spline(const double& T)
    {
        return exp(alglib::spline1dcalc(log_Vext_spline, T));
    }

    const double
    TransitionMetrics::get_false_vacuum_fraction(const double& T)
    {
        double Vext = get_extended_volume_from_spline(T);
        return exp(-Vext);
    }

    const double 
    TransitionMetrics::get_nucleation_rate(const double& T)
    {
        auto integrand = [this](double Tdash) 
        {
            double Pf = use_pf_in_nt_integrand ? get_false_vacuum_fraction(Tdash) : 1.0;
            double gamma = decay_rate.get_gamma(Tdash);
            double hubble = get_hubble_rate(Tdash);
            double dtdT = get_dtdT(Tdash);
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
            double dtdT = get_dtdT(Tdash);
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
        double dtdT = get_dtdT(T1);
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
    TransitionMetrics::get_duration(const double& T)
    {
        auto integrand = [this](double Tdash) 
        {
            double H = get_hubble_rate(Tdash);
            double dtdT = get_dtdT(Tdash);
            return dtdT*H;
        };

        if(T >= t_max) { return 0.0; }
        
        double result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, t_max, T, 5, 1e-5);
        return vw * result;
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

        if(type == MilestoneType::PERCOLATION && output.status == MilestoneStatus::YES)
        {
            bool action_gradient_negative_at_t_min = decay_rate.get_action_deriv(t_min) > 0.0;
            if (action_gradient_negative_at_t_min)
            {
                output.nucleation_type = NucleationType::EXPONENTIAL;
            } else 
            {
                output.nucleation_type = NucleationType::SIMULTANEOUS;
            }
        }

        return output;
    }

} // namespace PhaseTracer