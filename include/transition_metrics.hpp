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

#ifndef PHASETRACER_TRANSITION_METRICS_HPP_
#define PHASETRACER_TRANSITION_METRICS_HPP_

#include <cmath>
#include <chrono>
#include <vector>
#include <optional>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <interpolation.h>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/tools/roots.hpp>

#include "property.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "false_vacuum_decay_rate.hpp"
#include "equation_of_state.hpp"

namespace PhaseTracer {

enum MilestoneStatus
{
    YES,
    FAST,
    NO,
    ERR
}; // enum MilestoneStatus

enum MilestoneType
{
    PERCOLATION,
    NUCLEATION,
    COMPLETION,
    ONSET
}; // enum MilestoneType

enum PrintSettings
{
    MINIMAL,
    STANDARD,
    VERBOSE
}; // enum PrintSettings

struct TransitionMilestone
{
    MilestoneType type;
    MilestoneStatus status;
    double temperature;

    double alpha;
    double betaH;
    double H;
    double we;
    double cs_plus;
    double cs_minus;
    double n;
    double Rs;
    double Rbar;
    double dt;

private:
    PrintSettings print_setting = PrintSettings::STANDARD;

public:
    TransitionMilestone() = default;
    TransitionMilestone(const MilestoneType& type_in)
    : type(type_in), status(MilestoneStatus::ERR), temperature(0.0) {}

    void set_print_setting(PrintSettings setting) {
        print_setting = setting;
    }

    PrintSettings get_print_setting() const {
        return print_setting;
    }

    const std::string format_double(double value, int precision = 6) const {
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(precision) << value;
        return oss.str();
    }

    const std::string format_status_string() const {
        switch (status)
        {
            case MilestoneStatus::YES:
                return "YES";
            case MilestoneStatus::FAST:
                return "FAST";
            case MilestoneStatus::NO:
                return "NO";
            default:
                return "ERR";
        }
    }

    const std::string format_milestone_string() const
    {
        std::string output = "  status = " + format_status_string() + "\n";
        output += "  temperature = " + std::to_string(temperature) + " GeV\n";
        
        if (print_setting == PrintSettings::MINIMAL) {
            return output;
        }

        if (print_setting == PrintSettings::STANDARD || print_setting == PrintSettings::VERBOSE) {
            output += "  alpha = " + std::to_string(alpha) + "\n";
            output += "  betaH = " + std::to_string(betaH) + "\n";
            output += "  H = " + format_double(H) + "\n";
        }

        if (print_setting == PrintSettings::VERBOSE) {
            output += "  we = " + std::to_string(we) + "\n";
            output += "  cs_plus = " + std::to_string(cs_plus) + "\n";
            output += "  cs_minus = " + std::to_string(cs_minus) + "\n";
            output += "  Rs = " + std::to_string(Rs) + "\n";
            output += "  Rbar = " + std::to_string(Rbar) + "\n";
            output += "  dt = " + std::to_string(dt) + "\n";
        }
        
        return output;
    }

    friend std::ostream &operator<<(std::ostream& o, const TransitionMilestone &milestone) 
    {
        o << milestone.format_milestone_string();
        return o;
    }

}; // struct TransitionMilestone

struct RadiiDistribution 
{
    double ref_temperature;

    std::vector<double> temperature_values;
    std::vector<double> radii_values;
    std::vector<double> dndR_values;
    std::vector<double> log_dndR_values;

    double peak_radius;
    double peak_nuc_temperature;

    RadiiDistribution(
        const double& ref_temperature_in, 
        const std::vector<double>& temperature_values_in, 
        const std::vector<double>& radii_values_in, 
        const std::vector<double>& dndR_values_in,
        const std::vector<double>& log_dndR_values_in) :
    ref_temperature(ref_temperature_in), temperature_values(temperature_values_in), radii_values(radii_values_in), dndR_values(dndR_values_in), log_dndR_values(log_dndR_values_in)
    {
        alglib::real_1d_array t_array, r_array, log_dndR_array;

        t_array.setcontent(temperature_values.size(), temperature_values.data());
        r_array.setcontent(radii_values.size(), radii_values.data());
        log_dndR_array.setcontent(log_dndR_values.size(), log_dndR_values.data());

        alglib::spline1dbuildcubic(r_array, t_array, temperature_spline);
        alglib::spline1dbuildcubic(r_array, log_dndR_array, log_dndR_spline);

    }

    const double 
    get_nucleation_temperature(const double& radius)
    {
        double temperature = alglib::spline1dcalc(temperature_spline, radius);
        return temperature;
    }

    const double 
    get_dndR(const double& radius)
    {
        double log_dndR = alglib::spline1dcalc(log_dndR_spline, radius);
        return exp(log_dndR);
    }

private:

    /** Spline to extract nucleation temp from given radius */
    alglib::spline1dinterpolant temperature_spline;

    /** Spline to dndR from given radius */
    alglib::spline1dinterpolant log_dndR_spline;

    /** Extracts peak radius and temperature of dndR curve */
    // double find_peak_radius();

}; // struct RadiiDistribution

class TransitionMetrics 
{

    const FalseVacuumDecayRate& decay_rate;

    const EquationOfState& eos;

    double t_min, t_max;

    alglib::spline1dinterpolant a2a1_integrand_spline;

    alglib::spline1dinterpolant log_Vext_spline;

    PROPERTY(double, total_number_temp_steps, 200);

    PROPERTY(bool, use_pf_in_nt_integrand, true);

    PROPERTY(bool, use_bag_dtdT, false);

    PROPERTY(double, vw, 0.577);

    PROPERTY(double, dof, 106.75);

    PROPERTY(double, newtonG, 1/((1.22 * 1e19)*(1.22 * 1e19)));

    PROPERTY(double, percolation_target, 0.71);

    PROPERTY(double, completion_target, 1e-6);

    PROPERTY(double, onset_target, 1 - 1e-6);

    PROPERTY(double, nucleation_target, 1.00);

    PROPERTY(double, temperature_abs_tol, 1e-8);

public :

    TransitionMilestone onset_milestone;
    TransitionMilestone percolation_milestone;
    TransitionMilestone completion_milestone;
    TransitionMilestone nucleation_milestone;

    TransitionMetrics(const FalseVacuumDecayRate& decay_rate_in, const EquationOfState& eos_in) :
    decay_rate(decay_rate_in), eos(eos_in), t_min(decay_rate_in.get_t_min()), t_max(decay_rate_in.get_t_max()) 
    {
        auto start_scale_factor_calculation = std::chrono::high_resolution_clock::now();
        make_scale_factor_ratio_integrand_spline();
        auto end_scale_factor_calculation = std::chrono::high_resolution_clock::now();
        auto duration_scale_factor_calculation = std::chrono::duration_cast<std::chrono::milliseconds>(end_scale_factor_calculation - start_scale_factor_calculation);
        LOG(info) << "scale_factor_ratio_integrand spline created in " << duration_scale_factor_calculation.count() << " ms" << std::endl;

        auto start_log_extended_volume_calculation = std::chrono::high_resolution_clock::now();
        compute_log_extended_volume_spline();
        auto end_log_extended_volume_calculation = std::chrono::high_resolution_clock::now();
        auto duration_log_extended_volume_calculation = std::chrono::duration_cast<std::chrono::milliseconds>(end_log_extended_volume_calculation - start_log_extended_volume_calculation);
        LOG(info) << "log_extended_volume spline created in " << duration_log_extended_volume_calculation.count() << " ms" << std::endl;
    }

    void compute_milestones() {
        onset_milestone = get_transition_milestone(MilestoneType::ONSET);
        percolation_milestone = get_transition_milestone(MilestoneType::PERCOLATION);
        completion_milestone = get_transition_milestone(MilestoneType::COMPLETION);
        nucleation_milestone = get_transition_milestone(MilestoneType::NUCLEATION);
    }

    const double get_hubble_rate(const double& T);

    const double get_dtdT(const double& T);

    const double get_atop_abottom(const double& Ttop, const double& Tbottom);

    const double get_extended_volume(const double& T);

    const double get_extended_volume_from_spline(const double& T);

    const double get_false_vacuum_fraction(const double& T);

    const double get_nucleation_rate(const double& T);

    const double get_bubble_density(const double& T);

    const double get_bubble_radius_integral(const double& T);

    const double get_duration(const double& T);

    const TransitionMilestone get_transition_milestone(const MilestoneType type);

    const RadiiDistribution get_radii_distribution(const double& temperature);

private:

    const double find_temperature(std::function<double(double)> target_function, double tol = 1e-8, boost::uintmax_t max_iter = 100);

    const bool valid_lower_bound(std::function<double(double)> target_function, double tol = 1e-8)
    {
        return target_function(t_min) < tol;
    }

    std::function<double(double)> get_target_function(const MilestoneType type);

    void make_scale_factor_ratio_integrand_spline();

    const double get_volume_term(const double& T1, const double& T2);

    const double extended_volume_integrand(const double& T1, const double& T2);

    const double bubble_radius_integrand(const double& T1, const double& T2);

    void compute_log_extended_volume_spline();
};


} // namespace PhaseTracer

#endif // PHASETRACER_TRANSITION_METRICS_HPP_