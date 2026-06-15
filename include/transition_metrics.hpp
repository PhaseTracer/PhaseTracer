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
#include <boost/math/tools/minima.hpp>
#include <boost/numeric/odeint.hpp>

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

enum NucleationType
{
    EXPONENTIAL,
    SIMULTANEOUS
}; // enum NucleationType

struct NucleationHistory
{
    NucleationType nucleation_type;

    std::vector<double> temperature;
    std::vector<double> time;
    std::vector<double> scale_factor;
    std::vector<double> conformal_time;

    double betaH; // generic betaH, found from T * d(S/T)dT, unused?
    double betaH_1; // evaluated from d(S(t))dt
    double betaH_2; // evaluated from d^2(S(t))/dt^2

    double T_m; // temperature for peak Gamma_m
    double Gamma_m; // Peak Gamma_m, is normalised by H^4

    double T_0; // temperature for beta_1
    double Gamma_0;

    double Rs_exp; // computed from vw and beta
    double Rs_sim; // computed from Gamma_m and beta

    double betaH_2_perc;
    double betaH_2_nuc;
    double betaH_2_m;
    double betaH_1_perc;
    double betaH_1_nuc;
    double Gamma_0_perc;
    double Gamma_0_nuc;
    double Gamma_M;
    double RsH_exp_perc;
    double RsH_exp_nuc;
    double RsH_sim_perc;
    double RsH_sim_nuc;
    double RsH_sim_m;

}; // struct NucleationHistory

struct TransitionMilestone
{
    MilestoneType type;
    MilestoneStatus status;
    NucleationType nucleation_type = NucleationType::EXPONENTIAL;
    double temperature;

    double alpha;
    double alpha_munu;
    double betaH;
    double betaH_eff;
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
        if(type == MilestoneType::PERCOLATION && status == MilestoneStatus::YES)
        {
            output += "  nucleation_type = " + std::string(nucleation_type == NucleationType::EXPONENTIAL ? "exponential" : "simultaneous") + "\n";
        }

        if (print_setting == PrintSettings::MINIMAL) {
            return output;
        }

        if (print_setting == PrintSettings::STANDARD || print_setting == PrintSettings::VERBOSE) {
            output += "  alpha = " + std::to_string(alpha) + "\n";
            output += "  alpha_munu = " + std::to_string(alpha_munu) + "\n";
            output += "  betaH = " + std::to_string(betaH) + "\n";
            output += "  betaH_eff = " + std::to_string(betaH_eff) + "\n";
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

struct LifetimeDistribution
{
    std::vector<double> temperature_values;
    std::vector<double> lifetime_values;
};

struct FriedmannSystem
{
    std::vector<double> log_time;
    std::vector<double> time;

    std::vector<double> e_t;
    std::vector<double> p_t;
    std::vector<double> w_t;
    std::vector<double> s_t;
    std::vector<double> e_f;
    std::vector<double> p_f;
    std::vector<double> w_f;
    std::vector<double> s_f;

    std::vector<double> T_f;
    std::vector<double> T_t;
    std::vector<double> hubble;
    std::vector<double> a;
    std::vector<double> gamma;

    std::vector<double> I_0;
    std::vector<double> I_1;
    std::vector<double> I_2;
    std::vector<double> I_3;

    std::vector<double> nucleation_rate;
    std::vector<double> number_density;

    void write(std::string filename) const
    {
        std::ofstream out(filename);
        out << "# time,T_f,T_t,hubble,a,gamma,P_f,N,n\n";
        for (std::size_t i = 0; i < time.size(); ++i)
        {
            out << time[i]   << ","
                << T_f[i]    << ","
                << T_t[i]    << ","
                << hubble[i] << ","
                << a[i]      << ","
                << gamma[i]  << ","
                << std::exp( - 4.0 * M_PI * 0.577*0.577*0.577 / 3.0 * I_3[i]) << ","
                << nucleation_rate[i] << ","
                << number_density[i] <<
                "\n";
        }
        out.close();
    }
};

class TransitionMetrics 
{

    const FalseVacuumDecayRate& decay_rate;

    const EquationOfState& eos;

    double t_min, t_max;

    /* Friedmann splines */
    mutable alglib::spline1dinterpolant reheating_spline; // T_true(T_false)
    mutable alglib::spline1dinterpolant scale_factor_spline; // a(T_false)
    mutable alglib::spline1dinterpolant hubble_rate_spline; // H(T_false)
    mutable alglib::spline1dinterpolant log_I_3_spline; // log_I_3(T_false)
    mutable alglib::spline1dinterpolant log_nucleation_rate_spline; // log_N(T_false)
    mutable alglib::spline1dinterpolant log_bubble_number_density_spline; // log_n(T_false)
    mutable bool friedmann_splines_computed = false;

    PROPERTY(double, total_number_temp_steps, 200);

    PROPERTY(double, volume_term_integration_steps, 1000);

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

    FriedmannSystem system;

    TransitionMilestone onset_milestone;
    TransitionMilestone percolation_milestone;
    TransitionMilestone completion_milestone;
    TransitionMilestone nucleation_milestone;

    NucleationHistory nucleation_history;

    TransitionMetrics(const FalseVacuumDecayRate& decay_rate_in, const EquationOfState& eos_in) :
    decay_rate(decay_rate_in), eos(eos_in), t_min(decay_rate_in.get_t_min()), t_max(decay_rate_in.get_t_max()) 
    {
        {
            auto t0 = std::chrono::high_resolution_clock::now();
            refine_temperature_bounds();
            auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0);
            LOG(debug) << "Refined temperature bounds. Time: " << dt.count() << " ms"; 
        }
        
        {
            auto t0 = std::chrono::high_resolution_clock::now();
            evolve_friedmann();
            auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0);
            LOG(debug) << "Solved Friedmann equations. Time: " << dt.count() << " ms"; 
        }

        {
            auto t0 = std::chrono::high_resolution_clock::now();
            calculate_distributions();
            auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0);
            LOG(debug) << "Calculated distributions. Time: " << dt.count() << " ms"; 
        }

        {
            auto t0 = std::chrono::high_resolution_clock::now();
            fit_friedmann_splines();
            auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0);
            LOG(debug) << "Fit splines to Friedmann solution. Time: " << dt.count() << " ms"; 
        }

        // Write FridmannSystem arrays to reheating.csv
        system.write("example/TestThermalParameters/reheating_data/reheating.csv");
    }

    void compute_milestones() 
    {
        onset_milestone = get_transition_milestone(MilestoneType::ONSET);
        percolation_milestone = get_transition_milestone(MilestoneType::PERCOLATION);
        completion_milestone = get_transition_milestone(MilestoneType::COMPLETION);
        nucleation_milestone = get_transition_milestone(MilestoneType::NUCLEATION);
    }

    void compute_nucleation_history(const double& t_min, const double& t_max);

    const double get_hubble_rate(const double& T_false) const;

    const double get_time_temperature_false(const double& T_false) const;

    const double get_scale_factor(const double& T_false) const;

    const double get_scale_factor_ratio(const double& Ttop, const double& Tbottom) const;

    const double get_false_vacuum_fraction(const double& T_false) const;

    const double get_t(const double& T) const;

    const double get_nucleation_rate(const double& T_false) const;

    const double get_bubble_density(const double& T_false) const;

    const double get_t_min() const { return t_min; }

    const double get_t_max() const { return t_max; }

    const TransitionMilestone get_transition_milestone(const MilestoneType type);

    const RadiiDistribution get_radii_distribution(const double& temperature);

    const LifetimeDistribution get_lifetime_distribution(const double& timescale, const double& lifetime_min_fraction = 1e-6);

private:

    struct TransitionCompleteException {};

    const double find_temperature(std::function<double(double)> target_function, double tol = 1e-8, boost::uintmax_t max_iter = 100);

    const bool valid_lower_bound(std::function<double(double)> target_function, double tol = 1e-8)
    {
        return target_function(t_min) < tol;
    }

    std::function<double(double)> get_target_function(const MilestoneType type);

    const void make_scale_factor_ratio_spline() const;

    void make_volume_term_integral_spline() const;

    void refine_temperature_bounds();

    double match_T_true(const double& e_true, double tol = 1e-8, boost::uintmax_t max_iter = 100);

    double match_T_false(const double& e_false, double tol = 1e-8, boost::uintmax_t max_iter = 100);

    const double get_hubble_rate(const double& true_vacuum_fraction, const double& e_false, const double& e_true) const;

    const double get_false_vacuum_fraction_from_I3(const double& I3) const;

    void evolve_friedmann();

    const void fit_friedmann_splines() const;

    void calculate_distributions();

    std::vector<double> calculate_lifetime_distribution(const double beta, const std::vector<double>& lifetime_grid) const;

    const double get_volume_term(const double& T1, const double& T2) const;

    const double extended_volume_integrand(const double& T1, const double& T2) const;

    const double bubble_radius_integrand(const double& T1, const double& T2) const;

    void compute_log_extended_volume_spline();

    alglib::real_1d_array cumulative_simpson(const std::function<double(double)>& integrand, const alglib::real_1d_array& x, double F_initial = 0.0) const;

    std::vector<double> cumulative_simpson(const std::function<double(double)>& integrand, const std::vector<double>& x, double F_initial = 0.0) const;

    void integrate_and_fit_spline(alglib::spline1dinterpolant& spline, const std::function<double(double)>& integrand, const alglib::real_1d_array& x, double F_initial = 0.0) const;

    void integrate_and_fit_spline(alglib::spline1dinterpolant& spline, const std::function<double(double)>& integrand, int steps, double F_initial = 0.0) const;

};


} // namespace PhaseTracer

#endif // PHASETRACER_TRANSITION_METRICS_HPP_