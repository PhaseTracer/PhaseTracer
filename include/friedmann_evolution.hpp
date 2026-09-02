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

#ifndef PHASETRACER_FRIEDMANN_EVOLUTION_HPP_
#define PHASETRACER_FRIEDMANN_EVOLUTION_HPP_

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

#include "scale.hpp"
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
    double betaH_1;
    double betaH_2;
    double T_m;
}; // struct NucleationHistory

struct TransitionMilestone
{
    MilestoneType type;
    MilestoneStatus status;
    NucleationType nucleation_type = NucleationType::EXPONENTIAL;
    double temperature;
    double reheating_temperature;

    double alpha;
    double alpha_munu;
    double betaH;
    double beta1H;
    double beta2H;
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
        output += "  temperature = " + std::to_string(temperature) + " " + scale.name() + "\n";
        if(type == MilestoneType::COMPLETION && status == MilestoneStatus::YES)
        {
            output += "  reheating temperature = " + std::to_string(reheating_temperature) + " " + scale.name() + "\n";
        }
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
            output += "  beta1H = " + std::to_string(beta1H) + "\n"; 
            output += "  beta2H = " + std::to_string(beta2H) + "\n";
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
    double timescale;
    double mean_lifetime;

    std::vector<double> chi_values;
    std::vector<double> lifetime_values;
    std::vector<double> distribution_values;

    // log(I_3) and log(I_2) splines let dh/dt be evaluated analytically
    // (dh/dt = -4*pi*vw^3 * h * I_2), instead of numerically differentiating
    // a spline of h(t). I_2 and I_3 are splined in log-space (as elsewhere in
    // this class, e.g. fit_friedmann_splines) because they grow by many
    // orders of magnitude and are poorly conditioned for a direct cubic
    // spline (interpolation/extrapolation can ring/diverge). T_false_spline
    // lets gamma be evaluated as decay_rate.get_gamma(T_false(t)), using
    // decay_rate's own dedicated spline in temperature (built directly from
    // bounce-action data) rather than a second, coarser interpolation of
    // gamma vs time. Both gamma and dh/dt can be extremely sharp functions
    // of temperature/time, so avoiding an extra layer of interpolation of
    // these already-sharp quantities is important for robustness.
    alglib::spline1dinterpolant log_I3_spline;
    alglib::spline1dinterpolant log_I2_spline;
    alglib::spline1dinterpolant T_false_spline;
    alglib::spline1dinterpolant scale_factor_spline;
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
    std::vector<double> action;

    std::vector<double> I_0;
    std::vector<double> I_1;
    std::vector<double> I_2;
    std::vector<double> I_3;

    std::vector<double> nucleation_rate;
    std::vector<double> number_density;
    std::vector<double> mean_bubble_radius;

    void write(std::string filename, double vw = 1./sqrt(3)) const
    {
        std::ofstream out(filename);
        out << "# time,T_f,T_t,e_f,e_t,p_f,p_t,hubble,a,gamma,h,N,n,Rbar\n";
        for (std::size_t i = 0; i < time.size(); ++i)
        {
            out << time[i]   << ","
                << T_f[i]    << ","
                << T_t[i]    << ","
                << e_f[i]    << ","
                << e_t[i]    << ","
                << p_f[i] << ","
                << p_t[i] << ","
                << hubble[i] << ","
                << a[i]      << ","
                << gamma[i]  << ","
                << std::exp( - 4.0 * M_PI * vw*vw*vw / 3.0 * I_3[i]) << ","
                << nucleation_rate[i] << ","
                << number_density[i] << ","
                << mean_bubble_radius[i] <<
                "\n";
        }
        out.close();
    }
};

class FriedmannEvolution
{

public :

    /** Set by solve() if the h(t_min) = 1. */
    bool early_exit = false;

    FriedmannSystem system;

    TransitionMilestone onset_milestone;
    TransitionMilestone percolation_milestone;
    TransitionMilestone completion_milestone;
    TransitionMilestone nucleation_milestone;

    NucleationHistory nucleation_history;

    /**
     * @brief Binds the decay rate and equation of state; solves nothing.
     *
     * Both arguments must outlive this object.
     */
    FriedmannEvolution(FalseVacuumDecayRate& decay_rate_in, EquationOfState& eos_in) :
    decay_rate(decay_rate_in), eos(eos_in)
    {}

    /**
     * @brief Refines the temperature bounds, solves the coupled F-JMAK system, 
     * and fits the splines the accessors read.
     */
    void solve();

    /** @brief Whether solve() has completed successfully. */
    bool is_solved() const { return solved; }

    void compute_milestones()
    {
        // require_solved("compute_milestones"); // already check inside get_transition_milestone
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

    const std::pair<double, double> get_action_expansion(const double& temperature) const;

    const double get_t(const double& T) const;

    const double get_nucleation_rate(const double& T_false) const;

    const double get_bubble_density(const double& T_false) const;

    const double get_mean_bubble_radius(const double& T_false) const;

    const double get_T_true(const double& T_true) const;

    const double get_t_min() const { return t_min; }

    const double get_t_max() const { return t_max; }

    const TransitionMilestone get_transition_milestone(const MilestoneType type);

    const RadiiDistribution get_radii_distribution(const double& temperature);

    const LifetimeDistribution get_lifetime_distribution(const double& timescale, const double& lifetime_min_fraction = 1e-6);

private:

    /** Reference to the false vacuum decay rate class */
    FalseVacuumDecayRate& decay_rate;

    /** Reference to equation of state class */
    EquationOfState& eos;

    double t_min = 0.0;
    double t_max = 0.0;

    /** Set by solve() once the Friedmann system has been solved and splined */
    bool solved = false;

    /* Friedmann splines */
    mutable alglib::spline1dinterpolant reheating_spline; // T_true(T_false)
    mutable alglib::spline1dinterpolant log_time_spline; // log(t)(T_false)
    mutable alglib::spline1dinterpolant scale_factor_spline; // a(T_false)
    mutable alglib::spline1dinterpolant hubble_rate_spline; // H(T_false)
    mutable alglib::spline1dinterpolant log_action_spline; // Action(T_false)
    mutable alglib::spline1dinterpolant log_I_3_spline; // log_I_3(T_false)
    mutable alglib::spline1dinterpolant log_nucleation_rate_spline; // log_N(T_false)
    mutable alglib::spline1dinterpolant log_bubble_number_density_spline; // log_n(T_false)
    mutable alglib::spline1dinterpolant log_mean_bubble_radius_spline; // log_Rbar(T_false)
    mutable bool friedmann_splines_computed = false;

    PROPERTY(double, volume_term_integration_steps, 1000);

    PROPERTY(bool, use_bag_dtdT, false);

    PROPERTY(double, vw, 0.577);

    PROPERTY(double, M_planck, PhaseTracer::scale() * 1.22e19)

    PROPERTY(double, newtonG, 1/(M_planck*M_planck));

    PROPERTY(double, percolation_target, 0.71);

    PROPERTY(double, completion_target, 1e-6);

    PROPERTY(double, onset_target, 1 - 1e-6);

    PROPERTY(double, nucleation_target, 1.00);

    PROPERTY(double, temperature_abs_tol, 1e-8);

    struct TransitionCompleteException {};

    struct FalseVacuumTrappingException {};

    struct IntegrationStalledException {};

    /** Throws std::logic_error if solve() has not been called */
    void require_solved(const char* caller) const
    {
        if (!solved)
        {
            throw std::logic_error(std::string("FriedmannEvolution::") + caller + " called before solve().");
        }
    }

    const double find_temperature(std::function<double(double)> target_function, double tol = 1e-8, boost::uintmax_t max_iter = 100);

    const bool valid_lower_bound(std::function<double(double)> target_function, double tol = 1e-8)
    {
        return target_function(t_min) < tol;
    }

    std::function<double(double)> get_target_function(const MilestoneType type);

    void refine_temperature_bounds();

    double match_T_true(const double& e_true, double tol = 1e-8, boost::uintmax_t max_iter = 100);

    double match_T_false(const double& e_false, double tol = 1e-8, boost::uintmax_t max_iter = 100);

    const double get_hubble_rate(const double& true_vacuum_fraction, const double& e_false, const double& e_true) const;

    const double get_false_vacuum_fraction_from_I3(const double& I3) const;

    const double get_d_false_vacuum_fraction_from_I3(const double& I3, const double& I3_dot) const;

    void evolve_friedmann();

    const void fit_friedmann_splines() const;

    double simpson_integrate(const std::function<double(double)>& integrand, const double& x_min, const double& x_max, const int& steps = 500) const;

};

} // namespace PhaseTracer

#endif // PHASETRACER_FRIEDMANN_EVOLUTION_HPP_