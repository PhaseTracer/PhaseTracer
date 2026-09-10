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
#include "transition_milestones.hpp"
#include "bubble_statistics.hpp"

namespace PhaseTracer {

/**
 * @struct FriedmannSystem
 * @namespace PhaseTracer
 * @brief Represents the state of the Friedmann system at different times.
 */
struct FriedmannSystem
{
    /** @brief Logarithm of time values. */
    std::vector<double> log_time;
    /** @brief Time values. */
    std::vector<double> time;

    /** @brief Energy density in the true vacuum. */
    std::vector<double> e_t;
    /** @brief Energy density in the false vacuum. */
    std::vector<double> p_t;
    /** @brief Enthalpy density in the true vacuum. */
    std::vector<double> w_t;
    /** @brief Entropy density in the true vacuum. */
    std::vector<double> s_t;
    /** @brief Energy density in the false vacuum. */
    std::vector<double> e_f;
    /** @brief Pressure in the false vacuum. */
    std::vector<double> p_f;
    /** @brief Enthalpy density in the false vacuum. */
    std::vector<double> w_f;
    /** @brief Entropy density in the false vacuum. */
    std::vector<double> s_f;

    /** @brief Temperature in the false vacuum. */
    std::vector<double> T_f;
    /** @brief Temperature in the true vacuum. */
    std::vector<double> T_t;

    /** @brief Hubble parameter. */
    std::vector<double> hubble;
    /** @brief Scale factor. */
    std::vector<double> a;
    /** @brief False vacuum decay rate. */
    std::vector<double> gamma;
    /** @brief Bounce action. */
    std::vector<double> action;

    /** @brief Integrals I_n used in the evolution equations. */
    std::vector<double> I_0;
    std::vector<double> I_1;
    std::vector<double> I_2;
    std::vector<double> I_3;

    /** @brief Nucleation rate. */
    std::vector<double> nucleation_rate;
    /** @brief Mean bubble number density. */
    std::vector<double> number_density;
    /** @brief Mean bubble radius. */
    std::vector<double> mean_bubble_radius;

    /** @brief Writes the Friedmann system data to a CSV file. 
     * @param filename The name of the output CSV file.
     * @param vw The wall velocity, used in the calculation of the true vacuum fraction (default is 1/sqrt(3)).
    */
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

/**
 * @class FriedmannEvolution
 * @namespace PhaseTracer
 * @brief Solves the coupled Friedmann and JMAK equations for a cosmological phase transition.
 */
class FriedmannEvolution
{
public :

    /** @brief Early exit condition if true vacuum fraction at t_min is zero. */
    bool early_exit = false;

    /** @brief The Friedmann system. */
    FriedmannSystem system;

    /** @brief Onset TransitionMilestone. */
    TransitionMilestone onset_milestone;
    /** @brief Percolation TransitionMilestone. */
    TransitionMilestone percolation_milestone;
    /** @brief Completion TransitionMilestone. */
    TransitionMilestone completion_milestone;
    /** @brief Nucleation TransitionMilestone. */
    TransitionMilestone nucleation_milestone;
    /** @brief Nucleation history. */
    NucleationHistory nucleation_history;

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

    /** @brief Computes all transition milestones. */
    void compute_milestones()
    {
        onset_milestone = get_transition_milestone(MilestoneType::ONSET);
        percolation_milestone = get_transition_milestone(MilestoneType::PERCOLATION);
        completion_milestone = get_transition_milestone(MilestoneType::COMPLETION);
        nucleation_milestone = get_transition_milestone(MilestoneType::NUCLEATION);
    }

    /** 
     * @brief Computes the nucleation history over the specified temperature range.
     * @param t_min The minimum temperature for the nucleation history.
     * @param t_max The maximum temperature for the nucleation history.
     *
    */
    void compute_nucleation_history(const double& t_min, const double& t_max);

    /**
     * @brief Computes the nucleation history over the temperature range defined by t_min and t_max.
     * @param T_false The false vacuum temperature at which to compute the hubble rate.
     * @return The Hubble rate at the specified false vacuum temperature.
     * @note Does not check T_false against the valid temperature range.
     * 
     * This function returns the Hubble rate at the specified false vacuum temperature. If the system has not yet 
     * been solved, it falls back to assuming the false vacuum fraction is one. Otherwise, it uses
     * the Hubble rate spline obtained from the solved Friedmann system.
      */
    const double get_hubble_rate(const double& T_false) const;

    /**
     * @brief Computes the false vacuum time-temperature relation dt/dT.
     * @param T_false The false vacuum temperature at which to compute the time-temperature relation.
     * @return The derivative dt/dT at the specified false vacuum temperature.
     */
    const double get_time_temperature_false(const double& T_false) const;

    /**
     * @brief Computes the scale factor a(T_false).
     * @param T_false The false vacuum temperature at which to compute the scale factor.
     * @return The scale factor at the specified false vacuum temperature.
     * @note Requires the Friedmann system to be solved and thescale factor spline to be available.
     */
    const double get_scale_factor(const double& T_false) const;

    /**
     * @brief Computes the ratio of scale factors a(Ttop)/a(Tbottom).
     * @param Ttop The false vacuum temperature in the numerator.
     * @param Tbottom The false vacuum temperature in the denominator.
     * @return The ratio of scale factors at the specified false vacuum temperatures.
     * @note Requires the Friedmann system to be solved and the scale factor spline to be available.
     */
    const double get_scale_factor_ratio(const double& Ttop, const double& Tbottom) const;

    /**
     * @brief Computes the false vacuum fraction at the specified false vacuum temperature.
     * @param T_false The false vacuum temperature at which to compute the false vacuum fraction.
     * @return The false vacuum fraction at the specified false vacuum temperature.
     * @note Requires the Friedmann system to be solved and the false vacuum fraction spline to be available.
     */
    const double get_false_vacuum_fraction(const double& T_false) const;

    /**
     * @brief Computes the nucleation rate at the specified false vacuum temperature.
     * @param T_false The false vacuum temperature at which to compute the false vacuum fraction.
     * @return The nucleation rate at the specified false vacuum temperature.
     * @note Requires the Friedmann system to be solved and the nucleation rate spline to be available.
     */
    const double get_nucleation_rate(const double& T_false) const;

    /**
     * @brief Computes the average bubble density at the specified false vacuum temperature.
     * @param T_false The false vacuum temperature at which to compute the false vacuum fraction.
     * @return The average bubble density at the specified false vacuum temperature.
     * @note Requires the Friedmann system to be solved and the bubble density spline to be available.
     */
    const double get_bubble_density(const double& T_false) const;

    /**
     * @brief Computes the average bubble radius at the specified false vacuum temperature.
     * @param T_false The false vacuum temperature at which to compute the false vacuum fraction.
     * @return The average bubble radius at the specified false vacuum temperature.
     * @note Requires the Friedmann system to be solved and the bubble radius spline to be available.
     */
    const double get_mean_bubble_radius(const double& T_false) const;

    /**
     * @brief Computes the time at the specified false vacuum temperature.
     * @param T_false The false vacuum temperature at which to compute the false vacuum fraction.
     * @return The time at the specified false vacuum temperature.
     * 
     * This function returns the comoving time, t, at a given false vacuum temperature. This is obtained from the 
     * Friedmann system if computed. Otherwise, it integrates the false vacuum time-temperature relationship
     * to gain an estimate.
     */
    const double get_t(const double& T) const;

    /**
     * @brief Computes the true vacuum temperature at the specified false vacuum temperature.
     * @param T_false The false vacuum temperature at which to compute the false vacuum fraction.
     * @return The true vacuum temperature at the specified false vacuum temperature.
     * @note Requires the Friedmann system to be solved and the true vacuum temperature spline to be available.
     */
    const double get_T_true(const double& T_true) const;

    /**
     * @brief Performs a second-order expansion of the bounce action, S(t).
     * @param T_false The false vacuum temperature at which to compute the false vacuum fraction.
     * @return A pair containing the coefficients (beta_1, beta_2) of the Taylor expansion.
     */
    const std::pair<double, double> get_action_expansion(const double& T_false) const;

    /** @brief Returns the final temperature of the Friedmann evolution */
    const double get_t_min() const { return t_min; }
    
    /** @brief Returns the initial temperature of the Friedmann evolution */
    const double get_t_max() const { return t_max; }

    /** @brief Determines the transition milestone for a given type.
     *  @param type The type of milestone to determine.
     *  @return The transition milestone corresponding to the specified type.
     *  @note This is public because it is called by a percolation temperature wrapper in ThermoFinder.
     */
    const TransitionMilestone get_transition_milestone(const MilestoneType type);

    /** @brief Returns the distribution of bubble radii at a given temperature.
     *  @param temperature The temperature at which to compute the radii distribution.
     *  @return The distribution of bubble radii at the specified temperature.
     *  @note WIP: This function is still under development and may not provide accurate results.
     */
    const RadiiDistribution get_radii_distribution(const double& temperature);

    /** @brief Returns the distribution of bubble lifetimes for a given timescale.
     *  @param timescale The characteristic timescale for the lifetime distribution, should be invariant upon rescaling.
     *  @param lifetime_min_fraction The minimum fraction of the timescale to consider for the lifetime distribution.
     *  @return The distribution of bubble lifetimes corresponding to the specified timescale.
     *  @note WIP: This function is still under development and may not provide accurate results.
     */
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

    /** @brief Spline for T_true(T_false) */
    mutable alglib::spline1dinterpolant reheating_spline;
    /** @brief Spline for log(t)(T_false) */
    mutable alglib::spline1dinterpolant log_time_spline;
    /** @brief Spline for the scale factor a(T_false) */
    mutable alglib::spline1dinterpolant scale_factor_spline;
    /** @brief Spline for the Hubble rate H(T_false) */
    mutable alglib::spline1dinterpolant hubble_rate_spline;
    /** @brief Spline for the action Action(T_false) */
    mutable alglib::spline1dinterpolant log_action_spline;
    /** @brief Spline for log_I_3(T_false) */
    mutable alglib::spline1dinterpolant log_I_3_spline;
    /** @brief Spline for the nucleation rate log_N(T_false) */
    mutable alglib::spline1dinterpolant log_nucleation_rate_spline;
    /** @brief Spline for the bubble number density log_n(T_false) */
    mutable alglib::spline1dinterpolant log_bubble_number_density_spline;
    /** @brief Spline for the mean bubble radius log_Rbar(T_false) */
    mutable alglib::spline1dinterpolant log_mean_bubble_radius_spline;
    /** @brief Whether the Friedmann splines have been computed */
    mutable bool friedmann_splines_computed = false;

    /** @brief Throws std::logic_error if solve() has not been called */
    void require_solved(const char* caller) const
    {
        if (!solved)
        {
            throw std::logic_error(std::string("FriedmannEvolution::") + caller + " called before solve().");
        }
    }

    /** @brief Refine the temperature bounds for the Friedmann evolution */
    void refine_temperature_bounds();

    /** @brief Get the Hubble rate for the given true vacuum fraction and energy densities */
    const double get_hubble_rate(const double& true_vacuum_fraction, const double& e_false, const double& e_true) const;

    /** @brief Get the false vacuum fraction from the given I3 */
    const double get_false_vacuum_fraction_from_I3(const double& I3) const;

    /** @brief Get the derivative of the false vacuum fraction with respect to I3 */
    const double get_d_false_vacuum_fraction_from_I3(const double& I3, const double& I3_dot) const;

    /** @brief Match the true vacuum temperature to the given true vacuum energy density */
    double match_T_true(const double& e_true, double tol = 1e-8, boost::uintmax_t max_iter = 100);

    /** @brief Match the false vacuum temperature to the given false vacuum energy density */
    double match_T_false(const double& e_false, double tol = 1e-8, boost::uintmax_t max_iter = 100);

    /** @brief Evolve the Friedmann equations */
    void evolve_friedmann();

    /** @brief Fit the Friedmann splines */
    const void fit_friedmann_splines() const;

    /** @brief Check if the lower bound is valid for the target function */
    const bool valid_lower_bound(std::function<double(double)> target_function, double tol = 1e-8)
    {
        return target_function(t_min) < tol;
    }

    /** @brief Find the temperature that satisfies the target function */
    const double find_temperature(std::function<double(double)> target_function, double tol = 1e-8, boost::uintmax_t max_iter = 100);

    /** @brief Get the target function for a given milestone type */
    std::function<double(double)> get_target_function(const MilestoneType type);

    /** @brief Perform Simpson integration for the given integrand over the specified range */
    double simpson_integrate(const std::function<double(double)>& integrand, const double& x_min, const double& x_max, const int& steps = 500) const;

    /** @brief Exception thrown when the transition is complete */
    struct TransitionCompleteException {};

    /** @brief Exception thrown when the false vacuum is trapped */
    struct FalseVacuumTrappingException {};

    /** @brief Exception thrown when the integration stalls */
    struct IntegrationStalledException {};

    /** @brief Number of integration steps for the volume term */
    PROPERTY(double, volume_term_integration_steps, 1000);
    /** @brief Whether to use the bag model for dT/dt */
    PROPERTY(bool, use_bag_dtdT, false);
    /** @brief Bubble wall velocity */
    PROPERTY(double, vw, 0.577);
    /** @brief Planck mass */
    PROPERTY(double, M_planck, PhaseTracer::scale() * 1.22e19);
    /** @brief Newton's gravitational constant */
    PROPERTY(double, newtonG, 1/(M_planck*M_planck));
    /** @brief Percolation target */
    PROPERTY(double, percolation_target, 0.71);
    /** @brief Completion target */
    PROPERTY(double, completion_target, 1e-6);
    /** @brief Onset target */
    PROPERTY(double, onset_target, 1 - 1e-6);
    /** @brief Nucleation target */
    PROPERTY(double, nucleation_target, 1.00);
    /** @brief Absolute tolerance for temperature */
    PROPERTY(double, temperature_abs_tol, 1e-8);

};

} // namespace PhaseTracer

#endif // PHASETRACER_FRIEDMANN_EVOLUTION_HPP_