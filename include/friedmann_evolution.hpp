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