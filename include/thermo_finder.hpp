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

#ifndef PHASETRACER_THERMO_FINDER_HPP_
#define PHASETRACER_THERMO_FINDER_HPP_

#include <cmath>
#include <vector>
#include <optional>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <memory>
#include <interpolation.h>

#include "property.hpp"
#include "scale.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "action_calculator.hpp"
#include "false_vacuum_decay_rate.hpp"
#include "equation_of_state.hpp"
#include "friedmann_evolution.hpp"

namespace PhaseTracer {

/**
 * @struct ThermalProfiles
 * @brief Structure representing the thermal profiles of a phase transition.
 * 
 * This structure contains vectors representing various thermal quantities as functions of temperature or time, 
 * including the Hubble rate, bounce action, false vacuum decay rate, nucleation rate, and bubble properties.
 */
struct ThermalProfiles
{
    /** @brief Vector of temperatures. */
    std::vector<double> temperature;
    /** @brief Vector of the time-temperature relation dtdT */
    std::vector<double> dtdT;
    /** @brief Vector of times corresponding to the temperatures. */
    std::vector<double> time;
    /** @brief Vector of Hubble rates corresponding to the temperatures. */
    std::vector<double> hubble_rate;
    /** @brief Vector of bounce action corresponding to the temperatures. */
    std::vector<double> bounce_action;
    /** @brief Vector of extended volume corresponding to the temperatures. */
    std::vector<double> extended_volume;
    /** @brief Vector of false vacuum decay rate corresponding to the temperatures. */
    std::vector<double> false_vacuum_decay_rate;
    /** @brief Vector of false vacuum fraction corresponding to the temperatures. */
    std::vector<double> false_vacuum_fraction;
    /** @brief Vector of the derivative of the false vacuum fraction corresponding to the temperatures. */
    std::vector<double> d_false_vacuum_fraction;
    /** @brief Vector of nucleation rate corresponding to the temperatures. */
    std::vector<double> nucleation_rate;
    /** @brief Vector of mean bubble separation corresponding to the temperatures. */
    std::vector<double> mean_bubble_separation;
    /** @brief Vector of mean bubble radius corresponding to the temperatures. */
    std::vector<double> mean_bubble_radius;

    ThermalProfiles() = default;

    /** @brief Writes the thermal profiles to a CSV file.
     *  @param filename The name of the CSV file to write to.
     */
    void 
    write(const std::string& filename) const
    {
        std::ofstream file(filename);
        file << "# T,dtdT,t,H,S3/T,Gamma,Vext,Pf,dPf,Nt,RsH,RbarH\n";
        for (size_t i = 0; i < temperature.size(); ++i) {
            file << std::scientific << std::setprecision(10)
                 << temperature[i] << ","
                 << dtdT[i] << ","
                 << time[i] << ","
                 << hubble_rate[i] << ","
                 << bounce_action[i] << ","
                 << false_vacuum_decay_rate[i] << ","
                 << extended_volume[i] << ","
                 << false_vacuum_fraction[i] << ","
                 << d_false_vacuum_fraction[i] << ","
                 << nucleation_rate[i] << ","
                 << mean_bubble_separation[i] << ","
                 << mean_bubble_radius[i] << "\n";
        }
        file.close();
    }
};

/**
 * @brief Structure representing a set of thermal parameters for a phase transition.
 *
 * This structure contains references to the action calculator, transition, false vacuum decay rate,
 * equation of state, and Friedmann evolution. It also stores the critical temperature, nucleation history,
 * transition milestones, and thermal profiles.
 */
struct ThermalParameterSet
{
    /** @brief Critical temperature of the phase transition. */
    double TC;

    /** @brief Nucleation history of the phase transition. */
    NucleationHistory nucleation_history;

    /** @brief Onset milestone of the phase transition. */
    TransitionMilestone onset;
    /** @brief Percolation milestone of the phase transition. */
    TransitionMilestone percolation;
    /** @brief Completion milestone of the phase transition. */
    TransitionMilestone completion;
    /** @brief Nucleation milestone of the phase transition. */
    TransitionMilestone nucleation;

    /** @brief Thermal profiles of the phase transition. */
    ThermalProfiles profiles;

    ThermalParameterSet
    (
        const Transition& t_in, 
        const ActionCalculator& ac_in,
        // FalseVacuumDecayRate settings
        int action_spline_evaluations_in = 50,
        int warm_start_chunk_size_in = 0,
        // EquationOfState settings
        int eos_spline_evaluations_in = 100,
        double eos_background_dof_in = 0.0,
        // FriedmannEvolution settings
        double percolation_target_in = 0.71,
        double completion_target_in = 1e-8,
        double onset_target_in = 1 - 1e-8,
        double nucleation_target_in = 1.00,
        bool use_bag_dtdT_in = false,
        double temperature_abs_tol_in = 1e-8,
        // CustomPrefactor
        FalseVacuumDecayRate::PrefactorFunction prefactor_in = {}
    );

    /** @brief Returns by reference the false vacuum decay rate associated with this thermal parameter set */
    FalseVacuumDecayRate& get_decay_rate() const { return *decay_rate; }

    /** @brief Returns by reference the equation of state associated with this thermal parameter set */
    EquationOfState& get_equation_of_state() const { return *eos; }

    /** @brief Returns by reference the Friedmann evolution associated with this thermal parameter set */
    FriedmannEvolution& get_friedmann_evolution() const { return *friedmann_evolution; }
    
    /** @brief Returns the critical temperature associated with this thermal parameter set */
    double get_TC() const { return TC; }

    /** @brief Pretty-print the thermal parameter set to an output stream.
     *  @param o The output stream to which the thermal parameter set will be printed.
     *  @param tps The thermal parameter set to be printed.
     *  @return The output stream with the printed thermal parameter set.
     */
    friend std::ostream &operator<<(std::ostream& o, const ThermalParameterSet &tps) 
    {
        o << "=== transition @ TC = " << tps.TC << " ===" << "\n";
        o << "MILESTONE : ONSET" << "\n";
        o << tps.onset;
        o << "MILESTONE : PERCOLATION" << "\n";
        o << tps.percolation;
        o << "MILESTONE : NUCLEATION" << "\n";
        o << tps.nucleation;
        o << "MILESTONE : COMPLETION" << "\n";
        o << tps.completion;
        return o;
    }

private:
    /** @brief Reference to the action calculator */
    const ActionCalculator& ac;
    
    /** @brief Pointer to the transition associated with this thermal parameter set */
    std::unique_ptr<Transition> transition;

    /** @brief Pointer to the false vacuum decay rate associated with this thermal parameter set */
    std::unique_ptr<FalseVacuumDecayRate> decay_rate;

    /** @brief Pointer to the equation of state associated with this thermal parameter set */
    std::unique_ptr<EquationOfState> eos;

    /** @brief Pointer to the Friedmann evolution associated with this thermal parameter set */
    std::unique_ptr<FriedmannEvolution> friedmann_evolution;
};

/**
 * @enum ValidateMethod
 * @brief Enumeration for the validation method used in the ThermoFinder class.
 * 
 * TEMP: Validate using temperature.
 * VEV: Validate using vacuum expectation value.
 * NONE: No validation.
 */
enum ValidateMethod
{
    TEMP,
    VEV,
    NONE
};

class ThermoFinder 
{
public :

    /** @brief Pretty-printer for all ThermalParameterSets in this object */
    friend std::ostream &operator<<(std::ostream &o, const ThermoFinder &a);

    ThermoFinder(const ActionCalculator& ac_in) : ac(ac_in) {
        LOG(warning) << "ThermoFinder initialised without TransitionFinder.";
    };

    ThermoFinder(const std::optional<TransitionFinder>& tf_in, const ActionCalculator& ac_in) : tf(tf_in), ac(ac_in) {};

    /** @brief Sets the prefactor function for the false vacuum decay rate. */
    void set_prefactor_function(FalseVacuumDecayRate::PrefactorFunction f) { prefactor_function = f; }

    /** @brief Indicates whether all transitions have already been calculated. */
    bool calculated_thermal_parameters = false;

    /** @brief Container for all transitions between any two phases. */
    std::vector<ThermalParameterSet> thermal_parameters;

    /** 
     * @brief Calculates (once) and returns every transition's thermal parameters.
     *  @return Reference to the vector of ThermalParameterSet objects.
     *  @note ThermalParameterSet is move-only, so the vector cannot be copied out.
     */
    const std::vector<ThermalParameterSet>& get_thermal_parameters();

    /** @brief Finds all thermal parameters for the transitions. */
    void find_thermal_parameters();

    /** 
     * @brief Retrieves the thermal parameter set for a specific transition.
     *  @param t The transition for which to retrieve the thermal parameter set.
     *  @return The corresponding ThermalParameterSet object.
     */
    ThermalParameterSet get_thermal_parameter_set(Transition t);

    /** 
     * @brief Computes the thermal profiles for the phase transition.
     *  @param decay_rate The false vacuum decay rate object.
     *  @param tm The Friedmann evolution object.
     *  @return The computed ThermalProfiles object.
     */
    ThermalProfiles compute_thermal_profiles(const FalseVacuumDecayRate& decay_rate, FriedmannEvolution& tm);

    /** 
     * @brief Adds thermal parameter values to a transition milestone.
     *  @param milestone The transition milestone to update.
     *  @param decay_rate The false vacuum decay rate object.
     *  @param eos The equation of state object.
     *  @param tm The Friedmann evolution object.
     */
    const void add_thermal_parameter_values
    (
        TransitionMilestone& milestone, 
        const FalseVacuumDecayRate& decay_rate, 
        const EquationOfState& eos, 
        FriedmannEvolution& tm
    );

    /** 
     * @brief Adds the reheating temperature to a transition milestone.
     *  @param milestone The transition milestone to update.
     *  @param tm The Friedmann evolution object.
     */
    const void add_reheating_temperature(
        TransitionMilestone& milestone, 
        FriedmannEvolution& tm
    );

    /** 
     * @brief Fills the nucleation history with relevant milestones and decay rates.
     *  @param history The nucleation history object to update.
     *  @param percolation The percolation transition milestone.
     *  @param nucleation The nucleation transition milestone.
     *  @param decay_rate The false vacuum decay rate object.
     *  @param tm The Friedmann evolution object.
     */
    void fill_nucleation_history
    (
        NucleationHistory& history, 
        TransitionMilestone& percolation, 
        TransitionMilestone& nucleation, 
        const FalseVacuumDecayRate& decay_rate, 
        FriedmannEvolution& tm
    );

    /** 
     * @brief Computes the alpha parameter at a given temperature.
     *  @param temperature The temperature at which to compute alpha.
     *  @param eos The equation of state object.
     *  @param use_munu Flag indicating whether to use mu/nu corrections.
     *  @return The computed alpha value.
     */
    const double get_alpha(const double& temperature, const EquationOfState& eos, bool use_munu = false);

    /** 
     * @brief Computes the beta/H parameter at a given temperature.
     *  @param temperature The temperature at which to compute beta/H.
     *  @param decay_rate The false vacuum decay rate object.
     *  @return The computed beta/H value.
     * 
     * This function computes the beta/H parameter based on the false vacuum decay rate at the specified temperature.
     * Note this uses the old approximation, assuming radiation domination. More reliable estimates can be obtained
     * using beta_eff, beta_1, or beta_2.
     */
    const double get_betaH(const double& temperature, const FalseVacuumDecayRate& decay_rate);

    /** 
     * @brief Computes the effective beta/H parameter.
     *  @param vw The bubble wall velocity.
     *  @param RsH The ratio of the bubble radius to the Hubble radius.
     *  @return The computed effective beta/H value.
     */
    const double get_betaH_eff(const double& vw, const double& RsH);

    /** 
     * @brief Computes the beta_1/H at a given temperature.
     *  @param temperature The temperature at which to compute beta_1/H.
     *  @param decay_rate The false vacuum decay rate object.
     *  @param tm The Friedmann evolution object.
     *  @return The computed beta_1H value.
     */
    const double get_betaH_1(const double& temperature, const FalseVacuumDecayRate& decay_rate, FriedmannEvolution& tm);

    /** 
     * @brief Computes the beta_2/H at a given temperature.
     *  @param temperature The temperature at which to compute beta_2/H.
     *  @param decay_rate The false vacuum decay rate object.
     *  @param tm The Friedmann evolution object.
     *  @return The computed beta_2/H value.
     */
    const double get_betaH_2(const double& temperature, const FalseVacuumDecayRate& decay_rate, FriedmannEvolution& tm);

    /** 
     * @brief Computes the Hubble parameter at a given temperature.
     *  @param temperature The temperature at which to compute the Hubble parameter.
     *  @param tm The Friedmann evolution object.
     *  @return The computed Hubble parameter.
     */
    const double get_H(const double& temperature, FriedmannEvolution& tm);

    /** 
     * @brief Computes the enthalpy/energy ratio at a given temperature.
     *  @param temperature The temperature at which to compute w_e.
     *  @param eos The equation of state object.
     *  @return The computed w_e value.
     */
    const double get_we(const double& temperature, const EquationOfState& eos);

    /** 
     * @brief Computes the speed of sound squared at a given temperature.
     *  @param temperature The temperature at which to compute the speed of sound squared.
     *  @param eos The equation of state object.
     *  @return A pair containing the adiabatic and isothermal speed of sound squared.
     */
    const std::pair<double, double> get_cs(const double& temperature, const EquationOfState& eos);

    /** 
     * @brief Computes the number density at a given temperature.
     *  @param temperature The temperature at which to compute the number density.
     *  @param tm The Friedmann evolution object.
     *  @return The computed number density.
     */
    const double get_n(const double& temperature, FriedmannEvolution& tm);

    /** 
     * @brief Computes the average bubble radius at a given temperature.
     *  @param temperature The temperature at which to compute the average bubble radius.
     *  @param tm The Friedmann evolution object.
     *  @return The computed average bubble radius.
     */
    const double get_Rbar(const double& temperature, FriedmannEvolution& tm);

    /** 
     * @brief Computes the time interval at a given temperature.
     *  @param temperature The temperature at which to compute the time interval.
     *  @param tm The Friedmann evolution object.
     *  @return The computed time interval.
     */
    const double get_dt(const double& temperature, FriedmannEvolution& tm);

    /** 
     * @brief Wrapper function to compute the percolation temperature.
     *  @param vw The bubble wall velocity.
     *  @param percolation_target The target percolation fraction.
     *  @param tm The Friedmann evolution object.
     *  @return The computed percolation temperature.
     */
    const double get_percolation_temperature_wrapper(const double& vw, const double& percolation_target, const FriedmannEvolution& tm);

    /** 
     * @brief Wrapper function to compute the bubble wall velocity.
     *  @param temperature The temperature at which to compute the bubble wall velocity.
     *  @param tm The Friedmann evolution object.
     *  @param eos The equation of state object.
     *  @return The computed bubble wall velocity.
     */
    const double get_vw_wrapper(const double& temperature, const FriedmannEvolution& tm, const EquationOfState& eos);

    /** 
     * @brief Revises the percolation temperature based on the given percolation milestone.
     *  @param percolation The percolation milestone to be revised.
     *  @param eos The equation of state object.
     *  @param tm The Friedmann evolution object.
     */
    const void revise_percolation_temperature(TransitionMilestone& percolation, const EquationOfState& eos, const FriedmannEvolution& tm);

private : 

    /**
     * @brief Optional TransitionFinder object used for finding phase transitions.
     * 
     * This object is used internally by the ThermoFinder to locate phase transitions.
     */
    std::optional<TransitionFinder> tf;

    /** @brief Reference to the ActionCalculator object used for computing actions. */
    const ActionCalculator& ac;

    // ======================= Settings for ThermoFinder =======================

    /** @brief Print settings for the onset of the phase transition. */
    PROPERTY(PrintSettings, onset_print_setting, PrintSettings::MINIMAL);

    /** @brief Print settings for the percolation of the phase transition. */
    PROPERTY(PrintSettings, percolation_print_setting, PrintSettings::STANDARD);

    /** @brief Print settings for the nucleation of the phase transition. */
    PROPERTY(PrintSettings, nucleation_print_setting, PrintSettings::STANDARD);

    /** @brief Print settings for the completion of the phase transition. */
    PROPERTY(PrintSettings, completion_print_setting, PrintSettings::MINIMAL);

    /** @brief Flag indicating whether to update the percolation temperature. */
    PROPERTY(bool, update_percolation_temperature, false);

    /** @brief Flag indicating whether to compute profiles for the phase transition. */
    PROPERTY(bool, compute_profiles, false);

    /** @brief Number of temperature evluations for profiles. */
    PROPERTY(double, n_temp_profiles, 250);

    /** @brief Bubble wall velocity. 
     * @note This is overwritten by vw at the percolation temp if updated.
    */
    PROPERTY(double, vw, 1/sqrt(3.0));

    /** @brief Threshold for the temperature. */
    PROPERTY(double, temperature_threshold, 1);

    /** @brief Threshold for the vacuum expectation value (VEV). */
    PROPERTY(double, vev_threshold, 1);

    /** @brief Sets the validation method for screening transitions */
    PROPERTY(ValidateMethod, default_validation_method, NONE);

    /** @brief A custom function for validating transitions. */
    PROPERTY(
        std::function<std::vector<Transition>(const std::vector<Transition>&)>,
        transition_filter, {}
    ); // TODO define a custom signature

    // =================== Settings for FalseVacuumDecayRate ===================

    /** @brief Number of spline evaluations for false vacuum decay rate. */
    PROPERTY(int, action_spline_evaluations, 50);

    /** @brief Number of action evaluations solved as one block. */
    PROPERTY(int, warm_start_chunk_size, 0)

    // unused: t_min, t_max

    // ===================== Settings for EquationOfState =====================

    /** @brief Number of spline evaluations for the equation of state. */
    PROPERTY(int, eos_spline_evaluations, 250);

    /** @brief Number of background degrees of freedom in EOS. */
    PROPERTY(double, eos_background_dof, 0.0);

    // unused: t_min, t_max, energy_norm

    // ==================== Settings for FriedmannEvolution ====================

    /** @brief Target value for the percolation of the phase transition. */
    PROPERTY(double, percolation_target, 0.71);

    /** @brief Target value for the completion of the phase transition. */
    PROPERTY(double, completion_target, 1e-6);

    /** @brief Target value for the onset of the phase transition. */
    PROPERTY(double, onset_target, 1 - 1e-6);

    /** @brief Target value for the nucleation of the phase transition. */
    PROPERTY(double, nucleation_target, 1.00);

    /** @brief Flag indicating whether to use the bag model time-temperature relationship */
    PROPERTY(bool, use_bag_dtdT, false);

    /** @brief Absolute tolerance for the temperature. */
    PROPERTY(double, temperature_abs_tol, 1e-6);

    // unused: volume_term_integration_steps, vw, M_planck, newtonG


    std::vector<Transition> default_transition_filter(const std::vector<Transition>& transitions);

    /** Optional custom decay-rate prefactor (e.g. a BubbleDetPrefactor). When
     *  empty the FalseVacuumDecayRate uses its default analytic prefactor. */
    FalseVacuumDecayRate::PrefactorFunction prefactor_function;

};

} // namespace PhaseTracer

#endif // PHASETRACER_THERMO_FINDER_HPP_