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

#ifndef PHASETRACER_TEMPORARY_HPP_
#define PHASETRACER_TEMPORARY_HPP_

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

struct ThermalProfiles
{
    std::vector<double> temperature;
    std::vector<double> dtdT;
    std::vector<double> time;
    std::vector<double> hubble_rate;
    std::vector<double> bounce_action;
    std::vector<double> extended_volume;
    std::vector<double> false_vacuum_decay_rate;
    std::vector<double> false_vacuum_fraction;
    std::vector<double> d_false_vacuum_fraction;
    std::vector<double> nucleation_rate;
    std::vector<double> mean_bubble_separation;
    std::vector<double> mean_bubble_radius;

    ThermalProfiles() = default;

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

struct ThermalParameterSet
{
    /** Own the ActionCalculator: FalseVacuumDecayRate holds it by reference, so
     *  binding it to a by-value constructor parameter would dangle the moment
     *  the constructor returned. Declared first so it outlives decay_rate. */
    ActionCalculator ac;

    /** Held behind unique_ptr so their addresses survive a move of this struct.
     *  friedmann_evolution holds const references to *decay_rate and *eos; with
     *  the objects stored by value those references would point into the
     *  moved-from instance. */
    std::unique_ptr<FalseVacuumDecayRate> decay_rate;
    std::unique_ptr<EquationOfState> eos;
    std::unique_ptr<FriedmannEvolution> friedmann_evolution;

    double TC;

    NucleationHistory nucleation_history;

    TransitionMilestone onset;
    TransitionMilestone percolation;
    TransitionMilestone completion;
    TransitionMilestone nucleation;

    ThermalProfiles profiles;

    ThermalParameterSet(
        Transition t_in, 
        ActionCalculator ac_in,
        double n_temp_action = 50,
        double n_temp_eos = 100,
        double vw = 1/sqrt(3.0),
        double background_dof = 106.75,
        double dof = 106.75,
        bool use_pf_in_nt_integrand = true,
        bool use_bag_dtdT = false,
        double percolation_target = 0.71,
        double completion_target = 1e-8,
        double onset_target = 1 - 1e-8,
        double nucleation_target = 1.00,
        double temperature_abs_tol = 1e-8,
        FalseVacuumDecayRate::PrefactorFunction prefactor = {}
    ) :
    ac(ac_in)
    {
        decay_rate = std::make_unique<FalseVacuumDecayRate>(t_in, ac);
        decay_rate->set_t_min(t_in.false_phase.T.front());
        decay_rate->set_t_max(t_in.TC);
        decay_rate->set_spline_evaluations(n_temp_action);
        if (prefactor) { decay_rate->set_prefactor_function(prefactor); }
        decay_rate->calculate();

        eos = std::make_unique<EquationOfState>(t_in);
        eos->set_n_temp(n_temp_eos);
        eos->set_background_dof(background_dof);
        eos->calculate();

        friedmann_evolution = std::make_unique<FriedmannEvolution>(*decay_rate, *eos);
        friedmann_evolution->set_vw(vw);
        friedmann_evolution->set_use_bag_dtdT(use_bag_dtdT);
        friedmann_evolution->set_percolation_target(percolation_target);
        friedmann_evolution->set_completion_target(completion_target);
        friedmann_evolution->set_onset_target(onset_target);
        friedmann_evolution->set_nucleation_target(nucleation_target);
        friedmann_evolution->set_temperature_abs_tol(temperature_abs_tol);

        // Solve only after every property is in place: vw, newtonG,
        // temperature_abs_tol and use_bag_dtdT are all read inside
        // refine_temperature_bounds() and evolve_friedmann().
        friedmann_evolution->solve();

        TC = decay_rate->get_t_max();

        friedmann_evolution->compute_milestones();

        friedmann_evolution->compute_nucleation_history(t_in.false_phase.T.front(), t_in.TC);
    }

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
};

class ThermoFinder {

    ActionCalculator ac;

    PROPERTY(PrintSettings, onset_print_setting, PrintSettings::MINIMAL);

    PROPERTY(PrintSettings, percolation_print_setting, PrintSettings::STANDARD);

    PROPERTY(PrintSettings, nucleation_print_setting, PrintSettings::STANDARD);

    PROPERTY(PrintSettings, completion_print_setting, PrintSettings::MINIMAL);

    PROPERTY(bool, compute_profiles, false);

    PROPERTY(bool, update_percolation_temperature, false);

    PROPERTY(double, n_temp_profiles, 250);

    PROPERTY(double, vw, 1/sqrt(3.0));

    PROPERTY(double, dof, 106.75);

    PROPERTY(double, background_dof, 106.75);

    PROPERTY(bool, use_pf_in_nt_integrand, true);

    PROPERTY(bool, use_bag_dtdT, false);

    PROPERTY(double, n_temp_eos, 100);

    PROPERTY(double, n_temp_action, 50);

    PROPERTY(double, n_temp_pf_nt, 200);

    PROPERTY(double, percolation_target, 0.71);

    PROPERTY(double, completion_target, 1e-6);

    PROPERTY(double, onset_target, 1 - 1e-6);

    PROPERTY(double, nucleation_target, 1.00);

    PROPERTY(double, temperature_abs_tol, 1e-6);

    /** Optional custom decay-rate prefactor (e.g. a BubbleDetPrefactor). When
     *  empty the FalseVacuumDecayRate uses its default analytic prefactor. */
    FalseVacuumDecayRate::PrefactorFunction prefactor_function;

public :

    ThermoFinder(ActionCalculator ac_in) : ac(ac_in) {};

    /** Install a custom decay-rate prefactor used for all subsequent
     *  get_thermal_parameter_set calls. */
    void set_prefactor_function(FalseVacuumDecayRate::PrefactorFunction f) { prefactor_function = f; }

    /** Whether already calculated all transitions */
    bool calculated_thermal_parameters = false;

    /** Container for all transitions between any two phases */
    std::vector<ThermalParameterSet> thermal_parameters;
    
    std::vector<ThermalParameterSet> find_thermal_parameters;

    ThermalParameterSet get_thermal_parameter_set(Transition t);

    const void add_thermal_parameter_values
    (
        TransitionMilestone& milestone, 
        const FalseVacuumDecayRate& decay_rate, 
        const EquationOfState& eos, 
        FriedmannEvolution& tm
    );

    void fill_nucleation_history
    (
        NucleationHistory& history, 
        TransitionMilestone& percolation, 
        TransitionMilestone& nucleation, 
        const FalseVacuumDecayRate& decay_rate, 
        FriedmannEvolution& tm
    );

    const double get_alpha(const double& temperature, const EquationOfState& eos, bool use_munu = false);

    const double get_betaH(const double& temperature, const FalseVacuumDecayRate& decay_rate);

    const double get_betaH_eff(const double& vw, const double& RsH);

    const double get_betaH_1(const double& temperature, const FalseVacuumDecayRate& decay_rate, FriedmannEvolution& tm);

    const double get_betaH_2(const double& temperature, const FalseVacuumDecayRate& decay_rate, FriedmannEvolution& tm);

    const double get_H(const double& temperature, FriedmannEvolution& tm);

    const double get_we(const double& temperature, const EquationOfState& eos);

    const std::pair<double, double> get_cs(const double& temperature, const EquationOfState& eos);

    const double get_n(const double& temperature, FriedmannEvolution& tm);

    const double get_Rbar(const double& temperature, FriedmannEvolution& tm);

    const double get_dt(const double& temperature, FriedmannEvolution& tm);

    const double get_percolation_temperature_wrapper(const double& vw, const double& percolation_target, const FriedmannEvolution& tm);

    const double get_vw_wrapper(const double& temperature, const FriedmannEvolution& tm, const EquationOfState& eos);

    const void revise_percolation_temperature(TransitionMilestone& percolation, const EquationOfState& eos, const FriedmannEvolution& tm);

};

} // namespace PhaseTracer

#endif // PHASETRACER_TEMPORARY_HPP_