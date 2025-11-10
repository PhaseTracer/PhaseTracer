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
#include <interpolation.h>

#include "property.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "action_calculator.hpp"
#include "false_vacuum_decay_rate.hpp"
#include "equation_of_state.hpp"
#include "transition_metrics.hpp"

namespace PhaseTracer {

struct ThermalParameterSet 
{
    FalseVacuumDecayRate decay_rate;
    EquationOfState eos;
    TransitionMetrics transition_metrics;
    double TC;

    TransitionMilestone onset;
    TransitionMilestone percolation;
    TransitionMilestone completion;
    TransitionMilestone nucleation;

    ThermalParameterSet(
        Transition t_in, 
        TransitionFinder tf_in,
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
        double temperature_abs_tol = 1e-8
    ) :
    decay_rate(t_in, tf_in, t_in.false_phase.T.front(), t_in.TC, n_temp_action),
    eos(t_in, n_temp_eos, background_dof), 
    transition_metrics(decay_rate, eos), 
    TC(decay_rate.get_t_max()) 
    {
        transition_metrics.set_vw(vw);
        transition_metrics.set_dof(dof);
        transition_metrics.set_use_pf_in_nt_integrand(use_pf_in_nt_integrand);
        transition_metrics.set_use_bag_dtdT(use_bag_dtdT);
        transition_metrics.set_percolation_target(percolation_target);
        transition_metrics.set_completion_target(completion_target);
        transition_metrics.set_onset_target(onset_target);
        transition_metrics.set_nucleation_target(nucleation_target);
        transition_metrics.set_temperature_abs_tol(temperature_abs_tol);
        
        transition_metrics.compute_milestones();
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

class Temporary {

    TransitionFinder tf;

    ActionCalculator ac;

    PROPERTY(PrintSettings, onset_print_setting, PrintSettings::MINIMAL);

    PROPERTY(PrintSettings, percolation_print_setting, PrintSettings::STANDARD);

    PROPERTY(PrintSettings, nucleation_print_setting, PrintSettings::STANDARD);

    PROPERTY(PrintSettings, completion_print_setting, PrintSettings::MINIMAL);

    PROPERTY(double, vw, 1/sqrt(3.0));

    PROPERTY(double, dof, 106.75);

    PROPERTY(double, background_dof, 106.75);

    PROPERTY(bool, use_pf_in_nt_integrand, true);

    PROPERTY(bool, use_bag_dtdT, false);

    PROPERTY(double, n_temp_eos, 100);

    PROPERTY(double, n_temp_action, 50);

    PROPERTY(double, n_temp_pf_nt, 200);

    PROPERTY(double, percolation_target, 0.71);

    PROPERTY(double, completion_target, 1e-8);

    PROPERTY(double, onset_target, 1 - 1e-8);

    PROPERTY(double, nucleation_target, 1.00);

    PROPERTY(double, temperature_abs_tol, 1e-8);

public : 

    Temporary(TransitionFinder tf_in, ActionCalculator ac_in) : 
    tf(tf_in), ac(ac_in) {};

    // std::vector<ThermalParameterSet> find_thermal_parameters;

    ThermalParameterSet get_thermal_parameter_set(Transition t);

    const void add_thermal_parameter_values(TransitionMilestone& milestone, const FalseVacuumDecayRate& decay_rate, const EquationOfState& eos, TransitionMetrics& tm);

    const double get_alpha(const double& temperature, const EquationOfState& eos);

    const double get_betaH(const double& temperature, const FalseVacuumDecayRate& decay_rate);

    const double get_H(const double& temperature, TransitionMetrics& tm);

    const double get_we(const double& temperature, const EquationOfState& eos);

    const std::pair<double, double> get_cs(const double& temperature, const EquationOfState& eos);

    const double get_n(const double& temperature, TransitionMetrics& tm);

    const double get_Rbar_integral(const double& temperature, TransitionMetrics& tm);

    const double get_dt(const double& temperature, TransitionMetrics& tm);

};

} // namespace PhaseTracer

#endif // PHASETRACER_TEMPORARY_HPP_