#ifndef HELPER_FUNCTIONS_HPP_INCLUDED
#define HELPER_FUNCTIONS_HPP_INCLUDED

#include "phasetracer.hpp"
#include "hydrograv_interface.hpp"

constexpr double TL_BACKGROUND_DOF = 60.15297440;
constexpr double PHYSICAL_BACKGROUND_DOF = 66.25;

inline void
configure_phase_finder(PhaseTracer::PhaseFinder& phase_finder)
{
    phase_finder.set_seed(0);
    phase_finder.set_check_hessian_singular(true);
    phase_finder.set_check_vacuum_at_high(false);
    phase_finder.find_phases();
}

inline void
configure_transition_finder(PhaseTracer::TransitionFinder& transition_finder)
{
    transition_finder.find_transitions();
}

inline void
configure_action_calculator(PhaseTracer::ActionCalculator& action_calculator)
{
    action_calculator.set_action_calculator(PhaseTracer::ActionMethod::PathDeformation);
    action_calculator.set_PD_xtol(1e-6);
    action_calculator.set_PD_phitol(1e-6);
    action_calculator.set_PD_deformation_npoints(150);
}

inline void 
configure_thermo_finder(PhaseTracer::ThermoFinder& thermo_finder)
{
    thermo_finder.set_default_validation_method(PhaseTracer::ValidateMethod::TEMP);
    thermo_finder.set_temperature_threshold(1);
    thermo_finder.set_warm_start_chunk_size(0);
    thermo_finder.set_percolation_target(1-0.28957); // TL convention
    thermo_finder.set_eos_background_dof(TL_BACKGROUND_DOF);
    thermo_finder.set_percolation_print_setting(PhaseTracer::PrintSettings::VERBOSE);
    thermo_finder.set_nucleation_print_setting(PhaseTracer::PrintSettings::MINIMAL);

    thermo_finder.find_thermal_parameters();
}

#endif // HELPER_FUNCTIONS_HPP_INCLUDED