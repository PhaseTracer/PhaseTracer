#ifndef HELPER_FUNCTIONS_HPP_INCLUDED
#define HELPER_FUNCTIONS_HPP_INCLUDED

#include "phasetracer.hpp"

constexpr double TL_BACKGROUND_DOF = 60.15297440;
constexpr double PHYSICAL_BACKGROUND_DOF = 66.25;

inline void
configure_phase_finder(PhaseTracer::PhaseFinder& phase_finder)
{
    phase_finder.set_seed(0);
    phase_finder.set_check_hessian_singular(false);
    phase_finder.set_hessian_singular_rel_tol(1e-4);
    phase_finder.set_check_vacuum_at_high(false);
    phase_finder.set_hessian_eig_max_rel_change(0.);
    phase_finder.set_check_midpoint_hessian(false);
    // phase_finder.set_t_high(500);
    phase_finder.find_phases();
}

inline void
configure_transition_finder(PhaseTracer::TransitionFinder& transition_finder)
{
    transition_finder.find_transitions();
}

inline void
configure_action_calculator(
    PhaseTracer::ActionCalculator& action_calculator, 
    const double& tol=1e-6,
    const double& fRatioConv=0.01,
    const int& deformation_npoints=150
)
{
    action_calculator.set_action_calculator(PhaseTracer::ActionMethod::PathDeformation);
    action_calculator.set_PD_xtol(tol);
    action_calculator.set_PD_phitol(tol);
    action_calculator.set_PD_fRatioConv(fRatioConv);
    action_calculator.set_PD_deformation_npoints(deformation_npoints);
}

inline void 
configure_thermo_finder(
    PhaseTracer::ThermoFinder& thermo_finder, 
    const int& action_steps=50, 
    const bool& warm_start=true
)
{
    auto custom_validation = [](const std::vector<PhaseTracer::Transition>& input)
    {
        std::vector<PhaseTracer::Transition> output;
        std::copy_if(input.begin(), input.end(), std::back_inserter(output),
            [](const PhaseTracer::Transition& t) {
                const auto& tv = t.true_vacuum;
                const auto& fv = t.false_vacuum;
                // first ensure the fv is at the origin, so both values less than 1e-3
                bool fv_at_origin = std::abs(fv[0]) < 1e-3 && std::abs(fv[1]) < 1e-3;
                // then check the temperature range is above 1
                bool valid_temp = t.TC - t.false_phase.T.front() > 1;
                // lastly ensure Tc is less than 150
                bool valid_tc = t.TC < 150;
                return fv_at_origin && valid_temp && valid_tc;
            });
        return output;
    };
    thermo_finder.set_transition_filter(custom_validation);
    
    // thermo_finder.set_default_validation_method(PhaseTracer::ValidateMethod::TEMP);
    // thermo_finder.set_temperature_threshold(1);

    int warm_start_chunk;
    if(warm_start) { warm_start_chunk=0; } else { warm_start_chunk=1; }
    thermo_finder.set_warm_start_chunk_size(warm_start_chunk);
    thermo_finder.set_action_spline_evaluations(action_steps);
    thermo_finder.set_update_percolation_temperature(true);
    thermo_finder.set_action_smoothing_window(1);
    thermo_finder.set_action_smoothing_order(3);
    thermo_finder.set_action_laurent_tail(true);
    thermo_finder.set_percolation_target(1-0.28957); // TL convention
    thermo_finder.set_eos_background_dof(TL_BACKGROUND_DOF);
    thermo_finder.set_percolation_print_setting(PhaseTracer::PrintSettings::VERBOSE);
    thermo_finder.set_nucleation_print_setting(PhaseTracer::PrintSettings::MINIMAL);

    thermo_finder.find_thermal_parameters();
}

inline void
configure_grav_wave_calculator(PhaseTracer::GravWaveCalculator& gravwave_calculator, const double& vw = 0.3)
{
    gravwave_calculator.set_gw_method(PhaseTracer::GravWaveMethod::FitFormulae);
    gravwave_calculator.set_default_milestone(PhaseTracer::MilestoneType::PERCOLATION);
    gravwave_calculator.set_vw(1.0);
    gravwave_calculator.set_SNR_f_max(1e1);
    gravwave_calculator.set_max_frequency(1e1);
    gravwave_calculator.set_SNR_f_min(1e-5);
    gravwave_calculator.set_min_frequency(1e-5);
    gravwave_calculator.set_num_frequency(500);
    // gravwave_calculator.set_g_eff(68.9);
    // gravwave_calculator.set_h_eff(68.9);
    // gravwave_calculator.set_use_legacy_gw_methods(true);
    // gravwave_calculator.set_T_threshold_bubble_collision(1e10);
    gravwave_calculator.calc_spectrums();
}

#endif // HELPER_FUNCTIONS_HPP_INCLUDED