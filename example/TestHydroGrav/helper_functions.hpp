#ifndef HELPER_FUNCTIONS_HPP_INCLUDED
#define HELPER_FUNCTIONS_HPP_INCLUDED

#include "phasetracer.hpp"
#include "models/xSM_MSbar.hpp"

std::vector<PhaseTracer::Transition>
get_scalar_to_higgs_transition(std::vector<PhaseTracer::Transition> input)
{
    std::vector<PhaseTracer::Transition> output;
    for (auto t : input) 
    {
        // vac = [h, s]
        auto true_vac = t.true_vacuum;
        auto false_vac = t.false_vacuum;
        auto changed = t.changed;

        if(changed[0] && changed[1])
        {
        double h_true = true_vac[0];
        double h_false = false_vac[0];
        double s_true = true_vac[1];
        double s_false = false_vac[1];
        if( (abs(h_true) > 5. && abs(s_true) < 1e-3) && (abs(s_false) > 5. && abs(h_false) < 1e-3) )
        {
            output.push_back(t);
        }     
        }
    }
    return output;
}

EffectivePotential::xSM_MSbar 
get_xSM_model_from_parameters(const double lambda_hs, const double lambda_s, const double ms, const double Q, const double xi)
{
    bool use_1L_EWSB_in_0L_mass = false;
    bool use_Goldstone_resum = true;
    bool tree_level_tadpoles = false;
    bool use_covariant_gauge = false;
    auto model = EffectivePotential::xSM_MSbar::from_tadpoles(
        lambda_hs, 
        lambda_s, 
        ms, 
        Q, 
        xi, 
        use_covariant_gauge, 
        use_1L_EWSB_in_0L_mass, 
        use_Goldstone_resum, 
        tree_level_tadpoles, 
        {}
    );
    model.set_daisy_method(EffectivePotential::DaisyMethod::Parwani);
    return model;
}

PhaseTracer::Transition
get_transition_from_parameters(const double lambda_hs, const double lambda_s, const double ms, const double Q, const double xi)
{
    auto model = get_xSM_model_from_parameters(lambda_hs, lambda_s, ms, Q, xi);
    
    PhaseTracer::PhaseFinder phase_finder(model);
    phase_finder.set_seed(0);
    phase_finder.set_check_hessian_singular(true);
    phase_finder.set_check_vacuum_at_high(false);
    phase_finder.find_phases();

    PhaseTracer::TransitionFinder transition_finder(phase_finder);
    transition_finder.find_transitions();
    
    auto filtered_transitions = get_scalar_to_higgs_transition
    (
        transition_finder.get_transitions()
    );
    auto transition = filtered_transitions[0];

    return transition;
}

void
configure_phase_finder(PhaseTracer::PhaseFinder& phase_finder)
{
    phase_finder.set_seed(0);
    phase_finder.set_check_hessian_singular(true);
    phase_finder.set_check_vacuum_at_high(false);
    phase_finder.find_phases();
}

void
configure_transition_finder(PhaseTracer::TransitionFinder& transition_finder)
{
    transition_finder.find_transitions();
}

void
configure_action_calculator(PhaseTracer::ActionCalculator& action_calculator)
{
    action_calculator.set_action_calculator(PhaseTracer::ActionMethod::PathDeformation);
    action_calculator.set_PD_xtol(1e-4);
    action_calculator.set_PD_phitol(1e-4);
}

void 
configure_thermo_finder(PhaseTracer::ThermoFinder& thermo_finder)
{
    auto custom_validation = [](const std::vector<PhaseTracer::Transition>& input)
    {
        std::vector<PhaseTracer::Transition> output;
        std::copy_if(input.begin(), input.end(), std::back_inserter(output),
            [](const PhaseTracer::Transition& t) {
                const auto& tv = t.true_vacuum;
                const auto& fv = t.false_vacuum;
                return t.changed[0] && t.changed[1] &&
                    std::abs(tv[0]) > 5.  && std::abs(tv[1]) < 1e-3 &&
                    std::abs(fv[1]) > 5.  && std::abs(fv[0]) < 1e-3;
            });
        return output;
    };
    thermo_finder.set_transition_filter(custom_validation);
    thermo_finder.find_thermal_parameters();
}

#endif // HELPER_FUNCTIONS_HPP_INCLUDED