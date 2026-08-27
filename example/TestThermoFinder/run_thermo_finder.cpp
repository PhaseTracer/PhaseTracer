#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include <sys/stat.h>
#include <filesystem>
#include <sstream>
#include <chrono>
#include <omp.h>
#include <array>
#include <memory>
#include <unistd.h>

#include "phasetracer.hpp"
#include "helper_functions.hpp"
#include "models/xSM_MSbar.hpp"

using json = nlohmann::json;

json readFile(std::string fileName){
    std::ifstream file;
    file.open(fileName);

    if(file.fail()){
        throw std::runtime_error("error loading model parameters file!");
    }

    json data = json::parse(file);

    file.close();

    return data;
}

/*
  Main
*/
int main(int argc, char* argv[]) {

    LOGGER(fatal);

    /*
        This file demonstrates the ThermoFinder class. The performs the 
        calculation of thermal parameters for a given transition, and uses the 
        FalseVacuumDecayRate, EquationOfState, and FriedmannEvolution classes
        to facilitate a consistent treatment of all thermal parameters.
    */
    double ms, lambda_s, lambda_hs, Q, xi;

    try 
    {
        std::string json_filename = "example/TestThermoFinder/model_params.json";
        json modelParams = readFile(json_filename);

        ms = modelParams["ms"].get<double>();
        lambda_s = modelParams["lambda_s"].get<double>();
        lambda_hs = modelParams["lambda_hs"].get<double>();
        Q = modelParams["Q"].get<double>();
        xi = modelParams["xi"].get<double>();
    } catch (...) {
        std::cerr << "Running with default values." << std::endl;
        ms = 125;
        lambda_s = 1.00;
        lambda_hs = 1.05;
        Q = 100;
        xi = 1;
    }

    /*
        After reading off the model parameters, we run the PhaseFinder and 
        TransitionFinder parts of the calculation. These are unchanged from
        PhaseTracer2.
    */
    auto model = get_xSM_model_from_parameters(lambda_hs, lambda_s, ms, Q, xi);

    PhaseTracer::PhaseFinder phase_finder(model);
    phase_finder.set_seed(0);
    phase_finder.set_check_hessian_singular(true);
    phase_finder.set_check_vacuum_at_high(false);
    phase_finder.find_phases();
    std::cout << phase_finder;

    PhaseTracer::TransitionFinder transition_finder(phase_finder);
    transition_finder.find_transitions();
    std::cout << transition_finder;

    /*
        Because we use the FalseVacuumDecayRate class, we must also create an 
        instance of the ActionCalculator class.
    */
    PhaseTracer::ActionCalculator action_calculator(phase_finder);
    action_calculator.set_action_calculator
    (
        PhaseTracer::ActionMethod::PathDeformation
    );
    action_calculator.set_PD_xtol(1e-4);
    action_calculator.set_PD_phitol(1e-4);

    PhaseTracer::ThermoFinder thermo_finder(transition_finder, action_calculator);
    thermo_finder.set_default_validation_method(PhaseTracer::ValidateMethod::VEV);

    /*
        The calculation of thermal parameters is computationally expensive. As 
        such, we recommend screening transitions before running. We provide two 
        default validation methods. 

        The user can choose to validate using the temperature interval for the 
        transition:
            T_max - T_min > T_threshold,
        which can be set using:
            set_default_validation_method(PhaseTracer::ValidateMethod::TEMP);
        and T_threshold can be adjusted using:
            set_temperature_threshold(1.0)

        Alternatively, the VEV at the critical temperature can be used:
            norm(true_vacuum - false_vacuum) > vev_threshold,
        which can be set using:
            set_default_validation_method(PhaseTracer::ValidateMethod::VEV);
        and vev_threshold can be adjusted using:
            set_vev_threshold(1.0)

        Conversely, users wishing to maintain full control can pass their own
        filter. This has to be a function with the signature
        std::vector<PhaseTracer::Transition>(const std::vector<PhaseTracer::Transition>&),
        and can be set using set_transition_filter. We provide an example below
        selecting the (0, s) -> (h, 0) transition.
    */
    auto custom_validation = [](const std::vector<PhaseTracer::Transition>& input) 
    {
        std::vector<PhaseTracer::Transition> output;
        for (auto t : input) 
        {
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
    };
    thermo_finder.set_transition_filter(custom_validation);

    /*
        A given set of ThermalParameters contains a wealth of information,
        including the four transition milestones: onset, nucleation, percolation,
        and completion. Each milestone comes equiped with a full set of thermal
        parameters, but often only the percolation is of use. As such, the user
        can suppress printing the extra information by adjusting the print 
        settings. MINIMAL only prints the milestone temperature and status (ie 
        whether it occurs), whereas STANDARD includes alpha, beta, Rs, and vw.
        VERBSOE extends this to the full set of ThermalParameters the class 
        calculates.
    */
    thermo_finder.set_onset_print_setting(PhaseTracer::PrintSettings::MINIMAL);
    thermo_finder.set_nucleation_print_setting(PhaseTracer::PrintSettings::MINIMAL);
    thermo_finder.set_percolation_print_setting(PhaseTracer::PrintSettings::VERBOSE);
    thermo_finder.set_completion_print_setting(PhaseTracer::PrintSettings::MINIMAL);

    /*
        Then ThermoFinder works similar to PhaseFinder and TransitionFinder, 
        with find_thermal_parameters() and get_thermal_parameters() being
        the main methods and work analogously to the other classes.
    */
    thermo_finder.find_thermal_parameters();

    std::cout << thermo_finder << std::endl;

    return 0;
}