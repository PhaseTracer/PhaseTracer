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

    if (argc > 1 && std::string(argv[1]) == "-d")
    {
        LOGGER(debug);
    }

    /*
        This is a test of the FriedmannEvolution class. First, we read off the 
        model parameters.
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

    PhaseTracer::TransitionFinder transition_finder(phase_finder);
    transition_finder.find_transitions();
    
    auto filtered_transitions = get_scalar_to_higgs_transition
    (
        transition_finder.get_transitions()
    );
    auto transition = filtered_transitions[0];

    if(transition.message == PhaseTracer::Message::SUCCESS) 
    {
        std::cout << "Transition found:\n";
        std::cout << transition;
    } else {
        std::cout << "No transition found.\n";
    }

    /*
        The function above automatically selects the transition we are interested
        in: the (0, s) -> (h, 0) transition. For all functionality below, we 
        strongly recommend filtering the transitions instead of blindly running
        them.
    */

    /*
        The purpose of this file is to use the new FriedmannEvolution class. To
        use this, we need to provide both the FalseVacuumDecayRate and 
        EquationOfState classes. Information on these can be found in their
        respective files. We initialise each below.
    */

    PhaseTracer::ActionCalculator action_calculator(phase_finder);
    action_calculator.set_action_calculator(PhaseTracer::ActionMethod::PathDeformation);
    action_calculator.set_PD_xtol(1e-6);
    action_calculator.set_PD_phitol(1e-6);

    PhaseTracer::FalseVacuumDecayRate decay_rate(transition, action_calculator);
    decay_rate.set_spline_evaluations(50);
    decay_rate.calculate();

    PhaseTracer::EquationOfState eos(transition);
    eos.set_n_temp(200);
    eos.calculate();

    /*
        Then, we can initialise the FriedmannEvolution class using the above
        members. Note, it is not necessary to run calculate() for either 
    */
    PhaseTracer::FriedmannEvolution friedmann_evolution(decay_rate, eos);

    /*
        As with the other classes, 
    */
    friedmann_evolution.set_vw(0.577);

    /*
        We solve the F-JKAK equations by calling solve().
    */
    friedmann_evolution.solve();

    /*
        Once this is solved, we can access the solution to the F-JMAK
        equations (called system).
    */
    friedmann_evolution.system.write("example/TestThermalParameters/data/friedmann_system.csv");

    /*
        FriedmannEvolution is primarily used as a tool for the ThermoFinder 
        class, and as such contains many useful functions. This includes methods 
        that evaluate splines fit to the system.
    */
    double Tref = 75.0;
    std::cout << "False vacuum fraction at T = 75.0 GeV = " <<
        friedmann_evolution.get_false_vacuum_fraction(Tref) << "\n";

    /*
        In addition, 
    */

    return 0;
}