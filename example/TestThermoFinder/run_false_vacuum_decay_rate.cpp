#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include <sys/stat.h>
#include <filesystem>
#include <sstream>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
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

  TODO
    * add write function to decay rate
*/
int main(int argc, char* argv[]) {

    LOGGER(fatal);
    
    if (argc > 1 && std::string(argv[1]) == "-d")
    {
        LOGGER(debug);
    }

    /*
        This is a test of the FalseVacuumDecayRate class. First, we read off the 
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
        To proceed, we can construct the FalseVacuumDecayRate class. This 
        requires a copy of the ActionCalculator class, as well as the transition 
        we are interested in.

        For the action calculator, we note the code has only been tested using 
        the built in PD method. In addition, we recommend increasing the 
        tolerances for a smoother fit.
    */
    PhaseTracer::ActionCalculator action_calculator(phase_finder);
    action_calculator.set_action_calculator
    (
        PhaseTracer::ActionMethod::PathDeformation
    );
    action_calculator.set_PD_xtol(1e-4);
    action_calculator.set_PD_phitol(1e-4);

    PhaseTracer::FalseVacuumDecayRate decay_rate(transition, action_calculator);

    /*
        FalseVacuumDecayRate features three alterable properties: the minimum 
        and maximum temperatures of the spline, and the number of action 
        evaluations. On construction, these are initialised with the 
        default values below.
    */
    decay_rate.set_t_min(transition.false_phase.T.front());
    decay_rate.set_t_max(transition.TC);
    decay_rate.set_spline_evaluations(50);

    /*
        We then call calculate() to solve the bounce action and fit the splines.
        
        This calculation supports parallelisation via OpenMP, and is the most
        expensive part of the entire PhaseTracer pipeline. Below, we map the 
        time taken with and without OpenMP. The number of threads can be set via 
        the OMP_NUM_THREADS environment variable.
    */

    auto start = std::chrono::high_resolution_clock::now();

    decay_rate.calculate();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

#ifdef _OPENMP
    std::cout << "Action calculation on " << omp_get_max_threads()
        << " thread(s) took " << elapsed.count() << " seconds.\n";

    /*
        We then set the number of threads to 1 and rerun the calculation
        again.
    */
    omp_set_num_threads(1);
    start = std::chrono::high_resolution_clock::now();

    decay_rate.calculate();

    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Action calculation on 1 thread took " << elapsed.count()
        << " seconds.\n";
#else
    std::cout << "Action calculation (serial build) took " << elapsed.count()
        << " seconds.\n";
#endif
        
    /*
        We can then access the computed action, prefactor, and decay rate 
        splines as a function of the false vacuum temperature. We use the write
        method to store these in a file.
    */
    decay_rate.write("example/TestThermoFinder/data/decay_rate.csv");

    /*
        The FalseVacuumDecayRate class also supports using a custom prefactor
        function. By default, we use A(T) = T^4 x (S_3(T) / 2 pi T)^3/2, but this
        can be changed as follows. In particular, this functionality supports
        interfacing PhaseTracer with BubbleDet.

        First, we define a custom prefactor function. This function takes the 
        temperature and bounce solution and returns the prefactor.
    */
    auto custom_prefactor = []
    (
        double temperature, 
        double action_on_T, 
        const PhaseTracer::ActionResult& bounce
    ) {
        double T_4 = temperature*temperature*temperature*temperature;
        double prefactor = 1.1 * T_4 * std::pow(action_on_T / (2.*M_PI), 3./2.);
        return prefactor;
    };

    /*
        Then, the new custom prefactor is called using the associated setter.
        Naturally, calculate must be called again to rebuild the splines with 
        the new prefactor.
    */
    decay_rate.set_prefactor_function(custom_prefactor);
    decay_rate.calculate();
    decay_rate.write("example/TestThermoFinder/data/decay_rate_w_prefactor.csv");

    return 0;
}