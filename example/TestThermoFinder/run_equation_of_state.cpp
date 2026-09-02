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
        This is a test of the EquationOfState class. First, we read off the 
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
        PhaseTracer2 and this helper function just wraps them for convenience. 
        The transition we are interested in is the (0, s) -> (h, 0) transition.
    */
    auto transition = get_transition_from_parameters(lambda_hs, lambda_s, ms, Q, xi);
    
    /*
        We can then initialise the EquationOfState class using the transition. 
    */
    PhaseTracer::EquationOfState equation_of_state(transition);

    /*
        EquationOfState features two alterable properties: the number of spline 
        evaluations of the potential, and the background degrees of freedom. 
        This adds an extra term to the EoS of the form 
                            p_bg = g_bg x M_PI^2 x T^4 / 90, 
        which allows for additional degrees of freedom in the thermal plasma not 
        already counted by the effective potential. We set it to zero in this 
        case as all relativistic species are present in the potential. 
        
        On construction, these are initialised with the default values below.
    */

    equation_of_state.set_n_temp(200.0);
    equation_of_state.set_background_dof(0.0);

    /*
        We then calculate the equation of state and write to a csv.
    */
    equation_of_state.calculate();

    equation_of_state.write("example/TestThermoFinder/data/eos.csv");

    /*
        By default, the EoS is normalised by subtracting the energy density of
        the true vacuum at the minimum temperature. If the true vacuum
        corresponds to the zero-temperature ground state (as in this case), this
        is fine. However, we may be considering an intermediate transition for
        which this is no longer the case. In this case, the correct normalisation
        can be set explicitly.
    */
    const double normalisation = - 1e-8; // negative by definition
    equation_of_state.set_energy_norm(normalisation);

    equation_of_state.calculate();
    equation_of_state.write("example/TestThermoFinder/data/eos_w_norm.csv");

    return 0;
}