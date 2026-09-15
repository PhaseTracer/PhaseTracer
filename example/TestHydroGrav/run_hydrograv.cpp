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
#include "hydrograv_interface.hpp"
#include "helper_functions.hpp"
#include "models/xSM_MSbar.hpp"

using json = nlohmann::json;

json readFile(std::string fileName)
{
    std::ifstream file;
    file.open(fileName);

    if(file.fail())
    {
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

    /*
        This file demonstrates the GW calculation using HydroGrav. By default, 
        the GravWaveCalculator class of PhaseTracer3 utilises fitting formulas 
        for the calculation of the GW spectrum. However, we support interfacing 
        directly to the recently released code HydroGrav that performs consistent 
        hydrodynamics and a calculation of the acoustic spectrum using the sound 
        shell model.

        In this file, we demonstrate how to use HydroGrav very explicitly. In a
        separate example, we demonstrate how this process is streamlined by
        hiding the HydroGrav calculation away as a backend for GravWaveCalculator.

        As we demonstrate using the new ThermoFinder functionality elsewhere, we
        quickly perform the calculation of themal parameters below.
    
        Adding a timer to benchmark the full execution time for the self consistent
        pipeline.
    */
    auto start = std::chrono::high_resolution_clock::now();

    LOGGER(fatal);

    double ms, lambda_s, lambda_hs, Q, xi;

    if (argc > 1 && std::string(argv[1]) == "-d")
    {
        LOGGER(debug);
    }

    try 
    {
        std::string json_filename = "example/TestHydroGrav/model_params.json";
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

    auto model = get_xSM_model_from_parameters(lambda_hs, lambda_s, ms, Q, xi);

    PhaseTracer::PhaseFinder phase_finder(model);
    configure_phase_finder(phase_finder);

    PhaseTracer::TransitionFinder transition_finder(phase_finder);
    configure_transition_finder(transition_finder);

    PhaseTracer::ActionCalculator action_calculator(phase_finder);
    configure_action_calculator(action_calculator);

    PhaseTracer::ThermoFinder thermo_finder(transition_finder, action_calculator);
    configure_thermo_finder(thermo_finder);

    const auto& thermal_parameter_sets = thermo_finder.get_thermal_parameters();
    if (thermal_parameter_sets.empty())
    {
        LOG(fatal) << "No thermal parameters found.";
        return 1;
    }

    const auto& thermal_parameters = thermal_parameter_sets[0];
    std::cout << thermal_parameters;

    /*
        To proceed, we first need extract usable inputs for HydroGrav from our 
        ThermalParameterSet. First, we can extract our equation of state and 
        milestone. We will use the percolation milestone.
    */
    auto& percolation = thermal_parameters.percolation;
    auto& eos = thermal_parameters.get_equation_of_state();

    /*
        PhaseTracer and HydroGrav each have their own representation of the 
        equation of state and transition parameters objects, and HydroGrav has 
        an additional parameter in the form of the Universe object. The
        PhaseTracer::HydroGravBridge namespace provides methods for translating
        between them.

        We start with the Universe and EquationOfState objects.
    */
    namespace bridge = PhaseTracer::HydroGravBridge;

    /*
        HydroGrav reports its progress through its own logger. This line maps 
        the HydroGrav logger to the default PhaseTracer LOGGER above.
    */
    bridge::install_logger();

    auto dof = 107.75;
    HydroGrav::PhaseTransition::Universe universe
        = bridge::to_hydrograv_universe(percolation, dof);

    HydroGrav::PhaseTransition::EquationOfState
    hydrograv_eos = bridge::to_hydrograv_eos(eos);

    /*
        Then, we can create the PTParams_Veff object. We could build this from
        the pieces above, but we also provide another translation method for 
        this job.
    */
    auto vw = 0.577;
    auto pt_params = bridge::to_hydrograv_pt_params(
        thermal_parameters,
        percolation,
        vw, 
        dof
    );

    /*
        With this, we can use HydroGrav to solve for both the fluid profiles
        and the GW spectrum.

        Note that pt_params must stay in scope for as long as the spectrum does:
        PowerSpec and the FluidProfile it carries keep a non-owning pointer back
        to it, and write() below dereferences that pointer.
    */
    auto kRs_values = logspace(-3, 3, 200);
    HydroGrav::Spectrum::PowerSpec spectrum = HydroGrav::Spectrum::GWSpec(kRs_values, pt_params);

    /*
        Then, we can access the kRs, frequency, and amplitude vectors, as well
        as the fluid profiles as follows. For more information on the output, 
        see the following reference 2606.27775.
    */
    auto& profile = spectrum.profile();

    spectrum.write("example/TestHydroGrav/data/spectrum.csv");
    profile.write("example/TestHydroGrav/data/profiles.csv");
    
    /*
        Printing the elapsed time for the full execution of the pipeline.
    */
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

    return 0;
}