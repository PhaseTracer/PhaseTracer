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

        In this file, we demonstrate how HydroGrav is seemlessly incorporated as
        a backend for the GravWaveCalculator class. This allows users to utilise
        this much more powerful calculation, without needing to depertment from
        the familiar PhaseTracer UI.

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

    try {
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

    /*
        The 'configure_' methods below are defined in the helper_functions 
        header. Refer to this file for the settings we use below.
    */

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
        With this established, we can now calculate the SSM spectrum using
        HydroGrav. We have implemented this within the GravWaveCalculator class,
        which we build below.
    */
    PhaseTracer::GravWaveCalculator gravwave_calculator(thermo_finder);

    /*
        To select HydroGrav, we change gw_method using the setter below. This 
        can be either 'FitFormulae' or 'SoundShell', where the former is the 
        default. The SSM calculation evaluates the spectrum on a grid of 
        dimensionless k * R_s values, which are then converted to values of f_0.
        As such, these should be checked. 
    */
    gravwave_calculator.set_gw_method(PhaseTracer::GravWaveMethod::SoundShell);
    gravwave_calculator.set_n_kRs_value(100);
    gravwave_calculator.set_min_kRs_value(1e-3);
    gravwave_calculator.set_max_kRs_value(1e3);

    /*
        We then call calc_spectrums as usual. Note that this won't evaluate the
        collision or turbulence spectra.
    */
    gravwave_calculator.calc_spectrums();

    /*
        Everything GravWaveCalculator offers then works as usual - the summed
        spectrum, the pretty-printer, and writing to text.
    */
    std::cout << gravwave_calculator;

    /*
        The fluid profiles used in the HydroGrav calculation are also stored in 
        the spectrum object, and will be empty for the FitFormulae calculation.
        These are accessed using spectrum.profile, which returns the new 
        FluidProfile object. This is a direct PhaseTracer translation of the 
        equivalent fluid profile object in HydroGrav, so we refer readers to 
        2606.27775 for more information.
    */
    const auto& spectrums = gravwave_calculator.get_spectrums();
    for (int ii = 0; ii < spectrums.size(); ii++) {
        gravwave_calculator.write_spectrum_to_text(
            spectrums[ii], "example/TestHydroGrav/data/spectrum_" + std::to_string(ii) + ".csv");
        spectrums[ii].profile.write_profile_to_text(
            "example/TestHydroGrav/data/profiles_" + std::to_string(ii) + ".csv");
    }

    /*
        Printing the elapsed time for the full execution of the pipeline.
    */
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

    return 0;
}