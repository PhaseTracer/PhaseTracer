#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

#include "phasetracer.hpp"
#include "../TestThermoFinder/helper_functions.hpp"

/*
  The run_thermo_finder pipeline (example/TestThermoFinder/run_thermo_finder.cpp),
  driven through PhaseTracer::Config and PhaseTracer::Runner.
*/
int main(int argc, char* argv[]) {

    LOGGER(fatal);
    if (argc > 1 && std::string(argv[1]) == "-d")
    {
        LOGGER(debug);
    }

    double ms = 125, lambda_s = 1.00, lambda_hs = 1.05, Q = 100, xi = 1;
    std::ifstream file("example/TestThermoFinder/model_params.json");
    if (file)
    {
        const auto params = nlohmann::json::parse(file);
        ms = params["ms"].get<double>();
        lambda_s = params["lambda_s"].get<double>();
        lambda_hs = params["lambda_hs"].get<double>();
        Q = params["Q"].get<double>();
        xi = params["xi"].get<double>();
    }
    else
    {
        std::cerr << "Running with default values." << std::endl;
    }

    auto model = get_xSM_model_from_parameters(lambda_hs, lambda_s, ms, Q, xi);

    auto config = PhaseTracer::Config();

    config.phase_finder.seed = 0;
    config.phase_finder.check_hessian_singular = true;
    config.phase_finder.check_vacuum_at_high = false;

    config.action.method = PhaseTracer::ActionMethod::PathDeformation;
    config.action.PD_xtol = 1e-4;
    config.action.PD_phitol = 1e-4;

    config.thermo_finder.action_spline_evaluations = 50;
    config.thermo_finder.warm_start_chunk_size = 0;
    config.thermo_finder.eos_spline_evaluations = 250;
    config.thermo_finder.eos_background_dof = 0.0;
    config.thermo_finder.percolation_target = 0.71;

    // Select the (0, s) -> (h, 0) transition
    config.thermo_finder.transition_filter = [](const std::vector<PhaseTracer::Transition>& input)
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

    config.thermo_finder.onset_print_setting = PhaseTracer::PrintSettings::MINIMAL;
    config.thermo_finder.nucleation_print_setting = PhaseTracer::PrintSettings::MINIMAL;
    config.thermo_finder.percolation_print_setting = PhaseTracer::PrintSettings::VERBOSE;
    config.thermo_finder.completion_print_setting = PhaseTracer::PrintSettings::MINIMAL;

    config.gravwave.min_frequency = 1e-4;
    config.gravwave.max_frequency = 1e0;
    config.gravwave.num_frequency = 500;

    config.pipeline.to_print = true;

    PhaseTracer::Runner runner(model, config);
    const auto status = runner.run();
    std::cout << status;

    if (!status)
    {
        return 1;
    }

    const auto& percolation = runner.get_thermal_parameters()[0].percolation;
    std::cout << "Percolation T = " << percolation.temperature
              << ", alpha = " << percolation.alpha
              << ", betaH = " << percolation.betaH << std::endl;
    std::cout << "Number of GW spectra = " << runner.get_spectra().size() << std::endl;

    return 0;
}
