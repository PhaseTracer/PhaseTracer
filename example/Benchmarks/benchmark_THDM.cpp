/**
  2HDM+singlet DM
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <chrono>

#include "benchmark_models/THDM.hpp"
#include "phasetracer.hpp"
#include "helper_functions.hpp"


int main() {
    // start a timer
    auto start_time = std::chrono::high_resolution_clock::now();

    const bool debug_mode = true;
    
    // Set level of screen  output
    if (debug_mode) {
      LOGGER(debug);
    } else {
      LOGGER(fatal);
    }
    
    // Construct our model
    EffectivePotential::THDM model;
    model.init_params(21.3949, 1573.171, 0.28894, 0.26237, 5.80958300958301, -2.17494, -2.24170, 1);

    PhaseTracer::PhaseFinder pf(model);
    configure_phase_finder(pf);
    std::cout << pf;
    
    PhaseTracer::TransitionFinder tf(pf);
    configure_transition_finder(tf);
    std::cout << tf;

    PhaseTracer::ActionCalculator ac(pf);
    configure_action_calculator(ac);

    PhaseTracer::ThermoFinder tm(tf, ac);
    configure_thermo_finder(tm);

    std::cout << tm;

    PhaseTracer::GravWaveCalculator gw(tm);

    gw.set_gw_method(PhaseTracer::GravWaveMethod::FitFormulae);

    // print time before calculating gravitational wave spectrums
    auto gw_start_time = std::chrono::high_resolution_clock::now();
    gw.calc_spectrums();

    std::cout << gw;
    auto gw_end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> gw_elapsed = gw_end_time - gw_start_time;
    std::cout << "Elapsed time for gravitational wave calculation: " << gw_elapsed.count() << " seconds" << std::endl;
    
    // stop the timer and print the elapsed time
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

    return 0;

}
