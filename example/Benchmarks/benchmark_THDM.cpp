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

struct Results
{
    int flag;

    double lambda3;
    double Tcrit;
    double Tnuc;
    double Tperc;
    double Treh;
    double alpha;
    double betaH;
    double snrLISA;
    double runtime;

    Results() : flag(0), lambda3(0), Tcrit(0), Tnuc(0), Tperc(0), Treh(0), alpha(0), betaH(0), snrLISA(0), runtime(0) {}

    Results(const int& f, const double& l3) : flag(f), lambda3(l3), Tcrit(0), Tnuc(0), Tperc(0), Treh(0), alpha(0), betaH(0), snrLISA(0), runtime(0) {}

    Results(const int& f, const double& l3, const double& tc, const double& tn, const double& tp, const double& tr, const double& a, const double& bH, const double& snr, const double& rt)
        : flag(f), lambda3(l3), Tcrit(tc), Tnuc(tn), Tperc(tp), Treh(tr), alpha(a), betaH(bH), snrLISA(snr), runtime(rt) {}

    friend std::ostream& operator<<(std::ostream& os, const Results& res)
    {
        os << res.flag << ","
           << res.lambda3 << ","
           << res.Tcrit << ","
           << res.Tnuc << ","
           << res.Tperc << ","
           << res.Treh << ","
           << res.alpha << ","
           << res.betaH << ","
           << res.snrLISA << ","
           << res.runtime;
        return os;
    }
};

Results
get_results(
    const double& lambda3, 
    const double& action_tol=1e-6,
    const double& fRatioConv=0.01,
    const int& deformation_npoints=150,
    const int& action_steps=50,
    const bool& warm_start=true,
    const bool& debug=false, 
    const bool& print=false)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    EffectivePotential::THDM model; // lambda3 BP = 5.80958300958301
    model.init_params(21.3949, 1573.171, 0.28894, 0.26237, lambda3, -2.17494, -2.24170, 1);

    PhaseTracer::PhaseFinder pf(model);
    try {
        configure_phase_finder(pf);
    } catch (...) {
        return Results(-1, lambda3);
    }
    
    PhaseTracer::TransitionFinder tf(pf);
    try {
        configure_transition_finder(tf);
    } catch (...) {
        return Results(-2, lambda3);
    }

    PhaseTracer::ActionCalculator ac(pf);
    try {
        configure_action_calculator(ac, action_tol, fRatioConv, deformation_npoints);
    } catch (...) {
        return Results(-3, lambda3);
    }

    PhaseTracer::ThermoFinder tm(tf, ac);
    try {
        configure_thermo_finder(tm, action_steps, warm_start);
    } catch (...) {
        return Results(-4, lambda3);
    }

    if(print) 
    { 
        std::cout << pf;
        std::cout << tf;
        std::cout << tm; 
    }

    const auto& tps = tm.get_thermal_parameters();

    if(tps.size() == 0) { return Results(-5, lambda3);}

    if(debug)
    {
        auto& decay_rate = tps[0].get_decay_rate();
        auto t_min = decay_rate.get_t_min();
        auto t_max = decay_rate.get_t_max();
        int n_steps = 200;

        std::ofstream output("example/Benchmarks/data/THDM_debug_action.csv");
        output << "temperature,action,prefactor,decay_rate\n";
        double dt = (t_max - t_min) / (n_steps - 1);
        for (size_t i = 0; i < n_steps; ++i) {
            double t = t_min + i * dt;
            output << t << "," << decay_rate.get_action(t) << "," << decay_rate.get_prefactor(t) << "," << decay_rate.get_gamma(t) << "\n";
        }
        output.close();

        auto& friedmann = tps[0].get_friedmann_evolution();
        auto system = friedmann.system;
        system.write("example/Benchmarks/data/THDM_debug_friedmann.csv");
    }

    PhaseTracer::GravWaveCalculator gw(tm);
    try {
        double vw = tm.get_vw();
        configure_grav_wave_calculator(gw, vw);
    } catch (...) {
        return Results(-6, lambda3);
    }

    if(print) { std::cout << gw; }

    const auto& gws = gw.get_spectrums();

    if(gws.size() == 0) { return Results(-7, lambda3);}

    if(debug)
    {
        gw.write_spectrum_to_text(0, "example/Benchmarks/data/THDM_debug_spectrum.csv");
    }

    const auto& tp = tps[0];

    const double Tc = tp.TC;
    const double Tn = tp.nucleation.temperature;
    const double Tp = tp.percolation.temperature;
    const double Treh = tp.percolation.reheating_temperature;
    const double alpha = tp.percolation.alpha_bar;
    const double betaH = tp.percolation.betaH_eff;
    const double snrLISA = gws[0].SNR[0];

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();

    return Results(0, lambda3, Tc, Tn, Tp, Treh, alpha, betaH, snrLISA, elapsed_ms);
}


int main(int argc, char* argv[]) 
{
    double tol, fRatioConv;
    int nDeformations, nPoints;
    bool warm;
    if (argc == 2)
    {
        LOGGER(debug);
        double lambda3 = std::stod(argv[1]);

        Results res = get_results(lambda3, 1e-8, 0.005, 150, 50, true, true, true);
        std::cout << "Result for lambda3 = " << lambda3 << ": " << res << std::endl;
        return 0;
    } else if (argc == 1)
    {
        LOGGER(fatal);
        tol = 1e-8;
        fRatioConv = 0.005;
        nDeformations = 150;
        nPoints = 50;
        warm = true;
        std::cout << "Running scan with settings: tol = " << tol << ", fRatioConv = " << fRatioConv << ", nDeformations = " << nDeformations << ", nPoints = " << nPoints << ", warm = " << warm << std::endl;
    } else if (argc == 6)
    {
        LOGGER(fatal);
        tol = std::stod(argv[1]); // 1e-6
        fRatioConv = std::stod(argv[2]); // 0.01
        nDeformations = std::stoi(argv[3]); // 150
        nPoints = std::stoi(argv[4]); // 50
        warm = std::stoi(argv[5]); // true or false
        std::cout << "Running scan with settings: tol = " << tol << ", fRatioConv = " << fRatioConv << ", nDeformations = " << nDeformations << ", nPoints = " << nPoints << ", warm = " << warm << std::endl;
    } else {
        std::cerr << "Usage: " << argv[0] << " tol fRatioConv nDeformations nPoints warm" << std::endl;
        return 1;
    }

    std::ofstream output("example/Benchmarks/data/THDM_results.csv");
    output << "# status,lambda3,Tc,Tn,Tp,Treh,alpha,betaH,snrLISA,runtime\n";
    output << "# tol=" << tol << ", fRatioConv=" << fRatioConv << ", nDeformations=" << nDeformations << ", nPoints=" << nPoints << ", warm=" << warm << "\n";

    
    const int n_points = 200;
    const double lambda3_min = 5.5; // 5.5
    const double lambda3_max = 5.875; // 6.0
    for (int i = 0; i < n_points; ++i) {
        double lambda3 = lambda3_min + i * (lambda3_max - lambda3_min) / (n_points - 1);
        Results res = get_results(lambda3, tol, fRatioConv, nDeformations, nPoints, warm);
        output << std::setprecision(15) << res << "\n";
        output.flush();

        std::cout << "Ran point " << i + 1 << " with lambda3 = " << lambda3 << std::endl;
        std::cout << "Result: " << res << std::endl;
    }

    output.close();
    return 0;
}
