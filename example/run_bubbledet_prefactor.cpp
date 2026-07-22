/**
  End-to-end test of the production BubbleDet prefactor.

  Builds the 1D analytic test model, installs a BubbleDetPrefactor on a
  FalseVacuumDecayRate, and compares get_gamma() against the default analytic
  prefactor. The 1D test model derives from Potential (no built-in fluctuation
  spectrum), so we supply a single Higgs scalar with W = d2V/dphi2 via the
  config's spectrum_provider override.

  Usage:
    run_bubbledet_prefactor <python_exe> <worker_script.py>
  Defaults target the phasetracer venv and the in-tree worker script.
*/

#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "models/1D_test_model.hpp"
#include "phasetracer.hpp"
#include "bubbledet_prefactor.hpp"

int main(int argc, char *argv[]) {
  LOGGER(fatal);

  const std::string python_exe =
      argc > 1 ? argv[1] : "/home/wsea0003/.venvs/phasetracer_venv/bin/python";
  const std::string worker_script =
      argc > 2 ? argv[2]
               : "/home/wsea0003/PhysicsCodes/deepphase_interface/newBranch/PhaseTracer/"
                 "example/BubbleDetBridge/bubbledet_worker.py";

  EffectivePotential::OneDimModel model;

  PhaseTracer::PhaseFinder pf(model);
  pf.find_phases();

  PhaseTracer::ActionCalculator ac(pf);

  PhaseTracer::TransitionFinder tf(pf, ac);
  tf.find_transitions();
  const auto transitions = tf.get_transitions();
  if (transitions.empty()) {
    std::cerr << "no transition found\n";
    return 1;
  }
  const auto trans = transitions.at(0);

  // Single field: the functor supplies the longitudinal Higgs itself
  // (W = d2V/dphi2), so no spectrum is needed for this bare test model.
  PhaseTracer::BubbleDetConfig config;
  config.python_executable = python_exe;
  config.worker_script = worker_script;
  config.dim = 3;
  config.thermal = true;

  PhaseTracer::BubbleDetPrefactor bd_prefactor(model, config);

  // Narrow temperature window with a handful of spline knots (each BubbleDet
  // call costs ~seconds, so keep it small for a test).
  const double t_min = 40.0;
  const double t_max = 50.0;
  const int n_knots = 4;

  using clock = std::chrono::high_resolution_clock;
  using std::chrono::duration;

  // Time construction: this is where the prefactors actually differ -- the
  // BubbleDet rate runs n_knots determinants in get_splines, the analytic rate
  // just evaluates a closed-form expression. get_gamma() itself (below) is only
  // a spline lookup for both.
  auto c0 = clock::now();
  PhaseTracer::FalseVacuumDecayRate rate_bd(trans, ac, t_min, t_max, n_knots, bd_prefactor);
  auto c1 = clock::now();
  PhaseTracer::FalseVacuumDecayRate rate_analytic(trans, ac, t_min, t_max, n_knots);
  auto c2 = clock::now();

  const double build_bd_ms = duration<double, std::milli>(c1 - c0).count();
  const double build_an_ms = duration<double, std::milli>(c2 - c1).count();
  std::cout << "# Construction (get_splines, " << n_knots << " knots):\n";
  std::cout << "#   BubbleDet : " << std::fixed << std::setprecision(3) << build_bd_ms << " ms\n";
  std::cout << "#   analytic  : " << build_an_ms << " ms\n\n";

  std::cout << "# BubbleDet vs analytic prefactor (1D test model)\n";
  std::cout << std::setw(8) << "T" << std::setw(18) << "gamma_BD"
            << std::setw(18) << "gamma_analytic" << std::setw(16) << "ratio"
            << std::setw(16) << "t_BD[us]" << std::setw(16) << "t_an[us]" << "\n";
  for (double T : {43.0, 46.0, 48.0}) {
    auto t0 = clock::now();
    const double g_bd = rate_bd.get_gamma(T);
    auto t1 = clock::now();
    const double g_an = rate_analytic.get_gamma(T);
    auto t2 = clock::now();

    const double us_bd = duration<double, std::micro>(t1 - t0).count();
    const double us_an = duration<double, std::micro>(t2 - t1).count();

    std::cout << std::fixed << std::setprecision(1) << std::setw(8) << T
              << std::scientific << std::setprecision(6) << std::setw(18) << g_bd
              << std::setw(18) << g_an << std::setw(16) << g_bd / g_an
              << std::fixed << std::setprecision(3) << std::setw(16) << us_bd
              << std::setw(16) << us_an << "\n";
  }

  return 0;
}
