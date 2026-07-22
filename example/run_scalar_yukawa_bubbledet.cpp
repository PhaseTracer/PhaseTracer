/**
  PhaseTracer recreation of the BubbleDet reference script
  (the "BubbleDet" cell of example/TestThermalParameters/test_helpers.ipynb).

  For each temperature it finds the two minima of the d=3 scalar Yukawa model,
  builds the O(3) bounce with PhaseTracer's ActionCalculator (the analogue of the
  notebook's CosmoTransitions SingleFieldInstanton), and computes the one-loop
  prefactor S1 = -log(A) with BubbleDet through PhaseTracer's BubbleDetPrefactor
  interface. The columns mirror the notebook:

      T   S0   S1   S1_approx

  with S0 the bounce action, S1 the BubbleDet determinant, and
  S1_approx = -log[ T^4 (S0 / 2 pi T)^{3/2} ] the analytic estimate. A fallback
  (BubbleDet failing) shows up as S1 == S1_approx.

  Usage:  run_scalar_yukawa_bubbledet [python_exe] [worker_script.py]
*/

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "models/scalar_yukawa_3d.hpp"
#include "phasetracer.hpp"
#include "bubbledet_prefactor.hpp"

// Real roots of a x^3 + b x^2 + c x + d = 0 (a != 0), via Cardano/trig.
static std::vector<double> real_cubic_roots(double a, double b, double c, double d) {
  const double p = b / a, q = c / a, r = d / a;
  const double P = q - p * p / 3.0;            // depressed cubic t^3 + P t + Q
  const double Q = 2.0 * p * p * p / 27.0 - p * q / 3.0 + r;
  const double shift = -p / 3.0;               // x = t - p/3
  const double disc = Q * Q / 4.0 + P * P * P / 27.0;

  std::vector<double> roots;
  if (disc > 1e-12) {
    // One real root.
    const double sq = std::sqrt(disc);
    roots.push_back(std::cbrt(-Q / 2.0 + sq) + std::cbrt(-Q / 2.0 - sq) + shift);
  } else {
    // Three real roots (trigonometric form).
    const double m = 2.0 * std::sqrt(-P / 3.0);
    double arg = (3.0 * Q) / (2.0 * P) * std::sqrt(-3.0 / P);
    arg = std::clamp(arg, -1.0, 1.0);
    const double theta = std::acos(arg) / 3.0;
    for (int k = 0; k < 3; ++k) {
      roots.push_back(m * std::cos(theta - 2.0 * M_PI * k / 3.0) + shift);
    }
  }
  return roots;
}

int main(int argc, char *argv[]) {
  LOGGER(warning);

  const std::string python_exe =
      argc > 1 ? argv[1] : "/home/wsea0003/.venvs/phasetracer_venv/bin/python";
  const std::string worker_script =
      argc > 2 ? argv[2]
               : "/home/wsea0003/PhysicsCodes/deepphase_interface/newBranch/PhaseTracer/"
                 "example/BubbleDetBridge/bubbledet_worker.py";

  EffectivePotential::ScalarYukawa3D model;

  // Trace the phases over the transition window so the ThermoFinder section
  // below has a proper Transition. (The prefactor loop itself drives
  // get_action_full with the vacua directly and does not need this.) Doing it
  // before constructing ActionCalculator means ac's PhaseFinder copy has the
  // phases, matching run_Tp_test.
  PhaseTracer::PhaseFinder pf(model);
  pf.set_seed(0);
  pf.set_check_vacuum_at_high(false);
  pf.set_t_low(1.0);
  pf.set_t_high(20.0);
  pf.find_phases();

  // ActionCalculator: num_dims = 3 -> O(3) bounce, alpha = 2, matching the
  // notebook's alpha = dim - 1.
  PhaseTracer::ActionCalculator ac(pf);
  ac.set_num_dims(3);

  // Single-field model: BubbleDetPrefactor supplies the Higgs fluctuation
  // (W = d2V, spin 0, dof 1, zero_modes "Higgs", thermal) on its own -- exactly
  // the notebook's ParticleConfig.
  PhaseTracer::BubbleDetConfig config;
  config.python_executable = python_exe;
  config.worker_script = worker_script;
  config.dim = 3;
  config.thermal = true;
  PhaseTracer::BubbleDetPrefactor bd(model, config);

  const double T_max = 8.4;
  const double T_min = 7.9;
  const int nT = 51;

  std::cout << "# d=3 scalar Yukawa, BubbleDet via PhaseTracer\n";
  std::cout << std::setw(8) << "T" << std::setw(14) << "S0" << std::setw(14) << "S1"
            << std::setw(14) << "S1_approx" << std::setw(16) << "S0+S1" << "\n";

  for (int i = 0; i < nT; ++i) {
    const double T = T_max - (T_max - T_min) * i / (nT - 1);
    const auto p = model.params_3d(T);

    // Minima: real roots of dV/dx = s3 + m3sq x + (g3/2) x^2 + (lam3/6) x^3.
    auto roots = real_cubic_roots(p.lam3 / 6.0, p.g3 / 2.0, p.m3sq, p.s3);
    if (roots.size() < 3) {
      continue; // need two minima + a barrier
    }
    // Sort by potential; the two lowest are the minima (highest is the barrier).
    std::sort(roots.begin(), roots.end(), [&](double xa, double xb) {
      return model.V(Eigen::Matrix<double, 1, 1>(xa), T) <
             model.V(Eigen::Matrix<double, 1, 1>(xb), T);
    });
    const double phi_true = roots[0];
    const double phi_false = roots[1];

    PhaseTracer::ActionResult res;
    try {
      res = ac.get_action_full(Eigen::Matrix<double, 1, 1>(phi_true),
                               Eigen::Matrix<double, 1, 1>(phi_false), T);
    } catch (const std::exception &e) {
      continue;
    }
    const double S0 = res.action;
    if (!std::isfinite(S0) || S0 <= 0.) {
      continue;
    }

    const double A = bd(T, S0 / T, res); // BubbleDet prefactor
    const double S1 = -std::log(A);

    const double A_approx = std::pow(T, 4) * std::pow(S0 / (2.0 * M_PI * T), 1.5);
    const double S1_approx = -std::log(A_approx);

    std::cout << std::fixed << std::setprecision(4) << std::setw(8) << T
              << std::setprecision(6) << std::setw(14) << S0 << std::setw(14) << S1
              << std::setw(14) << S1_approx << std::setw(16) << (S0 + S1) << "\n";
  }

  // ---- Percolation temperature: default vs BubbleDet prefactor ------------
  // Solve for the percolation milestone with ThermoFinder (as in run_Tp_test),
  // once with the default analytic prefactor and once with BubbleDet, to see how
  // the determinant shifts Tp. Both runs share the same bounce action, so the
  // shift is purely the prefactor's effect.
  //
  // NOTE: FalseVacuumDecayRate forms the rate exponent as S/T, the 4d-thermal
  // convention; for this already-3d action that normalisation is not physical,
  // but it is identical in both runs, so the default-vs-BubbleDet comparison is
  // still meaningful.
  PhaseTracer::TransitionFinder tf(pf);
  tf.find_transitions();
  const auto transitions = tf.get_transitions();
  if (transitions.empty()) {
    std::cerr << "\nNo transition found; skipping percolation comparison.\n";
    return 0;
  }
  const auto trans = transitions[0];

  auto percolation_milestone = [&](bool use_bubbledet) {
    PhaseTracer::ThermoFinder tm(ac);
    tm.set_vw(1.0 / std::sqrt(3.0));
    tm.set_dof(3);
    tm.set_background_dof(3);
    if (use_bubbledet) {
      tm.set_prefactor_function(bd);
    }
    return tm.get_thermal_parameter_set(trans).percolation;
  };

  std::cout << "\n# Percolation temperature (TC = " << trans.TC << ")\n";
  for (bool use_bd : {false, true}) {
    const std::string label = use_bd ? "BubbleDet" : "analytic ";
    try {
      const auto perc = percolation_milestone(use_bd);
      std::cout << "  " << label << " : Tp = " << std::fixed << std::setprecision(5)
                << perc.temperature << "  (betaH = " << perc.betaH << ")\n";
    } catch (const std::exception &e) {
      std::cout << "  " << label << " : failed (" << e.what() << ")\n";
    } catch (...) {
      std::cout << "  " << label << " : failed (unknown exception)\n";
    }
  }

  return 0;
}
