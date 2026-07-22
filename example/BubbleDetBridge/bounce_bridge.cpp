// ============================================================================
//  bounce_bridge.cpp
//
//  Validation-harness bridge: exposes a PhaseTracer bounce (action + profile)
//  and the underlying effective potential to Python, so BubbleDet can be run on
//  the *PhaseTracer* solution and cross-checked against CosmoTransitions.
//
//  This is step 6.1 of the BubbleDet interface plan -- it de-risks the
//  units/dimension mapping before any in-C++ BubbleDet plumbing is committed.
//
//  Python module `pt_bounce`:
//
//      import pt_bounce
//      solver = pt_bounce.TestModelBounceSolver(num_dims=3)
//      solver.TC()                       -> critical temperature
//      solver.t_min()                    -> lowest temperature with a profile
//      b = solver.action_full(T)         -> BounceResult
//      b.action, b.R, b.Phi, b.dPhi      -> bounce action and profile arrays
//      b.phi_true, b.phi_false           -> the two vacua at T
//      solver.V(x, T) / dV / d2V         -> potential and derivatives (1 field)
//
//  Only the 1D analytic test model is wired up here on purpose: it has an exact
//  potential, so any disagreement is in the interface, not the physics.
// ============================================================================

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <Eigen/Core>

#include <memory>
#include <stdexcept>
#include <vector>

#include "models/1D_test_model.hpp"
#include "phasetracer.hpp"

namespace py = pybind11;

// ----------------------------------------------------------------------------
//  Plain-data result handed back to Python.
// ----------------------------------------------------------------------------
struct BounceResult {
  double T = 0.0;
  double action = 0.0; // S (NOT divided by T)
  double phi_true = 0.0;
  double phi_false = 0.0;
  Eigen::VectorXd R;
  Eigen::VectorXd Phi;
  Eigen::VectorXd dPhi;
};

// ----------------------------------------------------------------------------
//  Drives the full PhaseTracer pipeline once for the 1D test model and serves
//  bounce solutions at requested temperatures.
// ----------------------------------------------------------------------------
class TestModelBounceSolver {
public:
  explicit TestModelBounceSolver(int num_dims)
      : pf(model), ac(pf), tf(pf, ac) {
    // Keep the bridge quiet unless something goes wrong.
    LOGGER(fatal);

    ac.set_num_dims(static_cast<size_t>(num_dims));

    pf.find_phases();
    tf.find_transitions();

    transitions = tf.get_transitions();
    if (transitions.empty()) {
      throw std::runtime_error("no transition found for the 1D test model");
    }
  }

  double TC() const { return transitions.at(0).TC; }

  double t_min() const { return transitions.at(0).false_phase.T.front(); }

  // Full bounce at temperature T between the (first) transition's phases.
  BounceResult action_full(double T) const {
    const auto &trans = transitions.at(0);

    const auto result = ac.get_action_full(trans.true_phase, trans.false_phase, T);
    const auto vacua = ac.get_vacua_at_T(trans.true_phase, trans.false_phase, T, 0);

    BounceResult out;
    out.T = T;
    out.action = result.action;
    out.phi_false = vacua.at(0)(0); // false_phase point at T
    out.phi_true = vacua.at(1)(0);  // true_phase point at T
    out.R = result.bubble_profile.R;
    out.Phi = result.bubble_profile.Phi;
    out.dPhi = result.bubble_profile.dPhi;
    return out;
  }

  // Single-field potential and derivatives (the model has one scalar).
  double V(double x, double T) const {
    return model.V(field(x), T);
  }
  double dV(double x, double T) const {
    return model.dV_dx(field(x), T)(0);
  }
  double d2V(double x, double T) const {
    return model.d2V_dx2(field(x), T)(0, 0);
  }

private:
  static Eigen::VectorXd field(double x) {
    Eigen::VectorXd phi(1);
    phi(0) = x;
    return phi;
  }

  EffectivePotential::OneDimModel model;
  PhaseTracer::PhaseFinder pf;
  PhaseTracer::ActionCalculator ac;
  PhaseTracer::TransitionFinder tf;
  std::vector<PhaseTracer::Transition> transitions;
};

// ----------------------------------------------------------------------------
//  Module definition
// ----------------------------------------------------------------------------
PYBIND11_MODULE(pt_bounce, m) {
  m.doc() = "PhaseTracer bounce/profile bridge for BubbleDet validation (1D test model)";

  py::class_<BounceResult>(m, "BounceResult")
      .def_readonly("T", &BounceResult::T)
      .def_readonly("action", &BounceResult::action, "Bounce action S (not divided by T).")
      .def_readonly("phi_true", &BounceResult::phi_true)
      .def_readonly("phi_false", &BounceResult::phi_false)
      .def_readonly("R", &BounceResult::R, "Radial grid of the bounce profile.")
      .def_readonly("Phi", &BounceResult::Phi, "Field value along the profile.")
      .def_readonly("dPhi", &BounceResult::dPhi, "dPhi/dR along the profile.");

  py::class_<TestModelBounceSolver>(m, "TestModelBounceSolver")
      .def(py::init<int>(), py::arg("num_dims") = 3,
           "Run the PhaseTracer pipeline for the 1D analytic test model.")
      .def("TC", &TestModelBounceSolver::TC, "Critical temperature.")
      .def("t_min", &TestModelBounceSolver::t_min,
           "Lowest temperature with a profile (false-phase floor).")
      .def("action_full", &TestModelBounceSolver::action_full, py::arg("T"),
           "Bounce action and profile at temperature T.")
      .def("V", &TestModelBounceSolver::V, py::arg("x"), py::arg("T"))
      .def("dV", &TestModelBounceSolver::dV, py::arg("x"), py::arg("T"))
      .def("d2V", &TestModelBounceSolver::d2V, py::arg("x"), py::arg("T"));
}
