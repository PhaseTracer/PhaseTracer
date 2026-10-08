// ====================================================================
// This file is part of PhaseTracer

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ====================================================================

#include "bindings.hpp"

using PhaseTracer::RunnerError;
using PhaseTracer::RunStatus;
using PhaseTracer::Stage;
using PhaseTracer::StatusCode;

void bind_status(py::module_ &m) {

  py::enum_<Stage>(m, "Stage", "Stages of the pipeline, in the order they run")
      .value("None_", Stage::None)
      .value("Config", Stage::Config)
      .value("PhaseFinder", Stage::PhaseFinder)
      .value("ActionCalculator", Stage::ActionCalculator)
      .value("TransitionFinder", Stage::TransitionFinder)
      .value("ThermoFinder", Stage::ThermoFinder)
      .value("GravWave", Stage::GravWave)
      .def("__str__", [](Stage s) { return std::string(PhaseTracer::to_string(s)); });

  py::enum_<StatusCode>(
      m, "StatusCode",
      "Outcome of a run. The No* codes are physics outcomes (e.g. no first-order transition);\n"
      "the *Failed codes mean a stage threw.")
      .value("Success", StatusCode::Success)
      .value("InvalidConfig", StatusCode::InvalidConfig)
      .value("PhaseFinderFailed", StatusCode::PhaseFinderFailed)
      .value("NoPhases", StatusCode::NoPhases)
      .value("ActionCalculatorFailed", StatusCode::ActionCalculatorFailed)
      .value("TransitionFinderFailed", StatusCode::TransitionFinderFailed)
      .value("NoTransitions", StatusCode::NoTransitions)
      .value("ThermoFinderFailed", StatusCode::ThermoFinderFailed)
      .value("NoThermalParameters", StatusCode::NoThermalParameters)
      .value("GravWaveFailed", StatusCode::GravWaveFailed)
      .value("NoSpectra", StatusCode::NoSpectra)
      .def("__str__", [](StatusCode c) { return std::string(PhaseTracer::to_string(c)); });

  py::class_<RunStatus>(m, "RunStatus", "Result of Runner.run() or Config.validate(); truthy on success")
      .def(py::init<>())
      .def_readonly("code", &RunStatus::code)
      .def_readonly("stage", &RunStatus::stage, "Stage at which the run stopped")
      .def_readonly("message", &RunStatus::message)
      .def_readonly("warnings", &RunStatus::warnings, "Non-fatal problems, e.g. transitions that failed")
      .def("ok", &RunStatus::ok)
      .def("__bool__", &RunStatus::ok)
      .def("__str__", [](const RunStatus &s) {
        std::string text = repr_from_stream(s);
        while (!text.empty() && text.back() == '\n') {
          text.pop_back();
        }
        return text;
      })
      .def("__repr__", [](const RunStatus &s) {
        std::string repr = "<RunStatus " + std::string(PhaseTracer::to_string(s.code));
        if (s.stage != Stage::None) {
          repr += std::string(" at ") + PhaseTracer::to_string(s.stage);
        }
        return repr + ">";
      });

  // RunnerError carries the RunStatus as `.status`
  PYBIND11_CONSTINIT static py::gil_safe_call_once_and_store<py::object> runner_error;
  runner_error.call_once_and_store_result([&]() -> py::object {
    return py::exception<RunnerError>(m, "RunnerError", PyExc_RuntimeError);
  });
  m.attr("RunnerError").attr("__doc__") =
      "Raised by Runner.run() instead of returning a failed RunStatus when\n"
      "config.pipeline.throw_on_error is set; the status is available as .status";

  py::register_exception_translator([](std::exception_ptr p) {
    try {
      if (p) {
        std::rethrow_exception(p);
      }
    } catch (const RunnerError &e) {
      py::object type = runner_error.get_stored();
      py::object instance = type(e.what());
      instance.attr("status") = py::cast(e.status());
      PyErr_SetObject(type.ptr(), instance.ptr());
    }
  });
}
