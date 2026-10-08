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

#include <vector>

using PhaseTracer::Config;
using PhaseTracer::Runner;

void bind_runner(py::module_ &m) {

  py::class_<Runner>(
      m, "Runner",
      "Runs the full PhaseTracer pipeline for a model:\n"
      "PhaseFinder -> ActionCalculator -> TransitionFinder -> ThermoFinder -> GravWaveCalculator.\n\n"
      "    runner = Runner(model, config)\n"
      "    status = runner.run()\n\n"
      "The Runner keeps the model alive. Stage objects and thermal parameter sets obtained from it\n"
      "raise RuntimeError once the Runner has been re-run; phases, transitions, milestones and\n"
      "spectra are copies and stay valid.")
      .def(py::init<EffectivePotential::Potential &, Config>(), py::arg("model"), py::arg("config") = Config(),
           py::keep_alive<1, 2>())
      .def("run", &Runner::run,
           py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect, py::gil_scoped_release>(),
           "Rebuild every stage from scratch and run up to config.pipeline.stop_after.\n"
           "Returns a RunStatus, or raises RunnerError if config.pipeline.throw_on_error is set.")
      .def_property_readonly("status", &Runner::status, "RunStatus of the last run()")
      .def_property(
          "config", [](Runner &r) -> Config & { return r.config(); },
          [](Runner &r, const Config &c) { r.config() = c; }, py::return_value_policy::reference_internal,
          "Settings used by the next run(); edit in place or assign a new Config")
      .def_property_readonly("run_id", &Runner::run_id, "Number of times run() has been called")
      .def("has", &Runner::has, py::arg("stage"), "Whether the stage object was built in the last run()")

      // stage objects (stale-safe handles; see bind_stages.cpp)
      .def("phase_finder", [](py::object self) { return make_handle(self, self.cast<Runner &>().phase_finder()); })
      .def("action_calculator",
           [](py::object self) { return make_handle(self, self.cast<Runner &>().action_calculator()); })
      .def("transition_finder",
           [](py::object self) { return make_handle(self, self.cast<Runner &>().transition_finder()); })
      .def("thermo_finder", [](py::object self) { return make_handle(self, self.cast<Runner &>().thermo_finder()); })
      .def("gravwave_calculator",
           [](py::object self) { return make_handle(self, self.cast<Runner &>().gravwave_calculator()); })

      // results
      .def("get_phases", [](const Runner &r) { return std::vector<PhaseTracer::Phase>(r.get_phases()); },
           "Copies of the phases found by PhaseFinder")
      .def("get_transitions",
           [](const Runner &r) { return std::vector<PhaseTracer::Transition>(r.get_transitions()); },
           "Copies of the transitions found by TransitionFinder")
      .def("get_thermal_parameters",
           [](py::object self) {
             auto tm = make_handle(self, self.cast<Runner &>().thermo_finder());
             std::vector<ThermalParameterSetHandle> sets;
             for (std::size_t i = 0; i < tm.get().get_thermal_parameters().size(); ++i) {
               sets.push_back(ThermalParameterSetHandle{tm, i});
             }
             return sets;
           },
           "Thermal parameters of each analysed transition (valid until the next run())")
      .def("get_spectra",
           [](const Runner &r) { return std::vector<PhaseTracer::GravWaveSpectrum>(r.get_spectra()); },
           "Copies of the gravitational wave spectra")
      .def("__repr__", [](const Runner &r) {
        return "<Runner run_id=" + std::to_string(r.run_id()) + " status=" +
               PhaseTracer::to_string(r.status().code) + ">";
      });
}
