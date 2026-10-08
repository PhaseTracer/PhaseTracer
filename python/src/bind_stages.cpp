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

#include <string>
#include <utility>
#include <vector>

// The stage objects are owned by a Runner and only reachable through it. Each is exposed as a
// Handle (see bindings.hpp), which keeps the Runner alive and whose methods first check that the
// Runner has not been re-run.

using namespace PhaseTracer;

namespace {

template <typename T>
std::string handle_repr(const Handle<T> &h) {
  return repr_from_stream(h.get());
}

} // namespace

void bind_stages(py::module_ &m) {

  using PF = Handle<PhaseFinder>;
  py::class_<PF>(m, "PhaseFinder", "The Runner's PhaseFinder (obtain with Runner.phase_finder())")
      .def("get_phases", [](const PF &h) { return std::vector<Phase>(h.get().get_phases()); })
      .def("find_minima_at_t", [](const PF &h, double T) { return h.get().find_minima_at_t(T); }, py::arg("T"),
           "All minima of the potential at temperature T")
      .def("get_phases_at_T", [](const PF &h, double T) { return h.get().get_phases_at_T(T); }, py::arg("T"),
           "Phases that exist at temperature T")
      .def("phase_at_T", [](const PF &h, const Phase &phase, double T) { return h.get().phase_at_T(phase, T); },
           py::arg("phase"), py::arg("T"), "Location of a phase at temperature T")
      .def("get_deepest_phase_at_T", [](const PF &h, double T) { return h.get().get_deepest_phase_at_T(T); },
           py::arg("T"))
      .def("__repr__", &handle_repr<PhaseFinder>);

  using AC = Handle<ActionCalculator>;
  py::class_<AC>(m, "ActionCalculator", "The Runner's ActionCalculator (obtain with Runner.action_calculator())")
      .def("get_action",
           [](const AC &h, const Phase &phase1, const Phase &phase2, double T) {
             py::gil_scoped_release release;
             return h.get().get_action(phase1, phase2, T);
           },
           py::arg("phase1"), py::arg("phase2"), py::arg("T"), "Bounce action S between two phases at T")
      .def("get_action",
           [](const AC &h, const Eigen::VectorXd &true_vacuum, const Eigen::VectorXd &false_vacuum, double T) {
             py::gil_scoped_release release;
             return h.get().get_action(true_vacuum, false_vacuum, T);
           },
           py::arg("true_vacuum"), py::arg("false_vacuum"), py::arg("T"), "Bounce action S between two vacua at T")
      .def("get_action_full",
           [](const AC &h, const Phase &phase1, const Phase &phase2, double T) {
             py::gil_scoped_release release;
             return h.get().get_action_full(phase1, phase2, T);
           },
           py::arg("phase1"), py::arg("phase2"), py::arg("T"), "Bounce action with its bubble profile")
      .def("get_action_full",
           [](const AC &h, const Eigen::VectorXd &true_vacuum, const Eigen::VectorXd &false_vacuum, double T) {
             py::gil_scoped_release release;
             return h.get().get_action_full(true_vacuum, false_vacuum, T);
           },
           py::arg("true_vacuum"), py::arg("false_vacuum"), py::arg("T"));

  using TF = Handle<TransitionFinder>;
  py::class_<TF>(m, "TransitionFinder", "The Runner's TransitionFinder (obtain with Runner.transition_finder())")
      .def("get_transitions", [](const TF &h) { return std::vector<Transition>(h.get().get_transitions()); })
      .def("__repr__", &handle_repr<TransitionFinder>);

  using DR = Handle<FalseVacuumDecayRate>;
  py::class_<DR>(m, "FalseVacuumDecayRate", "Splined bounce action and decay rate of one transition")
      .def("get_action", [](const DR &h, double T) { return h.get().get_action(T); }, py::arg("T"),
           "Bounce action S(T)")
      .def("get_action_deriv", [](const DR &h, double T) { return h.get().get_action_deriv(T); }, py::arg("T"))
      .def("get_action_double_deriv", [](const DR &h, double T) { return h.get().get_action_double_deriv(T); },
           py::arg("T"))
      .def("get_gamma", [](const DR &h, double T) { return h.get().get_gamma(T); }, py::arg("T"),
           "Decay rate per unit volume Gamma(T)")
      .def("get_prefactor", [](const DR &h, double T) { return h.get().get_prefactor(T); }, py::arg("T"))
      .def_property_readonly("t_min", [](const DR &h) { return h.get().get_t_min(); })
      .def_property_readonly("t_max", [](const DR &h) { return h.get().get_t_max(); })
      .def("get_bubble_profile", [](const DR &h, double T) { return h.get().get_bubble_profile(T); }, py::arg("T"))
      .def("write", [](const DR &h, const std::string &filename, int n_steps) { h.get().write(filename, n_steps); },
           py::arg("filename"), py::arg("n_steps") = 100, "Write T, action, prefactor and rate to a CSV file");

  using EOS = Handle<EquationOfState>;
  py::class_<EOS>(m, "EquationOfState", "Equation of state of one transition; getters return (false, true) phase values")
      .def("get_energy", [](const EOS &h, double T) { return h.get().get_energy(T); }, py::arg("T"))
      .def("get_pressure", [](const EOS &h, double T) { return h.get().get_pressure(T); }, py::arg("T"))
      .def("get_enthalpy", [](const EOS &h, double T) { return h.get().get_enthalpy(T); }, py::arg("T"))
      .def("get_entropy", [](const EOS &h, double T) { return h.get().get_entropy(T); }, py::arg("T"))
      .def("get_sound_speed", [](const EOS &h, double T) { return h.get().get_sound_speed(T); }, py::arg("T"))
      .def_property_readonly("t_min", [](const EOS &h) { return h.get().get_t_min(); })
      .def_property_readonly("t_max", [](const EOS &h) { return h.get().get_t_max(); })
      .def("write", [](const EOS &h, const std::string &filename) { h.get().write(filename); }, py::arg("filename"));

  using TPS = ThermalParameterSetHandle;
  py::class_<TPS>(m, "ThermalParameterSet", "Thermal parameters of one transition (from Runner.get_thermal_parameters())")
      .def_property_readonly("TC", [](const TPS &h) { return h.get().TC; })
      .def_property_readonly("onset", [](const TPS &h) { return h.get().onset; })
      .def_property_readonly("nucleation", [](const TPS &h) { return h.get().nucleation; })
      .def_property_readonly("percolation", [](const TPS &h) { return h.get().percolation; })
      .def_property_readonly("completion", [](const TPS &h) { return h.get().completion; })
      .def_property_readonly("profiles", [](const TPS &h) { return h.get().profiles; },
                             "Thermal profiles (filled when thermo_finder.compute_profiles is set)")
      .def_property_readonly("nucleation_history", [](const TPS &h) { return h.get().nucleation_history; })
      .def("decay_rate", [](const TPS &h) { return sibling_handle(h.thermo_finder, h.get().get_decay_rate()); },
           "Splined action and decay rate of this transition")
      .def("equation_of_state",
           [](const TPS &h) { return sibling_handle(h.thermo_finder, h.get().get_equation_of_state()); },
           "Equation of state of this transition")
      .def("__repr__", [](const TPS &h) { return repr_from_stream(h.get()); });

  using TM = Handle<ThermoFinder>;
  py::class_<TM>(m, "ThermoFinder", "The Runner's ThermoFinder (obtain with Runner.thermo_finder())")
      .def("get_thermal_parameters",
           [](const TM &h) {
             std::vector<TPS> sets;
             for (std::size_t i = 0; i < h.get().get_thermal_parameters().size(); ++i) {
               sets.push_back(TPS{h, i});
             }
             return sets;
           })
      .def("get_failure_messages", [](const TM &h) { return h.get().get_failure_messages(); },
           "Why transitions that passed the filter have no thermal parameters")
      .def("__repr__", &handle_repr<ThermoFinder>);

  using GW = Handle<GravWaveCalculator>;
  py::class_<GW>(m, "GravWaveCalculator", "The Runner's GravWaveCalculator (obtain with Runner.gravwave_calculator())")
      .def("get_spectrums", [](const GW &h) { return std::vector<GravWaveSpectrum>(h.get().get_spectrums()); })
      .def("get_total_spectrum", [](const GW &h) { return h.get().get_total_spectrum(); })
      .def("get_SNR", [](const GW &h, const TransitionMilestone &milestone) { return h.get().get_SNR(milestone); },
           py::arg("milestone"), "[LISA, Taiji] signal-to-noise ratio for a milestone")
      .def("write_spectrum_to_text",
           [](const GW &h, const std::string &filename) { h.get().write_spectrum_to_text(filename); },
           py::arg("filename"), "Write the total spectrum")
      .def("write_spectrum_to_text",
           [](const GW &h, int i, const std::string &filename) { h.get().write_spectrum_to_text(i, filename); },
           py::arg("i"), py::arg("filename"), "Write the i-th spectrum")
      .def("__repr__", &handle_repr<GravWaveCalculator>);
}
