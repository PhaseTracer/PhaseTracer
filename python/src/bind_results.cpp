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

// Plain-data results. They are returned to Python by copy, so they stay valid after the
// Runner is re-run or destroyed.

void bind_results(py::module_ &m) {

  // ------------------------------------------------------------------ phases

  py::enum_<PhaseTracer::phase_end_descriptor>(m, "PhaseEnd", "Why a phase stopped being traced")
      .value("REACHED_T_STOP", PhaseTracer::REACHED_T_STOP)
      .value("FORBIDDEN_OR_BOUNDS", PhaseTracer::FORBIDDEN_OR_BOUNDS)
      .value("HESSIAN_SINGULAR", PhaseTracer::HESSIAN_SINGULAR)
      .value("HESSIAN_NOT_POSITIVE_DEFINITE", PhaseTracer::HESSIAN_NOT_POSITIVE_DEFINITE)
      .value("JUMP_INDICATED_END", PhaseTracer::JUMP_INDICATED_END);

  py::class_<PhaseTracer::Point>(m, "Point", "A minimum of the potential at one temperature")
      .def_readonly("x", &PhaseTracer::Point::x)
      .def_readonly("potential", &PhaseTracer::Point::potential)
      .def_readonly("t", &PhaseTracer::Point::t)
      .def("__repr__", &repr_from_stream<PhaseTracer::Point>);

  py::class_<PhaseTracer::Phase>(m, "Phase", "A minimum traced in temperature; T is ascending")
      .def_readonly("key", &PhaseTracer::Phase::key)
      .def_readonly("X", &PhaseTracer::Phase::X)
      .def_readonly("T", &PhaseTracer::Phase::T)
      .def_readonly("dXdT", &PhaseTracer::Phase::dXdT)
      .def_readonly("V", &PhaseTracer::Phase::V)
      .def_readonly("redundant", &PhaseTracer::Phase::redundant)
      .def_readonly("end_low", &PhaseTracer::Phase::end_low)
      .def_readonly("end_high", &PhaseTracer::Phase::end_high)
      .def("contains_t", &PhaseTracer::Phase::contains_t, py::arg("T"))
      .def("__repr__", &repr_from_stream<PhaseTracer::Phase>);

  // ------------------------------------------------------------------ transitions

  py::enum_<PhaseTracer::Message>(m, "Message", "Outcome of a critical-temperature search")
      .value("SUCCESS", PhaseTracer::SUCCESS)
      .value("NON_OVERLAPPING_T", PhaseTracer::NON_OVERLAPPING_T)
      .value("ERROR", PhaseTracer::ERROR);

  py::class_<PhaseTracer::Transition>(m, "Transition", "A transition between two phases at its critical temperature TC")
      .def_readonly("message", &PhaseTracer::Transition::message)
      .def_readonly("TC", &PhaseTracer::Transition::TC)
      .def_readonly("true_phase", &PhaseTracer::Transition::true_phase)
      .def_readonly("false_phase", &PhaseTracer::Transition::false_phase)
      .def_readonly("true_vacuum", &PhaseTracer::Transition::true_vacuum)
      .def_readonly("false_vacuum", &PhaseTracer::Transition::false_vacuum)
      .def_readonly("gamma", &PhaseTracer::Transition::gamma)
      .def_readonly("changed", &PhaseTracer::Transition::changed)
      .def_readonly("delta_potential", &PhaseTracer::Transition::delta_potential)
      .def_readonly("key", &PhaseTracer::Transition::key)
      .def_readonly("id", &PhaseTracer::Transition::id)
      .def_readonly("TN", &PhaseTracer::Transition::TN)
      .def_readonly("true_vacuum_TN", &PhaseTracer::Transition::true_vacuum_TN)
      .def_readonly("false_vacuum_TN", &PhaseTracer::Transition::false_vacuum_TN)
      .def_readonly("subcritical", &PhaseTracer::Transition::subcritical)
      .def("__repr__", &repr_from_stream<PhaseTracer::Transition>);

  // ------------------------------------------------------------------ bounce action

  py::class_<PhaseTracer::Profile1D>(m, "Profile1D", "Bubble profile: field Phi and its derivative dPhi against radius R")
      .def_readonly("R", &PhaseTracer::Profile1D::R)
      .def_readonly("Phi", &PhaseTracer::Profile1D::Phi)
      .def_readonly("dPhi", &PhaseTracer::Profile1D::dPhi);

  py::class_<PhaseTracer::ActionResult>(m, "ActionResult", "Bounce action and the solution it came from")
      .def_readonly("action", &PhaseTracer::ActionResult::action)
      .def_readonly("bubble_profile", &PhaseTracer::ActionResult::bubble_profile)
      .def_readonly("phi_for_profile", &PhaseTracer::ActionResult::phi_for_profile)
      .def_readonly("tunneling_path", &PhaseTracer::ActionResult::tunneling_path);

  // ------------------------------------------------------------------ thermal parameters

  py::enum_<PhaseTracer::MilestoneStatus>(m, "MilestoneStatus", "Whether a milestone was reached")
      .value("YES", PhaseTracer::YES)
      .value("FAST", PhaseTracer::FAST)
      .value("NO", PhaseTracer::NO)
      .value("ERR", PhaseTracer::ERR);

  py::enum_<PhaseTracer::NucleationType>(m, "NucleationType")
      .value("EXPONENTIAL", PhaseTracer::EXPONENTIAL)
      .value("SIMULTANEOUS", PhaseTracer::SIMULTANEOUS);

  py::class_<PhaseTracer::TransitionMilestone>(
      m, "TransitionMilestone", "Thermal parameters at one milestone (onset, nucleation, percolation or completion)")
      .def_readonly("type", &PhaseTracer::TransitionMilestone::type)
      .def_readonly("status", &PhaseTracer::TransitionMilestone::status)
      .def_readonly("nucleation_type", &PhaseTracer::TransitionMilestone::nucleation_type)
      .def_readonly("temperature", &PhaseTracer::TransitionMilestone::temperature)
      .def_readonly("reheating_temperature", &PhaseTracer::TransitionMilestone::reheating_temperature)
      .def_readonly("vw", &PhaseTracer::TransitionMilestone::vw)
      .def_readonly("alpha", &PhaseTracer::TransitionMilestone::alpha)
      .def_readonly("alpha_bar", &PhaseTracer::TransitionMilestone::alpha_bar)
      .def_readonly("alpha_munu", &PhaseTracer::TransitionMilestone::alpha_munu)
      .def_readonly("g_eff", &PhaseTracer::TransitionMilestone::g_eff)
      .def_readonly("h_eff", &PhaseTracer::TransitionMilestone::h_eff)
      .def_readonly("betaH", &PhaseTracer::TransitionMilestone::betaH)
      .def_readonly("beta1H", &PhaseTracer::TransitionMilestone::beta1H)
      .def_readonly("beta2H", &PhaseTracer::TransitionMilestone::beta2H)
      .def_readonly("betaH_eff", &PhaseTracer::TransitionMilestone::betaH_eff)
      .def_readonly("H", &PhaseTracer::TransitionMilestone::H)
      .def_readonly("we", &PhaseTracer::TransitionMilestone::we)
      .def_readonly("cs_plus", &PhaseTracer::TransitionMilestone::cs_plus)
      .def_readonly("cs_minus", &PhaseTracer::TransitionMilestone::cs_minus)
      .def_readonly("n", &PhaseTracer::TransitionMilestone::n)
      .def_readonly("Rs", &PhaseTracer::TransitionMilestone::Rs)
      .def_readonly("Rbar", &PhaseTracer::TransitionMilestone::Rbar)
      .def_readonly("dt", &PhaseTracer::TransitionMilestone::dt)
      .def("__repr__", &repr_from_stream<PhaseTracer::TransitionMilestone>);

  py::class_<PhaseTracer::ThermalProfiles>(m, "ThermalProfiles", "Thermal history on a temperature grid")
      .def_readonly("temperature", &PhaseTracer::ThermalProfiles::temperature)
      .def_readonly("dtdT", &PhaseTracer::ThermalProfiles::dtdT)
      .def_readonly("time", &PhaseTracer::ThermalProfiles::time)
      .def_readonly("hubble_rate", &PhaseTracer::ThermalProfiles::hubble_rate)
      .def_readonly("bounce_action", &PhaseTracer::ThermalProfiles::bounce_action)
      .def_readonly("extended_volume", &PhaseTracer::ThermalProfiles::extended_volume)
      .def_readonly("false_vacuum_decay_rate", &PhaseTracer::ThermalProfiles::false_vacuum_decay_rate)
      .def_readonly("false_vacuum_fraction", &PhaseTracer::ThermalProfiles::false_vacuum_fraction)
      .def_readonly("d_false_vacuum_fraction", &PhaseTracer::ThermalProfiles::d_false_vacuum_fraction)
      .def_readonly("nucleation_rate", &PhaseTracer::ThermalProfiles::nucleation_rate)
      .def_readonly("mean_bubble_separation", &PhaseTracer::ThermalProfiles::mean_bubble_separation)
      .def_readonly("mean_bubble_radius", &PhaseTracer::ThermalProfiles::mean_bubble_radius);

  py::class_<PhaseTracer::NucleationHistory>(m, "NucleationHistory")
      .def_readonly("nucleation_type", &PhaseTracer::NucleationHistory::nucleation_type)
      .def_readonly("betaH_1", &PhaseTracer::NucleationHistory::betaH_1)
      .def_readonly("betaH_2", &PhaseTracer::NucleationHistory::betaH_2)
      .def_readonly("T_m", &PhaseTracer::NucleationHistory::T_m);

  // ------------------------------------------------------------------ gravitational waves

  py::class_<PhaseTracer::FluidProfile>(m, "FluidProfile", "Self-similar fluid profile (SoundShell only)")
      .def_readonly("xi", &PhaseTracer::FluidProfile::xi)
      .def_readonly("v", &PhaseTracer::FluidProfile::v)
      .def_readonly("w", &PhaseTracer::FluidProfile::w)
      .def_readonly("lambda", &PhaseTracer::FluidProfile::lambda)
      .def_readonly("T", &PhaseTracer::FluidProfile::T)
      .def_readonly("xi_min", &PhaseTracer::FluidProfile::xi_min)
      .def_readonly("xi_max", &PhaseTracer::FluidProfile::xi_max)
      .def_readonly("cs_plus_sq", &PhaseTracer::FluidProfile::cs_plus_sq)
      .def_readonly("cs_minus_sq", &PhaseTracer::FluidProfile::cs_minus_sq)
      .def_readonly("mode", &PhaseTracer::FluidProfile::mode)
      .def_readonly("shock_converged", &PhaseTracer::FluidProfile::shock_converged)
      .def("empty", &PhaseTracer::FluidProfile::empty)
      .def("mode_str", &PhaseTracer::FluidProfile::mode_str)
      .def("__repr__", &repr_from_stream<PhaseTracer::FluidProfile>);

  py::class_<PhaseTracer::GravWaveSpectrum>(
      m, "GravWaveSpectrum", "Gravitational wave spectrum h^2 Omega(f); SNR is [LISA, Taiji]")
      .def_readonly("Tref", &PhaseTracer::GravWaveSpectrum::Tref)
      .def_readonly("Treh", &PhaseTracer::GravWaveSpectrum::Treh)
      .def_readonly("vw", &PhaseTracer::GravWaveSpectrum::vw)
      .def_readonly("alpha", &PhaseTracer::GravWaveSpectrum::alpha)
      .def_readonly("alpha_fit", &PhaseTracer::GravWaveSpectrum::alpha_fit)
      .def_readonly("g_eff", &PhaseTracer::GravWaveSpectrum::g_eff)
      .def_readonly("h_eff", &PhaseTracer::GravWaveSpectrum::h_eff)
      .def_readonly("cs", &PhaseTracer::GravWaveSpectrum::cs)
      .def_readonly("beta_H", &PhaseTracer::GravWaveSpectrum::beta_H)
      .def_readonly("peak_frequency", &PhaseTracer::GravWaveSpectrum::peak_frequency)
      .def_readonly("peak_amplitude", &PhaseTracer::GravWaveSpectrum::peak_amplitude)
      .def_readonly("frequency", &PhaseTracer::GravWaveSpectrum::frequency)
      .def_readonly("sound_wave", &PhaseTracer::GravWaveSpectrum::sound_wave)
      .def_readonly("turbulence", &PhaseTracer::GravWaveSpectrum::turbulence)
      .def_readonly("bubble_collision", &PhaseTracer::GravWaveSpectrum::bubble_collision)
      .def_readonly("total_amplitude", &PhaseTracer::GravWaveSpectrum::total_amplitude)
      .def_readonly("lisa_noise", &PhaseTracer::GravWaveSpectrum::lisa_noise)
      .def_readonly("taiji_noise", &PhaseTracer::GravWaveSpectrum::taiji_noise)
      .def_readonly("SNR", &PhaseTracer::GravWaveSpectrum::SNR)
      .def_readonly("method", &PhaseTracer::GravWaveSpectrum::method)
      .def_readonly("kRs", &PhaseTracer::GravWaveSpectrum::kRs)
      .def_readonly("profile", &PhaseTracer::GravWaveSpectrum::profile)
      .def_readonly("dtau", &PhaseTracer::GravWaveSpectrum::dtau)
      .def("__repr__", &repr_from_stream<PhaseTracer::GravWaveSpectrum>);
}
