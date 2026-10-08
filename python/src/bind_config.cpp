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

#include <optional>

#include "config.hpp"

namespace sev = boost::log::trivial;

void bind_config(py::module_ &m) {

  // ------------------------------------------------------------------ enums used by Config

  py::enum_<sev::severity_level>(m, "LogLevel", "Severity threshold of PhaseTracer's log output")
      .value("trace", sev::trace)
      .value("debug", sev::debug)
      .value("info", sev::info)
      .value("warning", sev::warning)
      .value("error", sev::error)
      .value("fatal", sev::fatal);

  py::enum_<PhaseTracer::ActionMethod>(m, "ActionMethod", "Method for the bounce action")
      .value("None_", PhaseTracer::ActionMethod::None)
      .value("BubbleProfiler", PhaseTracer::ActionMethod::BubbleProfiler)
      .value("PathDeformation", PhaseTracer::ActionMethod::PathDeformation)
      .value("All", PhaseTracer::ActionMethod::All);

  py::enum_<nlopt::algorithm>(m, "NLoptAlgorithm", "Minimisation algorithm used by PhaseFinder")
      .value("LN_SBPLX", nlopt::LN_SBPLX)
      .value("LN_NELDERMEAD", nlopt::LN_NELDERMEAD)
      .value("LN_COBYLA", nlopt::LN_COBYLA)
      .value("LN_BOBYQA", nlopt::LN_BOBYQA)
      .value("LN_NEWUOA", nlopt::LN_NEWUOA)
      .value("LN_PRAXIS", nlopt::LN_PRAXIS)
      .value("GN_DIRECT", nlopt::GN_DIRECT)
      .value("GN_DIRECT_L", nlopt::GN_DIRECT_L)
      .value("GN_CRS2_LM", nlopt::GN_CRS2_LM)
      .value("GN_ISRES", nlopt::GN_ISRES)
      .value("GN_ESCH", nlopt::GN_ESCH);

  py::enum_<PhaseTracer::PrintSettings>(m, "PrintSettings", "How much of a milestone ThermoFinder prints")
      .value("MINIMAL", PhaseTracer::MINIMAL)
      .value("STANDARD", PhaseTracer::STANDARD)
      .value("VERBOSE", PhaseTracer::VERBOSE);

  py::enum_<PhaseTracer::ValidateMethod>(m, "ValidateMethod", "Default screening of transitions in ThermoFinder")
      .value("TEMP", PhaseTracer::TEMP)
      .value("VEV", PhaseTracer::VEV)
      .value("NONE", PhaseTracer::NONE);

  py::enum_<PhaseTracer::MilestoneType>(m, "MilestoneType", "Transition milestone")
      .value("PERCOLATION", PhaseTracer::PERCOLATION)
      .value("NUCLEATION", PhaseTracer::NUCLEATION)
      .value("COMPLETION", PhaseTracer::COMPLETION)
      .value("ONSET", PhaseTracer::ONSET);

  py::enum_<PhaseTracer::GravWaveMethod>(m, "GravWaveMethod", "Backend for the gravitational wave spectrum")
      .value("FitFormulae", PhaseTracer::GravWaveMethod::FitFormulae)
      .value("SoundShell", PhaseTracer::GravWaveMethod::SoundShell);

  // ------------------------------------------------------------------ sub-structs
  // List fields (lower_bounds, upper_bounds, guess_points) are converted by copy: assign a whole
  // new list, since appending to the returned list does not change the config.

  py::class_<PhaseTracer::PhaseFinderConfig>(m, "PhaseFinderConfig", "Settings for PhaseFinder. Empty lower_bounds/upper_bounds keep the PhaseFinder defaults (+-1600 per field).")
      .def(py::init<>())
      .def_readwrite("x_abs_identical", &PhaseTracer::PhaseFinderConfig::x_abs_identical)
      .def_readwrite("x_rel_identical", &PhaseTracer::PhaseFinderConfig::x_rel_identical)
      .def_readwrite("x_abs_jump", &PhaseTracer::PhaseFinderConfig::x_abs_jump)
      .def_readwrite("x_rel_jump", &PhaseTracer::PhaseFinderConfig::x_rel_jump)
      .def_readwrite("find_min_x_tol_rel", &PhaseTracer::PhaseFinderConfig::find_min_x_tol_rel)
      .def_readwrite("find_min_x_tol_abs", &PhaseTracer::PhaseFinderConfig::find_min_x_tol_abs)
      .def_readwrite("find_min_algorithm", &PhaseTracer::PhaseFinderConfig::find_min_algorithm)
      .def_readwrite("find_min_max_f_eval", &PhaseTracer::PhaseFinderConfig::find_min_max_f_eval)
      .def_readwrite("find_min_min_step", &PhaseTracer::PhaseFinderConfig::find_min_min_step)
      .def_readwrite("find_min_max_time", &PhaseTracer::PhaseFinderConfig::find_min_max_time)
      .def_readwrite("find_min_trace_abs_step", &PhaseTracer::PhaseFinderConfig::find_min_trace_abs_step)
      .def_readwrite("find_min_locate_abs_step", &PhaseTracer::PhaseFinderConfig::find_min_locate_abs_step)
      .def_readwrite("n_test_points", &PhaseTracer::PhaseFinderConfig::n_test_points)
      .def_readwrite("lower_bounds", &PhaseTracer::PhaseFinderConfig::lower_bounds)
      .def_readwrite("upper_bounds", &PhaseTracer::PhaseFinderConfig::upper_bounds)
      .def_readwrite("t_low", &PhaseTracer::PhaseFinderConfig::t_low)
      .def_readwrite("t_high", &PhaseTracer::PhaseFinderConfig::t_high)
      .def_readwrite("dt_start_rel", &PhaseTracer::PhaseFinderConfig::dt_start_rel)
      .def_readwrite("dt_min_rel_split_phase", &PhaseTracer::PhaseFinderConfig::dt_min_rel_split_phase)
      .def_readwrite("t_jump_rel", &PhaseTracer::PhaseFinderConfig::t_jump_rel)
      .def_readwrite("dt_max_abs", &PhaseTracer::PhaseFinderConfig::dt_max_abs)
      .def_readwrite("dt_max_rel", &PhaseTracer::PhaseFinderConfig::dt_max_rel)
      .def_readwrite("dt_min_rel", &PhaseTracer::PhaseFinderConfig::dt_min_rel)
      .def_readwrite("dt_min_abs", &PhaseTracer::PhaseFinderConfig::dt_min_abs)
      .def_readwrite("n_ew_scalars", &PhaseTracer::PhaseFinderConfig::n_ew_scalars)
      .def_readwrite("v", &PhaseTracer::PhaseFinderConfig::v)
      .def_readwrite("check_vacuum_at_low", &PhaseTracer::PhaseFinderConfig::check_vacuum_at_low)
      .def_readwrite("check_vacuum_at_high", &PhaseTracer::PhaseFinderConfig::check_vacuum_at_high)
      .def_readwrite("check_dx_min_dt", &PhaseTracer::PhaseFinderConfig::check_dx_min_dt)
      .def_readwrite("hessian_singular_rel_tol", &PhaseTracer::PhaseFinderConfig::hessian_singular_rel_tol)
      .def_readwrite("linear_algebra_rel_tol", &PhaseTracer::PhaseFinderConfig::linear_algebra_rel_tol)
      .def_readwrite("check_hessian_singular", &PhaseTracer::PhaseFinderConfig::check_hessian_singular)
      .def_readwrite("hessian_eig_max_rel_change", &PhaseTracer::PhaseFinderConfig::hessian_eig_max_rel_change)
      .def_readwrite("check_midpoint_hessian", &PhaseTracer::PhaseFinderConfig::check_midpoint_hessian)
      .def_readwrite("seed", &PhaseTracer::PhaseFinderConfig::seed)
      .def_readwrite("trace_max_iter", &PhaseTracer::PhaseFinderConfig::trace_max_iter)
      .def_readwrite("guess_points", &PhaseTracer::PhaseFinderConfig::guess_points)
      .def_readwrite("check_merge_phase_gaps", &PhaseTracer::PhaseFinderConfig::check_merge_phase_gaps)
      .def_readwrite("dt_merge_phases", &PhaseTracer::PhaseFinderConfig::dt_merge_phases)
      .def_readwrite("dx_merge_phases", &PhaseTracer::PhaseFinderConfig::dx_merge_phases)
      .def("__repr__", &repr_from_stream<PhaseTracer::PhaseFinderConfig>)
      .def("__copy__", [](const PhaseTracer::PhaseFinderConfig &c) { return PhaseTracer::PhaseFinderConfig(c); })
      .def("__deepcopy__", [](const PhaseTracer::PhaseFinderConfig &c, py::dict) { return PhaseTracer::PhaseFinderConfig(c); }, py::arg("memo"));

  py::class_<PhaseTracer::TransitionFinderConfig>(m, "TransitionFinderConfig", "Settings for TransitionFinder")
      .def(py::init<>())
      .def_readwrite("n_ew_scalars", &PhaseTracer::TransitionFinderConfig::n_ew_scalars)
      .def_readwrite("separation", &PhaseTracer::TransitionFinderConfig::separation)
      .def_readwrite("assume_only_one_transition", &PhaseTracer::TransitionFinderConfig::assume_only_one_transition)
      .def_readwrite("TC_tol_rel", &PhaseTracer::TransitionFinderConfig::TC_tol_rel)
      .def_readwrite("max_iter", &PhaseTracer::TransitionFinderConfig::max_iter)
      .def_readwrite("change_rel_tol", &PhaseTracer::TransitionFinderConfig::change_rel_tol)
      .def_readwrite("change_abs_tol", &PhaseTracer::TransitionFinderConfig::change_abs_tol)
      .def_readwrite("Tnuc_step", &PhaseTracer::TransitionFinderConfig::Tnuc_step)
      .def_readwrite("Tnuc_tol_rel", &PhaseTracer::TransitionFinderConfig::Tnuc_tol_rel)
      .def_readwrite("check_subcritical_transitions", &PhaseTracer::TransitionFinderConfig::check_subcritical_transitions)
      .def("__repr__", &repr_from_stream<PhaseTracer::TransitionFinderConfig>)
      .def("__copy__", [](const PhaseTracer::TransitionFinderConfig &c) { return PhaseTracer::TransitionFinderConfig(c); })
      .def("__deepcopy__", [](const PhaseTracer::TransitionFinderConfig &c, py::dict) { return PhaseTracer::TransitionFinderConfig(c); }, py::arg("memo"));

  py::class_<PhaseTracer::ActionCalculatorConfig>(m, "ActionCalculatorConfig", "Settings for ActionCalculator (bounce action)")
      .def(py::init<>())
      .def_readwrite("method", &PhaseTracer::ActionCalculatorConfig::method)
      .def_readwrite("num_dims", &PhaseTracer::ActionCalculatorConfig::num_dims)
      .def_readwrite("BP_use_perturbative", &PhaseTracer::ActionCalculatorConfig::BP_use_perturbative)
      .def_readwrite("BP_initial_step_size", &PhaseTracer::ActionCalculatorConfig::BP_initial_step_size)
      .def_readwrite("BP_interpolation_points_fraction", &PhaseTracer::ActionCalculatorConfig::BP_interpolation_points_fraction)
      .def_readwrite("PD_xtol", &PhaseTracer::ActionCalculatorConfig::PD_xtol)
      .def_readwrite("PD_phitol", &PhaseTracer::ActionCalculatorConfig::PD_phitol)
      .def_readwrite("PD_thin_cutoff", &PhaseTracer::ActionCalculatorConfig::PD_thin_cutoff)
      .def_readwrite("PD_npoints", &PhaseTracer::ActionCalculatorConfig::PD_npoints)
      .def_readwrite("PD_rmin", &PhaseTracer::ActionCalculatorConfig::PD_rmin)
      .def_readwrite("PD_rmax", &PhaseTracer::ActionCalculatorConfig::PD_rmax)
      .def_readwrite("PD_max_iter", &PhaseTracer::ActionCalculatorConfig::PD_max_iter)
      .def_readwrite("PD_nb", &PhaseTracer::ActionCalculatorConfig::PD_nb)
      .def_readwrite("PD_kb", &PhaseTracer::ActionCalculatorConfig::PD_kb)
      .def_readwrite("PD_save_all_steps", &PhaseTracer::ActionCalculatorConfig::PD_save_all_steps)
      .def_readwrite("PD_v2min", &PhaseTracer::ActionCalculatorConfig::PD_v2min)
      .def_readwrite("PD_step_maxiter", &PhaseTracer::ActionCalculatorConfig::PD_step_maxiter)
      .def_readwrite("PD_path_maxiter", &PhaseTracer::ActionCalculatorConfig::PD_path_maxiter)
      .def_readwrite("PD_V_spline_samples", &PhaseTracer::ActionCalculatorConfig::PD_V_spline_samples)
      .def_readwrite("PD_extend_to_minima", &PhaseTracer::ActionCalculatorConfig::PD_extend_to_minima)
      .def_readwrite("PD_deformation_npoints", &PhaseTracer::ActionCalculatorConfig::PD_deformation_npoints)
      .def_readwrite("PD_fRatioConv", &PhaseTracer::ActionCalculatorConfig::PD_fRatioConv)
      .def_readwrite("PD_warm_start_fRatioConv", &PhaseTracer::ActionCalculatorConfig::PD_warm_start_fRatioConv)
      .def("__repr__", &repr_from_stream<PhaseTracer::ActionCalculatorConfig>)
      .def("__copy__", [](const PhaseTracer::ActionCalculatorConfig &c) { return PhaseTracer::ActionCalculatorConfig(c); })
      .def("__deepcopy__", [](const PhaseTracer::ActionCalculatorConfig &c, py::dict) { return PhaseTracer::ActionCalculatorConfig(c); }, py::arg("memo"));

  py::class_<PhaseTracer::ThermoFinderConfig>(m, "ThermoFinderConfig", "Settings for ThermoFinder, including those it forwards to FalseVacuumDecayRate, EquationOfState and FriedmannEvolution")
      .def(py::init<>())
      .def_readwrite("onset_print_setting", &PhaseTracer::ThermoFinderConfig::onset_print_setting)
      .def_readwrite("percolation_print_setting", &PhaseTracer::ThermoFinderConfig::percolation_print_setting)
      .def_readwrite("nucleation_print_setting", &PhaseTracer::ThermoFinderConfig::nucleation_print_setting)
      .def_readwrite("completion_print_setting", &PhaseTracer::ThermoFinderConfig::completion_print_setting)
      .def_readwrite("update_percolation_temperature", &PhaseTracer::ThermoFinderConfig::update_percolation_temperature)
      .def_readwrite("compute_profiles", &PhaseTracer::ThermoFinderConfig::compute_profiles)
      .def_readwrite("n_temp_profiles", &PhaseTracer::ThermoFinderConfig::n_temp_profiles)
      .def_readwrite("vw", &PhaseTracer::ThermoFinderConfig::vw)
      .def_readwrite("temperature_threshold", &PhaseTracer::ThermoFinderConfig::temperature_threshold)
      .def_readwrite("vev_threshold", &PhaseTracer::ThermoFinderConfig::vev_threshold)
      .def_readwrite("default_validation_method", &PhaseTracer::ThermoFinderConfig::default_validation_method)
      .def_readwrite("transition_filter", &PhaseTracer::ThermoFinderConfig::transition_filter, "Callable[[list[Transition]], list[Transition]] selecting the transitions to analyse; None for default_validation_method")
      .def_readwrite("action_spline_evaluations", &PhaseTracer::ThermoFinderConfig::action_spline_evaluations)
      .def_readwrite("warm_start_chunk_size", &PhaseTracer::ThermoFinderConfig::warm_start_chunk_size)
      .def_readwrite("action_smoothing_window", &PhaseTracer::ThermoFinderConfig::action_smoothing_window)
      .def_readwrite("action_smoothing_order", &PhaseTracer::ThermoFinderConfig::action_smoothing_order)
      .def_readwrite("action_laurent_tail", &PhaseTracer::ThermoFinderConfig::action_laurent_tail)
      .def_readwrite("prefactor_function", &PhaseTracer::ThermoFinderConfig::prefactor_function, "Callable[[float T, float S_over_T, ActionResult], float] for the decay-rate prefactor; None for the default")
      .def_readwrite("eos_spline_evaluations", &PhaseTracer::ThermoFinderConfig::eos_spline_evaluations)
      .def_readwrite("eos_background_dof", &PhaseTracer::ThermoFinderConfig::eos_background_dof)
      .def_readwrite("percolation_target", &PhaseTracer::ThermoFinderConfig::percolation_target)
      .def_readwrite("completion_target", &PhaseTracer::ThermoFinderConfig::completion_target)
      .def_readwrite("onset_target", &PhaseTracer::ThermoFinderConfig::onset_target)
      .def_readwrite("nucleation_target", &PhaseTracer::ThermoFinderConfig::nucleation_target)
      .def_readwrite("use_bag_dtdT", &PhaseTracer::ThermoFinderConfig::use_bag_dtdT)
      .def_readwrite("temperature_abs_tol", &PhaseTracer::ThermoFinderConfig::temperature_abs_tol)
      .def("__repr__", &repr_from_stream<PhaseTracer::ThermoFinderConfig>)
      .def("__copy__", [](const PhaseTracer::ThermoFinderConfig &c) { return PhaseTracer::ThermoFinderConfig(c); })
      .def("__deepcopy__", [](const PhaseTracer::ThermoFinderConfig &c, py::dict) { return PhaseTracer::ThermoFinderConfig(c); }, py::arg("memo"));

  py::class_<PhaseTracer::GravWaveConfig>(m, "GravWaveConfig", "Settings for GravWaveCalculator. g_eff, h_eff and vw of 0 mean: use the milestone value.")
      .def(py::init<>())
      .def_readwrite("gw_method", &PhaseTracer::GravWaveConfig::gw_method)
      .def_readwrite("use_legacy_gw_methods", &PhaseTracer::GravWaveConfig::use_legacy_gw_methods)
      .def_readwrite("min_kRs_value", &PhaseTracer::GravWaveConfig::min_kRs_value)
      .def_readwrite("max_kRs_value", &PhaseTracer::GravWaveConfig::max_kRs_value)
      .def_readwrite("n_kRs_value", &PhaseTracer::GravWaveConfig::n_kRs_value)
      .def_readwrite("min_frequency", &PhaseTracer::GravWaveConfig::min_frequency)
      .def_readwrite("max_frequency", &PhaseTracer::GravWaveConfig::max_frequency)
      .def_readwrite("num_frequency", &PhaseTracer::GravWaveConfig::num_frequency)
      .def_readwrite("num_frequency_ssm", &PhaseTracer::GravWaveConfig::num_frequency_ssm)
      .def_readwrite("T_threshold_bubble_collision", &PhaseTracer::GravWaveConfig::T_threshold_bubble_collision)
      .def_readwrite("h_dVdT", &PhaseTracer::GravWaveConfig::h_dVdT)
      .def_readwrite("h_dSdT", &PhaseTracer::GravWaveConfig::h_dSdT)
      .def_readwrite("np_dSdT", &PhaseTracer::GravWaveConfig::np_dSdT)
      .def_readwrite("default_milestone", &PhaseTracer::GravWaveConfig::default_milestone)
      .def_readwrite("include_col_and_turb_in_ssm", &PhaseTracer::GravWaveConfig::include_col_and_turb_in_ssm)
      .def_readwrite("g_0", &PhaseTracer::GravWaveConfig::g_0)
      .def_readwrite("h_0", &PhaseTracer::GravWaveConfig::h_0)
      .def_readwrite("g_eff", &PhaseTracer::GravWaveConfig::g_eff)
      .def_readwrite("h_eff", &PhaseTracer::GravWaveConfig::h_eff)
      .def_readwrite("omega_hsq_neutrino", &PhaseTracer::GravWaveConfig::omega_hsq_neutrino)
      .def_readwrite("D", &PhaseTracer::GravWaveConfig::D)
      .def_readwrite("vw", &PhaseTracer::GravWaveConfig::vw)
      .def_readwrite("epsilon", &PhaseTracer::GravWaveConfig::epsilon)
      .def_readwrite("run_time_LISA", &PhaseTracer::GravWaveConfig::run_time_LISA)
      .def_readwrite("run_time_Taiji", &PhaseTracer::GravWaveConfig::run_time_Taiji)
      .def_readwrite("use_legacy_LISA_noise", &PhaseTracer::GravWaveConfig::use_legacy_LISA_noise)
      .def_readwrite("SNR_f_min", &PhaseTracer::GravWaveConfig::SNR_f_min)
      .def_readwrite("SNR_f_max", &PhaseTracer::GravWaveConfig::SNR_f_max)
      .def_readwrite("SNR_steps_per_decade", &PhaseTracer::GravWaveConfig::SNR_steps_per_decade)
      .def("__repr__", &repr_from_stream<PhaseTracer::GravWaveConfig>)
      .def("__copy__", [](const PhaseTracer::GravWaveConfig &c) { return PhaseTracer::GravWaveConfig(c); })
      .def("__deepcopy__", [](const PhaseTracer::GravWaveConfig &c, py::dict) { return PhaseTracer::GravWaveConfig(c); }, py::arg("memo"));

  py::class_<PhaseTracer::PipelineConfig>(m, "PipelineConfig", "Settings for the pipeline as a whole")
      .def(py::init<>())
      .def_readwrite("stop_after", &PhaseTracer::PipelineConfig::stop_after, "Last stage to run")
      .def_readwrite("throw_on_error", &PhaseTracer::PipelineConfig::throw_on_error, "Raise RunnerError instead of returning a failed RunStatus")
      .def_readwrite("log_level", &PhaseTracer::PipelineConfig::log_level, "If set (a LogLevel), the global log level is changed before running")
      .def_readwrite("to_print", &PhaseTracer::PipelineConfig::to_print, "Print each stage after it is evaluated")
      .def("__repr__", &repr_from_stream<PhaseTracer::PipelineConfig>)
      .def("__copy__", [](const PhaseTracer::PipelineConfig &c) { return PhaseTracer::PipelineConfig(c); })
      .def("__deepcopy__", [](const PhaseTracer::PipelineConfig &c, py::dict) { return PhaseTracer::PipelineConfig(c); }, py::arg("memo"));

  // ------------------------------------------------------------------ Config

  py::class_<PhaseTracer::Config>(
      m, "Config",
      "All settings for a Runner, grouped per class. Defaults equal the class defaults.\n"
      "Sub-structs are returned by reference, so config.phase_finder.seed = 0 edits in place.")
      .def(py::init<>())
      .def_readwrite("phase_finder", &PhaseTracer::Config::phase_finder)
      .def_readwrite("transition_finder", &PhaseTracer::Config::transition_finder)
      .def_readwrite("action", &PhaseTracer::Config::action)
      .def_readwrite("thermo_finder", &PhaseTracer::Config::thermo_finder)
      .def_readwrite("gravwave", &PhaseTracer::Config::gravwave)
      .def_readwrite("pipeline", &PhaseTracer::Config::pipeline)
      .def("validate", &PhaseTracer::Config::validate,
           "Sanity checks that need no model; returns a RunStatus listing every problem")
      .def("__repr__", &repr_from_stream<PhaseTracer::Config>)
      .def("__copy__", [](const PhaseTracer::Config &c) { return PhaseTracer::Config(c); })
      .def("__deepcopy__", [](const PhaseTracer::Config &c, py::dict) { return PhaseTracer::Config(c); },
           py::arg("memo"));
}
