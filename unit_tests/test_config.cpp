#include <sstream>
#include <string>

#include "catch/catch.hpp"
#include "models/1D_test_model.hpp"
#include "config.hpp"
#include "logger.hpp"

// Config defaults must match the PROPERTY defaults of the classes they mirror
#define CHECK_DEFAULT(obj, cfg, field) CHECK((cfg).field == (obj).get_##field())

TEST_CASE("Config defaults match class defaults", "[Config]") {

  LOGGER(fatal);

  EffectivePotential::OneDimModel model;
  const PhaseTracer::Config config;

  PhaseTracer::PhaseFinder pf(model);
  PhaseTracer::ActionCalculator ac(pf);
  PhaseTracer::TransitionFinder tf(pf);
  PhaseTracer::ThermoFinder tm(tf, ac);
  PhaseTracer::GravWaveCalculator gw(tm);

  SECTION("PhaseFinder") {
    const auto& c = config.phase_finder;
    CHECK_DEFAULT(pf, c, x_abs_identical);
    CHECK_DEFAULT(pf, c, x_rel_identical);
    CHECK_DEFAULT(pf, c, x_abs_jump);
    CHECK_DEFAULT(pf, c, x_rel_jump);
    CHECK_DEFAULT(pf, c, find_min_x_tol_rel);
    CHECK_DEFAULT(pf, c, find_min_x_tol_abs);
    CHECK_DEFAULT(pf, c, find_min_algorithm);
    CHECK_DEFAULT(pf, c, find_min_max_f_eval);
    CHECK_DEFAULT(pf, c, find_min_min_step);
    CHECK_DEFAULT(pf, c, find_min_max_time);
    CHECK_DEFAULT(pf, c, find_min_trace_abs_step);
    CHECK_DEFAULT(pf, c, find_min_locate_abs_step);
    CHECK_DEFAULT(pf, c, n_test_points);
    CHECK_DEFAULT(pf, c, t_low);
    CHECK_DEFAULT(pf, c, t_high);
    CHECK_DEFAULT(pf, c, dt_start_rel);
    CHECK_DEFAULT(pf, c, dt_min_rel_split_phase);
    CHECK_DEFAULT(pf, c, t_jump_rel);
    CHECK_DEFAULT(pf, c, dt_max_abs);
    CHECK_DEFAULT(pf, c, dt_max_rel);
    CHECK_DEFAULT(pf, c, dt_min_rel);
    CHECK_DEFAULT(pf, c, dt_min_abs);
    CHECK_DEFAULT(pf, c, n_ew_scalars);
    CHECK_DEFAULT(pf, c, v);
    CHECK_DEFAULT(pf, c, check_vacuum_at_low);
    CHECK_DEFAULT(pf, c, check_vacuum_at_high);
    CHECK_DEFAULT(pf, c, check_dx_min_dt);
    CHECK_DEFAULT(pf, c, hessian_singular_rel_tol);
    CHECK_DEFAULT(pf, c, linear_algebra_rel_tol);
    CHECK_DEFAULT(pf, c, seed);
    CHECK_DEFAULT(pf, c, check_hessian_singular);
    CHECK_DEFAULT(pf, c, hessian_eig_max_rel_change);
    CHECK_DEFAULT(pf, c, check_midpoint_hessian);
    CHECK_DEFAULT(pf, c, trace_max_iter);
    CHECK(c.guess_points.size() == pf.get_guess_points().size());
    CHECK_DEFAULT(pf, c, check_merge_phase_gaps);
    CHECK_DEFAULT(pf, c, dt_merge_phases);
    CHECK_DEFAULT(pf, c, dx_merge_phases);
    // Bounds are filled by the PhaseFinder constructor; empty in Config means "keep them"
    CHECK(c.lower_bounds.empty());
    CHECK(c.upper_bounds.empty());
  }

  SECTION("TransitionFinder") {
    const auto& c = config.transition_finder;
    CHECK_DEFAULT(tf, c, n_ew_scalars);
    CHECK_DEFAULT(tf, c, separation);
    CHECK_DEFAULT(tf, c, assume_only_one_transition);
    CHECK_DEFAULT(tf, c, TC_tol_rel);
    CHECK_DEFAULT(tf, c, max_iter);
    CHECK_DEFAULT(tf, c, change_rel_tol);
    CHECK_DEFAULT(tf, c, change_abs_tol);
    CHECK_DEFAULT(tf, c, Tnuc_step);
    CHECK_DEFAULT(tf, c, Tnuc_tol_rel);
    CHECK_DEFAULT(tf, c, check_subcritical_transitions);
  }

  SECTION("ActionCalculator") {
    const auto& c = config.action;
    CHECK(c.method == ac.get_action_calculator());
    CHECK_DEFAULT(ac, c, num_dims);
    CHECK_DEFAULT(ac, c, BP_use_perturbative);
    CHECK_DEFAULT(ac, c, BP_initial_step_size);
    CHECK_DEFAULT(ac, c, BP_interpolation_points_fraction);
    CHECK_DEFAULT(ac, c, PD_xtol);
    CHECK_DEFAULT(ac, c, PD_phitol);
    CHECK_DEFAULT(ac, c, PD_thin_cutoff);
    CHECK_DEFAULT(ac, c, PD_npoints);
    CHECK_DEFAULT(ac, c, PD_rmin);
    CHECK_DEFAULT(ac, c, PD_rmax);
    CHECK_DEFAULT(ac, c, PD_max_iter);
    CHECK_DEFAULT(ac, c, PD_nb);
    CHECK_DEFAULT(ac, c, PD_kb);
    CHECK_DEFAULT(ac, c, PD_save_all_steps);
    CHECK_DEFAULT(ac, c, PD_v2min);
    CHECK_DEFAULT(ac, c, PD_step_maxiter);
    CHECK_DEFAULT(ac, c, PD_path_maxiter);
    CHECK_DEFAULT(ac, c, PD_V_spline_samples);
    CHECK_DEFAULT(ac, c, PD_extend_to_minima);
    CHECK_DEFAULT(ac, c, PD_deformation_npoints);
    CHECK_DEFAULT(ac, c, PD_fRatioConv);
    CHECK_DEFAULT(ac, c, PD_warm_start_fRatioConv);
  }

  SECTION("ThermoFinder") {
    const auto& c = config.thermo_finder;
    CHECK_DEFAULT(tm, c, onset_print_setting);
    CHECK_DEFAULT(tm, c, percolation_print_setting);
    CHECK_DEFAULT(tm, c, nucleation_print_setting);
    CHECK_DEFAULT(tm, c, completion_print_setting);
    CHECK_DEFAULT(tm, c, update_percolation_temperature);
    CHECK_DEFAULT(tm, c, compute_profiles);
    CHECK_DEFAULT(tm, c, n_temp_profiles);
    CHECK_DEFAULT(tm, c, vw);
    CHECK_DEFAULT(tm, c, temperature_threshold);
    CHECK_DEFAULT(tm, c, vev_threshold);
    CHECK_DEFAULT(tm, c, default_validation_method);
    CHECK(static_cast<bool>(c.transition_filter) == static_cast<bool>(tm.get_transition_filter()));
    CHECK_DEFAULT(tm, c, action_spline_evaluations);
    CHECK_DEFAULT(tm, c, warm_start_chunk_size);
    CHECK_DEFAULT(tm, c, action_smoothing_window);
    CHECK_DEFAULT(tm, c, action_smoothing_order);
    CHECK_DEFAULT(tm, c, action_laurent_tail);
    CHECK_DEFAULT(tm, c, eos_spline_evaluations);
    CHECK_DEFAULT(tm, c, eos_background_dof);
    CHECK_DEFAULT(tm, c, percolation_target);
    CHECK_DEFAULT(tm, c, completion_target);
    CHECK_DEFAULT(tm, c, onset_target);
    CHECK_DEFAULT(tm, c, nucleation_target);
    CHECK_DEFAULT(tm, c, use_bag_dtdT);
    CHECK_DEFAULT(tm, c, temperature_abs_tol);
  }

  SECTION("GravWaveCalculator") {
    const auto& c = config.gravwave;
    CHECK_DEFAULT(gw, c, gw_method);
    CHECK_DEFAULT(gw, c, use_legacy_gw_methods);
    CHECK_DEFAULT(gw, c, min_kRs_value);
    CHECK_DEFAULT(gw, c, max_kRs_value);
    CHECK_DEFAULT(gw, c, n_kRs_value);
    CHECK_DEFAULT(gw, c, min_frequency);
    CHECK_DEFAULT(gw, c, max_frequency);
    CHECK_DEFAULT(gw, c, num_frequency);
    CHECK_DEFAULT(gw, c, num_frequency_ssm);
    CHECK_DEFAULT(gw, c, T_threshold_bubble_collision);
    CHECK_DEFAULT(gw, c, h_dVdT);
    CHECK_DEFAULT(gw, c, h_dSdT);
    CHECK_DEFAULT(gw, c, np_dSdT);
    CHECK_DEFAULT(gw, c, default_milestone);
    CHECK_DEFAULT(gw, c, include_col_and_turb_in_ssm);
    CHECK_DEFAULT(gw, c, g_0);
    CHECK_DEFAULT(gw, c, h_0);
    CHECK_DEFAULT(gw, c, g_eff);
    CHECK_DEFAULT(gw, c, h_eff);
    CHECK_DEFAULT(gw, c, omega_hsq_neutrino);
    CHECK_DEFAULT(gw, c, D);
    CHECK_DEFAULT(gw, c, vw);
    CHECK_DEFAULT(gw, c, epsilon);
    CHECK_DEFAULT(gw, c, run_time_LISA);
    CHECK_DEFAULT(gw, c, run_time_Taiji);
    CHECK_DEFAULT(gw, c, use_legacy_LISA_noise);
    CHECK_DEFAULT(gw, c, SNR_f_min);
    CHECK_DEFAULT(gw, c, SNR_f_max);
    CHECK_DEFAULT(gw, c, SNR_steps_per_decade);
  }
}

#undef CHECK_DEFAULT

TEST_CASE("apply() pushes Config settings onto objects", "[Config]") {

  LOGGER(fatal);

  EffectivePotential::OneDimModel model;
  PhaseTracer::Config config;

  PhaseTracer::PhaseFinder pf(model);
  PhaseTracer::ActionCalculator ac(pf);
  PhaseTracer::TransitionFinder tf(pf);
  PhaseTracer::ThermoFinder tm(tf, ac);
  PhaseTracer::GravWaveCalculator gw(tm);

  SECTION("PhaseFinder") {
    const auto default_lower = pf.get_lower_bounds();
    config.phase_finder.seed = 7;
    config.phase_finder.t_high = 321.;
    config.phase_finder.check_hessian_singular = false;
    config.phase_finder.find_min_algorithm = nlopt::LN_COBYLA;
    config.phase_finder.upper_bounds = {500.};
    PhaseTracer::apply(config.phase_finder, pf);
    CHECK(pf.get_seed() == 7);
    CHECK(pf.get_t_high() == 321.);
    CHECK_FALSE(pf.get_check_hessian_singular());
    CHECK(pf.get_find_min_algorithm() == nlopt::LN_COBYLA);
    CHECK(pf.get_upper_bounds() == std::vector<double>{500.});
    // empty lower_bounds leaves the constructor default untouched
    CHECK(pf.get_lower_bounds() == default_lower);
  }

  SECTION("TransitionFinder") {
    config.transition_finder.TC_tol_rel = 1e-9;
    config.transition_finder.assume_only_one_transition = false;
    PhaseTracer::apply(config.transition_finder, tf);
    CHECK(tf.get_TC_tol_rel() == 1e-9);
    CHECK_FALSE(tf.get_assume_only_one_transition());
  }

  SECTION("ActionCalculator") {
    config.action.PD_xtol = 1e-8;
    config.action.PD_deformation_npoints = 150;
    PhaseTracer::apply(config.action, ac);
    CHECK(ac.get_PD_xtol() == 1e-8);
    CHECK(ac.get_PD_deformation_npoints() == 150);
  }

  SECTION("ThermoFinder") {
    config.thermo_finder.percolation_target = 1 - 0.28957;
    config.thermo_finder.action_laurent_tail = true;
    config.thermo_finder.percolation_print_setting = PhaseTracer::VERBOSE;
    config.thermo_finder.transition_filter = [](const std::vector<PhaseTracer::Transition>& t) { return t; };
    PhaseTracer::apply(config.thermo_finder, tm);
    CHECK(tm.get_percolation_target() == 1 - 0.28957);
    CHECK(tm.get_action_laurent_tail());
    CHECK(tm.get_percolation_print_setting() == PhaseTracer::VERBOSE);
    CHECK(static_cast<bool>(tm.get_transition_filter()));
  }

  SECTION("GravWaveCalculator") {
    config.gravwave.num_frequency = 123;
    config.gravwave.vw = 0.9;
    config.gravwave.default_milestone = PhaseTracer::NUCLEATION;
    PhaseTracer::apply(config.gravwave, gw);
    CHECK(gw.get_num_frequency() == 123);
    CHECK(gw.get_vw() == 0.9);
    CHECK(gw.get_default_milestone() == PhaseTracer::NUCLEATION);
  }

#ifndef BUILD_WITH_HG
  SECTION("SoundShell without HydroGrav throws") {
    config.gravwave.gw_method = PhaseTracer::GravWaveMethod::SoundShell;
    CHECK_THROWS(PhaseTracer::apply(config.gravwave, gw));
  }
#endif
}

TEST_CASE("Config::validate", "[Config]") {

  PhaseTracer::Config config;

  SECTION("Defaults are valid") {
    const auto status = config.validate();
    CHECK(status.ok());
    CHECK(static_cast<bool>(status));
    CHECK(status.message.empty());
  }

  SECTION("Bad frequency grid is reported") {
    config.gravwave.min_frequency = 1.;
    config.gravwave.max_frequency = 1e-3;
    const auto status = config.validate();
    CHECK_FALSE(status.ok());
    CHECK(status.code == PhaseTracer::StatusCode::InvalidConfig);
    CHECK(status.stage == PhaseTracer::Stage::Config);
    CHECK(status.message.find("min_frequency") != std::string::npos);
  }

  SECTION("All problems are listed") {
    config.phase_finder.t_low = 100.;
    config.phase_finder.t_high = 50.;
    config.thermo_finder.action_smoothing_window = 4;
    config.thermo_finder.percolation_target = 0.;
    const auto status = config.validate();
    CHECK_FALSE(status.ok());
    CHECK(status.message.find("t_low") != std::string::npos);
    CHECK(status.message.find("action_smoothing_window") != std::string::npos);
    CHECK(status.message.find("percolation_target") != std::string::npos);
  }

  SECTION("Mismatched bounds are reported") {
    config.phase_finder.lower_bounds = {-10., -10.};
    config.phase_finder.upper_bounds = {10.};
    CHECK_FALSE(config.validate().ok());
  }

  SECTION("stop_after must be a computation stage") {
    config.pipeline.stop_after = PhaseTracer::Stage::Config;
    CHECK_FALSE(config.validate().ok());
  }
}

TEST_CASE("RunStatus printing and RunnerError", "[Config]") {

  PhaseTracer::RunStatus status;
  status.code = PhaseTracer::StatusCode::NoTransitions;
  status.stage = PhaseTracer::Stage::TransitionFinder;
  status.message = "no transitions found";
  status.warnings = {"something odd"};

  std::ostringstream ss;
  ss << status;
  CHECK(ss.str() == "[NoTransitions] at stage TransitionFinder: no transitions found\n  warning: something odd\n");

  try {
    throw PhaseTracer::RunnerError(status);
  } catch (const PhaseTracer::RunnerError& e) {
    CHECK(e.status().code == PhaseTracer::StatusCode::NoTransitions);
    CHECK(std::string(e.what()) == ss.str());
  }
}
