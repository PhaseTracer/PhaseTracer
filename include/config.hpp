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

#ifndef PHASETRACER_CONFIG_HPP_
#define PHASETRACER_CONFIG_HPP_

/**
 * Settings for every stage of the PhaseTracer pipeline, collected in one struct.
 *
 * Each sub-struct mirrors the PROPERTY settings of one class, with identical names,
 * types and defaults (kept in sync by unit_tests/test_config.cpp). Use apply() to
 * push a sub-struct onto the corresponding object.
 *
 *   auto config = PhaseTracer::Config();
 *   config.phase_finder.seed = 0;
 *   config.action.PD_xtol = 1e-8;
 *   config.thermo_finder.percolation_target = 0.71;
 */

#include <cmath>
#include <functional>
#include <optional>
#include <ostream>
#include <vector>

#include <Eigen/Core>
#include <boost/cstdint.hpp>
#include "nlopt.hpp"

#include "logger.hpp"
#include "run_status.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "action_calculator.hpp"
#include "thermo_finder.hpp"
#include "gravwave_calculator.hpp"

namespace PhaseTracer {

/** @brief Settings for PhaseFinder (include/phase_finder.hpp). */
struct PhaseFinderConfig
{
    double x_abs_identical = 1.;
    double x_rel_identical = 1.e-3;
    double x_abs_jump = 0.5;
    double x_rel_jump = 1.e-2;

    double find_min_x_tol_rel = 0.0001;
    double find_min_x_tol_abs = 0.0001;
    nlopt::algorithm find_min_algorithm = nlopt::LN_SBPLX;
    double find_min_max_f_eval = 1000000;
    double find_min_min_step = 1.e-4;
    double find_min_max_time = 5.;
    double find_min_trace_abs_step = 1.;
    double find_min_locate_abs_step = 1.;
    size_t n_test_points = 100;

    std::vector<double> lower_bounds = {};
    std::vector<double> upper_bounds = {};

    double t_low = 0.;
    double t_high = 1000.;
    double dt_start_rel = 0.01;
    double dt_min_rel_split_phase = 0.001;
    double t_jump_rel = 0.005;
    double dt_max_abs = 50.;
    double dt_max_rel = 0.25;
    double dt_min_rel = 1.e-7;
    double dt_min_abs = 1.e-10;

    size_t n_ew_scalars = 0;
    double v = 246.;

    bool check_vacuum_at_low = true;
    bool check_vacuum_at_high = true;
    bool check_dx_min_dt = true;
    double hessian_singular_rel_tol = 1.e-2;
    double linear_algebra_rel_tol = 1.e-3;
    bool check_hessian_singular = true;
    double hessian_eig_max_rel_change = 0.5;
    bool check_midpoint_hessian = true;

    int seed = -1;
    unsigned int trace_max_iter = 100000;
    std::vector<Eigen::VectorXd> guess_points = {};

    bool check_merge_phase_gaps = false;
    double dt_merge_phases = 5.;
    double dx_merge_phases = 70.;
};

/** @brief Settings for TransitionFinder (include/transition_finder.hpp). */
struct TransitionFinderConfig
{
    int n_ew_scalars = -1;
    double separation = 1.;
    bool assume_only_one_transition = true;
    double TC_tol_rel = 1.e-4;
    boost::uintmax_t max_iter = 100;
    double change_rel_tol = 1.e-3;
    double change_abs_tol = 1.e-3;
    double Tnuc_step = 1.;
    double Tnuc_tol_rel = 1.e-3;
    bool check_subcritical_transitions = false;
};

/** @brief Settings for ActionCalculator (include/action_calculator.hpp). */
struct ActionCalculatorConfig
{
    ActionMethod method = ActionMethod::PathDeformation;
    size_t num_dims = 3;

    bool BP_use_perturbative = false;
    double BP_initial_step_size = 1.e-2;
    double BP_interpolation_points_fraction = 1.0;

    double PD_xtol = 1e-4;
    double PD_phitol = 1e-4;
    double PD_thin_cutoff = .01;
    double PD_npoints = 500;
    double PD_rmin = 1e-4;
    double PD_rmax = 1e4;
    boost::uintmax_t PD_max_iter = 100;

    size_t PD_nb = 10;
    size_t PD_kb = 3;
    bool PD_save_all_steps = false;
    double PD_v2min = 0.0;
    size_t PD_step_maxiter = 500;
    size_t PD_path_maxiter = 20;
    size_t PD_V_spline_samples = 140;
    bool PD_extend_to_minima = true;
    size_t PD_deformation_npoints = 300;
    double PD_fRatioConv = .02;
    double PD_warm_start_fRatioConv = .002;
};

/** @brief Settings for ThermoFinder (include/thermo_finder.hpp), including those it forwards. */
struct ThermoFinderConfig
{
    using TransitionFilter = std::function<std::vector<Transition>(const std::vector<Transition>&)>;

    PrintSettings onset_print_setting = PrintSettings::MINIMAL;
    PrintSettings percolation_print_setting = PrintSettings::STANDARD;
    PrintSettings nucleation_print_setting = PrintSettings::STANDARD;
    PrintSettings completion_print_setting = PrintSettings::MINIMAL;

    bool update_percolation_temperature = false;
    bool compute_profiles = false;
    double n_temp_profiles = 250;
    double vw = 1/sqrt(3.0);

    double temperature_threshold = 1;
    double vev_threshold = 1;
    ValidateMethod default_validation_method = NONE;
    /** Custom transition filter; empty means use default_validation_method. */
    TransitionFilter transition_filter = {};

    int action_spline_evaluations = 50;
    int warm_start_chunk_size = 0;
    int action_smoothing_window = 0;
    int action_smoothing_order = 3;
    bool action_laurent_tail = false;
    /** Custom prefactor (e.g. BubbleDetPrefactor); empty means keep the default. */
    FalseVacuumDecayRate::PrefactorFunction prefactor_function = {};

    int eos_spline_evaluations = 250;
    double eos_background_dof = 0.0;

    double percolation_target = 0.71;
    double completion_target = 1e-6;
    double onset_target = 1 - 1e-6;
    double nucleation_target = 1.00;
    bool use_bag_dtdT = false;
    double temperature_abs_tol = 1e-6;
};

/** @brief Settings for GravWaveCalculator (include/gravwave_calculator.hpp). */
struct GravWaveConfig
{
    GravWaveMethod gw_method = GravWaveMethod::FitFormulae;
    bool use_legacy_gw_methods = false;

    double min_kRs_value = 1e-3;
    double max_kRs_value = 1e3;
    int n_kRs_value = 200;
    double min_frequency = 1e-4;
    double max_frequency = 1e1;
    int num_frequency = 500;
    int num_frequency_ssm = 100;
    double T_threshold_bubble_collision = 10;

    double h_dVdT = 1e-2;
    double h_dSdT = 1e-1;
    double np_dSdT = 5;

    MilestoneType default_milestone = MilestoneType::PERCOLATION;
    bool include_col_and_turb_in_ssm = false;

    double g_0 = 2.0;
    double h_0 = 3.91;
    double g_eff = 0.0;
    double h_eff = 0.0;
    double omega_hsq_neutrino = 2.473e-5;
    double D = 1.0;
    double vw = 0.0;
    double epsilon = 0.1;

    double run_time_LISA = 4;
    double run_time_Taiji = 3;
    bool use_legacy_LISA_noise = false;
    double SNR_f_min = 1e-5;
    double SNR_f_max = 1e-1;
    double SNR_steps_per_decade = 200;
};

/** @brief Settings for the pipeline as a whole. */
struct PipelineConfig
{
    /** Last stage to run */
    Stage stop_after = Stage::GravWave;
    /** Throw RunnerError instead of returning a failed RunStatus. */
    bool throw_on_error = false;
    /** If set, the global log level is changed to this before running. */
    std::optional<boost::log::trivial::severity_level> log_level = {};
    /** Whether or not to print classes during evaluation. */
    bool to_print = true;
};

/** @brief All settings for a full PhaseTracer run. */
struct Config
{
    PhaseFinderConfig phase_finder;
    TransitionFinderConfig transition_finder;
    ActionCalculatorConfig action;
    ThermoFinderConfig thermo_finder;
    GravWaveConfig gravwave;
    PipelineConfig pipeline;

    /**
     * @brief Sanity checks that need no model.
     * @return Success, or InvalidConfig at Stage::Config with every problem listed in the message.
     */
    RunStatus validate() const;
};

void apply(const PhaseFinderConfig& c, PhaseFinder& pf);
void apply(const TransitionFinderConfig& c, TransitionFinder& tf);
void apply(const ActionCalculatorConfig& c, ActionCalculator& ac);
void apply(const ThermoFinderConfig& c, ThermoFinder& tm);
void apply(const GravWaveConfig& c, GravWaveCalculator& gw);

std::ostream& operator<<(std::ostream& o, const PhaseFinderConfig& c);
std::ostream& operator<<(std::ostream& o, const TransitionFinderConfig& c);
std::ostream& operator<<(std::ostream& o, const ActionCalculatorConfig& c);
std::ostream& operator<<(std::ostream& o, const ThermoFinderConfig& c);
std::ostream& operator<<(std::ostream& o, const GravWaveConfig& c);
std::ostream& operator<<(std::ostream& o, const PipelineConfig& c);
std::ostream& operator<<(std::ostream& o, const Config& c);

} // namespace PhaseTracer

#endif // PHASETRACER_CONFIG_HPP_
