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

#include "config.hpp"

#include <sstream>
#include <string>

namespace PhaseTracer {

// =========================== apply ===========================

void apply(const PhaseFinderConfig& c, PhaseFinder& pf)
{
    pf.set_x_abs_identical(c.x_abs_identical);
    pf.set_x_rel_identical(c.x_rel_identical);
    pf.set_x_abs_jump(c.x_abs_jump);
    pf.set_x_rel_jump(c.x_rel_jump);

    pf.set_find_min_x_tol_rel(c.find_min_x_tol_rel);
    pf.set_find_min_x_tol_abs(c.find_min_x_tol_abs);
    pf.set_find_min_algorithm(c.find_min_algorithm);
    pf.set_find_min_max_f_eval(c.find_min_max_f_eval);
    pf.set_find_min_min_step(c.find_min_min_step);
    pf.set_find_min_max_time(c.find_min_max_time);
    pf.set_find_min_trace_abs_step(c.find_min_trace_abs_step);
    pf.set_find_min_locate_abs_step(c.find_min_locate_abs_step);
    pf.set_n_test_points(c.n_test_points);

    if (!c.lower_bounds.empty()) { pf.set_lower_bounds(c.lower_bounds); }
    if (!c.upper_bounds.empty()) { pf.set_upper_bounds(c.upper_bounds); }

    pf.set_t_low(c.t_low);
    pf.set_t_high(c.t_high);
    pf.set_dt_start_rel(c.dt_start_rel);
    pf.set_dt_min_rel_split_phase(c.dt_min_rel_split_phase);
    pf.set_t_jump_rel(c.t_jump_rel);
    pf.set_dt_max_abs(c.dt_max_abs);
    pf.set_dt_max_rel(c.dt_max_rel);
    pf.set_dt_min_rel(c.dt_min_rel);
    pf.set_dt_min_abs(c.dt_min_abs);

    pf.set_n_ew_scalars(c.n_ew_scalars);
    pf.set_v(c.v);

    pf.set_check_vacuum_at_low(c.check_vacuum_at_low);
    pf.set_check_vacuum_at_high(c.check_vacuum_at_high);
    pf.set_check_dx_min_dt(c.check_dx_min_dt);
    pf.set_hessian_singular_rel_tol(c.hessian_singular_rel_tol);
    pf.set_linear_algebra_rel_tol(c.linear_algebra_rel_tol);
    pf.set_check_hessian_singular(c.check_hessian_singular);
    pf.set_hessian_eig_max_rel_change(c.hessian_eig_max_rel_change);
    pf.set_check_midpoint_hessian(c.check_midpoint_hessian);

    pf.set_seed(c.seed);
    pf.set_trace_max_iter(c.trace_max_iter);
    pf.set_guess_points(c.guess_points);

    pf.set_check_merge_phase_gaps(c.check_merge_phase_gaps);
    pf.set_dt_merge_phases(c.dt_merge_phases);
    pf.set_dx_merge_phases(c.dx_merge_phases);
}

void apply(const TransitionFinderConfig& c, TransitionFinder& tf)
{
    tf.set_n_ew_scalars(c.n_ew_scalars);
    tf.set_separation(c.separation);
    tf.set_assume_only_one_transition(c.assume_only_one_transition);
    tf.set_TC_tol_rel(c.TC_tol_rel);
    tf.set_max_iter(c.max_iter);
    tf.set_change_rel_tol(c.change_rel_tol);
    tf.set_change_abs_tol(c.change_abs_tol);
    tf.set_Tnuc_step(c.Tnuc_step);
    tf.set_Tnuc_tol_rel(c.Tnuc_tol_rel);
    tf.set_check_subcritical_transitions(c.check_subcritical_transitions);
}

void apply(const ActionCalculatorConfig& c, ActionCalculator& ac)
{
    ac.set_action_calculator(c.method);
    ac.set_num_dims(c.num_dims);

    ac.set_BP_use_perturbative(c.BP_use_perturbative);
    ac.set_BP_initial_step_size(c.BP_initial_step_size);
    ac.set_BP_interpolation_points_fraction(c.BP_interpolation_points_fraction);

    ac.set_PD_xtol(c.PD_xtol);
    ac.set_PD_phitol(c.PD_phitol);
    ac.set_PD_thin_cutoff(c.PD_thin_cutoff);
    ac.set_PD_npoints(c.PD_npoints);
    ac.set_PD_rmin(c.PD_rmin);
    ac.set_PD_rmax(c.PD_rmax);
    ac.set_PD_max_iter(c.PD_max_iter);

    ac.set_PD_nb(c.PD_nb);
    ac.set_PD_kb(c.PD_kb);
    ac.set_PD_save_all_steps(c.PD_save_all_steps);
    ac.set_PD_v2min(c.PD_v2min);
    ac.set_PD_step_maxiter(c.PD_step_maxiter);
    ac.set_PD_path_maxiter(c.PD_path_maxiter);
    ac.set_PD_V_spline_samples(c.PD_V_spline_samples);
    ac.set_PD_extend_to_minima(c.PD_extend_to_minima);
    ac.set_PD_deformation_npoints(c.PD_deformation_npoints);
    ac.set_PD_fRatioConv(c.PD_fRatioConv);
    ac.set_PD_warm_start_fRatioConv(c.PD_warm_start_fRatioConv);
}

void apply(const ThermoFinderConfig& c, ThermoFinder& tm)
{
    tm.set_onset_print_setting(c.onset_print_setting);
    tm.set_percolation_print_setting(c.percolation_print_setting);
    tm.set_nucleation_print_setting(c.nucleation_print_setting);
    tm.set_completion_print_setting(c.completion_print_setting);

    tm.set_update_percolation_temperature(c.update_percolation_temperature);
    tm.set_compute_profiles(c.compute_profiles);
    tm.set_n_temp_profiles(c.n_temp_profiles);
    tm.set_vw(c.vw);

    tm.set_temperature_threshold(c.temperature_threshold);
    tm.set_vev_threshold(c.vev_threshold);
    tm.set_default_validation_method(c.default_validation_method);
    tm.set_transition_filter(c.transition_filter);

    tm.set_action_spline_evaluations(c.action_spline_evaluations);
    tm.set_warm_start_chunk_size(c.warm_start_chunk_size);
    tm.set_action_smoothing_window(c.action_smoothing_window);
    tm.set_action_smoothing_order(c.action_smoothing_order);
    tm.set_action_laurent_tail(c.action_laurent_tail);
    if (c.prefactor_function) { tm.set_prefactor_function(c.prefactor_function); }

    tm.set_eos_spline_evaluations(c.eos_spline_evaluations);
    tm.set_eos_background_dof(c.eos_background_dof);

    tm.set_percolation_target(c.percolation_target);
    tm.set_completion_target(c.completion_target);
    tm.set_onset_target(c.onset_target);
    tm.set_nucleation_target(c.nucleation_target);
    tm.set_use_bag_dtdT(c.use_bag_dtdT);
    tm.set_temperature_abs_tol(c.temperature_abs_tol);
}

void apply(const GravWaveConfig& c, GravWaveCalculator& gw)
{
    gw.set_gw_method(c.gw_method);
    gw.set_use_legacy_gw_methods(c.use_legacy_gw_methods);

    gw.set_min_kRs_value(c.min_kRs_value);
    gw.set_max_kRs_value(c.max_kRs_value);
    gw.set_n_kRs_value(c.n_kRs_value);
    gw.set_min_frequency(c.min_frequency);
    gw.set_max_frequency(c.max_frequency);
    gw.set_num_frequency(c.num_frequency);
    gw.set_num_frequency_ssm(c.num_frequency_ssm);
    gw.set_T_threshold_bubble_collision(c.T_threshold_bubble_collision);

    gw.set_h_dVdT(c.h_dVdT);
    gw.set_h_dSdT(c.h_dSdT);
    gw.set_np_dSdT(c.np_dSdT);

    gw.set_default_milestone(c.default_milestone);
    gw.set_include_col_and_turb_in_ssm(c.include_col_and_turb_in_ssm);

    gw.set_g_0(c.g_0);
    gw.set_h_0(c.h_0);
    gw.set_g_eff(c.g_eff);
    gw.set_h_eff(c.h_eff);
    gw.set_omega_hsq_neutrino(c.omega_hsq_neutrino);
    gw.set_D(c.D);
    gw.set_vw(c.vw);
    gw.set_epsilon(c.epsilon);

    gw.set_run_time_LISA(c.run_time_LISA);
    gw.set_run_time_Taiji(c.run_time_Taiji);
    gw.set_use_legacy_LISA_noise(c.use_legacy_LISA_noise);
    gw.set_SNR_f_min(c.SNR_f_min);
    gw.set_SNR_f_max(c.SNR_f_max);
    gw.set_SNR_steps_per_decade(c.SNR_steps_per_decade);
}

// =========================== validate ===========================

RunStatus Config::validate() const
{
    std::vector<std::string> problems;
    auto require = [&problems](bool condition, const std::string& what) {
        if (!condition) { problems.push_back(what); }
    };

    const auto& pf = phase_finder;
    require(pf.t_low >= 0., "phase_finder.t_low must be >= 0");
    require(pf.t_low < pf.t_high, "phase_finder.t_low must be < phase_finder.t_high");
    require(pf.lower_bounds.empty() || pf.upper_bounds.empty() || pf.lower_bounds.size() == pf.upper_bounds.size(),
            "phase_finder.lower_bounds and upper_bounds must have the same size");
    if (pf.lower_bounds.size() == pf.upper_bounds.size())
    {
        for (size_t i = 0; i < pf.lower_bounds.size(); ++i)
        {
            require(pf.lower_bounds[i] < pf.upper_bounds[i],
                    "phase_finder.lower_bounds[" + std::to_string(i) + "] must be < upper_bounds[" + std::to_string(i) + "]");
        }
    }
    require(pf.n_test_points > 0, "phase_finder.n_test_points must be > 0");

    require(transition_finder.TC_tol_rel > 0., "transition_finder.TC_tol_rel must be > 0");

    require(action.num_dims > 0, "action.num_dims must be > 0");
    require(action.PD_xtol > 0., "action.PD_xtol must be > 0");
    require(action.PD_phitol > 0., "action.PD_phitol must be > 0");
    require(action.PD_rmin < action.PD_rmax, "action.PD_rmin must be < action.PD_rmax");

    const auto& tm = thermo_finder;
    require(tm.vw > 0. && tm.vw <= 1., "thermo_finder.vw must be in (0, 1]");
    require(tm.action_spline_evaluations >= 2, "thermo_finder.action_spline_evaluations must be >= 2");
    require(tm.warm_start_chunk_size >= 0, "thermo_finder.warm_start_chunk_size must be >= 0");
    require(tm.action_smoothing_window <= 1
            || (tm.action_smoothing_window % 2 == 1 && tm.action_smoothing_order >= 0
                && tm.action_smoothing_order < tm.action_smoothing_window),
            "thermo_finder.action_smoothing_window must be odd and larger than action_smoothing_order");
    require(tm.eos_spline_evaluations >= 2, "thermo_finder.eos_spline_evaluations must be >= 2");
    require(tm.eos_background_dof >= 0., "thermo_finder.eos_background_dof must be >= 0");
    require(tm.completion_target > 0. && tm.completion_target < 1., "thermo_finder.completion_target must be in (0, 1)");
    require(tm.percolation_target > 0. && tm.percolation_target < 1., "thermo_finder.percolation_target must be in (0, 1)");
    require(tm.onset_target > 0. && tm.onset_target < 1., "thermo_finder.onset_target must be in (0, 1)");
    require(tm.completion_target < tm.percolation_target && tm.percolation_target < tm.onset_target,
            "thermo_finder targets must satisfy completion_target < percolation_target < onset_target");
    require(tm.nucleation_target > 0., "thermo_finder.nucleation_target must be > 0");
    require(tm.temperature_abs_tol > 0., "thermo_finder.temperature_abs_tol must be > 0");

    const auto& gw = gravwave;
    require(gw.min_frequency > 0., "gravwave.min_frequency must be > 0");
    require(gw.min_frequency < gw.max_frequency, "gravwave.min_frequency must be < gravwave.max_frequency");
    require(gw.num_frequency >= 2, "gravwave.num_frequency must be >= 2");
    require(gw.min_kRs_value > 0. && gw.min_kRs_value < gw.max_kRs_value,
            "gravwave.min_kRs_value must be in (0, max_kRs_value)");
    require(gw.n_kRs_value >= 2, "gravwave.n_kRs_value must be >= 2");
    require(gw.SNR_f_min > 0. && gw.SNR_f_min < gw.SNR_f_max, "gravwave.SNR_f_min must be in (0, SNR_f_max)");
    require(gw.SNR_steps_per_decade > 0., "gravwave.SNR_steps_per_decade must be > 0");
    require(gw.vw >= 0. && gw.vw <= 1., "gravwave.vw must be in [0, 1] (0 = use the milestone's value)");
    require(gw.g_eff >= 0. && gw.h_eff >= 0., "gravwave.g_eff and h_eff must be >= 0 (0 = use the milestone's value)");

    require(pipeline.stop_after != Stage::None && pipeline.stop_after != Stage::Config,
            "pipeline.stop_after must be a computation stage");

    RunStatus status;
    if (!problems.empty())
    {
        status.code = StatusCode::InvalidConfig;
        status.stage = Stage::Config;
        for (size_t i = 0; i < problems.size(); ++i)
        {
            status.message += (i ? "; " : "") + problems[i];
        }
    }
    return status;
}

// =========================== printing ===========================

namespace {

const char* name(ActionMethod m)
{
    switch (m)
    {
        case ActionMethod::None:            return "None";
        case ActionMethod::BubbleProfiler:  return "BubbleProfiler";
        case ActionMethod::PathDeformation: return "PathDeformation";
        case ActionMethod::All:             return "All";
    }
    return "Unknown";
}

const char* name(PrintSettings p)
{
    switch (p)
    {
        case MINIMAL:  return "MINIMAL";
        case STANDARD: return "STANDARD";
        case VERBOSE:  return "VERBOSE";
    }
    return "Unknown";
}

const char* name(ValidateMethod v)
{
    switch (v)
    {
        case TEMP: return "TEMP";
        case VEV:  return "VEV";
        case NONE: return "NONE";
    }
    return "Unknown";
}

const char* name(MilestoneType m)
{
    switch (m)
    {
        case PERCOLATION: return "PERCOLATION";
        case NUCLEATION:  return "NUCLEATION";
        case COMPLETION:  return "COMPLETION";
        case ONSET:       return "ONSET";
    }
    return "Unknown";
}

std::string vec_str(const std::vector<double>& v)
{
    if (v.empty()) { return "default"; }
    std::ostringstream ss;
    ss << "[";
    for (size_t i = 0; i < v.size(); ++i) { ss << (i ? ", " : "") << v[i]; }
    ss << "]";
    return ss.str();
}

const char* set_str(bool is_set) { return is_set ? "custom" : "default"; }

} // namespace

#define PT_CONFIG_PRINT(field) o << "  " #field " = " << c.field << "\n"

std::ostream& operator<<(std::ostream& o, const PhaseFinderConfig& c)
{
    o << "[phase_finder]\n";
    PT_CONFIG_PRINT(x_abs_identical);
    PT_CONFIG_PRINT(x_rel_identical);
    PT_CONFIG_PRINT(x_abs_jump);
    PT_CONFIG_PRINT(x_rel_jump);
    PT_CONFIG_PRINT(find_min_x_tol_rel);
    PT_CONFIG_PRINT(find_min_x_tol_abs);
    o << "  find_min_algorithm = " << nlopt::algorithm_name(c.find_min_algorithm) << "\n";
    PT_CONFIG_PRINT(find_min_max_f_eval);
    PT_CONFIG_PRINT(find_min_min_step);
    PT_CONFIG_PRINT(find_min_max_time);
    PT_CONFIG_PRINT(find_min_trace_abs_step);
    PT_CONFIG_PRINT(find_min_locate_abs_step);
    PT_CONFIG_PRINT(n_test_points);
    o << "  lower_bounds = " << vec_str(c.lower_bounds) << "\n";
    o << "  upper_bounds = " << vec_str(c.upper_bounds) << "\n";
    PT_CONFIG_PRINT(t_low);
    PT_CONFIG_PRINT(t_high);
    PT_CONFIG_PRINT(dt_start_rel);
    PT_CONFIG_PRINT(dt_min_rel_split_phase);
    PT_CONFIG_PRINT(t_jump_rel);
    PT_CONFIG_PRINT(dt_max_abs);
    PT_CONFIG_PRINT(dt_max_rel);
    PT_CONFIG_PRINT(dt_min_rel);
    PT_CONFIG_PRINT(dt_min_abs);
    PT_CONFIG_PRINT(n_ew_scalars);
    PT_CONFIG_PRINT(v);
    PT_CONFIG_PRINT(check_vacuum_at_low);
    PT_CONFIG_PRINT(check_vacuum_at_high);
    PT_CONFIG_PRINT(check_dx_min_dt);
    PT_CONFIG_PRINT(hessian_singular_rel_tol);
    PT_CONFIG_PRINT(linear_algebra_rel_tol);
    PT_CONFIG_PRINT(check_hessian_singular);
    PT_CONFIG_PRINT(hessian_eig_max_rel_change);
    PT_CONFIG_PRINT(check_midpoint_hessian);
    PT_CONFIG_PRINT(seed);
    PT_CONFIG_PRINT(trace_max_iter);
    o << "  guess_points = " << c.guess_points.size() << " point(s)\n";
    PT_CONFIG_PRINT(check_merge_phase_gaps);
    PT_CONFIG_PRINT(dt_merge_phases);
    PT_CONFIG_PRINT(dx_merge_phases);
    return o;
}

std::ostream& operator<<(std::ostream& o, const TransitionFinderConfig& c)
{
    o << "[transition_finder]\n";
    PT_CONFIG_PRINT(n_ew_scalars);
    PT_CONFIG_PRINT(separation);
    PT_CONFIG_PRINT(assume_only_one_transition);
    PT_CONFIG_PRINT(TC_tol_rel);
    PT_CONFIG_PRINT(max_iter);
    PT_CONFIG_PRINT(change_rel_tol);
    PT_CONFIG_PRINT(change_abs_tol);
    PT_CONFIG_PRINT(Tnuc_step);
    PT_CONFIG_PRINT(Tnuc_tol_rel);
    PT_CONFIG_PRINT(check_subcritical_transitions);
    return o;
}

std::ostream& operator<<(std::ostream& o, const ActionCalculatorConfig& c)
{
    o << "[action]\n";
    o << "  method = " << name(c.method) << "\n";
    PT_CONFIG_PRINT(num_dims);
    PT_CONFIG_PRINT(BP_use_perturbative);
    PT_CONFIG_PRINT(BP_initial_step_size);
    PT_CONFIG_PRINT(BP_interpolation_points_fraction);
    PT_CONFIG_PRINT(PD_xtol);
    PT_CONFIG_PRINT(PD_phitol);
    PT_CONFIG_PRINT(PD_thin_cutoff);
    PT_CONFIG_PRINT(PD_npoints);
    PT_CONFIG_PRINT(PD_rmin);
    PT_CONFIG_PRINT(PD_rmax);
    PT_CONFIG_PRINT(PD_max_iter);
    PT_CONFIG_PRINT(PD_nb);
    PT_CONFIG_PRINT(PD_kb);
    PT_CONFIG_PRINT(PD_save_all_steps);
    PT_CONFIG_PRINT(PD_v2min);
    PT_CONFIG_PRINT(PD_step_maxiter);
    PT_CONFIG_PRINT(PD_path_maxiter);
    PT_CONFIG_PRINT(PD_V_spline_samples);
    PT_CONFIG_PRINT(PD_extend_to_minima);
    PT_CONFIG_PRINT(PD_deformation_npoints);
    PT_CONFIG_PRINT(PD_fRatioConv);
    PT_CONFIG_PRINT(PD_warm_start_fRatioConv);
    return o;
}

std::ostream& operator<<(std::ostream& o, const ThermoFinderConfig& c)
{
    o << "[thermo_finder]\n";
    o << "  onset_print_setting = " << name(c.onset_print_setting) << "\n";
    o << "  percolation_print_setting = " << name(c.percolation_print_setting) << "\n";
    o << "  nucleation_print_setting = " << name(c.nucleation_print_setting) << "\n";
    o << "  completion_print_setting = " << name(c.completion_print_setting) << "\n";
    PT_CONFIG_PRINT(update_percolation_temperature);
    PT_CONFIG_PRINT(compute_profiles);
    PT_CONFIG_PRINT(n_temp_profiles);
    PT_CONFIG_PRINT(vw);
    PT_CONFIG_PRINT(temperature_threshold);
    PT_CONFIG_PRINT(vev_threshold);
    o << "  default_validation_method = " << name(c.default_validation_method) << "\n";
    o << "  transition_filter = " << set_str(static_cast<bool>(c.transition_filter)) << "\n";
    PT_CONFIG_PRINT(action_spline_evaluations);
    PT_CONFIG_PRINT(warm_start_chunk_size);
    PT_CONFIG_PRINT(action_smoothing_window);
    PT_CONFIG_PRINT(action_smoothing_order);
    PT_CONFIG_PRINT(action_laurent_tail);
    o << "  prefactor_function = " << set_str(static_cast<bool>(c.prefactor_function)) << "\n";
    PT_CONFIG_PRINT(eos_spline_evaluations);
    PT_CONFIG_PRINT(eos_background_dof);
    PT_CONFIG_PRINT(percolation_target);
    PT_CONFIG_PRINT(completion_target);
    PT_CONFIG_PRINT(onset_target);
    PT_CONFIG_PRINT(nucleation_target);
    PT_CONFIG_PRINT(use_bag_dtdT);
    PT_CONFIG_PRINT(temperature_abs_tol);
    return o;
}

std::ostream& operator<<(std::ostream& o, const GravWaveConfig& c)
{
    o << "[gravwave]\n";
    o << "  gw_method = " << to_string(c.gw_method) << "\n";
    PT_CONFIG_PRINT(use_legacy_gw_methods);
    PT_CONFIG_PRINT(min_kRs_value);
    PT_CONFIG_PRINT(max_kRs_value);
    PT_CONFIG_PRINT(n_kRs_value);
    PT_CONFIG_PRINT(min_frequency);
    PT_CONFIG_PRINT(max_frequency);
    PT_CONFIG_PRINT(num_frequency);
    PT_CONFIG_PRINT(num_frequency_ssm);
    PT_CONFIG_PRINT(T_threshold_bubble_collision);
    PT_CONFIG_PRINT(h_dVdT);
    PT_CONFIG_PRINT(h_dSdT);
    PT_CONFIG_PRINT(np_dSdT);
    o << "  default_milestone = " << name(c.default_milestone) << "\n";
    PT_CONFIG_PRINT(include_col_and_turb_in_ssm);
    PT_CONFIG_PRINT(g_0);
    PT_CONFIG_PRINT(h_0);
    PT_CONFIG_PRINT(g_eff);
    PT_CONFIG_PRINT(h_eff);
    PT_CONFIG_PRINT(omega_hsq_neutrino);
    PT_CONFIG_PRINT(D);
    PT_CONFIG_PRINT(vw);
    PT_CONFIG_PRINT(epsilon);
    PT_CONFIG_PRINT(run_time_LISA);
    PT_CONFIG_PRINT(run_time_Taiji);
    PT_CONFIG_PRINT(use_legacy_LISA_noise);
    PT_CONFIG_PRINT(SNR_f_min);
    PT_CONFIG_PRINT(SNR_f_max);
    PT_CONFIG_PRINT(SNR_steps_per_decade);
    return o;
}

std::ostream& operator<<(std::ostream& o, const PipelineConfig& c)
{
    o << "[pipeline]\n";
    o << "  stop_after = " << to_string(c.stop_after) << "\n";
    PT_CONFIG_PRINT(throw_on_error);
    PT_CONFIG_PRINT(to_print);
    o << "  log_level = ";
    if (c.log_level) { o << *c.log_level; } else { o << "unchanged"; }
    o << "\n";
    return o;
}

#undef PT_CONFIG_PRINT

std::ostream& operator<<(std::ostream& o, const Config& c)
{
    return o << c.phase_finder << c.transition_finder << c.action
             << c.thermo_finder << c.gravwave << c.pipeline;
}

} // namespace PhaseTracer
