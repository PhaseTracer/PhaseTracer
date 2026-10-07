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

#include "runner.hpp"

#include <iostream>
#include <optional>
#include <utility>

#include "logger.hpp"

namespace PhaseTracer {

namespace {

const TransitionMilestone& milestone_of(const ThermalParameterSet& tps, MilestoneType type)
{
    switch (type)
    {
        case ONSET:      return tps.onset;
        case NUCLEATION: return tps.nucleation;
        case COMPLETION: return tps.completion;
        case PERCOLATION:
        default:         return tps.percolation;
    }
}

const char* milestone_name(MilestoneType type)
{
    switch (type)
    {
        case ONSET:       return "onset";
        case NUCLEATION:  return "nucleation";
        case COMPLETION:  return "completion";
        case PERCOLATION: return "percolation";
    }
    return "unknown";
}

std::logic_error not_constructed(Stage s)
{
    return std::logic_error(std::string("Runner: ") + to_string(s)
                            + " was not constructed in the last run(); check status()");
}

} // namespace

Runner::Runner(EffectivePotential::Potential& model, Config config)
    : model_(model), config_(std::move(config)) {}

// =========================== stage access ===========================

PhaseFinder& Runner::phase_finder()
{
    if (!pf_) { throw not_constructed(Stage::PhaseFinder); }
    return *pf_;
}

const PhaseFinder& Runner::phase_finder() const
{
    if (!pf_) { throw not_constructed(Stage::PhaseFinder); }
    return *pf_;
}

ActionCalculator& Runner::action_calculator()
{
    if (!ac_) { throw not_constructed(Stage::ActionCalculator); }
    return *ac_;
}

TransitionFinder& Runner::transition_finder()
{
    if (!tf_) { throw not_constructed(Stage::TransitionFinder); }
    return *tf_;
}

const TransitionFinder& Runner::transition_finder() const
{
    if (!tf_) { throw not_constructed(Stage::TransitionFinder); }
    return *tf_;
}

ThermoFinder& Runner::thermo_finder()
{
    if (!tm_) { throw not_constructed(Stage::ThermoFinder); }
    return *tm_;
}

GravWaveCalculator& Runner::gravwave_calculator()
{
    if (!gw_) { throw not_constructed(Stage::GravWave); }
    return *gw_;
}

const GravWaveCalculator& Runner::gravwave_calculator() const
{
    if (!gw_) { throw not_constructed(Stage::GravWave); }
    return *gw_;
}

bool Runner::has(Stage s) const
{
    switch (s)
    {
        case Stage::PhaseFinder:      return static_cast<bool>(pf_);
        case Stage::ActionCalculator: return static_cast<bool>(ac_);
        case Stage::TransitionFinder: return static_cast<bool>(tf_);
        case Stage::ThermoFinder:     return static_cast<bool>(tm_);
        case Stage::GravWave:         return static_cast<bool>(gw_);
        default:                      return false;
    }
}

// =========================== helpers ===========================

void Runner::reset()
{
    gw_.reset();
    tm_.reset();
    tf_.reset();
    ac_.reset();
    pf_.reset();
    status_ = RunStatus{};
}

template <class F>
bool Runner::stage(Stage s, StatusCode fail_code, F&& f)
{
    try {
        f();
        return true;
    } catch (const std::exception& e) {
        return fail(s, fail_code, e.what());
    } catch (...) {
        return fail(s, fail_code, "unknown exception");
    }
}

bool Runner::fail(Stage s, StatusCode code, std::string message)
{
    status_.code = code;
    status_.stage = s;
    status_.message = std::move(message);
    LOG(debug) << "Runner stopped: " << status_;
    return false;
}

RunStatus Runner::finish()
{
    if (config_.pipeline.throw_on_error && !status_.ok())
    {
        throw RunnerError(status_);
    }
    return status_;
}

// =========================== run ===========================

RunStatus Runner::run()
{
    reset();
    const auto& pipeline = config_.pipeline;

    if (pipeline.log_level)
    {
        boost::log::core::get()->set_filter(boost::log::trivial::severity >= *pipeline.log_level);
    }

    status_ = config_.validate();
    if (!status_.ok()) { return finish(); }

    // PhaseFinder
    if (!stage(Stage::PhaseFinder, StatusCode::PhaseFinderFailed, [&] { pf_ = std::make_unique<PhaseFinder>(model_); })
        || !stage(Stage::PhaseFinder, StatusCode::InvalidConfig, [&] { apply(config_.phase_finder, *pf_); })
        || !stage(Stage::PhaseFinder, StatusCode::PhaseFinderFailed, [&] { pf_->find_phases(); }))
    {
        return finish();
    }
    if (pf_->get_phases().empty())
    {
        fail(Stage::PhaseFinder, StatusCode::NoPhases, "PhaseFinder found no phases");
        return finish();
    }
    if (pipeline.to_print) { std::cout << *pf_; }
    if (done(Stage::PhaseFinder)) { return finish(); }

    // ActionCalculator: copies the PhaseFinder, so must be built after find_phases()
    if (!stage(Stage::ActionCalculator, StatusCode::ActionCalculatorFailed, [&] { ac_ = std::make_unique<ActionCalculator>(*pf_); })
        || !stage(Stage::ActionCalculator, StatusCode::InvalidConfig, [&] { apply(config_.action, *ac_); }))
    {
        return finish();
    }
    if (done(Stage::ActionCalculator)) { return finish(); }

    // TransitionFinder: never given the ActionCalculator, since nucleation comes from ThermoFinder
    if (!stage(Stage::TransitionFinder, StatusCode::TransitionFinderFailed, [&] { tf_ = std::make_unique<TransitionFinder>(*pf_); })
        || !stage(Stage::TransitionFinder, StatusCode::InvalidConfig, [&] { apply(config_.transition_finder, *tf_); })
        || !stage(Stage::TransitionFinder, StatusCode::TransitionFinderFailed, [&] { tf_->find_transitions(); }))
    {
        return finish();
    }
    if (tf_->get_transitions().empty())
    {
        fail(Stage::TransitionFinder, StatusCode::NoTransitions, "TransitionFinder found no transitions");
        return finish();
    }
    if (pipeline.to_print) { std::cout << *tf_; }
    if (done(Stage::TransitionFinder)) { return finish(); }

    // ThermoFinder: holds its own copy of the (already evaluated) TransitionFinder
    if (!stage(Stage::ThermoFinder, StatusCode::ThermoFinderFailed,
               [&] { tm_ = std::make_unique<ThermoFinder>(std::optional<TransitionFinder>(*tf_), *ac_); })
        || !stage(Stage::ThermoFinder, StatusCode::InvalidConfig, [&] { apply(config_.thermo_finder, *tm_); })
        || !stage(Stage::ThermoFinder, StatusCode::ThermoFinderFailed, [&] { tm_->find_thermal_parameters(); }))
    {
        return finish();
    }
    status_.warnings = tm_->get_failure_messages();
    if (tm_->get_thermal_parameters().empty())
    {
        const auto n_failed = status_.warnings.size();
        fail(Stage::ThermoFinder, StatusCode::NoThermalParameters,
             n_failed == 0 ? "no transition passed the transition filter"
                           : "thermal parameters failed for all " + std::to_string(n_failed)
                                 + " transition(s) that passed the filter (see warnings)");
        return finish();
    }
    if (pipeline.to_print) { std::cout << *tm_; }
    if (done(Stage::ThermoFinder)) { return finish(); }

    // GravWave: calc_spectrums() throws on an empty spectrum list, so check a milestone was reached first
    const auto milestone = config_.gravwave.default_milestone;
    bool any_milestone = false;
    for (const auto& tps : tm_->get_thermal_parameters())
    {
        any_milestone = any_milestone || milestone_of(tps, milestone).status == MilestoneStatus::YES;
    }
    if (!any_milestone)
    {
        fail(Stage::GravWave, StatusCode::NoSpectra,
             std::string("no transition reached the ") + milestone_name(milestone) + " milestone");
        return finish();
    }
    if (!stage(Stage::GravWave, StatusCode::GravWaveFailed, [&] { gw_ = std::make_unique<GravWaveCalculator>(*tm_); })
        || !stage(Stage::GravWave, StatusCode::InvalidConfig, [&] { apply(config_.gravwave, *gw_); })
        || !stage(Stage::GravWave, StatusCode::GravWaveFailed, [&] { gw_->calc_spectrums(); }))
    {
        return finish();
    }
    if (gw_->get_spectrums().empty())
    {
        fail(Stage::GravWave, StatusCode::NoSpectra, "GravWaveCalculator produced no spectra");
        return finish();
    }
    if (pipeline.to_print) { std::cout << *gw_; }

    return finish();
}

} // namespace PhaseTracer
