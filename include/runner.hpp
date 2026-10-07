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

#ifndef PHASETRACER_RUNNER_HPP_
#define PHASETRACER_RUNNER_HPP_

/**
 * Hands-off interface to the full PhaseTracer pipeline
 * (PhaseFinder -> ActionCalculator -> TransitionFinder -> ThermoFinder -> GravWaveCalculator).
 *
 *   auto config = PhaseTracer::Config();
 *   config.phase_finder.seed = 0;
 *   PhaseTracer::Runner runner(model, config);
 *   auto status = runner.run();
 *   if (!status) std::cerr << status;
 *   const auto& transitions = runner.get_transitions();
 *
 * The Runner owns the stage objects but no copies of their results: the get_* wrappers
 * return references into the stage objects, valid until the next run() or until the
 * Runner is destroyed.
 */

#include <memory>
#include <string>
#include <vector>

#include "potential.hpp"
#include "config.hpp"
#include "run_status.hpp"
#include "phase_finder.hpp"
#include "action_calculator.hpp"
#include "transition_finder.hpp"
#include "thermo_finder.hpp"
#include "gravwave_calculator.hpp"

namespace PhaseTracer {

class Runner
{
public:
    /** @brief The model is held by reference and must outlive the Runner. */
    Runner(EffectivePotential::Potential& model, Config config = {});
    Runner(EffectivePotential::Potential&&, Config = {}) = delete;

    Runner(const Runner&) = delete;
    Runner& operator=(const Runner&) = delete;

    /**
     * @brief Rebuilds every stage from scratch and runs up to config().pipeline.stop_after.
     * @return Success, or the stage and reason the run stopped.
     * @throws RunnerError instead of returning a failed status if config().pipeline.throw_on_error is set.
     */
    RunStatus run();

    /** @brief Status of the last run(). */
    const RunStatus& status() const { return status_; }

    /** @brief Settings used by the next run(). */
    Config& config() { return config_; }
    const Config& config() const { return config_; }

    // Stage objects; each throws std::logic_error if that stage was not constructed in the last run().
    PhaseFinder& phase_finder();
    ActionCalculator& action_calculator();
    TransitionFinder& transition_finder();
    ThermoFinder& thermo_finder();
    GravWaveCalculator& gravwave_calculator();
    const PhaseFinder& phase_finder() const;
    const TransitionFinder& transition_finder() const;
    const GravWaveCalculator& gravwave_calculator() const;

    /** @brief Whether the stage object was constructed in the last run(). */
    bool has(Stage s) const;

    // Results; wrappers for the stage getters, throwing std::logic_error if the stage was not constructed.
    /** @brief PhaseFinder::get_phases() */
    const std::vector<Phase>& get_phases() const { return phase_finder().get_phases(); }
    /** @brief TransitionFinder::get_transitions() */
    const std::vector<Transition>& get_transitions() const { return transition_finder().get_transitions(); }
    /** @brief ThermoFinder::get_thermal_parameters() */
    const std::vector<ThermalParameterSet>& get_thermal_parameters() { return thermo_finder().get_thermal_parameters(); }
    /** @brief GravWaveCalculator::get_spectrums() */
    const std::vector<GravWaveSpectrum>& get_spectra() const { return gravwave_calculator().get_spectrums(); }

private:
    EffectivePotential::Potential& model_;
    Config config_;
    RunStatus status_;

    // Declared in construction order, so they are destroyed in reverse: each stage holds
    // references or pointers to the ones before it.
    std::unique_ptr<PhaseFinder> pf_;
    std::unique_ptr<ActionCalculator> ac_;
    std::unique_ptr<TransitionFinder> tf_;
    std::unique_ptr<ThermoFinder> tm_;
    std::unique_ptr<GravWaveCalculator> gw_;

    /** @brief Destroys every stage object and clears the status. */
    void reset();

    /** @brief Runs f, turning any exception into a failed status_. Returns false on failure. */
    template <class F>
    bool stage(Stage s, StatusCode fail_code, F&& f);

    /** @brief Records a failed status_ and returns false. */
    bool fail(Stage s, StatusCode code, std::string message);

    /** @brief Whether the last completed stage is the one to stop after. */
    bool done(Stage s) const { return config_.pipeline.stop_after == s; }

    /** @brief Returns status_, or throws it as a RunnerError if requested. */
    RunStatus finish();
};

} // namespace PhaseTracer

#endif // PHASETRACER_RUNNER_HPP_
