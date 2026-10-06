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

#ifndef PHASETRACER_RUN_STATUS_HPP_
#define PHASETRACER_RUN_STATUS_HPP_

#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace PhaseTracer {

/** @brief Stages of the full PhaseTracer pipeline, in the order they are run. */
enum class Stage
{
    None,
    Config,
    PhaseFinder,
    ActionCalculator,
    TransitionFinder,
    ThermoFinder,
    GravWave
};

/**
 * @brief Outcome of a pipeline run.
 */
enum class StatusCode
{
    Success,
    InvalidConfig,
    PhaseFinderFailed,
    NoPhases,
    TransitionFinderFailed,
    NoTransitions,
    ThermoFinderFailed,
    NoThermalParameters,
    GravWaveFailed,
    NoSpectra
};

const char* to_string(Stage s);
const char* to_string(StatusCode c);

/** @brief Status returned by a pipeline run or by Config::validate(). */
struct RunStatus
{
    StatusCode code = StatusCode::Success;
    Stage stage = Stage::None;
    std::string message;
    std::vector<std::string> warnings;

    bool ok() const { return code == StatusCode::Success; }
    explicit operator bool() const { return ok(); }
};

std::ostream& operator<<(std::ostream& o, const RunStatus& s);

/** @brief Thrown instead of returning a failed RunStatus when throw_on_error is set. */
class RunnerError : public std::runtime_error
{
public:
    explicit RunnerError(RunStatus s);
    const RunStatus& status() const { return status_; }

private:
    RunStatus status_;
};

} // namespace PhaseTracer

#endif // PHASETRACER_RUN_STATUS_HPP_
