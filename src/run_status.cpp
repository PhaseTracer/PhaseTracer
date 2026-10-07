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

#include "run_status.hpp"

#include <sstream>

namespace PhaseTracer {

const char* to_string(Stage s)
{
    switch (s)
    {
        case Stage::None:             return "None";
        case Stage::Config:           return "Config";
        case Stage::PhaseFinder:      return "PhaseFinder";
        case Stage::ActionCalculator: return "ActionCalculator";
        case Stage::TransitionFinder: return "TransitionFinder";
        case Stage::ThermoFinder:     return "ThermoFinder";
        case Stage::GravWave:         return "GravWave";
    }
    return "Unknown";
}

const char* to_string(StatusCode c)
{
    switch (c)
    {
        case StatusCode::Success:                return "Success";
        case StatusCode::InvalidConfig:          return "InvalidConfig";
        case StatusCode::PhaseFinderFailed:      return "PhaseFinderFailed";
        case StatusCode::NoPhases:               return "NoPhases";
        case StatusCode::ActionCalculatorFailed: return "ActionCalculatorFailed";
        case StatusCode::TransitionFinderFailed: return "TransitionFinderFailed";
        case StatusCode::NoTransitions:          return "NoTransitions";
        case StatusCode::ThermoFinderFailed:     return "ThermoFinderFailed";
        case StatusCode::NoThermalParameters:    return "NoThermalParameters";
        case StatusCode::GravWaveFailed:         return "GravWaveFailed";
        case StatusCode::NoSpectra:              return "NoSpectra";
    }
    return "Unknown";
}

std::ostream& operator<<(std::ostream& o, const RunStatus& s)
{
    o << "[" << to_string(s.code) << "]";
    if (s.stage != Stage::None) { o << " at stage " << to_string(s.stage); }
    if (!s.message.empty()) { o << ": " << s.message; }
    for (const auto& w : s.warnings) { o << "\n  warning: " << w; }
    o << "\n";
    return o;
}

namespace {
std::string describe(const RunStatus& s)
{
    std::ostringstream ss;
    ss << s;
    return ss.str();
}
} // namespace

RunnerError::RunnerError(RunStatus s) : std::runtime_error(describe(s)), status_(std::move(s)) {}

} // namespace PhaseTracer
