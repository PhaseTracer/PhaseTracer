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

#ifndef PHASETRACER_TRANSITION_MILESTONES_HPP_
#define PHASETRACER_TRANSITION_MILESTONES_HPP_

#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>

#include "scale.hpp"

namespace PhaseTracer {

/**
 * @enum MilestoneStatus
 * @brief Enumeration of possible statuses for a transition milestone.
 * 
 * This enumeration defines the possible statuses that a transition milestone can have, indicating whether the milestone 
 * has or has not been reached. There are also status for a fast transition, and an error flag.
 */
enum MilestoneStatus
{
    YES,
    FAST,
    NO,
    ERR
}; // enum MilestoneStatus

/**
 * @enum MilestoneType
 * @brief Enumeration of possible types of transition milestones.
 */
enum MilestoneType
{
    PERCOLATION,
    NUCLEATION,
    COMPLETION,
    ONSET
}; // enum MilestoneType

/**
 * @enum PrintSettings
 * @brief Enumeration of print settings for transition milestones.
 */
enum PrintSettings
{
    MINIMAL,
    STANDARD,
    VERBOSE
}; // enum PrintSettings

/**
 * @enum NucleationType
 * @brief Enumeration of possible nucleation types.
 */
enum NucleationType
{
    EXPONENTIAL,
    SIMULTANEOUS
}; // enum NucleationType

/**
 * @struct NucleationHistory
 * @brief Structure representing the nucleation history.
 */
struct NucleationHistory
{
    /** @brief Type of nucleation. */
    NucleationType nucleation_type;

    /** @brief First beta parameter. */
    double betaH_1;

    /** @brief Second beta parameter. */
    double betaH_2;

    /** @brief Nucleation temperature. */
    double T_m;
}; // struct NucleationHistory

/**
 * @struct TransitionMilestone
 * @brief Structure representing a transition milestone.
 * 
 * This stucture contains information about a specific transition milestone. This includes its type, status, 
 * temperature, and other relevant thermal parameters. It also provides functionality for pretty printing the milestone 
 * information.
 */
struct TransitionMilestone
{
    /** @brief Type of the milestone. */
    MilestoneType type;

    /** @brief Status of the milestone. */
    MilestoneStatus status;

    /** @brief Type of nucleation. */
    NucleationType nucleation_type = NucleationType::EXPONENTIAL;

    /** @brief Temperature of the milestone. */
    double temperature;

    /** @brief Reheating temperature, if applicable. */
    double reheating_temperature;

    /** @brief Transition strength. */
    double alpha;

    /** @brief Transition strength (mu nu prescription). */
    double alpha_munu;

    /** @brief Beta/H in the usual approximation. */
    double betaH;

    /** @brief Full Beta/H from the first time derivative */
    double beta1H;

    /** @brief Full Beta_2/H from the second time derivative */
    double beta2H;

    /** @brief Effective timescale defined from mean bubble separation. */
    double betaH_eff;

    /** @brief Hubble rate */
    double H;

    /** @brief Enthalpy to energy density ratio. */
    double we;

    /** @brief Sound speed in the plus phase. */
    double cs_plus;

    /** @brief Sound speed in the minus phase. */
    double cs_minus;

    /** @brief Average bubble number density. */
    double n;

    /** @brief Average bubble separation. */
    double Rs;

    /** @brief Average bubble radius. */
    double Rbar;

    /** @brief Duration between T_C and the milestone temperature. */
    double dt;

private:

    /** @brief Print settings for milestone. */
    PROPERTY(PrintSettings, print_setting, PrintSettings::STANDARD);

    /** @brief Helper to format doubles. 
     * 
     * This formats a double into scientific notation with a specified precision. It is used for pretty printing the 
     * transition milestone information,
     * 
     * @param value The double value to format.
     * @param precision The number of decimal places to include in the formatted string (default is 6).
     * @return A string representation of the double value with the specified precision.
     *
    */
    const std::string format_double(double value, int precision = 6) const {
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(precision) << value;
        return oss.str();
    }

    /** @brief Helper to format the milestone status as a string. */
    const std::string format_status_string() const {
        switch (status)
        {
            case MilestoneStatus::YES:
                return "YES";
            case MilestoneStatus::FAST:
                return "FAST";
            case MilestoneStatus::NO:
                return "NO";
            default:
                return "ERR";
        }
    }

    /** @brief Helper to format the milestone as a string. */
    const std::string format_milestone_string() const
    {
        std::string output = "  status = " + format_status_string() + "\n";
        output += "  temperature = " + std::to_string(temperature) + " " + scale.name() + "\n";
        if(type == MilestoneType::COMPLETION && status == MilestoneStatus::YES)
        {
            output += "  reheating temperature = " + std::to_string(reheating_temperature) + " " + scale.name() + "\n";
        }
        if(type == MilestoneType::PERCOLATION && status == MilestoneStatus::YES)
        {
            output += "  nucleation_type = " + std::string(nucleation_type == NucleationType::EXPONENTIAL ? "exponential" : "simultaneous") + "\n";
        }

        if (print_setting == PrintSettings::MINIMAL) {
            return output;
        }

        if (print_setting == PrintSettings::STANDARD || print_setting == PrintSettings::VERBOSE) {
            output += "  alpha = " + std::to_string(alpha) + "\n";
            output += "  alpha_munu = " + std::to_string(alpha_munu) + "\n";
            output += "  betaH = " + std::to_string(betaH) + "\n";
            output += "  betaH_eff = " + std::to_string(betaH_eff) + "\n";
            output += "  H = " + format_double(H) + "\n";
        }

        if (print_setting == PrintSettings::VERBOSE) {
            output += "  beta1H = " + std::to_string(beta1H) + "\n"; 
            output += "  beta2H = " + std::to_string(beta2H) + "\n";
            output += "  we = " + std::to_string(we) + "\n";
            output += "  cs_plus = " + std::to_string(cs_plus) + "\n";
            output += "  cs_minus = " + std::to_string(cs_minus) + "\n";
            output += "  Rs = " + std::to_string(Rs) + "\n";
            output += "  Rbar = " + std::to_string(Rbar) + "\n";
            output += "  dt = " + std::to_string(dt) + "\n";
        }
        
        return output;
    }

public:
    TransitionMilestone() = default;
    TransitionMilestone(const MilestoneType& type_in)
    : type(type_in), status(MilestoneStatus::ERR), temperature(0.0) {}

    /** @brief Pretty print the milestone. */
    friend std::ostream &operator<<(std::ostream& o, const TransitionMilestone &milestone) 
    {
        o << milestone.format_milestone_string();
        return o;
    }

}; // struct TransitionMilestone

} // namespace PhaseTracer

#endif // PHASETRACER_TRANSITION_MILESTONES_HPP_