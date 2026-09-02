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

#ifndef PHASETRACER_BUBBLE_STATISTICS_HPP_
#define PHASETRACER_BUBBLE_STATISTICS_HPP_

#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <interpolation.h>

#include "scale.hpp"

namespace PhaseTracer {

/**
 * @struct RadiiDistribution
 * @namespace PhaseTracer
 * @brief Represents the distribution of bubble radii at different temperatures.
 * 
 * @note This class is currently under development.
 */
struct RadiiDistribution 
{
    double ref_temperature;

    std::vector<double> temperature_values;
    std::vector<double> radii_values;
    std::vector<double> dndR_values;
    std::vector<double> log_dndR_values;

    double peak_radius;
    double peak_nuc_temperature;

    RadiiDistribution(
        const double& ref_temperature_in, 
        const std::vector<double>& temperature_values_in, 
        const std::vector<double>& radii_values_in, 
        const std::vector<double>& dndR_values_in,
        const std::vector<double>& log_dndR_values_in) :
    ref_temperature(ref_temperature_in), temperature_values(temperature_values_in), radii_values(radii_values_in), dndR_values(dndR_values_in), log_dndR_values(log_dndR_values_in)
    {
        alglib::real_1d_array t_array, r_array, log_dndR_array;

        t_array.setcontent(temperature_values.size(), temperature_values.data());
        r_array.setcontent(radii_values.size(), radii_values.data());
        log_dndR_array.setcontent(log_dndR_values.size(), log_dndR_values.data());

        alglib::spline1dbuildcubic(r_array, t_array, temperature_spline);
        alglib::spline1dbuildcubic(r_array, log_dndR_array, log_dndR_spline);
    }

    const double 
    get_nucleation_temperature(const double& radius)
    {
        double temperature = alglib::spline1dcalc(temperature_spline, radius);
        return temperature;
    }

    const double 
    get_dndR(const double& radius)
    {
        double log_dndR = alglib::spline1dcalc(log_dndR_spline, radius);
        return exp(log_dndR);
    }

private:

    /** Spline to extract nucleation temp from given radius */
    alglib::spline1dinterpolant temperature_spline;

    /** Spline to dndR from given radius */
    alglib::spline1dinterpolant log_dndR_spline;

    /** Extracts peak radius and temperature of dndR curve */
    // double find_peak_radius();

}; // struct RadiiDistribution

/**
 * @struct LifetimeDistribution
 * @namespace PhaseTracer
 * @brief Represents the lifetime distribution of bubbles in a phase transition.
 * 
 * @note This class is currently under development.
 */
struct LifetimeDistribution
{
    /** @brief Characteristic timescale of the lifetime distribution. */
    double timescale;

    /** @brief Mean lifetime of bubbles. */
    double mean_lifetime;

    /** @brief Dimensionless parameter chi = beta * t. */
    std::vector<double> chi_values;

    /** @brief Lifetime values corresponding to chi_values. */
    std::vector<double> lifetime_values;

    /** @brief Probability distribution values corresponding to lifetime_values. */
    std::vector<double> distribution_values;

    /** @brief Spline for log(I_3) interpolation. */
    alglib::spline1dinterpolant log_I3_spline;

    /** @brief Spline for log(I_2) interpolation. */
    alglib::spline1dinterpolant log_I2_spline;

    /** @brief Spline for T_false interpolation. */
    alglib::spline1dinterpolant T_false_spline;

    /** @brief Spline for scale factor interpolation. */
    alglib::spline1dinterpolant scale_factor_spline;

}; // struct LifetimeDistribution

} // namespace PhaseTracer

#endif // PHASETRACER_BUBBLE_STATISTICS_HPP_