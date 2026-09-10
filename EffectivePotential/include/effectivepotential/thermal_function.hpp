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

#ifndef EFFECTIVEPOTENTIAL_THERMAL_FUNCTION_HPP_
#define EFFECTIVEPOTENTIAL_THERMAL_FUNCTION_HPP_

namespace EffectivePotential {

/**
 * @brief functions J_B(x) and J_F(x), with x = m^2 / T^2.
 * @param x The argument of the thermal function, x = m^2 / T^2.
 * @return The value of the thermal function at the given x.
 */
double J_B(double x);
double J_F(double x);

/**
 * @brief First derivatives dJ_B/dx and dJ_F/dx.
 * @param x The argument of the thermal function, x = m^2 / T^2.
 * @return The value of the derivative of the thermal function at the given x.
 */
double J_B_diff(double x);
double J_F_diff(double x);

} // namespace EffectivePotential

#endif // EFFECTIVEPOTENTIAL_THERMAL_FUNCTION_HPP_