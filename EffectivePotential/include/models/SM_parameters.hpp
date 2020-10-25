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

#ifndef POTENTIAL_SM_PARAMETERS_HPP_INCLUDED
#define POTENTIAL_SM_PARAMETERS_HPP_INCLUDED

#include <cmath>
#include "pow.hpp"


namespace SM
{
    // SM parameters
    const double v = 246.221;
    const double mh = 125.2;
    const double mtop = 173.2;
    const double mZ = 91.1876;
    const double mW = 80.385;
    const double g = 2. * mW / v;
    const double gp = std::sqrt(square(2.* mZ / v) - square(g));
    const double yt = std::sqrt(2.) * mtop / v;
    const double yt_sq = square(yt);

} // namespace SM

#endif
