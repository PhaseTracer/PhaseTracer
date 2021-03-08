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
    const double v = 247.4539585898549;
    const double mh = 125.2;
//
    const double g = 0.6477079010403815;
    const double gp = 0.3585657670932486;
    const double yt = 0.9346604329164533;
    const double yb = 0.01547231441771844;
    const double ytau = 0.01001420531547397;

    const double mZ = 0.5*std::sqrt(square(g)+square(gp))*v;
    const double mW = 0.5*g*v;

    const double yt_sq = square(yt);
    const double yb_sq = square(yt);
    const double ytau_sq = square(ytau);

    const double mtop = 173.2;




} // namespace SM

#endif
