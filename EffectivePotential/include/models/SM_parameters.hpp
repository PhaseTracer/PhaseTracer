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
    const double v = 247.4660195894076;
    const double mh = 125;
//
    const double g = 0.6477243040858367;
    const double gp = 0.3585619364547326;
    const double yt = 0.9342037748936957;
    const double yb = 0.0154730320228006;
    const double ytau = 0.01001372009234519;

    const double mZ = 0.5*std::sqrt(square(g)+square(gp))*v;
    const double mW = 0.5*g*v;

    const double yt_sq = square(yt);
    const double yb_sq = square(yb);
    const double ytau_sq = square(ytau);

    const double mtop = yt*v/sqrt(2);
    const double mb = yb*v/sqrt(2);
    const double mtau = ytau*v/sqrt(2);
    




} // namespace SM

#endif
