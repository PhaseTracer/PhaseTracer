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
//    const double v = 245.5724484517838;
//    const double mh = 125;
////
//    const double g = 0.6508563183329309;
//    const double gp = 0.3576314692605589;
//    const double yt = 0.9962622513632488;
//    const double yb = 0.01644374858393946;
//    const double ytau = 0.01023322576907503;

    const double v = 247.4554441277994;
    const double mh = 125;
//
    const double g = 0.6477108911331643;
    const double gp = 0.3585642676423744;
    const double yt = 0.9341420067357675;
    const double yb = 0.01547368604435034;
    const double ytau = 0.01001414937355005;



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
