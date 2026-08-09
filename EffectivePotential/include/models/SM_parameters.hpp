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

namespace SM {
// SM parameters
//    const double v = 245.5724484517838;
//    const double mh = 125;
////
//    const double g = 0.6508563183329309;
//    const double gp = 0.3576314692605589;
//    const double yt = 0.9962622513632488;
//    const double yb = 0.01644374858393946;
//    const double ytau = 0.01023322576907503;

const double v = 247.4544243292407;
const double mh = 125.25;
//
const double g = 0.6477096097526751;
const double gp = 0.3585644903737741;
const double yt = 0.9341361105006658;
const double yb = 0.0154737491181327;
const double ytau = 0.01001419077147658;

const double mZ = 0.5 * std::sqrt(square(g) + square(gp)) * v;
const double mW = 0.5 * g * v;

const double yt_sq = square(yt);
const double yb_sq = square(yb);
const double ytau_sq = square(ytau);

const double mtop = yt * v / sqrt(2);
const double mb = yb * v / sqrt(2);
const double mtau = ytau * v / sqrt(2);

// Light quark MS-bar masses (PDG 2022) and derived Yukawa couplings
// m_c, m_s, m_u, m_d are approximate values; corrections to the potential
// from these small Yukawas are negligible at the EW scale.
const double mc    = 1.27;      // GeV (charm)
const double ms_q  = 0.096;     // GeV (strange)
const double mu_q  = 0.0022;    // GeV (up)
const double md_q  = 0.0047;    // GeV (down)
const double mmu   = 0.10566;   // GeV (muon)
const double mel   = 0.000511;  // GeV (electron)

const double yc   = std::sqrt(2.) * mc    / v;
const double ys   = std::sqrt(2.) * ms_q  / v;
const double yu   = std::sqrt(2.) * mu_q  / v;
const double yd   = std::sqrt(2.) * md_q  / v;
const double ymu  = std::sqrt(2.) * mmu   / v;
const double ye   = std::sqrt(2.) * mel   / v;

const double yc_sq  = square(yc);
const double ys_sq  = square(ys);
const double yu_sq  = square(yu);
const double yd_sq  = square(yd);
const double ymu_sq = square(ymu);
const double ye_sq  = square(ye);

// Strong coupling at M_Z (alpha_s(M_Z) = 0.1181, PDG 2022)
const double alpha_s = 0.1181;
const double gs      = std::sqrt(4. * M_PI * alpha_s);
const double gs_sq   = square(gs);

} // namespace SM

#endif
