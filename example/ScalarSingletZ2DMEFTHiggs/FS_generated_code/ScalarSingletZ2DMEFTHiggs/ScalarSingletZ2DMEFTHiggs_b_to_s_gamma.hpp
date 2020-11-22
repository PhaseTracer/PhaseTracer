// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================


/**
 * @file ScalarSingletZ2DMEFTHiggs_b_to_s_gamma.hpp
 *
 * This file was generated with FlexibleSUSY 2.5.0 and SARAH 4.14.3 .
 */

#ifndef ScalarSingletZ2DMEFTHiggs_BToSGamma_H
#define ScalarSingletZ2DMEFTHiggs_BToSGamma_H

namespace softsusy {
class QedQcd;
}

namespace flexiblesusy {
class ScalarSingletZ2DMEFTHiggs_mass_eigenstates;

namespace ScalarSingletZ2DMEFTHiggs_b_to_s_gamma {
    std::array<std::complex<double>, 4> calculate_b_to_s_gamma(
            const ScalarSingletZ2DMEFTHiggs_mass_eigenstates&, const softsusy::QedQcd&);
}
} // namespace flexiblesusy

#endif
