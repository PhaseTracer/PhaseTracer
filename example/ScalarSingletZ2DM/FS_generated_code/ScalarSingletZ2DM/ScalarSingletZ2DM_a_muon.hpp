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

// File generated at Sat 24 Oct 2020 17:07:56

/**
 * @file ScalarSingletZ2DM_a_muon.hpp
 *
 * This file was generated at Sat 24 Oct 2020 17:07:56 with FlexibleSUSY
 * 2.4.2 and SARAH 4.14.3 .
 */

#ifndef ScalarSingletZ2DM_A_MUON_H
#define ScalarSingletZ2DM_A_MUON_H

namespace softsusy {
   class QedQcd;
} // namespace softsusy

namespace flexiblesusy {
class ScalarSingletZ2DM_mass_eigenstates;

namespace ScalarSingletZ2DM_a_muon {
/**
* @fn calculate_a_muon
* @brief Calculates \f$a_\mu = (g-2)_\mu/2\f$ of the muon.
*/
double calculate_a_muon(const ScalarSingletZ2DM_mass_eigenstates& model, const softsusy::QedQcd& qedqcd);

/**
* @fn calculate_a_muon_uncertainty
* @brief Calculates \f$\Delta a_\mu\f$ of the muon.
*/
double calculate_a_muon_uncertainty(const ScalarSingletZ2DM_mass_eigenstates& model, const softsusy::QedQcd& qedqcd);
} // namespace ScalarSingletZ2DM_a_muon
} // namespace flexiblesusy

#endif
