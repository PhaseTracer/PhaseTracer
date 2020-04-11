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

// File generated at Sat 11 Apr 2020 12:54:09

/**
 * @file THDMIISNMSSMBCsimple_a_muon.hpp
 *
 * This file was generated at Sat 11 Apr 2020 12:54:09 with FlexibleSUSY
 * 2.4.2 and SARAH 4.14.3 .
 */

#ifndef THDMIISNMSSMBCsimple_A_MUON_H
#define THDMIISNMSSMBCsimple_A_MUON_H

namespace softsusy {
   class QedQcd;
} // namespace softsusy

namespace flexiblesusy {
class THDMIISNMSSMBCsimple_mass_eigenstates;

namespace THDMIISNMSSMBCsimple_a_muon {
/**
* @fn calculate_a_muon
* @brief Calculates \f$a_\mu = (g-2)_\mu/2\f$ of the muon.
*/
double calculate_a_muon(const THDMIISNMSSMBCsimple_mass_eigenstates& model, const softsusy::QedQcd& qedqcd);

/**
* @fn calculate_a_muon_uncertainty
* @brief Calculates \f$\Delta a_\mu\f$ of the muon.
*/
double calculate_a_muon_uncertainty(const THDMIISNMSSMBCsimple_mass_eigenstates& model, const softsusy::QedQcd& qedqcd);
} // namespace THDMIISNMSSMBCsimple_a_muon
} // namespace flexiblesusy

#endif
