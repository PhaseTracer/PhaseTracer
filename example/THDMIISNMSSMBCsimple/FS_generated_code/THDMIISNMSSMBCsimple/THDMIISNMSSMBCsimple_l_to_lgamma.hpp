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

// File generated at Sat 11 Apr 2020 12:54:08

/**
 * @file THDMIISNMSSMBCsimple_edm.hpp
 *
 * This file was generated at Sat 11 Apr 2020 12:54:08 with FlexibleSUSY
 * 2.4.2 and SARAH 4.14.3 .
 */

#ifndef THDMIISNMSSMBCsimple_MuToEGamma_H
#define THDMIISNMSSMBCsimple_MuToEGamma_H

#include "lowe.h"
#include "physical_input.hpp"

namespace flexiblesusy {
class THDMIISNMSSMBCsimple_mass_eigenstates;

namespace THDMIISNMSSMBCsimple_l_to_lgamma {
template <typename FIn, typename FOut, typename T1, typename T2>
double lepton_total_decay_width(
      T1 const&, T2 const&, 
      const THDMIISNMSSMBCsimple_mass_eigenstates&, const softsusy::QedQcd&);

}
} // namespace flexiblesusy

#endif
