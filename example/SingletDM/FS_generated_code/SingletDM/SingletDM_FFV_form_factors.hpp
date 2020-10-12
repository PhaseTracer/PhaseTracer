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
 * @file SingletDM_FFV_form_factors.hpp
 *
 * This file was generated at Tue 18 Aug 2020 22:47:57 with FlexibleSUSY
 * 2.4.2 and SARAH 4.14.3 .
 */

#ifndef SingletDM_FFVFormFactors_H
#define SingletDM_FFVFormFactors_H

#include <valarray>

#include "cxx_qft/SingletDM_qft.hpp"

namespace flexiblesusy {

class SingletDM_mass_eigenstates;

namespace SingletDM_FFV_form_factors {
std::valarray<std::complex<double>> calculate_Fe_Fe_VP_form_factors (
   int generationIndex1, int generationIndex2,
   const SingletDM_mass_eigenstates& model, bool discard_SM_contributions);
}

} // namespace flexiblesusy

#endif

