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
 * @file cxx_qft/THDMIISNMSSMBCsimple_vertices.cpp
 *
 * This file was generated at Sat 11 Apr 2020 12:54:09 with FlexibleSUSY
 * 2.4.2 and SARAH 4.14.3 .
 */

#include "THDMIISNMSSMBCsimple_context_base.hpp"
#include "THDMIISNMSSMBCsimple_vertices.hpp"

#include "concatenate.hpp"
#include "wrappers.hpp"

#define INPUTPARAMETER(p) context.model.get_input().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

namespace flexiblesusy {
namespace THDMIISNMSSMBCsimple_cxx_diagrams {
namespace detail {

ChiralVertex unit_charge(const context_base& context)
{
   std::array<int, 1> electron_indices = { 0 };
   std::array<int, 0> photon_indices = {};
   std::array<int, 2> indices = concatenate(photon_indices, electron_indices, electron_indices);

      const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   const std::complex<double> right = 0.7745966692414834*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   return {left, right};
}

} // namespace detail
} // namespace THDMIISNMSSMBCsimple_cxx_diagrams
} // namespace flexiblesusy
