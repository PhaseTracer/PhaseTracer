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
 * @file ScalarSingletZ2DMMhInputMsInput_l_to_lgamma.cpp
 *
 * This file was generated with FlexibleSUSY 2.6.1 and SARAH 4.14.3 .
 */

#include <valarray>
#include <complex>

#include "ScalarSingletZ2DMMhInputMsInput_l_to_lgamma.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_mass_eigenstates.hpp"

#include "cxx_qft/ScalarSingletZ2DMMhInputMsInput_qft.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_FFV_form_factors.hpp"

#include "lowe.h"
#include "physical_input.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

using namespace ScalarSingletZ2DMMhInputMsInput_cxx_diagrams;
using namespace ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields;
using namespace ScalarSingletZ2DMMhInputMsInput_FFV_form_factors;

namespace ScalarSingletZ2DMMhInputMsInput_l_to_lgamma {

template <typename FIn, typename FOut, typename T1, typename T2>
double lepton_total_decay_width(
      T1 const& indices1, T2 const& indices2, 
      const ScalarSingletZ2DMMhInputMsInput_mass_eigenstates& model, const softsusy::QedQcd& qedqcd) {

   context_base context {model};

   const auto leptonInMass = context.mass<FIn>(indices1);
   const auto leptonOutMass = context.mass<FOut>(indices2);
   const auto r = Power2(leptonOutMass/leptonInMass);
   const auto GF2 = Power2(qedqcd.displayFermiConstant());
   const double Alpha = Power2(FIn::electric_charge * unit_charge(context))/(4*Pi);

   // eq. 10.6 of http://pdg.lbl.gov/2018/reviews/rpp2018-rev-standard-model.pdf
   return
      // tree-level (with massless outgoing lepton)
      GF2*Power5(leptonInMass)/(192.0*Power3(Pi))
      // outgoing lepton mass correction
      * (1 - 8*r + 8*Power3(r) - Power4(r) - 12*r*r*std::log(r));
      // higher order correction
      //* (1.0 + 0.5 * Alpha * (6.25-Sqr(Pi))/Pi)
}



}
} // namespace flexiblesusy
