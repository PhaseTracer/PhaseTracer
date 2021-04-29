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
 * @file ScalarSingletZ2DMMhInputMsInput_b_to_s_gamma.cpp
 *
 * This file was generated with FlexibleSUSY 2.5.0 and SARAH 4.14.3 .
 */

#include <array>
#include <complex>
#include <iostream>

#include "ScalarSingletZ2DMMhInputMsInput_b_to_s_gamma.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_mass_eigenstates.hpp"

#include "cxx_qft/ScalarSingletZ2DMMhInputMsInput_qft.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_FFV_form_factors.hpp"

#include "lowe.h"
#include "wrappers.hpp"

#include <fstream>

#define MODELPARAMETER(p) context.model.get_##p()

namespace flexiblesusy {

using namespace ScalarSingletZ2DMMhInputMsInput_cxx_diagrams;
using namespace ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields;
using namespace ScalarSingletZ2DMMhInputMsInput_FFV_form_factors;

namespace ScalarSingletZ2DMMhInputMsInput_b_to_s_gamma {

/**
 * @brief Writes the Wilson coefficients into a json file to make it usable
 *        for the wilson and flavio package
 */
void write_wilsoncoeffs(const std::complex<double>& C7NP_bs, const std::complex<double>& C7pNP_bs,
   const std::complex<double>& C8NP_bs, const std::complex<double>& C8pNP_bs,
   const double& matching_scale)
{
   std::ofstream wc_json;
   wc_json.open ("WC_ScalarSingletZ2DMMhInputMsInput.json");
   wc_json << "{\n\"eft\": \"WET\",\n";
   wc_json << "\t\"basis\": \"flavio\",\n";
   wc_json << "\t\"scale\": \"" << matching_scale << "\",\n";
   wc_json << "\t\"values\": { \n";
   wc_json << "\t\t\"C7_bs\": { \n";
   wc_json << "\t\t\t\"Re\": " << Re(C7NP_bs) << ",\n";
   wc_json << "\t\t\t\"Im\": " << Im(C7NP_bs) << "\n";
   wc_json << "\t\t},\n";
   wc_json << "\t\t\"C7p_bs\": { \n";
   wc_json << "\t\t\t\"Re\": " << Re(C7pNP_bs) << ",\n";
   wc_json << "\t\t\t\"Im\": " << Im(C7pNP_bs) << "\n";
   wc_json << "\t\t},\n";
   wc_json << "\t\t\"C8_bs\": { \n";
   wc_json << "\t\t\t\"Re\": " << Re(C8NP_bs) << ",\n";
   wc_json << "\t\t\t\"Im\": " << Im(C8NP_bs) << "\n";
   wc_json << "\t\t},\n";
   wc_json << "\t\t\"C8p_bs\": { \n";
   wc_json << "\t\t\t\"Re\": " << Re(C8pNP_bs) << ",\n";
   wc_json << "\t\t\t\"Im\": " << Im(C8pNP_bs) << "\n";
   wc_json << "\t\t}\n";
   wc_json << "\t}\n";
   wc_json << "}";
   wc_json.close();
}

/**
 * @brief Calculates the Wilson coefficients C7, C7p, C8 and C8p
 *        which are required to compute the branching ratio
 *        BR(B->Xs gamma)
 */
std::array<std::complex<double>, 4> calculate_b_to_s_gamma(
   const ScalarSingletZ2DMMhInputMsInput_mass_eigenstates& model, const softsusy::QedQcd& qedqcd)
{
   context_base context { model };

   const auto complex_CKM = qedqcd.get_complex_ckm();
   const std::complex<double> Vtb (complex_CKM(2, 2));
   std::complex<double> Vts;
   // avoid division by zero
   if (is_zero(std::abs(complex_CKM(2, 1)))) {
       Vts = std::complex<double>(-0.0404063888996, -0.00072107415577);
   } else {
       Vts = complex_CKM(2, 1);
   }
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);

   return {0, 0, 0, 0};
}

} // namespace ScalarSingletZ2DMMhInputMsInput_b_to_s_gamma
} // namespace flexiblesusy
