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

// File generated at Tue 17 Nov 2020 15:32:52

#include "ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_input_parameters::get() const
{
   Eigen::ArrayXd pars(6);

   pars(0) = muH2Input;
   pars(1) = LamSHInput;
   pars(2) = LamSInput;
   pars(3) = muS2Input;
   pars(4) = QEWSB;
   pars(5) = Qin;

   return pars;
}

void ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_input_parameters::set(const Eigen::ArrayXd& pars)
{
   muH2Input = pars(0);
   LamSHInput = pars(1);
   LamSInput = pars(2);
   muS2Input = pars(3);
   QEWSB = pars(4);
   Qin = pars(5);

}

std::ostream& operator<<(std::ostream& ostr, const ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_input_parameters& input)
{
   ostr << "muH2Input = " << INPUT(muH2Input) << ", ";
   ostr << "LamSHInput = " << INPUT(LamSHInput) << ", ";
   ostr << "LamSInput = " << INPUT(LamSInput) << ", ";
   ostr << "muS2Input = " << INPUT(muS2Input) << ", ";
   ostr << "QEWSB = " << INPUT(QEWSB) << ", ";
   ostr << "Qin = " << INPUT(Qin) << ", ";

   return ostr;
}

} // namespace flexiblesusy
