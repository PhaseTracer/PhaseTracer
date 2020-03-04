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

// File generated at Thu 7 Nov 2019 18:51:39

#include "THDMIISNMSSMBCsimple_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd THDMIISNMSSMBCsimple_input_parameters::get() const
{
   Eigen::ArrayXd pars(10);

   pars(0) = lambdaNMSSM;
   pars(1) = kappaNMSSM;
   pars(2) = AlambdaNMSSM;
   pars(3) = AkappaNMSSM;
   pars(4) = AtopNMSSM;
   pars(5) = mstopL;
   pars(6) = mstopR;
   pars(7) = MEWSB;
   pars(8) = TanBeta;
   pars(9) = vSIN;

   return pars;
}

void THDMIISNMSSMBCsimple_input_parameters::set(const Eigen::ArrayXd& pars)
{
   lambdaNMSSM = pars(0);
   kappaNMSSM = pars(1);
   AlambdaNMSSM = pars(2);
   AkappaNMSSM = pars(3);
   AtopNMSSM = pars(4);
   mstopL = pars(5);
   mstopR = pars(6);
   MEWSB = pars(7);
   TanBeta = pars(8);
   vSIN = pars(9);

}

std::ostream& operator<<(std::ostream& ostr, const THDMIISNMSSMBCsimple_input_parameters& input)
{
   ostr << "lambdaNMSSM = " << INPUT(lambdaNMSSM) << ", ";
   ostr << "kappaNMSSM = " << INPUT(kappaNMSSM) << ", ";
   ostr << "AlambdaNMSSM = " << INPUT(AlambdaNMSSM) << ", ";
   ostr << "AkappaNMSSM = " << INPUT(AkappaNMSSM) << ", ";
   ostr << "AtopNMSSM = " << INPUT(AtopNMSSM) << ", ";
   ostr << "mstopL = " << INPUT(mstopL) << ", ";
   ostr << "mstopR = " << INPUT(mstopR) << ", ";
   ostr << "MEWSB = " << INPUT(MEWSB) << ", ";
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "vSIN = " << INPUT(vSIN) << ", ";

   return ostr;
}

} // namespace flexiblesusy
