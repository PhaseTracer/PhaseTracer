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

// File generated at Sat 11 Apr 2020 12:52:55

#ifndef THDMIISNMSSMBCsimple_INPUT_PARAMETERS_H
#define THDMIISNMSSMBCsimple_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct THDMIISNMSSMBCsimple_input_parameters {
   double lambdaNMSSM{};
   double kappaNMSSM{};
   double AlambdaNMSSM{};
   double AkappaNMSSM{};
   double AtopNMSSM{};
   double mstopL{};
   double mstopR{};
   double MEWSB{};
   double TanBeta{};
   double vSIN{};


   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const THDMIISNMSSMBCsimple_input_parameters&);

} // namespace flexiblesusy

#endif
