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

#ifndef ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INPUT_PARAMETERS_H
#define ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_input_parameters {
   double muH2Input{};
   double LamSHInput{};
   double LamSInput{};
   double muS2Input{};
   double QEWSB{};
   double Qin{};


   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_input_parameters&);

} // namespace flexiblesusy

#endif
