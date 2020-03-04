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

// File generated at Thu 7 Nov 2019 18:53:49

#include "THDMIISNMSSMBCsimple_info.hpp"

#include <iostream>

namespace flexiblesusy {

namespace THDMIISNMSSMBCsimple_info {
   const double normalization_g1 = 0.7745966692414834;
   const double normalization_g2 = 1;
   const double normalization_g3 = 1;

   const std::array<int, NUMBER_OF_PARTICLES> particle_multiplicities = {1, 3, 3,
      3, 2, 3, 3, 3, 1, 1, 1};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_names = {"VG", "Fv"
      , "hh", "Ah", "Hm", "Fd", "Fu", "Fe", "VWm", "VP", "VZ"};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_latex_names = {
      "g", "\\nu", "h", "A^0", "H^-", "d", "u", "e", "W^-", "\\gamma", "Z"};

   const std::array<std::string, NUMBER_OF_PARAMETERS> parameter_names = {"g1",
      "g2", "g3", "Lambda7", "Lambda1", "Lambda4", "Lambda3", "Lambda2", "Lambda5"
      , "Lambda6", "Lambda8", "Yd(0,0)", "Yd(0,1)", "Yd(0,2)", "Yd(1,0)",
      "Yd(1,1)", "Yd(1,2)", "Yd(2,0)", "Yd(2,1)", "Yd(2,2)", "Ye(0,0)", "Ye(0,1)",
      "Ye(0,2)", "Ye(1,0)", "Ye(1,1)", "Ye(1,2)", "Ye(2,0)", "Ye(2,1)", "Ye(2,2)",
      "Yu(0,0)", "Yu(0,1)", "Yu(0,2)", "Yu(1,0)", "Yu(1,1)", "Yu(1,2)", "Yu(2,0)",
      "Yu(2,1)", "Yu(2,2)", "M123", "M5", "M112", "M222", "M332", "v1", "v2", "vS"
      };

   const std::array<std::string, NUMBER_OF_MIXINGS> particle_mixing_names = {
      "ZH(0,0)", "ZH(0,1)", "ZH(0,2)", "ZH(1,0)", "ZH(1,1)", "ZH(1,2)", "ZH(2,0)",
      "ZH(2,1)", "ZH(2,2)", "ZA(0,0)", "ZA(0,1)", "ZA(0,2)", "ZA(1,0)", "ZA(1,1)",
      "ZA(1,2)", "ZA(2,0)", "ZA(2,1)", "ZA(2,2)", "ZP(0,0)", "ZP(0,1)", "ZP(1,0)",
      "ZP(1,1)", "Re(Vd(0,0))", "Im(Vd(0,0))", "Re(Vd(0,1))", "Im(Vd(0,1))",
      "Re(Vd(0,2))", "Im(Vd(0,2))", "Re(Vd(1,0))", "Im(Vd(1,0))", "Re(Vd(1,1))",
      "Im(Vd(1,1))", "Re(Vd(1,2))", "Im(Vd(1,2))", "Re(Vd(2,0))", "Im(Vd(2,0))",
      "Re(Vd(2,1))", "Im(Vd(2,1))", "Re(Vd(2,2))", "Im(Vd(2,2))", "Re(Ud(0,0))",
      "Im(Ud(0,0))", "Re(Ud(0,1))", "Im(Ud(0,1))", "Re(Ud(0,2))", "Im(Ud(0,2))",
      "Re(Ud(1,0))", "Im(Ud(1,0))", "Re(Ud(1,1))", "Im(Ud(1,1))", "Re(Ud(1,2))",
      "Im(Ud(1,2))", "Re(Ud(2,0))", "Im(Ud(2,0))", "Re(Ud(2,1))", "Im(Ud(2,1))",
      "Re(Ud(2,2))", "Im(Ud(2,2))", "Re(Vu(0,0))", "Im(Vu(0,0))", "Re(Vu(0,1))",
      "Im(Vu(0,1))", "Re(Vu(0,2))", "Im(Vu(0,2))", "Re(Vu(1,0))", "Im(Vu(1,0))",
      "Re(Vu(1,1))", "Im(Vu(1,1))", "Re(Vu(1,2))", "Im(Vu(1,2))", "Re(Vu(2,0))",
      "Im(Vu(2,0))", "Re(Vu(2,1))", "Im(Vu(2,1))", "Re(Vu(2,2))", "Im(Vu(2,2))",
      "Re(Uu(0,0))", "Im(Uu(0,0))", "Re(Uu(0,1))", "Im(Uu(0,1))", "Re(Uu(0,2))",
      "Im(Uu(0,2))", "Re(Uu(1,0))", "Im(Uu(1,0))", "Re(Uu(1,1))", "Im(Uu(1,1))",
      "Re(Uu(1,2))", "Im(Uu(1,2))", "Re(Uu(2,0))", "Im(Uu(2,0))", "Re(Uu(2,1))",
      "Im(Uu(2,1))", "Re(Uu(2,2))", "Im(Uu(2,2))", "Re(Ve(0,0))", "Im(Ve(0,0))",
      "Re(Ve(0,1))", "Im(Ve(0,1))", "Re(Ve(0,2))", "Im(Ve(0,2))", "Re(Ve(1,0))",
      "Im(Ve(1,0))", "Re(Ve(1,1))", "Im(Ve(1,1))", "Re(Ve(1,2))", "Im(Ve(1,2))",
      "Re(Ve(2,0))", "Im(Ve(2,0))", "Re(Ve(2,1))", "Im(Ve(2,1))", "Re(Ve(2,2))",
      "Im(Ve(2,2))", "Re(Ue(0,0))", "Im(Ue(0,0))", "Re(Ue(0,1))", "Im(Ue(0,1))",
      "Re(Ue(0,2))", "Im(Ue(0,2))", "Re(Ue(1,0))", "Im(Ue(1,0))", "Re(Ue(1,1))",
      "Im(Ue(1,1))", "Re(Ue(1,2))", "Im(Ue(1,2))", "Re(Ue(2,0))", "Im(Ue(2,0))",
      "Re(Ue(2,1))", "Im(Ue(2,1))", "Re(Ue(2,2))", "Im(Ue(2,2))", "ZZ(0,0)",
      "ZZ(0,1)", "ZZ(1,0)", "ZZ(1,1)"};

   const std::array<std::string, NUMBER_OF_INPUT_PARAMETERS> input_parameter_names
       = {"lambdaNMSSM", "kappaNMSSM", "AlambdaNMSSM", "AkappaNMSSM", "AtopNMSSM",
      "mstopL", "mstopR", "MEWSB", "TanBeta", "vSIN"};

   const std::array<std::string, NUMBER_OF_EXTRA_PARAMETERS> extra_parameter_names
       = {};

   const std::string model_name = "THDMIISNMSSMBCsimple";

void print(std::ostream& ostr)
{
   ostr
      << "Model information\n"
      << "=================\n"
      << "Model name:                  " << model_name << '\n'
      << "Is a low-energy model:       "
      << (is_low_energy_model ? "yes" : "no") << '\n'
      << "Is a supersymmetric model:   "
      << (is_supersymmetric_model ? "yes" : "no") << '\n'
      << "Is a FlexibleEFTHiggs model: "
      << (is_FlexibleEFTHiggs ? "yes" : "no") << '\n'
      << "Number of multiplets:        " << NUMBER_OF_PARTICLES << '\n'
      << "Number of parameters:        " << NUMBER_OF_PARAMETERS << '\n'
      ;

   ostr << "\n"
      "Multiplets:                  ";
   for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
      ostr << particle_names[i]
           << '[' << particle_multiplicities[i] << ']';
      if (i + 1 < NUMBER_OF_PARTICLES)
         ostr << ", ";
   }

   ostr << "\n\n"
      "Parameters:                  ";
   for (int i = 0; i < NUMBER_OF_PARAMETERS; i++) {
      ostr << parameter_names[i];
      if (i + 1 < NUMBER_OF_PARAMETERS)
         ostr << ", ";
   }

   ostr << "\n\n"
      "Input parameters:            ";
   for (int i = 0; i < NUMBER_OF_INPUT_PARAMETERS; i++) {
      ostr << input_parameter_names[i];
      if (i + 1 < NUMBER_OF_INPUT_PARAMETERS)
         ostr << ", ";
   }

   ostr << '\n';
}

} // namespace THDMIISNMSSMBCsimple_info

} // namespace flexiblesusy

