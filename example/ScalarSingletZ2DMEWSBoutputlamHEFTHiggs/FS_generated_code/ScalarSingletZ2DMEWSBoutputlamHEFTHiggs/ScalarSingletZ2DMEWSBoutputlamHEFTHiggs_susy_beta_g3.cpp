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

// File generated at Tue 17 Nov 2020 15:32:48

#include "ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g3.
 *
 * @return 1-loop beta function
 */
double ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_parameters::calc_beta_g3_1_loop(const Susy_traces& susy_traces) const
{


   double beta_g3;

   beta_g3 = Re(-7*oneOver16PiSqr*Cube(g3));


   return beta_g3;
}

/**
 * Calculates the 2-loop beta function of g3.
 *
 * @return 2-loop beta function
 */
double ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_parameters::calc_beta_g3_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g3;

   beta_g3 = Re(-0.1*twoLoop*Cube(g3)*(20*traceYdAdjYd + 20*traceYuAdjYu - 11*
      Sqr(g1) - 45*Sqr(g2) + 260*Sqr(g3)));


   return beta_g3;
}

/**
 * Calculates the 3-loop beta function of g3.
 *
 * @return 3-loop beta function
 */
double ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_parameters::calc_beta_g3_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g3;

   beta_g3 = 0;


   return beta_g3;
}

/**
 * Calculates the 4-loop beta function of g3.
 *
 * @return 4-loop beta function
 */
double ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_parameters::calc_beta_g3_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g3;

   beta_g3 = 0;


   return beta_g3;
}

/**
 * Calculates the 5-loop beta function of g3.
 *
 * @return 5-loop beta function
 */
double ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_parameters::calc_beta_g3_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g3;

   beta_g3 = 0;


   return beta_g3;
}

} // namespace flexiblesusy
