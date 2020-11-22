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


#include "ScalarSingletZ2DMMhInput_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of LamSH.
 *
 * @return 1-loop beta function
 */
double ScalarSingletZ2DMMhInput_susy_parameters::calc_beta_LamSH_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_LamSH;

   beta_LamSH = Re(0.1*LamSH*(60*LamH + 60*LamS + 40*LamSH + 60*traceYdAdjYd +
      20*traceYeAdjYe + 60*traceYuAdjYu - 9*Sqr(g1) - 45*Sqr(g2)));


   return oneLoop * beta_LamSH;
}

/**
 * Calculates the 2-loop beta function of LamSH.
 *
 * @return 2-loop beta function
 */
double ScalarSingletZ2DMMhInput_susy_parameters::calc_beta_LamSH_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_LamSH;

   beta_LamSH = Re(-0.0025*LamSH*(14400*LamH*LamSH + 14400*LamS*LamSH + 14400*
      LamH*traceYdAdjYd + 4800*LamSH*traceYdAdjYd + 5400*traceYdAdjYdYdAdjYd +
      8400*traceYdAdjYuYuAdjYd + 4800*LamH*traceYeAdjYe + 1600*LamSH*
      traceYeAdjYe + 1800*traceYeAdjYeYeAdjYe + 14400*LamH*traceYuAdjYu + 4800*
      LamSH*traceYuAdjYu + 5400*traceYuAdjYuYuAdjYu - 1671*Quad(g1) + 3625*Quad
      (g2) - 2880*LamH*Sqr(g1) - 240*LamSH*Sqr(g1) - 500*traceYdAdjYd*Sqr(g1) -
      1500*traceYeAdjYe*Sqr(g1) - 1700*traceYuAdjYu*Sqr(g1) - 14400*LamH*Sqr(g2
      ) - 1200*LamSH*Sqr(g2) - 4500*traceYdAdjYd*Sqr(g2) - 1500*traceYeAdjYe*
      Sqr(g2) - 4500*traceYuAdjYu*Sqr(g2) - 450*Sqr(g1)*Sqr(g2) - 16000*
      traceYdAdjYd*Sqr(g3) - 16000*traceYuAdjYu*Sqr(g3) + 6000*Sqr(LamH) +
      12000*Sqr(LamS) + 4200*Sqr(LamSH)));


   return twoLoop * beta_LamSH;
}

/**
 * Calculates the 3-loop beta function of LamSH.
 *
 * @return 3-loop beta function
 */
double ScalarSingletZ2DMMhInput_susy_parameters::calc_beta_LamSH_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamSH;

   beta_LamSH = 0;


   return threeLoop * beta_LamSH;
}

/**
 * Calculates the 4-loop beta function of LamSH.
 *
 * @return 4-loop beta function
 */
double ScalarSingletZ2DMMhInput_susy_parameters::calc_beta_LamSH_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamSH;

   beta_LamSH = 0;


   return fourLoop * beta_LamSH;
}

/**
 * Calculates the 5-loop beta function of LamSH.
 *
 * @return 5-loop beta function
 */
double ScalarSingletZ2DMMhInput_susy_parameters::calc_beta_LamSH_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamSH;

   beta_LamSH = 0;


   return fiveLoop * beta_LamSH;
}

} // namespace flexiblesusy
