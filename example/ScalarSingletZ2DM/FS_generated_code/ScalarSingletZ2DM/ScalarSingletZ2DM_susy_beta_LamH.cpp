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

// File generated at Sat 24 Oct 2020 17:07:47

#include "ScalarSingletZ2DM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of LamH.
 *
 * @return 1-loop beta function
 */
double ScalarSingletZ2DM_susy_parameters::calc_beta_LamH_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_LamH;

   beta_LamH = Re(0.01*oneOver16PiSqr*(1200*LamH*traceYdAdjYd - 1200*
      traceYdAdjYdYdAdjYd + 400*LamH*traceYeAdjYe - 400*traceYeAdjYeYeAdjYe +
      1200*LamH*traceYuAdjYu - 1200*traceYuAdjYuYuAdjYu + 27*Quad(g1) + 225*
      Quad(g2) - 180*LamH*Sqr(g1) - 900*LamH*Sqr(g2) + 90*Sqr(g1)*Sqr(g2) +
      1200*Sqr(LamH) + 100*Sqr(LamSH)));


   return beta_LamH;
}

/**
 * Calculates the 2-loop beta function of LamH.
 *
 * @return 2-loop beta function
 */
double ScalarSingletZ2DM_susy_parameters::calc_beta_LamH_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYdYdAdjYdYdAdjYd = TRACE_STRUCT.
      traceYdAdjYdYdAdjYdYdAdjYd;
   const double traceYdAdjYdYdAdjYuYuAdjYd = TRACE_STRUCT.
      traceYdAdjYdYdAdjYuYuAdjYd;
   const double traceYdAdjYuYuAdjYdYdAdjYd = TRACE_STRUCT.
      traceYdAdjYuYuAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYuYuAdjYd = TRACE_STRUCT.
      traceYdAdjYuYuAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYeYeAdjYe = TRACE_STRUCT.
      traceYeAdjYeYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYuYuAdjYu = TRACE_STRUCT.
      traceYuAdjYuYuAdjYuYuAdjYu;


   double beta_LamH;

   beta_LamH = Re(0.001*twoLoop*(-3000*LamH*traceYdAdjYdYdAdjYd + 60000*
      traceYdAdjYdYdAdjYdYdAdjYd - 24000*traceYdAdjYdYdAdjYuYuAdjYd - 42000*
      LamH*traceYdAdjYuYuAdjYd + 12000*traceYdAdjYuYuAdjYdYdAdjYd - 12000*
      traceYdAdjYuYuAdjYuYuAdjYd - 1000*LamH*traceYeAdjYeYeAdjYe + 20000*
      traceYeAdjYeYeAdjYeYeAdjYe - 3000*LamH*traceYuAdjYuYuAdjYu + 60000*
      traceYuAdjYuYuAdjYuYuAdjYu - 78000*Cube(LamH) - 4000*Cube(LamSH) - 3411*
      Power6(g1) + 38125*Power6(g2) + 9435*LamH*Quad(g1) + 900*traceYdAdjYd*
      Quad(g1) - 4500*traceYeAdjYe*Quad(g1) - 3420*traceYuAdjYu*Quad(g1) - 9125
      *LamH*Quad(g2) - 4500*traceYdAdjYd*Quad(g2) - 1500*traceYeAdjYe*Quad(g2)
      - 4500*traceYuAdjYu*Quad(g2) + 2500*LamH*traceYdAdjYd*Sqr(g1) + 1600*
      traceYdAdjYdYdAdjYd*Sqr(g1) + 7500*LamH*traceYeAdjYe*Sqr(g1) - 4800*
      traceYeAdjYeYeAdjYe*Sqr(g1) + 8500*LamH*traceYuAdjYu*Sqr(g1) - 3200*
      traceYuAdjYuYuAdjYu*Sqr(g1) - 7225*Quad(g2)*Sqr(g1) + 22500*LamH*
      traceYdAdjYd*Sqr(g2) + 7500*LamH*traceYeAdjYe*Sqr(g2) + 22500*LamH*
      traceYuAdjYu*Sqr(g2) - 8385*Quad(g1)*Sqr(g2) + 5850*LamH*Sqr(g1)*Sqr(g2)
      + 5400*traceYdAdjYd*Sqr(g1)*Sqr(g2) + 6600*traceYeAdjYe*Sqr(g1)*Sqr(g2) +
      12600*traceYuAdjYu*Sqr(g1)*Sqr(g2) + 80000*LamH*traceYdAdjYd*Sqr(g3) -
      64000*traceYdAdjYdYdAdjYd*Sqr(g3) + 80000*LamH*traceYuAdjYu*Sqr(g3) -
      64000*traceYuAdjYuYuAdjYu*Sqr(g3) - 72000*traceYdAdjYd*Sqr(LamH) - 24000*
      traceYeAdjYe*Sqr(LamH) - 72000*traceYuAdjYu*Sqr(LamH) + 10800*Sqr(g1)*Sqr
      (LamH) + 54000*Sqr(g2)*Sqr(LamH) - 5000*LamH*Sqr(LamSH)));


   return beta_LamH;
}

/**
 * Calculates the 3-loop beta function of LamH.
 *
 * @return 3-loop beta function
 */
double ScalarSingletZ2DM_susy_parameters::calc_beta_LamH_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamH;

   beta_LamH = 0;


   return beta_LamH;
}

/**
 * Calculates the 4-loop beta function of LamH.
 *
 * @return 4-loop beta function
 */
double ScalarSingletZ2DM_susy_parameters::calc_beta_LamH_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamH;

   beta_LamH = 0;


   return beta_LamH;
}

/**
 * Calculates the 5-loop beta function of LamH.
 *
 * @return 5-loop beta function
 */
double ScalarSingletZ2DM_susy_parameters::calc_beta_LamH_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamH;

   beta_LamH = 0;


   return beta_LamH;
}

} // namespace flexiblesusy
