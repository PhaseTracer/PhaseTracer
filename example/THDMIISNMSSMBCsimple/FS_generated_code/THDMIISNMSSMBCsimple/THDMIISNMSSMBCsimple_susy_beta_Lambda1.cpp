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

// File generated at Sat 11 Apr 2020 12:52:48

#include "THDMIISNMSSMBCsimple_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda1.
 *
 * @return 1-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda1_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_Lambda1;

   beta_Lambda1 = Re(0.01*oneOver16PiSqr*(-400*Lambda3*Lambda4 + 1200*Lambda1*
      traceYdAdjYd - 1200*traceYdAdjYdYdAdjYd + 400*Lambda1*traceYeAdjYe - 400*
      traceYeAdjYeYeAdjYe + 27*Quad(g1) + 225*Quad(g2) - 180*Lambda1*Sqr(g1) -
      900*Lambda1*Sqr(g2) + 90*Sqr(g1)*Sqr(g2) + 1200*Sqr(Lambda1) + 400*Sqr(
      Lambda3) + 200*Sqr(Lambda4) + 200*Sqr(Lambda5)));


   return beta_Lambda1;
}

/**
 * Calculates the 2-loop beta function of Lambda1.
 *
 * @return 2-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda1_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYdAdjYdYdAdjYdYdAdjYd = TRACE_STRUCT.
      traceYdAdjYdYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYdYdAdjYd = TRACE_STRUCT.
      traceYdAdjYuYuAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYeYeAdjYe = TRACE_STRUCT.
      traceYeAdjYeYeAdjYeYeAdjYe;


   double beta_Lambda1;

   beta_Lambda1 = Re(0.001*twoLoop*(20000*Lambda1*Lambda3*Lambda4 - 3000*
      Lambda1*traceYdAdjYdYdAdjYd + 60000*traceYdAdjYdYdAdjYdYdAdjYd - 9000*
      Lambda1*traceYdAdjYuYuAdjYd + 12000*traceYdAdjYuYuAdjYdYdAdjYd - 1000*
      Lambda1*traceYeAdjYeYeAdjYe + 20000*traceYeAdjYeYeAdjYeYeAdjYe + 24000*
      Lambda3*Lambda4*traceYuAdjYu + 4000*Lambda1*AbsSqr(Lambda7) - 8000*
      Lambda3*AbsSqr(Lambda7) - 16000*Lambda5*AbsSqr(Lambda7) - 78000*Cube(
      Lambda1) - 16000*Cube(Lambda3) + 12000*Cube(Lambda4) - 8000*Cube(Lambda5)
      - 3537*Power6(g1) + 36375*Power6(g2) + 9765*Lambda1*Quad(g1) + 1800*
      Lambda3*Quad(g1) - 900*Lambda4*Quad(g1) + 900*traceYdAdjYd*Quad(g1) -
      4500*traceYeAdjYe*Quad(g1) - 6375*Lambda1*Quad(g2) + 15000*Lambda3*Quad(
      g2) - 7500*Lambda4*Quad(g2) - 4500*traceYdAdjYd*Quad(g2) - 1500*
      traceYeAdjYe*Quad(g2) - 4800*Lambda3*Lambda4*Sqr(g1) + 2500*Lambda1*
      traceYdAdjYd*Sqr(g1) + 1600*traceYdAdjYdYdAdjYd*Sqr(g1) + 7500*Lambda1*
      traceYeAdjYe*Sqr(g1) - 4800*traceYeAdjYeYeAdjYe*Sqr(g1) - 7575*Quad(g2)*
      Sqr(g1) - 24000*Lambda3*Lambda4*Sqr(g2) + 22500*Lambda1*traceYdAdjYd*Sqr(
      g2) + 7500*Lambda1*traceYeAdjYe*Sqr(g2) - 8595*Quad(g1)*Sqr(g2) + 5850*
      Lambda1*Sqr(g1)*Sqr(g2) + 3000*Lambda4*Sqr(g1)*Sqr(g2) + 5400*
      traceYdAdjYd*Sqr(g1)*Sqr(g2) + 6600*traceYeAdjYe*Sqr(g1)*Sqr(g2) + 80000*
      Lambda1*traceYdAdjYd*Sqr(g3) - 64000*traceYdAdjYdYdAdjYd*Sqr(g3) - 72000*
      traceYdAdjYd*Sqr(Lambda1) - 24000*traceYeAdjYe*Sqr(Lambda1) + 10800*Sqr(
      g1)*Sqr(Lambda1) + 54000*Sqr(g2)*Sqr(Lambda1) - 20000*Lambda1*Sqr(Lambda3
      ) + 24000*Lambda4*Sqr(Lambda3) - 24000*traceYuAdjYu*Sqr(Lambda3) + 4800*
      Sqr(g1)*Sqr(Lambda3) + 24000*Sqr(g2)*Sqr(Lambda3) - 12000*Lambda1*Sqr(
      Lambda4) - 32000*Lambda3*Sqr(Lambda4) - 12000*traceYuAdjYu*Sqr(Lambda4) +
      2400*Sqr(g1)*Sqr(Lambda4) + 6000*Sqr(g2)*Sqr(Lambda4) - 10000*Lambda1*Sqr
      (Lambda5)));


   return beta_Lambda1;
}

/**
 * Calculates the 3-loop beta function of Lambda1.
 *
 * @return 3-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda1_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda1;

   beta_Lambda1 = 0;


   return beta_Lambda1;
}

/**
 * Calculates the 4-loop beta function of Lambda1.
 *
 * @return 4-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda1_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda1;

   beta_Lambda1 = 0;


   return beta_Lambda1;
}

/**
 * Calculates the 5-loop beta function of Lambda1.
 *
 * @return 5-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda1_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda1;

   beta_Lambda1 = 0;


   return beta_Lambda1;
}

} // namespace flexiblesusy
