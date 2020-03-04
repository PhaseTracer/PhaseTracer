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

// File generated at Thu 7 Nov 2019 18:51:26

#include "THDMIISNMSSMBCsimple_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda3.
 *
 * @return 1-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda3_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambda3;

   beta_Lambda3 = Re(0.01*oneOver16PiSqr*(600*Lambda1*Lambda3 + 600*Lambda2*
      Lambda3 - 200*Lambda1*Lambda4 - 200*Lambda2*Lambda4 + 200*Lambda5*Lambda6
       + 600*Lambda3*traceYdAdjYd + 200*Lambda3*traceYeAdjYe + 600*Lambda3*
      traceYuAdjYu + 400*AbsSqr(Lambda7) + 27*Quad(g1) + 225*Quad(g2) - 180*
      Lambda3*Sqr(g1) - 900*Lambda3*Sqr(g2) + 90*Sqr(g1)*Sqr(g2) + 400*Sqr(
      Lambda3) + 200*Sqr(Lambda4)));


   return beta_Lambda3;
}

/**
 * Calculates the 2-loop beta function of Lambda3.
 *
 * @return 2-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda3_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYuYuAdjYdYdAdjYd = TRACE_STRUCT.
      traceYdAdjYuYuAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYuYuAdjYd = TRACE_STRUCT.
      traceYdAdjYuYuAdjYuYuAdjYd;


   double beta_Lambda3;

   beta_Lambda3 = Re(0.001*twoLoop*(16000*Lambda1*Lambda3*Lambda4 + 16000*
      Lambda2*Lambda3*Lambda4 - 8000*Lambda3*Lambda5*Lambda6 - 36000*Lambda1*
      Lambda3*traceYdAdjYd + 12000*Lambda1*Lambda4*traceYdAdjYd - 13500*Lambda3
      *traceYdAdjYdYdAdjYd - 9000*Lambda3*traceYdAdjYuYuAdjYd - 24000*Lambda4*
      traceYdAdjYuYuAdjYd + 12000*traceYdAdjYuYuAdjYdYdAdjYd + 12000*
      traceYdAdjYuYuAdjYuYuAdjYd - 12000*Lambda1*Lambda3*traceYeAdjYe + 4000*
      Lambda1*Lambda4*traceYeAdjYe - 4500*Lambda3*traceYeAdjYeYeAdjYe - 36000*
      Lambda2*Lambda3*traceYuAdjYu + 12000*Lambda2*Lambda4*traceYuAdjYu - 13500
      *Lambda3*traceYuAdjYuYuAdjYu - 8000*Lambda1*AbsSqr(Lambda7) - 8000*
      Lambda2*AbsSqr(Lambda7) - 4000*Lambda3*AbsSqr(Lambda7) - 24000*Lambda5*
      AbsSqr(Lambda7) - 24000*Lambda6*AbsSqr(Lambda7) - 12000*Cube(Lambda3) +
      12000*Cube(Lambda4) - 3537*Power6(g1) + 36375*Power6(g2) + 1350*Lambda1*
      Quad(g1) + 1350*Lambda2*Quad(g1) + 8865*Lambda3*Quad(g1) - 900*Lambda4*
      Quad(g1) + 450*traceYdAdjYd*Quad(g1) - 2250*traceYeAdjYe*Quad(g1) - 1710*
      traceYuAdjYu*Quad(g1) + 11250*Lambda1*Quad(g2) + 11250*Lambda2*Quad(g2) -
      13875*Lambda3*Quad(g2) - 7500*Lambda4*Quad(g2) - 2250*traceYdAdjYd*Quad(
      g2) - 750*traceYeAdjYe*Quad(g2) - 2250*traceYuAdjYu*Quad(g2) + 7200*
      Lambda1*Lambda3*Sqr(g1) + 7200*Lambda2*Lambda3*Sqr(g1) - 2400*Lambda1*
      Lambda4*Sqr(g1) - 2400*Lambda2*Lambda4*Sqr(g1) + 1250*Lambda3*
      traceYdAdjYd*Sqr(g1) + 3750*Lambda3*traceYeAdjYe*Sqr(g1) + 4250*Lambda3*
      traceYuAdjYu*Sqr(g1) - 7575*Quad(g2)*Sqr(g1) + 36000*Lambda1*Lambda3*Sqr(
      g2) + 36000*Lambda2*Lambda3*Sqr(g2) - 18000*Lambda1*Lambda4*Sqr(g2) -
      18000*Lambda2*Lambda4*Sqr(g2) + 12000*Lambda3*Lambda4*Sqr(g2) + 11250*
      Lambda3*traceYdAdjYd*Sqr(g2) + 3750*Lambda3*traceYeAdjYe*Sqr(g2) + 11250*
      Lambda3*traceYuAdjYu*Sqr(g2) - 8595*Quad(g1)*Sqr(g2) + 1500*Lambda1*Sqr(
      g1)*Sqr(g2) + 1500*Lambda2*Sqr(g1)*Sqr(g2) + 2850*Lambda3*Sqr(g1)*Sqr(g2)
      + 3000*Lambda4*Sqr(g1)*Sqr(g2) + 2700*traceYdAdjYd*Sqr(g1)*Sqr(g2) + 3300
      *traceYeAdjYe*Sqr(g1)*Sqr(g2) + 6300*traceYuAdjYu*Sqr(g1)*Sqr(g2) + 40000
      *Lambda3*traceYdAdjYd*Sqr(g3) + 40000*Lambda3*traceYuAdjYu*Sqr(g3) -
      15000*Lambda3*Sqr(Lambda1) + 4000*Lambda4*Sqr(Lambda1) - 15000*Lambda3*
      Sqr(Lambda2) + 4000*Lambda4*Sqr(Lambda2) - 36000*Lambda1*Sqr(Lambda3) -
      36000*Lambda2*Sqr(Lambda3) + 4000*Lambda4*Sqr(Lambda3) - 12000*
      traceYdAdjYd*Sqr(Lambda3) - 4000*traceYeAdjYe*Sqr(Lambda3) - 12000*
      traceYuAdjYu*Sqr(Lambda3) + 1200*Sqr(g1)*Sqr(Lambda3) + 6000*Sqr(g2)*Sqr(
      Lambda3) - 14000*Lambda1*Sqr(Lambda4) - 14000*Lambda2*Sqr(Lambda4) -
      16000*Lambda3*Sqr(Lambda4) - 6000*traceYdAdjYd*Sqr(Lambda4) - 2000*
      traceYeAdjYe*Sqr(Lambda4) - 6000*traceYuAdjYu*Sqr(Lambda4) + 2400*Sqr(g1)
      *Sqr(Lambda4) + 6000*Sqr(g2)*Sqr(Lambda4) - 1000*Lambda3*Sqr(Lambda5) -
      4000*Lambda6*Sqr(Lambda5) - 1000*Lambda3*Sqr(Lambda6) - 4000*Lambda5*Sqr(
      Lambda6)));


   return beta_Lambda3;
}

/**
 * Calculates the 3-loop beta function of Lambda3.
 *
 * @return 3-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda3_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda3;

   beta_Lambda3 = 0;


   return beta_Lambda3;
}

/**
 * Calculates the 4-loop beta function of Lambda3.
 *
 * @return 4-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda3_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda3;

   beta_Lambda3 = 0;


   return beta_Lambda3;
}

/**
 * Calculates the 5-loop beta function of Lambda3.
 *
 * @return 5-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda3_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda3;

   beta_Lambda3 = 0;


   return beta_Lambda3;
}

} // namespace flexiblesusy
