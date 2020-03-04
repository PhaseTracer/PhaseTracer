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

// File generated at Thu 7 Nov 2019 18:51:27

#include "THDMIISNMSSMBCsimple_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda2.
 *
 * @return 1-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda2_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda2;

   beta_Lambda2 = Re(0.01*oneOver16PiSqr*(-400*Lambda3*Lambda4 + 1200*Lambda2*
      traceYuAdjYu - 1200*traceYuAdjYuYuAdjYu + 27*Quad(g1) + 225*Quad(g2) -
      180*Lambda2*Sqr(g1) - 900*Lambda2*Sqr(g2) + 90*Sqr(g1)*Sqr(g2) + 1200*Sqr
      (Lambda2) + 400*Sqr(Lambda3) + 200*Sqr(Lambda4) + 200*Sqr(Lambda6)));


   return beta_Lambda2;
}

/**
 * Calculates the 2-loop beta function of Lambda2.
 *
 * @return 2-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda2_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYuYuAdjYuYuAdjYd = TRACE_STRUCT.
      traceYdAdjYuYuAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYuYuAdjYu = TRACE_STRUCT.
      traceYuAdjYuYuAdjYuYuAdjYu;


   double beta_Lambda2;

   beta_Lambda2 = Re(0.001*twoLoop*(20000*Lambda2*Lambda3*Lambda4 + 24000*
      Lambda3*Lambda4*traceYdAdjYd - 9000*Lambda2*traceYdAdjYuYuAdjYd + 12000*
      traceYdAdjYuYuAdjYuYuAdjYd + 8000*Lambda3*Lambda4*traceYeAdjYe - 3000*
      Lambda2*traceYuAdjYuYuAdjYu + 60000*traceYuAdjYuYuAdjYuYuAdjYu + 4000*
      Lambda2*AbsSqr(Lambda7) - 8000*Lambda3*AbsSqr(Lambda7) - 16000*Lambda6*
      AbsSqr(Lambda7) - 78000*Cube(Lambda2) - 16000*Cube(Lambda3) + 12000*Cube(
      Lambda4) - 8000*Cube(Lambda6) - 3537*Power6(g1) + 36375*Power6(g2) + 9765
      *Lambda2*Quad(g1) + 1800*Lambda3*Quad(g1) - 900*Lambda4*Quad(g1) - 3420*
      traceYuAdjYu*Quad(g1) - 6375*Lambda2*Quad(g2) + 15000*Lambda3*Quad(g2) -
      7500*Lambda4*Quad(g2) - 4500*traceYuAdjYu*Quad(g2) - 4800*Lambda3*Lambda4
      *Sqr(g1) + 8500*Lambda2*traceYuAdjYu*Sqr(g1) - 3200*traceYuAdjYuYuAdjYu*
      Sqr(g1) - 7575*Quad(g2)*Sqr(g1) - 24000*Lambda3*Lambda4*Sqr(g2) + 22500*
      Lambda2*traceYuAdjYu*Sqr(g2) - 8595*Quad(g1)*Sqr(g2) + 5850*Lambda2*Sqr(
      g1)*Sqr(g2) + 3000*Lambda4*Sqr(g1)*Sqr(g2) + 12600*traceYuAdjYu*Sqr(g1)*
      Sqr(g2) + 80000*Lambda2*traceYuAdjYu*Sqr(g3) - 64000*traceYuAdjYuYuAdjYu*
      Sqr(g3) - 72000*traceYuAdjYu*Sqr(Lambda2) + 10800*Sqr(g1)*Sqr(Lambda2) +
      54000*Sqr(g2)*Sqr(Lambda2) - 20000*Lambda2*Sqr(Lambda3) + 24000*Lambda4*
      Sqr(Lambda3) - 24000*traceYdAdjYd*Sqr(Lambda3) - 8000*traceYeAdjYe*Sqr(
      Lambda3) + 4800*Sqr(g1)*Sqr(Lambda3) + 24000*Sqr(g2)*Sqr(Lambda3) - 12000
      *Lambda2*Sqr(Lambda4) - 32000*Lambda3*Sqr(Lambda4) - 12000*traceYdAdjYd*
      Sqr(Lambda4) - 4000*traceYeAdjYe*Sqr(Lambda4) + 2400*Sqr(g1)*Sqr(Lambda4)
      + 6000*Sqr(g2)*Sqr(Lambda4) - 10000*Lambda2*Sqr(Lambda6)));


   return beta_Lambda2;
}

/**
 * Calculates the 3-loop beta function of Lambda2.
 *
 * @return 3-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda2_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda2;

   beta_Lambda2 = 0;


   return beta_Lambda2;
}

/**
 * Calculates the 4-loop beta function of Lambda2.
 *
 * @return 4-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda2_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda2;

   beta_Lambda2 = 0;


   return beta_Lambda2;
}

/**
 * Calculates the 5-loop beta function of Lambda2.
 *
 * @return 5-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda2_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda2;

   beta_Lambda2 = 0;


   return beta_Lambda2;
}

} // namespace flexiblesusy
