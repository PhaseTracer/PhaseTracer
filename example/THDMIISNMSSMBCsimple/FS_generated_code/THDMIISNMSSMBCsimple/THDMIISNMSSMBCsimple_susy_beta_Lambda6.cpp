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

// File generated at Sat 11 Apr 2020 12:52:49

#include "THDMIISNMSSMBCsimple_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda6.
 *
 * @return 1-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda6_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambda6;

   beta_Lambda6 = Re(0.1*oneOver16PiSqr*(40*Lambda3*Lambda5 - 20*Lambda4*
      Lambda5 + 60*Lambda2*Lambda6 + 80*Lambda6*Lambda8 + 60*Lambda6*
      traceYuAdjYu + 80*AbsSqr(Lambda7) - 9*Lambda6*Sqr(g1) - 45*Lambda6*Sqr(g2
      ) + 40*Sqr(Lambda6)));


   return beta_Lambda6;
}

/**
 * Calculates the 2-loop beta function of Lambda6.
 *
 * @return 2-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda6_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda6;

   beta_Lambda6 = Re(0.0025*twoLoop*(3200*Lambda3*Lambda4*Lambda5 + 800*Lambda3
      *Lambda4*Lambda6 - 6400*Lambda3*Lambda5*Lambda6 + 3200*Lambda4*Lambda5*
      Lambda6 - 9600*Lambda3*Lambda5*traceYdAdjYd + 4800*Lambda4*Lambda5*
      traceYdAdjYd - 1800*Lambda6*traceYdAdjYuYuAdjYd - 3200*Lambda3*Lambda5*
      traceYeAdjYe + 1600*Lambda4*Lambda5*traceYeAdjYe - 14400*Lambda2*Lambda6*
      traceYuAdjYu - 5400*Lambda6*traceYuAdjYuYuAdjYu - 9600*Lambda2*AbsSqr(
      Lambda7) - 12800*Lambda3*AbsSqr(Lambda7) - 3200*Lambda4*AbsSqr(Lambda7) -
      8000*Lambda5*AbsSqr(Lambda7) - 15200*Lambda6*AbsSqr(Lambda7) - 25600*
      Lambda8*AbsSqr(Lambda7) - 9600*traceYdAdjYd*AbsSqr(Lambda7) - 3200*
      traceYeAdjYe*AbsSqr(Lambda7) - 4400*Cube(Lambda6) + 360*Lambda5*Quad(g1)
      + 1737*Lambda6*Quad(g1) + 3000*Lambda5*Quad(g2) - 3075*Lambda6*Quad(g2) +
      1920*Lambda3*Lambda5*Sqr(g1) - 960*Lambda4*Lambda5*Sqr(g1) + 2880*Lambda2
      *Lambda6*Sqr(g1) + 1700*Lambda6*traceYuAdjYu*Sqr(g1) + 480*AbsSqr(Lambda7
      )*Sqr(g1) + 9600*Lambda3*Lambda5*Sqr(g2) - 4800*Lambda4*Lambda5*Sqr(g2) +
      14400*Lambda2*Lambda6*Sqr(g2) + 4500*Lambda6*traceYuAdjYu*Sqr(g2) + 2400*
      AbsSqr(Lambda7)*Sqr(g2) + 450*Lambda6*Sqr(g1)*Sqr(g2) + 16000*Lambda6*
      traceYuAdjYu*Sqr(g3) - 6000*Lambda6*Sqr(Lambda2) - 3200*Lambda5*Sqr(
      Lambda3) - 800*Lambda6*Sqr(Lambda3) - 3200*Lambda5*Sqr(Lambda4) - 800*
      Lambda6*Sqr(Lambda4) - 3200*Lambda3*Sqr(Lambda5) + 1600*Lambda4*Sqr(
      Lambda5) - 800*Lambda6*Sqr(Lambda5) - 14400*Lambda2*Sqr(Lambda6) - 19200*
      Lambda8*Sqr(Lambda6) - 4800*traceYuAdjYu*Sqr(Lambda6) + 240*Sqr(g1)*Sqr(
      Lambda6) + 1200*Sqr(g2)*Sqr(Lambda6) - 16000*Lambda6*Sqr(Lambda8)));


   return beta_Lambda6;
}

/**
 * Calculates the 3-loop beta function of Lambda6.
 *
 * @return 3-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda6_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda6;

   beta_Lambda6 = 0;


   return beta_Lambda6;
}

/**
 * Calculates the 4-loop beta function of Lambda6.
 *
 * @return 4-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda6_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda6;

   beta_Lambda6 = 0;


   return beta_Lambda6;
}

/**
 * Calculates the 5-loop beta function of Lambda6.
 *
 * @return 5-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda6_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda6;

   beta_Lambda6 = 0;


   return beta_Lambda6;
}

} // namespace flexiblesusy
