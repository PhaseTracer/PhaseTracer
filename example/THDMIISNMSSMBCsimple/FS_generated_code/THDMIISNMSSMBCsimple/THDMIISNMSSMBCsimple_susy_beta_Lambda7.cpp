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
 * Calculates the 1-loop beta function of Lambda7.
 *
 * @return 1-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda7_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambda7;

   beta_Lambda7 = Re(-0.1*Lambda7*oneOver16PiSqr*(-20*Lambda3 - 20*Lambda4 - 40
      *Lambda5 - 40*Lambda6 - 40*Lambda8 - 30*traceYdAdjYd - 10*traceYeAdjYe -
      30*traceYuAdjYu + 9*Sqr(g1) + 45*Sqr(g2)));


   return beta_Lambda7;
}

/**
 * Calculates the 2-loop beta function of Lambda7.
 *
 * @return 2-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda7_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda7;

   beta_Lambda7 = Re(0.0025*Lambda7*twoLoop*(-2400*Lambda1*Lambda3 - 2400*
      Lambda2*Lambda3 - 2400*Lambda3*Lambda4 - 4800*Lambda1*Lambda5 - 6400*
      Lambda3*Lambda5 - 1600*Lambda4*Lambda5 - 4800*Lambda2*Lambda6 - 6400*
      Lambda3*Lambda6 - 1600*Lambda4*Lambda6 - 4000*Lambda5*Lambda6 - 12800*
      Lambda5*Lambda8 - 12800*Lambda6*Lambda8 - 2400*Lambda3*traceYdAdjYd -
      2400*Lambda4*traceYdAdjYd - 4800*Lambda5*traceYdAdjYd - 2700*
      traceYdAdjYdYdAdjYd - 9000*traceYdAdjYuYuAdjYd - 800*Lambda3*traceYeAdjYe
       - 800*Lambda4*traceYeAdjYe - 1600*Lambda5*traceYeAdjYe - 900*
      traceYeAdjYeYeAdjYe - 2400*Lambda3*traceYuAdjYu - 2400*Lambda4*
      traceYuAdjYu - 4800*Lambda6*traceYuAdjYu - 2700*traceYuAdjYuYuAdjYu +
      4000*AbsSqr(Lambda7) + 1377*Quad(g1) - 6075*Quad(g2) + 960*Lambda3*Sqr(g1
      ) + 960*Lambda4*Sqr(g1) + 240*Lambda5*Sqr(g1) + 240*Lambda6*Sqr(g1) + 250
      *traceYdAdjYd*Sqr(g1) + 750*traceYeAdjYe*Sqr(g1) + 850*traceYuAdjYu*Sqr(
      g1) + 4800*Lambda3*Sqr(g2) + 4800*Lambda4*Sqr(g2) + 1200*Lambda5*Sqr(g2)
      + 1200*Lambda6*Sqr(g2) + 2250*traceYdAdjYd*Sqr(g2) + 750*traceYeAdjYe*Sqr
      (g2) + 2250*traceYuAdjYu*Sqr(g2) + 450*Sqr(g1)*Sqr(g2) + 8000*
      traceYdAdjYd*Sqr(g3) + 8000*traceYuAdjYu*Sqr(g3) + 600*Sqr(Lambda1) + 600
      *Sqr(Lambda2) + 2400*Sqr(Lambda4) - 3800*Sqr(Lambda5) - 3800*Sqr(Lambda6)
      - 9600*Sqr(Lambda8)));


   return beta_Lambda7;
}

/**
 * Calculates the 3-loop beta function of Lambda7.
 *
 * @return 3-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda7_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda7;

   beta_Lambda7 = 0;


   return beta_Lambda7;
}

/**
 * Calculates the 4-loop beta function of Lambda7.
 *
 * @return 4-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda7_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda7;

   beta_Lambda7 = 0;


   return beta_Lambda7;
}

/**
 * Calculates the 5-loop beta function of Lambda7.
 *
 * @return 5-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda7_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda7;

   beta_Lambda7 = 0;


   return beta_Lambda7;
}

} // namespace flexiblesusy
