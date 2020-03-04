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
 * Calculates the 1-loop beta function of Lambda5.
 *
 * @return 1-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda5_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_Lambda5;

   beta_Lambda5 = Re(0.1*oneOver16PiSqr*(60*Lambda1*Lambda5 + 40*Lambda3*
      Lambda6 - 20*Lambda4*Lambda6 + 80*Lambda5*Lambda8 + 60*Lambda5*
      traceYdAdjYd + 20*Lambda5*traceYeAdjYe + 80*AbsSqr(Lambda7) - 9*Lambda5*
      Sqr(g1) - 45*Lambda5*Sqr(g2) + 40*Sqr(Lambda5)));


   return beta_Lambda5;
}

/**
 * Calculates the 2-loop beta function of Lambda5.
 *
 * @return 2-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda5_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_Lambda5;

   beta_Lambda5 = Re(0.0025*twoLoop*(800*Lambda3*Lambda4*Lambda5 + 3200*Lambda3
      *Lambda4*Lambda6 - 6400*Lambda3*Lambda5*Lambda6 + 3200*Lambda4*Lambda5*
      Lambda6 - 14400*Lambda1*Lambda5*traceYdAdjYd - 5400*Lambda5*
      traceYdAdjYdYdAdjYd - 1800*Lambda5*traceYdAdjYuYuAdjYd - 4800*Lambda1*
      Lambda5*traceYeAdjYe - 1800*Lambda5*traceYeAdjYeYeAdjYe - 9600*Lambda3*
      Lambda6*traceYuAdjYu + 4800*Lambda4*Lambda6*traceYuAdjYu - 9600*Lambda1*
      AbsSqr(Lambda7) - 12800*Lambda3*AbsSqr(Lambda7) - 3200*Lambda4*AbsSqr(
      Lambda7) - 15200*Lambda5*AbsSqr(Lambda7) - 8000*Lambda6*AbsSqr(Lambda7) -
      25600*Lambda8*AbsSqr(Lambda7) - 9600*traceYuAdjYu*AbsSqr(Lambda7) - 4400*
      Cube(Lambda5) + 1737*Lambda5*Quad(g1) + 360*Lambda6*Quad(g1) - 3075*
      Lambda5*Quad(g2) + 3000*Lambda6*Quad(g2) + 2880*Lambda1*Lambda5*Sqr(g1) +
      1920*Lambda3*Lambda6*Sqr(g1) - 960*Lambda4*Lambda6*Sqr(g1) + 500*Lambda5*
      traceYdAdjYd*Sqr(g1) + 1500*Lambda5*traceYeAdjYe*Sqr(g1) + 480*AbsSqr(
      Lambda7)*Sqr(g1) + 14400*Lambda1*Lambda5*Sqr(g2) + 9600*Lambda3*Lambda6*
      Sqr(g2) - 4800*Lambda4*Lambda6*Sqr(g2) + 4500*Lambda5*traceYdAdjYd*Sqr(g2
      ) + 1500*Lambda5*traceYeAdjYe*Sqr(g2) + 2400*AbsSqr(Lambda7)*Sqr(g2) +
      450*Lambda5*Sqr(g1)*Sqr(g2) + 16000*Lambda5*traceYdAdjYd*Sqr(g3) - 6000*
      Lambda5*Sqr(Lambda1) - 800*Lambda5*Sqr(Lambda3) - 3200*Lambda6*Sqr(
      Lambda3) - 800*Lambda5*Sqr(Lambda4) - 3200*Lambda6*Sqr(Lambda4) - 14400*
      Lambda1*Sqr(Lambda5) - 19200*Lambda8*Sqr(Lambda5) - 4800*traceYdAdjYd*Sqr
      (Lambda5) - 1600*traceYeAdjYe*Sqr(Lambda5) + 240*Sqr(g1)*Sqr(Lambda5) +
      1200*Sqr(g2)*Sqr(Lambda5) - 3200*Lambda3*Sqr(Lambda6) + 1600*Lambda4*Sqr(
      Lambda6) - 800*Lambda5*Sqr(Lambda6) - 16000*Lambda5*Sqr(Lambda8)));


   return beta_Lambda5;
}

/**
 * Calculates the 3-loop beta function of Lambda5.
 *
 * @return 3-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda5_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda5;

   beta_Lambda5 = 0;


   return beta_Lambda5;
}

/**
 * Calculates the 4-loop beta function of Lambda5.
 *
 * @return 4-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda5_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda5;

   beta_Lambda5 = 0;


   return beta_Lambda5;
}

/**
 * Calculates the 5-loop beta function of Lambda5.
 *
 * @return 5-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda5_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda5;

   beta_Lambda5 = 0;


   return beta_Lambda5;
}

} // namespace flexiblesusy
