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
 * Calculates the 1-loop beta function of Lambda8.
 *
 * @return 1-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda8_1_loop(const Susy_traces& susy_traces) const
{


   double beta_Lambda8;

   beta_Lambda8 = Re(2*oneOver16PiSqr*(2*AbsSqr(Lambda7) + Sqr(Lambda5) + Sqr(
      Lambda6) + 10*Sqr(Lambda8)));


   return beta_Lambda8;
}

/**
 * Calculates the 2-loop beta function of Lambda8.
 *
 * @return 2-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda8_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambda8;

   beta_Lambda8 = Re(0.8*twoLoop*(-40*Lambda5*AbsSqr(Lambda7) - 40*Lambda6*
      AbsSqr(Lambda7) - 60*Lambda8*AbsSqr(Lambda7) - 15*traceYdAdjYd*AbsSqr(
      Lambda7) - 5*traceYeAdjYe*AbsSqr(Lambda7) - 15*traceYuAdjYu*AbsSqr(
      Lambda7) - 10*Cube(Lambda5) - 10*Cube(Lambda6) - 300*Cube(Lambda8) + 6*
      AbsSqr(Lambda7)*Sqr(g1) + 30*AbsSqr(Lambda7)*Sqr(g2) - 25*Lambda8*Sqr(
      Lambda5) - 15*traceYdAdjYd*Sqr(Lambda5) - 5*traceYeAdjYe*Sqr(Lambda5) + 3
      *Sqr(g1)*Sqr(Lambda5) + 15*Sqr(g2)*Sqr(Lambda5) - 25*Lambda8*Sqr(Lambda6)
      - 15*traceYuAdjYu*Sqr(Lambda6) + 3*Sqr(g1)*Sqr(Lambda6) + 15*Sqr(g2)*Sqr(
      Lambda6)));


   return beta_Lambda8;
}

/**
 * Calculates the 3-loop beta function of Lambda8.
 *
 * @return 3-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda8_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda8;

   beta_Lambda8 = 0;


   return beta_Lambda8;
}

/**
 * Calculates the 4-loop beta function of Lambda8.
 *
 * @return 4-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda8_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda8;

   beta_Lambda8 = 0;


   return beta_Lambda8;
}

/**
 * Calculates the 5-loop beta function of Lambda8.
 *
 * @return 5-loop beta function
 */
double THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Lambda8_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda8;

   beta_Lambda8 = 0;


   return beta_Lambda8;
}

} // namespace flexiblesusy
