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

// File generated at Sat 11 Apr 2020 12:52:50

#include "THDMIISNMSSMBCsimple_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Yd.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Yd_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (oneOver16PiSqr*(-0.25*Yd*(-12*traceYdAdjYd - 4*traceYeAdjYe + Sqr
      (g1) + 9*Sqr(g2) + 32*Sqr(g3)) + 1.5*(Yd*Yd.adjoint()*Yd) + 0.5*(Yd*Yu.
      adjoint()*Yu))).real();


   return beta_Yd;
}

/**
 * Calculates the 2-loop beta function of Yd.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Yd_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (twoLoop*(-0.0016666666666666668*Yd*(600*Lambda3*Lambda4 + 4050*
      traceYdAdjYdYdAdjYd + 1350*traceYdAdjYuYuAdjYd + 1350*traceYeAdjYeYeAdjYe
       - 600*AbsSqr(Lambda7) + 113*Quad(g1) + 3150*Quad(g2) + 64800*Quad(g3) -
      375*traceYdAdjYd*Sqr(g1) - 1125*traceYeAdjYe*Sqr(g1) - 3375*traceYdAdjYd*
      Sqr(g2) - 1125*traceYeAdjYe*Sqr(g2) + 810*Sqr(g1)*Sqr(g2) - 12000*
      traceYdAdjYd*Sqr(g3) - 1240*Sqr(g1)*Sqr(g3) - 5400*Sqr(g2)*Sqr(g3) - 900*
      Sqr(Lambda1) - 600*Sqr(Lambda3) - 600*Sqr(Lambda4) - 300*Sqr(Lambda5)) +
      0.0125*(-480*Lambda1 - 540*traceYdAdjYd - 180*traceYeAdjYe + 187*Sqr(g1)
      + 675*Sqr(g2) + 1280*Sqr(g3))*(Yd*Yd.adjoint()*Yd) + 0.004166666666666667
      *(-480*Lambda3 + 960*Lambda4 - 540*traceYuAdjYu - 53*Sqr(g1) + 495*Sqr(g2
      ) + 1280*Sqr(g3))*(Yd*Yu.adjoint()*Yu) + 1.5*(Yd*Yd.adjoint()*Yd*Yd.
      adjoint()*Yd) - 0.25*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 0.25*(Yd*Yu.
      adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_Yd;
}

/**
 * Calculates the 3-loop beta function of Yd.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Yd_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return beta_Yd;
}

/**
 * Calculates the 4-loop beta function of Yd.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Yd_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return beta_Yd;
}

/**
 * Calculates the 5-loop beta function of Yd.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> THDMIISNMSSMBCsimple_susy_parameters::calc_beta_Yd_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return beta_Yd;
}

} // namespace flexiblesusy
