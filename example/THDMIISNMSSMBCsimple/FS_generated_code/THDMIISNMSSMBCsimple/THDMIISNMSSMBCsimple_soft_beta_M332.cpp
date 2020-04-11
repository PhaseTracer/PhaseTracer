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

// File generated at Sat 11 Apr 2020 12:52:54

#include "THDMIISNMSSMBCsimple_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

namespace {

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(double n, const Eigen::MatrixBase<Derived>& m)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m - Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(double n, const Eigen::MatrixBase<Derived>& m)
{
   return - m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

} // anonymous namespace

/**
 * Calculates the 1-loop beta function of M332.
 *
 * @return 1-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M332_1_loop(const Soft_traces& soft_traces) const
{


   double beta_M332;

   beta_M332 = Re(4*oneOver16PiSqr*(Lambda5*M112 + Lambda6*M222 + 2*Lambda8*
      M332 + AbsSqr(M123) + AbsSqr(M5)));


   return beta_M332;
}

/**
 * Calculates the 2-loop beta function of M332.
 *
 * @return 2-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M332_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_M332;

   beta_M332 = Re(0.4*twoLoop*(-60*Lambda5*M112*traceYdAdjYd - 20*Lambda5*M112*
      traceYeAdjYe - 60*Lambda6*M222*traceYuAdjYu - 40*M112*AbsSqr(Lambda7) -
      40*M222*AbsSqr(Lambda7) - 20*M332*AbsSqr(Lambda7) - 30*Lambda5*AbsSqr(
      M123) - 30*Lambda6*AbsSqr(M123) - 40*Lambda8*AbsSqr(M123) - 30*
      traceYdAdjYd*AbsSqr(M123) - 10*traceYeAdjYe*AbsSqr(M123) - 30*
      traceYuAdjYu*AbsSqr(M123) - 120*Lambda8*AbsSqr(M5) - 40*Lambda7*M5*Conj(
      M123) - 40*M123*Conj(Lambda7)*Conj(M5) + 12*Lambda5*M112*Sqr(g1) + 12*
      Lambda6*M222*Sqr(g1) + 12*AbsSqr(M123)*Sqr(g1) + 60*Lambda5*M112*Sqr(g2)
      + 60*Lambda6*M222*Sqr(g2) + 60*AbsSqr(M123)*Sqr(g2) - 20*M112*Sqr(Lambda5
      ) - 5*M332*Sqr(Lambda5) - 20*M222*Sqr(Lambda6) - 5*M332*Sqr(Lambda6) -
      100*M332*Sqr(Lambda8)));


   return beta_M332;
}

/**
 * Calculates the 3-loop beta function of M332.
 *
 * @return 3-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M332_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M332;

   beta_M332 = 0;


   return beta_M332;
}

/**
 * Calculates the 4-loop beta function of M332.
 *
 * @return 4-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M332_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M332;

   beta_M332 = 0;


   return beta_M332;
}

/**
 * Calculates the 5-loop beta function of M332.
 *
 * @return 5-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M332_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M332;

   beta_M332 = 0;


   return beta_M332;
}

} // namespace flexiblesusy
