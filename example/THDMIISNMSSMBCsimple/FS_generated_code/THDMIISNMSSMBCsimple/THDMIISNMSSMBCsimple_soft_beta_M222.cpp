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

// File generated at Thu 7 Nov 2019 18:51:37

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
 * Calculates the 1-loop beta function of M222.
 *
 * @return 1-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M222_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_M222;

   beta_M222 = Re(0.1*oneOver16PiSqr*(40*Lambda3*M112 - 20*Lambda4*M112 + 60*
      Lambda2*M222 + 20*Lambda6*M332 + 60*M222*traceYuAdjYu + 20*AbsSqr(M123) -
      9*M222*Sqr(g1) - 45*M222*Sqr(g2)));


   return beta_M222;
}

/**
 * Calculates the 2-loop beta function of M222.
 *
 * @return 2-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M222_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_M222;

   beta_M222 = Re(0.0025*twoLoop*(3200*Lambda3*Lambda4*M112 + 800*Lambda3*
      Lambda4*M222 - 9600*Lambda3*M112*traceYdAdjYd + 4800*Lambda4*M112*
      traceYdAdjYd - 1800*M222*traceYdAdjYuYuAdjYd - 3200*Lambda3*M112*
      traceYeAdjYe + 1600*Lambda4*M112*traceYeAdjYe - 14400*Lambda2*M222*
      traceYuAdjYu - 5400*M222*traceYuAdjYuYuAdjYu - 1600*M112*AbsSqr(Lambda7)
      + 800*M222*AbsSqr(Lambda7) - 3200*M332*AbsSqr(Lambda7) - 2400*Lambda2*
      AbsSqr(M123) - 3200*Lambda3*AbsSqr(M123) - 800*Lambda4*AbsSqr(M123) -
      3200*Lambda6*AbsSqr(M123) - 2400*traceYdAdjYd*AbsSqr(M123) - 800*
      traceYeAdjYe*AbsSqr(M123) - 1600*Lambda6*AbsSqr(M5) - 1600*Lambda7*M5*
      Conj(M123) - 1600*M123*Conj(Lambda7)*Conj(M5) + 360*M112*Quad(g1) + 1737*
      M222*Quad(g1) + 3000*M112*Quad(g2) - 3075*M222*Quad(g2) + 1920*Lambda3*
      M112*Sqr(g1) - 960*Lambda4*M112*Sqr(g1) + 2880*Lambda2*M222*Sqr(g1) +
      1700*M222*traceYuAdjYu*Sqr(g1) + 120*AbsSqr(M123)*Sqr(g1) + 9600*Lambda3*
      M112*Sqr(g2) - 4800*Lambda4*M112*Sqr(g2) + 14400*Lambda2*M222*Sqr(g2) +
      4500*M222*traceYuAdjYu*Sqr(g2) + 600*AbsSqr(M123)*Sqr(g2) + 450*M222*Sqr(
      g1)*Sqr(g2) + 16000*M222*traceYuAdjYu*Sqr(g3) - 6000*M222*Sqr(Lambda2) -
      3200*M112*Sqr(Lambda3) - 800*M222*Sqr(Lambda3) - 3200*M112*Sqr(Lambda4) -
      800*M222*Sqr(Lambda4) - 400*M222*Sqr(Lambda6) - 1600*M332*Sqr(Lambda6)));


   return beta_M222;
}

/**
 * Calculates the 3-loop beta function of M222.
 *
 * @return 3-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M222_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M222;

   beta_M222 = 0;


   return beta_M222;
}

/**
 * Calculates the 4-loop beta function of M222.
 *
 * @return 4-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M222_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M222;

   beta_M222 = 0;


   return beta_M222;
}

/**
 * Calculates the 5-loop beta function of M222.
 *
 * @return 5-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M222_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M222;

   beta_M222 = 0;


   return beta_M222;
}

} // namespace flexiblesusy
