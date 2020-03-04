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
 * Calculates the 1-loop beta function of M112.
 *
 * @return 1-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M112_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_M112;

   beta_M112 = Re(0.1*oneOver16PiSqr*(60*Lambda1*M112 + 40*Lambda3*M222 - 20*
      Lambda4*M222 + 20*Lambda5*M332 + 60*M112*traceYdAdjYd + 20*M112*
      traceYeAdjYe + 20*AbsSqr(M123) - 9*M112*Sqr(g1) - 45*M112*Sqr(g2)));


   return beta_M112;
}

/**
 * Calculates the 2-loop beta function of M112.
 *
 * @return 2-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M112_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_M112;

   beta_M112 = Re(0.0025*twoLoop*(800*Lambda3*Lambda4*M112 + 3200*Lambda3*
      Lambda4*M222 - 14400*Lambda1*M112*traceYdAdjYd - 5400*M112*
      traceYdAdjYdYdAdjYd - 1800*M112*traceYdAdjYuYuAdjYd - 4800*Lambda1*M112*
      traceYeAdjYe - 1800*M112*traceYeAdjYeYeAdjYe - 9600*Lambda3*M222*
      traceYuAdjYu + 4800*Lambda4*M222*traceYuAdjYu + 800*M112*AbsSqr(Lambda7)
      - 1600*M222*AbsSqr(Lambda7) - 3200*M332*AbsSqr(Lambda7) - 2400*Lambda1*
      AbsSqr(M123) - 3200*Lambda3*AbsSqr(M123) - 800*Lambda4*AbsSqr(M123) -
      3200*Lambda5*AbsSqr(M123) - 2400*traceYuAdjYu*AbsSqr(M123) - 1600*Lambda5
      *AbsSqr(M5) - 1600*Lambda7*M5*Conj(M123) - 1600*M123*Conj(Lambda7)*Conj(
      M5) + 1737*M112*Quad(g1) + 360*M222*Quad(g1) - 3075*M112*Quad(g2) + 3000*
      M222*Quad(g2) + 2880*Lambda1*M112*Sqr(g1) + 1920*Lambda3*M222*Sqr(g1) -
      960*Lambda4*M222*Sqr(g1) + 500*M112*traceYdAdjYd*Sqr(g1) + 1500*M112*
      traceYeAdjYe*Sqr(g1) + 120*AbsSqr(M123)*Sqr(g1) + 14400*Lambda1*M112*Sqr(
      g2) + 9600*Lambda3*M222*Sqr(g2) - 4800*Lambda4*M222*Sqr(g2) + 4500*M112*
      traceYdAdjYd*Sqr(g2) + 1500*M112*traceYeAdjYe*Sqr(g2) + 600*AbsSqr(M123)*
      Sqr(g2) + 450*M112*Sqr(g1)*Sqr(g2) + 16000*M112*traceYdAdjYd*Sqr(g3) -
      6000*M112*Sqr(Lambda1) - 800*M112*Sqr(Lambda3) - 3200*M222*Sqr(Lambda3) -
      800*M112*Sqr(Lambda4) - 3200*M222*Sqr(Lambda4) - 400*M112*Sqr(Lambda5) -
      1600*M332*Sqr(Lambda5)));


   return beta_M112;
}

/**
 * Calculates the 3-loop beta function of M112.
 *
 * @return 3-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M112_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M112;

   beta_M112 = 0;


   return beta_M112;
}

/**
 * Calculates the 4-loop beta function of M112.
 *
 * @return 4-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M112_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M112;

   beta_M112 = 0;


   return beta_M112;
}

/**
 * Calculates the 5-loop beta function of M112.
 *
 * @return 5-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M112_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M112;

   beta_M112 = 0;


   return beta_M112;
}

} // namespace flexiblesusy
