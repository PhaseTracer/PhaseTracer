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

// File generated at Thu 7 Nov 2019 18:51:36

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
 * Calculates the 1-loop beta function of M123.
 *
 * @return 1-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M123_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_M123;

   beta_M123 = Re(0.1*oneOver16PiSqr*(20*Lambda3*M123 + 20*Lambda4*M123 + 20*
      Lambda5*M123 + 20*Lambda6*M123 + 40*Lambda7*M5 + 30*M123*traceYdAdjYd +
      10*M123*traceYeAdjYe + 30*M123*traceYuAdjYu - 9*M123*Sqr(g1) - 45*M123*
      Sqr(g2)));


   return beta_M123;
}

/**
 * Calculates the 2-loop beta function of M123.
 *
 * @return 2-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M123_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_M123;

   beta_M123 = Re(0.0025*twoLoop*(-2400*Lambda1*Lambda3*M123 - 2400*Lambda2*
      Lambda3*M123 - 2400*Lambda3*Lambda4*M123 - 2400*Lambda1*Lambda5*M123 -
      3200*Lambda3*Lambda5*M123 - 800*Lambda4*Lambda5*M123 - 2400*Lambda2*
      Lambda6*M123 - 3200*Lambda3*Lambda6*M123 - 800*Lambda4*Lambda6*M123 -
      2400*Lambda5*Lambda6*M123 - 3200*Lambda5*Lambda8*M123 - 3200*Lambda6*
      Lambda8*M123 - 4800*Lambda5*Lambda7*M5 - 4800*Lambda6*Lambda7*M5 - 6400*
      Lambda7*Lambda8*M5 - 2400*Lambda3*M123*traceYdAdjYd - 2400*Lambda4*M123*
      traceYdAdjYd - 2400*Lambda5*M123*traceYdAdjYd - 2700*M123*
      traceYdAdjYdYdAdjYd - 6600*M123*traceYdAdjYuYuAdjYd - 800*Lambda3*M123*
      traceYeAdjYe - 800*Lambda4*M123*traceYeAdjYe - 800*Lambda5*M123*
      traceYeAdjYe - 900*M123*traceYeAdjYeYeAdjYe - 2400*Lambda3*M123*
      traceYuAdjYu - 2400*Lambda4*M123*traceYuAdjYu - 2400*Lambda6*M123*
      traceYuAdjYu - 2700*M123*traceYuAdjYuYuAdjYu - 10400*M123*AbsSqr(Lambda7)
      + 1377*M123*Quad(g1) - 6075*M123*Quad(g2) + 960*Lambda3*M123*Sqr(g1) +
      960*Lambda4*M123*Sqr(g1) + 120*Lambda5*M123*Sqr(g1) + 120*Lambda6*M123*
      Sqr(g1) + 250*M123*traceYdAdjYd*Sqr(g1) + 750*M123*traceYeAdjYe*Sqr(g1) +
      850*M123*traceYuAdjYu*Sqr(g1) + 4800*Lambda3*M123*Sqr(g2) + 4800*Lambda4*
      M123*Sqr(g2) + 600*Lambda5*M123*Sqr(g2) + 600*Lambda6*M123*Sqr(g2) + 2250
      *M123*traceYdAdjYd*Sqr(g2) + 750*M123*traceYeAdjYe*Sqr(g2) + 2250*M123*
      traceYuAdjYu*Sqr(g2) + 450*M123*Sqr(g1)*Sqr(g2) + 8000*M123*traceYdAdjYd*
      Sqr(g3) + 8000*M123*traceYuAdjYu*Sqr(g3) + 600*M123*Sqr(Lambda1) + 600*
      M123*Sqr(Lambda2) + 2400*M123*Sqr(Lambda4) - 200*M123*Sqr(Lambda5) - 200*
      M123*Sqr(Lambda6) + 1600*M123*Sqr(Lambda8)));


   return beta_M123;
}

/**
 * Calculates the 3-loop beta function of M123.
 *
 * @return 3-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M123_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M123;

   beta_M123 = 0;


   return beta_M123;
}

/**
 * Calculates the 4-loop beta function of M123.
 *
 * @return 4-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M123_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M123;

   beta_M123 = 0;


   return beta_M123;
}

/**
 * Calculates the 5-loop beta function of M123.
 *
 * @return 5-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M123_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M123;

   beta_M123 = 0;


   return beta_M123;
}

} // namespace flexiblesusy
