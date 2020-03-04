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
 * Calculates the 1-loop beta function of M5.
 *
 * @return 1-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M5_1_loop(const Soft_traces& soft_traces) const
{


   double beta_M5;

   beta_M5 = Re(12*oneOver16PiSqr*(Lambda8*M5 + M123*Conj(Lambda7)));


   return beta_M5;
}

/**
 * Calculates the 2-loop beta function of M5.
 *
 * @return 2-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M5_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_M5;

   beta_M5 = Re(-0.6*twoLoop*(-20*M5*AbsSqr(Lambda7) + 60*Lambda5*M123*Conj(
      Lambda7) + 60*Lambda6*M123*Conj(Lambda7) + 80*Lambda8*M123*Conj(Lambda7)
      + 60*M123*traceYdAdjYd*Conj(Lambda7) + 20*M123*traceYeAdjYe*Conj(Lambda7)
      + 60*M123*traceYuAdjYu*Conj(Lambda7) - 24*M123*Conj(Lambda7)*Sqr(g1) -
      120*M123*Conj(Lambda7)*Sqr(g2) + 15*M5*Sqr(Lambda5) + 15*M5*Sqr(Lambda6)
      + 220*M5*Sqr(Lambda8)));


   return beta_M5;
}

/**
 * Calculates the 3-loop beta function of M5.
 *
 * @return 3-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M5_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M5;

   beta_M5 = 0;


   return beta_M5;
}

/**
 * Calculates the 4-loop beta function of M5.
 *
 * @return 4-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M5_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M5;

   beta_M5 = 0;


   return beta_M5;
}

/**
 * Calculates the 5-loop beta function of M5.
 *
 * @return 5-loop beta function
 */
double THDMIISNMSSMBCsimple_soft_parameters::calc_beta_M5_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M5;

   beta_M5 = 0;


   return beta_M5;
}

} // namespace flexiblesusy
