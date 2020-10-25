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

// File generated at Sat 24 Oct 2020 17:07:49

#include "ScalarSingletZ2DM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of v.
 *
 * @return 1-loop beta function
 */
double ScalarSingletZ2DM_soft_parameters::calc_beta_v_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_v;

   beta_v = Re(0.2*oneOver16PiSqr*v*(-15*traceYdAdjYd - 5*traceYeAdjYe - 15*
      traceYuAdjYu + 3*Sqr(g1) + 15*Sqr(g2)));


   return beta_v;
}

/**
 * Calculates the 2-loop beta function of v.
 *
 * @return 2-loop beta function
 */
double ScalarSingletZ2DM_soft_parameters::calc_beta_v_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_v;

   beta_v = Re(-0.00125*twoLoop*v*(-5400*traceYdAdjYdYdAdjYd + 1200*
      traceYdAdjYuYuAdjYd - 1800*traceYeAdjYeYeAdjYe - 5400*traceYuAdjYuYuAdjYu
       + 1221*Quad(g1) - 9475*Quad(g2) + 1220*traceYdAdjYd*Sqr(g1) + 1740*
      traceYeAdjYe*Sqr(g1) + 2420*traceYuAdjYu*Sqr(g1) + 8100*traceYdAdjYd*Sqr(
      g2) + 2700*traceYeAdjYe*Sqr(g2) + 8100*traceYuAdjYu*Sqr(g2) - 450*Sqr(g1)
      *Sqr(g2) + 16000*traceYdAdjYd*Sqr(g3) + 16000*traceYuAdjYu*Sqr(g3) + 1200
      *Sqr(LamH) + 200*Sqr(LamSH)));


   return beta_v;
}

/**
 * Calculates the 3-loop beta function of v.
 *
 * @return 3-loop beta function
 */
double ScalarSingletZ2DM_soft_parameters::calc_beta_v_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v;

   beta_v = 0;


   return beta_v;
}

/**
 * Calculates the 4-loop beta function of v.
 *
 * @return 4-loop beta function
 */
double ScalarSingletZ2DM_soft_parameters::calc_beta_v_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v;

   beta_v = 0;


   return beta_v;
}

/**
 * Calculates the 5-loop beta function of v.
 *
 * @return 5-loop beta function
 */
double ScalarSingletZ2DM_soft_parameters::calc_beta_v_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v;

   beta_v = 0;


   return beta_v;
}

} // namespace flexiblesusy
