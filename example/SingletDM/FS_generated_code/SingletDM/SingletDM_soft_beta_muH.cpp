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

// File generated at Tue 18 Aug 2020 22:46:30

#include "SingletDM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of muH.
 *
 * @return 1-loop beta function
 */
double SingletDM_soft_parameters::calc_beta_muH_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_muH;

   beta_muH = Re(0.1*oneOver16PiSqr*(60*LamH*muH + 10*LamSH*muS + 60*muH*
      traceYdAdjYd + 20*muH*traceYeAdjYe + 60*muH*traceYuAdjYu - 9*muH*Sqr(g1)
      - 45*muH*Sqr(g2)));


   return beta_muH;
}

/**
 * Calculates the 2-loop beta function of muH.
 *
 * @return 2-loop beta function
 */
double SingletDM_soft_parameters::calc_beta_muH_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_muH;

   beta_muH = Re(0.0025*twoLoop*(-14400*LamH*muH*traceYdAdjYd - 5400*muH*
      traceYdAdjYdYdAdjYd - 8400*muH*traceYdAdjYuYuAdjYd - 4800*LamH*muH*
      traceYeAdjYe - 1800*muH*traceYeAdjYeYeAdjYe - 14400*LamH*muH*traceYuAdjYu
       - 5400*muH*traceYuAdjYuYuAdjYu + 1671*muH*Quad(g1) - 3625*muH*Quad(g2) +
      2880*LamH*muH*Sqr(g1) + 500*muH*traceYdAdjYd*Sqr(g1) + 1500*muH*
      traceYeAdjYe*Sqr(g1) + 1700*muH*traceYuAdjYu*Sqr(g1) + 14400*LamH*muH*Sqr
      (g2) + 4500*muH*traceYdAdjYd*Sqr(g2) + 1500*muH*traceYeAdjYe*Sqr(g2) +
      4500*muH*traceYuAdjYu*Sqr(g2) + 450*muH*Sqr(g1)*Sqr(g2) + 16000*muH*
      traceYdAdjYd*Sqr(g3) + 16000*muH*traceYuAdjYu*Sqr(g3) - 6000*muH*Sqr(LamH
      ) - 200*muH*Sqr(LamSH) - 800*muS*Sqr(LamSH)));


   return beta_muH;
}

/**
 * Calculates the 3-loop beta function of muH.
 *
 * @return 3-loop beta function
 */
double SingletDM_soft_parameters::calc_beta_muH_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_muH;

   beta_muH = 0;


   return beta_muH;
}

/**
 * Calculates the 4-loop beta function of muH.
 *
 * @return 4-loop beta function
 */
double SingletDM_soft_parameters::calc_beta_muH_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_muH;

   beta_muH = 0;


   return beta_muH;
}

/**
 * Calculates the 5-loop beta function of muH.
 *
 * @return 5-loop beta function
 */
double SingletDM_soft_parameters::calc_beta_muH_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_muH;

   beta_muH = 0;


   return beta_muH;
}

} // namespace flexiblesusy
