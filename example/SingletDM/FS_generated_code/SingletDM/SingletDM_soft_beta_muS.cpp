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

// File generated at Tue 18 Aug 2020 22:46:29

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
 * Calculates the 1-loop beta function of muS.
 *
 * @return 1-loop beta function
 */
double SingletDM_soft_parameters::calc_beta_muS_1_loop(const Soft_traces& soft_traces) const
{


   double beta_muS;

   beta_muS = Re(2*(2*LamSH*muH + 3*LamS*muS)*oneOver16PiSqr);


   return beta_muS;
}

/**
 * Calculates the 2-loop beta function of muS.
 *
 * @return 2-loop beta function
 */
double SingletDM_soft_parameters::calc_beta_muS_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_muS;

   beta_muS = Re(0.4*twoLoop*(-60*LamSH*muH*traceYdAdjYd - 20*LamSH*muH*
      traceYeAdjYe - 60*LamSH*muH*traceYuAdjYu + 12*LamSH*muH*Sqr(g1) + 60*
      LamSH*muH*Sqr(g2) - 75*muS*Sqr(LamS) - 20*muH*Sqr(LamSH) - 5*muS*Sqr(
      LamSH)));


   return beta_muS;
}

/**
 * Calculates the 3-loop beta function of muS.
 *
 * @return 3-loop beta function
 */
double SingletDM_soft_parameters::calc_beta_muS_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_muS;

   beta_muS = 0;


   return beta_muS;
}

/**
 * Calculates the 4-loop beta function of muS.
 *
 * @return 4-loop beta function
 */
double SingletDM_soft_parameters::calc_beta_muS_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_muS;

   beta_muS = 0;


   return beta_muS;
}

/**
 * Calculates the 5-loop beta function of muS.
 *
 * @return 5-loop beta function
 */
double SingletDM_soft_parameters::calc_beta_muS_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_muS;

   beta_muS = 0;


   return beta_muS;
}

} // namespace flexiblesusy
