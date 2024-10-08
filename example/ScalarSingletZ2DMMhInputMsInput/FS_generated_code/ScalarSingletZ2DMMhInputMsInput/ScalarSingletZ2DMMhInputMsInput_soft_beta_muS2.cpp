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


#include "ScalarSingletZ2DMMhInputMsInput_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of muS2.
 *
 * @return 1-loop beta function
 */
double ScalarSingletZ2DMMhInputMsInput_soft_parameters::calc_beta_muS2_1_loop(const Soft_traces& soft_traces) const
{


   double beta_muS2;

   beta_muS2 = Re(2*(2*LamSH*muH2 + 3*LamS*muS2));


   return oneLoop * beta_muS2;
}

/**
 * Calculates the 2-loop beta function of muS2.
 *
 * @return 2-loop beta function
 */
double ScalarSingletZ2DMMhInputMsInput_soft_parameters::calc_beta_muS2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_muS2;

   beta_muS2 = Re(0.4*(-60*LamSH*muH2*traceYdAdjYd - 20*LamSH*muH2*traceYeAdjYe
       - 60*LamSH*muH2*traceYuAdjYu + 12*LamSH*muH2*Sqr(g1) + 60*LamSH*muH2*Sqr
      (g2) - 75*muS2*Sqr(LamS) - 20*muH2*Sqr(LamSH) - 5*muS2*Sqr(LamSH)));


   return twoLoop * beta_muS2;
}

/**
 * Calculates the 3-loop beta function of muS2.
 *
 * @return 3-loop beta function
 */
double ScalarSingletZ2DMMhInputMsInput_soft_parameters::calc_beta_muS2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_muS2;

   beta_muS2 = 0;


   return threeLoop * beta_muS2;
}

/**
 * Calculates the 4-loop beta function of muS2.
 *
 * @return 4-loop beta function
 */
double ScalarSingletZ2DMMhInputMsInput_soft_parameters::calc_beta_muS2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_muS2;

   beta_muS2 = 0;


   return fourLoop * beta_muS2;
}

/**
 * Calculates the 5-loop beta function of muS2.
 *
 * @return 5-loop beta function
 */
double ScalarSingletZ2DMMhInputMsInput_soft_parameters::calc_beta_muS2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_muS2;

   beta_muS2 = 0;


   return fiveLoop * beta_muS2;
}

} // namespace flexiblesusy
