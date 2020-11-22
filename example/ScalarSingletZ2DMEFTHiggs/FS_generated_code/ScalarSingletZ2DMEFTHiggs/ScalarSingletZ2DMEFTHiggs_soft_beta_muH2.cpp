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


#include "ScalarSingletZ2DMEFTHiggs_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of muH2.
 *
 * @return 1-loop beta function
 */
double ScalarSingletZ2DMEFTHiggs_soft_parameters::calc_beta_muH2_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_muH2;

   beta_muH2 = Re(0.1*(60*LamH*muH2 + 10*LamSH*muS2 + 60*muH2*traceYdAdjYd + 20
      *muH2*traceYeAdjYe + 60*muH2*traceYuAdjYu - 9*muH2*Sqr(g1) - 45*muH2*Sqr(
      g2)));


   return oneLoop * beta_muH2;
}

/**
 * Calculates the 2-loop beta function of muH2.
 *
 * @return 2-loop beta function
 */
double ScalarSingletZ2DMEFTHiggs_soft_parameters::calc_beta_muH2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_muH2;

   beta_muH2 = Re(0.0025*(-14400*LamH*muH2*traceYdAdjYd - 5400*muH2*
      traceYdAdjYdYdAdjYd - 8400*muH2*traceYdAdjYuYuAdjYd - 4800*LamH*muH2*
      traceYeAdjYe - 1800*muH2*traceYeAdjYeYeAdjYe - 14400*LamH*muH2*
      traceYuAdjYu - 5400*muH2*traceYuAdjYuYuAdjYu + 1671*muH2*Quad(g1) - 3625*
      muH2*Quad(g2) + 2880*LamH*muH2*Sqr(g1) + 500*muH2*traceYdAdjYd*Sqr(g1) +
      1500*muH2*traceYeAdjYe*Sqr(g1) + 1700*muH2*traceYuAdjYu*Sqr(g1) + 14400*
      LamH*muH2*Sqr(g2) + 4500*muH2*traceYdAdjYd*Sqr(g2) + 1500*muH2*
      traceYeAdjYe*Sqr(g2) + 4500*muH2*traceYuAdjYu*Sqr(g2) + 450*muH2*Sqr(g1)*
      Sqr(g2) + 16000*muH2*traceYdAdjYd*Sqr(g3) + 16000*muH2*traceYuAdjYu*Sqr(
      g3) - 6000*muH2*Sqr(LamH) - 200*muH2*Sqr(LamSH) - 800*muS2*Sqr(LamSH)));


   return twoLoop * beta_muH2;
}

/**
 * Calculates the 3-loop beta function of muH2.
 *
 * @return 3-loop beta function
 */
double ScalarSingletZ2DMEFTHiggs_soft_parameters::calc_beta_muH2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_muH2;

   beta_muH2 = 0;


   return threeLoop * beta_muH2;
}

/**
 * Calculates the 4-loop beta function of muH2.
 *
 * @return 4-loop beta function
 */
double ScalarSingletZ2DMEFTHiggs_soft_parameters::calc_beta_muH2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_muH2;

   beta_muH2 = 0;


   return fourLoop * beta_muH2;
}

/**
 * Calculates the 5-loop beta function of muH2.
 *
 * @return 5-loop beta function
 */
double ScalarSingletZ2DMEFTHiggs_soft_parameters::calc_beta_muH2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_muH2;

   beta_muH2 = 0;


   return fiveLoop * beta_muH2;
}

} // namespace flexiblesusy
