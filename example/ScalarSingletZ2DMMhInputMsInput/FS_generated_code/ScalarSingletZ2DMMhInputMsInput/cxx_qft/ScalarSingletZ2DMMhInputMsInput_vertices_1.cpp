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


/**
 * @file cxx_qft/ScalarSingletZ2DMMhInputMsInput_vertices.cpp
 *
 * This file was generated with FlexibleSUSY 2.5.0 and SARAH 4.14.3 .
 */

#include "ScalarSingletZ2DMMhInputMsInput_context_base.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_input_parameters.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_vertices.hpp"

#include "concatenate.hpp"
#include "wrappers.hpp"

#define INPUTPARAMETER(p) context.model.get_input_parameters().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

namespace flexiblesusy {
namespace ScalarSingletZ2DMMhInputMsInput_cxx_diagrams {
namespace detail {

ChiralVertex VertexImpl<fields::Ah, typename fields::bar<fields::Fe>::type, fields::Fe>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2));

   return {left, right};
}

ChiralVertex VertexImpl<fields::hh, typename fields::bar<fields::Fe>::type, fields::Fe>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2));

   return {left, right};
}

MomentumDifferenceVertex VertexImpl<fields::Hp, typename fields::conj<fields::Hp>::type, fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*(-3.872983346207417*g1*Cos(ThetaW) - 5*g2*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::Ah>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::hh>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   const std::complex<double> right = 0.7745966692414834*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, typename fields::conj<fields::Hp>::type, fields::Fv>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ue = MODELPARAMETER(Ue);

   const std::complex<double> left = -SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,gt2));

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Hp, fields::Fe>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ue = MODELPARAMETER(Ue);

   const std::complex<double> left = 0;

   const std::complex<double> right = -SUM(j1,0,2,Conj(Ye(j1,gt1))*Ue(gt2,j1));

   return {left, right};
}

} // namespace detail
} // namespace ScalarSingletZ2DMMhInputMsInput_cxx_diagrams
} // namespace flexiblesusy
