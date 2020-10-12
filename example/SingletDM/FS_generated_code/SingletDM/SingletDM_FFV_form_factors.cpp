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
 * @file SingletDM_FFV_form_factors.cpp
 *
 * This file was generated at Tue 18 Aug 2020 22:47:57 with FlexibleSUSY
 * 2.4.2 and SARAH 4.14.3 .
 */

#include <complex>
#include <valarray>

#include "SingletDM_mass_eigenstates.hpp"
#include "concatenate.hpp"
#include "cxx_qft/SingletDM_qft.hpp"

#include "wrappers.hpp"

using namespace flexiblesusy;
using namespace SingletDM_cxx_diagrams;
using namespace SingletDM_cxx_diagrams::fields;

namespace {

static constexpr double oneOver32PiSqr = 0.5*oneOver16PiSqr;

/**
 * @class FFV_SSF
 * @brief A template that calculate contributions to the FFV form
 *        factors of a given particles in a one loop diagram
 *        specified by a vector emitters and an exchange particle.
 * @tparam Args Specifies in order the field of which to
 *              calculate the electric dipole moment,
 *              the photon emitter and the exchange particle
 *              in a one-loop diagram where the photon emitter
 *              is a scalar and the exchange particle a fermion.
 *
 * This template evaluates the contribution to the electric
 * dipole moment of a one loop diagram with fields given by
 * \a Args.
 */
template <class Fj, class Fi, class V, class S1, class S2, class F>
struct FFV_SSF {
   static std::valarray<std::complex<double>>
   value(const typename field_indices<Fj>::type& indices_in,
         const typename field_indices<Fi>::type& indices_out,
         context_base const& context,
         bool discard_SM_contributions);
};

/**
* @class FFV_FFS
* @brief A template that calculate contributions to the FFV form
*        factors of a given particle in a one loop diagram
*        specified by a vector emitters and an exchange particle.
* @tparam Args Specifies in order the field of which to
*              calculate the electric dipole moment,
*              the photon emitter and the exchange particle
*              in a one-loop diagram where the photon emitter
*              is a fermion and the exchange particle a scalar.
*
* This template evaluates the contribution to the electric
* dipole moment of a one loop diagram with fields given by
* \a Args.
*/
template <class Fj, class Fi, class V, class F1, class F2, class S>
struct FFV_FFS {
   static std::valarray<std::complex<double>>
   value(const typename field_indices<Fj>::type& indices_in,
         const typename field_indices<Fi>::type& indices_out,
         context_base const& context,
         bool discard_SM_contributions);
};

} // anonymous namespace

namespace flexiblesusy {
namespace SingletDM_FFV_form_factors {
std::valarray<std::complex<double>> calculate_Fe_Fe_VP_form_factors (
   int generationIndex1, int generationIndex2,
   const SingletDM_mass_eigenstates& model, bool discard_SM_contributions) {

   context_base context {model};
   std::array<int, 1> indices1 = {generationIndex1 };
   std::array<int, 1> indices2 = {generationIndex2 };

   std::valarray<std::complex<double>> val {0.0, 0.0, 0.0, 0.0};

   val += std::complex<double> {1., 0.} * FFV_FFS<Fe,Fe,VP,Fe,Fe,Ah>::value(indices1, indices2, context, discard_SM_contributions);
   val += std::complex<double> {1., 0.} * FFV_FFS<Fe,Fe,VP,Fe,Fe,hh>::value(indices1, indices2, context, discard_SM_contributions);
   val += std::complex<double> {1., 0.} * FFV_SSF<Fe,Fe,VP,typename conj<Hp>::type,typename conj<Hp>::type,Fv>::value(indices1, indices2, context, discard_SM_contributions);

   return val;
}


}
} // namespace flexiblesusy

namespace {
/**
* @defgroup LoopFunctions Loop functions
* @brief The loop functions necessary for the Fe_I -> Fe_J gamma one-loop calculations.
*
* These are OneLoopFunctionA(), OneLoopFunctionB()
* as specified in arXiv:0808.1819
*/

// function from eq. 15 of hep-ph/9510309
double OneLoopFunctionA(double r)
{
   if (is_zero(1.0 - r)) {
      return 1.5;
   } else if (is_zero(r)) {
      return 2.0;
   } else {
      return (2.0 - 9.0 * r + 18.0 * r * r - 11.0 * r * r * r +
              6.0 * r * r * r * std::log(r)) /
             pow(1.0 - r, 4);
   }
}

double OneLoopFunctionB(double r)
{
   const double y = r - 1.0;
   if (is_zero(r)) {
      return 2.0;
   } else if (std::abs(y) < 0.23) {
      // error around x=1 is <= 10^-12 on an intel i7
      return (1.0000000000000000000 -
              0.4000000000000000000  * y +
              0.2000000000000000000  * y * y -
              0.11428571428571428571 * y * y * y +
              0.07142857142857142857 * y * y * y * y -
              0.04761904761904761905 * y * y * y * y * y +
              0.03333333333333333333 * y * y * y * y * y * y -
              0.02424242424242424242 * y * y * y * y * y * y * y +
              0.0181818181818181818  * y * y * y * y * y * y * y * y -
              0.01398601398601398601 * y * y * y * y * y * y * y * y * y +
              0.01098901098901098901 * y * y * y * y * y * y * y * y * y * y -
              0.0087912087912087912  * y * y * y * y * y * y * y * y * y * y * y +
              0.00714285714285714286 * y * y * y * y * y * y * y * y * y * y * y * y -
              0.0058823529411764706  * y * y * y * y * y * y * y * y * y * y * y * y * y +
              0.0049019607843137255  * y * y * y * y * y * y * y * y * y * y * y * y * y * y -
              0.0041279669762641899  * y * y * y * y * y * y * y * y * y * y * y * y * y * y * y);
   } else {
      return 2. *
             (1.0 - 6.0 * r + 3.0 * r * r + 2.0 * r * r * r -
              6.0 * r * r * std::log(r)) /
             pow(1.0 - r, 4);
   }
}

double OneLoopFunctionC(double r)
{
   const double y = r - 1.0;
   if (is_zero(r)) {
      return 3.0;
   } else if (std::abs(y) < 0.185) {
      // error around x=1 is <= 10^-13 on an intel i7
      return (1.0000000000000000000 -
              0.50000000000000000000 * y +
              0.30000000000000000000 * y * y -
              0.2000000000000000000  * y * y * y +
              0.14285714285714285714 * y * y * y * y -
              0.10714285714285714286 * y * y * y * y * y +
              0.08333333333333333333 * y * y * y * y * y * y -
              0.06666666666666666667 * y * y * y * y * y * y * y +
              0.05454545454545454545 * y * y * y * y * y * y * y * y -
              0.0454545454545454545  * y * y * y * y * y * y * y * y * y +
              0.0384615384615384615  * y * y * y * y * y * y * y * y * y * y -
              0.03296703296703296703 * y * y * y * y * y * y * y * y * y * y * y +
              0.0285714285714285714  * y * y * y * y * y * y * y * y * y * y * y * y -
              0.02500000000000000000 * y * y * y * y * y * y * y * y * y * y * y * y * y +
              0.0220588235294117647  * y * y * y * y * y * y * y * y * y * y * y * y * y * y -
              0.0196078431372549020  * y * y * y * y * y * y * y * y * y * y * y * y * y * y * y);
   } else {
      return 3. * (1.0 - r * r + 2.0 * r * std::log(r)) / pow(1.0 - r, 3);
   }
}

double OneLoopFunctionD(double r)
{
   if (is_zero(1.0 - r)) {
      return -9.0 / 2.0;
   } else {
      return (16.0 - 45.0 * r + 36.0 * r * r - 7.0 * r * r * r +
              6.0 * (2.0 - 3.0 * r) * std::log(r)) /
             pow(1.0 - r, 4);
   }
}

double OneLoopFunctionE(double r)
{
   const double y = r - 1.0;
   if (is_zero(r)) {
      return 12.0;
   } else if (std::abs(y) < 0.21) {
      // error around x=1 is <= 10^-12 on an intel i7
      return (1.0000000000000000000 -
              0.60000000000000000000  * y +
              0.40000000000000000000  * y * y -
              0.28571428571428571429  * y * y * y +
              0.21428571428571428571  * y * y * y * y -
              0.16666666666666666667  * y * y * y * y * y +
              0.13333333333333333333  * y * y * y * y * y * y -
              0.10909090909090909091  * y * y * y * y * y * y * y +
              0.090909090909090909091 * y * y * y * y * y * y * y * y -
              0.076923076923076923077 * y * y * y * y * y * y * y * y * y +
              0.065934065934065934066 * y * y * y * y * y * y * y * y * y * y -
              0.057142857142857142857 * y * y * y * y * y * y * y * y * y * y * y +
              0.050000000000000000000 * y * y * y * y * y * y * y * y * y * y * y * y -
              0.044117647058823529412 * y * y * y * y * y * y * y * y * y * y * y * y * y +
              0.039215686274509803922 * y * y * y * y * y * y * y * y * y * y * y * y * y * y -
              0.035087719298245614035 * y * y * y * y * y * y * y * y * y * y * y * y * y * y * y);
   } else {
      return 2. *
             (2.0 + 3.0 * r - 6.0 * r * r + r * r * r + 6.0 * r * std::log(r)) /
             pow(1.0 - r, 4);
   }
}

double OneLoopFunctionF(double r)
{
   const double y = r - 1.0;
   if (std::abs(y) < 0.155) {
      // error around x=1 is <= 10^-13 on an intel i7
      return (1.0 - 
              0.75 * y + 
              0.6 * y * y -
              0.50000000000000000000 * y * y * y +
              0.4285714285714285714  * y * y * y * y -
              0.37500000000000000000 * y * y * y * y * y +
              0.33333333333333333333 * y * y * y * y * y * y -
              0.3000000000000000000  * y * y * y * y * y * y * y +
              0.2727272727272727273  * y * y * y * y * y * y * y * y -
              0.2500000000000000000  * y * y * y * y * y * y * y * y * y +
              0.23076923076923076923 * y * y * y * y * y * y * y * y * y * y -
              0.21428571428571428571 * y * y * y * y * y * y * y * y * y * y * y +
              0.2000000000000000000  * y * y * y * y * y * y * y * y * y * y * y * y -
              0.1875000000000000000  * y * y * y * y * y * y * y * y * y * y * y * y * y +
              0.1764705882352941176  * y * y * y * y * y * y * y * y * y * y * y * y * y * y -
              0.16666666666666666667 * y * y * y * y * y * y * y * y * y * y * y * y * y * y * y);
   } else {
      return 3. / 2. * (-3.0 + 4.0 * r - r * r - 2.0 * std::log(r)) /
             pow(1.0 - r, 3);
   }
}

// emit massless vector boson from the internal scalar line
template <class Fj, class Fi, class V, class SA, class SB, class F>
std::valarray<std::complex<double>> FFV_SSF<Fj, Fi, V, SA, SB, F>::value(
   const typename field_indices<Fj>::type& indices_in,
   const typename field_indices<Fi>::type& indices_out,
   context_base const& context,
   bool discard_SM_contributions)
{

   static_assert(
      std::is_same<SA, SB>::type::value,
      "Internal scalars in the FFV_SSF instantiation must be of the same type."
   );

   using VertexFBarFjSBar = Vertex<typename F::lorentz_conjugate, typename SA::lorentz_conjugate, Fj>;
   using VertexFiBarFS    = Vertex<typename Fi::lorentz_conjugate, SB, F>;
   using VertexSBarSVBar  = Vertex<typename SB::lorentz_conjugate, SA, typename V::lorentz_conjugate>;

   // masses of external fermions
   const auto mj = context.mass<Fj>(indices_in);
   const auto mi = context.mass<Fi>(indices_out);

   // loop over all possible particle generations attached to both vertices
   std::valarray<std::complex<double>> res {0.0, 0.0, 0.0, 0.0};
   for (const auto& indexIn: index_range<VertexFBarFjSBar>()) {
      for (const auto& indexOut: index_range<VertexFiBarFS>()) {

         // cycle if generations of external fermions  are different then requested   
         const auto jFieldIndices = VertexFBarFjSBar::template indices_of_field<2>(indexIn);
         const auto iFieldIndices = VertexFiBarFS::template indices_of_field<0>(indexOut);
         if (jFieldIndices != indices_in || iFieldIndices != indices_out)
            continue;

         // match indices of the fermion in the loop
         const auto fermionFieldIndicesIn = VertexFBarFjSBar::template indices_of_field<0>(indexIn);
         const auto fermionFieldIndicesOut = VertexFiBarFS::template indices_of_field<2>(indexOut);
         if (fermionFieldIndicesIn != fermionFieldIndicesOut)
            continue;

         // match indices of the scalar in the loop
         const auto scalarFieldIndicesIn = VertexFBarFjSBar::template indices_of_field<1>(indexIn);
         const auto scalarIndicesOut = VertexFiBarFS::template indices_of_field<1>(indexOut);
         if (scalarFieldIndicesIn != scalarIndicesOut) 
            continue;

         if (discard_SM_contributions) {
            if (isSMField<SA>(scalarFieldIndicesIn) && isSMField<F>(fermionFieldIndicesIn)) {
               continue;
            }
         }

         const auto vertexIn = VertexFBarFjSBar::evaluate(indexIn, context);
         const auto vertexOut = VertexFiBarFS::evaluate(indexOut, context);

         const auto indexEmit = concatenate(scalarFieldIndicesIn, scalarFieldIndicesIn);
         const auto vertexEmit = VertexSBarSVBar::evaluate(indexEmit, context);

         const auto mS = context.mass<SA>(scalarFieldIndicesIn);
         const auto mF = context.mass<F>(fermionFieldIndicesIn);
         const auto x = pow(mF/mS, 2);

         // TODO: check the sign convention of this coupling
         std::complex<double> vector_boson_coupling {-vertexEmit.value(1,0)};

         // eq. 15 of hep-ph/9510309 (possibly with different sign)
         const std::complex<double> A1L =
            - 1./18. * vertexOut.right() * vertexIn.left() * OneLoopFunctionA(x);
         // eq. 16 of hep-ph/9510309 (possibly with different sign)
         const std::complex<double> A2L = 
            - vertexOut.left() * vertexIn.right() * OneLoopFunctionB(x)/12.
            - vertexOut.left()* vertexIn.left() * mF/mj * OneLoopFunctionC(x)/3.
            - mi/mj * vertexOut.right() * vertexIn.left() * OneLoopFunctionB(x)/12.; 

         // eq. 15 & 16 of hep-ph/9510309 after replacement L <-> R (possibly with different sign)
         const std::complex<double> A1R = 
            - 1./18. * vertexOut.left() * vertexIn.right() * OneLoopFunctionA(x);
         const std::complex<double> A2R = 
            - vertexOut.right() * vertexIn.left() * OneLoopFunctionB(x)/12. 
            - vertexOut.right()* vertexIn.right() * mF/mj * OneLoopFunctionC(x)/3.
            - mi/mj * vertexOut.left() * vertexIn.right() * OneLoopFunctionB(x)/12.; 

         const std::complex<double> massFactor = pow(mS,-2);

         res += oneOver32PiSqr * vector_boson_coupling * massFactor
            * std::valarray<std::complex<double>> {A1L, A1R, A2L, A2R};
      }
   }

   return res;
}

// emit massless vector boson from the internal fermion line
template <class Fj, class Fi, class V, class FA, class FB, class S>
std::valarray<std::complex<double>> FFV_FFS<Fj, Fi, V, FA, FB, S>::value(
   const typename field_indices<Fj>::type& indices_in,
   const typename field_indices<Fi>::type& indices_out,
   context_base const& context,
   bool discard_SM_contributions)
{

   static_assert(
      std::is_same<FA, FB>::type::value, 
      "Internal fermions in the FFV_FFS instantiation must be of the same type."
   );

   using VertexFBarFjSBar = Vertex<typename S::lorentz_conjugate, typename FA::lorentz_conjugate, Fj>;
   using VertexFiBarFS    = Vertex<typename Fi::lorentz_conjugate, FB, S>;
   using VertexFBarFVBar  = Vertex<typename FB::lorentz_conjugate, FA, typename V::lorentz_conjugate>;
   
   // masses of external fermions
   const auto mj = context.mass<Fj>(indices_in);
   const auto mi = context.mass<Fi>(indices_out);

   // loop over all possible particle generations attached to both vertices
   std::valarray<std::complex<double>> res {0.0, 0.0, 0.0, 0.0};
   for (const auto& indexIn: index_range<VertexFBarFjSBar>()) {
      for (const auto& indexOut: index_range<VertexFiBarFS>()) {

         // cycle if generations of external fermions are different then requested   
         const auto jFieldIndices = VertexFBarFjSBar::template indices_of_field<2>(indexIn);
         const auto iFieldIndices = VertexFiBarFS::template indices_of_field<0>(indexOut);
         if (jFieldIndices != indices_in || iFieldIndices != indices_out)
            continue;

         const auto fermionFieldIndicesIn = VertexFBarFjSBar::template indices_of_field<1>(indexIn);
         const auto fermionFieldIndicesOut = VertexFiBarFS::template indices_of_field<1>(indexOut);
         if (fermionFieldIndicesIn != fermionFieldIndicesOut )
            continue;

         const auto scalarFieldIndicesIn = VertexFBarFjSBar::template indices_of_field<0>(indexIn);
         const auto scalarIndicesOut = VertexFiBarFS::template indices_of_field<2>(indexOut);
         if (scalarFieldIndicesIn != scalarIndicesOut) 
            continue;

         if (discard_SM_contributions) {
            if (isSMField<S>(scalarFieldIndicesIn) && isSMField<FA>(fermionFieldIndicesIn)) {
               continue;
            }
         }

         const auto vertexIn = VertexFBarFjSBar::evaluate(indexIn, context);
         const auto vertexOut = VertexFiBarFS::evaluate(indexOut, context);
         
         const auto indexEmit = concatenate(fermionFieldIndicesIn, fermionFieldIndicesIn);
         const auto vertexEmit = VertexFBarFVBar::evaluate(indexEmit, context);

         const auto mF = context.mass<FA>(fermionFieldIndicesIn);
         const auto mS = context.mass<S>(scalarFieldIndicesIn);
         const auto x = pow(mF/mS, 2);

         std::complex<double> vector_boson_coupling {vertexEmit.left()};

         // eq. 18 of hep-ph/9510309 (possibly with different sign)
         const std::complex<double> A1L =
            - 1./18. * vertexOut.right() * vertexIn.left() * OneLoopFunctionD(x);
         // eq. 19 of hep-ph/9510309 (possibly with different sign)
         const std::complex<double> A2L = 
            - vertexOut.left() * vertexIn.right() * OneLoopFunctionE(x)/12.0 
            - vertexOut.left()* vertexIn.left() * mF/mj * OneLoopFunctionF(x) * 2./3.
            - mi/mj * vertexOut.right() * vertexIn.left() * OneLoopFunctionE(x)/12.0;

         // eq. 18 & 18 of hep-ph/9510309 after replacement L <-> R (possibly with different sign)
         const std::complex<double> A1R = 
            - 1./18. * vertexOut.left() * vertexIn.right() * OneLoopFunctionD(x);
         const std::complex<double> A2R = 
            - vertexOut.right() * vertexIn.left() * OneLoopFunctionE(x)/12.0 
            - vertexOut.right()* vertexIn.right() * mF/mj * OneLoopFunctionF(x) * 2./3.
            - mi/mj * vertexOut.left() * vertexIn.right() * OneLoopFunctionE(x)/12.0; 

         const std::complex<double> massFactor = pow(mS,-2);

         res += oneOver32PiSqr * vector_boson_coupling * massFactor
            * std::valarray<std::complex<double>> {A1L, A1R, A2L, A2R};
      }
   }

   return res;
}

} // anonymous namespace
