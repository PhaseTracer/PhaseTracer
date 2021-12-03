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
 * @file ScalarSingletZ2DMMhInputMsInput_edm.cpp
 *
 * This file was generated with FlexibleSUSY 2.6.1 and SARAH 4.14.3 .
 */

#include "ScalarSingletZ2DMMhInputMsInput_edm.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_mass_eigenstates.hpp"

#include "cxx_qft/ScalarSingletZ2DMMhInputMsInput_qft.hpp"

#include "wrappers.hpp"
#include "numerics2.hpp"

#define INPUTPARAMETER(p) context.model.get_input().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

using namespace flexiblesusy;
using namespace ScalarSingletZ2DMMhInputMsInput_cxx_diagrams;

namespace {
static constexpr double oneOver16PiSquared = 0.0063325739776461107152;

/**
 * @class EDMVertexCorrectionSF
 * @brief A template that calculate contributions to the EDM
 *        of a given particle in a one loop diagram specified
 *        by a photon emitter and an exchange particle.
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
template<class EDMField, class PhotonEmitter, class ExchangeParticle>
struct EDMVertexCorrectionSF {
   static double value(const typename field_indices<EDMField>::type& indices,
                       context_base& context);
};

/**
* @class EDMVertexCorrectionFS
* @brief A template that calculate contributions to the EDM
*        of a given particle in a one loop diagram specified
*        by a photon emitter and an exchange particle.
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
template<class EDMField, class PhotonEmitter, class ExchangeParticle>
struct EDMVertexCorrectionFS {
   static double value(const typename field_indices<EDMField>::type& indices,
                       context_base& context);
};
} // anonymous namespace

namespace flexiblesusy {
namespace ScalarSingletZ2DMMhInputMsInput_edm {

}
} // namespace flexiblesusy

namespace {
/**
* @defgroup LoopFunctions Loop functions
* @brief The loop functions necessary for edm one-loop calculations.
*
* These are OneLoopFunctionA(), OneLoopFunctionB()
* as specified in arXiv:0808.1819
*/

double OneLoopFunctionA(double r)
{
   if (is_zero(r))
      return std::numeric_limits<double>::infinity();

   const double d = r - 1.0;

   if (std::abs(d) < 0.15) {
      return (-0.33333333333333333333 +
              0.25000000000000000000  * d -
              0.20000000000000000000  * d * d +
              0.16666666666666666667  * d * d * d -
              0.14285714285714285714  * d * d * d * d +
              0.12500000000000000000  * d * d * d * d * d -
              0.11111111111111111111  * d * d * d * d * d * d +
              0.10000000000000000000  * d * d * d * d * d * d * d -
              0.090909090909090909091 * d * d * d * d * d * d * d * d +
              0.083333333333333333333 * d * d * d * d * d * d * d * d * d -
              0.076923076923076923077 * d * d * d * d * d * d * d * d * d * d);
   }

   return 1.0 / (2.0 * d * d) * (3.0 - r - 2.0 * std::log(r) / d);
}

double OneLoopFunctionB(double r)
{
   if (is_zero(r))
      return 1.0/2.0;

   const double d = r - 1.0;

   if (std::abs(d) < 0.15)
      return (0.16666666666666666667 -
              0.083333333333333333333  * d +
              0.050000000000000000000  * d * d -
              0.033333333333333333333  * d * d * d +
              0.023809523809523809524  * d * d * d * d -
              0.017857142857142857143  * d * d * d * d * d +
              0.013888888888888888889  * d * d * d * d * d * d -
              0.011111111111111111111  * d * d * d * d * d * d * d +
              0.0090909090909090909091 * d * d * d * d * d * d * d * d -
              0.0075757575757575757576 * d * d * d * d * d * d * d * d * d +
              0.0064102564102564102564 * d * d * d * d * d * d * d * d * d * d);

   return 1.0 / (2.0 * d * d) * (1.0 + r - 2.0 * r * std::log(r) / d);
}

template<class EDMField, class PhotonEmitter, class ExchangeField>
double EDMVertexCorrectionFS<
EDMField, PhotonEmitter, ExchangeField
>::value(const typename field_indices<EDMField>::type& indices, context_base& context)
{
   double res = 0.0;

   using FermionVertex = Vertex<
                         EDMField,
                         ExchangeField,
                         PhotonEmitter
                         >;

   for (const auto& index: index_range<FermionVertex>()) {
      const auto edmFieldIndices = FermionVertex::template indices_of_field<0>(index);

      if (edmFieldIndices != indices)
         continue;

      const auto photonEmitterIndices = FermionVertex::template indices_of_field<2>(index);
      const auto exchangeFieldIndices = FermionVertex::template indices_of_field<1>(index);
      const auto vertex = FermionVertex::evaluate(index, context);

      const auto photonEmitterMass = context.mass<PhotonEmitter>(photonEmitterIndices);
      const auto exchangeFieldMass = context.mass<ExchangeField>(exchangeFieldIndices);

      const double photonEmitterCharge =
         - PhotonEmitter::electric_charge * unit_charge(context);

      constexpr double numericFactor = oneOver16PiSquared;
      const double massFactor = photonEmitterMass/(exchangeFieldMass * exchangeFieldMass);
      const double couplingFactor = (std::conj(vertex.right()) * vertex.left()).imag();

      const double massRatioSquared = Sqr(photonEmitterMass/exchangeFieldMass);
      const double loopFactor = photonEmitterCharge * OneLoopFunctionA(massRatioSquared);
      const double contribution = numericFactor * massFactor * couplingFactor * loopFactor;

      res += contribution;
   }

   return res;
}

template<class EDMField, class PhotonEmitter, class ExchangeField>
double EDMVertexCorrectionSF<
EDMField, PhotonEmitter, ExchangeField
>::value(const typename field_indices<EDMField>::type& indices, context_base& context)
{
   double res = 0.0;

   using FermionVertex = Vertex<
                         EDMField,
                         ExchangeField,
                         PhotonEmitter
                         >;

   for (const auto& index: index_range<FermionVertex>()) {
      const auto edmFieldIndices = FermionVertex::template indices_of_field<0>(index);

      if (edmFieldIndices != indices)
         continue;

      const auto photonEmitterIndices = FermionVertex::template indices_of_field<2>(index);
      const auto exchangeFieldIndices = FermionVertex::template indices_of_field<1>(index);
      const auto vertex = FermionVertex::evaluate(index, context);

      const auto photonEmitterMass = context.mass<PhotonEmitter>(photonEmitterIndices);
      const auto exchangeFieldMass = context.mass<ExchangeField>(exchangeFieldIndices);

      const double photonEmitterCharge =
         - PhotonEmitter::electric_charge * unit_charge(context);

      constexpr double numericFactor = oneOver16PiSquared;
      const double massFactor = exchangeFieldMass/(photonEmitterMass * photonEmitterMass);
      const double couplingFactor = (std::conj(vertex.right()) * vertex.left()).imag();

      const double massRatioSquared = Sqr(exchangeFieldMass / photonEmitterMass);
      const double loopFactor = photonEmitterCharge * OneLoopFunctionB(massRatioSquared);
      const double contribution = numericFactor * massFactor * couplingFactor * loopFactor;

      res += contribution;
   }

   return res;
}

} // anonymous namespace
