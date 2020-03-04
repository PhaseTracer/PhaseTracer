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

// File generated at Thu 7 Nov 2019 18:53:58

/**
 * @file cxx_qft/THDMIISNMSSMBCsimple_context_base.hpp
 *
 * This file was generated at Thu 7 Nov 2019 18:53:58 with FlexibleSUSY
 * 2.4.1 and SARAH 4.14.2 .
 */

#ifndef THDMIISNMSSMBCsimple_CXXQFT_CONTEXT_BASE_H
#define THDMIISNMSSMBCsimple_CXXQFT_CONTEXT_BASE_H

#include "THDMIISNMSSMBCsimple_mass_eigenstates.hpp"

#include "THDMIISNMSSMBCsimple_fields.hpp"
#include "THDMIISNMSSMBCsimple_mass_eigenstates.hpp"

namespace flexiblesusy {
namespace THDMIISNMSSMBCsimple_cxx_diagrams {

   struct context_base {
      THDMIISNMSSMBCsimple_mass_eigenstates model; ///< The model object.

      template <class Field>
      double mass(const typename field_indices<Field>::type& indices) const
      {
         using CleanField =
            typename fields::remove_lorentz_conjugation<Field>::type;
         return mass_impl<CleanField>(indices);
      }

      context_base(const THDMIISNMSSMBCsimple_mass_eigenstates& m) : model(m) {}
      context_base(const context_base&) = default;
      context_base(context_base&&) = default;

      context_base& operator=(const context_base&) = default;
      context_base& operator=(context_base&&) = default;

      virtual ~context_base(void) = default;

   private:
      template <class Field>
      double
      mass_impl(const typename field_indices<Field>::type& indices) const;
   };

   template<> inline
double context_base::mass_impl<fields::VG>(const std::array<int, 0>& indices) const
{ return model.get_MVG(); }

template<> inline
double context_base::mass_impl<fields::gG>(const std::array<int, 0>& indices) const
{ return model.get_MVG(); }

template<> inline
double context_base::mass_impl<fields::Fv>(const std::array<int, 1>& indices) const
{ return model.get_MFv(indices[0]); }

template<> inline
double context_base::mass_impl<fields::VP>(const std::array<int, 0>& indices) const
{ return model.get_MVP(); }

template<> inline
double context_base::mass_impl<fields::VZ>(const std::array<int, 0>& indices) const
{ return model.get_MVZ(); }

template<> inline
double context_base::mass_impl<fields::gP>(const std::array<int, 0>& indices) const
{ return model.get_MVP(); }

template<> inline
double context_base::mass_impl<fields::gZ>(const std::array<int, 0>& indices) const
{ return model.get_MVZ(); }

template<> inline
double context_base::mass_impl<fields::gWm>(const std::array<int, 0>& indices) const
{ return model.get_MVWm(); }

template<> inline
double context_base::mass_impl<fields::gWmC>(const std::array<int, 0>& indices) const
{ return model.get_MVWm(); }

template<> inline
double context_base::mass_impl<fields::hh>(const std::array<int, 1>& indices) const
{ return model.get_Mhh(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Ah>(const std::array<int, 1>& indices) const
{ return model.get_MAh(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Hm>(const std::array<int, 1>& indices) const
{ return model.get_MHm(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Fd>(const std::array<int, 1>& indices) const
{ return model.get_MFd(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Fu>(const std::array<int, 1>& indices) const
{ return model.get_MFu(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Fe>(const std::array<int, 1>& indices) const
{ return model.get_MFe(indices[0]); }

template<> inline
double context_base::mass_impl<fields::VWm>(const std::array<int, 0>& indices) const
{ return model.get_MVWm(); }

} // namespace THDMIISNMSSMBCsimple_cxx_diagrams
} // namespace flexiblesusy

#include "THDMIISNMSSMBCsimple_vertices.hpp"

#endif
