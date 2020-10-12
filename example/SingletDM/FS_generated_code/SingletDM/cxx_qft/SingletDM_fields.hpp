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

// File generated at Tue 18 Aug 2020 22:47:58

/**
 * @file cxx_qft/SingletDM_fields.hpp
 *
 * This file was generated at Tue 18 Aug 2020 22:47:58 with FlexibleSUSY
 * 2.4.2 and SARAH 4.14.3 .
 */

#ifndef SingletDM_CXXQFT_FIELDS_H
#define SingletDM_CXXQFT_FIELDS_H

#include <array>

#include <boost/array.hpp>
#include <boost/version.hpp>

#include <boost/mpl/erase.hpp>
#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>

#include <boost/fusion/adapted/boost_array.hpp>
#include <boost/fusion/adapted/mpl.hpp>

#if BOOST_VERSION >= 105800
#include <boost/fusion/include/move.hpp>
#else
#include <boost/fusion/include/copy.hpp>
#endif

namespace flexiblesusy
{
namespace SingletDM_cxx_diagrams
{
   // Declare a type that can hold the field indices for any given field
   template <class Field>
   struct field_indices {
      using type = std::array<int, Field::numberOfFieldIndices>;
   };

   namespace detail
   {
   template <class Field>
   struct number_of_field_indices {
      static constexpr int value =
         std::tuple_size<typename field_indices<Field>::type>::value;
      using type = boost::mpl::int_<value>;
   };

   template <class FieldSequence>
   struct total_number_of_field_indices {
      using type = typename boost::mpl::fold<
         FieldSequence, boost::mpl::int_<0>,
         boost::mpl::plus<boost::mpl::_1,
                          number_of_field_indices<boost::mpl::_2>>>::type;
      static constexpr int value = type::value;
   };
   } // namespace detail

   template <class Field>
   typename std::enable_if<Field::numberOfGenerations != 1, bool>::type
   isSMField(const typename field_indices<Field>::type& indices)
   {
      boost::array<bool, Field::numberOfGenerations> sm_flags;

#if BOOST_VERSION >= 105800
      boost::fusion::move(typename Field::sm_flags(), sm_flags);
#else
      boost::fusion::copy(typename Field::sm_flags(), sm_flags);
#endif

      return sm_flags[indices.front()];
   }

   template <class Field>
   typename std::enable_if<Field::numberOfGenerations == 1, bool>::type
   isSMField(const typename field_indices<Field>::type&)
   {
      return boost::mpl::at_c<typename Field::sm_flags, 0>::type::value;
   }

   namespace fields
   {
   template <class Field>
   struct bar {
      using index_bounds = typename Field::index_bounds;
      static constexpr int numberOfGenerations = Field::numberOfGenerations;
      using sm_flags = typename Field::sm_flags;
      static constexpr int numberOfFieldIndices = Field::numberOfFieldIndices;
      static constexpr double electric_charge = -Field::electric_charge;
      using lorentz_conjugate = Field;

      using type = bar<Field>;
   };

   template <class Field>
   struct conj {
      using index_bounds = typename Field::index_bounds;
      static constexpr int numberOfGenerations = Field::numberOfGenerations;
      using sm_flags = typename Field::sm_flags;
      static constexpr int numberOfFieldIndices = Field::numberOfFieldIndices;
      static constexpr double electric_charge = -Field::electric_charge;
      using lorentz_conjugate = Field;

      using type = conj<Field>;
   };

   // Double Lorentz conjugation
   template <class Field>
   struct bar<bar<Field>> {
      using type = Field;
   };
   template <class Field>
   struct conj<conj<Field>> {
      using type = Field;
   };

   // Remove Lorentz conjugation
   template <class Field>
   struct remove_lorentz_conjugation {
      using type = Field;
   };

   template <class Field>
   struct remove_lorentz_conjugation<bar<Field>> {
      using type = Field;
   };

   template <class Field>
   struct remove_lorentz_conjugation<conj<Field>> {
      using type = Field;
   };

   struct VG {
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electric_charge = 0;
   using lorentz_conjugate = VG;
};

struct gG {
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electric_charge = 0;
   using lorentz_conjugate = typename bar<gG>::type;
};

struct Hp {
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electric_charge = -1;
   using lorentz_conjugate = typename conj<Hp>::type;
};

struct ss {
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, false>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electric_charge = 0;
   using lorentz_conjugate = ss;
};

struct Fv {
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 3>
   >;
   static constexpr int numberOfGenerations = 3;
   using sm_flags = boost::mpl::vector_c<bool, true, true, true>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electric_charge = 0;
   using lorentz_conjugate = typename bar<Fv>::type;
};

struct Ah {
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electric_charge = 0;
   using lorentz_conjugate = Ah;
};

struct hh {
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electric_charge = 0;
   using lorentz_conjugate = hh;
};

struct VP {
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electric_charge = 0;
   using lorentz_conjugate = VP;
};

struct VZ {
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electric_charge = 0;
   using lorentz_conjugate = VZ;
};

struct gP {
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electric_charge = 0;
   using lorentz_conjugate = typename bar<gP>::type;
};

struct gZ {
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electric_charge = 0;
   using lorentz_conjugate = typename bar<gZ>::type;
};

struct gWp {
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electric_charge = 1;
   using lorentz_conjugate = typename bar<gWp>::type;
};

struct gWpC {
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electric_charge = -1;
   using lorentz_conjugate = typename bar<gWpC>::type;
};

struct Fd {
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 3>
   >;
   static constexpr int numberOfGenerations = 3;
   using sm_flags = boost::mpl::vector_c<bool, true, true, true>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electric_charge = -0.3333333333333333;
   using lorentz_conjugate = typename bar<Fd>::type;
};

struct Fu {
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 3>
   >;
   static constexpr int numberOfGenerations = 3;
   using sm_flags = boost::mpl::vector_c<bool, true, true, true>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electric_charge = 0.6666666666666666;
   using lorentz_conjugate = typename bar<Fu>::type;
};

struct Fe {
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 3>
   >;
   static constexpr int numberOfGenerations = 3;
   using sm_flags = boost::mpl::vector_c<bool, true, true, true>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electric_charge = -1;
   using lorentz_conjugate = typename bar<Fe>::type;
};

struct VWp {
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electric_charge = 1;
   using lorentz_conjugate = typename conj<VWp>::type;
};

// Named fields
using Photon = VP;
using Electron = Fe;

// Fields that are their own Lorentz conjugates.
template<> struct conj<VG> { using type = VG; };
template<> struct conj<ss> { using type = ss; };
template<> struct conj<Ah> { using type = Ah; };
template<> struct conj<hh> { using type = hh; };
template<> struct conj<VP> { using type = VP; };
template<> struct conj<VZ> { using type = VZ; };

using scalars = boost::mpl::vector<Hp, ss, Ah, hh>;
using fermions = boost::mpl::vector<Fv, Fd, Fu, Fe>;
using vectors = boost::mpl::vector<VG, VP, VZ, VWp>;
using ghosts = boost::mpl::vector<gG, gP, gZ, gWp, gWpC>;
   } // namespace fields

   using fields::bar;
   using fields::conj;
   using fields::remove_lorentz_conjugation;
} // namespace SingletDM_cxx_diagrams
} // namespace flexiblesusy

#endif
