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
 * @file cxx_qft/ScalarSingletZ2DMMhInputMsInput_fields.hpp
 *
 * This file was generated with FlexibleSUSY 2.6.1 and SARAH 4.14.3 .
 */

#ifndef ScalarSingletZ2DMMhInputMsInput_CXXQFT_FIELDS_H
#define ScalarSingletZ2DMMhInputMsInput_CXXQFT_FIELDS_H

#include <array> /**< Only for field_indices and isSMField*/

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

namespace flexiblesusy {
namespace ScalarSingletZ2DMMhInputMsInput_cxx_diagrams {
   /** @brief Declare a type that can hold the field indices for any given field.
    * @todo Field should have nested numberOfFieldIndices.
    * @returns `::type` evaluates to `std::array<int,int N>`.
    */
   template <class Field>
   struct field_indices {
      /** @b using keyword here is type alias
       * (see https://en.cppreference.com/w/cpp/language/type_alias).
       * It allows to use @code{.cpp}typename field_indices<Field>::type@endcode
       * to get @code{.cpp}std::array<int, Field::numberOfFieldIndices>@endcode
       * type in the code.
       */
      using type = std::array<int, Field::numberOfFieldIndices>;
   };

   namespace detail {
      /** @brief Declare a metafunction that gives the int value for Field
       * during compilation time or `int_<int N>` typename of this index.
       * @returns `::value` gives compile time integer value.
       * @returns `::type` gives `int_<int N>` typename.
       */
      template <class Field>
      struct number_of_field_indices {
         static constexpr int value =
            std::tuple_size<typename field_indices<Field>::type>::value;
         using type = boost::mpl::int_<value>;
      };
      /*
       * @param FieldSequence a Forward Sequence, e.g vector<type1,typ2,...>
       * (see https://www.boost.org/doc/libs/1_69_0/libs/mpl/doc/refmanual/forward-sequence.html).
       * @returns `::type` which accomodates `mpl::int_<N>` with `int N` being the
       * sum of @b number of indices of fields inside FieldSequence.
       * @returns `::value` which represents `N` where `N` being the `int`
       * sum of @b number of indices of fields inside FieldSequence
       */
      template <class FieldSequence>
      struct total_number_of_field_indices {
         using type = typename boost::mpl::fold<
            FieldSequence,
            boost::mpl::int_<0>,
            boost::mpl::plus<boost::mpl::_1,
                             number_of_field_indices<boost::mpl::_2>
                            >
            >::type;
         static constexpr int value = type::value;
      };
   } // namespace detail

   /* @note Is defined by compiler only if the number of generation for a field
    * is not 1.
    */
   template <class Field>
   std::enable_if_t<Field::numberOfGenerations != 1, bool>
   isSMField(const typename field_indices<Field>::type& indices) {
      boost::array<bool, Field::numberOfGenerations> sm_flags;

      #if BOOST_VERSION >= 105800
      boost::fusion::move(typename Field::sm_flags(), sm_flags); /**< Making copy of the information from Field to sm_flags */
      #else
      boost::fusion::copy(typename Field::sm_flags(), sm_flags);
      #endif

      return sm_flags[indices.front()]; /**< front() accesses the first element. Either all are SM indices or none. */
   }

   template <class Field>
   std::enable_if_t<Field::numberOfGenerations == 1, bool>
   isSMField(const typename field_indices<Field>::type&) {
      return boost::mpl::at_c<typename Field::sm_flags, 0>::type::value;
   }

   namespace fields
   {
   enum class ParticleType {
      scalar,
      fermion,
      vector,
      ghost
   };

   template<typename Field>
   struct is_massless {
      static constexpr bool value = Field::massless;
   };
   template<typename Field>
   constexpr bool is_massless_v = is_massless<Field>::value;

   enum class ParticleColorRep {
      singlet,
      triplet,
      anti_triplet,
      sextet,
      octet
   };

   template<typename Field>
   struct is_singlet {
      static constexpr bool value =
         Field::color_rep == ParticleColorRep::singlet;
   };
   template<typename Field>
   constexpr bool is_singlet_v = is_singlet<Field>::value;

   template<typename Field>
   struct is_triplet {
      static constexpr bool value = Field::color_rep == ParticleColorRep::triplet;
   };
   template<typename Field>
   constexpr bool is_triplet_v = is_triplet<Field>::value;

   template<typename Field>
   struct is_anti_triplet {
      static constexpr bool value =
         Field::color_rep == ParticleColorRep::anti_triplet;
   };
   template<typename Field>
   constexpr bool is_anti_triplet_v = is_anti_triplet<Field>::value;

   template<typename Field>
   struct is_octet {
      static constexpr bool value = Field::color_rep == ParticleColorRep::octet;
   };
   template<typename Field>
   constexpr bool is_octet_v = is_octet<Field>::value;

   template<typename Field>
   constexpr std::enable_if_t<
      is_triplet<Field>::value, ParticleColorRep
      >
   color_conj() {
      return ParticleColorRep::anti_triplet;
   }
   template<typename Field>
   constexpr std::enable_if_t<
      is_anti_triplet<Field>::value, ParticleColorRep
      >
   color_conj() {
      return ParticleColorRep::triplet;
   }
   template<typename Field>
   constexpr std::enable_if_t<
      !is_triplet<Field>::value && !is_anti_triplet<Field>::value, ParticleColorRep
      >
   color_conj() {
      return Field::color_rep;
   }

   template<class Field>
   struct bar {
      using index_bounds = typename Field::index_bounds;
      using sm_flags = typename Field::sm_flags;
      using lorentz_conjugate = Field;
      using type = bar<Field>;

      static constexpr int numberOfGenerations = Field::numberOfGenerations;
      static constexpr int numberOfFieldIndices = Field::numberOfFieldIndices;
      static constexpr double electric_charge = -Field::electric_charge;
      static constexpr auto particle_type = Field::particle_type;
      static constexpr auto color_rep = color_conj<Field>();
      static constexpr auto massless = Field::massless;
   };

   template<class Field>
   struct conj {
      using index_bounds = typename Field::index_bounds;
      using sm_flags = typename Field::sm_flags;
      using lorentz_conjugate = Field;
      using type = conj<Field>;

      static constexpr int numberOfGenerations = Field::numberOfGenerations;
      static constexpr int numberOfFieldIndices = Field::numberOfFieldIndices;
      static constexpr double electric_charge = -Field::electric_charge;
      static constexpr auto particle_type = Field::particle_type;
      static constexpr auto color_rep = color_conj<Field>();
      static constexpr auto massless = Field::massless;
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
   static constexpr auto particle_type = ParticleType::vector;
   static constexpr auto color_rep = ParticleColorRep::octet;
   static constexpr auto massless = true;
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
   static constexpr auto particle_type = ParticleType::ghost;
   static constexpr auto color_rep = ParticleColorRep::octet;
   static constexpr auto massless = true;
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
   static constexpr auto particle_type = ParticleType::scalar;
   static constexpr auto color_rep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
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
   static constexpr auto particle_type = ParticleType::scalar;
   static constexpr auto color_rep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
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
   static constexpr auto particle_type = ParticleType::fermion;
   static constexpr auto color_rep = ParticleColorRep::singlet;
   static constexpr auto massless = true;
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
   static constexpr auto particle_type = ParticleType::scalar;
   static constexpr auto color_rep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
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
   static constexpr auto particle_type = ParticleType::scalar;
   static constexpr auto color_rep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
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
   static constexpr auto particle_type = ParticleType::vector;
   static constexpr auto color_rep = ParticleColorRep::singlet;
   static constexpr auto massless = true;
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
   static constexpr auto particle_type = ParticleType::vector;
   static constexpr auto color_rep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
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
   static constexpr auto particle_type = ParticleType::ghost;
   static constexpr auto color_rep = ParticleColorRep::singlet;
   static constexpr auto massless = true;
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
   static constexpr auto particle_type = ParticleType::ghost;
   static constexpr auto color_rep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
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
   static constexpr auto particle_type = ParticleType::ghost;
   static constexpr auto color_rep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
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
   static constexpr auto particle_type = ParticleType::ghost;
   static constexpr auto color_rep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
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
   static constexpr auto particle_type = ParticleType::fermion;
   static constexpr auto color_rep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
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
   static constexpr auto particle_type = ParticleType::fermion;
   static constexpr auto color_rep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
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
   static constexpr auto particle_type = ParticleType::fermion;
   static constexpr auto color_rep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
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
   static constexpr auto particle_type = ParticleType::vector;
   static constexpr auto color_rep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
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

   template <class Field>
   struct is_scalar : public std::false_type {};

   template <class Field>
   struct is_scalar<bar<Field> > : public is_scalar<Field> {};

   template <class Field>
   struct is_scalar<conj<Field> > : public is_scalar<Field> {};

   template <class Field>
   struct is_fermion : public std::false_type {};

   template <class Field>
   struct is_fermion<bar<Field> > : public is_fermion<Field> {};

   template <class Field>
   struct is_fermion<conj<Field> > : public is_fermion<Field> {};

   template <class Field>
   struct is_vector : public std::false_type {};

   template <class Field>
   struct is_vector<bar<Field> > : public is_vector<Field> {};

   template <class Field>
   struct is_vector<conj<Field> > : public is_vector<Field> {};

   template <class Field>
   struct is_ghost : public std::false_type {};

   template <class Field>
   struct is_ghost<bar<Field> > : public is_ghost<Field> {};

   template <class Field>
   struct is_ghost<conj<Field> > : public is_ghost<Field> {};

template<>
struct is_vector<ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields::VG > : public std::true_type {};
template<>
struct is_ghost<ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields::gG > : public std::true_type {};
template<>
struct is_scalar<ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields::Hp > : public std::true_type {};
template<>
struct is_scalar<ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields::ss > : public std::true_type {};
template<>
struct is_fermion<ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields::Fv > : public std::true_type {};
template<>
struct is_scalar<ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields::Ah > : public std::true_type {};
template<>
struct is_scalar<ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields::hh > : public std::true_type {};
template<>
struct is_vector<ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields::VP > : public std::true_type {};
template<>
struct is_vector<ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields::VZ > : public std::true_type {};
template<>
struct is_ghost<ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields::gP > : public std::true_type {};
template<>
struct is_ghost<ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields::gZ > : public std::true_type {};
template<>
struct is_ghost<ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields::gWp > : public std::true_type {};
template<>
struct is_ghost<ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields::gWpC > : public std::true_type {};
template<>
struct is_fermion<ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields::Fd > : public std::true_type {};
template<>
struct is_fermion<ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields::Fu > : public std::true_type {};
template<>
struct is_fermion<ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields::Fe > : public std::true_type {};
template<>
struct is_vector<ScalarSingletZ2DMMhInputMsInput_cxx_diagrams::fields::VWp > : public std::true_type {};


   } // namespace fields

   using fields::bar;
   using fields::conj;
   using fields::remove_lorentz_conjugation;
} // namespace ScalarSingletZ2DMMhInputMsInput_cxx_diagrams

} // namespace flexiblesusy

#endif
