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
 * @file cxx_qft/ScalarSingletZ2DMMhInputMsInput_npointfunctions_wilsoncoeffs.hpp
 *
 * This file was generated with FlexibleSUSY 2.5.0 and SARAH 4.14.3 .
 */

#ifndef ScalarSingletZ2DMMhInputMsInput_CXXQFT_NPOINTFUNCTIONS_H
#define ScalarSingletZ2DMMhInputMsInput_CXXQFT_NPOINTFUNCTIONS_H

#include "ScalarSingletZ2DMMhInputMsInput_qft.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_vertices.hpp"

#include <Eigen/Core>

#include <boost/mpl/map.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/apply.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/zip_view.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/unpack_args.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/front_inserter.hpp>

#include <boost/fusion/include/map.hpp>
#include <boost/fusion/adapted/mpl.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <boost/fusion/functional/adapter/fused.hpp>

#include <boost/functional/value_factory.hpp>

namespace flexiblesusy {
namespace ScalarSingletZ2DMMhInputMsInput_cxx_diagrams {
namespace npointfunctions {
namespace detail {
   struct metafunction_class_mpl_map {
      template<class ...Args>
      struct apply {
         using type = boost::mpl::map<Args...>;
      };
   };

   struct metafunction_class_fusion_map {
      template<class ...Args>
      struct apply {
         using type = boost::fusion::map<Args...>;
      };
   };

  template<class GenericKeys, template<typename> class F, class Wilsoncoeffs_Length, class StateVector
           >
   class accumulate_generic_impl {
      std::array<std::complex<double>, Wilsoncoeffs_Length::value> &value;
      StateVector state;

      public:
      template<class ...Args>
      accumulate_generic_impl(
         std::array<std::complex<double>, Wilsoncoeffs_Length::value> &v,
         Args && ...args ) :
      value( v ), state( std::forward<Args>( args )... ) {
      } // End of constructor

      template<class InsertionData>
      void operator()( InsertionData ) {
         using namespace boost::mpl;

         using field_insertions = typename at_c<InsertionData, 0>::type;
         static constexpr int combinatorial_factor = at_c<InsertionData, 1>::type::value;
         static constexpr auto colour_factor = at_c<InsertionData, 2>::type::value;

         using pairs = typename transform<GenericKeys,
                                          field_insertions,
                                          quote2<pair>,
                                          front_inserter<vector<>>
                                          >::type;

         using GenericFieldMap = typename apply<unpack_args<metafunction_class_mpl_map>,
                                                pairs
                                                >::type;

         const auto value_boost = boost::fusion::fused<boost::value_factory<F<GenericFieldMap>>
                                                       >{}( state )();
         const auto value_length = value.size();

         for ( std::size_t i=0; i<value_length; i++ ) {
            value.at(i) += colour_factor * combinatorial_factor * value_boost.at(i);
         }
      } // End of operator()
   }; // End of class accumulate_generic_impl

  template<class GenericFieldKeys,
           template<typename> class F,
           class Wilsoncoeffs_Length,
           class ...Args
           >
   accumulate_generic_impl<GenericFieldKeys, F, Wilsoncoeffs_Length,
      typename boost::fusion::result_of::make_vector<Args...>::type>
   make_accumulate_generic_impl(
      std::array<std::complex<double>,
      Wilsoncoeffs_Length::value> &v,
      Args &&... args
   ) {
      return accumulate_generic_impl<GenericFieldKeys, F, Wilsoncoeffs_Length,
         typename boost::fusion::result_of::make_vector<Args...>::type
      >{ v, std::forward<Args>( args )... };
   } // End of function make_accumulate_generic_impl

   template<class GenericFieldKey>
   struct field_vector {
      using type = typename GenericFieldKey::field_vector;
   };

   template<std::intmax_t Num, std::intmax_t Denom = 1>
   struct ratio_helper {
      static constexpr double value = static_cast<double>(Num) / static_cast<double>(Denom);
   };

   template<class RealPart, class ImPart>
   struct complex_helper {
      static constexpr std::complex<double> value {RealPart::value, ImPart::value};
   };
} // End of namespace detail

struct lorentz_scalar {}; // dummy struct
struct lorentz_left {}; // dummy struct
struct lorentz_right {}; // dummy struct
struct lorentz_momentum_diff {
   int minuend_index;
   int subtrahend_index;
}; // dummy struct
struct lorentz_inverse_metric {}; // dummy struct
struct triple_vector {}; // dummy struct

struct context_with_vertices : context_base {
   using context_base::context_base;

   template<class ...Fields>
   std::complex<double> vertex(lorentz_scalar,
      const typename Vertex<Fields...>::indices_type &indices) const {
      return Vertex<Fields...>::evaluate( indices, model ).value();
   }

   template<class ...Fields>
   std::complex<double> vertex(lorentz_left,
      const typename Vertex<Fields...>::indices_type &indices ) const {
      return Vertex<Fields...>::evaluate( indices, model ).left();
   }

   template<class ...Fields>
   std::complex<double> vertex(lorentz_right,
      const typename Vertex<Fields...>::indices_type &indices ) const {
      return Vertex<Fields...>::evaluate( indices, model ).right();
   }

   template<class ...Fields>
   std::complex<double> vertex(lorentz_momentum_diff lmd_indices,
      const typename Vertex<Fields...>::indices_type &indices ) const {
      return Vertex<Fields...>::evaluate( indices, model ).value(
      lmd_indices.minuend_index, lmd_indices.subtrahend_index );
   }

   template<class ...Fields>
   std::complex<double> vertex(lorentz_inverse_metric,
      const typename Vertex<Fields...>::indices_type &indices ) const {
      return Vertex<Fields...>::evaluate( indices, model ).value();
   }

   template<class ...Fields>
   std::complex<double> vertex(triple_vector,
      const typename Vertex<Fields...>::indices_type &indices ) const {
      // FIXME: hardcoded odd permutation to agree with SPheno
      return Vertex<Fields...>::evaluate( indices, model ).value( TripleVectorVertex::odd_permutation{} );
   }

   double scale( void ) const {
      return dynamic_cast<ScalarSingletZ2DMMhInputMsInput_mass_eigenstates const&>(model).get_scale();
   }
}; // End of context_with_vertices

template<int NumberOfExternalIndices, int NumberOfExternalMomenta>
class correlation_function_context :
public context_with_vertices {
private:
   const std::array<int,NumberOfExternalIndices> &ext_indices;
   const std::array<Eigen::Vector4d,NumberOfExternalMomenta> &ext_momenta;

public:
   correlation_function_context(
      const ScalarSingletZ2DMMhInputMsInput_mass_eigenstates &m,
      const std::array<int, NumberOfExternalIndices> &ei,
      const std::array<Eigen::Vector4d, NumberOfExternalMomenta> &em) :
   context_with_vertices( m ), ext_indices( ei ), ext_momenta( em ) {
   }

  using context_with_vertices::mass;
  using context_with_vertices::vertex;

  const std::array<int, NumberOfExternalIndices> &
  external_indices( void ) const { return ext_indices; }

  int external_indices( int index ) const { return ext_indices.at( index ); }

  const std::array<Eigen::Vector4d, NumberOfExternalMomenta> &
  external_momenta( void ) const { return ext_momenta; }

  Eigen::Vector4d external_momenta( int index ) const
  { return ext_momenta.at( index ); }
}; // End of correlation_function_context

template<class GenericFieldMap>
class field_index_map {
   using _1 = boost::mpl::_1;
   using pairs = typename boost::mpl::transform<
      GenericFieldMap,
      boost::fusion::pair< boost::mpl::first<_1>, field_indices<boost::mpl::second<_1>> >,
      boost::mpl::front_inserter<boost::mpl::vector<>>
   >::type;
public:
   using type = typename boost::mpl::apply<
      boost::mpl::unpack_args<detail::metafunction_class_fusion_map>,
      pairs
   >::type;
};

template<class GenericFieldMap>
class index_map_interface {
   using map = typename field_index_map<GenericFieldMap>::type;
   const map &index_map_;

   public:
   index_map_interface( const map &map ) : index_map_( map ) {}

   const map &index_map( void ) const { return index_map_; }
};

template<
   class GenericFieldKeys,
   class GenericFieldInsertions,
   class CombinatorialFactors,
   class ColourFactors,
   class Wilsoncoeffs_Length,
   template<typename> class F,
   class ...Args
>
std::array<std::complex<double>, Wilsoncoeffs_Length::value> accumulate_generic( Args &&...args ) {
   using insertion_data = boost::mpl::zip_view<
      boost::mpl::vector<
         GenericFieldInsertions,
         CombinatorialFactors,
         ColourFactors
      >
   >;

   std::array<std::complex<double>, Wilsoncoeffs_Length::value> value;
   value.fill(0.0);

   boost::mpl::for_each<insertion_data>(
      detail::make_accumulate_generic_impl<
         GenericFieldKeys,
         F,
         Wilsoncoeffs_Length
         >( value, std::forward<Args>( args )...
      )
   );

   return value;
} // End of accumulate_generic
}
}
}

#endif
