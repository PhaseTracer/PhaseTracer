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

// File generated at Tue 17 Nov 2020 16:11:26

#ifndef ScalarSingletZ2DM_UTILITIES_H
#define ScalarSingletZ2DM_UTILITIES_H

#include "ScalarSingletZ2DM_mass_eigenstates.hpp"
#include "ScalarSingletZ2DM_info.hpp"
#include "wrappers.hpp"

#include <Eigen/Core>
#include <array>
#include <string>
#include <vector>
#include <valarray>
#include <utility>

namespace softsusy {
class QedQcd;
}

namespace flexiblesusy {

struct ScalarSingletZ2DM_observables;
class Physical_input;

class ScalarSingletZ2DM_parameter_getter {
private:
   std::vector<std::string> get_mass_names(const std::string& head = "") const {
      using namespace ScalarSingletZ2DM_info;
      std::vector<std::string> masses;
      for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
         for (int m = 0; m < particle_multiplicities[i]; m++) {
            masses.push_back(
               head + "M" + particle_names[i] +
               (particle_multiplicities[i] == 1 ? "" : "("
                + std::to_string(static_cast<long long>(m)) + ")"));
         }
      }

      masses.shrink_to_fit();

      return masses;
   }

   std::array<std::string, ScalarSingletZ2DM_info::NUMBER_OF_MIXINGS> get_mixing_names() const {
      return ScalarSingletZ2DM_info::particle_mixing_names;
   }

public:
   /// returns DR-bar parameters
   Eigen::ArrayXd get_parameters(const ScalarSingletZ2DM_mass_eigenstates& model) {
      return model.get();
   }

   /// returns names of DR-bar parameters
   std::array<std::string, ScalarSingletZ2DM_info::NUMBER_OF_PARAMETERS> get_parameter_names() const {
      return ScalarSingletZ2DM_info::parameter_names;
   }

   /// returns names of particles
   std::array<std::string, ScalarSingletZ2DM_info::NUMBER_OF_PARTICLES> get_particle_names() const {
      return ScalarSingletZ2DM_info::particle_names;
   }

   /// returns names of DR-bar masses
   std::vector<std::string> get_DRbar_mass_names() const { return get_mass_names(); }

   /// returns names of pole masses
   std::vector<std::string> get_pole_mass_names() const { return get_mass_names("Pole"); }

   /// returns names of DR-bar mixing matrices
   std::array<std::string, ScalarSingletZ2DM_info::NUMBER_OF_MIXINGS> get_DRbar_mixing_names() const {
      return get_mixing_names();
   }

   /// returns names of pole mixing matrices
   std::array<std::string, ScalarSingletZ2DM_info::NUMBER_OF_MIXINGS> get_pole_mixing_names() const {
      auto mixing_names = get_mixing_names();
      for (auto& n: mixing_names)
         n = std::string("Pole") + n;
      return mixing_names;
   }

   /// returns names of input parameters
   std::array<std::string, ScalarSingletZ2DM_info::NUMBER_OF_INPUT_PARAMETERS> get_input_parameter_names() const {
      return ScalarSingletZ2DM_info::input_parameter_names;
   }

   /// returns names of input parameters
   std::array<std::string, ScalarSingletZ2DM_info::NUMBER_OF_EXTRA_PARAMETERS> get_extra_parameter_names() const {
      return ScalarSingletZ2DM_info::extra_parameter_names;
   }

   /// returns number of particle masses
   decltype(ScalarSingletZ2DM_info::NUMBER_OF_MASSES) get_number_of_masses() const {
      return ScalarSingletZ2DM_info::NUMBER_OF_MASSES;
   }
};

class ScalarSingletZ2DM_spectrum_plotter {
public:
   ScalarSingletZ2DM_spectrum_plotter() = default;
   explicit ScalarSingletZ2DM_spectrum_plotter(const ScalarSingletZ2DM_mass_eigenstates&);
   void extract_spectrum(const ScalarSingletZ2DM_mass_eigenstates&);
   void write_to_file(const std::string&) const;

private:
   struct TParticle {
      std::string name;
      std::string latex_name;
      std::valarray<double> masses;
      TParticle(const std::string& name_, const std::string& latex_name_,
                const std::valarray<double>& masses_)
         : name(name_)
         , latex_name(latex_name_)
         , masses(masses_)
         {}
   };
   using TSpectrum = std::vector<TParticle>;
   TSpectrum spectrum{};
   double scale{0.};
   int width{16};

   void write_spectrum(const TSpectrum&, std::ofstream&) const;
};

namespace ScalarSingletZ2DM_database {

/// append parameter point to database
void to_database(
   const std::string&,
   const ScalarSingletZ2DM_mass_eigenstates&,
   const softsusy::QedQcd* qedqcd = nullptr,
   const Physical_input* physical_input = nullptr,
   const ScalarSingletZ2DM_observables* observables = nullptr);

/// fill model from an entry of the database
ScalarSingletZ2DM_mass_eigenstates from_database(
   const std::string&,
   long long,
   softsusy::QedQcd* qedqcd = nullptr,
   Physical_input* physical_input = nullptr,
   ScalarSingletZ2DM_observables* observables = nullptr);

} // namespace ScalarSingletZ2DM_database

} // namespace flexiblesusy

#endif
