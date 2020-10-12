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

// File generated at Tue 18 Aug 2020 22:47:51

#ifndef SingletDM_INFO_H
#define SingletDM_INFO_H

#include "problems.hpp"

#include <array>
#include <iosfwd>
#include <string>

namespace flexiblesusy {

namespace SingletDM_info {
   enum Particles : int { VG, Hp, ss, Fv, Ah, hh, Fd, Fu, Fe, VWp, VP, VZ,
      NUMBER_OF_PARTICLES };

   enum Masses : int { MVG, MHp, Mss, MFv_1, MFv_2, MFv_3, MAh, Mhh, MFd_1, MFd_2,
      MFd_3, MFu_1, MFu_2, MFu_3, MFe_1, MFe_2, MFe_3, MVWp, MVP, MVZ,
      NUMBER_OF_MASSES };

   enum Parameters : int { g1, g2, g3, LamS, LamSH, LamH, Yu0_0, Yu0_1, Yu0_2,
      Yu1_0, Yu1_1, Yu1_2, Yu2_0, Yu2_1, Yu2_2, Yd0_0, Yd0_1, Yd0_2, Yd1_0, Yd1_1,
      Yd1_2, Yd2_0, Yd2_1, Yd2_2, Ye0_0, Ye0_1, Ye0_2, Ye1_0, Ye1_1, Ye1_2, Ye2_0,
      Ye2_1, Ye2_2, muS, muH, v, NUMBER_OF_PARAMETERS };

   enum Mixings : int { ReVd0_0, ImVd0_0, ReVd0_1, ImVd0_1, ReVd0_2, ImVd0_2,
      ReVd1_0, ImVd1_0, ReVd1_1, ImVd1_1, ReVd1_2, ImVd1_2, ReVd2_0, ImVd2_0,
      ReVd2_1, ImVd2_1, ReVd2_2, ImVd2_2, ReUd0_0, ImUd0_0, ReUd0_1, ImUd0_1,
      ReUd0_2, ImUd0_2, ReUd1_0, ImUd1_0, ReUd1_1, ImUd1_1, ReUd1_2, ImUd1_2,
      ReUd2_0, ImUd2_0, ReUd2_1, ImUd2_1, ReUd2_2, ImUd2_2, ReVu0_0, ImVu0_0,
      ReVu0_1, ImVu0_1, ReVu0_2, ImVu0_2, ReVu1_0, ImVu1_0, ReVu1_1, ImVu1_1,
      ReVu1_2, ImVu1_2, ReVu2_0, ImVu2_0, ReVu2_1, ImVu2_1, ReVu2_2, ImVu2_2,
      ReUu0_0, ImUu0_0, ReUu0_1, ImUu0_1, ReUu0_2, ImUu0_2, ReUu1_0, ImUu1_0,
      ReUu1_1, ImUu1_1, ReUu1_2, ImUu1_2, ReUu2_0, ImUu2_0, ReUu2_1, ImUu2_1,
      ReUu2_2, ImUu2_2, ReVe0_0, ImVe0_0, ReVe0_1, ImVe0_1, ReVe0_2, ImVe0_2,
      ReVe1_0, ImVe1_0, ReVe1_1, ImVe1_1, ReVe1_2, ImVe1_2, ReVe2_0, ImVe2_0,
      ReVe2_1, ImVe2_1, ReVe2_2, ImVe2_2, ReUe0_0, ImUe0_0, ReUe0_1, ImUe0_1,
      ReUe0_2, ImUe0_2, ReUe1_0, ImUe1_0, ReUe1_1, ImUe1_1, ReUe1_2, ImUe1_2,
      ReUe2_0, ImUe2_0, ReUe2_1, ImUe2_1, ReUe2_2, ImUe2_2, ZZ0_0, ZZ0_1, ZZ1_0,
      ZZ1_1, NUMBER_OF_MIXINGS };

   enum Input_parameters : int { HiggsIN, LamSHInput, LamSInput, muSInput, QEWSB,
      Qin, NUMBER_OF_INPUT_PARAMETERS };

   enum Extra_parameters : int { NUMBER_OF_EXTRA_PARAMETERS };

   extern const double normalization_g1;
   extern const double normalization_g2;
   extern const double normalization_g3;

   extern const std::array<int, NUMBER_OF_PARTICLES> particle_multiplicities;
   extern const std::array<std::string, NUMBER_OF_PARTICLES> particle_names;
   extern const std::array<std::string, NUMBER_OF_PARTICLES> particle_latex_names;
   extern const std::array<std::string, NUMBER_OF_PARAMETERS> parameter_names;
   extern const std::array<std::string, NUMBER_OF_MIXINGS> particle_mixing_names;
   extern const std::array<std::string, NUMBER_OF_INPUT_PARAMETERS> input_parameter_names;
   extern const std::array<std::string, NUMBER_OF_EXTRA_PARAMETERS> extra_parameter_names;
   extern const std::string model_name;
   constexpr bool is_low_energy_model = false;
   constexpr bool is_supersymmetric_model = false;
   constexpr bool is_FlexibleEFTHiggs = true;

   void print(std::ostream&);

   class SingletDM_particle_names : public Names {
   public:
      virtual ~SingletDM_particle_names() = default;
      virtual const std::string& get(int index) const override {
         return particle_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARTICLES;
      }
   };

   class SingletDM_parameter_names : public Names {
   public:
      virtual ~SingletDM_parameter_names() = default;
      virtual const std::string& get(int index) const override {
         return parameter_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARAMETERS;
      }
   };

   const SingletDM_particle_names  particle_names_getter{};
   const SingletDM_parameter_names parameter_names_getter{};

} // namespace SingletDM_info

} // namespace flexiblesusy

#endif
