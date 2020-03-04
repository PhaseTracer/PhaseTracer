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

// File generated at Thu 7 Nov 2019 18:53:46

#ifndef THDMIISNMSSMBCsimple_PHYSICAL_H
#define THDMIISNMSSMBCsimple_PHYSICAL_H

#include <Eigen/Core>

#include <iosfwd>

namespace flexiblesusy {

struct THDMIISNMSSMBCsimple_physical {
   void clear();
   void convert_to_hk();   ///< converts pole masses to HK convention
   void convert_to_slha(); ///< converts pole masses to SLHA convention
   Eigen::ArrayXd get() const; ///< returns array with all masses and mixings
   void set(const Eigen::ArrayXd&); ///< set all masses and mixings
   Eigen::ArrayXd get_masses() const; ///< returns array with all masses
   void set_masses(const Eigen::ArrayXd&); ///< set all masses
   void print(std::ostream&) const;

   double MVG{};
   Eigen::Array<double,3,1> MFv{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> Mhh{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MAh{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,2,1> MHm{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,3,1> MFd{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFu{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFe{Eigen::Array<double,3,1>::Zero()};
   double MVWm{};
   double MVP{};
   double MVZ{};

   Eigen::Matrix<double,3,3> ZH{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> ZA{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,2,2> ZP{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Vd{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ud{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Vu{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Uu{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ve{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ue{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<double,2,2> ZZ{Eigen::Matrix<double,2,2>::Zero()};

};

std::ostream& operator<<(std::ostream&, const THDMIISNMSSMBCsimple_physical&);

} // namespace flexiblesusy

#endif
