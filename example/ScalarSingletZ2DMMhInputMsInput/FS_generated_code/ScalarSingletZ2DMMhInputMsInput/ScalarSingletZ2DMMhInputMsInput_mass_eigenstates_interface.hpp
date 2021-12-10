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
 * @file ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_interface.hpp
 *
 * @brief Contains the mass eigenstates interface class definition
 *
 * This file was generated with FlexibleSUSY 2.6.1 and SARAH 4.14.3 .
 */

#ifndef ScalarSingletZ2DMMhInputMsInput_MASS_EIGENSTATES_INTERFACE_H
#define ScalarSingletZ2DMMhInputMsInput_MASS_EIGENSTATES_INTERFACE_H

#include "ScalarSingletZ2DMMhInputMsInput_physical.hpp"

#include <memory>

#include <Eigen/Core>

namespace flexiblesusy {

class Problems;
struct ScalarSingletZ2DMMhInputMsInput_input_parameters;

/**
 * @class ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_interface
 * @brief Interface definition for model parameters, masses and mixings
 */
class ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_interface {
public:
   virtual ~ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_interface() {}

   virtual std::unique_ptr<ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_interface> clone() const = 0;

   virtual void calculate_tree_level_mass_spectrum() = 0;
   virtual void calculate_pole_mass_spectrum() = 0;
   virtual void calculate_mass_spectrum() = 0;

   virtual int solve_ewsb_equations_tree_level() = 0;
   virtual int solve_ewsb_equations() = 0;

   virtual Eigen::ArrayXd get_tree_level_masses() const = 0;
   virtual Eigen::ArrayXd get_tree_level_masses_and_mixings() const = 0;
   virtual const ScalarSingletZ2DMMhInputMsInput_input_parameters& get_input_parameters() const = 0;
   virtual ScalarSingletZ2DMMhInputMsInput_input_parameters& get_input_parameters() = 0;
   virtual Eigen::ArrayXd get_extra_parameters() const = 0;
   virtual const ScalarSingletZ2DMMhInputMsInput_physical& get_physical() const = 0;
   virtual ScalarSingletZ2DMMhInputMsInput_physical& get_physical() = 0;
   virtual const Problems& get_problems() const = 0;
   virtual Problems& get_problems() = 0;
   virtual void set_tree_level_masses(const Eigen::ArrayXd&) = 0;
   virtual void set_tree_level_masses_and_mixings(const Eigen::ArrayXd&) = 0;
   virtual void set_extra_parameters(const Eigen::ArrayXd&) = 0;
   virtual void set_physical(const ScalarSingletZ2DMMhInputMsInput_physical&) = 0;
   virtual void clear_problems() = 0;

   virtual double get_g1() const = 0;
   virtual double get_g2() const = 0;
   virtual double get_g3() const = 0;
   virtual double get_LamS() const = 0;
   virtual double get_LamSH() const = 0;
   virtual double get_LamH() const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_Yu() const = 0;
   virtual double get_Yu(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_Yd() const = 0;
   virtual double get_Yd(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_Ye() const = 0;
   virtual double get_Ye(int i, int k) const = 0;
   virtual double get_muS2() const = 0;
   virtual double get_muH2() const = 0;
   virtual double get_v() const = 0;
   virtual void set_g1(double g1_) = 0;
   virtual void set_g2(double g2_) = 0;
   virtual void set_g3(double g3_) = 0;
   virtual void set_LamS(double LamS_) = 0;
   virtual void set_LamSH(double LamSH_) = 0;
   virtual void set_LamH(double LamH_) = 0;
   virtual void set_Yu(const Eigen::Matrix<double,3,3>& Yu_) = 0;
   virtual void set_Yu(int i, int k, const double& value) = 0;
   virtual void set_Yd(const Eigen::Matrix<double,3,3>& Yd_) = 0;
   virtual void set_Yd(int i, int k, const double& value) = 0;
   virtual void set_Ye(const Eigen::Matrix<double,3,3>& Ye_) = 0;
   virtual void set_Ye(int i, int k, const double& value) = 0;
   virtual void set_muS2(double muS2_) = 0;
   virtual void set_muH2(double muH2_) = 0;
   virtual void set_v(double v_) = 0;
   virtual double get_MVG() const = 0;
   virtual double get_MHp() const = 0;
   virtual double get_Mss() const = 0;
   virtual const Eigen::Array<double,3,1>& get_MFv() const = 0;
   virtual double get_MFv(int i) const = 0;
   virtual double get_MAh() const = 0;
   virtual double get_Mhh() const = 0;
   virtual const Eigen::Array<double,3,1>& get_MFd() const = 0;
   virtual double get_MFd(int i) const = 0;
   virtual const Eigen::Array<double,3,1>& get_MFu() const = 0;
   virtual double get_MFu(int i) const = 0;
   virtual const Eigen::Array<double,3,1>& get_MFe() const = 0;
   virtual double get_MFe(int i) const = 0;
   virtual double get_MVWp() const = 0;
   virtual double get_MVP() const = 0;
   virtual double get_MVZ() const = 0;

   virtual const Eigen::Matrix<std::complex<double>,3,3>& get_Vd() const = 0;
   virtual std::complex<double> get_Vd(int i, int k) const = 0;
   virtual const Eigen::Matrix<std::complex<double>,3,3>& get_Ud() const = 0;
   virtual std::complex<double> get_Ud(int i, int k) const = 0;
   virtual const Eigen::Matrix<std::complex<double>,3,3>& get_Vu() const = 0;
   virtual std::complex<double> get_Vu(int i, int k) const = 0;
   virtual const Eigen::Matrix<std::complex<double>,3,3>& get_Uu() const = 0;
   virtual std::complex<double> get_Uu(int i, int k) const = 0;
   virtual const Eigen::Matrix<std::complex<double>,3,3>& get_Ve() const = 0;
   virtual std::complex<double> get_Ve(int i, int k) const = 0;
   virtual const Eigen::Matrix<std::complex<double>,3,3>& get_Ue() const = 0;
   virtual std::complex<double> get_Ue(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,2,2>& get_ZZ() const = 0;
   virtual double get_ZZ(int i, int k) const = 0;



   virtual double get_mass_matrix_VG() const = 0;
   virtual void calculate_MVG() = 0;
   virtual double get_mass_matrix_Hp() const = 0;
   virtual void calculate_MHp() = 0;
   virtual double get_mass_matrix_ss() const = 0;
   virtual void calculate_Mss() = 0;
   virtual Eigen::Matrix<double,3,3> get_mass_matrix_Fv() const = 0;
   virtual void calculate_MFv() = 0;
   virtual double get_mass_matrix_Ah() const = 0;
   virtual void calculate_MAh() = 0;
   virtual double get_mass_matrix_hh() const = 0;
   virtual void calculate_Mhh() = 0;
   virtual Eigen::Matrix<double,3,3> get_mass_matrix_Fd() const = 0;
   virtual void calculate_MFd() = 0;
   virtual Eigen::Matrix<double,3,3> get_mass_matrix_Fu() const = 0;
   virtual void calculate_MFu() = 0;
   virtual Eigen::Matrix<double,3,3> get_mass_matrix_Fe() const = 0;
   virtual void calculate_MFe() = 0;
   virtual double get_mass_matrix_VWp() const = 0;
   virtual void calculate_MVWp() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_VPVZ() const = 0;
   virtual void calculate_MVPVZ() = 0;
   virtual double get_ewsb_eq_hh_1() const = 0;
   virtual double ThetaW() const = 0;
   virtual double VEV() const = 0;
};

} // namespace flexiblesusy

#endif
