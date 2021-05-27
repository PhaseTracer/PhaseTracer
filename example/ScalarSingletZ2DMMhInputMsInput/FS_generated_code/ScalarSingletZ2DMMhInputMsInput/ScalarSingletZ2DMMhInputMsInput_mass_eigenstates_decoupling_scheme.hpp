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
 * @file ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_decoupling_scheme.hpp
 *
 * @brief Defines model class for Stöckinger/Kotlarski decoupling scheme.
 *
 * This file was generated with FlexibleSUSY 2.5.0 and SARAH 4.14.3 .
 */

#ifndef ScalarSingletZ2DMMhInputMsInput_MASS_EIGENSTATES_DECOUPLING_SCHEME_H
#define ScalarSingletZ2DMMhInputMsInput_MASS_EIGENSTATES_DECOUPLING_SCHEME_H

#include "ScalarSingletZ2DMMhInputMsInput_info.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_physical.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_soft_parameters.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_interface.hpp"
#include "problems.hpp"

#include <iosfwd>
#include <memory>

#include <Eigen/Core>

#define SUPER(p) ScalarSingletZ2DMMhInputMsInput_soft_parameters::p

namespace flexiblesusy {

namespace standard_model {
class Standard_model;
}

struct ScalarSingletZ2DMMhInputMsInput_input_parameters;
class ScalarSingletZ2DMMhInputMsInput_mass_eigenstates;

/**
 * @class ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_decoupling_scheme
 *
 * @brief model class with routines for determing masses and mixings
 * and EWSB in the Stöckinger/Kotlarski decoupling scheme
 */
class ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_decoupling_scheme
   : private ScalarSingletZ2DMMhInputMsInput_soft_parameters
   , public ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_interface
{
public:
   explicit ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_decoupling_scheme(const ScalarSingletZ2DMMhInputMsInput_input_parameters& input_ = ScalarSingletZ2DMMhInputMsInput_input_parameters());
   explicit ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_decoupling_scheme(const ScalarSingletZ2DMMhInputMsInput_mass_eigenstates&);
   ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_decoupling_scheme(const ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_decoupling_scheme&) = default;
   ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_decoupling_scheme(ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_decoupling_scheme&&) = default;
   virtual ~ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_decoupling_scheme() = default;
   ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_decoupling_scheme& operator=(const ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_decoupling_scheme&) = default;
   ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_decoupling_scheme& operator=(ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_decoupling_scheme&&) = default;
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW

   std::unique_ptr<ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_interface> clone() const override;

   /// number of EWSB equations
   static const int number_of_ewsb_equations = 1;

   void check_pole_masses_for_tachyons();
   void do_force_output(bool);
   bool do_force_output() const;
   void fill_from(const standard_model::Standard_model&);
   void fill_from(const ScalarSingletZ2DMMhInputMsInput_mass_eigenstates&);
   void reorder_tree_level_masses();
   void reorder_pole_masses();
   void print(std::ostream&) const override;
   void set_precision(double);
   double get_precision() const;
   void clear() override;

   // mass_eigenstates_interface functions

   void calculate_tree_level_mass_spectrum() override;
   void calculate_pole_mass_spectrum() override;
   void calculate_mass_spectrum() override;
   int solve_ewsb_equations_tree_level() override;
   int solve_ewsb_equations() override;
   Eigen::ArrayXd get_tree_level_masses() const override;
   Eigen::ArrayXd get_tree_level_masses_and_mixings() const override;
   const ScalarSingletZ2DMMhInputMsInput_input_parameters& get_input_parameters() const override;
   ScalarSingletZ2DMMhInputMsInput_input_parameters& get_input_parameters() override;
   Eigen::ArrayXd get_extra_parameters() const override;
   const ScalarSingletZ2DMMhInputMsInput_physical& get_physical() const override;
   ScalarSingletZ2DMMhInputMsInput_physical& get_physical() override;
   const Problems& get_problems() const override;
   Problems& get_problems() override;
   void set_tree_level_masses(const Eigen::ArrayXd&) override;
   void set_tree_level_masses_and_mixings(const Eigen::ArrayXd&) override;
   void set_extra_parameters(const Eigen::ArrayXd&) override;
   void set_physical(const ScalarSingletZ2DMMhInputMsInput_physical&) override;
   void clear_problems() override;

   double get_g1() const override { return SUPER(g1); }
   double get_g2() const override { return SUPER(g2); }
   double get_g3() const override { return SUPER(g3); }
   double get_LamS() const override { return SUPER(LamS); }
   double get_LamSH() const override { return SUPER(LamSH); }
   double get_LamH() const override { return SUPER(LamH); }
   const Eigen::Matrix<double,3,3>& get_Yu() const override { return SUPER(Yu); }
   double get_Yu(int i, int k) const override { return SUPER(Yu(i,k)); }
   const Eigen::Matrix<double,3,3>& get_Yd() const override { return SUPER(Yd); }
   double get_Yd(int i, int k) const override { return SUPER(Yd(i,k)); }
   const Eigen::Matrix<double,3,3>& get_Ye() const override { return SUPER(Ye); }
   double get_Ye(int i, int k) const override { return SUPER(Ye(i,k)); }
   double get_muS2() const override { return SUPER(muS2); }
   double get_muH2() const override { return SUPER(muH2); }
   double get_v() const override { return SUPER(v); }

   void set_g1(double g1_) override { g1 = g1_; }
   void set_g2(double g2_) override { g2 = g2_; }
   void set_g3(double g3_) override { g3 = g3_; }
   void set_LamS(double LamS_) override { LamS = LamS_; }
   void set_LamSH(double LamSH_) override { LamSH = LamSH_; }
   void set_LamH(double LamH_) override { LamH = LamH_; }
   void set_Yu(const Eigen::Matrix<double,3,3>& Yu_) override { Yu = Yu_; }
   void set_Yu(int i, int k, const double& value) override { Yu(i,k) = value; }
   void set_Yd(const Eigen::Matrix<double,3,3>& Yd_) override { Yd = Yd_; }
   void set_Yd(int i, int k, const double& value) override { Yd(i,k) = value; }
   void set_Ye(const Eigen::Matrix<double,3,3>& Ye_) override { Ye = Ye_; }
   void set_Ye(int i, int k, const double& value) override { Ye(i,k) = value; }
   void set_muS2(double muS2_) override { muS2 = muS2_; }
   void set_muH2(double muH2_) override { muH2 = muH2_; }
   void set_v(double v_) override { v = v_; }

   double get_MVG() const override { return MVG; }
   double get_MHp() const override { return MHp; }
   double get_Mss() const override { return Mss; }
   const Eigen::Array<double,3,1>& get_MFv() const override { return MFv; }
   double get_MFv(int i) const override { return MFv(i); }
   double get_MAh() const override { return MAh; }
   double get_Mhh() const override { return Mhh; }
   const Eigen::Array<double,3,1>& get_MFd() const override { return MFd; }
   double get_MFd(int i) const override { return MFd(i); }
   const Eigen::Array<double,3,1>& get_MFu() const override { return MFu; }
   double get_MFu(int i) const override { return MFu(i); }
   const Eigen::Array<double,3,1>& get_MFe() const override { return MFe; }
   double get_MFe(int i) const override { return MFe(i); }
   double get_MVWp() const override { return MVWp; }
   double get_MVP() const override { return MVP; }
   double get_MVZ() const override { return MVZ; }

   

   
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vd() const override { return Vd; }
   std::complex<double> get_Vd(int i, int k) const override { return Vd(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ud() const override { return Ud; }
   std::complex<double> get_Ud(int i, int k) const override { return Ud(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vu() const override { return Vu; }
   std::complex<double> get_Vu(int i, int k) const override { return Vu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Uu() const override { return Uu; }
   std::complex<double> get_Uu(int i, int k) const override { return Uu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ve() const override { return Ve; }
   std::complex<double> get_Ve(int i, int k) const override { return Ve(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ue() const override { return Ue; }
   std::complex<double> get_Ue(int i, int k) const override { return Ue(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZZ() const override { return ZZ; }
   double get_ZZ(int i, int k) const override { return ZZ(i,k); }




   double get_mass_matrix_VG() const override;
   void calculate_MVG() override;
   double get_mass_matrix_Hp() const override;
   void calculate_MHp() override;
   double get_mass_matrix_ss() const override;
   void calculate_Mss() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Fv() const override;
   void calculate_MFv() override;
   double get_mass_matrix_Ah() const override;
   void calculate_MAh() override;
   double get_mass_matrix_hh() const override;
   void calculate_Mhh() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Fd() const override;
   void calculate_MFd() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Fu() const override;
   void calculate_MFu() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Fe() const override;
   void calculate_MFe() override;
   double get_mass_matrix_VWp() const override;
   void calculate_MVWp() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_VPVZ() const override;
   void calculate_MVPVZ() override;

   double get_ewsb_eq_hh_1() const override;

   double ThetaW() const override;
   double VEV() const override;


private:
   bool force_output{false};              ///< switch to force output of pole masses
   double precision{1.e-4};               ///< mass eigenstate precision
   ScalarSingletZ2DMMhInputMsInput_physical physical{}; ///< contains the pole masses and mixings
   Problems problems{ScalarSingletZ2DMMhInputMsInput_info::model_name,
                     &ScalarSingletZ2DMMhInputMsInput_info::particle_names_getter,
                     &ScalarSingletZ2DMMhInputMsInput_info::parameter_names_getter}; ///< problems

   void clear_tree_level_parameters();
   void copy_tree_level_masses_to_pole_masses();

   // DR-bar masses
   double MVG{};
   double MHp{};
   double Mss{};
   Eigen::Array<double,3,1> MFv{Eigen::Array<double,3,1>::Zero()};
   double MAh{};
   double Mhh{};
   Eigen::Array<double,3,1> MFd{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFu{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFe{Eigen::Array<double,3,1>::Zero()};
   double MVWp{};
   double MVP{};
   double MVZ{};

   // DR-bar mixing matrices
   Eigen::Matrix<std::complex<double>,3,3> Vd{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ud{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Vu{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Uu{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ve{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ue{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<double,2,2> ZZ{Eigen::Matrix<double,2,2>::Zero()};

   // phases

   // extra parameters

};

std::ostream& operator<<(std::ostream&, const ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_decoupling_scheme&);

} // namespace flexiblesusy

#undef SUPER

#endif
