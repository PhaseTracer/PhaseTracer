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
 * @file ScalarSingletZ2DMMhInputMsInput_mass_eigenstates.hpp
 *
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated with FlexibleSUSY 2.6.1 and SARAH 4.14.3 .
 */

#ifndef ScalarSingletZ2DMMhInputMsInput_MASS_EIGENSTATES_H
#define ScalarSingletZ2DMMhInputMsInput_MASS_EIGENSTATES_H

#include "ScalarSingletZ2DMMhInputMsInput_info.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_physical.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_soft_parameters.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_interface.hpp"
#include "loop_corrections.hpp"
#include "threshold_corrections.hpp"
#include "problems.hpp"

#include <iosfwd>
#include <memory>
#include <string>

#include <Eigen/Core>

#define SUPER(p) ScalarSingletZ2DMMhInputMsInput_soft_parameters::p

namespace flexiblesusy {

class ScalarSingletZ2DMMhInputMsInput_ewsb_solver_interface;
/**
 * @class ScalarSingletZ2DMMhInputMsInput_mass_eigenstates
 * @brief model class with routines for determing masses and mixinga and EWSB
 */
class ScalarSingletZ2DMMhInputMsInput_mass_eigenstates
   : public ScalarSingletZ2DMMhInputMsInput_soft_parameters
   , public ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_interface
{
public:
   explicit ScalarSingletZ2DMMhInputMsInput_mass_eigenstates(const ScalarSingletZ2DMMhInputMsInput_input_parameters& input_ = ScalarSingletZ2DMMhInputMsInput_input_parameters());
   ScalarSingletZ2DMMhInputMsInput_mass_eigenstates(const ScalarSingletZ2DMMhInputMsInput_mass_eigenstates&) = default;
   ScalarSingletZ2DMMhInputMsInput_mass_eigenstates(ScalarSingletZ2DMMhInputMsInput_mass_eigenstates&&) = default;
   virtual ~ScalarSingletZ2DMMhInputMsInput_mass_eigenstates() = default;
   ScalarSingletZ2DMMhInputMsInput_mass_eigenstates& operator=(const ScalarSingletZ2DMMhInputMsInput_mass_eigenstates&) = default;
   ScalarSingletZ2DMMhInputMsInput_mass_eigenstates& operator=(ScalarSingletZ2DMMhInputMsInput_mass_eigenstates&&) = default;
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW

   std::unique_ptr<ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_interface> clone() const override;

   /// number of EWSB equations
   static const int number_of_ewsb_equations = 1;

   void calculate_DRbar_masses();
   void calculate_pole_masses();
   void check_pole_masses_for_tachyons();
   virtual void clear() override;
   void clear_DRbar_parameters();
   Eigen::ArrayXd get_DRbar_masses() const;
   Eigen::ArrayXd get_DRbar_masses_and_mixings() const;
   void do_calculate_sm_pole_masses(bool);
   bool do_calculate_sm_pole_masses() const;
   void do_calculate_bsm_pole_masses(bool);
   bool do_calculate_bsm_pole_masses() const;
   void do_force_output(bool);
   bool do_force_output() const;
   void reorder_DRbar_masses();
   void reorder_pole_masses();
   void set_ewsb_iteration_precision(double);
   void set_ewsb_loop_order(int);
   void set_loop_corrections(const Loop_corrections&);
   const Loop_corrections& get_loop_corrections() const;
   void set_threshold_corrections(const Threshold_corrections&);
   const Threshold_corrections& get_threshold_corrections() const;
   void set_DRbar_masses(const Eigen::ArrayXd&);
   void set_DRbar_masses_and_mixings(const Eigen::ArrayXd&);
   void set_pole_mass_loop_order(int);
   int get_pole_mass_loop_order() const;
   double get_ewsb_iteration_precision() const;
   double get_ewsb_loop_order() const;
   void set_ewsb_solver(const std::shared_ptr<ScalarSingletZ2DMMhInputMsInput_ewsb_solver_interface>&);
   int solve_ewsb_tree_level();
   int solve_ewsb_one_loop();
   int solve_ewsb();            ///< solve EWSB at ewsb_loop_order level
   int solve_ewsb_tree_level_custom();
   
   virtual void calculate_spectrum();
   std::string name() const;
   void run_to(double scale, double eps = -1.0) override;
   void print(std::ostream&) const override;
   void set_precision(double);
   double get_precision() const;

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

   double CphhHpconjHp() const;
   double CpbargWpgZHp() const;
   double CpbargZgWpconjHp() const;
   double CpbargWpCgZconjHp() const;
   double CpbargZgWpCHp() const;
   double CpconjHpVPVWp() const;
   double CpconjHpVWpVZ() const;
   double CpAhAhHpconjHp() const;
   double CphhhhHpconjHp() const;
   double CpHpHpconjHpconjHp() const;
   double CpHpssssconjHp() const;
   std::complex<double> CpAhconjHpVWp() const;
   double CphhconjHpVWp() const;
   double CpHpconjHpVP() const;
   double CpHpconjHpVZ() const;
   double CpHpconjHpconjVWpVWp() const;
   std::complex<double> CpHpconjHpVZVZ() const;
   std::complex<double> CpbarFdFuconjHpPR(int gI1, int gI2) const;
   std::complex<double> CpbarFdFuconjHpPL(int gI1, int gI2) const;
   double CpbarFeFvconjHpPR(int , int ) const;
   std::complex<double> CpbarFeFvconjHpPL(int gI1, int gI2) const;
   double Cphhssss() const;
   double CpAhAhssss() const;
   double Cphhhhssss() const;
   double Cpssssssss() const;
   double CpAhAhhh() const;
   std::complex<double> CpbargWpgWpAh() const;
   std::complex<double> CpbargWpCgWpCAh() const;
   double CpAhAhAhAh() const;
   double CpAhAhhhhh() const;
   std::complex<double> CpAhhhVZ() const;
   std::complex<double> CpAhHpconjVWp() const;
   double CpAhAhconjVWpVWp() const;
   std::complex<double> CpAhAhVZVZ() const;
   std::complex<double> CpbarFdFdAhPR(int gI1, int gI2) const;
   std::complex<double> CpbarFdFdAhPL(int gI1, int gI2) const;
   std::complex<double> CpbarFeFeAhPR(int gI1, int gI2) const;
   std::complex<double> CpbarFeFeAhPL(int gI1, int gI2) const;
   std::complex<double> CpbarFuFuAhPR(int gI1, int gI2) const;
   std::complex<double> CpbarFuFuAhPL(int gI1, int gI2) const;
   double Cphhhhhh() const;
   double CphhVZVZ() const;
   double CphhconjVWpVWp() const;
   double CpbargWpgWphh() const;
   double CpbargWpCgWpChh() const;
   double CpbargZgZhh() const;
   double Cphhhhhhhh() const;
   double CphhHpconjVWp() const;
   double CphhhhconjVWpVWp() const;
   std::complex<double> CphhhhVZVZ() const;
   std::complex<double> CpbarFdFdhhPR(int gI1, int gI2) const;
   std::complex<double> CpbarFdFdhhPL(int gI1, int gI2) const;
   std::complex<double> CpbarFeFehhPR(int gI1, int gI2) const;
   std::complex<double> CpbarFeFehhPL(int gI1, int gI2) const;
   std::complex<double> CpbarFuFuhhPR(int gI1, int gI2) const;
   std::complex<double> CpbarFuFuhhPL(int gI1, int gI2) const;
   std::complex<double> CpVGVGVG() const;
   std::complex<double> CpbargGgGVG() const;
   double CpbarFdFdVGPL(int gI1, int gI2) const;
   double CpbarFdFdVGPR(int gI1, int gI2) const;
   double CpbarFuFuVGPL(int gI1, int gI2) const;
   double CpbarFuFuVGPR(int gI1, int gI2) const;
   double CpVGVGVGVG1() const;
   double CpVGVGVGVG2() const;
   double CpVGVGVGVG3() const;
   double CpHpconjVWpVP() const;
   double CpbargWpgWpVP() const;
   double CpbargWpCgWpCVP() const;
   std::complex<double> CpHpconjHpVPVP() const;
   double CpconjVWpVPVWp() const;
   double CpbarFdFdVPPL(int gI1, int gI2) const;
   double CpbarFdFdVPPR(int gI1, int gI2) const;
   double CpbarFeFeVPPL(int gI1, int gI2) const;
   double CpbarFeFeVPPR(int gI1, int gI2) const;
   double CpbarFuFuVPPL(int gI1, int gI2) const;
   double CpbarFuFuVPPR(int gI1, int gI2) const;
   double CpconjVWpVPVPVWp3() const;
   double CpconjVWpVPVPVWp1() const;
   double CpconjVWpVPVPVWp2() const;
   double CpHpconjVWpVZ() const;
   double CpbargWpgWpVZ() const;
   double CpbargWpCgWpCVZ() const;
   double CpconjVWpVWpVZ() const;
   double CpbarFdFdVZPL(int gI1, int gI2) const;
   double CpbarFdFdVZPR(int gI1, int gI2) const;
   double CpbarFeFeVZPL(int gI1, int gI2) const;
   double CpbarFeFeVZPR(int gI1, int gI2) const;
   double CpbarFuFuVZPL(int gI1, int gI2) const;
   double CpbarFuFuVZPR(int gI1, int gI2) const;
   double CpbarFvFvVZPL(int gI1, int gI2) const;
   double CpbarFvFvVZPR(int , int ) const;
   double CpconjVWpVWpVZVZ1() const;
   double CpconjVWpVWpVZVZ2() const;
   double CpconjVWpVWpVZVZ3() const;
   double CpbargPgWpconjVWp() const;
   double CpbargWpCgPconjVWp() const;
   double CpbargWpCgZconjVWp() const;
   double CpbargZgWpconjVWp() const;
   std::complex<double> CpbarFdFuconjVWpPL(int gI1, int gI2) const;
   double CpbarFdFuconjVWpPR(int , int ) const;
   std::complex<double> CpbarFeFvconjVWpPL(int gI1, int gI2) const;
   double CpbarFeFvconjVWpPR(int , int ) const;
   double CpconjVWpconjVWpVWpVWp2() const;
   double CpconjVWpconjVWpVWpVWp1() const;
   double CpconjVWpconjVWpVWpVWp3() const;
   std::complex<double> CpbarUFdFdAhPL(int gO2, int gI1) const;
   std::complex<double> CpbarUFdFdAhPR(int gO1, int gI1) const;
   std::complex<double> CpbarUFdFdhhPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdhhPR(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFdVGPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVGPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFdVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFdVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFuconjHpPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFuconjHpPR(int gO1, int gI2) const;
   double CpbarUFdFuconjVWpPR(int , int ) const;
   std::complex<double> CpbarUFdFuconjVWpPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuAhPL(int gO2, int gI1) const;
   std::complex<double> CpbarUFuFuAhPR(int gO1, int gI1) const;
   std::complex<double> CpbarUFuFdHpPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFdHpPR(int gO1, int gI2) const;
   double CpbarUFuFdVWpPR(int , int ) const;
   std::complex<double> CpbarUFuFdVWpPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuhhPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuhhPR(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVGPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVGPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFeAhPL(int gO2, int gI1) const;
   std::complex<double> CpbarUFeFeAhPR(int gO1, int gI1) const;
   std::complex<double> CpbarUFeFehhPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFeFehhPR(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFeVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFeFeVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFeVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFeFeVZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFvconjHpPL(int gO2, int gI2) const;
   double CpbarUFeFvconjHpPR(int , int ) const;
   double CpbarUFeFvconjVWpPR(int , int ) const;
   double CpbarUFeFvconjVWpPL(int gO1, int gI2) const;
   double CpbarFvFeHpPL(int , int ) const;
   std::complex<double> CpbarFvFeHpPR(int gO1, int gI2) const;
   double CpbarFvFeVWpPR(int , int ) const;
   std::complex<double> CpbarFvFeVWpPL(int gO1, int gI2) const;
   std::complex<double> CpbarFuFdHpPL(int gO2, int gI2) const;
   std::complex<double> CpbarFuFdHpPR(int gO1, int gI2) const;
   double CpbarFuFdVWpPR(int , int ) const;
   std::complex<double> CpbarFuFdVWpPL(int gO1, int gI2) const;
   std::complex<double> self_energy_Hp_1loop(double p ) const;
   std::complex<double> self_energy_ss_1loop(double p ) const;
   std::complex<double> self_energy_Ah_1loop(double p ) const;
   std::complex<double> self_energy_hh_1loop(double p ) const;
   std::complex<double> self_energy_VG_1loop(double p ) const;
   std::complex<double> self_energy_VP_1loop(double p ) const;
   std::complex<double> self_energy_VZ_1loop(double p ) const;
   std::complex<double> self_energy_VWp_1loop(double p ) const;
   std::complex<double> self_energy_Fd_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_1(double p) const;
   std::complex<double> self_energy_Fd_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_PR(double p) const;
   std::complex<double> self_energy_Fd_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_PL(double p) const;
   std::complex<double> self_energy_Fu_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_1(double p) const;
   std::complex<double> self_energy_Fu_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PR(double p) const;
   std::complex<double> self_energy_Fu_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PL(double p) const;
   std::complex<double> self_energy_Fe_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_1(double p) const;
   std::complex<double> self_energy_Fe_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_PR(double p) const;
   std::complex<double> self_energy_Fe_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_PL(double p) const;
   std::complex<double> self_energy_Fv_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fv_1loop_1(double p) const;
   std::complex<double> self_energy_Fv_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fv_1loop_PR(double p) const;
   std::complex<double> self_energy_Fv_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fv_1loop_PL(double p) const;
   std::complex<double> self_energy_Fd_1loop_1_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_1_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fd_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_PR_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fd_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_PL_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fe_1loop_1_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_1_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fe_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_PR_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fe_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_PL_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fu_1loop_1_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_1_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fu_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PR_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fu_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PL_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fu_1loop_1_heavy(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_1_heavy(double p) const;
   std::complex<double> self_energy_Fu_1loop_PR_heavy(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PR_heavy(double p) const;
   std::complex<double> self_energy_Fu_1loop_PL_heavy(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PL_heavy(double p) const;
   std::complex<double> tadpole_hh_1loop() const;


   /// calculates the tadpoles at current loop order
   void tadpole_equations(double[number_of_ewsb_equations]) const;
   /// calculates the tadpoles at current loop order
   Eigen::Matrix<double,number_of_ewsb_equations,1> tadpole_equations() const;
   /// calculates the tadpoles divided by VEVs at current loop order
   Eigen::Matrix<double,number_of_ewsb_equations,1> tadpole_equations_over_vevs() const;







   void calculate_MVG_pole();
   void calculate_Mss_pole();
   void calculate_MFv_pole();
   void calculate_Mhh_pole();
   void calculate_MVP_pole();
   void calculate_MVZ_pole();
   void calculate_MFd_pole();
   void calculate_MFu_pole();
   void calculate_MFe_pole();
   void calculate_MVWp_pole();
   double calculate_MVWp_pole(double);
   double calculate_MVZ_pole(double);

   double calculate_MFv_DRbar(double, int) const;
   double calculate_MFe_DRbar(double, int) const;
   double calculate_MFu_DRbar(double, int) const;
   double calculate_MFd_DRbar(double, int) const;
   double calculate_MVP_DRbar(double) const;
   double calculate_MVZ_DRbar(double) const;
   double calculate_MVWp_DRbar(double) const;

   double ThetaW() const override;
   double VEV() const override;


private:
   int ewsb_loop_order{4};           ///< loop order for EWSB
   int pole_mass_loop_order{4};      ///< loop order for pole masses
   bool calculate_sm_pole_masses{false};  ///< switch to calculate the pole masses of the Standard Model particles
   bool calculate_bsm_pole_masses{true};  ///< switch to calculate the pole masses of the BSM particles
   bool force_output{false};              ///< switch to force output of pole masses
   double precision{1.e-4};               ///< RG running precision
   double ewsb_iteration_precision{1.e-5};///< precision goal of EWSB solution
   ScalarSingletZ2DMMhInputMsInput_physical physical{}; ///< contains the pole masses and mixings
   mutable Problems problems{ScalarSingletZ2DMMhInputMsInput_info::model_name,
                             &ScalarSingletZ2DMMhInputMsInput_info::particle_names_getter,
                             &ScalarSingletZ2DMMhInputMsInput_info::parameter_names_getter}; ///< problems
   Loop_corrections loop_corrections{}; ///< used pole mass corrections
   std::shared_ptr<ScalarSingletZ2DMMhInputMsInput_ewsb_solver_interface> ewsb_solver{};
   Threshold_corrections threshold_corrections{}; ///< used threshold corrections

   int get_number_of_ewsb_iterations() const;
   int get_number_of_mass_iterations() const;
   void copy_DRbar_masses_to_pole_masses();

   // Passarino-Veltman loop functions
   double A0(double) const noexcept;
   double B0(double, double, double) const noexcept;
   double B1(double, double, double) const noexcept;
   double B00(double, double, double) const noexcept;
   double B22(double, double, double) const noexcept;
   double H0(double, double, double) const noexcept;
   double F0(double, double, double) const noexcept;
   double G0(double, double, double) const noexcept;

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

std::ostream& operator<<(std::ostream&, const ScalarSingletZ2DMMhInputMsInput_mass_eigenstates&);

} // namespace flexiblesusy

#undef SUPER

#endif
