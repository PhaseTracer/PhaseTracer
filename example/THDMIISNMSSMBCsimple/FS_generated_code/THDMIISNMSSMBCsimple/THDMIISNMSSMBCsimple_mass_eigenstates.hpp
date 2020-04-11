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

// File generated at Sat 11 Apr 2020 12:54:02

/**
 * @file THDMIISNMSSMBCsimple_mass_eigenstates.hpp
 *
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Sat 11 Apr 2020 12:54:02 with FlexibleSUSY
 * 2.4.2 (git commit: a94199e5620b8684f5d30d0eece5757a5a72c4a4) and SARAH 4.14.3 .
 */

#ifndef THDMIISNMSSMBCsimple_MASS_EIGENSTATES_H
#define THDMIISNMSSMBCsimple_MASS_EIGENSTATES_H

#include "THDMIISNMSSMBCsimple_info.hpp"
#include "THDMIISNMSSMBCsimple_physical.hpp"
#include "THDMIISNMSSMBCsimple_soft_parameters.hpp"
#include "loop_corrections.hpp"
#include "threshold_corrections.hpp"
#include "problems.hpp"

#include <iosfwd>
#include <memory>
#include <string>

#include <Eigen/Core>

namespace flexiblesusy {

class THDMIISNMSSMBCsimple_ewsb_solver_interface;
/**
 * @class THDMIISNMSSMBCsimple_mass_eigenstates
 * @brief model class with routines for determing masses and mixinga and EWSB
 */
class THDMIISNMSSMBCsimple_mass_eigenstates : public THDMIISNMSSMBCsimple_soft_parameters {
public:
   explicit THDMIISNMSSMBCsimple_mass_eigenstates(const THDMIISNMSSMBCsimple_input_parameters& input_ = THDMIISNMSSMBCsimple_input_parameters());
   THDMIISNMSSMBCsimple_mass_eigenstates(const THDMIISNMSSMBCsimple_mass_eigenstates&) = default;
   THDMIISNMSSMBCsimple_mass_eigenstates(THDMIISNMSSMBCsimple_mass_eigenstates&&) = default;
   virtual ~THDMIISNMSSMBCsimple_mass_eigenstates() = default;
   THDMIISNMSSMBCsimple_mass_eigenstates& operator=(const THDMIISNMSSMBCsimple_mass_eigenstates&) = default;
   THDMIISNMSSMBCsimple_mass_eigenstates& operator=(THDMIISNMSSMBCsimple_mass_eigenstates&&) = default;

   /// number of EWSB equations
   static const int number_of_ewsb_equations = 3;

   void calculate_DRbar_masses();
   void calculate_pole_masses();
   void check_pole_masses_for_tachyons();
   virtual void clear() override;
   void clear_DRbar_parameters();
   Eigen::ArrayXd get_DRbar_masses() const;
   Eigen::ArrayXd get_DRbar_masses_and_mixings() const;
   Eigen::ArrayXd get_extra_parameters() const;
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
   void set_extra_parameters(const Eigen::ArrayXd&);
   void set_pole_mass_loop_order(int);
   int get_pole_mass_loop_order() const;
   void set_physical(const THDMIISNMSSMBCsimple_physical&);
   double get_ewsb_iteration_precision() const;
   double get_ewsb_loop_order() const;
   const THDMIISNMSSMBCsimple_physical& get_physical() const;
   THDMIISNMSSMBCsimple_physical& get_physical();
   const Problems& get_problems() const;
   Problems& get_problems();
   void set_ewsb_solver(const std::shared_ptr<THDMIISNMSSMBCsimple_ewsb_solver_interface>&);
   int solve_ewsb_tree_level();
   int solve_ewsb_one_loop();
   int solve_ewsb();            ///< solve EWSB at ewsb_loop_order level

   void calculate_spectrum();
   void clear_problems();
   std::string name() const;
   void run_to(double scale, double eps = -1.0) override;
   void print(std::ostream& out = std::cerr) const override;
   void set_precision(double);
   double get_precision() const;


   double get_MVG() const { return MVG; }
   const Eigen::Array<double,3,1>& get_MFv() const { return MFv; }
   double get_MFv(int i) const { return MFv(i); }
   const Eigen::Array<double,3,1>& get_Mhh() const { return Mhh; }
   double get_Mhh(int i) const { return Mhh(i); }
   const Eigen::Array<double,3,1>& get_MAh() const { return MAh; }
   double get_MAh(int i) const { return MAh(i); }
   const Eigen::Array<double,2,1>& get_MHm() const { return MHm; }
   double get_MHm(int i) const { return MHm(i); }
   const Eigen::Array<double,3,1>& get_MFd() const { return MFd; }
   double get_MFd(int i) const { return MFd(i); }
   const Eigen::Array<double,3,1>& get_MFu() const { return MFu; }
   double get_MFu(int i) const { return MFu(i); }
   const Eigen::Array<double,3,1>& get_MFe() const { return MFe; }
   double get_MFe(int i) const { return MFe(i); }
   double get_MVWm() const { return MVWm; }
   double get_MVP() const { return MVP; }
   double get_MVZ() const { return MVZ; }

   
   Eigen::Array<double,1,1> get_MChargedHiggs() const;

   Eigen::Array<double,2,1> get_MPseudoscalarHiggs() const;

   const Eigen::Matrix<double,3,3>& get_ZH() const { return ZH; }
   double get_ZH(int i, int k) const { return ZH(i,k); }
   const Eigen::Matrix<double,3,3>& get_ZA() const { return ZA; }
   double get_ZA(int i, int k) const { return ZA(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZP() const { return ZP; }
   double get_ZP(int i, int k) const { return ZP(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vd() const { return Vd; }
   std::complex<double> get_Vd(int i, int k) const { return Vd(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ud() const { return Ud; }
   std::complex<double> get_Ud(int i, int k) const { return Ud(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vu() const { return Vu; }
   std::complex<double> get_Vu(int i, int k) const { return Vu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Uu() const { return Uu; }
   std::complex<double> get_Uu(int i, int k) const { return Uu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ve() const { return Ve; }
   std::complex<double> get_Ve(int i, int k) const { return Ve(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ue() const { return Ue; }
   std::complex<double> get_Ue(int i, int k) const { return Ue(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZZ() const { return ZZ; }
   double get_ZZ(int i, int k) const { return ZZ(i,k); }




   double get_mass_matrix_VG() const;
   void calculate_MVG();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fv() const;
   void calculate_MFv();
   Eigen::Matrix<double,3,3> get_mass_matrix_hh() const;
   void calculate_Mhh();
   Eigen::Matrix<double,3,3> get_mass_matrix_Ah() const;
   void calculate_MAh();
   Eigen::Matrix<double,2,2> get_mass_matrix_Hm() const;
   void calculate_MHm();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fd() const;
   void calculate_MFd();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fu() const;
   void calculate_MFu();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fe() const;
   void calculate_MFe();
   double get_mass_matrix_VWm() const;
   void calculate_MVWm();
   Eigen::Matrix<double,2,2> get_mass_matrix_VPVZ() const;
   void calculate_MVPVZ();

   double get_ewsb_eq_hh_1() const;
   double get_ewsb_eq_hh_2() const;
   double get_ewsb_eq_hh_3() const;

   std::complex<double> CpbargWmgWmUhh(int gO1) const;
   std::complex<double> CpbargWmCgWmCUhh(int gO1) const;
   std::complex<double> CpbargZgZUhh(int gO1) const;
   std::complex<double> CpUhhVZVZ(int gO2) const;
   std::complex<double> CpUhhconjVWmVWm(int gO2) const;
   std::complex<double> CpUhhUhhconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpUhhUhhVZVZ(int gO1, int gO2) const;
   std::complex<double> CpUhhUhhHmconjHm(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUhhHmconjHm(int gO2, int gI2, int gI1) const;
   std::complex<double> CpAhAhUhhUhh(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUhhUhh(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpAhAhUhh(int gI1, int gI2, int gO2) const;
   std::complex<double> CpAhhhUhh(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhhhUhh(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFdFdUhhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFdFdUhhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarFeFeUhhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFeFeUhhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarFuFuUhhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFuFuUhhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpUhhHmconjVWm(int gO2, int gI2) const;
   std::complex<double> CpAhUhhVZ(int gI2, int gO2) const;
   std::complex<double> CpbargWmgWmUAh(int gO1) const;
   std::complex<double> CpbargWmCgWmCUAh(int gO1) const;
   std::complex<double> CpUAhUAhconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpUAhUAhVZVZ(int gO1, int gO2) const;
   std::complex<double> CpUAhUAhHmconjHm(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhHmconjHm(int gO2, int gI2, int gI1) const;
   std::complex<double> CpAhAhUAhUAh(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpUAhUAhhhhh(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpAhAhUAh(int gI1, int gI2, int gO2) const;
   std::complex<double> CpAhUAhhh(int gI2, int gO2, int gI1) const;
   std::complex<double> CpUAhhhhh(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarFdFdUAhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFdFdUAhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarFeFeUAhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFeFeUAhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarFuFuUAhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFuFuUAhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpUAhHmconjVWm(int gO2, int gI2) const;
   std::complex<double> CpUAhhhVZ(int gO2, int gI2) const;
   std::complex<double> CpbargWmgZUHm(int gO2) const;
   std::complex<double> CpbargZgWmconjUHm(int gO1) const;
   std::complex<double> CpbargWmCgZconjUHm(int gO1) const;
   std::complex<double> CpbargZgWmCUHm(int gO2) const;
   std::complex<double> CpconjUHmVPVWm(int gO2) const;
   std::complex<double> CpconjUHmVWmVZ(int gO2) const;
   std::complex<double> CpUHmconjUHmconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpUHmconjUHmVZVZ(int gO1, int gO2) const;
   std::complex<double> CpHmUHmconjHmconjUHm(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpAhHmconjUHm(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhHmconjUHm(int gI2, int gI1, int gO2) const;
   std::complex<double> CpAhAhUHmconjUHm(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUHmconjUHm(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpbarFuFdconjUHmPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFuFdconjUHmPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarFvFeconjUHmPR(int gI1, int gI2, int gO2) const;
   double CpbarFvFeconjUHmPL(int , int , int ) const;
   std::complex<double> CpHmconjUHmVP(int gI2, int gO2) const;
   std::complex<double> CpHmconjUHmVZ(int gI2, int gO2) const;
   std::complex<double> CpAhconjUHmVWm(int gI2, int gO2) const;
   std::complex<double> CphhconjUHmVWm(int gI2, int gO2) const;
   std::complex<double> CpVGVGVG() const;
   std::complex<double> CpbargGgGVG() const;
   double CpbarFdFdVGPL(int gI1, int gI2) const;
   double CpbarFdFdVGPR(int gI1, int gI2) const;
   double CpbarFuFuVGPL(int gI1, int gI2) const;
   double CpbarFuFuVGPR(int gI1, int gI2) const;
   double CpVGVGVGVG1() const;
   double CpVGVGVGVG2() const;
   double CpVGVGVGVG3() const;
   double CpbargWmgWmVP() const;
   double CpbargWmCgWmCVP() const;
   double CpconjVWmVPVWm() const;
   std::complex<double> CpHmconjHmVPVP(int gI1, int gI2) const;
   double CpHmconjHmVP(int gI2, int gI1) const;
   double CpbarFdFdVPPL(int gI1, int gI2) const;
   double CpbarFdFdVPPR(int gI1, int gI2) const;
   double CpbarFeFeVPPL(int gI1, int gI2) const;
   double CpbarFeFeVPPR(int gI1, int gI2) const;
   double CpbarFuFuVPPL(int gI1, int gI2) const;
   double CpbarFuFuVPPR(int gI1, int gI2) const;
   std::complex<double> CpHmconjVWmVP(int gI2) const;
   double CpconjVWmVPVPVWm3() const;
   double CpconjVWmVPVPVWm1() const;
   double CpconjVWmVPVPVWm2() const;
   double CpbargWmgWmVZ() const;
   double CpbargWmCgWmCVZ() const;
   double CpconjVWmVWmVZ() const;
   std::complex<double> CpHmconjHmVZVZ(int gI1, int gI2) const;
   double CpHmconjHmVZ(int gI2, int gI1) const;
   std::complex<double> CpAhAhVZVZ(int gI1, int gI2) const;
   std::complex<double> CphhhhVZVZ(int gI1, int gI2) const;
   std::complex<double> CpAhhhVZ(int gI2, int gI1) const;
   double CpbarFdFdVZPL(int gI1, int gI2) const;
   double CpbarFdFdVZPR(int gI1, int gI2) const;
   double CpbarFeFeVZPL(int gI1, int gI2) const;
   double CpbarFeFeVZPR(int gI1, int gI2) const;
   double CpbarFuFuVZPL(int gI1, int gI2) const;
   double CpbarFuFuVZPR(int gI1, int gI2) const;
   double CpbarFvFvVZPL(int gI1, int gI2) const;
   double CpbarFvFvVZPR(int , int ) const;
   std::complex<double> CpHmconjVWmVZ(int gI2) const;
   std::complex<double> CphhVZVZ(int gI2) const;
   double CpconjVWmVWmVZVZ1() const;
   double CpconjVWmVWmVZVZ2() const;
   double CpconjVWmVWmVZVZ3() const;
   double CpbargPgWmconjVWm() const;
   double CpbargWmCgPconjVWm() const;
   double CpbargWmCgZconjVWm() const;
   double CpbargZgWmconjVWm() const;
   std::complex<double> CpHmconjHmconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CpAhHmconjVWm(int gI2, int gI1) const;
   std::complex<double> CphhHmconjVWm(int gI2, int gI1) const;
   std::complex<double> CpAhAhconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CphhhhconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CpbarFuFdconjVWmPL(int gI1, int gI2) const;
   double CpbarFuFdconjVWmPR(int , int ) const;
   std::complex<double> CpbarFvFeconjVWmPL(int gI1, int gI2) const;
   double CpbarFvFeconjVWmPR(int , int ) const;
   std::complex<double> CphhconjVWmVWm(int gI2) const;
   double CpconjVWmconjVWmVWmVWm2() const;
   double CpconjVWmconjVWmVWmVWm1() const;
   double CpconjVWmconjVWmVWmVWm3() const;
   std::complex<double> CpbarUFdFuHmPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFdFuHmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFdFdAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUFdFdAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUFdFdhhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFdFdhhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFdFdVGPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVGPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFdVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFdVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVZPL(int gO1, int gI2) const;
   double CpbarUFdFuVWmPR(int , int ) const;
   std::complex<double> CpbarUFdFuVWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFdconjHmPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFuFdconjHmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFuFuAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUFuFuAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUFuFuhhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFuFuhhPR(int gO1, int gI2, int gI1) const;
   double CpbarUFuFdconjVWmPR(int , int ) const;
   std::complex<double> CpbarUFuFdconjVWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVGPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVGPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFvHmPL(int gO2, int gI2, int gI1) const;
   double CpbarUFeFvHmPR(int , int , int ) const;
   std::complex<double> CpbarUFeFeAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUFeFeAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUFeFehhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFeFehhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFeFeVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFeFeVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFeVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFeFeVZPL(int gO1, int gI2) const;
   double CpbarUFeFvVWmPR(int , int ) const;
   double CpbarUFeFvVWmPL(int gO1, int gI2) const;
   double CpbarFvFeconjHmPL(int , int , int ) const;
   std::complex<double> CpbarFvFeconjHmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFdFuHmPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFdFuHmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFdFdAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarFdFdAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarFdFdhhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFdFdhhPR(int gO1, int gI2, int gI1) const;
   double CpbarFdFuVWmPR(int , int ) const;
   std::complex<double> CpbarFdFuVWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarFeFvHmPL(int gO2, int gI2, int gI1) const;
   double CpbarFeFvHmPR(int , int , int ) const;
   std::complex<double> CpbarFeFeAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarFeFeAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarFeFehhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFeFehhPR(int gO1, int gI2, int gI1) const;
   double CpbarFeFvVWmPR(int , int ) const;
   std::complex<double> CpbarFeFvVWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarFuFdconjHmPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFuFdconjHmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFuFuAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarFuFuAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarFuFuhhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFuFuhhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> self_energy_hh_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_hh_1loop(double p) const;
   std::complex<double> self_energy_Ah_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Ah_1loop(double p) const;
   std::complex<double> self_energy_Hm_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_Hm_1loop(double p) const;
   std::complex<double> self_energy_VG_1loop(double p ) const;
   std::complex<double> self_energy_VP_1loop(double p ) const;
   std::complex<double> self_energy_VZ_1loop(double p ) const;
   std::complex<double> self_energy_VWm_1loop(double p ) const;
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
   std::complex<double> tadpole_hh_1loop(int gO1) const;


   /// calculates the tadpoles at current loop order
   void tadpole_equations(double[number_of_ewsb_equations]) const;
   /// calculates the tadpoles at current loop order
   Eigen::Matrix<double,number_of_ewsb_equations,1> tadpole_equations() const;
   /// calculates the tadpoles divided by VEVs at current loop order
   Eigen::Matrix<double,number_of_ewsb_equations,1> tadpole_equations_over_vevs() const;







   void calculate_MVG_pole();
   void calculate_MFv_pole();
   void calculate_MVP_pole();
   void calculate_MVZ_pole();
   void calculate_Mhh_pole();
   void calculate_MAh_pole();
   void calculate_MHm_pole();
   void calculate_MFd_pole();
   void calculate_MFu_pole();
   void calculate_MFe_pole();
   void calculate_MVWm_pole();
   double calculate_MVWm_pole(double);
   double calculate_MVZ_pole(double);

   double calculate_MFv_DRbar(double, int) const;
   double calculate_MFe_DRbar(double, int) const;
   double calculate_MFu_DRbar(double, int) const;
   double calculate_MFd_DRbar(double, int) const;
   double calculate_MVP_DRbar(double);
   double calculate_MVZ_DRbar(double);
   double calculate_MVWm_DRbar(double);

   double v() const;
   double Betax() const;
   double ThetaW() const;
   double VEV() const;


private:
   int ewsb_loop_order{4};           ///< loop order for EWSB
   int pole_mass_loop_order{4};      ///< loop order for pole masses
   bool calculate_sm_pole_masses{false};  ///< switch to calculate the pole masses of the Standard Model particles
   bool calculate_bsm_pole_masses{true};  ///< switch to calculate the pole masses of the BSM particles
   bool force_output{false};              ///< switch to force output of pole masses
   double precision{1.e-4};               ///< RG running precision
   double ewsb_iteration_precision{1.e-5};///< precision goal of EWSB solution
   THDMIISNMSSMBCsimple_physical physical{}; ///< contains the pole masses and mixings
   Problems problems{THDMIISNMSSMBCsimple_info::model_name,
                     &THDMIISNMSSMBCsimple_info::particle_names_getter,
                     &THDMIISNMSSMBCsimple_info::parameter_names_getter}; ///< problems
   Loop_corrections loop_corrections{}; ///< used pole mass corrections
   std::shared_ptr<THDMIISNMSSMBCsimple_ewsb_solver_interface> ewsb_solver{};
   Threshold_corrections threshold_corrections{}; ///< used threshold corrections

   int get_number_of_ewsb_iterations() const;
   int get_number_of_mass_iterations() const;
   int solve_ewsb_tree_level_custom();
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

   // DR-bar mixing matrices
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

   // phases

   // extra parameters

};

std::ostream& operator<<(std::ostream&, const THDMIISNMSSMBCsimple_mass_eigenstates&);

} // namespace flexiblesusy

#endif
