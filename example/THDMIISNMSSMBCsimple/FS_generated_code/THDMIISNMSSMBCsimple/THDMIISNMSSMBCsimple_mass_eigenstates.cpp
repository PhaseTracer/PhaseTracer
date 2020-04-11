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
 * @file THDMIISNMSSMBCsimple_mass_eigenstates.cpp
 * @brief implementation of the THDMIISNMSSMBCsimple model class
 *
 * Contains the definition of the THDMIISNMSSMBCsimple model class methods
 * which solve EWSB and calculate pole masses and mixings from MSbar
 * parameters.
 *
 * This file was generated at Sat 11 Apr 2020 12:54:02 with FlexibleSUSY
 * 2.4.2 (git commit: a94199e5620b8684f5d30d0eece5757a5a72c4a4) and SARAH 4.14.3 .
 */

#include "THDMIISNMSSMBCsimple_mass_eigenstates.hpp"
#include "THDMIISNMSSMBCsimple_ewsb_solver_interface.hpp"
#include "config.h"
#include "eigen_utils.hpp"
#include "ewsb_solver.hpp"
#include "wrappers.hpp"
#include "linalg2.hpp"
#include "numerics2.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "pv.hpp"
#include "raii.hpp"
#include "functors.hpp"

#ifdef ENABLE_THREADS
#include "thread_pool.hpp"
#endif

#ifdef ENABLE_TWO_SCALE_SOLVER
#include "THDMIISNMSSMBCsimple_two_scale_ewsb_solver.hpp"
#endif






#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <algorithm>
#include <stdexcept>

namespace flexiblesusy {

#define STRINGIFY(s) XSTRINGIFY(s)
#define XSTRINGIFY(s) #s
#define CLASSNAME THDMIISNMSSMBCsimple_mass_eigenstates

#define PHYSICAL(parameter) physical.parameter
#define INPUT(parameter) model.get_input().parameter
#define LOCALINPUT(parameter) input.parameter
#define MODELPARAMETER(parameter) model.get_##parameter()
#define EXTRAPARAMETER(parameter) model.get_##parameter()

#define HIGGS_2LOOP_CORRECTION_AT_AS       loop_corrections.higgs_at_as
#define HIGGS_2LOOP_CORRECTION_AB_AS       loop_corrections.higgs_ab_as
#define HIGGS_2LOOP_CORRECTION_AT_AT       loop_corrections.higgs_at_at
#define HIGGS_2LOOP_CORRECTION_ATAU_ATAU   loop_corrections.higgs_atau_atau
#define TOP_POLE_QCD_CORRECTION            loop_corrections.top_qcd
#define HIGGS_3LOOP_CORRECTION_AT_AS_AS    loop_corrections.higgs_at_as_as
#define HIGGS_3LOOP_CORRECTION_AB_AS_AS    loop_corrections.higgs_ab_as_as
#define HIGGS_3LOOP_SCHEME                 loop_corrections.higgs_3L_scheme
#define HIGGS_3LOOP_CORRECTION_AT_AT_AS    loop_corrections.higgs_at_at_as
#define HIGGS_3LOOP_CORRECTION_AT_AT_AT    loop_corrections.higgs_at_at_at
#define HIGGS_4LOOP_CORRECTION_AT_AS_AS_AS loop_corrections.higgs_at_as_as_as

CLASSNAME::CLASSNAME(const THDMIISNMSSMBCsimple_input_parameters& input_)
   : THDMIISNMSSMBCsimple_soft_parameters(input_)
#if defined(ENABLE_TWO_SCALE_SOLVER)
   , ewsb_solver(new THDMIISNMSSMBCsimple_ewsb_solver<Two_scale>())
#endif
{
}

void CLASSNAME::do_calculate_sm_pole_masses(bool flag)
{
   calculate_sm_pole_masses = flag;
}

bool CLASSNAME::do_calculate_sm_pole_masses() const
{
   return calculate_sm_pole_masses;
}

void CLASSNAME::do_calculate_bsm_pole_masses(bool flag)
{
   calculate_bsm_pole_masses = flag;
}

bool CLASSNAME::do_calculate_bsm_pole_masses() const
{
   return calculate_bsm_pole_masses;
}

void CLASSNAME::do_force_output(bool flag)
{
   force_output = flag;
}

bool CLASSNAME::do_force_output() const
{
   return force_output;
}

void CLASSNAME::set_ewsb_loop_order(int loop_order)
{
   ewsb_loop_order = loop_order;
   if (ewsb_solver) {
      ewsb_solver->set_loop_order(ewsb_loop_order);
   }
}

void CLASSNAME::set_loop_corrections(const Loop_corrections& loop_corrections_)
{
   loop_corrections = loop_corrections_;
}

const Loop_corrections& CLASSNAME::get_loop_corrections() const
{
   return loop_corrections;
}

void CLASSNAME::set_threshold_corrections(const Threshold_corrections& tc)
{
   threshold_corrections = tc;
}

const Threshold_corrections& CLASSNAME::get_threshold_corrections() const
{
   return threshold_corrections;
}

int CLASSNAME::get_number_of_ewsb_iterations() const
{
   return static_cast<int>(std::abs(-log10(ewsb_iteration_precision) * 10));
}

int CLASSNAME::get_number_of_mass_iterations() const
{
   return static_cast<int>(std::abs(-log10(precision) * 10));
}

void CLASSNAME::set_precision(double precision_)
{
   precision = precision_;
   ewsb_iteration_precision = precision_;
   if (ewsb_solver) {
      ewsb_solver->set_precision(precision_);
   }
}

void CLASSNAME::set_pole_mass_loop_order(int loop_order)
{
   pole_mass_loop_order = loop_order;
}

int CLASSNAME::get_pole_mass_loop_order() const
{
   return pole_mass_loop_order;
}

void CLASSNAME::set_ewsb_iteration_precision(double precision)
{
   ewsb_iteration_precision = precision;
   if (ewsb_solver) {
      ewsb_solver->set_precision(precision);
   }
}

double CLASSNAME::get_ewsb_iteration_precision() const
{
   return ewsb_iteration_precision;
}

double CLASSNAME::get_precision() const
{
   return precision;
}

double CLASSNAME::get_ewsb_loop_order() const
{
   return ewsb_loop_order;
}

const THDMIISNMSSMBCsimple_physical& CLASSNAME::get_physical() const
{
   return physical;
}

THDMIISNMSSMBCsimple_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const THDMIISNMSSMBCsimple_physical& physical_)
{
   physical = physical_;
}

const Problems& CLASSNAME::get_problems() const
{
   return problems;
}

Problems& CLASSNAME::get_problems()
{
   return problems;
}

void CLASSNAME::set_ewsb_solver(const std::shared_ptr<THDMIISNMSSMBCsimple_ewsb_solver_interface>& solver)
{
   ewsb_solver = solver;
}

/**
 * Method which calculates the tadpoles at the current loop order.
 *
 * @param tadpole array of tadpole
 */
void CLASSNAME::tadpole_equations(double tadpole[number_of_ewsb_equations]) const
{
   const auto tadpole_(tadpole_equations());
   std::copy(tadpole_.data(), tadpole_.data() + number_of_ewsb_equations, tadpole);
}

/**
 * Method which calculates the tadpoles at the current loop order.
 *
 * @return array of tadpoles
 */
Eigen::Matrix<double,CLASSNAME::number_of_ewsb_equations,1> CLASSNAME::tadpole_equations() const
{
   Eigen::Matrix<double,number_of_ewsb_equations,1> tadpole(
      Eigen::Matrix<double,number_of_ewsb_equations,1>::Zero());

   tadpole[0] = get_ewsb_eq_hh_1();
   tadpole[1] = get_ewsb_eq_hh_2();
   tadpole[2] = get_ewsb_eq_hh_3();

   if (ewsb_loop_order > 0) {
      tadpole[0] -= Re(tadpole_hh_1loop(0));
      tadpole[1] -= Re(tadpole_hh_1loop(1));
      tadpole[2] -= Re(tadpole_hh_1loop(2));

      if (ewsb_loop_order > 1) {

      }
   }

   return tadpole;
}

/**
 * This function returns the vector of tadpoles, each divided by the
 * corresponding VEV.  Thus, the returned tadpoles have the dimension
 * GeV^2 each.
 *
 * @return vector of tadpoles
 */
Eigen::Matrix<double,CLASSNAME::number_of_ewsb_equations,1> CLASSNAME::tadpole_equations_over_vevs() const
{
   auto tadpole = tadpole_equations();

   tadpole[0] /= v1;
   tadpole[1] /= v2;
   tadpole[2] /= vS;


   return tadpole;
}

int CLASSNAME::solve_ewsb_tree_level_custom()
{
   int error = EWSB_solver::SUCCESS;



   return error;
}

int CLASSNAME::solve_ewsb_tree_level()
{
   if (!ewsb_solver) {
      throw SetupError(STRINGIFY(CLASSNAME) "::solve_ewsb_tree_level: "
                       "no EWSB solver set");
   }

   const int old_loop_order = ewsb_solver->get_loop_order();
   const auto save_loop_order = make_raii_guard(
      [this, old_loop_order] () {
         this->ewsb_solver->set_loop_order(old_loop_order);
      });

   const int old_iterations = ewsb_solver->get_number_of_iterations();
   const auto save_iterations = make_raii_guard(
      [this, old_iterations] () {
         this->ewsb_solver->set_number_of_iterations(old_iterations);
      });

   const double old_precision = ewsb_solver->get_precision();
   const auto save_precision = make_raii_guard(
      [this, old_precision] () {
         this->ewsb_solver->set_precision(old_precision);
      });

   ewsb_solver->set_loop_order(0);
   ewsb_solver->set_number_of_iterations(get_number_of_ewsb_iterations());
   ewsb_solver->set_precision(ewsb_iteration_precision);

   return ewsb_solver->solve(*this);
}

int CLASSNAME::solve_ewsb_one_loop()
{
   if (!ewsb_solver) {
      throw SetupError(STRINGIFY(CLASSNAME) "::solve_ewsb_one_loop: "
                       "no EWSB solver set");
   }

   const int old_loop_order = ewsb_solver->get_loop_order();
   const auto save_loop_order = make_raii_guard(
      [this, old_loop_order] () {
         this->ewsb_solver->set_loop_order(old_loop_order);
      });

   const int old_iterations = ewsb_solver->get_number_of_iterations();
   const auto save_iterations = make_raii_guard(
      [this, old_iterations] () {
         this->ewsb_solver->set_number_of_iterations(old_iterations);
      });

   const double old_precision = ewsb_solver->get_precision();
   const auto save_precision = make_raii_guard(
      [this, old_precision] () {
         this->ewsb_solver->set_precision(old_precision);
      });

   ewsb_solver->set_loop_order(1);
   ewsb_solver->set_number_of_iterations(get_number_of_ewsb_iterations());
   ewsb_solver->set_precision(ewsb_iteration_precision);

   return ewsb_solver->solve(*this);
}

int CLASSNAME::solve_ewsb()
{
   if (!ewsb_solver) {
      throw SetupError(STRINGIFY(CLASSNAME) "::solve_ewsb: "
                       "no EWSB solver set");
   }

   VERBOSE_MSG("\t\tSolving THDMIISNMSSMBCsimple EWSB at " << ewsb_loop_order << "-loop order");

   const int old_loop_order = ewsb_solver->get_loop_order();
   const auto save_loop_order = make_raii_guard(
      [this, old_loop_order] () {
         this->ewsb_solver->set_loop_order(old_loop_order);
      });

   const int old_iterations = ewsb_solver->get_number_of_iterations();
   const auto save_iterations = make_raii_guard(
      [this, old_iterations] () {
         this->ewsb_solver->set_number_of_iterations(old_iterations);
      });

   const double old_precision = ewsb_solver->get_precision();
   const auto save_precision = make_raii_guard(
      [this, old_precision] () {
         this->ewsb_solver->set_precision(old_precision);
      });

   ewsb_solver->set_loop_order(ewsb_loop_order);
   ewsb_solver->set_number_of_iterations(get_number_of_ewsb_iterations());
   ewsb_solver->set_precision(ewsb_iteration_precision);

   return ewsb_solver->solve(*this);
}

void CLASSNAME::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "THDMIISNMSSMBCsimple\n"
           "========================================\n";
   THDMIISNMSSMBCsimple_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level MSbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHm = " << MHm.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MVWm = " << MVWm << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';

   ostr << "----------------------------------------\n"
           "tree-level MSbar mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "ZH = " << ZH << '\n';
   ostr << "ZA = " << ZA << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "Vd = " << Vd << '\n';
   ostr << "Ud = " << Ud << '\n';
   ostr << "Vu = " << Vu << '\n';
   ostr << "Uu = " << Uu << '\n';
   ostr << "Ve = " << Ve << '\n';
   ostr << "Ue = " << Ue << '\n';
   ostr << "ZZ = " << ZZ << '\n';

   physical.print(ostr);
}

/**
 * wrapper routines for passarino Veltman functions
 * @note: They take squared arguments!
 */

double CLASSNAME::A0(double m) const noexcept
{
   return passarino_veltman::ReA0(m, Sqr(get_scale()));
}

double CLASSNAME::B0(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReB0(p, m1, m2, Sqr(get_scale()));
}

double CLASSNAME::B1(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReB1(p, m1, m2, Sqr(get_scale()));
}

double CLASSNAME::B00(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReB00(p, m1, m2, Sqr(get_scale()));
}

double CLASSNAME::B22(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReB22(p, m1, m2, Sqr(get_scale()));
}

double CLASSNAME::H0(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReH0(p, m1, m2, Sqr(get_scale()));
}

double CLASSNAME::F0(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReF0(p, m1, m2, Sqr(get_scale()));
}

double CLASSNAME::G0(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReG0(p, m1, m2, Sqr(get_scale()));
}

/**
 * routine which finds the MSbar mass eigenstates and mixings.
 */
void CLASSNAME::calculate_DRbar_masses()
{
   const auto save_M112_raii = make_raii_save(M112);
   const auto save_M222_raii = make_raii_save(M222);
   const auto save_M332_raii = make_raii_save(M332);

   const bool has_no_ewsb_flag = problems.no_ewsb();
   const auto save_ewsb_flag = make_raii_guard(
      [this, has_no_ewsb_flag] () {
         if (has_no_ewsb_flag) {
            this->problems.flag_no_ewsb_tree_level();
         } else {
            this->problems.unflag_no_ewsb_tree_level();
         }
      }
   );
   problems.unflag_no_ewsb_tree_level();
   solve_ewsb_tree_level();
#ifdef ENABLE_VERBOSE
   if (problems.no_ewsb()) {
      WARNING("solving EWSB at 0-loop order failed");
   }
#endif

   calculate_MVPVZ();
   calculate_MVWm();
   calculate_MFe();
   calculate_MFu();
   calculate_MFd();
   calculate_MHm();
   calculate_MAh();
   calculate_Mhh();
   calculate_MFv();
   calculate_MVG();

}

/**
 * routine which finds the pole mass eigenstates and mixings.
 */
void CLASSNAME::calculate_pole_masses()
{
#ifdef ENABLE_THREADS
   Thread_pool tp(std::min(std::thread::hardware_concurrency(), 11u));

   if (calculate_bsm_pole_masses) {
      tp.run_task([this] () { calculate_MAh_pole(); });
      tp.run_task([this] () { calculate_Mhh_pole(); });
      tp.run_task([this] () { calculate_MHm_pole(); });
   }

   if (calculate_sm_pole_masses) {
      tp.run_task([this] () { calculate_MVG_pole(); });
      tp.run_task([this] () { calculate_MFv_pole(); });
      tp.run_task([this] () { calculate_MVP_pole(); });
      tp.run_task([this] () { calculate_MVZ_pole(); });
      tp.run_task([this] () { calculate_MFd_pole(); });
      tp.run_task([this] () { calculate_MFu_pole(); });
      tp.run_task([this] () { calculate_MFe_pole(); });
      tp.run_task([this] () { calculate_MVWm_pole(); });
   }

#else
   if (calculate_bsm_pole_masses) {
      calculate_MAh_pole();
      calculate_Mhh_pole();
      calculate_MHm_pole();
   }

   if (calculate_sm_pole_masses) {
      calculate_MVG_pole();
      calculate_MFv_pole();
      calculate_MVP_pole();
      calculate_MVZ_pole();
      calculate_MFd_pole();
      calculate_MFu_pole();
      calculate_MFe_pole();
      calculate_MVWm_pole();
   }

#endif
}

void CLASSNAME::copy_DRbar_masses_to_pole_masses()
{
   PHYSICAL(MVG) = MVG;
   PHYSICAL(MFv) = MFv;
   PHYSICAL(Mhh) = Mhh;
   PHYSICAL(ZH) = ZH;
   PHYSICAL(MAh) = MAh;
   PHYSICAL(ZA) = ZA;
   PHYSICAL(MHm) = MHm;
   PHYSICAL(ZP) = ZP;
   PHYSICAL(MFd) = MFd;
   PHYSICAL(Vd) = Vd;
   PHYSICAL(Ud) = Ud;
   PHYSICAL(MFu) = MFu;
   PHYSICAL(Vu) = Vu;
   PHYSICAL(Uu) = Uu;
   PHYSICAL(MFe) = MFe;
   PHYSICAL(Ve) = Ve;
   PHYSICAL(Ue) = Ue;
   PHYSICAL(MVWm) = MVWm;
   PHYSICAL(MVP) = MVP;
   PHYSICAL(MVZ) = MVZ;

}

/**
 * reorders MSbar masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_DRbar_masses()
{
   move_goldstone_to(0, MVZ, MAh, ZA);
   move_goldstone_to(0, MVWm, MHm, ZP);

}

/**
 * reorders pole masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_pole_masses()
{
   move_goldstone_to(0, MVZ, PHYSICAL(MAh), PHYSICAL(ZA));
   move_goldstone_to(0, MVWm, PHYSICAL(MHm), PHYSICAL(ZP));

}

/**
 * Checks the pole masses for tachyons
 */
void CLASSNAME::check_pole_masses_for_tachyons()
{
   if (PHYSICAL(Mhh).tail<3>().minCoeff() < 0.) problems.flag_pole_tachyon(THDMIISNMSSMBCsimple_info::hh);
   if (PHYSICAL(MAh).tail<2>().minCoeff() < 0.) problems.flag_pole_tachyon(THDMIISNMSSMBCsimple_info::Ah);
   if (PHYSICAL(MHm).tail<1>().minCoeff() < 0.) problems.flag_pole_tachyon(THDMIISNMSSMBCsimple_info::Hm);

}

/**
 * calculates spectrum for model once the MSbar parameters at
 * at low energies are known
 */
void CLASSNAME::calculate_spectrum()
{
   calculate_DRbar_masses();
   if (pole_mass_loop_order > 0)
      calculate_pole_masses();

   // move goldstone bosons to the front
   reorder_DRbar_masses();
   if (pole_mass_loop_order == 0)
      copy_DRbar_masses_to_pole_masses();
   else
      reorder_pole_masses();

   check_pole_masses_for_tachyons();

   if (problems.have_problem() && !force_output) {
      clear_DRbar_parameters();
      physical.clear();
   }
}

void CLASSNAME::clear_DRbar_parameters()
{
   MVG = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   Mhh = Eigen::Matrix<double,3,1>::Zero();
   ZH = Eigen::Matrix<double,3,3>::Zero();
   MAh = Eigen::Matrix<double,3,1>::Zero();
   ZA = Eigen::Matrix<double,3,3>::Zero();
   MHm = Eigen::Matrix<double,2,1>::Zero();
   ZP = Eigen::Matrix<double,2,2>::Zero();
   MFd = Eigen::Matrix<double,3,1>::Zero();
   Vd = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ud = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   Vu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Uu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   Ve = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ue = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MVWm = 0.;
   MVP = 0.;
   MVZ = 0.;



}

void CLASSNAME::clear_problems()
{
   problems.clear();
}

void CLASSNAME::clear()
{
   THDMIISNMSSMBCsimple_soft_parameters::clear();
   clear_DRbar_parameters();
   physical.clear();
   problems.clear();
}

void CLASSNAME::set_DRbar_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MFv(0) = pars(1);
   MFv(1) = pars(2);
   MFv(2) = pars(3);
   Mhh(0) = pars(4);
   Mhh(1) = pars(5);
   Mhh(2) = pars(6);
   MAh(0) = pars(7);
   MAh(1) = pars(8);
   MAh(2) = pars(9);
   MHm(0) = pars(10);
   MHm(1) = pars(11);
   MFd(0) = pars(12);
   MFd(1) = pars(13);
   MFd(2) = pars(14);
   MFu(0) = pars(15);
   MFu(1) = pars(16);
   MFu(2) = pars(17);
   MFe(0) = pars(18);
   MFe(1) = pars(19);
   MFe(2) = pars(20);
   MVWm = pars(21);
   MVP = pars(22);
   MVZ = pars(23);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses() const
{
   Eigen::ArrayXd pars(24);

   pars(0) = MVG;
   pars(1) = MFv(0);
   pars(2) = MFv(1);
   pars(3) = MFv(2);
   pars(4) = Mhh(0);
   pars(5) = Mhh(1);
   pars(6) = Mhh(2);
   pars(7) = MAh(0);
   pars(8) = MAh(1);
   pars(9) = MAh(2);
   pars(10) = MHm(0);
   pars(11) = MHm(1);
   pars(12) = MFd(0);
   pars(13) = MFd(1);
   pars(14) = MFd(2);
   pars(15) = MFu(0);
   pars(16) = MFu(1);
   pars(17) = MFu(2);
   pars(18) = MFe(0);
   pars(19) = MFe(1);
   pars(20) = MFe(2);
   pars(21) = MVWm;
   pars(22) = MVP;
   pars(23) = MVZ;

   return pars;
}

void CLASSNAME::set_DRbar_masses_and_mixings(const Eigen::ArrayXd& pars)
{
   set_DRbar_masses(pars);

   ZH(0,0) = pars(24);
   ZH(0,1) = pars(25);
   ZH(0,2) = pars(26);
   ZH(1,0) = pars(27);
   ZH(1,1) = pars(28);
   ZH(1,2) = pars(29);
   ZH(2,0) = pars(30);
   ZH(2,1) = pars(31);
   ZH(2,2) = pars(32);
   ZA(0,0) = pars(33);
   ZA(0,1) = pars(34);
   ZA(0,2) = pars(35);
   ZA(1,0) = pars(36);
   ZA(1,1) = pars(37);
   ZA(1,2) = pars(38);
   ZA(2,0) = pars(39);
   ZA(2,1) = pars(40);
   ZA(2,2) = pars(41);
   ZP(0,0) = pars(42);
   ZP(0,1) = pars(43);
   ZP(1,0) = pars(44);
   ZP(1,1) = pars(45);
   Vd(0,0) = std::complex<double>(pars(46), pars(47));
   Vd(0,1) = std::complex<double>(pars(48), pars(49));
   Vd(0,2) = std::complex<double>(pars(50), pars(51));
   Vd(1,0) = std::complex<double>(pars(52), pars(53));
   Vd(1,1) = std::complex<double>(pars(54), pars(55));
   Vd(1,2) = std::complex<double>(pars(56), pars(57));
   Vd(2,0) = std::complex<double>(pars(58), pars(59));
   Vd(2,1) = std::complex<double>(pars(60), pars(61));
   Vd(2,2) = std::complex<double>(pars(62), pars(63));
   Ud(0,0) = std::complex<double>(pars(64), pars(65));
   Ud(0,1) = std::complex<double>(pars(66), pars(67));
   Ud(0,2) = std::complex<double>(pars(68), pars(69));
   Ud(1,0) = std::complex<double>(pars(70), pars(71));
   Ud(1,1) = std::complex<double>(pars(72), pars(73));
   Ud(1,2) = std::complex<double>(pars(74), pars(75));
   Ud(2,0) = std::complex<double>(pars(76), pars(77));
   Ud(2,1) = std::complex<double>(pars(78), pars(79));
   Ud(2,2) = std::complex<double>(pars(80), pars(81));
   Vu(0,0) = std::complex<double>(pars(82), pars(83));
   Vu(0,1) = std::complex<double>(pars(84), pars(85));
   Vu(0,2) = std::complex<double>(pars(86), pars(87));
   Vu(1,0) = std::complex<double>(pars(88), pars(89));
   Vu(1,1) = std::complex<double>(pars(90), pars(91));
   Vu(1,2) = std::complex<double>(pars(92), pars(93));
   Vu(2,0) = std::complex<double>(pars(94), pars(95));
   Vu(2,1) = std::complex<double>(pars(96), pars(97));
   Vu(2,2) = std::complex<double>(pars(98), pars(99));
   Uu(0,0) = std::complex<double>(pars(100), pars(101));
   Uu(0,1) = std::complex<double>(pars(102), pars(103));
   Uu(0,2) = std::complex<double>(pars(104), pars(105));
   Uu(1,0) = std::complex<double>(pars(106), pars(107));
   Uu(1,1) = std::complex<double>(pars(108), pars(109));
   Uu(1,2) = std::complex<double>(pars(110), pars(111));
   Uu(2,0) = std::complex<double>(pars(112), pars(113));
   Uu(2,1) = std::complex<double>(pars(114), pars(115));
   Uu(2,2) = std::complex<double>(pars(116), pars(117));
   Ve(0,0) = std::complex<double>(pars(118), pars(119));
   Ve(0,1) = std::complex<double>(pars(120), pars(121));
   Ve(0,2) = std::complex<double>(pars(122), pars(123));
   Ve(1,0) = std::complex<double>(pars(124), pars(125));
   Ve(1,1) = std::complex<double>(pars(126), pars(127));
   Ve(1,2) = std::complex<double>(pars(128), pars(129));
   Ve(2,0) = std::complex<double>(pars(130), pars(131));
   Ve(2,1) = std::complex<double>(pars(132), pars(133));
   Ve(2,2) = std::complex<double>(pars(134), pars(135));
   Ue(0,0) = std::complex<double>(pars(136), pars(137));
   Ue(0,1) = std::complex<double>(pars(138), pars(139));
   Ue(0,2) = std::complex<double>(pars(140), pars(141));
   Ue(1,0) = std::complex<double>(pars(142), pars(143));
   Ue(1,1) = std::complex<double>(pars(144), pars(145));
   Ue(1,2) = std::complex<double>(pars(146), pars(147));
   Ue(2,0) = std::complex<double>(pars(148), pars(149));
   Ue(2,1) = std::complex<double>(pars(150), pars(151));
   Ue(2,2) = std::complex<double>(pars(152), pars(153));
   ZZ(0,0) = pars(154);
   ZZ(0,1) = pars(155);
   ZZ(1,0) = pars(156);
   ZZ(1,1) = pars(157);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses_and_mixings() const
{
   Eigen::ArrayXd pars(get_DRbar_masses());

   pars.conservativeResize(158);

   pars(24) = ZH(0,0);
   pars(25) = ZH(0,1);
   pars(26) = ZH(0,2);
   pars(27) = ZH(1,0);
   pars(28) = ZH(1,1);
   pars(29) = ZH(1,2);
   pars(30) = ZH(2,0);
   pars(31) = ZH(2,1);
   pars(32) = ZH(2,2);
   pars(33) = ZA(0,0);
   pars(34) = ZA(0,1);
   pars(35) = ZA(0,2);
   pars(36) = ZA(1,0);
   pars(37) = ZA(1,1);
   pars(38) = ZA(1,2);
   pars(39) = ZA(2,0);
   pars(40) = ZA(2,1);
   pars(41) = ZA(2,2);
   pars(42) = ZP(0,0);
   pars(43) = ZP(0,1);
   pars(44) = ZP(1,0);
   pars(45) = ZP(1,1);
   pars(46) = Re(Vd(0,0));
   pars(47) = Im(Vd(0,0));
   pars(48) = Re(Vd(0,1));
   pars(49) = Im(Vd(0,1));
   pars(50) = Re(Vd(0,2));
   pars(51) = Im(Vd(0,2));
   pars(52) = Re(Vd(1,0));
   pars(53) = Im(Vd(1,0));
   pars(54) = Re(Vd(1,1));
   pars(55) = Im(Vd(1,1));
   pars(56) = Re(Vd(1,2));
   pars(57) = Im(Vd(1,2));
   pars(58) = Re(Vd(2,0));
   pars(59) = Im(Vd(2,0));
   pars(60) = Re(Vd(2,1));
   pars(61) = Im(Vd(2,1));
   pars(62) = Re(Vd(2,2));
   pars(63) = Im(Vd(2,2));
   pars(64) = Re(Ud(0,0));
   pars(65) = Im(Ud(0,0));
   pars(66) = Re(Ud(0,1));
   pars(67) = Im(Ud(0,1));
   pars(68) = Re(Ud(0,2));
   pars(69) = Im(Ud(0,2));
   pars(70) = Re(Ud(1,0));
   pars(71) = Im(Ud(1,0));
   pars(72) = Re(Ud(1,1));
   pars(73) = Im(Ud(1,1));
   pars(74) = Re(Ud(1,2));
   pars(75) = Im(Ud(1,2));
   pars(76) = Re(Ud(2,0));
   pars(77) = Im(Ud(2,0));
   pars(78) = Re(Ud(2,1));
   pars(79) = Im(Ud(2,1));
   pars(80) = Re(Ud(2,2));
   pars(81) = Im(Ud(2,2));
   pars(82) = Re(Vu(0,0));
   pars(83) = Im(Vu(0,0));
   pars(84) = Re(Vu(0,1));
   pars(85) = Im(Vu(0,1));
   pars(86) = Re(Vu(0,2));
   pars(87) = Im(Vu(0,2));
   pars(88) = Re(Vu(1,0));
   pars(89) = Im(Vu(1,0));
   pars(90) = Re(Vu(1,1));
   pars(91) = Im(Vu(1,1));
   pars(92) = Re(Vu(1,2));
   pars(93) = Im(Vu(1,2));
   pars(94) = Re(Vu(2,0));
   pars(95) = Im(Vu(2,0));
   pars(96) = Re(Vu(2,1));
   pars(97) = Im(Vu(2,1));
   pars(98) = Re(Vu(2,2));
   pars(99) = Im(Vu(2,2));
   pars(100) = Re(Uu(0,0));
   pars(101) = Im(Uu(0,0));
   pars(102) = Re(Uu(0,1));
   pars(103) = Im(Uu(0,1));
   pars(104) = Re(Uu(0,2));
   pars(105) = Im(Uu(0,2));
   pars(106) = Re(Uu(1,0));
   pars(107) = Im(Uu(1,0));
   pars(108) = Re(Uu(1,1));
   pars(109) = Im(Uu(1,1));
   pars(110) = Re(Uu(1,2));
   pars(111) = Im(Uu(1,2));
   pars(112) = Re(Uu(2,0));
   pars(113) = Im(Uu(2,0));
   pars(114) = Re(Uu(2,1));
   pars(115) = Im(Uu(2,1));
   pars(116) = Re(Uu(2,2));
   pars(117) = Im(Uu(2,2));
   pars(118) = Re(Ve(0,0));
   pars(119) = Im(Ve(0,0));
   pars(120) = Re(Ve(0,1));
   pars(121) = Im(Ve(0,1));
   pars(122) = Re(Ve(0,2));
   pars(123) = Im(Ve(0,2));
   pars(124) = Re(Ve(1,0));
   pars(125) = Im(Ve(1,0));
   pars(126) = Re(Ve(1,1));
   pars(127) = Im(Ve(1,1));
   pars(128) = Re(Ve(1,2));
   pars(129) = Im(Ve(1,2));
   pars(130) = Re(Ve(2,0));
   pars(131) = Im(Ve(2,0));
   pars(132) = Re(Ve(2,1));
   pars(133) = Im(Ve(2,1));
   pars(134) = Re(Ve(2,2));
   pars(135) = Im(Ve(2,2));
   pars(136) = Re(Ue(0,0));
   pars(137) = Im(Ue(0,0));
   pars(138) = Re(Ue(0,1));
   pars(139) = Im(Ue(0,1));
   pars(140) = Re(Ue(0,2));
   pars(141) = Im(Ue(0,2));
   pars(142) = Re(Ue(1,0));
   pars(143) = Im(Ue(1,0));
   pars(144) = Re(Ue(1,1));
   pars(145) = Im(Ue(1,1));
   pars(146) = Re(Ue(1,2));
   pars(147) = Im(Ue(1,2));
   pars(148) = Re(Ue(2,0));
   pars(149) = Im(Ue(2,0));
   pars(150) = Re(Ue(2,1));
   pars(151) = Im(Ue(2,1));
   pars(152) = Re(Ue(2,2));
   pars(153) = Im(Ue(2,2));
   pars(154) = ZZ(0,0);
   pars(155) = ZZ(0,1);
   pars(156) = ZZ(1,0);
   pars(157) = ZZ(1,1);


   return pars;
}

void CLASSNAME::set_extra_parameters(const Eigen::ArrayXd& pars)
{

}

Eigen::ArrayXd CLASSNAME::get_extra_parameters() const
{
   return Eigen::ArrayXd();

}

std::string CLASSNAME::name() const
{
   return "THDMIISNMSSMBCsimple";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   THDMIISNMSSMBCsimple_soft_parameters::run_to(scale, eps);
}


Eigen::Array<double,1,1> CLASSNAME::get_MChargedHiggs() const
{
   Eigen::Array<double,1,1> MHm_goldstone;
   MHm_goldstone(0) = MVWm;

   return remove_if_equal(MHm, MHm_goldstone);
}

Eigen::Array<double,2,1> CLASSNAME::get_MPseudoscalarHiggs() const
{
   Eigen::Array<double,1,1> MAh_goldstone;
   MAh_goldstone(0) = MVZ;

   return remove_if_equal(MAh, MAh_goldstone);
}




double CLASSNAME::get_mass_matrix_VG() const
{

   const double mass_matrix_VG = Re(0);

   return mass_matrix_VG;
}

void CLASSNAME::calculate_MVG()
{

   const auto mass_matrix_VG = get_mass_matrix_VG();
   MVG = mass_matrix_VG;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fv() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fv;

   mass_matrix_Fv(0,0) = 0;
   mass_matrix_Fv(0,1) = 0;
   mass_matrix_Fv(0,2) = 0;
   mass_matrix_Fv(1,1) = 0;
   mass_matrix_Fv(1,2) = 0;
   mass_matrix_Fv(2,2) = 0;

   Symmetrize(mass_matrix_Fv);

   return mass_matrix_Fv;
}

void CLASSNAME::calculate_MFv()
{

   MFv.setConstant(0);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_hh() const
{

   Eigen::Matrix<double,3,3> mass_matrix_hh;

   mass_matrix_hh(0,0) = M112 + 1.5*Lambda1*Sqr(v1) + 0.5*Lambda3*Sqr(v2) + 0.5
      *Lambda5*Sqr(vS);
   mass_matrix_hh(0,1) = Lambda3*v1*v2 - 0.35355339059327373*M123*vS -
      0.35355339059327373*vS*Conj(M123) + 0.25*Lambda7*Sqr(vS) + 0.25*Conj(
      Lambda7)*Sqr(vS);
   mass_matrix_hh(0,2) = -0.35355339059327373*M123*v2 + Lambda5*v1*vS + 0.5*
      Lambda7*v2*vS + 0.5*v2*vS*Conj(Lambda7) - 0.35355339059327373*v2*Conj(
      M123);
   mass_matrix_hh(1,1) = M222 + 0.5*Lambda3*Sqr(v1) + 1.5*Lambda2*Sqr(v2) + 0.5
      *Lambda6*Sqr(vS);
   mass_matrix_hh(1,2) = -0.35355339059327373*M123*v1 + 0.5*Lambda7*v1*vS +
      Lambda6*v2*vS + 0.5*v1*vS*Conj(Lambda7) - 0.35355339059327373*v1*Conj(
      M123);
   mass_matrix_hh(2,2) = M332 + 0.5*Lambda7*v1*v2 - 0.7071067811865475*M5*vS +
      0.5*v1*v2*Conj(Lambda7) - 0.7071067811865475*vS*Conj(M5) + 0.5*Lambda5*
      Sqr(v1) + 0.5*Lambda6*Sqr(v2) + 3*Lambda8*Sqr(vS);

   Symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH, eigenvalue_error);
   problems.flag_bad_mass(THDMIISNMSSMBCsimple_info::hh, eigenvalue_error >
      precision * Abs(Mhh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);
#endif
   normalize_to_interval(ZH);


   if (Mhh.minCoeff() < 0.) {
      problems.flag_running_tachyon(THDMIISNMSSMBCsimple_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Ah() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = M112 + 0.5*Lambda1*Sqr(v1) + 0.3872983346207417*g1*g2*
      Cos(ThetaW())*Sin(ThetaW())*Sqr(v1) + 0.5*Lambda3*Sqr(v2) + 0.5*Lambda5*
      Sqr(vS) + 0.25*Sqr(g2)*Sqr(v1)*Sqr(Cos(ThetaW())) + 0.15*Sqr(g1)*Sqr(v1)*
      Sqr(Sin(ThetaW()));
   mass_matrix_Ah(0,1) = 0.35355339059327373*M123*vS + 0.35355339059327373*vS*
      Conj(M123) - 0.3872983346207417*g1*g2*v1*v2*Cos(ThetaW())*Sin(ThetaW()) -
      0.25*Lambda7*Sqr(vS) - 0.25*Conj(Lambda7)*Sqr(vS) - 0.25*v1*v2*Sqr(g2)*
      Sqr(Cos(ThetaW())) - 0.15*v1*v2*Sqr(g1)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(0,2) = 0.35355339059327373*M123*v2 + 0.5*Lambda7*v2*vS + 0.5*
      v2*vS*Conj(Lambda7) + 0.35355339059327373*v2*Conj(M123);
   mass_matrix_Ah(1,1) = M222 + 0.5*Lambda3*Sqr(v1) + 0.5*Lambda2*Sqr(v2) +
      0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(v2) + 0.5*
      Lambda6*Sqr(vS) + 0.25*Sqr(g2)*Sqr(v2)*Sqr(Cos(ThetaW())) + 0.15*Sqr(g1)*
      Sqr(v2)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(1,2) = 0.35355339059327373*M123*v1 + 0.5*Lambda7*v1*vS + 0.5*
      v1*vS*Conj(Lambda7) + 0.35355339059327373*v1*Conj(M123);
   mass_matrix_Ah(2,2) = M332 - 0.5*Lambda7*v1*v2 + 0.7071067811865475*M5*vS -
      0.5*v1*v2*Conj(Lambda7) + 0.7071067811865475*vS*Conj(M5) + 0.5*Lambda5*
      Sqr(v1) + 0.5*Lambda6*Sqr(v2) + Lambda8*Sqr(vS);

   Symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA, eigenvalue_error);
   problems.flag_bad_mass(THDMIISNMSSMBCsimple_info::Ah, eigenvalue_error >
      precision * Abs(MAh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);
#endif
   normalize_to_interval(ZA);


   if (MAh.minCoeff() < 0.) {
      problems.flag_running_tachyon(THDMIISNMSSMBCsimple_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Hm() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Hm;

   mass_matrix_Hm(0,0) = M112 + 0.5*Lambda1*Sqr(v1) + 0.25*Sqr(g2)*Sqr(v1) +
      0.5*Lambda3*Sqr(v2) - 0.5*Lambda4*Sqr(v2) + 0.5*Lambda5*Sqr(vS);
   mass_matrix_Hm(0,1) = -0.5*Lambda4*v1*v2 + 0.7071067811865475*vS*Conj(M123)
      - 0.25*v1*v2*Sqr(g2) - 0.5*Conj(Lambda7)*Sqr(vS);
   mass_matrix_Hm(1,0) = -0.5*Lambda4*v1*v2 + 0.7071067811865475*M123*vS - 0.25
      *v1*v2*Sqr(g2) - 0.5*Lambda7*Sqr(vS);
   mass_matrix_Hm(1,1) = M222 + 0.5*Lambda3*Sqr(v1) - 0.5*Lambda4*Sqr(v1) + 0.5
      *Lambda2*Sqr(v2) + 0.25*Sqr(g2)*Sqr(v2) + 0.5*Lambda6*Sqr(vS);

   return mass_matrix_Hm;
}

void CLASSNAME::calculate_MHm()
{
   const auto mass_matrix_Hm(get_mass_matrix_Hm());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Hm, MHm, ZP, eigenvalue_error);
   problems.flag_bad_mass(THDMIISNMSSMBCsimple_info::Hm, eigenvalue_error >
      precision * Abs(MHm(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Hm, MHm, ZP);
#endif
   normalize_to_interval(ZP);


   if (MHm.minCoeff() < 0.) {
      problems.flag_running_tachyon(THDMIISNMSSMBCsimple_info::Hm);
   }

   MHm = AbsSqrt(MHm);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fd() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fd;

   mass_matrix_Fd(0,0) = -0.7071067811865475*v1*Yd(0,0);
   mass_matrix_Fd(0,1) = -0.7071067811865475*v1*Yd(1,0);
   mass_matrix_Fd(0,2) = -0.7071067811865475*v1*Yd(2,0);
   mass_matrix_Fd(1,0) = -0.7071067811865475*v1*Yd(0,1);
   mass_matrix_Fd(1,1) = -0.7071067811865475*v1*Yd(1,1);
   mass_matrix_Fd(1,2) = -0.7071067811865475*v1*Yd(2,1);
   mass_matrix_Fd(2,0) = -0.7071067811865475*v1*Yd(0,2);
   mass_matrix_Fd(2,1) = -0.7071067811865475*v1*Yd(1,2);
   mass_matrix_Fd(2,2) = -0.7071067811865475*v1*Yd(2,2);

   return mass_matrix_Fd;
}

void CLASSNAME::calculate_MFd()
{
   const auto mass_matrix_Fd(get_mass_matrix_Fd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud, eigenvalue_error);
   problems.flag_bad_mass(THDMIISNMSSMBCsimple_info::Fd, eigenvalue_error >
      precision * Abs(MFd(0)));
#else
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fu() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fu;

   mass_matrix_Fu(0,0) = -0.7071067811865475*v2*Yu(0,0);
   mass_matrix_Fu(0,1) = -0.7071067811865475*v2*Yu(1,0);
   mass_matrix_Fu(0,2) = -0.7071067811865475*v2*Yu(2,0);
   mass_matrix_Fu(1,0) = -0.7071067811865475*v2*Yu(0,1);
   mass_matrix_Fu(1,1) = -0.7071067811865475*v2*Yu(1,1);
   mass_matrix_Fu(1,2) = -0.7071067811865475*v2*Yu(2,1);
   mass_matrix_Fu(2,0) = -0.7071067811865475*v2*Yu(0,2);
   mass_matrix_Fu(2,1) = -0.7071067811865475*v2*Yu(1,2);
   mass_matrix_Fu(2,2) = -0.7071067811865475*v2*Yu(2,2);

   return mass_matrix_Fu;
}

void CLASSNAME::calculate_MFu()
{
   const auto mass_matrix_Fu(get_mass_matrix_Fu());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu, eigenvalue_error);
   problems.flag_bad_mass(THDMIISNMSSMBCsimple_info::Fu, eigenvalue_error >
      precision * Abs(MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fe() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fe;

   mass_matrix_Fe(0,0) = -0.7071067811865475*v1*Ye(0,0);
   mass_matrix_Fe(0,1) = -0.7071067811865475*v1*Ye(1,0);
   mass_matrix_Fe(0,2) = -0.7071067811865475*v1*Ye(2,0);
   mass_matrix_Fe(1,0) = -0.7071067811865475*v1*Ye(0,1);
   mass_matrix_Fe(1,1) = -0.7071067811865475*v1*Ye(1,1);
   mass_matrix_Fe(1,2) = -0.7071067811865475*v1*Ye(2,1);
   mass_matrix_Fe(2,0) = -0.7071067811865475*v1*Ye(0,2);
   mass_matrix_Fe(2,1) = -0.7071067811865475*v1*Ye(1,2);
   mass_matrix_Fe(2,2) = -0.7071067811865475*v1*Ye(2,2);

   return mass_matrix_Fe;
}

void CLASSNAME::calculate_MFe()
{
   const auto mass_matrix_Fe(get_mass_matrix_Fe());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue, eigenvalue_error);
   problems.flag_bad_mass(THDMIISNMSSMBCsimple_info::Fe, eigenvalue_error >
      precision * Abs(MFe(0)));
#else
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue);
#endif

}

double CLASSNAME::get_mass_matrix_VWm() const
{

   const double mass_matrix_VWm = Re(0.25*Sqr(g2)*(Sqr(v1) + Sqr(v2)));

   return mass_matrix_VWm;
}

void CLASSNAME::calculate_MVWm()
{

   const auto mass_matrix_VWm = get_mass_matrix_VWm();
   MVWm = mass_matrix_VWm;

   if (MVWm < 0.) {
      problems.flag_running_tachyon(THDMIISNMSSMBCsimple_info::VWm);
   }

   MVWm = AbsSqrt(MVWm);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_VPVZ() const
{

   Eigen::Matrix<double,2,2> mass_matrix_VPVZ;

   mass_matrix_VPVZ(0,0) = 0.15*Sqr(g1)*Sqr(v1) + 0.15*Sqr(g1)*Sqr(v2);
   mass_matrix_VPVZ(0,1) = -0.19364916731037085*g1*g2*Sqr(v1) -
      0.19364916731037085*g1*g2*Sqr(v2);
   mass_matrix_VPVZ(1,1) = 0.25*Sqr(g2)*Sqr(v1) + 0.25*Sqr(g2)*Sqr(v2);

   Symmetrize(mass_matrix_VPVZ);

   return mass_matrix_VPVZ;
}

void CLASSNAME::calculate_MVPVZ()
{
   const auto mass_matrix_VPVZ(get_mass_matrix_VPVZ());
   Eigen::Array<double,2,1> MVPVZ;


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_VPVZ, MVPVZ, ZZ, eigenvalue_error);
#else

   fs_diagonalize_hermitian(mass_matrix_VPVZ, MVPVZ, ZZ);
#endif
   ZZ.transposeInPlace();
   normalize_to_interval(ZZ);


   MVPVZ = AbsSqrt(MVPVZ);

   MVP = 0.;
   MVZ = MVPVZ(1);
}



double CLASSNAME::get_ewsb_eq_hh_1() const
{
   
   double result = Re(M112*v1 - 0.35355339059327373*M123*v2*vS -
      0.35355339059327373*v2*vS*Conj(M123) + 0.5*Lambda1*Cube(v1) + 0.5*Lambda3*v1
      *Sqr(v2) + 0.5*Lambda5*v1*Sqr(vS) + 0.25*Lambda7*v2*Sqr(vS) + 0.25*v2*Conj(
      Lambda7)*Sqr(vS));

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   
   double result = Re(M222*v2 - 0.35355339059327373*M123*v1*vS -
      0.35355339059327373*v1*vS*Conj(M123) + 0.5*Lambda2*Cube(v2) + 0.5*Lambda3*v2
      *Sqr(v1) + 0.25*Lambda7*v1*Sqr(vS) + 0.5*Lambda6*v2*Sqr(vS) + 0.25*v1*Conj(
      Lambda7)*Sqr(vS));

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_3() const
{
   
   double result = Re(-0.35355339059327373*M123*v1*v2 + M332*vS + 0.5*Lambda7*v1*
      v2*vS + 0.5*v1*v2*vS*Conj(Lambda7) - 0.35355339059327373*v1*v2*Conj(M123) +
      Lambda8*Cube(vS) + 0.5*Lambda5*vS*Sqr(v1) + 0.5*Lambda6*vS*Sqr(v2) -
      0.35355339059327373*M5*Sqr(vS) - 0.35355339059327373*Conj(M5)*Sqr(vS));

   return result;
}



std::complex<double> CLASSNAME::CpbargWmgWmUhh(int gO1) const
{
   
   const std::complex<double> result = -0.25*(v1*KroneckerDelta(0,gO1) + v2*
      KroneckerDelta(1,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbargWmCgWmCUhh(int gO1) const
{
   
   const std::complex<double> result = -0.25*(v1*KroneckerDelta(0,gO1) + v2*
      KroneckerDelta(1,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbargZgZUhh(int gO1) const
{
   
   const std::complex<double> result = -0.025*(v1*KroneckerDelta(0,gO1) + v2*
      KroneckerDelta(1,gO1))*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1)
      + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhVZVZ(int gO2) const
{
   
   const std::complex<double> result = 0.05*(v1*KroneckerDelta(0,gO2) + v2*
      KroneckerDelta(1,gO2))*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1)
      + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjVWmVWm(int gO2) const
{
   
   const std::complex<double> result = 0.5*(v1*KroneckerDelta(0,gO2) + v2*
      KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjVWmVWm(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0
      ,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhVZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*(7.745966692414834*g1*
      g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5
      *Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhHmconjHm(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*(2*KroneckerDelta(2,gO1)*KroneckerDelta
      (2,gO2)*(ZP(gI1,1)*(Conj(Lambda7)*ZP(gI2,0) - Lambda6*ZP(gI2,1)) + ZP(gI1,0)
      *(-(Lambda5*ZP(gI2,0)) + Lambda7*ZP(gI2,1))) + KroneckerDelta(1,gO1)*(
      Lambda4*KroneckerDelta(0,gO2)*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) -
      2*KroneckerDelta(1,gO2)*((Lambda3 - Lambda4)*ZP(gI1,0)*ZP(gI2,0) + Lambda2*
      ZP(gI1,1)*ZP(gI2,1))) + KroneckerDelta(0,gO1)*(Lambda4*KroneckerDelta(1,gO2)
      *(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) - 2*KroneckerDelta(0,gO2)*(
      Lambda1*ZP(gI1,0)*ZP(gI2,0) + (Lambda3 - Lambda4)*ZP(gI1,1)*ZP(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhHmconjHm(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.5*(KroneckerDelta(1,gO2)*(ZP(gI1,0)*(2*(-
      Lambda3 + Lambda4)*v2*ZP(gI2,0) + Lambda4*v1*ZP(gI2,1)) + ZP(gI1,1)*(Lambda4
      *v1*ZP(gI2,0) - 2*Lambda2*v2*ZP(gI2,1))) + KroneckerDelta(0,gO2)*(ZP(gI1,1)*
      (Lambda4*v2*ZP(gI2,0) + 2*(-Lambda3 + Lambda4)*v1*ZP(gI2,1)) + ZP(gI1,0)*(-2
      *Lambda1*v1*ZP(gI2,0) + Lambda4*v2*ZP(gI2,1))) - KroneckerDelta(2,gO2)*(ZP(
      gI1,1)*((1.4142135623730951*M123 - 2*Lambda7*vS)*ZP(gI2,0) + 2*Lambda6*vS*ZP
      (gI2,1)) + ZP(gI1,0)*(2*Lambda5*vS*ZP(gI2,0) + (-2*vS*Conj(Lambda7) +
      1.4142135623730951*Conj(M123))*ZP(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUhhUhh(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.5*(-(KroneckerDelta(1,gO1)*(2*
      KroneckerDelta(1,gO2)*(Lambda3*ZA(gI1,0)*ZA(gI2,0) + Lambda2*ZA(gI1,1)*ZA(
      gI2,1) + Lambda6*ZA(gI1,2)*ZA(gI2,2)) + (Lambda7 + Conj(Lambda7))*(-(
      KroneckerDelta(0,gO2)*ZA(gI1,2)*ZA(gI2,2)) + KroneckerDelta(2,gO2)*(ZA(gI1,2
      )*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,2))))) - KroneckerDelta(2,gO1)*(
      KroneckerDelta(2,gO2)*(-(ZA(gI1,1)*((Lambda7 + Conj(Lambda7))*ZA(gI2,0) - 2*
      Lambda6*ZA(gI2,1))) + ZA(gI1,0)*(2*Lambda5*ZA(gI2,0) - (Lambda7 + Conj(
      Lambda7))*ZA(gI2,1)) + 4*Lambda8*ZA(gI1,2)*ZA(gI2,2)) + (Lambda7 + Conj(
      Lambda7))*(KroneckerDelta(1,gO2)*(ZA(gI1,2)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,2))
      + KroneckerDelta(0,gO2)*(ZA(gI1,2)*ZA(gI2,1) + ZA(gI1,1)*ZA(gI2,2)))) -
      KroneckerDelta(0,gO1)*(2*KroneckerDelta(0,gO2)*(Lambda1*ZA(gI1,0)*ZA(gI2,0)
      + Lambda3*ZA(gI1,1)*ZA(gI2,1) + Lambda5*ZA(gI1,2)*ZA(gI2,2)) + (Lambda7 +
      Conj(Lambda7))*(-(KroneckerDelta(1,gO2)*ZA(gI1,2)*ZA(gI2,2)) +
      KroneckerDelta(2,gO2)*(ZA(gI1,2)*ZA(gI2,1) + ZA(gI1,1)*ZA(gI2,2)))));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUhhUhh(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.5*(-(KroneckerDelta(2,gO1)*(
      KroneckerDelta(1,gO2)*(ZH(gI1,2)*((Lambda7 + Conj(Lambda7))*ZH(gI2,0) + 2*
      Lambda6*ZH(gI2,1)) + ((Lambda7 + Conj(Lambda7))*ZH(gI1,0) + 2*Lambda6*ZH(gI1
      ,1))*ZH(gI2,2)) + KroneckerDelta(0,gO2)*(ZH(gI1,2)*(2*Lambda5*ZH(gI2,0) + (
      Lambda7 + Conj(Lambda7))*ZH(gI2,1)) + (2*Lambda5*ZH(gI1,0) + (Lambda7 + Conj
      (Lambda7))*ZH(gI1,1))*ZH(gI2,2)) + KroneckerDelta(2,gO2)*(ZH(gI1,1)*((
      Lambda7 + Conj(Lambda7))*ZH(gI2,0) + 2*Lambda6*ZH(gI2,1)) + ZH(gI1,0)*(2*
      Lambda5*ZH(gI2,0) + (Lambda7 + Conj(Lambda7))*ZH(gI2,1)) + 12*Lambda8*ZH(gI1
      ,2)*ZH(gI2,2)))) - KroneckerDelta(1,gO1)*(KroneckerDelta(2,gO2)*(ZH(gI1,2)*(
      (Lambda7 + Conj(Lambda7))*ZH(gI2,0) + 2*Lambda6*ZH(gI2,1)) + ((Lambda7 +
      Conj(Lambda7))*ZH(gI1,0) + 2*Lambda6*ZH(gI1,1))*ZH(gI2,2)) + 2*
      KroneckerDelta(1,gO2)*(Lambda3*ZH(gI1,0)*ZH(gI2,0) + 3*Lambda2*ZH(gI1,1)*ZH(
      gI2,1) + Lambda6*ZH(gI1,2)*ZH(gI2,2)) + KroneckerDelta(0,gO2)*(2*Lambda3*ZH(
      gI1,1)*ZH(gI2,0) + 2*Lambda3*ZH(gI1,0)*ZH(gI2,1) + (Lambda7 + Conj(Lambda7))
      *ZH(gI1,2)*ZH(gI2,2))) - KroneckerDelta(0,gO1)*(KroneckerDelta(2,gO2)*(ZH(
      gI1,2)*(2*Lambda5*ZH(gI2,0) + (Lambda7 + Conj(Lambda7))*ZH(gI2,1)) + (2*
      Lambda5*ZH(gI1,0) + (Lambda7 + Conj(Lambda7))*ZH(gI1,1))*ZH(gI2,2)) + 2*
      KroneckerDelta(0,gO2)*(3*Lambda1*ZH(gI1,0)*ZH(gI2,0) + Lambda3*ZH(gI1,1)*ZH(
      gI2,1) + Lambda5*ZH(gI1,2)*ZH(gI2,2)) + KroneckerDelta(1,gO2)*(2*Lambda3*ZH(
      gI1,1)*ZH(gI2,0) + 2*Lambda3*ZH(gI1,0)*ZH(gI2,1) + (Lambda7 + Conj(Lambda7))
      *ZH(gI1,2)*ZH(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUhh(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.25*(-(KroneckerDelta(1,gO2)*(4*Lambda2*v2
      *ZA(gI1,1)*ZA(gI2,1) + ZA(gI1,2)*((1.4142135623730951*M123 + 2*Lambda7*vS +
      2*vS*Conj(Lambda7) + 1.4142135623730951*Conj(M123))*ZA(gI2,0) - 2*(Lambda7*
      v1 - 2*Lambda6*v2 + v1*Conj(Lambda7))*ZA(gI2,2)) + ZA(gI1,0)*(4*Lambda3*v2*
      ZA(gI2,0) + (1.4142135623730951*M123 + 2*Lambda7*vS + 2*vS*Conj(Lambda7) +
      1.4142135623730951*Conj(M123))*ZA(gI2,2)))) - KroneckerDelta(0,gO2)*(4*
      Lambda1*v1*ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,2)*((1.4142135623730951*M123 + 2*
      Lambda7*vS + 2*vS*Conj(Lambda7) + 1.4142135623730951*Conj(M123))*ZA(gI2,1) +
      2*(2*Lambda5*v1 - Lambda7*v2 - v2*Conj(Lambda7))*ZA(gI2,2)) + ZA(gI1,1)*(4*
      Lambda3*v1*ZA(gI2,1) + (1.4142135623730951*M123 + 2*Lambda7*vS + 2*vS*Conj(
      Lambda7) + 1.4142135623730951*Conj(M123))*ZA(gI2,2))) - KroneckerDelta(2,gO2
      )*(ZA(gI1,1)*((1.4142135623730951*M123 - 2*Lambda7*vS - 2*vS*Conj(Lambda7) +
      1.4142135623730951*Conj(M123))*ZA(gI2,0) + 4*Lambda6*vS*ZA(gI2,1) + 2*v1*(
      Lambda7 + Conj(Lambda7))*ZA(gI2,2)) + ZA(gI1,0)*(4*Lambda5*vS*ZA(gI2,0) + (
      1.4142135623730951*M123 - 2*Lambda7*vS - 2*vS*Conj(Lambda7) +
      1.4142135623730951*Conj(M123))*ZA(gI2,1) + 2*v2*(Lambda7 + Conj(Lambda7))*ZA
      (gI2,2)) + 2*ZA(gI1,2)*(v2*(Lambda7 + Conj(Lambda7))*ZA(gI2,0) + v1*(Lambda7
       + Conj(Lambda7))*ZA(gI2,1) + (1.4142135623730951*M5 + 4*Lambda8*vS +
      1.4142135623730951*Conj(M5))*ZA(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpAhhhUhh(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.25)*(
      KroneckerDelta(1,gO2)*((1.4142135623730951*M123 - 2*Lambda7*vS + 2*vS*Conj(
      Lambda7) - 1.4142135623730951*Conj(M123))*ZA(gI2,0)*ZH(gI1,2) + ZA(gI2,2)*((
      1.4142135623730951*M123 + 2*Lambda7*vS - 2*vS*Conj(Lambda7) -
      1.4142135623730951*Conj(M123))*ZH(gI1,0) + 2*v1*(Lambda7 - Conj(Lambda7))*ZH
      (gI1,2))) + KroneckerDelta(0,gO2)*((1.4142135623730951*M123 - 2*Lambda7*vS +
      2*vS*Conj(Lambda7) - 1.4142135623730951*Conj(M123))*ZA(gI2,1)*ZH(gI1,2) + ZA
      (gI2,2)*((1.4142135623730951*M123 + 2*Lambda7*vS - 2*vS*Conj(Lambda7) -
      1.4142135623730951*Conj(M123))*ZH(gI1,1) + 2*v2*(Lambda7 - Conj(Lambda7))*ZH
      (gI1,2))) + KroneckerDelta(2,gO2)*(ZA(gI2,1)*((1.4142135623730951*M123 - 2*
      Lambda7*vS + 2*vS*Conj(Lambda7) - 1.4142135623730951*Conj(M123))*ZH(gI1,0) +
      2*v1*(-Lambda7 + Conj(Lambda7))*ZH(gI1,2)) + ZA(gI2,0)*((1.4142135623730951*
      M123 - 2*Lambda7*vS + 2*vS*Conj(Lambda7) - 1.4142135623730951*Conj(M123))*ZH
      (gI1,1) + 2*v2*(-Lambda7 + Conj(Lambda7))*ZH(gI1,2)) + 2*ZA(gI2,2)*(v2*(
      Lambda7 - Conj(Lambda7))*ZH(gI1,0) + v1*(Lambda7 - Conj(Lambda7))*ZH(gI1,1)
      + 1.4142135623730951*(M5 - Conj(M5))*ZH(gI1,2))));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUhh(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.25*(-(KroneckerDelta(1,gO2)*(4*ZH(gI1,1)*
      (Lambda3*v1*ZH(gI2,0) + 3*Lambda2*v2*ZH(gI2,1) + Lambda6*vS*ZH(gI2,2)) + ZH(
      gI1,2)*(-((1.4142135623730951*M123 - 2*Lambda7*vS - 2*vS*Conj(Lambda7) +
      1.4142135623730951*Conj(M123))*ZH(gI2,0)) + 4*Lambda6*vS*ZH(gI2,1) + 2*(
      Lambda7*v1 + 2*Lambda6*v2 + v1*Conj(Lambda7))*ZH(gI2,2)) + ZH(gI1,0)*(4*
      Lambda3*v2*ZH(gI2,0) + 4*Lambda3*v1*ZH(gI2,1) - (1.4142135623730951*M123 - 2
      *Lambda7*vS - 2*vS*Conj(Lambda7) + 1.4142135623730951*Conj(M123))*ZH(gI2,2))
      )) - KroneckerDelta(0,gO2)*(4*ZH(gI1,0)*(3*Lambda1*v1*ZH(gI2,0) + Lambda3*v2
      *ZH(gI2,1) + Lambda5*vS*ZH(gI2,2)) + ZH(gI1,2)*(4*Lambda5*vS*ZH(gI2,0) - (
      1.4142135623730951*M123 - 2*Lambda7*vS - 2*vS*Conj(Lambda7) +
      1.4142135623730951*Conj(M123))*ZH(gI2,1) + 2*(2*Lambda5*v1 + Lambda7*v2 + v2
      *Conj(Lambda7))*ZH(gI2,2)) + ZH(gI1,1)*(4*Lambda3*v2*ZH(gI2,0) + 4*Lambda3*
      v1*ZH(gI2,1) - (1.4142135623730951*M123 - 2*Lambda7*vS - 2*vS*Conj(Lambda7)
      + 1.4142135623730951*Conj(M123))*ZH(gI2,2))) - KroneckerDelta(2,gO2)*(ZH(gI1
      ,1)*(-((1.4142135623730951*M123 - 2*Lambda7*vS - 2*vS*Conj(Lambda7) +
      1.4142135623730951*Conj(M123))*ZH(gI2,0)) + 4*Lambda6*vS*ZH(gI2,1) + 2*(
      Lambda7*v1 + 2*Lambda6*v2 + v1*Conj(Lambda7))*ZH(gI2,2)) + ZH(gI1,0)*(4*
      Lambda5*vS*ZH(gI2,0) - (1.4142135623730951*M123 - 2*Lambda7*vS - 2*vS*Conj(
      Lambda7) + 1.4142135623730951*Conj(M123))*ZH(gI2,1) + 2*(2*Lambda5*v1 +
      Lambda7*v2 + v2*Conj(Lambda7))*ZH(gI2,2)) + 2*ZH(gI1,2)*((2*Lambda5*v1 +
      Lambda7*v2 + v2*Conj(Lambda7))*ZH(gI2,0) + (Lambda7*v1 + 2*Lambda6*v2 + v1*
      Conj(Lambda7))*ZH(gI2,1) - (1.4142135623730951*M5 - 12*Lambda8*vS +
      1.4142135623730951*Conj(M5))*ZH(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdUhhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.7071067811865475*KroneckerDelta(0,gO2)*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gI2,j1))*Vd(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdUhhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = 0.7071067811865475*KroneckerDelta(0,gO1)*
      SUM(j2,0,2,Conj(Vd(gI2,j2))*SUM(j1,0,2,Conj(Ud(gI1,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeUhhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.7071067811865475*KroneckerDelta(0,gO2)*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gI2,j1))*Ve(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeUhhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = 0.7071067811865475*KroneckerDelta(0,gO1)*
      SUM(j2,0,2,Conj(Ve(gI2,j2))*SUM(j1,0,2,Conj(Ue(gI1,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuUhhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.7071067811865475*KroneckerDelta(1,gO2)*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gI2,j1))*Vu(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuUhhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = 0.7071067811865475*KroneckerDelta(1,gO1)*
      SUM(j2,0,2,Conj(Vu(gI2,j2))*SUM(j1,0,2,Conj(Uu(gI1,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhHmconjVWm(int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.5*(-(g2*KroneckerDelta(0,gO2)*ZP(gI2,0))
      + g2*KroneckerDelta(1,gO2)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpAhUhhVZ(int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.1)*(5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*(KroneckerDelta(0,gO2)*ZA(
      gI2,0) - KroneckerDelta(1,gO2)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbargWmgWmUAh(int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.25)*(v1*
      KroneckerDelta(0,gO1) - v2*KroneckerDelta(1,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbargWmCgWmCUAh(int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.25)*(v1*
      KroneckerDelta(0,gO1) - v2*KroneckerDelta(1,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjVWmVWm(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0
      ,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhVZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*(7.745966692414834*g1*
      g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5
      *Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhHmconjHm(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*(-2*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*(ZP(gI1,1)*(Conj(Lambda7)*ZP(gI2,0) + Lambda6*ZP(gI2,1
      )) + ZP(gI1,0)*(Lambda5*ZP(gI2,0) + Lambda7*ZP(gI2,1))) - KroneckerDelta(1,
      gO1)*(Lambda4*KroneckerDelta(0,gO2)*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,
      1)) + 2*KroneckerDelta(1,gO2)*((Lambda3 - Lambda4)*ZP(gI1,0)*ZP(gI2,0) +
      Lambda2*ZP(gI1,1)*ZP(gI2,1))) - KroneckerDelta(0,gO1)*(Lambda4*
      KroneckerDelta(1,gO2)*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) + 2*
      KroneckerDelta(0,gO2)*(Lambda1*ZP(gI1,0)*ZP(gI2,0) + (Lambda3 - Lambda4)*ZP(
      gI1,1)*ZP(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhHmconjHm(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.5)*(Lambda4*v2*
      KroneckerDelta(0,gO2)*(ZP(gI1,1)*ZP(gI2,0) - ZP(gI1,0)*ZP(gI2,1)) + Lambda4*
      v1*KroneckerDelta(1,gO2)*(ZP(gI1,1)*ZP(gI2,0) - ZP(gI1,0)*ZP(gI2,1)) +
      KroneckerDelta(2,gO2)*((1.4142135623730951*M123 + 2*Lambda7*vS)*ZP(gI1,1)*ZP
      (gI2,0) - (2*vS*Conj(Lambda7) + 1.4142135623730951*Conj(M123))*ZP(gI1,0)*ZP(
      gI2,1)));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUAhUAh(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.5*(-(KroneckerDelta(2,gO1)*(
      KroneckerDelta(1,gO2)*(ZA(gI1,2)*((Lambda7 + Conj(Lambda7))*ZA(gI2,0) + 2*
      Lambda6*ZA(gI2,1)) + ((Lambda7 + Conj(Lambda7))*ZA(gI1,0) + 2*Lambda6*ZA(gI1
      ,1))*ZA(gI2,2)) + KroneckerDelta(0,gO2)*(ZA(gI1,2)*(2*Lambda5*ZA(gI2,0) + (
      Lambda7 + Conj(Lambda7))*ZA(gI2,1)) + (2*Lambda5*ZA(gI1,0) + (Lambda7 + Conj
      (Lambda7))*ZA(gI1,1))*ZA(gI2,2)) + KroneckerDelta(2,gO2)*(ZA(gI1,1)*((
      Lambda7 + Conj(Lambda7))*ZA(gI2,0) + 2*Lambda6*ZA(gI2,1)) + ZA(gI1,0)*(2*
      Lambda5*ZA(gI2,0) + (Lambda7 + Conj(Lambda7))*ZA(gI2,1)) + 12*Lambda8*ZA(gI1
      ,2)*ZA(gI2,2)))) - KroneckerDelta(1,gO1)*(KroneckerDelta(2,gO2)*(ZA(gI1,2)*(
      (Lambda7 + Conj(Lambda7))*ZA(gI2,0) + 2*Lambda6*ZA(gI2,1)) + ((Lambda7 +
      Conj(Lambda7))*ZA(gI1,0) + 2*Lambda6*ZA(gI1,1))*ZA(gI2,2)) + 2*
      KroneckerDelta(1,gO2)*(Lambda3*ZA(gI1,0)*ZA(gI2,0) + 3*Lambda2*ZA(gI1,1)*ZA(
      gI2,1) + Lambda6*ZA(gI1,2)*ZA(gI2,2)) + KroneckerDelta(0,gO2)*(2*Lambda3*ZA(
      gI1,1)*ZA(gI2,0) + 2*Lambda3*ZA(gI1,0)*ZA(gI2,1) + (Lambda7 + Conj(Lambda7))
      *ZA(gI1,2)*ZA(gI2,2))) - KroneckerDelta(0,gO1)*(KroneckerDelta(2,gO2)*(ZA(
      gI1,2)*(2*Lambda5*ZA(gI2,0) + (Lambda7 + Conj(Lambda7))*ZA(gI2,1)) + (2*
      Lambda5*ZA(gI1,0) + (Lambda7 + Conj(Lambda7))*ZA(gI1,1))*ZA(gI2,2)) + 2*
      KroneckerDelta(0,gO2)*(3*Lambda1*ZA(gI1,0)*ZA(gI2,0) + Lambda3*ZA(gI1,1)*ZA(
      gI2,1) + Lambda5*ZA(gI1,2)*ZA(gI2,2)) + KroneckerDelta(1,gO2)*(2*Lambda3*ZA(
      gI1,1)*ZA(gI2,0) + 2*Lambda3*ZA(gI1,0)*ZA(gI2,1) + (Lambda7 + Conj(Lambda7))
      *ZA(gI1,2)*ZA(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhhhhh(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*(-(KroneckerDelta(1,gO1)*(2*
      KroneckerDelta(1,gO2)*(Lambda3*ZH(gI1,0)*ZH(gI2,0) + Lambda2*ZH(gI1,1)*ZH(
      gI2,1) + Lambda6*ZH(gI1,2)*ZH(gI2,2)) + (Lambda7 + Conj(Lambda7))*(-(
      KroneckerDelta(0,gO2)*ZH(gI1,2)*ZH(gI2,2)) + KroneckerDelta(2,gO2)*(ZH(gI1,2
      )*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2))))) - KroneckerDelta(2,gO1)*(
      KroneckerDelta(2,gO2)*(-(ZH(gI1,1)*((Lambda7 + Conj(Lambda7))*ZH(gI2,0) - 2*
      Lambda6*ZH(gI2,1))) + ZH(gI1,0)*(2*Lambda5*ZH(gI2,0) - (Lambda7 + Conj(
      Lambda7))*ZH(gI2,1)) + 4*Lambda8*ZH(gI1,2)*ZH(gI2,2)) + (Lambda7 + Conj(
      Lambda7))*(KroneckerDelta(1,gO2)*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2))
      + KroneckerDelta(0,gO2)*(ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2)))) -
      KroneckerDelta(0,gO1)*(2*KroneckerDelta(0,gO2)*(Lambda1*ZH(gI1,0)*ZH(gI2,0)
      + Lambda3*ZH(gI1,1)*ZH(gI2,1) + Lambda5*ZH(gI1,2)*ZH(gI2,2)) + (Lambda7 +
      Conj(Lambda7))*(-(KroneckerDelta(1,gO2)*ZH(gI1,2)*ZH(gI2,2)) +
      KroneckerDelta(2,gO2)*(ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2)))));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUAh(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.25)*(
      KroneckerDelta(1,gO2)*((1.4142135623730951*M123 + 2*Lambda7*vS - 2*vS*Conj(
      Lambda7) - 1.4142135623730951*Conj(M123))*ZA(gI1,0)*ZA(gI2,2) + ZA(gI1,2)*((
      1.4142135623730951*M123 + 2*Lambda7*vS - 2*vS*Conj(Lambda7) -
      1.4142135623730951*Conj(M123))*ZA(gI2,0) + 2*v1*(-Lambda7 + Conj(Lambda7))*
      ZA(gI2,2))) + KroneckerDelta(0,gO2)*((1.4142135623730951*M123 + 2*Lambda7*vS
       - 2*vS*Conj(Lambda7) - 1.4142135623730951*Conj(M123))*ZA(gI1,1)*ZA(gI2,2) +
      ZA(gI1,2)*((1.4142135623730951*M123 + 2*Lambda7*vS - 2*vS*Conj(Lambda7) -
      1.4142135623730951*Conj(M123))*ZA(gI2,1) + 2*v2*(-Lambda7 + Conj(Lambda7))*
      ZA(gI2,2))) + KroneckerDelta(2,gO2)*(ZA(gI1,1)*((1.4142135623730951*M123 + 2
      *Lambda7*vS - 2*vS*Conj(Lambda7) - 1.4142135623730951*Conj(M123))*ZA(gI2,0)
      + 2*v1*(-Lambda7 + Conj(Lambda7))*ZA(gI2,2)) + ZA(gI1,0)*((
      1.4142135623730951*M123 + 2*Lambda7*vS - 2*vS*Conj(Lambda7) -
      1.4142135623730951*Conj(M123))*ZA(gI2,1) + 2*v2*(-Lambda7 + Conj(Lambda7))*
      ZA(gI2,2)) - 2*ZA(gI1,2)*(v2*(Lambda7 - Conj(Lambda7))*ZA(gI2,0) + v1*(
      Lambda7 - Conj(Lambda7))*ZA(gI2,1) + 1.4142135623730951*(-M5 + Conj(M5))*ZA(
      gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpAhUAhhh(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = 0.25*(-(KroneckerDelta(1,gO2)*((
      1.4142135623730951*M123 - 2*Lambda7*vS - 2*vS*Conj(Lambda7) +
      1.4142135623730951*Conj(M123))*ZA(gI2,0)*ZH(gI1,2) + 4*ZA(gI2,1)*(Lambda3*v1
      *ZH(gI1,0) + Lambda2*v2*ZH(gI1,1) + Lambda6*vS*ZH(gI1,2)) + ZA(gI2,2)*((
      1.4142135623730951*M123 + 2*Lambda7*vS + 2*vS*Conj(Lambda7) +
      1.4142135623730951*Conj(M123))*ZH(gI1,0) + 2*v1*(Lambda7 + Conj(Lambda7))*ZH
      (gI1,2)))) - KroneckerDelta(0,gO2)*((1.4142135623730951*M123 - 2*Lambda7*vS
      - 2*vS*Conj(Lambda7) + 1.4142135623730951*Conj(M123))*ZA(gI2,1)*ZH(gI1,2) +
      4*ZA(gI2,0)*(Lambda1*v1*ZH(gI1,0) + Lambda3*v2*ZH(gI1,1) + Lambda5*vS*ZH(gI1
      ,2)) + ZA(gI2,2)*((1.4142135623730951*M123 + 2*Lambda7*vS + 2*vS*Conj(
      Lambda7) + 1.4142135623730951*Conj(M123))*ZH(gI1,1) + 2*v2*(Lambda7 + Conj(
      Lambda7))*ZH(gI1,2))) - KroneckerDelta(2,gO2)*(ZA(gI2,1)*((
      1.4142135623730951*M123 + 2*Lambda7*vS + 2*vS*Conj(Lambda7) +
      1.4142135623730951*Conj(M123))*ZH(gI1,0) + 2*v1*(Lambda7 + Conj(Lambda7))*ZH
      (gI1,2)) + ZA(gI2,0)*((1.4142135623730951*M123 + 2*Lambda7*vS + 2*vS*Conj(
      Lambda7) + 1.4142135623730951*Conj(M123))*ZH(gI1,1) + 2*v2*(Lambda7 + Conj(
      Lambda7))*ZH(gI1,2)) + 2*ZA(gI2,2)*((2*Lambda5*v1 - Lambda7*v2 - v2*Conj(
      Lambda7))*ZH(gI1,0) - (Lambda7*v1 - 2*Lambda6*v2 + v1*Conj(Lambda7))*ZH(gI1,
      1) + (1.4142135623730951*M5 + 4*Lambda8*vS + 1.4142135623730951*Conj(M5))*ZH
      (gI1,2))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhhh(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.25)*(
      KroneckerDelta(1,gO2)*((1.4142135623730951*M123 - 2*Lambda7*vS + 2*vS*Conj(
      Lambda7) - 1.4142135623730951*Conj(M123))*ZH(gI1,0)*ZH(gI2,2) + ZH(gI1,2)*((
      1.4142135623730951*M123 - 2*Lambda7*vS + 2*vS*Conj(Lambda7) -
      1.4142135623730951*Conj(M123))*ZH(gI2,0) + 2*v1*(-Lambda7 + Conj(Lambda7))*
      ZH(gI2,2))) + KroneckerDelta(0,gO2)*((1.4142135623730951*M123 - 2*Lambda7*vS
       + 2*vS*Conj(Lambda7) - 1.4142135623730951*Conj(M123))*ZH(gI1,1)*ZH(gI2,2) +
      ZH(gI1,2)*((1.4142135623730951*M123 - 2*Lambda7*vS + 2*vS*Conj(Lambda7) -
      1.4142135623730951*Conj(M123))*ZH(gI2,1) + 2*v2*(-Lambda7 + Conj(Lambda7))*
      ZH(gI2,2))) + KroneckerDelta(2,gO2)*(ZH(gI1,1)*((1.4142135623730951*M123 + 2
      *Lambda7*vS - 2*vS*Conj(Lambda7) - 1.4142135623730951*Conj(M123))*ZH(gI2,0)
      + 2*v1*(Lambda7 - Conj(Lambda7))*ZH(gI2,2)) + ZH(gI1,0)*((1.4142135623730951
      *M123 + 2*Lambda7*vS - 2*vS*Conj(Lambda7) - 1.4142135623730951*Conj(M123))*
      ZH(gI2,1) + 2*v2*(Lambda7 - Conj(Lambda7))*ZH(gI2,2)) + 2*ZH(gI1,2)*(v2*(
      Lambda7 - Conj(Lambda7))*ZH(gI2,0) + v1*(Lambda7 - Conj(Lambda7))*ZH(gI2,1)
      + 1.4142135623730951*(M5 - Conj(M5))*ZH(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdUAhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*KroneckerDelta(0,gO2)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gI2,j1))*Vd
      (gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdUAhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *KroneckerDelta(0,gO1)*SUM(j2,0,2,Conj(Vd(gI2,j2))*SUM(j1,0,2,Conj(Ud(gI1,j1
      ))*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeUAhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*KroneckerDelta(0,gO2)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gI2,j1))*Ve
      (gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeUAhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *KroneckerDelta(0,gO1)*SUM(j2,0,2,Conj(Ve(gI2,j2))*SUM(j1,0,2,Conj(Ue(gI1,j1
      ))*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuUAhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*KroneckerDelta(1,gO2)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gI2,j1))*Vu
      (gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuUAhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *KroneckerDelta(1,gO1)*SUM(j2,0,2,Conj(Vu(gI2,j2))*SUM(j1,0,2,Conj(Uu(gI1,j1
      ))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhHmconjVWm(int gO2, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*g2*(
      KroneckerDelta(0,gO2)*ZP(gI2,0) + KroneckerDelta(1,gO2)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhVZ(int gO2, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.1)*(5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*(KroneckerDelta(0,gO2)*ZH(
      gI2,0) - KroneckerDelta(1,gO2)*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbargWmgZUHm(int gO2) const
{
   
   const std::complex<double> result = 0.05*g2*(v1*KroneckerDelta(0,gO2) - v2*
      KroneckerDelta(1,gO2))*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW
      ()));

   return result;
}

std::complex<double> CLASSNAME::CpbargZgWmconjUHm(int gO1) const
{
   
   const std::complex<double> result = -0.05*g2*(v1*KroneckerDelta(0,gO1) - v2*
      KroneckerDelta(1,gO1))*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW
      ()));

   return result;
}

std::complex<double> CLASSNAME::CpbargWmCgZconjUHm(int gO1) const
{
   
   const std::complex<double> result = 0.05*g2*(v1*KroneckerDelta(0,gO1) - v2*
      KroneckerDelta(1,gO1))*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW
      ()));

   return result;
}

std::complex<double> CLASSNAME::CpbargZgWmCUHm(int gO2) const
{
   
   const std::complex<double> result = -0.05*g2*(v1*KroneckerDelta(0,gO2) - v2*
      KroneckerDelta(1,gO2))*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW
      ()));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHmVPVWm(int gO2) const
{
   
   const std::complex<double> result = -0.3872983346207417*g1*g2*Cos(ThetaW())*(v1
      *KroneckerDelta(0,gO2) - v2*KroneckerDelta(1,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHmVWmVZ(int gO2) const
{
   
   const std::complex<double> result = 0.3872983346207417*g1*g2*(v1*KroneckerDelta
      (0,gO2) - v2*KroneckerDelta(1,gO2))*Sin(ThetaW());

   return result;
}

std::complex<double> CLASSNAME::CpUHmconjUHmconjVWmVWm(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0
      ,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUHmconjUHmVZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*(-7.745966692414834*g1
      *g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) +
      5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpHmUHmconjHmconjUHm(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -(KroneckerDelta(1,gO1)*(Lambda3*
      KroneckerDelta(0,gO2)*ZP(gI1,0)*ZP(gI2,1) + KroneckerDelta(1,gO2)*(Lambda3*
      ZP(gI1,0)*ZP(gI2,0) + 2*Lambda2*ZP(gI1,1)*ZP(gI2,1)))) - KroneckerDelta(0,
      gO1)*(Lambda3*KroneckerDelta(1,gO2)*ZP(gI1,1)*ZP(gI2,0) + KroneckerDelta(0,
      gO2)*(2*Lambda1*ZP(gI1,0)*ZP(gI2,0) + Lambda3*ZP(gI1,1)*ZP(gI2,1)));

   return result;
}

std::complex<double> CLASSNAME::CpAhHmconjUHm(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.5)*(
      KroneckerDelta(1,gO2)*(Lambda4*v2*ZA(gI2,0) + Lambda4*v1*ZA(gI2,1) + (
      1.4142135623730951*M123 + 2*Lambda7*vS)*ZA(gI2,2))*ZP(gI1,0) -
      KroneckerDelta(0,gO2)*(Lambda4*v2*ZA(gI2,0) + Lambda4*v1*ZA(gI2,1) + (2*vS*
      Conj(Lambda7) + 1.4142135623730951*Conj(M123))*ZA(gI2,2))*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CphhHmconjUHm(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = 0.5*(KroneckerDelta(1,gO2)*(ZH(gI2,0)*(
      Lambda4*v2*ZP(gI1,0) + 2*(-Lambda3 + Lambda4)*v1*ZP(gI1,1)) + ZH(gI2,1)*(
      Lambda4*v1*ZP(gI1,0) - 2*Lambda2*v2*ZP(gI1,1)) - ZH(gI2,2)*((
      1.4142135623730951*M123 - 2*Lambda7*vS)*ZP(gI1,0) + 2*Lambda6*vS*ZP(gI1,1)))
      + KroneckerDelta(0,gO2)*(ZH(gI2,1)*(2*(-Lambda3 + Lambda4)*v2*ZP(gI1,0) +
      Lambda4*v1*ZP(gI1,1)) + ZH(gI2,0)*(-2*Lambda1*v1*ZP(gI1,0) + Lambda4*v2*ZP(
      gI1,1)) - ZH(gI2,2)*(2*Lambda5*vS*ZP(gI1,0) + (-2*vS*Conj(Lambda7) +
      1.4142135623730951*Conj(M123))*ZP(gI1,1))));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUHmconjUHm(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.5*(-(KroneckerDelta(0,gO1)*(2*
      KroneckerDelta(0,gO2)*(Lambda1*ZA(gI1,0)*ZA(gI2,0) + (Lambda3 - Lambda4)*ZA(
      gI1,1)*ZA(gI2,1) + Lambda5*ZA(gI1,2)*ZA(gI2,2)) + KroneckerDelta(1,gO2)*(
      Lambda4*ZA(gI1,1)*ZA(gI2,0) + Lambda4*ZA(gI1,0)*ZA(gI2,1) + 2*Lambda7*ZA(gI1
      ,2)*ZA(gI2,2)))) - KroneckerDelta(1,gO1)*(2*KroneckerDelta(1,gO2)*((Lambda3
      - Lambda4)*ZA(gI1,0)*ZA(gI2,0) + Lambda2*ZA(gI1,1)*ZA(gI2,1) + Lambda6*ZA(
      gI1,2)*ZA(gI2,2)) + KroneckerDelta(0,gO2)*(Lambda4*ZA(gI1,1)*ZA(gI2,0) +
      Lambda4*ZA(gI1,0)*ZA(gI2,1) + 2*Conj(Lambda7)*ZA(gI1,2)*ZA(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUHmconjUHm(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.5*(KroneckerDelta(0,gO1)*(-2*
      KroneckerDelta(0,gO2)*(Lambda1*ZH(gI1,0)*ZH(gI2,0) + (Lambda3 - Lambda4)*ZH(
      gI1,1)*ZH(gI2,1) + Lambda5*ZH(gI1,2)*ZH(gI2,2)) + KroneckerDelta(1,gO2)*(
      Lambda4*ZH(gI1,1)*ZH(gI2,0) + Lambda4*ZH(gI1,0)*ZH(gI2,1) + 2*Lambda7*ZH(gI1
      ,2)*ZH(gI2,2))) + KroneckerDelta(1,gO1)*(-2*KroneckerDelta(1,gO2)*((Lambda3
      - Lambda4)*ZH(gI1,0)*ZH(gI2,0) + Lambda2*ZH(gI1,1)*ZH(gI2,1) + Lambda6*ZH(
      gI1,2)*ZH(gI2,2)) + KroneckerDelta(0,gO2)*(Lambda4*ZH(gI1,1)*ZH(gI2,0) +
      Lambda4*ZH(gI1,0)*ZH(gI2,1) + 2*Conj(Lambda7)*ZH(gI1,2)*ZH(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdconjUHmPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -(KroneckerDelta(0,gO2)*SUM(j2,0,2,SUM(j1,0
      ,2,Conj(Yd(j1,j2))*Ud(gI2,j1))*Vu(gI1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdconjUHmPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = -(KroneckerDelta(1,gO1)*SUM(j2,0,2,Conj(Vd(
      gI2,j2))*SUM(j1,0,2,Conj(Uu(gI1,j1))*Yu(j1,j2))));

   return result;
}

std::complex<double> CLASSNAME::CpbarFvFeconjUHmPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -(KroneckerDelta(0,gO2)*SUM(j1,0,2,Conj(Ye(
      j1,gI1))*Ue(gI2,j1)));

   return result;
}

double CLASSNAME::CpbarFvFeconjUHmPL(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpHmconjUHmVP(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 2,-0.3872983346207417*g1*Cos(
      ThetaW())*ZP(gI2,gO2),0) + IF(gI2 < 2,-0.5*g2*Sin(ThetaW())*ZP(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpHmconjUHmVZ(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 2,-0.5*g2*Cos(ThetaW())*ZP(gI2,gO2
      ),0) + IF(gI2 < 2,0.3872983346207417*g1*Sin(ThetaW())*ZP(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpAhconjUHmVWm(int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*g2*(
      KroneckerDelta(0,gO2)*ZA(gI2,0) + KroneckerDelta(1,gO2)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CphhconjUHmVWm(int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.5*g2*(KroneckerDelta(0,gO2)*ZH(gI2,0) -
      KroneckerDelta(1,gO2)*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVGVGVG() const
{
   
   const std::complex<double> result = std::complex<double>(0,-1)*g3;

   return result;
}

std::complex<double> CLASSNAME::CpbargGgGVG() const
{
   
   const std::complex<double> result = std::complex<double>(0,-1)*g3;

   return result;
}

double CLASSNAME::CpbarFdFdVGPL(int gI1, int gI2) const
{
   
   const double result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpbarFdFdVGPR(int gI1, int gI2) const
{
   
   const double result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpbarFuFuVGPL(int gI1, int gI2) const
{
   
   const double result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpbarFuFuVGPR(int gI1, int gI2) const
{
   
   const double result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpVGVGVGVG1() const
{
   
   const double result = -16*Sqr(g3);

   return result;
}

double CLASSNAME::CpVGVGVGVG2() const
{
   
   const double result = 0;

   return result;
}

double CLASSNAME::CpVGVGVGVG3() const
{
   
   const double result = 16*Sqr(g3);

   return result;
}

double CLASSNAME::CpbargWmgWmVP() const
{
   
   const double result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbargWmCgWmCVP() const
{
   
   const double result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWmVPVWm() const
{
   
   const double result = g2*Sin(ThetaW());

   return result;
}

std::complex<double> CLASSNAME::CpHmconjHmVPVP(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(7.745966692414834*g1*g2*Sin(2*ThetaW(
      )) + 3*Sqr(g1) + Cos(2*ThetaW())*(3*Sqr(g1) - 5*Sqr(g2)) + 5*Sqr(g2))*(ZP(
      gI1,0)*ZP(gI2,0) + ZP(gI1,1)*ZP(gI2,1));

   return result;
}

double CLASSNAME::CpHmconjHmVP(int gI2, int gI1) const
{
   
   const double result = -0.5*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos(
      ThetaW()) + g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFdFdVPPL(int gI1, int gI2) const
{
   
   const double result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(
      0.7745966692414834*g1*Cos(ThetaW()) - 3*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFdFdVPPR(int gI1, int gI2) const
{
   
   const double result = 0.2581988897471611*g1*Cos(ThetaW())*KroneckerDelta(gI1,
      gI2);

   return result;
}

double CLASSNAME::CpbarFeFeVPPL(int gI1, int gI2) const
{
   
   const double result = 0.5*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos(
      ThetaW()) + g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFeFeVPPR(int gI1, int gI2) const
{
   
   const double result = 0.7745966692414834*g1*Cos(ThetaW())*KroneckerDelta(gI1,
      gI2);

   return result;
}

double CLASSNAME::CpbarFuFuVPPL(int gI1, int gI2) const
{
   
   const double result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(
      0.7745966692414834*g1*Cos(ThetaW()) + 3*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFuFuVPPR(int gI1, int gI2) const
{
   
   const double result = -0.5163977794943222*g1*Cos(ThetaW())*KroneckerDelta(gI1,
      gI2);

   return result;
}

std::complex<double> CLASSNAME::CpHmconjVWmVP(int gI2) const
{
   
   const std::complex<double> result = -0.3872983346207417*g1*g2*Cos(ThetaW())*(v1
      *ZP(gI2,0) - v2*ZP(gI2,1));

   return result;
}

double CLASSNAME::CpconjVWmVPVPVWm3() const
{
   
   const double result = -2*Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWmVPVPVWm1() const
{
   
   const double result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWmVPVPVWm2() const
{
   
   const double result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbargWmgWmVZ() const
{
   
   const double result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpbargWmCgWmCVZ() const
{
   
   const double result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWmVWmVZ() const
{
   
   const double result = -(g2*Cos(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpHmconjHmVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(-7.745966692414834*g1*g2*Sin(2*ThetaW
      ()) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)))*(ZP(
      gI1,0)*ZP(gI2,0) + ZP(gI1,1)*ZP(gI2,1));

   return result;
}

double CLASSNAME::CpHmconjHmVZ(int gI2, int gI1) const
{
   
   const double result = 0.1*KroneckerDelta(gI1,gI2)*(-5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(7.745966692414834*g1*g2*Sin(2*ThetaW(
      )) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)))*(ZA(
      gI1,0)*ZA(gI2,0) + ZA(gI1,1)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CphhhhVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(7.745966692414834*g1*g2*Sin(2*ThetaW(
      )) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)))*(ZH(
      gI1,0)*ZH(gI2,0) + ZH(gI1,1)*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpAhhhVZ(int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*(g2*Cos(ThetaW(
      )) + 0.7745966692414834*g1*Sin(ThetaW()))*(ZA(gI2,0)*ZH(gI1,0) - ZA(gI2,1)*
      ZH(gI1,1));

   return result;
}

double CLASSNAME::CpbarFdFdVZPL(int gI1, int gI2) const
{
   
   const double result = 0.16666666666666666*KroneckerDelta(gI1,gI2)*(3*g2*Cos(
      ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFdFdVZPR(int gI1, int gI2) const
{
   
   const double result = -0.2581988897471611*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW
      ());

   return result;
}

double CLASSNAME::CpbarFeFeVZPL(int gI1, int gI2) const
{
   
   const double result = 0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) -
      0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFeFeVZPR(int gI1, int gI2) const
{
   
   const double result = -0.7745966692414834*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW
      ());

   return result;
}

double CLASSNAME::CpbarFuFuVZPL(int gI1, int gI2) const
{
   
   const double result = 0.03333333333333333*KroneckerDelta(gI1,gI2)*(-15*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFuFuVZPR(int gI1, int gI2) const
{
   
   const double result = 0.5163977794943222*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW(
      ));

   return result;
}

double CLASSNAME::CpbarFvFvVZPL(int gI1, int gI2) const
{
   
   const double result = -0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFvFvVZPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpHmconjVWmVZ(int gI2) const
{
   
   const std::complex<double> result = 0.3872983346207417*g1*g2*Sin(ThetaW())*(v1*
      ZP(gI2,0) - v2*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CphhVZVZ(int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()))*(v1*ZH(gI2,0) + v2*ZH(gI2,1));

   return result;
}

double CLASSNAME::CpconjVWmVWmVZVZ1() const
{
   
   const double result = -2*Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWmVWmVZVZ2() const
{
   
   const double result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWmVWmVZVZ3() const
{
   
   const double result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpbargPgWmconjVWm() const
{
   
   const double result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbargWmCgPconjVWm() const
{
   
   const double result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbargWmCgZconjVWm() const
{
   
   const double result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpbargZgWmconjVWm() const
{
   
   const double result = g2*Cos(ThetaW());

   return result;
}

std::complex<double> CLASSNAME::CpHmconjHmconjVWmVWm(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2)*(ZP(gI1,0)*ZP(gI2,0) + ZP(gI1,1
      )*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpAhHmconjVWm(int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*g2*(ZA(gI2,0)*
      ZP(gI1,0) + ZA(gI2,1)*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CphhHmconjVWm(int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.5*g2*(ZH(gI2,0)*ZP(gI1,0) - ZH(gI2,1)*ZP
      (gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhconjVWmVWm(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2)*(ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,1
      )*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CphhhhconjVWmVWm(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2)*(ZH(gI1,0)*ZH(gI2,0) + ZH(gI1,1
      )*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdconjVWmPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(Vd(
      gI2,j1))*Vu(gI1,j1));

   return result;
}

double CLASSNAME::CpbarFuFdconjVWmPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvFeconjVWmPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gI1 < 3,-0.7071067811865475*g2*Conj(Ve(
      gI2,gI1)),0);

   return result;
}

double CLASSNAME::CpbarFvFeconjVWmPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CphhconjVWmVWm(int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2)*(v1*ZH(gI2,0) + v2*ZH(gI2,1));

   return result;
}

double CLASSNAME::CpconjVWmconjVWmVWmVWm2() const
{
   
   const double result = -Sqr(g2);

   return result;
}

double CLASSNAME::CpconjVWmconjVWmVWmVWm1() const
{
   
   const double result = 2*Sqr(g2);

   return result;
}

double CLASSNAME::CpconjVWmconjVWmVWmVWm3() const
{
   
   const double result = -Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFuHmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-(SUM(j2,0,2,Conj(Vu(gI2,j2))*Yd
      (gO2,j2))*ZP(gI1,0)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFuHmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(SUM(j1,0,2,Conj(Yu(j1,gO1))*Uu
      (gI2,j1))*ZP(gI1,1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Vd(gI1,j2))*Yd(gO2,j2))*ZA(gI2,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,-
      0.7071067811865475)*SUM(j1,0,2,Conj(Yd(j1,gO1))*Ud(gI1,j1))*ZA(gI2,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*SUM(j2,0,2,
      Conj(Vd(gI2,j2))*Yd(gO2,j2))*ZH(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.7071067811865475*SUM(j1,0,2,
      Conj(Yd(j1,gO1))*Ud(gI2,j1))*ZH(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVGPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-(g3*Ud(gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVGPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-(g3*Conj(Vd(gI2,gO1))),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVPPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.2581988897471611*g1*Cos(ThetaW
      ())*Ud(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVPPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.12909944487358055*g1*Conj(Vd(
      gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 3,0.5*g2*Conj(Vd(gI2,gO1))*Sin(ThetaW(
      )),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVZPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.2581988897471611*g1*Sin(
      ThetaW())*Ud(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVZPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(Vd(gI2,gO1))*Cos(
      ThetaW()),0) + IF(gI2 < 3,0.12909944487358055*g1*Conj(Vd(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

double CLASSNAME::CpbarUFdFuVWmPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFuVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(Vu(
      gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFdconjHmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-(SUM(j2,0,2,Conj(Vd(gI2,j2))*Yu
      (gO2,j2))*ZP(gI1,1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFdconjHmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(SUM(j1,0,2,Conj(Yd(j1,gO1))*Ud
      (gI2,j1))*ZP(gI1,0)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Vu(gI1,j2))*Yu(gO2,j2))*ZA(gI2,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,-
      0.7071067811865475)*SUM(j1,0,2,Conj(Yu(j1,gO1))*Uu(gI1,j1))*ZA(gI2,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*SUM(j2,0,2,
      Conj(Vu(gI2,j2))*Yu(gO2,j2))*ZH(gI1,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.7071067811865475*SUM(j1,0,2,
      Conj(Yu(j1,gO1))*Uu(gI2,j1))*ZH(gI1,1),0);

   return result;
}

double CLASSNAME::CpbarUFuFdconjVWmPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFdconjVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(Vd(
      gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVGPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-(g3*Uu(gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVGPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-(g3*Conj(Vu(gI2,gO1))),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVPPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.5163977794943222*g1*Cos(
      ThetaW())*Uu(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVPPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.12909944487358055*g1*Conj(Vu(
      gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 3,-0.5*g2*Conj(Vu(gI2,gO1))*Sin(ThetaW
      ()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVZPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.5163977794943222*g1*Sin(ThetaW
      ())*Uu(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVZPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.5*g2*Conj(Vu(gI2,gO1))*Cos(
      ThetaW()),0) + IF(gI2 < 3,0.12909944487358055*g1*Conj(Vu(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFvHmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-(Ye(gO2,gI2)*ZP(gI1,0)),0);

   return result;
}

double CLASSNAME::CpbarUFeFvHmPR(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Ve(gI1,j2))*Ye(gO2,j2))*ZA(gI2,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,-
      0.7071067811865475)*SUM(j1,0,2,Conj(Ye(j1,gO1))*Ue(gI1,j1))*ZA(gI2,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFehhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*SUM(j2,0,2,
      Conj(Ve(gI2,j2))*Ye(gO2,j2))*ZH(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFehhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.7071067811865475*SUM(j1,0,2,
      Conj(Ye(j1,gO1))*Ue(gI2,j1))*ZH(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVPPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.7745966692414834*g1*Cos(ThetaW
      ())*Ue(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVPPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.3872983346207417*g1*Conj(Ve(
      gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 3,0.5*g2*Conj(Ve(gI2,gO1))*Sin(ThetaW(
      )),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVZPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.7745966692414834*g1*Sin(
      ThetaW())*Ue(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVZPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(Ve(gI2,gO1))*Cos(
      ThetaW()),0) + IF(gI2 < 3,-0.3872983346207417*g1*Conj(Ve(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

double CLASSNAME::CpbarUFeFvVWmPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

double CLASSNAME::CpbarUFeFvVWmPL(int gO1, int gI2) const
{
   
   const double result = IF(gI2 < 3,-0.7071067811865475*g2*KroneckerDelta(gI2,gO1)
      ,0);

   return result;
}

double CLASSNAME::CpbarFvFeconjHmPL(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvFeconjHmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(SUM(j1,0,2,Conj(Ye(j1,gO1))*Ue(gI2,j1))*
      ZP(gI1,0));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFuHmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(SUM(j2,0,2,Conj(Vu(gI2,j2))*SUM(j1,0,2,
      Conj(Ud(gO2,j1))*Yd(j1,j2)))*ZP(gI1,0));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFuHmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(
      gI2,j1))*Vd(gO1,j2))*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j2,0,2,Conj(Vd(gI1,j2))*SUM(j1,0,2,Conj(Ud(gO2,j1))*Yd(j1,j2)))*ZA(gI2,
      0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gI1,j1))*Vd(gO1,j2))*ZA(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.7071067811865475*SUM(j2,0,2,Conj(Vd(gI2,
      j2))*SUM(j1,0,2,Conj(Ud(gO2,j1))*Yd(j1,j2)))*ZH(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Yd(j1,j2))*Ud(gI2,j1))*Vd(gO1,j2))*ZH(gI1,0);

   return result;
}

double CLASSNAME::CpbarFdFuVWmPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFuVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(Vu(
      gI2,j1))*Vd(gO1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFvHmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(SUM(j1,0,2,Conj(Ue(gO2,j1))*Ye(j1,gI2))*
      ZP(gI1,0));

   return result;
}

double CLASSNAME::CpbarFeFvHmPR(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j2,0,2,Conj(Ve(gI1,j2))*SUM(j1,0,2,Conj(Ue(gO2,j1))*Ye(j1,j2)))*ZA(gI2,
      0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gI1,j1))*Ve(gO1,j2))*ZA(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFehhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.7071067811865475*SUM(j2,0,2,Conj(Ve(gI2,
      j2))*SUM(j1,0,2,Conj(Ue(gO2,j1))*Ye(j1,j2)))*ZH(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFehhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Ye(j1,j2))*Ue(gI2,j1))*Ve(gO1,j2))*ZH(gI1,0);

   return result;
}

double CLASSNAME::CpbarFeFvVWmPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFvVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.7071067811865475*g2*Ve(gO1,
      gI2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdconjHmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(SUM(j2,0,2,Conj(Vd(gI2,j2))*SUM(j1,0,2,
      Conj(Uu(gO2,j1))*Yu(j1,j2)))*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdconjHmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(
      gI2,j1))*Vu(gO1,j2))*ZP(gI1,0));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j2,0,2,Conj(Vu(gI1,j2))*SUM(j1,0,2,Conj(Uu(gO2,j1))*Yu(j1,j2)))*ZA(gI2,
      1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gI1,j1))*Vu(gO1,j2))*ZA(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.7071067811865475*SUM(j2,0,2,Conj(Vu(gI2,
      j2))*SUM(j1,0,2,Conj(Uu(gO2,j1))*Yu(j1,j2)))*ZH(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Yu(j1,j2))*Uu(gI2,j1))*Vu(gO1,j2))*ZH(gI1,1);

   return result;
}


std::complex<double> CLASSNAME::self_energy_hh_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -(B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*CpbargWmCgWmCUhh(gO1)*
      CpbargWmCgWmCUhh(gO2));
   result += -(B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*CpbargWmgWmUhh(gO1)*CpbargWmgWmUhh(
      gO2));
   result += -(B0(Sqr(p),Sqr(MVZ),Sqr(MVZ))*CpbargZgZUhh(gO1)*CpbargZgZUhh(gO2));
   result += 2*(-1 + 2*B0(Sqr(p),Sqr(MVWm),Sqr(MVWm)))*Conj(CpUhhconjVWmVWm(gO2))*
      CpUhhconjVWmVWm(gO1);
   result += (-1 + 2*B0(Sqr(p),Sqr(MVZ),Sqr(MVZ)))*Conj(CpUhhVZVZ(gO2))*CpUhhVZVZ(
      gO1);
   result += 2*CpUhhUhhconjVWmVWm(gO1,gO2)*(2*A0(Sqr(MVWm)) - Sqr(MVWm));
   result += CpUhhUhhVZVZ(gO1,gO2)*(2*A0(Sqr(MVZ)) - Sqr(MVZ));
   result += -SUM(gI1,0,1,A0(Sqr(MHm(gI1)))*CpUhhUhhHmconjHm(gO1,gO2,gI1,gI1));
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MHm(gI1)),Sqr(MHm(gI2)))*Conj(
      CpUhhHmconjHm(gO2,gI2,gI1))*CpUhhHmconjHm(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUhhUhh(gI1,gI1,gO1,gO2));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUhhUhh(gI1,gI1,gO1,gO2));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MAh(gI1)),Sqr(MAh(gI2)))*
      Conj(CpAhAhUhh(gI1,gI2,gO2))*CpAhAhUhh(gI1,gI2,gO1)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(MAh(gI2)))*Conj(
      CpAhhhUhh(gI2,gI1,gO2))*CpAhhhUhh(gI2,gI1,gO1)));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(Mhh(gI2)))*
      Conj(CphhhhUhh(gI1,gI2,gO2))*CphhhhUhh(gI1,gI2,gO1)));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFdFdUhhPL(gI1,gI2,gO2))*
      CpbarFdFdUhhPL(gI1,gI2,gO1) + Conj(CpbarFdFdUhhPR(gI1,gI2,gO2))*
      CpbarFdFdUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFeFeUhhPL(gI1,gI2,gO2))*
      CpbarFeFeUhhPL(gI1,gI2,gO1) + Conj(CpbarFeFeUhhPR(gI1,gI2,gO2))*
      CpbarFeFeUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFuFuUhhPL(gI1,gI2,gO2))*
      CpbarFuFuUhhPL(gI1,gI2,gO1) + Conj(CpbarFuFuUhhPR(gI1,gI2,gO2))*
      CpbarFuFuUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2)))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(
      gI2)))*(Conj(CpbarFdFdUhhPR(gI1,gI2,gO2))*CpbarFdFdUhhPL(gI1,gI2,gO1) + Conj
      (CpbarFdFdUhhPL(gI1,gI2,gO2))*CpbarFdFdUhhPR(gI1,gI2,gO1))*MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(
      gI2)))*(Conj(CpbarFeFeUhhPR(gI1,gI2,gO2))*CpbarFeFeUhhPL(gI1,gI2,gO1) + Conj
      (CpbarFeFeUhhPL(gI1,gI2,gO2))*CpbarFeFeUhhPR(gI1,gI2,gO1))*MFe(gI2)));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(
      gI2)))*(Conj(CpbarFuFuUhhPR(gI1,gI2,gO2))*CpbarFuFuUhhPL(gI1,gI2,gO1) + Conj
      (CpbarFuFuUhhPL(gI1,gI2,gO2))*CpbarFuFuUhhPR(gI1,gI2,gO1))*MFu(gI2)));
   result += 2*SUM(gI2,0,1,Conj(CpUhhHmconjVWm(gO2,gI2))*CpUhhHmconjVWm(gO1,gI2)*
      F0(Sqr(p),Sqr(MHm(gI2)),Sqr(MVWm)));
   result += SUM(gI2,0,2,Conj(CpAhUhhVZ(gI2,gO2))*CpAhUhhVZ(gI2,gO1)*F0(Sqr(p),Sqr
      (MAh(gI2)),Sqr(MVZ)));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_hh_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = i; k < 3; k++)
         self_energy(i, k) = self_energy_hh_1loop(p, i, k);

   Hermitianize(self_energy);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Ah_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -(B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*CpbargWmCgWmCUAh(gO1)*
      CpbargWmCgWmCUAh(gO2));
   result += -(B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*CpbargWmgWmUAh(gO1)*CpbargWmgWmUAh(
      gO2));
   result += 2*CpUAhUAhconjVWmVWm(gO1,gO2)*(2*A0(Sqr(MVWm)) - Sqr(MVWm));
   result += CpUAhUAhVZVZ(gO1,gO2)*(2*A0(Sqr(MVZ)) - Sqr(MVZ));
   result += -SUM(gI1,0,1,A0(Sqr(MHm(gI1)))*CpUAhUAhHmconjHm(gO1,gO2,gI1,gI1));
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MHm(gI1)),Sqr(MHm(gI2)))*Conj(
      CpUAhHmconjHm(gO2,gI2,gI1))*CpUAhHmconjHm(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUAhUAh(gI1,gI1,gO1,gO2));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CpUAhUAhhhhh(gO1,gO2,gI1,gI1));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MAh(gI1)),Sqr(MAh(gI2)))*
      Conj(CpAhAhUAh(gI1,gI2,gO2))*CpAhAhUAh(gI1,gI2,gO1)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(MAh(gI2)))*Conj(
      CpAhUAhhh(gI2,gO2,gI1))*CpAhUAhhh(gI2,gO1,gI1)));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(Mhh(gI2)))*
      Conj(CpUAhhhhh(gO2,gI1,gI2))*CpUAhhhhh(gO1,gI1,gI2)));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFdFdUAhPL(gI1,gI2,gO2))*
      CpbarFdFdUAhPL(gI1,gI2,gO1) + Conj(CpbarFdFdUAhPR(gI1,gI2,gO2))*
      CpbarFdFdUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFeFeUAhPL(gI1,gI2,gO2))*
      CpbarFeFeUAhPL(gI1,gI2,gO1) + Conj(CpbarFeFeUAhPR(gI1,gI2,gO2))*
      CpbarFeFeUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFuFuUAhPL(gI1,gI2,gO2))*
      CpbarFuFuUAhPL(gI1,gI2,gO1) + Conj(CpbarFuFuUAhPR(gI1,gI2,gO2))*
      CpbarFuFuUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2)))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(
      gI2)))*(Conj(CpbarFdFdUAhPR(gI1,gI2,gO2))*CpbarFdFdUAhPL(gI1,gI2,gO1) + Conj
      (CpbarFdFdUAhPL(gI1,gI2,gO2))*CpbarFdFdUAhPR(gI1,gI2,gO1))*MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(
      gI2)))*(Conj(CpbarFeFeUAhPR(gI1,gI2,gO2))*CpbarFeFeUAhPL(gI1,gI2,gO1) + Conj
      (CpbarFeFeUAhPL(gI1,gI2,gO2))*CpbarFeFeUAhPR(gI1,gI2,gO1))*MFe(gI2)));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(
      gI2)))*(Conj(CpbarFuFuUAhPR(gI1,gI2,gO2))*CpbarFuFuUAhPL(gI1,gI2,gO1) + Conj
      (CpbarFuFuUAhPL(gI1,gI2,gO2))*CpbarFuFuUAhPR(gI1,gI2,gO1))*MFu(gI2)));
   result += 2*SUM(gI2,0,1,Conj(CpUAhHmconjVWm(gO2,gI2))*CpUAhHmconjVWm(gO1,gI2)*
      F0(Sqr(p),Sqr(MHm(gI2)),Sqr(MVWm)));
   result += SUM(gI2,0,2,Conj(CpUAhhhVZ(gO2,gI2))*CpUAhhhVZ(gO1,gI2)*F0(Sqr(p),Sqr
      (Mhh(gI2)),Sqr(MVZ)));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Ah_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = i; k < 3; k++)
         self_energy(i, k) = self_energy_Ah_1loop(p, i, k);

   Hermitianize(self_energy);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Hm_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -(B0(Sqr(p),Sqr(MVWm),Sqr(MVZ))*CpbargWmgZUHm(gO2)*CpbargZgWmconjUHm(
      gO1));
   result += -(B0(Sqr(p),Sqr(MVZ),Sqr(MVWm))*CpbargWmCgZconjUHm(gO1)*
      CpbargZgWmCUHm(gO2));
   result += 2*(-1 + 2*B0(Sqr(p),0,Sqr(MVWm)))*Conj(CpconjUHmVPVWm(gO2))*
      CpconjUHmVPVWm(gO1);
   result += 2*(-1 + 2*B0(Sqr(p),Sqr(MVWm),Sqr(MVZ)))*Conj(CpconjUHmVWmVZ(gO2))*
      CpconjUHmVWmVZ(gO1);
   result += 2*CpUHmconjUHmconjVWmVWm(gO1,gO2)*(2*A0(Sqr(MVWm)) - Sqr(MVWm));
   result += CpUHmconjUHmVZVZ(gO1,gO2)*(2*A0(Sqr(MVZ)) - Sqr(MVZ));
   result += -SUM(gI1,0,1,A0(Sqr(MHm(gI1)))*CpHmUHmconjHmconjUHm(gI1,gO1,gI1,gO2))
      ;
   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MHm(gI1)),Sqr(MAh(gI2)))*Conj(
      CpAhHmconjUHm(gI2,gI1,gO2))*CpAhHmconjUHm(gI2,gI1,gO1)));
   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MHm(gI1)),Sqr(Mhh(gI2)))*Conj(
      CphhHmconjUHm(gI2,gI1,gO2))*CphhHmconjUHm(gI2,gI1,gO1)));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUHmconjUHm(gI1,gI1,gO1,gO2))
      ;
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUHmconjUHm(gI1,gI1,gO1,gO2))
      ;
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFuFdconjUHmPL(gI1,gI2,gO2))*
      CpbarFuFdconjUHmPL(gI1,gI2,gO1) + Conj(CpbarFuFdconjUHmPR(gI1,gI2,gO2))*
      CpbarFuFdconjUHmPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFd(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFvFeconjUHmPL(gI1,gI2,gO2))*
      CpbarFvFeconjUHmPL(gI1,gI2,gO1) + Conj(CpbarFvFeconjUHmPR(gI1,gI2,gO2))*
      CpbarFvFeconjUHmPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFe(gI2)))));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFd(
      gI2)))*(Conj(CpbarFuFdconjUHmPR(gI1,gI2,gO2))*CpbarFuFdconjUHmPL(gI1,gI2,gO1
      ) + Conj(CpbarFuFdconjUHmPL(gI1,gI2,gO2))*CpbarFuFdconjUHmPR(gI1,gI2,gO1))*
      MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFe(
      gI2)))*(Conj(CpbarFvFeconjUHmPR(gI1,gI2,gO2))*CpbarFvFeconjUHmPL(gI1,gI2,gO1
      ) + Conj(CpbarFvFeconjUHmPL(gI1,gI2,gO2))*CpbarFvFeconjUHmPR(gI1,gI2,gO1))*
      MFe(gI2)));
   result += SUM(gI2,0,1,Conj(CpHmconjUHmVP(gI2,gO2))*CpHmconjUHmVP(gI2,gO1)*F0(
      Sqr(p),Sqr(MHm(gI2)),0));
   result += SUM(gI2,0,1,Conj(CpHmconjUHmVZ(gI2,gO2))*CpHmconjUHmVZ(gI2,gO1)*F0(
      Sqr(p),Sqr(MHm(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,2,Conj(CpAhconjUHmVWm(gI2,gO2))*CpAhconjUHmVWm(gI2,gO1)*F0(
      Sqr(p),Sqr(MAh(gI2)),Sqr(MVWm)));
   result += SUM(gI2,0,2,Conj(CphhconjUHmVWm(gI2,gO2))*CphhconjUHmVWm(gI2,gO1)*F0(
      Sqr(p),Sqr(Mhh(gI2)),Sqr(MVWm)));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_Hm_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = i; k < 2; k++)
         self_energy(i, k) = self_energy_Hm_1loop(p, i, k);

   Hermitianize(self_energy);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_VG_1loop(double p ) const
{
   std::complex<double> result;

   result += 3*AbsSqr(CpbargGgGVG())*B00(Sqr(p),Sqr(MVG),Sqr(MVG));
   result += -(AbsSqr(CpVGVGVG())*(15*B00(Sqr(p),0,0) + Sqr(p) + 6*B0(Sqr(p),0,0)*
      Sqr(p)));
   result += 0;
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVGPL(gI1,gI2)) + AbsSqr(
      CpbarFdFdVGPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(CpbarFdFdVGPL(gI1,
      gI2))*CpbarFdFdVGPR(gI1,gI2))));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuVGPL(gI1,gI2)) + AbsSqr(
      CpbarFuFuVGPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFu(gI1)),Sqr(MFu(gI2)))*MFu(gI1)*MFu(gI2)*Re(Conj(CpbarFuFuVGPL(gI1,
      gI2))*CpbarFuFuVGPR(gI1,gI2))));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_VP_1loop(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpbargWmCgWmCVP())*B00(Sqr(p),Sqr(MVWm),Sqr(MVWm));
   result += AbsSqr(CpbargWmgWmVP())*B00(Sqr(p),Sqr(MVWm),Sqr(MVWm));
   result += -(A0(Sqr(MVWm))*(CpconjVWmVPVPVWm1() + CpconjVWmVPVPVWm2() + 4*
      CpconjVWmVPVPVWm3()));
   result += 2*CpconjVWmVPVPVWm3()*Sqr(MVWm);
   result += -0.6666666666666666*AbsSqr(CpconjVWmVPVWm())*(3*A0(Sqr(MVWm)) + 15*
      B00(Sqr(p),Sqr(MVWm),Sqr(MVWm)) - 6*Sqr(MVWm) + Sqr(p) + 3*B0(Sqr(p),Sqr(
      MVWm),Sqr(MVWm))*(Sqr(MVWm) + 2*Sqr(p)));
   result += SUM(gI1,0,1,A0(Sqr(MHm(gI1)))*CpHmconjHmVPVP(gI1,gI1));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpHmconjHmVP(gI2,gI1))*B00(Sqr(p),
      Sqr(MHm(gI1)),Sqr(MHm(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVPPL(gI1,gI2)) + AbsSqr(
      CpbarFdFdVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(CpbarFdFdVPPL(gI1,
      gI2))*CpbarFdFdVPPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFeVPPL(gI1,gI2)) + AbsSqr(
      CpbarFeFeVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFe(gI1)),Sqr(MFe(gI2)))*MFe(gI1)*MFe(gI2)*Re(Conj(CpbarFeFeVPPL(gI1,
      gI2))*CpbarFeFeVPPR(gI1,gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuVPPL(gI1,gI2)) + AbsSqr(
      CpbarFuFuVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFu(gI1)),Sqr(MFu(gI2)))*MFu(gI1)*MFu(gI2)*Re(Conj(CpbarFuFuVPPL(gI1,
      gI2))*CpbarFuFuVPPR(gI1,gI2))));
   result += 2*SUM(gI2,0,1,AbsSqr(CpHmconjVWmVP(gI2))*B0(Sqr(p),Sqr(MVWm),Sqr(MHm(
      gI2))));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_VZ_1loop(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpbargWmCgWmCVZ())*B00(Sqr(p),Sqr(MVWm),Sqr(MVWm));
   result += AbsSqr(CpbargWmgWmVZ())*B00(Sqr(p),Sqr(MVWm),Sqr(MVWm));
   result += -(A0(Sqr(MVWm))*(4*CpconjVWmVWmVZVZ1() + CpconjVWmVWmVZVZ2() +
      CpconjVWmVWmVZVZ3()));
   result += 2*CpconjVWmVWmVZVZ1()*Sqr(MVWm);
   result += -0.6666666666666666*AbsSqr(CpconjVWmVWmVZ())*(3*A0(Sqr(MVWm)) + 15*
      B00(Sqr(p),Sqr(MVWm),Sqr(MVWm)) - 6*Sqr(MVWm) + Sqr(p) + 3*B0(Sqr(p),Sqr(
      MVWm),Sqr(MVWm))*(Sqr(MVWm) + 2*Sqr(p)));
   result += SUM(gI1,0,1,A0(Sqr(MHm(gI1)))*CpHmconjHmVZVZ(gI1,gI1));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpHmconjHmVZ(gI2,gI1))*B00(Sqr(p),
      Sqr(MHm(gI1)),Sqr(MHm(gI2)))));
   result += 0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhVZVZ(gI1,gI1));
   result += 0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhVZVZ(gI1,gI1));
   result += -4*SUM(gI1,0,2,SUM(gI2,0,2,AbsSqr(CpAhhhVZ(gI2,gI1))*B00(Sqr(p),Sqr(
      MAh(gI2)),Sqr(Mhh(gI1)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVZPL(gI1,gI2)) + AbsSqr(
      CpbarFdFdVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(CpbarFdFdVZPL(gI1,
      gI2))*CpbarFdFdVZPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFeVZPL(gI1,gI2)) + AbsSqr(
      CpbarFeFeVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFe(gI1)),Sqr(MFe(gI2)))*MFe(gI1)*MFe(gI2)*Re(Conj(CpbarFeFeVZPL(gI1,
      gI2))*CpbarFeFeVZPR(gI1,gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuVZPL(gI1,gI2)) + AbsSqr(
      CpbarFuFuVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFu(gI1)),Sqr(MFu(gI2)))*MFu(gI1)*MFu(gI2)*Re(Conj(CpbarFuFuVZPL(gI1,
      gI2))*CpbarFuFuVZPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFvFvVZPL(gI1,gI2)) + AbsSqr(
      CpbarFvFvVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFv(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFv(gI1)),Sqr(MFv(gI2)))*MFv(gI1)*MFv(gI2)*Re(Conj(CpbarFvFvVZPL(gI1,
      gI2))*CpbarFvFvVZPR(gI1,gI2))));
   result += 2*SUM(gI2,0,1,AbsSqr(CpHmconjVWmVZ(gI2))*B0(Sqr(p),Sqr(MVWm),Sqr(MHm(
      gI2))));
   result += SUM(gI2,0,2,AbsSqr(CphhVZVZ(gI2))*B0(Sqr(p),Sqr(MVZ),Sqr(Mhh(gI2))));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_VWm_1loop(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpbargPgWmconjVWm())*B00(Sqr(p),Sqr(MVWm),Sqr(MVP));
   result += AbsSqr(CpbargWmCgPconjVWm())*B00(Sqr(p),Sqr(MVP),Sqr(MVWm));
   result += AbsSqr(CpbargWmCgZconjVWm())*B00(Sqr(p),Sqr(MVZ),Sqr(MVWm));
   result += AbsSqr(CpbargZgWmconjVWm())*B00(Sqr(p),Sqr(MVWm),Sqr(MVZ));
   result += -(A0(Sqr(MVWm))*(CpconjVWmconjVWmVWmVWm1() + 4*
      CpconjVWmconjVWmVWmVWm2() + CpconjVWmconjVWmVWmVWm3()));
   result += 0;
   result += 2*CpconjVWmconjVWmVWmVWm2()*Sqr(MVWm);
   result += -(AbsSqr(CpconjVWmVPVWm())*(A0(Sqr(MVWm)) + 10*B00(Sqr(p),Sqr(MVWm),0
      ) - 2*Sqr(MVWm) + 0.6666666666666666*Sqr(p) + B0(Sqr(p),Sqr(MVWm),0)*(Sqr(
      MVWm) + 4*Sqr(p))));
   result += -0.5*A0(Sqr(MVZ))*(4*CpconjVWmVWmVZVZ1() + CpconjVWmVWmVZVZ2() +
      CpconjVWmVWmVZVZ3()) + CpconjVWmVWmVZVZ1()*Sqr(MVZ);
   result += -(AbsSqr(CpconjVWmVWmVZ())*(A0(Sqr(MVWm)) + A0(Sqr(MVZ)) + 10*B00(Sqr
      (p),Sqr(MVZ),Sqr(MVWm)) - 2*(Sqr(MVWm) + Sqr(MVZ) - 0.3333333333333333*Sqr(p
      )) + B0(Sqr(p),Sqr(MVZ),Sqr(MVWm))*(Sqr(MVWm) + Sqr(MVZ) + 4*Sqr(p))));
   result += SUM(gI1,0,1,A0(Sqr(MHm(gI1)))*CpHmconjHmconjVWmVWm(gI1,gI1));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,2,AbsSqr(CpAhHmconjVWm(gI2,gI1))*B00(Sqr(p),
      Sqr(MAh(gI2)),Sqr(MHm(gI1)))));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,2,AbsSqr(CphhHmconjVWm(gI2,gI1))*B00(Sqr(p),
      Sqr(Mhh(gI2)),Sqr(MHm(gI1)))));
   result += 0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhconjVWmVWm(gI1,gI1));
   result += 0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhconjVWmVWm(gI1,gI1));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFdconjVWmPL(gI1,gI2)) +
      AbsSqr(CpbarFuFdconjVWmPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFd(gI2)))
      + 4*B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFd(gI2)))*MFd(gI2)*MFu(gI1)*Re(Conj(
      CpbarFuFdconjVWmPL(gI1,gI2))*CpbarFuFdconjVWmPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFvFeconjVWmPL(gI1,gI2)) + AbsSqr
      (CpbarFvFeconjVWmPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFe(gI2))) + 4*B0
      (Sqr(p),Sqr(MFv(gI1)),Sqr(MFe(gI2)))*MFe(gI2)*MFv(gI1)*Re(Conj(
      CpbarFvFeconjVWmPL(gI1,gI2))*CpbarFvFeconjVWmPR(gI1,gI2))));
   result += SUM(gI2,0,1,AbsSqr(CpHmconjVWmVP(gI2))*B0(Sqr(p),0,Sqr(MHm(gI2))));
   result += SUM(gI2,0,1,AbsSqr(CpHmconjVWmVZ(gI2))*B0(Sqr(p),Sqr(MVZ),Sqr(MHm(gI2
      ))));
   result += SUM(gI2,0,2,AbsSqr(CphhconjVWmVWm(gI2))*B0(Sqr(p),Sqr(MVWm),Sqr(Mhh(
      gI2))));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_Fd_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MHm(gI1)))*Conj(
      CpbarUFdFuHmPL(gO2,gI2,gI1))*CpbarUFdFuHmPR(gO1,gI2,gI1)*MFu(gI2)));
   result += SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarUFdFdAhPL(gO2,gI1,gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarUFdFdhhPL(gO2,gI2,gI1))*CpbarUFdFdhhPR(gO1,gI2,gI1)*MFd(gI2)));
   result += -5.333333333333333*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),0))*
      Conj(CpbarUFdFdVGPR(gO2,gI2))*CpbarUFdFdVGPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),0))*Conj(
      CpbarUFdFdVPPR(gO2,gI2))*CpbarUFdFdVPPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFdFdVZPR(gO2,gI2))*CpbarUFdFdVZPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWm)))*Conj(
      CpbarUFdFuVWmPR(gO2,gI2))*CpbarUFdFuVWmPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fd_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fd_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHm(gI1)))*
      Conj(CpbarUFdFuHmPR(gO2,gI2,gI1))*CpbarUFdFuHmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFdFdAhPR(gO2,gI1,gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFdFdhhPR(gO2,gI2,gI1))*CpbarUFdFdhhPR(gO1,gI2,gI1)));
   result += -1.3333333333333333*SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),0))*
      Conj(CpbarUFdFdVGPL(gO2,gI2))*CpbarUFdFdVGPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),0))*Conj(CpbarUFdFdVPPL(
      gO2,gI2))*CpbarUFdFdVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFdFdVZPL(gO2,gI2))*CpbarUFdFdVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWm)))*Conj(
      CpbarUFdFuVWmPL(gO2,gI2))*CpbarUFdFuVWmPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fd_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fd_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHm(gI1)))*
      Conj(CpbarUFdFuHmPL(gO2,gI2,gI1))*CpbarUFdFuHmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFdFdAhPL(gO2,gI1,gI2))*CpbarUFdFdAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFdFdhhPL(gO2,gI2,gI1))*CpbarUFdFdhhPL(gO1,gI2,gI1)));
   result += -1.3333333333333333*SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),0))*
      Conj(CpbarUFdFdVGPR(gO2,gI2))*CpbarUFdFdVGPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),0))*Conj(CpbarUFdFdVPPR(
      gO2,gI2))*CpbarUFdFdVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFdFdVZPR(gO2,gI2))*CpbarUFdFdVZPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWm)))*Conj(
      CpbarUFdFuVWmPR(gO2,gI2))*CpbarUFdFuVWmPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fd_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MHm(gI1)))*Conj(
      CpbarUFuFdconjHmPL(gO2,gI2,gI1))*CpbarUFuFdconjHmPR(gO1,gI2,gI1)*MFd(gI2)));
   result += SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarUFuFuAhPL(gO2,gI1,gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarUFuFuhhPL(gO2,gI2,gI1))*CpbarUFuFuhhPR(gO1,gI2,gI1)*MFu(gI2)));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm)))*Conj(
      CpbarUFuFdconjVWmPR(gO2,gI2))*CpbarUFuFdconjVWmPL(gO1,gI2)*MFd(gI2));
   result += -5.333333333333333*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),0))*
      Conj(CpbarUFuFuVGPR(gO2,gI2))*CpbarUFuFuVGPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),0))*Conj(
      CpbarUFuFuVPPR(gO2,gI2))*CpbarUFuFuVPPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFuFuVZPR(gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHm(gI1)))*
      Conj(CpbarUFuFdconjHmPR(gO2,gI2,gI1))*CpbarUFuFdconjHmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFuFuAhPR(gO2,gI1,gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFuFuhhPR(gO2,gI2,gI1))*CpbarUFuFuhhPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm)))*Conj(
      CpbarUFuFdconjVWmPL(gO2,gI2))*CpbarUFuFdconjVWmPL(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),0))*
      Conj(CpbarUFuFuVGPL(gO2,gI2))*CpbarUFuFuVGPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),0))*Conj(CpbarUFuFuVPPL(
      gO2,gI2))*CpbarUFuFuVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFuFuVZPL(gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHm(gI1)))*
      Conj(CpbarUFuFdconjHmPL(gO2,gI2,gI1))*CpbarUFuFdconjHmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFuFuAhPL(gO2,gI1,gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFuFuhhPL(gO2,gI2,gI1))*CpbarUFuFuhhPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm)))*Conj(
      CpbarUFuFdconjVWmPR(gO2,gI2))*CpbarUFuFdconjVWmPR(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),0))*
      Conj(CpbarUFuFuVGPR(gO2,gI2))*CpbarUFuFuVGPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),0))*Conj(CpbarUFuFuVPPR(
      gO2,gI2))*CpbarUFuFuVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFuFuVZPR(gO2,gI2))*CpbarUFuFuVZPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MHm(gI1)))*Conj(
      CpbarUFeFvHmPL(gO2,gI2,gI1))*CpbarUFeFvHmPR(gO1,gI2,gI1)*MFv(gI2)));
   result += SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarUFeFeAhPL(gO2,gI1,gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarUFeFehhPL(gO2,gI2,gI1))*CpbarUFeFehhPR(gO1,gI2,gI1)*MFe(gI2)));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFe(gI2)),0))*Conj(
      CpbarUFeFeVPPR(gO2,gI2))*CpbarUFeFeVPPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFeFeVZPR(gO2,gI2))*CpbarUFeFeVZPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWm)))*Conj(
      CpbarUFeFvVWmPR(gO2,gI2))*CpbarUFeFvVWmPL(gO1,gI2)*MFv(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fe_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHm(gI1)))*
      Conj(CpbarUFeFvHmPR(gO2,gI2,gI1))*CpbarUFeFvHmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFeFeAhPR(gO2,gI1,gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFeFehhPR(gO2,gI2,gI1))*CpbarUFeFehhPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),0))*Conj(CpbarUFeFeVPPL(
      gO2,gI2))*CpbarUFeFeVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFeFeVZPL(gO2,gI2))*CpbarUFeFeVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWm)))*Conj(
      CpbarUFeFvVWmPL(gO2,gI2))*CpbarUFeFvVWmPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fe_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHm(gI1)))*
      Conj(CpbarUFeFvHmPL(gO2,gI2,gI1))*CpbarUFeFvHmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFeFeAhPL(gO2,gI1,gI2))*CpbarUFeFeAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFeFehhPL(gO2,gI2,gI1))*CpbarUFeFehhPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),0))*Conj(CpbarUFeFeVPPR(
      gO2,gI2))*CpbarUFeFeVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFeFeVZPR(gO2,gI2))*CpbarUFeFeVZPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWm)))*Conj(
      CpbarUFeFvVWmPR(gO2,gI2))*CpbarUFeFvVWmPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fe_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fv_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MHm(gI1)))*Conj(
      CpbarFvFeconjHmPL(gO2,gI2,gI1))*CpbarFvFeconjHmPR(gO1,gI2,gI1)*MFe(gI2)));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWm)))*Conj(
      CpbarFvFeconjVWmPR(gO2,gI2))*CpbarFvFeconjVWmPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZ)))*Conj(
      CpbarFvFvVZPR(gO2,gI2))*CpbarFvFvVZPL(gO1,gI2)*MFv(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fv_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fv_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fv_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MHm(gI1)))*
      Conj(CpbarFvFeconjHmPR(gO2,gI2,gI1))*CpbarFvFeconjHmPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWm)))*Conj(
      CpbarFvFeconjVWmPL(gO2,gI2))*CpbarFvFeconjVWmPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZ)))*Conj(
      CpbarFvFvVZPL(gO2,gI2))*CpbarFvFvVZPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fv_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fv_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fv_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MHm(gI1)))*
      Conj(CpbarFvFeconjHmPL(gO2,gI2,gI1))*CpbarFvFeconjHmPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWm)))*Conj(
      CpbarFvFeconjVWmPR(gO2,gI2))*CpbarFvFeconjVWmPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZ)))*Conj(
      CpbarFvFvVZPR(gO2,gI2))*CpbarFvFvVZPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fv_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fv_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fd_1loop_1_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MHm(gI1)))*Conj(
      CpbarFdFuHmPL(gO2,gI2,gI1))*CpbarFdFuHmPR(gO1,gI2,gI1)*MFu(gI2)));
   result += SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarFdFdAhPL(gO2,gI1,gI2))*CpbarFdFdAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarFdFdhhPL(gO2,gI2,gI1))*CpbarFdFdhhPR(gO1,gI2,gI1)*MFd(gI2)));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ)))*Conj(
      CpbarFdFdVZPR(gO2,gI2))*CpbarFdFdVZPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWm)))*Conj(
      CpbarFdFuVWmPR(gO2,gI2))*CpbarFdFuVWmPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fd_1loop_1_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_1_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fd_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHm(gI1)))*
      Conj(CpbarFdFuHmPR(gO2,gI2,gI1))*CpbarFdFuHmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFdFdAhPR(gO2,gI1,gI2))*CpbarFdFdAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFdFdhhPR(gO2,gI2,gI1))*CpbarFdFdhhPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ)))*Conj(
      CpbarFdFdVZPL(gO2,gI2))*CpbarFdFdVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWm)))*Conj(
      CpbarFdFuVWmPL(gO2,gI2))*CpbarFdFuVWmPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fd_1loop_PR_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_PR_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fd_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHm(gI1)))*
      Conj(CpbarFdFuHmPL(gO2,gI2,gI1))*CpbarFdFuHmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFdFdAhPL(gO2,gI1,gI2))*CpbarFdFdAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFdFdhhPL(gO2,gI2,gI1))*CpbarFdFdhhPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ)))*Conj(
      CpbarFdFdVZPR(gO2,gI2))*CpbarFdFdVZPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWm)))*Conj(
      CpbarFdFuVWmPR(gO2,gI2))*CpbarFdFuVWmPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fd_1loop_PL_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_PL_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_1_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MHm(gI1)))*Conj(
      CpbarFeFvHmPL(gO2,gI2,gI1))*CpbarFeFvHmPR(gO1,gI2,gI1)*MFv(gI2)));
   result += SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarFeFeAhPL(gO2,gI1,gI2))*CpbarFeFeAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarFeFehhPL(gO2,gI2,gI1))*CpbarFeFehhPR(gO1,gI2,gI1)*MFe(gI2)));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ)))*Conj(
      CpbarFeFeVZPR(gO2,gI2))*CpbarFeFeVZPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWm)))*Conj(
      CpbarFeFvVWmPR(gO2,gI2))*CpbarFeFvVWmPL(gO1,gI2)*MFv(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fe_1loop_1_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_1_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHm(gI1)))*
      Conj(CpbarFeFvHmPR(gO2,gI2,gI1))*CpbarFeFvHmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFeFeAhPR(gO2,gI1,gI2))*CpbarFeFeAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFeFehhPR(gO2,gI2,gI1))*CpbarFeFehhPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ)))*Conj(
      CpbarFeFeVZPL(gO2,gI2))*CpbarFeFeVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWm)))*Conj(
      CpbarFeFvVWmPL(gO2,gI2))*CpbarFeFvVWmPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fe_1loop_PR_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_PR_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHm(gI1)))*
      Conj(CpbarFeFvHmPL(gO2,gI2,gI1))*CpbarFeFvHmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFeFeAhPL(gO2,gI1,gI2))*CpbarFeFeAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFeFehhPL(gO2,gI2,gI1))*CpbarFeFehhPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ)))*Conj(
      CpbarFeFeVZPR(gO2,gI2))*CpbarFeFeVZPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWm)))*Conj(
      CpbarFeFvVWmPR(gO2,gI2))*CpbarFeFvVWmPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fe_1loop_PL_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_PL_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_1_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MHm(gI1)))*Conj(
      CpbarFuFdconjHmPL(gO2,gI2,gI1))*CpbarFuFdconjHmPR(gO1,gI2,gI1)*MFd(gI2)));
   result += SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarFuFuAhPL(gO2,gI1,gI2))*CpbarFuFuAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarFuFuhhPL(gO2,gI2,gI1))*CpbarFuFuhhPR(gO1,gI2,gI1)*MFu(gI2)));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm)))*Conj(
      CpbarFuFdconjVWmPR(gO2,gI2))*CpbarFuFdconjVWmPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),0))*Conj(CpbarFuFuVPPR
      (gO2,gI2))*CpbarFuFuVPPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ)))*Conj(
      CpbarFuFuVZPR(gO2,gI2))*CpbarFuFuVZPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_1_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_1_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHm(gI1)))*
      Conj(CpbarFuFdconjHmPR(gO2,gI2,gI1))*CpbarFuFdconjHmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFuFuAhPR(gO2,gI1,gI2))*CpbarFuFuAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFuFuhhPR(gO2,gI2,gI1))*CpbarFuFuhhPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm)))*Conj(
      CpbarFuFdconjVWmPL(gO2,gI2))*CpbarFuFdconjVWmPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),0))*Conj(CpbarFuFuVPPL(
      gO2,gI2))*CpbarFuFuVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ)))*Conj(
      CpbarFuFuVZPL(gO2,gI2))*CpbarFuFuVZPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_PR_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PR_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHm(gI1)))*
      Conj(CpbarFuFdconjHmPL(gO2,gI2,gI1))*CpbarFuFdconjHmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFuFuAhPL(gO2,gI1,gI2))*CpbarFuFuAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFuFuhhPL(gO2,gI2,gI1))*CpbarFuFuhhPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm)))*Conj(
      CpbarFuFdconjVWmPR(gO2,gI2))*CpbarFuFdconjVWmPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),0))*Conj(CpbarFuFuVPPR(
      gO2,gI2))*CpbarFuFuVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ)))*Conj(
      CpbarFuFuVZPR(gO2,gI2))*CpbarFuFuVZPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_PL_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PL_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_1_heavy(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MHm(gI1)))*Conj(
      CpbarUFuFdconjHmPL(gO2,gI2,gI1))*CpbarUFuFdconjHmPR(gO1,gI2,gI1)*MFd(gI2)));
   result += SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarUFuFuAhPL(gO2,gI1,gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarUFuFuhhPL(gO2,gI2,gI1))*CpbarUFuFuhhPR(gO1,gI2,gI1)*MFu(gI2)));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm)))*Conj(
      CpbarUFuFdconjVWmPR(gO2,gI2))*CpbarUFuFdconjVWmPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),0))*Conj(
      CpbarUFuFuVPPR(gO2,gI2))*CpbarUFuFuVPPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFuFuVZPR(gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_1_heavy(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_1_heavy(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_PR_heavy(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHm(gI1)))*
      Conj(CpbarUFuFdconjHmPR(gO2,gI2,gI1))*CpbarUFuFdconjHmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFuFuAhPR(gO2,gI1,gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFuFuhhPR(gO2,gI2,gI1))*CpbarUFuFuhhPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm)))*Conj(
      CpbarUFuFdconjVWmPL(gO2,gI2))*CpbarUFuFdconjVWmPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),0))*Conj(CpbarUFuFuVPPL(
      gO2,gI2))*CpbarUFuFuVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFuFuVZPL(gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_PR_heavy(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PR_heavy(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_PL_heavy(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHm(gI1)))*
      Conj(CpbarUFuFdconjHmPL(gO2,gI2,gI1))*CpbarUFuFdconjHmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFuFuAhPL(gO2,gI1,gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFuFuhhPL(gO2,gI2,gI1))*CpbarUFuFuhhPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm)))*Conj(
      CpbarUFuFdconjVWmPR(gO2,gI2))*CpbarUFuFdconjVWmPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),0))*Conj(CpbarUFuFuVPPR(
      gO2,gI2))*CpbarUFuFuVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFuFuVZPR(gO2,gI2))*CpbarUFuFuVZPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_PL_heavy(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PL_heavy(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::tadpole_hh_1loop(int gO1) const
{
   std::complex<double> result;

   result += A0(Sqr(MVWm))*CpbargWmCgWmCUhh(gO1);
   result += A0(Sqr(MVWm))*CpbargWmgWmUhh(gO1);
   result += A0(Sqr(MVZ))*CpbargZgZUhh(gO1);
   result += 2*CpUhhconjVWmVWm(gO1)*(2*A0(Sqr(MVWm)) - Sqr(MVWm));
   result += CpUhhVZVZ(gO1)*(2*A0(Sqr(MVZ)) - Sqr(MVZ));
   result += -SUM(gI1,0,1,A0(Sqr(MHm(gI1)))*CpUhhHmconjHm(gO1,gI1,gI1));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUhh(gI1,gI1,gO1));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUhh(gI1,gI1,gO1));
   result += 6*SUM(gI1,0,2,A0(Sqr(MFd(gI1)))*(CpbarFdFdUhhPL(gI1,gI1,gO1) +
      CpbarFdFdUhhPR(gI1,gI1,gO1))*MFd(gI1));
   result += 2*SUM(gI1,0,2,A0(Sqr(MFe(gI1)))*(CpbarFeFeUhhPL(gI1,gI1,gO1) +
      CpbarFeFeUhhPR(gI1,gI1,gO1))*MFe(gI1));
   result += 6*SUM(gI1,0,2,A0(Sqr(MFu(gI1)))*(CpbarFuFuUhhPL(gI1,gI1,gO1) +
      CpbarFuFuUhhPR(gI1,gI1,gO1))*MFu(gI1));

   return result * oneOver16PiSqr;
}













void CLASSNAME::calculate_MVG_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MVG) = 0.;
}

void CLASSNAME::calculate_MFv_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MFv).setConstant(0.);
}

void CLASSNAME::calculate_MVP_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MVP) = 0.;
}

void CLASSNAME::calculate_MVZ_pole()
{
   if (!force_output && problems.is_running_tachyon(THDMIISNMSSMBCsimple_info::
      VZ))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVZ));
   const double p = MVZ;
   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(THDMIISNMSSMBCsimple_info::VZ);

   PHYSICAL(MVZ) = AbsSqrt(mass_sqr);
}

void CLASSNAME::calculate_Mhh_pole()
{
   if (!force_output && problems.is_running_tachyon(THDMIISNMSSMBCsimple_info::
      hh))
      return;

   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(Mhh) old_Mhh(Mhh), new_Mhh(Mhh);

   do {
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_hh());

      for (int es = 0; es < 3; ++es) {
         const double p = Abs(old_Mhh(es));
         Eigen::Matrix<double,3,3> self_energy = Re(self_energy_hh_1loop(p));
         const Eigen::Matrix<double,3,3> M_loop(M_tree - self_energy);
         Eigen::Array<double,3,1> eigen_values;
         Eigen::Matrix<double,3,3> mix_ZH;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZH,
               eigenvalue_error);
            problems.flag_bad_mass(THDMIISNMSSMBCsimple_info::hh,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else

            fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZH);
         #endif
            normalize_to_interval(mix_ZH);

         PHYSICAL(Mhh(es)) = SignedAbsSqrt(eigen_values(es));
         if (es == 0)
            PHYSICAL(ZH) = mix_ZH;
      }

      new_Mhh = PHYSICAL(Mhh);
      diff = MaxRelDiff(new_Mhh, old_Mhh);
      old_Mhh = new_Mhh;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(THDMIISNMSSMBCsimple_info::hh);
   else
      problems.unflag_no_pole_mass_convergence(THDMIISNMSSMBCsimple_info::hh);
}

void CLASSNAME::calculate_MAh_pole()
{
   if (!force_output && problems.is_running_tachyon(THDMIISNMSSMBCsimple_info::
      Ah))
      return;

   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(MAh) old_MAh(MAh), new_MAh(MAh);

   do {
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Ah());

      for (int es = 0; es < 3; ++es) {
         const double p = Abs(old_MAh(es));
         Eigen::Matrix<double,3,3> self_energy = Re(self_energy_Ah_1loop(p));
         const Eigen::Matrix<double,3,3> M_loop(M_tree - self_energy);
         Eigen::Array<double,3,1> eigen_values;
         Eigen::Matrix<double,3,3> mix_ZA;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZA,
               eigenvalue_error);
            problems.flag_bad_mass(THDMIISNMSSMBCsimple_info::Ah,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else

            fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZA);
         #endif
            normalize_to_interval(mix_ZA);

         PHYSICAL(MAh(es)) = SignedAbsSqrt(eigen_values(es));
         if (es == 1)
            PHYSICAL(ZA) = mix_ZA;
      }

      new_MAh = PHYSICAL(MAh);
      diff = MaxRelDiff(new_MAh, old_MAh);
      old_MAh = new_MAh;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(THDMIISNMSSMBCsimple_info::Ah);
   else
      problems.unflag_no_pole_mass_convergence(THDMIISNMSSMBCsimple_info::Ah);
}

void CLASSNAME::calculate_MHm_pole()
{
   if (!force_output && problems.is_running_tachyon(THDMIISNMSSMBCsimple_info::
      Hm))
      return;

   // diagonalization with medium precision
   const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Hm());

   for (int es = 0; es < 2; ++es) {
      const double p = Abs(MHm(es));
      Eigen::Matrix<double,2,2> self_energy = Re(self_energy_Hm_1loop(p));
      const Eigen::Matrix<double,2,2> M_loop(M_tree - self_energy);
      Eigen::Array<double,2,1> eigen_values;
      Eigen::Matrix<double,2,2> mix_ZP;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZP,
            eigenvalue_error);
         problems.flag_bad_mass(THDMIISNMSSMBCsimple_info::Hm, eigenvalue_error
             > precision * Abs(eigen_values(0)));
      #else

         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZP);
      #endif
         normalize_to_interval(mix_ZP);

      PHYSICAL(MHm(es)) = SignedAbsSqrt(eigen_values(es));
      if (es == 1)
         PHYSICAL(ZP) = mix_ZP;
   }
}

void CLASSNAME::calculate_MFd_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fd());
   for (int es = 0; es < 3; ++es) {
      const double p = Abs(MFd(es));
      const Eigen::Matrix<double,3,3> self_energy_1  = Re(
         self_energy_Fd_1loop_1(p));
      const Eigen::Matrix<double,3,3> self_energy_PL = Re(
         self_energy_Fd_1loop_PL(p));
      const Eigen::Matrix<double,3,3> self_energy_PR = Re(
         self_energy_Fd_1loop_PR(p));
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Vd) mix_Vd;
      decltype(Ud) mix_Ud;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_Vd, mix_Ud, eigenvalue_error);
      problems.flag_bad_mass(THDMIISNMSSMBCsimple_info::Fd, eigenvalue_error >
         precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_Vd, mix_Ud);
   #endif
      if (es == 0) {
         PHYSICAL(Vd) = mix_Vd;
         PHYSICAL(Ud) = mix_Ud;
      }
      PHYSICAL(MFd(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MFu_pole()
{
   // diagonalization with medium precision
   double qcd_1l = 0.;

   {
      const double currentScale = get_scale();
      qcd_1l = 0.008443431970194815*(-4. + 3.*Log(Sqr(MFu(2))/Sqr(currentScale)
         ))*Sqr(g3);
   }

   double qcd_2l = 0.;

   if (pole_mass_loop_order > 1 && TOP_POLE_QCD_CORRECTION > 0) {
      const double currentScale = get_scale();
      qcd_2l = 2.2278607323533713e-6*Quad(g3)*(-2372.129769909197 + 1452.*Log(
         Sqr(MFu(2))/Sqr(currentScale)) - 396.*Sqr(Log(Sqr(MFu(2))/Sqr(
         currentScale))));
   }

   double qcd_3l = 0.;

   if (pole_mass_loop_order > 2 && TOP_POLE_QCD_CORRECTION > 1) {
      const double currentScale = get_scale();
      qcd_3l = 1.0450439184830051e-10*Power6(g3)*(-1.1451958668162137e7 +
         320760.*Cube(Log(Sqr(MFu(2))/Sqr(currentScale))) + 4.935665892982448e6
         *Log(Sqr(MFu(2))/Sqr(currentScale)) - 1.7253e6*Sqr(Log(Sqr(MFu(2))/Sqr
         (currentScale))));
   }

   double qcd_4l = 0.;

   if (pole_mass_loop_order > 3 && TOP_POLE_QCD_CORRECTION > 2) {
      const double currentScale = get_scale();
      qcd_4l = -1.6081297554549208e-9*Power8(g3)*(211681.74421123447 - 5638.*
         Cube(Log(Sqr(MFu(2))/Sqr(currentScale))) - 104673.38261571848*Log(Sqr(
         MFu(2))/Sqr(currentScale)) + 825.*Quad(Log(Sqr(MFu(2))/Sqr(
         currentScale))) + 22162.91142653778*Sqr(Log(Sqr(MFu(2))/Sqr(
         currentScale))));
   }

   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fu());
   for (int es = 0; es < 3; ++es) {
      const double p = Abs(MFu(es));
      Eigen::Matrix<double,3,3> self_energy_1;
      Eigen::Matrix<double,3,3> self_energy_PL;
      Eigen::Matrix<double,3,3> self_energy_PR;
      for (int i1 = 0; i1 < 3; ++i1) {
         for (int i2 = 0; i2 < 3; ++i2) {
            if (i1 == 2 && i2 == 2) {
               self_energy_1(i1,i2)  = Re(self_energy_Fu_1loop_1_heavy(p,i1,i2)
                  );
               self_energy_PL(i1,i2) = Re(self_energy_Fu_1loop_PL_heavy(p,i1,i2
                  ));
               self_energy_PR(i1,i2) = Re(self_energy_Fu_1loop_PR_heavy(p,i1,i2
                  ));
            } else {
               self_energy_1(i1,i2)  = Re(self_energy_Fu_1loop_1(p,i1,i2));
               self_energy_PL(i1,i2) = Re(self_energy_Fu_1loop_PL(p,i1,i2));
               self_energy_PR(i1,i2) = Re(self_energy_Fu_1loop_PR(p,i1,i2));
            }
         }
      }
      Eigen::Matrix<double,3,3> delta_M(- self_energy_PR * M_tree - M_tree *
         self_energy_PL - self_energy_1);
      delta_M(2,2) -= M_tree(2,2) * (qcd_1l + qcd_2l + qcd_3l + qcd_4l);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Vu) mix_Vu;
      decltype(Uu) mix_Uu;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_Vu, mix_Uu, eigenvalue_error);
      problems.flag_bad_mass(THDMIISNMSSMBCsimple_info::Fu, eigenvalue_error >
         precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_Vu, mix_Uu);
   #endif
      if (es == 0) {
         PHYSICAL(Vu) = mix_Vu;
         PHYSICAL(Uu) = mix_Uu;
      }
      PHYSICAL(MFu(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MFe_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fe());
   for (int es = 0; es < 3; ++es) {
      const double p = Abs(MFe(es));
      const Eigen::Matrix<double,3,3> self_energy_1  = Re(
         self_energy_Fe_1loop_1(p));
      const Eigen::Matrix<double,3,3> self_energy_PL = Re(
         self_energy_Fe_1loop_PL(p));
      const Eigen::Matrix<double,3,3> self_energy_PR = Re(
         self_energy_Fe_1loop_PR(p));
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Ve) mix_Ve;
      decltype(Ue) mix_Ue;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_Ve, mix_Ue, eigenvalue_error);
      problems.flag_bad_mass(THDMIISNMSSMBCsimple_info::Fe, eigenvalue_error >
         precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_Ve, mix_Ue);
   #endif
      if (es == 0) {
         PHYSICAL(Ve) = mix_Ve;
         PHYSICAL(Ue) = mix_Ue;
      }
      PHYSICAL(MFe(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MVWm_pole()
{
   if (!force_output && problems.is_running_tachyon(THDMIISNMSSMBCsimple_info::
      VWm))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVWm));
   const double p = MVWm;
   const double self_energy = Re(self_energy_VWm_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(THDMIISNMSSMBCsimple_info::VWm);

   PHYSICAL(MVWm) = AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWm_pole(double p)
{
   if (!force_output && problems.is_running_tachyon(THDMIISNMSSMBCsimple_info::
      VWm))
      return 0.;

   const double self_energy = Re(self_energy_VWm_1loop(p));
   const double mass_sqr = Sqr(MVWm) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(THDMIISNMSSMBCsimple_info::VWm);

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVZ_pole(double p)
{
   if (!force_output && problems.is_running_tachyon(THDMIISNMSSMBCsimple_info::
      VZ))
      return 0.;

   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = Sqr(MVZ) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(THDMIISNMSSMBCsimple_info::VZ);

   return AbsSqrt(mass_sqr);
}



double CLASSNAME::calculate_MFv_DRbar(double, int) const
{
   return 0.0;
}

double CLASSNAME::calculate_MFe_DRbar(double m_sm_msbar, int idx) const
{
   const double p = m_sm_msbar;
   const double self_energy_1  = Re(self_energy_Fe_1loop_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fe_1loop_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fe_1loop_PR_heavy_rotated(p,
      idx, idx));
   const double drbar_conversion = 1;
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;
   const double delta_mf_1loop = - self_energy_1/m_sm_drbar - self_energy_PL -
      self_energy_PR;

   const double m_susy_drbar = m_sm_drbar / (1.0 + delta_mf_1loop);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MFu_DRbar(double m_pole, int idx) const
{
   const double p = m_pole;
   const double self_energy_1  = Re(self_energy_Fu_1loop_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fu_1loop_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fu_1loop_PR_heavy_rotated(p,
      idx, idx));

   const double currentScale = get_scale();
   double qcd_1l = 0., qcd_2l = 0., qcd_3l = 0., qcd_4l = 0.;
   double atas_S_2l = 0., atas_LR_2l = 0., atat_S_2l = 0., atat_LR_2l = 0.;

   qcd_1l = 0.008443431970194815*(-4. + 3.*Log(Sqr(MFu(idx))/Sqr(currentScale))
      )*Sqr(g3);

   if (get_thresholds() > 1 && threshold_corrections.mt > 1) {
      const double q_2l = 2.2278607323533713e-6*Quad(g3)*(2372.129769909197 -
         1452.*Log(Sqr(MFu(idx))/Sqr(currentScale)) + 396.*Sqr(Log(Sqr(MFu(idx)
         )/Sqr(currentScale))));

      qcd_2l = -q_2l + qcd_1l * qcd_1l;
   }

   if (get_thresholds() > 2 && threshold_corrections.mt > 2) {
      qcd_3l = 1.0450439184830051e-10*Power6(g3)*(-8.404731799492894e6 + 48600.
         *Cube(Log(Sqr(MFu(idx))/Sqr(currentScale))) + 1.112325741480515e6*Log(
         Sqr(MFu(idx))/Sqr(currentScale)) - 208980.*Sqr(Log(Sqr(MFu(idx))/Sqr(
         currentScale))));
   }

   if (get_thresholds() > 3 && threshold_corrections.mt > 3) {
      qcd_4l = 1.6544544809207004e-13*Power8(g3)*(-1.5015630809970136e9 +
         3.1428e6*Cube(Log(Sqr(MFu(idx))/Sqr(currentScale))) +
         4.409910543513119e8*Log(Sqr(MFu(idx))/Sqr(currentScale)) - 826200.*
         Quad(Log(Sqr(MFu(idx))/Sqr(currentScale))) - 1.7811910793512717e7*Sqr(
         Log(Sqr(MFu(idx))/Sqr(currentScale))));
   }

   const double m_susy_drbar = m_pole + self_energy_1 + atas_S_2l + atat_S_2l +
      m_pole * (self_energy_PL + self_energy_PR + qcd_1l + qcd_2l + qcd_3l +
      qcd_4l + atas_LR_2l + atat_LR_2l);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MFd_DRbar(double m_sm_msbar, int idx) const
{
   const double p = m_sm_msbar;
   const double self_energy_1  = Re(self_energy_Fd_1loop_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fd_1loop_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fd_1loop_PR_heavy_rotated(p,
      idx, idx));
   const double m_tree = MFd(idx);
   const double drbar_conversion = 1;
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;
   const double delta_mb_1loop = - self_energy_1/m_tree - self_energy_PL -
      self_energy_PR;
   double qcd_2l = 0.;

   const double m_susy_drbar = m_sm_drbar / (1.0 + delta_mb_1loop + qcd_2l);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MVP_DRbar(double)
{
   return 0.0;
}

double CLASSNAME::calculate_MVZ_DRbar(double m_pole)
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(THDMIISNMSSMBCsimple_info::VZ);
      return m_pole;
   }

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWm_DRbar(double m_pole)
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VWm_1loop(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(THDMIISNMSSMBCsimple_info::VWm);
      return m_pole;
   }

   return AbsSqrt(mass_sqr);
}



double CLASSNAME::v() const
{

   return Sqrt(Sqr(v1) + Sqr(v2));
}

double CLASSNAME::Betax() const
{

   return ArcSin(Abs(ZP(0,1)));
}

double CLASSNAME::ThetaW() const
{

   return ArcCos(Abs(ZZ(0,0)));
}

double CLASSNAME::VEV() const
{

   return Sqrt(Sqr(v1) + Sqr(v2));
}



std::ostream& operator<<(std::ostream& ostr, const CLASSNAME& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
