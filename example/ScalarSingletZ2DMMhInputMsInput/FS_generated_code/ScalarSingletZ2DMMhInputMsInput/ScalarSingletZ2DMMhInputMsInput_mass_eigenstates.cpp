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
 * @file ScalarSingletZ2DMMhInputMsInput_mass_eigenstates.cpp
 * @brief implementation of the ScalarSingletZ2DMMhInputMsInput model class
 *
 * Contains the definition of the ScalarSingletZ2DMMhInputMsInput model class methods
 * which solve EWSB and calculate pole masses and mixings from MSbar
 * parameters.
 *
 * This file was generated with FlexibleSUSY 2.6.1 and SARAH 4.14.3 .
 */

#include "ScalarSingletZ2DMMhInputMsInput_mass_eigenstates.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_ewsb_solver_interface.hpp"
#include "config.h"
#include "eigen_utils.hpp"
#include "ewsb_solver.hpp"
#include "wrappers.hpp"
#include "linalg2.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "loop_libraries/loop_library.hpp"
#include "raii.hpp"

#ifdef ENABLE_THREADS
#include "thread_pool.hpp"
#endif

#ifdef ENABLE_TWO_SCALE_SOLVER
#include "ScalarSingletZ2DMMhInputMsInput_two_scale_ewsb_solver.hpp"
#endif






#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <algorithm>

namespace flexiblesusy {

#define STRINGIFY(s) XSTRINGIFY(s)
#define XSTRINGIFY(s) #s
#define CLASSNAME ScalarSingletZ2DMMhInputMsInput_mass_eigenstates

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

CLASSNAME::CLASSNAME(const ScalarSingletZ2DMMhInputMsInput_input_parameters& input_)
   : ScalarSingletZ2DMMhInputMsInput_soft_parameters(input_)
#if defined(ENABLE_TWO_SCALE_SOLVER)
   , ewsb_solver(new ScalarSingletZ2DMMhInputMsInput_ewsb_solver<Two_scale>())
#endif
{
}

std::unique_ptr<ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_interface> CLASSNAME::clone() const
{
   return std::make_unique<ScalarSingletZ2DMMhInputMsInput_mass_eigenstates>(*this);
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

const ScalarSingletZ2DMMhInputMsInput_physical& CLASSNAME::get_physical() const
{
   return physical;
}

ScalarSingletZ2DMMhInputMsInput_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const ScalarSingletZ2DMMhInputMsInput_physical& physical_)
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

void CLASSNAME::set_ewsb_solver(const std::shared_ptr<ScalarSingletZ2DMMhInputMsInput_ewsb_solver_interface>& solver)
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

   if (ewsb_loop_order > 0) {
      tadpole[0] -= Re(tadpole_hh_1loop());

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

   tadpole[0] /= v;


   return tadpole;
}

int CLASSNAME::solve_ewsb_tree_level_custom()
{
   int error = EWSB_solver::SUCCESS;

   
   const double old_muH2 = muH2;

   muH2 = Re(-0.5*LamH*Sqr(v));

   const bool is_finite = IsFinite(muH2);

   if (!is_finite) {
      muH2 = old_muH2;
      error = EWSB_solver::FAIL;
   }

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

   VERBOSE_MSG("\t\tSolving ScalarSingletZ2DMMhInputMsInput EWSB at " << ewsb_loop_order << "-loop order");

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
           "ScalarSingletZ2DMMhInputMsInput\n"
           "========================================\n";
   ScalarSingletZ2DMMhInputMsInput_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level MSbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MHp = " << MHp << '\n';
   ostr << "Mss = " << Mss << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MAh = " << MAh << '\n';
   ostr << "Mhh = " << Mhh << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MVWp = " << MVWp << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';

   ostr << "----------------------------------------\n"
           "tree-level MSbar mixing matrices:\n"
           "----------------------------------------\n";
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
   return Loop_library::get().A0(m, Sqr(get_scale())).real();
}

double CLASSNAME::B0(double p, double m1, double m2) const noexcept
{
   return Loop_library::get().B0(p, m1, m2, Sqr(get_scale())).real();
}

double CLASSNAME::B1(double p, double m1, double m2) const noexcept
{
   return Loop_library::get().B1(p, m1, m2, Sqr(get_scale())).real();
}

double CLASSNAME::B00(double p, double m1, double m2) const noexcept
{
   return Loop_library::get().B00(p, m1, m2, Sqr(get_scale())).real();
}

double CLASSNAME::B22(double p, double m1, double m2) const noexcept
{
   const double scl2 = Sqr(get_scale());
   return (
      Loop_library::get().B00(p, m1, m2, scl2) - Loop_library::get().A0(m1, scl2)/4.0 -
      Loop_library::get().A0(m2, scl2)/4.0
   ).real();
}

double CLASSNAME::H0(double p, double m1, double m2) const noexcept
{
   const double scl2 = Sqr(get_scale());
   return 4.0*Loop_library::get().B00(p, m1, m2, scl2).real() + CLASSNAME::G0(p, m1, m2);
}

double CLASSNAME::F0(double p, double m1, double m2) const noexcept
{
   const double scl2 = Sqr(get_scale());
   return (
      Loop_library::get().A0(m1, scl2) - 2.0*Loop_library::get().A0(m2, scl2)
      - (2.0*p + 2.0*m1 - m2) * Loop_library::get().B0(p, m1, m2, scl2)
   ).real();
}

double CLASSNAME::G0(double p, double m1, double m2) const noexcept
{
   const double scl2 = Sqr(get_scale());
   return (
      (p - m1 - m2) * Loop_library::get().B0(p, m1, m2, scl2)
      - Loop_library::get().A0(m1, scl2) - Loop_library::get().A0(m2, scl2)
   ).real();
}

/**
 * routine which finds the MSbar mass eigenstates and mixings.
 */
void CLASSNAME::calculate_DRbar_masses()
{
   const auto save_muH2_raii = make_raii_save(muH2);

   solve_ewsb_tree_level_custom();

   calculate_MVPVZ();
   calculate_MVWp();
   calculate_MFe();
   calculate_MFu();
   calculate_MFd();
   calculate_Mhh();
   calculate_MAh();
   calculate_MFv();
   calculate_Mss();
   calculate_MHp();
   calculate_MVG();

}

/**
 * routine which finds the pole mass eigenstates and mixings.
 */
void CLASSNAME::calculate_pole_masses()
{
#ifdef ENABLE_THREADS
   Thread_pool tp(std::min(std::thread::hardware_concurrency(), 10u));

   if (calculate_bsm_pole_masses) {
      tp.run_task([this] () { calculate_Mss_pole(); });
   }

   if (calculate_sm_pole_masses) {
      tp.run_task([this] () { calculate_MVG_pole(); });
      tp.run_task([this] () { calculate_MFv_pole(); });
      tp.run_task([this] () { calculate_Mhh_pole(); });
      tp.run_task([this] () { calculate_MVP_pole(); });
      tp.run_task([this] () { calculate_MVZ_pole(); });
      tp.run_task([this] () { calculate_MFd_pole(); });
      tp.run_task([this] () { calculate_MFu_pole(); });
      tp.run_task([this] () { calculate_MFe_pole(); });
      tp.run_task([this] () { calculate_MVWp_pole(); });
   }

#else
   if (calculate_bsm_pole_masses) {
      calculate_Mss_pole();
   }

   if (calculate_sm_pole_masses) {
      calculate_MVG_pole();
      calculate_MFv_pole();
      calculate_Mhh_pole();
      calculate_MVP_pole();
      calculate_MVZ_pole();
      calculate_MFd_pole();
      calculate_MFu_pole();
      calculate_MFe_pole();
      calculate_MVWp_pole();
   }

#endif
}

void CLASSNAME::copy_DRbar_masses_to_pole_masses()
{
   PHYSICAL(MVG) = MVG;
   PHYSICAL(MHp) = MHp;
   PHYSICAL(Mss) = Mss;
   PHYSICAL(MFv) = MFv;
   PHYSICAL(MAh) = MAh;
   PHYSICAL(Mhh) = Mhh;
   PHYSICAL(MFd) = MFd;
   PHYSICAL(Vd) = Vd;
   PHYSICAL(Ud) = Ud;
   PHYSICAL(MFu) = MFu;
   PHYSICAL(Vu) = Vu;
   PHYSICAL(Uu) = Uu;
   PHYSICAL(MFe) = MFe;
   PHYSICAL(Ve) = Ve;
   PHYSICAL(Ue) = Ue;
   PHYSICAL(MVWp) = MVWp;
   PHYSICAL(MVP) = MVP;
   PHYSICAL(MVZ) = MVZ;
   PHYSICAL(ZZ) = ZZ;

}

/**
 * reorders MSbar masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_DRbar_masses()
{

}

/**
 * reorders pole masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_pole_masses()
{

}

/**
 * Checks the pole masses for tachyons
 */
void CLASSNAME::check_pole_masses_for_tachyons()
{
   if (PHYSICAL(Mss) < 0.) { problems.flag_pole_tachyon(ScalarSingletZ2DMMhInputMsInput_info::ss); }
   if (PHYSICAL(Mhh) < 0.) { problems.flag_pole_tachyon(ScalarSingletZ2DMMhInputMsInput_info::hh); }
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
   MHp = 0.;
   Mss = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MAh = 0.;
   Mhh = 0.;
   MFd = Eigen::Matrix<double,3,1>::Zero();
   Vd = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ud = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   Vu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Uu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   Ve = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ue = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MVWp = 0.;
   MVP = 0.;
   MVZ = 0.;



}

void CLASSNAME::clear_problems()
{
   problems.clear();
}

void CLASSNAME::clear()
{
   ScalarSingletZ2DMMhInputMsInput_soft_parameters::clear();
   clear_DRbar_parameters();
   physical.clear();
   problems.clear();
}

void CLASSNAME::set_DRbar_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MHp = pars(1);
   Mss = pars(2);
   MFv(0) = pars(3);
   MFv(1) = pars(4);
   MFv(2) = pars(5);
   MAh = pars(6);
   Mhh = pars(7);
   MFd(0) = pars(8);
   MFd(1) = pars(9);
   MFd(2) = pars(10);
   MFu(0) = pars(11);
   MFu(1) = pars(12);
   MFu(2) = pars(13);
   MFe(0) = pars(14);
   MFe(1) = pars(15);
   MFe(2) = pars(16);
   MVWp = pars(17);
   MVP = pars(18);
   MVZ = pars(19);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses() const
{
   Eigen::ArrayXd pars(20);

   pars(0) = MVG;
   pars(1) = MHp;
   pars(2) = Mss;
   pars(3) = MFv(0);
   pars(4) = MFv(1);
   pars(5) = MFv(2);
   pars(6) = MAh;
   pars(7) = Mhh;
   pars(8) = MFd(0);
   pars(9) = MFd(1);
   pars(10) = MFd(2);
   pars(11) = MFu(0);
   pars(12) = MFu(1);
   pars(13) = MFu(2);
   pars(14) = MFe(0);
   pars(15) = MFe(1);
   pars(16) = MFe(2);
   pars(17) = MVWp;
   pars(18) = MVP;
   pars(19) = MVZ;

   return pars;
}

void CLASSNAME::set_DRbar_masses_and_mixings(const Eigen::ArrayXd& pars)
{
   set_DRbar_masses(pars);

   Vd(0,0) = std::complex<double>(pars(20), pars(21));
   Vd(0,1) = std::complex<double>(pars(22), pars(23));
   Vd(0,2) = std::complex<double>(pars(24), pars(25));
   Vd(1,0) = std::complex<double>(pars(26), pars(27));
   Vd(1,1) = std::complex<double>(pars(28), pars(29));
   Vd(1,2) = std::complex<double>(pars(30), pars(31));
   Vd(2,0) = std::complex<double>(pars(32), pars(33));
   Vd(2,1) = std::complex<double>(pars(34), pars(35));
   Vd(2,2) = std::complex<double>(pars(36), pars(37));
   Ud(0,0) = std::complex<double>(pars(38), pars(39));
   Ud(0,1) = std::complex<double>(pars(40), pars(41));
   Ud(0,2) = std::complex<double>(pars(42), pars(43));
   Ud(1,0) = std::complex<double>(pars(44), pars(45));
   Ud(1,1) = std::complex<double>(pars(46), pars(47));
   Ud(1,2) = std::complex<double>(pars(48), pars(49));
   Ud(2,0) = std::complex<double>(pars(50), pars(51));
   Ud(2,1) = std::complex<double>(pars(52), pars(53));
   Ud(2,2) = std::complex<double>(pars(54), pars(55));
   Vu(0,0) = std::complex<double>(pars(56), pars(57));
   Vu(0,1) = std::complex<double>(pars(58), pars(59));
   Vu(0,2) = std::complex<double>(pars(60), pars(61));
   Vu(1,0) = std::complex<double>(pars(62), pars(63));
   Vu(1,1) = std::complex<double>(pars(64), pars(65));
   Vu(1,2) = std::complex<double>(pars(66), pars(67));
   Vu(2,0) = std::complex<double>(pars(68), pars(69));
   Vu(2,1) = std::complex<double>(pars(70), pars(71));
   Vu(2,2) = std::complex<double>(pars(72), pars(73));
   Uu(0,0) = std::complex<double>(pars(74), pars(75));
   Uu(0,1) = std::complex<double>(pars(76), pars(77));
   Uu(0,2) = std::complex<double>(pars(78), pars(79));
   Uu(1,0) = std::complex<double>(pars(80), pars(81));
   Uu(1,1) = std::complex<double>(pars(82), pars(83));
   Uu(1,2) = std::complex<double>(pars(84), pars(85));
   Uu(2,0) = std::complex<double>(pars(86), pars(87));
   Uu(2,1) = std::complex<double>(pars(88), pars(89));
   Uu(2,2) = std::complex<double>(pars(90), pars(91));
   Ve(0,0) = std::complex<double>(pars(92), pars(93));
   Ve(0,1) = std::complex<double>(pars(94), pars(95));
   Ve(0,2) = std::complex<double>(pars(96), pars(97));
   Ve(1,0) = std::complex<double>(pars(98), pars(99));
   Ve(1,1) = std::complex<double>(pars(100), pars(101));
   Ve(1,2) = std::complex<double>(pars(102), pars(103));
   Ve(2,0) = std::complex<double>(pars(104), pars(105));
   Ve(2,1) = std::complex<double>(pars(106), pars(107));
   Ve(2,2) = std::complex<double>(pars(108), pars(109));
   Ue(0,0) = std::complex<double>(pars(110), pars(111));
   Ue(0,1) = std::complex<double>(pars(112), pars(113));
   Ue(0,2) = std::complex<double>(pars(114), pars(115));
   Ue(1,0) = std::complex<double>(pars(116), pars(117));
   Ue(1,1) = std::complex<double>(pars(118), pars(119));
   Ue(1,2) = std::complex<double>(pars(120), pars(121));
   Ue(2,0) = std::complex<double>(pars(122), pars(123));
   Ue(2,1) = std::complex<double>(pars(124), pars(125));
   Ue(2,2) = std::complex<double>(pars(126), pars(127));
   ZZ(0,0) = pars(128);
   ZZ(0,1) = pars(129);
   ZZ(1,0) = pars(130);
   ZZ(1,1) = pars(131);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses_and_mixings() const
{
   Eigen::ArrayXd pars(get_DRbar_masses());

   pars.conservativeResize(132);

   pars(20) = Re(Vd(0,0));
   pars(21) = Im(Vd(0,0));
   pars(22) = Re(Vd(0,1));
   pars(23) = Im(Vd(0,1));
   pars(24) = Re(Vd(0,2));
   pars(25) = Im(Vd(0,2));
   pars(26) = Re(Vd(1,0));
   pars(27) = Im(Vd(1,0));
   pars(28) = Re(Vd(1,1));
   pars(29) = Im(Vd(1,1));
   pars(30) = Re(Vd(1,2));
   pars(31) = Im(Vd(1,2));
   pars(32) = Re(Vd(2,0));
   pars(33) = Im(Vd(2,0));
   pars(34) = Re(Vd(2,1));
   pars(35) = Im(Vd(2,1));
   pars(36) = Re(Vd(2,2));
   pars(37) = Im(Vd(2,2));
   pars(38) = Re(Ud(0,0));
   pars(39) = Im(Ud(0,0));
   pars(40) = Re(Ud(0,1));
   pars(41) = Im(Ud(0,1));
   pars(42) = Re(Ud(0,2));
   pars(43) = Im(Ud(0,2));
   pars(44) = Re(Ud(1,0));
   pars(45) = Im(Ud(1,0));
   pars(46) = Re(Ud(1,1));
   pars(47) = Im(Ud(1,1));
   pars(48) = Re(Ud(1,2));
   pars(49) = Im(Ud(1,2));
   pars(50) = Re(Ud(2,0));
   pars(51) = Im(Ud(2,0));
   pars(52) = Re(Ud(2,1));
   pars(53) = Im(Ud(2,1));
   pars(54) = Re(Ud(2,2));
   pars(55) = Im(Ud(2,2));
   pars(56) = Re(Vu(0,0));
   pars(57) = Im(Vu(0,0));
   pars(58) = Re(Vu(0,1));
   pars(59) = Im(Vu(0,1));
   pars(60) = Re(Vu(0,2));
   pars(61) = Im(Vu(0,2));
   pars(62) = Re(Vu(1,0));
   pars(63) = Im(Vu(1,0));
   pars(64) = Re(Vu(1,1));
   pars(65) = Im(Vu(1,1));
   pars(66) = Re(Vu(1,2));
   pars(67) = Im(Vu(1,2));
   pars(68) = Re(Vu(2,0));
   pars(69) = Im(Vu(2,0));
   pars(70) = Re(Vu(2,1));
   pars(71) = Im(Vu(2,1));
   pars(72) = Re(Vu(2,2));
   pars(73) = Im(Vu(2,2));
   pars(74) = Re(Uu(0,0));
   pars(75) = Im(Uu(0,0));
   pars(76) = Re(Uu(0,1));
   pars(77) = Im(Uu(0,1));
   pars(78) = Re(Uu(0,2));
   pars(79) = Im(Uu(0,2));
   pars(80) = Re(Uu(1,0));
   pars(81) = Im(Uu(1,0));
   pars(82) = Re(Uu(1,1));
   pars(83) = Im(Uu(1,1));
   pars(84) = Re(Uu(1,2));
   pars(85) = Im(Uu(1,2));
   pars(86) = Re(Uu(2,0));
   pars(87) = Im(Uu(2,0));
   pars(88) = Re(Uu(2,1));
   pars(89) = Im(Uu(2,1));
   pars(90) = Re(Uu(2,2));
   pars(91) = Im(Uu(2,2));
   pars(92) = Re(Ve(0,0));
   pars(93) = Im(Ve(0,0));
   pars(94) = Re(Ve(0,1));
   pars(95) = Im(Ve(0,1));
   pars(96) = Re(Ve(0,2));
   pars(97) = Im(Ve(0,2));
   pars(98) = Re(Ve(1,0));
   pars(99) = Im(Ve(1,0));
   pars(100) = Re(Ve(1,1));
   pars(101) = Im(Ve(1,1));
   pars(102) = Re(Ve(1,2));
   pars(103) = Im(Ve(1,2));
   pars(104) = Re(Ve(2,0));
   pars(105) = Im(Ve(2,0));
   pars(106) = Re(Ve(2,1));
   pars(107) = Im(Ve(2,1));
   pars(108) = Re(Ve(2,2));
   pars(109) = Im(Ve(2,2));
   pars(110) = Re(Ue(0,0));
   pars(111) = Im(Ue(0,0));
   pars(112) = Re(Ue(0,1));
   pars(113) = Im(Ue(0,1));
   pars(114) = Re(Ue(0,2));
   pars(115) = Im(Ue(0,2));
   pars(116) = Re(Ue(1,0));
   pars(117) = Im(Ue(1,0));
   pars(118) = Re(Ue(1,1));
   pars(119) = Im(Ue(1,1));
   pars(120) = Re(Ue(1,2));
   pars(121) = Im(Ue(1,2));
   pars(122) = Re(Ue(2,0));
   pars(123) = Im(Ue(2,0));
   pars(124) = Re(Ue(2,1));
   pars(125) = Im(Ue(2,1));
   pars(126) = Re(Ue(2,2));
   pars(127) = Im(Ue(2,2));
   pars(128) = ZZ(0,0);
   pars(129) = ZZ(0,1);
   pars(130) = ZZ(1,0);
   pars(131) = ZZ(1,1);


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
   return "ScalarSingletZ2DMMhInputMsInput";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   ScalarSingletZ2DMMhInputMsInput_soft_parameters::run_to(scale, eps);
}

void CLASSNAME::calculate_tree_level_mass_spectrum()
{
   calculate_DRbar_masses();
}

void CLASSNAME::calculate_pole_mass_spectrum()
{
   calculate_pole_masses();
}

void CLASSNAME::calculate_mass_spectrum()
{
   calculate_spectrum();
}

int CLASSNAME::solve_ewsb_equations_tree_level()
{
   return solve_ewsb_tree_level();
}

int CLASSNAME::solve_ewsb_equations()
{
   return solve_ewsb();
}

const ScalarSingletZ2DMMhInputMsInput_input_parameters& CLASSNAME::get_input_parameters() const
{
   return get_input();
}

ScalarSingletZ2DMMhInputMsInput_input_parameters& CLASSNAME::get_input_parameters()
{
   return get_input();
}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses() const
{
   return get_DRbar_masses();
}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses_and_mixings() const
{
   return get_DRbar_masses_and_mixings();
}

void CLASSNAME::set_tree_level_masses(const Eigen::ArrayXd& m)
{
   set_DRbar_masses(m);
}

void CLASSNAME::set_tree_level_masses_and_mixings(const Eigen::ArrayXd& m)
{
   set_DRbar_masses_and_mixings(m);
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

double CLASSNAME::get_mass_matrix_Hp() const
{

   const double mass_matrix_Hp = Re(muH2 + 0.25*(2*LamH + Sqr(g2))*Sqr(v));

   return mass_matrix_Hp;
}

void CLASSNAME::calculate_MHp()
{

   const auto mass_matrix_Hp = get_mass_matrix_Hp();
   MHp = mass_matrix_Hp;

   if (MHp < 0.) {
      problems.flag_running_tachyon(ScalarSingletZ2DMMhInputMsInput_info::Hp);
   }

   MHp = AbsSqrt(MHp);
}

double CLASSNAME::get_mass_matrix_ss() const
{

   const double mass_matrix_ss = Re(muS2 + 0.5*LamSH*Sqr(v));

   return mass_matrix_ss;
}

void CLASSNAME::calculate_Mss()
{

   const auto mass_matrix_ss = get_mass_matrix_ss();
   Mss = mass_matrix_ss;

   if (Mss < 0.) {
      problems.flag_running_tachyon(ScalarSingletZ2DMMhInputMsInput_info::ss);
   }

   Mss = AbsSqrt(Mss);
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

double CLASSNAME::get_mass_matrix_Ah() const
{

   const double mass_matrix_Ah = Re(0.25*(4*muH2 + 2*LamH*Sqr(v) +
      1.5491933384829668*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(v) + Sqr(g2)*Sqr
      (v)*Sqr(Cos(ThetaW())) + 0.6*Sqr(g1)*Sqr(v)*Sqr(Sin(ThetaW()))));

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{

   const auto mass_matrix_Ah = get_mass_matrix_Ah();
   MAh = mass_matrix_Ah;

   if (MAh < 0.) {
      problems.flag_running_tachyon(ScalarSingletZ2DMMhInputMsInput_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

double CLASSNAME::get_mass_matrix_hh() const
{

   const double mass_matrix_hh = Re(muH2 + 1.5*LamH*Sqr(v));

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{

   const auto mass_matrix_hh = get_mass_matrix_hh();
   Mhh = mass_matrix_hh;

   if (Mhh < 0.) {
      problems.flag_running_tachyon(ScalarSingletZ2DMMhInputMsInput_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fd() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fd;

   mass_matrix_Fd(0,0) = 0.7071067811865475*v*Yd(0,0);
   mass_matrix_Fd(0,1) = 0.7071067811865475*v*Yd(1,0);
   mass_matrix_Fd(0,2) = 0.7071067811865475*v*Yd(2,0);
   mass_matrix_Fd(1,0) = 0.7071067811865475*v*Yd(0,1);
   mass_matrix_Fd(1,1) = 0.7071067811865475*v*Yd(1,1);
   mass_matrix_Fd(1,2) = 0.7071067811865475*v*Yd(2,1);
   mass_matrix_Fd(2,0) = 0.7071067811865475*v*Yd(0,2);
   mass_matrix_Fd(2,1) = 0.7071067811865475*v*Yd(1,2);
   mass_matrix_Fd(2,2) = 0.7071067811865475*v*Yd(2,2);

   return mass_matrix_Fd;
}

void CLASSNAME::calculate_MFd()
{
   const auto mass_matrix_Fd(get_mass_matrix_Fd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud, eigenvalue_error);
   problems.flag_bad_mass(ScalarSingletZ2DMMhInputMsInput_info::Fd,
      eigenvalue_error > precision * Abs(MFd(0)));
#else
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fu() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fu;

   mass_matrix_Fu(0,0) = -0.7071067811865475*v*Yu(0,0);
   mass_matrix_Fu(0,1) = -0.7071067811865475*v*Yu(1,0);
   mass_matrix_Fu(0,2) = -0.7071067811865475*v*Yu(2,0);
   mass_matrix_Fu(1,0) = -0.7071067811865475*v*Yu(0,1);
   mass_matrix_Fu(1,1) = -0.7071067811865475*v*Yu(1,1);
   mass_matrix_Fu(1,2) = -0.7071067811865475*v*Yu(2,1);
   mass_matrix_Fu(2,0) = -0.7071067811865475*v*Yu(0,2);
   mass_matrix_Fu(2,1) = -0.7071067811865475*v*Yu(1,2);
   mass_matrix_Fu(2,2) = -0.7071067811865475*v*Yu(2,2);

   return mass_matrix_Fu;
}

void CLASSNAME::calculate_MFu()
{
   const auto mass_matrix_Fu(get_mass_matrix_Fu());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu, eigenvalue_error);
   problems.flag_bad_mass(ScalarSingletZ2DMMhInputMsInput_info::Fu,
      eigenvalue_error > precision * Abs(MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fe() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fe;

   mass_matrix_Fe(0,0) = 0.7071067811865475*v*Ye(0,0);
   mass_matrix_Fe(0,1) = 0.7071067811865475*v*Ye(1,0);
   mass_matrix_Fe(0,2) = 0.7071067811865475*v*Ye(2,0);
   mass_matrix_Fe(1,0) = 0.7071067811865475*v*Ye(0,1);
   mass_matrix_Fe(1,1) = 0.7071067811865475*v*Ye(1,1);
   mass_matrix_Fe(1,2) = 0.7071067811865475*v*Ye(2,1);
   mass_matrix_Fe(2,0) = 0.7071067811865475*v*Ye(0,2);
   mass_matrix_Fe(2,1) = 0.7071067811865475*v*Ye(1,2);
   mass_matrix_Fe(2,2) = 0.7071067811865475*v*Ye(2,2);

   return mass_matrix_Fe;
}

void CLASSNAME::calculate_MFe()
{
   const auto mass_matrix_Fe(get_mass_matrix_Fe());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue, eigenvalue_error);
   problems.flag_bad_mass(ScalarSingletZ2DMMhInputMsInput_info::Fe,
      eigenvalue_error > precision * Abs(MFe(0)));
#else
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue);
#endif

}

double CLASSNAME::get_mass_matrix_VWp() const
{

   const double mass_matrix_VWp = Re(0.25*Sqr(g2)*Sqr(v));

   return mass_matrix_VWp;
}

void CLASSNAME::calculate_MVWp()
{

   const auto mass_matrix_VWp = get_mass_matrix_VWp();
   MVWp = mass_matrix_VWp;

   if (MVWp < 0.) {
      problems.flag_running_tachyon(ScalarSingletZ2DMMhInputMsInput_info::VWp);
   }

   MVWp = AbsSqrt(MVWp);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_VPVZ() const
{

   Eigen::Matrix<double,2,2> mass_matrix_VPVZ;

   mass_matrix_VPVZ(0,0) = 0.15*Sqr(g1)*Sqr(v);
   mass_matrix_VPVZ(0,1) = -0.19364916731037085*g1*g2*Sqr(v);
   mass_matrix_VPVZ(1,1) = 0.25*Sqr(g2)*Sqr(v);

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
   
   double result = Re(muH2*v + 0.5*LamH*Cube(v));

   return result;
}



double CLASSNAME::CphhHpconjHp() const
{
   
   const double result = -(LamH*v);

   return result;
}

double CLASSNAME::CpbargWpgZHp() const
{
   
   const double result = 0.05*g2*v*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin
      (ThetaW()));

   return result;
}

double CLASSNAME::CpbargZgWpconjHp() const
{
   
   const double result = 0.25*g2*v*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(
      ThetaW()));

   return result;
}

double CLASSNAME::CpbargWpCgZconjHp() const
{
   
   const double result = 0.05*g2*v*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin
      (ThetaW()));

   return result;
}

double CLASSNAME::CpbargZgWpCHp() const
{
   
   const double result = 0.25*g2*v*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(
      ThetaW()));

   return result;
}

double CLASSNAME::CpconjHpVPVWp() const
{
   
   const double result = 0.3872983346207417*g1*g2*v*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpconjHpVWpVZ() const
{
   
   const double result = -0.3872983346207417*g1*g2*v*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpAhAhHpconjHp() const
{
   
   const double result = -LamH;

   return result;
}

double CLASSNAME::CphhhhHpconjHp() const
{
   
   const double result = -LamH;

   return result;
}

double CLASSNAME::CpHpHpconjHpconjHp() const
{
   
   const double result = -2*LamH;

   return result;
}

double CLASSNAME::CpHpssssconjHp() const
{
   
   const double result = -LamSH;

   return result;
}

std::complex<double> CLASSNAME::CpAhconjHpVWp() const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*g2;

   return result;
}

double CLASSNAME::CphhconjHpVWp() const
{
   
   const double result = 0.5*g2;

   return result;
}

double CLASSNAME::CpHpconjHpVP() const
{
   
   const double result = 0.5*(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW(
      )));

   return result;
}

double CLASSNAME::CpHpconjHpVZ() const
{
   
   const double result = 0.5*(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW(
      )));

   return result;
}

double CLASSNAME::CpHpconjHpconjVWpVWp() const
{
   
   const double result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpHpconjHpVZVZ() const
{
   
   const std::complex<double> result = 0.1*(-7.745966692414834*g1*g2*Cos(ThetaW())
      *Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW()))
      );

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFuconjHpPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = -SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(
      gI2,j1))*Vd(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFuconjHpPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -SUM(j2,0,2,Conj(Vu(gI2,j2))*SUM(j1,0,2,
      Conj(Ud(gI1,j1))*Yd(j1,j2)));

   return result;
}

double CLASSNAME::CpbarFeFvconjHpPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFvconjHpPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -SUM(j1,0,2,Conj(Ue(gI1,j1))*Ye(j1,gI2));

   return result;
}

double CLASSNAME::Cphhssss() const
{
   
   const double result = -(LamSH*v);

   return result;
}

double CLASSNAME::CpAhAhssss() const
{
   
   const double result = -LamSH;

   return result;
}

double CLASSNAME::Cphhhhssss() const
{
   
   const double result = -LamSH;

   return result;
}

double CLASSNAME::Cpssssssss() const
{
   
   const double result = -6*LamS;

   return result;
}

double CLASSNAME::CpAhAhhh() const
{
   
   const double result = -(LamH*v);

   return result;
}

std::complex<double> CLASSNAME::CpbargWpgWpAh() const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.25)*v*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbargWpCgWpCAh() const
{
   
   const std::complex<double> result = std::complex<double>(0,0.25)*v*Sqr(g2);

   return result;
}

double CLASSNAME::CpAhAhAhAh() const
{
   
   const double result = -3*LamH;

   return result;
}

double CLASSNAME::CpAhAhhhhh() const
{
   
   const double result = -LamH;

   return result;
}

std::complex<double> CLASSNAME::CpAhhhVZ() const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.1)*(5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpAhHpconjVWp() const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*g2;

   return result;
}

double CLASSNAME::CpAhAhconjVWpVWp() const
{
   
   const double result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhVZVZ() const
{
   
   const std::complex<double> result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834*g2
      *Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gI2,j1))*Vd(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j2,0,2,Conj(Vd(gI2,j2))*SUM(j1,0,2,Conj(Ud(gI1,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gI2,j1))*Ve(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j2,0,2,Conj(Ve(gI2,j2))*SUM(j1,0,2,Conj(Ue(gI1,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gI2,j1))*Vu(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j2,0,2,Conj(Vu(gI2,j2))*SUM(j1,0,2,Conj(Uu(gI1,j1))*Yu(j1,j2)));

   return result;
}

double CLASSNAME::Cphhhhhh() const
{
   
   const double result = -3*LamH*v;

   return result;
}

double CLASSNAME::CphhVZVZ() const
{
   
   const double result = 0.5*v*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(
      ThetaW()));

   return result;
}

double CLASSNAME::CphhconjVWpVWp() const
{
   
   const double result = 0.5*v*Sqr(g2);

   return result;
}

double CLASSNAME::CpbargWpgWphh() const
{
   
   const double result = -0.25*v*Sqr(g2);

   return result;
}

double CLASSNAME::CpbargWpCgWpChh() const
{
   
   const double result = -0.25*v*Sqr(g2);

   return result;
}

double CLASSNAME::CpbargZgZhh() const
{
   
   const double result = -0.25*v*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(
      ThetaW()));

   return result;
}

double CLASSNAME::Cphhhhhhhh() const
{
   
   const double result = -3*LamH;

   return result;
}

double CLASSNAME::CphhHpconjVWp() const
{
   
   const double result = -0.5*g2;

   return result;
}

double CLASSNAME::CphhhhconjVWpVWp() const
{
   
   const double result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CphhhhVZVZ() const
{
   
   const std::complex<double> result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834*g2
      *Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdhhPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Yd(j1,j2))*Ud(gI2,j1))*Vd(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdhhPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(Vd(gI2,
      j2))*SUM(j1,0,2,Conj(Ud(gI1,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFehhPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Ye(j1,j2))*Ue(gI2,j1))*Ve(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFehhPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(Ve(gI2,
      j2))*SUM(j1,0,2,Conj(Ue(gI1,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuhhPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Yu(j1,j2))*Uu(gI2,j1))*Vu(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuhhPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.7071067811865475*SUM(j2,0,2,Conj(Vu(gI2,
      j2))*SUM(j1,0,2,Conj(Uu(gI1,j1))*Yu(j1,j2)));

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

double CLASSNAME::CpHpconjVWpVP() const
{
   
   const double result = 0.3872983346207417*g1*g2*v*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpbargWpgWpVP() const
{
   
   const double result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbargWpCgWpCVP() const
{
   
   const double result = -(g2*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpHpconjHpVPVP() const
{
   
   const std::complex<double> result = 0.1*(g2*Sin(ThetaW())*(7.745966692414834*g1
      *Cos(ThetaW()) + 5*g2*Sin(ThetaW())) + 3*Sqr(g1)*Sqr(Cos(ThetaW())));

   return result;
}

double CLASSNAME::CpconjVWpVPVWp() const
{
   
   const double result = -(g2*Sin(ThetaW()));

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

double CLASSNAME::CpconjVWpVPVPVWp3() const
{
   
   const double result = -2*Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWpVPVPVWp1() const
{
   
   const double result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWpVPVPVWp2() const
{
   
   const double result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpHpconjVWpVZ() const
{
   
   const double result = -0.3872983346207417*g1*g2*v*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbargWpgWpVZ() const
{
   
   const double result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpbargWpCgWpCVZ() const
{
   
   const double result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWpVWpVZ() const
{
   
   const double result = g2*Cos(ThetaW());

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

double CLASSNAME::CpconjVWpVWpVZVZ1() const
{
   
   const double result = -2*Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWpVWpVZVZ2() const
{
   
   const double result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWpVWpVZVZ3() const
{
   
   const double result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpbargPgWpconjVWp() const
{
   
   const double result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbargWpCgPconjVWp() const
{
   
   const double result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbargWpCgZconjVWp() const
{
   
   const double result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpbargZgWpconjVWp() const
{
   
   const double result = -(g2*Cos(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFuconjVWpPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(Vu(
      gI2,j1))*Vd(gI1,j1));

   return result;
}

double CLASSNAME::CpbarFdFuconjVWpPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFvconjVWpPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.7071067811865475*g2*Ve(gI1,
      gI2),0);

   return result;
}

double CLASSNAME::CpbarFeFvconjVWpPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

double CLASSNAME::CpconjVWpconjVWpVWpVWp2() const
{
   
   const double result = -Sqr(g2);

   return result;
}

double CLASSNAME::CpconjVWpconjVWpVWpVWp1() const
{
   
   const double result = 2*Sqr(g2);

   return result;
}

double CLASSNAME::CpconjVWpconjVWpVWpVWp3() const
{
   
   const double result = -Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPL(int gO2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Vd(gI1,j2))*Yd(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPR(int gO1, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,-
      0.7071067811865475)*SUM(j1,0,2,Conj(Yd(j1,gO1))*Ud(gI1,j1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdhhPL(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*SUM(j2,0,2,
      Conj(Vd(gI2,j2))*Yd(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdhhPR(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*SUM(j1,0,2,
      Conj(Yd(j1,gO1))*Ud(gI2,j1)),0);

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

std::complex<double> CLASSNAME::CpbarUFdFuconjHpPL(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-SUM(j2,0,2,Conj(Vu(gI2,j2))*Yd(
      gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFuconjHpPR(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-SUM(j1,0,2,Conj(Yu(j1,gO1))*Uu(
      gI2,j1)),0);

   return result;
}

double CLASSNAME::CpbarUFdFuconjVWpPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFuconjVWpPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(Vu(
      gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPL(int gO2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Vu(gI1,j2))*Yu(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPR(int gO1, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,-
      0.7071067811865475)*SUM(j1,0,2,Conj(Yu(j1,gO1))*Uu(gI1,j1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFdHpPL(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-SUM(j2,0,2,Conj(Vd(gI2,j2))*Yu(
      gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFdHpPR(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-SUM(j1,0,2,Conj(Yd(j1,gO1))*Ud(
      gI2,j1)),0);

   return result;
}

double CLASSNAME::CpbarUFuFdVWpPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFdVWpPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(Vd(
      gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuhhPL(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*SUM(j2,0,2,
      Conj(Vu(gI2,j2))*Yu(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuhhPR(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.7071067811865475*SUM(j1,0,2,
      Conj(Yu(j1,gO1))*Uu(gI2,j1)),0);

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

std::complex<double> CLASSNAME::CpbarUFeFeAhPL(int gO2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Ve(gI1,j2))*Ye(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPR(int gO1, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,-
      0.7071067811865475)*SUM(j1,0,2,Conj(Ye(j1,gO1))*Ue(gI1,j1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFehhPL(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*SUM(j2,0,2,
      Conj(Ve(gI2,j2))*Ye(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFehhPR(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*SUM(j1,0,2,
      Conj(Ye(j1,gO1))*Ue(gI2,j1)),0);

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

std::complex<double> CLASSNAME::CpbarUFeFvconjHpPL(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-Ye(gO2,gI2),0);

   return result;
}

double CLASSNAME::CpbarUFeFvconjHpPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

double CLASSNAME::CpbarUFeFvconjVWpPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

double CLASSNAME::CpbarUFeFvconjVWpPL(int gO1, int gI2) const
{
   
   const double result = IF(gI2 < 3,-0.7071067811865475*g2*KroneckerDelta(gI2,gO1)
      ,0);

   return result;
}

double CLASSNAME::CpbarFvFeHpPL(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvFeHpPR(int gO1, int gI2) const
{
   
   const std::complex<double> result = -SUM(j1,0,2,Conj(Ye(j1,gO1))*Ue(gI2,j1));

   return result;
}

double CLASSNAME::CpbarFvFeVWpPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvFeVWpPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(Ve(
      gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdHpPL(int gO2, int gI2) const
{
   
   const std::complex<double> result = -SUM(j2,0,2,Conj(Vd(gI2,j2))*SUM(j1,0,2,
      Conj(Uu(gO2,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdHpPR(int gO1, int gI2) const
{
   
   const std::complex<double> result = -SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(
      gI2,j1))*Vu(gO1,j2));

   return result;
}

double CLASSNAME::CpbarFuFdVWpPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdVWpPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(Vd(
      gI2,j1))*Vu(gO1,j1));

   return result;
}


std::complex<double> CLASSNAME::self_energy_Hp_1loop(double p ) const
{
   std::complex<double> result;

   result += -0.5*A0(Sqr(MAh))*CpAhAhHpconjHp();
   result += -(B0(Sqr(p),Sqr(MVZ),Sqr(MVWp))*CpbargWpCgZconjHp()*CpbargZgWpCHp());
   result += -(B0(Sqr(p),Sqr(MVWp),Sqr(MVZ))*CpbargWpgZHp()*CpbargZgWpconjHp());
   result += 2*AbsSqr(CpconjHpVPVWp())*(-1 + 2*B0(Sqr(p),0,Sqr(MVWp)));
   result += 2*AbsSqr(CpconjHpVWpVZ())*(-1 + 2*B0(Sqr(p),Sqr(MVWp),Sqr(MVZ)));
   result += -0.5*A0(Sqr(Mhh))*CphhhhHpconjHp();
   result += AbsSqr(CphhHpconjHp())*B0(Sqr(p),Sqr(MHp),Sqr(Mhh));
   result += -(A0(Sqr(MHp))*CpHpHpconjHpconjHp());
   result += -0.5*A0(Sqr(Mss))*CpHpssssconjHp();
   result += AbsSqr(CpAhconjHpVWp())*F0(Sqr(p),Sqr(MAh),Sqr(MVWp));
   result += AbsSqr(CphhconjHpVWp())*F0(Sqr(p),Sqr(Mhh),Sqr(MVWp));
   result += AbsSqr(CpHpconjHpVP())*F0(Sqr(p),Sqr(MHp),0);
   result += AbsSqr(CpHpconjHpVZ())*F0(Sqr(p),Sqr(MHp),Sqr(MVZ));
   result += 4*A0(Sqr(MVWp))*CpHpconjHpconjVWpVWp() - 2*CpHpconjHpconjVWpVWp()*Sqr
      (MVWp);
   result += CpHpconjHpVZVZ()*(2*A0(Sqr(MVZ)) - Sqr(MVZ));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFuconjHpPL(gI1,gI2)) +
      AbsSqr(CpbarFdFuconjHpPR(gI1,gI2)))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFu(gI2)))))
      ;
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFvconjHpPL(gI1,gI2)) + AbsSqr(
      CpbarFeFvconjHpPR(gI1,gI2)))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFv(gI2)))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFu(
      gI2)))*(Conj(CpbarFdFuconjHpPR(gI1,gI2))*CpbarFdFuconjHpPL(gI1,gI2) + Conj(
      CpbarFdFuconjHpPL(gI1,gI2))*CpbarFdFuconjHpPR(gI1,gI2))*MFu(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFv(
      gI2)))*(Conj(CpbarFeFvconjHpPR(gI1,gI2))*CpbarFeFvconjHpPL(gI1,gI2) + Conj(
      CpbarFeFvconjHpPL(gI1,gI2))*CpbarFeFvconjHpPR(gI1,gI2))*MFv(gI2)));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_ss_1loop(double p ) const
{
   std::complex<double> result;

   result += -0.5*A0(Sqr(MAh))*CpAhAhssss();
   result += -0.5*A0(Sqr(Mhh))*Cphhhhssss();
   result += AbsSqr(Cphhssss())*B0(Sqr(p),Sqr(Mss),Sqr(Mhh));
   result += -(A0(Sqr(MHp))*CpHpssssconjHp());
   result += -0.5*A0(Sqr(Mss))*Cpssssssss();

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_Ah_1loop(double p ) const
{
   std::complex<double> result;

   result += -0.5*A0(Sqr(MAh))*CpAhAhAhAh();
   result += AbsSqr(CpAhAhhh())*B0(Sqr(p),Sqr(Mhh),Sqr(MAh));
   result += -0.5*A0(Sqr(Mhh))*CpAhAhhhhh();
   result += -(A0(Sqr(MHp))*CpAhAhHpconjHp());
   result += -0.5*A0(Sqr(Mss))*CpAhAhssss();
   result += -(B0(Sqr(p),Sqr(MVWp),Sqr(MVWp))*Sqr(CpbargWpCgWpCAh()));
   result += -(B0(Sqr(p),Sqr(MVWp),Sqr(MVWp))*Sqr(CpbargWpgWpAh()));
   result += AbsSqr(CpAhhhVZ())*F0(Sqr(p),Sqr(Mhh),Sqr(MVZ));
   result += 2*AbsSqr(CpAhHpconjVWp())*F0(Sqr(p),Sqr(MHp),Sqr(MVWp));
   result += 4*A0(Sqr(MVWp))*CpAhAhconjVWpVWp() - 2*CpAhAhconjVWpVWp()*Sqr(MVWp);
   result += CpAhAhVZVZ()*(2*A0(Sqr(MVZ)) - Sqr(MVZ));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdAhPL(gI1,gI2)) + AbsSqr(
      CpbarFdFdAhPR(gI1,gI2)))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFeAhPL(gI1,gI2)) + AbsSqr(
      CpbarFeFeAhPR(gI1,gI2)))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuAhPL(gI1,gI2)) + AbsSqr(
      CpbarFuFuAhPR(gI1,gI2)))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2)))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(
      gI2)))*(Conj(CpbarFdFdAhPR(gI1,gI2))*CpbarFdFdAhPL(gI1,gI2) + Conj(
      CpbarFdFdAhPL(gI1,gI2))*CpbarFdFdAhPR(gI1,gI2))*MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(
      gI2)))*(Conj(CpbarFeFeAhPR(gI1,gI2))*CpbarFeFeAhPL(gI1,gI2) + Conj(
      CpbarFeFeAhPL(gI1,gI2))*CpbarFeFeAhPR(gI1,gI2))*MFe(gI2)));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(
      gI2)))*(Conj(CpbarFuFuAhPR(gI1,gI2))*CpbarFuFuAhPL(gI1,gI2) + Conj(
      CpbarFuFuAhPL(gI1,gI2))*CpbarFuFuAhPR(gI1,gI2))*MFu(gI2)));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_hh_1loop(double p ) const
{
   std::complex<double> result;

   result += 0.5*AbsSqr(CpAhAhhh())*B0(Sqr(p),Sqr(MAh),Sqr(MAh));
   result += -0.5*A0(Sqr(MAh))*CpAhAhhhhh();
   result += -(B0(Sqr(p),Sqr(MVWp),Sqr(MVWp))*Sqr(CpbargWpCgWpChh()));
   result += -(B0(Sqr(p),Sqr(MVWp),Sqr(MVWp))*Sqr(CpbargWpgWphh()));
   result += -(B0(Sqr(p),Sqr(MVZ),Sqr(MVZ))*Sqr(CpbargZgZhh()));
   result += 2*AbsSqr(CphhconjVWpVWp())*(-1 + 2*B0(Sqr(p),Sqr(MVWp),Sqr(MVWp)));
   result += 0.5*AbsSqr(Cphhhhhh())*B0(Sqr(p),Sqr(Mhh),Sqr(Mhh));
   result += -0.5*A0(Sqr(Mhh))*Cphhhhhhhh();
   result += -(A0(Sqr(MHp))*CphhhhHpconjHp());
   result += -0.5*A0(Sqr(Mss))*Cphhhhssss();
   result += AbsSqr(CphhHpconjHp())*B0(Sqr(p),Sqr(MHp),Sqr(MHp));
   result += 0.5*AbsSqr(Cphhssss())*B0(Sqr(p),Sqr(Mss),Sqr(Mss));
   result += AbsSqr(CphhVZVZ())*(-1 + 2*B0(Sqr(p),Sqr(MVZ),Sqr(MVZ)));
   result += AbsSqr(CpAhhhVZ())*F0(Sqr(p),Sqr(MAh),Sqr(MVZ));
   result += 2*AbsSqr(CphhHpconjVWp())*F0(Sqr(p),Sqr(MHp),Sqr(MVWp));
   result += 4*A0(Sqr(MVWp))*CphhhhconjVWpVWp() - 2*CphhhhconjVWpVWp()*Sqr(MVWp);
   result += CphhhhVZVZ()*(2*A0(Sqr(MVZ)) - Sqr(MVZ));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdhhPL(gI1,gI2)) + AbsSqr(
      CpbarFdFdhhPR(gI1,gI2)))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFehhPL(gI1,gI2)) + AbsSqr(
      CpbarFeFehhPR(gI1,gI2)))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuhhPL(gI1,gI2)) + AbsSqr(
      CpbarFuFuhhPR(gI1,gI2)))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2)))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(
      gI2)))*(Conj(CpbarFdFdhhPR(gI1,gI2))*CpbarFdFdhhPL(gI1,gI2) + Conj(
      CpbarFdFdhhPL(gI1,gI2))*CpbarFdFdhhPR(gI1,gI2))*MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(
      gI2)))*(Conj(CpbarFeFehhPR(gI1,gI2))*CpbarFeFehhPL(gI1,gI2) + Conj(
      CpbarFeFehhPL(gI1,gI2))*CpbarFeFehhPR(gI1,gI2))*MFe(gI2)));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(
      gI2)))*(Conj(CpbarFuFuhhPR(gI1,gI2))*CpbarFuFuhhPL(gI1,gI2) + Conj(
      CpbarFuFuhhPL(gI1,gI2))*CpbarFuFuhhPR(gI1,gI2))*MFu(gI2)));

   return result * oneOver16PiSqr;
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

   result += AbsSqr(CpbargWpCgWpCVP())*B00(Sqr(p),Sqr(MVWp),Sqr(MVWp));
   result += AbsSqr(CpbargWpgWpVP())*B00(Sqr(p),Sqr(MVWp),Sqr(MVWp));
   result += -(A0(Sqr(MVWp))*(CpconjVWpVPVPVWp1() + CpconjVWpVPVPVWp2() + 4*
      CpconjVWpVPVPVWp3()));
   result += -4*AbsSqr(CpHpconjHpVP())*B00(Sqr(p),Sqr(MHp),Sqr(MHp));
   result += A0(Sqr(MHp))*CpHpconjHpVPVP();
   result += 2*AbsSqr(CpHpconjVWpVP())*B0(Sqr(p),Sqr(MVWp),Sqr(MHp));
   result += 2*CpconjVWpVPVPVWp3()*Sqr(MVWp);
   result += -0.6666666666666666*AbsSqr(CpconjVWpVPVWp())*(3*A0(Sqr(MVWp)) + 15*
      B00(Sqr(p),Sqr(MVWp),Sqr(MVWp)) - 6*Sqr(MVWp) + Sqr(p) + 3*B0(Sqr(p),Sqr(
      MVWp),Sqr(MVWp))*(Sqr(MVWp) + 2*Sqr(p)));
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

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_VZ_1loop(double p ) const
{
   std::complex<double> result;

   result += 0.5*A0(Sqr(MAh))*CpAhAhVZVZ();
   result += -4*AbsSqr(CpAhhhVZ())*B00(Sqr(p),Sqr(MAh),Sqr(Mhh));
   result += AbsSqr(CpbargWpCgWpCVZ())*B00(Sqr(p),Sqr(MVWp),Sqr(MVWp));
   result += AbsSqr(CpbargWpgWpVZ())*B00(Sqr(p),Sqr(MVWp),Sqr(MVWp));
   result += -(A0(Sqr(MVWp))*(4*CpconjVWpVWpVZVZ1() + CpconjVWpVWpVZVZ2() +
      CpconjVWpVWpVZVZ3()));
   result += 0.5*A0(Sqr(Mhh))*CphhhhVZVZ();
   result += AbsSqr(CphhVZVZ())*B0(Sqr(p),Sqr(MVZ),Sqr(Mhh));
   result += -4*AbsSqr(CpHpconjHpVZ())*B00(Sqr(p),Sqr(MHp),Sqr(MHp));
   result += A0(Sqr(MHp))*CpHpconjHpVZVZ();
   result += 2*AbsSqr(CpHpconjVWpVZ())*B0(Sqr(p),Sqr(MVWp),Sqr(MHp));
   result += 2*CpconjVWpVWpVZVZ1()*Sqr(MVWp);
   result += -0.6666666666666666*AbsSqr(CpconjVWpVWpVZ())*(3*A0(Sqr(MVWp)) + 15*
      B00(Sqr(p),Sqr(MVWp),Sqr(MVWp)) - 6*Sqr(MVWp) + Sqr(p) + 3*B0(Sqr(p),Sqr(
      MVWp),Sqr(MVWp))*(Sqr(MVWp) + 2*Sqr(p)));
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

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_VWp_1loop(double p ) const
{
   std::complex<double> result;

   result += 0.5*A0(Sqr(MAh))*CpAhAhconjVWpVWp();
   result += -4*AbsSqr(CpAhHpconjVWp())*B00(Sqr(p),Sqr(MAh),Sqr(MHp));
   result += AbsSqr(CpbargPgWpconjVWp())*B00(Sqr(p),Sqr(MVWp),Sqr(MVP));
   result += AbsSqr(CpbargWpCgPconjVWp())*B00(Sqr(p),Sqr(MVP),Sqr(MVWp));
   result += AbsSqr(CpbargWpCgZconjVWp())*B00(Sqr(p),Sqr(MVZ),Sqr(MVWp));
   result += AbsSqr(CpbargZgWpconjVWp())*B00(Sqr(p),Sqr(MVWp),Sqr(MVZ));
   result += -(A0(Sqr(MVWp))*(CpconjVWpconjVWpVWpVWp1() + 4*
      CpconjVWpconjVWpVWpVWp2() + CpconjVWpconjVWpVWpVWp3()));
   result += 0;
   result += AbsSqr(CphhconjVWpVWp())*B0(Sqr(p),Sqr(MVWp),Sqr(Mhh));
   result += 0.5*A0(Sqr(Mhh))*CphhhhconjVWpVWp();
   result += -4*AbsSqr(CphhHpconjVWp())*B00(Sqr(p),Sqr(Mhh),Sqr(MHp));
   result += A0(Sqr(MHp))*CpHpconjHpconjVWpVWp();
   result += AbsSqr(CpHpconjVWpVP())*B0(Sqr(p),0,Sqr(MHp));
   result += AbsSqr(CpHpconjVWpVZ())*B0(Sqr(p),Sqr(MVZ),Sqr(MHp));
   result += 2*CpconjVWpconjVWpVWpVWp2()*Sqr(MVWp);
   result += -(AbsSqr(CpconjVWpVPVWp())*(A0(Sqr(MVWp)) + 10*B00(Sqr(p),Sqr(MVWp),0
      ) - 2*Sqr(MVWp) + 0.6666666666666666*Sqr(p) + B0(Sqr(p),Sqr(MVWp),0)*(Sqr(
      MVWp) + 4*Sqr(p))));
   result += -0.5*A0(Sqr(MVZ))*(4*CpconjVWpVWpVZVZ1() + CpconjVWpVWpVZVZ2() +
      CpconjVWpVWpVZVZ3()) + CpconjVWpVWpVZVZ1()*Sqr(MVZ);
   result += -(AbsSqr(CpconjVWpVWpVZ())*(A0(Sqr(MVWp)) + A0(Sqr(MVZ)) + 10*B00(Sqr
      (p),Sqr(MVZ),Sqr(MVWp)) - 2*(Sqr(MVWp) + Sqr(MVZ) - 0.3333333333333333*Sqr(p
      )) + B0(Sqr(p),Sqr(MVZ),Sqr(MVWp))*(Sqr(MVWp) + Sqr(MVZ) + 4*Sqr(p))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFuconjVWpPL(gI1,gI2)) +
      AbsSqr(CpbarFdFuconjVWpPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFu(gI2)))
      + 4*B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFu(gI2)))*MFd(gI1)*MFu(gI2)*Re(Conj(
      CpbarFdFuconjVWpPL(gI1,gI2))*CpbarFdFuconjVWpPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFvconjVWpPL(gI1,gI2)) + AbsSqr
      (CpbarFeFvconjVWpPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFv(gI2))) + 4*B0
      (Sqr(p),Sqr(MFe(gI1)),Sqr(MFv(gI2)))*MFe(gI1)*MFv(gI2)*Re(Conj(
      CpbarFeFvconjVWpPL(gI1,gI2))*CpbarFeFvconjVWpPR(gI1,gI2))));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_Fd_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh))*Conj(CpbarUFdFdAhPL(gO2
      ,gI1))*CpbarUFdFdAhPR(gO1,gI1)*MFd(gI1));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh))*Conj(CpbarUFdFdhhPL(gO2
      ,gI2))*CpbarUFdFdhhPR(gO1,gI2)*MFd(gI2));
   result += -5.333333333333333*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),0))*
      Conj(CpbarUFdFdVGPR(gO2,gI2))*CpbarUFdFdVGPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),0))*Conj(
      CpbarUFdFdVPPR(gO2,gI2))*CpbarUFdFdVPPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFdFdVZPR(gO2,gI2))*CpbarUFdFdVZPL(gO1,gI2)*MFd(gI2));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MHp))*Conj(CpbarUFdFuconjHpPL
      (gO2,gI2))*CpbarUFdFuconjHpPR(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFdFuconjVWpPR(gO2,gI2))*CpbarUFdFuconjVWpPL(gO1,gI2)*MFu(gI2));

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

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh))*Conj(
      CpbarUFdFdAhPR(gO2,gI1))*CpbarUFdFdAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh))*Conj(
      CpbarUFdFdhhPR(gO2,gI2))*CpbarUFdFdhhPR(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),0))*
      Conj(CpbarUFdFdVGPL(gO2,gI2))*CpbarUFdFdVGPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),0))*Conj(CpbarUFdFdVPPL(
      gO2,gI2))*CpbarUFdFdVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFdFdVZPL(gO2,gI2))*CpbarUFdFdVZPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHp))*Conj(
      CpbarUFdFuconjHpPR(gO2,gI2))*CpbarUFdFuconjHpPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFdFuconjVWpPL(gO2,gI2))*CpbarUFdFuconjVWpPL(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh))*Conj(
      CpbarUFdFdAhPL(gO2,gI1))*CpbarUFdFdAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh))*Conj(
      CpbarUFdFdhhPL(gO2,gI2))*CpbarUFdFdhhPL(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),0))*
      Conj(CpbarUFdFdVGPR(gO2,gI2))*CpbarUFdFdVGPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),0))*Conj(CpbarUFdFdVPPR(
      gO2,gI2))*CpbarUFdFdVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFdFdVZPR(gO2,gI2))*CpbarUFdFdVZPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHp))*Conj(
      CpbarUFdFuconjHpPL(gO2,gI2))*CpbarUFdFuconjHpPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFdFuconjVWpPR(gO2,gI2))*CpbarUFdFuconjVWpPR(gO1,gI2));

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

   result += SUM(gI1,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh))*Conj(CpbarUFuFuAhPL(gO2
      ,gI1))*CpbarUFuFuAhPR(gO1,gI1)*MFu(gI1));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MHp))*Conj(CpbarUFuFdHpPL(gO2
      ,gI2))*CpbarUFuFdHpPR(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFuFdVWpPR(gO2,gI2))*CpbarUFuFdVWpPL(gO1,gI2)*MFd(gI2));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh))*Conj(CpbarUFuFuhhPL(gO2
      ,gI2))*CpbarUFuFuhhPR(gO1,gI2)*MFu(gI2));
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

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh))*Conj(
      CpbarUFuFuAhPR(gO2,gI1))*CpbarUFuFuAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHp))*Conj(
      CpbarUFuFdHpPR(gO2,gI2))*CpbarUFuFdHpPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFuFdVWpPL(gO2,gI2))*CpbarUFuFdVWpPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh))*Conj(
      CpbarUFuFuhhPR(gO2,gI2))*CpbarUFuFuhhPR(gO1,gI2));
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

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh))*Conj(
      CpbarUFuFuAhPL(gO2,gI1))*CpbarUFuFuAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHp))*Conj(
      CpbarUFuFdHpPL(gO2,gI2))*CpbarUFuFdHpPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFuFdVWpPR(gO2,gI2))*CpbarUFuFdVWpPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh))*Conj(
      CpbarUFuFuhhPL(gO2,gI2))*CpbarUFuFuhhPL(gO1,gI2));
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

   result += SUM(gI1,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh))*Conj(CpbarUFeFeAhPL(gO2
      ,gI1))*CpbarUFeFeAhPR(gO1,gI1)*MFe(gI1));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh))*Conj(CpbarUFeFehhPL(gO2
      ,gI2))*CpbarUFeFehhPR(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFe(gI2)),0))*Conj(
      CpbarUFeFeVPPR(gO2,gI2))*CpbarUFeFeVPPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFeFeVZPR(gO2,gI2))*CpbarUFeFeVZPL(gO1,gI2)*MFe(gI2));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MHp))*Conj(CpbarUFeFvconjHpPL
      (gO2,gI2))*CpbarUFeFvconjHpPR(gO1,gI2)*MFv(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFeFvconjVWpPR(gO2,gI2))*CpbarUFeFvconjVWpPL(gO1,gI2)*MFv(gI2));

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

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh))*Conj(
      CpbarUFeFeAhPR(gO2,gI1))*CpbarUFeFeAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh))*Conj(
      CpbarUFeFehhPR(gO2,gI2))*CpbarUFeFehhPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),0))*Conj(CpbarUFeFeVPPL(
      gO2,gI2))*CpbarUFeFeVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFeFeVZPL(gO2,gI2))*CpbarUFeFeVZPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHp))*Conj(
      CpbarUFeFvconjHpPR(gO2,gI2))*CpbarUFeFvconjHpPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFeFvconjVWpPL(gO2,gI2))*CpbarUFeFvconjVWpPL(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh))*Conj(
      CpbarUFeFeAhPL(gO2,gI1))*CpbarUFeFeAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh))*Conj(
      CpbarUFeFehhPL(gO2,gI2))*CpbarUFeFehhPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),0))*Conj(CpbarUFeFeVPPR(
      gO2,gI2))*CpbarUFeFeVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFeFeVZPR(gO2,gI2))*CpbarUFeFeVZPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHp))*Conj(
      CpbarUFeFvconjHpPL(gO2,gI2))*CpbarUFeFvconjHpPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFeFvconjVWpPR(gO2,gI2))*CpbarUFeFvconjVWpPR(gO1,gI2));

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

   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MHp))*Conj(CpbarFvFeHpPL(gO2,
      gI2))*CpbarFvFeHpPR(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWp)))*Conj(
      CpbarFvFeVWpPR(gO2,gI2))*CpbarFvFeVWpPL(gO1,gI2)*MFe(gI2));
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

   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MHp))*Conj(CpbarFvFeHpPR
      (gO2,gI2))*CpbarFvFeHpPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWp)))*Conj(
      CpbarFvFeVWpPL(gO2,gI2))*CpbarFvFeVWpPL(gO1,gI2));
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

   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MHp))*Conj(CpbarFvFeHpPL
      (gO2,gI2))*CpbarFvFeHpPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWp)))*Conj(
      CpbarFvFeVWpPR(gO2,gI2))*CpbarFvFeVWpPR(gO1,gI2));
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

   result += SUM(gI1,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh))*Conj(CpbarFdFdAhPL(gO2,
      gI1))*CpbarFdFdAhPR(gO1,gI1)*MFd(gI1));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh))*Conj(CpbarFdFdhhPL(gO2,
      gI2))*CpbarFdFdhhPR(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ)))*Conj(
      CpbarFdFdVZPR(gO2,gI2))*CpbarFdFdVZPL(gO1,gI2)*MFd(gI2));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MHp))*Conj(CpbarFdFuconjHpPL(
      gO2,gI2))*CpbarFdFuconjHpPR(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWp)))*Conj(
      CpbarFdFuconjVWpPR(gO2,gI2))*CpbarFdFuconjVWpPL(gO1,gI2)*MFu(gI2));

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

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh))*Conj(CpbarFdFdAhPR
      (gO2,gI1))*CpbarFdFdAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh))*Conj(CpbarFdFdhhPR
      (gO2,gI2))*CpbarFdFdhhPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ)))*Conj(
      CpbarFdFdVZPL(gO2,gI2))*CpbarFdFdVZPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHp))*Conj(
      CpbarFdFuconjHpPR(gO2,gI2))*CpbarFdFuconjHpPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWp)))*Conj(
      CpbarFdFuconjVWpPL(gO2,gI2))*CpbarFdFuconjVWpPL(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh))*Conj(CpbarFdFdAhPL
      (gO2,gI1))*CpbarFdFdAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh))*Conj(CpbarFdFdhhPL
      (gO2,gI2))*CpbarFdFdhhPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ)))*Conj(
      CpbarFdFdVZPR(gO2,gI2))*CpbarFdFdVZPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHp))*Conj(
      CpbarFdFuconjHpPL(gO2,gI2))*CpbarFdFuconjHpPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWp)))*Conj(
      CpbarFdFuconjVWpPR(gO2,gI2))*CpbarFdFuconjVWpPR(gO1,gI2));

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

   result += SUM(gI1,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh))*Conj(CpbarFeFeAhPL(gO2,
      gI1))*CpbarFeFeAhPR(gO1,gI1)*MFe(gI1));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh))*Conj(CpbarFeFehhPL(gO2,
      gI2))*CpbarFeFehhPR(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ)))*Conj(
      CpbarFeFeVZPR(gO2,gI2))*CpbarFeFeVZPL(gO1,gI2)*MFe(gI2));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MHp))*Conj(CpbarFeFvconjHpPL(
      gO2,gI2))*CpbarFeFvconjHpPR(gO1,gI2)*MFv(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWp)))*Conj(
      CpbarFeFvconjVWpPR(gO2,gI2))*CpbarFeFvconjVWpPL(gO1,gI2)*MFv(gI2));

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

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh))*Conj(CpbarFeFeAhPR
      (gO2,gI1))*CpbarFeFeAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh))*Conj(CpbarFeFehhPR
      (gO2,gI2))*CpbarFeFehhPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ)))*Conj(
      CpbarFeFeVZPL(gO2,gI2))*CpbarFeFeVZPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHp))*Conj(
      CpbarFeFvconjHpPR(gO2,gI2))*CpbarFeFvconjHpPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWp)))*Conj(
      CpbarFeFvconjVWpPL(gO2,gI2))*CpbarFeFvconjVWpPL(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh))*Conj(CpbarFeFeAhPL
      (gO2,gI1))*CpbarFeFeAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh))*Conj(CpbarFeFehhPL
      (gO2,gI2))*CpbarFeFehhPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ)))*Conj(
      CpbarFeFeVZPR(gO2,gI2))*CpbarFeFeVZPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHp))*Conj(
      CpbarFeFvconjHpPL(gO2,gI2))*CpbarFeFvconjHpPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWp)))*Conj(
      CpbarFeFvconjVWpPR(gO2,gI2))*CpbarFeFvconjVWpPR(gO1,gI2));

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

   result += SUM(gI1,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh))*Conj(CpbarFuFuAhPL(gO2,
      gI1))*CpbarFuFuAhPR(gO1,gI1)*MFu(gI1));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MHp))*Conj(CpbarFuFdHpPL(gO2,
      gI2))*CpbarFuFdHpPR(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWp)))*Conj(
      CpbarFuFdVWpPR(gO2,gI2))*CpbarFuFdVWpPL(gO1,gI2)*MFd(gI2));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh))*Conj(CpbarFuFuhhPL(gO2,
      gI2))*CpbarFuFuhhPR(gO1,gI2)*MFu(gI2));
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

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh))*Conj(CpbarFuFuAhPR
      (gO2,gI1))*CpbarFuFuAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHp))*Conj(CpbarFuFdHpPR
      (gO2,gI2))*CpbarFuFdHpPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWp)))*Conj(
      CpbarFuFdVWpPL(gO2,gI2))*CpbarFuFdVWpPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh))*Conj(CpbarFuFuhhPR
      (gO2,gI2))*CpbarFuFuhhPR(gO1,gI2));
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

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh))*Conj(CpbarFuFuAhPL
      (gO2,gI1))*CpbarFuFuAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHp))*Conj(CpbarFuFdHpPL
      (gO2,gI2))*CpbarFuFdHpPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWp)))*Conj(
      CpbarFuFdVWpPR(gO2,gI2))*CpbarFuFdVWpPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh))*Conj(CpbarFuFuhhPL
      (gO2,gI2))*CpbarFuFuhhPL(gO1,gI2));
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

   result += SUM(gI1,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh))*Conj(CpbarUFuFuAhPL(gO2
      ,gI1))*CpbarUFuFuAhPR(gO1,gI1)*MFu(gI1));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MHp))*Conj(CpbarUFuFdHpPL(gO2
      ,gI2))*CpbarUFuFdHpPR(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFuFdVWpPR(gO2,gI2))*CpbarUFuFdVWpPL(gO1,gI2)*MFd(gI2));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh))*Conj(CpbarUFuFuhhPL(gO2
      ,gI2))*CpbarUFuFuhhPR(gO1,gI2)*MFu(gI2));
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

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh))*Conj(
      CpbarUFuFuAhPR(gO2,gI1))*CpbarUFuFuAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHp))*Conj(
      CpbarUFuFdHpPR(gO2,gI2))*CpbarUFuFdHpPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFuFdVWpPL(gO2,gI2))*CpbarUFuFdVWpPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh))*Conj(
      CpbarUFuFuhhPR(gO2,gI2))*CpbarUFuFuhhPR(gO1,gI2));
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

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh))*Conj(
      CpbarUFuFuAhPL(gO2,gI1))*CpbarUFuFuAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHp))*Conj(
      CpbarUFuFdHpPL(gO2,gI2))*CpbarUFuFdHpPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFuFdVWpPR(gO2,gI2))*CpbarUFuFdVWpPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh))*Conj(
      CpbarUFuFuhhPL(gO2,gI2))*CpbarUFuFuhhPL(gO1,gI2));
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

std::complex<double> CLASSNAME::tadpole_hh_1loop() const
{
   std::complex<double> result;

   result += -0.5*A0(Sqr(MAh))*CpAhAhhh();
   result += A0(Sqr(MVWp))*CpbargWpCgWpChh();
   result += A0(Sqr(MVWp))*CpbargWpgWphh();
   result += A0(Sqr(MVZ))*CpbargZgZhh();
   result += -0.5*A0(Sqr(Mhh))*Cphhhhhh();
   result += -(A0(Sqr(MHp))*CphhHpconjHp());
   result += -0.5*A0(Sqr(Mss))*Cphhssss();
   result += 4*A0(Sqr(MVWp))*CphhconjVWpVWp() - 2*CphhconjVWpVWp()*Sqr(MVWp);
   result += CphhVZVZ()*(2*A0(Sqr(MVZ)) - Sqr(MVZ));
   result += 6*SUM(gI1,0,2,A0(Sqr(MFd(gI1)))*(CpbarFdFdhhPL(gI1,gI1) +
      CpbarFdFdhhPR(gI1,gI1))*MFd(gI1));
   result += 2*SUM(gI1,0,2,A0(Sqr(MFe(gI1)))*(CpbarFeFehhPL(gI1,gI1) +
      CpbarFeFehhPR(gI1,gI1))*MFe(gI1));
   result += 6*SUM(gI1,0,2,A0(Sqr(MFu(gI1)))*(CpbarFuFuhhPL(gI1,gI1) +
      CpbarFuFuhhPR(gI1,gI1))*MFu(gI1));

   return result * oneOver16PiSqr;
}













void CLASSNAME::calculate_MVG_pole()
{
   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(MVG) old_MVG(MVG), new_MVG(MVG);

   do {
      PHYSICAL(MVG) = 0.;

      new_MVG = PHYSICAL(MVG);
      diff = MaxRelDiff(new_MVG, old_MVG);
      old_MVG = new_MVG;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::VG);
   else
      problems.unflag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::VG);
}

void CLASSNAME::calculate_Mss_pole()
{
   if (!force_output && problems.is_running_tachyon(
      ScalarSingletZ2DMMhInputMsInput_info::ss))
      return;

   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(Mss) old_Mss(Mss), new_Mss(Mss);

   do {
      const double M_tree(get_mass_matrix_ss());
      const double p = old_Mss;
      double self_energy = Re(self_energy_ss_1loop(p));
      const double mass_sqr = M_tree - self_energy;

      PHYSICAL(Mss) = SignedAbsSqrt(mass_sqr);

      new_Mss = PHYSICAL(Mss);
      diff = MaxRelDiff(new_Mss, old_Mss);
      old_Mss = new_Mss;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::ss);
   else
      problems.unflag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::ss);
}

void CLASSNAME::calculate_MFv_pole()
{
   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(MFv) old_MFv(MFv), new_MFv(MFv);

   do {
      PHYSICAL(MFv).setConstant(0.);

      new_MFv = PHYSICAL(MFv);
      diff = MaxRelDiff(new_MFv, old_MFv);
      old_MFv = new_MFv;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::Fv);
   else
      problems.unflag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::Fv);
}

void CLASSNAME::calculate_Mhh_pole()
{
   if (!force_output && problems.is_running_tachyon(
      ScalarSingletZ2DMMhInputMsInput_info::hh))
      return;

   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(Mhh) old_Mhh(Mhh), new_Mhh(Mhh);

   do {
      const double M_tree(get_mass_matrix_hh());
      const double p = old_Mhh;
      double self_energy = Re(self_energy_hh_1loop(p));
      const double mass_sqr = M_tree - self_energy;

      PHYSICAL(Mhh) = SignedAbsSqrt(mass_sqr);

      new_Mhh = PHYSICAL(Mhh);
      diff = MaxRelDiff(new_Mhh, old_Mhh);
      old_Mhh = new_Mhh;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::hh);
   else
      problems.unflag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::hh);
}

void CLASSNAME::calculate_MVP_pole()
{
   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(MVP) old_MVP(MVP), new_MVP(MVP);

   do {
      PHYSICAL(MVP) = 0.;

      new_MVP = PHYSICAL(MVP);
      diff = MaxRelDiff(new_MVP, old_MVP);
      old_MVP = new_MVP;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::VP);
   else
      problems.unflag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::VP);
}

void CLASSNAME::calculate_MVZ_pole()
{
   if (!force_output && problems.is_running_tachyon(
      ScalarSingletZ2DMMhInputMsInput_info::VZ))
      return;

   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(MVZ) old_MVZ(MVZ), new_MVZ(MVZ);

   do {
      const double M_tree(Sqr(MVZ));
      const double p = old_MVZ;
      const double self_energy = Re(self_energy_VZ_1loop(p));
      const double mass_sqr = M_tree - self_energy;

      if (mass_sqr < 0.) {
         problems.flag_pole_tachyon(ScalarSingletZ2DMMhInputMsInput_info::VZ);
      }

      PHYSICAL(MVZ) = AbsSqrt(mass_sqr);

      new_MVZ = PHYSICAL(MVZ);
      diff = MaxRelDiff(new_MVZ, old_MVZ);
      old_MVZ = new_MVZ;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::VZ);
   else
      problems.unflag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::VZ);
}

void CLASSNAME::calculate_MFd_pole()
{
   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(MFd) old_MFd(MFd), new_MFd(MFd);

   do {
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fd());
      for (int es = 0; es < 3; ++es) {
         const double p = Abs(old_MFd(es));
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
         problems.flag_bad_mass(ScalarSingletZ2DMMhInputMsInput_info::Fd,
            eigenvalue_error > precision * Abs(eigen_values(0)));
      #else
         fs_svd(M_loop, eigen_values, mix_Vd, mix_Ud);
      #endif
         if (es == 0) {
            PHYSICAL(Vd) = mix_Vd;
            PHYSICAL(Ud) = mix_Ud;
         }
         PHYSICAL(MFd(es)) = Abs(eigen_values(es));
      }

      new_MFd = PHYSICAL(MFd);
      diff = MaxRelDiff(new_MFd, old_MFd);
      old_MFd = new_MFd;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::Fd);
   else
      problems.unflag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::Fd);
}

void CLASSNAME::calculate_MFu_pole()
{
   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(MFu) old_MFu(MFu), new_MFu(MFu);

   do {
      double qcd_1l = 0.;

      {
         const double currentScale = get_scale();
         qcd_1l = 0.008443431970194815*(-4. + 3.*Log(Sqr(MFu(2))/Sqr(
            currentScale)))*Sqr(g3);
      }

      double qcd_2l = 0.;

      if (pole_mass_loop_order > 1 && TOP_POLE_QCD_CORRECTION > 0) {
         const double currentScale = get_scale();
         qcd_2l = 2.2278607323533713e-6*Quad(g3)*(-2372.129769909197 + 1452.*
            Log(Sqr(MFu(2))/Sqr(currentScale)) - 396.*Sqr(Log(Sqr(MFu(2))/Sqr(
            currentScale))));
      }

      double qcd_3l = 0.;

      if (pole_mass_loop_order > 2 && TOP_POLE_QCD_CORRECTION > 1) {
         const double currentScale = get_scale();
         qcd_3l = 1.0450439184830051e-10*Power6(g3)*(-1.1451958668162137e7 +
            320760.*Cube(Log(Sqr(MFu(2))/Sqr(currentScale))) +
            4.935665892982448e6*Log(Sqr(MFu(2))/Sqr(currentScale)) - 1.7253e6*
            Sqr(Log(Sqr(MFu(2))/Sqr(currentScale))));
      }

      double qcd_4l = 0.;

      if (pole_mass_loop_order > 3 && TOP_POLE_QCD_CORRECTION > 2) {
         const double currentScale = get_scale();
         qcd_4l = -1.6081297554549208e-9*Power8(g3)*(211681.74421123447 - 5638.
            *Cube(Log(Sqr(MFu(2))/Sqr(currentScale))) - 104673.38261571848*Log(
            Sqr(MFu(2))/Sqr(currentScale)) + 825.*Quad(Log(Sqr(MFu(2))/Sqr(
            currentScale))) + 22162.91142653778*Sqr(Log(Sqr(MFu(2))/Sqr(
            currentScale))));
      }

      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fu());
      for (int es = 0; es < 3; ++es) {
         const double p = Abs(old_MFu(es));
         Eigen::Matrix<double,3,3> self_energy_1;
         Eigen::Matrix<double,3,3> self_energy_PL;
         Eigen::Matrix<double,3,3> self_energy_PR;
         for (int i1 = 0; i1 < 3; ++i1) {
            for (int i2 = 0; i2 < 3; ++i2) {
               if (i1 == 2 && i2 == 2) {
                  self_energy_1(i1,i2)  = Re(self_energy_Fu_1loop_1_heavy(p,i1,
                     i2));
                  self_energy_PL(i1,i2) = Re(self_energy_Fu_1loop_PL_heavy(p,i1
                     ,i2));
                  self_energy_PR(i1,i2) = Re(self_energy_Fu_1loop_PR_heavy(p,i1
                     ,i2));
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
         problems.flag_bad_mass(ScalarSingletZ2DMMhInputMsInput_info::Fu,
            eigenvalue_error > precision * Abs(eigen_values(0)));
      #else
         fs_svd(M_loop, eigen_values, mix_Vu, mix_Uu);
      #endif
         if (es == 0) {
            PHYSICAL(Vu) = mix_Vu;
            PHYSICAL(Uu) = mix_Uu;
         }
         PHYSICAL(MFu(es)) = Abs(eigen_values(es));
      }

      new_MFu = PHYSICAL(MFu);
      diff = MaxRelDiff(new_MFu, old_MFu);
      old_MFu = new_MFu;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::Fu);
   else
      problems.unflag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::Fu);
}

void CLASSNAME::calculate_MFe_pole()
{
   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(MFe) old_MFe(MFe), new_MFe(MFe);

   do {
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fe());
      for (int es = 0; es < 3; ++es) {
         const double p = Abs(old_MFe(es));
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
         problems.flag_bad_mass(ScalarSingletZ2DMMhInputMsInput_info::Fe,
            eigenvalue_error > precision * Abs(eigen_values(0)));
      #else
         fs_svd(M_loop, eigen_values, mix_Ve, mix_Ue);
      #endif
         if (es == 0) {
            PHYSICAL(Ve) = mix_Ve;
            PHYSICAL(Ue) = mix_Ue;
         }
         PHYSICAL(MFe(es)) = Abs(eigen_values(es));
      }

      new_MFe = PHYSICAL(MFe);
      diff = MaxRelDiff(new_MFe, old_MFe);
      old_MFe = new_MFe;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::Fe);
   else
      problems.unflag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::Fe);
}

void CLASSNAME::calculate_MVWp_pole()
{
   if (!force_output && problems.is_running_tachyon(
      ScalarSingletZ2DMMhInputMsInput_info::VWp))
      return;

   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(MVWp) old_MVWp(MVWp), new_MVWp(MVWp);

   do {
      const double M_tree(Sqr(MVWp));
      const double p = old_MVWp;
      const double self_energy = Re(self_energy_VWp_1loop(p));
      const double mass_sqr = M_tree - self_energy;

      if (mass_sqr < 0.) {
         problems.flag_pole_tachyon(ScalarSingletZ2DMMhInputMsInput_info::VWp);
      }

      PHYSICAL(MVWp) = AbsSqrt(mass_sqr);

      new_MVWp = PHYSICAL(MVWp);
      diff = MaxRelDiff(new_MVWp, old_MVWp);
      old_MVWp = new_MVWp;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::VWp);
   else
      problems.unflag_no_pole_mass_convergence(
         ScalarSingletZ2DMMhInputMsInput_info::VWp);
}

double CLASSNAME::calculate_MVWp_pole(double p)
{
   if (!force_output && problems.is_running_tachyon(
      ScalarSingletZ2DMMhInputMsInput_info::VWp))
      return 0.;

   const double self_energy = Re(self_energy_VWp_1loop(p));
   const double mass_sqr = Sqr(MVWp) - self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(ScalarSingletZ2DMMhInputMsInput_info::VWp);
   }

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVZ_pole(double p)
{
   if (!force_output && problems.is_running_tachyon(
      ScalarSingletZ2DMMhInputMsInput_info::VZ))
      return 0.;

   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = Sqr(MVZ) - self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(ScalarSingletZ2DMMhInputMsInput_info::VZ);
   }

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

double CLASSNAME::calculate_MVP_DRbar(double) const
{
   return 0.0;
}

double CLASSNAME::calculate_MVZ_DRbar(double m_pole) const
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(ScalarSingletZ2DMMhInputMsInput_info::VZ);
         return m_pole;
   }

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWp_DRbar(double m_pole) const
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VWp_1loop(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(ScalarSingletZ2DMMhInputMsInput_info::VWp);
         return m_pole;
   }

   return AbsSqrt(mass_sqr);
}



double CLASSNAME::ThetaW() const
{

   return ArcCos(Abs(ZZ(0,0)));
}

double CLASSNAME::VEV() const
{

   return v;
}



std::ostream& operator<<(std::ostream& ostr, const CLASSNAME& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
