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

// File generated at Tue 17 Nov 2020 16:11:27

#include "ScalarSingletZ2DM_two_scale_low_scale_constraint.hpp"
#include "ScalarSingletZ2DM_two_scale_model.hpp"
#include "ScalarSingletZ2DM_info.hpp"
#include "ScalarSingletZ2DM_weinberg_angle.hpp"
#include "config.h"
#include "wrappers.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "raii.hpp"
#include "root_finder.hpp"
#include "threshold_loop_functions.hpp"


#include <algorithm>
#include <cmath>
#include <limits>

namespace flexiblesusy {

#define DERIVEDPARAMETER(p) model->p()
#define EXTRAPARAMETER(p) model->get_##p()
#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define PHASE(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETAPARAMETER1(l,p) beta_functions_##l##L.get_##p()
#define BETA(p) beta_##p
#define BETA1(l,p) beta_##l##L_##p
#define LowEnergyGaugeCoupling(i) new_g##i
#define LowEnergyConstant(p) Electroweak_constants::p
#define MZPole qedqcd.displayPoleMZ()
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define SCALE model->get_scale()
#define MODEL model
#define MODELCLASSNAME ScalarSingletZ2DM<Two_scale>
#define MWMSbar mW_run
#define MWDRbar mW_run
#define MZMSbar mZ_run
#define MZDRbar mZ_run
#define EDRbar e_run
#define EMSbar e_run
#define CKM ckm
#define PMNS pmns
#define THETAW theta_w
#define THRESHOLD static_cast<int>(model->get_thresholds())
#define ALPHA_EM_DRBAR alpha_em_drbar
#define CALCULATE_DRBAR_MASSES() model->calculate_DRbar_masses()

ScalarSingletZ2DM_low_scale_constraint<Two_scale>::ScalarSingletZ2DM_low_scale_constraint(
   ScalarSingletZ2DM<Two_scale>* model_, const softsusy::QedQcd& qedqcd_)
   : model(model_)
   , qedqcd(qedqcd_)
{
   initialize();
}

void ScalarSingletZ2DM_low_scale_constraint<Two_scale>::apply()
{
   check_model_ptr();

   


   model->calculate_DRbar_masses();
   update_scale();
   qedqcd.run_to(scale, 1.0e-5);
   calculate_DRbar_gauge_couplings();
   calculate_running_SM_masses();

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);

   MODEL->set_v(Re((2*MZMSbar)/Sqrt(0.6*Sqr(g1) + Sqr(g2))));
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
   MODEL->set_g1(new_g1);
   MODEL->set_g2(new_g2);
   MODEL->set_g3(new_g3);

}

const Eigen::Matrix<std::complex<double>,3,3>& ScalarSingletZ2DM_low_scale_constraint<Two_scale>::get_ckm()
{
   return ckm;
}

const Eigen::Matrix<std::complex<double>,3,3>& ScalarSingletZ2DM_low_scale_constraint<Two_scale>::get_pmns()
{
   return pmns;
}

double ScalarSingletZ2DM_low_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double ScalarSingletZ2DM_low_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

void ScalarSingletZ2DM_low_scale_constraint<Two_scale>::set_model(Model* model_)
{
   model = cast_model<ScalarSingletZ2DM<Two_scale>*>(model_);
}

void ScalarSingletZ2DM_low_scale_constraint<Two_scale>::set_sm_parameters(
   const softsusy::QedQcd& qedqcd_)
{
   qedqcd = qedqcd_;
}

const softsusy::QedQcd& ScalarSingletZ2DM_low_scale_constraint<Two_scale>::get_sm_parameters() const
{
   return qedqcd;
}

void ScalarSingletZ2DM_low_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = nullptr;
   qedqcd = softsusy::QedQcd();
   ckm.setIdentity();
   pmns.setIdentity();
   upQuarksDRbar.setZero();
   downQuarksDRbar.setZero();
   downLeptonsDRbar.setZero();
   neutrinoDRbar.setZero();
   mW_run = 0.;
   mZ_run = 0.;
   AlphaS = 0.;
   e_run = 0.;
   ThetaWDRbar = 0.;
   new_g1 = 0.;
   new_g2 = 0.;
   new_g3 = 0.;
}

void ScalarSingletZ2DM_low_scale_constraint<Two_scale>::initialize()
{
   check_model_ptr();

   initial_scale_guess = qedqcd.displayPoleMZ();

   scale = initial_scale_guess;

   ckm = qedqcd.get_complex_ckm();
   pmns = qedqcd.get_complex_pmns();
   upQuarksDRbar.setZero();
   downQuarksDRbar.setZero();
   downLeptonsDRbar.setZero();
   neutrinoDRbar.setZero();
   mW_run = 0.;
   mZ_run = 0.;
   AlphaS = 0.;
   e_run = 0.;
   ThetaWDRbar = 0.;
   new_g1 = 0.;
   new_g2 = 0.;
   new_g3 = 0.;
}

void ScalarSingletZ2DM_low_scale_constraint<Two_scale>::update_scale()
{
   check_model_ptr();

   scale = qedqcd.displayPoleMZ();


}

void ScalarSingletZ2DM_low_scale_constraint<Two_scale>::calculate_threshold_corrections()
{
   check_model_ptr();

   if (qedqcd.get_scale() != get_scale())
      throw SetupError("Error: low-energy data"
          " set is not defined at the same scale as the low-energy"
          " constraint.  You need to run the low-energy data set to this"
          " scale!");

   const double alpha_em = qedqcd.displayAlpha(softsusy::ALPHA);
   const double alpha_s  = qedqcd.displayAlpha(softsusy::ALPHAS);
   const double mw_pole  = qedqcd.displayPoleMW();
   const double mz_pole  = qedqcd.displayPoleMZ();

   double delta_alpha_em = 0.;
   double delta_alpha_s  = 0.;

   if (model->get_thresholds() && model->get_threshold_corrections().alpha_em > 0)
      delta_alpha_em = calculate_delta_alpha_em(alpha_em);

   if (model->get_thresholds() && model->get_threshold_corrections().alpha_s > 0)
      delta_alpha_s  = calculate_delta_alpha_s(alpha_s);

   const double alpha_em_drbar = alpha_em / (1.0 - delta_alpha_em);
   const double alpha_s_drbar  = alpha_s  / (1.0 - delta_alpha_s);
   const double e_drbar        = Sqrt(4.0 * Pi * alpha_em_drbar);

   // interface variables
   mZ_run = mz_pole;
   mW_run = mw_pole;

   if (model->get_thresholds() && model->get_threshold_corrections().mz > 0)
      mZ_run = model->calculate_MVZ_DRbar(mz_pole);

   if (model->get_thresholds() && model->get_threshold_corrections().mw > 0)
      mW_run = model->calculate_MVWp_DRbar(mw_pole);

   AlphaS = alpha_s_drbar;
   e_run = e_drbar;
   ThetaWDRbar = calculate_theta_w();
}

double ScalarSingletZ2DM_low_scale_constraint<Two_scale>::calculate_theta_w()
{
   check_model_ptr();

   double theta_w = std::asin(Electroweak_constants::sinThetaW);

   ScalarSingletZ2DM_weinberg_angle::Sm_parameters sm_pars;
   sm_pars.fermi_constant = qedqcd.displayFermiConstant();
   sm_pars.mw_pole = qedqcd.displayPoleMW();
   sm_pars.mz_pole = qedqcd.displayPoleMZ();
   sm_pars.mt_pole = qedqcd.displayPoleMt();
   sm_pars.alpha_s = calculate_alpha_s_SM5_at(qedqcd, qedqcd.displayPoleMt());

   const int number_of_iterations =
       std::max(20, static_cast<int>(std::abs(-log10(MODEL->get_precision()) * 10)
          ));

   ScalarSingletZ2DM_weinberg_angle weinberg(MODEL, sm_pars);
   weinberg.set_number_of_loops(MODEL->get_threshold_corrections().sin_theta_w);
   weinberg.set_number_of_iterations(number_of_iterations);

   try {
      const auto result = weinberg.calculate();
      THETAW = ArcSin(result.first);

      if (MODEL->get_thresholds() && MODEL->get_threshold_corrections().
         sin_theta_w > 0)
         qedqcd.setPoleMW(result.second);

      MODEL->get_problems().unflag_no_sinThetaW_convergence();
   } catch (const Error& e) {
      VERBOSE_MSG(e.what_detailed());
      MODEL->get_problems().flag_no_sinThetaW_convergence();
   }

   return theta_w;
}

void ScalarSingletZ2DM_low_scale_constraint<Two_scale>::calculate_DRbar_gauge_couplings()
{
   check_model_ptr();
   calculate_threshold_corrections();

   new_g1 = 1.2909944487358056*EDRbar*Sec(ThetaWDRbar);
   new_g2 = EDRbar*Csc(ThetaWDRbar);
   new_g3 = 3.5449077018110318*Sqrt(AlphaS);

   if (IsFinite(new_g1)) {
      model->get_problems().unflag_non_perturbative_parameter(
         ScalarSingletZ2DM_info::g1);
   } else {
      model->get_problems().flag_non_perturbative_parameter(
         ScalarSingletZ2DM_info::g1, new_g1, get_scale());
      new_g1 = Electroweak_constants::g1;
   }

   if (IsFinite(new_g2)) {
      model->get_problems().unflag_non_perturbative_parameter(
         ScalarSingletZ2DM_info::g2);
   } else {
      model->get_problems().flag_non_perturbative_parameter(
         ScalarSingletZ2DM_info::g2, new_g2, get_scale());
      new_g2 = Electroweak_constants::g2;
   }
}

double ScalarSingletZ2DM_low_scale_constraint<Two_scale>::calculate_delta_alpha_em(double alphaEm) const
{
   check_model_ptr();

   const double currentScale = model->get_scale();
   const auto MFu = MODELPARAMETER(MFu);

   const double delta_alpha_em_SM = -0.28294212105225836*alphaEm*FiniteLog(Abs(MFu
      (2)/currentScale));

   const double delta_alpha_em = 0;

   return delta_alpha_em + delta_alpha_em_SM;

}

double ScalarSingletZ2DM_low_scale_constraint<Two_scale>::calculate_delta_alpha_s(double alphaS) const
{
   check_model_ptr();

   const double currentScale = model->get_scale();
   const auto MFu = MODELPARAMETER(MFu);

   const double delta_alpha_s_SM = -0.1061032953945969*alphaS*FiniteLog(Abs(MFu(2)
      /currentScale));

   const double delta_alpha_s = 0;

   const double delta_alpha_s_1loop = delta_alpha_s + delta_alpha_s_SM;
   double delta_alpha_s_2loop = 0.;
   double delta_alpha_s_3loop = 0.;
   double delta_alpha_s_4loop = 0.;

   return delta_alpha_s_1loop + delta_alpha_s_2loop + delta_alpha_s_3loop +
      delta_alpha_s_4loop;

}

double ScalarSingletZ2DM_low_scale_constraint<Two_scale>::calculate_alpha_s_SM5_at(
   softsusy::QedQcd qedqcd_tmp, double scale) const
{
   qedqcd_tmp.run_to(scale); // running in SM(5)
   return qedqcd_tmp.displayAlpha(softsusy::ALPHAS);
}

void ScalarSingletZ2DM_low_scale_constraint<Two_scale>::calculate_running_SM_masses()
{
   check_model_ptr();

   upQuarksDRbar.setZero();
   upQuarksDRbar(0,0) = qedqcd.displayMass(softsusy::mUp);
   upQuarksDRbar(1,1) = qedqcd.displayMass(softsusy::mCharm);
   upQuarksDRbar(2,2) = qedqcd.displayPoleMt();

   downQuarksDRbar.setZero();
   downQuarksDRbar(0,0) = qedqcd.displayMass(softsusy::mDown);
   downQuarksDRbar(1,1) = qedqcd.displayMass(softsusy::mStrange);
   downQuarksDRbar(2,2) = qedqcd.displayMass(softsusy::mBottom);

   downLeptonsDRbar.setZero();
   downLeptonsDRbar(0,0) = qedqcd.displayPoleMel();
   downLeptonsDRbar(1,1) = qedqcd.displayPoleMmuon();
   downLeptonsDRbar(2,2) = qedqcd.displayPoleMtau();

   neutrinoDRbar.setZero();
   neutrinoDRbar(0,0) = qedqcd.displayNeutrinoPoleMass(1);
   neutrinoDRbar(1,1) = qedqcd.displayNeutrinoPoleMass(2);
   neutrinoDRbar(2,2) = qedqcd.displayNeutrinoPoleMass(3);

   if (model->get_thresholds() && model->get_threshold_corrections().mt > 0) {
      upQuarksDRbar(2,2) = MODEL->calculate_MFu_DRbar(qedqcd.displayPoleMt(), 2);
   }

   if (model->get_thresholds() && model->get_threshold_corrections().mb > 0) {
      downQuarksDRbar(2,2) = MODEL->calculate_MFd_DRbar(qedqcd.displayMass(softsusy::mBottom), 2);
   }

   if (model->get_thresholds()) {
      downLeptonsDRbar(0,0) = MODEL->calculate_MFe_DRbar(qedqcd.displayMass(softsusy::mElectron), 0);
      downLeptonsDRbar(1,1) = MODEL->calculate_MFe_DRbar(qedqcd.displayMass(softsusy::mMuon), 1);
   }

   if (model->get_thresholds() && model->get_threshold_corrections().mtau > 0) {
      downLeptonsDRbar(2,2) = MODEL->calculate_MFe_DRbar(qedqcd.displayMass(softsusy::mTau), 2);
   }
}

void ScalarSingletZ2DM_low_scale_constraint<Two_scale>::calculate_DRbar_yukawa_couplings()
{
   calculate_running_SM_masses();
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
}

void ScalarSingletZ2DM_low_scale_constraint<Two_scale>::calculate_Yu_DRbar()
{
   check_model_ptr();

   const auto v = MODELPARAMETER(v);
   MODEL->set_Yu((-((1.4142135623730951*upQuarksDRbar)/v).transpose()).real());

}

void ScalarSingletZ2DM_low_scale_constraint<Two_scale>::calculate_Yd_DRbar()
{
   check_model_ptr();

   const auto v = MODELPARAMETER(v);
   MODEL->set_Yd((((1.4142135623730951*downQuarksDRbar)/v).transpose()).real());

}

void ScalarSingletZ2DM_low_scale_constraint<Two_scale>::calculate_Ye_DRbar()
{
   check_model_ptr();

   const auto v = MODELPARAMETER(v);
   MODEL->set_Ye((((1.4142135623730951*downLeptonsDRbar)/v).transpose()).real());

}

void ScalarSingletZ2DM_low_scale_constraint<Two_scale>::check_model_ptr() const
{
   if (!model)
      throw SetupError("ScalarSingletZ2DM_low_scale_constraint<Two_scale>: "
                       "model pointer is zero!");
}

} // namespace flexiblesusy
