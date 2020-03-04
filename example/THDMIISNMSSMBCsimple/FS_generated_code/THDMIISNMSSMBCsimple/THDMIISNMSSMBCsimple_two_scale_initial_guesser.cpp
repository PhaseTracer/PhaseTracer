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

// File generated at Thu 7 Nov 2019 18:53:56

#include "THDMIISNMSSMBCsimple_two_scale_initial_guesser.hpp"
#include "THDMIISNMSSMBCsimple_two_scale_model.hpp"
#include "lowe.h"
#include "error.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"

#include <Eigen/Core>

namespace flexiblesusy {

#define DERIVEDPARAMETER(p) model->p()
#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define PHASE(p) model->get_##p()
#define LowEnergyConstant(p) Electroweak_constants::p
#define MODEL model

THDMIISNMSSMBCsimple_initial_guesser<Two_scale>::THDMIISNMSSMBCsimple_initial_guesser(
   THDMIISNMSSMBCsimple<Two_scale>* model_,
   const softsusy::QedQcd& qedqcd_,
   const THDMIISNMSSMBCsimple_low_scale_constraint<Two_scale>& low_constraint_,
   const THDMIISNMSSMBCsimple_susy_scale_constraint<Two_scale>& susy_constraint_,
   const THDMIISNMSSMBCsimple_high_scale_constraint<Two_scale>& high_constraint_
)
   : model(model_)
   , qedqcd(qedqcd_)
   , low_constraint(low_constraint_)
   , susy_constraint(susy_constraint_)
   , high_constraint(high_constraint_)
{
   if (!model)
      throw SetupError("THDMIISNMSSMBCsimple_initial_guesser: Error: pointer to model"
                       " THDMIISNMSSMBCsimple<Two_scale> must not be zero");
}

/**
 * Guesses the DR-bar model parameters by calling
 * guess_susy_parameters() and guess_soft_parameters() .
 */
void THDMIISNMSSMBCsimple_initial_guesser<Two_scale>::guess()
{
   guess_susy_parameters();
   guess_soft_parameters();
}

/**
 * Guesses the SUSY parameters (gauge, Yukawa couplings) at
 * \f$m_\text{top}^\text{pole}\f$ from the Standard Model gauge
 * couplings and fermion masses.  Threshold corrections are ignored.
 * The user-defined initial guess at the low-scale
 * (InitialGuessAtLowScale) is applied here:
 *
 * \code{.cpp}
   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto vSIN = INPUTPARAMETER(vSIN);

   MODEL->set_v1(Re(LowEnergyConstant(vev)/Sqrt(1 + Sqr(TanBeta))));
   MODEL->set_v2(Re((TanBeta*LowEnergyConstant(vev))/Sqrt(1 + Sqr(TanBeta))));
   MODEL->set_vS(Re(vSIN));
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();

 * \endcode
 */
void THDMIISNMSSMBCsimple_initial_guesser<Two_scale>::guess_susy_parameters()
{
   softsusy::QedQcd leAtMt(qedqcd);
   const double mtpole = leAtMt.displayPoleMt();

   mu_guess = leAtMt.displayMass(softsusy::mUp);
   mc_guess = leAtMt.displayMass(softsusy::mCharm);
   mt_guess = model->get_thresholds() > 0 && model->get_threshold_corrections().mt > 0 ?
      leAtMt.displayMass(softsusy::mTop) - 30.0 :
      leAtMt.displayPoleMt();
   md_guess = leAtMt.displayMass(softsusy::mDown);
   ms_guess = leAtMt.displayMass(softsusy::mStrange);
   mb_guess = leAtMt.displayMass(softsusy::mBottom);
   me_guess = model->get_thresholds() > 0 ?
      leAtMt.displayMass(softsusy::mElectron) :
      leAtMt.displayPoleMel();
   mm_guess = model->get_thresholds() > 0 ?
      leAtMt.displayMass(softsusy::mMuon) :
      leAtMt.displayPoleMmuon();
   mtau_guess = leAtMt.displayMass(softsusy::mTau);

   calculate_running_SM_masses();

   // guess gauge couplings at mt
   const auto alpha_sm(leAtMt.guess_alpha_SM5(mtpole));

   MODEL->set_g1(Sqrt(4. * Pi * alpha_sm(0)));
   MODEL->set_g2(Sqrt(4. * Pi * alpha_sm(1)));
   MODEL->set_g3(Sqrt(4. * Pi * alpha_sm(2)));


   model->set_scale(mtpole);

   // apply user-defined initial guess at the low scale
   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto vSIN = INPUTPARAMETER(vSIN);

   MODEL->set_v1(Re(LowEnergyConstant(vev)/Sqrt(1 + Sqr(TanBeta))));
   MODEL->set_v2(Re((TanBeta*LowEnergyConstant(vev))/Sqrt(1 + Sqr(TanBeta))));
   MODEL->set_vS(Re(vSIN));
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();

}

void THDMIISNMSSMBCsimple_initial_guesser<Two_scale>::calculate_DRbar_yukawa_couplings()
{
   calculate_running_SM_masses();
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
}

void THDMIISNMSSMBCsimple_initial_guesser<Two_scale>::calculate_running_SM_masses()
{
   upQuarksDRbar.setZero();
   upQuarksDRbar(0,0) = mu_guess;
   upQuarksDRbar(1,1) = mc_guess;
   upQuarksDRbar(2,2) = mt_guess;

   downQuarksDRbar.setZero();
   downQuarksDRbar(0,0) = md_guess;
   downQuarksDRbar(1,1) = ms_guess;
   downQuarksDRbar(2,2) = mb_guess;

   downLeptonsDRbar.setZero();
   downLeptonsDRbar(0,0) = me_guess;
   downLeptonsDRbar(1,1) = mm_guess;
   downLeptonsDRbar(2,2) = mtau_guess;
}

/**
 * Calculates the Yukawa couplings Yu of the up-type quarks
 * from the Standard Model up-type quark masses (ignoring threshold
 * corrections).
 */
void THDMIISNMSSMBCsimple_initial_guesser<Two_scale>::calculate_Yu_DRbar()
{
   const auto v2 = MODELPARAMETER(v2);
   MODEL->set_Yu((-((1.4142135623730951*upQuarksDRbar)/v2).transpose()).real());

}

/**
 * Calculates the Yukawa couplings Yd of the down-type
 * quarks from the Standard Model down-type quark masses (ignoring
 * threshold corrections).
 */
void THDMIISNMSSMBCsimple_initial_guesser<Two_scale>::calculate_Yd_DRbar()
{
   const auto v1 = MODELPARAMETER(v1);
   MODEL->set_Yd((-((1.4142135623730951*downQuarksDRbar)/v1).transpose()).real());

}

/**
 * Calculates the Yukawa couplings Ye of the leptons
 * from the Standard Model down-type lepton masses (ignoring threshold
 * corrections).
 */
void THDMIISNMSSMBCsimple_initial_guesser<Two_scale>::calculate_Ye_DRbar()
{
   const auto v1 = MODELPARAMETER(v1);
   MODEL->set_Ye((-((1.4142135623730951*downLeptonsDRbar)/v1).transpose()).real())
      ;

}

/**
 * Guesses the soft-breaking parameters.  At first it runs to the
 * guess of the high-scale (HighScaleFirstGuess) and imposes the
 * high-scale constraint (HighScaleInput):
 *
 * \code{.cpp}
   

 * \endcode
 *
 * Afterwards, it runs to the low-scale guess (LowScaleFirstGuess) and
 * solves the EWSB conditions at the tree-level.  Finally the DR-bar
 * mass spectrum is calculated.
 */
void THDMIISNMSSMBCsimple_initial_guesser<Two_scale>::guess_soft_parameters()
{
   const double low_scale_guess = low_constraint.get_initial_scale_guess();
   const double high_scale_guess = high_constraint.get_initial_scale_guess();

   model->run_to(high_scale_guess, running_precision);

   // apply high-scale constraint
   high_constraint.set_model(model);
   high_constraint.apply();

   // apply user-defined initial guess at the high scale
   


   model->run_to(low_scale_guess, running_precision);

   // apply EWSB constraint
   model->solve_ewsb_tree_level();

   // calculate tree-level spectrum
   model->calculate_DRbar_masses();
}

} // namespace flexiblesusy
