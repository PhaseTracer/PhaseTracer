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

// File generated at Sat 24 Oct 2020 17:07:56

#include "ScalarSingletZ2DM_two_scale_initial_guesser.hpp"
#include "ScalarSingletZ2DM_two_scale_model.hpp"
#include "ScalarSingletZ2DM_standard_model_two_scale_matching.hpp"
#include "standard_model_two_scale_model.hpp"
#include "lowe.h"
#include "error.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"

#include <Eigen/Core>

namespace flexiblesusy {

#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define SMPARAMETER(p) eft->get_##p()
#define PHASE(p) model->get_##p()
#define LowEnergyConstant(p) Electroweak_constants::p
#define MODEL model

ScalarSingletZ2DM_standard_model_initial_guesser<Two_scale>::ScalarSingletZ2DM_standard_model_initial_guesser(
   ScalarSingletZ2DM<Two_scale>* model_,
   standard_model::StandardModel<Two_scale>* eft_,
   const softsusy::QedQcd& qedqcd_,
   const standard_model::Standard_model_low_scale_constraint<Two_scale>& low_constraint_,
   const ScalarSingletZ2DM_susy_scale_constraint<Two_scale>& susy_constraint_,
   const ScalarSingletZ2DM_high_scale_constraint<Two_scale>& high_constraint_
)
   : Initial_guesser()
   , model(model_)
   , eft(eft_)
   , qedqcd(qedqcd_)
   , low_constraint(low_constraint_)
   , susy_constraint(susy_constraint_)
   , high_constraint(high_constraint_)
{
   if (!model)
      throw SetupError("ScalarSingletZ2DM_initial_guesser: Error: pointer to model"
                       " ScalarSingletZ2DM<Two_scale> must not be zero");
}

ScalarSingletZ2DM_standard_model_initial_guesser<Two_scale>::~ScalarSingletZ2DM_standard_model_initial_guesser()
{
}

/**
 * Guesses the DR-bar model parameters by calling
 * guess_susy_parameters() and guess_soft_parameters() .
 */
void ScalarSingletZ2DM_standard_model_initial_guesser<Two_scale>::guess()
{
   guess_eft_parameters();
   guess_model_parameters();
}

/**
 * Guesses the SUSY parameters (gauge, Yukawa couplings) at
 * \f$m_\text{top}^\text{pole}\f$ from the Standard Model gauge
 * couplings and fermion masses.  Threshold corrections are ignored.
 * The user-defined initial guess at the low-scale
 * (InitialGuessAtLowScale) is applied here:
 *
 * \code{.cpp}
   

 * \endcode
 */
void ScalarSingletZ2DM_standard_model_initial_guesser<Two_scale>::guess_eft_parameters()
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

   eft->set_g1(sqrt(4.0 * Pi * alpha_sm(0)));
   eft->set_g2(sqrt(4.0 * Pi * alpha_sm(1)));
   eft->set_g3(sqrt(4.0 * Pi * alpha_sm(2)));
   eft->set_scale(mtpole);

   eft->set_v(Electroweak_constants::vev);
   eft->set_Yu(ZEROMATRIX(3,3));
   eft->set_Yd(ZEROMATRIX(3,3));
   eft->set_Ye(ZEROMATRIX(3,3));

   eft->set_Yu(0, 0, -Sqrt(2.)* mu_guess/ eft->get_v());
   eft->set_Yu(1, 1, -Sqrt(2.)* mc_guess/ eft->get_v());
   eft->set_Yu(2, 2, -Sqrt(2.)* mt_guess/ eft->get_v());

   eft->set_Yd(0, 0, Sqrt(2.)* md_guess/ eft->get_v());
   eft->set_Yd(1, 1, Sqrt(2.)* ms_guess/ eft->get_v());
   eft->set_Yd(2, 2, Sqrt(2.)* mb_guess/ eft->get_v());

   eft->set_Ye(0, 0, Sqrt(2.)* me_guess/ eft->get_v());
   eft->set_Ye(1, 1, Sqrt(2.)* mm_guess/ eft->get_v());
   eft->set_Ye(2, 2, Sqrt(2.)* mtau_guess/ eft->get_v());

   eft->set_Lambdax(0.12604);
   eft->solve_ewsb_tree_level();
}

void ScalarSingletZ2DM_standard_model_initial_guesser<Two_scale>::calculate_DRbar_yukawa_couplings()
{
   calculate_running_SM_masses();
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
}

void ScalarSingletZ2DM_standard_model_initial_guesser<Two_scale>::calculate_running_SM_masses()
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
void ScalarSingletZ2DM_standard_model_initial_guesser<Two_scale>::calculate_Yu_DRbar()
{

}

/**
 * Calculates the Yukawa couplings Yd of the down-type
 * quarks from the Standard Model down-type quark masses (ignoring
 * threshold corrections).
 */
void ScalarSingletZ2DM_standard_model_initial_guesser<Two_scale>::calculate_Yd_DRbar()
{

}

/**
 * Calculates the Yukawa couplings Ye of the leptons
 * from the Standard Model down-type lepton masses (ignoring threshold
 * corrections).
 */
void ScalarSingletZ2DM_standard_model_initial_guesser<Two_scale>::calculate_Ye_DRbar()
{

}

/**
 * Guesses the full model parameters.  At first it runs the effective SM to the
 * guess of the high-scale (SUSYScaleFirstGuess), and matches gauge and Yukawas
 * to the full theory. Then, initial guess and constrints are imposed at this scale:
 * \code{.cpp}
   const auto LamSHInput = INPUTPARAMETER(LamSHInput);
   const auto LamSInput = INPUTPARAMETER(LamSInput);
   const auto muS2Input = INPUTPARAMETER(muS2Input);
   const auto muH2Input = INPUTPARAMETER(muH2Input);

   MODEL->set_v(Re(LowEnergyConstant(vev)));
   MODEL->set_LamSH(Re(LamSHInput));
   MODEL->set_LamS(Re(LamSInput));
   MODEL->set_muS2(Re(muS2Input));
   MODEL->set_muH2(Re(muH2Input));

 * \endcode
 * After that, the model is run to the high scale, initial guess and constraints
 * are imposed (HighScaleInput):
 *
 * \code{.cpp}
   

 * \endcode
 *
 * Afterwards, it runs to the SUSY-scale guess (SUSYScaleFirstGuess) and
 * solves the EWSB conditions at the tree-level.  Finally the DR-bar
 * mass spectrum is calculated.
 */
void ScalarSingletZ2DM_standard_model_initial_guesser<Two_scale>::guess_model_parameters()
{
   const double susy_scale_guess = susy_constraint.get_initial_scale_guess();
   const double high_scale_guess = high_constraint.get_initial_scale_guess();
   const auto scale_getter = [this] () { return susy_constraint.get_scale(); };

   model->set_scale(susy_scale_guess);

   // apply susy-scale first guess
   {
   const auto LamSHInput = INPUTPARAMETER(LamSHInput);
   const auto LamSInput = INPUTPARAMETER(LamSInput);
   const auto muS2Input = INPUTPARAMETER(muS2Input);
   const auto muH2Input = INPUTPARAMETER(muH2Input);

   MODEL->set_v(Re(LowEnergyConstant(vev)));
   MODEL->set_LamSH(Re(LamSHInput));
   MODEL->set_LamS(Re(LamSInput));
   MODEL->set_muS2(Re(muS2Input));
   MODEL->set_muH2(Re(muH2Input));

   }

   eft->run_to(susy_scale_guess, running_precision);
   eft->calculate_DRbar_masses();

   //get gauge and Yukawa couplings from effective theory
   ScalarSingletZ2DM_standard_model_matching_up<Two_scale> matching_up;
   matching_up.set_models(eft, model);
   matching_up.set_scale(scale_getter);
   matching_up.match_tree_level();

   model->run_to(susy_scale_guess, running_precision);

   // apply susy-scale constraint
   susy_constraint.set_model(model);
   susy_constraint.apply();

   // run to high scale
   model->run_to(high_scale_guess, running_precision);

   // apply user-defined initial guess at the high scale
   {
   

   }

   // apply high-scale constraint
   high_constraint.set_model(model);
   high_constraint.apply();

   model->run_to(susy_scale_guess, running_precision);

   // apply EWSB constraint
   model->solve_ewsb_tree_level();

   // calculate tree-level spectrum
   model->calculate_DRbar_masses();
}

} // namespace flexiblesusy
