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


#include "ScalarSingletZ2DMMhInput_two_scale_susy_scale_constraint.hpp"
#include "ScalarSingletZ2DMMhInput_two_scale_model.hpp"
#include "config.h"
#include "wrappers.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "ew_input.hpp"
#include "minimizer.hpp"
#include "raii.hpp"
#include "root_finder.hpp"
#include "threshold_loop_functions.hpp"

#include <cmath>

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
#define LowEnergyConstant(p) Electroweak_constants::p
#define MZPole qedqcd.displayPoleMZ()
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define SCALE model->get_scale()
#define THRESHOLD static_cast<int>(model->get_thresholds())
#define MODEL model
#define MODELCLASSNAME ScalarSingletZ2DMMhInput<Two_scale>

ScalarSingletZ2DMMhInput_susy_scale_constraint<Two_scale>::ScalarSingletZ2DMMhInput_susy_scale_constraint(
   ScalarSingletZ2DMMhInput<Two_scale>* model_, const softsusy::QedQcd& qedqcd_)
   : model(model_)
   , qedqcd(qedqcd_)
{
   initialize();
}

void ScalarSingletZ2DMMhInput_susy_scale_constraint<Two_scale>::apply()
{
   check_model_ptr();

   


   model->calculate_DRbar_masses();
   update_scale();

   // apply user-defined susy scale constraints
   

   {
      auto solver_0 = [&](const Eigen::Matrix<double,1,1>& x) {
         MODEL->set_LamH(x(0));

         MODEL->calculate_DRbar_masses();
         const auto MhInput = INPUTPARAMETER(MhInput);
         MODEL->calculate_Mhh_pole();

         Eigen::Matrix<double,1,1> f;
         f(0) = MhInput - Pole(Mhh);

         return f;
      };

      const auto LamH = MODELPARAMETER(LamH);

      Eigen::VectorXd start_point(1);
      start_point << LamH;
      Root_finder<1> root_finder(solver_0, 100, 1.0e-2);
      const int status = root_finder.find_root(start_point);
      VERBOSE_MSG("\troot finder status: " << gsl_strerror(status));
   }

   MODEL->solve_ewsb();

}

double ScalarSingletZ2DMMhInput_susy_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double ScalarSingletZ2DMMhInput_susy_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const ScalarSingletZ2DMMhInput_input_parameters& ScalarSingletZ2DMMhInput_susy_scale_constraint<Two_scale>::get_input_parameters() const
{
   check_model_ptr();

   return model->get_input();
}

ScalarSingletZ2DMMhInput<Two_scale>* ScalarSingletZ2DMMhInput_susy_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void ScalarSingletZ2DMMhInput_susy_scale_constraint<Two_scale>::set_model(Model* model_)
{
   model = cast_model<ScalarSingletZ2DMMhInput<Two_scale>*>(model_);
}

void ScalarSingletZ2DMMhInput_susy_scale_constraint<Two_scale>::set_sm_parameters(
   const softsusy::QedQcd& qedqcd_)
{
   qedqcd = qedqcd_;
}

const softsusy::QedQcd& ScalarSingletZ2DMMhInput_susy_scale_constraint<Two_scale>::get_sm_parameters() const
{
   return qedqcd;
}

void ScalarSingletZ2DMMhInput_susy_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = nullptr;
   qedqcd = softsusy::QedQcd();
}

void ScalarSingletZ2DMMhInput_susy_scale_constraint<Two_scale>::initialize()
{
   check_model_ptr();

   const auto QEWSB = INPUTPARAMETER(QEWSB);

   initial_scale_guess = QEWSB;

   scale = initial_scale_guess;
}

void ScalarSingletZ2DMMhInput_susy_scale_constraint<Two_scale>::update_scale()
{
   check_model_ptr();

   const auto QEWSB = INPUTPARAMETER(QEWSB);

   scale = QEWSB;


}

void ScalarSingletZ2DMMhInput_susy_scale_constraint<Two_scale>::check_model_ptr() const
{
   if (!model)
      throw SetupError("ScalarSingletZ2DMMhInput_susy_scale_constraint<Two_scale>: "
                       "model pointer is zero!");
}

} // namespace flexiblesusy
