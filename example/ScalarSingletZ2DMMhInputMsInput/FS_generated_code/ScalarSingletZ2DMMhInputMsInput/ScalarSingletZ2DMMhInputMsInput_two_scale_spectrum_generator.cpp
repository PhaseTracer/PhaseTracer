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


#include "ScalarSingletZ2DMMhInputMsInput_two_scale_spectrum_generator.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_input_parameters.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_two_scale_convergence_tester.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_two_scale_ewsb_solver.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_two_scale_high_scale_constraint.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_two_scale_initial_guesser.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_two_scale_low_scale_constraint.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_two_scale_susy_scale_constraint.hpp"

#include "error.hpp"
#include "lowe.h"
#include "numerics2.hpp"
#include "two_scale_running_precision.hpp"
#include "two_scale_solver.hpp"

#include <limits>

namespace flexiblesusy {

double ScalarSingletZ2DMMhInputMsInput_spectrum_generator<Two_scale>::get_pole_mass_scale() const
{
   return settings.get(Spectrum_generator_settings::pole_mass_scale) != 0. ?
      settings.get(Spectrum_generator_settings::pole_mass_scale) :
      susy_scale;
}

/**
 * @brief Run's the RG solver with the given input parameters
 *
 * This function sets up the RG solver using a high-scale, susy-scale
 * and low-scale constraint.  Afterwards the solver is run until
 * convergence is reached or an error occours.  Finally the particle
 * spectrum (pole masses) is calculated.
 *
 * @param qedqcd Standard Model input parameters
 * @param input model input parameters
 */
void ScalarSingletZ2DMMhInputMsInput_spectrum_generator<Two_scale>::run_except(const softsusy::QedQcd& qedqcd,
                                const ScalarSingletZ2DMMhInputMsInput_input_parameters& input)
{
   VERBOSE_MSG("Solving BVP using two-scale solver");

   problems.set_bvp_solver_problems({ BVP_solver_problems("TwoScaleSolver") });

   auto& model = this->model;
   model.clear();
   model.set_input_parameters(input);
   model.do_calculate_sm_pole_masses(
      settings.get(Spectrum_generator_settings::calculate_sm_masses));
   model.do_calculate_bsm_pole_masses(
      settings.get(Spectrum_generator_settings::calculate_bsm_masses));
   model.do_force_output(
      settings.get(Spectrum_generator_settings::force_output));
   model.set_loops(
      settings.get(Spectrum_generator_settings::beta_loop_order));
   model.set_thresholds(
      settings.get(
         Spectrum_generator_settings::threshold_corrections_loop_order));
   model.set_zero_threshold(
      settings.get(Spectrum_generator_settings::beta_zero_threshold));

   ScalarSingletZ2DMMhInputMsInput_ewsb_solver<Two_scale> ewsb_solver;
   model.set_ewsb_solver(
      std::make_shared<ScalarSingletZ2DMMhInputMsInput_ewsb_solver<Two_scale> >(ewsb_solver));

   ScalarSingletZ2DMMhInputMsInput_high_scale_constraint<Two_scale> high_scale_constraint(&model);
   ScalarSingletZ2DMMhInputMsInput_susy_scale_constraint<Two_scale> susy_scale_constraint(&model, qedqcd);
   ScalarSingletZ2DMMhInputMsInput_low_scale_constraint<Two_scale>  low_scale_constraint(&model, qedqcd);

   high_scale_constraint.initialize();
   susy_scale_constraint.initialize();
   low_scale_constraint .initialize();

   low_scale_constraint.set_SM_like_Higgs_index(
      settings.get(Spectrum_generator_settings::eft_higgs_index));

   ScalarSingletZ2DMMhInputMsInput_convergence_tester<Two_scale> convergence_tester(
      &model, settings.get(Spectrum_generator_settings::precision));

   if (settings.get(Spectrum_generator_settings::pole_mass_scale) != 0.) {
      convergence_tester.set_scale_getter(
         [this](){ return settings.get(
               Spectrum_generator_settings::pole_mass_scale);
         });
   }

   if (settings.get(Spectrum_generator_settings::max_iterations) > 0) {
      convergence_tester.set_max_iterations(
         settings.get(Spectrum_generator_settings::max_iterations));
   }

   ScalarSingletZ2DMMhInputMsInput_initial_guesser<Two_scale> initial_guesser(&model, qedqcd,
                                                          low_scale_constraint,
                                                          susy_scale_constraint,
                                                          high_scale_constraint);

   Two_scale_increasing_precision precision(
      10.0, settings.get(Spectrum_generator_settings::precision));

   RGFlow<Two_scale> solver;
   solver.set_convergence_tester(&convergence_tester);
   solver.set_running_precision(&precision);
   solver.set_initial_guesser(&initial_guesser);
   solver.add(&low_scale_constraint, &model);
   solver.add(&high_scale_constraint, &model);
   solver.add(&susy_scale_constraint, &model);

   high_scale = susy_scale = low_scale = 0.;
   reached_precision = std::numeric_limits<double>::infinity();

   solver.solve();

   // impose low-scale constraint one last time
   model.run_to(low_scale_constraint.get_scale());
   low_scale_constraint.apply();

   high_scale = high_scale_constraint.get_scale();
   susy_scale = susy_scale_constraint.get_scale();
   low_scale  = low_scale_constraint.get_scale();
   reached_precision = convergence_tester.get_current_accuracy();

   calculate_spectrum();

   // copy calculated W pole mass
   model.get_physical().MVWp
      = low_scale_constraint.get_sm_parameters().displayPoleMW();

   // run to output scale (if scale > 0)
   if (!is_zero(parameter_output_scale))
      model.run_to(parameter_output_scale);
}

/**
 * Create a text file which contains the values of all model
 * parameters at all scales between the low-scale and the high-scale.
 *
 * @param filename name of output file
 */
void ScalarSingletZ2DMMhInputMsInput_spectrum_generator<Two_scale>::write_running_couplings(
   const std::string& filename) const
{
   ScalarSingletZ2DMMhInputMsInput_spectrum_generator_interface<Two_scale>::write_running_couplings(
      filename, get_low_scale(), get_high_scale());
}

void ScalarSingletZ2DMMhInputMsInput_spectrum_generator<Two_scale>::calculate_spectrum()
{
   model.run_to(get_pole_mass_scale());
   model.solve_ewsb();
   model.calculate_spectrum();
}

} // namespace flexiblesusy
