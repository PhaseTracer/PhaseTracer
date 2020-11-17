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

// File generated at Tue 17 Nov 2020 16:11:28

#include "config.h"

#include "ScalarSingletZ2DM_input_parameters.hpp"
#include "ScalarSingletZ2DM_observables.hpp"
#include "ScalarSingletZ2DM_slha_io.hpp"
#include "ScalarSingletZ2DM_spectrum_generator.hpp"
#include "ScalarSingletZ2DM_utilities.hpp"

#ifdef ENABLE_TWO_SCALE_SOLVER
#include "ScalarSingletZ2DM_two_scale_spectrum_generator.hpp"
#endif

#include "physical_input.hpp"
#include "spectrum_generator_settings.hpp"
#include "lowe.h"
#include "command_line_options.hpp"

#include <iostream>
#include <cstdlib>

/**
 * @brief Runs the spectrum generator of type \a solver_type
 * @tparam solver_type solver type
 * @param slha_io SLHA input
 * @param spectrum_generator_settings
 * @param slha_output_file output file for SLHA output
 * @param database_output_file output file for SQLite database
 * @param spectrum_file output file for the mass spectrum
 * @param rgflow_file output file for the RG flow
 * @return value of spectrum_generator::get_exit_code()
 */
template <class solver_type>
int run_solver(flexiblesusy::ScalarSingletZ2DM_slha_io& slha_io,
               const flexiblesusy::Spectrum_generator_settings& spectrum_generator_settings,
               const std::string& slha_output_file,
               const std::string& database_output_file,
               const std::string& spectrum_file,
               const std::string& rgflow_file)
{
   using namespace flexiblesusy;

   Physical_input physical_input; // extra non-SLHA physical input
   softsusy::QedQcd qedqcd;
   ScalarSingletZ2DM_input_parameters input;

   try {
      slha_io.fill(qedqcd);
      slha_io.fill(input);
      slha_io.fill(physical_input);
   } catch (const Error& error) {
      ERROR(error.what_detailed());
      return EXIT_FAILURE;
   }

   ScalarSingletZ2DM_spectrum_generator<solver_type> spectrum_generator;
   spectrum_generator.set_settings(spectrum_generator_settings);
   spectrum_generator.set_parameter_output_scale(
      slha_io.get_parameter_output_scale());

   spectrum_generator.run(qedqcd, input);

   auto models = spectrum_generator.get_models_slha();
   const auto& problems = spectrum_generator.get_problems();

   ScalarSingletZ2DM_scales scales;
   scales.HighScale = spectrum_generator.get_high_scale();
   scales.SUSYScale = spectrum_generator.get_susy_scale();
   scales.LowScale  = spectrum_generator.get_low_scale();
   scales.pole_mass_scale = spectrum_generator.get_pole_mass_scale();

   ScalarSingletZ2DM_observables observables;
   if (spectrum_generator_settings.get(Spectrum_generator_settings::calculate_observables))
      observables = calculate_observables(std::get<0>(models), qedqcd, physical_input, scales.pole_mass_scale);

   const bool show_result = !problems.have_problem() ||
      spectrum_generator_settings.get(Spectrum_generator_settings::force_output);
   // SLHA output
   if (!slha_output_file.empty()) {
      slha_io.set_spinfo(problems);
      slha_io.set_input(input);
      if (show_result) {
         slha_io.set_print_imaginary_parts_of_majorana_mixings(
            spectrum_generator_settings.get(
               Spectrum_generator_settings::force_positive_masses));
         slha_io.set_spectrum(models);
         slha_io.set_extra(std::get<0>(models), scales, observables);
      }

      slha_io.write_to(slha_output_file);
   }

   if (!database_output_file.empty() && show_result) {
      ScalarSingletZ2DM_database::to_database(
         database_output_file, std::get<0>(models), &qedqcd,
         &physical_input, &observables);
   }

   if (!spectrum_file.empty())
      spectrum_generator.write_spectrum(spectrum_file);

   if (!rgflow_file.empty())
      spectrum_generator.write_running_couplings(rgflow_file);

   return spectrum_generator.get_exit_code();
}

/**
 * @brief Runs the spectrum generator
 *
 * Reads the solver type from \a spectrum_generator_settings and calls
 * run_solver() with the corresponding solver type.
 *
 * @param slha_io SLHA input
 * @param spectrum_generator_settings
 * @param slha_output_file output file for SLHA output
 * @param database_output_file output file for SQLite database
 * @param spectrum_file output file for the mass spectrum
 * @param rgflow_file output file for the RG flow
 * @return return value of run_solver<>()
 */
int run(
   flexiblesusy::ScalarSingletZ2DM_slha_io& slha_io,
   const flexiblesusy::Spectrum_generator_settings& spectrum_generator_settings,
   const std::string& slha_output_file,
   const std::string& database_output_file,
   const std::string& spectrum_file,
   const std::string& rgflow_file)
{
   using namespace flexiblesusy;

   int exit_code = 0;
   const int solver_type
      = static_cast<int>(spectrum_generator_settings.get(
                            Spectrum_generator_settings::solver));

   switch (solver_type) {
   case 0:
#ifdef ENABLE_TWO_SCALE_SOLVER
   case 1:
      exit_code = run_solver<Two_scale>(
         slha_io, spectrum_generator_settings, slha_output_file,
         database_output_file, spectrum_file, rgflow_file);
      if (!exit_code || solver_type != 0) break;
#endif

   default:
      if (solver_type != 0) {
         ERROR("unknown solver type: " << solver_type);
         exit_code = -1;
      }
      break;
   }

   return exit_code;
}

int main(int argc, char* argv[])
{
   using namespace flexiblesusy;

   Command_line_options options(argc, argv);
   if (options.must_print_model_info())
      ScalarSingletZ2DM_info::print(std::cout);
   if (options.must_exit())
      return options.status();

   const std::string database_output_file(options.get_database_output_file());
   const std::string rgflow_file(options.get_rgflow_file());
   const std::string slha_input_source(options.get_slha_input_file());
   const std::string slha_output_file(options.get_slha_output_file());
   const std::string spectrum_file(options.get_spectrum_file());
   ScalarSingletZ2DM_slha_io slha_io;
   Spectrum_generator_settings spectrum_generator_settings;

   if (slha_input_source.empty()) {
      ERROR("No SLHA input source given!\n"
            "   Please provide one via the option --slha-input-file=");
      return EXIT_FAILURE;
   }

   try {
      slha_io.read_from_source(slha_input_source);
      slha_io.fill(spectrum_generator_settings);
   } catch (const Error& error) {
      ERROR(error.what_detailed());
      return EXIT_FAILURE;
   }

   const int exit_code
      = run(slha_io, spectrum_generator_settings, slha_output_file,
            database_output_file, spectrum_file, rgflow_file);

   return exit_code;
}
