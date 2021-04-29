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


#include "config.h"

#include "ScalarSingletZ2DMMhInputMsInput_input_parameters.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_observables.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_slha_io.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_spectrum_generator.hpp"

#ifdef ENABLE_TWO_SCALE_SOLVER
#include "ScalarSingletZ2DMMhInputMsInput_two_scale_spectrum_generator.hpp"
#endif

#include "array_view.hpp"
#include "command_line_options.hpp"
#include "lowe.h"
#include "logger.hpp"
#include "physical_input.hpp"

#include <iostream>

namespace flexiblesusy {

void print_usage()
{
   std::cout <<
      "Usage: run_cmd_line_ScalarSingletZ2DMMhInputMsInput.x [options]\n"
      "Options:\n"
      "  --MhInput=<value>\n"
      "  --LamSHInput=<value>\n"
      "  --LamSInput=<value>\n"
      "  --MsInput=<value>\n"
      "  --QEWSB=<value>\n"
      "  --Qin=<value>\n"

      "  --solver-type=<value>             an integer corresponding\n"
      "                                    to the solver type to use\n"
      "  --loop-library=<value>            an integer corresponding to the used\n"
      "                                    realization of loop library\n"
      "  --help,-h                         print this help message"
             << std::endl;
}

void set_command_line_parameters(const Dynamic_array_view<char*>& args,
                                 ScalarSingletZ2DMMhInputMsInput_input_parameters& input,
                                 int& solver_type,
                                 int& loop_library)
{
   for (int i = 1; i < args.size(); ++i) {
      const std::string option = args[i];

      if(Command_line_options::get_parameter_value(option, "--MhInput=", input.MhInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--LamSHInput=", input.LamSHInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--LamSInput=", input.LamSInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MsInput=", input.MsInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--QEWSB=", input.QEWSB))
         continue;

      if(Command_line_options::get_parameter_value(option, "--Qin=", input.Qin))
         continue;

      
      if (Command_line_options::get_parameter_value(
             option, "--solver-type=", solver_type))
         continue;

      if (Command_line_options::get_parameter_value(
             option, "--loop-library=", loop_library))
         continue;

      if (option == "--help" || option == "-h") {
         print_usage();
         exit(EXIT_SUCCESS);
      }

      ERROR("Unrecognized command line option: " << option);
      exit(EXIT_FAILURE);
   }
}

template<class solver_type>
int run_solver(int loop_library, const ScalarSingletZ2DMMhInputMsInput_input_parameters& input)
{
   Physical_input physical_input;
   softsusy::QedQcd qedqcd;

   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::precision, 1.0e-4);
   settings.set(Spectrum_generator_settings::loop_library, loop_library);

   ScalarSingletZ2DMMhInputMsInput_spectrum_generator<solver_type> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   ScalarSingletZ2DMMhInputMsInput_scales scales;
   scales.HighScale = spectrum_generator.get_high_scale();
   scales.SUSYScale = spectrum_generator.get_susy_scale();
   scales.LowScale  = spectrum_generator.get_low_scale();
   scales.pole_mass_scale = spectrum_generator.get_pole_mass_scale();

   auto models = spectrum_generator.get_models_slha();

   const auto observables = calculate_observables(
      std::get<0>(models), qedqcd, physical_input, scales.pole_mass_scale);

   // SLHA output
   ScalarSingletZ2DMMhInputMsInput_slha_io slha_io;
   slha_io.fill(models, qedqcd, scales, observables, settings);
   slha_io.write_to_stream(std::cout);

   return spectrum_generator.get_exit_code();
}

int run(int solver_type, int loop_library, const ScalarSingletZ2DMMhInputMsInput_input_parameters& input)
{
   int exit_code = 0;

   switch (solver_type) {
   case 0:
#ifdef ENABLE_TWO_SCALE_SOLVER
   case 1:
      exit_code = run_solver<Two_scale>(loop_library,input);
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

} // namespace flexiblesusy


int main(int argc, char* argv[])
{
   using namespace flexiblesusy;

   ScalarSingletZ2DMMhInputMsInput_input_parameters input;
   int solver_type = 0;
   int loop_library = 0;
   set_command_line_parameters(make_dynamic_array_view(&argv[0], argc), input,
                               solver_type,loop_library);

   const int exit_code = run(solver_type,loop_library,input);

   return exit_code;
}
