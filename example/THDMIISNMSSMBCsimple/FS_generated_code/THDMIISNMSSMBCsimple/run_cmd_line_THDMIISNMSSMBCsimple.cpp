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

// File generated at Sat 11 Apr 2020 12:54:09

#include "config.h"

#include "THDMIISNMSSMBCsimple_input_parameters.hpp"
#include "THDMIISNMSSMBCsimple_observables.hpp"
#include "THDMIISNMSSMBCsimple_slha_io.hpp"
#include "THDMIISNMSSMBCsimple_spectrum_generator.hpp"

#ifdef ENABLE_TWO_SCALE_SOLVER
#include "THDMIISNMSSMBCsimple_two_scale_spectrum_generator.hpp"
#endif

#include "command_line_options.hpp"
#include "lowe.h"
#include "logger.hpp"
#include "physical_input.hpp"

#include <iostream>
#include <cstring>

namespace flexiblesusy {

void print_usage()
{
   std::cout <<
      "Usage: run_cmd_line_THDMIISNMSSMBCsimple.x [options]\n"
      "Options:\n"
      "  --lambdaNMSSM=<value>\n"
      "  --kappaNMSSM=<value>\n"
      "  --AlambdaNMSSM=<value>\n"
      "  --AkappaNMSSM=<value>\n"
      "  --AtopNMSSM=<value>\n"
      "  --mstopL=<value>\n"
      "  --mstopR=<value>\n"
      "  --MEWSB=<value>\n"
      "  --TanBeta=<value>\n"
      "  --vSIN=<value>\n"

      "  --solver-type=<value>             an integer corresponding\n"
      "                                    to the solver type to use\n"
      "  --help,-h                         print this help message"
             << std::endl;
}

void set_command_line_parameters(const Dynamic_array_view<char*>& args,
                                 THDMIISNMSSMBCsimple_input_parameters& input,
                                 int& solver_type)
{
   for (int i = 1; i < args.size(); ++i) {
      const auto option = args[i];

      if(Command_line_options::get_parameter_value(option, "--lambdaNMSSM=", input.lambdaNMSSM))
         continue;

      if(Command_line_options::get_parameter_value(option, "--kappaNMSSM=", input.kappaNMSSM))
         continue;

      if(Command_line_options::get_parameter_value(option, "--AlambdaNMSSM=", input.AlambdaNMSSM))
         continue;

      if(Command_line_options::get_parameter_value(option, "--AkappaNMSSM=", input.AkappaNMSSM))
         continue;

      if(Command_line_options::get_parameter_value(option, "--AtopNMSSM=", input.AtopNMSSM))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mstopL=", input.mstopL))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mstopR=", input.mstopR))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MEWSB=", input.MEWSB))
         continue;

      if(Command_line_options::get_parameter_value(option, "--TanBeta=", input.TanBeta))
         continue;

      if(Command_line_options::get_parameter_value(option, "--vSIN=", input.vSIN))
         continue;

      
      if (Command_line_options::get_parameter_value(
             option, "--solver-type=", solver_type))
         continue;

      if (strcmp(option,"--help") == 0 || strcmp(option,"-h") == 0) {
         print_usage();
         exit(EXIT_SUCCESS);
      }

      ERROR("Unrecognized command line option: " << option);
      exit(EXIT_FAILURE);
   }
}

template<class solver_type>
int run_solver(const THDMIISNMSSMBCsimple_input_parameters& input)
{
   Physical_input physical_input;
   softsusy::QedQcd qedqcd;

   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::precision, 1.0e-4);

   THDMIISNMSSMBCsimple_spectrum_generator<solver_type> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   auto model = std::get<0>(spectrum_generator.get_models_slha());

   THDMIISNMSSMBCsimple_scales scales;
   scales.HighScale = spectrum_generator.get_high_scale();
   scales.SUSYScale = spectrum_generator.get_susy_scale();
   scales.LowScale  = spectrum_generator.get_low_scale();
   scales.pole_mass_scale = spectrum_generator.get_pole_mass_scale();

   const auto observables = calculate_observables(
      model, qedqcd, physical_input, scales.pole_mass_scale);

   // SLHA output
   SLHAea::Coll slhaea(THDMIISNMSSMBCsimple_slha_io::fill_slhaea(
                          model, qedqcd, scales, observables));

   std::cout << slhaea;

   return spectrum_generator.get_exit_code();
}

int run(int solver_type, const THDMIISNMSSMBCsimple_input_parameters& input)
{
   int exit_code = 0;

   switch (solver_type) {
   case 0:
#ifdef ENABLE_TWO_SCALE_SOLVER
   case 1:
      exit_code = run_solver<Two_scale>(input);
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

   THDMIISNMSSMBCsimple_input_parameters input;
   int solver_type = 0;
   set_command_line_parameters(make_dynamic_array_view(&argv[0], argc), input,
                               solver_type);

   const int exit_code = run(solver_type, input);

   return exit_code;
}
