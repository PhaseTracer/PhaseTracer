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
#include "ScalarSingletZ2DM_model_slha.hpp"
#include "ScalarSingletZ2DM_spectrum_generator.hpp"

#ifdef ENABLE_TWO_SCALE_SOLVER
#include "ScalarSingletZ2DM_two_scale_spectrum_generator.hpp"
#endif

#include "command_line_options.hpp"
#include "scan.hpp"
#include "lowe.h"
#include "logger.hpp"

#include <iostream>
#include <cstring>

#define INPUTPARAMETER(p) input.p

namespace flexiblesusy {

void print_usage()
{
   std::cout <<
      "Usage: scan_ScalarSingletZ2DM.x [options]\n"
      "Options:\n"
      "  --LamHInput=<value>\n"
      "  --LamSHInput=<value>\n"
      "  --LamSInput=<value>\n"
      "  --muS2Input=<value>\n"
      "  --QEWSB=<value>\n"
      "  --Qin=<value>\n"

      "  --solver-type=<value>             an integer corresponding\n"
      "                                    to the solver type to use\n"
      "  --help,-h                         print this help message"
             << std::endl;
}

void set_command_line_parameters(const Dynamic_array_view<char*>& args,
                                 ScalarSingletZ2DM_input_parameters& input,
                                 int& solver_type)
{
   for (int i = 1; i < args.size(); ++i) {
      const auto option = args[i];

      if(Command_line_options::get_parameter_value(option, "--LamHInput=", input.LamHInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--LamSHInput=", input.LamSHInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--LamSInput=", input.LamSInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--muS2Input=", input.muS2Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--QEWSB=", input.QEWSB))
         continue;

      if(Command_line_options::get_parameter_value(option, "--Qin=", input.Qin))
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

struct ScalarSingletZ2DM_scan_result {
   Spectrum_generator_problems problems;
   double higgs{0.};
};

template <class solver_type>
ScalarSingletZ2DM_scan_result run_parameter_point(const softsusy::QedQcd& qedqcd,
   ScalarSingletZ2DM_input_parameters& input)
{
   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::precision, 1.0e-4);

   ScalarSingletZ2DM_spectrum_generator<solver_type> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   const auto model = std::get<0>(spectrum_generator.get_models_slha());
   const auto& pole_masses = model.get_physical_slha();

   ScalarSingletZ2DM_scan_result result;
   result.problems = spectrum_generator.get_problems();
   result.higgs = pole_masses.Mhh;

   return result;
}

void scan(int solver_type, ScalarSingletZ2DM_input_parameters& input,
          const std::vector<double>& range)
{
   softsusy::QedQcd qedqcd;

   for (const auto p: range) {
      INPUTPARAMETER(LamHInput) = p;

      ScalarSingletZ2DM_scan_result result;
      switch (solver_type) {
      case 0:
#ifdef ENABLE_TWO_SCALE_SOLVER
      case 1:
         result = run_parameter_point<Two_scale>(qedqcd, input);
         if (!result.problems.have_problem() || solver_type != 0) break;
#endif

      default:
         if (solver_type != 0) {
            ERROR("unknown solver type: " << solver_type);
            exit(EXIT_FAILURE);
         }
      }

      const int error = result.problems.have_problem();
      std::cout << "  "
                << std::setw(12) << std::left << p << ' '
                << std::setw(12) << std::left << result.higgs << ' '
                << std::setw(12) << std::left << error;
      if (error) {
         std::cout << "\t# " << result.problems;
      }
      std::cout << '\n';
   }
}

} // namespace flexiblesusy


int main(int argc, char* argv[])
{
   using namespace flexiblesusy;

   ScalarSingletZ2DM_input_parameters input;
   int solver_type = 1;
   set_command_line_parameters(make_dynamic_array_view(&argv[0], argc), input,
                               solver_type);

   std::cout << "# "
             << std::setw(12) << std::left << "LamHInput" << ' '
             << std::setw(12) << std::left << "Mhh/GeV" << ' '
             << std::setw(12) << std::left << "error"
             << '\n';

   const std::vector<double> range(float_range(0., 100., 10));

   scan(solver_type, input, range);

   return 0;
}
