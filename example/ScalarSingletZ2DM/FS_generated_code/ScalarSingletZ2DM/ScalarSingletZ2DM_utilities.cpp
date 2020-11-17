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

// File generated at Tue 17 Nov 2020 16:11:26

#include "ScalarSingletZ2DM_utilities.hpp"
#include "ScalarSingletZ2DM_input_parameters.hpp"
#include "ScalarSingletZ2DM_observables.hpp"
#include "error.hpp"
#include "logger.hpp"
#include "physical_input.hpp"
#include "database.hpp"
#include "wrappers.hpp"
#include "lowe.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <iterator>

#define PHYSICAL(p) model.get_physical().p
#define MODELPARAMETER(p) model.get_##p()

namespace flexiblesusy {

namespace utilities {

template <typename Iterable>
void append(std::vector<std::string>& a, const Iterable& b)
{
   a.insert(a.end(), begin(b), end(b));
}

void append(Eigen::ArrayXd& a, const Eigen::ArrayXd& b)
{
   const auto a_rows = a.rows();
   a.conservativeResize(a_rows + b.rows());
   a.block(a_rows, 0, b.rows(), 1) = b;
}

} // namespace utilities

namespace {

std::valarray<double> to_valarray(double v)
{
   return std::valarray<double>(&v, 1);
}

template <class Scalar, int M, int N>
std::valarray<double> to_valarray(const Eigen::Array<Scalar, M, N>& v)
{
   return std::valarray<double>(v.data(), v.size());
}

} // anonymous namespace

ScalarSingletZ2DM_spectrum_plotter::ScalarSingletZ2DM_spectrum_plotter(const ScalarSingletZ2DM_mass_eigenstates& model)
{
   extract_spectrum(model);
}

void ScalarSingletZ2DM_spectrum_plotter::extract_spectrum(const ScalarSingletZ2DM_mass_eigenstates& model)
{
   spectrum.clear();
   scale = model.get_scale();

   spectrum.emplace_back(TParticle("ss", "S", to_valarray(PHYSICAL(Mss))));

   if (model.do_calculate_sm_pole_masses()) {
      spectrum.emplace_back(TParticle("Ah", "A^0", to_valarray(PHYSICAL(MAh))));
      spectrum.emplace_back(TParticle("Fd", "d", to_valarray(PHYSICAL(MFd))));
      spectrum.emplace_back(TParticle("Fe", "e", to_valarray(PHYSICAL(MFe))));
      spectrum.emplace_back(TParticle("Fu", "u", to_valarray(PHYSICAL(MFu))));
      spectrum.emplace_back(TParticle("Fv", "\\nu", to_valarray(PHYSICAL(MFv))));
      spectrum.emplace_back(TParticle("hh", "h", to_valarray(PHYSICAL(Mhh))));
      spectrum.emplace_back(TParticle("Hp", "H^+", to_valarray(PHYSICAL(MHp))));
      spectrum.emplace_back(TParticle("VG", "g", to_valarray(PHYSICAL(MVG))));
      spectrum.emplace_back(TParticle("VP", "\\gamma", to_valarray(PHYSICAL(MVP))));
      spectrum.emplace_back(TParticle("VWp", "W^+", to_valarray(PHYSICAL(MVWp))));
      spectrum.emplace_back(TParticle("VZ", "Z", to_valarray(PHYSICAL(MVZ))));

   }
}

void ScalarSingletZ2DM_spectrum_plotter::write_to_file(const std::string& file_name) const
{
   if (spectrum.empty())
      return;

   std::ofstream filestr(file_name, std::ios::out);
   VERBOSE_MSG("ScalarSingletZ2DM_spectrum_plotter::write_to_file: opening file: "
               << file_name);
   if (filestr.fail()) {
      ERROR("ScalarSingletZ2DM_spectrum_plotter::write_to_file: can't open file "
            << file_name);
      return;
   }

   filestr << "### pole masses (Q = " << scale << " GeV)\n";
   write_spectrum(spectrum, filestr);

   filestr.close();
   VERBOSE_MSG("ScalarSingletZ2DM_spectrum_plotter::write_to_file: file written: "
               << file_name);
}

void ScalarSingletZ2DM_spectrum_plotter::write_spectrum(const TSpectrum& spectrum, std::ofstream& filestr) const
{
   if (!filestr.good()) {
      ERROR("ScalarSingletZ2DM_spectrum_plotter::write_spectrum: "
            "file stream is corrupted");
      return;
   }

   for (std::size_t s = 0; s < spectrum.size(); ++s) {
      const std::string& name = spectrum[s].name;
      const std::string& latex_name = spectrum[s].latex_name;
      const std::valarray<double>& masses = spectrum[s].masses;
      const std::size_t multiplicity = masses.size();

      filestr << std::left << "# " << name << '\n';

      for (std::size_t i = 0; i < multiplicity; ++i) {
         const std::string lname("$" + latex_name + "$");
         std::stringstream lname_with_index;
         lname_with_index << "$" << latex_name;
         if (multiplicity > 1)
            lname_with_index << "_{" << (i+1) << "}";
         lname_with_index  << "$";

         filestr << std::left << std::setw(width) << s
                 << std::left << std::setw(width) << masses[i]
                 << std::left << std::setw(width) << name
                 << std::left << std::setw(2*width) << lname
                 << std::left << lname_with_index.str()
                 << '\n';
      }
   }
}

namespace ScalarSingletZ2DM_database {

/**
 * write mass eigenstates to database
 *
 * @param file_name database file name
 * @param model mass eigenstates
 * @param qedqcd pointer to low-enregy data. If zero, the low-enregy
 *    data will not be written.
 * @param physical_input pointer to physical non-SLHA input parameters
 * @param observables pointer to observables struct. If zero, the
 *    observables will not be written.
 */
void to_database(
   const std::string& file_name, const ScalarSingletZ2DM_mass_eigenstates& model,
   const softsusy::QedQcd* qedqcd, const Physical_input* physical_input,
   const ScalarSingletZ2DM_observables* observables)
{
   using utilities::append;

   std::vector<std::string> names(9,"");
   Eigen::ArrayXd values(9);

   // fill settings
   names[0] = "RGE_loop_order";
   names[1] = "Threshold_loop_order";
   names[2] = "EWSB_loop_order";
   names[3] = "Pole_mass_loop_order";
   names[4] = "BVP_precision_goal";
   names[5] = "EWSB_precision_goal";
   names[6] = "Problems";
   names[7] = "Warnings";
   names[8] = "Q";
   values(0) = model.get_loops();
   values(1) = model.get_thresholds();
   values(2) = model.get_ewsb_loop_order();
   values(3) = model.get_pole_mass_loop_order();
   values(4) = model.get_precision();
   values(5) = model.get_ewsb_iteration_precision();
   values(6) = model.get_problems().have_problem();
   values(7) = model.get_problems().have_warning();
   values(8) = model.get_scale();

   // fill input parameters
   append(names, ScalarSingletZ2DM_parameter_getter().get_input_parameter_names());
   append(values, model.get_input().get());

   // fill DR-bar parameters
   append(names, ScalarSingletZ2DM_parameter_getter().get_parameter_names());
   append(values, model.get());

   // fill extra parameters
   append(names, ScalarSingletZ2DM_parameter_getter().get_extra_parameter_names());
   append(values, model.get_extra_parameters());

   // fill DR-bar masses and mixings
   append(names, ScalarSingletZ2DM_parameter_getter().get_DRbar_mass_names());
   append(names, ScalarSingletZ2DM_parameter_getter().get_DRbar_mixing_names());
   append(values, model.get_DRbar_masses_and_mixings());

   // fill pole masses and mixings
   append(names, ScalarSingletZ2DM_parameter_getter().get_pole_mass_names());
   append(names, ScalarSingletZ2DM_parameter_getter().get_pole_mixing_names());
   append(values, model.get_physical().get());

   // fill low-energy data (optional)
   if (qedqcd) {
      append(names, softsusy::QedQcd::display_input_parameter_names());
      append(values, qedqcd->display_input());
   }

   // fill extra physical input (optional)
   if (physical_input) {
      append(names, Physical_input::get_names());
      append(values, physical_input->get());
   }

   // fill observables (optional)
   if (ScalarSingletZ2DM_observables::NUMBER_OF_OBSERVABLES > 0 && observables) {
      append(names, ScalarSingletZ2DM_observables::get_names());
      append(values, observables->get());
   }

   try {
      database::Database db(file_name);
      db.insert("Point", names, values);
   } catch(const flexiblesusy::Error& e) {
      ERROR(e.what_detailed());
   }
}

namespace {
Eigen::ArrayXd extract_entry(const std::string& file_name, long long entry)
{
   database::Database db(file_name);
   Eigen::ArrayXd values;

   try {
      values = db.extract("Point", entry);
   } catch (const flexiblesusy::Error& e) {
      ERROR(e.what_detailed());
   }

   return values;
}
} // anonymous namespace

/**
 * read mass eigenstates from database
 *
 * @param file_name database file name
 * @param entry entry number (0 = first entry)
 * @param qedqcd pointer to low-energy data.  If zero, the low-energy
 *    data structure will not be filled
 * @param physical_input pointer to physical non-SLHA input.  If zero,
 *    the physical_input data structure will not be filled.
 * @param observables pointer to observables.  If zero, the observables
 *    data structure will not be filled
 *
 * @return mass eigenstates
 */
ScalarSingletZ2DM_mass_eigenstates from_database(
   const std::string& file_name, long long entry, softsusy::QedQcd* qedqcd,
   Physical_input* physical_input, ScalarSingletZ2DM_observables* observables)
{
   using utilities::append;

   ScalarSingletZ2DM_mass_eigenstates model;
   const auto number_of_parameters = model.get_number_of_parameters();
   const auto number_of_masses = ScalarSingletZ2DM_parameter_getter().get_number_of_masses();
   const auto number_of_mixings = ScalarSingletZ2DM_info::NUMBER_OF_MIXINGS;
   const auto number_of_input_parameters = ScalarSingletZ2DM_info::NUMBER_OF_INPUT_PARAMETERS;
   const auto number_of_extra_parameters = ScalarSingletZ2DM_info::NUMBER_OF_EXTRA_PARAMETERS;
   const auto number_of_low_energy_input_parameters =
      (qedqcd ? softsusy::NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS : 0);
   const auto number_of_extra_physical_input_parameters =
      (physical_input ? Physical_input::NUMBER_OF_INPUT_PARAMETERS : 0);
   const auto number_of_observables =
      (observables ? ScalarSingletZ2DM_observables::NUMBER_OF_OBSERVABLES : 0);
   const auto total_entries = 9 + number_of_input_parameters + number_of_extra_parameters
      + number_of_parameters + 2*(number_of_masses + number_of_mixings)
      + number_of_low_energy_input_parameters
      + number_of_extra_physical_input_parameters + number_of_observables;

   const Eigen::ArrayXd values(extract_entry(file_name, entry));

   if (values.rows() < total_entries) {
      ERROR("data set " << entry << " extracted from " << file_name
            << " contains " << values.rows() << " entries."
            " Expected number of entries at least: " << total_entries);
      return model;
   }

   int offset = 0;

   // restore settings
   model.set_loops(values(offset++));
   model.set_thresholds(values(offset++));
   model.set_ewsb_loop_order(values(offset++));
   model.set_pole_mass_loop_order(values(offset++));
   model.set_precision(values(offset++));
   model.set_ewsb_iteration_precision(values(offset++));
   offset++; // Problems
   offset++; // Warnings
   model.set_scale(values(offset++));

   // restore input parameters
   model.get_input().set(values.block(offset, 0, number_of_input_parameters, 1));
   offset += number_of_input_parameters;

   // restore DR-bar parameters
   model.set(values.block(offset, 0, number_of_parameters, 1));
   offset += number_of_parameters;

   // restore extra parameters
   model.set_extra_parameters(values.block(offset, 0, number_of_extra_parameters, 1));
   offset += number_of_extra_parameters;

   // restore DR-bar masses and mixings
   model.set_DRbar_masses_and_mixings(values.block(offset, 0, number_of_masses + number_of_mixings, 1));
   offset += number_of_masses + number_of_mixings;

   // restore pole masses and mixings
   model.get_physical().set(values.block(offset, 0, number_of_masses + number_of_mixings, 1));
   offset += number_of_masses + number_of_mixings;

   // fill low-energy data (optional)
   if (qedqcd) {
      qedqcd->set_input(values.block(offset, 0, number_of_low_energy_input_parameters, 1));
      offset += number_of_low_energy_input_parameters;
   }

   // fill extra physical input parameters (optional)
   if (physical_input) {
      physical_input->set(values.block(offset, 0, number_of_extra_physical_input_parameters, 1));
      offset += number_of_extra_physical_input_parameters;
   }

   // fill observables (optional)
   if (number_of_observables > 0 && observables) {
      observables->set(values.block(offset, 0, number_of_observables, 1));
      offset += number_of_observables;
   }

   if (offset != total_entries)
      throw SetupError("from_database: offset (" + ToString(offset) +
                       ") != total_entries (" + ToString(total_entries) + ").");

   return model;
}

} // namespace ScalarSingletZ2DM_database

} // namespace flexiblesusy
