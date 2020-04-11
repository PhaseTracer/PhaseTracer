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

// File generated at Sat 11 Apr 2020 12:54:04

#ifndef THDMIISNMSSMBCsimple_SPECTRUM_GENERATOR_INTERFACE_H
#define THDMIISNMSSMBCsimple_SPECTRUM_GENERATOR_INTERFACE_H

#include "THDMIISNMSSMBCsimple_mass_eigenstates.hpp"
#include "THDMIISNMSSMBCsimple_model.hpp"
#include "THDMIISNMSSMBCsimple_model_slha.hpp"
#include "THDMIISNMSSMBCsimple_utilities.hpp"

#include "error.hpp"
#include "coupling_monitor.hpp"
#include "logger.hpp"
#include "lowe.h"
#include "spectrum_generator_problems.hpp"
#include "spectrum_generator_settings.hpp"
#include "loop_corrections.hpp"

#include <string>
#include <tuple>

namespace flexiblesusy {

struct THDMIISNMSSMBCsimple_input_parameters;

template <class T>
class THDMIISNMSSMBCsimple_spectrum_generator_interface {
public:
   virtual ~THDMIISNMSSMBCsimple_spectrum_generator_interface() = default;

   std::tuple<THDMIISNMSSMBCsimple<T>> get_models() const
   { return std::make_tuple(model); }
   std::tuple<THDMIISNMSSMBCsimple_slha<THDMIISNMSSMBCsimple<T>>> get_models_slha() const
   { return std::make_tuple(THDMIISNMSSMBCsimple_slha<THDMIISNMSSMBCsimple<T> >(model, settings.get(Spectrum_generator_settings::force_positive_masses) == 0.)); }

   const THDMIISNMSSMBCsimple<T>& get_model() const
   { return model; }
   THDMIISNMSSMBCsimple<T>& get_model()
   { return model; }
   THDMIISNMSSMBCsimple_slha<THDMIISNMSSMBCsimple<T>> get_model_slha() const
   { return THDMIISNMSSMBCsimple_slha<THDMIISNMSSMBCsimple<T>>(model, settings.get(Spectrum_generator_settings::force_positive_masses) == 0.); }

   Spectrum_generator_problems get_problems() const { return problems; }
   int get_exit_code() const { return problems.have_problem(); }
   double get_reached_precision() const { return reached_precision; }
   const Spectrum_generator_settings& get_settings() const { return settings; }
   void set_parameter_output_scale(double s) { parameter_output_scale = s; }
   void set_settings(const Spectrum_generator_settings&);

   void run(const softsusy::QedQcd&, const THDMIISNMSSMBCsimple_input_parameters&);
   void write_running_couplings(const std::string& filename, double, double) const;
   void write_spectrum(const std::string& filename = "THDMIISNMSSMBCsimple_spectrum.dat") const;

protected:
   THDMIISNMSSMBCsimple<T> model;
   Spectrum_generator_problems problems;
   Spectrum_generator_settings settings;
   double parameter_output_scale{0.}; ///< output scale for running parameters
   double reached_precision{std::numeric_limits<double>::infinity()}; ///< the precision that was reached

   void translate_exception_to_problem(THDMIISNMSSMBCsimple<T>& model);
   virtual void run_except(const softsusy::QedQcd&, const THDMIISNMSSMBCsimple_input_parameters&) = 0;
};

/**
 * Setup spectrum generator from a Spectrum_generator_settings object.
 *
 * @param settings_ spectrum generator settings
 */
template <class T>
void THDMIISNMSSMBCsimple_spectrum_generator_interface<T>::set_settings(
   const Spectrum_generator_settings& settings_)
{
   settings = settings_;
   model.set_pole_mass_loop_order(settings.get(Spectrum_generator_settings::pole_mass_loop_order));
   model.set_ewsb_loop_order(settings.get(Spectrum_generator_settings::ewsb_loop_order));
   model.set_loop_corrections(settings.get_loop_corrections());
   model.set_threshold_corrections(settings.get_threshold_corrections());
}

/**
 * @brief Run's the RG solver with the given input parameters
 *
 * This function calls run_except() from the derived class and
 * translates an emitted exception into an problem code.
 *
 * @param qedqcd_ Standard Model input parameters
 * @param input model input parameters
 */
template <class T>
void THDMIISNMSSMBCsimple_spectrum_generator_interface<T>::run(
   const softsusy::QedQcd& qedqcd_, const THDMIISNMSSMBCsimple_input_parameters& input)
{
   softsusy::QedQcd qedqcd = qedqcd_;

   try {
      qedqcd.to(qedqcd.displayPoleMZ());
      this->run_except(qedqcd, input);
   } catch (...) {
      this->translate_exception_to_problem(model);
   }

   problems.set_model_problems({ model.get_problems() });
}

/**
 * Create a text file which contains the values of all model
 * parameters at all scales between the low-scale and the high-scale.
 *
 * @param filename name of output file
 * @param start lowest scale
 * @param stop highest scale
 */
template <class T>
void THDMIISNMSSMBCsimple_spectrum_generator_interface<T>::write_running_couplings(
   const std::string& filename,
   double start, double stop) const
{
   THDMIISNMSSMBCsimple_mass_eigenstates tmp_model(model);
   try {
      tmp_model.run_to(start);
   } catch (const Error& error) {
      ERROR("write_running_couplings: running to scale "
            << start << " failed: " << error.what_detailed());
      return;
   }

   THDMIISNMSSMBCsimple_parameter_getter parameter_getter;
   Coupling_monitor<THDMIISNMSSMBCsimple_mass_eigenstates, THDMIISNMSSMBCsimple_parameter_getter>
      coupling_monitor(tmp_model, parameter_getter);

   coupling_monitor.run(start, stop, 100, true);
   coupling_monitor.write_to_file(filename);
}

/**
 * Write spectrum (pole masses) to a text file
 *
 * @param filename output file name
 */
template <class T>
void THDMIISNMSSMBCsimple_spectrum_generator_interface<T>::write_spectrum(
   const std::string& filename) const
{
   THDMIISNMSSMBCsimple_spectrum_plotter plotter(model);
   plotter.write_to_file(filename);
}

/**
 * Flags problems in the given model class from the current pending
 * exception.
 *
 * This function assumes that there is an active exception.
 *
 * @param model model class
 */
template <class T>
void THDMIISNMSSMBCsimple_spectrum_generator_interface<T>::translate_exception_to_problem(THDMIISNMSSMBCsimple<T>& model)
{
   try {
      throw;
   } catch (const NoConvergenceError&) {
      problems.flag_no_convergence();
   } catch (const NonPerturbativeRunningError& error) {
      model.get_problems().flag_no_perturbative();
      model.get_problems().flag_non_perturbative_parameter(
         error.get_parameter_index(), error.get_parameter_value(),
         error.get_scale());
   } catch (const NonPerturbativeRunningQedQcdError& error) {
      model.get_problems().flag_no_perturbative();
      model.get_problems().flag_thrown(error.what_detailed());
   } catch (const NoSinThetaWConvergenceError&) {
      model.get_problems().flag_no_sinThetaW_convergence();
   } catch (const Error& error) {
      model.get_problems().flag_thrown(error.what_detailed());
   } catch (const std::exception& error) {
      model.get_problems().flag_thrown(error.what());
   }
}

} // namespace flexiblesusy

#endif
