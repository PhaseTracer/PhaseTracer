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

// File generated at Sat 24 Oct 2020 17:07:54

#ifndef ScalarSingletZ2DM_SPECTRUM_GENERATOR_INTERFACE_H
#define ScalarSingletZ2DM_SPECTRUM_GENERATOR_INTERFACE_H

#include "ScalarSingletZ2DM_mass_eigenstates.hpp"
#include "ScalarSingletZ2DM_model.hpp"
#include "ScalarSingletZ2DM_model_slha.hpp"
#include "ScalarSingletZ2DM_utilities.hpp"

#include "error.hpp"
#include "coupling_monitor.hpp"
#include "logger.hpp"
#include "lowe.h"
#include "spectrum_generator_problems.hpp"
#include "spectrum_generator_settings.hpp"
#include "standard_model.hpp"
#include "loop_corrections.hpp"

#include <string>
#include <tuple>

namespace flexiblesusy {

struct ScalarSingletZ2DM_input_parameters;

template <class T>
class ScalarSingletZ2DM_spectrum_generator_interface {
public:
   virtual ~ScalarSingletZ2DM_spectrum_generator_interface() = default;

   std::tuple<ScalarSingletZ2DM<T>, standard_model::StandardModel<T>> get_models() const
   { return std::make_tuple(model, eft); }
   std::tuple<ScalarSingletZ2DM_slha<ScalarSingletZ2DM<T>>, standard_model::StandardModel<T>> get_models_slha() const
   { return std::make_tuple(ScalarSingletZ2DM_slha<ScalarSingletZ2DM<T>>(model, settings.get(Spectrum_generator_settings::force_positive_masses) == 0.), eft); }

   const ScalarSingletZ2DM<T>& get_model() const
   { return model; }
   ScalarSingletZ2DM<T>& get_model()
   { return model; }
   ScalarSingletZ2DM_slha<ScalarSingletZ2DM<T>> get_model_slha() const
   { return ScalarSingletZ2DM_slha<ScalarSingletZ2DM<T>>(model, settings.get(Spectrum_generator_settings::force_positive_masses) == 0.); }

   const standard_model::StandardModel<T>& get_sm() const
   { return eft; }
   standard_model::StandardModel<T>& get_sm()
   { return eft; }

   Spectrum_generator_problems get_problems() const { return problems; }
   int get_exit_code() const { return problems.have_problem(); }
   double get_reached_precision() const { return reached_precision; }
   const Spectrum_generator_settings& get_settings() const { return settings; }
   void set_parameter_output_scale(double s) { parameter_output_scale = s; }
   void set_settings(const Spectrum_generator_settings&);

   void run(const softsusy::QedQcd&, const ScalarSingletZ2DM_input_parameters&);
   void write_running_couplings(const std::string& filename, double, double) const;
   void write_spectrum(const std::string& filename = "ScalarSingletZ2DM_spectrum.dat") const;

protected:
   ScalarSingletZ2DM<T> model;
   standard_model::StandardModel<T> eft{};
   Spectrum_generator_problems problems;
   Spectrum_generator_settings settings;
   double parameter_output_scale{0.}; ///< output scale for running parameters
   double reached_precision{std::numeric_limits<double>::infinity()}; ///< the precision that was reached

   void translate_exception_to_problem(ScalarSingletZ2DM<T>& model);
   virtual void run_except(const softsusy::QedQcd&, const ScalarSingletZ2DM_input_parameters&) = 0;
};

/**
 * Setup spectrum generator from a Spectrum_generator_settings object.
 *
 * @param settings_ spectrum generator settings
 */
template <class T>
void ScalarSingletZ2DM_spectrum_generator_interface<T>::set_settings(
   const Spectrum_generator_settings& settings_)
{
   settings = settings_;
   model.set_pole_mass_loop_order(settings.get(Spectrum_generator_settings::pole_mass_loop_order));
   model.set_ewsb_loop_order(settings.get(Spectrum_generator_settings::ewsb_loop_order));
   model.set_loop_corrections(settings.get_loop_corrections());
   model.set_threshold_corrections(settings.get_threshold_corrections());
   eft.set_pole_mass_loop_order(settings.get(Spectrum_generator_settings::pole_mass_loop_order));
   eft.set_ewsb_loop_order(settings.get(Spectrum_generator_settings::ewsb_loop_order));
   eft.set_loop_corrections(settings.get_loop_corrections());
   eft.set_threshold_corrections(settings.get_threshold_corrections());
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
void ScalarSingletZ2DM_spectrum_generator_interface<T>::run(
   const softsusy::QedQcd& qedqcd_, const ScalarSingletZ2DM_input_parameters& input)
{
   softsusy::QedQcd qedqcd = qedqcd_;

   try {
      qedqcd.to(qedqcd.displayPoleMZ());
      this->run_except(qedqcd, input);
   } catch (...) {
      this->translate_exception_to_problem(model);
   }

   problems.set_model_problems({ model.get_problems(), eft.get_problems() });
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
void ScalarSingletZ2DM_spectrum_generator_interface<T>::write_running_couplings(
   const std::string& filename,
   double start, double stop) const
{
   ScalarSingletZ2DM_mass_eigenstates tmp_model(model);
   try {
      tmp_model.run_to(start);
   } catch (const Error& error) {
      ERROR("write_running_couplings: running to scale "
            << start << " failed: " << error.what_detailed());
      return;
   }

   ScalarSingletZ2DM_parameter_getter parameter_getter;
   Coupling_monitor<ScalarSingletZ2DM_mass_eigenstates, ScalarSingletZ2DM_parameter_getter>
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
void ScalarSingletZ2DM_spectrum_generator_interface<T>::write_spectrum(
   const std::string& filename) const
{
   ScalarSingletZ2DM_spectrum_plotter plotter(model);
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
void ScalarSingletZ2DM_spectrum_generator_interface<T>::translate_exception_to_problem(ScalarSingletZ2DM<T>& model)
{
   try {
      throw;
   } catch (const NoConvergenceError&) {
      problems.flag_no_convergence();
   } catch (const NonPerturbativeRunningError& error) {
      model.get_problems().flag_no_perturbative();
      eft.get_problems().flag_no_perturbative();
      model.get_problems().flag_non_perturbative_parameter(
         error.get_parameter_index(), error.get_parameter_value(),
         error.get_scale());
   } catch (const NonPerturbativeRunningQedQcdError& error) {
      model.get_problems().flag_no_perturbative();
      model.get_problems().flag_thrown(error.what_detailed());
      eft.get_problems().flag_no_perturbative();
      eft.get_problems().flag_thrown(error.what_detailed());
   } catch (const NoSinThetaWConvergenceError&) {
      model.get_problems().flag_no_sinThetaW_convergence();
      eft.get_problems().flag_no_sinThetaW_convergence();
   } catch (const Error& error) {
      model.get_problems().flag_thrown(error.what_detailed());
      eft.get_problems().flag_thrown(error.what_detailed());
   } catch (const std::exception& error) {
      model.get_problems().flag_thrown(error.what());
      eft.get_problems().flag_thrown(error.what());
   }
}

} // namespace flexiblesusy

#endif
