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


#ifndef ScalarSingletZ2DMMhInputMsInput_SPECTRUM_GENERATOR_INTERFACE_H
#define ScalarSingletZ2DMMhInputMsInput_SPECTRUM_GENERATOR_INTERFACE_H

#include "ScalarSingletZ2DMMhInputMsInput_info.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_mass_eigenstates.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_model.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_model_slha.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_utilities.hpp"

#include "error.hpp"
#include "coupling_monitor.hpp"
#include "logger.hpp"
#include "lowe.h"
#include "spectrum_generator_problems.hpp"
#include "spectrum_generator_settings.hpp"

#include "loop_libraries/loop_library.hpp"

#include <string>
#include <tuple>

namespace flexiblesusy {

struct ScalarSingletZ2DMMhInputMsInput_input_parameters;

template <class T>
class ScalarSingletZ2DMMhInputMsInput_spectrum_generator_interface {
public:
   virtual ~ScalarSingletZ2DMMhInputMsInput_spectrum_generator_interface() = default;

   std::tuple<ScalarSingletZ2DMMhInputMsInput<T>> get_models() const
   { return std::make_tuple(model); }
   std::tuple<ScalarSingletZ2DMMhInputMsInput<T>> get_models_slha() const
   { return std::make_tuple(ScalarSingletZ2DMMhInputMsInput<T>(model, settings.get(Spectrum_generator_settings::force_positive_masses) == 0.)); }

   const ScalarSingletZ2DMMhInputMsInput<T>& get_model() const
   { return model; }
   ScalarSingletZ2DMMhInputMsInput<T>& get_model()
   { return model; }
   ScalarSingletZ2DMMhInputMsInput<T> get_model_slha() const
   { return ScalarSingletZ2DMMhInputMsInput<T>(model, settings.get(Spectrum_generator_settings::force_positive_masses) == 0.); }

   Spectrum_generator_problems get_problems() const { return problems; }
   int get_exit_code() const { return problems.have_problem(); }
   double get_reached_precision() const { return reached_precision; }
   const Spectrum_generator_settings& get_settings() const { return settings; }
   void set_parameter_output_scale(double s) { parameter_output_scale = s; }
   void set_settings(const Spectrum_generator_settings&);

   void run(const softsusy::QedQcd&, const ScalarSingletZ2DMMhInputMsInput_input_parameters&);
   void write_running_couplings(const std::string& filename, double, double) const;
   void write_spectrum(const std::string& filename = "ScalarSingletZ2DMMhInputMsInput_spectrum.dat") const;

protected:
   ScalarSingletZ2DMMhInputMsInput<T> model;
   Spectrum_generator_problems problems;
   Spectrum_generator_settings settings;
   double parameter_output_scale{0.}; ///< output scale for running parameters
   double reached_precision{std::numeric_limits<double>::infinity()}; ///< the precision that was reached

   void translate_exception_to_problem(ScalarSingletZ2DMMhInputMsInput<T>& model);
   virtual void run_except(const softsusy::QedQcd&, const ScalarSingletZ2DMMhInputMsInput_input_parameters&) = 0;
};

/**
 * Setup spectrum generator from a Spectrum_generator_settings object.
 *
 * @param settings_ spectrum generator settings
 */
template <class T>
void ScalarSingletZ2DMMhInputMsInput_spectrum_generator_interface<T>::set_settings(
   const Spectrum_generator_settings& settings_)
{
   settings = settings_;
   Loop_library::set(settings.get(Spectrum_generator_settings::loop_library));

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
void ScalarSingletZ2DMMhInputMsInput_spectrum_generator_interface<T>::run(
   const softsusy::QedQcd& qedqcd_, const ScalarSingletZ2DMMhInputMsInput_input_parameters& input)
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
void ScalarSingletZ2DMMhInputMsInput_spectrum_generator_interface<T>::write_running_couplings(
   const std::string& filename,
   double start, double stop) const
{
   ScalarSingletZ2DMMhInputMsInput_mass_eigenstates tmp_model(model);
   try {
      tmp_model.run_to(start);
   } catch (const Error& error) {
      ERROR("write_running_couplings: running to scale "
            << start << " failed: " << error.what_detailed());
      return;
   }

   // returns parameters at given scale
   auto data_getter = [&tmp_model](double scale) {
      tmp_model.run_to(scale);
      return ScalarSingletZ2DMMhInputMsInput_parameter_getter::get_parameters(tmp_model);
   };

   std::vector<std::string> parameter_names(
      std::cbegin(ScalarSingletZ2DMMhInputMsInput_info::parameter_names),
      std::cend(ScalarSingletZ2DMMhInputMsInput_info::parameter_names));

   Coupling_monitor coupling_monitor(data_getter, parameter_names);
   coupling_monitor.run(start, stop, 100, true);
   coupling_monitor.write_to_file(filename);
}

/**
 * Write spectrum (pole masses) to a text file
 *
 * @param filename output file name
 */
template <class T>
void ScalarSingletZ2DMMhInputMsInput_spectrum_generator_interface<T>::write_spectrum(
   const std::string& filename) const
{
   ScalarSingletZ2DMMhInputMsInput_spectrum_plotter plotter(model);
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
void ScalarSingletZ2DMMhInputMsInput_spectrum_generator_interface<T>::translate_exception_to_problem(ScalarSingletZ2DMMhInputMsInput<T>& model)
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
