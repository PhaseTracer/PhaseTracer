// ====================================================================
// This file is part of PhaseTracer

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ====================================================================

#include "bindings.hpp"

#include "models/1D_test_model.hpp"
#include "models/2D_test_model.hpp"
#include "models/Z2_scalar_singlet_model.hpp"

using EffectivePotential::OneLoopPotential;
using EffectivePotential::Potential;

void bind_models(py::module_ &m) {

  py::class_<EffectivePotential::OneDimModel, Potential>(
      m, "OneDimModel",
      "One-field polynomial test model with an analytic critical temperature:\n"
      "V = (c T^2 - m2) phi^2 + kappa phi^3 + lambda phi^4 (see 1D_test_model.hpp)")
      .def(py::init<>())
      .def_property("m2", &EffectivePotential::OneDimModel::get_m2, &EffectivePotential::OneDimModel::set_m2)
      .def_property("kappa", &EffectivePotential::OneDimModel::get_kappa, &EffectivePotential::OneDimModel::set_kappa)
      .def_property("lambda_", &EffectivePotential::OneDimModel::get_lambda,
                    &EffectivePotential::OneDimModel::set_lambda, "The quartic coupling (`lambda` is a Python keyword)")
      .def_property("c", &EffectivePotential::OneDimModel::get_c, &EffectivePotential::OneDimModel::set_c)
      .def("get_TC_from_expression", &EffectivePotential::OneDimModel::get_TC_from_expression,
           "Analytic critical temperature")
      .def("get_true_vacuum_from_expression", &EffectivePotential::OneDimModel::get_true_vacuum_from_expression,
           "Analytic true vacuum at the critical temperature")
      .def("get_false_vacuum_from_expression", &EffectivePotential::OneDimModel::get_false_vacuum_from_expression,
           "Analytic false vacuum at the critical temperature");

  py::class_<EffectivePotential::TwoDimModel, OneLoopPotential>(
      m, "TwoDimModel", "Two-field one-loop test model (model1 of CosmoTransitions)")
      .def(py::init<>());

  py::class_<EffectivePotential::Z2ScalarSingletModel, Potential>(
      m, "Z2ScalarSingletModel",
      "High-temperature expansion of the Z2-symmetric real scalar singlet extension of the SM,\n"
      "with analytic results from arXiv:1611.02073")
      .def(py::init<>())
      .def("set_m_s", &EffectivePotential::Z2ScalarSingletModel::set_m_s, py::arg("m_s"), "Singlet mass (default 27)")
      .def("set_lambda_hs", &EffectivePotential::Z2ScalarSingletModel::set_lambda_hs, py::arg("lambda_hs"),
           "Higgs-singlet portal coupling (default 0.25)")
      .def("get_TC_from_expression", &EffectivePotential::Z2ScalarSingletModel::get_TC_from_expression,
           "Analytic critical temperature")
      .def("get_vs_from_expression", &EffectivePotential::Z2ScalarSingletModel::get_vs_from_expression,
           "Analytic singlet vev at the critical temperature")
      .def("get_vh_from_expression", &EffectivePotential::Z2ScalarSingletModel::get_vh_from_expression,
           "Analytic Higgs vev at the critical temperature");
}
