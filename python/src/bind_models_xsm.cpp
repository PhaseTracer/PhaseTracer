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

#include <vector>

#include "models/xSM_HT.hpp"
#include "models/xSM_MSbar.hpp"
#include "models/xSM_OSlike.hpp"

using EffectivePotential::OneLoopPotential;
using EffectivePotential::Potential;
using EffectivePotential::xSM_base;
using EffectivePotential::xSM_HT;
using EffectivePotential::xSM_MSbar;
using EffectivePotential::xSM_OSlike;

// Real scalar singlet extensions of the Standard Model (xSM)
void bind_models_xsm(py::module_ &m) {

  py::class_<xSM_base, OneLoopPotential>(
      m, "xSM_base", "Common base of xSM_MSbar and xSM_OSlike (not constructible): SM inputs and the field content")
      .def("set_SM_parameters", &xSM_base::set_SM_parameters, py::arg("SM_parameters"),
           "Seven SM inputs [mh, v, gp, g, yt_sq, yb_sq, ytau_sq]");

  py::class_<xSM_MSbar, xSM_base>(
      m, "xSM_MSbar", "Z2-symmetric xSM in the MS-bar scheme (xSM_MSbar.hpp). Use from_tadpoles to build a consistent model.")
      .def(py::init<double, double, double>(), py::arg("lambda_hs"), py::arg("lambda_s"), py::arg("ms"))
      .def_static("from_tadpoles", &xSM_MSbar::from_tadpoles, py::arg("lambda_hs"), py::arg("lambda_s"),
                  py::arg("ms"), py::arg("Q"), py::arg("xi"), py::arg("use_covariant_gauge") = false,
                  py::arg("use_1L_EWSB_in_0L_mass") = false, py::arg("use_Goldstone_resum") = true,
                  py::arg("use_tree_level_tadpole") = false, py::arg("SM_parameters") = std::vector<double>{},
                  "Model whose Lagrangian parameters satisfy the tree-level or one-loop tadpole conditions\n"
                  "at renormalisation scale Q, in gauge xi")
      .def("apply_tree_level", &xSM_MSbar::apply_tree_level, "Fix muh_sq, mus_sq and lambda_h from tree-level EWSB")
      .def("apply_one_loop", &xSM_MSbar::apply_one_loop, py::arg("tol") = 0.1,
           "Iterate the one-loop tadpole conditions; returns whether they converged")
      .def("iterate_one_loop", &xSM_MSbar::iterate_one_loop)
      .def("get_muh_sq", &xSM_MSbar::get_muh_sq)
      .def("get_mus_sq", &xSM_MSbar::get_mus_sq)
      .def("get_lambda_h", &xSM_MSbar::get_lambda_h)
      .def("get_lambda_s", &xSM_MSbar::get_lambda_s)
      .def("get_lambda_hs", &xSM_MSbar::get_lambda_hs)
      .def("get_v_tree_s", &xSM_MSbar::get_v_tree_s, "Tree-level singlet vev")
      .def_property("use_analytic_gradient", &xSM_MSbar::get_use_analytic_gradient,
                    &xSM_MSbar::set_use_analytic_gradient)
      .def_readonly("iteration_converged", &xSM_MSbar::iteration_converged,
                    "Whether from_tadpoles' one-loop iteration converged");

  py::class_<xSM_OSlike, xSM_base>(
      m, "xSM_OSlike", "Z2-symmetric xSM with on-shell-like renormalisation (xSM_OSlike.hpp)")
      .def(py::init<double, double, double>(), py::arg("lambda_hs"), py::arg("lambda_s"), py::arg("ms"))
      // the class hides the base getter with a non-const one that returns its own scale
      .def("get_renormalization_scale", [](xSM_OSlike &model) { return model.get_renormalization_scale(); });

  py::class_<xSM_HT, Potential>(
      m, "xSM_HT", "High-temperature expansion of the Z2-symmetric xSM, with analytic results (xSM_HT.hpp)")
      .def(py::init<double, double, double>(), py::arg("lambda_hs"), py::arg("lambda_s"), py::arg("m_s"))
      .def("get_TC_from_expression", &xSM_HT::get_TC_from_expression, "Analytic critical temperature")
      .def("get_vs_from_expression", &xSM_HT::get_vs_from_expression,
           "Analytic singlet vev at the critical temperature")
      .def("get_vh_from_expression", &xSM_HT::get_vh_from_expression,
           "Analytic Higgs vev at the critical temperature");
}
