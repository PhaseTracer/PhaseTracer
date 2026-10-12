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

#include <string>

#include "models/conformal_dark_u1.hpp"
#include "models/scalar_yukawa_3d.hpp"
#include "models/supercool_model.hpp"

using EffectivePotential::ConformalDarkU1;
using EffectivePotential::OneLoopPotential;
using EffectivePotential::Potential;
using EffectivePotential::ScalarYukawa3D;
using EffectivePotential::SuperCoolModel;

// Other beyond-the-SM models
void bind_models_bsm(py::module_ &m) {

  py::class_<ConformalDarkU1, OneLoopPotential>(
      m, "ConformalDarkU1",
      "Conformal dark U(1) model of arXiv:2502.19478 (conformal_dark_u1.hpp).\n"
      "The constructor raises RuntimeError('bad_Higgs_quartic') if (g, y, v) give an inconsistent Higgs quartic.")
      .def(py::init<double, double, double>(), py::arg("g"), py::arg("y"), py::arg("v"));

  py::class_<SuperCoolModel, OneLoopPotential>(
      m, "SuperCoolModel", "Model with strong supercooling (supercool_model.hpp)")
      .def(py::init<double, double, double, double>(), py::arg("kappa"), py::arg("mu_sq"), py::arg("lambda_"),
           py::arg("mu0_sq"))
      .def(py::init<std::string>(), py::arg("input_file"), "Read kappa, mu_sq, lambda and mu0_sq from a file")
      .def_property("mh", &SuperCoolModel::get_mh, &SuperCoolModel::set_mh)
      .def_property("mu_sq", &SuperCoolModel::get_mu_sq, &SuperCoolModel::set_mu_sq)
      .def_property("kappa", &SuperCoolModel::get_kappa, &SuperCoolModel::set_kappa)
      .def_property("lambda_", &SuperCoolModel::get_lambda, &SuperCoolModel::set_lambda,
                    "The quartic coupling (`lambda` is a Python keyword)")
      .def_property("mu0_sq", &SuperCoolModel::get_mu0_sq, &SuperCoolModel::set_mu0_sq)
      .def_property("vSq", &SuperCoolModel::get_vSq, &SuperCoolModel::set_vSq)
      .def_property("yt", &SuperCoolModel::get_yt, &SuperCoolModel::set_yt)
      .def_property("g", &SuperCoolModel::get_g, &SuperCoolModel::set_g)
      .def_property("g1", &SuperCoolModel::get_g1, &SuperCoolModel::set_g1)
      .def_property("raddof", &SuperCoolModel::get_raddof, &SuperCoolModel::set_raddof)
      .def_property("bUseBoltzmannSuppression", &SuperCoolModel::get_bUseBoltzmannSuppression,
                    &SuperCoolModel::set_bUseBoltzmannSuppression)
      .def_property("MAX_BOLTZ_EXP", &SuperCoolModel::get_MAX_BOLTZ_EXP, &SuperCoolModel::set_MAX_BOLTZ_EXP);

  py::class_<ScalarYukawa3D, Potential>(
      m, "ScalarYukawa3D", "Scalar-Yukawa toy model in three dimensions (scalar_yukawa_3d.hpp)")
      .def(py::init<>())
      .def_property("s", &ScalarYukawa3D::get_s, &ScalarYukawa3D::set_s)
      .def_property("msq", &ScalarYukawa3D::get_msq, &ScalarYukawa3D::set_msq)
      .def_property("g", &ScalarYukawa3D::get_g, &ScalarYukawa3D::set_g)
      .def_property("lam", &ScalarYukawa3D::get_lam, &ScalarYukawa3D::set_lam)
      .def_property("mf", &ScalarYukawa3D::get_mf, &ScalarYukawa3D::set_mf)
      .def_property("y", &ScalarYukawa3D::get_y, &ScalarYukawa3D::set_y);
}
