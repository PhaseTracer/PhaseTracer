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

#include <memory>

// THDM.hpp has `using namespace std;` at global scope, so it is kept in a file of its own
#include "models/THDM.hpp"

using EffectivePotential::OneLoopPotential;

// Two-Higgs-doublet model of the TransitionListener2 benchmark
void bind_models_thdm(py::module_ &m) {

     using EffectivePotential::THDM;
     py::class_<THDM, OneLoopPotential>(
          m, "THDM",
          "CP-conserving 2HDM of the TransitionListener2 benchmark (models/THDM.hpp),\n"
          "with BSMPT-style counterterms. Fields (phi1, phi2) in GeV.")
          .def(py::init([](double tan_beta, double m12_sq, double lambda1, double lambda2, double lambda3,
                           double lambda4, double lambda5, int yukawa_type, double RG_scale) {
                    auto model = std::make_unique<THDM>();
                    model->init_params(tan_beta, m12_sq, lambda1, lambda2, lambda3, lambda4, lambda5, yukawa_type,
                                       RG_scale);
                    return model;
               }),
               py::arg("tan_beta"), py::arg("m12_sq"), py::arg("lambda1"), py::arg("lambda2"), py::arg("lambda3"),
               py::arg("lambda4"), py::arg("lambda5"), py::arg("yukawa_type") = 2, py::arg("RG_scale") = -1.,
               "m12_sq in GeV^2; yukawa_type 1-4 (I, II, lepton specific, flipped); RG_scale <= 0 means v")
          .def("init_params", &THDM::init_params, py::arg("tan_beta"), py::arg("m12_sq"), py::arg("lambda1"),
               py::arg("lambda2"), py::arg("lambda3"), py::arg("lambda4"), py::arg("lambda5"),
               py::arg("yukawa_type") = 2, py::arg("RG_scale") = -1.)
          .def("get_counterterms", &THDM::get_counterterms, "Counterterms in TransitionListener's normalisation")
          .def("get_diagnostics", &THDM::get_diagnostics, "Derived couplings, for cross-checks with TransitionListener");
}
