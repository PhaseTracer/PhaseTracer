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

#include "models/RSS.hpp"
#include "models/RSS_HT.hpp"

using EffectivePotential::OneLoopPotential;
using EffectivePotential::RSS;
using EffectivePotential::RSS_HT;

// Real scalar singlet extension of the SM without Z2 symmetry (RSS)
void bind_models_rss(py::module_ &m) {

  py::class_<RSS, OneLoopPotential>(
      m, "RSS", "Real scalar singlet extension of the SM without Z2 symmetry, at one loop (RSS.hpp)")
      .def(py::init<>(), "Model without parameters; call setParameters before use")
      .def(py::init<double, double, double, double, double, double, double, double, double, double, double, double>(),
           py::arg("muh_sq"), py::arg("mus_sq"), py::arg("lh"), py::arg("ls"), py::arg("c1"), py::arg("c2"),
           py::arg("b3"), py::arg("muh_sq_tree_EWSB"), py::arg("mus_sq_tree_EWSB"), py::arg("vs"), py::arg("ms"),
           py::arg("theta"))
      .def("setParameters",
           py::overload_cast<double, double, double, double, double, double, double, double, double, double, double,
                             double>(&RSS::setParameters),
           py::arg("muh_sq"), py::arg("mus_sq"), py::arg("lh"), py::arg("ls"), py::arg("c1"), py::arg("c2"),
           py::arg("b3"), py::arg("muh_sq_tree_EWSB"), py::arg("mus_sq_tree_EWSB"), py::arg("vs"), py::arg("ms"),
           py::arg("theta"))
      .def("printParameters", &RSS::printParameters, py::call_guard<py::scoped_ostream_redirect>())
      .def("get_EW_VEV", &RSS::get_EW_VEV, "Electroweak vacuum (h, s) at zero temperature")
      .def("set_useGSResummation", &RSS::set_useGSResummation, py::arg("flag"), "Goldstone resummation on or off")
      .def_property("bUseBoltzmannSuppression", &RSS::get_bUseBoltzmannSuppression,
                    &RSS::set_bUseBoltzmannSuppression)
      .def_property("MAX_BOLTZ_EXP", &RSS::get_MAX_BOLTZ_EXP, &RSS::set_MAX_BOLTZ_EXP)
      .def_property("DEBUG", &RSS::get_DEBUG, &RSS::set_DEBUG);

  py::class_<RSS_HT, RSS>(m, "RSS_HT", "High-temperature expansion of the RSS model (RSS_HT.hpp)")
      .def(py::init<>(), "Model without parameters; call setParameters before use")
      .def(py::init<double, double, double, double, double, double, double, double, double, double>(),
           py::arg("muh_sq"), py::arg("mus_sq"), py::arg("lh"), py::arg("ls"), py::arg("c1"), py::arg("c2"),
           py::arg("b3"), py::arg("vs"), py::arg("ms"), py::arg("theta"));
}
