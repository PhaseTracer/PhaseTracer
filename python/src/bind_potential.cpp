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

#include "potential.hpp"
#include "one_loop_potential.hpp"

using EffectivePotential::DaisyMethod;
using EffectivePotential::OneLoopPotential;
using EffectivePotential::Potential;

namespace {

/**
 * Trampoline for Potential. V and get_n_scalars must be supplied by the Python subclass;
 * everything else falls back to the C++ default (numerical derivatives, no symmetries,
 * nothing forbidden) unless overridden.
 */
class PyPotential : public Potential {
public:
  using Potential::Potential;

  double V(Eigen::VectorXd phi, double T) const override {
    PYBIND11_OVERRIDE_PURE(double, Potential, V, phi, T);
  }
  size_t get_n_scalars() const override {
    PYBIND11_OVERRIDE_PURE(size_t, Potential, get_n_scalars, );
  }
  bool forbidden(Eigen::VectorXd phi) const override {
    PYBIND11_OVERRIDE(bool, Potential, forbidden, phi);
  }
  std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
    PYBIND11_OVERRIDE(std::vector<Eigen::VectorXd>, Potential, apply_symmetry, phi);
  }
  std::vector<std::vector<int>> get_symmetry_axes() const override {
    PYBIND11_OVERRIDE(std::vector<std::vector<int>>, Potential, get_symmetry_axes, );
  }
  std::vector<Eigen::VectorXd> get_low_t_phases() const override {
    PYBIND11_OVERRIDE(std::vector<Eigen::VectorXd>, Potential, get_low_t_phases, );
  }
  Eigen::VectorXd dV_dx(Eigen::VectorXd phi, double T) const override {
    PYBIND11_OVERRIDE(Eigen::VectorXd, Potential, dV_dx, phi, T);
  }
  Eigen::VectorXd d2V_dxdt(Eigen::VectorXd phi, double T) const override {
    PYBIND11_OVERRIDE(Eigen::VectorXd, Potential, d2V_dxdt, phi, T);
  }
  Eigen::MatrixXd d2V_dx2(Eigen::VectorXd phi, double T) const override {
    PYBIND11_OVERRIDE(Eigen::MatrixXd, Potential, d2V_dx2, phi, T);
  }
};

/** Trampoline for OneLoopPotential: V0 and get_n_scalars must be supplied, the rest have defaults. */
class PyOneLoopPotential : public OneLoopPotential {
public:
  using OneLoopPotential::OneLoopPotential;

  double V0(Eigen::VectorXd phi) const override {
    PYBIND11_OVERRIDE_PURE(double, OneLoopPotential, V0, phi);
  }
  size_t get_n_scalars() const override {
    PYBIND11_OVERRIDE_PURE(size_t, OneLoopPotential, get_n_scalars, );
  }
  bool forbidden(Eigen::VectorXd phi) const override {
    PYBIND11_OVERRIDE(bool, OneLoopPotential, forbidden, phi);
  }
  std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
    PYBIND11_OVERRIDE(std::vector<Eigen::VectorXd>, OneLoopPotential, apply_symmetry, phi);
  }
  std::vector<std::vector<int>> get_symmetry_axes() const override {
    PYBIND11_OVERRIDE(std::vector<std::vector<int>>, OneLoopPotential, get_symmetry_axes, );
  }
  std::vector<Eigen::VectorXd> get_low_t_phases() const override {
    PYBIND11_OVERRIDE(std::vector<Eigen::VectorXd>, OneLoopPotential, get_low_t_phases, );
  }
  std::vector<double> get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const override {
    PYBIND11_OVERRIDE(std::vector<double>, OneLoopPotential, get_scalar_masses_sq, phi, xi);
  }
  std::vector<double> get_fermion_masses_sq(Eigen::VectorXd phi) const override {
    PYBIND11_OVERRIDE(std::vector<double>, OneLoopPotential, get_fermion_masses_sq, phi);
  }
  std::vector<double> get_vector_masses_sq(Eigen::VectorXd phi) const override {
    PYBIND11_OVERRIDE(std::vector<double>, OneLoopPotential, get_vector_masses_sq, phi);
  }
  std::vector<double> get_ghost_masses_sq(Eigen::VectorXd phi, double xi) const override {
    PYBIND11_OVERRIDE(std::vector<double>, OneLoopPotential, get_ghost_masses_sq, phi, xi);
  }
  std::vector<double> get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const override {
    PYBIND11_OVERRIDE(std::vector<double>, OneLoopPotential, get_scalar_debye_sq, phi, xi, T);
  }
  std::vector<double> get_scalar_thermal_sq(double T) const override {
    PYBIND11_OVERRIDE(std::vector<double>, OneLoopPotential, get_scalar_thermal_sq, T);
  }
  std::vector<double> get_vector_debye_sq(Eigen::VectorXd phi, double T) const override {
    PYBIND11_OVERRIDE(std::vector<double>, OneLoopPotential, get_vector_debye_sq, phi, T);
  }
  std::vector<double> get_scalar_dofs() const override {
    PYBIND11_OVERRIDE(std::vector<double>, OneLoopPotential, get_scalar_dofs, );
  }
  std::vector<double> get_fermion_dofs() const override {
    PYBIND11_OVERRIDE(std::vector<double>, OneLoopPotential, get_fermion_dofs, );
  }
  std::vector<double> get_vector_dofs() const override {
    PYBIND11_OVERRIDE(std::vector<double>, OneLoopPotential, get_vector_dofs, );
  }
  std::vector<double> get_ghost_dofs() const override {
    PYBIND11_OVERRIDE(std::vector<double>, OneLoopPotential, get_ghost_dofs, );
  }
  double counter_term(Eigen::VectorXd phi, double T) const override {
    PYBIND11_OVERRIDE(double, OneLoopPotential, counter_term, phi, T);
  }
};

} // namespace

void bind_potential(py::module_ &m) {

  py::enum_<DaisyMethod>(m, "DaisyMethod", "Treatment of the thermal (daisy) masses")
      .value("NoDaisy", DaisyMethod::None)
      .value("ArnoldEspinosa", DaisyMethod::ArnoldEspinosa)
      .value("Parwani", DaisyMethod::Parwani);

  py::class_<Potential, PyPotential>(
      m, "Potential",
      "Finite-temperature effective potential.\n\n"
      "Subclass this in Python to define a model: call super().__init__() and provide\n"
      "V(phi, T) and get_n_scalars(). The gradient and Hessian default to numerical\n"
      "derivatives unless dV_dx / d2V_dx2 are overridden.\n\n"
      "A Python V holds the GIL on every call, so PhaseTracer's OpenMP loops run\n"
      "serially; call phasetracer.set_num_threads(1) to avoid idle threads.")
      .def(py::init<>())
      .def("V", &Potential::V, py::arg("phi"), py::arg("T"), "Potential at field values phi and temperature T")
      .def("__call__", &Potential::operator(), py::arg("phi"), py::arg("T"))
      .def("get_n_scalars", &Potential::get_n_scalars, "Number of scalar fields")
      .def("forbidden", &Potential::forbidden, py::arg("phi"), "Whether these field values are excluded")
      .def("apply_symmetry", &Potential::apply_symmetry, py::arg("phi"), "Symmetry partners of a field point")
      .def("get_symmetry_axes", &Potential::get_symmetry_axes, "Z2 symmetry axes of the potential")
      .def("get_low_t_phases", &Potential::get_low_t_phases, "Expected vacua at low temperature")
      .def("dV_dx", &Potential::dV_dx, py::arg("phi"), py::arg("T"), "Gradient of the potential")
      .def("d2V_dxdt", &Potential::d2V_dxdt, py::arg("phi"), py::arg("T"), "Temperature derivative of the gradient")
      .def("d2V_dx2", &Potential::d2V_dx2, py::arg("phi"), py::arg("T"), "Hessian of the potential")
      .def("set_h_4", &Potential::set_h_4, py::arg("h_4"),
           "Use fourth-order (True) or second-order (False) numerical derivatives")
      .def("get_h_4", &Potential::get_h_4)
      .def_property("h", &Potential::get_h, &Potential::set_h, "Step size of numerical derivatives")
      .def_property("field_scale", &Potential::get_field_scale, &Potential::set_field_scale)
      .def_property("temperature_scale", &Potential::get_temperature_scale, &Potential::set_temperature_scale);

  py::class_<OneLoopPotential, Potential, PyOneLoopPotential>(
      m, "OneLoopPotential",
      "One-loop thermal effective potential built from field-dependent masses.\n\n"
      "Subclass this in Python: call super().__init__() and provide V0(phi), get_n_scalars()\n"
      "and the mass and dof methods the model needs. V(phi, T) is assembled from V0, the\n"
      "Coleman-Weinberg and thermal corrections and the daisy resummation.")
      .def(py::init<>())
      .def("V0", &OneLoopPotential::V0, py::arg("phi"), "Tree-level potential")
      .def("V1", py::overload_cast<Eigen::VectorXd, double>(&OneLoopPotential::V1, py::const_),
           py::arg("phi"), py::arg("T") = 0., "Zero-temperature one-loop correction")
      .def("V1T", py::overload_cast<Eigen::VectorXd, double>(&OneLoopPotential::V1T, py::const_),
           py::arg("phi"), py::arg("T"), "Finite-temperature one-loop correction")
      .def("VHT", &OneLoopPotential::VHT, py::arg("phi"), py::arg("T"), "High-temperature expansion")
      .def("daisy", py::overload_cast<Eigen::VectorXd, double>(&OneLoopPotential::daisy, py::const_),
           py::arg("phi"), py::arg("T"), "Daisy correction")
      .def("counter_term", &OneLoopPotential::counter_term, py::arg("phi"), py::arg("T"))
      .def("d2V0_dx2", &OneLoopPotential::d2V0_dx2, py::arg("phi"), "Hessian of the tree-level potential")
      .def("get_scalar_masses_sq", &OneLoopPotential::get_scalar_masses_sq, py::arg("phi"), py::arg("xi"))
      .def("get_fermion_masses_sq", &OneLoopPotential::get_fermion_masses_sq, py::arg("phi"))
      .def("get_vector_masses_sq", &OneLoopPotential::get_vector_masses_sq, py::arg("phi"))
      .def("get_ghost_masses_sq", &OneLoopPotential::get_ghost_masses_sq, py::arg("phi"), py::arg("xi"))
      .def("get_scalar_debye_sq", &OneLoopPotential::get_scalar_debye_sq,
           py::arg("phi"), py::arg("xi"), py::arg("T"))
      .def("get_scalar_thermal_sq", &OneLoopPotential::get_scalar_thermal_sq, py::arg("T"))
      .def("get_vector_debye_sq", &OneLoopPotential::get_vector_debye_sq, py::arg("phi"), py::arg("T"))
      .def("get_scalar_dofs", &OneLoopPotential::get_scalar_dofs)
      .def("get_fermion_dofs", &OneLoopPotential::get_fermion_dofs)
      .def("get_vector_dofs", &OneLoopPotential::get_vector_dofs)
      .def("get_ghost_dofs", &OneLoopPotential::get_ghost_dofs)
      .def("get_tree_scalar_masses_sq", &OneLoopPotential::get_tree_scalar_masses_sq, py::arg("phi"))
      .def("get_1l_scalar_masses_sq", &OneLoopPotential::get_1l_scalar_masses_sq, py::arg("phi"), py::arg("T"))
      .def("get_renormalization_scale", &OneLoopPotential::get_renormalization_scale)
      .def("set_renormalization_scale", &OneLoopPotential::set_renormalization_scale, py::arg("Q"))
      .def_property("renormalization_scale", &OneLoopPotential::get_renormalization_scale,
                    &OneLoopPotential::set_renormalization_scale)
      .def("get_xi", &OneLoopPotential::get_xi)
      .def("set_xi", &OneLoopPotential::set_xi, py::arg("xi"))
      .def_property("xi", &OneLoopPotential::get_xi, &OneLoopPotential::set_xi)
      .def("get_daisy_method", &OneLoopPotential::get_daisy_method)
      .def("set_daisy_method", &OneLoopPotential::set_daisy_method, py::arg("method"))
      .def_property("daisy_method", &OneLoopPotential::get_daisy_method, &OneLoopPotential::set_daisy_method);
}
