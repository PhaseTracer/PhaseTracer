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

#ifndef PHASETRACER_PYTHON_BINDINGS_HPP_
#define PHASETRACER_PYTHON_BINDINGS_HPP_

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/iostream.h>
#include <pybind11/stl.h>

#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>

// not phasetracer.hpp: its plotting helpers are non-inline functions defined in the headers
#include "runner.hpp"

namespace py = pybind11;

void bind_potential(py::module_ &m);
void bind_models(py::module_ &m);
void bind_status(py::module_ &m);
void bind_config(py::module_ &m);
void bind_results(py::module_ &m);
void bind_stages(py::module_ &m);
void bind_runner(py::module_ &m);

/** @brief Builds a Python __repr__ from a type's operator<<. */
template <typename T>
std::string repr_from_stream(const T &value) {
  std::ostringstream out;
  out << value;
  return out.str();
}

/**
 * @brief Python-side reference to an object owned by a Runner.
 *
 * It holds a reference to the Python Runner object, so the Runner outlives it. Runner::run()
 * rebuilds every stage object, so a handle also remembers the run it came from and refuses to be
 * used after the Runner has run again, instead of dereferencing freed memory.
 */
template <typename T>
struct Handle {
  py::object owner;
  const PhaseTracer::Runner *runner;
  std::size_t run_id;
  T *obj;

  T &get() const {
    if (runner->run_id() != run_id) {
      throw std::runtime_error("stale handle: the Runner was re-run after this object was obtained");
    }
    return *obj;
  }
};

/** @brief Handle to obj, owned by the Runner wrapped by the Python object `runner`. */
template <typename T>
Handle<T> make_handle(const py::object &runner, T &obj) {
  const auto &r = runner.cast<const PhaseTracer::Runner &>();
  return Handle<T>{runner, &r, r.run_id(), &obj};
}

/** @brief Handle to another object owned by the same Runner as `parent`. */
template <typename T, typename U>
Handle<T> sibling_handle(const Handle<U> &parent, T &obj) {
  return Handle<T>{parent.owner, parent.runner, parent.run_id, &obj};
}

/** @brief Handle to one ThermalParameterSet, which is identified by its index in the ThermoFinder. */
struct ThermalParameterSetHandle {
  Handle<PhaseTracer::ThermoFinder> thermo_finder;
  std::size_t index;

  const PhaseTracer::ThermalParameterSet &get() const {
    return thermo_finder.get().thermal_parameters.at(index);
  }
};

#endif // PHASETRACER_PYTHON_BINDINGS_HPP_
