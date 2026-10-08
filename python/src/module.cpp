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

#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sev = boost::log::trivial;

namespace {

void set_log_level(sev::severity_level level) {
  boost::log::core::get()->set_filter(sev::severity >= level);
}

sev::severity_level parse_log_level(const std::string &name) {
  if (name == "trace") return sev::trace;
  if (name == "debug") return sev::debug;
  if (name == "info") return sev::info;
  if (name == "warning") return sev::warning;
  if (name == "error") return sev::error;
  if (name == "fatal") return sev::fatal;
  throw std::invalid_argument("unknown log level '" + name +
                              "' (use trace, debug, info, warning, error or fatal)");
}

} // namespace

PYBIND11_MODULE(_phasetracer, m) {
  m.doc() = "Python interface to PhaseTracer; import it as `phasetracer`";

  // order matters: types must be registered before they appear in signatures and defaults
  bind_potential(m);
  bind_status(m);
  bind_config(m);
  bind_results(m);
  bind_stages(m);
  bind_runner(m);
  py::module_ models = m.def_submodule("models", "Effective potentials shipped with PhaseTracer");
  bind_models(models);

  m.def("set_log_level", &set_log_level, py::arg("level"), "Set the global log level (a LogLevel)");
  m.def("set_log_level", [](const std::string &name) { set_log_level(parse_log_level(name)); }, py::arg("level"),
        "Set the global log level by name: trace, debug, info, warning, error or fatal");

  m.def("set_num_threads",
        [](int n) {
#ifdef _OPENMP
          omp_set_num_threads(n);
#else
          (void)n;
#endif
        },
        py::arg("n"), "Number of OpenMP threads used by the action calculations");
  m.def("get_max_threads", []() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
  });

  // quiet by default, as in the C++ examples
  set_log_level(sev::fatal);
}
