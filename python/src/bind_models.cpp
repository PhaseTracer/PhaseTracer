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

// The models are bound one family per file, so that each model header is included in a single
// translation unit. A new family file needs its bind_models_<family>() declared in bindings.hpp and
// called here; CMake picks the file up when it is re-run.
void bind_models(py::module_ &m) {
  bind_models_test(m);
  bind_models_xsm(m);
  bind_models_rss(m);
  bind_models_bsm(m);
  bind_models_thdm(m);
}
