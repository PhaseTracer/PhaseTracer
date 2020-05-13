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

#ifndef PHASETRACER_OVERLOAD_HPP_INCLUDED
#define PHASETRACER_OVERLOAD_HPP_INCLUDED

// TODO do these work for boost log?

#include <eigen3/Eigen/Core>
#include <ostream>
#include <iterator>
#include <vector>
#include <algorithm>

namespace PhaseTracer {

/** Overload stream for vectors */
template <typename T>
std::ostream& operator<< (std::ostream& o, std::vector<T> v) {
  o << std::boolalpha;
  o << "[";
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(o, ", "));
  o << "\b\b]";
  return o;
}

/** Overload for writing our field coordinates */
inline std::ostream& operator<< (std::ostream& o, Eigen::VectorXd v) {
  static Eigen::IOFormat format(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", ", ", "", "", "[", "]");
  o << v.format(format);
  return o;
}

}  // namespace PhaseTracer

#endif
