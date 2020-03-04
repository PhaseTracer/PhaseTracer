#ifndef PHASETRACER_OVERLOAD_HPP_INCLUDED
#define PHASETRACER_OVERLOAD_HPP_INCLUDED

// TODO do these work for boost log?

#include <Eigen/Core>
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
