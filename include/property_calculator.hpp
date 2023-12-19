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

#ifndef PHASETRACER_PROPERTYCALCULATOR_HPP_INCLUDED
#define PHASETRACER_PROPERTYCALCULATOR_HPP_INCLUDED

#include <ostream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <boost/cstdint.hpp>
#include "nlopt.hpp"


#include "logger.hpp"
#include "transition_finder.hpp"

namespace PhaseTracer {

class PropertyCalculator {
public:
  explicit PropertyCalculator(TransitionFinder& tf_) :
    tf(tf_) {
    std::cout << "PropertyCalculator starts " << std::endl;
  }
  virtual ~PropertyCalculator() = default;
  
private:
  TransitionFinder& tf;

  double get_action(const Eigen::VectorXd& vacuum_1, const Eigen::VectorXd& vacuum_2, double T) const;
  


};





}  // namespace PhaseTracer

#endif
