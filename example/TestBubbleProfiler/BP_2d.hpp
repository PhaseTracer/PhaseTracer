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

#ifndef POTENTIAL_BP_2d_HPP_INCLUDED
#define POTENTIAL_BP_2d_HPP_INCLUDED

/**
   The 2d example in BubbleProfiler :
   " ./bin/run_cmd_line_potential.x \
      --potential "(x^2 + y^2)*(1.8*(x - 1)^2 + 0.2*(y-1)^2 - 0.4)" \
      --field x --field y --local-minimum 0 0 "
 
*/

#include "potential.hpp"
#include <vector>
#include <cmath>

namespace EffectivePotential {

class BP_2d: public Potential {
 public:
  size_t get_n_scalars() const override {return 2;}

  double V(Eigen::VectorXd phi, double T) const override {
    double x = phi[0];
    double y = phi[1];
    return (x*x + y*y)*(1.8*pow(x - 1, 2) + 0.2*pow(y-1,2) - 0.4);
  }
};

}  // namespace EffectivePotential

#endif
