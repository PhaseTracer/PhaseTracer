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

/**
   The 2d example in BubbleProfiler :
   " ./bin/run_cmd_line_potential.x \
      --potential "(x^2 + y^2)*(1.8*(x - 1)^2 + 0.2*(y-1)^2 - 0.4)" \
      --field x --field y --local-minimum 0 0 "
 
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"
#include "shooting.hpp"
#include "path_deformation.hpp"

class BP_2d_potential: public EffectivePotential::Potential {
 public:
  size_t get_n_scalars() const override {return 2;}

  double V(Eigen::VectorXd phi, double T) const override {
    double x = phi[0];
    double y = phi[1];
    return (x*x + y*y)*(1.8*pow(x - 1, 2) + 0.2*pow(y-1,2) - 0.4);
  }
};


int main(int argc, char* argv[]) {
  
  LOGGER(debug);
  
  BP_2d_potential p;
  
  PhaseTracer::PathDeformation pd(p);
  
  
  std::vector<Eigen::VectorXd> path_pts;
  path_pts.push_back(Eigen::VectorXd(2));
  path_pts.push_back(Eigen::VectorXd(2));
  path_pts[0] << 1.046, 1.66;
  path_pts[1] << 0, 0;
  auto a = pd.full_tunneling(path_pts);
  
  LOG(debug)<< "Action = "<< std::setprecision(10) << a.action;
  LOG(debug)<< "fRatio = "<< std::setprecision(10) << a.fRatio;
  
}
