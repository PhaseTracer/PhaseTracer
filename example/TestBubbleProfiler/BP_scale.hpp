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

#ifndef POTENTIAL_BP_SCALE_HPP_INCLUDED
#define POTENTIAL_BP_SCALE_HPP_INCLUDED

/**
   The scale example in BubbleProfiler
*/

#include "potential.hpp"
#include <vector>
#include <cmath>

namespace EffectivePotential {

class BP_scale: public Potential {
 public:
  
  BP_scale(double E_,
            double alpha_,
            double scale_){
    E = E_;
    alpha = alpha_;
    scale = scale_;
  };
  
  size_t get_n_scalars() const override {return 1;}

  double V(Eigen::VectorXd phi_, double T) const override {
    double phi = phi_[0]/scale;
    return -pow(scale, 4) * E * (-alpha * phi * phi * phi * phi
                                 + phi * phi * phi
                                 + 0.5 * (4. * alpha - 3.) * phi * phi);
  }
  
  Eigen::VectorXd dV_dx(Eigen::VectorXd phi_, double T) const override {
    double phi = phi_[0]/scale;
    Eigen::VectorXd d(1);
    d(0) = -pow(scale, 3) * E * (-4. * alpha * phi * phi * phi
                                 + 3. * phi * phi
                                 + (4. * alpha - 3.) * phi);
    return d;
  }

  Eigen::MatrixXd d2V_dx2(Eigen::VectorXd phi_, double T) const  override {
    double phi = phi_[0]/scale;
    Eigen::MatrixXd d(1, 1);
    d(0, 0) = -pow(scale, 2) * E * (-12. * alpha * phi * phi
                                 + 6. * phi + (4. * alpha - 3.));
    return d;
  }
  
 private:
  PROPERTY(double, E, 1.)
  PROPERTY(double, alpha, 0.6)
  PROPERTY(double, scale, 2.)
};

}  // namespace EffectivePotential

#endif
