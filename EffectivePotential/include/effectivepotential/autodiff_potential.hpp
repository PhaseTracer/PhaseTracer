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

#ifndef POTENTIAL_AUTODIFF_POTENTIAL_HPP_INCLUDED
#define POTENTIAL_AUTODIFF_POTENTIAL_HPP_INCLUDED

#include <vector>
#include <Eigen/Core>

#include <autodiff/reverse.hpp>
#include <autodiff/reverse/eigen.hpp>

#include "property.hpp"
#include "potential.hpp"

namespace EffectivePotential {

class AutoDiffPotential : public Potential {
 public:
  AutoDiffPotential() = default;
  /** Potential, possibly at finite-temperature implemented from dual */
  virtual autodiff::var V(Eigen::VectorXvar phi, autodiff::var T) const = 0;
  /** Potential, possibly at finite-temperature implemented from dual */
  double V(Eigen::VectorXd phi, double T) const;
  /** The derivative of the gradient of potential with respect to temperature */
  virtual Eigen::VectorXd d2V_dxdt(Eigen::VectorXd phi, double T) const;
  /** The Hessian matrix of the one-loop potential at finite temperature */
  virtual Eigen::MatrixXd d2V_dx2(Eigen::VectorXd phi, double T) const;
};

}  // namespace EffectivePotential

#endif
