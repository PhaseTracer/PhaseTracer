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

#include "autodiff_potential.hpp"

namespace EffectivePotential {

double AutoDiffPotential::V(Eigen::VectorXd phi, double T) const {
  Eigen::VectorXvar phi_var = phi; 
  autodiff::var T_var = T;
  autodiff::var V_var = V(phi_var, T_var);
  return autodiff::val(V_var);
}

Eigen::MatrixXd AutoDiffPotential::d2V_dx2(Eigen::VectorXd phi, double T) const {
  Eigen::VectorXvar phi_var = phi; 
  autodiff::var T_var = T;
  autodiff::var u = V(phi_var, T_var);
  Eigen::VectorXd g;  // the gradient vector
  return autodiff::hessian(u, phi_var, g);
}

Eigen::VectorXd AutoDiffPotential::d2V_dxdt(Eigen::VectorXd phi, double T) const {
  // make vector containing both phi and T
  const int n_scalars = get_n_scalars();
  Eigen::VectorXvar phi_T_var(n_scalars + 1);
  phi_T_var.head(n_scalars) = phi;
  phi_T_var[n_scalars] = T; 

  // compute hessian
  autodiff::var u = V(phi_T_var.head(n_scalars), phi_T_var[n_scalars]);
  Eigen::VectorXd g;  // the gradient vector
  auto H = autodiff::hessian(u, phi_T_var, g);

  // return only relevant row of hessian
  return H.row(H.rows() - 1);
}

}  // namespace EffectivePotential
