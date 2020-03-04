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

#include "potential.hpp"
#include "pow.hpp"

namespace EffectivePotential {

Eigen::VectorXd Potential::d2V_dxdt(Eigen::VectorXd phi, double T) const {
  Eigen::VectorXd gradient = Eigen::VectorXd::Zero(phi.size());

  for (int ii = 0; ii < phi.size(); ++ii) {
    Eigen::VectorXd phi_shifted = phi;
    for (int jj = 0; jj < n_h_xy.size(); ++jj) {
      phi_shifted(ii) = phi(ii) + n_h_xy[jj] * h;
      for (int kk = 0; kk < n_h_xy.size(); ++kk) {
        const double T_shifted = T + n_h_xy[kk] * h;
        gradient(ii) += V(phi_shifted, T_shifted) * coeff_xy[jj] * coeff_xy[kk] / square(h);
      }
    }
  }
  return gradient;
}

Eigen::MatrixXd Potential::d2V_dx2(Eigen::VectorXd phi, double T) const {
  Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(phi.size(), phi.size());

  // diagonal elements
  for (int ii = 0; ii < phi.size(); ++ii) {
    Eigen::VectorXd phi_shifted = phi;
    for (int jj = 0; jj < n_h_xx.size(); ++jj) {
      phi_shifted(ii) = phi(ii) + n_h_xx[jj] * h;
      hessian(ii, ii) += V(phi_shifted, T) * coeff_xx[jj] / square(h);
    }
  }

  // off-diagonal elements
  for (int ii = 0; ii < phi.size(); ++ii) {
    for (int jj = 0; jj < ii; ++jj) {
      Eigen::VectorXd phi_shifted = phi;
      for (int kk = 0; kk < n_h_xy.size(); ++kk) {
        phi_shifted(ii) = phi(ii) + n_h_xy[kk] * h;
        for (int ll = 0; ll < n_h_xy.size(); ++ll) {
          phi_shifted(jj) = phi(jj) + n_h_xy[ll] * h;
          hessian(ii, jj) += V(phi_shifted, T) * coeff_xy[kk] * coeff_xy[ll] / square(h);
        }
      }
      hessian(jj, ii) = hessian(ii, jj);
    }
  }
  return hessian;
}

}  // namespace EffectivePotential
