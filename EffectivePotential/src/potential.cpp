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

Eigen::VectorXd Potential::dV_dx(Eigen::VectorXd phi, double T) const {
  Eigen::VectorXd gradient = Eigen::VectorXd::Zero(phi.size());
  
  for (int ii = 0; ii < phi.size(); ++ii) {
    Eigen::VectorXd phi_shifted = phi;
    for (int jj = 0; jj < n_h_xy.size(); ++jj) {
      phi_shifted(ii) = phi(ii) + n_h_xy[jj] * h;
      gradient(ii) += V(phi_shifted, T) * coeff_xy[jj] / h;
    }
  }
  return gradient;
}

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

alglib::spline1dinterpolant Potential::make_cubic_spline(alglib::real_1d_array x, alglib::real_1d_array y) {
    alglib::spline1dinterpolant spline_;
    alglib::spline1dbuildcubic(x, y, spline_);
    return spline_;
}

void Potential::solveBetas(std::vector<double> x0, double t0, double t_start, double t_end) {
    using state_type = std::vector<double>;

    state_type x = x0;
    std::vector<double> t_vec;
    std::vector<state_type> x_vec;

    // Stepper with constant step size
    using stepper_type = boost::numeric::odeint::runge_kutta4<state_type>; // Runge-Kutta 4th order with constant step size
    stepper_type stepper;

    double dt = 1; // Define a constant step size

    // Integrate backward from t0 to t_start
    boost::numeric::odeint::integrate_const(
        stepper,
        [this](const state_type& x, state_type& dxdt, double t) { Betas(x, dxdt, t); },
        x, t0, t_start, -dt,
        [&](const state_type& x, double t) {
            t_vec.push_back(t);
            x_vec.push_back(x);

        });

    // Reverse the vectors for the downward integration
    std::reverse(t_vec.begin(), t_vec.end());
    std::reverse(x_vec.begin(), x_vec.end());

    // Reinitialize for forward integration
    x = x0;

    boost::numeric::odeint::integrate_const(
        stepper,
        [this](const state_type& x, state_type& dxdt, double t) { Betas(x, dxdt, t); },
        x, t0, t_end, dt,
        [&](const state_type& x, double t) {
            if (t == t0) return; // Skip the initial point
            t_vec.push_back(t);
            x_vec.push_back(x);

        });

    // Since we used a constant step size, no need to re-sample
    // Proceed directly to creating the final spline on the uniform grid
    for (int j = 0; j < x0.size(); j++) {
        alglib::real_1d_array t_array, x_array;
        t_array.setlength(t_vec.size());
        x_array.setlength(t_vec.size());

        for (int i = 0; i < t_vec.size(); i++) {
            t_array[i] = t_vec[i];
            x_array[i] = x_vec[i][j];
        }

        // Perform cubic spline interpolation on the original data
        auto RGE = make_cubic_spline(t_array, x_array);

        RGEs.push_back(RGE);
    }
}

}  // namespace EffectivePotential
