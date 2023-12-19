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

void Potential::solveBetas(std::vector<double> x0, double t0, double t_start, double t_end, double dt){
    // Define the type for odeint
    using state_type = std::vector<double>;
    // Define the state variables
    state_type x;
    std::vector<double> t_vec;   // Store the values of t
    std::vector<state_type> x_vec;  // Store the values of x
    
    t_vec.push_back(t0);
    x_vec.push_back(x0);
    
    // Define the stepper
    x = x0;
    boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper_down;    // Solve the ODEs and print the results
    for (double t = t0; t >= t_start; t -= dt)
    {
        stepper_down.do_step(
                [this](const state_type& x, state_type& dxdt, double t) {Betas(x, dxdt, t);},
                x, t, -dt);  // Solve the ODEs
        t_vec.push_back(t-dt);
        x_vec.push_back(x);
    }
    std::reverse(t_vec.begin(), t_vec.end());
    std::reverse(x_vec.begin(), x_vec.end());
    
    // Define the stepper
    x=x0;
    boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper_up;    // Solve the ODEs and print the results
    for (double t = t0; t <= t_end; t += dt)
    {
        stepper_up.do_step(
             [this](const state_type& x, state_type& dxdt, double t) {Betas(x, dxdt, t);}, 
             x, t, dt);  // Solve the ODEs
        if (t == t0) continue;
        t_vec.push_back(t+dt);
        x_vec.push_back(x);
    }
    
//    TODO: checking the result. For example, x= -inf when t_start = 10
//    TODO: dynamically adjust the range
//    for (const auto& element : t_vec) {
//        std::cout << element << " ";
//    }
//    std::cout << std::endl;
//    for (const auto& element : x_vec) {
//        std::cout << element[8] << " ";
//    }
//    std::cout << std::endl;
    
    int n = t_vec.size();
    alglib::real_1d_array t_array;
    t_array.setlength(n);
    alglib::real_1d_array x_array;
    x_array.setlength(n);
    
    for (int j = 0; j < x0.size(); j++) {
      for (int i = 0; i < t_vec.size(); i++) {
        t_array[i] = t_vec[i];
        x_array[i] = x_vec[i][j];
      }
      auto RGE = make_cubic_spline(t_array, x_array);
      RGEs.push_back(RGE);
    }
}

}  // namespace EffectivePotential
