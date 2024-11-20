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

#ifndef POTENTIAL_POTENTIAL_HPP_INCLUDED
#define POTENTIAL_POTENTIAL_HPP_INCLUDED

#include <vector>
#include <Eigen/Core>
#include <boost/numeric/odeint.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <interpolation.h>

#include "property.hpp"

namespace EffectivePotential {

class Potential {
 public:
  Potential() {
    set_h_4(h_4);
  }
  /** Potential, possibly at finite-temperature */
  virtual double V(Eigen::VectorXd phi, double T) const = 0;
  /** Number of scalar fields */
  virtual size_t get_n_scalars() const = 0;
  /** Field values that are forbidden */
  virtual bool forbidden(Eigen::VectorXd phi) const { return false; };
  /** Apply symmetry operation on field */
  virtual std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const { return {}; };
  /** The gradient of potential*/
  virtual Eigen::VectorXd dV_dx(Eigen::VectorXd phi, double T) const;
  /** The derivative of the gradient of potential with respect to temperature */
  virtual Eigen::VectorXd d2V_dxdt(Eigen::VectorXd phi, double T) const;
  /** The Hessian matrix of the one-loop potential at finite temperature */
  virtual Eigen::MatrixXd d2V_dx2(Eigen::VectorXd phi, double T) const;

  /** Functor that returns potential */
  double operator () (Eigen::VectorXd phi, double T) const { return V(phi, T); }

  /** Set approximation order of derivatives */
  void set_h_4(bool h_4_) {
    h_4 = h_4_;
    n_h_xy = h_4 ? n_h_xy_4 : n_h_xy_2;
    coeff_xy = h_4 ? coeff_xy_4 : coeff_xy_2;
    n_h_xx = h_4 ? n_h_xx_4 : n_h_xx_2;
    coeff_xx = h_4 ? coeff_xx_4 : coeff_xx_2;
  }

  /* For 3D EFT*/
  virtual void Betas(const std::vector<double>& x, std::vector<double>& dxdt, const double t){};
  
 protected:
  /** The step-size in all numerical derivatives */
  PROTECTED_PROPERTY(double, h, 0.001);

  /** Whether to use a fourth-order approximation */
  PROTECTED_PROPERTY_CUSTOM_SETTER(bool, h_4, false);

  /** Steps in units of h for diagonal elements of Hessian */
  const std::vector<double> n_h_xx_4 = {-2., -1., 0., 1., 2.};
  const std::vector<double> n_h_xx_2 = {-1., 0., 1.};
  PROTECTED_PROPERTY_CUSTOM_SETTER(std::vector<double>, n_h_xx, {});

  /** Steps in units of h for off-diagonal elements of Hessian and gradient */
  const std::vector<double> n_h_xy_4 = {-2., -1., 1., 2.};
  const std::vector<double> n_h_xy_2 = {-1., 1.};
  PROTECTED_PROPERTY_CUSTOM_SETTER(std::vector<double>, n_h_xy, {});

  /** Convolution mask for diagonal elements of Hessian and gradient */
  const std::vector<double> coeff_xx_4 = {
      -1. / 12.,
      16. / 12. ,
      -30. / 12.,
      16. / 12.,
      -1. / 12.};
  const std::vector<double> coeff_xx_2 = {1., -2., 1.};
  PROTECTED_PROPERTY_CUSTOM_SETTER(std::vector<double>, coeff_xx, {});

  /** Convolution mask for off-diagonal elements of Hessian and gradient */
  const std::vector<double> coeff_xy_4 = {
    1. / 12.,
    -8. / 12.,
    8. / 12.,
    -1. / 12.};
  const std::vector<double> coeff_xy_2 = {-0.5, 0.5};
  PROTECTED_PROPERTY_CUSTOM_SETTER(std::vector<double>, coeff_xy, {});
  
  /*For 3d EFT */
  const double EulerGamma = 0.5772156649;
  const double Glaisher = 1.2824271291006226369;
  std::vector<alglib::spline1dinterpolant> RGEs;
  alglib::spline1dinterpolant make_cubic_spline(alglib::real_1d_array x, alglib::real_1d_array y);
  void solveBetas(std::vector<double> x0, double t0=100., double t_start = 20.0, double t_end = 5000.0, double dt = 1);
  
};

}  // namespace EffectivePotential

#endif
