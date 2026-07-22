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
  Single-field scalar Yukawa theory in d = 3, Example 1 of
  Phys.Rev.D 104 (2021) 096015 (arXiv:2108.04377). The finite-temperature
  potential is obtained from leading-order high-temperature dimensional
  reduction, so V(phi, T) below is the 3d effective potential and the relevant
  bounce is the O(3) (d = 3) instanton.

  This is the C++ port of the BubbleDet reference script (the "BubbleDet" cell of
  example/TestThermalParameters/test_helpers.ipynb).
*/

#ifndef EFFECTIVEPOTENTIAL_SCALAR_YUKAWA_3D_HPP_
#define EFFECTIVEPOTENTIAL_SCALAR_YUKAWA_3D_HPP_

#include <cmath>

#include <Eigen/Core>

#include "potential.hpp"
#include "property.hpp"
#include "pow.hpp"

namespace EffectivePotential {

class ScalarYukawa3D : public Potential {
public:
  double V(Eigen::VectorXd phi, double T) const override {
    const double x = phi[0];
    const Params3d p = params_3d(T);
    return p.s3 * x +
           0.5 * p.m3sq * square(x) +
           1. / 6. * p.g3 * cube(x) +
           1. / 24. * p.lam3 * pow_4(x);
  }

  size_t get_n_scalars() const override { return 1; }

  Eigen::VectorXd dV_dx(Eigen::VectorXd phi, double T) const override {
    const double x = phi[0];
    const Params3d p = params_3d(T);
    Eigen::VectorXd d(1);
    d(0) = p.s3 + p.m3sq * x + 0.5 * p.g3 * square(x) + 1. / 6. * p.lam3 * cube(x);
    return d;
  }

  Eigen::MatrixXd d2V_dx2(Eigen::VectorXd phi, double T) const override {
    const double x = phi[0];
    const Params3d p = params_3d(T);
    Eigen::MatrixXd d(1, 1);
    d(0, 0) = p.m3sq + p.g3 * x + 0.5 * p.lam3 * square(x);
    return d;
  }

  /** Coefficients of the 3d potential as polynomials in x. Public so a driver
   *  can find the minima directly (as the reference script does with np.roots),
   *  matching the notebook's `params_3d`. */
  struct Params3d {
    double s3, m3sq, g3, lam3;
  };
  Params3d params_3d(double T) const {
    return Params3d{
        (s + 1. / 24. * (g + 4. * y * mf) * square(T)) / std::sqrt(T),
        msq + 1. / 24. * (lam + 4. * square(y)) * square(T),
        std::sqrt(T) * g,
        T * lam};
  }

private:
  // Zero-temperature parameters (arXiv:2108.04377, Example 1).
  PROPERTY(double, s, 0.0)
  PROPERTY(double, msq, -1.0)
  PROPERTY(double, g, 0.30)
  PROPERTY(double, lam, 0.1)
  PROPERTY(double, mf, -0.2)
  PROPERTY(double, y, 0.3)
};

} // namespace EffectivePotential

#endif // EFFECTIVEPOTENTIAL_SCALAR_YUKAWA_3D_HPP_
