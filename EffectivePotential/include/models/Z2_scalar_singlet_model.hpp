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

#ifndef POTENTIAL_Z2_SCALAR_SINGLET_MODEL_HPP_INCLUDED
#define POTENTIAL_Z2_SCALAR_SINGLET_MODEL_HPP_INCLUDED

/**
   Z2 symmetric real scalar singlet extension of the Standard Model
   https://arxiv.org/pdf/1611.02073.pdf
*/

#include <vector>

#include "potential.hpp"
#include "pow.hpp"

namespace EffectivePotential {

class Z2ScalarSingletModel : public Potential {
public:
  void set_m_s(double m_s_) { m_s = m_s_; }
  void set_lambda_hs(double lambda_hs_) { lambda_hs = lambda_hs_; }

  double V(Eigen::VectorXd phi, double T) const override {
    const double mus_sq = square(m_s) - lambda_hs * square(v) / 2.;
    const double muh_sq_T = muh_sq + square(T) / 48. *
                                         (9. * g_sq + 3. * gp_sq + 12. * yt_sq + 24. * lambda_h + 2. * lambda_hs);
    const double mus_sq_T = mus_sq + square(T) / 12. * (2. * lambda_hs + 3. * lambda_s);
    return 0.5 * muh_sq_T * square(phi[0]) +
           0.25 * lambda_h * pow_4(phi[0]) +
           0.25 * lambda_hs * square(phi[0]) * square(phi[1]) +
           0.5 * mus_sq_T * square(phi[1]) +
           0.25 * lambda_s * pow_4(phi[1]);
  }

  size_t get_n_scalars() const override { return 2; }

  // Using expressions in arXiv:1611.02073
  double get_TC_from_expression() const {
    const double cs = 1. / 12. * (2. * lambda_hs + 3. * lambda_s);
    const double ch = 1. / 48. *
                      (9. * g_sq + 3. * gp_sq + 12. * yt_sq + 24. * lambda_h + 2. * lambda_hs);
    const double mus_sq = square(m_s) - lambda_hs * square(v) / 2.;
    const double TC_sq = -(lambda_s * ch * muh_sq - lambda_h * cs * mus_sq +
                           std::sqrt(lambda_s * lambda_h) * std::abs(cs * muh_sq - ch * mus_sq)) /
                         (lambda_s * square(ch) - lambda_h * square(cs));
    return std::sqrt(TC_sq);
  }

  double get_vs_from_expression() const {
    const double TC = get_TC_from_expression();
    const double mus_sq = square(m_s) - lambda_hs * square(v) / 2.;
    const double mus_sq_T = mus_sq + square(TC) / 12. * (2. * lambda_hs + 3. * lambda_s);
    return std::sqrt(-mus_sq_T / lambda_s);
  }

  double get_vh_from_expression() const {
    const double TC = get_TC_from_expression();
    const double muh_sq_T = muh_sq + square(TC) / 48. *
                                         (9. * g_sq + 3. * gp_sq + 12. * yt_sq + 24. * lambda_h + 2. * lambda_hs);
    return std::sqrt(-muh_sq_T / lambda_h);
  }

  bool check() const {
    if (lambda_hs < 2. * std::sqrt(lambda_s * lambda_h)) {
      return false;
    }

    const double mus_sq = square(m_s) - lambda_hs * square(v) / 2.;
    if (mus_sq > 0) {
      return false;
    }

    const double TC = get_TC_from_expression();
    if (std::isnan(TC)) {
      return false;
    }
    const double cs = 1. / 12. * (2. * lambda_hs + 3. * lambda_s);
    const double ch = 1. / 48. *
                      (9. * g_sq + 3. * gp_sq + 12. * yt_sq + 24. * lambda_h + 2. * lambda_hs);
    if (lambda_hs * (mus_sq + cs * square(TC)) > 2. * lambda_s * (muh_sq + ch * square(TC))) {
      return false;
    }

    return true;
  }

  std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
    auto phi1 = phi;
    phi1[0] = -phi[0];
    auto phi2 = phi;
    phi2[1] = -phi[1];
    return {phi1, phi2};
  };

private:
  const double v = 246.;
  const double mh = 125.;
  const double g_sq = 0.652 * 0.652;
  const double gp_sq = 0.352 * 0.352;
  const double yt_sq = 0.995 * 0.995;
  const double muh_sq = -0.5 * mh * mh;
  const double lambda_h = -muh_sq / square(v);
  const double lambda_s = 0.1;

  double lambda_hs = 0.25;
  double m_s = 27;
};

} // namespace EffectivePotential

#endif
