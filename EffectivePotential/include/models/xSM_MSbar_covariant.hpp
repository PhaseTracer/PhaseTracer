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

#ifndef XSM_MSBAR_COVARIANT_HPP_INCLUDED
#define XSM_MSBAR_COVARIANT_HPP_INCLUDED

/**
 * xSM model in covariant gauge
 */

#include <cmath>
#include <vector>

#include "SM_parameters.hpp"
#include "xSM_MSbar.hpp"
#include "pow.hpp"


namespace EffectivePotential {

/** @brief Real part of square root */
double sqrtr(double x) { return (x > 0.) ? std::sqrt(x) : 0.; }

class xSM_MSbar_covariant : public xSM_MSbar {
 public:
 
  std::vector<double> get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const override {
    return get_scalar_debye_sq(phi, xi, 0.);
  }
    
  std::vector<double> get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const override {
    const double mgb_sq = muh_sq + lambda_h * square(phi[0]) + 0.5 * lambda_hs * square(phi[1]);

    const double mode1 = (muh_sq + 0.5 * lambda_hs * square(phi[1]) + lambda_h * square(phi[0]))
                          * xi * square(SM::g) * square(phi[0]);

    const double mode2 = (muh_sq + 0.5 * lambda_hs * square(phi[1]) + lambda_h * square(phi[0]))
                          * xi * (square(SM::g) + square(SM::gp)) * square(phi[0]);

    const double pi_11 = (0.5 * lambda_h + 0.0625 * square(SM::gp) 
                          + 3. * 0.0625 * square(SM::g) + 0.25 * SM::yt_sq) * square(T);
    const double pi_22 = 0.25 * lambda_s * square(T); 

    const double m1p_sq = 0.5 * (mgb_sq + sqrtr(square(mgb_sq) - mode1)) + pi_11; 
    const double m1m_sq = 0.5 * (mgb_sq - sqrtr(square(mgb_sq) - mode1)) + pi_11; 
    const double m2p_sq = 0.5 * (mgb_sq + sqrtr(square(mgb_sq) - mode2)) + pi_11; 
    const double m2m_sq = 0.5 * (mgb_sq - sqrtr(square(mgb_sq) - mode2)) + pi_11; 

    const double m11_sq = muh_sq + 0.5 * lambda_hs * square(phi[1]) + 3. * lambda_h * square(phi[0]) + pi_11;
    const double m22_sq = mus_sq + 0.5 * lambda_hs * square(phi[0]) + 3. * lambda_s * square(phi[1]) + pi_22;
    const double m12_sq = lambda_hs * phi[1] * phi[0];

    const double trace = m11_sq + m22_sq;
    const double det = m11_sq * m22_sq - square(m12_sq);

    const double mh_sq = 0.5 * (trace - sqrtr(square(trace) -4. * det));
    const double ms_sq = 0.5 * (trace + sqrtr(square(trace) -4. * det));

    return {m1p_sq, m1m_sq, m2p_sq, m2m_sq, mh_sq, ms_sq};
  }
  
  // Charged Goldstone and plus-minus mode (2, 2)
  // Pseudo-scalar Goldsotne and plus-minus mode (1, 1)
  // Higgs and scalar singlet (1, 1,)
  std::vector<double> get_scalar_dofs() const override { return {2., 2., 1., 1., 1., 1.}; }

  // W, Z but not photon
  std::vector<double> get_vector_dofs() const override {
    return {6., 3.};
  }

  // W, Z but not photon
  std::vector<double> get_vector_debye_sq(Eigen::VectorXd phi, double T) const override {
    const double mw_sq = 11. / 6. * square(SM::g) * square(T)
                         + 0.25 * square(SM::g) * square(phi[0]);
    const double mz_sq = 11. / 6. * (pow_4(SM::g) + pow_4(SM::gp)) / (square(SM::g) + square(SM::gp)) * square(T)
                         + 0.25 * (square(SM::g) + square(SM::gp)) * square(phi[0]);
    return {mw_sq, mz_sq};
  }

  std::vector<double> get_vector_masses_sq(Eigen::VectorXd phi) const override {
    return get_vector_debye_sq(phi, 0.);
  }
};

}  // namespace EffectivePotential

#endif
