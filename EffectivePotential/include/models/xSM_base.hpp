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

#ifndef POTENTIAL_XSM_BASE_HPP_INCLUDED
#define POTENTIAL_XSM_BASE_HPP_INCLUDED

/**
 * Common parts of xSM models
 * See arXiv:2208.01319  [hep-ph] for details
 */

#include <algorithm>
#include <utility>
#include <vector>
#include <Eigen/Eigenvalues>

#include "one_loop_potential.hpp"
#include "pow.hpp"
#include "SM_parameters.hpp"


namespace EffectivePotential {

class xSM_base : public OneLoopPotential {
 public:

  double V0(Eigen::VectorXd phi) const override {
    return 0.5 * muh_sq * square(phi[0]) +
           0.25 * lambda_h * pow_4(phi[0]) +
           0.25 * lambda_hs * square(phi[0]) * square(phi[1]) +
           0.5 * mus_sq * square(phi[1]) +
           0.25 * lambda_s * pow_4(phi[1]);
  }
  
  /**
   * Thermal scalar masses of form c * T^2 etc for high-temperature expansion of potential
   */
  std::vector<double> get_scalar_thermal_sq(double T) const override {
    const double c_h = (9. * square(SM_g) +
                        3. * square(SM_gp) +
                        2. * (6. * SM_yt_sq + 6. * SM_yb_sq +
                              2. * SM_ytau_sq + 12. * lambda_h + lambda_hs)) / 48.;
    const double c_s = (2. * lambda_hs + 3. * lambda_s) / 12.;
    return {c_h * square(T), c_s * square(T)};
  }

  // Physical Higgs bosons and Goldstone bosons
  std::vector<double> get_scalar_dofs() const override {
    return {1., 1., 1., 1., 1};
  }
  
  // W, Z, photon
  std::vector<double> get_vector_masses_sq(Eigen::VectorXd phi) const override {
    return get_vector_debye_sq(phi, 0.);
  }

  std::vector<double> get_ghost_masses_sq(Eigen::VectorXd phi, double xi) const override {
    auto const vector_masses_sq = get_vector_debye_sq(phi, 0.);
    return {xi*vector_masses_sq[0], xi*vector_masses_sq[1], xi*vector_masses_sq[2]};
  }
  std::vector<double> get_ghost_dofs() const override {
    return {2., 1., 1.};
  }
  
  // W, Z, photon
  std::vector<double> get_vector_debye_sq(Eigen::VectorXd phi, double T) const override {
    const double h_sq = square(phi[0]);
    const double T_sq = square(T);
    const double MW_T_sq = 0.25 * square(SM_g) * h_sq;
    const double MZ_T_sq = 0.25 * (square(SM_g) + square(SM_gp)) * h_sq;
    const double Mphoton_T_sq = 0.;
    
    const double MW_L_sq = 0.25 * square(SM_g) * h_sq + 11. / 6. * square(SM_g) * T_sq;
		const double a_L = (square(SM_g) + square(SM_gp)) * (3. * h_sq + 22. * T_sq);
		const double b_L = std::sqrt(9. * square(square(SM_g) + square(SM_gp)) * square(h_sq)
                     + 132. * square(square(SM_g) - square(SM_gp)) * h_sq * T_sq
                     + 484. * square(square(SM_g) - square(SM_gp)) * pow_4(T));
		const double MZ_L_sq = (a_L + b_L) / 24.;
		const double Mphoton_L_sq = (a_L - b_L) / 24.;
    
    // Mphoton_sq must be put at the end, as it will not be used in the OSlike scheme.
    return {MW_L_sq, MZ_L_sq, Mphoton_L_sq, MW_T_sq, MZ_T_sq, Mphoton_T_sq};
  }

  // W, Z, photon
  std::vector<double> get_vector_dofs() const override {
    return {2., 1., 1., 4., 2, 2};
  }

  // top quark
  std::vector<double> get_fermion_masses_sq(Eigen::VectorXd phi) const override {
    return {0.5 * SM_yt_sq * square(phi[0]),
            0.5 * SM_yb_sq * square(phi[0]),
            0.5 * SM_ytau_sq * square(phi[0])};
  }

  // top quark
  std::vector<double> get_fermion_dofs() const override {
    return {12., 12., 4.};
  }

  size_t get_n_scalars() const override {
    return 2;
  }

  std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
    auto phi1 = phi;
    phi1[0] = - phi[0];
    auto phi2 = phi;
    phi2[1] = - phi[1];
    return {phi1, phi2};
  };
  
  void set_SM_parameters(std::vector<double> SM_parameters){
    SM_mh = SM_parameters[0];
    SM_v = SM_parameters[1];
    SM_gp = SM_parameters[2];
    SM_g = SM_parameters[3];
    SM_yt_sq = SM_parameters[4];
    SM_yb_sq = SM_parameters[5];
    SM_ytau_sq = SM_parameters[6];
  }
  
 protected:
 
   // Lagrangian parameters
  double ms;
  double lambda_hs;
  double muh_sq;
  double lambda_h;
  double mus_sq;
  double lambda_s;
  
  double SM_mh = SM::mh;
  double SM_v = SM::v;
  double SM_gp = SM::gp;
  double SM_g = SM::g;
  double SM_yt_sq = SM::yt_sq;
  double SM_yb_sq = SM::yb_sq;
  double SM_ytau_sq = SM::ytau_sq;
  
};

}  // namespace EffectivePotential

#endif
