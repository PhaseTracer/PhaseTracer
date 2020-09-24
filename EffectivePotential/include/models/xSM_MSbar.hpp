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

#ifndef POTENTIAL_XSM_MSbar_HPP_INCLUDED
#define POTENTIAL_XSM_MSbar_HPP_INCLUDED

/**
 * Z2 symmetric real scalar singlet extension of the Standard Model
 *
 * Conventions match 1808.01098.
 */

#include <vector>

#include "xSM_MSbar_solver.hpp"
#include "one_loop_potential.hpp"
#include "pow.hpp"
#include "SM_parameters.hpp"


namespace EffectivePotential {

class xSM_MSbar : public OneLoopPotential {
 public:
  /**
   * @brief Make an xSM model from Lagrangian parameters
   */
  xSM_MSbar(double lambda_hs_,
            double muh_sq_,
            double lambda_h_,
            double mus_sq_,
            double lambda_s_,
            double muh_sq_tree_EWSB_,
            double mus_sq_tree_EWSB_):
    lambda_hs(lambda_hs_), muh_sq(muh_sq_),
    lambda_h(lambda_h_), mus_sq(mus_sq_), lambda_s(lambda_s_),
    muh_sq_tree_EWSB(muh_sq_tree_EWSB_), mus_sq_tree_EWSB(mus_sq_tree_EWSB_) {}

  /**
   * @brief Make an xSM model using 1808.01098 conventions
   *
   * @param lambda_hs Quartic coupling
   * @param Q Renormalization scale
   * @param tree_level Whether to use tree-level tadpole conditions
   */
  xSM_MSbar(double lambda_hs_, double Q, bool tree_level = true) : lambda_hs(lambda_hs_) {
    // NB factors of 2 and 4 etc due to different conventions in solver

    // Match choices in 1808.01098
    const double ms = 0.5 * SM::mh;
    const double lambda_s_min = 2. / square(SM::mh * SM::v) *
      square(square(ms) - 0.5 * lambda_hs * square(SM::v));
    lambda_s = lambda_s_min + 0.1;

    // Solve constraints on Higgs, singlet masses etc
    xSM_MSbar_parameters_solver solver(ms, 0, 0.25 * lambda_hs, 0, 0, 0.25 * lambda_s);
    solver.set_renormalization_scale(Q);
    const bool valid = solver.solving_parameters();
    if (!valid) {
      throw std::runtime_error("Invalid when solving parameters");
    }

    // Extract parameters in appropriate convention
    muh_sq = 2. * (tree_level ? solver.get_mu_h_Sq_tree() : solver.get_mu_h_Sq());
    lambda_h = 4. * (tree_level ? solver.get_lambda_h_tree() : solver.get_lambda_h());
    mus_sq = 2. * (tree_level ? solver.get_mu_s_Sq_tree() : solver.get_mu_s_Sq());
    muh_sq_tree_EWSB = tree_level ? 2. * solver.get_mu_h_Sq_tree_EWSB() : muh_sq;
    mus_sq_tree_EWSB = tree_level ? 2. * solver.get_mu_s_Sq_tree_EWSB() : mus_sq;

    // Set renormalization scheme consistently with tadpoles
    set_renormalization_scale(Q);
  }

  double get_v_tree_s() const {
    return std::sqrt(-mus_sq / lambda_s);
  }

  double V0(Eigen::VectorXd phi) const override {
    return 0.5 * muh_sq * square(phi[0]) +
           0.25 * lambda_h * pow_4(phi[0]) +
           0.25 * lambda_hs * square(phi[0]) * square(phi[1]) +
           0.5 * mus_sq * square(phi[1]) +
           0.25 * lambda_s * pow_4(phi[1]);
  }

  std::vector<double> get_scalar_thermal_sq(double T) const override {
    const double c_h = (9. * square(SM::g) +
                        3. * square(SM::gp) +
                        2. * (6. * SM::yt_sq + 12. * lambda_h + lambda_hs)) / 48.;
    const double c_s = (2. * lambda_hs + 3. * lambda_s) / 12.;
    return {c_h * square(T), c_s * square(T)};
  }

  std::vector<double> get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const override {
    const double h = phi[0];
    const double s = phi[1];

    // Higgs and singlet. Construct elements of mass matrix
    const double Mhh = muh_sq_tree_EWSB + 3. * lambda_h * square(h) + 0.5 * lambda_hs * square(s);
    const double Mss = mus_sq_tree_EWSB + 3. * lambda_s * square(s) + 0.5 * lambda_hs * square(h);
    const double Mhs = lambda_hs * h * s;
    // diagonalize mass matrix
    const double A = 0.5 * (Mhh+Mss);
    const double B = std::sqrt(0.25 * square(Mhh - Mss) + square(Mhs));
    const double mH1 = A - B;
    const double mH2 = A + B;

    // Goldstone bosons
    const double mg = muh_sq_tree_EWSB + lambda_h * square(h) + 0.5 * lambda_hs * square(s);

    return {mH1, mH2, mg};
  }

  // Physical Higgs bosons and Goldstone bosons
  std::vector<double> get_scalar_dofs() const override {
    return {1., 1., 3.};
  }

  // W, Z, photon
  std::vector<double> get_vector_masses_sq(Eigen::VectorXd phi) const override {
    const double h_sq = square(phi[0]);
    const double MW_sq = 0.25 * square(SM::g) * h_sq;
    const double MZ_sq = 0.25 * (square(SM::g) + square(SM::gp)) * h_sq;
    const double MG_sq = 0;
    return {MW_sq, MZ_sq, MG_sq};
  }

  // W, Z, photon
  std::vector<double> get_vector_dofs() const override {
    return {6., 3., 2.};
  }

  // top quark
  std::vector<double> get_fermion_masses_sq(Eigen::VectorXd phi) const override {
    return {0.5 * SM::yt_sq * square(phi[0])};
  }

  // top quark
  std::vector<double> get_fermion_dofs() const override {
    return {12.};
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

 private:
  // Lagrangian parameters
  double lambda_hs;
  double muh_sq;
  double lambda_h;
  double mus_sq;
  double lambda_s;
  // For consistency in one-loop potential
  double muh_sq_tree_EWSB;
  double mus_sq_tree_EWSB;
};

}  // namespace EffectivePotential

#endif
