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

#include "xSM_base.hpp"
#include "logger.hpp"

#include <utility>
#include <vector>
#include <iostream>

namespace EffectivePotential {

class xSM_MSbar : public xSM_base {
 public:
  /**
   * @brief Make an xSM model from Lagrangian parameters
   */
  xSM_MSbar(double lambda_hs_,
            double lambda_s_,
            double ms_) {
      lambda_hs = lambda_hs_;
      lambda_s = lambda_s_;
      ms = ms_;
    }

  /**
   * @brief Make an xSM model using tree or one-loop tadpole constraints
   */
  static xSM_MSbar from_tadpoles(double lambda_hs, double lambda_s, double ms,
                                 double Q, double xi,
                                 bool use_covariant_gauge = false,
                                 bool use_1L_EWSB_in_0L_mass = false,
                                 bool use_Goldstone_resum = true,
                                 bool use_tree_level_tadpole = false,
                                 std::vector<double> SM_parameters = {}) {
    // Make model that we'll begin setting options etc on
    xSM_MSbar model(lambda_hs, lambda_s, ms);

    // Set SM parameters if present
    if (SM_parameters.size() == 7) {
      model.set_SM_parameters(SM_parameters);
    }

    // Set some options
    model.set_use_1L_EWSB_in_0L_mass(use_1L_EWSB_in_0L_mass);
    model.set_use_Goldstone_resum(use_Goldstone_resum);
    model.set_use_covariant_gauge(use_covariant_gauge);
    model.set_renormalization_scale(Q);

    // Special care about gauge. Covariant gauge only has xi dependence
    // in the masses; no explicit dependence in the potential.
    // So the ordinary xi is set to zero, as this one appears expicitly in the
    // potential. And the user input one is saved to a new parameter
    if (!use_covariant_gauge) {
      model.set_xi(xi);
    } else {
      model.xi_covariant_internal = xi;
      model.set_xi(0.);
    }

    // Apply the relevant tadpole conditions
    if (use_tree_level_tadpole) {
      model.apply_tree_level();
    } else {
      model.iteration_converged = model.apply_one_loop();
    }
    return model;
  }

  /**
   * Apply tree-level Higgs VEV, and Higgs and singlet mass to fix
   * three Lagrangian parameters.
   *
   * The singlet mass and quartic are fixed following 1808.01098.
   */
  void apply_tree_level() {
    double mhh2 = square(SM_mh);
    double mss2 = square(ms);

    // TODO: will this affect calculation of lambda_h?
    if (mss2 > mhh2) {
      std::swap(mhh2, mss2);
    }

    // Apply SM vacuum and Higgs and singlet masses to constraint three parameters
    lambda_h = mhh2 / (2. * square(SM_v));
    muh_sq = -lambda_h * square(SM_v);
    mus_sq = mss2 - 0.5 * lambda_hs * square(SM_v);
    muh_sq_use_0L_EWSB = muh_sq;
  }

  /**
   * Apply one-level Higgs VEV, and Higgs and singlet mass to fix
   * three Lagrangian parameters.
   *
   * The singlet mass and quartic are fixed following 1808.01098.
   *
   * This is an iterative solver that stops once the absolute change in
   * Lagrangian parameters is small.
   */
  bool apply_one_loop(double tol = 0.1) {
    apply_tree_level();
    muh_sq_use_0L_EWSB = muh_sq;  // TODO: this can be deleted?

    size_t ii = 0;
    while (true) {
      ++ii;
      double lambda_h_prev = lambda_h;
      double lambda_s_prev = lambda_s;
      double muh_sq_prev = muh_sq;
      double mus_sq_prev = mus_sq;

      iterate_one_loop();

      const double dmuh = std::abs(muh_sq - muh_sq_prev);
      const double dmus = std::abs(mus_sq - mus_sq_prev);
      const double dlambda_h = std::abs(lambda_h - lambda_h_prev);
      const double dlambda_s = std::abs(lambda_s - lambda_s_prev);

      const bool converged = (dmuh < tol) && (dlambda_h < tol) && (dmus < tol) && (dlambda_s < tol);
      if (converged) {
        return true;
      }
      if (ii > 100000) {
        LOG(fatal) << "1l iterations did not converge";
        return false;
      }
    }
  }

  Eigen::MatrixXd d2V1_dx2(Eigen::VectorXd phi) const {
    Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(phi.size(), phi.size());

    // diagonal elements
    for (int ii = 0; ii < phi.size(); ++ii) {
      Eigen::VectorXd phi_shifted = phi;
      for (int jj = 0; jj < n_h_xx.size(); ++jj) {
        phi_shifted(ii) = phi(ii) + n_h_xx[jj] * h;
        hessian(ii, ii) += V1(phi_shifted) * coeff_xx[jj] / square(h);
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
            hessian(ii, jj) += V1(phi_shifted) * coeff_xy[kk] * coeff_xy[ll] / square(h);
          }
        }
        hessian(jj, ii) = hessian(ii, jj);
      }
    }
    return hessian;
  }

  Eigen::VectorXd dV1_dx(Eigen::VectorXd phi) const {
    Eigen::VectorXd jacobian = Eigen::VectorXd::Zero(phi.size());

    for (int ii = 0; ii < phi.size(); ++ii) {
      Eigen::VectorXd f = phi;
      Eigen::VectorXd b = phi;
      f(ii) = phi(ii) + 0.5 * h;
      b(ii) = phi(ii) - 0.5 * h;
      jacobian(ii) = (V1(f) - V1(b)) / h;
    }

    return jacobian;
  }

  Eigen::VectorXd dV0_dx(Eigen::VectorXd phi) const {
    Eigen::VectorXd dV0dx = Eigen::VectorXd::Zero(phi.size());
    for (int ii = 0; ii < phi.size(); ++ii) {
      Eigen::VectorXd f = phi;
      Eigen::VectorXd b = phi;
      f(ii) = phi(ii) + 0.5 * h;
      b(ii) = phi(ii) - 0.5 * h;
      dV0dx(ii) = (V0(f) - V0(b)) / h;
    }
    return dV0dx;
  }

  Eigen::VectorXd dV1T_dx(Eigen::VectorXd phi, double T) const {
    Eigen::VectorXd dV1Tdx = Eigen::VectorXd::Zero(phi.size());
    for (int ii = 0; ii < phi.size(); ++ii) {
      Eigen::VectorXd f = phi;
      Eigen::VectorXd b = phi;
      f(ii) = phi(ii) + 0.5 * h;
      b(ii) = phi(ii) - 0.5 * h;
      dV1Tdx(ii) = (V1T(f, T) - V1T(b, T)) / h;
    }
    return dV1Tdx;
  }

  double dV1T_dT(Eigen::VectorXd phi, double T) const {
    const double eT =  0.001;
    return (V1T(phi, T+eT) - V1T(phi, T) )/ eT;
  }

  double ddaisy_dT(Eigen::VectorXd phi, double T) const {
    const double eT =  0.001;
    return (daisy(phi, T+eT) - daisy(phi, T) )/ eT;
  }

  /**
   * Single iteration of one-loop tadpole solver.
   *
   * Uses numerical derivatives of Coleman-Weinberg potential.
   */
  void iterate_one_loop() {
    Eigen::Vector2d vacuum;
    vacuum << SM_v, 0.;
    const auto jacobian = dV1_dx(vacuum);
    const auto hessian = d2V1_dx2(vacuum);

    double mhh2 = square(SM_mh);
    double mss2 = square(ms);

    // TODO: will this affect calculation of lambda_h?
    if (mhh2 < mss2) {
      std::swap(mhh2, mss2);
    }

    // Apply SM vacuum and Higgs and singlet masses to constraint three parameters
    lambda_h = (mhh2 + jacobian(0) / SM_v - hessian(0, 0)) / (2. * square(SM_v));
    muh_sq = -0.5 * mhh2 - 1.5 *jacobian(0) / SM_v + 0.5 * hessian(0, 0);
    mus_sq = mss2 - 0.5 * lambda_hs * square(SM_v) - hessian(1, 1);

    // Calculate muh_sq using tree level EWSB, for masses in CW potential
    muh_sq_use_0L_EWSB = - lambda_h * square(SM_v);
  }

  double get_v_tree_s() const {
    return std::sqrt(-mus_sq / lambda_s);
  }

  /**
   * Tree-level scalar masses including xi-dependence.
   *
   * These masses enter the Coleman-Weinberg potential.
   */
  std::vector<double> get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const override {
    return get_scalar_debye_sq(phi, xi, 0.);
  }

  /** @brief Real part of square root */
  double real_sqrt(double x) const {
    return (x > 0.) ? std::sqrt(x) : 0.;
  }

  /** @brief Scalar Debye masses careful treatment of covariant gauge etc */
  std::vector<double> get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const override {
    const double h = phi[0];
    const double s = phi[1];

    const double chosen_muh_sq = use_1L_EWSB_in_0L_mass ? muh_sq : muh_sq_use_0L_EWSB;
    const auto thermal_sq = get_scalar_thermal_sq(T);

    // Mass matrix elements
    const double m11_sq = chosen_muh_sq + 0.5 * lambda_hs * square(s) + 3. * lambda_h * square(h);
    const double m22_sq = mus_sq + 0.5 * lambda_hs * square(h) + 3. * lambda_s * square(s);
    const double m12_sq = lambda_hs * s * h;

    // Resummed Goldstone contributions
    const auto fm_sq = get_fermion_masses_sq(phi);
    const auto vm_sq = get_vector_masses_sq(phi);
    const double q_sq = square(get_renormalization_scale());
    const double sum = 1. / (16. * M_PI * M_PI) * (
                      + 3. * lambda_h * (q_sq * xlogx(m11_sq / q_sq) - m11_sq)
                      + 0.5 * lambda_hs * (q_sq * xlogx(m22_sq / q_sq) - m22_sq)
                      - 6. * SM_yt_sq * (q_sq * xlogx(fm_sq[0] / q_sq) - fm_sq[0])
                      - 6. * SM_yb_sq * (q_sq * xlogx(fm_sq[1] / q_sq) - fm_sq[1])  //  TODO: Need check
                      - 2. * SM_ytau_sq * (q_sq * xlogx(fm_sq[2] / q_sq) - fm_sq[2])  //  TODO: Need check
                      + 1.5 * square(SM_g) * (q_sq * xlogx(vm_sq[0] / q_sq) - 1. / 3. * vm_sq[0])
                      + 0.75 * (square(SM_g) + square(SM_gp)) * (q_sq * xlogx(vm_sq[1] / q_sq) - 1. / 3. * vm_sq[1]));

    // Goldstone mass
    double mg_sq = chosen_muh_sq + lambda_h * square(h) + 0.5 * lambda_hs * square(s) + (use_Goldstone_resum ? sum : 0.);

    // CP even Higgs thermal temperature masses
    Eigen::MatrixXd MTH2 = Eigen::MatrixXd::Zero(2, 2);
    MTH2(0, 0) = m11_sq + thermal_sq[0];
    MTH2(1, 1) = m22_sq + thermal_sq[1];
    // Mixing between Higgs and singlet
    MTH2(0, 1) = MTH2(1, 0) = m12_sq;
    // Get eigenvalues
    const Eigen::VectorXd mH_sq = MTH2.eigenvalues().real();

    if (!use_covariant_gauge) {
      // This is the ordinary R_\xi gauge

      // Goldstone finite temperature masses
      const double mg0_sq = mg_sq + thermal_sq[0] + 0.25 * xi * (square(SM_g * h) + square(SM_gp * h));
      const double mgpm_sq = mg_sq + thermal_sq[0] + 0.25 * xi * square(SM_g * h);

      // Vector for all scalars, including two mass degenerate charged goldstones
      return {mH_sq(0), mH_sq(1), mg0_sq, mgpm_sq, mgpm_sq};
    } else {
      // This is the covariant gauge. The parameter xi_covariant_internal
      // plays the role of xi. We must have xi = 0 so that there is no
      // explicit xi dependence in the potential

      const double mode1 = (chosen_muh_sq + 0.5 * lambda_hs * square(s) + lambda_h * square(h))
                            * xi_covariant_internal * square(SM::g) * square(h);

      const double mode2 = (chosen_muh_sq + 0.5 * lambda_hs * square(s) + lambda_h * square(h))
                            * xi_covariant_internal * (square(SM::g) + square(SM::gp)) * square(h);

      const double m1p_sq = 0.5 * (mg_sq + real_sqrt(square(mg_sq) - mode1)) + thermal_sq[0];
      const double m1m_sq = 0.5 * (mg_sq - real_sqrt(square(mg_sq) - mode1)) + thermal_sq[0];
      const double m2p_sq = 0.5 * (mg_sq + real_sqrt(square(mg_sq) - mode2)) + thermal_sq[0];
      const double m2m_sq = 0.5 * (mg_sq - real_sqrt(square(mg_sq) - mode2)) + thermal_sq[0];

      return {m1p_sq, m1m_sq, m2p_sq, m2m_sq, mH_sq(0), mH_sq(1)};
    }
  }

  // Physical Higgs bosons and Goldstone bosons
  std::vector<double> get_scalar_dofs() const override {
    if (use_covariant_gauge) {
      return {2., 2., 1., 1., 1., 1.};
    } else {
      return  {1., 1., 1., 1., 1};
    }
  }

  bool iteration_converged = false;
  double get_muh_sq() const { return muh_sq; }
  double get_mus_sq() const { return mus_sq; }
  double get_lambda_h() const { return lambda_h; }
  double get_lambda_s() const { return lambda_s; }
  double get_lambda_hs() const { return lambda_hs; }

 protected:
  /** Whether to use special tadpole constraints in masses entering Coleman-Weinberg potential */
  void set_use_1L_EWSB_in_0L_mass(bool use_1L_EWSB_in_0L_mass_) { use_1L_EWSB_in_0L_mass = use_1L_EWSB_in_0L_mass_; }
  void set_use_Goldstone_resum(bool use_Goldstone_resum_) { use_Goldstone_resum = use_Goldstone_resum_; }
  void set_use_covariant_gauge(bool use_covariant_gauge_) { use_covariant_gauge = use_covariant_gauge_; }

  // For consistency in one-loop potential
  double muh_sq_use_0L_EWSB;
  bool use_1L_EWSB_in_0L_mass{false};
  bool use_Goldstone_resum{true};
  bool use_covariant_gauge{false};

  // hack for covariant gauge
  double xi_covariant_internal{0.};
};

}  // namespace EffectivePotential

#endif
