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
 * See arXiv:2208.01319  [hep-ph] for details
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
      model.iteration_converged = true;
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
    muh_sq_use_0L_EWSB = muh_sq;

    size_t ii = 0;
    while (true) {
      ++ii;
      double lambda_h_prev = lambda_h;
      double lambda_s_prev = lambda_s;
      double muh_sq_prev = muh_sq;
      double mus_sq_prev = mus_sq;
      //        std::cout << "=========="<< std::endl;
      //        std::cout << "muh_sq_prev=" << muh_sq << std::endl;
      //        std::cout << "lambda_h_prev=" << lambda_h << std::endl;
      iterate_one_loop();
      //        std::cout << "muh_sq=" << muh_sq << std::endl;
      //        std::cout << "lambda_h=" << lambda_h << std::endl;
      const double dmuh = std::abs(muh_sq - muh_sq_prev);
      const double dmus = std::abs(mus_sq - mus_sq_prev);
      const double dlambda_h = std::abs(lambda_h - lambda_h_prev);
      const double dlambda_s = std::abs(lambda_s - lambda_s_prev);

      const bool converged = (dmuh < tol) && (dlambda_h < tol) && (dmus < tol) && (dlambda_s < tol);
      if (converged) {
        return true;
      }
      if (ii > 1000) {
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
    const double eT = 0.001;
    return (V1T(phi, T + eT) - V1T(phi, T)) / eT;
  }

  double ddaisy_dT(Eigen::VectorXd phi, double T) const {
    const double eT = 0.001;
    return (daisy(phi, T + eT) - daisy(phi, T)) / eT;
  }

  /**
   * Single iteration of one-loop tadpole solver.
   *
   * Uses numerical derivatives of Coleman-Weinberg potential.
   */
  void iterate_one_loop() {
    Eigen::Vector2d vacuum;
    vacuum << SM_v, 0.;
    double mhh2 = square(SM_mh);
    double mss2 = square(ms);

    if (mhh2 < mss2) {
      std::swap(mhh2, mss2);
    }

    //      std::cout<< "jacobian(0)=" << jacobian(0) << std::endl;
    //      std::cout<< "hessian(0, 0)=" << hessian(0, 0) << std::endl;
    // Apply SM vacuum and Higgs and singlet masses to constraint three parameters
    const auto jacobian = dV1_dx(vacuum);
    muh_sq = -lambda_h * square(SM_v) - jacobian(0) / SM_v;
    //    muh_sq = -0.5 * mhh2 - 1.5 *jacobian(0) / SM_v + 0.5 * hessian(0, 0);
    const auto jacobian_updated = dV1_dx(vacuum);
    const auto hessian = d2V1_dx2(vacuum);
    lambda_h = (mhh2 + jacobian_updated(0) / SM_v - hessian(0, 0)) / (2. * square(SM_v));

    mus_sq = mss2 - 0.5 * lambda_hs * square(SM_v) - hessian(1, 1);

    // Calculate muh_sq using tree level EWSB, for masses in CW potential
    muh_sq_use_0L_EWSB = -lambda_h * square(SM_v);
  }

  double get_v_tree_s() const {
    if (mus_sq < 0)
      return std::sqrt(-mus_sq / lambda_s);
    else
      return 0.;
  }

  /**
   * Tree-level scalar masses including xi-dependence.
   *
   * These masses enter the Coleman-Weinberg potential.
   */
  std::vector<double> get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const override {
    return get_scalar_debye_sq(phi, xi, 0.);
  }

  /**
   * @brief Use the analytic gradient (default) or the finite-difference one.
   */
  void set_use_analytic_gradient(bool use_analytic_gradient_) {
    use_analytic_gradient = use_analytic_gradient_;
  }
  bool get_use_analytic_gradient() const { return use_analytic_gradient; }

  /** @brief Real part of square root */
  double real_sqrt(double x) const {
    return (x > 0.) ? std::sqrt(x) : 0.;
  }

  /**
   * @brief Resummed Goldstone self-energy entering the Goldstone mass.
   */
  double goldstone_resum_sum(Eigen::VectorXd phi) const {
    const double h = phi[0];
    const double s = phi[1];
    const double chosen_muh_sq = use_1L_EWSB_in_0L_mass ? muh_sq : muh_sq_use_0L_EWSB;

    const double m11_sq = chosen_muh_sq + 0.5 * lambda_hs * square(s) + 3. * lambda_h * square(h);
    const double m22_sq = mus_sq + 0.5 * lambda_hs * square(h) + 3. * lambda_s * square(s);

    const auto fm_sq = get_fermion_masses_sq(phi);
    const auto vm_sq = get_vector_masses_sq(phi);
    const double q_sq = square(get_renormalization_scale());

    return 1. / (16. * M_PI * M_PI) * (+3. * lambda_h * (q_sq * xlogx(m11_sq / q_sq) - m11_sq) + 0.5 * lambda_hs * (q_sq * xlogx(m22_sq / q_sq) - m22_sq) - 6. * SM_yt_sq * (q_sq * xlogx(fm_sq[0] / q_sq) - fm_sq[0]) - 6. * SM_yb_sq * (q_sq * xlogx(fm_sq[1] / q_sq) - fm_sq[1]) - 2. * SM_ytau_sq * (q_sq * xlogx(fm_sq[2] / q_sq) - fm_sq[2]) + 1.5 * square(SM_g) * (q_sq * xlogx(vm_sq[0] / q_sq) - 1. / 3. * vm_sq[0]) + 0.75 * (square(SM_g) + square(SM_gp)) * (q_sq * xlogx(vm_sq[1] / q_sq) - 1. / 3. * vm_sq[1]));
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
    const double sum = goldstone_resum_sum(phi);

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

      const double mode1 = mg_sq * xi_covariant_internal * square(SM::g) * square(h);
      const double mode2 = mg_sq * xi_covariant_internal * (square(SM::g) + square(SM::gp)) * square(h);

      const double m1p_sq = 0.5 * (mg_sq + real_sqrt(square(mg_sq) - mode1)) + thermal_sq[0];
      const double m1m_sq = 0.5 * (mg_sq - real_sqrt(square(mg_sq) - mode1)) + thermal_sq[0];
      const double m2p_sq = 0.5 * (mg_sq + real_sqrt(square(mg_sq) - mode2)) + thermal_sq[0];
      const double m2m_sq = 0.5 * (mg_sq - real_sqrt(square(mg_sq) - mode2)) + thermal_sq[0];

      //        std::cout << "xi=" << xi << std::endl;
      //        std::cout << "xi_covariant_internal=" << xi_covariant_internal << std::endl;
      //        std::cout << ".................................." << std::endl;
      //        std::cout << ".......h=" << h << ", s =" << s << std::endl;
      //        std::cout << ".......m1p_sq=" << m1p_sq << ", m1m_sq =" << m1m_sq << std::endl;
      //        std::cout << " mg_sq = " << mg_sq << ", mode1 = " << mode1 << std::endl;
      //        std::cout << "mg_sq=" << mg_sq << std::endl;
      //        std::cout << "real_sqrt(square(mg_sq) - mode1)=" << real_sqrt(square(mg_sq) - mode1) << std::endl;
      //        std::cout << "m1m_sq=" << m1m_sq << std::endl;
      //        std::cout << "real_sqrt(square(mg_sq) - mode1)=" << real_sqrt(square(mg_sq) - mode1) << std::endl;
      //        std::cout << "mode1=" << mode1 << std::endl;
      //        std::cout << "xi_covariant_internal * square(SM::g) * square(h)=" << xi_covariant_internal * square(SM::g) * square(h) << std::endl;
      return {m1p_sq, m1m_sq, m2p_sq, m2m_sq, mH_sq(0), mH_sq(1)};
    }
  }

  // Physical Higgs bosons and Goldstone bosons
  std::vector<double> get_scalar_dofs() const override {
    if (use_covariant_gauge) {
      return {2., 2., 1., 1., 1., 1.};
    } else {
      return {1., 1., 1., 1., 1};
    }
  }

  /**
   * Full SM fermion content: t, b, tau, c, s, u, d, mu, e, nu_e, nu_mu, nu_tau.
   *
   * Light quarks and charged leptons couple to the Higgs via Yukawa interactions
   * and develop field-dependent masses m_f^2 = y_f^2 h^2 / 2.  Neutrinos are
   * treated as massless (SM Dirac limit, left-handed only).  Including these
   * species ensures the one-loop thermal potential captures the full g* T^4
   * radiation background without the need for a separate background_dof term.
   */
  std::vector<double> get_fermion_masses_sq(Eigen::VectorXd phi) const override {
    const double h_sq = square(phi[0]);
    return {
      0.5 * SM_yt_sq  * h_sq,  // top
      0.5 * SM_yb_sq  * h_sq,  // bottom
      0.5 * SM_ytau_sq * h_sq, // tau
      0.5 * SM_yc_sq  * h_sq,  // charm
      0.5 * SM_ys_sq  * h_sq,  // strange
      0.5 * SM_yu_sq  * h_sq,  // up
      0.5 * SM_yd_sq  * h_sq,  // down
      0.5 * SM_ymu_sq * h_sq,  // muon
      0.5 * SM_ye_sq  * h_sq,  // electron
      0.,                       // nu_e  (LH Weyl: nu_L + anti-nu_R)
      0.,                       // nu_mu
      0.,                       // nu_tau
    };
  }

  /** Degrees of freedom matching get_fermion_masses_sq order above. */
  std::vector<double> get_fermion_dofs() const override {
    return {
      12.,  // top:     N_c=3, 2 spins, particle+antiparticle
      12.,  // bottom
       4.,  // tau:     2 spins × 2
      12.,  // charm
      12.,  // strange
      12.,  // up
      12.,  // down
       4.,  // muon
       4.,  // electron
       2.,  // nu_e:    LH only × 2 (nu_L + anti-nu_R)
       2.,  // nu_mu
       2.,  // nu_tau
    };
  }

  /**
   * Vector sector: EW bosons (W, Z, gamma) plus 8 gluons.
   *
   * Gluons are SU(3) gauge bosons with no coupling to the Higgs field; their
   * T=0 (transverse) mass is zero and they pick up a field-independent Debye
   * mass m_D^2 = g_s^2 T^2 (N_c/3 + N_f/6) = 2 g_s^2 T^2 (N_c=3, N_f=6) for
   * the longitudinal (screened) modes.  The 16 transverse modes are massless at
   * all field values and contribute the full bosonic radiation g T^4 term.
   *
   * Layout: {MW_L(2), MZ_L(1), Mphoton_L(1), MW_T(4), MZ_T(2), Mphoton_T(2),
   *          Mgluon_L(8), Mgluon_T(16)}
   */
  std::vector<double> get_vector_dofs() const override {
    return {2., 1., 1., 4., 2., 2., 8., 16.};
  }

  std::vector<double> get_vector_debye_sq(Eigen::VectorXd phi, double T) const override {
    const double h_sq  = square(phi[0]);
    const double T_sq  = square(T);

    // ---- EW sector (identical to xSM_base) ----
    const double MW_T_sq      = 0.25 * square(SM_g) * h_sq;
    const double MZ_T_sq      = 0.25 * (square(SM_g) + square(SM_gp)) * h_sq;
    const double Mphoton_T_sq = 0.;

    const double MW_L_sq = 0.25 * square(SM_g) * h_sq + 11. / 6. * square(SM_g) * T_sq;
    const double a_L = (square(SM_g) + square(SM_gp)) * (3. * h_sq + 22. * T_sq);
    const double b_L = std::sqrt(9.  * square(square(SM_g) + square(SM_gp)) * square(h_sq) +
                                 132. * square(square(SM_g) - square(SM_gp)) * h_sq  * T_sq +
                                 484. * square(square(SM_g) - square(SM_gp)) * pow_4(T));
    const double MZ_L_sq      = (a_L + b_L) / 24.;
    const double Mphoton_L_sq = (a_L - b_L) / 24.;

    // ---- QCD gluon sector ----
    // Longitudinal Debye mass: m_D^2 = g_s^2 T^2 (N_c/3 + N_f/6) with N_c=3, N_f=6
    const double Mgluon_L_sq = 2. * SM_gs_sq * T_sq;
    // Transverse gluons are massless at every field point
    const double Mgluon_T_sq = 0.;

    return {MW_L_sq, MZ_L_sq, Mphoton_L_sq, MW_T_sq, MZ_T_sq, Mphoton_T_sq,
            Mgluon_L_sq, Mgluon_T_sq};
  }

  /** log|m^2/Q^2|, or 0 where the mass is below the cutoff xlogx uses. */
  double log_mass_ratio(double m_sq, double q_sq) const {
    if (std::abs(m_sq) <= std::numeric_limits<double>::min()) {
      return 0.;
    }
    return std::log(std::abs(m_sq / q_sq));
  }

  /** Gradient of the tree-level potential. */
  Eigen::VectorXd grad_V0(Eigen::VectorXd phi) const {
    const double h = phi[0];
    const double s = phi[1];
    Eigen::VectorXd d(2);
    d << muh_sq * h + lambda_h * h * square(h) + 0.5 * lambda_hs * h * square(s),
         mus_sq * s + lambda_s * s * square(s) + 0.5 * lambda_hs * square(h) * s;
    return d;
  }

  /** Gradients matching get_fermion_masses_sq: m^2 = y^2 h^2 / 2, so d/dh = y^2 h. */
  std::vector<Eigen::VectorXd> grad_fermion_masses_sq(Eigen::VectorXd phi) const {
    const double h = phi[0];
    const double y_sq[12] = {SM_yt_sq, SM_yb_sq, SM_ytau_sq, SM_yc_sq, SM_ys_sq,
                             SM_yu_sq, SM_yd_sq, SM_ymu_sq, SM_ye_sq, 0., 0., 0.};
    std::vector<Eigen::VectorXd> d;
    d.reserve(12);
    for (int i = 0; i < 12; ++i) {
      Eigen::VectorXd v(2);
      v << y_sq[i] * h, 0.;
      d.push_back(v);
    }
    return d;
  }

  /** Gradients matching get_vector_debye_sq. Only h enters. */
  std::vector<Eigen::VectorXd> grad_vector_debye_sq(Eigen::VectorXd phi, double T) const {
    const double h = phi[0];
    const double T_sq = square(T);
    const double g_sq = square(SM_g);
    const double gp_sq = square(SM_gp);
    const double sum_sq = g_sq + gp_sq;
    const double dif_sq = g_sq - gp_sq;

    const double d_MW_T = 0.5 * g_sq * h;
    const double d_MZ_T = 0.5 * sum_sq * h;
    const double d_MW_L = 0.5 * g_sq * h;
    const double d_a_L = 6. * sum_sq * h;

    // b_L = sqrt(inside); inside is quartic in h.
    const double inside = 9. * square(sum_sq) * pow_4(h) +
                          132. * square(dif_sq) * square(h) * T_sq +
                          484. * square(dif_sq) * pow_4(T);
    const double b_L = std::sqrt(inside);
    const double d_inside = 36. * square(sum_sq) * h * square(h) +
                            264. * square(dif_sq) * h * T_sq;
    // b_L vanishes only at h = T = 0, where the photon and Z longitudinal modes
    // become degenerate and the square root is not differentiable.
    const double d_b_L = (b_L > 0.) ? d_inside / (2. * b_L) : 0.;

    const auto along_h = [](double dh) {
      Eigen::VectorXd v(2);
      v << dh, 0.;
      return v;
    };

    return {along_h(d_MW_L),
            along_h((d_a_L + d_b_L) / 24.),
            along_h((d_a_L - d_b_L) / 24.),
            along_h(d_MW_T),
            along_h(d_MZ_T),
            along_h(0.),
            along_h(0.),   // gluon Debye mass is field independent
            along_h(0.)};
  }

  /** Gradients matching get_ghost_masses_sq: xi times the T=0 vector masses. */
  std::vector<Eigen::VectorXd> grad_ghost_masses_sq(Eigen::VectorXd phi, double xi_) const {
    const auto dv = grad_vector_debye_sq(phi, 0.);
    return {xi_ * dv[0], xi_ * dv[1], xi_ * dv[2]};
  }

  /**
   * Scalar masses and their gradients together, in the ordinary R_xi gauge.
   */
  void scalar_debye_sq_with_grad(Eigen::VectorXd phi, double xi_, double T,
                                 std::vector<double> &masses_sq,
                                 std::vector<Eigen::VectorXd> &grads) const {
    const double h = phi[0];
    const double s = phi[1];
    const double chosen_muh_sq = use_1L_EWSB_in_0L_mass ? muh_sq : muh_sq_use_0L_EWSB;
    const auto thermal_sq = get_scalar_thermal_sq(T);  // field independent

    const double m11_sq = chosen_muh_sq + 0.5 * lambda_hs * square(s) + 3. * lambda_h * square(h);
    const double m22_sq = mus_sq + 0.5 * lambda_hs * square(h) + 3. * lambda_s * square(s);
    const double m12_sq = lambda_hs * s * h;

    Eigen::VectorXd d_m11(2), d_m22(2), d_m12(2);
    d_m11 << 6. * lambda_h * h, lambda_hs * s;
    d_m22 << lambda_hs * h, 6. * lambda_s * s;
    d_m12 << lambda_hs * s, lambda_hs * h;

    // Goldstone resummation. With L(m) = m log|m/Q^2| - m the derivative is
    // L'(m) = log|m/Q^2|, and for the vector terms, which carry -m/3 instead of
    // -m, it is log|m/Q^2| + 2/3. Unlike in V1 the log is not multiplied by m^2,
    // so it diverges as a mass goes to zero; for the fermion and vector terms
    // d(m^2)/dh vanishes proportionally to h and tames it.
    Eigen::VectorXd d_sum = Eigen::VectorXd::Zero(2);
    if (use_Goldstone_resum) {
      const auto fm_sq = get_fermion_masses_sq(phi);
      const auto vm_sq = get_vector_masses_sq(phi);
      const auto d_fm = grad_fermion_masses_sq(phi);
      const auto d_vm = grad_vector_debye_sq(phi, 0.);
      const double q_sq = square(get_renormalization_scale());

      d_sum += (3. * lambda_h * log_mass_ratio(m11_sq, q_sq)) * d_m11;
      d_sum += (0.5 * lambda_hs * log_mass_ratio(m22_sq, q_sq)) * d_m22;
      d_sum -= (6. * SM_yt_sq * log_mass_ratio(fm_sq[0], q_sq)) * d_fm[0];
      d_sum -= (6. * SM_yb_sq * log_mass_ratio(fm_sq[1], q_sq)) * d_fm[1];
      d_sum -= (2. * SM_ytau_sq * log_mass_ratio(fm_sq[2], q_sq)) * d_fm[2];
      d_sum += (1.5 * square(SM_g) * (log_mass_ratio(vm_sq[0], q_sq) + 2. / 3.)) * d_vm[0];
      d_sum += (0.75 * (square(SM_g) + square(SM_gp)) *
                (log_mass_ratio(vm_sq[1], q_sq) + 2. / 3.)) * d_vm[1];
      d_sum /= (16. * M_PI * M_PI);
    }

    const double mg_sq = chosen_muh_sq + lambda_h * square(h) + 0.5 * lambda_hs * square(s) +
                         (use_Goldstone_resum ? goldstone_resum_sum(phi) : 0.);
    Eigen::VectorXd d_mg(2);
    d_mg << 2. * lambda_h * h, lambda_hs * s;
    d_mg += d_sum;

    // CP-even sector, closed form.
    const double A = m11_sq + thermal_sq[0];
    const double B = m22_sq + thermal_sq[1];
    const double C = m12_sq;
    const double D = std::sqrt(square(A - B) + 4. * square(C));

    const Eigen::VectorXd d_half_sum = 0.5 * (d_m11 + d_m22);
    Eigen::VectorXd d_half_D = Eigen::VectorXd::Zero(2);
    if (D > 0.) {
      d_half_D = ((A - B) * (d_m11 - d_m22) + 4. * C * d_m12) / (2. * D);
    }

    const double mHp_sq = 0.5 * (A + B) + 0.5 * D;
    const double mHm_sq = 0.5 * (A + B) - 0.5 * D;

    const double mg0_sq = mg_sq + thermal_sq[0] +
                          0.25 * xi_ * (square(SM_g * h) + square(SM_gp * h));
    const double mgpm_sq = mg_sq + thermal_sq[0] + 0.25 * xi_ * square(SM_g * h);

    Eigen::VectorXd d_mg0 = d_mg;
    d_mg0(0) += 0.5 * xi_ * (square(SM_g) + square(SM_gp)) * h;
    Eigen::VectorXd d_mgpm = d_mg;
    d_mgpm(0) += 0.5 * xi_ * square(SM_g) * h;

    masses_sq = {mHp_sq, mHm_sq, mg0_sq, mgpm_sq, mgpm_sq};
    grads = {d_half_sum + d_half_D, d_half_sum - d_half_D, d_mg0, d_mgpm, d_mgpm};
  }

  /**
   * @brief Analytic gradient of the full effective potential.
   */
  Eigen::VectorXd dV_dx(Eigen::VectorXd phi, double T) const override {

    if (use_covariant_gauge || !use_analytic_gradient) {
      return EffectivePotential::Potential::dV_dx(phi, T);
    }

    const double xi_ = get_xi();

    const auto fermion_masses_sq = get_fermion_masses_sq(phi);
    const auto d_fermion = grad_fermion_masses_sq(phi);
    const auto fermion_dofs = get_fermion_dofs();

    const auto ghost_masses_sq = (xi_ != 0.) ? get_ghost_masses_sq(phi, xi_) : std::vector<double>{};
    const auto d_ghost = (xi_ != 0.) ? grad_ghost_masses_sq(phi, xi_) : std::vector<Eigen::VectorXd>{};
    const auto ghost_dofs = (xi_ != 0.) ? get_ghost_dofs() : std::vector<double>{};

    const auto scalar_dofs = get_scalar_dofs();
    const auto vector_dofs = get_vector_dofs();

    // Scalar masses come back paired with their own derivatives; see
    // scalar_debye_sq_with_grad for why they are not read from
    // get_scalar_debye_sq.
    std::vector<double> scalar_sq, scalar_debye_sq;
    std::vector<Eigen::VectorXd> d_scalar, d_scalar_debye;

    Eigen::VectorXd grad = grad_V0(phi);

    // Fermions and ghosts always use the ordinary masses, in every branch.
    const auto add_fermions_and_ghosts = [&](bool thermal) {
      grad += dV1_term(fermion_masses_sq, fermion_dofs, d_fermion, -1., 1.5);
      if (xi_ != 0.) {
        grad += dV1_term(ghost_masses_sq, ghost_dofs, d_ghost, -1., 1.5);
      }
      if (thermal) {
        grad += dV1T_term(fermion_masses_sq, fermion_dofs, d_fermion, +1., T, true);
        if (xi_ != 0.) {
          grad += dV1T_term(ghost_masses_sq, ghost_dofs, d_ghost, -1., T, false);
        }
      }
    };

    if (T > 0) {
      switch (get_daisy_method()) {
      case DaisyMethod::None: {
        scalar_debye_sq_with_grad(phi, xi_, 0., scalar_sq, d_scalar);
        const auto vector_sq = get_vector_masses_sq(phi);
        const auto d_vector = grad_vector_debye_sq(phi, 0.);
        grad += dV1_term(scalar_sq, scalar_dofs, d_scalar, +1., 1.5);
        grad += dV1_term(vector_sq, vector_dofs, d_vector, +1., 5. / 6.);
        grad += dV1T_term(scalar_sq, scalar_dofs, d_scalar, +1., T, false);
        grad += dV1T_term(vector_sq, vector_dofs, d_vector, +1., T, false);
        add_fermions_and_ghosts(true);
        break;
      }
      case DaisyMethod::ArnoldEspinosa: {
        scalar_debye_sq_with_grad(phi, xi_, 0., scalar_sq, d_scalar);
        scalar_debye_sq_with_grad(phi, xi_, T, scalar_debye_sq, d_scalar_debye);
        const auto vector_sq = get_vector_masses_sq(phi);
        const auto d_vector = grad_vector_debye_sq(phi, 0.);
        const auto vector_debye_sq = get_vector_debye_sq(phi, T);
        const auto d_vector_debye = grad_vector_debye_sq(phi, T);
        grad += ddaisy_term(scalar_sq, scalar_debye_sq, scalar_dofs, d_scalar, d_scalar_debye, T);
        grad += ddaisy_term(vector_sq, vector_debye_sq, vector_dofs, d_vector, d_vector_debye, T);
        grad += dV1_term(scalar_sq, scalar_dofs, d_scalar, +1., 1.5);
        grad += dV1_term(vector_sq, vector_dofs, d_vector, +1., 5. / 6.);
        grad += dV1T_term(scalar_sq, scalar_dofs, d_scalar, +1., T, false);
        grad += dV1T_term(vector_sq, vector_dofs, d_vector, +1., T, false);
        add_fermions_and_ghosts(true);
        break;
      }
      case DaisyMethod::Parwani: {
        scalar_debye_sq_with_grad(phi, xi_, T, scalar_debye_sq, d_scalar_debye);
        const auto vector_debye_sq = get_vector_debye_sq(phi, T);
        const auto d_vector_debye = grad_vector_debye_sq(phi, T);
        grad += dV1_term(scalar_debye_sq, scalar_dofs, d_scalar_debye, +1., 1.5);
        grad += dV1_term(vector_debye_sq, vector_dofs, d_vector_debye, +1., 5. / 6.);
        grad += dV1T_term(scalar_debye_sq, scalar_dofs, d_scalar_debye, +1., T, false);
        grad += dV1T_term(vector_debye_sq, vector_dofs, d_vector_debye, +1., T, false);
        add_fermions_and_ghosts(true);
        break;
      }
      default:
        throw std::runtime_error("unknown daisy method");
      }
    } else {
      scalar_debye_sq_with_grad(phi, xi_, 0., scalar_sq, d_scalar);
      const auto vector_sq = get_vector_masses_sq(phi);
      const auto d_vector = grad_vector_debye_sq(phi, 0.);
      grad += dV1_term(scalar_sq, scalar_dofs, d_scalar, +1., 1.5);
      grad += dV1_term(vector_sq, vector_dofs, d_vector, +1., 5. / 6.);
      add_fermions_and_ghosts(false);
    }

    return grad;
  }

  // std::vector<double> get_4d_params() const {
  //   return {g_sq, gp_sq, v, m_s, muh_sq, lambda_h, lambda_hs, lambda_s, yt_sq};
  // }

  std::map<std::string, double> get_4d_parameter_map(double T) const {

    std::map<std::string, double> param_map;

    param_map["v0"] = 247.4544243292407;
    param_map["Mt"] = 173.0;
    param_map["MW"] = 80.379;
    param_map["MZ"] = 91.1876;
    param_map["mh1"] = 125.25;
    param_map["mh2"] = ms;
    param_map["lHS"] = lambda_hs;
    param_map["lSS"] = lambda_s;
    param_map["RGScale"] = 10.;
    param_map["g3"] = 1.2279920495357861;

    return param_map;
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
  bool use_analytic_gradient{true};

  // hack for covariant gauge
  double xi_covariant_internal{0.};

  // Additional SM Yukawa couplings (light quarks and leptons) and QCD coupling
  double SM_yc_sq  = SM::yc_sq;
  double SM_ys_sq  = SM::ys_sq;
  double SM_yu_sq  = SM::yu_sq;
  double SM_yd_sq  = SM::yd_sq;
  double SM_ymu_sq = SM::ymu_sq;
  double SM_ye_sq  = SM::ye_sq;
  double SM_gs_sq  = SM::gs_sq;
};

} // namespace EffectivePotential

#endif
