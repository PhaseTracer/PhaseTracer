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

#include <algorithm>
#include <utility>
#include <vector>
#include <Eigen/Eigenvalues>
#include <iostream>

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
            double lambda_s_,
            double ms_):
    lambda_hs(lambda_hs_), lambda_s(lambda_s_), ms(ms_) {}

  /**
   * @brief Make an xSM model using tree or one-loop tadpole constraints
   */
  static xSM_MSbar from_tadpoles(double lambda_hs, double lambda_s, double ms, 
                                 double Q, double xi, bool tree_level, bool tree_ewsb = false) {
    xSM_MSbar model(lambda_hs, lambda_s, ms);
    model.set_tree_ewsb(tree_ewsb);
    model.set_renormalization_scale(Q);
    model.set_xi(xi);
    if (tree_level) {
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
    double mhh2 = square(SM::mh);
    double mss2 = square(ms);
    
    // TODO: will this affect calculation of lambda_h?
    if (mss2 > mhh2) {
      std::swap(mhh2, mss2);
    }

    // Apply SM vacuum and Higgs and singlet masses to constraint three parameters
    lambda_h = mhh2 / (2. * square(SM::v));
    muh_sq = -lambda_h * square(SM::v);
    mus_sq = mss2 - 0.5 * lambda_hs * square(SM::v);
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
//      std::cout << "converged = " << converged << std::endl;
      if (converged) {
        return true;
      }
      if (ii > 1000){
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

  /**
   * Single iteration of one-loop tadpole solver.
   *
   * Uses numerical derivatives of Coleman-Weinberg potential.
   */
  void iterate_one_loop() {
    Eigen::Vector2d vacuum;
    vacuum << SM::v, 0.;
    const auto jacobian = dV1_dx(vacuum);
    const auto hessian = d2V1_dx2(vacuum);

    const double mh_sq = square(SM::mh);
    const double ms_sq = square(ms);
    
//    const double a = mh_sq + ms_sq;
//    const double discriminant = square(mh_sq - ms_sq) - 4. * square(hessian(1, 0));
//    if (discriminant < 0) {
//      throw std::runtime_error("Could not solve 1l tadpoles");
//    }
//    const double b = sqrt(discriminant);
//    double mhh2 = 0.5 * (a + b);
//    double mss2 = 0.5 * (a - b);
    
    double mhh2 = mh_sq;
    double mss2 = ms_sq;
    
    // TODO: will this affect calculation of lambda_h?
    if (mhh2 < mss2) {
      std::swap(mhh2, mss2);
    }

    // Apply SM vacuum and Higgs and singlet masses to constraint three parameters

    lambda_h = (mhh2 + jacobian(0) / SM::v - hessian(0, 0)) / (2. * square(SM::v));
    muh_sq = -0.5 * mhh2 - 1.5 *jacobian(0) / SM::v + 0.5 * hessian(0, 0);
    mus_sq = mss2 - 0.5 * lambda_hs * square(SM::v) - hessian(1, 1);

    muh_sq_tree_ewsb = - lambda_h * square(SM::v);
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

  /**
   * Thermal scalar masses of form c * T^2 etc for high-temperature expansion of potential
   */
  std::vector<double> get_scalar_thermal_sq(double T) const override {
    const double c_h = (9. * square(SM::g) +
                        3. * square(SM::gp) +
                        2. * (6. * SM::yt_sq + 12. * lambda_h + lambda_hs)) / 48.;
    const double c_s = (2. * lambda_hs + 3. * lambda_s) / 12.;
    return {c_h * square(T), c_s * square(T)};
  }

  /**
   * Tree-level scalar masses including xi-dependence.
   *
   * These masses enter the Coleman-Weinberg potential.
   */
  std::vector<double> get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const override {
    return get_scalar_debye_sq(phi, xi, 0.);
  }

  std::vector<double> get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const override {
    const double h = phi[0];
    const double s = phi[1];
    const auto thermal_sq = get_scalar_thermal_sq(T);
    
    // TODO: add thermal_sq here and use get_vector_debye_sq?
    const double mhh2 = (tree_ewsb ? muh_sq_tree_ewsb : muh_sq) + 3. * lambda_h * square(h) + 0.5 * lambda_hs * square(s);
    const double mgg2 = (tree_ewsb ? muh_sq_tree_ewsb : muh_sq) + lambda_h * square(h) + 0.5 * lambda_hs * square(s);
    const double mss2 = mus_sq + 3. * lambda_s * square(s) + 0.5 * lambda_hs * square(h);
    
    
    // resummed Goldstone contributions
    const auto fm2 = get_fermion_masses_sq(phi);
    const auto vm2 = get_vector_masses_sq(phi);
    const double Qsq = square(get_renormalization_scale());
    const double sum = 1. / (16. * M_PI * M_PI) * (
                       3.  * lambda_h * (Qsq*xlogx(mhh2/Qsq) - mhh2)
                      +0.5 * lambda_hs * (Qsq*xlogx(mss2/Qsq) - mss2)
                      -6.  * SM::yt_sq * (Qsq*xlogx(fm2[0]/Qsq) - fm2[0])
                      +1.5 * square(SM::g) * (Qsq*xlogx(vm2[0]/Qsq) - 1./3.*vm2[0])
                      +0.75* (square(SM::g)+square(SM::gp)) * (Qsq*xlogx(vm2[1]/Qsq) - 1./3.*vm2[1])
                      );
                      
    // Higgs and Goldstone diagonals
    // Eigen::MatrixXd M2 = Eigen::MatrixXd::Zero(5, 5);
    // M2(0, 0) = mgg2 + thermal_sq[0] + (tree_ewsb ? sum : 0);
    // M2(1, 1) = M2(0, 0);
    // M2(2, 2) = mhh2 + thermal_sq[0];
    // M2(3, 3) = M2(0, 0);

    // // Singlet mass diagonal
    // M2(4, 4) = mss2 + thermal_sq[1];

    // // Mixing between Higgs and singlet
    // M2(2, 4) = M2(4, 2) = lambda_hs * h * s;

    // // xi-dependence
    // M2(0, 0) += 0.25 * xi * square(SM::g * h);
    // M2(1, 1) += 0.25 * xi * square(SM::g * h);
    // M2(3, 3) += 0.25 * xi * (square(SM::g * h) + square(SM::gp * h));
    // Goldstone finite temperature masses
    double mTG02 =   mgg2 + thermal_sq[0] + (tree_ewsb ? sum : 0);
    double mTGpm2 = mTG02; // 2 degrees of freedom or two degenerate copies
    // xi-dependence
    mTG02 += 0.25 * xi * (square(SM::g * h) + square(SM::gp * h));
    mTGpm2 += 0.25 * xi * square(SM::g * h);
    // CP even Higgs thermal temperature masses
    Eigen::MatrixXd MTH2 = Eigen::MatrixXd::Zero(2, 2); 
    MTH2(0,0) = mhh2 + thermal_sq[0];
    MTH2(1,1) = mss2 + thermal_sq[1];
    // get eigenvalues
    const Eigen::VectorXd mH_sq = MTH2.eigenvalues().real();
    // vector for all scalars, including two mass degenerate charged goldstones
    std::vector<double> m_sq_vector{MTH2(0), MTH2(1), mTG02, mTGpm2, mTGpm2};
    // mass order
    std::sort(m_sq_vector.begin(), m_sq_vector.end());
    return m_sq_vector;
  }

  // Physical Higgs bosons and Goldstone bosons
  std::vector<double> get_scalar_dofs() const override {
    return {1., 1., 1., 1., 1};
  }

  // W, Z, photon
  std::vector<double> get_vector_masses_sq(Eigen::VectorXd phi) const override {
    return get_vector_debye_sq(phi, 0.);
  }

  // W, Z, photon
  std::vector<double> get_vector_debye_sq(Eigen::VectorXd phi, double T) const override {
    const double h_sq = square(phi[0]);
    const double T_sq = square(T);
    const double MW_sq = 0.25 * square(SM::g) * h_sq + 11. / 6. * square(SM::g) * T_sq;

		const double a = (square(SM::g) + square(SM::gp)) * (3. * h_sq + 22. * T_sq);
		const double b = std::sqrt(9. * square(square(SM::g) + square(SM::gp)) * square(h_sq) 
                     + 132. * square(square(SM::g) - square(SM::gp)) * h_sq * T_sq
                     + 484. * square(square(SM::g) - square(SM::gp)) * pow_4(T));

		const double MZ_sq = (a + b) / 24.;
		const double Mphoton_sq = (a - b) / 24.;

    return {MW_sq, MZ_sq, Mphoton_sq};
  }

  // W, Z, photon
  std::vector<double> get_vector_dofs() const override {
    return {6., 3., 2.}; // TODO: check photon dof
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

  /** Whether to use special tadpole constraints in masses entering Coleman-Weinberg potential */
  void set_tree_ewsb(bool tree_ewsb_) { tree_ewsb = tree_ewsb_; }

  bool iteration_converged = false;
  
 protected:
  double ms;
  // Lagrangian parameters
  double lambda_hs;
  double muh_sq;
  double lambda_h;
  double mus_sq;
  double lambda_s;
  // For consistency in one-loop potential
  double muh_sq_tree_ewsb;
  bool tree_ewsb{false};

};

}  // namespace EffectivePotential

#endif
