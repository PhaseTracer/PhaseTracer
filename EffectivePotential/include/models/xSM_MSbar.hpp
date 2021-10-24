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

namespace EffectivePotential {

class xSM_MSbar : public xSM_base {
 public:
  /**
   * @brief Make an xSM model from Lagrangian parameters
   */
  xSM_MSbar(double lambda_hs_,
            double lambda_s_,
            double ms_){
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
                                 std::vector<double> SM_parameters={}) {
    xSM_MSbar model(lambda_hs, lambda_s, ms);
    if (SM_parameters.size() == 7)
      model.set_SM_parameters(SM_parameters);
    model.set_use_1L_EWSB_in_0L_mass(use_1L_EWSB_in_0L_mass);
    model.set_use_Goldstone_resum(use_Goldstone_resum);
    model.set_use_covariant_gauge(use_covariant_gauge);
    model.set_renormalization_scale(Q);
    model.set_xi(xi);
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
    muh_sq_use_0L_EWSB = muh_sq; // TODO: this can be deleted?
    
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
  double real_sqrt(double x) const{
    return (x > 0.) ? std::sqrt(x) : 0.;
  }

  std::vector<double> get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const override {
    const double h = phi[0];
    const double s = phi[1];
    if (not use_covariant_gauge) {
      const auto thermal_sq = get_scalar_thermal_sq(T);
    
      const double mhh2 = (use_1L_EWSB_in_0L_mass ? muh_sq : muh_sq_use_0L_EWSB) + 3. * lambda_h * square(h) + 0.5 * lambda_hs * square(s);
      const double mgg2 = (use_1L_EWSB_in_0L_mass ? muh_sq : muh_sq_use_0L_EWSB) + lambda_h * square(h) + 0.5 * lambda_hs * square(s);
      const double mss2 = mus_sq + 3. * lambda_s * square(s) + 0.5 * lambda_hs * square(h);
      
      // resummed Goldstone contributions
      const auto fm2 = get_fermion_masses_sq(phi);
      const auto vm2 = get_vector_masses_sq(phi);
      const double Qsq = square(get_renormalization_scale());
      const double sum = 1. / (16. * M_PI * M_PI) * (
                         3.  * lambda_h * (Qsq*xlogx(mhh2/Qsq) - mhh2)
                        +0.5 * lambda_hs * (Qsq*xlogx(mss2/Qsq) - mss2)
                        -6.  * SM_yt_sq * (Qsq*xlogx(fm2[0]/Qsq) - fm2[0])
                        -6.  * SM_yb_sq * (Qsq*xlogx(fm2[1]/Qsq) - fm2[1]) // TODO: Need check
                        -2.  * SM_ytau_sq * (Qsq*xlogx(fm2[2]/Qsq) - fm2[2]) // TODO: Need check
                        +1.5 * square(SM_g) * (Qsq*xlogx(vm2[0]/Qsq) - 1./3.*vm2[0])
                        +0.75* (square(SM_g)+square(SM_gp)) * (Qsq*xlogx(vm2[1]/Qsq) - 1./3.*vm2[1])
                        );

      // Goldstone finite temperature masses
      double mTG02 =   mgg2 + thermal_sq[0] + (use_Goldstone_resum ? sum : 0);
      double mTGpm2 = mTG02; // 2 degrees of freedom or two degenerate copies
      // xi-dependence
      mTG02 += 0.25 * xi * (square(SM_g * h) + square(SM_gp * h));
      mTGpm2 += 0.25 * xi * square(SM_g * h);
      // CP even Higgs thermal temperature masses
      Eigen::MatrixXd MTH2 = Eigen::MatrixXd::Zero(2, 2); 
      MTH2(0,0) = mhh2 + thermal_sq[0];
      MTH2(1,1) = mss2 + thermal_sq[1];
      // Mixing between Higgs and singlet
      MTH2(0, 1) = MTH2(1, 0) = lambda_hs * h * s;
      // get eigenvalues
      const Eigen::VectorXd mH_sq = MTH2.eigenvalues().real();
      // vector for all scalars, including two mass degenerate charged goldstones
      return {mH_sq(0), mH_sq(1), mTG02, mTGpm2, mTGpm2};
    } else {
      const double mgb_sq = muh_sq + lambda_h * square(phi[0]) + 0.5 * lambda_hs * square(phi[1]);

      const double mode1 = (muh_sq + 0.5 * lambda_hs * square(phi[1]) + lambda_h * square(phi[0]))
                            * xi * square(SM::g) * square(phi[0]);

      const double mode2 = (muh_sq + 0.5 * lambda_hs * square(phi[1]) + lambda_h * square(phi[0]))
                            * xi * (square(SM::g) + square(SM::gp)) * square(phi[0]);

      const double pi_11 = (0.5 * lambda_h + 0.0625 * square(SM::gp) 
                            + 3. * 0.0625 * square(SM::g) + 0.25 * SM::yt_sq) * square(T);
      const double pi_22 = 0.25 * lambda_s * square(T); 

      const double m1p_sq = 0.5 * (mgb_sq + real_sqrt(square(mgb_sq) - mode1)) + pi_11; 
      const double m1m_sq = 0.5 * (mgb_sq - real_sqrt(square(mgb_sq) - mode1)) + pi_11; 
      const double m2p_sq = 0.5 * (mgb_sq + real_sqrt(square(mgb_sq) - mode2)) + pi_11; 
      const double m2m_sq = 0.5 * (mgb_sq - real_sqrt(square(mgb_sq) - mode2)) + pi_11; 

      const double m11_sq = muh_sq + 0.5 * lambda_hs * square(phi[1]) + 3. * lambda_h * square(phi[0]) + pi_11;
      const double m22_sq = mus_sq + 0.5 * lambda_hs * square(phi[0]) + 3. * lambda_s * square(phi[1]) + pi_22;
      const double m12_sq = lambda_hs * phi[1] * phi[0];

      const double trace = m11_sq + m22_sq;
      const double det = m11_sq * m22_sq - square(m12_sq);

      const double mh_sq = 0.5 * (trace - real_sqrt(square(trace) -4. * det));
      const double ms_sq = 0.5 * (trace + real_sqrt(square(trace) -4. * det));

      return {m1p_sq, m1m_sq, m2p_sq, m2m_sq, mh_sq, ms_sq};
    }
  }
  
  // Physical Higgs bosons and Goldstone bosons
  std::vector<double> get_scalar_dofs() const override {
    if (not use_covariant_gauge) {
      return {1., 1., 1., 1., 1};
    } else {
      return {2., 2., 1., 1., 1., 1.};
    }
    
  }

  bool iteration_converged = false;  
  double get_muh_sq() {return muh_sq;}
  double get_mus_sq() {return mus_sq;}
  double get_lambda_h() {return lambda_h;}
  double get_lambda_s() {return lambda_s;}
  double get_lambda_hs() {return lambda_hs;}
    
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

};

}  // namespace EffectivePotential

#endif
