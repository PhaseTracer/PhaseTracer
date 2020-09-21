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
   Z2 symmetric real scalar singlet extension of the Standard Model
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "xSM_MSbar_solver.hpp"
#include "one_loop_potential.hpp"
#include "pow.hpp"
#include "SM_parameters.hpp"


namespace EffectivePotential {

class xSM_MSbar : public OneLoopPotential {
 public:
  xSM_MSbar(double m_s_,
            double lambda_hs_,
            double muh_sq_,
            double lambda_h_,
            double mus_sq_,
            double lambda_s_,
            double muh_sq_tree_EWSB_,
            double mus_sq_tree_EWSB_):
    m_s(m_s_), lambda_hs(lambda_hs_), muh_sq(muh_sq_),
    lambda_h(lambda_h_), mus_sq(mus_sq_), lambda_s(lambda_s_),
    muh_sq_tree_EWSB(muh_sq_tree_EWSB_), mus_sq_tree_EWSB(mus_sq_tree_EWSB_) {}

  double V0(Eigen::VectorXd phi) const override {
    return 0.5 * muh_sq * square(phi[0]) +
           0.25 * lambda_h * pow_4(phi[0]) +
           0.25 * lambda_hs * square(phi[0]) * square(phi[1]) +
           0.5 * mus_sq * square(phi[1]) +
           0.25 * lambda_s * pow_4(phi[1]);
  }

  // Higgs
  std::vector<double> get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const override {
    const double h = phi[0];
    const double s = phi[1];
    // M is the mass matrix. h is higgs direction, s is singlet
    const double Mhh = muh_sq_tree_EWSB + 3. * lambda_h * square(h) + 0.5 * lambda_hs * square(s);
    const double Mss = mus_sq_tree_EWSB + 3. * lambda_s * square(s) + 0.5 * lambda_hs * square(h);
    const double Mhs = lambda_hs * h * s;
    // diagonalization
    const double A = 0.5 * (Mhh+Mss);
    const double B = std::sqrt(0.25 * square(Mhh - Mss) + square(Mhs));
    const double mH1 = A - B;
    const double mH2 = A + B;
    
    // NG
    const double mg = muh_sq_tree_EWSB + lambda_h * square(h) + 0.5 * lambda_hs * square(s);
    
    return {mH1, mH2, mg};
  }

  std::vector<double> get_scalar_dofs() const override { return {1., 1., 3.}; }

  // W, Z 
  std::vector<double> get_vector_masses_sq(Eigen::VectorXd phi) const override {
    const double h_sq = square(phi[0]);
    const double mW_ = 0.25 * square(SM::g) * h_sq;
    const double mZ_ = 0.25 * (square(SM::g) + square(SM::gp)) * h_sq;
    const double mPh_Sq = 0;
    return {mW_, mZ_, mPh_Sq};
  }

  std::vector<double> get_vector_dofs() const override { return {6., 3., 2.}; }
  
  // top
  std::vector<double> get_fermion_masses_sq(Eigen::VectorXd phi) const override{
    return {0.5 * SM::yt_sq * square(phi[0])};
  }
  // top, bottom and tau
  std::vector<double> get_fermion_dofs() const override {
    return {12.};
  }

  size_t get_n_scalars() const override { return 2; }

  std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
    auto phi1 = phi;
    phi1[0] = - phi[0];
    auto phi2 = phi;
    phi2[1] = - phi[1];
    return {phi1, phi2};
  };

 private:
  double m_s;
  double lambda_hs;
  double muh_sq;
  double lambda_h;
  double mus_sq;
  double lambda_s;
  
  double muh_sq_tree_EWSB;
  double mus_sq_tree_EWSB;  
};

xSM_MSbar make_xSM(double lambda_hs, double Q) {
  const double ms = 0.5 * SM::mh;

  const double lambda_s = 2. / square(SM::mh * SM::v) * square(square(ms) 
                          - 0.5 * lambda_hs * square(SM::v)) + 0.1;

  // Set ms, c1, c2, b3, vs, lambda_s
  xSM_MSbar_parameters_solver solver(ms, 0, 0.25 * lambda_hs, 0, 0, 0.25 * lambda_s);
  solver.set_renormalization_scale(Q);

  // Solve mu_h, lambda_h and mu_s
  const auto valid = solver.solving_parameters();
  if (!valid) {
    throw std::runtime_error("Invalid when solving parameters");
  }

  double mu_h_Sq = 2. * solver.get_mu_h_Sq();
  double lambda_h = 4. * solver.get_lambda_h();
  double mu_s_Sq = 2. * solver.get_mu_s_Sq();

  // Construct our model
  EffectivePotential::xSM_MSbar model(ms, lambda_hs, mu_h_Sq,
                                      lambda_h, mu_s_Sq, lambda_s,
                                      mu_h_Sq,
                                      mu_s_Sq);
  model.set_renormalization_scale(Q);
  return model;
}

}  // namespace EffectivePotential

#endif
