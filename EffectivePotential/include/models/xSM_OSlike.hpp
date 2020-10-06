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

#ifndef POTENTIAL_XSM_OS_LIKE_HPP_INCLUDED
#define POTENTIAL_XSM_OS_LIKE_HPP_INCLUDED

/**
   Z2 symmetric real scalar singlet extension of the Standard Model
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <boost/math/tools/roots.hpp>

#include "one_loop_potential.hpp"
#include "pow.hpp"
#include "SM_parameters.hpp"

namespace EffectivePotential {

struct TerminationCondition  {
  bool operator() (double min, double max)  {
    return abs(min - max) <= 1E-3;
  }
};

class xSM_OSlike : public OneLoopPotential {
 public:
  void set_m_s(double m_s_) {
    m_s = m_s_;
    set_mus_and_lambdas();
  }
  void set_lambda_hs(double lambda_hs_) {
    lambda_hs = lambda_hs_;
    set_mus_and_lambdas();
  }

  double VCW(Eigen::VectorXd phi, double mu){
    double correction = 0;
    Q = mu;
    const std::vector<double> scalar_masses_sq = get_scalar_masses_sq(phi, 0);
    const std::vector<double> fermion_masses_sq = get_fermion_masses_sq(phi);
    const std::vector<double> vector_masses_sq = get_vector_masses_sq(phi);
    
    static const auto scalar_dofs = get_scalar_dofs();
    static const auto fermion_dofs = get_fermion_dofs();
    static const auto vector_dofs = get_vector_dofs();
    
    // scalar correction
    for (size_t i = 0; i < scalar_masses_sq.size(); ++i) {
      const double x = scalar_masses_sq[i] / square(Q);
      correction += scalar_dofs[i] * scalar_masses_sq[i] *
                    (square(Q) * xlogx(x) - scalar_masses_sq[i] * 3. / 2.);
    }

    // fermion correction
    for (size_t i = 0; i < fermion_masses_sq.size(); ++i) {
      const double x = fermion_masses_sq[i] / square(Q);
      correction -= fermion_dofs[i] * fermion_masses_sq[i] *
                    (square(Q) * xlogx(x) - fermion_masses_sq[i] * 3. / 2.);
    }

    // vector correction
    for (size_t i = 0; i < vector_masses_sq.size(); ++i) {
      const double x = vector_masses_sq[i] / square(Q);
      correction += vector_dofs[i] * vector_masses_sq[i] *
                    (square(Q) * xlogx(x) - vector_masses_sq[i] * 5. / 6.);
    }
    
    return correction / (64. * M_PI * M_PI);
  }
  
  double dVCW(double mu){
    // first derivate
    Eigen::VectorXd phi_1(2);
    phi_1 << v+eps, 0;
    Eigen::VectorXd phi_2(2);
    phi_2 << v-eps, 0;
    
    double d1V = (VCW(phi_1, mu) - VCW(phi_2, mu))/(2*eps);
    return d1V;
  }
  
  double d2VCW(double mu){
    // second derivate
    const std::vector<double> n_h_xx = {-1., 0., 1.};
    const std::vector<double> coeff_xx = {1., -2., 1.};
    
    double d2V=0;
    Eigen::VectorXd x(2);
    x << v, 0;
    Eigen::VectorXd phi_shifted=x;
    for (int jj = 0; jj < n_h_xx.size(); ++jj) {
      phi_shifted[0] = x[0] + n_h_xx[jj] * eps;
      d2V += VCW(phi_shifted, mu) * coeff_xx[jj] / square(eps);
    }
    return d2V;
  }

  double d2VCWds2(double mu){
    // second derivate
    const std::vector<double> n_h_xx = {-1., 0., 1.};
    const std::vector<double> coeff_xx = {1., -2., 1.};
    
    double d2V=0;
    Eigen::VectorXd x(2);
    x << v, 0;
    Eigen::VectorXd phi_shifted=x;
    for (int jj = 0; jj < n_h_xx.size(); ++jj) {
      phi_shifted[1] = x[1] + n_h_xx[jj] * eps;
      d2V += VCW(phi_shifted, mu) * coeff_xx[jj] / square(eps);
    }
    return d2V;
  }
  
  double renormalization_condition_1(double mu){
    return d2VCW(mu) - dVCW(mu)/v;
  }
  

  void solve_Q(){
    
    const auto f = [this](double mu) {return this->renormalization_condition_1(mu);};
    const auto result = boost::math::tools::bisect(f, 10., 500., TerminationCondition());
    const double mu = (result.first + result.second) * 0.5;
    Q = mu;
    
    Eigen::VectorXd x(2);
    x << v, 0;
    const std::vector<double> scalar_masses_sq = get_scalar_masses_sq(x, 0);
    scalar_masses_sq_EW[2] = scalar_masses_sq[2];
    
//    std::cout << "scalar_masses_sq    = "<< scalar_masses_sq[0] << std::endl;
//    std::cout << "scalar_masses_sq    = "<< scalar_masses_sq[1] << std::endl;
//    std::cout << "scalar_masses_sq    = "<< scalar_masses_sq[2] << std::endl;
//    std::cout << "scalar_masses_sq_EW    = "<< scalar_masses_sq_EW[0] << std::endl;
//    std::cout << "scalar_masses_sq_EW    = "<< scalar_masses_sq_EW[1] << std::endl;
//    std::cout << "scalar_masses_sq_EW    = "<< scalar_masses_sq_EW[2] << std::endl;
    
//    delta_muSqH = d2VCW(Q);
//    delta_muSqS = d2VCWds2(Q);
//    std::cout << "Q    = "<< Q << std::endl;
//    std::cout << "delta_muSqH    = "<< delta_muSqH << std::endl;
//    std::cout << "delta_muSqS    = "<< delta_muSqS << std::endl;
//    std::cout << "renormalization condition 1    = "<< renormalization_condition_1(Q) << std::endl;
  }
  
  double delta_muSqH = 0;
  double delta_muSqS = 0;
  
  double counter_term(Eigen::VectorXd phi, double T) const override{
//    return -0.5*delta_muSqH*square(phi[0]) - 0.5*delta_muSqS*square(phi[1]);
    return 0;
  }
  
  double V0(Eigen::VectorXd phi) const override {
    return 0.5 * muh_sq * square(phi[0]) +
           0.25 * lambda_h * pow_4(phi[0]) +
           0.25 * lambda_hs * square(phi[0]) * square(phi[1]) +
           0.5 * mus_sq * square(phi[1]) +
           0.25 * lambda_s * pow_4(phi[1]);
  }

  double V1(std::vector<double> scalar_masses_sq,
            std::vector<double> fermion_masses_sq,
            std::vector<double> vector_masses_sq) const override {
    double correction = 0;

    static const auto scalar_dofs = get_scalar_dofs();
    static const auto fermion_dofs = get_fermion_dofs();
    static const auto vector_dofs = get_vector_dofs();
    
    // scalar correction
    for (size_t i = 0; i < scalar_masses_sq.size(); ++i) {
      const double x = scalar_masses_sq[i] / scalar_masses_sq_EW[i];
      correction += scalar_dofs[i] * scalar_masses_sq[i] *
                    (scalar_masses_sq_EW[i] * xlogx(x) - scalar_masses_sq[i] * 3. / 2.);
      correction += scalar_dofs[i] * 2. * scalar_masses_sq[i] * scalar_masses_sq_EW[i];
    }

    // fermion correction
    for (size_t i = 0; i < fermion_masses_sq.size(); ++i) {
      const double x = fermion_masses_sq[i] / fermion_masses_sq_EW[i];
      correction -= fermion_dofs[i] * fermion_masses_sq[i] *
                    (fermion_masses_sq_EW[i] * xlogx(x) - fermion_masses_sq[i] * 3. / 2.);
      correction -= fermion_dofs[i] * 2. * fermion_masses_sq[i] * fermion_masses_sq_EW[i];
    }

    // vector correction
    for (size_t i = 0; i < vector_masses_sq.size(); ++i) {
      const double x = vector_masses_sq[i] / vector_masses_sq_EW[i];
      correction += vector_dofs[i] * vector_masses_sq[i] *
                    (vector_masses_sq_EW[i] * xlogx(x) - vector_masses_sq[i] * 3. / 2.);
      correction += vector_dofs[i] * 2. * vector_masses_sq[i] * vector_masses_sq_EW[i];
    }
    return correction / (64. * M_PI * M_PI);                  
  }

  // Higgs
  std::vector<double> get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const override{
    const double h = phi[0];
    const double s = phi[1];
    // M is the mass matrix. h is higgs direction, s is singlet
    const double Mhh = muh_sq + 3*lambda_h*square(h) + 0.5*lambda_hs*square(s)
                      + square(T) * (3./16.*g_sq + 1./16.*gp_sq + 0.25*yt_sq
                      + 0.5*lambda_h + 1./24.*lambda_hs);
    const double Mss = mus_sq + 3*lambda_s*square(s) + 0.5*lambda_hs*square(h)
                      + square(T) * (1./6.*lambda_hs + 0.25*lambda_s);
    const double Mhs = lambda_hs*h*s;
    // diagonalization
    const double A = 0.5 * (Mhh+Mss);
    const double B = std::sqrt(0.25 * square(Mhh - Mss) + square(Mhs));
    const double mH1 = A-B;
    const double mH2 = A+B;
    
    // resummed NG contributions
    const double mg = muh_sq + lambda_h*square(h) + 0.5*lambda_hs*square(s)
                      + square(T) * (3./16.*g_sq + 1./16.*gp_sq + 0.25*yt_sq
                      + 0.5*lambda_h + 1./24.*lambda_hs);
    const auto fm2 = get_fermion_masses_sq(phi);
    const auto vm2 = get_vector_masses_sq(phi);
    
    const double Q2 = Q*Q;
    const double sum = 1. / (16. * M_PI * M_PI) * (
                       3.  * lambda_h     * (Q2*xlogx(mH2/Q2)     - mH2)
                      +0.5 * lambda_hs    * (Q2*xlogx(mH1/Q2)     - mH1)
                      -6.  * yt_sq        * (Q2*xlogx(fm2[0]/Q2)  - fm2[0])
                      +1.5 * g_sq         * (Q2*xlogx(vm2[0]/Q2) - 1./3.*vm2[0])
                      +0.75* (g_sq+gp_sq) * (Q2*xlogx(vm2[1]/Q2) - 1./3.*vm2[1])
                      );
    const double Mg = mg + sum;
    
    return {mH1, mH2, Mg};
  }
  std::vector<double> get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const override {
    return get_scalar_debye_sq(phi, xi, 0.);
  }
  std::vector<double> get_scalar_dofs() const override { return {1., 1., 3.}; }

  // W, Z
  std::vector<double> get_vector_debye_sq(Eigen::VectorXd phi, double T) const override{
    const double h_sq = square(phi[0]);
    const double A = (g_sq + gp_sq) * (11./12.*square(T) + 0.125 * h_sq);
    const double B = 0.125 * sqrt(square(g_sq - gp_sq) *
                               (16 * square(11./6.) * pow_4(T) + 8. * (11./6.) * square(T) * h_sq) +
                               square(g_sq + gp_sq) * square(h_sq));

    const double W_debye = g_sq * (0.25 * h_sq + 11./6. * square(T));
    const double Z_debye = A + B;
    const double g_debye = A - B;
    return {W_debye, Z_debye, g_debye};
  }
  std::vector<double> get_vector_masses_sq(Eigen::VectorXd phi) const override{
    return get_vector_debye_sq(phi, 0.);
  }
  std::vector<double> get_vector_dofs() const override { return {6., 3., 3.}; }
  
  // top
  std::vector<double> get_fermion_masses_sq(Eigen::VectorXd phi) const override{
    return {0.5*yt_sq*square(phi[0])};
  }
  // top, bottom and tau
  std::vector<double> get_fermion_dofs() const override {
    return {12.};
  }

  size_t get_n_scalars() const override {return 2;}

  std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
    auto phi1 = phi;
    phi1[0] = - phi[0];
    auto phi2 = phi;
    phi2[1] = - phi[1];
    return {phi1,phi2};
  };


 private:
  
  const double v = SM::v;
  const double mh = SM::mh;
  const double mtop = SM::mtop;
  const double mZ = SM::mZ;
  const double mW = SM::mW;
  const double g = SM::g;
  const double g_sq = g * g;
  const double gp = SM::gp;
  const double gp_sq = gp * gp;
  const double yt_sq = SM::yt_sq;
  
  double Q=mtop;       // renormalization scale
  double eps = 0.001;  // for solving renormalization scale
  
  const double muh_sq = -0.5 * mh * mh;
  const double lambda_h = -muh_sq / square(v);

  std::vector<double> scalar_masses_sq_EW = {mh*mh/4., mh*mh, 0.0};
  const std::vector<double> vector_masses_sq_EW = {mW*mW, mZ*mZ};
  const std::vector<double> fermion_masses_sq_EW = {mtop*mtop};
  
  double lambda_hs = 0.25;
  double m_s = 65;
  
  double mus_sq = square(m_s) - lambda_hs * square(v) / 2.;
  double lambda_s = lambda_h*square(mus_sq)/square(muh_sq) +0.1;
  
  void set_mus_and_lambdas(){
    mus_sq = square(m_s) - lambda_hs * square(v) / 2.;
    lambda_s = lambda_h*square(mus_sq)/square(muh_sq) +0.1;
  }
  
};

}  // namespace EffectivePotential

#endif
