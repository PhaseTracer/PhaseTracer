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
   See arXiv:2208.01319  [hep-ph] for details
*/

#include <boost/math/tools/roots.hpp>
#include <Eigen/Eigenvalues>
#include "xSM_base.hpp"

namespace EffectivePotential {

struct TerminationCondition {
  bool operator()(double min, double max) {
    return abs(min - max) <= 1E-3;
  }
};

class xSM_OSlike : public xSM_base {
public:
  xSM_OSlike(double lambda_hs_,
             double lambda_s_,
             double ms_) {
    lambda_hs = lambda_hs_;
    lambda_s = lambda_s_;
    ms = ms_;

    double mhh2 = square(SM_mh);
    double mss2 = square(ms);
    lambda_h = mhh2 / (2. * square(SM_v));
    muh_sq = -lambda_h * square(SM_v);
    mus_sq = mss2 - 0.5 * lambda_hs * square(SM_v);

    solve_renormalization_scale();
  }

  double VCW(Eigen::VectorXd phi, double mu) {
    double correction = 0;
    Q = mu;
    const std::vector<double> scalar_masses_sq = get_scalar_masses_sq(phi);
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

  double dVCW(double mu) {
    // first derivate
    Eigen::VectorXd phi_1(2);
    phi_1 << SM_v + eps, 0;
    Eigen::VectorXd phi_2(2);
    phi_2 << SM_v - eps, 0;

    double d1V = (VCW(phi_1, mu) - VCW(phi_2, mu)) / (2 * eps);
    return d1V;
  }

  double d2VCW(double mu) {
    // second derivate
    const std::vector<double> n_h_xx = {-1., 0., 1.};
    const std::vector<double> coeff_xx = {1., -2., 1.};

    double d2V = 0;
    Eigen::VectorXd x(2);
    x << SM_v, 0;
    Eigen::VectorXd phi_shifted = x;
    for (int jj = 0; jj < n_h_xx.size(); ++jj) {
      phi_shifted[0] = x[0] + n_h_xx[jj] * eps;
      d2V += VCW(phi_shifted, mu) * coeff_xx[jj] / square(eps);
    }
    return d2V;
  }

  double d2VCWds2(double mu) {
    // second derivate
    const std::vector<double> n_h_xx = {-1., 0., 1.};
    const std::vector<double> coeff_xx = {1., -2., 1.};

    double d2V = 0;
    Eigen::VectorXd x(2);
    x << SM_v, 0;
    Eigen::VectorXd phi_shifted = x;
    for (int jj = 0; jj < n_h_xx.size(); ++jj) {
      phi_shifted[1] = x[1] + n_h_xx[jj] * eps;
      d2V += VCW(phi_shifted, mu) * coeff_xx[jj] / square(eps);
    }
    return d2V;
  }

  double renormalization_condition_1(double mu) {
    return d2VCW(mu) - dVCW(mu) / SM_v;
  }

  void solve_renormalization_scale() {

    const auto f = [this](double mu) { return this->renormalization_condition_1(mu); };
    const auto result = boost::math::tools::bisect(f, 10., 500., TerminationCondition());
    const double mu = (result.first + result.second) * 0.5;
    Q = mu;

    Eigen::VectorXd EWVEV(2);
    EWVEV << SM_v, 0;
    scalar_masses_sq_EW = get_scalar_masses_sq(EWVEV);
    fermion_masses_sq_EW = get_fermion_masses_sq(EWVEV);
    vector_masses_sq_EW = get_vector_masses_sq(EWVEV);

    //    for (int ii=0; ii < scalar_masses_sq.size(); ++ii){
    //      std::cout << "scalar_masses_EW["<< ii << "] = "<< std::sqrt(std::abs(scalar_masses_sq_EW[ii])) << std::endl;
    //    }
    //    for (int ii=0; ii < vector_masses_sq_EW.size(); ++ii){
    //      std::cout << "vector_masses_EW["<< ii << "] = "<< std::sqrt(std::abs(vector_masses_sq_EW[ii])) << std::endl;
    //    }

    //    delta_muSqH = d2VCW(Q);
    //    delta_muSqS = d2VCWds2(Q);
    //    std::cout << "Q    = "<< Q << std::endl;
    //    std::cout << "delta_muSqH    = "<< delta_muSqH << std::endl;
    //    std::cout << "delta_muSqS    = "<< delta_muSqS << std::endl;
    //    std::cout << "renormalization condition 1    = "<< renormalization_condition_1(Q) << std::endl;
  }

  double V1(std::vector<double> scalar_masses_sq,
            std::vector<double> fermion_masses_sq,
            std::vector<double> vector_masses_sq,
            std::vector<double> ghost_masses_sq) const override {
    double correction = 0;

    static const auto scalar_dofs = get_scalar_dofs();
    static const auto fermion_dofs = get_fermion_dofs();
    static const auto vector_dofs = get_vector_dofs();

    // scalar correction
    for (size_t i = 0; i < scalar_masses_sq.size(); ++i) {
      const double x = scalar_masses_sq[i] / scalar_masses_sq_EW[i];
      correction += scalar_dofs[i] * scalar_masses_sq[i] *
                    (scalar_masses_sq_EW[i] * xlogx(x) - scalar_masses_sq[i] * 3. / 2.);
      correction += scalar_dofs[i] * (2. * scalar_masses_sq[i] * scalar_masses_sq_EW[i] - 0.5 * scalar_masses_sq_EW[i] * scalar_masses_sq_EW[i]);
    }

    // fermion correction
    for (size_t i = 0; i < fermion_masses_sq.size(); ++i) {
      const double x = fermion_masses_sq[i] / fermion_masses_sq_EW[i];
      correction -= fermion_dofs[i] * fermion_masses_sq[i] *
                    (fermion_masses_sq_EW[i] * xlogx(x) - fermion_masses_sq[i] * 3. / 2.);
      correction -= fermion_dofs[i] * (2. * fermion_masses_sq[i] * fermion_masses_sq_EW[i] - 0.5 * fermion_masses_sq_EW[i] * fermion_masses_sq_EW[i]);
    }

    // vector correction
    // deleted photon contribution as it has no effect
    for (size_t i = 0; i < vector_masses_sq.size(); ++i) {
      if (i != 2 and i != 5) {
        const double x = vector_masses_sq[i] / vector_masses_sq_EW[i];
        correction += vector_dofs[i] * vector_masses_sq[i] *
                      (vector_masses_sq_EW[i] * xlogx(x) - vector_masses_sq[i] * 3. / 2.);
        correction += vector_dofs[i] * (2. * vector_masses_sq[i] * vector_masses_sq_EW[i] - 0.5 * vector_masses_sq_EW[i] * vector_masses_sq_EW[i]);
      }
    }
    return correction / (64. * M_PI * M_PI);
  }

  // Higgs
  std::vector<double> get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const override {
    const double h = phi[0];
    const double s = phi[1];
    const auto thermal_sq = get_scalar_thermal_sq(T);

    const double mhh2 = muh_sq + 3. * lambda_h * square(h) + 0.5 * lambda_hs * square(s);
    const double mgg2 = muh_sq + lambda_h * square(h) + 0.5 * lambda_hs * square(s);
    const double mss2 = mus_sq + 3. * lambda_s * square(s) + 0.5 * lambda_hs * square(h);

    // resummed Goldstone contributions
    const auto fm2 = get_fermion_masses_sq(phi);
    const auto vm2 = get_vector_masses_sq(phi);
    const double Qsq = square(Q);
    const double sum = 1. / (16. * M_PI * M_PI) * (3. * lambda_h * (Qsq * xlogx(mhh2 / Qsq) - mhh2) + 0.5 * lambda_hs * (Qsq * xlogx(mss2 / Qsq) - mss2) - 6. * SM_yt_sq * (Qsq * xlogx(fm2[0] / Qsq) - fm2[0]) - 6. * SM_yb_sq * (Qsq * xlogx(fm2[1] / Qsq) - fm2[1])
                                                   - 2. * SM_ytau_sq * (Qsq * xlogx(fm2[2] / Qsq) - fm2[2])
                                                   + 1.5 * square(SM_g) * (Qsq * xlogx(vm2[0] / Qsq) - 1. / 3. * vm2[0]) + 0.75 * (square(SM_g) + square(SM_gp)) * (Qsq * xlogx(vm2[1] / Qsq) - 1. / 3. * vm2[1]));

    // Goldstone finite temperature masses
    double mTG02 = mgg2 + thermal_sq[0] + sum;
    double mTGpm2 = mTG02; // 2 degrees of freedom or two degenerate copies
    // CP even Higgs thermal temperature masses
    Eigen::MatrixXd MTH2 = Eigen::MatrixXd::Zero(2, 2);
    MTH2(0, 0) = mhh2 + thermal_sq[0];
    MTH2(1, 1) = mss2 + thermal_sq[1];
    // Mixing between Higgs and singlet
    MTH2(0, 1) = MTH2(1, 0) = lambda_hs * h * s;
    // get eigenvalues
    const Eigen::VectorXd mH_sq = MTH2.eigenvalues().real();
    // vector for all scalars, including two mass degenerate charged goldstones
    std::vector<double> m_sq_vector{mH_sq(0), mH_sq(1), mTG02, mTGpm2, mTGpm2};
    // mass order
    std::sort(m_sq_vector.begin(), m_sq_vector.end());
    return m_sq_vector;
  }
  std::vector<double> get_scalar_masses_sq(Eigen::VectorXd phi, double xi = 0) const override {
    return get_scalar_debye_sq(phi, xi, 0.);
  }

  double get_renormalization_scale() {
    return Q;
  }

  double get_muh_sq() const { return muh_sq; }
  double get_mus_sq() const { return mus_sq; }
  double get_lambda_h() const { return lambda_h; }
  double get_lambda_s() const { return lambda_s; }
  double get_lambda_hs() const { return lambda_hs; }

private:
  double Q = 173.;    // renormalization scale
  double eps = 0.001; // for solving renormalization scale
  std::vector<double> scalar_masses_sq_EW = {0., 0., 0., 0., 0., 0.};
  std::vector<double> vector_masses_sq_EW = {0., 0.};
  std::vector<double> fermion_masses_sq_EW = {0., 0., 0.};
};

} // namespace EffectivePotential

#endif
