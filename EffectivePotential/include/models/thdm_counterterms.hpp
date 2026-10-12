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

#ifndef POTENTIAL_THDM_COUNTERTERMS_HPP_INCLUDED
#define POTENTIAL_THDM_COUNTERTERMS_HPP_INCLUDED

/*
  Coleman-Weinberg counterterms for the CP-conserving 2HDM, following the BSMPT
  on-shell prescription (Eqs. (3.59)-(3.63) of arXiv:1803.02846).

  This is a port of TransitionListener's `counterterms/cw.py` and
  `counterterms/twohdm_solver.py`. Seven renormalisation conditions are imposed
  at the tree-level vacuum in the 8-component field basis: the two tadpoles, the
  three CP-even Hessian entries, one CP-odd entry and one charged entry. They
  fix (dm11, dm22, dm12, dl1, dl2, dl3, dl5); dl4 is not renormalised.

  Goldstone infrared divergences are handled the way BSMPT does: eigenvalues
  below a threshold are set to exactly zero and then skipped in the single-log
  terms, while the double sums use the regulated kernel `fbase`.
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <memory>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

#include "thdm_curvature_tensors.hpp"

namespace EffectivePotential {
namespace THDM_ct {

using THDM_tensors::NCT;
using THDM_tensors::NG;
using THDM_tensors::NH;
using THDM_tensors::NL;
using THDM_tensors::NQ;

using cdouble = std::complex<double>;

/** Coleman-Weinberg subtraction constants, matching TL's CW_CONSTANTS. */
constexpr double CW_SCALAR = 3. / 2.;
constexpr double CW_FERMION = 3. / 2.;
constexpr double CW_GAUGE = 5. / 6.;

/** Eigenvalues closer to zero than this are treated as exactly massless. */
constexpr double SCALAR_THRESHOLD = 1e-5;
constexpr double FERMION_THRESHOLD = 1e-10;

struct Result {
  double dm11_sq = 0.;
  double dm22_sq = 0.;
  double dm12_sq = 0.;
  double dlambda1 = 0.;
  double dlambda2 = 0.;
  double dlambda3 = 0.;
  double dlambda4 = 0.;
  double dlambda5 = 0.;
  /** Largest violation of the renormalisation conditions; should be ~0. */
  double residual = 0.;
  /** Coleman-Weinberg derivatives in the physical basis, for diagnostics. */
  Eigen::VectorXd cw_gradient;
  Eigen::MatrixXd cw_hessian;
};

namespace detail {

/** log(m^2/mu^2) - c + 1/2, regularised at m^2 = 0. */
inline double log_term(double m_sq, double scale_sq, double c) {
  if (m_sq == 0.) {
    return -c + 0.5;
  }
  return std::log(m_sq / scale_sq) - c + 0.5;
}

/** BSMPT's Class_Potential_Origin::fbase. */
inline double fbase(double m_sq_a, double m_sq_b, double scale) {
  if (m_sq_a == 0. && m_sq_b == 0.) {
    return 1.;
  }
  const double tol = 1e-5;
  const double log_scale = 2. * std::log(scale);
  const double log_a = (m_sq_a != 0.) ? std::log(m_sq_a) - log_scale : 0.;
  if (std::abs(m_sq_a - m_sq_b) > tol) {
    const double log_b = (m_sq_b != 0.) ? std::log(m_sq_b) - log_scale : 0.;
    if (m_sq_a == 0.) {
      return log_b;
    }
    if (m_sq_b == 0.) {
      return log_a;
    }
    return (log_a * m_sq_a - log_b * m_sq_b) / (m_sq_a - m_sq_b);
  }
  return 1. + log_a;
}

/**
 * One fermion sector (quarks or leptons) in the mass eigenbasis.
 *
 * `c21[(a*n + b)*NH + i]` and `d22[(a*NH + i)*NH + j]` are the rotated
 * two-fermion/one-Higgs and two-fermion/two-Higgs couplings; only the
 * fermion-diagonal part of the latter is ever needed.
 */
struct FermionSector {
  int n = 0;
  std::vector<double> mass_sq;
  std::vector<cdouble> c21;
  std::vector<cdouble> d22;
};

/**
 * Diagonalise M^dagger M for a Yukawa curvature tensor and build the rotated
 * couplings. `y[(i*n + j)*NH + k]` is the curvature tensor, `rot_h` the Higgs
 * rotation with rot_h[a*NH + i] the i-th component of the a-th eigenvector.
 */
inline FermionSector build_fermion_sector(const std::vector<cdouble> &y, int n,
                                          const std::array<double, NH> &vev,
                                          const std::vector<double> &rot_h) {
  FermionSector sector;
  sector.n = n;

  // M_ij = Y_ijk v_k
  Eigen::MatrixXcd mij = Eigen::MatrixXcd::Zero(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      cdouble total(0., 0.);
      for (int k = 0; k < NH; ++k) {
        total += y[(i * n + j) * NH + k] * vev[k];
      }
      mij(i, j) = total;
    }
  }

  const Eigen::MatrixXcd mass = mij.adjoint() * mij;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(mass);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Failed to diagonalise the fermion mass matrix");
  }

  sector.mass_sq.assign(n, 0.);
  for (int a = 0; a < n; ++a) {
    const double value = solver.eigenvalues()(a);
    sector.mass_sq[a] = (std::abs(value) < FERMION_THRESHOLD) ? 0. : value;
  }

  // rot[a][i] is component i of eigenvector a, i.e. the transpose of Eigen's
  // eigenvector matrix, matching TL's `eigvecs.T`.
  std::vector<cdouble> rot(static_cast<size_t>(n) * n);
  for (int a = 0; a < n; ++a) {
    for (int i = 0; i < n; ++i) {
      rot[a * n + i] = solver.eigenvectors()(i, a);
    }
  }

  // Lambda3_ijk = conj(Y_ilk) M_lj + conj(M_il) Y_ljk
  std::vector<cdouble> lambda3(static_cast<size_t>(n) * n * NH, cdouble(0., 0.));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < NH; ++k) {
        cdouble total(0., 0.);
        for (int l = 0; l < n; ++l) {
          total += std::conj(y[(i * n + l) * NH + k]) * mij(l, j);
          total += std::conj(mij(i, l)) * y[(l * n + j) * NH + k];
        }
        lambda3[(i * n + j) * NH + k] = total;
      }
    }
  }

  // Lambda4_ijkm = conj(Y_ilk) Y_ljm + conj(Y_ilm) Y_ljk
  std::vector<cdouble> lambda4(static_cast<size_t>(n) * n * NH * NH,
                               cdouble(0., 0.));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < NH; ++k) {
        for (int mm = 0; mm < NH; ++mm) {
          cdouble total(0., 0.);
          for (int l = 0; l < n; ++l) {
            total += std::conj(y[(i * n + l) * NH + k]) * y[(l * n + j) * NH + mm];
            total += std::conj(y[(i * n + l) * NH + mm]) * y[(l * n + j) * NH + k];
          }
          lambda4[((i * n + j) * NH + k) * NH + mm] = total;
        }
      }
    }
  }

  // c21_pqr = rot_pa rot_qb rot_h_ri Lambda3_abi (rotations are not conjugated,
  // following BSMPT).
  std::vector<cdouble> stage1(static_cast<size_t>(n) * n * NH, cdouble(0., 0.));
  for (int p = 0; p < n; ++p) {
    for (int b = 0; b < n; ++b) {
      for (int i = 0; i < NH; ++i) {
        cdouble total(0., 0.);
        for (int a = 0; a < n; ++a) {
          total += rot[p * n + a] * lambda3[(a * n + b) * NH + i];
        }
        stage1[(p * n + b) * NH + i] = total;
      }
    }
  }
  std::vector<cdouble> stage2(static_cast<size_t>(n) * n * NH, cdouble(0., 0.));
  for (int p = 0; p < n; ++p) {
    for (int q = 0; q < n; ++q) {
      for (int i = 0; i < NH; ++i) {
        cdouble total(0., 0.);
        for (int b = 0; b < n; ++b) {
          total += rot[q * n + b] * stage1[(p * n + b) * NH + i];
        }
        stage2[(p * n + q) * NH + i] = total;
      }
    }
  }
  sector.c21.assign(static_cast<size_t>(n) * n * NH, cdouble(0., 0.));
  for (int p = 0; p < n; ++p) {
    for (int q = 0; q < n; ++q) {
      for (int r = 0; r < NH; ++r) {
        cdouble total(0., 0.);
        for (int i = 0; i < NH; ++i) {
          total += rot_h[r * NH + i] * stage2[(p * n + q) * NH + i];
        }
        sector.c21[(p * n + q) * NH + r] = total;
      }
    }
  }

  // d22_aij = rot_ab rot_ac rot_h_im rot_h_jn Lambda4_bcmn
  sector.d22.assign(static_cast<size_t>(n) * NH * NH, cdouble(0., 0.));
  std::vector<cdouble> t1(static_cast<size_t>(n) * NH * NH);
  std::array<std::array<cdouble, NH>, NH> t2{};
  std::array<std::array<cdouble, NH>, NH> t3{};
  for (int a = 0; a < n; ++a) {
    for (int c = 0; c < n; ++c) {
      for (int mm = 0; mm < NH; ++mm) {
        for (int nn = 0; nn < NH; ++nn) {
          cdouble total(0., 0.);
          for (int b = 0; b < n; ++b) {
            total += rot[a * n + b] * lambda4[((b * n + c) * NH + mm) * NH + nn];
          }
          t1[(c * NH + mm) * NH + nn] = total;
        }
      }
    }
    for (int mm = 0; mm < NH; ++mm) {
      for (int nn = 0; nn < NH; ++nn) {
        cdouble total(0., 0.);
        for (int c = 0; c < n; ++c) {
          total += rot[a * n + c] * t1[(c * NH + mm) * NH + nn];
        }
        t2[mm][nn] = total;
      }
    }
    for (int i = 0; i < NH; ++i) {
      for (int nn = 0; nn < NH; ++nn) {
        cdouble total(0., 0.);
        for (int mm = 0; mm < NH; ++mm) {
          total += rot_h[i * NH + mm] * t2[mm][nn];
        }
        t3[i][nn] = total;
      }
    }
    for (int i = 0; i < NH; ++i) {
      for (int j = 0; j < NH; ++j) {
        cdouble total(0., 0.);
        for (int nn = 0; nn < NH; ++nn) {
          total += rot_h[j * NH + nn] * t3[i][nn];
        }
        sector.d22[(a * NH + i) * NH + j] = total;
      }
    }
  }

  return sector;
}

} // namespace detail

/**
 * Coleman-Weinberg gradient and Hessian in the physical basis, and the
 * counterterms that cancel them.
 *
 * @param p     model parameters in GeV units
 * @param yukawa_type 1..4
 * @param scale MS-bar renormalisation scale mu
 */
inline Result compute(const THDM_tensors::Parameters &p, int yukawa_type,
                      double scale) {
  using namespace detail;

  const double scale_sq = scale * scale;

  std::array<double, NH> vev{};
  vev[4] = p.v1;
  vev[6] = p.v2;

  // ------------------------------------------------------------------
  // Curvature tensors
  // ------------------------------------------------------------------
  THDM_tensors::Tensor2 h2{};
  auto h4 = std::unique_ptr<THDM_tensors::Tensor4>(new THDM_tensors::Tensor4());
  auto gauge =
      std::unique_ptr<THDM_tensors::GaugeTensor>(new THDM_tensors::GaugeTensor());
  auto quark =
      std::unique_ptr<THDM_tensors::QuarkTensor>(new THDM_tensors::QuarkTensor());
  auto lepton =
      std::unique_ptr<THDM_tensors::LeptonTensor>(new THDM_tensors::LeptonTensor());

  THDM_tensors::fill_higgs_l2(p, h2);
  THDM_tensors::fill_higgs_l4(p, *h4);
  THDM_tensors::fill_gauge(p, *gauge);
  THDM_tensors::fill_quark(p, yukawa_type, *quark);
  THDM_tensors::fill_lepton(p, yukawa_type, *lepton);

  // ------------------------------------------------------------------
  // Scalar and gauge mass matrices. The cubic Higgs curvature vanishes
  // identically for the 2HDM, so only the quadratic and quartic terms appear.
  // ------------------------------------------------------------------
  Eigen::MatrixXd mass_h = Eigen::MatrixXd::Zero(NH, NH);
  for (int i = 0; i < NH; ++i) {
    for (int j = 0; j < NH; ++j) {
      double total = h2[i][j];
      for (int k = 0; k < NH; ++k) {
        for (int l = 0; l < NH; ++l) {
          total += 0.5 * (*h4)[i][j][k][l] * vev[k] * vev[l];
        }
      }
      mass_h(i, j) = total;
    }
  }

  Eigen::MatrixXd mass_g = Eigen::MatrixXd::Zero(NG, NG);
  for (int a = 0; a < NG; ++a) {
    for (int b = 0; b < NG; ++b) {
      double total = 0.;
      for (int i = 0; i < NH; ++i) {
        for (int j = 0; j < NH; ++j) {
          total += (*gauge)[a][b][i][j] * vev[i] * vev[j];
        }
      }
      mass_g(a, b) = 0.5 * total;
    }
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver_h(mass_h);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver_g(mass_g);
  if (solver_h.info() != Eigen::Success || solver_g.info() != Eigen::Success) {
    throw std::runtime_error("Failed to diagonalise the boson mass matrices");
  }

  // rot_h[a*NH + i] is component i of the a-th eigenvector (TL's `eigvecs.T`).
  std::vector<double> rot_h(static_cast<size_t>(NH) * NH);
  for (int a = 0; a < NH; ++a) {
    for (int i = 0; i < NH; ++i) {
      rot_h[a * NH + i] = solver_h.eigenvectors()(i, a);
    }
  }
  std::vector<double> rot_g(static_cast<size_t>(NG) * NG);
  for (int a = 0; a < NG; ++a) {
    for (int i = 0; i < NG; ++i) {
      rot_g[a * NG + i] = solver_g.eigenvectors()(i, a);
    }
  }

  std::array<double, NH> mass_h_sq{};
  for (int a = 0; a < NH; ++a) {
    const double value = solver_h.eigenvalues()(a);
    mass_h_sq[a] = (std::abs(value) < SCALAR_THRESHOLD) ? 0. : value;
  }
  std::array<double, NG> mass_g_sq{};
  for (int a = 0; a < NG; ++a) {
    const double value = solver_g.eigenvalues()(a);
    mass_g_sq[a] = (std::abs(value) < SCALAR_THRESHOLD) ? 0. : value;
  }

  // ------------------------------------------------------------------
  // Rotated bosonic couplings
  // ------------------------------------------------------------------
  // LambdaGauge3_abi = G_abij v_j
  std::array<std::array<std::array<double, NH>, NG>, NG> lambda_g3{};
  for (int a = 0; a < NG; ++a) {
    for (int b = 0; b < NG; ++b) {
      for (int i = 0; i < NH; ++i) {
        double total = 0.;
        for (int j = 0; j < NH; ++j) {
          total += (*gauge)[a][b][i][j] * vev[j];
        }
        lambda_g3[a][b][i] = total;
      }
    }
  }

  // LambdaHiggs3_ijk = h4_ijkl v_l
  auto lambda_h3 = std::unique_ptr<std::array<std::array<std::array<double, NH>, NH>, NH>>(
      new std::array<std::array<std::array<double, NH>, NH>, NH>());
  for (int i = 0; i < NH; ++i) {
    for (int j = 0; j < NH; ++j) {
      for (int k = 0; k < NH; ++k) {
        double total = 0.;
        for (int l = 0; l < NH; ++l) {
          total += (*h4)[i][j][k][l] * vev[l];
        }
        (*lambda_h3)[i][j][k] = total;
      }
    }
  }

  // Generic three-index rotation: C_pqr = R1_pa R1_qb R2_ri L_abi.
  const auto rotate3 = [](const std::vector<double> &r1, int n1,
                          const std::vector<double> &r2, int n2,
                          const std::vector<double> &lambda) {
    // stage1_pbi = R1_pa L_abi
    std::vector<double> stage1(static_cast<size_t>(n1) * n1 * n2, 0.);
    for (int pp = 0; pp < n1; ++pp) {
      for (int b = 0; b < n1; ++b) {
        for (int i = 0; i < n2; ++i) {
          double total = 0.;
          for (int a = 0; a < n1; ++a) {
            total += r1[pp * n1 + a] * lambda[(a * n1 + b) * n2 + i];
          }
          stage1[(pp * n1 + b) * n2 + i] = total;
        }
      }
    }
    // stage2_pqi = R1_qb stage1_pbi
    std::vector<double> stage2(static_cast<size_t>(n1) * n1 * n2, 0.);
    for (int pp = 0; pp < n1; ++pp) {
      for (int q = 0; q < n1; ++q) {
        for (int i = 0; i < n2; ++i) {
          double total = 0.;
          for (int b = 0; b < n1; ++b) {
            total += r1[q * n1 + b] * stage1[(pp * n1 + b) * n2 + i];
          }
          stage2[(pp * n1 + q) * n2 + i] = total;
        }
      }
    }
    // out_pqr = R2_ri stage2_pqi
    std::vector<double> out(static_cast<size_t>(n1) * n1 * n2, 0.);
    for (int pp = 0; pp < n1; ++pp) {
      for (int q = 0; q < n1; ++q) {
        for (int r = 0; r < n2; ++r) {
          double total = 0.;
          for (int i = 0; i < n2; ++i) {
            total += r2[r * n2 + i] * stage2[(pp * n1 + q) * n2 + i];
          }
          out[(pp * n1 + q) * n2 + r] = total;
        }
      }
    }
    return out;
  };

  std::vector<double> flat_g3(static_cast<size_t>(NG) * NG * NH);
  for (int a = 0; a < NG; ++a) {
    for (int b = 0; b < NG; ++b) {
      for (int i = 0; i < NH; ++i) {
        flat_g3[(a * NG + b) * NH + i] = lambda_g3[a][b][i];
      }
    }
  }
  std::vector<double> flat_h3(static_cast<size_t>(NH) * NH * NH);
  for (int i = 0; i < NH; ++i) {
    for (int j = 0; j < NH; ++j) {
      for (int k = 0; k < NH; ++k) {
        flat_h3[(i * NH + j) * NH + k] = (*lambda_h3)[i][j][k];
      }
    }
  }

  const std::vector<double> c_gauge21 = rotate3(rot_g, NG, rot_h, NH, flat_g3);
  const std::vector<double> c_higgs3 = rotate3(rot_h, NH, rot_h, NH, flat_h3);

  // Boson-diagonal four-index couplings, D_aij = R_ab R_ac Rh_im Rh_jn T_bcmn.
  const auto rotate4_diag = [](const std::vector<double> &r1, int n1,
                               const std::vector<double> &r2, int n2,
                               const std::vector<double> &tensor) {
    std::vector<double> out(static_cast<size_t>(n1) * n2 * n2, 0.);
    std::vector<double> t1(static_cast<size_t>(n1) * n2 * n2);
    std::vector<double> t2(static_cast<size_t>(n2) * n2);
    std::vector<double> t3(static_cast<size_t>(n2) * n2);
    for (int a = 0; a < n1; ++a) {
      for (int c = 0; c < n1; ++c) {
        for (int mm = 0; mm < n2; ++mm) {
          for (int nn = 0; nn < n2; ++nn) {
            double total = 0.;
            for (int b = 0; b < n1; ++b) {
              total += r1[a * n1 + b] *
                       tensor[((b * n1 + c) * n2 + mm) * n2 + nn];
            }
            t1[(c * n2 + mm) * n2 + nn] = total;
          }
        }
      }
      for (int mm = 0; mm < n2; ++mm) {
        for (int nn = 0; nn < n2; ++nn) {
          double total = 0.;
          for (int c = 0; c < n1; ++c) {
            total += r1[a * n1 + c] * t1[(c * n2 + mm) * n2 + nn];
          }
          t2[mm * n2 + nn] = total;
        }
      }
      for (int i = 0; i < n2; ++i) {
        for (int nn = 0; nn < n2; ++nn) {
          double total = 0.;
          for (int mm = 0; mm < n2; ++mm) {
            total += r2[i * n2 + mm] * t2[mm * n2 + nn];
          }
          t3[i * n2 + nn] = total;
        }
      }
      for (int i = 0; i < n2; ++i) {
        for (int j = 0; j < n2; ++j) {
          double total = 0.;
          for (int nn = 0; nn < n2; ++nn) {
            total += r2[j * n2 + nn] * t3[i * n2 + nn];
          }
          out[(a * n2 + i) * n2 + j] = total;
        }
      }
    }
    return out;
  };

  std::vector<double> flat_gauge(static_cast<size_t>(NG) * NG * NH * NH);
  for (int a = 0; a < NG; ++a) {
    for (int b = 0; b < NG; ++b) {
      for (int i = 0; i < NH; ++i) {
        for (int j = 0; j < NH; ++j) {
          flat_gauge[((a * NG + b) * NH + i) * NH + j] = (*gauge)[a][b][i][j];
        }
      }
    }
  }
  std::vector<double> flat_h4(static_cast<size_t>(NH) * NH * NH * NH);
  for (int i = 0; i < NH; ++i) {
    for (int j = 0; j < NH; ++j) {
      for (int k = 0; k < NH; ++k) {
        for (int l = 0; l < NH; ++l) {
          flat_h4[((i * NH + j) * NH + k) * NH + l] = (*h4)[i][j][k][l];
        }
      }
    }
  }

  const std::vector<double> d_gauge22 =
      rotate4_diag(rot_g, NG, rot_h, NH, flat_gauge);
  const std::vector<double> d_higgs4 =
      rotate4_diag(rot_h, NH, rot_h, NH, flat_h4);

  // ------------------------------------------------------------------
  // Fermion sectors
  // ------------------------------------------------------------------
  std::vector<cdouble> flat_quark(static_cast<size_t>(NQ) * NQ * NH);
  for (int i = 0; i < NQ; ++i) {
    for (int j = 0; j < NQ; ++j) {
      for (int k = 0; k < NH; ++k) {
        flat_quark[(i * NQ + j) * NH + k] = (*quark)[i][j][k];
      }
    }
  }
  std::vector<cdouble> flat_lepton(static_cast<size_t>(NL) * NL * NH);
  for (int i = 0; i < NL; ++i) {
    for (int j = 0; j < NL; ++j) {
      for (int k = 0; k < NH; ++k) {
        flat_lepton[(i * NL + j) * NH + k] = (*lepton)[i][j][k];
      }
    }
  }

  const FermionSector quarks =
      build_fermion_sector(flat_quark, NQ, vev, rot_h);
  const FermionSector leptons =
      build_fermion_sector(flat_lepton, NL, vev, rot_h);

  // ------------------------------------------------------------------
  // Gradient
  // ------------------------------------------------------------------
  Eigen::VectorXd grad_mass_basis = Eigen::VectorXd::Zero(NH);
  for (int i = 0; i < NH; ++i) {
    double total = 0.;
    for (int a = 0; a < NG; ++a) {
      const double m_sq = mass_g_sq[a];
      if (m_sq != 0.) {
        total += 1.5 * c_gauge21[(a * NG + a) * NH + i] * m_sq *
                 log_term(m_sq, scale_sq, CW_GAUGE);
      }
    }
    for (int a = 0; a < NH; ++a) {
      const double m_sq = mass_h_sq[a];
      if (m_sq != 0.) {
        total += 0.5 * c_higgs3[(a * NH + a) * NH + i] * m_sq *
                 log_term(m_sq, scale_sq, CW_SCALAR);
      }
    }
    for (int a = 0; a < quarks.n; ++a) {
      const double m_sq = quarks.mass_sq[a];
      if (m_sq != 0.) {
        total -= 3.0 * quarks.c21[(a * quarks.n + a) * NH + i].real() * m_sq *
                 log_term(m_sq, scale_sq, CW_FERMION);
      }
    }
    for (int a = 0; a < leptons.n; ++a) {
      const double m_sq = leptons.mass_sq[a];
      if (m_sq != 0.) {
        total -= 1.0 * leptons.c21[(a * leptons.n + a) * NH + i].real() * m_sq *
                 log_term(m_sq, scale_sq, CW_FERMION);
      }
    }
    grad_mass_basis(i) = total;
  }

  const double eps_loop = 1. / (16. * M_PI * M_PI);
  Eigen::MatrixXd rot_h_mat(NH, NH);
  for (int a = 0; a < NH; ++a) {
    for (int i = 0; i < NH; ++i) {
      rot_h_mat(a, i) = rot_h[a * NH + i];
    }
  }
  const Eigen::VectorXd grad_phys =
      (rot_h_mat.transpose() * grad_mass_basis) * eps_loop;

  // ------------------------------------------------------------------
  // Hessian
  // ------------------------------------------------------------------
  Eigen::MatrixXd storage = Eigen::MatrixXd::Zero(NH, NH);
  for (int i = 0; i < NH; ++i) {
    for (int j = 0; j < NH; ++j) {
      double total = 0.;

      double block = 0.;
      for (int a = 0; a < NG; ++a) {
        for (int b = 0; b < NG; ++b) {
          const double coup1 = c_gauge21[(a * NG + b) * NH + i];
          const double coup2 = c_gauge21[(b * NG + a) * NH + j];
          block += coup1 * coup2 *
                   (fbase(mass_g_sq[a], mass_g_sq[b], scale) - CW_GAUGE + 0.5);
        }
        if (mass_g_sq[a] != 0.) {
          block += d_gauge22[(a * NH + i) * NH + j] * mass_g_sq[a] *
                   log_term(mass_g_sq[a], scale_sq, CW_GAUGE);
        }
      }
      total += 1.5 * block;

      block = 0.;
      for (int a = 0; a < NH; ++a) {
        for (int b = 0; b < NH; ++b) {
          const double coup1 = c_higgs3[(a * NH + b) * NH + i];
          const double coup2 = c_higgs3[(b * NH + a) * NH + j];
          block += coup1 * coup2 *
                   (fbase(mass_h_sq[a], mass_h_sq[b], scale) - CW_SCALAR + 0.5);
        }
        if (mass_h_sq[a] != 0.) {
          block += d_higgs4[(a * NH + i) * NH + j] * mass_h_sq[a] *
                   log_term(mass_h_sq[a], scale_sq, CW_SCALAR);
        }
      }
      total += 0.5 * block;

      for (const auto *sector : {&quarks, &leptons}) {
        const double weight = (sector == &quarks) ? 3.0 : 1.0;
        const int n = sector->n;
        double fermion_block = 0.;
        for (int a = 0; a < n; ++a) {
          for (int b = 0; b < n; ++b) {
            const cdouble coup = sector->c21[(a * n + b) * NH + i] *
                                 sector->c21[(b * n + a) * NH + j];
            fermion_block +=
                coup.real() *
                (fbase(sector->mass_sq[a], sector->mass_sq[b], scale) -
                 CW_FERMION + 0.5);
          }
          if (sector->mass_sq[a] != 0.) {
            fermion_block += sector->d22[(a * NH + i) * NH + j].real() *
                             sector->mass_sq[a] *
                             log_term(sector->mass_sq[a], scale_sq, CW_FERMION);
          }
        }
        total -= weight * fermion_block;
      }

      storage(i, j) = total;
    }
  }
  storage = 0.5 * (storage + storage.transpose()).eval();

  Eigen::MatrixXd hess_phys =
      (rot_h_mat.transpose() * storage * rot_h_mat) * eps_loop;
  hess_phys = 0.5 * (hess_phys + hess_phys.transpose()).eval();

  // ------------------------------------------------------------------
  // Solve the renormalisation conditions
  // ------------------------------------------------------------------
  THDM_tensors::CTMatrix a_raw{};
  THDM_tensors::fill_ct_matrix(p.v1, p.v2, a_raw);
  Eigen::MatrixXd A(NCT, NCT);
  for (int i = 0; i < NCT; ++i) {
    for (int j = 0; j < NCT; ++j) {
      A(i, j) = a_raw[i][j];
    }
  }

  // Same ordering as TL's CONDITION_SPECS: two tadpoles, then the CP-even
  // (1,1), (2,2) and (1,2) entries, one CP-odd and one charged entry.
  Eigen::VectorXd b(NCT);
  b << -grad_phys(4), -grad_phys(6), -hess_phys(4, 4), -hess_phys(6, 6),
      -hess_phys(4, 6), -hess_phys(5, 5), -hess_phys(0, 0);

  const Eigen::VectorXd delta = A.colPivHouseholderQr().solve(b);

  Result result;
  result.dm11_sq = delta(0);
  result.dm22_sq = delta(1);
  result.dm12_sq = delta(2);
  result.dlambda1 = delta(3);
  result.dlambda2 = delta(4);
  result.dlambda3 = delta(5);
  result.dlambda4 = 0.;
  result.dlambda5 = delta(6);
  result.residual = (A * delta - b).cwiseAbs().maxCoeff();
  result.cw_gradient = grad_phys;
  result.cw_hessian = hess_phys;
  return result;
}

} // namespace THDM_ct
} // namespace EffectivePotential

#endif // POTENTIAL_THDM_COUNTERTERMS_HPP_INCLUDED
