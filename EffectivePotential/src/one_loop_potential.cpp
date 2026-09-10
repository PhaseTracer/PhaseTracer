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

#include <cmath>
#include <algorithm>
#include <Eigen/Eigenvalues>

#include "one_loop_potential.hpp"
#include "thermal_function.hpp"
#include "pow.hpp"
#include "logger.hpp"

namespace EffectivePotential {

double xlogx(double x) {
  const double abs_x = std::abs(x);
  if (abs_x <= std::numeric_limits<double>::min()) {
    return 0.;
  } else {
    return x * std::log(abs_x);
  }
}

std::vector<double> OneLoopPotential::get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const {
  if (xi != 0.) {
    throw std::runtime_error("Default implementation of scalar masses assumes xi = 0.");
  }
  const auto scalar_dofs = get_scalar_dofs();
  const bool identical = std::all_of(scalar_dofs.begin(), scalar_dofs.end(),
                                     [scalar_dofs](double e) { return e == scalar_dofs.front(); });
  if (!identical) {
    throw std::runtime_error("Default implementation of scalar masses assumes scalar dof are equal.");
  }
  const Eigen::MatrixXd hessian = d2V0_dx2(phi);
  const Eigen::VectorXd m_sq = hessian.eigenvalues().real();
  std::vector<double> m_sq_vector(m_sq.data(), m_sq.data() + m_sq.rows() * m_sq.cols());
  std::sort(m_sq_vector.begin(), m_sq_vector.end());
  return m_sq_vector;
}

std::vector<double> OneLoopPotential::get_1l_scalar_masses_sq(Eigen::VectorXd phi, double T) const {
  const Eigen::MatrixXd hessian = d2V_dx2(phi, T);
  const Eigen::VectorXd m_sq = hessian.eigenvalues().real();
  std::vector<double> m_sq_vector(m_sq.data(), m_sq.data() + m_sq.rows() * m_sq.cols());
  std::sort(m_sq_vector.begin(), m_sq_vector.end());
  return m_sq_vector;
}

std::vector<double> OneLoopPotential::get_tree_scalar_masses_sq(Eigen::VectorXd phi) const {
  const Eigen::MatrixXd hessian = d2V0_dx2(phi);
  const Eigen::VectorXd m_sq = hessian.eigenvalues().real();
  std::vector<double> m_sq_vector(m_sq.data(), m_sq.data() + m_sq.rows() * m_sq.cols());
  std::sort(m_sq_vector.begin(), m_sq_vector.end());
  return m_sq_vector;
}

Eigen::MatrixXd OneLoopPotential::d2V0_dx2(Eigen::VectorXd phi) const {
  Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(phi.size(), phi.size());

  // diagonal elements
  for (int ii = 0; ii < phi.size(); ++ii) {
    Eigen::VectorXd phi_shifted = phi;
    for (int jj = 0; jj < n_h_xx.size(); ++jj) {
      phi_shifted(ii) = phi(ii) + n_h_xx[jj] * h;
      hessian(ii, ii) += V0(phi_shifted) * coeff_xx[jj] / square(h);
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
          hessian(ii, jj) += V0(phi_shifted) * coeff_xy[kk] * coeff_xy[ll] / square(h);
        }
      }
      hessian(jj, ii) = hessian(ii, jj);
    }
  }
  return hessian;
}

Eigen::VectorXd OneLoopPotential::d2V_dxdt(Eigen::VectorXd phi, double T) const {
  Eigen::VectorXd gradient = Eigen::VectorXd::Zero(phi.size());

  for (int ii = 0; ii < phi.size(); ++ii) {
    Eigen::VectorXd phi_shifted = phi;
    for (int jj = 0; jj < n_h_xy.size(); ++jj) {
      phi_shifted(ii) = phi(ii) + n_h_xy[jj] * h;
      for (int kk = 0; kk < n_h_xy.size(); ++kk) {
        const double T_shifted = T + n_h_xy[kk] * h;
        gradient(ii) += V(phi_shifted, T_shifted) * coeff_xy[jj] * coeff_xy[kk] / square(h);
      }
    }
  }
  return gradient;
}

std::vector<double> OneLoopPotential::get_scalar_dofs() const {
  return std::vector<double>(get_n_scalars(), 1.);
}

std::vector<ParticleSpec> OneLoopPotential::get_fluctuation_spectrum(double T) const {
  std::vector<ParticleSpec> spectrum;
  const double xi_ = get_xi();
  // The default uses zero-temperature field-dependent masses; models with an
  // explicit thermal/3d spectrum can use T when overriding.
  (void)T;

  // Append one ParticleSpec per dof entry, with a mass closure that recomputes
  // the corresponding field-dependent mass-squared on demand. Any T or xi
  // dependence is captured here so the returned closures take only phi.
  const auto append = [&spectrum](const std::string &tag, double spin,
                                  FieldStatistics statistics,
                                  const std::vector<double> &dofs,
                                  std::function<std::vector<double>(const Eigen::VectorXd &)> masses) {
    for (size_t i = 0; i < dofs.size(); ++i) {
      ParticleSpec s;
      s.name = tag + "_" + std::to_string(i);
      s.spin = spin;
      s.statistics = statistics;
      s.dof = dofs[i];
      s.zero_mode = ZeroModeType::None;
      s.mass_sq = [masses, i](const Eigen::VectorXd &phi) { return masses(phi).at(i); };
      spectrum.push_back(std::move(s));
    }
  };

  append("scalar", 0.0, FieldStatistics::Boson, get_scalar_dofs(),
         [this, xi_](const Eigen::VectorXd &phi) { return get_scalar_masses_sq(phi, xi_); });
  append("fermion", 0.5, FieldStatistics::Fermion, get_fermion_dofs(),
         [this](const Eigen::VectorXd &phi) { return get_fermion_masses_sq(phi); });
  append("vector", 1.0, FieldStatistics::Boson, get_vector_dofs(),
         [this](const Eigen::VectorXd &phi) { return get_vector_masses_sq(phi); });

  return spectrum;
}

namespace {

/** Number of fields a species-gradient list refers to, or 0 if it is empty. */
int gradient_dim(const std::vector<Eigen::VectorXd> &d_masses_sq) {
  return d_masses_sq.empty() ? 0 : static_cast<int>(d_masses_sq.front().size());
}

void check_gradient_sizes(const std::vector<double> &masses_sq,
                          const std::vector<double> &dofs,
                          const std::vector<Eigen::VectorXd> &d_masses_sq,
                          const char *what) {
  if (dofs.size() != masses_sq.size() || d_masses_sq.size() != masses_sq.size()) {
    throw std::runtime_error(std::string(what) +
                             " masses, dofs and mass gradients do not match");
  }
}

} // namespace

Eigen::VectorXd OneLoopPotential::dV1_term(
    const std::vector<double> &masses_sq, const std::vector<double> &dofs,
    const std::vector<Eigen::VectorXd> &d_masses_sq, double sign, double c) const {
  check_gradient_sizes(masses_sq, dofs, d_masses_sq, "dV1_term");
  Eigen::VectorXd grad = Eigen::VectorXd::Zero(gradient_dim(d_masses_sq));

  const double q_sq = square(renormalization_scale);
  for (size_t i = 0; i < masses_sq.size(); ++i) {
    const double m_sq = masses_sq[i];
    // d/dm^2 [ m^4 (log|m^2/Q^2| - c) ] = m^2 (2 log|m^2/Q^2| + 1 - 2c).
    // The log diverges as m^2 -> 0 but is multiplied by m^2, so the product
    // vanishes; cut off at the same point xlogx does, to match V1 exactly.
    if (std::abs(m_sq) <= std::numeric_limits<double>::min()) {
      continue;
    }
    const double factor = m_sq * (2. * std::log(std::abs(m_sq / q_sq)) + 1. - 2. * c);
    grad += (sign * dofs[i] * factor) * d_masses_sq[i];
  }
  return grad / (64. * M_PI * M_PI);
}

Eigen::VectorXd OneLoopPotential::dV1T_term(
    const std::vector<double> &masses_sq, const std::vector<double> &dofs,
    const std::vector<Eigen::VectorXd> &d_masses_sq, double sign, double T,
    bool fermionic) const {
  check_gradient_sizes(masses_sq, dofs, d_masses_sq, "dV1T_term");
  Eigen::VectorXd grad = Eigen::VectorXd::Zero(gradient_dim(d_masses_sq));

  if (T == 0.) {
    return grad;
  }

  const double T_sq = square(T);
  for (size_t i = 0; i < masses_sq.size(); ++i) {
    const double x = masses_sq[i] / T_sq;
    const double dJ = fermionic ? J_F_diff(x) : J_B_diff(x);
    grad += (sign * dofs[i] * dJ) * d_masses_sq[i];
  }
  // V1T = T^4/(2 pi^2) sum n J(m^2/T^2); the 1/T^2 from the argument leaves T^2.
  return grad * T_sq / (2. * square(M_PI));
}

Eigen::VectorXd OneLoopPotential::ddaisy_term(
    const std::vector<double> &masses_sq, const std::vector<double> &debye_sq,
    const std::vector<double> &dofs,
    const std::vector<Eigen::VectorXd> &d_masses_sq,
    const std::vector<Eigen::VectorXd> &d_debye_sq, double T) const {
  check_gradient_sizes(masses_sq, dofs, d_masses_sq, "ddaisy_term (ordinary)");
  check_gradient_sizes(debye_sq, dofs, d_debye_sq, "ddaisy_term (Debye)");
  Eigen::VectorXd grad = Eigen::VectorXd::Zero(gradient_dim(d_masses_sq));

  // d/dm^2 [ max(0, m^2)^{3/2} ] = 1.5 sqrt(m^2) for m^2 > 0, else 0. Continuous
  // at the clamp, but with an infinite second derivative there.
  for (size_t i = 0; i < masses_sq.size(); ++i) {
    if (debye_sq[i] > 0.) {
      grad += (1.5 * dofs[i] * std::sqrt(debye_sq[i])) * d_debye_sq[i];
    }
    if (masses_sq[i] > 0.) {
      grad -= (1.5 * dofs[i] * std::sqrt(masses_sq[i])) * d_masses_sq[i];
    }
  }
  return grad * -T / (12. * M_PI);
}

double OneLoopPotential::V1(std::vector<double> scalar_masses_sq,
                            std::vector<double> fermion_masses_sq,
                            std::vector<double> vector_masses_sq,
                            std::vector<double> ghost_masses_sq) const {
  double correction = 0;

  const auto scalar_dofs = get_scalar_dofs();
  const auto fermion_dofs = get_fermion_dofs();
  const auto vector_dofs = get_vector_dofs();
  const auto ghost_dofs = get_ghost_dofs();

  if (scalar_dofs.size() != scalar_masses_sq.size()) {
    throw std::runtime_error("Scalar dofs and masses do not match");
  }

  if (fermion_dofs.size() != fermion_masses_sq.size()) {
    throw std::runtime_error("Fermion dofs and masses do not match");
  }

  if (vector_dofs.size() != vector_masses_sq.size()) {
    throw std::runtime_error("Vector dofs and masses do not match");
  }

  // scalar correction
  for (size_t i = 0; i < scalar_masses_sq.size(); ++i) {
    const double x = scalar_masses_sq[i] / square(renormalization_scale);
    correction += scalar_dofs[i] * scalar_masses_sq[i] *
                  (square(renormalization_scale) * xlogx(x) - scalar_masses_sq[i] * 3. / 2.);
  }

  // fermion correction
  for (size_t i = 0; i < fermion_masses_sq.size(); ++i) {
    const double x = fermion_masses_sq[i] / square(renormalization_scale);
    correction -= fermion_dofs[i] * fermion_masses_sq[i] *
                  (square(renormalization_scale) * xlogx(x) - fermion_masses_sq[i] * 3. / 2.);
  }

  // vector correction
  for (size_t i = 0; i < vector_masses_sq.size(); ++i) {
    const double x = vector_masses_sq[i] / square(renormalization_scale);
    correction += vector_dofs[i] * vector_masses_sq[i] *
                  (square(renormalization_scale) * xlogx(x) - vector_masses_sq[i] * 5. / 6.);
  }

  // ghost correction
  if (xi != 0.) {
    for (size_t i = 0; i < ghost_masses_sq.size(); ++i) {
      const double x = ghost_masses_sq[i] / square(renormalization_scale);
      correction -= ghost_dofs[i] * ghost_masses_sq[i] *
                    (square(renormalization_scale) * xlogx(x) - ghost_masses_sq[i] * 3. / 2.);
    }
  }

  return correction / (64. * M_PI * M_PI);
}

double OneLoopPotential::V1(Eigen::VectorXd phi, double T) const {
  if (T != 0 && daisy_method != DaisyMethod::Parwani) {
    throw std::runtime_error("V1 dependent on T only with Parwani daisy method");
  }

  switch (daisy_method) {
  // Arnold-Espinosa method does not alter one-loop zero-temperature potential
  // hence here it is the same as no daisy corrections
  case DaisyMethod::None:
  case DaisyMethod::ArnoldEspinosa:
    return V1(get_scalar_masses_sq(phi, xi),
              get_fermion_masses_sq(phi),
              get_vector_masses_sq(phi),
              get_ghost_masses_sq(phi, xi));
  case DaisyMethod::Parwani:
    return V1(get_scalar_debye_sq(phi, xi, T),
              get_fermion_masses_sq(phi),
              get_vector_debye_sq(phi, T),
              get_ghost_masses_sq(phi, xi));
  default:
    throw std::runtime_error("unknown daisy method");
  }
}

double OneLoopPotential::V1T(std::vector<double> scalar_masses_sq,
                             std::vector<double> fermion_masses_sq,
                             std::vector<double> vector_masses_sq,
                             std::vector<double> ghost_masses_sq, double T) const {
  double correction = 0;

  if (T == 0.)
    return 0.;

  const auto scalar_dofs = get_scalar_dofs();
  const auto fermion_dofs = get_fermion_dofs();
  const auto vector_dofs = get_vector_dofs();
  const auto ghost_dofs = get_ghost_dofs();

  if (scalar_dofs.size() != scalar_masses_sq.size()) {
    throw std::runtime_error("Scalar dofs and masses do not match");
  }

  if (fermion_dofs.size() != fermion_masses_sq.size()) {
    throw std::runtime_error("Fermion dofs and masses do not match");
  }

  if (vector_dofs.size() != vector_masses_sq.size()) {
    throw std::runtime_error("Vector dofs and masses do not match");
  }

  if (T == 0.)
    LOG(debug) << "Correction: " << correction << " initially";

  // scalar correction
  for (size_t i = 0; i < scalar_masses_sq.size(); ++i) {
    correction += scalar_dofs[i] * J_B(scalar_masses_sq[i] / square(T));
  }

  if (T == 0.)
    LOG(debug) << "Correction: " << correction << " after scalars";

  // fermion correction
  for (size_t i = 0; i < fermion_masses_sq.size(); ++i) {
    correction += fermion_dofs[i] * J_F(fermion_masses_sq[i] / square(T));
  }

  if (T == 0.)
    LOG(debug) << "Correction: " << correction << " after fermions";

  // vector correction
  for (size_t i = 0; i < vector_masses_sq.size(); ++i) {
    correction += vector_dofs[i] * J_B(vector_masses_sq[i] / square(T));
  }

  if (T == 0.)
    LOG(debug) << "Correction: " << correction << " after vectors";

  // ghost correction
  if (xi != 0.) {
    for (size_t i = 0; i < ghost_masses_sq.size(); ++i) {
      correction -= ghost_dofs[i] * J_B(ghost_masses_sq[i] / square(T));
    }
  }

  if (T == 0.)
    LOG(debug) << "Correction: " << correction << " after gauge-dependent vector contributions";

  return correction * pow_4(T) / (2. * square(M_PI));
}

double OneLoopPotential::V1T(Eigen::VectorXd phi, double T) const {

  switch (daisy_method) {
  // Arnold-Espinosa method does not alter one-loop finite-temperature potential
  // hence here it is the same as no daisy corrections
  case DaisyMethod::None:
  case DaisyMethod::ArnoldEspinosa:
    return V1T(get_scalar_masses_sq(phi, xi),
               get_fermion_masses_sq(phi),
               get_vector_masses_sq(phi),
               get_ghost_masses_sq(phi, xi), T);
  case DaisyMethod::Parwani:
    return V1T(get_scalar_debye_sq(phi, xi, T),
               get_fermion_masses_sq(phi),
               get_vector_debye_sq(phi, T),
               get_ghost_masses_sq(phi, xi), T);
  default:
    throw std::runtime_error("unknown daisy method");
  }
}

double OneLoopPotential::VHT(Eigen::VectorXd phi, double T) const {
  const auto thermal_sq = get_scalar_thermal_sq(T);
  double V = V0(phi);
  for (int i = 0; i < thermal_sq.size(); i++) {
    V += 0.5 * thermal_sq[i] * square(phi(i));
  }
  return V;
}

double OneLoopPotential::daisy(std::vector<double> scalar_masses_sq,
                               std::vector<double> scalar_debye_sq,
                               std::vector<double> vector_masses_sq,
                               std::vector<double> vector_debye_sq, double T) const {
  if (daisy_method != DaisyMethod::ArnoldEspinosa) {
    throw std::runtime_error("Explicit daisy term only present in Arnold-Espinosa");
  }

  double correction = 0.;

  const auto scalar_dofs = get_scalar_dofs();
  const auto vector_dofs = get_vector_dofs();

  for (size_t i = 0; i < scalar_debye_sq.size(); ++i) {
    correction += scalar_dofs[i] * (std::pow(std::max(0., scalar_debye_sq[i]), 1.5) - std::pow(std::max(0., scalar_masses_sq[i]), 1.5));
  }

  // hack - i know only first 3 are longitudinal in xSM model
  for (size_t i = 0; i < vector_masses_sq.size(); ++i) {
    correction += vector_dofs[i] * (std::pow(std::max(0., vector_debye_sq[i]), 1.5) - std::pow(std::max(0., vector_masses_sq[i]), 1.5));
  }

  return correction * -T / (12. * M_PI);
}

double OneLoopPotential::daisy(Eigen::VectorXd phi, double T) const {
  const auto scalar_masses_sq = get_scalar_masses_sq(phi, xi);
  const auto vector_masses_sq = get_vector_masses_sq(phi);
  const auto scalar_debye_sq = get_scalar_debye_sq(phi, xi, T);
  const auto vector_debye_sq = get_vector_debye_sq(phi, T);
  return daisy(scalar_masses_sq, scalar_debye_sq, vector_masses_sq, vector_debye_sq, T);
}

double OneLoopPotential::V(Eigen::VectorXd phi, double T) const {

  const auto fermion_masses_sq = get_fermion_masses_sq(phi);

  const auto ghost_masses_sq = (xi != 0.) ? get_ghost_masses_sq(phi, xi) : std::vector<double>{};

  if (T > 0) {
    switch (daisy_method) {
    case DaisyMethod::None: {
      const auto scalar_masses_sq = get_scalar_masses_sq(phi, xi);
      const auto vector_masses_sq = get_vector_masses_sq(phi);
      return V0(phi) + V1(scalar_masses_sq, fermion_masses_sq, vector_masses_sq, ghost_masses_sq) + V1T(scalar_masses_sq, fermion_masses_sq, vector_masses_sq, ghost_masses_sq, T) + counter_term(phi, T);
    }
    case DaisyMethod::ArnoldEspinosa: {
      const auto scalar_masses_sq = get_scalar_masses_sq(phi, xi);
      const auto vector_masses_sq = get_vector_masses_sq(phi);
      const auto scalar_debye_sq = get_scalar_debye_sq(phi, xi, T);
      const auto vector_debye_sq = get_vector_debye_sq(phi, T);
      return V0(phi) + daisy(scalar_masses_sq, scalar_debye_sq, vector_masses_sq, vector_debye_sq, T) + V1(scalar_masses_sq, fermion_masses_sq, vector_masses_sq, ghost_masses_sq) + V1T(scalar_masses_sq, fermion_masses_sq, vector_masses_sq, ghost_masses_sq, T) + counter_term(phi, T);
    }
    case DaisyMethod::Parwani: {
      const auto scalar_debye_sq = get_scalar_debye_sq(phi, xi, T);
      const auto vector_debye_sq = get_vector_debye_sq(phi, T);
      return V0(phi) + V1(scalar_debye_sq, fermion_masses_sq, vector_debye_sq, ghost_masses_sq) + V1T(scalar_debye_sq, fermion_masses_sq, vector_debye_sq, ghost_masses_sq, T) + counter_term(phi, T);
    }
    default:
      throw std::runtime_error("unknown daisy method");
    }
  } else {
    const auto scalar_masses_sq = get_scalar_masses_sq(phi, xi);
    const auto vector_masses_sq = get_vector_masses_sq(phi);
    return V0(phi) + V1(scalar_masses_sq, fermion_masses_sq, vector_masses_sq, ghost_masses_sq) + counter_term(phi, T);
  }
}

} // namespace EffectivePotential
