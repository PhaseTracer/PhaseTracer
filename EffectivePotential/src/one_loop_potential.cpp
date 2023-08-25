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
#include <iostream>

#include "one_loop_potential.hpp"
#include "thermal_function.hpp"
#include "pow.hpp"

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
    [scalar_dofs](double e){ return e == scalar_dofs.front(); });
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

/*Eigen::VectorXd OneLoopPotential::d2V_dxdt_old(Eigen::VectorXd phi, double T) const {
  Eigen::VectorXd gradient = Eigen::VectorXd::Zero(phi.size());

  for (int ii = 0; ii < phi.size(); ++ii) {
    Eigen::VectorXd phi_shifted = phi;
    for (int jj = 0; jj < n_h_xy.size(); ++jj) {
      phi_shifted(ii) = phi(ii) + n_h_xy[jj] * h;
      const auto scalar_masses_sq = get_scalar_masses_sq(phi_shifted, xi);
      const auto fermion_masses_sq = get_fermion_masses_sq(phi_shifted);
      const auto vector_masses_sq = get_vector_masses_sq(phi_shifted);
      for (int kk = 0; kk < n_h_xy.size(); ++kk) {
        const double T_shifted = T + n_h_xy[kk] * h;
        const double V1T_ = V1T(scalar_masses_sq, fermion_masses_sq, vector_masses_sq, T_shifted);
        const auto scalar_debye_sq = get_scalar_debye_sq(phi_shifted, xi, T_shifted);
        const auto vector_debye_sq = get_vector_debye_sq(phi_shifted, T_shifted);
        const double daisy_ = daisy(scalar_masses_sq, scalar_debye_sq, vector_masses_sq, vector_debye_sq, T);
        const double counter_term_ = counter_term(phi, T);
        gradient(ii) += (V1T_ + daisy_ + counter_term_) * coeff_xy[jj] * coeff_xy[kk] / square(h);
      }
    }
  }
  return gradient;
}*/

/** Now handles all supported DaisyMethods. */
Eigen::VectorXd OneLoopPotential::d2V_dxdt(Eigen::VectorXd phi, double T) const
{
    Eigen::VectorXd gradient = Eigen::VectorXd::Zero(phi.size());
    Eigen::VectorXd phi_shifted;
    double V1_, V1T_, daisy_, counter_term_;

    for(int i = 0; i < phi.size(); ++i)
    {
        phi_shifted = phi;

        for(int j = 0; j < n_h_xy.size(); ++j)
        {
            phi_shifted(i) = phi(i) + n_h_xy[j]*h;

            const auto scalar_masses_sq = get_scalar_masses_sq(phi_shifted, xi);
            const auto fermion_masses_sq = get_fermion_masses_sq(phi_shifted);
            const auto vector_masses_sq = get_vector_masses_sq(phi_shifted);

            for(int k = 0; k < n_h_xy.size(); ++k)
            {
                const double T_shifted = T + n_h_xy[k]*h;

                counter_term_ = counter_term(phi_shifted, T_shifted);

                switch(daisy_method)
                {
                    case DaisyMethod::None:
                    {
                        V1T_ = V1T(scalar_masses_sq, fermion_masses_sq, vector_masses_sq, T);
                        gradient(i) += (V1T_ + counter_term_)*coeff_xy[j]*coeff_xy[k] / (h*h);
                        break;
                    }
                    case DaisyMethod::ArnoldEspinosa:
                    {
                        const auto scalar_debye_sq = get_scalar_debye_sq(phi_shifted, xi, T_shifted);
                        const auto vector_debye_sq = get_vector_debye_sq(phi_shifted, T_shifted);
                        
                        V1T_ = V1T(scalar_masses_sq, fermion_masses_sq, vector_masses_sq, T);
                        daisy_ = daisy(scalar_masses_sq, scalar_debye_sq, vector_masses_sq, vector_debye_sq,
                            T_shifted);
                        gradient(i) += (V1T_ + daisy_ + counter_term_)*coeff_xy[j]*coeff_xy[k] / (h*h);
                        break;
                    }
                    case DaisyMethod::Parwani:
                    {
                        const auto scalar_debye_sq = get_scalar_debye_sq(phi_shifted, xi, T_shifted);
                        const auto vector_debye_sq = get_vector_debye_sq(phi_shifted, T_shifted);
						
						/*if(T == 0.)
						{
							std::cout << "phi: " << phi_shifted << " T: " << T_shifted << std::endl;
							std::cout << "Scalar debye:";
							for(auto i: scalar_debye_sq)
								std::cout << " " << i;
							std::cout << std::endl;
							
							std::cout << "Vector debye:";
							for(auto i: vector_debye_sq)
								std::cout << " " << i;
							std::cout << std::endl;
						}*/
						
                        V1_= V1(scalar_debye_sq, fermion_masses_sq, vector_debye_sq);
                        V1T_ = V1T(scalar_debye_sq, fermion_masses_sq, vector_debye_sq, T);
						
						/*if(T == 0.)
						{
							std::cout << "V1_: " << V1_ << std::endl;
							std::cout << "V1T_: " << V1T_ << std::endl;
							std::cout << "Gradient += " << (V1_ + V1T_ + counter_term_)*coeff_xy[j]*coeff_xy[k] / (h*h) << std::endl;
						}*/
						
                        gradient(i) += (V1_ + V1T_ + counter_term_)*coeff_xy[j]*coeff_xy[k] / (h*h);
                        break;
                    }
                    default:
                        throw std::runtime_error("unknown daisy method");
                }
            }
        }
    }

    return gradient;
}

std::vector<double> OneLoopPotential::get_scalar_dofs() const {
  const std::vector<double> dof(get_n_scalars(), 1.);
  return dof;
}

double OneLoopPotential::V1(std::vector<double> scalar_masses_sq,
                            std::vector<double> fermion_masses_sq,
                            std::vector<double> vector_masses_sq) const {
  double correction = 0;

  const auto scalar_dofs = get_scalar_dofs();
  const auto fermion_dofs = get_fermion_dofs();
  const auto vector_dofs = get_vector_dofs();
  
  /*if(true)
	  {
		  std::cout << "Scalars:";
		  for(auto i: scalar_masses_sq)
			  std::cout << " " << i;
		  std::cout << std::endl;
		  
		  std::cout << "Fermions:";
		  for(auto i: fermion_masses_sq)
			  std::cout << " " << i;
		  std::cout << std::endl;
		  
		  std::cout << "Bosons:";
		  for(auto i: vector_masses_sq)
			  std::cout << " " << i;
		  std::cout << std::endl;
	  }*/

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

  // gauge dependent vector correction
  if (xi != 0.) {
    for (size_t i = 0; i < vector_masses_sq.size(); ++i) {
      const double x = xi * vector_masses_sq[i] / square(renormalization_scale);
      correction -= vector_dofs[i] / 3. * xi * vector_masses_sq[i] *
        (square(renormalization_scale) * xlogx(x) - xi * vector_masses_sq[i] * 3. / 2.);
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
                get_vector_masses_sq(phi));
    case DaisyMethod::Parwani:
      return V1(get_scalar_debye_sq(phi, xi, T),
                get_fermion_masses_sq(phi),
                get_vector_debye_sq(phi, T));
    default:
      throw std::runtime_error("unknown daisy method");
  }
}

double OneLoopPotential::V1T(std::vector<double> scalar_masses_sq,
                             std::vector<double> fermion_masses_sq,
                             std::vector<double> vector_masses_sq, double T) const {
  double correction = 0;
  
  if(T == 0.)
	  return 0.;

  const auto scalar_dofs = get_scalar_dofs();
  const auto fermion_dofs = get_fermion_dofs();
  const auto vector_dofs = get_vector_dofs();

  if (scalar_dofs.size() != scalar_masses_sq.size()) {
    throw std::runtime_error("Scalar dofs and masses do not match");
  }

  if (fermion_dofs.size() != fermion_masses_sq.size()) {
    throw std::runtime_error("Fermion dofs and masses do not match");
  }

  if (vector_dofs.size() != vector_masses_sq.size()) {
    throw std::runtime_error("Vector dofs and masses do not match");
  }
  
  if(T == 0.)
    std::cout << "Correction: " << correction << " initially." << std::endl;

  // scalar correction
  for (size_t i = 0; i < scalar_masses_sq.size(); ++i) {
    correction += scalar_dofs[i] * J_B(scalar_masses_sq[i] / square(T));
  }
  
  if(T == 0.)
    std::cout << "Correction: " << correction << " after scalars." << std::endl;

  // fermion correction
  for (size_t i = 0; i < fermion_masses_sq.size(); ++i) {
    correction += fermion_dofs[i] * J_F(fermion_masses_sq[i] / square(T));
  }
  
  if(T == 0.)
    std::cout << "Correction: " << correction << " after fermions." << std::endl;

  // vector correction
  for (size_t i = 0; i < vector_masses_sq.size(); ++i) {
    correction += vector_dofs[i] * J_B(vector_masses_sq[i] / square(T));
  }
  
  if(T == 0.)
    std::cout << "Correction: " << correction << " after vectors." << std::endl;

  // gauge dependent vector correction
  if (xi != 0.) {
    for (size_t i = 0; i < vector_masses_sq.size(); ++i) {
      correction -= vector_dofs[i] / 3. * J_B(xi * vector_masses_sq[i] / square(T));
    }
  }
  
  if(T == 0.)
    std::cout << "Correction: " << correction << " after gauge-dependent vector contributions." << std::endl;

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
                 get_vector_masses_sq(phi), T);
    case DaisyMethod::Parwani:
      return V1T(get_scalar_debye_sq(phi, xi, T),
                 get_fermion_masses_sq(phi),
                 get_vector_debye_sq(phi, T), T);
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
    correction += scalar_dofs[i] * (std::pow(std::max(0., scalar_debye_sq[i]), 1.5)
                    - std::pow(std::max(0., scalar_masses_sq[i]), 1.5));
  }

  for (size_t i = 0; i < vector_debye_sq.size(); ++i) {
    correction += vector_dofs[i] / 3. * (std::pow(std::max(0., vector_debye_sq[i]), 1.5)
                    - std::pow(std::max(0., vector_masses_sq[i]), 1.5));
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
  const auto scalar_masses_sq = get_scalar_masses_sq(phi, xi);
  const auto fermion_masses_sq = get_fermion_masses_sq(phi);
  const auto vector_masses_sq = get_vector_masses_sq(phi);

  if (T > 0) {
    const auto scalar_debye_sq = get_scalar_debye_sq(phi, xi, T);
    const auto vector_debye_sq = get_vector_debye_sq(phi, T);
    switch (daisy_method) {
      case DaisyMethod::None:
        return V0(phi)
               + V1(scalar_masses_sq, fermion_masses_sq, vector_masses_sq)
               + V1T(scalar_masses_sq, fermion_masses_sq, vector_masses_sq, T)
               + counter_term(phi, T);
      case DaisyMethod::ArnoldEspinosa:
        return V0(phi)
               + daisy(scalar_masses_sq, scalar_debye_sq, vector_masses_sq, vector_debye_sq, T)
               + V1(scalar_masses_sq, fermion_masses_sq, vector_masses_sq)
               + V1T(scalar_masses_sq, fermion_masses_sq, vector_masses_sq, T)
               + counter_term(phi, T);
      case DaisyMethod::Parwani:
        return V0(phi)
               + V1(scalar_debye_sq, fermion_masses_sq, vector_debye_sq)
               + V1T(scalar_debye_sq, fermion_masses_sq, vector_debye_sq, T)
               + counter_term(phi, T);
      default:
        throw std::runtime_error("unknown daisy method");
    }
  } else {
    return V0(phi) + V1(scalar_masses_sq, fermion_masses_sq, vector_masses_sq);
  }
}

}  // namespace EffectivePotential
