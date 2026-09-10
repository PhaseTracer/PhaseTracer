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

#ifndef EFFECTIVEPOTENTIAL_ONE_LOOP_POTENTIAL_HPP_
#define EFFECTIVEPOTENTIAL_ONE_LOOP_POTENTIAL_HPP_

#include <vector>

#include <Eigen/Core>

#include "property.hpp"
#include "potential.hpp"

namespace EffectivePotential {

/** x * log(x) that safely treats x = 0 */
double xlogx(double);

/** Method for daisy corrections */
enum class DaisyMethod { None,
                         ArnoldEspinosa,
                         Parwani };

class OneLoopPotential : public Potential {
public:
 // virtual ~OneLoopPotential() = default;
  virtual double V0(Eigen::VectorXd phi) const = 0;
  /** Functions for squared field dependent masses, depending on:
      a vector of fields and for scalars a xi gauge parameter.  Note the
      latter will be unused for potentials implemented in a fixed gauge. */
  virtual std::vector<double> get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const;
  virtual std::vector<double> get_fermion_masses_sq(Eigen::VectorXd phi) const { return {}; }
  virtual std::vector<double> get_vector_masses_sq(Eigen::VectorXd phi) const { return {}; }
  virtual std::vector<double> get_ghost_masses_sq(Eigen::VectorXd phi, double xi) const { return {}; }
  virtual std::vector<double> get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const { return {}; }
  virtual std::vector<double> get_scalar_thermal_sq(double T) const { return {}; }
  virtual std::vector<double> get_vector_debye_sq(Eigen::VectorXd phi, double T) const { return {}; }
  virtual std::vector<double> get_scalar_dofs() const;
  virtual std::vector<double> get_fermion_dofs() const { return {}; }
  virtual std::vector<double> get_vector_dofs() const { return {}; }
  virtual std::vector<double> get_ghost_dofs() const { return {}; }

  /**
   * Synthesise the fluctuation spectrum (for the one-loop bounce determinant)
   * from the existing field-dependent mass and dof accessors. Scalars, fermions
   * and vectors are included; ghosts are omitted as BubbleDet handles gauge
   * fields directly. Zero-mode assignment is left to the consumer (the field
   * that tunnels carries the Higgs zero mode), so every species is returned
   * with ZeroModeType::None here. Override in a model to refine.
   */
  std::vector<ParticleSpec> get_fluctuation_spectrum(double T) const override;

  /** The Hessian matrix of tree-level potential */
  Eigen::MatrixXd d2V0_dx2(Eigen::VectorXd phi) const;
  /** The derivative of the gradient of potential with respect to temperature */
  Eigen::VectorXd d2V_dxdt(Eigen::VectorXd phi, double T) const;
  /** Finite-temperature effective potential */
  double V(Eigen::VectorXd phi, double T) const;
  /** Zero-temperature one-loop correction */
  virtual double V1(std::vector<double> scalar_masses_sq,
                    std::vector<double> fermion_masses_sq,
                    std::vector<double> vector_masses_sq,
                    std::vector<double> ghost_masses_sq) const;
  virtual double V1(Eigen::VectorXd phi, double T = 0.) const;
  /** Finite-temperature one-loop correction */
  double V1T(std::vector<double> scalar_masses_sq,
             std::vector<double> fermion_masses_sq,
             std::vector<double> vector_masses_sq,
             std::vector<double> ghost_masses_sq, double T) const;
  double V1T(Eigen::VectorXd phi, double T) const;
  /** Daisy corrections to potential */
  double daisy(std::vector<double> scalar_masses_sq,
               std::vector<double> scalar_debye_sq,
               std::vector<double> vector_masses_sq,
               std::vector<double> vector_debye_sq, double T) const;
  double daisy(Eigen::VectorXd phi, double T) const;

  /** Counter-term to potential */
  virtual double counter_term(Eigen::VectorXd phi, double T) const { return 0; }

  /**
   * Contribution to dV1/dphi from one species group.
   * @param sign +1 for scalars and vectors, -1 for fermions and ghosts.
   * @param c    3/2 for scalars, fermions and ghosts; 5/6 for vectors.
   */
  Eigen::VectorXd dV1_term(const std::vector<double> &masses_sq,
                           const std::vector<double> &dofs,
                           const std::vector<Eigen::VectorXd> &d_masses_sq,
                           double sign, double c) const;

  /**
   * Contribution to dV1T/dphi from one species group.
   * @param sign  +1 for scalars, fermions and vectors, -1 for ghosts. Note the
   *              fermionic minus sign is already carried by the tabulated J_F,
   *              which is why fermions take +1 here but -1 in dV1_term.
   * @param fermionic select J_F rather than J_B.
   */
  Eigen::VectorXd dV1T_term(const std::vector<double> &masses_sq,
                            const std::vector<double> &dofs,
                            const std::vector<Eigen::VectorXd> &d_masses_sq,
                            double sign, double T, bool fermionic) const;

  /**
   * Contribution to d(daisy)/dphi from one species group (Arnold-Espinosa).
   * Both the ordinary and the Debye masses, and both their derivatives, are
   * needed because the term is a difference of the two.
   */
  Eigen::VectorXd ddaisy_term(const std::vector<double> &masses_sq,
                              const std::vector<double> &debye_sq,
                              const std::vector<double> &dofs,
                              const std::vector<Eigen::VectorXd> &d_masses_sq,
                              const std::vector<Eigen::VectorXd> &d_debye_sq,
                              double T) const;

  /** High-temperature expansion of potential */
  double VHT(Eigen::VectorXd phi, double T) const;

  /** Tree-level scalar masses */
  std::vector<double> get_tree_scalar_masses_sq(Eigen::VectorXd phi) const;

  /** One-loop scalar masses */
  std::vector<double> get_1l_scalar_masses_sq(Eigen::VectorXd phi, double T) const;

  /** The renormalization scale for the effective potential */
  void set_renormalization_scale(double Q) { renormalization_scale = Q; }
  double get_renormalization_scale() const { return renormalization_scale; }

  /** The gauge parameter \xi */
  void set_xi(double xi_in) { xi = xi_in; }
  double get_xi() const { return xi; }

  /** Treatment of the thermal masses */
  void set_daisy_method(DaisyMethod dm) { daisy_method = dm; }
  DaisyMethod get_daisy_method() const { return daisy_method; }

private:
  /** The renormalization scale for the effective potential */
  double renormalization_scale = 246.;
  /** The gauge parameter \xi */
  double xi = 0.;
  /** Treatment of the thermal masses */
  DaisyMethod daisy_method = DaisyMethod::ArnoldEspinosa;
};

} // namespace EffectivePotential

#endif // EFFECTIVEPOTENTIAL_ONE_LOOP_POTENTIAL_HPP_
