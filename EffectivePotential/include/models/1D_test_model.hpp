/**
  A one dimensional example model with known analytic solution.
*/

#include "potential.hpp"
#include "property.hpp"
#include "pow.hpp"

namespace EffectivePotential {

class OneDimModel : public Potential {
 public:
  double V(Eigen::VectorXd x, double T) const override {
    return (c * square(T) - m2) * square(x[0]) + kappa * cube(x[0]) + lambda * pow_4(x[0]);
  }
  size_t get_n_scalars() const override { return 1; }
  bool forbidden(Eigen::VectorXd x) const override { return x[0] < -0.1; }
  
  Eigen::VectorXd d2V_dxdt(Eigen::VectorXd phi, double T) const override {
    return 4. * c * T * phi;
  }

  Eigen::MatrixXd d2V_dx2(Eigen::VectorXd phi, double T) const  override {
    Eigen::MatrixXd d(1, 1);
    d(0, 0) = 2. * (c * square(T) - m2) + 6. * kappa * phi[0] + 12. * lambda * square(phi[0]);
    return d;
  }

  double get_TC_from_expression() const {
    return 0.5 * std::sqrt(square(kappa) / (c * lambda) + 4. * m2 / c);
  }

  double get_true_vacuum_from_expression() const {
    return -0.5 * kappa / lambda;
  }

  double get_false_vacuum_from_expression() const {
    return 0.;
  }

 private:
  PROPERTY(double, m2, 100.)
  PROPERTY(double, kappa, -10.)
  PROPERTY(double, lambda, 0.1)
  PROPERTY(double, c, 0.1)
};

}  // namespace EffectivePotential
