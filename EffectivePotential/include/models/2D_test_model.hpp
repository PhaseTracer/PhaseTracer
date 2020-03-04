#ifndef POTENTIAL_TEST_MODEL_HPP_INCLUDED
#define POTENTIAL_TEST_MODEL_HPP_INCLUDED

/**
   Make dummy model for testing.

  This is `model1` from the examples in `CosmoTransitions`.
*/

#include <vector>

#include "one_loop_potential.hpp"
#include "pow.hpp"


namespace EffectivePotential {

class TwoDimModel : public OneLoopPotential {
 public:

  double V0(Eigen::VectorXd phi) const override {
    return 0.25 * l1 * square(square(phi[0]) - square(v))
           + 0.25 * l2 * square(square(phi[1]) - square(v))
           - square(mu) * phi[0] * phi[1];
  }

  std::vector<double> get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const override {
    const double a = l1 * (3. * square(phi[0]) - square(v));
    const double b = l2 * (3. * square(phi[1]) - square(v));
    const double A = 0.5 * (a + b);
    const double B = std::sqrt(0.25 * square(a - b) + pow_4(mu));
    const double mb_sq = y1 * (square(phi[0]) + square(phi[1])) + y2 * phi[0] * phi[1];
    return {A + B, A - B,mb_sq};
  }

  std::vector<double> get_scalar_dofs() const override { return {1., 1., 30.}; }
  size_t get_n_scalars() const override { return 2; }

  std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
    return {-phi};
  };

 private:
  const double v = 246.;
  const double m1 = 120.;
  const double m2 = 50.;
  const double mu = 25.;
  const double l1 = 0.5 * square(m1 / v);
  const double l2 = 0.5 * square(m2 / v);
  const double y1 = 0.1;
  const double y2 = 0.15;
};

}  // namespace EffectivePotential

#endif
