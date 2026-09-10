#include <cmath>
#include <algorithm>
#include <string>

#include "catch/catch.hpp"
#include "models/xSM_MSbar.hpp"
#include "logger.hpp"

/*
  The analytic gradient of the xSM_MSbar one-loop potential must agree with a
  finite difference of the same potential. The finite difference is the baseline
  here, not the truth: it carries O(h^4) truncation error and a roundoff floor of
  roughly eps*|V|/h, so the tolerances below reflect its accuracy, not the
  analytic gradient's.
*/

namespace {

EffectivePotential::xSM_MSbar make_model() {
  auto model = EffectivePotential::xSM_MSbar::from_tadpoles(
      1.05, 1.00, 125., 100., 1.0, false, false, true, false, {});
  // Fourth-order stencil, for the most accurate finite-difference baseline.
  model.set_h_4(true);
  return model;
}

// Below this the gradient is numerically indistinguishable from zero and the
// finite difference is dominated by its own roundoff.
constexpr double gradient_noise_floor = 1.e-3;

} // namespace

TEST_CASE("Analytic gradient matches finite differences", "[xSM_MSbar][gradient]") {
  LOGGER(fatal);

  auto model = make_model();

  const std::vector<double> hs = {0., 0.7, 5., 25., 60., 123.4, 180., 246., 330., 400.};
  const std::vector<double> ss = {0., 0.7, 5., 30., 80., 150., 240., 330.};
  const std::vector<double> Ts = {0., 2., 15., 40., 70., 88.3, 110., 150.};

  const std::vector<std::pair<std::string, EffectivePotential::DaisyMethod>> methods = {
      {"None", EffectivePotential::DaisyMethod::None},
      {"ArnoldEspinosa", EffectivePotential::DaisyMethod::ArnoldEspinosa},
      {"Parwani", EffectivePotential::DaisyMethod::Parwani}};

  for (const auto &method : methods) {
    model.set_daisy_method(method.second);

    for (double h : hs) {
      for (double s : ss) {
        for (double T : Ts) {
          Eigen::VectorXd phi(2);
          phi << h, s;

          model.set_use_analytic_gradient(true);
          const Eigen::VectorXd analytic = model.dV_dx(phi, T);
          // Explicitly qualified so this is the base-class finite difference
          // rather than the override under test.
          const Eigen::VectorXd numerical =
              model.EffectivePotential::Potential::dV_dx(phi, T);

          INFO("daisy = " << method.first << ", h = " << h << ", s = " << s
                          << ", T = " << T);

          if (numerical.norm() < gradient_noise_floor) {
            // Where the true gradient vanishes, only require that the analytic
            // one vanishes too; the finite difference here is pure roundoff.
            CHECK(analytic.norm() < gradient_noise_floor);
          } else {
            CHECK((analytic - numerical).norm() / numerical.norm() ==
                  Approx(0.).margin(1.e-4));
          }
        }
      }
    }
  }
}

TEST_CASE("Analytic and finite-difference gradients are selectable", "[xSM_MSbar][gradient]") {
  LOGGER(fatal);

  auto model = make_model();
  model.set_daisy_method(EffectivePotential::DaisyMethod::Parwani);

  Eigen::VectorXd phi(2);
  phi << 123.4, 80.;
  const double T = 70.;

  model.set_use_analytic_gradient(false);
  const Eigen::VectorXd off = model.dV_dx(phi, T);
  model.set_use_analytic_gradient(true);
  const Eigen::VectorXd on = model.dV_dx(phi, T);

  CHECK(model.get_use_analytic_gradient());
  // Same quantity by two routes, so they must agree to finite-difference
  // accuracy while not being the identical code path.
  CHECK((on - off).norm() / off.norm() == Approx(0.).margin(1.e-4));
}

TEST_CASE("Covariant gauge falls back to finite differences", "[xSM_MSbar][gradient]") {
  LOGGER(fatal);

  // The covariant-gauge branch clamps its discriminant, so the gradient is
  // kinked and the analytic path deliberately does not cover it.
  auto model = EffectivePotential::xSM_MSbar::from_tadpoles(
      1.05, 1.00, 125., 100., 1.0, true, false, true, false, {});
  model.set_h_4(true);
  model.set_daisy_method(EffectivePotential::DaisyMethod::Parwani);

  Eigen::VectorXd phi(2);
  phi << 123.4, 80.;
  const double T = 70.;

  const Eigen::VectorXd from_override = model.dV_dx(phi, T);
  const Eigen::VectorXd from_base = model.EffectivePotential::Potential::dV_dx(phi, T);

  // Identical because the override delegates to the base implementation.
  CHECK(from_override(0) == from_base(0));
  CHECK(from_override(1) == from_base(1));
}
