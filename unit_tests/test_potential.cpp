#include "catch/catch.hpp"
#include "models/2D_test_model.hpp"
#include "models/Z2_scalar_singlet_model.hpp"

TEST_CASE("Check potential for our two-dimensional test model", "[2DTestModel]") {

  // Construct our model
  EffectivePotential::TwoDimModel model;

  // Points and temperature at which to test it
  Eigen::VectorXd x(model.get_n_scalars());
  x << 0., 0.;

  Eigen::VectorXd y(model.get_n_scalars());
  y << 100., 100.;

  const double T = 100.;
  const double rel_tol = 1.e-2;

  SECTION("Check potential at tree-level")

  CHECK(model.V0(x) == Approx(127840049.999999985).epsilon(rel_tol));
  CHECK(model.V0(y) == Approx(82830862.3471478522).epsilon(rel_tol));

  SECTION("Check potential at one-loop and zero-temperature")

  CHECK(model.V(x, 0.) == Approx(127525505.6793820262).epsilon(rel_tol));
  CHECK(model.V(y, 0.) == Approx(80201817.8519562483).epsilon(rel_tol));

  SECTION("Check potential at one-loop and finite-temperature")

  CHECK(model.V(x, T) == Approx(-226420010.2194074988).epsilon(rel_tol));
  CHECK(model.V(y, T) == Approx(-241315446.8260441422).epsilon(rel_tol));

  SECTION("Check gradient of potential")

  CHECK(model.d2V_dxdt(y, T)(0) == Approx(7088.2212877274).epsilon(rel_tol));
  CHECK(model.d2V_dxdt(y, T)(1) == Approx(6612.627245903).epsilon(rel_tol));

  SECTION("Check Hessian of potential")

  CHECK(model.d2V_dx2(y, T)(0, 1) == Approx(-693.2822294235).epsilon(rel_tol));
  CHECK(model.d2V_dx2(y, T)(0, 0) == Approx(-3077.9320335388).epsilon(rel_tol));
}

TEST_CASE("Check potential for Z2 scalar singlet model", "[Z2ScalarSingletModel]") {

  // Construct our model
  EffectivePotential::Z2ScalarSingletModel model;
  model.set_m_s(90);
  model.set_lambda_hs(0.4);

  // Points and temperature at which to test it
  Eigen::VectorXd x(model.get_n_scalars());
  x << 0., 0.;

  Eigen::VectorXd y(model.get_n_scalars());
  y << 100., 100.;

  const double T = 100.;
  const double rel_tol = 1.e-2;
  const double abs_tol = 1.;

  SECTION("Check potential at zero-temperature")

  CHECK(model.V(x, 0.) == Approx(0.).margin(abs_tol));
  CHECK(model.V(y, 0.) == Approx(-43351047.755965367).epsilon(rel_tol));

  SECTION("Check potential at finite-temperature")

  CHECK(model.V(x, T) == Approx(0.).margin(abs_tol));
  CHECK(model.V(y, T) == Approx(-17959066.345264066).epsilon(rel_tol));

  SECTION("Check gradient of potential")

  CHECK(model.d2V_dxdt(y, T)(0) == Approx(8323.4602203369).epsilon(rel_tol));
  CHECK(model.d2V_dxdt(y, T)(1) == Approx(1833.3332672119).epsilon(rel_tol));

  SECTION("Check Hessian of potential")

  CHECK(model.d2V_dx2(y, T)(0, 1) == Approx(3999.9988250732).epsilon(rel_tol));
  CHECK(model.d2V_dx2(y, T)(0, 0) == Approx(2222.1667480469).epsilon(rel_tol));
}
