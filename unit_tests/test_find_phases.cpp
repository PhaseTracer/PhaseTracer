#include "catch/catch.hpp"
#include "models/2D_test_model.hpp"
#include "models/Z2_scalar_singlet_model.hpp"
#include "phase_finder.hpp"
#include "logger.hpp"


TEST_CASE("Compute phases for a two-dimensional model", "[2DTestModel]") {

  LOGGER(fatal);

  // Construct our model
  EffectivePotential::TwoDimModel model;

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
  pf.set_seed(1);
  pf.find_phases();
  auto phases = pf.get_phases();

  const double rel_tol = 1.e-3;
  const double abs_tol = 1.;

  SECTION("Check number of phases")

  REQUIRE(phases.size() == 3);

  SECTION("Check the temperatures at which phases begin and end")

  CHECK(phases[0].T.front() == Approx(223.5400009155).epsilon(rel_tol));
  CHECK(phases[0].T.back() == Approx(1000.).epsilon(rel_tol));

  CHECK(phases[1].T.front() == Approx(77.6090621948).epsilon(rel_tol));
  CHECK(phases[1].T.back() == Approx(223.0857467651).epsilon(rel_tol));

  CHECK(phases[2].T.front() == Approx(0.).margin(abs_tol));
  CHECK(phases[2].T.back() == Approx(117.2025680542).epsilon(rel_tol));
}

TEST_CASE("Check phases for Z2 scalar singlet model", "[Z2ScalarSingletModel]") {

  LOGGER(fatal);

  // Construct our model
  EffectivePotential::Z2ScalarSingletModel model;
  model.set_m_s(90);
  model.set_lambda_hs(0.4);

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
  pf.set_seed(1);
  pf.find_phases();
  auto phases = pf.get_phases();

  const double rel_tol = 1.e-3;
  const double abs_tol = 1.;

  SECTION("Check number of phases")

  REQUIRE(phases.size() == 3);

  SECTION("Check the temperatures at which phases begin and end")

  CHECK(phases[0].T.front() == Approx(210.0191879272).epsilon(rel_tol));
  CHECK(phases[0].T.back() == Approx(1000.).epsilon(rel_tol));

  CHECK(phases[1].T.front() == Approx(0.).margin(abs_tol));
  CHECK(phases[1].T.back() == Approx(208.4613037109).epsilon(rel_tol));

  CHECK(phases[2].T.front() == Approx(0.).margin(abs_tol));
  CHECK(phases[2].T.back() == Approx(120.9191894531).epsilon(rel_tol));
}
