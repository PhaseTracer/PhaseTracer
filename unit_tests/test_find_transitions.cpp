#include "catch/catch.hpp"
#include "models/1D_test_model.hpp"
#include "models/2D_test_model.hpp"
#include "models/Z2_scalar_singlet_model.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "action_calculator.hpp"
#include "logger.hpp"

TEST_CASE("Compute transitions for a one-dimensional model", "[1DTestModel]") {

  LOGGER(fatal);

  // Construct our model
  EffectivePotential::OneDimModel model;

  // Since we have analytic results, let's turn up our precision and check we can
  // get very close
  const double rel_tol = 1.e-8;
  const double abs_tol = 1.e-5;

  // Make PhaseFinder object and find the transitions
  PhaseTracer::PhaseFinder pf(model);
  pf.set_seed(1);
  pf.set_find_min_x_tol_rel(1.e-8);
  pf.set_find_min_x_tol_abs(1.e-8);
  pf.find_phases();
  PhaseTracer::TransitionFinder tf(pf);
  tf.set_TC_tol_rel(1e-16);
  tf.find_transitions();

  auto transitions = tf.get_transitions();

  SECTION("Check number of transitions")

  REQUIRE(transitions.size() == 1);

  SECTION("Check the critical temperatures")

  CHECK(transitions[0].TC == Approx(model.get_TC_from_expression()).epsilon(rel_tol));

  SECTION("Check the true vacuum")

  CHECK(transitions[0].true_vacuum[0] == Approx(model.get_true_vacuum_from_expression()).epsilon(rel_tol));

  SECTION("Check the false vacuum")

  CHECK(transitions[0].false_vacuum[0] == Approx(model.get_false_vacuum_from_expression()).margin(abs_tol));
}

TEST_CASE("Compute nucleation for a one-dimensional model", "[1DTestModel]") {

  LOGGER(fatal);

  // Construct our model
  EffectivePotential::OneDimModel model;

  const double rel_tol = 1.e-3;
  const double abs_tol = 1.e-3;

  // Make PhaseFinder object and find the transitions
  PhaseTracer::PhaseFinder pf(model);
  pf.set_seed(1);
  pf.set_find_min_x_tol_rel(1.e-8);
  pf.set_find_min_x_tol_abs(1.e-8);
  pf.find_phases();

  PhaseTracer::ActionCalculator ac(model);
  PhaseTracer::TransitionFinder tf(pf, ac);
  tf.set_TC_tol_rel(1e-16);
  tf.find_transitions();

  auto transitions = tf.get_transitions();

  SECTION("Check number of transitions")

  REQUIRE(transitions.size() == 1);

  SECTION("Check the nucleation temperature")

  CHECK(transitions[0].TN == Approx(57.4280983265).epsilon(rel_tol));

  SECTION("Check the true vacuum at nucleation")

  CHECK(std::abs(transitions[0].true_vacuum_TN[0]) == Approx(0.).margin(abs_tol));

  SECTION("Check the false vacuum at nucleation")

  CHECK(std::abs(transitions[0].false_vacuum_TN[0]) == Approx(53.5392255174).epsilon(rel_tol));
}

TEST_CASE("Compute nucleation for a two-dimensional model", "[2DTestModel]") {

  LOGGER(fatal);

  // Construct our model
  EffectivePotential::TwoDimModel model;

  // Make PhaseFinder object and find the transitions
  PhaseTracer::PhaseFinder pf(model);
  pf.set_seed(1);
  pf.find_phases();

  PhaseTracer::ActionCalculator ac(model);
  PhaseTracer::TransitionFinder tf(pf, ac);
  tf.find_transitions();
  auto transitions = tf.get_transitions();

  const double rel_tol = 1.e-3;

  SECTION("Check number of transitions")

  REQUIRE(transitions.size() == 1);

  SECTION("Check the nucleation temperature")

  CHECK(transitions[0].TN == Approx(84.1901920318).epsilon(rel_tol));

  SECTION("Check the true vacuum at nucleation")

  CHECK(std::abs(transitions[0].true_vacuum_TN[0]) == Approx(231.122003632).epsilon(rel_tol));
  CHECK(std::abs(transitions[0].true_vacuum_TN[1]) == Approx(136.7545906901).epsilon(rel_tol));

  SECTION("Check the false vacuum at nucleation")

  CHECK(std::abs(transitions[0].false_vacuum_TN[0]) == Approx(286.4206084948).epsilon(rel_tol));
  CHECK(std::abs(transitions[0].false_vacuum_TN[1]) == Approx(382.254323455).epsilon(rel_tol));
}

TEST_CASE("Compute transitions for a two-dimensional model", "[2DTestModel]") {

  LOGGER(fatal);

  // Construct our model
  EffectivePotential::TwoDimModel model;

  // Make PhaseFinder object and find the transitions
  PhaseTracer::PhaseFinder pf(model);
  pf.set_seed(1);
  pf.find_phases();
  PhaseTracer::TransitionFinder tf(pf);
  tf.find_transitions();
  auto transitions = tf.get_transitions();

  const double rel_tol = 1.e-3;
  const double abs_tol = 1.;

  SECTION("Check number of transitions")

  REQUIRE(transitions.size() == 2);

  SECTION("Check the critical temperatures")

  CHECK(transitions[0].TC == Approx(109.4084075679).epsilon(rel_tol));
  CHECK(transitions[1].TC == transitions[0].TC);

  SECTION("Check the true vacuum")

  CHECK(std::abs(transitions[0].true_vacuum[0]) == Approx(263.4882058492).epsilon(rel_tol));
  CHECK(std::abs(transitions[0].true_vacuum[1]) == Approx(314.6546897515).epsilon(rel_tol));
  CHECK(transitions[0].true_vacuum[0] * transitions[0].true_vacuum[1] > 0);

  SECTION("Check the false vacuum")

  CHECK(std::abs(transitions[0].false_vacuum[0]) == Approx(220.0218876547).epsilon(rel_tol));
  CHECK(std::abs(transitions[0].false_vacuum[1]) == Approx(150.0149642738).epsilon(rel_tol));
  CHECK(transitions[0].false_vacuum[0] * transitions[0].false_vacuum[1] < 0);
}

TEST_CASE("Check transitions for Z2 scalar singlet model", "[Z2ScalarSingletModel]") {

  LOGGER(fatal);

  // Construct our model
  EffectivePotential::Z2ScalarSingletModel model;

  // Since we have analytic results, let's turn up our precision and check we can
  // get very close
  const double rel_tol = 1.e-8;
  const double abs_tol = 1.e-5;

  const double bins_lambda_hs = 30;
  const double bins_ms = 30;

  SECTION("Loop over many model parameters and check each one")

  for (double ii = 0; ii < bins_lambda_hs; ii++) {
    for (double jj = 0; jj < bins_ms; jj++) {

      // Set model parameters
      const double lambda_hs = 2. / bins_lambda_hs * ii;
      const double m_s = 250. / bins_ms * jj;
      model.set_m_s(m_s);
      model.set_lambda_hs(lambda_hs);

      if (!model.check()) {
        continue;
      }

      // Make PhaseFinder object and find the phases
      PhaseTracer::PhaseFinder pf(model);
      pf.set_seed(1);
      pf.set_find_min_x_tol_rel(1.e-8);
      pf.set_find_min_x_tol_abs(1.e-8);
      pf.find_phases();
      PhaseTracer::TransitionFinder tf(pf);
      tf.set_TC_tol_rel(1e-16);
      tf.find_transitions();
      auto transitions = tf.get_transitions();

      // Perform checks
      REQUIRE(transitions.size() == 1);
      CHECK(transitions[0].TC == Approx(model.get_TC_from_expression()).epsilon(rel_tol));
      CHECK(std::abs(transitions[0].true_vacuum[0]) == Approx(std::abs(model.get_vh_from_expression())).epsilon(rel_tol));
      CHECK(transitions[0].true_vacuum[1] == Approx(0.).margin(abs_tol));
      CHECK(std::abs(transitions[0].false_vacuum[1]) == Approx(std::abs(model.get_vs_from_expression())).epsilon(rel_tol));
      CHECK(transitions[0].false_vacuum[0] == Approx(0.).margin(abs_tol));
    }
  }
}
