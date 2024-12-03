#include "catch/catch.hpp"
#include "models/1D_test_model.hpp"
#include "models/2D_test_model.hpp"
#include "models/Z2_scalar_singlet_model.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "action_calculator.hpp"
#include "gravwave_calculator.hpp"
#include "logger.hpp"

TEST_CASE("Compute gravitational wave spectrum for a one-dimensional model", "[1DTestModel]") {

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

  PhaseTracer::GravWaveCalculator gc(tf);
  gc.calc_spectrums();
  auto spec = gc.get_spectrums();

  SECTION("Check number of spectra")

  REQUIRE(spec.size() == 1);

  SECTION("Check value of spectra")

  REQUIRE(spec[0].Tref == Approx(57.3603638993).epsilon(rel_tol));
  REQUIRE(spec[0].alpha == Approx(0.0013946917).epsilon(rel_tol));
  REQUIRE(spec[0].beta_H == Approx(6873.221784433).epsilon(rel_tol));
  REQUIRE(spec[0].peak_frequency == Approx(0.1191659843).epsilon(rel_tol));

  Catch::StringMaker<double>::precision = 30;
  REQUIRE(spec[0].peak_amplitude == Approx(3.8981023e-23).epsilon(rel_tol));
}

TEST_CASE("Compute gravitational wave spectrum for a two-dimensional model", "[2DTestModel]") {

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

  PhaseTracer::GravWaveCalculator gc(tf);
  gc.calc_spectrums();
  auto spec = gc.get_spectrums();

  const double rel_tol = 1.e-3;

  SECTION("Check number of spectra")

  REQUIRE(spec.size() == 1);

  SECTION("Check value of spectra")

  REQUIRE(spec[0].Tref == Approx(84.1901920318).epsilon(rel_tol));
  REQUIRE(spec[0].alpha == Approx(0.0970871041).epsilon(rel_tol));
  REQUIRE(spec[0].beta_H == Approx(1781.3022934312).epsilon(rel_tol));
  REQUIRE(spec[0].peak_frequency == Approx(0.0452179242).epsilon(rel_tol));

  Catch::StringMaker<double>::precision = 30;
  REQUIRE(spec[0].peak_amplitude == Approx(1.21499917010577e-16).epsilon(rel_tol));

}