#include <stdexcept>
#include <type_traits>
#include <vector>

#include "catch/catch.hpp"
#include "models/1D_test_model.hpp"
#include "runner.hpp"
#include "logger.hpp"

static_assert(!std::is_constructible_v<PhaseTracer::Runner, EffectivePotential::OneDimModel&&, PhaseTracer::Config>,
              "Runner must not accept a temporary model");
static_assert(!std::is_copy_constructible_v<PhaseTracer::Runner>, "Runner must not be copyable");

namespace {
PhaseTracer::Config one_dim_config()
{
  PhaseTracer::Config config;
  config.phase_finder.seed = 1;
  config.phase_finder.find_min_x_tol_rel = 1.e-8;
  config.phase_finder.find_min_x_tol_abs = 1.e-8;
  config.transition_finder.TC_tol_rel = 1e-16;
  config.pipeline.stop_after = PhaseTracer::Stage::TransitionFinder;
  config.pipeline.to_print = false;
  return config;
}
} // namespace

TEST_CASE("Runner reproduces the manual pipeline up to TransitionFinder", "[Runner]") {

  LOGGER(fatal);

  EffectivePotential::OneDimModel model;
  PhaseTracer::Runner runner(model, one_dim_config());
  CHECK(runner.run_id() == 0);

  const auto status = runner.run();
  REQUIRE(status.ok());
  CHECK(runner.run_id() == 1);

  CHECK_FALSE(runner.get_phases().empty());
  REQUIRE(runner.get_transitions().size() == 1);
  CHECK(runner.get_transitions()[0].TC == Approx(model.get_TC_from_expression()).epsilon(1.e-8));

  CHECK(runner.has(PhaseTracer::Stage::PhaseFinder));
  CHECK(runner.has(PhaseTracer::Stage::ActionCalculator));
  CHECK(runner.has(PhaseTracer::Stage::TransitionFinder));
  CHECK_FALSE(runner.has(PhaseTracer::Stage::ThermoFinder));
  CHECK_FALSE(runner.has(PhaseTracer::Stage::GravWave));
  CHECK_THROWS_AS(runner.get_thermal_parameters(), std::logic_error);
  CHECK_THROWS_AS(runner.get_spectra(), std::logic_error);

  SECTION("Wrappers return references into the stage objects") {
    CHECK(&runner.get_phases() == &runner.phase_finder().get_phases());
    CHECK(&runner.get_transitions() == &runner.transition_finder().get_transitions());
  }

  SECTION("run() rebuilds the stages") {
    const double TC = runner.get_transitions()[0].TC;
    runner.config().transition_finder.TC_tol_rel = 1e-12;
    REQUIRE(runner.run().ok());
    CHECK(runner.run_id() == 2);
    REQUIRE(runner.get_transitions().size() == 1);
    CHECK(runner.get_transitions()[0].TC == Approx(TC).epsilon(1.e-8));
  }
}

TEST_CASE("Runner reports failures", "[Runner]") {

  LOGGER(fatal);

  EffectivePotential::OneDimModel model;
  auto config = one_dim_config();

  SECTION("Stage accessors throw before run()") {
    PhaseTracer::Runner runner(model, config);
    CHECK_THROWS_AS(runner.phase_finder(), std::logic_error);
    CHECK_THROWS_AS(runner.get_transitions(), std::logic_error);
  }

  SECTION("Invalid config is returned as a status") {
    config.gravwave.min_frequency = 1.;
    config.gravwave.max_frequency = 1e-3;
    PhaseTracer::Runner runner(model, config);
    const auto status = runner.run();
    CHECK(status.code == PhaseTracer::StatusCode::InvalidConfig);
    CHECK(status.stage == PhaseTracer::Stage::Config);
    CHECK_FALSE(runner.has(PhaseTracer::Stage::PhaseFinder));
  }

  SECTION("Invalid config is thrown with throw_on_error") {
    config.gravwave.min_frequency = 1.;
    config.gravwave.max_frequency = 1e-3;
    config.pipeline.throw_on_error = true;
    PhaseTracer::Runner runner(model, config);
    CHECK_THROWS_AS(runner.run(), PhaseTracer::RunnerError);
    CHECK(runner.status().code == PhaseTracer::StatusCode::InvalidConfig);
  }

  SECTION("A filter rejecting every transition gives NoThermalParameters") {
    config.pipeline.stop_after = PhaseTracer::Stage::GravWave;
    config.thermo_finder.transition_filter = [](const std::vector<PhaseTracer::Transition>&) {
      return std::vector<PhaseTracer::Transition>{};
    };
    PhaseTracer::Runner runner(model, config);
    const auto status = runner.run();
    CHECK(status.code == PhaseTracer::StatusCode::NoThermalParameters);
    CHECK(status.stage == PhaseTracer::Stage::ThermoFinder);
    CHECK(runner.has(PhaseTracer::Stage::ThermoFinder));
    CHECK(runner.get_thermal_parameters().empty());
    CHECK_FALSE(runner.has(PhaseTracer::Stage::GravWave));
  }
}
