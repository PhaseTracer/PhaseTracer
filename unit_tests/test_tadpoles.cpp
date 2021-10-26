#include "catch/catch.hpp"
#include "models/xSM_MSbar.hpp"
#include "models/SM_parameters.hpp"
#include "logger.hpp"

std::vector<double> tadpole_result(bool use_tree_level_tadpole, bool use_covariant_gauge) {
  // Construct our model
  const double lambda_hs = 0.3;
  const double Q = SM::mtop;
  const double xi = 0.5;
  Eigen::ArrayXd physical(2);
  physical << SM::v, 0.;

  double ms = 0.5 * SM::mh;
  double lambda_s_min = 2. / square(SM::mh * SM::v) *
      square(square(ms) - 0.5 * lambda_hs * square(SM::v));
  double lambda_s = lambda_s_min + 0.1;

  auto model = EffectivePotential::xSM_MSbar::from_tadpoles(lambda_hs, lambda_s, ms, Q, xi,
              use_covariant_gauge, false, false, use_tree_level_tadpole);

  // Find singlet masses
  const auto mass_sq = use_tree_level_tadpole ? model.get_tree_scalar_masses_sq(physical) : model.get_1l_scalar_masses_sq(physical, 0.);
  const double mh_out = std::sqrt(mass_sq[1]);
  const double ms_out = std::sqrt(mass_sq[0]);

  return {mh_out, ms_out};
}

TEST_CASE("Z2 scalar singlet massses from tadpoles", "[TadpoleSolver]") {
  LOGGER(fatal);

  SECTION("covariant gauge") {
    for (bool use_tree_level_tadpole : {false, true}) {
      const auto cov = tadpole_result(use_tree_level_tadpole, true);
      CHECK(cov[0] == Approx(SM::mh));
      CHECK(cov[1] == Approx(0.5 * SM::mh));
    }
  }

  SECTION("R_xi gauge") {
    for (bool use_tree_level_tadpole : {false, true}) {
      const auto rxi = tadpole_result(use_tree_level_tadpole, false);
      CHECK(rxi[0] == Approx(SM::mh));
      CHECK(rxi[1] == Approx(0.5 * SM::mh));
    }
  }
}
