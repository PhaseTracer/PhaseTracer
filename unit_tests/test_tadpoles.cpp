#include "catch/catch.hpp"
#include "models/xSM_MSbar.hpp"
#include "models/SM_parameters.hpp"
#include "logger.hpp"

std::vector<double> tadpole_result(bool use_covariant_gauge,
                                   bool use_1L_EWSB_in_0L_mass,
                                   bool use_Goldstone_resum,
                                   bool use_tree_level_tadpole) {
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
                use_covariant_gauge, use_1L_EWSB_in_0L_mass, use_Goldstone_resum, use_tree_level_tadpole);

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
      for (bool use_Goldstone_resum : {false, true}) {
        for (bool use_1L_EWSB_in_0L_mass : {false, true}) {
          if (use_Goldstone_resum==use_1L_EWSB_in_0L_mass) continue;
          INFO("1L in 0L mass = " + std::to_string(use_1L_EWSB_in_0L_mass)
                + ". goldstone = " + std::to_string(use_Goldstone_resum)
                + ". tree tadpole = " + std::to_string(use_tree_level_tadpole));
          const auto rxi = tadpole_result(true, use_1L_EWSB_in_0L_mass, use_Goldstone_resum, use_tree_level_tadpole);
          CHECK(rxi[0] == Approx(SM::mh));
          CHECK(rxi[1] == Approx(0.5 * SM::mh));
        }
      }
    }
  }

  SECTION("R_xi gauge") {
    for (bool use_tree_level_tadpole : {false, true}) {
      for (bool use_Goldstone_resum : {false, true}) {
        for (bool use_1L_EWSB_in_0L_mass : {false, true}) {
          INFO("1L in 0L mass = " + std::to_string(use_1L_EWSB_in_0L_mass)
                + ". goldstone = " + std::to_string(use_Goldstone_resum)
                + ". tree tadpole = " + std::to_string(use_tree_level_tadpole));
          const auto rxi = tadpole_result(false, use_1L_EWSB_in_0L_mass, use_Goldstone_resum, use_tree_level_tadpole);
          CHECK(rxi[0] == Approx(SM::mh));
          CHECK(rxi[1] == Approx(0.5 * SM::mh));
        }
      }
    }
  }
}
