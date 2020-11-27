#include "catch/catch.hpp"
#include "models/xSM_MSbar.hpp"
#include "models/SM_parameters.hpp"
#include "logger.hpp"


TEST_CASE("Z2 scalar singlet massses from tadpoles", "[TadpoleSolver]") {

  LOGGER(fatal);

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
  
  // Construct models with tree-level and one-loop tadpoles
  auto tree = EffectivePotential::xSM_MSbar::from_tadpoles(lambda_hs, lambda_s, ms, Q, xi, true);
  auto one_loop = EffectivePotential::xSM_MSbar::from_tadpoles(lambda_hs, lambda_s, ms, Q, xi, false);

  // Find singlet masses
  const auto mass_sq_1l = one_loop.get_1l_scalar_masses_sq(physical, 0.);
  const double mh_1l = std::sqrt(mass_sq_1l[1]);
  const double ms_1l = std::sqrt(mass_sq_1l[0]);
  const auto mass_sq_tree = tree.get_tree_scalar_masses_sq(physical);
  const double mh_tree = std::sqrt(mass_sq_tree[1]);
  const double ms_tree = std::sqrt(mass_sq_tree[0]); 

  // Check them to constraints
  CHECK(mh_1l == Approx(SM::mh));
  CHECK(mh_tree == Approx(SM::mh));
  CHECK(ms_1l == Approx(0.5 * SM::mh));
  CHECK(ms_tree == Approx(0.5 * SM::mh));
}
