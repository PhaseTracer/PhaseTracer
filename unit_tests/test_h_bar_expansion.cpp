#include "catch/catch.hpp"
#include "models/xSM_MSbar.hpp"
#include "models/SM_parameters.hpp"
#include "h_bar_expansion.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"

std::vector<double> h_bar_expansion(bool use_covariant_gauge, double xi) {
  // Construct our model
  const double lambda_hs = 0.3;
  const double Q = 173.03;

  double ms = 0.5 * SM::mh;
  double lambda_s_min = 2. / square(SM::mh * SM::v) *
                        square(square(ms) - 0.5 * lambda_hs * square(SM::v));
  double lambda_s = lambda_s_min + 0.1;

  // Construct our model
  auto model = EffectivePotential::xSM_MSbar::from_tadpoles(lambda_hs, lambda_s, ms, Q, xi,
                                                            use_covariant_gauge, false, false, true);
  model.set_daisy_method(EffectivePotential::DaisyMethod::None);

  // Make HBarExpansion object and find the phases
  PhaseTracer::HbarExpansion hb(model);
  hb.set_seed(0);
  Eigen::ArrayXd pseudo(2);
  pseudo << 0., model.get_v_tree_s();
  hb.add_pseudo_phase(pseudo);
  hb.find_phases();

  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(hb);
  tf.set_TC_tol_rel(1e-12);
  tf.find_transitions();

  // Find gamma using high-temperature expansion
  PhaseTracer::HTExpansion ht(model);
  ht.set_seed(0);
  const double TC = tf.get_transitions()[0].TC;

  const auto ht_minima = ht.find_minima_at_t(TC);

  // Use minima with greatest Higgs
  double delta = 0.;
  for (const auto &m : ht_minima) {
    delta = std::max(delta, std::abs(m.x(0)));
  }

  const double gamma = delta / TC;

  return {gamma, TC};
}

TEST_CASE("h-bar expansion method", "[HBarExpansion]") {
  LOGGER(fatal);

  SECTION("covariant gauge") {
    for (double xi : {0., 0.5, 1., 2.}) {
      const auto cov = h_bar_expansion(true, xi);
      CHECK(cov[0] == Approx(1.9653276239));
      CHECK(cov[1] == Approx(94.6132178548));
    }
  }

  SECTION("R_xi gauge") {
    for (double xi : {0., 0.5, 1., 2.}) {
      const auto rxi = h_bar_expansion(false, xi);
      CHECK(rxi[0] == Approx(1.9653276251));
      CHECK(rxi[1] == Approx(94.6132177564));
    }
  }
}
