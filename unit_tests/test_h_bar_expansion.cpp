#include "catch/catch.hpp"
#include "models/xSM_MSbar.hpp"
#include "models/SM_parameters.hpp"
#include "h_bar_expansion.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"


TEST_CASE("h-bar expansion method", "[HBarExpansion]") {

  LOGGER(fatal);

  // Construct our model
  const double lambda_hs = 0.3;
  const double Q = 173.03;
  const double xi = 0.;
  const bool tree_level_tadpoles = true;
  const bool tree_ewsb = false;

  // Construct our model
  EffectivePotential::xSM_MSbar model(lambda_hs, Q, tree_level_tadpoles);
  model.set_daisy_method(EffectivePotential::DaisyMethod::None);
  model.set_xi(xi);
  model.set_tree_ewsb(tree_ewsb);

  // Make HBarExpansion object and find the phases
  PhaseTracer::HbarExpansion hb(model);
  hb.set_seed(0);
  Eigen::ArrayXd pseudo(2);
  pseudo << 0., model.get_v_tree_s();
  hb.add_pseudo_phase(pseudo);
  std::cerr << "find phases" << std::endl;
  hb.find_phases();

  
}
