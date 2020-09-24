/**
  Example of h-bar expansion using xSM model.
*/

#include <iostream>

#include "models/xSM_MSbar.hpp"
#include "models/SM_parameters.hpp"
#include "transition_finder.hpp"
#include "h_bar_expansion.hpp"
#include "phase_plotter.hpp"
#include "potential_plotter.hpp"
#include "potential_line_plotter.hpp"
#include "logger.hpp"


int main(int argc, char* argv[]) {
  LOGGER(debug);

  const double lambda_hs = atof(argv[1]);
  const double Q = atof(argv[2]);
  const double xi = atof(argv[3]);
  const bool tree_level_tadpoles = atoi(argv[4]);
  std::cout << "lambda_hs = " << lambda_hs << std::endl
            << "Q = " << Q << std::endl
            << "xi = " << xi << std::endl
            << "tree-level tadpoles = " << tree_level_tadpoles << std::endl;

  // Construct our model
  auto model = EffectivePotential::xSM_MSbar(lambda_hs, Q, tree_level_tadpoles);
  model.set_daisy_method(EffectivePotential::DaisyMethod::None);
  model.set_xi(xi);

  // Make PhaseFinder object and find the phases
  PhaseTracer::HbarExpansion hb(model);
  hb.set_seed(0);
  Eigen::ArrayXd pseudo(2);
  pseudo << 0., model.get_v_tree_s();
  hb.add_pseudo_phase(pseudo);
  hb.find_phases();
  std::cout << hb;

  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(hb);
  tf.set_TC_tol_rel(1e-8);
  tf.find_transitions();
  std::cout << std::setprecision(15) << tf;

  // Find gamma using high-temperature expansion
  PhaseTracer::HTExpansion ht(model);
  ht.set_seed(0);
  const double TC = tf.get_transitions()[0].TC;
  const auto ht_minima = ht.find_minima_at_t(TC);

  // Use minima with greatest Higgs
  double delta = 0.;
  for (const auto& m : ht_minima) {
    delta = std::max(delta, std::abs(m.x(0)));
  }

  const double gamma = delta / TC;
  std::cout << "gamma_HT = " << gamma << std::endl;

  // Get Higgs mass
  const auto minima_tree = hb.find_minima_at_t(0.);
  PhaseTracer::PhaseFinder pf(model);
  const auto minima_1l = pf.find_minima_at_t(0.);
  std::cout << "v_tree = " << std::abs(minima_tree[0].x[0]) << std::endl;
  std::cout << "v_1l = " << std::abs(minima_1l[0].x[0]) << std::endl;  

  Eigen::ArrayXd physical(2);
  physical << SM::v, 0.;
  const auto mass_sq_1l = model.get_1l_scalar_masses_sq(physical, 0.);
  const auto mass_sq_tree = model.get_tree_scalar_masses_sq(physical);
  std::cout << "mh_tree = " << std::sqrt(mass_sq_tree[1]) << std::endl;
  std::cout << "mh_1l = " << std::sqrt(mass_sq_1l[1]) << std::endl;

  return 0;
}
