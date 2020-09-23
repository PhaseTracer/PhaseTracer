/**
  Example of h-bar expansion using xSM model.
*/

#include <iostream>

#include "models/xSM_MSbar.hpp"
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
  std::cout << "lambda_hs = " << lambda_hs << std::endl
            << "Q = " <<  Q << std::endl;

  // Construct our model
  auto model = EffectivePotential::xSM_MSbar(lambda_hs, Q, true);
  model.set_daisy_method(EffectivePotential::DaisyMethod::None);

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

  return 0;
}
