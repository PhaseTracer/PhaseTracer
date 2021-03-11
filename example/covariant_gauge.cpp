/**
  Example of covariant gauge using xSM model.
*/

#include <iostream>
#include <iomanip>

#include "models/xSM_MSbar_covariant.hpp"
#include "models/SM_parameters.hpp"
#include "transition_finder.hpp"
#include "phase_finder.hpp"
#include "logger.hpp"


int main(int argc, char* argv[]) {
  LOGGER(debug);

  const double lambda_hs = atof(argv[1]);
  const double Q = atof(argv[2]);
  const double xi = atof(argv[3]);
  const bool tree_level_tadpoles = atoi(argv[4]);
  const bool use_1L_EWSB_in_0L_mass = atoi(argv[5]);
  const bool physical_vacuum = true;

  std::cout << "lambda_hs = " << lambda_hs << std::endl
            << "Q = " << Q << std::endl
            << "xi = " << xi << std::endl
            << "tree-level tadpoles = " << tree_level_tadpoles << std::endl
            << "use 1-level ewsb in tree-level masses = " << use_1L_EWSB_in_0L_mass << std::endl;

  double ms = 0.5 * SM::mh;
  double lambda_s_min = 2. / square(SM::mh * SM::v) *
      square(square(ms) - 0.5 * lambda_hs * square(SM::v));
  double lambda_s = lambda_s_min + 0.1;

  // Construct our model
  auto model = EffectivePotential::xSM_MSbar_covariant::from_tadpoles(lambda_hs, lambda_s, ms, Q, xi, tree_level_tadpoles);

  model.set_daisy_method(EffectivePotential::DaisyMethod::None);
  model.set_use_1L_EWSB_in_0L_mass(use_1L_EWSB_in_0L_mass);

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder hb(model);
  hb.find_phases();
  std::cout << hb;

  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(hb);
  tf.set_TC_tol_rel(1e-12);
  tf.find_transitions();
  auto t = tf.get_transitions();
  std::sort(t.begin(), t.end(), [](const PhaseTracer::Transition& a, const PhaseTracer::Transition& b) {
    return a.gamma < b.gamma;
  });
  std::cout << std::setprecision(15) << t.back();

  return 0;
}
