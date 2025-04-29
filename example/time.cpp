/**
  Time calls to PhaseTracer
*/

#include "models/2D_test_model.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"
#include "time.hpp"

int main(int argc, char *argv[]) {

  LOGGER(fatal)

  // Construct model
  EffectivePotential::TwoDimModel model;

  START_TIMER

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
  pf.find_phases();

  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(pf);
  tf.find_transitions();

  STOP_TIMER

  return 0;
}
