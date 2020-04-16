/**
  Example program for PhaseTracer.

  We scan a parameter in a two field model and show the phases and phase transitions
  for each value of the parameter.
*/

#include <iostream>
#include <string>
#include <vector>

#include "THDMIISNMSSMBCsimple.hpp"
#include "transition_finder.hpp"
#include "phase_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"

int main() {

  LOGGER(debug);

  // Construct our model
  EffectivePotential::THDMIISNMSSMBC model;

  // Set the input parameters for NMSSM
  std::vector<double> x(9);
  x[0] = 0.60057271091346375869; // lambda
  x[1] = 0.17539425572262870578; // kappa
  x[2] = 174.03961453074532528;  // Alambda
  x[3] = -25.183599724051667579; // Akappa
  x[4] = -24.598460281920118575; // Atop
  x[5] = 5857.4938270437187384;  // mstopL
  x[6] = 1619.4273820538846849;  // mstopR
  x[7] = 245.73342610981194412;  // vSIN
  x[8] = 2.5513456154594420511;  // TanBeta
  model.set_input(x);

  // Print masses of CP-even Higgs
  std::cout << "mh1 = " << model.get_mh()[0] << std::endl;

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
  pf.set_seed(3);
  pf.find_phases();
  std::cout << pf;

  // Make TransitionFinder object and find the transitions
  // This takes the phase finder object as an argument. We must have already
  // populated the phases by e.g. the find_phases method, as above
  PhaseTracer::TransitionFinder tf(pf);
  tf.find_transitions();
  std::cout << tf;

  // Print the data in a particular format for plotting
  PhaseTracer::phase_plotter(tf, "THDMIISNMSSMBCsimple");
}
