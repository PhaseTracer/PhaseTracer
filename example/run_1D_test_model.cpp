/**
  1D example program for PhaseTracer.
*/

#include <iostream>

#include "models/1D_test_model.hpp" // Located in effective-potential/include/models
#include "phasetracer.hpp"

int main(int argc, char *argv[]) {

  const bool debug_mode = argc > 1 && strcmp(argv[1], "-d") == 0;

  // Set level of screen  output
  if (debug_mode) {
    LOGGER(debug);
  } else {
    LOGGER(fatal);
  }

  // Construct model
  EffectivePotential::OneDimModel model;

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
  pf.find_phases();
  std::cout << pf;

  // Make ActionCalculator object
  PhaseTracer::ActionCalculator ac(pf);

  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(pf, ac);
  tf.find_transitions();
  std::cout << tf;

  // Make GravWaveCalculator object
  PhaseTracer::GravWaveCalculator gc(tf);
  gc.set_T_threshold_bubble_collision(1e10);

  const auto sps = gc.calc_spectrums();
  for (size_t ii = 0; ii < sps.size(); ii++) {
    std::cout << sps[ii];
    if (debug_mode) {
      gc.write_spectrum_to_text(sps[ii], "GW_spectrum_1D_test_model" + std::to_string(ii) + ".txt");
    }
  }

  if (debug_mode) {
    PhaseTracer::spectrum_plotter(gc, "1D_test");
  }

  if (debug_mode) {
    // Get actions between the two phases
    std::ofstream outputFile("action_1D_test_model.txt");
    LOGGER(fatal);
    const auto trans = tf.get_transitions();
    auto phase1 = trans[0].true_phase;
    auto phase2 = trans[0].false_phase;
    double min_T = std::max(trans[0].true_phase.T[0], trans[0].false_phase.T[0]);
    for (double Ttest = min_T; Ttest < trans[0].TC; Ttest += 0.1) {
      auto s = tf.get_action(phase1, phase2, Ttest);
      auto vacua = tf.get_vacua_at_T(phase1, phase2, Ttest);
      outputFile << vacua[0][0] << " " << vacua[1][0] << " " << Ttest << " " << s << " " << s / Ttest << std::endl;
    }
    outputFile.close();

    double action = ac.get_action(trans[0].true_vacuum_TN, trans[0].false_vacuum_TN, trans[0].TN);
    std::cout << "action = " << std::setprecision(15) << action << std::endl;

    std::ofstream file1("profile1D_for_1d_example.txt");
    file1 << "R,Phi,dPhi\n";
    PhaseTracer::Profile1D profile = ac.get_bubble_profile();
    for (int i = 0; i < profile.R.size(); ++i) {
      file1 << std::setprecision(10) << profile.R[i] << "," << profile.Phi[i] << "," << profile.dPhi[i] << std::endl;
    }
    file1.close();
  }

  return 0;
}
