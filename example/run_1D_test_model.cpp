/**
  1D example program for PhaseTracer.
*/

#include <iostream>

#include "models/1D_test_model.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"
#include "gravwave_calculator.hpp"

int main(int argc, char* argv[]) {

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
  PhaseTracer::ActionCalculator ac(model);
  
  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(pf,ac);
  tf.find_transitions();
  std::cout << tf;
  
  

//

//    const auto phases = pf.get_phases();
//    auto phase1 = phases[0];
//    auto phase2 = phases[1];
//
//  std::cout << tf.get_action(phase2, phase1, 56.9986)  << std::endl;
//
//    Eigen::VectorXd true_vacuum(1);
//    Eigen::VectorXd false_vacuum(1);
//  true_vacuum[0] = 9.18412e-06;
//  false_vacuum[0] = 54.2887;
//  std::cout << tf.get_action(true_vacuum, false_vacuum, 56.9986) << std::endl;
//
//  false_vacuum[0] = 54.28871;
//  std::cout << tf.get_action(true_vacuum, false_vacuum, 56.9986) << std::endl;
//
  
  std::ofstream outputFile("action_1D_test_model.txt", std::ios_base::app);
  const auto phases = pf.get_phases();
  auto phase1 = phases[1];
  auto phase2 = phases[0];
  for (double Ttest=31.623288422919078; Ttest < 57.975077465648184; Ttest += 0.5){
      auto s = tf.get_action(phase1, phase2, Ttest);
      auto vacua = tf.get_vacua_at_T(phase1, phase2, Ttest);
      outputFile << vacua[0][0] << " " << vacua[1][0] << " " << Ttest << " " << s << " " << s/Ttest  << std::endl;
  }
  outputFile.close();
  
  PhaseTracer::GravWaveCalculator gc(tf);
  auto gw = gc.cal_spectrums();
  std::cout << gc;
  for (int ii=0; ii < gw.size(); ii++){
    gc.write_spectrum_to_text(gw[ii],"GW_results"+std::to_string(ii)+".txt");
  }
  
  gc.write_spectrum_to_text(gc.get_total_spectrum(),"GW_results.txt");
  
//  gc.Set_parameters(0.3,500.,0.9,40.);
//  auto results = gc.GW_total_spectrum();
//  gc.Write_to_csv(results, "GW_results.csv");
  
//  Eigen::VectorXd d(1);
//  d(0) = 24.5;
//  std::cout << model.dV_dx(d,56.8716) << std::endl;
  
//  if (debug_mode) {
//    PhaseTracer::phase_plotter(tf, "1D_test_model");
//  }

  return 0;
}
