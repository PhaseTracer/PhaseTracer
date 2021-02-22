/**
  Example program for PhaseTracer.

*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "models/xSM_HT.hpp" 
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"

const std::string output_file_name = "xSM_HT_Results.txt";

int main(int argc, char* argv[]) {

  const bool debug_mode = argc == 1;
  std::ofstream output_file;
  
  double bins_lambda_hs;
  double lambda_hs;
  if (debug_mode){
    LOGGER(debug);
    bins_lambda_hs = 1;
    lambda_hs = 0.31;
  } else {
    LOGGER(fatal);
    bins_lambda_hs = 100;
    output_file.open(output_file_name);
  }

  double lambda_s = 0.0;
  double ms = SM::mh/2.;
  
  
  lambda_hs = 0.4;
  ms = 60;
  lambda_s =  0.16;
  
  for (double ii = 0; ii < bins_lambda_hs; ii++) {
    if (not debug_mode){
      lambda_hs = 0.4 / bins_lambda_hs * ii+0.2;
    }
    
    // Construct our model
    EffectivePotential::xSM_HT model(lambda_hs, lambda_s, ms);
    if (debug_mode)
      std::cout << model.check() << std::endl;
      
    // Make PhaseFinder object and find the phases
    PhaseTracer::PhaseFinder pf(model);
      
    pf.set_seed(0);
    pf.set_check_hessian_singular(false);
    pf.find_phases();
    if (debug_mode)
      std::cout << pf;

    // Make TransitionFinder object and find the transitions
    PhaseTracer::TransitionFinder tf(pf);
    tf.find_transitions();
    if (debug_mode)
      std::cout << tf;
      
    // Write results to data file
    if (debug_mode) {
      std::cout << lambda_hs << "\t"
              << model.get_TC_from_expression() << "\t"
              << model.get_vh_from_expression() << "\t"
              << model.get_vs_from_expression() << "\t";
    } else {
//      output_file << lambda_hs << "\t"
//              << model.get_TC_from_expression() << "\t"
//              << model.get_vh_from_expression() << "\t"
//              << model.get_vs_from_expression() << "\t";
      
      auto t = tf.get_transitions();
      if (t.size()==0){
        std::cout << "lambda_hs = " << lambda_hs << " found 0 transition!" << std::endl;
        output_file << lambda_hs << "\t" << 0 << "\t";
        output_file << -1 << "\t" << 0 << "\t"
                    << 0 << "\t" << 0 << "\t"
                    << 0 << "\t";
        output_file << lambda_s << "\t" << ms;
        output_file << std::endl;
        continue;
      }
      int jj=0;
      double gamme_max=0;
      for (int i=0; i<t.size(); i++) {
        double gamma = t[i].gamma;
        if (gamme_max < gamma){
          jj = i;
          gamme_max = gamma;
        }
      }
      output_file << lambda_hs << "\t" << t.size() << "\t";
      output_file << t[jj].TC << "\t" << t[jj].true_vacuum[0] << "\t"
                  << t[jj].true_vacuum[1] << "\t" << t[jj].false_vacuum[0] << "\t"
                  << t[jj].false_vacuum[1] << "\t";
      output_file << lambda_s << "\t" << ms;
      output_file << std::endl;
      
    }
      
    auto trans = tf.get_transitions();
    if (debug_mode) std::cout << trans.size() << "\t";
    for (auto &t : tf.get_transitions()) {
      if (debug_mode) {
        std::cout << t.TC << "\t" << t.true_vacuum[0] << "\t"
                  << t.true_vacuum[1] << "\t" << t.false_vacuum[0] << "\t"
                  << t.false_vacuum[1] << "\t";
//      }else{
//        output_file << t.TC << "\t" << t.true_vacuum[0] << "\t"
//                  << t.true_vacuum[1] << "\t" << t.false_vacuum[0] << "\t"
//                  << t.false_vacuum[1] << "\t";
      }
    }
    if (debug_mode) {
      std::cout << std::endl;
      // Print the data in a particular format for plotting
      PhaseTracer::phase_plotter(tf, "xSM_HT");
    } else {
      output_file << std::endl;
    }
    
  }
  if (not debug_mode) output_file.close();
  return 0;
}
