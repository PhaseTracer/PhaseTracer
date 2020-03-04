/**
  Example program for PhaseTracer.

  We scan two parameters in a Z2 real scalar singlet extension of
  the Standard Model and compare the critical temperatures with
  the expressions in arXiv:1611.02073.
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "models/Z2_scalar_singlet_model.hpp" // Located in effective-potential/include/models
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"

const std::string output_file_name = "Z2ScalarSingletModel_Results.txt";

int main() {

  LOGGER(fatal);

  // Construct our model
  EffectivePotential::Z2ScalarSingletModel model;

  std::ofstream output_file;
  output_file.open(output_file_name);
  output_file << "l_hs \tm_s \tT_c^EX \tv_h^EX \tv_s^EX "
                 "\tT1_c^PT \tv1_h^T \tv1_s^T \tv1_h^F \tv1_s^F"
//                 "\tT2_c^PT \tv2_h^T \tv2_s^T \tv2_h^F \tv2_s^F"
              << std::endl;

  const double bins_lambda_hs = 300;
  const double bins_ms = 300;

  for (double ii = 0; ii < bins_lambda_hs; ii++) {
    for (double jj = 0; jj < bins_ms; jj++) {
      const double lambda_hs = 2. / bins_lambda_hs * ii;
      const double m_s = 250. / bins_ms * jj;

      model.set_m_s(m_s);
      model.set_lambda_hs(lambda_hs);

      if (!model.check()) {
        continue;
      }

//      std::cout << "lambda_hs = " << lambda_hs << ". m_s = " << m_s << std::endl;

      // Make PhaseFinder object and find the phases
      PhaseTracer::PhaseFinder pf(model);
      pf.set_seed(0);
      pf.set_check_hessian_singular(false);
      pf.find_phases();

      // Make TransitionFinder object and find the transitions
      PhaseTracer::TransitionFinder tf(pf);
      tf.find_transitions();

      if (pf.get_phases().size()!=2 or tf.get_transitions().size() !=1) {
        std::cout << "n_p=" << pf.get_phases().size() << ", n_t=" << tf.get_transitions().size() << ", lambda_hs = " << lambda_hs << ", m_s = " << m_s << std::endl;
      }
      
      // Write results to data file
      output_file << lambda_hs << "\t" << m_s << "\t"
                  << model.get_TC_from_expression() << "\t"
                  << model.get_vh_from_expression() << "\t"
                  << model.get_vs_from_expression() << "\t";
      for (auto &t : tf.get_transitions()) {
        output_file << t.TC << "\t" << t.true_vacuum[0] << "\t"
                    << t.true_vacuum[1] << "\t" << t.false_vacuum[0] << "\t"
                    << t.false_vacuum[1] << "\t";
      }
      output_file << std::endl;
    }
  }
  output_file.close();
  return 0;
}
