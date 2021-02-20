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
  
  const bool debug_mode = argc == 1;
  
  double bins_lambda_hs;
  double lambda_hs;
  std::ofstream output_file;
  
  if (debug_mode){
    LOGGER(debug);
    bins_lambda_hs = 1;
    lambda_hs = 0.31;
  } else {
    LOGGER(fatal);
    bins_lambda_hs = 100;
    output_file.open("PRM.txt");
  }
  
  const double Q = SM::mtop;
  const double xi = 0;
  const bool tree_level_tadpoles = false;
  const bool tree_ewsb = false;
  const bool physical_vacuum = true;

  std::cout << "lambda_hs = " << lambda_hs << std::endl
            << "Q = " << Q << std::endl
            << "xi = " << xi << std::endl
            << "tree-level tadpoles = " << tree_level_tadpoles << std::endl
            << "tree-level ewsb parameters = " << tree_ewsb << std::endl;
  
  for (double ii = 0; ii < bins_lambda_hs; ii++) {
    if (not debug_mode){
      lambda_hs = 0.4 / bins_lambda_hs * ii+0.2;
    }
  
    double ms = 0.5 * SM::mh;
    double lambda_s_min = 2. / square(SM::mh * SM::v) *
        square(square(ms) - 0.5 * lambda_hs * square(SM::v));
    double lambda_s = lambda_s_min + 0.1;
    
    // Construct our model
    auto model = EffectivePotential::xSM_MSbar::from_tadpoles(lambda_hs, lambda_s, ms, Q, xi, tree_level_tadpoles, tree_ewsb);

    model.set_daisy_method(EffectivePotential::DaisyMethod::ArnoldEspinosa);
//    model.set_daisy_method(EffectivePotential::DaisyMethod::Parwani);

    // Make PhaseFinder object and find the phases
    PhaseTracer::HbarExpansion hb(model);
    hb.set_seed(0);
    Eigen::ArrayXd pseudo(2);
    pseudo << 0., model.get_v_tree_s();
    hb.add_pseudo_phase(pseudo);
    hb.find_phases();

    // Make TransitionFinder object and find the transitions
    PhaseTracer::TransitionFinder tf(hb);
    tf.set_TC_tol_rel(1e-12);
    tf.find_transitions();

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
    
    
    if (debug_mode) {
      std::cout << hb;
      std::cout << std::setprecision(15) << tf;
      std::cout << "gamma_HT = " << gamma << std::endl;
    } else {
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
      output_file << lambda_s << "\t" << ms << "\t" << gamma;
      output_file << std::endl;
    }

//    // Get Higgs mass
//    const auto minima_tree = hb.find_minima_at_t(0.);
//    PhaseTracer::PhaseFinder pf(model);
//    const auto minima_1l = pf.find_minima_at_t(0.);
//    std::cout << "v_tree = " << std::abs(minima_tree[0].x[0]) << std::endl;
//    std::cout << "v_1l = " << std::abs(minima_1l[0].x[0]) << std::endl;
//
//    Eigen::ArrayXd physical(2);
//    physical << SM::v, 0.;
//    const auto mass_sq_1l = physical_vacuum ? model.get_1l_scalar_masses_sq(physical, 0.) :
//                                            model.get_1l_scalar_masses_sq(minima_1l[0].x, 0.);
//    const auto mass_sq_tree = physical_vacuum ? model.get_tree_scalar_masses_sq(physical):
//                                              model.get_tree_scalar_masses_sq(minima_tree[0].x);
//    std::cout << "mh_tree = " << std::sqrt(mass_sq_tree[1]) << std::endl;
//    std::cout << "mh_1l = " << std::sqrt(mass_sq_1l[1]) << std::endl;
  
  }
  
  output_file.close();
  return 0;
}
