/**
  Example program for PhaseTracer.

  We scan a parameter in a two field model and show the phases and phase transitions
  for each value of the parameter.
*/

#include <iostream>
#include <string>
#include <vector>

#include "ScalarSingletZ2DMMhInput_withSingletVEVinPT.hpp"
#include "transition_finder.hpp"
#include "phase_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"
#include "potential_plotter.hpp"

#include "SM_parameters.hpp"

int main() {

  LOGGER(debug);

  // Construct our model
  EffectivePotential::ScalarSingletZ2DMMhInput_withSingletVEVinPT model;

  model.set_use_1L_EWSB_in_0L_mass(false);
  model.set_use_Goldstone_resum(false);
  
  //PA: use example SLHA values for now
  const double mtop =  SM::mtop;
  double MhInput = SM::mh;
  double lambda_hs = 0.25;
  double lambda_s =  0.1;
  double mus2 = -3000;
  
  std::vector<double> x(9);
  x[0] = mtop; // Qin = 173.0
  x[1] = mtop; // QEWSB = 173.0
  x[2] = MhInput;  // MhInput
  x[3] = mus2; // muSInput
  x[4] = lambda_s; // LamSInput
  x[5] = lambda_hs;  // LamSHInput
  
  model.set_input(x);

  // Print masses of CP-even Higgs
  std::cout << std::setprecision(16);
  std::cout << "mhh = " << model.get_mh()[0] << std::endl;
  std::cout << "mss = " << model.get_mh()[1] << std::endl;
  std::cout << std::endl;

        Eigen::VectorXd test(2);
        test <<  SM::v, 0;
//        std::cout << "Numerically derivatives of the full potential at EW VEV:" << std::endl;
//        auto d2Vdh2 = model.d2V_dx2(x,0);
//        std::cout << std::setprecision(16);
//        std::cout << "Sqrt[d^2V/dh^2] = "<< std::sqrt(abs(d2Vdh2(0,0))) << std::endl;
//        std::cout << "Sqrt[d^2V/ds^2] = "<< std::sqrt(abs(d2Vdh2(1,1))) << std::endl;
      
//      PhaseTracer::potential_plotter(model, 0, "potential", -10., 300, 1, -5., 200., 1);
//      PhaseTracer::potential_plotter(model, 142.35, "potential", 0., 160, 0.2, -2., 160., 0.2);
    
    
    double Vtree = model.V0(test);
    double VCW = model.V1(test);
    std::cout << "Vtree = "<< Vtree << std::endl;
    std::cout << "VCW = "<< VCW << std::endl;
    
    return 0;
    // Make PhaseFinder object and find the phases
    PhaseTracer::PhaseFinder pf(model);
    auto minima_at_t_low = pf.get_minima_at_t_low();
    
    std::cout << "minima_at_t_low:" << std::endl;
    for (auto minimum : minima_at_t_low) std::cout << minimum << std::endl;


  
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
  PhaseTracer::phase_plotter(tf, "ScalarSingletZ2DMMhInput_withSingletVEVinPT");
}
