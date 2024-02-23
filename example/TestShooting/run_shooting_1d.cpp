/**
 The scale example BubbleProfiler
 ./run_BP_scale 1. 0.6 200.
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"
#include "shooting.hpp"
#include "path_deformation.hpp"


class Example: public PhaseTracer::PotentialForShooting {
public:
  Example(double k1_, double k2_, double k3_) :
    k1(k1_), k2(k2_), k3(k3_) {}
  
  double V(double phi) const override {
    return k1*pow(phi,4) - k2*pow(phi,3) + k3*pow(phi,2);
//    return 0.25*pow(phi,4) - 0.4*pow(phi,3) + 0.1 * pow(phi,2);
  }

  
  double dV(double phi) const override {
    return 4.*k1*pow(phi,3) - 3.*k2*pow(phi,2) + 2.*k3*phi;
//    return phi*(phi-.2)*(phi-1);
  }
  
  double d2V(double phi) const override {
    return 12.*k1*pow(phi,2) - 6.*k2*phi + 2.*k3;
//    return (phi-.47)*(phi-1) + phi*(phi-1) + phis*(phi-.47);
//    return (phi-.2)*(phi-1) + phi*(phi-1) + phi*(phi-.2);
  }
  
private:
  double k1;
  double k2;
  double k3;
};

int main(int argc, char* argv[]) {

  LOGGER(debug);
  
//  Example p(0.25,0.49,0.235);
//  Example p(0.25,0.4,0.1);
//
//  PhaseTracer::Shooting s(p);
//
//  auto profile = s.findProfile(0,1);
  
//  std::ofstream file("test_data.txt");
//  for (int jj=0; jj< profile.R.size(); jj++){
//    std::cout << "R, phi, dphi = " << profile.R[jj] << ", " << profile.Phi(jj) << ", "  << profile.dPhi(jj) << std::endl;
//    file << profile.R[jj] << ", " << profile.Phi(jj) << ", "  << profile.dPhi(jj) << std::endl;
//  }
//  file.close();
//
//  auto action = s.calAction(profile);
//  std::cout << "action = " << action << std::endl;
//
  PhaseTracer::PathDeformation pd(0);
  auto a = pd.fullTunneling(0);
  
  
  
//  double E, alpha, scale;
//
//  if ( argc == 4 ) {
//    E = atof(argv[1]);
//    alpha = atof(argv[2]);
//    scale = atof(argv[3]);
//  } else {
//    std::cout << "Use ./run_BP_scale E alpha scale" << std::endl;
//    return 0;
//  }
//
//  LOGGER(debug);
//
//  // Construct our model
//  EffectivePotential::BP_scale model(E, alpha, scale);
//
//  // Make PhaseFinder object and find the phases
//  PhaseTracer::PhaseFinder pf(model);
////  pf.set_seed(0);
//  pf.set_check_hessian_singular(false);
//  pf.set_check_vacuum_at_high(false);
//  pf.set_guess_points({Eigen::VectorXd::Zero(1)});
//
//  pf.find_phases();
//  std::cout << pf;
//
//  // Make TransitionFinder object and find the transitions
//  PhaseTracer::TransitionFinder tf(pf);
////  tf.find_transitions();
////  std::cout << tf;
//
//  const auto phases = pf.get_phases();
//  auto phase1 = phases[1];
//  auto phase2 = phases[0];
//  auto s = tf.get_action(phase1, phase2, 0);
//  auto vacua = tf.get_vacua_at_T(phase1, phase2, 0);
//  std::cout << vacua[0][0] << " " << vacua[1][0] << " " << s << std::endl;
//
  return 0;
}
