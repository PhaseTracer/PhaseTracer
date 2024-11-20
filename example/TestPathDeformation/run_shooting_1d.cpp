/**
 The 1d shooting example in CosmoTransitions
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"
#include "shooting.hpp"


class Example: public PhaseTracer::PotentialForShooting {
public:
  Example(double k1_, double k2_, double k3_) :
    k1(k1_), k2(k2_), k3(k3_) {}
  
  double V(double phi) const override {
    return k1*pow(phi,4) - k2*pow(phi,3) + k3*pow(phi,2);
  }

  
  double dV(double phi) const override {
    return 4.*k1*pow(phi,3) - 3.*k2*pow(phi,2) + 2.*k3*phi;
  }
  
  double d2V(double phi) const override {
    return 12.*k1*pow(phi,2) - 6.*k2*phi + 2.*k3;
  }
  
private:
  double k1;
  double k2;
  double k3;
};

int main(int argc, char* argv[]) {

  LOGGER(fatal);
//  LOGGER(debug);
  
  double k2_min = 0.4;  // Thick-walled
  double k3_min = 0.1;
  double k2_max = 0.49; // Thin-walled
  double k3_max = 0.235;
  
  int num_values = 100;
  std::vector<double> k2_values;
  std::vector<double> k3_values;
  double k1 = 0.25;

  std::ofstream file("data_shooting_1d.txt");
  
  for (int i = 0; i < num_values; ++i) {
    double k2 = k2_min + i * (k2_max - k2_min) / (num_values - 1);
    double k3 = k3_min + i * (k3_max - k3_min) / (num_values - 1);

    Example p(k1,k2,k3);
    PhaseTracer::Shooting s(p, 2);
    auto profile = s.findProfile(0,1);
    auto action = s.calAction(profile);

    std::cout << "k2= "<< k2 << ", action = " << std::setprecision(10) << action << std::endl;

    file << k2 << "\t" << std::setprecision(10) << action << std::endl;
  }
  
  file.close();
  return 0;
}
