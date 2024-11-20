// ====================================================================
// This file is part of PhaseTracer

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ====================================================================

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "DRalgo_ah.hpp" 
#include "phasetracer.hpp"

std::string toString(std::vector<double> in, std::vector<double> out) {
  std::stringstream data_str;
  for (auto i : in    ) data_str << i << "\t";
  for (auto i : out   ) data_str << i << "\t";
  return data_str.str();;
}

int main(int argc, char* argv[]) {

  std::ofstream output_file;
  output_file.open("output.txt");

  bool debug_mode = false;
  double M, gsq, lam;
  int potential, matching;
  
  if ( argc == 1 ) {
    debug_mode = true;
    M = 100.;
    gsq = 0.42;
    lam =  0.005;
    potential = 1;
    matching = 1;
  } else if ( argc >= 6 ) {
    M = atof(argv[1]);
    gsq = atof(argv[2]);
    lam = atof(argv[3]);
    potential = atoi(argv[4]);
    matching = atoi(argv[5]);
  } else {
    std::cout << "Use ./run_DRalgo_ah <M> <gsq> <lam> <potential> <running>" << std::endl;
    return 0;
  }
  
  if (debug_mode){
    LOGGER(debug);
    std::cout << "M = " << M << std::endl
              << "gsq = " << gsq << std::endl
              << "lam = " << lam << std::endl;
  } else {
    LOGGER(fatal);
  }
    
  // Construct our model
  EffectivePotential::DR_ah model(gsq, M, lam);
  model.set_matching_flag(matching);
  model.set_potential_flag(potential);
  std::vector<double> in = {gsq, M, lam};

  // This makes the plot data found in fig1 of the PhaseTracer2 manual.
  // model.plot_data();

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
  
  pf.set_seed(0);
  pf.set_check_hessian_singular(false);
  pf.set_t_low(50);
  pf.set_t_high(1000);
  pf.set_lower_bounds({-1e5});
  pf.set_upper_bounds({1e5});
  
  try {
    pf.find_phases();
  } catch (...) {
    std::cout << "M = " << M << ",\t"
              << "gsq = " << gsq << ",\t"
              << "lam = " << lam << "\t"
              << "encounters bug!" << std::endl;
    std::vector<double> out = {-1, 0, 0, 0, 0, 0};
    output_file << toString(in, out) << std::endl;
    return 0;
  }
  std::cout << pf;

  // Make ActionCalculator object
  PhaseTracer::ActionCalculator ac(model);

  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(pf, ac);
  tf.find_transitions();

  std::cout << tf;

  auto t = tf.get_transitions();
  if (t.size()==0){
    std::cout << "M = " << M << ",\t"
              << "gsq = " << gsq << ",\t"
              << "lam = " << lam << "\t"
              << "found 0 transition!" << std::endl;
    std::vector<double> out = {-2, 0, 0, 0, 0, 0};
    output_file << toString(in, out) << std::endl;
    return 0;
  }

  PhaseTracer::GravWaveCalculator gc(tf);
  gc.set_min_frequency(1e-4);
  gc.set_max_frequency(1e+1);
  gc.calc_spectrums();

  std::cout << gc;

  auto gw = gc.get_spectrums();
  
  std::vector<double> out = {(float)t.size(), t[0].TC, t[0].true_vacuum[0], t[0].false_vacuum[0]};

  output_file << toString(in, out) << std::endl;
  output_file.close();
  return 0;
}
