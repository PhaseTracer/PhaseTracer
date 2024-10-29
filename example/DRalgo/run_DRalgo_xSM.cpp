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

#include "DRalgo_xSM.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "action_calculator.hpp"
#include "gravwave_calculator.hpp"
#include "logger.hpp"

std::string toString(std::vector<double> in, std::vector<double> out) {
  std::stringstream data_str;
  for (auto i : in    ) data_str << i << "\t";
  for (auto i : out   ) data_str << i << "\t";
  return data_str.str();
}

int main( int argc, char* argv[]) {

  std::ofstream output_file;
  output_file.open("output.txt");

  LOGGER(fatal)

  double Ms, lambda_s, lambda_hs;
  bool running;
  int potential, matching;

   if ( argc == 1 ) {
    Ms = 135;
    lambda_s = 1;
    lambda_hs =  1.05;
    running = true;
    potential = 1;
    matching = 1;
  }
   else if ( argc == 4 ) {
    Ms = atof(argv[1]);
    lambda_s = atof(argv[2]);
    lambda_hs =  atof(argv[3]);
    running = true;
    potential = 1;
    matching = 1;
  } else if ( argc == 7 ) {
    Ms = atof(argv[1]);
    lambda_s = atof(argv[2]);
    lambda_hs =  atof(argv[3]);
    // Takes either 0/1 or true/false inputs for running
    running = (std::string(argv[4]) == "1" || std::string(argv[4]) == "true");
    potential = atoi(argv[5]);
    if ( potential < 0 || potential > 2 ) {
      std::cerr << "Error: Incorrect input for potential. Use 0, 1, or 2." << std::endl;
      return 0;
    }
    matching = atoi(argv[6]);
    if ( matching < 0 || matching > 2 ) {
      std::cerr << "Error: Incorrect input for matching. Use 0, 1, or 2." << std::endl;
      return 0;
    }
  } else {
    std::cout << "Incorrect input. Use 'run_DRalgo_xsm <Ms> <lambda_s> <lambda_hs> <running> <potential> <matching>'." << std::endl;
    return 0;
  }

  // Initialise Model
  EffectivePotential::DR_xsm model(Ms, lambda_hs, lambda_s);
  model.set_matching_flag( matching );
  model.set_potential_flag( potential );
  model.set_running_flag( running );
  std::vector<double> in = {Ms, lambda_s, lambda_hs};

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
  
  pf.set_seed(0);
  pf.set_check_hessian_singular(false);
  pf.set_check_vacuum_at_high(false);
  pf.set_t_high(1000);
  pf.set_t_low(15);

  try {
    pf.find_phases();
  } catch (...) {
    std::vector<double> out = {-1, 0, 0, 0, 0, 0, 0, 0};
    output_file << toString(in, out) << std::endl;
    return 0;
  }
  std::cout << pf;

  // Check there's multiple phases
  auto p1 = pf.get_phases();
  if ( p1.size() == 0 ) { return 0; }

  // Make ActionCalculator object
  PhaseTracer::ActionCalculator ac(model);
  
  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(pf, ac);
  tf.find_transitions();

  std::cout << tf;

  auto t = tf.get_transitions();
  if (t.size()==0){
    std::vector<double> out = {-2, 0, 0, 0, 0, 0, 0, 0};
    output_file << toString(in, out) << std::endl;
    return 0;
  }

  if ( isnan(t[0].TN) ){
    std::vector<double> out = {-3, t[0].TC, t[0].TN, 0, 0, 0, 0, };
    output_file << toString(in, out) << std::endl;
    output_file.close();
    return 0;
  }

  // Make GravWave Object
  PhaseTracer::GravWaveCalculator gc(tf);
  gc.set_min_frequency(1e-4);
  gc.set_max_frequency(1e+1);
  gc.calc_spectrums();

  std::cout << gc;

  auto gw = gc.get_spectrums();

  if ( isnan(gw[0].beta_H) || gw[0].beta_H < 1 || gw[0].beta_H > 1e100){
    std::vector<double> out = {-5, t[0].TC, t[0].TN, (t[0].TC - t[0].TN)/t[0].TC, gw[0].beta_H, 0, 0, 0};
    output_file << toString(in, out) << std::endl;
    output_file.close();
    return 0;
  }

  std::vector<double> out = {(float)t.size(), t[0].TC, t[0].TN, (t[0].TC - t[0].TN)/t[0].TC, gw[0].beta_H, gw[0].peak_amplitude, gw[0].peak_frequency, gw[0].SNR[0]};
  output_file << toString(in, out) << std::endl;
  output_file.close();
  return 0;

}