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

#ifndef PHASETRACER_PHASE_PLOTTER_INCLUDED
#define PHASETRACER_PHASE_PLOTTER_INCLUDED

#include <boost/filesystem.hpp>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>

#include "logger.hpp"
#include "transition_finder.hpp"

namespace PhaseTracer {

void phase_plotter(PhaseTracer::TransitionFinder tf, std::string prefix = "model") {
  std::ofstream output_file;
  output_file.open(prefix + ".dat");

  auto phases = tf.get_phases();
  auto transitions = tf.get_transitions();

  for (auto &p : phases) {
    output_file << "# phase" << std::endl << std::endl;
    for (int i = 0; i < p.X.size(); i++) {
      output_file << p.T[i] << " " << p.V[i];
      for (int j = 0; j < p.X[i].size(); j++) {
        output_file << " " << p.X[i][j];
      }
      output_file << std::endl;
    }
  }

  for (auto &t : transitions) {
    output_file << "# transition" << std::endl << std::endl;
    output_file << t.TC;
    for (int j = 0; j < t.true_vacuum.size(); j++) {
      output_file << " " << t.true_vacuum[j];
    }
    for (int j = 0; j < t.false_vacuum.size(); j++) {
      output_file << " " << t.false_vacuum[j];
    }
    output_file << " " << t.key << std::endl;
  }
  output_file.close();

  const boost::filesystem::path this_file(__FILE__);
  const auto this_dir = this_file.parent_path();
  const boost::filesystem::path file("../make_plots/phase_plotter.py");
  const auto program = this_dir / file;
  std::string command = "python -W ignore " + program.string() + " " + prefix;

  LOG(debug) << "Executing " << command;
  const auto result = std::system(command.c_str());
}

}  // namespace PhaseTracer

#endif
