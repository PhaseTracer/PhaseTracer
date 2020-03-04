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

#ifndef PHASETRACER_POTENTIAL_PLOTTER_INCLUDED
#define PHASETRACER_POTENTIAL_PLOTTER_INCLUDED

#include <boost/filesystem.hpp>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>  

#include "logger.hpp"
#include "potential.hpp"

namespace PhaseTracer {

void potential_plotter(EffectivePotential::Potential &P,
                       double T,
                       std::string prefix = "model",
                       double xmin = -100.,
                       double xmax = 100.,
                       double dx = 0.1,
                       double ymin = -100.,
                       double ymax = 100.,
                       double dy = 0.1) {

  const auto dat_name = prefix + "_potential.dat";
  std::ofstream output_file;
  output_file.open(dat_name);
  output_file << std::fixed << std::setprecision(20);

  Eigen::VectorXd X(2);

  output_file << "#\t x\t y\t V" <<std::endl;
  for (double xx = xmin; xx < xmax; xx += dx) {
    for (double yy = ymin; yy < ymax; yy += dy) {
      X << xx, yy;
      output_file << xx << "\t" << yy << "\t"
                  << P.V(X, T)
                  << std::endl;
      }
  }
    
    output_file.close();

  const boost::filesystem::path this_file(__FILE__);
  const auto this_dir = this_file.parent_path();
  const boost::filesystem::path file("../make_plots/potential_plotter.gnu");
  const auto program = this_dir / file;

  const auto png_name = prefix + "_potential.png";
  std::string vars = "png_name = '" + png_name + "'; dat_name = '" + dat_name + "'; T = " + std::to_string(T) + ";";
  std::string command = "gnuplot -e \"" + vars + "\" " + program.string();

  LOG(debug) << "Executing " << command;
  const auto result = std::system(command.c_str());
}

}  // namespace PhaseTracer

#endif
