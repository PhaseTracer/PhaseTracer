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

#ifndef PHASETRACER_PHASE_PLOTTER_HPP_
#define PHASETRACER_PHASE_PLOTTER_HPP_

#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#include <boost/filesystem.hpp>

#include "logger.hpp"
#include "transition_finder.hpp"
#include "serialize.hpp"

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>

namespace PhaseTracer {

void phase_plotter(const PhaseTracer::TransitionFinder &tf, std::string folder_name = "", std::string file_basename = "model", bool make_plot = true) {

  std::ofstream output_file;

  if (folder_name == "") {
    output_file.open(file_basename + ".dat");
  } else {
    boost::filesystem::create_directories(folder_name);
    output_file.open(folder_name + "/" + file_basename + ".dat");
  }

  if (!output_file) {
    throw std::runtime_error("Could not open file for writing");
  }

  output_file << serialize(tf).str();
  output_file.close();

  if (make_plot) {
    const boost::filesystem::path this_file(__FILE__);
    const auto this_dir = this_file.parent_path();
    const boost::filesystem::path file("../make_plots/phase_plotter.py");
    const auto program = this_dir / file;
    std::string command = "python3 -W ignore \"" + program.string() + "\" \"" + folder_name + "\" \"" + file_basename + "\"";

    LOG(debug) << "Executing " << command;
    const auto result = std::system(command.c_str());
  }
}

} // namespace PhaseTracer

#endif // PHASETRACER_PHASE_PLOTTER_HPP_
