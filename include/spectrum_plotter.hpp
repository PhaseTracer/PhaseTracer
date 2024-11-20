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

#ifndef PHASETRACER_SPECTRUM_PLOTTER_INCLUDED
#define PHASETRACER_SPECTRUM_PLOTTER_INCLUDED

#include <boost/filesystem.hpp>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>

#include "logger.hpp"
#include "transition_finder.hpp"
#include "gravwave_calculator.hpp"

namespace PhaseTracer {

void spectrum_plotter(PhaseTracer::GravWaveCalculator& gc, std::string prefix = "model", std::string repo = "") {

  auto spectrums = gc.get_spectrums();
  int no_spectrums = spectrums.size();

  for (int ii=0; ii < no_spectrums; ii++){
    gc.write_spectrum_to_text(spectrums[ii], repo + "GW_" + prefix + "_" + std::to_string(ii) + ".txt");
    double peak_freq = spectrums[ii].peak_frequency;

    const boost::filesystem::path this_file(__FILE__);
    const auto this_dir = this_file.parent_path();
    const boost::filesystem::path file("../make_plots/spectrum_plotter.py");
    const auto program = this_dir / file;
    std::string command = "python3 -W ignore " + program.string() + " " + prefix + " " + std::to_string(ii) + " " + std::to_string(peak_freq) + " " + repo;

    LOG(debug) << "Executing " << command;
    const auto result = std::system(command.c_str());
  }

}

}  // namespace PhaseTracer

#endif