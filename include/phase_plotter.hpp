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

#include <boost/filesystem.hpp>

#include "logger.hpp"
#include "transition_finder.hpp"

//#include <sys/stat.h>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>

namespace PhaseTracer {

void phase_plotter(const PhaseTracer::TransitionFinder& tf, std::string folderName = "", std::string fileName = "model", bool bPlot = true) {
  std::ofstream output_file;
  // Increase the output precision, otherwise some of the output phase temperatures are rounded to identical values. This causes issues with
  // the arrow_plot in phase_plotter.py.
  output_file << std::setprecision(10);

  if(folderName == "")
  {
    output_file.open(fileName + ".dat");
  }
  else
  {
    /*if (!boost::starts_with(folderName, "output/"))
    {
      folderName = "output/" + folderName;
    }*/
    //char* folderNameCharArray = &folderName[0];
    //int check = mkdir(folderNameCharArray, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    //std::cout << "check = " << check << std::endl;
    boost::filesystem::create_directories(folderName);
    output_file.open(folderName + "/" + fileName + ".dat");
  }

  if(!output_file)
  {
    std::cout << "Could not open file!" << std::endl;
    return;
  }

  auto phases = tf.get_phases();
  auto transitions = tf.get_transitions();
  auto transition_paths = tf.get_transition_paths();

  for (auto &p : phases) {
    output_file << "# phase " << p.key << std::endl;
    //output_file << p.key << " " << std::endl;
    for (int i = 0; i < p.X.size(); i++) {
      output_file << p.T[i] << " " << p.V[i];
      for (int j = 0; j < p.X[i].size(); j++) {
        output_file << " " << p.X[i][j];
      }
      output_file << std::endl;
    }
    output_file << std::endl;
  }

  for (auto &t : transitions) {
    output_file << "# transition" << std::endl;
    output_file << t.false_phase.key << " " << t.true_phase.key << " ";
    output_file << t.TC;
    for (int j = 0; j < t.true_vacuum.size(); j++) {
      output_file << " " << t.true_vacuum[j];
    }
    for (int j = 0; j < t.false_vacuum.size(); j++) {
      output_file << " " << t.false_vacuum[j];
    }
    output_file << " " << t.key;
    output_file << " " << t.id;
    output_file << " " << int(t.subcritical) << std::endl << std::endl;
  }

  //std::cout << "Writing " << transition_paths.size() << " transition paths to the output file..." << std::endl;

  for (auto &tp : transition_paths) {
    //std::cout << "Transition path: " << tp << std::endl;
    //std::cout << "Size of path: " << tp.transitions.size() << std::endl;
    
    if (tp.transitions.size() == 0) {
      // Prior to 25/05/21 we didn't print a transition path if no transitions were involved. However, it is still useful
      // to know which phase the Universe stayed in, directly from the transition path.
      //continue;

      output_file << "# transition-path" << std::endl;
      // A negative number represents a phase rather than a transition.
      output_file << "-" << tp.phases[0] << std::endl << std::endl;
      continue;
    }


    output_file << "# transition-path" << std::endl;
    output_file << tp.transitions[0].transitionIndex;
    for (int i = 1; i < tp.transitions.size(); ++i) {
      output_file << " " << tp.transitions[i].transitionIndex;
    }
    output_file << std::endl << std::endl;
  }
  output_file.close();

  if(!bPlot)
  {
    return;
  }

  const boost::filesystem::path this_file(__FILE__);
  const auto this_dir = this_file.parent_path();
  const boost::filesystem::path file("../make_plots/phase_plotter.py");
  const auto program = this_dir / file;
  std::string command = "python3 -W ignore \"" + program.string() + "\" \"" + folderName + "\" \"" + fileName + "\"";

  LOG(debug) << "Executing " << command;
  const auto result = std::system(command.c_str());
}

}  // namespace PhaseTracer

#endif  // PHASETRACER_PHASE_PLOTTER_HPP_