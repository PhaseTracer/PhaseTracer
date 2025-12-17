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

#ifndef PHASETRACER_SERIALIZE_HPP_
#define PHASETRACER_SERIALIZE_HPP_

#include <string>
#include <sstream>

#include "transition_finder.hpp"


namespace PhaseTracer {

std::stringstream serialize(const PhaseTracer::TransitionFinder &tf) {

  std::stringstream ss;
  ss << std::setprecision(10);

  auto phases = tf.get_phases();
  auto transitions = tf.get_transitions();
  auto transition_paths = tf.get_transition_paths();

  for (auto &p : phases) {
    ss << "# phase " << p.key << "\n";
    for (int i = 0; i < p.X.size(); i++) {
      ss << p.T[i] << " " << p.V[i];
      for (int j = 0; j < p.X[i].size(); j++) {
        ss << " " << p.X[i][j];
      }
      ss << "\n";
    }
    ss << "\n";
  }

  for (auto &t : transitions) {
    ss << "# transition\n"
       << t.false_phase.key << " " << t.true_phase.key << " " << t.TC;

    for (int j = 0; j < t.true_vacuum.size(); j++) {
      ss << " " << t.true_vacuum[j];
    }
    
    for (int j = 0; j < t.false_vacuum.size(); j++) {
      ss << " " << t.false_vacuum[j];
    }
    
    ss << " " << t.key
       << " " << t.id
       << " " << int(t.subcritical)
       << "\n\n";
  }

  for (auto &tp : transition_paths) {

    if (tp.transitions.size() == 0) {
      // A negative number represents a phase rather than a transition.
      ss << "# transition-path\n"
         << "-" << tp.phases[0] << "\n\n";
      continue;
    }

    ss << "# transition-path\n"
       << tp.transitions[0].transitionIndex;
       
    for (int i = 1; i < tp.transitions.size(); ++i) {
      ss << " " << tp.transitions[i].transitionIndex;
    }
    
    ss << "\n\n";
  }

  return ss;
}

} // namespace PhaseTracer

#endif // PHASETRACER_SERIALIZE_HPP_
