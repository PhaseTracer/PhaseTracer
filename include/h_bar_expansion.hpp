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

#ifndef PHASETRACER_H_BAR_EXPANSION_HPP_INCLUDED
#define PHASETRACER_H_BAR_EXPANSION_HPP_INCLUDED

#include "phase_finder.hpp"

namespace PhaseTracer {

class HbarExpansion : public PhaseFinder {
 using PhaseFinder::PhaseFinder;
 public:
  std::function<double(Eigen::VectorXd)> make_objective(double T) const override;
  void find_phases() override;
  virtual Point phase_at_T(const Phase& phase, double T) const override;
  void add_symmetric_phase(Eigen::ArrayXd symmetric_phase) { symmetric_phases.push_back(symmetric_phase); };
 private:
  std::vector<Eigen::ArrayXd> symmetric_phases;
};

}  // namespace PhaseTracer

#endif
