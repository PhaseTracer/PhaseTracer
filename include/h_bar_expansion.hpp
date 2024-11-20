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

#ifndef PHASETRACER_H_BAR_EXPANSION_HPP_
#define PHASETRACER_H_BAR_EXPANSION_HPP_

#include <vector>

#include "phase_finder.hpp"
#include "one_loop_potential.hpp"

namespace PhaseTracer {

/** @brief Gauge-invariant h-bar expansion method */
class HbarExpansion : public PhaseFinder {
 public:
  /** @brief H-bar expansion needs a one-loop potential */
  HbarExpansion(EffectivePotential::OneLoopPotential &potential) : PhaseFinder(potential), P1l(potential) {};
  void find_phases() override;
  Point phase_at_T(const Phase& phase, double T) const override;
  /**
   * @brief h-bar expansion may consider one minima and one turning point or maxima
   *
   * Refer to latter as pseudo-phase and must add it manually.
   */
  void add_pseudo_phase(Eigen::ArrayXd pseudo_phase) { pseudo_phases.push_back(pseudo_phase); }

 private:
  std::function<double(Eigen::VectorXd)> make_objective(double T) const override;
  std::vector<Eigen::ArrayXd> pseudo_phases;
  EffectivePotential::OneLoopPotential &P1l;
};

/** @brief Trace high-temperature expansion of potential */
class HTExpansion : public PhaseFinder {
 public:
  /** @brief High-temperature expansion needs a one-loop potential */
  HTExpansion(EffectivePotential::OneLoopPotential &potential) : PhaseFinder(potential), P1l(potential) {};
 private:
  std::function<double(Eigen::VectorXd)> make_objective(double T) const override;
  EffectivePotential::OneLoopPotential &P1l;
};

}  // namespace PhaseTracer

#endif  // PHASETRACER_H_BAR_EXPANSION_HPP_
