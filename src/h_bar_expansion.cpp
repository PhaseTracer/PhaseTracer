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

#include "h_bar_expansion.hpp"
#include "one_loop_potential.hpp"


namespace PhaseTracer {

std::function<double(Eigen::VectorXd)> HTExpansion::make_objective(double T) const {
  // make objective tracing high-temperature expansion of potential
  std::function<double(Eigen::VectorXd)> objective = [this, T](Eigen::VectorXd x) {
    return P1l.VHT(x, T);
  };

  return objective;
}

std::function<double(Eigen::VectorXd)> HbarExpansion::make_objective(double /* T */) const {
  // make objective tracing tree-level part of potential
  std::function<double(Eigen::VectorXd)> objective = [this](Eigen::VectorXd x) {
    return P1l.V0(x);
  };

  return objective;
}

Point HbarExpansion::phase_at_T(const Phase& phase, double T) const {
  // ensure that we see full potential rather than tree-level part that
  // was traced
  Point new_;
  new_.x = phase.X.back();
  new_.t = T;
  new_.potential = P(new_.x, new_.t);  // full potential
  return new_;
}

void HbarExpansion::find_phases() {
  // tree-level minima - this does not depend on temperature
  const std::vector<Point> minima = get_minima_at_t_low();
  int key = 0;
  const auto zero = Eigen::VectorXd::Zero(n_scalars);

  for (const auto& m : minima) {
    Phase phase;
    phase.key = key++;
    phase.T = {t_low, t_high};
    phase.X = {m.x, m.x};
    phase.dXdT = {zero, zero};
    phase.V = {P.V(m.x, t_low), P.V(m.x, t_high)};
    phase.end_low = REACHED_T_STOP;
    phase.end_high = REACHED_T_STOP;
    phases.push_back(phase);
  }

  // add pseudo phases - might not correspond to a minima
  // of the potential at tree-level, but used for computation
  // for gauge invariance, should be a minima, maxima or saddle point

  for (const auto& p : pseudo_phases) {
    Phase phase;
    phase.key = key++;
    phase.T = {t_low, t_high};
    phase.X = {p, p};
    phase.dXdT = {zero, zero};
    phase.V = {P.V(p, t_low), P.V(p, t_high)};
    phase.end_low = REACHED_T_STOP;
    phase.end_high = REACHED_T_STOP;
    phases.push_back(phase);
  }
}

}  // namespace PhaseTracer
