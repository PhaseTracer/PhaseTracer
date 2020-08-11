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

#include <cmath>
#include <boost/math/tools/roots.hpp>

#include "transition_finder.hpp"
#include "logger.hpp"

namespace PhaseTracer {

std::vector<Transition> TransitionFinder::find_transition(Phase phase1, Phase phase2, double T1, double T2) const {
  if (T1 > T2) {
    LOG(debug) << "Phases do not overlap in temperature - no critical temperature";
    return {{NON_OVERLAPPING_T}};
  }

  try {
    const auto f = [this, phase1, phase2](double T) {return this->pf.delta_potential_at_T(phase1, phase2, T);};

    const double root_bits = 1. - std::log2(TC_tol_rel);
    boost::math::tools::eps_tolerance<double> stop(root_bits);
    boost::uintmax_t non_const_max_iter = max_iter;
    const auto result = boost::math::tools::toms748_solve(f, T1, T2, stop, non_const_max_iter);
    const double TC = (result.first + result.second) * 0.5;

    LOG(debug) << "Found critical temperature = " << TC;

    const bool ordered = f(T1) < 0.;
    if (!ordered) {
      std::swap(phase1, phase2);
    }
    const auto phase1_at_critical = pf.phase_at_T(phase1, TC);
    const auto phase2_at_critical = pf.phase_at_T(phase2, TC);
    const auto delta_potential = phase1_at_critical.potential - phase2_at_critical.potential;
    const auto true_vacua = pf.symmetric_partners(phase1_at_critical.x);
    const auto false_vacua =  pf.symmetric_partners(phase2_at_critical.x);
    std::vector<Transition> unique_transitions;
    // Fix false_vacuum, and loop all possible symmetric partners of true_vacuum
    const auto false_vacuum = false_vacua[0];
    for (const auto true_vacuum : true_vacua) {
      bool duplicate = false;
      const auto fv = pf.symmetric_partners(false_vacuum);
      const auto tv = pf.symmetric_partners(true_vacuum);
      LOG(trace) << "False vacuum: " << false_vacuum;
      LOG(trace) << "and true vacuum: " << true_vacuum;
      for (const auto tran : unique_transitions) {
        for (size_t i=0; i< fv.size(); ++i) {
          if ((fv[i]-tran.false_vacuum).norm() < change_abs_tol + change_rel_tol * fv[i].norm()
            && (tv[i]-tran.true_vacuum).norm() < change_abs_tol + change_rel_tol * tv[i].norm()) {
            LOG(trace) << "are duplicate";
            duplicate = true;
            break;
           }
        }
      }
      if (!duplicate) {
        LOG(trace) << "are not duplicate";
        const auto gamma_ = gamma(true_vacuum, false_vacuum, TC);
        const auto changed_ = changed(true_vacuum, false_vacuum);
        unique_transitions.push_back({SUCCESS, TC, phase1, phase2,
          true_vacuum, false_vacuum, gamma_, changed_,
          delta_potential, unique_transitions.size()});
      }
    }
    return unique_transitions;
  } catch (const std::exception& e) {
    LOG(debug) << e.what() << " - probably no sign change between T = " << T1 << " and " << T2;
    return {{ERROR}};
  }
}

std::vector<Transition> TransitionFinder::divide_and_find_transition(const Phase& phase1, const Phase& phase2, double T1, double T2) const {
  std::vector<Transition> roots;
  for (double T = T1; T <= T2; T += separation) {
    const auto bs = find_transition(phase1, phase2, T, std::min(T + separation, T2));
    for (const auto &b : bs) {
      roots.push_back(b);
    }
  }
  return roots;
}

double TransitionFinder::gamma(const Eigen::VectorXd& true_vacuum, const Eigen::VectorXd& false_vacuum, const double TC) const {
  const int b = true_vacuum.size() + 1;
  const int items = (b + (n_ew_scalars % b)) % b;
  return (true_vacuum-false_vacuum).head(items).norm()/TC;
}

std::vector<bool> TransitionFinder::changed(const Eigen::VectorXd& true_vacuum, const Eigen::VectorXd& false_vacuum) const {
  auto delta = true_vacuum-false_vacuum;
  std::vector<bool> changed_;
  for (unsigned int i = 0; i < delta.size(); i += 1) {
    changed_.push_back(std::abs(delta[i])> change_abs_tol + change_rel_tol * std::max(true_vacuum.norm(), false_vacuum.norm()));
  }
  return changed_;
}

bool TransitionFinder::phases_overlaped(const Phase& phase1, const Phase& phase2, double T) const {
  return pf.identical_within_tol(pf.phase_at_T(phase1, T).x, pf.phase_at_T(phase2, T).x);
}

std::vector<double> TransitionFinder::get_un_overlapped_T_range(const Phase& phase1, const Phase& phase2, double T1, double T2) const {
  std::vector<double> T_range;

  if (!(phases_overlaped(phase1, phase2, T1) || phases_overlaped(phase1, phase2, T2))) {
    T_range.push_back(T1);
    T_range.push_back(T2);
    return T_range;
  }

  if ( phases_overlaped(phase1, phase2, T1) && phases_overlaped(phase1, phase2, T2) ) {
    LOG(fatal) << "Phases " << phase1.key << " and " << phase2.key << " are redundant";
    return T_range;
  }

  double T_same = phases_overlaped(phase1, phase2, T1) ? T1 : T2;
  double T_diff = phases_overlaped(phase1, phase2, T2) ? T1 : T2;
  T_range.push_back(T_diff);

  while (std::abs(T_diff-T_same) > TC_tol_rel) {
    double T_mid = 0.5*(T_same+T_diff);
    if (phases_overlaped(phase1, phase2, T_mid)) {
      T_same = T_mid;
    } else {
      T_diff = T_mid;
    }
  }
  T_range.push_back(T_diff);
  T_range.push_back(T_same);

  return T_range;
}

void TransitionFinder::find_transitions() {
  if (calculated_transitions) {
    return;
  }

  transitions.empty();
  calculated_transitions = true;
  const auto phases = pf.get_phases();

  for (const auto &phase1 : phases) {
    for (const auto &phase2 : phases) {
      // Don't investigate n1 -> n1 or n1 -> n2 as well as n2 -> n1
      if (phase1.key >= phase2.key) {
        continue;
      }

      LOG(debug) << "Finding critical temperatures between phases "
                 << phase1.key << " and " << phase2.key;


      double tmax = std::min(phase1.T.back(), phase2.T.back());
      double tmin = std::max(phase1.T.front(), phase2.T.front());

      if (tmin > tmax) {
        LOG(debug) << "Phases do not overlap in temperature - no critical temperature";
        continue;
      }

      auto T_range = get_un_overlapped_T_range(phase1, phase2, tmin, tmax);

      if (T_range.size() == 0) continue;
      if (T_range.size() == 3) {
        const auto TC = T_range[2];
        LOG(debug) << "Phases " << phase1.key << " and " << phase2.key << " cross at T=" << TC;
        const auto phase1_at_critical = pf.phase_at_T(phase1, TC);
        const auto phase2_at_critical = pf.phase_at_T(phase2, TC);
        const double gamma_ = 0.;
        auto changed_ = changed(phase1_at_critical.x, phase2_at_critical.x);
        Transition f = {SUCCESS, TC, phase1, phase2,
                        phase1_at_critical.x, phase2_at_critical.x,
                        gamma_, changed_,
                        phase1_at_critical.potential - phase2_at_critical.potential};
        transitions.push_back(f);
      }

      const auto found = assume_only_one_transition ?
        find_transition(phase1, phase2, T_range[0], T_range[1]) :
        divide_and_find_transition(phase1, phase2, T_range[0], T_range[1]);

      for (const auto &f : found) {
        if (f.message == SUCCESS) {
            transitions.push_back(f);
        }
      }
    }
  }

  // Finally sort transitions by critical temperature
  std::sort(transitions.begin(), transitions.end(),
    [](const Transition &a, const Transition &b) {return a.TC < b.TC;});
}

std::ostream& operator << (std::ostream& o, const TransitionFinder& a) {
  if (!a.calculated_transitions) {
    o << "no transitions calculated yet" << std::endl << std::endl;
  } else if (a.transitions.empty()) {
    o << "found no transitions" << std::endl << std::endl;
  } else {
    o << "found " << a.transitions.size() << " transition";
    if (a.transitions.size() > 1) {
      o << "s";
    }
    o << std::endl << std::endl;
    for (const auto &t : a.transitions) {
      o << t << std::endl;
    }
  }
  return o;
}
}  // namespace PhaseTracer
