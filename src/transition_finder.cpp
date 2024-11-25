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

#include "logger.hpp"
#include "transition_finder.hpp"

#include <boost/math/tools/roots.hpp>

namespace PhaseTracer {

double TransitionFinder::find_critical_temperature(const Phase &phase1, const Phase &phase2, double T1, double T2) const {
  const auto f = [this, phase1, phase2](double T) { return this->pf.delta_potential_at_T(phase1, phase2, T); };
  const double root_bits = 1. - std::log2(TC_tol_rel);
  boost::math::tools::eps_tolerance<double> stop(root_bits);
  boost::uintmax_t non_const_max_iter = max_iter;
  const auto result = boost::math::tools::toms748_solve(f, T1, T2, stop, non_const_max_iter);
  return (result.first + result.second) * 0.5;
}

std::vector<Transition> TransitionFinder::find_transition(Phase phase1, Phase phase2, double T1, double T2, size_t currentID) const {
  if (T1 > T2) {
    LOG(debug) << "Phases do not overlap in temperature - no critical temperature";
    return {{NON_OVERLAPPING_T}};
  }

  double TC;

  try {
    TC = find_critical_temperature(phase1, phase2, T1, T2);
  } catch (const std::exception &e) {
    LOG(debug) << e.what() << " - probably no sign change between T = " << T1 << " and " << T2;
    return {{ERROR}};
  }

  LOG(debug) << "Found critical temperature = " << TC;

  const bool ordered = pf.delta_potential_at_T(phase1, phase2, T1) < 0.;
  if (!ordered) {
    std::swap(phase1, phase2);
  }

  const auto phase1_at_critical = pf.phase_at_T(phase1, TC);
  const auto phase2_at_critical = pf.phase_at_T(phase2, TC);
  const auto delta_potential = phase1_at_critical.potential - phase2_at_critical.potential;
  const auto true_vacua = pf.symmetric_partners(phase1_at_critical.x);
  const auto false_vacua = pf.symmetric_partners(phase2_at_critical.x);
  std::vector<Transition> unique_transitions;
  // Fix false_vacuum, and loop all possible symmetric partners of true_vacuum
  const auto false_vacuum = false_vacua[0];
  auto true_vacuum = true_vacua[0];
  for (size_t i_unique = 0; i_unique < true_vacua.size(); i_unique++) {
    true_vacuum = true_vacua[i_unique];
    bool duplicate = false;
    const auto fv = pf.symmetric_partners(false_vacuum);
    const auto tv = pf.symmetric_partners(true_vacuum);
    LOG(trace) << "False vacuum: " << false_vacuum;
    LOG(trace) << "and true vacuum: " << true_vacuum;
    for (const auto tran : unique_transitions) {
      for (size_t i = 0; i < fv.size(); ++i) {
        if ((fv[i] - tran.false_vacuum).norm() < change_abs_tol + change_rel_tol * fv[i].norm() && (tv[i] - tran.true_vacuum).norm() < change_abs_tol + change_rel_tol * tv[i].norm()) {
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
      unique_transitions.push_back({SUCCESS, TC, phase1, phase2, true_vacuum, false_vacuum, gamma_, changed_, delta_potential, std::numeric_limits<double>::quiet_NaN(), {}, {}, i_unique, false, currentID++});
    }
  }

  if (calculate_action) {
    size_t i_selected = 0;
    if (unique_transitions.size() > 1) {
      double min_action = std::numeric_limits<double>::max();
      for (size_t i_unique = 0; i_unique < unique_transitions.size(); i_unique++) {
        double Ttry = T1 + 0.9 * (TC - T1);
        double try_action = get_action(phase1, phase2, Ttry, i_unique);
        LOG(debug) << "Action at " << Ttry << " for transtion " << i_unique << " is " << try_action;
        if (try_action < min_action) {
          min_action = try_action;
          i_selected = i_unique;
        }
      }
      if (min_action < std::numeric_limits<double>::max()) {
        LOG(debug) << "Select the symmetric partner " << i_selected << ".";
      } else {
        LOG(debug) << "Can not select the symmetric partner.";
        // TODO
      }
    }

    double TN = get_Tnuc(phase1, phase2, i_selected, TC, T1);
    std::vector<Transition> selected_transition;
    selected_transition.push_back(unique_transitions[i_selected]);
    selected_transition[0].TN = TN;
    if (not std::isnan(TN)) {
      const auto vacua = get_vacua_at_T(phase1, phase2, TN, i_selected);
      selected_transition[0].true_vacuum_TN = vacua[0];
      selected_transition[0].false_vacuum_TN = vacua[1];
    }
    return selected_transition;
  } else {
    return unique_transitions;
  }
}

double TransitionFinder::get_Tnuc(const Phase &phase1, const Phase &phase2, size_t i_unique, double T_begin, double T_end) const {

  if (T_begin < T_end) {
    LOG(fatal) << "T_begin < T_end, so swith the values. ";
    T_begin = T_begin + T_end;
    T_end = T_begin - T_end;
    T_begin = T_begin - T_end;
  }

  LOG(debug) << "Find Tnuc between " << phase1.key << " and " << phase2.key << " in [" << T_begin << ", " << T_end << "]. ";

  const auto nucleation_criteria = [this, phase1, phase2, i_unique](double Ttry) {
    return this->get_action(phase1, phase2, Ttry, i_unique) / Ttry - 140.;
  };

  double Tnuc = std::numeric_limits<double>::quiet_NaN();
  // If action at T_begin is NaN, find the largest valid T_begin
  double nc = nucleation_criteria(T_begin);
  while (std::isnan(nc) or nc > 10000) {
    T_begin -= Tnuc_step;
    if (T_begin < T_end)
      return Tnuc;
    nc = nucleation_criteria(T_begin);
  }

  if (nc < 0) {
    LOG(debug) << "The tunneling possibility at T_begin satisfys the nucleation condition.";
    return T_begin;
  }

  nc = nucleation_criteria(T_end);
  // If the tunneling possibility at T_end is small, find T_end from T_begin
  if (nc > 0 or std::isnan(nc)) {
    double T_end_ = T_end;
    while (true) {
      T_end = T_begin - Tnuc_step;
      while (nucleation_criteria(T_end) < 0)
        break;
      while (T_end < T_end_)
        return Tnuc;
      T_begin = T_end;
    }
  }

  try {
    const double root_bits = 1. - std::log2(Tnuc_tol_rel);
    boost::math::tools::eps_tolerance<double> stop(root_bits);
    boost::uintmax_t non_const_max_iter = max_iter;
    const auto result = boost::math::tools::bisect(nucleation_criteria, T_end, T_begin, stop, non_const_max_iter);
    Tnuc = (result.first + result.second) * 0.5;
    LOG(debug) << "Found nucleation temperature = " << Tnuc;
  } catch (char *str) {
    std::cout << str << std::endl; // TODO
    // find the first
  }
  return Tnuc;
}

std::vector<Transition> TransitionFinder::divide_and_find_transition(const Phase &phase1, const Phase &phase2, double T1, double T2, size_t currentID) const {

  std::vector<Transition> roots;
  for (double T = T1; T <= T2; T += separation) {
    const auto bs = find_transition(phase1, phase2, T, std::min(T + separation, T2), currentID + roots.size());
    for (const auto &b : bs) {
      roots.push_back(b);
    }
  }
  return roots;
}

std::vector<Eigen::VectorXd> TransitionFinder::get_vacua_at_T(const Phase &phase1, const Phase &phase2, double T, size_t i_unique) const {
  const auto phase1_at_T = pf.phase_at_T(phase1, T);
  const auto phase2_at_T = pf.phase_at_T(phase2, T);
  const auto true_vacua_at_T = pf.symmetric_partners(phase1_at_T.x);
  const auto false_vacua_at_T = pf.symmetric_partners(phase2_at_T.x);

  return {false_vacua_at_T[0], true_vacua_at_T[i_unique]};
}

double TransitionFinder::get_action(const Eigen::VectorXd &vacuum_1, const Eigen::VectorXd &vacuum_2, double T) const {
  return ac.get_action(vacuum_1, vacuum_2, T);
}

double TransitionFinder::get_action(const Phase &phase1, const Phase &phase2, double T, size_t i_unique) const {
  const auto vacua = get_vacua_at_T(phase1, phase2, T, i_unique);
  return ac.get_action(vacua[0], vacua[1], T);
}

std::vector<double> TransitionFinder::get_action(const Phase &phase1, const Phase &phase2, std::vector<double> T_list, size_t i_unique) const {
  std::vector<double> action_list;
  for (const auto Ti : T_list) {
    double action;
    try {
      action = get_action(phase1, phase2, Ti, i_unique);
    } catch (char *str) {
      LOG(warning) << str << std::endl;
      action = std::numeric_limits<double>::max();
    }
    action_list.push_back(action);
  }
  return action_list;
}

void TransitionFinder::write_action_to_text(const Phase &phase1, const Phase &phase2, std::vector<double> T_list, const std::string &filename, size_t i_unique) const {
  std::vector<double> action_list = get_action(phase1, phase2, T_list, i_unique);
  std::ofstream outFile(filename);
  for (size_t i = 0; i < T_list.size(); ++i) {
    outFile << T_list[i] << " " << action_list[i] << std::endl;
  }
  outFile.close();
}

void TransitionFinder::write_action_to_text(const Transition &tran, double T_min, double T_max, size_t n_step, const std::string &filename, size_t i_unique) const {
  std::vector<double> T_list;
  for (double Ti = T_min; Ti <= T_max; Ti += (T_max - T_min) / n_step) {
    T_list.push_back(Ti);
  }
  auto phase1 = tran.true_phase;
  auto phase2 = tran.false_phase;
  write_action_to_text(phase1, phase2, T_list, filename, i_unique);
}

void TransitionFinder::write_action_to_text(const Transition &tran, const std::string &filename, size_t n_step, size_t i_unique) const {
  double T_min = std::max(tran.true_phase.T[0], tran.false_phase.T[0]);
  double T_max = tran.TC;
  write_action_to_text(tran, T_min, T_max, n_step, filename, i_unique);
}

double TransitionFinder::gamma(const Eigen::VectorXd &true_vacuum, const Eigen::VectorXd &false_vacuum, const double TC) const {
  const int b = true_vacuum.size() + 1;
  const int items = (b + (n_ew_scalars % b)) % b;
  return (true_vacuum - false_vacuum).head(items).norm() / TC;
}

std::vector<bool> TransitionFinder::changed(const Eigen::VectorXd &true_vacuum, const Eigen::VectorXd &false_vacuum) const {
  auto delta = true_vacuum - false_vacuum;
  std::vector<bool> changed_;
  for (unsigned int i = 0; i < delta.size(); i += 1) {
    changed_.push_back(std::abs(delta[i]) > change_abs_tol + change_rel_tol * std::max(true_vacuum.norm(), false_vacuum.norm()));
  }
  return changed_;
}

bool TransitionFinder::phases_overlaped(const Phase &phase1, const Phase &phase2, double T) const {
  return pf.identical_within_tol(pf.phase_at_T(phase1, T).x, pf.phase_at_T(phase2, T).x);
}

std::vector<double> TransitionFinder::get_un_overlapped_T_range(const Phase &phase1, const Phase &phase2, double T1, double T2) const {
  std::vector<double> T_range;

  if (!(phases_overlaped(phase1, phase2, T1) || phases_overlaped(phase1, phase2, T2))) {
    T_range.push_back(T1);
    T_range.push_back(T2);
    return T_range;
  }

  if (phases_overlaped(phase1, phase2, T1) && phases_overlaped(phase1, phase2, T2)) {
    LOG(fatal) << "Phases " << phase1.key << " and " << phase2.key << " are redundant";
    return T_range;
  }

  double T_same = phases_overlaped(phase1, phase2, T1) ? T1 : T2;
  double T_diff = phases_overlaped(phase1, phase2, T2) ? T1 : T2;
  T_range.push_back(T_diff);

  while (std::abs(T_diff - T_same) > TC_tol_rel) {
    double T_mid = 0.5 * (T_same + T_diff);
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

      if (T_range.size() == 0)
        continue;
      if (T_range.size() == 3) {
        const auto TC = T_range[2];
        const auto TN = T_range[2];
        LOG(debug) << "Phases " << phase1.key << " and " << phase2.key << " cross at T=" << TC;
        const auto phase1_at_critical = pf.phase_at_T(phase1, TC);
        const auto phase2_at_critical = pf.phase_at_T(phase2, TC);
        const double gamma_ = 0.;
        auto changed_ = changed(phase1_at_critical.x, phase2_at_critical.x);
        Transition f = {SUCCESS, TC, phase1, phase2, phase1_at_critical.x, phase2_at_critical.x, gamma_, changed_, phase1_at_critical.potential - phase2_at_critical.potential, TC, {}, {}, 0, false, transitions.size()};
        transitions.push_back(f);
      }

      const auto found = assume_only_one_transition ? find_transition(phase1, phase2, T_range[0], T_range[1], transitions.size()) : divide_and_find_transition(phase1, phase2, T_range[0], T_range[1], transitions.size());

      for (const auto &f : found) {
        if (f.message == SUCCESS) {
          transitions.push_back(f);
        }
      }
    }
  }

  if (check_subcritical_transitions) {
    LOG(debug) << "Checking for subcritical transitions...";
    append_subcritical_transitions();
  } else {
    LOG(debug) << "Ignoring subcritical transitions.";
  }

  // Finally sort transitions by critical temperature
  std::sort(transitions.begin(), transitions.end(),
            [](const Transition &a, const Transition &b) { return a.TC < b.TC; });
}

void TransitionFinder::append_subcritical_transitions() {
  const auto phases = pf.get_phases();

  std::vector<bool> isTransitionedTo(phases.size(), false);

  for (int i = 0; i < transitions.size(); ++i) {
    // Valid because we update the phase keys to their position in the phases list (at the end of merge_phase_gaps).
    isTransitionedTo[transitions[i].true_phase.key] = true;
  }

  // Loop through phases and check if any phase has lower energy than any other phase when it first appears (i.e. at
  // its Tmax).
  for (int i = 0; i < phases.size(); ++i) {
    double Tmax = phases[i].T.back();
    double energyAtTmax = phases[i].V.back();

    // Used to determine what phases need to be searched again to capture subcritical transitions *from* the new phase.
    int firstSubcriticalIndex = -1;

    for (int j = 0; j < phases.size(); ++j) {
      if (j == i) {
        continue;
      }

      LOG(debug) << "Checking subcritical transition " << phases[i].key << " -> " << phases[j].key;

      // The last parameter is whether we should consider transitions from the new phase. We want to check for this if
      // this new phase has lower energy than some other phase, because there is a subcritical transition to it so there
      // could be a transition from it, following that.
      if (firstSubcriticalIndex == -1 && checkSubcriticalTransition(phases, i, j, Tmax, energyAtTmax,
                                                                    firstSubcriticalIndex > -1, isTransitionedTo)) {
        firstSubcriticalIndex = j;
        LOG(debug) << "firstSubcriticalIndex j = " << firstSubcriticalIndex << " for phase i = " << i;
      }
    }

    // Search over all phases prior to the first subcritical transition to this phase. In the previous search we didn't
    // look for transitions from this phase in this range of phases.
    for (int j = 0; j < firstSubcriticalIndex; ++j) {
      if (j == i) {
        continue;
      }

      // This time we always want to check for a subcritical transition from this phase.
      checkSubcriticalTransition(phases, i, j, Tmax, energyAtTmax, true, isTransitionedTo);
    }
  }
}

bool TransitionFinder::checkSubcriticalTransition(const std::vector<PhaseTracer::Phase> &phases, int i, int j,
                                                  double Tmax, double energyAtTmax, bool checkFromNewPhase, std::vector<bool> &isTransitionedTo) {
  // If the phases do not overlap in temperature, there cannot be a transition between them.
  // if(abs(phasejAtTmax.t - Tmax) > 1)
  if (phases[j].T.back() < Tmax || phases[j].T.front() > Tmax) {
    LOG(debug) << "No overlap.";
    return false;
  }

  // Avoid double counting subcritical transitions between phases that appear at the same temperature.
  // if(i > j && phases[j].T.back() == Tmax)
  //{
  //  return false;
  //}

  PhaseTracer::Point phasejAtTmax = pf.phase_at_T(phases[j], Tmax);

  Eigen::VectorXd trueVacuum = phases[i].X.back();
  Eigen::VectorXd falseVacuum = phasejAtTmax.x;
  double dist = (trueVacuum - falseVacuum).norm();
  double gamma = (trueVacuum - falseVacuum).norm() / Tmax;
  std::vector<bool> vevChanged = changed(trueVacuum, falseVacuum);
  double deltaPotential = energyAtTmax - phasejAtTmax.potential;

  // If this second phase has a higher energy, then we have a subcritical transition from j -> i. Otherwise, we have
  // a subcritical transition from i -> j.
  // NOTE: changed on 16/11/2021 to code further below, because:
  //    - If phase i splits from phase j at T' (i.e. Tmax), with both phase continuing to exist below T', the precisely
  //        sampled potential at (phi_i, T') might *erroneously* be lower than the potential at (phi_j, T') which is
  //        linearly interpolated from surrounding samples in phase j.
  //    - If phase j was identified as the deeper phase, it should actually have a lower energy at T' but this may
  //        not be reflected after linear interpolation of potential energy samples.
  //    - So if phase i and phase j are very close together at T' (i.e. phase i's maximum temperature, Tmax), we should
  //        not add a subcritical transition j -> i. Instead, we should add (if deemed relevant) the subcritical
  //        transition i -> j.
  /*if (phasejAtTmax.potential > energyAtTmax)
  {
    LOG(debug) << "Adding subcritical transition " << j << " -> " << i;
    // TODO: need to handle symmetric transitions!
    //transitions.push_back({PhaseTracer::SUCCESS, -Tmax, phases[i], phases[j], trueVacuum, falseVacuum,
    //  -gamma, vevChanged, deltaPotential, transitions.size(), true});
    transitions.push_back({PhaseTracer::SUCCESS, Tmax, phases[i], phases[j], trueVacuum, falseVacuum,
      gamma, vevChanged, deltaPotential, 0, true, transitions.size()});

    return true;
  }
  else if(checkFromNewPhase)
  {
    LOG(debug) << "Adding subcritical transition " << i << " -> " << j;
    transitions.push_back({PhaseTracer::SUCCESS, Tmax, phases[j], phases[i], falseVacuum, trueVacuum,
      -gamma, vevChanged, -deltaPotential, 0, true, transitions.size()});
  }*/

  double distTol = 0.05 * pf.get_potential().get_field_scale();

  if (phasejAtTmax.potential > energyAtTmax) {
    LOG(debug) << "Correct energies for subcritical transition.";

    if (dist > distTol) {
      LOG(debug) << "Adding subcritical transition " << j << " -(" << Tmax << ")-> " << i;
      transitions.push_back({PhaseTracer::SUCCESS, Tmax, phases[i], phases[j], trueVacuum, falseVacuum, gamma, vevChanged, deltaPotential, Tmax, {}, {}, 0, true, transitions.size()});
      isTransitionedTo[i] = true;
    } else {
      LOG(debug) << "Not adding subcritical transition " << j << " -(" << Tmax << ")-> " << i << " from suspected phase splitting.";
    }
  } else {
    LOG(debug) << "Incorrect energies for subcritical transition: " << energyAtTmax << " vs. " << phasejAtTmax.potential;
  }

  if (isTransitionedTo[i] || checkFromNewPhase) {
    if (phasejAtTmax.potential < energyAtTmax || dist < distTol) {
      transitions.push_back({PhaseTracer::SUCCESS, Tmax, phases[j], phases[i], falseVacuum, trueVacuum, -gamma, vevChanged, -deltaPotential, Tmax, {}, {}, 0, true, transitions.size()});
      isTransitionedTo[j] = true;
    }
  }

  return checkFromNewPhase;
}

void TransitionFinder::find_transition_paths(const EffectivePotential::Potential &model, bool knownHighTPhase) {
  transition_paths = TransitionGraph::getPhaseHistory(pf, *this, model, knownHighTPhase);
}

std::ostream &operator<<(std::ostream &o, const TransitionFinder &a) {
  if (!a.calculated_transitions) {
    o << "no transitions calculated yet" << std::endl
      << std::endl;
  } else if (a.transitions.empty()) {
    o << "found no transitions" << std::endl
      << std::endl;
  } else {
    o << "found " << a.transitions.size() << " transition";
    if (a.transitions.size() > 1) {
      o << "s";
    }
    o << std::endl
      << std::endl;
    for (const auto &t : a.transitions) {
      o << t << std::endl;
    }
  }
  return o;
}
} // namespace PhaseTracer
