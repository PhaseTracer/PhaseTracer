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

std::vector<Transition> TransitionFinder::symmetric_partners(const Phase &phase1, const Phase &phase2, double TC, size_t currentID) const {

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

    for (const auto tran : unique_transitions) {
      for (size_t i = 0; i < fv.size(); ++i) {
        if ((fv[i] - tran.false_vacuum).norm() < change_abs_tol + change_rel_tol * fv[i].norm() && (tv[i] - tran.true_vacuum).norm() < change_abs_tol + change_rel_tol * tv[i].norm()) {
          LOG(trace) << "False vacuum: " << false_vacuum << " and true vacuum: " << true_vacuum << " are duplicate";
          duplicate = true;
          break;
        }
      }
    }
    if (!duplicate) {
      LOG(trace) << "False vacuum: " << false_vacuum << " and true vacuum: " << true_vacuum << " are not duplicate";
      const auto gamma_ = gamma(true_vacuum, false_vacuum, TC);
      const auto changed_ = changed(true_vacuum, false_vacuum);
      unique_transitions.push_back({SUCCESS, TC, phase1, phase2, true_vacuum, false_vacuum, gamma_, changed_, delta_potential, i_unique, currentID++});
    }
  }
  return unique_transitions;
}

std::vector<Transition> TransitionFinder::find_transition(const Phase &phase1, const Phase &phase2, double T1, double T2, size_t currentID) const {
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

  std::vector<Transition> unique_transitions = symmetric_partners(phase1, phase2, TC, currentID);

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
      if (std::isfinite(min_action)) {
        LOG(debug) << "Select the symmetric partner " << i_selected << ".";
      } else {
        LOG(debug) << "Cannot select the symmetric partner. min action = " << min_action;
        throw std::runtime_error("Non-finite action when selecting the symmetric partner");
      }
    }

    auto [TN, action_curve] = get_Tnuc(phase1, phase2, i_selected, TC, T1);
    if (!std::isnan(TN)) {
      const auto vacua = get_vacua_at_T(phase1, phase2, TN, i_selected);
      unique_transitions[i_selected].set_nucleation(TN, vacua[0], vacua[1]);
      unique_transitions[i_selected].set_action_curve(action_curve);
    }
    return {unique_transitions[i_selected]};
  } else {
    return unique_transitions;
  }
}

PolynomialFitterEigen TransitionFinder::get_action_curve(const Phase &phase1, const Phase &phase2, size_t i_unique, double TC, double T_end) const {

  PolynomialFitterEigen action_curve;

  // Get the nodes
  std::vector<double> T_list;
  double Tmax = TC - Tnuc_tol_rel;             // At TC action = infinity
  double Tmin = std::max(T_end, Tnuc_tol_rel); // We need S/T, so Tmin != 0
  if (Tmax <= Tmin) {
    return action_curve;
  }
  for (double Ti = Tmin; Ti <= Tmax; Ti += (Tmax - Tmin) / action_curve_nodes) {
    T_list.push_back(Ti);
  }
  std::vector<double> action_list = get_action(phase1, phase2, T_list, i_unique);

  // Filter the nodes
  std::vector<double> x;
  std::vector<double> y;
  x.reserve(action_list.size());
  y.reserve(action_list.size());
  for (size_t i = 0; i < action_list.size(); ++i) {
    double S3T = action_list[i] / T_list[i];
    if (S3T > 1 and S3T < 1000) {
      x.push_back(T_list[i]);
      y.push_back(action_list[i] * pow(TC - T_list[i], 2));
    }
  }
  x.shrink_to_fit();
  y.shrink_to_fit();

  // Perform fit
  action_curve.fit(x, y, TC, action_curve_order);
  action_curve.cal_MSE(x, y);

  std::ofstream output_file;
  output_file.open("action_curve_output.txt");
  output_file << T_end << "\t" << TC << "\t";
  Eigen::VectorXd ki = action_curve.getCoefficients();
  for (int i = 0; i < ki.size(); i++) {
    output_file << std::setprecision(15) << ki[i] << "\t";
  }
  output_file << action_curve.get_MSE() << std::endl;
  for (size_t i = 0; i < x.size(); ++i) {
    output_file << x[i] << "\t" << y[i] / pow(TC - x[i], 2) << std::endl;
  }
  output_file.close();

  return action_curve;
}

std::pair<double, PolynomialFitterEigen> TransitionFinder::get_Tnuc(const Phase &phase1, const Phase &phase2, size_t i_unique, double T_begin, double T_end) const {

  if (T_begin < T_end) {
    LOG(trace) << "T_begin < T_end, so swap the values";
    std::swap(T_begin, T_end);
  }

  LOG(debug) << "Find nucleation temperature between " << phase1.key << " and " << phase2.key << " in [" << T_begin << ", " << T_end << "]";

  PolynomialFitterEigen action_curve = get_action_curve(phase1, phase2, i_unique, T_begin, T_end);

  const auto nucleation_criteria = [this, phase1, phase2, i_unique, action_curve](double Ttry) {
    if (action_curve.get_success()) {
      return action_curve.get_S3T(Ttry) - 140.;
    } else {
      return this->get_action(phase1, phase2, Ttry, i_unique) / Ttry - 140.;
    }
  };

  // check criteria at start
  double nc = nucleation_criteria(T_begin);

  while (std::isnan(nc)) {
    LOG(debug) << "action was nan - stepping down in temperature";
    T_begin -= Tnuc_step;
    if (T_begin < T_end) {
      LOG(fatal) << "could not find nucleation temperature";
      return std::make_pair(std::numeric_limits<double>::quiet_NaN(), action_curve);
    }
    nc = nucleation_criteria(T_begin);
  }

  if (nc < 0) {
    LOG(debug) << "The tunneling possibility at T_begin satisfies the nucleation condition";
    return std::make_pair(T_begin, action_curve);
  }

  // check criteria at end

  nc = nucleation_criteria(T_end);

  if (nc > 0 || std::isnan(nc)) {
    LOG(debug) << "action was nan or no nucleation at end - stepping down in temperature";
    double T_end_ = T_end;
    while (true) {
      T_end = T_begin - Tnuc_step;
      if (T_end < T_end_) {
        LOG(fatal) << "could not find nucleation temperature";
        return std::make_pair(std::numeric_limits<double>::quiet_NaN(), action_curve);
      }
      if (nucleation_criteria(T_end) < 0) {
        break;
      }
      T_begin = T_end;
    }
  }

  // bisect for nucleation temperature

  const double root_bits = 1. - std::log2(Tnuc_tol_rel);
  boost::math::tools::eps_tolerance<double> stop(root_bits);
  boost::uintmax_t non_const_max_iter = max_iter;
  const auto result = boost::math::tools::toms748_solve(nucleation_criteria, T_end, T_begin, stop, non_const_max_iter);
  const double Tnuc = (result.first + result.second) * 0.5;
  LOG(debug) << "Found nucleation temperature = " << Tnuc;

  return std::make_pair(Tnuc, action_curve);
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
  return {phase2_at_T.x, true_vacua_at_T[i_unique]};
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

std::vector<double> TransitionFinder::get_un_overlapped_T_range(const Phase &phase1, const Phase &phase2, double T1, double T2) const {
  const bool o1 = pf.phases_overlap(phase1, phase2, T1);
  const bool o2 = pf.phases_overlap(phase1, phase2, T2);

  // no overlap - non-overlap is entire range

  if (!(o1 || o2)) {
    return {T1, T2};
  }

  // complete overlap - non-overlap is empty

  if (o1 && o2) {
    LOG(fatal) << "Phases " << phase1.key << " and " << phase2.key << " are redundant";
    return {};
  }

  // overlap at one end - return temperature at which distinct, and bracket for
  // point at which they overlap

  const double distinct = o1 ? T2 : T1;
  const double overlap = o1 ? T1 : T2;

  double a = distinct;
  double b = overlap;

  while (std::abs(a - b) > TC_tol_rel * 0.5 * (a + b)) {
    double T = 0.5 * (a + b);
    if (pf.phases_overlap(phase1, phase2, T)) {
      b = T;
    } else {
      a = T;
    }
  }

  return {distinct, a, b};
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

      LOG(debug) << "Finding transition between phases "
                 << phase1.key << " and " << phase2.key;

      double tmax = std::min(phase1.T.back(), phase2.T.back());
      double tmin = std::max(phase1.T.front(), phase2.T.front());

      if (tmin > tmax) {
        LOG(debug) << "Phases do not coexist in temperature - no transition";
        continue;
      }

      auto T_range = get_un_overlapped_T_range(phase1, phase2, tmin, tmax);

      if (T_range.size() == 0) {
        LOG(debug) << "Phases entirely overlap - no transition";
        continue;
      }

      if (T_range.size() == 3) {
        const double TC = T_range[2];
        LOG(debug) << "Phases join at T = " << TC;
        const auto phase_at_critical = pf.phase_at_T(phase1, TC);
        const std::vector<bool> changed_(pf.get_n_scalars(), false);
        const Transition f = {SUCCESS, TC, phase1, phase2, phase_at_critical.x, phase_at_critical.x, 0., changed_, 0., 0, transitions.size()};
        transitions.push_back(f);
      }

      LOG(debug) << "Phases co-exist - looking for transitions between " << T_range[0] << " and " << T_range[1];

      const bool ordered = pf.delta_potential_at_T(phase1, phase2, T_range[0]) < 0.;
      const Phase *tv = ordered ? &phase1 : &phase2;
      const Phase *fv = ordered ? &phase2 : &phase1;

      const std::vector<Transition> found = assume_only_one_transition ? find_transition(*tv, *fv, T_range[0], T_range[1], transitions.size()) : divide_and_find_transition(phase1, phase2, T_range[0], T_range[1], transitions.size());

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
      Transition t{SUCCESS, Tmax, phases[i], phases[j], trueVacuum, falseVacuum, gamma, vevChanged, deltaPotential, 0, transitions.size()};
      t.set_subcritical(true);
      transitions.push_back(t);
      isTransitionedTo[i] = true;
    } else {
      LOG(debug) << "Not adding subcritical transition " << j << " -(" << Tmax << ")-> " << i << " from suspected phase splitting.";
    }
  } else {
    LOG(debug) << "Incorrect energies for subcritical transition: " << energyAtTmax << " vs. " << phasejAtTmax.potential;
  }

  if (isTransitionedTo[i] || checkFromNewPhase) {
    if (phasejAtTmax.potential < energyAtTmax || dist < distTol) {
      Transition t{SUCCESS, Tmax, phases[j], phases[i], falseVacuum, trueVacuum, -gamma, vevChanged, -deltaPotential, 0, transitions.size()};
      t.set_subcritical(true);
      transitions.push_back(t);
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
    return o;
  }

  if (a.transitions.empty()) {
    o << "found no transitions";

    if (a.check_subcritical_transitions) {
      o << " and no subcritical transitions";
    } else {
      o << ", but did not check for subcritical transitions";
    }

    o << std::endl
      << std::endl;

    return o;
  }

  o << "found " << a.transitions.size() << " transition";
  if (a.transitions.size() > 1) {
    o << "s";
  }

  o << std::endl
    << std::endl;

  if (a.check_subcritical_transitions) {
    o << "checked for subcritical transitions";
  } else {
    o << "did not check for subcritical transitions";
  }

  o << std::endl
    << std::endl;

  for (const auto &t : a.transitions) {
    o << t << std::endl;
  }

  return o;
}
} // namespace PhaseTracer
