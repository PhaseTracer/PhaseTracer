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

#ifndef PHASETRACER_CRITICAL_TEMPERATURES_HPP_INCLUDED
#define PHASETRACER_CRITICAL_TEMPERATURES_HPP_INCLUDED

#include <ostream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <boost/cstdint.hpp>

#include "phase_finder.hpp"
#include "overload.hpp"
#include "potential.hpp"
#include "action_calculator.hpp"

namespace PhaseTracer {

/** Information about root-finding for a transition */
enum Message {SUCCESS, NON_OVERLAPPING_T, ERROR};

struct Transition {
  /** Data about a particular transition */
  Message message;
  double TC;
  Phase true_phase;
  Phase false_phase;
  Eigen::VectorXd true_vacuum;
  Eigen::VectorXd false_vacuum;
  double gamma;
  std::vector<bool> changed;
  double delta_potential;
  double TN;
  size_t key;

  /** Pretty-printer for single transition */
  friend std::ostream& operator << (std::ostream& o, const Transition& a) {
    if (a.message == SUCCESS) {
      o << "=== transition from phase " << a.false_phase.key;
      if (a.key == 0) {
        o << " to phase " << a.true_phase.key << " ===" << std::endl;
      } else {
        o << " to symmetric partner " << a.key
          << " of phase " << a.true_phase.key << " ===" << std::endl;
      }
      o << "changed = " << a.changed << std::endl
        << "TC = " << a.TC << std::endl
        << "false vacuum (TC) = " << a.false_vacuum << std::endl
        << "true vacuum (TC) = " << a.true_vacuum << std::endl
        << "gamma (TC) = " << a.gamma << std::endl
        << "delta potential (TC) = " << a.delta_potential << std::endl
        << "TN = " << a.TN << std::endl;
    } else {
      o << "=== failure. message =  " << a.message << " ===" << std::endl;
    }
    return o;
  }

};

class TransitionFinder {
 public:
  explicit TransitionFinder(PhaseFinder& pf_) :
    pf(pf_), ac(pf_.P){
      calculate_action = false;
    }
  explicit TransitionFinder(PhaseFinder& pf_, ActionCalculator ac_) :
    pf(pf_), ac(ac_){
      calculate_action = true;
    }
  virtual ~TransitionFinder() = default;

  /** Find all transitions between all phases */
  void find_transitions();

  /** Retrieve all transitions between all phases */
  std::vector<Transition> get_transitions() const { return transitions; }

  std::vector<Eigen::VectorXd> get_vacua_at_T(Phase phase1, Phase phase2, double T, size_t i_unique=0)const;
  
  double get_action(Phase phase1, Phase phase2, double T, size_t i_unique=0) const;
  
  double get_action(Eigen::VectorXd vacuum_1, Eigen::VectorXd vacuum_2, double T) const;
  
  double get_Tnuc(Phase phase1, Phase phase2, size_t i_unique, double T_begin, double T_end) const;
  
  /** Retrieve all phases */
  std::vector<Phase> get_phases() const { return pf.get_phases(); }

  /** Pretty-printer for set of transitions in this object */
  friend std::ostream& operator << (std::ostream& o, const TransitionFinder& a);

  /** Object with phases and potential */
  PhaseFinder &pf;
  
 private:
  
  ActionCalculator ac;
  
  /** Whether already calculated all transitions */
  bool calculated_transitions = false;

  /** Container for all transitions between any two phases */
  std::vector<Transition> transitions;

  /** Find transitions between two phases between two temperatures */
  std::vector<Transition> find_transition(Phase p1, Phase p2, double T1, double T2) const;

  /** Find many transitions between two phases at a particular resolution */
  std::vector<Transition> divide_and_find_transition(const Phase& phase1, const Phase& phase2, double T1, double T2) const;

  /** Check whether two phase are overlapped at T*/
  bool phases_overlaped(const Phase& phase1, const Phase& phase2, double T) const;

  /** Find un-overlapped temperature region between the two phases*/
  std::vector<double> get_un_overlapped_T_range(const Phase& phase1, const Phase& phase2, double T1, double T2) const;

  /** Strength of phase transition for first N fields */
  double gamma(const Eigen::VectorXd& true_vacuum, const Eigen::VectorXd& false_vacuum, const double TC) const;

  /** Note which VEVs changed */
  std::vector<bool> changed(const Eigen::VectorXd& true_vacuum, const Eigen::VectorXd& false_vacuum) const;
  
  /** Number of scalar fields that could break electroweak symmetry */
  PROPERTY(int, n_ew_scalars, -1)

  /** Minimum separation between critical temperatures */
  PROPERTY(double, separation, 1.)
  /** Assume at most one critical temperature between two phases */
  PROPERTY(bool, assume_only_one_transition, true)
  /** Relative precision in critical temperature */
  PROPERTY(double, TC_tol_rel, 1.e-4)
  /** Maximum number of iterations when finding a critical temperature */
  PROPERTY(boost::uintmax_t, max_iter, 100)
  /** Relative tolerance for judging whether a field is changed during a transition */
  PROPERTY(double, change_rel_tol, 1.e-3)
  /** Absolute tolerance for judging whether a field is changed during a transition */
  PROPERTY(double, change_abs_tol, 1.e-3)
  PROPERTY(double, Tnuc_step, 1.)
  
  PROPERTY(bool, calculate_action, false)
  
};

}  // namespace PhaseTracer

#endif
