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
#include "transition_graph_util.hpp"

namespace TransitionGraph
{
    struct Path;
}

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
  size_t key;
  bool subcritical;
  size_t id;

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
      o << "false vacuum = " << a.false_vacuum << std::endl
        << "true vacuum = " << a.true_vacuum << std::endl
        << "changed = " << a.changed << std::endl
        << "TC = " << a.TC << std::endl
        << "gamma = " << a.gamma << std::endl
        << "delta potential = " << a.delta_potential << std::endl;
      // TODO: Ideally we would only print this property if we check for subcritical transitions.
      // TODO: Unfortunately this is a property of the TransitionFinder and the Transition doesn't have knowledge of
      // TODO: this. We could store it here but that seems wasteful.
      //if(check_subcritical_transitions)
      {
        o << "subcritical = " << a.subcritical << std::endl;
      }
    } else {
      o << "=== failure. message =  " << a.message << " ===" << std::endl;
    }
    return o;
  }

};

class TransitionFinder {
 public:
  explicit TransitionFinder(PhaseFinder& pf_) : pf(pf_) {}
  virtual ~TransitionFinder() = default;

  /** Find all transitions between all phases */
  void find_transitions();

  /** Called from find_transitions; adds subcritical transitions to the transitions vector */
  void append_subcritical_transitions();

  /** Called from append_subcritical_transitions; checks whether there is a subcritical transition between two phases. */
  bool checkSubcriticalTransition(const std::vector<PhaseTracer::Phase>& phases, int i, int j, double Tmax,
    double energyAtTmax, bool checkFromNewPhase, std::vector<bool>& isTransitionedTo);

  /** Called from find_transitions; checks whether the transition should be kept or rejected. */
  bool validateTransition(const Transition& transition) const;

  /** Find all transition paths  */
  void find_transition_paths(const EffectivePotential::Potential& model, bool knownHighTPhase);

  /** Retrieve all transitions between all phases */
  std::vector<Transition> get_transitions() const { return transitions; }

  /** Retrieve all transition paths */
  std::vector<TransitionGraph::Path> get_transition_paths() const { return transition_paths; }

  /** Retrieve all phases */
  std::vector<Phase> get_phases() const { return pf.get_phases(); }

  /** Pretty-printer for set of transitions in this object */
  friend std::ostream& operator << (std::ostream& o, const TransitionFinder& a);

 private:
  /** Object with phases and potential */
  PhaseFinder &pf;

  /** Whether already calculated all transitions */
  bool calculated_transitions = false;

  /** Container for all transitions between any two phases */
  std::vector<Transition> transitions;

  /** Container for all transition paths. */
  std::vector<TransitionGraph::Path> transition_paths;

  /** Find transitions between two phases between two temperatures */
  std::vector<Transition> find_transition(Phase p1, Phase p2, double T1, double T2, size_t currentID) const;

  /** Find many transitions between two phases at a particular resolution */
  std::vector<Transition> divide_and_find_transition(const Phase& phase1, const Phase& phase2, double T1, double T2, size_t currentID) const;

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
  /**
    * Whether we should check for subcritical transitions, i.e. transitions between phases when one phase is strictly of
    * lower energy than the other. This can occur when a new phase appears near a phase that is not of the highest energy
    * at that temperature.
    */
  PROPERTY(bool, check_subcritical_transitions, false)
};

}  // namespace PhaseTracer

#endif
