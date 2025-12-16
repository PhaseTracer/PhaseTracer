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

#ifndef PHASETRACER_PHASE_FINDER_HPP_
#define PHASETRACER_PHASE_FINDER_HPP_

#include <ctime>
#include <ostream>
#include <vector>

#include <eigen3/Eigen/Core>
#include "nlopt.hpp"

#include "potential.hpp"
#include "overload.hpp"
#include "property.hpp"

namespace PhaseTracer {

/** Descriptors for possible locations of a minima */
enum minima_descriptor { ORIGIN,
                         CONSISTENT_VACUUM,
                         OTHER };

/** Descriptor for ends of phase */
enum end_descriptor { HIGH,
                      LOW,
                      BOTH };

/** Descriptors for ways of ending tracing a minimum */
enum phase_end_descriptor { REACHED_T_STOP,
                            FORBIDDEN_OR_BOUNDS,
                            HESSIAN_SINGULAR,
                            HESSIAN_NOT_POSITIVE_DEFINITE,
                            JUMP_INDICATED_END };
std::ostream &operator<<(std::ostream &o, const phase_end_descriptor &d);

struct Point {
  Eigen::VectorXd x;
  double potential;
  double t;

  /** Pretty-printer for single point */
  friend std::ostream &operator<<(std::ostream &o, const Point &p) {
    o << "x = " << p.x << ". T = " << p.t << ". Potential = " << p.potential;
    return o;
  }
};

struct Phase {
  size_t key;
  std::vector<Eigen::VectorXd> X;
  std::vector<double> T; // Ascending order
  std::vector<Eigen::VectorXd> dXdT;
  std::vector<double> V;
  bool redundant = false;
  phase_end_descriptor end_low;
  phase_end_descriptor end_high;

  /** Pretty-printer for single phase */
  friend std::ostream &operator<<(std::ostream &o, const Phase &p) {
    o << "=== phase key = " << p.key << " ===" << std::endl
      << "Maximum temperature = " << p.T.back() << std::endl
      << "Minimum temperature = " << p.T.front() << std::endl
      << "Field at tmax = " << p.X.back() << std::endl
      << "Field at tmin = " << p.X.front() << std::endl
      << "Potential at tmax = " << p.V.back() << std::endl
      << "Potential at tmin = " << p.V.front() << std::endl
      << "Ended at tmax = " << p.end_high << std::endl
      << "Ended at tmin = " << p.end_low << std::endl;
    return o;
  }

  bool contains_t(const double t) const {
    return t >= T.front() && t <= T.back();
  }
};

struct PhaseMerge {
  size_t fromPhase;
  size_t toPhase;
  double temperature;
  bool rejected = false;

  friend std::ostream &operator<<(std::ostream &o, const PhaseMerge &pm) {
    o << (pm.rejected ? "[REJECTED] " : "") << "Merge phases " << pm.fromPhase << " -> " << pm.toPhase << " at T = " << pm.temperature;
    return o;
  }

  friend bool operator<(const PhaseMerge &a, const PhaseMerge &b) {
    return a.temperature < b.temperature;
  }

  friend bool operator>(const PhaseMerge &a, const PhaseMerge &b) {
    return a.temperature > b.temperature;
  }
};

class PhaseFinder {
public:
  //! Find different phases as functions of temperature
  /*!
   * @param potential The total finite temperature effective potential object
   */
  explicit PhaseFinder(EffectivePotential::Potential &potential);
  virtual ~PhaseFinder() = default;

  /**
     For the test-points, use a list of guesses supplemented by points chosen randomly
     over a uniform interval in field space
  */
  std::vector<Eigen::VectorXd> generate_test_points() const;
  std::vector<Point> find_minima_at_t(double T) const;

  virtual void find_phases();

  double delta_potential_at_T(const Phase &phase1, const Phase &phase2, double T) const {
    LOG(debug) << "Delta potential";
    return phase_at_T(phase1, T).potential - phase_at_T(phase2, T).potential;
  }

  virtual Point phase_at_T(const Phase &phase, double T) const;

  /** Pretty-printer for a collection of phases */
  friend std::ostream &operator<<(std::ostream &o, const PhaseFinder &pf) {
    auto phases = pf.phases;

    o << "found " << phases.size() << " phase";
    if (phases.size() != 1) {
      o << "s";
    }
    o << std::endl
      << std::endl;

    for (auto &p : phases) {
      o << p << std::endl;
    }

    return o;
  }

  /** Set random seeds for algorithms for finding minima of potential */
  void set_seed(int _seed) {
    seed = _seed;
    if (seed >= 0) {
      std::srand(seed);
      nlopt::srand(seed);
    } else {
      std::srand(std::time(0));
      nlopt::srand_time();
    }
  }

  /** Generate symmetric partners for a point */
  std::vector<Eigen::VectorXd> symmetric_partners(const Eigen::VectorXd &a) const;

  /** Check that two minima are identical to within a particular tolerance */
  bool identical_within_tol(const Eigen::VectorXd &a, const Eigen::VectorXd &b) const;

  /** Check whether two phases overlap at T */
  bool phases_overlap(const Phase &phase1, const Phase &phase2, double T) const;

  /** return minima at T_low*/
  std::vector<Point> get_minima_at_t_low();
  /** return minima at T_high*/
  std::vector<Point> get_minima_at_t_high();

  /** Phases at T */
  std::vector<Phase> get_phases_at_T(double T);

  /** The deepest phase at T */
  Phase get_deepest_phase_at_T(double T);

  /** Allow the potential to be visible to other classes */
  const EffectivePotential::Potential &get_potential() const { return P; };

  /** The potential */
  EffectivePotential::Potential &P;

protected:
  /**
     Find local minima at a particular temperature. The overloads define
     different ways of passing an initial step.
  */
  virtual std::function<double(Eigen::VectorXd)> make_objective(double T) const;
  Point find_min(const Eigen::VectorXd &X, double T, Eigen::VectorXd step) const;
  Point find_min(const Eigen::VectorXd &X, double T, double step) const;
  Point find_min(const Eigen::VectorXd &X, double T) const;

  /** Check for a jump discontinuity between two phases*/
  bool jump(const Eigen::VectorXd &a, const Eigen::VectorXd &b) const;
  /** Get descriptor of a particular minimum */
  minima_descriptor get_minima_descriptor(const Point &minima) const;
  /**
    @overload
    Get descriptor of the deepest minimum of a set at a particular temperature
  */
  minima_descriptor get_minima_descriptor(const std::vector<Point> &minima, double T) const;
  /** Get deepest minimum of a set at a particular temperature */
  Point get_deepest_minima(const std::vector<Point> &minima, double T) const;
  /** Check whether the origin is the unique minima */
  bool origin_unique_minima(const std::vector<Point> &minima) const;

  static double wrap_nlopt(const std::vector<double> &x,
                           std::vector<double> &grad, void *data) {
    std::function<double(Eigen::VectorXd)> potential = *static_cast<std::function<double(Eigen::VectorXd)> *>(data);
    const auto phi = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(const_cast<double *>(x.data()), x.size());
    return potential(phi);
  }

  phase_end_descriptor trace_minimum(Point start, double tstop,
                                     double dtstart_, std::vector<Eigen::VectorXd> *X,
                                     std::vector<double> *T, std::vector<Eigen::VectorXd> *dXdT, std::vector<double> *V,
                                     Point *jumped) const;

  std::vector<Point> minima_at_t_low;
  std::vector<Point> minima_at_t_high;

  /** Expected change in minimum with temperature, dx/dt */
  Eigen::VectorXd dx_min_dt(const Eigen::VectorXd &X, double T) const;
  /** @overload Precomputed Hessian matrix */
  Eigen::VectorXd dx_min_dt(const Eigen::MatrixXd &hessian, const Eigen::VectorXd &X, double T) const;

  /** Check whether two phases are redundant */
  std::tuple<bool, bool> redundant(const Phase &phase1, const Phase &phase2, end_descriptor end = BOTH) const;

  /** Check whether point belongs to a known phase */
  bool belongs_known_phase(const Point &point) const;

  /** Electroweak vacuum at zero temperature */
  bool consistent_vacuum(const Eigen::VectorXd &x) const;

  /** Check whether point is out of boundary */
  bool out_of_bounds(const Eigen::VectorXd &x) const;

  /** Remove/combine identical phases */
  void remove_redundant();

  /** Split phases that overlap over some subset of their shared temperature range. */
  void split_overlapping_phases(Phase &phase1, Phase &phase2, bool isRedundantLow, bool isRedundantHigh);

  /** Merge phases separated by a negigible gap in field and temperature. */
  void merge_phase_gaps();

  bool should_merge_phases(const Phase &phase1, const Phase &phase2);

  int find_deepest_phase(const std::vector<PhaseMerge> &merges, const std::vector<int> &relevantMerges);

  void perform_phase_merge(const PhaseMerge &merge);

  /** Check whether Hessian is singular */
  bool hessian_singular(const Eigen::VectorXd &X, double T) const;
  /** @overload Precomputed Hessian matrix */
  bool hessian_singular(const Eigen::MatrixXd &hessian, const Eigen::VectorXd &X, double T) const;

  /** Check whether Hessian is positive definite  */
  bool hessian_positive_definite(const Eigen::VectorXd &X, double T) const;
  /** @overload Precomputed Hessian matrix */
  bool hessian_positive_definite(const Eigen::MatrixXd &hessian, const Eigen::VectorXd &X, double T) const;

  /** Default bound on fields */
protected:
  const double bound = 1600.;

  /** Absolute error below which field values are considered identical */
  PROPERTY(double, x_abs_identical, 1.)
  /** Relative error below which field values are considered identical */
  PROPERTY(double, x_rel_identical, 1.e-3)
  /** Absolute change in field that is considered a jump */
  PROPERTY(double, x_abs_jump, 0.5)
  /** Relative change in field that is considered a jump */
  PROPERTY(double, x_rel_jump, 1.e-2)
  /** Relative precision for finding minima */
  PROPERTY(double, find_min_x_tol_rel, 0.0001)
  /** Absolute precision for finding minima */
  PROPERTY(double, find_min_x_tol_abs, 0.0001)
  /** Algorithm for finding minima */
  PROPERTY(nlopt::algorithm, find_min_algorithm, nlopt::LN_SBPLX)
  /** Maximum number of function evaluations when finding minimum */
  PROPERTY(double, find_min_max_f_eval, 1000000)
  /** Minimum step size for finding minima*/
  PROPERTY(double, find_min_min_step, 1.e-4)
  /** Timeout for finding a minima */
  PROPERTY(double, find_min_max_time, 5.)
  /** Initial absolute step size when tracing a minimum */
  PROPERTY(double, find_min_trace_abs_step, 1.)
  /** Initial absolute step size when finding a minimum */
  PROPERTY(double, find_min_locate_abs_step, 1.)
  /** Number of guesses with which to find all minima before tracing them */
  PROPERTY(size_t, n_test_points, 100)
  /** Lower bounds on fields */
  PROPERTY(std::vector<double>, lower_bounds, {})
  /** Upper bounds on fields */
  PROPERTY(std::vector<double>, upper_bounds, {})
  /** Lowest temperature to consider */
  PROTECTED_PROPERTY(double, t_low, 0.)
  /** Highest temperature to consider */
  PROTECTED_PROPERTY(double, t_high, 1000.)
  /** The starting step-size relative to t_high - t_low */
  PROPERTY(double, dt_start_rel, 0.01)
  /**
    The minimum temperature interval (relative to t_high - t_low) with which to
    resolve the splitting of a phase.
    */
  PROPERTY(double, dt_min_rel_split_phase, 0.001)
  /**
     The jump in temperature from the end of one phase to the
     temperature at which we try to trace a new phase. If this is too
     large, intermediate phases may be skipped.
     This is relative to t_high - t_low.
  */
  PROPERTY(double, t_jump_rel, 0.005)
  /** The largest absolute step-size in temperature */
  PROPERTY(double, dt_max_abs, 50.)
  /** The largest relative step-size in temperature */
  PROPERTY(double, dt_max_rel, 0.25)
  /** The smallest relative step-size in temperature */
  PROPERTY(double, dt_min_rel, 1.e-7)
  /** The smallest absolute step-size in temperature */
  PROPERTY(double, dt_min_abs, 1.e-10)
  /** Container for the phases */
  PROTECTED_PROPERTY(std::vector<Phase>, phases, {})
  /** Number of scalar fields */
  PROTECTED_PROPERTY_CUSTOM_SETTER(size_t, n_scalars, 0)
  /** Number of scalar fields that could break electroweak symmetry */
  PROPERTY(size_t, n_ew_scalars, 0)
  /** Zero-temperature vacuum expectation value of electroweak charged scalars */
  PROPERTY(double, v, 246.)
  /** Check whether deepest vacuum at low temperature agrees with that expected */
  PROPERTY(bool, check_vacuum_at_low, true)
  /** Check whether there is a unique vacuum at high temperature */
  PROPERTY(bool, check_vacuum_at_high, true)
  /** Check whether dx_min_dt is valid - useful for 3dEFT */
  PROPERTY(bool, check_dx_min_dt, true)
  /** Tolerance for checking whether Hessian was singular */
  PROPERTY(double, hessian_singular_rel_tol, 1.e-2)
  /** Tolerance for checking solutions to linear algebra */
  PROPERTY(double, linear_algebra_rel_tol, 1.e-3)
  /** Random seed (if negative, use system clock) */
  PROPERTY_CUSTOM_SETTER(int, seed, -1)
  /** Whether to check singularity of the Hessian matrix */
  PROPERTY(bool, check_hessian_singular, true)
  /** Maximum number of iterations when tracing a minimum */
  PROPERTY(unsigned int, trace_max_iter, 100000)
  /** Minimum length of a phase in temperature */
  // Discard short phase checking, because it may cause endless loop.
  // PROPERTY(double, phase_min_length, 0.5)
  /** Guesses for locations of minima */
  PROPERTY(std::vector<Eigen::VectorXd>, guess_points, {})

  PROPERTY(bool, check_merge_phase_gaps, false)
  PROPERTY(double, dt_merge_phases, 5.)
  PROPERTY(double, dx_merge_phases, 70.)
};

} // namespace PhaseTracer

#endif // PHASETRACER_PHASE_FINDER_HPP_
