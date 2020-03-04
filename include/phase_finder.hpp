#ifndef PHASETRACER_PHASE_FINDER_HPP_INCLUDED
#define PHASETRACER_PHASE_FINDER_HPP_INCLUDED

#include <ostream>
#include <vector>
#include <Eigen/Core>
#include <ctime>

#include "nlopt.hpp"
#include "potential.hpp"
#include "overload.hpp"
#include "property.hpp"


namespace PhaseTracer {

/** Descriptors for possible locations of a minima */
enum minima_descriptor {ORIGIN, CONSISTENT_VACUUM, OTHER};

/** Descriptor for ends of phase */
enum end_descriptor {HIGH, LOW, BOTH};

/** Descriptors for ways of ending tracing a minimum */
enum phase_end_descriptor {REACHED_T_STOP, FORBIDDEN_OR_BOUNDS, HESSIAN_SINGULAR, HESSIAN_NOT_POSITIVE_DEFINITE, JUMP_INDICATED_END};
std::ostream& operator << (std::ostream& o, const phase_end_descriptor& d);

struct Point {
  Eigen::VectorXd x;
  double potential;
  double t;

  /** Pretty-printer for single point */
  friend std::ostream& operator << (std::ostream& o, const Point& p) {
    o << "x = " << p.x << ". T = " << p.t << ". Potential = " << p.potential;
    return o;
  }
};

struct Phase {
  size_t key;
  std::vector<Eigen::VectorXd> X;
  std::vector<double> T;  // Ascending order
  std::vector<Eigen::VectorXd> dXdT;
  std::vector<double> V;
  bool redundant = false;
  phase_end_descriptor end_low;
  phase_end_descriptor end_high;

  /** Pretty-printer for single phase */
  friend std::ostream& operator << (std::ostream& o, const Phase& p) {
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

class PhaseFinder {
 public:
  //! Find different phases as functions of temperature
  /*!
   * @param potential The total finite temperature effective potential object
   */
  PhaseFinder(EffectivePotential::Potential &potential);
  virtual ~PhaseFinder() = default;

  /**
     For the test-points, use a list of guesses supplemented by points chosen randomly
     over a uniform interval in field space
  */
  std::vector<Eigen::VectorXd> generate_test_points() const;
  std::vector<Point> find_minima_at_t(std::vector<Eigen::VectorXd> test_points, double T) const;

  void find_phases();

  double delta_potential_at_T(const Phase *phase1, const Phase *phase2, double T) const {
    return phase_at_T(phase1, T).potential - phase_at_T(phase2, T).potential;
  }

  Point phase_at_T(const Phase *phase, double T) const;

  /** Pretty-printer for a collection of phases */
  friend std::ostream& operator << (std::ostream& o, const PhaseFinder& pf) {
    auto phases = pf.phases;

    o << "found " << phases.size() << " phase";
    if (phases.size() != 1) {
      o << "s";
    }
    o << std::endl << std::endl;

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

  /** Generate symmetric partner for a point */
  std::vector<Eigen::VectorXd> symmetric_partners(const Eigen::VectorXd a) const;
  
  /** Check that two minima are identical to within a particular tolerance */
  bool identical_within_tol(const Eigen::VectorXd a, const Eigen::VectorXd b) const;

 private:
  EffectivePotential::Potential &P;

  /**
     Find local minima at a particular temperature. The overloads define
     different ways of passing an initial step.
  */
  Point find_min(const Eigen::VectorXd X, double T, double step) const;
  Point find_min(const Eigen::VectorXd X, double T, Eigen::VectorXd step) const;
  Point find_min(const Eigen::VectorXd X, double T) const;

  /** Check for a jump discontinuity between two phases*/
  bool jump(const Eigen::VectorXd a, const Eigen::VectorXd b) const;
  /** Get descriptor of a particular minimum */
  minima_descriptor get_minima_descriptor(const Point minima) const;
  /**
    @overload
    Get descriptor of the deepest minimum of a set at a particular temperature
  */
  minima_descriptor get_minima_descriptor(const std::vector<Point> minima, double T) const;
  /** Get deepest minimum of a set at a particular temperature */
  Point get_deepest_minima(const std::vector<Point> minima, double T) const;
  /** Check whether the origin is the unique minima */
  bool origin_unique_minima(const std::vector<Point> minima) const;

  static double wrap_nlopt(const std::vector<double> &x,
                           std::vector<double> &grad, void *data) {
    std::function<double(Eigen::VectorXd)> potential = *static_cast<std::function<double(Eigen::VectorXd)>*>(data);
    const auto phi = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(const_cast<double*>(x.data()), x.size());
    return potential(phi);
  }

  phase_end_descriptor trace_minimum(Point start, double tstop,
                                    double dtstart_, std::vector<Eigen::VectorXd> *X,
                                    std::vector<double> *T, std::vector<Eigen::VectorXd> *dXdT, std::vector<double> *V,
                                    Point *jumped) const;

  /** Expected change in minimum with temperature, dx/dt */
  Eigen::VectorXd dx_min_dt(Eigen::VectorXd X, double T) const;
  /** @overload Precomputed Hessian matrix */
  Eigen::VectorXd dx_min_dt(Eigen::MatrixXd hessian, Eigen::VectorXd X, double T) const;

  /** Check whether two phases are redundant */
  bool redundant(const Phase *phase1, const Phase *phase2, end_descriptor end = BOTH) const;

  /** Check whether point belongs to a known phase */
  bool belongs_known_phase(Point point) const;

  /** Electroweak vacuum at zero temperature */
  bool consistent_vacuum(Eigen::VectorXd x) const;

  /** Check whether point is out of boundary */
  bool out_of_bounds(Eigen::VectorXd x) const;

  /** Remove/combine identical phases */
  void remove_redundant();

  /** Check whether Hessian is singular */
  bool hessian_singular(Eigen::VectorXd X, double T) const;
  /** @overload Precomputed Hessian matrix */
  bool hessian_singular(Eigen::MatrixXd hessian, Eigen::VectorXd X, double T) const;

  /** Check whether Hessian is positive definite  */
  bool hessian_positive_definite(Eigen::VectorXd X, double T) const;
  /** @overload Precomputed Hessian matrix */
  bool hessian_positive_definite(Eigen::MatrixXd hessian, Eigen::VectorXd X, double T) const;

  /** Default bound on fields */
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
  PROPERTY(double, t_low, 0.)
  /** Highest temperature to consider */
  PROPERTY(double, t_high, 1000.)
  /** The starting step-size relative to t_high - t_low */
  PROPERTY(double, dt_start_rel, 0.01)
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
  PROPERTY(std::vector<Phase>, phases, {})
  /** Number of scalar fields */
  PROPERTY_CUSTOM_SETTER(size_t, n_scalars, 0)
  /** Number of scalar fields that could break electroweak symmetry */
  PROPERTY(size_t, n_ew_scalars, 0)
  /** Zero-temperature vacuum expectation value of electroweak charged scalars */
  PROPERTY(double, v, 246.)
  /** Check whether deepest vacuum at low temperature agrees with that expected */
  PROPERTY(bool, check_vacuum_at_low, true)
  /** Check whether there is a unique vacuum at high temperature */
  PROPERTY(bool, check_vacuum_at_high, true)
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
  PROPERTY(double, phase_min_length, 0.5)
  /** Guesses for locations of minima */
  PROPERTY(std::vector<Eigen::VectorXd>, guess_points, {})
};

}  // namespace PhaseTracer

#endif
