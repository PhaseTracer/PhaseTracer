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

#ifndef PHASETRACER_TRANSITION_FINDER_HPP_
#define PHASETRACER_TRANSITION_FINDER_HPP_

#include <algorithm>
#include <ostream>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <boost/cstdint.hpp>

#include "phase_finder.hpp"
#include "overload.hpp"
#include "potential.hpp"
#include "action_calculator.hpp"
#include "transition_graph_util.hpp"

namespace TransitionGraph {
struct Path;
}

namespace PhaseTracer {

/** Information about root-finding for a transition */
enum Message { SUCCESS,
               NON_OVERLAPPING_T,
               ERROR };

// A polynomial fitting function created by DeepSeek.
class PolynomialFitterEigen {
public:
  void fit(const std::vector<double> &T_list_, const std::vector<double> &S_list_, double TC_, int degree) {
    TC = TC_;
    T_list = T_list_;
    S_list = S_list_;
    fit_flag = true;

    // Filter the nodes
    std::vector<double> T_select;
    std::vector<double> S_select;
    T_select.reserve(T_list.size());
    S_select.reserve(T_list.size());
    for (size_t i = 0; i < T_list.size(); ++i) {
      double S3T = S_list[i] / T_list[i];
      if (S3T > 1 and S3T < 1000) {
        T_select.push_back(T_list[i]);
        S_select.push_back(S_list[i]);
      }
    }
    T_select.shrink_to_fit();
    S_select.shrink_to_fit();

    // TODO: this is only for test
    MSE_pre = fit_(T_select, S_select, degree);
    LOG(debug) << "MSE of action fit (pre selection) = " << MSE_pre;

    select(T_select, S_select);
    if (T_select.size() < degree * 2)
      return;

    MSE = fit_(T_select, S_select, degree);
    LOG(debug) << "MSE of action fit = " << MSE;
//    if (MSE < 10)
      success = true;
  }

  double fit_(const std::vector<double> T_select, const std::vector<double> S_select, int degree) {

    int n = T_select.size();
    std::vector<double> x = T_select;
    std::vector<double> y;
    y.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      y.push_back(S_select[i] * pow(TC - T_select[i], 2));
    }
    y.shrink_to_fit();

    Eigen::MatrixXd X(n, degree + 1);
    Eigen::VectorXd Y(n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j <= degree; j++) {
        X(i, j) = std::pow(x[i], j);
      }
      Y(i) = y[i];
    }
    coefficients = (X.transpose() * X).ldlt().solve(X.transpose() * Y);

    fit_range_min = *std::min_element(T_select.begin(), T_select.end());
    fit_range_max = *std::max_element(T_select.begin(), T_select.end());

    // Calculate mean square error
    double MSE_ = 0;
    S_predict.clear();
    S_predict.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      S_predict.push_back(predict(T_select[i]));
      MSE_ += pow(S_select[i] / T_select[i] - S_predict[i] / pow(TC - T_select[i], 2) / T_select[i], 2);
    }
    S_predict.shrink_to_fit();
    MSE_ = MSE_ / n;
    return MSE_;
  }

  void select(std::vector<double> &x, std::vector<double> &y) {

    if (x.size() != y.size() || x.size() < 3) {
      throw std::runtime_error("Need at least 3 points for action fit");
    }
    if (x[0] < x.back()) {
      std::reverse(x.begin(), x.end());
      std::reverse(y.begin(), y.end());
    }
    std::vector<int> selected_indices;
    for (size_t i = 0; i < y.size() - 1; ++i) {
      if (y[i + 1] - y[i] < 0) {
        selected_indices.push_back(i);
        selected_indices.push_back(i + 1);
        break;
      }
    }
    if (selected_indices.empty()) {
      return;
    }

    double y_prev = y[selected_indices.back()];
    double diff_prev = y[selected_indices[1]] - y[selected_indices[0]];
    for (size_t i = selected_indices.back() + 1; i < y.size(); ++i) {
      double diff = y[i] - y_prev;
      if ((diff < 0) && (std::abs(diff) < std::abs(diff_prev * 2))) {
        selected_indices.push_back(i);
        diff_prev = diff;
        y_prev = y[i];
      }
    }
    std::vector<double> x_selected;
    std::vector<double> y_selected;
    for (int idx : selected_indices) {
      x_selected.push_back(x[idx]);
      y_selected.push_back(y[idx]);
    }
    x = std::move(x_selected);
    y = std::move(y_selected);
  }

  double predict(double x) const {
    if (x < fit_range_min or x > fit_range_max) {
      LOG(warning) << "Input value exceeds the action fitting function's domain.";
    }
    double result = 0.0;
    for (int i = 0; i < coefficients.size(); i++) {
      result += coefficients[i] * std::pow(x, i);
    }
    return result;
  }
  std::vector<double> predict(const std::vector<double> &x) const {
    std::vector<double> result;
    result.reserve(x.size());
    for (int i = 0; i < x.size(); i++) {
      result.push_back(predict(x[i]));
    }
    return result;
  }

  double derivative(double x) const {
    double result = 0.0;
    for (int i = 1; i < coefficients.size(); i++) {
      result += i * coefficients[i] * std::pow(x, i - 1);
    }
    return result;
  }
  double secondDerivative(double x) const {
    double result = 0.0;
    for (int i = 2; i < coefficients.size(); i++) {
      result += i * (i - 1) * coefficients[i] * std::pow(x, i - 2);
    }
    return result;
  }

  double findLocalMinimum(double initial_guess = 0.0, double tolerance = 1e-4, int max_iterations = 100) const {
    double x = initial_guess;
    for (int i = 0; i < max_iterations; i++) {
      double f = derivative(x);
      double f_prime = secondDerivative(x);
      if (std::abs(f_prime) < 1e-12 || f_prime < 0) {
        x = x - 0.1 * f;
        continue;
      }
      double delta = f / f_prime;
      x = x - delta;
      if (std::abs(delta) < tolerance && f_prime > 0) {
        break;
      }
    }
    if (secondDerivative(x) > 0) {
      return x;
    } else {
      return std::numeric_limits<double>::quiet_NaN();
    }
  }

  const Eigen::VectorXd &getCoefficients() const {
    return coefficients;
  }

  const double get_MSE() const {
    return MSE;
  }

  const double get_MSE_pre() const {
    return MSE_pre;
  }

  const double get_S3T(double x) const {
    double S3T = predict(x) / pow(TC - x, 2) / x;
    return std::min(S3T, 1E30);
  }

  const double get_fit_flag() const {
    return fit_flag;
  }

  const bool get_success() const {
    return success;
  }

  const std::vector<double> get_T_list() const {
    return T_list;
  }

  const std::vector<double> get_S_list() const {
    return S_list;
  }

private:
  Eigen::VectorXd coefficients;
  double TC;
  std::vector<double> T_list;
  std::vector<double> S_list;
  std::vector<double> S_predict;
  double fit_range_min;
  double fit_range_max;
  double MSE = 1E10;
  double MSE_pre = 1E10; // TODO: this is only for test
  bool success = false;
  bool fit_flag = false;
};

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
  size_t id;

  double TN;
  Eigen::VectorXd true_vacuum_TN;
  Eigen::VectorXd false_vacuum_TN;
  bool subcritical = false;

  PolynomialFitterEigen action_curve;
  double TP;
  
  Transition(Message message) : message(message) {};

  Transition(Message message,
             double TC,
             Phase true_phase,
             Phase false_phase,
             Eigen::VectorXd true_vacuum,
             Eigen::VectorXd false_vacuum,
             double gamma,
             std::vector<bool> changed,
             double delta_potential,
             size_t key,
             size_t id) : message(message),
                          TC(TC),
                          true_phase(true_phase),
                          false_phase(false_phase),
                          true_vacuum(true_vacuum),
                          false_vacuum(false_vacuum),
                          gamma(gamma),
                          changed(changed),
                          delta_potential(delta_potential),
                          key(key),
                          id(id) {}

  void set_subcritical(bool subcritical_) {
    subcritical = subcritical_;
  }

  void set_nucleation(double T, Eigen::VectorXd true_vacuum, Eigen::VectorXd false_vacuum) {
    calculate_action = true;
    TN = T;
    true_vacuum_TN = true_vacuum;
    false_vacuum_TN = false_vacuum;
  }
  void set_percolation(double T ){
    calculate_percolation = true;
    TP=T;
  }

  void set_action_curve(PolynomialFitterEigen action_curve_) {
    action_curve = action_curve_;
  }
  /** Pretty-printer for single transition */
  friend std::ostream &operator<<(std::ostream &o, const Transition &a) {
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
        << "delta potential (TC) = " << a.delta_potential << std::endl;

      if (a.calculate_action) {
        if (a.action_curve.get_fit_flag()) {
          if (a.action_curve.get_success()) {
            o << "Action curve fitting succeeded with MSE = "
              << a.action_curve.get_MSE() << std::endl;
          } else {
            o << "Action curve fitting failed, MSE = "
              << a.action_curve.get_MSE() << std::endl;
          }
        }
        o << "TN = " << a.TN << std::endl
          << "false vacuum (TN) = " << a.true_vacuum_TN << std::endl
          << "true vacuum (TN) = " << a.false_vacuum_TN << std::endl;
        if (a.calculate_percolation)
          o << "TP = " << a.TP << std::endl;
      } else {
        o << "did not calculate action or check nucleation" << std::endl;
      }

      o << "transition was "
        << (a.subcritical ? "subcritical" : "not subcritical") << std::endl;

    } else {
      o << "=== failure. message =  " << a.message << " ===" << std::endl;
    }
    return o;
  }

private:
  bool calculate_action = false;
  bool calculate_percolation = false;
};

class TransitionFinder {
public:
  explicit TransitionFinder(PhaseFinder &pf_) : pf(pf_), ac(pf_.P) {
    calculate_action = false;
  }
  explicit TransitionFinder(PhaseFinder &pf_, ActionCalculator ac_) : pf(pf_), ac(ac_) {
    calculate_action = true;
  }
  virtual ~TransitionFinder() = default;

  /** Find all transitions between all phases */
  void find_transitions();

  /** Called from find_transitions; adds subcritical transitions to the transitions vector */
  void append_subcritical_transitions();

  /** Called from append_subcritical_transitions; checks whether there is a subcritical transition between two phases. */
  bool checkSubcriticalTransition(const std::vector<PhaseTracer::Phase> &phases, int i, int j, double Tmax,
                                  double energyAtTmax, bool checkFromNewPhase, std::vector<bool> &isTransitionedTo);

  /** Find all transition paths  */
  void find_transition_paths(const EffectivePotential::Potential &model, bool knownHighTPhase);

  /** Retrieve all transitions between all phases */
  std::vector<Transition> get_transitions() const { return transitions; }

  std::vector<Eigen::VectorXd> get_vacua_at_T(const Phase &phase1, const Phase &phase2, double T, size_t i_unique = 0) const;

  double get_action(const Eigen::VectorXd &vacuum_1, const Eigen::VectorXd &vacuum_2, double T) const;

  double get_action(const Phase &phase1, const Phase &phase2, double T, size_t i_unique = 0) const;

  std::vector<double> get_action(const Phase &phase1, const Phase &phase2, std::vector<double> T_list, size_t i_unique = 0) const;

  void write_action_to_text(const Phase &phase1, const Phase &phase2, std::vector<double> T_list, const std::string &filename, size_t i_unique = 0) const;

  void write_action_to_text(const Transition &tran, double T_min, double T_max, size_t n_step, const std::string &filename, size_t i_unique = 0) const;

  void write_action_to_text(const Transition &tran, const std::string &filename, size_t n_step = 50, size_t i_unique = 0) const;

  PolynomialFitterEigen get_action_curve(const Phase &phase1, const Phase &phase2, size_t i_unique, double T_begin, double T_end) const;

  std::pair<double, PolynomialFitterEigen> get_Tnuc(const Phase &phase1, const Phase &phase2, size_t i_unique, double T_begin, double T_end) const;

  
  double nucleation_rate(const Phase &phase1, const Phase &phase2, size_t i_unique, double T, PolynomialFitterEigen action_curve) const {
    double S3_over_T;
    if (action_curve.get_success()){
      S3_over_T = action_curve.get_S3T(T);
    }
    else{
      S3_over_T = get_action(phase1, phase2, T, i_unique)/T;
    }
    
    return pow(T, 4) * std::exp(-S3_over_T);
  }
  
  double false_vacuum_fraction_integrand(const Phase &phase1, const Phase &phase2, size_t i_unique, double T, double T_min, double vw, PolynomialFitterEigen action_curve) const {
    double nucl_rate = nucleation_rate(phase1, phase2, i_unique, T, action_curve);
    return -4./3. * M_PI * std::pow(vw, 3) * nucl_rate * pow(1./(T_min * T_min)-1./(T * T), 3) * 1/pow(T, 3);
  }
  
  double get_false_vacuum_fraction(const Phase &phase1, const Phase &phase2, size_t i_unique, double init_T, double end_T,  int num_T_list, PolynomialFitterEigen action_curve) const {
      const double G = 6.7088e-39;
      const double C = std::sqrt(8 * pow(M_PI, 3) * G * dof / 90);

      Eigen::VectorXd vec(num_T_list);
      vec.setLinSpaced(num_T_list, end_T, init_T);

      alglib::real_1d_array T_list;
      alglib::real_1d_array integrand;

      T_list.setlength(num_T_list);
      integrand.setlength(num_T_list);

      for (size_t i = 0; i < num_T_list; ++i) {
          T_list[i] = vec[i];
          integrand[i] = false_vacuum_fraction_integrand(phase1, phase2, i_unique, T_list[i], end_T, vw, action_curve);
      }

      alglib::spline1dinterpolant spline;
      alglib::spline1dbuildcubic(T_list, integrand, spline);

      double S = alglib::spline1dintegrate(spline, T_list[num_T_list - 1]);
      return std::exp(S / (8 * pow(C, 4)));
      //return S;
  }
  
  double get_percolation_temperature(const Phase &phase1, const Phase &phase2, size_t i_unique, double init_T, double end_T, PolynomialFitterEigen action_curve) const {
      double target = 0.7;
      double false_vacuum_init = get_false_vacuum_fraction(phase1, phase2, i_unique, init_T, init_T - Tperc_tol_rel, num_T_list, action_curve) - target;
      double false_vacuum_end = get_false_vacuum_fraction(phase1, phase2, i_unique, init_T, end_T, num_T_list, action_curve) - target;
      double Tp;
      double init_T_fix = init_T;

      if (false_vacuum_init * false_vacuum_end > 0) {
        LOG(error) << "false_vacuum_fraction at T_init = " << false_vacuum_init;
        LOG(error) << "false_vacuum_end at T_end = " << false_vacuum_init;
        LOG(error) << "Error: f(a) and f(b) must have opposite signs!";
        return 0;
      }

      while ((init_T - end_T) > Tperc_tol_rel) {
          double mid_T = (init_T + end_T) / 2.0;
          double fmid = get_false_vacuum_fraction(phase1, phase2, i_unique, init_T_fix, mid_T,  num_T_list, action_curve) - target;

          if (fabs(fmid) < Tperc_tol_rel) {
              return mid_T;
          }

          if (false_vacuum_end * fmid < 0) {
              init_T = mid_T;
              false_vacuum_init = fmid;
          } else {
              end_T = mid_T;
              false_vacuum_end = fmid;
          }
      }

      Tp = (init_T + end_T) / 2.0;
      return Tp;
  }
  
  
  /** Retrieve all transition paths */
  std::vector<TransitionGraph::Path> get_transition_paths() const { return transition_paths; }

  /** Retrieve all phases */
  std::vector<Phase> get_phases() const { return pf.get_phases(); }

  /** Pretty-printer for set of transitions in this object */
  friend std::ostream &operator<<(std::ostream &o, const TransitionFinder &a);

  /** Object with phases and potential */
  PhaseFinder &pf;

private:
  ActionCalculator ac;

  /** Whether already calculated all transitions */
  bool calculated_transitions = false;

  /** Container for all transitions between any two phases */
  std::vector<Transition> transitions;

  /** Container for all transition paths. */
  std::vector<TransitionGraph::Path> transition_paths;

  /** Find transitions between two phases between two temperatures */
  std::vector<Transition> find_transition(const Phase &phase1, const Phase &phase2, double T1, double T2, size_t currentID) const;

  std::vector<Transition> symmetric_partners(const Phase &phase1, const Phase &phase2, double TC, size_t currentID) const;

  /** Find critical temperature between two phases */
  double find_critical_temperature(const Phase &phase1, const Phase &phase2, double T1, double T2) const;

  /** Find many transitions between two phases at a particular resolution */
  std::vector<Transition> divide_and_find_transition(const Phase &phase1, const Phase &phase2, double T1, double T2, size_t currentID) const;

  /** Find un-overlapped temperature region between the two phases*/
  std::vector<double> get_un_overlapped_T_range(const Phase &phase1, const Phase &phase2, double T1, double T2) const;

  /** Strength of phase transition for first N fields */
  double gamma(const Eigen::VectorXd &true_vacuum, const Eigen::VectorXd &false_vacuum, const double TC) const;

  /** Note which VEVs changed */
  std::vector<bool> changed(const Eigen::VectorXd &true_vacuum, const Eigen::VectorXd &false_vacuum) const;

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
  /** Relative precision in nucleation temperature */
  PROPERTY(double, Tnuc_tol_rel, 1.e-3)
  PROPERTY(bool, fit_action_curve, true)
  /** Number of nodes for getting action curve */
  PROPERTY(double, action_curve_nodes, 30)
  /** Order of polynomial fitting for action curve */
  PROPERTY(double, action_curve_order, 6)
  PROPERTY(double, vw, 0.3)
  PROPERTY(double, dof, 106.75)
  PROPERTY(double, num_T_list, 300)
  PROPERTY(double, Tperc_tol_rel, 1e-4)
  
  PROPERTY(bool, calculate_action, false)

  PROPERTY(bool, calculate_percolation, false)
  
  /**
   * Whether we should check for subcritical transitions, i.e. transitions between phases when one phase is strictly of
   * lower energy than the other. This can occur when a new phase appears near a phase that is not of the highest energy
   * at that temperature.
   */
  PROPERTY(bool, check_subcritical_transitions, false)
};

} // namespace PhaseTracer

#endif // PHASETRACER_TRANSITION_FINDER_HPP_
