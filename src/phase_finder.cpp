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

#include <Eigen/Eigenvalues>
#include <boost/math/special_functions/sign.hpp>
#include <limits.h>
#include <iostream>

#include "phase_finder.hpp"
#include "logger.hpp"
#include "pow.hpp"

namespace PhaseTracer {

std::ostream& operator << (std::ostream& o, const phase_end_descriptor& d) {
  switch (d) {
    case REACHED_T_STOP:
      o << "Reached tstop";
      break;
    case FORBIDDEN_OR_BOUNDS:
      o << "Reached forbidden or out of bounds region";
      break;
    case HESSIAN_SINGULAR:
      o << "Hessian was singular";
      break;
    case HESSIAN_NOT_POSITIVE_DEFINITE:
      o << "Hessian was not positive definite";
      break;
    case JUMP_INDICATED_END:
      o << "Jump in fields indicated end of phase";
      break;
    default:
      throw std::runtime_error("Unknown type");
  }
  return o;
}

PhaseFinder::PhaseFinder(EffectivePotential::Potential &potential) : P(potential) {
  n_scalars = P.get_n_scalars();
  upper_bounds = std::vector<double>(n_scalars, bound);
  lower_bounds = std::vector<double>(n_scalars, -bound);
  set_seed(seed);
}

const EffectivePotential::Potential& PhaseFinder::get_potential() const
{
  return P;
}

std::vector<Eigen::VectorXd> PhaseFinder::generate_test_points() const {
  std::vector<Eigen::VectorXd> test_points(guess_points);
  const auto a = Eigen::Map<const Eigen::ArrayXd, Eigen::Unaligned>(lower_bounds.data(), lower_bounds.size());
  const auto b = Eigen::Map<const Eigen::ArrayXd, Eigen::Unaligned>(upper_bounds.data(), upper_bounds.size());
  while (test_points.size() < n_test_points) {
    const Eigen::VectorXd point = 0.5 * (b - a) * Eigen::ArrayXd::Random(n_scalars) + 0.5 * (b + a);
    test_points.push_back(point);
  }
  return test_points;
}

std::vector<Point> PhaseFinder::find_minima_at_t(double T) const {
  std::vector<Point> minima;
  const std::vector<Eigen::VectorXd> test_points = generate_test_points();

  for (const auto& p : test_points) {
    const auto polished = find_min(p, T, find_min_locate_abs_step);
    bool duplicate = false;

    for (const auto& m : minima) {
      if (identical_within_tol(polished.x, m.x)) {
        duplicate = true;
        break;
      }
    }

    if (!duplicate) {
      minima.push_back(polished);
    }
  }
  return minima;
}

bool PhaseFinder::consistent_vacuum(const Eigen::VectorXd& x) const {
  if (n_ew_scalars == 0) {
    return true;
  }
  const double found = x.segment(0, n_ew_scalars).norm();
  return std::abs(found - v) < x_abs_identical;
}

std::vector<Eigen::VectorXd> PhaseFinder::symmetric_partners(const Eigen::VectorXd& a) const {
  std::vector<Eigen::VectorXd> partners;
  partners.push_back(a);
  for (size_t i=0; i < P.apply_symmetry(a).size(); i++) {
    const size_t n = partners.size();
    for (size_t j=0; j< n; j++) {
      const auto x = partners[j];
      const auto x_ = P.apply_symmetry(x)[i];
      partners.push_back(x_);
    }
  }
  return partners;
}

bool PhaseFinder::identical_within_tol(const Eigen::VectorXd& a, const Eigen::VectorXd& b) const {
  double min_distance = (a-b).norm();
  for (const auto b_ : symmetric_partners(b)) {
    min_distance = std::min(min_distance, (a-b_).norm());
  }
  return min_distance < x_abs_identical + x_rel_identical * std::max(a.norm(), b.norm());
}

bool PhaseFinder::jump(const Eigen::VectorXd& a, const Eigen::VectorXd& b) const {
  return (a - b).norm() > x_abs_jump + x_rel_jump * std::max(a.norm(), b.norm());
}

Point PhaseFinder::get_deepest_minima(const std::vector<Point>& minima, double T) const {
  std::vector<double> potential;
  for (const auto &m : minima) {
    potential.push_back(m.potential);
  }

  const int deepest_index = std::min_element(potential.begin(), potential.end()) - potential.begin();
  const Point deepest = minima[deepest_index];

  LOG(debug) << "At T = "
             << T
             << ", the deepest minimum = "
             << deepest.x;

  return deepest;
}

minima_descriptor PhaseFinder::get_minima_descriptor(const Point& minima) const {
  if (identical_within_tol(minima.x, Eigen::VectorXd::Zero(n_scalars))) {
    return ORIGIN;
  }
  if (consistent_vacuum(minima.x)) {
    return CONSISTENT_VACUUM;
  }
  return OTHER;
}

minima_descriptor PhaseFinder::get_minima_descriptor(const std::vector<Point>& minima, double T) const {
  const Point deepest = get_deepest_minima(minima, T);
  return get_minima_descriptor(deepest);
}

bool PhaseFinder::origin_unique_minima(const std::vector<Point>& minima) const {
  for (const auto &m : minima) {
    if (!identical_within_tol(m.x, Eigen::VectorXd::Zero(n_scalars))) {
      return false;
    }
  }
  return true;
}

bool PhaseFinder::belongs_known_phase(const Point& point) const {
  for (const auto& phase : phases) {
    if (!phase.contains_t(point.t)) {
      continue;
    }

    const Eigen::VectorXd x = phase_at_T(phase, point.t).x;
    if (identical_within_tol(x, point.x)) {
      return true;
    }
  }
  return false;
}

std::vector<Point> PhaseFinder::get_minima_at_t_low() {
  if (minima_at_t_low.empty()) {
    LOG(debug) << "Check potential at T = t_low = " << t_low;
    minima_at_t_low = find_minima_at_t(t_low);
  }
  return minima_at_t_low;
}

std::vector<Point> PhaseFinder::get_minima_at_t_high() {
  if (minima_at_t_high.empty()) {
    LOG(debug) << "Check potential at T = t_high = " << t_high;
    minima_at_t_high = find_minima_at_t(t_high);
  }
  return minima_at_t_high;
}

std::vector<Phase> PhaseFinder::get_phases_at_T(double T) {
  if (phases.size() == 0) find_phases();
  std::vector<Phase> phases_at_T;
  for (const auto& pi: phases){
    if (T >= pi.T.front() and T <= pi.T.back()) phases_at_T.push_back(pi);
  }
  return phases_at_T;
}

Phase PhaseFinder::get_deepest_phase_at_T(double T) {
  const auto phases_at_T = get_phases_at_T(T);
  if (phases_at_T.size()==0)
    throw std::runtime_error("There is no phase at T = "+std::to_string(T));

  size_t i_deepest = 0 ;
  for (size_t i = 0; i<phases_at_T.size(); i++){
    if (phase_at_T(phases_at_T[i],T).potential < phase_at_T(phases_at_T[i_deepest],T).potential)
      i_deepest = i;
  }
  return phases_at_T.at(i_deepest);
}

void PhaseFinder::find_phases() {
  LOG(debug) << "Find global minima from " << guess_points.size() << " guesses";

  for (const auto& g : guess_points) {
    if (g.size() != n_scalars) {
      throw std::invalid_argument("Guesses have wrong size compared to number of fields");
    }
  }

  LOG(debug) << "Check potential at T = t_low = " << t_low;
  minima_at_t_low = find_minima_at_t(t_low);

  if (check_vacuum_at_low) {
      const minima_descriptor location = get_minima_descriptor(minima_at_t_low, t_low);
    if (location == ORIGIN) {
      throw std::runtime_error("No minimum lower than the origin at t_low");
    } else if (location != CONSISTENT_VACUUM) {
      throw std::runtime_error("The deepest minimum at t_low is not the electroweak VEV");
    }
  }
  LOG(debug) << "Check potential at T = t_high = " << t_high;
  minima_at_t_high = find_minima_at_t(t_high);
  if (check_vacuum_at_high && !origin_unique_minima(minima_at_t_high)) {
    throw std::runtime_error("Found minimum other than the origin at T = t_high");
  }

  // Concatenate minima at high and low temperature
  std::vector<Point> points(minima_at_t_low);
  points.insert(points.end(), minima_at_t_high.begin(), minima_at_t_high.end());

  LOG(debug) << "Start phase tracing with " << points.size() << " guesses";
  for (const auto &p : points) {
    LOG(debug) << p;
  }

  while (!points.empty()) {
    // Fetch point
    auto point = points.back();
    points.pop_back();

    LOG(debug) << "Tracing starts at " << point;

    if (point.t < t_low || point.t > t_high) {
      LOG(debug) << "T = " << point.t << " outside t_low and t_high";
      continue;
    }

    if (P.forbidden(point.x)) {
      LOG(debug) << point << " is forbidden";
      continue;
    }

    if (belongs_known_phase(point)) {
      LOG(debug) << point << " belongs to a known phase";
      continue;
    }

    if (out_of_bounds(point.x)) {
      LOG(debug) << point << " out of bounds";
      continue;
    }

    // Trace

    LOG(debug) << "Tracing down from " << point;
    const double dt_start = (t_high - t_low) * dt_start_rel;
    std::vector<Eigen::VectorXd> X_down;
    std::vector<double> T_down;
    std::vector<Eigen::VectorXd> dXdT_down;
    std::vector<double> V_down;
    Point jumped_down;
    const auto end_low = trace_minimum(point, t_low, -dt_start, &X_down, &T_down, &dXdT_down, &V_down, &jumped_down);

    LOG(debug) << "Tracing up from " << point;
    std::vector<Eigen::VectorXd> X_up;
    std::vector<double> T_up;
    std::vector<Eigen::VectorXd> dXdT_up;
    std::vector<double> V_up;
    Point jumped_up;
    const auto end_high = trace_minimum(point, t_high, dt_start, &X_up, &T_up, &dXdT_up, &V_up, &jumped_up);

    LOG(debug) << "Combine results from tracing up and down";
    std::vector<Eigen::VectorXd> X = X_down;
    std::vector<double> T = T_down;
    std::vector<Eigen::VectorXd> dXdT = dXdT_down;
    std::vector<double> V = V_down;

    // Reverse results from tracing down - as dt was negative
    std::reverse(X.begin(), X.end());
    std::reverse(T.begin(), T.end());
    std::reverse(dXdT.begin(), dXdT.end());
    std::reverse(V.begin(), V.end());

    // Append result of tracing up (don't double count point at which tracing began)
    X.insert(X.end(), X_up.begin() + 1, X_up.end());
    T.insert(T.end(), T_up.begin() + 1, T_up.end());
    dXdT.insert(dXdT.end(), dXdT_up.begin() + 1, dXdT_up.end());
    V.insert(V.end(), V_up.begin() + 1, V_up.end());

//  Ignore short phase may cause endless loop. 
//    if (std::abs(T.front() - T.back()) > phase_min_length) {
      Phase new_;
      new_.key = phases.size();
      new_.X = X;
      new_.T = T;
      new_.dXdT = dXdT;
      new_.V = V;
      new_.end_low = end_low;
      new_.end_high = end_high;
      phases.push_back(new_);
      LOG(debug) << "Added new phase:" << std::endl << new_;
//    } else {
//      LOG(warning) << "Did not add short phase";
//    }

    if (end_high == JUMP_INDICATED_END || end_high == HESSIAN_SINGULAR || end_high == HESSIAN_NOT_POSITIVE_DEFINITE) {
      Point top;
      if (end_high == JUMP_INDICATED_END) {
        top = jumped_up;
      } else {
        top = find_min(X.back(), T.back() + t_jump_rel * (t_high - t_low));
      }
      if (hessian_positive_definite(top.x, top.t)) {
        LOG(debug) << "Where does end nearest t_high go? Appending " << top;
        points.push_back(top);
      }
    }

    if (end_low == JUMP_INDICATED_END || end_low == HESSIAN_SINGULAR || end_low == HESSIAN_NOT_POSITIVE_DEFINITE) {
      Point bottom;
      if (end_low == JUMP_INDICATED_END) {
        bottom = jumped_down;
      } else {
        bottom = find_min(X.front(), T.front() - t_jump_rel * (t_high - t_low));
      }
      if (hessian_positive_definite(bottom.x, bottom.t)) {
        LOG(debug) << "Where does end nearest t_low go? Appending " << bottom;
        points.push_back(bottom);
      }
    }
  }
  // Remove any redundant phases - we only checked that the guess didn't belong
  // to a known phase. The phases might overlap
  remove_redundant();

  if (check_merge_phase_gaps) {
    merge_phase_gaps();
  }
  
  LOG(debug) << "Finished finding phases";
}

phase_end_descriptor PhaseFinder::trace_minimum(Point start, double tstop,
                                               double dt_start,
                                               std::vector<Eigen::VectorXd> *X,
                                               std::vector<double> *T,
                                               std::vector<Eigen::VectorXd> *dXdT,
                                               std::vector<double> *V,
                                               Point *jumped) const {

  const double time_scale = std::abs(t_high - t_low);
  const double dt_min = std::max(dt_min_rel * time_scale, dt_min_abs);
  const double dt_max = std::min(time_scale * dt_max_rel, dt_max_abs);
  double dt = dt_start;
  const double sign_dt = boost::math::sign(dt);
  bool significant_dx = false;
  bool hessian_singular_ = false;
  bool hessian_positive_definite_ = true;

  Eigen::VectorXd x0 = start.x;
  double t0 = start.t;
  double v0 = start.potential;

  Eigen::MatrixXd h0 = P.d2V_dx2(x0, t0);
  if (check_hessian_singular) {
    hessian_singular_ = hessian_singular(h0, x0, t0);
  }
  hessian_positive_definite_ = hessian_positive_definite(h0, x0, t0);

  // Append start of phase
  Eigen::VectorXd dxdt0 = dx_min_dt(h0, x0, t0);
  X->push_back(x0);
  T->push_back(t0);
  dXdT->push_back(dxdt0);
  V->push_back(v0);

  // This loop traces temperature in steps of dt
  for (unsigned int iter = 0; iter <= trace_max_iter; iter++) {

    LOG(trace) << "Tracing at T = " << t0 << " and dT = " << dt;

    if (iter == trace_max_iter) {
      throw std::runtime_error("Exceeded maximum number of iterations tracing minimum");
    }

    // Insignificant change in temperature
    const bool insignificant_dt = std::abs(dt) <= dt_min;

    // If the Hessian is singular, reached a transition
    if (check_hessian_singular && insignificant_dt && hessian_singular_) {
      LOG(debug) << "Hessian singular indicated that phase ended at about T = "
                 << T->back();
      return HESSIAN_SINGULAR;
    }

    // If the Hessian is not positive definite, reach an end of a phase
    if (insignificant_dt && !hessian_positive_definite_) {
      LOG(debug) << "Hessian not positive semi-definite indicated that phase ended at about T = "
                 << T->back();
      return HESSIAN_NOT_POSITIVE_DEFINITE;
    }

    // Insignificant change in temperature but significant change in field
    if (insignificant_dt && significant_dx) {
      LOG(debug) << "Jump in fields indicated that phase ended. Jump"
                 << " from " << X->back()
                 << " to " << jumped->x
                 << " from T = " << T->back()
                 << " to " << jumped->t;
      return JUMP_INDICATED_END;
    }

    const bool reached_tstop = dt > 0 ? T->back() >= tstop : T->back() <= tstop;

    if (reached_tstop) {
      LOG(trace) << "Phase did not end but traced until desired temperature at T = "
                 << tstop;
      return REACHED_T_STOP;
    }

    // Get the field values etc after incrementing the temperature by dt
    // but don't exceed tstop
    const double t1 = dt > 0 ? std::min(t0 + dt, tstop) : std::max(t0 + dt, tstop);
    const Eigen::VectorXd dx0 = dxdt0 * dt;
    const Eigen::VectorXd guess = x0 + dx0;

    // Find minimum at new temperature. If guess isn't acceptable, start from
    // previous point in phase with substantial step
    Point polish;
    if (out_of_bounds(guess) || P.forbidden(guess)) {
      polish = find_min(x0, t1, dx0);
    } else {
      polish = find_min(guess, t1);
    }

    const Eigen::VectorXd x1 = polish.x;
    const double v1 = polish.potential;
    const Eigen::MatrixXd h1 = P.d2V_dx2(x1, t1);
    const Eigen::VectorXd dxdt1 = dx_min_dt(h1, x1, t1);

    // If Hessian singular, phase may be ending
    if (check_hessian_singular) {
      hessian_singular_ = hessian_singular(h1, x1, t1);
      if (hessian_singular_) {
        dt *= 0.5;
        LOG(trace) << "Hessian singular at T = " << t1
                   << ". Reducing step-size to dT = " << dt;
        continue;
      }
    }

    // If Hessian is not positive definite, may be ending
    hessian_positive_definite_ = hessian_positive_definite(h1, x1, t1);
    if (!hessian_positive_definite_) {
      dt *= 0.5;
      LOG(trace) << "Hessian is not positive definite at T = " << t1
                 << ". Reducing step-size to dT = " << dt;
      continue;
    }

    // Check the change in minimum field values when we changed temperature by dt
    // to that expected from derivatives
    significant_dx = (jump(x0 + dx0, x1) || jump(x1 - dxdt1 * dt, x0)) && jump(x0, x1);

    // If the change was significant, reduce step size in temperature
    if (significant_dx) {
      dt *= 0.5;
      // Record minima that was jumped to
      jumped->x = x1;
      jumped->t = t1;
      jumped->potential = v1;
      LOG(trace) << "Possible jump from " << x0
                 << " to " << jumped->x
                 << " from T = " << t0 << " to " << jumped->t
                 << ". Reducing step-size to dT = " << dt;
      continue;
    }

    // If the minimum was out of bounds or forbidden, reduce step size in temperature
    if (out_of_bounds(x1) || P.forbidden(x1)) {
      dt *= 0.5;
      if (std::abs(dt) <= dt_min) {
        LOG(warning) << "Traced phase until reached a forbidden region or bound."
                     << " Fields = " << x1
                     << " and T = " << t1;
        return FORBIDDEN_OR_BOUNDS;
      }
      continue;
    }

    // Tweak the temperature step to get desired precision (but insure that
    // it remains between maximum and minimum allowed)

    if (identical_within_tol(guess , x1)) {
      dt = sign_dt * std::min(std::abs(dt) * 2., dt_max);
    } else {
      dt = sign_dt * std::max(std::abs(dt) * 0.5, dt_min);
    }

    // Add points and continue loop

    X->push_back(x1);
    T->push_back(t1);
    dXdT->push_back(dxdt1);
    V->push_back(v1);

    x0 = x1;
    t0 = t1;
    dxdt0 = dxdt1;
  }

  throw std::runtime_error("This should be unreachable");
}

bool PhaseFinder::hessian_singular(const Eigen::VectorXd& X, double T) const {
  return hessian_singular(P.d2V_dx2(X, T), X, T);
}

bool PhaseFinder::hessian_singular(const Eigen::MatrixXd& hessian, const Eigen::VectorXd& X, double T) const {
  const double t_min = hessian.eigenvalues().cwiseAbs().minCoeff();
  const double zero_t_min = P.d2V_dx2(X, 0.).eigenvalues().cwiseAbs().minCoeff();
  return std::abs(t_min) < hessian_singular_rel_tol * std::abs(zero_t_min);
}

bool PhaseFinder::hessian_positive_definite(const Eigen::VectorXd& X, double T) const {
  return hessian_positive_definite(P.d2V_dx2(X, T), X, T);
}

bool PhaseFinder::hessian_positive_definite(const Eigen::MatrixXd& hessian, const Eigen::VectorXd& X, double T) const {
  auto eivals = hessian.eigenvalues();
  for (int i = 0; i < eivals.size(); i++) {
    if (eivals[i].imag() != 0. || eivals[i].real() < 0.) {
      return false;
    }
  }
  return true;
}

Eigen::VectorXd PhaseFinder::dx_min_dt(const Eigen::VectorXd& X, double T) const {
  return dx_min_dt(P.d2V_dx2(X, T), X, T);
}

Eigen::VectorXd PhaseFinder::dx_min_dt(const Eigen::MatrixXd& hessian, const Eigen::VectorXd& X, double T) const {
  const Eigen::VectorXd b = -P.d2V_dxdt(X, T);
  const Eigen::VectorXd dxdt = hessian.colPivHouseholderQr().solve(b);
  const bool check = b.isApprox(hessian * dxdt, linear_algebra_rel_tol);

  if( check_dx_min_dt ) {
    if (!check) {
      throw std::runtime_error("Failed to find dxdt");
    }
  }

  return dxdt;
}

Point PhaseFinder::find_min(const Eigen::VectorXd& guess, double T) const {
  return find_min(guess, T, Eigen::VectorXd::Constant(n_scalars, find_min_trace_abs_step));
}

Point PhaseFinder::find_min(const Eigen::VectorXd& guess, double T, double abs_step) const {
  return find_min(guess, T, Eigen::VectorXd::Constant(n_scalars, abs_step));
}

std::function<double(Eigen::VectorXd)> PhaseFinder::make_objective(double T) const {
  std::function<double(Eigen::VectorXd)> objective = [this, T](Eigen::VectorXd x) {
    return this->P.V(x, T);
  };
  return objective;
}

Point PhaseFinder::find_min(const Eigen::VectorXd& guess, double T, Eigen::VectorXd step) const {
  if (out_of_bounds(guess)) {
    LOG(fatal) << "guess = " << guess << " was out of bounds";
    throw std::runtime_error("guess for nlopt was out of bounds");
  }

  for (int i = 0; i < step.size(); i++) {
    if (std::abs(step[i]) <= find_min_min_step) {
      step[i] = step[i] >= 0 ? find_min_min_step : -find_min_min_step;
      LOG(warning) << "Increasing step size. step = " << step;
    }
  }

  nlopt::opt opt(find_min_algorithm, n_scalars);
  opt.set_xtol_rel(find_min_x_tol_rel);
  opt.set_xtol_abs(find_min_x_tol_abs);
  opt.set_maxtime(find_min_max_time);
  opt.set_maxeval(find_min_max_f_eval);
  std::vector<double> step_vector(step.data(), step.data() + step.rows() * step.cols());

  opt.set_initial_step(step_vector);

  auto objective = make_objective(T);
  opt.set_min_objective(wrap_nlopt, &objective);
  std::vector<double> minima(guess.data(), guess.data() + guess.rows() * guess.cols());
  double potential_at_minima;
  nlopt::result result;

  try {
    result = opt.optimize(minima, potential_at_minima);
  } catch (...) {
    LOG(fatal) << "nlopt exception at T = " << T
               << " and guess = "
               << guess
               << ". step = " << step;
    throw;
  }

  const Eigen::VectorXd minima_VectorXd = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(minima.data(), minima.size());

  if (result < 0) {
    LOG(fatal) << "nlopt failed at T = " << T
               << " and guess = "
               << guess
               << ". error code = " << result
               << ". result = " << minima_VectorXd
               << ". step = " << step;
    std::runtime_error("nlopt failed");
  } else if (result == 5 || result == 6) {
    LOG(warning) << "nlopt timed out/reached maximum function evaluations at T = " << T
                 << " and guess = "
                 << guess
                 << ". error code = " << result
                 << ". result = " << minima_VectorXd
                 << ". step = " << step;
  } else {
    LOG(trace) << "nlopt successful at T = " << T
               << " and guess = "
               << guess
               << ". error code = " << result
               << ". result = " << minima_VectorXd;
  }

  return {minima_VectorXd, potential_at_minima, T};
}

Point PhaseFinder::phase_at_T(const Phase& phase, double T) const {
  if (T >= phase.T.back()) {
    return {phase.X.back(), phase.V.back(), T};
  }
  if (T <= phase.T.front()) {
    return {phase.X.front(), phase.V.front(), T};
  }

  size_t n = std::lower_bound(phase.T.begin(), phase.T.end(), T) - phase.T.begin();
  if (n > 0) {
    n = phase.T[n] - T < T - phase.T[n - 1] ? n : n - 1;
  }

  const double dt_start = T - phase.T[n];
  std::vector<Eigen::VectorXd> X;
  std::vector<double> T_;
  std::vector<Eigen::VectorXd> dXdT;
  std::vector<double> V;
  Point jumped;
  const Point start = {phase.X[n], phase.V[n], phase.T[n]};
  auto end = trace_minimum(start, T, dt_start, &X, &T_, &dXdT, &V, &jumped);
  if (end != REACHED_T_STOP) {
    LOG(warning) << "Expected to reach tstop but end = " << end;
  }
  return {X.back(), V.back(), T};
}

std::tuple<bool, bool> PhaseFinder::redundant(const Phase& phase1, const Phase& phase2, end_descriptor end) const {
  const double tmax_1 = phase1.T.back();
  const double tmin_1 = phase1.T.front();
  const double tmax_2 = phase2.T.back();
  const double tmin_2 = phase2.T.front();

  const double tmax = std::min(tmax_1, tmax_2);
  const double tmin = std::max(tmin_1, tmin_2);

  if (tmin > tmax) {
    return std::make_tuple(false, false);
  }

  Eigen::VectorXd x1, x2;

  LOG(debug) << "Checking whether phases " << phase1.key << " and " << phase2.key << " are redundant at endpoint(s) " << end;

  bool low = false;
  bool high = false;

  if (end != LOW) {
    // Check the high temperature endpoint.
    x1 = phase_at_T(phase1, tmax).x;
    x2 = phase_at_T(phase2, tmax).x;

    high = identical_within_tol(x1, x2);
  }

  if (end != HIGH) {
    // Check the low temperature endpoint.
    x1 = phase_at_T(phase1, tmin).x;
    x2 = phase_at_T(phase2, tmin).x;
    
    low = identical_within_tol(x1, x2);
  }

  return std::make_tuple(low, high);

  /*if (end != LOW) {
    x1 = phase_at_T(phase1, tmax).x;
    x2 = phase_at_T(phase2, tmax).x;

    if (!identical_within_tol(x1, x2)) {
      return false, false;
    } else if (end == HIGH) {
      LOG(debug) << "Phase " << phase1.key << " is identical to Phase " << phase2.key << " at Tmax: " << tmax;
      return false, true;
    }
  }

  if (end != HIGH) {
    x1 = phase_at_T(phase1, tmin).x;
    x2 = phase_at_T(phase2, tmin).x;

    if (!identical_within_tol(x1, x2)) {
      return false;
    } else if (end == LOW) {
      LOG(debug) << "Phase " << phase1.key << " is identical to Phase " << phase2.key << " at Tmin: " << tmin;
      return true;
    }
  }

  return true;*/
}

void PhaseFinder::merge_phase_gaps() {
  std::vector<PhaseMerge> merges;

  LOG(debug) << "Searching for small gaps between phases...";

  // Loop through the phases and determine if these phases should be merged into the same phase. The phase tracing can
  // treat a second-order phase transition as two separate phases, with a small temperature gap between the phases. This
  // prevents some valid transitions paths from being found. Merging these equivalent phases prevents such issues.
  for (size_t i = 0; i < phases.size(); ++i) {
    // Use this as a flag for whether a phase should be deleted after merging has been performed.
    phases[i].redundant = false;
    
    for (size_t j = 0; j < phases.size(); ++j) {
      if (i != j && should_merge_phases(phases[i], phases[j])) {
        // Always merge the higher temperature phase into the lower temperature phase. This convention was chosen
        // because we observe phases that split as the Universe cools, and do not observe phases that merge as the
        // Universe cools. If the phase splits, we can merge it into the low temperature phase that is deeper just below
        // the splitting temperature, assuming they are not exactly degenerate.
        /*if (phases[i].T.back() < phases[j].T.back()) {
          merges.push_back({j, i, phases[j].T.front(), false});
        } else {
          merges.push_back({i, j, phases[i].T.front(), false});
        }*/
        // We guarantee that the merge should be i->j in should_merge_phases.
        merges.push_back({i, j, phases[i].T.front(), false});
      }
    }
  }

  // Check whether there are any phase splittings, where a single phase splits into two phases as the Universe cools.
  // This can be found by two PhaseMerges having the same 'from' phase and having very similar temperatures.
  std::vector<int> splitMultiplicity(phases.size());
  for (size_t i = 0; i < merges.size(); ++i) {
    ++splitMultiplicity[merges[i].fromPhase];
    LOG(debug) << "Identified possible merge: " << merges[i];
  }

  for (size_t i = 0; i < splitMultiplicity.size(); ++i) {
    if (splitMultiplicity[i] > 1) {
      // Need to handle this case! Find the deepest toPhase and only keep that merge. The rest should be converted to
      // subcritical transitions (which will be handled automatically by append_subcritical_transitions).
      LOG(debug) << "Detected split multiplicity " << splitMultiplicity[i] << " > 1 for Phase " << i;

      // Search through the merges to find those that have this phase (i) as the from phase. Also find the lowest energy
      // phase out of those.
      std::vector<int> relevantMerges;
      int lowestEnergyPhaseIndex;
      double lowestEnergy = std::numeric_limits<double>::max();
      double energy;
      for (size_t j = 0; j < merges.size(); ++j) {
        if (merges[j].fromPhase == i) {
          relevantMerges.push_back(j);
          LOG(debug) << "Relevant merge: " << relevantMerges.back();
        }
      }

      // Note that if find_deepest_phase fails, it returns -1, in which case we reject all relevantMerges.
      lowestEnergyPhaseIndex = find_deepest_phase(merges, relevantMerges);

      if (lowestEnergyPhaseIndex >= 0) {
        LOG(debug) << "The deepest 'to' phase is Phase " << lowestEnergyPhaseIndex;
        LOG(debug) << "Removing all other relevant merges...";
      }
      else {
        LOG(debug) << "Unable to determine the deepest 'to' phase, removing all relevant merges...";
      }

      // Now loop through other merges in relevantMerges and flag for deletion those that aren't the lowest energy one.
      for (size_t j = 0; j < relevantMerges.size(); ++j) {
        if (merges[relevantMerges[j]].toPhase != lowestEnergyPhaseIndex) {
          merges[relevantMerges[j]].rejected = true;
        }
      }
    }
  }

  // Remove the rejected merges using the erase-remove idiom.
  merges.erase(std::remove_if(merges.begin(), merges.end(), [&](const PhaseMerge& pm) {return pm.rejected;}),
    merges.end());

  // We need to sort the merges in descending temperature order so that we don't merge into a phase that will be removed
  // later. For instance if we had merges 1->2 and 2->3, doing them in the order 2->3 then 1->2 would not save the data
  // from phase 1 in phase 3, the only phase remaining at the end.
  std::sort(merges.begin(), merges.end());
  std::reverse(merges.begin(), merges.end());

  for (size_t i = 0; i < merges.size(); ++i) {
    LOG(debug) << "Performing phase merge: " << merges[i];
    perform_phase_merge(merges[i]);
  }

  // Remove the merged fromPhases using the erase-remove idiom.
  phases.erase(std::remove_if(phases.begin(), phases.end(), [&](const Phase& p) {return p.redundant;}), phases.end());

  // Update phase keys to match their position in the array.
  // TODO: will this break any assumptions throughout the code?
  for (size_t i = 0; i < phases.size(); ++i) {
    phases[i].key = i;
  }
}

bool PhaseFinder::should_merge_phases(const Phase& phase1, const Phase& phase2) {
  if (phase1.key == phase2.key) {
    return false;
  }

  double deltaT = phase1.T.front() - phase2.T.back();
  double deltaPhi = (phase1.X.front() - phase2.X.back()).norm();

  // If phase1 is the higher temperature phase, it's minimum temperature is just above the maximum temperature of
  // phase2, and the phases are close to each in field-space at their almost shared temperature, then we could merge
  // these phases. If this is not true, we can return false now.
  if (!(deltaT >= 0. && deltaT < dt_merge_phases && deltaPhi < dx_merge_phases)) {
    return false;
  }

  // Before saying that we should merge these phases, we need to check if there are any phases existing between the
  // two merge candidates. If a phase is present between the merge candidates, the merging should not take place. It
  // is not correct to merge across another phase. Check at both ends of the temperature interval where the gap exists.
  for (const Phase& phaseMid : phases) {
    if (phaseMid.contains_t(phase1.T.front())) {
        Eigen::VectorXd x = phase_at_T(phaseMid, phase1.T.front()).x;

        if ((phase1.X.front() - x).norm() < deltaPhi && (phase2.X.back() - x).norm() < deltaPhi) {
          return false;
        }
    }

    if (phaseMid.contains_t(phase2.T.back())) {
      Eigen::VectorXd x = phase_at_T(phaseMid, phase2.T.back()).x;

      if ((phase1.X.front() - x).norm() < deltaPhi && (phase2.X.back() - x).norm() < deltaPhi) {
        return false;
      }
    }
  }

  return true;
}

int PhaseFinder::find_deepest_phase(const std::vector<PhaseMerge>& merges, const std::vector<int>& relevantMerges) {
  // Search through the high-T energy of the toPhases in relevantMerges to determine which phase would be transitioned
  // to. The expected scenario is that a phase splits into multiple phases as the Universe cools, and this splitting
  // should occur at a particular temperature. Thus we expect the low-T phases to all have a very similar maximum
  // temperature and it is then safe to evaluate all low-T phases at the lowest maximum temperature (where they should
  // all coexist).

  double minTemp = std::numeric_limits<double>::max();
  int minTempIndex = -1;
  double toPhaseMaxTemp;

  for (size_t i = 0; i < relevantMerges.size(); ++i) {
    toPhaseMaxTemp = phases[merges[relevantMerges[i]].toPhase].T.back();

    if (toPhaseMaxTemp < minTemp) {
      if (minTempIndex >= 0 && toPhaseMaxTemp < minTemp*0.95) {
        LOG(debug) << "When finding deepest phase for merging, encountered sizable gap between maximum temperature of "
          "'to' phases: Phase " << minTempIndex << " has Tmax=" << phases[minTempIndex].T.back() << " and Phase " <<
          merges[relevantMerges[i]].toPhase << " has Tmax=" << toPhaseMaxTemp;
      }

      // Make sure that all relevant phases exist at this temperature. If a phase doesn't exist at this new lowest
      // temperature, keep the old temperature and ignore this merge. The toPhase for this merge is assumed to be
      // not directly accessible from the fromPhase through a SOPT, since a SOPT would occur to a higher temperature
      // toPhase first.
      for (size_t j = 0; j < relevantMerges.size(); ++j) {
        if (phases[merges[relevantMerges[j]].toPhase].T.front() > minTemp) {
          LOG(debug) << "When finding deepest phase for merging, encountered toPhase (" <<
            merges[relevantMerges[i]].toPhase << ") with a maximum temperature (Tmax=" << toPhaseMaxTemp << ") where "
            "another toPhase (" << merges[relevantMerges[j]].toPhase << ") no longer exists.";
          continue;
        }
      }

      minTemp = toPhaseMaxTemp;
      minTempIndex = merges[relevantMerges[i]].toPhase;
    }
  }

  if (minTempIndex == -1) {
    // TODO: Hopefully this won't happen! Needs testing.
    std::cerr << "Unable to find deepest phase to merge to from Phase " << merges[relevantMerges[0]].fromPhase <<
        std::endl;
    return -1;
  }

  // We now have the temperature at which we should sample the toPhases, to compare their energies. It is now a matter
  // of finding the lowest energy toPhase at that sample temperature.
  double minEnergy = std::numeric_limits<double>::max();
  // This is the index in the phases array corresponding to the lowest energy phase out of the toPhases from the
  // relevantMerges.
  int minEnergyIndex = -1;
  double toPhaseEnergy;
  //Phase& toPhase;

  // Find which toPhase has the lowest energy at minTemp.
  for (size_t i = 0; i < relevantMerges.size(); ++i) {
    Phase* toPhase = &phases[merges[relevantMerges[i]].toPhase];

    // Ignore any phases that don't exist at this temperature (we are already guaranteed that T.front() is smaller).
    if (toPhase->T.back() < minTemp) {
      continue;
    }
    
    toPhaseEnergy = phase_at_T(*toPhase, minTemp).potential;

    if (toPhaseEnergy < minEnergy) {
      minEnergy = toPhaseEnergy;
      minEnergyIndex = merges[relevantMerges[i]].toPhase;
    }
  }

  return minEnergyIndex;
}

void PhaseFinder::perform_phase_merge(const PhaseMerge& merge) {
  if (phases[merge.fromPhase].redundant) {
    LOG(debug) << "While performing merge " << merge << ", fromPhase " << merge.fromPhase << " is already redundant!";
  }

  if (phases[merge.toPhase].redundant) {
    LOG(debug) << "While performing merge " << merge << ", toPhase " << merge.toPhase << " is redundant!";
  }

  Phase& fromPhase = *&phases[merge.fromPhase];
  Phase& toPhase = *&phases[merge.toPhase];

  toPhase.X.insert(toPhase.X.end(), fromPhase.X.begin(),
    fromPhase.X.end());
  toPhase.T.insert(toPhase.T.end(), fromPhase.T.begin(),
    fromPhase.T.end());
  toPhase.dXdT.insert(toPhase.dXdT.end(), fromPhase.dXdT.begin(),
    fromPhase.dXdT.end());
  toPhase.V.insert(toPhase.V.end(), fromPhase.V.begin(),
    fromPhase.V.end());
  
  toPhase.end_high = fromPhase.end_low;
  
  fromPhase.redundant = true;
}

void PhaseFinder::remove_redundant() {
  // Loop as we want to do as many passes as necessary until phases don't change
  bool changed = true;
  while (changed) {
    changed = false;

    for (auto &phase1 : phases) {
      for (auto &phase2 : phases) {
        // Don't investigate transitions more than once
        if (phase1.key >= phase2.key) {
          continue;
        }

        // Don't re-check any redundant phases
        if (phase1.redundant || phase2.redundant) {
          continue;
        }

        LOG(debug) << "Checking redundancy between phases " << phase1.key << " and " << phase2.key;

        bool isRedundantLow, isRedundantHigh;
        std::tie(isRedundantLow, isRedundantHigh) = redundant(phase1, phase2);

        // If the phases have the same field location at both the high and low temperatures that they share, mark one as redundant.
        if (isRedundantLow && isRedundantHigh) {
          LOG(debug) << "Phases " << phase1.key << " and " << phase2.key
                     << " are redundant";
          changed = true;
          const auto low_key = phase1.T.front() <= phase2.T.front() ? phase1.key : phase2.key;
          const auto high_key = phase1.T.back() >= phase2.T.back() ? phase1.key : phase2.key;

          if (low_key == high_key) {
            LOG(debug) << "One phase completely contained the other";
            const auto redundant_key = phase1.key == low_key ? phase2.key : phase1.key;
            phases[redundant_key].redundant = true;
            LOG(debug) << "Kept phase " << low_key
                       << " and marked redundant " << redundant_key;

          } else {
            LOG(debug) << "One phase did not completely contain the other";

            // Add additional points from high phase to low phase

            const double tmax = std::min(phase1.T.back(), phase2.T.back());
            const size_t n = std::upper_bound(phases[high_key].T.begin(), phases[high_key].T.end(), tmax) - phases[high_key].T.begin();
            phases[low_key].X.insert(phases[low_key].X.end(), phases[high_key].X.begin() + n, phases[high_key].X.end());
            phases[low_key].T.insert(phases[low_key].T.end(), phases[high_key].T.begin() + n, phases[high_key].T.end());
            phases[low_key].dXdT.insert(phases[low_key].dXdT.end(), phases[high_key].dXdT.begin() + n, phases[high_key].dXdT.end());
            phases[low_key].V.insert(phases[low_key].V.end(), phases[high_key].V.begin() + n, phases[high_key].V.end());

            phases[high_key].redundant = true;

            LOG(debug) << "Merged them into phase " << low_key << " and marked "
                       << high_key << " redundant";
          }
        }
        else if (isRedundantLow || isRedundantHigh) {
            split_overlapping_phases(phase1, phase2, isRedundantLow, isRedundantHigh);
        }
      }
    }
  }

  std::vector<int> prunedKeys;

  LOG(debug) << "Tracking keys of redundant phases...";

  // Track the keys of the redundant phases so we know how to decrement the surviving phases' keys after pruning.
  for(const auto &phase : phases) {
    if(phase.redundant) {
      prunedKeys.push_back(phase.key);
    }
  }

  LOG(debug) << "Pruning redundant phases...";

  // Finally, prune redundant phases
  phases.erase(std::remove_if(phases.begin(), phases.end(), [](Phase phase) {return phase.redundant;}), phases.end());

  LOG(debug) << "Determing amount to decrement remaining keys...";
  std::vector<int> decrementAmount(phases.size(), 0);

  // Determine how much each phase key needs to be decremented based on the keys of removed phases.
  for(int i = 0; i < prunedKeys.size(); ++i) {
    for(int j = 0; j < phases.size(); ++j) {
      if(phases[j].key > prunedKeys[i]) {
        ++decrementAmount[j];
      }
    }
  }

  LOG(debug) << "Decrementing keys...";

  // Decrement phase keys that come after any removed phase keys.
  for(int i = 0; i < phases.size(); ++i) {
    phases[i].key -= decrementAmount[i];
  }

  LOG(debug) << "Done.";
}

void PhaseFinder::split_overlapping_phases(Phase& phase1, Phase& phase2, bool isRedundantLow, bool isRedundantHigh) {
  if (isRedundantLow && isRedundantHigh) {
    LOG(debug) << "Attempted to split overlapping phases that are redundant at both temperature endpoints";
    return;
  }

  if (isRedundantLow) {
    std::cerr << "Detected two phases merging as temperature decreases! This is not yet supported and requires further"
      " investigation." << std::endl;
    return;
  }

  // Do the splitting. Determine the shared temperature interval, then bisect it to find where the phases split.
  double Tmax = std::min(phase1.T.back(), phase2.T.back());
  double Tmin = std::max(phase1.T.front(), phase2.T.front());
  double Tmid;
  Point point1, point2, lastNonIdenticalPoint1, lastNonIdenticalPoint2;
  bool identical;
  double lastNonIdenticalT = Tmin;

  LOG(debug) << "Splitting overlapping phases " << phase1.key << " and " << phase2.key;
  LOG(debug) << "Tmin: " << Tmin << ", Tmax: " << Tmax;

  // TODO: this was an attempt to get the optimiser to not jump to the nearby deeper phase, but it still doesn't work...
  double minimisationStepSize = get_find_min_trace_abs_step();
  set_find_min_trace_abs_step(0.01*minimisationStepSize);

  while (Tmax - Tmin > dt_min_rel_split_phase)
  {
    Tmid = 0.5*(Tmin + Tmax);
    point1 = phase_at_T(phase1, Tmid);
    point2 = phase_at_T(phase2, Tmid);
    identical = identical_within_tol(point1.x*10.0, point2.x*10.0);
    LOG(debug) << "Tmid: " << Tmid << ", dist: " << (point1.x - point2.x).norm() << ", identical: " << identical;

    if (identical)
    {
      Tmax = Tmid;
    }
    else
    {
      Tmin = Tmid;
      lastNonIdenticalT = Tmid;
      lastNonIdenticalPoint1.x = point1.x;
      lastNonIdenticalPoint1.t = point1.t;
      lastNonIdenticalPoint1.potential = point1.potential;
      lastNonIdenticalPoint2.x = point2.x;
      lastNonIdenticalPoint2.t = point2.t;
      lastNonIdenticalPoint2.potential = point2.potential;
    }
  }

  set_find_min_trace_abs_step(minimisationStepSize);

  LOG(debug) << "Last non-identical T: " << lastNonIdenticalT;

  /*if (Tmid != lastNonIdenticalT)
  {
    point1 = phase_at_T(phase1, lastNonIdenticalT);
    point2 = phase_at_T(phase2, lastNonIdenticalT);
  }*/

  // We now have a good estimate of where the phase splitting occurs. Now we want to split the phases.
  // Currently we have two phases that overlap at high temperature and diverge at some lower temperature.
  // Just like we handle phase splitting when we merge phase gaps, we want to determine which phase is
  // deeper just below the split point.
  Phase& shallowPhase = lastNonIdenticalPoint1.potential < lastNonIdenticalPoint2.potential ? phase2 : phase1;
  Phase& deeperPhase = lastNonIdenticalPoint1.potential < lastNonIdenticalPoint2.potential ? phase1 : phase2;
  Point& shallowPoint = lastNonIdenticalPoint1.potential < lastNonIdenticalPoint2.potential ? lastNonIdenticalPoint2 : lastNonIdenticalPoint1;
  Point& deeperPoint = lastNonIdenticalPoint1.potential < lastNonIdenticalPoint2.potential ? lastNonIdenticalPoint1 : lastNonIdenticalPoint2;

  LOG(debug) << "Shallow point: " << shallowPoint;
  LOG(debug) << "Deep point:    " << deeperPoint;

  // This method isn't sufficient to guarantee we pick the correct shallow and deeper phase. Due to numerical errors in
  // both the phase locations and the potential at those locations, it is possible to make the wrong choice. We should
  // double-check at a slightly lower temperature to ensure we have made the correct choice.

  LOG(debug) << "Identified phase " << shallowPhase.key << " as the shallow phase.";
  LOG(debug) << "Identified phase " << deeperPhase.key << " as the deep phase.";

  // The deeper phase should be left mostly unchanged. However, we want to add this new phase sample to the phases's
  // stored data. This is necessary to avoid spurious transitions between these two phases directly below the splitting.
  // The shallower phase will have precise data at the split point, including the potential energy, whereas the deeper
  // phase would otherwise rely on interpolation between samples surrounding this split point. This can lead to the
  // shallow phase appearing to have lower energy for a short temperature interval (which is not correct).

  // Find the insertion point in the deeper phase's data arrays, where this new sample should be inserted.
  int i = 0;

  LOG(debug) << "============================================================================";
  LOG(debug) << "Comparing against lastNonIdenticalT = " << lastNonIdenticalT;

  while (deeperPhase.T[i] < lastNonIdenticalT)
  {
    LOG(debug) << "T[" << i << "] = " << deeperPhase.T[i];
    ++i;
  }

  // i now contains the first index above the split point. This is precisely the element we want to insert before.
  // Insert the new deeper phase sample at this position in the data arrays of the deeper phase.
  deeperPhase.T.insert(deeperPhase.T.begin()+i, deeperPoint.t);
  deeperPhase.X.insert(deeperPhase.X.begin()+i, deeperPoint.x);
  deeperPhase.V.insert(deeperPhase.V.begin()+i, deeperPoint.potential);
  deeperPhase.dXdT.insert(deeperPhase.dXdT.begin()+i, dx_min_dt(deeperPoint.x, deeperPoint.t));

  LOG(debug) << "Adding point {T: " << deeperPhase.T[i] << ", X: " << deeperPhase.X[i] << ", V: " <<
    deeperPhase.V[i] << ", dXdT: " << deeperPhase.dXdT[i] << "} to the deep phase.";
  
  // The shallower phase should have all of the overlapping data removed. The split will be handled by merge_phase_gaps.
  i = 0;

  // Find the split point in the shallow phase's data arrays.
  while (shallowPhase.T[i] < lastNonIdenticalT)
  {
    ++i;
  }

  // i now contains the first index above the split point. Subtract 1 to find the last index below the split point.
  --i;

  // Remove all data stored in the shallow phase above this split point. This is the overlap region with the deeper phase.
  shallowPhase.T.resize(i+1);
  shallowPhase.X.resize(i+1);
  shallowPhase.V.resize(i+1);
  shallowPhase.dXdT.resize(i+1);

  // Add the closest sampled point below the split to the shallow phase. This provides the phase with accurate
  // knowledge of its high temperature endpoint, which might otherwise be significantly below the split.
  //Point closestPointBelowSplit = point1.potential < point2.potential ? point1 : point2;
  shallowPhase.T.push_back(shallowPoint.t);
  shallowPhase.X.push_back(shallowPoint.x);
  shallowPhase.V.push_back(shallowPoint.potential);
  shallowPhase.dXdT.push_back(dx_min_dt(shallowPoint.x, shallowPoint.t));

  LOG(debug) << "Adding point {T: " << shallowPhase.T.back() << ", X: " << shallowPhase.X.back() << ", V: " <<
    shallowPhase.V.back() << ", dXdT: " << shallowPhase.dXdT.back() << "} to the shallow phase.";

  // TODO: should update the HIGH phase end descriptor for the phase.
}

bool PhaseFinder::out_of_bounds(const Eigen::VectorXd& x) const {
  std::vector<double> vector_x(x.data(), x.data() + x.rows() * x.cols());
  return (vector_x <= lower_bounds) || (vector_x >= upper_bounds);
}


}  // namespace PhaseTracer
