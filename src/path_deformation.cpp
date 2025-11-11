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

#include "path_deformation.hpp"
#include <memory>

namespace PhaseTracer {

double max_norm(std::vector<Eigen::VectorXd> vec) {
  double max_sq = 0.;
  for (const auto &v : vec) {
    max_sq = std::max(max_sq, v.squaredNorm());
  }
  return std::sqrt(max_sq);
}

SplinePath::SplinePath(EffectivePotential::Potential &potential,
                       double T_,
                       std::vector<Eigen::VectorXd> pts_,
                       bool extend_to_minima_,
                       bool reeval_distances_) : P(potential), T(T_), pts(pts_), nphi(potential.get_n_scalars()), num_nodes(pts_.size()),
                                                 extend_to_minima(extend_to_minima_), reeval_distances(reeval_distances_) {

  // 1. Find derivs
  std::vector<Eigen::VectorXd> dpts = _pathDeriv(pts);
  // 2. Extend the path
  if (extend_to_minima) {
    double xmin = find_loc_min_w_guess(pts[0], dpts[0]); // TODO: This may return global mim instead of local minimum
    xmin = std::min(xmin, 0.0);
    int nx = static_cast<int>(std::ceil(std::abs(xmin) - .5)) + 1;
    if (nx > 1) {
      double stepSize = -xmin / (nx - 1);
      std::vector<Eigen::VectorXd> new_pts;
      for (int i = 0; i < nx; ++i) {
        double x = xmin + i * stepSize;
        Eigen::VectorXd pt_ext = pts[0] + x * dpts[0];
        new_pts.push_back(pt_ext);
      }
      new_pts.insert(new_pts.end(), pts.begin() + 1, pts.end());
      pts = new_pts;
    }

    xmin = find_loc_min_w_guess(pts.back(), dpts.back());
    xmin = std::max(xmin, 0.0);
    nx = static_cast<int>(std::ceil(std::abs(xmin) - 0.5)) + 1;
    if (nx > 1) {
      double stepSize = -xmin / (nx - 1);
      std::vector<Eigen::VectorXd> new_pts2;
      for (int i = 0; i < nx; ++i) {
        double x = xmin + (nx - i - 1) * stepSize;
        Eigen::VectorXd pt_ext = pts.back() + x * dpts.back();
        new_pts2.push_back(pt_ext);
      }
      pts.pop_back();
      pts.insert(pts.end(), new_pts2.begin(), new_pts2.end());
    }
    dpts = _pathDeriv(pts);
    num_nodes = pts.size();
  }
  // 3. Find knot positions and fit the spline
  std::vector<double> squared_sums;
  for (const auto &vec : dpts) {
    squared_sums.push_back(vec.norm());
  }

  std::vector<double> pdist = cumulative_trapezoidal_integration(squared_sums);
  pdist.insert(pdist.begin(), 0.0);
  length = pdist[pdist.size() - 1];
  set_path_tck(pdist);

  // 4. Re-evaluate the distance to each point.
  if (reeval_distances) {
    std::vector<double> pdist_(pdist.size());
    typedef double state_type;
    state_type x = 0.0;
    double t_start = 0.0;
    double dt = pdist.back() * 1E-4;
    using error_stepper_type =
        boost::numeric::odeint::runge_kutta_dopri5<state_type>;
    using controlled_stepper_type =
        boost::numeric::odeint::controlled_runge_kutta<error_stepper_type>;
    controlled_stepper_type stepper = make_controlled(0., pdist.back() * 1E-8, error_stepper_type());
    for (size_t i = 1; i < pdist.size(); ++i) {
      boost::numeric::odeint::integrate_adaptive(stepper, [this](const state_type &y, state_type &dydr, double r) { dpdx(y, dydr, r); }, x, pdist[i - 1], pdist[i], dt);
      pdist_[i] = x;
    }
    pdist = pdist_;
    length = pdist[pdist.size() - 1];
    set_path_tck(pdist);
  }

  // Make the potential spline.
  alglib::real_1d_array p_arr;
  alglib::real_1d_array v_arr;
  p_arr.setlength(V_spline_samples);
  v_arr.setlength(V_spline_samples);
  for (size_t i = 0; i < V_spline_samples; ++i) {
    // extend 20% beyond
    p_arr[i] = -0.2 * length + i * (1.4 * length) / (V_spline_samples - 1);
    v_arr[i] = P.V(vecp(p_arr[i]), T);
  }
  alglib::spline1dbuildcubic(p_arr, v_arr, _V_tck);
}

void SplinePath::update_path(std::vector<Eigen::VectorXd> pts_new, bool full_rebuild) {
  if (full_rebuild || extend_to_minima) {
    // Do full reconstruction for significant changes
    pts = pts_new;
    num_nodes = pts.size();
    
    std::vector<Eigen::VectorXd> dpts = _pathDeriv(pts);
    
    std::vector<double> squared_sums;
    for (const auto &vec : dpts) {
      squared_sums.push_back(vec.norm());
    }
    
    std::vector<double> pdist = cumulative_trapezoidal_integration(squared_sums);
    pdist.insert(pdist.begin(), 0.0);
    length = pdist[pdist.size() - 1];
    set_path_tck(pdist);
    
    // Update potential spline
    alglib::real_1d_array p_arr;
    alglib::real_1d_array v_arr;
    p_arr.setlength(V_spline_samples);
    v_arr.setlength(V_spline_samples);
    for (size_t i = 0; i < V_spline_samples; ++i) {
      p_arr[i] = -0.2 * length + i * (1.4 * length) / (V_spline_samples - 1);
      v_arr[i] = P.V(vecp(p_arr[i]), T);
    }
    alglib::spline1dbuildcubic(p_arr, v_arr, _V_tck);
  } else {
    // Fast update: just refit splines without extension/re-evaluation
    pts = pts_new;
    num_nodes = pts.size();
    
    std::vector<Eigen::VectorXd> dpts = _pathDeriv(pts);
    std::vector<double> squared_sums;
    for (const auto &vec : dpts) {
      squared_sums.push_back(vec.norm());
    }
    
    std::vector<double> pdist = cumulative_trapezoidal_integration(squared_sums);
    pdist.insert(pdist.begin(), 0.0);
    length = pdist[pdist.size() - 1];
    set_path_tck(pdist);
    
    // Update potential spline with fewer samples for speed
    size_t fast_samples = std::min(V_spline_samples / 2, static_cast<int>(50));
    alglib::real_1d_array p_arr;
    alglib::real_1d_array v_arr;
    p_arr.setlength(fast_samples);
    v_arr.setlength(fast_samples);
    for (size_t i = 0; i < fast_samples; ++i) {
      p_arr[i] = -0.2 * length + i * (1.4 * length) / (fast_samples - 1);
      v_arr[i] = P.V(vecp(p_arr[i]), T);
    }
    alglib::spline1dbuildcubic(p_arr, v_arr, _V_tck);
  }
}

std::vector<Eigen::VectorXd> SplinePath::_pathDeriv(const std::vector<Eigen::VectorXd> phi) { // rename phi
  int num_phi = phi.size();
  std::vector<Eigen::VectorXd> dphi(num_phi);
  if (num_phi < 2) {
    throw std::runtime_error("The number of points that describe the path must be larger than 1.");
  } else if (num_phi == 2) {
    dphi[0] = dphi[1] = phi[1] - phi[0];
  } else if (num_phi < 5) {
    // 1st/2nd order calculation
    dphi[0] = -1.5 * phi[0] + 2.0 * phi[1] - 0.5 * phi[2];
    dphi[num_phi - 1] = 1.5 * phi[num_phi - 1] - 2.0 * phi[num_phi - 2] + 0.5 * phi[num_phi - 3];
    for (int i = 1; i < num_phi - 1; ++i) {
      dphi[i] = 0.5 * (phi[i + 1] - phi[i - 1]);
    }
  } else {
    dphi = deriv14_const_dx(phi);
  }
  return dphi;
}

double SplinePath::find_loc_min_w_guess(Eigen::VectorXd p0, Eigen::VectorXd dp0, double guess) {

  std::shared_ptr<std::function<double(const std::vector<double> &)>> V_lin = std::make_shared<std::function<double(const std::vector<double> &)>>([this, p0, dp0](const std::vector<double> &x) {
    return this->P.V(p0 + x[0] * dp0, T);
  });

  auto V_lin_func = [](const std::vector<double> &x, std::vector<double> &grad, void *data) -> double {
    auto v_lin = reinterpret_cast<std::shared_ptr<std::function<double(const std::vector<double> &)>> *>(data);
    return (*(*v_lin))(x);
  };

  nlopt::opt optimizer(nlopt::LN_SBPLX, 1);
  optimizer.set_min_objective(V_lin_func, &V_lin);
  optimizer.set_xtol_abs(1e-6);

  std::vector<double> xmin(1, guess);
  double fmin;
  optimizer.optimize(xmin, fmin);
  return xmin[0];
}

void SplinePath::set_path_tck(std::vector<double> pdist) {
  alglib::real_1d_array x_arr;
  alglib::real_1d_array y_arr;
  x_arr.setlength(num_nodes);
  y_arr.setlength(num_nodes);
  _path_tck.clear();
  for (size_t i = 0; i < num_nodes; i++)
    x_arr[i] = pdist[i];
  for (size_t j = 0; j < nphi; j++) {
    alglib::spline1dinterpolant spline;
    for (size_t i = 0; i < num_nodes; i++)
      y_arr[i] = pts[i](j);
    alglib::spline1dbuildcubic(x_arr, y_arr, spline);
    _path_tck.push_back(spline);
  }
}

void SplinePath::dpdx(const double &x, double &dpdx_, double r) const {
  double norm_sq = 0.;
  double s, ds, d2s;
  for (int i = 0; i < nphi; i++) {
    alglib::spline1ddiff(_path_tck[i], x, s, ds, d2s);
    norm_sq += ds * ds;
  }
  dpdx_ = std::sqrt(norm_sq);
}

Eigen::VectorXd SplinePath::vecp(double x) const {
  Eigen::VectorXd vec(nphi);
  for (int i = 0; i < nphi; i++) {
    vec[i] = alglib::spline1dcalc(_path_tck[i], x);
  }
  return vec;
}

double SplinePath::V(double x) const {
  return alglib::spline1dcalc(_V_tck, x);
}

double SplinePath::dV(double x) const {
  double s, ds, d2s;
  alglib::spline1ddiff(_V_tck, x, s, ds, d2s);
  return ds;
}

double SplinePath::d2V(double x) const {
  double s, ds, d2s;
  alglib::spline1ddiff(_V_tck, x, s, ds, d2s);
  return d2s;
}

bool PathDeformation::deformPath(std::vector<double> dphidr) {
  num_steps = 0;
  phi_list.clear();
  F_list.clear();
  // Clear gradient cache when starting new deformation
  cached_gradients.clear();
  cached_phi_for_gradients.clear();
  
  // convert phi to a set of path lengths
  std::vector<double> dL;
  for (int i = 0; i < num_nodes - 1; i++) {
    dL.push_back((phi_node[i + 1] - phi_node[i]).norm());
  }
  t_node.clear();
  t_node.push_back(0);
  for (int i = 0; i < dL.size(); i++) {
    t_node.push_back(t_node[i] + dL[i]);
  }
  totalLength_node = t_node.back();
  for (int i = 0; i < t_node.size(); i++) {
    t_node[i] /= totalLength_node;
  }
  t_node[0] = 1e-50;

  // create the starting spline
  std::vector<double> t0;
  for (int i = 0; i < kb - 1; ++i) {
    t0.push_back(0.0);
  }
  for (int i = 0; i < nb; ++i) {
    t0.push_back(i / (nb - 1.));
  }
  for (int i = 0; i < kb - 1; ++i) {
    t0.push_back(1.0);
  }

  auto result = Nbspld2(t0, t_node, kb);
  X_node = std::get<0>(result);
  dX_node = std::get<1>(result);
  d2X_node = std::get<2>(result);
  
  // Cache SVD decomposition since X_node won't change during iterations
  X_node_svd.compute(X_node, Eigen::ComputeThinU | Eigen::ComputeThinV);
  svd_cached = true;

  Eigen::VectorXd phi0 = phi_node[0];
  Eigen::VectorXd phi1 = phi_node.back();
  beta_node.clear();
  for (int j = 0; j < nphi; ++j) {
    Eigen::VectorXd phi_delta(num_nodes);
    for (int i = 0; i < num_nodes; ++i) {
      Eigen::VectorXd phii = phi_node[i];
      double phi_lin = phi0[j] + (phi1[j] - phi0[j]) * t_node[i];
      phi_delta[i] = phii[j] - phi_lin;
    }
    // Use cached SVD if available
    if (svd_cached) {
      beta_node.push_back(X_node_svd.solve(phi_delta));
    } else {
      beta_node.push_back(X_node.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(phi_delta));
    }
  }

  double max_v2 = 0;
  for (int i = 0; i < num_nodes; ++i) {
    max_v2 = std::max(max_v2, P.dV_dx(phi_node[i], T).norm());
  }
  v2min *= max_v2 * totalLength_node / nb;

  v2_node.clear();
  for (int i = 0; i < num_nodes; ++i) {
    v2_node.push_back(std::max(v2min, dphidr[i] * dphidr[i]));
  }

  // Deform the path
  double minfRatio = std::numeric_limits<double>::infinity();
  size_t minfRatio_index = 0;
  std::vector<Eigen::VectorXd> minfRatio_beta;
  std::vector<Eigen::VectorXd> minfRatio_phi;
  double stepsize = startstep;
  bool deformation_converged = false;
  bool step_reversed;
  double fRatio;
  phi_prev.clear();
  F_prev.clear();
  while (true) {
    num_steps++;
    step(stepsize, step_reversed, fRatio);
    // TODO: add a user-defined callback function

    if (fRatio < fRatioConv || (num_steps == 1 && fRatio < converge_0 * fRatioConv)) {
      LOG(debug) << "Path deformation converged. "
                 << num_steps << " steps. fRatio = " << fRatio;
      deformation_converged = true;
      break;
    }

    if (minfRatio > fRatio) {
      minfRatio = fRatio;
      minfRatio_beta = beta_node;
      minfRatio_index = num_steps;
      minfRatio_phi = phi_node;
    }

    if (fRatio > fRatioIncrease * minfRatio && !step_reversed) {
      throw std::runtime_error("Deformation not converging");
    }

    if (num_steps >= step_maxiter) {
      throw std::runtime_error("Maximum number of deformation iterations reached");
    }
  }
  return deformation_converged;
}

void PathDeformation::step(double &lastStep, bool &step_reversed, double &fRatio) {

  std::vector<Eigen::VectorXd> _phi = phi_node;
  std::vector<Eigen::VectorXd> F, dV;
  forces(F, dV);
  double F_max = max_norm(F);
  double dV_max = max_norm(dV);
  double fRatio1 = F_max / dV_max;
  for (int ii = 0; ii < dX_node.rows(); ++ii)
    F[ii] *= totalLength_node / dV_max;

  // see how big the stepsize should be
  double stepsize = lastStep;
  step_reversed = false;
  std::vector<double> FdotFlast(dX_node.rows());
  if (reverseCheck < 1 and (not F_prev.empty())) {
    for (int ii = 0; ii < dX_node.rows(); ++ii)
      FdotFlast[ii] = F[ii].dot(F_prev[ii]);
    int numNegative = std::count_if(FdotFlast.begin(), FdotFlast.end(),
                                    [](double val) { return val < 0; });
    if (numNegative > num_nodes * reverseCheck) {
      if (stepsize > minstep) {
        step_reversed = true;
        _phi = phi_prev;
        F = F_prev;
        LOG(trace) << "step reversed";
        stepsize = lastStep / stepDecrease;
      }
    } else {
      stepsize = lastStep * stepIncrease;
    }
  }

  if (stepsize > maxstep)
    stepsize = maxstep;
  if (stepsize < minstep)
    stepsize = minstep;

  // Save the state before the step
  phi_prev = _phi;
  F_prev = F;
  if (save_all_steps) {
    phi_list.push_back(_phi);
    F_list.push_back(F);
  }

  for (int ii = 0; ii < dX_node.rows(); ++ii)
    _phi[ii] += F[ii] * stepsize;

  // fit to the spline
  std::vector<Eigen::VectorXd> phi_lin(dX_node.rows());
  Eigen::VectorXd delta_phi = _phi.back() - _phi[0];
  Eigen::VectorXd phi0 = _phi[0];
  for (int ii = 0; ii < dX_node.rows(); ++ii) {
    phi_lin[ii] = phi0 + delta_phi * t_node[ii];
    _phi[ii] -= phi_lin[ii];
  }

  //    std::cout << "betas:" << beta << std::endl;
  for (size_t ii = 0; ii < nphi; ++ii) {
    Eigen::VectorXd ii_phi(num_nodes);
    for (size_t jj = 0; jj < num_nodes; ++jj) {
      ii_phi(jj) = _phi[jj](ii);
    }
    // Use cached SVD if available
    if (svd_cached) {
      beta_node[ii] = X_node_svd.solve(ii_phi);
    } else {
      beta_node[ii] = X_node.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(ii_phi);
    }
  }

  double fRatio2 = 0.;
  for (int ii = 0; ii < dX_node.rows(); ++ii) {
    Eigen::VectorXd ii_phi(nphi);
    for (int jj = 0; jj < nphi; ++jj) {
      ii_phi[jj] = 0.;
      for (int kk = 0; kk < dX_node.cols(); ++kk) {
        ii_phi[jj] += beta_node[jj][kk] * X_node(ii, kk);
      }
    }
    _phi[ii] = ii_phi + phi_lin[ii];

    Eigen::VectorXd Ffit = (_phi[ii] - phi_prev[ii]) / stepsize;
    fRatio2 = std::max(fRatio2, Ffit.norm() / totalLength_node);
  }
  phi_node = _phi;
  lastStep = stepsize;

  LOG(trace) << "step: " << num_steps
             << "; stepsize: " << stepsize
             << "; fRatio1: " << fRatio1
             << "; fRatio2: " << fRatio2;
  fRatio = checkAfterFit ? fRatio2 : fRatio1;
}

// Calculate the normal force and potential gradient on the path
void PathDeformation::forces(std::vector<Eigen::VectorXd> &F_norm, std::vector<Eigen::VectorXd> &dV) {
  F_norm.resize(dX_node.rows());
  dV.resize(dX_node.rows());
  
  // Check if we can reuse cached gradients
  bool use_cache = (cached_gradients.size() == num_nodes && 
                    cached_phi_for_gradients.size() == num_nodes);
  if (use_cache) {
    for (int ii = 0; ii < num_nodes; ++ii) {
      if (!phi_node[ii].isApprox(cached_phi_for_gradients[ii], 1e-10)) {
        use_cache = false;
        break;
      }
    }
  }
  
  // Compute or reuse gradients
  if (!use_cache) {
    cached_gradients.resize(num_nodes);
    cached_phi_for_gradients.resize(num_nodes);
    for (int ii = 0; ii < num_nodes; ++ii) {
      cached_gradients[ii] = P.dV_dx(phi_node[ii], T);
      cached_phi_for_gradients[ii] = phi_node[ii];
    }
  }

  for (int ii = 0; ii < num_nodes; ++ii) {
    Eigen::VectorXd dphi(nphi);
    Eigen::VectorXd d2phi(nphi);
    for (int jj = 0; jj < nphi; ++jj) {
      dphi[jj] = 0.;
      d2phi[jj] = 0.;
      for (int kk = 0; kk < nb; ++kk) {
        dphi[jj] += beta_node[jj][kk] * dX_node(ii, kk);
        d2phi[jj] += beta_node[jj][kk] * d2X_node(ii, kk);
      }
      dphi[jj] += phi_node.back()[jj] - phi_node[1][jj]; // Note this is phi[1], not phi[0]
    }

    double dphi_ = dphi.norm();
    double dphi_sq = dphi_ * dphi_;
    Eigen::VectorXd dphids = dphi / dphi_;
    Eigen::VectorXd d2phids2 = (d2phi - dphi * (dphi.dot(d2phi)) / dphi_sq) / dphi_sq;

    dV[ii] = cached_gradients[ii];
    Eigen::VectorXd dV_perp = dV[ii] - dV[ii].dot(dphids) * dphids;
    F_norm[ii] = d2phids2 * v2_node[ii] - dV_perp;
  }
  if (fix_start) {
    F_norm[0] = Eigen::VectorXd::Zero(nphi);
    ;
  }
  if (fix_end) {
    F_norm.back() = Eigen::VectorXd::Zero(nphi);
  }
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> PathDeformation::Nbspld2(std::vector<double> t, std::vector<double> x, int k) const {
  int kmax = k;
  if (kmax >= t.size() - 2) {
    throw std::runtime_error("Input error in Nbspl: require that k < len(t)-2");
  }

  int size_x = static_cast<int>(x.size());

  Eigen::MatrixXd N = Eigen::MatrixXd::Ones(size_x, t.size() - 1);
  Eigen::MatrixXd dN = Eigen::MatrixXd::Zero(size_x, t.size() - 1);
  Eigen::MatrixXd d2N = Eigen::MatrixXd::Zero(size_x, t.size() - 1);

  for (int i = 0; i < N.rows(); ++i) {
    for (int j = 0; j < N.cols(); ++j) {
      N(i, j) = 1.0 * ((x[i] > t[j]) && (x[i] <= t[j + 1]));
    }
  }

  for (int i = 1; i <= kmax; ++i) {
    Eigen::VectorXd dt(t.size() - i);
    Eigen::VectorXd _dt(t.size() - i);

    for (int j = 0; j < t.size() - i; ++j) {
      dt(j) = t[j + i] - t[j];
      _dt(j) = (dt(j) != 0) ? 1.0 / dt(j) : 0.0;
    }
    for (int j = 0; j < size_x; ++j) {
      for (int l = 0; l < dt.size() - 1; ++l) {

        d2N(j, l) = d2N(j, l) * (x[j] - t[l]) * _dt[l] - d2N(j, 1 + l) * (x[j] - t[1 + i + l]) * _dt[1 + l] + 2. * dN(j, l) * _dt[l] - 2. * dN(j, 1 + l) * _dt[1 + l];

        dN(j, l) = dN(j, l) * (x[j] - t[l]) * _dt[l] - dN(j, 1 + l) * (x[j] - t[1 + i + l]) * _dt[1 + l] + N(j, l) * _dt[l] - N(j, 1 + l) * _dt[1 + l];

        N(j, l) = N(j, l) * (x[j] - t[l]) * _dt[l] - N(j, 1 + l) * (x[j] - t[1 + i + l]) * _dt[1 + l];
      }
    }
  }
  N.conservativeResize(N.rows(), N.cols() - 3);
  dN.conservativeResize(dN.rows(), dN.cols() - 3);
  d2N.conservativeResize(d2N.rows(), d2N.cols() - 3);

  return {N, dN, d2N};
}

FullTunneling PathDeformation::full_tunneling(std::vector<Eigen::VectorXd> path_pts) {
  FullTunneling ft;
  std::vector<double> phi_1d, dphi_1d;
  
  // Create SplinePath once and reuse
  std::unique_ptr<SplinePath> path_ptr;
  bool first_iteration = true;
  
  for (int num_iter = 1; num_iter <= path_maxiter; num_iter++) {
    LOG(debug) << "Starting tunneling step " << num_iter;
    
    if (first_iteration) {
      // First iteration: full construction
      path_ptr = std::make_unique<SplinePath>(P, T, path_pts, extend_to_minima, true);
      first_iteration = false;
    } else {
      // Subsequent iterations: fast update (no extension, no re-evaluation)
      path_ptr->update_path(path_pts, false);
    }
    
    PhaseTracer::Shooting tobj(*path_ptr, num_dims - 1);
    tobj.set_xtol(xtol);
    tobj.set_phitol(phitol);
    tobj.set_thin_cutoff(thin_cutoff);
    tobj.set_rmin(rmin);
    tobj.set_rmax(rmax);
    tobj.set_max_iter(max_iter);

    auto profile = tobj.findProfile(path_ptr->get_path_length(), 0.);
    if (profile.R.size() == 2) {
      if (profile.R.isApprox(tobj.profile_inf.R, tobj.get_xtol())) {
        ft.action = std::numeric_limits<double>::max();
        bounce_action = std::numeric_limits<double>::max();
        LOG(debug) << "Action is infinity.";
        return ft;
      } else if (profile.R.isApprox(tobj.profile_zero.R, tobj.get_xtol())) {
        ft.action = 0;
        bounce_action = 0;
        LOG(debug) << "Action is too small.";
        return ft;
      }
    }
    num_nodes = profile.Phi.size();
    tobj.evenlySpacedPhi(profile, &phi_1d, &dphi_1d, num_nodes, false);
    dphi_1d[0] = 0.;
    dphi_1d.back() = 0.;
    bounce_action = tobj.calAction(profile);
    //      std::cout << "action = " << std::setprecision(10) <<  bounce_action << std::endl;

    phi_node.resize(num_nodes);
    for (size_t ii = 0; ii < num_nodes; ii++) {
      phi_node[ii] = path_ptr->vecp(phi_1d[ii]);
    }

    bool converged;
    try {
      converged = deformPath(dphi_1d);
    } catch (...) {
      converged = false;
    }
    path_pts = phi_node;

    ft.saved_steps.push_back(phi_list);

    if (converged and num_steps < 2) {
      breakLoop = true;
      break;
    }
  }

  if (!breakLoop) {
    LOG(warning) << "Reached maxiter in full_tunneling. No convergence.";
  }

  // Calculate the ratio of max perpendicular force to max gradient.
  // Make sure that we go back a step and use the forces on the path, not the
  // most recently deformed path.
  // bool converged = deformPath(dphi_1d);
  path_pts = phi_node;

  std::vector<Eigen::VectorXd> F, dV;
  forces(F, dV);
  double F_max = max_norm(F);
  double dV_max = max_norm(dV);
  double fRatio = F_max / dV_max;

  // Final profile calculation - do full rebuild for accuracy
  path_ptr->update_path(path_pts, true);
  PhaseTracer::Shooting tobj(*path_ptr, num_dims - 1);
  tobj.set_xtol(xtol);
  tobj.set_phitol(phitol);
  tobj.set_thin_cutoff(thin_cutoff);
  tobj.set_rmin(rmin);
  tobj.set_rmax(rmax);
  tobj.set_max_iter(max_iter);

  auto profile = tobj.findProfile(path_ptr->get_path_length(), 0.);
  auto action = tobj.calAction(profile);
  bounce_action = action;

  std::vector<Eigen::VectorXd> phi_for_profile1D;
  for (size_t ii = 0; ii < profile.Phi.size(); ii++) {
    phi_for_profile1D.push_back(path_ptr->vecp(profile.Phi[ii]));
  }
  ft.phi_for_profile1D = phi_for_profile1D;

  ft.fRatio = fRatio;
  ft.phi = phi_node;
  ft.profile1D = profile;
  ft.action = action;

  LOG(debug) << "Tunneling step converged. ";
  return ft;
}

} // namespace PhaseTracer
