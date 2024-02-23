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

#ifndef PHASETRACER_PATH_DEFORMATION_HPP_INCLUDED
#define PHASETRACER_PATH_DEFORMATION_HPP_INCLUDED

#include <cmath>
#include <ostream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <boost/cstdint.hpp>
#include <nlopt.hpp>
//#include <boost/math/special_functions/bessel.hpp>
//#include <boost/numeric/odeint.hpp>
//#include <boost/math/quadrature/gauss_kronrod.hpp>
//#include <boost/math/tools/minima.hpp>
//#include <gsl/gsl_sf_bessel.h>


#include "logger.hpp"
#include "transition_finder.hpp"

namespace PhaseTracer {


//struct Profile1D {
//  Eigen::VectorXd R;
//  Eigen::VectorXd Phi;
//  Eigen::VectorXd dPhi;
//  double Rerr;
//};

/* Calculates dy/dx to fourth-order using finite differences. */
std::vector<Eigen::VectorXd> deriv14_const_dx(const std::vector<Eigen::VectorXd>& y, double dx = 1.) {
  const int ny = y.size();
  const double factor = 1.0 / (12.0*dx);
  std::vector<Eigen::VectorXd> dy(ny, Eigen::VectorXd::Zero(y[0].size()));
  for (int i = 2; i < ny - 2; ++i) {
      dy[i] = factor * (-y[i + 2] + 8.0 * y[i + 1] - 8.0 * y[i - 1] + y[i - 2]);
  }
  
  dy[0] = factor * (-25*y[0] + 48*y[1] - 36*y[2] + 16*y[3] - 3*y[4]);
  dy[1] = factor * (-3 *y[0] - 10*y[1] + 18*y[2] - 6 *y[3] +   y[4]);
  dy[ny-2] = factor * (3*y[ny-1] + 10*y[ny-2] - 18*y[ny-3] + 6*y[ny-4] - y[ny-5]);
  dy[ny-1] = factor * (25*y[ny-1] - 48*y[ny-2] + 36*y[ny-3] - 16*y[ny-4] + 3*y[ny-5]);
  return dy;
}


/* Fit a spline to a path in field space, and find the potential on that path */
class SplinePath: public PotentialForShooting {

public:
//  explicit SplinePath(TransitionFinder& tf_,
//                      std::vector<Eigen::VectorXd> pts_,
//                      bool extend_to_minima_ = true,
//                      bool reeval_distances_ = true) :
//    tf(tf_), pts(pts_), V_spline_samples(V_spline_samples_),
//    extend_to_minima(extend_to_minima_), reeval_distances(reeval_distances_) {
  explicit SplinePath(double a,
                      std::vector<Eigen::VectorXd> pts_,
                      bool extend_to_minima_ = true,
                      bool reeval_distances_ = true) :
    pts(pts_), npoints(pts_.size()),
    extend_to_minima(extend_to_minima_), reeval_distances(reeval_distances_) {
      // 1. Find derivs
      std::vector<Eigen::VectorXd> dpts = _pathDeriv(pts);
      std::cout << "SplinePath : " << std::endl;
      for (const auto& v : dpts) {
          std::cout << "  " << v.transpose() << std::endl;
      }
      
      // 2. Extend the path
      if (extend_to_minima){
        double xmin = find_loc_min_w_guess(pts[0], dpts[0]); // TODO: This may return global mim
        xmin = std::min(xmin, 0.0);
        int nx = static_cast<int>(std::ceil(std::abs(xmin)-.5)) + 1;
        double stepSize = xmin / nx;
        std::vector<Eigen::VectorXd> new_pts;
        for (int i = 0; i < nx; ++i) {
            double x = xmin + i * stepSize;
            Eigen::VectorXd pt_ext = pts[0] + x * dpts[0];
            new_pts.push_back(pt_ext);
        }
        new_pts.insert(new_pts.end(), pts.begin() + 1, pts.end());
        pts = new_pts;
        
//        std::cout << "pts= " << pts << std::endl;
        xmin = find_loc_min_w_guess(pts[npoints-1], dpts[npoints-1]);
//        std::cout << "xmin= " << xmin << std::endl;
        xmin = std::max(xmin, 0.0);
        nx = static_cast<int>(std::ceil(std::abs(xmin) - 0.5)) + 1;
        stepSize = xmin / nx;
        std::vector<Eigen::VectorXd> new_pts2;
        for (int i = 0; i < nx; ++i) {
          // TODO this is not checked
          double x = xmin + (nx-i-1) * stepSize;
          Eigen::VectorXd pt_ext = pts[npoints-1] + x * dpts[npoints-1];
          new_pts2.push_back(pt_ext);
        }
//        std::cout << "new_pts2= " << new_pts2 << std::endl;
//        std::cout << "pts= " << pts << std::endl;
        pts.pop_back();
//        std::cout << "pts= " << pts << std::endl;
        pts.insert(pts.end(), new_pts2.begin(), new_pts2.end());
//        std::cout << "pts= " << pts << std::endl;
      }
      // 3. Find knot positions and fit the spline
      std::cout << "dpts = " << dpts << std::endl;
      std::vector<double> squared_sums;
      squared_sums.reserve(dpts.size());
      for (const auto& vec : dpts)
      {
          squared_sums.push_back(vec.norm());
      }
//      std::cout << "squared_sums = " << squared_sums << std::endl;
      std::vector<double> pdist = cumulative_trapezoidal_integration(squared_sums);
      pdist.insert(pdist.begin(), 0.0);
//      std::cout << "pdist = ";
//      for (const auto& value : pdist)
//      {
//          std::cout << value << " ";
//      }
//      std::cout << std::endl;
      length = pdist[pdist.size()-1];
      
      get_path_tck(pdist);
//      alglib::real_1d_array x_arr;
//      alglib::real_1d_array y_arr;
//      x_arr.setlength(np);
//      y_arr.setlength(np);
//      _path_tck.clear();
//      for (size_t i = 0; i < np; i++) x_arr[i] = pdist[i];
//      for (size_t j = 0; j < nphi; j++){
//        alglib::spline1dinterpolant spline;
//        for (size_t i = 0; i < np; i++) y_arr[i] = pts[i](j);
//        alglib::spline1dbuildcubic(x_arr, y_arr, spline);
//        _path_tck.push_back(spline);
//      }
      
      // 4. Re-evaluate the distance to each point.
      if (reeval_distances){
        std::vector<double> pdist_(pdist.size());
        typedef double state_type;
        state_type x = 0.0;
        double t_start = 0.0;
        double dt = pdist[pdist.size()-1]*1E-4; // TODO this need be changed.
        using error_stepper_type =
           boost::numeric::odeint::runge_kutta_dopri5<state_type>;
        using controlled_stepper_type =
           boost::numeric::odeint::controlled_runge_kutta<error_stepper_type>;
        controlled_stepper_type stepper
           = make_controlled(0., pdist[pdist.size()-1]*1E-8, error_stepper_type());
        for (size_t i = 0; i < pdist.size(); ++i) // TODO this may be improved
        {
          boost::numeric::odeint::integrate_const(stepper,
            [this](const state_type& y, state_type& dydr, double r) {dpdx(y, dydr, r);},
            x, t_start, pdist[i], dt);
          pdist_[i] = x;
        }
        pdist = pdist_;
        length = pdist[pdist.size()-1];
        get_path_tck(pdist);
        
        
        std::cout << "pdist = ";
        for (const auto& value : pdist)
        {
            std::cout << value << " ";
        }
        std::cout << std::endl;
        
        std::cout << "x = 0.5 --->";
        for (const auto& value : vecp(0.5))
        {
            std::cout << value << " ";
        }
        std::cout << std::endl;
        std::cout << "V(0.5) = " << p_V(vecp(0.5)) << std::endl;
      }

      // Make the potential spline.
      alglib::real_1d_array p_arr;
      alglib::real_1d_array v_arr;
      p_arr.setlength(V_spline_samples);
      v_arr.setlength(V_spline_samples);
      for (size_t i = 0; i < V_spline_samples; ++i)
      {
        // extend 20% beyond
        p_arr[i] = -0.2*length + i * (1.4*length) / (V_spline_samples - 1);
        v_arr[i] = p_V(vecp(p_arr[i]));
      }
      alglib::spline1dbuildcubic(p_arr, v_arr, _V_tck);
    }
  
  virtual ~SplinePath() = default;
  
  /* Calculates to 4th order if len(phi) >= 5, otherwise 1st/2nd order. */
  std::vector<Eigen::VectorXd> _pathDeriv(const std::vector<Eigen::VectorXd> phi) { // rename phi
    std::vector<Eigen::VectorXd> dphi(npoints);
    if (npoints < 2){
      throw std::runtime_error("The number of points that describe the path must be larger than 1.");
    } else if (npoints == 2){
      dphi[0] = dphi[1] = phi[1] - phi[0];
    } else if (npoints < 5){
      // 1st/2nd order calculation
      dphi[0] = -1.5 * phi[0] + 2.0 * phi[1] - 0.5 * phi[2];
      dphi[npoints - 1] = 1.5 * phi[npoints - 1] - 2.0 * phi[npoints - 2] + 0.5 * phi[npoints - 3];
      for (int i = 1; i < npoints - 1; ++i) {
          dphi[i] = 0.5 * (phi[i + 1] - phi[i - 1]);
      }
    } else{
      dphi = deriv14_const_dx(phi);
    }
    return dphi;
  }
  

  double find_loc_min_w_guess(Eigen::VectorXd p0, Eigen::VectorXd dp0, double guess = 0.){

    std::shared_ptr<std::function<double(const std::vector<double>&)>> V_lin = std::make_shared<std::function<double(const std::vector<double>&)>>([this, p0, dp0](const std::vector<double>& x) {
      return this->p_V(p0 + x[0] * dp0);
    });
    
    auto V_lin_func = [](const std::vector<double>& x, std::vector<double>& grad, void* data) -> double {
      auto v_lin = reinterpret_cast<std::shared_ptr<std::function<double(const std::vector<double>&)>>*>(data);
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
  
  template<typename T>
  std::vector<T> cumulative_trapezoidal_integration(const std::vector<T>& values)
  {
      std::vector<T> integration_result;
      integration_result.reserve(values.size());
      T sum = T(0);
      for (size_t i = 0; i < values.size() - 1; ++i)
      {
          sum += (values[i] + values[i + 1]) / 2.0;
          integration_result.push_back(sum);
      }
      return integration_result;
  }
  
  void get_path_tck(std::vector<double> pdist){
    alglib::real_1d_array x_arr;
    alglib::real_1d_array y_arr;
    x_arr.setlength(npoints);
    y_arr.setlength(npoints);
    _path_tck.clear();
    for (size_t i = 0; i < npoints; i++) x_arr[i] = pdist[i];
    for (size_t j = 0; j < nphi; j++){
      alglib::spline1dinterpolant spline;
      for (size_t i = 0; i < npoints; i++) y_arr[i] = pts[i](j);
      alglib::spline1dbuildcubic(x_arr, y_arr, spline);
      _path_tck.push_back(spline);
    }
  }
  
  void dpdx(const double& x, double& dpdx_, double r){
    Eigen::VectorXd vecx(nphi);
    for (int i=0; i<nphi; i++){
      double s, ds, d2s;
      alglib::spline1ddiff(_path_tck[i], x, s, ds, d2s);
      vecx[i] = ds;
    }
    dpdx_ = vecx.norm();
  }
  
  Eigen::VectorXd vecp(double x){
    Eigen::VectorXd vec(nphi);
    for (int i=0; i<nphi; i++){
      vec[i] = alglib::spline1dcalc(_path_tck[i], x);
    }
    return vec;
  }
  
  double V(double x) const override {
    return alglib::spline1dcalc(_V_tck, x);
  }
  
  double dV(double x) const override {
    double s, ds, d2s;
    alglib::spline1ddiff(_V_tck, x, s, ds, d2s);
    return ds;
  }
  
  double d2V(double x) const override {
    double s, ds, d2s;
    alglib::spline1ddiff(_V_tck, x, s, ds, d2s);
    return d2s;
  }
  
  
  double p_V(Eigen::VectorXd phi){
    double x=phi[0];
    double y=phi[0];
    double r1 = x*x+c*y*y;
    double r2 = c*pow(x-1.,2) + pow(y-1.,2);
    double r3 = fx*(0.25*pow(x,4) - pow(x,3)/3.);
    r3 += fy*(0.25*pow(y,4) - pow(y,3)/3.);
    return r1*r2 + r3;
  }
  Eigen::VectorXd p_dV(Eigen::VectorXd phi){  // TODO this must be replcaed, // this is not used
    double x=phi[0];
    double y=phi[0];
    double r1 = x*x+c*y*y;
    double r2 = c*pow(x-1.,2) + pow(y-1.,2);
    double dr1dx = 2*x;
    double dr1dy = 2*c*y;
    double dr2dx = 2*c*(x-1.);
    double dr2dy = 2*(y-1.);
    double dVdx = r1*dr2dx + dr1dx*r2 + fx*x*x*(x-1.);
    double dVdy = r1*dr2dy + dr1dy*r2 + fy*y*y*(y-1.);
    Eigen::VectorXd rval(2);
    rval << dVdx, dVdy;
    return rval;
  }
  
  double get_path_length(){return length;}
  
private:
  
//  TransitionFinder& tf;
  std::vector<Eigen::VectorXd> pts;
  int V_spline_samples = 140; // larger than 5
  bool extend_to_minima;
  bool reeval_distances;
  
  double length;
  
  size_t nphi=2;
  size_t npoints;
  
  std::vector<alglib::spline1dinterpolant> _path_tck;
  alglib::spline1dinterpolant _V_tck;
  
  double c=5.;
  double fx=0.;
  double fy=2.;
  
};


class PathDeformation {
public:
//  explicit PathDeformation(TransitionFinder& tf_) :
//    tf(tf_) {
//      // TODO
//    }
//  TransitionFinder& tf;
  
  explicit PathDeformation(double a)
    {
      double b=a;
      // TODO
    }
  
  virtual ~PathDeformation() = default;
  
  
  double c=5.;
  double fx=0.;
  double fy=2.;
  
  double p_V(Eigen::VectorXd phi){
    double x=phi[0];
    double y=phi[0];
    double r1 = x*x+c*y*y;
    double r2 = c*pow(x-1.,2) + pow(y-1.,2);
    double r3 = fx*(0.25*pow(x,4) - pow(x,3)/3.);
    r3 += fy*(0.25*pow(y,4) - pow(y,3)/3.);
    return r1*r2 + r3;
  }
  
  Eigen::VectorXd p_dV(Eigen::VectorXd phi){  // TODO this must be replcaed, // this is not used
    double x=phi[0];
    double y=phi[0];
    double r1 = x*x+c*y*y;
    double r2 = c*pow(x-1.,2) + pow(y-1.,2);
    double dr1dx = 2*x;
    double dr1dy = 2*c*y;
    double dr2dx = 2*c*(x-1.);
    double dr2dy = 2*(y-1.);
    double dVdx = r1*dr2dx + dr1dx*r2 + fx*x*x*(x-1.);
    double dVdy = r1*dr2dy + dr1dy*r2 + fy*y*y*(y-1.);
    Eigen::VectorXd rval(2);
    rval << dVdx, dVdy;
    return rval;
  }
  
  void set_deformation(std::vector<Eigen::VectorXd> phi_, std::vector<double> dphidr){
    phi = phi_; // TODO
    
    size_t nb=10;
    size_t kb=3;
    double v2min=0.0;
    bool fix_start=false;
    bool fix_end=false;
    bool save_all_steps=false;
    
    
    // convert phi to a set of path lengths
    std::vector<double> dL;
    for (int i = 0; i < phi.size() - 1; i++) {
        dL.push_back((phi[i + 1] - phi[i]).norm());
    }
    std::vector<double> _t;
    _t.push_back(0);
    for (int i = 0; i < dL.size(); i++) {
        _t.push_back(_t[i] + dL[i]);
    }
    double totalLength = _t.back();
    for (int i = 0; i < _t.size(); i++) {
        _t[i] /= totalLength;
    }
    _t[0] = 1e-50;
    
    // create the starting spline
    std::vector<double> t0;
    t0.push_back(0.0);
    t0.push_back(0.0);
    for (int i = 0; i < nb; ++i) {
        t0.push_back( i / (nb - 1.) );
    }
    t0.push_back(1.0);
    t0.push_back(1.0);
    
//    for (int i = 0; i < t0.size(); ++i) {
//         std::cout << "t0[" << i << "]: " << t0[i] << std::endl;
//     }
    auto result = Nbspld2(t0, _t, kb);
    X = std::get<0>(result);
    dX = std::get<1>(result);
    d2X = std::get<2>(result);
    
    Eigen::VectorXd phi0 = phi[0];
    Eigen::VectorXd phi1 = phi[phi.size()-1];
    beta.clear();
    for (int j = 0; j < phi0.size(); ++j) {
      Eigen::VectorXd phi_delta(phi.size());
      for (int i = 0; i < phi.size(); ++i) {
        Eigen::VectorXd phii = phi[i];
        phi_delta[i] = phii[j] - (phi0[j] + (phi1[j] - phi0[j])* _t[i]);
      }
      beta.push_back(X.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(phi_delta));
//      std::cout << " phi_delta = " << phi_delta << std::endl;
//      std::cout << " beta[j] = " << beta[j] << std::endl;
    }
    
    
    
    double max_v2 = 0;
    for (int i = 0; i < phi.size(); ++i) {
      max_v2 = std::max(max_v2, p_dV(phi[i]).norm());
    }
    v2min *=  max_v2 * totalLength/nb;
    
    v2.clear();
    for (int i = 0; i < dphidr.size(); ++i) {
      v2.push_back(std::max(v2min, dphidr[i]*dphidr[i]));
    }
    
    // Deform the path
    double startstep=2e-3;
    double fRatioConv=.02;
    double converge_0=5.;
    double fRatioIncrease=5.;
    size_t maxiter=500;
    
//    minfRatio = np.inf
//    minfRatio_index = 0
//    minfRatio_beta = None
//    minfRatio_phi = None
    double stepsize = startstep;
//    deformation_converged = False
    while(true){
      num_steps++;
      step(stepsize);
      if (num_steps >= maxiter){
        LOG(debug) << "Maximum number of deformation iterations reached.";
        break;
      }
    }
    

  }
  
  size_t nphi=2;
  size_t num_steps=0;
  Eigen::MatrixXd X;
  Eigen::MatrixXd dX;
  Eigen::MatrixXd d2X;
  std::vector<Eigen::VectorXd> beta;
  std::vector<Eigen::VectorXd> phi;
  std::vector<double> v2;
  
  void step(size_t lastStep ){
    forces();
  }
  
  // Calculate the normal force and potential gradient on the path
  void forces(){
    // X, dX, d2X
    
    
    std::cout << "Matrix N.rows:" << dX.rows() << std::endl;
    std::cout << "Matrix N.cols:" << dX.cols() << std::endl;
    
    std::vector<Eigen::VectorXd> dphi(dX.rows());
    std::vector<Eigen::VectorXd> d2phi(dX.rows());
    std::vector<double> dphi_sq(dX.rows());
    std::vector<Eigen::VectorXd> dphids(dX.rows());
    std::vector<Eigen::VectorXd> d2phids2(dX.rows());
    
    std::vector<Eigen::VectorXd> dV(dX.rows());
    std::vector<Eigen::VectorXd> F_norm(dX.rows());
    for (int ii=0; ii<dX.rows(); ++ii){
      Eigen::VectorXd dphi_kk(nphi);
      Eigen::VectorXd d2phi_kk(nphi);
      for(int jj=0; jj<nphi; ++jj){
        dphi_kk[jj] = 0.;
        d2phi_kk[jj] = 0.;
        for (int kk=0; kk<dX.cols(); ++kk){
          dphi_kk[jj] += beta[jj][kk] * dX(ii,kk);
          d2phi_kk[jj] += beta[jj][kk] * d2X(ii,kk);
        }
        dphi_kk[jj] += phi.back()[jj] - phi[0][jj];
      }
      
      dphi[ii] = dphi_kk;
      d2phi[ii] = d2phi_kk;
      double dphi_ = dphi_kk.norm();
      dphi_sq[ii] = dphi_*dphi_;
      dphids[ii] = dphi_kk / dphi_;
      double sum_dphi_d2phi = 1;
      d2phids2[ii] = (d2phi_kk - dphi_kk*(dphi_kk.dot(d2phi_kk))/dphi_sq[ii])/dphi_sq[ii];
      
      dV[ii] = p_dV(phi[ii]);
      Eigen::VectorXd dV_perp = dV[ii] - dV[ii].dot(dphids[ii]) * dphids[ii];
      F_norm[ii] = d2phids2[ii] * v2[ii] - dV_perp;
//      std::cout << "phi[ii]:" << phi[ii] << std::endl;
      std::cout << "F_norm[ii]:" << F_norm[ii] << std::endl;
    }
    std::exit(0);
  }
  
  
  
  
  std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> Nbspld2(std::vector<double> t, std::vector<double> x, int k = 3) {
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
//    std::cout << "Matrix N:" << std::endl << N << std::endl;
//      std::cout << "Matrix N.rows:" << N.rows() << std::endl;
//    std::cout << "Matrix N.cols:" << N.cols() << std::endl;
    
      for (int i = 1; i <= kmax; ++i) {
        Eigen::VectorXd dt(t.size() - i);
        Eigen::VectorXd _dt(t.size() - i);
          
        for (int j = 0; j < t.size() - i; ++j) {
          dt(j) = t[j + i] - t[j];
          _dt(j) = (dt(j) != 0) ? 1.0 / dt(j) : 0.0;
        }
        for (int j = 0; j < size_x; ++j) {
          for (int l = 0; l < dt.size() - 1; ++l) {
            
            d2N(j,l) = d2N(j,l) * (x[j] - t[l]) * _dt[l]
                      -d2N(j,1+l) * (x[j] - t[1+i+l])* _dt[1+l]
                      +2.*dN(j,l) * _dt[l] - 2.*dN(j,1+l) * _dt[1+l];
            
            dN(j,l) = dN(j,l) * (x[j] - t[l]) * _dt[l]
                      -dN(j,1+l) * (x[j] - t[1+i+l])* _dt[1+l]
                      +N(j,l) * _dt[l] - N(j,1+l) * _dt[1+l];
            
            N(j,l) = N(j,l) * (x[j] - t[l]) * _dt[l]
                    -N(j,1+l) * (x[j] - t[1+i+l]) * _dt[1+l];
          }
        }
      }
    N.conservativeResize(N.rows(), N.cols() - 3);
    dN.conservativeResize(dN.rows(), dN.cols() - 3);
    d2N.conservativeResize(d2N.rows(), d2N.cols() - 3);
    
//    std::cout << "Matrix d2N:" << std::endl;
//    for (int ii = 0; ii < 3; ++ii) {
//        for (int j = 0; j < N.cols(); ++j) {
//            std::cout << d2N(ii, j) << " ";
//        }
//        std::cout << std::endl;
//    }
//    for (int ii = N.rows() - 3; ii < N.rows(); ++ii) {
//        for (int j = 0; j < N.cols(); ++j) {
//            std::cout << d2N(ii, j) << " ";
//        }
//        std::cout << std::endl;
//    }
    
    return {N, dN, d2N};
  }
  
  
  /* Calculate the instanton solution in multiple field dimension */
  double fullTunneling(double delta_phi){
    std::vector<Eigen::VectorXd> path_pts;
    path_pts.push_back(Eigen::VectorXd(2));
    path_pts.push_back(Eigen::VectorXd(2));
    path_pts[0] << 1, 1;
    path_pts[1] << 0, 0;
    
    int maxiter = 1;
    int V_spline_samples = 100;
    bool extend_to_minima = true;
    
    for (int num_iter = 1; num_iter <= maxiter; num_iter++) {
      LOG(debug) <<  "Starting tunneling step " << num_iter;
      SplinePath path(0, path_pts);
      
      std::cout << "path.V(0) = " << path.V(0) << std::endl;
      std::cout << "path.V(0.5) = " << path.V(0.5) << std::endl;
      std::cout << "path.V(1) = " << path.V(1) << std::endl;
      
      std::cout << "path.dV(0.5) = " << path.dV(0.5) << std::endl;
      std::cout << "path.d2V(0.5) = " << path.d2V(0.5) << std::endl;
      
      PhaseTracer::Shooting tobj(path);
      
      std::cout << "path.get_path_length() = " << path.get_path_length() << std::endl;
      
      auto profile = tobj.findProfile(path.get_path_length(),0.);
      
      std::vector<double> phi, dphi;
      tobj.evenlySpacedPhi(profile, &phi, &dphi, profile.Phi.size(), false);
      dphi[0]=0.;
      dphi[dphi.size()-1]=0.;
      auto action = tobj.calAction(profile);
      std::cout << "action = " << action << std::endl;
      std::vector<Eigen::VectorXd> pts(phi.size());
      for (size_t ii=0; ii<phi.size(); ii++ ){
        pts[ii] = path.vecp(phi[ii]);
      }
      // TODO add callback
      
      set_deformation(pts, dphi);
      
    }
    
    
    
    return 0;
  }
  
private:

  double phi_absMin;
  double phi_metaMin;
  
};

}  // namespace PhaseTracer

#endif
