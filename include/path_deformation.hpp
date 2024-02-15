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
      std::cout << "path.V(1) = " << path.V(
                                            1) << std::endl;
      
      PhaseTracer::Shooting s(path);
      
      auto profile = s.findProfile(0,path.get_path_length());
      
      auto action = s.calAction(profile);
      std::cout << "action = " << action << std::endl;
      
      
    }
    
    
    
    return 0;
  }
  
private:

  double phi_absMin;
  double phi_metaMin;
  
};

}  // namespace PhaseTracer

#endif
