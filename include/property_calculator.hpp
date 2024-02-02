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

#ifndef PHASETRACER_PROPERTYCALCULATOR_HPP_INCLUDED
#define PHASETRACER_PROPERTYCALCULATOR_HPP_INCLUDED

#include <cmath>
#include <ostream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <boost/cstdint.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/numeric/odeint.hpp>
#include <gsl/gsl_sf_bessel.h>
#include "nlopt.hpp"


#include "logger.hpp"
#include "transition_finder.hpp"

namespace PhaseTracer {


class CubicInterpFunction {
public:
  CubicInterpFunction(double y0, double dy0, double y1, double dy1, double c = 0)
    : y0(y0), dy0(dy0), y1(y1), dy1(dy1), c(c) {
      y_1 = y0 + dy0 / 3.0;
      y_2 = y1 - dy1 / 3.0;
      y_3 = y1;
    }
  double operator()(double t) const {
    double mt = 1. - t;
    return y0 * pow(mt, 3) + 3.0 * y_1 * mt * mt * t + 3.0 * y_2 * mt * t * t + y_3 * pow(t, 3) - c;
  }
private:
  double y0;
  double dy0;
  double y1;
  double dy1;
  double y_1;
  double y_2;
  double y_3;
  double c;
};

class PropertyCalculator {
public:
  explicit PropertyCalculator(TransitionFinder& tf_) :
    tf(tf_) {
    std::cout << "PropertyCalculator starts " << std::endl;
    phi_eps = phi_eps * fabs(phi_absMin - phi_metaMin); // TODO: this should not be forgot
  }
  virtual ~PropertyCalculator() = default;
  

  TransitionFinder& tf;

//  double get_action(const Eigen::VectorXd& vacuum_1, const Eigen::VectorXd& vacuum_2, double T) const{
//
//  }
  
  double V(double phi){
    return 0.25*pow(phi,4) - 0.49*pow(phi,3) + 0.235 * pow(phi,2);
  }
  double dV(double phi){
    return phi*(phi-.47)*(phi-1);
  }
  
  double d2V(double phi){
    return (phi-.47)*(phi-1) + phi*(phi-1) + phi*(phi-.47);
  }

  
  double phi_bar= 0.837171431429;
  double phi_absMin= 1.0;
  double phi_metaMin= 0.0;
  /*The approximate radial scale of the instanton*/
  double rscale = 1.6677089434383516;
 
  
  /*Calculates `dV/dphi` at ``phi = phi_absMin + delta_phi``.*/
  double dV_from_absMin(double delta_phi){
    double phi = phi_absMin + delta_phi;
    double dV_ = dV(phi);
    if (phi_eps > 0){
      double dV_eps = d2V(phi) * delta_phi;
      double blend_factor = exp(-pow((delta_phi/phi_eps),2));
      dV_ = dV_eps*blend_factor + dV_*(1-blend_factor);
    }
      
    return dV_;
  }
  
  /*Find `phi(r)` given `phi(r=0)`, assuming a quadratic potential.*/
  void exactSolution(double r, double phi0, double dV_, double d2V_,
                     double* phi_r, double* dphi_r){
    double beta = sqrt(fabs(d2V_));
    double beta_r = beta*r;
    double nu = 0.5*(alpha-1.);
    double phi = 0., dphi = 0.;
    if (beta_r < 1e-2){
      int s = d2V_>0 ? +1 : -1;
      for (int k=1; k<4; k++){
        double temp = pow(0.5*beta_r, 2.*k-2.) * pow(s,k) / (std::tgamma(k+1.)*std::tgamma(k+1.+nu));
        phi += temp;
        dphi += temp * (2.*k);
      }
      phi *= 0.25 * std::tgamma(nu+1) * pow(r,2) * dV_ * s;
      dphi *= 0.25 * std::tgamma(nu+1) * r * dV_ * s;
      phi += phi0;
      std::cout << "here1" << " phi = " << phi << std::endl;
    }else if (d2V_>0){
      // TODO: ignore the warnings if there is
      phi = (std::tgamma(nu+1)*pow(0.5*beta_r,-nu) * boost::math::cyl_bessel_i(nu, beta_r)-1) * dV_/d2V_;
      dphi = -nu*(pow(0.5*beta_r,-nu) / r) * boost::math::cyl_bessel_i(nu, beta_r);
      dphi += pow(0.5*beta_r,-nu) * 0.5*beta
              * (boost::math::cyl_bessel_i(nu-1, beta_r)+boost::math::cyl_bessel_i(nu+1, beta_r));
      dphi *= std::tgamma(nu+1) * dV_/d2V_;
      phi += phi0;
      std::cout << "here2" << " phi = " << phi << std::endl;
    }else{
      phi = (std::tgamma(nu+1)*pow(0.5*beta_r,-nu) * boost::math::cyl_bessel_j(nu, beta_r)-1) * dV_/d2V_;
      dphi = -nu*(pow(0.5*beta_r,-nu) / r) * boost::math::cyl_bessel_j(nu, beta_r);
      dphi += pow(0.5*beta_r,-nu) * 0.5*beta
              * (boost::math::cyl_bessel_j(nu-1, beta_r)-boost::math::cyl_bessel_j(nu+1, beta_r));
      dphi *= std::tgamma(nu+1) * dV_/d2V_;
      phi += phi0;
      std::cout << "here3" << std::endl;
    }
    *phi_r = phi;
    *dphi_r = dphi;
  }
  
  /*Finds the initial conditions for integration.*/
  void initialConditions(double delta_phi0, double rmin, double delta_phi_cutoff,
                         double* r0, double* phi_r0, double* dphi_r0){
    // Finds the value `r0` such that `phi(r0) = phi_cutoff`.
    // If there is no such value, it returns the intial conditions at `rmin`.
    double phi0 = phi_absMin + delta_phi0;
    double dV_ = dV_from_absMin(delta_phi0);
    double d2V_ = d2V(phi0);
    
    *r0 = rmin;
    
    std::cout << "rmin, phi0, dV_, d2V_ = " << rmin << "    " << phi0
    << "    " << dV_ << "    " << d2V_ << std::endl;
    
    const auto deltaPhiDiff = [this, phi0, dV_, d2V_, phi_r0, dphi_r0, delta_phi_cutoff](double rtry) {
      this->exactSolution(rtry, phi0, dV_, d2V_, phi_r0, dphi_r0);
      return fabs(*phi_r0 - phi_absMin) - fabs(delta_phi_cutoff);
    };
    if ( deltaPhiDiff(rmin) > 0) return;
    if ( (*dphi_r0)*delta_phi0 < 0 ) return; // wrong direction
 
    std::cout << "phi_r0, dphi_r0 = " << *phi_r0 << "    " << *dphi_r0 << std::endl;
    
    // Find the smallest r0 such that delta_phi_r0 > delta_phi_cutoff
    std::cout << "This part has not been checcked!" << std::endl;
    double r = rmin;
    double rlast = r;
    while(std::isfinite(r)){
      rlast = r;
      r *= 10.;
      if ( deltaPhiDiff(r) > 0) break;
    }
    
    // Find phi - phi_absMin = delta_phi_cutoff exactly
    
    const auto result = boost::math::tools::bisect(deltaPhiDiff, rlast, r, boost::math::tools::eps_tolerance<double>(), max_iter);
    *r0 = (result.first + result.second) * 0.5;
    exactSolution(*r0, phi0, dV_, d2V_, phi_r0, dphi_r0);
    return;
  }
  
  void equationOfMotion(const std::vector<double>& y, std::vector<double>& dydr, const double r){
    dydr[0] = y[1];
    dydr[1] = dV(y[0]) - alpha * y[1] /r ;
  }
  
  std::vector<double> dY(const std::vector<double> y, const double r){
    std::vector<double> dydr(2);
    equationOfMotion(y,dydr,r);
    return dydr;
  }
  
  int integrateProfile(double r0, std::vector<double> y0,
                        double* rf, std::vector<double>* yf, std::vector<double> epsabs, std::vector<double> epsfrac){
    
    int convergence_type;
    int ysign = y0[0] > phi_metaMin ? 1 : -1;
    double dr0 = 0.00016677089434383516;
    double rmax = 16717.451579307573;
    double drmin = 1.6677089434383517e-06;

    using state_type = std::vector<double>;
    using error_stepper_type =
       boost::numeric::odeint::runge_kutta_dopri5<state_type>;
    using controlled_stepper_type =
       boost::numeric::odeint::controlled_runge_kutta<error_stepper_type>;
    state_type y1 = y0;
    double r1 = r0;
    double dr = dr0;
    state_type dydr1(2), dydr0(2);
    double eps_abs = 0.0001; // TODO
    double eps_rel = 1.0e-04; // TODO
    controlled_stepper_type stepper
       = make_controlled(eps_rel, eps_abs,
                         error_stepper_type());
    while (true) {
      y0 = y1;
      r0 = r1;
//      std::cout << "-- r = " << r1 << ", phi = " << y1[0] << ", dphi/dr = " << y1[1] << std::endl;
      while ( stepper.try_step(
            [this](const state_type& y, state_type& dydr, double r) {equationOfMotion(y, dydr, r);},
                               y1, r1, dr) ) {};
//      std::cout << "-- r = " << r0 << ", phi = " << y0[0] << ", dphi/dr = " << y0[1] << std::endl;
      std::cout << "== r = " << r1 << ", phi = " << y1[0] << ", dphi/dr = " << y1[1] << std::endl;
      if (r1>rmax){
        throw std::runtime_error("r > rmax");
      } else if (dr < drmin){
        throw std::runtime_error("dr < drmin");
      } else if (fabs(y1[0]-phi_metaMin)< 3.*epsabs[0] and fabs(y1[1])< 3.*epsabs[1] ){
        convergence_type = 0; // converged
        *rf = r1;
        *yf = y1;
        break;
      } else if ( (y1[0]-phi_metaMin)*ysign < 0 or y1[1]*ysign > 0 ){
        dydr0 = dY(y0, r0);
        dydr1 = dY(y1, r1);
        CubicInterpFunction CF0(y0[0], dr*dydr0[0], y1[0], dr*dydr1[0], phi_metaMin);
        CubicInterpFunction CF1(y0[1], dr*dydr0[1], y1[1], dr*dydr1[1]);
        double x;
        if (y1[1]*ysign > 0){
          // Extrapolate to where dphi(r) = 0
          const auto result = boost::math::tools::bisect(CF1, 0., 1., boost::math::tools::eps_tolerance<double>(), max_iter);
          x = (result.first + result.second) * 0.5;
          convergence_type = -1; // "undershoot"
        } else {
          // Extrapolate to where phi(r) = phi_metaMin
          const auto result = boost::math::tools::bisect(CF0, 0., 1., boost::math::tools::eps_tolerance<double>(), max_iter);
          x = (result.first + result.second) * 0.5;
          convergence_type = 1; // "overshoot"
        }
        std::cout << "x = " << x << std::endl;
        *rf = r0+dr*x;
        *yf = {phi_metaMin + CF0(x), CF1(x)} ;
        break;
      }

    }
    std::cout << "convergence_type = " << convergence_type << std::endl;
    std::cout << "rf = " << *rf << ", yf = " << *yf << std::endl;
    
    if (fabs((*yf)[0]-phi_metaMin)< 3*epsabs[0] and fabs((*yf)[1]) < 3*epsabs[1] ) {
      convergence_type = 0;
    }
    
    return convergence_type;
  }
  
  struct Profile1D {
    Eigen::VectorXd R;
    Eigen::VectorXd Phi;
    Eigen::VectorXd dPhi;
    double Rerr;
  };
  
  Profile1D integrateAndSaveProfile(Eigen::VectorXd R, std::vector<double> y0,
                        double dr0, std::vector<double> epsabs, std::vector<double> epsfrac, double drmin){
    Profile1D profile;
    int N = R.size();
    double r0 = R[0];
    Eigen::VectorXd Phi(N);
    Eigen::VectorXd dPhi(N);
    Phi(0) = y0[0];
    dPhi(0) = y0[1];

    std::vector<double> y1 = y0;
    double r1 = r0;
    double dr = dr0;
    double Rerr = NAN;
    
    using state_type = std::vector<double>;
    using error_stepper_type =
       boost::numeric::odeint::runge_kutta_dopri5<state_type>;
    using controlled_stepper_type =
       boost::numeric::odeint::controlled_runge_kutta<error_stepper_type>;
    double eps_abs = 0.0001; // TODO
    double eps_rel = 1.0e-04; // TODO
    controlled_stepper_type stepper
       = make_controlled(eps_rel, eps_abs,
                         error_stepper_type());
    int ii = 1;
    while (ii<N){
      y0 = y1;
      r0 = r1;
//      std::cout << "-- r = " << r1 << ", phi = " << y1[0] << ", dphi/dr = " << y1[1] << std::endl;
      while ( stepper.try_step(
            [this](const state_type& y, state_type& dydr, double r) {equationOfMotion(y, dydr, r);},
                               y1, r1, dr) ) {};
      if (dr < drmin) {
        r1 = r0 + drmin;
        y1[0] = y0[0] + (y1[0]-y0[0])*drmin/dr;
        y1[1] = y0[1] + (y1[1]-y0[1])*drmin/dr;
        dr = drmin;
        if ( not std::isnan(Rerr) ) Rerr = r1; // TODO
      }

      if (r0 < R[ii] and R[ii] <= r1) {
        std::vector<double> dydr0 = dY(y0, r0);
        std::vector<double> dydr1 = dY(y1, r1);
        CubicInterpFunction CF0(y0[0], dr*dydr0[0], y1[0], dr*dydr1[0]);
        CubicInterpFunction CF1(y0[1], dr*dydr0[1], y1[1], dr*dydr1[1]);
        while ( ii < N and r0 < R[ii] and R[ii] <= r1 ) {
          double x = (R[ii]-r0)/dr;
          Phi[ii] = CF0(x);
          dPhi[ii] = CF1(x);
          ii ++;
        }
      }
    }
    profile.R = R;
    profile.Phi = Phi;
    profile.dPhi = dPhi;
    profile.Rerr = Rerr;
    
    //    throw std::runtime_error("test");
    return profile;
  }
  
  Profile1D findProfile(double xguess=NAN, int max_interior_pts=0){
    
    // Set r parameters
    rmin *= rscale;
    double dr0 = rmin;
    double drmin = 0.01*rmin;
    rmax *= rscale;
    // Set the phi parameters
    double delta_phi = phi_metaMin - phi_absMin;
    std::vector<double> epsabs = {fabs(delta_phi*phitol), fabs(delta_phi*phitol/rscale)};
    std::vector<double> epsfrac = {phitol, phitol};
    double delta_phi_cutoff = thinCutoff * delta_phi;
    // Set x parameters
    double xmin = xtol*10;
    double xmax = INF;
    double x;
    if (std::isnan(xguess)){ // TODO
      x =  -log(fabs((phi_bar-phi_absMin)/(phi_metaMin-phi_absMin)));
    } else {
      x = xguess;
    }
    double xincrease = 5.0;
    
    double r0, rf=NAN;
    std::vector<double> y0(2);
    double delta_phi0;
    while (true){
      delta_phi0 = exp(-x)*delta_phi;
      double r0_, phi0, dphi0;
      initialConditions(delta_phi0, rmin, delta_phi_cutoff, &r0_, &phi0, &dphi0);
      std::cout << "phi_r0 = " << phi0 << std::endl;
      
      if (not std::isfinite(r0_) or not std::isfinite(x)){
        if (std::isnan(rf)) {
          LOG(fatal) << "Failed to retrieve initial conditions on the first try.";
          std::terminate();
        }
        break;
      }
      
      r0 = r0_;
      y0[0] = phi0;
      y0[1] = dphi0;
      std::vector<double> yf(2);
      int shooting_type = integrateProfile(r0, y0, &rf, &yf, epsabs, epsfrac);

      if (shooting_type == 0){
        std::cout << "Converged" << std::endl;
        break;
      } else if (shooting_type == -1){ // undershoot
        xmin = x;
        x = std::isinf(xmax) ? x*xincrease : 0.5*(xmin+xmax);
      } else if (shooting_type == 1){ // overshoot
        xmax = x;
        x = .5*(xmin+xmax);
      }
      
      if ((xmax-xmin) < xtol){
        std::cout << "Reached xtol" << std::endl;
        break;
      }
    
    }
    
    // Getting the points along the path
    Eigen::VectorXd R = Eigen::VectorXd::LinSpaced(npoints, r0, rf);
    Profile1D profile_final = integrateAndSaveProfile(R, y0, dr0, epsfrac, epsabs, drmin);
    
    
    
    // Make points interior to the bubble
    Profile1D profile_int;
    
    if (max_interior_pts <= 0) max_interior_pts = int(R.size()/2);
    if (max_interior_pts>0){
      double dx0 = R[1]-R[0];
      
      int n;
      if (R[0] / dx0 <= max_interior_pts){
        n = std::ceil(R[0] / dx0);
        profile_int.R = Eigen::VectorXd::LinSpaced(n, 0.0, R[0]).head(n-1);
      } else {
        n = max_interior_pts;
        double a = (R[0]/dx0 - n) * 2./(n*(n+1.));
        profile_int.R.resize(n);
        for (int i = 0; i < n; i++) {
            int N = n - i;
          profile_int.R[i] = R[0] - dx0 * (N + 0.5 * a * N * (N + 1.));
        }
        profile_int.R[0] = 0.0;
      }
      int int_size = profile_int.R.size();
      profile_int.Phi.resize(int_size);
      profile_int.dPhi.resize(int_size);
      profile_int.Phi[0] = phi_absMin + delta_phi0;
      profile_int.dPhi[0] = 0.0;
      
      double dV_ = dV_from_absMin(delta_phi0);
      double d2V_ = d2V(profile_int.Phi[0]);
      
      for (int i=1; i < int_size; i++ ){
        exactSolution(profile_int.R[i], profile_int.Phi[0], dV_, d2V_, profile_int.Phi.data() + i, profile_int.dPhi.data() + i);
      }
    }
    Profile1D profile;
    int all_size = profile_int.R.size() + profile_final.R.size();
    profile.R.resize(all_size);
    profile.Phi.resize(all_size);
    profile.dPhi.resize(all_size);
    profile.R << profile_int.R, profile_final.R;
    profile.Phi << profile_int.Phi, profile_final.Phi;
    profile.dPhi << profile_int.dPhi, profile_final.dPhi;
    profile.Rerr = profile_final.Rerr;
    
    return profile;
  
  }

private:
  
//  const double NAN = std::numeric_limits<double>::quiet_NaN();
  const double INF = std::numeric_limits<double>::infinity();
  
  PROPERTY(double, phi_eps, 1e-3)
  PROPERTY(int, alpha, 2)
  
  PROPERTY(double, xtol, 1e-4)
  PROPERTY(double, phitol, 1e-4)
  PROPERTY(double, thinCutoff, .01)
  PROPERTY(double, npoints, 500)
  PROPERTY(double, rmin, 1e-4)
  PROPERTY(double, rmax, 1e4)
  
  PROPERTY(boost::uintmax_t, max_iter, 100)
  
};





}  // namespace PhaseTracer

#endif
