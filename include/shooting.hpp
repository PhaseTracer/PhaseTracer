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

#ifndef PHASETRACER_SHOOTING_HPP_INCLUDED
#define PHASETRACER_SHOOTING_HPP_INCLUDED

#include <cmath>
#include <ostream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <boost/cstdint.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/tools/minima.hpp>
#include <gsl/gsl_sf_bessel.h>

#include "potential.hpp"
#include "property.hpp"
#include "logger.hpp"


namespace PhaseTracer {


class PotentialForShooting{
public:
  virtual double V(double phi) const = 0;
  virtual double dV(double phi) const = 0;
  virtual double d2V(double phi) const = 0;
};

class CubicInterpFunction {
public:
  CubicInterpFunction(double y0, double dy0, double y1, double dy1, double c_ = 0)
    : c(c_) {
      x0 = y0;
      x1 = y0 + dy0 / 3.0;
      x2 = y1 - dy1 / 3.0;
      x3 = y1;
    }
  double operator()(double t) const {
    double mt = 1. - t;
    return x0 * pow(mt, 3) + 3.0 * x1 * mt * mt * t + 3.0 * x2 * mt * t * t + x3 * pow(t, 3) - c;
  }
private:
  double x0;
  double x1;
  double x2;
  double x3;
  double c;
};

struct Profile1D {
  Eigen::VectorXd R;
  Eigen::VectorXd Phi;
  Eigen::VectorXd dPhi;
  double Rerr;
};


class Shooting {
public:
  explicit Shooting(PotentialForShooting& ps_) :
    ps(ps_) {
      // TODO Add checking, this olny works for 1d
    }
  virtual ~Shooting() = default;
    
  /*Calculates `dV/dphi` at ``phi = phi_absMin + delta_phi``.*/
  double dV_from_absMin(double delta_phi);
  /* Find edge of the potential barrier. */
  void findBarrierLocation();
  /* Find the characteristic length scale for tunneling over the potential barrier. */
  void findRScale();
  /*Find `phi(r)` given `phi(r=0)`, assuming a quadratic potential.*/
  void exactSolution(double r, double phi0, double dV_, double d2V_,
                     double* phi_r, double* dphi_r);
  /*Finds the initial conditions, phi(r0) = phi_cutoff, for integration.*/
  void initialConditions(double delta_phi0, double rmin, double delta_phi_cutoff,
                         double* r0, double* phi_r0, double* dphi_r0);
  Eigen::Vector2d equationOfMotion(const Eigen::Vector2d y, const double r){
    Eigen::Vector2d dydr;
    dydr[0] = y[1];
    dydr[1] = ps.dV(y[0]) - alpha * y[1] /r ;
    return dydr;
  }
  /* Integrate the bubble wall equation */
  int integrateProfile(double r0, std::vector<double> y0, double* rf, std::vector<double>* yf,
      double dr0, std::vector<double> epsabs, std::vector<double> epsfrac, double drmin, double rmax);
  /* Integrate the bubble profile, saving the output in an array */
  Profile1D integrateAndSaveProfile(Eigen::VectorXd R, std::vector<double> y0,
      double dr0, std::vector<double> epsabs, std::vector<double> epsfrac, double drmin);
  /* Calculate the bubble profile */
  Profile1D findProfile(double metaMin, double absMin, double xguess=NAN, int max_interior_pts=0);
  /* Calculate the Euclidean action for the instanton */
  double calAction(Profile1D profile);
  
  /* Get linearly spaced phi */
  void evenlySpacedPhi(Profile1D pf, std::vector<double>* p,std::vector<double>* dp,
                       size_t npoints=100, bool fixAbs=true);
  
private:
  
  PotentialForShooting& ps;
  
  double phi_absMin;
  double phi_metaMin;
  double phi_bar;
  double rscale; // approximate radial scale of the instanton
  
  PROPERTY(double, phi_eps_rel, 1e-3)
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
