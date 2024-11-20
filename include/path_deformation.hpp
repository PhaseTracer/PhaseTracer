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

#ifndef PHASETRACER_PATH_DEFORMATION_HPP_
#define PHASETRACER_PATH_DEFORMATION_HPP_

#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>
#include <tuple>
#include <vector>

#include <boost/cstdint.hpp>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <nlopt.hpp>

#include "logger.hpp"
#include "shooting.hpp"

namespace PhaseTracer {



template<typename T>
std::vector<T> cumulative_trapezoidal_integration(const std::vector<T>& values)
{
    std::vector<T> integration_result;
    integration_result.clear();
    T sum = T(0);
    for (size_t i = 0; i < values.size() - 1; ++i)
    {
        sum += (values[i] + values[i + 1]) / 2.0;
        integration_result.push_back(sum);
    }
    return integration_result;
}

/* Fit a spline to a path in field space, and find the potential on that path */
class SplinePath: public PotentialForShooting {

public:
  explicit SplinePath(EffectivePotential::Potential &potential,
                      double T_,
                      std::vector<Eigen::VectorXd> pts_,
                      bool extend_to_minima_ = true,
                      bool reeval_distances_ = true);
  
  virtual ~SplinePath() = default;
  
  /* Calculates to 4th order if len(phi) >= 5, otherwise 1st/2nd order. */
  std::vector<Eigen::VectorXd> _pathDeriv(const std::vector<Eigen::VectorXd> phi);
  

  double find_loc_min_w_guess(Eigen::VectorXd p0, Eigen::VectorXd dp0, double guess = 0.);
  
  void get_path_tck(std::vector<double> pdist);
  
  void dpdx(const double& x, double& dpdx_, double r);
  
  Eigen::VectorXd vecp(double x);
  
  double V(double x) const override;
  
  double dV(double x) const override;
  
  double d2V(double x) const override;
  
  double get_path_length(){return length;}
  
private:
  
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
  
//  TransitionFinder& tf;
  std::vector<Eigen::VectorXd> pts;
  int V_spline_samples = 140; // larger than 5
  bool extend_to_minima;
  bool reeval_distances;
  
  double length;
  
  std::vector<alglib::spline1dinterpolant> _path_tck;
  alglib::spline1dinterpolant _V_tck;
  
  EffectivePotential::Potential &P;
  double T;
  
  size_t nphi;
  size_t num_nodes;
  
};


/* ************************************************************ */


struct FullTunneling {
  double action;
  double fRatio;
  Profile1D profile1D;
  std::vector<Eigen::VectorXd> phi;
  std::vector<std::vector<std::vector<Eigen::VectorXd>>> saved_steps;
};

class PathDeformation {
public:
  
  explicit PathDeformation(EffectivePotential::Potential &potential): P(potential), nphi(potential.get_n_scalars()) {}
  
  virtual ~PathDeformation() = default;
  
  bool deformPath(std::vector<double> dphidr);
  void step(double &lastStep, bool &step_reversed, double &fRatio);
  
  // Calculate the normal force and potential gradient on the path
  void forces(std::vector<Eigen::VectorXd>& F_norm, std::vector<Eigen::VectorXd>& dV);
  
  std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> Nbspld2(std::vector<double> t, std::vector<double> x, int k = 3);
  
  /* Calculate the instanton solution in multiple field dimension */
  FullTunneling full_tunneling(std::vector<Eigen::VectorXd> path_pts);
  
  double get_action(){return action_temp;}
  
  
private:

  EffectivePotential::Potential &P;
  
  double phi_absMin;
  double phi_metaMin;
  
  PROPERTY(double, T, 0);
  
  /* Number of dimensions */
  PROPERTY(size_t, num_dims, 3)
  
  /* Number of basis splines to use */
  PROPERTY(size_t, nb, 10);
  /* Order of basis splines */
  PROPERTY(size_t, kb, 3);
  /* Get each step saved */
  PROPERTY(bool, save_all_steps, false);
  /* The smallest the square of dphidr is allowed to be */
  PROPERTY(double, v2min, 0.0);
  /* Maximum number of steps to take in a deformation */
  PROPERTY(size_t, step_maxiter, 500);
  /* Maximum number of allowed deformation iterations */
  PROPERTY(size_t, path_maxiter, 20);
  /* Number of samples to take along the path to create the spline
   interpolation functions */
  PROPERTY(size_t, V_spline_samples, 100);
  /* Flag to extend the path to minimums*/
  PROPERTY(bool, extend_to_minima, true);
  
  /** Pass through the shooting settings **/
  /* The precision of field values after taking the logarithm */
  PROPERTY(double, xtol, 1e-4)
  /* The fractional error tolerance in integration*/
  PROPERTY(double, phitol, 1e-4)
  /* The cut off for finding the initial conditions for integration */
  PROPERTY(double, thin_cutoff, .01)
  /* Number of points to return in the profile */
  PROPERTY(double, npoints, 500)
  /* The smallest starting radius */
  PROPERTY(double, rmin, 1e-4)
  /* The maximum allowed integration distance */
  PROPERTY(double, rmax, 1e4)
  /* The maximum number of points to be positioned during the integration process. */
  PROPERTY(boost::uintmax_t, max_iter, 100)
  
  
  size_t num_nodes;
  size_t nphi;
  size_t num_steps=0;
  
  double startstep=2e-3;
  double fRatioConv=.02;
  double converge_0=5.;
  double fRatioIncrease=5.;
  
  double maxstep=.1;
  double minstep=1e-4;
  double reverseCheck=.15;
  double stepIncrease=1.5;
  double stepDecrease=5.;
  bool checkAfterFit=true;

  Eigen::MatrixXd X_node;
  Eigen::MatrixXd dX_node;
  Eigen::MatrixXd d2X_node;
  std::vector<double> t_node;
  std::vector<Eigen::VectorXd> beta_node;
  std::vector<Eigen::VectorXd> phi_node;
  std::vector<double> v2_node;
  double totalLength_node;
  
  bool breakLoop = false;
  bool fix_start=false;
  bool fix_end=false;
  std::vector<Eigen::VectorXd> phi_prev;
  std::vector<Eigen::VectorXd> F_prev;
  double action_temp=std::numeric_limits<double>::quiet_NaN();
  
  std::vector<std::vector<Eigen::VectorXd>> phi_list;
  std::vector<std::vector<Eigen::VectorXd>> F_list;

};

}  // namespace PhaseTracer

#endif  // PHASETRACER_PATH_DEFORMATION_HPP_