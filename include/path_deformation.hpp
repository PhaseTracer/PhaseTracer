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
#include <boost/math/special_functions/bessel.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/tools/minima.hpp>
#include <gsl/gsl_sf_bessel.h>


#include "logger.hpp"
#include "transition_finder.hpp"

namespace PhaseTracer {


//struct Profile1D {
//  Eigen::VectorXd R;
//  Eigen::VectorXd Phi;
//  Eigen::VectorXd dPhi;
//  double Rerr;
//};

/* Fit a spline to a path in field space, and find the potential on that path */
//class SplinePath:



class PathDeformation {
public:
  explicit PathDeformation(TransitionFinder& tf_) :
    tf(tf_) {
      // TODO
    }
  virtual ~PathDeformation() = default;
  
  double V(double phi){ // TODO this must be replcaed
    return 0.25*pow(phi,4) - 0.49*pow(phi,3) + 0.235 * pow(phi,2);
//    return 0.25*pow(phi,4) - 0.4*pow(phi,3) + 0.1 * pow(phi,2);
  }
  double dV(double phi){  // TODO this must be replcaed
    return phi*(phi-.47)*(phi-1);
//    return phi*(phi-.2)*(phi-1);
  }
  
  double d2V(double phi){  // TODO this must be replcaed
    return (phi-.47)*(phi-1) + phi*(phi-1) + phi*(phi-.47);
//    return (phi-.2)*(phi-1) + phi*(phi-1) + phi*(phi-.2);
  }
  
  /* Calculate the instanton solution in multiple field dimension */
  double fullTunneling(double delta_phi){
    std::vector<Eigen::VectorXd> path_pts;
    path_pts.push_back(Eigen::VectorXd(2));
    path_pts.push_back(Eigen::VectorXd(2));
    path_pts[0] << 1, 1;
    path_pts[1] << 0, 0;
    
    
    int maxiter = 1;
    
    for (int num_iter = 1; num_iter <= maxiter; num_iter++) {
      LOG(debug) <<  "Starting tunneling step " << num_iter;
    }
    
    
    
    return 0;
  }
  
private:
  
  TransitionFinder& tf;
  
  double phi_absMin;
  double phi_metaMin;
  
};

}  // namespace PhaseTracer

#endif
