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

#ifndef ACTION_CALCULATOR_HPP_INCLUDED
#define ACTION_CALCULATOR_HPP_INCLUDED

//#include <cmath>
#include <ostream>
//#include <vector>
//#include <Eigen/Core>
//#include <algorithm>
//#include <boost/cstdint.hpp>
//#include <boost/math/special_functions/bessel.hpp>
//#include <boost/numeric/odeint.hpp>
//#include <boost/math/quadrature/gauss_kronrod.hpp>
//#include <boost/math/tools/minima.hpp>
//#include <gsl/gsl_sf_bessel.h>


#include "bubble_profiler.hpp"
#include "path_deformation.hpp"

namespace PhaseTracer {


class ActionCalculator{
public:
  explicit ActionCalculator(EffectivePotential::Potential &potential_) :
    potential(potential_){
  }
  virtual ~ActionCalculator() = default;
  
  double get_action(const Eigen::VectorXd& vacuum_1, const Eigen::VectorXd& vacuum_2, double T) const{
    Eigen::VectorXd true_vacuum = vacuum_1;
    Eigen::VectorXd false_vacuum = vacuum_2;
    if ( potential.V(true_vacuum,T) > potential.V(false_vacuum,T) ) true_vacuum.swap(false_vacuum);
    
    double action_BP=std::numeric_limits<double>::quiet_NaN();
    if (use_BubbleProfiler){
      V_BubbleProfiler V_BP(potential); // perturbative_profiler only accept non-const potential
      V_BP.set_T(T); // This is necessary!!
      
      if (potential.get_n_scalars() == 1  && !use_perturbative) {

        double false_min = false_vacuum[0];
        double true_min = true_vacuum[0];
        auto barrier_ = V_BP.find_one_dimensional_barrier( true_vacuum, false_vacuum, T);
        double barrier = barrier_[0];

        LOG(debug) << "Calculate action(BP) at T=" << T << ", with false, true, barrier = " << false_min << ", " << true_min << ", " << barrier;

        try{
          BubbleProfiler::Shooting one_dim;
          one_dim.solve(V_BP, false_min,
                        true_min,barrier, num_dims, BubbleProfiler::Shooting::Solver_options::Compute_action);
          action_BP = one_dim.get_euclidean_action();
        }catch (const std::exception& e) {
          LOG(warning) << "At T=" << T << ", between[" << false_min << "] and [" << true_min << "]: "   << e.what() << std::endl;
        }
      } else {
        LOG(debug) << "Calculate action(BP) at T=" << T << ", between [" << false_vacuum.transpose().format(Eigen::IOFormat(4, Eigen::DontAlignCols, " ", " ")) << "] and [" << true_vacuum.transpose().format(Eigen::IOFormat(4, Eigen::DontAlignCols, " ", " ")) << "]";
        BubbleProfiler::RK4_perturbative_profiler profiler;

  //      profiler.set_domain_start(input.domain_start);
  //      profiler.set_domain_end(input.domain_end);
        profiler.set_initial_step_size(1.e-2);
        profiler.set_interpolation_points_fraction(1.0);


        profiler.set_false_vacuum_loc(false_vacuum);
        profiler.set_true_vacuum_loc(true_vacuum);
        profiler.set_number_of_dimensions(num_dims);
        auto root_finder = std::make_shared<BubbleProfiler::GSL_root_finder<Eigen::Dynamic> >();
        profiler.set_root_finder(root_finder);
        std::shared_ptr<BubbleProfiler::Profile_guesser> guesser;
        guesser = std::make_shared<BubbleProfiler::Kink_profile_guesser>();
        profiler.set_initial_guesser(guesser);
        auto convergence_tester = std::make_shared<BubbleProfiler::Relative_convergence_tester>(
                                  1.e-3, 1.e-3);
        profiler.set_convergence_tester(convergence_tester);

        try{
          profiler.calculate_bubble_profile(V_BP);
          action_BP = profiler.get_euclidean_action();
        }catch (const std::exception& e) {
          LOG(warning) << "At T=" << T <<  ", between [" << false_vacuum.transpose().format(Eigen::IOFormat(4, Eigen::DontAlignCols, " ", " ")) << "] and [" << true_vacuum.transpose().format(Eigen::IOFormat(4, Eigen::DontAlignCols, " ", " ")) << "]: "   << e.what() << std::endl;
        }
      }
      LOG(debug) << " S(BP) = " << action_BP << std::endl;

    }
    double action_PD=std::numeric_limits<double>::quiet_NaN();
    if (use_PathDeformation){
      LOG(debug) << "Calculate action(PD) at T=" << T << ", between [" << false_vacuum.transpose().format(Eigen::IOFormat(4, Eigen::DontAlignCols, " ", " ")) << "] and [" << true_vacuum.transpose().format(Eigen::IOFormat(4, Eigen::DontAlignCols, " ", " ")) << "]";
      try{
        std::vector<Eigen::VectorXd> path_pts;
        path_pts.push_back(true_vacuum);
        path_pts.push_back(false_vacuum);
        PathDeformation pd(potential);
        auto a_pd = pd.fullTunneling(path_pts);
      }catch (const std::exception& e) {
        LOG(warning) << "At T=" << T <<  ", between [" << false_vacuum.transpose().format(Eigen::IOFormat(4, Eigen::DontAlignCols, " ", " ")) << "] and [" << true_vacuum.transpose().format(Eigen::IOFormat(4, Eigen::DontAlignCols, " ", " ")) << "]: "   << e.what() << std::endl;
      }
      
      action_PD = 1.;
    }
    
    
    LOG(debug) << " S(PD) = " << action_PD << std::endl;

    return action_BP;

  }
  
private:
  EffectivePotential::Potential &potential;
  
  PROPERTY(size_t, num_dims, 3)
  
  PROPERTY(bool, use_BubbleProfiler, true)
  PROPERTY(bool, use_perturbative, false)
  
  PROPERTY(bool, use_PathDeformation, true)
  
  
  
  
  
};


}  // namespace PhaseTracer

#endif