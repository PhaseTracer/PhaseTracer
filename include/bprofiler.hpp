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

#ifndef PHASETRACER_BPROFILER_HPP_INCLUDED
#define PHASETRACER_BPROFILER_HPP_INCLUDED

#include <ostream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <boost/cstdint.hpp>
#include "nlopt.hpp"


#include "logger.hpp"
// Include .hpp of BubbleProfiler
#include "BubbleProfiler/include/potential.hpp"
#include "BubbleProfiler/include/action.hpp"
#include "BubbleProfiler/include/error.hpp"
#include "BubbleProfiler/include/logging_manager.hpp"
#include "BubbleProfiler/include/shooting.hpp"
#include "BubbleProfiler/include/nlopt_optimizer.hpp"

namespace PhaseTracer {



class V_BubbleProfiler : public BubbleProfiler::Potential {
public:
  explicit V_BubbleProfiler(PhaseFinder& pf_) :
    pf(pf_),
    n_fields(pf_.P.get_n_scalars()) {
    origin = Eigen::VectorXd::Zero(n_fields);
    origin_translation = origin;
    basis_transform = Eigen::MatrixXd::Identity(n_fields, n_fields);
  }
  virtual ~V_BubbleProfiler() = default;
  
  PhaseFinder pf;
  
  virtual V_BubbleProfiler * clone() const override {
     return new V_BubbleProfiler(*this);
  };
  
  double operator()(const Eigen::VectorXd& coords) const override{
    return pf.P.V(coords,T);
  }
  
  double partial(const Eigen::VectorXd& coords, int i) const override{
    auto const dV = pf.P.dV_dx(coords,T);
    return dV.coeff(i);
  }
  double partial(const Eigen::VectorXd& coords, int i, int j) const override{
    auto const d2V =  pf.P.d2V_dx2(coords,T);
    return d2V.coeff(i, j);
  }
  std::size_t get_number_of_fields() const override{
    return pf.P.get_n_scalars();
  }

  void translate_origin(const Eigen::VectorXd& translation) override{
      origin_translation = translation;
  }
  void apply_basis_change(const Eigen::MatrixXd& new_basis) override{
      basis_transform = basis_transform * (new_basis.transpose());
  }
  void add_constant_term(double constant) override{
      constant_term += constant;
  }
  
  void set_T(double T_) const { const_cast<double&>(T) = T_; }
  
  
  Eigen::VectorXd find_one_dimensional_barrier(
     const Eigen::VectorXd& true_vacuum_loc,
     const Eigen::VectorXd& false_vacuum_loc,
     double TT) const
  {
     if (n_fields != 1) {
       LOG(fatal) << ("automatically locating potential barrier only "
                      "supported for single field case");
     }

     const auto v = [this,TT](const Eigen::VectorXd& coords) {
      return pf.P.V(coords,TT);
     };

     BubbleProfiler::NLopt_optimizer optimizer(v, n_fields);

     // Don't need much precision for location of barrier
     optimizer.set_xtol_rel(1.e-2);
     optimizer.set_ftol_rel(1.e-2);

     optimizer.set_extremum_type(BubbleProfiler::NLopt_optimizer::Extremum_type::MAX);
     optimizer.set_lower_bounds(std::min(true_vacuum_loc(0),
                                         false_vacuum_loc(0)));
     optimizer.set_upper_bounds(std::max(true_vacuum_loc(0),
                                         false_vacuum_loc(0)));
     optimizer.set_max_time(10);

     Eigen::VectorXd initial_guess(0.5 * (true_vacuum_loc + false_vacuum_loc));
     const auto status = optimizer.optimize(initial_guess);

     if (!BubbleProfiler::optimization_succeeded(status)) {
        LOG(fatal) << "Error: unable to locate barrier. NLOPT status = "
                  << status;
        exit(EXIT_FAILURE);
     }

     return optimizer.get_extremum_location();
  }
  
  double get_action(const Eigen::VectorXd& vacuum_1, const Eigen::VectorXd& vacuum_2, double T) const{
    
    double action;
    
    Eigen::VectorXd true_vacuum = vacuum_1;
    Eigen::VectorXd false_vacuum = vacuum_2;
    if ( pf.P.V(true_vacuum,T) > pf.P.V(false_vacuum,T) ) true_vacuum.swap(false_vacuum);
//    LOG(debug) << "true vacuum";
    
    if (n_fields == 1) {
    
      double false_min = false_vacuum[0];
      double true_min = true_vacuum[0];
      auto barrier_ = find_one_dimensional_barrier( true_vacuum, false_vacuum, T);
      double barrier = barrier_[0];
      
      LOG(debug) << "Calculate action at " << T << ", with false, true, barrier = " << false_min << ", " << true_min << ", " << barrier;
      
      const auto potential = [this, T] (double phi) {
        const Eigen::VectorXd vector_phi = Eigen::VectorXd::Constant(1, phi);
        return pf.P.V(vector_phi,T);
      };

  //    std::cout << "potential 1  = " << potential(true_min) << std::endl;
  //    V_BP.set_T(T);
  //    std::cout << "potential 1  = " << V_BP(true_vacuum) << std::endl;
      
      const auto potential_first = [this, T] (double phi) {
        const Eigen::VectorXd vector_phi = Eigen::VectorXd::Constant(1, phi);
        const auto dV_dx = pf.P.dV_dx(vector_phi,T);
        return dV_dx.coeff(0);
      };

      const auto potential_second = [this, T] (double phi) {
        const Eigen::VectorXd vector_phi = Eigen::VectorXd::Constant(1, phi);
        const auto d2V_dx2 = pf.P.d2V_dx2(vector_phi,T);
        return d2V_dx2.coeff(0, 0);
      };

      try{
        BubbleProfiler::Shooting one_dim;
        one_dim.solve(potential, potential_first, potential_second,
                      false_min, true_min+1E-6, barrier,
                      4, BubbleProfiler::Shooting::Solver_options::Compute_action);
        action =one_dim.get_euclidean_action();
      }catch (const std::exception& e) {
        LOG(warning) << "At T=" << T << ", between[" << false_min << "] and [" << true_min << "]: "   << e.what();
      }
    } else {
      LOG(fatal) << "Action calculation for n_scalars != 1 is not ready!";
    }
    
    LOG(debug) << "S = " << action << ", S/T = " << action/T << std::endl;
    
    return action;
  }
  
private:
  mutable double T;
  
  std::size_t n_fields{0};
  double constant_term{0.};
  Eigen::VectorXd origin{};
  Eigen::VectorXd origin_translation{};
  Eigen::VectorXd mu{};
  Eigen::MatrixXd basis_transform{};
  
};





}  // namespace PhaseTracer

#endif
