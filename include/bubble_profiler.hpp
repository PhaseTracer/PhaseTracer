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

#ifndef PHASETRACER_BPROFILER_HPP_
#define PHASETRACER_BPROFILER_HPP_

#include <algorithm>
#include <ostream>
#include <vector>

#include <boost/cstdint.hpp>
#include <Eigen/Core>
#include "nlopt.hpp"

#include "logger.hpp"
// Include .hpp of BubbleProfiler
#include "BubbleProfiler/include/potential.hpp"
#include "BubbleProfiler/include/action.hpp"
#include "BubbleProfiler/include/error.hpp"
#include "BubbleProfiler/include/logging_manager.hpp"
#include "BubbleProfiler/include/shooting.hpp"
#include "BubbleProfiler/include/nlopt_optimizer.hpp"
#include "BubbleProfiler/include/perturbative_profiler.hpp"
#include "BubbleProfiler/include/generic_perturbative_profiler.hpp"
#include "BubbleProfiler/include/shooting_profile_guesser.hpp"
#include "BubbleProfiler/include/gsl_root_finder.hpp"
#include "BubbleProfiler/include/relative_convergence_tester.hpp"

namespace PhaseTracer {



class V_BubbleProfiler : public BubbleProfiler::Potential {
public:
  explicit V_BubbleProfiler(EffectivePotential::Potential &p_) :
    P(p_),
    n_fields(P.get_n_scalars()) {
    origin = Eigen::VectorXd::Zero(n_fields);
    origin_translation = origin;
    basis_transform = Eigen::MatrixXd::Identity(n_fields, n_fields);
  }
  virtual ~V_BubbleProfiler() = default;
  
  EffectivePotential::Potential &P;
  
  virtual V_BubbleProfiler * clone() const override {
     return new V_BubbleProfiler(*this);
  };
  
  double operator()(const Eigen::VectorXd& coords) const override{
    Eigen::VectorXd transformed_coords =
       (basis_transform * coords) + origin_translation;
    return P.V(transformed_coords,T);
  }
  
  double partial(const Eigen::VectorXd& coords, int i) const override{
    Eigen::VectorXd transformed_coords =
       (basis_transform * coords) + origin_translation;
    auto const dV = P.dV_dx(transformed_coords,T);
    return dV.coeff(i);
  }
  double partial(const Eigen::VectorXd& coords, int i, int j) const override{
    Eigen::VectorXd transformed_coords =
       (basis_transform * coords) + origin_translation;
    auto const d2V =  P.d2V_dx2(transformed_coords,T);
    return d2V.coeff(i, j);
  }
  std::size_t get_number_of_fields() const override{
    return P.get_n_scalars();
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
      Eigen::VectorXd transformed_coords =
          (basis_transform * coords) + origin_translation;
      return P.V(transformed_coords,TT);
     };

     BubbleProfiler::NLopt_optimizer optimizer(v, n_fields);

     // Don't need much precision for location of barrier
     optimizer.set_xtol_rel(1.e-5); // TODO: need lower than the one used to calcualte dV_dx and d2V_dx2
     optimizer.set_ftol_rel(1.e-5);

     optimizer.set_extremum_type(BubbleProfiler::NLopt_optimizer::Extremum_type::MAX);
     optimizer.set_lower_bounds(std::min(true_vacuum_loc(0),
                                         false_vacuum_loc(0)));
     optimizer.set_upper_bounds(std::max(true_vacuum_loc(0),
                                         false_vacuum_loc(0)));
     optimizer.set_max_time(1E6);

     Eigen::VectorXd initial_guess(0.5 * (true_vacuum_loc + false_vacuum_loc));
     const auto status = optimizer.optimize(initial_guess);

     if (!BubbleProfiler::optimization_succeeded(status)) {
        LOG(fatal) << "Error: unable to locate barrier. NLOPT status = "
                  << status;
        exit(EXIT_FAILURE);
     }

     return optimizer.get_extremum_location();
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

#endif  // PHASETRACER_BPROFILER_HPP_