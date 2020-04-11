// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Sat 11 Apr 2020 12:54:08

/**
 * @file THDMIISNMSSMBCsimple_two_scale_ewsb_solver.hpp
 *
 * @brief contains class for solving EWSB when two-scale algorithm is used
 *
 * This file was generated at Sat 11 Apr 2020 12:54:08 with FlexibleSUSY
 * 2.4.2 (git commit: a94199e5620b8684f5d30d0eece5757a5a72c4a4) and SARAH 4.14.3 .
 */

#ifndef THDMIISNMSSMBCsimple_TWO_SCALE_EWSB_SOLVER_H
#define THDMIISNMSSMBCsimple_TWO_SCALE_EWSB_SOLVER_H

#include "THDMIISNMSSMBCsimple_ewsb_solver.hpp"
#include "THDMIISNMSSMBCsimple_ewsb_solver_interface.hpp"
#include "error.hpp"

#include <Eigen/Core>

namespace flexiblesusy {

class EWSB_solver;
class Two_scale;

class THDMIISNMSSMBCsimple_mass_eigenstates;

template<>
class THDMIISNMSSMBCsimple_ewsb_solver<Two_scale> : public THDMIISNMSSMBCsimple_ewsb_solver_interface {
public:
   THDMIISNMSSMBCsimple_ewsb_solver() = default;
   THDMIISNMSSMBCsimple_ewsb_solver(const THDMIISNMSSMBCsimple_ewsb_solver&) = default;
   THDMIISNMSSMBCsimple_ewsb_solver(THDMIISNMSSMBCsimple_ewsb_solver&&) = default;
   virtual ~THDMIISNMSSMBCsimple_ewsb_solver() {}
   THDMIISNMSSMBCsimple_ewsb_solver& operator=(const THDMIISNMSSMBCsimple_ewsb_solver&) = default;
   THDMIISNMSSMBCsimple_ewsb_solver& operator=(THDMIISNMSSMBCsimple_ewsb_solver&&) = default;

   virtual void set_loop_order(int l) override { loop_order = l; }
   virtual void set_number_of_iterations(int n) override { number_of_iterations = n; }
   virtual void set_precision(double p) override { precision = p; }

   virtual int get_loop_order() const override { return loop_order; }
   virtual int get_number_of_iterations() const override { return number_of_iterations; }
   virtual double get_precision() const override { return precision; }

   virtual int solve(THDMIISNMSSMBCsimple_mass_eigenstates&) override;
private:
   static const int number_of_ewsb_equations = 3;
   using EWSB_vector_t = Eigen::Matrix<double,number_of_ewsb_equations,1>;

   class EEWSBStepFailed : public Error {
   public:
      EEWSBStepFailed() : Error("Could not perform EWSB step") {}
      virtual ~EEWSBStepFailed() = default;
   };

   int number_of_iterations{100}; ///< maximum number of iterations
   int loop_order{2};             ///< loop order to solve EWSB at
   double precision{1.e-5};       ///< precision goal

   void set_ewsb_solution(THDMIISNMSSMBCsimple_mass_eigenstates&, const EWSB_solver*);
   template <typename It> void set_best_ewsb_solution(THDMIISNMSSMBCsimple_mass_eigenstates&, It, It);

   int solve_tree_level(THDMIISNMSSMBCsimple_mass_eigenstates&);
   int solve_iteratively(THDMIISNMSSMBCsimple_mass_eigenstates&);
   int solve_iteratively_at(THDMIISNMSSMBCsimple_mass_eigenstates&, int);
   int solve_iteratively_with(THDMIISNMSSMBCsimple_mass_eigenstates&, EWSB_solver*, const EWSB_vector_t&);

   EWSB_vector_t initial_guess(const THDMIISNMSSMBCsimple_mass_eigenstates&) const;
   EWSB_vector_t tadpole_equations(const THDMIISNMSSMBCsimple_mass_eigenstates&) const;
   EWSB_vector_t ewsb_step(const THDMIISNMSSMBCsimple_mass_eigenstates&) const;
};

} // namespace flexiblesusy

#endif
