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


/**
 * @file ScalarSingletZ2DMMhInputMsInput_ewsb_solver_interface.hpp
 *
 * @brief contains interface for EWSB solver
 *
 * This file was generated with FlexibleSUSY 2.5.0 and SARAH 4.14.3 .
 */

#ifndef ScalarSingletZ2DMMhInputMsInput_EWSB_SOLVER_INTERFACE_H
#define ScalarSingletZ2DMMhInputMsInput_EWSB_SOLVER_INTERFACE_H

namespace flexiblesusy {

class ScalarSingletZ2DMMhInputMsInput_mass_eigenstates;

/**
 * @class ScalarSingletZ2DMMhInputMsInput_ewsb_solver_interface
 * @brief interface for EWSB solvers to be used to solve the EWSB equations
 */
class ScalarSingletZ2DMMhInputMsInput_ewsb_solver_interface {
public:
   virtual ~ScalarSingletZ2DMMhInputMsInput_ewsb_solver_interface() {}

   virtual void set_loop_order(int) = 0;
   virtual void set_number_of_iterations(int) = 0;
   virtual void set_precision(double) = 0;

   virtual int get_loop_order() const = 0;
   virtual int get_number_of_iterations() const = 0;
   virtual double get_precision() const = 0;

   virtual int solve(ScalarSingletZ2DMMhInputMsInput_mass_eigenstates&) = 0;
};

} // namespace flexiblesusy

#endif
