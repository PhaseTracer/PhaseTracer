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
 * @file THDMIISNMSSMBCsimple_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Sat 11 Apr 2020 12:54:08 with FlexibleSUSY
 * 2.4.2 (git commit: a94199e5620b8684f5d30d0eece5757a5a72c4a4) and SARAH 4.14.3 .
 */

#ifndef THDMIISNMSSMBCsimple_TWO_SCALE_H
#define THDMIISNMSSMBCsimple_TWO_SCALE_H

#include "THDMIISNMSSMBCsimple_model.hpp"
#include "THDMIISNMSSMBCsimple_mass_eigenstates.hpp"

#include "model.hpp"

namespace flexiblesusy {

class Two_scale;
/**
 * @class THDMIISNMSSMBCsimple<Two_scale>
 * @brief model class with routines for determining masses and mixings and EWSB
 */
template<>
class THDMIISNMSSMBCsimple<Two_scale> : public Model, public THDMIISNMSSMBCsimple_mass_eigenstates {
public:
   explicit THDMIISNMSSMBCsimple(const THDMIISNMSSMBCsimple_input_parameters& input_ = THDMIISNMSSMBCsimple_input_parameters());
   THDMIISNMSSMBCsimple(const THDMIISNMSSMBCsimple&) = default;
   THDMIISNMSSMBCsimple(THDMIISNMSSMBCsimple&&) = default;
   virtual ~THDMIISNMSSMBCsimple() = default;
   THDMIISNMSSMBCsimple& operator=(const THDMIISNMSSMBCsimple&) = default;
   THDMIISNMSSMBCsimple& operator=(THDMIISNMSSMBCsimple&&) = default;

   // interface functions
   virtual void calculate_spectrum() override;
   virtual void clear_problems() override;
   virtual std::string name() const override;
   virtual void run_to(double scale, double eps = -1.0) override;
   virtual void print(std::ostream& out = std::cerr) const override;
   virtual void set_precision(double) override;
};

std::ostream& operator<<(std::ostream&, const THDMIISNMSSMBCsimple<Two_scale>&);

} // namespace flexiblesusy

#endif
