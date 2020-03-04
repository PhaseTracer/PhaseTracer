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

// File generated at Thu 7 Nov 2019 18:53:56

#ifndef THDMIISNMSSMBCsimple_TWO_SCALE_HIGH_SCALE_CONSTRAINT_H
#define THDMIISNMSSMBCsimple_TWO_SCALE_HIGH_SCALE_CONSTRAINT_H

#include "THDMIISNMSSMBCsimple_high_scale_constraint.hpp"
#include "THDMIISNMSSMBCsimple_input_parameters.hpp"
#include "single_scale_constraint.hpp"

namespace flexiblesusy {

template <class T>
class THDMIISNMSSMBCsimple;

class Two_scale;

template<>
class THDMIISNMSSMBCsimple_high_scale_constraint<Two_scale> : public Single_scale_constraint {
public:
   THDMIISNMSSMBCsimple_high_scale_constraint() = default;
   THDMIISNMSSMBCsimple_high_scale_constraint(THDMIISNMSSMBCsimple<Two_scale>*);
   virtual ~THDMIISNMSSMBCsimple_high_scale_constraint() = default;
   virtual void apply() override;
   virtual double get_scale() const override;
   virtual std::string name() const override { return "THDMIISNMSSMBCsimple high-scale constraint"; }
   virtual void set_model(Model*) override;

   void clear();
   double get_initial_scale_guess() const;
   const THDMIISNMSSMBCsimple_input_parameters& get_input_parameters() const;
   THDMIISNMSSMBCsimple<Two_scale>* get_model() const;
   void initialize();
   void set_scale(double); ///< fix unification scale (0 = unfixed)

protected:
   void update_scale();
   bool check_non_perturbative();

private:
   double scale{0.};
   double initial_scale_guess{0.};
   THDMIISNMSSMBCsimple<Two_scale>* model{nullptr};

   void check_model_ptr() const;
};

} // namespace flexiblesusy

#endif
