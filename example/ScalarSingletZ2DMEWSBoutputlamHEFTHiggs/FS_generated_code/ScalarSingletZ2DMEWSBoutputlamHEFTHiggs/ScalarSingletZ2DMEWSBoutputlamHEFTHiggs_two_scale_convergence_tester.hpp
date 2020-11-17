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

// File generated at Tue 17 Nov 2020 15:32:58

#ifndef ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_TWO_SCALE_CONVERGENCE_TESTER_H
#define ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_TWO_SCALE_CONVERGENCE_TESTER_H

#include "ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_convergence_tester.hpp"
#include "ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_two_scale_model.hpp"

#include "convergence_tester_drbar.hpp"

namespace flexiblesusy {

class Two_scale;

template<>
class ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_convergence_tester<Two_scale> : public Convergence_tester_DRbar<ScalarSingletZ2DMEWSBoutputlamHEFTHiggs<Two_scale> > {
public:
   using Scale_getter = Convergence_tester_DRbar<ScalarSingletZ2DMEWSBoutputlamHEFTHiggs<Two_scale>>::Scale_getter;

   ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_convergence_tester(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs<Two_scale>*, double, const Scale_getter& sg = Scale_getter());
   virtual ~ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_convergence_tester() = default;

protected:
   virtual double max_rel_diff() const;
};

} // namespace flexiblesusy

#endif
