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

// File generated at Tue 17 Nov 2020 15:32:51

#include "ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters.hpp"
#include "config.h"
#ifdef ENABLE_THREADS
#include "global_thread_pool.hpp"
#endif
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces
#define TRACE_STRUCT_TYPE Soft_traces
#define CALCULATE_TRACES(l) calc_soft_traces(l);

const int ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters::numberOfParameters;

ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters::ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters(const ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_input_parameters& input_)
   : ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_parameters(input_)
{
   set_number_of_parameters(numberOfParameters);
}

ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters::ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters(
   const ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_parameters& susy_model
   , double muS2_, double muH2_, double v_
)
   : ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_parameters(susy_model)
   , muS2(muS2_), muH2(muH2_), v(v_)
{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters::calc_beta(int loops) const
{
   double beta_muS2 = 0.;
   double beta_muH2 = 0.;
   double beta_v = 0.;

   if (loops > 0) {
      const auto TRACE_STRUCT = CALCULATE_TRACES(loops);

      beta_muS2 += calc_beta_muS2_1_loop(TRACE_STRUCT);
      beta_muH2 += calc_beta_muH2_1_loop(TRACE_STRUCT);
      beta_v += calc_beta_v_1_loop(TRACE_STRUCT);

      if (loops > 1) {
         beta_muS2 += calc_beta_muS2_2_loop(TRACE_STRUCT);
         beta_muH2 += calc_beta_muH2_2_loop(TRACE_STRUCT);
         beta_v += calc_beta_v_2_loop(TRACE_STRUCT);

         if (loops > 2) {
         #ifdef ENABLE_THREADS
            {


            }
         #else
         #endif

            if (loops > 3) {

               if (loops > 4) {

               }
            }
         }
      }
   }


   const ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_parameters susy_betas(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_parameters::calc_beta(loops));

   return ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters(susy_betas, beta_muS2, beta_muH2, beta_v);
}

ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters::calc_beta() const
{
   return calc_beta(get_loops());
}

void ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters::clear()
{
   ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_parameters::clear();

   muS2 = 0.;
   muH2 = 0.;
   v = 0.;

}

Eigen::ArrayXd ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters::get() const
{
   Eigen::ArrayXd pars(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(33) = muS2;
   pars(34) = muH2;
   pars(35) = v;


   return pars;
}

void ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters::print(std::ostream& ostr) const
{
   ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_parameters::print(ostr);
   ostr << "soft parameters at Q = " << get_scale() << ":\n";
   ostr << "muS2 = " << muS2 << '\n';
   ostr << "muH2 = " << muH2 << '\n';
   ostr << "v = " << v << '\n';

}

void ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_parameters::set(pars);

   muS2 = pars(33);
   muH2 = pars(34);
   v = pars(35);

}

ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters::Soft_traces ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters::calc_soft_traces(int loops) const
{
   Soft_traces soft_traces;

   if (loops > 0) {
      

      TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());

   }

   if (loops > 1) {
      TRACE_STRUCT.traceYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint()).trace()
         );
      TRACE_STRUCT.traceYdAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace()
         );
      TRACE_STRUCT.traceYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace()
         );
      TRACE_STRUCT.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint()).trace()
         );

   }

   if (loops > 2) {

   }

   return soft_traces;
}

std::ostream& operator<<(std::ostream& ostr, const ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters& soft_pars)
{
   soft_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
