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

// File generated at Sat 11 Apr 2020 12:52:53

#include "THDMIISNMSSMBCsimple_soft_parameters.hpp"
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

const int THDMIISNMSSMBCsimple_soft_parameters::numberOfParameters;

THDMIISNMSSMBCsimple_soft_parameters::THDMIISNMSSMBCsimple_soft_parameters(const THDMIISNMSSMBCsimple_input_parameters& input_)
   : THDMIISNMSSMBCsimple_susy_parameters(input_)
{
   set_number_of_parameters(numberOfParameters);
}

THDMIISNMSSMBCsimple_soft_parameters::THDMIISNMSSMBCsimple_soft_parameters(
   const THDMIISNMSSMBCsimple_susy_parameters& susy_model
   , double M123_, double M5_, double M112_, double M222_, double M332_, double
   v1_, double v2_, double vS_
)
   : THDMIISNMSSMBCsimple_susy_parameters(susy_model)
   , M123(M123_), M5(M5_), M112(M112_), M222(M222_), M332(M332_), v1(v1_), v2(v2_)
   , vS(vS_)
{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd THDMIISNMSSMBCsimple_soft_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

THDMIISNMSSMBCsimple_soft_parameters THDMIISNMSSMBCsimple_soft_parameters::calc_beta(int loops) const
{
   double beta_M123 = 0.;
   double beta_M5 = 0.;
   double beta_M112 = 0.;
   double beta_M222 = 0.;
   double beta_M332 = 0.;
   double beta_v1 = 0.;
   double beta_v2 = 0.;
   double beta_vS = 0.;

   if (loops > 0) {
      const auto TRACE_STRUCT = CALCULATE_TRACES(loops);

      beta_M123 += calc_beta_M123_1_loop(TRACE_STRUCT);
      beta_M5 += calc_beta_M5_1_loop(TRACE_STRUCT);
      beta_M112 += calc_beta_M112_1_loop(TRACE_STRUCT);
      beta_M222 += calc_beta_M222_1_loop(TRACE_STRUCT);
      beta_M332 += calc_beta_M332_1_loop(TRACE_STRUCT);
      beta_v1 += calc_beta_v1_1_loop(TRACE_STRUCT);
      beta_v2 += calc_beta_v2_1_loop(TRACE_STRUCT);
      beta_vS += calc_beta_vS_1_loop(TRACE_STRUCT);

      if (loops > 1) {
         beta_M123 += calc_beta_M123_2_loop(TRACE_STRUCT);
         beta_M5 += calc_beta_M5_2_loop(TRACE_STRUCT);
         beta_M112 += calc_beta_M112_2_loop(TRACE_STRUCT);
         beta_M222 += calc_beta_M222_2_loop(TRACE_STRUCT);
         beta_M332 += calc_beta_M332_2_loop(TRACE_STRUCT);
         beta_v1 += calc_beta_v1_2_loop(TRACE_STRUCT);
         beta_v2 += calc_beta_v2_2_loop(TRACE_STRUCT);
         beta_vS += calc_beta_vS_2_loop(TRACE_STRUCT);

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


   const THDMIISNMSSMBCsimple_susy_parameters susy_betas(THDMIISNMSSMBCsimple_susy_parameters::calc_beta(loops));

   return THDMIISNMSSMBCsimple_soft_parameters(susy_betas, beta_M123, beta_M5, beta_M112, beta_M222, beta_M332, beta_v1, beta_v2, beta_vS);
}

THDMIISNMSSMBCsimple_soft_parameters THDMIISNMSSMBCsimple_soft_parameters::calc_beta() const
{
   return calc_beta(get_loops());
}

void THDMIISNMSSMBCsimple_soft_parameters::clear()
{
   THDMIISNMSSMBCsimple_susy_parameters::clear();

   M123 = 0.;
   M5 = 0.;
   M112 = 0.;
   M222 = 0.;
   M332 = 0.;
   v1 = 0.;
   v2 = 0.;
   vS = 0.;

}

Eigen::ArrayXd THDMIISNMSSMBCsimple_soft_parameters::get() const
{
   Eigen::ArrayXd pars(THDMIISNMSSMBCsimple_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(38) = M123;
   pars(39) = M5;
   pars(40) = M112;
   pars(41) = M222;
   pars(42) = M332;
   pars(43) = v1;
   pars(44) = v2;
   pars(45) = vS;


   return pars;
}

void THDMIISNMSSMBCsimple_soft_parameters::print(std::ostream& ostr) const
{
   THDMIISNMSSMBCsimple_susy_parameters::print(ostr);
   ostr << "soft parameters at Q = " << get_scale() << ":\n";
   ostr << "M123 = " << M123 << '\n';
   ostr << "M5 = " << M5 << '\n';
   ostr << "M112 = " << M112 << '\n';
   ostr << "M222 = " << M222 << '\n';
   ostr << "M332 = " << M332 << '\n';
   ostr << "v1 = " << v1 << '\n';
   ostr << "v2 = " << v2 << '\n';
   ostr << "vS = " << vS << '\n';

}

void THDMIISNMSSMBCsimple_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   THDMIISNMSSMBCsimple_susy_parameters::set(pars);

   M123 = pars(38);
   M5 = pars(39);
   M112 = pars(40);
   M222 = pars(41);
   M332 = pars(42);
   v1 = pars(43);
   v2 = pars(44);
   vS = pars(45);

}

THDMIISNMSSMBCsimple_soft_parameters::Soft_traces THDMIISNMSSMBCsimple_soft_parameters::calc_soft_traces(int loops) const
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

std::ostream& operator<<(std::ostream& ostr, const THDMIISNMSSMBCsimple_soft_parameters& soft_pars)
{
   soft_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
