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

// File generated at Sat 11 Apr 2020 12:52:47

#include "THDMIISNMSSMBCsimple_susy_parameters.hpp"
#include "config.h"
#ifdef ENABLE_THREADS
#include "global_thread_pool.hpp"
#endif
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define CLASSNAME THDMIISNMSSMBCsimple_susy_parameters
#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces
#define TRACE_STRUCT_TYPE Susy_traces
#define CALCULATE_TRACES(l) calc_susy_traces(l);

const int THDMIISNMSSMBCsimple_susy_parameters::numberOfParameters;

THDMIISNMSSMBCsimple_susy_parameters::THDMIISNMSSMBCsimple_susy_parameters(const THDMIISNMSSMBCsimple_input_parameters& input_)
   : input(input_)
{
   set_number_of_parameters(numberOfParameters);
}

THDMIISNMSSMBCsimple_susy_parameters::THDMIISNMSSMBCsimple_susy_parameters(
   double scale_, int loops_, int thresholds_,
   const THDMIISNMSSMBCsimple_input_parameters& input_
   , double g1_, double g2_, double g3_, double Lambda7_, double Lambda1_, double
   Lambda4_, double Lambda3_, double Lambda2_, double Lambda5_, double Lambda6_
   , double Lambda8_, const Eigen::Matrix<double,3,3>& Yd_, const Eigen::Matrix
   <double,3,3>& Ye_, const Eigen::Matrix<double,3,3>& Yu_
)
   : Beta_function()
   , g1(g1_), g2(g2_), g3(g3_), Lambda7(Lambda7_), Lambda1(Lambda1_), Lambda4(
   Lambda4_), Lambda3(Lambda3_), Lambda2(Lambda2_), Lambda5(Lambda5_), Lambda6(
   Lambda6_), Lambda8(Lambda8_), Yd(Yd_), Ye(Ye_), Yu(Yu_)
   , input(input_)
{
   set_number_of_parameters(numberOfParameters);
   set_scale(scale_);
   set_loops(loops_);
   set_thresholds(thresholds_);
}

Eigen::ArrayXd THDMIISNMSSMBCsimple_susy_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

THDMIISNMSSMBCsimple_susy_parameters THDMIISNMSSMBCsimple_susy_parameters::calc_beta(int loops) const
{
   double beta_g1 = 0.;
   double beta_g2 = 0.;
   double beta_g3 = 0.;
   double beta_Lambda7 = 0.;
   double beta_Lambda1 = 0.;
   double beta_Lambda4 = 0.;
   double beta_Lambda3 = 0.;
   double beta_Lambda2 = 0.;
   double beta_Lambda5 = 0.;
   double beta_Lambda6 = 0.;
   double beta_Lambda8 = 0.;
   Eigen::Matrix<double,3,3> beta_Yd = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_Ye = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_Yu = Eigen::Matrix<double,3,3>::Zero();

   if (loops > 0) {
      const auto TRACE_STRUCT = CALCULATE_TRACES(loops);

      beta_g1 += calc_beta_g1_1_loop(TRACE_STRUCT);
      beta_g2 += calc_beta_g2_1_loop(TRACE_STRUCT);
      beta_g3 += calc_beta_g3_1_loop(TRACE_STRUCT);
      beta_Lambda7 += calc_beta_Lambda7_1_loop(TRACE_STRUCT);
      beta_Lambda1 += calc_beta_Lambda1_1_loop(TRACE_STRUCT);
      beta_Lambda4 += calc_beta_Lambda4_1_loop(TRACE_STRUCT);
      beta_Lambda3 += calc_beta_Lambda3_1_loop(TRACE_STRUCT);
      beta_Lambda2 += calc_beta_Lambda2_1_loop(TRACE_STRUCT);
      beta_Lambda5 += calc_beta_Lambda5_1_loop(TRACE_STRUCT);
      beta_Lambda6 += calc_beta_Lambda6_1_loop(TRACE_STRUCT);
      beta_Lambda8 += calc_beta_Lambda8_1_loop(TRACE_STRUCT);
      beta_Yd += calc_beta_Yd_1_loop(TRACE_STRUCT);
      beta_Ye += calc_beta_Ye_1_loop(TRACE_STRUCT);
      beta_Yu += calc_beta_Yu_1_loop(TRACE_STRUCT);

      if (loops > 1) {
         beta_g1 += calc_beta_g1_2_loop(TRACE_STRUCT);
         beta_g2 += calc_beta_g2_2_loop(TRACE_STRUCT);
         beta_g3 += calc_beta_g3_2_loop(TRACE_STRUCT);
         beta_Lambda7 += calc_beta_Lambda7_2_loop(TRACE_STRUCT);
         beta_Lambda1 += calc_beta_Lambda1_2_loop(TRACE_STRUCT);
         beta_Lambda4 += calc_beta_Lambda4_2_loop(TRACE_STRUCT);
         beta_Lambda3 += calc_beta_Lambda3_2_loop(TRACE_STRUCT);
         beta_Lambda2 += calc_beta_Lambda2_2_loop(TRACE_STRUCT);
         beta_Lambda5 += calc_beta_Lambda5_2_loop(TRACE_STRUCT);
         beta_Lambda6 += calc_beta_Lambda6_2_loop(TRACE_STRUCT);
         beta_Lambda8 += calc_beta_Lambda8_2_loop(TRACE_STRUCT);
         beta_Yd += calc_beta_Yd_2_loop(TRACE_STRUCT);
         beta_Ye += calc_beta_Ye_2_loop(TRACE_STRUCT);
         beta_Yu += calc_beta_Yu_2_loop(TRACE_STRUCT);

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


   return THDMIISNMSSMBCsimple_susy_parameters(get_scale(), loops, get_thresholds(), input,
                    beta_g1, beta_g2, beta_g3, beta_Lambda7, beta_Lambda1, beta_Lambda4, beta_Lambda3, beta_Lambda2, beta_Lambda5, beta_Lambda6, beta_Lambda8, beta_Yd, beta_Ye, beta_Yu);
}

THDMIISNMSSMBCsimple_susy_parameters THDMIISNMSSMBCsimple_susy_parameters::calc_beta() const
{
   return calc_beta(get_loops());
}

void THDMIISNMSSMBCsimple_susy_parameters::clear()
{
   reset();
   g1 = 0.;
   g2 = 0.;
   g3 = 0.;
   Lambda7 = 0.;
   Lambda1 = 0.;
   Lambda4 = 0.;
   Lambda3 = 0.;
   Lambda2 = 0.;
   Lambda5 = 0.;
   Lambda6 = 0.;
   Lambda8 = 0.;
   Yd = Eigen::Matrix<double,3,3>::Zero();
   Ye = Eigen::Matrix<double,3,3>::Zero();
   Yu = Eigen::Matrix<double,3,3>::Zero();

}



Eigen::ArrayXd THDMIISNMSSMBCsimple_susy_parameters::get() const
{
   Eigen::ArrayXd pars(numberOfParameters);

   pars(0) = g1;
   pars(1) = g2;
   pars(2) = g3;
   pars(3) = Lambda7;
   pars(4) = Lambda1;
   pars(5) = Lambda4;
   pars(6) = Lambda3;
   pars(7) = Lambda2;
   pars(8) = Lambda5;
   pars(9) = Lambda6;
   pars(10) = Lambda8;
   pars(11) = Yd(0,0);
   pars(12) = Yd(0,1);
   pars(13) = Yd(0,2);
   pars(14) = Yd(1,0);
   pars(15) = Yd(1,1);
   pars(16) = Yd(1,2);
   pars(17) = Yd(2,0);
   pars(18) = Yd(2,1);
   pars(19) = Yd(2,2);
   pars(20) = Ye(0,0);
   pars(21) = Ye(0,1);
   pars(22) = Ye(0,2);
   pars(23) = Ye(1,0);
   pars(24) = Ye(1,1);
   pars(25) = Ye(1,2);
   pars(26) = Ye(2,0);
   pars(27) = Ye(2,1);
   pars(28) = Ye(2,2);
   pars(29) = Yu(0,0);
   pars(30) = Yu(0,1);
   pars(31) = Yu(0,2);
   pars(32) = Yu(1,0);
   pars(33) = Yu(1,1);
   pars(34) = Yu(1,2);
   pars(35) = Yu(2,0);
   pars(36) = Yu(2,1);
   pars(37) = Yu(2,2);


   return pars;
}

void THDMIISNMSSMBCsimple_susy_parameters::print(std::ostream& ostr) const
{
   ostr << "susy parameters at Q = " << get_scale() << ":\n";
   ostr << "g1 = " << g1 << '\n';
   ostr << "g2 = " << g2 << '\n';
   ostr << "g3 = " << g3 << '\n';
   ostr << "Lambda7 = " << Lambda7 << '\n';
   ostr << "Lambda1 = " << Lambda1 << '\n';
   ostr << "Lambda4 = " << Lambda4 << '\n';
   ostr << "Lambda3 = " << Lambda3 << '\n';
   ostr << "Lambda2 = " << Lambda2 << '\n';
   ostr << "Lambda5 = " << Lambda5 << '\n';
   ostr << "Lambda6 = " << Lambda6 << '\n';
   ostr << "Lambda8 = " << Lambda8 << '\n';
   ostr << "Yd = " << Yd << '\n';
   ostr << "Ye = " << Ye << '\n';
   ostr << "Yu = " << Yu << '\n';

}

void THDMIISNMSSMBCsimple_susy_parameters::set(const Eigen::ArrayXd& pars)
{
   g1 = pars(0);
   g2 = pars(1);
   g3 = pars(2);
   Lambda7 = pars(3);
   Lambda1 = pars(4);
   Lambda4 = pars(5);
   Lambda3 = pars(6);
   Lambda2 = pars(7);
   Lambda5 = pars(8);
   Lambda6 = pars(9);
   Lambda8 = pars(10);
   Yd(0,0) = pars(11);
   Yd(0,1) = pars(12);
   Yd(0,2) = pars(13);
   Yd(1,0) = pars(14);
   Yd(1,1) = pars(15);
   Yd(1,2) = pars(16);
   Yd(2,0) = pars(17);
   Yd(2,1) = pars(18);
   Yd(2,2) = pars(19);
   Ye(0,0) = pars(20);
   Ye(0,1) = pars(21);
   Ye(0,2) = pars(22);
   Ye(1,0) = pars(23);
   Ye(1,1) = pars(24);
   Ye(1,2) = pars(25);
   Ye(2,0) = pars(26);
   Ye(2,1) = pars(27);
   Ye(2,2) = pars(28);
   Yu(0,0) = pars(29);
   Yu(0,1) = pars(30);
   Yu(0,2) = pars(31);
   Yu(1,0) = pars(32);
   Yu(1,1) = pars(33);
   Yu(1,2) = pars(34);
   Yu(2,0) = pars(35);
   Yu(2,1) = pars(36);
   Yu(2,2) = pars(37);

}

const THDMIISNMSSMBCsimple_input_parameters& THDMIISNMSSMBCsimple_susy_parameters::get_input() const
{
   return input;
}

THDMIISNMSSMBCsimple_input_parameters& THDMIISNMSSMBCsimple_susy_parameters::get_input()
{
   return input;
}

void THDMIISNMSSMBCsimple_susy_parameters::set_input_parameters(const THDMIISNMSSMBCsimple_input_parameters& input_)
{
   input = input_;
}

THDMIISNMSSMBCsimple_susy_parameters::Susy_traces THDMIISNMSSMBCsimple_susy_parameters::calc_susy_traces(int loops) const
{
   Susy_traces susy_traces;

   if (loops > 0) {
      

      TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());
      TRACE_STRUCT.traceYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint()).trace()
         );
      TRACE_STRUCT.traceYdAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace()
         );
      TRACE_STRUCT.traceYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace()
         );
      TRACE_STRUCT.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint()).trace()
         );

   }

   if (loops > 1) {
      TRACE_STRUCT.traceYdAdjYdYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint()*
         Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYdAdjYdYdAdjYuYuAdjYd = Re((Yd*Yd.adjoint()*Yd*Yu.adjoint()*
         Yu*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYdAdjYuYuAdjYdYdAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint()*
         Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYdAdjYuYuAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yu.adjoint()*
         Yu*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYeYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint()*
         Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYuYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint()*
         Yu*Yu.adjoint()).trace());

   }

   if (loops > 2) {

   }

   return susy_traces;
}

std::ostream& operator<<(std::ostream& ostr, const THDMIISNMSSMBCsimple_susy_parameters& susy_pars)
{
   susy_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
