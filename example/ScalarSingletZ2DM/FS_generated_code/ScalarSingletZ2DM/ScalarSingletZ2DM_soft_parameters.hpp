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

// File generated at Tue 17 Nov 2020 16:11:20

#ifndef ScalarSingletZ2DM_soft_parameters_H
#define ScalarSingletZ2DM_soft_parameters_H

#include "ScalarSingletZ2DM_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class ScalarSingletZ2DM_soft_parameters : public ScalarSingletZ2DM_susy_parameters {
public:
   explicit ScalarSingletZ2DM_soft_parameters(const ScalarSingletZ2DM_input_parameters& input_ = ScalarSingletZ2DM_input_parameters());
   ScalarSingletZ2DM_soft_parameters(const ScalarSingletZ2DM_susy_parameters& , double muS2_, double muH2_, double v_);
   ScalarSingletZ2DM_soft_parameters(const ScalarSingletZ2DM_soft_parameters&) = default;
   ScalarSingletZ2DM_soft_parameters(ScalarSingletZ2DM_soft_parameters&&) = default;
   virtual ~ScalarSingletZ2DM_soft_parameters() = default;
   ScalarSingletZ2DM_soft_parameters& operator=(const ScalarSingletZ2DM_soft_parameters&) = default;
   ScalarSingletZ2DM_soft_parameters& operator=(ScalarSingletZ2DM_soft_parameters&&) = default;

   virtual Eigen::ArrayXd beta() const override;
   virtual Eigen::ArrayXd get() const override;
   virtual void print(std::ostream&) const override;
   virtual void set(const Eigen::ArrayXd&) override;

   ScalarSingletZ2DM_soft_parameters calc_beta() const;
   ScalarSingletZ2DM_soft_parameters calc_beta(int) const;
   virtual void clear() override;

   void set_muS2(double muS2_) { muS2 = muS2_; }
   void set_muH2(double muH2_) { muH2 = muH2_; }
   void set_v(double v_) { v = v_; }

   double get_muS2() const { return muS2; }
   double get_muH2() const { return muH2; }
   double get_v() const { return v; }


protected:
   double muS2{};
   double muH2{};
   double v{};


private:
   static const int numberOfParameters = 36;

   struct Soft_traces {
      double traceYdAdjYd{};
      double traceYeAdjYe{};
      double traceYuAdjYu{};
      double traceYdAdjYdYdAdjYd{};
      double traceYdAdjYuYuAdjYd{};
      double traceYeAdjYeYeAdjYe{};
      double traceYuAdjYuYuAdjYu{};

   };
   Soft_traces calc_soft_traces(int) const;

   double calc_beta_muS2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_muS2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_muS2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_muS2_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_muS2_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_muH2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_muH2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_muH2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_muH2_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_muH2_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_5_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const ScalarSingletZ2DM_soft_parameters&);

} // namespace flexiblesusy

#undef TRACE_STRUCT_TYPE

#endif
