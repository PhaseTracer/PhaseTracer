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

// File generated at Tue 18 Aug 2020 22:46:29

#ifndef SingletDM_soft_parameters_H
#define SingletDM_soft_parameters_H

#include "SingletDM_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class SingletDM_soft_parameters : public SingletDM_susy_parameters {
public:
   explicit SingletDM_soft_parameters(const SingletDM_input_parameters& input_ = SingletDM_input_parameters());
   SingletDM_soft_parameters(const SingletDM_susy_parameters& , double muS_, double muH_, double v_);
   SingletDM_soft_parameters(const SingletDM_soft_parameters&) = default;
   SingletDM_soft_parameters(SingletDM_soft_parameters&&) = default;
   virtual ~SingletDM_soft_parameters() = default;
   SingletDM_soft_parameters& operator=(const SingletDM_soft_parameters&) = default;
   SingletDM_soft_parameters& operator=(SingletDM_soft_parameters&&) = default;

   virtual Eigen::ArrayXd beta() const override;
   virtual Eigen::ArrayXd get() const override;
   virtual void print(std::ostream&) const override;
   virtual void set(const Eigen::ArrayXd&) override;

   SingletDM_soft_parameters calc_beta() const;
   SingletDM_soft_parameters calc_beta(int) const;
   virtual void clear() override;

   void set_muS(double muS_) { muS = muS_; }
   void set_muH(double muH_) { muH = muH_; }
   void set_v(double v_) { v = v_; }

   double get_muS() const { return muS; }
   double get_muH() const { return muH; }
   double get_v() const { return v; }


protected:
   double muS{};
   double muH{};
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

   double calc_beta_muS_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_muS_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_muS_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_muS_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_muS_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_muH_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_muH_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_muH_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_muH_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_muH_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_5_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const SingletDM_soft_parameters&);

} // namespace flexiblesusy

#undef TRACE_STRUCT_TYPE

#endif
