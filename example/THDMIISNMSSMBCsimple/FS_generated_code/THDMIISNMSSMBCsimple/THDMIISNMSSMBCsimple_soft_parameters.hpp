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

// File generated at Thu 7 Nov 2019 18:51:36

#ifndef THDMIISNMSSMBCsimple_soft_parameters_H
#define THDMIISNMSSMBCsimple_soft_parameters_H

#include "THDMIISNMSSMBCsimple_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class THDMIISNMSSMBCsimple_soft_parameters : public THDMIISNMSSMBCsimple_susy_parameters {
public:
   explicit THDMIISNMSSMBCsimple_soft_parameters(const THDMIISNMSSMBCsimple_input_parameters& input_ = THDMIISNMSSMBCsimple_input_parameters());
   THDMIISNMSSMBCsimple_soft_parameters(const THDMIISNMSSMBCsimple_susy_parameters& , double M123_, double M5_, double M112_, double M222_, double M332_, double
   v1_, double v2_, double vS_);
   THDMIISNMSSMBCsimple_soft_parameters(const THDMIISNMSSMBCsimple_soft_parameters&) = default;
   THDMIISNMSSMBCsimple_soft_parameters(THDMIISNMSSMBCsimple_soft_parameters&&) = default;
   virtual ~THDMIISNMSSMBCsimple_soft_parameters() = default;
   THDMIISNMSSMBCsimple_soft_parameters& operator=(const THDMIISNMSSMBCsimple_soft_parameters&) = default;
   THDMIISNMSSMBCsimple_soft_parameters& operator=(THDMIISNMSSMBCsimple_soft_parameters&&) = default;

   virtual Eigen::ArrayXd beta() const override;
   virtual Eigen::ArrayXd get() const override;
   virtual void print(std::ostream&) const override;
   virtual void set(const Eigen::ArrayXd&) override;

   THDMIISNMSSMBCsimple_soft_parameters calc_beta() const;
   THDMIISNMSSMBCsimple_soft_parameters calc_beta(int) const;
   virtual void clear() override;

   void set_M123(double M123_) { M123 = M123_; }
   void set_M5(double M5_) { M5 = M5_; }
   void set_M112(double M112_) { M112 = M112_; }
   void set_M222(double M222_) { M222 = M222_; }
   void set_M332(double M332_) { M332 = M332_; }
   void set_v1(double v1_) { v1 = v1_; }
   void set_v2(double v2_) { v2 = v2_; }
   void set_vS(double vS_) { vS = vS_; }

   double get_M123() const { return M123; }
   double get_M5() const { return M5; }
   double get_M112() const { return M112; }
   double get_M222() const { return M222; }
   double get_M332() const { return M332; }
   double get_v1() const { return v1; }
   double get_v2() const { return v2; }
   double get_vS() const { return vS; }


protected:
   double M123{};
   double M5{};
   double M112{};
   double M222{};
   double M332{};
   double v1{};
   double v2{};
   double vS{};


private:
   static const int numberOfParameters = 46;

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

   double calc_beta_M123_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M123_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M123_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M123_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M123_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M5_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M5_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M5_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M5_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M5_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M112_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M112_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M112_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M112_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M112_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M222_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M222_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M222_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M222_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M222_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M332_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M332_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M332_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M332_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M332_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v1_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v1_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v1_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v1_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v1_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v2_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v2_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vS_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vS_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vS_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vS_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vS_5_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const THDMIISNMSSMBCsimple_soft_parameters&);

} // namespace flexiblesusy

#undef TRACE_STRUCT_TYPE

#endif
