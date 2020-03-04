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

// File generated at Thu 7 Nov 2019 18:51:22

#ifndef THDMIISNMSSMBCsimple_susy_parameters_H
#define THDMIISNMSSMBCsimple_susy_parameters_H

#include "betafunction.hpp"
#include "THDMIISNMSSMBCsimple_input_parameters.hpp"

#include <iosfwd>
#include <string>
#include <vector>
#include <Eigen/Core>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Susy_traces

class THDMIISNMSSMBCsimple_susy_parameters : public Beta_function {
public:
   explicit THDMIISNMSSMBCsimple_susy_parameters(const THDMIISNMSSMBCsimple_input_parameters& input_ = THDMIISNMSSMBCsimple_input_parameters());
   THDMIISNMSSMBCsimple_susy_parameters(double scale_, int loops_, int thresholds_, const THDMIISNMSSMBCsimple_input_parameters& input_, double g1_, double g2_, double g3_, double Lambda7_, double Lambda1_, double
   Lambda4_, double Lambda3_, double Lambda2_, double Lambda5_, double Lambda6_
   , double Lambda8_, const Eigen::Matrix<double,3,3>& Yd_, const Eigen::Matrix
   <double,3,3>& Ye_, const Eigen::Matrix<double,3,3>& Yu_);
   THDMIISNMSSMBCsimple_susy_parameters(const THDMIISNMSSMBCsimple_susy_parameters&) = default;
   THDMIISNMSSMBCsimple_susy_parameters(THDMIISNMSSMBCsimple_susy_parameters&&) = default;
   virtual ~THDMIISNMSSMBCsimple_susy_parameters() = default;
   THDMIISNMSSMBCsimple_susy_parameters& operator=(const THDMIISNMSSMBCsimple_susy_parameters&) = default;
   THDMIISNMSSMBCsimple_susy_parameters& operator=(THDMIISNMSSMBCsimple_susy_parameters&&) = default;

   virtual Eigen::ArrayXd beta() const override;
   virtual Eigen::ArrayXd get() const override;
   virtual void print(std::ostream&) const;
   virtual void set(const Eigen::ArrayXd&) override;
   const THDMIISNMSSMBCsimple_input_parameters& get_input() const;
   THDMIISNMSSMBCsimple_input_parameters& get_input();
   void set_input_parameters(const THDMIISNMSSMBCsimple_input_parameters&);

   THDMIISNMSSMBCsimple_susy_parameters calc_beta() const;
   THDMIISNMSSMBCsimple_susy_parameters calc_beta(int) const;
   virtual void clear();

   void set_g1(double g1_) { g1 = g1_; }
   void set_g2(double g2_) { g2 = g2_; }
   void set_g3(double g3_) { g3 = g3_; }
   void set_Lambda7(double Lambda7_) { Lambda7 = Lambda7_; }
   void set_Lambda1(double Lambda1_) { Lambda1 = Lambda1_; }
   void set_Lambda4(double Lambda4_) { Lambda4 = Lambda4_; }
   void set_Lambda3(double Lambda3_) { Lambda3 = Lambda3_; }
   void set_Lambda2(double Lambda2_) { Lambda2 = Lambda2_; }
   void set_Lambda5(double Lambda5_) { Lambda5 = Lambda5_; }
   void set_Lambda6(double Lambda6_) { Lambda6 = Lambda6_; }
   void set_Lambda8(double Lambda8_) { Lambda8 = Lambda8_; }
   void set_Yd(const Eigen::Matrix<double,3,3>& Yd_) { Yd = Yd_; }
   void set_Yd(int i, int k, const double& value) { Yd(i,k) = value; }
   void set_Ye(const Eigen::Matrix<double,3,3>& Ye_) { Ye = Ye_; }
   void set_Ye(int i, int k, const double& value) { Ye(i,k) = value; }
   void set_Yu(const Eigen::Matrix<double,3,3>& Yu_) { Yu = Yu_; }
   void set_Yu(int i, int k, const double& value) { Yu(i,k) = value; }

   double get_g1() const { return g1; }
   double get_g2() const { return g2; }
   double get_g3() const { return g3; }
   double get_Lambda7() const { return Lambda7; }
   double get_Lambda1() const { return Lambda1; }
   double get_Lambda4() const { return Lambda4; }
   double get_Lambda3() const { return Lambda3; }
   double get_Lambda2() const { return Lambda2; }
   double get_Lambda5() const { return Lambda5; }
   double get_Lambda6() const { return Lambda6; }
   double get_Lambda8() const { return Lambda8; }
   const Eigen::Matrix<double,3,3>& get_Yd() const { return Yd; }
   double get_Yd(int i, int k) const { return Yd(i,k); }
   const Eigen::Matrix<double,3,3>& get_Ye() const { return Ye; }
   double get_Ye(int i, int k) const { return Ye(i,k); }
   const Eigen::Matrix<double,3,3>& get_Yu() const { return Yu; }
   double get_Yu(int i, int k) const { return Yu(i,k); }



protected:
   double g1{};
   double g2{};
   double g3{};
   double Lambda7{};
   double Lambda1{};
   double Lambda4{};
   double Lambda3{};
   double Lambda2{};
   double Lambda5{};
   double Lambda6{};
   double Lambda8{};
   Eigen::Matrix<double,3,3> Yd{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> Ye{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> Yu{Eigen::Matrix<double,3,3>::Zero()};

   THDMIISNMSSMBCsimple_input_parameters input{};

private:
   static const int numberOfParameters = 38;

   struct Susy_traces {
      double traceYdAdjYd{};
      double traceYeAdjYe{};
      double traceYuAdjYu{};
      double traceYdAdjYdYdAdjYd{};
      double traceYdAdjYuYuAdjYd{};
      double traceYeAdjYeYeAdjYe{};
      double traceYuAdjYuYuAdjYu{};
      double traceYdAdjYdYdAdjYdYdAdjYd{};
      double traceYdAdjYuYuAdjYdYdAdjYd{};
      double traceYeAdjYeYeAdjYeYeAdjYe{};
      double traceYdAdjYdYdAdjYuYuAdjYd{};
      double traceYdAdjYuYuAdjYuYuAdjYd{};
      double traceYuAdjYuYuAdjYuYuAdjYu{};

   };
   Susy_traces calc_susy_traces(int) const;

   double calc_beta_g1_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g1_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g1_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g1_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g1_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g3_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g3_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g3_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g3_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g3_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda7_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda7_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda7_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda7_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda7_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda1_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda1_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda1_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda1_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda1_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda4_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda4_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda4_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda4_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda4_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda3_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda3_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda3_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda3_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda3_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda2_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda2_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda5_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda5_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda5_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda5_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda5_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda6_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda6_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda6_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda6_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda6_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda8_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda8_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda8_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda8_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda8_5_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yd_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yd_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yd_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yd_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yd_5_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_5_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_5_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const THDMIISNMSSMBCsimple_susy_parameters&);

#undef TRACE_STRUCT_TYPE

} // namespace flexiblesusy

#endif
