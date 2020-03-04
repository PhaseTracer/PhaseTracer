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

#ifndef THDMIISNMSSMBCsimple_TWO_SCALE_LOW_SCALE_CONSTRAINT_H
#define THDMIISNMSSMBCsimple_TWO_SCALE_LOW_SCALE_CONSTRAINT_H

#include "THDMIISNMSSMBCsimple_low_scale_constraint.hpp"
#include "THDMIISNMSSMBCsimple_input_parameters.hpp"
#include "single_scale_constraint.hpp"
#include "lowe.h"
#include <Eigen/Core>

namespace flexiblesusy {

template <class T>
class THDMIISNMSSMBCsimple;

class Two_scale;

template<>
class THDMIISNMSSMBCsimple_low_scale_constraint<Two_scale> : public Single_scale_constraint {
public:
   THDMIISNMSSMBCsimple_low_scale_constraint() = default;
   THDMIISNMSSMBCsimple_low_scale_constraint(THDMIISNMSSMBCsimple<Two_scale>*, const softsusy::QedQcd&);
   virtual ~THDMIISNMSSMBCsimple_low_scale_constraint() = default;
   virtual void apply() override;
   virtual double get_scale() const override;
   virtual std::string name() const override { return "THDMIISNMSSMBCsimple low-scale constraint"; }
   virtual void set_model(Model*) override;

   void clear();
   const Eigen::Matrix<std::complex<double>,3,3>& get_ckm();
   const Eigen::Matrix<std::complex<double>,3,3>& get_pmns();
   double get_initial_scale_guess() const;
   void initialize();
   const softsusy::QedQcd& get_sm_parameters() const;
   void set_sm_parameters(const softsusy::QedQcd&);

private:
   double scale{0.};
   double initial_scale_guess{0.};
   THDMIISNMSSMBCsimple<Two_scale>* model{nullptr};
   softsusy::QedQcd qedqcd{};
   Eigen::Matrix<std::complex<double>,3,3> ckm{Eigen::Matrix<std::complex<double>,3,3>::Identity()};
   Eigen::Matrix<std::complex<double>,3,3> pmns{Eigen::Matrix<std::complex<double>,3,3>::Identity()};
   Eigen::Matrix<std::complex<double>,3,3> upQuarksDRbar{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> downQuarksDRbar{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> downLeptonsDRbar{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<double,3,3> neutrinoDRbar{Eigen::Matrix<double,3,3>::Zero()};
   double mW_run{0.};
   double mZ_run{0.};
   double AlphaS{0.};
   double e_run{0.};
   double ThetaWDRbar{0.};
   double new_g1{0.}, new_g2{0.}, new_g3{0.};

   double calculate_theta_w();
   void calculate_threshold_corrections();
   void calculate_DRbar_gauge_couplings();
   void calculate_DRbar_yukawa_couplings();
   void calculate_Yu_DRbar();
   void calculate_Yd_DRbar();
   void calculate_Ye_DRbar();
   void calculate_running_SM_masses();
   double calculate_delta_alpha_em(double) const;
   double calculate_delta_alpha_s(double) const;
   double calculate_alpha_s_SM5_at(softsusy::QedQcd, double) const;
   void check_model_ptr() const;
   void update_scale();
};

} // namespace flexiblesusy

#endif
