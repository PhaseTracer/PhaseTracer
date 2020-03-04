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

#include "THDMIISNMSSMBCsimple_two_scale_high_scale_constraint.hpp"
#include "THDMIISNMSSMBCsimple_two_scale_model.hpp"
#include "THDMIISNMSSMBCsimple_info.hpp"
#include "config.h"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "error.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "raii.hpp"
#include "root_finder.hpp"
#include "threshold_loop_functions.hpp"
#include "numerics2.hpp"

#ifdef ENABLE_HIMALAYA
#include "HierarchyCalculator.hpp"
#include "version.hpp"
#endif

#include <cmath>
#include <cerrno>
#include <cstring>

namespace flexiblesusy {

#define DERIVEDPARAMETER(p) model->p()
#define EXTRAPARAMETER(p) model->get_##p()
#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define PHASE(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETAPARAMETER1(l,p) beta_functions_##l##L.get_##p()
#define BETA(p) beta_##p
#define BETA1(l,p) beta_##l##L_##p
#define LowEnergyConstant(p) Electroweak_constants::p
#define MZPole Electroweak_constants::MZ
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define SCALE model->get_scale()
#define THRESHOLD static_cast<int>(model->get_thresholds())
#define MODEL model
#define MODELCLASSNAME THDMIISNMSSMBCsimple<Two_scale>

#if defined(ENABLE_HIMALAYA) && Himalaya_VERSION_MAJOR >= 2
#define FSHimalayaMh23L [&] () {                                        \
      MODEL->calculate_DRbar_masses();                                  \
                                                                        \
      himalaya::Parameters pars;                                        \
      const auto g1 = MODELPARAMETER(g1); \
      const auto g2 = MODELPARAMETER(g2); \
      const auto g3 = MODELPARAMETER(g3); \
      const auto v2 = MODELPARAMETER(v2); \
      const auto v1 = MODELPARAMETER(v1); \
      const auto Yu = MODELPARAMETER(Yu); \
      const auto Yd = MODELPARAMETER(Yd); \
      const auto Ye = MODELPARAMETER(Ye); \
      const auto Vu = MODELPARAMETER(Vu); \
      const auto Ud = MODELPARAMETER(Ud); \
      const auto Uu = MODELPARAMETER(Uu); \
      const auto Ve = MODELPARAMETER(Ve); \
      const auto Ue = MODELPARAMETER(Ue); \
      const auto MAh = MODELPARAMETER(MAh); \
       \
       \
      pars.scale = MODELPARAMETER(scale); \
      pars.mu = Re(Mu); \
      pars.g1 = Re(g1); \
      pars.g2 = Re(g2); \
      pars.g3 = Re(g3); \
      pars.vd = Re(v1); \
      pars.vu = Re(v2); \
      pars.mq2 = Re(Vu); \
      pars.md2 = Re(Ud); \
      pars.mu2 = Re(Uu); \
      pars.ml2 = Re(Ve); \
      pars.me2 = Re(Ue); \
      pars.Au(2,2) = Re(TrilinearUp); \
      pars.Ad(2,2) = Re(TrilinearDown); \
      pars.Ae(2,2) = Re(TrilinearLepton); \
      pars.Yu = Re(Yu); \
      pars.Yd = Re(Yd); \
      pars.Ye = Re(Ye); \
      pars.M1 = 0; \
      pars.M2 = 0; \
      pars.MG = MGluino; \
      pars.MA = MAh; \
       \
      const double msbar_scheme = 1; \
      const double lambda_3L_eft = 1; \
      const double lambda_3L_uncertainty = 0; \
       \
                                                                        \
      double lambda_3L = 0.;                                            \
                                                                        \
      try {                                                             \
         const bool verbose = false;                                    \
         himalaya::HierarchyCalculator hc(pars, verbose);               \
                                                                        \
         const auto ho = hc.calculateDMh3L(false);                      \
                                                                        \
         lambda_3L =                                                    \
            lambda_3L_eft * (                                           \
               ho.getDLambda(3)                                         \
               + msbar_scheme*ho.getDLambdaDRbarPrimeToMSbarShift(3)    \
               + lambda_3L_uncertainty*ho.getDLambdaUncertainty(3)      \
            );                                                          \
                                                                        \
         VERBOSE_MSG("Himalaya top (hierarchy, Dlambda_3L) = ("         \
                     << ho.getSuitableHierarchy() << ", "               \
                     << lambda_3L <<")");                               \
      } catch (const std::exception& e) {                               \
         model->get_problems().flag_bad_mass(THDMIISNMSSMBCsimple_info::hh); \
         WARNING(e.what());                                             \
         VERBOSE_MSG(pars);                                             \
      }                                                                 \
                                                                        \
      return lambda_3L;                                                 \
   }()
#else
#define FSHimalayaMh23L [] () {                                         \
      throw HimalayaError("The 3-loop corrections to lambda "           \
                          "require Himalaya 2.0.0 (or higher), but "    \
                          "FlexibleSUSY has not been configured with "  \
                          "this Himalaya version!");                    \
      return 0.;                                                        \
   }()
#endif

THDMIISNMSSMBCsimple_high_scale_constraint<Two_scale>::THDMIISNMSSMBCsimple_high_scale_constraint(
   THDMIISNMSSMBCsimple<Two_scale>* model_)
   : model(model_)
{
   initialize();
}

void THDMIISNMSSMBCsimple_high_scale_constraint<Two_scale>::apply()
{
   check_model_ptr();

   


   update_scale();

   const auto AtopNMSSM = INPUTPARAMETER(AtopNMSSM);
   const auto mstopL = INPUTPARAMETER(mstopL);
   const auto mstopR = INPUTPARAMETER(mstopR);
   const auto lambdaNMSSM = INPUTPARAMETER(lambdaNMSSM);
   const auto kappaNMSSM = INPUTPARAMETER(kappaNMSSM);
   const auto vSIN = INPUTPARAMETER(vSIN);
   const auto AlambdaNMSSM = INPUTPARAMETER(AlambdaNMSSM);
   const auto AkappaNMSSM = INPUTPARAMETER(AkappaNMSSM);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);

   MODEL->set_Lambda1(Re(0.25*(0.6*Sqr(g1) + Sqr(g2))));
   MODEL->set_Lambda2(Re(0.25*(0.6*Sqr(g1) + Sqr(g2)) + (0.037995443865876666*Sqr(
      AtopNMSSM)*(1 - (0.08333333333333333*Sqr(AtopNMSSM))/(mstopL*mstopR))*Sqr(Yu
      (2,2)))/(mstopL*mstopR)));
   MODEL->set_Lambda3(Re(-0.5*Sqr(g2) + 0.25*(-0.6*Sqr(g1) + Sqr(g2)) + Sqr(
      lambdaNMSSM)));
   MODEL->set_Lambda4(Re(-0.5*Sqr(g2) + Sqr(lambdaNMSSM)));
   MODEL->set_Lambda5(Re(Sqr(lambdaNMSSM)));
   MODEL->set_Lambda6(Re(Sqr(lambdaNMSSM)));
   MODEL->set_Lambda7(Re(-(kappaNMSSM*lambdaNMSSM)));
   MODEL->set_Lambda8(Re(Sqr(kappaNMSSM)));
   MODEL->set_vS(Re(vSIN));
   MODEL->set_M123(Re(AlambdaNMSSM*lambdaNMSSM));
   MODEL->set_M5(Re(-(AkappaNMSSM*kappaNMSSM)));


   check_non_perturbative();
}

bool THDMIISNMSSMBCsimple_high_scale_constraint<Two_scale>::check_non_perturbative()
{
   bool problem = false;

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto Lambda1 = MODELPARAMETER(Lambda1);
   const auto Lambda4 = MODELPARAMETER(Lambda4);
   const auto Lambda3 = MODELPARAMETER(Lambda3);
   const auto Lambda2 = MODELPARAMETER(Lambda2);
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto Lambda8 = MODELPARAMETER(Lambda8);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yu = MODELPARAMETER(Yu);

   if (MaxAbsValue(g1) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::g1, MaxAbsValue(g1), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::g1);
   }
   if (MaxAbsValue(g2) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::g2, MaxAbsValue(g2), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::g2);
   }
   if (MaxAbsValue(g3) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::g3, MaxAbsValue(g3), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::g3);
   }
   if (MaxAbsValue(Lambda7) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Lambda7, MaxAbsValue(Lambda7), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Lambda7);
   }
   if (MaxAbsValue(Lambda1) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Lambda1, MaxAbsValue(Lambda1), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Lambda1);
   }
   if (MaxAbsValue(Lambda4) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Lambda4, MaxAbsValue(Lambda4), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Lambda4);
   }
   if (MaxAbsValue(Lambda3) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Lambda3, MaxAbsValue(Lambda3), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Lambda3);
   }
   if (MaxAbsValue(Lambda2) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Lambda2, MaxAbsValue(Lambda2), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Lambda2);
   }
   if (MaxAbsValue(Lambda5) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Lambda5, MaxAbsValue(Lambda5), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Lambda5);
   }
   if (MaxAbsValue(Lambda6) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Lambda6, MaxAbsValue(Lambda6), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Lambda6);
   }
   if (MaxAbsValue(Lambda8) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Lambda8, MaxAbsValue(Lambda8), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Lambda8);
   }
   if (MaxAbsValue(Yd(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yd0_0, MaxAbsValue(Yd(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yd0_0);
   }

   if (MaxAbsValue(Yd(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yd0_1, MaxAbsValue(Yd(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yd0_1);
   }

   if (MaxAbsValue(Yd(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yd0_2, MaxAbsValue(Yd(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yd0_2);
   }

   if (MaxAbsValue(Yd(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yd1_0, MaxAbsValue(Yd(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yd1_0);
   }

   if (MaxAbsValue(Yd(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yd1_1, MaxAbsValue(Yd(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yd1_1);
   }

   if (MaxAbsValue(Yd(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yd1_2, MaxAbsValue(Yd(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yd1_2);
   }

   if (MaxAbsValue(Yd(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yd2_0, MaxAbsValue(Yd(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yd2_0);
   }

   if (MaxAbsValue(Yd(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yd2_1, MaxAbsValue(Yd(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yd2_1);
   }

   if (MaxAbsValue(Yd(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yd2_2, MaxAbsValue(Yd(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yd2_2);
   }
   if (MaxAbsValue(Ye(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Ye0_0, MaxAbsValue(Ye(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Ye0_0);
   }

   if (MaxAbsValue(Ye(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Ye0_1, MaxAbsValue(Ye(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Ye0_1);
   }

   if (MaxAbsValue(Ye(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Ye0_2, MaxAbsValue(Ye(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Ye0_2);
   }

   if (MaxAbsValue(Ye(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Ye1_0, MaxAbsValue(Ye(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Ye1_0);
   }

   if (MaxAbsValue(Ye(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Ye1_1, MaxAbsValue(Ye(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Ye1_1);
   }

   if (MaxAbsValue(Ye(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Ye1_2, MaxAbsValue(Ye(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Ye1_2);
   }

   if (MaxAbsValue(Ye(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Ye2_0, MaxAbsValue(Ye(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Ye2_0);
   }

   if (MaxAbsValue(Ye(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Ye2_1, MaxAbsValue(Ye(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Ye2_1);
   }

   if (MaxAbsValue(Ye(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Ye2_2, MaxAbsValue(Ye(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Ye2_2);
   }
   if (MaxAbsValue(Yu(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yu0_0, MaxAbsValue(Yu(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yu0_0);
   }

   if (MaxAbsValue(Yu(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yu0_1, MaxAbsValue(Yu(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yu0_1);
   }

   if (MaxAbsValue(Yu(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yu0_2, MaxAbsValue(Yu(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yu0_2);
   }

   if (MaxAbsValue(Yu(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yu1_0, MaxAbsValue(Yu(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yu1_0);
   }

   if (MaxAbsValue(Yu(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yu1_1, MaxAbsValue(Yu(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yu1_1);
   }

   if (MaxAbsValue(Yu(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yu1_2, MaxAbsValue(Yu(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yu1_2);
   }

   if (MaxAbsValue(Yu(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yu2_0, MaxAbsValue(Yu(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yu2_0);
   }

   if (MaxAbsValue(Yu(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yu2_1, MaxAbsValue(Yu(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yu2_1);
   }

   if (MaxAbsValue(Yu(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yu2_2, MaxAbsValue(Yu(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIISNMSSMBCsimple_info::Yu2_2);
   }


   return problem;
}

double THDMIISNMSSMBCsimple_high_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double THDMIISNMSSMBCsimple_high_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const THDMIISNMSSMBCsimple_input_parameters& THDMIISNMSSMBCsimple_high_scale_constraint<Two_scale>::get_input_parameters() const
{
   return model->get_input();
}

THDMIISNMSSMBCsimple<Two_scale>* THDMIISNMSSMBCsimple_high_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void THDMIISNMSSMBCsimple_high_scale_constraint<Two_scale>::set_model(Model* model_)
{
   model = cast_model<THDMIISNMSSMBCsimple<Two_scale>*>(model_);
}

void THDMIISNMSSMBCsimple_high_scale_constraint<Two_scale>::set_scale(double s)
{
   scale = s;
}

void THDMIISNMSSMBCsimple_high_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = nullptr;
}

void THDMIISNMSSMBCsimple_high_scale_constraint<Two_scale>::initialize()
{
   check_model_ptr();

   const auto mstopL = INPUTPARAMETER(mstopL);
   const auto mstopR = INPUTPARAMETER(mstopR);

   initial_scale_guess = Sqrt(mstopL*mstopR);

   scale = initial_scale_guess;
}

void THDMIISNMSSMBCsimple_high_scale_constraint<Two_scale>::update_scale()
{
   check_model_ptr();

   const auto mstopL = INPUTPARAMETER(mstopL);
   const auto mstopR = INPUTPARAMETER(mstopR);

   scale = Sqrt(mstopL*mstopR);


}

void THDMIISNMSSMBCsimple_high_scale_constraint<Two_scale>::check_model_ptr() const
{
   if (!model)
      throw SetupError("THDMIISNMSSMBCsimple_high_scale_constraint<Two_scale>: "
                       "model pointer is zero!");
}

} // namespace flexiblesusy
