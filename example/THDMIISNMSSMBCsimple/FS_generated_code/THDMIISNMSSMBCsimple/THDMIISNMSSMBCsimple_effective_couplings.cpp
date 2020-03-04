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

// File generated at Thu 7 Nov 2019 18:53:57

#include "THDMIISNMSSMBCsimple_effective_couplings.hpp"

#include "effective_couplings.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

using namespace effective_couplings;

#define INPUTPARAMETER(parameter) model.get_input().parameter
#define MODELPARAMETER(parameter) model.get_##parameter()
#define DERIVEDPARAMETER(parameter) model.##parameter()
#define PHASE(parameter) model.get_##parameter()
#define PHYSICAL(parameter) model.get_physical().parameter

THDMIISNMSSMBCsimple_effective_couplings::THDMIISNMSSMBCsimple_effective_couplings(
   const THDMIISNMSSMBCsimple_mass_eigenstates& model_,
   const softsusy::QedQcd& qedqcd_,
   const Physical_input& input_)
   : model(model_), qedqcd(qedqcd_), physical_input(input_)
   , rg_improve(true), include_qcd_corrections(true)
   , ZH(PHYSICAL(ZH)), ZA(PHYSICAL(ZA)), ZP(PHYSICAL(ZP)), Vd(PHYSICAL(Vd)), Ud(
      PHYSICAL(Ud)), Vu(PHYSICAL(Vu)), Uu(PHYSICAL(Uu)), Ve(PHYSICAL(Ve)), Ue(
      PHYSICAL(Ue)), ZZ(PHYSICAL(ZZ))
   , eff_CphhVPVP(Eigen::Array<std::complex<double>,3,1>::Zero()), eff_CphhVGVG(
      Eigen::Array<std::complex<double>,3,1>::Zero()), eff_CpAhVPVP(Eigen::Array<
      std::complex<double>,3,1>::Zero()), eff_CpAhVGVG(Eigen::Array<std::complex<
      double>,3,1>::Zero())
{
}

void THDMIISNMSSMBCsimple_effective_couplings::calculate_effective_couplings()
{
   const standard_model::Standard_model sm(initialise_SM());

   const double scale = model.get_scale();
   const Eigen::ArrayXd saved_parameters(model.get());

   const double saved_mt = PHYSICAL(MFu(2));
   PHYSICAL(MFu(2)) = qedqcd.displayPoleMt();

   const auto Mhh = PHYSICAL(Mhh);
   for (int gO1 = 0; gO1 < 3; ++gO1) {
      run_SM_strong_coupling_to(sm, 0.5 * Mhh(gO1));
      calculate_eff_CphhVPVP(gO1);
      run_SM_strong_coupling_to(sm, Mhh(gO1));
      calculate_eff_CphhVGVG(gO1);
   }

   const auto MAh = PHYSICAL(MAh);
   for (int gO1 = 1; gO1 < 3; ++gO1) {
      run_SM_strong_coupling_to(sm, 0.5 * MAh(gO1));
      calculate_eff_CpAhVPVP(gO1);
      run_SM_strong_coupling_to(sm, MAh(gO1));
      calculate_eff_CpAhVGVG(gO1);
   }

   PHYSICAL(MFu(2)) = saved_mt;
   model.set_scale(scale);
   model.set(saved_parameters);

}

void THDMIISNMSSMBCsimple_effective_couplings::set_model(const THDMIISNMSSMBCsimple_mass_eigenstates& model_)
{
   model = model_;
   copy_mixing_matrices_from_model();
}

void THDMIISNMSSMBCsimple_effective_couplings::copy_mixing_matrices_from_model()
{
   ZH = PHYSICAL(ZH);
   ZA = PHYSICAL(ZA);
   ZP = PHYSICAL(ZP);
   Vd = PHYSICAL(Vd);
   Ud = PHYSICAL(Ud);
   Vu = PHYSICAL(Vu);
   Uu = PHYSICAL(Uu);
   Ve = PHYSICAL(Ve);
   Ue = PHYSICAL(Ue);
   ZZ = PHYSICAL(ZZ);

}

standard_model::Standard_model THDMIISNMSSMBCsimple_effective_couplings::initialise_SM() const
{
   standard_model::Standard_model sm;

   sm.set_loops(2);
   sm.set_thresholds(2);
   sm.set_physical_input(physical_input);

   sm.initialise_from_input(qedqcd);

   return sm;
}

void THDMIISNMSSMBCsimple_effective_couplings::run_SM_strong_coupling_to(standard_model::Standard_model sm, double m)
{
   sm.run_to(m);

   model.set_g3(sm.get_g3());

}

std::complex<double> THDMIISNMSSMBCsimple_effective_couplings::scalar_scalar_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      if (m_loop > m_decay) {
         result = 1 + 0.06754745576155852*Sqr(g3);
      }
   }

   return result;
}

std::complex<double> THDMIISNMSSMBCsimple_effective_couplings::scalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) * scalar_diphoton_fermion_loop(
         m_decay, m_loop);

   }

   return result;
}

std::complex<double> THDMIISNMSSMBCsimple_effective_couplings::pseudoscalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         pseudoscalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

double THDMIISNMSSMBCsimple_effective_couplings::number_of_active_flavours(double m) const
{
   if (m < qedqcd.displayMbMb()) {
      return 4.0;
   } else if (m < qedqcd.displayPoleMt()) {
      return 5.0;
   } else {
      return 6.0;
   }
}

double THDMIISNMSSMBCsimple_effective_couplings::scalar_scaling_factor(double m) const
{
   const double Nf = number_of_active_flavours(m);
   const double mtpole = qedqcd.displayPoleMt();
   const double l = Log(Sqr(m) / Sqr(mtpole));

   const auto g3 = MODELPARAMETER(g3);

   const double nlo_qcd = 0.025330295910584444*(23.75 - 1.1666666666666667*Nf)*Sqr
      (g3);
   const double nnlo_qcd = 0.000641623890917771*Quad(g3)*(370.1956513893174 +
      2.375*l + (-47.18640261449638 + 0.6666666666666666*l)*Nf +
      0.9017702481178881*Sqr(Nf));
   const double nnnlo_qcd = 0.000016252523020247696*Power6(g3)*(467.683620788 +
      122.440972222*l + 10.9409722222*Sqr(l));

   return Sqrt(1.0 + nlo_qcd + nnlo_qcd + nnnlo_qcd);
}

double THDMIISNMSSMBCsimple_effective_couplings::pseudoscalar_scaling_factor(double m) const
{
   const double Nf = number_of_active_flavours(m);
   const double mtpole = qedqcd.displayPoleMt();
   const double l = Log(Sqr(m) / Sqr(mtpole));

   const auto g3 = MODELPARAMETER(g3);

   const double nlo_qcd = 0.025330295910584444*(24.25 - 1.1666666666666667*Nf)*Sqr
      (g3);
   const double nnlo_qcd = 0.000641623890917771*(171.54400563089382 + 5*l)*Quad(g3
      );
   const double nnnlo_qcd = 0;

   return Sqrt(1.0 + nlo_qcd + nnlo_qcd + nnnlo_qcd);
}

double THDMIISNMSSMBCsimple_effective_couplings::get_hhVPVP_partial_width(int gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CphhVPVP(gO1));
}

double THDMIISNMSSMBCsimple_effective_couplings::get_hhVGVG_partial_width(int gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CphhVGVG(gO1));
}

double THDMIISNMSSMBCsimple_effective_couplings::get_AhVPVP_partial_width(int gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CpAhVPVP(gO1));
}

double THDMIISNMSSMBCsimple_effective_couplings::get_AhVGVG_partial_width(int gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CpAhVGVG(gO1));
}

std::complex<double> THDMIISNMSSMBCsimple_effective_couplings::CphhconjVWmVWm(int gI2) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);

   const std::complex<double> result = 0.5*Sqr(g2)*(v1*ZH(gI2,0) + v2*ZH(gI2,1));

   return result;
}

std::complex<double> THDMIISNMSSMBCsimple_effective_couplings::CpbarFdFdAhPL(int gO2, int gI1, int gI2) const
{
   const auto Yd = MODELPARAMETER(Yd);

   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j2,0,2,Conj(Vd(gI1,j2))*SUM(j1,0,2,Conj(Ud(gO2,j1))*Yd(j1,j2)))*ZA(gI2,
      0);

   return result;
}

std::complex<double> THDMIISNMSSMBCsimple_effective_couplings::CpbarFdFdhhPL(int gO2, int gI2, int gI1) const
{
   const auto Yd = MODELPARAMETER(Yd);

   const std::complex<double> result = 0.7071067811865475*SUM(j2,0,2,Conj(Vd(gI2,
      j2))*SUM(j1,0,2,Conj(Ud(gO2,j1))*Yd(j1,j2)))*ZH(gI1,0);

   return result;
}

std::complex<double> THDMIISNMSSMBCsimple_effective_couplings::CpbarFeFeAhPL(int gO2, int gI1, int gI2) const
{
   const auto Ye = MODELPARAMETER(Ye);

   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j2,0,2,Conj(Ve(gI1,j2))*SUM(j1,0,2,Conj(Ue(gO2,j1))*Ye(j1,j2)))*ZA(gI2,
      0);

   return result;
}

std::complex<double> THDMIISNMSSMBCsimple_effective_couplings::CpbarFeFehhPL(int gO2, int gI2, int gI1) const
{
   const auto Ye = MODELPARAMETER(Ye);

   const std::complex<double> result = 0.7071067811865475*SUM(j2,0,2,Conj(Ve(gI2,
      j2))*SUM(j1,0,2,Conj(Ue(gO2,j1))*Ye(j1,j2)))*ZH(gI1,0);

   return result;
}

std::complex<double> THDMIISNMSSMBCsimple_effective_couplings::CpbarFuFuAhPL(int gO2, int gI1, int gI2) const
{
   const auto Yu = MODELPARAMETER(Yu);

   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j2,0,2,Conj(Vu(gI1,j2))*SUM(j1,0,2,Conj(Uu(gO2,j1))*Yu(j1,j2)))*ZA(gI2,
      1);

   return result;
}

std::complex<double> THDMIISNMSSMBCsimple_effective_couplings::CpbarFuFuhhPL(int gO2, int gI2, int gI1) const
{
   const auto Yu = MODELPARAMETER(Yu);

   const std::complex<double> result = 0.7071067811865475*SUM(j2,0,2,Conj(Vu(gI2,
      j2))*SUM(j1,0,2,Conj(Uu(gO2,j1))*Yu(j1,j2)))*ZH(gI1,1);

   return result;
}

std::complex<double> THDMIISNMSSMBCsimple_effective_couplings::CphhHmconjHm(int gt1, int gt2, int gt3) const
{
   const auto Lambda1 = MODELPARAMETER(Lambda1);
   const auto Lambda2 = MODELPARAMETER(Lambda2);
   const auto Lambda3 = MODELPARAMETER(Lambda3);
   const auto Lambda4 = MODELPARAMETER(Lambda4);
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto M123 = MODELPARAMETER(M123);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto vS = MODELPARAMETER(vS);

   const std::complex<double> result = 0.5*(ZH(gt1,1)*(ZP(gt2,0)*(2*(-Lambda3 +
      Lambda4)*v2*ZP(gt3,0) + Lambda4*v1*ZP(gt3,1)) + ZP(gt2,1)*(Lambda4*v1*ZP(gt3
      ,0) - 2*Lambda2*v2*ZP(gt3,1))) + ZH(gt1,0)*(ZP(gt2,1)*(Lambda4*v2*ZP(gt3,0)
      + 2*(-Lambda3 + Lambda4)*v1*ZP(gt3,1)) + ZP(gt2,0)*(-2*Lambda1*v1*ZP(gt3,0)
      + Lambda4*v2*ZP(gt3,1))) - ZH(gt1,2)*(ZP(gt2,1)*(-2*vS*Conj(Lambda7)*ZP(gt3,
      0) + 1.4142135623730951*Conj(M123)*ZP(gt3,0) + 2*Lambda6*vS*ZP(gt3,1)) + ZP(
      gt2,0)*(2*Lambda5*vS*ZP(gt3,0) + (1.4142135623730951*M123 - 2*Lambda7*vS)*ZP
      (gt3,1))));

   return result;
}

std::complex<double> THDMIISNMSSMBCsimple_effective_couplings::CpAhHmconjHm(int gt1, int gt2, int gt3) const
{
   const auto Lambda4 = MODELPARAMETER(Lambda4);
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto M123 = MODELPARAMETER(M123);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto vS = MODELPARAMETER(vS);

   const std::complex<double> result = std::complex<double>(0,0.5)*(Lambda4*v2*ZA(
      gt1,0)*(ZP(gt2,1)*ZP(gt3,0) - ZP(gt2,0)*ZP(gt3,1)) + Lambda4*v1*ZA(gt1,1)*(
      ZP(gt2,1)*ZP(gt3,0) - ZP(gt2,0)*ZP(gt3,1)) + ZA(gt1,2)*(2*vS*Conj(Lambda7)*
      ZP(gt2,1)*ZP(gt3,0) + 1.4142135623730951*Conj(M123)*ZP(gt2,1)*ZP(gt3,0) - (
      1.4142135623730951*M123 + 2*Lambda7*vS)*ZP(gt2,0)*ZP(gt3,1)));

   return result;
}

void THDMIISNMSSMBCsimple_effective_couplings::calculate_eff_CphhVPVP(int gO1)
{
   const auto MHm = PHYSICAL(MHm);
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const auto MFe = PHYSICAL(MFe);
   const auto MVWm = PHYSICAL(MVWm);
   const auto decay_mass = PHYSICAL(Mhh)(gO1);
   const auto decay_scale = 0.25 * Sqr(decay_mass);
   const auto saved_ZH = ZH;
   ZH = PHYSICAL(ZH);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (int gI1 = 1; gI1 < 2; ++gI1) {
      result += 0.5 * CphhHmconjHm(gO1, gI1, gI1) * vev * AS0(decay_scale / Sqr(
         MHm(gI1))) / Sqr(MHm(gI1));
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += 0.3333333333333333 * scalar_fermion_qcd_factor(decay_mass, MFd(gI1
         )) * CpbarFdFdhhPL(gI1, gI1, gO1) * vev * AS12(decay_scale / Sqr(MFd(gI1)
         )) / MFd(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += 1.3333333333333333 * scalar_fermion_qcd_factor(decay_mass, MFu(gI1
         )) * CpbarFuFuhhPL(gI1, gI1, gO1) * vev * AS12(decay_scale / Sqr(MFu(gI1)
         )) / MFu(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFeFehhPL(gI1, gI1, gO1) * vev * AS12(decay_scale / Sqr(MFe(
         gI1))) / MFe(gI1);
   }
   result += -0.5 * CphhconjVWmVWm(gO1) * vev * AS1(decay_scale / Sqr(MVWm)) / Sqr
      (MVWm);


   result *= 0.1892681907127351 * physical_input.get(Physical_input::alpha_em_0) *
      Sqrt(qedqcd.displayFermiConstant());

   ZH = saved_ZH;
   eff_CphhVPVP(gO1) = result;

}

void THDMIISNMSSMBCsimple_effective_couplings::calculate_eff_CphhVGVG(int gO1)
{
   const auto g3 = MODELPARAMETER(g3);
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const double alpha_s = 0.07957747154594767*Sqr(g3);
   const auto decay_mass = PHYSICAL(Mhh)(gO1);
   const auto decay_scale = 0.25 * Sqr(decay_mass);
   const auto saved_ZH = ZH;
   ZH = PHYSICAL(ZH);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFdFdhhPL(gI1, gI1, gO1) * vev * AS12(decay_scale / Sqr(MFd(
         gI1))) / MFd(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFuFuhhPL(gI1, gI1, gO1) * vev * AS12(decay_scale / Sqr(MFu(
         gI1))) / MFu(gI1);
   }
   result *= 0.75;

   if (include_qcd_corrections) {
      result *= scalar_scaling_factor(decay_mass);
   }


   result *= 0.12617879380849006 * alpha_s * Sqrt(qedqcd.displayFermiConstant());

   ZH = saved_ZH;
   eff_CphhVGVG(gO1) = result;

}

void THDMIISNMSSMBCsimple_effective_couplings::calculate_eff_CpAhVPVP(int gO1)
{
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const auto MFe = PHYSICAL(MFe);
   const auto decay_mass = PHYSICAL(MAh)(gO1);
   const auto decay_scale = 0.25 * Sqr(decay_mass);
   const auto saved_ZA = ZA;
   ZA = PHYSICAL(ZA);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += 0.3333333333333333 * pseudoscalar_fermion_qcd_factor(decay_mass,
         MFd(gI1)) * CpbarFdFdAhPL(gI1, gI1, gO1) * vev * AP12(decay_scale / Sqr(
         MFd(gI1))) / MFd(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += 1.3333333333333333 * pseudoscalar_fermion_qcd_factor(decay_mass,
         MFu(gI1)) * CpbarFuFuAhPL(gI1, gI1, gO1) * vev * AP12(decay_scale / Sqr(
         MFu(gI1))) / MFu(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFeFeAhPL(gI1, gI1, gO1) * vev * AP12(decay_scale / Sqr(MFe(
         gI1))) / MFe(gI1);
   }
   result *= 2.0;


   result *= 0.1892681907127351 * physical_input.get(Physical_input::alpha_em_0) *
      Sqrt(qedqcd.displayFermiConstant());

   ZA = saved_ZA;
   eff_CpAhVPVP(gO1) = result;

}

void THDMIISNMSSMBCsimple_effective_couplings::calculate_eff_CpAhVGVG(int gO1)
{
   const auto g3 = MODELPARAMETER(g3);
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const double alpha_s = 0.07957747154594767*Sqr(g3);
   const auto decay_mass = PHYSICAL(MAh)(gO1);
   const auto decay_scale = 0.25 * Sqr(decay_mass);
   const auto saved_ZA = ZA;
   ZA = PHYSICAL(ZA);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFdFdAhPL(gI1, gI1, gO1) * vev * AP12(decay_scale / Sqr(MFd(
         gI1))) / MFd(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFuFuAhPL(gI1, gI1, gO1) * vev * AP12(decay_scale / Sqr(MFu(
         gI1))) / MFu(gI1);
   }
   result *= 1.5;

   if (include_qcd_corrections) {
      result *= pseudoscalar_scaling_factor(decay_mass);
   }


   result *= 0.12617879380849006 * alpha_s * Sqrt(qedqcd.displayFermiConstant());

   ZA = saved_ZA;
   eff_CpAhVGVG(gO1) = result;

}


} // namespace flexiblesusy
