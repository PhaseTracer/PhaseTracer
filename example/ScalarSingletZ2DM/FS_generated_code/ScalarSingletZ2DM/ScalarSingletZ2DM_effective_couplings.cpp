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

// File generated at Tue 17 Nov 2020 16:11:28

#include "ScalarSingletZ2DM_effective_couplings.hpp"

#include "effective_couplings.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

using namespace effective_couplings;

#define INPUTPARAMETER(parameter) model.get_input().parameter
#define MODELPARAMETER(parameter) model.get_##parameter()
#define DERIVEDPARAMETER(parameter) model.##parameter()
#define PHASE(parameter) model.get_##parameter()
#define PHYSICAL(parameter) model.get_physical().parameter

ScalarSingletZ2DM_effective_couplings::ScalarSingletZ2DM_effective_couplings(
   const ScalarSingletZ2DM_mass_eigenstates& model_,
   const softsusy::QedQcd& qedqcd_,
   const Physical_input& input_)
   : model(model_), qedqcd(qedqcd_), physical_input(input_)
   , rg_improve(true), include_qcd_corrections(true)
   , Vd(PHYSICAL(Vd)), Ud(PHYSICAL(Ud)), Vu(PHYSICAL(Vu)), Uu(PHYSICAL(Uu)), Ve(
      PHYSICAL(Ve)), Ue(PHYSICAL(Ue)), ZZ(PHYSICAL(ZZ))
   , eff_CphhVPVP(0), eff_CphhVGVG(0), eff_CpAhVPVP(0), eff_CpAhVGVG(0)
{
}

void ScalarSingletZ2DM_effective_couplings::calculate_effective_couplings()
{
   const standard_model::Standard_model sm(initialise_SM());

   const double scale = model.get_scale();
   const Eigen::ArrayXd saved_parameters(model.get());

   const double saved_mt = PHYSICAL(MFu(2));
   PHYSICAL(MFu(2)) = qedqcd.displayPoleMt();

   const auto Mhh = PHYSICAL(Mhh);
   run_SM_strong_coupling_to(sm, 0.5 * Mhh);
   calculate_eff_CphhVPVP();
   run_SM_strong_coupling_to(sm, Mhh);
   calculate_eff_CphhVGVG();

   PHYSICAL(MFu(2)) = saved_mt;
   model.set_scale(scale);
   model.set(saved_parameters);

}

void ScalarSingletZ2DM_effective_couplings::set_model(const ScalarSingletZ2DM_mass_eigenstates& model_)
{
   model = model_;
   copy_mixing_matrices_from_model();
}

void ScalarSingletZ2DM_effective_couplings::copy_mixing_matrices_from_model()
{
   Vd = PHYSICAL(Vd);
   Ud = PHYSICAL(Ud);
   Vu = PHYSICAL(Vu);
   Uu = PHYSICAL(Uu);
   Ve = PHYSICAL(Ve);
   Ue = PHYSICAL(Ue);
   ZZ = PHYSICAL(ZZ);

}

standard_model::Standard_model ScalarSingletZ2DM_effective_couplings::initialise_SM() const
{
   standard_model::Standard_model sm;

   sm.set_loops(2);
   sm.set_thresholds(2);
   sm.set_physical_input(physical_input);

   sm.initialise_from_input(qedqcd);

   return sm;
}

void ScalarSingletZ2DM_effective_couplings::run_SM_strong_coupling_to(standard_model::Standard_model sm, double m)
{
   sm.run_to(m);

   model.set_g3(sm.get_g3());

}

std::complex<double> ScalarSingletZ2DM_effective_couplings::scalar_scalar_qcd_factor(double m_decay, double m_loop) const
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

std::complex<double> ScalarSingletZ2DM_effective_couplings::scalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) * scalar_diphoton_fermion_loop(
         m_decay, m_loop);

   }

   return result;
}

std::complex<double> ScalarSingletZ2DM_effective_couplings::pseudoscalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         pseudoscalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

double ScalarSingletZ2DM_effective_couplings::number_of_active_flavours(double m) const
{
   if (m < qedqcd.displayMbMb()) {
      return 4.0;
   } else if (m < qedqcd.displayPoleMt()) {
      return 5.0;
   } else {
      return 6.0;
   }
}

double ScalarSingletZ2DM_effective_couplings::scalar_scaling_factor(double m) const
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

double ScalarSingletZ2DM_effective_couplings::pseudoscalar_scaling_factor(double m) const
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

double ScalarSingletZ2DM_effective_couplings::get_hhVPVP_partial_width() const
{
   const double mass = PHYSICAL(Mhh);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CphhVPVP);
}

double ScalarSingletZ2DM_effective_couplings::get_hhVGVG_partial_width() const
{
   const double mass = PHYSICAL(Mhh);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CphhVGVG);
}

double ScalarSingletZ2DM_effective_couplings::get_AhVPVP_partial_width() const
{
   const double mass = PHYSICAL(MAh);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CpAhVPVP);
}

double ScalarSingletZ2DM_effective_couplings::get_AhVGVG_partial_width() const
{
   const double mass = PHYSICAL(MAh);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CpAhVGVG);
}

std::complex<double> ScalarSingletZ2DM_effective_couplings::CpbarFdFdAhPL(int gI1, int gI2) const
{
   const auto Yd = MODELPARAMETER(Yd);

   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j2,0,2,Conj(Vd(gI2,j2))*SUM(j1,0,2,Conj(Ud(gI1,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> ScalarSingletZ2DM_effective_couplings::CpbarFeFeAhPL(int gI1, int gI2) const
{
   const auto Ye = MODELPARAMETER(Ye);

   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j2,0,2,Conj(Ve(gI2,j2))*SUM(j1,0,2,Conj(Ue(gI1,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> ScalarSingletZ2DM_effective_couplings::CpbarFuFuAhPL(int gI1, int gI2) const
{
   const auto Yu = MODELPARAMETER(Yu);

   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j2,0,2,Conj(Vu(gI2,j2))*SUM(j1,0,2,Conj(Uu(gI1,j1))*Yu(j1,j2)));

   return result;
}

double ScalarSingletZ2DM_effective_couplings::CphhconjVWpVWp() const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);

   const double result = 0.5*v*Sqr(g2);

   return result;
}

std::complex<double> ScalarSingletZ2DM_effective_couplings::CpbarFdFdhhPL(int gI1, int gI2) const
{
   const auto Yd = MODELPARAMETER(Yd);

   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(Vd(gI2,
      j2))*SUM(j1,0,2,Conj(Ud(gI1,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> ScalarSingletZ2DM_effective_couplings::CpbarFeFehhPL(int gI1, int gI2) const
{
   const auto Ye = MODELPARAMETER(Ye);

   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(Ve(gI2,
      j2))*SUM(j1,0,2,Conj(Ue(gI1,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> ScalarSingletZ2DM_effective_couplings::CpbarFuFuhhPL(int gI1, int gI2) const
{
   const auto Yu = MODELPARAMETER(Yu);

   const std::complex<double> result = 0.7071067811865475*SUM(j2,0,2,Conj(Vu(gI2,
      j2))*SUM(j1,0,2,Conj(Uu(gI1,j1))*Yu(j1,j2)));

   return result;
}

void ScalarSingletZ2DM_effective_couplings::calculate_eff_CphhVPVP()
{
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const auto MFe = PHYSICAL(MFe);
   const auto MVWp = PHYSICAL(MVWp);
   const auto decay_mass = PHYSICAL(Mhh);
   const auto decay_scale = 0.25 * Sqr(decay_mass);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += 0.3333333333333333 * scalar_fermion_qcd_factor(decay_mass, MFd(gI1
         )) * CpbarFdFdhhPL(gI1, gI1) * vev * AS12(decay_scale / Sqr(MFd(gI1))) /
         MFd(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += 1.3333333333333333 * scalar_fermion_qcd_factor(decay_mass, MFu(gI1
         )) * CpbarFuFuhhPL(gI1, gI1) * vev * AS12(decay_scale / Sqr(MFu(gI1))) /
         MFu(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFeFehhPL(gI1, gI1) * vev * AS12(decay_scale / Sqr(MFe(gI1)))
         / MFe(gI1);
   }
   result += -0.5 * CphhconjVWpVWp() * vev * AS1(decay_scale / Sqr(MVWp)) / Sqr(
      MVWp);


   result *= 0.1892681907127351 * physical_input.get(Physical_input::alpha_em_0) *
      Sqrt(qedqcd.displayFermiConstant());

   eff_CphhVPVP = result;

}

void ScalarSingletZ2DM_effective_couplings::calculate_eff_CphhVGVG()
{
   const auto g3 = MODELPARAMETER(g3);
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const double alpha_s = 0.07957747154594767*Sqr(g3);
   const auto decay_mass = PHYSICAL(Mhh);
   const auto decay_scale = 0.25 * Sqr(decay_mass);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFdFdhhPL(gI1, gI1) * vev * AS12(decay_scale / Sqr(MFd(gI1)))
         / MFd(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFuFuhhPL(gI1, gI1) * vev * AS12(decay_scale / Sqr(MFu(gI1)))
         / MFu(gI1);
   }
   result *= 0.75;

   if (include_qcd_corrections) {
      result *= scalar_scaling_factor(decay_mass);
   }


   result *= 0.12617879380849006 * alpha_s * Sqrt(qedqcd.displayFermiConstant());

   eff_CphhVGVG = result;

}

void ScalarSingletZ2DM_effective_couplings::calculate_eff_CpAhVPVP()
{
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const auto MFe = PHYSICAL(MFe);
   const auto decay_mass = PHYSICAL(MAh);
   const auto decay_scale = 0.25 * Sqr(decay_mass);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += 0.3333333333333333 * pseudoscalar_fermion_qcd_factor(decay_mass,
         MFd(gI1)) * CpbarFdFdAhPL(gI1, gI1) * vev * AP12(decay_scale / Sqr(MFd(
         gI1))) / MFd(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += 1.3333333333333333 * pseudoscalar_fermion_qcd_factor(decay_mass,
         MFu(gI1)) * CpbarFuFuAhPL(gI1, gI1) * vev * AP12(decay_scale / Sqr(MFu(
         gI1))) / MFu(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFeFeAhPL(gI1, gI1) * vev * AP12(decay_scale / Sqr(MFe(gI1)))
         / MFe(gI1);
   }
   result *= 2.0;


   result *= 0.1892681907127351 * physical_input.get(Physical_input::alpha_em_0) *
      Sqrt(qedqcd.displayFermiConstant());

   eff_CpAhVPVP = result;

}

void ScalarSingletZ2DM_effective_couplings::calculate_eff_CpAhVGVG()
{
   const auto g3 = MODELPARAMETER(g3);
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const double alpha_s = 0.07957747154594767*Sqr(g3);
   const auto decay_mass = PHYSICAL(MAh);
   const auto decay_scale = 0.25 * Sqr(decay_mass);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFdFdAhPL(gI1, gI1) * vev * AP12(decay_scale / Sqr(MFd(gI1)))
         / MFd(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFuFuAhPL(gI1, gI1) * vev * AP12(decay_scale / Sqr(MFu(gI1)))
         / MFu(gI1);
   }
   result *= 1.5;

   if (include_qcd_corrections) {
      result *= pseudoscalar_scaling_factor(decay_mass);
   }


   result *= 0.12617879380849006 * alpha_s * Sqrt(qedqcd.displayFermiConstant());

   eff_CpAhVGVG = result;

}


} // namespace flexiblesusy
