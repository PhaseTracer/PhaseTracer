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


#include "ScalarSingletZ2DMMhInputMsInput_observables.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_mass_eigenstates.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_a_muon.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_edm.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_l_to_lgamma.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_b_to_s_gamma.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_f_to_f_conversion.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_effective_couplings.hpp"
#include "config.h"
#include "eigen_utils.hpp"
#include "numerics2.hpp"
#include "wrappers.hpp"
#include "lowe.h"
#include "physical_input.hpp"

#ifdef ENABLE_GM2Calc
#include "gm2calc_interface.hpp"
#endif

#define MODEL model
#define AMU a_muon
#define AMUUNCERTAINTY a_muon_uncertainty
#define AMUGM2CALC a_muon_gm2calc
#define AMUGM2CALCUNCERTAINTY a_muon_gm2calc_uncertainty
#define EDM0(p) edm_ ## p
#define EDM1(p,idx) edm_ ## p ## _ ## idx
#define LToLGamma0(pIn, pOut, spec) pIn ## _to_ ## pOut ## _ ## spec
#define LToLGamma1(pIn,idxIn,pOut,idxOut,spec) pIn ## _to_ ## pOut ## _ ## spec
#define FToFConversion1(pIn,idxIn,pOut,idxOut,nuclei,qedqcd) pIn ## _to_ ## pOut ## _in_ ## nuclei
#define BSGAMMA b_to_s_gamma
#define EFFCPHIGGSPHOTONPHOTON eff_cp_higgs_photon_photon
#define EFFCPHIGGSGLUONGLUON eff_cp_higgs_gluon_gluon
#define EFFCPPSEUDOSCALARPHOTONPHOTON eff_cp_pseudoscalar_photon_photon
#define EFFCPPSEUDOSCALARGLUONGLUON eff_cp_pseudoscalar_gluon_gluon

#define ALPHA_S_MZ qedqcd.displayAlpha(softsusy::ALPHAS)
#define MWPole qedqcd.displayPoleMW()
#define MZPole qedqcd.displayPoleMZ()
#define MTPole qedqcd.displayPoleMt()
#define MBMB qedqcd.displayMbMb()
#define MTauPole qedqcd.displayPoleMtau()
#define MMPole qedqcd.displayPoleMmuon()

namespace flexiblesusy {

const int ScalarSingletZ2DMMhInputMsInput_observables::NUMBER_OF_OBSERVABLES;

ScalarSingletZ2DMMhInputMsInput_observables::ScalarSingletZ2DMMhInputMsInput_observables()
   : a_muon(0)
   , eff_cp_higgs_photon_photon(0)
   , eff_cp_higgs_gluon_gluon(0)

{
}

Eigen::ArrayXd ScalarSingletZ2DMMhInputMsInput_observables::get() const
{
   Eigen::ArrayXd vec(ScalarSingletZ2DMMhInputMsInput_observables::NUMBER_OF_OBSERVABLES);

   vec(0) = a_muon;
   vec(1) = Re(eff_cp_higgs_photon_photon);
   vec(2) = Im(eff_cp_higgs_photon_photon);
   vec(3) = Re(eff_cp_higgs_gluon_gluon);
   vec(4) = Im(eff_cp_higgs_gluon_gluon);

   return vec;
}

std::vector<std::string> ScalarSingletZ2DMMhInputMsInput_observables::get_names()
{
   std::vector<std::string> names(ScalarSingletZ2DMMhInputMsInput_observables::NUMBER_OF_OBSERVABLES);

   names[0] = "a_muon";
   names[1] = "Re(eff_cp_higgs_photon_photon)";
   names[2] = "Im(eff_cp_higgs_photon_photon)";
   names[3] = "Re(eff_cp_higgs_gluon_gluon)";
   names[4] = "Im(eff_cp_higgs_gluon_gluon)";

   return names;
}

void ScalarSingletZ2DMMhInputMsInput_observables::clear()
{
   a_muon = 0.;
   eff_cp_higgs_photon_photon = std::complex<double>(0.,0.);
   eff_cp_higgs_gluon_gluon = std::complex<double>(0.,0.);

}

void ScalarSingletZ2DMMhInputMsInput_observables::set(const Eigen::ArrayXd& vec)
{
   assert(vec.rows() == ScalarSingletZ2DMMhInputMsInput_observables::NUMBER_OF_OBSERVABLES);

   a_muon = vec(0);
   eff_cp_higgs_photon_photon = std::complex<double>(vec(1), vec(2));
   eff_cp_higgs_gluon_gluon = std::complex<double>(vec(3), vec(4));

}

ScalarSingletZ2DMMhInputMsInput_observables calculate_observables(ScalarSingletZ2DMMhInputMsInput_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input,
                                              double scale)
{
   auto model_at_scale = model;

   if (scale > 0.) {
      try {
         model_at_scale.run_to(scale);
      } catch (const Error& e) {
         model.get_problems().flag_thrown(e.what_detailed());
         return ScalarSingletZ2DMMhInputMsInput_observables();
      } catch (const std::exception& e) {
         model.get_problems().flag_thrown(e.what());
         return ScalarSingletZ2DMMhInputMsInput_observables();
      }
   }

   return calculate_observables(model_at_scale, qedqcd, physical_input);
}

ScalarSingletZ2DMMhInputMsInput_observables calculate_observables(ScalarSingletZ2DMMhInputMsInput_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input)
{
   ScalarSingletZ2DMMhInputMsInput_observables observables;

   try {
      ScalarSingletZ2DMMhInputMsInput_effective_couplings effective_couplings(model, qedqcd, physical_input);
      effective_couplings.calculate_effective_couplings();

      observables.AMU = ScalarSingletZ2DMMhInputMsInput_a_muon::calculate_a_muon(MODEL, qedqcd);
      observables.EFFCPHIGGSPHOTONPHOTON = effective_couplings.get_eff_CphhVPVP();
      observables.EFFCPHIGGSGLUONGLUON = effective_couplings.get_eff_CphhVGVG();
   } catch (const Error& e) {
      model.get_problems().flag_thrown(e.what_detailed());
   } catch (const std::exception& e) {
      model.get_problems().flag_thrown(e.what());
   }

   return observables;
}

} // namespace flexiblesusy
