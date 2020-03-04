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

// File generated at Thu 7 Nov 2019 18:53:58

/**
 * @file THDMIISNMSSMBCsimple_a_muon.cpp
 *
 * This file was generated at Thu 7 Nov 2019 18:53:58 with FlexibleSUSY
 * 2.4.1 and SARAH 4.14.2 .
 */

#include "THDMIISNMSSMBCsimple_a_muon.hpp"
#include "THDMIISNMSSMBCsimple_mass_eigenstates.hpp"

#include "cxx_qft/THDMIISNMSSMBCsimple_qft.hpp"
#include "THDMIISNMSSMBCsimple_FFV_form_factors.hpp"

#include "lowe.h"
#include "wrappers.hpp"
#include "numerics2.hpp"
#include "logger.hpp"

using namespace flexiblesusy;
using namespace THDMIISNMSSMBCsimple_cxx_diagrams;

using Muon = fields::Fe;

namespace {

double get_QED_2L(context_base&, const softsusy::QedQcd&);

double get_MSUSY(const THDMIISNMSSMBCsimple_mass_eigenstates& model)
{
   return Min(model.get_MHm().tail<1>().minCoeff());
}

void run_to_MSUSY(THDMIISNMSSMBCsimple_mass_eigenstates& model)
{
   const double precision_goal = model.get_precision();
   const int Nmax = static_cast<int>(-log10(precision_goal) * 10);
   int it = 0;
   double precision = 0.;
   double MSUSY_old = 0., MSUSY_new = 0.;

   VERBOSE_MSG("THDMIISNMSSMBCsimple_a_muon: iterate to run to MSUSY ...");

   do {
      MSUSY_new = get_MSUSY(model);

      if (MSUSY_new <= 0.)
         return;

      VERBOSE_MSG("THDMIISNMSSMBCsimple_a_muon:    running to new MSUSY = " << MSUSY_new);

      model.run_to(MSUSY_new);
      model.calculate_DRbar_masses();

      precision = MaxRelDiff(MSUSY_old, MSUSY_new);
      MSUSY_old = MSUSY_new;
   } while (precision > precision_goal && ++it < Nmax);

   VERBOSE_MSG("THDMIISNMSSMBCsimple_a_muon: iteration finished. MSUSY = " << MSUSY_new);

   if (precision > precision_goal) {
      WARNING("THDMIISNMSSMBCsimple_a_muon: run_to_MSUSY(): "
              "Iteration did not converge after " << Nmax
              << " iterations (precision goal: " << precision_goal << ")");
   }
}

double muonPhysicalMass(const softsusy::QedQcd& qedqcd)
{
   return qedqcd.displayPoleMmuon();
}

double calculate_a_muon_impl(const THDMIISNMSSMBCsimple_mass_eigenstates& model, const softsusy::QedQcd& qedqcd)
{
   VERBOSE_MSG("THDMIISNMSSMBCsimple_a_muon: calculating a_mu at Q = " << model.get_scale());

   context_base context{ model };
   
   using namespace THDMIISNMSSMBCsimple_cxx_diagrams::fields;

   const std::valarray<std::complex<double>> form_factors {0., 0., 0., 0.};

   if (!is_zero((form_factors[2] + form_factors[3]).imag())) {
      ERROR("Error in the g-2 calculation! Form factor F2 should be real");
      return std::numeric_limits<double>::quiet_NaN();
   }

   double val = 
      // vector form factor
      1./2.*(form_factors[2] + form_factors[3]).real()
      // reinstate the mass that was factored out in definition of form factor 
      * context.mass<Muon>({1}) 
      // factor out e/(2*muon pole mass)
      / (unit_charge(context_base{model})/(2.*muonPhysicalMass(qedqcd))) 
      // definition of photon momentum flow for g-2 is opposite than for form factors
      * (-1.);

   // add 2-loop QED logarithms
   val *= 1. + get_QED_2L(context, qedqcd);

   return val;
}

/// generates array with N scales from mean/factor to mean*factor
template <int N>
std::array<double,N> generate_scales(double mean, double factor)
{
   static_assert(N > 1, "N must be greater than 1!");

   const double start = mean / factor, stop = mean * factor;
   std::array<double,N> scales;

   scales[0] = start;

   for (int i = 1; i < (N-1); i++)
      scales[i] = std::exp(std::log(start) + (std::log(stop) - std::log(start))*i / N);

   scales[N-1] = stop;

   return scales;
}

/// returns minimum and maximum a_mu when scale is varied by a factor 2
std::pair<double,double> vary_scale(const THDMIISNMSSMBCsimple_mass_eigenstates& model, const softsusy::QedQcd& qedqcd)
{
   auto scales = generate_scales<7>(model.get_scale(), 2.);

   std::transform(scales.begin(), scales.end(), scales.begin(),
                  [&model,&qedqcd] (double scale) {
                     double amu = 0.;
                     try {
                        auto m = model;
                        m.run_to(scale);
                        m.get_physical().clear();
                        m.calculate_DRbar_masses();
                        m.solve_ewsb();
                        m.calculate_MFe_pole();
                        amu = calculate_a_muon_impl(m, qedqcd);
                     }
                     catch(const Error& e) {
                        ERROR("THDMIISNMSSMBCsimple_a_muon: scale variation: " << e.what_detailed());
                     }
                     return amu;
                  });

   const auto minmax = std::minmax_element(scales.cbegin(), scales.cend());

   return std::make_pair(*(minmax.first), *(minmax.second));
}

} // anonymous namespace

double THDMIISNMSSMBCsimple_a_muon::calculate_a_muon(const THDMIISNMSSMBCsimple_mass_eigenstates& model_, const softsusy::QedQcd& qedqcd)
{
   THDMIISNMSSMBCsimple_mass_eigenstates model(model_);

   VERBOSE_MSG("THDMIISNMSSMBCsimple_a_muon: starting calculation of a_mu ...");

   try {
      run_to_MSUSY(model);
      model.get_physical().clear();
   } catch (const Error& e) {
      ERROR("THDMIISNMSSMBCsimple_a_muon:" << e.what_detailed());
      return std::numeric_limits<double>::quiet_NaN();
   }

   double m_muon_pole = muonPhysicalMass(qedqcd);

   if (m_muon_pole == 0.0) {
      model.solve_ewsb();
      model.calculate_MFe_pole();
   }

   return calculate_a_muon_impl(model, qedqcd);
}

double THDMIISNMSSMBCsimple_a_muon::calculate_a_muon_uncertainty(const THDMIISNMSSMBCsimple_mass_eigenstates& model_, const softsusy::QedQcd& qedqcd)
{
   THDMIISNMSSMBCsimple_mass_eigenstates model(model_);

   VERBOSE_MSG("THDMIISNMSSMBCsimple_a_muon: starting calculation of a_mu uncertainty ...");

   try {
      run_to_MSUSY(model);
      model.get_physical().clear();
   } catch (const Error& e) {
      ERROR("THDMIISNMSSMBCsimple_a_muon uncertainty: " << e.what_detailed());
      return std::numeric_limits<double>::quiet_NaN();
   }

   const auto delta_amu_scale_minmax = vary_scale(model, qedqcd);
   const auto delta_amu_scale = std::abs(delta_amu_scale_minmax.second - delta_amu_scale_minmax.first);

   return delta_amu_scale;
}

namespace {
double get_QED_2L(context_base& context, const softsusy::QedQcd& qedqcd)
{
   const double MSUSY = Abs(get_MSUSY(context.model));
   const double m_muon = muonPhysicalMass(qedqcd);
   const double alpha_em = Sqr(Muon::electric_charge * unit_charge(context))/(4*Pi);
   const double qed_2L = alpha_em/(4*Pi) * 16 * FiniteLog(m_muon/MSUSY);

   return qed_2L;
}

} // anonymous namespace
