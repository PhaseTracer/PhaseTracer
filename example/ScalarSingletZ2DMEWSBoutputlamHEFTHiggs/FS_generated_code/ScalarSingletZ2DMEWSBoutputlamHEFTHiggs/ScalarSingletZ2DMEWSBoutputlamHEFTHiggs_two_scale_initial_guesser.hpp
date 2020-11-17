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

// File generated at Tue 17 Nov 2020 15:32:58

#ifndef ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_STANDARD_MODEL_TWO_SCALE_INITIAL_GUESSER_H
#define ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_STANDARD_MODEL_TWO_SCALE_INITIAL_GUESSER_H

#include "ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_initial_guesser.hpp"
#include "ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_two_scale_susy_scale_constraint.hpp"
#include "ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_two_scale_high_scale_constraint.hpp"
#include "standard_model_two_scale_low_scale_constraint.hpp"
#include "initial_guesser.hpp"
#include "lowe.h"

#include <sstream>
#include <Eigen/Core>

namespace flexiblesusy {

class Two_scale;

template <class T>
class ScalarSingletZ2DMEWSBoutputlamHEFTHiggs;

template <class T>
class StandardModel;

template <class T>
class ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_standard_model_initial_guesser;

/**
 * @class ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_standard_model_initial_guesser<Two_scale>
 * @brief initial guesser for the ScalarSingletZ2DMEWSBoutputlamHEFTHiggs tower
 */

template<>
class ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_standard_model_initial_guesser<Two_scale> : public Initial_guesser {
public:
   ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_standard_model_initial_guesser(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs<Two_scale>*,
                               standard_model::StandardModel<Two_scale>*,
                               const softsusy::QedQcd&,
                               const standard_model::Standard_model_low_scale_constraint<Two_scale>&,
                               const ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_scale_constraint<Two_scale>&,
                               const ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_high_scale_constraint<Two_scale>&);
   virtual ~ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_standard_model_initial_guesser();
   virtual void guess(); ///< initial guess

   void set_running_precision(double p) { running_precision = p; }

private:
   ScalarSingletZ2DMEWSBoutputlamHEFTHiggs<Two_scale>* model{nullptr}; ///< pointer to model class
   standard_model::StandardModel<Two_scale>* eft{nullptr}; ///< pointer to effective low energy model
   softsusy::QedQcd qedqcd{}; ///< Standard Model low-energy data
   double mu_guess{0.}; ///< guessed DR-bar mass of up-quark
   double mc_guess{0.}; ///< guessed DR-bar mass of charm-quark
   double mt_guess{0.}; ///< guessed DR-bar mass of top-quark
   double md_guess{0.}; ///< guessed DR-bar mass of down-quark
   double ms_guess{0.}; ///< guessed DR-bar mass of strange-quark
   double mb_guess{0.}; ///< guessed DR-bar mass of bottom-quark
   double me_guess{0.}; ///< guessed DR-bar mass of electron
   double mm_guess{0.}; ///< guessed DR-bar mass of muon
   double mtau_guess{0.}; ///< guessed DR-bar mass of tau
   double running_precision{1.0e-3}; ///< Runge-Kutta RG running precision
   standard_model::Standard_model_low_scale_constraint<Two_scale> low_constraint{};
   ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_scale_constraint<Two_scale> susy_constraint{};
   ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_high_scale_constraint<Two_scale> high_constraint{};
   Eigen::Matrix<std::complex<double>,3,3> upQuarksDRbar{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> downQuarksDRbar{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> downLeptonsDRbar{Eigen::Matrix<std::complex<double>,3,3>::Zero()};

   void guess_eft_parameters();
   void guess_model_parameters();
   void calculate_DRbar_yukawa_couplings();
   void calculate_Yu_DRbar();
   void calculate_Yd_DRbar();
   void calculate_Ye_DRbar();
   void calculate_running_SM_masses();
};

} // namespace flexiblesusy

#endif
