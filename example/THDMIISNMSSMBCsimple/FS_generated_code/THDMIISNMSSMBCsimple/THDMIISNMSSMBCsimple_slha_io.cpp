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

// File generated at Thu 7 Nov 2019 18:53:49

#include "THDMIISNMSSMBCsimple_slha_io.hpp"
#include "THDMIISNMSSMBCsimple_input_parameters.hpp"
#include "logger.hpp"
#include "wrappers.hpp"
#include "numerics2.hpp"
#include "config.h"

#include <array>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

#define Pole(p) physical.p
#define PHYSICAL(p) model.get_physical().p
#define PHYSICAL_SLHA(p) model.get_physical_slha().p
#define LOCALPHYSICAL(p) physical.p
#define MODELPARAMETER(p) model.get_##p()
#define INPUTPARAMETER(p) input.p
#define EXTRAPARAMETER(p) model.get_##p()
#define DEFINE_PHYSICAL_PARAMETER(p) decltype(LOCALPHYSICAL(p)) p;
#define LowEnergyConstant(p) Electroweak_constants::p

namespace flexiblesusy {

THDMIISNMSSMBCsimple_slha_io::THDMIISNMSSMBCsimple_slha_io()
   : slha_io()
   , print_imaginary_parts_of_majorana_mixings(false)
{
}

void THDMIISNMSSMBCsimple_slha_io::clear()
{
   slha_io.clear();
}

void THDMIISNMSSMBCsimple_slha_io::set_print_imaginary_parts_of_majorana_mixings(bool flag)
{
   print_imaginary_parts_of_majorana_mixings = flag;
}

/**
 * Stores the EXTPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void THDMIISNMSSMBCsimple_slha_io::set_extpar(const THDMIISNMSSMBCsimple_input_parameters& input)
{

}

/**
 * Stores the IMMINPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void THDMIISNMSSMBCsimple_slha_io::set_imminpar(const THDMIISNMSSMBCsimple_input_parameters& input)
{

}

/**
 * Stores the IMEXTPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void THDMIISNMSSMBCsimple_slha_io::set_imextpar(const THDMIISNMSSMBCsimple_input_parameters& input)
{

}

/**
 * Stores the MODSEL input parameters in the SLHA object.
 *
 * @param modsel struct of MODSEL parameters
 */
void THDMIISNMSSMBCsimple_slha_io::set_modsel(const SLHA_io::Modsel& modsel)
{
   slha_io.set_modsel(modsel);
}

/**
 * Stores the MINPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void THDMIISNMSSMBCsimple_slha_io::set_minpar(const THDMIISNMSSMBCsimple_input_parameters& input)
{
   std::ostringstream minpar;

   minpar << "Block MINPAR\n";
   minpar << FORMAT_ELEMENT(1, input.lambdaNMSSM, "lambdaNMSSM");
   minpar << FORMAT_ELEMENT(2, input.kappaNMSSM, "kappaNMSSM");
   minpar << FORMAT_ELEMENT(3, input.AlambdaNMSSM, "AlambdaNMSSM");
   minpar << FORMAT_ELEMENT(4, input.AkappaNMSSM, "AkappaNMSSM");
   minpar << FORMAT_ELEMENT(5, input.AtopNMSSM, "AtopNMSSM");
   minpar << FORMAT_ELEMENT(6, input.mstopL, "mstopL");
   minpar << FORMAT_ELEMENT(7, input.mstopR, "mstopR");
   minpar << FORMAT_ELEMENT(8, input.MEWSB, "MEWSB");
   minpar << FORMAT_ELEMENT(10, input.TanBeta, "TanBeta");
   minpar << FORMAT_ELEMENT(13, input.vSIN, "vSIN");
   slha_io.set_block(minpar);

}

/**
 * Stores all input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void THDMIISNMSSMBCsimple_slha_io::set_input(const THDMIISNMSSMBCsimple_input_parameters& input)
{
   set_minpar(input);
   set_extpar(input);
   set_imminpar(input);
   set_imextpar(input);


}

/**
 * Stores the additional physical input (FlexibleSUSYInput block) in
 * the SLHA object.
 *
 * @param input class of input
 */
void THDMIISNMSSMBCsimple_slha_io::set_physical_input(const Physical_input& input)
{
   slha_io.set_physical_input(input);
}

/**
 * Stores the settings (FlexibleSUSY block) in the SLHA object.
 *
 * @param settings class of settings
 */
void THDMIISNMSSMBCsimple_slha_io::set_settings(const Spectrum_generator_settings& settings)
{
   slha_io.set_settings(settings);
}

/**
 * Stores the SMINPUTS input parameters in the SLHA object.
 *
 * @param qedqcd class of Standard Model parameters
 */
void THDMIISNMSSMBCsimple_slha_io::set_sminputs(const softsusy::QedQcd& qedqcd)
{
   slha_io.set_sminputs(qedqcd);
}

/**
 * Stores the spectrum generator information in the SPINFO block in
 * the SLHA object.
 *
 * @param problems struct with parameter point problems
 */
void THDMIISNMSSMBCsimple_slha_io::set_spinfo(const Spectrum_generator_problems& problems)
{
   set_spinfo(problems.get_problem_strings(), problems.get_warning_strings());
}

/**
 * Stores the spectrum generator information in the SPINFO block in
 * the SLHA object.
 *
 * @param problems struct with parameter point problems
 */
void THDMIISNMSSMBCsimple_slha_io::set_spinfo(const Problems& problems)
{
   set_spinfo(problems.get_problem_strings(), problems.get_warning_strings());
}

/**
 * Stores the given problems and warnings in the SPINFO block in the
 * SLHA object.
 *
 * @param problems vector of problem strings
 * @param warnings vector of warning strings
 */
void THDMIISNMSSMBCsimple_slha_io::set_spinfo(
   const std::vector<std::string>& problems,
   const std::vector<std::string>& warnings)
{
   std::ostringstream spinfo;
   spinfo << "Block SPINFO\n"
          << FORMAT_SPINFO(1, PKGNAME)
          << FORMAT_SPINFO(2, FLEXIBLESUSY_VERSION);

   for (const auto& s: warnings)
      spinfo << FORMAT_SPINFO(3, s);

   for (const auto& s: problems)
      spinfo << FORMAT_SPINFO(4, s);

   spinfo << FORMAT_SPINFO(5, THDMIISNMSSMBCsimple_info::model_name)
          << FORMAT_SPINFO(9, SARAH_VERSION);

   slha_io.set_block(spinfo, SLHA_io::front);
}

/**
 * Stores the particle masses in the SLHA object.
 *
 * @param physical struct of physical parameters
 *
 * @param write_sm_masses flag to indicate if Standard Model
 *    particle masses should be written as well
 */
void THDMIISNMSSMBCsimple_slha_io::set_mass(const THDMIISNMSSMBCsimple_physical& physical,
                                   bool write_sm_masses)
{
   std::ostringstream mass;

   mass << "Block MASS\n"
      << FORMAT_MASS(24, LOCALPHYSICAL(MVWm), "VWm")
      << FORMAT_MASS(37, LOCALPHYSICAL(MHm(1)), "Hm(2)")
      << FORMAT_MASS(25, LOCALPHYSICAL(Mhh(0)), "hh(1)")
      << FORMAT_MASS(35, LOCALPHYSICAL(Mhh(1)), "hh(2)")
      << FORMAT_MASS(45, LOCALPHYSICAL(Mhh(2)), "hh(3)")
      << FORMAT_MASS(36, LOCALPHYSICAL(MAh(1)), "Ah(2)")
      << FORMAT_MASS(46, LOCALPHYSICAL(MAh(2)), "Ah(3)")
   ;

   if (write_sm_masses) {
      mass
         << FORMAT_MASS(21, LOCALPHYSICAL(MVG), "VG")
         << FORMAT_MASS(12, LOCALPHYSICAL(MFv(0)), "Fv(1)")
         << FORMAT_MASS(14, LOCALPHYSICAL(MFv(1)), "Fv(2)")
         << FORMAT_MASS(16, LOCALPHYSICAL(MFv(2)), "Fv(3)")
         << FORMAT_MASS(1, LOCALPHYSICAL(MFd(0)), "Fd(1)")
         << FORMAT_MASS(3, LOCALPHYSICAL(MFd(1)), "Fd(2)")
         << FORMAT_MASS(5, LOCALPHYSICAL(MFd(2)), "Fd(3)")
         << FORMAT_MASS(2, LOCALPHYSICAL(MFu(0)), "Fu(1)")
         << FORMAT_MASS(4, LOCALPHYSICAL(MFu(1)), "Fu(2)")
         << FORMAT_MASS(6, LOCALPHYSICAL(MFu(2)), "Fu(3)")
         << FORMAT_MASS(11, LOCALPHYSICAL(MFe(0)), "Fe(1)")
         << FORMAT_MASS(13, LOCALPHYSICAL(MFe(1)), "Fe(2)")
         << FORMAT_MASS(15, LOCALPHYSICAL(MFe(2)), "Fe(3)")
         << FORMAT_MASS(22, LOCALPHYSICAL(MVP), "VP")
         << FORMAT_MASS(23, LOCALPHYSICAL(MVZ), "VZ")
      ;
   }

   slha_io.set_block(mass);

}

/**
 * Stores the mixing matrices in the SLHA object.
 *
 * @param physical struct of physical parameters
 *
 * @param write_sm_mixing_matrics flag to indicate if Standard Model
 *    particle mixing matrices should be written as well
 */
void THDMIISNMSSMBCsimple_slha_io::set_mixing_matrices(const THDMIISNMSSMBCsimple_physical& physical,
                                              bool write_sm_mixing_matrics)
{
   slha_io.set_block("PSEUDOSCALARMIX", LOCALPHYSICAL(ZA), "ZA");
   slha_io.set_block("SCALARMIX", LOCALPHYSICAL(ZH), "ZH");
   slha_io.set_block("CHARGEMIX", LOCALPHYSICAL(ZP), "ZP");

   if (write_sm_mixing_matrics) {
      slha_io.set_block("UULMIX", LOCALPHYSICAL(Vu), "Vu");
      slha_io.set_block("UDLMIX", LOCALPHYSICAL(Vd), "Vd");
      slha_io.set_block("UURMIX", LOCALPHYSICAL(Uu), "Uu");
      slha_io.set_block("UDRMIX", LOCALPHYSICAL(Ud), "Ud");
      slha_io.set_block("UELMIX", LOCALPHYSICAL(Ve), "Ve");
      slha_io.set_block("UERMIX", LOCALPHYSICAL(Ue), "Ue");
   }

   if (print_imaginary_parts_of_majorana_mixings) {
   }

}

void THDMIISNMSSMBCsimple_slha_io::set_ckm(
   const Eigen::Matrix<std::complex<double>,3,3>& ckm_matrix,
   double scale)
{
   slha_io.set_block("VCKM"  , ckm_matrix.real(), "Re(CKM)", scale);
   slha_io.set_block("IMVCKM", ckm_matrix.imag(), "Im(CKM)", scale);
}

void THDMIISNMSSMBCsimple_slha_io::set_pmns(
   const Eigen::Matrix<std::complex<double>,3,3>& pmns_matrix,
   double scale)
{
   slha_io.set_block("VPMNS"  , pmns_matrix.real(), "Re(PMNS)", scale);
   slha_io.set_block("IMVPMNS", pmns_matrix.imag(), "Im(PMNS)", scale);
}

void THDMIISNMSSMBCsimple_slha_io::set_model_parameters(const standard_model::Standard_model& model)
{
   {
      std::ostringstream block;
      block << "Block SMGAUGE Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (model.get_g1() * standard_model_info::normalization_g1), "gY")
            << FORMAT_ELEMENT(2, (model.get_g2()), "g2")
            << FORMAT_ELEMENT(3, (model.get_g3()), "g3")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("SMYu", ToMatrix(model.get_Yu()), "Yu", model.get_scale());
   slha_io.set_block("SMYd", ToMatrix(model.get_Yd()), "Yd", model.get_scale());
   slha_io.set_block("SMYe", ToMatrix(model.get_Ye()), "Ye", model.get_scale());
   {
      std::ostringstream block;
      block << "Block SMSM Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (model.get_mu2()), "mu2")
            << FORMAT_ELEMENT(2, (model.get_Lambdax()), "Lambdax")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block SMHMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(3, (model.get_v()), "v")
      ;
      slha_io.set_block(block);
   }
}

void THDMIISNMSSMBCsimple_slha_io::set_mass(const standard_model::Standard_model_physical& physical)
{
   std::ostringstream mass;

   mass << "Block SMMASS\n"
      << FORMAT_MASS(24, physical.MVWp, "VWp")
      << FORMAT_MASS(21, physical.MVG, "VG")
      << FORMAT_MASS(12, physical.MFv(0), "Fv(1)")
      << FORMAT_MASS(14, physical.MFv(1), "Fv(2)")
      << FORMAT_MASS(16, physical.MFv(2), "Fv(3)")
      << FORMAT_MASS(25, physical.Mhh, "hh")
      << FORMAT_MASS(1, physical.MFd(0), "Fd(1)")
      << FORMAT_MASS(3, physical.MFd(1), "Fd(2)")
      << FORMAT_MASS(5, physical.MFd(2), "Fd(3)")
      << FORMAT_MASS(2, physical.MFu(0), "Fu(1)")
      << FORMAT_MASS(4, physical.MFu(1), "Fu(2)")
      << FORMAT_MASS(6, physical.MFu(2), "Fu(3)")
      << FORMAT_MASS(11, physical.MFe(0), "Fe(1)")
      << FORMAT_MASS(13, physical.MFe(1), "Fe(2)")
      << FORMAT_MASS(15, physical.MFe(2), "Fe(3)")
      << FORMAT_MASS(22, physical.MVP, "VP")
      << FORMAT_MASS(23, physical.MVZ, "VZ")
      ;

   slha_io.set_block(mass);
}

void THDMIISNMSSMBCsimple_slha_io::set_mixing_matrices(const standard_model::Standard_model_physical& physical)
{
   slha_io.set_block("SMUULMIX", physical.Vu, "Vu");
   slha_io.set_block("SMUDLMIX", physical.Vd, "Vd");
   slha_io.set_block("SMUURMIX", physical.Uu, "Uu");
   slha_io.set_block("SMUDRMIX", physical.Ud, "Ud");
   slha_io.set_block("SMUELMIX", physical.Ve, "Ve");
   slha_io.set_block("SMUERMIX", physical.Ue, "Ue");
}

void THDMIISNMSSMBCsimple_slha_io::set_spectrum(const standard_model::Standard_model& model)
{
   const auto& physical = model.get_physical();

   set_model_parameters(model);
   set_mass(physical);
   set_mixing_matrices(physical);
}

/**
 * Write SLHA object to given output.  If output == "-", then the SLHA
 * object is written to std::cout.  Otherwise, output is interpreted
 * as a file name
 *
 * @param output "-" for cout, or file name
 */
void THDMIISNMSSMBCsimple_slha_io::write_to(const std::string& output) const
{
   if (output == "-")
      write_to_stream(std::cout);
   else
      write_to_file(output);
}

/**
 * Read (DR-bar) model parameter output scale from MODSEL entry 12
 */
double THDMIISNMSSMBCsimple_slha_io::get_parameter_output_scale() const
{
   return slha_io.get_modsel().parameter_output_scale;
}

/**
 * Read SLHA object from file
 *
 * @param file_name file name
 */
void THDMIISNMSSMBCsimple_slha_io::read_from_file(const std::string& file_name)
{
   slha_io.read_from_file(file_name);
}

/**
 * Read SLHA object from source
 *
 * calls SLHA_io::read_from_source()
 *
 * @param source source name
 */
void THDMIISNMSSMBCsimple_slha_io::read_from_source(const std::string& source)
{
   slha_io.read_from_source(source);
}

/**
 * Read SLHA object from stream
 *
 * @param istr stream name
 */
void THDMIISNMSSMBCsimple_slha_io::read_from_stream(std::istream& istr)
{
   slha_io.read_from_stream(istr);
}

/**
 * Fill struct of model input parameters from SLHA object (MINPAR,
 * EXTPAR and IMEXTPAR blocks)
 *
 * @param input struct of model input parameters
 */
void THDMIISNMSSMBCsimple_slha_io::fill(THDMIISNMSSMBCsimple_input_parameters& input) const
{
   SLHA_io::Tuple_processor minpar_processor = [&input, this] (int key, double value) {
      return fill_minpar_tuple(input, key, value);
   };

   SLHA_io::Tuple_processor extpar_processor = [&input, this] (int key, double value) {
      return fill_extpar_tuple(input, key, value);
   };

   SLHA_io::Tuple_processor imminpar_processor = [&input, this] (int key, double value) {
      return fill_imminpar_tuple(input, key, value);
   };

   SLHA_io::Tuple_processor imextpar_processor = [&input, this] (int key, double value) {
      return fill_imextpar_tuple(input, key, value);
   };

   slha_io.read_block("MINPAR", minpar_processor);
   slha_io.read_block("EXTPAR", extpar_processor);
   slha_io.read_block("IMMINPAR", imminpar_processor);
   slha_io.read_block("IMEXTPAR", imextpar_processor);


}

/**
 * Reads DR-bar parameters from a SLHA output file.
 *
 * @param model model class to be filled
 */
void THDMIISNMSSMBCsimple_slha_io::fill_drbar_parameters(THDMIISNMSSMBCsimple_mass_eigenstates& model) const
{
   model.set_g1(slha_io.read_entry("gauge", 1) * 1.2909944487358056);
   model.set_g2(slha_io.read_entry("gauge", 2));
   model.set_g3(slha_io.read_entry("gauge", 3));
   {
      Eigen::Matrix<double,3,3> Yu;
      slha_io.read_block("Yu", Yu);
      model.set_Yu(Yu);
   }
   {
      Eigen::Matrix<double,3,3> Yd;
      slha_io.read_block("Yd", Yd);
      model.set_Yd(Yd);
   }
   {
      Eigen::Matrix<double,3,3> Ye;
      slha_io.read_block("Ye", Ye);
      model.set_Ye(Ye);
   }
   model.set_Lambda1(slha_io.read_entry("HMIX", 31));
   model.set_Lambda2(slha_io.read_entry("HMIX", 32));
   model.set_Lambda3(slha_io.read_entry("HMIX", 33));
   model.set_Lambda4(slha_io.read_entry("HMIX", 34));
   model.set_Lambda5(slha_io.read_entry("HMIX", 35));
   model.set_Lambda6(slha_io.read_entry("HMIX", 36));
   model.set_Lambda7(slha_io.read_entry("HMIX", 37));
   model.set_Lambda8(slha_io.read_entry("HMIX", 38));
   model.set_M112(slha_io.read_entry("HMIX", 20));
   model.set_M222(slha_io.read_entry("HMIX", 21));
   model.set_M332(slha_io.read_entry("HMIX", 23));
   model.set_M123(slha_io.read_entry("HMIX", 24));
   model.set_M5(slha_io.read_entry("HMIX", 25));
   model.set_v1(slha_io.read_entry("HMIX", 102));
   model.set_v2(slha_io.read_entry("HMIX", 103));


   model.set_scale(read_scale());
}

/**
 * Reads DR-bar parameters, pole masses and mixing matrices (in
 * Haber-Kane convention) from a SLHA output file.
 *
 * @param model model class to be filled
 */
void THDMIISNMSSMBCsimple_slha_io::fill(THDMIISNMSSMBCsimple_mass_eigenstates& model) const
{
   fill_drbar_parameters(model);

   THDMIISNMSSMBCsimple_physical physical_hk;
   fill_physical(physical_hk);
   physical_hk.convert_to_hk();
   model.get_physical() = physical_hk;
}

/**
 * Fill struct of extra physical input parameters from SLHA object
 * (FlexibleSUSYInput block)
 *
 * @param input struct of physical non-SLHA input parameters
 */
void THDMIISNMSSMBCsimple_slha_io::fill(Physical_input& input) const
{
   slha_io.fill(input);
}

/**
 * Fill struct of spectrum generator settings from SLHA object
 * (FlexibleSUSY block)
 *
 * @param settings struct of spectrum generator settings to be filled
 */
void THDMIISNMSSMBCsimple_slha_io::fill(Spectrum_generator_settings& settings) const
{
   slha_io.fill(settings);
}

void THDMIISNMSSMBCsimple_slha_io::fill_minpar_tuple(THDMIISNMSSMBCsimple_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   case 1: input.lambdaNMSSM = value; break;
   case 2: input.kappaNMSSM = value; break;
   case 3: input.AlambdaNMSSM = value; break;
   case 4: input.AkappaNMSSM = value; break;
   case 5: input.AtopNMSSM = value; break;
   case 6: input.mstopL = value; break;
   case 7: input.mstopR = value; break;
   case 8: input.MEWSB = value; break;
   case 10: input.TanBeta = value; break;
   case 13: input.vSIN = value; break;
   default: WARNING("Unrecognized entry in block MINPAR: " << key); break;
   }

}

void THDMIISNMSSMBCsimple_slha_io::fill_extpar_tuple(THDMIISNMSSMBCsimple_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   default: WARNING("Unrecognized entry in block EXTPAR: " << key); break;
   }

}

void THDMIISNMSSMBCsimple_slha_io::fill_imminpar_tuple(THDMIISNMSSMBCsimple_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   default: WARNING("Unrecognized entry in block IMMINPAR: " << key); break;
   }

}

void THDMIISNMSSMBCsimple_slha_io::fill_imextpar_tuple(THDMIISNMSSMBCsimple_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   default: WARNING("Unrecognized entry in block IMEXTPAR: " << key); break;
   }

}

/**
 * Reads pole masses and mixing matrices from a SLHA output file to be filled.
 */
void THDMIISNMSSMBCsimple_slha_io::fill_physical(THDMIISNMSSMBCsimple_physical& physical) const
{
   {
      DEFINE_PHYSICAL_PARAMETER(ZH);
      slha_io.read_block("SCALARMIX", ZH);
      LOCALPHYSICAL(ZH) = ZH;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZA);
      slha_io.read_block("PSEUDOSCALARMIX", ZA);
      LOCALPHYSICAL(ZA) = ZA;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZP);
      slha_io.read_block("CHARGEMIX", ZP);
      LOCALPHYSICAL(ZP) = ZP;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Vu);
      slha_io.read_block("UULMIX", Vu);
      LOCALPHYSICAL(Vu) = Vu;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Vd);
      slha_io.read_block("UDLMIX", Vd);
      LOCALPHYSICAL(Vd) = Vd;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Uu);
      slha_io.read_block("UURMIX", Uu);
      LOCALPHYSICAL(Uu) = Uu;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Ud);
      slha_io.read_block("UDRMIX", Ud);
      LOCALPHYSICAL(Ud) = Ud;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Ve);
      slha_io.read_block("UELMIX", Ve);
      LOCALPHYSICAL(Ve) = Ve;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Ue);
      slha_io.read_block("UERMIX", Ue);
      LOCALPHYSICAL(Ue) = Ue;
   }

   LOCALPHYSICAL(MVG) = slha_io.read_entry("MASS", 21);
   LOCALPHYSICAL(MFv)(0) = slha_io.read_entry("MASS", 12);
   LOCALPHYSICAL(MFv)(1) = slha_io.read_entry("MASS", 14);
   LOCALPHYSICAL(MFv)(2) = slha_io.read_entry("MASS", 16);
   LOCALPHYSICAL(MVP) = slha_io.read_entry("MASS", 22);
   LOCALPHYSICAL(MVZ) = slha_io.read_entry("MASS", 23);
   LOCALPHYSICAL(Mhh)(0) = slha_io.read_entry("MASS", 25);
   LOCALPHYSICAL(Mhh)(1) = slha_io.read_entry("MASS", 35);
   LOCALPHYSICAL(Mhh)(2) = slha_io.read_entry("MASS", 45);
   LOCALPHYSICAL(MAh)(1) = slha_io.read_entry("MASS", 36);
   LOCALPHYSICAL(MAh)(2) = slha_io.read_entry("MASS", 46);
   LOCALPHYSICAL(MHm)(1) = slha_io.read_entry("MASS", 37);
   LOCALPHYSICAL(MFd)(0) = slha_io.read_entry("MASS", 1);
   LOCALPHYSICAL(MFd)(1) = slha_io.read_entry("MASS", 3);
   LOCALPHYSICAL(MFd)(2) = slha_io.read_entry("MASS", 5);
   LOCALPHYSICAL(MFu)(0) = slha_io.read_entry("MASS", 2);
   LOCALPHYSICAL(MFu)(1) = slha_io.read_entry("MASS", 4);
   LOCALPHYSICAL(MFu)(2) = slha_io.read_entry("MASS", 6);
   LOCALPHYSICAL(MFe)(0) = slha_io.read_entry("MASS", 11);
   LOCALPHYSICAL(MFe)(1) = slha_io.read_entry("MASS", 13);
   LOCALPHYSICAL(MFe)(2) = slha_io.read_entry("MASS", 15);
   LOCALPHYSICAL(MVWm) = slha_io.read_entry("MASS", 24);

}

/**
 * Reads the renormalization scales from all DR-bar parameter blocks.
 * If blocks with different scales are found the last scale is
 * returned and a warning is printed.
 *
 * @return common renormalization scale
 */
double THDMIISNMSSMBCsimple_slha_io::read_scale() const
{
   static const std::array<std::string, 5> drbar_blocks =
      { "gauge", "Yu", "Yd", "Ye", "HMIX" };

   double scale = 0.;

   for (const auto& block: drbar_blocks) {
      const double block_scale = slha_io.read_scale(block);
      if (!is_zero(block_scale)) {
         if (!is_zero(scale) && !is_equal(scale, block_scale))
            WARNING("DR-bar parameters defined at different scales");
         scale = block_scale;
      }
   }

   return scale;
}

} // namespace flexiblesusy
