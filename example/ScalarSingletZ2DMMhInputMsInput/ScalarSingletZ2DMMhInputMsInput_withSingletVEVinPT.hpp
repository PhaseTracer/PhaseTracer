#ifndef PHASETRACER_SDMMHINPUTMSINPUT_MODEL_HPP_INCLUDED
#define PHASETRACER_SDMMHINPUTMSINPUT_MODEL_HPP_INCLUDED

/**
  ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT

 */

#include <cmath>
#include <vector>
#include <iomanip>
#include <Eigen/Eigenvalues>

#include "one_loop_potential.hpp"
#include "pow.hpp"

#include "ScalarSingletZ2DMMhInputMsInput_input_parameters.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_slha_io.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_spectrum_generator.hpp"
#include "ScalarSingletZ2DMMhInputMsInput_two_scale_spectrum_generator.hpp"


namespace fs = flexiblesusy;
typedef fs::ScalarSingletZ2DMMhInputMsInput_spectrum_generator<fs::Two_scale> SpectrumGenerator;
typedef fs::ScalarSingletZ2DMMhInputMsInput_slha Model;
typedef fs::ScalarSingletZ2DMMhInputMsInput_input_parameters Input;
typedef fs::Spectrum_generator_settings Settings;


namespace EffectivePotential {

class ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT : public OneLoopPotential {
 public:
  ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT() {
    set_xi(1.); // TODO better way to set this?
  }

  void set_input(std::vector<double> x);

  // set data members of this class from FS Model object
  void set_VH_pars_from_FS(double renormalization_scale) {
    // Fetch parameters from FS model
    gp = model.get_g1() * sqrt(3. / 5.);
    g = model.get_g2();
  
    muH2 = model.get_muH2();
    muS2 = model.get_muS2();
    lambda_h = model.get_LamH() *0.5;
    lambda_s = model.get_LamS();
    lambda_hs = model.get_LamSH();
    
    if (debug){
      std::cout << "muH2=" << muH2 << std::endl;
      std::cout << "lambda_h=" << lambda_h << std::endl;
      std::cout << "muS2=" << muS2 << std::endl;
      std::cout << "lambda_s=" << lambda_s << std::endl;
      std::cout << "lambda_hs=" << lambda_hs << std::endl;
      std::cout << "pole mass Mh = "  << model.get_Mhh_pole_slha() << std::endl;
      std::cout << "running mass mh = "  << model.get_Mhh() << std::endl;
      std::cout << "pole mass Ms = "  << model.get_Mss_pole_slha() << std::endl;
      std::cout << "running mass ms = "  << model.get_Mss() << std::endl;
      std::cout << "running mass for MZ = " << model.get_MVZ() << std::endl;
      std::cout << "running mass for MW = " << model.get_MVWp() << std::endl;
    }
								
    set_renormalization_scale(renormalization_scale);
//    std::cout << "renormalization_scale= " << renormalization_scale << std::endl;
    // Calculate Debye coefficients
    yt = model.get_Yu(2, 2);
    yb = model.get_Yd(2, 2);
    ytau = model.get_Ye(2, 2);
    
  }
 
  std::vector<double> get_scalar_thermal_sq(double T) const override {
    const double c_h = (9. * square(g) +
                        3. * square(gp) +
                        2. * (6. * square(yt) + 6. * square(yb) +
                              2. * square(ytau) + 12. * lambda_h + lambda_hs)) / 48.;
    const double c_s = (2. * lambda_hs + 3. * lambda_s) / 12.;
    return {c_h * square(T), c_s * square(T)};
  }
  
  double V0(Eigen::VectorXd phi) const override;
  std::vector<double> get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const override;
  std::vector<double> get_vector_masses_sq(Eigen::VectorXd phi) const override;
  std::vector<double> get_fermion_masses_sq(Eigen::VectorXd phi) const override;
  std::vector<double> get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const override;
  std::vector<double> get_vector_debye_sq(Eigen::VectorXd phi, double T) const override;

  // W, Z and photon
  std::vector<double> get_vector_dofs() const override {
    return {2., 1., 1., 4., 2., 2.};
  }

  // 1 complex doublet plus 1 real scalar singlet = 5
  // --> Physical Higgs, Singlet, G^0, G^\pm 
  std::vector<double> get_scalar_dofs() const override {
    return {1., 1., 1., 1., 1.};
  }

  // top, bottom and tau
  std::vector<double> get_fermion_dofs() const override {
    return {12., 12., 4.};
  }

  size_t get_n_scalars() const override { return 2; }
  
  std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
    auto phi1 = phi;
    phi1[0] = - phi[0];
    auto phi2 = phi;
    phi2[1] = - phi[1];
    return {phi1, phi2};
  };

  // Allows running to a different scale, for checking scale dependence
  // TODO;  do we leave this as a public method, this changes the
  // renormalsiation scale permanently
  // Could instead make a copy then run and return that
  // and/or only have a public methiod that doe sthe whole thing of
  // varying the scale and recaculating the effective piotential 
  void Run_pars_to(double scale, double tol) {
    model.run_to(scale,tol);
    set_VH_pars_from_FS(scale); // includes setting renormalisation scale
  }

  // Whether to use special tadpole constraints in masses entering Coleman-Weinberg potential 
  void set_use_1L_EWSB_in_0L_mass(bool use_1L_EWSB_in_0L_mass_) { use_1L_EWSB_in_0L_mass = use_1L_EWSB_in_0L_mass_; } 
  void set_use_Goldstone_resum(bool use_Goldstone_resum_) { use_Goldstone_resum = use_Goldstone_resum_; } 
  void set_use_tree_level_tadpoles(bool use_tree_level_tadpoles_) { use_tree_level_tadpoles = use_tree_level_tadpoles_; } 
  void set_use_tree_level_beta(bool use_tree_level_beta_) {use_tree_level_beta = use_tree_level_beta_;}
  
  // For debuging
  void set_debug(bool debug_) { debug=debug_; }
  
    /// this code is just for testing masses delete and replace with proper
  /// tests if possible when done
  double get_EW_VEV() {return model.get_v();}
  double get_g() {return g;}
  double get_gp() {return gp;}
  double get_yt() {return yt;}
  double get_yb() {return yb;}
  double get_ytau() {return ytau;}
  
  double get_ms() {return model.get_Mss();}
  double get_muh_sq() {return muH2;}
  double get_mus_sq() {return muS2;}
  double get_lambda_h() {return lambda_h;}
  double get_lambda_s() {return lambda_s;}
  double get_lambda_hs() {return lambda_hs;}
  
  double get_v_tree_s() const {
    if (muS2>0)
      return 0;
    else
      return std::sqrt(-muS2 / lambda_s);
  }
  
 private:
  Model model;
  //TODO: Make it an input
  double mt_pole = 173;
  // paramters in tree level Higgs potential
  double lambda_h, lambda_s, lambda_hs, muH2, muS2;
  // other parameters that enter at the loop level 
  double gp, g, yt, ytau, yb;

  // flag for using tree-level or one=-loop EWSB conditions in tree-level masses
  bool use_1L_EWSB_in_0L_mass{false};
  bool use_Goldstone_resum{false};
  bool use_tree_level_tadpoles{false};
  bool use_tree_level_beta{false};

  // For debuging
  bool debug = false;
};

void ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT::set_input(std::vector<double> x) {
  Input input;
  input.Qin = x[0];
  input.QEWSB = x[1];
  input.MhInput = x[2];
  input.MsInput = x[3];
  input.LamSInput = x[4];
  input.LamSHInput = x[5];

  // Arrange settings
  Settings settings;
  settings.set(Settings::precision, 1.e-8);
  settings.set(Settings::calculate_sm_masses, 1);
  settings.set(Settings::loop_library, 0);
  if (use_tree_level_tadpoles) {
    settings.set(Settings::pole_mass_loop_order, 0);
    settings.set(Settings::ewsb_loop_order, 0);
  } else {
    settings.set(Settings::pole_mass_loop_order, 1);
    settings.set(Settings::ewsb_loop_order, 1);
  }

  if (use_tree_level_beta) {
    settings.set(Settings::threshold_corrections_loop_order,0);
    settings.set(Settings::beta_loop_order, 0);
  }
  
  SpectrumGenerator spectrum_generator;
  spectrum_generator.set_settings(settings);

  // Run FS spectrum generator
  softsusy::QedQcd qedqcd;
  qedqcd.setPoleMt(mt_pole);
  spectrum_generator.run(qedqcd, input);

  // Fetch model from FS
  model = std::get<0>(spectrum_generator.get_models_slha());

  // Check that model is physical
  if (spectrum_generator.get_problems().have_problem()) {
    throw std::runtime_error("Unphysical spectrum");
  }

  // Fetch parameters from FS model anduse them to  set data members for this class
  set_VH_pars_from_FS(model.get_scale()); // includes setting renormalisation scale
}

double ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT::V0(Eigen::VectorXd phi) const {
  //Note: \phi_h has \sqrt{1/2} pulled out
  const double h = phi[0] ;
  const double s = phi[1] ;
  const double V0 =
    0.5 * muH2 * square(h) +
    0.25 * lambda_h * pow_4(h) +
    0.25 * lambda_hs * square(h) * square (s) +
    0.5 * muS2 * square(s) +
    0.25 * lambda_s * pow_4(s);
  return V0;
}

std::vector<double> ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT::get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const {
    const double h = phi[0];
    const double s = phi[1];
    const auto thermal_sq = get_scalar_thermal_sq(T);
    
    // Copy so that we can set tree level EWSB needed in mhh2 and mgg2
    auto model_copy = model;
    // PA: Should I use _custom version here and earlier actually?  Check FS model_file and code details. 
    model_copy.solve_ewsb_tree_level();
    const double muh_sq_use_0L_EWSB = model_copy.get_muH2();
    
    const double mhh2 = (use_1L_EWSB_in_0L_mass ? muH2 : muh_sq_use_0L_EWSB) + 3. * lambda_h * square(h) + 0.5 * lambda_hs * square(s);
    const double mgg2 = (use_1L_EWSB_in_0L_mass ? muH2 : muh_sq_use_0L_EWSB) + lambda_h * square(h) + 0.5 * lambda_hs * square(s);
    const double mss2 = muS2 + 3. * lambda_s * square(s) + 0.5 * lambda_hs * square(h);
    
    // resummed Goldstone contributions
    const auto fm2 = get_fermion_masses_sq(phi);
    const auto vm2 = get_vector_masses_sq(phi);
    
    const double Qsq = square( get_renormalization_scale() );
    const double sum = 1. / (16. * M_PI * M_PI) * (
                       3.  * lambda_h * (Qsq*xlogx(mhh2/Qsq) - mhh2)
                      +0.5 * lambda_hs * (Qsq*xlogx(mss2/Qsq) - mss2)
		                  -6.  * square(yt) * (Qsq*xlogx(fm2[0]/Qsq) - fm2[0])
                      -6.  * square(yb) * (Qsq*xlogx(fm2[1]/Qsq) - fm2[1])
                      -2.  * square(ytau) * (Qsq*xlogx(fm2[2]/Qsq) - fm2[2])
                      +1.5 * square(g) * (Qsq*xlogx(vm2[0]/Qsq) - 1./3.*vm2[0])
                      +0.75* (square(g)+square(gp)) * (Qsq*xlogx(vm2[1]/Qsq) - 1./3.*vm2[1])
                      );
             
    // Goldstone finite temperature masses
    double mTG02 =   mgg2 + thermal_sq[0] + ( ( not use_1L_EWSB_in_0L_mass and use_Goldstone_resum )? sum : 0 );
    double mTGpm2 = mTG02; // 2 degrees of freedom or two degenerate copies
    // xi-dependence
    mTG02 += 0.25 * xi * (square(g * h) + square(gp * h));
    mTGpm2 += 0.25 * xi * square(g * h);
    // CP even Higgs thermal temperature masses
    Eigen::MatrixXd MTH2 = Eigen::MatrixXd::Zero(2, 2); 
    MTH2(0,0) = mhh2 + thermal_sq[0];
    MTH2(1,1) = mss2 + thermal_sq[1];
    // Mixing between Higgs and singlet
    MTH2(0, 1) = MTH2(1, 0) = lambda_hs * h * s;
    
    // get eigenvalues
    const Eigen::VectorXd mH_sq = MTH2.eigenvalues().real();
    
    // vector for all scalars, including two mass degenerate charged goldstones
    std::vector<double> m_sq_vector{mH_sq(0), mH_sq(1), mTG02, mTGpm2, mTGpm2};
    // mass order
    std::sort(m_sq_vector.begin(), m_sq_vector.end());
    return m_sq_vector;
  } 


std::vector<double> ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT::get_vector_debye_sq(Eigen::VectorXd phi, double T) const {
  //squared gauge couplings
  const double g_sq = square(g);
  const double gp_sq = square(gp);   
  const double h_sq = square(phi[0]);
  const double T_sq = T * T;
  const double MW_T_sq = 0.25 * g_sq * h_sq;
  const double MZ_T_sq = 0.25 * (g_sq + gp_sq) * h_sq;
  const double Mphoton_T_sq = 0.;


  const double MW_L_sq = 0.25 * g_sq * h_sq + 11. / 6. * g_sq * T_sq;
	const double a_L = (g_sq + gp_sq) * (3. * h_sq + 22. * T_sq);
	const double b_L = std::sqrt(9. * square(g_sq + gp_sq) * square(h_sq)
                     + 132. * square(g_sq - gp_sq) * h_sq * T_sq
                     + 484. * square(g_sq - gp_sq) * T_sq * T_sq);
	const double MZ_L_sq = (a_L + b_L) / 24.;
	const double Mphoton_L_sq = (a_L - b_L) / 24.;

  return {MW_L_sq, MZ_L_sq, Mphoton_L_sq, MW_T_sq, MZ_T_sq, Mphoton_T_sq};
}

std::vector<double> ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT::get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const {
  return get_scalar_debye_sq(phi, xi, 0.);
}

std::vector<double> ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT::get_vector_masses_sq(Eigen::VectorXd phi) const {
  return get_vector_debye_sq(phi, 0.);
}

std::vector<double> ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT::get_fermion_masses_sq(Eigen::VectorXd phi) const {
  auto model_copy = model;
  model_copy.solve_ewsb_tree_level();
  model_copy.set_v(phi[0]);

  model_copy.calculate_MFu();
  model_copy.calculate_MFd();
  model_copy.calculate_MFe();

  return {square(model_copy.get_MFu()[2]), square(model_copy.get_MFd()[2]), square(model_copy.get_MFe()[2])};
}

}  // namespace EffectivePotential

#endif
