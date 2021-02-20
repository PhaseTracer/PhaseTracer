#ifndef PHASETRACER_SDMEFTHIGGS_MODEL_HPP_INCLUDED
#define PHASETRACER_SDMEFTHIGGS_MODEL_HPP_INCLUDED

/**
  ScalarSingletZ2DMEFTHiggs
*/

#include <cmath>
#include <vector>

#include "one_loop_potential.hpp"
#include "pow.hpp"

#include "ScalarSingletZ2DMEFTHiggs_input_parameters.hpp"
#include "ScalarSingletZ2DMEFTHiggs_slha_io.hpp"
#include "ScalarSingletZ2DMEFTHiggs_spectrum_generator.hpp"
#include "ScalarSingletZ2DMEFTHiggs_two_scale_spectrum_generator.hpp"

namespace fs = flexiblesusy;
typedef fs::ScalarSingletZ2DMEFTHiggs_spectrum_generator<fs::Two_scale> SpectrumGenerator;
typedef fs::ScalarSingletZ2DMEFTHiggs_slha Model;
typedef fs::ScalarSingletZ2DMEFTHiggs_input_parameters Input;
typedef fs::Spectrum_generator_settings Settings;


namespace EffectivePotential {

class ScalarSingletZ2DMEFTHiggs : public OneLoopPotential {
 public:
  ScalarSingletZ2DMEFTHiggs() {
    set_xi(1.); // TODO better way to set this?
  }

  void set_input(std::vector<double> x);

  double V0(Eigen::VectorXd phi) const override;
  std::vector<double> get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const override;
  std::vector<double> get_vector_masses_sq(Eigen::VectorXd phi) const override;
  std::vector<double> get_fermion_masses_sq(Eigen::VectorXd phi) const override;
  std::vector<double> get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const override;
  std::vector<double> get_vector_debye_sq(Eigen::VectorXd phi, double T) const override;

  // W, Z and photon
  // should the photon (3rd entry( be 3 or 2?
  std::vector<double> get_vector_dofs() const override {
    return {6., 3., 3.};
  }

  // 1 complex doublet plus 1 real scalar singlet = 5
  // --> Physical Higgs, Singlet, G^0, G^\pm 
  std::vector<double> get_scalar_dofs() const override {
    return {1., 1., 1., 2.};
  }

  // top, bottom and tau
  std::vector<double> get_fermion_dofs() const override {
    return {12., 12., 4.};
  }

  size_t get_n_scalars() const override { return 2; }
  
  std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
    phi[0] = - phi[0];
    phi[1] = - phi[1];
    return {phi};
  };
    
  Eigen::Array<double, 2, 1> get_mh() { return {model.get_Mhh_pole_slha(), model.get_Mss_pole_slha()}; }

 private:
  Model model;
  //PA: do we really want mt to be const?  why?
  const double mt = 173.1;
  // Higss potential parameters
  double lambda_h, lambda_s, lambda_hs, muH2, muS2;
  // PA: what about the Higgs VEV do we need that here?
  // PA: seems odd to have gauge couplings here, but not Yukawas
  double gp, g;
  double c_h, c_s;
};

void ScalarSingletZ2DMEFTHiggs::set_input(std::vector<double> x) {
  Input input;
  input.Qin = x[0];
  input.QEWSB = x[1];
  input.LamHInput = x[2];
  input.muS2Input = x[3];
  input.LamSInput = x[4];
  input.LamSHInput = x[5];

  // Arrange settings
  Settings settings;
  settings.set(Settings::precision, 1.e-4);
  settings.set(Settings::calculate_sm_masses, 1);
  SpectrumGenerator spectrum_generator;
  spectrum_generator.set_settings(settings);

  // Run FS spectrum generator
  softsusy::QedQcd qedqcd;
  qedqcd.setPoleMt(mt);
  spectrum_generator.run(qedqcd, input);

  // Fetch model from FS
  model = std::get<0>(spectrum_generator.get_models_slha());

  // Check that model is physical
  if (spectrum_generator.get_problems().have_problem()) {
    throw std::runtime_error("Unphysical spectrum");
  }

  // Fetch parameters from FS model
  gp = model.get_g1() * sqrt(3. / 5.);
  g = model.get_g2();
  
  muH2 = model.get_muH2();
  muS2 = model.get_muS2();
  lambda_h = model.get_LamH();
  lambda_s = model.get_LamS();
  lambda_hs = model.get_LamSH();
  
  std::cout << "muH2=" << muH2 << std::endl;
  std::cout << "lambda_h=" << lambda_h << std::endl;
  std::cout << "muS2=" << muS2 << std::endl;
  std::cout << "lambda_s=" << lambda_s << std::endl;
  std::cout << "lambda_hs=" << lambda_hs << std::endl;
  std::cout << "pole mass Mh = "  << model.get_Mhh_pole_slha() << std::endl;
  std::cout << "running mass mh = "  << model.get_Mhh() << std::endl;
  std::cout << "pole mass Ms = "  << model.get_Mss_pole_slha() << std::endl;
  std::cout << "running mass ms = "  << model.get_Mss() << std::endl;
  set_renormalization_scale(model.get_scale());
  
  
  // Calculate Debye coefficients
  //PA: Do we include all 3 3rd gen fermions /  xSM_MSbar does not 
  const double yt = model.get_Yu(2, 2);
  const double yb = model.get_Yd(2, 2);
  const double ytau = model.get_Ye(2, 2);

  // PA: I took these from Yang's code already and they match Lachlan's
  // PA: up to factor 4 differences for quartics presumably from coupling defs
  // PA: but i should independently check them. 
  c_h = 1. / 48. *  ( 9. * square(g) + 3. * square(gp)
		      + 12. * square(yt) + 4. * square(yb) + 4. * square(ytau)
		      + 24. * lambda_h + 2. * lambda_hs );

  c_s =  (2. * lambda_hs + 3. * lambda_s) / 12.;
}

double ScalarSingletZ2DMEFTHiggs::V0(Eigen::VectorXd phi) const {
  //PA: Do we need  * M_SQRT1_2 for the H field like in THDMIISNMSSMBC?
  //PA: I guess it depends how I write the poetntial 
  const double h = phi[0] ;
  const double s = phi[1] ;

  const double V0 =
    0.5 * muH2 * square(h) + 0.25 * lambda_h * pow_4(h) +
    0.25 * lambda_hs * square(h) * square (s) +
    0.5 * muS2 * square(s) + 0.25 * lambda_s * pow_4(s);
  return V0;
}

std::vector<double> ScalarSingletZ2DMEFTHiggs::get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const {
  auto model_copy = model;
  // we use tree-levekl EWSB conditions for tree-level Higgs masses
  // There amy some complications in this for xSM_MSbar and comparison
  // ppaer work
  model_copy.solve_ewsb_tree_level();
  model_copy.set_v(phi[0]);

  // CP-even Higgs in \xi = 1
  //PA:  I don't see the need fpor calculating the mass at this step.
  model_copy.calculate_Mhh();
  // get mass matrix for Higgs - in this case it is just the tree mass^2
  double matrix_hh = model_copy.get_mass_matrix_hh();
  double matrix_ss = model_copy.get_mass_matrix_ss();
  // add debeye masses
  matrix_hh = matrix_hh  += c_h * square(T);
  matrix_ss = matrix_ss  += c_s * square(T);  

  // G^0 golstone mass, Do we need goldstones, I assume so
  // PA: as above I don't see the need fpor calculating the mass at this step.
  model_copy.calculate_MAh();
  double matrix_Ah = model_copy.get_mass_matrix_Ah();
  // PA: double check it is correct to just add same as hh case
  // when auto-generating this we probbaly define them all debeye coeffs
  // indepenedently even though somme are the same
  matrix_Ah += c_h * square(T);
  
  // G^+- golstone mass, Do we need goldstones, I assume so
  // PA: as above I don't see the need fpor calculating the mass at this step.
  model_copy.calculate_MHp();
  double matrix_Hp = model_copy.get_mass_matrix_Hp();
  // PA: double check it is correct to just add same as hh case
  // when auto-generating this we probbaly define them all debeye coeffs
  // indepenedently even though somme are the same
  matrix_Hp += c_h * square(T);
  
  std::vector<double> scalar_debye_sq;
  scalar_debye_sq.push_back(matrix_hh);
  scalar_debye_sq.push_back(matrix_ss);
  scalar_debye_sq.push_back(matrix_Ah);
  scalar_debye_sq.push_back(matrix_Hp);

  return scalar_debye_sq;
}

std::vector<double> ScalarSingletZ2DMEFTHiggs::get_vector_debye_sq(Eigen::VectorXd phi, double T) const {
  
  const double h_sq = square(phi[0]);
  const double T2 = T * T;
  //squared gauge couplings
  const double g_sq = square(g);
  const double gp_sq = square(gp); 
  const double MW_sq = 0.25 * square(g) * h_sq;
  // Use Debeye coefficients fromlachlan for now, these are different to thdm,
  // TODO: cross check literature and understand difference.
  const double MW_sq_T = MW_sq + 11.0/6.0 * square(g) * square(T);
  /// Z and photon thermal corrections come from diagonalsing 2 by 2 
  const double a = ( g_sq + gp_sq ) * ( 3 * h_sq + 22 * T2 );
  const double b = std::sqrt( 9.0 * square(g_sq + gp_sq) * square(h_sq)
			     + 132.0 * square(g_sq - gp_sq) * h_sq * T2
			     + 484.0 * square(g_sq - gp_sq) * T2 * T2);
  
  const double mZSq_T = (a + b) / 24.0;
  const double mPhotonSq_T = (a - b) / 24.0;
  // const double MZ_sq = 0.25 * (square(g) + square(gp)) * h_sq;
  // const double MG_sq = 0;
  
  // const double vev_sq = square(phi[0]) + square(phi[1]);

  // const double A = (square(g) + square(gp)) * (square(T) + 0.125 * vev_sq);
  // const double B = 0.125 * sqrt(square(square(g) - square(gp)) *
  //                              (64. * pow_4(T) + 16. * square(T) * vev_sq) +
  //                              square(square(g) + square(gp)) * square(vev_sq));

  // const double W_debye = square(g) * (0.25 * vev_sq + 2. * square(T));
  // const double Z_debye = A + B;
  // const double g_debye = A - B;

  return {MW_sq_T, mZSq_T,  mPhotonSq_T};
}

std::vector<double> ScalarSingletZ2DMEFTHiggs::get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const {
  return get_scalar_debye_sq(phi, xi, 0.);
}

std::vector<double> ScalarSingletZ2DMEFTHiggs::get_vector_masses_sq(Eigen::VectorXd phi) const {
  return get_vector_debye_sq(phi, 0.);
}

std::vector<double> ScalarSingletZ2DMEFTHiggs::get_fermion_masses_sq(Eigen::VectorXd phi) const {
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