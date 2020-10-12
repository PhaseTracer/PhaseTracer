#ifndef PHASETRACER_SDM_MODEL_HPP_INCLUDED
#define PHASETRACER_SDM_MODEL_HPP_INCLUDED

/**
  SingletDM
*/

#include <cmath>
#include <vector>

#include "one_loop_potential.hpp"
#include "pow.hpp"

#include "SingletDM_input_parameters.hpp"
#include "SingletDM_slha_io.hpp"
#include "SingletDM_spectrum_generator.hpp"
#include "SingletDM_two_scale_spectrum_generator.hpp"


namespace fs = flexiblesusy;
typedef fs::SingletDM_spectrum_generator<fs::Two_scale> SpectrumGenerator;
typedef fs::SingletDM_slha<fs::SingletDM<fs::Two_scale>> Model;
typedef fs::SingletDM_input_parameters Input;
typedef fs::Spectrum_generator_settings Settings;


namespace EffectivePotential {

class SingletDM : public OneLoopPotential {
 public:
  SingletDM() {
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
  std::vector<double> get_vector_dofs() const override {
    return {6., 3., 3.};
  }

  // 3 CP-even, 3 CP-odd and 2 charged Higgs
  std::vector<double> get_scalar_dofs() const override {
    return {1., 1., 3.};
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
  const double mt = 173.1;
  double m1, m2, m3, m4, m5;
  double l1, l2, l3, l4, l5, l6, l7, l8;
  double g1, g2;
  double gp, g;
  double c_h, c_s;
};

void SingletDM::set_input(std::vector<double> x) {
  Input input;
  input.Qin = x[0];
  input.QEWSB = x[1];
  input.HiggsIN = x[2];
  input.muSInput = x[3];
  input.LamSInput = x[4];
  input.LamSHInput = x[5];
//  input.K2input = x[6];
//  input.vSInput = x[7];
//  x[8];

  // Arrange settings
  Settings settings;
  settings.set(Settings::precision, 1.e-4);
  settings.set(Settings::calculate_sm_masses, 0);
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

//  m1 = model.get_M112();
//  m2 = model.get_M222();
//  m3 = model.get_M332();
//  m4 = model.get_M123();
//  m5 = model.get_M5();
//
//  l1 = model.get_Lambda1();
//  l2 = model.get_Lambda2();
//  l3 = model.get_Lambda3() - model.get_Lambda4();
//  l4 = model.get_Lambda4();
//  l5 = model.get_Lambda5();
//  l6 = model.get_Lambda6();
//  l7 = model.get_Lambda7();
//  l8 = model.get_Lambda8();

  gp = model.get_g1() * sqrt(3. / 5.);
  g = model.get_g2();
  
  set_renormalization_scale(model.get_scale());
  
  
  // Calculate Debye coefficients

  const double yt = model.get_Yu(2, 2);
  const double yb = model.get_Yd(2, 2);
  const double ytau = model.get_Ye(2, 2);

  c_h = 1. / 48. * (3. * square(gp)
                     + 9. * square(g)
                     + 12. * square(yt));

  c_s = 0;
}

double SingletDM::V0(Eigen::VectorXd phi) const {
  const double h = phi[0] ;
  const double S = phi[1] ;

  return 0.5;
}

std::vector<double> SingletDM::get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const {
//  auto model_copy = model;
//  model_copy.solve_ewsb_tree_level();
//  model_copy.set_v(phi[0]);
//  model_copy.set_vS(phi[1]);
//
//  // CP-even Higgs in \xi = 1
//  model_copy.calculate_Mhh();
//  Eigen::Matrix<double, 3, 3> matrix_hh = model_copy.get_mass_matrix_hh();
//  Eigen::Array<double, 3, 1> hh_debye;
//  Eigen::Matrix<double, 3, 3> Zh;
//  matrix_hh(0, 0) += c_hd * square(T);
//  matrix_hh(1, 1) += c_hu * square(T);
//  matrix_hh(2, 2) += c_s * square(T);
//  fs::fs_diagonalize_hermitian(matrix_hh, hh_debye, Zh);
//
//  // CP-odd Higgs in \xi = 1
//  model_copy.calculate_MAh();
//  Eigen::Matrix<double, 3, 3> matrix_Ah = model_copy.get_mass_matrix_Ah();
//  Eigen::Array<double, 3, 1> Ah_debye;
//  Eigen::Matrix<double, 3, 3> ZA;
//  matrix_Ah(0, 0) += c_hd * square(T);
//  matrix_Ah(1, 1) += c_hu * square(T);
//  matrix_Ah(2, 2) += c_s * square(T);
//  fs::fs_diagonalize_hermitian(matrix_Ah, Ah_debye, ZA);
//
//  // Charged Higgs in \xi = 1
//  model_copy.calculate_MHm();
//  Eigen::Matrix<double, 2, 2> matrix_Hm = model_copy.get_mass_matrix_Hm();
//  Eigen::Array<double, 2, 1> Hm_debye;
//  Eigen::Matrix<double, 2, 2> ZP;
//  matrix_Hm(0, 0) += c_hd * square(T);
//  matrix_Hm(1, 1) += c_hu * square(T);
//  fs::fs_diagonalize_hermitian(matrix_Hm, Hm_debye, ZP);

  std::vector<double> scalar_debye_sq;
  scalar_debye_sq.push_back(0);
  scalar_debye_sq.push_back(0);
  scalar_debye_sq.push_back(0);
  return scalar_debye_sq;
}

std::vector<double> SingletDM::get_vector_debye_sq(Eigen::VectorXd phi, double T) const {
  const double vev_sq = square(phi[0]) + square(phi[1]);

  const double A = (square(g) + square(gp)) * (square(T) + 0.125 * vev_sq);
  const double B = 0.125 * sqrt(square(square(g) - square(gp)) *
                               (64. * pow_4(T) + 16. * square(T) * vev_sq) +
                               square(square(g) + square(gp)) * square(vev_sq));

  const double W_debye = square(g) * (0.25 * vev_sq + 2. * square(T));
  const double Z_debye = A + B;
  const double g_debye = A - B;

  return {W_debye, Z_debye, g_debye};
}

std::vector<double> SingletDM::get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const {
  return get_scalar_debye_sq(phi, xi, 0.);
}

std::vector<double> SingletDM::get_vector_masses_sq(Eigen::VectorXd phi) const {
  return get_vector_debye_sq(phi, 0.);
}

std::vector<double> SingletDM::get_fermion_masses_sq(Eigen::VectorXd phi) const {
  auto model_copy = model;
  model_copy.solve_ewsb_tree_level();
  model_copy.set_v(phi[0]);
//  model_copy.set_vs(phi[1]);

  model_copy.calculate_MFu();
  model_copy.calculate_MFd();
  model_copy.calculate_MFe();

  return {square(model_copy.get_MFu()[2]), square(model_copy.get_MFd()[2]), square(model_copy.get_MFe()[2])};
}

}  // namespace EffectivePotential

#endif
