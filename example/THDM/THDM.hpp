// ====================================================================
// This file is part of PhaseTracer

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef POTENTIAL_THDM_HPP_INCLUDED
#define POTENTIAL_THDM_HPP_INCLUDED

/*
   2HDM
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "one_loop_potential.hpp"
#include "pow.hpp"
using namespace std;


namespace EffectivePotential {

class THDM : public OneLoopPotential {
 public:
  
  void init_params(double tanb_, double m12_, double lam1_, double lam2_, double lam3_, double lam4_, double lam5_, double RG_scale_=246., double delta_ = 1e-1) {
  	tanb = tanb_;
  	m12 = m12_;
  	lam1 = lam1_;
  	lam2 = lam2_;
  	lam3 = lam3_;
  	lam4 = lam4_;
  	lam5 = lam5_;
  	lam345 = lam3 + lam4 + lam5;
  	cosb_sq = 1. / (1 + pow(tanb, 2));
  	sinb_sq = 1. - cosb_sq;
  	sin2b = 2. / (tanb + 1./tanb);
  	RG_scale = RG_scale_;	  
  	delta = delta_;
  	CT_term = get_counter_term(RG_scale, delta); 	
  }
  
  size_t get_n_scalars() const override {return 2;}
  
  vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
    auto phi1 = phi;
    phi1[0] = - phi[0];
    phi1[1] = - phi[1];
    return {phi1};
  };
//  v0-----
  double V0(Eigen::VectorXd phi) const override {
    const double h1 = phi[0];
    const double h2 = phi[1];
    
    double r = 1./2 * pow(m12, 2) * tanb * pow(h1 - h2 / tanb, 2) - vhsq / 4 * (lam1 * pow(h1, 2) + lam2 * pow(h2 * tanb, 2)) / (1 + pow(tanb, 2))
    - vhsq / 4 * lam345 * (pow(h1 * tanb, 2) +  pow(h2, 2)) / (1 + pow(tanb,2)) + 1./8 * pow(h1, 4) * lam1 + 1./8 * pow(h2, 4) * lam2 + 1./4 * pow(h1, 2) * pow(h2, 2) * lam345;
    return r;
  }
//  vct--------
  double counter_term(Eigen::VectorXd phi, double T) const override {
    const double h1 = phi[0];
    const double h2 = phi[1];
    double r = CT_term[0] * pow(h1, 2) + CT_term[1] * pow(h2, 2) + CT_term[2] * pow(h1, 4) + CT_term[3] * pow(h2, 4) + CT_term[4] * pow(h1, 2) * pow(h2, 2);
    return r;
    
  }
  
  vector<double> get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const override{
    const double h1 = phi[0];
    const double h2 = phi[1];
    
    const double T11 = pow(T, 2) / 24. * (9 * pow(Y1, 2) / 2 + 3 * pow(Y2, 2) / 2 + 6 * pow(Yb, 2) / cosb_sq + 6 * lam1 + 4 * lam3 + 2 * lam4);
    const double T22 = pow(T, 2) / 24. * (9 * pow(Y1, 2) / 2 + 3 * pow(Y2, 2) / 2 + 6 * pow(Yt, 2) / sinb_sq + 6 * lam2 + 4 * lam3 + 2 * lam4);
    // mass matrix of Higgs
    const double M11 = 3. / 2 * lam1 * pow(h1, 2) + 1. / 2 * lam345 * pow(h2, 2) + pow(m12, 2) * tanb - 1. / 2 * lam1 * vhsq * cosb_sq -
             1. / 2 * lam345 * vhsq * sinb_sq ;
    const double M22 = 3. / 2 * lam2 * pow(h2, 2) + 1. / 2 * lam345 * pow(h1, 2) + pow(m12, 2) / tanb - 1. / 2 * lam2 * vhsq * sinb_sq -
             1. / 2 * lam345 * vhsq * cosb_sq ;
    const double M12 = lam345 * h1 * h2 - pow(m12, 2);
    const double M11_T = M11 + T11;
    const double M22_T = M22 + T22;

    const double mhT_sq = 1. / 2 * (M11_T + M22_T - sqrt(pow(M11_T, 2) + 4 * pow(M12, 2) + pow(M22_T, 2) - 2 * M11_T * M22_T));
    const double mHT_sq = 1. / 2 * (M11_T + M22_T + sqrt(pow(M11_T, 2) + 4 * pow(M12, 2) + pow(M22_T, 2) - 2 * M11_T * M22_T));
    // mass matrix of A
    const double MA11 = 1. / 2 * lam1 * pow(h1, 2) + pow(m12, 2) * tanb - 1. / 2 * lam1 * vhsq * cosb_sq - 1. / 2 * lam345 * vhsq * sinb_sq + 1. / 2 * (lam3 + lam4 - lam5) * 				 pow(h2, 2);  
    const double MA22 = 1. / 2 * lam2 * pow(h2, 2) + pow(m12, 2) / tanb - 1. / 2 * lam2 * vhsq * sinb_sq - 1. / 2 * lam345 * vhsq * cosb_sq + 1. / 2 * (lam3 + lam4 - lam5) * 				 pow(h1, 2); 
    const double MA12 = lam5 * h1 * h2 - pow(m12, 2);
    const double MA11_T = MA11 + T11;
    const double MA22_T = MA22 + T22;
    
    const double mAT_sq = 1. / 2  * (MA11_T + MA22_T + sqrt(pow(MA11_T, 2) + 4 * pow(MA12, 2) + pow(MA22_T, 2) - 2 * MA11_T * MA22_T));
    const double mG0T_sq = 1. / 2  * (MA11_T + MA22_T - sqrt(pow(MA11_T, 2) + 4 * pow(MA12, 2) + pow(MA22_T, 2) - 2 * MA11_T * MA22_T));
    // mass matrix of Hpm
    const double MC11 = 1. / 2 * lam1 * pow(h1, 2) + pow(m12, 2) * tanb - 1. / 2 * lam1 * vhsq * cosb_sq - 1. / 2 * lam345 * vhsq *  sinb_sq + 1. / 2 * lam3 * pow(h2, 2);
    const double MC22 = 1. / 2 * lam2 * pow(h2, 2) + pow(m12, 2) / tanb - 1. / 2 * lam2 * vhsq * sinb_sq - 1. / 2 * lam345 * vhsq *  cosb_sq + 1. / 2 * lam3 * pow(h1, 2);
    const double MC12 = 1. / 2 * (lam4 + lam5) * h2 * h1 - pow(m12, 2);
    const double MC11_T = MC11 + T11;
    const double MC22_T = MC22 + T22;
    
    const double mHpmT_sq = 1. / 2  * (MC11_T + MC22_T + sqrt(pow(MC11_T, 2) + 4 * pow(MC12, 2) + pow(MC22_T, 2) - 2 * MC11_T * MC22_T));
    const double mGpmT_sq = 1. / 2  * (MC11_T + MC22_T - sqrt(pow(MC11_T, 2) + 4 * pow(MC12, 2) + pow(MC22_T, 2) - 2 * MC11_T * MC22_T));
    
    
    return {mhT_sq, mHT_sq, mAT_sq, mG0T_sq, mGpmT_sq, mHpmT_sq, };
  }
  
 vector<double> get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const override {
    return get_scalar_debye_sq(phi, xi, 0.);    

  }
//  ni of scalar
 vector<double> get_scalar_dofs() const override { return {1., 1., 1., 1., 2., 2.}; }

  // mass  matrix of W, Z 
 vector<double> get_vector_debye_sq(Eigen::VectorXd phi, double T) const override{
    const double h1 = phi[0];
    const double h2 = phi[1];
    const double mwLT_sq = 1. / 4 * pow(Y1, 2) * (pow(h1, 2) + pow(h2, 2)) + 2 * pow(Y1, 2) * pow(T, 2);
    const double Delta = 1. / 64 * pow((pow(Y1, 2) + pow(Y2, 2)), 2) * pow(pow(h1, 2) + pow(h2, 2) + 8 * pow(T, 2), 2) - pow(Y1 * Y2 * T, 2) * (pow(h1, 2) + pow(h2, 2) + 4 * pow(T, 2));
    const double mzLT_sq = 1. / 8 * (pow(Y1, 2) + pow(Y2, 2)) * (pow(h1, 2) + pow(h2, 2)) + (pow(Y1, 2) + pow(Y2, 2)) * pow(T, 2) + sqrt(Delta);
    const double mgamaLT_sq = 1. / 8 * (pow(Y1, 2) + pow(Y2, 2)) * (pow(h1, 2) + pow(h2, 2)) + (pow(Y1, 2) + pow(Y2, 2)) * pow(T, 2) - sqrt(Delta);
    const double mwTT_sq = 1. / 4 * pow(Y1, 2) * (pow(h1, 2) + pow(h2, 2));
    const double mzTT_sq = 1. / 4 * (pow(Y1, 2) + pow(Y2, 2)) * (pow(h1, 2) + pow(h2, 2));
    const double mgamaTT_sq = 0;

    return {mwTT_sq, mzTT_sq, mgamaTT_sq, mwLT_sq, mzLT_sq, mgamaLT_sq};
  }
  
 vector<double> get_vector_masses_sq(Eigen::VectorXd phi) const override {
    return get_vector_debye_sq(phi, 0.); 
  }

//ni of w, z, photon
 vector<double> get_vector_dofs() const override { return {4., 2., 1., 2., 1., 1.};}


// matrix of top and bottom  
 vector<double> get_fermion_masses_sq(Eigen::VectorXd phi) const override {
    return {0.5 * square(phi[1] * Yt ) / sinb_sq, 0.5 * square(phi[0] * Yb ) / cosb_sq};
  }
// ni of top and bottom
 vector<double> get_fermion_dofs() const override {
    return {12., 12.};
  }
  
 vector<double> get_tree_minimum() const {
    return {vh*sqrt(cosb_sq), vh*sqrt(sinb_sq)};
  }
 
 //function to obtain the Vcw
 vector<double> get_scalarsq_withoutgoldstone(Eigen::VectorXd phi){
    const double h1 = phi[0];
    const double h2 = phi[1];
    
    // mass matrix of Higgs
    const double M11 = 3. / 2 * lam1 * pow(h1, 2) + 1. / 2 * lam345 * pow(h2, 2) + pow(m12, 2) * tanb - 1. / 2 * lam1 * vhsq * cosb_sq -
             1. / 2 * lam345 * vhsq * sinb_sq ;
    const double M22 = 3. / 2 * lam2 * pow(h2, 2) + 1. / 2 * lam345 * pow(h1, 2) + pow(m12, 2) / tanb - 1. / 2 * lam2 * vhsq * sinb_sq -
             1. / 2 * lam345 * vhsq * cosb_sq ;
    const double M12 = lam345 * h1 * h2 - pow(m12, 2);

    const double mh_sq = 1. / 2 * (M11 + M22 - sqrt(pow(M11, 2) + 4 * pow(M12, 2) + pow(M22, 2) - 2 * M11 * M22));
    const double mH_sq = 1. / 2 * (M11 + M22 + sqrt(pow(M11, 2) + 4 * pow(M12, 2) + pow(M22, 2) - 2 * M11 * M22));
    // mass matrix of A
    const double MA11 = 1. / 2 * lam1 * pow(h1, 2) + pow(m12, 2) * tanb - 1. / 2 * lam1 * vhsq * cosb_sq - 1. / 2 * lam345 * vhsq * sinb_sq + 1. / 2 * (lam3 + lam4 - lam5) * 				 pow(h2, 2);  
    const double MA22 = 1. / 2 * lam2 * pow(h2, 2) + pow(m12, 2) / tanb - 1. / 2 * lam2 * vhsq * sinb_sq - 1. / 2 * lam345 * vhsq * cosb_sq + 1. / 2 * (lam3 + lam4 - lam5) * 				 pow(h1, 2); 
    const double MA12 = lam5 * h1 * h2 - pow(m12, 2);

    const double mA_sq = 1. / 2  * (MA11 + MA22 + sqrt(pow(MA11, 2) + 4 * pow(MA12, 2) + pow(MA22, 2) - 2 * MA11 * MA22));
    // const double mG0T_sq = 1. / 2  * (MA11_T + MA22_T - sqrt(pow(MA11_T, 2) + 4 * pow(MA12, 2) + pow(MA22_T, 2) - 2 * MA11_T * MA22_T));
    // mass matrix of Hpm
    const double MC11 = 1. / 2 * lam1 * pow(h1, 2) + pow(m12, 2) * tanb - 1. / 2 * lam1 * vhsq * cosb_sq - 1. / 2 * lam345 * vhsq *  sinb_sq + 1. / 2 * lam3 * pow(h2, 2);
    const double MC22 = 1. / 2 * lam2 * pow(h2, 2) + pow(m12, 2) / tanb - 1. / 2 * lam2 * vhsq * sinb_sq - 1. / 2 * lam345 * vhsq *  cosb_sq + 1. / 2 * lam3 * pow(h1, 2);
    const double MC12 = 1. / 2 * (lam4 + lam5) * h2 * h1 - pow(m12, 2);
    
    const double mHpm_sq = 1. / 2  * (MC11 + MC22 + sqrt(pow(MC11, 2) + 4 * pow(MC12, 2) + pow(MC22, 2) - 2 * MC11 * MC22));
    //const double mGpmT_sq = 1. / 2  * (MC11_T + MC22_T - sqrt(pow(MC11_T, 2) + 4 * pow(MC12, 2) + pow(MC22_T, 2) - 2 * MC11_T * MC22_T));
    
    
    return {mh_sq, mH_sq, mA_sq, mHpm_sq};
 
 }
 
 vector<double> get_scalar_dofs_withoutgoldstone(){ return {1., 1., 1., 2.}; }
 
 vector<double>goldstone_mass(double h1, double h2){
    const double MA11 = 1. / 2 * lam1 * pow(h1, 2) + pow(m12, 2) * tanb - 1. / 2 * lam1 * vhsq * cosb_sq - 1. / 2 * lam345 * vhsq * sinb_sq + 1. / 2 * (lam3 + lam4 - lam5) * 				 pow(h2, 2);  
    const double MA22 = 1. / 2 * lam2 * pow(h2, 2) + pow(m12, 2) / tanb - 1. / 2 * lam2 * vhsq * sinb_sq - 1. / 2 * lam345 * vhsq * cosb_sq + 1. / 2 * (lam3 + lam4 - lam5) * 				 pow(h1, 2); 
    const double MA12 = lam5 * h1 * h2 - pow(m12, 2);
    
    //const double mAT_sq = 1. / 2  * (MA11_T + MA22_T + sqrt(pow(MA11_T, 2) + 4 * pow(MA12, 2) + pow(MA22_T, 2) - 2 * MA11_T * MA22_T));
    const double mG0_sq = 1. / 2  * (MA11 + MA22 - sqrt(pow(MA11, 2) + 4 * pow(MA12, 2) + pow(MA22, 2) - 2 * MA11 * MA22));
    // mass matrix of Hpm
    const double MC11 = 1. / 2 * lam1 * pow(h1, 2) + pow(m12, 2) * tanb - 1. / 2 * lam1 * vhsq * cosb_sq - 1. / 2 * lam345 * vhsq *  sinb_sq + 1. / 2 * lam3 * pow(h2, 2);
    const double MC22 = 1. / 2 * lam2 * pow(h2, 2) + pow(m12, 2) / tanb - 1. / 2 * lam2 * vhsq * sinb_sq - 1. / 2 * lam345 * vhsq *  cosb_sq + 1. / 2 * lam3 * pow(h1, 2);
    const double MC12 = 1. / 2 * (lam4 + lam5) * h2 * h1 - pow(m12, 2);
    
    //const double mHpmT_sq = 1. / 2  * (MC11_T + MC22_T + sqrt(pow(MC11_T, 2) + 4 * pow(MC12, 2) + pow(MC22_T, 2) - 2 * MC11_T * MC22_T));
    const double mGpm_sq = 1. / 2  * (MC11 + MC22 - sqrt(pow(MC11, 2) + 4 * pow(MC12, 2) + pow(MC22, 2) - 2 * MC11 * MC22));
    return {mG0_sq, mGpm_sq};
}    
 
 vector<double>goldstone_derivative(double h1, double h2, double delta){
   vector<double> mg_h1h2 = goldstone_mass(h1, h2);
   vector<double> mg_dh1h2 = goldstone_mass(h1 + delta, h2);
   vector<double> mg_h1dh2 = goldstone_mass(h1, h2 + delta);
   double mG0_sq1 = mg_h1h2[0];
   double mGpm_sq1 = mg_h1h2[1];
   
   double mG0_sq2 = mg_dh1h2[0];
   double mGpm_sq2 = mg_dh1h2[1];
   
   double mG0_sq3 = mg_h1dh2[0];
   double mGpm_sq3 = mg_h1dh2[1];
   
   double dG0dx = (mG0_sq2 - mG0_sq1) / delta;
   double dG0dy = (mG0_sq3 - mG0_sq1) / delta;
   double dGpmdx = (mGpm_sq2 - mGpm_sq1) / delta;
   double dGpmdy = (mGpm_sq3 - mGpm_sq1) / delta;
   return {dG0dx, dG0dy, dGpmdx, dGpmdy};
 }
 
 //Vcw_withoutgoldstone
 double xlogx_withoutgoldtone(double x) {
  const double abs_x = std::abs(x);
  if (abs_x <= std::numeric_limits<double>::min()) {
    return 0.;
  } else {
    return x * std::log(abs_x);
  }
}
 
 double Vcw_withoutgoldstone(Eigen::VectorXd phi, double renormalization_scale, double xi){
  double correction = 0;
  
  const auto scalar_masses_sq = THDM::get_scalarsq_withoutgoldstone(phi);
  const auto fermion_masses_sq = THDM::get_fermion_masses_sq(phi);
  const auto vector_masses_sq = THDM::get_vector_masses_sq(phi);

  const auto scalar_dofs = THDM::get_scalar_dofs_withoutgoldstone();
  const auto fermion_dofs = THDM::get_fermion_dofs();
  const auto vector_dofs = THDM::get_vector_dofs();

  if (scalar_dofs.size() != scalar_masses_sq.size()) {
    throw std::runtime_error("Scalar dofs and masses do not match");
  }

  if (fermion_dofs.size() != fermion_masses_sq.size()) {
    throw std::runtime_error("Fermion dofs and masses do not match");
  }

  if (vector_dofs.size() != vector_masses_sq.size()) {
    throw std::runtime_error("Vector dofs and masses do not match");
  }

  // scalar correction
  for (size_t i = 0; i < scalar_masses_sq.size(); ++i) {
    const double x = scalar_masses_sq[i] / square(renormalization_scale);
    correction += scalar_dofs[i] * scalar_masses_sq[i] *
      (square(renormalization_scale) * THDM::xlogx_withoutgoldtone(x) - scalar_masses_sq[i] * 3. / 2.);
  }

  // fermion correction
  for (size_t i = 0; i < fermion_masses_sq.size(); ++i) {
    const double x = fermion_masses_sq[i] / square(renormalization_scale);
    correction -= fermion_dofs[i] * fermion_masses_sq[i] *
      (square(renormalization_scale) * THDM::xlogx_withoutgoldtone(x) - fermion_masses_sq[i] * 3. / 2.);
  }

  // vector correction
  for (size_t i = 0; i < vector_masses_sq.size(); ++i) {
    const double x = vector_masses_sq[i] / square(renormalization_scale);
    correction += vector_dofs[i] * vector_masses_sq[i] *
      (square(renormalization_scale) * THDM::xlogx_withoutgoldtone(x) - vector_masses_sq[i] * 5. / 6.);
  }

  // gauge dependent vector correction
  // hack - i know only first 3 are longitudinal in xSM model
  if (xi != 0.) {
    for (size_t i = 0; i < vector_masses_sq.size(); ++i) {
      if (i < 3) {
        const double x = xi * vector_masses_sq[i] / square(renormalization_scale);
        correction -= vector_dofs[i] * xi * vector_masses_sq[i] *
          (square(renormalization_scale) * THDM::xlogx_withoutgoldtone(x) - xi * vector_masses_sq[i] * 3. / 2.);
      }
    }
  }

  return correction / (64. * M_PI * M_PI);
}

//obtain counter term
 
 vector<double> get_counter_term(double renormalization_scale, double delta){
       double vh1 = vh * sqrt(cosb_sq);
       double vh2 = vh * sqrt(sinb_sq);
       
       Eigen::VectorXd EW_vev(2);
       EW_vev << vh1, vh2;
       
       Eigen::VectorXd EW_h1_vev(2);
       EW_h1_vev << vh1 + delta, vh2;
       
       Eigen::VectorXd EW_h1m_vev(2);
       EW_h1m_vev << vh1 - delta, vh2;
       
       Eigen::VectorXd EW_h2_vev(2);
       EW_h2_vev << vh1, vh2 + delta;
       
       Eigen::VectorXd EW_h2m_vev(2);
       EW_h2m_vev << vh1, vh2 - delta;
       
       Eigen::VectorXd EW_h1h2_vev(2);
       EW_h1h2_vev << vh1 + delta, vh2 + delta;
       
       Eigen::VectorXd EW_h1mh2m_vev(2);
       EW_h1mh2m_vev << vh1 - delta, vh2 - delta;
       
       Eigen::VectorXd EW_h1h2m_vev(2);
       EW_h1h2m_vev << vh1 + delta, vh2 - delta;
       
       Eigen::VectorXd EW_h1mh2_vev(2);
       EW_h1mh2_vev << vh1 - delta, vh2 + delta;
       
       vector<double> diff_G = THDM::goldstone_derivative(vh1, vh2, delta);
       double dG0dx = diff_G[0];
       double dG0dy = diff_G[1]; 
       double dGpmdx = diff_G[2];
       double dGpmdy = diff_G[3];
             
       Eigen::Matrix<double, 5, 5> Coe;
       Coe << 2. * vh1, 0, 4 * pow(vh1, 3), 0, 2 * vh1 * pow(vh2, 2),
              0, 2 * vh2, 0, 4 * pow(vh2, 3), 2 * vh2 * pow(vh1, 2),
              2, 0, 12 * pow(vh1, 2), 0, 2 * pow(vh2, 2),
              0, 2, 0, 12 * pow(vh2, 2), 2 * pow(vh1, 2),
              0, 0, 0, 0, 4 * vh1 * vh2;
       
       double xi = 0;
       double b1 = -1 * (THDM::Vcw_withoutgoldstone(EW_h1_vev, RG_scale, xi) - THDM::Vcw_withoutgoldstone(EW_vev, RG_scale, xi)) / delta;
       double b2 = -1 * (THDM::Vcw_withoutgoldstone(EW_h2_vev, RG_scale, xi) - THDM::Vcw_withoutgoldstone(EW_vev, RG_scale, xi)) / delta;
       double b3 = -1 * (THDM::Vcw_withoutgoldstone(EW_h1_vev, RG_scale, xi) - 2 * THDM::Vcw_withoutgoldstone(EW_vev, RG_scale, xi) + THDM::Vcw_withoutgoldstone(EW_h1m_vev, RG_scale, xi)) / pow(delta, 2) - 1. / (32. * M_PI * M_PI) * std::log(125. * 125. / vhsq) * (dG0dx * dG0dx + 2 * dGpmdx * dGpmdx);	
       double b4 = -1 * (THDM::Vcw_withoutgoldstone(EW_h2_vev, RG_scale, xi) - 2 * THDM::Vcw_withoutgoldstone(EW_vev, RG_scale, xi) + THDM::Vcw_withoutgoldstone(EW_h2m_vev, RG_scale, xi)) / pow(delta, 2) - 1. / (32. * M_PI * M_PI) * std::log(125. * 125. / vhsq) * (dG0dy * dG0dy + 2 * dGpmdy * dGpmdy);
       double b5 = -1 * (THDM::Vcw_withoutgoldstone(EW_h1h2_vev, RG_scale, xi) - THDM::Vcw_withoutgoldstone(EW_h1h2m_vev, RG_scale, xi) - THDM::Vcw_withoutgoldstone(EW_h1mh2_vev, RG_scale, xi)+ THDM::Vcw_withoutgoldstone(EW_h1mh2m_vev, RG_scale, xi)) / pow(2 * delta, 2) - 1. / (32. * M_PI * M_PI) * std::log(125. * 125. / vhsq) * (dG0dx * dG0dy + 2 * dGpmdx * dGpmdy);
       
       Eigen::Matrix<double, 5, 1> b;
       b <<b1, b2, b3, b4, b5;
       Eigen::Matrix<double, 5, 1> solution = Coe.colPivHouseholderQr().solve(b);
       std::vector<double> counter_term={solution(0,0), solution(1,0), solution(2,0), solution(3,0), solution(4,0)};   
       return counter_term;                     
}
 
  
  

  
 private:
  
  const double vh = 246.0;
  const double vhsq = vh * vh;
  const double alpha = sqrt(4.00 * M_PI/128.0);
  const double sin2w = 0.2315;
  const double sw = sqrt(sin2w);
  const double cw = sqrt(1.0 - sin2w);
  const double Y1 = alpha / sw;
  const double Y2 = alpha / cw;
  const double Yt = sqrt(2.0) * 172.5 / vh;
  const double Yb = sqrt(2.0) * 4.75 / vh;
  
  double tanb;
  double cosb_sq;
  double sinb_sq;
  double sin2b;
  double RG_scale;
  double delta;
  double m12;
  double lam1;
  double lam2;
  double lam3;
  double lam4;
  double lam5;
  double lam345;
  std::vector<double> CT_term;
  

};

}  // namespace EffectivePotential

#endif
