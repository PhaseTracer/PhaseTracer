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

// ====================================================================

// This file was created by Claude Opus 5.0, provided with the equivalent
// file from TransitionListener2 (TL2). We claim no authorship of the code 
// in this file, and all credit should be given to the authors of TL2. 
// arXiv key: 
// ====================================================================


#ifndef POTENTIAL_THDM_HPP_INCLUDED
#define POTENTIAL_THDM_HPP_INCLUDED

#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "logger.hpp"
#include "one_loop_potential.hpp"
#include "pow.hpp"

#include "thdm_counterterms.hpp"
using namespace std;


namespace EffectivePotential {

class THDM : public OneLoopPotential {
 public:
  
  /**
   * @param m12_sq_      m_{12}^2 in GeV^2 (note: the squared mass, matching
   *                     TransitionListener's `m12_sq_GeV2`).
   * @param yukawa_type_ 1 = Type I, 2 = Type II, 3 = lepton specific, 4 = flipped.
   * @param RG_scale_    MS-bar renormalisation scale; <= 0 means "use v".
   */
  void init_params(double tanb_, double m12_sq_, double lam1_, double lam2_, double lam3_, double lam4_, double lam5_,
                   int yukawa_type_ = 2, double RG_scale_ = -1.) {
  	tanb = tanb_;
  	m12_sq = m12_sq_;
  	lam1 = lam1_;
  	lam2 = lam2_;
  	lam3 = lam3_;
  	lam4 = lam4_;
  	lam5 = lam5_;
  	lam345 = lam3 + lam4 + lam5;
  	cosb_sq = 1. / (1 + pow(tanb, 2));
  	sinb_sq = 1. - cosb_sq;
  	sin2b = 2. / (tanb + 1./tanb);

  	if (yukawa_type_ < 1 || yukawa_type_ > 4) {
  	  throw std::runtime_error("Invalid Yukawa type: expected an integer in [1, 4]");
  	}
  	yukawa_type = yukawa_type_;
  	// Type I and III (lepton specific) couple the down quarks to Phi_2;
  	// Type I and IV (flipped) couple the charged leptons to Phi_2.
  	bottom_uses_phi2 = (yukawa_type == 1 || yukawa_type == 3);
  	lepton_uses_phi2 = (yukawa_type == 1 || yukawa_type == 4);

  	// The vevs carry the sign of tan(beta), as in TransitionListener.
  	cosb = sqrt(cosb_sq);
  	sinb = (tanb < 0.) ? -sqrt(sinb_sq) : sqrt(sinb_sq);
  	v1 = vh * cosb;
  	v2 = vh * sinb;

  	// Yukawa couplings are defined by m_f = y_f <phi_i>, i.e. without a
  	// factor of sqrt(2), so that m_f^2(phi) = (y_f phi_i)^2.
  	Yt = mt / v2;
  	Yb = mb / (bottom_uses_phi2 ? v2 : v1);
  	Ytau = mtau / (lepton_uses_phi2 ? v2 : v1);

  	// Scalar Debye coefficients: Pi_i = CTempC_i T^2 (BSMPT conventions).
  	// Note the tau Yukawa is deliberately absent, matching TransitionListener.
  	const double coeff_gauge = 3. * (3. * pow(Y1, 2) + pow(Y2, 2));
  	const double c_t = sqrt(2.) * Yt;
  	const double c_b = sqrt(2.) * Yb;
  	CTempC1 = (12. * lam1 + 8. * lam3 + 4. * lam4 + coeff_gauge) / 48.;
  	CTempC2 = (12. * lam2 + 8. * lam3 + 4. * lam4 + coeff_gauge + 12. * c_t * c_t) / 48.;
  	if (bottom_uses_phi2) {
  	  CTempC2 += 12. * c_b * c_b / 48.;
  	} else {
  	  CTempC1 += 12. * c_b * c_b / 48.;
  	}

  	// Tree-level mass parameters, fixed by the EWSB (tadpole) conditions.
  	// V0 below writes them out inline; these are kept for the counterterm
  	// solver and for diagnostics.
  	m11_sq = m12_sq * tanb - 0.5 * vhsq * sinb_sq * lam345 - 0.5 * vhsq * cosb_sq * lam1;
  	m22_sq = m12_sq / tanb - 0.5 * vhsq * cosb_sq * lam345 - 0.5 * vhsq * sinb_sq * lam2;

  	RG_scale = (RG_scale_ > 0.) ? RG_scale_ : vh;
  	// The base class keeps its own renormalisation scale, used by V1; keep the
  	// two in step rather than relying on its default happening to match.
  	set_renormalization_scale(RG_scale);

  	THDM_tensors::Parameters tensor_params;
  	tensor_params.lambda1 = lam1;
  	tensor_params.lambda2 = lam2;
  	tensor_params.lambda3 = lam3;
  	tensor_params.lambda4 = lam4;
  	tensor_params.lambda5 = lam5;
  	tensor_params.m11sq = m11_sq;
  	tensor_params.m22sq = m22_sq;
  	tensor_params.m12sq = m12_sq;
  	tensor_params.v1 = v1;
  	tensor_params.v2 = v2;
  	tensor_params.Cg = Y1;
  	tensor_params.Cgs = Y2;
  	ct = THDM_ct::compute(tensor_params, yukawa_type, RG_scale);

  	if (ct.residual > 1e-6) {
  	  LOG(warning) << "Renormalised Coleman-Weinberg conditions deviate by "
  	               << ct.residual << " in the 2HDM counterterm solve";
  	}
  }
  
  size_t get_n_scalars() const override {return 2;}

  /**
   * The counterterms in TransitionListener's normalisation, i.e. as they appear
   * in V_ct = 1/2 dm11 h1^2 + 1/2 dm22 h2^2 - dm12 h1 h2
   *        + 1/8 dl1 h1^4 + 1/8 dl2 h2^4 + 1/4 (dl3 + dl4 + dl5) h1^2 h2^2.
   */
  std::map<std::string, double> get_counterterms() const {
    return {
      {"delta_m11_sq", ct.dm11_sq},
      {"delta_m22_sq", ct.dm22_sq},
      {"delta_m12_sq", ct.dm12_sq},
      {"delta_lambda1", ct.dlambda1},
      {"delta_lambda2", ct.dlambda2},
      {"delta_lambda3", ct.dlambda3},
      {"delta_lambda4", ct.dlambda4},
      {"delta_lambda5", ct.dlambda5},
      {"residual", ct.residual},
    };
  }

  /** Derived quantities, for cross-checking against TransitionListener. */
  std::map<std::string, double> get_diagnostics() const {
    return {
      {"g2", Y1}, {"g1", Y2},
      {"yt", Yt}, {"yb", Yb}, {"ytau", Ytau},
      {"v1", v1}, {"v2", v2},
      {"m11_sq", m11_sq}, {"m22_sq", m22_sq}, {"m12_sq", m12_sq},
      {"CTempC1", CTempC1}, {"CTempC2", CTempC2},
      {"renormScaleSq", RG_scale * RG_scale},
    };
  }
  
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
    
    double r = 1./2 * m12_sq * tanb * pow(h1 - h2 / tanb, 2) - vhsq / 4 * (lam1 * pow(h1, 2) + lam2 * pow(h2 * tanb, 2)) / (1 + pow(tanb, 2))
    - vhsq / 4 * lam345 * (pow(h1 * tanb, 2) +  pow(h2, 2)) / (1 + pow(tanb,2)) + 1./8 * pow(h1, 4) * lam1 + 1./8 * pow(h2, 4) * lam2 + 1./4 * pow(h1, 2) * pow(h2, 2) * lam345;
    return r;
  }
//  vct--------
  double counter_term(Eigen::VectorXd phi, double T) const override {
    const double h1 = phi[0];
    const double h2 = phi[1];
    return 0.5 * ct.dm11_sq * pow(h1, 2) + 0.5 * ct.dm22_sq * pow(h2, 2)
           - ct.dm12_sq * h1 * h2
           + 1. / 8 * ct.dlambda1 * pow(h1, 4) + 1. / 8 * ct.dlambda2 * pow(h2, 4)
           + 1. / 4 * (ct.dlambda3 + ct.dlambda4 + ct.dlambda5) * pow(h1, 2) * pow(h2, 2);
  }
  
  vector<double> get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const override{
    const double h1 = phi[0];
    const double h2 = phi[1];
    
    // Debye masses, added to the diagonal of each 2x2 block before diagonalising.
    const double T11 = CTempC1 * pow(T, 2);
    const double T22 = CTempC2 * pow(T, 2);
    // mass matrix of Higgs
    const double M11 = 3. / 2 * lam1 * pow(h1, 2) + 1. / 2 * lam345 * pow(h2, 2) + m12_sq * tanb - 1. / 2 * lam1 * vhsq * cosb_sq -
             1. / 2 * lam345 * vhsq * sinb_sq ;
    const double M22 = 3. / 2 * lam2 * pow(h2, 2) + 1. / 2 * lam345 * pow(h1, 2) + m12_sq / tanb - 1. / 2 * lam2 * vhsq * sinb_sq -
             1. / 2 * lam345 * vhsq * cosb_sq ;
    const double M12 = lam345 * h1 * h2 - m12_sq;
    const double M11_T = M11 + T11;
    const double M22_T = M22 + T22;

    const double mhT_sq = 1. / 2 * (M11_T + M22_T - sqrt(pow(M11_T, 2) + 4 * pow(M12, 2) + pow(M22_T, 2) - 2 * M11_T * M22_T));
    const double mHT_sq = 1. / 2 * (M11_T + M22_T + sqrt(pow(M11_T, 2) + 4 * pow(M12, 2) + pow(M22_T, 2) - 2 * M11_T * M22_T));
    // mass matrix of A
    const double MA11 = 1. / 2 * lam1 * pow(h1, 2) + m12_sq * tanb - 1. / 2 * lam1 * vhsq * cosb_sq - 1. / 2 * lam345 * vhsq * sinb_sq + 1. / 2 * (lam3 + lam4 - lam5) * 				 pow(h2, 2);  
    const double MA22 = 1. / 2 * lam2 * pow(h2, 2) + m12_sq / tanb - 1. / 2 * lam2 * vhsq * sinb_sq - 1. / 2 * lam345 * vhsq * cosb_sq + 1. / 2 * (lam3 + lam4 - lam5) * 				 pow(h1, 2); 
    const double MA12 = lam5 * h1 * h2 - m12_sq;
    const double MA11_T = MA11 + T11;
    const double MA22_T = MA22 + T22;
    
    const double mAT_sq = 1. / 2  * (MA11_T + MA22_T + sqrt(pow(MA11_T, 2) + 4 * pow(MA12, 2) + pow(MA22_T, 2) - 2 * MA11_T * MA22_T));
    const double mG0T_sq = 1. / 2  * (MA11_T + MA22_T - sqrt(pow(MA11_T, 2) + 4 * pow(MA12, 2) + pow(MA22_T, 2) - 2 * MA11_T * MA22_T));
    // mass matrix of Hpm
    const double MC11 = 1. / 2 * lam1 * pow(h1, 2) + m12_sq * tanb - 1. / 2 * lam1 * vhsq * cosb_sq - 1. / 2 * lam345 * vhsq *  sinb_sq + 1. / 2 * lam3 * pow(h2, 2);
    const double MC22 = 1. / 2 * lam2 * pow(h2, 2) + m12_sq / tanb - 1. / 2 * lam2 * vhsq * sinb_sq - 1. / 2 * lam345 * vhsq *  cosb_sq + 1. / 2 * lam3 * pow(h1, 2);
    const double MC12 = 1. / 2 * (lam4 + lam5) * h2 * h1 - m12_sq;
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

//ni of w, z, photon (transverse: 4, 2, 2; longitudinal: 2, 1, 1 -> 12 in total)
 vector<double> get_vector_dofs() const override { return {4., 2., 2., 2., 1., 1.};}


// top, bottom and tau; which doublet the down-type fermions see is set by the
// Yukawa type.
 vector<double> get_fermion_masses_sq(Eigen::VectorXd phi) const override {
    const double phi_bottom = bottom_uses_phi2 ? phi[1] : phi[0];
    const double phi_lepton = lepton_uses_phi2 ? phi[1] : phi[0];
    return {square(phi[1] * Yt), square(phi_bottom * Yb), square(phi_lepton * Ytau)};
  }
// ni of top, bottom and tau
 vector<double> get_fermion_dofs() const override {
    return {12., 12., 4.};
  }

 vector<double> get_tree_minimum() const {
    return {v1, v2};
  }
 
 private:
  
  // SM inputs, matching TransitionListener's constants.py (BSMPT defaults).
  const double vh = 246.21965079413735;
  const double vhsq = vh * vh;
  const double mW = 80.385;
  const double mZ = 91.1876;
  const double mt = 172.5;
  const double mb = 4.92;
  const double mtau = 1.77682;
  /** SU(2) gauge coupling g. */
  const double Y1 = 2.0 * mW / vh;
  /** U(1)_Y gauge coupling g'. */
  const double Y2 = sqrt(4.0 * mZ * mZ / vhsq - Y1 * Y1);

  // Yukawas and Debye coefficients depend on tan(beta), so they are set in
  // init_params rather than here.
  double Yt;
  double Yb;
  double Ytau;
  double CTempC1;
  double CTempC2;

  double tanb;
  double cosb;
  double sinb;
  double cosb_sq;
  double sinb_sq;
  double sin2b;
  double v1;
  double v2;
  int yukawa_type;
  bool bottom_uses_phi2;
  bool lepton_uses_phi2;
  double RG_scale;
  double m11_sq;
  double m22_sq;
  double m12_sq;
  double lam1;
  double lam2;
  double lam3;
  double lam4;
  double lam5;
  double lam345;
  THDM_ct::Result ct;


};

}  // namespace EffectivePotential

#endif
