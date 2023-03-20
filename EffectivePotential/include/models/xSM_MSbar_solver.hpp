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


/**
 Determine parameters \mu_h^2, \mu_s^2, \lambda_h and \lambda_s
 for real scalar singlet extension of the Standard Model
 by solving the one-loop corrected tadpole conditions and Higgs masses
 in the \bar{MS} scheme.
 
 To avoid massless Goldstone in V_CW, the vacuum is defined by the
 1-loop tadpole condition, which is same to arXiv:1808.01089.
 
 Notation:
 V_tree(h,s) = \mu_h^2 h^2 + \mu_s^2 s^2 + \lambda_h h^4 + \lambda_s s^4
          + c_1 h^2 s + c_2 h^2 s^2 + b_3 s^3
 
*/

#ifndef XSM_MSBAR_PARAMETERS_SOLVER_HPP_INCLUDED
#define XSM_MSBAR_PARAMETERS_SOLVER_HPP_INCLUDED

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Core>
#include "pow.hpp"
#include "SM_parameters.hpp"

double logx(double x) {
  const double abs_x = std::abs(x);
  if (abs_x <= std::numeric_limits<double>::min()) {
    return 0.;
  } else {
    return std::log(abs_x);
  }
}

class xSM_MSbar_parameters_solver {
 public:
  xSM_MSbar_parameters_solver(double ms_,
                              double c1_,
                              double c2_,
                              double b3_,
                              double vs_,
                              double lambda_s_):
  ms(ms_), c1(c1_), c2(c2_), b3(b3_), vs(vs_), lambda_s(lambda_s_){
    vev << v, vs;
    Z2_symmetry = vs == 0;
  }
  bool solving_parameters(){
    auto const valid = DetermineParameters();
    return valid;
  }
  void set_renormalization_scale(double Q_){
    Q = Q_;
  }
  
  double get_mu_h_Sq(){return mu_h_Sq_iteration;}
  double get_lambda_h(){return lambda_h_iteration;}
  double get_mu_s_Sq(){return mu_s_Sq_iteration;}
  double get_lambda_s(){return lambda_s_iteration;}

  double get_mu_h_Sq_tree_EWSB(){return mu_h_Sq_tree_EWSB;}
  double get_mu_s_Sq_tree_EWSB(){return mu_s_Sq_tree_EWSB;}

  double get_mu_h_Sq_tree(){return mu_h_Sq_tree;}
  double get_lambda_h_tree(){return lambda_h_tree;}
  double get_mu_s_Sq_tree(){return mu_s_Sq_tree;}
  double get_lambda_s_tree(){return lambda_s_tree;}
  
 private:
  
  // SM parameters
  const double v = SM::v;
  const double mh = SM::mh;
  const double mtop = SM::mtop;
  const double mZ = SM::mZ;
  const double mW = SM::mW;
  const double g = SM::g;
  const double g_sq = g * g;
  const double gp = SM::gp;
  const double gp_sq = gp * gp;
  const double yt_sq = SM::yt_sq;
  
  // renormalization scale
  double Q=173.03;
  size_t print_level = 0;
  bool Z2_symmetry = false;
  
  // input parameters
  // TODO PA: we shoudl really chnage the names of these paremeters
  // to either something more intuitive or something widely used in literature 
  double ms;
  double c1;
  double c2;
  double b3;
  double vs;
  double lambda_s;

  // output parameters
  double lambda_h_iteration;
  double lambda_s_iteration;
  double mu_h_Sq_iteration;
  double mu_s_Sq_iteration;

  // mu^2 satisfy tree-level EWSB condition with lambda_iteration
  double mu_h_Sq_tree_EWSB;
  double mu_s_Sq_tree_EWSB;
  
  double lambda_h_tree;
  double lambda_s_tree;
  double mu_h_Sq_tree;
  double mu_s_Sq_tree;

  double tol=1e-5;
  Eigen::VectorXd vev = Eigen::VectorXd::Zero(2);
  
  const double zero = 0;
  
  bool DetermineParameters(){
    bool valid = DetermineParametersAtTreeLevel();
    
    if (not valid) {
      if (print_level == 1)
        std::cout << "Invalid at tree-level" << std::endl;
      return false;
    }
    
    const auto mhsq_tree = massEigenvalues(vev);
    if (print_level > 1 ) {
      std::cout << "mu_h_Sq_tree:" << mu_h_Sq_tree << std::endl;
      std::cout << "mu_s_Sq_tree:" << mu_s_Sq_tree << std::endl;
      std::cout << "lambda_h_tree:" << lambda_h_tree << std::endl;
      std::cout << "lambda_s_tree:" << lambda_s_tree << std::endl;
      std::cout << "mh_1_tree:" << sqrt(std::abs(mhsq_tree[0])) << std::endl;
      std::cout << "mh_2_tree:" << sqrt(std::abs(mhsq_tree[1])) << std::endl;
    }
    
    // Initialize the Lagrangian parameters
    lambda_h_iteration = lambda_h_tree;
    lambda_s_iteration = lambda_s_tree;
    mu_h_Sq_iteration = mu_h_Sq_tree;
    mu_s_Sq_iteration = mu_s_Sq_tree;
    
    double lambda_h_prev = lambda_h_iteration;
    double lambda_s_prev = lambda_s_iteration;
    double mu_h_Sq_prev = mu_h_Sq_iteration;
    double mu_s_Sq_prev = mu_s_Sq_iteration;
    
    bool hasConverged = false;
    int i = 0;
    const int iterMax = 50;
    
    while (not hasConverged){
      if (i > iterMax){
        if (print_level == 1)
          std::cout << "Failed to converge after" << iterMax << "iterations." << std::endl;
        return false;
      }
      
      if (print_level >2){
        std::cout << i << ", "<< mu_h_Sq_iteration << ", "<< mu_s_Sq_iteration << ", " << lambda_h_iteration << ", " << lambda_s_iteration << ", " << std::endl;
        const auto mhsq = massEigenvalues_full(vev);
        std::cout << "mh_1:" << sqrt(std::abs(mhsq[0])) << std::endl;
        std::cout << "mh_2:" << sqrt(std::abs(mhsq[1])) << std::endl;
        
      }
      i++;
      
      lambda_h_prev = lambda_h_iteration;
      lambda_s_prev = lambda_s_iteration;
      mu_h_Sq_prev = mu_h_Sq_iteration;
      mu_s_Sq_prev = mu_s_Sq_iteration;
      
      valid = DetermineParametersAtOneLoop();
      
      if (not valid) {
        if (print_level >= 1)
          std::cout << "Invalid at one-loop after iteration:" << i << std::endl;
        return false;
      }
      
      const double dmuh = std::abs(mu_h_Sq_iteration - mu_h_Sq_prev);
      const double dmus = std::abs(mu_s_Sq_iteration - mu_s_Sq_prev);
      const double dlh = std::abs(lambda_h_iteration - lambda_h_prev);
      const double dls = std::abs(lambda_s_iteration - lambda_s_prev);
      if (print_level > 2 ){
        std::cout << "d_muh:" << dmuh << std::endl;
        std::cout << "d_mus:" << dmus << std::endl;
        std::cout << "d_lh:" << dlh << std::endl;
        std::cout << "d_ls:" << dls << std::endl;
      }
      
      hasConverged = (dmuh < tol) and (dlh < tol) and (dmus < tol) and (dls < tol);
    }
    if (print_level >= 1 ){
      std::cout << "Converged at"<< std::endl;
      std::cout << "  muh:" << mu_h_Sq_iteration << std::endl;
      std::cout << "  mus:" << mu_s_Sq_iteration << std::endl;
      std::cout << "  lh:" << lambda_h_iteration << std::endl;
      std::cout << "  ls:" << lambda_s_iteration << std::endl;
      const auto mhsq = massEigenvalues_full(vev);
      std::cout << "  mh_1:" << sqrt(std::abs(mhsq[0])) << std::endl;
      std::cout << "  mh_2:" << sqrt(std::abs(mhsq[1])) << std::endl;
      
      std::cout << "=================== "<< std::endl;
      Eigen::VectorXd x(2);
      x << 246.221, 0;
      std::cout << "d2Vdh2:" << sqrt(d2Vdh2(x)) << std::endl;
      std::cout << "d2Vds2:" << sqrt(d2Vds2(x)) << std::endl;
      std::cout << "d2V0dh2:" << sqrt(std::abs(d2V0dh2(x))) << std::endl;
      std::cout << "d2V0ds2:" << sqrt(std::abs(d2V0ds2(x))) << std::endl;
      std::cout << "d2VCWdh2:" << sqrt(std::abs(d2VCWdh2(x))) << std::endl;
      std::cout << "d2VCWds2:" << sqrt(std::abs(d2VCWds2(x))) << std::endl;
      std::cout << "d2VCWds2_test:" << sqrt(std::abs(d2VCWds2_test(x))) << std::endl;
      
      std::cout << "V1:" << VCW(x) << std::endl;
      
      std::cout << "Numerically derivatives of the CW potential at EW VEV:" << std::endl;
       // derivate
       const std::vector<double> n_h_xx = {-1., 0., 1.};
       const double h = 0.001;
       const std::vector<double> coeff_xx = {1., -2., 1.};
            
       Eigen::MatrixXd hessianCW = Eigen::MatrixXd::Zero(x.size(), x.size());
       for (int ii = 0; ii < x.size(); ++ii) {
         Eigen::VectorXd phi_shifted = x;
         for (int jj = 0; jj < n_h_xx.size(); ++jj) {
           phi_shifted(ii) = x(ii) + n_h_xx[jj] * h;
//           std::cout << "VCW [" << ii << ", " << jj <<"] =" << std::setprecision(16) << VCW(phi_shifted) << std::endl;
           hessianCW(ii, ii) += VCW(phi_shifted) * coeff_xx[jj] / square(h);
         }
       }
       
       std::cout << "Sqrt[d^2VCW/dh^2] = "<< std::sqrt(abs(hessianCW(0,0))) << std::endl;
       std::cout << "Sqrt[d^2VCW/ds^2] = "<< std::sqrt(abs(hessianCW(1,1))) << std::endl;
      
    }
    return valid;
  }
  
  bool DetermineParametersAtTreeLevel(){
    const double mhSq = square(mh);
    const double msSq = square(ms);
    const double mhs = 2*c1*v + 4*c2*v*vs;
    
    const double a = mhSq + msSq;
    const double discriminant = square(mhSq - msSq) - 4*square(mhs);
    
    // TODO
    if (discriminant < 0)
      return false;
    const double b = discriminant < 0 ? 0 : sqrt(discriminant);
    
    const double mhh = mhSq > msSq ? 0.5*(a + b) : 0.5*(a - b);
    const double mss = mhSq > msSq ? 0.5*(a - b) : 0.5*(a + b);
    
    // Apply tree-level EWSB to constrain muh and mus, and use the Higgs and singlet masses to constrain lh and ls.
    lambda_h_tree  = mhh / (8*square(v));
    mu_h_Sq_tree = -2*lambda_h_tree*square(v) - c1*vs - c2*square(vs);
    if (not Z2_symmetry) {
      lambda_s_tree = (mss + c1*square(v)/vs - 3*b3*vs) / (8*square(vs));
      mu_s_Sq_tree = -2*lambda_s_tree*square(vs) - 0.5*c1*square(v)/vs
      - c2*square(v) - 1.5*b3*vs;
    } else {
      lambda_s_tree = lambda_s;
      mu_s_Sq_tree = 0.5*mss - c2*square(v);
    }
    mu_h_Sq_tree_EWSB = mu_h_Sq_tree;
    mu_s_Sq_tree_EWSB = mu_s_Sq_tree;
        
    return true;
  }
  
  bool DetermineParametersAtOneLoop(){
    const double dVCWdh_ = dVCWdh(vev);
    const double dVCWds_ = dVCWds(vev);
    const double d2VCWdh2_ = d2VCWdh2(vev);
    const double d2VCWds2_ = d2VCWds2(vev);
    const double d2VCWdhds_ = d2VCWdhds(vev);
    
    const double mhSq = square(mh);
    const double msSq = square(ms);
    
    const double a = mhSq + msSq;
    const double discriminant = square(mhSq - msSq) - 4*square(2*c1*v + 4*c2*v*vs                         + d2VCWdhds_);
    //TODO
    if (discriminant < 0)
      return false;
    const double b = discriminant < 0 ? 0 : sqrt(discriminant);
    
    const double mhh = mhSq > msSq ? 0.5*(a + b) : 0.5*(a - b);
    const double mss = mhSq > msSq ? 0.5*(a - b) : 0.5*(a + b);
    
    // Apply one-loop EWSB to constrain muh and mus, and use the Higgs and singlet masses to constrain lh and ls.
    lambda_h_iteration = (mhh + 1./v*dVCWdh_ - d2VCWdh2_) / (8*v*v);
    mu_h_Sq_iteration = -0.25*mhh - 0.75/v*dVCWdh_ + 0.25*d2VCWdh2_
                        -c1*vs -c2*vs*vs;
    if (not Z2_symmetry){
      lambda_s_iteration = (mss + c1*v*v/vs - 3*b3*vs + 1/vs*dVCWds_ - d2VCWds2_)
                          / (8.*vs*vs);
      mu_s_Sq_iteration = -0.25*mss - 0.75*(c1*v*v/vs + b3*vs + dVCWds_/vs - d2VCWds2_/3.)
                          - c2*v*v;
    } else {
      lambda_s_iteration = lambda_s ;
      mu_s_Sq_iteration = 0.5*mss - c2*square(v) - 0.5*d2VCWds2_;
//      mu_s_Sq_iteration = -0.25*mss - c2*square(v) + 0.25*d2VCWds2_;
    }
//    std::cout << "!! d2VCWds2_ = " << d2VCWds2_ <<std::endl;
//    std::cout << "!! mu_s_Sq_iteration = " << mu_s_Sq_iteration <<std::endl;
    
    // Store the quadratic parameters calculated from the quadratic parameters using tree-level EWSB. These are used
    // for masses that enter the Coleman-Weinberg potential.
    mu_h_Sq_tree_EWSB = - 2*lambda_h_iteration*square(v) - c1*vs -c2*vs*vs;
    if (not Z2_symmetry) {
      mu_s_Sq_tree_EWSB = -1./(2.*vs)*(4.*lambda_s_iteration*vs*vs*vs + c1*v*v + 2.*c2*v*v*vs
                                       + 3.*b3*vs*vs);
    } else {
//      mu_s_Sq_tree_EWSB = 0.5*mss - c2*square(v);
      mu_s_Sq_tree_EWSB = 0.5*mss - c2*square(v);
    }
    
//    std::cout << "mu_s_Sq_iteration" << mu_h_Sq_iteration << "\t";
//    std::cout << mu_s_Sq_iteration <<"\t";
//    std::cout << lambda_h_iteration <<"\t";
//    std::cout << lambda_s_iteration <<std::endl;
//    std::cout << "mu_h_Sq_tree_EWSB" << mu_h_Sq_tree_EWSB << "\t";
//    std::cout << mu_s_Sq_tree_EWSB <<"\t";
//    std::cout << std::endl;

    return true;
  }
  
  
  //The eigenvalues of the one-loop mass matrix, in descending order of mass. These are used to check if the masses are reproduced after constraining the parameters in the potential.
  std::vector<double> massEigenvalues_full(Eigen::VectorXd phi){
    const double mhh = d2Vdh2(phi);
    const double mhs = d2Vdhds(phi);
    const double mss = d2Vds2(phi);
    const double a = mhh + mss;
    const double b = sqrt(square(mhh - mss) + 4*square(mhs));
//    std::cout << "massEigenvalues_full = " << mhh << "\t"<< mhs << "\t"<< mss << std::endl;
    
    return {0.5*(a + b), 0.5*(a - b)};
  }
  
  // The eigenvalues of the tree-level mass matrix, in descending order of mass. These are used as masses in the Coleman-Weinberg potential.
  std::vector<double> massEigenvalues(Eigen::VectorXd phi){
    const double mhh = d2V0dh2(phi);
    const double mhs = d2V0dhds(phi);
    const double mss = d2V0ds2(phi);
    const double a = mhh + mss;
    const double b = sqrt(square(mhh - mss) + 4*square(mhs));
//    std::cout << "massEigenvalues Mss = " << std::setprecision(16) << mss << std::endl;
//    std::cout << "massEigenvalues = " << mhh << "\t"<< mhs << "\t"<< mss << std::endl;
    return {0.5*(a + b), 0.5*(a - b)};
  }
  
  std::vector<double> d_massEigenvalues_dh(Eigen::VectorXd phi){
    const double mhh = d2V0dh2(phi);
    const double mhs = d2V0dhds(phi);
    const double mss = d2V0ds2(phi);
  
    const double dmhhdh = d3V0dh3(phi);
    const double dmssdh = d3V0dhds2(phi);
    const double dmhsdh = d3V0dh2ds(phi);
  
    const double a = dmhhdh + dmssdh;
    const double b = ((dmhhdh - dmssdh)*(mhh - mss) + 4*dmhsdh*mhs) / sqrt(square(mhh - mss) + 4*square(mhs));
  
    return {0.5*(a + b), 0.5*(a - b)};
  }
  
  std::vector<double> d_massEigenvalues_ds(Eigen::VectorXd phi){
    const double mhh = d2V0dh2(phi);
    const double mhs = d2V0dhds(phi);
    const double mss = d2V0ds2(phi);
  
    const double dmhhds = d3V0dh2ds(phi);
    const double dmssds = d3V0ds3(phi);
    const double dmhsds = d3V0dhds2(phi);
  
    const double a = dmhhds + dmssds;
    const double b = ((dmhhds - dmssds) * (mhh - mss) + 4 * dmhsds * mhs) / sqrt(square(mhh - mss) + 4*square(mhs));
  
    return {0.5*(a + b), 0.5*(a - b)};
  }
  
  std::vector<double> d2_massEigenvalues_dh2(Eigen::VectorXd phi){
    const double mhh = d2V0dh2(phi);
    const double mhs = d2V0dhds(phi);
    const double mss = d2V0ds2(phi);

    const double dmhhdh = d3V0dh3(phi);
    const double dmssdh = d3V0dhds2(phi);
    const double dmhsdh = d3V0dh2ds(phi);
    
    const double d2mhhdh2 = d4V0dh4(phi);
    const double d2mssdh2 = d4V0dh2ds2(phi);
    const double d2mhsdh2 = d4V0dh3ds(phi);
  
    const double a = d2mhhdh2 + d2mssdh2;
    const double b = (((d2mhhdh2 - d2mssdh2)*(mhh - mss)
                    + square(dmhhdh - dmssdh) + 4*d2mhsdh2*mhs
                    + 4*square(dmhsdh))*(square(mhh - mss)
                    + 4*mhs) - square((dmhhdh - dmssdh)*(mhh - mss)
                    + 4*dmhsdh*mhs)) / pow((square(mhh - mss)
                    + 4*square(mhs)),1.5);
    return {0.5*(a + b), 0.5*(a - b)};
  }
  
  std::vector<double> d2_massEigenvalues_ds2(Eigen::VectorXd phi){
    const double mhh = d2V0dh2(phi);
    const double mhs = d2V0dhds(phi);
    const double mss = d2V0ds2(phi);
    
    const double dmhhds = d3V0dh2ds(phi);
    const double dmssds = d3V0ds3(phi);
    const double dmhsds = d3V0dhds2(phi);
    
    const double d2mhhds2 = d4V0dh2ds2(phi);
    const double d2mssds2 = d4V0ds4(phi);
    const double d2mhsds2 = d4V0dhds3(phi);
    
    const double a = d2mhhds2 + d2mssds2;
    const double b = (((d2mhhds2 - d2mssds2)*(mhh - mss)
                    + square(dmhhds - dmssds) + 4*d2mhsds2*mhs
                    + 4*square(dmhsds))*(square(mhh - mss)
                    + 4*mhs) - square((dmhhds - dmssds)*(mhh - mss)
                    + 4*dmhsds*mhs)) / pow((square(mhh - mss)
                    + 4*square(mhs)),1.5);
    return {0.5*(a + b), 0.5*(a - b)};
  }
  
  std::vector<double> d2_massEigenvalues_dhds(Eigen::VectorXd phi){
    const double mhh = d2V0dh2(phi);
    const double mhs = d2V0dhds(phi);
    const double mss = d2V0ds2(phi);
    
    const double dmhhdh = d3V0dh3(phi);
    const double dmssdh = d3V0dhds2(phi);
    const double dmhsdh = d3V0dh2ds(phi);
    
    const double dmhhds = d3V0dh2ds(phi);
    const double dmssds = d3V0ds3(phi);
    const double dmhsds = d3V0dhds2(phi);
  
    const double d2mhhdhds = d4V0dh3ds(phi);
    const double d2mssdhds = d4V0dhds3(phi);
    const double d2mhsdhds = d4V0dh2ds2(phi);
    
    const double a = d2mhhdhds + d2mssdhds;
    const double b = (((d2mhhdhds - d2mssdhds)*(mhh - mss)
                    + (dmhhdh - dmssdh)*(dmhhds - dmssds)
                    + 4*d2mhsdhds*mhs + 4*dmhsdh*dmhsds)*(square(mhh - mss)
                    + 4*mhs) - ((dmhhdh - dmssdh)*(mhh - mss)
                    + 4*dmhsdh*mhs) * ((dmhhds - dmssds)*(mhh - mss)
                    + 4*dmhsds*mhs)) / pow((square(mhh - mss)
                    + 4*square(mhs)),1.5);
//    std::cout << "d2_massEigenvalues_dhds=" << a << ", "<< b << std::endl;
    return {0.5*(a + b), 0.5*(a - b)};
  }
  
  /* ****************************************
  *******************************************
  Derivatives of the full one-loop potential up to second order.
  ******************************************
  ****************************************** */
  double dVdh(Eigen::VectorXd phi){
    const double dV0dh_ = dV0dh(phi);
    const double dVCWdh_ = dVCWdh(phi);
    return dV0dh_ + dVCWdh_;
  }
  double dVds(Eigen::VectorXd phi){
    const double dV0ds_ = dV0ds(phi);
    const double dVCWds_ = dVCWds(phi);
    return dV0ds_ + dVCWds_;
  }
  double d2Vdh2(Eigen::VectorXd phi){
    const double d2V0dh2_ = d2V0dh2(phi);
    const double d2VCWdh2_ = d2VCWdh2(phi);
//    std::cout << "d2Vdh2 = " << d2V0dh2_ << "\t"<< d2VCWdh2_ << std::endl;
    return d2V0dh2_ + d2VCWdh2_;
  }
  double d2Vdhds(Eigen::VectorXd phi){
    const double d2V0dhds_ = d2V0dhds(phi);
    const double d2VCWdhds_ = d2VCWdhds(phi);
    return d2V0dhds_ + d2VCWdhds_;
  }
  double d2Vds2(Eigen::VectorXd phi){
    const double d2V0ds2_ = d2V0ds2(phi);
    const double d2VCWds2_ = d2VCWds2(phi);
    return d2V0ds2_ + d2VCWds2_;
  }
  
  /* ****************************************
  *******************************************
  Derivatives of the tree-level potential up to fourth order using tree-level EWSB to calculate the quadratic
  parameters from the quartic parameters (i.e. muh from lh and mus from ls). The 'subscript' 'm' in V0m signifies
  these forms are for use in masses that enter the Coleman-Weinberg potential, where we want the parameters to
  be related through EWSB at tree-level only.
  ******************************************
  ****************************************** */
  
  double dV0mdh(Eigen::VectorXd phi){
    const auto h = phi[0];
    const auto s = phi[1];
    return 2*mu_h_Sq_tree_EWSB*h + 4*lambda_h_iteration*cube(h) + 2*c1*h*s + 2*c2*h*square(s);
  }
  double dV0mds(Eigen::VectorXd phi){
    const auto h = phi[0];
    const auto s = phi[1];
    return 2*mu_s_Sq_tree_EWSB*s + 4*lambda_s_iteration*cube(s) + c1*square(h) + 2*c2*square(h)*s + 3*b3*square(s);
  }
  // Second derivative of tree level potential
  double d2V0mdh2(Eigen::VectorXd phi){
    const auto h = phi[0];
    const auto s = phi[1];
    return 2*mu_h_Sq_tree_EWSB + 12*lambda_h_iteration*square(h) + 2*c1*s + 2*c2*square(s);
  }
  double d2V0mds2(Eigen::VectorXd phi){
    const auto h = phi[0];
    const auto s = phi[1];
    return 2*mu_s_Sq_tree_EWSB + 12*lambda_s_iteration*square(s) + 2*c2*square(h) + 6*b3*s;
  }
  double d2V0mdhds(Eigen::VectorXd phi){
    const auto h = phi[0];
    const auto s = phi[1];
    return 2*c1*h + 4*c2*h*s;
  }
  double d3V0mdh3(Eigen::VectorXd phi){
    const auto h = phi[0];
    return 24*lambda_h_iteration*h;
  }
  double d3V0mdh2ds(Eigen::VectorXd phi){
    const auto s = phi[1];
    return 2*c1 + 4*c2*s;
  }
  double d3V0mdhds2(Eigen::VectorXd phi){
    const auto h = phi[0];
    return 4*c2*h;
  }
  double d3V0mds3(Eigen::VectorXd phi){
    const auto s = phi[1];
    return 24*lambda_s_iteration*s + 6*b3;
  }
  double d4V0mdh4(Eigen::VectorXd phi){ return 24*lambda_h_iteration;}
  double d4V0mdh3ds(Eigen::VectorXd phi){ return 0.;}
  double d4V0mdh2ds2(Eigen::VectorXd phi){ return 4.*c2;}
  double d4V0mdhds3(Eigen::VectorXd phi){ return 0.;}
  double d4V0mds4(Eigen::VectorXd phi){ return 24*lambda_s_iteration;}
  
  /* ****************************************
  *******************************************
  Derivatives of the tree-level potential up to fourth order.
  ******************************************
  ****************************************** */
  double dV0dh(Eigen::VectorXd phi){
    const auto h = phi[0];
    const auto s = phi[1];
    return 2*mu_h_Sq_iteration*h + 4*lambda_h_iteration*cube(h) + 2*c1*h*s + 2*c2*h*square(s);
  }
  double dV0ds(Eigen::VectorXd phi){
    const auto h = phi[0];
    const auto s = phi[1];
    return 2*mu_s_Sq_iteration*s + 4*lambda_s_iteration*cube(s) + c1*square(h) + 2*c2*square(h)*s + 3*b3*square(s);
  }
  // Second derivative of tree level potential
  double d2V0dh2(Eigen::VectorXd phi){
    const auto h = phi[0];
    const auto s = phi[1];
    return 2*mu_h_Sq_iteration + 12*lambda_h_iteration*square(h) + 2*c1*s + 2*c2*square(s);
  }
  double d2V0ds2(Eigen::VectorXd phi){
    const auto h = phi[0];
    const auto s = phi[1];
//    std::cout << "massEigenvalues lambda_hs = " << 4*c2 << std::endl;
//    std::cout << "massEigenvalues mus_sq_tree_EWSB = " << 2*mu_s_Sq_iteration << std::endl;
    return 2*mu_s_Sq_iteration + 12*lambda_s_iteration*square(s) + 2*c2*square(h) + 6*b3*s;
  }
  double d2V0dhds(Eigen::VectorXd phi){
    const auto h = phi[0];
    const auto s = phi[1];
    return 2*c1*h + 4*c2*h*s;
  }
  double d3V0dh3(Eigen::VectorXd phi){
    const auto h = phi[0];
    return 24*lambda_h_iteration*h;
  }
  double d3V0dh2ds(Eigen::VectorXd phi){
    const auto s = phi[1];
    return 2*c1 + 4*c2*s;
  }
  double d3V0dhds2(Eigen::VectorXd phi){
    const auto h = phi[0];
    return 4*c2*h;
  }
  double d3V0ds3(Eigen::VectorXd phi){
    const auto s = phi[1];
    return 24*lambda_s_iteration*s + 6*b3;
  }
  double d4V0dh4(Eigen::VectorXd phi){ return 24*lambda_h_iteration;}
  double d4V0dh3ds(Eigen::VectorXd phi){ return 0.;}
  double d4V0dh2ds2(Eigen::VectorXd phi){ return 4.*c2;}
  double d4V0dhds3(Eigen::VectorXd phi){ return 0.;}
  double d4V0ds4(Eigen::VectorXd phi){ return 24*lambda_s_iteration;}
  
  /* ****************************************
  *******************************************
  Derivatives of the Coleman-Weinberg potential up to second order.
  ******************************************
  ****************************************** */
  double dVCWdh(Eigen::VectorXd phi){
    const auto mb = boson_massSq(phi,0);
    const auto dmbdh = d_boson_massSq_dh(phi);
    
    double y=0;
    for (size_t i = 0; i < mb.size(); ++i) {
      y += boson_dofs[i] * mb[i] * dmbdh[i]
      * (logx(std::abs(mb[i]/square(Q))+zero) - boson_cs[i] + 0.5);
    }
    
    const auto mf = fermion_massSq(phi);
    const auto dmfdh = d_fermion_massSq_dh(phi);
    for (size_t i = 0; i < mf.size(); ++i) {
      y -= fermion_dofs[i] * mf[i] * dmfdh[i]
      * (logx(std::abs(mf[i]/square(Q))+zero) - 1.5 + 0.5);
    }
    
    return y/(32*M_PI*M_PI);
  }
  double dVCWds(Eigen::VectorXd phi){
    const auto mb = boson_massSq(phi,0);
    const auto dmbds = d_boson_massSq_ds(phi);
    
    double y=0;
    for (size_t i = 0; i < mb.size(); ++i) {
      y += boson_dofs[i] * mb[i] * dmbds[i]
      * (logx(std::abs(mb[i]/square(Q))+zero) - boson_cs[i] + 0.5);
    }
    
    const auto mf = fermion_massSq(phi);
    const auto dmfds = d_fermion_massSq_ds(phi);
    for (size_t i = 0; i < mf.size(); ++i) {
      y -= fermion_dofs[i] * mf[i] * dmfds[i]
      * (logx(std::abs(mf[i]/square(Q))+zero) - 1.5 + 0.5);
    }
    
    return y/(32*M_PI*M_PI);
  }
  
  
  double VCW(Eigen::VectorXd phi){
    const auto mb = boson_massSq(phi,0);
    double y=0;
    for (size_t i = 0; i < mb.size(); ++i) {
//      if (i<=3) continue;
      y += boson_dofs[i] * mb[i] * mb[i]
           * (logx(mb[i]/square(Q)) - boson_cs[i]);
    }
    
    const auto mf = fermion_massSq(phi);
    for (size_t i = 0; i < mf.size(); ++i) {
      y -= fermion_dofs[i] * mf[i] * mf[i]
      * (logx(mf[i]/square(Q)) - 1.5);
    }
    return y/(64*M_PI*M_PI);
    
  }

  double d2VCWds2_test(Eigen::VectorXd phi){
    const auto mb = boson_massSq(phi,0);
    const auto dmbds = d_boson_massSq_ds(phi);
    const auto d2mbds2 = d2_boson_massSq_ds2(phi);
    double y=0;
//        std::cout << std::endl
//        <<"mb=" << mb[0]
//        <<"\t" << mb[1]
//        <<"\t" << mb[2]
//        <<"\t" << mb[3]
//        <<"\t" << mb[4]
//        <<"\t" << mb[5]
//        << std::endl;
    for (size_t i = 0; i < mb.size(); ++i) {
      if (i<=4) continue;
      y += boson_dofs[i] * ( (square(dmbds[i]) + mb[i] * d2mbds2[i])
           * (logx(std::abs(mb[i]/square(Q))+zero) - boson_cs[i]+0.5)
           + square(dmbds[i]));
    }
    
//    const auto mf = fermion_massSq(phi);
//    const auto dmfds = d_fermion_massSq_ds(phi);
//    const auto d2mfds2 = d2_fermion_massSq_ds2(phi);
//    for (size_t i = 0; i < mf.size(); ++i) {
//      y -= fermion_dofs[i] * ( (square(dmfds[i]) + mf[i] * d2mfds2[i])
//           * (logx(std::abs(mf[i]/square(Q))+zero) - 1.5+0.5)
//           + square(dmfds[i]));
//    }
    return y/(32*M_PI*M_PI);
  }
  
  double d2VCWdh2(Eigen::VectorXd phi){
    const auto mb = boson_massSq(phi,0);
//    std::cout << "!!!!!! d2VCWdh2 mb=" << mb[0]
//    <<"\t" << mb[1]
//    <<"\t" << mb[2]
//    <<"\t" << mb[3]
//    <<"\t" << mb[4]
//    <<"\t" << mb[5]
//    << std::endl;
    const auto dmbdh = d_boson_massSq_dh(phi);
//    std::cout << "!!!!!! d2VCWdh2 dmbdh=" << dmbdh[0]
//    <<"\t" << dmbdh[1]
//    <<"\t" << dmbdh[2]
//    <<"\t" << dmbdh[3]
//    <<"\t" << dmbdh[4]
//    <<"\t" << dmbdh[5]
//    << std::endl;
    const auto d2mbdh2 = d2_boson_massSq_dh2(phi);
//    std::cout << "!!!!!! d2mbdh2 d2mbdh2=" << d2mbdh2[0]
//    <<"\t" << d2mbdh2[1]
//    <<"\t" << d2mbdh2[2]
//    <<"\t" << d2mbdh2[3]
//    <<"\t" << d2mbdh2[4]
//    <<"\t" << d2mbdh2[5]
//    << std::endl;
    double y=0;
    for (size_t i = 0; i < mb.size(); ++i) {
      y += boson_dofs[i] * ( (square(dmbdh[i]) + mb[i] * d2mbdh2[i])
           * (logx(std::abs(mb[i]/square(Q))+zero) - boson_cs[i]+0.5)
           + square(dmbdh[i]));
//      std::cout << "!!!!!-----! y=" <<
//      logx(std::abs(mb[i]/square(Q))+zero)
//      <<"\t" << square(Q)
//       << std::endl;
    }
    
    const auto mf = fermion_massSq(phi);
    const auto dmfdh = d_fermion_massSq_dh(phi);
    const auto d2mfdh2 = d2_fermion_massSq_dh2(phi);
    for (size_t i = 0; i < mf.size(); ++i) {
      y -= fermion_dofs[i] * ( (square(dmfdh[i]) + mf[i] * d2mfdh2[i])
           * (logx(std::abs(mf[i]/square(Q))+zero) - 1.5+0.5)
           + square(dmfdh[i]));
//      std::cout << "!!!!!! y=" << y << std::endl;
    }
    return y/(32*M_PI*M_PI);
  }
  double d2VCWds2(Eigen::VectorXd phi){
    const auto mb = boson_massSq(phi,0);
    const auto dmbds = d_boson_massSq_ds(phi);
    const auto d2mbds2 = d2_boson_massSq_ds2(phi);
    double y=0;
    
    for (size_t i = 0; i < mb.size(); ++i) {
      y += boson_dofs[i] * ( (square(dmbds[i]) + mb[i] * d2mbds2[i])
           * (logx(std::abs(mb[i]/square(Q))+zero) - boson_cs[i]+0.5)
           + square(dmbds[i]));
    }
    
    const auto mf = fermion_massSq(phi);
    const auto dmfds = d_fermion_massSq_ds(phi);
    const auto d2mfds2 = d2_fermion_massSq_ds2(phi);
    for (size_t i = 0; i < mf.size(); ++i) {
      y -= fermion_dofs[i] * ( (square(dmfds[i]) + mf[i] * d2mfds2[i])
           * (logx(std::abs(mf[i]/square(Q))+zero) - 1.5+0.5)
           + square(dmfds[i]));
    }
    return y/(32*M_PI*M_PI);
  }
  double d2VCWdhds(Eigen::VectorXd phi){
    const auto mb = boson_massSq(phi,0);
    const auto dmbds = d_boson_massSq_ds(phi);
    const auto dmbdh = d_boson_massSq_dh(phi);
    const auto d2mbdhds = d2_boson_massSq_dhds(phi);
    double y=0;
    for (size_t i = 0; i < mb.size(); ++i) {
      y += boson_dofs[i] * ( (dmbdh[i]*dmbds[i] + mb[i] * d2mbdhds[i])
           * (logx(std::abs(mb[i]/square(Q))+zero) - boson_cs[i]+0.5)
           + dmbdh[i]*dmbds[i]);
    }
    
    const auto mf = fermion_massSq(phi);
    const auto dmfdh = d_fermion_massSq_dh(phi);
    const auto dmfds = d_fermion_massSq_ds(phi);
    const auto d2mfdhds = d2_fermion_massSq_dhds(phi);
    for (size_t i = 0; i < mf.size(); ++i) {
      y -= fermion_dofs[i] * ( (dmfds[i]*dmfdh[i] + mf[i] * d2mfdhds[i])
           * (logx(std::abs(mf[i]/square(Q))+zero) - 1.5+0.5)
           + dmfds[i]*dmfdh[i]);
    }
    return y/(32*M_PI*M_PI);
  }
  
  
  /* ****************************************
  *******************************************
  Boson masses and their derivatives.
  ******************************************
  ****************************************** */
  std::vector<double> boson_dofs = {3, 6 ,3, 2, 1, 1};
  std::vector<double> boson_cs = {1.5, 5./6. ,5./6., 5./6., 1.5, 1.5};
  std::vector<double> boson_massSq(Eigen::VectorXd phi, double T){
    const double h = phi[0];
    const double s = phi[1];
    const double h_sq = square(h);
    const double s2 = square(s);
  
    const double mW_Sq = 0.25 * g_sq * h_sq;
    const double mZ_Sq = 0.25 * (g_sq + gp_sq) * h_sq;
    const double mPh_Sq = 0;
    
    const auto m = massEigenvalues(phi);
    const double mp = m[0];
    const double mm = m[1];
    const double mgb_Sq = 2*mu_h_Sq_iteration + 4*lambda_h_iteration*h_sq + 2*c1*s + 2*c2*s2;
    
    return {mgb_Sq, mW_Sq, mZ_Sq, mPh_Sq, mp, mm};
  }
 
  std::vector<double> d_boson_massSq_dh(Eigen::VectorXd phi){
    const double h = phi[0];
    
    const double dmW = g_sq*h/2;
    const double dmZ = (g_sq + gp_sq)*h/2;
    
    const double dmgb = 8*lambda_h_iteration*h;
    const double dmPh = 0;
    
    const auto dm = d_massEigenvalues_dh(phi);
    const double dmp = dm[0];
    const double dmm = dm[1];
    return {dmgb, dmW, dmZ, dmPh, dmp, dmm};
  }

  std::vector<double> d_boson_massSq_ds(Eigen::VectorXd phi){
    const double s = phi[1];
    
    const double dmW = 0;
    const double dmZ = 0;
    const double dmgb = 2*c1 + 4*c2*s;
    const double dmPh = 0;
    
    const auto dm = d_massEigenvalues_ds(phi);
    const double dmp = dm[0];
    const double dmm = dm[1];
    return {dmgb, dmW, dmZ, dmPh, dmp, dmm};
  }
  std::vector<double> d2_boson_massSq_dh2(Eigen::VectorXd phi){
    const double d2mW = g_sq/2;
    const double d2mZ = (g_sq + gp_sq)/2;
    const double d2mgb = 8*lambda_h_iteration;
    const double d2mPh = 0;
    
    const auto d2m = d2_massEigenvalues_dh2(phi);
    const double d2mp = d2m[0];
    const double d2mm = d2m[1];
    return {d2mgb, d2mW, d2mZ, d2mPh, d2mp, d2mm};
  }
  std::vector<double> d2_boson_massSq_ds2(Eigen::VectorXd phi){
    const double d2mW = 0;
    const double d2mZ = 0;
    const double d2mgb = 4*c2;
    const double d2mPh = 0;
    
    const auto d2m = d2_massEigenvalues_ds2(phi);
    const double d2mp = d2m[0];
    const double d2mm = d2m[1];
    return {d2mgb, d2mW, d2mZ, d2mPh, d2mp, d2mm};
  }
  std::vector<double> d2_boson_massSq_dhds(Eigen::VectorXd phi){
    const double d2mW = 0;
    const double d2mZ = 0;
    const double d2mgb = 0;
    const double d2mPh = 0;
    
    const auto d2m = d2_massEigenvalues_dhds(phi);
    const double d2mp = d2m[0];
    const double d2mm = d2m[1];
    return {d2mgb, d2mW, d2mZ, d2mPh, d2mp, d2mm};
  }
  
  /* ****************************************
  *******************************************
  Fermion masses and their derivatives.
  *******************************************
  ******************************************* */
  std::vector<double> fermion_dofs = {12};
  std::vector<double> fermion_massSq(Eigen::VectorXd phi){
    return {0.5 * yt_sq * square(phi[0])};
  }
  
  std::vector<double> d_fermion_massSq_dh(Eigen::VectorXd phi){
    return {yt_sq * phi[0]};
  }
  std::vector<double> d_fermion_massSq_ds(Eigen::VectorXd phi){
    return {0};
  }
  std::vector<double> d2_fermion_massSq_dh2(Eigen::VectorXd phi){
    return {yt_sq};
  }
  std::vector<double> d2_fermion_massSq_ds2(Eigen::VectorXd phi){
    return {0};
  }
  std::vector<double> d2_fermion_massSq_dhds(Eigen::VectorXd phi){
    return {0};
  }
};

#endif
