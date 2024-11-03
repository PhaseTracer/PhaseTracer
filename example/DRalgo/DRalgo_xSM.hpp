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

#ifndef POTENTIAL_DRALGO_xSM_HPP_INCLUDED
#define POTENTIAL_DRALGO_xSM_MODEL_HPP_INCLUDED

/**
   The xSMin  DRalgo
   https://arxiv.org/pdf/22xx.xxxxx.pdf

  Sqrt -> sqrt, pow -> pow, 2 -> 2.
  M_PI -> M_PI, log-> log
   /.{\[CurlyPhi]^2 -> Hsq, \[CurlyPhi]^4 -> Hsq^2, s^2->Ssq, s^4->Ssq^2,lambdaH->lamH,lambdaHS->lamHS, lambdaS->lamS,\[Mu]Ssq->muSsq,\[Mu]Hsq->muHsq,g1^2 -> g1sq, g2^2 -> g2sq, g13d^2->g13dsq,g23d^2->g23dsq,g33d^2->g33dsq,lambdaH3d->lamH3d, lambdaHS3d->lamHS3d,lambdaS3d->lamS3d}
 

*/

#include <vector>
#include <cmath>
#include <complex>
#include <Eigen/Eigenvalues>
#include <interpolation.h>

#include "pow.hpp"
#include "potential.hpp"
// #include "DRalgo_xSM_nnlo.hpp"
#include "models/SM_parameters.hpp"

namespace EffectivePotential {

class DR_xsm: public Potential {
  public :

    DR_xsm(double Ms_, double lamHS_input_, double lamS_input_) {
    
      bool print = false;

      // First set inputs
      Ms = Ms_;
      lamHS_input = lamHS_input_;
      lamS_input = lamS_input_;
      if(print){
        std::cout << "Input Ms (physical scalar mass) = " << Ms << std::endl;
        std::cout << "Input lamHS = " << lamHS_input << std::endl;
        std::cout << "Input lamS = " << lamS_input << std::endl;
      }

      // Then, solve for mphiSq_input, msSq_input, lamH_input
      set_Params();
      
      if(print){
        std::cout << "msSq = " << msSq_input << std::endl;
        std::cout << "mphiSq = " << mphiSq_input << std::endl;
        std::cout << "lamH = " << lamH_input << std::endl;
      }

      // Check points
      if( !check_points(print) ) { 
        std::cout << "Failed!" << std::endl;
        exit(EXIT_FAILURE);
      }
    
      std::vector<double> x0 = {g1sq_input, g2sq_input, g3sq_input, b1_input, mphiSq_input, msSq_input, a1_input, b3_input, lamH_input, lamHS_input, lamS_input, yt_input};
  
      solveBetas(x0, 246., 1., 5000.);

    }

    size_t get_n_scalars() const override { return 2;}

    std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
      auto phi1 = phi;
      phi1[0] = - phi[0];
      auto phi2 = phi;
      phi2[1] = - phi[1];
      return {phi1,phi2};
    };

    void set_Params() {
      
      msSq_input = pow(Ms, 2) - 0.5 * lamHS_input * pow(v, 2);
      lamH_input = 0.5 * MhSq/pow(v, 2);
      mphiSq_input = - lamH_input * pow(v, 2);

    }

    bool check_points(bool print_) {

      if ( mphiSq_input > 0){
        std::cout << "Higgs mass (mphiSq) is positive." << std::endl;
        return false;
      }

      if ( msSq_input > 0){
        std::cout << "Scalar mass (msSq) is positive." << std::endl;
        return false;
      }

      if(print_) {
        std::cout << "- pow(mphiSq_input,2)/lamH_input = " << - pow(mphiSq_input,2)/lamH_input << std::endl;
        std::cout << "- pow(msSq_input,2)/lamS_input = " << - pow(msSq_input,2)/lamS_input << std::endl;
      }

      return true;
    }

    double get_TC_from_expression() const {
      const double cs = 1. / 12. * (2. * lamHS_input + 3. * lamS_input);
      const double ch = 1. / 48. * (9. * g2sq_input + 3. * g1sq_input + 12. * yt_input*yt_input + 24. * lamH_input + 2. * lamHS_input);
      const double TC_sq = -(lamS_input * ch * mphiSq_input - lamH_input * cs * msSq_input + std::sqrt(lamS_input * lamH_input) * std::abs(cs * mphiSq_input - ch * msSq_input)) / (lamS_input * square(ch) - lamH_input * square(cs));
      return std::sqrt(TC_sq);
      }

    double V(Eigen::VectorXd phi, double T) const override {
    
      const std::vector<double> par = get_3d_parameters(T);
    
      std::complex<double> g1sq(par[0],0);
      std::complex<double> g2sq(par[1],0);
      std::complex<double> g3sq(par[2],0);
      std::complex<double> b1(par[3],0);
      std::complex<double> mphiSq(par[4],0);
      std::complex<double> msSq(par[5],0);
      std::complex<double> a1(par[6],0);
      std::complex<double> b3(par[7],0);
      std::complex<double> lamH(par[8],0);
      std::complex<double> lamHS(par[9],0);
      std::complex<double> lamS(par[10],0);
    
      std::complex<double> h(phi[0]/sqrt(T + 1e-15),0);
      std::complex<double> rs(phi[1]/sqrt(T + 1e-15),0);

      std::complex<double> veffLO = 0.25*pow(h,4)*lamH + 0.5*pow(h,2)*mphiSq + 0.25*a1*pow(h,2)*rs + 0.25*pow(h,2)*lamHS*pow(rs,2) + 0.5*msSq*pow(rs,2) + 0.3333333333333333*b3*pow(rs,3) + 0.25*lamS*pow(rs,4);

      if ( potential_flag == 0 ) {
        return T * veffLO.real();
      }

      std::complex<double> veffNLO = -1./(12. * M_PI) * ( 0. 
        + 3. * pow( pow(h,2)*lamH + mphiSq + 0.5*a1*rs + 0.5*lamHS*pow(rs,2) , 3./2.)
        + pow( 0.25*(6.*pow(h,2)*lamH + pow(h,2)*lamHS + 2.*mphiSq + 2.*msSq + a1*rs + 4.*b3*rs + lamHS*pow(rs,2) + 6.*lamS*pow(rs,2) - sqrt(pow(-6.*pow(h,2)*lamH - 1.*pow(h,2)*lamHS - 2.*mphiSq - 2.*msSq - 1.*a1*rs - 4.*b3*rs - 1.*lamHS*pow(rs,2) - 6.*lamS*pow(rs,2),2) - 4.*(-1.*pow(a1,2)*pow(h,2) + 6.*pow(h,4)*lamH*lamHS + 2.*pow(h,2)*lamHS*mphiSq + 12.*pow(h,2)*lamH*msSq + 4.*mphiSq*msSq + 24.*b3*pow(h,2)*lamH*rs - 3.*a1*pow(h,2)*lamHS*rs + 8.*b3*mphiSq*rs + 2.*a1*msSq*rs + 4.*a1*b3*pow(rs,2) - 3.*pow(h,2)*pow(lamHS,2)*pow(rs,2) + 36.*pow(h,2)*lamH*lamS*pow(rs,2) + 12.*lamS*mphiSq*pow(rs,2) + 2.*lamHS*msSq*pow(rs,2) + 4.*b3*lamHS*pow(rs,3) + 6.*a1*lamS*pow(rs,3) + 6.*lamHS*lamS*pow(rs,4)))), 3./2.)
        + pow( 0.25*(6.*pow(h,2)*lamH + pow(h,2)*lamHS + 2.*mphiSq + 2.*msSq + a1*rs + 4.*b3*rs + lamHS*pow(rs,2) + 6.*lamS*pow(rs,2) + sqrt(pow(-6.*pow(h,2)*lamH - 1.*pow(h,2)*lamHS - 2.*mphiSq - 2.*msSq - 1.*a1*rs - 4.*b3*rs - 1.*lamHS*pow(rs,2) - 6.*lamS*pow(rs,2),2) - 4.*(-1.*pow(a1,2)*pow(h,2) + 6.*pow(h,4)*lamH*lamHS + 2.*pow(h,2)*lamHS*mphiSq + 12.*pow(h,2)*lamH*msSq + 4.*mphiSq*msSq + 24.*b3*pow(h,2)*lamH*rs - 3.*a1*pow(h,2)*lamHS*rs + 8.*b3*mphiSq*rs + 2.*a1*msSq*rs + 4.*a1*b3*pow(rs,2) - 3.*pow(h,2)*pow(lamHS,2)*pow(rs,2) + 36.*pow(h,2)*lamH*lamS*pow(rs,2) + 12.*lamS*mphiSq*pow(rs,2) + 2.*lamHS*msSq*pow(rs,2) + 4.*b3*lamHS*pow(rs,3) + 6.*a1*lamS*pow(rs,3) + 6.*lamHS*lamS*pow(rs,4)))), 3./2.)
      )
      - 1./(6. * M_PI) * ( 2. * pow( 0.25*g2sq*pow(h,2) , 3./2.) + pow( 0.25*(g1sq + g2sq)*pow(h,2), 3./2.))
      + 1./(12. * M_PI) * ( 4. * pow( mphiSq, 3./2.) +  pow( msSq, 3./2.)); // this ensure correct normalisation

      if ( potential_flag == 1 ) {
        return T * veffLO.real() + T * veffNLO.real();
      }

      // DR_xsm_nnlo Vnnlo;
      // std::complex<double> veffNNLO = Vnnlo.V2(phi, T, par);

      // if ( potential == 2 ) {
      // return T * veffLO.real() + T * veffNLO.real() + T * veffNNLO.real();
      // }

      if ( potential_flag == 2 ) {
      return 0.;
      }

      return 0;
    }

    void Betas( const std::vector<double>& x, std::vector<double>& dxdt, const double t) override {
    
      double g1sq = x[0];
      double g2sq = x[1];
      double g3sq = x[2];
      double b1 = x[3];
      double mphiSq = x[4];
      double msSq = x[5];
      double a1 = x[6];
      double b3 = x[7];
      double lamH = x[8];
      double lamHS = x[9];
      double lamS = x[10];
      double yt = x[11];
      dxdt[0] = 1./t * 0.25541381709839317*pow(g1sq,2);
      dxdt[1] = 1./t * 0.0612148817839124*pow(g2sq,2);
      dxdt[2] = 1./t * 0.012665147955292222*pow(g3sq,2);
      dxdt[3] = 1./t * 0.012665147955292222*(a1*mphiSq + b3*msSq);
      dxdt[4] = 1./t * 0.0031662869888230555*(pow(a1,2) + 2.*lamHS*msSq - 3.*mphiSq*(g1sq + 3.*g2sq - 8.*lamH - 4.*pow(yt,2)));
      dxdt[5] = 1./t * 0.006332573977646111*(pow(a1,2) + 4.*pow(b3,2) + 4.*lamHS*mphiSq + 6.*lamS*msSq);
      dxdt[6] = 1./t * 0.0031662869888230555*(8.*b3*lamHS + a1*(-3.*g1sq - 9.*g2sq + 24.*lamH + 8.*lamHS + 12.*pow(yt,2)));
      dxdt[7] = 1./t * 0.018997721932938333*(a1*lamHS + 6.*b3*lamS);
      dxdt[8] = 1./t * 0.0007915717472057639*(3.*pow(g1sq,2) + 9.*pow(g2sq,2) + 6.*g1sq*(g2sq - 4.*lamH) - 72.*g2sq*lamH + 4.*(48.*pow(lamH,2) + pow(lamHS,2) + 24.*lamH*pow(yt,2) - 12.*pow(yt,4)));
      dxdt[9] = 1./t * 0.0031662869888230555*lamHS*(-3.*g1sq - 9.*g2sq + 24.*lamH + 8.*lamHS + 12.*(lamS + pow(yt,2)));
      dxdt[10] = 1./t * 0.012665147955292222*(pow(lamHS,2) + 9.*pow(lamS,2));
      dxdt[11] = 1./t * 0.0005277144981371759*yt*(-17.*g1sq - 27.*g2sq - 96.*g3sq + 54.*pow(yt,2));
    
    }

    std::vector<double> get_3d_parameters(double T) const {
    
      double Gamma = scaleFactor * T;
      double scaleFactor3;
      double scaleFactor3dUS;
      double Lb = 2.*log(scaleFactor) + 2. * EulerGamma - 2. * log(4 * M_PI);
      double Lf = Lb + 4. * log(2.);
      double g1sq, g2sq, g3sq, b1, mphiSq, msSq, a1,b3, lamH, lamHS, lamS, yt;
    
      if ( running_flag ) {
        g1sq = alglib::spline1dcalc(RGEs[0], Gamma);
        g2sq = alglib::spline1dcalc(RGEs[1], Gamma);
        g3sq = alglib::spline1dcalc(RGEs[2], Gamma);
        b1 = alglib::spline1dcalc(RGEs[3], Gamma);
        mphiSq = alglib::spline1dcalc(RGEs[4], Gamma);
        msSq = alglib::spline1dcalc(RGEs[5], Gamma);
        a1 = alglib::spline1dcalc(RGEs[6], Gamma);
        b3 = alglib::spline1dcalc(RGEs[7], Gamma);
        lamH = alglib::spline1dcalc(RGEs[8], Gamma);
        lamHS = alglib::spline1dcalc(RGEs[9], Gamma);
        lamS = alglib::spline1dcalc(RGEs[10], Gamma);
        yt = alglib::spline1dcalc(RGEs[11], Gamma);
      } else {
        g1sq = g1sq_input;
        g2sq = g2sq_input;
        g3sq = g3sq_input;
        b1 = b1_input;
        mphiSq = mphiSq_input;
        msSq = msSq_input;
        a1 = a1_input;
        b3 = b3_input;
        lamH = lamH_input;
        lamHS = lamHS_input;
        lamS = lamS_input;
        yt = yt_input;
      }
      //---------------------------------------------------------------------------------
      double g1sq3d = g1sq*T - 0.0010554289962743518*pow(g1sq,2)*(Lb + 120.*Lf)*T;
      double g2sq3d = g2sq*T + 0.0010554289962743518*pow(g2sq,2)*(4. + 43.*Lb - 72.*Lf)*T;
      double g3sq3d = g3sq*T + 0.006332573977646111*pow(g3sq,2)*(1. + 11.*Lb - 12.*Lf)*T;
      double a13d = 0.0015831434944115277*sqrt(T)*(-8.*b3*lamHS*Lb + a1*(631.6546816697189 + 3.*g1sq*Lb + 9.*g2sq*Lb - 8.*(3.*lamH + lamHS)*Lb - 12.*Lf*pow(yt,2)));
      double b33d = 0.5*(2.*b3 - 0.018997721932938333*(a1*lamHS + 6.*b3*lamS)*Lb)*sqrt(T);
      double lamH3d = 0.00039578587360288194*T*(2526.6187266788756*lamH + (pow(g1sq,2) + 2.*g1sq*g2sq + 3.*pow(g2sq,2))*(2. - 3.*Lb) - 4.*(48.*pow(lamH,2) + pow(lamHS,2))*Lb + 48.*Lf*pow(yt,4) + 24.*lamH*(g1sq*Lb + 3.*g2sq*Lb - 4.*Lf*pow(yt,2)));
      double lamHS3d = 0.0015831434944115277*lamHS*T*(631.6546816697189 + (3.*g1sq + 9.*g2sq - 4.*(6.*lamH + 2.*lamHS + 3.*lamS))*Lb - 12.*Lf*pow(yt,2));
      double lamS3d = lamS*T - 0.006332573977646111*(pow(lamHS,2) + 9.*pow(lamS,2))*Lb*T;
      //---------------------------------------------------------------------------------
      double b13d_LO = (0.08333333333333333*(12.*b1 + (a1 + b3)*pow(T,2)))/sqrt(T);
      double b13d_NLO = (0.00013192862453429398*(-4.*b3*(12.*Lb*msSq - 28.89405671404654*lamS*pow(T,2) + (2.*lamHS - 3.*lamS)*Lb*pow(T,2)) + a1*(2.*pow(T,2)*(g1sq - 40.3410850710698*g2sq + 14.44702835702327*(-1.*g1sq + 2.*lamHS) + 3.*Lf*pow(yt,2)) - 1.*Lb*(48.*mphiSq + pow(T,2)*(9.*g1sq + 27.*g2sq + 24.*lamH - 10.*lamHS + 18.*pow(yt,2)))) + 24.*(-1.*a13d*(g1sq3d + 3.*g2sq3d - 2.*lamHS3d) + 4.*b33d*lamS3d)*sqrt(T)*log(scaleFactor3/scaleFactor)))/sqrt(T);
    
      double b13d = b13d_LO + b13d_NLO;
      //---------------------------------------------------------------------------------
      double lambdaVLL_1 = -0.24063781115055222*pow(g2sq,2)*T;
      double lambdaVLL_2 = -0.22797266319526*g2sq*g3sq*T;
      double lambdaVLL_3 = -0.22797266319526*pow(g3sq,2)*T;
      double lambdaVLL_4 = 0.;
      double lambdaVLL_5 = 0.;
      double lambdaVLL_6 = -0.07599088773175333*sqrt(g1sq)*pow(g3sq,1.5)*T;
      double lambdaVLL_7 = -0.2786332550164289*g1sq*g3sq*T;
      double lambdaVLL_8 = -0.13931662750821444*g1sq*g2sq*T;
      double lambdaVLL_9 = -1.5915869263817226*pow(g1sq,2)*T;
      double lambdaVL_1 = 0.012665147955292222*g2sq*lamHS*T;
      double lambdaVL_2 = -0.025330295910584444*g3sq*T*pow(yt,2);
      double lambdaVL_3 = 0.012665147955292222*g1sq*lamHS*T;
      double lambdaVL_4 = 0.0005277144981371759*g2sq*T*(947.4820225045784 + 3.*g1sq + 72.*lamH + g2sq*(123. + 43.*Lb - 72.*Lf) - 36.*pow(yt,2));
      double lambdaVL_5 = -0.0026041666666666665*sqrt(g1sq)*sqrt(g2sq)*T*(192. + 0.10132118364233778*(48.*lamH + g1sq*(118. - 1.*Lb) + g2sq*(60. + 43.*Lb) - 6.*g1sq*(-1. + Lf) - 9.*(12.666666666666666*g1sq + 8.*g2sq)*Lf - 12.*(-2. + Lf)*pow(yt,2) + 12.*Lf*pow(yt,2)));
      double lambdaVL_6 = 0.005208333333333333*g1sq*T*(96. + 0.10132118364233778*(9.*g2sq + 72.*lamH + g1sq*(67. - 1.*Lb - 48.*(-1. + Lf) - 66.*Lf) + 18.*(1.7777777777777777*(-2. + Lf) - 1.*Lf)*pow(yt,2) - 16.*Lf*pow(yt,2) + 2.*(-3.*g1sq*(-1. + Lf) + (-2. + Lf)*pow(yt,2))));
      double lambdaVVSL_1 = 0.006332573977646111*a1*g2sq*sqrt(T);
      double lambdaVVSL_2 = 0.006332573977646111*a1*g1sq*sqrt(T);
      //---------------------------------------------------------------------------------
      double MusqSU3_LO = 4.*g3sq*pow(T,2);
      double MusqSU2_LO = 3.8333333333333335*g2sq*pow(T,2);
      double MusqU1_LO = 5.166666666666667*g1sq*pow(T,2);
      double MusqSU3_NLO = -0.017414578438526805*g1sq*g3sq*pow(T,2) - 0.04274487434911125*g2sq*g3sq*pow(T,2) - 0.33001039755056605*pow(g3sq,2)*pow(T,2) - 0.006332573977646111*g3sq*pow(T,2)*pow(yt,2) + 0.13931662750821444*pow(g3sq,2)*pow(T,2)*log(0.07957747154594767*scaleFactor) + 0.18997721932938333*pow(g3sq,2)*pow(T,2)*log(3.141592653589793*T) - 0.18997721932938333*pow(g3sq,2)*pow(T,2)*log(scaleFactor*T);
      double MusqSU2_NLO = 0.012665147955292222*g2sq*mphiSq - 0.008707289219263403*g1sq*g2sq*pow(T,2) - 0.06816774467240017*pow(g2sq,2)*pow(T,2) - 0.11398633159763*g2sq*g3sq*pow(T,2) + 0.006332573977646111*g2sq*lamH*pow(T,2) + 0.0005277144981371759*g2sq*lamHS*pow(T,2) - 0.0015831434944115277*g2sq*pow(T,2)*pow(yt,2) + 0.030255631226531417*pow(g2sq,2)*pow(T,2)*log(scaleFactor) + 0.26491267806486235*pow(g2sq,2)*pow(T,2)*log(T) - 0.26491267806486235*pow(g2sq,2)*pow(T,2)*log(scaleFactor*T);
      double MusqU1_NLO = 0.012665147955292222*g1sq*mphiSq - 0.14268965380560672*pow(g1sq,2)*pow(T,2) - 0.026121867657790208*g1sq*g2sq*pow(T,2) - 0.13931662750821444*g1sq*g3sq*pow(T,2) + 0.006332573977646111*g1sq*lamH*pow(T,2) + 0.0005277144981371759*g1sq*lamHS*pow(T,2) - 0.005804859479508935*g1sq*pow(T,2)*pow(yt,2) - 0.0527714498137176*pow(g1sq,2)*pow(T,2)*log(scaleFactor) + 1.2665147955292222*pow(g1sq,2)*pow(T,2)*log(3.141592653589793*T) + 0.0003518096654247839*pow(g1sq,2)*pow(T,2)*log(12.566370614359172*T) - 1.2668666051946469*pow(g1sq,2)*pow(T,2)*log(scaleFactor*T);
    
      double MusqSU3 = MusqSU3_LO/3 + MusqSU3_NLO/3;
      double MusqSU2 = MusqSU2_LO/3 + MusqSU2_NLO/3;
      double MusqU1 = MusqU1_LO/3 + MusqU1_NLO/3;
      //---------------------------------------------------------------------------------
      scaleFactor3 = sqrt(MusqU1);
      scaleFactor3dUS = sqrt(g1sq) * scaleFactor3;
      //---------------------------------------------------------------------------------
      double mphiSq3d_LO = mphiSq + 0.020833333333333332*pow(T,2)*(3.*g1sq + 9.*g2sq + 24.*lamH + 2.*lamHS + 12.*pow(yt,2));
      double msSq3d_LO = msSq + 0.08333333333333333*(2.*lamHS + 3.*lamS)*pow(T,2);
      double mphiSq3d_NLO = 0.000021988104089048995*(-72.*pow(a1,2)*Lb + 216.*g1sq*Lb*mphiSq + 648.*g2sq*Lb*mphiSq - 1728.*lamH*Lb*mphiSq - 144.*lamHS*Lb*msSq + 272.6937977487443*pow(g1sq,2)*pow(T,2) + 596.116276066047*g1sq*g2sq*pow(T,2) - 1038.313945378327*pow(g2sq,2)*pow(T,2) - 968.1860417056755*g1sq*lamH*pow(T,2) - 2904.558125117026*g2sq*lamH*pow(T,2) + 4160.744166822702*pow(lamH,2)*pow(T,2) + 173.36434028427922*pow(lamHS,2)*pow(T,2) - 507.*pow(g1sq,2)*Lb*pow(T,2) + 216.*g1sq*g2sq*Lb*pow(T,2) - 1395.*pow(g2sq,2)*Lb*pow(T,2) - 216.*g1sq*lamH*Lb*pow(T,2) - 648.*g2sq*lamH*Lb*pow(T,2) + 9.*g1sq*lamHS*Lb*pow(T,2) + 27.*g2sq*lamHS*Lb*pow(T,2) - 72.*lamH*lamHS*Lb*pow(T,2) + 12.*pow(lamHS,2)*Lb*pow(T,2) - 36.*lamHS*lamS*Lb*pow(T,2) + 180.*pow(g1sq,2)*Lf*pow(T,2) + 324.*pow(g2sq,2)*Lf*pow(T,2) - 864.*Lf*mphiSq*pow(yt,2) - 66.*g1sq*pow(T,2)*pow(yt,2) - 54.*g2sq*pow(T,2)*pow(yt,2) - 576.*g3sq*pow(T,2)*pow(yt,2) + 47.*g1sq*Lb*pow(T,2)*pow(yt,2) + 189.*g2sq*Lb*pow(T,2)*pow(yt,2) - 192.*g3sq*Lb*pow(T,2)*pow(yt,2) - 648.*lamH*Lb*pow(T,2)*pow(yt,2) + 55.*g1sq*Lf*pow(T,2)*pow(yt,2) - 27.*g2sq*Lf*pow(T,2)*pow(yt,2) + 768.*g3sq*Lf*pow(T,2)*pow(yt,2) - 216.*lamH*Lf*pow(T,2)*pow(yt,2) - 36.*lamHS*Lf*pow(T,2)*pow(yt,2) + 108.*Lb*pow(T,2)*pow(yt,4) + 18.*log(scaleFactor3/scaleFactor)*(5.*pow(g1sq3d,2) - 39.*pow(g2sq3d,2) + 6.*g1sq3d*(3.*g2sq3d - 8.*lamH3d) - 48.*g2sq3d*(3.*lamH3d + 2.*lambdaVL_4) + 8.*(24.*pow(lamH3d,2) + pow(lamHS3d,2) - 48.*g3sq3d*lambdaVL_2 + 8.*pow(lambdaVL_2,2) + 3.*pow(lambdaVL_4,2) + 6.*pow(lambdaVL_5,2) + pow(lambdaVL_6,2))));
      double msSq3d_NLO = -0.00026385724906858796*(12.*pow(a1,2)*Lb + 48.*pow(b3,2)*Lb + 48.*lamHS*Lb*mphiSq + 72.*lamS*Lb*msSq + 26.894056714046535*g1sq*lamHS*pow(T,2) + 80.68217014213961*g2sq*lamHS*pow(T,2) - 57.78811342809307*pow(lamHS,2)*pow(T,2) - 173.36434028427922*pow(lamS,2)*pow(T,2) + 9.*g1sq*lamHS*Lb*pow(T,2) + 27.*g2sq*lamHS*Lb*pow(T,2) + 24.*lamH*lamHS*Lb*pow(T,2) - 10.*pow(lamHS,2)*Lb*pow(T,2) + 12.*lamHS*lamS*Lb*pow(T,2) - 18.*pow(lamS,2)*Lb*pow(T,2) + 18.*lamHS*Lb*pow(T,2)*pow(yt,2) - 6.*lamHS*Lf*pow(T,2)*pow(yt,2) + 12.*log(scaleFactor3/scaleFactor)*(2.*g1sq3d*lamHS3d - 4.*pow(lamHS3d,2) - 12.*pow(lamS3d,2) - 3.*pow(lambdaVL_1,2) + 6.*g2sq3d*(lamHS3d + 2.*lambdaVL_1) - 1.*pow(lambdaVL_3,2)));
    
      double mphiSq3d = mphiSq3d_LO + mphiSq3d_NLO;
      double msSq3d = msSq3d_LO + msSq3d_NLO;
      //---------------------------------------------------------------------------------
      if ( matching_flag == 0 ) {
        return {g1sq*T, g2sq*T, g3sq*T, b1/sqrt(T), mphiSq3d_LO, msSq3d_LO, a1*sqrt(T), b3*sqrt(T), lamH*T, lamHS*T, lamS*T};
      } else if ( matching_flag == 1 ) {
        return {g1sq3d, g2sq3d, g3sq3d, b13d_NLO, mphiSq3d, msSq3d, a13d, b33d, lamH3d, lamHS3d, lamS3d};
      }
      //---------------------------------------------------------------------------------
      double g1sq3dUS = g1sq3d;
      double g2sq3dUS = g2sq3d - (0.013262911924324612*pow(g2sq3d,2))/sqrt(MusqSU2);
      double g3sq3dUS = g3sq3d - (0.019894367886486918*pow(g3sq3d,2))/sqrt(MusqSU3);
      double a13dUS = -0.0008289319952702883*(48.*((3.*lambdaVL_4*lambdaVVSL_1)/sqrt(MusqSU2) + (lambdaVL_6*lambdaVVSL_2)/sqrt(MusqU1)) + a13d*(-1206.3715789784806 + (3.*pow(lambdaVVSL_1,2))/pow(MusqSU2,1.5) + pow(lambdaVVSL_2,2)/pow(MusqU1,1.5)));
      double b33dUS = -0.0024867959858108648*((36.*lambdaVL_1*lambdaVVSL_1)/sqrt(MusqSU2) + (12.*lambdaVL_3*lambdaVVSL_2)/sqrt(MusqU1) + b33d*(-402.1238596594935 + (3.*pow(lambdaVVSL_1,2))/pow(MusqSU2,1.5) + pow(lambdaVVSL_2,2)/pow(MusqU1,1.5)));
      double lamH3dUS = lamH3d - 0.009947183943243459*((8.*pow(lambdaVL_2,2))/sqrt(MusqSU3) + (3.*pow(lambdaVL_4,2))/sqrt(MusqSU2) + (4.*pow(lambdaVL_5,2))/(sqrt(MusqSU2) + sqrt(MusqU1)) + pow(lambdaVL_6,2)/sqrt(MusqU1));
      double lamHS3dUS = 0.0016578639905405765*(603.1857894892403*lamHS3d - (36.*lambdaVL_1*lambdaVL_4)/sqrt(MusqSU2) - (12.*lambdaVL_3*lambdaVL_6)/sqrt(MusqU1) + (3.*lamHS3d*pow(lambdaVVSL_1,2))/pow(MusqSU2,1.5) + (lamHS3d*pow(lambdaVVSL_2,2))/pow(MusqU1,1.5));
      double lamS3dUS = 0.003315727981081153*(301.59289474462014*lamS3d - (9.*pow(lambdaVL_1,2))/sqrt(MusqSU2) - (3.*pow(lambdaVL_3,2))/sqrt(MusqU1) + (3.*lamS3d*pow(lambdaVVSL_1,2))/pow(MusqSU2,1.5) + (lamS3d*pow(lambdaVVSL_2,2))/pow(MusqU1,1.5));
      //---------------------------------------------------------------------------------
      double b13dUS = b13d - 0.039788735772973836*(3.*sqrt(MusqSU2)*lambdaVVSL_1 + sqrt(MusqU1)*lambdaVVSL_2) - 0.0008289319952702883*b13d*((3.*pow(lambdaVVSL_1,2))/pow(MusqSU2,1.5) + pow(lambdaVVSL_2,2)/pow(MusqU1,1.5));
      //---------------------------------------------------------------------------------
      double mphiSq3dUS_LO = mphiSq3d - 0.039788735772973836*(8.*sqrt(MusqSU3)*lambdaVL_2 + 3.*sqrt(MusqSU2)*lambdaVL_4 + sqrt(MusqU1)*lambdaVL_6);
      double msSq3dUS_LO = msSq3d - 0.019894367886486918*(6.*sqrt(MusqSU2)*lambdaVL_1 + (3.*pow(lambdaVVSL_1,2))/sqrt(MusqSU2) + (2.*MusqU1*lambdaVL_3 + pow(lambdaVVSL_2,2))/sqrt(MusqU1));
      double mphiSq3dUS_NLO = 0.0007915717472057639*(48.*g3sq3d*lambdaVL_2 + 32.*log((0.5*scaleFactor3)/sqrt(MusqSU3))*(6.*g3sq3d - 1.*lambdaVL_2)*lambdaVL_2 - 16.*pow(lambdaVL_2,2) + 12.*g2sq3d*lambdaVL_4 - 6.*pow(lambdaVL_4,2) - 6.*log((0.5*scaleFactor3)/sqrt(MusqSU2))*(pow(g2sq3d,2) - 8.*g2sq3d*lambdaVL_4 + 2.*pow(lambdaVL_4,2)) - 12.*pow(lambdaVL_5,2) - 24.*log((scaleFactor3)/(sqrt(MusqSU2) + sqrt(MusqU1)))*pow(lambdaVL_5,2) - 2.*pow(lambdaVL_6,2) - 4.*log((0.5*scaleFactor3)/sqrt(MusqU1))*pow(lambdaVL_6,2) + 5.*lambdaVL_4*lambdaVLL_1 + (24.*sqrt(MusqSU2)*lambdaVL_2*lambdaVLL_2)/sqrt(MusqSU3) + (24.*sqrt(MusqSU3)*lambdaVL_4*lambdaVLL_2)/sqrt(MusqSU2) + 26.666666666666668*lambdaVL_2*lambdaVLL_3 + (8.*sqrt(MusqU1)*lambdaVL_2*lambdaVLL_7)/sqrt(MusqSU3) + (8.*sqrt(MusqSU3)*lambdaVL_6*lambdaVLL_7)/sqrt(MusqU1) + (3.*sqrt(MusqU1)*lambdaVL_4*lambdaVLL_8)/sqrt(MusqSU2) + (3.*sqrt(MusqSU2)*lambdaVL_6*lambdaVLL_8)/sqrt(MusqU1) + lambdaVL_6*lambdaVLL_9);
      double msSq3dUS_NLO = 0.00026385724906858796*(36.*g2sq3d*(1. + 4.*log((0.5*scaleFactor3)/sqrt(MusqSU2)))*lambdaVL_1 - 18.*(1. + 2.*log((0.5*scaleFactor3)/sqrt(MusqSU2)))*pow(lambdaVL_1,2) - 6.*(1. + 2.*log((0.5*scaleFactor3)/sqrt(MusqU1)))*pow(lambdaVL_3,2) + (3.*lambdaVL_1*(5.*sqrt(MusqSU2)*lambdaVLL_1 + 24.*sqrt(MusqSU3)*lambdaVLL_2 + 3.*sqrt(MusqU1)*lambdaVLL_8))/sqrt(MusqSU2) + (3.*lambdaVL_3*(8.*sqrt(MusqSU3)*lambdaVLL_7 + 3.*sqrt(MusqSU2)*lambdaVLL_8 + sqrt(MusqU1)*lambdaVLL_9))/sqrt(MusqU1) + (0.125*(3.*pow(MusqU1,1.5)*pow(lambdaVVSL_1,2) + pow(MusqSU2,1.5)*pow(lambdaVVSL_2,2))*(-50.26548245743669*msSq3d*sqrt(MusqSU2)*sqrt(MusqU1) + 6.*MusqSU2*sqrt(MusqU1)*lambdaVL_1 + 3.*sqrt(MusqU1)*pow(lambdaVVSL_1,2) + sqrt(MusqSU2)*(2.*MusqU1*lambdaVL_3 + pow(lambdaVVSL_2,2))))/(pow(MusqSU2,2)*pow(MusqU1,2)));
      double mphiSq3dUS_beta = 0.00039578587360288194*(5.*pow(g1sq3dUS,2) + 18.*g1sq3dUS*g2sq3dUS - 51.*pow(g2sq3dUS,2) - 48.*(g1sq3dUS + 3.*g2sq3dUS)*lamH3dUS + 192.*pow(lamH3dUS,2) + 8.*pow(lamHS3dUS,2));
      double msSq3dUS_beta = 0.006332573977646111*(-1.*(g1sq3dUS + 3.*g2sq3dUS - 2.*lamHS3dUS)*lamHS3dUS + 6.*pow(lamS3dUS,2));
    
      double mphiSq3dUS = mphiSq3dUS_LO + mphiSq3dUS_NLO + mphiSq3dUS_beta * log(scaleFactor3dUS / scaleFactor3);
      double msSq3dUS = msSq3dUS_LO + msSq3dUS_NLO + msSq3dUS_beta * log(scaleFactor3dUS / scaleFactor3);
      //---------------------------------------------------------------------------------
      if ( matching_flag == 2 ) {
        return {g1sq3dUS, g2sq3dUS, g3sq3dUS, b13dUS, mphiSq3dUS, msSq3dUS, a13dUS, b33dUS, lamH3dUS, lamHS3dUS, lamS3dUS};
      }
      return {0.};
    }

  private :

    PROPERTY(int, potential_flag, 1);
    PROPERTY(int, matching_flag, 1);
    PROPERTY(bool, running_flag, true);

    const double v = SM::v;
    const double Mh = SM::mh; // Captial for physical mass.
    const double MhSq = Mh*Mh;
    const double g1 = SM::gp;
    const double g2 = SM::g;

    double Ms;
  
    double g1sq_input = g1*g1;
    double g2sq_input = g2*g2;
    double g3sq_input = 0.1183 * 0.1183;
    double b1_input = 0;
    double mphiSq_input;
    double msSq_input;
    double a1_input = 0;
    double b3_input = 0;
    double lamH_input;
    double lamHS_input;
    double lamS_input;
    double yt_input = SM::yt;
  
    const double scaleFactor = M_PI;

};

} // namespace EffectivePotential

#endif
