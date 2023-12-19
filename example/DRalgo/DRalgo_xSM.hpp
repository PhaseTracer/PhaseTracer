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

#ifndef POTENTIAL_DRALGO_XSM_MODEL_HPP_INCLUDED
#define POTENTIAL_DRALGO_XSM_MODEL_HPP_INCLUDED

/**
   The xSMin  DRalgo
   https://arxiv.org/pdf/22xx.xxxxx.pdf

  Sqrt -> sqrt, pow -> pow, 2 -> 2.
  M_PI -> M_PI, log-> log
   /.{\[CurlyPhi]^2 -> Hsq, \[CurlyPhi]^4 -> Hsq^2, s^2->Ssq, s^4->Ssq^2,lambdaH->lamH,lambdaHS->lamHS, lambdaS->lamS,\[Mu]Ssq->muSsq,\[Mu]Hsq->muHsq,g1^2 -> g1sq, g2^2 -> g2sq, g13d^2->g13dsq,g23d^2->g23dsq,g33d^2->g33dsq,lambdaH3d->lamH3d, lambdaHS3d->lamHS3d,lambdaS3d->lamS3d}
 

*/

#include <vector>
#include <cmath>
#include <interpolation.h>

#include "pow.hpp"
#include "potential.hpp"
#include "pow.hpp"

namespace EffectivePotential {

class DR_xSM: public Potential {
 public:

  DR_xSM(double lamdaHS_input_,
         double lamdaS_input_,
         double mS_input_){
      lamdaHS_input = lamdaHS_input_;
      lamdaS_input = lamdaS_input_;
      mS_input = mS_input_;
      ms_sq = mS_input * mS_input;
      muSsq_init = -ms_sq + 0.5 * lamdaHS_input * vev_higgs_sq;
<<<<<<< HEAD
      solveBetas({g1sq_init, g2sq_init, g3sq_init, lamdaHS_input, lamdaH_init, lamdaS_input, Yt_init, muHsq_init, muSsq_init},246.);
    std::cout << " Before RGE g1sq = " << g1sq_init << std::endl;
    std::cout << " Before RGE g2sq = " << g2sq_init << std::endl;
    std::cout << " Before RGE g3sq = " << g3sq_init << std::endl;
    std::cout << " Before RGE Yt   = " << Yt_init << std::endl;
    std::cout << " Before RGE muSsq = " << muSsq_init << std::endl;
    std::cout << " Before RGE muHsq = " << muHsq_init << std::endl;
    std::cout << " Before RGE lamdaH = " << lamdaH_init << std::endl;
    std::cout << " Before RGE lamdaS = " << lamdaS_input << std::endl;
    std::cout << " Before RGE lamdaHS = " << lamdaHS_input << std::endl;
=======
      solveBetas({g1sq_init, g2sq_init, g3sq_init, lamdaHS_input, lamdaH_init, lamdaS_input, Yt_init, muHsq_init, muSsq_init},246.); // The order is same to `Betas'
//    std::cout << " Before RGE muHsq = " << muHsq_init << std::endl;
>>>>>>> 749ac667f7180d29d07853feac7c185fa63958c2
//    std::cout << " Before RGE muSsq = " << muSsq_init << std::endl;
    }
  
  size_t get_n_scalars() const override {return 2;}

  std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
    auto phi1 = phi;
    phi1[0] = - phi[0];
    auto phi2 = phi;
    phi2[1] = - phi[1];
    return {phi1,phi2};
  };

  double V(Eigen::VectorXd phi, double T) const override {

    const std::vector<double> par = DRstep(T);

    std::complex<double> g1sq(par[0],0);
    std::complex<double> g2sq(par[1],0);
    std::complex<double> lamH(par[2],0);
    std::complex<double> lamHS(par[3],0);
    std::complex<double> lamS(par[4],0);
    std::complex<double> muHsq(par[5],0);
    std::complex<double> muSsq(par[6],0);

//    std::cout << std::endl;
//    std::cout << "par = " << par[0] << ", " <<  par[1] << ", " << par[2] << std::endl;
//    std::cout << "par = " << par[3] << ", " <<  par[4] << ", " << par[5] << std::endl;
//    std::cout << "par = " << par[6] << std::endl;

//        std::cout << " 3dUS g1sq = " << g1sq << std::endl;
//        std::cout << " 3dUS g2sq = " << g2sq << std::endl;
//        std::cout << " 3dUS lamHS = " << lamHS << std::endl;
//        std::cout << " 3dUS lamH = " << lamH << std::endl;
//        std::cout << " 3dUS lamS = " << lamS << std::endl;
//        std::cout << " 3dUS muHsq = " << muHsq << std::endl;
//        std::cout << " 3dUS muSsq = " << muSsq << std::endl;
    
    std::complex<double> Hsq(phi[0] * phi[0] ,0);
    std::complex<double> Ssq(phi[1] * phi[1] ,0);
    
    std::complex<double> veffLO = (pow(Hsq,2.)*lamH)/4. - (Hsq*muHsq)/2.
             + (Hsq*lamHS*Ssq)/4. - (muSsq*Ssq)/2. + (lamS*pow(Ssq,2.))/4.;
    
    std::complex<double> veffNLO = 2.*(-0.020833333333333332*pow(g2sq*Hsq,1.5)/M_PI - pow((g1sq +
              g2sq)*Hsq,1.5)/(96.*M_PI)) - pow(Hsq*lamH - muHsq + (lamHS*Ssq)/2.,1.5)/(4.*M_PI) -
              pow(6.*Hsq*lamH + Hsq*lamHS - 2.*muHsq - 2.*muSsq + lamHS*Ssq + 6.*lamS*Ssq -
              sqrt(pow(Hsq*(-6.*lamH + lamHS) + 2.*muHsq - 2.*muSsq,2.) +
              2.*(7.*Hsq*pow(lamHS,2.) + 12.*lamS*(-3.*Hsq*lamH + muHsq - muSsq) +
              2.*lamHS*(3.*Hsq*(lamH + lamS) - muHsq + muSsq))*Ssq +
              pow(lamHS - 6.*lamS,2.)*pow(Ssq,2.)),1.5)/(96.*M_PI) -
              pow(6.*Hsq*lamH + Hsq*lamHS - 2.*muHsq - 2.*muSsq + lamHS*Ssq + 6.*lamS*Ssq +
              sqrt(pow(Hsq*(-6.*lamH + lamHS) + 2.*muHsq - 2.*muSsq,2.) +
              2.*(7.*Hsq*pow(lamHS,2.) + 12.*lamS*(-3.*Hsq*lamH + muHsq - muSsq) +
              2.*lamHS*(3.*Hsq*(lamH + lamS) - muHsq + muSsq))*Ssq +
              pow(lamHS - 6.*lamS,2.)*pow(Ssq,2.)),1.5)/(96.*M_PI);
    
//    std::cout << "veffLO = " << veffLO.real() << std::endl;
//    std::cout << "veffNLO = " << veffNLO.real() << std::endl;
    
    return veffLO.real() + veffNLO.real();
  }
  
  void Betas(const std::vector<double>& x, std::vector<double>& dxdt, const double t) override {
    double g1sq = x[0];
    double g2sq = x[1];
    double g3sq = x[2];
    double lamHS = x[3];
    double lamH = x[4];
    double lamS = x[5];
    double yt1 = x[6];
    double muHsq = x[7];
    double muSsq = x[8];
    dxdt[0] = 1./t*(43*pow(g1sq,2))/(144.*pow(M_PI,2));
    dxdt[1] = 1./t*(-35*pow(g2sq,2))/(48.*pow(M_PI,2));
    dxdt[2] = 1./t*(-29*pow(g3sq,2))/(24.*pow(M_PI,2));
    dxdt[3] = 1./t*(lamHS*(-3*g1sq - 9*g2sq + 4*(6*lamH + 2*lamHS + 3*lamS + 3*pow(yt1,2))))/(32.*pow(M_PI,2));
    dxdt[4] = 1./t*(3*pow(g1sq,2) + 9*pow(g2sq,2) + 6*g1sq*(g2sq - 4*lamH) - 72*g2sq*lamH + 4*pow(lamHS,2) -
               48*pow(yt1,4) + 96*lamH*(2*lamH + pow(yt1,2)))/(128.*pow(M_PI,2));
    dxdt[5] = 1./t*(pow(lamHS,2) + 9*pow(lamS,2))/(8.*pow(M_PI,2));
    dxdt[6] = 1./t*(yt1*(-17*g1sq - 27*g2sq - 96*g3sq + 54*pow(yt1,2)))/(192.*pow(M_PI,2));
    dxdt[7] = 1./t*(2*lamHS*muSsq - 3*muHsq*(g1sq + 3*g2sq - 4*(2*lamH + pow(yt1,2))))/(32.*pow(M_PI,2));
    dxdt[8] = 1./t*(2*lamHS*muHsq + 3*lamS*muSsq)/(8.*pow(M_PI,2));
  }
  
  std::vector<double> DRstep(double T) const {
    
    double Gamma = scaleFactor * T;
    double Lb = 2.*log(scaleFactor) + 2. * EulerGamma - 2. * log(4 * M_PI);
    double Lf = Lb + 4. * log(2.);
    
    double g1sq = alglib::spline1dcalc(RGEs[0], Gamma);
    double g2sq = alglib::spline1dcalc(RGEs[1], Gamma);
    double g3sq = alglib::spline1dcalc(RGEs[2], Gamma);

    double lamHS = alglib::spline1dcalc(RGEs[3], Gamma);
    double lamH  = alglib::spline1dcalc(RGEs[4], Gamma);
    double lamS  = alglib::spline1dcalc(RGEs[5], Gamma);
    
    double yt1   = alglib::spline1dcalc(RGEs[6], Gamma);
    double muHsq = alglib::spline1dcalc(RGEs[7], Gamma);
    double muSsq = alglib::spline1dcalc(RGEs[8], Gamma);
//    std::cout << std::endl;
//    std::cout << " Gamma = " << Gamma << std::endl;
//    std::cout << " After RGE g1sq = " << g1sq << std::endl;
//    std::cout << " After RGE g2sq = " << g2sq << std::endl;
//    std::cout << " After RGE g3sq = " << g3sq << std::endl;
//    std::cout << " After RGE lamHS = " << lamHS << std::endl;
//    std::cout << " After RGE lamH = " << lamH << std::endl;
//    std::cout << " After RGE lamS = " << lamS << std::endl;
//    std::cout << " After RGE yt1 = " << yt1<< std::endl;
//    std::cout << " After RGE muHsq = " << muHsq << std::endl;
//    std::cout << " After RGE muSsq = " << muSsq << std::endl;
    
    // ---------------------------------------------------------------------------------
    // Couplings
    double g13dsq = g1sq*T - (pow(g1sq,2)*(3*Lb + 40*Lf)*T)/(288.*pow(M_PI,2));
    double g23dsq = g2sq*T + (pow(g2sq,2)*(4 + 43*Lb - 8*Lf)*T)/(96.*pow(M_PI,2));
    double g33dsq = g3sq*T + (pow(g3sq,2)*(3 + 33*Lb - 4*Lf)*T)/(48.*pow(M_PI,2));
    
//    std::cout << " g13dsq = " << g13dsq << std::endl;
//    std::cout << " g23dsq = " << g23dsq << std::endl;
//    std::cout << " g33dsq = " << g33dsq << std::endl;
    
    double lamH3d = (T*(pow(g2sq,2)*(6 - 9*Lb) + pow(g1sq,2)*(2 - 3*Lb) + 72*g2sq*lamH*Lb - 4*pow(lamHS,2)*Lb +\
                    g1sq*(g2sq*(4 - 6*Lb) + 24*lamH*Lb) +\
                    48*Lf*pow(yt1,4) - 32*lamH*(6*lamH*Lb - 8*pow(M_PI,2) + 3*Lf*pow(yt1,2))))/(256.*pow(M_PI,2));
    double lamHS3d = (lamHS*T*((3*g1sq + 9*g2sq - 4*(6*lamH + 2*lamHS + 3*lamS))*Lb + 64*pow(M_PI,2) - 12*Lf*pow(yt1,2)))/(64.*pow(M_PI,2));
    double lamS3d = lamS*T - ((pow(lamHS,2) + 9*pow(lamS,2))*Lb*T)/(16.*pow(M_PI,2));
    // ---------------------------------------------------------------------------------
    // TemporalScalarCouplings
    double lambdaVLL_1 = (-353*pow(g1sq,2)*T)/(216.*pow(M_PI,2));
    double lambdaVLL_2 = -0.041666666666666664*(g1sq*g2sq*T)/pow(M_PI,2);
    double lambdaVLL_3 = (13*pow(g2sq,2)*T)/(24.*pow(M_PI,2));
    double lambdaVLL_4 = (-11*g1sq*g3sq*T)/(36.*pow(M_PI,2));
    double lambdaVLL_5 = -0.25*(g2sq*g3sq*T)/pow(M_PI,2);
    double lambdaVLL_6 = -0.08333333333333333*(sqrt(g1sq)*pow(sqrt(g3sq),3)*T)/pow(M_PI,2);
    double lambdaVLL_7 = (7*pow(g3sq,2)*T)/(12.*pow(M_PI,2));
    double lambdaVL_1 = -0.25*(g3sq*T*pow(yt1,2))/pow(M_PI,2);
    double lambdaVL_2 = -0.0008680555555555555*(sqrt(g1sq)*sqrt(g2sq)*T*(-3*g2sq*(-4 + 43*Lb - 8*Lf) + g1sq*(-52 + 3*Lb + 40*Lf) - 72*(2*lamH + 8*pow(M_PI,2) + pow(yt1,2))))/ pow(M_PI,2);
    double lambdaVL_3 = (g2sq*T*(3*g1sq + g2sq*(59 + 43*Lb - 8*Lf) + 12*(6*lamH + 8*pow(M_PI,2) - 3*pow(yt1,2))))/(192.*pow(M_PI,2));
    double lambdaVL_4 = (g1sq*T*(g1sq*(43 - 3*Lb - 40*Lf) + 3*(9*g2sq + 72*lamH + 96*pow(M_PI,2) - 68*pow(yt1,2))))/(576.*pow(M_PI,2));
    double lambdaVL_5 = (g1sq*lamHS*T)/(8.*pow(M_PI,2));
    double lambdaVL_6 = (g2sq*lamHS*T)/(8.*pow(M_PI,2));
    // ---------------------------------------------------------------------------------
    // 3d DebyeMass
    double MusqSU2_LO = (7*g2sq*pow(T,2))/6.;
    double MusqSU3_LO = (4*g3sq*pow(T,2))/3.;
    double MusqU1_LO = (13*g1sq*pow(T,2))/18.;
    
    double MusqSU2_NLO = (g2sq*(-3*g1sq*pow(T,2) + g2sq*(283 + 602*Lb - 112*Lf)*pow(T,2) -
                      6*(24*muHsq + 24*g3sq*pow(T,2) + pow(T,2)*(-12*lamH - lamHS +
                      3*pow(yt1,2)))))/(1152.*pow(M_PI,2));
    double MusqSU3_NLO = -0.001736111111111111*(g3sq*pow(T,2)*(11*g1sq + 27*g2sq -
                      16*g3sq*(13 + 33*Lb - 4*Lf) + 36*pow(yt1,2)))/pow(M_PI,2);
    double MusqU1_NLO = -0.000048225308641975306*(g1sq*(2592*muHsq + pow(T,2)*
                    (2*g1sq*(175 + 78*Lb + 1040*Lf) + 18*(9*g2sq + 176*g3sq - 72*lamH - 6*lamHS + 66*pow(yt1,2)))))/pow(M_PI,2);
    
    double MusqSU2 = MusqSU2_LO + MusqSU2_NLO;
    double MusqSU3 = MusqSU3_LO + MusqSU3_NLO;
    double MusqU1 = MusqU1_LO + MusqU1_NLO;
    
    // ---------------------------------------------------------------------------------
    // 3d ScalarMass
    double MuHsq3d_LO = (-3*g1sq*pow(T,2) - 9*g2sq*pow(T,2) - 2*(-24*muHsq + 12*lamH*pow(T,2) + lamHS*pow(T,2) + 6*pow(T,2)*pow(yt1,2)))/48.;
    double MuSsq3d_LO = muSsq - ((2*lamHS + 3*lamS)*pow(T,2))/12.;

//        std::cout << " MuHsq3d_LO = " << MuHsq3d_LO << std::endl;
//        std::cout << " MuSsq3d_LO = " << MuSsq3d_LO << std::endl;
//    std::cout << " g1sq = " << g1sq << std::endl;
//    std::cout << " g2sq = " << g2sq << std::endl;
//    std::cout << " g3sq = " << g3sq << std::endl;
//    std::cout << " Lb = " << Lb << std::endl;
//    std::cout << " Lf = " << Lf << std::endl;
//    std::cout << " ====================== "<< std::endl;
//    std::cout << " muHsq = " << muHsq << std::endl;
//    std::cout << " muSsq = " << muSsq << std::endl;
//    std::cout << " lamHS = " << lamHS << std::endl;
//    std::cout << " lamH = " << lamH << std::endl;
//    std::cout << " lamS = " << lamS << std::endl;
//    std::cout << " yt1 = " << yt1 << std::endl;
//    std::cout << " ====================== "<< std::endl;
//    std::cout << " g13dsq = " << g13dsq << std::endl;
//    std::cout << " g23dsq = " << g23dsq << std::endl;
//    std::cout << " g33dsq = " << g33dsq << std::endl;
//    std::cout << " lamH3d = " << lamH3d << std::endl;
//    std::cout << " lamHS3d = " << lamHS3d << std::endl;
//    std::cout << " lamS3d = " << lamS3d << std::endl;
//    std::cout << " ====================== "<< std::endl;
//    std::cout << " lambdaVL_1 = " << lambdaVL_1 << std::endl;
//    std::cout << " lambdaVL_2 = " << lambdaVL_2 << std::endl;
//    std::cout << " lambdaVL_3 = " << lambdaVL_3 << std::endl;
//    std::cout << " lambdaVL_4 = " << lambdaVL_4 << std::endl;
//    std::cout << " lambdaVL_5 = " << lambdaVL_5 << std::endl;
//    std::cout << " lambdaVL_6 = " << lambdaVL_6 << std::endl;
//
//    std::cout << " ====================== "<< std::endl;
//    std::cout << " lambdaVLL_1 = " << lambdaVLL_1 << std::endl;
//    std::cout << " lambdaVLL_2 = " << lambdaVLL_2 << std::endl;
//    std::cout << " lambdaVLL_3 = " << lambdaVLL_3 << std::endl;
//    std::cout << " lambdaVLL_4 = " << lambdaVLL_4 << std::endl;
//    std::cout << " lambdaVLL_5 = " << lambdaVLL_5 << std::endl;
//    std::cout << " lambdaVLL_6 = " << lambdaVLL_6 << std::endl;
//    std::cout << " lambdaVLL_7 = " << lambdaVLL_7 << std::endl;
//
//    std::cout << " ====================== "<< std::endl;
//    std::cout << " EulerGamma = " << EulerGamma << std::endl;
//    std::cout << " Glaisher = " << Glaisher << std::endl;

    double MuHsq3d_NLO = (-9*(9*g2sq*(-24*Lb*muHsq + pow(T,2)*(lamHS*Lb + (-2 + 7*Lb - Lf)*pow(yt1,2) +
              8*lamH*(1 + 6*EulerGamma - 3*Lb - 72*log(Glaisher)))) -
              4*(16*g3sq*(3 + Lb - 4*Lf)*pow(T,2)*pow(yt1,2) + 3*Lf*(-24*muHsq +
              (6*lamH + lamHS)*pow(T,2))*pow(yt1,2) -
              Lb*(12*(12*lamH*muHsq + lamHS*muSsq) + pow(T,2)*(lamHS*(-6*lamH + lamHS - 3*lamS) -
              54*lamH*pow(yt1,2) + 9*pow(yt1,4))) +
              6*(24*pow(lamH,2) + pow(lamHS,2))*pow(T,2)*(EulerGamma - 12*log(Glaisher))) +
              pow(g2sq,2)*pow(T,2)*(175 + 243*EulerGamma - 177*Lb + 12*Lf - 2916*log(Glaisher))) +
              3*g1sq*(216*Lb*muHsq - pow(T,2)*(9*lamHS*Lb + (-66 + 47*Lb + 55*Lf)*pow(yt1,2) +
              72*lamH*(1 + 6*EulerGamma - 3*Lb - 72*log(Glaisher))) +
              54*g2sq*pow(T,2)*(1 + 5*EulerGamma - 4*Lb - 60*log(Glaisher))) +
              pow(g1sq,2)*pow(T,2)*(-43 + 189*EulerGamma + 81*Lb - 60*Lf - 2268*log(Glaisher)) -
              54*log(scaleFactor3/scaleFactor)*(5*pow(g13dsq,2) - 39*pow(g23dsq,2) +
              6*g13dsq*(3*g23dsq - 8*lamH3d) - 48*g23dsq*(3*lamH3d + 2*lambdaVL_3) +
              8*(24*pow(lamH3d,2) + pow(lamHS3d,2) - 48*g33dsq*lambdaVL_1 + 8*pow(lambdaVL_1,2) +
              6*pow(lambdaVL_2,2) + 3*pow(lambdaVL_3,2) + pow(lambdaVL_4,2))))/(13824.*pow(M_PI,2));
    
    double MuSsq3d_NLO = (-48*lamHS*Lb*muHsq - 72*lamS*Lb*muSsq - 6*g2sq*lamHS*pow(T,2) -
              36*EulerGamma*g2sq*lamHS*pow(T,2) +
              24*EulerGamma*pow(lamHS,2)*pow(T,2) + 72*EulerGamma*pow(lamS,2)*pow(T,2) +
              27*g2sq*lamHS*Lb*pow(T,2) + 24*lamH*lamHS*Lb*pow(T,2) -
              10*pow(lamHS,2)*Lb*pow(T,2) + 12*lamHS*lamS*Lb*pow(T,2) - 18*pow(lamS,2)*Lb*pow(T,2) +
              18*lamHS*Lb*pow(T,2)*pow(yt1,2) -
              6*lamHS*Lf*pow(T,2)*pow(yt1,2) + 432*g2sq*lamHS*pow(T,2)*log(Glaisher) -
              288*pow(lamHS,2)*pow(T,2)*log(Glaisher) -
              864*pow(lamS,2)*pow(T,2)*log(Glaisher) + g1sq*lamHS*pow(T,2)*(-2 - 12*EulerGamma +
              9*Lb + 144*log(Glaisher)) +
              12*log(scaleFactor3 / scaleFactor)*(2*g13dsq*lamHS3d - 4*pow(lamHS3d,2) -
              12*pow(lamS3d,2) - pow(lambdaVL_5,2) - 3*pow(lambdaVL_6,2) + 6*g23dsq*(lamHS3d
              + 2*lambdaVL_6)))/(384.*pow(M_PI,2));
    
    double MuHsq3d = MuHsq3d_LO + MuHsq3d_NLO;
    double MuSsq3d = MuSsq3d_LO + MuSsq3d_NLO;
    
//    std::cout << " MuHsq3d_LO = " << MuHsq3d_LO << std::endl;
//    std::cout << " MuHsq3d_NLO = " << MuHsq3d_NLO << std::endl;
//
//    std::cout << " MuSsq3d_LO = " << MuSsq3d_LO << std::endl;
//    std::cout << " MuSsq3d_NLO = " << MuSsq3d_NLO << std::endl;
    
    // ---------------------------------------------------------------------------------
    // 3d US Couplings
//    std::cout << " lamH3d = " << lamH3d << std::endl;
//    std::cout << " (8*pow(lambdaVL_1,2))/sqrt(MusqSU3) = " << (8*pow(lambdaVL_1,2))/sqrt(MusqSU3) << std::endl;
//    std::cout << " (4*pow(lambdaVL_2,2))/(sqrt(MusqSU2) + sqrt(MusqU1)) = " << (4*pow(lambdaVL_2,2))/(sqrt(MusqSU2) + sqrt(MusqU1))<< std::endl;
//    std::cout << " (3*pow(lambdaVL_3,2))/sqrt(MusqSU2) = " << (3*pow(lambdaVL_3,2))/sqrt(MusqSU2) << std::endl;
//    std::cout << " pow(lambdaVL_4,2)/sqrt(MusqU1) = " << pow(lambdaVL_4,2)/sqrt(MusqU1) << std::endl;
    double LambdaH3dUS = lamH3d - ((8*pow(lambdaVL_1,2))/sqrt(MusqSU3) +
              (4*pow(lambdaVL_2,2))/(sqrt(MusqSU2) + sqrt(MusqU1)) +
              (3*pow(lambdaVL_3,2))/sqrt(MusqSU2) +
              pow(lambdaVL_4,2)/sqrt(MusqU1))/(32.*M_PI);
    double LambdaHS3dUS = lamHS3d - ((lambdaVL_4*lambdaVL_5)/sqrt(MusqU1) + (3*lambdaVL_3*lambdaVL_6)/sqrt(MusqSU2))/(16.*M_PI);
    double LambdaS3dUS = lamS3d - (pow(lambdaVL_5,2)/sqrt(MusqU1) + (3*pow(lambdaVL_6,2))/sqrt(MusqSU2))/(32.*M_PI);
    double g13dUSsq = g13dsq;
    double g23dUSsq = g23dsq - pow(g23dsq,2)/(24.*M_PI*sqrt(MusqSU2));
    // ---------------------------------------------------------------------------------
    // 3d US Masses
    double MuHsq3dUS_LO = MuHsq3d + (8*sqrt(MusqSU3)*lambdaVL_1 + 3*sqrt(MusqSU2)*lambdaVL_3 + sqrt(MusqU1)*lambdaVL_4)/(8.*M_PI);
    double MuSsq3dUS_LO = MuSsq3d + (sqrt(MusqU1)*lambdaVL_5 + 3*sqrt(MusqSU2)*lambdaVL_6)/(8.*M_PI);
    
//    std::cout << " MuHsq3d = " << MuHsq3d << std::endl;
    
    double MuHsq3dUS_NLO = -0.0078125*(48*g33dsq*lambdaVL_1 +
            32*log(scaleFactor3*T/(2.*sqrt(MusqSU3)))*(6*g33dsq - lambdaVL_1)*lambdaVL_1 -
            16*pow(lambdaVL_1,2) - 12*pow(lambdaVL_2,2) -
            24*log(scaleFactor3*T/(sqrt(MusqSU2) + sqrt(MusqU1)))*pow(lambdaVL_2,2) +
            12*g23dsq*lambdaVL_3 - 6*pow(lambdaVL_3,2) -
            6*log(scaleFactor3*T/(2.*sqrt(MusqSU2)))*(pow(g23dsq,2) - 8*g23dsq*lambdaVL_3 +
            2*pow(lambdaVL_3,2)) - 2*pow(lambdaVL_4,2) -
            4*log(scaleFactor3*T/(2.*sqrt(MusqU1)))*pow(lambdaVL_4,2) + lambdaVL_4*lambdaVLL_1 +
            (3*sqrt(MusqU1)*lambdaVL_3*lambdaVLL_2)/sqrt(MusqSU2) +
            (3*sqrt(MusqSU2)*lambdaVL_4*lambdaVLL_2)/sqrt(MusqU1) + 15*lambdaVL_3*lambdaVLL_3 +
            (8*sqrt(MusqU1)*lambdaVL_1*lambdaVLL_4)/sqrt(MusqSU3) +
            (8*sqrt(MusqSU3)*lambdaVL_4*lambdaVLL_4)/sqrt(MusqU1) +
            (24*sqrt(MusqSU2)*lambdaVL_1*lambdaVLL_5)/sqrt(MusqSU3) +
            (24*sqrt(MusqSU3)*lambdaVL_3*lambdaVLL_5)/sqrt(MusqSU2) +
            80*lambdaVL_1*lambdaVLL_7)/pow(M_PI,2);
      double MuSsq3dUS_NLO = -0.0078125*
            (-2*(1 + 2*log(scaleFactor3*T/(2.*sqrt(MusqU1))))*pow(lambdaVL_5,2) + 12*g23dsq*(1 +
            4*log(scaleFactor3*T/(2.*sqrt(MusqSU2))))*lambdaVL_6 -
            6*(1 + 2*log(scaleFactor3*T/(2.*sqrt(MusqSU2))))*pow(lambdaVL_6,2) +
            (lambdaVL_5*(sqrt(MusqU1)*lambdaVLL_1 + 3*sqrt(MusqSU2)*lambdaVLL_2 +
            8*sqrt(MusqSU3)*lambdaVLL_4))/
            sqrt(MusqU1) + (3*lambdaVL_6*(sqrt(MusqU1)*lambdaVLL_2 + 5*sqrt(MusqSU2)*lambdaVLL_3 +
            8*sqrt(MusqSU3)*lambdaVLL_5))/sqrt(MusqSU2))/pow(M_PI,2);
    double MuHsq3dUS_beta = (-5*pow(g13dUSsq,2) + 51*pow(g23dUSsq,2) +
                            144*g23dUSsq*LambdaH3dUS + 6*g13dUSsq*(-3*g23dUSsq + 8*LambdaH3dUS) -
                            8*(24*pow(LambdaH3dUS,2) + pow(LambdaHS3dUS,2)))/(256.*pow(M_PI,2));
    double MuSsq3dUS_beta = ((g13dUSsq + 3*g23dUSsq - 2*LambdaHS3dUS)*LambdaHS3dUS -
                            6*pow(LambdaS3dUS,2))/(16.*pow(M_PI,2));
    double MuHsq3dUS = MuHsq3dUS_LO + MuHsq3dUS_NLO + MuHsq3dUS_beta * log(scaleFactor3dUS / scaleFactor3);
    double MuSsq3dUS = MuSsq3dUS_LO + MuSsq3dUS_NLO + MuSsq3dUS_beta * log(scaleFactor3dUS / scaleFactor3);

    
//    std::cout << " LambdaH3dUS = " << LambdaH3dUS << std::endl;
//    std::cout << " LambdaS3dUS = " << LambdaS3dUS << std::endl;
//    std::cout << " LambdaHS3dUS = " << LambdaHS3dUS << std::endl;
//    std::cout << " MuHsq3dUS_LO = " << MuHsq3dUS_LO << std::endl;
//    std::cout << " MuSsq3dUS_LO = " << MuSsq3dUS_LO << std::endl;
//    std::cout << " MuHsq3dUS_beta = " << MuHsq3dUS_beta << std::endl;
//    std::cout << " scaleFactor3dUS = " << scaleFactor3dUS << std::endl;
//    std::cout << " scaleFactor3 = " << scaleFactor3 << std::endl;
//    std::cout << " scaleFactor = " << scaleFactor << std::endl;
    
    return {g13dUSsq, g23dUSsq, LambdaH3dUS, LambdaHS3dUS, LambdaS3dUS, MuHsq3dUS, MuSsq3dUS};
  };
  

 private:
  
  double lamdaHS_input = 1.6;
  double lamdaS_input = 1.0;
  double mS_input = 160;
  
  const double mZsq = 91.2*91.2;
  const double mWsq = 80.4 * 80.4;
  const double mt = 172.4;
  const double mh = 125.;
  const double mh_sq = mh*mh;
  const double vev_higgs = 246.;
  const double vev_higgs_sq = vev_higgs*vev_higgs;
  
  const double g1sq_init = 4. * mWsq / vev_higgs_sq;
  const double g2sq_init = 4. * (mZsq - mWsq) / vev_higgs_sq;
  const double g3sq_init = pow(4 * M_PI * 0.1179 / (1 - 0.1179 / (4 * M_PI)),2);
  const double Yt_init = sqrt(2.) * mt / vev_higgs;

  double ms_sq = mS_input * mS_input;
  double lamdaH_init = 0.5 * mh_sq / vev_higgs_sq;
  double muHsq_init = 0.5 * mh_sq;
  double muSsq_init = -ms_sq + 0.5 * lamdaHS_input * vev_higgs_sq;
  
  double scaleFactor = M_PI;
  double scaleFactor3 = g1sq_init;
  double scaleFactor3dUS = g1sq_init*g1sq_init;
};

}  // namespace EffectivePotential

#endif
