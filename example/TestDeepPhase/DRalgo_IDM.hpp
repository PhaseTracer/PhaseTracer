// ====================================================================
// This file was made by DRtracer for use within PhaseTracer. 
// ====================================================================

#ifndef POTENTIAL_DRALGO_IDM_HPP_INCLUDED
#define POTENTIAL_DRALGO_IDM_MODEL_HPP_INCLUDED

#include <vector>
#include <cmath>
#include <complex>
#include <Eigen/Eigenvalues>
#include <interpolation.h>

#include "pow.hpp"
#include "potential.hpp"
#include "helperIncludes/complex_operators.hpp"

namespace EffectivePotential {

class DR_idm: public Potential {
  public :

    DR_idm( std::vector<double> inputParams ) {
    
      std::vector<double> x0 = lagrangianParams(inputParams);

      std::vector<std::string> names = {"gYsq", "gWsq", "gYsq", "m11", "m22", "m12", "l1", "l2", "l3", "l4", "l5", "yt"};

      for ( int i = 0; i < x0.size(); i += 1 ) { std::cout << names[i] << " : " << x0[i] << "\n"; }
    
      std::cout << "solving beta functions" << std::endl;
      solveBetas(x0, 91.1880, 20., 1000., 1.);
      std::cout << "solved beta functions" << std::endl;
    }

    std::vector<double> lagrangianParams(std::vector<double> inputParams) {

      // some hard-coded SM params at Q = mZ
      double GF = 1.166e-5;
      double thetaEW = asin(sqrt(0.23129));
      double mW = 80.3692;
      double mZ = 91.1880;
      double alphaS = 0.1180;
      double mt = 172.56;

      // defined from above
      double v = 1./sqrt(sqrt(2)*GF); double vsq = v*v;
      double gS = sqrt(4*M_PI*alphaS); double gSsq = gS*gS;
      double gW = 2 * mW / v; double gWsq = gW*gW;
      double gY = tan(thetaEW)*gW; double gYsq = gY*gY;

      // extract inputs
      double m12 = 0.0; double m12sq = m12*m12;
      double mh = inputParams[0]; double mhsq = mh*mh;
      double mH = inputParams[1]; double mHsq = mH*mH;
      double mA = inputParams[2]; double mAsq = mA*mA;
      double mHpm = inputParams[3]; double mHpmsq = mHpm*mHpm;
      double lambda2 = inputParams[4];
      double lambdaL = inputParams[5];

      double yt = sqrt(2) * mt/v;

      double l1 = 2 * (0.5 * mhsq / vsq);
      double l2 = 2 * (lambda2);
      double l5 = (mHsq - mAsq) / vsq;
      double l4 = -2.0 * (mHpmsq - mAsq)/vsq + l5;
      double l3 = 2 * lambdaL - l4 - l5;
      
      double m11 = - l1 * vsq;
      double m22 = mHpmsq - lambdaL * vsq / 2.0;

      if( !check_quartics(l1, l2, l3, l4, l5) ) {
        throw std::runtime_error("IDM quartic couplings do not satisfy bounded-from-below conditions.");
      }

      // if( !check_masses(mH, mHpm, mA, mW, mZ) ) {
      //   throw std::runtime_error("IDM mass parameters do not satisfy tree-level unitarity conditions.");
      // }

      return {gYsq, gWsq, gSsq, m11, m22, m12, l1, l2, l3, l4, l5, yt};
    }

    bool check_quartics(const double& l1, const double& l2, const double& l3, const double& l4, const double& l5) const {
      if ( l1 < 0 ) { return false; }
      if ( l2 < 0 ) { return false; }
      if ( l3 + sqrt(l1*l2) < 0 ) { return false; }
      if ( l3 + l4 - abs(l5) + sqrt(l1*l2) < 0 ) { return false; }
      return true;
    }

    bool check_masses(const double& mH, const double& mHpm, const double& mA, const double& mW, const double& mZ) const {
      if ( mH + mHpm < mW ) { return false; }
      if ( mA + mHpm < mW ) { return false; }
      if ( mH + mA < mZ ) { return false; }
      if ( 2 * mHpm < mZ ) { return false; }
      return true;
    }

    size_t get_n_scalars() const override { return 1;}

    std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
      auto phi1 = phi;
      phi1[0] = - phi[0];
      // auto phi2 = phi;
      // phi2[0] = - phi[0];
      return {phi1};
    };

    double V(Eigen::VectorXd phi, double T) const override {
    
      const std::vector<double> par = get_3d_parameters(T);
    
      std::complex<double> gYsq(par[0],0);
      std::complex<double> gWsq(par[1],0);
      std::complex<double> gSsq(par[2],0);
      std::complex<double> m11(par[3],0);
      std::complex<double> m22(par[4],0);
      std::complex<double> m12(par[5],0);
      std::complex<double> l1(par[6],0);
      std::complex<double> l2(par[7],0);
      std::complex<double> l3(par[8],0);
      std::complex<double> l4(par[9],0);
      std::complex<double> l5(par[10],0);
    
      std::complex<double> h1(phi[0]/sqrt(T + 1e-15),0.0);
      // std::complex<double> h2(phi[1]/sqrt(T + 1e-15),0.0);
      std::complex<double> h2(0.0,0.0);

      std::complex<double> h1sq = h1*h1;
      std::complex<double> h2sq = h2*h2;
    
      std::complex<double> veffLO = (
        0.5 * m11 * h1sq 
        + 0.5 * m22 * h2sq 
        - m12 * h1 * h2 
        + 0.125 * l1 * h1sq * h1sq
        + 0.125 * l2 * h2sq * h2sq
        + 0.25 * (l3+l4+l5) * h1sq * h2sq
      );

      std::complex<double> eigen_1 = 0.25*(pow(h1,2)*(l1 + l3) + pow(h2,2)*(l2 + l3) + 2.*m11 + 2.*m22 - 1.*sqrt(4.*pow(m11,2) + 4.*pow(h1*h2*(l4 + l5) - 2.*pow(m12,2),2) - 4.*m11*(pow(h2,2)*(l2 - 1.*l3) + pow(h1,2)*(-1.*l1 + l3) + 2.*m22) + pow(pow(h2,2)*(l2 - 1.*l3) + pow(h1,2)*(-1.*l1 + l3) + 2.*m22,2)));
      std::complex<double> eigen_2 = 0.25*(pow(h1,2)*(l1 + l3) + pow(h2,2)*(l2 + l3) + 2.*m11 + 2.*m22 + sqrt(4.*pow(m11,2) + 4.*pow(h1*h2*(l4 + l5) - 2.*pow(m12,2),2) - 4.*m11*(pow(h2,2)*(l2 - 1.*l3) + pow(h1,2)*(-1.*l1 + l3) + 2.*m22) + pow(pow(h2,2)*(l2 - 1.*l3) + pow(h1,2)*(-1.*l1 + l3) + 2.*m22,2)));
      std::complex<double> eigen_3 = 0.25*(pow(h1,2)*l1 + pow(h2,2)*l2 + pow(h1,2)*l3 + pow(h2,2)*l3 + pow(h1,2)*l4 + pow(h2,2)*l4 - 1.*pow(h1,2)*l5 - 1.*pow(h2,2)*l5 + 2.*m11 + 2.*m22 - 1.*sqrt(pow(pow(h1,2)*(l1 + l3 + l4 - 1.*l5) + pow(h2,2)*(l2 + l3 + l4 - 1.*l5) + 2.*m11 + 2.*m22,2) + 4.*(pow(h1,4)*l1*l5 + 2.*pow(h1,2)*l5*m11 - 8.*h1*h2*l5*pow(m12,2) + 4.*pow(m12,4) - 1.*(pow(h1,2)*l1 + pow(h2,2)*(l3 + l4) + 2.*m11)*(pow(h2,2)*l2 + pow(h1,2)*(l3 + l4) + 2.*m22) + pow(h2,2)*l5*(pow(h2,2)*l2 + pow(h1,2)*(2.*(l3 + l4) + 3.*l5) + 2.*m22))));
      std::complex<double> eigen_4 = 0.25*(pow(h1,2)*l1 + pow(h2,2)*l2 + pow(h1,2)*l3 + pow(h2,2)*l3 + pow(h1,2)*l4 + pow(h2,2)*l4 - 1.*pow(h1,2)*l5 - 1.*pow(h2,2)*l5 + 2.*m11 + 2.*m22 + sqrt(pow(pow(h1,2)*(l1 + l3 + l4 - 1.*l5) + pow(h2,2)*(l2 + l3 + l4 - 1.*l5) + 2.*m11 + 2.*m22,2) + 4.*(pow(h1,4)*l1*l5 + 2.*pow(h1,2)*l5*m11 - 8.*h1*h2*l5*pow(m12,2) + 4.*pow(m12,4) - 1.*(pow(h1,2)*l1 + pow(h2,2)*(l3 + l4) + 2.*m11)*(pow(h2,2)*l2 + pow(h1,2)*(l3 + l4) + 2.*m22) + pow(h2,2)*l5*(pow(h2,2)*l2 + pow(h1,2)*(2.*(l3 + l4) + 3.*l5) + 2.*m22))));
      std::complex<double> eigen_5 = 0.25*(3.*pow(h1,2)*l1 + 3.*pow(h2,2)*l2 + pow(h1,2)*l3 + pow(h2,2)*l3 + pow(h1,2)*l4 + pow(h2,2)*l4 + pow(h1,2)*l5 + pow(h2,2)*l5 + 2.*m11 + 2.*m22 - 1.*sqrt(pow(pow(h1,2)*(3.*l1 + l3 + l4 + l5) + pow(h2,2)*(3.*l2 + l3 + l4 + l5) + 2.*m11 + 2.*m22,2) - 4.*(3.*pow(h1,4)*l1*l3 + 3.*pow(h2,4)*l2*l3 - 3.*pow(h1,2)*pow(h2,2)*pow(l3,2) + 3.*pow(h1,4)*l1*l4 + 3.*pow(h2,4)*l2*l4 - 6.*pow(h1,2)*pow(h2,2)*l3*l4 - 3.*pow(h1,2)*pow(h2,2)*pow(l4,2) + 3.*pow(h2,4)*l2*l5 - 3.*pow(h1,2)*pow(h2,2)*pow(l5,2) + 2.*pow(h1,2)*l3*m11 + 2.*pow(h1,2)*l4*m11 + pow(h1,2)*l5*(3.*pow(h1,2)*l1 - 6.*pow(h2,2)*(l3 + l4) + 2.*m11) + 8.*h1*h2*(l3 + l4 + l5)*pow(m12,2) - 4.*pow(m12,4) + 2.*pow(h2,2)*(l3 + l4 + l5)*m22 + (3.*pow(h1,2)*l1 + 2.*m11)*(3.*pow(h2,2)*l2 + 2.*m22))));
      std::complex<double> eigen_6 = 0.25*(3.*pow(h1,2)*l1 + 3.*pow(h2,2)*l2 + pow(h1,2)*l3 + pow(h2,2)*l3 + pow(h1,2)*l4 + pow(h2,2)*l4 + pow(h1,2)*l5 + pow(h2,2)*l5 + 2.*m11 + 2.*m22 + sqrt(pow(pow(h1,2)*(3.*l1 + l3 + l4 + l5) + pow(h2,2)*(3.*l2 + l3 + l4 + l5) + 2.*m11 + 2.*m22,2) - 4.*(3.*pow(h1,4)*l1*l3 + 3.*pow(h2,4)*l2*l3 - 3.*pow(h1,2)*pow(h2,2)*pow(l3,2) + 3.*pow(h1,4)*l1*l4 + 3.*pow(h2,4)*l2*l4 - 6.*pow(h1,2)*pow(h2,2)*l3*l4 - 3.*pow(h1,2)*pow(h2,2)*pow(l4,2) + 3.*pow(h2,4)*l2*l5 - 3.*pow(h1,2)*pow(h2,2)*pow(l5,2) + 2.*pow(h1,2)*l3*m11 + 2.*pow(h1,2)*l4*m11 + pow(h1,2)*l5*(3.*pow(h1,2)*l1 - 6.*pow(h2,2)*(l3 + l4) + 2.*m11) + 8.*h1*h2*(l3 + l4 + l5)*pow(m12,2) - 4.*pow(m12,4) + 2.*pow(h2,2)*(l3 + l4 + l5)*m22 + (3.*pow(h1,2)*l1 + 2.*m11)*(3.*pow(h2,2)*l2 + 2.*m22))));

      std::complex<double> scalarNLO = - 1/(12.*M_PI) * ( 2.*pow(eigen_1,1.5) + 2.*pow(eigen_2,1.5) + pow(eigen_3,1.5) + pow(eigen_4,1.5) + pow(eigen_5,1.5) + pow(eigen_6,1.5));

      std::complex<double> counter_1 = 0.5*(m11 + m22 - 1.*sqrt(pow(m11,2) + 4.*pow(m12,4) - 2.*m11*m22 + pow(m22,2)));
      std::complex<double> counter_2 = 0.5*(m11 + m22 + sqrt(pow(m11,2) + 4.*pow(m12,4) - 2.*m11*m22 + pow(m22,2)));

      std::complex<double> counterTerm = 1/(12.*M_PI) * ( 4.*pow(counter_1, 1.5) + 4.*pow(counter_2, 1.5) );

      std::complex<double> wMassTerm = 0.25 * gWsq * (h1sq + h2sq);
      std::complex<double> zMassTerm = 0.25 * (gWsq + gYsq) * (h1sq + h2sq);

      std::complex<double> gaugeNLO = - 1/(6.*M_PI) * ( 2.*pow(wMassTerm,1.5) + pow(zMassTerm,1.5) );

      std::complex<double> veffNLO = scalarNLO + gaugeNLO + counterTerm;
      
      return T * veffLO.real() + T * veffNLO.real();
    
    }

    void Betas( const std::vector<double>& x, std::vector<double>& dxdt, const double t) override {
    
      double gYsq = x[0];
      double gWsq = x[1];
      double gSsq = x[2];
      double m11 = x[3];
      double m22 = x[4];
      double m12 = x[5];
      double l1 = x[6];
      double l2 = x[7];
      double l3 = x[8];
      double l4 = x[9];
      double l5 = x[10];
      double yt = x[11];
      dxdt[0] = 1./t * 0.08865603568704555*pow(gYsq,2);
      dxdt[1] = 1./t * -0.037995443865876666*pow(gWsq,2);
      dxdt[2] = 1./t * -0.08865603568704555*pow(gSsq,2);
      dxdt[3] = 1./t * 0.0031662869888230555*(4.*(2.*l3 + l4)*m22 - 3.*m11*(3.*gWsq + gYsq - 4.*(l1 + pow(yt,2))));
      dxdt[4] = 1./t * 0.0031662869888230555*(8.*l3*m11 + 4.*l4*m11 - 3.*(3.*gWsq + gYsq - 4.*l2)*m22);
      dxdt[5] = 1./t * 0.0031662869888230555*m12*(-9.*gWsq - 3.*gYsq + 4.*l3 + 8.*l4 + 12.*l5 + 6.*pow(yt,2));
      dxdt[6] = 1./t * 0.0015831434944115277*(9.*pow(gWsq,2) + 3.*pow(gYsq,2) + 6.*gWsq*(gYsq - 6.*l1) - 12.*gYsq*l1 + 8.*(6.*pow(l1,2) + 2.*pow(l3,2) + 2.*l3*l4 + pow(l4,2) + pow(l5,2) + 6.*l1*pow(yt,2) - 6.*pow(yt,4)));
      dxdt[7] = 1./t * 0.0015831434944115277*(9.*pow(gWsq,2) + 3.*pow(gYsq,2) + 6.*gWsq*(gYsq - 6.*l2) - 12.*gYsq*l2 + 8.*(6.*pow(l2,2) + 2.*pow(l3,2) + 2.*l3*l4 + pow(l4,2) + pow(l5,2)));
      dxdt[8] = 1./t * 0.0015831434944115277*(9.*pow(gWsq,2) + 3.*pow(gYsq,2) - 12.*gYsq*l3 - 6.*gWsq*(gYsq + 6.*l3) + 8.*(3.*l1*l3 + 3.*l2*l3 + 2.*pow(l3,2) + l1*l4 + l2*l4 + pow(l4,2) + pow(l5,2) + 3.*l3*pow(yt,2)));
      dxdt[9] = 1./t * 0.006332573977646111*(3.*gWsq*(gYsq - 3.*l4) - 3.*gYsq*l4 + 8.*pow(l5,2) + 2.*l4*(l1 + l2 + 4.*l3 + 2.*l4 + 3.*pow(yt,2)));
      dxdt[10] = 1./t * 0.006332573977646111*l5*(-9.*gWsq - 3.*gYsq + 2.*(l1 + l2 + 4.*l3 + 6.*l4 + 3.*pow(yt,2)));
      dxdt[11] = 1./t * 0.0005277144981371759*yt*(-96.*gSsq - 27.*gWsq - 17.*gYsq + 54.*pow(yt,2));
    }

    std::vector<double> get_4d_parameters(double T) const {

      double Gamma = scaleFactor*T;
      std::vector<double> params;

      for (int i = 0; i < RGEs.size(); i += 1 ){
        params.push_back(alglib::spline1dcalc(RGEs[i], Gamma));
      }

      return params;
    }

    std::vector<double> get_3d_parameters(double T) const {
    
      double Gamma = scaleFactor * T;
      double scaleFactor3;
      double scaleFactor3dUS;
      double Lb = 2. * log(scaleFactor) + 2. * EulerGamma - 2. * log(4. * M_PI);
      double Lf = Lb + 4. * log(2.);
    
      double gYsq = alglib::spline1dcalc(RGEs[0], Gamma);
      double gWsq = alglib::spline1dcalc(RGEs[1], Gamma);
      double gSsq = alglib::spline1dcalc(RGEs[2], Gamma);
      double m11 = alglib::spline1dcalc(RGEs[3], Gamma);
      double m22 = alglib::spline1dcalc(RGEs[4], Gamma);
      double m12 = alglib::spline1dcalc(RGEs[5], Gamma);
      double l1 = alglib::spline1dcalc(RGEs[6], Gamma);
      double l2 = alglib::spline1dcalc(RGEs[7], Gamma);
      double l3 = alglib::spline1dcalc(RGEs[8], Gamma);
      double l4 = alglib::spline1dcalc(RGEs[9], Gamma);
      double l5 = alglib::spline1dcalc(RGEs[10], Gamma);
      double yt = alglib::spline1dcalc(RGEs[11], Gamma);
      //---------------------------------------------------------------------------------
      double gYsq3d = gYsq*T - 0.0021108579925487037*pow(gYsq,2)*(Lb + 20.*Lf)*T;
      double gWsq3d = gWsq*T + 0.0021108579925487037*pow(gWsq,2)*(2. + 21.*Lb - 12.*Lf)*T;
      double gSsq3d = gSsq*T + 0.006332573977646111*pow(gSsq,2)*(1. + 11.*Lb - 4.*Lf)*T;
      double l13d = -0.0007915717472057639*T*(-12.*gYsq*l1*Lb + pow(gYsq,2)*(-2. + 3.*Lb) + pow(gWsq,2)*(-6. + 9.*Lb) + gWsq*(-36.*l1*Lb + gYsq*(-4. + 6.*Lb)) + 8.*(6.*pow(l1,2)*Lb + (2.*pow(l3,2) + 2.*l3*l4 + pow(l4,2) + pow(l5,2))*Lb - 6.*Lf*pow(yt,4) + l1*(-157.91367041742973 + 6.*Lf*pow(yt,2))));
      double l23d = -0.0007915717472057639*(-1263.3093633394378*l2 - 12.*gYsq*l2*Lb + 8.*(6.*pow(l2,2) + 2.*pow(l3,2) + 2.*l3*l4 + pow(l4,2) + pow(l5,2))*Lb + pow(gYsq,2)*(-2. + 3.*Lb) + pow(gWsq,2)*(-6. + 9.*Lb) + gWsq*(-36.*l2*Lb + gYsq*(-4. + 6.*Lb)))*T;
      double l33d = -0.0007915717472057639*T*(-12.*gYsq*l3*Lb + pow(gYsq,2)*(-2. + 3.*Lb) + pow(gWsq,2)*(-6. + 9.*Lb) + gWsq*(gYsq*(4. - 6.*Lb) - 36.*l3*Lb) + 8.*(-157.91367041742973*l3 + 2.*pow(l3,2)*Lb + pow(l4,2)*Lb + l1*(3.*l3 + l4)*Lb + l2*(3.*l3 + l4)*Lb + pow(l5,2)*Lb + 3.*l3*Lf*pow(yt,2)));
      double l43d = -0.0031662869888230555*T*(-315.82734083485946*l4 - 3.*gYsq*l4*Lb + 2.*l1*l4*Lb + 2.*l2*l4*Lb + 8.*l3*l4*Lb + 4.*pow(l4,2)*Lb + 8.*pow(l5,2)*Lb + gWsq*(-9.*l4*Lb + gYsq*(-2. + 3.*Lb)) + 6.*l4*Lf*pow(yt,2));
      double l53d = 0.0031662869888230555*l5*T*(315.82734083485946 + 9.*gWsq*Lb + 3.*gYsq*Lb - 2.*(l1 + l2 + 4.*l3 + 6.*l4)*Lb - 6.*Lf*pow(yt,2));
      //---------------------------------------------------------------------------------
      double lambdaVLL_1 = 0.025330295910584444*pow(gSsq,2)*T;
      double lambdaVLL_2 = -0.07599088773175333*gSsq*gWsq*T;
      double lambdaVLL_3 = 0.025330295910584444*pow(gWsq,2)*T;
      double lambdaVLL_4 = -0.025330295910584444*pow(gSsq,1.5)*sqrt(gYsq)*T;
      double lambdaVLL_5 = -0.09287775167214296*gSsq*gYsq*T;
      double lambdaVLL_6 = -0.025330295910584444*gWsq*gYsq*T;
      double lambdaVLL_7 = -0.5094203955350871*pow(gYsq,2)*T;
      double lambdaVL_1 = 0.;
      double lambdaVL_2 = 0.;
      double lambdaVL_3 = 0.;
      double lambdaVL_4 = -0.025330295910584444*gSsq*T*pow(yt,2);
      double lambdaVL_5 = -0.0005277144981371759*gYsq*T*(-9.*gWsq + gYsq*(-39. + 2.*Lb + 40.*Lf) - 4.*(236.8705056261446 + 9.*l1 + 6.*l3 + 3.*l4 - 17.*pow(yt,2)));
      double lambdaVL_6 = 0.0005277144981371759*sqrt(gWsq)*sqrt(gYsq)*T*(gWsq*(5. + 21.*Lb - 12.*Lf) - 1.*gYsq*(-21. + Lb + 20.*Lf) + 12.*(78.95683520871486 + l1 + l4 + pow(yt,2)));
      double lambdaVL_7 = 0.0005277144981371759*gWsq*T*(gWsq*(73. + 42.*Lb - 24.*Lf) + 3.*(gYsq + 4.*(78.95683520871486 + 3.*l1 + 2.*l3 + l4 - 3.*pow(yt,2))));
      double lambdaVL_8 = 0.;
      double lambdaVL_9 = 0.;
      double lambdaVL_10 = 0.;
      double lambdaVL_11 = 0.;
      double lambdaVL_12 = 0.;
      double lambdaVL_13 = -0.0005277144981371759*gYsq*(-9.*gWsq - 4.*(236.8705056261446 + 9.*l2 + 6.*l3 + 3.*l4) + gYsq*(-39. + 2.*Lb + 40.*Lf))*T;
      double lambdaVL_14 = 0.0005277144981371759*sqrt(gWsq)*sqrt(gYsq)*(12.*(78.95683520871486 + l2 + l4) + gWsq*(5. + 21.*Lb - 12.*Lf) - 1.*gYsq*(-21. + Lb + 20.*Lf))*T;
      double lambdaVL_15 = 0.0005277144981371759*gWsq*(3.*(gYsq + 4.*(78.95683520871486 + 3.*l2 + 2.*l3 + l4)) + gWsq*(73. + 42.*Lb - 24.*Lf))*T;
      //---------------------------------------------------------------------------------
      double MusqSU3_LO = 2.*gSsq*pow(T,2);
      double MusqSU2_LO = 2.*gWsq*pow(T,2);
      double MusqU1_LO = 2.*gYsq*pow(T,2);
      double MusqSU3_NLO = -0.02746081017515274*pow(gSsq,2)*pow(T,2) - 0.01424829144970375*gSsq*gWsq*pow(T,2) - 0.005804859479508935*gSsq*gYsq*pow(T,2) - 0.006332573977646111*gSsq*pow(T,2)*pow(yt,2) + 0.13931662750821444*pow(gSsq,2)*pow(T,2)*log(0.07957747154594767*scaleFactor) - 0.037995443865876666*pow(gSsq,2)*pow(T,2)*log(3.141592653589793*T) + 0.037995443865876666*pow(gSsq,2)*pow(T,2)*log(scaleFactor*T);
      double MusqSU2_NLO = 0.012665147955292222*gWsq*m11 + 0.012665147955292222*gWsq*m22 - 0.037995443865876666*gSsq*gWsq*pow(T,2) - 0.22824546686883768*pow(gWsq,2)*pow(T,2) - 0.0015831434944115277*gWsq*gYsq*pow(T,2) + 0.0031662869888230555*gWsq*l1*pow(T,2) + 0.0031662869888230555*gWsq*l2*pow(T,2) + 0.004221715985097407*gWsq*l3*pow(T,2) + 0.0021108579925487037*gWsq*l4*pow(T,2) - 0.0015831434944115277*gWsq*pow(T,2)*pow(yt,2) + 0.04080992118927494*pow(gWsq,2)*pow(T,2)*log(scaleFactor) - 0.03518096654247839*pow(gWsq,2)*pow(T,2)*log(T) + 0.03518096654247839*pow(gWsq,2)*pow(T,2)*log(scaleFactor*T);
      double MusqU1_NLO = 0.012665147955292222*gYsq*m11 + 0.012665147955292222*gYsq*m22 - 0.04643887583607148*gSsq*gYsq*pow(T,2) - 0.004749430483234583*gWsq*gYsq*pow(T,2) + 0.15648431719011083*pow(gYsq,2)*pow(T,2) + 0.0031662869888230555*gYsq*l1*pow(T,2) + 0.0031662869888230555*gYsq*l2*pow(T,2) + 0.004221715985097407*gYsq*l3*pow(T,2) + 0.0021108579925487037*gYsq*l4*pow(T,2) - 0.005804859479508935*gYsq*pow(T,2)*pow(yt,2) - 0.03518096654247839*pow(gYsq,2)*pow(T,2)*log(scaleFactor) + 0.1421311048316127*pow(gYsq,2)*pow(T,2)*log(T) - 0.1421311048316127*pow(gYsq,2)*pow(T,2)*log(scaleFactor*T);
    
      double MusqSU3 = MusqSU3_LO + MusqSU3_NLO;
      double MusqSU2 = MusqSU2_LO + MusqSU2_NLO;
      double MusqU1 = MusqU1_LO + MusqU1_NLO;
      //---------------------------------------------------------------------------------
      scaleFactor3 = sqrt(MusqU1);
      scaleFactor3dUS = sqrt(gYsq)*scaleFactor3;
      //---------------------------------------------------------------------------------
      double m113d_LO = m11 + 0.020833333333333332*pow(T,2)*(9.*gWsq + 3.*gYsq + 4.*(3.*l1 + 2.*l3 + l4 + 3.*pow(yt,2)));
      double m223d_LO = m22 + 0.020833333333333332*(9.*gWsq + 3.*gYsq + 12.*l2 + 8.*l3 + 4.*l4)*pow(T,2);
      double m123d_LO = m12;
      double l5I3d = 0, l6R3d = 0, l6I3d = 0, l7R3d = 0, l7I3d = 0;
      double m113d_NLO = -0.000021988104089048995*(-216.*gYsq*Lb*m11 + 864.*l1*Lb*m11 + 576.*l3*Lb*m22 + 288.*l4*Lb*m22 - 246.03488281981413*pow(gYsq,2)*pow(T,2) + 484.09302085283775*gYsq*l1*pow(T,2) - 1040.1860417056755*pow(l1,2)*pow(T,2) + 322.72868056855845*gYsq*l3*pow(T,2) - 693.4573611371169*pow(l3,2)*pow(T,2) + 161.36434028427922*gYsq*l4*pow(T,2) - 693.4573611371169*l3*l4*pow(T,2) - 693.4573611371169*pow(l4,2)*pow(T,2) - 1040.1860417056755*pow(l5,2)*pow(T,2) + 150.*pow(gYsq,2)*Lb*pow(T,2) + 108.*gYsq*l1*Lb*pow(T,2) + 72.*gYsq*l3*Lb*pow(T,2) + 144.*l1*l3*Lb*pow(T,2) + 144.*l2*l3*Lb*pow(T,2) - 48.*pow(l3,2)*Lb*pow(T,2) + 36.*gYsq*l4*Lb*pow(T,2) + 72.*l1*l4*Lb*pow(T,2) + 72.*l2*l4*Lb*pow(T,2) - 48.*l3*l4*Lb*pow(T,2) - 120.*pow(l4,2)*Lb*pow(T,2) - 216.*pow(l5,2)*Lb*pow(T,2) - 60.*pow(gYsq,2)*Lf*pow(T,2) - 9.*pow(gWsq,2)*(-113.58785446279089 - 84.*Lb + 12.*Lf)*pow(T,2) + 864.*Lf*m11*pow(yt,2) + 576.*gSsq*pow(T,2)*pow(yt,2) + 66.*gYsq*pow(T,2)*pow(yt,2) + 192.*gSsq*Lb*pow(T,2)*pow(yt,2) - 47.*gYsq*Lb*pow(T,2)*pow(yt,2) + 324.*l1*Lb*pow(T,2)*pow(yt,2) - 768.*gSsq*Lf*pow(T,2)*pow(yt,2) - 55.*gYsq*Lf*pow(T,2)*pow(yt,2) + 108.*l1*Lf*pow(T,2)*pow(yt,2) + 144.*l3*Lf*pow(T,2)*pow(yt,2) + 72.*l4*Lf*pow(T,2)*pow(yt,2) - 108.*Lb*pow(T,2)*pow(yt,4) - 9.*gWsq*(pow(T,2)*(66.23514178511634*gYsq - 161.3643402842792*l1 - 107.57622685618614*l3 - 53.78811342809307*l4 - 6.*pow(yt,2) - 3.*Lf*pow(yt,2)) + 3.*Lb*(24.*m11 + pow(T,2)*(8.*gYsq - 12.*l1 - 8.*l3 - 4.*l4 + 7.*pow(yt,2)))) + 18.*log(scaleFactor3/scaleFactor*T)*(33.*pow(gWsq3d,2) - 7.*pow(gYsq3d,2) + 8.*gYsq3d*(3.*l13d + 2.*l33d + l43d) + 6.*gWsq3d*(-3.*gYsq3d + 4.*(3.*l13d + 2.*l33d + l43d + 4.*lambdaVL_7)) - 8.*(6.*pow(l13d,2) + 4.*pow(l33d,2) + 4.*l33d*l43d + 4.*pow(l43d,2) + 6.*pow(l5I3d,2) + 6.*pow(l53d,2) + 18.*pow(l6I3d,2) + 18.*pow(l6R3d,2) + 6.*pow(l7I3d,2) + 6.*pow(l7R3d,2) + 3.*pow(lambdaVL_1,2) + 6.*pow(lambdaVL_2,2) + pow(lambdaVL_3,2) - 48.*gSsq3d*lambdaVL_4 + 8.*pow(lambdaVL_4,2) + pow(lambdaVL_5,2) + 6.*pow(lambdaVL_6,2) + 3.*pow(lambdaVL_7,2) + 8.*pow(lambdaVL_8,2) + pow(lambdaVL_10,2) + 3.*pow(lambdaVL_11,2) + 6.*pow(lambdaVL_12,2))));
      double m223d_NLO = 0.00006596431226714699*(-96.*l4*Lb*m11 + 216.*gWsq*Lb*m22 + 72.*gYsq*Lb*m22 - 288.*l2*Lb*m22 - 340.76356338837263*pow(gWsq,2)*pow(T,2) + 198.70542535534904*gWsq*gYsq*pow(T,2) + 82.01162760660472*pow(gYsq,2)*pow(T,2) - 484.09302085283775*gWsq*l2*pow(T,2) - 161.36434028427922*gYsq*l2*pow(T,2) + 346.72868056855845*pow(l2,2)*pow(T,2) - 161.36434028427922*gWsq*l4*pow(T,2) - 53.78811342809307*gYsq*l4*pow(T,2) + 231.15245371237228*pow(l4,2)*pow(T,2) + 346.72868056855845*pow(l5,2)*pow(T,2) - 16.*pow(l3,2)*(-14.447028357023267 - 1.*Lb)*pow(T,2) - 252.*pow(gWsq,2)*Lb*pow(T,2) + 72.*gWsq*gYsq*Lb*pow(T,2) - 50.*pow(gYsq,2)*Lb*pow(T,2) - 108.*gWsq*l2*Lb*pow(T,2) - 36.*gYsq*l2*Lb*pow(T,2) - 36.*gWsq*l4*Lb*pow(T,2) - 12.*gYsq*l4*Lb*pow(T,2) - 24.*l1*l4*Lb*pow(T,2) - 24.*l2*l4*Lb*pow(T,2) + 40.*pow(l4,2)*Lb*pow(T,2) + 72.*pow(l5,2)*Lb*pow(T,2) + 36.*pow(gWsq,2)*Lf*pow(T,2) + 20.*pow(gYsq,2)*Lf*pow(T,2) - 36.*l4*Lb*pow(T,2)*pow(yt,2) + 12.*l4*Lf*pow(T,2)*pow(yt,2) - 8.*l3*Lb*(24.*m11 + pow(T,2)*(9.*gWsq + 3.*gYsq + 6.*l1 + 6.*l2 - 2.*l4 + 9.*pow(yt,2))) + 8.*l3*pow(T,2)*(-40.3410850710698*gWsq - 13.447028357023267*gYsq + 3.*(9.631352238015513*l4 + Lf*pow(yt,2))) - 6.*log(scaleFactor3/scaleFactor*T)*(33.*pow(gWsq3d,2) - 7.*pow(gYsq3d,2) + 8.*gYsq3d*(3.*l23d + 2.*l33d + l43d) - 8.*(6.*pow(l23d,2) + 4.*pow(l33d,2) + 4.*l33d*l43d + 4.*pow(l43d,2) + 6.*pow(l5I3d,2) + 6.*pow(l53d,2) + 6.*pow(l6I3d,2) + 6.*pow(l6R3d,2) + 18.*pow(l7I3d,2) + 18.*pow(l7R3d,2) + 3.*pow(lambdaVL_1,2) + 6.*pow(lambdaVL_2,2) + pow(lambdaVL_3,2) + 8.*pow(lambdaVL_8,2) - 48.*gSsq3d*lambdaVL_9 + 8.*pow(lambdaVL_9,2) + pow(lambdaVL_10,2) + 3.*pow(lambdaVL_11,2) + 6.*pow(lambdaVL_12,2) + pow(lambdaVL_13,2) + 6.*pow(lambdaVL_14,2) + 3.*pow(lambdaVL_15,2)) + 6.*gWsq3d*(-3.*gYsq3d + 4.*(3.*l23d + 2.*l33d + l43d + 4.*lambdaVL_15))));
      double m123d_NLO = 0.0015831434944115277*(m12*(9.*gWsq*Lb + 3.*gYsq*Lb - 2.*(2.*l3*Lb + 4.*l4*Lb + 6.*l5*Lb + 3.*Lf*pow(yt,2))) - 2.*log(scaleFactor3/scaleFactor*T)*(-3.*gYsq3d*l6R3d + 6.*l13d*l6R3d + 6.*l33d*l6R3d + 6.*l43d*l6R3d + 6.*l53d*l6R3d + 6.*l5I3d*(l6I3d + l7I3d) - 3.*gYsq3d*l7R3d + 6.*l23d*l7R3d + 6.*l33d*l7R3d + 6.*l43d*l7R3d + 6.*l53d*l7R3d - 48.*gSsq3d*lambdaVL_8 + 8.*lambdaVL_4*lambdaVL_8 + 8.*lambdaVL_8*lambdaVL_9 + lambdaVL_5*lambdaVL_10 + 3.*lambdaVL_7*lambdaVL_11 - 3.*gWsq3d*(3.*l6R3d + 3.*l7R3d + 4.*lambdaVL_11) - 6.*lambdaVL_6*lambdaVL_12 + lambdaVL_10*lambdaVL_13 - 6.*lambdaVL_12*lambdaVL_14 + 3.*lambdaVL_11*lambdaVL_15));
    
      double m113d = m113d_LO + m113d_NLO;
      double m223d = m223d_LO + m223d_NLO;
      double m123d = m123d_LO + m123d_NLO;
      
      return {gYsq, gWsq, gSsq, m113d, m223d, m123d, l13d, l23d, l33d, l43d, l53d};
    }

  private :
  
    double gYsq_input;
    double gWsq_input;
    double gSsq_input;
    double m11_input;
    double m22_input;
    double m12_input;
    double l1_input;
    double l2_input;
    double l3_input;
    double l4_input;
    double l5_input;
    double yt_input;
  
    double scaleFactor = M_PI;

};

} // namespace EffectivePotential

#endif