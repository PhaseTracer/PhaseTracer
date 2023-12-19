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

#ifndef POTENTIAL_DRALGO_AH_MODEL_HPP_INCLUDED
#define POTENTIAL_DRALGO_AH_MODEL_HPP_INCLUDED

/**
   The Abelian Higgs Model in  DRalgo
   https://arxiv.org/pdf/2205.08815.pdf

   /. {\[Phi]^2 -> phisq, \[Phi]^4 -> phisq^2, \[Lambda] -> lam ,
   g1^2 -> g1sq, g1^4 -> g1sq^2, Y\[Phi]^2 -> Yphisq,
   Y\[Phi]^4 -> Yphisq^2, \[Mu]3US -> mu3US}
 
*/

#include <vector>
#include <cmath>
#include <interpolation.h>

#include "pow.hpp"
#include "potential.hpp"
#include "pow.hpp"

namespace EffectivePotential {

class DR_ah: public Potential {
 public:

  DR_ah(double M_input,
        double gsq_input,
        double lam_input){
        double muhsq_input = -0.5*square(M_input);
      solveBetas({gsq_input,muhsq_input,lam_input}); // The order is same to `Betas'
    }

  size_t get_n_scalars() const override {return 1;}

  std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
    auto phi1 = phi;
    phi1[0] = - phi[0];
    return {phi1};
  };

  double V(Eigen::VectorXd phi, double T) const override {

    const std::vector<double> par = DRstep(T);
//    std::cout << "par = " << par[0] << ", " <<  par[1] << ", " << par[2] << std::endl;
    std::complex<double> msq(par[1],0);
    std::complex<double> lam(par[2],0);
    std::complex<double> g1sq(par[0],0);

    std::complex<double> phisq(phi[0] * phi[0] / T,0);
    
    auto mu3US = g1sq;
    
    std::complex<double> veffLO = (msq*phisq)/2. + (lam*pow(phisq,2))/4.;
    
    std::complex<double> veffNLO = -pow(msq + lam*phisq,1.5)/(12.*M_PI) -
                      pow(msq + 3.*lam*phisq,1.5)/(12.*M_PI) -
                      pow(g1sq*phisq*Yphisq,1.5)/(6.*M_PI);

    std::complex<double> veffNNLO = (3.*lam*(msq + lam*phisq))/(64.*pow(M_PI,2)) +
            (lam*sqrt(msq + lam*phisq)*sqrt(msq + 3.*lam*phisq))/(32.*pow(M_PI,2)) +
            (3.*lam*(msq + 3.*lam*phisq))/(64.*pow(M_PI,2)) +
            (g1sq*sqrt(msq + lam*phisq)*Yphisq*sqrt(g1sq*phisq*Yphisq))/
            (16.*pow(M_PI,2)) + (g1sq*sqrt(msq + 3.*lam*phisq)*Yphisq*
            sqrt(g1sq*phisq*Yphisq))/(16.*pow(M_PI,2)) -
            (3.*pow(lam,2)*phisq*(0.5 + log(mu3US/(3.*sqrt(msq + 3.*lam*phisq)))))/
            (16.*pow(M_PI,2)) - (pow(lam,2)*phisq*
            (0.5 + log(mu3US/(2.*sqrt(msq + lam*phisq) + sqrt(msq + 3.*lam*phisq)))))/
            (16.*pow(M_PI,2)) + ((g1sq*phisq*sqrt(msq + lam*phisq)*
            sqrt(msq + 3.*lam*phisq)*Yphisq)/(16.*pow(M_PI,2)) +
            (sqrt(msq + 3.*lam*phisq)*sqrt(g1sq*phisq*Yphisq)*
            (-2.*lam*phisq - g1sq*phisq*Yphisq))/(16.*pow(M_PI,2)) +
            (sqrt(msq + lam*phisq)*sqrt(g1sq*phisq*Yphisq)*
            (2.*lam*phisq - g1sq*phisq*Yphisq))/(16.*pow(M_PI,2)) +
            (pow(lam,2)*pow(phisq,2)*
            (0.5 + log(mu3US/(sqrt(msq + lam*phisq) + sqrt(msq + 3.*lam*phisq)))))/
            (4.*pow(M_PI,2)) - ((pow(msq + lam*phisq,2) -
            2.*(msq + lam*phisq)*(msq + 3.*lam*phisq) + pow(msq + 3.*lam*phisq,2) -
            2.*g1sq*phisq*(msq + lam*phisq)*Yphisq -
            2.*g1sq*phisq*(msq + 3.*lam*phisq)*Yphisq +
            pow(g1sq,2)*pow(phisq,2)*pow(Yphisq,2))*
            (0.5 + log(mu3US/
            (sqrt(msq + lam*phisq) + sqrt(msq + 3.*lam*phisq) +
            sqrt(g1sq*phisq*Yphisq)))))/(16.*pow(M_PI,2)))/(2.*phisq) +
            pow(g1sq,2)*phisq*pow(Yphisq,2)*
            (1/(32.*pow(M_PI,2)) - (0.5 +
            log(mu3US/(sqrt(msq + 3.*lam*phisq) + 2.*sqrt(g1sq*phisq*Yphisq))))/
            (16.*pow(M_PI,2)) + ((g1sq*phisq*sqrt(msq + 3.*lam*phisq)*Yphisq*
            sqrt(g1sq*phisq*Yphisq))/(8.*pow(M_PI,2)) +
            (g1sq*phisq*Yphisq*(msq + 3.*lam*phisq - 2.*g1sq*phisq*Yphisq))/
            (16.*pow(M_PI,2)) - (pow(msq + 3.*lam*phisq,2)*
            (0.5 + log(mu3US/sqrt(msq + 3.*lam*phisq))))/(16.*pow(M_PI,2)) +
            (pow(-msq - 3.*lam*phisq + g1sq*phisq*Yphisq,2)*
            (0.5 + log(mu3US/(sqrt(msq + 3.*lam*phisq) + sqrt(g1sq*phisq*Yphisq))))
            )/(8.*pow(M_PI,2)) -
            (pow(-msq - 3.*lam*phisq + 2.*g1sq*phisq*Yphisq,2)*
            (0.5 + log(mu3US/
            (sqrt(msq + 3.*lam*phisq) + 2.*sqrt(g1sq*phisq*Yphisq)))))/
            (16.*pow(M_PI,2)))/(4.*pow(g1sq,2)*pow(Yphisq,2)*pow(phisq,2)));

//    std::cout << std::endl << "veffLO=" << veffLO << "veffNLO=" << veffNLO << "veffNNLO=" << veffNNLO << std::endl;
    return veffLO.real() + veffNLO.real() + veffNNLO.real();
  }
  
  void Betas(const std::vector<double>& x, std::vector<double>& dxdt, const double t) override{
    double gsq   = x[0];
    double muhsq = x[1];
    double lam   = x[2];
    dxdt[0] = 1/t * ( Yphisq * square(gsq)) / (24.*square(M_PI));
    dxdt[1] = 1/t * ( -(3.*Yphisq * gsq * muhsq) / (8.*square(M_PI))
                      + (lam * muhsq) / (2.*square(M_PI)) );
    dxdt[2] = 1/t * ((3.*square(Yphisq) * square(gsq))/ (8.*square(M_PI))
                      -(3.*Yphisq * gsq*lam )/ (4.*square(M_PI))
                      +(5.*square(lam) )/ (4.*square(M_PI)) );
  }

  std::vector<double> DRstep(double T) const {
    double scale = scaleFactor * M_PI * T;
    
    double gsq   = alglib::spline1dcalc(RGEs[0], scale);
    double muhsq = alglib::spline1dcalc(RGEs[1], scale);
    double lamh  = alglib::spline1dcalc(RGEs[2], scale);
    
    double Lb = -2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T);
    
    double bargsq3 = gsq*T - (pow(gsq,2)*T*Yphisq*(-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)))/
    (48.*pow(M_PI,2));
    
    double barlamh3 = (T*(6.*gsq*lamh*Yphisq*(-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)) +
                      pow(gsq,2)*pow(Yphisq,2)*
                      (2 - 3.*(-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T))) +
                      2.*lamh*(8.*pow(M_PI,2) - 5.*lamh*
                      (-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)))))/(16.*pow(M_PI,2))\
                      - (pow(gsq,2)*pow(T,2)*pow(Yphisq,2)*
                      pow(24.*lamh + 48.*pow(M_PI,2) -
                      gsq*Yphisq*(-4 - 2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)),2))/
                      (18432.*pow(M_PI,5)*sqrt((gsq*pow(T,2)*Yphisq)/3. +
                      (gsq*Yphisq*(72.*muhsq + pow(T,2)*
                      (24.*lamh - 2.*gsq*Yphisq*
                      (-7 - 2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)))))/
                      (288.*pow(M_PI,2))));
    double mu3US=bargsq3;
    
    double barmuhsq3 = muhsq + (pow(T,2)*(4.*lamh + 3.*gsq*Yphisq))/12. +
                      (gsq*T*Yphisq*(24.*lamh + 48.*pow(M_PI,2) -
                      gsq*Yphisq*(-4 - 2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)))*
                      ((pow(gsq,2)*T*pow(Yphisq,2))/pow(M_PI,2) -
                      (gsq*T*Yphisq*(24.*lamh + 48.*pow(M_PI,2) -
                      gsq*Yphisq*(-4 - 2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T))))/
                      (12.*pow(M_PI,2)) + (gsq*T*Yphisq*log(2)*
                      (24.*lamh + 48.*pow(M_PI,2) -
                      gsq*Yphisq*(-4 - 2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T))))/
                      (6.*pow(M_PI,2))))/(3072.*pow(M_PI,4)) -
                      (gsq*T*Yphisq*(24.*lamh + 48.*pow(M_PI,2) -
                      gsq*Yphisq*(-4 - 2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)))*
                      sqrt((gsq*pow(T,2)*Yphisq)/3. +
                      (gsq*Yphisq*(72.*muhsq + pow(T,2)*
                      (24.*lamh - 2.*gsq*Yphisq*
                      (-7 - 2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)))))/
                      (288.*pow(M_PI,2))))/(192.*pow(M_PI,3)) +
                      ((pow(Yphisq,2)*(gsq*T - (pow(gsq,2)*T*Yphisq*
                      (-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)))/(48.*pow(M_PI,2)))\
                      - 2.*Yphisq*(gsq*T - (pow(gsq,2)*T*Yphisq*
                      (-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)))/(48.*pow(M_PI,2)))*
                      ((T*(6.*gsq*lamh*Yphisq*(-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)) +
                      pow(gsq,2)*pow(Yphisq,2)*
                      (2 - 3.*(-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T))) +
                      2.*lamh*(8.*pow(M_PI,2) -
                      5.*lamh*(-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)))))/
                      (16.*pow(M_PI,2)) - (pow(gsq,2)*pow(T,2)*pow(Yphisq,2)*
                      pow(24.*lamh + 48.*pow(M_PI,2) -
                      gsq*Yphisq*(-4 - 2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)),2))
                      /(18432.*pow(M_PI,5)*
                      sqrt((gsq*pow(T,2)*Yphisq)/3. +
                      (gsq*Yphisq*(72.*muhsq +
                      pow(T,2)*(24.*lamh -
                      2.*gsq*Yphisq*
                      (-7 - 2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)))))/
                      (288.*pow(M_PI,2))))) +
                      2.*pow((T*(6.*gsq*lamh*Yphisq*
                      (-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)) +
                      pow(gsq,2)*pow(Yphisq,2)*
                      (2 - 3.*(-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T))) +
                      2.*lamh*(8.*pow(M_PI,2) -
                      5.*lamh*(-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)))))/
                      (16.*pow(M_PI,2)) - (pow(gsq,2)*pow(T,2)*pow(Yphisq,2)*
                      pow(24.*lamh + 48.*pow(M_PI,2) -
                      gsq*Yphisq*(-4 - 2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)),2))
                      /(18432.*pow(M_PI,5)*
                      sqrt((gsq*pow(T,2)*Yphisq)/3. +
                      (gsq*Yphisq*(72.*muhsq +
                      pow(T,2)*(24.*lamh -
                      2.*gsq*Yphisq*
                      (-7 - 2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)))))/
                      (288.*pow(M_PI,2)))),2))*
                      log(mu3US/sqrt((gsq*pow(T,2)*Yphisq)/3. +
                      (gsq*Yphisq*(72.*muhsq +
                      pow(T,2)*(24.*lamh -
                      2.*gsq*Yphisq*(-7 - 2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T))
                      )))/(288.*pow(M_PI,2)))))/(4.*pow(M_PI,2)) +
                      (pow(gsq,2)*pow(T,2)*pow(Yphisq,2)*
                      (-8 - 108.*EulerGamma + 1296.*log(Glaisher) +
                      69.*(-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T))) +
                      12.*gsq*Yphisq*(2.*lamh*pow(T,2)*(1 + 6.*EulerGamma - 72.*log(Glaisher)) +
                      (9.*muhsq - 6.*lamh*pow(T,2))*
                      (-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T))) +
                      24.*lamh*(-6.*lamh*pow(T,2)*(EulerGamma - 12.*log(Glaisher)) +
                      (-6.*muhsq + lamh*pow(T,2))*
                      (-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T))) +
                      18.*((pow(gsq,2)*pow(T,2)*pow(Yphisq,2)*
                      pow(24.*lamh + 48.*pow(M_PI,2) -
                      gsq*Yphisq*(-4 - 2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)),2))/
                      (576.*pow(M_PI,4)) + 8.*pow(Yphisq,2)*
                      pow(gsq*T - (pow(gsq,2)*T*Yphisq*
                      (-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)))/(48.*pow(M_PI,2)),
                      2) - (T*Yphisq*(gsq*T -
                      (pow(gsq,2)*T*Yphisq*
                      (-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)))/
                      (48.*pow(M_PI,2)))*
                      (6.*gsq*lamh*Yphisq*(-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)) +
                      pow(gsq,2)*pow(Yphisq,2)*
                      (2 - 3.*(-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T))) +
                      2.*lamh*(8.*pow(M_PI,2) -
                      5.*lamh*(-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)))))/
                      pow(M_PI,2) + (pow(T,2)*
                      pow(6.*gsq*lamh*Yphisq*
                      (-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T)) +
                      pow(gsq,2)*pow(Yphisq,2)*
                      (2 - 3.*(-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T))) +
                      2.*lamh*(8.*pow(M_PI,2) -
                      5.*lamh*(-2.*(-EulerGamma + log(4.*M_PI)) + 2.*log(scale/T))),2))/
                      (16.*pow(M_PI,4)))*log(sqrt((gsq*pow(T,2)*Yphisq)/3. +
                      (gsq*Yphisq*(72.*muhsq +
                      pow(T,2)*(24.*lamh -
                      2.*gsq*Yphisq*(-7 - 2.*(-EulerGamma + log(4.*M_PI)) +
                      2.*log(scale/T)))))/(288.*pow(M_PI,2)))/scale))/
                      (576.*pow(M_PI,2));
    
    return {bargsq3, barmuhsq3, barlamh3};
  };
  

 private:
  
  static constexpr double Yphisq = 1;
  const double scaleFactor = 1;

};

}  // namespace EffectivePotential

#endif
