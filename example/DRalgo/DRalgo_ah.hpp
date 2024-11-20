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

#ifndef POTENTIAL_DRALGO_ah_HPP_INCLUDED
#define POTENTIAL_DRALGO_ah_HPP_INCLUDED

#include <vector>
#include <cmath>
#include <complex>
#include <Eigen/Eigenvalues>
#include <interpolation.h>

#include "pow.hpp"
#include "potential.hpp"

namespace EffectivePotential {

class DR_ah: public Potential {
  public :

    DR_ah(double gsq_input_, double M_input_, double lam_input_) {
    
      gsq_input = gsq_input_;
      m_input = -0.5*square(M_input_);
      lam_input = lam_input_;
    
      std::vector<double> x0 = {gsq_input, m_input, lam_input};
    
      solveBetas(x0);
    
    }

    size_t get_n_scalars() const override { return 1;}

    std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
      auto phi1 = phi;
      phi1[0] = - phi[0];
      return {phi1};
    };

    void debug_print() {

    std::cout << "Input parameters:\n";
    std::cout << "\tgsq(Q=100) = " << gsq_input << "\n";
    std::cout << "\tm(Q=100) = " << m_input << "\n";
    std::cout << "\tlam(Q=100) = " << lam_input << "\n";

    std::cout << "Q = 20 parameters:\n";
    std::cout << "\tgsq(Q=20) = " << alglib::spline1dcalc(RGEs[0], 20) << "\n";
    std::cout << "\tm(Q=200) = " << alglib::spline1dcalc(RGEs[1], 20) << "\n";
    std::cout << "\tlam(Q=200) = " << alglib::spline1dcalc(RGEs[2], 20) << "\n";

    std::cout << "Q = 5000 parameters:\n";
    std::cout << "\tgsq(Q=5000) = " << alglib::spline1dcalc(RGEs[0], 5000) << "\n";
    std::cout << "\tm(Q=5000) = " << alglib::spline1dcalc(RGEs[1], 5000) << "\n";
    std::cout << "\tlam(Q=5000) = " << alglib::spline1dcalc(RGEs[2], 5000) << "\n";

    std::vector<double> par_1000 = get_3d_parameters(1000);
    std::complex<double> gsq3d_1000(par_1000[0],0);
    std::complex<double> msq3d_1000(par_1000[1],0);
    std::complex<double> lam3d_1000(par_1000[2],0);
    std::cout << "T = 1000 3d parameters:\n";
    std::cout << "\tgsq3d(T=1000) = " << gsq3d_1000.real() << "\n";
    std::cout << "\tm3d(T=1000) = " << msq3d_1000.real() << "\n";
    std::cout << "\tlam3d(T=1000) = " << lam3d_1000.real() << "\n";

    std::vector<double> par_100 = get_3d_parameters(100);
    std::complex<double> gsq3d_100(par_100[0],0);
    std::complex<double> msq3d_100(par_100[1],0);
    std::complex<double> lam3d_100(par_100[2],0);
    std::cout << "T = 100 3d parameters:\n";
    std::cout << "\tgsq3d(T=100) = " << gsq3d_100.real() << "\n";
    std::cout << "\tm3d(T=100) = " << msq3d_100.real() << "\n";
    std::cout << "\tlam3d(T=100) = " << lam3d_100.real() << "\n";

    }

    void plot_data() {

      std::string path = "my_stuff/PT2_plots/sec_7_3/data/T_";

      for ( double T = 50; T < 360; T += 10 ) {

        std::ofstream output_file;
        output_file.open(path + std::to_string(T) + ".txt"); // replace this path

        for ( double h = -30; h < 301; h += 1 ){
          Eigen::VectorXd phi(1);
          phi[0] = h;

          output_file << h << "\t" << V(phi, T) << std::endl;
        }

        output_file.close();
      }
    }

    double V(Eigen::VectorXd phi, double T) const override {
    
      const std::vector<double> par = get_3d_parameters(T);
    
      std::complex<double> gsq(par[0],0);
      std::complex<double> m(par[1],0);
      std::complex<double> lam(par[2],0);

      std::complex<double> h(phi[0]/sqrt(T + 1e-15),0);

      std::complex<double> Yphisq(1., 0);
      std::complex<double> phisq = pow(h, 2.);
      auto mu3US = gsq;

      std::complex<double> veffLO = 0.25*pow(h,4.)*lam + 0.5*pow(h,2.)*m;

      if ( potential_flag == 0 ) {
        return T * veffLO.real();
      }

      std::complex<double> veffNLO = - pow( 3.*pow(h,2.)*lam + m, 1.5)/(12. * M_PI) - pow( pow(h,2.)*lam + m, 1.5)/(12. * M_PI) - pow( gsq*pow(h,2.), 1.5)/(6. * M_PI)  + pow(m, 1.5)/(6. * M_PI);

      if ( potential_flag == 1 ) {
        return T * veffLO.real() + T * veffNLO.real();
      }

      std::complex<double> veffNNLO = (3.*lam*(m + lam*phisq))/(64.*pow(M_PI,2)) + (lam*sqrt(m + lam*phisq)*sqrt(m + 3.*lam*phisq))/(32.*pow(M_PI,2)) + (3.*lam*(m + 3.*lam*phisq))/(64.*pow(M_PI,2)) + (gsq*sqrt(m + lam*phisq)*Yphisq*sqrt(gsq*phisq*Yphisq))/(16.*pow(M_PI,2)) + (gsq*sqrt(m + 3.*lam*phisq)*Yphisq*sqrt(gsq*phisq*Yphisq))/(16.*pow(M_PI,2)) - (3.*pow(lam,2)*phisq*(0.5 + log(mu3US/(3.*sqrt(m + 3.*lam*phisq)))))/(16.*pow(M_PI,2)) - (pow(lam,2)*phisq*(0.5 + log(mu3US/(2.*sqrt(m + lam*phisq) + sqrt(m + 3.*lam*phisq)))))/(16.*pow(M_PI,2)) + ((gsq*phisq*sqrt(m + lam*phisq)*sqrt(m + 3.*lam*phisq)*Yphisq)/(16.*pow(M_PI,2)) + (sqrt(m + 3.*lam*phisq)*sqrt(gsq*phisq*Yphisq)* (-2.*lam*phisq - gsq*phisq*Yphisq))/(16.*pow(M_PI,2)) + (sqrt(m + lam*phisq)*sqrt(gsq*phisq*Yphisq)*(2.*lam*phisq - gsq*phisq*Yphisq))/(16.*pow(M_PI,2)) + (pow(lam,2)*pow(phisq,2)*(0.5 + log(mu3US/(sqrt(m + lam*phisq) + sqrt(m + 3.*lam*phisq)))))/(4.*pow(M_PI,2)) - ((pow(m + lam*phisq,2) - 2.*(m + lam*phisq)*(m + 3.*lam*phisq) + pow(m + 3.*lam*phisq,2) - 2.*gsq*phisq*(m + lam*phisq)*Yphisq -2.*gsq*phisq*(m + 3.*lam*phisq)*Yphisq + pow(gsq,2)*pow(phisq,2)*pow(Yphisq,2))*(0.5 + log(mu3US/(sqrt(m + lam*phisq) + sqrt(m + 3.*lam*phisq) +sqrt(gsq*phisq*Yphisq)))))/(16.*pow(M_PI,2)))/(2.*phisq) +pow(gsq,2)*phisq*pow(Yphisq,2)*(1/(32.*pow(M_PI,2)) - (0.5 + log(mu3US/(sqrt(m + 3.*lam*phisq) + 2.*sqrt(gsq*phisq*Yphisq))))/(16.*pow(M_PI,2)) + ((gsq*phisq*sqrt(m + 3.*lam*phisq)*Yphisq*sqrt(gsq*phisq*Yphisq))/(8.*pow(M_PI,2)) +(gsq*phisq*Yphisq*(m + 3.*lam*phisq - 2.*gsq*phisq*Yphisq))/(16.*pow(M_PI,2)) - (pow(m + 3.*lam*phisq,2)*(0.5 + log(mu3US/sqrt(m + 3.*lam*phisq))))/(16.*pow(M_PI,2)) +(pow(-m - 3.*lam*phisq + gsq*phisq*Yphisq,2)*(0.5 + log(mu3US/(sqrt(m + 3.*lam*phisq) + sqrt(gsq*phisq*Yphisq)))))/(8.*pow(M_PI,2)) -(pow(-m - 3.*lam*phisq + 2.*gsq*phisq*Yphisq,2)*(0.5 + log(mu3US/(sqrt(m + 3.*lam*phisq) + 2.*sqrt(gsq*phisq*Yphisq)))))/(16.*pow(M_PI,2)))/(4.*pow(gsq,2)*pow(Yphisq,2)*pow(phisq,2)));
      
      if ( potential_flag == 2 ) {
        return T * veffLO.real() + T * veffNLO.real() + T * veffNNLO.real();
      }
      
      return 0.;
    }

    void Betas( const std::vector<double>& x, std::vector<double>& dxdt, const double t) override {
    
      double gsq = x[0];
      double m = x[1];
      double lam = x[2];
      dxdt[0] = 1./t * 0.004221715985097407*pow(gsq,2);
      dxdt[1] = 1./t * -0.012665147955292222*(3.*gsq - 4.*lam)*m;
      dxdt[2] = 1./t * 0.012665147955292222*(3.*pow(gsq,2) - 6.*gsq*lam + 10.*pow(lam,2));
    
    }

    std::vector<double> get_3d_parameters(double T) const {
    
      double Gamma = scaleFactor * T;
      double scaleFactor3;
      double scaleFactor3dUS;
      double Lb = 2. * log(scaleFactor) + 2. * EulerGamma - 2. * log(4. * M_PI);
      double Lf = Lb + 4. * log(2.);
    
      double gsq = alglib::spline1dcalc(RGEs[0], Gamma);
      double m = alglib::spline1dcalc(RGEs[1], Gamma);
      double lam = alglib::spline1dcalc(RGEs[2], Gamma);
      //---------------------------------------------------------------------------------
      double gsq3d = gsq*T - 0.0021108579925487037*pow(gsq,2)*Lb*T;
      double lam3d = 0.006332573977646111*(pow(gsq,2)*(2. - 3.*Lb) + 6.*gsq*lam*Lb + 2.*lam*(78.95683520871486 - 5.*lam*Lb))*T;
      //---------------------------------------------------------------------------------
      double lambdaVLL_1 = 0.10132118364233778*pow(gsq,2)*T;
      double lambdaVL_1 = 0.004221715985097407*gsq*(473.7410112522892 + 24.*lam - 1.*gsq*(-4. + Lb))*T;
      //---------------------------------------------------------------------------------
      double MusqU1_LO = 0.3333333333333333*gsq*pow(T,2);
      double MusqU1_NLO = 0.025330295910584444*gsq*m + 0.004113055116159166*pow(gsq,2)*pow(T,2) + 0.008443431970194815*gsq*lam*pow(T,2) + 0.0014072386616991357*pow(gsq,2)*pow(T,2)*log(12.566370614359172*T) - 0.0014072386616991357*pow(gsq,2)*pow(T,2)*log(scaleFactor*T);
    
      double MusqU1 = MusqU1_LO + MusqU1_NLO;
      //---------------------------------------------------------------------------------
      scaleFactor3 = sqrt(MusqU1);
      scaleFactor3dUS = sqrt(gsq) * scaleFactor3;
      //---------------------------------------------------------------------------------
      double m3d_LO = m + 0.08333333333333333*(3.*gsq + 4.*lam)*pow(T,2);
      double m3d_NLO = 0.00017590483271239196*(pow(gsq,2)*(252.04651042641888 + 69.*Lb)*pow(T,2) + 12.*gsq*(-26.894056714046535*lam*pow(T,2) + Lb*(9.*m - 6.*lam*pow(T,2))) + 24.*lam*(14.44702835702327*lam*pow(T,2) + Lb*(-6.*m + lam*pow(T,2))) + 18.*log(scaleFactor3/Gamma)*(8.*pow(gsq3d,2) - 16.*gsq3d*lam3d + 16.*pow(lam3d,2) + pow(lambdaVL_1,2)));
    
      double m3d = m3d_LO + m3d_NLO;
      //---------------------------------------------------------------------------------
      if ( matching_flag == 0 ) {
        return {gsq*T, m3d_LO, lam*T};
      } else if ( matching_flag == 1 ) {
        return {gsq3d, m3d, lam3d};
      }
      //---------------------------------------------------------------------------------
      double gsq3dUS = gsq3d;
      double lam3dUS = lam3d - (0.009947183943243459*pow(lambdaVL_1,2))/sqrt(MusqU1);
      //---------------------------------------------------------------------------------
      double m3dUS_LO = m3d - 0.039788735772973836*sqrt(MusqU1)*lambdaVL_1;
      double m3dUS_NLO = 0.0007915717472057639*lambdaVL_1*(-2.*lambdaVL_1 - 4.*log((0.5*scaleFactor3)/sqrt(MusqU1))*lambdaVL_1 + lambdaVLL_1);
      double m3dUS_beta = 0.025330295910584444*(pow(gsq3dUS,2) - 2.*gsq3dUS*lam3dUS + 2.*pow(lam3dUS,2));
    
      double m3dUS = m3dUS_LO + m3dUS_NLO + m3dUS_beta * log(scaleFactor3dUS / scaleFactor3);
      //---------------------------------------------------------------------------------
      if ( matching_flag == 2 ) {
        return {gsq3dUS, m3dUS, lam3dUS};
      }
      return {0.};
    }

  private :

    PROPERTY(int, potential_flag, 1);
    PROPERTY(int, matching_flag, 1);
    PROPERTY(bool, running_flag, true);
  
    double gsq_input;
    double m_input;
    double lam_input;
  
    double scaleFactor = M_PI;
};

} // namespace EffectivePotential

#endif
