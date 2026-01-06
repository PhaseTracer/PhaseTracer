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

#ifndef POTENTIAL_RSS_HPP_INCLUDED
#define POTENTIAL_RSS_HPP_INCLUDED

/**
 * Set check_vacuum_at_high to false, otherwise we'll get errors for minima away from the origin at T>>Tc.
 * The origin is not a minimum at high temperature due to explicit Z2-symmetry-breaking terms.
 */

/**
 * Real scalar singlet extension of the Standard Model, not assuming Z2 symmetry.
 */

#include "one_loop_potential.hpp"
#include "pow.hpp"
#include "SM_parameters.hpp"

#include <vector>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
#include <iomanip>
#include <math.h>

namespace EffectivePotential {
class RSS : public OneLoopPotential {
public:
  /**
   * @brief Make a real scalar singlet model (with no Z2 symmetry) from Lagrangian parameters.
   */
  RSS(double muh_sq_, double mus_sq_, double lh_, double ls_, double c1_, double c2_, double b3_, double muh_sq_tree_EWSB_, double mus_sq_tree_EWSB_, double vs_, double ms_, double theta_) : muh_sq(muh_sq_), mus_sq(mus_sq_), lh(lh_), ls(ls_), c1(c1_), c2(c2_), b3(b3_),
                                                                                                                                                                                               muh_sq_tree_EWSB(muh_sq_tree_EWSB_), mus_sq_tree_EWSB(mus_sq_tree_EWSB_), vs(vs_), theta_vev(theta_) {
    ms2 = ms_ * ms_;

    cos2Theta_vev = std::cos(2 * theta_vev);
    sin2Theta_vev = std::sin(2 * theta_vev);
    cosSqTheta_vev = std::cos(theta_vev);
    cosSqTheta_vev *= cosSqTheta_vev;
    sinSqTheta_vev = 1 - cosSqTheta_vev;

    init();
  }

  RSS(const std::vector<double> &params) {
    init();
    setParameters(params);
  }

  /**
   * Default constructor. Use at your own risk. There are no checks when using the potential that it has been
   * properly initialised with parameter values. If using this default constructor, it is assumed you will use
   * setParameters afterwards.
   */
  RSS() {
    init();
  }

  void init() {
    setCouplingsAndScale();

    bIgnoreGoldstone = false;

    goldstoneStepSize = 1;
  }

  void setCouplingsAndScale() {
    // From PDG (04/06/2021).
    double mW = 80.379;
    double mZ = 91.1876;
    double mt = 162.5; // Using running mass rather than pole mass.
    double alpha = 1. / 137.036;
    mh2 = 125.1 * 125.1;

    double e = std::sqrt(4 * M_PI * alpha);
    double sinWeinbergAngle = std::sqrt(1 - pow(mW / mZ, 2));
    double cosWeinbergAngle = mW / mZ;

    double gp = e / cosWeinbergAngle;
    gp2 = gp * gp;
    double g = e / sinWeinbergAngle;
    g2 = g * g;

    vh = 2 * mW / g;

    double yt = std::sqrt(2) * mt / vh;
    yt2 = yt * yt;

    set_renormalization_scale(mZ);

    field_scale = vh;
    temperature_scale = vh;
  }

  /** Used to consider a new potential without creating a new model object. */
  void setParameters(double muh_sq_, double mus_sq_, double lh_, double ls_, double c1_, double c2_, double b3_,
                     double muh_sq_tree_EWSB_, double mus_sq_tree_EWSB_, double vs_, double ms_, double theta_) {
    muh_sq = muh_sq_;
    mus_sq = mus_sq_;
    lh = lh_;
    ls = ls_;
    c1 = c1_;
    c2 = c2_;
    b3 = b3_;
    muh_sq_tree_EWSB = muh_sq_tree_EWSB_;
    mus_sq_tree_EWSB = mus_sq_tree_EWSB_;
    vs = vs_;
    ms2 = ms_ * ms_;
    theta_vev = theta_;

    cos2Theta_vev = std::cos(2 * theta_vev);
    sin2Theta_vev = std::sin(2 * theta_vev);
    cosSqTheta_vev = std::cos(theta_vev);
    cosSqTheta_vev *= cosSqTheta_vev;
    sinSqTheta_vev = 1 - cosSqTheta_vev;
  }

  virtual void setParameters(std::vector<double> params) {
    if (params.size() != 12) {
      std::cerr << "Called RSS::setParameters with the wrong number of parameters: got " << params.size() << ", expected 12." << std::endl;
      return;
    }

    muh_sq = params[5];
    mus_sq = params[6];
    lh = params[7];
    ls = params[8];
    c1 = params[0];
    c2 = params[9];
    b3 = params[1];
    muh_sq_tree_EWSB = params[10];
    mus_sq_tree_EWSB = params[11];
    vs = params[3];
    ms2 = params[4] * params[4];
    theta_vev = params[2];

    cos2Theta_vev = std::cos(2 * theta_vev);
    sin2Theta_vev = std::sin(2 * theta_vev);
    cosSqTheta_vev = std::cos(theta_vev);
    cosSqTheta_vev *= cosSqTheta_vev;
    sinSqTheta_vev = 1 - cosSqTheta_vev;
  }

  virtual void printParameters() const {
    auto precisionPrev = std::cout.precision(std::numeric_limits<double>::max_digits10);

    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "muh_sq: " << muh_sq << std::endl;
    std::cout << "mus_sq: " << mus_sq << std::endl;
    std::cout << "lh: " << lh << std::endl;
    std::cout << "ls: " << ls << std::endl;
    std::cout << "c1: " << c1 << std::endl;
    std::cout << "c2: " << c2 << std::endl;
    std::cout << "b3: " << b3 << std::endl;
    std::cout << "muh_sq_tree_EWSB: " << muh_sq_tree_EWSB << std::endl;
    std::cout << "mus_sq_tree_EWSB: " << mus_sq_tree_EWSB << std::endl;
    std::cout << "vs: " << vs << std::endl;
    std::cout << "ms: " << std::sqrt(ms2) << std::endl;
    std::cout << "theta: " << theta_vev << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;

    std::cout.precision(precisionPrev);
  }

  std::vector<Eigen::VectorXd> get_low_t_phases() const override {
    return {get_EW_VEV(), apply_symmetry(get_EW_VEV())[0]};
  }

  Eigen::VectorXd get_EW_VEV() const {
    Eigen::VectorXd vev(2);
    vev << vh, vs;

    return vev;
  }

  double get_vh() {
    return vh;
  }

  void set_useGSResummation(bool flag) {
    bUseGoldstoneResummation = flag;
  }

  double V0(Eigen::VectorXd phi) const override {
    double h = phi[0];
    double s = phi[1];
    double h2 = h * h;
    double s2 = s * s;

    /*std::cout << "vh: " << vh << std::endl;
    std::cout << "muh_sq*h2 " << muh_sq*h2 << std::endl;
    std::cout << "mus_sq*s2 " << mus_sq*s2 << std::endl;
    std::cout << "lh*h2*h2  " << lh*h2*h2  << std::endl;
    std::cout << "ls*s2*s2  " << ls*s2*s2  << std::endl;
    std::cout << "c1*h2*s   " << c1*h2*s   << std::endl;
    std::cout << "c2*h2*s2  " << c2*h2*s2  << std::endl;
    std::cout << "b3*s2*s   " << b3*s2*s   << std::endl;*/

    return muh_sq * h2 + mus_sq * s2 + lh * h2 * h2 + ls * s2 * s2 + c1 * h2 * s + c2 * h2 * s2 + b3 * s2 * s;
  }

  double get_Boltzmann_suppression(double massSq, double TSq) const {
    return bUseBoltzmannSuppression ? (TSq == 0. || abs(massSq / TSq) > MAX_BOLTZ_EXP ? 0. : exp(-abs(massSq / TSq))) : 1.;
  }

  std::vector<double> get_scalar_thermal_sq(double T) const override {
    const double T2 = T * T;
    const double debyeHiggs = 3. / 16. * g2 + gp2 / 16. + yt2 / 4. + 2 * lh + c2 / 6.;
    const double debyeSinglet = 2. / 3. * c2 + ls;
    const double debyeGoldstone = debyeHiggs;

    return {debyeHiggs * T2, debyeSinglet * T2, debyeGoldstone * T2};
  }

  std::vector<double> get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const override {
    double h = phi[0];
    double s = phi[1];
    double h2 = h * h;
    double s2 = s * s;

    std::vector<double> scalar_eigs = get_scalar_eigenvalues(phi, xi);

    // double goldstone = 2*muh_sq_tree_EWSB + 4*lh*h2 + 2*c1*s + 2*c2*s2;
    double goldstone = get_goldstone_massSq(phi, xi);

    return {scalar_eigs[0], scalar_eigs[1], goldstone};
  }

  std::vector<double> get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const override {
    double h = phi[0];
    double s = phi[1];
    double h2 = h * h;
    double s2 = s * s;

    std::vector<double> scalar_eigs = get_scalar_eigenvalues_thermal(phi, xi, T);

    double goldstone = get_goldstone_mass_thermal(phi, xi, T);

    // double goldstone = 2*muh_sq_tree_EWSB + 4*lh*h2 + 2*c1*s + 2*c2*s2;
    // goldstone += (3./16.*g2 + gp2/16. + yt2/4. + 2*lh + c2/6.)*T*T;

    return {scalar_eigs[0], scalar_eigs[1], goldstone};
  }

  std::vector<double> get_scalar_eigenvalues(Eigen::VectorXd phi, double xi) const {
    double mhh = d2V0mdh2(phi);
    double mhs = d2V0mdhds(phi);
    double mss = d2V0mds2(phi);
    double theta = get_mixing_angle_supplied(mhh, mhs, mss);

    double cosSqTheta = std::cos(theta);
    cosSqTheta *= cosSqTheta;
    double sinSqTheta = 1 - cosSqTheta;
    double sin2Theta = std::sin(2 * theta);

    double mh = mhh * cosSqTheta + mhs * sin2Theta + mss * sinSqTheta;
    double ms = mhh * sinSqTheta - mhs * sin2Theta + mss * cosSqTheta;

    return {mh, ms};
  }

  std::vector<double> get_scalar_eigenvalues_thermal(Eigen::VectorXd phi, double xi, double T) const {
    std::vector<double> thermal_corrections = get_scalar_thermal_sq(T);

    double mhh = d2V0mdh2(phi); // + thermal_corrections[0];
    const double mhs = d2V0mdhds(phi);
    double mss = d2V0mds2(phi); // + thermal_corrections[1];
    const double T2 = T * T;

    // mhh += bUseBoltzmannSuppression ? (T2 == 0. || mhh/T2 > MAX_BOLTZ_EXP ? 0. : thermal_corrections[0]*exp(-mhh/T2)) : thermal_corrections[0];
    // mss += bUseBoltzmannSuppression ? (T2 == 0. || mss/T2 > MAX_BOLTZ_EXP ? 0. : thermal_corrections[1]*exp(-mss/T2)) : thermal_corrections[1];
    mhh += thermal_corrections[0] * get_Boltzmann_suppression(mhh, T2);
    mss += thermal_corrections[1] * get_Boltzmann_suppression(mss, T2);

    const double theta = get_mixing_angle_supplied(mhh, mhs, mss);
    const double cosSqTheta = square(std::cos(theta));
    const double sinSqTheta = 1 - cosSqTheta;
    const double sin2Theta = std::sin(2 * theta);

    const double mh = mhh * cosSqTheta + mhs * sin2Theta + mss * sinSqTheta;
    const double ms = mhh * sinSqTheta - mhs * sin2Theta + mss * cosSqTheta;

    return {mh, ms};
  }

  double get_goldstone_massSq(Eigen::VectorXd phi, double xi) const {
    return get_goldstone_mass_thermal(phi, xi, 0.0);
  }

  double get_goldstone_mass_thermal(Eigen::VectorXd phi, double xi, double T) const {
    if (bIgnoreGoldstone) {
      return 0.0;
    }

    std::vector<double> thermal_corrections;

    if (T > 0) {
      thermal_corrections = get_scalar_thermal_sq(T);
    } else {
      thermal_corrections = {0.0, 0.0};
    }

    double hCopy = phi[0];
    double hShift = 1e-8 * field_scale;

    if (std::abs(phi[0]) < hShift) {
      phi[0] = hShift;
    }

    double dVCWdh = 0.0;

    if (bUseGoldstoneResummation) {
      Eigen::VectorXd phi_left(2);
      Eigen::VectorXd phi_right(2);

      phi_left << phi[0] - goldstoneStepSize, phi[1];
      phi_right << phi[0] + goldstoneStepSize, phi[1];

      // Temporarily ignore the Goldstone contribution when calculating the Goldstone resummation.
      bool prevIgnoreGoldstone = bIgnoreGoldstone;
      bIgnoreGoldstone = true;
      dVCWdh = (V1(phi_right, 0.0) - V1(phi_left, 0.0)) / (2 * goldstoneStepSize);
      bIgnoreGoldstone = prevIgnoreGoldstone;
    }

    double mass = 1 / phi[0] * (dV0mdh(phi) + dVCWdh);

    phi[0] = hCopy;

    // thermal_corrections[0] = bUseBoltzmannSuppression ? (T2 == 0. || mass/T2 > MAX_BOLTZ_EXP ? 0. : thermal_corrections[0]*exp(-mass/T2))
    //	: thermalCorrections[0];
    return mass + thermal_corrections[0] * get_Boltzmann_suppression(mass, T * T);
  }

  double get_mixing_angle_supplied(double mhh, double mhs, double mss) const {
    return 0.5 * std::atan2(2 * mhs, mhh - mss);
  }

  std::vector<double> get_scalar_dofs() const override {
    return {1., 1., 3.};
  }

  std::vector<double> get_vector_masses_sq(Eigen::VectorXd phi) const override {
    const double h2 = phi[0] * phi[0];
    const double mWSq_T = 0.25 * g2 * h2;
    const double mWSq_L = mWSq_T;
    const double mZSq_T = 0.25 * (g2 + gp2) * h2;
    const double mZSq_L = mZSq_T;
    const double mPhSq_L = 0.;

    return {mWSq_L, mWSq_T, mZSq_L, mZSq_T, mPhSq_L};
  }

  std::vector<double> get_vector_debye_sq(Eigen::VectorXd phi, double T) const override {
    const double h2 = phi[0] * phi[0];
    const double T2 = T * T;

    const double mWSq_T = 0.25 * g2 * h2;
    const double mWSq_L = mWSq_T + 11. / 6. * g2 * T2 * get_Boltzmann_suppression(0.25 * g2 * h2, T2);
    const double mZSq_T = 0.25 * (g2 + gp2) * h2;

    double a, b;

    if (bUseBoltzmannSuppression && T2 > 0) {
      const double TSqForZ = T2 * get_Boltzmann_suppression(mZSq_T, T2);
      a = (g2 + gp2) * (3 * h2 + 22 * TSqForZ);
      b = std::sqrt(9 * square((g2 + gp2) * h2) + 132 * square(g2 - gp2) * h2 * TSqForZ + 484 * square((g2 - gp2) * TSqForZ));
    } else {
      a = (g2 + gp2) * (3 * h2 + 22 * T2);
      b = std::sqrt(9 * square((g2 + gp2) * h2) + 132 * square(g2 - gp2) * h2 * T2 + 484 * square((g2 - gp2) * T2));
    }

    const double mZSq_L = (a + b) / 24.;
    const double mPhSq_L = (a - b) / 24.;

    return {mWSq_L, mWSq_T, mZSq_L, mZSq_T, mPhSq_L};
  }

  std::vector<double> get_vector_dofs() const override {
    return {2., 4., 1., 2., 1.};
  }

  /*std::vector<double> get_fermion_masses_sq(Eigen::VectorXd phi) const override
  {
          return {0.5*yt2*phi[0]*phi[0]};
  }*/

  std::vector<double> get_fermion_masses_sq(Eigen::VectorXd phi) const override {
    const double hSq = phi[0] * phi[0] / (vh * vh);

    // For the top quark we can use the Yuakawa which was set already
    const double topSq = 0.5 * yt2 * phi[0] * phi[0];
    // const double topSq = 162.5*162.5*hSq;
    //  Mass values from PDG.
    const double upSq = 0.00216 * 0.00216 * hSq;
    const double downSq = 0.00467 * 0.00467 * hSq;
    const double strangeSq = 0.0934 * 0.0934 * hSq;
    const double charmSq = 1.27 * 1.27 * hSq;
    const double bottomSq = 4.18 * 4.18 * hSq;
    const double muonSq = 0.10566 * 0.10566 * hSq;
    const double tauonSq = 1.777 * 1.777 * hSq;

    return {topSq, upSq, downSq, strangeSq, charmSq, bottomSq, muonSq, tauonSq};
  }

  std::vector<double> get_fermion_dofs() const override {
    return {12., 12., 12., 12., 12., 12., 4., 4.};
  }

  size_t get_n_scalars() const override {
    return 2;
  }

  std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
    auto negH = phi;
    negH[0] = -phi[0];

    return {negH};
  }

  /* ----------------------------------------------------------------------------------------------------
   * Derivatives of the tree-level potential up to second order using tree-level EWSB to calculate the
   * quadratic parameters from the quartic parameters (i.e. mu{h,s}_sq from l{h,s}). The 'subscript' 'm'
   * in V0m signifies these forms are for use in masses that enter the Coleman-Weinberg potential, where
   * we want the parameters to be related through EWSB at tree-level only.
   * --------------------------------------------------------------------------------------------------*/

  double dV0mdh(Eigen::VectorXd phi) const {
    const double h = phi[0];
    const double s = phi[1];

    return 2 * muh_sq_tree_EWSB * h + 4 * lh * h * h * h + 2 * c1 * h * s + 2 * c2 * h * s * s;
  }

  double dV0mds(Eigen::VectorXd phi) const {
    const double h = phi[0];
    const double s = phi[1];

    return 2 * mus_sq_tree_EWSB * s + 4 * ls * s * s * s + c1 * h * h + 2 * c2 * h * h * s + 3 * b3 * s * s;
  }

  double d2V0mdh2(Eigen::VectorXd phi) const {
    const double h = phi[0];
    const double s = phi[1];

    return 2 * muh_sq_tree_EWSB + 12 * lh * h * h + 2 * c1 * s + 2 * c2 * s * s;
  }

  double d2V0mds2(Eigen::VectorXd phi) const {
    const double h = phi[0];
    const double s = phi[1];

    return 2 * mus_sq_tree_EWSB + 12 * ls * s * s + 2 * c2 * h * h + 6 * b3 * s;
  }

  double d2V0mdhds(Eigen::VectorXd phi) const {
    const double h = phi[0];
    const double s = phi[1];

    return 2 * c1 * h + 4 * c2 * h * s;
  }

  /* ----------------------------------------------------------------------------------------------------
   * Utility functions for determining the maximum temperature to sample during phase tracing.
   * --------------------------------------------------------------------------------------------------*/

  // ds is the finite difference step size.
  bool hasUniqueMinimum(double T, const std::vector<double> &samples, double ds) {
    Eigen::VectorXd point1(2), point2(2);
    point1 << 0, samples[0] + ds;
    point2 << 0, samples[0] - ds;
    double oneOnTwoDS = 1 / (2 * ds);

    double dVds = (V(point1, T) - V(point2, T)) * oneOnTwoDS;
    bool prevWasPositive = dVds > 0;
    bool foundRoot = false;

    // Loop through the samples and check for roots by finding where the sign of the derivative changes.
    for (int i = 1; i < samples.size(); ++i) {
      point1[1] = samples[i] + ds;
      point2[1] = samples[i] - ds;
      dVds = (V(point1, T) - V(point2, T)) * oneOnTwoDS;
      bool positive = dVds > 0;

      // Check for a sign change.
      if (positive != prevWasPositive) {
        // If this isn't the first sign change, we have multiple roots.
        if (foundRoot) {
          return false;
        }

        foundRoot = true;
      }

      prevWasPositive = positive;
    }

    return true;
  }

  double getMaxTemp(double startTemp, double Tmax, int numSamples, double precision, double finiteDifferenceStepSize) {
    double smax = 1500.;
    double smin = -1500.;
    std::vector<double> samples = std::vector<double>(numSamples);
    double ds = (smax - smin) / (numSamples - 1);

    for (int i = 0; i < numSamples; ++i) {
      samples[i] = smin + i * ds;
    }

    double T = startTemp;

    if (hasUniqueMinimum(T, samples, finiteDifferenceStepSize)) {
      return T;
    }

    T *= 2;

    while (T < Tmax) {
      if (hasUniqueMinimum(T, samples, finiteDifferenceStepSize)) {
        break;
      }

      T *= 2;
    }

    if (T >= Tmax) {
      return Tmax;
    }

    double Tstep = 0.25 * T;
    T -= Tstep;
    Tstep *= 0.5;
    bool lastWasUnique = true;

    while (Tstep > precision) {
      if (hasUniqueMinimum(T, samples, finiteDifferenceStepSize)) {
        T -= Tstep;
        lastWasUnique = true;
      } else {
        T += Tstep;
        lastWasUnique = false;
      }

      Tstep *= 0.5;
    }

    if (lastWasUnique) {
      return T + Tstep;
    } else {
      return T + 2 * Tstep;
    }
  }

protected:
  double muh_sq;
  double mus_sq;
  double lh;
  double ls;
  double c1;
  double c2;
  double b3;
  double muh_sq_tree_EWSB;
  double mus_sq_tree_EWSB;
  double vh;
  double vs;
  double mh2;
  double ms2;
  double theta_vev;

  double cosSqTheta_vev;
  double sinSqTheta_vev;
  double cos2Theta_vev;
  double sin2Theta_vev;

  double goldstoneStepSize;
  bool bUseGoldstoneResummation;

  double g2;
  double gp2;
  double yt2;

  // This flag is used to ignore the Goldstone contribution to the Coleman-Weinberg potential, specifically for the
  // Goldstone one-loop self-energy to avoid the IR divergence.
  mutable bool bIgnoreGoldstone;

  PROPERTY(bool, bUseBoltzmannSuppression, false)
  PROPERTY(double, MAX_BOLTZ_EXP, 12.)
  PROPERTY(bool, DEBUG, false)
}; // class RSS
} // namespace EffectivePotential

#endif
