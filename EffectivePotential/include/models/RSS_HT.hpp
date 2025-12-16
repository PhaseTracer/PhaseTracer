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

#ifndef POTENTIAL_RSS_HT_HPP_INCLUDED
#define POTENTIAL_RSS_HT_HPP_INCLUDED

/**
 * High temperature expansion of the real scalar singlet extension of the Standard Model, not assuming Z2 symmetry.
 */

#include "one_loop_potential.hpp"
#include "pow.hpp"
#include "SM_parameters.hpp"

#include <vector>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <iomanip>
#include <math.h>

namespace EffectivePotential {
class RSS_HT : public RSS {
public:
  /**
   * @brief Make a real scalar singlet model (with no Z2 symmetry) from Lagrangian parameters.
   */
  RSS_HT(double muh_sq_, double mus_sq_, double lh_, double ls_, double c1_, double c2_, double b3_,
         double vs_, double ms_, double theta_) : RSS(muh_sq_, mus_sq_, lh_, ls_, c1_, c2_, b3_, 0., 0., vs_, ms_, theta_) {
    ms2 = ms_ * ms_;

    cos2Theta_vev = std::cos(2 * theta_vev);
    sin2Theta_vev = std::sin(2 * theta_vev);
    cosSqTheta_vev = std::cos(theta_vev);
    cosSqTheta_vev *= cosSqTheta_vev;
    sinSqTheta_vev = 1 - cosSqTheta_vev;

    init();
  }

  /**
   * Default constructor. Use at your own risk. There are no checks when using the potential that it has been
   * properly initialised with parameter values. If using this default constructor, it is assumed you will use
   * setParameters afterwards.
   */
  RSS_HT() {
    init();
  }

  /** Used to consider a new potential without creating a new model object. */
  void setParameters(double muh_sq_, double mus_sq_, double lh_, double ls_, double c1_, double c2_, double b3_, double vs_, double ms_, double theta_) {
    muh_sq = muh_sq_;
    mus_sq = mus_sq_;
    lh = lh_;
    ls = ls_;
    c1 = c1_;
    c2 = c2_;
    b3 = b3_;
    vs = vs_;
    ms2 = ms_ * ms_;
    theta_vev = theta_;

    cos2Theta_vev = std::cos(2 * theta_vev);
    sin2Theta_vev = std::sin(2 * theta_vev);
    cosSqTheta_vev = std::cos(theta_vev);
    cosSqTheta_vev *= cosSqTheta_vev;
    sinSqTheta_vev = 1 - cosSqTheta_vev;
  }

  void setParameters(std::vector<double> params) override {
    if (params.size() != 10) {
      std::cerr << "Called RSS_HT::setParameters with the wrong number of parameters: got " << params.size() << ", expected 10." << std::endl;
      return;
    }

    muh_sq = params[5];
    mus_sq = params[6];
    lh = params[7];
    ls = params[8];
    c1 = params[0];
    c2 = params[9];
    b3 = params[1];
    vs = params[3];
    ms2 = params[4] * params[4];
    theta_vev = params[2];

    cos2Theta_vev = std::cos(2 * theta_vev);
    sin2Theta_vev = std::sin(2 * theta_vev);
    cosSqTheta_vev = std::cos(theta_vev);
    cosSqTheta_vev *= cosSqTheta_vev;
    sinSqTheta_vev = 1 - cosSqTheta_vev;
  }

  void printParameters() const override {
    auto precisionPrev = std::cout.precision(std::numeric_limits<double>::max_digits10);

    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "muh_sq: " << muh_sq << std::endl;
    std::cout << "mus_sq: " << mus_sq << std::endl;
    std::cout << "lh: " << lh << std::endl;
    std::cout << "ls: " << ls << std::endl;
    std::cout << "c1: " << c1 << std::endl;
    std::cout << "c2: " << c2 << std::endl;
    std::cout << "b3: " << b3 << std::endl;
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

  double V(Eigen::VectorXd phi, double T) const override {
    double h = phi[0];
    double s = phi[1];
    double h2 = h * h;
    double s2 = s * s;
    std::vector<double> debye = get_scalar_thermal_sq(T);

    return (muh_sq + 0.5 * debye[0]) * h2 + (mus_sq + 0.5 * debye[1]) * s2 + lh * h2 * h2 + ls * s2 * s2 + c1 * h2 * s + c2 * h2 * s2 + b3 * s2 * s;
  }
}; // class RSS_HT
} // namespace EffectivePotential

#endif
