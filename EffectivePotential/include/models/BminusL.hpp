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

// This is the effective potential for a minimal confromal B-L model (lps =0)
// and its singlet extension (lps != 0)  implemented for https://arxiv.org/abs/2511.10288
// Model written and provided by Zhao-Yang Zhang, Satyabrata Datta and Peter Athron
// Uses simple treatment of resumming logs locally around v
// Added to illustrate a model that is relevant over many different input mass scales

#ifndef POTENTIAL_BminusL_HPP_INCLUDED
#define POTENTIAL_BminusL_HPP_INCLUDED
#include "potential.hpp"
#include "property.hpp"
#include "pow.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include "one_loop_potential.hpp"
#include "thermal_function.hpp"
#include "boost/filesystem.hpp"

namespace EffectivePotential {

class BminusL : public Potential {
private:
  double lps = 0.;
  double gbl = 0.;
  double lr1 = 0.;
  double lr2 = 0.;
  double lr3 = 0.;
  double vphi = 1.;

  void init(double lps_, double vphi_, double gbl_, double lr1_, double lr2_, double lr3_) {
    lps = lps_;
    vphi = vphi_;
    gbl = gbl_;
    lr1 = lr1_;
    lr2 = lr2_;
    lr3 = lr3_;
    field_scale = vphi;
    temperature_scale = vphi;

    // The step size for deruvatives should scale with the general field / mass scale
    // which is set by the VEV.
    const double small_scale = vphi * 1e-3;
    set_h(0.001 * small_scale);
  }

public:
  BminusL(double lps_, double vphi_, double gbl_, double lr1_, double lr2_, double lr3_) {
    init(lps_, vphi_, gbl_, lr1_, lr2_, lr3_);
  }

  BminusL(std::string inputFileName) {
    std::ifstream inputFile(inputFileName);

    if (!inputFile) {
      std::cerr << "Cannot open the file: " << inputFileName << std::endl;
      std::cerr << "Absolute path: " << boost::filesystem::complete(inputFileName) << std::endl;
      return;
    }

    std::string line;
    std::vector<double> data;

    std::cout.precision(std::numeric_limits<double>::max_digits10);

    std::getline(inputFile, line);
    if (line.size() == 0) {
      std::cerr << "Parameter point line is empty!" << std::endl;
      return;
    }

    // Keep only the first 6 values, in case there is other data stored in the file.
    int dataValues = 6;
    data = split(line, ' ', dataValues);
    inputFile.close();
    init(data[0], data[1], data[2], data[3], data[4], data[5]);
  }

  BminusL(const std::vector<double> &point) {
    if (point.size() < 6) {
      throw std::runtime_error("BminusL expects at least 6 input parameters.");
    }
    init(point[0], point[1], point[2], point[3], point[4], point[5]);
  }

  // Add a getter for the VEV so we can use this to re-scale settings
  double get_vphi() const { return vphi; }

  double get_raddof() const override { return 106.0; }
  double get_minimum_temperature() const override { return 0.01; }

  double V(Eigen::VectorXd phi, double T) const override {
    const double T2 = T * T;
    const double vphi2 = vphi * vphi;
    const double phi_sq = phi[0] * phi[0];
    const double T4 = T2 * T2;

    const double IPi2 = 0.10132118364233777144;    //= 1. / (M_PI * M_PI);
    const double I12PiSqrt2 = 0.01875658991993971; // = 1. / (12. * M_PI * std::sqrt(2.));
    const double twoI3Pi = 0.2122065907891938;     // 2. / (3. * M_PI);
    const double I96 = 0.01041666666666667;        // 1. / 96.;
    const double I12 = 0.08333333333333333;        // 1. / 12.;
    const double sixIPi2 = 0.6079271018540266;

    const double B_CW = sixIPi2 * (lps * lps * I96 + gbl * gbl * gbl * gbl - (lr1 * lr1 * lr1 * lr1 + lr2 * lr2 * lr2 * lr2 + lr3 * lr3 * lr3 * lr3) * I96);

    const double VzeroT = 0.25 * B_CW * (0.5 * vphi2 * phi_sq * xlogx(phi_sq / vphi2) - 0.25 * phi_sq * phi_sq);

    if (T == 0.0)
      return VzeroT;
    else
      return VzeroT + IPi2 * T4 * J_B((0.5 * lps * phi_sq) / (T2)) + 1.5 * IPi2 * T4 * J_B((4 * gbl * gbl * phi_sq) / (T2)) + IPi2 * T4 * J_F((lr1 * lr1 * phi_sq) / (2 * T2)) + IPi2 * T4 * J_F((lr2 * lr2 * phi_sq) / (2 * T2)) + IPi2 * T4 * J_F((lr3 * lr3 * phi_sq) / (2 * T2)) - I12PiSqrt2 * std::pow(T2, 0.5) * std::pow(lps, 1.5) * (std::pow(phi_sq + I12 * T2, 1.5) - abs(phi_sq * std::pow((phi_sq), 0.5))) - twoI3Pi * std::pow(T2, 0.5) * std::pow(gbl, 3) * (std::pow(phi_sq + T2, 1.5) - std::pow(abs(phi_sq), 1.5));
  }
  size_t get_n_scalars() const override { return 1; }
  std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
    return {-phi};
  };
  bool forbidden(Eigen::VectorXd x) const override { return false; }

  // Based on https://stackoverflow.com/questions/236129/how-do-i-iterate-over-the-words-of-a-string
  template <typename DataType>
  void split(const std::string &inputString, char delimiter, DataType result, int keepCount) {
    std::istringstream inputStringStream(inputString);
    std::string substring;
    int i = 0;

    while (std::getline(inputStringStream, substring, delimiter) && i++ < keepCount) {
      *result++ = std::stod(substring);
    }
  }

  std::vector<double> split(const std::string &inputString, char delimiter, int keepCount) {
    std::vector<double> elements;
    split(inputString, delimiter, std::back_inserter(elements), keepCount);
    return elements;
  }
}; // class BLModel
} // namespace EffectivePotential

#endif
