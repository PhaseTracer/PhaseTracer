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
  Toy model with with supercooling.
*/

#include "one_loop_potential.hpp"
#include "property.hpp"
#include "pow.hpp"
#include "logger.hpp"
#include <math.h>
#include <iostream>
#include <fstream>
#include "boost/filesystem.hpp"

namespace EffectivePotential {
class SuperCoolModel : public OneLoopPotential {
private:
  void init(double kappa_, double mu_sq_, double lambda_, double mu0_sq_) {
    mu_sq = mu_sq_;
    kappa = kappa_;
    lambda = lambda_;
    mu0_sq = mu0_sq_;
    vSq = pow(get_renormalization_scale(), 2);
  }

public:
  SuperCoolModel(double kappa_, double mu_sq, double lambda_, double mu0_sq_) {
    init(mu_sq, kappa_, lambda_, mu0_sq_);
  }

  SuperCoolModel(std::string inputFileName) {
    std::ifstream inputFile(inputFileName);

    if (!inputFile) {
      std::cerr << "Cannot open the file: " << inputFileName << std::endl;
      std::cerr << "Absolute path: " << boost::filesystem::complete(inputFileName) << std::endl;
      return;
    }

    std::string line;
    std::vector<double> data;

    std::cout.precision(std::numeric_limits<double>::max_digits10);

    int lineIndex = 0;

    int dataValues = 4;

    std::getline(inputFile, line);
    if (line.size() == 0) {
      std::cerr << "Parameter point line is empty!" << std::endl;
      return;
    }

    // Keep only the first 3 values, in case there is other data stored in the file.
    data = split(line, ' ', dataValues);

    inputFile.close();

    init(data[0], data[1], data[2], data[3]);

    LOG(debug) << "Params: " << std::endl;
    LOG(debug) << mu_sq << " " << kappa << " " << lambda << " " << mu0_sq << std::endl;
  }

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

  /*double V(Eigen::VectorXd phi, double T) const override
  {
          // Subtract off the radiation energy density from the light particles that are not included explicitly in the one-loop corrections.
          return OneLoopPotential::V(phi, T) - M_PI*M_PI/90.*raddof*pow(T, 4);
  }*/

  double V0(Eigen::VectorXd phi) const override {
    return -0.5 * mu_sq * pow(phi[0], 2) + kappa * pow(phi[0], 3) / 3 + 0.25 * lambda * pow(phi[0], 4);
  }

  size_t get_n_scalars() const override {
    return 1;
  }

  // bool forbidden(Eigen::VectorXd phi) const override { return phi[0] < -5.; } // todo: why?

  std::vector<double> get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const override {
    const double higgsSq = 3. * lambda * pow(phi[0], 2) + 2. * kappa * phi[0] - mu0_sq;
    // LOG(debug) << "higgs mass: " << higgsSq << std::endl;
    return {higgsSq};
  }

  std::vector<double> get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const override {
    const double higgsSqZeroT = get_scalar_masses_sq(phi, xi)[0];
    const double T2 = T * T;
    const double T2forH = bUseBoltzmannSuppression ? (T2 == 0. || higgsSqZeroT / T2 > MAX_BOLTZ_EXP ? 0. : T2 * exp(-higgsSqZeroT / T2)) : T2;
    const double higgsSq = higgsSqZeroT + T2forH * (lambda / 4. + g * g + (g * g + g1 * g1) / 16. + yt * yt / 4.);
    return {higgsSq};
  }

  std::vector<double> get_scalar_dofs() const override {
    return {1};
  }

  std::vector<double> get_fermion_masses_sq(Eigen::VectorXd phi) const override {
    const double hSq = phi[0] * phi[0] / vSq;

    // Mass values from PDG.
    const double topSq = 162.5 * 162.5 * hSq;
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

  std::vector<double> get_vector_masses_sq(Eigen::VectorXd phi) const override {
    const double MW2_T = pow(0.5 * g * phi[0], 2);
    const double MW2_L = MW2_T;
    const double MZ2_T = (pow(g, 2) + pow(g1, 2)) / 4. * pow(phi[0], 2);
    const double MZ2_L = MZ2_T;
    const double MPh2_L = 0.;

    // LOG(debug) << "bosons: " << MW2_T << " " << MW2_L << " " << MZ2_T << " " << MZ2_L << " " << MPh2_L << std::endl;

    return {MW2_T, MW2_L, MZ2_T, MZ2_L, MPh2_L};
  }

  std::vector<double> get_vector_debye_sq(Eigen::VectorXd phi, double T) const override {
    const std::vector<double> massesZeroT = get_vector_masses_sq(phi);

    const double h2 = phi[0] * phi[0];
    const double T2 = T * T;

    const double MW2_T = pow(0.5 * g * phi[0], 2);
    const double T2forW = bUseBoltzmannSuppression ? (T2 == 0. || massesZeroT[1] / T2 > MAX_BOLTZ_EXP ? 0. : T2 * exp(-massesZeroT[1] / T2)) : T2;
    const double MW2_L = MW2_T + 11. / 6. * g * g * T2forW;
    // const double T2forZ = T2 * (bUseBoltzmannSuppression && T2 > 0. ? exp(-massesZeroT[3]/T2) : 1.);
    const double T2forZ = bUseBoltzmannSuppression ? (T2 == 0. || massesZeroT[3] / T2 > MAX_BOLTZ_EXP ? 0. : T2 * exp(-massesZeroT[3] / T2)) : T2;
    const double a = (g * g + g1 * g1) * (3 * h2 + 22 * T2forZ);
    const double b = std::sqrt(9 * pow(g * g + g1 * g1, 2) * h2 * h2 + 44 * T2forZ * pow(g * g - g1 * g1, 2) * (3 * h2 + 11 * T2forZ));
    const double MZ2_T = (pow(g, 2) + pow(g1, 2)) / 4. * pow(phi[0], 2);
    const double MZ2_L = (a + b) / 24.;
    const double MPh2_L = (a - b) / 24. + (T == 0. ? 1e-10 : 0.);

    // LOG(debug) << "bosons (debye) (x=" << phi[0] << ", T=" << T << "): " << MW2_T << " " << MW2_L << " " <<
    //	MZ2_T << " " << MZ2_L << " " << MPh2_L << std::endl;

    return {MW2_T, MW2_L, MZ2_T, MZ2_L, MPh2_L};
  }

  std::vector<double> get_vector_dofs() const override {
    return {4., 2., 2., 1., 1.};
  }

private:
  /*double mu2() const
  {
          return  0.5 * (pow(mh, 2) + get_renormalization_scale() * kappa);
  }

  double lambda() const
  {
          return 0.5 / pow(get_renormalization_scale(), 2) * (pow(mh, 2) - get_renormalization_scale() * kappa);
  }*/

  PROPERTY(double, mh, 125.)
  PROPERTY(double, mu_sq, 0)
  PROPERTY(double, kappa, 0)
  PROPERTY(double, lambda, 0)
  PROPERTY(double, mu0_sq, 0)
  PROPERTY(double, vSq, 0)
  PROPERTY(double, yt, 0.9946)
  PROPERTY(double, g, 0.6535)
  PROPERTY(double, g1, 0.35)
  PROPERTY(double, raddof, 25.75)
  PROPERTY(bool, bUseBoltzmannSuppression, false)
  PROPERTY(double, MAX_BOLTZ_EXP, 12.);
};

} // namespace EffectivePotential