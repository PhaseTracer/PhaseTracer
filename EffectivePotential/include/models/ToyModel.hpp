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

#ifndef POTENTIAL_ToyModel_HPP_INCLUDED
#define POTENTIAL_ToyModel_HPP_INCLUDED

#include "potential.hpp"
#include "property.hpp"
#include "pow.hpp"
#include <iostream>
#include <fstream>
#include "boost/filesystem.hpp"

namespace EffectivePotential {
class ToyModel : public Potential {
private:
  double A, AonV, v, l, E, D, T0Sq;

  void init(double AonV_, double v_, double D_, double E_) {
    v = v_;
    field_scale = v;
    temperature_scale = v;
    D = D_;
    E = E_;
    setAonV(AonV_);
  }

public:
  ToyModel(double AonV_, double v_) {
    // D and E set to the values from arxiv:1611.05853.
    init(AonV_, v_, 0.4424835167, 0.0625);
  }

  ToyModel(double AonV_, double v_, double D_, double E_) {
    init(AonV_, v_, D_, E_);
  }

  ToyModel(std::string inputFileName) {
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

    // Keep only the first 4 values, in case there is other data stored in the file.
    data = split(line, ' ', dataValues);

    inputFile.close();

    init(data[0], data[1], data[2], data[3]);
  }

  void setAonV(double AonV_) {
    AonV = AonV_;
    A = AonV * v;
    l = 0.125 + 1.5 * AonV;
    T0Sq = (l - 3. * AonV) / (2. * D) * v * v;
  }

  double V(Eigen::VectorXd phi, double T) const override {
    const double x = phi[0];
    const double xSq = x * x;

    return D * (T * T - T0Sq) * xSq - (E * T + A) * xSq * x + 0.25 * l * xSq * xSq;
  }

  size_t get_n_scalars() const override { return 1; }

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
}; // class ToyModel
} // namespace EffectivePotential

#endif
