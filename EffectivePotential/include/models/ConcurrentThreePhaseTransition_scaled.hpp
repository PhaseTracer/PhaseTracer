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

#ifndef POTENTIAL_ConcurrentThreePhaseTransition_HPP_INCLUDED
#define POTENTIAL_ConcurrentThreePhaseTransition_HPP_INCLUDED

#include "potential.hpp"
#include "property.hpp"
#include "pow.hpp"

namespace EffectivePotential {
class ConcurrentThreePhaseTransition_scaled : public Potential {
public:
  double V(Eigen::VectorXd phi, double T) const override {
    const double x = phi[0];
    const double x2 = x * x;
    const double x3 = x2 * x;
    const double x4 = x3 * x;
    const double x5 = x4 * x;
    const double x6 = x5 * x;
    const double T2 = T * T;
    const double T3 = T2 * T;
    const double T4 = T3 * T;
    const double T5 = T4 * T;

    // T1=80, T2=100
    return (5 * (-967680000000 + 7994880000 * T + 1706880000 * T2 + 3072000 * T3 - 688800 * T4 - 5043 * T5) * x) / 246016 + (5 * (-8156160000 - 5817600000 * T - 15456000 * T2 + 4879200 * T3 + 45510 * T4) * x2) / 246016 + (5 * (6717440000 + 35840000 * T - 17184000 * T2 - 218080 * T3) * x3) / 246016 + (5 * (-33024000 + 30086400 * T + 585120 * T2) * x4) / 246016 + (5 * (-20951040 - 833280 * T) * x5) / 246016 + 10 * x6;
  }

  size_t get_n_scalars() const override { return 1; }

  bool forbidden(Eigen::VectorXd x) const override { return false; }
}; // class ConcurrentThreePhaseTransition_scaled
} // namespace EffectivePotential

#endif
