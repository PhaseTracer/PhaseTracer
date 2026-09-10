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

#include "thermal_function.hpp"
#include "thermal_function_tables.hpp"

namespace EffectivePotential {

namespace {

class ThermalSpline {
public:
  ThermalSpline(const alglib::real_1d_array &x, const alglib::real_1d_array &y,
                double min_x_, double max_x_)
      : min_x(min_x_), max_x(max_x_) {
    alglib::spline1dbuildcubic(x, y, spline);
    value_at_min_x = alglib::spline1dcalc(spline, min_x);
  }

  double value(double x) const 
  {
    if (x < min_x) 
    {
      return value_at_min_x;
    } else if (x > max_x) 
    {
      return 0.;
    }
    return alglib::spline1dcalc(spline, x);
  }

  double deriv(double x) const 
  {
    if (x < min_x || x > max_x) 
    {
      return 0.;
    }
    double s, ds, d2s;
    alglib::spline1ddiff(spline, x, s, ds, d2s);
    return ds;
  }

private:
  alglib::spline1dinterpolant spline;
  double min_x;
  double max_x;
  double value_at_min_x;
};

const ThermalSpline &boson_spline() {
  static const ThermalSpline s(J_B_X_DATA, J_B_Y_DATA, -3.72402637, 1.41e3);
  return s;
}

const ThermalSpline &fermion_spline() {
  static const ThermalSpline s(J_F_X_DATA, J_F_Y_DATA, -6.82200203, 1.35e3);
  return s;
}

} // namespace

double J_B(double x) { return boson_spline().value(x); }

double J_F(double x) { return fermion_spline().value(x); }

double J_B_diff(double x) { return boson_spline().deriv(x); }

double J_F_diff(double x) { return fermion_spline().deriv(x); }

} // namespace EffectivePotential
