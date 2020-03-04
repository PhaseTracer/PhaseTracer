#include "thermal_function.hpp"
#include "thermal_function_tables.hpp"

namespace EffectivePotential {

alglib::spline1dinterpolant make_cubic_spline(alglib::real_1d_array x, alglib::real_1d_array y) {
  alglib::spline1dinterpolant spline;
  alglib::spline1dbuildcubic(x, y, spline);
  return spline;
}

double J_B(double x) {
  constexpr double min_x = -3.72402637;
  constexpr double max_x = 1.41e3;

  static const auto spline = make_cubic_spline(J_B_X_DATA, J_B_Y_DATA);
  static const double spline_at_min_x = alglib::spline1dcalc(spline, min_x);

  if (x < min_x) {
    return spline_at_min_x;
  } else if (x > max_x) {
    return 0.;
  } else {
    return alglib::spline1dcalc(spline, x);
  }
}

double J_F(double x) {
  constexpr double min_x = -6.82200203;
  constexpr double max_x = 1.35e3;

  static const auto spline = make_cubic_spline(J_F_X_DATA, J_F_Y_DATA);
  static const double spline_at_min_x = alglib::spline1dcalc(spline, min_x);

  if (x < min_x) {
    return spline_at_min_x;
  } else if (x > max_x) {
    return 0.;
  } else {
    return alglib::spline1dcalc(spline, x);
  }
}

}  // namespace EffectivePotential
