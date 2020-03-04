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

#ifndef POTENTIAL_POW_HPP_INCLUDED
#define POTENTIAL_POW_HPP_INCLUDED

/**
   Optimized pow(double, int) etc
*/

inline double square(double x) {
  return x * x;
}

inline double cube(double x) {
  return x * x * x;
}

inline double pow_4(double x) {
  x *= x;
  return x * x;
}

inline double pow_int(double x, int n) {

  if (n < 0) {
    return pow_int(1. / x, -n);
  }

  double result = 1.;  

  while (n) {
     if ((n & 1) != 0) {
        result *= x;
     }
     n >>= 1;
     x *= x;
  }

  return result;
}

#endif
