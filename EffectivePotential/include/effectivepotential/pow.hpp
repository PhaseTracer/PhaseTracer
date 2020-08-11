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
   Optimized pow(T, int) etc
*/

template<class T>
T square(T x) {
  return x * x;
}

template<class T>
T cube(T x) {
  return x * x * x;
}

template<class T>
T pow_4(T x) {
  x *= x;
  return x * x;
}

template<class T>
T pow_int(T x, int n) {

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
