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
