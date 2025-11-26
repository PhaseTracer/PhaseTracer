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

#ifndef COMPLEXOPERATORS_h
#define COMPLEXOPERATORS_h

#include <vector>
#include <cmath>
#include <complex>

std::complex<double> operator * (const int & a, const std::complex<double> & b){
  if (a != 0){
      return std::complex<double>(b.real() * a, b.imag() * a);
  }else{
      return std::complex<double>(0, 0);
  }
}

std::complex<double> operator / (const std::complex<double> & b, const int & a ){
  if (a != 0){
      return std::complex<double>(b.real() / a, b.imag() / a);
  }else{
      // std::cout<<"Divide by zero error"<<std::endl;
      return std::complex<double>(0, 0);
  }
}

std::complex<double> operator - (const int & a, const std::complex<double> & b){
  return std::complex<double>(a-b.real(), -b.imag());
}

std::complex<double> operator + (const int & a, const std::complex<double> & b){
  return std::complex<double>(a+b.real(), b.imag());
}

#endif