// Copyright © 2020 Nico Petri, 03172 Guben, Germany
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files
// (the “Software”), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to permit
// persons to whom the Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
//
//  AlcR22.cpp
//  Created by Nico Petri on 01.11.20.

#include <assert.h>
#include <cmath>
#include "include/AlcR22.h"

inline const double AlcR22::maxDensity() const { return 999.9688158097358; }

inline const double AlcR22::minDensity() const { return 771.9323112798857; }
  
const double AlcR22::density(const double Mas, const double T) const {
  return firstPol(Mas) + secondPol(T) + thirdPol(Mas, T);
}

const double AlcR22::densityOfWater(const double T) const {
  return _A_[0] + secondPol(T);
}

const double AlcR22::firstPol(const double Mas) const {
    
  // mass fraction value range exceeded
  assert((Mas > -1.0e-14) && (Mas < (100.0 + 1.0e-14)));
  
  double _sum = 0.0;
  
  // A(i) * p^(i-1) | i = 2..12
  if (Mas > 1.0e-14) {
    for (int k = 1; k < 12; k++) {
      _sum += _A_[k] * pow(Mas, static_cast<double>(k));
    }
  }
  return _A_[0] + _sum;
}
  
const double AlcR22::secondPol(const double T) const {
  
  // temperatur value range exceeded
  assert((T > (-20.0 - 1.0e-12)) && (T < (40.0 + 1.0e-12)));
  
  double _sum = 0.0;
  
  if (abs(T-20.0) > 1.0e-14) {
    // B(i) * (t-20°C)^i | i = 1..6
    for (int k = 0; k < 6; k++) {
      _sum += _B_[k] * pow(T - 20.0, static_cast<double>(k) + 1.0);
    }
  }
  return _sum;
}
  
  
const double AlcR22::thirdPol(const double Mas, const double T) const {
  
  // mass fraction value range exceeded
  assert((Mas > -1.0e-12) && (Mas < (100.0 + 1.0e-12)));
  // temperature value range exceeded
  assert((T > (-20.0 - 1.0e-12)) && (T < (40.0 + 1.0e-12)));
  
  double _sum = 0.0;
  int size_a[5] = { 11, 10, 9, 4, 2 };
  
  if ((abs(T - 20.0) > 1.0e-14) && (Mas > 0.0)) {
    // C(i,k) * p^(k) * (t - 20°C)^i
    for (int i = 0; i < 5; i++) {
      for ( int k = 0; k < size_a[i]; k++) {
        _sum += _C_[i][k] * pow(Mas, static_cast<double>(k) + 1.0) * pow(T - 20.0, static_cast<double>(i) + 1.0);
      }
    }
  }
  return _sum;
}
