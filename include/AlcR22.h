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
//  AlcR22.h
//  Created by Nico Petri on 01.11.20.


#ifndef ALC_R22_H
#define ALC_R22_H

/**
 * Class for the reference implementation of the density formula for alcohol-water mixtures (AWM).
 * according to OIML from the year 1973 (OIMLR22).
 */
class AlcR22 {
  
public:
    
  /**
   * Maximum possible density is that of water at 3.96913770 °C in kg/㎥.
   * @return As already described.
   */
  inline const double maxDensity() const;
  
  
  /**
   * Minimum possible density is pure alcohol at 40.0 °C in kg/㎥.
   * @return As already described.
   */
  inline const double minDensity() const;
  
  /**
   * Calculates the density of a water alcohol mixture.
   *
   * @param T Temperatur in °C
   * @param Mas The mass concentration.
   * @return Density in kg/㎥
   */
  const double density(const double Mas, const double T) const;

  /**
   * Calculates the density of water at specified temperature.
   *
   * @param T Temperatur in °C
   * @return Density in kg/㎥
   */
  const double densityOfWater(const double T) const;
  
protected:

  /**
   * Coefficients of the 1st polynomial.
   */
  const double _A_[12] = { 998.20123, -192.9769495, 389.1238958,
                          -1668.103923, 13522.15441 ,-88292.78388,
                           306287.4042, -613838.1234, 747017.2998,
                           -547846.1354, 223446.0334, -39032.85426 };
  
  /**
   * Coefficients of the 2nd polynomial.
   */
  const double _B_[6] = { -2.0618513e-1, -5.2682542e-3, 3.6130013e-5,
                           -3.8957702e-7, 7.169354e-9, -9.9739231e-11 };
  /**
   * Coefficients of the 3rd polynomial.
   */
  const double _C_ [5][11] = {{ 1.693443461530087e-1, -1.046914743455169e+1, 7.196353469546523e+1,
                               -7.047478054272792e+2, 3.924090430035045e+3, -1.210164659068747e+4,
                                2.248646550400788e+4, -2.605562982188164e+4, 1.852373922069467e+4,
                               -7.420201433430137e+3, 1.285617841998974e+3 },
    
                              { -1.19301300505701e-2, 2.517399633803461e-1, -2.170575700563993,
                                 1.353034988843029e+1, -5.029988758537014e+1, 1.09635566657757e+2,
                                -1.422753946421155e+2, 1.08043594285623e+2, -4.414153236817392e+1,
                                 7.442971530188783 },
                              { -6.802995733503803e-4, 1.876837790289664e-2, -0.2002561813734156,
                                 1.02299296671922, -2.895696483903638, 4.810060584300675,
                                -4.672147440794683, 2.458043105903461, -5.411227621436812e-1},
                              { 4.075376675622027e-6, -8.76305857347111e-6, 6.515031360099368e-6,
                               -1.51578483698721e-6 },
                              {-2.788074354782409e-08, 1.345612883493354e-08 }};

  /**
   * Calculates the first part of the polynomial, which corresponds to the density of the alcohol-water mixture at 20 °C.
   *
   * @param Mas The mass concentration.
   * @return As already described.
   */
  const double firstPol(const double Mas) const;
  
  /**
   * Calculates the second part of the polynomial.
   *
   * @param T The temperature in °C
   * @return As already described.
   */
  const double secondPol(const double T) const;
  
  
  /**
   * Calculates the third part of the polynomial.
   *
   * @param Mas The mass concentration.
   * @param T The temperature in °C
   * @return As already described.
   */
  const double thirdPol(const double Mas, const double T) const;
  
};

#endif /* ALC_R22_H */
