
/*
 * Copyright (C) Bernard Gingold, 2020-2026 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*
 prob, a C++ code which handles various discrete and continuous probability density functions (PDF).

For a discrete variable X, PDF(X) is the probability that the value X will occur; for a continuous variable, PDF(X) is the probability density of X, that is, the probability of a value between X and X+dX is PDF(X) * dX.

The corresponding cumulative density functions or "CDF"'s are also handled. For a discrete or continuous variable, CDF(X) is the probability that the variable takes on a value less than or equal to X.

In some cases, the inverse of the CDF can easily be computed. If


        X = CDF_INV ( P )
      

then we are asserting that the value X has a cumulative probability density function of P, in other words, the probability that the variable is less than or equal to X is P. If the CDF cannot be analytically inverted, there are simple ways to try to estimate the inverse. Depending on the PDF, these methods may be rapid and accurate, or not.

For most distributions, the mean or "average value" or "expected value" is also available. For a discrete variable, MEAN is simply the sum of the products X * PDF(X); for a continuous variable, MEAN is the integral of X * PDF(X) over the range. For the distributions covered here, the means are known beforehand, and no summation or integration is required.

For most distributions, the variance is available. For a discrete variable, the variance is the sum of the products ( X - MEAN )^2 * PDF(X); for a continuous variable, the variance is the integral of ( X - MEAN )^2 * PDF(X) over the range. The square root of the variance is known as the standard deviation. For the distributions covered here, the variances are often known beforehand, and no summation or integration is required.

For many of the distributions, it is possible to repeatedly request "samples", that is, a pseudorandom sequence of realizations of the PDF. These samples are always associated with an integer seed, which controls the calculation. Using the same seed as input will guarantee the same sample value on output. Ultimately, a random number generator must be invoked internally. In most cases, the current code will call a routine called R8_RANDOM or I4_RANDOM, each of which in turn calls a routine called R8_UNIFORM_01. You may prefer a different random number generator for this purpose.
Licensing:

The information on this web page is distributed under the MIT license. 
*/

#ifndef __GMS_UNIFORM_DISTRIBUTION_H__ 
#define __GMS_UNIFORM_DISTRIBUTION_H__ 030920260837

#include <cstdint>
#include <valarray>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif 
#include "GMS_config.h"

namespace gms
{

namespace math
{

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
double r8_uniform_01(std::int32_t &seed)
{
//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      unif = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
  std::int32_t k;
  double r;
  k = seed / 127773;
  seed = 16807 * ( seed - k * 127773 ) - k * 2836;
  if (__builtin_expect(seed<0,0))
  {
    seed = seed + 2147483647;
  }
  r = ( double ) ( seed ) * 4.656612875E-10;
  return r;
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
double uniform_01_sample(std::int32_t &seed)
{
//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_01_SAMPLE is a random number generator.
//
//  Discussion:
//
//    SEED = SEED * (7^5) mod (2^31 - 1)
//    UNIFORM_01_SAMPLE = SEED * / ( 2^31 - 1 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &SEED, the integer "seed" used to generate
//    the output random number, and updated in preparation for the
//    next one.  SEED should not be zero.
//
//    Output, double UNIFORM_01_SAMPLE, a random value between 0 and 1.
//
//  Local Parameters:
//
//    Local, int IA = 7^5.
//
//    Local, int IB = 2^15.
//
//    Local, int IB16 = 2^16.
//
//    Local, int IP = 2^31 - 1.
//
#ifdef _OPENMP
  static int ia = 16807;
  static int ib15 = 32768;
  static int ib16 = 65536;
  static int ip = 2147483647;
#pragma omp threadprivate(ia)
#pragma omp threadprivate(ib15)
#pragma omp threadprivate(ib16)
#pragma omp threadprivate(ip)
#else 
  thread_local static int ia = 16807;
  thread_local static int ib15 = 32768;
  thread_local static int ib16 = 65536;
  thread_local static int ip = 2147483647; 
#endif
  int iprhi;
  int ixhi;
  int k;
  int leftlo;
  int loxa;
  double temp;
//
//  Don't let SEED be 0.
//
  if (__builtin_expect(seed==0,0) )
  {
    seed = ip;
  }
//
//  Get the 15 high order bits of SEED.
//
  ixhi = seed / ib16;
//
//  Get the 16 low bits of SEED and form the low product.
//
  loxa = ( seed - ixhi * ib16 ) * ia;
//
//  Get the 15 high order bits of the low product.
//
  leftlo = loxa / ib16;
//
//  Form the 31 highest bits of the full product.
//
  iprhi = ixhi * ia + leftlo;
//
//  Get overflow past the 31st bit of full product.
//
  k = iprhi / ib15;
//
//  Assemble all the parts and presubtract IP.  The parentheses are essential.
//
  seed = ( ( ( loxa - leftlo * ib16 ) - ip ) +
    ( iprhi - k * ib15 ) * ib16 ) + k;
//
//  Add IP back in if necessary.
//
  if ( seed < 0 )
  {
    seed = seed + ip;
  }
//
//  Multiply by 1 / (2^31-1).
//
  temp = ( ( double ) seed ) * 4.656612875E-10;
  return temp;
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
void uniform_01_order_sample(std::valarray<double> &samples,std::int32_t &seed)
{
    double * __restrict__ p_samples = &samples.operator[](0);
    const std::size_t n_samples     = samples.size();
    std::size_t i;
    double u;
    double v;
    v = 1.0;
    for(i = n_samples-1ull; 0ull <= i; --i)
    {   
        const double t_i{static_cast<double>(i)};
        u = r8_uniform_01(seed);
        v = std::pow(u,1.0/t_i);
        p_samples[i] = v;
    }
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_ALWAYS_INLINE__
static inline
double uniform_cdf_inv(const double cdf,const double a,const double b)
{
//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_CDF_INV inverts the Uniform CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    A < B.
//
//    Output, double UNIFORM_CDF_INV, the corresponding argument.
//
    if(__builtin_expect(cdf<0.0,0) || __builtin_expect(1.0<cdf,0))
      { return (-1.0);}
    return (a+(b-a)*cdf);
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_ALWAYS_INLINE__
static inline
double uniform_sample(const double a,const double b,std::int32_t &seed)
{
    double cdf = r8_uniform_01(seed);
    double x   = uniform_cdf_inv(cdf,a,b);
    return (x);
}

} // math

} // gms

#endif /*__GMS_UNIFORM_DISTRIBUTION_H__*/