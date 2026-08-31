// Copyright (C) 2006-2025 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.


/** @file tr1/hypergeometric.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{tr1/cmath}
 */

//
// ISO C++ 14882 TR1: 5.2  Special functions
//

// Written by Edward Smith-Rowland based:
//   (1) Handbook of Mathematical Functions,
//       ed. Milton Abramowitz and Irene A. Stegun,
//       Dover Publications,
//       Section 6, pp. 555-566
//   (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl

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
#ifndef __GMS_HYPERGEOMETRIC_FUNC_H__
#define __GMS_HYPERGEOMETRIC_FUNC_H__ 310820260628

#include <cstdint>
#include <limits>
#include <cmath>
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_HYPERGEOMETRIC_FUNC_MAJOR = 1;
     static const unsigned int GMS_HYPERGEOMETRIC_FUNC_MINOR = 1;
     static const unsigned int GMS_HYPERGEOMETRIC_FUNC_MICRO = 0;
     static const unsigned int GMS_HYPERGEOMETRIC_FUNC_FULLVER =
       1000U*GMS_HYPERGEOMETRIC_FUNC_MAJOR+100U*GMS_HYPERGEOMETRIC_FUNC_MINOR+
       10U*GMS_HYPERGEOMETRIC_FUNC_MICRO;
     static const char GMS_HYPERGEOMETRIC_FUNC_CREATION_DATE[] = "31-08-2026 06:31AM +00200 (MON 31 AUG 2026 GMT+2)";
     static const char GMS_HYPERGEOMETRIC_FUNC_BUILD_DATE[]    = __DATE__; 
     static const char GMS_HYPERGEOMETRIC_FUNC_BUILD_TIME[]    = __TIME__;
     static const char GMS_HYPERGEOMETRIC_FUNC_SYNOPSIS[]      = "Adapted from the 'hypergeometric.tcc' source implementation.";

}

namespace gms 
{

namespace math 
{

 
   /**
     *   @brief This routine returns the confluent hypergeometric function
     *          by series expansion.
     * 
     *   @f[
     *     _1F_1(a;c;x) = \frac{\Gamma(c)}{\Gamma(a)}
     *                      \sum_{n=0}^{\infty}
     *                      \frac{\Gamma(a+n)}{\Gamma(c+n)}
     *                      \frac{x^n}{n!}
     *   @f]
     * 
     *   If a and b are integers and a < 0 and either b > 0 or b < a
     *   then the series is a polynomial with a finite number of
     *   terms.  If b is an integer and b <= 0 the confluent
     *   hypergeometric function is undefined.
     *
     *   @param  __a  The "numerator" parameter.
     *   @param  __c  The "denominator" parameter.
     *   @param  __x  The argument of the confluent hypergeometric function.
     *   @return  The confluent hypergeometric function, or NaN due to convergence error
     */
    template<typename _Tp>
    __ATTR_ALWAYS_INLINE__
    static inline 
    _Tp
    conf_hyperg_series(_Tp __a, _Tp __c, _Tp __x) noexcept(true)
    {
      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();

      _Tp __term = _Tp(1);
      _Tp __Fac = _Tp(1);
      const unsigned int __max_iter = 100000;
      unsigned int __i;
      for (__i = 0; __i < __max_iter; ++__i)
        {
          __term *= (__a + _Tp(__i)) * __x
                  / ((__c + _Tp(__i)) * _Tp(1 + __i));
          if (std::abs(__term) < __eps)
            {
              break;
            }
          __Fac += __term;
        }
      if (__i == __max_iter) {return (std::numeric_limits<_Tp>::quiet_NaN());}
          
      return __Fac;
    }

     /**
     *  @brief  Return the hypogeometric function @f$ _2F_1(a,b;c;x) @f$
     *          by an iterative procedure described in
     *          Luke, Algorithms for the Computation of Mathematical Functions.
     *
     *  Like the case of the 2F1 rational approximations, these are 
     *  probably guaranteed to converge for x < 0, barring gross    
     *  numerical instability in the pre-asymptotic regime.         
     */
    template<typename _Tp>
    _Tp
    __ATTR_ALWAYS_INLINE__
    static inline 
    conf_hyperg_luke(_Tp __a, _Tp __c, _Tp __xin) 
    {
      const _Tp __big = std::pow(std::numeric_limits<_Tp>::max(), _Tp(0.16L));
      const int __nmax = 20000;
      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp __x  = -__xin;
      const _Tp __x3 = __x * __x * __x;
      const _Tp __t0 = __a / __c;
      const _Tp __t1 = (__a + _Tp(1)) / (_Tp(2) * __c);
      const _Tp __t2 = (__a + _Tp(2)) / (_Tp(2) * (__c + _Tp(1)));
      _Tp __F = _Tp(1);
      _Tp __prec;

      _Tp __Bnm3 = _Tp(1);
      _Tp __Bnm2 = _Tp(1) + __t1 * __x;
      _Tp __Bnm1 = _Tp(1) + __t2 * __x * (_Tp(1) + __t1 / _Tp(3) * __x);

      _Tp __Anm3 = _Tp(1);
      _Tp __Anm2 = __Bnm2 - __t0 * __x;
      _Tp __Anm1 = __Bnm1 - __t0 * (_Tp(1) + __t2 * __x) * __x
                 + __t0 * __t1 * (__c / (__c + _Tp(1))) * __x * __x;

      int __n = 3;
      while(1)
        {
          _Tp __npam1 = _Tp(__n - 1) + __a;
          _Tp __npcm1 = _Tp(__n - 1) + __c;
          _Tp __npam2 = _Tp(__n - 2) + __a;
          _Tp __npcm2 = _Tp(__n - 2) + __c;
          _Tp __tnm1  = _Tp(2 * __n - 1);
          _Tp __tnm3  = _Tp(2 * __n - 3);
          _Tp __tnm5  = _Tp(2 * __n - 5);
          _Tp __F1 =  (_Tp(__n - 2) - __a) / (_Tp(2) * __tnm3 * __npcm1);
          _Tp __F2 =  (_Tp(__n) + __a) * __npam1
                   / (_Tp(4) * __tnm1 * __tnm3 * __npcm2 * __npcm1);
          _Tp __F3 = -__npam2 * __npam1 * (_Tp(__n - 2) - __a)
                   / (_Tp(8) * __tnm3 * __tnm3 * __tnm5
                   * (_Tp(__n - 3) + __c) * __npcm2 * __npcm1);
          _Tp __E  = -__npam1 * (_Tp(__n - 1) - __c)
                   / (_Tp(2) * __tnm3 * __npcm2 * __npcm1);

          _Tp __An = (_Tp(1) + __F1 * __x) * __Anm1
                   + (__E + __F2 * __x) * __x * __Anm2 + __F3 * __x3 * __Anm3;
          _Tp __Bn = (_Tp(1) + __F1 * __x) * __Bnm1
                   + (__E + __F2 * __x) * __x * __Bnm2 + __F3 * __x3 * __Bnm3;
          _Tp __r = __An / __Bn;

          __prec = std::abs((__F - __r) / __F);
          __F = __r;

          if (__prec < __eps || __n > __nmax)
            break;

          if (std::abs(__An) > __big || std::abs(__Bn) > __big)
            {
              __An   /= __big;
              __Bn   /= __big;
              __Anm1 /= __big;
              __Bnm1 /= __big;
              __Anm2 /= __big;
              __Bnm2 /= __big;
              __Anm3 /= __big;
              __Bnm3 /= __big;
            }
          else if (std::abs(__An) < _Tp(1) / __big
                || std::abs(__Bn) < _Tp(1) / __big)
            {
              __An   *= __big;
              __Bn   *= __big;
              __Anm1 *= __big;
              __Bnm1 *= __big;
              __Anm2 *= __big;
              __Bnm2 *= __big;
              __Anm3 *= __big;
              __Bnm3 *= __big;
            }

          ++__n;
          __Bnm3 = __Bnm2;
          __Bnm2 = __Bnm1;
          __Bnm1 = __Bn;
          __Anm3 = __Anm2;
          __Anm2 = __Anm1;
          __Anm1 = __An;
        }

      if (__n >= __nmax) {return (std::numeric_limits<_Tp>::quiet_NaN());}

      return __F;
    }

    

    /**
     *   @brief Return the hypogeometric function @f$ _2F_1(a,b;c;x) @f$
     *   by series expansion.
     * 
     *   The hypogeometric function is defined by
     *   @f[
     *     _2F_1(a,b;c;x) = \frac{\Gamma(c)}{\Gamma(a)\Gamma(b)}
     *                      \sum_{n=0}^{\infty}
     *                      \frac{\Gamma(a+n)\Gamma(b+n)}{\Gamma(c+n)}
     *                      \frac{x^n}{n!}
     *   @f]
     * 
     *   This works and it's pretty fast.
     *
     *   @param  __a  The first @a numerator parameter.
     *   @param  __a  The second @a numerator parameter.
     *   @param  __c  The @a denominator parameter.
     *   @param  __x  The argument of the confluent hypergeometric function.
     *   @return  The confluent hypergeometric function.
     */
    template<typename _Tp>
    __ATTR_ALWAYS_INLINE__
    static inline 
    _Tp
    hyperg_series(_Tp __a, _Tp __b, _Tp __c, _Tp __x) noexcept(true)
    {
      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();

      _Tp __term = _Tp(1);
      _Tp __Fabc = _Tp(1);
      const unsigned int __max_iter = 100000;
      unsigned int __i;
      for (__i = 0; __i < __max_iter; ++__i)
        {
          __term *= (__a + _Tp(__i)) * (__b + _Tp(__i)) * __x
                  / ((__c + _Tp(__i)) * _Tp(1 + __i));
          if (std::abs(__term) < __eps)
            {
              break;
            }
          __Fabc += __term;
        }
      if (__i == __max_iter) {return (std::numeric_limits<_Tp>::quiet_NaN());}

      return __Fabc;
    }

    /**
     *   @brief  Return the hypogeometric function @f$ _2F_1(a,b;c;x) @f$
     *           by an iterative procedure described in
     *           Luke, Algorithms for the Computation of Mathematical Functions.
     */
    template<typename _Tp>
    __ATTR_ALWAYS_INLINE__
    static inline 
    _Tp
    hyperg_luke(_Tp __a, _Tp __b, _Tp __c, _Tp __xin)
    {
      const _Tp __big = std::pow(std::numeric_limits<_Tp>::max(), _Tp(0.16L));
      const int __nmax = 20000;
      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp __x  = -__xin;
      const _Tp __x3 = __x * __x * __x;
      const _Tp __t0 = __a * __b / __c;
      const _Tp __t1 = (__a + _Tp(1)) * (__b + _Tp(1)) / (_Tp(2) * __c);
      const _Tp __t2 = (__a + _Tp(2)) * (__b + _Tp(2))
                     / (_Tp(2) * (__c + _Tp(1)));

      _Tp __F = _Tp(1);

      _Tp __Bnm3 = _Tp(1);
      _Tp __Bnm2 = _Tp(1) + __t1 * __x;
      _Tp __Bnm1 = _Tp(1) + __t2 * __x * (_Tp(1) + __t1 / _Tp(3) * __x);

      _Tp __Anm3 = _Tp(1);
      _Tp __Anm2 = __Bnm2 - __t0 * __x;
      _Tp __Anm1 = __Bnm1 - __t0 * (_Tp(1) + __t2 * __x) * __x
                 + __t0 * __t1 * (__c / (__c + _Tp(1))) * __x * __x;

      int __n = 3;
      while (1)
        {
          const _Tp __npam1 = _Tp(__n - 1) + __a;
          const _Tp __npbm1 = _Tp(__n - 1) + __b;
          const _Tp __npcm1 = _Tp(__n - 1) + __c;
          const _Tp __npam2 = _Tp(__n - 2) + __a;
          const _Tp __npbm2 = _Tp(__n - 2) + __b;
          const _Tp __npcm2 = _Tp(__n - 2) + __c;
          const _Tp __tnm1  = _Tp(2 * __n - 1);
          const _Tp __tnm3  = _Tp(2 * __n - 3);
          const _Tp __tnm5  = _Tp(2 * __n - 5);
          const _Tp __n2 = __n * __n;
          const _Tp __F1 = (_Tp(3) * __n2 + (__a + __b - _Tp(6)) * __n
                         + _Tp(2) - __a * __b - _Tp(2) * (__a + __b))
                         / (_Tp(2) * __tnm3 * __npcm1);
          const _Tp __F2 = -(_Tp(3) * __n2 - (__a + __b + _Tp(6)) * __n
                         + _Tp(2) - __a * __b) * __npam1 * __npbm1
                         / (_Tp(4) * __tnm1 * __tnm3 * __npcm2 * __npcm1);
          const _Tp __F3 = (__npam2 * __npam1 * __npbm2 * __npbm1
                         * (_Tp(__n - 2) - __a) * (_Tp(__n - 2) - __b))
                         / (_Tp(8) * __tnm3 * __tnm3 * __tnm5
                         * (_Tp(__n - 3) + __c) * __npcm2 * __npcm1);
          const _Tp __E  = -__npam1 * __npbm1 * (_Tp(__n - 1) - __c)
                         / (_Tp(2) * __tnm3 * __npcm2 * __npcm1);

          _Tp __An = (_Tp(1) + __F1 * __x) * __Anm1
                   + (__E + __F2 * __x) * __x * __Anm2 + __F3 * __x3 * __Anm3;
          _Tp __Bn = (_Tp(1) + __F1 * __x) * __Bnm1
                   + (__E + __F2 * __x) * __x * __Bnm2 + __F3 * __x3 * __Bnm3;
          const _Tp __r = __An / __Bn;

          const _Tp __prec = std::abs((__F - __r) / __F);
          __F = __r;

          if (__prec < __eps || __n > __nmax)
            break;

          if (std::abs(__An) > __big || std::abs(__Bn) > __big)
            {
              __An   /= __big;
              __Bn   /= __big;
              __Anm1 /= __big;
              __Bnm1 /= __big;
              __Anm2 /= __big;
              __Bnm2 /= __big;
              __Anm3 /= __big;
              __Bnm3 /= __big;
            }
          else if (std::abs(__An) < _Tp(1) / __big
                || std::abs(__Bn) < _Tp(1) / __big)
            {
              __An   *= __big;
              __Bn   *= __big;
              __Anm1 *= __big;
              __Bnm1 *= __big;
              __Anm2 *= __big;
              __Bnm2 *= __big;
              __Anm3 *= __big;
              __Bnm3 *= __big;
            }

          ++__n;
          __Bnm3 = __Bnm2;
          __Bnm2 = __Bnm1;
          __Bnm1 = __Bn;
          __Anm3 = __Anm2;
          __Anm2 = __Anm1;
          __Anm1 = __An;
        }

      if (__n >= __nmax) {return (std::numeric_limits<_Tp>::quiet_NaN());}

      return __F;
    }
    
    /**
     *   @brief  Return the confluent hypogeometric function
     *           @f$ _1F_1(a;c;x) @f$.
     * 
     *   @todo  Handle b == nonpositive integer blowup - return NaN.
     *
     *   @param  __a  The @a numerator parameter.
     *   @param  __c  The @a denominator parameter.
     *   @param  __x  The argument of the confluent hypergeometric function.
     *   @return  The confluent hypergeometric function.
     */
    template<typename _Tp>
    __ATTR_ALWAYS_INLINE__
    static inline 
    _Tp
    conf_hyperg(_Tp __a, _Tp __c, _Tp __x)
    {
      const _Tp __c_nint = static_cast<int>(__c + _Tp(0.5L));
      if (__isnan(__a) || __isnan(__c) || __isnan(__x))
        return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__c_nint == __c && __c_nint <= 0)
        return std::numeric_limits<_Tp>::infinity();
      else if (__a == _Tp(0))
        return _Tp(1);
      else if (__c == __a)
        return std::exp(__x);
      else if (__x < _Tp(0))
        return conf_hyperg_luke(__a, __c, __x);
      else
        return conf_hyperg_series(__a, __c, __x);
    }

    /**
     *   @brief Return \f$log(\Gamma(x))\f$ by the Lanczos method.
     *          This method dominates all others on the positive axis I think.
     *
     *   @param __x The argument of the log of the gamma function.
     *   @return  The logarithm of the gamma function.
     */
    template<typename _Tp>
    __ATTR_ALWAYS_INLINE__
    static inline 
    _Tp
    log_gamma_lanczos(_Tp __x)
    {
      constexpr _Tp __xm1 = __x - _Tp(1);

      constexpr _Tp __lanczos_cheb_7[9] = {
       _Tp( 0.99999999999980993227684700473478L),
       _Tp( 676.520368121885098567009190444019L),
       _Tp(-1259.13921672240287047156078755283L),
       _Tp( 771.3234287776530788486528258894L),
       _Tp(-176.61502916214059906584551354L),
       _Tp( 12.507343278686904814458936853L),
       _Tp(-0.13857109526572011689554707L),
       _Tp( 9.984369578019570859563e-6L),
       _Tp( 1.50563273514931155834e-7L)
      };

      static const _Tp __LOGROOT2PI
          = _Tp(0.9189385332046727417803297364056176L);

      _Tp __sum = __lanczos_cheb_7[0];
      for(unsigned int __k = 1; __k < 9; ++__k)
        __sum += __lanczos_cheb_7[__k] / (__xm1 + __k);

      const _Tp __term1 = (__xm1 + _Tp(0.5L))
                        * std::log((__xm1 + _Tp(7.5L))
                       / static_cast<_Tp>(2.7182818284590452353602874713526625L));
      const _Tp __term2 = __LOGROOT2PI + std::log(__sum);
      const _Tp __result = __term1 + (__term2 - _Tp(7));

      return __result;
    }

     /**
     *   @brief Return the sign of \f$ \Gamma(x) \f$.
     *          At nonpositive integers zero is returned.
     *
     *   @param __x The argument of the gamma function.
     *   @return  The sign of the gamma function.
     */
    template<typename _Tp>
    __ATTR_ALWAYS_INLINE__
    static inline 
    _Tp
    log_gamma_sign(_Tp __x)
    {
      if (__x > _Tp(0))
        return _Tp(1);
      else
        {
          const _Tp __sin_fact
                  = std::sin(static_cast<_Tp>(3.1415926535897932384626433832795029L) * __x);
          if (__sin_fact > _Tp(0))
            return (1);
          else if (__sin_fact < _Tp(0))
            return -_Tp(1);
          else
            return _Tp(0);
        }
    }


    /**
     *   @brief Return \f$ log(|\Gamma(x)|) \f$.
     *          This will return values even for \f$ x < 0 \f$.
     *          To recover the sign of \f$ \Gamma(x) \f$ for
     *          any argument use @a __log_gamma_sign.
     *
     *   @param __x The argument of the log of the gamma function.
     *   @return  The logarithm of the gamma function.
     */
    template<typename _Tp>
    __ATTR_ALWAYS_INLINE__
    static inline 
    _Tp
    log_gamma(_Tp __x)
    {
      if (__x > _Tp(0.5L))
        return __log_gamma_lanczos(__x);
      else
        {
          const _Tp __sin_fact
                 = std::abs(std::sin(static_cast<_Tp>(3.1415926535897932384626433832795029L) * __x));
          if (__sin_fact == _Tp(0)) {return (std::numeric_limits<_Tp>::quiet_NaN());}
             
          return __numeric_constants<_Tp>::__lnpi()
                     - std::log(__sin_fact)
                     - log_gamma_lanczos(_Tp(1) - __x);
        }
    }

     /**
     *   @brief Return the hypogeometric function @f$ _2F_1(a,b;c;x) @f$.
     *
     *   The hypogeometric function is defined by
     *   @f[
     *     _2F_1(a,b;c;x) = \frac{\Gamma(c)}{\Gamma(a)\Gamma(b)}
     *                      \sum_{n=0}^{\infty}
     *                      \frac{\Gamma(a+n)\Gamma(b+n)}{\Gamma(c+n)}
     *                      \frac{x^n}{n!}
     *   @f]
     *
     *   @param  __a  The first @a numerator parameter.
     *   @param  __a  The second @a numerator parameter.
     *   @param  __c  The @a denominator parameter.
     *   @param  __x  The argument of the confluent hypergeometric function.
     *   @return  The confluent hypergeometric function.
     */
    template<typename _Tp>
    __ATTR_ALWAYS_INLINE__ 
    static inline 
    _Tp
    hyperg(_Tp __a, _Tp __b, _Tp __c, _Tp __x)
    {
      const _Tp __a_nint = static_cast<int>(__a + _Tp(0.5L));
      const _Tp __b_nint = static_cast<int>(__b + _Tp(0.5L));
      const _Tp __c_nint = static_cast<int>(__c + _Tp(0.5L));
      const _Tp __toler = _Tp(1000) * std::numeric_limits<_Tp>::epsilon();
      if (std::abs(__x) >= _Tp(1)) { return (std::numeric_limits<_Tp>::quiet_NaN());}
      else if (__isnan(__a) || __isnan(__b)
            || __isnan(__c) || __isnan(__x))
        return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__c_nint == __c && __c_nint <= _Tp(0))
        return std::numeric_limits<_Tp>::infinity();
      else if (std::abs(__c - __b) < __toler || std::abs(__c - __a) < __toler)
        return std::pow(_Tp(1) - __x, __c - __a - __b);
      else if (__a >= _Tp(0) && __b >= _Tp(0) && __c >= _Tp(0)
            && __x >= _Tp(0) && __x < _Tp(0.995L))
        return hyperg_series(__a, __b, __c, __x);
      else if (std::abs(__a) < _Tp(10) && std::abs(__b) < _Tp(10))
        {
          //  For integer a and b the hypergeometric function is a
          //  finite polynomial.
          if (__a < _Tp(0)  &&  std::abs(__a - __a_nint) < __toler)
            return hyperg_series(__a_nint, __b, __c, __x);
          else if (__b < _Tp(0)  &&  std::abs(__b - __b_nint) < __toler)
            return hyperg_series(__a, __b_nint, __c, __x);
          else if (__x < -_Tp(0.25L))
            return hyperg_luke(__a, __b, __c, __x);
          else if (__x < _Tp(0.5L))
            return hyperg_series(__a, __b, __c, __x);
          else if (std::abs(__c) > _Tp(10))
              return hyperg_series(__a, __b, __c, __x);
        }
      else
        return hyperg_luke(__a, __b, __c, __x);
    }


}

}

#endif /*__GMS_HYPERGEOMETRIC_FUNC_H__*/