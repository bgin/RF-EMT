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

#ifndef __GMS_INTEGRANDS_HELPERS_H__
#define __GMS_INTEGRANDS_HELPERS_H__ 130820260553

#include <cstdint>
#include "GMS_config.h"


namespace file_info 
{

     static const unsigned int GMS_INTEGRANDS_HELPERS_MAJOR = 1;
     static const unsigned int GMS_INTEGRANDS_HELPERS_MINOR = 1;
     static const unsigned int GMS_INTEGRANDS_HELPERS_MICRO = 0;
     static const unsigned int GMS_INTEGRANDS_HELPERS_FULLVER =
       1000U*GMS_INTEGRANDS_HELPERS_MAJOR+100U*GMS_INTEGRANDS_HELPERS_MINOR+
       10U*GMS_INTEGRANDS_HELPERS_MICRO;
     static const char GMS_INTEGRANDS_HELPERS_CREATION_DATE[] = "13-08-2026 06:05AM +00200 (THR 13 AUG 2026 GMT+2)";
     static const char GMS_INTEGRANDS_HELPERS_BUILD_DATE[]    = __DATE__; 
     static const char GMS_INTEGRANDS_HELPERS_BUILD_TIME[]    = __TIME__;
     static const char GMS_INTEGRANDS_HELPERS_SYNOPSIS[]      = "Integrands of the functional form used for fading channel and BER,BEP,SER,SEP quadrature computation(QUADPACK)\
                                                                 and by the polynomial fit of the tabulated data.\
                                                                 This file contains a helper functions for the assesment of inverse exponential argument range.\
	                                                             Based on the M.K Simon, M.S. Alouini: Digital Communication over Fading Channels 1st ed\
															     ISBN-13 978-0471317791";

}

#include <cstdint>
#include "GMS_config.h"


namespace gms
{

namespace fading_channel
{


/*
  Chapter 4, formula: 4.1
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t 
integrand_4_1_gauss_Q_arg_checker(double * __restrict__,
                                  const std::int32_t,
                                  const std::int32_t,
                                  const double,
                                  const double,
                                  const double);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t
integrand_4_2_gauss_Q_arg_checker(double * __restrict__,
                                  const std::int32_t,
                                  const std::int32_t,
                                  const double,
                                  const double,
                                  const double,
                                  const double);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t
integrand_4_6_gauss_Q_sin_arg_checker(double * __restrict__,
                                      const std::int32_t,
                                      const std::int32_t,
                                      const double,
                                      const double,
                                      const double,
                                      const double,
                                      const double,
                                      const double,
                                      const double,
                                      const double);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t
integrand_4_6_gauss_Q_cos_arg_checker(double * __restrict__,
                                      const std::int32_t,
                                      const std::int32_t,
                                      const double,
                                      const double,
                                      const double,
                                      const double,
                                      const double,
                                      const double,
                                      const double,
                                      const double);



} //math

} // gms

#endif /*__GMS_INTEGRANDS_HELPERS_H__*/