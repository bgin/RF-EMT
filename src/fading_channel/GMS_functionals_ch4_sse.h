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

#ifndef __GMS_FUNCTIONALS_CH4_SSE_H__
#define __GMS_FUNCTIONALS_CH4_SSE_H__ 1108020260108

#include <cstdint>
#include "GMS_config.h"


namespace file_info 
{

     static const unsigned int GMS_FUNCTIONALS_CH4_SSE_MAJOR = 1;
     static const unsigned int GMS_FUNCTIONALS_CH4_SSE_MINOR = 1;
     static const unsigned int GMS_FUNCTIONALS_CH4_SSE_MICRO = 0;
     static const unsigned int GMS_FUNCTIONALS_CH4_SSE_FULLVER =
       1000U*GMS_FUNCTIONALS_CH4_SSE_MAJOR+100U*GMS_FUNCTIONALS_CH4_SSE_MINOR+
       10U*GMS_FUNCTIONALS_CH4_SSE_MICRO;
     static const char GMS_FUNCTIONALS_CH4_SSE_CREATION_DATE[] = "11-08-2026 01:07AM +00200 (TUE 11 AUG 2026 GMT+2)";
     static const char GMS_FUNCTIONALS_CH4_SSE_BUILD_DATE[]    = __DATE__; 
     static const char GMS_FUNCTIONALS_CH4_SSE_BUILD_TIME[]    = __TIME__;
     static const char GMS_FUNCTIONALS_CH4_SSE_SYNOPSIS[]      = "Numerical integration (computing of Functionals) of the integrands (ch.4) used for the fading channel PDF/CDF and BER,BEP,SER,SEP quadrature computation(QUADPACK).\
	                                                          Based on the M.K Simon, M.S. Alouini: Digital Communication over Fading Channels 1st ed\
															  ISBN-13 978-0471317791";

}


namespace gms 
{

namespace fading_channel
{

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t 
integrate_4_1_gauss_Q_func_sse(double * __restrict__ __ATTR_ALIGN__(16), /* in_buf*/
                               double * __restrict__ __ATTR_ALIGN__(16), /* out_buf*/
                               double * __restrict__, /* functional values*/
                               const std::int32_t, /* n_func_args*/
                               const std::int32_t, /* n_gauss_q_vals*/
                               const double, /* lo_x*/
                               const double, /* hi_x*/
                               const double);/* hi_y*/



} // fading_channel

} // gms

#endif /*__GMS_FUNCTIONALS_CH4_SSE_H__*/