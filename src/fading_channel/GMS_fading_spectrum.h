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
#ifndef __GMS_FADING_SPECTRUM_H__
#define __GMS_FADING_SPECTRUM_H__ 280720260558


#include <cstdint>
#include <cmath> // bessel J0
#include "GMS_config.h"
#if (USE_OPENMP) == 1
#include <omp.h>
#endif
#include "GMS_cephes_sin_cos.h"


namespace file_info 
{

     static const unsigned int GMS_FADING_SPECTRUM_MAJOR = 1;
     static const unsigned int GMS_FADING_SPECTRUM_MINOR = 1;
     static const unsigned int GMS_FADING_SPECTRUM_MICRO = 0;
     static const unsigned int GMS_FADING_SPECTRUM_FULLVER =
       1000U*GMS_FADING_SPECTRUM_MAJOR+100U*GMS_FADING_SPECTRUM_MINOR+
       10U*GMS_FADING_SPECTRUM_MICRO;
     static const char GMS_FADING_SPECTRUM_CREATION_DATE[] = "28-07-2026 05:58PM +00200 (MON 28 JUL 2026 GMT+2)";
     static const char GMS_FADING_SPECTRUM_BUILD_DATE[]    = __DATE__; 
     static const char GMS_FADING_SPECTRUM_BUILD_TIME[]    = __TIME__;
     static const char GMS_FADING_SPECTRUM_SYNOPSIS[]      = "Fading [Radio]Channel analytical autocorrelation and PSD derived filters.";

}

namespace gms
{

namespace fading_channel
{

 
/*
    OpeMP auto-vectorization will probably preclude the optmization
	of cephes scalar code, hence the usage of cmath functions will
	be needed for the OpenMP path
*/
/* Frequency selective fading*/

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif 
__ATTR_ALWAYS_INLINE__
static inline 
float 
rectangular_autocorr(const float f_d, // doppler frequency spread
                     const float T_s) // symbol time-period
{
    constexpr float C6283185307179586476925{6.283185307179586476925f};
	const float sinc_arg{C6283185307179586476925*f_d*T_s};
#if (USE_OPENMP) == 1
	const float sinc_val{std::sin(sinc_arg)};
#else 
    const float sinc_val{gms::math::ceph_sinf(sinc_arg)};
#endif 
	return (sinc_val/sinc_arg);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif 
__ATTR_ALWAYS_INLINE__
static inline 
float
rect_normalized_psd(const float f_d)
{
	return (1.0f/(f_d+f_d));
}
///////////////////////////////////////////////////////////
#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif 
__ATTR_ALWAYS_INLINE__
static inline 
float
gaussian_autocorr(const float f_d,
                  const float T_s)
{
	constexpr float C3141592653589793238463{3.141592653589793238463f};
	const float exp_arg{C3141592653589793238463*f_d*T_s};
#if (USE_OPENMP) == 1
	const float exp_val{std::exp(-(exp_arg*exp_arg))};
#else 
    const float exp_val{gms::math::ceph_expf(-(exp_arg*exp_arg))};
#endif
    return (exp_val);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif 
__ATTR_ALWAYS_INLINE__
static inline 
float
gaussian_normalized_psd(const float f,
                        const float f_d)
{
	constexpr float C1772453850905516027298{1.772453850905516027298f};
	const float sec_factor{1.0f/(C1772453850905516027298*f_d)};
	const float exp_arg{f/f_d};
#if (USE_OPENMP) == 1
	const float exp_val{std::exp(-(exp_arg*exp_arg))};
#else 
    const float exp_val{gms::math::ceph_expf(-(exp_arg*exp_arg))};
#endif 
    return (exp_val*sec_factor);
}
//////////////////////////////////////////////////////////////////////////
#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif 
__ATTR_ALWAYS_INLINE__
static inline 
float
land_mobile_autocorr(const float f_d,
                     const float T_s)
{
	constexpr float C6283185307179586476925{6.283185307179586476925f};
    const float bes_j0_arg{C6283185307179586476925*f_d*T_s};
	const float bes_j0_val{std::cyl_bessel_j(0,bes_j0_arg)};
	return (bes_j0_val);
}

/* Precondition: std::abs(f)<=f_d */
#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif 
__ATTR_ALWAYS_INLINE__
static inline 
float
land_mobile_normalized_psd(const float f,
                           const float f_d) 
{
	constexpr float C9869604401089358618834{9.869604401089358618834f};
	const float f_sqr{f*f};
	const float f_d_sqr{f_d*f_d};
#if (USE_OPENMP) == 1
	const float inv_sqr{C9869604401089358618834*std::sqrt(f_sqr-f_d_sqr)};
#else 
    const float inv_sqr{C9869604401089358618834*gms::math::sqrtf(f_sqr-f_d_sqr)};
#endif 
    return (1.0f/inv_sqr);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif 
__ATTR_ALWAYS_INLINE__
static inline 
float
butterworth_1st_order(const float f_d,
                      const float T_s)
{
	constexpr float CN6283185307179586476925{-6.283185307179586476925f};
	const float exp_arg{CN6283185307179586476925/(f_d*T_s)};
#if (USE_OPENMP) == 1
	const float exp_val{std::exp(exp_arg)};
#else 
    const float exp_val{gms::math::ceph_expf(exp_arg)};
#endif 
    return (exp_val);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif 
__ATTR_ALWAYS_INLINE__
static inline 
float
butterworth_1st_order_psd(const float f,
                          const float f_d) 
{
	constexpr float C3141592653589793238463{3.141592653589793238463f};
	const float sec_term{1.0f+(f/f_d)};
	const float first_term{C3141592653589793238463*f};
	const float sqr_term{sec_term*sec_term};
	return (1.0f/(first_term*sqr_term));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif 
__ATTR_ALWAYS_INLINE__
static inline 
float
butterworth_2nd_order(const float f_d,
                      const float T_s) 
{
	constexpr float C3141592653589793238463{3.141592653589793238463f};
	constexpr float inv_sqrt2{0.707106781186547524401f};
	const float funcs_arg{C3141592653589793238463*f_d*T_s*inv_sqrt2};
#if (USE_OPENMP) == 1
    const float funcs_val{std::exp(-funcs_arg)*(std::cos(funcs_arg)+std::sin(funcs_arg))};
#else 
    const float exp_val{gms::math::ceph_expf(-funcs_arg)};
	const float cos_sin_val{gms::math::ceph_cosf(funcs_arg)+gms::math::sinf(funcs_arg)};
	const float funcs_val{exp_val*cos_sin_val};
#endif 
    return (funcs_val);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif 
__ATTR_ALWAYS_INLINE__
static inline 
float
butterworth_2nd_order_psd(const float f,
                          const float f_d)
{
	const float ratio{f/f_d};
	const float rat_sqr2{ratio*ratio};
	const float rat_pow4{rat_sqr2*rat_sqr2};
	const float term{1.0f+16.0f*rat_pow4};
	return (1.0f/term);
}

} // fading channel

} // gms

#endif /*__GMS_FADING_SPECTRUM_H__*/