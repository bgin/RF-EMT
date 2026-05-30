/* Copyright (C) J. Racine, 1995-2001 
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

/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef __GMS_PDF_KERNELS_SSE_H__
#define __GMS_PDF_KERNELS_SSE_H__ 230520260617

namespace file_info 
{

     const unsigned int GMS_PDF_KERNELS_SSE_MAJOR = 1;
     const unsigned int GMS_PDF_KERNELS_SSE_MINOR = 0;
     const unsigned int GMS_PDF_KERNELS_SSE_MICRO = 0;
     const unsigned int GMS_PDF_KERNELS_SSE_FULLVER =
       1000U*GMS_PDF_KERNELS_SSE_MAJOR+100U*GMS_PDF_KERNELS_SSE_MINOR+
       10U*GMS_PDF_KERNELS_SSE_MICRO;
     const char * const GMS_PDF_KERNELS_SSE_CREATION_DATE = "23-05-2026 08:31 +00200 (SAT 23 MAY 2026 GMT+2)";
     const char * const GMS_PDF_KERNELS_SSE_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_PDF_KERNELS_SSE_SYNOPSIS      = "The Density estimation kernel functions (SSE version).Original implementation by J. Raccine modified by \
                                                              Bernard Gingold.";
}

#include <cstdint>
#include <cmath>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_simd_utils.h"

#if !defined(PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE)
#define PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE 1
#endif

namespace gms
{

namespace math 
{

namespace np_standalone_funcs
{

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_gauss2_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
     __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = {
                                                 0.398942280401432677939946059934f,
                                                 0.398942280401432677939946059934f,
                                                 0.398942280401432677939946059934f,
                                                 0.398942280401432677939946059934f,
                                                 -0.5f,-0.5f,-0.5f,-0.5f,
                                                 0.0f,0.0f,0.0f,0.0f,
                                                 0.0f,0.0f,0.0f,0.0f};
    if constexpr(use_prefetching==true)
    {
         _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0); 
    }
    const __m128 vC0398942280401432677939946059934{_mm_load_ps(&prefetched_constants[0])};
    const __m128 vC05{_mm_load_ps(&prefetched_constants[4])};
#else
    const __m128 vC0398942280401432677939946059934{_mm_set1_ps(0.398942280401432677939946059934f)};
    const __m128 vC05{_mm_set1_ps(-0.5f)};
#endif 
    const __m128 z_to_pow2{_mm_mul_ps(z,z)};
    const __m128 vexpf{_mm_exp_ps(_mm_mul_ps(vC05,z_to_pow2))};
    return (_mm_mul_ps(vC0398942280401432677939946059934,vexpf));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_gauss4_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
     __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = {
                                                 0.398942280401432677939946059934f,
                                                 0.398942280401432677939946059934f,
                                                 0.398942280401432677939946059934f,
                                                 0.398942280401432677939946059934f,
                                                 -0.5f,-0.5f,-0.5f,-0.5f,
                                                 1.5f,1.5f,1.5f,1.5f,
                                                 0.5f,0.5f,0.5f,0.5f};
    if constexpr(use_prefetching==true)
    {
      _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0); 
    }
    const __m128 vC0398942280401432677939946059934{_mm_load_ps(&prefetched_constants[0])};
    const __m128 vnC05{_mm_load_ps(&prefetched_constants[4])};
    const __m128 vC15{_mm_load_ps(&prefetched_constants[8])};
    const __m128 vC05{_mm_load_ps(&prefetched_constants[12])};
#else
    const __m128 vC0398942280401432677939946059934{_mm_set1_ps(0.398942280401432677939946059934f)};
    const __m128 vnC05{_mm_set1_ps(-0.5f)};
    const __m128 vC15{_mm_set1_ps(1.5f)};
    const __m128 vC05{_mm_set1_ps(0.5f)};
#endif 
    const __m128 z_to_pow2{_mm_mul_ps(z,z)};
    const __m128 vexpf_arg{_mm_mul_ps(vnC05,z_to_pow2)};
    const __m128 vexpf{_mm_exp_ps(vexpf_arg)};
    const __m128 vleft_term{_mm_mul_ps(vC0398942280401432677939946059934,_mm_sub_ps(vC15,_mm_mul_ps(vC05,z_to_pow2)))};
    return (_mm_mul_ps(vleft_term,vexpf));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_gauss6_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
     __ATTR_ALIGN__(16) constexpr float prefetched_constants[32] = {
                                                 0.398942280401432677939946059934f,
                                                 0.398942280401432677939946059934f,
                                                 0.398942280401432677939946059934f,
                                                 0.398942280401432677939946059934f,
                                                 0.125f,0.125f,0.125f,0.125f,
                                                 1.25f,1.25f,1.25f,1.25f,
                                                 1.875f,1.875f,1.875f,1.875f
                                                 -0.5f,-0.5f,-0.5f,-0.5f,
                                                 0.0f,0.0f,0.0f,0.0f,
                                                 0.0f,0.0f,0.0f,0.0f,
                                                 0.0f,0.0f,0.0f,0.f};
    if constexpr(use_prefetching==true)
    {
       _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0); 
       _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
    }
    const __m128 vC0398942280401432677939946059934{_mm_load_ps(&prefetched_constants[0])};
    const __m128 vnC0125{_mm_load_ps(&prefetched_constants[4])};
    const __m128 vC125{_mm_load_ps(&prefetched_constants[8])};
    const __m128 vC1875{_mm_load_ps(&prefetched_constants[12])};
    const __m128 vnC05{_mm_load_ps(&prefetched_constants[16])};
#else
    const __m128 vC0398942280401432677939946059934{_mm_set1_ps(0.398942280401432677939946059934f)};
    const __m128 vnC0125{_mm_set1_ps(0.125f)};
    const __m128 vC125{_mm_set1_ps(1.25f)};
    const __m128 vC1875{_mm_set1_ps(1.875f)};
    const __m128 vn05{_mm_set1_ps(-0.5f)};
#endif 
    const __m128 z_to_pow2{_mm_mul_ps(z,z)};
    const __m128 term1{_mm_fmsub_ps(z_to_pow2,vnC0125,vC125)};
    const __m128 vexpf{_mm_exp_ps(_mm_mul_ps(vnC05,z_to_pow2))};
    return (_mm_mul_ps(vC0398942280401432677939946059934,_mm_mul_ps(_mm_fmadd_ps(term1,z_to_pow2,vC1875),vexpf)));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_gauss8_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
     __ATTR_ALIGN__(16) constexpr float prefetched_constants[32] = {
                                                 0.398942280401432677939946059934f,
                                                 0.398942280401432677939946059934f,
                                                 0.398942280401432677939946059934f,
                                                 0.398942280401432677939946059934f,
                                                 0.4375f,0.4375f,0.4375f,0.4375f,
                                                 0.02083333333f,0.02083333333f,
                                                 0.02083333333f,0.02083333333f,
                                                 -2.1875f,-2.1875f,-2.1875f,-2.1875f,
                                                 2.1875f,2.1875f,2.1875f,2.1875f,
                                                 -0.5f,-0.5f,-0.5f,-0.5f,
                                                 0.0f,0.0f,0.0f,0.0f,
                                                 0.0f,0.0f,0.0f,0.0f};
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0); 
        _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
    }
    const __m128 vC0398942280401432677939946059934{_mm_load_ps(&prefetched_constants[0])};
    const __m128 vC4375{_mm_load_ps(&prefetched_constants[4])};
    const __m128 vC02083333333{_mm_load_ps(&prefetched_constants[8])};
    const __m128 vnC21875{_mm_load_ps(&prefetched_constants[12])};
    const __m128 vC21875{_mm_load_ps(&prefetched_constants[16])};
    const __m128 vnC05{_mm_load_ps(&prefetched_constants[20])};
#else 
    const __m128 vC0398942280401432677939946059934{_mm_set1_ps(0.398942280401432677939946059934f)};
    const __m128 vC4375{_mm_set1_ps(0.4375f)};
    const __m128 vC0208333333{_mm_set1_ps(0.02083333333f)};
    const __m128 vnC21875{_mm_set1_ps(-2.1875f)};
    const __m128 vC21875{_mm_set1_ps(2.1875f)};
    const __m128 vnC05{_mm_set1_ps(-0.5f)};
#endif 
    const __m128 z_to_pow2{_mm_mul_ps(z,z)};
    const __m128 term1{_mm_sub_ps(vC4375,_mm_mul_ps(z_to_pow2,vC02083333333))};
    const __m128 vexpf{_mm_exp_ps(_mm_mul_ps(vnC05,z_to_pow2))};
    const __m128 term2{_mm_fmadd_ps(_mm_fmadd_ps(term1,z_to_pow2,vnC21875),z_to_pow2,vC21875)};
    return (_mm_mul_ps(vC0398942280401432677939946059934,_mm_mul_ps(term2,vexpf)));
}

/*
   For Epanechnikov kernels the input check i.e. z^2 < 5.0 was removed 
   in order to remove the branch-misprediction.
*/

/*
    The input range is (z*z < 5.0)?result:0.0
*/

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_epanechnikov2_sse_ps(const __m128 z) 
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = {
                                                0.33541019662496845446f,
                                                0.33541019662496845446f,
                                                0.33541019662496845446f,
                                                0.33541019662496845446f,
                                                0.067082039324993690892f,
                                                0.067082039324993690892f,
                                                0.067082039324993690892f,
                                                0.067082039324993690892f,
                                                0.0f,0.0f,0.0f,0.0f,
                                                0.0f,0.0f,0.0f,0.0f};
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);       
    }                           
    const __m128 vC033541019662496845446{_mm_load_ps(&prefetched_constants[0])};
    const __m128 vC0067082039324993690892{_mm_load_ps(&prefetched_constants[4])};
#else 
    const __m128 vC033541019662496845446{_mm_set1_ps(0.33541019662496845446f)};
    const __m128 vC0067082039324993690892{_mm_set1_ps(0.067082039324993690892f)};
#endif 
    const __m128 z_to_pow2{_mm_mul_ps(z,z)};
    return (_mm_sub_ps(vC033541019662496845446,_mm_mul_ps(vC0067082039324993690892,z_to_pow2)));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_epanechnikov4_sse_ps(const __m128 z) 
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = {
                                    0.008385254916f,
                                    0.008385254916f,
                                    0.008385254916f,
                                    0.008385254916f,
                                    7.0f,7.0f,7.0f,7.0f,
                                    -15.0f,-15.0f,-15.0f,-15.0f,
                                    -5.0f,-5.0f,-5.0f,-5.0f};
    if constexpr(use_prefetching==true)
    {                                         
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0); 
    }
    const __m128 vC0008385254916{_mm_load_ps(&prefetched_constants[0])};
    const __m128 vC7{_mm_load_ps(&prefetched_constants[4])};
    const __m128 vnC15{_mm_load_ps(&prefetched_constants[8])};
    const __m128 vnC5{_mm_load_ps(&prefetched_constants[12])};
#else 
    const __m128 vC0008385254916{_mm_set1_ps(0.008385254916f)};
    const __m128 vC7{_mm_set1_ps(7.0f)};
    const __m128 vnC15{_mm_set1_ps(-15.0f)};
    const __m128 vnC5{_mm_set1_ps(-5.0f)};
#endif 
    const __m128 z_to_pow2{_mm_mul_ps(z,z)};
    const __m128 term1{_mm_fmadd_ps(z_to_pow2,vC7,vnC15)};
    const __m128 term2{_mm_add_ps(vnC5,z_to_pow2)};
    return (_mm_mul_ps(vC0008385254916,_mm_mul_ps(term1,term2)));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_epanechnikov6_sse_ps(const __m128 z) 
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[32] = {
                                        0.33541019662496845446f,
                                        0.33541019662496845446f,
                                        0.33541019662496845446f,
                                        0.33541019662496845446f,
                                        2.734375f,2.734375f,
                                        2.734375f,2.734375f,
                                        -3.28125f,-3.28125f,
                                        -3.28125f,-3.28125f,
                                        0.721875f,0.721875f,
                                        0.721875f,0.721875f,
                                        1.0f,1.0f,1.0f,1.0f,
                                        0.2f,0.2f,0.2f,0.2f,
                                        0.0f,0.0f,0.0f,0.0f,
                                        0.0f,0.0f,0.0f,0.0f};
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
    }
    const __m128 vC033541019662496845446{_mm_load_ps(&prefetched_constants[0])};
    const __m128 vC2734375{_mm_load_ps(&prefetched_constants[4])};
    const __m128 vnC328125{_mm_load_ps(&prefetched_constants[8])};
    const __m128 vC0721875{_mm_load_ps(&prefetched_constants[12])};
    const __m128 vC1{_mm_load_ps(&prefetched_constants[16])};
    const __m128 vC02{_mm_load_ps(&prefetched_constants[20])};
#else 
    const __m128 vC033541019662496845446{_mm_set1_ps(0.33541019662496845446f)};
    const __m128 vC2734375{_mm_set1_ps(2.734375f)};
    const __m128 vnC328125{_mm_set1_ps(-3.28125f)};
    const __m128 vC0721875{_mm_set1_ps(0.721875f)};
    const __m128 vC1{_mm_set1_ps(1.0f)};
    const __m128 vC02{_mm_set1_ps(0.2f)};
#endif 
    const __m128 z_to_pow2{_mm_mul_ps(z,z)};
    const __m128 term1{_mm_sub_ps(vC1,_mm_mul_ps(vC02,z_to_pow2))};
    const __m128 term2{_mm_fmadd_ps(_mm_fmadd_ps(z_to_pow2,vC0721875,vnC328125),z_to_pow2,vC2734375)};
    return (_mm_mul_ps(vC033541019662496845446,_mm_mul_ps(term1,term2)));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_epanechnikov8_sse_ps(const __m128 z) 
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[32] = {
                                     0.33541019662496845446f,
                                     0.33541019662496845446f,
                                     0.33541019662496845446f,
                                     0.33541019662496845446f,
                                     3.5888671875f,3.5888671875f,
                                     3.5888671875f,3.5888671875f,
                                     -7.8955078125f,-7.8955078125f,
                                     -7.8955078125f,-7.8955078125f,
                                     4.1056640625f,4.1056640625f,
                                     4.1056640625f,4.1056640625f,
                                     0.5865234375f,0.5865234375f,
                                     0.5865234375f,0.5865234375f,
                                     1.0f,1.0f,1.0f,1.0f,
                                     0.2f,0.2f,0.2f,0.2f,
                                     0.0f,0.0f,0.0f,0.0f};
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
    }
    const __m128 vC033541019662496845446{_mm_load_ps(&prefetched_constants[0])};
    const __m128 vC35888671875{_mm_load_ps(&prefetched_constants[4])};
    const __m128 vnC78955078125{_mm_load_ps(&prefetched_constants[8])};
    const __m128 vC41056640625{_mm_load_ps(&prefetched_constants[12])};
    const __m128 vC05865234375{_mm_load_ps(&prefetched_constants[16])};
    const __m128 vC1{_mm_load_ps(&prefetched_constants[20])};
    const __m128 vC02{_mm_load_ps(&prefetched_constants[24])};
#else 
    const __m128 vC033541019662496845446{_mm_set1_ps(0.33541019662496845446f)};
    const __m128 vC35888671875{_mm_set1_ps(3.5888671875f)};
    const __m128 vnC78955078125{_mm_set1_ps(-7.8955078125f)};
    const __m128 vC41056640625{_mm_set1_ps(4.1056640625f)};
    const __m128 vC05865234375{_mm_set1_ps(0.5865234375f)};
    const __m128 vC1{_mm_set1_ps(1.0f)};
    const __m128 vC02{_mm_set1_ps(0.2f)};
#endif 
    const __m128 z_to_pow2{_mm_mul_ps(z,z)};
    const __m128 term1{_mm_sub_ps(vC1,_mm_mul_ps(vC02,z_to_pow2))};
    const __m128 term2{_mm_sub_ps(vC41056640625,_mm_mul_ps(vC05865234375,z_to_pow2))};
    const __m128 term3{_mm_fmadd_ps(_mm_fmadd_ps(term2,z_to_pow2,vnC78955078125),z_to_pow2,vC35888671875)};
    return (_mm_mul_ps(vC033541019662496845446,_mm_mul_ps(term3,term1)));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_rect_sse_ps(const __m128 z) 
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = {  
                                        0.5f,0.5f,0.5f,0.5f,
                                        0.0f,0.0f,0.0f,0.0f,
                                        1.0f,1.0f,1.0f,1.0f,
                                        0.0f,0.0f,0.0f,0.0f};
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }
    const __m128 vC05{_mm_load_ps(&prefetched_constants[0])};
    const __m128 vC0{_mm_load_ps(&prefetched_constants[4])};
    const __m128 vC1{_mm_load_ps(&prefetched_constants[8])};
#else 
    const __m128 vC05{_mm_set1_ps(0.5f)};
    const __m128 vC0{_mm_setzero_ps()};
    const __m128 vC1{_mm_set1_ps(1.0f)};
#endif 
    const __m128 z_to_pow2{_mm_mul_ps(z,z)};
    const __mmask8 z_to_pow2_lt_1{_mm_cmp_ps_mask(z_to_pow2,vC1,_CMP_LT_OQ)};
    return (_mm_mask_blend_ps(z_to_pow2_lt_1,vC0,vC05));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_uaa_sse_ps(const bool   same_cat,
                     const __m128 vlambda,
                     const __m128 vc)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = { 
                                      1.0f,1.0f,1.0f,1.0f,
                                      0.0f,0.0f,0.0f,0.0f,
                                      0.0f,0.0f,0.0f,0.0f,
                                      0.0f,0.0f,0.0f,0.0f};
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }
    const __m128 vC1{_mm_load_ps(&prefetched_constants[0])};
#else 
    const __m128 vC1{_mm_set1_ps(1.0f)};
#endif
    return ((same_cat==true) ? _mm_sub_ps(vC1,vlambda) : _mm_div_ps(vC1,_mm_sub_ps(vc,vC1)));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_score_uaa_sse_ps(const bool same_cat,
                           const __m128 vc) 
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = { 
                                      1.0f,1.0f,1.0f,1.0f,
                                     -1.0f,-1.0f,-1.0f,-1.0f,
                                      0.0f,0.0f,0.0f,0.0f,
                                      0.0f,0.0f,0.0f,0.0f};
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }
    const __m128 vC1{_mm_load_ps(&prefetched_constants[0])};
    const __m128 vnC1{_mm_load_ps(&prefetched_constants[4])};
#else 
    const __m128 vC1{_mm_set1_ps(1.0f)};
    const __m128 vnC1{_mm_set1_ps(-1.0f)};
#endif
    return ((same_cat==true) ? vnC1 : _mm_div_ps(vC1,_mm_sub_ps(vc,vC1)));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_uli_racine_sse_ps(const bool same_cat,
                            const __m128 vlambda) 
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = { 
                                      1.0f,1.0f,1.0f,1.0f,
                                      0.0f,0.0f,0.0f,0.0f,
                                      0.0f,0.0f,0.0f,0.0f,
                                      0.0f,0.0f,0.0f,0.0f};
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }
    const __m128 vC1{_mm_load_ps(&prefetched_constants[0])};
#else 
    const __m128 vC1{_mm_set1_ps(1.0f)};
#endif
    return ((same_cat==true)?vC1:vlambda);
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_unli_racine_sse_ps(const bool same_cat,
                             const __m128 vlambda,
                             const __m128 vc) 
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = { 
                                      1.0f,1.0f,1.0f,1.0f,
                                      0.0f,0.0f,0.0f,0.0f,
                                      0.0f,0.0f,0.0f,0.0f,
                                      0.0f,0.0f,0.0f,0.0f};
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }
    const __m128 vC1{_mm_load_ps(&prefetched_constants[0])};
#else 
    const __m128 vC1{_mm_set1_ps(1.0f)};
#endif
    return ((same_cat==true) ? vC1 : _mm_div_ps(vC1,_mm_fmadd_ps(_mm_sub_ps(vc,vC1),vlambda,vC1)));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_score_unli_racine_sse_ps(const bool same_cat)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = { 
                                      1.0f,1.0f,1.0f,1.0f,
                                      0.0f,0.0f,0.0f,0.0f,
                                      0.0f,0.0f,0.0f,0.0f,
                                      0.0f,0.0f,0.0f,0.0f};
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }
    const __m128 vC1{_mm_load_ps(&prefetched_constants[0])};
#else 
    const __m128 vC1{_mm_set1_ps(1.0f)};
    
#endif
    return ((same_cat==true) ? _mm_setzero_ps() : vC1);
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_score_unli_racine_sse_ps(const bool same_cat,
                                   const __m128 vlambda,
                                   const __m128 vc)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = { 
                                      1.0f,1.0f,1.0f,1.0f,
                                      0.0f,0.0f,0.0f,0.0f,
                                      0.0f,0.0f,0.0f,0.0f,
                                      0.0f,0.0f,0.0f,0.0f};
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }
    const __m128 vC1{_mm_load_ps(&prefetched_constants[0])};
#else 
    const __m128 vC1{_mm_set1_ps(1.0f)};
#endif   
    const __m128 tmp1{_mm_sub_ps(vc,vC1)};
    const __m128 vinorm{_mm_div_ps(vC1,_mm_fmadd_ps(tmp1,vlambda,vC1))};
    const __m128 tmp2{_mm_sub_ps(vC1,vc)};
    const __m128 true_br{_mm_mul_ps(tmp2,_mm_mul_ps(vinorm,vinorm))};
    const __m128 false_br{_mm_mul_ps(vinorm,_mm_fmadd_ps(_mm_mul_ps(vlambda,tmp2),vinorm,vC1))};
    return ((same_cat==true)?true_br:false_br);
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_owang_van_ryzin_sse_ps(const __m128 vx,const __m128 vy,
                                 const __m128 vlambda)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = { 
                                      1.0f,1.0f,1.0f,1.0f,
                                      0.5f,0.5f,0.5f,0.5f,
                                      0.0f,0.0f,0.0f,0.0f,
                                      0.0f,0.0f,0.0f,0.0f};
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }
    const __m128 vC1{_mm_load_ps(&prefetched_constants[0])};
    const __m128 vC05{_mm_load_ps(&prefetched_constants[4])};
#else 
    const __m128 vC1{_mm_set1_ps(1.0f)};
    const __m128 vC05{_mm_set1_ps(0.5f)};
#endif      
    const __m128 exponent{_mm_floor_ps(gms::math::xmm4r4_abs(_mm_sub_ps(vx,vy)))};
    const __mmask8 vx_eq_vy{_mm_cmp_ps_mask(vx,vy,_CMP_EQ_OQ)};
    const __m128 vtmp{_mm_sub_ps(vC1,vlambda)};
    const __m128 vlambda_to_powxy{_mm_pow_ps(vlambda,exponent)};
    return (_mm_mask_blend_ps(vx_eq_vy,_mm_mul_ps(vlambda_to_powxy,_mm_mul_ps(vtmp,vC05)),vtmp));
}
                                 
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_score_owang_van_ryzin_sse_ps(const __m128 vx,const __m128 vy,
                                       const __m128 vlambda)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = { 
                                     -1.0f,-1.0f,-1.0f,-1.0f,
                                      0.5f,0.5f,0.5f,0.5f,
                                      2.0f,2.0f,2.0f,2.0f,
                                      0.0f,0.0f,0.0f,0.0f};
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }
    const __m128 vnC1{_mm_load_ps(&prefetched_constants[0])};
    const __m128 vC05{_mm_load_ps(&prefetched_constants[4])};
    const __m128 vC2{_mm_load_ps(&prefetched_constants[8])};
#else 
    const __m128 vnC1{_mm_set1_ps(-1.0f)};
    const __m128 vC05{_mm_set1_ps(0.5f)};
    const __m128 vC2{_mm_set1_ps(2.0f)}
#endif   
    const __m128   exponent{_mm_floor_ps(gms::math::xmm4r4_abs(_mm_sub_ps(vx,vy)))};
    const __m128   vlambda_to_xypow{_mm_mul_ps(vC05,_mm_pow_ps(vlambda,exponent))};
    const __mmask8 vx_eq_vy{_mm_cmp_ps_mask(vx,vy,_CMP_EQ_OQ)};
    const __m128   last_term{_mm_sub_ps(_mm_div_ps(exponent,vlambda),vC2)};
    return (_mm_mask_blend_ps(vx_eq_vy,_mm_mul_ps(vlambda_to_xypow,last_term),vnC1));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
__ATTR_ALWAYS_INLINE__
static inline
__m128 np_oli_racine_sse_ps(const __m128 vx,const __m128 vy,
                            const __m128 vlambda)
{
    const __m128   exponent{_mm_floor_ps(gms::math::xmm4r4_abs(_mm_sub_ps(vx,vy)))};
    return (_mm_pow_ps(vlambda,exponent));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_score_oli_racine_sse_ps(const __m128 vx,const __m128 vy,
                                  const __m128 vlambda)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = { 
                                      1.0f,1.0f,1.0f,1.0f,
                                      0.0f,0.0f,0.0f,0.0f,
                                      0.0f,0.0f,0.0f,0.0f,
                                      0.0f,0.0f,0.0f,0.0f};
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }
    const __m128 vC1{_mm_load_ps(&prefetched_constants[0])};
#else 
    const __m128 vC1{_mm_set1_ps(1.0f)};
#endif   
    const __m128 tmp1{gms::math::xmm4r4_abs(_mm_sub_ps(vx,vy))};
    const __m128 exponent{_mm_floor_ps(gms::math::xmm4r4_abs(_mm_sub_ps(vx,vy)))};
    return (_mm_mul_ps(tmp1,_mm_pow_ps(vlambda,_mm_sub_ps(exponent,vC1))));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_score_onli_racine_sse_ps(const __m128 vx,const __m128 vy,
                                  const __m128 vlambda)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = { 
                                      0.0f,0.0f,0.0f,0.0f,
                                      1.0f,1.0f,1.0f,1.0f,
                                      2.0f,2.0f,2.0f,2.0f,
                                     -2.0f,-2.0f,-2.0f,-2.0f};
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }
    const __m128 vC0{_mm_load_ps(&prefetched_constants[0])};
    const __m128 vC1{_mm_load_ps(&prefetched_constants[4])};
    const __m128 vC2{_mm_load_ps(&prefetched_constants[8])};
    const __m128 vCn2{_mm_load_ps(&prefetched_constants(12))};
#else 
    const __m128 vC0{_mm_setzero_ps()};
    const __m128 vC1{_mm_set1_ps(1.0f)};
    const __m128 vC2{_mm_set1_ps(2.0f)};
    const __m128 vCn2{_mm_set1_ps(-2.0f)};
#endif       
    const __m128   cxy{_mm_floor_ps(gms::math::xmm4r4_abs(_mm_sub_ps(vx,vy)))};
    const __mmask8 lambda_neq_0{_mm_cmp_ps_mask(vlambda,vC0,_CMP_NEQ_OQ)};
    const __m128   lambda_to_cxypow{_mm_pow_ps(vlambda,_mm_sub_ps(cxy,vC1))};
    const __mmask8 cxy_neq_0{_mm_cmp_ps_mask(cxy,vC0,_CMP_NEQ_OQ)};
    const __m128   term1{_mm_sub_ps(vC1,_mm_mul_ps(vlambda,vlambda))};
    const __m128   term2{_mm_fmsub_ps(cxy,term1,_mm_add_ps(vlambda,vlambda))};
    // Must be reimplemented by using kunpckbw or similiar intrinsic!!
    return (_mm_mask_blend_ps(_kor_mask8(cxy_neq_0,lambda_neq_0),vCn2,_mm_mul_ps(lambda_to_cxypow,term2)));
}

// the simple convolution kernels
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_econvol_gauss2_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = { 
                                      0.28209479177387814348f,
                                      0.28209479177387814348f,
                                      0.28209479177387814348f,
                                      0.28209479177387814348f,
                                     -0.25f,-0.25f,-0.25f,-0.25f,
                                      0.0f,0.0f,0.0f,0.0f,
                                      0.0f,0.0f,0.0f,0.0f};
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }
    const __m128 vC028209479177387814348{_mm_load_ps(&prefetched_constants[0])};
    const __m128 vCn025{_mm_load_ps(&prefetched_constants[4])};
#else 
    const __m128 vC028209479177387814348{_mm_set1_ps(0.28209479177387814348f)};
    const __m128 vCn025{_mm_set1_ps(-0.25f)};
#endif   
    const __m128 z_to_pow2{_mm_mul_ps(z,z)};   
    const __m128 vexpf{_mm_exp_ps(_mm_mul_ps(vCn025,z_to_pow2))};
    return (_mm_mul_ps(vC028209479177387814348,vexpf));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_econvol_gauss4_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = 
    {
                                         0.0044077311214668459918f,
                                         0.0044077311214668459918f,
                                         0.0044077311214668459918f,
                                         0.0044077311214668459918f,
                                        -0.25f,-0.25f,-0.25f,-0.25f,
                                         108.0f,108.0f,108.0f,108.0f,
                                         28.0f,28.0f,28.0f,28.0f
    };
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }
    const __m128 vC00044077311214668459918{_mm_load_ps(&prefetched_constants[0])};
    const __m128 vC025{_mm_load_ps(&prefetched_constants[4])};
    const __m128 vC108{_mm_load_ps(&prefetched_constants[8])};
    const __m128 vC28{_mm_load_ps(&prefetched_constants[12])};
#else 
    const __m128 vC00044077311214668459918{_mm_set1_ps(0.0044077311214668459918f)};
    const __m128 vC025{_mm_set1_ps(0.25f)};
    const __m128 vC108{_mm_set1_ps(108.0f)};
    const __m128 vC28{_mm_set1_ps(28.0f)};
#endif
    const __m128 z_to_pow2{_mm_mul_ps(z,z)}; 
    const __m128 tmp1{_mm_sub_ps(z_to_pow2,vC28)};
    const __m128 vexpf{_mm_exp_ps(_mm_mul_ps(vC025,z_to_pow2))};
    return (_mm_mul_ps(vC00044077311214668459918,_mm_mul_ps(vexpf,_mm_fmadd_ps(tmp1,z_to_pow2,vC108))));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_econvol_gauss6_sse_ps(const __m128 z)
{ 
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[32] = 
    {
                                      0.00001721769969f,
                                      0.00001721769969f,
                                      0.00001721769969f,
                                      0.00001721769969f,
                                     -0.25f,-0.25f,-0.25f,-0.25f,
                                      36240.0f,36240.0f,
                                      36240.0f,36240.0f,
                                     -19360.0f,-19360.0f
                                     -19360.0f,-19360.0f,
                                      2312.0f,2312.0f,
                                      2312.0f,2312.0f,
                                      -88.0f,-88.0f,-88.0f,-88.0f
                                      0.0f,0.0f,0.0f,0.0f,
                                      0.0f,0.0f,0.0f,0.0f
    };
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
    }
        const __m128 vC000001721769969{_mm_load_ps(&prefetched_constants[0])};
        const __m128 vCN025{_mm_load_ps(&prefetched_constants[4])};
        const __m128 vC3624{_mm_load_ps(&prefetched_constants[8])};
        const __m128 vCN19360{_mm_load_ps(&prefetched_constants[12])};
        const __m128 vC2312{_mm_load_ps(&prefetched_constants[16])};
        const __m128 vCN88{_mm_load_ps(&prefetched_constants[20])};
    
#else 
        const __m128 vC000001721769969{_mm_set1_ps(0.00001721769969f)};
        const __m128 vCN025{_mm_set1_ps(-0.25f)};
        const __m128 vC3624{_mm_set1_ps(3624.0f)};
        const __m128 vCN19360{_mm_set1_ps(19360.0f)};
        const __m128 vC2312{_mm_set1_ps(2312.0f)};
        const __m128 vCN88{_mm_set1_ps(-88.0f)};
#endif 
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 factor1{_mm_mul_ps(vCN025,z_to_pow2)};
        const __m128 factor2{_mm_add_ps(vCN88,z_to_pow2)};
        const __m128 term_exp{_mm_exp_ps(factor1)};
        const __m128 horner_series{_mm_fmadd_ps(
                                            _mm_fmadd_ps(
                                                     _mm_fmadd_ps(factor2,z_to_pow2,vC2312),
                                                                             z_to_pow2,vCN19360),
                                                                                      z_to_pow2,vC3624)};
        return (_mm_mul_ps(vC000001721769969,_mm_mul_ps(term_exp,horner_series)));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_econvol_gauss8_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[32] = 
    {
                                       0.2989183974E-7f,0.2989183974E-7,
                                       0.2989183974E-7f,0.2989183974E-7f,
                                      -0.25f,-0.25f,-0.25f,-0.25f,
                                       25018560.0f,25018560.0f,
                                       25018560.0f,25018560.0f,
                                      -20462400.0f,-20462400.0f,
                                      -20462400.0f,-20462400.0f,
                                       4202352.0f,4202352.0f,
                                       4202352.0f,4202352.0f,
                                      -331680.0f,-331680.0f,
                                      -331680.0f,-331680.0f,
                                       11604.0f,11604.0f,
                                       11604.0f,11604.0f,
                                      -180.0f,-180.0f,-180.0f,-180.0f
    };
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
    }
        const __m128 vC02989183974{_mm_load_ps(&prefetched_constants[0])};
        const __m128 vCN025{_mm_load_ps(&prefetched_constants[4])};
        const __m128 vC25018560{_mm_load_ps(&prefetched_constants[8])};
        const __m128 vCN20462400{_mm_load_ps(&prefetched_constants[12])};
        const __m128 vC4202352{_mm_load_ps(&prefetched_constants[16])};
        const __m128 vCN331680{_mm_load_ps(&prefetched_constants[20])};
        const __m128 vC11604{_mm_load_ps(&prefetched_constants[24])};
        const __m128 vCN180{_mm_load_ps(&prefetched_constants[28])};
#else 
        const __m128 vC02989183974{_mm_set1_ps(0.2989183974E-7f)};
        const __m128 vCN025{_mm_set1_ps(-0.25f)};
        const __m128 vC25018560{_mm_set1_ps(25018560.0f)};
        const __m128 vCN20462400{_mm_set1_ps(-20462400.0f)};
        const __m128 vC4202352{_mm_set1_ps(4202352.0f)};
        const __m128 vCN331680{_mm_set1_ps(-331680.0f)};
        const __m128 vC11604{_mm_set1_ps(11604.0f)};
        const __m128 vCN180{_mm_set1_ps(-180.0f)};
#endif 
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 factor1{_mm_mul_ps(vCN025,z_to_pow2)};
        const __m128 factor2{_mm_add_ps(vCN180,z_to_pow2)};
        const __m128 term_exp{_mm_exp_ps(factor1)};
        const __m128 horner_series{_mm_fmadd_ps(
                                        _mm_fmadd_ps(
                                            _mm_fmadd_ps(
                                                _mm_fmadd_ps(
                                                    _mm_fmadd_ps(factor2,z_to_pow2,vC11604),
                                                                            z_to_pow2,vCN331680),
                                                                                    z_to_pow2,vC4202352),
                                                                                            z_to_pow2,vCN20462400),
                                                                                                z_to_pow2,vC25018560)};
        return (_mm_mul_ps(vC02989183974,_mm_mul_ps(term_exp,horner_series)));
}

/*
   Accordingly to original implementation these kernel may be inaccurate.
   // These kernels are in error jracine 27 5 2009
   // need to edit both kernel.c and this file
   // note we require an additional test for negative (z) as well 
*/
/*
***WARNING***
   The check for z>=SQRT(20) was removed in order to minimize
   the if-conversion nesting level.
   Further optimized by splitting into z < 0  and z > 0
   This version computes the case: z < 0
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_econvol_epan2_z_lt_0_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[32] = 
    {   
                       5.579734404642339E-9f,5.579734404642339E-9f,
                       5.579734404642339E-9f,5.579734404642339E-9f,
                       26883.0f,26883.0f,26883.0f,26883.0f,
                       2688300.0f,2688300.0f,2688300.0f,2688300.0f,
                       12022443.0f,12022443.0f,12022443.0f,12022443.0f,
                       48089773.0f,48089773.0f,48089773.0f,48089773.0f,
                       0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
                       0.0f,0.0f,0.0f,0.0f
    };
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
    }
        const __m128 vC5579734404642339{_mm_load_ps(&prefetched_constants[0])};
        const __m128 vC26883{_mm_load_ps(&prefetched_constants[4])};
        const __m128 vC2688300{_mm_load_ps(&prefetched_constants[8])};
        const __m128 vC12022443{_mm_load_ps(&prefetched_constants[12])};
        const __m128 vC48089773{_mm_load_ps(&prefetched_constants[16])};
#else 
        const __m128 vC5579734404642339{_mm_set1_ps(5.579734404642339E-9f)};
        const __m128 vC26883{_mm_set1_ps(26883.0f)};
        const __m128 vC2688300{_mm_set1_ps(2688300.0f)};
        const __m128 vC12022443{_mm_set1_ps(12022443.0f)};
        const __m128 vC48089773{_mm_set1_ps(48089773.0f)};
#endif 
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 z_to_pow3{_mm_mul_ps(z_to_pow2,z)};
        const __m128 z_to_pow5{_mm_mul_ps(z_to_pow3,z_to_pow2)};
        const __m128 horner_series{_mm_fmsub_ps(vC26883,z_to_pow5,
                                            _mm_fmsub_ps(vC2688300,z_to_pow3,
                                                     _mm_fmadd_ps(vC12022443,z_to_pow2,vC48089773)))};
        return (_mm_mul_ps(vC5579734404642339,horner_series));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_econvol_epan2_z_gt_0_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[32] = 
    {   
                      -5.579734404642339E-9f,-5.579734404642339E-9f,
                      -5.579734404642339E-9f,-5.579734404642339E-9f,
                       26883.0f,26883.0f,26883.0f,26883.0f,
                       2688300.0f,2688300.0f,2688300.0f,2688300.0f,
                       12022443.0f,12022443.0f,12022443.0f,12022443.0f,
                       48089773.0f,48089773.0f,48089773.0f,48089773.0f,
                       0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
                       0.0f,0.0f,0.0f,0.0f
    };
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
    }
        const __m128 vCN5579734404642339{_mm_load_ps(&prefetched_constants[0])};
        const __m128 vC26883{_mm_load_ps(&prefetched_constants[4])};
        const __m128 vC2688300{_mm_load_ps(&prefetched_constants[8])};
        const __m128 vC12022443{_mm_load_ps(&prefetched_constants[12])};
        const __m128 vC48089773{_mm_load_ps(&prefetched_constants[16])};
#else 
        const __m128 vCN5579734404642339{_mm_set1_ps(-5.579734404642339E-9f)};
        const __m128 vC26883{_mm_set1_ps(26883.0f)};
        const __m128 vC2688300{_mm_set1_ps(2688300.0f)};
        const __m128 vC12022443{_mm_set1_ps(12022443.0f)};
        const __m128 vC48089773{_mm_set1_ps(48089773.0f)};
#endif 
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 z_to_pow3{_mm_mul_ps(z_to_pow2,z)};
        const __m128 z_to_pow5{_mm_mul_ps(z_to_pow3,z_to_pow2)};
        const __m128 horner_series{_mm_fmsub_ps(vC26883,z_to_pow5,
                                            _mm_fmadd_ps(vC2688300,z_to_pow3,
                                                     _mm_fmsub_ps(vC12022443,z_to_pow2,vC48089773)))};
        return (_mm_mul_ps(vCN5579734404642339,horner_series));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_econvol_epan4_z_lt_0_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[32] = 
    {  
                               3.756009615384615e-9f,3.756009615384615e-9f
                               3.756009615384615e-9f,3.756009615384615e-9f,
                               1456.0f,1456.0f,1456.0f,1456.0f,
                               124800.0f,124800.0f,124800.0f,124800.0f,
                               5491200.0f,5491200.0f,5491200.0f,5491200.0f,
                               15627432.0f,15627432.0f,15627432.0f,15627432.0f,
                               24960000.0f,24960000.0f,24960000.0f,24960000.0f,
                               111624513.0f,111624513.0f,111624513.0f,111624513.0f,
                               148832684.0f,148832684.0f,148832684.0f,148832684.0f
    };
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
    }
        const __m128 vC3756009615384615{_mm_load_ps(&prefetched_constants[0])};
        const __m128 vC1456{_mm_load_ps(&prefetched_constants[4])};
        const __m128 vC124800{_mm_load_ps(&prefetched_constants[8])};
        const __m128 vC5491200{_mm_load_ps(&prefetched_constants[12])};
        const __m128 vC15627432{_mm_load_ps(&prefetched_constants[16])};
        const __m128 vC24960000{_mm_load_ps(&prefetched_constants[20])};
        const __m128 vC111624513{_mm_load_ps(&prefetched_constants[24])};
        const __m128 vC148832684{_mm_load_ps(&prefetched_constants[28])};
#else 
        const __m128 vC3756009615384615{_mm_set1_ps(3.756009615384615e-9f)};
        const __m128 vC1456{_mm_set1_ps(1456.0f)};
        const __m128 vC124800{_mm_set1_ps(124800.0f)};
        const __m128 vC5491200{_mm_set1_ps(5491200.0f)};
        const __m128 vC15627432{_mm_set1_ps(15627432.0f)};
        const __m128 vC24960000{_mm_set1_ps(24960000.0f)};
        const __m128 vC111624513{_mm_set1_ps(111624513.0f)};
        const __m128 vC148832684{_mm_set1_ps(148832684.0f)};
#endif 
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 z_to_pow3{_mm_mul_ps(z_to_pow2,z)};
        const __m128 z_to_pow4{_mm_mul_ps(z_to_pow2,z_to_pow2)};
        const __m128 z_to_pow5{_mm_mul_ps(z_to_pow3,z_to_pow2)};
        const __m128 z_to_pow7{_mm_mul_ps(z_to_pow5,z_to_pow2)};
        const __m128 z_to_pow9{_mm_mul_ps(z_to_pow7,z_to_pow2)};
        const __m128 horner_series{_mm_fmsub_ps(vC1456,z_to_pow9,
                                      _mm_fmadd_ps(vC124800,z_to_pow7,
                                          _mm_fmadd_ps(vC5491200,z_to_pow5,
                                            _mm_fmadd_ps(vC15627432,z_to_pow4,
                                                _mm_fmsub_ps(vC24960000,z_to_pow3,
                                                    _mm_fmadd_ps(vC111624513,z_to_pow2,vC148832684))))))};
        return (_mm_mul_ps(vC3756009615384615,horner_series));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching,bool use_horner_scheme>
__m128 np_econvol_epan4_z_gt_0_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[32] = 
    {  
                              -3.756009615384615e-9f,-3.756009615384615e-9f
                              -3.756009615384615e-9f,-3.756009615384615e-9f,
                               1456.0f,1456.0f,1456.0f,1456.0f,
                               124800.0f,124800.0f,124800.0f,124800.0f,
                               5491200.0f,5491200.0f,5491200.0f,5491200.0f,
                               15627432.0f,15627432.0f,15627432.0f,15627432.0f,
                               24960000.0f,24960000.0f,24960000.0f,24960000.0f,
                               111624513.0f,111624513.0f,111624513.0f,111624513.0f,
                               148832684.0f,148832684.0f,148832684.0f,148832684.0f
    };
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
    }
        const __m128 vCN3756009615384615{_mm_load_ps(&prefetched_constants[0])};
        const __m128 vC1456{_mm_load_ps(&prefetched_constants[4])};
        const __m128 vC124800{_mm_load_ps(&prefetched_constants[8])};
        const __m128 vC5491200{_mm_load_ps(&prefetched_constants[12])};
        const __m128 vC15627432{_mm_load_ps(&prefetched_constants[16])};
        const __m128 vC24960000{_mm_load_ps(&prefetched_constants[20])};
        const __m128 vC111624513{_mm_load_ps(&prefetched_constants[24])};
        const __m128 vC148832684{_mm_load_ps(&prefetched_constants[28])};
#else 
        const __m128 vCN3756009615384615{_mm_set1_ps(-3.756009615384615e-9f)};
        const __m128 vC1456{_mm_set1_ps(1456.0f)};
        const __m128 vC124800{_mm_set1_ps(124800.0f)};
        const __m128 vC5491200{_mm_set1_ps(5491200.0f)};
        const __m128 vC15627432{_mm_set1_ps(15627432.0f)};
        const __m128 vC24960000{_mm_set1_ps(24960000.0f)};
        const __m128 vC111624513{_mm_set1_ps(111624513.0f)};
        const __m128 vC148832684{_mm_set1_ps(148832684.0f)};
#endif 
    if constexpr(use_horner_scheme==true)
    {
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 z_to_pow3{_mm_mul_ps(z_to_pow2,z)};
        const __m128 z_to_pow4{_mm_mul_ps(z_to_pow2,z_to_pow2)};
        const __m128 z_to_pow5{_mm_mul_ps(z_to_pow3,z_to_pow2)};
        const __m128 z_to_pow7{_mm_mul_ps(z_to_pow5,z_to_pow2)};
        const __m128 z_to_pow9{_mm_mul_ps(z_to_pow7,z_to_pow2)};
        const __m128 horner_series{_mm_fmsub_ps(vC1456,z_to_pow9,
                                      _mm_fmadd_ps(vC124800,z_to_pow7,
                                          _mm_fmadd_ps(vC5491200,z_to_pow5,
                                            _mm_fmsub_ps(vC15627432,z_to_pow4,
                                                _mm_fmadd_ps(vC24960000,z_to_pow3,
                                                    _mm_fmsub_ps(vC111624513,z_to_pow2,vC148832684))))))};
        return (_mm_mul_ps(vCN3756009615384615,horner_series));
    }
    else 
    {
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 factor_zpow2{_mm_fmadd_ps(vC111624513,z_to_pow2,vC148832684)};
        const __m128 z_to_pow3{_mm_mul_ps(z_to_pow2,z)};
        const __m128 factor_zpow3{_mm_mul_ps(vC24960000,z_to_pow3)};
        const __m128 z_to_pow4{_mm_mul_ps(z_to_pow2,z_to_pow2)};
        const __m128 factor_zpow4{_mm_mul_ps(vC15627432,z_to_pow4)};
        const __m128 z_to_pow5{_mm_mul_ps(z_to_pow3,z_to_pow2)};
        const __m128 factor_zpow5{_mm_mul_ps(vC5491200,z_to_pow5)};
        const __m128 z_to_pow7{_mm_mul_ps(z_to_pow5,z_to_pow2)};
        const __m128 factor_zpow7{_mm_mul_ps(vC124800,z_to_pow7)};
        const __m128 z_to_pow9{_mm_mul_ps(z_to_pow7,z_to_pow2)};
        const __m128 factor_zpow9{_mm_mul_ps(vC1456,z_to_pow9)};
        const __m128 partial_factor_1{_mm_add_ps(_mm_sub_ps(factor_zpow9,factor_zpow7),factor_zpow5)};
        const __m128 partial_factor_2{_mm_add_ps(_mm_sub_ps(factor_zpow4,factor_zpow3),factor_zpow2)};
        return (_mm_mul_ps(vCN3756009615384615,_mm_sub_ps(partial_factor_1,partial_factor_2)));
    }
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_econvol_epan6_z_lt_0_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[48] = 
    { 
                                  9.390024038461537E-11f,9.390024038461537E-11,
                                  9.390024038461537E-11f,9.390024038461537E-11f,
                                  2079.0f,2079.0f,2079.0f,2079.0f,
                                  206388.0f,206388.0f,206388.0f,206388.0f,
                                  8867040.0f,8867040.0f,8867040.0f,8867040.0f,
                                  255528000.0f,255528000.0f,255528000.0f,255528000.0f,
                                  515705252.0f,515705252.0f,515705252.0f,515705252.0f,
                                  1681680000.0f,1681680000.0f,1681680000.0f,1681680000.0f,
                                  4922641042.0f,4922641042.0f,4922641042.0f,4922641042.0f,
                                  3057600000.0f,3057600000.0f,3057600000.0f,3057600000.0f,
                                  13674002896.0f,13674002896.0f,13674002896.0f,13674002896.0f,
                                  9015826085.0f,9015826085.0f,9015826085.0f,9015826085.0f,
                                  0.0f,0.0f,0.0f,0.0f
    };

#else 

#endif 
}



} // np_standalone_funcs

} // math

} // gms







#endif /*__GMS_PDF_KERNELS_SSE_H__*/