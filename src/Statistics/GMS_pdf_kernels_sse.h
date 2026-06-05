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
    const __mmask8 mask_or{cxy_neq_0 | lambda_neq_0};
    const __m128   term2{_mm_fmsub_ps(cxy,term1,_mm_add_ps(vlambda,vlambda))};
    return (_mm_mask_blend_ps(mask_or,vCn2,_mm_mul_ps(lambda_to_cxypow,term2)));
    //return (_mm_mask_blend_ps(_kor_mask8(cxy_neq_0,lambda_neq_0),vCn2,_mm_mul_ps(lambda_to_cxypow,term2)));
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
template<bool use_prefetching, bool use_horner_scheme>
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
    if constexpr (use_horner_scheme==true)
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
                                            _mm_fmadd_ps(vC15627432,z_to_pow4,
                                                _mm_fmsub_ps(vC24960000,z_to_pow3,
                                                    _mm_fmadd_ps(vC111624513,z_to_pow2,vC148832684))))))};
        return (_mm_mul_ps(vC3756009615384615,horner_series));
    }
    else 
    {
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 factor_zpow2{_mm_fmadd_ps(vC111624513,z_to_pow2,vC148832684)};
        const __m128 z_to_pow3{_mm_mul_ps(z_to_pow2,z)};
        const __m128 factor_zpow3{_mm_mul_ps(vC24960000,z_to_pow3)};
        const __m128 z_to_pow4{_mm_mul_ps(z_to_pow3,z)};
        const __m128 factor_zpow4{_mm_mul_ps(vC15627432,z_to_pow4)};
        const __m128 z_to_pow5{_mm_mul_ps(z_to_pow4,z)};
        const __m128 factor_zpow5{_mm_mul_ps(vC5491200,z_to_pow5)};
        const __m128 z_to_pow7{_mm_mul_ps(z_to_pow5,z_to_pow2)};
        const __m128 factor_zpow7{vC124800,z_to_pow7};
        const __m128 z_to_pow9{_mm_mul_ps(z_to_pow7,z_to_pow2)};
        const __m128 factor_zpow9{_mm_mul_ps(vC1456,z_to_pow9)};
        const __m128 partial_factor_1{_mm_sub_ps(factor_zpow4,_mm_sub_ps(factor_zpow3,factor_zpow2))};
        const __m128 partial_factor_2{_mm_add_ps(_mm_sub_ps(factor_zpow9,factor_zpow7),factor_zpow5)};
        return (_mm_mul_ps(vC3756009615384615,_mm_add_ps(partial_factor_1,partial_factor_2)));
    }
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
template<bool use_prefetching, bool use_horner_scheme>
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
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[32],_MM_HINT_T0);
    }    
        const __m128 vC9390024038461537{_mm_load_ps(&prefetched_constants[0])};
        const __m128 vC2079{_mm_load_ps(&prefetched_constants[4])};
        const __m128 vC206388{_mm_load_ps(&prefetched_constants[8])};
        const __m128 vC8867040{_mm_load_ps(&prefetched_constants[12])};
        const __m128 vC255528000{_mm_load_ps(&prefetched_constants[16])};
        const __m128 vC515705252{_mm_load_ps(&prefetched_constants[20])};
        const __m128 vC1681680000{_mm_load_ps(&prefetched_constants[24])};
        const __m128 vC4922641042{_mm_load_ps(&prefetched_constants[28])};
        const __m128 vC3057600000{_mm_load_ps(&prefetched_constants[32])};
        const __m128 vC13674002896{_mm_load_ps(&prefetched_constants[36])};
        const __m128 vC9015826085{_mm_load_ps(&prefetched_constants[40])};
#else 
        const __m128 vC9390024038461537{_mm_set1_ps(9.390024038461537E-11f)};
        const __m128 vC2079{_mm_set1_ps(2079.0f)};
        const __m128 vC206388{_mm_set1_ps(206388.0f)};
        const __m128 vC8867040{_mm_set1_ps(8867040.0f)};
        const __m128 vC255528000{_mm_set1_ps(255528000.0f)};
        const __m128 vC515705252{_mm_set1_ps(515705252.0f)};
        const __m128 vC1681680000{_mm_set1_ps(1681680000.0f)};
        const __m128 vC4922641042{_mm_set1_ps(4922641042.0f)};
        const __m128 vC3057600000{_mm_set1_ps(3057600000.0f)};
        const __m128 vC13674002896{_mm_set1_ps(13674002896.0f)};
        const __m128 vC9015826085{_mm_load_ps(9015826085.0f)};
#endif 
    if constexpr(use_horner_scheme==true)
    {
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 z_to_pow3{_mm_mul_ps(z_to_pow2,z)};
        const __m128 z_to_pow4{_mm_mul_ps(z_to_pow2,z_to_pow2)};
        const __m128 z_to_pow5{_mm_mul_ps(z_to_pow3,z_to_pow2)};
        const __m128 z_to_pow6{_mm_mul_ps(z_to_pow5,z)};
        const __m128 z_to_pow7{_mm_mul_ps(z_to_pow5,z_to_pow2)};
        const __m128 z_to_pow9{_mm_mul_ps(z_to_pow7,z_to_pow2)};
        const __m128 z_to_pow11{_mm_mul_ps(z_to_pow9,z_to_pow2)};
        const __m128 z_to_pow13{_mm_mul_ps(z_to_pow11,z_to_pow2)};
        // Long micro-ops dependency chain!!
        const __m128 horner_series{_mm_fmsub_ps(vC2079,z_to_pow13,
                                      _mm_fmadd_ps(vC206388,z_to_pow11,
                                         _mm_fmsub_ps(vC8867040,z_to_pow9,
                                            _mm_fmsub_ps(vC255528000,z_to_pow7,
                                               _mm_fmadd_ps(vC515705252,z_to_pow6,
                                                  _mm_fmadd_ps(vC1681680000,z_to_pow5,
                                                     _mm_fmsub_ps(vC4922641042,z_to_pow4,
                                                        _mm_fmsub_ps(vC3057600000,z_to_pow3,
                                                           _mm_fmadd_ps(vC13674002896,z_to_pow2,vC9015826085)))))))))};
        return (_mm_mul_ps(vC9390024038461537,horner_series));
    }
    else 
    {
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 factor_zpow2{_mm_fmadd_ps(vC13674002896,z_to_pow2,vC9015826085)};
        const __m128 z_to_pow3{_mm_mul_ps(z_to_pow2,z)};
        const __m128 factor_zpow3{_mm_mul_ps(vC3057600000,z_to_pow3)};
        const __m128 z_to_pow4{_mm_mul_ps(z_to_pow2,z_to_pow2)};
        const __m128 factor_zpow4{_mm_mul_ps(vC4922641042,z_to_pow4)};
        const __m128 z_to_pow5{_mm_mul_ps(z_to_pow4,z)};
        const __m128 factor_zpow5{_mm_mul_ps(vC1681680000,z_to_pow5)};
        const __m128 z_to_pow6{_mm_mul_ps(z_to_pow5,z)};
        const __m128 factor_zpow6{_mm_mul_ps(vC515705252,z_to_pow6)};
        const __m128 z_to_pow7{_mm_mul_ps(z_to_pow5,z_to_pow2)};
        const __m128 factor_zpow7{_mm_mul_ps(vC255528000,z_to_pow7)};
        const __m128 z_to_pow9{_mm_mul_ps(z_to_pow7,z_to_pow2)};
        const __m128 factor_zpow9{_mm_mul_ps(vC8867040,z_to_pow9)};
        const __m128 z_to_pow11{_mm_mul_ps(z_to_pow9,z_to_pow2)};
        const __m128 factor_zpow11{_mm_mul_ps(vC206388,z_to_pow11)};
        const __m128 z_to_pow13{_mm_mul_ps(z_to_pow11,z_to_pow2)};
        const __m128 factor_zpow13{_mm_mul_ps(vC2079,z_to_pow13)};
        const __m128 partial_factor_1{_mm_add_ps(_mm_sub_ps(factor_zpow13,factor_zpow11),factor_zpow9)};
        const __m128 partial_factor_2{_mm_add_ps(_mm_sub_ps(factor_zpow7,factor_zpow6),factor_zpow5)};
        const __m128 partial_factor_3{_mm_sub_ps(_mm_sub_ps(factor_zpow4,factor_zpow3),factor_zpow2)};
        return (_mm_mul_ps(vC9390024038461537,_mm_add_ps(_mm_sub_ps(partial_factor_1,partial_factor_2),partial_factor_3)));
    }
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching, bool use_horner_scheme>
__m128 np_econvol_epan6_z_gt_0_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[48] = 
    { 
                                 -9.390024038461537E-11f,-9.390024038461537E-11,
                                 -9.390024038461537E-11f,-9.390024038461537E-11f,
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
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[32],_MM_HINT_T0);
    }    
        const __m128 vCN9390024038461537{_mm_load_ps(&prefetched_constants[0])};
        const __m128 vC2079{_mm_load_ps(&prefetched_constants[4])};
        const __m128 vC206388{_mm_load_ps(&prefetched_constants[8])};
        const __m128 vC8867040{_mm_load_ps(&prefetched_constants[12])};
        const __m128 vC255528000{_mm_load_ps(&prefetched_constants[16])};
        const __m128 vC515705252{_mm_load_ps(&prefetched_constants[20])};
        const __m128 vC1681680000{_mm_load_ps(&prefetched_constants[24])};
        const __m128 vC4922641042{_mm_load_ps(&prefetched_constants[28])};
        const __m128 vC3057600000{_mm_load_ps(&prefetched_constants[32])};
        const __m128 vC13674002896{_mm_load_ps(&prefetched_constants[36])};
        const __m128 vC9015826085{_mm_load_ps(&prefetched_constants[40])};
#else 
        const __m128 vCN9390024038461537{_mm_set1_ps(-9.390024038461537E-11f)};
        const __m128 vC2079{_mm_set1_ps(2079.0f)};
        const __m128 vC206388{_mm_set1_ps(206388.0f)};
        const __m128 vC8867040{_mm_set1_ps(8867040.0f)};
        const __m128 vC255528000{_mm_set1_ps(255528000.0f)};
        const __m128 vC515705252{_mm_set1_ps(515705252.0f)};
        const __m128 vC1681680000{_mm_set1_ps(1681680000.0f)};
        const __m128 vC4922641042{_mm_set1_ps(4922641042.0f)};
        const __m128 vC3057600000{_mm_set1_ps(3057600000.0f)};
        const __m128 vC13674002896{_mm_set1_ps(13674002896.0f)};
        const __m128 vC9015826085{_mm_load_ps(9015826085.0f)};
#endif 
    if constexpr(use_horner_scheme==true)
    {
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 z_to_pow3{_mm_mul_ps(z_to_pow2,z)};
        const __m128 z_to_pow4{_mm_mul_ps(z_to_pow2,z_to_pow2)};
        const __m128 z_to_pow5{_mm_mul_ps(z_to_pow3,z_to_pow2)};
        const __m128 z_to_pow6{_mm_mul_ps(z_to_pow5,z)};
        const __m128 z_to_pow7{_mm_mul_ps(z_to_pow5,z_to_pow2)};
        const __m128 z_to_pow9{_mm_mul_ps(z_to_pow7,z_to_pow2)};
        const __m128 z_to_pow11{_mm_mul_ps(z_to_pow9,z_to_pow2)};
        const __m128 z_to_pow13{_mm_mul_ps(z_to_pow11,z_to_pow2)};
        // Long micro-ops dependency chain!!
        const __m128 horner_series{_mm_fmsub_ps(vC2079,z_to_pow13,
                                      _mm_fmadd_ps(vC206388,z_to_pow11,
                                         _mm_fmsub_ps(vC8867040,z_to_pow9,
                                            _mm_fmadd_ps(vC255528000,z_to_pow7,
                                               _mm_fmadd_ps(vC515705252,z_to_pow6,
                                                  _mm_fmsub_ps(vC1681680000,z_to_pow5,
                                                     _mm_fmsub_ps(vC4922641042,z_to_pow4,
                                                        _mm_fmadd_ps(vC3057600000,z_to_pow3,
                                                           _mm_fsub_ps(vC13674002896,z_to_pow2,vC9015826085)))))))))};
        return (_mm_mul_ps(vCN9390024038461537,horner_series));
    }
    else 
    {
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 factor_zpow2{_mm_fmsub_ps(vC13674002896,z_to_pow2,vC9015826085)};
        const __m128 z_to_pow3{_mm_mul_ps(z_to_pow2,z)};
        const __m128 factor_zpow3{_mm_mul_ps(vC3057600000,z_to_pow3)};
        const __m128 z_to_pow4{_mm_mul_ps(z_to_pow2,z_to_pow2)};
        const __m128 factor_zpow4{_mm_mul_ps(vC4922641042,z_to_pow4)};
        const __m128 z_to_pow5{_mm_mul_ps(z_to_pow4,z)};
        const __m128 factor_zpow5{_mm_mul_ps(vC1681680000,z_to_pow5)};
        const __m128 z_to_pow6{_mm_mul_ps(z_to_pow5,z)};
        const __m128 factor_zpow6{_mm_mul_ps(vC515705252,z_to_pow6)};
        const __m128 z_to_pow7{_mm_mul_ps(z_to_pow5,z_to_pow2)};
        const __m128 factor_zpow7{_mm_mul_ps(vC255528000,z_to_pow7)};
        const __m128 z_to_pow9{_mm_mul_ps(z_to_pow7,z_to_pow2)};
        const __m128 factor_zpow9{_mm_mul_ps(vC8867040,z_to_pow9)};
        const __m128 z_to_pow11{_mm_mul_ps(z_to_pow9,z_to_pow2)};
        const __m128 factor_zpow11{_mm_mul_ps(vC206388,z_to_pow11)};
        const __m128 z_to_pow13{_mm_mul_ps(z_to_pow11,z_to_pow2)};
        const __m128 factor_zpow13{_mm_mul_ps(vC2079,z_to_pow13)};
        const __m128 partial_factor_1{_mm_add_ps(_mm_sub_ps(factor_zpow13,factor_zpow11),factor_zpow9)};
        const __m128 partial_factor_2{_mm_add_ps(_mm_add_ps(factor_zpow7,factor_zpow6),factor_zpow5)};
        const __m128 partial_factor_3{_mm_add_ps(_mm_sub_ps(factor_zpow4,factor_zpow3),factor_zpow2)};
        return (_mm_mul_ps(vCN9390024038461537,_mm_sub_ps(_mm_sub_ps(partial_factor_1,partial_factor_2),partial_factor_3)));
    }
}

/*
   ***Warning***
   True branch of this function: 
   double np_econvol_uaa(const int same_cat, const double lambda, const int c){
          const double dcm1 = (double)(c-1);
          const double oml = 1.0-lambda;
          return (same_cat)?(oml*oml+lambda*lambda/dcm1):(lambda/dcm1*(2.0*oml+((double)(c-2))*lambda/dcm1));
}
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
inline 
__m128 np_econvol_uaa_br_true_sse_ps(const __m128 lambda,const __m128 c) 
{
    __ATTR_ALIGN__(16) constexpr float constants[4] = {1.0f,1.0f,1.0f,1.0f};
    const __m128 vC1{_mm_load_ps(&constants[0])};
    const __m128 oml{_mm_sub_ps(vC1,lambda)};
    const __m128 dcm1{_mm_sub_ps(c,vC1)};   
    const __m128 factor_1{_mm_div_ps(lambda,dcm1)};
    return (_mm_fmadd_ps(oml,oml,_mm_mul_ps(lambda,factor_1)));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
inline 
__m128 
np_econvol_uaa_br_false_sse_ps(const __m128 lambda,const __m128 c)
{
    __ATTR_ALIGN__(16) constexpr float constants[8] = {1.0f,1.0f,1.0f,1.0f,
                                                       2.0f,2.0f,2.0f,2.0f};
    const __m128 vC1{_mm_load_ps(&constants[0])};
    const __m128 oml{_mm_sub_ps(vC1,lambda)};
    const __m128 dcm1{_mm_sub_ps(c,vC1)};
    const __m128 vC2{_mm_load_ps(&constants[4])};
    const __m128 factor_1{_mm_div_ps(lambda,dcm1)};
    const __m128 factor_2{_mm_fmadd_ps(vC2,oml,_mm_mul_ps(_mm_sub_ps(c,vC2),factor_1))};
    return (_mm_mul_ps(factor_1,factor_2));                                                   
}

/*
   True branch
   double np_econvol_uli_racine(const int same_cat, const double lambda, const int c){
  return (same_cat)?(1.0 + (double)(c-1)*lambda*lambda):(lambda*(2.0+(double)(c-2)*lambda));
}
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
inline 
__m128 
np_econvol_uli_racine_true_br_sse_ps(const __m128 lambda,const __m128 c)
{
    __ATTR_ALIGN__(16) constexpr float constants[4] = {1.0f,1.0f,1.0f,1.0f};
    const __m128 vC1{_mm_load_ps(&constants[0])};
    const __m128 factor_1{_mm_mul_ps(lambda,lambda)};
    return (_mm_add_ps(vC1,_mm_mul_ps(_mm_sub_ps(c,vC1),factor_1)));
}

/* 
    False branch of above-pasted code
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
inline 
__m128 
np_econvol_uli_racine_false_br_sse_ps(const __m128 lambda,const __m128 c)
{
    __ATTR_ALIGN__(16) constexpr float constants[4] = {2.0f,2.0f,2.0f,2.0f};
    const __m128 vC2{_mm_load_ps(&constants[0])};
    return (_mm_mul_ps(lambda,_mm_fmadd_ps(lambda,_mm_sub_ps(c,vC2),vC2)));
}

/*
   True branch of:
   double np_econvol_unli_racine(const int same_cat, const double lambda, const int c){
   const double inorm = 1.0/((c-1.0)*lambda + 1.0);
   return ((same_cat)?(1.0 + (double)(c-1)*lambda*lambda):(lambda*(2.0+(double)(c-2)*lambda)))*inorm*inorm;
}
*/

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
inline 
__m128 
np_econvol_unli_racine_true_br_sse_ps(const __m128 lambda,const __m128 c)
{
    __ATTR_ALIGN__(16) constexpr float constants[4] = {1.0f,1.0f,1.0f,1.0f};
    const __m128 vC1{_mm_load_ps(&constants[0])};
    const __m128 factor_1{_mm_sub_ps(c,vC1)};
    return (_mm_fmadd_ps(lambda,_mm_mul_ps(lambda,factor_1),vC1));
}

/*
    False branch of above-pasted code
*/

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
inline 
__m128 
np_econvol_unli_racine_true_br_sse_ps(const __m128 lambda,const __m128 c)
{
    __ATTR_ALIGN__(16) constexpr float constants[8] = {1.0f,1.0f,1.0f,1.0f,
                                                       2.0f,2.0f,2.0f,2.0f};
    const __m128 vC1{_mm_load_ps(&constants[0])};
    const __m128 factor_1{_mm_fmadd_ps(_mm_sub_ps(c,vC1),lambda,vC1)};
    const __m128 vC2{_mm_load_ps(&constants[4])};
    const __m128 factor_2{_mm_mul_ps(lambda,_mm_fmadd_ps(lambda,_mm_sub_ps(c,vC2),vC2))};
    const __m128 inorm{_mm_rcp14_ps(factor_1)};
    return (_mm_mul_ps(factor_2,_mm_mul_ps(inorm,inorm)));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
inline
__m128 
np_econvol_onli_racine_sse_ps(const __m128 x,const __m128 y,
                              const __m128 lambda)
{
    __ATTR_ALIGN__(16) constexpr float constants[4] = {1.0f,1.0f,1.0f,1.0f};
    const __m128 vC1{_mm_load_ps(&constants[0])};
    const __m128 cxy{_mm_floor_ps(gms::math::xmm4r4_abs(_mm_sub_ps(x,y)))};
    const __m128 l2{_mm_mul_ps(lambda,lambda)};
    const __m128 lnorm{_mm_div_ps(_mm_sub_ps(vC1,lambda),_mm_sub_ps(vC1,lambda))};
    const __m128 lambda_to_cxypow{_mm_pow_ps(lambda,cxy)};
    const __m128 factor_1{_mm_add_ps(_mm_div_ps(_mm_add_ps(vC1,l2),_mm_sub_ps(vC1,l2)),cxy)};
    const __m128 factor_2{_mm_mul_ps(lnorm,lnorm)};
    return (_mm_mul_ps(factor_2,_mm_mul_ps(lambda_to_cxypow,factor_1)));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_econvol_owang_van_ryzin_sse_ps(const __m128 x,const __m128 y,
                                         const __m128 lambda)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = 
    { 
                        0.5f,0.5f,0.5f,0.5f,
                        1.0f,1.0f,1.0f,1.0f,
                        2.0f,2.0f,2.0f,2.0f,
                        0.0f,0.0f,0.0f,0.0f
    };
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }
       const __m128 vC05{_mm_load_ps(&prefetched_constants[0])};
       const __m128 vC1{_mm_load_ps(&prefetched_constants[4])};
       const __m128 vC2{_mm_load_ps(&prefetched_constants[8])};
#else 
       const __m128 vC05{_mm_set1_ps(0.5f)};
       const __m128 vC1{_mm_set1_ps(1.0f)};
       const __m128 vC2{_mm_set1_ps(2.0f)};
#endif 
       /*
           Intermixed the false and true branch operations
           due to unknown probablity density function (PDF)
           of the comparison result: x == y.
           Could have assumed the uniform density thus assigning
           a value of 1/2 to the PDF, but the real distribution
           must be measured rigouristically.
       */
       const __m128 one_m_lambda{_mm_sub_ps(vC1,lambda)};
       const __m128 lambda_sqr{_mm_mul_ps(lambda,lambda)};
       const __m128 cxy{_mm_floor_ps(gms::math::xmm4r4_abs(_mm_sub_ps(x,y)))};
       const __m128 lnorm{_mm_mul_ps(vC05,one_m_lambda)};
       const __m128 lambda_to_pow_cxy{_mm_pow_ps(lambda,cxy)};
       const __m128 one_m_lambda_sqr{_mm_sub_ps(vC1,lambda_sqr)};
       const __m128 left_term_br_true{_mm_mul_ps(vC05,_mm_mul_ps(one_m_lambda,one_m_lambda))};
       const __m128 left_term_br_false{_mm_mul_ps(_mm_mul_ps(lnorm,lnorm),lambda_to_pow_cxy)};
       const __m128 right_term_br_true{_mm_add_ps(vC1,_mm_rcp14_ps(one_m_lambda_sqr))};
       const __m128 two_div_one_m_l2{_mm_div_ps(vC2,_mm_sub_ps(vC1,lambda_sqr))};
       const __m128 true_br_term{_mm_mul_ps(left_term_br_true,right_term_br_true)};
       const __mmask8 x_eq_y{_mm_cmp_ps_mask(x,y,_CMP_EQ_OQ)};
       const __m128 right_term_br_false{_mm_add_ps(vC1,_mm_add_ps(cxy,two_div_one_m_l2))};
       const __m128 false_br_term{_mm_mul_ps(left_term_br_false,right_term_br_false)};
       return (_mm_mask_blend_ps(x_eq_y,false_br_term,true_br_term));
}

// Derivative kernels

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_deriv_gauss2_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = 
    {
                        0.398942280401432677939946059934f,
                        0.398942280401432677939946059934f,
                        0.398942280401432677939946059934f,
                        0.398942280401432677939946059934f,
                       -0.5f,-0.5f,-0.5f,-0.5f,
                        0.0f,0.0f,0.0f,0.0f,
                        0.0f,0.0f,0.0f,0.0f
    };
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }
        const __m128 vC0398942280401432677939946059934{_mm_load_ps(&prefetched_constants[0])};
        const __m128 vCN05{_mm_load_ps(&prefetched_constants[4])};
#else 
        const __m128 vC0398942280401432677939946059934{_mm_set1_ps(0.398942280401432677939946059934f)};
        const __m128 vCN05{_mm_set1_ps(-0.5f)};
#endif 
        const __m128 neg_z{gms::math::negate_xmm4r4(z)};
        const __m128 exp_arg{_mm_mul_ps(vCN05,_mm_mul_ps(z,z))};
        const __m128 factor_1{_mm_mul_ps(neg_z,vC0398942280401432677939946059934)};
        /*
            ***WARNING***
            The call of SVML exp kernel may not be inlined, thus causing a far jump
            to the address space/range of SVML machine code implementation -- possible
            L1I miss amd maybe a L2 miss also and the memory load can not be excluded.
        */
        const __m128 exp_term{_mm_exp_ps(exp_arg)};
        return (_mm_mul_ps(factor_1,exp_term));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_deriv_gauss4_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = 
    {
                       -0.398942280401432677939946059934f,
                       -0.398942280401432677939946059934f,
                       -0.398942280401432677939946059934f,
                       -0.398942280401432677939946059934f,
                        2.5f,2.5f,2.5f,2.5f,
                        0.5f,0.5f,0.5f,0.5f,
                       -0.5f,-0.5f,-0.5f,-0.5f
    };
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }
        const __m128 vCN0398942280401432677939946059934{_mm_load_ps(&prefetched_constants[0])};
        const __m128 vC25{_mm_load_ps(&prefetched_constants[4])};
        const __m128 vC05{_mm_load_ps(&prefetched_constants[8])};
        const __m128 vCN05{_mm_load_ps(&prefetched_constants[12])};
#else 
        const __m128 vCN0398942280401432677939946059934{_mm_set1_ps(-0.398942280401432677939946059934f)};
        const __m128 vC25{_mm_set1_ps(2.5f)};
        const __m128 vC05{_mm_set1_ps(0.5f)};
        const __m128 vCN05{_mm_set1_ps(-0.5f)};
#endif 
        const __m128 factor_1{_mm_mul_ps(vCN0398942280401432677939946059934,z)};
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 exp_arg{_mm_mul_ps(vCN05,z_to_pow2)};
        const __m128 factor_2{_mm_sub_ps(vC25,_mm_mul_ps(vC05,z_to_pow2))};
         /*
            ***WARNING***
            The call of SVML exp kernel may not be inlined, thus causing a far jump
            to the address space/range of SVML machine code implementation -- possible
            L1I miss amd maybe a L2 miss also and the memory load can not be excluded.
        */
        const __m128 exp_term{_mm_exp_ps(exp_arg)};
        return (_mm_mul_ps(factor_1,_mm_mul_ps(factor_2,exp_term)));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_deriv_gauss6_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = 
    {
                               -0.049867785050179084743f,
                               -0.049867785050179084743f,
                               -0.049867785050179084743f,
                               -0.049867785050179084743f,
                               -0.5f,-0.5f,-0.5f,-0.5f,
                                35.0f,35.0f,35.0f,35.0f,
                               -14.0f,-14.0f-14.0f,-14.0f
    };
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }   
        const __m128 vCN0049867785050179084743{_mm_load_ps(&prefetched_constants[0])};
        const __m128 vCN05{_mm_load_ps(&prefetched_constants[4])};
        const __m128 vC35{_mm_load_ps(&prefetched_constants[8])};
        const __m128 vCN14{_mm_load_ps(&prefetched_constants[12])};
#else 
        const __m128 vCN0049867785050179084743{_mm_set1_ps(-0.049867785050179084743f)};
        const __m128 vCN05{_mm_set1_ps(-0.5f)};
        const __m128 vC35{_mm_set1_ps(35.0f)};
        const __m128 vCN14{_mm_set1_ps(-14.0f)};
#endif 
        const __m128 factor_1{_mm_mul_ps(vCN0049867785050179084743,z)};
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 factor_2{_mm_fmadd_ps(_mm_add_ps(vCN14,z_to_pow2),z_to_pow2,vC35)};
         /*
            ***WARNING***
            The call of SVML exp kernel may not be inlined, thus causing a far jump
            to the address space/range of SVML machine code implementation -- possible
            L1I miss amd maybe a L2 miss also and the memory load can not be excluded.
        */
        const __m128 exp_arg{_mm_mul_ps(vCN05,z_to_pow2)};
        const __m128 exp_term{_mm_exp_ps(exp_arg)};
        return (_mm_mul_ps(factor_1,_mm_mul_ps(exp_term,factor_2)));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif  
template<bool use_prefetching>
__m128 np_deriv_gauss8_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[32] = 
    {
                       -0.398942280401432677939946059934f,
                       -0.398942280401432677939946059934f,
                       -0.398942280401432677939946059934f,
                       -0.398942280401432677939946059934f,
                        6.5625f,6.5625f,6.5625f,6.5625f,
                       -3.9375f,-3.9375f,-3.9375f,-3.9375f,
                        0.5625f,0.5625f,0.5625f,0.5625f,
                        0.02083333333f,0.02083333333f,
                        0.02083333333f,0.02083333333f,
                       -0.5f,-0.5f,-0.5f,-0.5f,
                        0.0f,0.0f,0.0f,0.0f,
                        0.0f,0.0f,0.0f,0.0f
    };
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
    } 
        const __m128 vCN0398942280401432677939946059934{_mm_load_ps(&prefetched_constants[0])};
        const __m128 vC65625{_mm_load_ps(&prefetched_constants[4])};
        const __m128 vCN39375{_mm_load_ps(&prefetched_constants[8])};
        const __m128 vC05625{_mm_load_ps(&prefetched_constants[12])};
        const __m128 vC002083333333{_mm_load_ps(&prefetched_constants[16])};
        const __m128 vCN05{_mm_load_ps(&prefetched_constants[20])};
#else 
        const __m128 vCN0398942280401432677939946059934{_mm_set1_ps(-0.398942280401432677939946059934)};
        const __m128 vC65625{_mm_set1_ps(6.5625f)};
        const __m128 vCN39375{_mm_set1_ps(3.9375f)};
        const __m128 vC05625{_mm_set1_ps(0.5625f)};
        const __m128 vC002083333333{_mm_set1_ps(0.02083333333f)};
        const __m128 vCN05{_mm_set1_ps(-0.5f)};
#endif
        const __m128 factor_1{_mm_mul_ps(vCN0398942280401432677939946059934,z)};
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 factor_2{_mm_sub_ps(vC05625,_mm_mul_ps(vC002083333333,z_to_pow2))};
        const __m128 exp_arg{_mm_mul_ps(vCN05,z_to_pow2)};
        const __m128 horner_series{_mm_fmadd_ps(z_to_pow2,
                                            _mm_fmadd_ps(factor_2,z_to_pow2,vCN39375),vC65625)};
         /*
            ***WARNING***
            The call of SVML exp kernel may not be inlined, thus causing a far jump
            to the address space/range of SVML machine code implementation -- possible
            L1I miss amd maybe a L2 miss also and the memory load can not be excluded.
        */
        const __m128 exp_term{_mm_exp_ps(exp_arg)};
        return (_mm_mul_ps(_mm_mul_ps(factor_1,horner_series),exp_term));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
inline 
__m128 np_deriv_epan2_sse_ps(const __m128 z)
{
       __ATTR_ALIGN__(16) constexpr float constants[4] = 
       {
                            -0.13416407864998738178f,
                            -0.13416407864998738178f,
                            -0.13416407864998738178f,
                            -0.13416407864998738178f    
       };
       return (_mm_mul_ps(_mm_load_ps(&constants[0]),z));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
inline 
__m128 np_deriv_epan4_sse_ps(const __m128 z)
{
       __ATTR_ALIGN__(16) constexpr float constants[8] = 
       {
                           2.347871374742824e-1f,
                           2.347871374742824e-1f,
                           2.347871374742824e-1f,
                           2.347871374742824e-1f,
                           8.385254921942804e-1f,
                           8.385254921942804e-1f,
                           8.385254921942804e-1f,
                           8.385254921942804e-1f
       };
       const __m128 z_to_pow2{_mm_mul_ps(z,z)};
       return (_mm_mul_ps(z,_mm_fmsub_ps(_mm_load_ps(&constants[0]),z_to_pow2,_mm_load_ps(&constants[4]))));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_prefetching>
__m128 np_deriv_epan6_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = 
    {
                          -2.567984320334919f,
                          -2.567984320334919f,
                          -2.567984320334919f,
                          -2.567984320334919f,
                           1.848948710641142f,
                           1.848948710641142f,
                           1.848948710641142f,
                           1.848948710641142f,
                           2.905490831007508e-1f,
                           2.905490831007508e-1f,
                           2.905490831007508e-1f,
                           2.905490831007508e-1f,
                           0.0f,0.0f,0.0f,0.0f
    };
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }  
        const __m128 vCN2567984320334919{_mm_load_ps(&prefetched_constants[0])};
        const __m128 vC1848948710641142{_mm_load_ps(&prefetched_constants[4])};
        const __m128 vC2905490831007508{_mm_load_ps(&prefetched_constants[8])};
#else 
        const __m128 vCN2567984320334919{_mm_set1_ps(-2.567984320334919f)};
        const __m128 vC1848948710641142{_mm_set1_ps(1.848948710641142f)};
        const __m128 vC2905490831007508{_mm_set1_ps(2.905490831007508e-1f)};
#endif 
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 factor_1{_mm_sub_ps(vC1848948710641142,_mm_mul_ps(vC2905490831007508,z_to_pow2))};
        const __m128 factor_2{_mm_fmadd_ps(factor_1,z_to_pow2,vCN2567984320334919)};
        return (_mm_mul_ps(z,factor_2));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_prefetching>
__m128 np_deriv_epan8_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = 
    {
                              -5.777964720753567f,
                              -5.777964720753567f,
                              -5.777964720753567f,
                              -5.777964720753567f,
                               7.626913431394709f,
                               7.626913431394709f,
                               7.626913431394709f,
                               7.626913431394709f,
                              -2.83285356023232f,
                              -2.83285356023232f,
                              -2.83285356023232f,
                              -2.83285356023232f,
                               3.147615066924801e-1f,
                               3.147615066924801e-1f,
                               3.147615066924801e-1f,
                               3.147615066924801e-1f
    };
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }  
        const __m128 vCN5777964720753567{_mm_load_ps(&prefetched_constants[0])};
        const __m128 vC7626913431394709{_mm_load_ps(&prefetched_constants[4])};
        const __m128 vCN283285356023232{_mm_load_ps(&prefetched_constants[8])};
        const __m128 vC3147615066924801{_mm_load_ps(&prefetched_constants[12])};
#else 
        const __m128 vCN5777964720753567{_mm_set1_ps(-5.777964720753567f)};
        const __m128 vC7626913431394709{_mm_set1_ps(7.626913431394709f)};
        const __m128 vCN283285356023232{_mm_set1_ps(2.83285356023232f)};
        const __m128 vC3147615066924801{_mm_set1_ps(3.147615066924801e-1f)};
#endif 
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 horner_series{_mm_fmadd_ps(
                                            _mm_fmadd_ps(
                                                     _mm_fmadd_ps(z_to_pow2,vC3147615066924801,vCN283285356023232),
                                                                                            z_to_pow2,vC7626913431394709),
                                                                                                       z_to_pow2,vCN5777964720753567)};
        return (_mm_mul_ps(z,horner_series));
}

// CDF kernels

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_prefetching>
__m128 np_cdf_gauss2_sse_ps(const __m128 z)
{
#if (PDF_KERNELS_SSE_OPTIMIZE_OUT_RIP_RODATA_STORE) == 1
    __ATTR_ALIGN__(16) constexpr float prefetched_constants[16] = 
    {
                           0.5f,0.5f,0.5f,0.5f,
                           0.7071067810f,0.7071067810f,
                           0.7071067810f,0.7071067810f,
                           0.0f,0.0f,0.0f,0.0f,
                           0.0f,0.0f,0.0f,0.0f
    };
    if constexpr(use_prefetching==true)
    {
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
    }  
        const __m128 vC05{_mm_load_ps(&prefetched_constants[0])};
        const __m128 vC07071067810{_mm_load_ps(&prefetched_constants[4])};
#else 
        const __m128 vC05{_mm_set1_ps(0.5f)};
        const __m128 vC07071067810{_mm_set1_ps(0.7071067810f)};
#endif
        return (_mm_fmadd_ps(vC05,_mm_erf_ps(_mm_mul_ps(vC07071067810,z),vC05)));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_prefetching,bool optimize_out_rip_store>
__m128 np_cdf_gauss4_sse_ps(const __m128 z)
{
    __m128 vC05;
    __m128 vC07071067810;
    __m128 vC01994711401;
    __m128 vCN05;
    if constexpr(static_cast<std::int32_t>(optimize_out_rip_store)==static_cast<std::int32_t>(true))
    {
        __ATTR_ALIGN__(16) constexpr const float prefetched_constants[16] = 
        {
                        0.5f,0.5f,0.5f,0.5f,
                        0.7071067810f,0.7071067810f,
                        0.7071067810f,0.7071067810f,
                        0.1994711401f,0.1994711401f,
                        0.1994711401f,0.1994711401f,
                       -0.5f,-0.5f,-0.5f,-0.5f
        };
        if constexpr(static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
        {
            _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
        } 
        vC05          = _mm_load_ps((const float*)&prefetched_constants[0]);
        vC07071067810 = _mm_load_ps((const float*)&prefetched_constants[4]);
        vC01994711401 = _mm_load_ps((const float*)&prefetched_constants[8]);
        vCN05         = _mm_load_ps((const float*)&prefetched_constants[12]);
    }
    else 
    {
        vC05          = _mm_set1_ps(0.5f);
        vC07071067810 = _mm_set1_ps(0.7071067810f);
        vC01994711401 = _mm_set1_ps(0.1994711401f);
        vCN05         = _mm_set1_ps(-0.5f);
    } 
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 exp_arg{_mm_mul_ps(vCN05,z_to_pow2)};
        const __m128 erf_arg{_mm_mul_ps(vC07071067810,z)};
        const __m128 exp_val{_mm_exp_ps(exp_arg)};
        const __m128 factor_1{_mm_mul_ps(vC01994711401,z)};
        const __m128 erf_val{_mm_erf_ps(erf_arg)};
        return (_mm_fmadd_ps(vC05,erf_val,_mm_fmadd_ps(factor_1,exp_val,vC05)));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_prefetching,bool optimize_out_rip_store>
__m128 np_cdf_gauss6_sse_ps(const __m128 z)
{
    __m128 vC05;
    __m128 vCN05; 
    __m128 vC07071067810f; 
    __m128 vC03490744952f; 
    __m128 vC004986778504f; 
    if constexpr(static_cast<std::int32_t>(optimize_out_rip_store)==static_cast<std::int32_t>(true))
    {
        __ATTR_ALIGN__(16) constexpr float prefetched_constants[32] = 
        {
                        0.5f,0.5f,0.5f,0.5f,
                       -0.5f,-0.5f,-0.5f,-0.5f,
                        0.7071067810f,0.7071067810f,
                        0.7071067810f,0.7071067810f,
                        0.3490744952f,0.3490744952f,
                        0.3490744952f,0.3490744952f,
                        0.04986778504f,0.04986778504f,
                        0.04986778504f,0.04986778504f,
                        0.0f,0.0f,0.0f,0.0f,
                        0.0f,0.0f,0.0f,0.0f,
                        0.0f,0.0f,0.0f,0.0f
        };
        if constexpr(static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
        {
            _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
            _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
        }
            vC05            = _mm_load_ps((const float*)&prefetched_constants[0]);
            vCN05           = _mm_load_ps((const float*)&prefetched_constants[4]);
            vC07071067810f  = _mm_load_ps((const float*)&prefetched_constants[8]);
            vC03490744952f  = _mm_load_ps((const float*)&prefetched_constants[12]);
            vC004986778504f = _mm_load_ps((const float*)&prefetched_constants[16]);
    }
    else 
    { 
            vC05            = _mm_set1_ps(0.5f);
            vCN05           = _mm_set1_ps(-0.5f);
            vC07071067810f  = _mm_set1_ps(0.7071067810f);
            vC03490744952f  = _mm_set1_ps(0.3490744952f);
            vC004986778504f = _mm_set1_ps(0.04986778504f);
    } 
            const __m128 z_to_pow2{_mm_mul_ps(z,z)};
            const __m128 erf_arg{_mm_mul_ps(vC07071067810f,z)};
            const __m128 erf_val{_mm_erf_ps(erf_arg)};
            const __m128 factor_1{_mm_sub_ps(vC03490744952f,_mm_mul_ps(vC004986778504f,z))};
            const __m128 exp_arg{_mm_mul_ps(vCN05,z_to_pow2)};
            const __m128 exp_val{_mm_exp_ps(exp_arg)};
            return (_mm_fmadd_ps(vC05,erf_val,_mm_fmadd_ps(z,_mm_mul_ps(exp_val,factor_1),vC05)));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_prefetching,bool optimize_out_rip_store>
__m128 np_cdf_gauss8_sse_ps(const __m128 z)
{
    __m128 vC05;            
    __m128 vC07071067810;   
    __m128 vCN05;           
    __m128 vC04737439578;   
    __m128 vCN01329807601;  
    __m128 vC0008311297511;
    if constexpr(static_cast<std::int32_t>(optimize_out_rip_store)==static_cast<std::int32_t>(true))
    {
        __ATTR_ALIGN__(16) constexpr const float prefetched_constants[32] = 
        {
                        0.5f,0.5f,0.5f,0.5f,
                        0.7071067810f,0.7071067810f,
                        0.7071067810f,0.7071067810f,
                       -0.5f,-0.5f,-0.5f,-0.5f,
                        0.4737439578f,0.4737439578f,
                        0.4737439578f,0.4737439578f,
                       -0.1329807601f,-0.1329807601f,
                       -0.1329807601f,-0.1329807601f,
                        0.008311297511f,0.008311297511f,
                        0.008311297511f,0.008311297511f,
                        0.0f,0.0f,0.0f,0.0f,
                        0.0f,0.0f,0.0f,0.0f
        };
        if constexpr(static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
        {
            _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
            _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
        }
        vC05            = _mm_load_ps((const float*)&prefetched_constants[0]);
        vC07071067810   = _mm_load_ps((const float*)&prefetched_constants[4]);
        vCN05           = _mm_load_ps((const float*)&prefetched_constants[8]);
        vC04737439578   = _mm_load_ps((const float*)&prefetched_constants[12]);
        vCN01329807601  = _mm_load_ps((const float*)&prefetched_constants[16]);
        vC0008311297511 = _mm_load_ps((const float*)&prefetched_constants[20]);
    }
    else 
    {
        vC05            = _mm_set1_ps(0.5f);
        vC07071067810   = _mm_set1_ps(0.7071067810f);
        vCN05           = _mm_set1_ps(-0.5f);
        vC04737439578   = _mm_set1_ps(0.4737439578f);
        vCN01329807601  = _mm_set1_ps(-0.1329807601f);
        vC0008311297511 = _mm_set1_ps(0.008311297511f);
    }
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 erf_arg{_mm_mul_ps(vC07071067810,z)};
        const __m128 erf_val{_mm_erf_ps(erf_arg)};
        const __m128 horner_series{_mm_fmadd_ps(
                                      _mm_fmadd_ps(z_to_pow2,vC0008311297511,vCN01329807601),
                                                                            z_to_pow2,vC04737439578)};
        const __m128 exp_arg{_mm_mul_ps(vCN05,z_to_pow2)};
        const __m128 exp_val{_mm_mul_ps(z,_mm_exp_ps(exp_val))};
        return (_mm_fmadd_ps(vC05,erf_val,_mm_fmadd_ps(exp_val,horner_series,vC05)));
}

/*
*** WARNING ***
    Removed the input checking in order to eliminate tje if-convertion overhead!!
    Below is the disassembly of the scalar original version, pay attention to 
    costly double branch input checking.
 np_cdf_epan2(double):
 vmovapd xmm1,xmm0
 vmovsd xmm0,QWORD PTR [rip+0x0]        # c <np_cdf_epan2(double)+0xc>
    R_X86_64_PC32 .LC2-0x4
 vcomisd xmm0,xmm1
 ja     28 <np_cdf_epan2(double)+0x28>
 vcomisd xmm1,QWORD PTR [rip+0x0]        # 1a <np_cdf_epan2(double)+0x1a>
    R_X86_64_PC32 .LC3-0x4
 vmovsd xmm0,QWORD PTR [rip+0x0]        # 22 <np_cdf_epan2(double)+0x22>
    R_X86_64_PC32 .LC0-0x4
 jbe    30 <np_cdf_epan2(double)+0x30>
 ret
 nop    DWORD PTR [rax]
 vxorpd xmm0,xmm0,xmm0
 ret
 nop    DWORD PTR [rax]
 vmulsd xmm2,xmm1,QWORD PTR [rip+0x0]        # 38 <np_cdf_epan2(double)+0x38>
    R_X86_64_PC32 .LC4-0x4
 vmovsd xmm0,QWORD PTR [rip+0x0]        # 40 <np_cdf_epan2(double)+0x40>
    R_X86_64_PC32 .LC5-0x4
 vmulsd xmm2,xmm2,xmm1
 vsubsd xmm0,xmm0,xmm2
 vmulsd xmm0,xmm0,xmm1
 vaddsd xmm0,xmm0,QWORD PTR [rip+0x0]        # 54 <np_cdf_epan2(double)+0x54>
    R_X86_64_PC32 .LC6-0x4
 ret
*/

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_prefetching,bool optimize_out_rip_store>
__m128 np_cdf_epan2_sse_ps(const __m128 z)
{
    __m128 vC03354101967;
    __m128 vC002236067978; 
    __m128 vC05;
    if constexpr(static_cast<std::int32_t>(optimize_out_rip_store)==static_cast<std::int32_t>(true))
    {
        __ATTR_ALIGN__(16) constexpr const float prefetched_constants[16] = 
        {
                                0.3354101967f,0.3354101967f,
                                0.3354101967f,0.3354101967f,
                                0.02236067978f,0.02236067978f,
                                0.02236067978f,0.02236067978f,
                                0.5f,0.5f,0.5f,0.5f,
                                0.0f,0.0f,0.0f,0.0f
        };
        if constexpr(static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
        {
            _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
        }
            vC03354101967   = _mm_load_ps((const float*)&prefetched_constants[0]);
            vC002236067978  = _mm_load_ps((const float*)&prefetched_constants[4]);
            vC05            = _mm_load_ps((const float*)&prefetched_constants[8]);
    }
    else
    {
            vC03354101967   = _mm_set1_ps(0.3354101967f);
            vC002236067978  =_mm_set1_ps(0.02236067978f);
            vC05            = _mm_set1_ps(0.5f);
    }
            const __m128 z_to_pow2{_mm_mul_ps(z,z)};
            const __m128 factor_1{_mm_sub_ps(vC03354101967,_mm_mul_ps(vC002236067978,z_to_pow2))};
            return (_mm_fmadd_ps(z,factor_1,vC05));
}


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_prefetching,bool optimize_out_rip_store>
__m128 np_cdf_epan4_sse_ps(const __m128 z)
{   
    __m128 vC05;
    __m128 vC06288941188f;
    __m128 vCN01397542486;
    __m128 vC001173935688;
    if constexpr(static_cast<std::int32_t>(optimize_out_rip_store)==static_cast<std::int32_t>(true))
    {
        __ATTR_ALIGN__(16) constexpr const float prefetched_constants[16] = 
        {
                           0.5f,0.5f,0.5f,0.5f,
                           0.6288941188f,0.6288941188f,
                           0.6288941188f,0.6288941188f,
                          -0.1397542486f,-0.1397542486f,
                          -0.1397542486f,-0.1397542486f,
                           0.01173935688f,0.01173935688f,
                           0.01173935688f,0.01173935688f
        };
        if constexpr(static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
        {
             _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
        }
        vC05            = _mm_load_ps((const float*)&prefetched_constants[0]);
        vC06288941188f  = _mm_load_ps((const float*)&prefetched_constants[4]);
        vCN01397542486  = _mm_load_ps((const float*)&prefetched_constants[8]);
        vC001173935688  = _mm_load_ps((const float*)&prefetched_constants[12]);
    }
    else 
    {
        vC05            = _mm_set1_ps(0.5f);
        vC06288941188f  = _mm_set1_ps(0.6288941188f);
        vCN01397542486  = _mm_set1_ps(-0.1397542486f);
        vC001173935688  = _mm_set1_ps(0.01173935688f);
    } 
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        return (_mm_fmadd_ps(
                         _mm_fmadd_ps(
                                  _mm_fmadd_ps(z_to_pow2,vC001173935688,vCN01397542486),
                                                                        z_to_pow2,vC06288941188f),
                                                                                                z,vC05));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_prefetching,bool optimize_out_rip_store>
__m128 np_cdf_epan6_sse_ps(const __m128 z)
{
    __m128 vC05;
    __m128 vC09171372566;
    __m128 vCN04279973864;
    __m128 vC009244743547;
    __m128 vC0006917835307;
    if constexpr(static_cast<std::int32_t>(optimize_out_rip_store)==static_cast<std::int32_t>(true))
    {
        __ATTR_ALIGN__(16) constexpr const float prefetched_constants[32] = 
        {
                        0.5f,0.5f,0.5f,0.5f,
                        0.9171372566f,0.9171372566f,
                        0.9171372566f,0.9171372566f,
                       -0.4279973864f,-0.4279973864f,
                       -0.4279973864f,-0.4279973864f,
                        0.09244743547f,0.09244743547f,
                        0.09244743547f,0.09244743547f,
                        0.006917835307f,0.006917835307f,
                        0.006917835307f,0.006917835307f,
                        0.0f,0.0f,0.0f,0.0f,
                        0.0f,0.0f,0.0f,0.0f,
                        0.0f,0.0f,0.0f,0.0f
        };
        if constexpr(static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
        {
            _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
            _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
        }
            vC05            = _mm_load_ps((const float*)&prefetched_constants[0]);
            vC09171372566   = _mm_load_ps((const float*)&prefetched_constants[4]);
            vCN04279973864  = _mm_load_ps((const float*)&prefetched_constants[8]);
            vC009244743547  = _mm_load_ps((const float*)&prefetched_constants[12]);
            vC0006917835307 = _mm_load_ps((const float*)&prefetched_constants[16]);
    }
    else 
    {
            vC05            = _mm_set1_ps(0.5f);
            vC09171372566   = _mm_set1_ps(0.9171372566f);
            vCN04279973864  = _mm_set1_ps(-0.4279973864f);
            vC009244743547  = _mm_set1_ps(0.09244743547f);
            vC0006917835307 = _mm_set1_ps(0.006917835307f);
    }
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 factor_1{_mm_sub_ps(vC009244743547,
                                     _mm_mul_ps(vC0006917835307,z_to_pow2))};
        return (_mm_fmadd_ps(
                         _mm_fmadd_ps(
                                  _mm_fmadd_ps(factor_1,z_to_pow2,vCN04279973864),
                                                        z_to_pow2,vC009244743547),z,vC05));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_prefetching,bool optimize_out_rip_store>
__m128 np_cdf_epan8_sse_ps(const __m128 z)
{
    __m128 vC05;
    __m128 vC1203742649;
    __m128 vCN09629941194;
    __m128 vC03813456714;
    __m128 vCN006744889424;
    __m128 vC0004371687590;
    if constexpr(static_cast<std::int32_t>(optimize_out_rip_store)==static_cast<std::int32_t>(true))
    {    
        __ATTR_ALIGN__(16) constexpr const float prefetched_constants[32] = 
        {
                        0.5f,0.5f,0.5f,0.5f,
                        1.203742649f,1.203742649f,
                        1.203742649f,1.203742649f,
                       -0.9629941194f,-0.9629941194f,
                       -0.9629941194f,-0.9629941194f,
                        0.3813456714f,0.3813456714f,
                        0.3813456714f,0.3813456714f,
                       -0.06744889424f,-0.06744889424f,
                       -0.06744889424f,-0.06744889424f,
                        0.004371687590f,0.004371687590f,
                        0.004371687590f,0.004371687590f,
                        0.0f,0.0f,0.0f,0.0f,
                        0.0f,0.0f,0.0f,0.0f
        };
        if constexpr(static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
        {
            _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
            _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
        }
            vC05            = _mm_load_ps((const float*)&prefetched_constants[0]);
            vC1203742649    = _mm_load_ps((const float*)&prefetched_constants[4]);
            vCN09629941194  = _mm_load_ps((const float*)&prefetched_constants[8]);
            vC03813456714   = _mm_load_ps((const float*)&prefetched_constants[12]);
            vCN006744889424 = _mm_load_ps((const float*)&prefetched_constants[16]);
            vC0004371687590 = _mm_load_ps((const float*)&prefetched_constants[20]);
    }
    else 
    {
            vC05            =_mm_set1_ps(0.5f);
            vC1203742649    = _mm_set1_ps(1.203742649f);
            vCN09629941194  = _mm_set1_ps(-0.9629941194f);
            vC03813456714   = _mm_set1_ps(0.3813456714f);
            vCN006744889424 = _mm_set1_ps(-0.06744889424f);
            vC0004371687590 = _mm_set1_ps(0.004371687590f);
    }
        const __m128 z_to_pow2{_mm_mul_ps(z,z)};
        const __m128 factor_1{_mm_fmadd_ps(z_to_pow2,vC0004371687590,vCN006744889424)};
        return (_mm_fmadd_ps(
                         _mm_fmadd_ps(
                                  _mm_fmadd_ps(
                                           _mm_fmadd_ps(factor_1,z_to_pow2,vC03813456714),
                                                                     z_to_pow2,vCN09629941194),
                                                                             z_to_pow2,vC1203742649),z,vC05));
}



} // np_standalone_funcs

} // math

} // gms







#endif /*__GMS_PDF_KERNELS_SSE_H__*/