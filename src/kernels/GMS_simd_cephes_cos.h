/* SIMD (SSE1+MMX or SSE2) implementation of sin, cos, exp and log

   Inspired by Intel Approximate Math library, and based on the
   corresponding algorithms of the cephes math library

   The default is to use the SSE1 version. If you define USE_SSE2 the
   the SSE2 intrinsics will be used in place of the MMX intrinsics. Do
   not expect any significant performance improvement with SSE2.
*/

/* Copyright (C) 2007  Julien Pommier

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  (this is the zlib license)
*/

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1989 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

/*MIT License
!Copyright (c) 2020 Bernard Gingold
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
*/

#ifndef __GMS_SIMD_CEPHES_COS_H__
#define __GMS_SIMD_CEPHES_COS_H__  100620260843

namespace file_info 
{

     static const unsigned int GMS_SIMD_CEPHES_COS_MAJOR = 1;
     static const unsigned int GMS_SIMD_CEPHES_COS_MINOR = 0;
     static const unsigned int GMS_SIMD_CEPHES_COS_MICRO = 0;
     static const unsigned int GMS_SIMD_CEPHES_COS_FULLVER =
       1000U*GMS_SIMD_CEPHES_COS_MAJOR+100U*GMS_SIMD_CEPHES_COS_MINOR+
       10U*GMS_SIMD_CEPHES_COS_MICRO;
     static const char GMS_SIMD_CEPHES_COS_CREATION_DATE[] = "10-06-2026 08:43AM +00200 (WED 10 JUN 2026 GMT+2)";
     static const char GMS_SIM_CEPHES_COS_BUILD_DATE[]     = __DATE__; 
     static const char GMS_SIMD_CEPHES_COS_BUILD_TIME[]    = __TIME__;
     static const char GMS_SIMD_CEPHES_COS_SYNOPSIS[]      = "Manual vectorization of CEPHES cosine(float) function.";

}

#include <immintrin.h>
#include <limits>
#include <cstdint>
#include "GMS_config.h"

namespace gms
{

namespace math 
{

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
__m128 _mm_ceph_cosf_ps(__m128 x) 
{

    const  __m128i inv_sign_mask{_mm_set1_epi32(~0x80000000)};
    const  __m128i sign_mask{_mm_set1_epi32(static_cast<std::int32_t>(0x80000000))};
    const  __m128i ONE{_mm_set1_epi32(1)};
    const  __m128i ONEINV{_mm_set1_epi32(~1)};
    const  __m128i TWO{_mm_set1_epi32(2)};
    const  __m128i FOUR{_mm_set1_epi32(4)};
    const  __m128  FOPI{_mm_set1_ps(1.27323954473516f)};
    const __m128 DP1N{_mm_set1_ps(-0.78515625f)};
	const __m128 DP2N{_mm_set1_ps(-2.4187564849853515625e-4f)};
	const __m128 DP3N{_mm_set1_ps(-3.77489497744594108e-8f)};
    const __m128 C2443315711809948E005{_mm_set1_ps(2.443315711809948E-005)};
    const __m128 CN1388731625493765E003{_mm_set1_ps(-1.388731625493765E-003)};
    const __m128 C4166664568298827E002{_mm_set1_ps(4.166664568298827E-002)};
    const __m128 C05{_mm_set1_ps(0.5f)};
    const __m128 C1{_mm_set1_ps(1.0f)};
    const __m128 CN19515295891E4{_mm_set1_ps(-1.9515295891E-4)};
    const __m128 C83321608736E3{_mm_set1_ps(8.3321608736E-3)};
    const __m128 CN16666654611E1{_mm_set1_ps(-1.6666654611E-1)};
    __m128  xmm1;
    __m128  xmm2;
    __m128  xmm3;
    __m128  y;
    __m128i emm0;
    __m128i emm2;
    
    x           = _mm_and_ps(x,*(__m128*)&inv_sign_mask);
    y           = _mm_mul_ps(x,FOPI);
    emm2        = _mm_cvttps_epi32(y);
    emm2        = _mm_add_epi32(emm2,ONE);
    emm2        = _mm_and_si128(emm2,ONEINV);
    y           = _mm_cvtepi32_ps(emm2);
    emm2        = _mm_sub_epi32(emm2,TWO);
    emm0        = _mm_andnot_si128(emm2,FOUR);
    emm0        = _mm_slli_epi32(emm0,29);
    emm2        = _mm_and_si128(emm2,TWO);
    emm2        = _mm_cmpeq_epi32(emm2,_mm_setzero_si128());
    __m128 sign = _mm_castsi128_ps(emm0);
    __m128 mask = _mm_castsi128_ps(emm2);
    xmm1        = _mm_mul_ps(y,DP1N);
    x           = _mm_add_ps(x,xmm1);
    xmm2        = _mm_mul_ps(y,DP2N);
    x           = _mm_add_ps(x,xmm2);
    xmm3        = _mm_mul_ps(y,DP3N);
    x           = _mm_add_ps(x,xmm3);
    __m128 z    = _mm_mul_ps(x,x);
    y           = C2443315711809948E005;
    y           = _mm_mul_ps(y,z);
    y           = _mm_add_ps(y,CN1388731625493765E003);
    y           = _mm_mul_ps(y,z);
    y           = _mm_add_ps(y,C4166664568298827E002);
    y           = _mm_mul_ps(y,z);
    y           = _mm_mul_ps(y,z);
    __m128 tmp  = _mm_mul_ps(z,C05);
    y           = _mm_sub_ps(y,tmp);
    y           = _mm_add_ps(y,C1);
    __m128 y2   = CN19515295891E4;
    y2          = _mm_mul_ps(y2,z);
    y2          = _mm_add_ps(y2,C83321608736E3);
    y2          = _mm_mul_ps(y2,z);
    y2          = _mm_add_ps(y2,CN16666654611E1);
    y2          = _mm_mul_ps(y2,z);
    y2          = _mm_mul_ps(y2,x);
    y2          = _mm_add_ps(y2,x);
    xmm3        = mask;
    y2          = _mm_and_ps(xmm3,y2);
    y           = _mm_andnot_ps(xmm3,y);
    y           = _mm_add_ps(y,y2);
    y           = _mm_xor_ps(y,sign);
    return (y);
}

} // math

}// gms









#endif /*__GMS_SIMD_CEPHES_COS_H*/