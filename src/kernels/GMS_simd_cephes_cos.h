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

namespace simd_cos_constants
{

__ATTR_ALIGN__(16)
__attribute__((section(".rodata"))) 
constexpr const static float prefetched_constants_4x_f32[48] =
{
    1.27323954473516f,1.27323954473516f,
    1.27323954473516f,1.27323954473516f,
   -0.78515625f,     -0.78515625f,
   -0.78515625f,     -0.78515625f,
   -2.4187564849853515625e-4f,-2.4187564849853515625e-4f,
   -2.4187564849853515625e-4f,-2.4187564849853515625e-4f,
   -3.77489497744594108e-8f,-3.77489497744594108e-8f,
   -3.77489497744594108e-8f,-3.77489497744594108e-8f,
    2.443315711809948E-005f,2.443315711809948E-005f,
    2.443315711809948E-005f,2.443315711809948E-005f,
   -1.388731625493765E-003f,-1.388731625493765E-003f,
   -1.388731625493765E-003f,-1.388731625493765E-003f,
    4.166664568298827E-002f,4.166664568298827E-002f,
    4.166664568298827E-002f,4.166664568298827E-002f,
    0.5f,0.5f,0.5f,0.5f,1.0f,1.0f,1.0f,1.0f,
   -1.9515295891E-4f,-1.9515295891E-4f,
   -1.9515295891E-4f,-1.9515295891E-4f,
    8.3321608736E-3f,8.3321608736E-3f,
    8.3321608736E-3f,8.3321608736E-3f,
   -1.6666654611E-1f,-1.6666654611E-1f,
   -1.6666654611E-1f,-1.6666654611E-1f
};

__ATTR_ALIGN__(16)
__attribute__((section(".rodata"))) 
constexpr const static std::int32_t prefetched_constants_4x_i32[32] = 
{
   ~0x80000000,~0x80000000,
   ~0x80000000,~0x80000000,
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   1,1,1,1,
   ~1,~1,~1,~1,
   2,2,2,2,
   4,4,4,4,
   0,0,0,0,
   0,0,0,0
};

__ATTR_ALIGN__(32)
__attribute__((section(".rodata"))) 
constexpr const static float prefetched_constants_8x_f32[96] =
{
   1.27323954473516f,1.27323954473516f,
   1.27323954473516f,1.27323954473516f,
   1.27323954473516f,1.27323954473516f,
   1.27323954473516f,1.27323954473516f,
  -0.78515625f,     -0.78515625f,
  -0.78515625f,     -0.78515625f,
  -0.78515625f,     -0.78515625f,
  -0.78515625f,     -0.78515625f,
  -2.4187564849853515625e-4f,-2.4187564849853515625e-4f,
  -2.4187564849853515625e-4f,-2.4187564849853515625e-4f,
  -2.4187564849853515625e-4f,-2.4187564849853515625e-4f,
  -2.4187564849853515625e-4f,-2.4187564849853515625e-4f,
  -3.77489497744594108e-8f,-3.77489497744594108e-8f,
  -3.77489497744594108e-8f,-3.77489497744594108e-8f,
  -3.77489497744594108e-8f,-3.77489497744594108e-8f,
  -3.77489497744594108e-8f,-3.77489497744594108e-8f,
   2.443315711809948E-005f,2.443315711809948E-005f,
   2.443315711809948E-005f,2.443315711809948E-005f,
   2.443315711809948E-005f,2.443315711809948E-005f,
   2.443315711809948E-005f,2.443315711809948E-005f,
  -1.388731625493765E-003f,-1.388731625493765E-003f,
  -1.388731625493765E-003f,-1.388731625493765E-003f,
  -1.388731625493765E-003f,-1.388731625493765E-003f,
  -1.388731625493765E-003f,-1.388731625493765E-003f,
   4.166664568298827E-002f,4.166664568298827E-002f,
   4.166664568298827E-002f,4.166664568298827E-002f,
   4.166664568298827E-002f,4.166664568298827E-002f,
   4.166664568298827E-002f,4.166664568298827E-002f,
   0.5f,0.5f,0.5f,0.5f, 0.5f,0.5f,0.5f,0.5f,
   1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,
  -1.9515295891E-4f,-1.9515295891E-4f,
  -1.9515295891E-4f,-1.9515295891E-4f,
  -1.9515295891E-4f,-1.9515295891E-4f,
  -1.9515295891E-4f,-1.9515295891E-4f,
   8.3321608736E-3f,8.3321608736E-3f,
   8.3321608736E-3f,8.3321608736E-3f,
   8.3321608736E-3f,8.3321608736E-3f,
   8.3321608736E-3f,8.3321608736E-3f,
  -1.6666654611E-1f,-1.6666654611E-1f,
  -1.6666654611E-1f,-1.6666654611E-1f,
  -1.6666654611E-1f,-1.6666654611E-1f,
  -1.6666654611E-1f,-1.6666654611E-1f
};

__ATTR_ALIGN__(32)
__attribute__((section(".rodata"))) 
constexpr const static std::int32_t prefetched_constants_8x_i32[48] =
{
   ~0x80000000,~0x80000000,
   ~0x80000000,~0x80000000,
   ~0x80000000,~0x80000000,
   ~0x80000000,~0x80000000,
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   1,1,1,1,1,1,1,1,
  ~1,~1,~1,~1,~1,~1,~1,~1,
   2,2,2,2,2,2,2,2,
   4,4,4,4,4,4,4,4
};

__ATTR_ALIGN__(64)
__attribute__((section(".rodata"))) 
constexpr const static float prefetched_constants_16x_f32[192] =
{
   1.27323954473516f,1.27323954473516f,
   1.27323954473516f,1.27323954473516f,
   1.27323954473516f,1.27323954473516f,
   1.27323954473516f,1.27323954473516f,
   1.27323954473516f,1.27323954473516f,
   1.27323954473516f,1.27323954473516f,
   1.27323954473516f,1.27323954473516f,
   1.27323954473516f,1.27323954473516f,
  -0.78515625f,     -0.78515625f,
  -0.78515625f,     -0.78515625f,
  -0.78515625f,     -0.78515625f,
  -0.78515625f,     -0.78515625f,
  -0.78515625f,     -0.78515625f,
  -0.78515625f,     -0.78515625f,
  -0.78515625f,     -0.78515625f,
  -0.78515625f,     -0.78515625f,
  -2.4187564849853515625e-4f,-2.4187564849853515625e-4f,
  -2.4187564849853515625e-4f,-2.4187564849853515625e-4f,
  -2.4187564849853515625e-4f,-2.4187564849853515625e-4f,
  -2.4187564849853515625e-4f,-2.4187564849853515625e-4f,
  -2.4187564849853515625e-4f,-2.4187564849853515625e-4f,
  -2.4187564849853515625e-4f,-2.4187564849853515625e-4f,
  -2.4187564849853515625e-4f,-2.4187564849853515625e-4f,
  -2.4187564849853515625e-4f,-2.4187564849853515625e-4f,
  -3.77489497744594108e-8f,-3.77489497744594108e-8f,
  -3.77489497744594108e-8f,-3.77489497744594108e-8f,
  -3.77489497744594108e-8f,-3.77489497744594108e-8f,
  -3.77489497744594108e-8f,-3.77489497744594108e-8f,
  -3.77489497744594108e-8f,-3.77489497744594108e-8f,
  -3.77489497744594108e-8f,-3.77489497744594108e-8f,
  -3.77489497744594108e-8f,-3.77489497744594108e-8f,
  -3.77489497744594108e-8f,-3.77489497744594108e-8f,
   2.443315711809948E-005f,2.443315711809948E-005f,
   2.443315711809948E-005f,2.443315711809948E-005f,
   2.443315711809948E-005f,2.443315711809948E-005f,
   2.443315711809948E-005f,2.443315711809948E-005f,
   2.443315711809948E-005f,2.443315711809948E-005f,
   2.443315711809948E-005f,2.443315711809948E-005f,
   2.443315711809948E-005f,2.443315711809948E-005f,
   2.443315711809948E-005f,2.443315711809948E-005f,
  -1.388731625493765E-003f,-1.388731625493765E-003f,
  -1.388731625493765E-003f,-1.388731625493765E-003f,
  -1.388731625493765E-003f,-1.388731625493765E-003f,
  -1.388731625493765E-003f,-1.388731625493765E-003f,
  -1.388731625493765E-003f,-1.388731625493765E-003f,
  -1.388731625493765E-003f,-1.388731625493765E-003f,
  -1.388731625493765E-003f,-1.388731625493765E-003f,
  -1.388731625493765E-003f,-1.388731625493765E-003f,
   4.166664568298827E-002f,4.166664568298827E-002f,
   4.166664568298827E-002f,4.166664568298827E-002f,
   4.166664568298827E-002f,4.166664568298827E-002f,
   4.166664568298827E-002f,4.166664568298827E-002f,
   4.166664568298827E-002f,4.166664568298827E-002f,
   4.166664568298827E-002f,4.166664568298827E-002f,
   4.166664568298827E-002f,4.166664568298827E-002f,
   4.166664568298827E-002f,4.166664568298827E-002f,
   0.5f,0.5f,0.5f,0.5f, 0.5f,0.5f,0.5f,0.5f,
   0.5f,0.5f,0.5f,0.5f, 0.5f,0.5f,0.5f,0.5f,
   1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,
   1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,
  -1.9515295891E-4f,-1.9515295891E-4f,
  -1.9515295891E-4f,-1.9515295891E-4f,
  -1.9515295891E-4f,-1.9515295891E-4f,
  -1.9515295891E-4f,-1.9515295891E-4f,
  -1.9515295891E-4f,-1.9515295891E-4f,
  -1.9515295891E-4f,-1.9515295891E-4f,
  -1.9515295891E-4f,-1.9515295891E-4f,
  -1.9515295891E-4f,-1.9515295891E-4f,
   8.3321608736E-3f,8.3321608736E-3f,
   8.3321608736E-3f,8.3321608736E-3f,
   8.3321608736E-3f,8.3321608736E-3f,
   8.3321608736E-3f,8.3321608736E-3f,
   8.3321608736E-3f,8.3321608736E-3f,
   8.3321608736E-3f,8.3321608736E-3f,
   8.3321608736E-3f,8.3321608736E-3f,
   8.3321608736E-3f,8.3321608736E-3f,
  -1.6666654611E-1f,-1.6666654611E-1f,
  -1.6666654611E-1f,-1.6666654611E-1f,
  -1.6666654611E-1f,-1.6666654611E-1f,
  -1.6666654611E-1f,-1.6666654611E-1f
  -1.6666654611E-1f,-1.6666654611E-1f,
  -1.6666654611E-1f,-1.6666654611E-1f,
  -1.6666654611E-1f,-1.6666654611E-1f,
  -1.6666654611E-1f,-1.6666654611E-1f
};

__ATTR_ALIGN__(64)
__attribute__((section(".rodata"))) 
constexpr const static std::int32_t prefetched_constants_16x_i32[96] =
{
   ~0x80000000,~0x80000000,
   ~0x80000000,~0x80000000,
   ~0x80000000,~0x80000000,
   ~0x80000000,~0x80000000,
   ~0x80000000,~0x80000000,
   ~0x80000000,~0x80000000,
   ~0x80000000,~0x80000000,
   ~0x80000000,~0x80000000,
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   static_cast<std::int32_t>(0x80000000),
   1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
  ~1,~1,~1,~1,~1,~1,~1,~1,~1,~1,~1,~1,~1,~1,~1,~1,
   2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
   4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4
};

} // simd_cos_constants

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
    const  __m128  DP1N{_mm_set1_ps(-0.78515625f)};
	 const  __m128  DP2N{_mm_set1_ps(-2.4187564849853515625e-4f)};
	 const  __m128  DP3N{_mm_set1_ps(-3.77489497744594108e-8f)};
    const  __m128  C2443315711809948E005{_mm_set1_ps(2.443315711809948E-005f)};
    const  __m128  CN1388731625493765E003{_mm_set1_ps(-1.388731625493765E-003f)};
    const  __m128  C4166664568298827E002{_mm_set1_ps(4.166664568298827E-002f)};
    const  __m128  C05{_mm_set1_ps(0.5f)};
    const  __m128  C1{_mm_set1_ps(1.0f)};
    const  __m128  CN19515295891E4{_mm_set1_ps(-1.9515295891E-4f)};
    const  __m128  C83321608736E3{_mm_set1_ps(8.3321608736E-3f)};
    const  __m128  CN16666654611E1{_mm_set1_ps(-1.6666654611E-1f)};
    __m128  xmm1;
    __m128  xmm2;
    __m128  xmm3;
    __m128  y;
    __m128i emm0;
    __m128i emm2;
    /*
        Unfortunately the achieved accurracy is between the 2ULP to 3ULP
        and it is completely dependent upon the radian's range of the input data.
    */
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
    // Long dependency chain which may stall pipeline on the 'y' register interdependent result computation
    // Shall be considering the Estrin's scheme (if possible).
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
    // Long dependency chain which may stall pipeline on the 'y2' register interdependent result computation
    // Shall be considering the Estrin's scheme (if possible).
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

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_prefetching,bool optimize_out_rip_store>
inline static __m128 _mm_ceph_cosf_ps_v2(__m128 x) 
{
   __m128i inv_sign_mask;
   __m128i sign_mask;
   __m128i ONE;
   __m128i ONEINV;
   __m128i TWO;
   __m128i FOUR;
   __m128  FOPI;
   __m128  DP1N;
	__m128  DP2N;
	__m128  DP3N;
   __m128  C2443315711809948E005;
   __m128  CN1388731625493765E003;
   __m128  C4166664568298827E002;
   __m128  C05;
   __m128  C1;
   __m128  CN19515295891E4;
   __m128  C83321608736E3;
   __m128  CN16666654611E1;
   __m128  xmm1;
   __m128  xmm2;
   __m128  xmm3;
   __m128  y;
   __m128i emm0;
   __m128i emm2;

   if constexpr(static_cast<std::int32_t>(optimize_out_rip_store)==static_cast<std::int32_t>(true))
   {
	   const float *        __restrict__ ptr_prefetched_constants_4x_f32 = &simd_cos_constants::prefetched_constants_4x_f32[0];
		const std::int32_t * __restrict__ ptr_prefetched_constants_4x_i32 = &simd_cos_constants::prefetched_constants_4x_i32[0];
      if constexpr(static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
      {
            _mm_prefetch((const char*)&ptr_prefetched_constants_4x_i32[0],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_4x_f32[0],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_4x_i32[16],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_4x_f32[16],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_4x_f32[32],_MM_HINT_T0);
      }
      inv_sign_mask          = _mm_load_epi32((const std::int32_t*)&ptr_prefetched_constants_4x_i32[0]);
      FOPI                   = _mm_load_ps((const float*)&ptr_prefetched_constants_4x_f32[0]);
      sign_mask              = _mm_load_epi32((const std::int32_t*)&ptr_prefetched_constants_4x_i32[4]);
      DP1N                   = _mm_load_ps((const float*)&ptr_prefetched_constants_4x_f32[4]);
      ONE                    = _mm_load_epi32((const std::int32_t*)&ptr_prefetched_constants_4x_i32[8]);
      DP2N                   = _mm_load_ps((const float*)&ptr_prefetched_constants_4x_f32[8]);
      ONEINV                 = _mm_load_epi32((const std::int32_t*)&ptr_prefetched_constants_4x_i32[12]);
      DP3N                   = _mm_load_ps((const float*)&ptr_prefetched_constants_4x_f32[12]);
      TWO                    = _mm_load_epi32((const std::int32_t*)&ptr_prefetched_constants_4x_i32[16]);
      C2443315711809948E005  = _mm_load_ps((const float*)&ptr_prefetched_constants_4x_f32[16]);
      FOUR                   = _mm_load_epi32((const std::int32_t*)&ptr_prefetched_constants_4x_i32[20]);
      CN1388731625493765E003 = _mm_load_ps((const float*)&ptr_prefetched_constants_4x_f32[20]);
      C4166664568298827E002  = _mm_load_ps((const float*)&ptr_prefetched_constants_4x_f32[24]);
      C05                    = _mm_load_ps((const float*)&ptr_prefetched_constants_4x_f32[28]);
      C1                     = _mm_load_ps((const float*)&ptr_prefetched_constants_4x_f32[32]);
      CN19515295891E4        = _mm_load_ps((const float*)&ptr_prefetched_constants_4x_f32[36]);
      C83321608736E3         = _mm_load_ps((const float*)&ptr_prefetched_constants_4x_f32[40]);
      CN16666654611E1        = _mm_load_ps((const float*)&ptr_prefetched_constants_4x_f32[44]);
   }
   else 
   {
      inv_sign_mask          = _mm_set1_epi32(~0x80000000);
      sign_mask              = _mm_set1_epi32(static_cast<std::int32_t>(0x80000000));
      ONE                    = _mm_set1_epi32(1);
      ONEINV                 = _mm_set1_epi32(~1);
      TWO                    = _mm_set1_epi32(2);
      FOUR                   = _mm_set1_epi32(4);
      FOPI                   = _mm_set1_ps(1.27323954473516f);
      DP1N                   = _mm_set1_ps(-0.78515625f);
	   DP2N                   = _mm_set1_ps(-2.4187564849853515625e-4f);
	   DP3N                   = _mm_set1_ps(-3.77489497744594108e-8f);
      C2443315711809948E005  = _mm_set1_ps(2.443315711809948E-005f);
      CN1388731625493765E003 = _mm_set1_ps(-1.388731625493765E-003f);
      C4166664568298827E002  = _mm_set1_ps(4.166664568298827E-002f);
      C05                    = _mm_set1_ps(0.5f);
      C1                     = _mm_set1_ps(1.0f);
      CN19515295891E4        = _mm_set1_ps(-1.9515295891E-4f);
      C83321608736E3         = _mm_set1_ps(8.3321608736E-3f);
      CN16666654611E1        = _mm_set1_ps(-1.6666654611E-1f);
   }

    /*
        Unfortunately the achieved accurracy is between the 2ULP to 3ULP
        and it is completely dependent upon the radian's range of the input data.
    */
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
    // Long dependency chain which may stall pipeline on the 'y' register interdependent result computation
    // Shall be considering the Estrin's scheme (if possible).
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
    // Long dependency chain which may stall pipeline on the 'y2' register interdependent result computation
    // Shall be considering the Estrin's scheme (if possible).
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

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX2
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx2")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
__m256 _mm256_ceph_cosf_ps(__m256 x) 
{

    const  __m256i inv_sign_mask{_mm256_set1_epi32(~0x80000000)};
    const  __m256i sign_mask{_mm256_set1_epi32(static_cast<std::int32_t>(0x80000000))};
    const  __m256i ONE{_mm256_set1_epi32(1)};
    const  __m256i ONEINV{_mm256_set1_epi32(~1)};
    const  __m256i TWO{_mm256_set1_epi32(2)};
    const  __m256i FOUR{_mm256_set1_epi32(4)};
    const  __m256  FOPI{_mm256_set1_ps(1.27323954473516f)};
    const  __m256  DP1N{_mm256_set1_ps(-0.78515625f)};
	 const  __m256  DP2N{_mm256_set1_ps(-2.4187564849853515625e-4f)};
	 const  __m256  DP3N{_mm256_set1_ps(-3.77489497744594108e-8f)};
    const  __m256  C2443315711809948E005{_mm256_set1_ps(2.443315711809948E-005f)};
    const  __m256  CN1388731625493765E003{_mm256_set1_ps(-1.388731625493765E-003f)};
    const  __m256  C4166664568298827E002{_mm256_set1_ps(4.166664568298827E-002f)};
    const  __m256  C05{_mm256_set1_ps(0.5f)};
    const  __m256  C1{_mm256_set1_ps(1.0f)};
    const  __m256  CN19515295891E4{_mm256_set1_ps(-1.9515295891E-4f)};
    const  __m256  C83321608736E3{_mm256_set1_ps(8.3321608736E-3f)};
    const  __m256  CN16666654611E1{_mm256_set1_ps(-1.6666654611E-1f)};
    __m256  xmm1;
    __m256  xmm2;
    __m256  xmm3;
    __m256  y;
    __m256i emm0;
    __m256i emm2;
    /*
        Unfortunately the achieved accurracy is between the 2ULP to 3ULP
        and it is completely dependent upon the radian's range of the input data.
    */
    x           = _mm256_and_ps(x,*(__m256*)&inv_sign_mask);
    y           = _mm256_mul_ps(x,FOPI);
    emm2        = _mm256_cvttps_epi32(y);
    emm2        = _mm256_add_epi32(emm2,ONE);
    emm2        = _mm256_and_si256(emm2,ONEINV);
    y           = _mm256_cvtepi32_ps(emm2);
    emm2        = _mm256_sub_epi32(emm2,TWO);
    emm0        = _mm256_andnot_si256(emm2,FOUR);
    emm0        = _mm256_slli_epi32(emm0,29);
    emm2        = _mm256_and_si256(emm2,TWO);
    emm2        = _mm256_cmpeq_epi32(emm2,_mm256_setzero_si256());
    __m256 sign = _mm256_castsi256_ps(emm0);
    __m256 mask = _mm256_castsi256_ps(emm2);
    xmm1        = _mm256_mul_ps(y,DP1N);
    x           = _mm256_add_ps(x,xmm1);
    xmm2        = _mm256_mul_ps(y,DP2N);
    x           = _mm256_add_ps(x,xmm2);
    xmm3        = _mm256_mul_ps(y,DP3N);
    x           = _mm256_add_ps(x,xmm3);
    __m256 z    = _mm256_mul_ps(x,x);
    // Long dependency chain which may stall pipeline on the 'y' register interdependent result computation
    // Shall be considering the Estrin's scheme (if possible).
    y           = C2443315711809948E005;
    y           = _mm256_mul_ps(y,z);
    y           = _mm256_add_ps(y,CN1388731625493765E003);
    y           = _mm256_mul_ps(y,z);
    y           = _mm256_add_ps(y,C4166664568298827E002);
    y           = _mm256_mul_ps(y,z);
    y           = _mm256_mul_ps(y,z);
    __m256 tmp  = _mm256_mul_ps(z,C05);
    y           = _mm256_sub_ps(y,tmp);
    y           = _mm256_add_ps(y,C1);
    __m256 y2   = CN19515295891E4;
    // Long dependency chain which may stall pipeline on the 'y2' register interdependent result computation
    // Shall be considering the Estrin's scheme (if possible).
    y2          = _mm256_mul_ps(y2,z);
    y2          = _mm256_add_ps(y2,C83321608736E3);
    y2          = _mm256_mul_ps(y2,z);
    y2          = _mm256_add_ps(y2,CN16666654611E1);
    y2          = _mm256_mul_ps(y2,z);
    y2          = _mm256_mul_ps(y2,x);
    y2          = _mm256_add_ps(y2,x);
    xmm3        = mask;
    y2          = _mm256_and_ps(xmm3,y2);
    y           = _mm256_andnot_ps(xmm3,y);
    y           = _mm256_add_ps(y,y2);
    y           = _mm256_xor_ps(y,sign);
    return (y);
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX2
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx2")
#endif
template<bool use_prefetching,bool optimize_out_rip_store>
inline static __m256 _mm256_ceph_cosf_ps_v2(__m256 x) 
{
   __m256i inv_sign_mask;
   __m256i sign_mask;
   __m256i ONE;
   __m256i ONEINV;
   __m256i TWO;
   __m256i FOUR;
   __m256  FOPI;
   __m256  DP1N;
	__m256  DP2N;
	__m256  DP3N;
   __m256  C2443315711809948E005;
   __m256  CN1388731625493765E003;
   __m256  C4166664568298827E002;
   __m256  C05;
   __m256  C1;
   __m256  CN19515295891E4;
   __m256  C83321608736E3;
   __m256  CN16666654611E1;
   __m256  xmm1;
   __m256  xmm2;
   __m256  xmm3;
   __m256  y;
   __m256i emm0;
   __m256i emm2;

   if constexpr(static_cast<std::int32_t>(optimize_out_rip_store)==static_cast<std::int32_t>(true))
   {
	   const float *        __restrict__ ptr_prefetched_constants_8x_f32 = &simd_cos_constants::prefetched_constants_8x_f32[0];
		const std::int32_t * __restrict__ ptr_prefetched_constants_8x_i32 = &simd_cos_constants::prefetched_constants_8x_i32[0];
      if constexpr(static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
      {
            _mm_prefetch((const char*)&ptr_prefetched_constants_8x_i32[0],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_8x_f32[0],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_8x_i32[16],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_8x_f32[16],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_8x_i32[32],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_8x_f32[32],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_8x_f32[48],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_8x_f32[64],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_8x_f32[80],_MM_HINT_T0);
      }
      inv_sign_mask          = _mm256_load_epi32((const std::int32_t*)&ptr_prefetched_constants_8x_i32[0]);
      FOPI                   = _mm256_load_ps((const float*)&ptr_prefetched_constants_8x_f32[0]);
      sign_mask              = _mm256_load_epi32((const std::int32_t*)&ptr_prefetched_constants_8x_i32[8]);
      DP1N                   = _mm256_load_ps((const float*)&ptr_prefetched_constants_8x_f32[8]);
      ONE                    = _mm256_load_epi32((const std::int32_t*)&ptr_prefetched_constants_8x_i32[16]);
      DP2N                   = _mm256_load_ps((const float*)&ptr_prefetched_constants_8x_f32[16]);
      ONEINV                 = _mm256_load_epi32((const std::int32_t*)&ptr_prefetched_constants_8x_i32[24]);
      DP3N                   = _mm256_load_ps((const float*)&ptr_prefetched_constants_8x_f32[24]);
      TWO                    = _mm256_load_epi32((const std::int32_t*)&ptr_prefetched_constants_8x_i32[32]);
      C2443315711809948E005  = _mm256_load_ps((const float*)&ptr_prefetched_constants_8x_f32[32]);
      FOUR                   = _mm256_load_epi32((const std::int32_t*)&ptr_prefetched_constants_8x_i32[40]);
      CN1388731625493765E003 = _mm256_load_ps((const float*)&ptr_prefetched_constants_8x_f32[40]);
      C4166664568298827E002  = _mm256_load_ps((const float*)&ptr_prefetched_constants_8x_f32[48]);
      C05                    = _mm256_load_ps((const float*)&ptr_prefetched_constants_8x_f32[56]);
      C1                     = _mm256_load_ps((const float*)&ptr_prefetched_constants_8x_f32[64]);
      CN19515295891E4        = _mm256_load_ps((const float*)&ptr_prefetched_constants_8x_f32[72]);
      C83321608736E3         = _mm256_load_ps((const float*)&ptr_prefetched_constants_8x_f32[80]);
      CN16666654611E1        = _mm256_load_ps((const float*)&ptr_prefetched_constants_8x_f32[88]);
   }
   else 
   {
      inv_sign_mask          = _mm256_set1_epi32(~0x80000000);
      sign_mask              = _mm256_set1_epi32(static_cast<std::int32_t>(0x80000000));
      ONE                    = _mm256_set1_epi32(1);
      ONEINV                 = _mm256_set1_epi32(~1);
      TWO                    = _mm256_set1_epi32(2);
      FOUR                   = _mm256_set1_epi32(4);
      FOPI                   = _mm256_set1_ps(1.27323954473516f);
      DP1N                   = _mm256_set1_ps(-0.78515625f);
	   DP2N                   = _mm256_set1_ps(-2.4187564849853515625e-4f);
	   DP3N                   = _mm256_set1_ps(-3.77489497744594108e-8f);
      C2443315711809948E005  = _mm256_set1_ps(2.443315711809948E-005f);
      CN1388731625493765E003 = _mm256_set1_ps(-1.388731625493765E-003f);
      C4166664568298827E002  = _mm256_set1_ps(4.166664568298827E-002f);
      C05                    = _mm256_set1_ps(0.5f);
      C1                     = _mm256_set1_ps(1.0f);
      CN19515295891E4        = _mm256_set1_ps(-1.9515295891E-4f);
      C83321608736E3         = _mm256_set1_ps(8.3321608736E-3f);
      CN16666654611E1        = _mm256_set1_ps(-1.6666654611E-1f);
   }

    /*
        Unfortunately the achieved accurracy is between the 2ULP to 3ULP
        and it is completely dependent upon the radian's range of the input data.
    */
    x           = _mm256_and_ps(x,*(__m256*)&inv_sign_mask);
    y           = _mm256_mul_ps(x,FOPI);
    emm2        = _mm256_cvttps_epi32(y);
    emm2        = _mm256_add_epi32(emm2,ONE);
    emm2        = _mm256_and_si256(emm2,ONEINV);
    y           = _mm256_cvtepi32_ps(emm2);
    emm2        = _mm256_sub_epi32(emm2,TWO);
    emm0        = _mm256_andnot_si256(emm2,FOUR);
    emm0        = _mm256_slli_epi32(emm0,29);
    emm2        = _mm256_and_si256(emm2,TWO);
    emm2        = _mm256_cmpeq_epi32(emm2,_mm256_setzero_si256());
    __m256 sign = _mm256_castsi256_ps(emm0);
    __m256 mask = _mm256_castsi256_ps(emm2);
    xmm1        = _mm256_mul_ps(y,DP1N);
    x           = _mm256_add_ps(x,xmm1);
    xmm2        = _mm256_mul_ps(y,DP2N);
    x           = _mm256_add_ps(x,xmm2);
    xmm3        = _mm256_mul_ps(y,DP3N);
    x           = _mm256_add_ps(x,xmm3);
    __m256 z    = _mm256_mul_ps(x,x);
    // Long dependency chain which may stall pipeline on the 'y' register interdependent result computation
    // Shall be considering the Estrin's scheme (if possible).
    y           = C2443315711809948E005;
    y           = _mm256_mul_ps(y,z);
    y           = _mm256_add_ps(y,CN1388731625493765E003);
    y           = _mm256_mul_ps(y,z);
    y           = _mm256_add_ps(y,C4166664568298827E002);
    y           = _mm256_mul_ps(y,z);
    y           = _mm256_mul_ps(y,z);
    __m256 tmp  = _mm256_mul_ps(z,C05);
    y           = _mm256_sub_ps(y,tmp);
    y           = _mm256_add_ps(y,C1);
    __m256 y2   = CN19515295891E4;
    // Long dependency chain which may stall pipeline on the 'y2' register interdependent result computation
    // Shall be considering the Estrin's scheme (if possible).
    y2          = _mm256_mul_ps(y2,z);
    y2          = _mm256_add_ps(y2,C83321608736E3);
    y2          = _mm256_mul_ps(y2,z);
    y2          = _mm256_add_ps(y2,CN16666654611E1);
    y2          = _mm256_mul_ps(y2,z);
    y2          = _mm256_mul_ps(y2,x);
    y2          = _mm256_add_ps(y2,x);
    xmm3        = mask;
    y2          = _mm256_and_ps(xmm3,y2);
    y           = _mm256_andnot_ps(xmm3,y);
    y           = _mm256_add_ps(y,y2);
    y           = _mm256_xor_ps(y,sign);
    return (y);
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
__m512 _mm512_ceph_cosf_ps(__m512 x) 
{

    const  __m512i inv_sign_mask{_mm512_set1_epi32(~0x80000000)};
    const  __m512i sign_mask{_mm512_set1_epi32(static_cast<std::int32_t>(0x80000000))};
    const  __m512i ONE{_mm512_set1_epi32(1)};
    const  __m512i ONEINV{_mm512_set1_epi32(~1)};
    const  __m512i TWO{_mm512_set1_epi32(2)};
    const  __m512i FOUR{_mm512_set1_epi32(4)};
    const  __m512  FOPI{_mm512_set1_ps(1.27323954473516f)};
    const  __m512  DP1N{_mm512_set1_ps(-0.78515625f)};
	 const  __m512  DP2N{_mm512_set1_ps(-2.4187564849853515625e-4f)};
	 const  __m512  DP3N{_mm512_set1_ps(-3.77489497744594108e-8f)};
    const  __m512  C2443315711809948E005{_mm512_set1_ps(2.443315711809948E-005f)};
    const  __m512  CN1388731625493765E003{_mm512_set1_ps(-1.388731625493765E-003f)};
    const  __m512  C4166664568298827E002{_mm512_set1_ps(4.166664568298827E-002f)};
    const  __m512  C05{_mm512_set1_ps(0.5f)};
    const  __m512  C1{_mm512_set1_ps(1.0f)};
    const  __m512  CN19515295891E4{_mm512_set1_ps(-1.9515295891E-4f)};
    const  __m512  C83321608736E3{_mm512_set1_ps(8.3321608736E-3f)};
    const  __m512  CN16666654611E1{_mm512_set1_ps(-1.6666654611E-1f)};
    __m512  zmm1;
    __m512  zmm2;
    __m512  zmm3;
    __m512  y;
    __m512i emm0;
    __m512i emm2;
    //__mmask16 emm2_eq_zero;
    /*
        Unfortunately the achieved accurracy is between the 2ULP to 3ULP
        and it is completely dependent upon the radian's range of the input data.
        The full precision with an error of 0ULP!! is attained only in the range
        of [-PI,+PI]
    */
    x               = _mm512_and_ps(x,*(__m512*)&inv_sign_mask);
    y               = _mm512_mul_ps(x,FOPI);
    emm2            = _mm512_cvttps_epi32(y);
    emm2            = _mm512_add_epi32(emm2,ONE);
    emm2            = _mm512_and_si512(emm2,ONEINV);
    y               = _mm512_cvtepi32_ps(emm2);
    emm2            = _mm512_sub_epi32(emm2,TWO);
    emm0            = _mm512_andnot_si512(emm2,FOUR);
    emm0            = _mm512_slli_epi32(emm0,29);
    emm2            = _mm512_and_si512(emm2,TWO);
    // Temporary workaround, shall be looking for more efficient implementation
    __m256i lo_emm2 = _mm512_extracti32x8_epi32(emm2,0);
    __m256i hi_emm2 = _mm512_extracti32x8_epi32(emm2,1);
    lo_emm2         = _mm256_cmpeq_epi32(lo_emm2,_mm256_setzero_si256());
    hi_emm2         = _mm256_cmpeq_epi32(hi_emm2,_mm256_setzero_si256());
    // Causing wrong values.
    //emm2_eq_zero = _mm512_cmp_epi32_mask(emm2,_mm512_setzero_si512(),_MM_CMPINT_EQ);
    //emm2        = _mm512_mask_mov_epi32(_mm512_setzero_si512(),emm2_eq_zero,ONE);
    emm2            = _mm512_inserti32x8(emm2,lo_emm2,0);
    emm2            = _mm512_inserti32x8(emm2,hi_emm2,1);
    __m512 sign     = _mm512_castsi512_ps(emm0);
    __m512 mask     = _mm512_castsi512_ps(emm2);
    zmm1            = _mm512_mul_ps(y,DP1N);
    x               = _mm512_add_ps(x,zmm1);
    zmm2            = _mm512_mul_ps(y,DP2N);
    x               = _mm512_add_ps(x,zmm2);
    zmm3            = _mm512_mul_ps(y,DP3N);
    x               = _mm512_add_ps(x,zmm3);
    __m512 z        = _mm512_mul_ps(x,x);
    // Long dependency chain which may stall pipeline on the 'y' register interdependent result computation
    // Shall be considering the Estrin's scheme (if possible).
    y               = C2443315711809948E005;
    y               = _mm512_mul_ps(y,z);
    y               = _mm512_add_ps(y,CN1388731625493765E003);
    y               = _mm512_mul_ps(y,z);
    y               = _mm512_add_ps(y,C4166664568298827E002);
    y               = _mm512_mul_ps(y,z);
    y               = _mm512_mul_ps(y,z);
    __m512 tmp      = _mm512_mul_ps(z,C05);
    y               = _mm512_sub_ps(y,tmp);
    y               = _mm512_add_ps(y,C1);
    __m512 y2       = CN19515295891E4;
    // Long dependency chain which may stall pipeline on the 'y2' register interdependent result computation
    // Shall be considering the Estrin's scheme (if possible).
    y2              = _mm512_mul_ps(y2,z);
    y2              = _mm512_add_ps(y2,C83321608736E3);
    y2              = _mm512_mul_ps(y2,z);
    y2              = _mm512_add_ps(y2,CN16666654611E1);
    y2              = _mm512_mul_ps(y2,z);
    y2              = _mm512_mul_ps(y2,x);
    y2              = _mm512_add_ps(y2,x);
    zmm3            = mask;
    y2              = _mm512_and_ps(zmm3,y2);
    y               = _mm512_andnot_ps(zmm3,y);
    y               = _mm512_add_ps(y,y2);
    y               = _mm512_xor_ps(y,sign);
    return (y);
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool optimize_out_rip_store>
inline static __m512 _mm512_ceph_cosf_ps_v2(__m512 x) 
{
   __m512i inv_sign_mask;
   __m512i sign_mask;
   __m512i ONE;
   __m512i ONEINV;
   __m512i TWO;
   __m512i FOUR;
   __m512  FOPI;
   __m512  DP1N;
	__m512  DP2N;
	__m512  DP3N;
   __m512  C2443315711809948E005;
   __m512  CN1388731625493765E003;
   __m512  C4166664568298827E002;
   __m512  C05;
   __m512  C1;
   __m512  CN19515295891E4;
   __m512  C83321608736E3;
   __m512  CN16666654611E1;
   __m512  zmm1;
   __m512  zmm2;
   __m512  zmm3;
   __m512  y;
   __m512i emm0;
   __m512i emm2;
    //__mmask16 emm2_eq_zero;
   if constexpr(static_cast<std::int32_t>(optimize_out_rip_store)==static_cast<std::int32_t>(true))
   {
	   const float *        __restrict__ ptr_prefetched_constants_16x_f32 = &simd_cos_constants::prefetched_constants_16x_f32[0];
		const std::int32_t * __restrict__ ptr_prefetched_constants_16x_i32 = &simd_cos_constants::prefetched_constants_16x_i32[0];
      /*
        ***IMPORTANT*** -- Disable this prefetching block after the first call (in the case of cold I1D) cache 
                           of the _mm512_cpehes_cos_ps_v2 kernel.
                           This may be done before the main loop entry.
      */
      if constexpr(static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
      {
            _mm_prefetch((const char*)&ptr_prefetched_constants_16x_i32[0],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[0],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_16x_i32[16],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[16],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_16x_i32[32],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[32],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_16x_i32[48],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[48],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_16x_i32[64],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[64],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_16x_i32[80],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[80],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[96],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[112],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[128],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[144],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[160],_MM_HINT_T0);
            _mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[176],_MM_HINT_T0);

      }
      inv_sign_mask          = _mm512_load_epi32((const std::int32_t*)&ptr_prefetched_constants_16x_i32[0]);
      FOPI                   = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[0]);
      sign_mask              = _mm512_load_epi32((const std::int32_t*)&ptr_prefetched_constants_16x_i32[16]);
      DP1N                   = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[16]);
      ONE                    = _mm512_load_epi32((const std::int32_t*)&ptr_prefetched_constants_16x_i32[32]);
      DP2N                   = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[32]);
      ONEINV                 = _mm512_load_epi32((const std::int32_t*)&ptr_prefetched_constants_16x_i32[48]);
      DP3N                   = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[48]);
      TWO                    = _mm512_load_epi32((const std::int32_t*)&ptr_prefetched_constants_16x_i32[64]);
      C2443315711809948E005  = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[64]);
      FOUR                   = _mm512_load_epi32((const std::int32_t*)&ptr_prefetched_constants_16x_i32[80]);
      CN1388731625493765E003 = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[80]);
      C4166664568298827E002  = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[96]);
      C05                    = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[112]);
      C1                     = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[128]);
      CN19515295891E4        = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[144]);
      C83321608736E3         = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[160]);
      CN16666654611E1        = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[176]);
   }
   else 
   {
      inv_sign_mask          = _mm512_set1_epi32(~0x80000000);
      sign_mask              = _mm512_set1_epi32(static_cast<std::int32_t>(0x80000000));
      ONE                    = _mm512_set1_epi32(1);
      ONEINV                 = _mm512_set1_epi32(~1);
      TWO                    = _mm512_set1_epi32(2);
      FOUR                   = _mm512_set1_epi32(4);
      FOPI                   = _mm512_set1_ps(1.27323954473516f);
      DP1N                   = _mm512_set1_ps(-0.78515625f);
	   DP2N                   = _mm512_set1_ps(-2.4187564849853515625e-4f);
	   DP3N                   = _mm512_set1_ps(-3.77489497744594108e-8f);
      C2443315711809948E005  = _mm512_set1_ps(2.443315711809948E-005f);
      CN1388731625493765E003 = _mm512_set1_ps(-1.388731625493765E-003f);
      C4166664568298827E002  = _mm512_set1_ps(4.166664568298827E-002f);
      C05                    = _mm512_set1_ps(0.5f);
      C1                     = _mm512_set1_ps(1.0f);
      CN19515295891E4        = _mm512_set1_ps(-1.9515295891E-4f);
      C83321608736E3         = _mm512_set1_ps(8.3321608736E-3f);
      CN16666654611E1        = _mm512_set1_ps(-1.6666654611E-1f);
   }
   
    /*
        Unfortunately the achieved accurracy is between the 2ULP to 3ULP
        and it is completely dependent upon the radian's range of the input data.
        The full precision with an error of 0ULP!! is attained only in the range
        of [-PI,+PI]
    */
    x               = _mm512_and_ps(x,*(__m512*)&inv_sign_mask);
    y               = _mm512_mul_ps(x,FOPI);
    emm2            = _mm512_cvttps_epi32(y);
    emm2            = _mm512_add_epi32(emm2,ONE);
    emm2            = _mm512_and_si512(emm2,ONEINV);
    y               = _mm512_cvtepi32_ps(emm2);
    emm2            = _mm512_sub_epi32(emm2,TWO);
    emm0            = _mm512_andnot_si512(emm2,FOUR);
    emm0            = _mm512_slli_epi32(emm0,29);
    emm2            = _mm512_and_si512(emm2,TWO);
    // Temporary workaround, shall be looking for more efficient implementation
    __m256i lo_emm2 = _mm512_extracti32x8_epi32(emm2,0);
    __m256i hi_emm2 = _mm512_extracti32x8_epi32(emm2,1);
    lo_emm2         = _mm256_cmpeq_epi32(lo_emm2,_mm256_setzero_si256());
    hi_emm2         = _mm256_cmpeq_epi32(hi_emm2,_mm256_setzero_si256());
    // Causing wrong values.
    //emm2_eq_zero = _mm512_cmp_epi32_mask(emm2,_mm512_setzero_si512(),_MM_CMPINT_EQ);
    //emm2        = _mm512_mask_mov_epi32(_mm512_setzero_si512(),emm2_eq_zero,ONE);
    emm2            = _mm512_inserti32x8(emm2,lo_emm2,0);
    emm2            = _mm512_inserti32x8(emm2,hi_emm2,1);
    __m512 sign     = _mm512_castsi512_ps(emm0);
    __m512 mask     = _mm512_castsi512_ps(emm2);
    zmm1            = _mm512_mul_ps(y,DP1N);
    x               = _mm512_add_ps(x,zmm1);
    zmm2            = _mm512_mul_ps(y,DP2N);
    x               = _mm512_add_ps(x,zmm2);
    zmm3            = _mm512_mul_ps(y,DP3N);
    x               = _mm512_add_ps(x,zmm3);
    __m512 z        = _mm512_mul_ps(x,x);
    // Long dependency chain which may stall pipeline on the 'y' register interdependent result computation
    // Shall be considering the Estrin's scheme (if possible).
    y               = C2443315711809948E005;
    y               = _mm512_mul_ps(y,z);
    y               = _mm512_add_ps(y,CN1388731625493765E003);
    y               = _mm512_mul_ps(y,z);
    y               = _mm512_add_ps(y,C4166664568298827E002);
    y               = _mm512_mul_ps(y,z);
    y               = _mm512_mul_ps(y,z);
    __m512 tmp      = _mm512_mul_ps(z,C05);
    y               = _mm512_sub_ps(y,tmp);
    y               = _mm512_add_ps(y,C1);
    __m512 y2       = CN19515295891E4;
    // Long dependency chain which may stall pipeline on the 'y2' register interdependent result computation
    // Shall be considering the Estrin's scheme (if possible).
    y2              = _mm512_mul_ps(y2,z);
    y2              = _mm512_add_ps(y2,C83321608736E3);
    y2              = _mm512_mul_ps(y2,z);
    y2              = _mm512_add_ps(y2,CN16666654611E1);
    y2              = _mm512_mul_ps(y2,z);
    y2              = _mm512_mul_ps(y2,x);
    y2              = _mm512_add_ps(y2,x);
    zmm3            = mask;
    y2              = _mm512_and_ps(zmm3,y2);
    y               = _mm512_andnot_ps(zmm3,y);
    y               = _mm512_add_ps(y,y2);
    y               = _mm512_xor_ps(y,sign);
    return (y);
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching>
inline static __m512 _mm512_ceph_cosf_ps_v3(__m512 x)
{
   __m512i inv_sign_mask;
   __m512i sign_mask;
   __m512i ONE;
   __m512i ONEINV;
   __m512i TWO;
   __m512i FOUR;
   __m512  FOPI;
   __m512  DP1N;
	__m512  DP2N;
	__m512  DP3N;
   __m512  C2443315711809948E005;
   __m512  CN1388731625493765E003;
   __m512  C4166664568298827E002;
   __m512  C05;
   __m512  C1;
   __m512  CN19515295891E4;
   __m512  C83321608736E3;
   __m512  CN16666654611E1;
   __m512  zmm1;
   __m512  zmm2;
   __m512  zmm3;
   __m512  y;
   __m512i emm0;
   __m512i emm2;
   const float *        __restrict__ ptr_prefetched_constants_16x_f32 = &simd_cos_constants::prefetched_constants_16x_f32[0];
	const std::int32_t * __restrict__ ptr_prefetched_constants_16x_i32 = &simd_cos_constants::prefetched_constants_16x_i32[0];
   const std::int32_t i_use_prefetching = static_cast<std::int32_t>(use_prefetching);
   const std::int32_t i_true            = static_cast<std::int32_t>(true);
   const bool cond                      = (i_use_prefetching==i_true);
   if constexpr (cond) {_mm_prefetch((const char*)&ptr_prefetched_constants_16x_i32[0],_MM_HINT_T0);}
   inv_sign_mask    = _mm512_load_epi32((const std::int32_t*)&ptr_prefetched_constants_16x_i32[0]);
   x                = _mm512_and_ps(x,*(__m512*)&inv_sign_mask);
   if constexpr (cond) {_mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[0],_MM_HINT_T0);}
   FOPI             = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[0]);
   y                = _mm512_mul_ps(x,FOPI);
   if constexpr (cond) {_mm_prefetch((const char*)&ptr_prefetched_constants_16x_i32[32],_MM_HINT_T0);}
   ONE              = _mm512_load_epi32((const std::int32_t*)&ptr_prefetched_constants_16x_i32[32]);
   emm2             = _mm512_cvttps_epi32(y);
   if constexpr (cond) {_mm_prefetch((const char*)&ptr_prefetched_constants_16x_i32[48],_MM_HINT_T0);}
   ONEINV           = _mm512_load_epi32((const std::int32_t*)&ptr_prefetched_constants_16x_i32[48]);
   emm2             = _mm512_add_epi32(emm2,ONE);
   emm2             = _mm512_and_si512(emm2,ONEINV);
   if constexpr (cond) {_mm_prefetch((const char*)&ptr_prefetched_constants_16x_i32[64],_MM_HINT_T0);}
   TWO              = _mm512_load_epi32((const std::int32_t*)&ptr_prefetched_constants_16x_i32[64]);
   y                = _mm512_cvtepi32_ps(emm2);
   if constexpr (cond) { _mm_prefetch((const char*)&ptr_prefetched_constants_16x_i32[80],_MM_HINT_T0);}
   FOUR             = _mm512_load_epi32((const std::int32_t*)&ptr_prefetched_constants_16x_i32[80]);
   emm2             = _mm512_sub_epi32(emm2,TWO);
   emm0             = _mm512_andnot_si512(emm2,FOUR);
   emm0             = _mm512_slli_epi32(emm0,29);
   emm2             = _mm512_and_si512(emm2,TWO);
   // Temporary workaround, shall be looking for more efficient implementation
   if constexpr (cond) {_mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[16],_MM_HINT_T0);}
   DP1N            = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[16]);
   __m256i lo_emm2 = _mm512_extracti32x8_epi32(emm2,0);
   __m256i hi_emm2 = _mm512_extracti32x8_epi32(emm2,1);
   if constexpr (cond) {_mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[32],_MM_HINT_T0);}
   DP2N            = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[32]);
   lo_emm2         = _mm256_cmpeq_epi32(lo_emm2,_mm256_setzero_si256());
   hi_emm2         = _mm256_cmpeq_epi32(hi_emm2,_mm256_setzero_si256());
    // Causing wrong values.
    //emm2_eq_zero = _mm512_cmp_epi32_mask(emm2,_mm512_setzero_si512(),_MM_CMPINT_EQ);
    //emm2        = _mm512_mask_mov_epi32(_mm512_setzero_si512(),emm2_eq_zero,ONE);
   if constexpr (cond) {_mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[48],_MM_HINT_T0);}
   DP3N            = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[48]);
   emm2            = _mm512_inserti32x8(emm2,lo_emm2,0);
   emm2            = _mm512_inserti32x8(emm2,hi_emm2,1);
   if constexpr (cond) {_mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[64],_MM_HINT_T0);}
   C2443315711809948E005  = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[64]);
   __m512 sign     = _mm512_castsi512_ps(emm0);
   __m512 mask     = _mm512_castsi512_ps(emm2);
   if constexpr (cond) {_mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[80],_MM_HINT_T0);}
   CN1388731625493765E003 = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[80]);
   zmm1            = _mm512_mul_ps(y,DP1N);
   x               = _mm512_add_ps(x,zmm1);
   if constexpr (cond) {_mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[96],_MM_HINT_T0);}
   C4166664568298827E002  = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[96]);
   zmm2            = _mm512_mul_ps(y,DP2N);
   x               = _mm512_add_ps(x,zmm2);
   if constexpr (cond) {_mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[112],_MM_HINT_T0);}
   C05             = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[112]);
   zmm3            = _mm512_mul_ps(y,DP3N);
   x               = _mm512_add_ps(x,zmm3);
   if constexpr (cond) {_mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[128],_MM_HINT_T0);}
   C1              = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[128]);
   __m512 z        = _mm512_mul_ps(x,x);
   // Long dependency chain which may stall pipeline on the 'y' register interdependent result computation
   // Shall be considering the Estrin's scheme (if possible).
   y               = C2443315711809948E005;
   y               = _mm512_mul_ps(y,z);
   y               = _mm512_add_ps(y,CN1388731625493765E003);
   y               = _mm512_mul_ps(y,z);
   y               = _mm512_add_ps(y,C4166664568298827E002);
   if constexpr (cond) {_mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[144],_MM_HINT_T0);}
   CN19515295891E4 = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[144]);
   y               = _mm512_mul_ps(y,z);
   y               = _mm512_mul_ps(y,z);
   __m512 tmp      = _mm512_mul_ps(z,C05);
   if constexpr (cond) {_mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[160],_MM_HINT_T0);}
   C83321608736E3  = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[160]);
   y               = _mm512_sub_ps(y,tmp);
   y               = _mm512_add_ps(y,C1);
   __m512 y2       = CN19515295891E4;
   if constexpr (cond) {_mm_prefetch((const char*)&ptr_prefetched_constants_16x_f32[176],_MM_HINT_T0);}
   CN16666654611E1 = _mm512_load_ps((const float*)&ptr_prefetched_constants_16x_f32[176]);
    // Long dependency chain which may stall pipeline on the 'y2' register interdependent result computation
    // Shall be considering the Estrin's scheme (if possible).
   y2              = _mm512_mul_ps(y2,z);
   y2              = _mm512_add_ps(y2,C83321608736E3);
   y2              = _mm512_mul_ps(y2,z);
   y2              = _mm512_add_ps(y2,CN16666654611E1);
   y2              = _mm512_mul_ps(y2,z);
   y2              = _mm512_mul_ps(y2,x);
   y2              = _mm512_add_ps(y2,x);
   zmm3            = mask;
   y2              = _mm512_and_ps(zmm3,y2);
   y               = _mm512_andnot_ps(zmm3,y);
   y               = _mm512_add_ps(y,y2);
   y               = _mm512_xor_ps(y,sign);
   return (y);
} 


} // math

}// gms









#endif /*__GMS_SIMD_CEPHES_COS_H*/