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

/* Copyright (C) Bernard Gingold, 2020-2026 
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

#ifndef __GMS_SIMD_CEPHES_SIN_H__
#define __GMS_SIMD_CEPHES_SIN_H__ 040820260616

namespace file_info 
{

     static const unsigned int GMS_SIMD_CEPHES_SIN_MAJOR = 1;
     static const unsigned int GMS_SIMD_CEPHES_SIN_MINOR = 0;
     static const unsigned int GMS_SIMD_CEPHES_SIN_MICRO = 0;
     static const unsigned int GMS_SIMD_CEPHES_SIN_FULLVER =
       1000U*GMS_SIMD_CEPHES_SIN_MAJOR+100U*GMS_SIMD_CEPHES_SIN_MINOR+
       10U*GMS_SIMD_CEPHES_SIN_MICRO;
     static const char GMS_SIMD_CEPHES_SIN_CREATION_DATE[] = "04-08-2026 06:18AM +00200 (TUE 04 AUG 2026 GMT+2)";
     static const char GMS_SIM_CEPHES_SIN_BUILD_DATE[]     = __DATE__; 
     static const char GMS_SIMD_CEPHES_SIN_BUILD_TIME[]    = __TIME__;
     static const char GMS_SIMD_CEPHES_SIN_SYNOPSIS[]      = "Manual vectorization of CEPHES sine(float) function.";

}

#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace gms
{

namespace math 
{

/* Much simplier implementation without the rodata prefetched storage.
   This was decided due to overcomplicated implementation by the Compiler
   which did insert an additional loads beside the usually expected pattern
   of the direct load-compute operands.*/

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
__m128 _mm_ceph_sinf_ps(__m128 x)
{
    __m128 xmm1;
    __m128 xmm2;
    __m128 xmm3;
    __m128 sign_bit, y;
	__m128i emm0, emm2;
	sign_bit = x;
	/* take the absolute value */
	const __m128 inv_sign_mask = _mm_castsi128_ps(_mm_set1_epi32(~0x80000000));
	x = _mm_and_ps(x, inv_sign_mask);
	/* extract the sign bit (upper one) */
	const __m128 sign_mask = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
	sign_bit = _mm_and_ps(sign_bit, sign_mask);

	/* scale by 4/Pi */
	const __m128 cephes_FOPI = _mm_set1_ps(1.27323954473516); // 4 / M_PI
	y = _mm_mul_ps(x, cephes_FOPI);

	/* store the integer part of y in mm0 */
	emm2 = _mm_cvttps_epi32(y);
	/* j=(j+1) & (~1) (see the cephes sources) */
	emm2 = _mm_add_epi32(emm2, _mm_set1_epi32(1));
	emm2 = _mm_and_si128(emm2, _mm_set1_epi32(~1));
	y = _mm_cvtepi32_ps(emm2);

	/* get the swap sign flag */
	emm0 = _mm_and_si128(emm2, _mm_set1_epi32(4));
	emm0 = _mm_slli_epi32(emm0, 29);
	/* get the polynom selection mask
	   there is one polynom for 0 <= x <= Pi/4
	   and another one for Pi/4<x<=Pi/2

	   Both branches will be computed.
	*/
	emm2 = _mm_and_si128(emm2, _mm_set1_epi32(2));
	emm2 = _mm_cmpeq_epi32(emm2, _mm_setzero_si128());

	__m128 swap_sign_bit = _mm_castsi128_ps(emm0);
	__m128 poly_mask = _mm_castsi128_ps(emm2);
	sign_bit = _mm_xor_ps(sign_bit, swap_sign_bit);

	/* The magic pass: "Extended precision modular arithmetic"
	   x = ((x - y * DP1) - y * DP2) - y * DP3; */
	xmm1 = _mm_set1_ps(-0.78515625);
	xmm2 = _mm_set1_ps(-2.4187564849853515625e-4);
	xmm3 = _mm_set1_ps(-3.77489497744594108e-8);
	xmm1 = _mm_mul_ps(y, xmm1);
	xmm2 = _mm_mul_ps(y, xmm2);
	xmm3 = _mm_mul_ps(y, xmm3);
	x = _mm_add_ps(x, xmm1);
	x = _mm_add_ps(x, xmm2);
	x = _mm_add_ps(x, xmm3);

	/* Evaluate the first polynom  (0 <= x <= Pi/4) */
	y = _mm_set1_ps(2.443315711809948E-005);
	__m128 z = _mm_mul_ps(x, x);

	y = _mm_mul_ps(y, z);
	y = _mm_add_ps(y, _mm_set1_ps(-1.388731625493765E-003));
	y = _mm_mul_ps(y, z);
	y = _mm_add_ps(y, _mm_set1_ps(4.166664568298827E-002));
	y = _mm_mul_ps(y, z);
	y = _mm_mul_ps(y, z);
	__m128 tmp = _mm_mul_ps(z, _mm_set1_ps(0.5));
	y = _mm_sub_ps(y, tmp);
	y = _mm_add_ps(y, _mm_set1_ps(1.0));

	/* Evaluate the second polynom  (Pi/4 <= x <= 0) */

	__m128 y2 = _mm_set1_ps(-1.9515295891E-4);
	y2 = _mm_mul_ps(y2, z);
	y2 = _mm_add_ps(y2, _mm_set1_ps(8.3321608736E-3));
	y2 = _mm_mul_ps(y2, z);
	y2 = _mm_add_ps(y2, _mm_set1_ps(-1.6666654611E-1));
	y2 = _mm_mul_ps(y2, z);
	y2 = _mm_mul_ps(y2, x);
	y2 = _mm_add_ps(y2, x);

	/* select the correct result from the two polynoms */
	xmm3 = poly_mask;
	y2 = _mm_and_ps(xmm3, y2); //, xmm3);
	y = _mm_andnot_ps(xmm3, y);
	y = _mm_add_ps(y, y2);
	/* update the sign */
	y = _mm_xor_ps(y, sign_bit);
	return y;
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
__m256 _mm256_ceph_sinf_ps(__m256 x)
{
    __m256 xmm1;
    __m256 xmm2;
    __m256 xmm3;
    __m256 sign_bit, y;
	__m256i emm0, emm2;
	sign_bit = x;
	/* take the absolute value */
	const __m256 inv_sign_mask = _mm256_castsi256_ps(_mm256_set1_epi32(~0x80000000));
	x = _mm256_and_ps(x, inv_sign_mask);
	/* extract the sign bit (upper one) */
	const __m256 sign_mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));
	sign_bit = _mm256_and_ps(sign_bit, sign_mask);

	/* scale by 4/Pi */
	const __m256 cephes_FOPI = _mm256_set1_ps(1.27323954473516); // 4 / M_PI
	y = _mm256_mul_ps(x, cephes_FOPI);

	/* store the integer part of y in mm0 */
	emm2 = _mm256_cvttps_epi32(y);
	/* j=(j+1) & (~1) (see the cephes sources) */
	emm2 = _mm256_add_epi32(emm2, _mm256_set1_epi32(1));
	emm2 = _mm256_and_si256(emm2, _mm256_set1_epi32(~1));
	y = _mm256_cvtepi32_ps(emm2);

	/* get the swap sign flag */
	emm0 = _mm256_and_si256(emm2, _mm256_set1_epi32(4));
	emm0 = _mm256_slli_epi32(emm0, 29);
	/* get the polynom selection mask
	   there is one polynom for 0 <= x <= Pi/4
	   and another one for Pi/4<x<=Pi/2

	   Both branches will be computed.
	*/
	emm2 = _mm256_and_si256(emm2, _mm256_set1_epi32(2));
	emm2 = _mm256_cmpeq_epi32(emm2, _mm256_setzero_si256());

	__m256 swap_sign_bit = _mm256_castsi256_ps(emm0);
	__m256 poly_mask = _mm256_castsi256_ps(emm2);
	sign_bit = _mm256_xor_ps(sign_bit, swap_sign_bit);

	/* The magic pass: "Extended precision modular arithmetic"
	   x = ((x - y * DP1) - y * DP2) - y * DP3; */
	xmm1 = _mm256_set1_ps(-0.78515625);
	xmm2 = _mm256_set1_ps(-2.4187564849853515625e-4);
	xmm3 = _mm256_set1_ps(-3.77489497744594108e-8);
	xmm1 = _mm256_mul_ps(y, xmm1);
	xmm2 = _mm256_mul_ps(y, xmm2);
	xmm3 = _mm256_mul_ps(y, xmm3);
	x = _mm256_add_ps(x, xmm1);
	x = _mm256_add_ps(x, xmm2);
	x = _mm256_add_ps(x, xmm3);

	/* Evaluate the first polynom  (0 <= x <= Pi/4) */
	y = _mm256_set1_ps(2.443315711809948E-005);
	__m256 z = _mm256_mul_ps(x, x);

	y = _mm256_mul_ps(y, z);
	y = _mm256_add_ps(y, _mm256_set1_ps(-1.388731625493765E-003));
	y = _mm256_mul_ps(y, z);
	y = _mm256_add_ps(y, _mm256_set1_ps(4.166664568298827E-002));
	y = _mm256_mul_ps(y, z);
	y = _mm256_mul_ps(y, z);
	__m256 tmp = _mm256_mul_ps(z, _mm256_set1_ps(0.5));
	y = _mm256_sub_ps(y, tmp);
	y = _mm256_add_ps(y, _mm256_set1_ps(1.0));

	/* Evaluate the second polynom  (Pi/4 <= x <= 0) */

	__m256 y2 = _mm256_set1_ps(-1.9515295891E-4);
	y2 = _mm256_mul_ps(y2, z);
	y2 = _mm256_add_ps(y2, _mm256_set1_ps(8.3321608736E-3));
	y2 = _mm256_mul_ps(y2, z);
	y2 = _mm256_add_ps(y2, _mm256_set1_ps(-1.6666654611E-1));
	y2 = _mm256_mul_ps(y2, z);
	y2 = _mm256_mul_ps(y2, x);
	y2 = _mm256_add_ps(y2, x);

	/* select the correct result from the two polynoms */
	xmm3 = poly_mask;
	y2 = _mm256_and_ps(xmm3, y2); //, xmm3);
	y = _mm256_andnot_ps(xmm3, y);
	y = _mm256_add_ps(y, y2);
	/* update the sign */
	y = _mm256_xor_ps(y, sign_bit);
	return y;
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
__m512 _mm512_ceph_sinf_ps(__m512 x)
{
    __m512 zmm1;
    __m512 zmm2;
    __m512 zmm3;
    __m512 sign_bit, y;
	__m512i emm0, emm2;
	sign_bit = x;
	/* take the absolute value */
	const __m512 inv_sign_mask = _mm512_castsi512_ps(_mm512_set1_epi32(~0x80000000));
	x = _mm512_and_ps(x, inv_sign_mask);
	/* extract the sign bit (upper one) */
	const __m512 sign_mask = _mm512_castsi512_ps(_mm512_set1_epi32(0x80000000));
	sign_bit = _mm512_and_ps(sign_bit, sign_mask);

	/* scale by 4/Pi */
	const __m512 cephes_FOPI = _mm512_set1_ps(1.27323954473516); // 4 / M_PI
	y = _mm512_mul_ps(x, cephes_FOPI);

	/* store the integer part of y in mm0 */
	emm2 = _mm512_cvttps_epi32(y);
	/* j=(j+1) & (~1) (see the cephes sources) */
	emm2 = _mm512_add_epi32(emm2, _mm512_set1_epi32(1));
	emm2 = _mm512_and_si512(emm2, _mm512_set1_epi32(~1));
	y = _mm512_cvtepi32_ps(emm2);

	/* get the swap sign flag */
	emm0 = _mm512_and_si512(emm2, _mm512_set1_epi32(4));
	emm0 = _mm512_slli_epi32(emm0, 29);
	/* get the polynom selection mask
	   there is one polynom for 0 <= x <= Pi/4
	   and another one for Pi/4<x<=Pi/2

	   Both branches will be computed.
	*/
	emm2 = _mm512_and_si512(emm2, _mm512_set1_epi32(2));
	//emm2 = _mm256_cmpeq_epi32(emm2, _mm256_setzero_si256());
    __m256i lo_emm2 = _mm512_extracti32x8_epi32(emm2,0);
    __m256i hi_emm2 = _mm512_extracti32x8_epi32(emm2,1);
    lo_emm2         = _mm256_cmpeq_epi32(lo_emm2,_mm256_setzero_si256());
    hi_emm2         = _mm256_cmpeq_epi32(hi_emm2,_mm256_setzero_si256());
    // Causing wrong values.
    //emm2_eq_zero = _mm512_cmp_epi32_mask(emm2,_mm512_setzero_si512(),_MM_CMPINT_EQ);
    //emm2        = _mm512_mask_mov_epi32(_mm512_setzero_si512(),emm2_eq_zero,ONE);
    emm2            = _mm512_inserti32x8(emm2,lo_emm2,0);
    emm2            = _mm512_inserti32x8(emm2,hi_emm2,1);
	__m512 swap_sign_bit = _mm512_castsi512_ps(emm0);
    zmm1 = _mm512_set1_ps(-0.78515625);
	__m512 poly_mask = _mm512_castsi512_ps(emm2);
    zmm2 = _mm512_set1_ps(-2.4187564849853515625e-4);
	sign_bit = _mm512_xor_ps(sign_bit, swap_sign_bit);

	/* The magic pass: "Extended precision modular arithmetic"
	   x = ((x - y * DP1) - y * DP2) - y * DP3; */
	zmm3 = _mm512_set1_ps(-3.77489497744594108e-8);
	zmm1 = _mm512_mul_ps(y, zmm1);
	zmm2 = _mm512_mul_ps(y, zmm2);
	zmm3 = _mm512_mul_ps(y, zmm3);
	x = _mm512_add_ps(x, zmm1);
	x = _mm512_add_ps(x, zmm2);
	x = _mm512_add_ps(x, zmm3);

	/* Evaluate the first polynom  (0 <= x <= Pi/4) */
	y = _mm512_set1_ps(2.443315711809948E-005);
	__m512 z = _mm512_mul_ps(x, x);

	y = _mm512_mul_ps(y, z);
	y = _mm512_add_ps(y, _mm512_set1_ps(-1.388731625493765E-003));
	y = _mm512_mul_ps(y, z);
	y = _mm512_add_ps(y, _mm512_set1_ps(4.166664568298827E-002));
	y = _mm512_mul_ps(y, z);
	y = _mm512_mul_ps(y, z);
	__m512 tmp = _mm512_mul_ps(z, _mm512_set1_ps(0.5));
	y = _mm512_sub_ps(y, tmp);
	y = _mm512_add_ps(y, _mm512_set1_ps(1.0));

	/* Evaluate the second polynom  (Pi/4 <= x <= 0) */

	__m512 y2 = _mm512_set1_ps(-1.9515295891E-4);
	y2 = _mm512_mul_ps(y2, z);
	y2 = _mm512_add_ps(y2, _mm512_set1_ps(8.3321608736E-3));
	y2 = _mm512_mul_ps(y2, z);
	y2 = _mm512_add_ps(y2, _mm512_set1_ps(-1.6666654611E-1));
	y2 = _mm512_mul_ps(y2, z);
	y2 = _mm512_mul_ps(y2, x);
	y2 = _mm512_add_ps(y2, x);

	/* select the correct result from the two polynoms */
	zmm3 = poly_mask;
	y2 = _mm512_and_ps(zmm3, y2); //, xmm3);
	y = _mm512_andnot_ps(zmm3, y);
	y = _mm512_add_ps(y, y2);
	/* update the sign */
	y = _mm512_xor_ps(y, sign_bit);
	return y;
}

} // math

} //gms

#endif /*__GMS_SIMD_CEPHES_SIN_H__*/