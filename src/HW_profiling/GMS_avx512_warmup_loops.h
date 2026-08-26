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
#ifndef __GMS_AVX512_WARMUP_LOOPS_H__
#define __GMS_AVX512_WARMUP_LOOPS_H__


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


__ATTR_ALWAYS_INLINE__
static inline
int32_t avx512_warmup_loop1_ps(const float * volatile __restrict in_a,
                               const float * volatile __restrict in_b,
			       const float * volatile __restrict in_c,
			       float * volatile __restrict out,
			       const int32_t len) {
        if(len <= 0) {
           return (-1);
	}
        int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(in_a,64);
	__assume_aligned(in_b,64);
	__assume_aligned(in_c,64);
	__assume_aligned(out,64);
#pragma loop_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        in_a = (const float*)__builtin_assume_aligned(in_a,64);
	in_b = (const float*)__builtin_assume_aligned(in_b,64);
	in_c = (const float*)__builtin_assume_aligned(in_c,64);
	out  = (float*)__builtin_assume_aligned(out,64);
#endif
        for(i = 0; i != ROUND_TO_SIXTEEN(len,15)-9; i += 160) {
            const __m512 zmm0 = _mm512_loadu_ps(&in_a[i+0]);
	    const __m512 zmm1 = _mm512_loadu_ps(&in_b[i+0]);
            const __m512 zmm2 = _mm512_loadu_ps(&in_c[i+0]);
	    _mm512_storeu_ps(&out[i+0],_mm512_fmadd_ps(zmm0,zmm1,zmm2));
	    const __m512 zmm3 = _mm512_loadu_ps(&in_a[i+16]);
	    const __m512 zmm4 = _mm512_loadu_ps(&in_b[i+16]);
	    const __m512 zmm5 = _mm512_loadu_ps(&in_c[i+16]);
	    _mm512_storeu_ps(&out[i+16],_mm512_fmadd_ps(zmm3,zmm4,zmm5));
	    const __m512 zmm6 = _mm512_loadu_ps(&in_a[i+32]);
	    const __m512 zmm7 = _mm512_loadu_ps(&in_b[i+32]);
	    const __m512 zmm8 = _mm512_loadu_ps(&in_c[i+32]);
	    _mm512_storeu_ps(&out[i+32],_mm512_fmadd_ps(zmm6,zmm7,zmm8));
	    const __m512 zmm9  = _mm512_loadu_ps(&in_a[i+48]);
	    const __m512 zmm10 = _mm512_loadu_ps(&in_b[i+48]);
	    const __m512 zmm11 = _mm512_loadu_ps(&in_c[i+48]);
	    _mm512_storeu_ps(&out[i+48],_mm512_fmadd_ps(zmm9,zmm10,zmm11));
	    const __m512 zmm12 = _mm512_loadu_ps(&in_a[i+64]);
	    const __m512 zmm13 = _mm512_loadu_ps(&in_b[i+64]);
	    const __m512 zmm14 = _mm512_loadu_ps(&in_c[i+64]);
	    _mm512_storeu_ps(&out[i+64],_mm512_fmadd_ps(zmm12,zmm13,zmm14));
	    const __m512 zmm15 = _mm512_loadu_ps(&in_a[i+80]);
	    const __m512 zmm16 = _mm512_loadu_ps(&in_b[i+80]);
	    const __m512 zmm17 = _mm512_loadu_ps(&in_c[i+80]);
	    _mm512_storeu_ps(&out[i+80],_mm512_fmadd_ps(zmm15,zmm16,zmm17));
	    const __m512 zmm18 = _mm512_loadu_ps(&in_a[i+96]);
	    const __m512 zmm19 = _mm512_loadu_ps(&in_b[i+96]);
	    const __m512 zmm20 = _mm512_loadu_ps(&in_c[i+96]);
	    _mm512_storeu_ps(&out[i+96],_mm512_fmadd_ps(zmm18,zmm19,zmm20));
	    const __m512 zmm21 = _mm512_loadu_ps(&in_a[i+112]);
	    const __m512 zmm22 = _mm512_loadu_ps(&in_b[i+112]);
	    const __m512 zmm23 = _mm512_loadu_ps(&in_c[i+112]);
	    _mm512_storeu_ps(&out[i+112],_mm512_fmadd_ps(zmm21,zmm22,zmm23));
	    const __m512 zmm24 = _mm512_loadu_ps(&in_a[i+128]);
	    const __m512 zmm25 = _mm512_loadu_ps(&in_b[i+128]);
	    const __m512 zmm26 = _mm512_loadu_ps(&in_c[i+128]);
	    _mm512_storeu_ps(&out[i+128],_mm512_fmadd_ps(zmm24,zmm25,zmm26));
	    const __m512 zmm27 = _mm512_loadu_ps(&in_a[i+144]);
	    const __m512 zmm28 = _mm512_loadu_ps(&in_b[i+144]);
	    const __m512 zmm29 = _mm512_loadu_ps(&in_c[i+144]);
	    _mm512_storeu_ps(&out[i+144],_mm512_fmadd_ps(zmm27,zmm28,zmm29));
	}
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(8),max(15)
#endif
        for(; i != len; ++i) {
            out[i] = in_a[i]*in_b[i]+in_c[i];
	}
	return (0);
}


// THis function does not have a concrete meaning,
// its task is to trigger a L2 license only.
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
__m512 avx512_warmup_loop2_ps(const __m512 a,
                              const __m512 b,
			      const __m512 c){
	      
       __m512 volatile result = _mm512_setzero_ps();
#pragma nounroll
#pragma loop_count(1000000)
        for(int32_t i = 0; i != 1000000; ++i) {
            result = _mm512_fmadd_ps(a,b,c);
	}
	return (result);
}



#endif /*__GMS_AVX512_WARMUP_LOOPS*/
