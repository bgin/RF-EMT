




/*
     This version is based on BLIS framework original implementation.
     Few changes were introduced in order to more easily adapt it 
     to 'RF-EMT' project.
     
   The original authors license
   BLIS
   An object-based framework for developing high-performance BLAS-like
   libraries.
   Copyright (C) 2016 - 2019, Advanced Micro Devices, Inc.
   Copyright (C) 2018 - 2020, The University of Texas at Austin. All rights reserved.
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:
    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    - Neither the name(s) of the copyright holder(s) nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
     
*/


#include <immintrin.h>
#include <omp.h>
#include "GMS_copy_avx_unroll16x.h"





		   void gms::math::scopy_u_ymm8r8_unroll16x(const int32_t n,
		                                 float * __restrict x,
						 const int32_t incx,
						 float * __restrict y,
						 const int32_t incy) {
                        //Early exit
		        if(__builtin_expect(0==n,0)) { return;}
			__ATTR_ALIGN__(32) __m256 xv[16];
			int32_t i = 0;
			
			if(__builtin_expect(1==incx,1) &&
			   __builtin_expect(1==incy,1)) {

                             for(i = 0; (i+127) < n; i += 128) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
				 _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
                                 xv[0] = _mm256_loadu_ps(&x[i+0]);
				 _mm256_storeu_ps(&y[i+0],xv[0]);
				 xv[1] = _mm256_loadu_ps(&x[i+8]);
				 _mm256_storeu_ps(&y[i+8],xv[1]);
				 xv[2] = _mm256_loadu_ps(&x[i+16]);
				 _mm256_storeu_ps(&y[i+16],xv[2]);
				 xv[3] = _mm256_loadu_ps(&x[i+24]);
				 _mm256_storeu_ps(&y[i+24],xv[3]);
				 xv[4] = _mm256_loadu_ps(&x[i+32]);
				 _mm256_storeu_ps(&y[i+32],xv[4]);
				 _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				 xv[5] = _mm256_loadu_ps(&x[i+40]);
				 _mm256_storeu_ps(&y[i+40],xv[5]);
				 xv[6] = _mm256_loadu_ps(&x[i+48]);
				 _mm256_storeu_ps(&y[i+48],xv[6]);
				 xv[7] = _mm256_loadu_ps(&x[i+56]);
				 _mm256_storeu_ps(&y[i+56],xv[7]);
				 xv[8] = _mm256_loadu_ps(&x[i+64]);
				 _mm256_storeu_ps(&y[i+64],xv[8]);
				 xv[9] = _mm256_loadu_ps(&x[i+72]);
				 _mm256_storeu_ps(&y[i+72],xv[9]);
				 _mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				 xv[10] = _mm256_loadu_ps(&x[i+80]);
				 _mm256_storeu_ps(&y[i+80],xv[10]);
				 xv[11] = _mm256_loadu_ps(&x[i+88]);
				 _mm256_storeu_ps(&y[i+88],xv[11]);
				  xv[12] = _mm256_loadu_ps(&x[i+96]);
				 _mm256_storeu_ps(&y[i+96],xv[12]);
				 xv[13] = _mm256_loadu_ps(&x[i+104]);
				 _mm256_storeu_ps(&y[i+104],xv[13]);
				 xv[14] = _mm256_loadu_ps(&x[i+112]);
				 _mm256_storeu_ps(&y[i+112],xv[14]);
				 xv[15] = _mm256_loadu_ps(&x[i+120]);
				 _mm256_storeu_ps(&y[i+120],xv[15]);
#else
                                 _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				 xv[0]  = _mm256_loadu_ps(&x[i+0]);
				 xv[1]  = _mm256_loadu_ps(&x[i+8]);
				 xv[2]  = _mm256_loadu_ps(&x[i+16]);
				 xv[3]  = _mm256_loadu_ps(&x[i+24]);
				 xv[4]  = _mm256_loadu_ps(&x[i+32]);
				 _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				 xv[5]  = _mm256_loadu_ps(&x[i+40]);
				 xv[6]  = _mm256_loadu_ps(&x[i+48]);
				 xv[7]  = _mm256_loadu_ps(&x[i+56]);
				 xv[8]  = _mm256_loadu_ps(&x[i+64]);
				 xv[9]  = _mm256_loadu_ps(&x[i+72]);
				 _mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				 xv[10] = _mm256_loadu_ps(&x[i+80]);
				 xv[11] = _mm256_loadu_ps(&x[i+88]);
				 xv[12] = _mm256_loadu_ps(&x[i+96]);
				 xv[13] = _mm256_loadu_ps(&x[i+104]);
				 xv[14] = _mm256_loadu_ps(&x[i+112]);
				 xv[15] = _mm256_loadu_ps(&x[i+120]);

				 _mm256_storeu_ps(&y[i+0],xv[0]);
				 _mm256_storeu_ps(&y[i+8],xv[1]);
				 _mm256_storeu_ps(&y[i+16],xv[2]);
				 _mm256_storeu_ps(&y[i+24],xv[3]);
				 _mm256_storeu_ps(&y[i+32],xv[4]);
				 _mm256_storeu_ps(&y[i+40],xv[5]);
				 _mm256_storeu_ps(&y[i+48],xv[6]);
				 _mm256_storeu_ps(&y[i+56],xv[7]);
				 _mm256_storeu_ps(&y[i+64],xv[8]);
				 _mm256_storeu_ps(&y[i+72],xv[9]);
				 _mm256_storeu_ps(&y[i+80],xv[10]);
				 _mm256_storeu_ps(&y[i+88],xv[11]);
				 _mm256_storeu_ps(&y[i+96],xv[12]);
				 _mm256_storeu_ps(&y[i+104],xv[13]);
				 _mm256_storeu_ps(&y[i+112],xv[14]);
				 _mm256_storeu_ps(&y[i+120],xv[15]);

#endif
			     }

                            for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
			         xv[0] = _mm256_loadu_ps(&x[i+0]);
				 _mm256_storeu_ps(&y[i+0],xv[0]);
				 xv[1] = _mm256_loadu_ps(&x[i+8]);
				 _mm256_storeu_ps(&y[i+8],xv[1]);
				 xv[2] = _mm256_loadu_ps(&x[i+16]);
				 _mm256_storeu_ps(&y[i+16],xv[2]);
				 xv[3] = _mm256_loadu_ps(&x[i+24]);
				 _mm256_storeu_ps(&y[i+24],xv[3]);
				 xv[4] = _mm256_loadu_ps(&x[i+32]);
				 _mm256_storeu_ps(&y[i+32],xv[4]);
				// _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				 xv[5] = _mm256_loadu_ps(&x[i+40]);
				 _mm256_storeu_ps(&y[i+40],xv[5]);
				 xv[6] = _mm256_loadu_ps(&x[i+48]);
				 _mm256_storeu_ps(&y[i+48],xv[6]);
				 xv[7] = _mm256_loadu_ps(&x[i+56]);
				 _mm256_storeu_ps(&y[i+56],xv[7]);
#else
                                // _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				 xv[0]  = _mm256_loadu_ps(&x[i+0]);
				 xv[1]  = _mm256_loadu_ps(&x[i+8]);
				 xv[2]  = _mm256_loadu_ps(&x[i+16]);
				 xv[3]  = _mm256_loadu_ps(&x[i+24]);
				 xv[4]  = _mm256_loadu_ps(&x[i+32]);
				// _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				 xv[5]  = _mm256_loadu_ps(&x[i+40]);
				 xv[6]  = _mm256_loadu_ps(&x[i+48]);
				 xv[7]  = _mm256_loadu_ps(&x[i+56]);

				 _mm256_storeu_ps(&y[i+0],xv[0]);
				 _mm256_storeu_ps(&y[i+8],xv[1]);
				 _mm256_storeu_ps(&y[i+16],xv[2]);
				 _mm256_storeu_ps(&y[i+24],xv[3]);
				 _mm256_storeu_ps(&y[i+32],xv[4]);
				 _mm256_storeu_ps(&y[i+40],xv[5]);
				 _mm256_storeu_ps(&y[i+48],xv[6]);
				 _mm256_storeu_ps(&y[i+56],xv[7]);

#endif
			    }

			    for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                 xv[0] = _mm256_loadu_ps(&x[i+0]);
				 _mm256_storeu_ps(&y[i+0],xv[0]);
				 xv[1] = _mm256_loadu_ps(&x[i+8]);
				 _mm256_storeu_ps(&y[i+8],xv[1]);
				 xv[2] = _mm256_loadu_ps(&x[i+16]);
				 _mm256_storeu_ps(&y[i+16],xv[2]);
				 xv[3] = _mm256_loadu_ps(&x[i+24]);
#else
                                 xv[0]  = _mm256_loadu_ps(&x[i+0]);
				 xv[1]  = _mm256_loadu_ps(&x[i+8]);
				 xv[2]  = _mm256_loadu_ps(&x[i+16]);
				 xv[3]  = _mm256_loadu_ps(&x[i+24]);
				 _mm256_storeu_ps(&y[i+0],xv[0]);
				 _mm256_storeu_ps(&y[i+8],xv[1]);
				 _mm256_storeu_ps(&y[i+16],xv[2]);
				 _mm256_storeu_ps(&y[i+24],xv[3]);
#endif

			    }

			    for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                 xv[0] = _mm256_loadu_ps(&x[i+0]);
				 _mm256_storeu_ps(&y[i+0],xv[0]);
				 xv[1] = _mm256_loadu_ps(&x[i+8]);
				 _mm256_storeu_ps(&y[i+8],xv[1]);
#else
                                 xv[0]  = _mm256_loadu_ps(&x[i+0]);
				 xv[1]  = _mm256_loadu_ps(&x[i+8]);
				 _mm256_storeu_ps(&y[i+0],xv[0]);
				 _mm256_storeu_ps(&y[i+8],xv[1]);
#endif

			    }

			    for(; (i+7) < n; i += 8) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                   xv[0] = _mm256_loadu_ps(&x[i+0]);
				   _mm256_storeu_ps(&y[i+0],xv[0]);
#else
                                   xv[0]  = _mm256_loadu_ps(&x[i+0]);
				   _mm256_storeu_ps(&y[i+0],xv[0]);
#endif

			    }

			    for(; (i+0) < n; i += 1) {
                                y[i] = x[i];
			    }
			    
			}
			else {

			    for(int32_t i = 0; i != n; ++i) {
                                *y = *x;
				 x += incx;
				 y += incy;
			    }

			}
		  }



		  /* COPY kernel single-precision aligned */
		
		   void gms::math::scopy_a_ymm8r8_unroll16x(const int32_t n,
		                                 float * __restrict __ATTR_ALIGN__(32) x,
						 const int32_t incx,
						 float * __restrict __ATTR_ALIGN__(32) y,
						 const int32_t incy) {
                        //Early exit
		        if(__builtin_expect(0==n,0)) { return;}
			__ATTR_ALIGN__(32) __m256 xv[16];
			int32_t i = 0;
			
			if(__builtin_expect(1==incx,1) &&
			   __builtin_expect(1==incy,1)) {
#if defined(__INTEL_COMPILER) || defined(__ICC)
                             __assume_aligned(x,32);
			     __assume_aligned(y,32);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                             x = (float*)__builtin_assume_aligned(x,32);
			     y = (float*)__builtin_assume_aligned(y,32);
#endif
                             for(i = 0; (i+127) < n; i += 128) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
				 _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
                                 xv[0] = _mm256_load_ps(&x[i+0]);
				 _mm256_store_ps(&y[i+0],xv[0]);
				 xv[1] = _mm256_load_ps(&x[i+8]);
				 _mm256_store_ps(&y[i+8],xv[1]);
				 xv[2] = _mm256_load_ps(&x[i+16]);
				 _mm256_store_ps(&y[i+16],xv[2]);
				 xv[3] = _mm256_load_ps(&x[i+24]);
				 _mm256_store_ps(&y[i+24],xv[3]);
				 xv[4] = _mm256_load_ps(&x[i+32]);
				 _mm256_store_ps(&y[i+32],xv[4]);
				 _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				 xv[5] = _mm256_load_ps(&x[i+40]);
				 _mm256_store_ps(&y[i+40],xv[5]);
				 xv[6] = _mm256_load_ps(&x[i+48]);
				 _mm256_store_ps(&y[i+48],xv[6]);
				 xv[7] = _mm256_loadu_ps(&x[i+56]);
				 _mm256_store_ps(&y[i+56],xv[7]);
				 xv[8] = _mm256_load_ps(&x[i+64]);
				 _mm256_store_ps(&y[i+64],xv[8]);
				 xv[9] = _mm256_load_ps(&x[i+72]);
				 _mm256_store_ps(&y[i+72],xv[9]);
				 _mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				 xv[10] = _mm256_load_ps(&x[i+80]);
				 _mm256_store_ps(&y[i+80],xv[10]);
				 xv[11] = _mm256_load_ps(&x[i+88]);
				 _mm256_store_ps(&y[i+88],xv[11]);
				  xv[12] = _mm256_load_ps(&x[i+96]);
				 _mm256_store_ps(&y[i+96],xv[12]);
				 xv[13] = _mm256_load_ps(&x[i+104]);
				 _mm256_store_ps(&y[i+104],xv[13]);
				 xv[14] = _mm256_load_ps(&x[i+112]);
				 _mm256_store_ps(&y[i+112],xv[14]);
				 xv[15] = _mm256_load_ps(&x[i+120]);
				 _mm256_store_ps(&y[i+120],xv[15]);
#else
                                 _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				 xv[0]  = _mm256_load_ps(&x[i+0]);
				 xv[1]  = _mm256_load_ps(&x[i+8]);
				 xv[2]  = _mm256_load_ps(&x[i+16]);
				 xv[3]  = _mm256_load_ps(&x[i+24]);
				 xv[4]  = _mm256_load_ps(&x[i+32]);
				 _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				 xv[5]  = _mm256_load_ps(&x[i+40]);
				 xv[6]  = _mm256_load_ps(&x[i+48]);
				 xv[7]  = _mm256_load_ps(&x[i+56]);
				 xv[8]  = _mm256_load_ps(&x[i+64]);
				 xv[9]  = _mm256_load_ps(&x[i+72]);
				 _mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				 xv[10] = _mm256_load_ps(&x[i+80]);
				 xv[11] = _mm256_load_ps(&x[i+88]);
				 xv[12] = _mm256_load_ps(&x[i+96]);
				 xv[13] = _mm256_load_ps(&x[i+104]);
				 xv[14] = _mm256_load_ps(&x[i+112]);
				 xv[15] = _mm256_load_ps(&x[i+120]);

				 _mm256_store_ps(&y[i+0],xv[0]);
				 _mm256_store_ps(&y[i+8],xv[1]);
				 _mm256_store_ps(&y[i+16],xv[2]);
				 _mm256_storeu_ps(&y[i+24],xv[3]);
				 _mm256_store_ps(&y[i+32],xv[4]);
				 _mm256_store_ps(&y[i+40],xv[5]);
				 _mm256_store_ps(&y[i+48],xv[6]);
				 _mm256_store_ps(&y[i+56],xv[7]);
				 _mm256_store_ps(&y[i+64],xv[8]);
				 _mm256_store_ps(&y[i+72],xv[9]);
				 _mm256_store_ps(&y[i+80],xv[10]);
				 _mm256_store_ps(&y[i+88],xv[11]);
				 _mm256_store_ps(&y[i+96],xv[12]);
				 _mm256_store_ps(&y[i+104],xv[13]);
				 _mm256_store_ps(&y[i+112],xv[14]);
				 _mm256_store_ps(&y[i+120],xv[15]);

#endif
			     }

                            for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
			         xv[0] = _mm256_load_ps(&x[i+0]);
				 _mm256_store_ps(&y[i+0],xv[0]);
				 xv[1] = _mm256_load_ps(&x[i+8]);
				 _mm256_store_ps(&y[i+8],xv[1]);
				 xv[2] = _mm256_load_ps(&x[i+16]);
				 _mm256_store_ps(&y[i+16],xv[2]);
				 xv[3] = _mm256_load_ps(&x[i+24]);
				 _mm256_store_ps(&y[i+24],xv[3]);
				 xv[4] = _mm256_load_ps(&x[i+32]);
				 _mm256_store_ps(&y[i+32],xv[4]);
				// _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				 xv[5] = _mm256_loadu_ps(&x[i+40]);
				 _mm256_store_ps(&y[i+40],xv[5]);
				 xv[6] = _mm256_load_ps(&x[i+48]);
				 _mm256_store_ps(&y[i+48],xv[6]);
				 xv[7] = _mm256_load_ps(&x[i+56]);
				 _mm256_store_ps(&y[i+56],xv[7]);
#else
                                // _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				 xv[0]  = _mm256_load_ps(&x[i+0]);
				 xv[1]  = _mm256_load_ps(&x[i+8]);
				 xv[2]  = _mm256_load_ps(&x[i+16]);
				 xv[3]  = _mm256_load_ps(&x[i+24]);
				 xv[4]  = _mm256_load_ps(&x[i+32]);
				// _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				 xv[5]  = _mm256_load_ps(&x[i+40]);
				 xv[6]  = _mm256_load_ps(&x[i+48]);
				 xv[7]  = _mm256_load_ps(&x[i+56]);

				 _mm256_store_ps(&y[i+0],xv[0]);
				 _mm256_store_ps(&y[i+8],xv[1]);
				 _mm256_store_ps(&y[i+16],xv[2]);
				 _mm256_store_ps(&y[i+24],xv[3]);
				 _mm256_store_ps(&y[i+32],xv[4]);
				 _mm256_store_ps(&y[i+40],xv[5]);
				 _mm256_store_ps(&y[i+48],xv[6]);
				 _mm256_store_ps(&y[i+56],xv[7]);

#endif
			    }

			    for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                 xv[0] = _mm256_load_ps(&x[i+0]);
				 _mm256_store_ps(&y[i+0],xv[0]);
				 xv[1] = _mm256_load_ps(&x[i+8]);
				 _mm256_store_ps(&y[i+8],xv[1]);
				 xv[2] = _mm256_load_ps(&x[i+16]);
				 _mm256_store_ps(&y[i+16],xv[2]);
				 xv[3] = _mm256_load_ps(&x[i+24]);
#else
                                 xv[0]  = _mm256_load_ps(&x[i+0]);
				 xv[1]  = _mm256_load_ps(&x[i+8]);
				 xv[2]  = _mm256_load_ps(&x[i+16]);
				 xv[3]  = _mm256_load_ps(&x[i+24]);
				 _mm256_store_ps(&y[i+0],xv[0]);
				 _mm256_store_ps(&y[i+8],xv[1]);
				 _mm256_store_ps(&y[i+16],xv[2]);
				 _mm256_store_ps(&y[i+24],xv[3]);
#endif

			    }

			    for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                 xv[0] = _mm256_load_ps(&x[i+0]);
				 _mm256_storeu_ps(&y[i+0],xv[0]);
				 xv[1] = _mm256_load_ps(&x[i+8]);
				 _mm256_store_ps(&y[i+8],xv[1]);
#else
                                 xv[0]  = _mm256_load_ps(&x[i+0]);
				 xv[1]  = _mm256_load_ps(&x[i+8]);
				 _mm256_store_ps(&y[i+0],xv[0]);
				 _mm256_store_ps(&y[i+8],xv[1]);
#endif

			    }

			    for(; (i+7) < n; i += 8) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                   xv[0] = _mm256_load_ps(&x[i+0]);
				   _mm256_store_ps(&y[i+0],xv[0]);
#else
                                   xv[0]  = _mm256_load_ps(&x[i+0]);
				   _mm256_store_ps(&y[i+0],xv[0]);
#endif

			    }

			    for(; (i+0) < n; i += 1) {
                                y[i] = x[i];
			    }
			    
			}
			else {

			    for(int32_t i = 0; i != n; ++i) {
                                *y = *x;
				 x += incx;
				 y += incy;
			    }

			}
		  }



		    /* COPY kernel single-precision aligned */
		  
		   void gms::math::scopy_a_ymm8r8_unroll16x_omp(const int32_t n,
		                                     float * __restrict __ATTR_ALIGN__(32) x,
						     const int32_t incx,
						     float * __restrict __ATTR_ALIGN__(32) y,
						     const int32_t incy) {
                        //Early exit
		        if(__builtin_expect(0==n,0)) { return;}
			
			int32_t i = 0;
			int32_t last_i = 0;
			
			if(__builtin_expect(1==incx,1) &&
			   __builtin_expect(1==incy,1)) {
#if defined(__INTEL_COMPILER) || defined(__ICC)
                             __assume_aligned(x,32);
			     __assume_aligned(y,32);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                             x = (float*)__builtin_assume_aligned(x,32);
			     y = (float*)__builtin_assume_aligned(y,32);
#endif
#pragma omp parallel for schedule(static,128) default(none) \
                             lastprivate(last_i) private(i) shared(n,x,y)
                             for(i = 0; (i+127) < n; i += 128) {
			         last_i = i;
				 _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				 _mm256_store_ps(&y[i+0],
				             _mm256_load_ps(&x[i+0]));
				 _mm256_store_ps(&y[i+8],
				             _mm256_load_ps(&x[i+8]));
				 _mm256_store_ps(&y[i+16],
				             _mm256_load_ps(&x[i+16]));
				 _mm256_store_ps(&y[i+24],
				             _mm256_load_ps(&x[i+24]));
				 _mm256_store_ps(&y[i+32],
				             _mm256_load_ps(&x[i+32]));
				 _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				 _mm256_store_ps(&y[i+40],
				             _mm256_load_ps(&x[i+40]));
				 _mm256_store_ps(&y[i+48],
				             _mm256_load_ps(&x[i+48]));
				 _mm256_store_ps(&y[i+56],
				             _mm256_load_ps(&x[i+56]));
				 _mm256_store_ps(&y[i+64],
				             _mm256_load_ps(&x[i+64]));
				 _mm256_store_ps(&y[i+72],
				             _mm256_load_ps(&x[i+72]));
				 _mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				 _mm256_store_ps(&y[i+80],
				             _mm256_load_ps(&x[i+80]));
				 _mm256_store_ps(&y[i+88],
				             _mm256_load_ps(&x[i+88]));
				 _mm256_store_ps(&y[i+96],
				             _mm256_load_ps(&x[i+96]));
				 _mm256_store_ps(&y[i+104],
				            _mm256_load_ps(&x[i+104]));
				 _mm256_store_ps(&y[i+112],
				             _mm256_load_ps(&x[i+112]));
				 _mm256_store_ps(&y[i+120],
				             _mm256_load_ps(&x[i+120]));
			                              
			     }

                            for(; (last_i+63) < n; last_i += 64) {

			       	// _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				 _mm256_store_ps(&y[last_i+0],
				             _mm256_load_ps(&x[last_i+0]));
				 _mm256_store_ps(&y[last_i+8],
				             _mm256_load_ps(&x[last_i+8]));
				 _mm256_store_ps(&y[last_i+16],
				             _mm256_load_ps(&x[last_i+16]));
				 _mm256_store_ps(&y[last_i+24],
				             _mm256_load_ps(&x[last_i+24]));
				 _mm256_store_ps(&y[last_i+32],
				             _mm256_load_ps(&x[last_i+32]));
				// _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				 _mm256_store_ps(&y[last_i+40],
				             _mm256_load_ps(&x[last_i+40]));
				 _mm256_store_ps(&y[last_i+48],
				             _mm256_load_ps(&x[last_i+48]));
				 _mm256_store_ps(&y[last_i+56],
				             _mm256_load_ps(&x[last_i+56]));
			    }

			    for(; (last_i+31) < n; last_i += 32) {

                               	 _mm256_store_ps(&y[last_i+0],
				             _mm256_load_ps(&x[last_i+0]));
				 _mm256_store_ps(&y[last_i+8],
				             _mm256_load_ps(&x[last_i+8]));
				 _mm256_store_ps(&y[last_i+16],
				             _mm256_load_ps(&x[last_i+16]));
				 _mm256_store_ps(&y[last_i+24],
				             _mm256_load_ps(&x[last_i+24]));
			    }

			    for(; (last_i+15) < n; last_i += 16) {
                               
                                 _mm256_store_ps(&y[last_i+0],
				             _mm256_load_ps(&x[last_i+0]));
				 _mm256_store_ps(&y[last_i+8],
				             _mm256_load_ps(&x[last_i+8]));
			    }

			    for(; (last_i+7) < n; last_i += 8) {
                           
                                 _mm256_store_ps(&y[last_i+0],
				             _mm256_load_ps(&x[last_i+0]));
			    }

			    for(; (last_i+0) < n; last_i += 1) {
                                y[last_i] = x[last_i];
			    }
			    
			}
			else {

			    for(int32_t i = 0; i != n; ++i) {
                                *y = *x;
				 x += incx;
				 y += incy;
			    }

			}
		  }


		     /* COPY kernel double-precision unaligned */
		   
                   void gms::math::dcopy_u_ymm4r8_unroll16x(const int32_t n,
		                                 double * __restrict x,
						 const int32_t incx,
						 double * __restrict y,
						 const int32_t incy) {
                         if(__builtin_expect(0==n,0)) { return;}

			 __ATTR_ALIGN__(32) __m256d xv[16];
			 int32_t i;

			 if(__builtin_expect(1==incx,1) &&
			    __builtin_expect(1==incy,1)) {

                              for(i = 0; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  _mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
				  xv[0] = _mm256_loadu_pd(&x[i+0]);
				  _mm256_storeu_pd(&y[i+0],xv[0]);
				  xv[1] = _mm256_loadu_pd(&x[i+4]);
				  _mm256_storeu_pd(&y[i+4],xv[1]);
				  xv[2] = _mm256_loadu_pd(&x[i+8]);
				  _mm256_storeu_pd(&y[i+8],xv[2]);
				  xv[3] = _mm256_loadu_pd(&x[i+12]);
				  _mm256_storeu_pd(&y[i+12],xv[3]);
				  xv[4] = _mm256_loadu_pd(&x[i+16]);
				  _mm256_storeu_pd(&y[i+16],xv[4]);
				  _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  xv[5] = _mm256_loadu_pd(&x[i+20]);
				  _mm256_storeu_pd(&y[i+20],xv[5]);
				  xv[6] = _mm256_loadu_pd(&x[i+24]);
				  _mm256_storeu_pd(&y[i+24],xv[6]);
				  xv[7] = _mm256_loadu_pd(&x[i+28]);
				  _mm256_storeu_pd(&y[i+28],xv[7]);
				  xv[8] = _mm256_loadu_pd(&x[i+32]);
				  _mm256_storeu_pd(&y[i+32],xv[8]);
				  _mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
				   xv[9] = _mm256_loadu_pd(&x[i+36]);
				  _mm256_storeu_pd(&y[i+36],xv[9]);
				  xv[10] = _mm256_loadu_pd(&x[i+40]);
				  _mm256_storeu_pd(&y[i+40],xv[10]);
				   xv[11] = _mm256_loadu_pd(&x[i+44]);
				  _mm256_storeu_pd(&y[i+44],xv[11]);
				  xv[12] = _mm256_loadu_pd(&x[i+48]);
				  _mm256_storeu_pd(&y[i+48],xv[12]);
				  _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				  xv[13] = _mm256_loadu_pd(&x[i+52]);
				  _mm256_storeu_pd(&y[i+52],xv[13]);
				  xv[14] = _mm256_loadu_pd(&x[i+56]);
				  _mm256_storeu_pd(&y[i+56],xv[14]);
				  xv[15] = _mm256_loadu_pd(&x[i+60]);
				  _mm256_storeu_pd(&y[i+60],xv[15]);
#else
                                  _mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
                                  xv[0] = _mm256_loadu_pd(&x[i+0]);
				  xv[1] = _mm256_loadu_pd(&x[i+4]);
				  xv[2] = _mm256_loadu_pd(&x[i+8]);
				  xv[3] = _mm256_loadu_pd(&x[i+12]);
				  xv[4] = _mm256_loadu_pd(&x[i+16]);
				  _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  xv[5] = _mm256_loadu_pd(&x[i+20]);
				  xv[6] = _mm256_loadu_pd(&x[i+24]);
				  xv[7] = _mm256_loadu_pd(&x[i+28]);
				  xv[8] = _mm256_loadu_pd(&x[i+32]);
				  _mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
				  xv[9] = _mm256_loadu_pd(&x[i+36]);
				  xv[10] = _mm256_loadu_pd(&x[i+40]);
				  xv[11] = _mm256_loadu_pd(&x[i+44]);
				  xv[12] = _mm256_loadu_pd(&x[i+48]);
				  _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				  xv[13] = _mm256_loadu_pd(&x[i+52]);
				  xv[14] = _mm256_loadu_pd(&x[i+56]);
				  xv[15] = _mm256_loadu_pd(&x[i+60]);

				  _mm256_storeu_pd(&y[i+0],xv[0]);
				  _mm256_storeu_pd(&y[i+4],xv[1]);
				  _mm256_storeu_pd(&y[i+8],xv[2]);
				  _mm256_storeu_pd(&y[i+12],xv[3]);
				  _mm256_storeu_pd(&y[i+16],xv[4]);
				  _mm256_storeu_pd(&y[i+20],xv[5]);
				  _mm256_storeu_pd(&y[i+24],xv[6]);
				  _mm256_storeu_pd(&y[i+28],xv[7]);
				  _mm256_storeu_pd(&y[i+32],xv[8]);
				  _mm256_storeu_pd(&y[i+36],xv[9]);
				  _mm256_storeu_pd(&y[i+40],xv[10]);
				  _mm256_storeu_pd(&y[i+44],xv[11]);
				  _mm256_storeu_pd(&y[i+48],xv[12]);
				  _mm256_storeu_pd(&y[i+52],xv[13]);
				  _mm256_storeu_pd(&y[i+56],xv[14]);
				  _mm256_storeu_pd(&y[i+60],xv[15]);
				  
				  
#endif

			      }

			      for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                   xv[0] = _mm256_loadu_pd(&x[i+0]);
				  _mm256_storeu_pd(&y[i+0],xv[0]);
				  xv[1] = _mm256_loadu_pd(&x[i+4]);
				  _mm256_storeu_pd(&y[i+4],xv[1]);
				  xv[2] = _mm256_loadu_pd(&x[i+8]);
				  _mm256_storeu_pd(&y[i+8],xv[2]);
				  xv[3] = _mm256_loadu_pd(&x[i+12]);
				  _mm256_storeu_pd(&y[i+12],xv[3]);
				  xv[4] = _mm256_loadu_pd(&x[i+16]);
				  _mm256_storeu_pd(&y[i+16],xv[4]);
				  //_mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  xv[5] = _mm256_loadu_pd(&x[i+20]);
				  _mm256_storeu_pd(&y[i+20],xv[5]);
				  xv[6] = _mm256_loadu_pd(&x[i+24]);
				  _mm256_storeu_pd(&y[i+24],xv[6]);
				  xv[7] = _mm256_loadu_pd(&x[i+28]);
				  _mm256_storeu_pd(&y[i+28],xv[7]);  
#else
                                  //_mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
                                  xv[0] = _mm256_loadu_pd(&x[i+0]);
				  xv[1] = _mm256_loadu_pd(&x[i+4]);
				  xv[2] = _mm256_loadu_pd(&x[i+8]);
				  xv[3] = _mm256_loadu_pd(&x[i+12]);
				  xv[4] = _mm256_loadu_pd(&x[i+16]);
				  //_mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  xv[5] = _mm256_loadu_pd(&x[i+20]);
				  xv[6] = _mm256_loadu_pd(&x[i+24]);
				  xv[7] = _mm256_loadu_pd(&x[i+28]);

				  _mm256_storeu_pd(&y[i+0],xv[0]);
				  _mm256_storeu_pd(&y[i+4],xv[1]);
				  _mm256_storeu_pd(&y[i+8],xv[2]);
				  _mm256_storeu_pd(&y[i+12],xv[3]);
				  _mm256_storeu_pd(&y[i+16],xv[4]);
				  _mm256_storeu_pd(&y[i+20],xv[5]);
				  _mm256_storeu_pd(&y[i+24],xv[6]);
				  _mm256_storeu_pd(&y[i+28],xv[7]);
#endif
			      }

			      for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  xv[0] = _mm256_loadu_pd(&x[i+0]);
				  _mm256_storeu_pd(&y[i+0],xv[0]);
				  xv[1] = _mm256_loadu_pd(&x[i+4]);
				  _mm256_storeu_pd(&y[i+4],xv[1]);
				  xv[2] = _mm256_loadu_pd(&x[i+8]);
				  _mm256_storeu_pd(&y[i+8],xv[2]);
				  xv[3] = _mm256_loadu_pd(&x[i+12]);
				  _mm256_storeu_pd(&y[i+12],xv[3]);
#else
                                  xv[0] = _mm256_loadu_pd(&x[i+0]);
				  xv[1] = _mm256_loadu_pd(&x[i+4]);
				  xv[2] = _mm256_loadu_pd(&x[i+8]);
				  xv[3] = _mm256_loadu_pd(&x[i+12]);

				  _mm256_storeu_pd(&y[i+0],xv[0]);
				  _mm256_storeu_pd(&y[i+4],xv[1]);
				  _mm256_storeu_pd(&y[i+8],xv[2]);
				  _mm256_storeu_pd(&y[i+12],xv[3]);
#endif
			      }

			      for(; (i+7) < n; i += 8) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  xv[0] = _mm256_loadu_pd(&x[i+0]);
				  _mm256_storeu_pd(&y[i+0],xv[0]);
				  xv[1] = _mm256_loadu_pd(&x[i+4]);
				  _mm256_storeu_pd(&y[i+4],xv[1]);
#else
                                  xv[0] = _mm256_loadu_pd(&x[i+0]);
				  xv[1] = _mm256_loadu_pd(&x[i+4]);

				  _mm256_storeu_pd(&y[i+0],xv[0]);
				  _mm256_storeu_pd(&y[i+4],xv[1]);
#endif
			      }

			      for(; (i+3) < n; i += 4) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                   xv[0] = _mm256_loadu_pd(&x[i+0]);
				  _mm256_storeu_pd(&y[i+0],xv[0]);
#else
                                   xv[0] = _mm256_loadu_pd(&x[i+0]);
				   _mm256_storeu_pd(&y[i+0],xv[0]);
#endif
			      }

			      for(; (i+0) < n; i += 1) {
                                    y[i] = x[i];
			      }
			 }
			 else {
                                for(i = 0; i != n; ++i) {
                                    *y = *x;
				    x += incx;
				    y += incy;
				}
			 }

		   }


		      /* COPY kernel double-precision aligned */
		 
                   void gms::math::dcopy_a_ymm4r8_unroll16x(const int32_t n,
		                                 double * __restrict __ATTR_ALIGN__(32) x,
						 const int32_t incx,
						 double * __restrict __ATTR_ALIGN__(32) y,
						 const int32_t incy) {
                         if(__builtin_expect(0==n,0)) { return;}

			 __ATTR_ALIGN__(32) __m256d xv[16];
			 int32_t i;

			 if(__builtin_expect(1==incx,1) &&
			    __builtin_expect(1==incy,1)) {
#if defined(__INTEL_COMPILER) || defined(__ICC)
                             __assume_aligned(x,32);
			     __assume_aligned(y,32);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                             x = (double*)__builtin_assume_aligned(x,32);
			     y = (double*)__builtin_assume_aligned(y,32);
#endif
                              for(i = 0; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  _mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
				  xv[0] = _mm256_load_pd(&x[i+0]);
				  _mm256_store_pd(&y[i+0],xv[0]);
				  xv[1] = _mm256_load_pd(&x[i+4]);
				  _mm256_store_pd(&y[i+4],xv[1]);
				  xv[2] = _mm256_load_pd(&x[i+8]);
				  _mm256_store_pd(&y[i+8],xv[2]);
				  xv[3] = _mm256_load_pd(&x[i+12]);
				  _mm256_store_pd(&y[i+12],xv[3]);
				  xv[4] = _mm256_load_pd(&x[i+16]);
				  _mm256_store_pd(&y[i+16],xv[4]);
				  _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  xv[5] = _mm256_load_pd(&x[i+20]);
				  _mm256_store_pd(&y[i+20],xv[5]);
				  xv[6] = _mm256_load_pd(&x[i+24]);
				  _mm256_store_pd(&y[i+24],xv[6]);
				  xv[7] = _mm256_load_pd(&x[i+28]);
				  _mm256_store_pd(&y[i+28],xv[7]);
				  xv[8] = _mm256_load_pd(&x[i+32]);
				  _mm256_store_pd(&y[i+32],xv[8]);
				  _mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
				   xv[9] = _mm256_load_pd(&x[i+36]);
				  _mm256_store_pd(&y[i+36],xv[9]);
				  xv[10] = _mm256_load_pd(&x[i+40]);
				  _mm256_store_pd(&y[i+40],xv[10]);
				   xv[11] = _mm256_load_pd(&x[i+44]);
				  _mm256_store_pd(&y[i+44],xv[11]);
				  xv[12] = _mm256_load_pd(&x[i+48]);
				  _mm256_store_pd(&y[i+48],xv[12]);
				  _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				  xv[13] = _mm256_load_pd(&x[i+52]);
				  _mm256_store_pd(&y[i+52],xv[13]);
				  xv[14] = _mm256_load_pd(&x[i+56]);
				  _mm256_store_pd(&y[i+56],xv[14]);
				  xv[15] = _mm256_load_pd(&x[i+60]);
				  _mm256_store_pd(&y[i+60],xv[15]);
#else
                                  _mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
                                  xv[0] = _mm256_load_pd(&x[i+0]);
				  xv[1] = _mm256_load_pd(&x[i+4]);
				  xv[2] = _mm256_load_pd(&x[i+8]);
				  xv[3] = _mm256_load_pd(&x[i+12]);
				  xv[4] = _mm256_load_pd(&x[i+16]);
				  _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  xv[5] = _mm256_load_pd(&x[i+20]);
				  xv[6] = _mm256_load_pd(&x[i+24]);
				  xv[7] = _mm256_load_pd(&x[i+28]);
				  xv[8] = _mm256_load_pd(&x[i+32]);
				  _mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
				  xv[9] = _mm256_load_pd(&x[i+36]);
				  xv[10] = _mm256_load_pd(&x[i+40]);
				  xv[11] = _mm256_load_pd(&x[i+44]);
				  xv[12] = _mm256_load_pd(&x[i+48]);
				  _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				  xv[13] = _mm256_load_pd(&x[i+52]);
				  xv[14] = _mm256_load_pd(&x[i+56]);
				  xv[15] = _mm256_load_pd(&x[i+60]);

				  _mm256_store_pd(&y[i+0],xv[0]);
				  _mm256_store_pd(&y[i+4],xv[1]);
				  _mm256_store_pd(&y[i+8],xv[2]);
				  _mm256_store_pd(&y[i+12],xv[3]);
				  _mm256_store_pd(&y[i+16],xv[4]);
				  _mm256_store_pd(&y[i+20],xv[5]);
				  _mm256_store_pd(&y[i+24],xv[6]);
				  _mm256_store_pd(&y[i+28],xv[7]);
				  _mm256_store_pd(&y[i+32],xv[8]);
				  _mm256_store_pd(&y[i+36],xv[9]);
				  _mm256_store_pd(&y[i+40],xv[10]);
				  _mm256_store_pd(&y[i+44],xv[11]);
				  _mm256_store_pd(&y[i+48],xv[12]);
				  _mm256_store_pd(&y[i+52],xv[13]);
				  _mm256_store_pd(&y[i+56],xv[14]);
				  _mm256_store_pd(&y[i+60],xv[15]);
				  
				  
#endif

			      }

			      for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                   xv[0] = _mm256_load_pd(&x[i+0]);
				  _mm256_store_pd(&y[i+0],xv[0]);
				  xv[1] = _mm256_load_pd(&x[i+4]);
				  _mm256_store_pd(&y[i+4],xv[1]);
				  xv[2] = _mm256_load_pd(&x[i+8]);
				  _mm256_store_pd(&y[i+8],xv[2]);
				  xv[3] = _mm256_load_pd(&x[i+12]);
				  _mm256_store_pd(&y[i+12],xv[3]);
				  xv[4] = _mm256_load_pd(&x[i+16]);
				  _mm256_store_pd(&y[i+16],xv[4]);
				  //_mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  xv[5] = _mm256_load_pd(&x[i+20]);
				  _mm256_store_pd(&y[i+20],xv[5]);
				  xv[6] = _mm256_load_pd(&x[i+24]);
				  _mm256_store_pd(&y[i+24],xv[6]);
				  xv[7] = _mm256_load_pd(&x[i+28]);
				  _mm256_store_pd(&y[i+28],xv[7]);  
#else
                                  //_mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
                                  xv[0] = _mm256_load_pd(&x[i+0]);
				  xv[1] = _mm256_load_pd(&x[i+4]);
				  xv[2] = _mm256_load_pd(&x[i+8]);
				  xv[3] = _mm256_load_pd(&x[i+12]);
				  xv[4] = _mm256_load_pd(&x[i+16]);
				  //_mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  xv[5] = _mm256_load_pd(&x[i+20]);
				  xv[6] = _mm256_load_pd(&x[i+24]);
				  xv[7] = _mm256_load_pd(&x[i+28]);

				  _mm256_store_pd(&y[i+0],xv[0]);
				  _mm256_store_pd(&y[i+4],xv[1]);
				  _mm256_store_pd(&y[i+8],xv[2]);
				  _mm256_store_pd(&y[i+12],xv[3]);
				  _mm256_store_pd(&y[i+16],xv[4]);
				  _mm256_store_pd(&y[i+20],xv[5]);
				  _mm256_store_pd(&y[i+24],xv[6]);
				  _mm256_store_pd(&y[i+28],xv[7]);
#endif
			      }

			      for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  xv[0] = _mm256_load_pd(&x[i+0]);
				  _mm256_store_pd(&y[i+0],xv[0]);
				  xv[1] = _mm256_load_pd(&x[i+4]);
				  _mm256_store_pd(&y[i+4],xv[1]);
				  xv[2] = _mm256_load_pd(&x[i+8]);
				  _mm256_store_pd(&y[i+8],xv[2]);
				  xv[3] = _mm256_load_pd(&x[i+12]);
				  _mm256_store_pd(&y[i+12],xv[3]);
#else
                                  xv[0] = _mm256_load_pd(&x[i+0]);
				  xv[1] = _mm256_load_pd(&x[i+4]);
				  xv[2] = _mm256_load_pd(&x[i+8]);
				  xv[3] = _mm256_load_pd(&x[i+12]);

				  _mm256_store_pd(&y[i+0],xv[0]);
				  _mm256_store_pd(&y[i+4],xv[1]);
				  _mm256_store_pd(&y[i+8],xv[2]);
				  _mm256_store_pd(&y[i+12],xv[3]);
#endif
			      }

			      for(; (i+7) < n; i += 8) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  xv[0] = _mm256_load_pd(&x[i+0]);
				  _mm256_store_pd(&y[i+0],xv[0]);
				  xv[1] = _mm256_load_pd(&x[i+4]);
				  _mm256_store_pd(&y[i+4],xv[1]);
#else
                                  xv[0] = _mm256_load_pd(&x[i+0]);
				  xv[1] = _mm256_load_pd(&x[i+4]);

				  _mm256_store_pd(&y[i+0],xv[0]);
				  _mm256_store_pd(&y[i+4],xv[1]);
#endif
			      }

			      for(; (i+3) < n; i += 4) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                   xv[0] = _mm256_load_pd(&x[i+0]);
				  _mm256_store_pd(&y[i+0],xv[0]);
#else
                                   xv[0] = _mm256_load_pd(&x[i+0]);
				   _mm256_store_pd(&y[i+0],xv[0]);
#endif
			      }

			      for(; (i+0) < n; i += 1) {
                                    y[i] = x[i];
			      }
			 }
			 else {
                                for(i = 0; i != n; ++i) {
                                    *y = *x;
				    x += incx;
				    y += incy;
				}
			 }

		   }



	                /* COPY kernel double-precision aligned */
		  
                   void gms::math::dcopy_a_ymm4r8_unroll16x_omp(const int32_t n,
		                                     double * __restrict __ATTR_ALIGN__(32) x,
						     const int32_t incx,
						     double * __restrict __ATTR_ALIGN__(32) y,
						     const int32_t incy) {
                         if(__builtin_expect(0==n,0)) { return;}

			 __ATTR_ALIGN__(32) __m256d xv[16];
                       
			 int32_t i;
			 int32_t last_i;
			 last_i = 0;

			 if(__builtin_expect(1==incx,1) &&
			    __builtin_expect(1==incy,1)) {
#if defined(__INTEL_COMPILER) || defined(__ICC)
                             __assume_aligned(x,32);
			     __assume_aligned(y,32);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                             x = (double*)__builtin_assume_aligned(x,32);
			     y = (double*)__builtin_assume_aligned(y,32);
#endif
#pragma omp parallel for schedule(static,128) default(none) \
                             lastprivate(last_i) private(i) shared(n,x,y)
                              for(i = 0; (i+63) < n; i += 64) {
                                  last_i = i;
                                  _mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
				  _mm256_store_pd(&y[i+0],
						  _mm256_load_pd(&x[i+0]));
				  _mm256_store_pd(&y[i+4],
						  _mm256_load_pd(&x[i+4]));
				  _mm256_store_pd(&y[i+8],
						  _mm256_load_pd(&x[i+8]));
				  _mm256_store_pd(&y[i+12],
						  _mm256_load_pd(&x[i+12]));
				  _mm256_store_pd(&y[i+16],
						  _mm256_load_pd(&x[i+16]));		  
				  _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  _mm256_store_pd(&y[i+20],
						  _mm256_load_pd(&x[i+20]));
				  _mm256_store_pd(&y[i+24],
						  _mm256_load_pd(&x[i+24]));
				  _mm256_store_pd(&y[i+28],
						  _mm256_load_pd(&x[i+28]));
				  _mm256_store_pd(&y[i+32],
						  _mm256_load_pd(&x[i+32]));		  
				  _mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
				  _mm256_store_pd(&y[i+36],
						  _mm256_load_pd(&x[i+36]));
				  _mm256_store_pd(&y[i+40],
						  _mm256_load_pd(&x[i+40]));
				  _mm256_store_pd(&y[i+44],
						  _mm256_load_pd(&x[i+44]));
				  _mm256_store_pd(&y[i+48],
						  _mm256_load_pd(&x[i+48]));		  
				  _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				  _mm256_store_pd(&y[i+52],
						  _mm256_load_pd(&x[i+52]));
				  _mm256_store_pd(&y[i+56],
						  _mm256_load_pd(&x[i+56]));
				  _mm256_store_pd(&y[i+60],
					          _mm256_load_pd(&x[i+60]));		  

			      }

			      for(; (last_i+31) < n; last_i += 32) {
                                  _mm256_store_pd(&y[last_i+0],
						  _mm256_load_pd(&x[last_i+0]));
				  _mm256_store_pd(&y[last_i+4],
						  _mm256_load_pd(&x[last_i+4]));
				  _mm256_store_pd(&y[last_i+8],
						  _mm256_load_pd(&x[last_i+8]));
				  _mm256_store_pd(&y[last_i+12],
						  _mm256_load_pd(&x[last_i+12]));
				  _mm256_store_pd(&y[last_i+16],
						  _mm256_load_pd(&x[last_i+16]));		  
				   _mm256_store_pd(&y[last_i+20],
						  _mm256_load_pd(&x[last_i+20]));
				  _mm256_store_pd(&y[last_i+24],
						  _mm256_load_pd(&x[last_i+24]));
				  _mm256_store_pd(&y[last_i+28],
						  _mm256_load_pd(&x[last_i+28]));
			      }

			      for(; (last_i+15) < n; last_i += 16) {
			          _mm256_store_pd(&y[last_i+0],
						  _mm256_load_pd(&x[last_i+0]));
				  _mm256_store_pd(&y[last_i+4],
						  _mm256_load_pd(&x[last_i+4]));
				  _mm256_store_pd(&y[last_i+8],
						  _mm256_load_pd(&x[last_i+8]));
				  _mm256_store_pd(&y[last_i+12],
						  _mm256_load_pd(&x[last_i+12]));

			      }

			      for(; (last_i+7) < n; last_i += 8) {
			          _mm256_store_pd(&y[last_i+0],
						  _mm256_load_pd(&x[last_i+0]));
				  _mm256_store_pd(&y[last_i+4],
						  _mm256_load_pd(&x[last_i+4]));

			      }

			      for(; (last_i+3) < n; last_i += 4) {
			             _mm256_store_pd(&y[last_i+0],
						  _mm256_load_pd(&x[last_i+0]));

			      }

			      for(; (last_i+0) < n; last_i += 1) {
                                    y[last_i] = x[last_i];
			      }
			 }
			 else {
                                for(i = 0; i != n; ++i) {
                                    *y = *x;
				    x += incx;
				    y += incy;
				}
			 }

		   }





		  

 
