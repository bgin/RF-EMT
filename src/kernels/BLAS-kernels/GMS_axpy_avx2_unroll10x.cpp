


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


#include <cstdint>
#include <immintrin.h>
#include <omp.h>
#include "GMS_axpy_avx2_unroll10x.h"
#include "GMS_config.h"


              /*
                   SAXPY kernel unaligned.
               */
	     
	      void gms::math::saxpy_u_ymm8r4_unroll10x(const int32_t n,
	                                    const float alpha,
	                                    float * __restrict  x,
					    const int32_t incx,
					    float * __restrict  y,
					    const int32_t incy) {
                   if(__builtin_expect(0==n,0) ||
		      __builtin_expect(alpha==0.0f,0)) { return;}
		   
		   __ATTR_ALIGN__(32) __m256 xv[10];
		   __ATTR_ALIGN__(32) __m256 yv[10];
		   __ATTR_ALIGN__(32) __m256 zv[10];
                   __m256 valpha;
		   float * __restrict x0 = NULL;
		   float * __restrict y0 = NULL;
		   int32_t i;
		   x0 = x;
		   y0 = y;

		   if(__builtin_expect(incx==1,1) &&
		      __builtin_expect(incy==1,1) {

                      valpha = _mm256_broadcast_ss(alpha);
		      // Unrolled 10 times in order exploit the ratio
		      // between the FMA latency to its throughput.
		      for(i = 0; (i+79) < n; i += 80) {
		          _mm_prefetch((const char*)&x0[i+80],_MM_HINT_T0);
                          xv[0] = _mm256_loadu_ps(&x0[i+0]);
			  yv[0] = _mm256_loadu_ps(&y0[i+0]);
			  zv[0] = _mm256_fmadd_ps(xv[0],valpha,yv[0]);
			  _mm256_storeu_ps(&y0[i+0],zv[0]);
			  xv[1] = _mm256_loadu_ps(&x0[i+8]);
			  yv[1] = _mm256_loadu_ps(&y0[i+8]);
			  zv[1] = _mm256_fmadd_ps(xv[1],valpha,yv[1]);
			  _mm256_storeu_ps(&y0[i+8],zv[1]);
			  xv[2] = _mm256_loadu_ps(&x0[i+16]);
			  yv[2] = _mm256_loadu_ps(&y0[i+16]);
			  zv[2] = _mm256_fmadd_ps(xv[2],valpha,yv[2]);
			  _mm256_storeu_ps(&y0[i+16],zv[2]);
			  xv[3] = _mm256_loadu_ps(&x0[i+24]);
			  yv[3] = _mm256_loadu_ps(&y0[i+24]);
			  zv[3] = _mm256_fmadd_ps(xv[3],valpha,yv[3]);
			  _mm256_storeu_ps(&y0[i+24],zv[3]);
			   xv[4] = _mm256_loadu_ps(&x0[i+32]);
			  yv[4] = _mm256_loadu_ps(&y0[i+32]);
			  zv[4] = _mm256_fmadd_ps(xv[4],valpha,yv[4]);
			  _mm256_storeu_ps(&y0[i+32],zv[4]);
			  xv[5] = _mm256_loadu_ps(&x0[i+40]);
			  yv[5] = _mm256_loadu_ps(&y0[i+40]);
			  zv[5] = _mm256_fmadd_ps(xv[5],valpha,yv[5]);
			  _mm256_storeu_ps(&y0[i+40],zv[5]);
			  xv[6] = _mm256_loadu_ps(&x0[i+48]);
			  yv[6] = _mm256_loadu_ps(&y0[i+48]);
			  zv[6] = _mm256_fmadd_ps(xv[6],valpha,yv[6]);
			  _mm256_storeu_ps(&y[i+48],zv[6]);
			  xv[7] = _mm256_loadu_ps(&x0[i+56]);
			  yv[7] = _mm256_loadu_ps(&y0[i+56]);
			  zv[7] = _mm256_fmadd_ps(xv[7],valpha,yv[7]);
			  _mm256_storeu_ps(&y0[i+56],zv[7]);
			  xv[8] = _mm256_loadu_ps(&x0[i+64]);
			  yv[8] = _mm256_loadu_ps(&y0[i+64]);
			  zv[8] = _mm256_fmadd_ps(xv[8],valpha,yv[8]);
			  _mm256_storeu_ps(&y0[i+64],zv[8]);
			  xv[9] = _mm256_loadu_ps(&x0[i+72]);
			  yv[9] = _mm256_loadu_ps(&y0[i+72]);
			  zv[9] = _mm256_fmadd_ps(xv[9],valpha,yv[9]);
			  _mm256_storeu_ps(&y0[i+72],zv[9]);
		      }
		      for(; (i+39) < n; i += 40) {

		          xv[0] = _mm256_loadu_ps(&x0[i+0]);
			  yv[0] = _mm256_loadu_ps(&y0[i+0]);
			  zv[0] = _mm256_fmadd_ps(xv[0],valpha,yv[0]);
			  _mm256_storeu_ps(&y0[i+0],zv[0]);
			  xv[1] = _mm256_loadu_ps(&x0[i+8]);
			  yv[1] = _mm256_loadu_ps(&y0[i+8]);
			  zv[1] = _mm256_fmadd_ps(xv[1],valpha,yv[1]);
			  _mm256_storeu_ps(&y0[i+8],zv[1]);
			  xv[2] = _mm256_loadu_ps(&x0[i+16]);
			  yv[2] = _mm256_loadu_ps(&y0[i+16]);
			  zv[2] = _mm256_fmadd_ps(xv[2],valpha,yv[2]);
			  _mm256_storeu_ps(&y0[i+16],zv[2]);
			  xv[3] = _mm256_loadu_ps(&x0[i+24]);
			  yv[3] = _mm256_loadu_ps(&y0[i+24]);
			  zv[3] = _mm256_fmadd_ps(xv[3],valpha,yv[3]);
			  _mm256_storeu_ps(&y0[i+24],zv[3]);
			  xv[4] = _mm256_loadu_ps(&x0[i+32]);
			  yv[4] = _mm256_loadu_ps(&y0[i+32]);
			  zv[4] = _mm256_fmadd_ps(xv[4],valpha,yv[4]);
			  _mm256_storeu_ps(&y0[i+32],zv[4]);
		      }
		      for(; (i+31) < n; i += 32) {

		          xv[0] = _mm256_loadu_ps(&x0[i+0]);
			  yv[0] = _mm256_loadu_ps(&y0[i+0]);
			  zv[0] = _mm256_fmadd_ps(xv[0],valpha,yv[0]);
			  _mm256_storeu_ps(&y0[i+0],zv[0]);
			  xv[1] = _mm256_loadu_ps(&x0[i+8]);
			  yv[1] = _mm256_loadu_ps(&y0[i+8]);
			  zv[1] = _mm256_fmadd_ps(xv[1],valpha,yv[1]);
			  _mm256_storeu_ps(&y0[i+8],zv[1]);
			  xv[2] = _mm256_loadu_ps(&x0[i+16]);
			  yv[2] = _mm256_loadu_ps(&y0[i+16]);
			  zv[2] = _mm256_fmadd_ps(xv[2],valpha,yv[2]);
			  _mm256_storeu_ps(&y0[i+16],zv[2]);
			  xv[3] = _mm256_loadu_ps(&x0[i+24]);
			  yv[3] = _mm256_loadu_ps(&y0[i+24]);
			  zv[3] = _mm256_fmadd_ps(xv[3],valpha,yv[3]);
			  _mm256_storeu_ps(&y0[i+24],zv[3]);
		      }
		      for(; (i+15) < n; i += 16) {

		          xv[0] = _mm256_loadu_ps(&x0[i+0]);
			  yv[0] = _mm256_loadu_ps(&y0[i+0]);
			  zv[0] = _mm256_fmadd_ps(xv[0],valpha,yv[0]);
			  _mm256_storeu_ps(&y0[i+0],zv[0]);
			  xv[1] = _mm256_loadu_ps(&x0[i+8]);
			  yv[1] = _mm256_loadu_ps(&y0[i+8]);
			  zv[1] = _mm256_fmadd_ps(xv[1],valpha,yv[1]);
			  _mm256_storeu_ps(&y0[i+8],zv[1]);
		      }
		      for(; (i+7) < n; i += 8) {

		          xv[0] = _mm256_loadu_ps(&x0[i+0]);
			  yv[0] = _mm256_loadu_ps(&y0[i+0]);
			  zv[0] = _mm256_fmadd_ps(xv[0],valpha,yv[0]);
			  _mm256_storeu_ps(&y0[i+0],zv[0]);
		      }
		      _mm256_vzeroupper();
		      for(; (i+0) < n; i += 1) {
                          y0[i] += alpha * x0[i];
		      }
		   }
		   else {

                       for(i = 0; i != n; ++i) {
		           const float tx0 = *x0;
                           *y0 += alpha * tx0;
			   x0 += incx;
			   y0 += incy;
		       }
		   }
	     }

            /*
                 SAXPY kernel aligned
              */
           
	     void gms::math::saxpy_a_ymm8r4_unroll10x(const int32_t n,
	                                    const float alpha,
	                                    float * __restrict __ATTR_ALIGN__(32) x,
					    const int32_t incx,
					    float * __restrict  __ATTR_ALIGN__(32) y,
					    const int32_t incy) {
                   if(__builtin_expect(0==n,0) || 
		      __builtin_expect(alpha==0.0f,0)) { return;}
		   
		   __ATTR_ALIGN__(32) __m256 xv[10];
		   __ATTR_ALIGN__(32) __m256 yv[10];
		   __ATTR_ALIGN__(32) __m256 zv[10];
                   __m256 valpha;
		   float * __restrict x0 = NULL;
		   float * __restrict y0 = NULL;
		   int32_t i;
		   x0 = x;
		   y0 = y;

		   if(__builtin_expect(incx==1,1) &&
		      __builtin_expect(incy==1,1)) {
#if defined(INTEL_COMPILER) || defined(__ICC)
                      __assume_aligned(x,32);
		      __assume_aligned(y,32);
		      __assume_aligned(x0,32);
		      __assume_aligned(y0,32);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC)
                      x  = (float*)__builtin_assume_aligned(x,32);
		      y  = (float*)__builtin_assume_aligned(y,32);
		      x0 = (float*)__builtin_assume_aligned(x0,32);
		      y0 = (float*)__builtin_assume_aligned(y0,32);
#endif
                      valpha = _mm256_broadcast_ss(alpha);
		      // Unrolled 10 times in order exploit the ratio
		      // between the FMA latency to its throughput.
		      for(i = 0; (i+79) < n; i += 80) {
		          _mm_prefetch((const char*)&x[i+80],_MM_HINT_T0);
                          xv[0] = _mm256_load_ps(&x0[i+0]);
			  yv[0] = _mm256_load_ps(&y0[i+0]);
			  zv[0] = _mm256_fmadd_ps(xv[0],valpha,yv[0]);
			  _mm256_store_ps(&y0[i+0],zv[0]);
			  xv[1] = _mm256_load_ps(&x0[i+8]);
			  yv[1] = _mm256_load_ps(&y0[i+8]);
			  zv[1] = _mm256_fmadd_ps(xv[1],valpha,yv[1]);
			  _mm256_store_ps(&y0[i+8],zv[1]);
			  xv[2] = _mm256_load_ps(&x0[i+16]);
			  yv[2] = _mm256_load_ps(&y0[i+16]);
			  zv[2] = _mm256_fmadd_ps(xv[2],valpha,yv[2]);
			  _mm256_store_ps(&y0[i+16],zv[2]);
			  xv[3] = _mm256_load_ps(&x0[i+24]);
			  yv[3] = _mm256_load_ps(&y0[i+24]);
			  zv[3] = _mm256_fmadd_ps(xv[3],valpha,yv[3]);
			  _mm256_store_ps(&y0[i+24],zv[3]);
			  xv[4] = _mm256_load_ps(&x0[i+32]);
			  yv[4] = _mm256_load_ps(&y0[i+32]);
			  zv[4] = _mm256_fmadd_ps(xv[4],valpha,yv[4]);
			  _mm256_store_ps(&y0[i+32],zv[4]);
			  xv[5] = _mm256_load_ps(&x0[i+40]);
			  yv[5] = _mm256_load_ps(&y0[i+40]);
			  zv[5] = _mm256_fmadd_ps(xv[5],valpha,yv[5]);
			  _mm256_store_ps(&y0[i+40],zv[5]);
			  xv[6] = _mm256_load_ps(&x0[i+48]);
			  yv[6] = _mm256_load_ps(&y0[i+48]);
			  zv[6] = _mm256_fmadd_ps(xv[6],valpha,yv[6]);
			  _mm256_store_ps(&y0[i+48],zv[6]);
			  xv[7] = _mm256_load_ps(&x0[i+56]);
			  yv[7] = _mm256_load_ps(&y0[i+56]);
			  zv[7] = _mm256_fmadd_ps(xv[7],valpha,yv[7]);
			  _mm256_store_ps(&y0[i+56],zv[7]);
			  xv[8] = _mm256_load_ps(&x0[i+64]);
			  yv[8] = _mm256_load_ps(&y0[i+64]);
			  zv[8] = _mm256_fmadd_ps(xv[8],valpha,yv[8]);
			  _mm256_store_ps(&y0[i+64],zv[8]);
			  xv[9] = _mm256_load_ps(&x0[i+72]);
			  yv[9] = _mm256_load_ps(&y0[i+72]);
			  zv[9] = _mm256_fmadd_ps(xv[9],valpha,yv[9]);
			  _mm256_store_ps(&y0[i+72],zv[9]);
		      }
		      for(; (i+39) < n; i += 40) {

		          xv[0] = _mm256_load_ps(&x0[i+0]);
			  yv[0] = _mm256_load_ps(&y0[i+0]);
			  zv[0] = _mm256_fmadd_ps(xv[0],valpha,yv[0]);
			  _mm256_store_ps(&y0[i+0],zv[0]);
			  xv[1] = _mm256_load_ps(&x0[i+8]);
			  yv[1] = _mm256_load_ps(&y0[i+8]);
			  zv[1] = _mm256_fmadd_ps(xv[1],valpha,yv[1]);
			  _mm256_store_ps(&y0[i+8],zv[1]);
			  xv[2] = _mm256_load_ps(&x0[i+16]);
			  yv[2] = _mm256_load_ps(&y0[i+16]);
			  zv[2] = _mm256_fmad_ps(xv[2],valpha,yv[2]);
			  _mm256_store_ps(&y0[i+16],zv[2]);
			  xv[3] = _mm256_load_ps(&x0[i+24]);
			  yv[3] = _mm256_load_ps(&y0[i+24]);
			  zv[3] = _mm256_fmadd_ps(xv[3],valpha,yv[3]);
			  _mm256_store_ps(&y0[i+24],zv[3]);
			  xv[4] = _mm256_load_ps(&x0[i+32]);
			  yv[4] = _mm256_load_ps(&y0[i+32]);
			  zv[4] = _mm256_fmadd_ps(xv[4],valpha,yv[4]);
			  _mm256_store_ps(&y0[i+32],zv[4]);
		      }
		      for(; (i+31) < n; i += 32) {

		          xv[0] = _mm256_load_ps(&x0[i+0]);
			  yv[0] = _mm256_load_ps(&y0[i+0]);
			  zv[0] = _mm256_fmadd_ps(xv[0],valpha,yv[0]);
			  _mm256_store_ps(&y0[i+0],zv[0]);
			  xv[1] = _mm256_load_ps(&x0[i+8]);
			  yv[1] = _mm256_load_ps(&y0[i+8]);
			  zv[1] = _mm256_fmadd_ps(xv[1],valpha,yv[1]);
			  _mm256_store_ps(&y0[i+8],zv[1]);
			  xv[2] = _mm256_load_ps(&x0[i+16]);
			  yv[2] = _mm256_load_ps(&y0[i+16]);
			  zv[2] = _mm256_fmadd_ps(xv[2],valpha,yv[2]);
			  _mm256_store_ps(&y0[i+16],zv[2]);
			  xv[3] = _mm256_load_ps(&x0[i+24]);
			  yv[3] = _mm256_load_ps(&y0[i+24]);
			  zv[3] = _mm256_fmadd_ps(xv[3],valpha,yv[3]);
			  _mm256_store_ps(&y0[i+24],zv[3]);
		      }
		      for(; (i+15) < n; i += 16) {

		          xv[0] = _mm256_load_ps(&x0[i+0]);
			  yv[0] = _mm256_load_ps(&y0[i+0]);
			  zv[0] = _mm256_fmadd_ps(xv[0],valpha,yv[0]);
			  _mm256_store_ps(&y0[i+0],zv[0]);
			  xv[1] = _mm256_load_ps(&x0[i+8]);
			  yv[1] = _mm256_load_ps(&y0[i+8]);
			  zv[1] = _mm256_fmadd_ps(xv[1],valpha,yv[1]);
			  _mm256_store_ps(&y0[i+8],zv[1]);
		      }
		      for(; (i+7) < n; i += 8) {

		          xv[0] = _mm256_load_ps(&x0[i+0]);
			  yv[0] = _mm256_load_ps(&y0[i+0]);
			  zv[0] = _mm256_fmadd_ps(xv[0],valpha,yv[0]);
			  _mm256_store_ps(&y0[i+0],zv[0]);
		      }
		      _mm256_vzeroupper();
		      for(; (i+0) < n; i += 1) {
                          y0[i] += alpha * x0[i];
		      }
		   }
		   else {

                       for(i = 0; i != n; ++i) {
		           const float tx0 = *x0;
                           *y0 += alpha * tx0;
			   x0 += incx;
			   y0 += incy;
		       }
		   }
	     }

	     /*
                 SAXPY kernel OpenMP parallelized (aligned version only)
              */
           
	     void gms::math::saxpy_a_ymm8r4_unroll10x_omp(const int32_t n,
	                                       const float alpha,
	                                       float * __restrict __ATTR_ALIGN__(32) x,
					       const int32_t incx,
					       float * __restrict  __ATTR_ALIGN__(32) y,
					       const int32_t incy) {
                   if(__builtin_expect(0==n,0) ||
		      __builtin_expect(alpha==0.0f,0)) { return;}
		   
		  
                   __m256 valpha;
		   float * __restrict x0 = NULL;
		   float * __restrict y0 = NULL;
		   int32_t i, last_i;
		   last_i = 0;
		   x0 = x;
		   y0 = y;

		   if(__builtin_expect(incx==1,1) &&
		      __builtin_expect(incy==1,1)) {
		        valpha = _mm256_broadcast_ss(alpha);
#if defined(INTEL_COMPILER) || defined(__ICC)
                      __assume_aligned(x,32);
		      __assume_aligned(y,32);
		      __assume_aligned(x0,32);
		      __assume_aligned(y0,32)
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC)
                      x  = (float*)__builtin_assume_aligned(x,32);
		      y  = (float*)__builtin_assume_aligned(y,32);
		      x0 = (float*)__builtin_assume_aligned(x0,32);
		      y0 = (float*)__builtin_assume_aligned(y0,32);
#endif
                    
		      // Unrolled 10 times in order exploit the ratio
		      // between the FMA latency to its throughput.
#pragma omp parallel for schedule(static,80) default(none) \
                      lastprivate(last_i) private(i) shared(valpha,n,x0,y0)
		      for(i = 0; (i+79) < n; i += 80) {
		          last_i = i;
			  _mm_prefetch((const char*)&x0[i+80],_MM_HINT_T0);
			  _mm256_store_ps(&y0[i+0],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[i+0]),valpha,
				                                _mm256_load_ps(&y0[i+0])));
			  _mm256_store_ps(&y0[i+8],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[i+8]),valpha,
				                                _mm256_load_ps(&y0[i+8])));
			  _mm256_store_ps(&y0[i+16],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[i+16]),valpha,
				                                _mm256_load_ps(&y0[i+16])));
			  _mm256_store_ps(&y0[i+24],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[i+24]),valpha,
				                                _mm256_load_ps(&y0[i+24])));
			  _mm256_store_ps(&y0[i+32],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[i+32]),valpha,
				                                _mm256_load_ps(&y0[i+32])));
			  _mm256_store_ps(&y0[i+40],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[i+40]),valpha,
				                                _mm256_load_ps(&y0[i+40])));
			  _mm256_store_ps(&y0[i+48],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[i+48]),valpha,
				                                _mm256_load_ps(&y0[i+48])));
			  _mm256_store_ps(&y0[i+56],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[i+56]),valpha,
				                                _mm256_load_ps(&y0[i+56])));
			  _mm256_store_ps(&y0[i+64],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[i+64]),valpha,
				                                _mm256_load_ps(&y0[i+64])));
			  _mm256_store_ps(&y0[i+72],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[i+72]),valpha,
				                                _mm256_load_ps(&y0[i+72])));	
                        
		      }
		      
		      for(; (last_i+39) < n; last_i += 40) {

		          _mm256_store_ps(&y0[last_i+0],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[last_i+0]),valpha,
				                                _mm256_load_ps(&y0[last_i+0])));
			  _mm256_store_ps(&y0[last_i+8],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[last_i+8]),valpha,
				                                _mm256_load_ps(&y0[last_i+8])));
			  _mm256_store_ps(&y0[last_i+16],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[last_i+16]),valpha,
				                                _mm256_load_ps(&y0[last_i+16])));
			  _mm256_store_ps(&y0[last_i+24],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[last_i+24]),valpha,
				                                _mm256_load_ps(&y0[last_i+24])));
			  _mm256_store_ps(&y0[last_i+32],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[last_i+32]),valpha,
				                                _mm256_load_ps(&y0[last_i+32])));
			
		       
		      }
		      for(; (last_i+31) < n; last_i += 32) {

		           _mm256_store_ps(&y0[last_i+0],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[last_i+0]),valpha,
				                                _mm256_load_ps(&y0[last_i+0])));
			  _mm256_store_ps(&y0[last_i+8],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[last_i+8]),valpha,
				                                _mm256_load_ps(&y0[last_i+8])));
			  _mm256_store_ps(&y0[last_i+16],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[last_i+16]),valpha,
				                                _mm256_load_ps(&y0[last_i+16])));
			  _mm256_store_ps(&y0[last_i+24],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[last_i+24]),valpha,
				                                _mm256_load_ps(&y0[last_i+24])));
		       
		      }
		      for(; (last_i+15) < n; last_i += 16) {
                          _mm256_store_ps(&y0[last_i+0],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[last_i+0]),valpha,
				                                _mm256_load_ps(&y0[last_i+0])));
			  _mm256_store_ps(&y0[last_i+8],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[last_i+8]),valpha,
				                                _mm256_load_ps(&y0[last_i+8])));
		       
		      }
		      for(; (last_i+7) < n; last_i += 8) {

		           _mm256_store_ps(&y0[last_i+0],
			         _mm256_fmadd_ps(
				           _mm256_load_ps(&x0[last_i+0]),valpha,
				                                _mm256_load_ps(&y0[last_i+0])));
		      }
		      _mm256_vzeroupper();
		      for(; (last_i+0) < n; last_i += 1) {
                          y0[last_i] += alpha * x0[last_i];
		      }
		   }
		   else {

                       for(i = 0; i != n; ++i) {
		           const float tx0 = *x0;
                           *y0 += alpha * tx0;
			   x0 += incx;
			   y0 += incy;
		       }
		   }
	     }


	      /*
                DAXPY kernel unaligned
              */
	  
	      void gms::math::daxpy_u_ymm4r8_unroll10x(const int32_t n,
	                                    const double alpha,
	                                    double * __restrict  x,
					    const int32_t incx,
					    double * __restrict  y,
					    const int32_t incy) {
                   if(__builtin_expect(0==n,0) ||
		      __builtin_expect(alpha==0.0,0)) { return;}
		   
		   __ATTR_ALIGN__(32) __m256d xv[10];
		   __ATTR_ALIGN__(32) __m256d yv[10];
		   __ATTR_ALIGN__(32) __m256d zv[10];
                   __m256d valpha;
		   double * __restrict x0 = NULL;
		   double * __restrict y0 = NULL;
		   int32_t i;
		   x0 = x;
		   y0 = y;

		   if(__builtin_expect(incx==1,1) &&
		      __builtin_expec(incy==1,1)) {

                      valpha = _mm256_broadcast_sd(alpha);
		      // Unrolled 10 times in order exploit the ratio
		      // between the FMA latency to its throughput.
		      for(i = 0; (i+39) < n; i += 40) {
                          xv[0] = _mm256_loadu_pd(&x0[i+0]);
			  yv[0] = _mm256_loadu_pd(&y0[i+0]);
			  zv[0] = _mm256_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm256_storeu_pd(&y0[i+0],zv[0]);
			  xv[1] = _mm256_loadu_pd(&x0[i+4]);
			  yv[1] = _mm256_loadu_pd(&y0[i+4]);
			  zv[1] = _mm256_fmadd_pd(xv[1],valpha,yv[1]);
			  _mm256_storeu_pd(&y0[i+4],zv[1]);
			  xv[2] = _mm256_loadu_pd(&x0[i+8]);
			  yv[2] = _mm256_loadu_pd(&y0[i+8]);
			  zv[2] = _mm256_fmadd_pd(xv[2],valpha,yv[2]);
			  _mm256_storeu_pd(&y0[i+8],zv[2]);
			  xv[3] = _mm256_loadu_pd(&x0[i+12]);
			  yv[3] = _mm256_loadu_pd(&y0[i+12]);
			  zv[3] = _mm256_fmadd_pd(xv[3],valpha,yv[3]);
			  _mm256_storeu_ps(&y0[i+12],zv[3]);
			  xv[4] = _mm256_loadu_pd(&x0[i+16]);
			  yv[4] = _mm256_loadu_pd(&y0[i+16]);
			  zv[4] = _mm256_fmadd_pd(xv[4],valpha,yv[4]);
			  _mm256_storeu_pd(&y0[i+16],zv[4]);
			  xv[5] = _mm256_loadu_pd(&x0[i+20]);
			  yv[5] = _mm256_loadu_pd(&y0[i+20]);
			  zv[5] = _mm256_fmadd_pd(xv[5],valpha,yv[5]);
			  _mm256_storeu_pd(&y0[i+20],zv[5]);
			  xv[6] = _mm256_loadu_pd(&x0[i+24]);
			  yv[6] = _mm256_loadu_pd(&y0[i+24]);
			  zv[6] = _mm256_fmadd_pd(xv[6],valpha,yv[6]);
			  _mm256_storeu_pd(&y0[i+24],zv[6]);
			  xv[7] = _mm256_loadu_pd(&x0[i+28]);
			  yv[7] = _mm256_loadu_pd(&y0[i+28]);
			  zv[7] = _mm256_fmadd_pd(xv[7],valpha,yv[7]);
			  _mm256_storeu_pd(&y0[i+28],zv[7]);
			  xv[8] = _mm256_loadu_pd(&x0[i+32]);
			  yv[8] = _mm256_loadu_pd(&y0[i+32]);
			  zv[8] = _mm256_fmadd_pd(xv[8],valpha,yv[8]);
			  _mm256_storeu_pd(&y0[i+32],zv[8]);
			  xv[9] = _mm256_loadu_pd(&x0[i+36]);
			  yv[9] = _mm256_loadu_pd(&y0[i+36]);
			  zv[9] = _mm256_fmadd_pd(xv[9],valpha,yv[9]);
			  _mm256_storeu_pd(&y0[i+36],zv[9]);
		      }
		      for(; (i+19) < n; i += 20) {

		          xv[0] = _mm256_loadu_pd(&x0[i+0]);
			  yv[0] = _mm256_loadu_pd(&y0[i+0]);
			  zv[0] = _mm256_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm256_storeu_pd(&y0[i+0],zv[0]);
			  xv[1] = _mm256_loadu_pd(&x0[i+4]);
			  yv[1] = _mm256_loadu_pd(&y0[i+4]);
			  zv[1] = _mm256_fmadd_pd(xv[1],valpha,yv[1]);
			  _mm256_storeu_pd(&y0[i+4],zv[1]);
			  xv[2] = _mm256_loadu_pd(&x0[i+8]);
			  yv[2] = _mm256_loadu_pd(&y0[i+8]);
			  zv[2] = _mm256_fmadd_pd(xv[2],valpha,yv[2]);
			  _mm256_storeu_pd(&y0[i+8],zv[2]);
			  xv[3] = _mm256_loadu_pd(&x0[i+12]);
			  yv[3] = _mm256_loadu_pd(&y0[i+12]);
			  zv[3] = _mm256_fmadd_pd(xv[3],valpha,yv[3]);
			  _mm256_storeu_ps(&y0[i+12],zv[3]);
			  xv[4] = _mm256_loadu_pd(&x0[i+16]);
			  yv[4] = _mm256_loadu_pd(&y0[i+16]);
			  zv[4] = _mm256_fmadd_pd(xv[4],valpha,yv[4]);
			  _mm256_storeu_pd(&y0[i+16],zv[4]);
		      }
		      for(; (i+15) < n; i += 16) {

		          xv[0] = _mm256_loadu_pd(&x0[i+0]);
			  yv[0] = _mm256_loadu_pd(&y0[i+0]);
			  zv[0] = _mm256_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm256_storeu_pd(&y0[i+0],zv[0]);
			  xv[1] = _mm256_loadu_pd(&x0[i+4]);
			  yv[1] = _mm256_loadu_pd(&y0[i+4]);
			  zv[1] = _mm256_fmadd_pd(xv[1],valpha,yv[1]);
			  _mm256_storeu_pd(&y0[i+4],zv[1]);
			  xv[2] = _mm256_loadu_pd(&x0[i+8]);
			  yv[2] = _mm256_loadu_pd(&y0[i+8]);
			  zv[2] = _mm256_fmadd_pd(xv[2],valpha,yv[2]);
			  _mm256_storeu_pd(&y0[i+8],zv[2]);
			  xv[3] = _mm256_loadu_pd(&x0[i+12]);
			  yv[3] = _mm256_loadu_pd(&y0[i+12]);
			  zv[3] = _mm256_fmadd_pd(xv[3],valpha,yv[3]);
			  _mm256_storeu_ps(&y0[i+12],zv[3]);
		      }
		      for(; (i+7) < n; i += 8) {

		          xv[0] = _mm256_loadu_pd(&x0[i+0]);
			  yv[0] = _mm256_loadu_pd(&y0[i+0]);
			  zv[0] = _mm256_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm256_storeu_pd(&y0[i+0],zv[0]);
			  xv[1] = _mm256_loadu_pd(&x0[i+4]);
			  yv[1] = _mm256_loadu_pd(&y0[i+4]);
			  zv[1] = _mm256_fmadd_pd(xv[1],valpha,yv[1]);
			  _mm256_storeu_pd(&y0[i+4],zv[1]);
		      }
		      for(; (i+3) < n; i += 4) {

		          xv[0] = _mm256_loadu_pd(&x0[i+0]);
			  yv[0] = _mm256_loadu_pd(&y0[i+0]);
			  zv[0] = _mm256_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm256_storeu_pd(&y0[i+0],zv[0]);
		      }
		      _mm256_vzeroupper();
		      for(; (i+0) < n; i += 1) {
                          y0[i] += alpha * x0[i];
		      }
		   }
		   else {

                       for(i = 0; i != n; ++i) {
		           const double tx0 = *x0;
                           *y0 += alpha * tx0;
			   x0 += incx;
			   y0 += incy;
		       }
		   }
	     }


	    /*
                DAXPY kernel aligned
            */
	     
	      void gms::math::daxpy_a_ymm4r8_unroll10x(const int32_t n,
	                                    const double alpha,
	                                    double * __restrict  __ATTR_ALIGN__(32) x,
					    const int32_t incx,
					    double * __restrict  __ATTR_ALIGN__(32) y,
					    const int32_t incy) {
                   if(__builtin_expect(0==n,0) ||
		      __builtin_expect(alpha==0.0,0)) { return;}
		   
		   __ATTR_ALIGN__(32) __m256d xv[10];
		   __ATTR_ALIGN__(32) __m256d yv[10];
		   __ATTR_ALIGN__(32) __m256d zv[10];
                   __m256d valpha;
		   double * __restrict x0 = NULL;
		   double * __restrict y0 = NULL;
		   int32_t i;
		   x0 = x;
		   y0 = y;

		   if(__builtin_expect(incx==1,1) &&
		      __builtin_expect(incy==1,1)) {
		      valpha = _mm256_broadcast_sd(alpha);
#if defined(INTEL_COMPILER) || defined(__ICC)
                      __assume_aligned(x,32);
		      __assume_aligned(y,32);
		      __assume_aligned(x0,32);
		      __assume_aligned(y0,32);
#pragma code_align(64)

#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC)
                      x  = (double*)__builtin_assume_aligned(x,32);
		      y  = (double*)__builtin_assume_aligned(y,32);
		      x0 = (double*)__builtin_assume_aligned(x0,32);
		      y0 = (double*)__builtin_assume_aligned(y0,32);
#endif
                      
		      // Unrolled 10 times in order exploit the ratio
		      // between the FMA latency to its throughput.
		      for(i = 0; (i+39) < n; i += 40) {
		          _mm_prefetch((const char*)&x[i+40],_MM_HINT_T0);
                          xv[0] = _mm256_load_pd(&x0[i+0]);
			  yv[0] = _mm256_load_pd(&y0[i+0]);
			  zv[0] = _mm256_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm256_store_pd(&y0[i+0],zv[0]);
			  xv[1] = _mm256_load_pd(&x0[i+4]);
			  yv[1] = _mm256_load_pd(&y0[i+4]);
			  zv[1] = _mm256_fmadd_pd(xv[1],valpha,yv[1]);
			  _mm256_store_pd(&y0[i+4],zv[1]);
			  xv[2] = _mm256_load_pd(&x0[i+8]);
			  yv[2] = _mm256_load_pd(&y0[i+8]);
			  zv[2] = _mm256_fmadd_pd(xv[2],valpha,yv[2]);
			  _mm256_store_pd(&y0[i+8],zv[2]);
			  xv[3] = _mm256_load_pd(&x0[i+12]);
			  yv[3] = _mm256_load_pd(&y0[i+12]);
			  zv[3] = _mm256_fmadd_pd(xv[3],valpha,yv[3]);
			  _mm256_store_ps(&y0[i+12],zv[3]);
			  xv[4] = _mm256_load_pd(&x0[i+16]);
			  yv[4] = _mm256_load_pd(&y0[i+16]);
			  zv[4] = _mm256_fmadd_pd(xv[4],valpha,yv[4]);
			  _mm256_store_pd(&y0[i+16],zv[4]);
			  xv[5] = _mm256_load_pd(&x0[i+20]);
			  yv[5] = _mm256_load_pd(&y0[i+20]);
			  zv[5] = _mm256_fmadd_pd(xv[5],valpha,yv[5]);
			  _mm256_store_pd(&y0[i+20],zv[5]);
			  xv[6] = _mm256_load_pd(&x0[i+24]);
			  yv[6] = _mm256_load_pd(&y0[i+24]);
			  zv[6] = _mm256_fmadd_pd(xv[6],valpha,yv[6]);
			  _mm256_store_pd(&y0[i+24],zv[6]);
			  xv[7] = _mm256_load_pd(&x0[i+28]);
			  yv[7] = _mm256_load_pd(&y0[i+28]);
			  zv[7] = _mm256_fmadd_pd(xv[7],valpha,yv[7]);
			  _mm256_storeu_pd(&y0[i+28],zv[7]);
			  xv[8] = _mm256_load_pd(&x0[i+32]);
			  yv[8] = _mm256_load_pd(&y0[i+32]);
			  zv[8] = _mm256_fmadd_pd(xv[8],valpha,yv[8]);
			  _mm256_store_pd(&y0[i+32],zv[8]);
			  xv[9] = _mm256_load_pd(&x0[i+36]);
			  yv[9] = _mm256_load_pd(&y0[i+36]);
			  zv[9] = _mm256_fmadd_pd(xv[9],valpha,yv[9]);
			  _mm256_store_pd(&y0[i+36],zv[9]);
		      }
		      for(; (i+19) < n; i += 20) {

		          xv[0] = _mm256_load_pd(&x0[i+0]);
			  yv[0] = _mm256_load_pd(&y0[i+0]);
			  zv[0] = _mm256_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm256_store_pd(&y0[i+0],zv[0]);
			  xv[1] = _mm256_load_pd(&x0[i+4]);
			  yv[1] = _mm256_load_pd(&y0[i+4]);
			  zv[1] = _mm256_fmadd_pd(xv[1],valpha,yv[1]);
			  _mm256_store_pd(&y0[i+4],zv[1]);
			  xv[2] = _mm256_load_pd(&x0[i+8]);
			  yv[2] = _mm256_load_pd(&y0[i+8]);
			  zv[2] = _mm256_fmadd_pd(xv[2],valpha,yv[2]);
			  _mm256_store_pd(&y0[i+8],zv[2]);
			  xv[3] = _mm256_load_pd(&x0[i+12]);
			  yv[3] = _mm256_load_pd(&y0[i+12]);
			  zv[3] = _mm256_fmadd_pd(xv[3],valpha,yv[3]);
			  _mm256_store_ps(&y0[i+12],zv[3]);
			  xv[4] = _mm256_load_pd(&x0[i+16]);
			  yv[4] = _mm256_load_pd(&y0[i+16]);
			  zv[4] = _mm256_fmadd_pd(xv[4],valpha,yv[4]);
			  _mm256_store_pd(&y0[i+16],zv[4]);
		      }
		      for(; (i+15) < n; i += 16) {

		          xv[0] = _mm256_load_pd(&x0[i+0]);
			  yv[0] = _mm256_load_pd(&y0[i+0]);
			  zv[0] = _mm256_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm256_store_pd(&y0[i+0],zv[0]);
			  xv[1] = _mm256_load_pd(&x0[i+4]);
			  yv[1] = _mm256_load_pd(&y0[i+4]);
			  zv[1] = _mm256_fmadd_pd(xv[1],valpha,yv[1]);
			  _mm256_store_pd(&y0[i+4],zv[1]);
			  xv[2] = _mm256_load_pd(&x0[i+8]);
			  yv[2] = _mm256_load_pd(&y0[i+8]);
			  zv[2] = _mm256_fmadd_pd(xv[2],valpha,yv[2]);
			  _mm256_store_pd(&y0[i+8],zv[2]);
			  xv[3] = _mm256_load_pd(&x0[i+12]);
			  yv[3] = _mm256_load_pd(&y0[i+12]);
			  zv[3] = _mm256_fmadd_pd(xv[3],valpha,yv[3]);
			  _mm256_store_ps(&y0[i+12],zv[3]);
		      }
		      for(; (i+7) < n; i += 8) {

		          xv[0] = _mm256_load_pd(&x0[i+0]);
			  yv[0] = _mm256_load_pd(&y0[i+0]);
			  zv[0] = _mm256_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm256_store_pd(&y0[i+0],zv[0]);
			  xv[1] = _mm256_load_pd(&x0[i+4]);
			  yv[1] = _mm256_load_pd(&y0[i+4]);
			  zv[1] = _mm256_fmadd_pd(xv[1],valpha,yv[1]);
			  _mm256_store_pd(&y0[i+4],zv[1]);
		      }
		      for(; (i+3) < n; i += 4) {

		          xv[0] = _mm256_load_pd(&x0[i+0]);
			  yv[0] = _mm256_load_pd(&y0[i+0]);
			  zv[0] = _mm256_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm256_store_pd(&y0[i+0],zv[0]);
		      }
		      _mm256_vzeroupper();
		      for(; (i+0) < n; i += 1) {
                          y0[i] += alpha * x0[i];
		      }
		   }
		   else {

                       for(i = 0; i != n; ++i) {
		           const double tx0 = *x0;
                           *y0 += alpha * tx0;
			   x0 += incx;
			   y0 += incy;
		       }
		   }
	     }



	       /*
                DAXPY kernel OpenMP parallelized aligned
            */
	    
	      void gms::math::daxpy_a_ymm4r8_unroll10x_omp(const int32_t n,
	                                        const double alpha,
	                                        double * __restrict  __ATTR_ALIGN__(32) x,
					        const int32_t incx,
					        double * __restrict  __ATTR_ALIGN__(32) y,
					        const int32_t incy) {
                   if(__builtin_expect(0==n,0) ||
		      __builtin_expect(alpha==0.0)) { return;}
		   
		   __ATTR_ALIGN__(32) __m256d xv[10];
		   __ATTR_ALIGN__(32) __m256d yv[10];
		   __ATTR_ALIGN__(32) __m256d zv[10];
                   __m256d valpha;
		   double * __restrict x0 = NULL;
		   double * __restrict y0 = NULL;
		   int32_t i, last_i;
		   last_i = 0;
		   x0 = x;
		   y0 = y;

		   if(__builtin_expect(incx==1,1) &&
		      __builtin_expect(incy==1,1)) {
		       valpha = _mm256_broadcast_sd(alpha);
#if defined(INTEL_COMPILER) || defined(__ICC)
                      __assume_aligned(x,32);
		      __assume_aligned(y,32);
		      __assume_aligned(x0,32);
		      __assume_aligned(y0,32);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC)
                      x = (double*)__builtin_assume_aligned(x,32);
		      y = (double*)__builtin_assume_aligned(y,32);
		      x0 = (double*)__builtin_assume_aligned(x0,32);
		      y0 = (double*)__builtin_assume_aligned(y0,32);
#endif
                      valpha = _mm256_broadcast_sd(alpha);
		      // Unrolled 10 times in order exploit the ratio
		      // between the FMA latency to its throughput.
#pragma omp parallel for schedule(static,40) default(none)  \
        lastprivate(last_i) private(i) shared(valpha,n,x0,y0)
		      for(i = 0; (i+39) < n; i += 40) {
		          last_i = i;
		          _mm256_store_pd(&y0[i+0],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[i+0]),valpha,
					                       _mm256_load_pd(&y[i+0])));
			 _mm256_store_pd(&y0[i+4],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[i+4]),valpha,
					                       _mm256_load_pd(&y[i+4])));
			 _mm256_store_pd(&y0[i+8],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[i+8]),valpha,
					                       _mm256_load_pd(&y[i+8])));
			 _mm256_store_pd(&y0[i+12],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[i+12]),valpha,
					                       _mm256_load_pd(&y[i+12])));
			 _mm256_store_pd(&y0[i+16],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[i+16]),valpha,
					                       _mm256_load_pd(&y[i+16])));
			 _mm256_store_pd(&y0[i+20],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[i+20]),valpha,
					                       _mm256_load_pd(&y[i+20])));
			 _mm256_store_pd(&y0[i+24],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[i+24]),valpha,
					                       _mm256_load_pd(&y[i+24])));
			 _mm256_store_pd(&y0[i+28],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[i+28]),valpha,
					                       _mm256_load_pd(&y[i+28])));
			 _mm256_store_pd(&y0[i+32],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[i+32]),valpha,
					                       _mm256_load_pd(&y[i+32])));
			 _mm256_store_pd(&y0[i+36],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[i+36]),valpha,
					                       _mm256_load_pd(&y[i+36])));
                         
		      }
		      for(; (last_i+19) < n; last_i += 20) {

		         _mm256_store_pd(&y0[last_i+0],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[last_i+0]),valpha,
					                       _mm256_load_pd(&y[last_i+0])));
			 _mm256_store_pd(&y0[last_i+4],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[last_i+4]),valpha,
					                       _mm256_load_pd(&y[last_i+4])));
			 _mm256_store_pd(&y0[last_i+8],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[last_i+8]),valpha,
					                       _mm256_load_pd(&y[last_i+8])));
			 _mm256_store_pd(&y0[last_i+12],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[last_i+12]),valpha,
					                       _mm256_load_pd(&y[last_i+12])));
			 _mm256_store_pd(&y0[last_i+16],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[last_i+16]),valpha,
					                       _mm256_load_pd(&y[last_i+16])));

		      }
		      for(; (last_i+15) < n; last_i += 16) {

		         _mm256_store_pd(&y0[last_i+0],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[last_i+0]),valpha,
					                       _mm256_load_pd(&y[last_i+0])));
			 _mm256_store_pd(&y0[last_i+4],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[last_i+4]),valpha,
					                       _mm256_load_pd(&y[last_i+4])));
			 _mm256_store_pd(&y0[last_i+8],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[last_i+8]),valpha,
					                       _mm256_load_pd(&y[last_i+8])));
			 _mm256_store_pd(&y0[last_i+12],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[last_i+12]),valpha,
					                       _mm256_load_pd(&y[last_i+12])));
		        
		      }
		      for(; (last_i+7) < n; last_i += 8) {

		         _mm256_store_pd(&y0[last_i+0],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[last_i+0]),valpha,
					                       _mm256_load_pd(&y[last_i+0])));
			 _mm256_store_pd(&y0[last_i+4],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[last_i+4]),valpha,
					                       _mm256_load_pd(&y[last_i+4])));
		        
		      }
		      for(; (last_i+3) < n; last_i += 4) {

		          _mm256_store_pd(&y0[last_i+0],
			          _mm256_fmadd_pd(
				            _mm256_load_pd(&x0[last_i+0]),valpha,
					                       _mm256_load_pd(&y[last_i+0])));
		      }
		      _mm256_vzeroupper();
		      for(; (last_i+0) < n; last_i += 1) {
                          y0[last_i] += alpha * x0[last_i];
		      }
		   }
		   else {

                       for(i = 0; i != n; ++i) {
		           const double tx0 = *x0;
                           *y0 += alpha * tx0;
			   x0 += incx;
			   y0 += incy;
		       }
		   }
	     }


	     


	     
	     
  
