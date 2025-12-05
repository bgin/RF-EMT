



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



#include "GMS_axpy_avx512_unroll10.h"

	        /*
                   SAXPY kernel unaligned.
               */
	          
                   void gms::math::saxpy_u_zmm16r4_unroll10x(const int32_t n,
		                                  const float alpha,
		                                  float * __restrict x,
                                                  const int32_t incx,
						  float * __restrict y,
						  const int32_t incy) {
                         if(__builtin_expect(0==n,0) ||
			    __builtin_expect(0.0==alpha,0)) { return;}
			    
			 __ATTR_ALIGN__(64) __m512 xv[10];
			 __ATTR_ALIGN__(64) __m512 yv[10];
			 __ATTR_ALIGN__(64) __m512 zv[10];
			 __m512 valpha;
			 float * __restrict x0 = NULL;
			 float * __restrict y0 = NULL;
			 int32_t i;
			 x0 = x;
			 y0 = y;

			 if(__builtin_expect(incx==1,1) &&
			    __builtin_expect(incy==1,1)) {

			    valpha = _mm512_broadcast_ss(alpha);
			     // Unrolled 10 times in order exploit the ratio
		             // between the FMA latency to its throughput.
			     for(i = 0; (i+159) < n; i += 160) {
                                 _mm_prefetch((const char*)&x0[i+0],_MM_HINT_T0);
				 xv[0] = _mm512_loadu_ps(&x0[i+0]);
				 yv[0] = _mm512_loadu_ps(&y0[i+0]);
				 zv[0] = _mm512_fmadd_ps(xv[0],valpha,yv[0]);
				 _mm512_storeu_ps(&y0[i+0],zv[0]);
				 xv[1] = _mm512_loadu_ps(&x0[i+16]);
				 yv[1] = _mm512_loadu_ps(&y0[i+16]);
				 zv[1] = _mm512_fmadd_ps(xv[1],valpha,yv[1]);
				 _mm512_storeu_ps(&y0[i+16],zv[1]);
				 _mm_prefetch((const char*)&x0[i+32],_MM_HINT_T0);
				 xv[2] = _mm512_loadu_ps(&x0[i+32]);
				 yv[2] = _mm512_loadu_ps(&y0[i+32]);
				 zv[2] = _mm512_fmadd_ps(xv[2],valpha,yv[2]);
				 _mm512_storeu_ps(&y0[i+32],zv[2]);
				 xv[3] = _mm512_loadu_ps(&x0[i+48]);
				 yv[3] = _mm512_loadu_ps(&y0[i+48]);
				 zv[3] = _mm512_fmadd_ps(xv[3],valpha,yv[3]);
				 _mm512_storeu_ps(&y0[i+48],zv[3]);
				 _mm_prefetch((const char*)&x0[i+64],_MM_HINT_T0);
				 xv[4] = _mm512_loadu_ps(&x0[i+64]);
				 yv[4] = _mm512_loadu_ps(&y0[i+64]);
				 zv[4] = _mm512_fmadd_ps(xv[4],valpha,yv[4]);
				 _mm512_storeu_ps(&y0[i+64],zv[4]);
				 xv[5] = _mm512_loadu_ps(&x0[i+80]);
				 yv[5] = _mm512_loadu_ps(&y0[i+80]);
				 zv[5] = _mm512_fmadd_ps(xv[5],valpha,yv[5]);
				 _mm512_storeu_ps(&y0[i+80],zv[5]);
				 _mm_prefetch((const char*)&x0[i+96],_MM_HINT_T0);
				 xv[6] = _mm512_loadu_ps(&x0[i+96]);
				 yv[6] = _mm512_loadu_ps(&y0[i+96]);
				 zv[6] = _mm512_fmadd_ps(xv[6],valpha,yv[6]);
				 _mm512_storeu_ps(&y0[i+96],zv[6]);
				 xv[7] = _mm512_loadu_ps(&x0[i+112]);
				 yv[7] = _mm512_loadu_ps(&y0[i+112]);
				 zv[7] = _mm512_fmadd_ps(xv[7],valpha,yv[7]);
				 _mm512_storeu_ps(&y0[i+112],zv[7]);
				 _mm_prefetch((const char*)&x0[i+128],_MM_HINT_T0);
				 xv[8] = _mm512_loadu_ps(&x0[i+128]);
				 yv[8] = _mm512_loadu_ps(&y0[i+128]);
				 zv[8] = _mm512_fmadd_ps(xv[8],valpha,yv[8]);
				 _mm512_storeu_ps(&y0[i+128],zv[8]);
				 xv[9] = _mm512_loadu_ps(&x0[i+144]);
				 yv[9] = _mm512_loadu_ps(&y0[i+144]);
				 zv[9] = _mm512_fmadd_ps(xv[9],valpha,yv[9]);
				 _mm512_storeu_ps(&y0[i+144],zv[9]);
			     }
			    for(; (i+79) < n; i += 80) {
                                 //_mm_prefetch((const char*)&x[i+0],_MM_HINT_T0);
				 xv[0] = _mm512_loadu_ps(&x0[i+0]);
				 yv[0] = _mm512_loadu_ps(&y0[i+0]);
				 zv[0] = _mm512_fmadd_ps(xv[0],valpha,yv[0]);
				 _mm512_storeu_ps(&y0[i+0],zv[0]);
				 xv[1] = _mm512_loadu_ps(&x0[i+16]);
				 yv[1] = _mm512_loadu_ps(&y0[i+16]);
				 zv[1] = _mm512_fmadd_ps(xv[1],valpha,yv[1]);
				 _mm512_storeu_ps(&y0[i+16],zv[1]);
				// _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				 xv[2] = _mm512_loadu_ps(&x0[i+32]);
				 yv[2] = _mm512_loadu_ps(&y0[i+32]);
				 zv[2] = _mm512_fmadd_ps(xv[2],valpha,yv[2]);
				 _mm512_storeu_ps(&y0[i+32],zv[2]);
				 xv[3] = _mm512_loadu_ps(&x0[i+48]);
				 yv[3] = _mm512_loadu_ps(&y0[i+48]);
				 zv[3] = _mm512_fmadd_ps(xv[3],valpha,yv[3]);
				 _mm512_storeu_ps(&y0[i+48],zv[3]);
				// _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				 xv[4] = _mm512_loadu_ps(&x0[i+64]);
				 yv[4] = _mm512_loadu_ps(&y0[i+64]);
				 zv[4] = _mm512_fmadd_ps(xv[4],valpha,yv[4]);
				 _mm512_storeu_ps(&y0[i+64],zv[4]);
			    }
			   for(; (i+63) < n; i += 64) {
                                 //_mm_prefetch((const char*)&x[i+0],_MM_HINT_T0);
				 xv[0] = _mm512_loadu_ps(&x0[i+0]);
				 yv[0] = _mm512_loadu_ps(&y0[i+0]);
				 zv[0] = _mm512_fmadd_ps(xv[0],valpha,yv[0]);
				 _mm512_storeu_ps(&y0[i+0],zv[0]);
				 xv[1] = _mm512_loadu_ps(&x0[i+16]);
				 yv[1] = _mm512_loadu_ps(&y0[i+16]);
				 zv[1] = _mm512_fmadd_ps(xv[1],valpha,yv[1]);
				 _mm512_storeu_ps(&y0[i+16],zv[1]);
				// _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				 xv[2] = _mm512_loadu_ps(&x0[i+32]);
				 yv[2] = _mm512_loadu_ps(&y0[i+32]);
				 zv[2] = _mm512_fmadd_ps(xv[2],valpha,yv[2]);
				 _mm512_storeu_ps(&y0[i+32],zv[2]);
				 xv[3] = _mm512_loadu_ps(&x0[i+48]);
				 yv[3] = _mm512_loadu_ps(&y0[i+48]);
				 zv[3] = _mm512_fmadd_ps(xv[3],valpha,yv[3]);
				 _mm512_storeu_ps(&y0[i+48],zv[3]);
			   }
			  for(; (i+31) < n; i += 32) {
                                //_mm_prefetch((const char*)&x[i+0],_MM_HINT_T0);
				 xv[0] = _mm512_loadu_ps(&x0[i+0]);
				 yv[0] = _mm512_loadu_ps(&y0[i+0]);
				 zv[0] = _mm512_fmadd_ps(xv[0],valpha,yv[0]);
				 _mm512_storeu_ps(&y0[i+0],zv[0]);
				 xv[1] = _mm512_loadu_ps(&x0[i+16]);
				 yv[1] = _mm512_loadu_ps(&y0[i+16]);
				 zv[1] = _mm512_fmadd_ps(xv[1],valpha,yv[1]);
				 _mm512_storeu_ps(&y0[i+16],zv[1]);
			  }
			 for(; (i+15) < n; i += 16) {
                               //_mm_prefetch((const char*)&x[i+0],_MM_HINT_T0);
				 xv[0] = _mm512_loadu_ps(&x0[i+0]);
				 yv[0] = _mm512_loadu_ps(&y0[i+0]);
				 zv[0] = _mm512_fmadd_ps(xv[0],valpha,yv[0]);
				 _mm512_storeu_ps(&y0[i+0],zv[0]);
 			 }
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
                   
                   void gms::math::saxpy_a_zmm16r4_unroll10x(const int32_t n,
		                                  const float alpha,
		                                  float * __restrict __ATTR_ALIGN__(64) x,
                                                  const int32_t incx,
						  float * __restrict __ATTR_ALIGN__(64) y,
						  const int32_t incy) {
                         if(__builtin_expect(0==n,0) ||
			    __builtin_expect(0.0f == alpha,0)) { return;}
			    
			 __ATTR_ALIGN__(64) __m512 xv[10];
			 __ATTR_ALIGN__(64) __m512 yv[10];
			 __ATTR_ALIGN__(64) __m512 zv[10];
			 __m512 valpha;
			 float * __restrict x0 = NULL;
			 float * __restrict y0 = NULL;
			 int32_t i;
			 x0 = x;
			 y0 = y;

			 if(__builtin_expect(incx==1,1) &&
			    __builtin_expect(incy==1,1)) {

			    valpha = _mm512_broadcast_ss(alpha);
			     // Unrolled 10 times in order exploit the ratio
		             // between the FMA latency to its throughput.
#if defined(INTEL_COMPILER) || defined(__ICC)
                      __assume_aligned(x,64);
		      __assume_aligned(y,64);
		      __assume_aligned(x0,64);
		      __assume_aligned(y0,64);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC)
                      x  = (float*)__builtin_assume_aligned(x,64);
		      y  = (float*)__builtin_assume_aligned(y,64);
		      x0 = (float*)__builtin_assume_aligned(x0,64);
		      y0 = (float*)__builtin_assume_aligned(y0,64);
#endif			     
			     for(i = 0; (i+159) < n; i += 160) {
                                 _mm_prefetch((const char*)&x0[i+0],_MM_HINT_T0);
				 xv[0] = _mm512_load_ps(&0[i+0]);
				 yv[0] = _mm512_load_ps(&y0[i+0]);
				 zv[0] = _mm512_fmadd_ps(xv[0],valpha,yv[0]);
				 _mm512_store_ps(&y0[i+0],zv[0]);
				 xv[1] = _mm512_load_ps(&x0[i+16]);
				 yv[1] = _mm512_load_ps(&y0[i+16]);
				 zv[1] = _mm512_fmadd_ps(xv[1],valpha,yv[1]);
				 _mm512_store_ps(&y0[i+16],zv[1]);
				 _mm_prefetch((const char*)&x0[i+32],_MM_HINT_T0);
				 xv[2] = _mm512_load_ps(&x0[i+32]);
				 yv[2] = _mm512_load_ps(&y0[i+32]);
				 zv[2] = _mm512_fmadd_ps(xv[2],valpha,yv[2]);
				 _mm512_store_ps(&y0[i+32],zv[2]);
				 xv[3] = _mm512_load_ps(&x0[i+48]);
				 yv[3] = _mm512_load_ps(&y0[i+48]);
				 zv[3] = _mm512_fmadd_ps(xv[3],valpha,yv[3]);
				 _mm512_store_ps(&y0[i+48],zv[3]);
				 _mm_prefetch((const char*)&x0[i+64],_MM_HINT_T0);
				 xv[4] = _mm512_load_ps(&x0[i+64]);
				 yv[4] = _mm512_load_ps(&y0[i+64]);
				 zv[4] = _mm512_fmadd_ps(xv[4],valpha,yv[4]);
				 _mm512_store_ps(&y0[i+64],zv[4]);
				 xv[5] = _mm512_load_ps(&x0[i+80]);
				 yv[5] = _mm512_load_ps(&y0[i+80]);
				 zv[5] = _mm512_fmadd_ps(xv[5],valpha,yv[5]);
				 _mm512_store_ps(&y0[i+80],zv[5]);
				 _mm_prefetch((const char*)&x0[i+96],_MM_HINT_T0);
				 xv[6] = _mm512_load_ps(&x0[i+96]);
				 yv[6] = _mm512_load_ps(&y0[i+96]);
				 zv[6] = _mm512_fmadd_ps(xv[6],valpha,yv[6]);
				 _mm512_store_ps(&y0[i+96],zv[6]);
				 xv[7] = _mm512_load_ps(&x0[i+112]);
				 yv[7] = _mm512_load_ps(&y0[i+112]);
				 zv[7] = _mm512_fmadd_ps(xv[7],valpha,yv[7]);
				 _mm512_store_ps(&y0[i+112],zv[7]);
				 _mm_prefetch((const char*)&x0[i+128],_MM_HINT_T0);
				 xv[8] = _mm512_load_ps(&x0[i+128]);
				 yv[8] = _mm512_load_ps(&0y[i+128]);
				 zv[8] = _mm512_fmadd_ps(xv[8],valpha,yv[8]);
				 _mm512_store_ps(&y0[i+128],zv[8]);
				 xv[9] = _mm512_load_ps(&x0[i+144]);
				 yv[9] = _mm512_load_ps(&y0[i+144]);
				 zv[9] = _mm512_fmadd_ps(xv[9],valpha,yv[9]);
				 _mm512_store_ps(&y0[i+144],zv[9]);
			     }
			    for(; (i+79) < n; i += 80) {
                                 //_mm_prefetch((const char*)&x[i+0],_MM_HINT_T0);
				 xv[0] = _mm512_load_ps(&x0[i+0]);
				 yv[0] = _mm512_load_ps(&y0[i+0]);
				 zv[0] = _mm512_fmadd_ps(xv[0],valpha,yv[0]);
				 _mm512_store_ps(&y0[i+0],zv[0]);
				 xv[1] = _mm512_load_ps(&x0[i+16]);
				 yv[1] = _mm512_load_ps(&y0[i+16]);
				 zv[1] = _mm512_fmadd_ps(xv[1],valpha,yv[1]);
				 _mm512_storeu_ps(&y0[i+16],zv[1]);
				// _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				 xv[2] = _mm512_load_ps(&x0[i+32]);
				 yv[2] = _mm512_load_ps(&y0[i+32]);
				 zv[2] = _mm512_fmadd_ps(xv[2],valpha,yv[2]);
				 _mm512_store_ps(&y0[i+32],zv[2]);
				 xv[3] = _mm512_load_ps(&x0[i+48]);
				 yv[3] = _mm512_load_ps(&y0[i+48]);
				 zv[3] = _mm512_fmadd_ps(xv[3],valpha,yv[3]);
				 _mm512_store_ps(&y0[i+48],zv[3]);
				// _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				 xv[4] = _mm512_load_ps(&x0[i+64]);
				 yv[4] = _mm512_load_ps(&y0[i+64]);
				 zv[4] = _mm512_fmadd_ps(xv[4],valpha,yv[4]);
				 _mm512_store_ps(&y0[i+64],zv[4]);
			    }
			   for(; (i+63) < n; i += 64) {
                                 //_mm_prefetch((const char*)&x[i+0],_MM_HINT_T0);
				 xv[0] = _mm512_load_ps(&x0[i+0]);
				 yv[0] = _mm512_load_ps(&y0[i+0]);
				 zv[0] = _mm512_fmadd_ps(xv[0],valpha,yv[0]);
				 _mm512_store_ps(&y0[i+0],zv[0]);
				 xv[1] = _mm512_load_ps(&x0[i+16]);
				 yv[1] = _mm512_load_ps(&y0[i+16]);
				 zv[1] = _mm512_fmadd_ps(xv[1],valpha,yv[1]);
				 _mm512_store_ps(&y0[i+16],zv[1]);
				// _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				 xv[2] = _mm512_load_ps(&x0[i+32]);
				 yv[2] = _mm512_load_ps(&y0[i+32]);
				 zv[2] = _mm512_fmadd_ps(xv[2],valpha,yv[2]);
				 _mm512_store_ps(&y0[i+32],zv[2]);
				 xv[3] = _mm512_load_ps(&x0[i+48]);
				 yv[3] = _mm512_load_ps(&y0[i+48]);
				 zv[3] = _mm512_fmadd_ps(xv[3],valpha,yv[3]);
				 _mm512_store_ps(&y0[i+48],zv[3]);
			   }
			  for(; (i+31) < n; i += 32) {
                                //_mm_prefetch((const char*)&x[i+0],_MM_HINT_T0);
				 xv[0] = _mm512_load_ps(&x0[i+0]);
				 yv[0] = _mm512_load_ps(&y0[i+0]);
				 zv[0] = _mm512_fmadd_ps(xv[0],valpha,yv[0]);
				 _mm512_store_ps(&y0[i+0],zv[0]);
				 xv[1] = _mm512_load_ps(&x0[i+16]);
				 yv[1] = _mm512_load_ps(&y0[i+16]);
				 zv[1] = _mm512_fmadd_ps(xv[1],valpha,yv[1]);
				 _mm512_store_ps(&y0[i+16],zv[1]);
			  }
			 for(; (i+15) < n; i += 16) {
                               //_mm_prefetch((const char*)&x[i+0],_MM_HINT_T0);
				 xv[0] = _mm512_load_ps(&x0[i+0]);
				 yv[0] = _mm512_load_ps(&y0[i+0]);
				 zv[0] = _mm512_fmadd_ps(xv[0],valpha,yv[0]);
				 _mm512_store_ps(&y0[i+0],zv[0]);
 			 }
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
                  
                   void gms::math::saxpy_a_zmm16r4_unroll10x_omp(const int32_t n,
		                                      const float alpha,
		                                      float * __restrict __ATTR_ALIGN__(64) x,
                                                      const int32_t incx,
						      float * __restrict __ATTR_ALIGN__(64) y,
						      const int32_t incy) {
                         if(__builtin_expect(0==n,0) ||
			    __builtin_expect(0.0f == alpha,0)) { return;}
			    
			 __ATTR_ALIGN__(64) __m512 xv[10];
			 __ATTR_ALIGN__(64) __m512 yv[10];
			 __ATTR_ALIGN__(64) __m512 zv[10];
			 __m512 valpha;
			 float * __restrict x0 = NULL;
			 float * __restrict y0 = NULL;
			 int32_t i, last_i;
			 last_i = 0;
			 x0 = x;
			 y0 = y;

			 if(__builtin_expect(incx==1,1) &&
			    __builtin_expect(incy==1,1)) {

			    valpha = _mm512_broadcast_ss(alpha);
			     // Unrolled 10 times in order exploit the ratio
		             // between the FMA latency to its throughput.
#if defined(INTEL_COMPILER) || defined(__ICC)
                      __assume_aligned(x,64);
		      __assume_aligned(y,64);
		      __assume_aligned(x0,64);
		      __assume_aligned(y0,64);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC)
                      x  = (float*)__builtin_assume_aligned(x,64);
		      y  = (float*)__builtin_assume_aligned(y,64);
		      x0 = (float*)__builtin_assume_aligned(x0,64);
		      y0 = (float*)__builtin_assume_aligned(y0,64);
#endif
#pragma omp parallel for schedule(static,160) default(none)  \
                             lastprivate(last_i) private(i) shared(valpha,n,x0,y0)
			     for(i = 0; (i+159) < n; i += 160) {
			         last_i = i;
                                 _mm_prefetch((const char*)&x0[i+0],_MM_HINT_T0);
				 _mm512_store_ps(&y0[i+0],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[i+0]),valpha,
						                        _mm512_load_ps(&y0[i+0])));
				 _mm512_store_ps(&y0[i+16],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[i+16]),valpha,
						                        _mm512_load_ps(&y0[i+16])));
			         
				 _mm_prefetch((const char*)&x0[i+32],_MM_HINT_T0);
				 _mm512_store_ps(&y0[i+32],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[i+32]),valpha,
						                        _mm512_load_ps(&y0[i+32])));
				 _mm512_store_ps(&y0[i+48],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[i+48]),valpha,
						                        _mm512_load_ps(&y0[i+48])));
									
			      	 _mm_prefetch((const char*)&x0[i+64],_MM_HINT_T0);
				 _mm512_store_ps(&y0[i+64],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[i+64]),valpha,
						                        _mm512_load_ps(&y0[i+64])));
				 _mm512_store_ps(&y0[i+80],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[i+80]),valpha,
						                        _mm512_load_ps(&y0[i+80])));
				 
				 _mm_prefetch((const char*)&x0[i+96],_MM_HINT_T0);
				 _mm512_store_ps(&y0[i+96],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[i+96]),valpha,
						                        _mm512_load_ps(&y0[i+96])));
				 _mm512_store_ps(&y0[i+112],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[i+112]),valpha,
						                        _mm512_load_ps(&y0[i+112])));
				
				 _mm_prefetch((const char*)&x0[i+128],_MM_HINT_T0);
				 _mm512_store_ps(&y0[i+128],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[i+128]),valpha,
						                        _mm512_load_ps(&y0[i+128])));
				 _mm512_store_ps(&y0[i+144],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[i+144]),valpha,
						                        _mm512_load_ps(&y0[i+144])));
			
			     }
			    for(; (last_i+79) < n; last_i += 80) {
                                 _mm512_store_ps(&y0[i+0],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[last_i+0]),valpha,
						                        _mm512_load_ps(&y0[last_i+0])));
				 _mm512_store_ps(&y0[last_i+16],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[last_i+16]),valpha,
						                        _mm512_load_ps(&y0[last_i+16])));
			         
				// _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				 _mm512_store_ps(&y0[last_i+32],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[last_i+32]),valpha,
						                        _mm512_load_ps(&y0[last_i+32])));
				 _mm512_store_ps(&y0[last_i+48],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[last_i+48]),valpha,
						                        _mm512_load_ps(&y0[last_i+48])));
									
			      	 //_mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				 _mm512_store_ps(&y0[last_i+64],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[last_i+64]),valpha,
						                        _mm512_load_ps(&y0[last_i+64])));
			    }
			   for(; (last_i+63) < n; last_i += 64) {
                                 _mm512_store_ps(&y0[last_i+0],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[last_i+0]),valpha,
						                        _mm512_load_ps(&y0[last_i+0])));
				 _mm512_store_ps(&y0[last_i+16],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[last_i+16]),valpha,
						                        _mm512_load_ps(&y0[lasT_i+16])));
			         
				// _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				 _mm512_store_ps(&y0[last_i+32],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[last_i+32]),valpha,
						                        _mm512_load_ps(&y0[last_i+32])));
				 _mm512_store_ps(&y0[last_i+48],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[last_i+48]),valpha,
						                        _mm512_load_ps(&y0[last_i+48])));
			   }
			  for(; (last_i+31) < n; last_i += 32) {
                                 _mm512_store_ps(&y0[last_i+0],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[last_i+0]),valpha,
						                        _mm512_load_ps(&y0[last_i+0])));
				 _mm512_store_ps(&y0[last_i+16],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[last_i+16]),valpha,
						                        _mm512_load_ps(&y0[last_i+16])));
			  }
			 for(; (last_i+15) < n; last_i += 16) {
                               _mm512_store_ps(&y0[last_i+0],
				           _mm512_fmadd_ps(
					             _mm512_load_ps(&x0[last_i+0]),valpha,
						                        _mm512_load_ps(&y0[last_i+0])));
 			 }
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
	          
                   void gms::math::daxpy_u_zmm8r8_unroll10x(const int32_t n,
	                                         const double alpha,
	                                         double * __restrict  x,
					         const int32_t incx,
					         double * __restrict  y,
					         const int32_t incy) {
                   if(__builtin_except(0==n,0) ||
		      __builtin_except(alpha==0.0,0)) { return;}
		   
		   __ATTR_ALIGN__(64) __m512d xv[10];
		   __ATTR_ALIGN__(64) __m512d yv[10];
		   __ATTR_ALIGN__(64) __m512d zv[10];
                   __m512d valpha;
		   double * __restrict x0 = NULL;
		   double * __restrict y0 = NULL;
		   int32_t i;
		   x0 = x;
		   y0 = y;

		   if(__builtin_expect(incx==1,1) &&
		      __builtin_expect(incy==1,1) {

                      valpha = _mm512_broadcast_sd(alpha);
		      // Unrolled 10 times in order exploit the ratio
		      // between the FMA latency to its throughput.
		      for(i = 0; (i+79) < n; i += 80) {
		          _mm_prefetch((const char*)&x[i+0],_MM_HINT_T0);
                          xv[0] = _mm512_loadu_pd(&x0[i+0]);
			  yv[0] = _mm512_loadu_pd(&y0[i+0]);
			  zv[0] = _mm512_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm512_storeu_pd(&y0[i+0],zv[0]);
			  xv[1] = _mm512_loadu_pd(&x0[i+8]);
			  yv[1] = _mm512_loadu_pd(&y0[i+8]);
			  zv[1] = _mm512_fmadd_pd(xv[1],valpha,yv[1]);
			  _mm512_storeu_pd(&y0[i+8],zv[1]);
			  _mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
			  xv[2] = _mm512_loadu_pd(&x0[i+16]);
			  yv[2] = _mm512_loadu_pd(&y0[i+16]);
			  zv[2] = _mm512_fmadd_pd(xv[2],valpha,yv[2]);
			  _mm512_storeu_pd(&y0[i+16],zv[2]);
			  xv[3] = _mm512_loadu_pd(&x0[i+24]);
			  yv[3] = _mm512_loadu_pd(&y0[i+24]);
			  zv[3] = _mm512_fmadd_pd(xv[3],valpha,yv[3]);
			  _mm512_storeu_pd(&y0[i+24],zv[3]);
			  _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
			  xv[4] = _mm512_loadu_pd(&x0[i+32]);
			  yv[4] = _mm512_loadu_pd(&y0[i+32]);
			  zv[4] = _mm512_fmadd_pd(xv[4],valpha,yv[4]);
			  _mm512_storeu_pd(&y0[i+32],zv[4]);
			  xv[5] = _mm512_loadu_pd(&x0[i+40]);
			  yv[5] = _mm512_loadu_pd(&y0[i+40]);
			  zv[5] = _mm512_fmadd_pd(xv[5],valpha,yv[5]);
			  _mm512_storeu_pd(&y0[i+40],zv[5]);
			  _mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
			  xv[6] = _mm512_loadu_pd(&x0[i+48]);
			  yv[6] = _mm512_loadu_pd(&y0[i+48]);
			  zv[6] = _mm512_fmadd_pd(xv[6],valpha,yv[6]);
			  _mm512_storeu_pd(&y0[i+48],zv[6]);
			  xv[7] = _mm512_loadu_pd(&x0[i+56]);
			  yv[7] = _mm512_loadu_pd(&y0[i+56]);
			  zv[7] = _mm512_fmadd_pd(xv[7],valpha,yv[7]);
			  _mm512_storeu_pd(&y0[i+56],zv[7]);
			   _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
			  xv[8] = _mm512_loadu_pd(&x0[i+64]);
			  yv[8] = _mm512_loadu_pd(&y0[i+64]);
			  zv[8] = _mm512_fmadd_pd(xv[8],valpha,yv[8]);
			  _mm256_storeu_pd(&y0[i+64],zv[8]);
			  xv[9] = _mm256_loadu_ps(&x0[i+72]);
			  yv[9] = _mm256_loadu_ps(&y0[i+72]);
			  zv[9] = _mm256_fmadd_ps(xv[9],valpha,yv[9]);
			  _mm256_storeu_ps(&y0[i+72],zv[9]);
		      }
		      for(; (i+39) < n; i += 40) {

		          //_mm_prefetch((const char*)&x[i+0],_MM_HINT_T0);
                          xv[0] = _mm512_loadu_pd(&x0[i+0]);
			  yv[0] = _mm512_loadu_pd(&y0[i+0]);
			  zv[0] = _mm512_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm512_storeu_pd(&y0[i+0],zv[0]);
			  xv[1] = _mm512_loadu_pd(&x0[i+8]);
			  yv[1] = _mm512_loadu_pd(&y0[i+8]);
			  zv[1] = _mm512_fmadd_pd(xv[1],valpha,yv[1]);
			  _mm512_storeu_pd(&y0[i+8],zv[1]);
			  //_mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
			  xv[2] = _mm512_loadu_pd(&x0[i+16]);
			  yv[2] = _mm512_loadu_pd(&y0[i+16]);
			  zv[2] = _mm512_fmadd_pd(xv[2],valpha,yv[2]);
			  _mm512_storeu_pd(&y0[i+16],zv[2]);
			  xv[3] = _mm512_loadu_pd(&x0[i+24]);
			  yv[3] = _mm512_loadu_pd(&y0[i+24]);
			  zv[3] = _mm512_fmadd_pd(xv[3],valpha,yv[3]);
			  _mm512_storeu_pd(&y0[i+24],zv[3]);
			 // _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
			  xv[4] = _mm512_loadu_pd(&x0[i+32]);
			  yv[4] = _mm512_loadu_pd(&y0[i+32]);
			  zv[4] = _mm512_fmadd_pd(xv[4],valpha,yv[4]);
			  _mm512_storeu_pd(&y0[i+32],zv[4]);
		         
		      }
		      for(; (i+31) < n; i += 32) {
                         
		           //_mm_prefetch((const char*)&x[i+0],_MM_HINT_T0);
                          xv[0] = _mm512_loadu_pd(&x0[i+0]);
			  yv[0] = _mm512_loadu_pd(&y0[i+0]);
			  zv[0] = _mm512_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm512_storeu_pd(&y0[i+0],zv[0]);
			  xv[1] = _mm512_loadu_pd(&x0[i+8]);
			  yv[1] = _mm512_loadu_pd(&y0[i+8]);
			  zv[1] = _mm512_fmadd_pd(xv[1],valpha,yv[1]);
			  _mm512_storeu_pd(&y0[i+8],zv[1]);
			  //_mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
			  xv[2] = _mm512_loadu_pd(&x0[i+16]);
			  yv[2] = _mm512_loadu_pd(&y0[i+16]);
			  zv[2] = _mm512_fmadd_pd(xv[2],valpha,yv[2]);
			  _mm512_storeu_pd(&y0[i+16],zv[2]);
			  xv[3] = _mm512_loadu_pd(&x0[i+24]);
			  yv[3] = _mm512_loadu_pd(&y0[i+24]);
			  zv[3] = _mm512_fmadd_pd(xv[3],valpha,yv[3]);
			  _mm512_storeu_pd(&y0[i+24],zv[3]);
		      }
		      for(; (i+15) < n; i += 16) {

		           //_mm_prefetch((const char*)&x[i+0],_MM_HINT_T0);
                          xv[0] = _mm512_loadu_pd(&x0[i+0]);
			  yv[0] = _mm512_loadu_pd(&y0[i+0]);
			  zv[0] = _mm512_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm512_storeu_pd(&y0[i+0],zv[0]);
			  xv[1] = _mm512_loadu_pd(&x0[i+8]);
			  yv[1] = _mm512_loadu_pd(&y0[i+8]);
			  zv[1] = _mm512_fmadd_pd(xv[1],valpha,yv[1]);
			  _mm512_storeu_pd(&y0[i+8],zv[1]);
		      }
		      for(; (i+7) < n; i += 8) {

		          xv[0] = _mm512_loadu_pd(&x0[i+0]);
			  yv[0] = _mm512_loadu_pd(&y0[i+0]);
			  zv[0] = _mm512_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm512_storeu_pd(&y0[i+0],zv[0]);
		      }
		      
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


	         
                   void gms::math::daxpy_a_zmm8r8_unroll10x(const int32_t n,
	                                         const double alpha,
	                                         double * __restrict __ATTR_ALIGN__(64) x,
					         const int32_t incx,
					         double * __restrict  __ATTR_ALIGN__(64) y,
					         const int32_t incy) {
                   if(__builtin_except(0==n,0) ||
		      __builtin_except(alpha==0.0,0)) { return;}
		   
		   __ATTR_ALIGN__(64) __m512d xv[10];
		   __ATTR_ALIGN__(64) __m512d yv[10];
		   __ATTR_ALIGN__(64) __m512d zv[10];
                   __m512d valpha;
		   double * __restrict x0 = NULL;
		   double * __restrict y0 = NULL;
		   int32_t i;
		   x0 = x;
		   y0 = y;

		   if(__builtin_expect(incx==1,1) &&
		      __builtin_expect(incy==1,1) {

                      valpha = _mm512_broadcast_sd(alpha);
		      // Unrolled 10 times in order exploit the ratio
		      // between the FMA latency to its throughput.
#if defined(INTEL_COMPILER) || defined(__ICC)
                      __assume_aligned(x,64);
		      __assume_aligned(y,64);
		      __assume_aligned(x0,64);
		      __assume_aligned(y0,64);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC)
                      x  = (double*)__builtin_assume_aligned(x,64);
		      y  = (double*)__builtin_assume_aligned(y,64);
		      x0 = (double*)__builtin_assume_aligned(x0,64);
		      y0 = (double*)__builtin_assume_aligned(y0,64);
#endif
		      for(i = 0; (i+79) < n; i += 80) {
		          _mm_prefetch((const char*)&x0[i+0],_MM_HINT_T0);
                          xv[0] = _mm512_load_pd(&x0[i+0]);
			  yv[0] = _mm512_load_pd(&y0[i+0]);
			  zv[0] = _mm512_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm512_store_pd(&y0[i+0],zv[0]);
			  xv[1] = _mm512_load_pd(&x0[i+8]);
			  yv[1] = _mm512_load_pd(&y0[i+8]);
			  zv[1] = _mm512_fmadd_pd(xv[1],valpha,yv[1]);
			  _mm512_store_pd(&y0[i+8],zv[1]);
			  _mm_prefetch((const char*)&x0[i+16],_MM_HINT_T0);
			  xv[2] = _mm512_load_pd(&x0[i+16]);
			  yv[2] = _mm512_load_pd(&y0[i+16]);
			  zv[2] = _mm512_fmadd_pd(xv[2],valpha,yv[2]);
			  _mm512_store_pd(&y0[i+16],zv[2]);
			  xv[3] = _mm512_load_pd(&x0[i+24]);
			  yv[3] = _mm512_load_pd(&y0[i+24]);
			  zv[3] = _mm512_fmadd_pd(xv[3],valpha,yv[3]);
			  _mm512_store_pd(&y0[i+24],zv[3]);
			  _mm_prefetch((const char*)&x0[i+32],_MM_HINT_T0);
			  xv[4] = _mm512_load_pd(&x0[i+32]);
			  yv[4] = _mm512_load_pd(&y0[i+32]);
			  zv[4] = _mm512_fmadd_pd(xv[4],valpha,yv[4]);
			  _mm512_store_pd(&y0[i+32],zv[4]);
			  xv[5] = _mm512_load_pd(&x0[i+40]);
			  yv[5] = _mm512_load_pd(&y0[i+40]);
			  zv[5] = _mm512_fmadd_pd(xv[5],valpha,yv[5]);
			  _mm512_store_pd(&y0[i+40],zv[5]);
			  _mm_prefetch((const char*)&x0[i+48],_MM_HINT_T0);
			  xv[6] = _mm512_load_pd(&x0[i+48]);
			  yv[6] = _mm512_load_pd(&y0[i+48]);
			  zv[6] = _mm512_fmadd_pd(xv[6],valpha,yv[6]);
			  _mm512_storeu_pd(&y0[i+48],zv[6]);
			  xv[7] = _mm512_load_pd(&x0[i+56]);
			  yv[7] = _mm512_load_pd(&y0[i+56]);
			  zv[7] = _mm512_fmadd_pd(xv[7],valpha,yv[7]);
			  _mm512_store_pd(&y0[i+56],zv[7]);
			   _mm_prefetch((const char*)&x0[i+64],_MM_HINT_T0);
			  xv[8] = _mm512_load_pd(&x0[i+64]);
			  yv[8] = _mm512_load_pd(&y0[i+64]);
			  zv[8] = _mm512_fmadd_pd(xv[8],valpha,yv[8]);
			  _mm256_store_pd(&y0[i+64],zv[8]);
			  xv[9] = _mm256_load_ps(&x0[i+72]);
			  yv[9] = _mm256_load_ps(&y0[i+72]);
			  zv[9] = _mm256_fmadd_ps(xv[9],valpha,yv[9]);
			  _mm256_store_ps(&y0[i+72],zv[9]);
		      }
		      for(; (i+39) < n; i += 40) {

		          //_mm_prefetch((const char*)&x[i+0],_MM_HINT_T0);
                          xv[0] = _mm512_load_pd(&x0[i+0]);
			  yv[0] = _mm512_load_pd(&y0[i+0]);
			  zv[0] = _mm512_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm512_store_pd(&y0[i+0],zv[0]);
			  xv[1] = _mm512_load_pd(&x0[i+8]);
			  yv[1] = _mm512_load_pd(&y0[i+8]);
			  zv[1] = _mm512_fmadd_pd(xv[1],valpha,yv[1]);
			  _mm512_store_pd(&y0[i+8],zv[1]);
			  //_mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
			  xv[2] = _mm512_load_pd(&x0[i+16]);
			  yv[2] = _mm512_load_pd(&y0[i+16]);
			  zv[2] = _mm512_fmadd_pd(xv[2],valpha,yv[2]);
			  _mm512_store_pd(&y0[i+16],zv[2]);
			  xv[3] = _mm512_load_pd(&x0[i+24]);
			  yv[3] = _mm512_load_pd(&y0[i+24]);
			  zv[3] = _mm512_fmadd_pd(xv[3],valpha,yv[3]);
			  _mm512_store_pd(&y0[i+24],zv[3]);
			 // _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
			  xv[4] = _mm512_load_pd(&x0[i+32]);
			  yv[4] = _mm512_load_pd(&y0[i+32]);
			  zv[4] = _mm512_fmadd_pd(xv[4],valpha,yv[4]);
			  _mm512_store_pd(&y0[i+32],zv[4]);
		         
		      }
		      for(; (i+31) < n; i += 32) {
                         
		           //_mm_prefetch((const char*)&x[i+0],_MM_HINT_T0);
                          xv[0] = _mm512_load_pd(&x0[i+0]);
			  yv[0] = _mm512_load_pd(&y0[i+0]);
			  zv[0] = _mm512_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm512_storeu_pd(&y0[i+0],zv[0]);
			  xv[1] = _mm512_load_pd(&x0[i+8]);
			  yv[1] = _mm512_load_pd(&y0[i+8]);
			  zv[1] = _mm512_fmadd_pd(xv[1],valpha,yv[1]);
			  _mm512_store_pd(&y0[i+8],zv[1]);
			  //_mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
			  xv[2] = _mm512_load_pd(&x0[i+16]);
			  yv[2] = _mm512_load_pd(&y0[i+16]);
			  zv[2] = _mm512_fmadd_pd(xv[2],valpha,yv[2]);
			  _mm512_store_pd(&y0[i+16],zv[2]);
			  xv[3] = _mm512_load_pd(&x0[i+24]);
			  yv[3] = _mm512_load_pd(&y0[i+24]);
			  zv[3] = _mm512_fmadd_pd(xv[3],valpha,yv[3]);
			  _mm512_store_pd(&y0[i+24],zv[3]);
		      }
		      for(; (i+15) < n; i += 16) {

		           //_mm_prefetch((const char*)&x[i+0],_MM_HINT_T0);
                          xv[0] = _mm512_load_pd(&x0[i+0]);
			  yv[0] = _mm512_load_pd(&y0[i+0]);
			  zv[0] = _mm512_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm512_store_pd(&y0[i+0],zv[0]);
			  xv[1] = _mm512_load_pd(&x0[i+8]);
			  yv[1] = _mm512_load_pd(&y0[i+8]);
			  zv[1] = _mm512_fmadd_pd(xv[1],valpha,yv[1]);
			  _mm512_store_pd(&y0[i+8],zv[1]);
		      }
		      for(; (i+7) < n; i += 8) {

		          xv[0] = _mm512_load_pd(&x0[i+0]);
			  yv[0] = _mm512_load_pd(&y0[i+0]);
			  zv[0] = _mm512_fmadd_pd(xv[0],valpha,yv[0]);
			  _mm512_store_pd(&y0[i+0],zv[0]);
		      }
		      
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



		  
                   void gms::math::daxpy_a_zmm8r8_unroll10x_omp(const int32_t n,
	                                             const double alpha,
	                                             double * __restrict __ATTR_ALIGN__(64) x,
					             const int32_t incx,
					             double * __restrict  __ATTR_ALIGN__(64) y,
					             const int32_t incy) {
                   if(__builtin_except(0==n,0) ||
		      __builtin_except(alpha==0.0,0)) { return;}
		   
		   __ATTR_ALIGN__(64) __m512d xv[10];
		   __ATTR_ALIGN__(64) __m512d yv[10];
		   __ATTR_ALIGN__(64) __m512d zv[10];
                   __m512d valpha;
		   double * __restrict x0 = NULL;
		   double * __restrict y0 = NULL;
		   int32_t i,last_i;
		   last_i = 0;
		   x0 = x;
		   y0 = y;

		   if(__builtin_expect(incx==1,1) &&
		      __builtin_expect(incy==1,1) {

                      valpha = _mm512_broadcast_sd(alpha);
		      // Unrolled 10 times in order exploit the ratio
		      // between the FMA latency to its throughput.
#if defined(INTEL_COMPILER) || defined(__ICC)
                      __assume_aligned(x,64);
		      __assume_aligned(y,64);
		      __assume_aligned(x0,64);
		      __assume_aligned(y0,64);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC)
                      x  = (double*)__builtin_assume_aligned(x,64);
		      y  = (double*)__builtin_assume_aligned(y,64);
		      x0 = (double*)__builtin_assume_aligned(x0,64);
		      y0 = (double*)__builtin_assume_aligned(y0,64);
#endif
#pragma omp parallel for schedule(static,160) default(none)  \
                      lastprivate(last_i) private(i) shared(valpha,n,x0,y0)
		      for(i = 0; (i+79) < n; i += 80) {
		          _mm_prefetch((const char*)&x0[i+0],_MM_HINT_T0);
			  _mm512_store_pd(&y0[i+0],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+0]),valpha,
					                  _mm512_load_pd(&y0[i+0])));
			  _mm512_store_pd(&y0[i+8],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+8]),valpha,
					                  _mm512_load_pd(&y0[i+8])));
					            
                          _mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
			  _mm512_store_pd(&y0[i+16],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+16]),valpha,
					                  _mm512_load_pd(&y0[i+16])));
			  _mm512_store_pd(&y0[i+24],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+24]),valpha,
					                  _mm512_load_pd(&y0[i+24])));
			  _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
			  _mm512_store_pd(&y[i+32],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+32]),valpha,
					                  _mm512_load_pd(&y0[i+32])));
			  _mm512_store_pd(&y[i+40],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+40]),valpha,
					                  _mm512_load_pd(&y0[i+40])));
			  _mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
			  _mm512_store_pd(&y0[i+48],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+48]),valpha,
					                  _mm512_load_pd(&y0[i+48])));
			  _mm512_store_pd(&y0[i+56],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+56]),valpha,
					                  _mm512_load_pd(&y0[i+56])));
			  _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
			  _mm512_store_pd(&y[i+64],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+64]),valpha,
					                  _mm512_load_pd(&y0[i+64])));
			  _mm512_store_pd(&y0[i+72],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+72]),valpha,
					                  _mm512_load_pd(&y0[i+72])));
			
		      }
		      for(; (i+39) < n; i += 40) {

		           _mm512_store_pd(&y0[i+0],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+0]),valpha,
					                  _mm512_load_pd(&y0[i+0])));
			  _mm512_store_pd(&y0[i+8],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+8]),valpha,
					                  _mm512_load_pd(&y0[i+8])));
					            
                          //_mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
			  _mm512_store_pd(&y0[i+16],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+16]),valpha,
					                  _mm512_load_pd(&y0[i+16])));
			  _mm512_store_pd(&y0[i+24],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+24]),valpha,
					                  _mm512_load_pd(&y0[i+24])));
			  //_mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
			  _mm512_store_pd(&y0[i+32],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+32]),valpha,
					                  _mm512_load_pd(&y0[i+32])));
		         
		      }
		      for(; (i+31) < n; i += 32) {
                         
		          _mm512_store_pd(&y0[i+0],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+0]),valpha,
					                  _mm512_load_pd(&y0[i+0])));
			  _mm512_store_pd(&y0[i+8],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+8]),valpha,
					                  _mm512_load_pd(&y0[i+8])));
					            
                          //_mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
			  _mm512_store_pd(&y[i+16],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+16]),valpha,
					                  _mm512_load_pd(&y0[i+16])));
			  _mm512_store_pd(&y0[i+24],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+24]),valpha,
					                  _mm512_load_pd(&y0[i+24])));
		      }
		      for(; (i+15) < n; i += 16) {

		           _mm512_store_pd(&y0[i+0],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+0]),valpha,
					                  _mm512_load_pd(&y0[i+0])));
			  _mm512_store_pd(&y0[i+8],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+8]),valpha,
					                  _mm512_load_pd(&y0[i+8])));
		      }
		      for(; (i+7) < n; i += 8) {

		         _mm512_store_pd(&y0[i+0],
				   _mm512_fmadd_pd(
				           _mm512_load_pd(&x0[i+0]),valpha,
					                  _mm512_load_pd(&y0[i+0])));
		      }
		      
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




	        


		   

  
