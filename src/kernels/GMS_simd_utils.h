

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

#ifndef __GMS_SIMD_UTILS_H__
#define __GMS_SIMD_UTILS_H__ 040120220918




namespace file_info 
{

const unsigned int GMS_SIMD_UTILS_MAJOR = 1U;
const unsigned int GMS_SIMD_UTILS_MINOR = 0U;
const unsigned int GMS_SIMD_UTILS_MICRO = 1U;
const unsigned int GMS_SIMD_UTILS_FULLVER =
       1000U*GMS_SIMD_UTILS_MAJOR+
       100U*GMS_SIMD_UTILS_MINOR +
       10U*GMS_SIMD_UTILS_MICRO;
const char * const GMS_SIMD_UTILS_CREATION_DATE = "04-01-2022 09:18 AM +00200 (TUE 04 JAN 2022 GMT+2)";
const char * const GMS_SIMD_UTILS_BUILD_DATE    = __DATE__ ":" __TIME__;
const char * const GMS_SIMD_UTILS_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
const char * const GMS_SIMD_UTILS_DESCRIPTION   = "Various SIMD utility functions.";


}



#include <immintrin.h>
#include "GMS_config.h"


namespace  gms {


          namespace  math {


            namespace 
			{

                          const __m128  _0PS     = _mm_set1_ps(0.0F);
			              const __m256d _0PD     = _mm256_setzero_pd();
			              const __m128   NZ128SP = _mm_set1_ps(-0.0F);
						  const __m128d  NZ128DP = _mm_set1_pd(-0.0);
			              const __m256   NZ256SP = _mm256_set1_ps(-0.0F);
			              const __m256d  NZ256DP = _mm256_set1_pd(-0.0);
			              const __m512   NZ512SP = _mm512_set1_ps(-0.0F);
			              const __m512d  NZ512DP = _mm512_set1_pd(-0.0F);
		      }


	              
        __ATTR_ALWAYS_INLINE__
		static inline
		__mmask8 isinf_zmm8r8(__m512d x) {
             union {
               __m512i u;
			   __m512d f;
			 } ieee754;
			 const __m512i c0 = _mm512_set1_epi64(0x7fffffff);
			 const __m512i c1 = _mm512_set1_epi64(0x7ff00000);
			 const __m512i _0 = _mm512_setzero_si512();
			 __m512i t0,t1;
			 __mmask8 b0,b1;
			 ieee754.f = x;
			 t0 = _mm512_and_epi64(_mm512_srli_epi64(ieee754.u,32),c0);
			 b0 = _mm512_cmp_epi64_mask(t0,c1,_MM_CMPINT_EQ);
			 b1 = _mm512_cmp_epi64_mask(ieee754.u,_0,_MM_CMPINT_EQ);
			 return (b0 && b1);
		}

/*
		     
                      __ATTR_ALWAYS_INLINE__
		      static inline
		      __mmask8 isnan_zmm8r8(__m512d x) {

		        union {
                           __m512i u;
			   __m512d f;
			 } ieee754;
			 const __m512i c0 = _mm512_set1_epi64(0x7fffffff);
			 const __m512i c1 = _mm512_set1_epi64(0x7ff00000);
			 const __m512i _0 = _mm512_setzero_si512();
			 __m512i t0,t1;
			 __mmask8 b0,b1;
			 ieee754.f = x;
			 t0 = _mm512_and_epi64(_mm512_srli_epi64(ieee754.u,32),c0);
			 b0 = _mm512_cmp_epi64_mask(ieee754.u,_0,_MM_CMPINT_NE);
			 t1 = _mm512_movm_epi64(b0);
			 b1 = _mm512_cmp_epi64_mask(t1,c1,_MM_CMPINT_LT);
			 b0 = _mm512_movm_epi64(b1);
			 return (b0+b1);
			
		      }
*/

		      // Load only 3 elements (lower) of XMM register.
		      // Single-precision

                     
            __ATTR_ALWAYS_INLINE__
		    static inline
		    __m128
		    xmm4r4_load_3u_avx512(const float * __restrict v) {
                            const __mmask8 k = 0x7;
                            return (_mm_mask_loadu_ps(_0PS,k,v));
		    }


		      
            __ATTR_ALWAYS_INLINE__
		    static inline
		    __m128
		    xmm4r4_load_3a_avx512(const float * __restrict __ATTR_ALIGN__(16) v) {
                            const __mmask8 k = 0x7;
                            return (_mm_mask_load_ps(_0PS,k,v));
		    }


		     
            __ATTR_ALWAYS_INLINE__
		    static inline
		    __m128
		    xmm4r4_load_3u_avx(const float * __restrict v) {

                          const __m128i k = _mm_set_epi32(0,-1,-1,-1);
			  return (_mm_maskload_ps(v,(__m128i)k));
		    }


		    
            __ATTR_ALWAYS_INLINE__
		    static inline
		    __m128
		    xmm4r4_load_3a_avx(const float * __restrict __ATTR_ALIGN__(16) v) {

                          const __m128i k = _mm_set_epi32(0,-1,-1,-1);
			  return (_mm_maskload_ps(v,(__m128i)k));
		    }

                    // Load only 3 elements (lower) of XMM register.
		    // Double-precision
		  
            __ATTR_ALWAYS_INLINE__
		    static inline
		    __m256d
		    ymm8r4_load_3u_avx512(const double * __restrict v) {

		          const __mmask8 k = 0x7;
			  return (_mm256_mask_loadu_pd(_0PD,k,v));
		    }


		    
            __ATTR_ALWAYS_INLINE__
		    static inline
		    __m256d
		    ymm8r4_load_3a_avx512(const double * __restrict v) {

		          const __mmask8 k = 0x7;
			  return (_mm256_mask_load_pd(_0PD,k,v));
		    }


/*	    
                      __ATTR_ALWAYS_INLINE__
		     static inline
		      __m256d
		      ymm8r4_load_3u_avx(const double * __restrict v) {

		          const __m256i k = _mm256_set_epi32(0,-1,-1,-1);
			  return (_mm256_maskload_pd(v,k));
		    }


		    
                      __ATTR_ALWAYS_INLINE__
		     
		     
		      static inline
		      __m256d
		      ymm8r4_load_3a_avx(const double * __restrict v) {

		          const __m256i k = _mm256_set_epi32(0,-1,-1,-1);
			  return (_mm256_maskload_pd(v,k));
		    }
*/		    
		  // Store only 3 elements (lower) of XMM register.

		    
            __ATTR_ALWAYS_INLINE__
		    static inline   
            void
		    xmm4r4_store_3u_avx512(float * __restrict v,
		                             const __m128 x) {

                          const __mmask8 k = 0x7;
			  return (_mm_mask_storeu_ps(v,k,x));
		   }


		    
            __ATTR_ALWAYS_INLINE__		     
		    static inline   
            void
		    xmm4r4_store_3a_avx512(float * __restrict  __ATTR_ALIGN__(16) v,
		                             const __m128 x) {

                          const __mmask8 k = 0x7;
			  return (_mm_mask_store_ps(v,k,x));
		   }


		    
             __ATTR_ALWAYS_INLINE__		    		     
		      static inline   
                      void
		      xmm4r4_store_3u_avx(float * __restrict v,
		                          const __m128 x) 
		  {

                          const __m128i k = _mm_set_epi32(0,-1,-1,-1);
			  return (_mm_maskstore_ps(v,(__m128i)k,x));
		   }


		    
            __ATTR_ALWAYS_INLINE__		     
		      static inline   
                      void
		      xmm4r4_store_3a_avx(float * __restrict __ATTR_ALIGN__(16) v,
		                          const __m128 x) {

                          const __m128i k = _mm_set_epi32(0,-1,-1,-1);
			  return (_mm_maskstore_ps(v,(__m128i)k,x));
		   }


		    
            __ATTR_ALWAYS_INLINE__		     
		      static inline   
                      void
		      ymm4r8_store_3u_avx512(double * __restrict v,
		                          const __m256d x) {

                          const __mmask8 k = 0x7;
			  return (_mm256_mask_storeu_pd(v,k,x));
		   }


		    
            __ATTR_ALWAYS_INLINE__		     
		      static inline   
                      void
		      ymm8r4_store_3a_avx512(double * __restrict __ATTR_ALIGN__(16) v,
		                             const __m256d x) {

                          const __mmask8 k = 0x7;
			  return (_mm256_mask_store_pd(v,k,x));
		   }


		   // The whole register negated
		   
		     
              __ATTR_ALWAYS_INLINE__		    
		      static inline 
             __m128
		      negate_xmm4r4(const __m128 v) 
			  {
                   const __m128  NZ128SP = _mm_set1_ps(-0.0F);
		           return (_mm_xor_ps(v,NZ128SP));
		      }
		   
                   
		    
            __ATTR_ALWAYS_INLINE__		     
		      static inline 
                      __m256
		      negate_ymm8r4(const __m256 v) {

		           return (_mm256_xor_ps(v,NZ256SP));
		    }
		    
		    
                    
            __ATTR_ALWAYS_INLINE__		     
		      static inline 
                      __m128d
		      negate_xmm2r8(const __m128d v) {

		           return (_mm_xor_pd(v,NZ128DP));
		    }
		   

		    
            __ATTR_ALWAYS_INLINE__		     
		      static inline 
                      __m256d
		      negate_ymm4r8(const __m256d v) {

		           return (_mm256_xor_pd(v,NZ256DP));
		    }


		    
            __ATTR_ALWAYS_INLINE__		     
		      static inline 
                      __m512
		      negate_zmm16r4(const __m512 v) {

                          return (_mm512_xor_ps(v,NZ512SP));
		    }


		    
            __ATTR_ALWAYS_INLINE__		     
		      static inline 
                      __m512d
		      negate_zmm8r8(const __m512d v) {

                          return (_mm512_xor_pd(v,NZ512DP));
		    }


		    // Dot product
		    
            __ATTR_ALWAYS_INLINE__		     
		      static inline 
                      __m128d
		      ymm8r4_dot(const __m256d x,
		                 const __m256d y) {

			  const __m256d t   = _mm256_mul_pd(x,y);
                          const __m256d slo = _mm256_hadd_pd(t, t);
                          const __m128d shi = _mm256_extractf128_pd(slo, 0x1);
			  return _mm_add_pd(shi, _mm256_castpd256_pd128(slo));
		   }
		    
                   // Getting a single value form SIMD register

		    
              __ATTR_ALWAYS_INLINE__		     
		      static inline 
		      float
		      ymm8r4_0_elem(__m256 a) {
                           return _mm_cvtss_f32(_mm256_castps256_ps128(a));
                      }
		      

		    
            __ATTR_ALWAYS_INLINE__		     
		      static inline 
                      float
		      ymm8r4_1_elem(__m256 a) {
                           __m128 t = _mm256_castps256_ps128(a);
                           return _mm_cvtss_f32(_mm_shuffle_ps(t,t,_MM_SHUFFLE(0,0,0,1)));
                      }

		    
            __ATTR_ALWAYS_INLINE__		     
		      static inline 
                     float
		     ymm8r4_2_elem(__m256 a) {
                           __m128 t = _mm256_castps256_ps128(a);
                           return _mm_cvtss_f32(_mm_shuffle_ps(t,t,_MM_SHUFFLE(0,0,0,2)));
                     }

		    
            __ATTR_ALWAYS_INLINE__		     
		      static inline 
                      float
		      ymm8r4_3_elem(__m256 a) {
                           __m128 t = _mm256_castps256_ps128(a);
                           return _mm_cvtss_f32(_mm_shuffle_ps(t,t,_MM_SHUFFLE(0,0,0,3)));
                     }

		    
            __ATTR_ALWAYS_INLINE__		     
		      static inline 
                      float
		      ymm8r4_4_elem(__m256 a) {
                         return _mm_cvtss_f32(_mm256_extractf128_ps(a,0x1));
                      }

		     
            __ATTR_ALWAYS_INLINE__
		      static inline 
                      float
		      ymm8r4_5_elem(__m256 a) {
                           __m128 t = _mm256_extractf128_ps(a,0x1);
                           return _mm_cvtss_f32(_mm_shuffle_ps(t,t,_MM_SHUFFLE(0,0,0,1)));
                     }
		     
		     
            __ATTR_ALWAYS_INLINE__
		       static inline 
                      float
		      ymm8r4_6_elem(__m256 a) {
                            __m128 t = _mm256_extractf128_ps(a,0x1);
                            return _mm_cvtss_f32(_mm_shuffle_ps(t,t,_MM_SHUFFLE(0,0,0,2)));
                      }

		    
            __ATTR_ALWAYS_INLINE__
		      static inline 
                      float
		      ymm8r4_7_elem(__m256 a) {
                            __m128 t = _mm256_extractf128_ps(a,0x1);
                            return _mm_cvtss_f32(_mm_shuffle_ps(t,t,_MM_SHUFFLE(0,0,0,3)));
                      }

		   
            __ATTR_ALWAYS_INLINE__
		      static inline 
                      double
		      ymm4r8_0_elem(__m256d a) {
                          return _mm_cvtsd_f64(_mm256_castpd256_pd128(a));
                      }

		     
            __ATTR_ALWAYS_INLINE__
		     static inline 
                      double
		      ymm4r8_1_elem(__m256d a) {
                            __m128d t = _mm256_castpd256_pd128(a);
                           return _mm_cvtsd_f64(_mm_shuffle_pd(t,t,_MM_SHUFFLE2(0,1)));
                      }

		     
            __ATTR_ALWAYS_INLINE__
		       static inline 
                      double
		      ymm4r8_2_elem(__m256d a) {
                            return _mm_cvtsd_f64(_mm256_extractf128_pd(a,0x1));
                      }

		    
            __ATTR_ALWAYS_INLINE__
		     static inline 
                      double
		      ymm4r8_3_elem(__m256d a) {
                           __m128d t = _mm256_extractf128_pd(a,0x1);
                           return _mm_cvtsd_f64(_mm_shuffle_pd(t,t,_MM_SHUFFLE2(0,1)));
                     }


		     // Horizontal summation.
                   
            __ATTR_ALWAYS_INLINE__
		      static inline 
		      float
		      ymm8r4_horizontal_sum(const __m256 a) {
                          __m256 sum    = _mm256_hadd_ps(a, a);
                          sum           = _mm256_hadd_ps(sum, sum);
                          __m128 r      = _mm_add_ps(_mm256_castps256_ps128(sum),
			                          _mm256_extractf128_ps(sum, 0x1));
                          return _mm_cvtss_f32(r);
                      }

		      
            __ATTR_ALWAYS_INLINE__
		      static inline 
                      double
		      ymm4r8_horizontal_sum(const __m256d a) {
                            __m256d sum     = _mm256_hadd_pd(a, a);
                            const __m128d r = _mm_add_sd(_mm256_castpd256_pd128(sum),
			                                 _mm256_extractf128_pd(sum, 0x1));
                             return _mm_cvtsd_f64(r);
                      }

		  
            __ATTR_ALWAYS_INLINE__
		       static inline 
                      double
		      ymm4r8_horizontal_prod(const __m256d a) {
                           const __m256d sum = _mm256_mul_pd(a, _mm256_shuffle_pd(a,a,0x5));
                           const __m128d sh = _mm256_extractf128_pd(sum, 0x1);
                           __m128d r        = _mm_mul_sd(sh, _mm256_castpd256_pd128(sum));
                           return _mm_cvtsd_f64(r);
                      }

			
             __ATTR_ALWAYS_INLINE__
		     static inline 
		      __m128
		      xmm4r4_abs(const __m128 x) 
			  {
                          const __m128 mask = _mm_set1_ps(-0.0f); 
                          return _mm_andnot_ps(mask, x);
              }
		     

		     
              __ATTR_ALWAYS_INLINE__
		       static inline 
             __m128d
		     xmm2r8_abs(const __m128d x) 
			 {
                           const __m128d mask = _mm_set1_pd(-0.0); 
                           return _mm_andnot_pd(mask, x); 
             }
		      
                   
             __ATTR_ALWAYS_INLINE__
		     static inline 
		      __m256
		      ymm8r4_abs(const __m256 x) 
			  {
                          const __m256 mask = _mm256_set1_ps(-0.0f); 
                          return _mm256_andnot_ps(mask, x);
              }
		     

		     
              __ATTR_ALWAYS_INLINE__
		       static inline 
             __m256d
		     ymm4r8_abs(const __m256d x) 
			 {
                           const __m256d mask = _mm256_set1_pd(-0.0); 
                           return _mm256_andnot_pd(mask, x); 
             }
                     

             __ATTR_ALWAYS_INLINE__
             static inline 
			 __m128i 
             xmm4i4_maskload(int const * __restrict mem, __m128i mask) 
			 {
                               return _mm_castps_si128(_mm_maskload_ps((float const *) mem, mask));
             }

                 




		     
    } // math


} // gms











#endif /*__GMS_SIMD_UTILS_H__*/
