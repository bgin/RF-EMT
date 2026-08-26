
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

#ifndef __GMS_AVX512C16F32_H__
#define __GMS_AVX512C16F32_H__ 220220201228


namespace file_info {

 const unsigned int gGMS_AVX512C16F32_MAJOR = 1U;
 const unsigned int gGMS_AVX512C16F32_MINOR = 1U;
 const unsigned int gGMS_AVX512C16F32_MICRO = 0U;
 const unsigned int gGMS_AVX512C16F32_FULLVER =
  1000U*gGMS_AVX512C16F32_MAJOR+100U*gGMS_AVX512C16F32_MINOR+10U*gGMS_AVX512C16F32_MICRO;
 const char * const pgGMS_AVX512C16F32_CREATION_DATE = "22-02-2020 12:28 +00200 (SAT 22 FEB 2020 12:28 GMT+2)";
 const char * const pgGMS_AVX512C16F32_BUILD_DATE    = __DATE__ " " __TIME__ ;
 const char * const pgGMS_AVX512C16F32_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
 const char * const pgGMS_AVX512C16F32_SYNOPSIS      = "AVX512 complex number class decomposed into real and imaginary parts stored as 16-tuple;"


}

#include <cstdint>
#include <immintrin.h>
#include <complex>
#include "GMS_config.h"

namespace gms {

         namespace  math {

#if !defined(ZMM16c4_SETPS_ONE)
     define ZMM16c4_SETPS_ONE(x) _mm512_set_ps(1.0f,1.0f,1.0f,1.0f,1.0f,
                                                     1.0f,1.0f,1.0f,1.0f,1.0f,
						     1.0f,1.0f,1.0f,1.0f,1.0f,
						     (x));
#endif

                 struct ZMM16c4 {

		      __m512 m_re;
		      __m512 m_im;

		      __ATTR_HOT__
		      __ATTR_ALIGN__(16)
		      ZMM16c4() {
                         m_re = _mm512_setzero_ps();
		         m_im = _mm512_setzero_ps();
		    }

		      __ATTR_HOT__
		      __ATTR_ALIGN__(16)
		      __ATTR_VECTORCALL__
		      ZMM16c4(const float * __restrict __ATTR_ALIGN__(64) re,
		                   const float * __restrict __ATTR_ALIGN__(64) im) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                           re = (const float*)__builtin_assume_aligned(re,64);
			   im = (const float*)__builtin_assume_aligned(im,64);
#elif defined __ICC || defined __INTEL_COMPILER
                           __assume_aligned(re,64);
			   __assume_aligned(im,64);
#endif
                           m_re = _mm512_load_ps(&re[0]);
			   m_im = _mm512_load_ps(&im[0]);
		    }

		      __ATTR_HOT__
		      __ATTR_ALIGN__(16)
		      ZMM16c4(const float re,
		              const float im) {
                           m_re = _mm512_set1_ps(re);
			   m_im = _mm512_set1_ps(im);
		    }

		      __ATTR_HOT__
		      __ATTR_ALIGN__(16)
		      ZMM16c4(const std::complex<float> c) {
                           m_re = _mm512_set1_ps(c.real());
			   m_im = _mm512_set1_ps(c.imag());
		    }

		      __ATTR_HOT__
		      __ATTR_ALIGN__(16)
		      ZMM16c4(const float re) {
                           m_re = _mm512_set1_ps(re);
			   m_im = _mm512_setzero_ps();
		    }

		    __ATTR_HOT__
		    __ATTR_ALIGN__(16)
		    __ATTR_VECTORCALL__
		    ZMM16c4(const float re0,
		                 const float re1,
				 const float re2,
				 const float re3,
				 const float re4,
				 const float re5,
				 const float re6,
				 const float re7,
				 const float re8,
				 const float re9,
				 const float re10,
				 const float re11,
				 const float re12,
				 const float re13,
				 const float re14,
				 const float re15) {
                          m_re = _mm512_setr_ps(re15,re14,re13,re12,
			                        re11,re10,re9,re8,
						re7,re6,re5,re4,
						re3,re2,re1,re0);
			  m_im = _mm512_setzero_ps();
		   }

		   __ATTR_HOT__
		   __ATTR_ALIGN__(16)
		   __ATTR_VECTORCALL__
		   ZMM16c4( const float re0,
		                 const float re1,
				 const float re2,
				 const float re3,
				 const float re4,
				 const float re5,
				 const float re6,
				 const float re7,
				 const float re8,
				 const float re9,
				 const float re10,
				 const float re11,
				 const float re12,
				 const float re13,
				 const float re14,
				 const float re15,
				 const float im0,
				 const float im1,
				 const float im2,
				 const float im3,
				 const float im4,
				 const float im5,
				 const float im6,
				 const float im7,
				 const float im8,
				 const float im9,
				 const float im10,
				 const float im11,
				 const float im12,
				 const float im13,
				 const float im14,
				 const float im15) {
                              m_re = _mm512_setr_ps(re15,re14,re13,re12,
			                            re11,re10,re9,re8,
						    re7,re6,re5,re4,
						    re3,re2,re1,re0);
			      m_im = _mm512_setr_ps(im15,im14,im13,im12,
			                            im11,im10,im9,im8,
						    im7,im6,im5,im4,
						    im3,im2,im1,im0);
          
		  }

		    __ATTR_HOT__
		    __ATTR_ALIGN__(16)
		    __ATTR_VECTORCALL__
		    ZMM16c4(const __m512 re,
		                 const __m512 im) {
                          m_re = re;
			  m_im = im;
                  }


		   __ATTR_HOT__
		   __ATTR_ALIGN__(16)
		   __ATTR_VECTORCALL__
		     ZMM16c4(const __m512 re) {
		                
                          m_re = re;
			  m_im = _mm512_setzero_ps();
                  }

		  __ATTR_HOT__
		  __ATTR_ALIGN__(16)
		  __ATTR_VECTORCALL__
		  ZMM16c4(const ZMM16c4 & x) {
                          m_re = x.m_re;
                          m_im = x.m_im;
		  }

		  __ATTR_HOT__
		  __ATTR_ALIGN__(16)
		  __ATTR_VECTORCALL__
		  ZMM16c4 &
		  load_a(const float * __restrict __ATTR_ALIGN__(64) re,
		         const float * __restrict __ATTR_ALIGN__(64) im) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                         re = (const float*)__builtin_assume_aligned(re,64);
			 im = (const float*)__builtin_assume_aligned(im,64);
#elif defined __ICC || defined __INTEL_COMPILER
                         __assume_aligned(re,64);
			 __assume_aligned(im,64);
#endif
                       m_re = _mm512_load_ps(&re[0]);
		       m_im = _mm512_load_ps(&im[0]);
		       return (*this);
		 }

		 __ATTR_HOT__
		 __ATTR_ALIGN__(16)
		 __ATTR_VECTORCALL__
		 ZMM16c4 &
		 load_u(const float * __restrict re,
		        const float * __restrict im) {
                      m_re = _mm512_loadu_ps(&re[0]);
                      m_im = _mm512_loadu_ps(&im[0]);
                      return (*this);
		 }

		 __ATTR_HOT__
		 __ATTR_ALIGN__(16)
		 __ATTR_VECTORCALL__
		 void store_a(float * __restrict __ATTR_ALIGN__(64) re,
		              float * __restrict __ATTR_ALIGN__(64) im) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                         re = (float*)__builtin_assume_aligned(re,64);
			 im = (float*)__builtin_assume_aligned(im,64);
#elif defined __ICC || defined __INTEL_COMPILER
                         __assume_aligned(re,64);
			 __assume_aligned(im,64);
#endif
                          _mm512_store_ps(&re[0],m_re);
                          _mm512_store_ps(&im[0],m_im);
		}

		__ATTR_HOT__
		__ATTR_ALIGN__(16)
		__ATTR_VECTORCALL__
		void store_u(float * __restrict re,
		             float * __restrict im) {
                      _mm512_store_ps(&re[0],m_re);
                      _mm512_store_ps(&im[0],m_im);
	        }

	      __ATTR_HOT__
	      __ATTR_ALIGN__(16)
	      __ATTR_VECTORCALL__
	      void stream_nt(float * __restrict re,
	                     float * __restrict im) {
                     _mm512_stream_ps(&re[0],m_re);
                     _mm512_stream_ps(&im[0],m_im);
                    _mm_sfence();
	        }

	     __ATTR_COLD__
	     __ATTR_ALIGN__(32)
	     float extract_1xf32(const int32_t pos) {
                     __attribute__((aligned(64))) float mem[32] = {};
                     store_a(&mem[0],&mem[16]);
                      return (mem[pos & 0x1F]);
                }

	     __ATTR_COLD__
	     __ATTR_ALIGN__(32)
	     std::pair<float,float>
	     extract_2xf32(const int32_t posx,
	                   const int32_t posy) {
                     __attribute__((aligned(64))) float re_mem[16] = {};
                     __attribute__((aligned(64))) float im_mem[16] = {};
                     store_a(&re_mem[0],&im_mem[0]);
                     return (std::make_pair(re_mem[posx & 0xF],im_mem[posy  & 0xF]));
                }
               
	    __ATTR_COLD__
	    __ATTR_ALIGN__(32)
	    ZMM16c4 &
	    insert_1xf32(const int32_t pos,
	                 const float value) {
                    __attribute__((aligned(64))) float mem[32] = {};

                   store_a(&mem[0],&mem[16]);
                   mem[idx & 0x1F] = value;
                   load_a(&mem[0],&mem[16]);
                   return (*this);
	        }

	    __ATTR_COLD__
	    __ATTR_ALIGN__(32)
	    ZMM16c4 &
	    insert_2xf32(const int32_t re_idx
	                 const int32_t im_idx
			 const float re,
			 const float im) {
                 __attribute__((aligned(64))) float mem_re[16] = {};
                 __attribute__((aligned(64))) float mem_im[16] = {};

                 store_a(&mem_re[0],&mem_im[0]);
                 mem_re[re_idx & 0xF] = re;
                 mem_im[im_idx & 0xF] = im;
                 load_a(&mem_re[0],&mem_im[0]);
                 return (*this);
	       }

	    __ATTR_COLD__
            __ATTR_ALIGN__(16)
	    __ATTR_VECTORCALL__
	    void concatenate_a(float * __restrict __ATTR_ALIGN__(64) out) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                  out = (float*)__builtin_assume_aligned(out,64);
#elif defined __ICC || defined __INTEL_COMPILER
                  __assume_aligned(out,64);
#endif
                  store_a(&out[0],&out[16]);
	       }

            __ATTR_COLD__
	    __ATTR_ALIGN__(16)
	    __ATTR_VECTORCALL__
	    void concatenate_u(float * __restrict out) {
                  store_u(&out[0],&out[16]);
	      }

	    __ATTR_COLD__
	    __ATTR_ALIGN__(16)
	    ZMM16c4 &
	    partial_loadu(const float * __restrict re,
	                  const int32_t n_re,
			  const float * __restrict im,
			  const int32_t n_im) {
                   m_re = _mm512_maskz_loadu_ps(__mmask16((1 << n_re)-1),re);
                   m_im = _mm512_maskz_loadu_ps(__mmask16((1 << n_im)-1),im);
                   return (*this);
	      }

	    __ATTR_COLD__
	    __ATTR_ALIGN__(16)
	    ZMM16c4 &
	    partial_loada(const float * __restrict re,
	                  const int32_t n_re,
			  const float * __restrict im,
			  const int32_t n_im) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                   re = (const float*)__builtin_assume_aligned(re,64);
		   im = (const float*)__builtin_assume_aligned(im,64);
#elif defined __ICC || defined __INTEL_COMPILER
                   __assume_aligned(re,64);
		   __assume_aligned(im,64);
#endif
                    m_re = _mm512_maskz_load_ps(__mmask16((1 << n_re)-1),re);
                    m_im = _mm512_maskz_load_ps(__mmask16((1 << n_im)-1),im);
                    return (*this);
	      }

	      __ATTR_COLD__
	      __ATTR_ALIGN__(16)
	      void partial_storeu(float * __restrict re,
                                       const int32_t n_re,
				       float * __restrict im,
				       const int32_t n_im) {
                 _mm512_mask_storeu_ps(&re[0],__mmask16((1 << n_re)-1),m_re);
                 _mm512_mask_storeu_ps(&im[0],__mmask16((1 << m_im)-1),m_im);
              }

	      __ATTR_COLD__
	      __ATTR_ALIGN__(16)
	      void partial_storea(float * __restrict __ATTR_ALIGN__(64) re,
	                          float * __restrict __ATTR_ALIGN__(64) im) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                   re = (float*)__builtin_assume_aligned(re,64);
		   im = (float*)__builtin_assume_aligned(im,64);
#elif defined __ICC || defined __INTEL_COMPILER
                   __assume_aligned(re,64);
		   __assume_aligned(im,64);
#endif
                    _mm512_mask_store_ps(&re[0],__mmask16((1 << n_re)-1),m_re);
                    _mm512_mask_store_ps(&im[0],__mmask16((1 << m_im)-1),m_im);
	      }

	      __ATTR_COLD__
	      __ATTR_ALIGN__(16)
	      __ATTR_VECTORCALL__
	      ZMM16c4 &
	             expand(const ZMM16c4 x,
                            const __mmask16 mask) {
                     m_re = _mm512_maskz_expand_ps(mask,x.m_re);
                     m_im = _mm512_maskz_expand_ps(mask.x.m_im);
                     return (*this);
              }

	      __ATTR_COLD__
	      __ATTR_ALIGN__(16)
	      __ATTR_VECTORCALL__
	      ZMM16c4 &
	            expand_loadu(const ZMM16c4 x,
                                 const __mmask16 mask,
				 const float * __restrict re,
				 const float * __restrict im) {
                     m_re = _mm512_mask_expandloadu_ps(x.m_re,mask,&re[0]);
                     m_im = _mm512_mask_expandloadu_ps(x.m_im,mask,&im[0]);
                     return (*this);
               }

	       __ATTR_COLD__
	       __ATTR_ALIGN__(16)
	       ZMM16c4 & permute(const __mmask16 mask,
                                      const int32_t imm) {
                    m_re = _mm512_mask_permute_ps(m_re,mask,m_im,imm);
                    m_im = _mm512_mask_permute_ps(m_im,mask,m_re,imm);
                    return (*this);
               }

	       __ATTR_COLD__
	       __ATTR_ALIGN__(16)
	       __ATTR_VECTORCALL__
	       __m256 re_8xf32() const {
                     return (_mm512_extractf32x8_ps(m_re,0));
	       }

	       __ATTR_COLD__
	       __ATTR_ALIGN__(16)
	       __ATTR_VECTORCALL__
	       __m256 im_8xf32() const {
                      return (_mm512_extract32x8_ps(m_im,0));
	       }

	       __ATTR_HOT__
	       __ATTR_ALIGN__(16)
	       __ATTR_VECTORCALL__
	       ZMM16c4 &
	              operator=(const ZMM16c4 x) {
                        if(this == &x) return (*this);
                        m_re = x.m_re;
                        m_im = x.m_im;
                        return (*this);
	       }

	       
	       
	  }__ATTR_ALIGN__(64);


	     __ATTR_HOT__
	     __ATTR_ALIGN__(16)
	     __ATTR_VECTORCALL__
	     static inline
	     ZMM16c4 conj(ZMM16c4 x) {
	       auto tmp = ~x;
		  return (ZMM16c4{x.m_re,tmp.m_im});
	     }

	     __ATTR__HOT__
	     __ATTR_ALIGN__(32)
	     __ATTR_VECTORCALL__
	     static inline
	     ZMM16c4 polar(const __m512 rho,
	                        const __m512 theta) {
                   const __m512 re_part =
                        _mm512_mul_ps(rho,_mm512_cos_ps(theta));
                   const __m512 im_part =
                        _mm512_mul_ps(rho,_mm512_sin_ps(theta));
                   return (ZMM16c4{re_part,im_part});
	     }

	     __ATTR_HOT__
	     __ATTR_ALIGN__(16)
	     __ATTR_VECTORCALL__
	     static inline
	     __m512 carg(const ZMM16c4 x) {
                  return (_mm512_atan2_ps(x.m_re,x.m_im));
	     }

	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     __ATTR_VECTORCALL__
	     static inline
	     ZMM16c4 csin(const ZMM16c4 x) {
                    const __m512 re_part =
                             _mm512_mul_ps(_mm512_sin_ps(x.m_re),
	                                   _mm512_cosh_ps(x.m_im));
                    const __m512 im_part =
                             _mm512_mul_ps(_mm512_cos_ps(x.m_re),
	                                   _mm512_sinh_ps(x.m_im));
                    return (ZMM16c4{re_part,im_part});
	     }

	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     __ATTR_VECTORCALL__
	     static inline
	     ZMM16c4 csinh(const ZMM16c4 x) {
                    const __m512 re_part =
                             _mm512_mul_ps(_mm512_sinh_ps(x.m_re),
	                                   _mm512_cos_ps(x.m_im));
                    const __m512 im_part =
                             _mm512_mul_ps(_mm512_cosh_ps(x.m_re),
	                                   _mm512_sin_ps(x.m_im));
                    return (ZMM16c4{re_part,im_part});
	     }

	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     __ATTR_VECTORCALL__
	     static inline
	     ZMM16c4 ccos(const ZMM16c4 x) {
                    const __m512 re_part =
                                _mm512_mul_ps(_mm512_cos_ps(x.m_re),
	                                      _mm512_cosh_ps(x.m_im));
                    const __m512 im_part =
                                _mm512_mul_part(_mm512_sin_ps(x.m_re),
	                                      _mm512_sinh_ps(x.m_im));
                    return (ZMM16c4{re_part,im_part});
	     }
	     
             __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     __ATTR_VECTORCALL__
	     static inline
	     ZMM16c4 ccosh(const ZMM16c4 x) {
                    const __m512 re_part =
                                _mm512_mul_ps(_mm512_cosh_ps(x.m_re),
	                                      _mm512_cos_ps(x.m_im));
                    const __m512 im_part =
                                _mm512_mul_part(_mm512_sinh_ps(x.m_re),
	                                      _mm512_sin_ps(x.m_im));
                    return (ZMM16c4{re_part,im_part});
	     }

	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     __ATTR_VECTORCALL__
	     static inline
	     ZMM16c4 cexp(const ZMM16c4 x) {
                    const __m512 re_part =
                                 _mm512_mul_ps(_mm512_exp_ps(x.m_re),
	                                       _mm512_cos_ps(x.m_im));
                    const __m512 im_part =
                                 _mm512_mul_ps(_mm512_exp_ps(x.m_re),
	                                       _mm512_sin_ps(x.m_im));
                    return (ZMM16c4{re_part,im_part});
	     }

	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     __ATTR_VECTORCALL__
	     static inline
	     __m512 cnorm(const ZMM16c4 x) {
                    const __m512 re_part =
                             _mm512_mul_ps(x.m_re,x.m_re);
                    const __m512 im_part =
                             _mm512_mul_ps(x.m_im,x.m_im);
                    return (_mm512_sqrt_ps(_mm512_add_ps(re_part,im_part)));
             }
		 
		 
             

	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     __ATTR_VECTORCALL__
	     static inline
	     ZMM16c4 cpow(const ZMM16c4 x,
                               const float n) {
                    const __m512 re_part =
                               _mm512_mul_ps(x.m_re,x.m_re);
                    const __m512 im_part =
                               _mm512_mul_ps(x.m_im,x.m_im);
                    const __m512 r =
                               _mm512_sqrt_ps(_mm512_add_ps(re_part,im_part));
                    const __m512 theta =
                               _mm512_atan_ps(_mm512_div_ps(x.m_im,x.m_re));
                    const __m512 vn = _mm512_set1_ps(n);
                    const __m512 pow_term = _mm512_pow_ps(r,vn);
                    const __m512 trig_arg = _mm512_mul_ps(vn,theta);
                    return (ZMM16c4{_mm512_mul_ps(pow_term,_mm512_cos_ps(trig_arg),
						  _mm512_mul_ps(pow_term,_mm512_sin_ps(trig_arg)))});
             }

	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     __ATTR_VECTORCALL__
	     static inline
	     ZMM16c4 clog(const ZMM16c4 x) {
                     const __m512 t1  = cabs(x);
                     const __m512 t2  = carg(x);
                     const __m512 re_part = _mm512_log_ps(t1);
                     return (ZMM16c4{re_part,t2});
	     }

	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     __ATTR_VECTORCALL__
	     static inline
	     ZMM16c4 csqrt(const ZMM16c4 x) {
                       const __m512 t = cabs(x);
                       const __m512 re_part =
                               _mm512_mul_ps(_mm512_set1_ps(0.5f),
	                                     _mm512_add_ps(t,x.m_re));
                       const __m512 im_part =
                               _mm512_mul_ps(_mm512_set1_ps(0.5f),
	                                     _mm512_sub_ps(t,x.m_re));
                       return (ZMM16c4{_mm512_sqrt_ps(re_part),
                                            _mm512_sqrt_ps(im_part)});
	     }

	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     __ATTR_VECTORCALL__
	     static inline
	     ZMM16c4 ctan(const ZMM16c4 x) {
                       return (csin(x)/ccos(x));
	     }

	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     __ATTR_VECTORCALL__
             static inline
	     ZMM16c4
             ctanh(const ZMM16c4 x) {
                        return (csinh(x)/ccosh(x));
             }

	     __ATTR_HOT__
	     __ATTR_ALIGN__(16)
	     __ATTR_VECTORCALL__
             static inline
	     ZMM16c4
             select(const ZMM16c4 x,
                    const ZMM16c4 y,
		    __mmask16 mask) {
                 return (ZMM16c4{_mm512_mask_blend_ps(mask,x.m_re,y.m_re),
                                     _mm512_mask_blend_ps(mask,x.m_im,y.m_im)});
              }

	 /*    __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     __ATTR_VECTORCALL__
             static inline
	    ZMM16c4
             cdiv_smith(const ZMM16c4 x,
                        const ZMM16c4 y) {
                    __m512 ratio,denom,re_part,im_part;
                    constexpr __mmask16 all_ones = 0xFFFF;
                    re_part = _mm512_setzero_ps();
                    im_part = _mm512_setzero_ps();
                    __mmask16 is_gte = _mm512_cmp_ps_mask(
                                _mm512_abs_ps(y.m_re),
			        _mm512_abs_ps(y.m_im),_CMP_GE_OQ);
                    ratio = _mm512_setzero_ps();
                    denom = _mm512_setzero_ps();
                    if(is_gte == all_ones) {
                        ratio = _mm512_div_ps(y.m_im,y.m_re);
	                denom = _mm512_add_ps(y.m_re,
	                         _mm512_mul_ps(ratio,y.m_im));
	                re_part = _mm512_div_ps(_mm512_add_ps(x.m_re,
	                          _mm512_mul_ps(x.m_im,ratio)),denom);
	                im_part = _mm512_div_ps(_mm512_sub_ps(x.m_im,
	                          _mm512_mul_ps(x.m_re,ratio)),denom);
	                return (ZMM16c4{re_part,im_part});
                     }
                     else {
                        ratio = _mm512_div_ps(y.m_re,y.m_im);
	                denom = _mm512_add_ps(y.m_im,_mm512_mul_ps(ratio,y.m_re));
	                re_part = _mm512_div_ps(_mm512_add_ps(
	                          _mm512_mul_ps(x.m_re,ratio)),denom);
	                im_part = _mm512_div_ps(_mm512_sub_ps(
	                          _mm512_mul_ps(x.m_im,ratio)),denom);
	                return (ZMM16c4{re_part,im_part});
                      }
               } */

	      __ATTR_HOT__
	      __ATTR_ALIGN__(16)
	      __ATTR_VECTORCALL__
              static inline
	     ZMM16c4
              operator+(const ZMM16c4 x,
                        const ZMM16c4 y) {
                  return (ZMM16c4{_mm512_add_ps(x.m_re,y.m_re),
                                       _mm512_add_ps(x.m_im,x.m_im)});
               }

	       __ATTR_HOT__
	       __ATTR_ALIGN__(16)
	       __ATTR_VECTORCALL__
	       static inline
	       ZMM16c4
	       operator+(const AVX512c16f32 x,
	                 const std::complex<float> y) {
                   return(ZMM16c4{_mm512_add_ps(x.m_re,_mm512_set1_ps(y.real())),
		                       _mm512_add_ps(x.m_im,_mm512_set1_ps(y.imag()))});
	       }

	       __ATTR_HOT__
	       __ATTR_ALIGN__(16)
	       __ATTR_VECTORCALL__
	       static inline
	      ZMM16c4
	       operator+(const std::complex<float> x,
	                 const ZMM16c4 y) {
                   return (ZMM16c4{_mm512_add_ps(_mm512_set1_ps(x.real()),y.m_re),
		                        _mm512_add_ps(_mm512_set1_ps(x.imag()),y.m_im)});
		}

	      __ATTR_HOT__
	      __ATTR_ALIGN__(16)
	      __ATTR_VECTORCALL__
              static inline ZMM16c4
              operator+(const ZMM16c4 x,
                        const __m512 v) {
                  return (ZMM16c4{_mm512_add_ps(x.m_re,v),
                                                    x.m_im});
               }

              __ATTR_HOT__
	      __ATTR_ALIGN__(16)
	      __ATTR_VECTORCALL__
              static inline ZMM16c4
              operator+(const __m512 v,
                        const ZMM16c4 x) {
                   return (ZMM16c4{_mm512_add_ps(v,x.m_re),x.m_im});
               }

              __ATTR_HOT__
	      __ATTR_ALIGN__(16)
	      __ATTR_VECTORCALL__
              static inline ZMM16c4
              operator+(const ZMM16c4 x,
                        const float s) {
                   return (x + ZMM16c4{s});
               }

              __ATTR_HOT__
	      __ATTR_ALIGN__(16)
	      __ATTR_VECTORCALL__
              static inline ZMM16c4
              operator+(const float s,
                       const ZMM16c4 x) {
                   return (ZMM16c4{s} + x);
               }

	     __ATTR_HOT__
	     __ATTR_ALIGN__(16)
	     __ATTR_VECTORCALL__
             static inline ZMM16c4
             operator+=(ZMM16c4 x,
                        const ZMM16c4 y) {
                  x = x + y;
                  return (x)
              }

	     __ATTR_HOT__
	     __ATTR_ALIGN__(16)
	     __ATTR_VECTORCALL__
	     static inline ZMM16c4
	     operator+=(ZMM16c4 x,
	                const std::complex<float> y) {
                 x = x + y;
		 return (x);
	     }

	     __ATTR_HOT__
	     __ATTR_ALIGN__(16)
	     __ATTR_VECTORCALL__
	     static inline ZMM16c4
	     operator+=(const std::complex<float> x,
	                ZMM16c4 y) {
                 y = y + x;
		 return (y);
	     }

             __ATTR_HOT__
	     __ATTR_ALIGN__(16)
	     __ATTR_VECTORCALL__
             static inline ZMM16c4
             operator+=(ZMM16c4 x,
                        const __m512 v) {
                 x = x + v;
                 return (x);
              }

             __ATTR_HOT__
	     __ATTR_ALIGN__(16)
	     __ATTR_VECTORCALL__
             static inline ZMM16c4
             operator+=(const __m512 v,
                        ZMM16c4 x) {
                x = v + x;
                return (x);
              }

             __ATTR_HOT__
	     __ATTR_ALIGN__(16)
	     __ATTR_VECTORCALL__
             static inline ZMM16c4
             operator+=(ZMM16c4 x,
                        const float s) {
                  x = x + ZMM16c4{s};
                  return (x)
              }

             __ATTR_HOT__
	     __ATTR_ALIGN__(16)
	     __ATTR_VECTORCALL__
             static inline ZMM16c4
             operator+=(const float s,
                       ZMM16c4 x) {
                x = ZMM16c42{s} + x;
                return (x);
              }

	     __ATTR_HOT__
	     __ATTR_ALIGN__(16)
	     __ATTR_VECTORCALL__
             static inline ZMM16c4
             operator-(const ZMM16c4 x,
                     const ZMM16c4 y) {
                  return (ZMM16c4{_mm512_sub_ps(x.m_re,y.m_re),
                                      _mm512_sub_ps(x.m_im,y.m_im)});
              }

	    __ATTR_HOT__
	    __ATTR_ALIGN__(16)
	    __ATTR_VECTORCALL__
	    static inline ZMM16c4
	    operator-(const ZMM16c4 x,
	              const std::complex<float> y) {
                 return(ZMM16c4{_mm512_sub_ps(x.m_re,_mm512_set1_ps(y.real())),
		                     _mm512_sub_ps(x.m_im,_mm512_set1_ps(y.imag()))});
	    }

	    __ATTR_HOT__
	    __ATTR_ALIGN__
	    __ATTR_VECTORCALL__
	    static inlineZMM16c4
	    operator-(const std::complex<float> x,
	              const ZMM16c4 y) {
                  return (ZMM16c4{_mm512_sub_ps(_mm512_set1_ps(x.real()),y.m_re),
		                     _mm512_sub_ps(_mm512_set1_ps(x.imag()),y.m_im)});
	    }
	    

            __ATTR_HOT__
	    __ATTR_ALIGN__(16)
	    __ATTR_VECTORCALL__
            static inline ZMM16c4
            operator-(const ZMM16c4 x,
                      const __m512 v) {
                return (ZMM16c4{_mm512_sub_ps(x.m_re,v),x.m_im});
             }

            __ATTR_HOT__
	    __ATTR_ALIGN__(16)
	    __ATTR_VECTORCALL__
            static inline ZMM16c4
            operator-(const __m512 v,
                     const ZMM16c4 x) {
                  return (ZMM16c4{_mm512_sub_ps(v,x.m_re),x.m_im});
             }

            __ATTR_HOT__
	    __ATTR_ALIGN__(16)
	    __ATTR_VECTORCALL__
            static inline ZMM16c4
            operator-(const ZMM16c4 x,
                     const float s) {
                 return (x - ZMM16c4{s});
             }

            __ATTR_HOT__
	    __ATTR_ALIGN__(16)
	    __ATTR_VECTORCALL__
            static inline ZMM16c4
            operator-(const float s,
                     const ZMM16c4 x) {
                  return (ZMM16c4{s} - x);
              }

	     __ATTR_HOT__
	    __ATTR_ALIGN__(16)
	    __ATTR_VECTORCALL__
            static inline ZMM16c4
	    operator-(const ZMM16c4 x) {
               const __m512 _0 = _mm512_setzero_ps();
               return (ZMM16c4{_mm512_sub_ps(_0,x.re),
	                       _mm512_sub_ps(_0,x.im)});
	    }

	    __ATTR_HOT__
	    __ATTR_ALIGN__(16)
	    __ATTR_VECTORCALL__
	    static inline ZMM16c4
	    operator-=(ZMM16c4 x,
	               const ZMM16c4 y) {
                x = x - y;
		return (x);
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline ZMM16c4
	   operator-=(ZMM16c4 x,
	              const std::complex<float> y) {
                x = x - y;
		return (x);
	    }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline ZMM16c4
	   operator-=(const std::complex<float> x,
	              ZMM16c4 y) {
               y = x - y;
	       return (y);
	    }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline ZMM16c4
	   operator-=(ZMM16c4 x,
	              const __m512 y) {
              x = x - y;
	      return (x);
	   }

	  __ATTR_HOT__
	  __ATTR_ALIGN__(16)
	  __ATTR_VECTORCALL__
	  static inline ZMM16c4
	  operator-=(const __m512 x,
	            ZMM16c4 y) {
              y = x - y;
	      return (y);
	  }

	 __ATTR_HOT__
	 __ATTR_ALIGN__(16)
	 __ATTR_VECTORCALL__
	 static inline ZMM16c4
	 operator-=(ZMM16c4 x,
	            const float s) {
            x = x - ZMM16c4{s};
	    return (x);
	 }

	 __ATTR_HOT__
	 __ATTR_ALIGN__(16)
	 __ATTR_VECTORCALL__
	 static inline ZMM16c4
	 operator-=(const float s,
	            ZMM16c4 x) {
            x = ZMM16c4{s}-x;
	    return (x);
	 }

	__ATTR_HOT__
	__ATTR_ALIGN__(16)
	__ATTR_VECTORCALL__
        static inline ZMM16c4
        operator*(const ZMM16c4 x,
                  const ZMM16c4 y) {
              const __m512 zmm0 = _mm512_mul_ps(x.m_re,y.m_re);
              const __m512 zmm1 = _mm512_mul_ps(x.m_im,y.m_im);
              const __m512 zmm2 = _mm512_mul_ps(x.m_im,y.m_re);
              const __m512 zmm3 = _mm512_mul_ps(x.m_re,y.m_im);
              return (ZMM16c4{_mm512_add_ps(zmm0,zmm1),
                                   _mm512_sub_ps(zmm2,zmm3)});
           }

	 __ATTR_HOT__
	 __ATTR_ALIGN__(16)
	 __ATTR_VECTORCALL__
	 static inline ZMM16c4
	 operator*(const ZMM16c4 x,
	           const std::complex<float> y) {
            const __m512 y_re = _mm512_set1_ps(y.real());
	    const __m512 y_im = _mm512_set1_ps(y.imag());
	    const __m512 zmm0 = _mm512_mul_ps(x.m_re,y_re);
	    const __m512 zmm1 = _mm512_mul_ps(x.m_im,y_im);
	    const __m512 zmm2 = _mm512_mul_ps(x.m_im,y_re);
	    const __m512 zmm3 = _mm512_mul_ps(x.m_re,y_im);
	    return (ZMM16c4{_mm512_add_ps(zmm0,zmm1),
                                 _mm512_sub_ps(zmm2,zmm3)});
	 }

        __ATTR_HOT__
	__ATTR_ALIGN__(16)
	__ATTR_VECTORCALL__
        static inline ZMM16c4
        operator*(const ZMM16c4 x,
                     const __m512 v) {
              return (ZMM16c4{_mm512_mul_ps(x.m_re,v),
                                  _mm512_mul_ps(x.m_im,v)});
          }

       __ATTR_HOT__
       __ATTR_ALIGN__(16)
       __ATTR_VECTORCALL__
       static inline ZMM16c4
       operator*(const __m512 v,
                 const ZMM16c4 x) {
            return (ZMM16c4{_mm512_mul_ps(v,x.m_re),
                                _mm512_mul_ps(v,x.m_im)});
          }

       __ATTR_HOT__
       __ATTR_ALIGN__(16)
       __ATTR_VECTORCALL__
       static inline ZMM16c4
       operator*(const ZMM16c4 x,
                 const float s) {
         const __m512 zmm0 = _mm512_set1_ps(s);
              return (ZMM16c4{_mm512_mul_ps(x.m_re,zmm0),
                                  _mm512_mul_ps(x.m_im,zmm0)});
          }

       __ATTR_HOT__
       __ATTR_ALIGN__(16)
       __ATTR_VECTORCALL__
       static inline ZMM16c4
       operator*=(ZMM16c4 x,
                  const ZMM16c4 y) {
              x = x * y;
              return (x);
         }

       __ATTR_HOT__
       __ATTR_ALIGN__(16)
       __ATTR_VECTORCALL__
       static inline ZMM16c4
       operator*=(ZMM16c4 x,
                  const __m512 v) {
            x = x * v;
            return (x);
          }

       __ATTR_HOT__
       __ATTR_ALIGN__(16)
       __ATTR_VECTORCALL__
       static inline ZMM16c4
       operator*=(const __m512 v,
                   ZMM16c4 x) {
              x = v * x;
              return (x);
          }

        __ATTR_HOT__
	__ATTR_ALIGN__(16)
	__ATTR_VECTORCALL__
        static inline ZMM16c4
        operator*=(ZMM16c4 x,
                   const float s) {
             x = x * s;
             return (x);
          }

        __ATTR_HOT__
	__ATTR_ALIGN__(16)
	__ATTR_VECTORCALL__
         static inline ZMM16c4
         operator*=(const float s,
                   ZMM16c4 x) {
             x = s * x;
             return (x);
          }

	__ATTR_HOT__
	__ATTR_ALIGN__(32)
	__ATTR_VECTORCALL__
        static inline ZMM16c4
        operator/(const ZMM16c4 x,
                  const ZMM16c4 y) {
//#if defined USE_SAFE_COMPLEX_DIVISION && (USE_SAFE_COMPLEX_DIVISION) == 1
//              return (cdiv_smith(x,y));
//#else
              const __m512 zmm0 = _mm512_mul_ps(x.m_re,y.m_re);
              const __m512 zmm1 = _mm512_mul_ps(x.m_im,y.m_im);
              const __m512 zmm2 = _mm512_mul_ps(x.m_im,y.m_re);
              const __m512 zmm3 = _mm512_mul_ps(x.m_re,y.m_im);
              const __m512 den  = _mm512_add_ps(_mm512_mul_ps(y.m_re,y.m_re),
                                                _mm512_mul_ps(y.m_im,y.m_im));
              const __m512 re_part = _mm512_add_ps(zmm0,zmm1);
              const __m512 im_part = _mm512_sub_ps(zmm2,zmm3);
              return (ZMM16c4{_mm512_div_ps(re_part,den),
                                  _mm512_div_ps(im_part,den)});
//#endif
           }

	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 __ATTR_VECTORCALL__
	 static inline  ZMM16c4
	 operator/(const  ZMM16c4 x,
	           const std::complex<float> y) {

              const __m512 y_re = _mm512_set1_ps(y.real());
	      const __m512 y_im = _mm512_set1_ps(y.imag());
	      const __m512 zmm0 = _mm512_mul_ps(x.m_re,y_re);
              const __m512 zmm1 = _mm512_mul_ps(x.m_im,y_im);
              const __m512 zmm2 = _mm512_mul_ps(x.m_im,y_re);
              const __m512 zmm3 = _mm512_mul_ps(x.m_re,y_im);
              const __m512 den  = _mm512_add_ps(_mm512_mul_ps(y_re,y_re),
                                                _mm512_mul_ps(y_im,y_im));
              const __m512 re_part = _mm512_add_ps(zmm0,zmm1);
              const __m512 im_part = _mm512_sub_ps(zmm2,zmm3);
              return ( ZMM16c4{_mm512_div_ps(re_part,den),
                                  _mm512_div_ps(im_part,den)});
	  }

	  __ATTR_HOT__
	  __ATTR_ALIGN__(32)
	  __ATTR_VECTORCALL__
	  static inline  ZMM16c4
	  operator/(const std::complex<float> x,
	            const  ZMM16c4 y) {

	      const __m512 x_re = _mm512_set1_ps(x.real());
	      const __m512 x_im = _mm512_set1_ps(x.imag());
	      const __m512 zmm0 = _mm512_mul_ps(x_re,y.m_re);
              const __m512 zmm1 = _mm512_mul_ps(x_im,y.m_im);
              const __m512 zmm2 = _mm512_mul_ps(x_im,y.m_re);
              const __m512 zmm3 = _mm512_mul_ps(x_re,y.m_im);
              const __m512 den  = _mm512_add_ps(_mm512_mul_ps(y.m_re,y.m_re),
                                                _mm512_mul_ps(y.m_im,y.m_im));
              const __m512 re_part = _mm512_add_ps(zmm0,zmm1);
              const __m512 im_part = _mm512_sub_ps(zmm2,zmm3);
              return ( ZMM16c4{_mm512_div_ps(re_part,den),
                                  _mm512_div_ps(im_part,den)});
	  }

	  __ATTR_HOT__
	  __ATTR_ALIGN__(16)
	  __ATTR_VECTORCALL__
          static inline  ZMM16c4
          operator/(const  ZMM16c4 x,
                     const __m512 v) {
                return ( ZMM16c4{_mm512_div_ps(x.m_re,v),
                                    _mm512_div_ps(x.m_im,v)});
           }

         __ATTR_HOT__
	 __ATTR_ALIGN__(16)
	 __ATTR_VECTORCALL__
         static inline  ZMM16c4
         operator/(const __m512 v,
                   const  ZMM16c4 x) {
                return ( ZMM16c4{_mm512_div_ps(v,x.m_re),
                                     _mm512_div_ps(v,x.m_im)});
           }

         __ATTR_HOT__
	 __ATTR_ALIGN__(16)
	 __ATTR_VECTORCALL__
        static inline ZMM16c4
        operator/(const ZMM16c4 x,
                     const float s) {
               const __m512 zmm0 = _mm52_set1_ps(s);
               return (ZMM16c4{_mm512_div_ps(x.m_re,zmm0),
                                    _mm512_div_ps(x.m_im,zmm0)});
           }

         __ATTR_HOT__
	 __ATTR_ALIGN__(16)
	 __ATTR_VECTORCALL__
        static inline ZMM16c4
        operator/(const float s,
                  const ZMM16c4 x) {
               const __m512 zmm0 = _mm512_set1_ps(s);
               return (ZMM16c4{_mm512_div_ps(zmm0,x.m_re),
                                    _mm512_div_ps(zmm0,x.m_im)});
           }

	__ATTR_HOT__
	__ATTR_ALIGN__(16)
	__ATTR_VECTORCALL__
        static inline ZMM16c4
        operator/=(AVX512c16f32 x,
                   const ZMM16c4 y) {
              x = x / y;
              return (x);
          }

	__ATTR_HOT__
	__ATTR_ALIGN__(16)
	__ATTR_VECTORCALL__
	static inline ZMM16c4
	operator/=(ZMM16c4 x,
	           const std::complex<float> y) {
              x = x / y
	      return (x);
	 }

	 __ATTR_HOT__
	 __ATTR_ALIGN__(16)
	 __ATTR_VECTORCALL__
	 static inline ZMM16c4
	 operator/=(const std::complex<float> x,
	            ZMM16c4 y) {
               y = x / y;
               return (y);
	 }

        __ATTR_HOT__
	__ATTR_ALIGN__(16)
	__ATTR_VECTORCALL__
        static inline ZMM16c4
        operator/=(ZMM16c4 x,
                      const __m512 v) {
              x = x / v;
              return (x);
          }

         __ATTR_HOT__
	 __ATTR_ALIGN__(16)
	 __ATTR_VECTORCALL__
         static inline ZMM16c4
         operator/=(const __m512 v,
                    ZMM16c4 x) {
             x = v / x;
             return (x);
          }

         __ATTR_HOT__
	 __ATTR_ALIGN__(16)
	 __ATTR_VECTORCALL__
         static inline ZMM16c4
         operator/=(ZMM16c4 x,
                    const float s) {
              x = x / s;
              return (x);
          }

	 __ATTR_HOT__
	 __ATTR_ALIGN__(16)
	 __ATTR_VECTORCALL__
         static inline ZMM16c4
         operator/=(const float s,
                    ZMM16c4 x) {
              x = s / x;
              return (x);
          }

          __ATTR_HOT__
	  __ATTR_ALIGN__(16)
	  __ATTR_VECTORCALL__
          static inline ZMM16c4
          operator~(ZMM16c4 x) {
                x.m_re = _mm512_sub_ps(_mm512_setzero_ps(),x.m_re);
                return (x);
           }

         __ATTR_HOT__
	 __ATTR_ALIGN__(16)
	 __ATTR_VECTORCALL__
         static inline MMASK16_RETTYPE
         operator==(const ZMM16c4 x,
                    const ZMM16c4 y) {
               const __mmask16 m1 =
                     _mm512_cmp_ps_mask(x.m_re,y.m_re,_CMP_EQ_OQ);
               const __mmask16 m2 =
                     _mm512_cmp_ps_mask(x.m_im,y.m_im,_CMP_EQ_OQ);
               return (std::make_pair(m1,m2));
           }

	 __ATTR_HOT__
	 __ATTR_ALIGN__(16)
	 __ATTR_VECTORCALL__
         static inline MMASK16_RETTYPE
         operator==(const ZMM16c4 x,
                    const std::complex<float> c) {
               const __mmask16 m1 =
                  _mm512_cmp_ps_mask(x.m_re,_mm512_set1_ps(c.real()),_CMP_EQ_OQ);
               const __mmask16 m2 =
                   _mm512_cmp_ps_mask(x.m_im,_mm512_set1_ps(c.imag()),_CMP_EQ_OQ);
               return (std::make_pair(m1,m2));
           }

         __ATTR_HOT__
	 __ATTR_ALIGN__(16)
	 __ATTR_VECTORCALL__
         static inline MMASK16_RETTYPE
         operator==(const std::complex<float> c,
                    const ZMM16c4 x) {
              const __mmask16 m1 =
                    _mm512_cmp_ps_mask(_mm512_set1_ps(c.real()),x.m_re,_CMP_EQ_OQ);
              const __mmask16 m2 =
                    _mm512_cmp_ps_mask(_mm512_set1_ps(c.imag()),x.m_im,_CMP_EQ_OQ);
              return (std::make_pair(m1,m2));
           }

          __ATTR_HOT__
	  __ATTR_ALIGN__(16)
	  __ATTR_VECTORCALL__
          static inline MMASK16_RETTYPE
          operator!=(const ZMM16c4 x,
                      const ZMM16c4 y) {
             const __mmask16 m1 =
                   _mm512_cmp_ps_mask(x.m_re,y.m_re,_CMP_NEQ_OQ);
             const __mmask16 m2 =
                   _mm512_cmp_ps_mask(x.m_im,y.m_im,_CMP_NEQ_OQ);
             return (std::make_pair(m1,m2));
           }

          __ATTR_HOT__
	  __ATTR_ALIGN__(16)
	  __ATTR_VECTORCALL__
          static inline MMASK16_RETTYPE
          operator!=(const ZMM16c4 x,
                      const std::complex<float> c) {
             const __mmask16 m1 =
                    _mm512_cmp_ps_mask(x.m_re,_mm512_set1_ps(c.real()),_CMP_NEQ_OQ);
             const __mmask16 m2 =
                    _mm512_cmp_ps_mask(x.m_im,_mm512_set1_ps(c.imag()),_CMP_NEQ_OQ);
             return (std::make_pair(m1,m2));
            }

           __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
           static inline MMASK16_RETTYPE
           operator!=(const std::complex<float> c,
                      const ZMM16c4 x) {
                const __mmask16 m1 =
                   _mm512_cmp_ps_mask(_mm512_set1_ps(c.real()),x.m_re,_CMP_NEQ_OQ);
                const __mmask16 m2 =
                   _mm512_cmp_ps_mask(_mm512_set1_ps(c.imag()),x.m_im,_CMP_NEQ_OQ);
                return (std::make_pair(m1,m2));
             }

           __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
           static inline MMASK16_RETTYPE
           operator>(const ZMM16c4 x,
                      const ZMM16c4 y) {
              const __mmask16 m1 =
                  _mm512_cmp_ps_mask(x.m_re,y.m_re,_CMP_GT_OQ);
              const __mmask16 m2 =
                  _mm512_cmp_ps_mask(x.m_im,y.m_im,_CMP_GT_OQ);
               return (std::make_pair(m1,m2));
             }

          __ATTR_HOT__
	  __ATTR_ALIGN__(16)
	  __ATTR_VECTORCALL__
          static inline MMASK16_RETTYPE
          operator>(const ZMM16c4 x,
                      const std::complex<float> c) {
              const __mmask16 m1 =
                   _mm512_cmp_ps_mask(x.m_re,_mm512_set1_ps(c.real()),_CMP_GT_OQ);
              const __mmask16 m2 =
                   _mm512_cmp_ps_mask(x.m_im,_mm512_set1_ps(c.imag()),_CMP_GT_OQ);
              return (std::make_pair(m1,m2));
             }

           __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
           static inline MMASK16_RETTYPE
           operator>(const std::complex<float> c,
                      const ZMM16c4 x) {
              const __mmask16 m1 =
                 _mm512_cmp_ps_mask(_mm512_set1_ps(c.real()),x.m_re,_CMP_GT_OQ);
              const __mmask16 m2 =
                 _mm512_cmp_ps_mask(_mm512_set1_ps(c.imag()),x.m_im,_CMP_GT_OQ);
              return (std::make_pair(m1,m2));
             }

          __ATTR_HOT__
	  __ATTR_ALIGN__(16)
	  __ATTR_VECTORCALL__
          static inline MMASK16_RETTYPE
          gms::math::operator<(const ZMM16c4 x,
                      const ZMM16c4 y) {
              const __mmask16 m1 =
                  _mm512_cmp_ps_mask(x.m_re,y.m_re,_CMP_LT_OQ);
              const __mmask16 m2 =
                  _mm512_cmp_ps_mask(x.m_im,y.m_im,_CMP_LT_OQ);
              return (std::make_pair(m1,m2));
             }

          __ATTR_HOT__
	  __ATTR_ALIGN__(16)
	  __ATTR_VECTORCALL__
          static inline MMASK16_RETTYPE
          operator<(const ZMM16c4 x,
                      const std::complex<float> c) {
               const __mmask16 m1 =
                   _mm512_cmp_ps_mask(x.m_re,_mm512_set1_ps(c.real()),_CMP_LT_OQ);
               const __mmask16 m2 =
                   _mm512_cmp_ps_mask(x.m_im,_mm512_set1_ps(c.imag()),_CMP_LT_OQ);
               return (std::make_pair(m1,m2));
              }

          __ATTR_HOT__
	  __ATTR_ALIGN__(16)
	  __ATTR_VECTORCALL__
          static inline MMASK16_RETTYPE
          operator<(const std::complex<float> c,
                               const ZMM16c4 x) {
               const __mmask16 m1 =
                    _mm512_cmp_ps_mask(_mm512_set1_ps(c.real()),x.m_re,_CMP_LT_OQ);
               const __mmask16 m2 =
                    _mm512_cmp_ps_mask(_mm512_set1_ps(c.imag()),x.m_im,_CMP_LT_OQ);
               return (std::make_pair(m1,m2));
               }

          __ATTR_HOT__
	  __ATTR_ALIGN__(16)
	  __ATTR_VECTORCALL__
          static inline MMASK16_RETTYPE
          operator>=(const ZMM16c4 x,
                      const ZMM16c4 y) {
               const __mmask16 m1 =
                    _mm512_cmp_ps_mask(x.m_re,y.m_re,_CMP_GE_OQ);
               const __mmask16 m2 =
                    _mm512_cmp_ps_mask(x.m_im,y.m_im,_CMP_GE_OQ);
                return (std::make_pair(m1,m2));
              }

        __ATTR_HOT__
	__ATTR_ALIGN__(16)
	__ATTR_VECTORCALL__
        static inline MMASK16_RETTYPE
        operator>=(const ZMM16c4 x,
                      const std::complex<float> c) {
            const __mmask16 m1 =
                 _mm512_cmp_ps_mask(x.m_re,_mm512_set1_ps(c.real()),_CMP_GE_OQ);
            const __mmask16 m2 =
                 _mm512_cmp_ps_mask(x.m_im,_mm512_set1_ps(c.imag()),_CMP_GE_OQ);
             return (std::make_pair(m1,m2));
          }

       __ATTR_HOT__
       __ATTR_ALIGN__(16)
       __ATTR_VECTORCALL__
       static inline MMASK16_RETTYPE
       gms::math::operator>=(const std::complex<float> c,
                             const ZMM16c4 x) {
            const __mmask16 m1 =
                  _mm512_cmp_ps_mask(_mm512_set1_ps(c.real()),x.m_re,_CMP_GE_OQ);
            const __mmask16 m2 =
                  _mm512_cmp_ps_mask(_mm512_set1_ps(c.imag()),x.m_im,_CMP_GE_OQ);
            return (std::make_pair(m1,m2));
        }

        __ATTR_HOT__
	__ATTR_ALIGN__(16)
	__ATTR_VECTORCALL__
        static inline MMASK16_RETTYPE
        operator<=(const ZMM16c4 x,
                      const ZMM16c4 y) {
              const __mmask16 m1 =
                 _mm512_cmp_ps_mask(x.m_re,y.m_re,_CMP_LE_OQ);
              const __mmask16 m2 =
                 _mm512_cmp_ps_mask(x.m_im,y.m_im,_CMP_LE_OQ);
              return (std::make_pair(m1,m2));
          }

         __ATTR_HOT__
	 __ATTR_ALIGN__(16)
	 __ATTR_VECTORCALL__
         static inline MMASK16_RETTYPE
         operator<=(const ZMM16c4 x,
                    const std::complex<float> c) {
              const __mmask16 m1 =
                   _mm512_cmp_ps_mask(x.m_re,_mm512_set1_ps(c.real()),_CMP_LE_OQ);
              const __mmask16 m2 =
                   _mm512_cmp_ps_mask(x.m_im,_mm512_set1_ps(c.imag()),_CMP_LE_OQ);
              return (std::make_pair(m1,m2));
         }

        __ATTR_HOT__
	__ATTR_ALIGN__(16)
	__ATTR_VECTORCALL__
        static inline MMASK16_RETTYPE
        operator<=(const std::complex<float> c,
                      const ZMM16c4 x) {
             const __mmask16 m1 =
               _mm512_cmp_ps_mask(_mm512_set1_ps(c.real()),x.m_re,_CMP_LE_OQ);
             const __mmask16 m2 =
               _mm512_cmp_ps_mask(_mm512_set1_ps(c.imag()),x.m_im,_CMP_LE_OQ);
             return (std::make_pair(m1,m2));
         }













     } // math




} // gms









#endif /*__GMS_AVX512C16F32_H__*/
