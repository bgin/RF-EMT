
#ifndef __GMS_SSEC4F32_H__
#define __GMS_SSEC4F32_H__ 310320201538

namespace file_info {

  const unsigned int gGMS_SSEC4F32_MAJOR = 1;
  const unsigned int gGMS_SSEC4F32_MINOR = 0;
  const unsigned int gGMS_SSEC4F32_MICRO = 0;
  const unsigned int gGMS_SSEC4F32_FULLVER =
       1000U*gGMS_SSEC4F32_MAJOR+100U*gGMS_SSEC4F32_MINOR+10U*gGMS_SSEC4F32_MICRO;
  const char * const pgGMS_SSEC4F32_CREATION_DATE = "31-03-2020 15:38 +00200 (TUE 31 MAR 2020 15:38 GMT+2)";
  const char * const pgGMS_SSEC4F32_BUILD_DATE    = __DATE__ ":" __TIME__;
  const char * const pgGMS_SSEC4F32_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
  const char * const pgGMS_SSEC4F32_SYNOPSIS      = "Complex class based on __m128 data type.";

}

#include <cstdint>
#include <immintrin.h>
#include <complex>
#include "GMS_config.h"

namespace gms {

          namespace math {

    
	struct  SSEc4f32 
	{

        __m128 m_re;
		__m128 m_im;
	
	       __ATTR_ALWAYS_INLINE__
		   SSEc4f32() 
		   {
               m_re = _mm_setzero_ps();
			   m_im = _mm_setzero_ps();
		   }

		    __ATTR_ALWAYS_INLINE__
		    SSEc4f32(const float * __restrict  re,
			    const float * __restrict  im) 
		   {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                 re = (const float*)__builtin_assume_aligned(re,16);
				 im = (const float*)__builtin_assume_aligned(im,16);
#elif defined __ICC || defined __INTEL_COMPILER
                __assume_aligned(re,16);
				__assume_aligned(im,16);
#endif
                m_re = _mm_load_ps(&re[0]);
				m_im = _mm_load_ps(&im[0]);
			}

		    __ATTR_ALWAYS_INLINE__
		    SSEc4f32(const float re,
		         const float im) 
		    {
               m_re = _mm_set1_ps(re);
			   m_im = _mm_set1_ps(im);
		    }

		    __ATTR_ALWAYS_INLINE__
		    SSEc4f32(const std::complex<float> c) 
			  {
                m_re = _mm_set1_ps(c.real());
			    m_im = _mm_set1_ps(c.imag());
		    }

		    __ATTR_ALWAYS_INLINE__
		    SSEc4f32(const float re) {
               m_re = _mm_set1_ps(re);
			   m_im = _mm_setzero_ps();
		    }

		    __ATTR_ALWAYS_INLINE__
		      SSEc4f32(const float re) {
                m_re = _mm_set1_ps(re);
			   m_im = _mm_setzero_ps();
		    }

		   
		    __ATTR_ALWAYS_INLINE__
		    SSEc4f32(const float re0,
		             const float re1,
				     const float re2,
				     const float re3)
			{
              m_re = _mm_setr_ps(re3,re2,re1,re0);
			  m_im = _mm_setzero_ps();
		   }

		   __ATTR_ALWAYS_INLINE__
		   SSEc4f32( const float re0,
		             const float re1,
				     const float re2,
				     const float re3,
			         const float im0,
				     const float im1,
				     const float im2,
				     const float im3)
		  {
                m_re = _mm_setr_ps(re3,re2,re1,re0);
			    m_im = _mm_setr_ps(im3,im2,im1,im0);
          
		  }

		    __ATTR_ALWAYS_INLINE__
		    SSEc4f32(const __m128 re,
		             const __m128 im) 
			{
                m_re = re;
			    m_im = im;
            }

		  __ATTR_ALWAYS_INLINE__
		  SSEc4f32(const SSEc4f32 x) 
		  {
                m_re = x.m_re;
                m_im = x.m_im;
		  }

		  __ATTR_ALWAYS_INLINE__
		  SSEc4f32 &
		  load_a(const float * __restrict  re,
		         const float * __restrict  im) 
		  {
#if defined __GNUC__ && !defined __INTEL_COMPILER
             re = (const float*)__builtin_assume_aligned(re,16);
			 im = (const float*)__builtin_assume_aligned(im,16);
#elif defined __ICC || defined __INTEL_COMPILER
            __assume_aligned(re,16);
			 __assume_aligned(im,16);
#endif
               m_re = _mm_load_ps(&re[0]);
		       m_im = _mm_load_ps(&im[0]);
		       return (*this);
		   }

		 __ATTR_ALWAYS_INLINE__
		 SSEc4f32 &
		 load_u(const float * __restrict re,
		        const float * __restrict im) {
                m_re = _mm_loadu_ps(&re[0]);
                m_im = _mm_loadu_ps(&im[0]);
                return (*this);
		 }

		 __ATTR_ALWAYS_INLINE__
		 void store_a(float * __restrict  re,
		              float * __restrict  im) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
             re = (float*)__builtin_assume_aligned(re,16);
			 im = (float*)__builtin_assume_aligned(im,16);
#elif defined __ICC || defined __INTEL_COMPILER
             __assume_aligned(re,16);
			 __assume_aligned(im,16);
#endif
            _mm_store_ps(&re[0],m_re);
            _mm_store_ps(&im[0],m_im);
		}

		
		
		
		void store_u(float * __restrict re,
		             float * __restrict im) {
                      _mm_store_ps(&re[0],m_re);
                      _mm_store_ps(&im[0],m_im);
	        }

	      
	      
	      
	      void stream_nt(float * __restrict re,
	                     float * __restrict im) {
                     _mm_stream_ps(&re[0],m_re);
                     _mm_stream_ps(&im[0],m_im);
                    _mm_sfence();
	        }

	     
	     
	     float extract_1xf32(const int32_t pos) 
		 {
            __attribute__((aligned(16))) float mem[8] = {};
                     store_a(&mem[0],&mem[4]);
                      return (mem[pos & 0x7]);
         }

	     
	     
	     std::pair<float,float>
	     extract_2xf32(const int32_t posx,
	                   const int32_t posy) {
                     __attribute__((aligned(16))) float re_mem[4] = {};
                     __attribute__((aligned(16))) float im_mem[4] = {};
                     store_a(&re_mem[0],&im_mem[0]);
                     return (std::make_pair(re_mem[posx & 0x3],im_mem[posy  & 0x3]));
                }

            
	    
	    SSEc4f32 &
	    insert_1xf32(const int32_t pos,
	                 const float value) {
                    __attribute__((aligned(16))) float mem[8] = {};

                   store_a(&mem[0],&mem[4]);
                   mem[idx & 0x7] = value;
                   load_a(&mem[0],&mem[4]);
                   return (*this);
	        }

	    
	    
	    SSEc4f32 &
	    insert_2xf32(const int32_t re_idx
	                 const int32_t im_idx
			 const float re,
			 const float im) {
                 __attribute__((aligned(16))) float mem_re[4] = {};
                 __attribute__((aligned(16))) float mem_im[4] = {};

                 store_a(&mem_re[0],&mem_im[0]);
                 mem_re[re_idx & 0x3] = re;
                 mem_im[im_idx & 0x3] = im;
                 load_a(&mem_re[0],&mem_im[0]);
                 return (*this);
	       }

	    
            
	    
	    void concatenate_a(float * __restrict  out) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                  out = (float*)__builtin_assume_aligned(out,16);
#elif defined __ICC || defined __INTEL_COMPILER
                  __assume_aligned(out,16);
#endif
                  store_a(&out[0],&out[4]);
	       }

            
	
	    void concatenate_u(float * __restrict out) {
                  store_u(&out[0],&out[4]);
	      }
	      
#if defined __AVX512F__

            
	    
	    SSEc4f32 &
	    partial_loadu(const float * __restrict re,
	                  const int32_t n_re,
			  const float * __restrict im,
			  const int32_t n_im) {
                   m_re = _mm_maskz_loadu_ps(__mmask8((1 << n_re)-1),re);
                   m_im = _mm_maskz_loadu_ps(__mmask8((1 << n_im)-1),im);
                   return (*this);
	      }

	    
	   
	    SSEc4f32 &
	    partial_loada(const float * __restrict re,
	                  const int32_t n_re,
			  const float * __restrict im,
			  const int32_t n_im) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                   re = (const float*)__builtin_assume_aligned(re,16);
		   im = (const float*)__builtin_assume_aligned(im,16);
#elif defined __ICC || defined __INTEL_COMPILER
                   __assume_aligned(re,16);
		   __assume_aligned(im,16);
#endif
                    m_re = _mm_maskz_load_ps(__mmask8((1 << n_re)-1),re);
                    m_im = _mm_maskz_load_ps(__mmask8((1 << n_im)-1),im);
                    return (*this);
	      }

	      
	      
	      void partial_storeu(float * __restrict re,
                                       const int32_t n_re,
				       float * __restrict im,
				       const int32_t n_im) {
                 _mm_mask_storeu_ps(&re[0],__mmask8((1 << n_re)-1),m_re);
                 _mm_mask_storeu_ps(&im[0],__mmask8((1 << m_im)-1),m_im);
              }

	      
	      
	      
	      SSEc4f32 &
	             expand(const SSEc4f32 x,
                            const __mmask8 mask) {
                     m_re = _mm_maskz_expand_ps(mask,x.m_re);
                     m_im = _mm_maskz_expand_ps(mask.x.m_im);
                     return (*this);
              }

	      
	     
	      
	      SSEc4f32 &
	            expand_loadu(const SSEc4f32 x,
                                 const __mmask8 mask,
				 const float * __restrict re,
				 const float * __restrict im) {
                     m_re = _mm_mask_expandloadu_ps(x.m_re,mask,&re[0]);
                     m_im = _mm_mask_expandloadu_ps(x.m_im,mask,&im[0]);
                     return (*this);
               }

	       
	       
	       SSEc4f32 & permute(const __mmask8 mask,
                                      const int32_t imm) {
                    m_re = _mm_mask_permute_ps(m_re,mask,m_im,imm);
                    m_im = _mm_mask_permute_ps(m_im,mask,m_re,imm);
                    return (*this);
               }
#endif
	       
	      
	       
	       SSEc4f32 &
	              operator=(const SSEc4f32 x) {
                        if(this == &x) return (*this);
                        m_re = x.m_re;
                        m_im = x.m_im;
                        return (*this);
	          }

		 
		     
	     }__ATTR_ALIGN__(16);



	     
	     
	     
	     static inline
	     SSEc4f32 conj(SSEc4f32 x) {
	          auto tmp = ~x;
		  return (SSEc4f32{x.m_re,tmp.m_im});
	     }

	     __ATTR__HOT__
	     
	     
	     static inline
	     SSEc4f32 polar(const __m128 rho,
	                    const __m128 theta) {
                   const __m128 re_part =
                        _mm_mul_ps(rho,_mm_cos_ps(theta));
                   const __m128 im_part =
                        _mm_mul_ps(rho,_mm_sin_ps(theta));
                   return (SSEc4f32{re_part,im_part});
	     }

	     
	     
	     
	     static inline
	     __m128 carg(const SSEc4f32 x) {
                  return (_mm_atan2_ps(x.m_re,x.m_im));
	     }

	     
	     
	     
	     static inline
	     SSEc4f32 csin(const SSEc4f32 x) {
                    const __m128 re_part =
                             _mm_mul_ps(_mm_sin_ps(x.m_re),
	                                   _mm_cosh_ps(x.m_im));
                    const __m128 im_part =
                             _mm_mul_ps(_mm_cos_ps(x.m_re),
	                                   _mm_sinh_ps(x.m_im));
                    return (SSEc4f32{re_part,im_part});
	     }

	     
	     
	     
	     static inline
	     SSEc4f32 csinh(const SSEc4f32 x) {
                    const __m128 re_part =
                             _mm_mul_ps(_mm_sinh_ps(x.m_re),
	                                   _mm_cos_ps(x.m_im));
                    const __m128 im_part =
                             _mm_mul_ps(_mm_cosh_ps(x.m_re),
	                                   _mm_sin_ps(x.m_im));
                    return (SSEc4f32{re_part,im_part});
	     }

	     
	     
	     
	     static inline
	     SSEc4f32 ccos(const SSEc4f32 x) {
                    const __m128 re_part =
                                _mm_mul_ps(_mm_cos_ps(x.m_re),
	                                      _mm_cosh_ps(x.m_im));
                    const __m128 im_part =
                                _mm_mul_part(_mm_sin_ps(x.m_re),
	                                      _mm_sinh_ps(x.m_im));
                    return (SSEc4f32{re_part,im_part});
	     }

	     
	     
	     
	     static inline
	     SSEc4f32 ccosh(const SSEc4f32 x) {
                    const __m128 re_part =
                                _mm_mul_ps(_mm_cosh_ps(x.m_re),
	                                      _mm_cos_ps(x.m_im));
                    const __m128 im_part =
                                _mm_mul_part(_mm_sinh_ps(x.m_re),
	                                      _mm_sin_ps(x.m_im));
                    return (SSEc4f32{re_part,im_part});
	     }

	     
	     
	     
	     static inline
	     SSEc4f32 cexp(const SSEc4f32 x) {
                    const __m128 re_part =
                                 _mm_mul_ps(_mm_exp_ps(x.m_re),
	                                       _mm_cos_ps(x.m_im));
                    const __m128 im_part =
                                 _mm_mul_ps(_mm_exp_ps(x.m_re),
	                                       _mm_sin_ps(x.m_im));
                    return (SSEc4f32{re_part,im_part});
	     }

	     
	     
	     
	     static inline
	     __m128 cabs(const SSEc4f32 x) {
                    const __m128 re_part =
                             _mm_mul_ps(x.m_re,x.m_re);
                    const __m128 im_part =
                             _mm_mul_ps(x.m_im,x.m_im);
                    return (_mm_sqrt_ps(_mm_add_ps(re_part,im_part)));
             }

	     
	     
	     
	     static inline
	     SSEc4f32 cpow(const SSEc4f32 x,
                               const float n) {
                    const __m128 re_part =
                               _mm_mul_ps(x.m_re,x.m_re);
                    const __m128 im_part =
                               _mm_mul_ps(x.m_im,x.m_im);
                    const __m128 r =
                               _mm_sqrt_ps(_mm_add_ps(re_part,im_part));
                    const __m128 theta =
                               _mm_atan_ps(_mm_div_ps(x.m_im,x.m_re));
                    const __m128 vn = _mm128_set1_ps(n);
                    const __m128 pow_term = _mm_pow_ps(r,vn);
                    const __m128 trig_arg = _mm_mul_ps(vn,theta);
                    return (SSEc4f32{_mm_mul_ps(pow_term,_mm_cos_ps(trig_arg),
                                     _mm_mul_ps(pow_term,_mm_sin_ps(trig_arg))}));
             }

	     
	     
	     
	     static inline
	     SSEc4f32 clog(const SSEc4f32 x) {
                     const __m128 t1  = cabs(x);
                     const __m128 t2  = carg(x);
                     const __m128 re_part = _mm_log_ps(t1);
                     return (SSEc4f32{re_part,t2});
	     }

	     
	     
	     
	     static inline
	     SSEc4f32 csqrt(const SSEc4f32 x) {
                       const __m128 t = cabs(x);
                       const __m128 re_part =
                               _mm_mul_ps(_mm_set1_ps(0.5f),
	                                     _mm_add_ps(t,x.m_re));
                       const __m128 im_part =
                               _mm_mul_ps(_mm_set1_ps(0.5f),
	                                     _mm_sub_ps(t,x.m_re));
                       return (SSEc4f32{_mm_sqrt_ps(re_part),
                                        _mm_sqrt_ps(im_part)});
	     }

	     
	     
	     
	     static inline
	     SSEc4f32 ctan(const SSEc4f32 x) {
                       return (ctan(x)/csin(x));
	     }

	     
	     
	     
             static inline
	     SSEc4f32
             ctanh(const SSEc4f32 x) {
                        return (csinh(x)/ccosh(x));
             }
#if defined __AVX512F__
	     
	     
	     
             static inline
	     SSEc4f32
             select(const SSEc4f32 x,
                    const SSEc4f32 y,
		    __mmask8 mask) {
                 return (SSEc4f32{_mm_mask_blend_ps(mask,x.m_re,y.m_re),
                                     _mm_mask_blend_ps(mask,x.m_im,y.m_im)});
              }
#endif

              
              
              
              static inline
	      __m128 abs_real(const __m128 x) {

			 static const __m128 mask = _mm_set1_pd(0x7FFFFFFF);
			 return (_mm_and_pd(x,mask));
	       }

             static inline
	         SSEc4f32
             cdiv_smith(const SSEc4f32 x,
                        const SSEc4f32 y) {
                    __m128 ratio,denom,re_part,im_part;
                    constexpr __mmask8 all_ones = 0xFF;
                    re_part = _mm_setzero_ps();
                    im_part = _mm_setzero_ps();
                    __mmask8 is_gte = _mm_cmp_ps_mask(
                                abs_real(y.m_re),
			        abs_real(y.m_im),_CMP_GE_OQ);
                    ratio = _mm_setzero_ps();
                    denom = _mm_setzero_ps();
                    if(is_gte == all_ones) {
                        ratio = _mm_div_ps(y.m_im,y.m_re);
	                denom = _mm_add_ps(y.m_re,
	                         _mm_mul_ps(ratio,y.m_im));
	                re_part = _mm_div_ps(_mm_add_ps(x.m_re,
	                          _mm_mul_ps(x.m_im,ratio)),denom);
	                im_part = _mm_div_ps(_mm_sub_ps(x.m_im,
	                          _mm_mul_ps(x.m_re,ratio)),denom);
	                return (SSEc4f32{re_part,im_part});
                     }
                     else {
                        ratio = _mm_div_ps(y.m_re,y.m_im);
	                denom = _mm_add_ps(y.m_im,_mm_mul_ps(ratio,y.m_re));
	                re_part = _mm_div_ps(_mm_add_ps(
	                          _mm_mul_ps(x.m_re,ratio)),denom);
	                im_part = _mm_div_ps(_mm_sub_ps(
	                          _mm_mul_ps(x.m_im,ratio)),denom);
	                return (SSEc4f32{re_part,im_part});
                      }
               }

                  
	      
	      
	      
              static inline
	      SSEc4f32
              operator+(const SSEc4f32 x,
                        const SSEc4f32 y) {
                  return (SSEc4f32{_mm_add_ps(x.m_re,y.m_re),
                                       _mm_add_ps(x.m_im,x.m_im)});
               }

	       
	      
	       
	       static inline
	       SSEc4f32
	       operator+(const SSEc4f32 x,
	                 const std::complex<float> y) {
                   return(SSEc4f32{_mm_add_ps(x.m_re,_mm_set1_ps(y.real())),
		                       _mm_add_ps(x.m_im,_mm_set1_ps(y.imag()))});
	       }

	       
	       
	       
	       static inline
	       SSEc4f32
	       operator+(const std::complex<float> x,
	                 const SSEc4f32 y) {
                   return (SSEc4f32{_mm_add_ps(_mm_set1_ps(x.real()),y.m_re),
		                        _mm_add_ps(_mm_set1_ps(x.imag()),y.m_im)});
		}

	      
	      
	      
              static inline SSEc4f32
              operator+(const SSEc4f32 x,
                        const __m128 v) {
                  return (SSEc4f32{_mm_add_ps(x.m_re,v),
                                              x.m_im});
               }

              
	      
	      
              static inline SSEc4f32
              operator+(const __m128 v,
                        const SSEc4f32 x) {
                   return (SSEc4f32{_mm_add_ps(v,x.m_re),x.m_im});
               }

              
	      
	      
              static inline SSEc4f32
              operator+(const SSEc4f32 x,
                        const float s) {
                   return (x + SSEc4f32{s});
               }

              
	      
	      
              static inline SSEc4f32
              operator+(const float s,
                       const SSEc4f32 x) {
                   return (SSEc4f32{s} + x);
               }

	     
	     
	     
             static inline SSEc4f32
             operator+=(SSEc4f32 x,
                        const SSEc4f32 y) {
                  x = x + y;
                  return (x)
              }

	     
	    
	     
	     static inline SSEc4f32
	     operator+=(SSEc4f32 x,
	                const std::complex<float> y) {
                 x = x + y;
		 return (x);
	     }

	     
	     
	     
	     static inline SSEc4f32
	     operator+=(const std::complex<float> x,
	                SSEc4f32 y) {
                 y = y + x;
		 return (y);
	     }

             
	     
	     
             static inline SSEc4f32
             operator+=(SSEc4f32 x,
                        const __m128 v) {
                 x = x + v;
                 return (x);
              }

             
	     
	     
             static inline SSEc4f32
             operator+=(const __m128 v,
                        SSEc4f32 x) {
                x = v + x;
                return (x);
              }

             
	     
	     
             static inline SSEc4f32
             operator+=(SSEc4f32 x,
                        const float s) {
                  x = x + SSEc4f32{s};
                  return (x)
              }

             
	     
	     
             static inline SSEc4f32
             operator+=(const float s,
                        SSEc4f32 x) {
                x = SSEc4f32{s} + x;
                return (x);
              }

	     
	     
	     
             static inline SSEc4f32
             operator-(const SSEc4f32 x,
                     const SSEc4f32 y) {
                  return (SSEc4f32{_mm_sub_ps(x.m_re,y.m_re),
                                      _mm_sub_ps(x.m_im,y.m_im)});
              }

	    
	    
	    
	    static inline SSEc4f32
	    operator-(const SSEc4f32 x,
	              const std::complex<float> y) {
                 return(SSEc4f32{_mm_sub_ps(x.m_re,_mm_set1_ps(y.real())),
		                     _mm_sub_ps(x.m_im,_mm_set1_ps(y.imag()))});
	    }

	    
	    
	    
	    static inline SSEc4f32
	    operator-(const std::complex<float> x,
	              const SSEc4f32 y) {
                  return (SSEc4f32{_mm_sub_ps(_mm_set1_ps(x.real()),y.m_re),
		                     _mm_sub_ps(_mm_set1_ps(x.imag()),y.m_im)});
	    }
	    

            
	    
	    
            static inline SSEc4f32
            operator-(const SSEc4f32 x,
                      const __m128 v) {
                return (SSEc4f32{_mm_sub_ps(x.m_re,v),x.m_im});
             }

            
	    
	    
            static inline SSEc4f32
            operator-(const __m128 v,
                     const SSEc4f32 x) {
                  return (SSEc4f32{_mm_sub_ps(v,x.m_re),x.m_im});
             }

            
	    
	    
            static inline SSEc4f32
            operator-(const SSEc4f32 x,
                     const float s) {
                 return (x - SSEc4f32{s});
             }

	    
	    
	    
            static inline SSEc4f32
            operator-(const float s,
                     const SSEc4f32 x) {
                  return (SSEc4f32{s} - x);
              }

	    
	    
	    
	    static inline SSEc4f32
	    operator-=(SSEc4f32 x,
	               const SSEc4f32 y) {
                x = x - y;
		return (x);
	     }

	   
	   
	   
	   static inline SSEc4f32
	   operator-=(SSEc4f32 x,
	              const std::complex<float> y) {
                x = x - y;
		return (x);
	    }

	   
	   
	   
	   static inline SSEc4f32
	   operator-=(const std::complex<float> x,
	              SSEc4f32 y) {
               y = x - y;
	       return (y);
	    }

	   
	   
	   
	   static inline SSEc4f32
	   operator-=(SSEc4f32 x,
	              const __m128 y) {
              x = x - y;
	      return (x);
	   }

	  
	  
	  
	  static inline SSEc4f32
	  operator-=(const __m128 x,
	             SSEc4f32 y) {
              y = x - y;
	      return (y);
	  }

	 
	 
	 
	 static inline SSEc4f32
	 operator-=(SSEc4f32 x,
	            const float s) {
            x = x - SSEc4f32{s};
	    return (x);
	 }

	 
	 
	 
	 static inline SSEc4f32
	 operator-=(const float s,
	            SSEc4f32 x) {
            x = SSEc4f32{s}-x;
	    return (x);
	 }

	
	
	
        static inline SSEc4f32
        operator*(const SSEc4f32 x,
                  const SSEc4f32 y) {
              const __m128 xmm0 = _mm_mul_ps(x.m_re,y.m_re);
              const __m128 xmm1 = _mm_mul_ps(x.m_im,y.m_im);
              const __m128 xmm2 = _mm_mul_ps(x.m_im,y.m_re);
              const __m128 xmm3 = _mm_mul_ps(x.m_re,y.m_im);
              return (SSEc4f32{_mm_add_ps(xmm0,xmm1),
                                   _mm_sub_ps(xmm2,xmm3)});
           }

	 
	 
	 
	 static inline SSEc4f32
	 operator*(const SSEc4f32 x,
	           const std::complex<float> y) {
            const __m128 y_re = _mm_set1_ps(y.real());
	    const __m128 y_im = _mm_set1_ps(y.imag());
	    const __m128 zmm0 = _mm_mul_ps(x.m_re,y_re);
	    const __m128 zmm1 = _mm_mul_ps(x.m_im,y_im);
	    const __m128 zmm2 = _mm_mul_ps(x.m_im,y_re);
	    const __m128 zmm3 = _mm_mul_ps(x.m_re,y_im);
	    return (SSEc4f32{_mm_add_ps(zmm0,zmm1),
                                 _mm_sub_ps(zmm2,zmm3)});
	 }

        
	
	
        static inline SSEc4f32
        operator*(const SSEc4f32 x,
                     const __m128 v) {
              return (SSEc4f32{_mm_mul_ps(x.m_re,v),
                                  _mm_mul_ps(x.m_im,v)});
          }

       
       
       
       static inline SSEc4f32
       operator*(const __m128 v,
                 const SSEc4f32 x) {
            return (SSEc4f32{_mm_mul_ps(v,x.m_re),
                                _mm_mul_ps(v,x.m_im)});
          }


       
       
       
       static inline SSEc4f32
       operator*(const SSEc4f32 x,
                 const float s) {
         const __m128 zmm0 = _mm_set1_ps(s);
              return (SSEc4f32{_mm_mul_ps(x.m_re,zmm0),
                                  _mm_mul_ps(x.m_im,zmm0)});
          }

       
       
       
       static inline SSEc4f32
       operator*=(SSEc4f32 x,
                  const SSEc4f32 y) {
              x = x * y;
              return (x);
         }

       
       
       
       static inline SSEc4f32
       operator*=(SSEc4f32 x,
                  const __m128 v) {
            x = x * v;
            return (x);
          }

       
       
       
       static inline SSEc4f32
       operator*=(const __m128 v,
                   SSEc4f32 x) {
              x = v * x;
              return (x);
          }

        
	
	
        static inline SSEc4f32
        operator*=(SSEc4f32 x,
                   const float s) {
             x = x * s;
             return (x);
          }

        
	
	
         static inline SSEc4f32
         operator*=(const float s,
                    SSEc4f32 x) {
             x = s * x;
             return (x);
          }

	
	
	
        static inline SSEc4f32
        operator/(const SSEc4f32 x,
                  const SSEc4f32 y) {
#if defined USE_SAFE_COMPLEX_DIVISION && (USE_SAFE_COMPLEX_DIVISION) == 1
              return (cdiv_smith(x,y));
#else
              const __m128 zmm0 = _mm_mul_ps(x.m_re,y.m_re);
              const __m128 zmm1 = _mm_mul_ps(x.m_im,y.m_im);
              const __m128 zmm2 = _mm_mul_ps(x.m_im,y.m_re);
              const __m128 zmm3 = _mm_mul_ps(x.m_re,y.m_im);
              const __m128 den  = _mm_add_ps(_mm_mul_ps(y.m_re,y.m_re),
                                                _mm_mul_ps(y.m_im,y.m_im));
              const __m128 re_part = _mm_add_ps(zmm0,zmm1);
              const __m128 im_part = _mm_sub_ps(zmm2,zmm3);
              return (SSEc4f32{_mm_div_ps(re_part,den),
                                  _mm_div_ps(im_part,den)});
#endif
           }

	 
	 
	 
	 static inline SSEc4f32
	 operator/(const SSEc4f32 x,
	           const std::complex<float> y) {

              const __m128 y_re = _mm_set1_ps(y.real());
	      const __m128 y_im = _mm_set1_ps(y.imag());
	      const __m128 zmm0 = _mm_mul_ps(x.m_re,y_re);
              const __m128 zmm1 = _mm_mul_ps(x.m_im,y_im);
              const __m128 zmm2 = _mm_mul_ps(x.m_im,y_re);
              const __m128 zmm3 = _mm_mul_ps(x.m_re,y_im);
              const __m128 den  = _mm_add_ps(_mm_mul_ps(y_re,y_re),
                                                _mm_mul_ps(y_im,y_im));
              const __m128 re_part = _mm_add_ps(zmm0,zmm1);
              const __m128 im_part = _mm_sub_ps(zmm2,zmm3);
              return (SSEc4f32{_mm_div_ps(re_part,den),
                                  _mm_div_ps(im_part,den)});
	  }


	  	  
	  
	  
	  static inline SSEc4f32
	  operator/(const std::complex<float> x,
	            const SSEc4f32 y) {

	      const __m128 x_re = _mm_set1_ps(x.real());
	      const __m128 x_im = _mm_set1_ps(x.imag());
	      const __m128 zmm0 = _mm_mul_ps(x_re,y.m_re);
              const __m128 zmm1 = _mm_mul_ps(x_im,y.m_im);
              const __m128 zmm2 = _mm_mul_ps(x_im,y.m_re);
              const __m128 zmm3 = _mm_mul_ps(x_re,y.m_im);
              const __m128 den  = _mm_add_ps(_mm_mul_ps(y.m_re,y.m_re),
                                                _mm_mul_ps(y.m_im,y.m_im));
              const __m128 re_part = _mm_add_ps(zmm0,zmm1);
              const __m128 im_part = _mm_sub_ps(zmm2,zmm3);
              return (SSEc4f32{_mm_div_ps(re_part,den),
                                  _mm_div_ps(im_part,den)});
	  }

	  
	  
	  
          static inline SSEc4f32
          operator/(const SSEc4f32 x,
                     const __m128 v) {
                return (SSEc4f32{_mm_div_ps(x.m_re,v),
                                    _mm_div_ps(x.m_im,v)});
           }

         
	 
	 
         static inline SSEc4f32
         operator/(const __m128 v,
                   const SSEc4f32 x) {
                return (SSEc4f32{_mm_div_ps(v,x.m_re),
                                     _mm_div_ps(v,x.m_im)});
           }

         
	 
	 
        static inline SSEc4f32
        operator/(const SSEc4f32 x,
                     const float s) {
               const __m128 zmm0 = _mm_set1_ps(s);
               return (SSEc4f32{_mm_div_ps(x.m_re,zmm0),
                                    _mm_div_ps(x.m_im,zmm0)});
           }

         
	 
	 
        static inline SSEc4f32
        operator/(const float s,
                  const SSEc4f32 x) {
               const __m128 zmm0 = _mm_set1_ps(s);
               return (SSEc4f32{_mm_div_ps(zmm0,x.m_re),
                                    _mm_div_ps(zmm0,x.m_im)});
           }

	
	
	
        static inline SSEc4f32
        operator/=(SSEc4f32 x,
                   const SSEc4f32 y) {
              x = x / y;
              return (x);
          }

	
	
	
	static inline AVX512c16f32
	operator/=(AVX512c16f32 x,
	           const std::complex<float> y) {
              x = x / y
	      return (x);
	 }

	 
	 
	 
	 static inline SSEc4f32
	 operator/=(const std::complex<float> x,
	            SSEc4f32 y) {
               y = x / y;
               return (y);
	 }

        
	
	
        static inline SSEc4f32
        operator/=(SSEc4f32 x,
                      const __m128 v) {
              x = x / v;
              return (x);
          }

         
	 
	 
         static inline SSEc4f32
         operator/=(const __m128 v,
                    SSEc4f32 x) {
             x = v / x;
             return (x);
          }

         
	 
	 
         static inline SSEc4f32
         operator/=(SSEc4f32 x,
                    const float s) {
              x = x / s;
              return (x);
          }

	 
	 
	 
         static inline SSEc4f32
         operator/=(const float s,
                    SSEc4f32 x) {
              x = s / x;
              return (x);
          }

          
	  
	  
          static inline SSEc4f32
          operator~(SSEc4f32 x) {
                x.m_re = _mm_sub_ps(_mm_setzero_ps(),x.m_re);
                return (x);
           }

#if defined __AVX512F__


	 
	 
	 
	 static inline std::pair<__mmask8,__mmask8>
         operator==(const SSEc4f32 x,
                    const SSEc4f32 y) {
               const __mmask8 m1 =
                     _mm_cmp_ps_mask(x.m_re,y.m_re,_CMP_EQ_OQ);
               const __mmask8 m2 =
                     _mm_cmp_ps_mask(x.m_im,y.m_im,_CMP_EQ_OQ);
               return (std::make_pair(m1,m2));
           }

	 
	 
	 
	 static inline std::pair<__mmask8,__mmask8>
         operator==(const SSEc4f32 x,
                    const std::complex<float> c) {
               const __mmask8 m1 =
                  _mm_cmp_ps_mask(x.m_re,_mm_set1_ps(c.real()),_CMP_EQ_OQ);
               const __mmask8 m2 =
                   _mm_cmp_ps_mask(x.m_im,_mm_set1_ps(c.imag()),_CMP_EQ_OQ);
               return (std::make_pair(m1,m2));
           }

         
	 
	 
	 static inline std::pair<__mmask8,__mmask8>
         operator==(const std::complex<float> c,
                    const SSEc4f32 x) {
              const __mmask8 m1 =
                    _mm_cmp_ps_mask(_mm_set1_ps(c.real()),x.m_re,_CMP_EQ_OQ);
              const __mmask8 m2 =
                    _mm_cmp_ps_mask(_mm_set1_ps(c.imag()),x.m_im,_CMP_EQ_OQ);
              return (std::make_pair(m1,m2));
           }

          
	  
	  
	  static inline std::pair<__mmask8,__mmask8>
          operator!=( const SSEc4f32 x,
                      const SSEc4f32 y) {
             const __mmask8 m1 =
                   _mm_cmp_ps_mask(x.m_re,y.m_re,_CMP_NEQ_OQ);
             const __mmask8 m2 =
                   _mm_cmp_ps_mask(x.m_im,y.m_im,_CMP_NEQ_OQ);
             return (std::make_pair(m1,m2));
           }

          
	  
	  
	  static inline std::pair<__mmask8,__mmask8>
          operator!=(const SSEc4f32 x,
                      const std::complex<float> c) {
             const __mmask8 m1 =
                    _mm_cmp_ps_mask(x.m_re,_mm_set1_ps(c.real()),_CMP_NEQ_OQ);
             const __mmask8 m2 =
                    _mm_cmp_ps_mask(x.m_im,_mm_set1_ps(c.imag()),_CMP_NEQ_OQ);
             return (std::make_pair(m1,m2));
            }

           
	   
	   
	   static inline std::pair<__mmask8,__mmask8>
           operator!=(const std::complex<float> c,
                      const SSEc4f32 x) {
                const __mmask8 m1 =
                   _mm_cmp_ps_mask(_mm_set1_ps(c.real()),x.m_re,_CMP_NEQ_OQ);
                const __mmask8 m2 =
                   _mm_cmp_ps_mask(_mm_set1_ps(c.imag()),x.m_im,_CMP_NEQ_OQ);
                return (std::make_pair(m1,m2));
             }

           
	   
	   
	   static inline std::pair<__mmask8,__mmask8>
           operator>(const SSEc4f32 x,
                      const SSEc4f32 y) {
              const __mmask8 m1 =
                  _mm_cmp_ps_mask(x.m_re,y.m_re,_CMP_GT_OQ);
              const __mmask8 m2 =
                  _mm_cmp_ps_mask(x.m_im,y.m_im,_CMP_GT_OQ);
               return (std::make_pair(m1,m2));
             }

          
	  
	  
	  static inline std::pair<__mmask8,__mmask8>
          operator>(const SSEc4f32 x,
                      const std::complex<float> c) {
              const __mmask8 m1 =
                   _mm_cmp_ps_mask(x.m_re,_mm_set1_ps(c.real()),_CMP_GT_OQ);
              const __mmask8 m2 =
                   _mm_cmp_ps_mask(x.m_im,_mm_set1_ps(c.imag()),_CMP_GT_OQ);
              return (std::make_pair(m1,m2));
             }

           
	   
	   
	   static inline std::pair<__mmask8,__mmask8>
           operator>(const std::complex<float> c,
                      const SSEc4f32 x) {
              const __mmask8 m1 =
                 _mm_cmp_ps_mask(_mm_set1_ps(c.real()),x.m_re,_CMP_GT_OQ);
              const __mmask8 m2 =
                 _mm_cmp_ps_mask(_mm_set1_ps(c.imag()),x.m_im,_CMP_GT_OQ);
              return (std::make_pair(m1,m2));
             }

          
	  
	  
	  static inline std::pair<__mmask8,__mmask8>
          operator<(const SSEc4f32 x,
                      const SSEc4f32 y) {
              const __mmask8 m1 =
                  _mm_cmp_ps_mask(x.m_re,y.m_re,_CMP_LT_OQ);
              const __mmask8 m2 =
                  _mm_cmp_ps_mask(x.m_im,y.m_im,_CMP_LT_OQ);
              return (std::make_pair(m1,m2));
             }

          
	  
	  
	    static inline std::pair<__mmask8,__mmask8>
          operator<(const AVX512c16f32 x,
                      const std::complex<float> c) {
               const __mmask8 m1 =
                   _mm_cmp_ps_mask(x.m_re,_mm_set1_ps(c.real()),_CMP_LT_OQ);
               const __mmask8 m2 =
                   _mm_cmp_ps_mask(x.m_im,_mm_set1_ps(c.imag()),_CMP_LT_OQ);
               return (std::make_pair(m1,m2));
              }

          
	  
	  
	    static inline std::pair<__mmask8,__mmask8>
          operator<(const std::complex<float> c,
                    const SSEc4f32 x) {
               const __mmask8 m1 =
                    _mm_cmp_ps_mask(_mm_set1_ps(c.real()),x.m_re,_CMP_LT_OQ);
               const __mmask8 m2 =
                    _mm_cmp_ps_mask(_mm_set1_ps(c.imag()),x.m_im,_CMP_LT_OQ);
               return (std::make_pair(m1,m2));
               }

          
	  
	  
	  static inline std::pair<__mmask8,__mmask8>
          operator>=(const SSEc4f32 x,
                      const SSEc4f32 y) {
               const __mmask8 m1 =
                    _mm_cmp_ps_mask(x.m_re,y.m_re,_CMP_GE_OQ);
               const __mmask8 m2 =
                    _mm_cmp_ps_mask(x.m_im,y.m_im,_CMP_GE_OQ);
                return (std::make_pair(m1,m2));
              }

        
	
	
	static inline std::pair<__mmask8,__mmask8>
        operator>=(const SSEc4f32 x,
                   const std::complex<float> c) {
            const __mmask8 m1 =
                 _mm_cmp_ps_mask(x.m_re,_mm_set1_ps(c.real()),_CMP_GE_OQ);
            const __mmask8 m2 =
                 _mm_cmp_ps_mask(x.m_im,_mm_set1_ps(c.imag()),_CMP_GE_OQ);
             return (std::make_pair(m1,m2));
          }

       
       
       
       static inline std::pair<__mmask8,__mmask8>
       operator>=(const std::complex<float> c,
                   const SSEc4f32 x) {
            const __mmask8 m1 =
                  _mm_cmp_ps_mask(_mm_set1_ps(c.real()),x.m_re,_CMP_GE_OQ);
            const __mmask8 m2 =
                  _mm_cmp_ps_mask(_mm_set1_ps(c.imag()),x.m_im,_CMP_GE_OQ);
            return (std::make_pair(m1,m2));
        }

        
	
	
	static inline std::pair<__mmask8,__mmask8>
        operator<=(const SSEc4f32 x,
                   const SSEc4f32 y) {
              const __mmask8 m1 =
                 _mm_cmp_ps_mask(x.m_re,y.m_re,_CMP_LE_OQ);
              const __mmask8 m2 =
                 _mm_cmp_ps_mask(x.m_im,y.m_im,_CMP_LE_OQ);
              return (std::make_pair(m1,m2));
          }

         
	 
	 
	 static inline std::pair<__mmask8,__mmask8>
         operator<=(const SSEc4f32 x,
                    const std::complex<float> c) {
              const __mmask8 m1 =
                   _mm_cmp_ps_mask(x.m_re,_mm_set1_ps(c.real()),_CMP_LE_OQ);
              const __mmask16 m2 =
                   _mm_cmp_ps_mask(x.m_im,_mm_set1_ps(c.imag()),_CMP_LE_OQ);
              return (std::make_pair(m1,m2));
         }

        
	
	
	static inline std::pair<__mmask8,__mmask8>
        operator<=(const std::complex<float> c,
                      const SSEc4f32 x) {
             const __mmask8 m1 =
               _mm_cmp_ps_mask(_mm_set1_ps(c.real()),x.m_re,_CMP_LE_OQ);
             const __mmask8 m2 =
               _mm_cmp_ps_mask(_mm_set1_ps(c.imag()),x.m_im,_CMP_LE_OQ);
             return (std::make_pair(m1,m2));
         }


#else

	  
	  
	  
          static inline
	  std::pair<int32_t,int32_t>
	  operator==(const SSEc4f32 x,
		     const SSEc4f32 y) {

	      const __m128 eqre,eqim;
	      eqre = _mm_setzero_ps();
	      eqim = _mm_setzero_ps();
	      eqre = _mm_cmp_ps(x.m_re,y.m_re,
					 _CMP_EQ_OQ);
	      eqim = _mm_cmp_ps(x.m_im,y.m_im,
					 _CMP_EQ_OQ);
              return (std::make_pair(_mm_testz_ps(eqre,eqre),
				     _mm_testz_ps(eqim,eqim)));
 	     }
	  
           
	   
	   
           static inline
	   std::pair<int32_t,int32_t>
	   operator==(const SSEc4f32 x,
		      const std::complex<float> y) {

                   const __m128 eqre,eqim;
		   eqre = _mm_setzero_ps();
		   eqim = _mm_setzero_ps();
		   eqre = _mm_cmp_ps(x.m_re,
				    _mm_set1_ps(y.real()),
						_CMP_EQ_OQ);
		   eqim = _mm_cmp_ps(x.m_im,
				     _mm_set1_ps(y.imag()),
						_CMP_EQ_OQ);
		   return (std::make_pair(_mm_testz_ps(eqre,eqre),
					 _mm_testz_ps(eqim,eqim))); 
	       }
	   
           
	   
	   
           static inline
	   std::pair<int32_t,int32_t>
	   operator==(const std::complex<float> x,
		      const SSEc4f32 y) {

                  const __m128 eqre,eqim;
		  eqre = _mm_setzero_ps();
		  eqim = _mm_setzero_ps();
		  eqre = _mm_cmp_ps(_mm_set1_ps(x.real()),
						y.m_re,
						_CMP_EQ_OQ);
		  eqim = _mm_cmp_ps(_mm_set1_ps(x.imag()),
						 y.m_im,
						 _CMP_EQ_OQ);
		  return (std::make_pair(_mm_testz_ps(eqre,eqre),
					 _mm_testz_ps(eqim,eqim))); 
		}
	   
            
	    
	    
            static inline
	    std::pair<int32_t,int32_t>
	    operator!=(const SSEc4f32 x,
		       const SSEc4f32 y) {

		   const __m128 eqre,eqim;
		   eqre = _mm_setzero_ps();
		   eqim = _mm_setzero_ps();
		   eqre = _mm_cmp_ps(x.m_re,y.m_re,
				        _CMP_NEQ_OQ);
		   eqim = _mm_cmp_ps(x.m_im,y.m_im,
				        _CMP_NEQ_OQ);
                   return (std::make_pair(_mm_testz_ps(eqre,eqre),
					  _mm_testz_ps(eqim,eqim)));
 		}
	    
             
	     
	     
             static inline
	     std::pair<int32_t,int32_t>
	     operator!=(const SSEc4f32 x,
			const std::complex<float> y) {

                     const __m128 eqre,eqim;
		     eqre = _mm_setzero_ps();
		     eqim = _mm_setzero_ps();
		     eqre = _mm_cmp_ps(x.m_re,
				       _mm_set1_ps(y.real()),
				       _CMP_NEQ_OQ);
		     eqim = _mm_cmp_ps(x.m_im,
				       _mm_set1_ps(y.imag()),
				       _CMP_NEQ_OQ);
		     return (std::make_pair(_mm_testz_ps(eqre,eqre),
					    _mm_testz_ps(eqim,eqim))); 
		 }
	     
               
	       
	       
               static inline
	       std::pair<int32_t,int32_t>
	       operator!=(const std::complex<float> x,
			  const SSEc4f32 y) {

                     const __m128 eqre,eqim;
		     eqre = _mm_setzero_ps();
		     eqim = _mm_setzero_ps();
		     eqre = _mm_cmp_ps(_mm_set1_ps(x.real()),
						   y.m_re,
						   _CMP_NEQ_OQ);
		     eqim = _mm_cmp_ps(_mm_set1_ps(x.imag()),
						   y.m_im,
						   _CMP_NEQ_OQ);
		     return (std::make_pair(_mm_testz_ps(eqre,eqre),
					    _mm_testz_ps(eqim,eqim))); 
		  }
	       
               
	       
	       
               static inline
	       std::pair<int32_t,int32_t>
	       operator>(const SSEc4f32 x,
			 const SSEc4f32 y) {

		     const __m128 eqre,eqim;
		     eqre = _mm_setzero_ps();
		     eqim = _mm_setzero_ps();
		     eqre = _mm_cmp_ps(x.m_re,y.m_re,
					     _CMP_GT_OQ);
		     eqim = _mm_cmp_ps(x.m_im,y.m_im,
			                     _CMP_GT_OQ);
                      return (std::make_pair(_mm_testz_ps(eqre,eqre),
					     _mm_testz_ps(eqim,eqim)));
 		   }
	       
                 
		 
		 
                 static inline
		 std::pair<int32_t,int32_t>
		 operator>(const SSEc4f32 x,
		           const std::complex<float> y) {

                        const __m128 eqre,eqim;
			eqre = _mm_setzero_ps();
			eqim = _mm_setzero_ps();
			eqre = _mm_cmp_ps(x.m_re,
				          _mm_set1_ps(y.real()),
						     _CMP_GT_OQ);
			eqim = _mm_cmp_ps(x.m_im,
					  _mm_set1_ps(y.imag()),
						      _CMP_GT_OQ);
			return (std::make_pair(_mm_testz_ps(eqre,eqre),
					       _mm_testz_ps(eqim,eqim))); 
		    }
		 
                   
		   
		   
                   static inline
		   std::pair<int32_t,int32_t>
		   operator>(const std::complex<float> x,
			     const SSEc4f32 y) {

                           const __m128 eqre,eqim;
			   eqre = _mm_setzero_ps();
			   eqim = _mm_setzero_ps();
			   eqre = _mm_cmp_ps(_mm_set1_ps(x.real()),
						         y.m_re,
						       _CMP_GT_OQ);
			   eqim = _mm_cmp_ps(_mm_set1_ps(x.imag()),
						         y.m_im,
							 _CMP_GT_OQ);
			   return (std::make_pair(_mm_testz_ps(eqre,eqre),
						  _mm_testz_ps(eqim,eqim))); 
		     }
		   
                   
		   
		  
                  static inline
		  std::pair<int32_t,int32_t>
		  operator<(const SSEc4f32 x,
			    const SSEc4f32 y) {

			const __m128 eqre,eqim;
			eqre = _mm_setzero_ps();
			eqim = _mm_setzero_ps();
			eqre = _mm_cmp_ps(x.m_re,y.m_re,
						 _CMP_LT_OQ);
			eqim = _mm_cmp_ps(x.m_im,y.m_im,
						 _CMP_LT_OQ);
                       return (std::make_pair(_mm_testz_ps(eqre,eqre),
					      _mm_testz_ps(eqim,eqim)));
 		      }
		   
                   
		   
		   
                   static inline
		   std::pair<int32_t,int32_t>
		   operator<(const SSEc4f32 x,
			     const std::complex<float> y) {

                         const __m128 eqre,eqim;
			 eqre = _mm_setzero_ps();
			 eqim = _mm_setzero_ps();
			 eqre = _mm_cmp_ps(x.m_re,
					   _mm_set1_ps(y.real()),
						      _CMP_LT_OQ);
			 eqim = _mm_cmp_ps(x.m_im,
					  _mm_set1_ps(y.imag()),
						      _CMP_LT_OQ);
		     return (std::make_pair(_mm_testz_ps(eqre,eqre),
					    _mm_testz_ps(eqim,eqim))); 
		      }
		   
                   
		   
		   
                   static inline
		   std::pair<int32_t,int32_t>
		   operator<(const std::complex<float> x,
			     const SSEc4f32 y) {

                           const __m128 eqre,eqim;
			   eqre = _mm_setzero_ps();
			   eqim = _mm_setzero_ps();
			   eqre = _mm_cmp_ps(_mm_set1_ps(x.real()),
						            y.m_re,
							 _CMP_LT_OQ);
			   eqim = _mm_cmp_ps(_mm_set1_ps(x.imag()),
						             y.m_im,
							 _CMP_LT_OQ);
			   return (std::make_pair(_mm_testz_ps(eqre,eqre),
						  _mm_testz_ps(eqim,eqim))); 
			}
		   
                   
		   
		   
                   static inline
		   std::pair<int32_t,int32_t>
		  operator>=(const SSEc4f32 x,
			     const SSEc4f32 y) {

			 const __m128 eqre,eqim;
			 eqre = _mm_setzero_ps();
			 eqim = _mm_setzero_ps();
			 eqre = _mm_cmp_ps(x.m_re,y.m_re,
						_CMP_GE_OQ);
			 eqim = _mm_cmp_ps(x.m_im,y.m_im,
						_CMP_GE_OQ);
                         return (std::make_pair(_mm_testz_ps(eqre,eqre),
						_mm_testz_ps(eqim,eqim)));
 			 }
		   
                     
		     
		     
                     static inline
		     std::pair<int32_t,int32_t>
		     operator>=(const SSEc4f32 x,
				const std::complex<float> y) {

                             const __m128 eqre,eqim;
			     eqre = _mm_setzero_ps();
			     eqim = _mm_setzero_ps();
			     eqre = _mm_cmp_ps(x.m_re,
						  _mm_set1_ps(y.real()),
								_CMP_GE_OQ);
			     eqim = _mm_cmp_ps(x.m_im,
						  _mm_set1_ps(y.imag()),
								_CMP_GE_OQ);
			     return (std::make_pair(_mm_testz_ps(eqre,eqre),
						    _mm_testz_ps(eqim,eqim))); 

		     }
		     
                   
		   
		   
                   static inline
		   std::pair<int32_t,int32_t>
		   operator>=(const std::complex<float> x,
			     const SSEc4f32 y) {

                          const __m128 eqre,eqim;
			  eqre = _mm_setzero_ps();
			  eqim = _mm_setzero_ps();
			  eqre = _mm_cmp_ps(_mm_set1_ps(x.real()),
						         y.m_re,
							  _CMP_GE_OQ);
			  eqim = _mm_cmp_ps(_mm_set1_ps(x.imag()),
						         y.m_im,
							 _CMP_GE_OQ);
			  return (std::make_pair(_mm_testz_ps(eqre,eqre),
						 _mm_testz_ps(eqim,eqim))); 
		      }
		   
                  
		  
		  
                  static inline
		  std::pair<int32_t,int32_t>
		  operator<=(const SSEc4f32 x,
			     const SSEc4f32 y) {

			  const __m128 eqre,eqim;
			  eqre = _mm_setzero_ps();
			  eqim = _mm_setzero_ps();
			  eqre = _mm_cmp_ps(x.m_re,y.m_re,
						    _CMP_LE_OQ);
			  eqim = _mm_cmp_ps(x.m_im,y.m_im,
						    _CMP_LE_OQ);
                          return (std::make_pair(_mm_testz_ps(eqre,eqre),
						 _mm_testz_ps(eqim,eqim)));
 		       }
		  
                     
		     
		     
                     static inline
		     std::pair<int32_t,int32_t>
		     operator<=(const SSEc4f32 x,
				const std::complex<float> y) {

                            const __m128 eqre,eqim;
			    eqre = _mm_setzero_ps();
			    eqim = _mm_setzero_ps();
			    eqre = _mm_cmp_ps(x.m_re,
					      _mm_set1_ps(y.real()),
							    _CMP_LE_OQ);
			    eqim = _mm_cmp_ps(x.m_im,
					      _mm_set1_ps(y.imag()),
							    _CMP_LE_OQ);
			    return (std::make_pair(_mm_testz_ps(eqre,eqre),
						   _mm_testz_ps(eqim,eqim))); 
			  }
 
                  
		  
		  
                  static inline
		  std::pair<int32_t,int32_t>
		  operator<=(const std::complex<float> x,
			     const SSEc4f32 y) {

                         const __m128 eqre,eqim;
			 eqre = _mm_setzero_ps();
			 eqim = _mm_setzero_ps();
			 eqre = _mm256_cmp_ps(_mm256_set1_ps(x.real()),
						                y.m_re,
							   _CMP_LE_OQ);
			 eqim = _mm_cmp_ps(_mm256_set1_ps(x.imag()),
			       		                         y.m_im,
							    _CMP_LE_OQ);
			 return (std::make_pair(_mm_testz_ps(eqre,eqre),
						_mm_testz_ps(eqim,eqim))); 
		         }
	 
#endif
	     
     }


}

#endif /*__GMS_SSEC4F32_H__*/
