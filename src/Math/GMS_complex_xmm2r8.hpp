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
#ifndef __GMS_COMPLEX_XMM2R8_HPP__
#define __GMS_COMPLEX_XMM2R8_HPP__ 181020231557


namespace file_version {

    const unsigned int GMS_COMPLEX_XMM2R8_MAJOR = 1U;
    const unsigned int GMS_COMPLEX_XMM2R8_MINOR = 0U;
    const unsigned int GMS_COMPLEX_XMM2R8_MICRO = 0U;
    const unsigned int GMS_COMPLEX_XMM2R8_FULLVER =
      1000U*GMS_COMPLEX_XMM2R8_MAJOR+
      100U*GMS_COMPLEX_XMM2R8_MINOR+
      10U*GMS_COMPLEX_XMM2R8_MICRO;
    const char * const GMS_COMPLEX_XMM2R8_CREATION_DATE = "22-10-2023 10:35 AM +00200 (SAT 22 OCT 2022 GMT+2)";
    const char * const GMS_COMPLEX_XMM2R8_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_COMPLEX_XMM2R8_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_COMPLEX_XMM2R8_DESCRIPTION   = "SSE optimized complex number implementation."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_simd_utils.hpp"

namespace  gms {


       namespace math {
       
       
                   struct __ATTR_ALIGN__(16) xmm2c8_t {
                   
                          __m128d re;
                          __m128d im;
                   };
                   
                   

                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cadd_xmm2c8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       const double * __restrict yre,
                                       const double * __restrict yim,
                                       double *       __restrict zre,
                                       double *       __restrict zim) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_loadu_pd(&xre[0]);
                        xmm1  = _mm_loadu_pd(&yre[0]);
                        _mm_storeu_pd(&zre[0], _mm_add_pd(xmm0,xmm1));
                        xmm2  = _mm_loadu_pd(&xim[0]);
                        xmm3  = _mm_loadu_pd(&yim[0]);
                        _mm_storeu_pd(&zim[0], _mm_add_pd(xmm2,xmm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cadd_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) xre,
                                       const double * __restrict __ATTR_ALIGN__(16) xim,
                                       const double * __restrict __ATTR_ALIGN__(16) yre,
                                       const double * __restrict __ATTR_ALIGN__(16) yim,
                                       double *       __restrict __ATTR_ALIGN__(16) zre,
                                       double *       __restrict __ATTR_ALIGN__(16) zim) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_load_pd(&xre[0]);
                        xmm1  = _mm_load_pd(&yre[0]);
                        _mm_store_pd(&zre[0], _mm_add_pd(xmm0,xmm1));
                        xmm2  = _mm_load_pd(&xim[0]);
                        xmm3  = _mm_load_pd(&yim[0]);
                        _mm_store_pd(&zim[0], _mm_add_pd(xmm2,xmm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cadd_xmm2c8(const __m128d xre,
                                     const __m128d xim,
                                     const __m128d yre,
                                     const __m128d yim,
                                     __m128d * __restrict zre,
                                     __m128d * __restrict zim) {
                     
                        register __m128d xmm0,xmm1;
                        xmm0  = _mm_add_pd(xre,yre);
                        *zre  = xmm0;
                        xmm1  = _mm_add_pd(xim,yim);
                        *zim  = xmm1;
                }
                
                  __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cadd_xmm2c8(const __m128d xre,
                                          const __m128d xim,
                                          const __m128d yre,
                                          const __m128d yim) {
                                     
                        xmm2c8_t cv;
                        register __m128d xmm0,xmm1;
                        xmm0   = _mm_add_pd(xre,yre);
                        cv.re  = xmm0;
                        xmm1   = _mm_add_pd(xim,yim);
                        cv.im  = xmm1;  
                        return (cv);            
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cadd_xmm2c8(const xmm2c8_t x,
                                          const xmm2c8_t y) {
                                     
                        xmm2c8_t cv;
                        register __m128d xmm0,xmm1;
                        xmm0   = _mm_add_pd(x.re,y.re);
                        cv.re  = xmm0;
                        xmm1   = _mm_add_pd(x.im,y.im);
                        cv.im  = xmm1;  
                        return (cv);            
               }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cadd_xmm2c8(const __m128d xre,
                                     const __m128d xim,
                                     const __m128d s,
                                     __m128d * __restrict     zre,
                                     __m128d * __restrict     zim) {

                        *zre = _mm_add_pd(xre,s);
                        *zim = _mm_setzero_pd();
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           xmm2c8_t cadd_xmm2c8(const __m128d xre,
                                          const __m128d xim,
                                          const __m128d s) {
                      
                      xmm2c8_t cv;
                      cv.re =  _mm_add_pd(xre,s);
                      cv.im =  _mm_setzero_pd();
                      return (cv);                       
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           xmm2c8_t cadd_xmm2c8(const xmm2c8_t x,
                                          const __m128d s) {
                      
                      xmm2c8_t cv;
                      cv.re =  _mm_add_pd(x.re,s);
                      cv.im =  _mm_setzero_pd();
                      return (cv);                       
                }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cadd_xmm2c8_uip(const double * __restrict xre,
                                         const double * __restrict xim,
                                         double *       __restrict zre,
                                         double *       __restrict zim) {
                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_loadu_pd(&xre[0]);
                        xmm1  = _mm_loadu_pd(&xim[0]);
                        xmm2  = _mm_loadu_pd(&zre[0]);
                        xmm3  = _mm_loadu_pd(&zim[0])
                        _mm_storeu_pd(&zre[0], _mm_add_pd(xmm2,xmm0));
                        _mm_storeu_pd(&zim[0], _mm_add_pd(xmm3,xmm1));
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cadd_xmm2c8_aip(const double * __restrict __ATTR_ALIGN__(16) xre,
                                         const double * __restrict __ATTR_ALIGN__(16) xim,
                                         double *       __restrict __ATTR_ALIGN__(16) zre,
                                         double *       __restrict __ATTR_ALIGN__(16) zim) {
                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_load_pd(&xre[0]);
                        xmm1  = _mm_load_pd(&xim[0]);
                        xmm2  = _mm_load_pd(&zre[0]);
                        xmm3  = _mm_load_pd(&zim[0])
                        _mm_store_pd(&zre[0], _mm_add_pd(xmm2,xmm0));
                        _mm_store_pd(&zim[0], _mm_add_pd(xmm3,xmm1));
              }


                ////////////////////////////////////////////////////////////////////

                  __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void csub_xmm2c8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       const double * __restrict yre,
                                       const double * __restrict yim,
                                       double *       __restrict zre,
                                       double *       __restrict zim) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_loadu_pd(&xre[0]);
                        xmm1  = _mm_loadu_pd(&yre[0]);
                        _mm_storeu_pd(&zre[0], _mm_sub_pd(xmm0,xmm1));
                        xmm2  = _mm_loadu_pd(&xim[0]);
                        xmm3  = _mm_loadu_pd(&yim[0]);
                        _mm_storeu_pd(&zim[0], _mm_sub_pd(xmm2,xmm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void csub_xmm2c8_a( const double * __restrict __ATTR_ALIGN__(16) xre,
                                       const double * __restrict __ATTR_ALIGN__(16) xim,
                                       const double * __restrict __ATTR_ALIGN__(16) yre,
                                       const double * __restrict __ATTR_ALIGN__(16) yim,
                                       double *       __restrict __ATTR_ALIGN__(16) zre,
                                       double *       __restrict __ATTR_ALIGN__(16) zim) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_load_pd(&xre[0]);
                        xmm1  = _mm_load_pd(&yre[0]);
                        _mm_store_pd(&zre[0], _mm_sub_pd(xmm0,xmm1));
                        xmm2  = _mm_load_pd(&xim[0]);
                        xmm3  = _mm_load_pd(&yim[0]);
                        _mm_store_pd(&zim[0], _mm_sub_pd(xmm2,xmm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void csub_xmm2c8(const __m128d xre,
                                     const __m128d xim,
                                     const __m128d yre,
                                     const __m128d yim,
                                     __m128d * __restrict  zre,
                                     __m128d * __restrict  zim) {
                     
                        register __m128d xmm0,xmm1;
                        xmm0  = _mm_sub_pd(xre,yre);
                        *zre  = xmm0;
                        xmm1  = _mm_sub_pd(xim,yim);
                        *zim  = xmm1;
                }
                
                
                
                
                
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t csub_xmm2c8(const __m128d xre,
                                          const __m128d xim,
                                          const __m128d yre,
                                          const __m128d yim) {
                                    
                        xmm2c8_t cv;
                        register __m128d xmm0,xmm1;
                        xmm0  = _mm_sub_pd(xre,yre);
                        cv.re  = xmm0;
                        xmm1  = _mm_sub_pd(xim,yim);
                        cv.im  = xmm1;
                        return (cv);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t csub_xmm2c8(const xmm2c8_t x,
                                          const xmm2c8_t y) {
                                    
                        xmm2c8_t cv;
                        register __m128d xmm0,xmm1;
                        xmm0  = _mm_sub_pd(x.re,y.re);
                        cv.re  = xmm0;
                        xmm1  = _mm_sub_pd(x.im,y.im);
                        cv.im  = xmm1;
                        return (cv);
                }
                
                


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void csub_xmm2c8(const __m128d xre,
                                     const __m128d xim,
                                     const __m128d s,
                                     __m128d * __restrict     zre,
                                     __m128d * __restrict     zim) {

                        *zre = _mm_sub_pd(xre,s);
                        *zim = _mm_setzero_pd();
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t csub_xmm2c8(const __m128d xre,
                                     const __m128d xim,
                                     const __m128d s) {
                                    
                        xmm2c8_t cv;
                        cv.re = _mm_sub_pd(xre,s);
                        cv.im = _mm_setzero_pd();
                        return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t csub_xmm2c8(const xmm2c8_t x,
                                          const __m128d s) {
                                    
                        xmm2c8_t cv;
                        cv.re = _mm_sub_pd(x.re,s);
                        cv.im = _mm_setzero_pd();
                        return (cv);
               }
               


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void csub_xmm2c8_uip( const double * __restrict xre,
                                         const double * __restrict xim,
                                         double *       __restrict zre,
                                         double *       __restrict zim) {
                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_loadu_pd(&xre[0]);
                        xmm1  = _mm_loadu_pd(&xim[0]);
                        xmm2  = _mm_loadu_pd(&zre[0]);
                        xmm3  = _mm_loadu_pd(&zim[0])
                        _mm_storeu_pd(&zre[0], _mm_sub_pd(xmm2,xmm0));
                        _mm_storeu_pd(&zim[0], _mm_sub_pd(xmm3,xmm1));
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void csub_xmm2c8_aip(const double * __restrict __ATTR_ALIGN__(16) xre,
                                         const double * __restrict __ATTR_ALIGN__(16) xim,
                                         double *       __restrict __ATTR_ALIGN__(16) zre,
                                         double *       __restrict __ATTR_ALIGN__(16) zim) {
                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_load_pd(&xre[0]);
                        xmm1  = _mm_load_pd(&xim[0]);
                        xmm2  = _mm_load_pd(&zre[0]);
                        xmm3  = _mm_load_pd(&zim[0])
                        _mm_store_pd(&zre[0], _mm_sub_pd(xmm2,xmm0));
                        _mm_store_pd(&zim[0], _mm_sub_pd(xmm3,xmm1));
              }


               ////////////////////////////////////////////////////////////////////////

                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmul_xmm2c8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       const double * __restrict yre,
                                       const double * __restrict yim,
                                       double *       __restrict zre,
                                       double *       __restrict zim) {

                           register __m128d xmm0,xmm1,xmm2,xmm3,xmm4,xmm5;
                           xmm0  = _mm_loadu_pd(&xre[0]);
                           xmm1  = _mm_loadu_pd(&yre[0]);
                           xmm2  = _mm_loadu_pd(&xim[0]);
                           xmm3  = _mm_loadu_pd(&yim[0]);
                           xmm4  = _mm_sub_pd(_mm_mul_pd(xmm0,xmm1),
                                                   _mm_mul_pd(xmm2,xmm3));
                           _mm_storeu_pd(&zre[0], xmm4);
                           xmm5  = _mm_mul_pd(_mm_mul_pd(xmm2,xmm1),
                                                   _mm_mul_pd(xmm0,xmm3));
                           _mm_storeu_pd(&zim[0], xmm5);
               }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmul_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) xre,
                                       const double * __restrict __ATTR_ALIGN__(16) xim,
                                       const double * __restrict __ATTR_ALIGN__(16) yre,
                                       const double * __restrict __ATTR_ALIGN__(16) yim,
                                       double *       __restrict __ATTR_ALIGN__(16) zre,
                                       double *       __restrict __ATTR_ALIGN__(16) zim) {

                           register __m128d xmm0,xmm1,xmm2,xmm3,xmm4,xmm5;
                           xmm0  = _mm_load_pd(&xre[0]);
                           xmm1  = _mm_load_pd(&yre[0]);
                           xmm2  = _mm_load_pd(&xim[0]);
                           xmm3  = _mm_load_pd(&yim[0]);
                           xmm4  = _mm_sub_pd(_mm_mul_pd(xmm0,xmm1),
                                                                        _mm_mul_pd(xmm2,xmm3));
                           _mm_store_pd(&zre[0], xmm4);
                           xmm5  = _mm_mul_pd(_mm_mul_pd(xmm2,xmm1),
                                                                        _mm_mul_pd(xmm0,xmm3));
                           _mm_store_pd(&zim[0], xmm5);
               }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmul_xmm2c8(const __m128d xre,
                                     const __m128d xim,
                                     const __m128d yre,
                                     const __m128d yim,
                                     __m128d * __restrict     zre,
                                     __m128d * __restrict     zim) {

                         register __m128d xmm0,xmm1;
                         xmm0 = _mm_sub_pd(_mm_mul_pd(xre,yre),
                                              _mm_mul_pd(xim,yim));
                         *zre  = xmm0;
                         xmm1 = _mm_mul_pd(_mm_mul_pd(xim,yre),
                                              _mm_mul_pd(xre,yim));
                         *zim  = xmm1;
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cmul_xmm2c8(const __m128d xre,
                                          const __m128d xim,
                                          const __m128d yre,
                                          const __m128d yim) {
                                     
                         xmm2c8_t cv
                         register __m128d xmm0,xmm1;
                         xmm0 = _mm_sub_pd(_mm_mul_pd(xre,yre),
                                              _mm_mul_pd(xim,yim));
                         cv.re  = xmm0;
                         xmm1 = _mm_mul_pd(_mm_mul_pd(xim,yre),
                                              _mm_mul_pd(xre,yim));
                         cv.im  = xmm1;
                         return (cv);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cmul_xmm2c8(const xmm2c8_t x,
                                          const xmm2c8_t y) {
                                     
                         xmm2c8_t cv
                         register __m128d xmm0,xmm1;
                         xmm0 = _mm_sub_pd(_mm_mul_pd(x.re,y.re),
                                              _mm_mul_pd(x.im,y.im));
                         cv.re  = xmm0;
                         xmm1 = _mm_mul_pd(_mm_mul_pd(x.im,y.re),
                                              _mm_mul_pd(x.re,y.im));
                         cv.im  = xmm1;
                         return (cv);
                }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmul_xmm2c8(const __m128d xre,
                                     const __m128d xim,
                                     const __m128d s,
                                     __m128d * __restrict   zre,
                                     __m128d * __restrict   zim) {

                        *zre = _mm_mul_pd(xre,s);
                        *zim = _mm_mul_pd(xim,s);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cmul_xmm2c8(const __m128d xre,
                                          const __m128d xim,
                                          const __m128d s) {
                                     
                        xmm2c8_t cv;
                        cv.re = _mm_mul_pd(xre,s);
                        cv.im = _mm_mul_pd(xim,s);
                        return (cv);
               }
               
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cmul_xmm2c8(const xmm2c8_t x,
                                          const __m128d s) {
                                     
                        xmm2c8_t cv;
                        cv.re = _mm_mul_pd(x.re,s);
                        cv.im = _mm_mul_pd(x.im,s);
                        return (cv);
               }
               


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmul_xmm2c8_uip(const double * __restrict xre,
                                         const double * __restrict xim,
                                         double *       __restrict zre,
                                         double *       __restrict zim) {

                           register __m128d xmm0,xmm1,xmm2,xmm3,xmm4,xmm5;
                           xmm0  = _mm_loadu_pd(&xre[0]);
                           xmm1  = _mm_loadu_pd(&zre[0]);
                           xmm2  = _mm_loadu_pd(&xim[0]);
                           xmm3  = _mm_loadu_pd(&zim[0]);
                           xmm4  = _mm_sub_pd(_mm_mul_pd(xmm0,xmm1),
                                                 _mm_mul_pd(xmm2,xmm3));
                           _mm_storeu_pd(&zre[0], xmm4);
                           xmm5  = _mm_mul_pd(_mm_mul_pd(xmm2,xmm1),
                                                 _mm_mul_pd(xmm0,xmm3));
                           _mm_storeu_pd(&zim[0], xmm5);
               }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmul_xmm2c8_aip(const double * __restrict __ATTR_ALIGN__(16) xre,
                                         const double * __restrict __ATTR_ALIGN__(16) xim,
                                         double *       __restrict __ATTR_ALIGN__(16) zre,
                                         double *       __restrict __ATTR_ALIGN__(16) zim) {

                           register __m128d xmm0,xmm1,xmm2,xmm3,xmm4,xmm5;
                           xmm0  = _mm_load_pd(&xre[0]);
                           xmm1  = _mm_load_pd(&zre[0]);
                           xmm2  = _mm_load_pd(&xim[0]);
                           xmm3  = _mm_load_pd(&zim[0]);
                           xmm4  = _mm_sub_pd(_mm_mul_pd(xmm0,xmm1),
                                                 _mm_mul_pd(xmm2,xmm3));
                           _mm_store_pd(&zre[0], xmm4);
                           xmm5  = _mm_mul_pd(_mm_mul_pd(xmm2,xmm1),
                                                 _mm_mul_pd(xmm0,xmm3));
                           _mm_store_pd(&zim[0], xmm5);
               }

                 ////////////////////////////////////////////////////////////////////////////

                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cdiv_xmm2c8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       const double * __restrict yre,
                                       const double * __restrict yim,
                                       double *       __restrict zre,
                                       double *       __restrict zim) {

                        register __m128d xmm0,xmm1,xmm2,xmm3; 
                        register __m128d xmm4,xmm5,xmm6;
                        xmm0  = _mm_loadu_pd(&xre[0]); //a
                        xmm1  = _mm_loadu_pd(&yim[0]); //d
                        xmm2  = _mm_loadu_pd(&xim[0]); //b
                        xmm3  = _mm_loadu_pd(&yre[0]); //c
                        xmm4  = _mm_fmadd_pd(xmm0,xmm3,
                                                _mm_mul_pd(xmm2,xmm1));
                        xmm5  = _mm_fmsub_pd(xmm2,xmm3,
                                                _mm_mul_pd(xmm0,xmm1));
                        xmm6  = _mm_fmadd_pd(xmm3,xmm3),
                                                _mm_mul_pd(xmm1,xmm1));
                        _mm_storeu_pd(&zre[0], _mm_div_pd(xmm4,xmm6));
                        _mm_storeu_pd(&zim[0], _mm_div_pd(xmm5,xmm6));
                }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cdiv_xmm2c8_a(const double * __restrict xre,
                                       const double * __restrict xim,
                                       const double * __restrict yre,
                                       const double * __restrict yim,
                                       double *       __restrict zre,
                                       double *       __restrict zim) {

                        register __m128d xmm0,xmm1,xmm2,xmm3; 
                        register __m128d xmm4,xmm5,xmm6;
                        xmm0  = _mm_load_pd(&xre[0]); //a
                        xmm1  = _mm_load_pd(&yim[0]); //d
                        xmm2  = _mm_load_pd(&xim[0]); //b
                        xmm3  = _mm_load_pd(&yre[0]); //c
                        xmm4  = _mm_fmadd_pd(xmm0,xmm3,
                                                _mm_mul_pd(xmm2,xmm1));
                        xmm5  = _mm_fmsub_pd(xmm2,xmm3,
                                                _mm_mul_pd(xmm0,xmm1));
                        xmm6  = _mm_fmadd_pd(xmm3,xmm3,
                                                _mm_mul_pd(xmm1,xmm1));
                        _mm_store_pd(&zre[0], _mm_div_pd(xmm4,xmm6));
                        _mm_store_pd(&zim[0], _mm_div_pd(xmm5,xmm6));
                }
              
                                         
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cdiv_xmm2c8(const __m128d xre,
                                     const __m128d xim,
                                     const __m128d yre,
                                     const __m128d yim,
                                     __m128d * __restrict zre,
                                     __m128d * __restrict zim) {

                      register __m128d xmm0,xmm1,xmm2;
                      xmm0 = _mm_fmadd_pd(xre,yre,
                                           _mm_mul_pd(xim,yim));
                      xmm1 = _mm_fmsub_pd(xim,yre,
                                           _mm_mul_pd(xre,yim));
                      xmm2 = _mm_fmadd_pd(xmm3,xmm3,
                                           _mm_mul_pd(xmm1,xmm1));
                      *zre  = _mm_div_pd(xmm0,xmm2);
                      *zim  = _mm_div_pd(xmm1,xmm2);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cdiv_xmm2c8(const __m128d xre,
                                          const __m128d xim,
                                          const __m128d yre,
                                          const __m128d yim) {
                                     
                      xmm2c8_t
                      register __m128d xmm0,xmm1,xmm2;
                      xmm0 = _mm_fmadd_pd(xre,yre,
                                           _mm_mul_pd(xim,yim));
                      xmm1 = _mm_fmsub_pd(x.im,y.re,
                                           _mm_mul_pd(xre,yim));
                      xmm2 = _mm_fmadd_pd(xmm3,xmm3,
                                           _mm_mul_pd(xmm1,xmm1));
                      cv.re  = _mm_div_pd(xmm0,xmm2);
                      cv.im  = _mm_div_pd(xmm1,xmm2);
                      return (cv);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cdiv_xmm2c8(const xmm2c8_t x,
                                          const xmm2c8_t y) {
                                     
                      xmm2c8_t
                      register __m128d xmm0,xmm1,xmm2;
                      xmm0 = _mm_fmadd_pd(x.re,y.re,
                                           _mm_mul_pd(x.im,y.im));
                      xmm1 = _mm_fmsub_pd(x.im,y.re,
                                           _mm_mul_pd(x.re,y.im));
                      xmm2 = _mm_fmadd_pd(xmm3,xmm3,
                                           _mm_mul_pd(xmm1,xmm1));
                      cv.re  = _mm_div_pd(xmm0,xmm2);
                      cv.im  = _mm_div_pd(xmm1,xmm2);
                      return (cv);
                }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cdiv_xmm2c8(const __m128d xre,
                                     const __m128d xim,
                                     const __m128d s,
                                     __m128d * __restrict zre,
                                     __m128d * __restrict zim) {

                        *zre = _mm_div_pd(xre,s);
                        *zim = _mm_div_pd(xim,s);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cdiv_xmm2c8(const __m128d xre,
                                          const __m128d xim,
                                          const __m128d s) {
                                     
                         xmm2c8_t cv;
                         cv.re = _mm_div_pd(xre,s);
                         cv.im = _mm_div_pd(xim,s);
                         return (cv);
               }
               
               
                 
               
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cdiv_xmm2c8(const xmm2c8_t x,
                                          const __m128d s) {
                                     
                         xmm2c8_t cv;
                         cv.re = _mm_div_pd(x.re,s);
                         cv.im = _mm_div_pd(x.im,s);
                         return (cv);
               }
               
               
               
               
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cdiv_xmm2c8_s(const __m128d s,
                                       const __m128d xre,
                                       const __m128d xim,
                                       __m128d * __restrict zre,
                                       __m128d * __restrict zim) {
                        
                        register __m128d t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm_setzero_pd();
                        cdiv_xmm2c8(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        *zre = tmpr;
                        *zim = tmpi;                     
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cdiv_xmm2c8_s(const __m128d s,
                                            const __m128d xre,
                                            const __m128d xim) {
                                       
                        xmm2c8_t cv;
                        register __m128d t0r,t0i;
                        t0r = s;
                        t0i = _mm_setzero_pd();
                        cdiv_xmm2c8(t0r,t0i,xre,xim,&cv.re,&cv.im);
                        return (cv);                     
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cdiv_xmm2c8_s(const __m128d s,
                                            const xmm2c8_t x) {
                                       
                        xmm2c8_t cv;
                        register __m128d t0r,t0i;
                        t0r = s;
                        t0i = _mm_setzero_pd();
                        cdiv_xmm2c8(t0r,t0i,x.re,x.im,&cv.re,&cv.im);
                        return (cv);                     
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cdiv_smith_xmm2c8_s(const __m128d s,
                                             const __m128d xre,
                                             const __m128d xim,
                                             __m128d * __restrict zre,
                                             __m128d * __restrict zim) {
                                             
                        register __m128d t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm_setzero_pd(); 
                        cdiv_smith_xmm2c8(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        *zre = tmpr;
                        *zim = tmpi;                   
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cdiv_smith_xmm2c8_s(const __m128d s,
                                                  const __m128d xre,
                                                  const __m128d xim) {
                                             
                        xmm2c8_t cv;                    
                        register __m128d t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm_setzero_pd(); 
                        cdiv_smith_xmm2c8(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        cv.re = tmpr;
                        cv.im = tmpi;  
                        return (cv);                 
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cdiv_smith_xmm2c8_s(const __m128d s,
                                                  const xmm2c8_t x) {
                                             
                        xmm2c8_t cv,t0;                    
                        t0.re = s;
                        t0.im = _mm_setzero_pd(); 
                        cv = cdiv_smith_xmm2c8(t0,x);
                        return (cv);                 
                 }
                 
                   


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cdiv_xmm2c8_uip(const double * __restrict xre,
                                         const double * __restrict xim,
                                         double *       __restrict zre,
                                         double *       __restrict zim) {

                        register __m128d xmm0,xmm1,xmm2,xmm3; 
                        register __m128d xmm4,xmm5,xmm6;
                        xmm0  = _mm_loadu_pd(&xre[0]); //a
                        xmm1  = _mm_loadu_pd(&zim[0]); //d
                        xmm2  = _mm_loadu_pd(&xim[0]); //b
                        xmm3  = _mm_loadu_pd(&zre[0]); //c
                        xmm4  = _mm_fmadd_pd(xmm0,xmm3,
                                                _mm_mul_pd(xmm2,xmm1));
                        xmm5  = _mm_fmsub_pd(xmm2,xmm3,
                                                _mm_mul_pd(xmm0,xmm1));
                        xmm6  = _mm_fmadd_pd(xmm3,xmm3,
                                                _mm_mul_pd(xmm1,xmm1));
                        _mm_storeu_pd(&zre[0], _mm_div_pd(xmm4,xmm6));
                        _mm_storeu_pd(&zim[0], _mm_div_pd(xmm5,xmm6));
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cdiv_xmm2c8_aip(const double * __restrict __ATTR_ALIGN__(16) xre,
                                         const double * __restrict __ATTR_ALIGN__(16) xim,
                                         double *       __restrict __ATTR_ALIGN__(16) zre,
                                         double *       __restrict __ATTR_ALIGN__(16) zim) {

                        register __m128d xmm0,xmm1,xmm2,xmm3; 
                        register __m128d xmm4,xmm5,xmm6;
                        xmm0  = _mm_load_pd(&xre[0]); //a
                        xmm1  = _mm_load_pd(&zim[0]); //d
                        xmm2  = _mm_load_pd(&xim[0]); //b
                        xmm3  = _mm_load_pd(&zre[0]); //c
                        xmm4  = _mm_fmadd_pd(xmm0,xmm3,
                                                _mm_mul_pd(xmm2,xmm1));
                        xmm5  = _mm_fmsub_pd(xmm2,xmm3,
                                                _mm_mul_pd(xmm0,xmm1));
                        xmm6  = _mm_fmadd_pd(xmm3,xmm3,
                                                _mm_mul_pd(xmm1,xmm1));
                        _mm_store_pd(&zre[0], _mm_div_pd(xmm4,xmm6));
                        _mm_store_pd(&zim[0], _mm_div_pd(xmm5,xmm6));
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cdiv_smith_xmm2c8_u(const double * __restrict xre,
                                             const double * __restrict xim,
                                             const double * __restrict yre,
                                             const double * __restrict yim,
                                             double *       __restrict zre,
                                             double *       __restrict zim) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        register __m128d r,den;
                        __mmask8 m = 0x0;
                        xmm0 = _mm_loadu_pd(&yre[0]); // c
                        xmm1 = _mm_loadu_pd(&yim[0]); // d
                        xmm2 = _mm_loadu_pd(&xre[0]); // a
                        xmm3 = _mm_loadu_pd(&xim[0]); // b
                        m    = _mm_cmp_pd_mask(_mm_abs_pd(xmm0),
                                                  _mm_abs_pd(xmm1),
                                                  _CMP_GE_OQ);
                        r    = _mm_mask_blend_pd(m,_mm_div_pd(xmm0,xmm1),
                                                      _mm_div_pd(xmm1,xmm0)); // r
                        den  = _mm_mask_blend_pd(m,_mm_fmadd_pd(r,xmm0,xmm1),
                                                      _mm_fmadd_pd(r,xmm1,xmm0));
                        _mm_storeu_pd(&zre[0], _mm_mask_blend_pd(m,
                                                _mm_div_pd(_mm_fmadd_pd(xmm2,r,xmm3),den),
                                                _mm_div_pd(_mm_fmadd_pd(xmm3,r,xmm2),den)));
                        _mm_storeu_pd(&zim[0], _mm_mask_blend_pd(m,
                                                _mm_div_pd(_mm_fmsub_pd(xmm3,r,xmm2),den),
                                                _mm_div_pd(_mm_sub_pd(xmm3,_mm_mul_pd(xmm2,r)),den)));
               }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cdiv_smith_xmm2c8_a(const double * __restrict xre,
                                             const double * __restrict xim,
                                             const double * __restrict yre,
                                             const double * __restrict yim,
                                             double *       __restrict zre,
                                             double *       __restrict zim) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        register __m128d r,den;
                        __mmask8 m = 0x0;
                        xmm0 = _mm_load_pd(&yre[0]); // c
                        xmm1 = _mm_load_pd(&yim[0]); // d
                        xmm2 = _mm_load_pd(&xre[0]); // a
                        xmm3 = _mm_load_pd(&xim[0]); // b
                        m    = _mm_cmp_pd_mask(_mm_abs_pd(xmm0),
                                                  _mm_abs_pd(xmm1),
                                                  _CMP_GE_OQ);
                        r    = _mm_mask_blend_pd(m,_mm_div_pd(xmm0,xmm1),
                                                      _mm_div_pd(xmm1,xmm0)); // r
                        den  = _mm_mask_blend_pd(m,_mm_fmadd_pd(r,xmm0,xmm1),
                                                      _mm_fmadd_pd(r,xmm1,xmm0));
                        _mm_storeu_pd(&zre[0], _mm_mask_blend_pd(m,
                                                _mm_div_pd(_mm_fmadd_pd(xmm2,r,xmm3),den),
                                                _mm_div_pd(_mm_fmadd_pd(xmm3,r,xmm2),den)));
                        _mm_storeu_pd(&zim[0], _mm_mask_blend_pd(m,
                                                _mm_div_pd(_mm_fmsub_pd(xmm3,r,xmm2),den),
                                                _mm_div_pd(_mm_sub_pd(xmm3,_mm_mul_pd(xmm2,r)),den)));
               }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cdiv_smith_xmm2c8(const __m128d xre,
                                           const __m128d xim,
                                           const __m128d yre,
                                           const __m128d yim,
                                           __m128d * __restrict zre,
                                           __m128d * __restrict zim) {

                        register __m128d r,den;
                        __mmask8 m = 0x0;
                        m    = _mm_cmp_pd_mask(_mm_abs_pd(yre),
                                                  _mm_abs_pd(yim),
                                                  _CMP_GE_OQ);
                        r    = _mm_mask_blend_pd(m,_mm_div_pd(yre,yim),
                                                      _mm_div_pd(yim,yre)); // r
                        den  = _mm_mask_blend_pd(m,_mm_fmadd_pd(r,yre,yim),
                                                      _mm_fmadd_pd(r,yim,yre));
                        *zre  =  _mm_mask_blend_pd(m,
                                                _mm_div_pd(_mm_fmadd_pd(xre,r,xim),den),
                                                _mm_div_pd(_mm_fmadd_pd(xim,r,xre),den));
                        *zim  =  _mm_mask_blend_pd(m,
                                                _mm_div_pd(_mm_fmsub_pd(xim,r,xre),den),
                                                _mm_div_pd(_mm_sub_pd(xim,_mm_mul_pd(xre,r)),den)));
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cdiv_smith_xmm2c8(const __m128d xre,
                                                const __m128d xim,
                                                const __m128d yre,
                                                const __m128d yim) {
                                           
                        xmm2c8_t cv
                        register __m128d r,den;
                        __mmask8 m = 0x0;
                        m    = _mm_cmp_pd_mask(_mm_abs_pd(yre),
                                                  _mm_abs_pd(yim),
                                                  _CMP_GE_OQ);
                        r    = _mm_mask_blend_pd(m,_mm_div_pd(yre,yim),
                                                      _mm_div_pd(yim,yre)); // r
                        den  = _mm_mask_blend_pd(m,_mm_fmadd_pd(r,yre,yim),
                                                      _mm_fmadd_pd(r,yim,yre));
                        cv.re  =  _mm_mask_blend_pd(m,
                                                _mm_div_pd(_mm_fmadd_pd(xre,r,xim),den),
                                                _mm_div_pd(_mm_fmadd_pd(xim,r,xre),den));
                        cv.im  =  _mm_mask_blend_pd(m,
                                                _mm_div_pd(_mm_fmsub_pd(xim,r,xre),den),
                                                _mm_div_pd(_mm_sub_pd(xim,_mm_mul_pd(xre,r)),den)));
                        return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cdiv_smith_xmm2c8(const xmm2c8_t x,
                                                const xmm2c8_t y) {
                                           
                        xmm2c8_t cv
                        register __m128d r,den;
                        __mmask8 m = 0x0;
                        m    = _mm_cmp_pd_mask(_mm_abs_pd(y.re),
                                                  _mm_abs_pd(y.im),
                                                  _CMP_GE_OQ);
                        r    = _mm_mask_blend_pd(m,_mm_div_pd(y.re,y.im),
                                                      _mm_div_pd(y.im,y.re)); // r
                        den  = _mm_mask_blend_pd(m,_mm_fmadd_pd(r,y.re,y.im),
                                                      _mm_fmadd_pd(r,y.im,y.re));
                        cv.re  =  _mm_mask_blend_pd(m,
                                                _mm_div_pd(_mm_fmadd_pd(x.re,r,x.im),den),
                                                _mm_div_pd(_mm_fmadd_pd(x.im,r,x.re),den));
                        cv.im  =  _mm_mask_blend_pd(m,
                                                _mm_div_pd(_mm_fmsub_pd(x.im,r,x.re),den),
                                                _mm_div_pd(_mm_sub_pd(x.im,_mm_mul_pd(x.re,r)),den)));
                        return (cv);
               }





                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cabs_xmm2c8_u(const double * __restrict re,
                                       const double * __restrict im,
                                       double * __restrict  cabs) {

                        register __m128d xmm0,xmm1,xmm2,xmm3,xmm4;
                        xmm0  = _mm_loadu_pd(&re[0]);
                        xmm1  = _mm_mul_pd(xmm0,xmm0);
                        xmm2  = _mm_loadu_pd(&im[0]);
                        xmm3  = _mm_mul_pd(xmm2,xmm2);
                        xmm4  = _mm_sqrt_pd(_mm_add_pd(xmm1,xmm3));
                        _mm_storeu_pd(&cabs[0],xmm4);
                 }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cabs_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) re,
                                       const double * __restrict __ATTR_ALIGN__(16) im,
                                       double * __restrict  __ATTR_ALIGN__(16) cabs) {

                        register __m128d xmm0,xmm1,xmm2,xmm3,xmm4;
                        xmm0  = _mm_load_pd(&re[0]);
                        xmm1  = _mm_mul_pd(xmm0,xmm0);
                        xmm2  = _mm_load_pd(&im[0]);
                        xmm3  = _mm_mul_pd(xmm2,xmm2);
                        xmm4  = _mm_sqrt_pd(_mm_add_pd(xmm1,xmm3));
                        _mm_store_pd(&cabs[0],xmm4);
                 }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   __m128d cabs_xmm2c8(const __m128d re,
                                       const __m128d im) {

                        register __m128d xmm0,xmm1,cabs;
                        xmm0 = _mm_mul_pd(re,re);
                        xmm1 = _mm_mul_pd(im,im);
                        cabs = _mm_sqrt_pd(_mm_add_pd(xmm0,xmm1));
                        return (cabs);
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   __m128d cabs_xmm2c8(const xmm2c8_t x) {

                        register __m128d xmm0,xmm1,cabs;
                        xmm0 = _mm_mul_pd(x.re,x.re);
                        xmm1 = _mm_mul_pd(x.im,x.im);
                        cabs = _mm_sqrt_pd(_mm_add_pd(xmm0,xmm1));
                        return (cabs);
                 }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void carg_xmm2c8_u(const double * __restrict re,
                                       const double * __restrict im,
                                       double * __restrict  carg) {

                        register __m128d xmm0,xmm1;
                        xmm0 = _mm_loadu_pd(&re[0]);
                        xmm1 = _mm_loadu_pd(&im[0]);
                        _mm_storeu_pd(&carg[0],_mm_atan2_pd(xmm0,xmm1));
                }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void carg_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) re,
                                       const double * __restrict __ATTR_ALIGN__(16) im,
                                       double * __restrict  __ATTR_ALIGN__(16) carg) {

                        register __m128d xmm0,xmm1;
                        xmm0 = _mm_load_pd(&re[0]);
                        xmm1 = _mm_load_pd(&im[0]);
                        _mm_store_pd(&carg[0],_mm_atan2_pd(xmm0,xmm1));
                }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   __m128d carg_xmm2c8(const __m128d re,
                                       const __m128d im) {

                       register __m128d carg;
                       carg = _mm_atan2_pd(re,im);
                       return (carg);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   __m128d carg_xmm2c8(xmm2c8_t x) {

                       register __m128d carg;
                       carg = _mm_atan2_pd(x.re,x.im);
                       return (carg);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           void clog_xmm2c8(const __m128d re,
	                             const __m128d im,
	                             __m128d * __restrict clogr,
	                             __m128d * __restrict clogi) {
	                
	                register __m128d t1,t2,ln;
	                t1  = cabs_xmm2c8(re,im);
	                t2  = carg_xmm2c8(re,im);
	                ln  = _mm_log_pd(t1);
	                *clogr = ln;
	                *clogi = t2;                    
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           void clog_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) pre,
	                               const double * __restrict __ATTR_ALIGN__(16) pim,
	                               double * __restrict clogr,
	                               double * __restrict clogi) {
	                
	                register __m128d re = _mm_load_pd(&pre[0]);
	                register __m128d im = _mm_load_pd(&pim[0]);
	                register __m128d t1,t2,ln;
	                t1  = cabs_xmm2c8(re,im);
	                t2  = carg_xmm2c8(re,im);
	                ln  = _mm_log_pd(t1);
	                _mm_store_pd(&clogr[0] ,ln);
	                _mm_store_pd(&clogi[0] ,t2);                    
	        }
	        
	        
	          __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           void clog_xmm2c8_u(const double * __restrict  pre,
	                               const double * __restrict  pim,
	                               double * __restrict clogr,
	                               double * __restrict clogi) {
	                
	                register __m128d re = _mm_loadu_pd(&pre[0]);
	                register __m128d im = _mm_loadu_pd(&pim[0]);
	                register __m128d t1,t2,ln;
	                t1  = cabs_xmm2c8(re,im);
	                t2  = carg_xmm2c8(re,im);
	                ln  = _mm_log_pd(t1);
	                _mm_storeu_pd(&clogr[0] ,ln);
	                _mm_storeu_pd(&clogi[0] ,t2);                    
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           xmm2c8_t clog_xmm2c8(const xmm2c8_t x){
	                                  
	                xmm2c8_t clog;                           
	                register __m128d t1,t2,ln;
	                t1  = cabs_xmm2c8(x.re,x.im);
	                t2  = carg_xmm2c8(x.re,x.im);
	                ln  = _mm_log_pd(t1);
	                clog.re = ln;
	                clog.im = t2;
	                return (clog);                    
	        }

                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cconj_xmm2c8_u(double * __restrict re,
                                        double * __restrict im) {

                        register __m128d c;
                        c = negate_xmm2r8(_mm_loadu_pd(&im[0]));
                        _mm_storeu_pd(&im[0],c);
               }  


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cconj_xmm2c8_a(double * __restrict __ATTR_ALIGN__(16) re,
                                        double * __restrict __ATTR_ALIGN__(16) im) {
                                        
                        register __m128d c;
                        c = negate_xmm2r8(_mm_load_pd(&im[0]));
                        _mm_store_pd(&im[0],c);
               }  


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cconj_xmm2c8(__m128d * __restrict re,
                                      __m128d * __restrict im) {
                         
                        register __m128d c;              
                        c = negate_xmm2r8(*im);
                        *im = c;
                   } 
                   
                   
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cconj_xmm2c8_v2(const __m128d xre,
                                         const __m128d xim,
                                         __m128d * __restrict yre,
                                         __m128d * __restrict yim) {
                         
                        //register __m128d c;              
                        //c = negate_xmm2c8(*im);
                        //*im = c;
                        *yre = xre; 
                        *yim = negate_xmm2r8(xim);
                   } 
                   
                   
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cconj_xmm2c8_v2(const __m128d xre,
                                              const __m128d xim) {                                              
                         
                        //register __m128d c;              
                        //c = negate_xmm2c8(*im);
                        //*im = c;
                        xmm2c8_t cv;
                        cv.re = xre; 
                        cv.im = negate_xmm2r8(xim);
                        return (cv);
                   } 
                   
                   
                
	          
                  
	           static inline
                   xmm2c8_t cconj_xmm2c8_v2(const xmm2c8_t x) {                                              
                         
                        //register __m128d c;              
                        //c = negate_xmm2c8(*im);
                        //*im = c;
                        xmm2c8_t cv;
                        cv.re = x.re; 
                        cv.im = negate_xmm2r8(x.im);
                        return (cv);
                   } 
                   
                   


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void ccos_xmm2c8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       double * __restrict  csre,
                                       double * __restrict  csim) {

                      register __m128d xmm0,xmm1,xmm2,xmm3;
                      xmm0  = _mm_loadu_pd(&xre[0]);
                      xmm1  = _mm_loadu_pd(&xim[0]);
                      xmm2  = _mm_mul_pd(_mm_cos_pd(xmm0),_mm_cosh_pd(xmm1));
                      _mm_storeu_pd(&csre[0],xmm2);
                      xmm3  = _mm_mul_pd(_mm_sin_pd(xmm0),_mm_sinh_pd(xmm1));
                      _mm_storeu_pd(&csim[0],xmm3);
               }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void ccos_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) xre,
                                       const double * __restrict __ATTR_ALIGN__(16) xim,
                                       double * __restrict  __ATTR_ALIGN__(16) csre,
                                       double * __restrict  __ATTR_ALIGN__(16) csim) {

                      register __m128d xmm0,xmm1,xmm2,xmm3;
                      xmm0  = _mm_load_pd(&xre[0]);
                      xmm1  = _mm_load_pd(&xim[0]);
                      xmm2  = _mm_mul_pd(_mm_cos_pd(xmm0),_mm_cosh_pd(xmm1));
                      _mm_store_pd(&csre[0],xmm2);
                      xmm3  = _mm_mul_pd(_mm_sin_pd(xmm0),_mm_sinh_pd(xmm1));
                      _mm_store_pd(&csim[0],xmm3);
               }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void ccos_xmm2c8(const __m128d xre,
                                     const __m128d xim,
                                     __m128d * __restrict csre,
                                     __m128d * __restrict csim) {

                      register __m128d xmm0,xmm1;
                      xmm0  = _mm_mul_pd(_mm_cos_pd(xre),_mm_cosh_pd(xim));
                      *csre = xmm0;
                      xmm1  = _mm_mul_pd(_mm_sin_pd(xre),_mm_sinh_pd(xim));
                      *csim = xmm1; 
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t ccos_xmm2c8(const __m128d xre,
                                          const __m128d xim) {
                                    
                      xmm2c8_t cv;
                      register __m128d xmm0,xmm1;
                      xmm0  = _mm_mul_pd(_mm_cos_pd(xre),_mm_cosh_pd(xim));
                      cv.re = xmm0;
                      xmm1  = _mm_mul_pd(_mm_sin_pd(xre),_mm_sinh_pd(xim));
                      cv.im = xmm1;
                      return (cv); 
               }
               
               
                  __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t ccos_xmm2c8(const xmm2c8_t x) {
                                    
                      xmm2c8_t cv;
                      register __m128d xmm0,xmm1;
                      xmm0  = _mm_mul_pd(_mm_cos_pd(x.re),_mm_cosh_pd(x.im));
                      cv.re = xmm0;
                      xmm1  = _mm_mul_pd(_mm_sin_pd(x.re),_mm_sinh_pd(x.im));
                      cv.im = xmm1;
                      return (cv); 
               }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void ccosh_xmm2c8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       double * __restrict  csre,
                                       double * __restrict  csim) {

                      register __m128d xmm0,xmm1,xmm2,xmm3;
                      xmm0  = _mm_loadu_pd(&xre[0]);
                      xmm1  = _mm_loadu_pd(&xim[0]);
                      xmm2  = _mm_mul_pd(_mm_cosh_pd(xmm0),_mm_cos_pd(xmm1));
                      _mm_storeu_pd(&csre[0],xmm2);
                      xmm3  = _mm_mul_pd(_mm_sinh_pd(xmm0),_mm_sin_pd(xmm1));
                      _mm_storeu_pd(&csim[0],xmm3);
               }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void ccosh_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) xre,
                                       const double * __restrict __ATTR_ALIGN__(16) xim,
                                       double * __restrict  __ATTR_ALIGN__(16) csre,
                                       double * __restrict  __ATTR_ALIGN__(16) csim) {

                      register __m128d xmm0,xmm1,xmm2,xmm3;
                      xmm0  = _mm_load_pd(&xre[0]);
                      xmm1  = _mm_load_pd(&xim[0]);
                      xmm2  = _mm_mul_pd(_mm_cosh_pd(xmm0),_mm_cos_pd(xmm1));
                      _mm_store_pd(&csre[0],xmm2);
                      xmm3  = _mm_mul_pd(_mm_sinh_pd(xmm0),_mm_sin_pd(xmm1));
                      _mm_store_pd(&csim[0],xmm3);
               }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void ccosh_xmm2c8(const __m128d xre,
                                     const __m128d xim,
                                     __m128d * __restrict csre,
                                     __m128d * __restrict csim) {

                      register __m128d xmm0,xmm1;
                      xmm0  = _mm_mul_pd(_mm_cosh_pd(xre),_mm_cos_pd(xim));
                      *csre = xmm0;
                      xmm1  = _mm_mul_pd(_mm_sinh_pd(xre),_mm_sin_pd(xim));
                      *csim = xmm1; 
               }
               
               
                  __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t ccosh_xmm2c8(const __m128d xre,
                                           const __m128d xim) {
                                          
                      xmm2c8_t cv;
                      register __m128d xmm0,xmm1;
                      xmm0  = _mm_mul_pd(_mm_cosh_pd(xre),_mm_cos_pd(xim));
                      cv.re = xmm0;
                      xmm1  = _mm_mul_pd(_mm_sinh_pd(xre),_mm_sin_pd(xim));
                      cv.im = xmm1; 
                      return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t ccosh_xmm2c8(const xmm2c8_t x) {
                                          
                      xmm2c8_t cv;
                      register __m128d xmm0,xmm1;
                      xmm0  = _mm_mul_pd(_mm_cosh_pd(x.re),_mm_cos_pd(x.im));
                      cv.re = xmm0;
                      xmm1  = _mm_mul_pd(_mm_sinh_pd(x.re),_mm_sin_pd(x.im));
                      cv.im = xmm1; 
                      return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           void cpow_xmm2c8(const __m128d xre,
	                             const __m128d xim,
	                             const double n,
	                             __m128d * __restrict powr,
	                             __m128d * __restrict powi) {
	                             
	                register __m128d xmm0,xmm1;
	                register __m128d r,tht;
	                register __m128d vn,pt;
	                register __m128d ta;
	                xmm0  = _mm_mul_pd(xre,xre);
	                vn    = _mm_set1_pd(n);
	                xmm1  = _mm_mul_pd(xim,xim);
	                r     = _mm_sqrt_pd(_mm_add_pd(xmm0,xmm1));
	                tht   = _mm_atan_pd(_mm_div_pd(xim,xre));
	                pt    = _mm_pow_pd(r,vn);
	                ta    = _mm_mul_pd(vn,tht);
	                *powr = _mm_mul_pd(pt,_mm_cos_pd(ta));
	                *powi = _mm_mul_pd(pt,_mm_sin_pd(ta));                    
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           void cpow_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) pxre,
	                               const double * __restrict __ATTR_ALIGN__(16) pxim,
	                               const double n,
	                               double * __restrict __ATTR_ALIGN__(16) ppowr
	                               double * __restrict __ATTR_ALIGN__(16) ppowi) {
	                  
	                register __m128d xre = _mm_load_pd(&pxre[0]);
	                register __m128d xim = _mm_load_pd(&pxim[0]);           
	                register __m128d xmm0,xmm1;
	                register __m128d r,tht;
	                register __m128d vn,pt;
	                register __m128d ta;
	                xmm0  = _mm_mul_pd(xre,xre);
	                vn    = _mm_set1_pd(n);
	                xmm1  = _mm_mul_pd(xim,xim);
	                r     = _mm_sqrt_pd(_mm_add_pd(xmm0,xmm1));
	                tht   = _mm_atan_pd(_mm_div_pd(xim,xre));
	                pt    = _mm_pow_pd(r,vn);
	                ta    = _mm_mul_pd(vn,tht);
	                _mm_store_pd(&ppowr[0] ,_mm_mul_pd(pt,_mm_cos_pd(ta)));
	                _mm_store_pd(&ppowi[0] ,_mm_mul_pd(pt,_mm_sin_pd(ta)));                    
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           void cpow_xmm2c8_u(const double * __restrict  pxre,
	                               const double * __restrict  pxim,
	                               const double n,
	                               double * __restrict  ppowr
	                               double * __restrict  ppowi) {
	                  
	                register __m128d xre = _mm_loadu_pd(&pxre[0]);
	                register __m128d xim = _mm_loadu_pd(&pxim[0]);           
	                register __m128d xmm0,xmm1;
	                register __m128d r,tht;
	                register __m128d vn,pt;
	                register __m128d ta;
	                xmm0  = _mm_mul_pd(xre,xre);
	                vn    = _mm_set1_pd(n);
	                xmm1  = _mm_mul_pd(xim,xim);
	                r     = _mm_sqrt_pd(_mm_add_pd(xmm0,xmm1));
	                tht   = _mm_atan_pd(_mm_div_pd(xim,xre));
	                pt    = _mm_pow_pd(r,vn);
	                ta    = _mm_mul_pd(vn,tht);
	                _mm_storeu_pd(&ppowr[0] ,_mm_mul_pd(pt,_mm_cos_pd(ta)));
	                _mm_storeu_pd(&ppowi[0] ,_mm_mul_pd(pt,_mm_sin_pd(ta)));                    
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           ZMM16c4_t cpow_xmm2c8(const ZMM16c4_t x,
	                                  const double n) {
	                   
	                ZMM16c4_t cp;        
	                register __m128d xmm0,xmm1;
	                register __m128d r,tht;
	                register __m128d vn,pt;
	                register __m128d ta;
	                xmm0  = _mm_mul_pd(x.re,x.re);
	                vn    = _mm_set1_pd(n);
	                xmm1  = _mm_mul_pd(x.im,x.im);
	                r     = _mm_sqrt_pd(_mm_add_pd(xmm0,xmm1));
	                tht   = _mm_atan_pd(_mm_div_pd(x.im,x.re));
	                pt    = _mm_pow_pd(r,vn);
	                ta    = _mm_mul_pd(vn,tht);
	                cp.re = _mm_mul_pd(pt,_mm_cos_pd(ta));
	                cp.im = _mm_mul_pd(pt,_mm_sin_pd(ta));      
	                return (cp);              
	       }
	       
#include <utility>

                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   ceq_xmm2c8_u(const double * __restrict xre,
                                      const double * __restrict xim,
                                      const double * __restrict yre,
                                      const double * __restrict yim) {
                                      
                      register __m128d xmm0,xmm1,xmm2,xmm3;
                      __mmask8 eqr;
                      __mmask8 eqi;
                      xmm0 = _mm_loadu_pd(&xre[0]);
                      xmm1 = _mm_loadu_pd(&yre[0]);
                      eqr  = _mm_cmp_pd_mask(xmm0,xmm1,_CMP_EQ_OQ);
                      xmm2 = _mm_loadu_pd(&xim[0]);
                      xmm3 = _mm_loadu_pd(&yim[0]);
                      eqi  = _mm_cmp_pd_mask(xmm2,xmm3,_CMP_EQ_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   ceq_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) xre,
                                      const double * __restrict __ATTR_ALIGN__(16) xim,
                                      const double * __restrict __ATTR_ALIGN__(16) yre,
                                      const double * __restrict __ATTR_ALIGN__(16) yim) {
                                      
                      register __m128d xmm0,xmm1,xmm2,xmm3;
                      __mmask8 eqr;
                      __mmask8 eqi;
                      xmm0 = _mm_load_pd(&xre[0]);
                      xmm1 = _mm_load_pd(&yre[0]);
                      eqr  = _mm_cmp_pd_mask(xmm0,xmm1,_CMP_EQ_OQ);
                      xmm2 = _mm_load_pd(&xim[0]);
                      xmm3 = _mm_load_pd(&yim[0]);
                      eqi  = _mm_cmp_pd_mask(xmm2,xmm3,_CMP_EQ_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   ceq_xmm2c8(      const __m128d xre,
                                    const __m128d xim,
                                    const __m128d yre,
                                    const __m128d yim) {
                                    
                          __mmask8 eqr;
                          __mmask8 eqi;
                         eqr = _mm_cmp_pd_mask(xre,yre,_CMP_EQ_OQ);
                         eqi = _mm_cmp_pd_mask(xim,yim,_CMP_EQ_OQ);
                         return std::make_pair(eqr,eqi);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   ceq_xmm2c8(      const xmm2c8_t x,
                                    const xmm2c8_t y) {
                                   
                          __mmask8 eqr;
                          __mmask8 eqi;
                         eqr = _mm_cmp_pd_mask(x.re,y.re,_CMP_EQ_OQ);
                         eqi = _mm_cmp_pd_mask(x.im,y.im,_CMP_EQ_OQ);
                         return std::make_pair(eqr,eqi);
              }
              


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cgt_xmm2c8_u(      const double * __restrict xre,
                                      const double * __restrict xim,
                                      const double * __restrict yre,
                                      const double * __restrict yim) {
                                     
                      register __m128d xmm0,xmm1,xmm2,xmm3;
                      __mmask8 eqr;
                      __mmask8 eqi;
                      xmm0 = _mm_loadu_pd(&xre[0]);
                      xmm1 = _mm_loadu_pd(&yre[0]);
                      eqr  = _mm_cmp_pd_mask(xmm0,xmm1,_CMP_GT_OQ);
                      xmm2 = _mm_loadu_pd(&xim[0]);
                      xmm3 = _mm_loadu_pd(&yim[0]);
                      eqi  = _mm_cmp_pd_mask(xmm2,xmm3,_CMP_GT_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cgt_xmm2c8_a(      const double * __restrict __ATTR_ALIGN__(16) xre,
                                      const double * __restrict __ATTR_ALIGN__(16) xim,
                                      const double * __restrict __ATTR_ALIGN__(16) yre,
                                      const double * __restrict __ATTR_ALIGN__(16) yim) {
                                      
                      register __m128d xmm0,xmm1,xmm2,xmm3;
                      __mmask8 eqr;
                      __mmask8 eqi;
                      xmm0 = _mm_load_pd(&xre[0]);
                      xmm1 = _mm_load_pd(&yre[0]);
                      eqr  = _mm_cmp_pd_mask(xmm0,xmm1,_CMP_GT_OQ);
                      xmm2 = _mm_load_pd(&xim[0]);
                      xmm3 = _mm_load_pd(&yim[0]);
                      eqi  = _mm_cmp_pd_mask(xmm2,xmm3,_CMP_GT_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cgt_xmm2c8(      const __m128d xre,
                                    const __m128d xim,
                                    const __m128d yre,
                                    const __m128d yim) {
                            
                         __mmask8 eqr;
                         __mmask8 eqi;      
                         eqr = _mm_cmp_pd_mask(xre,yre,_CMP_GT_OQ);
                         eqi = _mm_cmp_pd_mask(xim,yim,_CMP_GT_OQ);
                         return std::make_pair(eqr,eqi);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cgt_xmm2c8(      const xmm2c8_t x,
                                    const xmm2c8_t y) {
                           
                         __mmask8 eqr;
                         __mmask8 eqi;         
                         eqr = _mm_cmp_pd_mask(x.re,y.re,_CMP_GT_OQ);
                         eqi = _mm_cmp_pd_mask(x.im,y.im,_CMP_GT_OQ);
                         return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   clt_xmm2c8_u(const double * __restrict xre,
                                      const double * __restrict xim,
                                      const double * __restrict yre,
                                      const double * __restrict yim){
                                      
                      register __m128d xmm0,xmm1,xmm2,xmm3;
                      __mmask8 eqr;
                      __mmask8 eqi; 
                      xmm0 = _mm_loadu_pd(&xre[0]);
                      xmm1 = _mm_loadu_pd(&yre[0]);
                      eqr  = _mm_cmp_pd_mask(xmm0,xmm1,_CMP_LT_OQ);
                      xmm2 = _mm_loadu_pd(&xim[0]);
                      xmm3 = _mm_loadu_pd(&yim[0]);
                      eqi  = _mm_cmp_pd_mask(xmm2,xmm3,_CMP_LT_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   clt_xmm2c8_a(      const double * __restrict __ATTR_ALIGN__(16) xre,
                                      const double * __restrict __ATTR_ALIGN__(16) xim,
                                      const double * __restrict __ATTR_ALIGN__(16) yre,
                                      const double * __restrict __ATTR_ALIGN__(16) yim) {
                                     
                      register __m128d xmm0,xmm1,xmm2,xmm3;
                      __mmask8 eqr;
                      __mmask8 eqi; 
                      xmm0 = _mm_load_pd(&xre[0]);
                      xmm1 = _mm_load_pd(&yre[0]);
                      eqr  = _mm_cmp_pd_mask(xmm0,xmm1,_CMP_LT_OQ);
                      xmm2 = _mm_load_pd(&xim[0]);
                      xmm3 = _mm_load_pd(&yim[0]);
                      eqi  = _mm_cmp_pd_mask(xmm2,xmm3,_CMP_LT_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   clt_xmm2c8(      const __m128d xre,
                                    const __m128d xim,
                                    const __m128d yre,
                                    const __m128d yim) {
                               
                         __mmask8 eqr;
                         __mmask8 eqi;      
                         eqr = _mm_cmp_pd_mask(xre,yre,_CMP_LT_OQ);
                         eqi = _mm_cmp_pd_mask(xim,yim,_CMP_LT_OQ);
                         return std::make_pair(eqr,eqi);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   clt_xmm2c8(     const xmm2c8_t x,
                                    const xmm2c8_t y) {
                                   
                         __mmask8 eqr;
                         __mmask8 eqi; 
                         eqr = _mm_cmp_pd_mask(x.re,y.re,_CMP_LT_OQ);
                         eqi = _mm_cmp_pd_mask(x.im,y.im,_CMP_LT_OQ);
                         return std::make_pair(eqr,eqi);
              }



                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cneq_xmm2c8_u(     const double * __restrict xre,
                                      const double * __restrict xim,
                                      const double * __restrict yre,
                                      const double * __restrict yim) {
                                      
                      register __m128d xmm0,xmm1,xmm2,xmm3;
                      __mmask8 eqr;
                      __mmask8 eqi; 
                      xmm0 = _mm_loadu_pd(&xre[0]);
                      xmm1 = _mm_loadu_pd(&yre[0]);
                      eqr  = _mm_cmp_pd_mask(xmm0,xmm1,_CMP_NEQ_OQ);
                      xmm2 = _mm_loadu_pd(&xim[0]);
                      xmm3 = _mm_loadu_pd(&yim[0]);
                      eqi  = _mm_cmp_pd_mask(xmm2,xmm3,_CMP_NEQ_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cneq_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) xre,
                                      const double * __restrict __ATTR_ALIGN__(16) xim,
                                      const double * __restrict __ATTR_ALIGN__(16) yre,
                                      const double * __restrict __ATTR_ALIGN__(16) yim) {
                                      
                      register __m128d xmm0,xmm1,xmm2,xmm3;
                      __mmask8 eqr;
                      __mmask8 eqi;
                      xmm0 = _mm_load_pd(&xre[0]);
                      xmm1 = _mm_load_pd(&yre[0]);
                      eqr  = _mm_cmp_pd_mask(xmm0,xmm1,_CMP_NEQ_OQ);
                      xmm2 = _mm_load_pd(&xim[0]);
                      xmm3 = _mm_load_pd(&yim[0]);
                      eqi  = _mm_cmp_pd_mask(xmm2,xmm3,_CMP_NEQ_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cneq_xmm2c8(     const __m128d xre,
                                    const __m128d xim,
                                    const __m128d yre,
                                    const __m128d yim) {
                                    
                         __mmask8 eqr;
                         __mmask8 eqi;
                         eqr = _mm_cmp_pd_mask(xre,yre,_CMP_NEQ_OQ);
                         eqi = _mm_cmp_pd_mask(xim,yim,_CMP_NEQ_OQ);
                         return std::make_pair(eqr,eqi);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cneq_xmm2c8(      const xmm2c8_t x,
                                     const xmm2c8_t y) {
                                    
                         __mmask8 eqr;
                         __mmask8 eqi;
                         eqr = _mm_cmp_pd_mask(x.re,y.re,_CMP_NEQ_OQ);
                         eqi = _mm_cmp_pd_mask(x.im,y.im,_CMP_NEQ_OQ);
                         return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cexp_xmm2c8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       double * __restrict cexpr,
                                       double * __restrict cexpi ) {

                        register const __m128d I = _mm_set1_pd(1.0);
                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_loadu_pd(&xre[0]);
                        xmm1  = _mm_loadu_pd(&xim[0]);
                        xmm2  = _mm_exp_pd(xmm0);
                        xmm3  = _mm_mul_pd(xmm2,_mm_cos_pd(xmm0));
                        _mm_storeu_pd(&cexpr[0],xmm3);
                        xmm4  = _mm_mul_pd(xmm2,_mm_mul_pd(_mm_sin_pd(xmm1),I));
                        _mm_storeu_pd(&cexpi[0],xmm4);
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cexp_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) xre,
                                       const double * __restrict __ATTR_ALIGN__(16) xim,
                                       double * __restrict __ATTR_ALIGN__(16) cexpr,
                                       double * __restrict __ATTR_ALIGN__(16) cexpi ) {

                        register const __m128d I = _mm_set1_pd(1.0);
                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_load_pd(&xre[0]);
                        xmm1  = _mm_load_pd(&xim[0]);
                        xmm2  = _mm_exp_pd(xmm0);
                        xmm3  = _mm_mul_pd(xmm2,_mm_cos_pd(xmm0));
                        _mm_store_pd(&cexpr[0],xmm3);
                        xmm4  = _mm_mul_pd(xmm2,_mm_mul_pd(_mm_sin_pd(xmm1),I));
                        _mm_store_pd(&cexpi[0],xmm4);
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cexp_xmm2c8(const __m128d xre,
                                     const __m128d xim,
                                     __m128d * __restrict cexpr,
                                     __m128d * __restrict cexpi) {

                        register const __m128d I = _mm_set1_pd(1.0);
                        register __m128d xmm0;
                        xmm0   = _mm_exp_pd(xre);
                        *cexpr = _mm_mul_pd(xmm0,_mm_cos_pd(xre));
                        *cexpi = _mm_mul_pd(xmm0,_mm_mul_pd(_mm_sin_pd(xim),I));
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cexp_xmm2c8(const __m128d xre,
                                          const __m128d xim) {
                                     
                        xmm2c8_t cv;
                        register const __m128d I = _mm_set1_pd(1.0);
                        register __m128d xmm0;
                        xmm0   = _mm_exp_pd(xre);
                        cv.re = _mm_mul_pd(xmm0,_mm_cos_pd(xre));
                        cv.im = _mm_mul_pd(xmm0,_mm_mul_pd(_mm_sin_pd(xim),I));
                        return (cv);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cexp_xmm2c8(const xmm2c8_t x) {
                                     
                        xmm2c8_t cv;
                        register const __m128d I = _mm_set1_pd(1.0);
                        register __m128d xmm0;
                        xmm0   = _mm_exp_pd(x.re);
                        cv.re = _mm_mul_pd(xmm0,_mm_cos_pd(x.re));
                        cv.im = _mm_mul_pd(xmm0,_mm_mul_pd(_mm_sin_pd(x.im),I));
                        return (cv);
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cpolar_xmm2c8_u(const double * __restrict rho,
                                         const double * __restrict tht,
                                         double * __restrict  re,
                                         double * __restrict  im) {

                         register __m128d xmm0,xmm1,xmm2,xmm3;
                         xmm0 = _mm_loadu_pd(&rho[0]);
                         xmm1 = _mm_loadu_pd(&tht[0]);
                         xmm2 = _mm_mul_pd(xmm0,_mm_cos_pd(xmm1)); //tht
                         _mm_storeu_pd(&re[0],xmm2);
                         xmm3 = _mm_mul_pd(xmm0,_mm_sin_pd(xmm1)); //tht
                         _mm_storeu_pd(&im[0],xmm3);
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cpolar_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) rho,
                                         const double * __restrict __ATTR_ALIGN__(16) tht,
                                         double * __restrict  __ATTR_ALIGN__(16) re,
                                         double * __restrict  __ATTR_ALIGN__(16) im) {

                         register __m128d xmm0,xmm1,xmm2,xmm3;
                         xmm0 = _mm_load_pd(&rho[0]);
                         xmm1 = _mm_load_pd(&tht[0]);
                         xmm2 = _mm_mul_pd(xmm0,_mm_cos_pd(xmm1)); //tht
                         _mm_store_pd(&re[0],xmm2);
                         xmm3 = _mm_mul_pd(xmm0,_mm_sin_pd(xmm1)); //tht
                         _mm_store_pd(&im[0],xmm3);
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cpolar_xmm2c8(const __m128d rho,
                                       const __m128d tht,
                                       __m128d * __restrict re,
                                       __m128d * __restrict im) {

                        register __m128d xmm0,xmm1;
                        xmm0 = _mm_mul_pd(rho,_mm_cos_pd(tht));
                        *re  = xmm0;
                        xmm1 = _mm_mul_pd(rho,_mm_sin_pd(tht));
                        *im  = xmm1;
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cpolar_xmm2c8(const __m128d rho,
                                            const __m128d tht) {
                                      
                        xmm2c8_t cv
                        register __m128d xmm0,xmm1;
                        xmm0 = _mm_mul_pd(rho,_mm_cos_pd(tht));
                        cv.re  = xmm0;
                        xmm1 = _mm_mul_pd(rho,_mm_sin_pd(tht));
                        cv.im  = xmm1;
                        return (cv);
              }
              


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void csqrt_xmm2c8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       double * __restrict wrkc,
                                       double * __restrict csqr,
                                       double * __restrict csqi) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        register __m128d half = _mm_set1_pd(0.5);
                        cabs_xmm2c8_u(xre,xim,wrkc);
                        xmm0  = _mm_loadu_pd(&xre[0]);
                        xmm1  = _mm_loadu_pd(&wrkc[0]);
                        xmm2  = _mm_mul_pd(half,_mm_add_pd(xmm1,xmm0));
                        _mm_storeu_pd(&csqr[0],_mm_sqrt_pd(xmm2));
                        xmm3  = _mm_mul_pd(half,_mm_sub_pd(xmm1,xmm0));
                        _mm_storeu_pd(&csqi[0],_mm_sqrt_pd(xmm3));
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void csqrt_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) xre,
                                       const double * __restrict __ATTR_ALIGN__(16) xim,
                                       double * __restrict __ATTR_ALIGN__(16) wrkc,
                                       double * __restrict __ATTR_ALIGN__(16) csqr,
                                       double * __restrict __ATTR_ALIGN__(16) csqi) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        register __m128d half = _mm_set1_pd(0.5);
                        cabs_xmm2c8_a(xre,xim,wrkc);
                        xmm0  = _mm_load_pd(&xre[0]);
                        xmm1  = _mm_load_pd(&wrkc[0]);
                        xmm2  = _mm_mul_pd(half,_mm_add_pd(xmm1,xmm0));
                        _mm_store_pd(&csqr[0],_mm_sqrt_pd(xmm2));
                        xmm3  = _mm_mul_pd(half,_mm_sub_pd(xmm1,xmm0));
                        _mm_store_pd(&csqi[0],_mm_sqrt_pd(xmm3));
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void csqrt_xmm2c8(const __m128d xre,
                                      const __m128d xim,
                                      __m128d * __restrict wrkc,
                                      __m128d * __restrict csqr,
                                      __m128d * __restrict csqi) {

                       register __m128d xmm0,xmm1;
                       register __m128d half = _mm_set1_pd(0.5); 
                       cabs_xmm2c8(xre,xim,wrkc);
                       xmm0  = _mm_mul_pd(half,_mm_add_pd(*wrkc,xre));
                       *csqr = xmm0;
                       xmm1  = _mm_mul_pd(half,_mm_sub_pd(*wrkc,xre));
                       *csqi = xmm1; 
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t csqrt_xmm2c8(const __m128d xre,
                                           const __m128d xim,
                                          __m128d * __restrict wrkc) {
                                          
                       xmm2c8_t cv;
                       register __m128d xmm0,xmm1;
                       register __m128d half = _mm_set1_pd(0.5); 
                       cabs_xmm2c8(xre,xim,wrkc);
                       xmm0  = _mm_mul_pd(half,_mm_add_pd(*wrkc,xre));
                       cv.re = xmm0;
                       xmm1  = _mm_mul_pd(half,_mm_sub_pd(*wrkc,xre));
                       cv.im = xmm1; 
                       return (cv);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t csqrt_xmm2c8(const xmm2c8_t x,
                                          __m128d * __restrict wrkc) {
                                          
                       xmm2c8_t cv;
                       register __m128d xmm0,xmm1;
                       register __m128d half = _mm_set1_pd(0.5); 
                       cabs_xmm2c8(x.re,x.im,wrkc);
                       xmm0  = _mm_mul_pd(half,_mm_add_pd(*wrkc,x.re));
                       cv.re = xmm0;
                       xmm1  = _mm_mul_pd(half,_mm_sub_pd(*wrkc,x.re));
                       cv.im = xmm1; 
                       return (cv);
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cnorm_prod_xmm2c8_u(const double * __restrict xre,
                                             const double * __restrict xim,
                                             const double * __restrict yre,
                                             const double * __restrict yim,
                                             double * __restrict zre,
                                             double * __restrict zim) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        register __m128d rep,imp;
                        xmm0 = _mm_loadu_pd(&xre[0]);
                        xmm1 = _mm_loadu_pd(&yre[0]);
                        xmm2 = _mm_loadu_pd(&xim[0]);
                        xmm3 = _mm_loadu_pd(&yim[0]);
                        rep  = _mm_fmsub_pd(xmm0,xmm1,
                                               _mm_mul_pd(xmm2,xmm3));
                        imp  = _mm_fmadd_pd(xmm2,xmm1,
                                               _mm_mul_pd(xmm0,xmm3));
                        xmm0 = _mm_mul_pd(rep,rep);
                        xmm1 = _mm_mul_pd(imp,imp);
                        xmm2 = _mm_sqrt_pd(_mm_add_pd(xmm0,xmm1));
                        _mm_storeu_pd(&zre[0], _mm_div_pd(rep,xmm2));
                        _mm_storeu_pd(&zim[0], _mm_div_pd(imp,xmm2));
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cnorm_prod_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) xre,
                                             const double * __restrict __ATTR_ALIGN__(16) xim,
                                             const double * __restrict __ATTR_ALIGN__(16) yre,
                                             const double * __restrict __ATTR_ALIGN__(16) yim,
                                             double * __restrict __ATTR_ALIGN__(16) zre,
                                             double * __restrict __ATTR_ALIGN__(16) zim) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        register __m128d rep,imp;
                        xmm0 = _mm_load_pd(&xre[0]);
                        xmm1 = _mm_load_pd(&yre[0]);
                        xmm2 = _mm_load_pd(&xim[0]);
                        xmm3 = _mm_load_pd(&yim[0]);
                        rep  = _mm_fmsub_pd(xmm0,xmm1,
                                               _mm_mul_pd(xmm2,xmm3));
                        imp  = _mm_fmadd_pd(xmm2,xmm1,
                                               _mm_mul_pd(xmm0,xmm3));
                        xmm0 = _mm_mul_pd(rep,rep);
                        xmm1 = _mm_mul_pd(imp,imp);
                        xmm2 = _mm_sqrt_pd(_mm_add_pd(xmm0,xmm1));
                        _mm_store_pd(&zre[0], _mm_div_pd(rep,xmm2));
                        _mm_store_pd(&zim[0], _mm_div_pd(imp,xmm2));
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cnorm_prod_xmm2c8(  const __m128d  xre,
                                             const __m128d  xim,
                                             const __m128d  yre,
                                             const __m128d  yim,
                                             __m128d * __restrict zre,
                                             __m128d * __restrict zim) {

                        register __m128d rep,imp,xmm0,xmm1,xmm2;
                        rep  = _mm_fmsub_pd(xre,yre,
                                               _mm_mul_pd(xim,yim));
                        imp  = _mm_fmadd_pd(xim,yre,
                                               _mm_mul_pd(xre,yim));
                        xmm0 = _mm_mul_pd(rep,rep);
                        xmm1 = _mm_mul_pd(imp,imp);
                        xmm2 = _mm_sqrt_pd(_mm_add_pd(xmm0,xmm1));
                        *zre = _mm_div_pd(rep,xmm2);
                        *zim = _mm_div_pd(imp,xmm2);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cnorm_prod_xmm2c8(  const __m128d  xre,
                                                  const __m128d  xim,
                                                  const __m128d  yre,
                                                  const __m128d  yim) {
                                             
                        xmm2c8_t cv;
                        register __m128d rep,imp,xmm0,xmm1,xmm2;
                        rep  = _mm_fmsub_pd(xre,yre,
                                               _mm_mul_pd(xim,yim));
                        imp  = _mm_fmadd_pd(xim,yre,
                                               _mm_mul_pd(xre,yim));
                        xmm0 = _mm_mul_pd(rep,rep);
                        xmm1 = _mm_mul_pd(imp,imp);
                        xmm2 = _mm_sqrt_pd(_mm_add_pd(xmm0,xmm1));
                        cv.re = _mm_div_pd(rep,xmm2);
                        cv.im = _mm_div_pd(imp,xmm2);
                        return (cv);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cnorm_prod_xmm2c8(  const xmm2c8_t x,
                                                  const xmm2c8_t y) {
                                             
                        xmm2c8_t cv;
                        register __m128d rep,imp,xmm0,xmm1,xmm2;
                        rep  = _mm_fmsub_pd(x.re,y.re,
                                               _mm_mul_pd(x.im,y.im));
                        imp  = _mm_fmadd_pd(x.im,y.re,
                                               _mm_mul_pd(x.re,y.im));
                        xmm0 = _mm_mul_pd(rep,rep);
                        xmm1 = _mm_mul_pd(imp,imp);
                        xmm2 = _mm_sqrt_pd(_mm_add_pd(xmm0,xmm1));
                        cv.re = _mm_div_pd(rep,xmm2);
                        cv.im = _mm_div_pd(imp,xmm2);
                        return (cv);
             }
             
#include "GMS_simd_utils.hpp"


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmean_prod_xmm2c8_u(const double * __restrict xre,
                                             const double * __restrict xim,
                                             const double * __restrict yre,
                                             const double * __restrict yim,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        register __m128d rep,imp;
                        __m128d t0,t1;
                        constexpr double inv2 = 0.5;
                        double sre,sim;
                        xmm0 = _mm_loadu_pd(&xre[0]);
                        xmm1 = _mm_loadu_pd(&yre[0]);
                        xmm2 = _mm_loadu_pd(&xim[0]);
                        xmm3 = _mm_loadu_pd(&yim[0]);
                        sre = 0.0;
                        rep  = _mm_fmsub_pd(xmm0,xmm1,
                                               _mm_mul_pd(xmm2,xmm3));
                        t0   = _mm_hadd_pd(rep,rep);//xmm8r4_horizontal_sum(rep);
                        sre  = *(double*)&t0;
                        *mre = sre*inv2;
                        sim  = 0.0;
                        imp  = _mm_fmadd_pd(xmm2,xmm1,
                                               _mm_mul_pd(xmm0,xmm3));
                        t1   = _mm_hadd_pd(imp,imp);
                        sim  = *(double*)&t1;
                        *mim = sim*inv2;
              }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmean_prod_xmm2c8_a( const double * __restrict __ATTR_ALIGN__(16) xre,
                                             const double * __restrict __ATTR_ALIGN__(16) xim,
                                             const double * __restrict __ATTR_ALIGN__(16) yre,
                                             const double * __restrict __ATTR_ALIGN__(16) yim,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        register __m128d rep,imp;
                        __m128d t0,t1;
                        constexpr double inv2 = 0.5;
                        double sre,sim;
                        xmm0 = _mm_load_pd(&xre[0]);
                        xmm1 = _mm_load_pd(&yre[0]);
                        xmm2 = _mm_load_pd(&xim[0]);
                        xmm3 = _mm_load_pd(&yim[0]);
                        sre = 0.0;
                        rep  = _mm_fmsub_pd(xmm0,xmm1,
                                               _mm_mul_pd(xmm2,xmm3));
                        t0   = _mm_hadd_pd(rep,rep);
                        sre  = *(double*)&t0;
                        *mre = sre*inv2;
                        sim  = 0.0;
                        imp  = _mm_fmadd_pd(xmm2,xmm1,
                                               _mm_mul_pd(xmm0,xmm3));
                        t1   = _mm_hadd_pd(imp,imp);
                        sim  = *(double*)&t1;
                        *mim = sim*inv2;
              } 


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmean_prod_xmm2c8( const __m128d xre,
                                           const __m128d xim,
                                           const __m128d yre,
                                           const __m128d yim,
                                           double * __restrict mre,
                                           double * __restrict min) {

                        register __m128d rep,imp;
                        __m128d t0,t1;
                        constexpr double inv16 = 0.5;
                        double sre,sim;
                        sre = 0.0;
                        rep  = _mm_fmsub_pd(xre,yre,
                                               _mm_mul_pd(xim,yim));
                        t0   = _mm_hadd_pd(rep,rep);
                        sre  = *(double*)&t0;
                        *mre = sre*inv2;
                        sim  = 0.0;
                        imp  = _mm_fmadd_pd(xim,yre,
                                               _mm_mul_pd(xre,yim));
                        t1   = _mm_hadd_pd(imp,imp);
                        sim  = *(double*)&t1;
                        *mim = sim*inv2;
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmean_prod_xmm2c8(const xmm2c8_t x,
                                           const xmm2c8_t y,
                                           double * __restrict mre,
                                           double * __restrict min) {

                        register __m128d rep,imp;
                        __m128d t0,t1
                        constexpr double inv16 = 0.5;
                        double sre,sim;
                        sre = 0.0;
                        rep  = _mm_fmsub_pd(x.re,y.re,
                                               _mm_mul_pd(x.im,y.im));
                        t0   = _mm_hadd_pd(rep,rep);
                        sre  = *(double*)&t0;
                        *mre = sre*inv2;
                        sim  = 0.0;
                        imp  = _mm_fmadd_pd(x.im,y.re,
                                               _mm_mul_pd(x.re,y.im));
                        t1   = _mm_hadd_pd(imp,imp);
                        sim  = *(double*)&t1;
                        *mim = sim*inv2;
             }


               
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmean_quot_xmm2c8_u(const double * __restrict xre,
                                             const double * __restrict xim,
                                             const double * __restrict yre,
                                             const double * __restrict yim,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        register __m128d rep,imp,den,rquot,iquot;
                        __m128d t0,t1;
                        constexpr double inv2 = 0.5;
                        double sre,sim;
                        xmm0 = _mm_loadu_pd(&xre[0]);
                        xmm1 = _mm_loadu_pd(&yre[0]);
                        xmm2 = _mm_loadu_pd(&xim[0]);
                        xmm3 = _mm_loadu_pd(&yim[0]);
                        sre  = 0.0;
                        rep  = _mm_fmsub_pd(xmm0,xmm1,
                                               _mm_mul_pd(xmm2,xmm3));
                        imp  = _mm_fmadd_pd(xmm2,xmm1,
                                               _mm_mul_pd(xmm0,xmm3));
                        sim  = 0.0;
                        den  = _mm_fmadd_pd(xmm1,xmm1,
                                               _mm_mul_pd(xmm3,xmm3));
                        rquot = _mm_div_pd(rep,den);
                        t0    = _mm_hadd_pd(rquot,rquot);
                        sre   = *(double*)&t0;
                        *mre  = sre*inv2;
                        iquot = _mm_div_pd(imp,den);
                        t1    = _mm_hadd_pd(iquot,iquot);
                        sim   = *(double*)&t1;
                        *mim  = sre*inv2;
              }  


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmean_quot_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) xre,
                                             const double * __restrict __ATTR_ALIGN__(16) xim,
                                             const double * __restrict __ATTR_ALIGN__(16) yre,
                                             const double * __restrict __ATTR_ALIGN__(16) yim,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        register __m128d rep,imp,den,rquot,iquot;
                        __m128d t0,t1;
                        constexpr double inv2 = 0.5;
                        double sre,sim;
                        xmm0 = _mm_load_pd(&xre[0]);
                        xmm1 = _mm_load_pd(&yre[0]);
                        xmm2 = _mm_load_pd(&xim[0]);
                        xmm3 = _mm_load_pd(&yim[0]);
                        sre  = 0.0f;
                        rep  = _mm_fmsub_pd(xmm0,xmm1,
                                               _mm_mul_pd(xmm2,xmm3));
                        imp  = _mm_fmadd_pd(xmm2,xmm1,
                                               _mm_mul_pd(xmm0,xmm3));
                        sim  = 0.0f;
                        den  = _mm_fmadd_pd(xmm1,xmm1,
                                               _mm_mul_pd(xmm3,xmm3));
                        rquot = _mm_div_pd(rep,den);
                        t0    = _mm_hadd_pd(rquot,rquot);
                        sre   = *(double*)&t0;
                        *mre  = sre*inv2;
                        iquot = _mm_div_pd(imp,den);
                        t1    = _mm_hadd_pd(iquot,iquot);
                        sim   = *(double*)&t1;
                        *mim  = sre*inv2;
              }  


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmean_quot_xmm2c8(  const __m128d xre,
                                             const __m128d xim,
                                             const __m128d yre,
                                             const __m128d yim,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m128d rep,imp,den,rquot,iquot;
                        __m128d t0,t1;
                        constexpr double inv2 = 0.5;
                        double sre,sim;
                        sre  = 0.0;
                        rep  = _mm_fmsub_pd(xre,yre,
                                               _mm_mul_pd(xim,yim));
                        imp  = _mm_fmadd_pd(xim,yre,
                                               _mm_mul_pd(xre,yim));
                        sim  = 0.0;
                        den  = _mm_fmadd_pd(yre,yre,
                                               _mm_mul_pd(yim,yim));
                        rquot = _mm_div_pd(rep,den);
                        t0    = _mm_hadd_pd(rquot,rquot);
                        sre   = *(double*)&t0;
                        *mre  = sre*inv2;
                        iquot = _mm_div_pd(imp,den);
                        t1    = _mm_hadd_pd(iquot,iquot);
                        sim   = *(double*)&t1;
                        *mim  = sre*inv2;
              }  
              
              
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmean_quot_xmm2c8(  const xmm2c8_t x,
                                             const xmm2c8_t y,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m128d rep,imp,den,rquot,iquot;
                        __m128 t0,t1;
                        constexpr double inv2 = 0.5;
                        double sre,sim;
                        sre  = 0.0;
                        rep  = _mm_fmsub_pd(x.re,y.re,
                                               _mm_mul_pd(x.im,y.im));
                        imp  = _mm_fmadd_pd(xim,yre,
                                               _mm_mul_pd(x.re,y.im));
                        sim  = 0.0;
                        den  = _mm_fmadd_pd(y.re,y.re,
                                               _mm_mul_pd(y.im,y.im));
                        rquot = _mm_div_pd(rep,den);
                        t0    = _mm_hadd_pd(rquot,rquot);
                        sre   = *(double*)&t0;
                        *mre  = sre*inv2;
                        iquot = _mm_div_pd(imp,den);
                        t1    = _mm_hadd_pd(iquot,iquot);
                        sim   = *(double*)&t1;
                        *mim  = sre*inv2;
              }  


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cnorm_cprod_xmm2c8_u(const double * __restrict xre,
                                              const double * __restrict xim,
                                              const double * __restrict yre,
                                              const double * __restrict yim,
                                              double * __restrict mre,
                                              double * __restrict mim ) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        register __m128d rep,imp,magc1,magc2,vcmag;
                        xmm0 = _mm_loadu_pd(&xre[0]);
                        xmm1 = _mm_loadu_pd(&yre[0]);
                        xmm2 = _mm_loadu_pd(&xim[0]);
                        xmm3 = _mm_loadu_pd(&yim[0]);
                        rep  = _mm_fmadd_pd(xmm0,xmm1,
                                               _mm_mul_pd(xmm2,xmm3));
                        magc1= _mm_mul_pd(rep,rep);
                        imp  = _mm_fmsub_pd(xmm2,xmm1,
                                               _mm_mul_pd(xmm0,xmm3));
                        magc2= _mm_mul_pd(imp,imp);
                        vcmag= _mm_sqrt_pd(_mm_add_pd(magc1,magc2));
                        _mm_storeu_pd(&mre[0], _mm_div_pd(rep,vcmag));
                        _mm_storeu_pd(&mim[0], _mm_div_pd(imp,vcmag));
             }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cnorm_cprod_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) xre,
                                              const double * __restrict __ATTR_ALIGN__(16) xim,
                                              const double * __restrict __ATTR_ALIGN__(16) yre,
                                              const double * __restrict __ATTR_ALIGN__(16) yim,
                                              double * __restrict __ATTR_ALIGN__(16) mre,
                                              double * __restrict __ATTR_ALIGN__(16) mim ) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        register __m128d rep,imp,magc1,magc2,vcmag;
                        xmm0 = _mm_load_pd(&xre[0]);
                        xmm1 = _mm_load_pd(&yre[0]);
                        xmm2 = _mm_load_pd(&xim[0]);
                        xmm3 = _mm_load_pd(&yim[0]);
                        rep  = _mm_fmad_pd(xmm0,xmm1,
                                               _mm_mul_pd(xmm2,xmm3));
                        magc1= _mm_mul_pd(rep,rep);
                        imp  = _mm_fmsub_pd(xmm2,xmm1,
                                               _mm_mul_pd(xmm0,xmm3));
                        magc2= _mm_mul_pd(imp,imp);
                        vcmag= _mm_sqrt_pd(_mm_add_pd(magc1,magc2));
                        _mm_store_pd(&mre[0], _mm_div_pd(rep,vcmag));
                        _mm_store_pd(&mim[0], _mm_div_pd(imp,vcmag));
             }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cnorm_cprod_xmm2c8(const __m128d xre,
                                            const __m128d xim,
                                            const __m128d yre,
                                            const __m128d yim,
                                            __m128d * __restrict mre,
                                            __m128d * __restrict mim) {

                        register __m128d rep,imp,magc1,magc2,vcmag;
                        rep  = _mm_fmad_pd(xre,yre,
                                               _mm_mul_pd(xim,yim));
                        magc1= _mm_mul_pd(rep,rep);
                        imp  = _mm_fmsub_pd(xim,yre,
                                               _mm_mul_pd(xre,yim));
                        magc2= _mm_mul_pd(imp,imp);
                        vcmag= _mm_sqrt_pd(_mm_add_pd(magc1,magc2));
                        *mre = _mm_div_pd(rep,vcmag);
                        *mim = _mm_div_pd(imp,vcmag);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cnorm_cprod_xmm2c8(const xmm2c8_t x,
                                            const xmm2c8_t y,
                                            __m128d * __restrict mre,
                                            __m128d * __restrict mim) {

                        register __m128d rep,imp,magc1,magc2,vcmag;
                        rep  = _mm_fmad_pd(x.re,y.re,
                                               _mm_mul_pd(x.im,y.im));
                        magc1= _mm_mul_pd(rep,rep);
                        imp  = _mm_fmsub_pd(x.im,y.re,
                                               _mm_mul_pd(x.re,y.im));
                        magc2= _mm_mul_pd(imp,imp);
                        vcmag= _mm_sqrt_pd(_mm_add_pd(magc1,magc2));
                        *mre = _mm_div_pd(rep,vcmag);
                        *mim = _mm_div_pd(imp,vcmag);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cnorm_cprod_xmm2c8(const __m128d xre,
                                                 const __m128d xim,
                                                 const __m128d yre,
                                                 const __m128d yim) {
                                               
                        xmm2c8_t cv;
                        register __m128d rep,imp,magc1,magc2,vcmag;
                        rep  = _mm_fmad_pd(xre,yre,
                                               _mm_mul_pd(xim,yim));
                        magc1= _mm_mul_pd(rep,rep);
                        imp  = _mm_fmsub_pd(xim,yre,
                                               _mm_mul_pd(xre,yim));
                        magc2= _mm_mul_pd(imp,imp);
                        vcmag= _mm_sqrt_pd(_mm_add_pd(magc1,magc2));
                        cv.re = _mm_div_pd(rep,vcmag);
                        cv.im = _mm_div_pd(imp,vcmag);
                        return (cv);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cnorm_cprod_xmm2c8(const xmm2c8_t x,
                                                 const xmm2c8_t y) {
                                               
                        xmm2c8_t cv;
                        register __m128d rep,imp,magc1,magc2,vcmag;
                        rep  = _mm_fmad_pd(x.re,y.re,
                                               _mm_mul_pd(x.im,y.im));
                        magc1= _mm_mul_pd(rep,rep);
                        imp  = _mm_fmsub_pd(x.im,y.re,
                                               _mm_mul_pd(x.re,y.im));
                        magc2= _mm_mul_pd(imp,imp);
                        vcmag= _mm_sqrt_pd(_mm_add_pd(magc1,magc2));
                        cv.re = _mm_div_pd(rep,vcmag);
                        cv.im = _mm_div_pd(imp,vcmag);
                        return (cv);
             }
             
             
             


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmean_cprod_xmm2c8_u(const double * __restrict xre,
                                              const double * __restrict xim,
                                              const double * __restrict yre,
                                              const double * __restrict yim,
                                              double * __restrict mre,
                                              double * __restrict mim) {
                      
                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        register __m128d re,im;
                        __m128d t0,t1;
                        constexpr double inv2 = 0.5;
                        double sre,sim;
                        xmm0 = _mm_loadu_pd(&xre[0]);
                        xmm1 = _mm_loadu_pd(&yre[0]);
                        xmm2 = _mm_loadu_pd(&xim[0]);
                        xmm3 = _mm_loadu_pd(&yim[0]);
                        re   = _mm_fmadd_pd(xmm0,xmm1,
                                               _mm_mul_pd(xmm2,xmm3));
                        t0   = _mm_hadd_pd(re,re);
                        sre  = *(double*)&t0;
                        *mre = sre*inv2;
                        im   = _mm_fmsub_pd(xmm2,xmm1,
                                               _mm_mul_pd(xmm0,xmm3));
                        t1   = _mm_hadd_pd(im,im);
                        sim  = *(double*)&t1;
                        *mim = sim*inv2;
             }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmean_cprod_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) xre,
                                              const double * __restrict __ATTR_ALIGN__(16) xim,
                                              const double * __restrict __ATTR_ALIGN__(16) yre,
                                              const double * __restrict __ATTR_ALIGN__(16) yim,
                                              double * __restrict mre,
                                              double * __restrict mim) {
                      
                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        register __m128d re,im;
                        __m128d t0,t1;
                        constexpr double inv2 = 0.5;
                        double sre,sim;
                        xmm0 = _mm_load_pd(&xre[0]);
                        xmm1 = _mm_load_pd(&yre[0]);
                        xmm2 = _mm_load_pd(&xim[0]);
                        xmm3 = _mm_load_pd(&yim[0]);
                        re   = _mm_fmadd_pd(xmm0,xmm1,
                                               _mm_mul_pd(xmm2,xmm3));
                        t0   = _mm_hadd_pd(re,re);
                        sre  = *(double*)&t0;
                        *mre = sre*inv2;
                        im   = _mm_fmsub_pd(xmm2,xmm1,
                                               _mm_mul_pd(xmm0,xmm3));
                        t1   = _mm_hadd_pd(im,im);
                        sim  = *(double*)&t1;
                        *mim = sim*inv2;
             }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmean_cprod_xmm2c8( const __m128d xre,
                                            const __m128d xim,
                                            const __m128d yre,
                                            const __m128d yim,
                                            double * __restrict mre,
                                            double * __restrict min) {

                        register __m128d re,im;
                        __m128d t0,t1;
                        constexpr double inv2 = 0.5;
                        double sre,sim;
                        re   = _mm_fmadd_pd(xre,yre,
                                               _mm_mul_pd(xim,yim));
                        t0   = _mm_hadd_pd(re,re);
                        sre  = *(double*)&t0;
                        *mre = sre*inv2;
                        im   = _mm_fmsub_pd(xim,yre,
                                               _mm_mul_pd(xre,yim));
                        t1   = _mm_hadd_pd(im,im);
                        sim  = *(double*)&t1;
                        *mim = sim*inv2;
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmean_cprod_xmm2c8(const xmm2c8_t x,
                                            const xmm2c8_t y,
                                            double * __restrict mre,
                                            double * __restrict min) {

                        register __m128d re,im;
                        __m128d t0,t1;
                        constexpr double inv2 = 0.5;
                        double sre,sim;
                        re   = _mm_fmadd_pd(x.re,y.re,
                                               _mm_mul_pd(x.im,y.im));
                        t0   = _mm_hadd_pd(re,re);
                        sre  = *(double*)&t0;
                        *mre = sre*inv2;
                        im   = _mm_fmsub_pd(x.im,y.re,
                                               _mm_mul_pd(x.re,y.im));
                        t1   = _mm_hadd_pd(im,im);
                        sim  = *(double*)&t1;
                        *mim = sim*inv2;
             }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void arith_cmean_xmm2c8_u(const double * __restrict xre,
                                              const double * __restrict xim,
                                              double * __restrict mre,
                                              double * __restrict min) {

                        register __m128d re,im;
                        __m128d t0,t1;
                        constexpr double inv16 = 0.5;
                        double sre,sim;
                        re   = _mm_loadu_pd(&xre[0]);
                        t0   = _mm_hadd_pd(re,re);
                        sre  = *(double*)&t0;
                        *mre = sre*inv2;
                        im   = _mm_loadu_pd(&xim[0]);
                        t1   = _mm_hadd_pd(im,im);
                        sim  = *(double*)&t1;
                        *mim = sim*inv2; 
             }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void arith_cmean_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) xre,
                                              const double * __restrict __ATTR_ALIGN__(16) xim,
                                              double * __restrict mre,
                                              double * __restrict min) {

                        register __m128d re,im;
                        __m128d t0,t1;
                        constexpr double inv16 = 0.5;
                        double sre,sim;
                        re   = _mm_load_pd(&xre[0]);
                        t0   = _mm_hadd_pd(re,re);
                        sre  = *(double*)&t0;
                        *mre = sre*inv2;
                        im   = _mm_load_pd(&xim[0]);
                        t1   = _mm_hadd_pd(im,im);
                        sim  = *(double*)&t1;
                        *mim = sim*inv2; 
             }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void arith_cmean_xmm2c8(  const __m128d xre,
                                              const __m128d xim,
                                              double * __restrict mre,
                                              double * __restrict min) {

                        constexpr double inv2 = 0.5;
                        __m128d t0,t1;
                        double sre,sim;
                        t0   = _mm_hadd_pd(xre,xre);
                        sre  = *(double*)&t0;
                        *mre = sre*inv2;
                        t1   = _mm_hadd_pd(xim,xim);
                        sim  = *(double*)&t1;
                        *mim = sim*inv2; 
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void arith_cmean_xmm2c8(  const xmm2c8_t x,
                                              double * __restrict mre,
                                              double * __restrict min) {

                        constexpr double inv2 = 0.5;
                        __m128d t0,t1;
                        double sre,sim;
                        t0   = _mm_hadd_pd(x.re,x.re);
                        sre  = *(double*)&t0;
                        *mre = sre*inv2;
                        t1   = _mm_hadd_pd(x.im,x.im);
                        sim  = *(double*)&t1;
                        *mim = sim*inv2; 
             }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cnormalize_xmm2c8_u( const double * __restrict xre,
                                              const double * __restrict xim,
                                              const double * __restrict yre,
                                              const double * __restrict yim,
                                              double * __restrict mre,
                                              double * __restrict mim ) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        register __m128d re,im,cvmag;
                        xmm0 = _mm_loadu_pd(&xre[0]);
                        xmm1 = _mm_loadu_pd(&yre[0]);
                        xmm2 = _mm_loadu_pd(&xim[0]);
                        xmm3 = _mm_loadu_pd(&yim[0]);
                        cvmag= _mm_sqrt_pd(_mm_fmadd_pd(xmm0,xmm1,
                                                              _mm_mul_pd(xmm2,xmm3)));
                        _mm_storeu_pd(&mre[0], _mm_div_pd(xmm0,cvmag));
                        _mm_storeu_pd(&mim[0], _mm_div_pd(xmm2,cvmag));
             }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cnormalize_xmm2c8_a( const double * __restrict __ATTR_ALIGN__(16) xre,
                                              const double * __restrict __ATTR_ALIGN__(16) xim,
                                              const double * __restrict __ATTR_ALIGN__(16) yre,
                                              const double * __restrict __ATTR_ALIGN__(16) yim,
                                              double * __restrict __ATTR_ALIGN__(16) mre,
                                              double * __restrict __ATTR_ALIGN__(16) mim ) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        register __m128d re,im,cvmag;
                        xmm0 = _mm_load_pd(&xre[0]);
                        xmm1 = _mm_load_pd(&yre[0]);
                        xmm2 = _mm_load_pd(&xim[0]);
                        xmm3 = _mm_load_pd(&yim[0]);
                        cvmag= _mm_sqrt_pd(_mm_fmadd_pd(xmm0,xmm1,
                                                              _mm_mul_pd(xmm2,xmm3)));
                        _mm_store_pd(&mre[0], _mm_div_pd(xmm0,cvmag));
                        _mm_store_pd(&mim[0], _mm_div_pd(xmm2,cvmag));
             }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cnormalize_xmm2c8( const __m128d xre,
                                            const __m128d xim,
                                            const __m128d yre,
                                            const __m128d yim,
                                            __m128d * __restrict mre,
                                            __m128d * __restrict mim ) {

                        register __m128d re,im,cvmag;
                        cvmag= _mm_sqrt_pd(_mm_fmadd_pd(xre,yre,
                                                    _mm_mul_pd(xim,yim)));
                        *mre = _mm_div_pd(xre,cvmag));
                        *mim =  _mm_div_pd(xim,cvmag));
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cnormalize_xmm2c8( const xmm2c8_t x,
                                            const xmm2c8_t y,
                                            __m128d * __restrict mre,
                                            __m128d * __restrict mim ) {

                        register __m128d re,im,cvmag;
                        cvmag= _mm_sqrt_pd(_mm_fmadd_pd(x.re,y.re,
                                                    _mm_mul_pd(x.im,y.im)));
                        *mre = _mm_div_pd(x.re,cvmag));
                        *mim =  _mm_div_pd(x.im,cvmag));
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cnormalize_xmm2c8( const __m128d xre,
                                                 const __m128d xim,
                                                 const __m128d yre,
                                                 const __m128d yim) {
                                            
                        xmm2c8_t cv;
                        register __m128d re,im,cvmag;
                        cvmag= _mm_sqrt_pd(_mm_fmadd_pd(xre,yre,
                                                    _mm_mul_pd(xim,yim)));
                        cv.re = _mm_div_pd(xre,cvmag));
                        cv.im =  _mm_div_pd(xim,cvmag));
                        return (cv);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   xmm2c8_t cnormalize_xmm2c8( const xmm2c8_t x,
                                                 const xmm2c8_t y,) {
                                            
                        xmm2c8_t cv;
                        register __m128d re,im,cvmag;
                        cvmag= _mm_sqrt_pd(_mm_fmadd_pd(x.re,y.re,
                                                    _mm_mul_pd(x.im,y.im)));
                        cv.re = _mm_div_pd(x.re,cvmag));
                        cv.im =  _mm_div_pd(x.im,cvmag));
                        return (cv);
             }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmagnitude_xmm2c8_u( const double * __restrict xre,
                                              const double * __restrict xim,
                                              const double * __restrict yre,
                                              const double * __restrict yim,
                                              double * __restrict mre) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        register __m128d cvmag;
                        xmm0 = _mm_loadu_pd(&xre[0]);
                        xmm1 = _mm_loadu_pd(&yre[0]);
                        xmm2 = _mm_loadu_pd(&xim[0]);
                        xmm3 = _mm_loadu_pd(&yim[0]);
                        cvmag= _mm_sqrt_pd(_mm_fmadd_pd(xmm0,xmm1,
                                                          _mm_mul_pd(xmm2,xmm3)));
                        _mm_storeu_pd(&mre[0], cvmag);
             }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmagnitude_xmm2c8_a( const double * __restrict __ATTR_ALIGN__(16) xre,
                                              const double * __restrict __ATTR_ALIGN__(16) xim,
                                              const double * __restrict __ATTR_ALIGN__(16) yre,
                                              const double * __restrict __ATTR_ALIGN__(16) yim,
                                              double * __restrict __ATTR_ALIGN__(16) mre) {

                        register __m128d xmm0,xmm1,xmm2,xmm3;
                        register __m128d cvmag;
                        xmm0 = _mm_load_pd(&xre[0]);
                        xmm1 = _mm_load_pd(&yre[0]);
                        xmm2 = _mm_load_pd(&xim[0]);
                        xmm3 = _mm_load_pd(&yim[0]);
                        cvmag= _mm_sqrt_pd(_mm_fmadd_pd(xmm0,xmm1,
                                                          _mm_mul_pd(xmm2,xmm3)));
                        _mm_store_pd(&mre[0], cvmag);
             }


                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmagnitude_xmm2c8(   const __m128d xre,
                                              const __m128d xim,
                                              const __m128d yre,
                                              const __m128d yim,
                                              __m128d * __restrict  mre) {

                        register __m128d cvmag;
                        cvmag= _mm_sqrt_pd(_mm_fmadd_pd(xre,yre,
                                                          _mm_mul_pd(xim,yim)));
                        *mre = cvmag;
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	        
	          
                  
	           static inline
                   void cmagnitude_xmm2c8(   const xmm2c8_t x,
                                              const xmm2c8_t y,
                                              __m128d * __restrict  mre) {

                        register __m128d cvmag;
                        cvmag= _mm_sqrt_pd(_mm_fmadd_pd(x.re,y.re,
                                                          _mm_mul_pd(x.im,y.im)));
                        *mre = cvmag;
             }


#include <complex>


                   __ATTR_ALWAYS_INLINE__
	        
	           
                   static inline
                   void copy_2xr4_c4_unroll16x(double * __restrict __ATTR_ALIGN__(64) xre,
                                               double * __restrict __ATTR_ALIGN__(64) xim,
                                               std::complex<double> * __restrict __ATTR_ALIGN__(64) vc,
                                               const int32_t n) {

                        if(__builtin_expect(0<=n,0)) { return;}
                        register double re0,im0,re1,im1,re2,im2;
                        register double re3,im3,re4,im4,re5,im5;
                        register double re6,im6,re7,im7;
                        int32_t j,m,m1;
                        m = n%16;
#if defined(__INTEL_COMPILER)
                        __assume_aligned(xre,64);
                        __assume_aligned(xim,64);
                        __assume_aligned(c,64);
#elif defined(__GNUC__) && (!defined __INTEL_COMPILER)
                        xre = (double*)__builtin_assume_aligned(xre,64);
                        xim = (double*)__builtin_assume_aligned(xim,);
                        vc   = (std::complex<double>*)__builtin_assume_aligned(vc,64);
#endif
                        if(n != 0) {
#pragma omp simd aligned(xre:64) aligned(xim:64) aligned(c:64)
                           for(j = 0; j != m; ++j) {
                               const double re = xre[j];
                               const double im = xim[j];
                               const std::complex<double> tj{re,im};
                               vc[i] = tj;
                           }
                            
                            if(n<16) {return;}
                       }

                       m1 = m+1;
#if defined(__INTEL_COMPILER)
#pragma code_align(32)
#endif
#pragma omp simd aligned(xre:64) aligned(xim:64) aligned(vc:64)
                       for(j = m1; j != n; j += 16) {
                           re0     = xre[i+0];
                           im0     = xim[i+0];
                           vc[i+0] = {re0,im0};
                           re1     = xre[i+1];
                           im1     = xim[i+1];
                           vc[i+1] = {re1,im1}; 
                           re2     = xre[i+2];
                           im2     = xim[i+2];
                           vc[i+2] = {re2,im2};
                           re3     = xre[i+3];
                           im3     = xim[i+3];
                           vc[i+3] = {re3,im3};
                           re4     = xre[i+4];
                           im4     = xim[i+4];
                           vc[i+4] = {re4,im4};
                           re5     = xre[i+5];
                           im5     = xim[i+5];
                           vc[i+5] = {re5,im5}; 
                           re6     = xre[i+6];
                           im6     = xim[i+6];
                           vc[i+6] = {re6,im6};
                           re7     = xre[i+7];
                           im7     = xim[i+7];
                           vc[i+7] = {re7,im7};
                           re0     = xre[i+8];
                           im0     = xim[i+8];
                           vc[i+8] = {re0,im0};
                           re1     = xre[i+9];
                           im1     = xim[i+9];
                           vc[i+9] = {re1,im1}; 
                           re2     = xre[i+10];
                           im2     = xim[i+10];
                           vc[i+10]= {re2,im2};
                           re3     = xre[i+11];
                           im3     = xim[i+11];
                           vc[i+11]= {re3,im3};
                           re4     = xre[i+12];
                           im4     = xim[i+12];
                           vc[i+12]= {re4,im4};
                           re5     = xre[i+13];
                           im5     = xim[i+13];
                           vc[i+13]= {re5,im5}; 
                           re6     = xre[i+14];
                           im6     = xim[i+14];
                           vc[i+14]= {re6,im6};
                           re7     = xre[i+15];
                           im7     = xim[i+15];
                           vc[i+15]= {re7,im7};
                      }
                 }


                /*   __ATTR_ALWAYS_INLINE__
	        
	           
                   static inline
                   void copy_2xr8_c8_unroll16x(double * __restrict xre,
                                               double * __restrict xim,
                                               std::complex<double> * __restrict vc,
                                               const int32_t n) {

                        if(__builtin_expect(0<=n,0)) { return;}
                        register double re0,im0,re1,im1,re2,im2;
                        register double re3,im3,re4,im4,re5,im5;
                        register double re6,im6,re7,im7;
                        int32_t j,m,m1;
                        m = n%16;
#if defined(__INTEL_COMPILER)
                        __assume_aligned(xre,64);
                        __assume_aligned(xim,64);
                        __assume_aligned(c,64);
#elif defined(__GNUC__) && (!defined __INTEL_COMPILER)
                        xre = (double*)__builtin_assume_aligned(xre,64);
                        xim = (double*)__builtin_assume_aligned(xim,64);
                        vc  = (std::complex<double>*)__builtin_assume_aligned(vc,64);
#endif
                        if(n != 0) {
#pragma omp simd aligned(xre:64) aligned(xim:64) aligned(c:64)
                           for(j = 0; j != m; ++j) {
                               const double re = xre[j];
                               const double im = xim[j];
                               const std::complex<double> tj{re,im};
                               vc[i] = tj;
                           }
                            
                            if(n<16) {return;}
                       }

                       m1 = m+1;
#if defined(__INTEL_COMPILER)
#pragma code_align(32)
#endif
#pragma omp simd aligned(xre:64) aligned(xim:64) aligned(vc:64)
                       for(j = m1; j != n; j += 16) {
                           re0     = xre[i+0];
                           im0     = xim[i+0];
                           vc[i+0] = {re0,im0};
                           re1     = xre[i+1];
                           im1     = xim[i+1];
                           vc[i+1] = {re1,im1}; 
                           re2     = xre[i+2];
                           im2     = xim[i+2];
                           vc[i+2] = {re2,im2};
                           re3     = xre[i+3];
                           im3     = xim[i+3];
                           vc[i+3] = {re3,im3};
                           re4     = xre[i+4];
                           im4     = xim[i+4];
                           vc[i+4] = {re4,im4};
                           re5     = xre[i+5];
                           im5     = xim[i+5];
                           vc[i+5] = {re5,im5}; 
                           re6     = xre[i+6];
                           im6     = xim[i+6];
                           vc[i+6] = {re6,im6};
                           re7     = xre[i+7];
                           im7     = xim[i+7];
                           vc[i+7] = {re7,im7};
                           re0     = xre[i+8];
                           im0     = xim[i+8];
                           vc[i+8] = {re0,im0};
                           re1     = xre[i+9];
                           im1     = xim[i+9];
                           vc[i+9] = {re1,im1}; 
                           re2     = xre[i+10];
                           im2     = xim[i+10];
                           vc[i+10]= {re2,im2};
                           re3     = xre[i+11];
                           im3     = xim[i+11];
                           vc[i+11]= {re3,im3};
                           re4     = xre[i+12];
                           im4     = xim[i+12];
                           vc[i+12]= {re4,im4};
                           re5     = xre[i+13];
                           im5     = xim[i+13];
                           vc[i+13]= {re5,im5}; 
                           re6     = xre[i+14];
                           im6     = xim[i+14];
                           vc[i+14]= {re6,im6};
                           re7     = xre[i+15];
                           im7     = xim[i+15];
                           vc[i+15]= {re7,im7};
                      }
                 }

*/
                   __ATTR_ALWAYS_INLINE__
	        
	           
                   static inline
                   void copy_c4_2xr4_unroll16x( double * __restrict __ATTR_ALIGN__(64) xre,
                                               double * __restrict __ATTR_ALIGN__(64) xim,
                                               std::complex<double> * __restrict __ATTR_ALIGN__(64) vc,
                                               const int32_t n) {

                        if(__builtin_expect(0<=n,0)) { return;}
                        register double re0,im0,re1,im1,re2,im2;
                        register double re3,im3,re4,im4,re5,im5;
                        register double re6,im6,re7,im7;
                        int32_t j,m,m1;
                        m = n%16;
#if defined(__INTEL_COMPILER)
                        __assume_aligned(xre,64);
                        __assume_aligned(xim,64);
                        __assume_aligned(c,64);
#elif defined(__GNUC__) && (!defined __INTEL_COMPILER)
                        xre = (double*)__builtin_assume_aligned(xre,64);
                        xim = (double*)__builtin_assume_aligned(xim,64);
                        vc   = (std::complex<double>*)__builtin_assume_aligned(vc,64);
#endif
                        if(n != 0) {
#pragma omp simd aligned(xre:64) aligned(xim:64) aligned(c:64)
                           for(j = 0; j != m; ++j) {
                               const double re = vc[j].real();
                               const double im = vc[j].imag();
                               xre[j]         = re;
                               xim[j]         = im;
                           }
                            
                            if(n<16) {return;}
                       }

                       m1 = m+1;
#if defined(__INTEL_COMPILER)
#pragma code_align(32)
#endif
#pragma omp simd aligned(xre:64) aligned(xim:64) aligned(vc:64)
                       for(j = m1; j != n; j += 16) {
                           re0    = vc[j+0].real();
                           xre[j+0] = re0;
                           im0    = vc[j+0].imag();
                           xim[j+0] = im0;
                           re1    = vc[j+1].real();
                           xre[j+1] = re1;
                           im1    = vc[j+1].imag();
                           xim[j+1] = im1;
                           re2    = vc[j+2].real();
                           xre[j+2] = re2;
                           im2    = vc[j+2].imag();
                           xim[j+2] = im2;
                           re3    = vc[j+3].real();
                           xre[j+3] = re3;
                           im3    = vc[j+3].imag();
                           xim[j+3] = im3;
                           re4    = vc[j+4].real();
                           xre[j+4] = re4;
                           im4    = vc[j+4].imag();
                           xim[j+4] = im4;
                           re5    = vc[j+5].real();
                           xre[j+5] = re5;
                           im5    = vc[j+5].imag();
                           xim[j+5] = im5;
                           re6    = vc[j+6].real();
                           xre[j+6] = re6;
                           im6    = vc[j+6].imag();
                           xim[j+6] = im6;
                           re7    = vc[j+7].real();
                           xre[j+7] = re7;
                           im7    = vc[j+7].imag();
                           xim[j+7] = im7;
                           re0    = vc[j+8].real();
                           xre[j+8] = re0;
                           im0    = vc[j+8].imag();
                           xim[j+8] = im0;
                           re1    = vc[j+9].real();
                           xre[j+9] = re1;
                           im1    = vc[j+9].imag();
                           xim[j+9] = im1;
                           re2    = vc[j+10].real();
                           xre[j+10] = re2;
                           im2    = vc[j+10].imag();
                           xim[j+10] = im2;
                           re3    = vc[j+11].real();
                           xre[j+11] = re3;
                           im3    = vc[j+11].imag();
                           xim[j+11] = im3;
                           re4    = vc[j+12].real();
                           xre[j+12] = re4;
                           im4    = vc[j+12].imag();
                           xim[j+12] = im4;
                           re5    = vc[j+13].real();
                           xre[j+13] = re5;
                           im5    = vc[j+13].imag();
                           xim[j+13] = im5;
                           re6    = vc[j+14].real();
                           xre[j+14] = re6;
                           im6    = vc[j+14].imag();
                           xim[j+14] = im6;
                           re7    = vc[j+15].real();
                           xre[j+15] = re7;
                           im7    = vc[j+15].imag();
                           xim[j+15] = im7;
                      }
                 }


               /*    __ATTR_ALWAYS_INLINE__
	        
	           
                   static inline
                   void copy_c8_2xr8_unroll16x( double * __restrict xre,
                                                double * __restrict xim,
                                               std::complex<double> * __restrict vc,
                                               const int32_t n) {

                        if(__builtin_expect(0<=n,0)) { return;}
                        register double re0,im0,re1,im1,re2,im2;
                        register double re3,im3,re4,im4,re5,im5;
                        register double re6,im6,re7,im7;
                        int32_t j,m,m1;
                        m = n%16;
#if defined(__INTEL_COMPILER)
                        __assume_aligned(xre,64);
                        __assume_aligned(xim,64);
                        __assume_aligned(c,64);
#elif defined(__GNUC__) && (!defined __INTEL_COMPILER)
                        xre = (double*)__builtin_assume_aligned(xre,64);
                        xim = (double*)__builtin_assume_aligned(xim,64);
                        vc   = (std::complex<double>*)__builtin_assume_aligned(vc,64);
#endif
                        if(n != 0) {
#pragma omp simd aligned(xre:64) aligned(xim:64) aligned(c:64)
                           for(j = 0; j != m; ++j) {
                               const double re = vc[j].real();
                               const double im = vc[j].imag();
                               xre[j]         = re;
                               xim[j]         = im;
                           }
                            
                            if(n<16) {return;}
                       }

                       m1 = m+1;
#if defined(__INTEL_COMPILER)
#pragma code_align(32)
#endif
#pragma omp simd aligned(xre:64) aligned(xim:64) aligned(vc:64)
                       for(j = m1; j != n; j += 16) {
                           re0    = vc[j+0].real();
                           xre[j+0] = re0;
                           im0    = vc[j+0].imag();
                           xim[j+0] = im0;
                           re1    = vc[j+1].real();
                           xre[j+1] = re1;
                           im1    = vc[j+1].imag();
                           xim[j+1] = im1;
                           re2    = vc[j+2].real();
                           xre[j+2] = re2;
                           im2    = vc[j+2].imag();
                           xim[j+2] = im2;
                           re3    = vc[j+3].real();
                           xre[j+3] = re3;
                           im3    = vc[j+3].imag();
                           xim[j+3] = im3;
                           re4    = vc[j+4].real();
                           xre[j+4] = re4;
                           im4    = vc[j+4].imag();
                           xim[j+4] = im4;
                           re5    = vc[j+5].real();
                           xre[j+5] = re5;
                           im5    = vc[j+5].imag();
                           xim[j+5] = im5;
                           re6    = vc[j+6].real();
                           xre[j+6] = re6;
                           im6    = vc[j+6].imag();
                           xim[j+6] = im6;
                           re7    = vc[j+7].real();
                           xre[j+7] = re7;
                           im7    = vc[j+7].imag();
                           xim[j+7] = im7;
                           re0    = vc[j+8].real();
                           xre[j+8] = re0;
                           im0    = vc[j+8].imag();
                           xim[j+8] = im0;
                           re1    = vc[j+9].real();
                           xre[j+9] = re1;
                           im1    = vc[j+9].imag();
                           xim[j+9] = im1;
                           re2    = vc[j+10].real();
                           xre[j+10] = re2;
                           im2    = vc[j+10].imag();
                           xim[j+10] = im2;
                           re3    = vc[j+11].real();
                           xre[j+11] = re3;
                           im3    = vc[j+11].imag();
                           xim[j+11] = im3;
                           re4    = vc[j+12].real();
                           xre[j+12] = re4;
                           im4    = vc[j+12].imag();
                           xim[j+12] = im4;
                           re5    = vc[j+13].real();
                           xre[j+13] = re5;
                           im5    = vc[j+13].imag();
                           xim[j+13] = im5;
                           re6    = vc[j+14].real();
                           xre[j+14] = re6;
                           im6    = vc[j+14].imag();
                           xim[j+14] = im6;
                           re7    = vc[j+15].real();
                           xre[j+15] = re7;
                           im7    = vc[j+15].imag();
                           xim[j+15] = im7;
                      }
                 }
               */






      } // math


} // gms















#endif /*__GMS_COMPLEX_XMM2R8_HPP__*/
