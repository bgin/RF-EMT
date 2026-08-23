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

#ifndef __GMS_STAT_ANTENNA_THEORY_ZMM16R4_H__
#define __GMS_STAT_ANTENNA_THEORY_ZMM16R4_H__ 060620231414

namespace file_version {

    const unsigned int GMS_STAT_ANTENNA_THEORY_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_STAT_ANTENNA_THEORY_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_STAT_ANTENNA_THEORY_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_STAT_ANTENNA_THEORY_ZMM16R4_FULLVER =
      1000U*GMS_STAT_ANTENNA_THEORY_ZMM16R4_MAJOR+
      100U*GMS_STAT_ANTENNA_THEORY_ZMM16R4_MINOR+
      10U*GMS_STAT_ANTENNA_THEORY_ZMM16R4_MICRO;
    const char * const GMS_STAT_ANTENNA_THEORY_ZMM16R4_CREATION_DATE = "04-06-2023 14:14 PM +00200 (SUN 04 JUN 2023 GMT+2)";
    const char * const GMS_STAT_ANTENNA_THEORY_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_STAT_ANTENNA_THEORY_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_STAT_ANTENNA_THEORY_ZMM16R4_DESCRIPTION   = "AVX512 (single) optimized statistical antenna theory kernels.";

}


#include <immintrin.h>
#include <complex>
#include <cstdint>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_complex_zmm16r4.hpp"
#include "GMS_simd_utils.hpp"
#include "GMS_cspint_quad.hpp"
#include "GMS_avint_quad.hpp"
#include "GMS_cephes.h"

namespace  gms {


      namespace radiolocation {
      
      
                    /* The whole content is based on:
                       Шифрин Я.С. - "Вопросы статистической теории антенн"
                    */
                    
                    
               /*
                   Work (input) arrays for integrator 'cspint'
                   Multi-threaded i.e. (two-threaded)
               */
               __ATTR_ALIGN__(64) struct CSPINT_DATA_2T {
               
                       float * __restrict  Ya1; 
                       float * __restrict  Ya2; 
                       float * __restrict  Ya3; 
                       float * __restrict  Ea;  
                       float * __restrict  WRKa; 
                       float * __restrict  Yb1;
                       float * __restrict  Yb2; 
                       float * __restrict  Yb3; 
                       float * __restrict  Eb; 
                       float * __restrict  WRKb;  
               };
               
                /*
                   Work (input) arrays for integrator 'cspint'
                   Single-threaded.
               */
                __ATTR_ALIGN__(64) struct CSPINT_DATA_1T {
               
                       float * __restrict  Y1; 
                       float * __restrict  Y2; 
                       float * __restrict  Y3; 
                       float * __restrict  E;  
                       float * __restrict  WRK; 
                       
               };
               
               
                   /*
                        Formula 1.1, p. 14
                        Integrator cspint single data vector ZMM.
                   */
                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   std::complex<float> 
                   fth_f11_zmm16r4_cspint_16e_a(const float * __restrict __ATTR_ALIGN__(64) pAz,
                                            const float * __restrict __ATTR_ALIGN__(64) pphiz,
                                            const float * __restrict __ATTR_ALIGN__(64) pz,
                                            const float tht,
                                            const float gam,
                                            const int32_t L2); 
                 
                 
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   std::complex<float> 
                   fth_f11_zmm16r4_cspint_16e_u(const float * __restrict  pAz,
                                            const float * __restrict  pphiz,
                                            const float * __restrict  pz,
                                            const float tht,
                                            const float gam,
                                            const int32_t L2); 
                 
                 
                   /*
                        Formula 1.1, p. 14
                        Integrator avint (irregular abscissas) single data vector ZMM.
                   */
                 
                 
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   std::complex<float> 
                   fth_f11_zmm16r4_avint_16e_a(const float * __restrict __ATTR_ALIGN__(64) pAz,
                                               const float * __restrict __ATTR_ALIGN__(64) pphiz,
                                               const float * __restrict __ATTR_ALIGN__(64) pz,
                                               const float tht,
                                               const float gam,
                                               const int32_t L2,
                                               int32_t & ierr,
                                               int32_t & ieri); 
                                               
                 
                 
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   std::complex<float> 
                   fth_f11_zmm16r4_avint_16e_u(const float * __restrict  pAz,
                                               const float * __restrict  pphiz,
                                               const float * __restrict  pz,
                                               const float tht,
                                               const float gam,
                                               const int32_t L2,
                                               int32_t & ierr,
                                               int32_t & ieri); 
                 
                 
                 
                 
                   /*
                        Formula 1.1, p. 14
                        Integrator cspint n-elements data arrays.
                        Single-threaded.
                   */
                   
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   std::complex<float>
                   fth_f11_zmm16r4_cspint_ne_a( const float * __restrict __ATTR_ALIGN__(64) pAz,
                                            const float * __restrict __ATTR_ALIGN__(64) pphiz,
                                            const float * __restrict __ATTR_ALIGN__(64) pz,
                                            float * __restrict __ATTR_ALIGN__(64) intr,
                                            float * __restrict __ATTR_ALIGN__(64) inti,
                                            CSPINT_DATA_1T & csd,
                                            const float tht,
                                            const float gam,
                                            const int32_t L2,
                                            const int32_t NTAB); 
                                            
                
                
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   std::complex<float>
                   fth_f11_zmm16r4_cspint_ne_u( const float * __restrict  pAz,
                                            const float * __restrict  pphiz,
                                            const float * __restrict  pz,
                                            float * __restrict  intr,
                                            float * __restrict  inti,
                                            CSPINT_DATA_1T & csd,
                                            const float tht,
                                            const float gam,
                                            const int32_t L2,
                                            const int32_t NTAB); 
                 /*
                        Formula 1.1, p. 14
                        Integrator avint (irregular abscissas) n-elements data arrays.
                        Single-threaded.
                   */
                   
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           
                   std::complex<float>
                   fth_f11_zmm16r4_avint_ne_a( const float * __restrict __ATTR_ALIGN__(64) pAz,
                                            const float * __restrict __ATTR_ALIGN__(64) pphiz,
                                            const float * __restrict __ATTR_ALIGN__(64) pz,
                                            float * __restrict __ATTR_ALIGN__(64) intr,
                                            float * __restrict __ATTR_ALIGN__(64) inti,
                                            const float tht,
                                            const float gam,
                                            const int32_t L2,
                                            const int32_t n,
                                            int32_t & ierr,
                                            int32_t & ieri); 
                                            
                                            
                
                
                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   std::complex<float>
                   fth_f11_zmm16r4_avint_ne_u( const float * __restrict pAz,
                                            const float * __restrict  pphiz,
                                            const float * __restrict  pz,
                                            float * __restrict  intr,
                                            float * __restrict  inti,
                                            const float tht,
                                            const float gam,
                                            const int32_t L2,
                                            const int32_t n,
                                            int32_t & ierr,
                                            int32_t & ieri); 
                
                
                
                
                /*
                    Formula 1.3, p. 15
                    The level of an amplitude.
                */
                
                   __ATTR_ALWAYS_INLINE__
	           
	           static inline
                   __m512 Bx_f13_zmm16r4(const __m512 A0,
                                      const __m512 A) {
                                      
                          register __m512 Bx,rat;
                          rat = _mm512_div_ps(A,A0);
                          Bx  = xlogf(rat);
                          return (Bx);                    
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	          
	           static inline
                   __m512 Bx_f13_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pA0,
                                        const float * __restrict __ATTR_ALIGN__(64) pA) {
                              
                          register __m512 A0 = _mm512_load_ps(&pA0[0]);
                          register __m512 A  = _mm512_load_ps(&pA[0]);        
                          register __m512 Bx,rat;
                          rat = _mm512_div_ps(A,A0);
                          Bx  = xlogf(rat);
                          return (Bx);                    
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	          
	           static inline
                   __m512 Bx_f13_zmm16r4_u(const float * __restrict  pA0,
                                        const float * __restrict  pA) {
                              
                          register __m512 A0 = _mm512_loadu_ps(&pA0[0]);
                          register __m512 A  = _mm512_loadu_ps(&pA[0]);        
                          register __m512 Bx,rat;
                          rat = _mm512_div_ps(A,A0);
                          Bx  = xlogf(rat);
                          return (Bx);                    
                 }
                 
                 
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   void Bx_f13_zmm16r4_u10x_a(const float * __restrict __ATTR_ALIGN__(64) pA0,
                                              const float * __restrict __ATTR_ALIGN__(64) pA,
                                              float * __restrict __ATTR_ALIGN__(64) pB,
                                              const int32_t n); 
                
                
                
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   void Bx_f13_zmm16r4_u10x_u(const float * __restrict pA0,
                                              const float * __restrict  pA,
                                              float * __restrict  pB,
                                              const int32_t n); 
                
                
                
                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   void Bx_f13_zmm16r4_u2x_u(const float * __restrict pA0,
                                              const float * __restrict  pA,
                                              float * __restrict  pB,
                                              const int32_t n); 
                                              
                 
                 
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   void Bx_f13_zmm16r4_u2x_a(const float * __restrict __ATTR_ALIGN__(64) pA0,
                                              const float * __restrict __ATTR_ALIGN__(64) pA,
                                              float * __restrict __ATTR_ALIGN__(64) pB,
                                              const int32_t n); 
                 
                 
                 /*
                      Helper formula serving the needs of 
                      formula 1.4, p. 15
                      The mathematical expectation of random 
                      variable B(x).
                      Integrator 'cspint'.
                      Short data vector (ZMM-size).
                 */
                 
                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   float Ex_Bx_zmm16r4_cspint_16e_a(const float * __restrict __ATTR_ALIGN__(64) pBx,
                                                    const float * __restrict __ATTR_ALIGN__(64) ppdf,
                                                    const float * __restrict __ATTR_ALIGN__(64) px,
                                                    const float a,
                                                    const float b); 
                
                
                
                
                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   float Ex_Bx_zmm16r4_cspint_16e_u(const float * __restrict  pBx,
                                                    const float * __restrict  ppdf,
                                                    const float * __restrict  px,
                                                    const float a,
                                                    const float b); 
                
                
                  /*
                      Helper formula serving the needs of 
                      formula 1.4, p. 15
                      The mathematical expectation of random 
                      variable B(x).
                      Integrator 'avint' irregular abscissas.
                      Short data vector (ZMM-size).
                 */
                 
                 
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   float Ex_Bx_zmm16r4_avint_16e_a(const float * __restrict __ATTR_ALIGN__(64) pBx,
                                                    const float * __restrict __ATTR_ALIGN__(64) ppdf,
                                                    const float * __restrict __ATTR_ALIGN__(64) px,
                                                    const float a,
                                                    const float b,
                                                    int32_t & ier);
                
                
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   float Ex_Bx_zmm16r4_avint_16e_u(const float * __restrict  pBx,
                                                    const float * __restrict  ppdf,
                                                    const float * __restrict  px,
                                                    const float a,
                                                    const float b,
                                                    int32_t & ier);
                
                 /*
                      Helper formula serving the needs of 
                      formula 1.4, p. 15
                      The mathematical expectation of random 
                      variable B(x).
                      Integrator 'cspint'.
                      Long data vector (arrays).
                 */
                 
                 
                 
                 
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   float Ex_Bx_zmm16r4_cspint_ne_a(const float * __restrict __ATTR_ALIGN__(64) pBx,
                                                   const float * __restrict __ATTR_ALIGN__(64) ppdf,
                                                   const float * __restrict __ATTR_ALIGN__(64) px,
                                                   float * __restrict __ATTR_ALIGN__(64)       intr,
                                                   CSPINT_DATA_1T csd,
                                                   const float a,
                                                   const float b,
                                                   const int32_t NTAB); 
                                                   
                  
                  
                  
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   float Ex_Bx_zmm16r4_cspint_ne_u(const float * __restrict  pBx,
                                                   const float * __restrict  ppdf,
                                                   const float * __restrict  px,
                                                   float * __restrict intr,
                                                   CSPINT_DATA_1T csd,
                                                   const float a,
                                                   const float b,
                                                   const int32_t NTAB); 
                  
                  
                   /*
                      Helper formula serving the needs of 
                      formula 1.4, p. 15
                      The mathematical expectation of random 
                      variable B(x).
                      Integrator 'avint' irregular abscissas.
                      Long data vector (arrays).
                 */
                 
                 
                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           
                   float Ex_Bx_zmm16r4_avint_ne_a(const float * __restrict __ATTR_ALIGN__(64) pBx,
                                                   const float * __restrict __ATTR_ALIGN__(64) ppdf,
                                                   const float * __restrict __ATTR_ALIGN__(64) px,
                                                   float * __restrict __ATTR_ALIGN__(64)       intr,
                                                   const int32_t n,
                                                   const float a,
                                                   const float b,
                                                   int32_t & ier); 
                                                   
                  
                  
                  
                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   float Ex_Bx_zmm16r4_avint_ne_u(const float * __restrict  pBx,
                                                   const float * __restrict  ppdf,
                                                   const float * __restrict  px,
                                                   float * __restrict intr,
                                                   const int32_t n,
                                                   const float a,
                                                   const float b,
                                                   int32_t & ier); 
                  
                  
                    
                 /*
                      Helper formula serving the needs of 
                      formula 1.4, p. 15
                      The mathematical expectation of phase random 
                      variable phi(x).
                      Integrator 'cspint'.
                      Short data vector (ZMM-size).
                 */
                 
                 
                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   float Ex_phix_zmm16r4_cspint_16e_a(const float * __restrict __ATTR_ALIGN__(64) pphix,
                                                      const float * __restrict __ATTR_ALIGN__(64) ppdf,
                                                      const float * __restrict __ATTR_ALIGN__(64) px,
                                                      const float a,
                                                      const float b); 
                                                      
                
                
                
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           
                   float Ex_phix_zmm16r4_cspint_16e_u(const float * __restrict  pphix,
                                                      const float * __restrict  ppdf,
                                                      const float * __restrict  px,
                                                      const float a,
                                                      const float b); 
                
                
                  /*
                      Helper formula serving the needs of 
                      formula 1.4, p. 15
                      The mathematical expectation of phase random 
                      variable phi(x).
                      Integrator 'avint' irregular abscissas.
                      Short data vector (ZMM-size).
                 */
                 
                 
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   float Ex_phix_zmm16r4_avint_16e_a(const float * __restrict __ATTR_ALIGN__(64) pphix,
                                                      const float * __restrict __ATTR_ALIGN__(64) ppdf,
                                                      const float * __restrict __ATTR_ALIGN__(64) px,
                                                      const float a,
                                                      const float b,
                                                      int32_t & ier); 
                
                
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   float Ex_phix_zmm16r4_avint_16e_u(const float * __restrict pphix,
                                                      const float * __restrict ppdf,
                                                      const float * __restrict  px,
                                                      const float a,
                                                      const float b,
                                                      int32_t & ier); 
                
                
                
                  /*
                      Helper formula serving the needs of 
                      formula 1.4, p. 15
                      The mathematical expectation of phase random 
                      variable phix(x).
                      Integrator 'cspint'.
                      Long data vector (arrays).
                 */
                 
                 
                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   float Ex_phix_zmm16r4_cspint_ne_a(const float * __restrict __ATTR_ALIGN__(64) pphix,
                                                   const float * __restrict __ATTR_ALIGN__(64) ppdf,
                                                   const float * __restrict __ATTR_ALIGN__(64) px,
                                                   float * __restrict __ATTR_ALIGN__(64)       intr,
                                                   CSPINT_DATA_1T csd,
                                                   const float a,
                                                   const float b,
                                                   const int32_t NTAB); 
                 
                 
                 
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   float Ex_phix_zmm16r4_cspint_ne_u(const float * __restrict pphix,
                                                   const float * __restrict  ppdf,
                                                   const float * __restrict  px,
                                                   float * __restrict intr,
                                                   CSPINT_DATA_1T csd,
                                                   const float a,
                                                   const float b,
                                                   const int32_t NTAB); 
                  
                  
                   /*
                      Helper formula serving the needs of 
                      formula 1.4, p. 15
                      The mathematical expectation of phase random 
                      variable phix(x).
                      Integrator 'avint' irregular abscissas.
                      Long data vector (arrays).
                 */
                 
                 
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   float Ex_phix_zmm16r4_avint_ne_a(const float * __restrict __ATTR_ALIGN__(64) pphix,
                                                   const float * __restrict __ATTR_ALIGN__(64) ppdf,
                                                   const float * __restrict __ATTR_ALIGN__(64) px,
                                                   float * __restrict __ATTR_ALIGN__(64)       intr,
                                                   const float a,
                                                   const float b,
                                                   const int32_t n,
                                                   int32_t ier); 
                  
                  
                  
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	        
                   float Ex_phix_zmm16r4_avint_ne_u(const float * __restrict  pphix,
                                                   const float * __restrict  ppdf,
                                                   const float * __restrict px,
                                                   float * __restrict  intr,
                                                   const float a,
                                                   const float b,
                                                   const int32_t n,
                                                   int32_t ier); 
                  
                  
                  /*
                      Formula 1.5, p. 16
                      Far-zone Electric-field
                  */
                  
                   __ATTR_ALWAYS_INLINE__
	          
	           static inline
                   void Efz_f15_zmm16r4(const __m512 fthr,
                                        const __m512 fthi,
                                        const __m512 Ethr,
                                        const __m512 Ethi,
                                        const __m512 Ephr,
                                        const __m512 Ephi,
                                        __m512 * __restrict Efztr,
                                        __m512 * __restrict Efzti,
                                        __m512 * __restrict Efzpr,
                                        __m512 * __restrict Efzpi) {
                                        
                        
                        register __m512 t0r,t0i,t1r,ti1;
                        cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                        *Efztr = t0r;
                        *Efzti = t0i;
                        cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                        *Efzpr = t1r;
                        *Efzpi = t1i;         
                 }
                 
                 
                 
                   __ATTR_ALWAYS_INLINE__
	          
	           static inline
                   void Efz_f15_zmm16r4_a(const  float  __restrict __ATTR_ALIGN__(64) pfthr,
                                          const  float  __restrict __ATTR_ALIGN__(64) pfthi,
                                          const  float  __restrict __ATTR_ALIGN__(64) pEthr,
                                          const  float  __restrict __ATTR_ALIGN__(64) pEthi,
                                          const  float  __restrict __ATTR_ALIGN__(64) pEphr,
                                          const  float  __restrict __ATTR_ALIGN__(64) pEphi,
                                          float  __restrict __ATTR_ALIGN__(64) Efztr,
                                          float  __restrict __ATTR_ALIGN__(64) Efzti,
                                          float  __restrict __ATTR_ALIGN__(64) Efzpr,
                                          float  __restrict __ATTR_ALIGN__(64) Efzpi) {
                                        
                        register __m512 fthr = _mm512_load_ps(&pfthr[0]);
                        register __m512 fthi = _mm512_load_ps(&pfthi[0]);
                        register __m512 Ethr = _mm512_load_ps(&pEthr[0]);
                        register __m512 Ethi = _mm512_load_ps(&pEthi[0]);
                        register __m512 Ephr = _mm512_load_ps(&pEphr[0]);
                        register __m512 Ephi = _mm512_load_ps(&pEphi[0]);
                        register __m512 t0r,t0i,t1r,ti1;
                        cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                        _mm512_store_ps(&Efztr[0] ,t0r);
                        _mm512_store_ps(&Efzti[0] ,t0i);
                        cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                        _mm512_store_ps(&Efzpr[0] ,t1r);
                        _mm512_store_ps(&Efzpi[0] ,t1i);         
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	          
	           static inline
                   void Efz_f15_zmm16r4_u(const  float  __restrict  pfthr,
                                          const  float  __restrict  pfthi,
                                          const  float  __restrict  pEthr,
                                          const  float  __restrict  pEthi,
                                          const  float  __restrict  pEphr,
                                          const  float  __restrict  pEphi,
                                          float  __restrict  Efztr,
                                          float  __restrict  Efzti,
                                          float  __restrict  Efzpr,
                                          float  __restrict  Efzpi) {
                                        
                        register __m512 fthr = _mm512_loadu_ps(&pfthr[0]);
                        register __m512 fthi = _mm512_loadu_ps(&pfthi[0]);
                        register __m512 Ethr = _mm512_loadu_ps(&pEthr[0]);
                        register __m512 Ethi = _mm512_loadu_ps(&pEthi[0]);
                        register __m512 Ephr = _mm512_loadu_ps(&pEphr[0]);
                        register __m512 Ephi = _mm512_loadu_ps(&pEphi[0]);
                        register __m512 t0r,t0i,t1r,ti1;
                        cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                        _mm512_storeu_ps(&Efztr[0] ,t0r);
                        _mm512_storeu_ps(&Efzti[0] ,t0i);
                        cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                        _mm512_storeu_ps(&Efzpr[0] ,t1r);
                        _mm512_storeu_ps(&Efzpi[0] ,t1i);         
                 }
                 
                 
                  /*
                      Formula 1.5, p. 16
                      Far-zone Electric-field, n-terms
                  */
                  
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   void Efz_f15_zmm16r4_unroll10x(const __m512 * __restrict __ATTR_ALIGN__(64) pfthr,
                                                  const __m512 * __restrict __ATTR_ALIGN__(64) pfthi,
                                                  const __m512 * __restrict __ATTR_ALIGN__(64) pEthr,
                                                  const __m512 * __restrict __ATTR_ALIGN__(64) pEthi,
                                                  const __m512 * __restrict __ATTR_ALIGN__(64) pEphr,
                                                  const __m512 * __restrict __ATTR_ALIGN__(64) pEphi,
                                                  __m512 * __restrict __ATTR_ALIGN__(64) Efztr,
                                                  __m512 * __restrict __ATTR_ALIGN__(64) Efzti,
                                                  __m512 * __restrict __ATTR_ALIGN__(64) Efzpr,
                                                  __m512 * __restrict __ATTR_ALIGN__(64) Efzpi,
                                                  const int32_t n);
                 
                 
                 
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   void Efz_f15_zmm16r4_unroll8x(const __m512 * __restrict __ATTR_ALIGN__(64) pfthr,
                                                  const __m512 * __restrict __ATTR_ALIGN__(64) pfthi,
                                                  const __m512 * __restrict __ATTR_ALIGN__(64) pEthr,
                                                  const __m512 * __restrict __ATTR_ALIGN__(64) pEthi,
                                                  const __m512 * __restrict __ATTR_ALIGN__(64) pEphr,
                                                  const __m512 * __restrict __ATTR_ALIGN__(64) pEphi,
                                                  __m512 * __restrict __ATTR_ALIGN__(64) Efztr,
                                                  __m512 * __restrict __ATTR_ALIGN__(64) Efzti,
                                                  __m512 * __restrict __ATTR_ALIGN__(64) Efzpr,
                                                  __m512 * __restrict __ATTR_ALIGN__(64) Efzpi,
                                                  const int32_t n); 
                 
                 
                 
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   void Efz_f15_zmm16r4_unroll4x(const __m512 * __restrict __ATTR_ALIGN__(64) pfthr,
                                                  const __m512 * __restrict __ATTR_ALIGN__(64) pfthi,
                                                  const __m512 * __restrict __ATTR_ALIGN__(64) pEthr,
                                                  const __m512 * __restrict __ATTR_ALIGN__(64) pEthi,
                                                  const __m512 * __restrict __ATTR_ALIGN__(64) pEphr,
                                                  const __m512 * __restrict __ATTR_ALIGN__(64) pEphi,
                                                  __m512 * __restrict __ATTR_ALIGN__(64) Efztr,
                                                  __m512 * __restrict __ATTR_ALIGN__(64) Efzti,
                                                  __m512 * __restrict __ATTR_ALIGN__(64) Efzpr,
                                                  __m512 * __restrict __ATTR_ALIGN__(64) Efzpi,
                                                  const int32_t n); 
                
                
              
                
                
                 
                 
      
      } // radiolocation


} // radiolocation


















#endif /*__GMS_STAT_ANTENNA_THEORY_ZMM16R4_H__*/
