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

#ifndef __GMS_STAT_ANTENNA_THEORY_ZMM8R8_H__
#define __GMS_STAT_ANTENNA_THEORY_ZMM8R8_H__ 040420251840

namespace file_version {

    const unsigned int GMS_STAT_ANTENNA_THEORY_ZMM8R8_MAJOR = 1U;
    const unsigned int GMS_STAT_ANTENNA_THEORY_ZMM8R8_MINOR = 0U;
    const unsigned int GMS_STAT_ANTENNA_THEORY_ZMM8R8_MICRO = 0U;
    const unsigned int GMS_STAT_ANTENNA_THEORY_ZMM8R8_FULLVER =
      1000U*GMS_STAT_ANTENNA_THEORY_ZMM8R8_MAJOR+
      100U*GMS_STAT_ANTENNA_THEORY_ZMM8R8_MINOR+
      10U*GMS_STAT_ANTENNA_THEORY_ZMM8R8_MICRO;
    const char * const GMS_STAT_ANTENNA_THEORY_ZMM8R8_CREATION_DATE = "04-04-2025 18:40 PM +00200 (FRI 04 APR 2025 GMT+2)";
    const char * const GMS_STAT_ANTENNA_THEORY_ZMM8R8_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_STAT_ANTENNA_THEORY_ZMM8R8_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_STAT_ANTENNA_THEORY_ZMM8R8_DESCRIPTION   = "AVX512 (double) optimized statistical antenna theory kernels.";

}


#include <immintrin.h>
#include <complex>
#include <cstdint>
#include "GMS_config.h"
#include "GMS_sleefsimddp.hpp"
#include "GMS_complex_zmm8r8.hpp"
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
               
                       double * __restrict  Ya1; 
                       double * __restrict  Ya2; 
                       double * __restrict  Ya3; 
                       double * __restrict  Ea;  
                       double * __restrict  WRKa; 
                       double * __restrict  Yb1;
                       double * __restrict  Yb2; 
                       double * __restrict  Yb3; 
                       double * __restrict  Eb; 
                       double * __restrict  WRKb;  
               };
               
                /*
                   Work (input) arrays for integrator 'cspint'
                   Single-threaded.
               */
                __ATTR_ALIGN__(64) struct CSPINT_DATA_1T {
               
                       double * __restrict  Y1; 
                       double * __restrict  Y2; 
                       double * __restrict  Y3; 
                       double * __restrict  E;  
                       double * __restrict  WRK; 
                       
               };
               
               
                   /*
                        Formula 1.1, p. 14
                        Integrator cspint single data vector ZMM.
                   */
                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   std::complex<double> 
                   fth_f11_zmm8r8_cspint_8e_a(const double * __restrict __ATTR_ALIGN__(64) pAz,
                                            const double * __restrict __ATTR_ALIGN__(64) pphiz,
                                            const double * __restrict __ATTR_ALIGN__(64) pz,
                                            const double tht,
                                            const double gam,
                                            const int32_t L2); 
                 
                 
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   std::complex<double> 
                   fth_f11_zmm8r8_cspint_8e_u(const double * __restrict  pAz,
                                            const double * __restrict  pphiz,
                                            const double * __restrict  pz,
                                            const double tht,
                                            const double gam,
                                            const int32_t L2); 
                 
                 
                   /*
                        Formula 1.1, p. 14
                        Integrator avint (irregular abscissas) single data vector ZMM.
                   */
                 
                 
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   std::complex<double> 
                   fth_f11_zmm8r8_avint_8e_a(const double * __restrict __ATTR_ALIGN__(64) pAz,
                                               const double * __restrict __ATTR_ALIGN__(64) pphiz,
                                               const double * __restrict __ATTR_ALIGN__(64) pz,
                                               const double tht,
                                               const double gam,
                                               const int32_t L2,
                                               int32_t & ierr,
                                               int32_t & ieri); 
                                               
                 
                 
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   std::complex<double> 
                   fth_f11_zmm8r8_avint_8e_u(const double * __restrict  pAz,
                                               const double * __restrict  pphiz,
                                               const double * __restrict  pz,
                                               const double tht,
                                               const double gam,
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
	         
                   std::complex<double>
                   fth_f11_zmm8r8_cspint_ne_a( const double * __restrict __ATTR_ALIGN__(64) pAz,
                                            const double * __restrict __ATTR_ALIGN__(64) pphiz,
                                            const double * __restrict __ATTR_ALIGN__(64) pz,
                                            double * __restrict __ATTR_ALIGN__(64) intr,
                                            double * __restrict __ATTR_ALIGN__(64) inti,
                                            CSPINT_DATA_1T & csd,
                                            const double tht,
                                            const double gam,
                                            const int32_t L2,
                                            const int32_t NTAB); 
                                            
                
                
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   std::complex<double>
                   fth_f11_zmm8r8_cspint_ne_u( const double * __restrict  pAz,
                                            const double * __restrict  pphiz,
                                            const double * __restrict  pz,
                                            double * __restrict  intr,
                                            double * __restrict  inti,
                                            CSPINT_DATA_1T & csd,
                                            const double tht,
                                            const double gam,
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
	           
                   std::complex<double>
                   fth_f11_zmm8r8_avint_ne_a( const double * __restrict __ATTR_ALIGN__(64) pAz,
                                            const double * __restrict __ATTR_ALIGN__(64) pphiz,
                                            const double * __restrict __ATTR_ALIGN__(64) pz,
                                            double * __restrict __ATTR_ALIGN__(64) intr,
                                            double * __restrict __ATTR_ALIGN__(64) inti,
                                            const double tht,
                                            const double gam,
                                            const int32_t L2,
                                            const int32_t n,
                                            int32_t & ierr,
                                            int32_t & ieri); 
                                            
                                            
                
                
                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   std::complex<double>
                   fth_f11_zmm8r8_avint_ne_u( const double * __restrict pAz,
                                            const double * __restrict  pphiz,
                                            const double * __restrict  pz,
                                            double * __restrict  intr,
                                            double * __restrict  inti,
                                            const double tht,
                                            const double gam,
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
                   __m512d Bx_f13_zmm8r8(const __m512d A0,
                                      const __m512d A) {
                                      
                           __m512d Bx,rat;
                          rat = _mm512_div_pd(A,A0);
                          Bx  = xlogf(rat);
                          return (Bx);                    
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	          
	           static inline
                   __m512d Bx_f13_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pA0,
                                        const double * __restrict __ATTR_ALIGN__(64) pA) {
                              
                           __m512d A0 = _mm512_load_pd(&pA0[0]);
                           __m512d A  = _mm512_load_pd(&pA[0]);        
                           __m512d Bx,rat;
                          rat = _mm512_div_pd(A,A0);
                          Bx  = xlogf(rat);
                          return (Bx);                    
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	          
	           static inline
                   __m512d Bx_f13_zmm8r8_u(const double * __restrict  pA0,
                                        const double * __restrict  pA) {
                              
                           __m512d A0 = _mm512_loadu_pd(&pA0[0]);
                           __m512d A  = _mm512_loadu_pd(&pA[0]);        
                           __m512d Bx,rat;
                          rat = _mm512_div_pd(A,A0);
                          Bx  = xlogf(rat);
                          return (Bx);                    
                 }
                 
                 
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   void Bx_f13_zmm8r8_u10x_a(const double * __restrict __ATTR_ALIGN__(64) pA0,
                                              const double * __restrict __ATTR_ALIGN__(64) pA,
                                              double * __restrict __ATTR_ALIGN__(64) pB,
                                              const int32_t n); 
                
                
                
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   void Bx_f13_zmm8r8_u10x_u(const double * __restrict pA0,
                                              const double * __restrict  pA,
                                              double * __restrict  pB,
                                              const int32_t n); 
                
                
                
                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   void Bx_f13_zmm8r8_u2x_u(const double * __restrict pA0,
                                              const double * __restrict  pA,
                                              double * __restrict  pB,
                                              const int32_t n); 
                                              
                 
                 
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   void Bx_f13_zmm8r8_u2x_a(const double * __restrict __ATTR_ALIGN__(64) pA0,
                                              const double * __restrict __ATTR_ALIGN__(64) pA,
                                              double * __restrict __ATTR_ALIGN__(64) pB,
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
	         
                   double Ex_Bx_zmm8r8_cspint_8e_a(const double * __restrict __ATTR_ALIGN__(64) pBx,
                                                    const double * __restrict __ATTR_ALIGN__(64) ppdf,
                                                    const double * __restrict __ATTR_ALIGN__(64) px,
                                                    const double a,
                                                    const double b); 
                
                
                
                
                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   double Ex_Bx_zmm8r8_cspint_8e_u(const double * __restrict  pBx,
                                                    const double * __restrict  ppdf,
                                                    const double * __restrict  px,
                                                    const double a,
                                                    const double b); 
                
                
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
	          
                   double Ex_Bx_zmm8r8_avint_8e_a(const double * __restrict __ATTR_ALIGN__(64) pBx,
                                                    const double * __restrict __ATTR_ALIGN__(64) ppdf,
                                                    const double * __restrict __ATTR_ALIGN__(64) px,
                                                    const double a,
                                                    const double b,
                                                    int32_t & ier);
                
                
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   double Ex_Bx_zmm8r8_avint_8e_u(const double * __restrict  pBx,
                                                    const double * __restrict  ppdf,
                                                    const double * __restrict  px,
                                                    const double a,
                                                    const double b,
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
	          
                   double Ex_Bx_zmm8r8_cspint_ne_a(const double * __restrict __ATTR_ALIGN__(64) pBx,
                                                   const double * __restrict __ATTR_ALIGN__(64) ppdf,
                                                   const double * __restrict __ATTR_ALIGN__(64) px,
                                                   double * __restrict __ATTR_ALIGN__(64)       intr,
                                                   CSPINT_DATA_1T csd,
                                                   const double a,
                                                   const double b,
                                                   const int32_t NTAB); 
                                                   
                  
                  
                  
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   double Ex_Bx_zmm8r8_cspint_ne_u(const double * __restrict  pBx,
                                                   const double * __restrict  ppdf,
                                                   const double * __restrict  px,
                                                   double * __restrict intr,
                                                   CSPINT_DATA_1T csd,
                                                   const double a,
                                                   const double b,
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
	           
                   double Ex_Bx_zmm8r8_avint_ne_a(const double * __restrict __ATTR_ALIGN__(64) pBx,
                                                   const double * __restrict __ATTR_ALIGN__(64) ppdf,
                                                   const double * __restrict __ATTR_ALIGN__(64) px,
                                                   double * __restrict __ATTR_ALIGN__(64)       intr,
                                                   const int32_t n,
                                                   const double a,
                                                   const double b,
                                                   int32_t & ier); 
                                                   
                  
                  
                  
                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   double Ex_Bx_zmm8r8_avint_ne_u(const double * __restrict  pBx,
                                                   const double * __restrict  ppdf,
                                                   const double * __restrict  px,
                                                   double * __restrict intr,
                                                   const int32_t n,
                                                   const double a,
                                                   const double b,
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
	          
                   double Ex_phix_zmm8r8_cspint_8e_a(const double * __restrict __ATTR_ALIGN__(64) pphix,
                                                      const double * __restrict __ATTR_ALIGN__(64) ppdf,
                                                      const double * __restrict __ATTR_ALIGN__(64) px,
                                                      const double a,
                                                      const double b); 
                                                      
                
                
                
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           
                   double Ex_phix_zmm8r8_cspint_8e_u(const double * __restrict  pphix,
                                                      const double * __restrict  ppdf,
                                                      const double * __restrict  px,
                                                      const double a,
                                                      const double b); 
                
                
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
	         
                   double Ex_phix_zmm8r8_avint_8e_a(const double * __restrict __ATTR_ALIGN__(64) pphix,
                                                      const double * __restrict __ATTR_ALIGN__(64) ppdf,
                                                      const double * __restrict __ATTR_ALIGN__(64) px,
                                                      const double a,
                                                      const double b,
                                                      int32_t & ier); 
                
                
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   double Ex_phix_zmm8r8_avint_8e_u(const double * __restrict pphix,
                                                      const double * __restrict ppdf,
                                                      const double * __restrict  px,
                                                      const double a,
                                                      const double b,
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
	          
                   double Ex_phix_zmm8r8_cspint_ne_a(const double * __restrict __ATTR_ALIGN__(64) pphix,
                                                   const double * __restrict __ATTR_ALIGN__(64) ppdf,
                                                   const double * __restrict __ATTR_ALIGN__(64) px,
                                                   double * __restrict __ATTR_ALIGN__(64)       intr,
                                                   CSPINT_DATA_1T csd,
                                                   const double a,
                                                   const double b,
                                                   const int32_t NTAB); 
                 
                 
                 
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   double Ex_phix_zmm8r8_cspint_ne_u(const double * __restrict pphix,
                                                   const double * __restrict  ppdf,
                                                   const double * __restrict  px,
                                                   double * __restrict intr,
                                                   CSPINT_DATA_1T csd,
                                                   const double a,
                                                   const double b,
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
	          
                   double Ex_phix_zmm8r8_avint_ne_a(const double * __restrict __ATTR_ALIGN__(64) pphix,
                                                   const double * __restrict __ATTR_ALIGN__(64) ppdf,
                                                   const double * __restrict __ATTR_ALIGN__(64) px,
                                                   double * __restrict __ATTR_ALIGN__(64)       intr,
                                                   const double a,
                                                   const double b,
                                                   const int32_t n,
                                                   int32_t ier); 
                  
                  
                  
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	        
                   double Ex_phix_zmm8r8_avint_ne_u(const double * __restrict  pphix,
                                                   const double * __restrict  ppdf,
                                                   const double * __restrict px,
                                                   double * __restrict  intr,
                                                   const double a,
                                                   const double b,
                                                   const int32_t n,
                                                   int32_t ier); 
                  
                  
                  /*
                      Formula 1.5, p. 16
                      Far-zone Electric-field
                  */
                  
                   __ATTR_ALWAYS_INLINE__
	          
	           static inline
                   void Efz_f15_zmm8r8(const __m512d fthr,
                                        const __m512d fthi,
                                        const __m512d Ethr,
                                        const __m512d Ethi,
                                        const __m512d Ephr,
                                        const __m512d Ephi,
                                        __m512d * __restrict Efztr,
                                        __m512d * __restrict Efzti,
                                        __m512d * __restrict Efzpr,
                                        __m512d * __restrict Efzpi) {
                                        
                        
                         __m512d t0r,t0i,t1r,ti1;
                        cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                        *Efztr = t0r;
                        *Efzti = t0i;
                        cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                        *Efzpr = t1r;
                        *Efzpi = t1i;         
                 }
                 
                 
                 
                   __ATTR_ALWAYS_INLINE__
	          
	           static inline
                   void Efz_f15_zmm8r8_a(const  double  __restrict __ATTR_ALIGN__(64) pfthr,
                                          const  double  __restrict __ATTR_ALIGN__(64) pfthi,
                                          const  double  __restrict __ATTR_ALIGN__(64) pEthr,
                                          const  double  __restrict __ATTR_ALIGN__(64) pEthi,
                                          const  double  __restrict __ATTR_ALIGN__(64) pEphr,
                                          const  double  __restrict __ATTR_ALIGN__(64) pEphi,
                                          double  __restrict __ATTR_ALIGN__(64) Efztr,
                                          double  __restrict __ATTR_ALIGN__(64) Efzti,
                                          double  __restrict __ATTR_ALIGN__(64) Efzpr,
                                          double  __restrict __ATTR_ALIGN__(64) Efzpi) {
                                        
                         __m512d fthr = _mm512_load_pd(&pfthr[0]);
                         __m512d fthi = _mm512_load_pd(&pfthi[0]);
                         __m512d Ethr = _mm512_load_pd(&pEthr[0]);
                         __m512d Ethi = _mm512_load_pd(&pEthi[0]);
                         __m512d Ephr = _mm512_load_pd(&pEphr[0]);
                         __m512d Ephi = _mm512_load_pd(&pEphi[0]);
                         __m512d t0r,t0i,t1r,ti1;
                        cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                        _mm512_store_pd(&Efztr[0] ,t0r);
                        _mm512_store_pd(&Efzti[0] ,t0i);
                        cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                        _mm512_store_pd(&Efzpr[0] ,t1r);
                        _mm512_store_pd(&Efzpi[0] ,t1i);         
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	          
	           static inline
                   void Efz_f15_zmm8r8_u(const  double  __restrict  pfthr,
                                          const  double  __restrict  pfthi,
                                          const  double  __restrict  pEthr,
                                          const  double  __restrict  pEthi,
                                          const  double  __restrict  pEphr,
                                          const  double  __restrict  pEphi,
                                          double  __restrict  Efztr,
                                          double  __restrict  Efzti,
                                          double  __restrict  Efzpr,
                                          double  __restrict  Efzpi) {
                                        
                         __m512d fthr = _mm512_loadu_pd(&pfthr[0]);
                         __m512d fthi = _mm512_loadu_pd(&pfthi[0]);
                         __m512d Ethr = _mm512_loadu_pd(&pEthr[0]);
                         __m512d Ethi = _mm512_loadu_pd(&pEthi[0]);
                         __m512d Ephr = _mm512_loadu_pd(&pEphr[0]);
                         __m512d Ephi = _mm512_loadu_pd(&pEphi[0]);
                         __m512d t0r,t0i,t1r,ti1;
                        cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                        _mm512_storeu_pd(&Efztr[0] ,t0r);
                        _mm512_storeu_pd(&Efzti[0] ,t0i);
                        cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                        _mm512_storeu_pd(&Efzpr[0] ,t1r);
                        _mm512_storeu_pd(&Efzpi[0] ,t1i);         
                 }
                 
                 
                  /*
                      Formula 1.5, p. 16
                      Far-zone Electric-field, n-terms
                  */
                  
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   void Efz_f15_zmm8r8_unroll10x(const __m512d * __restrict __ATTR_ALIGN__(64) pfthr,
                                                  const __m512d * __restrict __ATTR_ALIGN__(64) pfthi,
                                                  const __m512d * __restrict __ATTR_ALIGN__(64) pEthr,
                                                  const __m512d * __restrict __ATTR_ALIGN__(64) pEthi,
                                                  const __m512d * __restrict __ATTR_ALIGN__(64) pEphr,
                                                  const __m512d * __restrict __ATTR_ALIGN__(64) pEphi,
                                                  __m512d * __restrict __ATTR_ALIGN__(64) Efztr,
                                                  __m512d * __restrict __ATTR_ALIGN__(64) Efzti,
                                                  __m512d * __restrict __ATTR_ALIGN__(64) Efzpr,
                                                  __m512d * __restrict __ATTR_ALIGN__(64) Efzpi,
                                                  const int32_t n);
                 
                 
                 
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   void Efz_f15_zmm8r8_unroll8x(const __m512d * __restrict __ATTR_ALIGN__(64) pfthr,
                                                  const __m512d * __restrict __ATTR_ALIGN__(64) pfthi,
                                                  const __m512d * __restrict __ATTR_ALIGN__(64) pEthr,
                                                  const __m512d * __restrict __ATTR_ALIGN__(64) pEthi,
                                                  const __m512d * __restrict __ATTR_ALIGN__(64) pEphr,
                                                  const __m512d * __restrict __ATTR_ALIGN__(64) pEphi,
                                                  __m512d * __restrict __ATTR_ALIGN__(64) Efztr,
                                                  __m512d * __restrict __ATTR_ALIGN__(64) Efzti,
                                                  __m512d * __restrict __ATTR_ALIGN__(64) Efzpr,
                                                  __m512d * __restrict __ATTR_ALIGN__(64) Efzpi,
                                                  const int32_t n); 
                 
                 
                 
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   void Efz_f15_zmm8r8_unroll4x(const __m512d * __restrict __ATTR_ALIGN__(64) pfthr,
                                                  const __m512d * __restrict __ATTR_ALIGN__(64) pfthi,
                                                  const __m512d * __restrict __ATTR_ALIGN__(64) pEthr,
                                                  const __m512d * __restrict __ATTR_ALIGN__(64) pEthi,
                                                  const __m512d * __restrict __ATTR_ALIGN__(64) pEphr,
                                                  const __m512d * __restrict __ATTR_ALIGN__(64) pEphi,
                                                  __m512d * __restrict __ATTR_ALIGN__(64) Efztr,
                                                  __m512d * __restrict __ATTR_ALIGN__(64) Efzti,
                                                  __m512d * __restrict __ATTR_ALIGN__(64) Efzpr,
                                                  __m512d * __restrict __ATTR_ALIGN__(64) Efzpi,
                                                  const int32_t n); 
                
                
              
                
                
                 
                 
      
      } // radiolocation


} // radiolocation


















#endif /*__GMS_STAT_ANTENNA_THEORY_ZMM8R8_H__*/
