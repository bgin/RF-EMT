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

#include <cmath>
#include "GMS_stat_antenna_theory_zmm8r8.h"


               
               
                   /*
                        Formula 1.1, p. 14
                        Integrator cspint single data vector ZMM.
                   */
                 
                   std::complex<double> 
                   fth_f11_zmm8r8_cspint_8e_a(const double * __restrict __ATTR_ALIGN__(64) pAz,
                                            const double * __restrict __ATTR_ALIGN__(64) pphiz,
                                            const double * __restrict __ATTR_ALIGN__(64) pz,
                                            const double tht,
                                            const double gam,
                                            const int32_t L2) {
                                            
                         __ATTR_ALIGN__(64) double intr[8] = {};
                         __ATTR_ALIGN__(64) double inti[8] = {}; 
                         __ATTR_ALIGN__(64) double Y1[8]   = {};
                         __ATTR_ALIGN__(64) double Y2[8]   = {};
                         __ATTR_ALIGN__(64) double Y3[8]   = {}; 
                         __ATTR_ALIGN__(64) double E[8]    = {};
                         __ATTR_ALIGN__(64) double WRK[8]  = {}; 
                         constexpr int32_t NTAB = 8;               
                         constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                          __m512d Az    = _mm512_load_pd(&pAz[0]);
                          __m512d phiz  = _mm512_load_pd(&pphiz[0]);
                          __m512d z     = _mm512_load_pd(&pz[0]);
                          __m512d stht,ear,eai,cer,cei;
                          __m512d k,vtht;
                         std::complex<double> fth;
                          double sumr,sumi;
                         vtht= _mm512_set1_pd(tht);
                         k   = _mm512_mul_pd(C314159265358979323846264338328,
                                                          _mm512_set1_pd(gam));
                         stht= xsin(tht);
                         ear = _mm512_setzero_pd();
                         eai = _mm512_fmadd_pd(stht,
                                           _mm512_mul_pd(k,z),phiz);
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         sumr= 0.0;
                         cer = _mm512_mul_pd(Az,cer);
                         _mm512_store_pd(&intr[0],cer);
                         sumi= 0.0;
                         cei = _mm512_mul_pd(Az,cei);
                         _mm512_store_pd(&inti[0],cei);
                         cspint(NTAB,&pz[0],&intr[0],-L2,L2,&Y1[0],&Y2[0],&Y3[0],&E[0],&WRK[0],sumr);
                         cspint(NTAB,&pz[0],&inti[0],-L2,L2,&Y1[0],&Y2[0],&Y3[0],&E[0],&WRK[0],sumi);
                         fth = {sumr,sumi};
                         return (fth);
                 }
                 
                 
                 
                   std::complex<double> 
                   fth_f11_zmm8r8_cspint_8e_u(const double * __restrict  pAz,
                                            const double * __restrict  pphiz,
                                            const double * __restrict  pz,
                                            const double tht,
                                            const double gam,
                                            const int32_t L2) {
                                            
                         __ATTR_ALIGN__(64) double intr[8] = {};
                         __ATTR_ALIGN__(64) double inti[8] = {}; 
                         __ATTR_ALIGN__(64) double Y1[8]   = {};
                         __ATTR_ALIGN__(64) double Y2[8]   = {};
                         __ATTR_ALIGN__(64) double Y3[8]   = {}; 
                         __ATTR_ALIGN__(64) double E[8]    = {};
                         __ATTR_ALIGN__(64) double WRK[8]  = {}; 
                         constexpr int32_t NTAB = 8;               
                         constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                          __m512d Az    = _mm512_loadu_pd(&pAz[0]);
                          __m512d phiz  = _mm512_loadu_pd(&pphiz[0]);
                          __m512d z     = _mm512_loadu_pd(&pz[0]);
                          __m512d stht,ear,eai,cer,cei;
                          __m512d k,vtht;
                         std::complex<double> fth;
                          double sumr,sumi;
                         vtht= _mm512_set1_pd(tht);
                         k   = _mm512_mul_pd(C314159265358979323846264338328,
                                                          _mm512_set1_pd(gam));
                         stht= xsin(tht);
                         ear = _mm512_setzero_pd();
                         eai = _mm512_fmadd_pd(stht,
                                           _mm512_mul_pd(k,z),phiz);
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         sumr= 0.0;
                         cer = _mm512_mul_pd(Az,cer);
                         _mm512_storeu_pd(&intr[0],cer);
                         sumi= 0.0;
                         cei = _mm512_mul_pd(Az,cei);
                         _mm512_storeu_pd(&inti[0],cei);
                         cspint(NTAB,&pz[0],&intr[0],-L2,L2,&Y1[0],&Y2[0],&Y3[0],&E[0],&WRK[0],sumr);
                         cspint(NTAB,&pz[0],&inti[0],-L2,L2,&Y1[0],&Y2[0],&Y3[0],&E[0],&WRK[0],sumi);
                         fth = {sumr,sumi};
                         return (fth);
                 }
                 
                 
                   /*
                        Formula 1.1, p. 14
                        Integrator avint (irregular abscissas) single data vector ZMM.
                   */
                 
                 
                  
                   std::complex<double> 
                   fth_f11_zmm8r8_avint_8e_a(const double * __restrict __ATTR_ALIGN__(64) pAz,
                                               const double * __restrict __ATTR_ALIGN__(64) pphiz,
                                               const double * __restrict __ATTR_ALIGN__(64) pz,
                                               const double tht,
                                               const double gam,
                                               const int32_t L2,
                                               int32_t & ierr,
                                               int32_t & ieri) {
                                            
                         __ATTR_ALIGN__(64) double intr[8] = {};
                         __ATTR_ALIGN__(64) double inti[8] = {}; 
                         
                                      
                         constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                          __m512d Az    = _mm512_load_pd(&pAz[0]);
                          __m512d phiz  = _mm512_load_pd(&pphiz[0]);
                          __m512d z     = _mm512_load_pd(&pz[0]);
                          __m512d stht,ear,eai,cer,cei;
                          __m512d k,vtht;
                         std::complex<double> fth;
                          double sumr,sumi;
                         int32_t err,eri;
                         vtht= _mm512_set1_pd(tht);
                         k   = _mm512_mul_pd(C314159265358979323846264338328,
                                                          _mm512_set1_pd(gam));
                         stht= xsin(tht);
                         ear = _mm512_setzero_pd();
                         eai = _mm512_fmadd_pd(stht,
                                           _mm512_mul_pd(k,z),phiz);
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         sumr= 0.0;
                         cer = _mm512_mul_pd(Az,cer);
                         _mm512_store_pd(&intr[0],cer);
                         sumi= 0.0;
                         cei = _mm512_mul_pd(Az,cei);
                         _mm512_store_pd(&inti[0],cei);
                         sumr = avint(&pz[0],&intr[0],16,-L2,L2,err);
                         sumi = avint(&pz[0],&inti[0],16,-L2,L2,eri);
                         ierr = err;
                         ieri = eri;
                         if(ierr==3 || ieri==3) {
                            const double NAN = std::numeric_limits<double>::quiet_NaN();
                            fth = {NAN,NAN};
                            return (fth);
                         }
                         fth = {sumr,sumi};
                         return (fth);
                 }
                 
                 
                   
                   std::complex<double> 
                   fth_f11_zmm8r8_avint_8e_u(const double * __restrict  pAz,
                                               const double * __restrict  pphiz,
                                               const double * __restrict  pz,
                                               const double tht,
                                               const double gam,
                                               const int32_t L2,
                                               int32_t & ierr,
                                               int32_t & ieri) {
                                            
                          double intr[8] = {};
                          double inti[8] = {}; 
                         
                                    
                         constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                          __m512d Az    = _mm512_loadu_pd(&pAz[0]);
                          __m512d phiz  = _mm512_loadu_pd(&pphiz[0]);
                          __m512d z     = _mm512_loadu_pd(&pz[0]);
                          __m512d stht,ear,eai,cer,cei;
                          __m512d k,vtht;
                         std::complex<double> fth;
                          double sumr,sumi;
                         int32_t err,eri;
                         vtht= _mm512_set1_pd(tht);
                         k   = _mm512_mul_pd(C314159265358979323846264338328,
                                                          _mm512_set1_pd(gam));
                         stht= xsin(tht);
                         ear = _mm512_setzero_pd();
                         eai = _mm512_fmadd_pd(stht,
                                           _mm512_mul_pd(k,z),phiz);
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         sumr= 0.0;
                         cer = _mm512_mul_pd(Az,cer);
                         _mm512_storeu_pd(&intr[0],cer);
                         sumi= 0.0;
                         cei = _mm512_mul_pd(Az,cei);
                         _mm512_storeu_pd(&inti[0],cei);
                         sumr = avint(&pz[0],&intr[0],16,-L2,L2,err);
                         sumi = avint(&pz[0],&inti[0],16,-L2,L2,eri);
                         ierr = err;
                         ieri = eri;
                         if(ierr==3 || ieri==3) {
                            const double NAN = std::numeric_limits<double>::quiet_NaN();
                            fth = {NAN,NAN};
                            return (fth);
                         }
                         fth = {sumr,sumi};
                         return (fth);
                 }
                 
                 
                 
                 
                   /*
                        Formula 1.1, p. 14
                        Integrator cspint n-elements data arrays.
                        Single-threaded.
                   */
                   
                  
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
                                            const int32_t NTAB) {
                                            
                        if(__builtin_expect(NTAB==16,0)) {
                           std::complex<double> fth = {0.0,0.0};
                           fth = f11_zmm8r8_cspint_8e_a(pAz,pphiz,pz,tht,gam,L2);
                           return (fth);
                        }     
                        
                        constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         __m512d vtht,stht,vk;
                         __m512d ear,cer,cei;
                         std::complex<double> fth;
                         double sumr,sumi; 
                         double sint,k;
                        int32_t i; 
                        vtht = _mm512_set1_pd(tht);
                        sint = sin(tht);
                        vk   = _mm512_mul_pd(C314159265358979323846264338328,
                                                          _mm512_set1_pd(gam));
                        ear  = _mm512_setzero_pd();
                        k    = C314159265358979323846264338328*gam;
                        stht = xsin(vtht);
                        for(i = 0; i != ROUND_TO_EIGHT(n,7); i += 8) {
                             _mm_prefetch((const char*)&pAz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pphiz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pz[i],_MM_HINT_T0);
                              __m512d t0 = _mm512_load_pd(&pz[i]);
                              __m512d t1 = _mm512_load_pd(&pphiz[i]);
                              __m512d eai= _mm512_fmadd_pd(stht,
                                                            _mm512_mul_pd(vk,t0),t1);
                              __m512d t2 = _mm512_load_pd(&pAz[i]);
                             cexp_zmm8r8(ear,eai,&cer,&cei);
                             cer = _mm512_mul_pd(cer,t2);
                             _mm512_store_pd(&intr[i],cer);
                             cei = _mm512_mul_pd(cei,t2);
                             _mm512_store_pd(&inti[i],cei);
                        }   
                        sumr = 0.0;
                        sumi = 0.0;
                        for(; i != n; ++i) {
                              double t0 = pAz[i];
                              double t1 = pphiz[i];
                              double t2 = pz[i];
                              double zz = k*t0;
                              double eai= sint*zz+t1;
                              std::complex<double> c = std::exp({0.0,eai});
                             intr[i] = c.real()*t0;
                             inti[i] = c.imag()*t0;
                        } 
                        cspint(NTAB,&pz[0],&intr[0],-L2,L2,&csd.Y1[0],
                               &csd.Y2[0],&csd.Y3[0],&csd.E[0],&csd.WRK[0],sumr); 
                        cspint(NTAB,&pz[0],&inti[0],-L2,L2, &csd.Y1[0],
                               &csd.Y2[0],&csd.Y3[0],&csd.E[0],&csd.WRK[0],sumi);  
                        fth = {sumr,sumi};
                        return (fth);            
                }
                
                
                 
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
                                            const int32_t NTAB) {
                                            
                        if(__builtin_expect(NTAB==16,0)) {
                           std::complex<double> fth = {0.0,0.0};
                           fth = f11_zmm8r8_cspint_8e_u(pAz,pphiz,pz,tht,gam,L2);
                           return (fth);
                        }     
                        
                        constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         __m512d vtht,stht,vk;
                         __m512d ear,cer,cei;
                         std::complex<double> fth;
                         double sumr,sumi; 
                         double sint,k;
                        int32_t i; 
                        vtht = _mm512_set1_pd(tht);
                        sint = sin(tht);
                        vk   = _mm512_mul_pd(C314159265358979323846264338328,
                                                          _mm512_set1_pd(gam));
                        ear  = _mm512_setzero_pd();
                        k    = C314159265358979323846264338328*gam;
                        stht = xsin(vtht);
                        for(i = 0; i != ROUND_TO_EIGHT(n,7); i += 8) {
                             _mm_prefetch((const char*)&pAz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pphiz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pz[i],_MM_HINT_T0);
                              __m512d t0 = _mm512_loadu_pd(&pz[i]);
                              __m512d t1 = _mm512_loadu_pd(&pphiz[i]);
                              __m512d eai= _mm512_fmadd_pd(stht,
                                                            _mm512_mul_pd(vk,t0),t1);
                              __m512d t2 = _mm512_loadu_pd(&pAz[i]);
                             cexp_zmm8r8(ear,eai,&cer,&cei);
                             cer = _mm512_mul_pd(cer,t2);
                             _mm512_storeu_pd(&intr[i],cer);
                             cei = _mm512_mul_pd(cei,t2);
                             _mm512_storeu_pd(&inti[i],cei);
                        }   
                        sumr = 0.0;
                        sumi = 0.0;
                        for(; i != n; ++i) {
                              double t0 = pAz[i];
                              double t1 = pphiz[i];
                              double t2 = pz[i];
                              double zz = k*t0;
                              double eai= sint*zz+t1;
                              std::complex<double> c = std::exp({0.0,eai});
                             intr[i] = c.real()*t0;
                             inti[i] = c.imag()*t0;
                        } 
                        cspint(NTAB,&pz[0],&intr[0],-L2,L2,&csd.Y1[0],
                               &csd.Y2[0],&csd.Y3[0],&csd.E[0],&csd.WRK[0],sumr); 
                        cspint(NTAB,&pz[0],&inti[0],-L2,L2, &csd.Y1[0],
                               &csd.Y2[0],&csd.Y3[0],&csd.E[0],&csd.WRK[0],sumi);  
                        fth = {sumr,sumi};
                        return (fth);            
                }
                
                 /*
                        Formula 1.1, p. 14
                        Integrator avint (irregular abscissas) n-elements data arrays.
                        Single-threaded.
                   */
                   
                  
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
                                            int32_t & ieri) {
                                            
                        if(__builtin_expect(NTAB==16,0)) {
                           std::complex<double> fth = {0.0,0.0};
                           fth = f11_zmm8r8_avint_8e_a(pAz,pphiz,pz,tht,gam,L2,ierr,ieri);
                           return (fth);
                        }     
                        
                        constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         __m512d vtht,stht,vk;
                         __m512d ear,cer,cei;
                         std::complex<double> fth;
                         double sumr,sumi; 
                         double sint,k;
                        int32_t i,err,eri; 
                        vtht = _mm512_set1_pd(tht);
                        sint = sin(tht);
                        vk   = _mm512_mul_pd(C314159265358979323846264338328,
                                                          _mm512_set1_pd(gam));
                        ear  = _mm512_setzero_pd();
                        k    = C314159265358979323846264338328*gam;
                        stht = xsin(vtht);
                        for(i = 0; i != ROUND_TO_EIGHT(n,7); i += 8) {
                             _mm_prefetch((const char*)&pAz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pphiz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pz[i],_MM_HINT_T0);
                              __m512d t0 = _mm512_load_pd(&pz[i]);
                              __m512d t1 = _mm512_load_pd(&pphiz[i]);
                              __m512d eai= _mm512_fmadd_pd(stht,
                                                            _mm512_mul_pd(vk,t0),t1);
                              __m512d t2 = _mm512_load_pd(&pAz[i]);
                             cexp_zmm8r8(ear,eai,&cer,&cei);
                             cer = _mm512_mul_pd(cer,t2);
                             _mm512_store_pd(&intr[i],cer);
                             cei = _mm512_mul_pd(cei,t2);
                             _mm512_store_pd(&inti[i],cei);
                        }   
                        sumr = 0.0;
                        sumi = 0.0;
                        for(; i != n; ++i) {
                              double t0 = pAz[i];
                              double t1 = pphiz[i];
                              double t2 = pz[i];
                              double zz = k*t0;
                              double eai= sint*zz+t1;
                              std::complex<double> c = std::exp({0.0,eai});
                             intr[i] = c.real()*t0;
                             inti[i] = c.imag()*t0;
                        } 
                        sumr = avint(&pz[0],&intr[0],n,-L2,L2,err); 
                        sumi = avint(&pz[0],&inti[0],n,-L2,L2,eri);
                        ierr = err;
                        ieri = eri;
                        if(ierr==3 || ieri==3) {
                            const double NAN = std::numeric_limits<double>::quiet_NaN();
                            fth = {NAN,NAN};
                            return (fth);
                         }
                        fth = {sumr,sumi};
                        return (fth);            
                }
                
                
                 
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
                                            int32_t & ieri) {
                                            
                        if(__builtin_expect(NTAB==16,0)) {
                           std::complex<double> fth = {0.0,0.0};
                           fth = f11_zmm8r8_avint_8e_u(pAz,pphiz,pz,tht,gam,L2,ierr,ieri);
                           return (fth);
                        }     
                        
                        constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         __m512d vtht,stht,vk;
                         __m512d ear,cer,cei;
                         std::complex<double> fth;
                         double sumr,sumi; 
                         double sint,k;
                        int32_t i,err,eri; 
                        vtht = _mm512_set1_pd(tht);
                        sint = sin(tht);
                        vk   = _mm512_mul_pd(C314159265358979323846264338328,
                                                          _mm512_set1_pd(gam));
                        ear  = _mm512_setzero_pd();
                        k    = C314159265358979323846264338328*gam;
                        stht = xsin(vtht);
                        for(i = 0; i != ROUND_TO_EIGHT(n,7); i += 8) {
                             _mm_prefetch((const char*)&pAz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pphiz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pz[i],_MM_HINT_T0);
                              __m512d t0 = _mm512_loadu_pd(&pz[i]);
                              __m512d t1 = _mm512_loadu_pd(&pphiz[i]);
                              __m512d eai= _mm512_fmadd_pd(stht,
                                                            _mm512_mul_pd(vk,t0),t1);
                              __m512d t2 = _mm512_load_pd(&pAz[i]);
                             cexp_zmm8r8(ear,eai,&cer,&cei);
                             cer = _mm512_mul_pd(cer,t2);
                             _mm512_storeu_pd(&intr[i],cer);
                             cei = _mm512_mul_pd(cei,t2);
                             _mm512_storeu_pd(&inti[i],cei);
                        }   
                        sumr = 0.0;
                        sumi = 0.0;
                        for(; i != n; ++i) {
                              double t0 = pAz[i];
                              double t1 = pphiz[i];
                              double t2 = pz[i];
                              double zz = k*t0;
                              double eai= sint*zz+t1;
                              std::complex<double> c = std::exp({0.0,eai});
                             intr[i] = c.real()*t0;
                             inti[i] = c.imag()*t0;
                        } 
                        sumr = avint(&pz[0],&intr[0],n,-L2,L2,err); 
                        sumi = avint(&pz[0],&inti[0],n,-L2,L2,eri);
                        ierr = err;
                        ieri = eri;
                        if(ierr==3 || ieri==3) {
                            const double NAN = std::numeric_limits<double>::quiet_NaN();
                            fth = {NAN,NAN};
                            return (fth);
                         }
                        fth = {sumr,sumi};
                        return (fth);            
                }
                
                
               
                          
                
                   void Bx_f13_zmm8r8_u10x_a(const double * __restrict __ATTR_ALIGN__(64) pA0,
                                              const double * __restrict __ATTR_ALIGN__(64) pA,
                                              double * __restrict __ATTR_ALIGN__(64) pB,
                                              const int32_t n) {
                                              
                        if(__builtin_expect(0==n,0)) { return;}
                         __m512d zmm0,zmm1,zmm2,zmm3;
                         __m512d zmm4,zmm5,zmm6,zmm7;
                         __m512d zmm8,zmm9,zmm10,zmm11;
                         __m512d zmm12,zmm13,zmm14,zmm15;
                        int32_t i;
                        
                        while(n && ((uintptr_t)pB & 63)) {
                              double t0;
                              double t1;
                              double bx;
                              double rat;
                              t0 = *pA0;
                              t1 = *pA;
                              rat= t1/t0;
                              bx = std::log(rat);
                             *pB++ = bx;
                              pA0++;
                              pA++;
                              n--;
                        }
                        
                        for(i = 0; (i+79) < n; i += 80) {
                             _mm_prefetch((const char*)&pA0[i+0],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+0], _MM_HINT_T0);
                             zmm0 = _mm512_load_pd(&pA0[i+0]);
                             zmm1 = _mm512_load_pd(&pA[i+0]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_store_pd(&pB[i+0],zmm3);
                             zmm4 = _mm512_load_pd(&pA0[i+8]);
                             zmm5 = _mm512_load_pd(&pA[i+8]);
                             zmm6 = _mm512_div_pd(zmm5,zmm4);
                             zmm7 = xlog(zmm6);
                             _mm512_store_pd(&pB[i+8],zmm7);
                             _mm_prefetch((const char*)&pA0[i+16],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+16], _MM_HINT_T0);
                             zmm8 = _mm512_load_pd(&pA0[i+16]);
                             zmm9 = _mm512_load_pd(&pA[i+16]);
                             zmm10 = _mm512_div_pd(zmm9,zmm8);
                             zmm11 = xlog(zmm10);
                             _mm512_store_pd(&pB[i+16],zmm10);
                             zmm12 = _mm512_load_pd(&pA0[i+24]);
                             zmm13 = _mm512_load_pd(&pA[i+24]);
                             zmm14 = _mm512_div_pd(zmm13,zmm12);
                             zmm15 = xlog(zmm14);
                             _mm512_store_pd(&pB[i+24],zmm15);
                             _mm_prefetch((const char*)&pA0[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+32], _MM_HINT_T0);
                             zmm0 = _mm512_load_pd(&pA0[i+32]);
                             zmm1 = _mm512_load_pd(&pA[i+32]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_store_pd(&pB[i+32],zmm3);
                             zmm4 = _mm512_load_pd(&pA0[i+40]);
                             zmm5 = _mm512_load_pd(&pA[i+40]);
                             zmm6 = _mm512_div_pd(zmm5,zmm4);
                             zmm7 = xlog(zmm6);
                             _mm512_store_pd(&pB[i+40],zmm7);
                             _mm_prefetch((const char*)&pA0[i+48],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+48], _MM_HINT_T0);
                             zmm8 = _mm512_load_pd(&pA0[i+48]);
                             zmm9 = _mm512_load_pd(&pA[i+48]);
                             zmm10= _mm512_div_pd(zmm9,zmm8);
                             zmm11= xlog(zmm10);
                             _mm512_store_pd(&pB[i+48],zmm11);
                             zmm12 = _mm512_load_pd(&pA0[i+56]);
                             zmm13 = _mm512_load_pd(&pA[i+56]);
                             zmm14= _mm512_div_pd(zmm13,zmm12);
                             zmm15= xlog(zmm14);
                             _mm512_store_pd(&pB[i+56],zmm15);
                             _mm_prefetch((const char*)&pA0[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+64], _MM_HINT_T0);
                             zmm0 = _mm512_load_pd(&pA0[i+64]);
                             zmm1 = _mm512_load_pd(&pA[i+64]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_store_pd(&pB[i+64],zmm3);
                             zmm4 = _mm512_load_pd(&pA0[i+72]);
                             zmm5 = _mm512_load_pd(&pA[i+72]);
                             zmm6 = _mm512_div_pd(zmm5,zmm4);
                             zmm7 = xlog(zmm6);
                             _mm512_store_pd(&pB[i+72],zmm7);
                        }   
                        
                        for(; (i+63) < n; i += 64) {
                             zmm0 = _mm512_load_pd(&pA0[i+0]);
                             zmm1 = _mm512_load_pd(&pA[i+0]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_store_pd(&pB[i+0],zmm3);
                             zmm4 = _mm512_load_pd(&pA0[i+8]);
                             zmm5 = _mm512_load_pd(&pA[i+8]);
                             zmm6 = _mm512_div_pd(zmm5,zmm4);
                             zmm7 = xlog(zmm6);
                             _mm512_store_pd(&pB[i+8],zmm7);
                             zmm8 = _mm512_load_pd(&pA0[i+16]);
                             zmm9 = _mm512_load_pd(&pA[i+16]);
                             zmm10 = _mm512_div_pd(zmm9,zmm8);
                             zmm11 = xlog(zmm10);
                             _mm512_store_pd(&pB[i+16],zmm10);
                             zmm12 = _mm512_load_pd(&pA0[i+24]);
                             zmm13 = _mm512_load_pd(&pA[i+24]);
                             zmm14 = _mm512_div_pd(zmm13,zmm12);
                             zmm15 = xlog(zmm14);
                             _mm512_store_pd(&pB[i+24],zmm15);
                             zmm0 = _mm512_load_pd(&pA0[i+32]);
                             zmm1 = _mm512_load_pd(&pA[i+32]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_store_pd(&pB[i+32],zmm3);
                             zmm4 = _mm512_load_pd(&pA0[i+40]);
                             zmm5 = _mm512_load_pd(&pA[i+40]);
                             zmm6 = _mm512_div_pd(zmm5,zmm4);
                             zmm7 = xlog(zmm6);
                             _mm512_store_pd(&pB[i+40],zmm7);
                             zmm8 = _mm512_load_pd(&pA0[i+48]);
                             zmm9 = _mm512_load_pd(&pA[i+48]);
                             zmm10= _mm512_div_pd(zmm9,zmm8);
                             zmm11= xlog(zmm10);
                             _mm512_store_pd(&pB[i+48],zmm11);
                             zmm12 = _mm512_load_pd(&pA0[i+56]);
                             zmm13 = _mm512_load_pd(&pA[i+56]);
                             zmm14= _mm512_div_pd(zmm13,zmm12);
                             zmm15= xlog(zmm14);
                             _mm512_store_pd(&pB[i+56],zmm15);
                        }
                        
                                                         
                        for(; (i+31) < n; i += 32) {
                             zmm0 = _mm512_load_pd(&pA0[i+0]);
                             zmm1 = _mm512_load_pd(&pA[i+0]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_store_pd(&pB[i+0],zmm3);
                             zmm4 = _mm512_load_pd(&pA0[i+8]);
                             zmm5 = _mm512_load_pd(&pA[i+8]);
                             zmm6 = _mm512_div_pd(zmm5,zmm4);
                             zmm7 = xlog(zmm6);
                             _mm512_store_pd(&pB[i+8],zmm7);
                             zmm8 = _mm512_load_pd(&pA0[i+16]);
                             zmm9 = _mm512_load_pd(&pA[i+16]);
                             zmm10 = _mm512_div_pd(zmm9,zmm8);
                             zmm11 = xlog(zmm10);
                             _mm512_store_pd(&pB[i+16],zmm10);
                             zmm12 = _mm512_load_pd(&pA0[i+24]);
                             zmm13 = _mm512_load_pd(&pA[i+24]);
                             zmm14 = _mm512_div_pd(zmm13,zmm12);
                             zmm15 = xlog(zmm14);
                             _mm512_store_pd(&pB[i+24],zmm15);
                        } 
                        
                        for(; (i+15) < n; i += 16) {
                             zmm0 = _mm512_load_pd(&pA0[i+0]);
                             zmm1 = _mm512_load_pd(&pA[i+0]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_store_pd(&pB[i+0],zmm3);
                             zmm4 = _mm512_load_pd(&pA0[i+8]);
                             zmm5 = _mm512_load_pd(&pA[i+8]);
                             zmm6 = _mm512_div_pd(zmm5,zmm4);
                             zmm7 = xlog(zmm6);
                             _mm512_store_pd(&pB[i+8],zmm7); 
                        } 
                        
                        for(; (i+7) < n; i += 8) {
                             zmm0 = _mm512_load_pd(&pA0[i+0]);
                             zmm1 = _mm512_load_pd(&pA[i+0]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_store_pd(&pB[i+0],zmm3);
                        }      
                        
                        for(; (i+0) < n; i += 1) {
                              double t0 = pA0[i];
                              double t1 = pA[i];
                              double rat= t1/t0;
                              double bx = std::log(rat);
                              pB[i] = bx;
                        }           
                }
                
                
                
                  
                   void Bx_f13_zmm8r8_u10x_u(const double * __restrict pA0,
                                              const double * __restrict  pA,
                                              double * __restrict  pB,
                                              const int32_t n) {
                                              
                        if(__builtin_expect(0==n,0)) { return;}
                         __m512d zmm0,zmm1,zmm2,zmm3;
                         __m512d zmm4,zmm5,zmm6,zmm7;
                         __m512d zmm8,zmm9,zmm10,zmm11;
                         __m512d zmm12,zmm13,zmm14,zmm15;
                        int32_t i;
                        
                        for(i = 0; (i+79) < n; i += 80) {
                             _mm_prefetch((const char*)&pA0[i+0],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+0], _MM_HINT_T0);
                             zmm0 = _mm512_loadu_pd(&pA0[i+0]);
                             zmm1 = _mm512_loadu_pd(&pA[i+0]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_storeu_pd(&pB[i+0],zmm3);
                             zmm4 = _mm512_loadu_pd(&pA0[i+8]);
                             zmm5 = _mm512_loadu_pd(&pA[i+8]);
                             zmm6 = _mm512_div_pd(zmm5,zmm4);
                             zmm7 = xlog(zmm6);
                             _mm512_storeu_pd(&pB[i+8],zmm7);
                             _mm_prefetch((const char*)&pA0[i+16],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+16], _MM_HINT_T0);
                             zmm8 = _mm512_loadu_pd(&pA0[i+16]);
                             zmm9 = _mm512_loadu_pd(&pA[i+16]);
                             zmm10 = _mm512_div_pd(zmm9,zmm8);
                             zmm11 = xlog(zmm10);
                             _mm512_storeu_pd(&pB[i+16],zmm10);
                             zmm12 = _mm512_loadu_pd(&pA0[i+24]);
                             zmm13 = _mm512_loadu_pd(&pA[i+24]);
                             zmm14 = _mm512_div_pd(zmm13,zmm12);
                             zmm15 = xlog(zmm14);
                             _mm512_storeu_pd(&pB[i+24],zmm15);
                             _mm_prefetch((const char*)&pA0[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+32], _MM_HINT_T0);
                             zmm0 = _mm512_loadu_pd(&pA0[i+32]);
                             zmm1 = _mm512_loadu_pd(&pA[i+32]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_storeu_pd(&pB[i+32],zmm3);
                             zmm4 = _mm512_loadu_pd(&pA0[i+40]);
                             zmm5 = _mm512_loadu_pd(&pA[i+40]);
                             zmm6 = _mm512_div_pd(zmm5,zmm4);
                             zmm7 = xlog(zmm6);
                             _mm512_storeu_pd(&pB[i+40],zmm7);
                             _mm_prefetch((const char*)&pA0[i+48],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+48], _MM_HINT_T0);
                             zmm8 = _mm512_loadu_pd(&pA0[i+48]);
                             zmm9 = _mm512_loadu_pd(&pA[i+48]);
                             zmm10= _mm512_div_pd(zmm9,zmm8);
                             zmm11= xlog(zmm10);
                             _mm512_storeu_pd(&pB[i+48],zmm11);
                             zmm12 = _mm512_loadu_pd(&pA0[i+56]);
                             zmm13 = _mm512_loadu_pd(&pA[i+56]);
                             zmm14= _mm512_div_pd(zmm13,zmm12);
                             zmm15= xlog(zmm14);
                             _mm512_store_pd(&pB[i+56],zmm15);
                             _mm_prefetch((const char*)&pA0[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+64], _MM_HINT_T0);
                             zmm0 = _mm512_loadu_pd(&pA0[i+64]);
                             zmm1 = _mm512_loadu_pd(&pA[i+64]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_storeu_pd(&pB[i+64],zmm3);
                             zmm4 = _mm512_loadu_pd(&pA0[i+72]);
                             zmm5 = _mm512_loadu_pd(&pA[i+72]);
                             zmm6 = _mm512_div_pd(zmm5,zmm4);
                             zmm7 = xlog(zmm6);
                             _mm512_storeu_pd(&pB[i+72],zmm7);
                        }   
                        
                        for(; (i+63) < n; i += 64) {
                             zmm0 = _mm512_load_pd(&pA0[i+0]);
                             zmm1 = _mm512_load_pd(&pA[i+0]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_store_pd(&pB[i+0],zmm3);
                             zmm4 = _mm512_load_pd(&pA0[i+8]);
                             zmm5 = _mm512_load_pd(&pA[i+8]);
                             zmm6 = _mm512_div_pd(zmm5,zmm4);
                             zmm7 = xlog(zmm6);
                             _mm512_store_pd(&pB[i+8],zmm7);
                             zmm8 = _mm512_load_pd(&pA0[i+16]);
                             zmm9 = _mm512_load_pd(&pA[i+16]);
                             zmm10 = _mm512_div_pd(zmm9,zmm8);
                             zmm11 = xlog(zmm10);
                             _mm512_store_pd(&pB[i+16],zmm10);
                             zmm12 = _mm512_load_pd(&pA0[i+24]);
                             zmm13 = _mm512_load_pd(&pA[i+24]);
                             zmm14 = _mm512_div_pd(zmm13,zmm12);
                             zmm15 = xlog(zmm14);
                             _mm512_store_pd(&pB[i+24],zmm15);
                             zmm0 = _mm512_load_pd(&pA0[i+32]);
                             zmm1 = _mm512_load_pd(&pA[i+32]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_store_pd(&pB[i+32],zmm3);
                             zmm4 = _mm512_load_pd(&pA0[i+40]);
                             zmm5 = _mm512_load_pd(&pA[i+40]);
                             zmm6 = _mm512_div_pd(zmm5,zmm4);
                             zmm7 = xlog(zmm6);
                             _mm512_store_pd(&pB[i+40],zmm7);
                             zmm8 = _mm512_load_pd(&pA0[i+48]);
                             zmm9 = _mm512_load_pd(&pA[i+48]);
                             zmm10= _mm512_div_pd(zmm9,zmm8);
                             zmm11= xlog(zmm10);
                             _mm512_store_pd(&pB[i+48],zmm11);
                             zmm12 = _mm512_load_pd(&pA0[i+56]);
                             zmm13 = _mm512_load_pd(&pA[i+56]);
                             zmm14= _mm512_div_pd(zmm13,zmm12);
                             zmm15= xlog(zmm14);
                             _mm512_store_pd(&pB[i+56],zmm15);
                        }
                        
                                                         
                        for(; (i+31) < n; i += 32) {
                             zmm0 = _mm512_load_pd(&pA0[i+0]);
                             zmm1 = _mm512_load_pd(&pA[i+0]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_store_pd(&pB[i+0],zmm3);
                             zmm4 = _mm512_load_pd(&pA0[i+8]);
                             zmm5 = _mm512_load_pd(&pA[i+8]);
                             zmm6 = _mm512_div_pd(zmm5,zmm4);
                             zmm7 = xlog(zmm6);
                             _mm512_store_pd(&pB[i+8],zmm7);
                             zmm8 = _mm512_load_pd(&pA0[i+16]);
                             zmm9 = _mm512_load_pd(&pA[i+16]);
                             zmm10 = _mm512_div_pd(zmm9,zmm8);
                             zmm11 = xlog(zmm10);
                             _mm512_store_pd(&pB[i+16],zmm10);
                             zmm12 = _mm512_load_pd(&pA0[i+24]);
                             zmm13 = _mm512_load_pd(&pA[i+24]);
                             zmm14 = _mm512_div_pd(zmm13,zmm12);
                             zmm15 = xlog(zmm14);
                             _mm512_store_pd(&pB[i+24],zmm15);
                        } 
                        
                        for(; (i+15) < n; i += 16) {
                             zmm0 = _mm512_load_pd(&pA0[i+0]);
                             zmm1 = _mm512_load_pd(&pA[i+0]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_store_pd(&pB[i+0],zmm3);
                             zmm4 = _mm512_load_pd(&pA0[i+8]);
                             zmm5 = _mm512_load_pd(&pA[i+8]);
                             zmm6 = _mm512_div_pd(zmm5,zmm4);
                             zmm7 = xlog(zmm6);
                             _mm512_store_pd(&pB[i+8],zmm7); 
                        } 
                        
                        for(; (i+7) < n; i += 8) {
                             zmm0 = _mm512_load_pd(&pA0[i+0]);
                             zmm1 = _mm512_load_pd(&pA[i+0]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_store_pd(&pB[i+0],zmm3);
                        }      
                        
                        for(; (i+0) < n; i += 1) {
                              double t0 = pA0[i];
                              double t1 = pA[i];
                              double rat= t1/t0;
                              double bx = std::log(rat);
                              pB[i] = bx;
                        }           
                }
                
                
                
               
                   void Bx_f13_zmm8r8_u2x_u(const double * __restrict pA0,
                                              const double * __restrict  pA,
                                              double * __restrict  pB,
                                              const int32_t n) {
                                              
                        if(__builtin_expect(0==n,0)) { return;}
                         __m512d zmm0,zmm1,zmm2,zmm3;
                         __m512d zmm4,zmm5,zmm6,zmm7;
                        int32_t i;
                        
                        for(i = 0; (i+15) < n; i += 16) {
                             _mm_prefetch((const char*)&pA0[i+0],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+0], _MM_HINT_T0);
                             zmm0 = _mm512_loadu_pd(&pA0[i+0]);
                             zmm1 = _mm512_loadu_pd(&pA[i+0]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_storeu_pd(&pB[i+0],zmm3);
                             zmm4 = _mm512_loadu_pd(&pA0[i+8]);
                             zmm5 = _mm512_loadu_pd(&pA[i+8]);
                             zmm6 = _mm512_div_pd(zmm5,zmm4);
                             zmm7 = xlog(zmm6);
                             _mm512_storeu_pd(&pB[i+8],zmm7);
                       }  
                       
                        for(; (i+7) < n; i += 8) {
                             zmm0 = _mm512_loadu_pd(&pA0[i+0]);
                             zmm1 = _mm512_loadu_pd(&pA[i+0]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_storeu_pd(&pB[i+0],zmm3);
                        }      
                        
                        for(; (i+0) < n; i += 1) {
                              double t0 = pA0[i];
                              double t1 = pA[i];
                              double rat= t1/t0;
                              double bx = std::log(rat);
                             pB[i] = bx;
                        }        
                        
                 }
                 
                 
                  
                   void Bx_f13_zmm8r8_u2x_a(const double * __restrict __ATTR_ALIGN__(64) pA0,
                                              const double * __restrict __ATTR_ALIGN__(64) pA,
                                              double * __restrict __ATTR_ALIGN__(64) pB,
                                              const int32_t n) {
                                              
                        if(__builtin_expect(0==n,0)) { return;}
                         __m512d zmm0,zmm1,zmm2,zmm3;
                         __m512d zmm4,zmm5,zmm6,zmm7;
                        int32_t i;
                        
                         while(n && ((uintptr_t)pB & 63)) {
                              double t0;
                              double t1;
                              double bx;
                              double rat;
                              t0 = *pA0;
                              t1 = *pA;
                              rat= t1/t0;
                              bx = std::log(rat);
                             *pB++ = bx;
                              pA0++;
                              pA++;
                              n--;
                        }
                        
                        for(i = 0; (i+15) < n; i += 16) {
                             _mm_prefetch((const char*)&pA0[i+0],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+0], _MM_HINT_T0);
                             zmm0 = _mm512_load_pd(&pA0[i+0]);
                             zmm1 = _mm512_load_pd(&pA[i+0]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_store_pd(&pB[i+0],zmm3);
                             zmm4 = _mm512_load_pd(&pA0[i+8]);
                             zmm5 = _mm512_load_pd(&pA[i+8]);
                             zmm6 = _mm512_div_pd(zmm5,zmm4);
                             zmm7 = xlog(zmm6);
                             _mm512_store_pd(&pB[i+8],zmm7);
                       }  
                       
                        for(; (i+7) < n; i += 8) {
                             zmm0 = _mm512_load_pd(&pA0[i+0]);
                             zmm1 = _mm512_load_pd(&pA[i+0]);
                             zmm2 = _mm512_div_pd(zmm1,zmm0);
                             zmm3 = xlog(zmm2);
                             _mm512_store_pd(&pB[i+0],zmm3);
                        }      
                        
                        for(; (i+0) < n; i += 1) {
                              double t0 = pA0[i];
                              double t1 = pA[i];
                              double rat= t1/t0;
                              double bx = std::log(rat);
                             pB[i] = bx;
                        }        
                        
                 }
                 
                 
                 /*
                      Helper formula serving the needs of 
                      formula 1.4, p. 15
                      The mathematical expectation of random 
                      variable B(x).
                      Integrator 'cspint'.
                      Short data vector (ZMM-size).
                 */
                 
                 
                   double Ex_Bx_zmm8r8_cspint_8e_a(const double * __restrict __ATTR_ALIGN__(64) pBx,
                                                    const double * __restrict __ATTR_ALIGN__(64) ppdf,
                                                    const double * __restrict __ATTR_ALIGN__(64) px,
                                                    const double a,
                                                    const double b) {
                                                    
                         __ATTR_ALIGN__(64) double intr[8] = {};
                         __ATTR_ALIGN__(64) double Y1[8]   = {};
                         __ATTR_ALIGN__(64) double Y2[8]   = {};
                         __ATTR_ALIGN__(64) double Y3[8]   = {}; 
                         __ATTR_ALIGN__(64) double E[8]    = {};
                         __ATTR_ALIGN__(64) double WRK[8]  = {}; 
                         constexpr int32_t NTAB = 8;   
                          __m512d Bx = _mm512_load_pd(&pBx[0]);
                          __m512d pdf= _mm512_load_pd(&ppdf[0]);
                          __m512d prod;
                         double sum;
                         prod = _mm512_mul_pd(Bx,pdf);
                         sum  = 0.0;
                         _mm512_store_pd(&intr[0],prod);
                         cspint(NTAB,&px[0],&intr[0],a,b,&Y1[0],
                                &Y2[0],&Y3[0],&E[0],&WRK[0],sum);
                         return (sum);                                          
                }
                
                
                
                
                
                 
                   double Ex_Bx_zmm8r8_cspint_8e_u(const double * __restrict  pBx,
                                                    const double * __restrict  ppdf,
                                                    const double * __restrict  px,
                                                    const double a,
                                                    const double b) {
                                                    
                         __ATTR_ALIGN__(64) double intr[8] = {};
                         __ATTR_ALIGN__(64) double Y1[8]   = {};
                         __ATTR_ALIGN__(64) double Y2[8]   = {};
                         __ATTR_ALIGN__(64) double Y3[8]   = {}; 
                         __ATTR_ALIGN__(64) double E[8]    = {};
                         __ATTR_ALIGN__(64) double WRK[8]  = {}; 
                         constexpr int32_t NTAB = 8;   
                          __m512d Bx = _mm512_loadu_pd(&pBx[0]);
                          __m512d pdf= _mm512_loadu_pd(&ppdf[0]);
                          __m512d prod;
                         double sum;
                         prod = _mm512_mul_pd(Bx,pdf);
                         sum  = 0.0;
                         _mm512_storeu_pd(&intr[0],prod);
                         cspint(NTAB,&px[0],&intr[0],a,b,&Y1[0],
                                &Y2[0],&Y3[0],&E[0],&WRK[0],sum);
                         return (sum);                                          
                }
                
                
                  /*
                      Helper formula serving the needs of 
                      formula 1.4, p. 15
                      The mathematical expectation of random 
                      variable B(x).
                      Integrator 'avint' irregular abscissas.
                      Short data vector (ZMM-size).
                 */
                 
                 
                 
                   double Ex_Bx_zmm8r8_avint_8e_a(const double * __restrict __ATTR_ALIGN__(64) pBx,
                                                    const double * __restrict __ATTR_ALIGN__(64) ppdf,
                                                    const double * __restrict __ATTR_ALIGN__(64) px,
                                                    const double a,
                                                    const double b,
                                                    int32_t & ier){
                                                  
                                                    
                         __ATTR_ALIGN__(64) double intr[8] = {};
                                                 
                          __m512d Bx = _mm512_load_pd(&pBx[0]);
                          __m512d pdf= _mm512_load_pd(&ppdf[0]);
                          __m512d prod;
                         double sum;
                         int32_t err;
                         prod = _mm512_mul_pd(Bx,pdf);
                         sum  = 0.0;
                         _mm512_store_pd(&intr[0],prod);
                         sum = avint(&px[0],&intr[0],16,a,b,err);
                         ier = err;
                         if(ier==3) {
                            sum = std::numeric_limits<double>::quiet_NaN();
                            return (sum);
                         }
                         return (sum);                                          
                }
                
                
                 
                   double Ex_Bx_zmm8r8_avint_8e_u(const double * __restrict  pBx,
                                                    const double * __restrict  ppdf,
                                                    const double * __restrict  px,
                                                    const double a,
                                                    const double b,
                                                    int32_t & ier){
                                                  
                                                    
                         __ATTR_ALIGN__(64) double intr[8] = {};
                                                 
                          __m512d Bx = _mm512_loadu_pd(&pBx[0]);
                          __m512d pdf= _mm512_loadu_pd(&ppdf[0]);
                          __m512d prod;
                         double sum;
                         int32_t err;
                         prod = _mm512_mul_pd(Bx,pdf);
                         sum  = 0.0;
                         _mm512_storeu_pd(&intr[0],prod);
                         sum = avint(&px[0],&intr[0],16,a,b,err);
                         ier = err;
                         if(ier==3) {
                            sum = std::numeric_limits<double>::quiet_NaN();
                            return (sum);
                         }
                         return (sum);                                          
                }
                
                 /*
                      Helper formula serving the needs of 
                      formula 1.4, p. 15
                      The mathematical expectation of random 
                      variable B(x).
                      Integrator 'cspint'.
                      Long data vector (arrays).
                 */
                 
                 
                 
                 
                  
                   double Ex_Bx_zmm8r8_cspint_ne_a(const double * __restrict __ATTR_ALIGN__(64) pBx,
                                                   const double * __restrict __ATTR_ALIGN__(64) ppdf,
                                                   const double * __restrict __ATTR_ALIGN__(64) px,
                                                   double * __restrict __ATTR_ALIGN__(64)       intr,
                                                   CSPINT_DATA_1T csd,
                                                   const double a,
                                                   const double b,
                                                   const int32_t NTAB) {
                                                   
                        if(__buitlin_expect(NTAB==16,0)) {
                             double sum = 0.0;
                            sum = Ex_Bx_zmm8r8_cspint_8e_a(pBx,ppdf,px,a,b);
                            return (sum);
                        }                
                        
                         __m512d Bx,pdf,prod;
                         double sum;
                        int32_t i;
                        
                        for(i = 0; i != ROUND_TO_EIGHT(n,7); i += 8) {
                             _mm_prefetch((const char*)&pBx[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&ppdf[i],_MM_HINT_T0);
                             Bx = _mm512_load_pd(&pBx[i]);
                             pdf= _mm512_load_pd(&ppdf[i]);
                             prod = _mm512_mul_pd(Bx,pdf);
                             _mm512_store_pd(&intr[i], prod); 
                        }  
                        sum = 0.0;
                        for(; i != n; ++i) {
                              double Bx = pBx[i];
                              double pdf= ppdf[i];
                              double prod = Bx*pdf;
                             intr[i] = prod;
                        }  
                        
                        cspint(NTAB,&px[0],&intr[0],a,b,&csd.Y1[0],
                                &csd.Y2[0],&csd.Y3[0],&csd.E[0],&csd.WRK[0],sum);
                         return (sum);             
                  }
                  
                  
                  
                 
                   double Ex_Bx_zmm8r8_cspint_ne_u(const double * __restrict  pBx,
                                                   const double * __restrict  ppdf,
                                                   const double * __restrict  px,
                                                   double * __restrict intr,
                                                   CSPINT_DATA_1T csd,
                                                   const double a,
                                                   const double b,
                                                   const int32_t NTAB) {
                                                   
                        if(__buitlin_expect(NTAB==16,0)) {
                             double sum = 0.0;
                            sum = Ex_Bx_zmm8r8_cspint_8e_u(pBx,ppdf,px,a,b);
                            return (sum);
                        }                
                        
                         __m512d Bx,pdf,prod;
                         double sum;
                        int32_t i;
                        
                        for(i = 0; i != ROUND_TO_EIGHT(n,7); i += 8) {
                             _mm_prefetch((const char*)&pBx[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&ppdf[i],_MM_HINT_T0);
                             Bx = _mm512_loadu_pd(&pBx[i]);
                             pdf= _mm512_loadu_pd(&ppdf[i]);
                             prod = _mm512_mul_pd(Bx,pdf);
                             _mm512_storeu_pd(&intr[i], prod); 
                        }  
                        sum = 0.0;
                        for(; i != n; ++i) {
                              double Bx = pBx[i];
                              double pdf= ppdf[i];
                              double prod = Bx*pdf;
                             intr[i] = prod;
                        }  
                        
                        cspint(NTAB,&px[0],&intr[0],a,b,&csd.Y1[0],
                                &csd.Y2[0],&csd.Y3[0],&csd.E[0],&csd.WRK[0],sum);
                         return (sum);             
                  }
                  
                  
                   /*
                      Helper formula serving the needs of 
                      formula 1.4, p. 15
                      The mathematical expectation of random 
                      variable B(x).
                      Integrator 'avint' irregular abscissas.
                      Long data vector (arrays).
                 */
                 
                 
                 
                   double Ex_Bx_zmm8r8_avint_ne_a(const double * __restrict __ATTR_ALIGN__(64) pBx,
                                                   const double * __restrict __ATTR_ALIGN__(64) ppdf,
                                                   const double * __restrict __ATTR_ALIGN__(64) px,
                                                   double * __restrict __ATTR_ALIGN__(64)       intr,
                                                   const int32_t n,
                                                   const double a,
                                                   const double b,
                                                   int32_t & ier) {
                                                   
                        if(__buitlin_expect(n==16,0)) {
                             double sum = 0.0;
                            sum = Ex_Bx_zmm8r8_avint_8e_a(pBx,ppdf,px,a,b,ier);
                            return (sum);
                        }                
                        
                         __m512d Bx,pdf,prod;
                         double sum;
                        int32_t i,err;
                        
                        for(i = 0; i != ROUND_TO_EIGHT(n,7); i += 8) {
                             _mm_prefetch((const char*)&pBx[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&ppdf[i],_MM_HINT_T0);
                             Bx = _mm512_load_pd(&pBx[i]);
                             pdf= _mm512_load_pd(&ppdf[i]);
                             prod = _mm512_mul_pd(Bx,pdf);
                             _mm512_store_pd(&intr[i], prod); 
                        }  
                        sum = 0.0;
                        for(; i != n; ++i) {
                              double Bx = pBx[i];
                              double pdf= ppdf[i];
                              double prod = Bx*pdf;
                             intr[i] = prod;
                        }  
                        
                        sum = avint(&px[0],&intr[0],n,a,b,err);
                        ier = err;
                        if(ier==3) {
                           sum = std::numeric_limits<double>::quiet_NaN();
                           return (sum);
                        }
                         return (sum);             
                  }
                  
                  
                  
                 
                   double Ex_Bx_zmm8r8_avint_ne_u(const double * __restrict  pBx,
                                                   const double * __restrict  ppdf,
                                                   const double * __restrict  px,
                                                   double * __restrict intr,
                                                   const int32_t n,
                                                   const double a,
                                                   const double b,
                                                   int32_t & ier) {
                                                   
                        if(__buitlin_expect(n==16,0)) {
                             double sum = 0.0;
                            sum = Ex_Bx_zmm8r8_avint_8e_u(pBx,ppdf,px,a,b,ier);
                            return (sum);
                        }                
                        
                         __m512d Bx,pdf,prod;
                         double sum;
                        int32_t i,err;
                        
                        for(i = 0; i != ROUND_TO_EIGHT(n,7); i += 8) {
                             _mm_prefetch((const char*)&pBx[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&ppdf[i],_MM_HINT_T0);
                             Bx = _mm512_loadu_pd(&pBx[i]);
                             pdf= _mm512_loadu_pd(&ppdf[i]);
                             prod = _mm512_mul_pd(Bx,pdf);
                             _mm512_storeu_pd(&intr[i], prod); 
                        }  
                        sum = 0.0;
                        for(; i != n; ++i) {
                              double Bx = pBx[i];
                              double pdf= ppdf[i];
                              double prod = Bx*pdf;
                             intr[i] = prod;
                        }  
                        
                        sum = avint(&px[0],&intr[0],n,a,b,err);
                        ier = err;
                        if(ier==3) {
                           sum = std::numeric_limits<double>::quiet_NaN();
                           return (sum);
                        }
                         return (sum);             
                  }
                  
                  
                    
                 /*
                      Helper formula serving the needs of 
                      formula 1.4, p. 15
                      The mathematical expectation of phase random 
                      variable phi(x).
                      Integrator 'cspint'.
                      Short data vector (ZMM-size).
                 */
                 
                 
                  
                   double Ex_phix_zmm8r8_cspint_8e_a(const double * __restrict __ATTR_ALIGN__(64) pphix,
                                                      const double * __restrict __ATTR_ALIGN__(64) ppdf,
                                                      const double * __restrict __ATTR_ALIGN__(64) px,
                                                      const double a,
                                                      const double b) {
                                                    
                         __ATTR_ALIGN__(64) double intr[8] = {};
                         __ATTR_ALIGN__(64) double Y1[8]   = {};
                         __ATTR_ALIGN__(64) double Y2[8]   = {};
                         __ATTR_ALIGN__(64) double Y3[8]   = {}; 
                         __ATTR_ALIGN__(64) double E[8]    = {};
                         __ATTR_ALIGN__(64) double WRK[8]  = {}; 
                         constexpr int32_t NTAB = 8;   
                          __m512d phix = _mm512_load_pd(&pphix[0]);
                          __m512d pdf  = _mm512_load_pd(&ppdf[0]);
                          __m512d prod;
                         double sum;
                         prod = _mm512_mul_pd(phix,pdf);
                         sum  = 0.0;
                         _mm512_store_pd(&intr[0],prod);
                         cspint(NTAB,&px[0],&intr[0],a,b,&Y1[0],
                                &Y2[0],&Y3[0],&E[0],&WRK[0],sum);
                         return (sum);                                          
                }
                
                
                
                  
                   double Ex_phix_zmm8r8_cspint_8e_u(const double * __restrict  pphix,
                                                      const double * __restrict  ppdf,
                                                      const double * __restrict  px,
                                                      const double a,
                                                      const double b) {
                                                    
                         __ATTR_ALIGN__(64) double intr[8] = {};
                         __ATTR_ALIGN__(64) double Y1[8]   = {};
                         __ATTR_ALIGN__(64) double Y2[8]   = {};
                         __ATTR_ALIGN__(64) double Y3[8]   = {}; 
                         __ATTR_ALIGN__(64) double E[8]    = {};
                         __ATTR_ALIGN__(64) double WRK[8]  = {}; 
                         constexpr int32_t NTAB = 8;   
                          __m512d phix = _mm512_loadu_pd(&pphix[0]);
                          __m512d pdf  = _mm512_loadu_pd(&ppdf[0]);
                          __m512d prod;
                         double sum;
                         prod = _mm512_mul_pd(phix,pdf);
                         sum  = 0.0;
                         _mm512_storeu_pd(&intr[0],prod);
                         cspint(NTAB,&px[0],&intr[0],a,b,&Y1[0],
                                &Y2[0],&Y3[0],&E[0],&WRK[0],sum);
                         return (sum);                                          
                }
                
                
                  /*
                      Helper formula serving the needs of 
                      formula 1.4, p. 15
                      The mathematical expectation of phase random 
                      variable phi(x).
                      Integrator 'avint' irregular abscissas.
                      Short data vector (ZMM-size).
                 */
                 
                 
                  
                   double Ex_phix_zmm8r8_avint_8e_a(const double * __restrict __ATTR_ALIGN__(64) pphix,
                                                      const double * __restrict __ATTR_ALIGN__(64) ppdf,
                                                      const double * __restrict __ATTR_ALIGN__(64) px,
                                                      const double a,
                                                      const double b,
                                                      int32_t & ier) {
                                                    
                         __ATTR_ALIGN__(64) double intr[8] = {};
                         constexpr int32_t n = 16;   
                          __m512d phix = _mm512_load_pd(&pphix[0]);
                          __m512d pdf  = _mm512_load_pd(&ppdf[0]);
                          __m512d prod;
                          double sum;
                         int32_t err;
                         err = -1;
                         prod = _mm512_mul_pd(phix,pdf);
                         sum  = 0.0;
                         _mm512_store_pd(&intr[0],prod);
                         sum = cspint(&px[0],&intr[0],n,a,b,err);
                         ier = err;
                         if(ier==3) {
                             sum = std::numeric_limits<double>::quiet_NaN();
                             return (sum);
                         }      
                         return (sum);                                          
                }
                
                
                 
                   double Ex_phix_zmm8r8_avint_8e_u(const double * __restrict pphix,
                                                      const double * __restrict ppdf,
                                                      const double * __restrict  px,
                                                      const double a,
                                                      const double b,
                                                      int32_t & ier) {
                                                    
                         double intr[8] = {};
                         constexpr int32_t n = 16;   
                          __m512d phix = _mm512_loadu_pd(&pphix[0]);
                          __m512d pdf  = _mm512_loadu_pd(&ppdf[0]);
                          __m512d prod;
                          double sum;
                         int32_t err;
                         err = -1;
                         prod = _mm512_mul_pd(phix,pdf);
                         sum  = 0.0;
                         _mm512_storeu_pd(&intr[0],prod);
                         sum = cspint(&px[0],&intr[0],n,a,b,err);
                         ier = err;
                         if(ier==3) {
                             sum = std::numeric_limits<double>::quiet_NaN();
                             return (sum);
                         }      
                         return (sum);                                          
                }
                
                
                
                  /*
                      Helper formula serving the needs of 
                      formula 1.4, p. 15
                      The mathematical expectation of phase random 
                      variable phix(x).
                      Integrator 'cspint'.
                      Long data vector (arrays).
                 */
                 
                 
                  
                   double Ex_phix_zmm8r8_cspint_ne_a(const double * __restrict __ATTR_ALIGN__(64) pphix,
                                                   const double * __restrict __ATTR_ALIGN__(64) ppdf,
                                                   const double * __restrict __ATTR_ALIGN__(64) px,
                                                   double * __restrict __ATTR_ALIGN__(64)       intr,
                                                   CSPINT_DATA_1T csd,
                                                   const double a,
                                                   const double b,
                                                   const int32_t NTAB) {
                                                   
                        if(__buitlin_expect(NTAB==16,0)) {
                             double sum = 0.0;
                            sum = Ex_phix_zmm8r8_cspint_8e_a(pBx,ppdf,px,a,b);
                            return (sum);
                        }                
                        
                         __m512d phix,pdf,prod;
                         double sum;
                        int32_t i;
                        
                        for(i = 0; i != ROUND_TO_EIGHT(n,7); i += 8) {
                             _mm_prefetch((const char*)&pphix[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&ppdf[i],_MM_HINT_T0);
                             phix = _mm512_load_pd(&pphix[i]);
                             pdf  = _mm512_load_pd(&ppdf[i]);
                             prod = _mm512_mul_pd(phix,pdf);
                             _mm512_store_pd(&intr[i], prod); 
                        }  
                        sum = 0.0;
                        for(; i != n; ++i) {
                              double phix = pphix[i];
                              double pdf  = ppdf[i];
                              double prod = phix*pdf;
                             intr[i] = prod;
                        }  
                        
                        cspint(NTAB,&px[0],&intr[0],a,b,&csd.Y1[0],
                                &csd.Y2[0],&csd.Y3[0],&csd.E[0],&csd.WRK[0],sum);
                         return (sum);             
                  }
                 
                 
                 
                  
                   double Ex_phix_zmm8r8_cspint_ne_u(const double * __restrict pphix,
                                                   const double * __restrict  ppdf,
                                                   const double * __restrict  px,
                                                   double * __restrict intr,
                                                   CSPINT_DATA_1T csd,
                                                   const double a,
                                                   const double b,
                                                   const int32_t NTAB) {
                                                   
                        if(__buitlin_expect(NTAB==16,0)) {
                             double sum = 0.0;
                            sum = Ex_phix_zmm8r8_cspint_8e_u(pBx,ppdf,px,a,b);
                            return (sum);
                        }                
                        
                         __m512d phix,pdf,prod;
                         double sum;
                        int32_t i;
                        
                        for(i = 0; i != ROUND_TO_EIGHT(n,7); i += 8) {
                             _mm_prefetch((const char*)&pphix[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&ppdf[i],_MM_HINT_T0);
                             phix = _mm512_loadu_pd(&pphix[i]);
                             pdf  = _mm512_loadu_pd(&ppdf[i]);
                             prod = _mm512_mul_pd(phix,pdf);
                             _mm512_storeu_pd(&intr[i], prod); 
                        }  
                        sum = 0.0;
                        for(; i != n; ++i) {
                              double phix = pphix[i];
                              double pdf  = ppdf[i];
                              double prod = phix*pdf;
                             intr[i] = prod;
                        }  
                        
                        cspint(NTAB,&px[0],&intr[0],a,b,&csd.Y1[0],
                                &csd.Y2[0],&csd.Y3[0],&csd.E[0],&csd.WRK[0],sum);
                         return (sum);             
                  }
                  
                  
                   /*
                      Helper formula serving the needs of 
                      formula 1.4, p. 15
                      The mathematical expectation of phase random 
                      variable phix(x).
                      Integrator 'avint' irregular abscissas.
                      Long data vector (arrays).
                 */
                 
                 
                 
                   double Ex_phix_zmm8r8_avint_ne_a(const double * __restrict __ATTR_ALIGN__(64) pphix,
                                                   const double * __restrict __ATTR_ALIGN__(64) ppdf,
                                                   const double * __restrict __ATTR_ALIGN__(64) px,
                                                   double * __restrict __ATTR_ALIGN__(64)       intr,
                                                   const double a,
                                                   const double b,
                                                   const int32_t n,
                                                   int32_t ier) {
                                                   
                        if(__buitlin_expect(n==16,0)) {
                             double sum = 0.0;
                            sum = Ex_phix_zmm8r8_avint_8e_a(pBx,ppdf,px,a,b);
                            return (sum);
                        }                
                        
                         __m512d phix,pdf,prod;
                         double sum;
                        int32_t i,err;
                        
                        for(i = 0; i != ROUND_TO_EIGHT(n,7); i += 8) {
                             _mm_prefetch((const char*)&pphix[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&ppdf[i],_MM_HINT_T0);
                             phix = _mm512_load_pd(&pphix[i]);
                             pdf  = _mm512_load_pd(&ppdf[i]);
                             prod = _mm512_mul_pd(phix,pdf);
                             _mm512_store_pd(&intr[i], prod); 
                        }  
                        sum = 0.0;
                        for(; i != n; ++i) {
                              double phix = pphix[i];
                              double pdf  = ppdf[i];
                              double prod = phix*pdf;
                             intr[i] = prod;
                        }  
                        
                        sum = cspint(&px[0],&intr[0],n,a,b,err);
                        ier = err;
                        if(ier==3) {
                           sum = std::numeric_limits<double>::quiet_NaN();
                           return (sum);
                        }
                        return (sum);             
                  }
                  
                  
                  
                  
                   double Ex_phix_zmm8r8_avint_ne_u(const double * __restrict  pphix,
                                                   const double * __restrict  ppdf,
                                                   const double * __restrict px,
                                                   double * __restrict  intr,
                                                   const double a,
                                                   const double b,
                                                   const int32_t n,
                                                   int32_t ier) {
                                                   
                        if(__buitlin_expect(n==16,0)) {
                             double sum = 0.0;
                            sum = Ex_phix_zmm8r8_avint_8e_u(pBx,ppdf,px,a,b);
                            return (sum);
                        }                
                        
                         __m512d phix,pdf,prod;
                         double sum;
                        int32_t i,err;
                        
                        for(i = 0; i != ROUND_TO_EIGHT(n,7); i += 8) {
                             _mm_prefetch((const char*)&pphix[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&ppdf[i],_MM_HINT_T0);
                             phix = _mm512_loadu_pd(&pphix[i]);
                             pdf  = _mm512_loadu_pd(&ppdf[i]);
                             prod = _mm512_mul_pd(phix,pdf);
                             _mm512_storeu_pd(&intr[i], prod); 
                        }  
                        sum = 0.0;
                        for(; i != n; ++i) {
                              double phix = pphix[i];
                              double pdf  = ppdf[i];
                              double prod = phix*pdf;
                             intr[i] = prod;
                        }  
                        
                        sum = cspint(&px[0],&intr[0],n,a,b,err);
                        ier = err;
                        if(ier==3) {
                           sum = std::numeric_limits<double>::quiet_NaN();
                           return (sum);
                        }
                        return (sum);             
                  }
                  
                  
                      
                                 
                            
                  /*
                      Formula 1.5, p. 16
                      Far-zone Electric-field, n-terms
                  */
                  
                  
	           static inline
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
                                                  const int32_t n) {
                                                  
                         if(__builtin_expect(n<=0,0)) { return;}
                          __m512d fthr,fthi,Ethr,Ethi;
                          __m512d Ephr,Ephi,t0r,t0i;
                          __m512d t1r,t1i;   
                         int32_t j,m,m1;
                         
                         m = n%10;
                         if(m!=0) {
                             for(j = 0; j != m; ++j) {
                                  fthr = pfthr[j];
                                  fthi = pfthi[j];
                                  Ethr = pEthr[j];
                                  Ethi = pEthi[j];
                                  Ephr = pEphr[j];
                                  Ephi = pEphi[j];
                                  cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                                  Efztr[j] = t0r;
                                  Eftzi[j] = t0i;
                                  cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                                  Efzpr[j] = t1r;
                                  Efzpi[j] = t1i;
                             }
                             if(n<10) { return;}
                                   
                         }  
                         
                         m1 = m+1;
                         for(j = m1; j != n; j += 10) {
                              fthr = pfthr[j];
                              fthi = pfthi[j];
                              Ethr = pEthr[j];
                              Ethi = pEthi[j];
                              Ephr = pEphr[j];
                              Ephi = pEphi[j];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j] = t0r;
                              Eftzi[j] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j] = t1r;
                              Efzpi[j] = t1i;
                              fthr = pfthr[j+1];
                              fthi = pfthi[j+1];
                              Ethr = pEthr[j+1];
                              Ethi = pEthi[j+1];
                              Ephr = pEphr[j+1];
                              Ephi = pEphi[j+1];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+1] = t0r;
                              Eftzi[j+1] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+1] = t1r;
                              Efzpi[j+1] = t1i;
                              fthr = pfthr[j+2];
                              fthi = pfthi[j+2];
                              Ethr = pEthr[j+2];
                              Ethi = pEthi[j+2];
                              Ephr = pEphr[j+2];
                              Ephi = pEphi[j+2];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+2] = t0r;
                              Eftzi[j+2] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+2] = t1r;
                              Efzpi[j+2] = t1i;
                              fthr = pfthr[j+3];
                              fthi = pfthi[j+3];
                              Ethr = pEthr[j+3];
                              Ethi = pEthi[j+3];
                              Ephr = pEphr[j+3];
                              Ephi = pEphi[j+3];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+3] = t0r;
                              Eftzi[j+3] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+3] = t1r;
                              Efzpi[j+3] = t1i;
                              fthr = pfthr[j+4];
                              fthi = pfthi[j+4];
                              Ethr = pEthr[j+4];
                              Ethi = pEthi[j+4];
                              Ephr = pEphr[j+4];
                              Ephi = pEphi[j+4];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+4] = t0r;
                              Eftzi[j+4] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+4] = t1r;
                              Efzpi[j+4] = t1i;
                              fthr = pfthr[j+5];
                              fthi = pfthi[j+5];
                              Ethr = pEthr[j+5];
                              Ethi = pEthi[j+5];
                              Ephr = pEphr[j+5];
                              Ephi = pEphi[j+5];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+5] = t0r;
                              Eftzi[j+5] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+5] = t1r;
                              Efzpi[j+5] = t1i;
                              fthr = pfthr[j+6];
                              fthi = pfthi[j+6];
                              Ethr = pEthr[j+6];
                              Ethi = pEthi[j+6];
                              Ephr = pEphr[j+6];
                              Ephi = pEphi[j+6];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+6] = t0r;
                              Eftzi[j+6] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+6] = t1r;
                              Efzpi[j+6] = t1i;
                              fthr = pfthr[j+7];
                              fthi = pfthi[j+7];
                              Ethr = pEthr[j+7];
                              Ethi = pEthi[j+7];
                              Ephr = pEphr[j+7];
                              Ephi = pEphi[j+7];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+7] = t0r;
                              Eftzi[j+7] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+7] = t1r;
                              Efzpi[j+7] = t1i;
                              fthr = pfthr[j+8];
                              fthi = pfthi[j+8];
                              Ethr = pEthr[j+8];
                              Ethi = pEthi[j+8];
                              Ephr = pEphr[j+8];
                              Ephi = pEphi[j+8];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+8] = t0r;
                              Eftzi[j+8] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+8] = t1r;
                              Efzpi[j+8] = t1i;
                              fthr = pfthr[j+9];
                              fthi = pfthi[j+9];
                              Ethr = pEthr[j+9];
                              Ethi = pEthi[j+9];
                              Ephr = pEphr[j+9];
                              Ephi = pEphi[j+9];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+9] = t0r;
                              Eftzi[j+9] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+9] = t1r;
                              Efzpi[j+9] = t1i;
                            
                         }                           
                }
                 
                 
                 
                
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
                                                  const int32_t n) {
                                                  
                         if(__builtin_expect(n<=0,0)) { return;}
                          __m512d fthr,fthi,Ethr,Ethi;
                          __m512d Ephr,Ephi,t0r,t0i;
                          __m512d t1r,t1i;   
                         int32_t j,m,m1;
                         
                         m = n%8;
                         if(m!=0) {
                             for(j = 0; j != m; ++j) {
                                  fthr = pfthr[j];
                                  fthi = pfthi[j];
                                  Ethr = pEthr[j];
                                  Ethi = pEthi[j];
                                  Ephr = pEphr[j];
                                  Ephi = pEphi[j];
                                  cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                                  Efztr[j] = t0r;
                                  Eftzi[j] = t0i;
                                  cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                                  Efzpr[j] = t1r;
                                  Efzpi[j] = t1i;
                             }
                             if(n<8) { return;}
                                   
                         }  
                         
                         m1 = m+1;
                         for(j = m1; j != n; j += 8) {
                              fthr = pfthr[j];
                              fthi = pfthi[j];
                              Ethr = pEthr[j];
                              Ethi = pEthi[j];
                              Ephr = pEphr[j];
                              Ephi = pEphi[j];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j] = t0r;
                              Eftzi[j] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j] = t1r;
                              Efzpi[j] = t1i;
                              fthr = pfthr[j+1];
                              fthi = pfthi[j+1];
                              Ethr = pEthr[j+1];
                              Ethi = pEthi[j+1];
                              Ephr = pEphr[j+1];
                              Ephi = pEphi[j+1];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+1] = t0r;
                              Eftzi[j+1] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+1] = t1r;
                              Efzpi[j+1] = t1i;
                              fthr = pfthr[j+2];
                              fthi = pfthi[j+2];
                              Ethr = pEthr[j+2];
                              Ethi = pEthi[j+2];
                              Ephr = pEphr[j+2];
                              Ephi = pEphi[j+2];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+2] = t0r;
                              Eftzi[j+2] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+2] = t1r;
                              Efzpi[j+2] = t1i;
                              fthr = pfthr[j+3];
                              fthi = pfthi[j+3];
                              Ethr = pEthr[j+3];
                              Ethi = pEthi[j+3];
                              Ephr = pEphr[j+3];
                              Ephi = pEphi[j+3];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+3] = t0r;
                              Eftzi[j+3] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+3] = t1r;
                              Efzpi[j+3] = t1i;
                              fthr = pfthr[j+4];
                              fthi = pfthi[j+4];
                              Ethr = pEthr[j+4];
                              Ethi = pEthi[j+4];
                              Ephr = pEphr[j+4];
                              Ephi = pEphi[j+4];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+4] = t0r;
                              Eftzi[j+4] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+4] = t1r;
                              Efzpi[j+4] = t1i;
                              fthr = pfthr[j+5];
                              fthi = pfthi[j+5];
                              Ethr = pEthr[j+5];
                              Ethi = pEthi[j+5];
                              Ephr = pEphr[j+5];
                              Ephi = pEphi[j+5];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+5] = t0r;
                              Eftzi[j+5] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+5] = t1r;
                              Efzpi[j+5] = t1i;
                              fthr = pfthr[j+6];
                              fthi = pfthi[j+6];
                              Ethr = pEthr[j+6];
                              Ethi = pEthi[j+6];
                              Ephr = pEphr[j+6];
                              Ephi = pEphi[j+6];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+6] = t0r;
                              Eftzi[j+6] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+6] = t1r;
                              Efzpi[j+6] = t1i;
                              fthr = pfthr[j+7];
                              fthi = pfthi[j+7];
                              Ethr = pEthr[j+7];
                              Ethi = pEthi[j+7];
                              Ephr = pEphr[j+7];
                              Ephi = pEphi[j+7];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+7] = t0r;
                              Eftzi[j+7] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+7] = t1r;
                              Efzpi[j+7] = t1i;
                                                
                         }                           
                }
                 
                 
                 
                
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
                                                  const int32_t n) {
                                                  
                         if(__builtin_expect(n<=0,0)) { return;}
                          __m512d fthr,fthi,Ethr,Ethi;
                          __m512d Ephr,Ephi,t0r,t0i;
                          __m512d t1r,t1i;   
                         int32_t j,m,m1;
                         
                         m = n%4;
                         if(m!=0) {
                             for(j = 0; j != m; ++j) {
                                  fthr = pfthr[j];
                                  fthi = pfthi[j];
                                  Ethr = pEthr[j];
                                  Ethi = pEthi[j];
                                  Ephr = pEphr[j];
                                  Ephi = pEphi[j];
                                  cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                                  Efztr[j] = t0r;
                                  Eftzi[j] = t0i;
                                  cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                                  Efzpr[j] = t1r;
                                  Efzpi[j] = t1i;
                             }
                             if(n<4) { return;}
                                   
                         }  
                         
                         m1 = m+1;
                         for(j = m1; j != n; j += 4) {
                              fthr = pfthr[j];
                              fthi = pfthi[j];
                              Ethr = pEthr[j];
                              Ethi = pEthi[j];
                              Ephr = pEphr[j];
                              Ephi = pEphi[j];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j] = t0r;
                              Eftzi[j] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j] = t1r;
                              Efzpi[j] = t1i;
                              fthr = pfthr[j+1];
                              fthi = pfthi[j+1];
                              Ethr = pEthr[j+1];
                              Ethi = pEthi[j+1];
                              Ephr = pEphr[j+1];
                              Ephi = pEphi[j+1];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+1] = t0r;
                              Eftzi[j+1] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+1] = t1r;
                              Efzpi[j+1] = t1i;
                              fthr = pfthr[j+2];
                              fthi = pfthi[j+2];
                              Ethr = pEthr[j+2];
                              Ethi = pEthi[j+2];
                              Ephr = pEphr[j+2];
                              Ephi = pEphi[j+2];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+2] = t0r;
                              Eftzi[j+2] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+2] = t1r;
                              Efzpi[j+2] = t1i;
                              fthr = pfthr[j+3];
                              fthi = pfthi[j+3];
                              Ethr = pEthr[j+3];
                              Ethi = pEthi[j+3];
                              Ephr = pEphr[j+3];
                              Ephi = pEphi[j+3];
                              cmul_zmm8r8(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+3] = t0r;
                              Eftzi[j+3] = t0i;
                              cmul_zmm8r8(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+3] = t1r;
                              Efzpi[j+3] = t1i;
                                                                      
                         }                           
                }
                
                
              
                
                
                 
                 
      
  
