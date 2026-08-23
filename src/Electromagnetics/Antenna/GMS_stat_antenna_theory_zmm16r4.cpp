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
#include "GMS_stat_antenna_theory_zmm16r4.h"


               
               
                   /*
                        Formula 1.1, p. 14
                        Integrator cspint single data vector ZMM.
                   */
                 
                   std::complex<float> 
                   fth_f11_zmm16r4_cspint_16e_a(const float * __restrict __ATTR_ALIGN__(64) pAz,
                                            const float * __restrict __ATTR_ALIGN__(64) pphiz,
                                            const float * __restrict __ATTR_ALIGN__(64) pz,
                                            const float tht,
                                            const float gam,
                                            const int32_t L2) {
                                            
                         __ATTR_ALIGN__(64) float intr[16] = {};
                         __ATTR_ALIGN__(64) float inti[16] = {}; 
                         __ATTR_ALIGN__(64) float Y1[16]   = {};
                         __ATTR_ALIGN__(64) float Y2[16]   = {};
                         __ATTR_ALIGN__(64) float Y3[16]   = {}; 
                         __ATTR_ALIGN__(64) float E[16]    = {};
                         __ATTR_ALIGN__(64) float WRK[16]  = {}; 
                         constexpr int32_t NTAB = 16;               
                         constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                          __m512 Az    = _mm512_load_ps(&pAz[0]);
                          __m512 phiz  = _mm512_load_pd(&pphiz[0]);
                          __m512 z     = _mm512_load_pd(&pz[0]);
                          __m512 stht,ear,eai,cer,cei;
                          __m512 k,vtht;
                         std::complex<float> fth;
                          float sumr,sumi;
                         vtht= _mm512_set1_ps(tht);
                         k   = _mm512_mul_ps(C314159265358979323846264338328,
                                                          _mm512_set1_ps(gam));
                         stht= xsinf(tht);
                         ear = _mm512_setzero_ps();
                         eai = _mm512_fmadd_pd(stht,
                                           _mm512_mul_pd(k,z),phiz);
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         sumr= 0.0f;
                         cer = _mm512_mul_ps(Az,cer);
                         _mm512_store_ps(&intr[0],cer);
                         sumi= 0.0f;
                         cei = _mm512_mul_ps(Az,cei);
                         _mm512_store_ps(&inti[0],cei);
                         cspint(NTAB,&pz[0],&intr[0],-L2,L2,&Y1[0],&Y2[0],&Y3[0],&E[0],&WRK[0],sumr);
                         cspint(NTAB,&pz[0],&inti[0],-L2,L2,&Y1[0],&Y2[0],&Y3[0],&E[0],&WRK[0],sumi);
                         fth = {sumr,sumi};
                         return (fth);
                 }
                 
                 
                 
                   std::complex<float> 
                   fth_f11_zmm16r4_cspint_16e_u(const float * __restrict  pAz,
                                            const float * __restrict  pphiz,
                                            const float * __restrict  pz,
                                            const float tht,
                                            const float gam,
                                            const int32_t L2) {
                                            
                         __ATTR_ALIGN__(64) float intr[16] = {};
                         __ATTR_ALIGN__(64) float inti[16] = {}; 
                         __ATTR_ALIGN__(64) float Y1[16]   = {};
                         __ATTR_ALIGN__(64) float Y2[16]   = {};
                         __ATTR_ALIGN__(64) float Y3[16]   = {}; 
                         __ATTR_ALIGN__(64) float E[16]    = {};
                         __ATTR_ALIGN__(64) float WRK[16]  = {}; 
                         constexpr int32_t NTAB = 16;               
                         constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                          __m512 Az    = _mm512_loadu_ps(&pAz[0]);
                          __m512 phiz  = _mm512_loadu_pd(&pphiz[0]);
                          __m512 z     = _mm512_loadu_pd(&pz[0]);
                          __m512 stht,ear,eai,cer,cei;
                          __m512 k,vtht;
                         std::complex<float> fth;
                          float sumr,sumi;
                         vtht= _mm512_set1_ps(tht);
                         k   = _mm512_mul_ps(C314159265358979323846264338328,
                                                          _mm512_set1_ps(gam));
                         stht= xsinf(tht);
                         ear = _mm512_setzero_ps();
                         eai = _mm512_fmadd_pd(stht,
                                           _mm512_mul_pd(k,z),phiz);
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         sumr= 0.0f;
                         cer = _mm512_mul_ps(Az,cer);
                         _mm512_storeu_ps(&intr[0],cer);
                         sumi= 0.0f;
                         cei = _mm512_mul_ps(Az,cei);
                         _mm512_storeu_ps(&inti[0],cei);
                         cspint(NTAB,&pz[0],&intr[0],-L2,L2,&Y1[0],&Y2[0],&Y3[0],&E[0],&WRK[0],sumr);
                         cspint(NTAB,&pz[0],&inti[0],-L2,L2,&Y1[0],&Y2[0],&Y3[0],&E[0],&WRK[0],sumi);
                         fth = {sumr,sumi};
                         return (fth);
                 }
                 
                 
                   /*
                        Formula 1.1, p. 14
                        Integrator avint (irregular abscissas) single data vector ZMM.
                   */
                 
                 
                  
                   std::complex<float> 
                   fth_f11_zmm16r4_avint_16e_a(const float * __restrict __ATTR_ALIGN__(64) pAz,
                                               const float * __restrict __ATTR_ALIGN__(64) pphiz,
                                               const float * __restrict __ATTR_ALIGN__(64) pz,
                                               const float tht,
                                               const float gam,
                                               const int32_t L2,
                                               int32_t & ierr,
                                               int32_t & ieri) {
                                            
                         __ATTR_ALIGN__(64) float intr[16] = {};
                         __ATTR_ALIGN__(64) float inti[16] = {}; 
                         
                                      
                         constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                          __m512 Az    = _mm512_load_ps(&pAz[0]);
                          __m512 phiz  = _mm512_load_pd(&pphiz[0]);
                          __m512 z     = _mm512_load_pd(&pz[0]);
                          __m512 stht,ear,eai,cer,cei;
                          __m512 k,vtht;
                         std::complex<float> fth;
                          float sumr,sumi;
                         int32_t err,eri;
                         vtht= _mm512_set1_ps(tht);
                         k   = _mm512_mul_ps(C314159265358979323846264338328,
                                                          _mm512_set1_ps(gam));
                         stht= xsinf(tht);
                         ear = _mm512_setzero_ps();
                         eai = _mm512_fmadd_pd(stht,
                                           _mm512_mul_pd(k,z),phiz);
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         sumr= 0.0f;
                         cer = _mm512_mul_ps(Az,cer);
                         _mm512_store_ps(&intr[0],cer);
                         sumi= 0.0f;
                         cei = _mm512_mul_ps(Az,cei);
                         _mm512_store_ps(&inti[0],cei);
                         sumr = avint(&pz[0],&intr[0],16,-L2,L2,err);
                         sumi = avint(&pz[0],&inti[0],16,-L2,L2,eri);
                         ierr = err;
                         ieri = eri;
                         if(ierr==3 || ieri==3) {
                            const float NAN = std::numeric_limits<float>::quiet_NaN();
                            fth = {NAN,NAN};
                            return (fth);
                         }
                         fth = {sumr,sumi};
                         return (fth);
                 }
                 
                 
                   
                   std::complex<float> 
                   fth_f11_zmm16r4_avint_16e_u(const float * __restrict  pAz,
                                               const float * __restrict  pphiz,
                                               const float * __restrict  pz,
                                               const float tht,
                                               const float gam,
                                               const int32_t L2,
                                               int32_t & ierr,
                                               int32_t & ieri) {
                                            
                          float intr[16] = {};
                          float inti[16] = {}; 
                         
                                    
                         constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                          __m512 Az    = _mm512_loadu_ps(&pAz[0]);
                          __m512 phiz  = _mm512_loadu_pd(&pphiz[0]);
                          __m512 z     = _mm512_loadu_pd(&pz[0]);
                          __m512 stht,ear,eai,cer,cei;
                          __m512 k,vtht;
                         std::complex<float> fth;
                          float sumr,sumi;
                         int32_t err,eri;
                         vtht= _mm512_set1_ps(tht);
                         k   = _mm512_mul_ps(C314159265358979323846264338328,
                                                          _mm512_set1_ps(gam));
                         stht= xsinf(tht);
                         ear = _mm512_setzero_ps();
                         eai = _mm512_fmadd_pd(stht,
                                           _mm512_mul_pd(k,z),phiz);
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         sumr= 0.0f;
                         cer = _mm512_mul_ps(Az,cer);
                         _mm512_storeu_ps(&intr[0],cer);
                         sumi= 0.0f;
                         cei = _mm512_mul_ps(Az,cei);
                         _mm512_storeu_ps(&inti[0],cei);
                         sumr = avint(&pz[0],&intr[0],16,-L2,L2,err);
                         sumi = avint(&pz[0],&inti[0],16,-L2,L2,eri);
                         ierr = err;
                         ieri = eri;
                         if(ierr==3 || ieri==3) {
                            const float NAN = std::numeric_limits<float>::quiet_NaN();
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
                                            const int32_t NTAB) {
                                            
                        if(__builtin_expect(NTAB==16,0)) {
                           std::complex<float> fth = {0.0f,0.0f};
                           fth = f11_zmm16r4_cspint_16e_a(pAz,pphiz,pz,tht,gam,L2);
                           return (fth);
                        }     
                        
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         __m512 vtht,stht,vk;
                         __m512 ear,cer,cei;
                         std::complex<float> fth;
                         float sumr,sumi; 
                         float sint,k;
                        int32_t i; 
                        vtht = _mm512_set1_ps(tht);
                        sint = cephes_sinf(tht);
                        vk   = _mm512_mul_ps(C314159265358979323846264338328,
                                                          _mm512_set1_ps(gam));
                        ear  = _mm512_setzero_ps();
                        k    = C314159265358979323846264338328*gam;
                        stht = xsinf(vtht);
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pAz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pphiz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pz[i],_MM_HINT_T0);
                              __m512 t0 = _mm512_load_ps(&pz[i]);
                              __m512 t1 = _mm512_load_ps(&pphiz[i]);
                              __m512 eai= _mm512_fmadd_ps(stht,
                                                            _mm512_mul_ps(vk,t0),t1);
                              __m512 t2 = _mm512_load_ps(&pAz[i]);
                             cexp_zmm16r4(ear,eai,&cer,&cei);
                             cer = _mm512_mul_ps(cer,t2);
                             _mm512_store_ps(&intr[i],cer);
                             cei = _mm512_mul_ps(cei,t2);
                             _mm512_store_ps(&inti[i],cei);
                        }   
                        sumr = 0.0f;
                        sumi = 0.0f;
                        for(; i != NTAB; ++i) {
                              float t0 = pAz[i];
                              float t1 = pphiz[i];
                              float t2 = pz[i];
                              float zz = k*t0;
                              float eai= sint*zz+t1;
                              std::complex<float> c = std::exp({0.0f,eai});
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
                                            const int32_t NTAB) {
                                            
                        if(__builtin_expect(NTAB==16,0)) {
                           std::complex<float> fth = {0.0f,0.0f};
                           fth = f11_zmm16r4_cspint_16e_u(pAz,pphiz,pz,tht,gam,L2);
                           return (fth);
                        }     
                        
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         __m512 vtht,stht,vk;
                         __m512 ear,cer,cei;
                         std::complex<float> fth;
                         float sumr,sumi; 
                         float sint,k;
                        int32_t i; 
                        vtht = _mm512_set1_ps(tht);
                        sint = cephes_sinf(tht);
                        vk   = _mm512_mul_ps(C314159265358979323846264338328,
                                                          _mm512_set1_ps(gam));
                        ear  = _mm512_setzero_ps();
                        k    = C314159265358979323846264338328*gam;
                        stht = xsinf(vtht);
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pAz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pphiz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pz[i],_MM_HINT_T0);
                              __m512 t0 = _mm512_loadu_ps(&pz[i]);
                              __m512 t1 = _mm512_loadu_ps(&pphiz[i]);
                              __m512 eai= _mm512_fmadd_ps(stht,
                                                            _mm512_mul_ps(vk,t0),t1);
                              __m512 t2 = _mm512_loadu_ps(&pAz[i]);
                             cexp_zmm16r4(ear,eai,&cer,&cei);
                             cer = _mm512_mul_ps(cer,t2);
                             _mm512_storeu_ps(&intr[i],cer);
                             cei = _mm512_mul_ps(cei,t2);
                             _mm512_storeu_ps(&inti[i],cei);
                        }   
                        sumr = 0.0f;
                        sumi = 0.0f;
                        for(; i != NTAB; ++i) {
                              float t0 = pAz[i];
                              float t1 = pphiz[i];
                              float t2 = pz[i];
                              float zz = k*t0;
                              float eai= sint*zz+t1;
                              std::complex<float> c = std::exp({0.0f,eai});
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
                                            int32_t & ieri) {
                                            
                        if(__builtin_expect(NTAB==16,0)) {
                           std::complex<float> fth = {0.0f,0.0f};
                           fth = f11_zmm16r4_avint_16e_a(pAz,pphiz,pz,tht,gam,L2,ierr,ieri);
                           return (fth);
                        }     
                        
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         __m512 vtht,stht,vk;
                         __m512 ear,cer,cei;
                         std::complex<float> fth;
                         float sumr,sumi; 
                         float sint,k;
                        int32_t i,err,eri; 
                        vtht = _mm512_set1_ps(tht);
                        sint = cephes_sinf(tht);
                        vk   = _mm512_mul_ps(C314159265358979323846264338328,
                                                          _mm512_set1_ps(gam));
                        ear  = _mm512_setzero_ps();
                        k    = C314159265358979323846264338328*gam;
                        stht = xsinf(vtht);
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pAz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pphiz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pz[i],_MM_HINT_T0);
                              __m512 t0 = _mm512_load_ps(&pz[i]);
                              __m512 t1 = _mm512_load_ps(&pphiz[i]);
                              __m512 eai= _mm512_fmadd_ps(stht,
                                                            _mm512_mul_ps(vk,t0),t1);
                              __m512 t2 = _mm512_load_ps(&pAz[i]);
                             cexp_zmm16r4(ear,eai,&cer,&cei);
                             cer = _mm512_mul_ps(cer,t2);
                             _mm512_store_ps(&intr[i],cer);
                             cei = _mm512_mul_ps(cei,t2);
                             _mm512_store_ps(&inti[i],cei);
                        }   
                        sumr = 0.0f;
                        sumi = 0.0f;
                        for(; i != NTAB; ++i) {
                              float t0 = pAz[i];
                              float t1 = pphiz[i];
                              float t2 = pz[i];
                              float zz = k*t0;
                              float eai= sint*zz+t1;
                              std::complex<float> c = std::exp({0.0f,eai});
                             intr[i] = c.real()*t0;
                             inti[i] = c.imag()*t0;
                        } 
                        sumr = avint(&pz[0],&intr[0],n,-L2,L2,err); 
                        sumi = avint(&pz[0],&inti[0],n,-L2,L2,eri);
                        ierr = err;
                        ieri = eri;
                        if(ierr==3 || ieri==3) {
                            const float NAN = std::numeric_limits<float>::quiet_NaN();
                            fth = {NAN,NAN};
                            return (fth);
                         }
                        fth = {sumr,sumi};
                        return (fth);            
                }
                
                
                 
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
                                            int32_t & ieri) {
                                            
                        if(__builtin_expect(NTAB==16,0)) {
                           std::complex<float> fth = {0.0f,0.0f};
                           fth = f11_zmm16r4_avint_16e_u(pAz,pphiz,pz,tht,gam,L2,ierr,ieri);
                           return (fth);
                        }     
                        
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         __m512 vtht,stht,vk;
                         __m512 ear,cer,cei;
                         std::complex<float> fth;
                         float sumr,sumi; 
                         float sint,k;
                        int32_t i,err,eri; 
                        vtht = _mm512_set1_ps(tht);
                        sint = cephes_sinf(tht);
                        vk   = _mm512_mul_ps(C314159265358979323846264338328,
                                                          _mm512_set1_ps(gam));
                        ear  = _mm512_setzero_ps();
                        k    = C314159265358979323846264338328*gam;
                        stht = xsinf(vtht);
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pAz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pphiz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pz[i],_MM_HINT_T0);
                              __m512 t0 = _mm512_loadu_ps(&pz[i]);
                              __m512 t1 = _mm512_loadu_ps(&pphiz[i]);
                              __m512 eai= _mm512_fmadd_ps(stht,
                                                            _mm512_mul_ps(vk,t0),t1);
                              __m512 t2 = _mm512_load_ps(&pAz[i]);
                             cexp_zmm16r4(ear,eai,&cer,&cei);
                             cer = _mm512_mul_ps(cer,t2);
                             _mm512_storeu_ps(&intr[i],cer);
                             cei = _mm512_mul_ps(cei,t2);
                             _mm512_storeu_ps(&inti[i],cei);
                        }   
                        sumr = 0.0f;
                        sumi = 0.0f;
                        for(; i != NTAB; ++i) {
                              float t0 = pAz[i];
                              float t1 = pphiz[i];
                              float t2 = pz[i];
                              float zz = k*t0;
                              float eai= sint*zz+t1;
                              std::complex<float> c = std::exp({0.0f,eai});
                             intr[i] = c.real()*t0;
                             inti[i] = c.imag()*t0;
                        } 
                        sumr = avint(&pz[0],&intr[0],n,-L2,L2,err); 
                        sumi = avint(&pz[0],&inti[0],n,-L2,L2,eri);
                        ierr = err;
                        ieri = eri;
                        if(ierr==3 || ieri==3) {
                            const float NAN = std::numeric_limits<float>::quiet_NaN();
                            fth = {NAN,NAN};
                            return (fth);
                         }
                        fth = {sumr,sumi};
                        return (fth);            
                }
                
                
               
                          
                
                   void Bx_f13_zmm16r4_u10x_a(const float * __restrict __ATTR_ALIGN__(64) pA0,
                                              const float * __restrict __ATTR_ALIGN__(64) pA,
                                              float * __restrict __ATTR_ALIGN__(64) pB,
                                              const int32_t n) {
                                              
                        if(__builtin_expect(0==n,0)) { return;}
                         __m512 zmm0,zmm1,zmm2,zmm3;
                         __m512 zmm4,zmm5,zmm6,zmm7;
                         __m512 zmm8,zmm9,zmm10,zmm11;
                         __m512 zmm12,zmm13,zmm14,zmm15;
                        int32_t i;
                        
                        for(i = 0; (i+159) < n; i += 160) {
                             _mm_prefetch((const char*)&pA0[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+32], _MM_HINT_T0);
                             zmm0 = _mm512_load_ps(&pA0[i+0]);
                             zmm1 = _mm512_load_ps(&pA[i+0]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_store_ps(&pB[i+0],zmm3);
                             zmm4 = _mm512_load_ps(&pA0[i+16]);
                             zmm5 = _mm512_load_ps(&pA[i+16]);
                             zmm6 = _mm512_div_ps(zmm5,zmm4);
                             zmm7 = xlogf(zmm6);
                             _mm512_store_ps(&pB[i+16],zmm7);
                             _mm_prefetch((const char*)&pA0[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+64], _MM_HINT_T0);
                             zmm8 = _mm512_load_ps(&pA0[i+32]);
                             zmm9 = _mm512_load_ps(&pA[i+32]);
                             zmm10 = _mm512_div_ps(zmm9,zmm8);
                             zmm11 = xlogf(zmm10);
                             _mm512_store_ps(&pB[i+32],zmm10);
                             zmm12 = _mm512_load_ps(&pA0[i+48]);
                             zmm13 = _mm512_load_ps(&pA[i+48]);
                             zmm14 = _mm512_div_ps(zmm13,zmm12);
                             zmm15 = xlogf(zmm14);
                             _mm512_store_ps(&pB[i+48],zmm15);
                             _mm_prefetch((const char*)&pA0[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+96], _MM_HINT_T0);
                             zmm0 = _mm512_load_ps(&pA0[i+64]);
                             zmm1 = _mm512_load_ps(&pA[i+64]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_store_ps(&pB[i+64],zmm3);
                             zmm4 = _mm512_load_ps(&pA0[i+80]);
                             zmm5 = _mm512_load_ps(&pA[i+80]);
                             zmm6 = _mm512_div_ps(zmm5,zmm4);
                             zmm7 = xlogf(zmm6);
                             _mm512_store_ps(&pB[i+80],zmm7);
                             _mm_prefetch((const char*)&pA0[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+128], _MM_HINT_T0);
                             zmm8 = _mm512_load_ps(&pA0[i+96]);
                             zmm9 = _mm512_load_ps(&pA[i+96]);
                             zmm10= _mm512_div_ps(zmm9,zmm8);
                             zmm11= xlogf(zmm10);
                             _mm512_store_ps(&pB[i+96],zmm11);
                             zmm12 = _mm512_load_ps(&pA0[i+112]);
                             zmm13 = _mm512_load_ps(&pA[i+112]);
                             zmm14= _mm512_div_ps(zmm13,zmm12);
                             zmm15= xlogf(zmm14);
                             _mm512_store_ps(&pB[i+112],zmm15);
                             _mm_prefetch((const char*)&pA0[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+160], _MM_HINT_T0);
                             zmm0 = _mm512_load_ps(&pA0[i+128]);
                             zmm1 = _mm512_load_ps(&pA[i+128]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_store_ps(&pB[i+128],zmm3);
                             zmm4 = _mm512_load_ps(&pA0[i+144]);
                             zmm5 = _mm512_load_ps(&pA[i+144]);
                             zmm6 = _mm512_div_ps(zmm5,zmm4);
                             zmm7 = xlogf(zmm6);
                             _mm512_store_ps(&pB[i+144],zmm7);
                        }   
                        
                        for(; (i+127) < n; i += 128) {
                             zmm0 = _mm512_load_ps(&pA0[i+0]);
                             zmm1 = _mm512_load_ps(&pA[i+0]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_store_ps(&pB[i+0],zmm3);
                             zmm4 = _mm512_load_ps(&pA0[i+16]);
                             zmm5 = _mm512_load_ps(&pA[i+16]);
                             zmm6 = _mm512_div_ps(zmm5,zmm4);
                             zmm7 = xlogf(zmm6);
                             _mm512_store_ps(&pB[i+16],zmm7);
                             zmm8 = _mm512_load_ps(&pA0[i+32]);
                             zmm9 = _mm512_load_ps(&pA[i+32]);
                             zmm10 = _mm512_div_ps(zmm9,zmm8);
                             zmm11 = xlogf(zmm10);
                             _mm512_store_ps(&pB[i+32],zmm10);
                             zmm12 = _mm512_load_ps(&pA0[i+48]);
                             zmm13 = _mm512_load_ps(&pA[i+48]);
                             zmm14 = _mm512_div_ps(zmm13,zmm12);
                             zmm15 = xlogf(zmm14);
                             _mm512_store_ps(&pB[i+48],zmm15);
                             zmm0 = _mm512_load_ps(&pA0[i+64]);
                             zmm1 = _mm512_load_ps(&pA[i+64]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_store_ps(&pB[i+64],zmm3);
                             zmm4 = _mm512_load_ps(&pA0[i+80]);
                             zmm5 = _mm512_load_ps(&pA[i+80]);
                             zmm6 = _mm512_div_ps(zmm5,zmm4);
                             zmm7 = xlogf(zmm6);
                             _mm512_store_ps(&pB[i+80],zmm7);
                             zmm8 = _mm512_load_ps(&pA0[i+96]);
                             zmm9 = _mm512_load_ps(&pA[i+96]);
                             zmm10= _mm512_div_ps(zmm9,zmm8);
                             zmm11= xlogf(zmm10);
                             _mm512_store_ps(&pB[i+96],zmm11);
                             zmm12 = _mm512_load_ps(&pA0[i+112]);
                             zmm13 = _mm512_load_ps(&pA[i+112]);
                             zmm14= _mm512_div_ps(zmm13,zmm12);
                             zmm15= xlogf(zmm14);
                             _mm512_store_ps(&pB[i+112],zmm15);
                        }  
                        
                        for(; (i+79) < n; i += 80) {
                             zmm0 = _mm512_load_ps(&pA0[i+0]);
                             zmm1 = _mm512_load_ps(&pA[i+0]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_store_ps(&pB[i+0],zmm3);
                             zmm4 = _mm512_load_ps(&pA0[i+16]);
                             zmm5 = _mm512_load_ps(&pA[i+16]);
                             zmm6 = _mm512_div_ps(zmm5,zmm4);
                             zmm7 = xlogf(zmm6);
                             _mm512_store_ps(&pB[i+16],zmm7);
                             zmm8 = _mm512_load_ps(&pA0[i+32]);
                             zmm9 = _mm512_load_ps(&pA[i+32]);
                             zmm10 = _mm512_div_ps(zmm9,zmm8);
                             zmm11 = xlogf(zmm10);
                             _mm512_store_ps(&pB[i+32],zmm10);
                             zmm12 = _mm512_load_ps(&pA0[i+48]);
                             zmm13 = _mm512_load_ps(&pA[i+48]);
                             zmm14 = _mm512_div_ps(zmm13,zmm12);
                             zmm15 = xlogf(zmm14);
                             _mm512_store_ps(&pB[i+48],zmm15);
                             zmm0 = _mm512_load_ps(&pA0[i+64]);
                             zmm1 = _mm512_load_ps(&pA[i+64]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_store_ps(&pB[i+64],zmm3);
                        }   
                        
                        for(; (i+63) < n; i += 64) {
                             zmm0 = _mm512_load_ps(&pA0[i+0]);
                             zmm1 = _mm512_load_ps(&pA[i+0]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_store_ps(&pB[i+0],zmm3);
                             zmm4 = _mm512_load_ps(&pA0[i+16]);
                             zmm5 = _mm512_load_ps(&pA[i+16]);
                             zmm6 = _mm512_div_ps(zmm5,zmm4);
                             zmm7 = xlogf(zmm6);
                             _mm512_store_ps(&pB[i+16],zmm7);
                             zmm8 = _mm512_load_ps(&pA0[i+32]);
                             zmm9 = _mm512_load_ps(&pA[i+32]);
                             zmm10 = _mm512_div_ps(zmm9,zmm8);
                             zmm11 = xlogf(zmm10);
                             _mm512_store_ps(&pB[i+32],zmm10);
                             zmm12 = _mm512_load_ps(&pA0[i+48]);
                             zmm13 = _mm512_load_ps(&pA[i+48]);
                             zmm14 = _mm512_div_ps(zmm13,zmm12);
                             zmm15 = xlogf(zmm14);
                             _mm512_store_ps(&pB[i+48],zmm15);
                        } 
                        
                        for(; (i+31) < n; i += 32) {
                             zmm0 = _mm512_load_ps(&pA0[i+0]);
                             zmm1 = _mm512_load_ps(&pA[i+0]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_store_ps(&pB[i+0],zmm3);
                             zmm4 = _mm512_load_ps(&pA0[i+16]);
                             zmm5 = _mm512_load_ps(&pA[i+16]);
                             zmm6 = _mm512_div_ps(zmm5,zmm4);
                             zmm7 = xlogf(zmm6);
                             _mm512_store_ps(&pB[i+16],zmm7); 
                        } 
                        
                        for(; (i+15) < n; i += 16) {
                             zmm0 = _mm512_load_ps(&pA0[i+0]);
                             zmm1 = _mm512_load_ps(&pA[i+0]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_store_ps(&pB[i+0],zmm3);
                        }      
                        
                        for(; (i+0) < n; i += 1) {
                              float t0 = pA0[i];
                              float t1 = pA[i];
                              float rat= t1/t0;
                              float bx = ceph_logf(rat);
                             pB[i] = bx;
                        }           
                }
                
                
                
                  
                   void Bx_f13_zmm16r4_u10x_u(const float * __restrict pA0,
                                              const float * __restrict  pA,
                                              float * __restrict  pB,
                                              const int32_t n) {
                                              
                        if(__builtin_expect(0==n,0)) { return;}
                         __m512 zmm0,zmm1,zmm2,zmm3;
                         __m512 zmm4,zmm5,zmm6,zmm7;
                         __m512 zmm8,zmm9,zmm10,zmm11;
                         __m512 zmm12,zmm13,zmm14,zmm15;
                        int32_t i;
                        
                        for(i = 0; (i+159) < n; i += 160) {
                             _mm_prefetch((const char*)&pA0[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+32], _MM_HINT_T0);
                             zmm0 = _mm512_loadu_ps(&pA0[i+0]);
                             zmm1 = _mm512_loadu_ps(&pA[i+0]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_storeu_ps(&pB[i+0],zmm3);
                             zmm4 = _mm512_loadu_ps(&pA0[i+16]);
                             zmm5 = _mm512_loadu_ps(&pA[i+16]);
                             zmm6 = _mm512_div_ps(zmm5,zmm4);
                             zmm7 = xlogf(zmm6);
                             _mm512_storeu_ps(&pB[i+16],zmm7);
                             _mm_prefetch((const char*)&pA0[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+64], _MM_HINT_T0);
                             zmm8 = _mm512_loadu_ps(&pA0[i+32]);
                             zmm9 = _mm512_loadu_ps(&pA[i+32]);
                             zmm10 = _mm512_div_ps(zmm9,zmm8);
                             zmm11 = xlogf(zmm10);
                             _mm512_storeu_ps(&pB[i+32],zmm10);
                             zmm12 = _mm512_loadu_ps(&pA0[i+48]);
                             zmm13 = _mm512_loadu_ps(&pA[i+48]);
                             zmm14 = _mm512_div_ps(zmm13,zmm12);
                             zmm15 = xlogf(zmm14);
                             _mm512_storeu_ps(&pB[i+48],zmm15);
                             _mm_prefetch((const char*)&pA0[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+96], _MM_HINT_T0);
                             zmm0 = _mm512_loadu_ps(&pA0[i+64]);
                             zmm1 = _mm512_loadu_ps(&pA[i+64]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_storeu_ps(&pB[i+64],zmm3);
                             zmm4 = _mm512_loadu_ps(&pA0[i+80]);
                             zmm5 = _mm512_loadu_ps(&pA[i+80]);
                             zmm6 = _mm512_div_ps(zmm5,zmm4);
                             zmm7 = xlogf(zmm6);
                             _mm512_storeu_ps(&pB[i+80],zmm7);
                             _mm_prefetch((const char*)&pA0[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+128], _MM_HINT_T0);
                             zmm8 = _mm512_loadu_ps(&pA0[i+96]);
                             zmm9 = _mm512_loadu_ps(&pA[i+96]);
                             zmm10= _mm512_div_ps(zmm9,zmm8);
                             zmm11= xlogf(zmm10);
                             _mm512_storeu_ps(&pB[i+96],zmm11);
                             zmm12 = _mm512_loadu_ps(&pA0[i+112]);
                             zmm13 = _mm512_loadu_ps(&pA[i+112]);
                             zmm14= _mm512_div_ps(zmm13,zmm12);
                             zmm15= xlogf(zmm14);
                             _mm512_storeu_ps(&pB[i+112],zmm15);
                             _mm_prefetch((const char*)&pA0[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+160], _MM_HINT_T0);
                             zmm0 = _mm512_loadu_ps(&pA0[i+128]);
                             zmm1 = _mm512_loadu_ps(&pA[i+128]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_storeu_ps(&pB[i+128],zmm3);
                             zmm4 = _mm512_loadu_ps(&pA0[i+144]);
                             zmm5 = _mm512_loadu_ps(&pA[i+144]);
                             zmm6 = _mm512_div_ps(zmm5,zmm4);
                             zmm7 = xlogf(zmm6);
                             _mm512_storeu_ps(&pB[i+144],zmm7);
                        }   
                        
                        for(; (i+127) < n; i += 128) {
                             zmm0 = _mm512_loadu_ps(&pA0[i+0]);
                             zmm1 = _mm512_loadu_ps(&pA[i+0]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_storeu_ps(&pB[i+0],zmm3);
                             zmm4 = _mm512_loadu_ps(&pA0[i+16]);
                             zmm5 = _mm512_loadu_ps(&pA[i+16]);
                             zmm6 = _mm512_div_ps(zmm5,zmm4);
                             zmm7 = xlogf(zmm6);
                             _mm512_storeu_ps(&pB[i+16],zmm7);
                             zmm8 = _mm512_loadu_ps(&pA0[i+32]);
                             zmm9 = _mm512_loadu_ps(&pA[i+32]);
                             zmm10 = _mm512_div_ps(zmm9,zmm8);
                             zmm11 = xlogf(zmm10);
                             _mm512_storeu_ps(&pB[i+32],zmm10);
                             zmm12 = _mm512_loadu_ps(&pA0[i+48]);
                             zmm13 = _mm512_loadu_ps(&pA[i+48]);
                             zmm14 = _mm512_div_ps(zmm13,zmm12);
                             zmm15 = xlogf(zmm14);
                             _mm512_storeu_ps(&pB[i+48],zmm15);
                             zmm0 = _mm512_loadu_ps(&pA0[i+64]);
                             zmm1 = _mm512_loadu_ps(&pA[i+64]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_storeu_ps(&pB[i+64],zmm3);
                             zmm4 = _mm512_loadu_ps(&pA0[i+80]);
                             zmm5 = _mm512_loadu_ps(&pA[i+80]);
                             zmm6 = _mm512_div_ps(zmm5,zmm4);
                             zmm7 = xlogf(zmm6);
                             _mm512_storeu_ps(&pB[i+80],zmm7);
                             zmm8 = _mm512_loadu_ps(&pA0[i+96]);
                             zmm9 = _mm512_loadu_ps(&pA[i+96]);
                             zmm10= _mm512_div_ps(zmm9,zmm8);
                             zmm11= xlogf(zmm10);
                             _mm512_storeu_ps(&pB[i+96],zmm11);
                             zmm12 = _mm512_loadu_ps(&pA0[i+112]);
                             zmm13 = _mm512_loadu_ps(&pA[i+112]);
                             zmm14= _mm512_div_ps(zmm13,zmm12);
                             zmm15= xlogf(zmm14);
                             _mm512_storeu_ps(&pB[i+112],zmm15);
                        }  
                        
                        for(; (i+79) < n; i += 80) {
                             zmm0 = _mm512_loadu_ps(&pA0[i+0]);
                             zmm1 = _mm512_loadu_ps(&pA[i+0]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_storeu_ps(&pB[i+0],zmm3);
                             zmm4 = _mm512_loadu_ps(&pA0[i+16]);
                             zmm5 = _mm512_loadu_ps(&pA[i+16]);
                             zmm6 = _mm512_div_ps(zmm5,zmm4);
                             zmm7 = xlogf(zmm6);
                             _mm512_storeu_ps(&pB[i+16],zmm7);
                             zmm8 = _mm512_loadu_ps(&pA0[i+32]);
                             zmm9 = _mm512_loadu_ps(&pA[i+32]);
                             zmm10 = _mm512_div_ps(zmm9,zmm8);
                             zmm11 = xlogf(zmm10);
                             _mm512_store_ps(&pB[i+32],zmm10);
                             zmm12 = _mm512_loadu_ps(&pA0[i+48]);
                             zmm13 = _mm512_loadu_ps(&pA[i+48]);
                             zmm14 = _mm512_div_ps(zmm13,zmm12);
                             zmm15 = xlogf(zmm14);
                             _mm512_storeu_ps(&pB[i+48],zmm15);
                             zmm0 = _mm512_loadu_ps(&pA0[i+64]);
                             zmm1 = _mm512_loadu_ps(&pA[i+64]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_storeu_ps(&pB[i+64],zmm3);
                        }   
                        
                        for(; (i+63) < n; i += 64) {
                             zmm0 = _mm512_loadu_ps(&pA0[i+0]);
                             zmm1 = _mm512_loadu_ps(&pA[i+0]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_storeu_ps(&pB[i+0],zmm3);
                             zmm4 = _mm512_loadu_ps(&pA0[i+16]);
                             zmm5 = _mm512_loadu_ps(&pA[i+16]);
                             zmm6 = _mm512_div_ps(zmm5,zmm4);
                             zmm7 = xlogf(zmm6);
                             _mm512_storeu_ps(&pB[i+16],zmm7);
                             zmm8 = _mm512_loadu_ps(&pA0[i+32]);
                             zmm9 = _mm512_loadu_ps(&pA[i+32]);
                             zmm10 = _mm512_div_ps(zmm9,zmm8);
                             zmm11 = xlogf(zmm10);
                             _mm512_storeu_ps(&pB[i+32],zmm10);
                             zmm12 = _mm512_loadu_ps(&pA0[i+48]);
                             zmm13 = _mm512_loadu_ps(&pA[i+48]);
                             zmm14 = _mm512_div_ps(zmm13,zmm12);
                             zmm15 = xlogf(zmm14);
                             _mm512_storeu_ps(&pB[i+48],zmm15);
                        } 
                        
                        for(; (i+31) < n; i += 32) {
                             zmm0 = _mm512_loadu_ps(&pA0[i+0]);
                             zmm1 = _mm512_loadu_ps(&pA[i+0]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_storeu_ps(&pB[i+0],zmm3);
                             zmm4 = _mm512_loadu_ps(&pA0[i+16]);
                             zmm5 = _mm512_loadu_ps(&pA[i+16]);
                             zmm6 = _mm512_div_ps(zmm5,zmm4);
                             zmm7 = xlogf(zmm6);
                             _mm512_storeu_ps(&pB[i+16],zmm7); 
                        } 
                        
                        for(; (i+15) < n; i += 16) {
                             zmm0 = _mm512_loadu_ps(&pA0[i+0]);
                             zmm1 = _mm512_loadu_ps(&pA[i+0]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_storeu_ps(&pB[i+0],zmm3);
                        }      
                        
                        for(; (i+0) < n; i += 1) {
                              float t0 = pA0[i];
                              float t1 = pA[i];
                              float rat= t1/t0;
                              float bx = ceph_logf(rat);
                             pB[i] = bx;
                        }           
                }
                
                
                
               
                   void Bx_f13_zmm16r4_u2x_u(const float * __restrict pA0,
                                              const float * __restrict  pA,
                                              float * __restrict  pB,
                                              const int32_t n) {
                                              
                        if(__builtin_expect(0==n,0)) { return;}
                         __m512 zmm0,zmm1,zmm2,zmm3;
                         __m512 zmm4,zmm5,zmm6,zmm7;
                        int32_t i;
                        
                        for(i = 0; (i+31) < n; i += 32) {
                             _mm_prefetch((const char*)&pA0[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+32], _MM_HINT_T0);
                             zmm0 = _mm512_loadu_ps(&pA0[i+0]);
                             zmm1 = _mm512_loadu_ps(&pA[i+0]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_storeu_ps(&pB[i+0],zmm3);
                             zmm4 = _mm512_loadu_ps(&pA0[i+16]);
                             zmm5 = _mm512_loadu_ps(&pA[i+16]);
                             zmm6 = _mm512_div_ps(zmm5,zmm4);
                             zmm7 = xlogf(zmm6);
                             _mm512_storeu_ps(&pB[i+16],zmm7);
                       }  
                       
                        for(; (i+15) < n; i += 16) {
                             zmm0 = _mm512_loadu_ps(&pA0[i+0]);
                             zmm1 = _mm512_loadu_ps(&pA[i+0]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_storeu_ps(&pB[i+0],zmm3);
                        }      
                        
                        for(; (i+0) < n; i += 1) {
                              float t0 = pA0[i];
                              float t1 = pA[i];
                              float rat= t1/t0;
                              float bx = ceph_logf(rat);
                             pB[i] = bx;
                        }        
                        
                 }
                 
                 
                  
                   void Bx_f13_zmm16r4_u2x_a(const float * __restrict __ATTR_ALIGN__(64) pA0,
                                              const float * __restrict __ATTR_ALIGN__(64) pA,
                                              float * __restrict __ATTR_ALIGN__(64) pB,
                                              const int32_t n) {
                                              
                        if(__builtin_expect(0==n,0)) { return;}
                         __m512 zmm0,zmm1,zmm2,zmm3;
                         __m512 zmm4,zmm5,zmm6,zmm7;
                        int32_t i;
                        
                        for(i = 0; (i+31) < n; i += 32) {
                             _mm_prefetch((const char*)&pA0[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pA[i+32], _MM_HINT_T0);
                             zmm0 = _mm512_load_ps(&pA0[i+0]);
                             zmm1 = _mm512_load_ps(&pA[i+0]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_store_ps(&pB[i+0],zmm3);
                             zmm4 = _mm512_load_ps(&pA0[i+16]);
                             zmm5 = _mm512_load_ps(&pA[i+16]);
                             zmm6 = _mm512_div_ps(zmm5,zmm4);
                             zmm7 = xlogf(zmm6);
                             _mm512_store_ps(&pB[i+16],zmm7);
                       }  
                       
                        for(; (i+15) < n; i += 16) {
                             zmm0 = _mm512_load_ps(&pA0[i+0]);
                             zmm1 = _mm512_load_ps(&pA[i+0]);
                             zmm2 = _mm512_div_ps(zmm1,zmm0);
                             zmm3 = xlogf(zmm2);
                             _mm512_store_ps(&pB[i+0],zmm3);
                        }      
                        
                        for(; (i+0) < n; i += 1) {
                              float t0 = pA0[i];
                              float t1 = pA[i];
                              float rat= t1/t0;
                              float bx = ceph_logf(rat);
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
                 
                 
                   float Ex_Bx_zmm16r4_cspint_16e_a(const float * __restrict __ATTR_ALIGN__(64) pBx,
                                                    const float * __restrict __ATTR_ALIGN__(64) ppdf,
                                                    const float * __restrict __ATTR_ALIGN__(64) px,
                                                    const float a,
                                                    const float b) {
                                                    
                         __ATTR_ALIGN__(64) float intr[16] = {};
                         __ATTR_ALIGN__(64) float Y1[16]   = {};
                         __ATTR_ALIGN__(64) float Y2[16]   = {};
                         __ATTR_ALIGN__(64) float Y3[16]   = {}; 
                         __ATTR_ALIGN__(64) float E[16]    = {};
                         __ATTR_ALIGN__(64) float WRK[16]  = {}; 
                         constexpr int32_t NTAB = 16;   
                          __m512 Bx = _mm512_load_ps(&pBx[0]);
                          __m512 pdf= _mm512_load_ps(&ppdf[0]);
                          __m512 prod;
                         float sum;
                         prod = _mm512_mul_ps(Bx,pdf);
                         sum  = 0.0f;
                         _mm512_store_ps(&intr[0],prod);
                         cspint(NTAB,&px[0],&intr[0],a,b,&Y1[0],
                                &Y2[0],&Y3[0],&E[0],&WRK[0],sum);
                         return (sum);                                          
                }
                
                
                
                
                
                 
                   float Ex_Bx_zmm16r4_cspint_16e_u(const float * __restrict  pBx,
                                                    const float * __restrict  ppdf,
                                                    const float * __restrict  px,
                                                    const float a,
                                                    const float b) {
                                                    
                         __ATTR_ALIGN__(64) float intr[16] = {};
                         __ATTR_ALIGN__(64) float Y1[16]   = {};
                         __ATTR_ALIGN__(64) float Y2[16]   = {};
                         __ATTR_ALIGN__(64) float Y3[16]   = {}; 
                         __ATTR_ALIGN__(64) float E[16]    = {};
                         __ATTR_ALIGN__(64) float WRK[16]  = {}; 
                         constexpr int32_t NTAB = 16;   
                          __m512 Bx = _mm512_loadu_ps(&pBx[0]);
                          __m512 pdf= _mm512_loadu_ps(&ppdf[0]);
                          __m512 prod;
                         float sum;
                         prod = _mm512_mul_ps(Bx,pdf);
                         sum  = 0.0f;
                         _mm512_storeu_ps(&intr[0],prod);
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
                 
                 
                 
                   float Ex_Bx_zmm16r4_avint_16e_a(const float * __restrict __ATTR_ALIGN__(64) pBx,
                                                    const float * __restrict __ATTR_ALIGN__(64) ppdf,
                                                    const float * __restrict __ATTR_ALIGN__(64) px,
                                                    const float a,
                                                    const float b,
                                                    int32_t & ier){
                                                  
                                                    
                         __ATTR_ALIGN__(64) float intr[16] = {};
                                                 
                          __m512 Bx = _mm512_load_ps(&pBx[0]);
                          __m512 pdf= _mm512_load_ps(&ppdf[0]);
                          __m512 prod;
                         float sum;
                         int32_t err;
                         prod = _mm512_mul_ps(Bx,pdf);
                         sum  = 0.0f;
                         _mm512_store_ps(&intr[0],prod);
                         sum = avint(&px[0],&intr[0],16,a,b,err);
                         ier = err;
                         if(ier==3) {
                            sum = std::numeric_limits<float>::quiet_NaN();
                            return (sum);
                         }
                         return (sum);                                          
                }
                
                
                 
                   float Ex_Bx_zmm16r4_avint_16e_u(const float * __restrict  pBx,
                                                    const float * __restrict  ppdf,
                                                    const float * __restrict  px,
                                                    const float a,
                                                    const float b,
                                                    int32_t & ier){
                                                  
                                                    
                         __ATTR_ALIGN__(64) float intr[16] = {};
                                                 
                          __m512 Bx = _mm512_loadu_ps(&pBx[0]);
                          __m512 pdf= _mm512_loadu_ps(&ppdf[0]);
                          __m512 prod;
                         float sum;
                         int32_t err;
                         prod = _mm512_mul_ps(Bx,pdf);
                         sum  = 0.0f;
                         _mm512_storeu_ps(&intr[0],prod);
                         sum = avint(&px[0],&intr[0],16,a,b,err);
                         ier = err;
                         if(ier==3) {
                            sum = std::numeric_limits<float>::quiet_NaN();
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
                 
                 
                 
                 
                  
                   float Ex_Bx_zmm16r4_cspint_ne_a(const float * __restrict __ATTR_ALIGN__(64) pBx,
                                                   const float * __restrict __ATTR_ALIGN__(64) ppdf,
                                                   const float * __restrict __ATTR_ALIGN__(64) px,
                                                   float * __restrict __ATTR_ALIGN__(64)       intr,
                                                   CSPINT_DATA_1T csd,
                                                   const float a,
                                                   const float b,
                                                   const int32_t NTAB) {
                                                   
                        if(__buitlin_expect(NTAB==16,0)) {
                             float sum = 0.0f;
                            sum = Ex_Bx_zmm16r4_cspint_16e_a(pBx,ppdf,px,a,b);
                            return (sum);
                        }                
                        
                         __m512 Bx,pdf,prod;
                         float sum;
                        int32_t i;
                        
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pBx[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&ppdf[i],_MM_HINT_T0);
                             Bx = _mm512_load_ps(&pBx[i]);
                             pdf= _mm512_load_ps(&ppdf[i]);
                             prod = _mm512_mul_pd(Bx,pdf);
                             _mm512_store_ps(&intr[i], prod); 
                        }  
                        sum = 0.0f;
                        for(; i != NTAB; ++i) {
                              float Bx = pBx[i];
                              float pdf= ppdf[i];
                              float prod = Bx*pdf;
                             intr[i] = prod;
                        }  
                        
                        cspint(NTAB,&px[0],&intr[0],a,b,&csd.Y1[0],
                                &csd.Y2[0],&csd.Y3[0],&csd.E[0],&csd.WRK[0],sum);
                         return (sum);             
                  }
                  
                  
                  
                 
                   float Ex_Bx_zmm16r4_cspint_ne_u(const float * __restrict  pBx,
                                                   const float * __restrict  ppdf,
                                                   const float * __restrict  px,
                                                   float * __restrict intr,
                                                   CSPINT_DATA_1T csd,
                                                   const float a,
                                                   const float b,
                                                   const int32_t NTAB) {
                                                   
                        if(__buitlin_expect(NTAB==16,0)) {
                             float sum = 0.0f;
                            sum = Ex_Bx_zmm16r4_cspint_16e_u(pBx,ppdf,px,a,b);
                            return (sum);
                        }                
                        
                         __m512 Bx,pdf,prod;
                         float sum;
                        int32_t i;
                        
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pBx[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&ppdf[i],_MM_HINT_T0);
                             Bx = _mm512_loadu_ps(&pBx[i]);
                             pdf= _mm512_loadu_ps(&ppdf[i]);
                             prod = _mm512_mul_pd(Bx,pdf);
                             _mm512_storeu_ps(&intr[i], prod); 
                        }  
                        sum = 0.0f;
                        for(; i != NTAB; ++i) {
                              float Bx = pBx[i];
                              float pdf= ppdf[i];
                              float prod = Bx*pdf;
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
                 
                 
                 
                   float Ex_Bx_zmm16r4_avint_ne_a(const float * __restrict __ATTR_ALIGN__(64) pBx,
                                                   const float * __restrict __ATTR_ALIGN__(64) ppdf,
                                                   const float * __restrict __ATTR_ALIGN__(64) px,
                                                   float * __restrict __ATTR_ALIGN__(64)       intr,
                                                   const int32_t n,
                                                   const float a,
                                                   const float b,
                                                   int32_t & ier) {
                                                   
                        if(__buitlin_expect(n==16,0)) {
                             float sum = 0.0f;
                            sum = Ex_Bx_zmm16r4_avint_16e_a(pBx,ppdf,px,a,b,ier);
                            return (sum);
                        }                
                        
                         __m512 Bx,pdf,prod;
                         float sum;
                        int32_t i,err;
                        
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pBx[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&ppdf[i],_MM_HINT_T0);
                             Bx = _mm512_load_ps(&pBx[i]);
                             pdf= _mm512_load_ps(&ppdf[i]);
                             prod = _mm512_mul_pd(Bx,pdf);
                             _mm512_store_ps(&intr[i], prod); 
                        }  
                        sum = 0.0f;
                        for(; i != NTAB; ++i) {
                              float Bx = pBx[i];
                              float pdf= ppdf[i];
                              float prod = Bx*pdf;
                             intr[i] = prod;
                        }  
                        
                        sum = avint(&px[0],&intr[0],n,a,b,err);
                        ier = err;
                        if(ier==3) {
                           sum = std::numeric_limits<float>::quiet_NaN();
                           return (sum);
                        }
                         return (sum);             
                  }
                  
                  
                  
                 
                   float Ex_Bx_zmm16r4_avint_ne_u(const float * __restrict  pBx,
                                                   const float * __restrict  ppdf,
                                                   const float * __restrict  px,
                                                   float * __restrict intr,
                                                   const int32_t n,
                                                   const float a,
                                                   const float b,
                                                   int32_t & ier) {
                                                   
                        if(__buitlin_expect(n==16,0)) {
                             float sum = 0.0f;
                            sum = Ex_Bx_zmm16r4_avint_16e_u(pBx,ppdf,px,a,b,ier);
                            return (sum);
                        }                
                        
                         __m512 Bx,pdf,prod;
                         float sum;
                        int32_t i,err;
                        
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pBx[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&ppdf[i],_MM_HINT_T0);
                             Bx = _mm512_loadu_ps(&pBx[i]);
                             pdf= _mm512_loadu_ps(&ppdf[i]);
                             prod = _mm512_mul_pd(Bx,pdf);
                             _mm512_storeu_ps(&intr[i], prod); 
                        }  
                        sum = 0.0f;
                        for(; i != NTAB; ++i) {
                              float Bx = pBx[i];
                              float pdf= ppdf[i];
                              float prod = Bx*pdf;
                             intr[i] = prod;
                        }  
                        
                        sum = avint(&px[0],&intr[0],n,a,b,err);
                        ier = err;
                        if(ier==3) {
                           sum = std::numeric_limits<float>::quiet_NaN();
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
                 
                 
                  
                   float Ex_phix_zmm16r4_cspint_16e_a(const float * __restrict __ATTR_ALIGN__(64) pphix,
                                                      const float * __restrict __ATTR_ALIGN__(64) ppdf,
                                                      const float * __restrict __ATTR_ALIGN__(64) px,
                                                      const float a,
                                                      const float b) {
                                                    
                         __ATTR_ALIGN__(64) float intr[16] = {};
                         __ATTR_ALIGN__(64) float Y1[16]   = {};
                         __ATTR_ALIGN__(64) float Y2[16]   = {};
                         __ATTR_ALIGN__(64) float Y3[16]   = {}; 
                         __ATTR_ALIGN__(64) float E[16]    = {};
                         __ATTR_ALIGN__(64) float WRK[16]  = {}; 
                         constexpr int32_t NTAB = 16;   
                          __m512 phix = _mm512_load_ps(&pphix[0]);
                          __m512 pdf  = _mm512_load_ps(&ppdf[0]);
                          __m512 prod;
                         float sum;
                         prod = _mm512_mul_ps(phix,pdf);
                         sum  = 0.0f;
                         _mm512_store_ps(&intr[0],prod);
                         cspint(NTAB,&px[0],&intr[0],a,b,&Y1[0],
                                &Y2[0],&Y3[0],&E[0],&WRK[0],sum);
                         return (sum);                                          
                }
                
                
                
                  
                   float Ex_phix_zmm16r4_cspint_16e_u(const float * __restrict  pphix,
                                                      const float * __restrict  ppdf,
                                                      const float * __restrict  px,
                                                      const float a,
                                                      const float b) {
                                                    
                         __ATTR_ALIGN__(64) float intr[16] = {};
                         __ATTR_ALIGN__(64) float Y1[16]   = {};
                         __ATTR_ALIGN__(64) float Y2[16]   = {};
                         __ATTR_ALIGN__(64) float Y3[16]   = {}; 
                         __ATTR_ALIGN__(64) float E[16]    = {};
                         __ATTR_ALIGN__(64) float WRK[16]  = {}; 
                         constexpr int32_t NTAB = 16;   
                          __m512 phix = _mm512_loadu_ps(&pphix[0]);
                          __m512 pdf  = _mm512_loadu_ps(&ppdf[0]);
                          __m512 prod;
                         float sum;
                         prod = _mm512_mul_ps(phix,pdf);
                         sum  = 0.0f;
                         _mm512_storeu_ps(&intr[0],prod);
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
                 
                 
                  
                   float Ex_phix_zmm16r4_avint_16e_a(const float * __restrict __ATTR_ALIGN__(64) pphix,
                                                      const float * __restrict __ATTR_ALIGN__(64) ppdf,
                                                      const float * __restrict __ATTR_ALIGN__(64) px,
                                                      const float a,
                                                      const float b,
                                                      int32_t & ier) {
                                                    
                         __ATTR_ALIGN__(64) float intr[16] = {};
                         constexpr int32_t n = 16;   
                          __m512 phix = _mm512_load_ps(&pphix[0]);
                          __m512 pdf  = _mm512_load_ps(&ppdf[0]);
                          __m512 prod;
                          float sum;
                         int32_t err;
                         err = -1;
                         prod = _mm512_mul_ps(phix,pdf);
                         sum  = 0.0f;
                         _mm512_store_ps(&intr[0],prod);
                         sum = cspint(&px[0],&intr[0],n,a,b,err);
                         ier = err;
                         if(ier==3) {
                             sum = std::numeric_limits<float>::quiet_NaN();
                             return (sum);
                         }      
                         return (sum);                                          
                }
                
                
                 
                   float Ex_phix_zmm16r4_avint_16e_u(const float * __restrict pphix,
                                                      const float * __restrict ppdf,
                                                      const float * __restrict  px,
                                                      const float a,
                                                      const float b,
                                                      int32_t & ier) {
                                                    
                         float intr[16] = {};
                         constexpr int32_t n = 16;   
                          __m512 phix = _mm512_loadu_ps(&pphix[0]);
                          __m512 pdf  = _mm512_loadu_ps(&ppdf[0]);
                          __m512 prod;
                          float sum;
                         int32_t err;
                         err = -1;
                         prod = _mm512_mul_ps(phix,pdf);
                         sum  = 0.0f;
                         _mm512_storeu_ps(&intr[0],prod);
                         sum = cspint(&px[0],&intr[0],n,a,b,err);
                         ier = err;
                         if(ier==3) {
                             sum = std::numeric_limits<float>::quiet_NaN();
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
                 
                 
                  
                   float Ex_phix_zmm16r4_cspint_ne_a(const float * __restrict __ATTR_ALIGN__(64) pphix,
                                                   const float * __restrict __ATTR_ALIGN__(64) ppdf,
                                                   const float * __restrict __ATTR_ALIGN__(64) px,
                                                   float * __restrict __ATTR_ALIGN__(64)       intr,
                                                   CSPINT_DATA_1T csd,
                                                   const float a,
                                                   const float b,
                                                   const int32_t NTAB) {
                                                   
                        if(__buitlin_expect(NTAB==16,0)) {
                             float sum = 0.0f;
                            sum = Ex_phix_zmm16r4_cspint_16e_a(pBx,ppdf,px,a,b);
                            return (sum);
                        }                
                        
                         __m512 phix,pdf,prod;
                         float sum;
                        int32_t i;
                        
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pphix[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&ppdf[i],_MM_HINT_T0);
                             phix = _mm512_load_ps(&pphix[i]);
                             pdf  = _mm512_load_ps(&ppdf[i]);
                             prod = _mm512_mul_pd(phix,pdf);
                             _mm512_store_ps(&intr[i], prod); 
                        }  
                        sum = 0.0f;
                        for(; i != NTAB; ++i) {
                              float phix = pphix[i];
                              float pdf  = ppdf[i];
                              float prod = phix*pdf;
                             intr[i] = prod;
                        }  
                        
                        cspint(NTAB,&px[0],&intr[0],a,b,&csd.Y1[0],
                                &csd.Y2[0],&csd.Y3[0],&csd.E[0],&csd.WRK[0],sum);
                         return (sum);             
                  }
                 
                 
                 
                  
                   float Ex_phix_zmm16r4_cspint_ne_u(const float * __restrict pphix,
                                                   const float * __restrict  ppdf,
                                                   const float * __restrict  px,
                                                   float * __restrict intr,
                                                   CSPINT_DATA_1T csd,
                                                   const float a,
                                                   const float b,
                                                   const int32_t NTAB) {
                                                   
                        if(__buitlin_expect(NTAB==16,0)) {
                             float sum = 0.0f;
                            sum = Ex_phix_zmm16r4_cspint_16e_u(pBx,ppdf,px,a,b);
                            return (sum);
                        }                
                        
                         __m512 phix,pdf,prod;
                         float sum;
                        int32_t i;
                        
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pphix[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&ppdf[i],_MM_HINT_T0);
                             phix = _mm512_loadu_ps(&pphix[i]);
                             pdf  = _mm512_loadu_ps(&ppdf[i]);
                             prod = _mm512_mul_pd(phix,pdf);
                             _mm512_storeu_ps(&intr[i], prod); 
                        }  
                        sum = 0.0f;
                        for(; i != NTAB; ++i) {
                              float phix = pphix[i];
                              float pdf  = ppdf[i];
                              float prod = phix*pdf;
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
                 
                 
                 
                   float Ex_phix_zmm16r4_avint_ne_a(const float * __restrict __ATTR_ALIGN__(64) pphix,
                                                   const float * __restrict __ATTR_ALIGN__(64) ppdf,
                                                   const float * __restrict __ATTR_ALIGN__(64) px,
                                                   float * __restrict __ATTR_ALIGN__(64)       intr,
                                                   const float a,
                                                   const float b,
                                                   const int32_t n,
                                                   int32_t ier) {
                                                   
                        if(__buitlin_expect(n==16,0)) {
                             float sum = 0.0f;
                            sum = Ex_phix_zmm16r4_avint_16e_a(pBx,ppdf,px,a,b);
                            return (sum);
                        }                
                        
                         __m512 phix,pdf,prod;
                         float sum;
                        int32_t i,err;
                        
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pphix[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&ppdf[i],_MM_HINT_T0);
                             phix = _mm512_load_ps(&pphix[i]);
                             pdf  = _mm512_load_ps(&ppdf[i]);
                             prod = _mm512_mul_pd(phix,pdf);
                             _mm512_store_ps(&intr[i], prod); 
                        }  
                        sum = 0.0f;
                        for(; i != NTAB; ++i) {
                              float phix = pphix[i];
                              float pdf  = ppdf[i];
                              float prod = phix*pdf;
                             intr[i] = prod;
                        }  
                        
                        sum = cspint(&px[0],&intr[0],n,a,b,err);
                        ier = err;
                        if(ier==3) {
                           sum = std::numeric_limits<float>::quiet_NaN();
                           return (sum);
                        }
                        return (sum);             
                  }
                  
                  
                  
                  
                   float Ex_phix_zmm16r4_avint_ne_u(const float * __restrict  pphix,
                                                   const float * __restrict  ppdf,
                                                   const float * __restrict px,
                                                   float * __restrict  intr,
                                                   const float a,
                                                   const float b,
                                                   const int32_t n,
                                                   int32_t ier) {
                                                   
                        if(__buitlin_expect(n==16,0)) {
                             float sum = 0.0f;
                            sum = Ex_phix_zmm16r4_avint_16e_u(pBx,ppdf,px,a,b);
                            return (sum);
                        }                
                        
                         __m512 phix,pdf,prod;
                         float sum;
                        int32_t i,err;
                        
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pphix[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&ppdf[i],_MM_HINT_T0);
                             phix = _mm512_loadu_ps(&pphix[i]);
                             pdf  = _mm512_loadu_ps(&ppdf[i]);
                             prod = _mm512_mul_pd(phix,pdf);
                             _mm512_storeu_ps(&intr[i], prod); 
                        }  
                        sum = 0.0f;
                        for(; i != NTAB; ++i) {
                              float phix = pphix[i];
                              float pdf  = ppdf[i];
                              float prod = phix*pdf;
                             intr[i] = prod;
                        }  
                        
                        sum = cspint(&px[0],&intr[0],n,a,b,err);
                        ier = err;
                        if(ier==3) {
                           sum = std::numeric_limits<float>::quiet_NaN();
                           return (sum);
                        }
                        return (sum);             
                  }
                  
                  
                      
                                 
                            
                  /*
                      Formula 1.5, p. 16
                      Far-zone Electric-field, n-terms
                  */
                  
                  
	           static inline
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
                                                  const int32_t n) {
                                                  
                         if(__builtin_expect(n<=0,0)) { return;}
                          __m512 fthr,fthi,Ethr,Ethi;
                          __m512 Ephr,Ephi,t0r,t0i;
                          __m512 t1r,t1i;   
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
                                  cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                                  Efztr[j] = t0r;
                                  Eftzi[j] = t0i;
                                  cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
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
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j] = t0r;
                              Eftzi[j] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j] = t1r;
                              Efzpi[j] = t1i;
                              fthr = pfthr[j+1];
                              fthi = pfthi[j+1];
                              Ethr = pEthr[j+1];
                              Ethi = pEthi[j+1];
                              Ephr = pEphr[j+1];
                              Ephi = pEphi[j+1];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+1] = t0r;
                              Eftzi[j+1] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+1] = t1r;
                              Efzpi[j+1] = t1i;
                              fthr = pfthr[j+2];
                              fthi = pfthi[j+2];
                              Ethr = pEthr[j+2];
                              Ethi = pEthi[j+2];
                              Ephr = pEphr[j+2];
                              Ephi = pEphi[j+2];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+2] = t0r;
                              Eftzi[j+2] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+2] = t1r;
                              Efzpi[j+2] = t1i;
                              fthr = pfthr[j+3];
                              fthi = pfthi[j+3];
                              Ethr = pEthr[j+3];
                              Ethi = pEthi[j+3];
                              Ephr = pEphr[j+3];
                              Ephi = pEphi[j+3];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+3] = t0r;
                              Eftzi[j+3] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+3] = t1r;
                              Efzpi[j+3] = t1i;
                              fthr = pfthr[j+4];
                              fthi = pfthi[j+4];
                              Ethr = pEthr[j+4];
                              Ethi = pEthi[j+4];
                              Ephr = pEphr[j+4];
                              Ephi = pEphi[j+4];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+4] = t0r;
                              Eftzi[j+4] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+4] = t1r;
                              Efzpi[j+4] = t1i;
                              fthr = pfthr[j+5];
                              fthi = pfthi[j+5];
                              Ethr = pEthr[j+5];
                              Ethi = pEthi[j+5];
                              Ephr = pEphr[j+5];
                              Ephi = pEphi[j+5];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+5] = t0r;
                              Eftzi[j+5] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+5] = t1r;
                              Efzpi[j+5] = t1i;
                              fthr = pfthr[j+6];
                              fthi = pfthi[j+6];
                              Ethr = pEthr[j+6];
                              Ethi = pEthi[j+6];
                              Ephr = pEphr[j+6];
                              Ephi = pEphi[j+6];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+6] = t0r;
                              Eftzi[j+6] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+6] = t1r;
                              Efzpi[j+6] = t1i;
                              fthr = pfthr[j+7];
                              fthi = pfthi[j+7];
                              Ethr = pEthr[j+7];
                              Ethi = pEthi[j+7];
                              Ephr = pEphr[j+7];
                              Ephi = pEphi[j+7];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+7] = t0r;
                              Eftzi[j+7] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+7] = t1r;
                              Efzpi[j+7] = t1i;
                              fthr = pfthr[j+8];
                              fthi = pfthi[j+8];
                              Ethr = pEthr[j+8];
                              Ethi = pEthi[j+8];
                              Ephr = pEphr[j+8];
                              Ephi = pEphi[j+8];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+8] = t0r;
                              Eftzi[j+8] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+8] = t1r;
                              Efzpi[j+8] = t1i;
                              fthr = pfthr[j+9];
                              fthi = pfthi[j+9];
                              Ethr = pEthr[j+9];
                              Ethi = pEthi[j+9];
                              Ephr = pEphr[j+9];
                              Ephi = pEphi[j+9];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+9] = t0r;
                              Eftzi[j+9] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+9] = t1r;
                              Efzpi[j+9] = t1i;
                            
                         }                           
                }
                 
                 
                 
                
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
                                                  const int32_t n) {
                                                  
                         if(__builtin_expect(n<=0,0)) { return;}
                          __m512 fthr,fthi,Ethr,Ethi;
                          __m512 Ephr,Ephi,t0r,t0i;
                          __m512 t1r,t1i;   
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
                                  cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                                  Efztr[j] = t0r;
                                  Eftzi[j] = t0i;
                                  cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
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
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j] = t0r;
                              Eftzi[j] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j] = t1r;
                              Efzpi[j] = t1i;
                              fthr = pfthr[j+1];
                              fthi = pfthi[j+1];
                              Ethr = pEthr[j+1];
                              Ethi = pEthi[j+1];
                              Ephr = pEphr[j+1];
                              Ephi = pEphi[j+1];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+1] = t0r;
                              Eftzi[j+1] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+1] = t1r;
                              Efzpi[j+1] = t1i;
                              fthr = pfthr[j+2];
                              fthi = pfthi[j+2];
                              Ethr = pEthr[j+2];
                              Ethi = pEthi[j+2];
                              Ephr = pEphr[j+2];
                              Ephi = pEphi[j+2];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+2] = t0r;
                              Eftzi[j+2] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+2] = t1r;
                              Efzpi[j+2] = t1i;
                              fthr = pfthr[j+3];
                              fthi = pfthi[j+3];
                              Ethr = pEthr[j+3];
                              Ethi = pEthi[j+3];
                              Ephr = pEphr[j+3];
                              Ephi = pEphi[j+3];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+3] = t0r;
                              Eftzi[j+3] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+3] = t1r;
                              Efzpi[j+3] = t1i;
                              fthr = pfthr[j+4];
                              fthi = pfthi[j+4];
                              Ethr = pEthr[j+4];
                              Ethi = pEthi[j+4];
                              Ephr = pEphr[j+4];
                              Ephi = pEphi[j+4];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+4] = t0r;
                              Eftzi[j+4] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+4] = t1r;
                              Efzpi[j+4] = t1i;
                              fthr = pfthr[j+5];
                              fthi = pfthi[j+5];
                              Ethr = pEthr[j+5];
                              Ethi = pEthi[j+5];
                              Ephr = pEphr[j+5];
                              Ephi = pEphi[j+5];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+5] = t0r;
                              Eftzi[j+5] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+5] = t1r;
                              Efzpi[j+5] = t1i;
                              fthr = pfthr[j+6];
                              fthi = pfthi[j+6];
                              Ethr = pEthr[j+6];
                              Ethi = pEthi[j+6];
                              Ephr = pEphr[j+6];
                              Ephi = pEphi[j+6];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+6] = t0r;
                              Eftzi[j+6] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+6] = t1r;
                              Efzpi[j+6] = t1i;
                              fthr = pfthr[j+7];
                              fthi = pfthi[j+7];
                              Ethr = pEthr[j+7];
                              Ethi = pEthi[j+7];
                              Ephr = pEphr[j+7];
                              Ephi = pEphi[j+7];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+7] = t0r;
                              Eftzi[j+7] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+7] = t1r;
                              Efzpi[j+7] = t1i;
                                                
                         }                           
                }
                 
                 
                 
                
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
                                                  const int32_t n) {
                                                  
                         if(__builtin_expect(n<=0,0)) { return;}
                          __m512 fthr,fthi,Ethr,Ethi;
                          __m512 Ephr,Ephi,t0r,t0i;
                          __m512 t1r,t1i;   
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
                                  cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                                  Efztr[j] = t0r;
                                  Eftzi[j] = t0i;
                                  cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
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
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j] = t0r;
                              Eftzi[j] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j] = t1r;
                              Efzpi[j] = t1i;
                              fthr = pfthr[j+1];
                              fthi = pfthi[j+1];
                              Ethr = pEthr[j+1];
                              Ethi = pEthi[j+1];
                              Ephr = pEphr[j+1];
                              Ephi = pEphi[j+1];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+1] = t0r;
                              Eftzi[j+1] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+1] = t1r;
                              Efzpi[j+1] = t1i;
                              fthr = pfthr[j+2];
                              fthi = pfthi[j+2];
                              Ethr = pEthr[j+2];
                              Ethi = pEthi[j+2];
                              Ephr = pEphr[j+2];
                              Ephi = pEphi[j+2];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+2] = t0r;
                              Eftzi[j+2] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+2] = t1r;
                              Efzpi[j+2] = t1i;
                              fthr = pfthr[j+3];
                              fthi = pfthi[j+3];
                              Ethr = pEthr[j+3];
                              Ethi = pEthi[j+3];
                              Ephr = pEphr[j+3];
                              Ephi = pEphi[j+3];
                              cmul_zmm16r4(fthr,fthi,Ethr,Ethi,&t0r,&t0i);
                              Efztr[j+3] = t0r;
                              Eftzi[j+3] = t0i;
                              cmul_zmm16r4(fthr,fthi,Ephr,Ephi,&t1r,&t1i);
                              Efzpr[j+3] = t1r;
                              Efzpi[j+3] = t1i;
                                                                      
                         }                           
                }
                
                
              
                
                
                 
                 
      
  
