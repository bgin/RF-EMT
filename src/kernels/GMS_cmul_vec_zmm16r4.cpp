



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


#include <immintrin.h>
#include "GMS_cmul_vec_zmm16r4.h"




                   void gms::math::cmul_zmm16r4_unroll_12x_u(const float * __restrict  xre,
                                                  const float * __restrict  xim,
                                                  const float * __restrict  yre,
                                                  const float * __restrict  yim,
                                                  float * __restrict        zre,
                                                  float * __restrict        zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) {return;}

                          __m512 zmm0,zmm1,zmm2,zmm3;
                          __m512 zmm4,zmm5,zmm6,zmm7;
                          __m512 zmm8,zmm9,zmm10,zmm11;
                          __m512 zmm12,zmm13,zmm14,zmm15;
                          __m512 zmm16,zmm17,zmm18,zmm19;
                          __m512 zmm20,zmm21,zmm22,zmm23;
                          __m512 zmm24,zmm25,zmm26,zmm27;
                          __m512 zmm28,zmm29;
                         int32_t i;

                         for(i = 0; (i+191) < n; i += 192) {
                            _mm_prefetch((const char *)&xre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+0],_MM_HINT_T0);
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&zim[i+0], zmm5);
                            _mm_prefetch((const char *)&xre[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+16],_MM_HINT_T0);
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_add_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_storeu_ps(&zim[i+16], zmm11);
                            _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+32],_MM_HINT_T0);
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_add_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&zim[i+32], zmm17);
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_storeu_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_add_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_storeu_ps(&zim[i+48], zmm23);
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+64],_MM_HINT_T0);
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             zmm28 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                                        _mm512_mul_ps(zmm26,zmm27));
                            _mm512_storeu_ps(&zre[i+64], zmm28);
                             zmm29  = _mm512_add_ps(_mm512_mul_ps(zmm26,zmm25),
                                                                        _mm512_mul_ps(zmm24,zmm27));
                            _mm512_storeu_ps(&zim[i+64], zmm29);
                            _mm_prefetch((const char *)&xre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+80],_MM_HINT_T0);
                            zmm1 = _mm512_loadu_ps(&xre[i+80]);
                            zmm2 = _mm512_loadu_ps(&yre[i+80]);
                            zmm3 = _mm512_loadu_ps(&xim[i+80]);
                            zmm4 = _mm512_loadu_ps(&yim[i+80]);
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_storeu_ps(&zre[i+80], zmm5);
                             zmm6  = _mm512_add_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_storeu_ps(&zim[i+80], zmm6);
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0);
                            zmm7 = _mm512_loadu_ps(&xre[i+96]);
                            zmm8 = _mm512_loadu_ps(&yre[i+96]);
                            zmm9 = _mm512_loadu_ps(&xim[i+96]);
                            zmm10= _mm512_loadu_ps(&yim[i+96]);
                            zmm11= _mm512_sub_ps(_mm512_mul_ps(zmm7,zmm8),
                                                       _mm512_mul_ps(zmm9,zmm10));
                            _mm512_storeu_ps(&zre[i+96], zmm11);
                             zmm12  = _mm512_add_ps(_mm512_mul_ps(zmm9,zmm8),
                                                           _mm512_mul_ps(zmm7,zmm10));
                            _mm512_storeu_ps(&zim[i+96], zmm12);
                            _mm_prefetch((const char *)&xre[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+112],_MM_HINT_T0);
                            zmm13 = _mm512_loadu_ps(&xre[i+112]);
                            zmm14 = _mm512_loadu_ps(&yre[i+112]);
                            zmm15 = _mm512_loadu_ps(&xim[i+112]);
                            zmm16 = _mm512_loadu_ps(&yim[i+112]);
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm13,zmm14),
                                                       _mm512_mul_ps(zmm15,zmm16));
                            _mm512_storeu_ps(&zre[i+112], zmm17);
                             zmm18  = _mm512_add_ps(_mm512_mul_ps(zmm15,zmm14),
                                                           _mm512_mul_ps(zmm13,zmm16));
                            _mm512_storeu_ps(&zim[i+112], zmm18);
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+128],_MM_HINT_T0);
                            zmm19 = _mm512_loadu_ps(&xre[i+128]);
                            zmm20 = _mm512_loadu_ps(&yre[i+128]);
                            zmm21 = _mm512_loadu_ps(&xim[i+128]);
                            zmm22 = _mm512_loadu_ps(&yim[i+128]);
                            zmm23 = _mm512_sub_ps(_mm512_mul_ps(zmm19,zmm20),
                                                       _mm512_mul_ps(zmm21,zmm22));
                            _mm512_storeu_ps(&zre[i+128], zmm23);
                             zmm24  = _mm512_add_ps(_mm512_mul_ps(zmm21,zmm20),
                                                           _mm512_mul_ps(zmm19,zmm22));
                            _mm512_storeu_ps(&zim[i+128], zmm24);
                            _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+144],_MM_HINT_T0);
                            zmm25 = _mm512_loadu_ps(&xre[i+144]);
                            zmm26 = _mm512_loadu_ps(&yre[i+144]);
                            zmm27 = _mm512_loadu_ps(&xim[i+144]);
                            zmm28 = _mm512_loadu_ps(&yim[i+144]);
                            zmm29 = _mm512_sub_ps(_mm512_mul_ps(zmm25,zmm26),
                                                       _mm512_mul_ps(zmm27,zmm28));
                            _mm512_storeu_ps(&zre[i+144], zmm29);
                             zmm1  = _mm512_add_ps(_mm512_mul_ps(zmm27,zmm26),
                                                           _mm512_mul_ps(zmm25,zmm28));
                            _mm512_storeu_ps(&zim[i+144], zmm1);
                            _mm_prefetch((const char *)&xre[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+160],_MM_HINT_T0);
                            zmm2 = _mm512_loadu_ps(&xre[i+160]);
                            zmm3 = _mm512_loadu_ps(&yre[i+160]);
                            zmm4 = _mm512_loadu_ps(&xim[i+160]);
                            zmm5 = _mm512_loadu_ps(&yim[i+160]);
                            zmm6 = _mm512_sub_ps(_mm512_mul_ps(zmm2,zmm3),
                                                       _mm512_mul_ps(zmm4,zmm5));
                            _mm512_storeu_ps(&zre[i+160], zmm6);
                             zmm7  = _mm512_add_ps(_mm512_mul_ps(zmm4,zmm3),
                                                           _mm512_mul_ps(zmm2,zmm5));
                            _mm512_storeu_ps(&zim[i+160], zmm7);
                            _mm_prefetch((const char *)&xre[i+176],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+176],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+176],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+176],_MM_HINT_T0);
                            zmm8 = _mm512_loadu_ps(&xre[i+176]);
                            zmm9 = _mm512_loadu_ps(&yre[i+176]);
                            zmm10= _mm512_loadu_ps(&xim[i+176]);
                            zmm11= _mm512_loadu_ps(&yim[i+176]);
                            zmm12= _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                       _mm512_mul_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&zre[i+176], zmm12);
                             zmm13 = _mm512_add_ps(_mm512_mul_ps(zmm10,zmm9),
                                                           _mm512_mul_ps(zmm8,zmm11));
                            _mm512_storeu_ps(&zim[i+176], zmm13);

                         }

                      for(; (i+159) < n; i += 160) {
                           
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&zim[i+0], zmm5);
                           
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_add_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_storeu_ps(&zim[i+16], zmm11);
                           
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_add_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&zim[i+32], zmm17);
                           
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_storeu_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_add_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_storeu_ps(&zim[i+48], zmm23);
                          
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             zmm28 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                                        _mm512_mul_ps(zmm26,zmm27));
                            _mm512_storeu_ps(&zre[i+64], zmm28);
                             zmm29  = _mm512_add_ps(_mm512_mul_ps(zmm26,zmm25),
                                                                        _mm512_mul_ps(zmm24,zmm27));
                            _mm512_storeu_ps(&zim[i+64], zmm29);
                           
                            zmm1 = _mm512_loadu_ps(&xre[i+80]);
                            zmm2 = _mm512_loadu_ps(&yre[i+80]);
                            zmm3 = _mm512_loadu_ps(&xim[i+80]);
                            zmm4 = _mm512_loadu_ps(&yim[i+80]);
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_storeu_ps(&zre[i+80], zmm5);
                             zmm6  = _mm512_add_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_storeu_ps(&zim[i+80], zmm6);
                           
                            zmm7 = _mm512_loadu_ps(&xre[i+96]);
                            zmm8 = _mm512_loadu_ps(&yre[i+96]);
                            zmm9 = _mm512_loadu_ps(&xim[i+96]);
                            zmm10= _mm512_loadu_ps(&yim[i+96]);
                            zmm11= _mm512_sub_ps(_mm512_mul_ps(zmm7,zmm8),
                                                       _mm512_mul_ps(zmm9,zmm10));
                            _mm512_storeu_ps(&zre[i+96], zmm11);
                             zmm12  = _mm512_add_ps(_mm512_mul_ps(zmm9,zmm8),
                                                           _mm512_mul_ps(zmm7,zmm10));
                            _mm512_storeu_ps(&zim[i+96], zmm12);
                           
                            zmm13 = _mm512_loadu_ps(&xre[i+112]);
                            zmm14 = _mm512_loadu_ps(&yre[i+112]);
                            zmm15 = _mm512_loadu_ps(&xim[i+112]);
                            zmm16 = _mm512_loadu_ps(&yim[i+112]);
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm13,zmm14),
                                                       _mm512_mul_ps(zmm15,zmm16));
                            _mm512_storeu_ps(&zre[i+112], zmm17);
                             zmm18  = _mm512_add_ps(_mm512_mul_ps(zmm15,zmm14),
                                                           _mm512_mul_ps(zmm13,zmm16));
                            _mm512_storeu_ps(&zim[i+112], zmm18);
                            
                            zmm19 = _mm512_loadu_ps(&xre[i+128]);
                            zmm20 = _mm512_loadu_ps(&yre[i+128]);
                            zmm21 = _mm512_loadu_ps(&xim[i+128]);
                            zmm22 = _mm512_loadu_ps(&yim[i+128]);
                            zmm23 = _mm512_sub_ps(_mm512_mul_ps(zmm19,zmm20),
                                                       _mm512_mul_ps(zmm21,zmm22));
                            _mm512_storeu_ps(&zre[i+128], zmm23);
                             zmm24  = _mm512_add_ps(_mm512_mul_ps(zmm21,zmm20),
                                                           _mm512_mul_ps(zmm19,zmm22));
                            _mm512_storeu_ps(&zim[i+128], zmm24);
                           
                            zmm25 = _mm512_loadu_ps(&xre[i+144]);
                            zmm26 = _mm512_loadu_ps(&yre[i+144]);
                            zmm27 = _mm512_loadu_ps(&xim[i+144]);
                            zmm28 = _mm512_loadu_ps(&yim[i+144]);
                            zmm29 = _mm512_sub_ps(_mm512_mul_ps(zmm25,zmm26),
                                                       _mm512_mul_ps(zmm27,zmm28));
                            _mm512_storeu_ps(&zre[i+144], zmm29);
                             zmm1  = _mm512_add_ps(_mm512_mul_ps(zmm27,zmm26),
                                                           _mm512_mul_ps(zmm25,zmm28));
                            _mm512_storeu_ps(&zim[i+144], zmm1);
                      }

                    for(; (i+127) < n; i += 128) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_add_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_storeu_ps(&zim[i+16], zmm11);
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_add_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&zim[i+32], zmm17);
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_storeu_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_add_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_storeu_ps(&zim[i+48], zmm23);
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             zmm28 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                                        _mm512_mul_ps(zmm26,zmm27));
                            _mm512_storeu_ps(&zre[i+64], zmm28);
                             zmm29  = _mm512_add_ps(_mm512_mul_ps(zmm26,zmm25),
                                                                        _mm512_mul_ps(zmm24,zmm27));
                            _mm512_storeu_ps(&zim[i+64], zmm29);
                            zmm1 = _mm512_loadu_ps(&xre[i+80]);
                            zmm2 = _mm512_loadu_ps(&yre[i+80]);
                            zmm3 = _mm512_loadu_ps(&xim[i+80]);
                            zmm4 = _mm512_loadu_ps(&yim[i+80]);
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_storeu_ps(&zre[i+80], zmm5);
                             zmm6  = _mm512_add_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_storeu_ps(&zim[i+80], zmm6);
                            zmm7 = _mm512_loadu_ps(&xre[i+96]);
                            zmm8 = _mm512_loadu_ps(&yre[i+96]);
                            zmm9 = _mm512_loadu_ps(&xim[i+96]);
                            zmm10= _mm512_loadu_ps(&yim[i+96]);
                            zmm11= _mm512_sub_ps(_mm512_mul_ps(zmm7,zmm8),
                                                       _mm512_mul_ps(zmm9,zmm10));
                            _mm512_storeu_ps(&zre[i+96], zmm11);
                             zmm12  = _mm512_add_ps(_mm512_mul_ps(zmm9,zmm8),
                                                           _mm512_mul_ps(zmm7,zmm10));
                            _mm512_storeu_ps(&zim[i+96], zmm12);
                            zmm13 = _mm512_loadu_ps(&xre[i+112]);
                            zmm14 = _mm512_loadu_ps(&yre[i+112]);
                            zmm15 = _mm512_loadu_ps(&xim[i+112]);
                            zmm16 = _mm512_loadu_ps(&yim[i+112]);
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm13,zmm14),
                                                       _mm512_mul_ps(zmm15,zmm16));
                            _mm512_storeu_ps(&zre[i+112], zmm17);
                             zmm18  = _mm512_add_ps(_mm512_mul_ps(zmm15,zmm14),
                                                           _mm512_mul_ps(zmm13,zmm16));
                            _mm512_storeu_ps(&zim[i+112], zmm18);

                    }

                    for(; (i+95) < n; i += 96) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_add_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_storeu_ps(&zim[i+16], zmm11);
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_add_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&zim[i+32], zmm17);
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_storeu_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_add_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_storeu_ps(&zim[i+48], zmm23);
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             zmm28 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                                        _mm512_mul_ps(zmm26,zmm27));
                            _mm512_storeu_ps(&zre[i+64], zmm28);
                             zmm29  = _mm512_add_ps(_mm512_mul_ps(zmm26,zmm25),
                                                                        _mm512_mul_ps(zmm24,zmm27));
                            _mm512_storeu_ps(&zim[i+64], zmm29);
                            zmm1 = _mm512_loadu_ps(&xre[i+80]);
                            zmm2 = _mm512_loadu_ps(&yre[i+80]);
                            zmm3 = _mm512_loadu_ps(&xim[i+80]);
                            zmm4 = _mm512_loadu_ps(&yim[i+80]);
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_storeu_ps(&zre[i+80], zmm5);
                             zmm6  = _mm512_add_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_storeu_ps(&zim[i+80], zmm6);
                   }

                    for(; (i+63) < n; i += 64) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_add_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_storeu_ps(&zim[i+16], zmm11);
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_add_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&zim[i+32], zmm17);
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_storeu_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_add_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_storeu_ps(&zim[i+48], zmm23);

                    }

                    for(; (i+31) < n; i += 32) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_add_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_storeu_ps(&zim[i+16], zmm11);
                    }

                    for(; (i+15) < n; i += 16) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&zim[i+0], zmm5);

                    }

                    for(; (i+0) < n; i += 1) {
                         zre[i] = (xre[i] * yre[i]) - (xim[i] * yim[i]);
		             zim[i] = (xim[i] * yre[i]) + (xre[i] * yim[i]);
                    }
               }


                 
                   void gms::math::cmul_zmm16r4_unroll_12x_a(const float * __restrict  __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict  __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict  __ATTR_ALIGN__(64) yre,
                                                  const float * __restrict  __ATTR_ALIGN__(64) yim,
                                                  float * __restrict        __ATTR_ALIGN__(64) zre,
                                                  float * __restrict        __ATTR_ALIGN__(64) zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) {return;}

                          __m512 zmm0,zmm1,zmm2,zmm3;
                          __m512 zmm4,zmm5,zmm6,zmm7;
                          __m512 zmm8,zmm9,zmm10,zmm11;
                          __m512 zmm12,zmm13,zmm14,zmm15;
                          __m512 zmm16,zmm17,zmm18,zmm19;
                          __m512 zmm20,zmm21,zmm22,zmm23;
                          __m512 zmm24,zmm25,zmm26,zmm27;
                          __m512 zmm28,zmm29;
                         int32_t i;

                         for(i = 0; (i+191) < n; i += 192) {
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_load_ps(&xre[i+16]);
                             zmm7  = _mm512_load_ps(&yre[i+16]);
                             zmm8  = _mm512_load_ps(&xim[i+16]);
                             zmm9  = _mm512_load_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_store_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_mul_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_store_ps(&zim[i+16], zmm11);
                             zmm12 = _mm512_load_ps(&xre[i+32]);
                             zmm13 = _mm512_load_ps(&yre[i+32]);
                             zmm14 = _mm512_load_ps(&xim[i+32]);
                             zmm15 = _mm512_load_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&zim[i+32], zmm17);
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0);
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_store_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_store_ps(&zim[i+48], zmm23);
                             zmm24 = _mm512_load_ps(&xre[i+64]);
                             zmm25 = _mm512_load_ps(&yre[i+64]);
                             zmm26 = _mm512_load_ps(&xim[i+64]);
                             zmm27 = _mm512_load_ps(&yim[i+64]);
                             zmm28 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                                        _mm512_mul_ps(zmm26,zmm27));
                            _mm512_store_ps(&zre[i+64], zmm28);
                             zmm29  = _mm512_mul_ps(_mm512_mul_ps(zmm26,zmm25),
                                                                        _mm512_mul_ps(zmm24,zmm27));
                            _mm512_store_ps(&zim[i+64], zmm29);
                            zmm1 = _mm512_load_ps(&xre[i+80]);
                            zmm2 = _mm512_load_ps(&yre[i+80]);
                            zmm3 = _mm512_load_ps(&xim[i+80]);
                            zmm4 = _mm512_load_ps(&yim[i+80]);
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_store_ps(&zre[i+80], zmm5);
                             zmm6  = _mm512_mul_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_storeu_ps(&zim[i+80], zmm6);
                            _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+144],_MM_HINT_T0);
                            zmm7 = _mm512_load_ps(&xre[i+96]);
                            zmm8 = _mm512_load_ps(&yre[i+96]);
                            zmm9 = _mm512_load_ps(&xim[i+96]);
                            zmm10= _mm512_load_ps(&yim[i+96]);
                            zmm11= _mm512_sub_ps(_mm512_mul_ps(zmm7,zmm8),
                                                       _mm512_mul_ps(zmm9,zmm10));
                            _mm512_store_ps(&zre[i+96], zmm11);
                             zmm12  = _mm512_mul_ps(_mm512_mul_ps(zmm9,zmm8),
                                                           _mm512_mul_ps(zmm7,zmm10));
                            _mm512_storeu_ps(&zim[i+96], zmm12);
                            zmm13 = _mm512_load_ps(&xre[i+112]);
                            zmm14 = _mm512_load_ps(&yre[i+112]);
                            zmm15 = _mm512_load_ps(&xim[i+112]);
                            zmm16 = _mm512_load_ps(&yim[i+112]);
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm13,zmm14),
                                                       _mm512_mul_ps(zmm15,zmm16));
                            _mm512_store_ps(&zre[i+112], zmm17;
                             zmm18  = _mm512_mul_ps(_mm512_mul_ps(zmm15,zmm14),
                                                           _mm512_mul_ps(zmm13,zmm16));
                            _mm512_store_ps(&zim[i+112], zmm18);
                            zmm19 = _mm512_load_ps(&xre[i+128]);
                            zmm20 = _mm512_load_ps(&yre[i+128]);
                            zmm21 = _mm512_load_ps(&xim[i+128]);
                            zmm22 = _mm512_load_ps(&yim[i+128]);
                            zmm23 = _mm512_sub_ps(_mm512_mul_ps(zmm19,zmm20),
                                                       _mm512_mul_ps(zmm21,zmm22));
                            _mm512_store_ps(&zre[i+128], zmm23);
                             zmm24  = _mm512_mul_ps(_mm512_mul_ps(zmm21,zmm20),
                                                           _mm512_mul_ps(zmm19,zmm22));
                            _mm512_store_ps(&zim[i+128], zmm24);
                            _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+192],_MM_HINT_T0);
                            zmm25 = _mm512_load_ps(&xre[i+144]);
                            zmm26 = _mm512_load_ps(&yre[i+144]);
                            zmm27 = _mm512_load_ps(&xim[i+144]);
                            zmm28 = _mm512_load_ps(&yim[i+144]);
                            zmm29 = _mm512_sub_ps(_mm512_mul_ps(zmm25,zmm26),
                                                       _mm512_mul_ps(zmm27,zmm28));
                            _mm512_store_ps(&zre[i+144], zmm29);
                             zmm1  = _mm512_mul_ps(_mm512_mul_ps(zmm27,zmm26),
                                                           _mm512_mul_ps(zmm25,zmm28));
                            _mm512_store_ps(&zim[i+144], zmm1);
                            zmm2 = _mm512_load_ps(&xre[i+160]);
                            zmm3 = _mm512_load_ps(&yre[i+160]);
                            zmm4 = _mm512_load_ps(&xim[i+160]);
                            zmm5 = _mm512_load_ps(&yim[i+160]);
                            zmm6 = _mm512_sub_ps(_mm512_mul_ps(zmm2,zmm3),
                                                       _mm512_mul_ps(zmm4,zmm5));
                            _mm512_store_ps(&zre[i+160], zmm6);
                             zmm7  = _mm512_mul_ps(_mm512_mul_ps(zmm4,zmm3),
                                                           _mm512_mul_ps(zmm2,zmm5));
                            _mm512_store_ps(&zim[i+160], zmm7);
                            zmm8 = _mm512_load_ps(&xre[i+176]);
                            zmm9 = _mm512_load_ps(&yre[i+176]);
                            zmm10= _mm512_load_ps(&xim[i+176]);
                            zmm11= _mm512_load_ps(&yim[i+176]);
                            zmm12= _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                       _mm512_mul_ps(zmm10,zmm11));
                            _mm512_store_ps(&zre[i+176], zmm12);
                             zmm13 = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                           _mm512_mul_ps(zmm8,zmm11));
                            _mm512_store_ps(&zim[i+176], zmm13);

                         }

                      for(; (i+159) < n; i += 160) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_load_ps(&xre[i+16]);
                             zmm7  = _mm512_load_ps(&yre[i+16]);
                             zmm8  = _mm512_load_ps(&xim[i+16]);
                             zmm9  = _mm512_load_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_store_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_mul_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_store_ps(&zim[i+16], zmm11);
                             zmm12 = _mm512_load_ps(&xre[i+32]);
                             zmm13 = _mm512_load_ps(&yre[i+32]);
                             zmm14 = _mm512_load_ps(&xim[i+32]);
                             zmm15 = _mm512_load_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&zim[i+32], zmm17);
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_store_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_store_ps(&zim[i+48], zmm23);
                             zmm24 = _mm512_load_ps(&xre[i+64]);
                             zmm25 = _mm512_load_ps(&yre[i+64]);
                             zmm26 = _mm512_load_ps(&xim[i+64]);
                             zmm27 = _mm512_load_ps(&yim[i+64]);
                             zmm28 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                                        _mm512_mul_ps(zmm26,zmm27));
                            _mm512_storeu_ps(&zre[i+64], zmm28);
                             zmm29  = _mm512_mul_ps(_mm512_mul_ps(zmm26,zmm25),
                                                                        _mm512_mul_ps(zmm24,zmm27));
                            _mm512_storeu_ps(&zim[i+64], zmm29);
                            zmm1 = _mm512_load_ps(&xre[i+80]);
                            zmm2 = _mm512_load_ps(&yre[i+80]);
                            zmm3 = _mm512_load_ps(&xim[i+80]);
                            zmm4 = _mm512_load_ps(&yim[i+80]);
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_store_ps(&zre[i+80], zmm5);
                             zmm6  = _mm512_mul_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_store_ps(&zim[i+80], zmm6);
                            zmm7 = _mm512_load_ps(&xre[i+96]);
                            zmm8 = _mm512_load_ps(&yre[i+96]);
                            zmm9 = _mm512_load_ps(&xim[i+96]);
                            zmm10= _mm512_load_ps(&yim[i+96]);
                            zmm11= _mm512_sub_ps(_mm512_mul_ps(zmm7,zmm8),
                                                       _mm512_mul_ps(zmm9,zmm10));
                            _mm512_store_ps(&zre[i+96], zmm11);
                             zmm12  = _mm512_mul_ps(_mm512_mul_ps(zmm9,zmm8),
                                                           _mm512_mul_ps(zmm7,zmm10));
                            _mm512_store_ps(&zim[i+96], zmm12);
                            zmm13 = _mm512_load_ps(&xre[i+112]);
                            zmm14 = _mm512_load_ps(&yre[i+112]);
                            zmm15 = _mm512_load_ps(&xim[i+112]);
                            zmm16 = _mm512_load_ps(&yim[i+112]);
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm13,zmm14),
                                                       _mm512_mul_ps(zmm15,zmm16));
                            _mm512_store_ps(&zre[i+112], zmm17;
                             zmm18  = _mm512_mul_ps(_mm512_mul_ps(zmm15,zmm14),
                                                           _mm512_mul_ps(zmm13,zmm16));
                            _mm512_store_ps(&zim[i+112], zmm18);
                            zmm19 = _mm512_load_ps(&xre[i+128]);
                            zmm20 = _mm512_load_ps(&yre[i+128]);
                            zmm21 = _mm512_load_ps(&xim[i+128]);
                            zmm22 = _mm512_load_ps(&yim[i+128]);
                            zmm23 = _mm512_sub_ps(_mm512_mul_ps(zmm19,zmm20),
                                                       _mm512_mul_ps(zmm21,zmm22));
                            _mm512_store_ps(&zre[i+128], zmm23);
                             zmm24  = _mm512_mul_ps(_mm512_mul_ps(zmm21,zmm20),
                                                           _mm512_mul_ps(zmm19,zmm22));
                            _mm512_store_ps(&zim[i+128], zmm24);
                            zmm25 = _mm512_load_ps(&xre[i+144]);
                            zmm26 = _mm512_load_ps(&yre[i+144]);
                            zmm27 = _mm512_load_ps(&xim[i+144]);
                            zmm28 = _mm512_load_ps(&yim[i+144]);
                            zmm29 = _mm512_sub_ps(_mm512_mul_ps(zmm25,zmm26),
                                                       _mm512_mul_ps(zmm27,zmm28));
                            _mm512_store_ps(&zre[i+144], zmm29);
                             zmm1  = _mm512_mul_ps(_mm512_mul_ps(zmm27,zmm26),
                                                           _mm512_mul_ps(zmm25,zmm28));
                            _mm512_store_ps(&zim[i+144], zmm1);

                      }

                    for(; (i+127) < n; i += 128) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_load_ps(&xre[i+16]);
                             zmm7  = _mm512_load_ps(&yre[i+16]);
                             zmm8  = _mm512_load_ps(&xim[i+16]);
                             zmm9  = _mm512_load_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_store_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_mul_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_store_ps(&zim[i+16], zmm11);
                             zmm12 = _mm512_load_ps(&xre[i+32]);
                             zmm13 = _mm512_load_ps(&yre[i+32]);
                             zmm14 = _mm512_load_ps(&xim[i+32]);
                             zmm15 = _mm512_load_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&zim[i+32], zmm17);
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_store_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_store_ps(&zim[i+48], zmm23);
                             zmm24 = _mm512_load_ps(&xre[i+64]);
                             zmm25 = _mm512_load_ps(&yre[i+64]);
                             zmm26 = _mm512_load_ps(&xim[i+64]);
                             zmm27 = _mm512_load_ps(&yim[i+64]);
                             zmm28 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                                        _mm512_mul_ps(zmm26,zmm27));
                            _mm512_store_ps(&zre[i+64], zmm28);
                             zmm29  = _mm512_mul_ps(_mm512_mul_ps(zmm26,zmm25),
                                                                        _mm512_mul_ps(zmm24,zmm27));
                            _mm512_store_ps(&zim[i+64], zmm29);
                            zmm1 = _mm512_load_ps(&xre[i+80]);
                            zmm2 = _mm512_load_ps(&yre[i+80]);
                            zmm3 = _mm512_load_ps(&xim[i+80]);
                            zmm4 = _mm512_load_ps(&yim[i+80]);
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_store_ps(&zre[i+80], zmm5);
                             zmm6  = _mm512_mul_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_store_ps(&zim[i+80], zmm6);
                            zmm7 = _mm512_load_ps(&xre[i+96]);
                            zmm8 = _mm512_load_ps(&yre[i+96]);
                            zmm9 = _mm512_load_ps(&xim[i+96]);
                            zmm10= _mm512_load_ps(&yim[i+96]);
                            zmm11= _mm512_sub_ps(_mm512_mul_ps(zmm7,zmm8),
                                                       _mm512_mul_ps(zmm9,zmm10));
                            _mm512_store_ps(&zre[i+96], zmm11);
                             zmm12  = _mm512_mul_ps(_mm512_mul_ps(zmm9,zmm8),
                                                           _mm512_mul_ps(zmm7,zmm10));
                            _mm512_storeu_ps(&zim[i+96], zmm12);
                            zmm13 = _mm512_load_ps(&xre[i+112]);
                            zmm14 = _mm512_load_ps(&yre[i+112]);
                            zmm15 = _mm512_load_ps(&xim[i+112]);
                            zmm16 = _mm512_load_ps(&yim[i+112]);
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm13,zmm14),
                                                       _mm512_mul_ps(zmm15,zmm16));
                            _mm512_store_ps(&zre[i+112], zmm17;
                             zmm18  = _mm512_mul_ps(_mm512_mul_ps(zmm15,zmm14),
                                                           _mm512_mul_ps(zmm13,zmm16));
                            _mm512_store_ps(&zim[i+112], zmm18);

                    }

                    for(; (i+95) < n; i += 96) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_load_ps(&xre[i+16]);
                             zmm7  = _mm512_load_ps(&yre[i+16]);
                             zmm8  = _mm512_load_ps(&xim[i+16]);
                             zmm9  = _mm512_load_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_store_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_mul_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_store_ps(&zim[i+16], zmm11);
                             zmm12 = _mm512_load_ps(&xre[i+32]);
                             zmm13 = _mm512_load_ps(&yre[i+32]);
                             zmm14 = _mm512_load_ps(&xim[i+32]);
                             zmm15 = _mm512_load_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&zim[i+32], zmm17);
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_store_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_storeu_ps(&zim[i+48], zmm23);
                             zmm24 = _mm512_load_ps(&xre[i+64]);
                             zmm25 = _mm512_load_ps(&yre[i+64]);
                             zmm26 = _mm512_load_ps(&xim[i+64]);
                             zmm27 = _mm512_load_ps(&yim[i+64]);
                             zmm28 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                                        _mm512_mul_ps(zmm26,zmm27));
                            _mm512_store_ps(&zre[i+64], zmm28);
                             zmm29  = _mm512_mul_ps(_mm512_mul_ps(zmm26,zmm25),
                                                                        _mm512_mul_ps(zmm24,zmm27));
                            _mm512_storeu_ps(&zim[i+64], zmm29);
                            zmm1 = _mm512_load_ps(&xre[i+80]);
                            zmm2 = _mm512_load_ps(&yre[i+80]);
                            zmm3 = _mm512_load_ps(&xim[i+80]);
                            zmm4 = _mm512_load_ps(&yim[i+80]);
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_store_ps(&zre[i+80], zmm5);
                             zmm6  = _mm512_mul_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_store_ps(&zim[i+80], zmm6);
                   }

                    for(; (i+63) < n; i += 64) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_load_ps(&xre[i+16]);
                             zmm7  = _mm512_load_ps(&yre[i+16]);
                             zmm8  = _mm512_load_ps(&xim[i+16]);
                             zmm9  = _mm512_load_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_store_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_mul_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_storeu_ps(&zim[i+16], zmm11);
                             zmm12 = _mm512_load_ps(&xre[i+32]);
                             zmm13 = _mm512_load_ps(&yre[i+32]);
                             zmm14 = _mm512_load_ps(&xim[i+32]);
                             zmm15 = _mm512_load_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&zim[i+32], zmm17);
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_store_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_store_ps(&zim[i+48], zmm23);

                    }

                    for(; (i+31) < n; i += 32) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_load_ps(&xre[i+16]);
                             zmm7  = _mm512_load_ps(&yre[i+16]);
                             zmm8  = _mm512_load_ps(&xim[i+16]);
                             zmm9  = _mm512_load_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_store_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_mul_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_store_ps(&zim[i+16], zmm11);
                    }

                    for(; (i+15) < n; i += 16) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&zim[i+0], zmm5);

                    }

                    for(; (i+0) < n; i += 1) {
                         zre[i] = (xre[i] * yre[i]) - (xim[i] * yim[i]);
		         zim[i] = (xim[i] * yre[i]) + (xre[i] * yim[i]);
                    }
               }

         
                 
                   void gms::math::cmul_zmm16r4_unroll_16x_u(const float * __restrict xre,
                                                  const float * __restrict xim, 
                                                  const __m512             vs,
                                                  float  * __restrict zre,
                                                  float  * __restrict zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                           __m512 zmm0,zmm1,zmm2,zmm3;
                           __m512 zmm4,zmm5,zmm6,zmm7;
                           __m512 zmm8,zmm9,zmm10,zmm11;
                           __m512 zmm12,zmm13,zmm14,zmm15;
                           __m512 zmm16,zmm17,zmm18,zmm19;
                           __m512 zmm20,zmm21,zmm22,zmm23;
                           __m512 zmm24,zmm25,zmm26,zmm27;
                           __m512 zmm28,zmm29,zmm30,zmm31;
                           __m512 zmmx    = vs;
                          const float * __restrict pzmm =  (const float*)&vs[0];
                          for(i = 0; (i+255) < n; i += 256) {
                              _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                               zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                               zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                               zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                               zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                               zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                               zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                               zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                               zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));
                               zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_mul_ps(zmm8,zmmx));
                               zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_mul_ps(zmm9,zmmx));
                              _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                               zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_mul_ps(zmm10,zmmx));
                               zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_mul_ps(zmm11,zmmx));
                               zmm12 = _mm512_loadu_ps(&xre[i+96]);
                              _mm512_storeu_ps(&zre[i+96], _mm512_mul_ps(zmm12,zmmx));
                               zmm13 = _mm512_loadu_ps(&xim[i+96]);
                              _mm512_storeu_ps(&zim[i+96], _mm512_mul_ps(zmm13,zmmx));
                               zmm14 = _mm512_loadu_ps(&xre[i+112]);
                              _mm512_storeu_ps(&zre[i+112], _mm512_mul_ps(zmm14,zmmx));
                               zmm15 = _mm512_loadu_ps(&xim[i+112]);
                              _mm512_storeu_ps(&zim[i+112], _mm512_mul_ps(zmm15,zmmx));
                               zmm16 = _mm512_loadu_ps(&xre[i+128]);
                              _mm512_storeu_ps(&zre[i+128], _mm512_mul_ps(zmm16,zmmx));
                               zmm17 = _mm512_loadu_ps(&xim[i+128]);
                              _mm512_storeu_ps(&zim[i+128], _mm512_mul_ps(zmm17,zmmx));
                              _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                               zmm18 = _mm512_loadu_ps(&xre[i+144]);
                              _mm512_storeu_ps(&zre[i+144], _mm512_mul_ps(zmm18,zmmx));
                               zmm19 = _mm512_loadu_ps(&xim[i+144]);
                              _mm512_storeu_ps(&zim[i+144], _mm512_mul_ps(zmm19,zmmx));
                               zmm20 = _mm512_loadu_ps(&xre[i+160]);
                              _mm512_storeu_ps(&zre[i+160], _mm512_mul_ps(zmm20,zmmx));
                               zmm21 = _mm512_loadu_ps(&xim[i+160]);
                              _mm512_storeu_ps(&zim[i+160], _mm512_mul_ps(zmm21,zmmx));
                               zmm22 = _mm512_loadu_ps(&xre[i+176]);
                              _mm512_storeu_ps(&zre[i+176], _mm512_mul_ps(zmm22,zmmx));
                               zmm23 = _mm512_loadu_ps(&xim[i+176]);
                              _mm512_storeu_ps(&zim[i+176], _mm512_mul_ps(zmm23,zmmx));
                               zmm24 = _mm512_loadu_ps(&xre[i+192]);
                              _mm512_storeu_ps(&zre[i+192], _mm512_mul_ps(zmm24,zmmx));
                               zmm25 = _mm512_loadu_ps(&xim[i+192]);
                              _mm512_storeu_ps(&zim[i+192], _mm512_mul_ps(zmm25,zmmx));
                              _mm_prefetch((const char *)&xre[i+255],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+255],_MM_HINT_T0);
                               zmm26 = _mm512_loadu_ps(&xre[i+208]);
                              _mm512_storeu_ps(&zre[i+208], _mm512_mul_ps(zmm26,zmmx));
                               zmm27 = _mm512_loadu_ps(&xim[i+208]);
                              _mm512_storeu_ps(&zim[i+208], _mm512_mul_ps(zmm27,zmmx));
                               zmm28 = _mm512_loadu_ps(&xre[i+224]);
                              _mm512_storeu_ps(&zre[i+224], _mm512_mul_ps(zmm28,zmmx));
                               zmm29 = _mm512_loadu_ps(&xim[i+224]);
                              _mm512_storeu_ps(&zim[i+224], _mm512_mul_ps(zmm29,zmmx));
                               zmm30 = _mm512_loadu_ps(&xre[i+240]);
                              _mm512_storeu_ps(&zre[i+240], _mm512_mul_ps(zmm30,zmmx));
                               zmm31 = _mm512_loadu_ps(&xim[i+240]);
                              _mm512_storeu_ps(&zim[i+240], _mm512_mul_ps(zmm31,zmmx)); 

                         }

                        for(; (i+191) < n; i += 192) {
                               zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                               zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                               zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                               zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                               zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                               zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                               zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                               zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));
                               zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_mul_ps(zmm8,zmmx));
                               zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_mul_ps(zmm9,zmmx));
                              zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_mul_ps(zmm10,zmmx));
                               zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_mul_ps(zmm11,zmmx));
                               zmm12 = _mm512_loadu_ps(&xre[i+96]);
                              _mm512_storeu_ps(&zre[i+96], _mm512_mul_ps(zmm12,zmmx));
                               zmm13 = _mm512_loadu_ps(&xim[i+96]);
                              _mm512_storeu_ps(&zim[i+96], _mm512_mul_ps(zmm13,zmmx));
                               zmm14 = _mm512_loadu_ps(&xre[i+112]);
                              _mm512_storeu_ps(&zre[i+112], _mm512_mul_ps(zmm14,zmmx));
                               zmm15 = _mm512_loadu_ps(&xim[i+112]);
                              _mm512_storeu_ps(&zim[i+112], _mm512_mul_ps(zmm15,zmmx));
                               zmm16 = _mm512_loadu_ps(&xre[i+128]);
                              _mm512_storeu_ps(&zre[i+128], _mm512_mul_ps(zmm16,zmmx));
                               zmm17 = _mm512_loadu_ps(&xim[i+128]);
                              _mm512_storeu_ps(&zim[i+128], _mm512_mul_ps(zmm17,zmmx));
                               zmm18 = _mm512_loadu_ps(&xre[i+144]);
                              _mm512_storeu_ps(&zre[i+144], _mm512_mul_ps(zmm18,zmmx));
                               zmm19 = _mm512_loadu_ps(&xim[i+144]);
                              _mm512_storeu_ps(&zim[i+144], _mm512_mul_ps(zmm19,zmmx));
                               zmm20 = _mm512_loadu_ps(&xre[i+160]);
                              _mm512_storeu_ps(&zre[i+160], _mm512_mul_ps(zmm20,zmmx));
                               zmm21 = _mm512_loadu_ps(&xim[i+160]);
                              _mm512_storeu_ps(&zim[i+160], _mm512_mul_ps(zmm21,zmmx));
                               zmm22 = _mm512_loadu_ps(&xre[i+176]);
                              _mm512_storeu_ps(&zre[i+176], _mm512_mul_ps(zmm22,zmmx));
                               zmm23 = _mm512_loadu_ps(&xim[i+176]);
                              _mm512_storeu_ps(&zim[i+176], _mm512_mul_ps(zmm23,zmmx));

                        }

                        for(; (i+127) < n; i += 128) {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                              zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                              zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                              zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                              zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                              zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));
                              zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_mul_ps(zmm8,zmmx));
                              zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_mul_ps(zmm9,zmmx));
                              zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_mul_ps(zmm10,zmmx));
                              zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_mul_ps(zmm11,zmmx));
                              zmm12 = _mm512_loadu_ps(&xre[i+96]);
                              _mm512_storeu_ps(&zre[i+96], _mm512_mul_ps(zmm12,zmmx));
                              zmm13 = _mm512_loadu_ps(&xim[i+96]);
                              _mm512_storeu_ps(&zim[i+96], _mm512_mul_ps(zmm13,zmmx));
                              zmm14 = _mm512_loadu_ps(&xre[i+112]);
                              _mm512_storeu_ps(&zre[i+112], _mm512_mul_ps(zmm14,zmmx));
                              zmm15 = _mm512_loadu_ps(&xim[i+112]);
                              _mm512_storeu_ps(&zim[i+112], _mm512_mul_ps(zmm15,zmmx));

                        }

                       for(; (i+63) < n; i += 64) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                             zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                             zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                             zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                             zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                             zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                             zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                             zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));

                       }

                      for(; (i+31) < n; i += 32) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                             zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                             zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                             zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                      }

                     for(; (i+15) < n; i += 16) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                             zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                     }
                    
                    
                    for(; (i+0) < n; i += 1) {
                          const float s = pzmm[i];
                          zre[i] = xre[i] * s;
                          zim[i] = zim[i] * s;
                    }
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void cmul_zmm16r4_unroll_16x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim, 
                                                  const __m512             vs,
                                                  float  * __restrict __ATTR_ALIGN__(64) zre,
                                                  float  * __restrict __ATTR_ALIGN__(64) zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                           __m512 zmm0,zmm1,zmm2,zmm3;
                           __m512 zmm4,zmm5,zmm6,zmm7;
                           __m512 zmm8,zmm9,zmm10,zmm11;
                           __m512 zmm12,zmm13,zmm14,zmm15;
                           __m512 zmm16,zmm17,zmm18,zmm19;
                           __m512 zmm20,zmm21,zmm22,zmm23;
                           __m512 zmm24,zmm25,zmm26,zmm27;
                           __m512 zmm28,zmm29,zmm30,zmm31;
                           __m512 zmmx    = vs;
                          const float * __restrict pzmm =  (const float*)&vs[0];
                          for(i = 0; (i+255) < n; i += 256) {
                              _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                              zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                              zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));
                              zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_mul_ps(zmm8,zmmx));
                              zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_mul_ps(zmm9,zmmx));
                              _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                              zmm10 = _mm512_load_ps(&xre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_mul_ps(zmm10,zmmx));
                              zmm11 = _mm512_load_ps(&xim[i+80]);
                              _mm512_store_ps(&zim[i+80], _mm512_mul_ps(zmm11,zmmx));
                              zmm12 = _mm512_load_ps(&xre[i+96]);
                              _mm512_store_ps(&zre[i+96], _mm512_mul_ps(zmm12,zmmx));
                              zmm13 = _mm512_load_ps(&xim[i+96]);
                              _mm512_store_ps(&zim[i+96], _mm512_mul_ps(zmm13,zmmx));
                              zmm14 = _mm512_load_ps(&xre[i+112]);
                              _mm512_store_ps(&zre[i+112], _mm512_mul_ps(zmm14,zmmx));
                              zmm15 = _mm512_load_ps(&xim[i+112]);
                              _mm512_store_ps(&zim[i+112], _mm512_mul_ps(zmm15,zmmx));
                              zmm16 = _mm512_load_ps(&xre[i+128]);
                              _mm512_store_ps(&zre[i+128], _mm512_mul_ps(zmm16,zmmx));
                              zmm17 = _mm512_load_ps(&xim[i+128]);
                              _mm512_store_ps(&zim[i+128], _mm512_mul_ps(zmm17,zmmx));
                              _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                              zmm18 = _mm512_load_ps(&xre[i+144]);
                              _mm512_store_ps(&zre[i+144], _mm512_mul_ps(zmm18,zmmx));
                              zmm19 = _mm512_load_ps(&xim[i+144]);
                              _mm512_store_ps(&zim[i+144], _mm512_mul_ps(zmm19,zmmx));
                              zmm20 = _mm512_load_ps(&xre[i+160]);
                              _mm512_store_ps(&zre[i+160], _mm512_mul_ps(zmm20,zmmx));
                              zmm21 = _mm512_load_ps(&xim[i+160]);
                              _mm512_store_ps(&zim[i+160], _mm512_mul_ps(zmm21,zmmx));
                              zmm22 = _mm512_load_ps(&xre[i+176]);
                              _mm512_store_ps(&zre[i+176], _mm512_mul_ps(zmm22,zmmx));
                              zmm23 = _mm512_load_ps(&xim[i+176]);
                              _mm512_store_ps(&zim[i+176], _mm512_mul_ps(zmm23,zmmx));
                              zmm24 = _mm512_load_ps(&xre[i+192]);
                              _mm512_store_ps(&zre[i+192], _mm512_mul_ps(zmm24,zmmx));
                              zmm25 = _mm512_load_ps(&xim[i+192]);
                              _mm512_storeu_ps(&zim[i+192], _mm512_mul_ps(zmm25,zmmx));
                              _mm_prefetch((const char *)&xre[i+255],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+255],_MM_HINT_T0);
                              zmm26 = _mm512_load_ps(&xre[i+208]);
                              _mm512_store_ps(&zre[i+208], _mm512_mul_ps(zmm26,zmmx));
                              zmm27 = _mm512_load_ps(&xim[i+208]);
                              _mm512_store_ps(&zim[i+208], _mm512_mul_ps(zmm27,zmmx));
                              zmm28 = _mm512_load_ps(&xre[i+224]);
                              _mm512_store_ps(&zre[i+224], _mm512_mul_ps(zmm28,zmmx));
                              zmm29 = _mm512_load_ps(&xim[i+224]);
                              _mm512_store_ps(&zim[i+224], _mm512_mul_ps(zmm29,zmmx));
                              zmm30 = _mm512_load_ps(&xre[i+240]);
                              _mm512_store_ps(&zre[i+240], _mm512_mul_ps(zmm30,zmmx));
                              zmm31 = _mm512_load_ps(&xim[i+240]);
                              _mm512_storeu_ps(&zim[i+240], _mm512_mul_ps(zmm31,zmmx)); 

                         }

                        for(; (i+191) < n; i += 192) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                              zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                              zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));
                              zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_mul_ps(zmm8,zmmx));
                              zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_mul_ps(zmm9,zmmx));
                              zmm10 = _mm512_load_ps(&xre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_mul_ps(zmm10,zmmx));
                              zmm11 = _mm512_load_ps(&xim[i+80]);
                              _mm512_store_ps(&zim[i+80], _mm512_mul_ps(zmm11,zmmx));
                              zmm12 = _mm512_load_ps(&xre[i+96]);
                              _mm512_store_ps(&zre[i+96], _mm512_mul_ps(zmm12,zmmx));
                              zmm13 = _mm512_load_ps(&xim[i+96]);
                              _mm512_store_ps(&zim[i+96], _mm512_mul_ps(zmm13,zmmx));
                              zmm14 = _mm512_load_ps(&xre[i+112]);
                              _mm512_store_ps(&zre[i+112], _mm512_mul_ps(zmm14,zmmx));
                              zmm15 = _mm512_load_ps(&xim[i+112]);
                              _mm512_store_ps(&zim[i+112], _mm512_mul_ps(zmm15,zmmx));
                              zmm16 = _mm512_load_ps(&xre[i+128]);
                              _mm512_store_ps(&zre[i+128], _mm512_mul_ps(zmm16,zmmx));
                              zmm17 = _mm512_load_ps(&xim[i+128]);
                              _mm512_store_ps(&zim[i+128], _mm512_mul_ps(zmm17,zmmx));
                              zmm18 = _mm512_load_ps(&xre[i+144]);
                              _mm512_store_ps(&zre[i+144], _mm512_mul_ps(zmm18,zmmx));
                              zmm19 = _mm512_load_ps(&xim[i+144]);
                              _mm512_store_ps(&zim[i+144], _mm512_mul_ps(zmm19,zmmx));
                              zmm20 = _mm512_load_ps(&xre[i+160]);
                              _mm512_store_ps(&zre[i+160], _mm512_mul_ps(zmm20,zmmx));
                              zmm21 = _mm512_load_ps(&xim[i+160]);
                              _mm512_store_ps(&zim[i+160], _mm512_mul_ps(zmm21,zmmx));
                              zmm22 = _mm512_load_ps(&xre[i+176]);
                              _mm512_store_ps(&zre[i+176], _mm512_mul_ps(zmm22,zmmx));
                              zmm23 = _mm512_load_ps(&xim[i+176]);
                              _mm512_store_ps(&zim[i+176], _mm512_mul_ps(zmm23,zmmx)); 

                        }

                        for(; (i+127) < n; i += 128) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                              zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                              zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));
                              zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_mul_ps(zmm8,zmmx));
                              zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_mul_ps(zmm9,zmmx));
                              zmm10 = _mm512_load_ps(&xre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_mul_ps(zmm10,zmmx));
                              zmm11 = _mm512_load_ps(&xim[i+80]);
                              _mm512_store_ps(&zim[i+80], _mm512_mul_ps(zmm11,zmmx));
                              zmm12 = _mm512_load_ps(&xre[i+96]);
                              _mm512_store_ps(&zre[i+96], _mm512_mul_ps(zmm12,zmmx));
                              zmm13 = _mm512_load_ps(&xim[i+96]);
                              _mm512_store_ps(&zim[i+96], _mm512_mul_ps(zmm13,zmmx));
                              zmm14 = _mm512_load_ps(&xre[i+112]);
                              _mm512_store_ps(&zre[i+112], _mm512_mul_ps(zmm14,zmmx));
                              zmm15 = _mm512_load_ps(&xim[i+112]);
                              _mm512_store_ps(&zim[i+112], _mm512_mul_ps(zmm15,zmmx));
                        }

                       for(; (i+63) < n; i += 64) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                              zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                              zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));
                       }

                      for(; (i+31) < n; i += 32) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                      }

                     for(; (i+15) < n; i += 16) {
                             zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                     }
                    
                    
                    for(; (i+0) < n; i += 1) {
                          const float s = pzmm[i];
                          zre[i] = xre[i] * s;
                          zim[i] = zim[i] * s;
                    }
                 }


                   
                   void gms::math::cmul_zmm16r4_unroll_10x_u(const float * __restrict xre,
                                                  const float * __restrict xim, 
                                                  const __m512             vs,
                                                  float  * __restrict zre,
                                                  float  * __restrict zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                           __m512 zmm0,zmm1,zmm2,zmm3;
                           __m512 zmm4,zmm5,zmm6,zmm7;
                           __m512 zmm8,zmm9,zmm10,zmm11;
                           __m512 zmm12,zmm13,zmm14,zmm15;
                           __m512 zmm16,zmm17,zmm18,zmm19;
                            __m512 zmmx    = vs;
                          const float * __restrict pzmm =  (const float*)&vs[0];
                         for(i = 0; (i+159) < n; i += 160) {
                              _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+32],_MM_HINT_T0);
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                              zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                              zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                              zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                              _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                              zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                              zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));
                              zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_mul_ps(zmm8,zmmx));
                              zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_mul_ps(zmm9,zmmx));
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                              zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_mul_ps(zmm10,zmmx));
                              zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_mul_ps(zmm11,zmmx));
                              zmm12 = _mm512_loadu_ps(&xre[i+96]);
                              _mm512_storeu_ps(&zre[i+96], _mm512_mul_ps(zmm12,zmmx));
                              zmm13 = _mm512_loadu_ps(&xim[i+96]);
                              _mm512_storeu_ps(&zim[i+96], _mm512_mul_ps(zmm13,zmmx));
                              _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                              zmm14 = _mm512_loadu_ps(&xre[i+112]);
                              _mm512_storeu_ps(&zre[i+112], _mm512_mul_ps(zmm14,zmmx));
                              zmm15 = _mm512_loadu_ps(&xim[i+112]);
                              _mm512_storeu_ps(&zim[i+112], _mm512_mul_ps(zmm15,zmmx));
                              zmm16 = _mm512_loadu_ps(&xre[i+128]);
                              _mm512_storeu_ps(&zre[i+128], _mm512_mul_ps(zmm16,zmmx));
                              zmm17 = _mm512_loadu_ps(&xim[i+128]);
                              _mm512_storeu_ps(&zim[i+128], _mm512_mul_ps(zmm17,zmmx));
                              zmm18 = _mm512_loadu_ps(&xre[i+144]);
                              _mm512_storeu_ps(&zre[i+144], _mm512_mul_ps(zmm18,zmmx));
                              zmm19 = _mm512_loadu_ps(&xim[i+144]);
                              _mm512_storeu_ps(&zim[i+144], _mm512_mul_ps(zmm19,zmmx));
                        }

                       for(; (i+127) < n; i += 128) {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                              zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                              zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                              zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                              zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                              zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));
                              zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_mul_ps(zmm8,zmmx));
                              zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_mul_ps(zmm9,zmmx));
                              zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_mul_ps(zmm10,zmmx));
                              zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_mul_ps(zmm11,zmmx));
                              zmm12 = _mm512_loadu_ps(&xre[i+96]);
                              _mm512_storeu_ps(&zre[i+96], _mm512_mul_ps(zmm12,zmmx));
                              zmm13 = _mm512_loadu_ps(&xim[i+96]);
                              _mm512_storeu_ps(&zim[i+96], _mm512_mul_ps(zmm13,zmmx));
                              zmm14 = _mm512_loadu_ps(&xre[i+112]);
                              _mm512_storeu_ps(&zre[i+112], _mm512_mul_ps(zmm14,zmmx));
                              zmm15 = _mm512_loadu_ps(&xim[i+112]);
                              _mm512_storeu_ps(&zim[i+112], _mm512_mul_ps(zmm15,zmmx));
                       }

                       for(; (i+79) < n; i += 80) {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                              zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                              zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                              zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                              zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                              zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));
                              zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_mul_ps(zmm8,zmmx));
                              zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_mul_ps(zmm9,zmmx));
                       }

                      for(; (i+63) < n; i += 64) {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                              zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                              zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                              zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                              zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                              zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx)); 
                      }

                     for(; (i+31) < n; i += 32) {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                              zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                     }

                    for(; (i+15) < n; i += 16) {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));    
                    }

                    for(; (i+0) < n; i += 1) {
                          const float s = pzmm[i];
                          zre[i] = xre[i] * s;
                          zim[i] = zim[i] * s;
                    }

                 }


                  
                   void gms::math::cmul_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim, 
                                                  const __m512             vs,
                                                  float  * __restrict __ATTR_ALIGN__(64) zre,
                                                  float  * __restrict __ATTR_ALIGN__(64) zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                           __m512 zmm0,zmm1,zmm2,zmm3;
                           __m512 zmm4,zmm5,zmm6,zmm7;
                           __m512 zmm8,zmm9,zmm10,zmm11;
                           __m512 zmm12,zmm13,zmm14,zmm15;
                           __m512 zmm16,zmm17,zmm18,zmm19;
                            __m512 zmmx    = vs;
                          const float * __restrict pzmm =  (const float*)&vs[0];
                          for(i = 0; (i+159) < n; i += 160) {
                              _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+32],_MM_HINT_T0);
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                              zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                              _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                               zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));
                              zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_mul_ps(zmm8,zmmx));
                              zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_mul_ps(zmm9,zmmx));
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                               zmm10 = _mm512_load_ps(&xre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_mul_ps(zmm10,zmmx));
                              zmm11 = _mm512_load_ps(&xim[i+80]);
                              _mm512_store_ps(&zim[i+80], _mm512_mul_ps(zmm11,zmmx));
                               zmm12 = _mm512_load_ps(&xre[i+96]);
                              _mm512_store_ps(&zre[i+96], _mm512_mul_ps(zmm12,zmmx));
                              zmm13 = _mm512_load_ps(&xim[i+96]);
                              _mm512_storeu_ps(&zim[i+96], _mm512_mul_ps(zmm13,zmmx));
                              _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                              zmm14 = _mm512_load_ps(&xre[i+112]);
                              _mm512_store_ps(&zre[i+112], _mm512_mul_ps(zmm14,zmmx));
                              zmm15 = _mm512_load_ps(&xim[i+112]);
                              _mm512_store_ps(&zim[i+112], _mm512_mul_ps(zmm15,zmmx));
                               zmm16 = _mm512_load_ps(&xre[i+128]);
                              _mm512_store_ps(&zre[i+128], _mm512_mul_ps(zmm16,zmmx));
                               zmm17 = _mm512_load_ps(&xim[i+128]);
                              _mm512_store_ps(&zim[i+128], _mm512_mul_ps(zmm17,zmmx));
                               zmm18 = _mm512_load_ps(&xre[i+144]);
                              _mm512_store_ps(&zre[i+144], _mm512_mul_ps(zmm18,zmmx));
                               zmm19 = _mm512_load_ps(&xim[i+144]);
                              _mm512_store_ps(&zim[i+144], _mm512_mul_ps(zmm19,zmmx));
                        }

                       for(; (i+127) < n; i += 128) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                              zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                               zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));
                              zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_mul_ps(zmm8,zmmx));
                              zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_mul_ps(zmm9,zmmx));
                              zmm10 = _mm512_load_ps(&xre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_mul_ps(zmm10,zmmx));
                              zmm11 = _mm512_load_ps(&xim[i+80]);
                              _mm512_store_ps(&zim[i+80], _mm512_mul_ps(zmm11,zmmx));
                               zmm12 = _mm512_load_ps(&xre[i+96]);
                              _mm512_store_ps(&zre[i+96], _mm512_mul_ps(zmm12,zmmx));
                              zmm13 = _mm512_load_ps(&xim[i+96]);
                              _mm512_storeu_ps(&zim[i+96], _mm512_mul_ps(zmm13,zmmx));
                              zmm14 = _mm512_load_ps(&xre[i+112]);
                              _mm512_store_ps(&zre[i+112], _mm512_mul_ps(zmm14,zmmx));
                              zmm15 = _mm512_load_ps(&xim[i+112]);
                              _mm512_store_ps(&zim[i+112], _mm512_mul_ps(zmm15,zmmx));
                       }

                       for(; (i+79) < n; i += 80) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                              zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                               zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));
                              zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_mul_ps(zmm8,zmmx));
                              zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_mul_ps(zmm9,zmmx));
                       }

                      for(; (i+63) < n; i += 64) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                              zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                               zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));
                      }

                     for(; (i+31) < n; i += 32) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                     }

                    for(; (i+15) < n; i += 16) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                    }

                    for(; (i+0) < n; i += 1) {
                          const float s = pzmm[i];
                          zre[i] = xre[i] * s;
                          zim[i] = zim[i] * s;
                    }

                 }


                
                   void gms::math::cmul_zmm16r4_unroll_6x_u(const float * __restrict xre,
                                                  const float * __restrict xim, 
                                                  const __m512             vs,
                                                  float  * __restrict zre,
                                                  float  * __restrict zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                           __m512 zmm0,zmm1,zmm2,zmm3;
                           __m512 zmm4,zmm5,zmm6,zmm7;
                           __m512 zmm8,zmm9,zmm10,zmm11;
                            __m512 zmmx    = vs;
                          const float * __restrict pzmm =  (const float*)&vs[0];
                          for(i = 0; (i+95) < n; i += 96) {
                              _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                               zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                               zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                              zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                               zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                              zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                              zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));
                              zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_mul_ps(zmm8,zmmx));
                              zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_mul_ps(zmm9,zmmx));
                              zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_mul_ps(zmm10,zmmx));
                              zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_mul_ps(zmm11,zmmx));
                            
                        }

                       for(; (i+63) < n; i += 64) {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                               zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                               zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                              zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                               zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                              zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                              zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));
                      }

                     for(; (i+31) < n; i += 32) {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                               zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                               zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                              zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));  
                     }

                    for(; (i+15) < n; i += 16) {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                               zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                    }

                    for(; (i+0) < n; i += 1) {
                          const float s = pzmm[i];
                          zre[i] = xre[i] * s;
                          zim[i] = zim[i] * s;
                    }

                 }


                  
                   void gms::math::cmul_zmm16r4_unroll_6x_a( const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim, 
                                                  const __m512             vs,
                                                  float  * __restrict __ATTR_ALIGN__(64) zre,
                                                  float  * __restrict __ATTR_ALIGN__(64) zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                           __m512 zmm0,zmm1,zmm2,zmm3;
                           __m512 zmm4,zmm5,zmm6,zmm7;
                           __m512 zmm8,zmm9,zmm10,zmm11;
                            __m512 zmmx    = vs;
                          const float * __restrict pzmm =  (const float*)&vs[0];
                          for(i = 0; (i+95) < n; i += 96) {
                              _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                               zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                              zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                              zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));
                               zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_mul_ps(zmm8,zmmx));
                               zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_mul_ps(zmm9,zmmx));
                              zmm10 = _mm512_load_ps(&xre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_mul_ps(zmm10,zmmx));
                               zmm11 = _mm512_load_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_mul_ps(zmm11,zmmx));
                        }

                        for(; (i+63) < n; i += 64) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                               zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm4,zmmx));
                              zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm5,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm6,zmmx));
                              zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_mul_ps(zmm7,zmmx));
                      }

                     for(; (i+31) < n; i += 32) {
                             zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                             zmm1 = _mm512_load_ps(&xim[i+0]);
                             _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                             zmm2 = _mm512_load_ps(&xre[i+16]);
                             _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm2,zmmx));
                             zmm3 = _mm512_load_ps(&xim[i+16]);
                             _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm3,zmmx));
                     }

                    for(; (i+15) < n; i += 16) {
                             zmm0 = _mm512_load_ps(&xre[i+0]);
                             _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmmx));
                             zmm1 = _mm512_load_ps(&xim[i+0]);
                             _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm1,zmmx));
                    }

                    for(; (i+0) < n; i += 1) {
                          const float s = pzmm[i];
                          zre[i] = xre[i] * s;
                          zim[i] = zim[i] * s;
                    }

                 }


                 
                   void gms::math::cmul_zmm16r4_unroll_16x_u(const float * __restrict  xre,
                                                  const float * __restrict  xim,
                                                  const float * __restrict  yre,
                                                  float * __restrict        zre,
                                                  float * __restrict        zim,
                                                  const int32_t n) {

                        if(__builtin_expect(0==n,0)) {return;}
                         __m512 zmm0,zmm1,zmm2,zmm3;
                         __m512 zmm4,zmm5,zmm6,zmm7;
                         __m512 zmm8,zmm9,zmm10,zmm11;
                         __m512 zmm12,zmm13,zmm14,zmm15;
                         __m512 zmm16,zmm17,zmm18,zmm19;
                         __m512 zmm20,zmm21,zmm22,zmm23;
                         __m512 zmm24,zmm25,zmm26,zmm27;
                         __m512 zmm28,zmm29,zmm30,zmm31;
                        int32_t i;
                        
                        for(i = 0; (i+255) < n; i += 256) {
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                            _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_loadu_ps(&xre[i+16]);
                            zmm4  = _mm512_loadu_ps(&yre[i+16]); //yre
                            _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_loadu_ps(&xim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                            zmm6  = _mm512_loadu_ps(&xre[i+32]);
                            zmm7  = _mm512_loadu_ps(&yre[i+32]); //yre
                            _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                            zmm9  = _mm512_loadu_ps(&xre[i+48]);
                            zmm10 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                            zmm11 = _mm512_loadu_ps(&xim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                            zmm12  = _mm512_loadu_ps(&xre[i+64]);
                            zmm13  = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_mul_ps(zmm12,zmm13));
                            zmm14  = _mm512_loadu_ps(&xim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_mul_ps(zmm14,zmm13));
                            zmm15  = _mm512_loadu_ps(&xre[i+80]);
                            zmm16  = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_mul_ps(zmm15,zmm16));
                            zmm17  = _mm512_loadu_ps(&xim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_mul_ps(zmm17,zmm16));
                            zmm18  = _mm512_loadu_ps(&xre[i+96]);
                            zmm19  = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_mul_ps(zmm18,zmm19));
                            zmm20  = _mm512_loadu_ps(&xim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_mul_ps(zmm20,zmm19));
                            zmm21  = _mm512_loadu_ps(&xre[i+112]);
                            zmm22  = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_mul_ps(zmm21,zmm22));
                            zmm23  = _mm512_loadu_ps(&xim[i+112]);
                            _mm512_storeu_ps(&zim[i+112], _mm512_mul_ps(zmm23,zmm22));
                            _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                            zmm24  = _mm512_loadu_ps(&xre[i+128]);
                            zmm25  = _mm512_loadu_ps(&yre[i+128]);
                            _mm512_storeu_ps(&zre[i+128], _mm512_mul_ps(zmm25,zmm24));
                            zmm26  = _mm512_loadu_ps(&xim[i+128]);
                            _mm512_storeu_ps(&zim[i+128], _mm512_mul_ps(zmm26,zmm24));
                            zmm27  = _mm512_loadu_ps(&xre[i+144]);
                            zmm28  = _mm512_loadu_ps(&yre[i+144]);
                            _mm512_storeu_ps(&zre[i+144], _mm512_mul_ps(zmm27,zmm28));
                            zmm29  = _mm512_loadu_ps(&xim[i+144]);
                            _mm512_storeu_ps(&zim[i+144], _mm512_mul_ps(zmm29,zmm28));
                            zmm30  = _mm512_loadu_ps(&xre[i+160]);
                            zmm31  = _mm512_loadu_ps(&yre[i+160]);
                            _mm512_storeu_ps(&zre[i+160], _mm512_mul_ps(zmm30,zmm31));
                            zmm1  = _mm512_loadu_ps(&xim[i+160]);
                            _mm512_storeu_ps(&zim[i+160], _mm512_mul_ps(zmm1,zmm31));
                            zmm2  = _mm512_loadu_ps(&xre[i+176]);
                            zmm3  = _mm512_loadu_ps(&yre[i+176]);
                            _mm512_storeu_ps(&zre[i+176], _mm512_mul_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xim[i+176]);
                            _mm512_storeu_ps(&zim[i+176], _mm512_mul_ps(zmm4,zmm3));
                            _mm_prefetch((const char *)&xre[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+255],_MM_HINT_T0); 
                            zmm5  = _mm512_loadu_ps(&xre[i+192]);
                            zmm6  = _mm512_loadu_ps(&yre[i+192]);
                            _mm512_storeu_ps(&zre[i+192], _mm512_mul_ps(zmm5,zmm6));
                            zmm7  = _mm512_loadu_ps(&xim[i+192]);
                            _mm512_storeu_ps(&zim[i+192], _mm512_mul_ps(zmm7,zmm6));
                            zmm8  = _mm512_loadu_ps(&xre[i+208]);
                            zmm9  = _mm512_loadu_ps(&yre[i+208]);
                            _mm512_storeu_ps(&zre[i+208], _mm512_mul_ps(zmm8,zmm9));
                            zmm10  = _mm512_loadu_ps(&xim[i+208]);
                            _mm512_storeu_ps(&zim[i+208], _mm512_mul_ps(zmm10,zmm9));
                            zmm11 = _mm512_loadu_ps(&xre[i+224]);
                            zmm12 = _mm512_loadu_ps(&yre[i+224]);
                            _mm512_storeu_ps(&zre[i+224], _mm512_mul_ps(zmm11,zmm12));
                            zmm13  = _mm512_loadu_ps(&xim[i+224]);
                            _mm512_storeu_ps(&zim[i+224], _mm512_mul_ps(zmm13,zmm12));
                            zmm14  = _mm512_loadu_ps(&xre[i+240]);
                            zmm15  = _mm512_loadu_ps(&yre[i+240]);
                            _mm512_storeu_ps(&zre[i+240], _mm512_mul_ps(zmm14,zmm15));
                            zmm16  = _mm512_loadu_ps(&xim[i+240]);
                            _mm512_storeu_ps(&zim[i+240], _mm512_mul_ps(zmm16,zmm15));
                        }

                       for(; (i+191) < n; i += 192) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                            _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_loadu_ps(&xre[i+16]);
                            zmm4  = _mm512_loadu_ps(&yre[i+16]); //yre
                            _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_loadu_ps(&xim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                            zmm6  = _mm512_loadu_ps(&xre[i+32]);
                            zmm7  = _mm512_loadu_ps(&yre[i+32]); //yre
                            _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                            zmm9  = _mm512_loadu_ps(&xre[i+48]);
                            zmm10 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                            zmm11 = _mm512_loadu_ps(&xim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));
                            zmm12  = _mm512_loadu_ps(&xre[i+64]);
                            zmm13  = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_mul_ps(zmm12,zmm13));
                            zmm14  = _mm512_loadu_ps(&xim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_mul_ps(zmm14,zmm13));
                            zmm15  = _mm512_loadu_ps(&xre[i+80]);
                            zmm16  = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_mul_ps(zmm15,zmm16));
                            zmm17  = _mm512_loadu_ps(&xim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_mul_ps(zmm17,zmm16));
                            zmm18  = _mm512_loadu_ps(&xre[i+96]);
                            zmm19  = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_mul_ps(zmm18,zmm19));
                            zmm20  = _mm512_loadu_ps(&xim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_mul_ps(zmm20,zmm19));
                            zmm21  = _mm512_loadu_ps(&xre[i+112]);
                            zmm22  = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_mul_ps(zmm21,zmm22));
                            zmm23  = _mm512_loadu_ps(&xim[i+112]);
                            _mm512_storeu_ps(&zim[i+112], _mm512_mul_ps(zmm23,zmm22));
                            zmm24  = _mm512_loadu_ps(&xre[i+128]);
                            zmm25  = _mm512_loadu_ps(&yre[i+128]);
                            _mm512_storeu_ps(&zre[i+128], _mm512_mul_ps(zmm25,zmm24));
                            zmm26  = _mm512_loadu_ps(&xim[i+128]);
                            _mm512_storeu_ps(&zim[i+128], _mm512_mul_ps(zmm26,zmm24));
                            zmm27  = _mm512_loadu_ps(&xre[i+144]);
                            zmm28  = _mm512_loadu_ps(&yre[i+144]);
                            _mm512_storeu_ps(&zre[i+144], _mm512_mul_ps(zmm27,zmm28));
                            zmm29  = _mm512_loadu_ps(&xim[i+144]);
                            _mm512_storeu_ps(&zim[i+144], _mm512_mul_ps(zmm29,zmm28));
                            zmm30  = _mm512_loadu_ps(&xre[i+160]);
                            zmm31  = _mm512_loadu_ps(&yre[i+160]);
                            _mm512_storeu_ps(&zre[i+160], _mm512_mul_ps(zmm30,zmm31));
                            zmm1  = _mm512_loadu_ps(&xim[i+160]);
                            _mm512_storeu_ps(&zim[i+160], _mm512_mul_ps(zmm1,zmm31));
                            zmm2  = _mm512_loadu_ps(&xre[i+176]);
                            zmm3  = _mm512_loadu_ps(&yre[i+176]);
                            _mm512_storeu_ps(&zre[i+176], _mm512_mul_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xim[i+176]);
                            _mm512_storeu_ps(&zim[i+176], _mm512_mul_ps(zmm4,zmm3));
                       }

                       for(; (i+127) < n; i += 128) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                            _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_loadu_ps(&xre[i+16]);
                            zmm4  = _mm512_loadu_ps(&yre[i+16]); //yre
                            _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_loadu_ps(&xim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                            zmm6  = _mm512_loadu_ps(&xre[i+32]);
                            zmm7  = _mm512_loadu_ps(&yre[i+32]); //yre
                            _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                            zmm9  = _mm512_loadu_ps(&xre[i+48]);
                            zmm10 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                            zmm11 = _mm512_loadu_ps(&xim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));
                            zmm12  = _mm512_loadu_ps(&xre[i+64]);
                            zmm13  = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_mul_ps(zmm12,zmm13));
                            zmm14  = _mm512_loadu_ps(&xim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_mul_ps(zmm14,zmm13));
                            zmm15  = _mm512_loadu_ps(&xre[i+80]);
                            zmm16  = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_mul_ps(zmm15,zmm16));
                            zmm17  = _mm512_loadu_ps(&xim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_mul_ps(zmm17,zmm16));
                            zmm18  = _mm512_loadu_ps(&xre[i+96]);
                            zmm19  = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_mul_ps(zmm18,zmm19));
                            zmm20  = _mm512_loadu_ps(&xim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_mul_ps(zmm20,zmm19));
                            zmm21  = _mm512_loadu_ps(&xre[i+112]);
                            zmm22  = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_mul_ps(zmm21,zmm22));
                            zmm23  = _mm512_loadu_ps(&xim[i+112]);
                            _mm512_storeu_ps(&zim[i+112], _mm512_mul_ps(zmm23,zmm22));
                       }

                       for(; (i+63) < n; i += 64) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                            _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_loadu_ps(&xre[i+16]);
                            zmm4  = _mm512_loadu_ps(&yre[i+16]); //yre
                            _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_loadu_ps(&xim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                            zmm6  = _mm512_loadu_ps(&xre[i+32]);
                            zmm7  = _mm512_loadu_ps(&yre[i+32]); //yre
                            _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                            zmm9  = _mm512_loadu_ps(&xre[i+48]);
                            zmm10 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                            zmm11 = _mm512_loadu_ps(&xim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));

                       }

                       for(; (i+31) < n; i += 32) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                            _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_loadu_ps(&xre[i+16]);
                            zmm4  = _mm512_loadu_ps(&yre[i+16]); //yre
                            _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_loadu_ps(&xim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));

                       }

                      for(; (i+15) < n; i += 16) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                            _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));

                      }

                      for(; (i+0) < n; i += 1) {
                            zre[i] = xre[i] * yre[i];
                            zim[i] = xim[i] * yre[i];
                      }
                  }


                 
                   void cgms::math::mul_zmm16r4_unroll_16x_a(const float * __restrict  __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict  __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict  __ATTR_ALIGN__(64) yre,
                                                  float * __restrict        __ATTR_ALIGN__(64) zre,
                                                  float * __restrict        __ATTR_ALIGN__(64) zim,
                                                  const int32_t n) {

                        if(__builtin_expect(0==n,0)) {return;}
                         __m512 zmm0,zmm1,zmm2,zmm3;
                         __m512 zmm4,zmm5,zmm6,zmm7;
                         __m512 zmm8,zmm9,zmm10,zmm11;
                         __m512 zmm12,zmm13,zmm14,zmm15;
                         __m512 zmm16,zmm17,zmm18,zmm19;
                         __m512 zmm20,zmm21,zmm22,zmm23;
                         __m512 zmm24,zmm25,zmm26,zmm27;
                         __m512 zmm28,zmm29,zmm30,zmm31;
                        int32_t i;
                        
                        for(i = 0; (i+255) < n; i += 256) {
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                            _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_load_ps(&xre[i+16]);
                            zmm4  = _mm512_load_ps(&yre[i+16]); //yre
                            _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_load_ps(&xim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                            zmm6  = _mm512_load_ps(&xre[i+32]);
                            zmm7  = _mm512_load_ps(&yre[i+32]); //yre
                            _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                            zmm9  = _mm512_load_ps(&xre[i+48]);
                            zmm10 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                            zmm11 = _mm512_load_ps(&xim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                            zmm12  = _mm512_load_ps(&xre[i+64]);
                            zmm13  = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_mul_ps(zmm12,zmm13));
                            zmm14  = _mm512_loadu_ps(&xim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_mul_ps(zmm14,zmm13));
                            zmm15  = _mm512_load_ps(&xre[i+80]);
                            zmm16  = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_mul_ps(zmm15,zmm16));
                            zmm17  = _mm512_load_ps(&xim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_mul_ps(zmm17,zmm16));
                            zmm18  = _mm512_load_ps(&xre[i+96]);
                            zmm19  = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_mul_ps(zmm18,zmm19));
                            zmm20  = _mm512_load_ps(&xim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_mul_ps(zmm20,zmm19));
                            zmm21  = _mm512_load_ps(&xre[i+112]);
                            zmm22  = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_mul_ps(zmm21,zmm22));
                            zmm23  = _mm512_loadu_ps(&xim[i+112]);
                            _mm512_store_ps(&zim[i+112], _mm512_mul_ps(zmm23,zmm22));
                            _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                            zmm24  = _mm512_load_ps(&xre[i+128]);
                            zmm25  = _mm512_load_ps(&yre[i+128]);
                            _mm512_store_ps(&zre[i+128], _mm512_mul_ps(zmm25,zmm24));
                            zmm26  = _mm512_load_ps(&xim[i+128]);
                            _mm512_store_ps(&zim[i+128], _mm512_mul_ps(zmm26,zmm24));
                            zmm27  = _mm512_load_ps(&xre[i+144]);
                            zmm28  = _mm512_load_ps(&yre[i+144]);
                            _mm512_store_ps(&zre[i+144], _mm512_mul_ps(zmm27,zmm28));
                            zmm29  = _mm512_load_ps(&xim[i+144]);
                            _mm512_store_ps(&zim[i+144], _mm512_mul_ps(zmm29,zmm28));
                            zmm30  = _mm512_load_ps(&xre[i+160]);
                            zmm31  = _mm512_load_ps(&yre[i+160]);
                            _mm512_store_ps(&zre[i+160], _mm512_mul_ps(zmm30,zmm31));
                            zmm1  = _mm512_load_ps(&xim[i+160]);
                            _mm512_store_ps(&zim[i+160], _mm512_mul_ps(zmm1,zmm31));
                            zmm2  = _mm512_load_ps(&xre[i+176]);
                            zmm3  = _mm512_load_ps(&yre[i+176]);
                            _mm512_store_ps(&zre[i+176], _mm512_mul_ps(zmm2,zmm3));
                            zmm4  = _mm512_loau_ps(&xim[i+176]);
                            _mm512_store_ps(&zim[i+176], _mm512_mul_ps(zmm4,zmm3));
                            _mm_prefetch((const char *)&xre[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+255],_MM_HINT_T0); 
                            zmm5  = _mm512_load_ps(&xre[i+192]);
                            zmm6  = _mm512_load_ps(&yre[i+192]);
                            _mm512_store_ps(&zre[i+192], _mm512_mul_ps(zmm5,zmm6));
                            zmm7  = _mm512_load_ps(&xim[i+192]);
                            _mm512_store_ps(&zim[i+192], _mm512_mul_ps(zmm7,zmm6));
                            zmm8  = _mm512_load_ps(&xre[i+208]);
                            zmm9  = _mm512_load_ps(&yre[i+208]);
                            _mm512_store_ps(&zre[i+208], _mm512_mul_ps(zmm8,zmm9));
                            zmm10  = _mm512_load_ps(&xim[i+208]);
                            _mm512_store_ps(&zim[i+208], _mm512_mul_ps(zmm10,zmm9));
                            zmm11 = _mm512_load_ps(&xre[i+224]);
                            zmm12 = _mm512_load_ps(&yre[i+224]);
                            _mm512_store_ps(&zre[i+224], _mm512_mul_ps(zmm11,zmm12));
                            zmm13  = _mm512_load_ps(&xim[i+224]);
                            _mm512_store_ps(&zim[i+224], _mm512_mul_ps(zmm13,zmm12));
                            zmm14  = _mm512_load_ps(&xre[i+240]);
                            zmm15  = _mm512_load_ps(&yre[i+240]);
                            _mm512_store_ps(&zre[i+240], _mm512_mul_ps(zmm14,zmm15));
                            zmm16  = _mm512_load_ps(&xim[i+240]);
                            _mm512_store_ps(&zim[i+240], _mm512_mul_ps(zmm16,zmm15));
                        }

                       for(; (i+191) < n; i += 192) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                            _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_load_ps(&xre[i+16]);
                            zmm4  = _mm512_load_ps(&yre[i+16]); //yre
                            _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_load_ps(&xim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                            zmm6  = _mm512_load_ps(&xre[i+32]);
                            zmm7  = _mm512_load_ps(&yre[i+32]); //yre
                            _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                            zmm9  = _mm512_load_ps(&xre[i+48]);
                            zmm10 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                            zmm11 = _mm512_load_ps(&xim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));
                            zmm12  = _mm512_load_ps(&xre[i+64]);
                            zmm13  = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_mul_ps(zmm12,zmm13));
                            zmm14  = _mm512_load_ps(&xim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_mul_ps(zmm14,zmm13));
                            zmm15  = _mm512_load_ps(&xre[i+80]);
                            zmm16  = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_mul_ps(zmm15,zmm16));
                            zmm17  = _mm512_loadu_ps(&xim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_mul_ps(zmm17,zmm16));
                            zmm18  = _mm512_load_ps(&xre[i+96]);
                            zmm19  = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_mul_ps(zmm18,zmm19));
                            zmm20  = _mm512_load_ps(&xim[i+96]);
                            _mm512_store_ps(&zim[i+96], _mm512_mul_ps(zmm20,zmm19));
                            zmm21  = _mm512_load_ps(&xre[i+112]);
                            zmm22  = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_mul_ps(zmm21,zmm22));
                            zmm23  = _mm512_load_ps(&xim[i+112]);
                            _mm512_store_ps(&zim[i+112], _mm512_mul_ps(zmm23,zmm22));
                            zmm24  = _mm512_load_ps(&xre[i+128]);
                            zmm25  = _mm512_load_ps(&yre[i+128]);
                            _mm512_store_ps(&zre[i+128], _mm512_mul_ps(zmm25,zmm24));
                            zmm26  = _mm512_load_ps(&xim[i+128]);
                            _mm512_store_ps(&zim[i+128], _mm512_mul_ps(zmm26,zmm24));
                            zmm27  = _mm512_load_ps(&xre[i+144]);
                            zmm28  = _mm512_load_ps(&yre[i+144]);
                            _mm512_store_ps(&zre[i+144], _mm512_mul_ps(zmm27,zmm28));
                            zmm29  = _mm512_load_ps(&xim[i+144]);
                            _mm512_store_ps(&zim[i+144], _mm512_mul_ps(zmm29,zmm28));
                            zmm30  = _mm512_load_ps(&xre[i+160]);
                            zmm31  = _mm512_load_ps(&yre[i+160]);
                            _mm512_store_ps(&zre[i+160], _mm512_mul_ps(zmm30,zmm31));
                            zmm1  = _mm512_load_ps(&xim[i+160]);
                            _mm512_store_ps(&zim[i+160], _mm512_mul_ps(zmm1,zmm31));
                            zmm2  = _mm512_load_ps(&xre[i+176]);
                            zmm3  = _mm512_load_ps(&yre[i+176]);
                            _mm512_store_ps(&zre[i+176], _mm512_mul_ps(zmm2,zmm3));
                            zmm4  = _mm512_load_ps(&xim[i+176]);
                            _mm512_store_ps(&zim[i+176], _mm512_mul_ps(zmm4,zmm3));
                       }

                       for(; (i+127) < n; i += 128) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                            _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_load_ps(&xre[i+16]);
                            zmm4  = _mm512_load_ps(&yre[i+16]); //yre
                            _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_load_ps(&xim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                            zmm6  = _mm512_load_ps(&xre[i+32]);
                            zmm7  = _mm512_load_ps(&yre[i+32]); //yre
                            _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                            zmm9  = _mm512_load_ps(&xre[i+48]);
                            zmm10 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                            zmm11 = _mm512_load_ps(&xim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));
                            zmm12  = _mm512_load_ps(&xre[i+64]);
                            zmm13  = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_mul_ps(zmm12,zmm13));
                            zmm14  = _mm512_load_ps(&xim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_mul_ps(zmm14,zmm13));
                            zmm15  = _mm512_load_ps(&xre[i+80]);
                            zmm16  = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_mul_ps(zmm15,zmm16));
                            zmm17  = _mm512_load_ps(&xim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_mul_ps(zmm17,zmm16));
                            zmm18  = _mm512_load_ps(&xre[i+96]);
                            zmm19  = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_mul_ps(zmm18,zmm19));
                            zmm20  = _mm512_load_ps(&xim[i+96]);
                            _mm512_store_ps(&zim[i+96], _mm512_mul_ps(zmm20,zmm19));
                            zmm21  = _mm512_load_ps(&xre[i+112]);
                            zmm22  = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_mul_ps(zmm21,zmm22));
                            zmm23  = _mm512_load_ps(&xim[i+112]);
                            _mm512_store_ps(&zim[i+112], _mm512_mul_ps(zmm23,zmm22));
                       }

                       for(; (i+63) < n; i += 64) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                            _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_load_ps(&xre[i+16]);
                            zmm4  = _mm512_load_ps(&yre[i+16]); //yre
                            _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_load_ps(&xim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                            zmm6  = _mm512_load_ps(&xre[i+32]);
                            zmm7  = _mm512_load_ps(&yre[i+32]); //yre
                            _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                            zmm9  = _mm512_load_ps(&xre[i+48]);
                            zmm10 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                            zmm11 = _mm512_load_ps(&xim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));

                       }

                       for(; (i+31) < n; i += 32) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                            _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_load_ps(&xre[i+16]);
                            zmm4  = _mm512_load_ps(&yre[i+16]); //yre
                            _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_load_ps(&xim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));

                       }

                      for(; (i+15) < n; i += 16) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                            _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));

                      }

                      for(; (i+0) < n; i += 1) {
                            zre[i] = xre[i] * yre[i];
                            zim[i] = xim[i] * yre[i];
                      }
                  }


                 
                   void gms::math::cmul_zmm16r4_unroll_10x_u(const float * __restrict  xre,
                                                  const float * __restrict  xim,
                                                  const float * __restrict  yre,
                                                  float * __restrict        zre,
                                                  float * __restrict        zim,
                                                  const int32_t n) {

                        if(__builtin_expect(0==n,0)) {return;}
                         __m512 zmm0,zmm1,zmm2,zmm3;
                         __m512 zmm4,zmm5,zmm6,zmm7;
                         __m512 zmm8,zmm9,zmm10,zmm11;
                         __m512 zmm12,zmm13,zmm14,zmm15;
                         __m512 zmm16,zmm17,zmm18,zmm19;
                         __m512 zmm20,zmm21,zmm22,zmm23;
                         __m512 zmm24,zmm25,zmm26,zmm27;
                         __m512 zmm28,zmm29,zmm30,zmm31;
                        int32_t i;
                        
                        for(i = 0; (i+159) < n; i += 160) {
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                            _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_loadu_ps(&xre[i+16]);
                            zmm4  = _mm512_loadu_ps(&yre[i+16]); //yre
                            _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_loadu_ps(&xim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                            zmm6  = _mm512_loadu_ps(&xre[i+32]);
                            zmm7  = _mm512_loadu_ps(&yre[i+32]); //yre
                            _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            zmm9  = _mm512_loadu_ps(&xre[i+48]);
                            zmm10 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                            zmm11 = _mm512_loadu_ps(&xim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));
                            zmm12  = _mm512_loadu_ps(&xre[i+64]);
                            zmm13  = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_mul_ps(zmm12,zmm13));
                            zmm14  = _mm512_loadu_ps(&xim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_mul_ps(zmm14,zmm13));
                            zmm15  = _mm512_loadu_ps(&xre[i+80]);
                            zmm16  = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_mul_ps(zmm15,zmm16));
                            zmm17  = _mm512_loadu_ps(&xim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_mul_ps(zmm17,zmm16));
                            _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+144],_MM_HINT_T0);
                            zmm18  = _mm512_loadu_ps(&xre[i+96]);
                            zmm19  = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_mul_ps(zmm18,zmm19));
                            zmm20  = _mm512_loadu_ps(&xim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_mul_ps(zmm20,zmm19));
                            zmm21  = _mm512_loadu_ps(&xre[i+112]);
                            zmm22  = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_mul_ps(zmm21,zmm22));
                            zmm23  = _mm512_loadu_ps(&xim[i+112]);
                            _mm512_storeu_ps(&zim[i+112], _mm512_mul_ps(zmm23,zmm22));
                            zmm24  = _mm512_loadu_ps(&xre[i+128]);
                            zmm25  = _mm512_loadu_ps(&yre[i+128]);
                            _mm512_storeu_ps(&zre[i+128], _mm512_mul_ps(zmm25,zmm24));
                            zmm26  = _mm512_loadu_ps(&xim[i+128]);
                            _mm512_storeu_ps(&zim[i+128], _mm512_mul_ps(zmm26,zmm24));
                            zmm27  = _mm512_loadu_ps(&xre[i+144]);
                            zmm28  = _mm512_loadu_ps(&yre[i+144]);
                            _mm512_storeu_ps(&zre[i+144], _mm512_mul_ps(zmm27,zmm28));
                            zmm29  = _mm512_loadu_ps(&xim[i+144]);
                            _mm512_storeu_ps(&zim[i+144], _mm512_mul_ps(zmm29,zmm28));
                          
                        }

                       for(; (i+127) < n; i += 128) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                            _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_loadu_ps(&xre[i+16]);
                            zmm4  = _mm512_loadu_ps(&yre[i+16]); //yre
                            _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_loadu_ps(&xim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                            zmm6  = _mm512_loadu_ps(&xre[i+32]);
                            zmm7  = _mm512_loadu_ps(&yre[i+32]); //yre
                            _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                            zmm9  = _mm512_loadu_ps(&xre[i+48]);
                            zmm10 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                            zmm11 = _mm512_loadu_ps(&xim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));
                            zmm12  = _mm512_loadu_ps(&xre[i+64]);
                            zmm13  = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_mul_ps(zmm12,zmm13));
                            zmm14  = _mm512_loadu_ps(&xim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_mul_ps(zmm14,zmm13));
                            zmm15  = _mm512_loadu_ps(&xre[i+80]);
                            zmm16  = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_mul_ps(zmm15,zmm16));
                            zmm17  = _mm512_loadu_ps(&xim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_mul_ps(zmm17,zmm16));
                            zmm18  = _mm512_loadu_ps(&xre[i+96]);
                            zmm19  = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_mul_ps(zmm18,zmm19));
                            zmm20  = _mm512_loadu_ps(&xim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_mul_ps(zmm20,zmm19));
                            zmm21  = _mm512_loadu_ps(&xre[i+112]);
                            zmm22  = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_mul_ps(zmm21,zmm22));
                            zmm23  = _mm512_loadu_ps(&xim[i+112]);
                            _mm512_storeu_ps(&zim[i+112], _mm512_mul_ps(zmm23,zmm22));
                       }

                       for(; (i+79) < n; i += 80) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                            _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_loadu_ps(&xre[i+16]);
                            zmm4  = _mm512_loadu_ps(&yre[i+16]); //yre
                            _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_loadu_ps(&xim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                            zmm6  = _mm512_loadu_ps(&xre[i+32]);
                            zmm7  = _mm512_loadu_ps(&yre[i+32]); //yre
                            _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                            zmm9  = _mm512_loadu_ps(&xre[i+48]);
                            zmm10 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                            zmm11 = _mm512_loadu_ps(&xim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));
                            zmm12  = _mm512_loadu_ps(&xre[i+64]);
                            zmm13  = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_mul_ps(zmm12,zmm13));
                            zmm14  = _mm512_loadu_ps(&xim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_mul_ps(zmm14,zmm13));
                          
                       }

                       for(; (i+63) < n; i += 64) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                            _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_loadu_ps(&xre[i+16]);
                            zmm4  = _mm512_loadu_ps(&yre[i+16]); //yre
                            _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_loadu_ps(&xim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                            zmm6  = _mm512_loadu_ps(&xre[i+32]);
                            zmm7  = _mm512_loadu_ps(&yre[i+32]); //yre
                            _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                            zmm9  = _mm512_loadu_ps(&xre[i+48]);
                            zmm10 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                            zmm11 = _mm512_loadu_ps(&xim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));

                       }

                       for(; (i+31) < n; i += 32) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                            _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_loadu_ps(&xre[i+16]);
                            zmm4  = _mm512_loadu_ps(&yre[i+16]); //yre
                            _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_loadu_ps(&xim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));

                       }

                      for(; (i+15) < n; i += 16) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                            _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));

                      }

                      for(; (i+0) < n; i += 1) {
                            zre[i] = xre[i] * yre[i];
                            zim[i] = xim[i] * yre[i];
                      }
                  }


                  
                   void gms::math::cmul_zmm16r4_unroll_10x_a(const float * __restrict  __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict  __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict  __ATTR_ALIGN__(64) yre,
                                                  float * __restrict        __ATTR_ALIGN__(64) zre,
                                                  float * __restrict        __ATTR_ALIGN__(64) zim,
                                                  const int32_t n) {

                        if(__builtin_expect(0==n,0)) {return;}
                         __m512 zmm0,zmm1,zmm2,zmm3;
                         __m512 zmm4,zmm5,zmm6,zmm7;
                         __m512 zmm8,zmm9,zmm10,zmm11;
                         __m512 zmm12,zmm13,zmm14,zmm15;
                         __m512 zmm16,zmm17,zmm18,zmm19;
                         __m512 zmm20,zmm21,zmm22,zmm23;
                         __m512 zmm24,zmm25,zmm26,zmm27;
                         __m512 zmm28,zmm29,zmm30,zmm31;
                        int32_t i;
                        
                        for(i = 0; (i+159) < n; i += 160) {
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                            _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_load_ps(&xre[i+16]);
                            zmm4  = _mm512_load_ps(&yre[i+16]); //yre
                            _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_load_ps(&xim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                            zmm6  = _mm512_load_ps(&xre[i+32]);
                            zmm7  = _mm512_load_ps(&yre[i+32]); //yre
                            _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            zmm9  = _mm512_load_ps(&xre[i+48]);
                            zmm10 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                            zmm11 = _mm512_load_ps(&xim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));
                            zmm12  = _mm512_load_ps(&xre[i+64]);
                            zmm13  = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_mul_ps(zmm12,zmm13));
                            zmm14  = _mm512_load_ps(&xim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_mul_ps(zmm14,zmm13));
                            zmm15  = _mm512_load_ps(&xre[i+80]);
                            zmm16  = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_mul_ps(zmm15,zmm16));
                            zmm17  = _mm512_load_ps(&xim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_mul_ps(zmm17,zmm16));
                            _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+144],_MM_HINT_T0);
                            zmm18  = _mm512_load_ps(&xre[i+96]);
                            zmm19  = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_mul_ps(zmm18,zmm19));
                            zmm20  = _mm512_load_ps(&xim[i+96]);
                            _mm512_store_ps(&zim[i+96], _mm512_mul_ps(zmm20,zmm19));
                            zmm21  = _mm512_load_ps(&xre[i+112]);
                            zmm22  = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_mul_ps(zmm21,zmm22));
                            zmm23  = _mm512_load_ps(&xim[i+112]);
                            _mm512_store_ps(&zim[i+112], _mm512_mul_ps(zmm23,zmm22));
                            zmm24  = _mm512_load_ps(&xre[i+128]);
                            zmm25  = _mm512_load_ps(&yre[i+128]);
                            _mm512_store_ps(&zre[i+128], _mm512_mul_ps(zmm25,zmm24));
                            zmm26  = _mm512_load_ps(&xim[i+128]);
                            _mm512_store_ps(&zim[i+128], _mm512_mul_ps(zmm26,zmm24));
                            zmm27  = _mm512_load_ps(&xre[i+144]);
                            zmm28  = _mm512_load_ps(&yre[i+144]);
                            _mm512_store_ps(&zre[i+144], _mm512_mul_ps(zmm27,zmm28));
                            zmm29  = _mm512_load_ps(&xim[i+144]);
                            _mm512_store_ps(&zim[i+144], _mm512_mul_ps(zmm29,zmm28));
                          
                        }

                       for(; (i+127) < n; i += 128) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                            _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_load_ps(&xre[i+16]);
                            zmm4  = _mm512_load_ps(&yre[i+16]); //yre
                            _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_load_ps(&xim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                            zmm6  = _mm512_load_ps(&xre[i+32]);
                            zmm7  = _mm512_load_ps(&yre[i+32]); //yre
                            _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                            zmm9  = _mm512_load_ps(&xre[i+48]);
                            zmm10 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                            zmm11 = _mm512_load_ps(&xim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));
                            zmm12  = _mm512_load_ps(&xre[i+64]);
                            zmm13  = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_mul_ps(zmm12,zmm13));
                            zmm14  = _mm512_load_ps(&xim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_mul_ps(zmm14,zmm13));
                            zmm15  = _mm512_load_ps(&xre[i+80]);
                            zmm16  = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_mul_ps(zmm15,zmm16));
                            zmm17  = _mm512_load_ps(&xim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_mul_ps(zmm17,zmm16));
                            zmm18  = _mm512_load_ps(&xre[i+96]);
                            zmm19  = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_mul_ps(zmm18,zmm19));
                            zmm20  = _mm512_load_ps(&xim[i+96]);
                            _mm512_store_ps(&zim[i+96], _mm512_mul_ps(zmm20,zmm19));
                            zmm21  = _mm512_load_ps(&xre[i+112]);
                            zmm22  = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_mul_ps(zmm21,zmm22));
                            zmm23  = _mm512_load_ps(&xim[i+112]);
                            _mm512_store_ps(&zim[i+112], _mm512_mul_ps(zmm23,zmm22));
                       }

                       for(; (i+79) < n; i += 80) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                            _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_load_ps(&xre[i+16]);
                            zmm4  = _mm512_load_ps(&yre[i+16]); //yre
                            _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_load_ps(&xim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                            zmm6  = _mm512_load_ps(&xre[i+32]);
                            zmm7  = _mm512_load_ps(&yre[i+32]); //yre
                            _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                            zmm9  = _mm512_load_ps(&xre[i+48]);
                            zmm10 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                            zmm11 = _mm512_load_ps(&xim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));
                            zmm12  = _mm512_load_ps(&xre[i+64]);
                            zmm13  = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_mul_ps(zmm12,zmm13));
                            zmm14  = _mm512_load_ps(&xim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_mul_ps(zmm14,zmm13));
                          
                       }

                       for(; (i+63) < n; i += 64) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                            _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_load_ps(&xre[i+16]);
                            zmm4  = _mm512_load_ps(&yre[i+16]); //yre
                            _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_load_ps(&xim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                            zmm6  = _mm512_load_ps(&xre[i+32]);
                            zmm7  = _mm512_load_ps(&yre[i+32]); //yre
                            _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                            zmm9  = _mm512_load_ps(&xre[i+48]);
                            zmm10 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                            zmm11 = _mm512_load_ps(&xim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));

                       }

                       for(; (i+31) < n; i += 32) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                            _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                            zmm3  = _mm512_load_ps(&xre[i+16]);
                            zmm4  = _mm512_load_ps(&yre[i+16]); //yre
                            _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                            zmm5  = _mm512_load_ps(&xim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));

                       }

                      for(; (i+15) < n; i += 16) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                            _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));

                      }

                      for(; (i+0) < n; i += 1) {
                            zre[i] = xre[i] * yre[i];
                            zim[i] = xim[i] * yre[i];
                      }
                  }


                 
                   void gms::math::cmul_zmm16r4_unroll_6x_u( const float * __restrict xre,
                                                  const float * __restrict xim, 
                                                  const float * __restrict yre
                                                  float  * __restrict zre,
                                                  float  * __restrict zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                           __m512 zmm0,zmm1,zmm2,zmm3;
                           __m512 zmm4,zmm5,zmm6,zmm7;
                           __m512 zmm8,zmm9,zmm10,zmm11;
                          for(i = 0; (i+95) < n; i += 96) {
                              _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                              zmm2  = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                              zmm3  = _mm512_loadu_ps(&xre[i+16]);
                              zmm4  = _mm512_loadu_ps(&yre[i+16]); //yre
                              _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                              zmm6  = _mm512_loadu_ps(&xre[i+32]);
                              zmm7  = _mm512_loadu_ps(&yre[i+32]); //yre
                              _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                              zmm8  = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                              zmm9  = _mm512_loadu_ps(&xre[i+48]);
                              zmm10 = _mm512_loadu_ps(&yre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                              zmm11 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));
                              zmm12  = _mm512_loadu_ps(&xre[i+64]);
                              zmm13  = _mm512_loadu_ps(&yre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_mul_ps(zmm12,zmm13));
                              zmm14  = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_mul_ps(zmm14,zmm13));
                              zmm15  = _mm512_loadu_ps(&xre[i+80]);
                              zmm16  = _mm512_loadu_ps(&yre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_mul_ps(zmm15,zmm16));
                              zmm17  = _mm512_loadu_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_mul_ps(zmm17,zmm16));                           
                        }

                       for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                              zmm2  = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                              zmm3  = _mm512_loadu_ps(&xre[i+16]);
                              zmm4  = _mm512_loadu_ps(&yre[i+16]); //yre
                              _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                              zmm6  = _mm512_loadu_ps(&xre[i+32]);
                              zmm7  = _mm512_loadu_ps(&yre[i+32]); //yre
                              _mm512_storeu_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                              zmm8  = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                              zmm9  = _mm512_loadu_ps(&xre[i+48]);
                              zmm10 = _mm512_loadu_ps(&yre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                              zmm11 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));
                      }

                     for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                              zmm2  = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                              zmm3  = _mm512_loadu_ps(&xre[i+16]);
                              zmm4  = _mm512_loadu_ps(&yre[i+16]); //yre
                              _mm512_storeu_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));   
                     }

                    for(; (i+15) < n; i += 16) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                              _mm512_storeu_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                              zmm2  = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));     
                    }

                    for(; (i+0) < n; i += 1) {
                         
                          zre[i] = xre[i] * yre[i];
                          zim[i] = xim[i] * yre[i];
                    }

                 }


                 
                   void gms::math::cmul_zmm16r4_unroll_6x_a( const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim, 
                                                  const float * __restrict __ATTR_ALIGN__(64) yre,
                                                  float  * __restrict __ATTR_ALIGN__(64) zre,
                                                  float  * __restrict __ATTR_ALIGN__(64) zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                           __m512 zmm0,zmm1,zmm2,zmm3;
                           __m512 zmm4,zmm5,zmm6,zmm7;
                           __m512 zmm8,zmm9,zmm10,zmm11;
                          for(i = 0; (i+95) < n; i += 96) {
                              _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                              zmm2  = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                              zmm3  = _mm512_load_ps(&xre[i+16]);
                              zmm4  = _mm512_load_ps(&yre[i+16]); //yre
                              _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                              zmm6  = _mm512_load_ps(&xre[i+32]);
                              zmm7  = _mm512_load_ps(&yre[i+32]); //yre
                              _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                              zmm8  = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                              zmm9  = _mm512_load_ps(&xre[i+48]);
                              zmm10 = _mm512_load_ps(&yre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                              zmm11 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));
                              zmm12  = _mm512_load_ps(&xre[i+64]);
                              zmm13  = _mm512_load_ps(&yre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_mul_ps(zmm12,zmm13));
                              zmm14  = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_mul_ps(zmm14,zmm13));
                              zmm15  = _mm512_load_ps(&xre[i+80]);
                              zmm16  = _mm512_load_ps(&yre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_mul_ps(zmm15,zmm16));
                              zmm17  = _mm512_load_ps(&xim[i+80]);
                              _mm512_store_ps(&zim[i+80], _mm512_mul_ps(zmm17,zmm16));                           
                        }

                       for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                              zmm2  = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                              zmm3  = _mm512_load_ps(&xre[i+16]);
                              zmm4  = _mm512_load_ps(&yre[i+16]); //yre
                              _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));
                              zmm6  = _mm512_load_ps(&xre[i+32]);
                              zmm7  = _mm512_load_ps(&yre[i+32]); //yre
                              _mm512_store_ps(&zre[i+32], _mm512_mul_ps(zmm6,zmm7));
                              zmm8  = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_mul_ps(zmm8,zmm7));
                              zmm9  = _mm512_load_ps(&xre[i+48]);
                              zmm10 = _mm512_load_ps(&yre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_mul_ps(zmm9,zmm10));
                              zmm11 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_mul_ps(zmm11,zmm10));
                      }

                     for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                              zmm2  = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));
                              zmm3  = _mm512_load_ps(&xre[i+16]);
                              zmm4  = _mm512_load_ps(&yre[i+16]); //yre
                              _mm512_store_ps(&zre[i+16], _mm512_mul_ps(zmm3,zmm4));
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_mul_ps(zmm5,zmm4));   
                     }

                    for(; (i+15) < n; i += 16) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                              _mm512_store_ps(&zre[i+0], _mm512_mul_ps(zmm0,zmm1));
                              zmm2  = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_mul_ps(zmm2,zmm1));     
                    }

                    for(; (i+0) < n; i += 1) {
                         
                          zre[i] = xre[i] * yre[i];
                          zim[i] = xim[i] * yre[i];
                    }

                 }


                /*
                    All of above, but in-place
                 */

                 
                   void gms::math::cmul_zmm16r4_12x_uip(const float * __restrict xre,
                                             const float * __restrict yre,
                                             float * __restrict yre,
                                             float * __restrict yim,
                                             const int32_t n) {

                        if(__builtin_expect(0==n,0)) {return;}

                          __m512 zmm0,zmm1,zmm2,zmm3;
                          __m512 zmm4,zmm5,zmm6,zmm7;
                          __m512 zmm8,zmm9,zmm10,zmm11;
                          __m512 zmm12,zmm13,zmm14,zmm15;
                          __m512 zmm16,zmm17,zmm18,zmm19;
                          __m512 zmm20,zmm21,zmm22,zmm23;
                          __m512 zmm24,zmm25,zmm26,zmm27;
                          __m512 zmm28,zmm29,zmm30,zmm31;
                         int32_t i;
                        for(i = 0; (i+191) < n; i += 192) {
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_storeu_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_loadu_ps(&xre[i+32]);
                            zmm9 = _mm512_loadu_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_loadu_ps(&xim[i+32]);
                            zmm11= _mm512_loadu_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_storeu_ps(&yim[i+32], zmm11);
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0);
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&yim[i+48], zmm15);
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]); // yre modified
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]); // yim modified
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm16,zmm17),
                                                                        _mm512_mul_ps(zmm18,zmm19));
                            _mm512_storeu_ps(&yre[i+64], zmm17);
                            zmm19  = _mm512_mul_ps(_mm512_mul_ps(zmm18,zmm17),
                                                                        _mm512_mul_ps(zmm16,zmm19));
                            _mm512_storeu_ps(&yim[i+64], zmm19);
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]); // yre modified
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]); // yim modified
                            zmm21 = _mm512_sub_ps(_mm512_mul_ps(zmm20,zmm21),
                                                       _mm512_mul_ps(zmm22,zmm23));
                            _mm512_storeu_ps(&yre[i+80], zmm21);
                            zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm22,zmm21),
                                                           _mm512_mul_ps(zmm20,zmm23));
                            _mm512_storeu_ps(&yim[i+80], zmm23);
                            _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+144],_MM_HINT_T0);
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&yre[i+96]); // yre modified
                            zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            zmm27 = _mm512_loadu_ps(&yim[i+96]); // yim modified
                            zmm25 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                       _mm512_mul_ps(zmm25,zmm27));
                            _mm512_storeu_ps(&yre[i+96], zmm25);
                             zmm27  = _mm512_mul_ps(_mm512_mul_ps(zmm26,zmm25),
                                                           _mm512_mul_ps(zmm24,zmm27));
                            _mm512_storeu_ps(&yim[i+96], zmm27);
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&yre[i+112]);  // yre modified
                            zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            zmm31 = _mm512_loadu_ps(&yim[i+112]);  // yim modified
                            zmm29 = _mm512_sub_ps(_mm512_mul_ps(zmm28,zmm29),
                                                       _mm512_mul_ps(zmm30,zmm31));
                            _mm512_storeu_ps(&yre[i+112], zmm29);
                            zmm31  = _mm512_mul_ps(_mm512_mul_ps(zmm30,zmm29),
                                                           _mm512_mul_ps(zmm28,zmm31));
                            _mm512_storeu_ps(&yim[i+112], zmm31);
                            zmm1 = _mm512_loadu_ps(&xre[i+128]);
                            zmm2 = _mm512_loadu_ps(&yre[i+128]); // yre modified
                            zmm3 = _mm512_loadu_ps(&xim[i+128]);
                            zmm4 = _mm512_loadu_ps(&yim[i+128]); // yim modified
                            zmm2 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_storeu_ps(&yre[i+128], zmm2);
                            zmm4  = _mm512_mul_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_storeu_ps(&yim[i+128], zmm4);
                            _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+192],_MM_HINT_T0);
                            zmm5 = _mm512_loadu_ps(&xre[i+144]);
                            zmm6 = _mm512_loadu_ps(&yre[i+144]); // yre modified
                            zmm7 = _mm512_loadu_ps(&xim[i+144]);
                            zmm8 = _mm512_loadu_ps(&yim[i+144]); // yim modified
                            zmm6 = _mm512_sub_ps(_mm512_mul_ps(zmm5,zmm6),
                                                       _mm512_mul_ps(zmm7,zmm8));
                            _mm512_storeu_ps(&yre[i+144], zmm6);
                             zmm8  = _mm512_mul_ps(_mm512_mul_ps(zmm7,zmm6),
                                                           _mm512_mul_ps(zmm5,zmm8));
                            _mm512_storeu_ps(&yim[i+144], zmm8);
                            zmm9 = _mm512_loadu_ps(&xre[i+160]);
                            zmm10= _mm512_loadu_ps(&yre[i+160]); // yre modified
                            zmm11= _mm512_loadu_ps(&xim[i+160]);
                            zmm12= _mm512_loadu_ps(&yim[i+160]); // yim modified
                            zmm10= _mm512_sub_ps(_mm512_mul_ps(zmm9,zmm10),
                                                       _mm512_mul_ps(zmm11,zmm12));
                            _mm512_storeu_ps(&yre[i+160], zmm10);
                            zmm12 = _mm512_mul_ps(_mm512_mul_ps(zmm11,zmm10),
                                                           _mm512_mul_ps(zmm9,zmm12));
                            _mm512_storeu_ps(&yim[i+160], zmm12);
                            zmm13 = _mm512_loadu_ps(&xre[i+176]);
                            zmm14 = _mm512_loadu_ps(&yre[i+176]); // yre modified
                            zmm15 = _mm512_loadu_ps(&xim[i+176]);
                            zmm16 = _mm512_loadu_ps(&yim[i+176]);  // yim modified
                            zmm14 = _mm512_sub_ps(_mm512_mul_ps(zmm13,zmm14),
                                                       _mm512_mul_ps(zmm15,zmm16));
                            _mm512_storeu_ps(&yre[i+176], zmm14);
                            zmm16 = _mm512_mul_ps(_mm512_mul_ps(zmm15,zmm14),
                                                           _mm512_mul_ps(zmm13,zmm16));
                            _mm512_storeu_ps(&yim[i+176], zmm16);
                      }

                        for(; (i+159) < n; i += 160) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_storeu_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_loadu_ps(&xre[i+32]);
                            zmm9 = _mm512_loadu_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_loadu_ps(&xim[i+32]);
                            zmm11= _mm512_loadu_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_storeu_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&yim[i+48], zmm15);
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]); // yre modified
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]); // yim modified
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm16,zmm17),
                                                                        _mm512_mul_ps(zmm18,zmm19));
                            _mm512_storeu_ps(&yre[i+64], zmm17);
                            zmm19  = _mm512_mul_ps(_mm512_mul_ps(zmm18,zmm17),
                                                                        _mm512_mul_ps(zmm16,zmm19));
                            _mm512_storeu_ps(&yim[i+64], zmm19);
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]); // yre modified
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]); // yim modified
                            zmm21 = _mm512_sub_ps(_mm512_mul_ps(zmm20,zmm21),
                                                       _mm512_mul_ps(zmm22,zmm23));
                            _mm512_storeu_ps(&yre[i+80], zmm21);
                            zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm22,zmm21),
                                                           _mm512_mul_ps(zmm20,zmm23));
                            _mm512_storeu_ps(&yim[i+80], zmm23);
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&yre[i+96]); // yre modified
                            zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            zmm27 = _mm512_loadu_ps(&yim[i+96]); // yim modified
                            zmm25 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                       _mm512_mul_ps(zmm25,zmm27));
                            _mm512_storeu_ps(&yre[i+96], zmm25);
                             zmm27  = _mm512_mul_ps(_mm512_mul_ps(zmm26,zmm25),
                                                           _mm512_mul_ps(zmm24,zmm27));
                            _mm512_storeu_ps(&yim[i+96], zmm27);
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&yre[i+112]);  // yre modified
                            zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            zmm31 = _mm512_loadu_ps(&yim[i+112]);  // yim modified
                            zmm29 = _mm512_sub_ps(_mm512_mul_ps(zmm28,zmm29),
                                                       _mm512_mul_ps(zmm30,zmm31));
                            _mm512_storeu_ps(&yre[i+112], zmm29);
                            zmm31  = _mm512_mul_ps(_mm512_mul_ps(zmm30,zmm29),
                                                           _mm512_mul_ps(zmm28,zmm31));
                            _mm512_storeu_ps(&yim[i+112], zmm31);
                            zmm1 = _mm512_loadu_ps(&xre[i+128]);
                            zmm2 = _mm512_loadu_ps(&yre[i+128]); // yre modified
                            zmm3 = _mm512_loadu_ps(&xim[i+128]);
                            zmm4 = _mm512_loadu_ps(&yim[i+128]); // yim modified
                            zmm2 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_storeu_ps(&yre[i+128], zmm2);
                            zmm4  = _mm512_mul_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_storeu_ps(&yim[i+128], zmm4);
                            zmm5 = _mm512_loadu_ps(&xre[i+144]);
                            zmm6 = _mm512_loadu_ps(&yre[i+144]); // yre modified
                            zmm7 = _mm512_loadu_ps(&xim[i+144]);
                            zmm8 = _mm512_loadu_ps(&yim[i+144]); // yim modified
                            zmm6 = _mm512_sub_ps(_mm512_mul_ps(zmm5,zmm6),
                                                       _mm512_mul_ps(zmm7,zmm8));
                            _mm512_storeu_ps(&yre[i+144], zmm6);
                             zmm8  = _mm512_mul_ps(_mm512_mul_ps(zmm7,zmm6),
                                                           _mm512_mul_ps(zmm5,zmm8));
                            _mm512_storeu_ps(&yim[i+144], zmm8);

                      }

                        for(; (i+127) < n; i += 128) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_storeu_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_loadu_ps(&xre[i+32]);
                            zmm9 = _mm512_loadu_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_loadu_ps(&xim[i+32]);
                            zmm11= _mm512_loadu_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_storeu_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&yim[i+48], zmm15);
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]); // yre modified
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]); // yim modified
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm16,zmm17),
                                                                        _mm512_mul_ps(zmm18,zmm19));
                            _mm512_storeu_ps(&yre[i+64], zmm17);
                            zmm19  = _mm512_mul_ps(_mm512_mul_ps(zmm18,zmm17),
                                                                        _mm512_mul_ps(zmm16,zmm19));
                            _mm512_storeu_ps(&yim[i+64], zmm19);
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]); // yre modified
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]); // yim modified
                            zmm21 = _mm512_sub_ps(_mm512_mul_ps(zmm20,zmm21),
                                                       _mm512_mul_ps(zmm22,zmm23));
                            _mm512_storeu_ps(&yre[i+80], zmm21);
                            zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm22,zmm21),
                                                           _mm512_mul_ps(zmm20,zmm23));
                            _mm512_storeu_ps(&yim[i+80], zmm23);
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&yre[i+96]); // yre modified
                            zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            zmm27 = _mm512_loadu_ps(&yim[i+96]); // yim modified
                            zmm25 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                       _mm512_mul_ps(zmm25,zmm27));
                            _mm512_storeu_ps(&yre[i+96], zmm25);
                             zmm27  = _mm512_mul_ps(_mm512_mul_ps(zmm26,zmm25),
                                                           _mm512_mul_ps(zmm24,zmm27));
                            _mm512_storeu_ps(&yim[i+96], zmm27);
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&yre[i+112]);  // yre modified
                            zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            zmm31 = _mm512_loadu_ps(&yim[i+112]);  // yim modified
                            zmm29 = _mm512_sub_ps(_mm512_mul_ps(zmm28,zmm29),
                                                       _mm512_mul_ps(zmm30,zmm31));
                            _mm512_storeu_ps(&yre[i+112], zmm29);
                            zmm31  = _mm512_mul_ps(_mm512_mul_ps(zmm30,zmm29),
                                                           _mm512_mul_ps(zmm28,zmm31));
                            _mm512_storeu_ps(&yim[i+112], zmm31);

                      }

                        for(; (i+95) < n; i += 96) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_storeu_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_loadu_ps(&xre[i+32]);
                            zmm9 = _mm512_loadu_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_loadu_ps(&xim[i+32]);
                            zmm11= _mm512_loadu_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_storeu_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&yim[i+48], zmm15);
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]); // yre modified
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]); // yim modified
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm16,zmm17),
                                                                        _mm512_mul_ps(zmm18,zmm19));
                            _mm512_storeu_ps(&yre[i+64], zmm17);
                            zmm19  = _mm512_mul_ps(_mm512_mul_ps(zmm18,zmm17),
                                                                        _mm512_mul_ps(zmm16,zmm19));
                            _mm512_storeu_ps(&yim[i+64], zmm19);
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]); // yre modified
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]); // yim modified
                            zmm21 = _mm512_sub_ps(_mm512_mul_ps(zmm20,zmm21),
                                                       _mm512_mul_ps(zmm22,zmm23));
                            _mm512_storeu_ps(&yre[i+80], zmm21);
                            zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm22,zmm21),
                                                           _mm512_mul_ps(zmm20,zmm23));
                            _mm512_storeu_ps(&yim[i+80], zmm23);

                      }

                        for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_storeu_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_loadu_ps(&xre[i+32]);
                            zmm9 = _mm512_loadu_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_loadu_ps(&xim[i+32]);
                            zmm11= _mm512_loadu_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_storeu_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&yim[i+48], zmm15);

                      }

                        for(; (i+31) < n; i += 32) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_storeu_ps(&yim[i+16], zmm7);

                      }

                        for(; (i+15) < n; i += 16) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);

                      }

                      for(; (i+0) < n; i += 1) {
                         yre[i] = (xre[i] * yre[i]) - (xim[i] * yim[i]);
		         yim[i] = (xim[i] * yre[i]) + (xre[i] * yim[i]);
                    }
                }


               
                   void gms::math::cmul_zmm16r4_12x_aip(const float * __restrict __ATTR_ALIGN__(64) xre,
                                             const float * __restrict __ATTR_ALIGN__(64) yre,
                                             float * __restrict __ATTR_ALIGN__(64) yre,
                                             float * __restrict __ATTR_ALIGN__(64) yim,
                                             const int32_t n) {

                        if(__builtin_expect(0==n,0)) {return;}

                          __m512 zmm0,zmm1,zmm2,zmm3;
                          __m512 zmm4,zmm5,zmm6,zmm7;
                          __m512 zmm8,zmm9,zmm10,zmm11;
                          __m512 zmm12,zmm13,zmm14,zmm15;
                          __m512 zmm16,zmm17,zmm18,zmm19;
                          __m512 zmm20,zmm21,zmm22,zmm23;
                          __m512 zmm24,zmm25,zmm26,zmm27;
                          __m512 zmm28,zmm29,zmm30,zmm31;
                         int32_t i;
                        for(i = 0; (i+191) < n; i += 192) {
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_store_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_store_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_load_ps(&xre[i+32]);
                            zmm9 = _mm512_load_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_load_ps(&xim[i+32]);
                            zmm11= _mm512_load_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_store_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_store_ps(&yim[i+32], zmm11);
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0);
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&yim[i+48], zmm15);
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]); // yre modified
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]); // yim modified
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm16,zmm17),
                                                                        _mm512_mul_ps(zmm18,zmm19));
                            _mm512_store_ps(&yre[i+64], zmm17);
                            zmm19  = _mm512_mul_ps(_mm512_mul_ps(zmm18,zmm17),
                                                                        _mm512_mul_ps(zmm16,zmm19));
                            _mm512_store_ps(&yim[i+64], zmm19);
                            zmm20 = _mm512_load_ps(&xre[i+80]);
                            zmm21 = _mm512_load_ps(&yre[i+80]); // yre modified
                            zmm22 = _mm512_load_ps(&xim[i+80]);
                            zmm23 = _mm512_load_ps(&yim[i+80]); // yim modified
                            zmm21 = _mm512_sub_ps(_mm512_mul_ps(zmm20,zmm21),
                                                       _mm512_mul_ps(zmm22,zmm23));
                            _mm512_store_ps(&yre[i+80], zmm21);
                            zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm22,zmm21),
                                                           _mm512_mul_ps(zmm20,zmm23));
                            _mm512_store_ps(&yim[i+80], zmm23);
                            _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+144],_MM_HINT_T0);
                            zmm24 = _mm512_load_ps(&xre[i+96]);
                            zmm25 = _mm512_load_ps(&yre[i+96]); // yre modified
                            zmm26 = _mm512_load_ps(&xim[i+96]);
                            zmm27 = _mm512_load_ps(&yim[i+96]); // yim modified
                            zmm25 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                       _mm512_mul_ps(zmm25,zmm27));
                            _mm512_store_ps(&yre[i+96], zmm25);
                             zmm27  = _mm512_mul_ps(_mm512_mul_ps(zmm26,zmm25),
                                                           _mm512_mul_ps(zmm24,zmm27));
                            _mm512_store_ps(&yim[i+96], zmm27);
                            zmm28 = _mm512_load_ps(&xre[i+112]);
                            zmm29 = _mm512_load_ps(&yre[i+112]);  // yre modified
                            zmm30 = _mm512_load_ps(&xim[i+112]);
                            zmm31 = _mm512_load_ps(&yim[i+112]);  // yim modified
                            zmm29 = _mm512_sub_ps(_mm512_mul_ps(zmm28,zmm29),
                                                       _mm512_mul_ps(zmm30,zmm31));
                            _mm512_store_ps(&yre[i+112], zmm29);
                            zmm31  = _mm512_mul_ps(_mm512_mul_ps(zmm30,zmm29),
                                                           _mm512_mul_ps(zmm28,zmm31));
                            _mm512_store_ps(&yim[i+112], zmm31);
                            zmm1 = _mm512_load_ps(&xre[i+128]);
                            zmm2 = _mm512_load_ps(&yre[i+128]); // yre modified
                            zmm3 = _mm512_load_ps(&xim[i+128]);
                            zmm4 = _mm512_load_ps(&yim[i+128]); // yim modified
                            zmm2 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_store_ps(&yre[i+128], zmm2);
                            zmm4  = _mm512_mul_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_storeu_ps(&yim[i+128], zmm4);
                            _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+192],_MM_HINT_T0);
                            zmm5 = _mm512_load_ps(&xre[i+144]);
                            zmm6 = _mm512_load_ps(&yre[i+144]); // yre modified
                            zmm7 = _mm512_load_ps(&xim[i+144]);
                            zmm8 = _mm512_load_ps(&yim[i+144]); // yim modified
                            zmm6 = _mm512_sub_ps(_mm512_mul_ps(zmm5,zmm6),
                                                       _mm512_mul_ps(zmm7,zmm8));
                            _mm512_store_ps(&yre[i+144], zmm6);
                             zmm8  = _mm512_mul_ps(_mm512_mul_ps(zmm7,zmm6),
                                                           _mm512_mul_ps(zmm5,zmm8));
                            _mm512_store_ps(&yim[i+144], zmm8);
                            zmm9 = _mm512_load_ps(&xre[i+160]);
                            zmm10= _mm512_load_ps(&yre[i+160]); // yre modified
                            zmm11= _mm512_load_ps(&xim[i+160]);
                            zmm12= _mm512_load_ps(&yim[i+160]); // yim modified
                            zmm10= _mm512_sub_ps(_mm512_mul_ps(zmm9,zmm10),
                                                       _mm512_mul_ps(zmm11,zmm12));
                            _mm512_store_ps(&yre[i+160], zmm10);
                            zmm12 = _mm512_mul_ps(_mm512_mul_ps(zmm11,zmm10),
                                                           _mm512_mul_ps(zmm9,zmm12));
                            _mm512_store_ps(&yim[i+160], zmm12);
                            zmm13 = _mm512_load_ps(&xre[i+176]);
                            zmm14 = _mm512_load_ps(&yre[i+176]); // yre modified
                            zmm15 = _mm512_load_ps(&xim[i+176]);
                            zmm16 = _mm512_load_ps(&yim[i+176]);  // yim modified
                            zmm14 = _mm512_sub_ps(_mm512_mul_ps(zmm13,zmm14),
                                                       _mm512_mul_ps(zmm15,zmm16));
                            _mm512_store_ps(&yre[i+176], zmm14);
                            zmm16 = _mm512_mul_ps(_mm512_mul_ps(zmm15,zmm14),
                                                           _mm512_mul_ps(zmm13,zmm16));
                            _mm512_store_ps(&yim[i+176], zmm16);
                      }

                        for(; (i+159) < n; i += 160) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_store_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_store_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_load_ps(&xre[i+32]);
                            zmm9 = _mm512_load_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_load_ps(&xim[i+32]);
                            zmm11= _mm512_load_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_store_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_store_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&yim[i+48], zmm15);
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]); // yre modified
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]); // yim modified
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm16,zmm17),
                                                                        _mm512_mul_ps(zmm18,zmm19));
                            _mm512_store_ps(&yre[i+64], zmm17);
                            zmm19  = _mm512_mul_ps(_mm512_mul_ps(zmm18,zmm17),
                                                                        _mm512_mul_ps(zmm16,zmm19));
                            _mm512_store_ps(&yim[i+64], zmm19);
                            zmm20 = _mm512_load_ps(&xre[i+80]);
                            zmm21 = _mm512_load_ps(&yre[i+80]); // yre modified
                            zmm22 = _mm512_load_ps(&xim[i+80]);
                            zmm23 = _mm512_load_ps(&yim[i+80]); // yim modified
                            zmm21 = _mm512_sub_ps(_mm512_mul_ps(zmm20,zmm21),
                                                       _mm512_mul_ps(zmm22,zmm23));
                            _mm512_store_ps(&yre[i+80], zmm21);
                            zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm22,zmm21),
                                                           _mm512_mul_ps(zmm20,zmm23));
                            _mm512_store_ps(&yim[i+80], zmm23);
                            zmm24 = _mm512_load_ps(&xre[i+96]);
                            zmm25 = _mm512_load_ps(&yre[i+96]); // yre modified
                            zmm26 = _mm512_load_ps(&xim[i+96]);
                            zmm27 = _mm512_load_ps(&yim[i+96]); // yim modified
                            zmm25 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                       _mm512_mul_ps(zmm25,zmm27));
                            _mm512_store_ps(&yre[i+96], zmm25);
                             zmm27  = _mm512_mul_ps(_mm512_mul_ps(zmm26,zmm25),
                                                           _mm512_mul_ps(zmm24,zmm27));
                            _mm512_store_ps(&yim[i+96], zmm27);
                            zmm28 = _mm512_loadups(&xre[i+112]);
                            zmm29 = _mm512_load_ps(&yre[i+112]);  // yre modified
                            zmm30 = _mm512_load_ps(&xim[i+112]);
                            zmm31 = _mm512_load_ps(&yim[i+112]);  // yim modified
                            zmm29 = _mm512_sub_ps(_mm512_mul_ps(zmm28,zmm29),
                                                       _mm512_mul_ps(zmm30,zmm31));
                            _mm512_store_ps(&yre[i+112], zmm29);
                            zmm31  = _mm512_mul_ps(_mm512_mul_ps(zmm30,zmm29),
                                                           _mm512_mul_ps(zmm28,zmm31));
                            _mm512_store_ps(&yim[i+112], zmm31);
                            zmm1 = _mm512_load_ps(&xre[i+128]);
                            zmm2 = _mm512_load_ps(&yre[i+128]); // yre modified
                            zmm3 = _mm512_load_ps(&xim[i+128]);
                            zmm4 = _mm512_load_ps(&yim[i+128]); // yim modified
                            zmm2 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_store_ps(&yre[i+128], zmm2);
                            zmm4  = _mm512_mul_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_store_ps(&yim[i+128], zmm4);
                            zmm5 = _mm512_load_ps(&xre[i+144]);
                            zmm6 = _mm512_load_ps(&yre[i+144]); // yre modified
                            zmm7 = _mm512_load_ps(&xim[i+144]);
                            zmm8 = _mm512_load_ps(&yim[i+144]); // yim modified
                            zmm6 = _mm512_sub_ps(_mm512_mul_ps(zmm5,zmm6),
                                                       _mm512_mul_ps(zmm7,zmm8));
                            _mm512_store_ps(&yre[i+144], zmm6);
                             zmm8  = _mm512_mul_ps(_mm512_mul_ps(zmm7,zmm6),
                                                           _mm512_mul_ps(zmm5,zmm8));
                            _mm512_store_ps(&yim[i+144], zmm8);

                      }

                        for(; (i+127) < n; i += 128) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_store_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_storeu_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_load_ps(&xre[i+32]);
                            zmm9 = _mm512_load_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_load_ps(&xim[i+32]);
                            zmm11= _mm512_load_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_store_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_store_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&yim[i+48], zmm15);
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]); // yre modified
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]); // yim modified
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm16,zmm17),
                                                                        _mm512_mul_ps(zmm18,zmm19));
                            _mm512_store_ps(&yre[i+64], zmm17);
                            zmm19  = _mm512_mul_ps(_mm512_mul_ps(zmm18,zmm17),
                                                                        _mm512_mul_ps(zmm16,zmm19));
                            _mm512_store_ps(&yim[i+64], zmm19);
                            zmm20 = _mm512_load_ps(&xre[i+80]);
                            zmm21 = _mm512_load_ps(&yre[i+80]); // yre modified
                            zmm22 = _mm512_load_ps(&xim[i+80]);
                            zmm23 = _mm512_load_ps(&yim[i+80]); // yim modified
                            zmm21 = _mm512_sub_ps(_mm512_mul_ps(zmm20,zmm21),
                                                       _mm512_mul_ps(zmm22,zmm23));
                            _mm512_store_ps(&yre[i+80], zmm21);
                            zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm22,zmm21),
                                                           _mm512_mul_ps(zmm20,zmm23));
                            _mm512_store_ps(&yim[i+80], zmm23);
                            zmm24 = _mm512_load_ps(&xre[i+96]);
                            zmm25 = _mm512_load_ps(&yre[i+96]); // yre modified
                            zmm26 = _mm512_load_ps(&xim[i+96]);
                            zmm27 = _mm512_load_ps(&yim[i+96]); // yim modified
                            zmm25 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                       _mm512_mul_ps(zmm25,zmm27));
                            _mm512_store_ps(&yre[i+96], zmm25);
                             zmm27  = _mm512_mul_ps(_mm512_mul_ps(zmm26,zmm25),
                                                           _mm512_mul_ps(zmm24,zmm27));
                            _mm512_store_ps(&yim[i+96], zmm27);
                            zmm28 = _mm512_load_ps(&xre[i+112]);
                            zmm29 = _mm512_load_ps(&yre[i+112]);  // yre modified
                            zmm30 = _mm512_load_ps(&xim[i+112]);
                            zmm31 = _mm512_load_ps(&yim[i+112]);  // yim modified
                            zmm29 = _mm512_sub_ps(_mm512_mul_ps(zmm28,zmm29),
                                                       _mm512_mul_ps(zmm30,zmm31));
                            _mm512_store_ps(&yre[i+112], zmm29);
                            zmm31  = _mm512_mul_ps(_mm512_mul_ps(zmm30,zmm29),
                                                           _mm512_mul_ps(zmm28,zmm31));
                            _mm512_store_ps(&yim[i+112], zmm31);

                      }

                        for(; (i+95) < n; i += 96) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_store_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_store_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_load_ps(&xre[i+32]);
                            zmm9 = _mm512_load_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_load_ps(&xim[i+32]);
                            zmm11= _mm512_load_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_store_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_store_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&yim[i+48], zmm15);
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]); // yre modified
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]); // yim modified
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm16,zmm17),
                                                                        _mm512_mul_ps(zmm18,zmm19));
                            _mm512_store_ps(&yre[i+64], zmm17);
                            zmm19  = _mm512_mul_ps(_mm512_mul_ps(zmm18,zmm17),
                                                                        _mm512_mul_ps(zmm16,zmm19));
                            _mm512_store_ps(&yim[i+64], zmm19);
                            zmm20 = _mm512_load_ps(&xre[i+80]);
                            zmm21 = _mm512_load_ps(&yre[i+80]); // yre modified
                            zmm22 = _mm512_load_ps(&xim[i+80]);
                            zmm23 = _mm512_load_ps(&yim[i+80]); // yim modified
                            zmm21 = _mm512_sub_ps(_mm512_mul_ps(zmm20,zmm21),
                                                       _mm512_mul_ps(zmm22,zmm23));
                            _mm512_store_ps(&yre[i+80], zmm21);
                            zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm22,zmm21),
                                                           _mm512_mul_ps(zmm20,zmm23));
                            _mm512_store_ps(&yim[i+80], zmm23);

                      }

                        for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_store_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_store_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_load_ps(&xre[i+32]);
                            zmm9 = _mm512_load_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_load_ps(&xim[i+32]);
                            zmm11= _mm512_load_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_store_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_store_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&yim[i+48], zmm15);

                      }

                        for(; (i+31) < n; i += 32) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_store_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_store_ps(&yim[i+16], zmm7);

                      }

                        for(; (i+15) < n; i += 16) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&yim[i+0], zmm3);

                      }

                      for(; (i+0) < n; i += 1) {
                         yre[i] = (xre[i] * yre[i]) - (xim[i] * yim[i]);
		         yim[i] = (xim[i] * yre[i]) + (xre[i] * yim[i]);
                    }
                }


                  
                   void gms::math::cmul_zmm16r4_10x_aip(const float * __restrict __ATTR_ALIGN__(64)  xre,
                                             const float * __restrict __ATTR_ALIGN__(64)  yre,
                                             float * __restrict __ATTR_ALIGN__(64)  yre,
                                             float * __restrict __ATTR_ALIGN__(64)  yim,
                                             const int32_t n) {

                        if(__builtin_expect(0==n,0)) {return;}

                          __m512 zmm0,zmm1,zmm2,zmm3;
                          __m512 zmm4,zmm5,zmm6,zmm7;
                          __m512 zmm8,zmm9,zmm10,zmm11;
                          __m512 zmm12,zmm13,zmm14,zmm15;
                          __m512 zmm16,zmm17,zmm18,zmm19;
                          __m512 zmm20,zmm21,zmm22,zmm23;
                          __m512 zmm24,zmm25,zmm26,zmm27;
                          __m512 zmm28,zmm29,zmm30,zmm31;
                         int32_t i;

                        for(i = 0; (i+159) < n; i += 160) {
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_store_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_store_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_load_ps(&xre[i+32]);
                            zmm9 = _mm512_load_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_load_ps(&xim[i+32]);
                            zmm11= _mm512_load_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_store_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_store_ps(&yim[i+32], zmm11);
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0);
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&yim[i+48], zmm15);
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]); // yre modified
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]); // yim modified
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm16,zmm17),
                                                                        _mm512_mul_ps(zmm18,zmm19));
                            _mm512_store_ps(&yre[i+64], zmm17);
                            zmm19  = _mm512_mul_ps(_mm512_mul_ps(zmm18,zmm17),
                                                                        _mm512_mul_ps(zmm16,zmm19));
                            _mm512_store_ps(&yim[i+64], zmm19);
                            zmm20 = _mm512_load_ps(&xre[i+80]);
                            zmm21 = _mm512_load_ps(&yre[i+80]); // yre modified
                            zmm22 = _mm512_load_ps(&xim[i+80]);
                            zmm23 = _mm512_load_ps(&yim[i+80]); // yim modified
                            zmm21 = _mm512_sub_ps(_mm512_mul_ps(zmm20,zmm21),
                                                       _mm512_mul_ps(zmm22,zmm23));
                            _mm512_store_ps(&yre[i+80], zmm21);
                            zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm22,zmm21),
                                                           _mm512_mul_ps(zmm20,zmm23));
                            _mm512_store_ps(&yim[i+80], zmm23);
                            _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+144],_MM_HINT_T0);
                            zmm24 = _mm512_load_ps(&xre[i+96]);
                            zmm25 = _mm512_load_ps(&yre[i+96]); // yre modified
                            zmm26 = _mm512_load_ps(&xim[i+96]);
                            zmm27 = _mm512_load_ps(&yim[i+96]); // yim modified
                            zmm25 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                       _mm512_mul_ps(zmm25,zmm27));
                            _mm512_store_ps(&yre[i+96], zmm25);
                             zmm27  = _mm512_mul_ps(_mm512_mul_ps(zmm26,zmm25),
                                                           _mm512_mul_ps(zmm24,zmm27));
                            _mm512_store_ps(&yim[i+96], zmm27);
                            zmm28 = _mm512_load_ps(&xre[i+112]);
                            zmm29 = _mm512_load_ps(&yre[i+112]);  // yre modified
                            zmm30 = _mm512_load_ps(&xim[i+112]);
                            zmm31 = _mm512_load_ps(&yim[i+112]);  // yim modified
                            zmm29 = _mm512_sub_ps(_mm512_mul_ps(zmm28,zmm29),
                                                       _mm512_mul_ps(zmm30,zmm31));
                            _mm512_store_ps(&yre[i+112], zmm29);
                            zmm31  = _mm512_mul_ps(_mm512_mul_ps(zmm30,zmm29),
                                                           _mm512_mul_ps(zmm28,zmm31));
                            _mm512_store_ps(&yim[i+112], zmm31);
                            zmm1 = _mm512_load_ps(&xre[i+128]);
                            zmm2 = _mm512_load_ps(&yre[i+128]); // yre modified
                            zmm3 = _mm512_load_ps(&xim[i+128]);
                            zmm4 = _mm512_load_ps(&yim[i+128]); // yim modified
                            zmm2 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_store_ps(&yre[i+128], zmm2);
                            zmm4  = _mm512_mul_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_storeu_ps(&yim[i+128], zmm4);
                            zmm5 = _mm512_load_ps(&xre[i+144]);
                            zmm6 = _mm512_load_ps(&yre[i+144]); // yre modified
                            zmm7 = _mm512_load_ps(&xim[i+144]);
                            zmm8 = _mm512_load_ps(&yim[i+144]); // yim modified
                            zmm6 = _mm512_sub_ps(_mm512_mul_ps(zmm5,zmm6),
                                                       _mm512_mul_ps(zmm7,zmm8));
                            _mm512_store_ps(&yre[i+144], zmm6);
                             zmm8  = _mm512_mul_ps(_mm512_mul_ps(zmm7,zmm6),
                                                           _mm512_mul_ps(zmm5,zmm8));
                            _mm512_store_ps(&yim[i+144], zmm8);

                       }

                      for(; (i+127) < n; i += 128) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_store_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_store_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_load_ps(&xre[i+32]);
                            zmm9 = _mm512_load_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_load_ps(&xim[i+32]);
                            zmm11= _mm512_load_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_store_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_store_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&yim[i+48], zmm15);
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]); // yre modified
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]); // yim modified
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm16,zmm17),
                                                                        _mm512_mul_ps(zmm18,zmm19));
                            _mm512_store_ps(&yre[i+64], zmm17);
                            zmm19  = _mm512_mul_ps(_mm512_mul_ps(zmm18,zmm17),
                                                                        _mm512_mul_ps(zmm16,zmm19));
                            _mm512_store_ps(&yim[i+64], zmm19);
                            zmm20 = _mm512_load_ps(&xre[i+80]);
                            zmm21 = _mm512_load_ps(&yre[i+80]); // yre modified
                            zmm22 = _mm512_load_ps(&xim[i+80]);
                            zmm23 = _mm512_load_ps(&yim[i+80]); // yim modified
                            zmm21 = _mm512_sub_ps(_mm512_mul_ps(zmm20,zmm21),
                                                       _mm512_mul_ps(zmm22,zmm23));
                            _mm512_store_ps(&yre[i+80], zmm21);
                            zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm22,zmm21),
                                                           _mm512_mul_ps(zmm20,zmm23));
                            _mm512_store_ps(&yim[i+80], zmm23);
                            zmm24 = _mm512_load_ps(&xre[i+96]);
                            zmm25 = _mm512_load_ps(&yre[i+96]); // yre modified
                            zmm26 = _mm512_load_ps(&xim[i+96]);
                            zmm27 = _mm512_load_ps(&yim[i+96]); // yim modified
                            zmm25 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                       _mm512_mul_ps(zmm25,zmm27));
                            _mm512_store_ps(&yre[i+96], zmm25);
                             zmm27  = _mm512_mul_ps(_mm512_mul_ps(zmm26,zmm25),
                                                           _mm512_mul_ps(zmm24,zmm27));
                            _mm512_store_ps(&yim[i+96], zmm27);
                            zmm28 = _mm512_load_ps(&xre[i+112]);
                            zmm29 = _mm512_load_ps(&yre[i+112]);  // yre modified
                            zmm30 = _mm512_load_ps(&xim[i+112]);
                            zmm31 = _mm512_load_ps(&yim[i+112]);  // yim modified
                            zmm29 = _mm512_sub_ps(_mm512_mul_ps(zmm28,zmm29),
                                                       _mm512_mul_ps(zmm30,zmm31));
                            _mm512_store_ps(&yre[i+112], zmm29);
                            zmm31  = _mm512_mul_ps(_mm512_mul_ps(zmm30,zmm29),
                                                           _mm512_mul_ps(zmm28,zmm31));
                            _mm512_store_ps(&yim[i+112], zmm31);

                      }

                       for(; (i+79) < n; i += 80) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_store_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_store_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_load_ps(&xre[i+32]);
                            zmm9 = _mm512_load_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_load_ps(&xim[i+32]);
                            zmm11= _mm512_load_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_store_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_store_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&yim[i+48], zmm15);
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]); // yre modified
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]); // yim modified
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm16,zmm17),
                                                                        _mm512_mul_ps(zmm18,zmm19));
                            _mm512_store_ps(&yre[i+64], zmm17);
                            zmm19  = _mm512_mul_ps(_mm512_mul_ps(zmm18,zmm17),
                                                                        _mm512_mul_ps(zmm16,zmm19));
                            _mm512_store_ps(&yim[i+64], zmm19);

                      }

                       for(; (i+63) < n; i += 64) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_store_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_store_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_load_ps(&xre[i+32]);
                            zmm9 = _mm512_load_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_load_ps(&xim[i+32]);
                            zmm11= _mm512_load_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_store_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_store_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&yim[i+48], zmm15);

                      }

                       for(; (i+31) < n; i += 32) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_store_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_store_ps(&yim[i+16], zmm7);

                      }

                       for(; (i+15) < n; i += 16) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&yim[i+0], zmm3);
                      }

                       for(; (i+0) < n; i += 1) {
                         yre[i] = (xre[i] * yre[i]) - (xim[i] * yim[i]);
		         yim[i] = (xim[i] * yre[i]) + (xre[i] * yim[i]);
                    }
                 }


                
                   void gms::math::cmul_zmm16r4_10x_uip(const float * __restrict xre,
                                             const float * __restrict yre,
                                             float * __restrict  yre,
                                             float * __restrict  yim,
                                             const int32_t n) {

                        if(__builtin_expect(0==n,0)) {return;}

                          __m512 zmm0,zmm1,zmm2,zmm3;
                          __m512 zmm4,zmm5,zmm6,zmm7;
                          __m512 zmm8,zmm9,zmm10,zmm11;
                          __m512 zmm12,zmm13,zmm14,zmm15;
                          __m512 zmm16,zmm17,zmm18,zmm19;
                          __m512 zmm20,zmm21,zmm22,zmm23;
                          __m512 zmm24,zmm25,zmm26,zmm27;
                          __m512 zmm28,zmm29,zmm30,zmm31;
                         int32_t i;

                        for(i = 0; (i+159) < n; i += 160) {
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_storeu_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_loadu_ps(&xre[i+32]);
                            zmm9 = _mm512_loadu_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_loadu_ps(&xim[i+32]);
                            zmm11= _mm512_loadu_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_storeu_ps(&yim[i+32], zmm11);
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0);
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&yim[i+48], zmm15);
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]); // yre modified
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]); // yim modified
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm16,zmm17),
                                                                        _mm512_mul_ps(zmm18,zmm19));
                            _mm512_storeu_ps(&yre[i+64], zmm17);
                            zmm19  = _mm512_mul_ps(_mm512_mul_ps(zmm18,zmm17),
                                                                        _mm512_mul_ps(zmm16,zmm19));
                            _mm512_storeu_ps(&yim[i+64], zmm19);
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]); // yre modified
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]); // yim modified
                            zmm21 = _mm512_sub_ps(_mm512_mul_ps(zmm20,zmm21),
                                                       _mm512_mul_ps(zmm22,zmm23));
                            _mm512_storeu_ps(&yre[i+80], zmm21);
                            zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm22,zmm21),
                                                           _mm512_mul_ps(zmm20,zmm23));
                            _mm512_storeu_ps(&yim[i+80], zmm23);
                            _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+144],_MM_HINT_T0);
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&yre[i+96]); // yre modified
                            zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            zmm27 = _mm512_loadu_ps(&yim[i+96]); // yim modified
                            zmm25 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                       _mm512_mul_ps(zmm25,zmm27));
                            _mm512_storeu_ps(&yre[i+96], zmm25);
                             zmm27  = _mm512_mul_ps(_mm512_mul_ps(zmm26,zmm25),
                                                           _mm512_mul_ps(zmm24,zmm27));
                            _mm512_storeu_ps(&yim[i+96], zmm27);
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&yre[i+112]);  // yre modified
                            zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            zmm31 = _mm512_loadu_ps(&yim[i+112]);  // yim modified
                            zmm29 = _mm512_sub_ps(_mm512_mul_ps(zmm28,zmm29),
                                                       _mm512_mul_ps(zmm30,zmm31));
                            _mm512_storeu_ps(&yre[i+112], zmm29);
                            zmm31  = _mm512_mul_ps(_mm512_mul_ps(zmm30,zmm29),
                                                           _mm512_mul_ps(zmm28,zmm31));
                            _mm512_storeu_ps(&yim[i+112], zmm31);
                            zmm1 = _mm512_loadu_ps(&xre[i+128]);
                            zmm2 = _mm512_loadu_ps(&yre[i+128]); // yre modified
                            zmm3 = _mm512_loadu_ps(&xim[i+128]);
                            zmm4 = _mm512_loadu_ps(&yim[i+128]); // yim modified
                            zmm2 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_storeu_ps(&yre[i+128], zmm2);
                            zmm4  = _mm512_mul_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_storeu_ps(&yim[i+128], zmm4);
                            zmm5 = _mm512_loadu_ps(&xre[i+144]);
                            zmm6 = _mm512_loadu_ps(&yre[i+144]); // yre modified
                            zmm7 = _mm512_loadu_ps(&xim[i+144]);
                            zmm8 = _mm512_loadu_ps(&yim[i+144]); // yim modified
                            zmm6 = _mm512_sub_ps(_mm512_mul_ps(zmm5,zmm6),
                                                       _mm512_mul_ps(zmm7,zmm8));
                            _mm512_storeu_ps(&yre[i+144], zmm6);
                             zmm8  = _mm512_mul_ps(_mm512_mul_ps(zmm7,zmm6),
                                                           _mm512_mul_ps(zmm5,zmm8));
                            _mm512_storeu_ps(&yim[i+144], zmm8);

                       }

                      for(; (i+127) < n; i += 128) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_storeu_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_loadu_ps(&xre[i+32]);
                            zmm9 = _mm512_loadu_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_loadu_ps(&xim[i+32]);
                            zmm11= _mm512_loadu_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_storeu_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&yim[i+48], zmm15);
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]); // yre modified
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]); // yim modified
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm16,zmm17),
                                                                        _mm512_mul_ps(zmm18,zmm19));
                            _mm512_storeu_ps(&yre[i+64], zmm17);
                            zmm19  = _mm512_mul_ps(_mm512_mul_ps(zmm18,zmm17),
                                                                        _mm512_mul_ps(zmm16,zmm19));
                            _mm512_storeu_ps(&yim[i+64], zmm19);
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]); // yre modified
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]); // yim modified
                            zmm21 = _mm512_sub_ps(_mm512_mul_ps(zmm20,zmm21),
                                                       _mm512_mul_ps(zmm22,zmm23));
                            _mm512_storeu_ps(&yre[i+80], zmm21);
                            zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm22,zmm21),
                                                           _mm512_mul_ps(zmm20,zmm23));
                            _mm512_storeu_ps(&yim[i+80], zmm23);
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&yre[i+96]); // yre modified
                            zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            zmm27 = _mm512_loadu_ps(&yim[i+96]); // yim modified
                            zmm25 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                       _mm512_mul_ps(zmm25,zmm27));
                            _mm512_storeu_ps(&yre[i+96], zmm25);
                             zmm27  = _mm512_mul_ps(_mm512_mul_ps(zmm26,zmm25),
                                                           _mm512_mul_ps(zmm24,zmm27));
                            _mm512_storeu_ps(&yim[i+96], zmm27);
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&yre[i+112]);  // yre modified
                            zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            zmm31 = _mm512_loadu_ps(&yim[i+112]);  // yim modified
                            zmm29 = _mm512_sub_ps(_mm512_mul_ps(zmm28,zmm29),
                                                       _mm512_mul_ps(zmm30,zmm31));
                            _mm512_storeu_ps(&yre[i+112], zmm29);
                            zmm31  = _mm512_mul_ps(_mm512_mul_ps(zmm30,zmm29),
                                                           _mm512_mul_ps(zmm28,zmm31));
                            _mm512_storeu_ps(&yim[i+112], zmm31);

                      }

                       for(; (i+79) < n; i += 80) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_storeu_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_loadu_ps(&xre[i+32]);
                            zmm9 = _mm512_loadu_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_loadu_ps(&xim[i+32]);
                            zmm11= _mm512_loadu_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_storeu_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&yim[i+48], zmm15);
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]); // yre modified
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]); // yim modified
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm16,zmm17),
                                                                        _mm512_mul_ps(zmm18,zmm19));
                            _mm512_storeu_ps(&yre[i+64], zmm17);
                            zmm19  = _mm512_mul_ps(_mm512_mul_ps(zmm18,zmm17),
                                                                        _mm512_mul_ps(zmm16,zmm19));
                            _mm512_storeu_ps(&yim[i+64], zmm19);

                      }

                       for(; (i+63) < n; i += 64) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_storeu_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_loadu_ps(&xre[i+32]);
                            zmm9 = _mm512_loadu_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_loadu_ps(&xim[i+32]);
                            zmm11= _mm512_loadu_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_storeu_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&yim[i+48], zmm15);

                      }

                       for(; (i+31) < n; i += 32) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_storeu_ps(&yim[i+16], zmm7);

                      }

                       for(; (i+15) < n; i += 16) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);
                      }

                       for(; (i+0) < n; i += 1) {
                         yre[i] = (xre[i] * yre[i]) - (xim[i] * yim[i]);
		         yim[i] = (xim[i] * yre[i]) + (xre[i] * yim[i]);
                    }
                 }


                  
                   void gms::math::cmul_zmm16r4_6x_uip( const float * __restrict  xre,
                                             const float * __restrict  yre,
                                             float * __restrict  yre,
                                             float * __restrict  yim,
                                             const int32_t n) {

                        if(__builtin_expect(0==n,0)) {return;}

                          __m512 zmm0,zmm1,zmm2,zmm3;
                          __m512 zmm4,zmm5,zmm6,zmm7;
                          __m512 zmm8,zmm9,zmm10,zmm11;
                          __m512 zmm12,zmm13,zmm14,zmm15;
                          __m512 zmm16,zmm17,zmm18,zmm19;
                          __m512 zmm20,zmm21,zmm22,zmm23;
                         int32_t i;

                        for(; (i+95) < n; i += 96) {
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_storeu_ps(&yim[i+16], zmm7);
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0);
                            zmm8 = _mm512_loadu_ps(&xre[i+32]);
                            zmm9 = _mm512_loadu_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_loadu_ps(&xim[i+32]);
                            zmm11= _mm512_loadu_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_storeu_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&yim[i+48], zmm15);
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]); // yre modified
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]); // yim modified
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm16,zmm17),
                                                                        _mm512_mul_ps(zmm18,zmm19));
                            _mm512_storeu_ps(&yre[i+64], zmm17);
                            zmm19  = _mm512_mul_ps(_mm512_mul_ps(zmm18,zmm17),
                                                                        _mm512_mul_ps(zmm16,zmm19));
                            _mm512_storeu_ps(&yim[i+64], zmm19);
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]); // yre modified
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]); // yim modified
                            zmm21 = _mm512_sub_ps(_mm512_mul_ps(zmm20,zmm21),
                                                       _mm512_mul_ps(zmm22,zmm23));
                            _mm512_storeu_ps(&yre[i+80], zmm21);
                            zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm22,zmm21),
                                                           _mm512_mul_ps(zmm20,zmm23));
                            _mm512_storeu_ps(&yim[i+80], zmm23);
                      }

                    for(; (i+79) < n; i += 80) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_storeu_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_loadu_ps(&xre[i+32]);
                            zmm9 = _mm512_loadu_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_loadu_ps(&xim[i+32]);
                            zmm11= _mm512_loadu_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_storeu_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeU_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeU_ps(&yim[i+48], zmm15);
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]); // yre modified
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]); // yim modified
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm16,zmm17),
                                                                        _mm512_mul_ps(zmm18,zmm19));
                            _mm512_storeu_ps(&yre[i+64], zmm17);
                            zmm19  = _mm512_mul_ps(_mm512_mul_ps(zmm18,zmm17),
                                                                        _mm512_mul_ps(zmm16,zmm19));
                            _mm512_storeu_ps(&yim[i+64], zmm19);
                    }

                    for(; (i+63) < n; i += 64) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_storeu_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_loadu_ps(&xre[i+32]);
                            zmm9 = _mm512_loadu_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_loadu_ps(&xim[i+32]);
                            zmm11= _mm512_loadu_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_storeu_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeU_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeU_ps(&yim[i+48], zmm15);

                   }

                    for(; (i+31) < n; i += 32) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_storeu_ps(&yim[i+16], zmm7);
                    }

                     for(; (i+15) < n; i += 16) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&yim[i+0], zmm3);
                    }

                     for(; (i+0) < n; i += 1) {
                         yre[i] = (xre[i] * yre[i]) - (xim[i] * yim[i]);
		         yim[i] = (xim[i] * yre[i]) + (xre[i] * yim[i]);
                    }

               }


                
                   void gms::math::cmul_zmm16r4_6x_aip( const float * __restrict  __ATTR_ALIGN__(64) xre,
                                             const float * __restrict  __ATTR_ALIGN__(64) yre,
                                             float * __restrict  __ATTR_ALIGN__(64) yre,
                                             float * __restrict  __ATTR_ALIGN__(64) yim,
                                             const int32_t n) {

                        if(__builtin_expect(0==n,0)) {return;}

                          __m512 zmm0,zmm1,zmm2,zmm3;
                          __m512 zmm4,zmm5,zmm6,zmm7;
                          __m512 zmm8,zmm9,zmm10,zmm11;
                          __m512 zmm12,zmm13,zmm14,zmm15;
                          __m512 zmm16,zmm17,zmm18,zmm19;
                          __m512 zmm20,zmm21,zmm22,zmm23;
                         int32_t i;

                        for(; (i+95) < n; i += 96) {
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_store_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_store_ps(&yim[i+16], zmm7);
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0);
                            zmm8 = _mm512_load_ps(&xre[i+32]);
                            zmm9 = _mm512_load_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_load_ps(&xim[i+32]);
                            zmm11= _mm512_load_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_store_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_store_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&yim[i+48], zmm15);
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]); // yre modified
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]); // yim modified
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm16,zmm17),
                                                                        _mm512_mul_ps(zmm18,zmm19));
                            _mm512_store_ps(&yre[i+64], zmm17);
                            zmm19  = _mm512_mul_ps(_mm512_mul_ps(zmm18,zmm17),
                                                                        _mm512_mul_ps(zmm16,zmm19));
                            _mm512_store_ps(&yim[i+64], zmm19);
                            zmm20 = _mm512_load_ps(&xre[i+80]);
                            zmm21 = _mm512_load_ps(&yre[i+80]); // yre modified
                            zmm22 = _mm512_load_ps(&xim[i+80]);
                            zmm23 = _mm512_load_ps(&yim[i+80]); // yim modified
                            zmm21 = _mm512_sub_ps(_mm512_mul_ps(zmm20,zmm21),
                                                       _mm512_mul_ps(zmm22,zmm23));
                            _mm512_store_ps(&yre[i+80], zmm21);
                            zmm23  = _mm512_mul_ps(_mm512_mul_ps(zmm22,zmm21),
                                                           _mm512_mul_ps(zmm20,zmm23));
                            _mm512_store_ps(&yim[i+80], zmm23);
                      }

                    for(; (i+79) < n; i += 80) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_store_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_store_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_load_ps(&xre[i+32]);
                            zmm9 = _mm512_load_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_load_ps(&xim[i+32]);
                            zmm11= _mm512_load_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_store_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_store_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&yim[i+48], zmm15);
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]); // yre modified
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]); // yim modified
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm16,zmm17),
                                                                        _mm512_mul_ps(zmm18,zmm19));
                            _mm512_store_ps(&yre[i+64], zmm17);
                            zmm19  = _mm512_mul_ps(_mm512_mul_ps(zmm18,zmm17),
                                                                        _mm512_mul_ps(zmm16,zmm19));
                            _mm512_store_ps(&yim[i+64], zmm19);
                    }

                    for(; (i+63) < n; i += 64) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_store_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_storeu_ps(&yim[i+16], zmm7);
                            zmm8 = _mm512_load_ps(&xre[i+32]);
                            zmm9 = _mm512_load_ps(&yre[i+32]); // yre modified
                            zmm10= _mm512_load_ps(&xim[i+32]);
                            zmm11= _mm512_load_ps(&yim[i+32]); // yim modified
                            zmm9 = _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                                        _mm512_mul_ps(zmm10,zmm11));
                            _mm512_store_ps(&yre[i+32], zmm9);
                            zmm11  = _mm512_mul_ps(_mm512_mul_ps(zmm10,zmm9),
                                                                        _mm512_mul_ps(zmm8,zmm11));
                            _mm512_store_ps(&yim[i+32], zmm11);
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]); // yre modified
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]); // yim modified
                            zmm13 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&yre[i+48], zmm13);
                            zmm15  = _mm512_mul_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&yim[i+48], zmm15);

                   }

                    for(; (i+31) < n; i += 32) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&yim[i+0], zmm3);
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);  // yre modified
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);  // yim modified
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm4,zmm5),
                                                                        _mm512_mul_ps(zmm6,zmm7));
                            _mm512_store_ps(&yre[i+16], zmm5);
                            zmm7 = _mm512_mul_ps(_mm512_mul_ps(zmm6,zmm5),
                                                                        _mm512_mul_ps(zmm7,zmm7));
                            _mm512_store_ps(&yim[i+16], zmm7);
                    }

                     for(; (i+15) < n; i += 16) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); // yre modified
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]); // yim modified
                            zmm1  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&yre[i+0], zmm1);
                             zmm3  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&yim[i+0], zmm3);
                    }

                     for(; (i+0) < n; i += 1) {
                         yre[i] = (xre[i] * yre[i]) - (xim[i] * yim[i]);
		         yim[i] = (xim[i] * yre[i]) + (xre[i] * yim[i]);
                    }

               }






                
















