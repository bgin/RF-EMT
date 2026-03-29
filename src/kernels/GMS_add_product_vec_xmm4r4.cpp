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
#include "GMS_add_product_vec_xmm4r4.h"


std::int32_t 
gms::math::
vaddpc_f32_1x128_1xf32_u16x_u(const float * __restrict__ vec_x,
                              float       * __restrict__ vec_y,
                              const float                c,
                              std::size_t                n) 
{
        if(__builtin_expect(n==0ull,0)) { return (-1);}

        while(((std::uintptr_t)&vec_y && 15) && n) 
        {
              const float x = *vec_x;
              const float y = *vec_y;
              const float z =  x+y*c;
              *vec_y        =  z;
               vec_x++;
               vec_y++;
               n--;
        }
        
        const __m128 vc{_mm_set1_ps(c)};
        std::size_t i;
        for(i = 0ull; (i+63ull) < n; i += 64ull) 
        {
            /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
            */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_storeu_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_storeu_ps(&vec_y[i+4ull],vz1);
                vx2 = _mm_loadu_ps(&vec_x[i+8ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+8ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_storeu_ps(&vec_y[i+8ull],vz2);
                vx3 = _mm_loadu_ps(&vec_x[i+12ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+12ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_storeu_ps(&vec_y[i+12ull],vz3);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_loadu_ps(&vec_x[i+16ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+16ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_storeu_ps(&vec_y[i+16ull],vz4);
                vx5 = _mm_loadu_ps(&vec_x[i+20ull]);
                vy5 = _mm_loadu_ps(&vec_y[i+20ull]);
                vz5 = _mm_fmadd_ps(vx5,vc,vy5);
                _mm_storeu_ps(&vec_y[i+20ull],vz5);
                vx6 = _mm_loadu_ps(&vec_x[i+24ull]);
                vy6 = _mm_loadu_ps(&vec_y[i+24ull]);
                vz6 = _mm_fmadd_ps(vx6,vc,vy6);
                _mm_storeu_ps(&vec_y[i+24ull],vz6);
                vx7 = _mm_loadu_ps(&vec_x[i+28ull]);
                vy7 = _mm_loadu_ps(&vec_y[i+28ull]);
                vz7 = _mm_fmadd_ps(vx7,vc,vy7);
                _mm_storeu_ps(&vec_y[i+28ull],vz7);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+32ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+32ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_storeu_ps(&vec_y[i+32ull],vz0);
                vx1 = _mm_loadu_ps(&vec_x[i+36ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+36ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_storeu_ps(&vec_y[i+36ull],vz1);
                vx2 = _mm_loadu_ps(&vec_x[i+40ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+40ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_storeu_ps(&vec_y[i+40ull],vz2);
                vx3 = _mm_loadu_ps(&vec_x[i+44ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+44ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_storeu_ps(&vec_y[i+44ull],vz3);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+48ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_loadu_ps(&vec_x[i+48ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+48ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_storeu_ps(&vec_y[i+48ull],vz4);
                vx5 = _mm_loadu_ps(&vec_x[i+52ull]);
                vy5 = _mm_loadu_ps(&vec_y[i+52ull]);
                vz5 = _mm_fmadd_ps(vx5,vc,vy5);
                _mm_storeu_ps(&vec_y[i+52ull],vz5); 
                vx6 = _mm_loadu_ps(&vec_x[i+56ull]);
                vy6 = _mm_loadu_ps(&vec_y[i+56ull]);
                vz6 = _mm_fmadd_ps(vx6,vc,vy6);
                _mm_storeu_ps(&vec_y[i+56ull],vz6);
                vx7 = _mm_loadu_ps(&vec_x[i+60ull]);
                vy7 = _mm_loadu_ps(&vec_y[i+60ull]);
                vz7 = _mm_fmadd_ps(vx7,vc,vy7);
                _mm_storeu_ps(&vec_y[i+60ull],vz7);

        }

        for(; (i+47ull) < n; i += 48ull) 
        {
               /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_storeu_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_storeu_ps(&vec_y[i+4ull],vz1);
                vx2 = _mm_loadu_ps(&vec_x[i+8ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+8ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_storeu_ps(&vec_y[i+8ull],vz2);
                vx3 = _mm_loadu_ps(&vec_x[i+12ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+12ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_storeu_ps(&vec_y[i+12ull],vz3);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_loadu_ps(&vec_x[i+16ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+16ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_storeu_ps(&vec_y[i+16ull],vz4);
                vx5 = _mm_loadu_ps(&vec_x[i+20ull]);
                vy5 = _mm_loadu_ps(&vec_y[i+20ull]);
                vz5 = _mm_fmadd_ps(vx5,vc,vy5);
                _mm_storeu_ps(&vec_y[i+20ull],vz5);
                vx6 = _mm_loadu_ps(&vec_x[i+24ull]);
                vy6 = _mm_loadu_ps(&vec_y[i+24ull]);
                vz6 = _mm_fmadd_ps(vx6,vc,vy6);
                _mm_storeu_ps(&vec_y[i+24ull],vz6);
                vx7 = _mm_loadu_ps(&vec_x[i+28ull]);
                vy7 = _mm_loadu_ps(&vec_y[i+28ull]);
                vz7 = _mm_fmadd_ps(vx7,vc,vy7);
                _mm_storeu_ps(&vec_y[i+28ull],vz7);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+32ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+32ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_storeu_ps(&vec_y[i+32ull],vz0);
                vx1 = _mm_loadu_ps(&vec_x[i+36ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+36ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_storeu_ps(&vec_y[i+36ull],vz1);
                vx2 = _mm_loadu_ps(&vec_x[i+40ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+40ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_storeu_ps(&vec_y[i+40ull],vz2);
                vx3 = _mm_loadu_ps(&vec_x[i+44ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+44ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_storeu_ps(&vec_y[i+44ull],vz3);
        }

        for(; (i+31ull) < n; i += 32ull) 
        {
               /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;

#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_storeu_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_storeu_ps(&vec_y[i+4ull],vz1);
                vx2 = _mm_loadu_ps(&vec_x[i+8ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+8ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_storeu_ps(&vec_y[i+8ull],vz2);
                vx3 = _mm_loadu_ps(&vec_x[i+12ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+12ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_storeu_ps(&vec_y[i+12ull],vz3);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_loadu_ps(&vec_x[i+16ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+16ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_storeu_ps(&vec_y[i+16ull],vz4);
                vx5 = _mm_loadu_ps(&vec_x[i+20ull]);
                vy5 = _mm_loadu_ps(&vec_y[i+20ull]);
                vz5 = _mm_fmadd_ps(vx5,vc,vy5);
                _mm_storeu_ps(&vec_y[i+20ull],vz5);
                vx6 = _mm_loadu_ps(&vec_x[i+24ull]);
                vy6 = _mm_loadu_ps(&vec_y[i+24ull]);
                vz6 = _mm_fmadd_ps(vx6,vc,vy6);
                _mm_storeu_ps(&vec_y[i+24ull],vz6);
                vx7 = _mm_loadu_ps(&vec_x[i+28ull]);
                vy7 = _mm_loadu_ps(&vec_y[i+28ull]);
                vz7 = _mm_fmadd_ps(vx7,vc,vy7);
                _mm_storeu_ps(&vec_y[i+28ull],vz7);
        }

        for(; (i+19ull) < n; i += 20ull) 
        {
              /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4; 
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_storeu_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_storeu_ps(&vec_y[i+4ull],vz1);
                vx2 = _mm_loadu_ps(&vec_x[i+8ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+8ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_storeu_ps(&vec_y[i+8ull],vz2);
                vx3 = _mm_loadu_ps(&vec_x[i+12ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+12ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_storeu_ps(&vec_y[i+12ull],vz3);
                vx4 = _mm_loadu_ps(&vec_x[i+16ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+16ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_storeu_ps(&vec_y[i+16ull],vz4);
        }

        for(; (i+7ull) < n; i += 8ull) 
        {
             /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                             
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_storeu_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_storeu_ps(&vec_y[i+4ull],vz1);
        }

        for(; (i+3ull) < n; i += 4ull) 
        {
             /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                                             
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_storeu_ps(&vec_y[i+0ull],vz0);
        }

        for(; (i+0ull) < n; i += 1ull) 
        {
             const float x{vec_x[i]};
             const float y{vec_y[i]};
             const float z{x+y*c};
             vec_y[i] = z;
        }

        return (0);
}


std::int32_t 
gms::math::
vaddpc_f32_1x128_1xf32_u16x_a(const float * __restrict__ vec_x,
                              float       * __restrict__ vec_y,
                              const float                c,
                              std::size_t                n) 
{
        if(__builtin_expect(n==0ull,0)) { return (-1);}
               
        const __m128 vc{_mm_set1_ps(c)};
        std::size_t i;
        for(i = 0ull; (i+63ull) < n; i += 64ull) 
        {
            /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
            */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_store_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_store_ps(&vec_y[i+4ull],vz1);
                vx2 = _mm_load_ps(&vec_x[i+8ull]);
                vy2 = _mm_load_ps(&vec_y[i+8ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_store_ps(&vec_y[i+8ull],vz2);
                vx3 = _mm_load_ps(&vec_x[i+12ull]);
                vy3 = _mm_load_ps(&vec_y[i+12ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_store_ps(&vec_y[i+12ull],vz3);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_load_ps(&vec_x[i+16ull]);
                vy4 = _mm_load_ps(&vec_y[i+16ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_store_ps(&vec_y[i+16ull],vz4);
                vx5 = _mm_load_ps(&vec_x[i+20ull]);
                vy5 = _mm_load_ps(&vec_y[i+20ull]);
                vz5 = _mm_fmadd_ps(vx5,vc,vy5);
                _mm_store_ps(&vec_y[i+20ull],vz5);
                vx6 = _mm_load_ps(&vec_x[i+24ull]);
                vy6 = _mm_load_ps(&vec_y[i+24ull]);
                vz6 = _mm_fmadd_ps(vx6,vc,vy6);
                _mm_store_ps(&vec_y[i+24ull],vz6);
                vx7 = _mm_load_ps(&vec_x[i+28ull]);
                vy7 = _mm_load_ps(&vec_y[i+28ull]);
                vz7 = _mm_fmadd_ps(vx7,vc,vy7);
                _mm_store_ps(&vec_y[i+28ull],vz7);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+32ull]);
                vy0 = _mm_load_ps(&vec_y[i+32ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_store_ps(&vec_y[i+32ull],vz0);
                vx1 = _mm_load_ps(&vec_x[i+36ull]);
                vy1 = _mm_load_ps(&vec_y[i+36ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_store_ps(&vec_y[i+36ull],vz1);
                vx2 = _mm_load_ps(&vec_x[i+40ull]);
                vy2 = _mm_load_ps(&vec_y[i+40ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_store_ps(&vec_y[i+40ull],vz2);
                vx3 = _mm_load_ps(&vec_x[i+44ull]);
                vy3 = _mm_load_ps(&vec_y[i+44ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_store_ps(&vec_y[i+44ull],vz3);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+48ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_load_ps(&vec_x[i+48ull]);
                vy4 = _mm_load_ps(&vec_y[i+48ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_store_ps(&vec_y[i+48ull],vz4);
                vx5 = _mm_load_ps(&vec_x[i+52ull]);
                vy5 = _mm_load_ps(&vec_y[i+52ull]);
                vz5 = _mm_fmadd_ps(vx5,vc,vy5);
                _mm_store_ps(&vec_y[i+52ull],vz5); 
                vx6 = _mm_load_ps(&vec_x[i+56ull]);
                vy6 = _mm_load_ps(&vec_y[i+56ull]);
                vz6 = _mm_fmadd_ps(vx6,vc,vy6);
                _mm_store_ps(&vec_y[i+56ull],vz6);
                vx7 = _mm_load_ps(&vec_x[i+60ull]);
                vy7 = _mm_load_ps(&vec_y[i+60ull]);
                vz7 = _mm_fmadd_ps(vx7,vc,vy7);
                _mm_store_ps(&vec_y[i+60ull],vz7);

        }

        for(; (i+47ull) < n; i += 48ull) 
        {
               /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_store_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_store_ps(&vec_y[i+4ull],vz1);
                vx2 = _mm_load_ps(&vec_x[i+8ull]);
                vy2 = _mm_load_ps(&vec_y[i+8ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_store_ps(&vec_y[i+8ull],vz2);
                vx3 = _mm_load_ps(&vec_x[i+12ull]);
                vy3 = _mm_load_ps(&vec_y[i+12ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_store_ps(&vec_y[i+12ull],vz3);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_load_ps(&vec_x[i+16ull]);
                vy4 = _mm_load_ps(&vec_y[i+16ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_store_ps(&vec_y[i+16ull],vz4);
                vx5 = _mm_load_ps(&vec_x[i+20ull]);
                vy5 = _mm_load_ps(&vec_y[i+20ull]);
                vz5 = _mm_fmadd_ps(vx5,vc,vy5);
                _mm_store_ps(&vec_y[i+20ull],vz5);
                vx6 = _mm_load_ps(&vec_x[i+24ull]);
                vy6 = _mm_load_ps(&vec_y[i+24ull]);
                vz6 = _mm_fmadd_ps(vx6,vc,vy6);
                _mm_store_ps(&vec_y[i+24ull],vz6);
                vx7 = _mm_load_ps(&vec_x[i+28ull]);
                vy7 = _mm_load_ps(&vec_y[i+28ull]);
                vz7 = _mm_fmadd_ps(vx7,vc,vy7);
                _mm_store_ps(&vec_y[i+28ull],vz7);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+32ull]);
                vy0 = _mm_load_ps(&vec_y[i+32ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_store_ps(&vec_y[i+32ull],vz0);
                vx1 = _mm_load_ps(&vec_x[i+36ull]);
                vy1 = _mm_load_ps(&vec_y[i+36ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_store_ps(&vec_y[i+36ull],vz1);
                vx2 = _mm_load_ps(&vec_x[i+40ull]);
                vy2 = _mm_load_ps(&vec_y[i+40ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_store_ps(&vec_y[i+40ull],vz2);
                vx3 = _mm_load_ps(&vec_x[i+44ull]);
                vy3 = _mm_load_ps(&vec_y[i+44ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_store_ps(&vec_y[i+44ull],vz3);
        }

        for(; (i+31ull) < n; i += 32ull) 
        {
               /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;

#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_store_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_store_ps(&vec_y[i+4ull],vz1);
                vx2 = _mm_load_ps(&vec_x[i+8ull]);
                vy2 = _mm_load_ps(&vec_y[i+8ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_store_ps(&vec_y[i+8ull],vz2);
                vx3 = _mm_load_ps(&vec_x[i+12ull]);
                vy3 = _mm_load_ps(&vec_y[i+12ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_store_ps(&vec_y[i+12ull],vz3);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_load_ps(&vec_x[i+16ull]);
                vy4 = _mm_load_ps(&vec_y[i+16ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_store_ps(&vec_y[i+16ull],vz4);
                vx5 = _mm_load_ps(&vec_x[i+20ull]);
                vy5 = _mm_load_ps(&vec_y[i+20ull]);
                vz5 = _mm_fmadd_ps(vx5,vc,vy5);
                _mm_store_ps(&vec_y[i+20ull],vz5);
                vx6 = _mm_load_ps(&vec_x[i+24ull]);
                vy6 = _mm_load_ps(&vec_y[i+24ull]);
                vz6 = _mm_fmadd_ps(vx6,vc,vy6);
                _mm_store_ps(&vec_y[i+24ull],vz6);
                vx7 = _mm_load_ps(&vec_x[i+28ull]);
                vy7 = _mm_load_ps(&vec_y[i+28ull]);
                vz7 = _mm_fmadd_ps(vx7,vc,vy7);
                _mm_store_ps(&vec_y[i+28ull],vz7);
        }

        for(; (i+19ull) < n; i += 20ull) 
        {
              /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4; 
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_store_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_store_ps(&vec_y[i+4ull],vz1);
                vx2 = _mm_load_ps(&vec_x[i+8ull]);
                vy2 = _mm_load_ps(&vec_y[i+8ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_store_ps(&vec_y[i+8ull],vz2);
                vx3 = _mm_load_ps(&vec_x[i+12ull]);
                vy3 = _mm_load_ps(&vec_y[i+12ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_store_ps(&vec_y[i+12ull],vz3);
                vx4 = _mm_load_ps(&vec_x[i+16ull]);
                vy4 = _mm_load_ps(&vec_y[i+16ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_store_ps(&vec_y[i+16ull],vz4);
        }

        for(; (i+7ull) < n; i += 8ull) 
        {
             /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                             
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_store_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_store_ps(&vec_y[i+4ull],vz1);
        }

        for(; (i+3ull) < n; i += 4ull) 
        {
             /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                                             
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_store_ps(&vec_y[i+0ull],vz0);
        }

        for(; (i+0ull) < n; i += 1ull) 
        {
             const float x{vec_x[i]};
             const float y{vec_y[i]};
             const float z{x+y*c};
             vec_y[i] = z;
        }

        return (0);
}


std::int32_t 
gms::math::
vaddpc_f32_1x128_1xf32_u10x_u(const float * __restrict__ vec_x,
                              float       * __restrict__ vec_y,
                              const float                c,
                              std::size_t                n) 
{
        if(__builtin_expect(n==0ull,0)) { return (-1);}

        while(((std::uintptr_t)&vec_y && 15) && n) 
        {
              const float x = *vec_x;
              const float y = *vec_y;
              const float z =  x+y*c;
              *vec_y        =  z;
               vec_x++;
               vec_y++;
               n--;
        }
        
        const __m128 vc{_mm_set1_ps(c)};
        std::size_t i;
        for(i = 0ull; (i+39ull) < n; i += 40ull) 
        {
            /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
            */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_storeu_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_storeu_ps(&vec_y[i+4ull],vz1);
                vx2 = _mm_loadu_ps(&vec_x[i+8ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+8ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_storeu_ps(&vec_y[i+8ull],vz2);
                vx3 = _mm_loadu_ps(&vec_x[i+12ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+12ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_storeu_ps(&vec_y[i+12ull],vz3);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_loadu_ps(&vec_x[i+16ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+16ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_storeu_ps(&vec_y[i+16ull],vz4);
                vx5 = _mm_loadu_ps(&vec_x[i+20ull]);
                vy5 = _mm_loadu_ps(&vec_y[i+20ull]);
                vz5 = _mm_fmadd_ps(vx5,vc,vy5);
                _mm_storeu_ps(&vec_y[i+20ull],vz5);
                vx6 = _mm_loadu_ps(&vec_x[i+24ull]);
                vy6 = _mm_loadu_ps(&vec_y[i+24ull]);
                vz6 = _mm_fmadd_ps(vx6,vc,vy6);
                _mm_storeu_ps(&vec_y[i+24ull],vz6);
                vx7 = _mm_loadu_ps(&vec_x[i+28ull]);
                vy7 = _mm_loadu_ps(&vec_y[i+28ull]);
                vz7 = _mm_fmadd_ps(vx7,vc,vy7);
                _mm_storeu_ps(&vec_y[i+28ull],vz7);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+32ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+32ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_storeu_ps(&vec_y[i+32ull],vz0);
                vx1 = _mm_loadu_ps(&vec_x[i+36ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+36ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_storeu_ps(&vec_y[i+36ull],vz1);
                
        }

        for(; (i+31ull) < n; i += 32ull) 
        {
               /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;

#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_storeu_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_storeu_ps(&vec_y[i+4ull],vz1);
                vx2 = _mm_loadu_ps(&vec_x[i+8ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+8ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_storeu_ps(&vec_y[i+8ull],vz2);
                vx3 = _mm_loadu_ps(&vec_x[i+12ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+12ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_storeu_ps(&vec_y[i+12ull],vz3);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_loadu_ps(&vec_x[i+16ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+16ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_storeu_ps(&vec_y[i+16ull],vz4);
                vx5 = _mm_loadu_ps(&vec_x[i+20ull]);
                vy5 = _mm_loadu_ps(&vec_y[i+20ull]);
                vz5 = _mm_fmadd_ps(vx5,vc,vy5);
                _mm_storeu_ps(&vec_y[i+20ull],vz5);
                vx6 = _mm_loadu_ps(&vec_x[i+24ull]);
                vy6 = _mm_loadu_ps(&vec_y[i+24ull]);
                vz6 = _mm_fmadd_ps(vx6,vc,vy6);
                _mm_storeu_ps(&vec_y[i+24ull],vz6);
                vx7 = _mm_loadu_ps(&vec_x[i+28ull]);
                vy7 = _mm_loadu_ps(&vec_y[i+28ull]);
                vz7 = _mm_fmadd_ps(vx7,vc,vy7);
                _mm_storeu_ps(&vec_y[i+28ull],vz7);
        }

        for(; (i+19ull) < n; i += 20ull) 
        {
              /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4; 
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_storeu_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_storeu_ps(&vec_y[i+4ull],vz1);
                vx2 = _mm_loadu_ps(&vec_x[i+8ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+8ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_storeu_ps(&vec_y[i+8ull],vz2);
                vx3 = _mm_loadu_ps(&vec_x[i+12ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+12ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_storeu_ps(&vec_y[i+12ull],vz3);
                vx4 = _mm_loadu_ps(&vec_x[i+16ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+16ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_storeu_ps(&vec_y[i+16ull],vz4);
        }

        for(; (i+7ull) < n; i += 8ull) 
        {
             /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                             
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_storeu_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_storeu_ps(&vec_y[i+4ull],vz1);
        }

        for(; (i+3ull) < n; i += 4ull) 
        {
             /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                                             
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_storeu_ps(&vec_y[i+0ull],vz0);
        }
        
        for(; (i+0ull) < n; i += 1ull) 
        {
             const float x{vec_x[i]};
             const float y{vec_y[i]};
             const float z{x+y*c};
             vec_y[i] = z;
        }

        return (0);
}


std::int32_t 
gms::math::
vaddpc_f32_1x128_1xf32_u10x_a(const float * __restrict__ vec_x,
                              float       * __restrict__ vec_y,
                              const float                c,
                              std::size_t                n) 
{
        if(__builtin_expect(n==0ull,0)) { return (-1);}

              
        const __m128 vc{_mm_set1_ps(c)};
        std::size_t i;
        for(i = 0ull; (i+39ull) < n; i += 40ull) 
        {
            /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
            */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_store_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_store_ps(&vec_y[i+4ull],vz1);
                vx2 = _mm_load_ps(&vec_x[i+8ull]);
                vy2 = _mm_load_ps(&vec_y[i+8ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_store_ps(&vec_y[i+8ull],vz2);
                vx3 = _mm_load_ps(&vec_x[i+12ull]);
                vy3 = _mm_load_ps(&vec_y[i+12ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_store_ps(&vec_y[i+12ull],vz3);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_load_ps(&vec_x[i+16ull]);
                vy4 = _mm_load_ps(&vec_y[i+16ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_store_ps(&vec_y[i+16ull],vz4);
                vx5 = _mm_load_ps(&vec_x[i+20ull]);
                vy5 = _mm_load_ps(&vec_y[i+20ull]);
                vz5 = _mm_fmadd_ps(vx5,vc,vy5);
                _mm_store_ps(&vec_y[i+20ull],vz5);
                vx6 = _mm_load_ps(&vec_x[i+24ull]);
                vy6 = _mm_load_ps(&vec_y[i+24ull]);
                vz6 = _mm_fmadd_ps(vx6,vc,vy6);
                _mm_store_ps(&vec_y[i+24ull],vz6);
                vx7 = _mm_load_ps(&vec_x[i+28ull]);
                vy7 = _mm_load_ps(&vec_y[i+28ull]);
                vz7 = _mm_fmadd_ps(vx7,vc,vy7);
                _mm_store_ps(&vec_y[i+28ull],vz7);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+32ull]);
                vy0 = _mm_load_ps(&vec_y[i+32ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_store_ps(&vec_y[i+32ull],vz0);
                vx1 = _mm_load_ps(&vec_x[i+36ull]);
                vy1 = _mm_load_ps(&vec_y[i+36ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_store_ps(&vec_y[i+36ull],vz1);
                
        }

        for(; (i+31ull) < n; i += 32ull) 
        {
               /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;

#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_store_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_store_ps(&vec_y[i+4ull],vz1);
                vx2 = _mm_load_ps(&vec_x[i+8ull]);
                vy2 = _mm_load_ps(&vec_y[i+8ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_store_ps(&vec_y[i+8ull],vz2);
                vx3 = _mm_load_ps(&vec_x[i+12ull]);
                vy3 = _mm_load_ps(&vec_y[i+12ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_store_ps(&vec_y[i+12ull],vz3);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_load_ps(&vec_x[i+16ull]);
                vy4 = _mm_load_ps(&vec_y[i+16ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_store_ps(&vec_y[i+16ull],vz4);
                vx5 = _mm_load_ps(&vec_x[i+20ull]);
                vy5 = _mm_load_ps(&vec_y[i+20ull]);
                vz5 = _mm_fmadd_ps(vx5,vc,vy5);
                _mm_store_ps(&vec_y[i+20ull],vz5);
                vx6 = _mm_load_ps(&vec_x[i+24ull]);
                vy6 = _mm_load_ps(&vec_y[i+24ull]);
                vz6 = _mm_fmadd_ps(vx6,vc,vy6);
                _mm_store_ps(&vec_y[i+24ull],vz6);
                vx7 = _mm_load_ps(&vec_x[i+28ull]);
                vy7 = _mm_load_ps(&vec_y[i+28ull]);
                vz7 = _mm_fmadd_ps(vx7,vc,vy7);
                _mm_store_ps(&vec_y[i+28ull],vz7);
        }

        for(; (i+19ull) < n; i += 20ull) 
        {
              /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4; 
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_store_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_store_ps(&vec_y[i+4ull],vz1);
                vx2 = _mm_load_ps(&vec_x[i+8ull]);
                vy2 = _mm_load_ps(&vec_y[i+8ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_store_ps(&vec_y[i+8ull],vz2);
                vx3 = _mm_load_ps(&vec_x[i+12ull]);
                vy3 = _mm_load_ps(&vec_y[i+12ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_store_ps(&vec_y[i+12ull],vz3);
                vx4 = _mm_load_ps(&vec_x[i+16ull]);
                vy4 = _mm_load_ps(&vec_y[i+16ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_store_ps(&vec_y[i+16ull],vz4);
        }

        for(; (i+7ull) < n; i += 8ull) 
        {
             /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                             
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_store_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_store_ps(&vec_y[i+4ull],vz1);
        }

        for(; (i+3ull) < n; i += 4ull) 
        {
             /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                                             
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_store_ps(&vec_y[i+0ull],vz0);
        }
        
        for(; (i+0ull) < n; i += 1ull) 
        {
             const float x{vec_x[i]};
             const float y{vec_y[i]};
             const float z{x+y*c};
             vec_y[i] = z;
        }

        return (0);
}


std::int32_t 
gms::math::
vaddpc_f32_1x128_1xf32_u6x_u( const float * __restrict__ vec_x,
                              float       * __restrict__ vec_y,
                              const float                c,
                              std::size_t                n) 
{
        if(__builtin_expect(n==0ull,0)) { return (-1);}

        while(((std::uintptr_t)&vec_y && 15) && n) 
        {
              const float x = *vec_x;
              const float y = *vec_y;
              const float z =  x+y*c;
              *vec_y        =  z;
               vec_x++;
               vec_y++;
               n--;
        }
        
        const __m128 vc{_mm_set1_ps(c)};
        std::size_t i;
        for(i = 0ull; (i+23ull) < n; i += 24ull) 
        {
            /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
            */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_storeu_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_storeu_ps(&vec_y[i+4ull],vz1);
                vx2 = _mm_loadu_ps(&vec_x[i+8ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+8ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_storeu_ps(&vec_y[i+8ull],vz2);
                vx3 = _mm_loadu_ps(&vec_x[i+12ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+12ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_storeu_ps(&vec_y[i+12ull],vz3);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_loadu_ps(&vec_x[i+16ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+16ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_storeu_ps(&vec_y[i+16ull],vz4);
                vx5 = _mm_loadu_ps(&vec_x[i+20ull]);
                vy5 = _mm_loadu_ps(&vec_y[i+20ull]);
                vz5 = _mm_fmadd_ps(vx5,vc,vy5);
                _mm_storeu_ps(&vec_y[i+20ull],vz5);
        }

       
        for(; (i+19ull) < n; i += 20ull) 
        {
              /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4; 
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_storeu_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_storeu_ps(&vec_y[i+4ull],vz1);
                vx2 = _mm_loadu_ps(&vec_x[i+8ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+8ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_storeu_ps(&vec_y[i+8ull],vz2);
                vx3 = _mm_loadu_ps(&vec_x[i+12ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+12ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_storeu_ps(&vec_y[i+12ull],vz3);
                vx4 = _mm_loadu_ps(&vec_x[i+16ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+16ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_storeu_ps(&vec_y[i+16ull],vz4);
        }

        for(; (i+7ull) < n; i += 8ull) 
        {
             /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                             
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_storeu_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_storeu_ps(&vec_y[i+4ull],vz1);
        }

        for(; (i+3ull) < n; i += 4ull) 
        {
             /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                                             
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_storeu_ps(&vec_y[i+0ull],vz0);
        }
        
        for(; (i+0ull) < n; i += 1ull) 
        {
             const float x{vec_x[i]};
             const float y{vec_y[i]};
             const float z{x+y*c};
             vec_y[i] = z;
        }

        return (0);
}


std::int32_t 
gms::math::
vaddpc_f32_1x128_1xf32_u6x_a( const float * __restrict__ vec_x,
                              float       * __restrict__ vec_y,
                              const float                c,
                              std::size_t                n) 
{
        if(__builtin_expect(n==0ull,0)) { return (-1);}
             
        const __m128 vc{_mm_set1_ps(c)};
        std::size_t i;
        for(i = 0ull; (i+23ull) < n; i += 24ull) 
        {
            /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
            */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_store_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_store_ps(&vec_y[i+4ull],vz1);
                vx2 = _mm_load_ps(&vec_x[i+8ull]);
                vy2 = _mm_load_ps(&vec_y[i+8ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_store_ps(&vec_y[i+8ull],vz2);
                vx3 = _mm_load_ps(&vec_x[i+12ull]);
                vy3 = _mm_load_ps(&vec_y[i+12ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_store_ps(&vec_y[i+12ull],vz3);
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_load_ps(&vec_x[i+16ull]);
                vy4 = _mm_load_ps(&vec_y[i+16ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_store_ps(&vec_y[i+16ull],vz4);
                vx5 = _mm_load_ps(&vec_x[i+20ull]);
                vy5 = _mm_load_ps(&vec_y[i+20ull]);
                vz5 = _mm_fmadd_ps(vx5,vc,vy5);
                _mm_store_ps(&vec_y[i+20ull],vz5);
        }

       
        for(; (i+19ull) < n; i += 20ull) 
        {
              /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4; 
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_store_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_store_ps(&vec_y[i+4ull],vz1);
                vx2 = _mm_load_ps(&vec_x[i+8ull]);
                vy2 = _mm_load_ps(&vec_y[i+8ull]);
                vz2 = _mm_fmadd_ps(vx2,vc,vy2);
                _mm_store_ps(&vec_y[i+8ull],vz2);
                vx3 = _mm_load_ps(&vec_x[i+12ull]);
                vy3 = _mm_load_ps(&vec_y[i+12ull]);
                vz3 = _mm_fmadd_ps(vx3,vc,vy3);
                _mm_store_ps(&vec_y[i+12ull],vz3);
                vx4 = _mm_load_ps(&vec_x[i+16ull]);
                vy4 = _mm_load_ps(&vec_y[i+16ull]);
                vz4 = _mm_fmadd_ps(vx4,vc,vy4);
                _mm_store_ps(&vec_y[i+16ull],vz4);
        }

        for(; (i+7ull) < n; i += 8ull) 
        {
             /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                             
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_store_ps(&vec_y[i+0ull],vz0);
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_fmadd_ps(vx1,vc,vy1);
                _mm_store_ps(&vec_y[i+4ull],vz1);
        }

        for(; (i+3ull) < n; i += 4ull) 
        {
             /*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                                             
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_fmadd_ps(vx0,vc,vy0);
                _mm_store_ps(&vec_y[i+0ull],vz0);
        }
        
        for(; (i+0ull) < n; i += 1ull) 
        {
             const float x{vec_x[i]};
             const float y{vec_y[i]};
             const float z{x+y*c};
             vec_y[i] = z;
        }

        return (0);
}


          /*
             Accumulates product of vector and vector
          */

std::int32_t 
gms::math::
vaddpv_f32_2x128_u16x_u(const float * __restrict__ vec_x,
                        const float * __restrict__ vec_y,
                        float       * __restrict__ vec_z,
                        std::size_t                    n) 
{
        if(__builtin_expect(n==0ull,0)) { return (-1);}

        while(((std::uintptr_t)&vec_z & 15) && n) 
        {
             const float x   = *vec_x;
             const float y   = *vec_y;
             const float z   = *vec_z;
             const float res = z+y*x;
             *vec_z          = res;
             vec_x++;
             vec_y++;
             vec_z++;
             n--; 
        }

        std::size_t i;
        for(i = 0ull; (i+63ull) < n; i += 64ull) 
        {
              /*
                VFMADD, Port=0-1, LAT=4, TPUT=1-0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_loadu_ps(&vec_z[i+0ull]);
                _mm_storeu_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_loadu_ps(&vec_z[i+4ull]);
                _mm_storeu_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
                vx2 = _mm_loadu_ps(&vec_x[i+8ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+8ull]);
                vz2 = _mm_loadu_ps(&vec_z[i+8ull]);
                _mm_storeu_ps(&vec_z[i+8ull], _mm_fmadd_ps(vx2,vy2,vz2));
                vx3 = _mm_loadu_ps(&vec_x[i+12ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+12ull]);
                vz3 = _mm_loadu_ps(&vec_z[i+12ull]);
                _mm_storeu_ps(&vec_z[i+12ull], _mm_fmadd_ps(vx3,vy3,vz3));
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_loadu_ps(&vec_x[i+16ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+16ull]);
                vz4 = _mm_loadu_ps(&vec_z[i+16ull]);
                _mm_storeu_ps(&vec_z[i+16ull], _mm_fmadd_ps(vx4,vy4,vz4));
                vx5 = _mm_loadu_ps(&vec_x[i+20ull]);
                vy5 = _mm_loadu_ps(&vec_y[i+20ull]);
                vz5 = _mm_loadu_ps(&vec_z[i+20ull]);
                _mm_storeu_ps(&vec_z[i+20ull], _mm_fmadd_ps(vx5,vy5,vz5));
                vx6 = _mm_loadu_ps(&vec_x[i+24ull]);
                vy6 = _mm_loadu_ps(&vec_y[i+24ull]);
                vz6 = _mm_loadu_ps(&vec_z[i+24ull]);
                _mm_storeu_ps(&vec_z[i+24ull], _mm_fmadd_ps(vx6,vy6,vz6));
                vx7 = _mm_loadu_ps(&vec_x[i+28ull]);
                vy7 = _mm_loadu_ps(&vec_y[i+28ull]);
                vz7 = _mm_loadu_ps(&vec_z[i+28ull]);
                _mm_storeu_ps(&vec_z[i+28ull], _mm_fmadd_ps(vx7,vy7,vz7));
 #if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+32ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+32ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+32ull]);
                vz0 = _mm_loadu_ps(&vec_z[i+32ull]);
                _mm_storeu_ps(&vec_z[i+32ull], _mm_fmadd_ps(vx0,vy0,vz0));  
                vx1 = _mm_loadu_ps(&vec_x[i+36ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+36ull]);
                vz1 = _mm_loadu_ps(&vec_z[i+36ull]);
                _mm_storeu_ps(&vec_z[i+36ull], _mm_fmadd_ps(vx1,vy1,vz1));  
                vx2 = _mm_loadu_ps(&vec_x[i+40ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+40ull]);
                vz2 = _mm_loadu_ps(&vec_z[i+40ull]);
                _mm_storeu_ps(&vec_z[i+40ull], _mm_fmadd_ps(vx2,vy2,vz2)); 
                vx3 = _mm_loadu_ps(&vec_x[i+44ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+44ull]);
                vz3 = _mm_loadu_ps(&vec_z[i+44ull]);
                _mm_storeu_ps(&vec_z[i+44ull], _mm_fmadd_ps(vx3,vy3,vz3));   
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+48ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+48ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_loadu_ps(&vec_x[i+48ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+48ull]);
                vz4 = _mm_loadu_ps(&vec_z[i+48ull]);
                _mm_storeu_ps(&vec_z[i+48ull], _mm_fmadd_ps(vx4,vy4,vz4));   
                vx5 = _mm_loadu_ps(&vec_x[i+52ull]);
                vy5 = _mm_loadu_ps(&vec_y[i+52ull]);
                vz5 = _mm_loadu_ps(&vec_z[i+52ull]);
                _mm_storeu_ps(&vec_z[i+52ull], _mm_fmadd_ps(vx5,vy5,vz5));  
                vx6 = _mm_loadu_ps(&vec_x[i+56ull]);
                vy6 = _mm_loadu_ps(&vec_y[i+56ull]);
                vz6 = _mm_loadu_ps(&vec_z[i+56ull]);
                _mm_storeu_ps(&vec_z[i+56ull], _mm_fmadd_ps(vx6,vy6,vz6)); 
                vx7 = _mm_loadu_ps(&vec_x[i+60ull]);
                vy7 = _mm_loadu_ps(&vec_y[i+60ull]);
                vz7 = _mm_loadu_ps(&vec_z[i+60ull]);
                _mm_storeu_ps(&vec_z[i+60ull], _mm_fmadd_ps(vx7,vy7,vz7)); 

        }

        for(; (i+47ull) < n; i += 48ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;

#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_loadu_ps(&vec_z[i+0ull]);
                _mm_storeu_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_loadu_ps(&vec_z[i+4ull]);
                _mm_storeu_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
                vx2 = _mm_loadu_ps(&vec_x[i+8ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+8ull]);
                vz2 = _mm_loadu_ps(&vec_z[i+8ull]);
                _mm_storeu_ps(&vec_z[i+8ull], _mm_fmadd_ps(vx2,vy2,vz2));
                vx3 = _mm_loadu_ps(&vec_x[i+12ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+12ull]);
                vz3 = _mm_loadu_ps(&vec_z[i+12ull]);
                _mm_storeu_ps(&vec_z[i+12ull], _mm_fmadd_ps(vx3,vy3,vz3));
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_loadu_ps(&vec_x[i+16ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+16ull]);
                vz4 = _mm_loadu_ps(&vec_z[i+16ull]);
                _mm_storeu_ps(&vec_z[i+16ull], _mm_fmadd_ps(vx4,vy4,vz4));
                vx5 = _mm_loadu_ps(&vec_x[i+20ull]);
                vy5 = _mm_loadu_ps(&vec_y[i+20ull]);
                vz5 = _mm_loadu_ps(&vec_z[i+20ull]);
                _mm_storeu_ps(&vec_z[i+20ull], _mm_fmadd_ps(vx5,vy5,vz5));
                vx6 = _mm_loadu_ps(&vec_x[i+24ull]);
                vy6 = _mm_loadu_ps(&vec_y[i+24ull]);
                vz6 = _mm_loadu_ps(&vec_z[i+24ull]);
                _mm_storeu_ps(&vec_z[i+24ull], _mm_fmadd_ps(vx6,vy6,vz6));
                vx7 = _mm_loadu_ps(&vec_x[i+28ull]);
                vy7 = _mm_loadu_ps(&vec_y[i+28ull]);
                vz7 = _mm_loadu_ps(&vec_z[i+28ull]);
                _mm_storeu_ps(&vec_z[i+28ull], _mm_fmadd_ps(vx7,vy7,vz7));
 #if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+32ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+32ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+32ull]);
                vz0 = _mm_loadu_ps(&vec_z[i+32ull]);
                _mm_storeu_ps(&vec_z[i+32ull], _mm_fmadd_ps(vx0,vy0,vz0));  
                vx1 = _mm_loadu_ps(&vec_x[i+36ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+36ull]);
                vz1 = _mm_loadu_ps(&vec_z[i+36ull]);
                _mm_storeu_ps(&vec_z[i+36ull], _mm_fmadd_ps(vx1,vy1,vz1));  
                vx2 = _mm_loadu_ps(&vec_x[i+40ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+40ull]);
                vz2 = _mm_loadu_ps(&vec_z[i+40ull]);
                _mm_storeu_ps(&vec_z[i+40ull], _mm_fmadd_ps(vx2,vy2,vz2)); 
                vx3 = _mm_loadu_ps(&vec_x[i+44ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+44ull]);
                vz3 = _mm_loadu_ps(&vec_z[i+44ull]);
                _mm_storeu_ps(&vec_z[i+44ull], _mm_fmadd_ps(vx3,vy3,vz3)); 
        }

        for(; (i+31ull) < n; i += 32ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;

#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_loadu_ps(&vec_z[i+0ull]);
                _mm_storeu_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_loadu_ps(&vec_z[i+4ull]);
                _mm_storeu_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
                vx2 = _mm_loadu_ps(&vec_x[i+8ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+8ull]);
                vz2 = _mm_loadu_ps(&vec_z[i+8ull]);
                _mm_storeu_ps(&vec_z[i+8ull], _mm_fmadd_ps(vx2,vy2,vz2));
                vx3 = _mm_loadu_ps(&vec_x[i+12ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+12ull]);
                vz3 = _mm_loadu_ps(&vec_z[i+12ull]);
                _mm_storeu_ps(&vec_z[i+12ull], _mm_fmadd_ps(vx3,vy3,vz3));
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_loadu_ps(&vec_x[i+16ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+16ull]);
                vz4 = _mm_loadu_ps(&vec_z[i+16ull]);
                _mm_storeu_ps(&vec_z[i+16ull], _mm_fmadd_ps(vx4,vy4,vz4));
                vx5 = _mm_loadu_ps(&vec_x[i+20ull]);
                vy5 = _mm_loadu_ps(&vec_y[i+20ull]);
                vz5 = _mm_loadu_ps(&vec_z[i+20ull]);
                _mm_storeu_ps(&vec_z[i+20ull], _mm_fmadd_ps(vx5,vy5,vz5));
                vx6 = _mm_loadu_ps(&vec_x[i+24ull]);
                vy6 = _mm_loadu_ps(&vec_y[i+24ull]);
                vz6 = _mm_loadu_ps(&vec_z[i+24ull]);
                _mm_storeu_ps(&vec_z[i+24ull], _mm_fmadd_ps(vx6,vy6,vz6));
                vx7 = _mm_loadu_ps(&vec_x[i+28ull]);
                vy7 = _mm_loadu_ps(&vec_y[i+28ull]);
                vz7 = _mm_loadu_ps(&vec_z[i+28ull]);
                _mm_storeu_ps(&vec_z[i+28ull], _mm_fmadd_ps(vx7,vy7,vz7));
        }

        for(; (i+19ull) < n; i += 20ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_loadu_ps(&vec_z[i+0ull]);
                _mm_storeu_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_loadu_ps(&vec_z[i+4ull]);
                _mm_storeu_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
                vx2 = _mm_loadu_ps(&vec_x[i+8ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+8ull]);
                vz2 = _mm_loadu_ps(&vec_z[i+8ull]);
                _mm_storeu_ps(&vec_z[i+8ull], _mm_fmadd_ps(vx2,vy2,vz2));
                vx3 = _mm_loadu_ps(&vec_x[i+12ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+12ull]);
                vz3 = _mm_loadu_ps(&vec_z[i+12ull]);
                _mm_storeu_ps(&vec_z[i+12ull], _mm_fmadd_ps(vx3,vy3,vz3));
                vx4 = _mm_loadu_ps(&vec_x[i+16ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+16ull]);
                vz4 = _mm_loadu_ps(&vec_z[i+16ull]);
                _mm_storeu_ps(&vec_z[i+16ull], _mm_fmadd_ps(vx4,vy4,vz4));
        }

        for(; (i+7ull) < n; i += 8ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                               
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_loadu_ps(&vec_z[i+0ull]);
                _mm_storeu_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_loadu_ps(&vec_z[i+4ull]);
                _mm_storeu_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
        }

        for(; (i+3ull) < n; i += 4ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                                       
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_loadu_ps(&vec_z[i+0ull]);
                _mm_storeu_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
        }

        for(; (i+0ull) < n; i += 1ull) 
        {
            const float x{vec_x[i]};
            const float y{vec_y[i]};
            const float z{vec_z[i]};
            const float res{z+x*y};
            vec_z[i]  = res;
        }

        return (0);
}


std::int32_t 
gms::math::
vaddpv_f32_2x128_u16x_a(const float * __restrict__ vec_x,
                        const float * __restrict__ vec_y,
                        float       * __restrict__ vec_z,
                        std::size_t                    n) 
{
        if(__builtin_expect(n==0ull,0)) { return (-1);}
        
        std::size_t i;
        for(i = 0ull; (i+63ull) < n; i += 64ull) 
        {
              /*
                VFMADD, Port=0-1, LAT=4, TPUT=1-0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_load_ps(&vec_z[i+0ull]);
                _mm_store_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_load_ps(&vec_z[i+4ull]);
                _mm_store_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
                vx2 = _mm_load_ps(&vec_x[i+8ull]);
                vy2 = _mm_load_ps(&vec_y[i+8ull]);
                vz2 = _mm_load_ps(&vec_z[i+8ull]);
                _mm_store_ps(&vec_z[i+8ull], _mm_fmadd_ps(vx2,vy2,vz2));
                vx3 = _mm_load_ps(&vec_x[i+12ull]);
                vy3 = _mm_load_ps(&vec_y[i+12ull]);
                vz3 = _mm_load_ps(&vec_z[i+12ull]);
                _mm_store_ps(&vec_z[i+12ull], _mm_fmadd_ps(vx3,vy3,vz3));
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_load_ps(&vec_x[i+16ull]);
                vy4 = _mm_load_ps(&vec_y[i+16ull]);
                vz4 = _mm_load_ps(&vec_z[i+16ull]);
                _mm_store_ps(&vec_z[i+16ull], _mm_fmadd_ps(vx4,vy4,vz4));
                vx5 = _mm_load_ps(&vec_x[i+20ull]);
                vy5 = _mm_load_ps(&vec_y[i+20ull]);
                vz5 = _mm_load_ps(&vec_z[i+20ull]);
                _mm_store_ps(&vec_z[i+20ull], _mm_fmadd_ps(vx5,vy5,vz5));
                vx6 = _mm_load_ps(&vec_x[i+24ull]);
                vy6 = _mm_load_ps(&vec_y[i+24ull]);
                vz6 = _mm_load_ps(&vec_z[i+24ull]);
                _mm_store_ps(&vec_z[i+24ull], _mm_fmadd_ps(vx6,vy6,vz6));
                vx7 = _mm_load_ps(&vec_x[i+28ull]);
                vy7 = _mm_load_ps(&vec_y[i+28ull]);
                vz7 = _mm_load_ps(&vec_z[i+28ull]);
                _mm_store_ps(&vec_z[i+28ull], _mm_fmadd_ps(vx7,vy7,vz7));
 #if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+32ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+32ull]);
                vy0 = _mm_load_ps(&vec_y[i+32ull]);
                vz0 = _mm_load_ps(&vec_z[i+32ull]);
                _mm_store_ps(&vec_z[i+32ull], _mm_fmadd_ps(vx0,vy0,vz0));  
                vx1 = _mm_load_ps(&vec_x[i+36ull]);
                vy1 = _mm_load_ps(&vec_y[i+36ull]);
                vz1 = _mm_load_ps(&vec_z[i+36ull]);
                _mm_store_ps(&vec_z[i+36ull], _mm_fmadd_ps(vx1,vy1,vz1));  
                vx2 = _mm_load_ps(&vec_x[i+40ull]);
                vy2 = _mm_load_ps(&vec_y[i+40ull]);
                vz2 = _mm_load_ps(&vec_z[i+40ull]);
                _mm_store_ps(&vec_z[i+40ull], _mm_fmadd_ps(vx2,vy2,vz2)); 
                vx3 = _mm_load_ps(&vec_x[i+44ull]);
                vy3 = _mm_load_ps(&vec_y[i+44ull]);
                vz3 = _mm_load_ps(&vec_z[i+44ull]);
                _mm_store_ps(&vec_z[i+44ull], _mm_fmadd_ps(vx3,vy3,vz3));   
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+48ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+48ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_load_ps(&vec_x[i+48ull]);
                vy4 = _mm_load_ps(&vec_y[i+48ull]);
                vz4 = _mm_load_ps(&vec_z[i+48ull]);
                _mm_store_ps(&vec_z[i+48ull], _mm_fmadd_ps(vx4,vy4,vz4));   
                vx5 = _mm_load_ps(&vec_x[i+52ull]);
                vy5 = _mm_load_ps(&vec_y[i+52ull]);
                vz5 = _mm_load_ps(&vec_z[i+52ull]);
                _mm_store_ps(&vec_z[i+52ull], _mm_fmadd_ps(vx5,vy5,vz5));  
                vx6 = _mm_load_ps(&vec_x[i+56ull]);
                vy6 = _mm_load_ps(&vec_y[i+56ull]);
                vz6 = _mm_load_ps(&vec_z[i+56ull]);
                _mm_store_ps(&vec_z[i+56ull], _mm_fmadd_ps(vx6,vy6,vz6)); 
                vx7 = _mm_load_ps(&vec_x[i+60ull]);
                vy7 = _mm_load_ps(&vec_y[i+60ull]);
                vz7 = _mm_load_ps(&vec_z[i+60ull]);
                _mm_store_ps(&vec_z[i+60ull], _mm_fmadd_ps(vx7,vy7,vz7)); 

        }

        for(; (i+47ull) < n; i += 48ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;

#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_load_ps(&vec_z[i+0ull]);
                _mm_store_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_load_ps(&vec_z[i+4ull]);
                _mm_store_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
                vx2 = _mm_load_ps(&vec_x[i+8ull]);
                vy2 = _mm_load_ps(&vec_y[i+8ull]);
                vz2 = _mm_load_ps(&vec_z[i+8ull]);
                _mm_store_ps(&vec_z[i+8ull], _mm_fmadd_ps(vx2,vy2,vz2));
                vx3 = _mm_load_ps(&vec_x[i+12ull]);
                vy3 = _mm_load_ps(&vec_y[i+12ull]);
                vz3 = _mm_load_ps(&vec_z[i+12ull]);
                _mm_store_ps(&vec_z[i+12ull], _mm_fmadd_ps(vx3,vy3,vz3));
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_load_ps(&vec_x[i+16ull]);
                vy4 = _mm_load_ps(&vec_y[i+16ull]);
                vz4 = _mm_load_ps(&vec_z[i+16ull]);
                _mm_store_ps(&vec_z[i+16ull], _mm_fmadd_ps(vx4,vy4,vz4));
                vx5 = _mm_load_ps(&vec_x[i+20ull]);
                vy5 = _mm_load_ps(&vec_y[i+20ull]);
                vz5 = _mm_load_ps(&vec_z[i+20ull]);
                _mm_store_ps(&vec_z[i+20ull], _mm_fmadd_ps(vx5,vy5,vz5));
                vx6 = _mm_load_ps(&vec_x[i+24ull]);
                vy6 = _mm_load_ps(&vec_y[i+24ull]);
                vz6 = _mm_load_ps(&vec_z[i+24ull]);
                _mm_store_ps(&vec_z[i+24ull], _mm_fmadd_ps(vx6,vy6,vz6));
                vx7 = _mm_load_ps(&vec_x[i+28ull]);
                vy7 = _mm_load_ps(&vec_y[i+28ull]);
                vz7 = _mm_load_ps(&vec_z[i+28ull]);
                _mm_store_ps(&vec_z[i+28ull], _mm_fmadd_ps(vx7,vy7,vz7));
 #if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+32ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+32ull]);
                vy0 = _mm_load_ps(&vec_y[i+32ull]);
                vz0 = _mm_load_ps(&vec_z[i+32ull]);
                _mm_store_ps(&vec_z[i+32ull], _mm_fmadd_ps(vx0,vy0,vz0));  
                vx1 = _mm_load_ps(&vec_x[i+36ull]);
                vy1 = _mm_load_ps(&vec_y[i+36ull]);
                vz1 = _mm_load_ps(&vec_z[i+36ull]);
                _mm_store_ps(&vec_z[i+36ull], _mm_fmadd_ps(vx1,vy1,vz1));  
                vx2 = _mm_load_ps(&vec_x[i+40ull]);
                vy2 = _mm_load_ps(&vec_y[i+40ull]);
                vz2 = _mm_load_ps(&vec_z[i+40ull]);
                _mm_store_ps(&vec_z[i+40ull], _mm_fmadd_ps(vx2,vy2,vz2)); 
                vx3 = _mm_load_ps(&vec_x[i+44ull]);
                vy3 = _mm_load_ps(&vec_y[i+44ull]);
                vz3 = _mm_load_ps(&vec_z[i+44ull]);
                _mm_store_ps(&vec_z[i+44ull], _mm_fmadd_ps(vx3,vy3,vz3)); 
        }

        for(; (i+31ull) < n; i += 32ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;

#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_load_ps(&vec_z[i+0ull]);
                _mm_store_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_load_ps(&vec_z[i+4ull]);
                _mm_store_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
                vx2 = _mm_load_ps(&vec_x[i+8ull]);
                vy2 = _mm_load_ps(&vec_y[i+8ull]);
                vz2 = _mm_load_ps(&vec_z[i+8ull]);
                _mm_store_ps(&vec_z[i+8ull], _mm_fmadd_ps(vx2,vy2,vz2));
                vx3 = _mm_load_ps(&vec_x[i+12ull]);
                vy3 = _mm_load_ps(&vec_y[i+12ull]);
                vz3 = _mm_load_ps(&vec_z[i+12ull]);
                _mm_store_ps(&vec_z[i+12ull], _mm_fmadd_ps(vx3,vy3,vz3));
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_load_ps(&vec_x[i+16ull]);
                vy4 = _mm_load_ps(&vec_y[i+16ull]);
                vz4 = _mm_load_ps(&vec_z[i+16ull]);
                _mm_store_ps(&vec_z[i+16ull], _mm_fmadd_ps(vx4,vy4,vz4));
                vx5 = _mm_load_ps(&vec_x[i+20ull]);
                vy5 = _mm_load_ps(&vec_y[i+20ull]);
                vz5 = _mm_load_ps(&vec_z[i+20ull]);
                _mm_store_ps(&vec_z[i+20ull], _mm_fmadd_ps(vx5,vy5,vz5));
                vx6 = _mm_load_ps(&vec_x[i+24ull]);
                vy6 = _mm_load_ps(&vec_y[i+24ull]);
                vz6 = _mm_load_ps(&vec_z[i+24ull]);
                _mm_store_ps(&vec_z[i+24ull], _mm_fmadd_ps(vx6,vy6,vz6));
                vx7 = _mm_load_ps(&vec_x[i+28ull]);
                vy7 = _mm_load_ps(&vec_y[i+28ull]);
                vz7 = _mm_load_ps(&vec_z[i+28ull]);
                _mm_store_ps(&vec_z[i+28ull], _mm_fmadd_ps(vx7,vy7,vz7));
        }

        for(; (i+19ull) < n; i += 20ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_load_ps(&vec_z[i+0ull]);
                _mm_store_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_load_ps(&vec_z[i+4ull]);
                _mm_store_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
                vx2 = _mm_load_ps(&vec_x[i+8ull]);
                vy2 = _mm_load_ps(&vec_y[i+8ull]);
                vz2 = _mm_load_ps(&vec_z[i+8ull]);
                _mm_store_ps(&vec_z[i+8ull], _mm_fmadd_ps(vx2,vy2,vz2));
                vx3 = _mm_load_ps(&vec_x[i+12ull]);
                vy3 = _mm_load_ps(&vec_y[i+12ull]);
                vz3 = _mm_load_ps(&vec_z[i+12ull]);
                _mm_store_ps(&vec_z[i+12ull], _mm_fmadd_ps(vx3,vy3,vz3));
                vx4 = _mm_load_ps(&vec_x[i+16ull]);
                vy4 = _mm_load_ps(&vec_y[i+16ull]);
                vz4 = _mm_load_ps(&vec_z[i+16ull]);
                _mm_store_ps(&vec_z[i+16ull], _mm_fmadd_ps(vx4,vy4,vz4));
        }

        for(; (i+7ull) < n; i += 8ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                               
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_load_ps(&vec_z[i+0ull]);
                _mm_store_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_load_ps(&vec_z[i+4ull]);
                _mm_store_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
        }

        for(; (i+3ull) < n; i += 4ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                                       
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_load_ps(&vec_z[i+0ull]);
                _mm_store_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
        }

        for(; (i+0ull) < n; i += 1ull) 
        {
            const float x{vec_x[i]};
            const float y{vec_y[i]};
            const float z{vec_z[i]};
            const float res{z+x*y};
            vec_z[i]  = res;
        }

        return (0);
}


std::int32_t 
gms::math::
vaddpv_f32_2x128_u10x_u(const float * __restrict__ vec_x,
                        const float * __restrict__ vec_y,
                        float       * __restrict__ vec_z,
                        std::size_t                    n) 
{
        if(__builtin_expect(n==0ull,0)) { return (-1);}

        while(((std::uintptr_t)&vec_z & 15) && n) 
        {
             const float x   = *vec_x;
             const float y   = *vec_y;
             const float z   = *vec_z;
             const float res = z+y*x;
             *vec_z          = res;
             vec_x++;
             vec_y++;
             vec_z++;
             n--; 
        }

        std::size_t i;
        for(i = 0ull; (i+39ull) < n; i += 40ull) 
        {
              /*
                VFMADD, Port=0-1, LAT=4, TPUT=1-0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_loadu_ps(&vec_z[i+0ull]);
                _mm_storeu_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_loadu_ps(&vec_z[i+4ull]);
                _mm_storeu_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
                vx2 = _mm_loadu_ps(&vec_x[i+8ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+8ull]);
                vz2 = _mm_loadu_ps(&vec_z[i+8ull]);
                _mm_storeu_ps(&vec_z[i+8ull], _mm_fmadd_ps(vx2,vy2,vz2));
                vx3 = _mm_loadu_ps(&vec_x[i+12ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+12ull]);
                vz3 = _mm_loadu_ps(&vec_z[i+12ull]);
                _mm_storeu_ps(&vec_z[i+12ull], _mm_fmadd_ps(vx3,vy3,vz3));
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_loadu_ps(&vec_x[i+16ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+16ull]);
                vz4 = _mm_loadu_ps(&vec_z[i+16ull]);
                _mm_storeu_ps(&vec_z[i+16ull], _mm_fmadd_ps(vx4,vy4,vz4));
                vx5 = _mm_loadu_ps(&vec_x[i+20ull]);
                vy5 = _mm_loadu_ps(&vec_y[i+20ull]);
                vz5 = _mm_loadu_ps(&vec_z[i+20ull]);
                _mm_storeu_ps(&vec_z[i+20ull], _mm_fmadd_ps(vx5,vy5,vz5));
                vx6 = _mm_loadu_ps(&vec_x[i+24ull]);
                vy6 = _mm_loadu_ps(&vec_y[i+24ull]);
                vz6 = _mm_loadu_ps(&vec_z[i+24ull]);
                _mm_storeu_ps(&vec_z[i+24ull], _mm_fmadd_ps(vx6,vy6,vz6));
                vx7 = _mm_loadu_ps(&vec_x[i+28ull]);
                vy7 = _mm_loadu_ps(&vec_y[i+28ull]);
                vz7 = _mm_loadu_ps(&vec_z[i+28ull]);
                _mm_storeu_ps(&vec_z[i+28ull], _mm_fmadd_ps(vx7,vy7,vz7));
 #if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+32ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+32ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+32ull]);
                vz0 = _mm_loadu_ps(&vec_z[i+32ull]);
                _mm_storeu_ps(&vec_z[i+32ull], _mm_fmadd_ps(vx0,vy0,vz0));  
                vx1 = _mm_loadu_ps(&vec_x[i+36ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+36ull]);
                vz1 = _mm_loadu_ps(&vec_z[i+36ull]);
                _mm_storeu_ps(&vec_z[i+36ull], _mm_fmadd_ps(vx1,vy1,vz1));  
        }

        for(; (i+31ull) < n; i += 32ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;

#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_loadu_ps(&vec_z[i+0ull]);
                _mm_storeu_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_loadu_ps(&vec_z[i+4ull]);
                _mm_storeu_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
                vx2 = _mm_loadu_ps(&vec_x[i+8ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+8ull]);
                vz2 = _mm_loadu_ps(&vec_z[i+8ull]);
                _mm_storeu_ps(&vec_z[i+8ull], _mm_fmadd_ps(vx2,vy2,vz2));
                vx3 = _mm_loadu_ps(&vec_x[i+12ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+12ull]);
                vz3 = _mm_loadu_ps(&vec_z[i+12ull]);
                _mm_storeu_ps(&vec_z[i+12ull], _mm_fmadd_ps(vx3,vy3,vz3));
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_loadu_ps(&vec_x[i+16ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+16ull]);
                vz4 = _mm_loadu_ps(&vec_z[i+16ull]);
                _mm_storeu_ps(&vec_z[i+16ull], _mm_fmadd_ps(vx4,vy4,vz4));
                vx5 = _mm_loadu_ps(&vec_x[i+20ull]);
                vy5 = _mm_loadu_ps(&vec_y[i+20ull]);
                vz5 = _mm_loadu_ps(&vec_z[i+20ull]);
                _mm_storeu_ps(&vec_z[i+20ull], _mm_fmadd_ps(vx5,vy5,vz5));
                vx6 = _mm_loadu_ps(&vec_x[i+24ull]);
                vy6 = _mm_loadu_ps(&vec_y[i+24ull]);
                vz6 = _mm_loadu_ps(&vec_z[i+24ull]);
                _mm_storeu_ps(&vec_z[i+24ull], _mm_fmadd_ps(vx6,vy6,vz6));
                vx7 = _mm_loadu_ps(&vec_x[i+28ull]);
                vy7 = _mm_loadu_ps(&vec_y[i+28ull]);
                vz7 = _mm_loadu_ps(&vec_z[i+28ull]);
                _mm_storeu_ps(&vec_z[i+28ull], _mm_fmadd_ps(vx7,vy7,vz7));
        }

        for(; (i+19ull) < n; i += 20ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_loadu_ps(&vec_z[i+0ull]);
                _mm_storeu_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_loadu_ps(&vec_z[i+4ull]);
                _mm_storeu_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
                vx2 = _mm_loadu_ps(&vec_x[i+8ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+8ull]);
                vz2 = _mm_loadu_ps(&vec_z[i+8ull]);
                _mm_storeu_ps(&vec_z[i+8ull], _mm_fmadd_ps(vx2,vy2,vz2));
                vx3 = _mm_loadu_ps(&vec_x[i+12ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+12ull]);
                vz3 = _mm_loadu_ps(&vec_z[i+12ull]);
                _mm_storeu_ps(&vec_z[i+12ull], _mm_fmadd_ps(vx3,vy3,vz3));
                vx4 = _mm_loadu_ps(&vec_x[i+16ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+16ull]);
                vz4 = _mm_loadu_ps(&vec_z[i+16ull]);
                _mm_storeu_ps(&vec_z[i+16ull], _mm_fmadd_ps(vx4,vy4,vz4));
        }

        for(; (i+7ull) < n; i += 8ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                               
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_loadu_ps(&vec_z[i+0ull]);
                _mm_storeu_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_loadu_ps(&vec_z[i+4ull]);
                _mm_storeu_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
        }

        for(; (i+3ull) < n; i += 4ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                                       
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_loadu_ps(&vec_z[i+0ull]);
                _mm_storeu_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
        }

        for(; (i+0ull) < n; i += 1ull) 
        {
            const float x{vec_x[i]};
            const float y{vec_y[i]};
            const float z{vec_z[i]};
            const float res{z+x*y};
            vec_z[i]  = res;
        }

        return (0);
}


std::int32_t 
gms::math::
vaddpv_f32_2x128_u10x_a(const float * __restrict__ vec_x,
                        const float * __restrict__ vec_y,
                        float       * __restrict__ vec_z,
                        std::size_t                    n) 
{
        if(__builtin_expect(n==0ull,0)) { return (-1);}
        
        std::size_t i;
        for(i = 0ull; (i+39ull) < n; i += 40ull) 
        {
              /*
                VFMADD, Port=0-1, LAT=4, TPUT=1-0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_load_ps(&vec_z[i+0ull]);
                _mm_store_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_load_ps(&vec_z[i+4ull]);
                _mm_store_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
                vx2 = _mm_load_ps(&vec_x[i+8ull]);
                vy2 = _mm_load_ps(&vec_y[i+8ull]);
                vz2 = _mm_load_ps(&vec_z[i+8ull]);
                _mm_store_ps(&vec_z[i+8ull], _mm_fmadd_ps(vx2,vy2,vz2));
                vx3 = _mm_load_ps(&vec_x[i+12ull]);
                vy3 = _mm_load_ps(&vec_y[i+12ull]);
                vz3 = _mm_load_ps(&vec_z[i+12ull]);
                _mm_store_ps(&vec_z[i+12ull], _mm_fmadd_ps(vx3,vy3,vz3));
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_load_ps(&vec_x[i+16ull]);
                vy4 = _mm_load_ps(&vec_y[i+16ull]);
                vz4 = _mm_load_ps(&vec_z[i+16ull]);
                _mm_store_ps(&vec_z[i+16ull], _mm_fmadd_ps(vx4,vy4,vz4));
                vx5 = _mm_load_ps(&vec_x[i+20ull]);
                vy5 = _mm_load_ps(&vec_y[i+20ull]);
                vz5 = _mm_load_ps(&vec_z[i+20ull]);
                _mm_store_ps(&vec_z[i+20ull], _mm_fmadd_ps(vx5,vy5,vz5));
                vx6 = _mm_load_ps(&vec_x[i+24ull]);
                vy6 = _mm_load_ps(&vec_y[i+24ull]);
                vz6 = _mm_load_ps(&vec_z[i+24ull]);
                _mm_store_ps(&vec_z[i+24ull], _mm_fmadd_ps(vx6,vy6,vz6));
                vx7 = _mm_load_ps(&vec_x[i+28ull]);
                vy7 = _mm_load_ps(&vec_y[i+28ull]);
                vz7 = _mm_load_ps(&vec_z[i+28ull]);
                _mm_store_ps(&vec_z[i+28ull], _mm_fmadd_ps(vx7,vy7,vz7));
 #if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+32ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+32ull]);
                vy0 = _mm_load_ps(&vec_y[i+32ull]);
                vz0 = _mm_load_ps(&vec_z[i+32ull]);
                _mm_store_ps(&vec_z[i+32ull], _mm_fmadd_ps(vx0,vy0,vz0));  
                vx1 = _mm_load_ps(&vec_x[i+36ull]);
                vy1 = _mm_load_ps(&vec_y[i+36ull]);
                vz1 = _mm_load_ps(&vec_z[i+36ull]);
                _mm_store_ps(&vec_z[i+36ull], _mm_fmadd_ps(vx1,vy1,vz1));  
        }

        for(; (i+31ull) < n; i += 32ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;

#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_load_ps(&vec_z[i+0ull]);
                _mm_store_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_load_ps(&vec_z[i+4ull]);
                _mm_store_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
                vx2 = _mm_load_ps(&vec_x[i+8ull]);
                vy2 = _mm_load_ps(&vec_y[i+8ull]);
                vz2 = _mm_load_ps(&vec_z[i+8ull]);
                _mm_store_ps(&vec_z[i+8ull], _mm_fmadd_ps(vx2,vy2,vz2));
                vx3 = _mm_load_ps(&vec_x[i+12ull]);
                vy3 = _mm_load_ps(&vec_y[i+12ull]);
                vz3 = _mm_load_ps(&vec_z[i+12ull]);
                _mm_store_ps(&vec_z[i+12ull], _mm_fmadd_ps(vx3,vy3,vz3));
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+16ull],_MM_HINT_T0);
#endif 
                vx4 = _mm_load_ps(&vec_x[i+16ull]);
                vy4 = _mm_load_ps(&vec_y[i+16ull]);
                vz4 = _mm_load_ps(&vec_z[i+16ull]);
                _mm_store_ps(&vec_z[i+16ull], _mm_fmadd_ps(vx4,vy4,vz4));
                vx5 = _mm_load_ps(&vec_x[i+20ull]);
                vy5 = _mm_load_ps(&vec_y[i+20ull]);
                vz5 = _mm_load_ps(&vec_z[i+20ull]);
                _mm_store_ps(&vec_z[i+20ull], _mm_fmadd_ps(vx5,vy5,vz5));
                vx6 = _mm_load_ps(&vec_x[i+24ull]);
                vy6 = _mm_load_ps(&vec_y[i+24ull]);
                vz6 = _mm_load_ps(&vec_z[i+24ull]);
                _mm_store_ps(&vec_z[i+24ull], _mm_fmadd_ps(vx6,vy6,vz6));
                vx7 = _mm_load_ps(&vec_x[i+28ull]);
                vy7 = _mm_load_ps(&vec_y[i+28ull]);
                vz7 = _mm_load_ps(&vec_z[i+28ull]);
                _mm_store_ps(&vec_z[i+28ull], _mm_fmadd_ps(vx7,vy7,vz7));
        }

        for(; (i+19ull) < n; i += 20ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_load_ps(&vec_z[i+0ull]);
                _mm_store_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_load_ps(&vec_z[i+4ull]);
                _mm_store_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
                vx2 = _mm_load_ps(&vec_x[i+8ull]);
                vy2 = _mm_load_ps(&vec_y[i+8ull]);
                vz2 = _mm_load_ps(&vec_z[i+8ull]);
                _mm_store_ps(&vec_z[i+8ull], _mm_fmadd_ps(vx2,vy2,vz2));
                vx3 = _mm_load_ps(&vec_x[i+12ull]);
                vy3 = _mm_load_ps(&vec_y[i+12ull]);
                vz3 = _mm_load_ps(&vec_z[i+12ull]);
                _mm_store_ps(&vec_z[i+12ull], _mm_fmadd_ps(vx3,vy3,vz3));
                vx4 = _mm_load_ps(&vec_x[i+16ull]);
                vy4 = _mm_load_ps(&vec_y[i+16ull]);
                vz4 = _mm_load_ps(&vec_z[i+16ull]);
                _mm_store_ps(&vec_z[i+16ull], _mm_fmadd_ps(vx4,vy4,vz4));
        }

        for(; (i+7ull) < n; i += 8ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                               
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_load_ps(&vec_z[i+0ull]);
                _mm_store_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_load_ps(&vec_z[i+4ull]);
                _mm_store_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
        }

        for(; (i+3ull) < n; i += 4ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                                       
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_load_ps(&vec_z[i+0ull]);
                _mm_store_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
        }

        for(; (i+0ull) < n; i += 1ull) 
        {
            const float x{vec_x[i]};
            const float y{vec_y[i]};
            const float z{vec_z[i]};
            const float res{z+x*y};
            vec_z[i]  = res;
        }

        return (0);
}


std::int32_t 
gms::math::
vaddpv_f32_2x128_u6x_u(const float * __restrict__ vec_x,
                        const float * __restrict__ vec_y,
                        float       * __restrict__ vec_z,
                        std::size_t                    n) 
{
        if(__builtin_expect(n==0ull,0)) { return (-1);}

        while(((std::uintptr_t)&vec_z & 15) && n) 
        {
             const float x   = *vec_x;
             const float y   = *vec_y;
             const float z   = *vec_z;
             const float res = z+y*x;
             *vec_z          = res;
             vec_x++;
             vec_y++;
             vec_z++;
             n--; 
        }

        std::size_t i;
        for(i = 0ull; (i+23ull) < n; i += 24ull) 
        {
              /*
                VFMADD, Port=0-1, LAT=4, TPUT=1-0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_loadu_ps(&vec_z[i+0ull]);
                _mm_storeu_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_loadu_ps(&vec_z[i+4ull]);
                _mm_storeu_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
                vx2 = _mm_loadu_ps(&vec_x[i+8ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+8ull]);
                vz2 = _mm_loadu_ps(&vec_z[i+8ull]);
                _mm_storeu_ps(&vec_z[i+8ull], _mm_fmadd_ps(vx2,vy2,vz2));
                vx3 = _mm_loadu_ps(&vec_x[i+12ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+12ull]);
                vz3 = _mm_loadu_ps(&vec_z[i+12ull]);
                _mm_storeu_ps(&vec_z[i+12ull], _mm_fmadd_ps(vx3,vy3,vz3));
                vx4 = _mm_loadu_ps(&vec_x[i+16ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+16ull]);
                vz4 = _mm_loadu_ps(&vec_z[i+16ull]);
                _mm_storeu_ps(&vec_z[i+16ull], _mm_fmadd_ps(vx4,vy4,vz4));
                vx5 = _mm_loadu_ps(&vec_x[i+20ull]);
                vy5 = _mm_loadu_ps(&vec_y[i+20ull]);
                vz5 = _mm_loadu_ps(&vec_z[i+20ull]);
                _mm_storeu_ps(&vec_z[i+20ull], _mm_fmadd_ps(vx5,vy5,vz5));
        }

        for(; (i+19ull) < n; i += 20ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_loadu_ps(&vec_z[i+0ull]);
                _mm_storeu_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_loadu_ps(&vec_z[i+4ull]);
                _mm_storeu_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
                vx2 = _mm_loadu_ps(&vec_x[i+8ull]);
                vy2 = _mm_loadu_ps(&vec_y[i+8ull]);
                vz2 = _mm_loadu_ps(&vec_z[i+8ull]);
                _mm_storeu_ps(&vec_z[i+8ull], _mm_fmadd_ps(vx2,vy2,vz2));
                vx3 = _mm_loadu_ps(&vec_x[i+12ull]);
                vy3 = _mm_loadu_ps(&vec_y[i+12ull]);
                vz3 = _mm_loadu_ps(&vec_z[i+12ull]);
                _mm_storeu_ps(&vec_z[i+12ull], _mm_fmadd_ps(vx3,vy3,vz3));
                vx4 = _mm_loadu_ps(&vec_x[i+16ull]);
                vy4 = _mm_loadu_ps(&vec_y[i+16ull]);
                vz4 = _mm_loadu_ps(&vec_z[i+16ull]);
                _mm_storeu_ps(&vec_z[i+16ull], _mm_fmadd_ps(vx4,vy4,vz4));
        }

        for(; (i+7ull) < n; i += 8ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                               
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_loadu_ps(&vec_z[i+0ull]);
                _mm_storeu_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_loadu_ps(&vec_x[i+4ull]);
                vy1 = _mm_loadu_ps(&vec_y[i+4ull]);
                vz1 = _mm_loadu_ps(&vec_z[i+4ull]);
                _mm_storeu_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
        }

        for(; (i+3ull) < n; i += 4ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                                       
                vx0 = _mm_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm_loadu_ps(&vec_z[i+0ull]);
                _mm_storeu_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
        }
        for(; (i+0ull) < n; i += 1ull) 
        {
            const float x{vec_x[i]};
            const float y{vec_y[i]};
            const float z{vec_z[i]};
            const float res{z+x*y};
            vec_z[i]  = res;
        }

        return (0);
}


std::int32_t 
gms::math::
vaddpv_f32_2x128_u6x_a( const float * __restrict__ vec_x,
                        const float * __restrict__ vec_y,
                        float       * __restrict__ vec_z,
                        std::size_t                    n) 
{
        if(__builtin_expect(n==0ull,0)) { return (-1);}
        
        std::size_t i;
        for(i = 0ull; (i+23ull) < n; i += 24ull) 
        {
              /*
                VFMADD, Port=0-1, LAT=4, TPUT=1-0.5
               */
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                __m128 vx5;
                __m128 vy5;
                __m128 vz5;
                __m128 vx6;
                __m128 vy6;
                __m128 vz6;
                __m128 vx7;
                __m128 vy7;
                __m128 vz7;
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_load_ps(&vec_z[i+0ull]);
                _mm_store_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_load_ps(&vec_z[i+4ull]);
                _mm_store_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
                vx2 = _mm_load_ps(&vec_x[i+8ull]);
                vy2 = _mm_load_ps(&vec_y[i+8ull]);
                vz2 = _mm_load_ps(&vec_z[i+8ull]);
                _mm_store_ps(&vec_z[i+8ull], _mm_fmadd_ps(vx2,vy2,vz2));
                vx3 = _mm_load_ps(&vec_x[i+12ull]);
                vy3 = _mm_load_ps(&vec_y[i+12ull]);
                vz3 = _mm_load_ps(&vec_z[i+12ull]);
                _mm_store_ps(&vec_z[i+12ull], _mm_fmadd_ps(vx3,vy3,vz3));
                vx4 = _mm_load_ps(&vec_x[i+16ull]);
                vy4 = _mm_load_ps(&vec_y[i+16ull]);
                vz4 = _mm_load_ps(&vec_z[i+16ull]);
                _mm_store_ps(&vec_z[i+16ull], _mm_fmadd_ps(vx4,vy4,vz4));
                vx5 = _mm_load_ps(&vec_x[i+20ull]);
                vy5 = _mm_load_ps(&vec_y[i+20ull]);
                vz5 = _mm_load_ps(&vec_z[i+20ull]);
                _mm_store_ps(&vec_z[i+20ull], _mm_fmadd_ps(vx5,vy5,vz5));
        }

        for(; (i+19ull) < n; i += 20ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                __m128 vx2;
                __m128 vy2;
                __m128 vz2;
                __m128 vx3;
                __m128 vy3;
                __m128 vz3;
                __m128 vx4;
                __m128 vy4;
                __m128 vz4;
                
#if (ADD_PRODUCT_VEC_XMM4R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_z[i+0ull],_MM_HINT_T0);
#endif 
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_load_ps(&vec_z[i+0ull]);
                _mm_store_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_load_ps(&vec_z[i+4ull]);
                _mm_store_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
                vx2 = _mm_load_ps(&vec_x[i+8ull]);
                vy2 = _mm_load_ps(&vec_y[i+8ull]);
                vz2 = _mm_load_ps(&vec_z[i+8ull]);
                _mm_store_ps(&vec_z[i+8ull], _mm_fmadd_ps(vx2,vy2,vz2));
                vx3 = _mm_load_ps(&vec_x[i+12ull]);
                vy3 = _mm_load_ps(&vec_y[i+12ull]);
                vz3 = _mm_load_ps(&vec_z[i+12ull]);
                _mm_store_ps(&vec_z[i+12ull], _mm_fmadd_ps(vx3,vy3,vz3));
                vx4 = _mm_load_ps(&vec_x[i+16ull]);
                vy4 = _mm_load_ps(&vec_y[i+16ull]);
                vz4 = _mm_load_ps(&vec_z[i+16ull]);
                _mm_store_ps(&vec_z[i+16ull], _mm_fmadd_ps(vx4,vy4,vz4));
        }

        for(; (i+7ull) < n; i += 8ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                __m128 vx1;
                __m128 vy1;
                __m128 vz1;
                               
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_load_ps(&vec_z[i+0ull]);
                _mm_store_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
                vx1 = _mm_load_ps(&vec_x[i+4ull]);
                vy1 = _mm_load_ps(&vec_y[i+4ull]);
                vz1 = _mm_load_ps(&vec_z[i+4ull]);
                _mm_store_ps(&vec_z[i+4ull], _mm_fmadd_ps(vx1,vy1,vz1)); 
        }

        for(; (i+3ull) < n; i += 4ull) 
        {
                __m128 vx0;
                __m128 vy0;
                __m128 vz0;
                                       
                vx0 = _mm_load_ps(&vec_x[i+0ull]);
                vy0 = _mm_load_ps(&vec_y[i+0ull]);
                vz0 = _mm_load_ps(&vec_z[i+0ull]);
                _mm_store_ps(&vec_z[i+0ull], _mm_fmadd_ps(vx0,vy0,vz0));
        }

        for(; (i+0ull) < n; i += 1ull) 
        {
            const float x{vec_x[i]};
            const float y{vec_y[i]};
            const float z{vec_z[i]};
            const float res{z+x*y};
            vec_z[i]  = res;
        }

        return (0);
}







