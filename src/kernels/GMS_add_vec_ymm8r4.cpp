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
#include "GMS_add_vec_ymm8r4.h"

std::int32_t 
gms::math::vadd_f32_2x256v_u16x_u(const float * __restrict__ vec_x,
                                  const float * __restrict__ vec_y,
                                  float       * __restrict__ vec_z,
                                  std::size_t          n)
{
            if(__builtin_expect(n==0ull,0)) { return -1;}

#if (ADD_VEC_YMM8R4_ADD_PEEL_LOOP) == 1
            while(((std::uintptr_t)&vec_z & 31) && n)
            {
                   const float x = *vec_x;
                   const float y = *vec_y;
                   const float z = x+y;
                   *vec_z        = z;
                   vec_x++;
                   vec_y++;
                   vec_z++;
                   n--;
            }
#endif 

            std::size_t i;
            for(i = 0ull; (i+127ull) < n; i += 128ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;
                __m256 vx4;
                __m256 vy4;
                __m256 vz4;
                __m256 vx5;
                __m256 vy5;
                __m256 vz5;
                __m256 vx6;
                __m256 vy6;
                __m256 vz6;
                __m256 vx7;
                __m256 vy7;
                __m256 vz7;

#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                vy1 = _mm256_loadu_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_loadu_ps(&vec_x[i+16ull]);
                vy2 = _mm256_loadu_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_storeu_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_loadu_ps(&vec_x[i+24ull]);
                vy3 = _mm256_loadu_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_storeu_ps(&vec_z[i+24ull],vz3);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif         
                vx4 = _mm256_loadu_ps(&vec_x[i+32ull]);
                vy4 = _mm256_loadu_ps(&vec_y[i+32ull]);
                vz4 = _mm256_add_ps(vx4,vy4);
                _mm256_storeu_ps(&vec_z[i+32ull],vz4);      
                vx5 = _mm256_loadu_ps(&vec_x[i+40ull]);
                vy5 = _mm256_loadu_ps(&vec_y[i+40ull]);
                vz5 = _mm256_add_ps(vx5,vy5);
                _mm256_storeu_ps(&vec_z[i+40ull],vz5);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+48ull],_MM_HINT_T0);
#endif              
                vx6 = _mm256_loadu_ps(&vec_x[i+48ull]);
                vy6 = _mm256_loadu_ps(&vec_y[i+48ull]);
                vz6 = _mm256_add_ps(vx6,vy6);
                _mm256_storeu_ps(&vec_z[i+48ull],vz6);
                vx7 = _mm256_loadu_ps(&vec_x[i+56ull]);
                vy7 = _mm256_loadu_ps(&vec_y[i+56ull]);
                vz7 = _mm256_add_ps(vx7,vy7);
                _mm256_storeu_ps(&vec_z[i+56ull],vz7);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+64ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+64ull],_MM_HINT_T0);
#endif      
                vx0 = _mm256_loadu_ps(&vec_x[i+64ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+64ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+64ull],vz0);  
                vx1 = _mm256_loadu_ps(&vec_x[i+72ull]);
                vy1 = _mm256_loadu_ps(&vec_y[i+72ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_storeu_ps(&vec_z[i+72ull],vz1);    
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+80ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+80ull],_MM_HINT_T0);
#endif                  
                vx2 = _mm256_loadu_ps(&vec_x[i+80ull]);
                vy2 = _mm256_loadu_ps(&vec_y[i+80ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_storeu_ps(&vec_z[i+80ull],vz2);          
                vx3 = _mm256_loadu_ps(&vec_x[i+88ull]);
                vy3 = _mm256_loadu_ps(&vec_y[i+88ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_storeu_ps(&vec_z[i+88ull],vz3);    
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+96ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+96ull],_MM_HINT_T0);
#endif         
                vx4 = _mm256_loadu_ps(&vec_x[i+96ull]);
                vy4 = _mm256_loadu_ps(&vec_y[i+96ull]);
                vz4 = _mm256_add_ps(vx4,vy4);
                _mm256_storeu_ps(&vec_z[i+96ull],vz4); 
                vx5 = _mm256_loadu_ps(&vec_x[i+104ull]);
                vy5 = _mm256_loadu_ps(&vec_y[i+104ull]);
                vz5 = _mm256_add_ps(vx5,vy5);
                _mm256_storeu_ps(&vec_z[i+104ull],vz5);    
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+112ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+112ull],_MM_HINT_T0);
#endif                           
                vx6 = _mm256_loadu_ps(&vec_x[i+112ull]);
                vy6 = _mm256_loadu_ps(&vec_y[i+112ull]);
                vz6 = _mm256_add_ps(vx6,vy6);
                _mm256_storeu_ps(&vec_z[i+112ull],vz6);         
                vx7 = _mm256_loadu_ps(&vec_x[i+120ull]);
                vy7 = _mm256_loadu_ps(&vec_y[i+120ull]);
                vz7 = _mm256_add_ps(vx7,vy7);
                _mm256_storeu_ps(&vec_z[i+120ull],vz7);                                 
            }

            for(i = 0ull; (i+95ull) < n; i += 96ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;
                __m256 vx4;
                __m256 vy4;
                __m256 vz4;
                __m256 vx5;
                __m256 vy5;
                __m256 vz5;
                __m256 vx6;
                __m256 vy6;
                __m256 vz6;
                __m256 vx7;
                __m256 vy7;
                __m256 vz7;

#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                vy1 = _mm256_loadu_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_loadu_ps(&vec_x[i+16ull]);
                vy2 = _mm256_loadu_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_storeu_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_loadu_ps(&vec_x[i+24ull]);
                vy3 = _mm256_loadu_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_storeu_ps(&vec_z[i+24ull],vz3);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif         
                vx4 = _mm256_loadu_ps(&vec_x[i+32ull]);
                vy4 = _mm256_loadu_ps(&vec_y[i+32ull]);
                vz4 = _mm256_add_ps(vx4,vy4);
                _mm256_storeu_ps(&vec_z[i+32ull],vz4);      
                vx5 = _mm256_loadu_ps(&vec_x[i+40ull]);
                vy5 = _mm256_loadu_ps(&vec_y[i+40ull]);
                vz5 = _mm256_add_ps(vx5,vy5);
                _mm256_storeu_ps(&vec_z[i+40ull],vz5);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+48ull],_MM_HINT_T0);
#endif              
                vx6 = _mm256_loadu_ps(&vec_x[i+48ull]);
                vy6 = _mm256_loadu_ps(&vec_y[i+48ull]);
                vz6 = _mm256_add_ps(vx6,vy6);
                _mm256_storeu_ps(&vec_z[i+48ull],vz6);
                vx7 = _mm256_loadu_ps(&vec_x[i+56ull]);
                vy7 = _mm256_loadu_ps(&vec_y[i+56ull]);
                vz7 = _mm256_add_ps(vx7,vy7);
                _mm256_storeu_ps(&vec_z[i+56ull],vz7);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+64ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+64ull],_MM_HINT_T0);
#endif      
                vx0 = _mm256_loadu_ps(&vec_x[i+64ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+64ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+64ull],vz0);  
                vx1 = _mm256_loadu_ps(&vec_x[i+72ull]);
                vy1 = _mm256_loadu_ps(&vec_y[i+72ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_storeu_ps(&vec_z[i+72ull],vz1);    
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+80ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+80ull],_MM_HINT_T0);
#endif                  
                vx2 = _mm256_loadu_ps(&vec_x[i+80ull]);
                vy2 = _mm256_loadu_ps(&vec_y[i+80ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_storeu_ps(&vec_z[i+80ull],vz2);          
                vx3 = _mm256_loadu_ps(&vec_x[i+88ull]);
                vy3 = _mm256_loadu_ps(&vec_y[i+88ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_storeu_ps(&vec_z[i+88ull],vz3);    
                                   
            }

            for(i = 0ull; (i+63ull) < n; i += 64ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;
                __m256 vx4;
                __m256 vy4;
                __m256 vz4;
                __m256 vx5;
                __m256 vy5;
                __m256 vz5;
                __m256 vx6;
                __m256 vy6;
                __m256 vz6;
                
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                vy1 = _mm256_loadu_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_loadu_ps(&vec_x[i+16ull]);
                vy2 = _mm256_loadu_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_storeu_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_loadu_ps(&vec_x[i+24ull]);
                vy3 = _mm256_loadu_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_storeu_ps(&vec_z[i+24ull],vz3);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif         
                vx4 = _mm256_loadu_ps(&vec_x[i+32ull]);
                vy4 = _mm256_loadu_ps(&vec_y[i+32ull]);
                vz4 = _mm256_add_ps(vx4,vy4);
                _mm256_storeu_ps(&vec_z[i+32ull],vz4);      
                vx5 = _mm256_loadu_ps(&vec_x[i+40ull]);
                vy5 = _mm256_loadu_ps(&vec_y[i+40ull]);
                vz5 = _mm256_add_ps(vx5,vy5);
                _mm256_storeu_ps(&vec_z[i+40ull],vz5);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+48ull],_MM_HINT_T0);
#endif              
                vx6 = _mm256_loadu_ps(&vec_x[i+48ull]);
                vy6 = _mm256_loadu_ps(&vec_y[i+48ull]);
                vz6 = _mm256_add_ps(vx6,vy6);
                _mm256_storeu_ps(&vec_z[i+48ull],vz6);
            }

            for(i = 0ull; (i+31ull) < n; i += 32ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;
                
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                vy1 = _mm256_loadu_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_loadu_ps(&vec_x[i+16ull]);
                vy2 = _mm256_loadu_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_storeu_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_loadu_ps(&vec_x[i+24ull]);
                vy3 = _mm256_loadu_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_storeu_ps(&vec_z[i+24ull],vz3);
               
            }

            for(i = 0ull; (i+15ull) < n; i += 16ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                vy1 = _mm256_loadu_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);        
            }
           
            for(i = 0ull; (i+7ull) < n; i += 8ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                               
                vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);                           
            }

            for(; (i+0ull) < n; i += 1ull)
            {
                const float x{vec_x[i]};
                const float y{vec_y[i]};
                const float z{x+y};
                vec_z[i] = z;
            }

            return (0);
}


std::int32_t 
gms::math::vadd_f32_2x256v_u16x_a(const float * __restrict__ vec_x,
                                  const float * __restrict__ vec_y,
                                  float       * __restrict__ vec_z,
                                  std::size_t          n)
{
            if(__builtin_expect(n==0ull,0)) { return -1;}


            std::size_t i;
             for(i = 0ull; (i+127ull) < n; i += 128ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;
                __m256 vx4;
                __m256 vy4;
                __m256 vz4;
                __m256 vx5;
                __m256 vy5;
                __m256 vz5;
                __m256 vx6;
                __m256 vy6;
                __m256 vz6;
                __m256 vx7;
                __m256 vy7;
                __m256 vz7;

#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                vy0 = _mm256_load_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                vy1 = _mm256_load_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_store_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_load_ps(&vec_x[i+16ull]);
                vy2 = _mm256_load_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_store_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_load_ps(&vec_x[i+24ull]);
                vy3 = _mm256_load_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_store_ps(&vec_z[i+24ull],vz3);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif         
                vx4 = _mm256_load_ps(&vec_x[i+32ull]);
                vy4 = _mm256_load_ps(&vec_y[i+32ull]);
                vz4 = _mm256_add_ps(vx4,vy4);
                _mm256_store_ps(&vec_z[i+32ull],vz4);      
                vx5 = _mm256_load_ps(&vec_x[i+40ull]);
                vy5 = _mm256_load_ps(&vec_y[i+40ull]);
                vz5 = _mm256_add_ps(vx5,vy5);
                _mm256_store_ps(&vec_z[i+40ull],vz5);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+48ull],_MM_HINT_T0);
#endif              
                vx6 = _mm256_load_ps(&vec_x[i+48ull]);
                vy6 = _mm256_load_ps(&vec_y[i+48ull]);
                vz6 = _mm256_add_ps(vx6,vy6);
                _mm256_store_ps(&vec_z[i+48ull],vz6);
                vx7 = _mm256_load_ps(&vec_x[i+56ull]);
                vy7 = _mm256_load_ps(&vec_y[i+56ull]);
                vz7 = _mm256_add_ps(vx7,vy7);
                _mm256_store_ps(&vec_z[i+56ull],vz7);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+64ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+64ull],_MM_HINT_T0);
#endif      
                vx0 = _mm256_load_ps(&vec_x[i+64ull]);
                vy0 = _mm256_load_ps(&vec_y[i+64ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+64ull],vz0);  
                vx1 = _mm256_load_ps(&vec_x[i+72ull]);
                vy1 = _mm256_load_ps(&vec_y[i+72ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_store_ps(&vec_z[i+72ull],vz1);    
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+80ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+80ull],_MM_HINT_T0);
#endif                  
                vx2 = _mm256_load_ps(&vec_x[i+80ull]);
                vy2 = _mm256_load_ps(&vec_y[i+80ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_store_ps(&vec_z[i+80ull],vz2);          
                vx3 = _mm256_load_ps(&vec_x[i+88ull]);
                vy3 = _mm256_load_ps(&vec_y[i+88ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_store_ps(&vec_z[i+88ull],vz3);    
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+96ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+96ull],_MM_HINT_T0);
#endif         
                vx4 = _mm256_load_ps(&vec_x[i+96ull]);
                vy4 = _mm256_load_ps(&vec_y[i+96ull]);
                vz4 = _mm256_add_ps(vx4,vy4);
                _mm256_store_ps(&vec_z[i+96ull],vz4); 
                vx5 = _mm256_load_ps(&vec_x[i+104ull]);
                vy5 = _mm256_load_ps(&vec_y[i+104ull]);
                vz5 = _mm256_add_ps(vx5,vy5);
                _mm256_store_ps(&vec_z[i+104ull],vz5);    
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+112ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+112ull],_MM_HINT_T0);
#endif                           
                vx6 = _mm256_load_ps(&vec_x[i+112ull]);
                vy6 = _mm256_load_ps(&vec_y[i+112ull]);
                vz6 = _mm256_add_ps(vx6,vy6);
                _mm256_store_ps(&vec_z[i+112ull],vz6);         
                vx7 = _mm256_load_ps(&vec_x[i+120ull]);
                vy7 = _mm256_load_ps(&vec_y[i+120ull]);
                vz7 = _mm256_add_ps(vx7,vy7);
                _mm256_store_ps(&vec_z[i+120ull],vz7);                                 
            }

            for(i = 0ull; (i+95ull) < n; i += 96ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;
                __m256 vx4;
                __m256 vy4;
                __m256 vz4;
                __m256 vx5;
                __m256 vy5;
                __m256 vz5;
                __m256 vx6;
                __m256 vy6;
                __m256 vz6;
                __m256 vx7;
                __m256 vy7;
                __m256 vz7;

#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                vy0 = _mm256_load_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                vy1 = _mm256_load_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_store_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_load_ps(&vec_x[i+16ull]);
                vy2 = _mm256_load_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_store_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_load_ps(&vec_x[i+24ull]);
                vy3 = _mm256_load_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_store_ps(&vec_z[i+24ull],vz3);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif         
                vx4 = _mm256_load_ps(&vec_x[i+32ull]);
                vy4 = _mm256_load_ps(&vec_y[i+32ull]);
                vz4 = _mm256_add_ps(vx4,vy4);
                _mm256_store_ps(&vec_z[i+32ull],vz4);      
                vx5 = _mm256_load_ps(&vec_x[i+40ull]);
                vy5 = _mm256_load_ps(&vec_y[i+40ull]);
                vz5 = _mm256_add_ps(vx5,vy5);
                _mm256_store_ps(&vec_z[i+40ull],vz5);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+48ull],_MM_HINT_T0);
#endif              
                vx6 = _mm256_load_ps(&vec_x[i+48ull]);
                vy6 = _mm256_load_ps(&vec_y[i+48ull]);
                vz6 = _mm256_add_ps(vx6,vy6);
                _mm256_store_ps(&vec_z[i+48ull],vz6);
                vx7 = _mm256_load_ps(&vec_x[i+56ull]);
                vy7 = _mm256_load_ps(&vec_y[i+56ull]);
                vz7 = _mm256_add_ps(vx7,vy7);
                _mm256_store_ps(&vec_z[i+56ull],vz7);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+64ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+64ull],_MM_HINT_T0);
#endif      
                vx0 = _mm256_load_ps(&vec_x[i+64ull]);
                vy0 = _mm256_load_ps(&vec_y[i+64ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+64ull],vz0);  
                vx1 = _mm256_load_ps(&vec_x[i+72ull]);
                vy1 = _mm256_load_ps(&vec_y[i+72ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_store_ps(&vec_z[i+72ull],vz1);    
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+80ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+80ull],_MM_HINT_T0);
#endif                  
                vx2 = _mm256_load_ps(&vec_x[i+80ull]);
                vy2 = _mm256_load_ps(&vec_y[i+80ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_store_ps(&vec_z[i+80ull],vz2);          
                vx3 = _mm256_load_ps(&vec_x[i+88ull]);
                vy3 = _mm256_load_ps(&vec_y[i+88ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_store_ps(&vec_z[i+88ull],vz3);    
                                   
            }

            for(i = 0ull; (i+63ull) < n; i += 64ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;
                __m256 vx4;
                __m256 vy4;
                __m256 vz4;
                __m256 vx5;
                __m256 vy5;
                __m256 vz5;
                __m256 vx6;
                __m256 vy6;
                __m256 vz6;
                
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                vy0 = _mm256_load_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                vy1 = _mm256_load_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_store_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_load_ps(&vec_x[i+16ull]);
                vy2 = _mm256_load_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_store_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_load_ps(&vec_x[i+24ull]);
                vy3 = _mm256_load_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_store_ps(&vec_z[i+24ull],vz3);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif         
                vx4 = _mm256_load_ps(&vec_x[i+32ull]);
                vy4 = _mm256_load_ps(&vec_y[i+32ull]);
                vz4 = _mm256_add_ps(vx4,vy4);
                _mm256_store_ps(&vec_z[i+32ull],vz4);      
                vx5 = _mm256_load_ps(&vec_x[i+40ull]);
                vy5 = _mm256_load_ps(&vec_y[i+40ull]);
                vz5 = _mm256_add_ps(vx5,vy5);
                _mm256_store_ps(&vec_z[i+40ull],vz5);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+48ull],_MM_HINT_T0);
#endif              
                vx6 = _mm256_load_ps(&vec_x[i+48ull]);
                vy6 = _mm256_load_ps(&vec_y[i+48ull]);
                vz6 = _mm256_add_ps(vx6,vy6);
                _mm256_store_ps(&vec_z[i+48ull],vz6);
            }

            for(i = 0ull; (i+31ull) < n; i += 32ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;
                
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                vy0 = _mm256_load_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                vy1 = _mm256_load_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_store_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_load_ps(&vec_x[i+16ull]);
                vy2 = _mm256_load_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_store_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_load_ps(&vec_x[i+24ull]);
                vy3 = _mm256_load_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_store_ps(&vec_z[i+24ull],vz3);
               
            }

            for(i = 0ull; (i+15ull) < n; i += 16ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                vy0 = _mm256_load_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                vy1 = _mm256_load_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_store_ps(&vec_z[i+8ull],vz1);        
            }
           
            for(i = 0ull; (i+7ull) < n; i += 8ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                               
                vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                vy0 = _mm256_load_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);                           
            }

            for(; (i+0ull) < n; i += 1ull)
            {
                const float x{vec_x[i]};
                const float y{vec_y[i]};
                const float z{x+y};
                vec_z[i] = z;
            }

            return (0);
}


std::int32_t 
gms::math::vadd_f32_2x256v_u10x_u(const float * __restrict__ vec_x,
                                  const float * __restrict__ vec_y,
                                  float       * __restrict__ vec_z,
                                  std::size_t          n)
{
            if(__builtin_expect(n==0ull,0)) { return -1;}

#if (ADD_VEC_YMM8R4_ADD_PEEL_LOOP) == 1
            while(((std::uintptr_t)&vec_z & 31) && n)
            {
                   const float x = *vec_x;
                   const float y = *vec_y;
                   const float z = x+y;
                   *vec_z        = z;
                   vec_x++;
                   vec_y++;
                   vec_z++;
                   n--;
            }
#endif 

            std::size_t i;
            for(i = 0ull; (i+79ull) < n; i += 80ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;
                __m256 vx4;
                __m256 vy4;
                __m256 vz4;
                __m256 vx5;
                __m256 vy5;
                __m256 vz5;
                __m256 vx6;
                __m256 vy6;
                __m256 vz6;
                __m256 vx7;
                __m256 vy7;
                __m256 vz7;

#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                vy1 = _mm256_loadu_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_loadu_ps(&vec_x[i+16ull]);
                vy2 = _mm256_loadu_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_storeu_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_loadu_ps(&vec_x[i+24ull]);
                vy3 = _mm256_loadu_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_storeu_ps(&vec_z[i+24ull],vz3);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif         
                vx4 = _mm256_loadu_ps(&vec_x[i+32ull]);
                vy4 = _mm256_loadu_ps(&vec_y[i+32ull]);
                vz4 = _mm256_add_ps(vx4,vy4);
                _mm256_storeu_ps(&vec_z[i+32ull],vz4);      
                vx5 = _mm256_loadu_ps(&vec_x[i+40ull]);
                vy5 = _mm256_loadu_ps(&vec_y[i+40ull]);
                vz5 = _mm256_add_ps(vx5,vy5);
                _mm256_storeu_ps(&vec_z[i+40ull],vz5);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+48ull],_MM_HINT_T0);
#endif              
                vx6 = _mm256_loadu_ps(&vec_x[i+48ull]);
                vy6 = _mm256_loadu_ps(&vec_y[i+48ull]);
                vz6 = _mm256_add_ps(vx6,vy6);
                _mm256_storeu_ps(&vec_z[i+48ull],vz6);
                vx7 = _mm256_loadu_ps(&vec_x[i+56ull]);
                vy7 = _mm256_loadu_ps(&vec_y[i+56ull]);
                vz7 = _mm256_add_ps(vx7,vy7);
                _mm256_storeu_ps(&vec_z[i+56ull],vz7);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+64ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+64ull],_MM_HINT_T0);
#endif      
                vx0 = _mm256_loadu_ps(&vec_x[i+64ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+64ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+64ull],vz0);  
                vx1 = _mm256_loadu_ps(&vec_x[i+72ull]);
                vy1 = _mm256_loadu_ps(&vec_y[i+72ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_storeu_ps(&vec_z[i+72ull],vz1);        
            }

            for(i = 0ull; (i+63ull) < n; i += 64ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;
                __m256 vx4;
                __m256 vy4;
                __m256 vz4;
                __m256 vx5;
                __m256 vy5;
                __m256 vz5;
                __m256 vx6;
                __m256 vy6;
                __m256 vz6;
                __m256 vx7;
                __m256 vy7;
                __m256 vz7;

#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                vy1 = _mm256_loadu_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_loadu_ps(&vec_x[i+16ull]);
                vy2 = _mm256_loadu_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_storeu_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_loadu_ps(&vec_x[i+24ull]);
                vy3 = _mm256_loadu_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_storeu_ps(&vec_z[i+24ull],vz3);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif         
                vx4 = _mm256_loadu_ps(&vec_x[i+32ull]);
                vy4 = _mm256_loadu_ps(&vec_y[i+32ull]);
                vz4 = _mm256_add_ps(vx4,vy4);
                _mm256_storeu_ps(&vec_z[i+32ull],vz4);      
                vx5 = _mm256_loadu_ps(&vec_x[i+40ull]);
                vy5 = _mm256_loadu_ps(&vec_y[i+40ull]);
                vz5 = _mm256_add_ps(vx5,vy5);
                _mm256_storeu_ps(&vec_z[i+40ull],vz5);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+48ull],_MM_HINT_T0);
#endif              
                vx6 = _mm256_loadu_ps(&vec_x[i+48ull]);
                vy6 = _mm256_loadu_ps(&vec_y[i+48ull]);
                vz6 = _mm256_add_ps(vx6,vy6);
                _mm256_storeu_ps(&vec_z[i+48ull],vz6);
                vx7 = _mm256_loadu_ps(&vec_x[i+56ull]);
                vy7 = _mm256_loadu_ps(&vec_y[i+56ull]);
                vz7 = _mm256_add_ps(vx7,vy7);
                _mm256_storeu_ps(&vec_z[i+56ull],vz7);    
                      
            }

            for(i = 0ull; (i+47ull) < n; i += 48ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;
                __m256 vx4;
                __m256 vy4;
                __m256 vz4;
                __m256 vx5;
                __m256 vy5;
                __m256 vz5;
    
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                vy1 = _mm256_loadu_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_loadu_ps(&vec_x[i+16ull]);
                vy2 = _mm256_loadu_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_storeu_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_loadu_ps(&vec_x[i+24ull]);
                vy3 = _mm256_loadu_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_storeu_ps(&vec_z[i+24ull],vz3);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif         
                vx4 = _mm256_loadu_ps(&vec_x[i+32ull]);
                vy4 = _mm256_loadu_ps(&vec_y[i+32ull]);
                vz4 = _mm256_add_ps(vx4,vy4);
                _mm256_storeu_ps(&vec_z[i+32ull],vz4);      
                vx5 = _mm256_loadu_ps(&vec_x[i+40ull]);
                vy5 = _mm256_loadu_ps(&vec_y[i+40ull]);
                vz5 = _mm256_add_ps(vx5,vy5);
                _mm256_storeu_ps(&vec_z[i+40ull],vz5);                
                           
            }

            for(i = 0ull; (i+31ull) < n; i += 32ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;
                   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                vy1 = _mm256_loadu_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_loadu_ps(&vec_x[i+16ull]);
                vy2 = _mm256_loadu_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_storeu_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_loadu_ps(&vec_x[i+24ull]);
                vy3 = _mm256_loadu_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_storeu_ps(&vec_z[i+24ull],vz3);      
                                        
            }

            for(i = 0ull; (i+15ull) < n; i += 16ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                         
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                vy1 = _mm256_loadu_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);        
                                             
            }

            for(i = 0ull; (i+7ull) < n; i += 8ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;                                         
                vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);   
                                                  
            }

            for(; (i+0ull) < n; i += 1ull)
            {
                const float x{vec_x[i]};
                const float y{vec_y[i]};
                const float z{x+y};
                vec_z[i] = z;
            }

            return (0);
}


std::int32_t 
gms::math::vadd_f32_2x256v_u10x_a(const float * __restrict__ vec_x,
                                  const float * __restrict__ vec_y,
                                  float       * __restrict__ vec_z,
                                  std::size_t          n)
{
            if(__builtin_expect(n==0ull,0)) { return -1;}


            std::size_t i;
            for(i = 0ull; (i+79ull) < n; i += 80ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;
                __m256 vx4;
                __m256 vy4;
                __m256 vz4;
                __m256 vx5;
                __m256 vy5;
                __m256 vz5;
                __m256 vx6;
                __m256 vy6;
                __m256 vz6;
                __m256 vx7;
                __m256 vy7;
                __m256 vz7;

#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                vy0 = _mm256_load_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                vy1 = _mm256_load_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_store_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_load_ps(&vec_x[i+16ull]);
                vy2 = _mm256_load_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_store_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_load_ps(&vec_x[i+24ull]);
                vy3 = _mm256_load_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_store_ps(&vec_z[i+24ull],vz3);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif         
                vx4 = _mm256_load_ps(&vec_x[i+32ull]);
                vy4 = _mm256_load_ps(&vec_y[i+32ull]);
                vz4 = _mm256_add_ps(vx4,vy4);
                _mm256_store_ps(&vec_z[i+32ull],vz4);      
                vx5 = _mm256_load_ps(&vec_x[i+40ull]);
                vy5 = _mm256_load_ps(&vec_y[i+40ull]);
                vz5 = _mm256_add_ps(vx5,vy5);
                _mm256_store_ps(&vec_z[i+40ull],vz5);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+48ull],_MM_HINT_T0);
#endif              
                vx6 = _mm256_load_ps(&vec_x[i+48ull]);
                vy6 = _mm256_load_ps(&vec_y[i+48ull]);
                vz6 = _mm256_add_ps(vx6,vy6);
                _mm256_store_ps(&vec_z[i+48ull],vz6);
                vx7 = _mm256_load_ps(&vec_x[i+56ull]);
                vy7 = _mm256_load_ps(&vec_y[i+56ull]);
                vz7 = _mm256_add_ps(vx7,vy7);
                _mm256_store_ps(&vec_z[i+56ull],vz7);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+64ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+64ull],_MM_HINT_T0);
#endif      
                vx0 = _mm256_load_ps(&vec_x[i+64ull]);
                vy0 = _mm256_load_ps(&vec_y[i+64ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+64ull],vz0);  
                vx1 = _mm256_load_ps(&vec_x[i+72ull]);
                vy1 = _mm256_load_ps(&vec_y[i+72ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_store_ps(&vec_z[i+72ull],vz1);        
            }

            for(i = 0ull; (i+63ull) < n; i += 64ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;
                __m256 vx4;
                __m256 vy4;
                __m256 vz4;
                __m256 vx5;
                __m256 vy5;
                __m256 vz5;
                __m256 vx6;
                __m256 vy6;
                __m256 vz6;
                __m256 vx7;
                __m256 vy7;
                __m256 vz7;

#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                vy0 = _mm256_load_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                vy1 = _mm256_load_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_store_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_load_ps(&vec_x[i+16ull]);
                vy2 = _mm256_load_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_store_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_load_ps(&vec_x[i+24ull]);
                vy3 = _mm256_load_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_store_ps(&vec_z[i+24ull],vz3);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif         
                vx4 = _mm256_load_ps(&vec_x[i+32ull]);
                vy4 = _mm256_load_ps(&vec_y[i+32ull]);
                vz4 = _mm256_add_ps(vx4,vy4);
                _mm256_store_ps(&vec_z[i+32ull],vz4);      
                vx5 = _mm256_load_ps(&vec_x[i+40ull]);
                vy5 = _mm256_load_ps(&vec_y[i+40ull]);
                vz5 = _mm256_add_ps(vx5,vy5);
                _mm256_store_ps(&vec_z[i+40ull],vz5);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+48ull],_MM_HINT_T0);
#endif              
                vx6 = _mm256_load_ps(&vec_x[i+48ull]);
                vy6 = _mm256_load_ps(&vec_y[i+48ull]);
                vz6 = _mm256_add_ps(vx6,vy6);
                _mm256_store_ps(&vec_z[i+48ull],vz6);
                vx7 = _mm256_load_ps(&vec_x[i+56ull]);
                vy7 = _mm256_load_ps(&vec_y[i+56ull]);
                vz7 = _mm256_add_ps(vx7,vy7);
                _mm256_store_ps(&vec_z[i+56ull],vz7);    
                      
            }

            for(i = 0ull; (i+47ull) < n; i += 48ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;
                __m256 vx4;
                __m256 vy4;
                __m256 vz4;
                __m256 vx5;
                __m256 vy5;
                __m256 vz5;
    
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                vy0 = _mm256_load_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                vy1 = _mm256_load_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_store_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_load_ps(&vec_x[i+16ull]);
                vy2 = _mm256_load_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_store_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_load_ps(&vec_x[i+24ull]);
                vy3 = _mm256_load_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_store_ps(&vec_z[i+24ull],vz3);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif         
                vx4 = _mm256_load_ps(&vec_x[i+32ull]);
                vy4 = _mm256_load_ps(&vec_y[i+32ull]);
                vz4 = _mm256_add_ps(vx4,vy4);
                _mm256_store_ps(&vec_z[i+32ull],vz4);      
                vx5 = _mm256_load_ps(&vec_x[i+40ull]);
                vy5 = _mm256_load_ps(&vec_y[i+40ull]);
                vz5 = _mm256_add_ps(vx5,vy5);
                _mm256_store_ps(&vec_z[i+40ull],vz5);                
                           
            }

            for(i = 0ull; (i+31ull) < n; i += 32ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;
                   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                vy0 = _mm256_load_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                vy1 = _mm256_load_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_store_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_load_ps(&vec_x[i+16ull]);
                vy2 = _mm256_load_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_store_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_load_ps(&vec_x[i+24ull]);
                vy3 = _mm256_load_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_store_ps(&vec_z[i+24ull],vz3);      
                                        
            }

            for(i = 0ull; (i+15ull) < n; i += 16ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                         
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                vy0 = _mm256_load_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                vy1 = _mm256_load_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_store_ps(&vec_z[i+8ull],vz1);        
                                             
            }

            for(i = 0ull; (i+7ull) < n; i += 8ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;                                         
                vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                vy0 = _mm256_load_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);   
                                                  
            }

            for(; (i+0ull) < n; i += 1ull)
            {
                const float x{vec_x[i]};
                const float y{vec_y[i]};
                const float z{x+y};
                vec_z[i] = z;
            }

            return (0);
}


std::int32_t 
gms::math::vadd_f32_2x256v_u6x_u( const float * __restrict__ vec_x,
                                  const float * __restrict__ vec_y,
                                  float       * __restrict__ vec_z,
                                  std::size_t          n)
{
            if(__builtin_expect(n==0ull,0)) { return -1;}

#if (ADD_VEC_YMM8R4_ADD_PEEL_LOOP) == 1
            while(((std::uintptr_t)&vec_z & 31) && n)
            {
                   const float x = *vec_x;
                   const float y = *vec_y;
                   const float z = x+y;
                   *vec_z        = z;
                   vec_x++;
                   vec_y++;
                   vec_z++;
                   n--;
            }
#endif 

            std::size_t i;
            for(i = 0ull; (i+47ull) < n; i += 48ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;
                __m256 vx4;
                __m256 vy4;
                __m256 vz4;
                __m256 vx5;
                __m256 vy5;
                __m256 vz5;
                
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                vy1 = _mm256_loadu_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_loadu_ps(&vec_x[i+16ull]);
                vy2 = _mm256_loadu_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_storeu_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_loadu_ps(&vec_x[i+24ull]);
                vy3 = _mm256_loadu_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_storeu_ps(&vec_z[i+24ull],vz3);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif         
                vx4 = _mm256_loadu_ps(&vec_x[i+32ull]);
                vy4 = _mm256_loadu_ps(&vec_y[i+32ull]);
                vz4 = _mm256_add_ps(vx4,vy4);
                _mm256_storeu_ps(&vec_z[i+32ull],vz4);      
                vx5 = _mm256_loadu_ps(&vec_x[i+40ull]);
                vy5 = _mm256_loadu_ps(&vec_y[i+40ull]);
                vz5 = _mm256_add_ps(vx5,vy5);
                _mm256_storeu_ps(&vec_z[i+40ull],vz5); 
                    
            }

           for(i = 0ull; (i+31ull) < n; i += 32ull)
           {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;

#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                vy1 = _mm256_loadu_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_loadu_ps(&vec_x[i+16ull]);
                vy2 = _mm256_loadu_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_storeu_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_loadu_ps(&vec_x[i+24ull]);
                vy3 = _mm256_loadu_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_storeu_ps(&vec_z[i+24ull],vz3);                
            
            }

            for(i = 0ull; (i+15ull) < n; i += 16ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
               
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                vy1 = _mm256_loadu_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);   
                                  
            }

            for(i = 0ull; (i+7ull) < n; i += 8ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                              
                vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                vy0 = _mm256_loadu_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);                                                
            }

           
            for(; (i+0ull) < n; i += 1ull)
            {
                const float x{vec_x[i]};
                const float y{vec_y[i]};
                const float z{x+y};
                vec_z[i] = z;
            }

            return (0);
}


std::int32_t 
gms::math::vadd_f32_2x256v_u6x_a( const float * __restrict__ vec_x,
                                  const float * __restrict__ vec_y,
                                  float       * __restrict__ vec_z,
                                  std::size_t          n)
{
            if(__builtin_expect(n==0ull,0)) { return -1;}


            std::size_t i;
             for(i = 0ull; (i+47ull) < n; i += 48ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;
                __m256 vx4;
                __m256 vy4;
                __m256 vz4;
                __m256 vx5;
                __m256 vy5;
                __m256 vz5;
                
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                vy0 = _mm256_load_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                vy1 = _mm256_load_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_store_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_load_ps(&vec_x[i+16ull]);
                vy2 = _mm256_load_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_store_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_load_ps(&vec_x[i+24ull]);
                vy3 = _mm256_load_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_store_ps(&vec_z[i+24ull],vz3);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+32ull],_MM_HINT_T0);
#endif         
                vx4 = _mm256_load_ps(&vec_x[i+32ull]);
                vy4 = _mm256_load_ps(&vec_y[i+32ull]);
                vz4 = _mm256_add_ps(vx4,vy4);
                _mm256_store_ps(&vec_z[i+32ull],vz4);      
                vx5 = _mm256_load_ps(&vec_x[i+40ull]);
                vy5 = _mm256_load_ps(&vec_y[i+40ull]);
                vz5 = _mm256_add_ps(vx5,vy5);
                _mm256_store_ps(&vec_z[i+40ull],vz5); 
                    
            }

           for(i = 0ull; (i+31ull) < n; i += 32ull)
           {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
                __m256 vx2;
                __m256 vy2;
                __m256 vz2;
                __m256 vx3;
                __m256 vy3;
                __m256 vz3;

#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                vy0 = _mm256_load_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                vy1 = _mm256_load_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_store_ps(&vec_z[i+8ull],vz1);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+16ull],_MM_HINT_T0);
#endif          
                vx2 = _mm256_load_ps(&vec_x[i+16ull]);
                vy2 = _mm256_load_ps(&vec_y[i+16ull]);
                vz2 = _mm256_add_ps(vx2,vy2);
                _mm256_store_ps(&vec_z[i+16ull],vz2);   
                vx3 = _mm256_load_ps(&vec_x[i+24ull]);
                vy3 = _mm256_load_ps(&vec_y[i+24ull]);
                vz3 = _mm256_add_ps(vx3,vy3);
                _mm256_store_ps(&vec_z[i+24ull],vz3);                
            
            }

            for(i = 0ull; (i+15ull) < n; i += 16ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                __m256 vx1;
                __m256 vy1;
                __m256 vz1;
               
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
                            _mm_prefetch((const char *)&vec_y[i+0ull],_MM_HINT_T0);
#endif         
                vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                vy0 = _mm256_load_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);   
                vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                vy1 = _mm256_load_ps(&vec_y[i+8ull]);
                vz1 = _mm256_add_ps(vx1,vy1);
                _mm256_store_ps(&vec_z[i+8ull],vz1);   
                                  
            }

            for(i = 0ull; (i+7ull) < n; i += 8ull)
            {
/*
                VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
*/
                __m256 vx0;
                __m256 vy0;
                __m256 vz0;
                              
                vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                vy0 = _mm256_load_ps(&vec_y[i+0ull]);
                vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);                                                
            }
           
            for(; (i+0ull) < n; i += 1ull)
            {
                const float x{vec_x[i]};
                const float y{vec_y[i]};
                const float z{x+y};
                vec_z[i] = z;
            }

            return (0);
}


std::int32_t 
gms::math::vadd_f32_1x256v_1xf32_u16x_u(const float * __restrict__ vec_x,
                                        const float                y,
                                        float       * __restrict__ vec_z,
                                        std::size_t          n) 
{
           if(__builtin_expect(n==0ull,0)) { return -1;}

#if (ADD_VEC_YMM8R4_ADD_PEEL_LOOP) == 1
            while(((std::uintptr_t)&vec_z & 31) && n)
            {
                   const float x = *vec_x;
                   const float z = x+y;
                  *vec_z         = z;
                   vec_x++;
                   vec_z++;
                   n--;
            }
#endif   

           const __m256 vy0{_mm256_set1_ps(y)};
           std::size_t i; 
           for(i = 0ull; (i+127ull) < n; i += 128ull) 
           {
                /*
                   VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
                */
                __m256 vx0;
                __m256 vz0;
                __m256 vx1;
                __m256 vz1;
                __m256 vx2;
                __m256 vz2;
                __m256 vx3;
                __m256 vz3;
                __m256 vx4;
                __m256 vz4;
                __m256 vx5;
                __m256 vz5;
                __m256 vx6;
                __m256 vz6;
                __m256 vx7;
                __m256 vz7;

#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
#endif                  
                 vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                 vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);                  
                 vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                 vz1 = _mm256_add_ps(vx1,vy0);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
#endif                               
                 vx2 = _mm256_loadu_ps(&vec_x[i+16ull]);
                 vz2 = _mm256_add_ps(vx2,vy0);
                _mm256_storeu_ps(&vec_z[i+16ull],vz2);                 
                 vx3 = _mm256_loadu_ps(&vec_x[i+24ull]);
                 vz3 = _mm256_add_ps(vx3,vy0);
                _mm256_storeu_ps(&vec_z[i+24ull],vz3);    
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
#endif             
                 vx4 = _mm256_loadu_ps(&vec_x[i+32ull]);
                 vz4 = _mm256_add_ps(vx4,vy0);
                _mm256_storeu_ps(&vec_z[i+32ull],vz4);                 
                 vx5 = _mm256_loadu_ps(&vec_x[i+40ull]);
                 vz5 = _mm256_add_ps(vx5,vy0);
                _mm256_storeu_ps(&vec_z[i+40ull],vz5);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
#endif
                 vx6 = _mm256_loadu_ps(&vec_x[i+48ull]);
                 vz6 = _mm256_add_ps(vx6,vy0);
                _mm256_storeu_ps(&vec_z[i+48ull],vz6);                  
                 vx7 = _mm256_loadu_ps(&vec_x[i+56ull]);
                 vz7 = _mm256_add_ps(vx7,vy0);
                _mm256_storeu_ps(&vec_z[i+56ull],vz7);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+64ull],_MM_HINT_T0);
#endif                              
                 vx0 = _mm256_loadu_ps(&vec_x[i+64ull]);
                 vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+64ull],vz0);                 
                 vx1 = _mm256_loadu_ps(&vec_x[i+72ull]);
                 vz1 = _mm256_add_ps(vx1,vy0);
                _mm256_storeu_ps(&vec_z[i+72ull],vz1);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+80ull],_MM_HINT_T0);
#endif             
                 vx2 = _mm256_loadu_ps(&vec_x[i+80ull]);
                 vz2 = _mm256_add_ps(vx2,vy0);
                _mm256_storeu_ps(&vec_z[i+80ull],vz2);                  
                 vx3 = _mm256_loadu_ps(&vec_x[i+88ull]);
                 vz3 = _mm256_add_ps(vx3,vy0);
                _mm256_storeu_ps(&vec_z[i+88ull],vz3);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+96ull],_MM_HINT_T0);
#endif                  
                 vx4 = _mm256_loadu_ps(&vec_x[i+96ull]);
                 vz4 = _mm256_add_ps(vx4,vy0);
                _mm256_storeu_ps(&vec_z[i+96ull],vz4);                 
                 vx5 = _mm256_loadu_ps(&vec_x[i+104ull]);
                 vz5 = _mm256_add_ps(vx5,vy0);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+104ull],_MM_HINT_T0);
#endif
                _mm256_storeu_ps(&vec_z[i+104ull],vz5);                
                 vx6 = _mm256_loadu_ps(&vec_x[i+112ull]);
                 vz6 = _mm256_add_ps(vx6,vy0);
                _mm256_storeu_ps(&vec_z[i+112ull],vz6);               
                 vx7 = _mm256_loadu_ps(&vec_x[i+120ull]);
                 vz7 = _mm256_add_ps(vx7,vy0);
                _mm256_storeu_ps(&vec_z[i+120ull],vz7);
           }

           for(i = 0ull; (i+95ull) < n; i += 96ull) 
           {
                /*
                   VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
                */
                __m256 vx0;
                __m256 vz0;
                __m256 vx1;
                __m256 vz1;
                __m256 vx2;
                __m256 vz2;
                __m256 vx3;
                __m256 vz3;
                __m256 vx4;
                __m256 vz4;
                __m256 vx5;
                __m256 vz5;
                __m256 vx6;
                __m256 vz6;
                __m256 vx7;
                __m256 vz7;

#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
#endif                  
                 vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                 vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);                  
                 vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                 vz1 = _mm256_add_ps(vx1,vy0);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
#endif                               
                 vx2 = _mm256_loadu_ps(&vec_x[i+16ull]);
                 vz2 = _mm256_add_ps(vx2,vy0);
                _mm256_storeu_ps(&vec_z[i+16ull],vz2);                 
                 vx3 = _mm256_loadu_ps(&vec_x[i+24ull]);
                 vz3 = _mm256_add_ps(vx3,vy0);
                _mm256_storeu_ps(&vec_z[i+24ull],vz3);    
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
#endif             
                 vx4 = _mm256_loadu_ps(&vec_x[i+32ull]);
                 vz4 = _mm256_add_ps(vx4,vy0);
                _mm256_storeu_ps(&vec_z[i+32ull],vz4);                 
                 vx5 = _mm256_loadu_ps(&vec_x[i+40ull]);
                 vz5 = _mm256_add_ps(vx5,vy0);
                _mm256_storeu_ps(&vec_z[i+40ull],vz5);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
#endif
                 vx6 = _mm256_loadu_ps(&vec_x[i+48ull]);
                 vz6 = _mm256_add_ps(vx6,vy0);
                _mm256_storeu_ps(&vec_z[i+48ull],vz6);                  
                 vx7 = _mm256_loadu_ps(&vec_x[i+56ull]);
                 vz7 = _mm256_add_ps(vx7,vy0);
                _mm256_storeu_ps(&vec_z[i+56ull],vz7);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+64ull],_MM_HINT_T0);
#endif                              
                 vx0 = _mm256_loadu_ps(&vec_x[i+64ull]);
                 vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+64ull],vz0);                 
                 vx1 = _mm256_loadu_ps(&vec_x[i+72ull]);
                 vz1 = _mm256_add_ps(vx1,vy0);
                _mm256_storeu_ps(&vec_z[i+72ull],vz1);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+80ull],_MM_HINT_T0);
#endif             
                 vx2 = _mm256_loadu_ps(&vec_x[i+80ull]);
                 vz2 = _mm256_add_ps(vx2,vy0);
                _mm256_storeu_ps(&vec_z[i+80ull],vz2);     
           }

           for(i = 0ull; (i+63ull) < n; i += 64ull) 
           {
                /*
                   VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
                */
                __m256 vx0;
                __m256 vz0;
                __m256 vx1;
                __m256 vz1;
                __m256 vx2;
                __m256 vz2;
                __m256 vx3;
                __m256 vz3;
                __m256 vx4;
                __m256 vz4;
                __m256 vx5;
                __m256 vz5;
                __m256 vx6;
                __m256 vz6;
                __m256 vx7;
                __m256 vz7;

#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
#endif                  
                 vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                 vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);                  
                 vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                 vz1 = _mm256_add_ps(vx1,vy0);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
#endif                               
                 vx2 = _mm256_loadu_ps(&vec_x[i+16ull]);
                 vz2 = _mm256_add_ps(vx2,vy0);
                _mm256_storeu_ps(&vec_z[i+16ull],vz2);                 
                 vx3 = _mm256_loadu_ps(&vec_x[i+24ull]);
                 vz3 = _mm256_add_ps(vx3,vy0);
                _mm256_storeu_ps(&vec_z[i+24ull],vz3);    
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
#endif             
                 vx4 = _mm256_loadu_ps(&vec_x[i+32ull]);
                 vz4 = _mm256_add_ps(vx4,vy0);
                _mm256_storeu_ps(&vec_z[i+32ull],vz4);                 
                 vx5 = _mm256_loadu_ps(&vec_x[i+40ull]);
                 vz5 = _mm256_add_ps(vx5,vy0);
                _mm256_storeu_ps(&vec_z[i+40ull],vz5);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
#endif
                 vx6 = _mm256_loadu_ps(&vec_x[i+48ull]);
                 vz6 = _mm256_add_ps(vx6,vy0);
                _mm256_storeu_ps(&vec_z[i+48ull],vz6);                  
                 vx7 = _mm256_loadu_ps(&vec_x[i+56ull]);
                 vz7 = _mm256_add_ps(vx7,vy0);
                _mm256_storeu_ps(&vec_z[i+56ull],vz7);  
           }

           for(i = 0ull; (i+47ull) < n; i += 48ull) 
           {
                /*
                   VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
                */
                __m256 vx0;
                __m256 vz0;
                __m256 vx1;
                __m256 vz1;
                __m256 vx2;
                __m256 vz2;
                __m256 vx3;
                __m256 vz3;
                __m256 vx4;
                __m256 vz4;
                __m256 vx5;
                __m256 vz5;
               
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
#endif                  
                 vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                 vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);                  
                 vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                 vz1 = _mm256_add_ps(vx1,vy0);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
#endif                               
                 vx2 = _mm256_loadu_ps(&vec_x[i+16ull]);
                 vz2 = _mm256_add_ps(vx2,vy0);
                _mm256_storeu_ps(&vec_z[i+16ull],vz2);                 
                 vx3 = _mm256_loadu_ps(&vec_x[i+24ull]);
                 vz3 = _mm256_add_ps(vx3,vy0);
                _mm256_storeu_ps(&vec_z[i+24ull],vz3);    
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
#endif             
                 vx4 = _mm256_loadu_ps(&vec_x[i+32ull]);
                 vz4 = _mm256_add_ps(vx4,vy0);
                _mm256_storeu_ps(&vec_z[i+32ull],vz4);                 
                 vx5 = _mm256_loadu_ps(&vec_x[i+40ull]);
                 vz5 = _mm256_add_ps(vx5,vy0);
                _mm256_storeu_ps(&vec_z[i+40ull],vz5);              
           }

           for(i = 0ull; (i+31ull) < n; i += 32ull) 
           {
                /*
                   VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
                */
                __m256 vx0;
                __m256 vz0;
                __m256 vx1;
                __m256 vz1;
                __m256 vx2;
                __m256 vz2;
                __m256 vx3;
                __m256 vz3;
               
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
#endif                  
                 vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                 vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);                  
                 vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                 vz1 = _mm256_add_ps(vx1,vy0);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
#endif                               
                 vx2 = _mm256_loadu_ps(&vec_x[i+16ull]);
                 vz2 = _mm256_add_ps(vx2,vy0);
                _mm256_storeu_ps(&vec_z[i+16ull],vz2);                 
                 vx3 = _mm256_loadu_ps(&vec_x[i+24ull]);
                 vz3 = _mm256_add_ps(vx3,vy0);
                _mm256_storeu_ps(&vec_z[i+24ull],vz3);                                 
           }

           for(i = 0ull; (i+15ull) < n; i += 16ull) 
           {
                /*
                   VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
                */
                __m256 vx0;
                __m256 vz0;
                __m256 vx1;
                __m256 vz1;
                               
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
#endif                  
                 vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                 vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);                  
                 vx1 = _mm256_loadu_ps(&vec_x[i+8ull]);
                 vz1 = _mm256_add_ps(vx1,vy0);
                _mm256_storeu_ps(&vec_z[i+8ull],vz1);                                                     
           }

           for(i = 0ull; (i+7ull) < n; i += 8ull) 
           {
                /*
                   VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
                */
                __m256 vx0;
                __m256 vz0;
                               
                 vx0 = _mm256_loadu_ps(&vec_x[i+0ull]);
                 vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_storeu_ps(&vec_z[i+0ull],vz0);                                                                       
           }

           for(; (i+0ull) < n; i += 1ull) 
           {
                const float x{vec_x[i]};
                const float z{x+y};
                vec_z[i] =  z;
           }

           return (0);
}


std::int32_t 
gms::math::vadd_f32_1x256v_1xf32_u16x_a(const float * __restrict__ vec_x,
                                        const float                y,
                                        float       * __restrict__ vec_z,
                                        std::size_t          n) 
{
           if(__builtin_expect(n==0ull,0)) { return -1;}

           const __m256 vy0{_mm256_set1_ps(y)};
           std::size_t i; 
           
           for(i = 0ull; (i+127ull) < n; i += 128ull) 
           {
                /*
                   VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
                */
                __m256 vx0;
                __m256 vz0;
                __m256 vx1;
                __m256 vz1;
                __m256 vx2;
                __m256 vz2;
                __m256 vx3;
                __m256 vz3;
                __m256 vx4;
                __m256 vz4;
                __m256 vx5;
                __m256 vz5;
                __m256 vx6;
                __m256 vz6;
                __m256 vx7;
                __m256 vz7;

#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
#endif                  
                 vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                 vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);                  
                 vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                 vz1 = _mm256_add_ps(vx1,vy0);
                _mm256_store_ps(&vec_z[i+8ull],vz1);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
#endif                               
                 vx2 = _mm256_load_ps(&vec_x[i+16ull]);
                 vz2 = _mm256_add_ps(vx2,vy0);
                _mm256_store_ps(&vec_z[i+16ull],vz2);                 
                 vx3 = _mm256_load_ps(&vec_x[i+24ull]);
                 vz3 = _mm256_add_ps(vx3,vy0);
                _mm256_store_ps(&vec_z[i+24ull],vz3);    
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
#endif             
                 vx4 = _mm256_load_ps(&vec_x[i+32ull]);
                 vz4 = _mm256_add_ps(vx4,vy0);
                _mm256_store_ps(&vec_z[i+32ull],vz4);                 
                 vx5 = _mm256_load_ps(&vec_x[i+40ull]);
                 vz5 = _mm256_add_ps(vx5,vy0);
                _mm256_store_ps(&vec_z[i+40ull],vz5);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
#endif
                 vx6 = _mm256_load_ps(&vec_x[i+48ull]);
                 vz6 = _mm256_add_ps(vx6,vy0);
                _mm256_store_ps(&vec_z[i+48ull],vz6);                  
                 vx7 = _mm256_load_ps(&vec_x[i+56ull]);
                 vz7 = _mm256_add_ps(vx7,vy0);
                _mm256_store_ps(&vec_z[i+56ull],vz7);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+64ull],_MM_HINT_T0);
#endif                              
                 vx0 = _mm256_load_ps(&vec_x[i+64ull]);
                 vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+64ull],vz0);                 
                 vx1 = _mm256_load_ps(&vec_x[i+72ull]);
                 vz1 = _mm256_add_ps(vx1,vy0);
                _mm256_store_ps(&vec_z[i+72ull],vz1);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+80ull],_MM_HINT_T0);
#endif             
                 vx2 = _mm256_load_ps(&vec_x[i+80ull]);
                 vz2 = _mm256_add_ps(vx2,vy0);
                _mm256_store_ps(&vec_z[i+80ull],vz2);                  
                 vx3 = _mm256_load_ps(&vec_x[i+88ull]);
                 vz3 = _mm256_add_ps(vx3,vy0);
                _mm256_store_ps(&vec_z[i+88ull],vz3);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+96ull],_MM_HINT_T0);
#endif                  
                 vx4 = _mm256_load_ps(&vec_x[i+96ull]);
                 vz4 = _mm256_add_ps(vx4,vy0);
                _mm256_store_ps(&vec_z[i+96ull],vz4);                 
                 vx5 = _mm256_load_ps(&vec_x[i+104ull]);
                 vz5 = _mm256_add_ps(vx5,vy0);
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+104ull],_MM_HINT_T0);
#endif
                _mm256_store_ps(&vec_z[i+104ull],vz5);                
                 vx6 = _mm256_load_ps(&vec_x[i+112ull]);
                 vz6 = _mm256_add_ps(vx6,vy0);
                _mm256_store_ps(&vec_z[i+112ull],vz6);               
                 vx7 = _mm256_loadu_ps(&vec_x[i+120ull]);
                 vz7 = _mm256_add_ps(vx7,vy0);
                _mm256_store_ps(&vec_z[i+120ull],vz7);
           }

           for(i = 0ull; (i+95ull) < n; i += 96ull) 
           {
                /*
                   VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
                */
                __m256 vx0;
                __m256 vz0;
                __m256 vx1;
                __m256 vz1;
                __m256 vx2;
                __m256 vz2;
                __m256 vx3;
                __m256 vz3;
                __m256 vx4;
                __m256 vz4;
                __m256 vx5;
                __m256 vz5;
                __m256 vx6;
                __m256 vz6;
                __m256 vx7;
                __m256 vz7;

#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
#endif                  
                 vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                 vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);                  
                 vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                 vz1 = _mm256_add_ps(vx1,vy0);
                _mm256_store_ps(&vec_z[i+8ull],vz1);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
#endif                               
                 vx2 = _mm256_load_ps(&vec_x[i+16ull]);
                 vz2 = _mm256_add_ps(vx2,vy0);
                _mm256_store_ps(&vec_z[i+16ull],vz2);                 
                 vx3 = _mm256_load_ps(&vec_x[i+24ull]);
                 vz3 = _mm256_add_ps(vx3,vy0);
                _mm256_store_ps(&vec_z[i+24ull],vz3);    
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
#endif             
                 vx4 = _mm256_load_ps(&vec_x[i+32ull]);
                 vz4 = _mm256_add_ps(vx4,vy0);
                _mm256_store_ps(&vec_z[i+32ull],vz4);                 
                 vx5 = _mm256_load_ps(&vec_x[i+40ull]);
                 vz5 = _mm256_add_ps(vx5,vy0);
                _mm256_store_ps(&vec_z[i+40ull],vz5);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
#endif
                 vx6 = _mm256_load_ps(&vec_x[i+48ull]);
                 vz6 = _mm256_add_ps(vx6,vy0);
                _mm256_store_ps(&vec_z[i+48ull],vz6);                  
                 vx7 = _mm256_load_ps(&vec_x[i+56ull]);
                 vz7 = _mm256_add_ps(vx7,vy0);
                _mm256_store_ps(&vec_z[i+56ull],vz7);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+64ull],_MM_HINT_T0);
#endif                              
                 vx0 = _mm256_load_ps(&vec_x[i+64ull]);
                 vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+64ull],vz0);                 
                 vx1 = _mm256_load_ps(&vec_x[i+72ull]);
                 vz1 = _mm256_add_ps(vx1,vy0);
                _mm256_store_ps(&vec_z[i+72ull],vz1);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+80ull],_MM_HINT_T0);
#endif             
                 vx2 = _mm256_load_ps(&vec_x[i+80ull]);
                 vz2 = _mm256_add_ps(vx2,vy0);
                _mm256_store_ps(&vec_z[i+80ull],vz2);     
           }

           for(i = 0ull; (i+63ull) < n; i += 64ull) 
           {
                /*
                   VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
                */
                __m256 vx0;
                __m256 vz0;
                __m256 vx1;
                __m256 vz1;
                __m256 vx2;
                __m256 vz2;
                __m256 vx3;
                __m256 vz3;
                __m256 vx4;
                __m256 vz4;
                __m256 vx5;
                __m256 vz5;
                __m256 vx6;
                __m256 vz6;
                __m256 vx7;
                __m256 vz7;

#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
#endif                  
                 vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                 vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);                  
                 vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                 vz1 = _mm256_add_ps(vx1,vy0);
                _mm256_store_ps(&vec_z[i+8ull],vz1);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
#endif                               
                 vx2 = _mm256_load_ps(&vec_x[i+16ull]);
                 vz2 = _mm256_add_ps(vx2,vy0);
                _mm256_store_ps(&vec_z[i+16ull],vz2);                 
                 vx3 = _mm256_load_ps(&vec_x[i+24ull]);
                 vz3 = _mm256_add_ps(vx3,vy0);
                _mm256_store_ps(&vec_z[i+24ull],vz3);    
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
#endif             
                 vx4 = _mm256_load_ps(&vec_x[i+32ull]);
                 vz4 = _mm256_add_ps(vx4,vy0);
                _mm256_store_ps(&vec_z[i+32ull],vz4);                 
                 vx5 = _mm256_load_ps(&vec_x[i+40ull]);
                 vz5 = _mm256_add_ps(vx5,vy0);
                _mm256_store_ps(&vec_z[i+40ull],vz5);  
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+48ull],_MM_HINT_T0);
#endif
                 vx6 = _mm256_load_ps(&vec_x[i+48ull]);
                 vz6 = _mm256_add_ps(vx6,vy0);
                _mm256_store_ps(&vec_z[i+48ull],vz6);                  
                 vx7 = _mm256_load_ps(&vec_x[i+56ull]);
                 vz7 = _mm256_add_ps(vx7,vy0);
                _mm256_store_ps(&vec_z[i+56ull],vz7);  
           }

           for(i = 0ull; (i+47ull) < n; i += 48ull) 
           {
                /*
                   VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
                */
                __m256 vx0;
                __m256 vz0;
                __m256 vx1;
                __m256 vz1;
                __m256 vx2;
                __m256 vz2;
                __m256 vx3;
                __m256 vz3;
                __m256 vx4;
                __m256 vz4;
                __m256 vx5;
                __m256 vz5;
               
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
#endif                  
                 vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                 vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);                  
                 vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                 vz1 = _mm256_add_ps(vx1,vy0);
                _mm256_store_ps(&vec_z[i+8ull],vz1);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
#endif                               
                 vx2 = _mm256_load_ps(&vec_x[i+16ull]);
                 vz2 = _mm256_add_ps(vx2,vy0);
                _mm256_store_ps(&vec_z[i+16ull],vz2);                 
                 vx3 = _mm256_load_ps(&vec_x[i+24ull]);
                 vz3 = _mm256_add_ps(vx3,vy0);
                _mm256_store_ps(&vec_z[i+24ull],vz3);    
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+32ull],_MM_HINT_T0);
#endif             
                 vx4 = _mm256_load_ps(&vec_x[i+32ull]);
                 vz4 = _mm256_add_ps(vx4,vy0);
                _mm256_store_ps(&vec_z[i+32ull],vz4);                 
                 vx5 = _mm256_load_ps(&vec_x[i+40ull]);
                 vz5 = _mm256_add_ps(vx5,vy0);
                _mm256_store_ps(&vec_z[i+40ull],vz5);              
           }

           for(i = 0ull; (i+31ull) < n; i += 32ull) 
           {
                /*
                   VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
                */
                __m256 vx0;
                __m256 vz0;
                __m256 vx1;
                __m256 vz1;
                __m256 vx2;
                __m256 vz2;
                __m256 vx3;
                __m256 vz3;
               
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
#endif                  
                 vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                 vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);                  
                 vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                 vz1 = _mm256_add_ps(vx1,vy0);
                _mm256_store_ps(&vec_z[i+8ull],vz1);   
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+16ull],_MM_HINT_T0);
#endif                               
                 vx2 = _mm256_load_ps(&vec_x[i+16ull]);
                 vz2 = _mm256_add_ps(vx2,vy0);
                _mm256_store_ps(&vec_z[i+16ull],vz2);                 
                 vx3 = _mm256_load_ps(&vec_x[i+24ull]);
                 vz3 = _mm256_add_ps(vx3,vy0);
                _mm256_store_ps(&vec_z[i+24ull],vz3);                                 
           }

           for(i = 0ull; (i+15ull) < n; i += 16ull) 
           {
                /*
                   VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
                */
                __m256 vx0;
                __m256 vz0;
                __m256 vx1;
                __m256 vz1;
                               
#if (ADD_VEC_YMM8R4_SOFT_PREFETCH) == 1
                _mm_prefetch((const char *)&vec_x[i+0ull],_MM_HINT_T0);
#endif                  
                 vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                 vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);                  
                 vx1 = _mm256_load_ps(&vec_x[i+8ull]);
                 vz1 = _mm256_add_ps(vx1,vy0);
                _mm256_store_ps(&vec_z[i+8ull],vz1);                                                     
           }

           for(i = 0ull; (i+7ull) < n; i += 8ull) 
           {
                /*
                   VADDPS/D, Port=0-1, LAT=4, TPUT=0.5
                */
                __m256 vx0;
                __m256 vz0;
                               
                 vx0 = _mm256_load_ps(&vec_x[i+0ull]);
                 vz0 = _mm256_add_ps(vx0,vy0);
                _mm256_store_ps(&vec_z[i+0ull],vz0);                                                                       
           }

           for(; (i+0ull) < n; i += 1ull) 
           {
                const float x{vec_x[i]};
                const float z{x+y};
                vec_z[i] =  z;
           }

           return (0);
}



