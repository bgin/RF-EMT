
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

#ifndef __GMS_ANTENNA_INTEGRANDS_CUH__
#define __GMS_ANTENNA_INTEGRANDS_CUH__



#include <cuda/std/complex>

using namespace cuda::std;

// Integrand: (2-42)
    
      struct Intf242 {
      
          __device__ complex<float>
          operator()(const complex<float> ikz,
                     const float th,
                     const float z,
                     const float k,
                     const float m,
                     const float I0) {
             
             constexpr float pi = 3.14159265358979323846264f;
             complex<float> earg;
             complex<float> cexp;
             complex<float> value;
             float cost,mpi2,kz,sarg;
             mpi2 = m*pi*0.5f;
             cost = cosf(th);
             kz   = k*z;
             sarg = I0*sinf(kz+mpi2);
             earg = ikz*cost;
             cexp = exp(earg);
             value= sarg*cexp;
             return (value);     
         }
    };
    
    
// Integrand: (2-69)
   
      struct Intf269 {
      
          __device__ float
          operator()(const float M,
                     const float u,
                     const float x,
                     const float L) {
            
            constexpr float pi = 3.14159265358979323846264f;
            float MM,uu,piu,_2xL;
            float chtrm,cost1,cost2,sqr;
            float value;
            MM   = M*M;
            piu  = pi*u;
            cost1= cosf(piu);
            uu   = u*u;
            _2xL = (2.0f*x)/L;
            sqr  = sqrtf(MM-uu);
            chtrm= coshf(pi*sqr);
            cost2= cosf(piu*_2xL);
            value= (chtrm-cost1)*cost2;
            return (value);
        }
   };

      
   

   
   













#endif /*__GMS_ANTENNA_INTEGRANDS_CUH__*/
