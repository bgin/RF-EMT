
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

#ifndef __GMS_CUDA_SVECT_CUH__
#define __GMS_CUDA_SVECT_CUH__


#include <cuda_runtime.h>

/*
   Cuda short vector implementation of mathematical
   operations by defining an overloaded operators.
*/


__forceinline__
__host__ __device__ float2
create_vecf2(const float x) {
    return (make_float2(x,x));
}

__forceinline__
__host__ __device__ float2
create_vecf2(const float x,
             const float y) {
    return (make_float2(x,y));
}

__forceinline__
__host__ __device__ float2
create_vecf2(const float3 v) {
     return (make_float2(v.x,v.y));
}

#define FLOAT4_HI 0

__forceinline__
__host__ __device__ float2
create_vecf2(const float4 v) {
#if (FLOAT4_HI) == 1
     return (make_float2(v.z,v.w)); // higher part
#else
     return (make_float2(v.x,v.y)); // lower part
endif
}

__forceinline__
__host__ __device__ float2
create_vecf2(const int2 v) {
      return (make_float2(float(v.x),float(v.y)));
}

__forceinline__
__host__ __device__ float2
create_vecf2(const uint2 v) {
       return (make_float2(float(v.x),float(v.y)));
}

__foreceinline__
__host__ __device__ int2
create_veci2(const int x) {
       return (make_int2(x,x));
}

__foreceinline__
__host__ __device__ int2
create_veci2(const int x,
             const int y) {
       return (make_int2(x,y));
}

__forceinline__
__host__ __device__ int2
create_veci2(const int2 v) {
       return (make_int2(v));
}

__forceinline__
__host__ __device__ int2
create_veci2(const float2 v) {
     return (make_int2(int(v.x),int(v.y)));
}

__forceinline__
__host__ __device__ int2
create_veci2(const float3 v) {
     return (make_int2(int(v.x),int(v.y)));
}

__forceinline__
__host__ __device__ int2
create_veci2(const float4 v) {
#if (FLOAT4_HI) == 1
     return (make_int2(int(v.z),int(v.w)));
#else
     return (make_int2(int(v.x),int(v.y)));
#endif
}

__forceinline__
__host__ __device__ float3
create_vecf3(const float x) {
      return (make_float3(x,x,x));
}

__forceinline__
__host__ __device__ float3
create_vecf3(const float x,
             const float y) {
      return (make_float3(x,y,0.0f));
}

__forceinline__
__host__ __device__ float3
create_vecf3(const float x,
             const float y,
	     const float z) {
      return (make_float3(x,y,z));
}

__forceinline__
__host__ __device__ float3
create_vecf3(const float2 v,
             const float  z) {
      return (make_float3(v.x,v.y,z));
}

__forceinline__
__host__ __device__ float3
create_vecf3(const float3 v) {
       return (make_float3(v.x,v.y,v.z));
}

__forceinline__
__host__ __device__ float3
create_vecf3(const float4 v) {
       return (make_float3(v.x,x.y,v.z));
}

__forceinline__
__host__ __device__ float4
create_vecf4(const float x) {
       return (make_float4(x,x,x,x));
}

__forceinline__
__host__ __device__ float4
create_vecf4(const float2 v,
             const float z,
	     const float w) {
        return (make_float4(v.x,v.y,z,w));
}

__forceinline__
__host__ __device__ float4
create_vecf4(const float3 v,
             const float w) {
        return (make_float4(v.x,v.y,v.z,w));
}

__forceinline__
__host__ __device__ float4
create_vecf4(const float4 v) {
       return (make_float4(v.x,v.y,v.z,v.w));
}

__forceinline__
__host__ __device__ float4
create_vecf4(const int4 v) {
       return (make_float4(float(v.x),float(v.y),
                           float(v.z),float(v.w)));
}


__forceinline__
__host__ __device__ float2
operator+(const float2 v1,
          const float2 v2) {
   return (make_float2(v1.x+v2.x,
                       v1.y+v2.y));
}

__forceinline__
__host__ __device__ float2
operator+(const float2 v,
          const float y) {
    return (make_float2(v.x+y,v.y+y));
}

__forceinline__
__host__ __device__ float2
operator+(const float y,
          const float2 v) {
     return (make_float2(v.x+y,v.y+y));
}

__forceinline__
__host__ __device__ float2
operator+=(float2 &v1,
           const float2 v2) {
   v1.x += v2.x;
   v1.y += v2.y;
}

__forceinline__
__host__ __device__ float2
operator+=(float2 &v1,
           const float x) {
   v1.x += x;
   v1.y += x;
}

__forceinline__
__host__ __device__ float2
operator+=(const float x,
           float2 &v1) {
   v1.x += x;
   v1.y += x;
}


//Skipping integral operators (no much usage of those operations)

__forceinline__
__host__ __device__ float3
operator+(const float3 v1,
          const float3 v2) {
      return (make_float3(v1.x+v2.x,
                          v1.y+v2.y,
			  v1.z+v2.z));
}

__forceinline__
__host__ __device__ float3
operator+(const float3 v1,
          const float  x) {
       return (make_float3(v1.x+x,
                           v1.y+x,
			   v1.z+x));
}

__forceinline__
__host__ __device__ float3
operator+(const float x,
          const float3 v1) {
       return (make_float3(v1.x+x,
                           v1.y+x,
			   v1.z+x));
}

__forceinline__
__host__ __device__ float3
operator+=(float3 &v1,
           float3 v2) {
   v1.x += v2.x;
   v1.y += v2.y;
   v1.z += v2.z;
}

__forceinline__
__host__ __device__ float3
operator+=(float3 &v1,
           const float x) {
   v1.x += x;
   v1.y += x;
   v1.z += x;
}

__forceinline__
__host__ __device__ float3
operator+=(const float x,
           float3 v1) {
   v1.x += x;
   v1.y += x;
   v1.z += x;
}

__forceinline__
__host__ __device__ float4
operator+(const float4 v1,
          const float4 v2) {
    return (make_float4(v1.x+v2.x,
                        v1.y+v2.y,
			v1.z+v2.z,
			v1.w+v2.w));
}

__forceinline__
__host__ __device__ float4
operator+(const float4 v1,
          const float x) {
     return (make_float4(v1.x+x,
                         v1.y+x,
			 v1.z+x,
			 v1.w+x));
}

__forceinline__
__host__ __device__ float4
operator+(const float x,
          const float4 v1) {
      return (make_float4(v1.x+x,
                         v1.y+x,
			 v1.z+x,
			 v1.w+x));
}

__forceinline__
__host__ __device__ float4
operator+=(float4 &v1,
           const float4 v2) {
   v1.x += v2.x;
   v1.y += v2.y;
   v1.z += v2.z;
   v1.w += v2.w;
}

__forceinline__
__host__ __device__ float4
operator+=(float4 &v1,
           const float x) {
   v1.x += x;
   v1.y += x;
   v1.z += x;
   v1.w += x;
}

__forceinline__
__host__ __device__ float4
operator+=(const float x,
           float4 &v1) {
   v1.x += x;
   v1.y += x;
   v1.z += x;
   v1.w += x;
}


//**********************************************************//
//**********************************************************//

__forceinline__
__host__ __device__ float2
operator-(const float2 v1,
          const float2 v2) {
   return (make_float2(v1.x-v2.x,
                       v1.y-v2.y));
}

__forceinline__
__host__ __device__ float2
operator-(const float2 v,
          const float y) {
    return (make_float2(v.x-y,v.y-y));
}

__forceinline__
__host__ __device__ float2
operator-(const float y,
          const float2 v) {
     return (make_float2(v.x-y,v.y-y));
}

__forceinline__
__host__ __device__ float2
operator-=(float2 &v1,
           const float2 v2) {
   v1.x -= v2.x;
   v1.y -= v2.y;
}

__forceinline__
__host__ __device__ float2
operator-=(float2 &v1,
           const float x) {
   v1.x -= x;
   v1.y -= x;
}

__forceinline__
__host__ __device__ float2
operator-=(const float x,
           float2 &v1) {
   v1.x -= x;
   v1.y -= x;
}


//Skipping integral operators (no much usage of those operations)

__forceinline__
__host__ __device__ float3
operator-(const float3 v1,
          const float3 v2) {
      return (make_float3(v1.x-v2.x,
                          v1.y-v2.y,
			  v1.z-v2.z));
}

__forceinline__
__host__ __device__ float3
operator-(const float3 v1,
          const float  x) {
       return (make_float3(v1.x-x,
                           v1.y-x,
			   v1.z-x));
}

__forceinline__
__host__ __device__ float3
operator-(const float x,
          const float3 v1) {
       return (make_float3(v1.x-x,
                           v1.y-x,
			   v1.z-x));
}

__forceinline__
__host__ __device__ float3
operator-=(float3 &v1,
           float3 v2) {
   v1.x -= v2.x;
   v1.y -= v2.y;
   v1.z -= v2.z;
}

__forceinline__
__host__ __device__ float3
operator-=(float3 &v1,
           const float x) {
   v1.x -= x;
   v1.y -= x;
   v1.z -= x;
}

__forceinline__
__host__ __device__ float3
operator-=(const float x,
           float3 v1) {
   v1.x -= x;
   v1.y -= x;
   v1.z -= x;
}

__forceinline__
__host__ __device__ float4
operator-(const float4 v1,
          const float4 v2) {
    return (make_float4(v1.x-v2.x,
                        v1.y-v2.y,
			v1.z-v2.z,
			v1.w-v2.w));
}

__forceinline__
__host__ __device__ float4
operator-(const float4 v1,
          const float x) {
     return (make_float4(v1.x-x,
                         v1.y-x,
			 v1.z-x,
			 v1.w-x));
}

__forceinline__
__host__ __device__ float4
operator-(const float x,
          const float4 v1) {
      return (make_float4(v1.x-x,
                         v1.y-x,
			 v1.z-x,
			 v1.w-x));
}

__forceinline__
__host__ __device__ float4
operator-=(float4 &v1,
           const float4 v2) {
   v1.x -= v2.x;
   v1.y -= v2.y;
   v1.z -= v2.z;
   v1.w -= v2.w;
}

__forceinline__
__host__ __device__ float4
operator-=(float4 &v1,
           const float x) {
   v1.x -= x;
   v1.y -= x;
   v1.z -= x;
   v1.w -= x;
}

__forceinline__
__host__ __device__ float4
operator-=(const float x,
           float4 &v1) {
   v1.x -= x;
   v1.y -= x;
   v1.z -= x;
   v1.w -= x;
}

//****************************************************//
//****************************************************//

__forceinline__
__host__ __device__ float2
operator*(const float2 v1,
          const float2 v2) {
   return (make_float2(v1.x*v2.x,
                       v1.y*v2.y));
}

__forceinline__
__host__ __device__ float2
operator*(const float2 v,
          const float y) {
    return (make_float2(v.x*y,v.y*y));
}

__forceinline__
__host__ __device__ float2
operator*(const float y,
          const float2 v) {
     return (make_float2(v.x*y,v.y*y));
}

__forceinline__
__host__ __device__ float2
operator*=(float2 &v1,
           const float2 v2) {
   v1.x *= v2.x;
   v1.y *= v2.y;
}

__forceinline__
__host__ __device__ float2
operator*=(float2 &v1,
           const float x) {
   v1.x *= x;
   v1.y *= x;
}

__forceinline__
__host__ __device__ float2
operator*=(const float x,
           float2 &v1) {
   v1.x *= x;
   v1.y *= x;
}


//Skipping integral operators (no much usage of those operations)

__forceinline__
__host__ __device__ float3
operator*(const float3 v1,
          const float3 v2) {
      return (make_float3(v1.x*v2.x,
                          v1.y*v2.y,
			  v1.z*v2.z));
}

__forceinline__
__host__ __device__ float3
operator*(const float3 v1,
          const float  x) {
       return (make_float3(v1.x*x,
                           v1.y*x,
			   v1.z*x));
}

__forceinline__
__host__ __device__ float3
operator*(const float x,
          const float3 v1) {
       return (make_float3(v1.x*x,
                           v1.y*x,
			   v1.z*x));
}

__forceinline__
__host__ __device__ float3
operator*=(float3 &v1,
           float3 v2) {
   v1.x *= v2.x;
   v1.y *= v2.y;
   v1.z *= v2.z;
}

__forceinline__
__host__ __device__ float3
operator*=(float3 &v1,
           const float x) {
   v1.x *= x;
   v1.y *= x;
   v1.z *= x;
}

__forceinline__
__host__ __device__ float3
operator*=(const float x,
           float3 v1) {
   v1.x *= x;
   v1.y *= x;
   v1.z *= x;
}

__forceinline__
__host__ __device__ float4
operator*(const float4 v1,
          const float4 v2) {
    return (make_float4(v1.x*v2.x,
                        v1.y*v2.y,
			v1.z*v2.z,
			v1.w*v2.w));
}

__forceinline__
__host__ __device__ float4
operator*(const float4 v1,
          const float x) {
     return (make_float4(v1.x*x,
                         v1.y*x,
			 v1.z*x,
			 v1.w*x));
}

__forceinline__
__host__ __device__ float4
operator*(const float x,
          const float4 v1) {
      return (make_float4(v1.x*x,
                         v1.y*x,
			 v1.z*x,
			 v1.w*x));
}

__forceinline__
__host__ __device__ float4
operator*=(float4 &v1,
           const float4 v2) {
   v1.x *= v2.x;
   v1.y *= v2.y;
   v1.z *= v2.z;
   v1.w *= v2.w;
}

__forceinline__
__host__ __device__ float4
operator*=(float4 &v1,
           const float x) {
   v1.x *= x;
   v1.y *= x;
   v1.z *= x;
   v1.w *= x;
}

__forceinline__
__host__ __device__ float4
operator*=(const float x,
           float4 &v1) {
   v1.x *= x;
   v1.y *= x;
   v1.z *= x;
   v1.w *= x;
}


//***************************************************//
//***************************************************//


__forceinline__
__host__ __device__ float2
operator/(const float2 v1,
          const float2 v2) {
   return (make_float2(v1.x/v2.x,
                       v1.y/v2.y));
}

__forceinline__
__host__ __device__ float2
operator/(const float2 v,
          const float y) {
    return (make_float2(v.x/y,v.y/y));
}

__forceinline__
__host__ __device__ float2
operator/(const float y,
          const float2 v) {
     return (make_float2(v.x/y,v.y/y));
}

__forceinline__
__host__ __device__ float2
operator/=(float2 &v1,
           const float2 v2) {
   v1.x /= v2.x;
   v1.y /= v2.y;
}

__forceinline__
__host__ __device__ float2
operator/=(float2 &v1,
           const float x) {
   v1.x /= x;
   v1.y /= x;
}

__forceinline__
__host__ __device__ float2
operator/=(const float x,
           float2 &v1) {
   v1.x /= x;
   v1.y /= x;
}


//Skipping integral operators (no much usage of those operations)

__forceinline__
__host__ __device__ float3
operator/(const float3 v1,
          const float3 v2) {
      return (make_float3(v1.x/v2.x,
                          v1.y/v2.y,
			  v1.z/v2.z));
}

__forceinline__
__host__ __device__ float3
operator/(const float3 v1,
          const float  x) {
       return (make_float3(v1.x/x,
                           v1.y/x,
			   v1.z/x));
}

__forceinline__
__host__ __device__ float3
operator/(const float x,
          const float3 v1) {
       return (make_float3(v1.x/x,
                           v1.y/x,
			   v1.z/x));
}

__forceinline__
__host__ __device__ float3
operator/=(float3 &v1,
           float3 v2) {
   v1.x /= v2.x;
   v1.y /= v2.y;
   v1.z /= v2.z;
}

__forceinline__
__host__ __device__ float3
operator/=(float3 &v1,
           const float x) {
   v1.x /= x;
   v1.y /= x;
   v1.z /= x;
}

__forceinline__
__host__ __device__ float3
operator/=(const float x,
           float3 v1) {
   v1.x /= x;
   v1.y /= x;
   v1.z /= x;
}

__forceinline__
__host__ __device__ float4
operator/(const float4 v1,
          const float4 v2) {
    return (make_float4(v1.x/v2.x,
                        v1.y/v2.y,
			v1.z/v2.z,
			v1.w/v2.w));
}

__forceinline__
__host__ __device__ float4
operator/(const float4 v1,
          const float x) {
     return (make_float4(v1.x/x,
                         v1.y/x,
			 v1.z/x,
			 v1.w/x));
}

__forceinline__
__host__ __device__ float4
operator/(const float x,
          const float4 v1) {
      return (make_float4(v1.x/x,
                         v1.y/x,
			 v1.z/x,
			 v1.w/x));
}

__forceinline__
__host__ __device__ float4
operator/=(float4 &v1,
           const float4 v2) {
   v1.x /= v2.x;
   v1.y /= v2.y;
   v1.z /= v2.z;
   v1.w /= v2.w;
}

__forceinline__
__host__ __device__ float4
operator/=(float4 &v1,
           const float x) {
   v1.x /= x;
   v1.y /= x;
   v1.z /= x;
   v1.w /= x;
}

__forceinline__
__host__ __device__ float4
operator/=(const float x,
           float4 &v1) {
   v1.x /= x;
   v1.y /= x;
   v1.z /= x;
   v1.w /= x;
}

//***************************************************//
//***************************************************//

__forceinline__
__host__ __device__ float2
dotp(const float2 v1,
     const float2 v2) {
  return (v1.x*v2.x+v1.y*v2.y);
}

__forceinline__
__host__ __device__ float3
dotp(const float3 v1,
     const float3 v2) {
   return (v1.x*v2.x+v1.y*v2.y+v1.z*v2.z);
}

__forceinline__
__host__ __device__ float4
dotp(const float4 v1,
     const float4 v2) {
   return (v1.x*v2.x+v1.y*v2.y+v1.z*v2.z+v1.w*v2.w);
}

//=================================================//
//================================================//

__forceinline__
__host__ __device__ float
vmag(const float2 v) {
  return (sqrtf(dotp(v,v)));
}

__forceinline__
__host__ __device__ float
vmag(const float3 v) {
  return (sqrtf(dotp(v,v)));
}

__forceinline__
__host__ __device__ float
vmag(const float4 v) {
  return (sqrtf(dotp(v,v)));
}

//***************************************************//
//**************************************************//

__forceinline__
__host__ __device__ float2
vabs(const float2 v) {
  return (make_float2(fabs(v.x),fabs(v.y)));
}

__forceinline__
__host__ __device__ float3
vabs(const float3 v) {
  return (make_float3(fabs(v.x),fabs(v.y),fabs(v.z)));
}

__forceinline__
__host__ __device__ float4
vabs(const float4 v) {
  return (make_float4(fabs(v.x),fabs(v.y),
                      fabs(v.z),fabs(v.w))));
}

//===================================================//
//===================================================//

__forceinline__
__host__ __device__ float3
crossp(const float3 v1,
       const float3 v2) {

   const float c0 = v1.y*v2.z-v1.z*v2.y;
   const float c1 = v1.z*v2.x-v1.x*v2.z;
   const float c2 = v1.x*v2.y-v1.y*v2.x;
   return (make_float3(c0,c1,c2));
}

//=====================================================//
//=====================================================//

__forceinline__
__host__ __device__ float2
vnormalize(const float2 v) {
      const float inv = rsqrtf(dotp(v,v));
      return (v*inv);
}

__forceinline__
__host__ __device__ float3
vnormalize(const float3 v) {
      const float inv = rsqrtf(dotp(v,v));
      return (v*inv);
}

__forceinline__
__host__ __device__ float4
vnormalize(const float4 v) {
      const float inv = rsqrtf(dotp(v,v));
      return (v*inv);
}

__forceinline__
__host__ __device__ float
dot3f4(const float4 v1,
       const float4 v2) {
   const float4 a = make_float4(v1.x,v1.y,v1.z,0.0f);
   const float4 b = make_float4(v2.x,v2.y,v2.z,0.0f);
   return (dotp(a,b));
}


//*****************************************************//
// Host/Device matrix 3x3 implementation
//*****************************************************//

const __host__ __device__ float4 zero = make_float4(0.0f,0.0f,0.0f,0.0f);
typedef struct __align__(16) {
   
  float4 row3[3];
} Mat3x3;

__forceinline__
__host__ __device__ void
mat3x3_set_zero(Mat3x3 &mat) {
     
     mat.row3[0] = zero;
     mat.row3[1] = zero;
     mat.row3[2] = zero;
}

__forceinline__
__host__ __device__ Mat3x3
mat3x3_set_zero() {
     Mat3x3 mat;
     mat.row3[0] = zero;
     mat.row3[1] = zero;
     mat.row3[2] = zero;
     return (mat);
}

__forceinline__ 
__host__ __device__ void
mat3x3_set_v1(const float4 a,
              const float4 b,
              const float4 c,
              Mat3x3 &mat) {
    
     mat.row3[0] = a;
     mat.row3[1] = b;
     mat.row3[2] = c;
}

__forceinline__ 
__host__ __device__ Mat3x3
mat3x3_set_v1(const float4 a,
              const float4 b,
              const float4 c) {
 
     Mat3x3 mat;
     mat.row3[0] = a;
     mat.row3[1] = b;
     mat.row3[2] = c;
     return (mat);
}



__forceinline__
__host__ __device__ void
mat3x3_identity(Mat3x3 &mat) {
  
     mat.row3[0] = make_float4(1.0f,0.0f,0.0f,0.0f);
     mat.row3[1] = make_float4(0.0f,1.0f,0.0f,0.0f);
     mat.row3[2] = make_float4(0.0f,0.0f,1.0f,0.0f);
}

__forceinline__
__host__ __device__ Mat3x3
mat3x3_identity() {
    
     Mat3x3 mat;
     mat.row3[0] = make_float4(1.0f,0.0f,0.0f,0.0f);
     mat.row3[1] = make_float4(0.0f,1.0f,0.0f,0.0f);
     mat.row3[2] = make_float4(0.0f,0.0f,1.0f,0.0f);
     return (mat);
}

__forceinline__
__host__ __device__ void
mat3x3_diagonal(const float x,
                const float y,
                const float z,
                Mat3x3 &mat) {
    
    mat.row3[0] = make_float4(x,0.f,0.f,0.f);
    mat.row3[1] = make_float4(0.f,y,0.f,0.f);
    mat.row3[2] = make_float4(0.f,0.f,z,0.f);  
}

__forceinline__
__host__ __device__ Mat3x3
mat3x3_diagonal(const float x,
                const float y,
                const float z) {
  
    Mat3x3 mat;
    mat.row3[0] = make_float4(x,0.f,0.f,0.f);
    mat.row3[1] = make_float4(0.f,y,0.f,0.f);
    mat.row3[2] = make_float4(0.f,0.f,z,0.f);  
    return (mat);
}

__forceinline__
__host__ __device__ void
mat3x3_transpose(const Mat3x3 &in) {

   Mat3x3 mat;
   mat.row3[0] = make_float4(in.row3[0].x,in.row3[1].x,in.row3[2].x,0.0f);
   mat.row3[1] = make_float4(in.row3[0].y,in.row3[1].y,in.row3[2].y,0.0f);
   mat.row3[2] = make_float4(in.row3[0].z,in.row3[2].z,in.row3[2].z,0.0f);
   return (mat);
}

__forceinline__
__host__ __device__ Mat3x3
mat3x3_transpose(Mat3x3 &mat,
                 const Mat3x3 &in) {
   
   mat.row3[0] = make_float4(in.row3[0].x,in.row3[1].x,in.row3[2].x,0.0f);
   mat.row3[1] = make_float4(in.row3[0].y,in.row3[1].y,in.row3[2].y,0.0f);
   mat.row3[2] = make_float4(in.row3[0].z,in.row3[2].z,in.row3[2].z,0.0f);
}

__forceinline__
__host__ __device__ void
mat3x3_mul_mat3x3(Mat3x3 &mat,
                  const Mat3x3 &m1,
                  const Mat3x3 &m2) {

    Mat3x3 tmp    = mat3x3_transpose(m2);
    for(int i = 0; i != 3; ++i) {
        mat.row3[i].x = dot3f4(m1.row3[i],tmp.row3[0]);
        mat.row3[i].y = dot3f4(m1.row3[i],tmp.row3[1]);
        mat.row3[i].z = dot3f4(m1.row3[i],tmp.row3[2]);
        mat.row3[i].w = 0.0f;
    }
}

__forceinline__
__host__ __device__ Mat3x3
mat3x3_mul_mat3x3(const Mat3x3 &m1,
                  const Mat3x3 &m2) {
     
    Mat3x3 mat;
    Mat3x3 tmp    = mat3x3_transpose(m2);
    for(int i = 0; i != 3; ++i) {
        mat.row3[i].x = dot3f4(m1.row3[i],tmp.row3[0]);
        mat.row3[i].y = dot3f4(m1.row3[i],tmp.row3[1]);
        mat.row3[i].z = dot3f4(m1.row3[i],tmp.row3[2]);
        mat.row3[i].w = 0.0f;
    }
   return (mat);
}

__forceinline__
__host__ __device__ void
mat3x3_mul_float4(float4 &v4,
                  const Mat3x3 &in,
                  const float4 v) {
   
    v4.x = dot3f4(in.row3[0],v);
    v4.y = dot3f4(in.row3[1],v);
    v4.z = dot3f4(in.row3[2],v);
    v4.w = 0.0f;
}

__forceinline__
__host__ __device__ float4
mat3x3_mul_float4(const Mat3x3 &in,
                  const float4 v) {
   
    float4 v4;
    v4.x = dot3f4(in.row3[0],v);
    v4.y = dot3f4(in.row3[1],v);
    v4.z = dot3f4(in.row3[2],v);
    v4.w = 0.0f;
    return (v4);
}

__forceinline__
__host__ __device__ void
mat3x3_mul_float(Mat3x3 &mat,
                 const Mat3x3 &in,
                 const float s) {
    
    mat.row3[0] = s*in.row3[0];
    mat.row3[1] = s*in.row3[1];
    mat.row3[2] = s*in.row3[2];
}

__forceinline__
__host__ __device__ Mat3x3
mat3x3_mul_float(const Mat3x3 &in,
                 const float s) {

     Mat3x3 mat;
     mat.row3[0] = s*in.row3[0];
     mat.row3[1] = s*in.row3[1];
     mat.row3[2] = s*in.row3[2];
     return (mat);
}















/*
       Example of potential usage
#include <thrust/device_vector.h>

#define BLOCKSIZE 256


int iDivUp(int a, int b){ return ((a % b) != 0) ? (a / b + 1) : (a / b); }


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


__global__ void add_float(float *d_a, float *d_b, float *d_c, unsigned int N) {

    const int tid = 4 * threadIdx.x + blockIdx.x * (4 * blockDim.x);

    if (tid < N) {

        float a1 = d_a[tid];
        float b1 = d_b[tid];

        float a2 = d_a[tid+1];
        float b2 = d_b[tid+1];

        float a3 = d_a[tid+2];
        float b3 = d_b[tid+2];

        float a4 = d_a[tid+3];
        float b4 = d_b[tid+3];

        float c1 = a1 + b1;
        float c2 = a2 + b2;
        float c3 = a3 + b3;
        float c4 = a4 + b4;

        d_c[tid] = c1;
        d_c[tid+1] = c2;
        d_c[tid+2] = c3;
        d_c[tid+3] = c4;

        //if ((tid < 1800) && (tid > 1790)) {
            //printf("%i %i %i %f %f %f\n", tid, threadIdx.x, blockIdx.x, a1, b1, c1);
            //printf("%i %i %i %f %f %f\n", tid+1, threadIdx.x, blockIdx.x, a2, b2, c2);
            //printf("%i %i %i %f %f %f\n", tid+2, threadIdx.x, blockIdx.x, a3, b3, c3);
            //printf("%i %i %i %f %f %f\n", tid+3, threadIdx.x, blockIdx.x, a4, b4, c4);
        //}

    }

}


__global__ void add_float2(float2 *d_a, float2 *d_b, float2 *d_c, unsigned int N) {

    const int tid = 2 * threadIdx.x + blockIdx.x * (2 * blockDim.x);

    if (tid < N) {

        float2 a1 = d_a[tid];
        float2 b1 = d_b[tid];

        float2 a2 = d_a[tid+1];
        float2 b2 = d_b[tid+1];

        float2 c1;
        c1.x = a1.x + b1.x;
        c1.y = a1.y + b1.y;

        float2 c2;
        c2.x = a2.x + b2.x;
        c2.y = a2.y + b2.y;

        d_c[tid] = c1;
        d_c[tid+1] = c2;

    }

}


__global__ void add_float4(float4 *d_a, float4 *d_b, float4 *d_c, unsigned int N) {

    const int tid = 1 * threadIdx.x + blockIdx.x * (1 * blockDim.x);

    if (tid < N/4) {

        float4 a1 = d_a[tid];
        float4 b1 = d_b[tid];

        float4 c1;
        c1.x = a1.x + b1.x;
        c1.y = a1.y + b1.y;
        c1.z = a1.z + b1.z;
        c1.w = a1.w + b1.w;

        d_c[tid] = c1;

    }

}


int main() {

    const int N = 4*10000000;

    const float a = 3.f;
    const float b = 5.f;

    // --- float

    thrust::device_vector<float> d_A(N, a);
    thrust::device_vector<float> d_B(N, b);
    thrust::device_vector<float> d_C(N);

    float time;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    add_float<<<iDivUp(N/4, BLOCKSIZE), BLOCKSIZE>>>(thrust::raw_pointer_cast(d_A.data()), thrust::raw_pointer_cast(d_B.data()), thrust::raw_pointer_cast(d_C.data()), N);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    printf("Elapsed time:  %3.1f ms \n", time); gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    thrust::host_vector<float> h_float = d_C;
    for (int i=0; i<N; i++) {
        if (h_float[i] != (a+b)) {
            printf("Error for add_float at %i: result is %f\n",i, h_float[i]);
            return -1;
        }
    }

    // --- float2

    thrust::device_vector<float> d_A2(N, a);
    thrust::device_vector<float> d_B2(N, b);
    thrust::device_vector<float> d_C2(N);

    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    add_float2<<<iDivUp(N/4, BLOCKSIZE), BLOCKSIZE>>>((float2*)thrust::raw_pointer_cast(d_A2.data()), (float2*)thrust::raw_pointer_cast(d_B2.data()), (float2*)thrust::raw_pointer_cast(d_C2.data()), N);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    printf("Elapsed time:  %3.1f ms \n", time); gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    thrust::host_vector<float> h_float2 = d_C2;
    for (int i=0; i<N; i++) {
        if (h_float2[i] != (a+b)) {
            printf("Error for add_float2 at %i: result is %f\n",i, h_float2[i]);
            return -1;
        }
    }

    // --- float4

    thrust::device_vector<float> d_A4(N, a);
    thrust::device_vector<float> d_B4(N, b);
    thrust::device_vector<float> d_C4(N);

    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    add_float4<<<iDivUp(N/4, BLOCKSIZE), BLOCKSIZE>>>((float4*)thrust::raw_pointer_cast(d_A4.data()), (float4*)thrust::raw_pointer_cast(d_B4.data()), (float4*)thrust::raw_pointer_cast(d_C4.data()), N);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    printf("Elapsed time:  %3.1f ms \n", time); gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    thrust::host_vector<float> h_float4 = d_C4;
    for (int i=0; i<N; i++) {
        if (h_float4[i] != (a+b)) {
            printf("Error for add_float4 at %i: result is %f\n",i, h_float4[i]);
            return -1;
        }
    }

    return 0;
}

*/














#endif /*__GMS_CUDA_SVECT_CUH__*/
