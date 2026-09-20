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

#ifndef __GMS_GPU_HELPERS_CUH__
#define __GMS_GPU_HELPERS_CUH__

/*
    Various GPU short helper functions (mainly computational).
    Partly adapted from various sources.
*/

// Atomic-add double precision (CUDA-C programming guide)
__device__ double
atomicAdd_r8_1(double * vaddr,
          const double value) {
   const unsigned long long int * __restrict__ to_int_rep = (unsigned long long int*)vaddr;
   unsigned long long int prev_vaddr = *to_int_rep;
   unsigned long long int t_vaddr_val = 0ULL;
   do {
        t_vaddr_val = prev_vaddr;
	prev_vaddr  = atomicCas(to_int_rep,t_vaddr_val,
	                        __double_as_longlong(value+
				                     __longlong_as_double(t_vaddr_val)));
   } while(t_vaddr_val != prev_vaddr);
   return (__longlong_as_double(prev_vaddr));
}


// Atomic-min single precision
__device__ float
atomicMin_r4_1(float * vaddr,
               const float value) {
    unsigned int * __restrict__ to_int_rep = (unsigned int*)vaddr;
    unsigned int prev_vaddr  = *to_int_rep;
    unsigned int t_vaddr_val = 0;
    do {
         t_vaddr_val = prev_vaddr;
	 prev_vaddr = atomicCas(to_int_rep,t_vaddr_val,
	                        __float_as_int(fminf(value,__int_as_float(prev_vaddr))));
	 
    }while(t_vaddr_val != prev_vadr);
    return (__int_as_float(prev_vaddr));
}


// Atomic-max single precision
__device__ float
atomicMax_r4_1(float * vaddr,
               const float value) {
    unsigned int * __restrict__ to_int_rep = (unsigned int*)vaddr;
    unsigned int prev_vaddr  = *to_int_rep;
    unsigned int t_vaddr_val = 0;
    do {
         t_vaddr_val = prev_vaddr;
	 prev_vaddr = atomicCas(to_int_rep,t_vaddr_val,
	                        __float_as_int(fmaxf(value,__int_as_float(prev_vaddr))));
	 
    }while(t_vaddr_val != prev_vadr);
    return (__int_as_float(prev_vaddr));
}


//Fast approximation of ERFC by Norbert Juff, NVIDIA.
__device__ __forceinline__
float erfc_r4_1(const float x) {
    float t = 0.0f;
    t =                       (float)-1.6488499458192755E-006;
    t = __internal_fmad(t, x, (float)2.9524665006554534E-005);
    t = __internal_fmad(t, x, (float)-2.3341951153749626E-004);
    t = __internal_fmad(t, x, (float)1.0424943374047289E-003);
    t = __internal_fmad(t, x, (float)-2.5501426008983853E-003);
    t = __internal_fmad(t, x, (float)3.1979939710877236E-004);
    t = __internal_fmad(t, x, (float)2.7605379075746249E-002);
    t = __internal_fmad(t, x, (float)-1.4827402067461906E-001);
    t = __internal_fmad(t, x, (float)-9.1844764013203406E-001);
    t = __internal_fmad(t, x, (float)-1.6279070384382459E+000);
    t = t * x;
    return exp2f(t);
}

/*
    Computation of global thread indexing
*/


// 1D grid of 1D blocks
__device__ __forceinline__
int globalIdx_1D_1D() {
    return (blockIdx.x*blockDim.x+threadIdx.x);
}

// 1D grid of 2D blocks
__device__ __forceinline__
int globalIdx_1D_2D() {
    return (blockIdx.x*blockDim.x*blockDim.y+
             threadIdx.y*blockDim.x+threadIdx.x);
}

// 1D grid of 3D blocks
__device__ __forceinline__
int globalIdx_1D_3D() {
    return (blockIdx.x*blockDim.x*blockDim.y*blockDim.z+
            threadIdx.z*blockDim.y*blockDim.x+
            threadIdx.y*blockDim.x+threadIdx.x);
}

// 2D grid of 1D blocks
__device__ __forceinline__
int globalIdx_2D_1D() {
     int blockId  = blockIdx.y*gridDim.x+blockIdx.x;
     int threadId = blockId*blockDim.x+threadIdx.x;
     return (threadId);
}

// 2D grid of 2D blocks
__device__ __forceinline__
int globalIdx_2D_2D() {
    int blockId  = blockIdx.x+blockIdx.y*gridDim.x;
    int threadId = blockId*(blockDim.x*blockDim.y)+
                   (threadIdx.y*blockDim.x)+threadIdx.x;
    return (threadId);
}

// 2D grid of 3D blocks
__device__ __forceinline__
int globalIdx_2D_3D() {
    int blockId  =  blockIdx.x+blockIdx.y*gridDim.x;
    int threadId =  blockId*(blockDim.x*blockDim.y*blockDim.z)+ 
                   (threadIdx.z* (blockDim.x*blockDim.y)) +
                   (threadIdx.y* blockDim.x)+threadIdx.x;
    return (threadId);
}

// 3D grid of 1D blocks
__device__ __forceinline__
int globalIdx_3D_1D() {
    int blockId  = blockIdx.x+blockIdx.y*gridDim.x+
                   gridDim.x*gridDim.y*blockIdx.z;
    int threadId = blockId*blockDim.x+threadIdx.x;
    return (threadId);
}

// 3D grid of 2D blocks
__device__ __forceinline__
int globalIdx_3D_2D() {
    int blockId  = blockIdx.x+blockIdx.y*gridDim.x+
                   gridDim.x*gridDim.y*blockIdx.z;
    int threadId = blockId*(blockDim.x* blockDim.y)+
                   (threadIdx.y*blockDim.x)+threadIdx.x;
    return (threadId);
}

// 3D grid of 3D blocks
__device__ __forceinline__
int globalIdx_3D_3D() {
     int blockId  = blockIdx.x+blockIdx.y*gridDim.x+
                    gridDim.x*gridDim.y*blockIdx.z;
     int threadId = blockId*(blockDim.x*blockDim.y*blockDim.z)+
                    (threadIdx.z * (blockDim.x * blockDim.y)) +
                    (threadIdx.y * blockDim.x)  + threadIdx.x;
     return (threadId);
}

#define PXL_GLOBAL_PTR   "r"

// Software prefetching.
inline __device__ 
void __prefetch_global_l1(const void* const ptr)
{
  asm("prefetch.global.L1 [%0];" : : PXL_GLOBAL_PTR(ptr));
}

inline __device__
void __prefetch_global_uniform(const void* const ptr)
{
  asm("prefetchu.L1 [%0];" : : PXL_GLOBAL_PTR(ptr));
}

inline __device__ 
void __prefetch_global_l2(const void* const ptr)
{
  asm("prefetch.global.L2 [%0];" : : PXL_GLOBAL_PTR(ptr));
}



/*

     Simple host device complex arithmetic auxiliary routines.
*/

#include <cuComplex.h>



__host__ __device__
__forceinline__ cuFloatComplex
create_complex_c4_1(const float re,
		  const float im) {
  return (make_cuFloatComplex(re,im));
}


__host__ __device__
__forceinline__ cuFloatDouble
create_complex_c8_1(const double re,
		  const double im) {
  return (make_cuDoubleComplex(re,im));
}


__host__ __device__
__forceinline__ float
real_c4_1(const cuFloatComplex c) {
    return (cuCrealf(c));
}


__host__ __device__
__forceinline__ double
real_c8_1(const cuDoubleComplex c) {
    return (cuCreal(c));
}


__host__ __device__
__forceinline__ float
imag_c4_1(const cuFloatComplex c) {
     return (cuCimagf(c));
}


__host__ __device__
__forceinline__ double
imag_c8_1(const cuDoubleComplex c) {
     return (cuCimag(c));
}


__host__ __device__
__forceinline__ cuFloatComplex
conjugate_c4_1(const cuFloatComplex c) {
      return (cuConjf(c));
}


__host__ __device__
__forceinline__ cuDoubleComplex
conjugate_c8_1(const cuDoubleComplex c) {
      return (cuConj(c));
}


__host__ __device__
__forceinline__ float
cabs_c4_1(const cuFloatComplex c) {
       return (cuCabsf(c));
}


__host__ __device__
__forceinline__ double
cabs_c8_1(const cuDoubleComplex c) {
       return (cuCabs(c));
}


__host__ __device__
__forceinline__ cuFloatComplex
cdiv_c4_1(const cuFloatComplex c1,
          const cuFloatComplex c2) {
       return (cuCdivf(c1,c2));
}


__host__ __device__
__forceinline__ cuDoubleComplex
cdiv_c8_1(const cuDoubleComplex c1,
          const cuDoubleComplex c2) {
       return (cuCdiv(c1,c2));
}


__host__ __device__
__forceinline__ cuFloatComplex
cmul_c4_1(const cuFloatComplex c1,
          const cuFloatComplex c2) {
       return (cuCmulf(c1,c2));
}


__host__ __device__
__forceinline__ cuDoubleComplex
cmul_c8_1(const cuDoubleComplex c1,
          const cuDoubleComplex c2) {
       return (cuCmul(c1,c2));
}


__host__ __device__
__forceinline__ cuFloatComplex
csub_c4_1(const cuFloatComplex c1,
          const cuFloatComplex c2) {
      return (cuCsubf(c1,c2));
}


__host__ __device__
__forceinline__ cuDoubleComplex
csub_c8_1(const cuDoubleComplex c1,
          const cuDoubleComplex c2) {
      return (cuCsub(c1,c2));
}


__host__ __device__
__forceinline__ cuFloatComplex
cadd_c4_1(const cuFloatComplex c1,
          const cuFloatComplex c2) {
      return (cuCaddf(c1,c2));
}


__host__ __device__
__forceinline__ cuDoubleComplex
cadd_c8_1(const cuDoubleComplex c1,
          const cuDoubleComplex c2) {
      return (cuCadd(c1,c2));
}


__host__ __device__
__forceinline__ float
cmag_c4_1(const cuFloatComplex c) {
    const float re=cuCrealf(c);
    const float im=cuCimagf(c);
    return (re*re+im*im);
}


__host__ __device__
__forceinline__ double
cmag_c8_1(const cuDoubleComplex c) {
    const double re=cuCreal(c);
    const double im=cuCimag(c);
    return (re*re+im*im);
}















#endif /*__GMS_GPU_HELPERS_CUH__*/
