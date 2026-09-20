

#include <cstdint>
#include "GMS_atomicops.h"
//#include "snrm2_macros.h"
// Modified by Bernard Gingold, contact: bgin@gmail.com
/*
    GLAS is licensed under the The MIT License:

Copyright (c) 2011 Hans Henrik Brandenborg Sørensen, DTU

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

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

// 
// Generic Euclidian norm kernel.
//
// - Designed for architectures with atomicCAS, i.e., sm_11 and higher.
// - Simple design for all capabilities
// - Handled by splitting matrix as
// --------------------------
// |                 |      |
// |                 |      |
// |                 |      |
// | left tiled part | rest | 
// |                 |      |
// |                 |      |
// |                 |      |
// --------------------------
// - Kernel is generic and works for all cases
// - Unrolling of inner loops
//     1) to reduce instruction count from loop
//     2) to collect shared loads into 128 bit transactions (compiler should do that)
//   * Level determined by WORKSIZE_n
//   * For WORKSIZE_n = 1:
//       Given by a constant: full unroll
//       Given by a variable: unroll level 4
__device__ unsigned int count = 0; 
__shared__ bool isLastBlockDone;


inline __device__ void block_reduce_r4(float * __restrict__ y, 
                                    volatile float * __restrict__ smem, 
                                    float sum, 
                                    const int32_t tid,
                                    const int32_t BLOCKSIZE_n )
{
	smem[tid] = sum;
	__syncthreads();
	if (BLOCKSIZE_n >= 1024) { if (tid < 512) { smem[tid] += smem[tid + 512]; } __syncthreads(); }
	if (BLOCKSIZE_n >= 512) { if (tid < 256) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
	if (BLOCKSIZE_n == 480) { if (tid < 224) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
	if (BLOCKSIZE_n == 448) { if (tid < 192) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
	if (BLOCKSIZE_n == 416) { if (tid < 160) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
	if (BLOCKSIZE_n == 384) { if (tid < 128) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
	if (BLOCKSIZE_n >= 256) { if (tid < 128) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
	if (BLOCKSIZE_n == 352) { if (tid < 96) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
	if (BLOCKSIZE_n == 224) { if (tid < 96) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
	if (BLOCKSIZE_n == 320) { if (tid < 64) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
	if (BLOCKSIZE_n == 192) { if (tid < 64) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
	if (BLOCKSIZE_n >= 128) { if (tid < 64) { smem[tid] += smem[tid + 64]; } __syncthreads(); }
	if (BLOCKSIZE_n == 288) { if (tid < 32) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
	if (BLOCKSIZE_n == 160) { if (tid < 32) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
	if (BLOCKSIZE_n == 96) { if (tid < 32) { smem[tid] += smem[tid + 64]; } __syncthreads(); }
	if (BLOCKSIZE_n >= 64) { if (tid < 32) { smem[tid] += smem[tid + 32]; } }
	if (BLOCKSIZE_n >= 32) { if (tid < 16) { smem[tid] += smem[tid + 16]; } }
	if (BLOCKSIZE_n >= 16) { if (tid < 8) { smem[tid] += smem[tid + 8]; } }
	if (BLOCKSIZE_n >= 8) { if (tid < 4) { smem[tid] += smem[tid + 4]; } }
	if (BLOCKSIZE_n >= 4) { if (tid < 2) { smem[tid] += smem[tid + 2]; } }
	if (tid == 0) atomicFloatAdd(y, smem[tid] + smem[tid + 1]); 
}


inline __device__ void block_reduce_r8(double * __restrict__ y, 
                                       volatile double * __restrict__ smem, 
                                       double sum, 
                                       const int32_t tid,
                                       const int32_t BLOCKSIZE_n )
{
	smem[tid] = sum;
	__syncthreads();
	if (BLOCKSIZE_n >= 1024) { if (tid < 512) { smem[tid] += smem[tid + 512]; } __syncthreads(); }
	if (BLOCKSIZE_n >= 512) { if (tid < 256) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
	if (BLOCKSIZE_n == 480) { if (tid < 224) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
	if (BLOCKSIZE_n == 448) { if (tid < 192) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
	if (BLOCKSIZE_n == 416) { if (tid < 160) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
	if (BLOCKSIZE_n == 384) { if (tid < 128) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
	if (BLOCKSIZE_n >= 256) { if (tid < 128) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
	if (BLOCKSIZE_n == 352) { if (tid < 96) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
	if (BLOCKSIZE_n == 224) { if (tid < 96) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
	if (BLOCKSIZE_n == 320) { if (tid < 64) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
	if (BLOCKSIZE_n == 192) { if (tid < 64) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
	if (BLOCKSIZE_n >= 128) { if (tid < 64) { smem[tid] += smem[tid + 64]; } __syncthreads(); }
	if (BLOCKSIZE_n == 288) { if (tid < 32) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
	if (BLOCKSIZE_n == 160) { if (tid < 32) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
	if (BLOCKSIZE_n == 96) { if (tid < 32) { smem[tid] += smem[tid + 64]; } __syncthreads(); }
	if (BLOCKSIZE_n >= 64) { if (tid < 32) { smem[tid] += smem[tid + 32]; } }
	if (BLOCKSIZE_n >= 32) { if (tid < 16) { smem[tid] += smem[tid + 16]; } }
	if (BLOCKSIZE_n >= 16) { if (tid < 8) { smem[tid] += smem[tid + 8]; } }
	if (BLOCKSIZE_n >= 8) { if (tid < 4) { smem[tid] += smem[tid + 4]; } }
	if (BLOCKSIZE_n >= 4) { if (tid < 2) { smem[tid] += smem[tid + 2]; } }
	if (tid == 0) atomicFloatAdd(y, smem[tid] + smem[tid + 1]); 
}


__global__ void snrm2_kernel_r4(const float * __restrict__ x, 
                             float * __restrict__ y, 
                             const int32_t n, 
                             const int32_t ntile,
                             const int32_t BLOCKSIZE_n, 
                             const int32_t WORKSIZE_m, 
                             const int32_t WORKSIZE_n)
{
	const int tid = threadIdx.x;
	const int idx = blockIdx.x * BLOCKSIZE_n * WORKSIZE_n;
	x += idx + tid;
	float sum = 0.0f;
	if (idx != ntile)
	{
        #pragma unroll
          for (int i = 0; i < WORKSIZE_n; i++ )				
	   {												
		sum += x[0] * x[0]; x += BLOCKSIZE_n;		
	   }
		//inner_loop(WORKSIZE_n);
	}
	else
	{
        #pragma unroll 4
           for (int i = 0; i < n-idx-tid; i += BLOCKSIZE_n)
	    {
		sum += x[i] * x[i];
	   }
		//inner_loop2(n - idx - tid);
	}
	__shared__ volatile float smem[BLOCKSIZE_n];
	block_reduce_r4(y, smem, sum, tid, BLOCKSIZE_n);
	if (tid == 0) 
	{
		unsigned int32_t value = atomicInc(&count, gridDim.x); 
		isLastBlockDone = (value == (gridDim.x - 1)); 
	}
	__syncthreads(); 
	if (isLastBlockDone) {
		y[0] = sqrt(y[0]);
		count = 0;
	}
}


__global__ void snrm2_kernel_r8(const double * __restrict__ x, 
                                double * __restrict__ y, 
                                const int32_t n, 
                                const int32_t ntile,
                                const int32_t BLOCKSIZE_n, 
                                const int32_t WORKSIZE_m, 
                                const int32_t WORKSIZE_n)
{
	const int tid = threadIdx.x;
	const int idx = blockIdx.x * BLOCKSIZE_n * WORKSIZE_n;
	x += idx + tid;
	double sum = 0.0;
	if (idx != ntile)
	{
        #pragma unroll
          for (int i = 0; i < WORKSIZE_n; i++ )				
	   {												
		sum += x[0] * x[0]; x += BLOCKSIZE_n;		
	   }
		//inner_loop(WORKSIZE_n);
	}
	else
	{
        #pragma unroll 4
           for (int i = 0; i < n-idx-tid; i += BLOCKSIZE_n)
	    {
		sum += x[i] * x[i];
	   }
		//inner_loop2(n - idx - tid);
	}
	__shared__ volatile float smem[BLOCKSIZE_n];
	block_reduce_r8(y, smem, sum, tid, BLOCKSIZE_n);
	if (tid == 0) 
	{
		unsigned int32_t value = atomicInc(&count, gridDim.x); 
		isLastBlockDone = (value == (gridDim.x - 1)); 
	}
	__syncthreads(); 
	if (isLastBlockDone) {
		y[0] = sqrt(y[0]);
		count = 0;
	}
}




__global__ void snrm2_kernel2_r4(const float* __restrict__ x, 
                                 float* __restrict__ y, 
                                 const int32_t n, 
                                 const int32_t ntile,
                                 const int32_t BLOCKSIZE_n, 
                                 const int32_t WORKSIZE_m, 
                                 const int32_t WORKSIZE_n)
{
	const int tid = threadIdx.x;
	const int idx = blockIdx.x * BLOCKSIZE_n * WORKSIZE_n;
	x += idx + tid;
	float sum = 0.0f;
	if (idx != ntile)
	{
		if (WORKSIZE_m == 1)
		{	// NB: full unroll
            #pragma unroll
                    for (int i = 0; i < WORKSIZE_n; i++ )				
	                 {												
		            sum += x[0] * x[0]; x += BLOCKSIZE_n;		
	            }
			//inner_loop(WORKSIZE_n);
		}
		if (WORKSIZE_m == 2)
		{
            #pragma unroll 2
                   for (int i = 0; i < WORKSIZE_n; i++ )				
	                 {												
		            sum += x[0] * x[0]; x += BLOCKSIZE_n;		
	            }
			//inner_loop(WORKSIZE_n);
		}
		if (WORKSIZE_m == 3)
		{
            #pragma unroll 3
                   for (int i = 0; i < WORKSIZE_n; i++ )				
	                 {												
		            sum += x[0] * x[0]; x += BLOCKSIZE_n;		
	            }
			//inner_loop(WORKSIZE_n);
		}
		if (WORKSIZE_m == 4)
		{
            #pragma unroll 4
                    for (int i = 0; i < WORKSIZE_n; i++ )				
	                 {												
		            sum += x[0] * x[0]; x += BLOCKSIZE_n;		
	            }
			//inner_loop(WORKSIZE_n);
		}
		if (WORKSIZE_m == 5)
		{
            #pragma unroll 5 
                    for (int i = 0; i < WORKSIZE_n; i++ )				
	                 {												
		            sum += x[0] * x[0]; x += BLOCKSIZE_n;		
	            }
			//inner_loop(WORKSIZE_n);
		}
		if (WORKSIZE_m == 6)
		{
            #pragma unroll 6
                    for (int i = 0; i < WORKSIZE_n; i++ )				
	                 {												
		            sum += x[0] * x[0]; x += BLOCKSIZE_n;		
	            }
			//inner_loop(WORKSIZE_n);
		}
		if (WORKSIZE_m == 7)
		{
            #pragma unroll 7
                     for (int i = 0; i < WORKSIZE_n; i++ )				
	                 {												
		            sum += x[0] * x[0]; x += BLOCKSIZE_n;		
	            }
			//inner_loop(WORKSIZE_n);
		}
		if (WORKSIZE_m == 8)
		{
            #pragma unroll 8
                     for (int i = 0; i < WORKSIZE_n; i++ )				
	                 {												
		            sum += x[0] * x[0]; x += BLOCKSIZE_n;		
	            }
			///inner_loop(WORKSIZE_n);
		}
	}
	else
	{
		if (WORKSIZE_m == 1)
		{
            #pragma unroll 4
                    for (int i = 0; i < (n-idx-tid); i += BLOCKSIZE_n) {
	
		        sum += x[i] * x[i];
	            }
			//inner_loop2(n - idx - tid);
		}
		if (WORKSIZE_m == 2)
		{
            #pragma unroll 2
                    for (int i = 0; i < (n-idx-tid); i += BLOCKSIZE_n) {
	
		        sum += x[i] * x[i];
	            }
			//inner_loop2(n - idx - tid);
		}
		if (WORKSIZE_m == 3)
		{
            #pragma unroll 3
                    for (int i = 0; i < (n-idx-tid); i += BLOCKSIZE_n) {
	
		        sum += x[i] * x[i];
	            }
			//inner_loop2(n - idx - tid);
		}
		if (WORKSIZE_m == 4)
		{
            #pragma unroll 4
                    for (int i = 0; i < (n-idx-tid); i += BLOCKSIZE_n) {
	
		        sum += x[i] * x[i];
	            }
			//inner_loop2(n - idx - tid);
		}
		if (WORKSIZE_m == 5)
		{
            #pragma unroll 5 
                     for (int i = 0; i < (n-idx-tid); i += BLOCKSIZE_n) {
	
		        sum += x[i] * x[i];
	            }
			//inner_loop2(n - idx - tid);
		}
		if (WORKSIZE_m == 6)
		{
            #pragma unroll 6
                      for (int i = 0; i < (n-idx-tid); i += BLOCKSIZE_n) {
	
		        sum += x[i] * x[i];
	            }
			//inner_loop2(n - idx - tid);
		}
		if (WORKSIZE_m == 7)
		{
            #pragma unroll 7
                      for (int i = 0; i < (n-idx-tid); i += BLOCKSIZE_n) {
	
		        sum += x[i] * x[i];
	            }
			//inner_loop2(n - idx - tid);
		}
		if (WORKSIZE_m == 8)
		{
            #pragma unroll 8
                     for (int i = 0; i < (n-idx-tid); i += BLOCKSIZE_n) {
	
		        sum += x[i] * x[i];
	            }
			//inner_loop2(n - idx - tid);
		}
	}
	__shared__ volatile float smem[BLOCKSIZE_n];
	block_reduce_r4(y, smem, sum, tid, BLOCKSIZE_n);
	if (tid == 0) 
	{
		unsigned int value = atomicInc(&count, gridDim.x); 
		isLastBlockDone = (value == (gridDim.x - 1)); 
	}
	__syncthreads(); 
	if (isLastBlockDone) {
		y[0] = sqrt(y[0]);
		count = 0;
	}
}


__global__ void snrm2_kernel2_r8(const double* __restrict__ x, 
                                 double* __restrict__ y, 
                                 const int32_t n, 
                                 const int32_t ntile,
                                 const int32_t BLOCKSIZE_n, 
                                 const int32_t WORKSIZE_m, 
                                 const int32_t WORKSIZE_n)
{
	const int tid = threadIdx.x;
	const int idx = blockIdx.x * BLOCKSIZE_n * WORKSIZE_n;
	x += idx + tid;
	double sum = 0.0;
	if (idx != ntile)
	{
		if (WORKSIZE_m == 1)
		{	// NB: full unroll
            #pragma unroll
                    for (int i = 0; i < WORKSIZE_n; i++ )				
	                 {												
		            sum += x[0] * x[0]; x += BLOCKSIZE_n;		
	            }
			//inner_loop(WORKSIZE_n);
		}
		if (WORKSIZE_m == 2)
		{
            #pragma unroll 2
                   for (int i = 0; i < WORKSIZE_n; i++ )				
	                 {												
		            sum += x[0] * x[0]; x += BLOCKSIZE_n;		
	            }
			//inner_loop(WORKSIZE_n);
		}
		if (WORKSIZE_m == 3)
		{
            #pragma unroll 3
                   for (int i = 0; i < WORKSIZE_n; i++ )				
	                 {												
		            sum += x[0] * x[0]; x += BLOCKSIZE_n;		
	            }
			//inner_loop(WORKSIZE_n);
		}
		if (WORKSIZE_m == 4)
		{
            #pragma unroll 4
                    for (int i = 0; i < WORKSIZE_n; i++ )				
	                 {												
		            sum += x[0] * x[0]; x += BLOCKSIZE_n;		
	            }
			//inner_loop(WORKSIZE_n);
		}
		if (WORKSIZE_m == 5)
		{
            #pragma unroll 5 
                    for (int i = 0; i < WORKSIZE_n; i++ )				
	                 {												
		            sum += x[0] * x[0]; x += BLOCKSIZE_n;		
	            }
			//inner_loop(WORKSIZE_n);
		}
		if (WORKSIZE_m == 6)
		{
            #pragma unroll 6
                    for (int i = 0; i < WORKSIZE_n; i++ )				
	                 {												
		            sum += x[0] * x[0]; x += BLOCKSIZE_n;		
	            }
			//inner_loop(WORKSIZE_n);
		}
		if (WORKSIZE_m == 7)
		{
            #pragma unroll 7
                     for (int i = 0; i < WORKSIZE_n; i++ )				
	                 {												
		            sum += x[0] * x[0]; x += BLOCKSIZE_n;		
	            }
			//inner_loop(WORKSIZE_n);
		}
		if (WORKSIZE_m == 8)
		{
            #pragma unroll 8
                     for (int i = 0; i < WORKSIZE_n; i++ )				
	                 {												
		            sum += x[0] * x[0]; x += BLOCKSIZE_n;		
	            }
			///inner_loop(WORKSIZE_n);
		}
	}
	else
	{
		if (WORKSIZE_m == 1)
		{
            #pragma unroll 4
                    for (int i = 0; i < (n-idx-tid); i += BLOCKSIZE_n) {
	
		        sum += x[i] * x[i];
	            }
			//inner_loop2(n - idx - tid);
		}
		if (WORKSIZE_m == 2)
		{
            #pragma unroll 2
                    for (int i = 0; i < (n-idx-tid); i += BLOCKSIZE_n) {
	
		        sum += x[i] * x[i];
	            }
			//inner_loop2(n - idx - tid);
		}
		if (WORKSIZE_m == 3)
		{
            #pragma unroll 3
                    for (int i = 0; i < (n-idx-tid); i += BLOCKSIZE_n) {
	
		        sum += x[i] * x[i];
	            }
			//inner_loop2(n - idx - tid);
		}
		if (WORKSIZE_m == 4)
		{
            #pragma unroll 4
                    for (int i = 0; i < (n-idx-tid); i += BLOCKSIZE_n) {
	
		        sum += x[i] * x[i];
	            }
			//inner_loop2(n - idx - tid);
		}
		if (WORKSIZE_m == 5)
		{
            #pragma unroll 5 
                     for (int i = 0; i < (n-idx-tid); i += BLOCKSIZE_n) {
	
		        sum += x[i] * x[i];
	            }
			//inner_loop2(n - idx - tid);
		}
		if (WORKSIZE_m == 6)
		{
            #pragma unroll 6
                      for (int i = 0; i < (n-idx-tid); i += BLOCKSIZE_n) {
	
		        sum += x[i] * x[i];
	            }
			//inner_loop2(n - idx - tid);
		}
		if (WORKSIZE_m == 7)
		{
            #pragma unroll 7
                      for (int i = 0; i < (n-idx-tid); i += BLOCKSIZE_n) {
	
		        sum += x[i] * x[i];
	            }
			//inner_loop2(n - idx - tid);
		}
		if (WORKSIZE_m == 8)
		{
            #pragma unroll 8
                     for (int i = 0; i < (n-idx-tid); i += BLOCKSIZE_n) {
	
		        sum += x[i] * x[i];
	            }
			//inner_loop2(n - idx - tid);
		}
	}
	__shared__ volatile double smem[BLOCKSIZE_n];
	block_reduce_r8(y, smem, sum, tid, BLOCKSIZE_n);
	if (tid == 0) 
	{
		unsigned int value = atomicInc(&count, gridDim.x); 
		isLastBlockDone = (value == (gridDim.x - 1)); 
	}
	__syncthreads(); 
	if (isLastBlockDone) {
		y[0] = sqrt(y[0]);
		count = 0;
	}
}


