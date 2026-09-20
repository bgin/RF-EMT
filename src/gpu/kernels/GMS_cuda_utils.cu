
#ifndef __GMS_CUDA_UTILS_CU__
#define __GMS_CUDA_UTILS_CU__

/*

     A collection of various short cuda utility
     functions both device and kernels.
     Part of them adapted from the SO answers.
*/


#include <cub/cub.cuh>
#include <cuda.h>
#include <cstdint>


//
// Block reduction kernel
// 

__global__
void
block_reduction_r4(const float * __restrict__ in,
                   float       * __restrict__ out,
                   const int32_t n) {
     
       int tid = blockIdx.x*blockDIm.x+threadIdx.x;
       typedef cub::BlockReduce<float,32> br; // block-size == 32
       __shared__ br::TempStorage ts;
       float sum = 0.0f;
       if(tid < n) sum = br(ts).Sum(in[tid]);
       if(threadIdx.x == 0) out[blockIdx.x] = sum;
}
//===========================================================================//
// Adapted from the Cuda-samples and modified by removal of templated code.
//===========================================================================//
/* Copyright (c) 2022, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>

// Utility class used to avoid the linker errors.
// Templating removed created two three versions for the following data types
// 1) float
// 2) int
// 3) cuComplex

struct SharedMemoryR4 {
 
   __device__ inline operator float *() {
         extern __shared__ float __smem[];
         return (float*)__smem;
  }

  __device__ inline operator const float *() const {
         extern __shared__ float __smem[];
         return (float*)__smem;
  }

};


struct SharedMemoryI4 {

   __device__ inline operator int *() {
        extern __shared__ int __smem[];
        return (int*)__smem;
  }

  __device__ inline operator const int *() const {
        extern __shared__ int __smem[];
        return (int*)__smem;
  }
};


struct SharedMemoryC4 {

    __device__ inline operator cuComplex *() {
            extern __shared__ cuComplex __smem[];
            return (cuComplex*)__smem;
    }

    __device__ inline operator const cuComplex *() const {
            extern __shared__ cuComplex __smem[];
            return (cuComplex*)__smem;
    }
};


__device__ __forceinline__
float warpReduceSum(uint32_t mask, 
                    float sum) {
  for(int32_t off = warpSize/2; off > 0; off /= 2) {
      sum += __shfl_down_sync(mask,sum,off);
  }
   return (sum);
}

  /*
    This version uses sequential addressing -- no divergence or bank conflicts.
*/

__global__
void reduce_r4_v1(float * __restrict__ in,
                  float * __restrict__ out,
                  const uint32_t n) {
 
   // Create thread block group
   namespace cg = cooperative_groups;
   cg::thread_block tb = cg::this_thread_block();
   float * __restrict__ shmem  = SharedMemoryR4();
   
   // Shmem loading
   uint32_t tid = threadIdx.x;
   uint32_t i   = blockIdx.x*blockDIm.x+threadIdx.x;
   shmem[tid]   = (i<n) ? in[i] : 0;
   cg::sync(tb);
   
   // Reduction in shared memory
   for(uint32_t s = blockDIm.x/2; s > 0; s >>= 1) {
       if(tid<s){
          shmem[tid] += shmem[tid+s];
       }
    cg::sync(tb);
   }
   // SToring result in global memory
   if(tid == 0) out[blockIdx.x] = shmem[0];
}

__global__
void reduce_i4_v1(int32_t * __restrict__ in,
                  int32_t * __restrict__ out,
                  const uint32_t n) {
 
   // Create thread block group
   namespace cg = cooperative_groups;
   cg::thread_block tb = cg::this_thread_block();
   int32_t * __restrict__ shmem  = SharedMemoryI4();
   
   // Shmem loading
   uint32_t tid = threadIdx.x;
   uint32_t i   = blockIdx.x*blockDIm.x+threadIdx.x;
   shmem[tid]   = (i<n) ? in[i] : 0;
   cg::sync(tb);
   
   // Reduction in shared memory
   for(uint32_t s = blockDIm.x/2; s > 0; s >>= 1) {
       if(tid<s){
          shmem[tid] += shmem[tid+s];
       }
    cg::sync(tb);
   }
   // SToring result in global memory
   if(tid == 0) out[blockIdx.x] = shmem[0];
}

__global__
void reduce_c4_v1(cuComplex * __restrict__ in,
                  cuComplex * __restrict__ out,
                  const uint32_t n) {
 
   // Create thread block group
   namespace cg = cooperative_groups;
   cg::thread_block tb = cg::this_thread_block();
   cuComplex * __restrict__ shmem  = SharedMemoryC4();
   
   // Shmem loading
   uint32_t tid = threadIdx.x;
   uint32_t i   = blockIdx.x*blockDIm.x+threadIdx.x;
   shmem[tid]   = (i<n) ? in[i] : 0;
   cg::sync(tb);
   
   // Reduction in shared memory
   for(uint32_t s = blockDIm.x/2; s > 0; s >>= 1) {
       if(tid<s){
          shmem[tid] += shmem[tid+s];
       }
    cg::sync(tb);
   }
   // SToring result in global memory
   if(tid == 0) out[blockIdx.x] = shmem[0];
}


/*
   This version uses n/2 threads --
    it performs the first level of reduction when reading from global memory. 
*/

__global__
void reduce_r4_v2(float * __restrict__ in,
                  float * __restrict__ out,
                  const uint32_t n) {
     
     // Thread block group instantiation
     namespace cg = cooperative_groups;
     cg::thread_block tb = cg::this_thread_block();
     float * __restrict__ shmem = SharedMemoryR4();
     
     // first level of reduction
     // global memory read, write to shared memory
     uint32_t tid = threadIdx.x;
     uint32_t i   = blockIdx.x*(blockDimx.x/2)+threadIdx.x;
     float sum    = (i<n) ? in[i] : 0.0f;
     if(i+blockDim.x<n) sum += in[i+blockDim.x];
     shmem[tid] = sum;
     cg::sync(tb); // barier before an entry to shared memory reduce ops
     for(uint32_t s = blockDim.x/2; s > 0; s >>= 1) {
         if(tid<s) {
             const float tmp = sum
             shmem[tid] = tmp + shmem[tid+s];
         }
     }
     cg::synx(tb);
     // Result is written back to global memory
     if(tid==0) out[blockIdx.x] = sum;
}


__global__
void reduce_i4_v2(int32_t * __restrict__ in,
                  int32_t * __restrict__ out,
                  const uint32_t n) {
     
     // Thread block group instantiation
     namespace cg = cooperative_groups;
     cg::thread_block tb = cg::this_thread_block();
     int32_t * __restrict__ shmem = SharedMemoryI4();
     
     // first level of reduction
     // global memory read, write to shared memory
     uint32_t tid = threadIdx.x;
     uint32_t i   = blockIdx.x*(blockDimx.x/2)+threadIdx.x;
     int32_t sum    = (i<n) ? in[i] : 0;
     if(i+blockDim.x<n) sum += in[i+blockDim.x];
     shmem[tid] = sum;
     cg::sync(tb); // barier before an entry to shared memory reduce ops
     for(uint32_t s = blockDim.x/2; s > 0; s >>= 1) {
         if(tid<s) {
             const int32_t tmp = sum
             shmem[tid] = tmp + shmem[tid+s];
         }
     }
     cg::synx(tb);
     // Result is written back to global memory
     if(tid==0) out[blockIdx.x] = sum;
}


__global__
void reduce_c4_v2(cuComplex * __restrict__ in,
                  cuComplex * __restrict__ out,
                  const uint32_t n) {
     
     // Thread block group instantiation
     namespace cg = cooperative_groups;
     cg::thread_block tb = cg::this_thread_block();
     cuComplex * __restrict__ shmem = SharedMemoryC4();
     
     // first level of reduction
     // global memory read, write to shared memory
     uint32_t tid = threadIdx.x;
     uint32_t i   = blockIdx.x*(blockDimx.x/2)+threadIdx.x;
     cuComplex sum    = (i<n) ? in[i] : 0;
     if(i+blockDim.x<n) sum += in[i+blockDim.x];
     shmem[tid] = sum;
     cg::sync(tb); // barier before an entry to shared memory reduce ops
     for(uint32_t s = blockDim.x/2; s > 0; s >>= 1) {
         if(tid<s) {
             const cuComplex tmp = sum
             shmem[tid] = tmp + shmem[tid+s];
         }
     }
     cg::synx(tb);
     // Result is written back to global memory
     if(tid==0) out[blockIdx.x] = sum;
}


/*
    This version uses the warp shuffle operation if available to reduce
    warp synchronization. When shuffle is not available the final warp's
    worth of work is unrolled to reduce looping overhead.

    Note, this kernel needs a minimum of 64*sizeof(T) bytes of shared memory.
    In other words if blockSize <= 32, allocate 64*sizeof(T) bytes.
    If blockSize > 32, allocate blockSize*sizeof(T) bytes.
*/

__global__
void reduce_r4_v3(float * __restrict__ in,
                  float * __restrict__ out,
                  const uint32_t n,
                  const uint32_t blocksize) {
  
      // Thread block group instantiation
     namespace cg = cooperative_groups;
     cg::thread_block tb = cg::this_thread_block();
     float * __restrict__ shmem = SharedMemoryR4();
      // first level of reduction
     // global memory read, write to shared memory
     uint32_t tid = threadIdx.x;
     uint32_t i   = blockIdx.x*(blockDim.x/2)+threadIdx.x;
     float sum    = (i<n) ? in[i] : 0.0f;
     if(i+blocksize<n) sum += in[i+blocksize];
     shmem[tid] = sum;
     cg::sync(tb);
     // Reduction in shared memory
     for(uint32_t s = blockDim.x/2; s > 0; s >>= 1) {
         if(tid<s) {
             const float tmp = sum
             shmem[tid] = tmp + shmem[tid+s];
         }
     }
     cg::sync(tb);
     
     cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(tb);
     if(tb.thread_rank() < 32) {
        // Fetch final intermediate sum from 2nd warp
        if(blocksize >= 64) sum += shmem[tid+32];
        // Final warp reduction (shuffle)
        for(int32_t off = tile32.size()/2; off > 0; off /= 2) {
             sum += tile32.shfl_down(sum,off);
        }
     }
     // Result is written to global memory
     if(tb.thread_rank() == 0) out[blockIdx.x] = sum;
}


__global__
void reduce_i4_v3(int32_t * __restrict__ in,
                  int32_t * __restrict__ out,
                  const uint32_t n,
                  const uint32_t blocksize) {
  
      // Thread block group instantiation
     namespace cg = cooperative_groups;
     cg::thread_block tb = cg::this_thread_block();
     int32_t * __restrict__ shmem = SharedMemoryI4();
      // first level of reduction
     // global memory read, write to shared memory
     uint32_t tid = threadIdx.x;
     uint32_t i   = blockIdx.x*(blockDim.x/2)+threadIdx.x;
     int32_t sum    = (i<n) ? in[i] : 0;
     if(i+blocksize<n) sum += in[i+blocksize];
     shmem[tid] = sum;
     cg::sync(tb);
     // Reduction in shared memory
     for(uint32_t s = blockDim.x/2; s > 0; s >>= 1) {
         if(tid<s) {
             const int32_t tmp = sum
             shmem[tid] = tmp + shmem[tid+s];
         }
     }
     cg::sync(tb);
     
     cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(tb);
     if(tb.thread_rank() < 32) {
        // Fetch final intermediate sum from 2nd warp
        if(blocksize >= 64) sum += shmem[tid+32];
        // Final warp reduction (shuffle)
        for(int32_t off = tile32.size()/2; off > 0; off /= 2) {
             sum += tile32.shfl_down(sum,off);
        }
     }
     // Result is written to global memory
     if(tb.thread_rank() == 0) out[blockIdx.x] = sum;
}


__global__
void reduce_c4_v3(cuComplex * __restrict__ in,
                  cuComplex * __restrict__ out,
                  const uint32_t n,
                  const uint32_t blocksize) {
  
      // Thread block group instantiation
     namespace cg = cooperative_groups;
     cg::thread_block tb = cg::this_thread_block();
     cuComplex * __restrict__ shmem = SharedMemoryC4();
      // first level of reduction
     // global memory read, write to shared memory
     uint32_t tid = threadIdx.x;
     uint32_t i   = blockIdx.x*(blockDim.x/2)+threadIdx.x;
     cuComplex sum    = (i<n) ? in[i] : {0.0f,0.0f};
     if(i+blocksize<n) sum += in[i+blocksize];
     shmem[tid] = sum;
     cg::sync(tb);
     // Reduction in shared memory
     for(uint32_t s = blockDim.x/2; s > 0; s >>= 1) {
         if(tid<s) {
             const cuComplex tmp = sum
             shmem[tid] = tmp + shmem[tid+s];
         }
     }
     cg::sync(tb);
     
     cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(tb);
     if(tb.thread_rank() < 32) {
        // Fetch final intermediate sum from 2nd warp
        if(blocksize >= 64) sum += shmem[tid+32];
        // Final warp reduction (shuffle)
        for(int32_t off = tile32.size()/2; off > 0; off /= 2) {
             sum += tile32.shfl_down(sum,off);
        }
     }
     // Result is written to global memory
     if(tb.thread_rank() == 0) out[blockIdx.x] = sum;
}

/*
    This version is completely unrolled, unless warp shuffle is available, then
    shuffle is used within a loop.  It uses a template parameter to achieve
    optimal code for any (power of 2) number of threads.  This requires a switch
    statement in the host code to handle all the different thread block sizes at
    compile time. When shuffle is available, it is used to reduce warp
   synchronization.
    Note, this kernel needs a minimum of 64*sizeof(T) bytes of shared memory.
    In other words if blockSize <= 32, allocate 64*sizeof(T) bytes.
    If blockSize > 32, allocate blockSize*sizeof(T) bytes.
*/

__global__
void reduce_r4_v4(float * __restrict__ in,
                  float * __restrict__ out,
                  const uint32_t n,
                  const uint32_t blocksize) {

     // Thread block group instantiation
     namespace cg = cooperative_groups;
     cg::thread_block tb = cg::this_thread_block();
     float * __restrict__ shmem = SharedMemoryR4();
      // first level of reduction
     // global memory read, write to shared memory
     uint32_t tid = threadIdx.x;
     uint32_t i   = blockIdx.x*(blocksize*2)+threadIdx.x;
     float sum = (i<n) ? in[i] : 0.0f;
     if(i+blocksize<n) sum += in[i+blocksize];
     shmem[tid] = sum; // partial sums
     cg::sync(tb);
     // reduction in shared memory
     if ((blocksize>=512) && (tid<256)) {
         const float tmp = sum
         shmem[tid] = tmp + shmem[tid+256];
      }
      cg::sync(tb);
     if ((blocksize>=256) && (tid<128)) {
         const float tmp = sum
         shmem[tid] = tmp + shmem[tid+128];
     }
     cg::sync(tb);
     if ((blocksize>=128) && (tid<64)) {
         const float tmp = sum
         shmem[tid] = tmp + shmem[tid+64];
     }
     cg::sync(tb);

     cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(tb);
     if(tb.thread_rank()<32) {
        if(blocksize>=64) sum += shmem[tid+32];
        for(int32_t off = tile32.size()/2; off > 0; offset /= 2) {
            sum += tile32.shfl_down(sum,off);
         }
     }
     if(tb.thread_rank()==0) out[blockIdx.x] = sum;
}


__global__
void reduce_i4_v4(int32_t * __restrict__ in,
                  int32_t * __restrict__ out,
                  const uint32_t n,
                  const uint32_t blocksize) {

     // Thread block group instantiation
     namespace cg = cooperative_groups;
     cg::thread_block tb = cg::this_thread_block();
     int32_t * __restrict__ shmem = SharedMemoryI4();
      // first level of reduction
     // global memory read, write to shared memory
     uint32_t tid = threadIdx.x;
     uint32_t i   = blockIdx.x*(blocksize*2)+threadIdx.x;
     int32_t sum = (i<n) ? in[i] : 0;
     if(i+blocksize<n) sum += in[i+blocksize];
     shmem[tid] = sum; // partial sums
     cg::sync(tb);
     // reduction in shared memory
     if ((blocksize>=512) && (tid<256)) {
         const int32_t tmp = sum
         shmem[tid] = tmp + shmem[tid+256];
      }
      cg::sync(tb);
     if ((blocksize>=256) && (tid<128)) {
         const int32_t tmp = sum
         shmem[tid] = tmp + shmem[tid+128];
     }
     cg::sync(tb);
     if ((blocksize>=128) && (tid<64)) {
         const int32_t tmp = sum
         shmem[tid] = tmp + shmem[tid+64];
     }
     cg::sync(tb);

     cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(tb);
     if(tb.thread_rank()<32) {
        if(blocksize>=64) sum += shmem[tid+32];
        for(int32_t off = tile32.size()/2; off > 0; offset /= 2) {
            sum += tile32.shfl_down(sum,off);
         }
     }
     if(tb.thread_rank()==0) out[blockIdx.x] = sum;
}


__global__
void reduce_c4_v4(cuComplex * __restrict__ in,
                  cuComplex * __restrict__ out,
                  const uint32_t n,
                  const uint32_t blocksize) {

     // Thread block group instantiation
     namespace cg = cooperative_groups;
     cg::thread_block tb = cg::this_thread_block();
     cuComplex * __restrict__ shmem = SharedMemoryC4();
      // first level of reduction
     // global memory read, write to shared memory
     uint32_t tid = threadIdx.x;
     uint32_t i   = blockIdx.x*(blocksize*2)+threadIdx.x;
     cuComplex sum = (i<n) ? in[i] : {0.0f,0.0f};
     if(i+blocksize<n) sum += in[i+blocksize];
     shmem[tid] = sum; // partial sums
     cg::sync(tb);
     // reduction in shared memory
     if ((blocksize>=512) && (tid<256)) {
         const cuComplex tmp = sum
         shmem[tid] = tmp + shmem[tid+256];
      }
      cg::sync(tb);
     if ((blocksize>=256) && (tid<128)) {
         const cuComplex tmp = sum
         shmem[tid] = tmp + shmem[tid+128];
     }
     cg::sync(tb);
     if ((blocksize>=128) && (tid<64)) {
         const cuComplex tmp = sum
         shmem[tid] = tmp + shmem[tid+64];
     }
     cg::sync(tb);

     cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(tb);
     if(tb.thread_rank()<32) {
        if(blocksize>=64) sum += shmem[tid+32];
        for(int32_t off = tile32.size()/2; off > 0; offset /= 2) {
            sum += tile32.shfl_down(sum,off);
         }
     }
     if(tb.thread_rank()==0) out[blockIdx.x] = sum;
}


/*
    This version adds multiple elements per thread sequentially.  This reduces
   the overall cost of the algorithm while keeping the work complexity O(n) and
   the step complexity O(log n). (Brent's Theorem optimization)
    Note, this kernel needs a minimum of 64*sizeof(T) bytes of shared memory.
    In other words if blockSize <= 32, allocate 64*sizeof(T) bytes.
    If blockSize > 32, allocate blockSize*sizeof(T) bytes.
*/


__global__ void reduce_r4_v5(float * __restrict__ in,
                            float * __restrict__ out,
                            uint32_t n,
                            uint32_t blocksize,
                            bool npow2 ) {
  // Handle to thread block group
  cg::thread_block cta = cg::this_thread_block();
  float * __restrict__ shmem = SharedMemoryR4();
  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  uint32_t tid = threadIdx.x;
  uint32_t gridSize = blockSize * gridDim.x;

  float sum = 0.0f;

  // we reduce multiple elements per thread.  The number is determined by the
  // number of active thread blocks (via gridDim).  More blocks will result
  // in a larger gridSize and therefore fewer elements per thread
  if (npow2) {
      uint32_t i = blockIdx.x * blockSize * 2 + threadIdx.x;
      gridSize = gridSize << 1;

    while (i < n) {
      sum += in[i];
      // ensure we don't read out of bounds -- this is optimized away for
      // powerOf2 sized arrays
      if ((i + blockSize) < n) {
           sum += in[i + blockSize];
      }
      i += gridSize;
    }
  } else {
    unsigned int i = blockIdx.x * blockSize + threadIdx.x;
    while (i < n) {
      sum += in[i];
      i += gridSize;
    }
  }

  // each thread puts its local sum into shared memory
  shmem[tid] = sum;
  cg::sync(cta);

  // do reduction in shared mem
  if ((blockSize >= 512) && (tid < 256)) {
      shmem[tid] = sum = sum + sdata[tid + 256];
  }

  cg::sync(cta);

  if ((blockSize >= 256) && (tid < 128)) {
      shmem[tid] = sum = sum + sdata[tid + 128];
  }

  cg::sync(cta);

  if ((blockSize >= 128) && (tid < 64)) {
     shmem[tid] = sum = sum + sdata[tid + 64];
  }

  cg::sync(cta);

  cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cta);

  if (cta.thread_rank() < 32) {
    // Fetch final intermediate sum from 2nd warp
    if (blocksize >= 64) sum += shmem[tid + 32];
    // Reduce final warp using shuffle
    for (int32_t off = tile32.size() / 2; off > 0; off /= 2) {
         sum += tile32.shfl_down(sum, off);
    }
  }

  // write result for this block to global mem
  if (cta.thread_rank() == 0) out[blockIdx.x] = sum;
}


__global__ void reduce_i4_v5(int32_t * __restrict__ in,
                             int32_t * __restrict__ out,
                             uint32_t n,
                             uint32_t blocksize,
                             bool npow2 ) {
  // Handle to thread block group
  cg::thread_block cta = cg::this_thread_block();
  int32_t * __restrict__ shmem = SharedMemoryI4();
  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  uint32_t tid = threadIdx.x;
  uint32_t gridSize = blockSize * gridDim.x;

  int32_t sum = 0;

  // we reduce multiple elements per thread.  The number is determined by the
  // number of active thread blocks (via gridDim).  More blocks will result
  // in a larger gridSize and therefore fewer elements per thread
  if (npow2) {
      uint32_t i = blockIdx.x * blockSize * 2 + threadIdx.x;
      gridSize = gridSize << 1;

    while (i < n) {
      sum += in[i];
      // ensure we don't read out of bounds -- this is optimized away for
      // powerOf2 sized arrays
      if ((i + blockSize) < n) {
           sum += in[i + blockSize];
      }
      i += gridSize;
    }
  } else {
    unsigned int i = blockIdx.x * blockSize + threadIdx.x;
    while (i < n) {
      sum += in[i];
      i += gridSize;
    }
  }

  // each thread puts its local sum into shared memory
  shmem[tid] = sum;
  cg::sync(cta);

  // do reduction in shared mem
  if ((blockSize >= 512) && (tid < 256)) {
      shmem[tid] = sum = sum + sdata[tid + 256];
  }

  cg::sync(cta);

  if ((blockSize >= 256) && (tid < 128)) {
      shmem[tid] = sum = sum + sdata[tid + 128];
  }

  cg::sync(cta);

  if ((blockSize >= 128) && (tid < 64)) {
     shmem[tid] = sum = sum + sdata[tid + 64];
  }

  cg::sync(cta);

  cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cta);

  if (cta.thread_rank() < 32) {
    // Fetch final intermediate sum from 2nd warp
    if (blocksize >= 64) sum += shmem[tid + 32];
    // Reduce final warp using shuffle
    for (int32_t off = tile32.size() / 2; off > 0; off /= 2) {
         sum += tile32.shfl_down(sum, off);
    }
  }

  // write result for this block to global mem
  if (cta.thread_rank() == 0) out[blockIdx.x] = sum;
}


__global__ void reduce_c4_v5(cuComplex * __restrict__ in,
                             cuComplex * __restrict__ out,
                             uint32_t n,
                             uint32_t blocksize,
                             bool npow2 ) {
  // Handle to thread block group
  cg::thread_block cta = cg::this_thread_block();
  cuComplex* __restrict__ shmem = SharedMemoryC4();
  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  uint32_t tid = threadIdx.x;
  uint32_t gridSize = blockSize * gridDim.x;

  cuComplex sum = {0.0f,0.0f};

  // we reduce multiple elements per thread.  The number is determined by the
  // number of active thread blocks (via gridDim).  More blocks will result
  // in a larger gridSize and therefore fewer elements per thread
  if (npow2) {
      uint32_t i = blockIdx.x * blockSize * 2 + threadIdx.x;
      gridSize = gridSize << 1;

    while (i < n) {
      sum += in[i];
      // ensure we don't read out of bounds -- this is optimized away for
      // powerOf2 sized arrays
      if ((i + blockSize) < n) {
           sum += in[i + blockSize];
      }
      i += gridSize;
    }
  } else {
    unsigned int i = blockIdx.x * blockSize + threadIdx.x;
    while (i < n) {
      sum += in[i];
      i += gridSize;
    }
  }

  // each thread puts its local sum into shared memory
  shmem[tid] = sum;
  cg::sync(cta);

  // do reduction in shared mem
  if ((blockSize >= 512) && (tid < 256)) {
      shmem[tid] = sum = sum + sdata[tid + 256];
  }

  cg::sync(cta);

  if ((blockSize >= 256) && (tid < 128)) {
      shmem[tid] = sum = sum + sdata[tid + 128];
  }

  cg::sync(cta);

  if ((blockSize >= 128) && (tid < 64)) {
     shmem[tid] = sum = sum + sdata[tid + 64];
  }

  cg::sync(cta);

  cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cta);

  if (cta.thread_rank() < 32) {
    // Fetch final intermediate sum from 2nd warp
    if (blocksize >= 64) sum += shmem[tid + 32];
    // Reduce final warp using shuffle
    for (int32_t off = tile32.size() / 2; off > 0; off /= 2) {
         sum += tile32.shfl_down(sum, off);
    }
  }

  // write result for this block to global mem
  if (cta.thread_rank() == 0) out[blockIdx.x] = sum;
}

/*

   template <class T>
void reduce(int size, int threads, int blocks, int whichKernel, T *d_idata,
            T *d_odata) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize =
      (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);

  // as kernel 9 - multi_warp_cg_reduce cannot work for more than 64 threads
  // we choose to set kernel 7 for this purpose.
  if (threads < 64 && whichKernel == 9)
  {
    whichKernel = 7;
  }

  // choose which of the optimized versions of reduction to launch
  switch (whichKernel) {
    case 0:
      reduce0<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
      break;

    case 1:
      reduce1<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
      break;

    case 2:
      reduce2<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
      break;

    case 3:
      reduce3<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
      break;

    case 4:
      switch (threads) {
        case 512:
          reduce4<T, 512>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;

        case 256:
          reduce4<T, 256>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;

        case 128:
          reduce4<T, 128>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;

        case 64:
          reduce4<T, 64>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;

        case 32:
          reduce4<T, 32>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;

        case 16:
          reduce4<T, 16>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;

        case 8:
          reduce4<T, 8>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;

        case 4:
          reduce4<T, 4>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;

        case 2:
          reduce4<T, 2>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;

        case 1:
          reduce4<T, 1>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;
      }

      break;

    case 5:
      switch (threads) {
        case 512:
          reduce5<T, 512>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;

        case 256:
          reduce5<T, 256>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;

        case 128:
          reduce5<T, 128>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;

        case 64:
          reduce5<T, 64>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;

        case 32:
          reduce5<T, 32>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;

        case 16:
          reduce5<T, 16>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;

        case 8:
          reduce5<T, 8>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;

        case 4:
          reduce5<T, 4>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;

        case 2:
          reduce5<T, 2>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;

        case 1:
          reduce5<T, 1>
              <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
          break;
      }

      break;

    case 6:
      if (isPow2(size)) {
        switch (threads) {
          case 512:
            reduce6<T, 512, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 256:
            reduce6<T, 256, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 128:
            reduce6<T, 128, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 64:
            reduce6<T, 64, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 32:
            reduce6<T, 32, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 16:
            reduce6<T, 16, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 8:
            reduce6<T, 8, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 4:
            reduce6<T, 4, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 2:
            reduce6<T, 2, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 1:
            reduce6<T, 1, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;
        }
      } else {
        switch (threads) {
          case 512:
            reduce6<T, 512, false>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 256:
            reduce6<T, 256, false>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 128:
            reduce6<T, 128, false>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 64:
            reduce6<T, 64, false>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 32:
            reduce6<T, 32, false>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 16:
            reduce6<T, 16, false>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 8:
            reduce6<T, 8, false>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 4:
            reduce6<T, 4, false>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 2:
            reduce6<T, 2, false>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 1:
            reduce6<T, 1, false>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;
        }
      }

      break;

    case 7:
      // For reduce7 kernel we require only blockSize/warpSize
      // number of elements in shared memory
      smemSize = ((threads / 32) + 1) * sizeof(T);
      if (isPow2(size)) {
        switch (threads) {
          case 1024:
            reduce7<T, 1024, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;
          case 512:
            reduce7<T, 512, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 256:
            reduce7<T, 256, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 128:
            reduce7<T, 128, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 64:
            reduce7<T, 64, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 32:
            reduce7<T, 32, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 16:
            reduce7<T, 16, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 8:
            reduce7<T, 8, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 4:
            reduce7<T, 4, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 2:
            reduce7<T, 2, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;

          case 1:
            reduce7<T, 1, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
            break;
        }
*/

//
//   Transposition adapted from NVIDIA samples
//

// Each block transposes/copies a tile of TILE_DIM x TILE_DIM elements
// using TILE_DIM x BLOCK_ROWS threads, so that each thread transposes
// TILE_DIM/BLOCK_ROWS elements.  TILE_DIM must be an integral multiple of
// BLOCK_ROWS

#define TILE_DIM 16
#define BLOCK_ROWS 16

// This sample assumes that MATRIX_SIZE_X = MATRIX_SIZE_Y
int MATRIX_SIZE_X = 1024;
int MATRIX_SIZE_Y = 1024;
int MUL_FACTOR = TILE_DIM;

#define FLOOR(a, b) (a - (a % b))

// Compute the tile size necessary to illustrate performance cases for SM20+
// hardware
int MAX_TILES = (FLOOR(MATRIX_SIZE_X, 512) * FLOOR(MATRIX_SIZE_Y, 512)) /
                (TILE_DIM * TILE_DIM);

// Number of repetitions used for timing.  Two sets of repetitions are
// performed: 1) over kernel launches and 2) inside the kernel over just the
// loads and stores

#define NUM_REPS 100

// -------------------------------------------------------
// Copies
// width and height must be integral multiples of TILE_DIM
// -------------------------------------------------------

__global__ void copy_r4(float * __restrict__ out, 
                     const float * __restrict__ in, 
                     int width, 
                     int height) {
  int32_t xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int32_t yIndex = blockIdx.y * TILE_DIM + threadIdx.y;

  int32_t index = xIndex + width * yIndex;

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
       out[index + i * width] = in[index + i * width];
  }
}

__global__ void copy_i4(int32_t * __restrict__ out, 
                        const int32_t * __restrict__ in, 
                        int32_t width, 
                        int32_t height) {
  int32_t xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int32_t yIndex = blockIdx.y * TILE_DIM + threadIdx.y;

  int32_t index = xIndex + width * yIndex;

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
       out[index + i * width] = in[index + i * width];
  }
}


__global__ void copy_c4(cuComplex * __restrict__ out, 
                        const cuComplex * __restrict__ in, 
                        int width, 
                        int height) {
  int32_t xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int32_t yIndex = blockIdx.y * TILE_DIM + threadIdx.y;

  int32_t index = xIndex + width * yIndex;

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
       out[index + i * width] = in[index + i * width];
  }
}


__global__ void copy_shared_mem_r4(float * __restrict__ out, 
                                   const float * __restrict__ in, 
                                   int32_t width,
                                   int32_t height) {
  // Handle to thread block group
  cg::thread_block cta = cg::this_thread_block();
  __shared__ float tile[TILE_DIM][TILE_DIM];

  int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;

  int index = xIndex + width * yIndex;

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
    if (xIndex < width && yIndex < height) {
      tile[threadIdx.y][threadIdx.x] = in[index];
    }
  }

  cg::sync(cta);

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
    if (xIndex < height && yIndex < width) {
      out[index] = tile[threadIdx.y][threadIdx.x];
    }
  }
}


__global__ void copy_shared_mem_i4(int32_t * __restrict__ out, 
                                   const int32_t * __restrict__ in, 
                                   int32_t width,
                                   int32_t height) {
  // Handle to thread block group
  cg::thread_block cta = cg::this_thread_block();
  __shared__ int32_t tile[TILE_DIM][TILE_DIM];

  int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;

  int index = xIndex + width * yIndex;

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
    if (xIndex < width && yIndex < height) {
      tile[threadIdx.y][threadIdx.x] = in[index];
    }
  }

  cg::sync(cta);

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
    if (xIndex < height && yIndex < width) {
      out[index] = tile[threadIdx.y][threadIdx.x];
    }
  }
}


__global__ void copy_shared_mem_c4(cuComplex * __restrict__ out, 
                                   const cuComplex * __restrict__ in, 
                                   int32_t width,
                                   int32_t height) {
  // Handle to thread block group
  cg::thread_block cta = cg::this_thread_block();
  __shared__ cuComplex tile[TILE_DIM][TILE_DIM];

  int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;

  int index = xIndex + width * yIndex;

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
    if (xIndex < width && yIndex < height) {
      tile[threadIdx.y][threadIdx.x] = in[index];
    }
  }

  cg::sync(cta);

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
    if (xIndex < height && yIndex < width) {
      out[index] = tile[threadIdx.y][threadIdx.x];
    }
  }
}

// Coalesced transpose with no bank conflicts

_global__ void transpose_r4_v1(float * __restrict__ out, 
                               const float * __restrict__ in, 
                               int32_t width,
                               int32_t height) {
  // Handle to thread block group
  cg::thread_block cta = cg::this_thread_block();
  __shared__ float tile[TILE_DIM][TILE_DIM + 1];

  int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;
  int index_in = xIndex + (yIndex)*width;

  xIndex = blockIdx.y * TILE_DIM + threadIdx.x;
  yIndex = blockIdx.x * TILE_DIM + threadIdx.y;
  int index_out = xIndex + (yIndex)*height;

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
    tile[threadIdx.y + i][threadIdx.x] = in[index_in + i * width];
  }

  cg::sync(cta);

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
    out[index_out + i * height] = tile[threadIdx.x][threadIdx.y + i];
  }
}


_global__ void transpose_i4_v1(int32_t * __restrict__ out, 
                               const int32_t * __restrict__ in, 
                               int32_t width,
                               int32_t height) {
  // Handle to thread block group
  cg::thread_block cta = cg::this_thread_block();
  __shared__ int32_t tile[TILE_DIM][TILE_DIM + 1];

  int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;
  int index_in = xIndex + (yIndex)*width;

  xIndex = blockIdx.y * TILE_DIM + threadIdx.x;
  yIndex = blockIdx.x * TILE_DIM + threadIdx.y;
  int index_out = xIndex + (yIndex)*height;

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
    tile[threadIdx.y + i][threadIdx.x] = in[index_in + i * width];
  }

  cg::sync(cta);

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
    out[index_out + i * height] = tile[threadIdx.x][threadIdx.y + i];
  }
}


_global__ void transpose_c4_v1(cuComplex * __restrict__ out, 
                               const cuComplex * __restrict__ in, 
                               int32_t width,
                               int32_t height) {
  // Handle to thread block group
  cg::thread_block cta = cg::this_thread_block();
  __shared__ float tile[TILE_DIM][TILE_DIM + 1];

  int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;
  int index_in = xIndex + (yIndex)*width;

  xIndex = blockIdx.y * TILE_DIM + threadIdx.x;
  yIndex = blockIdx.x * TILE_DIM + threadIdx.y;
  int index_out = xIndex + (yIndex)*height;

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
    tile[threadIdx.y + i][threadIdx.x] = in[index_in + i * width];
  }

  cg::sync(cta);

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
    out[index_out + i * height] = tile[threadIdx.x][threadIdx.y + i];
  }
}


__global__ void transpose_r4_v2(float * __restrict__ out, 
                                const float * __restrict__ in, 
                                int32_t width,
                                int32_t height) {
  // Handle to thread block group
  cg::thread_block cta = cg::this_thread_block();
  __shared__ float tile[TILE_DIM][TILE_DIM + 1];

  int blockIdx_x, blockIdx_y;

  // do diagonal reordering
  if (width == height) {
    blockIdx_y = blockIdx.x;
    blockIdx_x = (blockIdx.x + blockIdx.y) % gridDim.x;
  } else {
    int bid = blockIdx.x + gridDim.x * blockIdx.y;
    blockIdx_y = bid % gridDim.y;
    blockIdx_x = ((bid / gridDim.y) + blockIdx_y) % gridDim.x;
  }

  // from here on the code is same as previous kernel except blockIdx_x replaces
  // blockIdx.x and similarly for y

  int xIndex = blockIdx_x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx_y * TILE_DIM + threadIdx.y;
  int index_in = xIndex + (yIndex)*width;

  xIndex = blockIdx_y * TILE_DIM + threadIdx.x;
  yIndex = blockIdx_x * TILE_DIM + threadIdx.y;
  int index_out = xIndex + (yIndex)*height;

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
    tile[threadIdx.y + i][threadIdx.x] = in[index_in + i * width];
  }

  cg::sync(cta);

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
       out[index_out + i * height] = tile[threadIdx.x][threadIdx.y + i];
  }
}


__global__ void transpose_i4_v2(int32_t * __restrict__ out, 
                                const int32_t * __restrict__ in, 
                                int32_t width,
                                int32_t height) {
  // Handle to thread block group
  cg::thread_block cta = cg::this_thread_block();
  __shared__ int32_t tile[TILE_DIM][TILE_DIM + 1];

  int blockIdx_x, blockIdx_y;

  // do diagonal reordering
  if (width == height) {
    blockIdx_y = blockIdx.x;
    blockIdx_x = (blockIdx.x + blockIdx.y) % gridDim.x;
  } else {
    int bid = blockIdx.x + gridDim.x * blockIdx.y;
    blockIdx_y = bid % gridDim.y;
    blockIdx_x = ((bid / gridDim.y) + blockIdx_y) % gridDim.x;
  }

  // from here on the code is same as previous kernel except blockIdx_x replaces
  // blockIdx.x and similarly for y

  int xIndex = blockIdx_x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx_y * TILE_DIM + threadIdx.y;
  int index_in = xIndex + (yIndex)*width;

  xIndex = blockIdx_y * TILE_DIM + threadIdx.x;
  yIndex = blockIdx_x * TILE_DIM + threadIdx.y;
  int index_out = xIndex + (yIndex)*height;

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
    tile[threadIdx.y + i][threadIdx.x] = in[index_in + i * width];
  }

  cg::sync(cta);

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
       out[index_out + i * height] = tile[threadIdx.x][threadIdx.y + i];
  }
}


__global__ void transpose_c4_v1(cuComplex * __restrict__ out, 
                                const cuComplex * __restrict__ in, 
                                int32_t width,
                                int32_t height) {
  // Handle to thread block group
  cg::thread_block cta = cg::this_thread_block();
  __shared__ cuComplex tile[TILE_DIM][TILE_DIM + 1];

  int blockIdx_x, blockIdx_y;

  // do diagonal reordering
  if (width == height) {
    blockIdx_y = blockIdx.x;
    blockIdx_x = (blockIdx.x + blockIdx.y) % gridDim.x;
  } else {
    int bid = blockIdx.x + gridDim.x * blockIdx.y;
    blockIdx_y = bid % gridDim.y;
    blockIdx_x = ((bid / gridDim.y) + blockIdx_y) % gridDim.x;
  }

  // from here on the code is same as previous kernel except blockIdx_x replaces
  // blockIdx.x and similarly for y

  int xIndex = blockIdx_x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx_y * TILE_DIM + threadIdx.y;
  int index_in = xIndex + (yIndex)*width;

  xIndex = blockIdx_y * TILE_DIM + threadIdx.x;
  yIndex = blockIdx_x * TILE_DIM + threadIdx.y;
  int index_out = xIndex + (yIndex)*height;

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
    tile[threadIdx.y + i][threadIdx.x] = in[index_in + i * width];
  }

  cg::sync(cta);

  for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
       out[index_out + i * height] = tile[threadIdx.x][threadIdx.y + i];
  }
}


__global__
void atomic_reduction(float  * __restrict__  out,
                      const float * __restrict__ in,
                      const int32_t n) {
     const int32_t tid = threadIdx.x;
     float partsum = 0.0f;
     int32_t i = blockIdx.x*blockDIm.x+tid;
     const int32_t stride = blockDIm.x*gridDIm.x;
     for(; i < n; i += stride) {
         partsum += in[i];
    }
    atomicAdd(&out[i],partsum);
}



#define R4_BUFF_LEN  256
#define C4_BUFF_LEN  128


__global__ void
device_copy_r4(float * __restrict__ to
               const float * __restrict__ from,
               const int32_t n) {
     float buf[R4_BUFF_LEN];
     int32_t tid = R4_BUFF_LEN*blockIdx.x*blockDim.x+threadIdx.x;
     const int32_t stride = R4_BUFF_LEN*blockDIm.x*gridDIm.x;
     int32_t i;
     for(i = tid; i < n-stride i += stride) {
         for(int32_t j = 0; j < R4_BUFF_LEN; ++j) {
             const int32_t idx = i+j*blockDim.x;
             buf[j] = from[idx];
         }
         for(int32_t j = 0; j < R4_BUFF_LEN; ++j) {
             const int32_t idx = i+j*blockDim.x;
             to[idx] = buf[j];
         }
     }
     
     for(int32_t j = 0; j < R4_BUFF_LEN; ++j) {
         for(int32_t j = 0; j < R4_BUFF_LEN; ++j) {
             const int32_t idx = i+j*blockDIm.x;
             if(idx<n) buf[j] = from[idx];
         }
         for(int32_t j = 0; j < R4_BUFF_LEN; ++j) {
             const int32_t idx = i+j*blockDim.x;
             if(idx<n) to[idx] = buf[j];
         }
     }
}


__global__ void
device_copy_c4(cuComplex * __restrict__ to,
               const cuComplex * __restrict__ from,
               const int32_t n) {

     cuComplex buf[C4_BUFF_LEN];
     int32_t tid = C4_BUFF_LEN*blockIdx.x*blockDim.x+threadIdx.x;
     const int32_t stride = C4_BUFF_LEN*blockDIm.x*gridDIm.x;
     int32_t i;
     for(i = tid; i < n-stride i += stride) {
         for(int32_t j = 0; j < C4_BUFF_LEN; ++j) {
             const int32_t idx = i+j*blockDim.x;
             buf[j] = from[idx];
         }
         for(int32_t j = 0; j < C4_BUFF_LEN; ++j) {
             const int32_t idx = i+j*blockDim.x;
             to[idx] = buf[j];
         }
     }
     
     for(int32_t j = 0; j < C4_BUFF_LEN; ++j) {
         for(int32_t j = 0; j < C4_BUFF_LEN; ++j) {
             const int32_t idx = i+j*blockDIm.x;
             if(idx<n) buf[j] = from[idx];
         }
         for(int32_t j = 0; j < C4_BUFF_LEN; ++j) {
             const int32_t idx = i+j*blockDim.x;
             if(idx<n) to[idx] = buf[j];
         }
     }
}

#endif /*__GMS_CUDA_UTILS_CU__*/
