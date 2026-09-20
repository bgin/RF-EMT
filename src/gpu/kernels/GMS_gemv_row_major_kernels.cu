
#include <cstdint>
#include "GMS_atomicops.cuh"
// Modified by Bernard Gingold, contact: bgin@gmail.com
/*
    GLAS is licensed under the The MIT License:
Copyright (c) 2011 Hans Henrik Brandenborg Sørensen, DTU
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
// Row-major matrix-vector multiplication kernel.
//
// - Designed for architectures with atomicCAS, i.e., sm_11 and higher.
// - Split into 4 separate device kernels.
// -----------------
// |         |     |
// |    1    |  2  |
// |         |     |
// |---------------|
// |    3    |  4  |
// -----------------
// - Valid for all cases.
// - Uses atomicCAS (sm_11-sm_13) or atomicAdd (>=sm_20) for grid level reduction.
// - Unrolling of inner loops
//     * Given by a constant: full unroll
//     * Given by a variable: unroll level 4
//     1) to reduce instruction count from loop
//     2) to collect shared loads unsigned into 128 bit transactions (compiler should do that)
//



inline __device__ void block_reduce_r4(float * __restrict__ sdata, 
                                       volatile float * __restrict__ smem, 
                                       float sum, 
                                       int tid,
                                       const int32_t BLOCKSIZE_n)
{
	sdata[tid] = sum;
	__syncthreads();
	if (BLOCKSIZE_n > 1)
	{
		smem = sdata;
		if (BLOCKSIZE_n >= 512) { if (tid < 256) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
		if (BLOCKSIZE_n >= 256) { if (tid < 128) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
		if (BLOCKSIZE_n == 224) { if (tid < 96) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
		if (BLOCKSIZE_n == 192) { if (tid < 64) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
		if (BLOCKSIZE_n >= 128) { if (tid < 64) { smem[tid] += smem[tid + 64]; } __syncthreads(); }
		if (BLOCKSIZE_n == 160) { if (tid < 32) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
		if (BLOCKSIZE_n == 96) { if (tid < 32) { smem[tid] += smem[tid + 64]; } __syncthreads(); }
		if (BLOCKSIZE_n >= 64) { if (tid < 32) { smem[tid] += smem[tid + 32]; } }
		if (BLOCKSIZE_n >= 32) { if (tid < 16) { smem[tid] += smem[tid + 16]; } }
		if (BLOCKSIZE_n >= 16) { if (tid < 8) { smem[tid] += smem[tid + 8]; } }
		if (BLOCKSIZE_n >= 8) { if (tid < 4) { smem[tid] += smem[tid + 4]; } }
		if (BLOCKSIZE_n >= 4) { if (tid < 2) { smem[tid] += smem[tid + 2]; } }
		if (BLOCKSIZE_n >= 2) { if (tid < 1) { smem[tid] += smem[tid + 1]; } }
	} 
}

inline __device__ void block_reduce_r8(double * __restrict__ sdata, 
                                       volatile double * __restrict__ smem, 
                                       double sum, 
                                       const int32_t tid,
                                       const int32_t BLOCKSIZE_n)
{
	sdata[tid] = sum;
	__syncthreads();
	if (BLOCKSIZE_n > 1)
	{
		smem = sdata;
		if (BLOCKSIZE_n >= 512) { if (tid < 256) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
		if (BLOCKSIZE_n >= 256) { if (tid < 128) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
		if (BLOCKSIZE_n == 224) { if (tid < 96) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
		if (BLOCKSIZE_n == 192) { if (tid < 64) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
		if (BLOCKSIZE_n >= 128) { if (tid < 64) { smem[tid] += smem[tid + 64]; } __syncthreads(); }
		if (BLOCKSIZE_n == 160) { if (tid < 32) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
		if (BLOCKSIZE_n == 96) { if (tid < 32) { smem[tid] += smem[tid + 64]; } __syncthreads(); }
		if (BLOCKSIZE_n >= 64) { if (tid < 32) { smem[tid] += smem[tid + 32]; } }
		if (BLOCKSIZE_n >= 32) { if (tid < 16) { smem[tid] += smem[tid + 16]; } }
		if (BLOCKSIZE_n >= 16) { if (tid < 8) { smem[tid] += smem[tid + 8]; } }
		if (BLOCKSIZE_n >= 8) { if (tid < 4) { smem[tid] += smem[tid + 4]; } }
		if (BLOCKSIZE_n >= 4) { if (tid < 2) { smem[tid] += smem[tid + 2]; } }
		if (BLOCKSIZE_n >= 2) { if (tid < 1) { smem[tid] += smem[tid + 1]; } }
	} 
}




inline __device__ void block_reduce_2d_r4(float * __restrict__ sdata, 
                                          volatile float * __restrict__ smem, 
                                          float sum, 
                                          const int32_t tid, 
                                          const int32_t tidx,
                                          const int32_t BLOCKSIZE_n)
{
	sdata[tid] = sum;
	__syncthreads();
	if (BLOCKSIZE_n > 1)
	{
		smem = sdata;
		if (BLOCKSIZE_n >= 512) { if (tidx < 256) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
		if (BLOCKSIZE_n >= 256) { if (tidx < 128) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
		if (BLOCKSIZE_n == 224) { if (tidx < 96) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
		if (BLOCKSIZE_n == 192) { if (tidx < 64) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
		if (BLOCKSIZE_n >= 128) { if (tidx < 64) { smem[tid] += smem[tid + 64]; } __syncthreads(); }
		if (BLOCKSIZE_n == 160) { if (tidx < 32) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
		if (BLOCKSIZE_n == 96) { if (tidx < 32) { smem[tid] += smem[tid + 64]; } __syncthreads(); }
		if (BLOCKSIZE_n >= 64) { if (tidx < 32) { smem[tid] += smem[tid + 32]; } }
		if (BLOCKSIZE_n >= 32) { if (tidx < 16) { smem[tid] += smem[tid + 16]; } }
		if (BLOCKSIZE_n >= 16) { if (tidx < 8) { smem[tid] += smem[tid + 8]; } }
		if (BLOCKSIZE_n >= 8) { if (tidx < 4) { smem[tid] += smem[tid + 4]; } }
		if (BLOCKSIZE_n >= 4) { if (tidx < 2) { smem[tid] += smem[tid + 2]; } }
		if (BLOCKSIZE_n >= 2) { if (tidx < 1) { smem[tid] += smem[tid + 1]; } }
	} 
}


inline __device__ void block_reduce_2d_r8(double * __restrict__ sdata, 
                                          volatile double * __restrict__ smem, 
                                          double sum, 
                                          const int32_t tid, 
                                          const int32_t tidx,
                                          const int32_t BLOCKSIZE_n)
{
	sdata[tid] = sum;
	__syncthreads();
	if (BLOCKSIZE_n > 1)
	{
		smem = sdata;
		if (BLOCKSIZE_n >= 512) { if (tidx < 256) { smem[tid] += smem[tid + 256]; } __syncthreads(); }
		if (BLOCKSIZE_n >= 256) { if (tidx < 128) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
		if (BLOCKSIZE_n == 224) { if (tidx < 96) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
		if (BLOCKSIZE_n == 192) { if (tidx < 64) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
		if (BLOCKSIZE_n >= 128) { if (tidx < 64) { smem[tid] += smem[tid + 64]; } __syncthreads(); }
		if (BLOCKSIZE_n == 160) { if (tidx < 32) { smem[tid] += smem[tid + 128]; } __syncthreads(); }
		if (BLOCKSIZE_n == 96) { if (tidx < 32) { smem[tid] += smem[tid + 64]; } __syncthreads(); }
		if (BLOCKSIZE_n >= 64) { if (tidx < 32) { smem[tid] += smem[tid + 32]; } }
		if (BLOCKSIZE_n >= 32) { if (tidx < 16) { smem[tid] += smem[tid + 16]; } }
		if (BLOCKSIZE_n >= 16) { if (tidx < 8) { smem[tid] += smem[tid + 8]; } }
		if (BLOCKSIZE_n >= 8) { if (tidx < 4) { smem[tid] += smem[tid + 4]; } }
		if (BLOCKSIZE_n >= 4) { if (tidx < 2) { smem[tid] += smem[tid + 2]; } }
		if (BLOCKSIZE_n >= 2) { if (tidx < 1) { smem[tid] += smem[tid + 1]; } }
	} 
}



__global__ void gemv_one_block_per_row_kernel_r4( float alpha, 
                                                  const float* __restrict__ A, 
                                                  const int32_t lda, 
                                                  const float* __restrict__ x, 
                                                  float beta, 
                                                  float* __restrict__ y, 
                                                  const int32_t m, 
                                                  const int32_t n, 
                                                  const int32_t mtile, 
                                                  const int32_t ntile,
                                                  const int32_t BLOCKSIZE_n, 
                                                  const int32_t WORKSIZE_m, 
                                                  const int32_t WORKSIZE_n )
{
	int tid = threadIdx.x;
	__shared__ float sdata[BLOCKSIZE_n];
	volatile float * __restrict__ smem;
	float sum = 0.0f;
	for(int i = 0; i < ntile; i += BLOCKSIZE_n)
	{
		sum += A[tid + i + lda * blockIdx.x] * x[tid + i];
	}
	if (tid + ntile <  n)
	{
		sum += A[tid + ntile + lda * blockIdx.x] * x[tid + ntile];
	}
	block_reduce_r4(sdata, smem, sum, tid, BLOCKSIZE_n);
	if (tid == 0) y[blockIdx.x] = alpha * sdata[0] + beta * y[blockIdx.x];
}


__global__ void gemv_one_block_per_row_kernel_r8( double alpha, 
                                                  const double* __restrict__ A, 
                                                  const int32_t lda, 
                                                  const double* __restrict__ x, 
                                                  double beta, 
                                                  double * __restrict__ y, 
                                                  const int32_t m, 
                                                  const int32_t n, 
                                                  const int32_t mtile, 
                                                  const int32_t ntile,
                                                  const int32_t BLOCKSIZE_n, 
                                                  const int32_t WORKSIZE_m, 
                                                  const int32_t WORKSIZE_n )
{
	int tid = threadIdx.x;
	__shared__ double sdata[BLOCKSIZE_n];
	volatile double * __restrict__ smem;
	float sum = 0.0;
	for(int i = 0; i < ntile; i += BLOCKSIZE_n)
	{
		sum += A[tid + i + lda * blockIdx.x] * x[tid + i];
	}
	if (tid + ntile <  n)
	{
		sum += A[tid + ntile + lda * blockIdx.x] * x[tid + ntile];
	}
	block_reduce_r8(sdata, smem, sum, tid, BLOCKSIZE_n);
	if (tid == 0) y[blockIdx.x] = alpha * sdata[0] + beta * y[blockIdx.x];
}






__global__ void sgemv_one_block_per_row_kernel2_r4(float alpha, 
                                                   const float* __restrict__ A, 
                                                   int lda, 
                                                   const float* __restrict__ x, 
                                                   float beta, 
                                                   float* __restrict__ y, 
                                                   const int32_t m, 
                                                   const int32_t n, 
                                                   const int32_t mtile, 
                                                   const int32_t ntile,
                                                   const int32_t BLOCKSIZE_n, 
                                                   const int32_t WORKSIZE_m, 
                                                   const int32_t WORKSIZE_n)
{
	int32_t tid = threadIdx.x;
	__shared__ float sdata[BLOCKSIZE_n];
	A += tid + lda * blockIdx.x;
	x += tid;
	volatile float *smem;
	float sum = 0.0f;
	for(int i = 0; i < ntile; i += BLOCKSIZE_n)
	{
		sum += A[0] * x[0]; A += BLOCKSIZE_n; x += BLOCKSIZE_n;
	}
	if (tid + ntile <  n)
	{
		sum += A[0] * x[0];
	}
	block_reduce_r4(sdata, smem, sum, tid, BLOCKSIZE_n);
	if (tid == 0) y[blockIdx.x] = alpha * sdata[0] + beta * y[blockIdx.x];
}


__global__ void sgemv_one_block_per_row_kernel2_r8(double alpha, 
                                                   const double* __restrict__ A, 
                                                   int lda, 
                                                   const double* __restrict__ x, 
                                                   double beta, 
                                                   double* __restrict__ y, 
                                                   const int32_t m, 
                                                   const int32_t n, 
                                                   const int32_t mtile, 
                                                   const int32_t ntile,
                                                   const int32_t BLOCKSIZE_n, 
                                                   const int32_t WORKSIZE_m, 
                                                   const int32_t WORKSIZE_n)
{
	int32_t tid = threadIdx.x;
	__shared__ double sdata[BLOCKSIZE_n];
	A += tid + lda * blockIdx.x;
	x += tid;
	volatile double *smem;
	float sum = 0.0;
	for(int i = 0; i < ntile; i += BLOCKSIZE_n)
	{
		sum += A[0] * x[0]; A += BLOCKSIZE_n; x += BLOCKSIZE_n;
	}
	if (tid + ntile <  n)
	{
		sum += A[0] * x[0];
	}
	block_reduce_r8(sdata, smem, sum, tid, BLOCKSIZE_n);
	if (tid == 0) y[blockIdx.x] = alpha * sdata[0] + beta * y[blockIdx.x];
}

/*template <int BLOCKSIZE_m, int BLOCKSIZE_n, int WORKSIZE_m, int WORKSIZE_n>
__launch_bounds__(BLOCKSIZE_n, MIN_BLOCKS_PER_MP)
__global__ void sgemv_several_blocks_several_rows_kernel(float alpha, const float* A, int lda, const float* x, float beta, float* y, int m, int n, int mtile, int ntile)
{
	int tid = threadIdx.x;
	int row = WORKSIZE_m * blockIdx.x;
	int col = BLOCKSIZE_n * WORKSIZE_n * blockIdx.y;
	__shared__ float sdata[BLOCKSIZE_n];
	A += tid + row * lda + col;
	x += tid + col;
	y += row;
	volatile float *smem;
	initialize_sums(WORKSIZE_m);
	if (row !=  mtile)
	{
		if (col != ntile)
		{
			//
			// The top-left fully tiled part of the Matrix
			//
            #pragma unroll
			inner_loop(BLOCKSIZE_n * WORKSIZE_n);
			reduce_and_atomicAdd_sums(BLOCKSIZE_n);
		}
		else
		{
			//
			// The top-right m-tiled part of the Matrix
			//
			int nrest = n - ntile;
			int nblks = (nrest / BLOCKSIZE_n) * BLOCKSIZE_n;
            #pragma unroll 4
			inner_loop(nblks);
			if (tid + nblks < nrest)
			{
				add_to_sums(lda, WORKSIZE_m);
			}
			reduce_and_atomicAdd_sums(BLOCKSIZE_n);
		}
	}
	else
	{
		if (col != ntile)
		{	
			//
			// The bottom-left n-tiled part of the Matrix
			//
            #pragma unroll
			inner_loop_check_bounds(row, m, BLOCKSIZE_n * WORKSIZE_n);
			reduce_and_atomicAdd_sums_check_bounds(row, m, BLOCKSIZE_n);
		}
		else
		{
			//
			// The bottom-right rest of the Matrix
			//
			int nrest = n - ntile;
			int nblks = (nrest / BLOCKSIZE_n) * BLOCKSIZE_n;
            #pragma unroll 4
			inner_loop_check_bounds(row, m, nblks);
			if (tid + nblks < nrest)
			{
				add_to_sums_check_bounds(row, m, lda, WORKSIZE_m);
			}
			reduce_and_atomicAdd_sums_check_bounds(row, m, BLOCKSIZE_n);
		}
	}
}

template <int BLOCKSIZE_m, int BLOCKSIZE_n, int WORKSIZE_m, int WORKSIZE_n>
__launch_bounds__(BLOCKSIZE_m * BLOCKSIZE_n, MIN_BLOCKS_PER_MP)
__global__ void sgemv_several_blocks_per_row_2d_kernel(float alpha, const float* A, int lda, const float* x, float beta, float* y, int m, int n, int mtile, int ntile)
{
	int tidx = threadIdx.x;
	int tidy = threadIdx.y;
	int tid = tidx + BLOCKSIZE_n * tidy;
	int row = BLOCKSIZE_m * WORKSIZE_m * blockIdx.y;
	int col = BLOCKSIZE_n * WORKSIZE_n * blockIdx.x;
	int therow = row + tidy * WORKSIZE_m;
	__shared__ float sdata[BLOCKSIZE_m * BLOCKSIZE_n];
	A += tidx + therow * lda + col;
	x += tidx + col;
	y += therow;
	volatile float *smem;
	initialize_sums(WORKSIZE_m);
	if (row !=  mtile)
	{
		if (col != ntile)
		{
			//
			// The top-left fully tiled part of the Matrix
			//
            #pragma unroll
			inner_loop(BLOCKSIZE_n * WORKSIZE_n);
			reduce_2d_and_atomicAdd_sums(BLOCKSIZE_n);
			//printf("HER tidx=%d tidy=%d tid=%d m=%d therow=%d sum[0]=%f\n", tidx, tidy, tid, m, therow, sum[0]);
		}
		else
		{
			//
			// The top-right m-tiled part of the Matrix
			//
			int nrest = n - ntile;
			int nblks = (nrest / BLOCKSIZE_n) * BLOCKSIZE_n;
            #pragma unroll 4
			inner_loop(nblks);
			if (tidx + nblks < nrest)
			{
				add_to_sums(lda, WORKSIZE_m);
			}
			reduce_2d_and_atomicAdd_sums(BLOCKSIZE_n);
		}
	}
	else if (therow < m)
	{
		if (col != ntile)
		{	
			//
			// The bottom-left n-tiled part of the Matrix
			//
            #pragma unroll
			inner_loop_check_bounds(therow, m, BLOCKSIZE_n * WORKSIZE_n);
			reduce_2d_and_atomicAdd_sums_check_bounds(row, m, BLOCKSIZE_n);
		}
		else
		{
			//
			// The bottom-right rest of the Matrix
			//
			int nrest = n - ntile;
			int nblks = (nrest / BLOCKSIZE_n) * BLOCKSIZE_n;
            #pragma unroll 4
			inner_loop_check_bounds(therow, m, nblks);
			if (tidx + nblks < nrest)
			{
				add_to_sums_check_bounds(therow, m, lda, WORKSIZE_m);
			}
			reduce_2d_and_atomicAdd_sums_check_bounds(therow, m, BLOCKSIZE_n);
		}
	}
}


template <int BLOCKSIZE_m, int BLOCKSIZE_n, int WORKSIZE_m, int WORKSIZE_n>
__launch_bounds__(BLOCKSIZE_m * BLOCKSIZE_n, MIN_BLOCKS_PER_MP)
__global__ void sgemv_xdim_block_per_row_2d_kernel(float alpha, const float* A, int lda, const float* x, float beta, float* y, int m, int n, int mtile, int ntile, int lds)
{
	int tidx = threadIdx.x;
	int tidy = threadIdx.y;
	int tid = tidx + BLOCKSIZE_n * tidy;
	int row = BLOCKSIZE_m * WORKSIZE_m * blockIdx.y;
	int col = BLOCKSIZE_n * WORKSIZE_n * blockIdx.x;
	int therow = row + tidy * WORKSIZE_m;
	extern __shared__ float sdata[];
	A += row * lda + col;
	y += therow;
	const float *x_s = x + tidx + col;
	/*
	  // If using the shared memory for x - no difference compared to cache
	float *x_s = sdata + lds * BLOCKSIZE_m * WORKSIZE_m;
	x += col;
	if (tid < n)
	{
		x_s[tid] = x[tid];
	}
	x_s += tidx;
	
	float *A_s = sdata + tidx + tidy * WORKSIZE_m * lds;
	if (row !=  mtile)
	{
		ntile = n * BLOCKSIZE_m * WORKSIZE_m;
	}
	else
	{
		ntile = n * (m - mtile);
	}
	if (lds == lda)
	{
		A += tid;
		float *A_s = sdata + tid;
		int nblks = (ntile / (BLOCKSIZE_m * BLOCKSIZE_n)) * (BLOCKSIZE_m * BLOCKSIZE_n);
        #pragma unroll 4
		for (int i = 0; i < nblks; i += BLOCKSIZE_m * BLOCKSIZE_n)
		{
			A_s[i] = A[i];
		}
		if (tid + nblks < ntile)
		{
			A_s[nblks] = A[nblks];
		}
	}
	else
	{
		A += tidx + tidy * WORKSIZE_m * lda;
        #pragma unroll 4
		for (int i = 0; i < n - tidx; i += BLOCKSIZE_n)
		{
			copy_A_sh(lds, lda, WORKSIZE_m);
		}
	}
	__syncthreads();
	volatile float *smem;
	initialize_sums(WORKSIZE_m);
	if (row !=  mtile)
	{
		//
		// The top-right m-tiled part of the Matrix
		//
		int nblks = (n / BLOCKSIZE_n) * BLOCKSIZE_n;
        #pragma unroll 4
		inner_loop_sh(nblks);
		if (tidx + nblks < n)
		{
			add_to_sums_sh(lds, WORKSIZE_m);
		}
		__syncthreads();
		reduce_2d_and_atomicAdd_sums(BLOCKSIZE_n);
	}
	else if (therow < m)
	{
		//
		// The bottom-right rest of the Matrix
		//
        #pragma unroll 4
		int nblks = (n / BLOCKSIZE_n) * BLOCKSIZE_n;
		inner_loop_check_bounds_sh(therow, m, nblks);
		if (tidx + nblks < n)
		{
			add_to_sums_check_bounds_sh(therow, m, lds, WORKSIZE_m);
		}
		__syncthreads();
		reduce_2d_and_atomicAdd_sums_check_bounds(therow, m, BLOCKSIZE_n);
	}
}
*/
