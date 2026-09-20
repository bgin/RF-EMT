
#include <cstdint>
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


// Generic vector copy kernel.
//
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
// - Kernels are generic and work for all cases
//

// Version 1: copy as floats
// -----------------------------------------
// Threads per block    : BLOCKSIZE_n
// Work per thread      : WORKSIZE_n
// Shared mem per block : 0
//

__global__ void copy_kernel_r4(const float* __restrict__ x, 
                               float* __restrict__ y, 
                               const int32_t n, 
                               const int32_t ntile,
                               const int32_t BLOCKSIZE_n,
                               const int32_t WORKSIZE_m,
                               const int32_t WORKSIZE_n )
{
	const int32_t tid = threadIdx.x;
	const int32_t idx = blockIdx.x * BLOCKSIZE_n * WORKSIZE_n;
	x += idx + tid;
	y += idx + tid;
	if (idx != ntile)
	{
        #pragma unroll
		for (int i = 0; i < WORKSIZE_n; i++ )
		{
			y[0] = x[0]; x += BLOCKSIZE_n; y += BLOCKSIZE_n;
		}
	}
	else
	{
        #pragma unroll 4
		for (int i = 0; i < n - idx - tid; i += BLOCKSIZE_n)
		{
			y[i] = x[i];
		}
	}
}


__global__ void copy_kernel_r8(const double* __restrict__ x, 
                               double* __restrict__ y, 
                               const int32_t n, 
                               const int32_t ntile,
                               const int32_t BLOCKSIZE_n,
                               const int32_t WORKSIZE_m,
                               const int32_t WORKSIZE_n )
{
	const int32_t tid = threadIdx.x;
	const int32_t idx = blockIdx.x * BLOCKSIZE_n * WORKSIZE_n;
	x += idx + tid;
	y += idx + tid;
	if (idx != ntile)
	{
        #pragma unroll
		for (int i = 0; i < WORKSIZE_n; i++ )
		{
			y[0] = x[0]; x += BLOCKSIZE_n; y += BLOCKSIZE_n;
		}
	}
	else
	{
        #pragma unroll 4
		for (int i = 0; i < n - idx - tid; i += BLOCKSIZE_n)
		{
			y[i] = x[i];
		}
	}
}

// Version 2: copy as doubles (requires WORKSIZE_n & 1 == 0)
// -----------------------------------------
// Threads per block    : BLOCKSIZE_n
// Work per thread      : WORKSIZE_n / 2
// Shared mem per block : 0
//
/*template <int BLOCKSIZE_n, int WORKSIZE_m, int WORKSIZE_n>
__launch_bounds__(BLOCKSIZE_n, MIN_BLOCKS_PER_MP)
__global__ void scopy_kernel2(const float* x, float* y, int n, int ntile)
{
	const int tid = threadIdx.x;
	const int idx = blockIdx.x * BLOCKSIZE_n * WORKSIZE_n;
	x += idx;
	y += idx;
	const double *xd = (double *) x + tid;
	double *yd = (double *) y + tid;
	if (idx != ntile)
	{
        #pragma unroll
		for (int i = 0; i < (WORKSIZE_n / 2); i++ )
		{
			yd[0] = xd[0]; xd += BLOCKSIZE_n; yd += BLOCKSIZE_n;
		}
	}
	else
	{
        #pragma unroll 4
		for (int i = 0; i < ((n - idx) >> 1) - tid; i += BLOCKSIZE_n)
		{
			yd[i] = xd[i];
		}
	}
	if (blockIdx.x == 0 && tid == 0 && n & 1) y[n - 1] = x[n - 1];
}
*/
