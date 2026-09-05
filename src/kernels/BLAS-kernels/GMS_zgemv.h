

#ifndef __GMS_ZGEMV_H__
#define __GMS_ZGEMV_H__

//=========================================================================
// Modified and optimized version of OpenBLAS zgemv and kernels zgemv_4
// Programmer: Bernard Gingold, contact: beniekg@gmail.com
// 07-02-2021 10:06 AM +00200
// Original copyright below
//========================================================================
/***************************************************************************
Copyright (c) 2014, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

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

#include <cstdint>
#include "GMS_config" // must define 'CONJ'

// This is important macro definition
// Set this to '1' in order to enable simd reduction in zgemv_function
// Vectorized reduction may cause (probably) performance reduction
// because of irregular stride advancing.
#if !defined(GMS_ZGEMV_EXPLICITLY_VECTORIZE)
#define GMS_ZGEMV_EXPLICITLY_VECTORIZE 0
#endif

// Kernels

#if !defined(ZGEMV_KERNEL_USE_PREFETCHT0)
    #define ZGEMV_KERNEL_USE_PREFETCHT0 1
#endif

__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
void zgemv_kernel_4x4(int32_t n,
                      double **__restrict ap,
		      double * __restrict x,
		      double * __restrict y);
         


__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
void zgemv_kernel_4x2(const int32_t n,
                      double ** __restrict ap,
		      double *  __restrict x,
		      double *  __restrict y); 

__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
void zgemv_kernel_4x1(const int32_t n,
                      double * __restrict ap,
		      double * __restrict x,
		      double * __restrict y); 
		      
__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
void add_y(const int32_t n,
           double * __restrict src,
	   double * __restrict dest,
	   const int32_t inc_dest,
	   const double alpha_r,
	   const double alpha_i); 
	   
	   


#define ZGEMV_NBMAX 1024

__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
void zgemv(const int32_t m,
           const int32_t n,
	   const double alpha_r,
	   const double alpha_i,
	   double * __restrict a,
	   const int32_t lda,
	   double * __restrict x,
	   const int32_t inc_x,
	   double * __restrict y,
	   const int32_t inc_y,
	   double * buffer); 






#endif /*__GMS_ZGEMV_H__*/
