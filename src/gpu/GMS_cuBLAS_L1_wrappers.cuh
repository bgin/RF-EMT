
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
#ifndef __GMS_CUBLAS_L1_WRAPPERS_CUH__
#define __GMS_CUBLAS_L1_WRAPPERS_CUH__


/*
     cuBLAS Level-1 wrappers 
*/

#include <stdint>



extern "C"
void
cuBLAS_Isamax_iface(const float * __restrict,
                    const int32_t,
                    const int32_t,
                    int32_t *     __restrict,
                    cudaError_t * __restrict,
                    int32_t * __restrict,
                    uint64_t * __restrict) 
                                        __attribute__((hot))
                                        __attribute__((aligned(32)));
extern "C"
void
cuBLAS_Icamax_iface(const cuComplex * __restrict,
                    const int32_t,
                    const int32_t,
                    int32_t *  __restrict,
                    cudaError_t * __restrict,
                    int32_t * __restrict,
                    uint64_t * __restrict)  
                                        __attribute__((hot))
                                        __attribute__((aligned(32)));


extern "C"
void
cuBLAS_Isamin_iface(const float * __restrict,
                    const int32_t,
                    const int32_t,
                    int32_t *     __restrict,
                    cudaError_t * __restrict,
                    int32_t * __restrict,
                    uint64_t * __restrict) 
                                        __attribute__((hot))
                                        __attribute__((aligned(32)));

extern "C"
void
cuBLAS_Icamin_iface(const cuComplex * __restrict,
                    const int32_t,
                    const int32_t,
                    int32_t *  __restrict,
                    cudaError_t * __restrict,
                    int32_t * __restrict,
                    uint64_t * __restrict)  
                                        __attribute__((hot))
                                        __attribute__((aligned(32)));

extern "C"
void
cuBLAS_Sasum_iface(const float * __restrict,
                   const int32_t,
                   const int32_t,
                   float * __restrict,
                   cudaError_t * __restrict,
                   int32_t * __restrict,
                   uint64_t * __restrict) 
                                        __attribute__((hot))
                                        __attribute__((aligned(32)));


extern "C"
void
cuBLAS_Scasum_iface(const cuComplex * __restrict,
                   const int32_t,
                   const int32_t,
                   cuComplex * __restrict,
                   cudaError_t * __restrict,
                   int32_t * __restrict,
                   uint64_t * __restrict) 
                                        __attribute__((hot))
                                        __attribute__((aligned(32)));




#endif /*__GMS_CUBLAS_L1_WRAPPERS_CUH__*/
