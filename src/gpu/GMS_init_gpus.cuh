
#ifndef __GMS_INIT_GPUS_CUH__
#define __GMS_INIT_GPUS_CUH__


#include <stdint.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>


//
//	Declarations
//

//
//	Print device capabilities to screen.
//
void devcaps_to_screen(const int32_t, 
                       int32_t * , const bool);


//
//	Print device capabilities to file.
//
void devcaps_to_file(const int32_t, const char *, 
		     int32_t * , const bool );

//
//   Taken and adapted slightly from the cuda-samples   
//
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

 // Defines for GPU Architecture types (using the SM version to determine
  // the # of cores per SM
int32_t convertSMVer2Cores(const int32_t,
                           const int32_t);

  // Defines for GPU Architecture types (using the SM version to determine
  // the GPU Arch name)
const char * convertSMVer2ArchName(const int32_t,
                                   const int32_t);



int32_t gpuDeviceInit(const int32_t);

// The device with highest Gflops throughput
int32_t gpuGetMaxGflopsDevId();

// General check for CUDA GPU SM capabilities
bool checkCudaCapabilities(const int32_t,
                           const int32_t);

// Finds the integrated GPU compute-capable
int32_t findIntegratedGPU();

// Finds the best CUDA device
int32_t findCudaDevice();

void setBlocksThreadNumber(const int32_t,
                           const int32_t,
                           const int32_t,
                           const int32_t,
                           int32_t &,
                           int32_t &);






//
// Run integer and real arithmentic test.
// Host code.
//
cudaError_t  gpu_vec_add_tests(const32_t int);

// 
// Device code.
// simple vector addition.
//
__global__ void kvec_add_int32(int32_t * __restrict,
			       const int32_t * __restrict,
			       const int32_t * __restrict);

__global__ void kvec_add_float(float * __restrict,
			       const float * __restrict ,
			       const float * __restrict);

__global__ void kvec_add_double(double * __restrict,
				const double * __restrict,
				const double * __restrict);




#endif /*__GMS_INIT_GPUS_CUH__*/
