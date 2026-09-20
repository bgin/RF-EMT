

#include <immintrin.h>
#include <stdio.h>
#include "GMS_init_gpus.cuh"
//
//	Implementation
//
// Few of functions adapted from NVIDIA Samples.
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

#define TO_KiB(value) ((double)(value)) / (1024.0)
#ifndef MIN
#define MIN(x, y) ((x < y) ? x : y)
#endif

void devcaps_to_screen(const int dev_num, int * ierr, const bool all_info) {



#if (LOG_ACTIVITY) == 1
	GMS_CUDA_LOG("At prolog of devcaps_to_screen");
#endif
	cudaDeviceProp devp;
	cudaError status;
	
	status = cudaGetDeviceProperties(&devp,dev_num);
	if(status != cudaSuccess) {
		fprintf(stderr,"cudaGetDeviceProperties failed with an error: %s\n",cudaGetErrorString(status));
		*ierr = -1;
		return;
	}
	// Check CUDA device no. x is installed
	if(devp.major == 9999 && devp.minor == 9999) {
		fprintf(stderr,"cudaGetDeviceProperties find out no CUDA device(s) ... terminating!!\n");
		exit_fatal("No CUDA device has been found",
			      "CuWRF/Run-time errors.txt"  );
	}
	printf("Begin dump of CUDA device no: %d\n", dev_num);
	printf("Device ASCII name:			   %s       \n",devp.name);
	printf("Device major compute caps:    %d       \n",devp.major);
	printf("Device minor compute caps:    %d       \n",devp.minor);
	printf("\n");
	printf("CUDA Device HW characteristics: \n");
	printf("Device clock rate:             %d Khz  .\n",devp.clockRate);
	printf("Device async engine count:     %d       \n",devp.asyncEngineCount);
	printf("Device L2 cache size:	        %d bytes.\n",devp.l2CacheSize);
	printf("Device multi-processor count:  %d       \n",devp.multiProcessorCount);
	printf("Device concurrent kernels:     %d       \n",devp.concurrentKernels);
	printf("Device compute mode:           %d       \n",devp.computeMode);
	printf("Device kernel timeout enabled: %d       \n",devp.kernelExecTimeoutEnabled);
	printf("Device PCI Device ID:          %d       \n",devp.pciDeviceID);
	printf("Device PCI Domain ID:		    %d       \n",devp.pciDomainID);
	printf("Device PCI Bus ID:             %d       \n",devp.pciBusID);
	printf("Device ECC enabled:		    %d       \n",devp.ECCEnabled);
	printf("Device host mem mapping:       %d       \n",devp.canMapHostMemory);
	printf("Device registers per block:    %d       \n",devp.regsPerBlock);
	printf("Device warp size:              %d       \n",devp.warpSize);
	printf("Device overlap:                %d       \n",devp.deviceOverlap);
	printf("\n");
	printf("CUDA Device memory characteristics: \n");
	printf("Device memory bus width:       %d       \n",devp.memoryBusWidth);
	printf("Device memory clock rate:      %d       \n",devp.memoryClockRate);
	printf("Device memory pitch:           %d bytes.\n",devp.memPitch);
	printf("Device shared memory/block:    %d bytes.\n",devp.sharedMemPerBlock);
	printf("Device total const memory:     %.9f KiB.\n",TO_KiB(devp.totalConstMem));
	printf("Device total global memory:    %.9f KiB.\n",TO_KiB(devp.totalGlobalMem));
	printf("Device unified addressing:     %d       \n",devp.unifiedAddressing);
	printf("\n");
	printf("CUDA Device compute characteristics: \n");
	printf("Device max threads/block:      %d    \n",devp.maxThreadsPerBlock);
	printf("Device max threads/dimension:        \n");
	printf("Dim 1 = %d threads.\n",devp.maxThreadsDim[0]);
	printf("Dim 2 = %d threads.\n",devp.maxThreadsDim[1]);
	printf("Dim 3 = %d threads.\n",devp.maxThreadsDim[2]);
	printf("Device max thread/multiproc:   %d    \n",devp.maxThreadsPerMultiProcessor);
	printf("Device max size of grid:             \n");
	printf("Grid dim 1 = %d \n",devp.maxGridSize[0]);
	printf("Grid dim 2 = %d \n",devp.maxGridSize[1]);
	printf("Grid dim 3 = %d \n",devp.maxGridSize[2]);

	if(all_info == true) {

	   printf("\n");
	   printf("CUDA Device texturing characteristics: \n");
	   printf("Device max surface 1D:              %d     \n",devp.maxSurface1D);
	   printf("Device max surface layered 1D:      %d     \n",devp.maxSurface1DLayered);
	   printf("Device max surface 2D:              %d     \n",devp.maxSurface2D);
	   printf("Device max surface layered 2D:      %d     \n",devp.maxSurface2DLayered);
	   printf("Device max surface 3D:              %d     \n",devp.maxSurface3D);
	   printf("Device max surface Cubemap:         %d     \n",devp.maxSurfaceCubemap);
	   printf("Device max surface Cubemap layered: %d     \n",devp.maxSurfaceCubemapLayered);
	   printf("Device max texture 1D:			   %d     \n",devp.maxTexture1D);
	   printf("Device max texture 1D layered:      %d     \n",devp.maxTexture1DLayered);
	   printf("Device max texture 1D linear:       %d     \n",devp.maxTexture1DLinear);
	   printf("Device max texture 1D mipmap:       %d     \n",devp.maxTexture1DMipmap);
	   printf("Device max texture 2D:              %d     \n",devp.maxTexture2D);
	   printf("Device max texture 2D Gather:       %d     \n",devp.maxTexture2DGather);
	   printf("Device max texture 2D Linear:       %d     \n",devp.maxTexture2DLinear);
	   printf("Device max texture 2D Layered:      %d     \n",devp.maxTexture2DLayered);
	   printf("Device max texture 2D Gather:       %d     \n",devp.maxTexture2DGather);
	   printf("Device max texture 2D Mipmap:       %d     \n",devp.maxTexture2DMipmap);
	   printf("Device max texture 3D:			   %d     \n",devp.maxTexture3D);
	   printf("Device max texture 3D Alt:          %d     \n",devp.maxTexture3DAlt);
	}
	printf("CUDA Device name: %s -- end of capabilities dump!! \n",devp.name);
	*ierr = 0;
}

void devcaps_to_file(const int dev_num, const char * fname, 
		     int * ierr,     const bool all_info) {



#if (LOG_ACTIVITY) == 1
	GMS_CUDA_LOG("At prolog of devcaps_to_file");
#endif
	cudaDeviceProp devp;
	cudaError status;
#if (GMS_CUDA_DEBUG_ON) == 1
	 GMS_CUDA_DEBUG_CHECK(cudaGetDeviceProperties(&devp,dev_num));
#else
	 GMS_CUDA_CHECK(cudaGetDeviceProperties(&devp,dev_num));
#endif
	// Check CUDA device no. x is installed
	if(devp.major == 9999 && devp.minor == 9999) {
		fprintf(stderr,"cudaGetDeviceProperties find out no CUDA device(s) ... terminating!!\n");
		exit_fatal("No CUDA device has been found",
			      "CuWRF/Run-time errors.txt"  );
	}
	FILE *fp = NULL;
	fp = fopen(fname,"a+");
	if(NULL != fp) {
	   fprintf(fp,"Begin dump of CUDA device no: %d\n", dev_num);
	   fprintf(fp,"Device ASCII name:			   %s       \n",devp.name);
	   fprintf(fp,"Device major compute caps:      %d       \n",devp.major);
	   fprintf(fp,"Device minor compute caps:      %d       \n",devp.minor);
	   fprintf(fp,"\n");
	   fprintf(fp,"CUDA Device HW characteristics: \n");
	   fprintf(fp,"Device clock rate:             %d Khz  .\n",devp.clockRate);
	   fprintf(fp,"Device async engine count:     %d       \n",devp.asyncEngineCount);
	   fprintf(fp,"Device L2 cache size:	       %d bytes.\n",devp.l2CacheSize);
	   fprintf(fp,"Device multi-processor count:  %d       \n",devp.multiProcessorCount);
	   fprintf(fp,"Device concurrent kernels:     %d       \n",devp.concurrentKernels);
	   fprintf(fp,"Device compute mode:           %d       \n",devp.computeMode);
	   fprintf(fp,"Device kernel timeout enabled: %d       \n",devp.kernelExecTimeoutEnabled);
	   fprintf(fp,"Device PCI Device ID:          %d       \n",devp.pciDeviceID);
	   fprintf(fp,"Device PCI Domain ID:		  %d       \n",devp.pciDomainID);
	   fprintf(fp,"Device PCI Bus ID:             %d       \n",devp.pciBusID);
	   fprintf(fp,"Device ECC enabled:		       %d       \n",devp.ECCEnabled);
	   fprintf(fp,"Device host mem mapping:       %d       \n",devp.canMapHostMemory);
	   fprintf(fp,"Device registers per block:    %d       \n",devp.regsPerBlock);
	   fprintf(fp,"Device warp size:              %d       \n",devp.warpSize);
	   fprintf(fp,"Device overlap:                %d       \n",devp.deviceOverlap);
	   fprintf(fp,"\n");
	   fprintf(fp,"CUDA Device memory characteristics: \n");
	   fprintf(fp,"Device memory bus width:       %d       \n",devp.memoryBusWidth);
	   fprintf(fp,"Device memory clock rate:      %d       \n",devp.memoryClockRate);
	   fprintf(fp,"Device memory pitch:           %d bytes.\n",devp.memPitch);
	   fprintf(fp,"Device shared memory/block:    %d bytes.\n",devp.sharedMemPerBlock);
	   fprintf(fp,"Device total const memory:     %.9f KiB.\n",TO_KiB(devp.totalConstMem));
	   fprintf(fp,"Device total global memory:    %.9f KiB.\n",TO_KiB(devp.totalGlobalMem));
	   fprintf(fp,"Device unified addressing:     %d       \n",devp.unifiedAddressing);
	   fprintf(fp,"\n");
	   fprintf(fp,"CUDA Device compute characteristics: \n");
	   fprintf(fp,"Device max threads/block:      %d    \n",devp.maxThreadsPerBlock);
	   fprintf(fp,"Device max threads/dimension:        \n");
	   fprintf(fp,"Dim 1 = %d threads.\n",devp.maxThreadsDim[0]);
	   fprintf(fp,"Dim 2 = %d threads.\n",devp.maxThreadsDim[1]);
	   fprintf(fp,"Dim 3 = %d threads.\n",devp.maxThreadsDim[2]);
	   fprintf(fp,"Device max thread/multiproc:   %d    \n",devp.maxThreadsPerMultiProcessor);
	   fprintf(fp,"Device max size of grid:             \n");
	   fprintf(fp,"Grid dim 1 = %d \n",devp.maxGridSize[0]);
	   fprintf(fp,"Grid dim 2 = %d \n",devp.maxGridSize[1]);
	   fprintf(fp,"Grid dim 3 = %d \n",devp.maxGridSize[2]);

	   if(all_info == true) {

	      fprintf(fp,"\n");
	      fprintf(fp,"CUDA Device texturing characteristics: \n");
	      fprintf(fp,"Device max surface 1D:              %d     \n",devp.maxSurface1D);
	      fprintf(fp,"Device max surface layered 1D:      %d     \n",devp.maxSurface1DLayered);
	      fprintf(fp,"Device max surface 2D:              %d     \n",devp.maxSurface2D);
	      fprintf(fp,"Device max surface layered 2D:      %d     \n",devp.maxSurface2DLayered);
	      fprintf(fp,"Device max surface 3D:              %d     \n",devp.maxSurface3D);
	      fprintf(fp,"Device max surface Cubemap:         %d     \n",devp.maxSurfaceCubemap);
	      fprintf(fp,"Device max surface Cubemap layered: %d     \n",devp.maxSurfaceCubemapLayered);
	      fprintf(fp,"Device max texture 1D:			  %d     \n",devp.maxTexture1D);
	      fprintf(fp,"Device max texture 1D layered:      %d     \n",devp.maxTexture1DLayered);
	      fprintf(fp,"Device max texture 1D linear:       %d     \n",devp.maxTexture1DLinear);
	      fprintf(fp,"Device max texture 1D mipmap:       %d     \n",devp.maxTexture1DMipmap);
	      fprintf(fp,"Device max texture 2D:              %d     \n",devp.maxTexture2D);
	      fprintf(fp,"Device max texture 2D Gather:       %d     \n",devp.maxTexture2DGather);
	      fprintf(fp,"Device max texture 2D Linear:       %d     \n",devp.maxTexture2DLinear);
	      fprintf(fp,"Device max texture 2D Layered:      %d     \n",devp.maxTexture2DLayered);
	      fprintf(fp,"Device max texture 2D Gather:       %d     \n",devp.maxTexture2DGather);
	      fprintf(fp,"Device max texture 2D Mipmap:       %d     \n",devp.maxTexture2DMipmap);
	      fprintf(fp,"Device max texture 3D:			  %d     \n",devp.maxTexture3D);
	      fprintf(fp,"Device max texture 3D Alt:          %d     \n",devp.maxTexture3DAlt);
	   }
	   fprintf(fp,"CUDA Device name: %s -- end of capabilities dump!! \n",devp.name);
	   *ierr = 0;
	}
	else {
		fprintf(stderr,"fopen returned invalid pointer: %p \n",fp);
		*ierr = -2;
		return;
	}
Error:
	 *ierr = -1;
	 return;
}

// Adapted from NVIDIA Samples


int32_t convertSMVer2Cores(const int32_t major,
                           const int32_t minor) {
#if (LOG_ACTIVITY) == 1
	GMS_CUDA_LOG("At prolog of convertSMVer2Cores.");
#endif
     typedef struct {
         int32_t SM;  // 0xMm (hexidecimal notation), M = SM Major version,
                  // and m = SM minor version
         int32_t Cores;
  } sSMtoCores;

  sSMtoCores nGpuArchCoresPerSM[] = {
      {0x30, 192},
      {0x32, 192},
      {0x35, 192},
      {0x37, 192},
      {0x50, 128},
      {0x52, 128},
      {0x53, 128},
      {0x60,  64},
      {0x61, 128},
      {0x62, 128},
      {0x70,  64},
      {0x72,  64},
      {0x75,  64},
      {0x80,  64},
      {0x86, 128},
      {0x87, 128},
      {-1, -1}};

    int32_t index = 0;

    while (nGpuArchCoresPerSM[index].SM != -1) {
      if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor)) {
          return nGpuArchCoresPerSM[index].Cores;
     }

     index++;
   }

  
     printf(
      "MapSMtoCores for SM %d.%d is undefined."
      "  Default to use %d Cores/SM\n",
        major, minor, nGpuArchCoresPerSM[index - 1].Cores);
     return (nGpuArchCoresPerSM[index - 1].Cores);
}

const char * convertSMVer2ArchName(const int32_t major,
                                   const int32_t minor) {
#if (LOG_ACTIVITY) == 1
	GMS_CUDA_LOG("At prolog of convertSMVer2ArchName.");
#endif
     typedef struct {
         int32_t SM;  // 0xMm (hexidecimal notation), M = SM Major version,
                  // and m = SM minor version
         const char* name;
     } sSMtoArchName;

    sSMtoArchName nGpuArchNameSM[] = {
      {0x30, "Kepler"},
      {0x32, "Kepler"},
      {0x35, "Kepler"},
      {0x37, "Kepler"},
      {0x50, "Maxwell"},
      {0x52, "Maxwell"},
      {0x53, "Maxwell"},
      {0x60, "Pascal"},
      {0x61, "Pascal"},
      {0x62, "Pascal"},
      {0x70, "Volta"},
      {0x72, "Xavier"},
      {0x75, "Turing"},
      {0x80, "Ampere"},
      {0x86, "Ampere"},
      {-1, "Graphics Device"}};

   int32_t index = 0;

  while (nGpuArchNameSM[index].SM != -1) {
    if (nGpuArchNameSM[index].SM == ((major << 4) + minor)) {
        return nGpuArchNameSM[index].name;
     }

    index++;
  }

  
  printf(
      "MapSMtoArchName for SM %d.%d is undefined."
      "  Default to use %s\n",
      major, minor, nGpuArchNameSM[index - 1].name);
     return (nGpuArchNameSM[index - 1].name);
}


 
int32_t gpuDeviceInit(const int32_t devID) {
  int32_t device_count;
  

    GMS_CUDA_DEBUG_CHECK(cudaGetDeviceCount(&device_count));
     
  if (device_count == 0) {
    fprintf(stderr,
            "gpuDeviceInit() CUDA error: "
            "no devices supporting CUDA.\n");
    exit(EXIT_FAILURE);
  }

  if (devID < 0) {
    devID = 0;
  }

  if (devID > device_count - 1) {
    fprintf(stderr, "\n");
    fprintf(stderr, ">> %d CUDA capable GPU device(s) detected. <<\n",
            device_count);
    fprintf(stderr,
            ">> gpuDeviceInit (-device=%d) is not a valid"
            " GPU device. <<\n",
            devID);
    fprintf(stderr, "\n");
    return -devID;
  }

  int computeMode = -1, major = 0, minor = 0;
  
   GMS_CUDA_DEBUG_CHECK(cudaDeviceGetAttribute(&computeMode, cudaDevAttrComputeMode, devID));
   GMS_CUDA_DEBUG_CHECK(cudaDeviceGetAttribute(&major, cudaDevAttrComputeCapabilityMajor, devID));
   GMS_CUDA_DEBUG_CHECK(cudaDeviceGetAttribute(&minor, cudaDevAttrComputeCapabilityMinor, devID));
  if (computeMode == cudaComputeModeProhibited) {
    fprintf(stderr,
            "Error: device is running in <Compute Mode "
            "Prohibited>, no threads can use cudaSetDevice().\n");
    return -1;
  }

  if (major < 1) {
    fprintf(stderr, "gpuDeviceInit(): GPU device does not support CUDA.\n");
    exit(EXIT_FAILURE);
  }

  GMS_CUDA_DEBUG_CHECK(cudaSetDevice(devID));
  printf("gpuDeviceInit() CUDA Device [%d]: \"%s\n", devID, _ConvertSMVer2ArchName(major, minor));

  return devID;  
  Error:
        return (-1);
}

int32_t gpuGetMaxGflopsDeviceId() {
      
  int32_t current_device = 0, sm_per_multiproc = 0;
  int32_t max_perf_device = 0;
  int32_t device_count = 0;
  int32_t devices_prohibited = 0;

  uint64_t max_compute_perf = 0;
  GMS_CUDA_DEBUG_CHECK(cudaGetDeviceCount(&device_count));

  if (device_count == 0) {
    fprintf(stderr,
            "gpuGetMaxGflopsDeviceId() CUDA error:"
            " no devices supporting CUDA.\n");
    exit(EXIT_FAILURE);
  }

  // Find the best CUDA capable GPU device
  current_device = 0;

  while (current_device < device_count) {
    int computeMode = -1, major = 0, minor = 0;
    GMS_CUDA_DEBUG_CHECK(cudaDeviceGetAttribute(&computeMode, cudaDevAttrComputeMode, current_device));
    GMS_CUDA_DEBUG_CHECK(cudaDeviceGetAttribute(&major, cudaDevAttrComputeCapabilityMajor, current_device));
    GMS_CUDA_DEBUG_CHECK(cudaDeviceGetAttribute(&minor, cudaDevAttrComputeCapabilityMinor, current_device));

    // If this GPU is not running on Compute Mode prohibited,
    // then we can add it to the list
    if (computeMode != cudaComputeModeProhibited) {
      if (major == 9999 && minor == 9999) {
        sm_per_multiproc = 1;
      } else {
        sm_per_multiproc =
            _ConvertSMVer2Cores(major,  minor);
      }
      int multiProcessorCount = 0, clockRate = 0;
      GMS_CUDA_DEBUG_CHECK(cudaDeviceGetAttribute(&multiProcessorCount, cudaDevAttrMultiProcessorCount, current_device));
      cudaError_t result = cudaDeviceGetAttribute(&clockRate, cudaDevAttrClockRate, current_device);
      if (result != cudaSuccess) {
        // If cudaDevAttrClockRate attribute is not supported we
        // set clockRate as 1, to consider GPU with most SMs and CUDA Cores.
        if(result == cudaErrorInvalidValue) {
          clockRate = 1;
        }
        else {
          fprintf(stderr, "CUDA error at %s:%d code=%d(%s) \n", __FILE__, __LINE__,
            static_cast<unsigned int>(result), _cudaGetErrorEnum(result));
          exit(EXIT_FAILURE);
        }
      }
      uint64_t compute_perf = (uint64_t)multiProcessorCount * sm_per_multiproc * clockRate;

      if (compute_perf > max_compute_perf) {
        max_compute_perf = compute_perf;
        max_perf_device = current_device;
      }
    } else {
      devices_prohibited++;
    }

    ++current_device;
  }

  if (devices_prohibited == device_count) {
    fprintf(stderr,
            "gpuGetMaxGflopsDeviceId() CUDA error:"
            " all devices have compute mode prohibited.\n");
    exit(EXIT_FAILURE);
  }

  return max_perf_device;   
  Error:
         return (-1); 
}

int32_t findIntegratedGPU() {
    
  int32_t current_device = 0;
  int32_t device_count = 0;
  int32_t devices_prohibited = 0;

  GMS_CUDA_DEBUG_CHECK(cudaGetDeviceCount(&device_count));

  if (device_count == 0) {
    fprintf(stderr, "CUDA error: no devices supporting CUDA.\n");
    exit(EXIT_FAILURE);
  }

 
  while (current_device < device_count) {
    int computeMode = -1, integrated = -1;
    GMS_CUDA_DEBUG_CHECK(cudaDeviceGetAttribute(&computeMode, cudaDevAttrComputeMode, current_device));
    GMS_CUDA_DEBUG_CHECK(cudaDeviceGetAttribute(&integrated, cudaDevAttrIntegrated, current_device));
    
    if (integrated && (computeMode != cudaComputeModeProhibited)) {
        GMS_CUDA_DEBUG_CHECK(cudaSetDevice(current_device));

      int major = 0, minor = 0;
      GMS_CUDA_DEBUG_CHECK(cudaDeviceGetAttribute(&major, cudaDevAttrComputeCapabilityMajor, current_device));
      GMS_CUDA_DEBUG_CHECK(cudaDeviceGetAttribute(&minor, cudaDevAttrComputeCapabilityMinor, current_device));
      printf("GPU Device %d: \"%s\" with compute capability %d.%d\n\n",
             current_device, convertSMVer2ArchName(major, minor), major, minor);

      return current_device;
    } else {
      devices_prohibited++;
    }

    current_device++;
  }

  if (devices_prohibited == device_count) {
    fprintf(stderr,
            "CUDA error:"
            " No GLES-CUDA Interop capable GPU found.\n");
    exit(EXIT_FAILURE);
  }

   return 0;
   Error:
          return (-1);
}

bool checkCudaCapabilities(const int32_t major,
                           const int32_t minor) {
  int32_t dev;
  int32_t major = 0, minor = 0;

  GMS_CUDA_DEBUG_CHECK(cudaGetDevice(&dev));
  GMS_CUDA_DEBUG_CHECK(cudaDeviceGetAttribute(&major, cudaDevAttrComputeCapabilityMajor, dev));
  gms_cuda_debug_check(cudaDeviceGetAttribute(&minor, cudaDevAttrComputeCapabilityMinor, dev));

  if ((major > major_version) ||
      (major == major_version &&
       minor >= minor_version)) {
    printf("  Device %d: <%16s >, Compute SM %d.%d detected\n", dev,
           _ConvertSMVer2ArchName(major, minor), major, minor);
    return true;
  } else {
    printf(
        "  No GPU device was found that can support "
        "CUDA compute capability %d.%d.\n",
        major_version, minor_version);
    return false;
    }
    Error:
           return (-1);
}

 static uint32_t npow2(uint32_t x) {
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return ++x;
}

void setBlocksThreadNumber(const int32_t arch,
                           const int32_t n,
                           const int32_t mblocks,
                           const int32_t mthreads,
                           int32_t & blocks,
                           int32_t & threads) {
 
      cudaDeviceProp prop;
      int32_t device;
      GMS_CUDA_DEBUG_CHECK(cudaGetDevice(&device));
      GMS_CUDA_DEBUG_CHECK(cudaGetDeviceProperties(&prop, device));

     if (arch < 3) {
         threads = (n < mthreads) ? npow2(n) : mthreads;
         blocks = (n + threads - 1) / threads;
     } else {
        threads = (n < mthreads * 2) ? npow2((n + 1) / 2) : mthreads;
        blocks = (n + (threads * 2 - 1)) / (threads * 2);
  }

  if ((float)threads * blocks >
      (float)prop.maxGridSize[0] * prop.maxThreadsPerBlock) {
      printf("n is too large!\n");
  }

  if (blocks > prop.maxGridSize[0]) {
    printf(
        "Grid size <%d> exceeds the device capability <%d>, set block size as "
        "%d (original %d)\n",
        blocks, prop.maxGridSize[0], threads * 2, threads);

    blocks /= 2;
    threads *= 2;
  }

  if (arch >= 6) {
    blocks = MIN(mblocks, blocks);
  } 
  return;
  Error:
           blocks  = -1;
           threads = -1;
}

#if !defined (CHECK_FAILURE_GENERATOR)
#define CHECK_FAILURE_GENERATOR(ierr,msg) \
 do {                                  \
     if((ierr) < 0) {                     \
         fprintf(stderr,"%s %d\n",(msg),(ierr));      \
         return;                        \
	 }                                  \
 }while(0);
#endif

#if !defined (MALLOC_FAILED_GENERATOR)
#define MALLOC_FAILED_GENERATOR(ptr,msg) \
	do {                              \
		 if(NULL == ptr) {             \
		    fprintf(stderr, "%s %p\n",(msg),(ptr)); \
			return;								 \
		 }										 \
	}while(0);
#endif

#if !defined (INIT_VECI32_GENERATOR) 
#define INIT_VECI32_GENERATOR     \
	for(int32_t i = 0; i != len; ++i) {  \
        h_i32a[i] = i;                \
        h_i32b[i] = i;                \
		h_i32c[i] = 0;                \
	}                                  
#endif

#if !defined (INIT_VECR4_GENERATOR)
#define INIT_VECR4_GENEREATOR     \
	for(int32_t i = 0; i != len; ++i) { \
		h_r4a[i] = (float)i;       \
		h_r4b[i] = (float)i;       \
		h_r4c[i] = 0.F;            \
	}
#endif

#if !defined (INIT_VECR8_GENERATOR)
#define INIT_VECR8_GENERATOR      \
	for(int32_t i = 0; i != len; ++i) { \
		h_r8a[i] = (double)i;       \
		h_r8b[i] = (double)i;       \
		h_r8c[i] = 0.0;            \
	}
#endif

// Random value initialization

#if !defined (INIT_RAND_VECI4_GENERATOR)
#define INIT_RAND_VECI4_GENERATOR             \
	for(int32_t i = 0; i != len; ++i) {       \
	   uint32_t seed = 0;			           \
	   int32_t ret = _rdrand32_step(&seed);   \
	   if(ret == 1){                        \
	      h_i32a[i] = seed;                  \
		  h_i32b[i] = seed;                  \
	   }                                    \
	   else{                                \
	        h_i32a[i] = i;				    \
			h_i32b[i] = i;				    \
	   }								        \
	   h_i32c[i] = 0;                        \
    }
#endif
		
#if !defined (INIT_RAND_VECR4_GENERATOR)
#define INIT_RAND_VECR4_GENERATOR            \
	for(int32_t i = 0; i != len; ++i) {      \
		uint32_t seed = 0;                  \
		int32_t ret = _rdrand32_step(&seed);  \
		if(1 == ret) {                      \
		   h_r4a[i] = (float)seed;           \
		   h_r4b[i] = (float)seed;           \
		}                                  \
		else {                             \
		    h_r4a[i] = (float)i;            \
            h_r4b[i] = (float)i;            \
		}                                 \
		h_r4c[i] = 0.F;                    \
	}
#endif

#if !defined (INIT_RAND_VECR8_GENERATOR)
#define INIT_RAND_VECR8_GENERATOR            \
	for(int32_t i = 0; i != len; ++i) {      \
	   uint64_t seed = 0ULL;                \
	   int32_t ret = _rdrand64_step(&seed);  \
	   if(1 == ret) {                      \
          h_r8a[i] = (double)seed;          \
		  h_r8b[i] = (double)seed;          \
	   }									  \
	   else {                             \
	         h_r8a[i] = (double)i;         \
			 h_r8b[i] = (double)i;         \
	   }                                 \
	    h_r8c[i] = 0.0;                   \
	}
#endif

#if !defined (FREE_GPU_MEMORY) 
#define FREE_GPU_MEMORY               \
	 cudaFree(dev_r8c);			    \
	 cudaFree(dev_r8b);				\
	 cudaFree(dev_r8a);				\
	 cudaFree(dev_r4c);		        \
	 cudaFree(dev_r4b);		        \
	 cudaFree(dev_r4a);			    \
	 cudaFree(dev_i32c);              \
	 cudaFree(dev_i32b);              \
	 cudaFree(dev_i32a);              
#endif

#if !defined (FREE_CPU_MEMORY)
#define FREE_CPU_MEMORY               \
     _aligned_free(h_r8c);            \
	 _mm_free(h_r8b);			\
	 _mm_free(h_r8a);		    \
	 _mm_free(h_r4c);			\
	 _mm_free(h_r4b);			\
	 _mm_free(h_r4a);		    \
	 _mm_free(h_i32c);		    \
	 _mm_free(h_i32b);           \
	 _mm_free(h_i32a);			
#endif

#if !defined (PRINT_VECI4)
#define PRINT_VECI4	  \
	for(int32_t i = 0; i != len; ++i)  \
	    printf("GPU computed int32_t vector: %d\n",h_i32c[i]);
#endif

#if !defined (PRINT_VECR4)
#define PRINT_VECR4     \
	for(int32_t i = 0; i != len; ++i) \
		 printf("GPU computed REAL4 vector: %.6f\n",h_r4c[i]);
#endif

#if !defined (PRINT_VECR8)
#define PRINT_VECR8    \
	for(int32_t i = 0; i != len; ++i) \
		printf("GPU computed REAL8 vector: %.16f\n",h_r8c[i]);
#endif






__global__ void kvec_add_int32(int32_t * __restrict c,
			       const int32_t * __restrict b,
			       const int32_t * __restrict a) {
	int i = threadIdx.x;
	c[i] = b[i] + a[i];
}

__global__ void kvec_add_real4(float * __restrict c,
			       const float * __restrict b,
			       const float * __restrict a) {
	int i = threadIdx.x;
	c[i] = b[i] + a[i];
}

__global__ void kvec_add_real8(double * __restrict c,
			        const double * __restrict b,
				const double * __restrict a) {
	int i = threadIdx.x;
	c[i] = b[i] + a[i];
}

cudaError_t gpu_vec_add_tests(const int dev_num) {


#if (LOG_ACTIVITY) == 1
	GMS_CUDA_LOG("At prolog of gpu_vec_add_tests.\n");
#endif
	cudaError status;
	int32_t * __restrict dev_i32a = NULL;
	int32_t * __restrict dev_i32b = NULL;
	int32_t * __restrict dev_i32c = NULL;
	float   * __restrict dev_r4a  = NULL;
	float   * __restrict dev_r4b  = NULL;
	float   * __restrict dev_r4c  = NULL;
	double  * __restrict dev_r8a  = NULL;
	double  * __restrict dev_r8b  = NULL;
	double  * __restrict dev_r8c  = NULL;
	// Host arrays.
	int32_t * __restrict h_i32a = NULL;
	int32_t * __restrict h_i32b = NULL;
	int32_t * __restrict h_i32c = NULL;
	float *   __restrict h_r4a  = NULL;
	float *   __restrict h_r4b  = NULL;
	float *   __restrict h_r4c  = NULL;
	double *   __restrict h_r8a  = NULL;
	double *   __restrict h_r8b  = NULL;
	double *   __restrict h_r8c  = NULL;
	int32_t  ierr = 0;
	const int32_t len = 1 << 16;

	status = cudaSetDevice(dev_num);
	if(status != cudaSuccess) {
		fprintf(stderr,"CUDA Device: -- FATAL-ERROR -- : [%s]\n",cudaGetErrorString(status));
		return;
	}
	
	alloc_int32_gpu(&dev_i32a[0],len,&ierr);
	CHECK_FAILURE_GENERATOR(ierr,"alloc_int32_gpu failed allocating: dev_i32a")

	alloc_int32_gpu(&dev_i32b[0],len,&ierr);
	CHECK_FAILURE_GENERATOR(ierr,"alloc_int32_gpu failed allocating: dev_i32b")

	alloc_int32_gpu(&dev_i32c[0],len,&ierr);
	CHECK_FAILURE_GENERATOR(ierr,"alloc_int32_gpu failed allocating: dev_i32c")

	alloc_float_gpu(&dev_r4a[0],len,&ierr);
	CHECK_FAILURE_GENERATOR(ierr,"alloc_float_gpu failed allocating: dev_r4a")

	alloc_float_gpu(&dev_r4b[0],len,&ierr);
	CHECK_FAILURE_GENERATOR(ierr,"alloc_float_gpu failed allocating: dev_r4b")

	alloc_float_gpu(&dev_r4c[0],len,&ierr);
	CHECK_FAILURE_GENERATOR(ierr,"alloc_float_gpu failed allocating: dev_r4c")

	alloc_double_gpu(&dev_r8a[0],len,&ierr);
	CHECK_FAILURE_GENERATOR(ierr,"alloc_double_gpu failed allocating: dev_r8a")

	alloc_double_gpu(&dev_r8b[0],len,&ierr);
	CHECK_FAILURE_GENERATOR(ierr,"alloc_double_gpu failed allocating: dev_r8b")

	alloc_double_gpu(&dev_r8c[0],len,&ierr);
	CHECK_FAILURE_GENERATOR(ierr,"alloc_double_gpu failed allocating: dev_r8c")

    // Allocate host arrays.
	h_i32a = (int32_t*)_mm_malloc(((size_t)len) * sizeof(int32_t),HOST_ALIGN64);
	MALLOC_FAILED_GENERATOR(h_i32a,"_mm_malloc failed to allocate: h_i32a")

	h_i32b = (int32_t*)_mm_malloc(((size_t)len) * sizeof(int32_t),HOST_ALIGN64);
	MALLOC_FAILED_GENERATOR(h_i32b,"_aligned_malloc failed to allocate: h_i32b")

	h_i32c = (int32_t*)_mm_malloc(((size_t)len) * sizeof(int32_t),HOST_ALIGN64);
	MALLOC_FAILED_GENERATOR(h_i32c,"_aligned_malloc failed to allocate: h_i32c")

	h_r4a = (float*)_mm_malloc(((size_t)len) * sizeof(float),HOST_ALIGN64);
	MALLOC_FAILED_GENERATOR(h_r4a,"_aligned_malloc failed to allocate: h_r4a")

	h_r4b = (float*)_mm_malloc(((size_t)len) *sizeof(float),HOST_ALIGN64);
	MALLOC_FAILED_GENERATOR(h_r4b,"_aligned_malloc failed to allocate: h_r4b")

	h_r4c = (float*)_mm_malloc(((size_t)len) * sizeof(float),HOST_ALIGN64);
	MALLOC_FAILED_GENERATOR(h_r4c,"_aligned_malloc failed to allocate: h_r4c")

	h_r8a = (double*)_mm_malloc(((size_t)len) * sizeof(double),HOST_ALIGN64);
	MALLOC_FAILED_GENERATOR(h_r8a,"_aligned_malloc failed to allocate: h_r8a")

	h_r8b = (double*)_mm_malloc(((size_t)len) * sizeof(double),HOST_ALIGN64);
	MALLOC_FAILED_GENERATOR(h_r8b,"_aligned_malloc failed to allocate: h_r8b")

	h_r8c = (double*)_mm_malloc(((size_t)len) * sizeof(double),HOST_ALIGN64);
	MALLOC_FAILED_GENERATOR(h_r8c,"_aligned_malloc failed to allocate: h_r8c")

	INIT_RAND_VECI4_GENERATOR
	INIT_RAND_VECR4_GENERATOR
	INIT_RAND_VECR8_GENERATOR

	// Copy host arrays to device arrays
	// Do not check for minor errors (non-critical)
	copy_int32_cpu_to_gpu(&dev_i32a[0],&h_i32a[0],len,&ierr);
	copy_int32_cpu_to_gpu(&dev_i32b[0],&h_i32b[0],len,&ierr);
	copy_int32_cpu_to_gpu(&dev_i32c[0],&h_i32c[0],len,&ierr);
	copy_float_cpu_to_gpu(&dev_r4a[0],&h_r4a[0],len,&ierr);
	copy_float_cpu_to_gpu(&dev_r4b[0],&h_r4b[0],len,&ierr);
	copy_float_cpu_to_gpu(&dev_r4c[0],&h_r4c[0],len,&ierr);
	copy_double_cpu_to_gpu(&dev_r8a[0],&h_r8a[0],len,&ierr);
	copy_double_cpu_to_gpu(&dev_r8b[0],&h_r8b[0],len,&ierr);
	copy_double_cpu_to_gpu(&dev_r8c[0],&h_r8c[0],len,&ierr);

	//
	// Launch simple kernels with single thread
	// per array element.
	//

	kvec_add_int32<<<1,len>>>(&dev_i32c[0],&dev_i32b[0],&dev_i32a[0]);
	kvec_add_real4<<<1,len>>>(&dev_r4c[0],&dev_r4b[0],&dev_r4a[0]);
	kvec_add_real8<<<1,len>>>(&dev_r8c[0],&dev_r8b[0],&dev_r8a[0]);
	
	GMS_CUDA_DEBUG_CHECK(cudaGetLastError());
        GMS_CUDA_DEBUG_CHECK(cudaDeviceSynchronize());

	//
	// Copy rresulting device arrays to host memory.
	//Do not check for minor errors (non-critical)
	//
	copy_int32_gpu_to_cpu(&dev_i32c[0],&h_i32c[0],len,&ierr);
	copy_float_gpu_to_cpu(&dev_r4c[0],&h_r4c[0],len,&ierr);
	copy_double_gpu_to_cpu(&dev_r8c[0],&h_r8c[0],len,&ierr);

	

	PRINT_VECI4
	PRINT_VECR4
	PRINT_VECR8

	GMS_CUDA_DEBUG_CHECK(cudaDeviceReset());
	
	FREE_CPU_MEMORY
Error:
	
	FREE_GPU_MEMORY
	FREE_CPU_MEMORY
	return;
}
