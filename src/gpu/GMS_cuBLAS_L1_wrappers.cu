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
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include "GMS_cuBLAS_L1_wrappers.cuh"
#include "GMS_gpu_config.cuh"
#if (PROFILE_HOST_TO_DEVICE) == 1
#include <immintrin.h> //rdtscp
#endif 
#include "GMS_cuda_memops.cuh"

static const uint64_t rdtscp_cost = 42; // Skylake uarch

void
cuBLAS_Isamax_iface(const float * __restrict h_ptr,
                    const int32_t n,
                    const int32_t incx,
                    int32_t * result,
                    cudaError_t * __restrict err,
                    int32_t * __restrict ierr
                    uint64_t * __restrict tsc_delta) { // for profiling usage 
                                  
    if(__builtin_expect(n<=1,0)) {
       *ierr = -1;
       return;
      }
    *err    = 0;
    *ierr   = 0;
    *result = 0;
#if (PROFILE_HOST_TO_DEVICE) == 1
     volatile uint64_t dummy1;
     volatile uint32_t dummy2;
     volatile uint64_t tsc_start,tsc_end;
     volatile uint32_t coreid;
#endif
    float * __restrict d_ptr = NULL;
    cublasHandle_t handle    = NULL;
    cudaStream_t   stream    = NULL;
    cudaError_t status;
    int32_t merr = 0;
    GMS_CUBLAS_STAT_CHECK(cublasCreate(&handle));
    GMS_CUDA_DEBUG_CHECK(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
    GMS_CUBLAS_STAT_CHECK(cublasSetStream(handle,stream));
    alloc_float_gpu(&d_ptr[0],(size_t)n,&merr);
    if(merr != 0) goto Error;
    GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptr,h_ptr,sizeof(float)*((size_t)n),
                         cudaMemcpyHostToDevice,stream));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      dummy1    = __rdtscp(&dummy2);
      tsc_start = __rdtscp(&coreid);
      __asm__("lfence");
#endif
     GMS_CUBLAS_STAT_CHECK(cublasIsamax(handle,n,d_ptr,incx,&result));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      tsc_end     = __rdtscp(&coreid);
      *tsc_delta  = tsc_end-tsc_start-rdtscp_cost;
      __asm__("lfence");
      printf("cublasIsamax executed in: %llu reference cycles\n",*tsc_delta);
#endif
      GMS_CUDA_DEBUG_CHECK(cudaStreamSynchronize(stream));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptr));
      GMS_CUBLAS_STAT_CHECK(cublasDestroy(handle));
      GMS_CUDA_DEBUG_CHECK(cudaStreamDestroy(stream));
      *ierr = 0;
      *err  = status;
      return;
      Error: 
            {
                if(d_ptr) cudaFree(d_ptr);
                cublasDestroy(handle);
                cudaStreamDestroy(stream);  
                ierr = -2;
                *err = status;  
                return;  
      }
}


cuBLAS_Icamax_iface(const cuComplex * __restrict h_ptr,
                    const int32_t n,
                    const int32_t incx,
                    int32_t * result,
                    cudaError_t * __restrict err,
                    int32_t * __restrict ierr
                    uint64_t * __restrict tsc_delta) { // for profiling usage 
                                  
    if(__builtin_expect(n<=1,0)) {
       *ierr = -1;
       return;
      }
     *err    = 0;
     *ierr   = 0;
     *result = 0;
#if (PROFILE_HOST_TO_DEVICE) == 1
     volatile uint64_t dummy1;
     volatile uint32_t dummy2;
     volatile uint64_t tsc_start,tsc_end;
     volatile uint32_t coreid;
#endif
    cuComplex * __restrict d_ptr = NULL;
    cublasHandle_t handle    = NULL;
    cudaStream_t   stream    = NULL;
    cudaError_t status;
    int32_t merr = 0;
    GMS_CUBLAS_STAT_CHECK(cublasCreate(&handle));
    GMS_CUDA_DEBUG_CHECK(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
    GMS_CUBLAS_STAT_CHECK(cublasSetStream(handle,stream));
    alloc_complex4_gpu(&d_ptr[0],(size_t)n,&merr);
    if(merr != 0) goto Error;
    GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptr,h_ptr,sizeof(cuComplex)*((size_t)n),
                         cudaMemcpyHostToDevice,stream));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      dummy1    = __rdtscp(&dummy2);
      tsc_start = __rdtscp(&coreid);
      __asm__("lfence");
#endif
     GMS_CUBLAS_STAT_CHECK(cublasIcamax(handle,n,d_ptr,incx,&result));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      tsc_end     = __rdtscp(&coreid);
      *tsc_delta  = tsc_end-tsc_start-rdtscp_cost;
      __asm__("lfence");
      printf("cublasIcamax executed in: %llu reference cycles\n",*tsc_delta);
#endif
      GMS_CUDA_DEBUG_CHECK(cudaStreamSynchronize(stream));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptr));
      GMS_CUBLAS_STAT_CHECK(cublasDestroy(handle));
      GMS_CUDA_DEBUG_CHECK(cudaStreamDestroy(stream));
      *ierr = 0;
      *err  = status;
       return;
      Error: 
            {
                if(d_ptr) cudaFree(d_ptr);
                cublasDestroy(handle);
                cudaStreamDestroy(stream);  
                ierr = -2;
                *err = status;  
                return;  
      }
}


void
cuBLAS_Isamin_iface(const float * __restrict h_ptr,
                    const int32_t n,
                    const int32_t incx,
                    int32_t *     __restrict result,
                    cudaError_t * __restrict err,
                    int32_t * __restrict ierr,
                    uint64_t * __restrict tsc_delta) {
  
     if(__builtin_expect(n<=1,0)) {
       *ierr = -1;
       return;
      }
    *err    = 0;
    *ierr   = 0;
    *result = 0;
#if (PROFILE_HOST_TO_DEVICE) == 1
     volatile uint64_t dummy1;
     volatile uint32_t dummy2;
     volatile uint64_t tsc_start,tsc_end;
     volatile uint32_t coreid;
#endif
    float * __restrict d_ptr = NULL;
    cublasHandle_t handle    = NULL;
    cudaStream_t   stream    = NULL;
    cudaError_t_t status;
    int32_t merr = 0;
    GMS_CUBLAS_STAT_CHECK(cublasCreate(&handle));
    GMS_CUDA_DEBUG_CHECK(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
    GMS_CUBLAS_STAT_CHECK(cublasSetStream(handle,stream));
    alloc_float_gpu(&d_ptr[0],(size_t)n,&merr);
    if(merr != 0) goto Error;
    GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptr,h_ptr,sizeof(float)*((size_t)n),
                         cudaMemcpyHostToDevice,stream));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      dummy1    = __rdtscp(&dummy2);
      tsc_start = __rdtscp(&coreid);
      __asm__("lfence");
#endif
     GMS_CUBLAS_STAT_CHECK(cublasIsamin(handle,n,d_ptr,incx,&result));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      tsc_end     = __rdtscp(&coreid);
      *tsc_delta  = tsc_end-tsc_start-rdtscp_cost;
      __asm__("lfence");
      printf("cublasIsamin executed in: %llu reference cycles\n",*tsc_delta);
#endif
      GMS_CUDA_DEBUG_CHECK(cudaStreamSynchronize(stream));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptr));
      GMS_CUBLAS_STAT_CHECK(cublasDestroy(handle));
      GMS_CUDA_DEBUG_CHECK(cudaStreamDestroy(stream));
      *ierr = 0;
      *err  = status;
      return;
      Error: 
            {
                if(d_ptr) cudaFree(d_ptr);
                cublasDestroy(handle);
                cudaStreamDestroy(stream);  
                ierr = -2;
                *err = status; 
                return;   
      }
} 


void
cuBLAS_Icamin_iface(const cuComplex * __restrict h_ptr,
                    const int32_t n,
                    const int32_t incx,
                    int32_t *  __restrict result,
                    cudaError_t * __restrict err,
                    int32_t * __restrict ierr,
                    uint64_t * __restrict tsc_delta) {
      
      if(__builtin_expect(n<=1,0)) {
       *ierr = -1;
       return;
      }
    *err    = 0;
    *ierr   = 0;
    *result = 0;
#if (PROFILE_HOST_TO_DEVICE) == 1
     volatile uint64_t dummy1;
     volatile uint32_t dummy2;
     volatile uint64_t tsc_start,tsc_end;
     volatile uint32_t coreid;
#endif
    cuComplex * __restrict d_ptr = NULL;
    cublasHandle_t handle    = NULL;
    cudaStream_t   stream    = NULL;
    cudaError_t status;
    int32_t merr = 0;
    GMS_CUBLAS_STAT_CHECK(cublasCreate(&handle));
    GMS_CUDA_DEBUG_CHECK(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
    GMS_CUBLAS_STAT_CHECK(cublasSetStream(handle,stream));
    alloc_complex4_gpu(&d_ptr[0],(size_t)n,&merr);
    if(merr != 0) goto Error;
    GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptr,h_ptr,sizeof(cuComplex)*((size_t)n),
                         cudaMemcpyHostToDevice,stream));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      dummy1    = __rdtscp(&dummy2);
      tsc_start = __rdtscp(&coreid);
      __asm__("lfence");
#endif
     GMS_CUBLAS_STAT_CHECK(cublasIcamin(handle,n,d_ptr,incx,&result));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      tsc_end     = __rdtscp(&coreid);
      *tsc_delta  = tsc_end-tsc_start-rdtscp_cost;
      __asm__("lfence");
      printf("cublasIcamin executed in: %llu reference cycles\n",*tsc_delta);
#endif
      GMS_CUDA_DEBUG_CHECK(cudaStreamSynchronize(stream));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptr));
      GMS_CUBLAS_STAT_CHECK(cublasDestroy(handle));
      GMS_CUDA_DEBUG_CHECK(cudaStreamDestroy(stream));
      *ierr = 0;
      *err  = status;
       return;
      Error: 
            {
                if(d_ptr) cudaFree(d_ptr);
                cublasDestroy(handle);
                cudaStreamDestroy(stream);  
                ierr = -2;
                *err = status;   
                return; 
      }
}


void
cuBLAS_Sasum_iface(const float * __restrict h_ptr,
                   const int32_t n,
                   const int32_t incx,
                   float * __restrict result,
                   cudaError_t * __restrict err,
                   int32_t * __restrict ierr,
                   uint64_t * __restrict tsc_delta) {

     if(__builtin_expect(n<=1,0)) {
       *ierr = -1;
       return;
      }
    *err    = 0;
    *ierr   = 0;
    *result = 0;
#if (PROFILE_HOST_TO_DEVICE) == 1
     volatile uint64_t dummy1;
     volatile uint32_t dummy2;
     volatile uint64_t tsc_start,tsc_end;
     volatile uint32_t coreid;
#endif
    float * __restrict d_ptr = NULL;
    cublasHandle_t handle    = NULL;
    cudaStream_t   stream    = NULL;
    cudaError_t status;
    int32_t merr = 0;
    GMS_CUBLAS_STAT_CHECK(cublasCreate(&handle));
    GMS_CUDA_DEBUG_CHECK(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
    GMS_CUBLAS_STAT_CHECK(cublasSetStream(handle,stream));
    alloc_float_gpu(&d_ptr[0],(size_t)n,&merr);
    if(merr != 0) goto Error;
    GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptr,h_ptr,sizeof(float)*((size_t)n),
                         cudaMemcpyHostToDevice,stream));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      dummy1    = __rdtscp(&dummy2);
      tsc_start = __rdtscp(&coreid);
      __asm__("lfence");
#endif
     GMS_CUBLAS_STAT_CHECK(cublasSasum(handle,n,d_ptr,incx,&result));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      tsc_end     = __rdtscp(&coreid);
      *tsc_delta  = tsc_end-tsc_start-rdtscp_cost;
      __asm__("lfence");
      printf("cublasSasum executed in: %llu reference cycles\n",*tsc_delta);
#endif
      GMS_CUDA_DEBUG_CHECK(cudaStreamSynchronize(stream));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptr));
      GMS_CUBLAS_STAT_CHECK(cublasDestroy(handle));
      GMS_CUDA_DEBUG_CHECK(cudaStreamDestroy(stream));
      *ierr = 0;
      *err  = status;
      return;
      Error: 
            {
                if(d_ptr) cudaFree(d_ptr);
                cublasDestroy(handle);
                cudaStreamDestroy(stream);  
                ierr = -2;
                *err = status;  
                return;  
      }
}


void
cuBLAS_Scasum_iface(const cuComplex * __restrict h_ptr,
                   const int32_t n,
                   const int32_t incx,
                   cuComplex * __restrict result,
                   cudaError_t * __restrict err,
                   int32_t * __restrict ierr,
                   uint64_t * __restrict tsc_delta) {
    
     if(__builtin_expect(n<=1,0)) {
       *ierr = -1;
       return;
      }
    *err    = 0;
    *ierr   = 0;
    *result = 0;
#if (PROFILE_HOST_TO_DEVICE) == 1
     volatile uint64_t dummy1;
     volatile uint32_t dummy2;
     volatile uint64_t tsc_start,tsc_end;
     volatile uint32_t coreid;
#endif
    cuComplex * __restrict d_ptr = NULL;
    cublasHandle_t handle    = NULL;
    cudaStream_t   stream    = NULL;
    cudaError_t_t status;
    int32_t merr = 0;
    GMS_CUBLAS_STAT_CHECK(cublasCreate(&handle));
    GMS_CUDA_DEBUG_CHECK(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
    GMS_CUBLAS_STAT_CHECK(cublasSetStream(handle,stream));
    alloc_complex4_gpu(&d_ptr[0],(size_t)n,&merr);
    if(merr != 0) goto Error;
    GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptr,h_ptr,sizeof(cuComplex)*((size_t)n),
                         cudaMemcpyHostToDevice,stream));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      dummy1    = __rdtscp(&dummy2);
      tsc_start = __rdtscp(&coreid);
      __asm__("lfence");
#endif
     GMS_CUBLAS_STAT_CHECK(cublasScasum(handle,n,d_ptr,incx,&result));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      tsc_end     = __rdtscp(&coreid);
      *tsc_delta  = tsc_end-tsc_start-rdtscp_cost;
      __asm__("lfence");
      printf("cublasScasum executed in: %llu reference cycles\n",*tsc_delta);
#endif
      GMS_CUDA_DEBUG_CHECK(cudaStreamSynchronize(stream));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptr));
      GMS_CUBLAS_STAT_CHECK(cublasDestroy(handle));
      GMS_CUDA_DEBUG_CHECK(cudaStreamDestroy(stream));
      *ierr = 0;
      *err  = status;
      return;
      Error: 
            {
                if(d_ptr) cudaFree(d_ptr));
                cublasDestroy(handle);
                cudaStreamDestroy(stream);  
                ierr = -2;
                *err = status;   
                return; 
      }
} 




  


