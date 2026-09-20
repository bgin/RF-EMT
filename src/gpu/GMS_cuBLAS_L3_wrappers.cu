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
#include "GMS_cuBLAS_L3_wrappers.cuh"
#include "GMS_gpu_config.cuh"
#if (PROFILE_HOST_TO_DEVICE) == 1
#include <immintrin.h> //rdtscp
#endif 
#include "GMS_cuda_memops.cuh"


    #define GMS_MIN(a,b) (((a)<(b))?(a):(b))

    #define GMS_MAX(a,b) (((a)>(b))?(a):(b))


void
cuBLAS_Sgemm_iface(cublasOperation_t transa,
                   cublasOperation_t trnasb,
                   const int32_t m,
                   const int32_t n,
                   const int32_t k,
                   float alpha,
                   const float * __restrict A,
                   const int32_t LDA,
                   const float * __restrict B,
                   const int32_t LDB,
                   float beta,
                   float * __restrict C,
                   const int32_t LDC,
                   cudaError_t * __restrict err,
                   int32_t * __restrict ierr,
                   uint64_t * __restrict tsc_delta)  {

    if(__builtin_expect(n<=1,0)) {
       *ierr = -1;
       return;
      }
    *err    = 0;
    *ierr   = 0;
     
#if (PROFILE_HOST_TO_DEVICE) == 1
     volatile uint64_t dummy1;
     volatile uint32_t dummy2;
     volatile uint64_t tsc_start,tsc_end;
     volatile uint32_t coreid;
#endif
     cublasHandle_t handle       = NULL;
     cudaStream_t   stream       = NULL;
     float * __restrict   d_ptrA = NULL;
     float * __restrict   d_ptrB = NULL;
     float * __restrict   d_ptrC = NULL;
     float alph = alpha;
     float bet  = beta;
     cudaError_t_t status;
     int32_t merr = 0;
     GMS_CUBLAS_STAT_CHECK(cublasCreate(&handle));
     GMS_CUDA_DEBUG_CHECK(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
     GMS_CUBLAS_STAT_CHECK(cublasSetStream(handle,stream));
     alloc_float_gpu(d_ptrA,(size_t)(m*k),&merr);
     if(merr != 0) goto Error;
     alloc_float_gpu(d_ptrB,(size_t)(n*k),&merr);
     if(merr != 0) goto Error;
     alloc_float_gpu(d_ptrC,(size_t)(m*n),&merr);
     if(merr != 0) goto Error;
     GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrA,A,sizeof(float)*(m*k),
                          cudaMemcpyHostToDevice,stream));
     GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrB,B,sizeof(float)*(n*k),
                          cudaMemcpyHostToDevice,stream));
#if (PROFILE_HOST_TO_DEVICE) == 1
      dummy1    = __rdtscp(&dummy2);
      __asm__("lfence");
      tsc_start = __rdtscp(&coreid);
      __asm__("lfence");
#endif
      GMS_CUBLAS_STAT_CHECK(cublasSgemm(handle,transa,transb,m,n,k,&alph,
                                        d_ptrA,LDA,d_ptrB,LDB,&bet,d_ptrC,LDC));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      tsc_end     = __rdtscp(&coreid);
      *tsc_delta  = tsc_end-tsc_start-rdtscp_cost;
      __asm__("lfence");
      printf("cublasSgemm executed in: %llu reference cycles\n",*tsc_delta);
#endif  
      GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrC,d_C,sizeof(float)*(m*n),
                                           cudaMemcpyDeviceToHost,stream));
      GMS_CUDA_DEBUG_CHECK(cudaStreamSynchronize(stream));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrA));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrB));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrC));
      GMS_CUBLAS_STAT_CHECK(cublasDestroy(handle));
      GMS_CUDA_DEBUG_CHECK(cudaStreamDestroy(stream));
      *ierr = 0;
      *err  = status;
      return; 
Error:   {  
           if(d_ptrA) cudaFree(d_ptrA);
           if(d_ptrB) cudaFree(d_ptrB);
           if(d_ptrC) cudaFree(d_ptrC);
           cublasDestroy(handle);
           cudaStreamDestroy(stream); 
           ierr = -2;
           *err = status;   
           return; 
    }
}


void
cuBLAS_Cgemm_iface(cublasOperation_t transa,
                   cublasOperation_t transb,
                   const int32_t m,
                   const int32_t n,
                   const int32_t k,
                   const cuComplex alpha,
                   const cuComplex * __restrict A,
                   const int32_t LDA,
                   const cuComplex * __restrict B,
                   const int32_t LDB,
                   cuComplex beta,
                   cuComplex * __restrict C,
                   const int32_t LDC,
                   cudaError_t * __restrict err,
                   int32_t * __restrict ierr,
                   uint64_t * __restrict tsc_delta) {

   if(__builtin_expect(n<=1,0)) {
       *ierr = -1;
       return;
      }
    *err    = 0;
    *ierr   = 0;
     
#if (PROFILE_HOST_TO_DEVICE) == 1
     volatile uint64_t dummy1;
     volatile uint32_t dummy2;
     volatile uint64_t tsc_start,tsc_end;
     volatile uint32_t coreid;
#endif
     cublasHandle_t handle       = NULL;
     cudaStream_t   stream       = NULL;
     cuComplex * __restrict   d_ptrA = NULL;
     cuComplex * __restrict   d_ptrB = NULL;
     cuComplex * __restrict   d_ptrC = NULL;
     cuComplex alph = alpha;
     cuComplex bet  = beta;
     cudaError_t_t status;
     int32_t merr = 0;
     GMS_CUBLAS_STAT_CHECK(cublasCreate(&handle));
     GMS_CUDA_DEBUG_CHECK(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
     GMS_CUBLAS_STAT_CHECK(cublasSetStream(handle,stream));
     alloc_complex4_gpu(d_ptrA,(size_t)(m*k),&merr);
     if(merr != 0) goto Error;
     alloc_complex4_gpu(d_ptrB,(size_t)(n*k),&merr);
     if(merr != 0) goto Error;
     alloc_complex4_gpu(d_ptrC,(size_t)(m*n),&merr);
     if(merr != 0) goto Error;
     GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrA,A,sizeof(cuComplex)*(m*k),
                          cudaMemcpyHostToDevice,stream));
     GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrB,B,sizeof(cuComplex)*(n*k),
                          cudaMemcpyHostToDevice,stream));
#if (PROFILE_HOST_TO_DEVICE) == 1
      dummy1    = __rdtscp(&dummy2);
      __asm__("lfence");
      tsc_start = __rdtscp(&coreid);
      __asm__("lfence");
#endif
      GMS_CUBLAS_STAT_CHECK(cublasCgemm(handle,transa,transb,m,n,k,&alph,
                                        d_ptrA,LDA,d_ptrB,LDB,&bet,d_ptrC,LDC));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      tsc_end     = __rdtscp(&coreid);
      *tsc_delta  = tsc_end-tsc_start-rdtscp_cost;
      __asm__("lfence");
      printf("cublasCgemm executed in: %llu reference cycles\n",*tsc_delta);
#endif  
      GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrC,d_C,sizeof(cuComplex)*(m*n),
                                           cudaMemcpyDeviceToHost,stream));
      GMS_CUDA_DEBUG_CHECK(cudaStreamSynchronize(stream));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrA));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrB));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrC));
      GMS_CUBLAS_STAT_CHECK(cublasDestroy(handle));
      GMS_CUDA_DEBUG_CHECK(cudaStreamDestroy(stream));
      *ierr = 0;
      *err  = status;
      return; 
Error:   {  
           if(d_ptrA) cudaFree(d_ptrA);
           if(d_ptrB) cudaFree(d_ptrB);
           if(d_ptrC) cudaFree(d_ptrC);
           cublasDestroy(handle);
           cudaStreamDestroy(stream); 
           ierr = -2;
           *err = status;   
           return; 
    }  
} 



void
cuBLAS_Cgemm3m_iface(cublasOperation_t transa,
                   cublasOperation_t transb,
                   const int32_t m,
                   const int32_t n,
                   const int32_t k,
                   const cuComplex alpha,
                   const cuComplex * __restrict A,
                   const int32_t LDA,
                   const cuComplex * __restrict B,
                   const int32_t LDB,
                   cuComplex beta,
                   cuComplex * __restrict C,
                   const int32_t LDC,
                   cudaError_t * __restrict err,
                   int32_t * __restrict ierr,
                   uint64_t * __restrict tsc_delta) {

   if(__builtin_expect(n<=1,0)) {
       *ierr = -1;
       return;
      }
    *err    = 0;
    *ierr   = 0;
     
#if (PROFILE_HOST_TO_DEVICE) == 1
     volatile uint64_t dummy1;
     volatile uint32_t dummy2;
     volatile uint64_t tsc_start,tsc_end;
     volatile uint32_t coreid;
#endif
     cublasHandle_t handle       = NULL;
     cudaStream_t   stream       = NULL;
     cuComplex * __restrict   d_ptrA = NULL;
     cuComplex * __restrict   d_ptrB = NULL;
     cuComplex * __restrict   d_ptrC = NULL;
     cuComplex alph = alpha;
     cuComplex bet  = beta;
     cudaError_t_t status;
     int32_t merr = 0;
     GMS_CUBLAS_STAT_CHECK(cublasCreate(&handle));
     GMS_CUDA_DEBUG_CHECK(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
     GMS_CUBLAS_STAT_CHECK(cublasSetStream(handle,stream));
     alloc_complex4_gpu(d_ptrA,(size_t)(m*k),&merr);
     if(merr != 0) goto Error;
     alloc_complex4_gpu(d_ptrB,(size_t)(n*k),&merr);
     if(merr != 0) goto Error;
     alloc_complex4_gpu(d_ptrC,(size_t)(m*n),&merr);
     if(merr != 0) goto Error;
     GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrA,A,sizeof(cuComplex)*(m*k),
                          cudaMemcpyHostToDevice,stream));
     GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrB,B,sizeof(cuComplex)*(n*k),
                          cudaMemcpyHostToDevice,stream));
#if (PROFILE_HOST_TO_DEVICE) == 1
      dummy1    = __rdtscp(&dummy2);
      __asm__("lfence");
      tsc_start = __rdtscp(&coreid);
      __asm__("lfence");
#endif
      GMS_CUBLAS_STAT_CHECK(cublasCgemm3m(handle,transa,transb,m,n,k,&alph,
                                        d_ptrA,LDA,d_ptrB,LDB,&bet,d_ptrC,LDC));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      tsc_end     = __rdtscp(&coreid);
      *tsc_delta  = tsc_end-tsc_start-rdtscp_cost;
      __asm__("lfence");
      printf("cublasCgemm3m executed in: %llu reference cycles\n",*tsc_delta);
#endif  
      GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrC,d_C,sizeof(cuComplex)*(m*n),
                                           cudaMemcpyDeviceToHost,stream));
      GMS_CUDA_DEBUG_CHECK(cudaStreamSynchronize(stream));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrA));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrB));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrC));
      GMS_CUBLAS_STAT_CHECK(cublasDestroy(handle));
      GMS_CUDA_DEBUG_CHECK(cudaStreamDestroy(stream));
      *ierr = 0;
      *err  = status;
      return; 
Error:   {  
           if(d_ptrA) cudaFree(d_ptrA);
           if(d_ptrB) cudaFree(d_ptrB);
           if(d_ptrC) cudaFree(d_ptrC);
           cublasDestroy(handle);
           cudaStreamDestroy(stream); 
           ierr = -2;
           *err = status;   
           return; 
    }  
} 



void
cuBLAS_Chemm_iface(cublasSideMode_t side,
                   cublasFillMode_t uplo,
                   const int32_t m,
                   const int32_t n,
                   const cuComplex alpha,
                   const cuComplex * __restrict A,
                   const int32_t LDA,
                   const cuComplex * __restrict B,
                   const int32_t LDB,
                   const cuComplex beta,
                   cuComplex * __restrict C,
                   const int32_t LDC,
                   cudaError_t * __restrict err,
                   int32_t * __restrict ierr,
                   uint64_t * __restrict tsc_delta) {\

    if(__builtin_expect(n<=1,0)) {
       *ierr = -1;
       return;
      }
    *err    = 0;
    *ierr   = 0;
     
#if (PROFILE_HOST_TO_DEVICE) == 1
     volatile uint64_t dummy1;
     volatile uint32_t dummy2;
     volatile uint64_t tsc_start,tsc_end;
     volatile uint32_t coreid;
#endif
     cublasHandle_t handle       = NULL;
     cudaStream_t   stream       = NULL;
     cuComplex * __restrict   d_ptrA = NULL;
     cuComplex * __restrict   d_ptrB = NULL;
     cuComplex * __restrict   d_ptrC = NULL;
     cuComplex alph = alpha;
     cuComplex bet  = beta;
     cudaError_t_t status;
     int32_t merr = 0;
     GMS_CUBLAS_STAT_CHECK(cublasCreate(&handle));
     GMS_CUDA_DEBUG_CHECK(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
     GMS_CUBLAS_STAT_CHECK(cublasSetStream(handle,stream));
     if(side==CUBLAS_SIDE_LEFT) {
        alloc_complex4_gpu(d_ptrA,(size_t)(LDA*m),&merr);
     }else{
        alloc_complex4_gpu(d_ptrA,(size_t)(LDA*n),&merr);
     {
     if(merr != 0) goto Error;
     alloc_complex4_gpu(d_ptrB,(size_t)(LDB*n),&merr);
     if(merr != 0) goto Error;
     alloc_complex4_gpu(d_ptrC,(size_t)(LDC*n),&merr);
     if(merr != 0) goto Error;
     if(side==CUBLAS_SIDE_LEFT) {
        GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrA,A,sizeof(cuComplex)*(LDA*m),
                          cudaMemcpyHostToDevice,stream));
     }else {
        GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrA,A,sizeof(cuComplex)*(LDA*n),
                          cudaMemcpyHostToDevice,stream));
     } 
     GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrB,B,sizeof(cuComplex)*(LDB*n),
                          cudaMemcpyHostToDevice,stream));
#if (PROFILE_HOST_TO_DEVICE) == 1
      dummy1    = __rdtscp(&dummy2);
      __asm__("lfence");
      tsc_start = __rdtscp(&coreid);
      __asm__("lfence");
#endif
      GMS_CUBLAS_STAT_CHECK(cublasChemm(handle,side,uplo,m,n,&alph,d_ptrA,LDA,
                                        d_ptrB,LDB,&bet,d_ptrC,LDC));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      tsc_end     = __rdtscp(&coreid);
      *tsc_delta  = tsc_end-tsc_start-rdtscp_cost;
      __asm__("lfence");
      printf("cublasChemm executed in: %llu reference cycles\n",*tsc_delta);
#endif  
      GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrC,d_C,sizeof(cuComplex)*(LDC*n),
                                           cudaMemcpyDeviceToHost,stream));
      GMS_CUDA_DEBUG_CHECK(cudaStreamSynchronize(stream));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrA));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrB));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrC));
      GMS_CUBLAS_STAT_CHECK(cublasDestroy(handle));
      GMS_CUDA_DEBUG_CHECK(cudaStreamDestroy(stream));
      *ierr = 0;
      *err  = status;
      return; 
Error:   {  
           if(d_ptrA) cudaFree(d_ptrA);
           if(d_ptrB) cudaFree(d_ptrB);
           if(d_ptrC) cudaFree(d_ptrC);
           cublasDestroy(handle);
           cudaStreamDestroy(stream); 
           ierr = -2;
           *err = status;   
           return; 
    }  


}


void
cuBLAS_Cher2k_iface(cublasOperation_t transa,
                    cublasOperation_t transb,
                    cublasFillMode_t  uplo,
                    const int32_t n,
                    const int32_t k,
                    const cuComplex alpha,
                    const cuComplex * __restrict A,
                    const int32_t LDA,
                    const cuComplex * __restrict B,
                    const int32_t LDB,
                    const cuComplex beta,
                    cuComplex * __restrict C,
                    const int32_t LDC,
                    cudaError_t * __restrict err,
                    int32_t * __restrict ierr,
                    uint64_t * __restrict tsc_delta) {

    if(__builtin_expect(n<=1,0)) {
       *ierr = -1;
       return;
      }
    *err    = 0;
    *ierr   = 0;
     
#if (PROFILE_HOST_TO_DEVICE) == 1
     volatile uint64_t dummy1;
     volatile uint32_t dummy2;
     volatile uint64_t tsc_start,tsc_end;
     volatile uint32_t coreid;
#endif
     cublasHandle_t handle       = NULL;
     cudaStream_t   stream       = NULL;
     cuComplex * __restrict   d_ptrA = NULL;
     cuComplex * __restrict   d_ptrB = NULL;
     cuComplex * __restrict   d_ptrC = NULL;
     cuComplex alph = alpha;
     cuComplex bet  = beta;
     cudaError_t_t status;
     int32_t merr = 0;
     GMS_CUBLAS_STAT_CHECK(cublasCreate(&handle));
     GMS_CUDA_DEBUG_CHECK(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
     GMS_CUBLAS_STAT_CHECK(cublasSetStream(handle,stream));
     if(transa==CUBLAS_OP_N) {
         alloc_complex4_gpu(d_ptrA,(size_t)(LDA*k),&merr);
     }else{
         alloc_complex4_gpu(d_ptrA,(size_t)(LDA*n),&merr);
     {
     if(merr != 0) goto Error;
     if(transb==CUBLAS_OP_N) {
        alloc_complex4_gpu(d_ptrB,(size_t)(LDB*k),&merr);
     }
      else {
         alloc_complex4_gpu(d_ptrB,(size_t)(LDB*n),&merr); 
     }
     if(merr != 0) goto Error;
     alloc_complex4_gpu(d_ptrC,(size_t)(LDC*n),&merr);
     if(merr != 0) goto Error;
     if(transa==CUBLAS_OP_N) {
        GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrA,A,sizeof(cuComplex)*(LDA*k),
                          cudaMemcpyHostToDevice,stream));
     }else {
         GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrA,A,sizeof(cuComplex)*(LDA*n),
                          cudaMemcpyHostToDevice,stream));
     }
     if(transb==CUBLAS_OP_N) {
         GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrB,B,sizeof(cuComplex)*(LDB*k),
                          cudaMemcpyHostToDevice,stream));
     }else {
         GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrB,B,sizeof(cuComplex)*(LDB*n),
                          cudaMemcpyHostToDevice,stream));
     }
#if (PROFILE_HOST_TO_DEVICE) == 1
      dummy1    = __rdtscp(&dummy2);
      __asm__("lfence");
      tsc_start = __rdtscp(&coreid);
      __asm__("lfence");
#endif
      GMS_CUBLAS_STAT_CHECK(cublasCher2k(handle,uplo,transa,n,k,&alph,d_ptrA,LDA,
                                         d_ptrB,LDB,&bet,d_ptrC,LDC));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      tsc_end     = __rdtscp(&coreid);
      *tsc_delta  = tsc_end-tsc_start-rdtscp_cost;
      __asm__("lfence");
      printf("cublasCher2k executed in: %llu reference cycles\n",*tsc_delta);
#endif  
      GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrC,d_C,sizeof(cuComplex)*(LDC*n),
                                           cudaMemcpyDeviceToHost,stream));
      GMS_CUDA_DEBUG_CHECK(cudaStreamSynchronize(stream));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrA));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrB));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrC));
      GMS_CUBLAS_STAT_CHECK(cublasDestroy(handle));
      GMS_CUDA_DEBUG_CHECK(cudaStreamDestroy(stream));
      *ierr = 0;
      *err  = status;
      return; 
Error:   {  
           if(d_ptrA) cudaFree(d_ptrA);
           if(d_ptrB) cudaFree(d_ptrB);
           if(d_ptrC) cudaFree(d_ptrC);
           cublasDestroy(handle);
           cudaStreamDestroy(stream); 
           ierr = -2;
           *err = status;   
           return; 
    }  

} 


void
cuBLAS_Ssymm_iface(cublasSideMode_t side,
                   cublasFillMode_t uplo,
                   int32_t m,
                   int32_t n,
                   const float alpha,
                   const float * __restrict A,
                   const int32_t LDA,
                   const float * __restrict B,
                   int32_t LDB,
                   const float beta,
                   float * __restrict C,
                   const int32_t LDC,
                   cudaError_t * __restrict,
                   int32_t * __restrict,
                   uint64_t * __restrict)  {
    
     if(__builtin_expect(n<=1,0)) {
       *ierr = -1;
       return;
      }
    *err    = 0;
    *ierr   = 0;
     
#if (PROFILE_HOST_TO_DEVICE) == 1
     volatile uint64_t dummy1;
     volatile uint32_t dummy2;
     volatile uint64_t tsc_start,tsc_end;
     volatile uint32_t coreid;
#endif
     cublasHandle_t handle       = NULL;
     cudaStream_t   stream       = NULL;
     float * __restrict   d_ptrA = NULL;
     float * __restrict   d_ptrB = NULL;
     float * __restrict   d_ptrC = NULL;
     float alph = alpha;
     float bet  = beta;
     cudaError_t_t status;
     int32_t merr = 0;
     GMS_CUBLAS_STAT_CHECK(cublasCreate(&handle));
     GMS_CUDA_DEBUG_CHECK(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
     GMS_CUBLAS_STAT_CHECK(cublasSetStream(handle,stream));
     if(side==CUBLAS_SIDE_LEFT) {
         alloc_float_gpu(d_ptrA,(size_t)(LDA*m),&merr);
     }else{
         alloc_float_gpu(d_ptrA,(size_t)(LDA*n),&merr);
     {
     if(merr != 0) goto Error;
     alloc_float_gpu(d_ptrB,(size_t)(LDB*n),&merr);
    
     if(merr != 0) goto Error;
     alloc_float_gpu(d_ptrC,(size_t)(LDC*n),&merr);
     if(merr != 0) goto Error;
     if(side==CUBLAS_SIDE_LEFT) {
        GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrA,A,sizeof(float)*(LDA*m),
                          cudaMemcpyHostToDevice,stream));
      }else {
        GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrA,A,sizeof(float)*(LDA*n),
                          cudaMemcpyHostToDevice,stream)); 
      }
     GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrB,B,sizeof(float)*(LDB*n),
                          cudaMemcpyHostToDevice,stream));
#if (PROFILE_HOST_TO_DEVICE) == 1
      dummy1    = __rdtscp(&dummy2);
      __asm__("lfence");
      tsc_start = __rdtscp(&coreid);
      __asm__("lfence");
#endif
      GMS_CUBLAS_STAT_CHECK(cublasSsymm(handle,side,uplo,m,n,&alph,d_ptrA,LDA,
                                         d_ptrB,LDB,&bet,d_ptrC,LDC));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      tsc_end     = __rdtscp(&coreid);
      *tsc_delta  = tsc_end-tsc_start-rdtscp_cost;
      __asm__("lfence");
      printf("cublasSsymm executed in: %llu reference cycles\n",*tsc_delta);
#endif  
      GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrC,d_C,sizeof(float)*(LDC*n),
                                           cudaMemcpyDeviceToHost,stream));
      GMS_CUDA_DEBUG_CHECK(cudaStreamSynchronize(stream));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrA));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrB));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrC));
      GMS_CUBLAS_STAT_CHECK(cublasDestroy(handle));
      GMS_CUDA_DEBUG_CHECK(cudaStreamDestroy(stream));
      *ierr = 0;
      *err  = status;
      return; 
Error:   {  
           if(d_ptrA) cudaFree(d_ptrA);
           if(d_ptrB) cudaFree(d_ptrB);
           if(d_ptrC) cudaFree(d_ptrC);
           cublasDestroy(handle);
           cudaStreamDestroy(stream); 
           ierr = -2;
           *err = status;   
           return; 
    }  


}


void
cuBLAS_Ssyrk_iface(cublasFillMode_t uplo,
                   cublasOperation_t trans,
                   const int32_t n,
                   const int32_t k,
                   const float alpha,
                   const float * __restrict A,
                   const int32_t LDA,
                   const float beta,
                   float * __restrict C,
                   const int32_t LDC,
                   cudaError_t * __restrict err,
                   int32_t * __restrict ierr,
                   uint64_t * __restrict tsc_delta) {

     if(__builtin_expect(n<=1,0)) {
       *ierr = -1;
       return;
      }
    *err    = 0;
    *ierr   = 0;
     
#if (PROFILE_HOST_TO_DEVICE) == 1
     volatile uint64_t dummy1;
     volatile uint32_t dummy2;
     volatile uint64_t tsc_start,tsc_end;
     volatile uint32_t coreid;
#endif
     cublasHandle_t handle       = NULL;
     cudaStream_t   stream       = NULL;
     float * __restrict   d_ptrA = NULL;
     float * __restrict   d_ptrC = NULL;
     float alph = alpha;
     float bet  = beta;
     cudaError_t_t status;
     int32_t merr = 0;
     GMS_CUBLAS_STAT_CHECK(cublasCreate(&handle));
     GMS_CUDA_DEBUG_CHECK(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
     GMS_CUBLAS_STAT_CHECK(cublasSetStream(handle,stream));
     if(trans==CUBLAS_OP_N) {
         alloc_float_gpu(d_ptrA,(size_t)(LDA*k),&merr);
     }else{
         alloc_float_gpu(d_ptrA,(size_t)(LDA*n),&merr);
     {
     if(merr != 0) goto Error;
     alloc_float_gpu(d_ptrC,(size_t)(LDC*n),&merr);
     if(merr != 0) goto Error;
     if(trans==CUBLAS_OP_N) {
        GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrA,A,sizeof(float)*(LDA*k),
                          cudaMemcpyHostToDevice,stream));
      }else {
        GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrA,A,sizeof(float)*(LDA*n),
                          cudaMemcpyHostToDevice,stream)); 
      }
     
#if (PROFILE_HOST_TO_DEVICE) == 1
      dummy1    = __rdtscp(&dummy2);
      __asm__("lfence");
      tsc_start = __rdtscp(&coreid);
      __asm__("lfence");
#endif
      GMS_CUBLAS_STAT_CHECK(cublasSsymm(handle,uplo,trans,n,k,&alph,d_ptrA,LDA,
                                        &bet,d_ptrC,LDC));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      tsc_end     = __rdtscp(&coreid);
      *tsc_delta  = tsc_end-tsc_start-rdtscp_cost;
      __asm__("lfence");
      printf("cublasSsyrk executed in: %llu reference cycles\n",*tsc_delta);
#endif  
      GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrC,d_C,sizeof(float)*(LDC*n),
                                           cudaMemcpyDeviceToHost,stream));
      GMS_CUDA_DEBUG_CHECK(cudaStreamSynchronize(stream));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrA));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrC));
      GMS_CUBLAS_STAT_CHECK(cublasDestroy(handle));
      GMS_CUDA_DEBUG_CHECK(cudaStreamDestroy(stream));
      *ierr = 0;
      *err  = status;
      return; 
Error:   {  
           if(d_ptrA) cudaFree(d_ptrA);
           if(d_ptrC) cudaFree(d_ptrC);
           cublasDestroy(handle);
           cudaStreamDestroy(stream); 
           ierr = -2;
           *err = status;   
           return; 
    }  


} 


void
cuBLAS_Csyrk_iface(cublasFillMode_t uplo,
                   cublasOperation_t trans,
                   const int32_t n,
                   const int32_t k,
                   const cuComplex alpha,
                   const cuComplex * __restrict A,
                   const int32_t LDA,
                   const cuComplex beta,
                   cuComplex * __restrict C,
                   const int32_t LDC,
                   cudaError_t * __restrict err,
                   int32_t * __restrict ierr,
                   uint64_t * __restrict tsc_delta) {

     if(__builtin_expect(n<=1,0)) {
       *ierr = -1;
       return;
      }
    *err    = 0;
    *ierr   = 0;
     
#if (PROFILE_HOST_TO_DEVICE) == 1
     volatile uint64_t dummy1;
     volatile uint32_t dummy2;
     volatile uint64_t tsc_start,tsc_end;
     volatile uint32_t coreid;
#endif
     cublasHandle_t handle       = NULL;
     cudaStream_t   stream       = NULL;
     cuComplex * __restrict   d_ptrA = NULL;
     cuComplex * __restrict   d_ptrC = NULL;
     cuComplex alph = alpha;
     cuComplex bet  = beta;
     cudaError_t_t status;
     int32_t merr = 0;
     GMS_CUBLAS_STAT_CHECK(cublasCreate(&handle));
     GMS_CUDA_DEBUG_CHECK(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
     GMS_CUBLAS_STAT_CHECK(cublasSetStream(handle,stream));
     if(trans==CUBLAS_OP_N) {
         alloc_complex4_gpu(d_ptrA,(size_t)(LDA*k),&merr);
     }else{
         alloc_complex4_gpu(d_ptrA,(size_t)(LDA*n),&merr);
     {
     if(merr != 0) goto Error;
     alloc_complex4_gpu(d_ptrC,(size_t)(LDC*n),&merr);
     if(merr != 0) goto Error;
     if(trans==CUBLAS_OP_N) {
        GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrA,A,sizeof(cuComplex)*(LDA*k),
                          cudaMemcpyHostToDevice,stream));
      }else {
        GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrA,A,sizeof(cuComplex)*(LDA*n),
                          cudaMemcpyHostToDevice,stream)); 
      }
     
#if (PROFILE_HOST_TO_DEVICE) == 1
      dummy1    = __rdtscp(&dummy2);
      __asm__("lfence");
      tsc_start = __rdtscp(&coreid);
      __asm__("lfence");
#endif
      GMS_CUBLAS_STAT_CHECK(cublasCsymm(handle,uplo,trans,n,k,&alph,d_ptrA,LDA,
                                        &bet,d_ptrC,LDC));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      tsc_end     = __rdtscp(&coreid);
      *tsc_delta  = tsc_end-tsc_start-rdtscp_cost;
      __asm__("lfence");
      printf("cublasCsyrk executed in: %llu reference cycles\n",*tsc_delta);
#endif  
      GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrC,d_C,sizeof(cuComplex)*(LDC*n),
                                           cudaMemcpyDeviceToHost,stream));
      GMS_CUDA_DEBUG_CHECK(cudaStreamSynchronize(stream));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrA));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrC));
      GMS_CUBLAS_STAT_CHECK(cublasDestroy(handle));
      GMS_CUDA_DEBUG_CHECK(cudaStreamDestroy(stream));
      *ierr = 0;
      *err  = status;
      return; 
Error:   {  
           if(d_ptrA) cudaFree(d_ptrA);
           if(d_ptrC) cudaFree(d_ptrC);
           cublasDestroy(handle);
           cudaStreamDestroy(stream); 
           ierr = -2;
           *err = status;   
           return; 
    }  


} 


void
cuBLAS_Strsm_iface(cublasSideMode_t side,
                   cublasFillMode_t uplo,
                   cublasOperation_t trans,
                   cublasDiagType_t diag,
                   const int32_t m,
                   const int32_t n,
                   const float alpha,
                   const float * __restrict A,
                   const int32_t LDA,
                   float * __restrict B,
                   const int32_t LDB,
                   cudaError_t * __restrict err,
                   int32_t * __restrict ierr,
                   uint64_t * __restrict tsc_delta) {

     if(__builtin_expect(n<=1,0)) {
       *ierr = -1;
       return;
      }
    *err    = 0;
    *ierr   = 0;
     
#if (PROFILE_HOST_TO_DEVICE) == 1
     volatile uint64_t dummy1;
     volatile uint32_t dummy2;
     volatile uint64_t tsc_start,tsc_end;
     volatile uint32_t coreid;
#endif
     cublasHandle_t handle       = NULL;
     cudaStream_t   stream       = NULL;
     float * __restrict   d_ptrA = NULL;
     float * __restrict   d_ptrB = NULL;
     float alph = alpha;
     float bet  = beta;
     cudaError_t_t status;
     int32_t merr = 0;
     GMS_CUBLAS_STAT_CHECK(cublasCreate(&handle));
     GMS_CUDA_DEBUG_CHECK(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
     GMS_CUBLAS_STAT_CHECK(cublasSetStream(handle,stream));
     if(side==CUBLAS_SIDE_LEFT) {
         alloc_float_gpu(d_ptrA,(size_t)(LDA*m),&merr);
     }else{
         alloc_float_gpu(d_ptrA,(size_t)(LDA*n),&merr);
     {
     if(merr != 0) goto Error;
     alloc_float_gpu(d_ptrB,(size_t)(LDB*n),&merr);
     if(merr != 0) goto Error;
     if(side==CUBLAS_SIDE_LEFT) {
        GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrA,A,sizeof(float)*(LDA*m),
                          cudaMemcpyHostToDevice,stream));
      }else {
        GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrA,A,sizeof(float)*(LDA*n),
                          cudaMemcpyHostToDevice,stream)); 
      }
     
#if (PROFILE_HOST_TO_DEVICE) == 1
      dummy1    = __rdtscp(&dummy2);
      __asm__("lfence");
      tsc_start = __rdtscp(&coreid);
      __asm__("lfence");
#endif
      GMS_CUBLAS_STAT_CHECK(cublasStrsm(handle,side,uplo,trans,diag,m,n,&alph,d_ptrA,LDA,
                                        &bet,d_ptrB,LDB));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      tsc_end     = __rdtscp(&coreid);
      *tsc_delta  = tsc_end-tsc_start-rdtscp_cost;
      __asm__("lfence");
      printf("cublasStrsm executed in: %llu reference cycles\n",*tsc_delta);
#endif  
      GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrB,d_B,sizeof(float)*(LDB*n),
                                           cudaMemcpyDeviceToHost,stream));
      GMS_CUDA_DEBUG_CHECK(cudaStreamSynchronize(stream));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrA));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrB));
      GMS_CUBLAS_STAT_CHECK(cublasDestroy(handle));
      GMS_CUDA_DEBUG_CHECK(cudaStreamDestroy(stream));
      *ierr = 0;
      *err  = status;
      return; 
Error:   {  
           if(d_ptrA) cudaFree(d_ptrA);
           if(d_ptrC) cudaFree(d_ptrB);
           cublasDestroy(handle);
           cudaStreamDestroy(stream); 
           ierr = -2;
           *err = status;   
           return; 
    }  


} 


void
cuBLAS_Ctrsm_iface(cublasSideMode_t side,
                   cublasFillMode_t uplo,
                   cublasOperation_t trans,
                   cublasDiagType_t diag,
                   const int32_t m,
                   const int32_t n,
                   const cuComplex alpha,
                   const cuComplex * __restrict A,
                   const int32_t LDA,
                   cuComplex * __restrict B,
                   const int32_t LDB,
                   cudaError_t * __restrict err,
                   int32_t * __restrict ierr,
                   uint64_t * __restrict tsc_delta) {

     if(__builtin_expect(n<=1,0)) {
       *ierr = -1;
       return;
      }
    *err    = 0;
    *ierr   = 0;
     
#if (PROFILE_HOST_TO_DEVICE) == 1
     volatile uint64_t dummy1;
     volatile uint32_t dummy2;
     volatile uint64_t tsc_start,tsc_end;
     volatile uint32_t coreid;
#endif
     cublasHandle_t handle       = NULL;
     cudaStream_t   stream       = NULL;
     cuComplex * __restrict   d_ptrA = NULL;
     cuComplex * __restrict   d_ptrB = NULL;
     cuComplex alph = alpha;
     cuComplex bet  = beta;
     cudaError_t_t status;
     int32_t merr = 0;
     GMS_CUBLAS_STAT_CHECK(cublasCreate(&handle));
     GMS_CUDA_DEBUG_CHECK(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
     GMS_CUBLAS_STAT_CHECK(cublasSetStream(handle,stream));
     if(side==CUBLAS_SIDE_LEFT) {
         alloc_complex4_gpu(d_ptrA,(size_t)(LDA*m),&merr);
     }else{
         alloc_complex4_gpu(d_ptrA,(size_t)(LDA*n),&merr);
     {
     if(merr != 0) goto Error;
     alloc_complex4_gpu(d_ptrB,(size_t)(LDB*n),&merr);
     if(merr != 0) goto Error;
     if(side==CUBLAS_SIDE_LEFT) {
        GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrA,A,sizeof(cuComplex)*(LDA*m),
                          cudaMemcpyHostToDevice,stream));
      }else {
        GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrA,A,sizeof(cuComplex)*(LDA*n),
                          cudaMemcpyHostToDevice,stream)); 
      }
     
#if (PROFILE_HOST_TO_DEVICE) == 1
      dummy1    = __rdtscp(&dummy2);
      __asm__("lfence");
      tsc_start = __rdtscp(&coreid);
      __asm__("lfence");
#endif
      GMS_CUBLAS_STAT_CHECK(cublasCtrsm(handle,side,uplo,trans,diag,m,n,&alph,d_ptrA,LDA,
                                        &bet,d_ptrB,LDB));
#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      tsc_end     = __rdtscp(&coreid);
      *tsc_delta  = tsc_end-tsc_start-rdtscp_cost;
      __asm__("lfence");
      printf("cublasCtrsm executed in: %llu reference cycles\n",*tsc_delta);
#endif  
      GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptrB,d_B,sizeof(cuComplex)*(LDB*n),
                                           cudaMemcpyDeviceToHost,stream));
      GMS_CUDA_DEBUG_CHECK(cudaStreamSynchronize(stream));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrA));
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptrB));
      GMS_CUBLAS_STAT_CHECK(cublasDestroy(handle));
      GMS_CUDA_DEBUG_CHECK(cudaStreamDestroy(stream));
      *ierr = 0;
      *err  = status;
      return; 
Error:   {  
           if(d_ptrA) cudaFree(d_ptrA);
           if(d_ptrC) cudaFree(d_ptrB);
           cublasDestroy(handle);
           cudaStreamDestroy(stream); 
           ierr = -2;
           *err = status;   
           return; 
    }  

} 
