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
#include <cuda.h>
#include <cuda_runtime.h>
#include "GMS_cuda_memops.cuh"


//
//	Implementation
//



void copy_int32_cpu_to_gpu(int32_t * __restrict d_ptr, 
			   const int32_t * __restrict h_ptr,
			   const size_t n,
			   int32_t * ierr )
 {
	if(*ierr <= 0) *ierr = 0;
	if( NULL == h_ptr || n <= 0){ //  Host error handling
	   *ierr = -1;
	   return; 
	}
	cudaError_t status;

#if (GMS_CUDA_DEBUG_ON) == 1

	GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(int32_t)));
	GMS_CUDA_DEBUG_CHECK(cudaMemcpy((void*)&d_ptr[0],(void*)&h_ptr[0],n * sizeof(int32_t),cudaMemcpyHostToDevice));
#else
	GMS_CUDA_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(int32_t)));
	GMS_CUDA_CHECK(cudaMemcpy((void*)&d_ptr[0],(void*)&h_ptr[0],n*sizeof(int32_t),cudaMemcpyHostToDevice));
#endif
	*ierr = 0;
        return;
Error:
      
	if(__builtin_expect(status==cudaErrorMemoryAllocation,0)){
	   fatal_gpu_error(__PRETTY_FUNCTION__,
		            status);
        }
        else{
	   if(d_ptr) cudaFree(d_ptr);
	      *ierr = -2;
             return;
       } 
}






void copy_float_cpu_to_gpu(float * __restrict d_ptr,
			   const float * __restrict h_ptr,
			   const size_t n,
			   int32_t * ierr ) {
	if(*ierr <= 0) *ierr = 0;
	if(NULL == h_ptr || n <= 0){
	   *ierr = -1;
	   return;
	}
	cudaError_t status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(float)));
	GMS_CUDA_DEBUG_CHECK(cudaMemcpy((void*)&d_ptr[0],(void*)&h_ptr[0],n*sizeof(float),
					cudaMemcpyHostToDevice));
#else
	GMS_CUDA_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(float)));
	GMS_CUDA_CHECK(cudaMemcpy((void*)&d_ptr[0],(void*)&h_ptr[0],n*sizeof(float),
					cudaMemcpyHostToDevice));
#endif
	*ierr = 0;
        return;
Error:
	if(__builtin_expect(status==cudaErrorMemoryAllocation,0)){
	   fatal_gpu_error(__PRETTY_FUNCTION__,
		            status);
        }
         else {
	   if(d_ptr) cudaFree(d_ptr);
	      *ierr = -2;
             return;
       } 
}






void copy_double_cpu_to_gpu(double * __restrict d_ptr,
			    const double * __restrict h_ptr,
			    const size_t n,
			    int32_t * ierr ) {
	if(*ierr <= 0) *ierr = 0;
	if(NULL == h_ptr ||
	   0 >= n     ) {
	    *ierr = -1;
	    return;
	}
	cudaError_t status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(double)));
	GMS_CUDA_DEBUG_CHECK(cudaMemcpy((void*)&d_ptr[0],(void*)&h_ptr[0],n*sizeof(double),
					cudaMemcpyHostToDevice));
#else
	GMS_CUDA_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(double)));
	GMS_CUDA_CHECK(cudaMemcpy((void*)&d_ptr[0],(void*)&h_ptr[0],n*sizeof(double),
					cudaMemcpyHostToDevice));
#endif
	*ierr = 0;
        return;
Error:
      if(__builtin_expect(status==cudaErrorMemoryAllocation,0)){
	 fatal_gpu_error(__PRETTY_FUNCTION__,
		            status);
        }
        else {
	   if(d_ptr) cudaFree(d_ptr);
	      *ierr = -2;
             return;
       } 
}

void copy_complex4_cpu_to_gpu(cuComplex * __restrict d_ptr,
                              cuComplex * __restrict h_ptr,
                              const size_t n,
                              int32_t * ierr) {
       if(*ierr <= 0) *ierr = 0;
       if(NULL == h_ptr ||
	   0 >= n     ) {
	    *ierr = -1;
	    return;
	}
	cudaError_t status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(cuComplex)));
	GMS_CUDA_DEBUG_CHECK(cudaMemcpy((void*)&d_ptr[0],(void*)&h_ptr[0],n*sizeof(cuComplex),
					cudaMemcpyHostToDevice));
#else
	GMS_CUDA_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(cuComplex)));
	GMS_CUDA_CHECK(cudaMemcpy((void*)&d_ptr[0],(void*)&h_ptr[0],n*sizeof(cuComplex),
					cudaMemcpyHostToDevice));
#endif
	*ierr = 0;
        return;
Error:
        if(__builtin_expect(status==cudaErrorMemoryAllocation,0)){
	fatal_gpu_error(__PRETTY_FUNCTION__,
		            status);
        }
         else {
	   if(d_ptr) cudaFree(d_ptr);
	      *ierr = -2;
             return;
       } 
}



//
//	Allocate memory on GPU.
//

void alloc_int32_gpu(int32_t * __restrict d_ptr,
                     const size_t n,
		     int32_t * ierr ) {
	if(*ierr < 0) *ierr = 0;
	if(0 >= n) {

		*ierr = -1;
		return;
	}
	cudaError_t status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(int32_t)));
#else
	GMS_CUDA_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(int32_t)));
#endif
	*ierr = 0;
        return;
Error:
	if(__builtin_expect(status==cudaErrorMemoryAllocation,0)){
	   fatal_gpu_error(__PRETTY_FUNCTION__,
		            status);
        }
        else {
	   if(d_ptr) cudaFree(d_ptr);
	      *ierr = -2;
             return;
       } 
}







void alloc_float_gpu(float * __restrict d_ptr,
		     const size_t n,
		     int32_t * ierr ) {
	if(*ierr < 0) *ierr = 0;
	if(0 >= n ) {

		*ierr = -1;
		return;
	}
	cudaError_t status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(float)));
#else
	GMS_CUDA_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(float)));
#endif
	*ierr = 0;
        return;
Error:
	if(__builtin_expect(status==cudaErrorMemoryAllocation,0)){
	   fatal_gpu_error(__PRETTY_FUNCTION__,
		            status);
        }
         else {
	   if(d_ptr) cudaFree(d_ptr);
	      *ierr = -2;
             return;
       } 
}





void alloc_double_gpu(double * __restrict d_ptr,
		      const size_t n,
		      int32_t * ierr ) {
	if(*ierr < 0) *ierr = 0;
	if(0 >= n ) {

		*ierr = -1;
		return;
	}
	cudaError_t status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(double)));
#else
	GMS_CUDA_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(double)));
#endif
	*ierr = 0;
        return;
Error:
	
        if(__builtin_expect(status==cudaErrorMemoryAllocation,0)){
	  fatal_gpu_error(__PRETTY_FUNCTION__,
		            status);
          }
          else {
	      if(d_ptr) cudaFree(d_ptr);
	      *ierr = -2;
             return;
       } 
}

void alloc_complex4_gpu(cuComplex * __restrict d_ptr,
                        const size_t n,
                        int32_t * ierr) {
    if(*ierr < 0) *ierr = 0;
    if(0 >= n ) {

		*ierr = -1;
		return;
	}
	cudaError_t status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(cuComplex)));
#else
	GMS_CUDA_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(cuComplex)));
#endif
	*ierr = 0;
        return;
Error:
       if(__builtin_expect(status==cudaErrorMemoryAllocation,0)){
	  fatal_gpu_error(__PRETTY_FUNCTION__,
		            status);
        }
        else {
	   if(d_ptr) cudaFree(d_ptr);
	      *ierr = -2;
             return;
       } 
}



//
// GPU to CPU memory copy routines
//


//
// Copy array  of int32_t from GPU to CPU.
//

void copy_int32_gpu_to_cpu(const int32_t * __restrict d_ptr,
			   int32_t * __restrict h_ptr,
			   const int32_t n,
			   int32_t * ierr ) {
	if(*ierr < 0) *ierr = 0;
	if(0 >= n) {

		*ierr = -1;
		return;
	}
	cudaError_t status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMemcpy(&h_ptr[0],&d_ptr[0],n*sizeof(int32_t),
					cudaMemcpyDeviceToHost));
	
#else
	GMS_CUDA_CHECK(cudaMemcpy(&h_ptr[0],&d_ptr[0],n*sizeof(int32_t),
					cudaMemcpyDeviceToHost));
	
#endif
        *ierr = 0;
         return;
Error:
	
	*ierr = -2;
	return;
}



void copy_float_gpu_to_cpu(const float * __restrict d_ptr,
			   float * __restrict h_ptr,
			   const int32_t n,
			   int32_t * ierr ) {
	if(*ierr < 0) *ierr = 0;
	if(0 >= n) {

		*ierr = -1;
		return;
	}
	cudaError_t status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMemcpy(&h_ptr[0],&d_ptr[0],n*sizeof(float),
					cudaMemcpyDeviceToHost));
	
#else
	GMS_CUDA_CHECK(cudaMemcpy(&h_ptr[0],&d_ptr[0],n*sizeof(float),
					cudaMemcpyDeviceToHost));
	
#endif
        *ierr = 0;
        return;
Error:
	
	*ierr = -2;
	return;
}



void copy_double_gpu_to_cpu(const double * __restrict d_ptr,
			    double * __restrict h_ptr,
			    const int32_t n,
			    int32_t * ierr ) {
	if(*ierr < 0) *ierr = 0;
	if(0 >= n) {

		*ierr = -1;
		return;
	}
	cudaError_t status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMemcpy(&h_ptr[0],&d_ptr[0],n*sizeof(double),
					cudaMemcpyDeviceToHost));
	
#else
	GMS_CUDA_CHECK(cudaMemcpy(&h_ptr[0],&d_ptr[0],n*sizeof(double),
			   cudaMemcpyDeviceToHost));
	
#endif
        *ierr = 0;
Error:
	
	*ierr = -2;
	 return;
}


void copy_complex4_gpu_to_cpu(cuComplex * __restrict d_ptr,
                              cuComplex * __restrict h_ptr,
                              const size_t n,
                              int32_t ierr) {
        if(*ierr < 0) *ierr = 0;
	if(0 >= n) {

		*ierr = -1;
		return;
	}
	cudaError_t status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMemcpy(&h_ptr[0],&d_ptr[0],n*sizeof(cuComplex),
					cudaMemcpyDeviceToHost));
	
#else
	GMS_CUDA_CHECK(cudaMemcpy(&h_ptr[0],&d_ptr[0],n*sizeof(cuComplex),
			   cudaMemcpyDeviceToHost));
	
#endif
         *ierr = 0;
         return;
Error:
	
	*ierr = -2;
	 return;
} 

























