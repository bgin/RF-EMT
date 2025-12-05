
#ifndef __GMS_AXPY_AVX2_UNROLL_10X_H__
#define __GMS_AXPY_AVX2_UNROLL_10X_H__

/*
     This version is based on BLIS framework original implementation.
     Few changes were introduced in order to more easily adapt it 
     to 'RF-EMT' project.
     
   The original authors license
   BLIS
   An object-based framework for developing high-performance BLAS-like
   libraries.
   Copyright (C) 2016 - 2019, Advanced Micro Devices, Inc.
   Copyright (C) 2018 - 2020, The University of Texas at Austin. All rights reserved.
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:
    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    - Neither the name(s) of the copyright holder(s) nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
     
*/

namespace file_version {

    const unsigned int gGMS_AXPY_AVX2_UNROLL10X_MAJOR = 1U;
    const unsigned int gGMS_AXPY_AVX2_UNROLL10X_MINOR = 0U;
    const unsigned int gGMS_AXPY_AVX2_UNROLL10X_MICRO = 0U;
    const unsigned int gGMS_AXPY_AVX2_UNROLL10X_FULLVER =
      1000U*gGMS_AXPY_AVX2_UNROLL10X_MAJOR+
      100U*gGMS_AXPY_AVX2_UNROLL10X_MINOR+
      10U*gGMS_AXPY_AVX2_UNROLL10X_MICRO;
    const char * const pgGMS_AXPY_AVX2_UNROLL10X_CREATION_DATE = "14-08-2021 10:39 AM +00200 (SAT 14 AUG 2021 GMT+2)";
    const char * const pgGMS_AXPY_AVX2_UNROLL10X_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const pgGMS_AXPY_AVX2_UNROLL10X_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_AXPY_AVX2_UNROLL10X_DESCRIPTION   = "AVX/AVX2 optimized AXPY kernels."

}

#include <cstdint>
#include <immintrin.h>
#include <omp.h>
#include "GMS_config.h"

namespace gms {

         namespace math {

              /*
                   SAXPY kernel unaligned.
           */
	      __ATTR_HOT__
	      __ATTR_ALIGN__(32)
	    
	      void saxpy_u_ymm8r4_unroll10x(const int32_t n,
	                                    const float alpha,
	                                    float * __restrict  x,
					    const int32_t incx,
					    float * __restrict  y,
					    const int32_t incy);
            /*
                 SAXPY kernel aligned
              */
           
	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	    void saxpy_a_ymm8r4_unroll10x(const int32_t n,
	                                    const float alpha,
	                                    float * __restrict __ATTR_ALIGN__(32) x,
					    const int32_t incx,
					    float * __restrict  __ATTR_ALIGN__(32) y,
					    const int32_t incy); 

	     /*
                 SAXPY kernel OpenMP parallelized (aligned version only)
              */
            
	     __ATTR_HOT__
	     __ATTR_ALIGN__(32)
	     void saxpy_a_ymm8r4_unroll10x_omp(const int32_t n,
	                                       const float alpha,
	                                       float * __restrict __ATTR_ALIGN__(32) x,
					       const int32_t incx,
					       float * __restrict  __ATTR_ALIGN__(32) y,
					       const int32_t incy); 

	      /*
                DAXPY kernel unaligned
              */
	   
	      __ATTR_HOT__
	      __ATTR_ALIGN__(32)
	     void daxpy_u_ymm4r8_unroll10x(const int32_t n,
	                                    const double alpha,
	                                    double * __restrict  x,
					    const int32_t incx,
					    double * __restrict  y,
					    const int32_t incy); 

	    /*
                DAXPY kernel aligned
            */
	   
	      __ATTR_HOT__
	      __ATTR_ALIGN__(32)
	      void daxpy_a_ymm4r8_unroll10x(const int32_t n,
	                                    const double alpha,
	                                    double * __restrict  __ATTR_ALIGN__(32) x,
					    const int32_t incx,
					    double * __restrict  __ATTR_ALIGN__(32) y,
					    const int32_t incy); 


	       /*
                DAXPY kernel OpenMP parallelized aligned
            */
	      
	      __ATTR_HOT__
	      __ATTR_ALIGN__(32)
	     void daxpy_a_ymm4r8_unroll10x_omp(const int32_t n,
	                                        const double alpha,
	                                        double * __restrict  __ATTR_ALIGN__(32) x,
					        const int32_t incx,
					        double * __restrict  __ATTR_ALIGN__(32) y,
					        const int32_t incy); 

	     


	     
	     
    } // math

} // gsm







#endif /*__GMS_AXPY_AVX2_UNROLL_10X_H__*/
