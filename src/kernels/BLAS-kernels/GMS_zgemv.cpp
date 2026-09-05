



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

#include "GMS_zgemv.h" // must define 'CONJ'



void zgemv_kernel_4x4(int32_t n,
                      double **__restrict ap,
		      double * __restrict x,
		      double * __restrict y) {
      register int32_t i = 0;
      if(((uintptr_t) & 0x1F) != 0ULL &&
         ((uintptr_t)y & 0x1F) != 0ULL) {

	    __asm__ __voltile__ (

                    "vzeroupper			 \n\t"

	            "vbroadcastsd	  (%2), %%ymm0                  \n\t"  // real part x0
	            "vbroadcastsd	 8(%2), %%ymm1                  \n\t"  // imag part x0
	            "vbroadcastsd	16(%2), %%ymm2                  \n\t"  // real part x1
	            "vbroadcastsd	24(%2), %%ymm3                  \n\t"  // imag part x1
	            "vbroadcastsd	32(%2), %%ymm4                  \n\t"  // real part x2
	            "vbroadcastsd	40(%2), %%ymm5                  \n\t"  // imag part x2
	            "vbroadcastsd	48(%2), %%ymm6                  \n\t"  // real part x3
	            "vbroadcastsd	56(%2), %%ymm7                  \n\t"  // imag part x3


	//	".align 16				        \n\t"
	           "1:				        \n\t"
#if (ZGEMV_KERNEL_USE_PREFETCT0) == 1
	           "prefetcht0      192(%4,%0,8)			\n\t"
#endif
	           "vmovups	(%4,%0,8), %%ymm8	        \n\t" // 2 complex values form a0
	           "vmovups      32(%4,%0,8), %%ymm9	        \n\t" // 2 complex values form a0
#if (ZGEMV_KERNEL_USE_PREFETCHT0) == 1
	           "prefetcht0      192(%5,%0,8)			\n\t"
#endif
	           "vmovups	(%5,%0,8), %%ymm10              \n\t" // 2 complex values form a1
	           "vmovups      32(%5,%0,8), %%ymm11              \n\t" // 2 complex values form a1

	           "vmulpd      %%ymm8 , %%ymm0, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	           "vmulpd      %%ymm8 , %%ymm1, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	           "vmulpd      %%ymm9 , %%ymm0, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	           "vmulpd      %%ymm9 , %%ymm1, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i
#if (ZGEMV_KERNEL_USE_PREFETCHT0) == 1
	           "prefetcht0      192(%6,%0,8)			\n\t"
#endif
	           "vmovups	(%6,%0,8), %%ymm8	        \n\t" // 2 complex values form a2
	           "vmovups      32(%6,%0,8), %%ymm9	        \n\t" // 2 complex values form a2

	           "vfmadd231pd      %%ymm10, %%ymm2, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	           "vfmadd231pd      %%ymm10, %%ymm3, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	           "vfmadd231pd      %%ymm11, %%ymm2, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	           "vfmadd231pd      %%ymm11, %%ymm3, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i
#if (ZGEMV_KERNEL_USE_PREFETCHT0) == 1
	           "prefetcht0      192(%7,%0,8)			\n\t"
#endif
	           "vmovups	(%7,%0,8), %%ymm10              \n\t" // 2 complex values form a3
	           "vmovups      32(%7,%0,8), %%ymm11              \n\t" // 2 complex values form a3

	           "vfmadd231pd      %%ymm8 , %%ymm4, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	           "vfmadd231pd      %%ymm8 , %%ymm5, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	           "vfmadd231pd      %%ymm9 , %%ymm4, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	           "vfmadd231pd      %%ymm9 , %%ymm5, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i

	           "vfmadd231pd      %%ymm10, %%ymm6, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	           "vfmadd231pd      %%ymm10, %%ymm7, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	           "vfmadd231pd      %%ymm11, %%ymm6, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	           "vfmadd231pd      %%ymm11, %%ymm7, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i
#if (ZGEMV_KERNEL_USE_PREFETCHT0) == 1
	           "prefetcht0      192(%3,%0,8)			\n\t"
#endif
	           "vmovups	  (%3,%0,8),  %%ymm10           \n\t"
	           "vmovups	32(%3,%0,8),  %%ymm11           \n\t"

#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
                  "vpermilpd      $0x5 , %%ymm13, %%ymm13               \n\t"
                  "vpermilpd      $0x5 , %%ymm15, %%ymm15               \n\t"
                  "vaddsubpd      %%ymm13, %%ymm12, %%ymm8              \n\t"
                  "vaddsubpd      %%ymm15, %%ymm14, %%ymm9              \n\t"
#else
                  "vpermilpd      $0x5 , %%ymm12, %%ymm12               \n\t"
                  "vpermilpd      $0x5 , %%ymm14, %%ymm14               \n\t"
                  "vaddsubpd      %%ymm12, %%ymm13, %%ymm8              \n\t"
                  "vaddsubpd      %%ymm14, %%ymm15, %%ymm9              \n\t"
                  "vpermilpd      $0x5 , %%ymm8 , %%ymm8                \n\t"
                  "vpermilpd      $0x5 , %%ymm9 , %%ymm9                \n\t"
#endif

                  "vaddpd         %%ymm8, %%ymm10, %%ymm12              \n\t"
                  "vaddpd         %%ymm9, %%ymm11, %%ymm13              \n\t"

	          "vmovups  %%ymm12,   (%3,%0,8)		        \n\t" // 2 complex values to y	
	          "vmovups  %%ymm13, 32(%3,%0,8)		        \n\t"	

                  "addq		$8 , %0	  	 	        \n\t"
	          "subq	        $4 , %1			        \n\t"		
	          "jnz		1b		        \n\t"
	          "vzeroupper			 \n\t"

	:
          "+r" (i),	// 0	
	  "+r" (n)  	// 1
	:
          "r" (x),      // 2
          "r" (y),      // 3
          "r" (ap[0]),  // 4
          "r" (ap[1]),  // 5
          "r" (ap[2]),  // 6
          "r" (ap[3])   // 7
	: "cc", 
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	   );
      }
       else {

             __asm__ __volatile__ (
              "vzeroupper			 \n\t"

	      "vbroadcastsd	  (%2), %%ymm0                  \n\t"  // real part x0
	      "vbroadcastsd	 8(%2), %%ymm1                  \n\t"  // imag part x0
	      "vbroadcastsd	16(%2), %%ymm2                  \n\t"  // real part x1
	      "vbroadcastsd	24(%2), %%ymm3                  \n\t"  // imag part x1
	      "vbroadcastsd	32(%2), %%ymm4                  \n\t"  // real part x2
	      "vbroadcastsd	40(%2), %%ymm5                  \n\t"  // imag part x2
	      "vbroadcastsd	48(%2), %%ymm6                  \n\t"  // real part x3
	      "vbroadcastsd	56(%2), %%ymm7                  \n\t"  // imag part x3


	//	".align 16				        \n\t"
	     "1:				        \n\t"
#if (ZGEMV_KERNEL_USE_PREFETCT0) == 1
	     "prefetcht0      192(%4,%0,8)			\n\t"
#endif
	     "vmovaps	(%4,%0,8), %%ymm8	        \n\t" // 2 complex values form a0
	     "vmovaps      32(%4,%0,8), %%ymm9	        \n\t" // 2 complex values form a0
#if (ZGEMV_KERNEL_USE_PREFETCT0) == 1
	     "prefetcht0      192(%5,%0,8)			\n\t"
#endif
	     "vmovaps	(%5,%0,8), %%ymm10              \n\t" // 2 complex values form a1
	     "vmovaps      32(%5,%0,8), %%ymm11              \n\t" // 2 complex values form a1

	     "vmulpd      %%ymm8 , %%ymm0, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	     "vmulpd      %%ymm8 , %%ymm1, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	     "vmulpd      %%ymm9 , %%ymm0, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	     "vmulpd      %%ymm9 , %%ymm1, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i
#if (ZGEMV_KERNEL_USE_PREFETCT0) == 1
	     "prefetcht0      192(%6,%0,8)			\n\t"
#endif
	     "vmovaps	(%6,%0,8), %%ymm8	        \n\t" // 2 complex values form a2
	     "vmovaps      32(%6,%0,8), %%ymm9	        \n\t" // 2 complex values form a2

	     "vfmadd231pd      %%ymm10, %%ymm2, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	     "vfmadd231pd      %%ymm10, %%ymm3, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	     "vfmadd231pd      %%ymm11, %%ymm2, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	     "vfmadd231pd      %%ymm11, %%ymm3, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i
#if (ZGEMV_KERNEL_USE_PREFETCT0) == 1
	     "prefetcht0      192(%7,%0,8)			\n\t"
#endif
	     "vmovaps	(%7,%0,8), %%ymm10              \n\t" // 2 complex values form a3
	     "vmovaps      32(%7,%0,8), %%ymm11              \n\t" // 2 complex values form a3

	     "vfmadd231pd      %%ymm8 , %%ymm4, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	     "vfmadd231pd      %%ymm8 , %%ymm5, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	     "vfmadd231pd      %%ymm9 , %%ymm4, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	     "vfmadd231pd      %%ymm9 , %%ymm5, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i

	     "vfmadd231pd      %%ymm10, %%ymm6, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	     "vfmadd231pd      %%ymm10, %%ymm7, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	     "vfmadd231pd      %%ymm11, %%ymm6, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	     "vfmadd231pd      %%ymm11, %%ymm7, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i
#if (ZGEMV_KERNEL_USE_PREFETCT0) == 1
	     "prefetcht0      192(%3,%0,8)			\n\t"
#endif
	     "vmovaps	  (%3,%0,8),  %%ymm10           \n\t"
	     "vmovaps	32(%3,%0,8),  %%ymm11           \n\t"

#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
             "vpermilpd      $0x5 , %%ymm13, %%ymm13               \n\t"
             "vpermilpd      $0x5 , %%ymm15, %%ymm15               \n\t"
             "vaddsubpd      %%ymm13, %%ymm12, %%ymm8              \n\t"
             "vaddsubpd      %%ymm15, %%ymm14, %%ymm9              \n\t"
#else
             "vpermilpd      $0x5 , %%ymm12, %%ymm12               \n\t"
             "vpermilpd      $0x5 , %%ymm14, %%ymm14               \n\t"
             "vaddsubpd      %%ymm12, %%ymm13, %%ymm8              \n\t"
             "vaddsubpd      %%ymm14, %%ymm15, %%ymm9              \n\t"
             "vpermilpd      $0x5 , %%ymm8 , %%ymm8                \n\t"
             "vpermilpd      $0x5 , %%ymm9 , %%ymm9                \n\t"
#endif

             "vaddpd         %%ymm8, %%ymm10, %%ymm12              \n\t"
             "vaddpd         %%ymm9, %%ymm11, %%ymm13              \n\t"

	      "vmovaps  %%ymm12,   (%3,%0,8)		        \n\t" // 2 complex values to y	
	      "vmovaps  %%ymm13, 32(%3,%0,8)		        \n\t"	

            "addq		$8 , %0	  	 	        \n\t"
	    "subq	        $4 , %1			        \n\t"		
	    "jnz		1b		        \n\t"
	    "vzeroupper			 \n\t"

	:
          "+r" (i),	// 0	
	  "+r" (n)  	// 1
	:
          "r" (x),      // 2
          "r" (y),      // 3
          "r" (ap[0]),  // 4
          "r" (ap[1]),  // 5
          "r" (ap[2]),  // 6
          "r" (ap[3])   // 7
	: "cc", 
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
          );
       }

}



void zgemv_kernel_4x2(const int32_t n,
                      double ** __restrict ap,
		      double *  __restrict x,
		      double *  __restrict y) {
         register int32_t i = 0;
         if(((uintptr_t) & 0x1F) != 0ULL &&
            ((uintptr_t)y & 0x1F) != 0ULL) {

	       __asm__ __volatile__ (

	              	"vzeroupper			 \n\t"

	                "vbroadcastsd	  (%2), %%ymm0                  \n\t"  // real part x0
	                "vbroadcastsd	 8(%2), %%ymm1                  \n\t"  // imag part x0
	                "vbroadcastsd	16(%2), %%ymm2                  \n\t"  // real part x1
	                "vbroadcastsd	24(%2), %%ymm3                  \n\t"  // imag part x1


	//	".align 16				        \n\t"
	               "1:				        \n\t"
#if (ZGEMV_KERNEL_USE_PREFETCT0) == 1
	               "prefetcht0      192(%4,%0,8)			\n\t"
#endif
	               "vmovups	(%4,%0,8), %%ymm8	        \n\t" // 2 complex values form a0
	               "vmovups      32(%4,%0,8), %%ymm9	        \n\t" // 2 complex values form a0
#if (ZGEMV_KERNEL_USE_PREFETCT0) == 1
	               "prefetcht0      192(%5,%0,8)			\n\t"
#endif
	               "vmovups	(%5,%0,8), %%ymm10              \n\t" // 2 complex values form a1
	               "vmovups      32(%5,%0,8), %%ymm11              \n\t" // 2 complex values form a1

	               "vmulpd      %%ymm8 , %%ymm0, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	               "vmulpd      %%ymm8 , %%ymm1, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	               "vmulpd      %%ymm9 , %%ymm0, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	               "vmulpd      %%ymm9 , %%ymm1, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i

	               "vfmadd231pd      %%ymm10, %%ymm2, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	               "vfmadd231pd      %%ymm10, %%ymm3, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	               "vfmadd231pd      %%ymm11, %%ymm2, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	               "vfmadd231pd      %%ymm11, %%ymm3, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i

#if (ZGEMV_KERNEL_USE_PREFETCT0) == 1
	               "prefetcht0      192(%3,%0,8)			\n\t"
#endif
	               "vmovups	  (%3,%0,8),  %%ymm10           \n\t"
	               "vmovups	32(%3,%0,8),  %%ymm11           \n\t"

#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
                       "vpermilpd      $0x5 , %%ymm13, %%ymm13               \n\t"
                       "vpermilpd      $0x5 , %%ymm15, %%ymm15               \n\t"
                       "vaddsubpd      %%ymm13, %%ymm12, %%ymm8              \n\t"
                       "vaddsubpd      %%ymm15, %%ymm14, %%ymm9              \n\t"
#else
                       "vpermilpd      $0x5 , %%ymm12, %%ymm12               \n\t"
                       "vpermilpd      $0x5 , %%ymm14, %%ymm14               \n\t"
                       "vaddsubpd      %%ymm12, %%ymm13, %%ymm8              \n\t"
                       "vaddsubpd      %%ymm14, %%ymm15, %%ymm9              \n\t"
                       "vpermilpd      $0x5 , %%ymm8 , %%ymm8                \n\t"
                       "vpermilpd      $0x5 , %%ymm9 , %%ymm9                \n\t"
#endif

                       "vaddpd         %%ymm8, %%ymm10, %%ymm12              \n\t"
                       "vaddpd         %%ymm9, %%ymm11, %%ymm13              \n\t"

	               "vmovups  %%ymm12,   (%3,%0,8)		        \n\t" // 2 complex values to y	
	               "vmovups  %%ymm13, 32(%3,%0,8)		        \n\t"	

                       "addq		$8 , %0	  	 	        \n\t"
	               "subq	        $4 , %1			        \n\t"		
	               "jnz		1b		        \n\t"
	               "vzeroupper			 \n\t"

	               :
                       "+r" (i),	// 0	
	               "+r" (n)  	// 1
	               :
                       "r" (x),      // 2
                       "r" (y),      // 3
                       "r" (ap[0]),  // 4
                       "r" (ap[1])   // 5
	              : "cc", 
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"

	       );
       }
       else {
                  __asm__  __volatile__ (
               	         "vzeroupper			 \n\t"

	                "vbroadcastsd	  (%2), %%ymm0                  \n\t"  // real part x0
	                "vbroadcastsd	 8(%2), %%ymm1                  \n\t"  // imag part x0
	                "vbroadcastsd	16(%2), %%ymm2                  \n\t"  // real part x1
	                "vbroadcastsd	24(%2), %%ymm3                  \n\t"  // imag part x1


	//	".align 16				        \n\t"
	               "1:				        \n\t"
#if (ZGEMV_KERNEL_USE_PREFETCT0) == 1
	               "prefetcht0      192(%4,%0,8)			\n\t"
#endif
	               "vmovaps	(%4,%0,8), %%ymm8	        \n\t" // 2 complex values form a0
	               "vmovaps      32(%4,%0,8), %%ymm9	        \n\t" // 2 complex values form a0
#if (ZGEMV_KERNEL_USE_PREFETCT0) == 1
	               "prefetcht0      192(%5,%0,8)			\n\t"
#endif
	               "vmovaps	(%5,%0,8), %%ymm10              \n\t" // 2 complex values form a1
	               "vmovaps      32(%5,%0,8), %%ymm11              \n\t" // 2 complex values form a1

	               "vmulpd      %%ymm8 , %%ymm0, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	               "vmulpd      %%ymm8 , %%ymm1, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	               "vmulpd      %%ymm9 , %%ymm0, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	               "vmulpd      %%ymm9 , %%ymm1, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i

	               "vfmadd231pd      %%ymm10, %%ymm2, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	               "vfmadd231pd      %%ymm10, %%ymm3, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	               "vfmadd231pd      %%ymm11, %%ymm2, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	               "vfmadd231pd      %%ymm11, %%ymm3, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i

#if (ZGEMV_KERNEL_USE_PREFETCHT0) == 1
	               "prefetcht0      192(%3,%0,8)			\n\t"
#endif
	               "vmovaps	  (%3,%0,8),  %%ymm10           \n\t"
	               "vmovaps	32(%3,%0,8),  %%ymm11           \n\t"

#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
                       "vpermilpd      $0x5 , %%ymm13, %%ymm13               \n\t"
                       "vpermilpd      $0x5 , %%ymm15, %%ymm15               \n\t"
                       "vaddsubpd      %%ymm13, %%ymm12, %%ymm8              \n\t"
                       "vaddsubpd      %%ymm15, %%ymm14, %%ymm9              \n\t"
#else
                       "vpermilpd      $0x5 , %%ymm12, %%ymm12               \n\t"
                       "vpermilpd      $0x5 , %%ymm14, %%ymm14               \n\t"
                       "vaddsubpd      %%ymm12, %%ymm13, %%ymm8              \n\t"
                       "vaddsubpd      %%ymm14, %%ymm15, %%ymm9              \n\t"
                       "vpermilpd      $0x5 , %%ymm8 , %%ymm8                \n\t"
                       "vpermilpd      $0x5 , %%ymm9 , %%ymm9                \n\t"
#endif

                       "vaddpd         %%ymm8, %%ymm10, %%ymm12              \n\t"
                       "vaddpd         %%ymm9, %%ymm11, %%ymm13              \n\t"

	               "vmovaps  %%ymm12,   (%3,%0,8)		        \n\t" // 2 complex values to y	
	               "vmovaps  %%ymm13, 32(%3,%0,8)		        \n\t"	

                       "addq		$8 , %0	  	 	        \n\t"
	               "subq	        $4 , %1			        \n\t"		
	               "jnz		1b		        \n\t"
	               "vzeroupper			 \n\t"

	               :
                       "+r" (i),	// 0	
	               "+r" (n)  	// 1
	               :
                       "r" (x),      // 2
                       "r" (y),      // 3
                       "r" (ap[0]),  // 4
                       "r" (ap[1])   // 5
	              : "cc", 
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
           );
       }
   
}



void zgemv_kernel_4x1(const int32_t n,
                      double * __restrict ap,
		      double * __restrict x,
		      double * __restrict y) {
         register int32_t i = 0;
	 if(((uintptr_t) & 0x1F) != 0ULL &&
            ((uintptr_t)y & 0x1F) != 0ULL) {

	      __asm__ __volatile__ (

	           	"vzeroupper			 \n\t"

	                "vbroadcastsd	  (%2), %%ymm0                  \n\t"  // real part x0
	                "vbroadcastsd	 8(%2), %%ymm1                  \n\t"  // imag part x0

	 	        ".align 16				        \n\t"
	                "1:				        \n\t"
#if (ZGEMV_KERNEL_USE_PREFETCHT0) == 1			
	                "prefetcht0      192(%4,%0,8)			\n\t"
#endif
	                "vmovaps	(%4,%0,8), %%ymm8	        \n\t" // 2 complex values form a0
	                "vmovaps      32(%4,%0,8), %%ymm9	        \n\t" // 2 complex values form a0

	                "vmulpd      %%ymm8 , %%ymm0, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	                "vmulpd      %%ymm8 , %%ymm1, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	                "vmulpd      %%ymm9 , %%ymm0, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	                "vmulpd      %%ymm9 , %%ymm1, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i
#if (ZGEMV_KERNEL_USE_PREFETCHT0) == 1
	                "prefetcht0      192(%3,%0,8)			\n\t"
#endif
	                "vmovaps	  (%3,%0,8),  %%ymm10           \n\t"
	                "vmovaps	32(%3,%0,8),  %%ymm11           \n\t"

#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
                        "vpermilpd      $0x5 , %%ymm13, %%ymm13               \n\t"
                        "vpermilpd      $0x5 , %%ymm15, %%ymm15               \n\t"
                        "vaddsubpd      %%ymm13, %%ymm12, %%ymm8              \n\t"
                        "vaddsubpd      %%ymm15, %%ymm14, %%ymm9              \n\t"
#else
                        "vpermilpd      $0x5 , %%ymm12, %%ymm12               \n\t"
                        "vpermilpd      $0x5 , %%ymm14, %%ymm14               \n\t"
                        "vaddsubpd      %%ymm12, %%ymm13, %%ymm8              \n\t"
                        "vaddsubpd      %%ymm14, %%ymm15, %%ymm9              \n\t"
                        "vpermilpd      $0x5 , %%ymm8 , %%ymm8                \n\t"
                        "vpermilpd      $0x5 , %%ymm9 , %%ymm9                \n\t"
#endif

                        "vaddpd         %%ymm8, %%ymm10, %%ymm12              \n\t"
                        "vaddpd         %%ymm9, %%ymm11, %%ymm13              \n\t"

	                "vmovaps  %%ymm12,   (%3,%0,8)		        \n\t" // 2 complex values to y	
	                "vmovaps  %%ymm13, 32(%3,%0,8)		        \n\t"	

                        "addq		$8 , %0	  	 	        \n\t"
	                "subq	        $4 , %1			        \n\t"		
	                "jnz		1b		        \n\t"
	                "vzeroupper			 \n\t"

	              :
                      "+r" (i),	// 0	
	              "+r" (n)  	// 1
	              :
                      "r" (x),      // 2
                      "r" (y),      // 3
                      "r" (ap)      // 4
	              : "cc", 
	              "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	               "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	              "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	              "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	              "memory"
	           );
	 } else {

                  __asm__ __volatile__ (

	           	"vzeroupper			 \n\t"

	                "vbroadcastsd	  (%2), %%ymm0                  \n\t"  // real part x0
	                "vbroadcastsd	 8(%2), %%ymm1                  \n\t"  // imag part x0

	 	        ".align 16				        \n\t"
	                "1:				        \n\t"
#if (ZGEMV_KERNEL_USE_PREFETCHT0) == 1			
	                "prefetcht0      192(%4,%0,8)			\n\t"
#endif
	                "vmovaps	(%4,%0,8), %%ymm8	        \n\t" // 2 complex values form a0
	                "vmovaps      32(%4,%0,8), %%ymm9	        \n\t" // 2 complex values form a0

	                "vmulpd      %%ymm8 , %%ymm0, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	                "vmulpd      %%ymm8 , %%ymm1, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	                "vmulpd      %%ymm9 , %%ymm0, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	                "vmulpd      %%ymm9 , %%ymm1, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i
#if (ZGEMV_KERNEL_USE_PREFETCHT0) == 1
	                "prefetcht0      192(%3,%0,8)			\n\t"
#endif
	                "vmovaps	  (%3,%0,8),  %%ymm10           \n\t"
	                "vmovaps	32(%3,%0,8),  %%ymm11           \n\t"

#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
                        "vpermilpd      $0x5 , %%ymm13, %%ymm13               \n\t"
                        "vpermilpd      $0x5 , %%ymm15, %%ymm15               \n\t"
                        "vaddsubpd      %%ymm13, %%ymm12, %%ymm8              \n\t"
                        "vaddsubpd      %%ymm15, %%ymm14, %%ymm9              \n\t"
#else
                        "vpermilpd      $0x5 , %%ymm12, %%ymm12               \n\t"
                        "vpermilpd      $0x5 , %%ymm14, %%ymm14               \n\t"
                        "vaddsubpd      %%ymm12, %%ymm13, %%ymm8              \n\t"
                        "vaddsubpd      %%ymm14, %%ymm15, %%ymm9              \n\t"
                        "vpermilpd      $0x5 , %%ymm8 , %%ymm8                \n\t"
                        "vpermilpd      $0x5 , %%ymm9 , %%ymm9                \n\t"
#endif

                        "vaddpd         %%ymm8, %%ymm10, %%ymm12              \n\t"
                        "vaddpd         %%ymm9, %%ymm11, %%ymm13              \n\t"

	                "vmovaps  %%ymm12,   (%3,%0,8)		        \n\t" // 2 complex values to y	
	                "vmovaps  %%ymm13, 32(%3,%0,8)		        \n\t"	

                        "addq		$8 , %0	  	 	        \n\t"
	                "subq	        $4 , %1			        \n\t"		
	                "jnz		1b		        \n\t"
	                "vzeroupper			 \n\t"

	              :
                      "+r" (i),	// 0	
	              "+r" (n)  	// 1
	              :
                      "r" (x),      // 2
                      "r" (y),      // 3
                      "r" (ap)      // 4
	              : "cc", 
	              "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	               "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	              "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	              "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	              "memory"
	           );
	 }
	  

}



void add_y(const int32_t n,
           double * __restrict src,
	   double * __restrict dest,
	   const int32_t inc_dest,
	   const double alpha_r,
	   const double alpha_i) {

        register int32_t i;
	if(inc_dest != 2) {
            double temp_r = 0.0;
	    double temp_i = 0.0;
	    for(i = 0; i != n; ++i) {
#if !defined(XCONJ) 
			temp_r = alpha_r * src[0] - alpha_i * src[1];
			temp_i = alpha_r * src[1] + alpha_i * src[0];
#else
			temp_r =  alpha_r * src[0] + alpha_i * src[1];
			temp_i = -alpha_r * src[1] + alpha_i * src[0];
#endif
                        *dest += temp_r;
			*(dest+1) += temp_i;
			src += 2;
			dest += inc_dest;
	    }
	    return;
	}
	i = 0;
	const bool isAligned32 = false;
	isAligned32 = ((uintptr_t) & 0x1F) != 0ULL && ((uintptr_t)y & 0x1F) != 0ULL;
        if(!isAligned32) {

	    	__asm__  __volatile__ (
	

	              "vzeroupper			 \n\t"

	              "vbroadcastsd	  (%4), %%ymm0                  \n\t"  // alpha_r
	              "vbroadcastsd	  (%5), %%ymm1                  \n\t"  // alpha_i

		".align 16				        \n\t"
	        "1:                                              \n\t"
#if (ZGEMV_KERNEL_USE_PREFETCHT0) == 1					
	        "prefetcht0      192(%2,%0,8)			\n\t"
#endif
	        "vmovups	(%2,%0,8), %%ymm8	        \n\t" // 2 complex values from src
	        "vmovups      32(%2,%0,8), %%ymm9	        \n\t" 

	        "vmulpd      %%ymm8 , %%ymm0, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	        "vmulpd      %%ymm8 , %%ymm1, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	        "vmulpd      %%ymm9 , %%ymm0, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	        "vmulpd      %%ymm9 , %%ymm1, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i
#if (ZGEMV_KERNEL_USE_PREFETCHT0) == 1	
	        "prefetcht0      192(%3,%0,8)			\n\t"
#endif
	        "vmovups	  (%3,%0,8),  %%ymm10           \n\t" // 2 complex values from dest
	        "vmovups	32(%3,%0,8),  %%ymm11           \n\t"

#if !defined(XCONJ)
               "vpermilpd      $0x5 , %%ymm13, %%ymm13               \n\t"
               "vpermilpd      $0x5 , %%ymm15, %%ymm15               \n\t"
               "vaddsubpd      %%ymm13, %%ymm12, %%ymm8              \n\t"
               "vaddsubpd      %%ymm15, %%ymm14, %%ymm9              \n\t"
#else
               "vpermilpd      $0x5 , %%ymm12, %%ymm12               \n\t"
               "vpermilpd      $0x5 , %%ymm14, %%ymm14               \n\t"
               "vaddsubpd      %%ymm12, %%ymm13, %%ymm8              \n\t"
               "vaddsubpd      %%ymm14, %%ymm15, %%ymm9              \n\t"
                "vpermilpd      $0x5 , %%ymm8 , %%ymm8                \n\t"
                "vpermilpd      $0x5 , %%ymm9 , %%ymm9                \n\t"
#endif

              "vaddpd         %%ymm8, %%ymm10, %%ymm12              \n\t"
              "vaddpd         %%ymm9, %%ymm11, %%ymm13              \n\t"

	      "vmovups  %%ymm12,   (%3,%0,8)		        \n\t" // 2 complex values to y	
	      "vmovups  %%ymm13, 32(%3,%0,8)		        \n\t"	

              "addq		$8 , %0	  	 	        \n\t"
	      "subq	        $4 , %1			        \n\t"		
	      "jnz		1b		        \n\t"
	      "vzeroupper			 \n\t"

	      :
              "+r" (i),	      // 0	
	      "+r" (n)  	      // 1
	      :
              "r" (src),          // 2
              "r" (dest),         // 3
              "r" (&alpha_r),     // 4
              "r" (&alpha_i)      // 5
	      : "cc", 
	      "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	      "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	      "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	      "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	      "memory"
	     );

	}
	 else {

	              "vzeroupper			 \n\t"

	              "vbroadcastsd	  (%4), %%ymm0                  \n\t"  // alpha_r
	              "vbroadcastsd	  (%5), %%ymm1                  \n\t"  // alpha_i

		".align 16				        \n\t"
	        "1:                                              \n\t"
#if (ZGEMV_KERNEL_USE_PREFETCHT0) == 1					
	        "prefetcht0      192(%2,%0,8)			\n\t"
#endif
	        "vmovaps	(%2,%0,8), %%ymm8	        \n\t" // 2 complex values from src
	        "vmovaps      32(%2,%0,8), %%ymm9	        \n\t" 

	        "vmulpd      %%ymm8 , %%ymm0, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	        "vmulpd      %%ymm8 , %%ymm1, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	        "vmulpd      %%ymm9 , %%ymm0, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	        "vmulpd      %%ymm9 , %%ymm1, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i
#if (ZGEMV_KERNEL_USE_PREFETCHT0) == 1	
	        "prefetcht0      192(%3,%0,8)			\n\t"
#endif
	        "vmovaps	  (%3,%0,8),  %%ymm10           \n\t" // 2 complex values from dest
	        "vmovaps	32(%3,%0,8),  %%ymm11           \n\t"

#if !defined(XCONJ)
               "vpermilpd      $0x5 , %%ymm13, %%ymm13               \n\t"
               "vpermilpd      $0x5 , %%ymm15, %%ymm15               \n\t"
               "vaddsubpd      %%ymm13, %%ymm12, %%ymm8              \n\t"
               "vaddsubpd      %%ymm15, %%ymm14, %%ymm9              \n\t"
#else
               "vpermilpd      $0x5 , %%ymm12, %%ymm12               \n\t"
               "vpermilpd      $0x5 , %%ymm14, %%ymm14               \n\t"
               "vaddsubpd      %%ymm12, %%ymm13, %%ymm8              \n\t"
               "vaddsubpd      %%ymm14, %%ymm15, %%ymm9              \n\t"
                "vpermilpd      $0x5 , %%ymm8 , %%ymm8                \n\t"
                "vpermilpd      $0x5 , %%ymm9 , %%ymm9                \n\t"
#endif

              "vaddpd         %%ymm8, %%ymm10, %%ymm12              \n\t"
              "vaddpd         %%ymm9, %%ymm11, %%ymm13              \n\t"

	      "vmovaps  %%ymm12,   (%3,%0,8)		        \n\t" // 2 complex values to y	
	      "vmovaps  %%ymm13, 32(%3,%0,8)		        \n\t"	

              "addq		$8 , %0	  	 	        \n\t"
	      "subq	        $4 , %1			        \n\t"		
	      "jnz		1b		        \n\t"
	      "vzeroupper			 \n\t"

	      :
              "+r" (i),	      // 0	
	      "+r" (n)  	      // 1
	      :
              "r" (src),          // 2
              "r" (dest),         // 3
              "r" (&alpha_r),     // 4
              "r" (&alpha_i)      // 5
	      : "cc", 
	      "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	      "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	      "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	      "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	      "memory"
	     );
	}
}

#include <cstring>
#if defined __GNUC__ && !defined __INTEL_COMPILER
#include <omp.h>
#endif

#define ZGEMV_NBMAX 1024


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
	   double * buffer) {

      if(m<1 || n<1) return(0);
      __ATTR_ALIGN__(32) double xbuffer[8] = {};
      __ATTR_ALIGN__(32) double *ap[4] = {};
      double * __restrict ybuffer = NULL;
      double * __restrict a_ptr = NULL;
      double * __restrict x_ptr = NULL;
      double * __restrict y_ptr = NULL;
      int32_t n1,m1,m2,m3,n2,lda4,i,t0,nb;
      char pad[18] = {};

      ybuffer = buffer;
      inc_x *= 2;
      inc_y *= 2;
      lda   *= 2;
      lda4   = 4*lda;
      n1 = n/4;
      n2 = n%4;
      t0 = m%4;
      m3 = t0;
      m1 = m-t0;
      m2 = (m%ZGEMV_NBMAX)-t0;
      y_ptr = y;
      nb = NBMAX;
      
      while(nb == NBMAX) {
           m1 -= nb;
	   if(m1<0) {
	      if(m2==0) break;
	      nb = m2;
	   }
	   a_ptr = a;
	   ap[0] = a_ptr;
	   ap[1] = a_ptr+lda;
	   ap[2] = ap[1]+lda;
	   ap[3] = ap[2]+lda;
	   x_ptr = x;
	   memset(ybuffer,0,nb*16);

	   if(inc_x == 2) {

	       for(i = 0; i != n1; ++i) {
                   zgemv_kernel_4x4(nb,ap,x_ptr,ybuffer);
		   ap[0] += lda4;
		   ap[1] += lda4;
		   ap[2] += lda4;
		   ap[3] += lda4;
		   a_ptr += lda4;
		   x_ptr += 8;
	       }

	       if(n2 & 2) {
                  zgemv_kernel_4x2(nb,ap,x_ptr,ybuffer);
		  x_ptr += 4;
		  a_ptr += (lda + lda);
	       }

	       if(n2 & 1) {
                  zgemv_kernel_4x1(nb,a_ptr,x_ptr,ybuffer);
	       }
	   }
	   else {

	          for(i = 0; i != n1; ++i) {
                      xbuffer[0] = x_ptr[0];
		      xbuffer[1] = x_ptr[1];
		      x_ptr += inc_x;
		      xbuffer[2] = x_ptr[0];
		      xbuffer[3] = x_ptr[1];
		      x_ptr += inc_x;	
		      xbuffer[4] = x_ptr[0];
		      xbuffer[5] = x_ptr[1];
		      x_ptr += inc_x;	
		      xbuffer[6] = x_ptr[0];
		      xbuffer[7] = x_ptr[1];
		      x_ptr += inc_x;
		      zgemv_kernel_4x4(nb,ap_xbuffer,ybuffer);
		      ap[0] += lda4;
		      ap[1] += lda4;
		      ap[2] += lda4;
		      ap[3] += lda4;
		      a_ptr += lda4;
		  }
		  for(i = 0; i != n2; ++i) {
                      xbuffer[0] = x_ptr[0];
		      xbuffer[1] = x_ptr[1];
		      x_ptr += inc_x;
		      zgemv_kernel_4x1(nb,a_ptr,xbuffer,ybuffer);
		      a_ptr += lda;
		  }
	   }

	   add_y(nb,ybuffer,y_ptr,inc_y,alpha_r,alpha_i);
	   a  += 2*nb;
	   y_ptr += nb*inc_y;
      }

      if(m3 == 0) return(0);
      if(m3 == 1) {
         a_ptr = a;
	 x_ptr = x;
	 double temp_r = 0.0;
	 double temp_i = 0.0;
	 if(lda==2 && inc_x==2) {
#if (GMS_ZGEMV_EXPLICITLY_VECTORIZE) == 1
#if defined __INTEL_COMPILER || defined __ICC
#pragma vectorlengthfor(double) reduction(+:temp_r,temp_i)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd reduction(+:temp_r,temp_i)
#endif
#endif
            for(i = 0; i != (n & -2); i += 2) {
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
		    temp_r += a_ptr[0] * x_ptr[0] - a_ptr[1] * x_ptr[1];
		    temp_i += a_ptr[0] * x_ptr[1] + a_ptr[1] * x_ptr[0];
		    temp_r += a_ptr[2] * x_ptr[2] - a_ptr[3] * x_ptr[3];
		    temp_i += a_ptr[2] * x_ptr[3] + a_ptr[3] * x_ptr[2];
#else
		    temp_r += a_ptr[0] * x_ptr[0] + a_ptr[1] * x_ptr[1];
		    temp_i += a_ptr[0] * x_ptr[1] - a_ptr[1] * x_ptr[0];
		    temp_r += a_ptr[2] * x_ptr[2] + a_ptr[3] * x_ptr[3];
		    temp_i += a_ptr[2] * x_ptr[3] - a_ptr[3] * x_ptr[2];
#endif

		    a_ptr += 4;
		    x_ptr += 4;
	    }
#if ( GMS_ZGEMV_EXPLICITLY_VECTORIZE) == 1
#if defined __INTEL_COMPILER || defined __ICC
#pragma vectorlengthfor(double) reduction(+:temp_r,temp_i)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd reduction(+:temp_r,temp_i)
#endif
#endif
	    for(; i != n; ++i) {
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
		    temp_r += a_ptr[0] * x_ptr[0] - a_ptr[1] * x_ptr[1];
		    temp_i += a_ptr[0] * x_ptr[1] + a_ptr[1] * x_ptr[0];
#else
		    temp_r += a_ptr[0] * x_ptr[0] + a_ptr[1] * x_ptr[1];
		    temp_i += a_ptr[0] * x_ptr[1] - a_ptr[1] * x_ptr[0];
#endif

		    a_ptr += 2;
		    x_ptr += 2;
	    }
	 }
	 else {
#if ( GMS_ZGEMV_EXPLICITLY_VECTORIZE) == 1
#if defined __INTEL_COMPILER || defined __ICC
#pragma vectorlengthfor(double) reduction(+:temp_r,temp_i)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd reduction(+:temp_r,temp_i)
#endif
#endif
	        for(i = 0; i != n; ++i) {
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
		    temp_r += a_ptr[0] * x_ptr[0] - a_ptr[1] * x_ptr[1];
		    temp_i += a_ptr[0] * x_ptr[1] + a_ptr[1] * x_ptr[0];
#else
		    temp_r += a_ptr[0] * x_ptr[0] + a_ptr[1] * x_ptr[1];
		    temp_i += a_ptr[0] * x_ptr[1] - a_ptr[1] * x_ptr[0];
#endif

		    a_ptr += lda;
		    x_ptr += inc_x;
		}
	 }
#if !defined(XCONJ) 
		y_ptr[0] += alpha_r * temp_r - alpha_i * temp_i;
		y_ptr[1] += alpha_r * temp_i + alpha_i * temp_r;
#else
		y_ptr[0] += alpha_r * temp_r + alpha_i * temp_i;
		y_ptr[1] -= alpha_r * temp_i - alpha_i * temp_r;
#endif
		return(0);	 
      }

      if(m3 == 2) {

         a_ptr = a;
	 x_ptr = x;
	 double temp_r0 = 0.0;
	 double temp_i0 = 0.0;
	 double temp_r1 = 0.0;
	 double temp_i1 = 0.0;
	 double temp_r2 = 0.0;
	 double temp_i2 = 0.0;

	 if(lda==6 && inc_x==2) {
#if ( GMS_ZGEMV_EXPLICITLY_VECTORIZE) == 1
#if defined __INTEL_COMPILER || defined __ICC
#pragma vectorlengthfor(double) reduction(+:temp_r0,temp_i0,temp_r1,temp_i1,temp_r2,temp_i2)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd reduction(+:temp_r0,temp_i0,temp_r1,temp_i1,temp_r2,temp_i2)
#endif
#endif
	     for(i = 0; i != n; ++i) {
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
		 temp_r0 += a_ptr[0] * x_ptr[0] - a_ptr[1] * x_ptr[1];
		 temp_i0 += a_ptr[0] * x_ptr[1] + a_ptr[1] * x_ptr[0];
		 temp_r1 += a_ptr[2] * x_ptr[0] - a_ptr[3] * x_ptr[1];
		 temp_i1 += a_ptr[2] * x_ptr[1] + a_ptr[3] * x_ptr[0];
		 temp_r2 += a_ptr[4] * x_ptr[0] - a_ptr[5] * x_ptr[1];
		 temp_i2 += a_ptr[4] * x_ptr[1] + a_ptr[5] * x_ptr[0];
#else
		 temp_r0 += a_ptr[0] * x_ptr[0] + a_ptr[1] * x_ptr[1];
		 temp_i0 += a_ptr[0] * x_ptr[1] - a_ptr[1] * x_ptr[0];
		 temp_r1 += a_ptr[2] * x_ptr[0] + a_ptr[3] * x_ptr[1];
		 temp_i1 += a_ptr[2] * x_ptr[1] - a_ptr[3] * x_ptr[0];
		 temp_r2 += a_ptr[4] * x_ptr[0] + a_ptr[5] * x_ptr[1];
		 temp_i2 += a_ptr[4] * x_ptr[1] - a_ptr[5] * x_ptr[0];
#endif

		 a_ptr += 6;
		 x_ptr += 2;               
	     }
	 }
	 else {
#if ( GMS_ZGEMV_EXPLICITLY_VECTORIZE) == 1
#if defined __INTEL_COMPILER || defined __ICC
#pragma vectorlengthfor(double) reduction(+:temp_r0,temp_i0,temp_r1,temp_i1,temp_r2,temp_i2)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd reduction(+:temp_r0,temp_i0,temp_r1,temp_i1,temp_r2,temp_i2)
#endif
#endif
	      for(i = 0; i != n; ++i) {
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )

		   temp_r0 += a_ptr[0] * x_ptr[0] - a_ptr[1] * x_ptr[1];
		   temp_i0 += a_ptr[0] * x_ptr[1] + a_ptr[1] * x_ptr[0];
		   temp_r1 += a_ptr[2] * x_ptr[0] - a_ptr[3] * x_ptr[1];
		   temp_i1 += a_ptr[2] * x_ptr[1] + a_ptr[3] * x_ptr[0];
		   temp_r2 += a_ptr[4] * x_ptr[0] - a_ptr[5] * x_ptr[1];
		   temp_i2 += a_ptr[4] * x_ptr[1] + a_ptr[5] * x_ptr[0];
#else
		   temp_r0 += a_ptr[0] * x_ptr[0] + a_ptr[1] * x_ptr[1];
		   temp_i0 += a_ptr[0] * x_ptr[1] - a_ptr[1] * x_ptr[0];
		   temp_r1 += a_ptr[2] * x_ptr[0] + a_ptr[3] * x_ptr[1];
		   temp_i1 += a_ptr[2] * x_ptr[1] - a_ptr[3] * x_ptr[0];
		   temp_r2 += a_ptr[4] * x_ptr[0] + a_ptr[5] * x_ptr[1];
		   temp_i2 += a_ptr[4] * x_ptr[1] - a_ptr[5] * x_ptr[0];
#endif

				a_ptr += lda;
				x_ptr += inc_x;
	      }
	 }

#if !defined(XCONJ) 
		y_ptr[0] += alpha_r * temp_r0 - alpha_i * temp_i0;
		y_ptr[1] += alpha_r * temp_i0 + alpha_i * temp_r0;
		y_ptr    += inc_y;
		y_ptr[0] += alpha_r * temp_r1 - alpha_i * temp_i1;
		y_ptr[1] += alpha_r * temp_i1 + alpha_i * temp_r1;
		y_ptr    += inc_y;
		y_ptr[0] += alpha_r * temp_r2 - alpha_i * temp_i2;
		y_ptr[1] += alpha_r * temp_i2 + alpha_i * temp_r2;
#else
		y_ptr[0] += alpha_r * temp_r0 + alpha_i * temp_i0;
		y_ptr[1] -= alpha_r * temp_i0 - alpha_i * temp_r0;
		y_ptr    += inc_y;
		y_ptr[0] += alpha_r * temp_r1 + alpha_i * temp_i1;
		y_ptr[1] -= alpha_r * temp_i1 - alpha_i * temp_r1;
		y_ptr    += inc_y;
		y_ptr[0] += alpha_r * temp_r2 + alpha_i * temp_i2;
		y_ptr[1] -= alpha_r * temp_i2 - alpha_i * temp_r2;
#endif
		return(0);	 
      }
}









