



//=========================================================================
// Modified and optimized version of OpenBLAS zscal and kernels zscal_kernel_x
// Programmer: Bernard Gingold, contact: beniekg@gmail.com
// 13-02-2021 14:52  +00200
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

#include "GMS_zscal.h" 


void zscal_kernel_8(const int32_t n,
		    double * __restrict alpha,
		    double * __restrict x) {
      bool isAligned = false;
      isAligned32 = ((uintptr_t)alpha & 0x1F) != 0ULL && ((uintptr_t)x & 0x1F) != 0ULL;
      if(!isAligned32) {

            __asm__ __volatile__ (
                        "vzeroupper                                         \n\t"
                 	"vbroadcastsd		(%2), %%ymm0		    \n\t"  // da_r	
	                "vbroadcastsd          8(%2), %%ymm1		    \n\t"  // da_i 	

	                "addq	$128, %1				    \n\t"

	                "vmovups	-128(%1), %%ymm4		    \n\t"
	                "vmovups	 -96(%1), %%ymm5		    \n\t"
	                "vmovups	 -64(%1), %%ymm6		    \n\t"
	                "vmovups	 -32(%1), %%ymm7		    \n\t"

	                "vpermilpd	$0x05 , %%ymm4, %%ymm12		    \n\t"
	                "vpermilpd	$0x05 , %%ymm5, %%ymm13		    \n\t"
	                "vpermilpd	$0x05 , %%ymm6, %%ymm14 	    \n\t"
	                "vpermilpd	$0x05 , %%ymm7, %%ymm15		    \n\t"

	                "subq	        $8 , %0			            \n\t"		
	                "jz	2f					    \n\t"

	                ".p2align 4				            \n\t"
	                "1:				            	    \n\t"

	//"prefetcht0     128(%1)				    \n\t"
	// ".align 2				            \n\t"

	                "vmulpd		%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	                "vmovups	   0(%1), %%ymm4		    \n\t"
	                "vmulpd		%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	                "vmovups	  32(%1), %%ymm5		    \n\t"
	                "vmulpd		%%ymm0, %%ymm6 , %%ymm10	    \n\t" 
	                "vmovups	  64(%1), %%ymm6		    \n\t"
	                "vmulpd		%%ymm0, %%ymm7 , %%ymm11	    \n\t" 
	                "vmovups	  96(%1), %%ymm7		    \n\t"

	                "vmulpd		%%ymm1, %%ymm12, %%ymm12	    \n\t" // da_i*x1 , da_i *x0
	                "vaddsubpd	%%ymm12 , %%ymm8 , %%ymm8	    \n\t"
	                "vmulpd		%%ymm1, %%ymm13, %%ymm13	    \n\t" 
	                "vaddsubpd	%%ymm13 , %%ymm9 , %%ymm9	    \n\t"
	                "vmulpd		%%ymm1, %%ymm14, %%ymm14	    \n\t" 
	                "vaddsubpd	%%ymm14 , %%ymm10, %%ymm10	    \n\t"
	                "vmulpd		%%ymm1, %%ymm15, %%ymm15	    \n\t" 
	                "vaddsubpd	%%ymm15 , %%ymm11, %%ymm11	    \n\t"

	                "vmovups	%%ymm8 , -128(%1)		    \n\t"
	                "vpermilpd	$0x05 , %%ymm4, %%ymm12		    \n\t"
	                "vmovups	%%ymm9 ,  -96(%1)		    \n\t"
	                "vpermilpd	$0x05 , %%ymm5, %%ymm13		    \n\t"
	                "vmovups	%%ymm10,  -64(%1)		    \n\t"
	                "vpermilpd	$0x05 , %%ymm6, %%ymm14 	    \n\t"
	                "vmovups	%%ymm11,  -32(%1)		    \n\t"
	                "vpermilpd	$0x05 , %%ymm7, %%ymm15		    \n\t"

	                "addq		$128 ,%1  	 	            \n\t"
	                "subq	        $8 , %0			            \n\t"		
	                "jnz		1b		             	    \n\t"

	                "2:				            	    \n\t"


	                "vmulpd		%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	                "vmulpd		%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	                "vmulpd		%%ymm0, %%ymm6 , %%ymm10	    \n\t" 
	                "vmulpd		%%ymm0, %%ymm7 , %%ymm11	    \n\t" 

	                "vmulpd		%%ymm1, %%ymm12, %%ymm12	    \n\t" // da_i*x1 , da_i *x0
	                "vaddsubpd	%%ymm12 , %%ymm8 , %%ymm8	    \n\t"
	                "vmulpd		%%ymm1, %%ymm13, %%ymm13	    \n\t" 
	                "vaddsubpd	%%ymm13 , %%ymm9 , %%ymm9	    \n\t"
	                "vmulpd		%%ymm1, %%ymm14, %%ymm14	    \n\t" 
	                "vaddsubpd	%%ymm14 , %%ymm10, %%ymm10	    \n\t"
	                "vmulpd		%%ymm1, %%ymm15, %%ymm15	    \n\t" 
	                "vaddsubpd	%%ymm15 , %%ymm11, %%ymm11	    \n\t"

	                "vmovups	%%ymm8 , -128(%1)		    \n\t"
	                "vmovups	%%ymm9 ,  -96(%1)		    \n\t"
	                "vmovups	%%ymm10,  -64(%1)		    \n\t"
	                "vmovups	%%ymm11,  -32(%1)		    \n\t"

	                "vzeroupper					    \n\t"

	               :
	               "+r" (n),  	// 0
                       "+r" (x)      // 1
	               :
                       "r" (alpha)   // 2
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
                        "vzeroupper                                         \n\t"
                 	"vbroadcastsd		(%2), %%ymm0		    \n\t"  // da_r	
	                "vbroadcastsd          8(%2), %%ymm1		    \n\t"  // da_i 	

	                "addq	$128, %1				    \n\t"

	                "vmovaps	-128(%1), %%ymm4		    \n\t"
	                "vmovaps	 -96(%1), %%ymm5		    \n\t"
	                "vmovaps	 -64(%1), %%ymm6		    \n\t"
	                "vmovaps	 -32(%1), %%ymm7		    \n\t"

	                "vpermilpd	$0x05 , %%ymm4, %%ymm12		    \n\t"
	                "vpermilpd	$0x05 , %%ymm5, %%ymm13		    \n\t"
	                "vpermilpd	$0x05 , %%ymm6, %%ymm14 	    \n\t"
	                "vpermilpd	$0x05 , %%ymm7, %%ymm15		    \n\t"

	                "subq	        $8 , %0			            \n\t"		
	                "jz	2f					    \n\t"

	                ".p2align 4				            \n\t"
	                "1:				            	    \n\t"

	//"prefetcht0     128(%1)				    \n\t"
	// ".align 2				            \n\t"

	                "vmulpd		%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	                "vmovaps	   0(%1), %%ymm4		    \n\t"
	                "vmulpd		%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	                "vmovaps	  32(%1), %%ymm5		    \n\t"
	                "vmulpd		%%ymm0, %%ymm6 , %%ymm10	    \n\t" 
	                "vmovaps	  64(%1), %%ymm6		    \n\t"
	                "vmulpd		%%ymm0, %%ymm7 , %%ymm11	    \n\t" 
	                "vmovaps	  96(%1), %%ymm7		    \n\t"

	                "vmulpd		%%ymm1, %%ymm12, %%ymm12	    \n\t" // da_i*x1 , da_i *x0
	                "vaddsubpd	%%ymm12 , %%ymm8 , %%ymm8	    \n\t"
	                "vmulpd		%%ymm1, %%ymm13, %%ymm13	    \n\t" 
	                "vaddsubpd	%%ymm13 , %%ymm9 , %%ymm9	    \n\t"
	                "vmulpd		%%ymm1, %%ymm14, %%ymm14	    \n\t" 
	                "vaddsubpd	%%ymm14 , %%ymm10, %%ymm10	    \n\t"
	                "vmulpd		%%ymm1, %%ymm15, %%ymm15	    \n\t" 
	                "vaddsubpd	%%ymm15 , %%ymm11, %%ymm11	    \n\t"

	                "vmovaps	%%ymm8 , -128(%1)		    \n\t"
	                "vpermilpd	$0x05 , %%ymm4, %%ymm12		    \n\t"
	                "vmovaps	%%ymm9 ,  -96(%1)		    \n\t"
	                "vpermilpd	$0x05 , %%ymm5, %%ymm13		    \n\t"
	                "vmovaps	%%ymm10,  -64(%1)		    \n\t"
	                "vpermilpd	$0x05 , %%ymm6, %%ymm14 	    \n\t"
	                "vmovaps	%%ymm11,  -32(%1)		    \n\t"
	                "vpermilpd	$0x05 , %%ymm7, %%ymm15		    \n\t"

	                "addq		$128 ,%1  	 	            \n\t"
	                "subq	        $8 , %0			            \n\t"		
	                "jnz		1b		             	    \n\t"

	                "2:				            	    \n\t"


	                "vmulpd		%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	                "vmulpd		%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	                "vmulpd		%%ymm0, %%ymm6 , %%ymm10	    \n\t" 
	                "vmulpd		%%ymm0, %%ymm7 , %%ymm11	    \n\t" 

	                "vmulpd		%%ymm1, %%ymm12, %%ymm12	    \n\t" // da_i*x1 , da_i *x0
	                "vaddsubpd	%%ymm12 , %%ymm8 , %%ymm8	    \n\t"
	                "vmulpd		%%ymm1, %%ymm13, %%ymm13	    \n\t" 
	                "vaddsubpd	%%ymm13 , %%ymm9 , %%ymm9	    \n\t"
	                "vmulpd		%%ymm1, %%ymm14, %%ymm14	    \n\t" 
	                "vaddsubpd	%%ymm14 , %%ymm10, %%ymm10	    \n\t"
	                "vmulpd		%%ymm1, %%ymm15, %%ymm15	    \n\t" 
	                "vaddsubpd	%%ymm15 , %%ymm11, %%ymm11	    \n\t"

	                "vmovaps	%%ymm8 , -128(%1)		    \n\t"
	                "vmovaps	%%ymm9 ,  -96(%1)		    \n\t"
	                "vmovaps	%%ymm10,  -64(%1)		    \n\t"
	                "vmovaps	%%ymm11,  -32(%1)		    \n\t"

	                "vzeroupper					    \n\t"

	               :
	               "+r" (n),  	// 0
                       "+r" (x)      // 1
	               :
                       "r" (alpha)   // 2
	               : "cc",
	               "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	               "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	               "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	               "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	               "memory"
	    );
      }
}



void zscal_kernel_8_zero_r(const int32_t n,
                           double * __restrict alpha,
			   double * __restrict x) {
         bool isAligned32 = false;
         isAligned32 = ((uintptr_t)alpha & 0x1F) != 0ULL && ((uintptr_t)x & 0x1F) != 0ULL;
	 if(!isAligned32) {
	         __asm__ __volatile__ (

		           "vzeroupper                                              \n\t "
                           "vxorpd	           %%ymm0, %%ymm0, %%ymm0	    \n\t"
	                   "vbroadcastsd          8(%2), %%ymm1		    \n\t"  // da_i

	                   "addq	$128, %1				    \n\t"

	                   "vmovups	-128(%1), %%ymm4		    \n\t"
	                   "vmovups	 -96(%1), %%ymm5		    \n\t"
	                   "vmovups	 -64(%1), %%ymm6		    \n\t"
	                   "vmovups	 -32(%1), %%ymm7		    \n\t"

	                   "vpermilpd	$0x05 , %%ymm4, %%ymm12		    \n\t"
	                   "vpermilpd	$0x05 , %%ymm5, %%ymm13		    \n\t"
	                   "vpermilpd	$0x05 , %%ymm6, %%ymm14 	    \n\t"
	                   "vpermilpd	$0x05 , %%ymm7, %%ymm15		    \n\t"

	                   "subq	        $8 , %0			            \n\t"		
	                   "jz	2f					    \n\t"

	                   ".p2align 4				            \n\t"
	                   "1:				            	    \n\t"

	//"prefetcht0     128(%1)				    \n\t" // Software prefetching not needed here.
	// ".align 2				            \n\t"
#if (GMS_MAN_PREFETCH) == 1			   
                           "prefetcht0     256(%1)				    \n\t"
#endif
	                   "vmovups	   0(%1), %%ymm4		    \n\t"
	                   "vmovups	  32(%1), %%ymm5		    \n\t"
	                   "vmovups	  64(%1), %%ymm6		    \n\t"
	                   "vmovups	  96(%1), %%ymm7		    \n\t"

	                   "vmulpd		%%ymm1, %%ymm12, %%ymm12	    \n\t" // da_i*x1 , da_i *x0
	                   "vaddsubpd	%%ymm12 , %%ymm0 , %%ymm8	    \n\t"
	                   "vmulpd		%%ymm1, %%ymm13, %%ymm13	    \n\t" 
	                   "vaddsubpd	%%ymm13 , %%ymm0 , %%ymm9	    \n\t"
	                   "vmulpd		%%ymm1, %%ymm14, %%ymm14	    \n\t" 
	                   "vaddsubpd	%%ymm14 , %%ymm0 , %%ymm10	    \n\t"
	                   "vmulpd		%%ymm1, %%ymm15, %%ymm15	    \n\t" 
	                   "vaddsubpd	%%ymm15 , %%ymm0 , %%ymm11	    \n\t"

	                   "vmovups	%%ymm8 , -128(%1)		    \n\t"
	                   "vpermilpd	$0x05 , %%ymm4, %%ymm12		    \n\t"
	                   "vmovups	%%ymm9 ,  -96(%1)		    \n\t"
	                   "vpermilpd	$0x05 , %%ymm5, %%ymm13		    \n\t"
	                   "vmovups	%%ymm10,  -64(%1)		    \n\t"
	                   "vpermilpd	$0x05 , %%ymm6, %%ymm14 	    \n\t"
	                   "vmovups	%%ymm11,  -32(%1)		    \n\t"
	                   "vpermilpd	$0x05 , %%ymm7, %%ymm15		    \n\t"

	                   "addq		$128 ,%1  	 	            \n\t"
	                   "subq	        $8 , %0			            \n\t"		
	                   "jnz		1b		             	    \n\t"

	                   "2:				            	    \n\t"

	                   "vmulpd		%%ymm1, %%ymm12, %%ymm12	    \n\t" // da_i*x1 , da_i *x0
	                   "vaddsubpd	%%ymm12 , %%ymm0 , %%ymm8	    \n\t"
	                   "vmulpd		%%ymm1, %%ymm13, %%ymm13	    \n\t" 
	                   "vaddsubpd	%%ymm13 , %%ymm0 , %%ymm9	    \n\t"
	                   "vmulpd		%%ymm1, %%ymm14, %%ymm14	    \n\t" 
	                   "vaddsubpd	%%ymm14 , %%ymm0 , %%ymm10	    \n\t"
	                   "vmulpd		%%ymm1, %%ymm15, %%ymm15	    \n\t" 
	                   "vaddsubpd	%%ymm15 , %%ymm0 , %%ymm11	    \n\t"

	                   "vmovups	%%ymm8 , -128(%1)		    \n\t"
	                   "vmovups	%%ymm9 ,  -96(%1)		    \n\t"
	                   "vmovups	%%ymm10,  -64(%1)		    \n\t"
	                   "vmovups	%%ymm11,  -32(%1)		    \n\t"

	                   "vzeroupper					    \n\t"

	                   :
	                   "+r" (n),  	// 0
                           "+r" (x)      // 1
	                   :
                           "r" (alpha)   // 2
	                   : "cc",
	                   "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	                   "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	                   "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	                   "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	                   "memory"
	            );

       }else {

                   __asm__ __volatile__ (

		           "vzeroupper                                              \n\t "
                           "vxorpd	           %%ymm0, %%ymm0, %%ymm0	    \n\t"
	                   "vbroadcastsd          8(%2), %%ymm1		    \n\t"  // da_i 	

	                   "addq	$128, %1				    \n\t"

	                   "vmovaps	-128(%1), %%ymm4		    \n\t"
	                   "vmovaps	 -96(%1), %%ymm5		    \n\t"
	                   "vmovaps	 -64(%1), %%ymm6		    \n\t"
	                   "vmovaps	 -32(%1), %%ymm7		    \n\t"

	                   "vpermilpd	$0x05 , %%ymm4, %%ymm12		    \n\t"
	                   "vpermilpd	$0x05 , %%ymm5, %%ymm13		    \n\t"
	                   "vpermilpd	$0x05 , %%ymm6, %%ymm14 	    \n\t"
	                   "vpermilpd	$0x05 , %%ymm7, %%ymm15		    \n\t"

	                   "subq	        $8 , %0			            \n\t"		
	                   "jz	2f					    \n\t"

	                   ".p2align 4				            \n\t"
	                   "1:				            	    \n\t"

	//"prefetcht0     128(%1)				    \n\t"
	// ".align 2				            \n\t"
#if (GMS_MAN_PREFETCH) == 1			   
                           "prefetcht0     256(%1)				    \n\t"
#endif
	                   "vmovaps	   0(%1), %%ymm4		    \n\t"
	                   "vmovaps	  32(%1), %%ymm5		    \n\t"
	                   "vmovaps	  64(%1), %%ymm6		    \n\t"
	                   "vmovaps	  96(%1), %%ymm7		    \n\t"

	                   "vmulpd		%%ymm1, %%ymm12, %%ymm12	    \n\t" // da_i*x1 , da_i *x0
	                   "vaddsubpd	%%ymm12 , %%ymm0 , %%ymm8	    \n\t"
	                   "vmulpd		%%ymm1, %%ymm13, %%ymm13	    \n\t" 
	                   "vaddsubpd	%%ymm13 , %%ymm0 , %%ymm9	    \n\t"
	                   "vmulpd		%%ymm1, %%ymm14, %%ymm14	    \n\t" 
	                   "vaddsubpd	%%ymm14 , %%ymm0 , %%ymm10	    \n\t"
	                   "vmulpd		%%ymm1, %%ymm15, %%ymm15	    \n\t" 
	                   "vaddsubpd	%%ymm15 , %%ymm0 , %%ymm11	    \n\t"

	                   "vmovaps	%%ymm8 , -128(%1)		    \n\t"
	                   "vpermilpd	$0x05 , %%ymm4, %%ymm12		    \n\t"
	                   "vmovaps	%%ymm9 ,  -96(%1)		    \n\t"
	                   "vpermilpd	$0x05 , %%ymm5, %%ymm13		    \n\t"
	                   "vmovaps	%%ymm10,  -64(%1)		    \n\t"
	                   "vpermilpd	$0x05 , %%ymm6, %%ymm14 	    \n\t"
	                   "vmovaps	%%ymm11,  -32(%1)		    \n\t"
	                   "vpermilpd	$0x05 , %%ymm7, %%ymm15		    \n\t"

	                   "addq		$128 ,%1  	 	            \n\t"
	                   "subq	        $8 , %0			            \n\t"		
	                   "jnz		1b		             	    \n\t"

	                   "2:				            	    \n\t"

	                   "vmulpd		%%ymm1, %%ymm12, %%ymm12	    \n\t" // da_i*x1 , da_i *x0
	                   "vaddsubpd	%%ymm12 , %%ymm0 , %%ymm8	    \n\t"
	                   "vmulpd		%%ymm1, %%ymm13, %%ymm13	    \n\t" 
	                   "vaddsubpd	%%ymm13 , %%ymm0 , %%ymm9	    \n\t"
	                   "vmulpd		%%ymm1, %%ymm14, %%ymm14	    \n\t" 
	                   "vaddsubpd	%%ymm14 , %%ymm0 , %%ymm10	    \n\t"
	                   "vmulpd		%%ymm1, %%ymm15, %%ymm15	    \n\t" 
	                   "vaddsubpd	%%ymm15 , %%ymm0 , %%ymm11	    \n\t"

	                   "vmovaps	%%ymm8 , -128(%1)		    \n\t"
	                   "vmovaps	%%ymm9 ,  -96(%1)		    \n\t"
	                   "vmovaps	%%ymm10,  -64(%1)		    \n\t"
	                   "vmovaps	%%ymm11,  -32(%1)		    \n\t"

	                   "vzeroupper					    \n\t"

	                   :
	                   "+r" (n),  	// 0
                           "+r" (x)      // 1
	                   :
                           "r" (alpha)   // 2
	                   : "cc",
	                   "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	                   "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	                   "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	                   "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	                   "memory"
	            );
          
    
       }
}



void zscal_kernel_8_zero_i(const int32_t n,
			   double * __restrict alpha,
			   double * __restrict x) {
      bool isAligned32 = false;
      isAligned32 = ((uintptr_t)alpha & 0x1F) != 0ULL && ((uintptr_t)x & 0x1F) != 0ULL;
      if(!isAligned32) {
          __asm__ __volatile__ (

	           	"vbroadcastsd		(%2), %%ymm0		    \n\t"  // da_r	

	                "addq	$128, %1				    \n\t"

	                "vmovups	-128(%1), %%ymm4		    \n\t"
	                "vmovups	 -96(%1), %%ymm5		    \n\t"
	                "vmovups	 -64(%1), %%ymm6		    \n\t"
	                "vmovups	 -32(%1), %%ymm7		    \n\t"


	                "subq	        $8 , %0			            \n\t"		
	                "jz	2f					    \n\t"

	                ".p2align 4				            \n\t"
	                "1:				            	    \n\t"

	//"prefetcht0     128(%1)				    \n\t"
	// ".align 2				            \n\t"
#if (GMS_MAN_PREFETCH) == 1			   
                           "prefetcht0     256(%1)				    \n\t"
#endif
	                "vmulpd		%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	                "vmovups	   0(%1), %%ymm4		    \n\t"
	                "vmulpd		%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	                "vmovups	  32(%1), %%ymm5		    \n\t"
	                "vmulpd		%%ymm0, %%ymm6 , %%ymm10	    \n\t" 
	                "vmovups	  64(%1), %%ymm6		    \n\t"
	                "vmulpd		%%ymm0, %%ymm7 , %%ymm11	    \n\t" 
	                "vmovups	  96(%1), %%ymm7		    \n\t"

	                "vmovups	%%ymm8 , -128(%1)		    \n\t"
	                "vmovups	%%ymm9 ,  -96(%1)		    \n\t"
	                "vmovups	%%ymm10,  -64(%1)		    \n\t"
	                "vmovups	%%ymm11,  -32(%1)		    \n\t"

	                "addq		$128 ,%1  	 	            \n\t"
	                "subq	        $8 , %0			            \n\t"		
	                "jnz		1b		             	    \n\t"

	                "2:				            	    \n\t"


	                "vmulpd		%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	                "vmulpd		%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	                "vmulpd		%%ymm0, %%ymm6 , %%ymm10	    \n\t" 
	                "vmulpd		%%ymm0, %%ymm7 , %%ymm11	    \n\t" 

	                "vmovups	%%ymm8 , -128(%1)		    \n\t"
	                "vmovups	%%ymm9 ,  -96(%1)		    \n\t"
	                "vmovups	%%ymm10,  -64(%1)		    \n\t"
	                "vmovups	%%ymm11,  -32(%1)		    \n\t"

	                "vzeroupper					    \n\t"

	                :
	                "+r" (n),  	// 0
                        "+r" (x)      // 1
	                :
                        "r" (alpha)   // 2
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

	           	"vbroadcastsd		(%2), %%ymm0		    \n\t"  // da_r	

	                "addq	$128, %1				    \n\t"

	                "vmovaps	-128(%1), %%ymm4		    \n\t"
	                "vmovaps	 -96(%1), %%ymm5		    \n\t"
	                "vmovaps	 -64(%1), %%ymm6		    \n\t"
	                "vmovaps	 -32(%1), %%ymm7		    \n\t"


	                "subq	        $8 , %0			            \n\t"		
	                "jz	2f					    \n\t"

	                ".p2align 4				            \n\t"
	                "1:				            	    \n\t"

	//"prefetcht0     128(%1)				    \n\t"
	// ".align 2				            \n\t"
#if (GMS_MAN_PREFETCH) == 1			   
                           "prefetcht0     256(%1)				    \n\t"
#endif
	                "vmulpd		%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	                "vmovaps	   0(%1), %%ymm4		    \n\t"
	                "vmulpd		%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	                "vmovaps	  32(%1), %%ymm5		    \n\t"
	                "vmulpd		%%ymm0, %%ymm6 , %%ymm10	    \n\t" 
	                "vmovaps	  64(%1), %%ymm6		    \n\t"
	                "vmulpd		%%ymm0, %%ymm7 , %%ymm11	    \n\t" 
	                "vmovaps	  96(%1), %%ymm7		    \n\t"

	                "vmovaps	%%ymm8 , -128(%1)		    \n\t"
	                "vmovaps	%%ymm9 ,  -96(%1)		    \n\t"
	                "vmovaps	%%ymm10,  -64(%1)		    \n\t"
	                "vmovaps	%%ymm11,  -32(%1)		    \n\t"

	                "addq		$128 ,%1  	 	            \n\t"
	                "subq	        $8 , %0			            \n\t"		
	                "jnz		1b		             	    \n\t"

	                "2:				            	    \n\t"


	                "vmulpd		%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	                "vmulpd		%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	                "vmulpd		%%ymm0, %%ymm6 , %%ymm10	    \n\t" 
	                "vmulpd		%%ymm0, %%ymm7 , %%ymm11	    \n\t" 

	                "vmovaps	%%ymm8 , -128(%1)		    \n\t"
	                "vmovaps	%%ymm9 ,  -96(%1)		    \n\t"
	                "vmovaps	%%ymm10,  -64(%1)		    \n\t"
	                "vmovaps	%%ymm11,  -32(%1)		    \n\t"

	                "vzeroupper					    \n\t"

	                :
	                "+r" (n),  	// 0
                        "+r" (x)      // 1
	                :
                        "r" (alpha)   // 2
	                : "cc",
	                "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	                "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	                "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	                "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	                "memory"
	              );
      }
}
			   


void zscal_kernel_8_zero(const int32_t n,
                         double * __restrict alpha,
			 double * __restrict x) {
       bool isAligned32 = false;
       isAligned32 =  ((uintptr_t)alpha & 0x1F) != 0ULL && ((uintptr_t)x & 0x1F) != 0ULL;
       if(!isAligned32) {

             __asm__ __volatile__ (

                 	"vxorpd	           %%ymm0, %%ymm0, %%ymm0	    \n\t"

	                "addq	$128, %1				    \n\t"

	                ".p2align 4				            \n\t"
	                "1:				            	    \n\t"

	           //"prefetcht0     128(%1)				    \n\t"
	// ".align 2				            \n\t"
#if (GMS_MAN_PREFETCH) == 1			   
                           "prefetcht0     256(%1)				    \n\t"
#endif
	                "vmovups	%%ymm0 , -128(%1)		    \n\t"
	                "vmovups	%%ymm0 ,  -96(%1)		    \n\t"
	                "vmovups	%%ymm0 ,  -64(%1)		    \n\t"
	                "vmovups	%%ymm0 ,  -32(%1)		    \n\t"

	                "addq		$128 ,%1  	 	            \n\t"
	                "subq	        $8 , %0			            \n\t"		
	                "jnz		1b		             	    \n\t"

	                "vzeroupper					    \n\t"

	                :
	                "+r" (n),  	// 0
                        "+r" (x)      // 1
	                :
                        "r" (alpha)   // 2
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

                 	"vxorpd	           %%ymm0, %%ymm0, %%ymm0	    \n\t"

	                "addq	$128, %1				    \n\t"

	                ".p2align 4				            \n\t"
	                "1:				            	    \n\t"

	           //"prefetcht0     128(%1)				    \n\t"
	// ".align 2				            \n\t"
#if (GMS_MAN_PREFETCH) == 1			   
                           "prefetcht0     256(%1)				    \n\t"
#endif
	                "vmovaps	%%ymm0 , -128(%1)		    \n\t"
	                "vmovaps	%%ymm0 ,  -96(%1)		    \n\t"
	                "vmovaps	%%ymm0 ,  -64(%1)		    \n\t"
	                "vmovaps	%%ymm0 ,  -32(%1)		    \n\t"

	                "addq		$128 ,%1  	 	            \n\t"
	                "subq	        $8 , %0			            \n\t"		
	                "jnz		1b		             	    \n\t"

	                "vzeroupper					    \n\t"

	                :
	                "+r" (n),  	// 0
                        "+r" (x)      // 1
	                :
                        "r" (alpha)   // 2
	                : "cc",
	                "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	                "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	                "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	                "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	                "memory"
	         );
       }
}

_
void zscal_kernel_inc_8(const int32_t n,
                        double * __restrict alpha,
			double * __restrict x,
			const int32_t inc_x){


	int32_t i;
	const int32_t inc_x2 = 2 * inc_x;
	const int32_t inc_x3 = inc_x2 + inc_x;
	double t0,t1,t2,t3;
	double da_r = alpha[0];
	double da_i = alpha[1];
        t0 = 0.0;
	t1 = 0.0;
	t2 = 0.0;
	t3 = 0.0;
#if defined __INTEL_COMPILER || defined __ICC
#pragma code_align(16)
#pragma simd vectorlengthfor(double)
#elif defined __GNUC__ || !defined __INTEL_COMPILER
#pragma omp simd
#endif
	for (i=0; i != n; i +=4) {
	
		t0 = da_r * x[0]         - da_i *x[1];	
		t1 = da_r * x[inc_x]     - da_i *x[inc_x  + 1];	
		t2 = da_r * x[inc_x2]    - da_i *x[inc_x2 + 1];	
		t3 = da_r * x[inc_x3]    - da_i *x[inc_x3 + 1];	

		x[1]               = da_i * x[0]       + da_r * x[1];
		x[inc_x  +1]       = da_i * x[inc_x]   + da_r * x[inc_x  +1];
		x[inc_x2 +1]       = da_i * x[inc_x2]  + da_r * x[inc_x2 +1];
		x[inc_x3 +1]       = da_i * x[inc_x3]  + da_r * x[inc_x3 +1];

		x[0]        = t0;
		x[inc_x]    = t1;
		x[inc_x2]   = t2;
		x[inc_x3]   = t3;

		x+=4*inc_x;

	}
}



void zscal(const int32_t n,
           const double da_r,
	   const double da_i,
           double * __restrict x,
	   const int32_t inc_x,
	   double * __restrict y,
	   const int32_t inc_y) {

      double alpha[2];
      double temp0;
      double temp1;
      int32_t i,j;
      i = 0;
      j = 0;
      if(inc_x != 1) {
         inc_x <<= 1;
	 if(da_r == 0.0) {
	    const int32_t n1 = n & -2;
	    if(da_i == 0.0) {
	    
	       while(j < n1) {
                  x[i] = 0.0;
		  x[i+1] = 0.0;
		  x[i+inc_x] = 0.0;
		  x[i+1+inc_x] = 0.0;
		  i += incx+inc_x;
		  j += 2;
	       }

	       while(j < n) {
                  x[i] = 0.0;
		  x[i+1] = 0.0;
		  i += inc_x;
		  j++;
	       }
      }
      else {
               	while(j < n1){
		    temp0  = -da_i * x[i+1];
		    x[i+1] =  da_i * x[i];
		    x[i]   =  temp0;
		    temp1  = -da_i * x[i+1+inc_x];
		    x[i+1+inc_x] =  da_i * x[i+inc_x];
		    x[i+inc_x]   =  temp1;
		    i += 2*inc_x ;
		    j+=2;
	       }

	       	while(j < n){
		    temp0 = -da_i * x[i+1];
		    x[i+1] =  da_i * x[i];
		    x[i]   =  temp0;
		    i += inc_x ;
		    j++;
	       }
	       
      }
  }
     else {

        if(da_i == 0.0) {
           const int32_t n1 = n & -2;

	   while(j < n1) {
               	temp0 =  da_r * x[i];
		x[i+1] =  da_r * x[i+1];
		x[i]  =  temp0;
		temp1 =  da_r * x[i+inc_x];
		x[i+1+inc_x] =  da_r * x[i+1+inc_x];
		x[i+inc_x]   =  temp1;
		i += 2*inc_x ;
		j+=2;
	   }

	   while(j < n) {
               	temp0 =  da_r * x[i];
		x[i+1] =  da_r * x[i+1];
		x[i]   =  temp0;
		i += inc_x ;
		j++;
	   }
	}
	else {
             const int32_t n1 = n & -8;
	     if(n1 > 0) {
                  alpha[0] = da_r;
		  alpha[1] = da_i;
		  zscal_kernel_inc_8(n1, alpha, x, inc_x);
		  j = n1;
		  i = n1 * inc_x;
	     }
	}
    }
}

	const int32_t n1 = n & -8;
	if (n1 > 0){
	     alpha[0] = da_r;
	     alpha[1] = da_i;
	     if(da_r == 0.0)
		    if(da_i == 0)
			   zscal_kernel_8_zero(n1,alpha,x);
		        else
			   zscal_kernel_8_zero_r(n1,alpha,x);
	     else
			if(da_i == 0)
			    zscal_kernel_8_zero_i(n1 , alpha , x);
			else
		          zscal_kernel_8(n1 , alpha , x);

		i = n1 << 1;
		j = n1;
	}


	if ( da_r == 0.0 ){
	      if ( da_i == 0.0 ){
		   while(j < n){
			x[i]=0.0;
			x[i+1]=0.0;
			i += 2 ;
			j++;
	          	}

		}
		else {
		  while(j < n){
		      temp0 = -da_i * x[i+1];
		      x[i+1] =  da_i * x[i];
		      x[i]   =  temp0;
		      i += 2 ;
		      j++;
                    }
	       }

	}
	else
	{

		if ( da_i == 0.0 )
		{

			while(j < n)
			{
			
					temp0        =  da_r * x[i];
					x[i+1]       =  da_r * x[i+1];
					x[i]         =  temp0;
					i += 2 ;
					j++;

			}

		}
		else
		{

			while(j < n)
			{

				temp0        =  da_r * x[i]   - da_i * x[i+1];
				x[i+1]       =  da_r * x[i+1] + da_i * x[i];
				x[i]         =  temp0;
				i += 2 ;
				j++;

			}

		}		

	}

	return(0);
}

}






