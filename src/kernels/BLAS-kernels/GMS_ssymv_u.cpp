




//=========================================================================
// Modified and optimized version of OpenBLAS ssymv_u and kernels kernel_ssymv_u
// Programmer: Bernard Gingold, contact: beniekg@gmail.com
// 20-02-2021 14:58  +00200
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

#if defined __GNUC__ && !defined __INTEL_COMPILER
#include <omp.h>
#endif
#include "GMS_ssymv_u.h" 


void ssymv_kernel_4x4(const int32_t n,
                      float * __restrict a0,
		      float * __restrict a1,
		      float * __restrict a2,
		      float * __restrict a3,
		      float * __restrict x,
		      float * __restrict y,
		      float * __restrict temp1,
		      float * __restrict temp2) {

        int32_t i = 0;
	__asm__ __volatile__ (
              	"vzeroupper				     \n\t"
	"vxorps		%%ymm0 , %%ymm0 , %%ymm0     \n\t"	// temp2[0]
	"vxorps		%%ymm1 , %%ymm1 , %%ymm1     \n\t"	// temp2[1]
	"vxorps		%%ymm2 , %%ymm2 , %%ymm2     \n\t"	// temp2[2]
	"vxorps		%%ymm3 , %%ymm3 , %%ymm3     \n\t"	// temp2[3]
	"vbroadcastss   (%8),    %%ymm4	             \n\t"	// temp1[0]
	"vbroadcastss  4(%8),    %%ymm5	             \n\t"	// temp1[1]
	"vbroadcastss  8(%8),    %%ymm6	             \n\t"	// temp1[1]
	"vbroadcastss 12(%8),    %%ymm7	             \n\t"	// temp1[1]
	"xorq           %0,%0                        \n\t"

	".p2align 4				     \n\t"
	"1:				     \n\t"

	"vmovups	(%3,%0,4), %%ymm9	           \n\t"  // 2 * y
	"vmovups	(%2,%0,4), %%ymm8	           \n\t"  // 2 * x

	"vmovups	(%4,%0,4), %%ymm12	           \n\t"  // 2 * a
	"vmovups	(%5,%0,4), %%ymm13	           \n\t"  // 2 * a
	"vmovups	(%6,%0,4), %%ymm14	           \n\t"  // 2 * a
	"vmovups	(%7,%0,4), %%ymm15	           \n\t"  // 2 * a

	"vfmadd231ps	%%ymm4, %%ymm12 , %%ymm9  \n\t"  // y     += temp1 * a
	"vfmadd231ps	%%ymm8, %%ymm12 , %%ymm0  \n\t"  // temp2 += x * a

	"vfmadd231ps	%%ymm5, %%ymm13 , %%ymm9  \n\t"  // y     += temp1 * a
	"vfmadd231ps	%%ymm8, %%ymm13 , %%ymm1  \n\t"  // temp2 += x * a

	"vfmadd231ps	%%ymm6, %%ymm14 , %%ymm9  \n\t"  // y     += temp1 * a
	"vfmadd231ps	%%ymm8, %%ymm14 , %%ymm2  \n\t"  // temp2 += x * a

	"vfmadd231ps	%%ymm7, %%ymm15 , %%ymm9  \n\t"  // y     += temp1 * a
	"vfmadd231ps	%%ymm8, %%ymm15 , %%ymm3  \n\t"  // temp2 += x * a

	"vmovups	%%ymm9 ,  (%3,%0,4)	      \n\t"
	"addq		$8 , %0	  	 	      \n\t"
	"subq		$8 , %1	  	 	      \n\t"

	"jnz		1b		      \n\t"

	"vmovss		  (%9), %%xmm4		      \n\t"
	"vmovss		 4(%9), %%xmm5		      \n\t"
	"vmovss		 8(%9), %%xmm6		      \n\t"
	"vmovss		12(%9), %%xmm7		      \n\t"

	"vextractf128 $0x01, %%ymm0 , %%xmm12	      \n\t"
	"vextractf128 $0x01, %%ymm1 , %%xmm13	      \n\t"
	"vextractf128 $0x01, %%ymm2 , %%xmm14	      \n\t"
	"vextractf128 $0x01, %%ymm3 , %%xmm15	      \n\t"

	"vaddps	        %%xmm0, %%xmm12, %%xmm0	      \n\t"
	"vaddps	        %%xmm1, %%xmm13, %%xmm1	      \n\t"
	"vaddps	        %%xmm2, %%xmm14, %%xmm2	      \n\t"
	"vaddps	        %%xmm3, %%xmm15, %%xmm3	      \n\t"

	"vhaddps        %%xmm0, %%xmm0, %%xmm0  \n\t"
	"vhaddps        %%xmm1, %%xmm1, %%xmm1  \n\t"
	"vhaddps        %%xmm2, %%xmm2, %%xmm2  \n\t"
	"vhaddps        %%xmm3, %%xmm3, %%xmm3  \n\t"

	"vhaddps        %%xmm0, %%xmm0, %%xmm0  \n\t"
	"vhaddps        %%xmm1, %%xmm1, %%xmm1  \n\t"
	"vhaddps        %%xmm2, %%xmm2, %%xmm2  \n\t"
	"vhaddps        %%xmm3, %%xmm3, %%xmm3  \n\t"

	"vaddsd		%%xmm4, %%xmm0, %%xmm0  \n\t"
	"vaddsd		%%xmm5, %%xmm1, %%xmm1  \n\t"
	"vaddsd		%%xmm6, %%xmm2, %%xmm2  \n\t"
	"vaddsd		%%xmm7, %%xmm3, %%xmm3  \n\t"

	"vmovss         %%xmm0 ,  (%9)		\n\t"	// save temp2
	"vmovss         %%xmm1 , 4(%9)		\n\t"	// save temp2
	"vmovss         %%xmm2 , 8(%9)		\n\t"	// save temp2
	"vmovss         %%xmm3 ,12(%9)		\n\t"	// save temp2
	"vzeroupper				     \n\t"

	:
          "+r" (i),	// 0	
	  "+r" (n)  	// 1
        :
          "r" (x),      // 2
          "r" (y),      // 3
          "r" (a0),	// 4
          "r" (a1),	// 5
          "r" (a2),	// 6
          "r" (a3),	// 8
          "r" (temp1),  // 8
          "r" (temp2)   // 9
	: "cc", 
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);
}



void ssymv_kernel_8x1(const int32_t n,
                      float * __restrict a0,
		      float * __restrict xp,
		      float * __restrict yp,
		      float * __restrict temp1,
		      float * __restric temp2) {
        float at0,at1,at2,at3;
        float temp = 0.0f;
	float t1   = *temp1;
	const char pad[8];
	const int32_t len = (n/4)*4;
#if defined __GNUC__ && !defined __INTEL_COMPILER
        xp = (float*)__builtin_assume_aligned(xp,32);
	yp = (float*)__builtin_assume_aligned(yp,32);
#elif defined __INTEL_COMPILER
        __assume_aligned(xp,32);
	__assume_aligned(yp,32);
#endif
#if defined __ICC || defined __INTEL_COMPILER
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(xp:32,yp:32) safelen(4)
#endif
	for (i=0; i<len; i+=4) {
	
		at0     = a0[i];
		at1     = a0[i+1];
		at2     = a0[i+2];
		at3     = a0[i+3];

		yp[i]   += t1    * at0;
		temp    += at0   * xp[i];
		yp[i+1] += t1    * at1;
		temp    += at1   * xp[i+1];

		yp[i+2] += t1    * at2;
		temp    += at2   * xp[i+2];
		yp[i+3] += t1    * at3;
		temp    += at3   * xp[i+3];

	}
	*temp2 = temp;
}


void ssymv_kernel_1x4(const int32_t from,
                      const int32_t to,
		      float * __restrict a0,
		      float * __restrict a1,
		      float * __restrict a2,
		      float * __restrict a3,
		      float * __restrict xp,
		      float * __restrict yp,
		      float * __restrict temp1,
		      float * __restrict temp2) {
        float tmp2[4] = {};
	float at0,at1,at2,at3;
	float tp0,tp1,tp2,tp3;
	float x;
	const char pad[4];
	int32_t i;
	tp0 = temp1[0];
	tp1 = temp1[1];
	tp2 = temp1[2];
	tp3 = temp1[3];
#if defined __GNUC__ && !defined __INTEL_COMPILER
        xp = (float*)__builtin_assume_aligned(xp,64);
	yp = (float*)__builtin_assume_aligned(yp,64);
#elif defined __INTEL_COMPILER
        __assume_aligned(xp,64);
	__assume_aligned(yp,64);
#endif
#if defined __ICC || defined __INTEL_COMPILER
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(xp:32,yp:32)
#endif
      	for(i=from; i<to; i++) {
	
		at0     = a0[i];
		at1     = a1[i];
		at2     = a2[i];
		at3     = a3[i];
		x       = xp[i];
		yp[i]   += tp0 * at0 + tp1 *at1 + tp2 * at2 + tp3 * at3;
		tmp2[0] += at0 * x;
		tmp2[1] += at1 * x;
		tmp2[2] += at2 * x;
		tmp2[3] += at3 * x;

	}

	temp2[0] += tmp2[0];
	temp2[1] += tmp2[1];
	temp2[2] += tmp2[2];
	temp2[3] += tmp2[3];  
}



int32_t ssymv_u(const int32_t m,
                const int32_t offset,
		const float alpha,
		float * __restrict a,
		const int32_t lda,
		float * __restrict x,
		const int32_t inc_x,
		float * __restrict y,
		const int32_t inc_y,
		float * __restrict buffer) {

      float tmp1[4] = {};
      float tmp2[4] = {};
      float temp1,temp2;
      float * __restrict __ATTR_ALIGN__(64) xp = NULL;
      float * __restrict __ATTR_ALIGN__(64) yp = NULL;
      float * __restrict a0 = NULL;
      float * __restrict a1 = NULL;
      float * __restrict a2 = NULL;
      float * __restrict a3 = NULL;
      int32_t i,j,ix,jx,iy,jy;
      int32_t j1,j2,m2;
      const int32_t m1 = m-offset;
      const int32_t mrange = m - m1;
#if defined __GNUC__ && !defined __INTEL_COMPILER
        a = (float*)__builtin_assume_aligned(a,64);
        x = (float*)__builtin_assume_aligned(x,64);
	y = (float*)__builtin_assume_aligned(y,64);
#elif defined __INTEL_COMPILER
        __assume_aligned(a,64);
        __assume_aligned(x,64);
	__assume_aligned(y,64);
#endif
      if ( (inc_x!=1) || (inc_y!=1) || (mrange<16) ) {
	

		jx = m1 * inc_x;
		jy = m1 * inc_y;

		for (j=m1; j<m; j++)
		{
			temp1 = alpha * x[jx];
			temp2 = 0.0f;
			iy = 0;
			ix = 0;
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd reduction(+:y,temp2) aligned(y:64,a:64)
#elif defined __ICC || defined __INTEL_COMPILER
#pragma vector aligned
#pragma vector always
#endif
			for (i=0; i<j; i++)
			{
				y[iy] += temp1 * a[j*lda+i];
				temp2 += a[j*lda+i] * x[ix];
				ix += inc_x;
				iy += inc_y;
			
			}
			y[jy] += temp1 * a[j*lda+j] + alpha * temp2;
			jx    += inc_x;
			jy    += inc_y;
		}
		return(0);
	}

	xp = x;
	yp = y;

	m2 = m - ( mrange % 4 );

	for (j=m1; j<m2; j+=4)
	{
		tmp1[0] = alpha * xp[j];
		tmp1[1] = alpha * xp[j+1];
		tmp1[2] = alpha * xp[j+2];
		tmp1[3] = alpha * xp[j+3];
		tmp2[0] = 0.0;
		tmp2[1] = 0.0;
		tmp2[2] = 0.0;
		tmp2[3] = 0.0;
		a0    = &a[j*lda];
		a1    = a0+lda;
		a2    = a1+lda;
		a3    = a2+lda;
		j1 = (j/8)*8;		
		if ( j1 )
			ssymv_kernel_4x4(j1, a0, a1, a2, a3, xp, yp, tmp1, tmp2);
		if ( j1 < j )
			ssymv_kernel_1x4(j1, j,  a0, a1, a2, a3, xp, yp, tmp1, tmp2);

		j2 = 0;
		for ( j1 = j ; j1 < j+4 ; j1++ )
		{
			temp1 = tmp1[j2];
			temp2 = tmp2[j2];
			a0    = &a[j1*lda];
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd reduction(+:yp,temp2)  aligned(xp:64,yp:64)
#elif defined __ICC || defined __INTEL_COMPILER
#pragma vector aligned
#pragma vector always
#endif			
			for ( i=j ; i<j1; i++ )
			{
				yp[i] += temp1 * a0[i];	
				temp2 += a0[i] * xp[i];
				
			}
			y[j1] += temp1 * a0[j1] + alpha * temp2;
			j2++;

		}

	}

	for ( ; j<m; j++)
	{
		temp1 = alpha * xp[j];
		temp2 = 0.0f;
		a0    = &a[j*lda];
		float at0;
		j1 = (j/8)*8;		

		if ( j1 )
			ssymv_kernel_8x1(j1, a0, xp, yp, &temp1, &temp2);

#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd reduction(+:yp,temp2)  aligned(xp:64,yp:64)
#elif defined __ICC || defined __INTEL_COMPILER
#pragma vector aligned
#pragma vector always
#endif				
		for (i=j1 ; i<j; i++)
		{
			at0     = a0[i];
			yp[i] += temp1 * at0;
			temp2 += at0 * xp[i];
			
		}

		yp[j] += temp1 * a0[j] + alpha * temp2;
	}

	return(0);
	
}






