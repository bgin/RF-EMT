
/*******************************************************************************
*
*
*  Copyright [2019 - 2023] [Intel Corporation]
* 
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  
*  You may obtain a copy of the License at
*  
*     http://www.apache.org/licenses/LICENSE-2.0 
*  
*  Unless required by applicable law or agreed to in writing, software 
*  distributed under the License is distributed on an "AS IS" BASIS, 
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and 
*  limitations under the License. 
*  
*  SPDX-License-Identifier: Apache-2.0 
*  
* 
*
*******************************************************************************/

/*MIT License
!Copyright (c) 2020 Bernard Gingold
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
*/

#include "GMS_matrix_inv_cholesky.h"

// complex multiplication: A * A'
#define GET_AxAH(re, im, out)\
{\
    out = _mm512_fmadd_ps(re,re,_mm512_mul_ps(im,im));\
}

// complex multiplication: A * B'
#define GET_AxBH(are, aim, bre, bim, outre, outim)\
{\
    outre = _mm512_fmadd_ps(are,bre,_mm512_mul_ps(aim,bim));\
    outim = _mm512_fmsub_ps(bre,aim,_mm512_mul_ps(are,bim));\
}

// complex multiplication: A * B
#define GET_AxB(are, aim, bre, bim, outre, outim)\
{\
    outre = _mm512_fmsub_ps(are,bre,_mm512_mul_ps(aim,bim));\
    outim = _mm512_fmadd_ps(bre,aim,_mm512_mul_ps(are,bim));\
}

// complex multiplication: A * real(B)
#define GET_AxRealB(are, aim, bre, outre, outim)\
{\
    outre = _mm512_mul_ps(are, bre);\
    outim = _mm512_mul_ps(aim, bre);\
}

// get G00
#define GET_G00(matGRe, matBRe, matD, matND)\
{\
    matD[0] = _mm512_rsqrt14_ps(matBRe[0][0]);\
    matND[0] = _mm512_sub_ps(vzero, matD[0]);\
}

#define GET_G00_SAFE(matGRe, matBRe, matD, matND)\
{\
    __mmask16 k0 = _mm512_cmp_ps_mask(matBRe[0][0],vzero,_CMP_LT_OQ);\
    matD[0] = _mm512_rsqrt14_ps(_mm512_mask_blend_ps(k0,matBRe[0][0],_mm512_mul_ps(vneg_one,matBRe[0][0])));\
    matND[0] = _mm512_sub_ps(vzero, matD[0]);\
}

// get column 0 of matrix G
#define GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, j)\
{\
    matGRe[j][0] = _mm512_mul_ps(matBRe[j][0], matD[0]);\
    matGIm[j][0] = _mm512_mul_ps(matBIm[j][0], matD[0]);\
}

// get G11
#define GET_G11(matGRe, matGIm, matBRe, matD, matND, temp0)\
{\
    GET_AxAH(matGRe[1][0], matGIm[1][0], temp0);\
    matGRe[1][1] = _mm512_sub_ps(matBRe[1][1], temp0);\
    matD[1] = _mm512_rsqrt14_ps(matGRe[1][1]);\
    matND[1] = _mm512_sub_ps(vzero, matD[1]);\
}

// get G11
#define GET_G11_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0)\
{\
    GET_AxAH(matGRe[1][0], matGIm[1][0], temp0);\
    matGRe[1][1] = _mm512_sub_ps(matBRe[1][1], temp0);\
    __mmask16 k1 = _mm512_cmp_ps_mask(matGRe[1][1],vzero,_CMP_LT_OQ);\
    matD[1] = _mm512_rsqrt14_ps(_mm512_mask_blend_ps(k1,matGRe[1][1],_mm512_mul_ps(vneg_one,matGRe[1][1])));\
    matND[1] = _mm512_sub_ps(vzero, matD[1]);\
}


// This version albeit ~4.66 times slower than VRSQRT14 version achieves 
// minimize an approximation error of 2^-14 propagated by the usage of 
// VRSQRT14 machine code instruction.
#define GET_G11_MIN_APPROX_ERROR(matGRe, matGIm, matBRe, matD, matND, temp0)\
{\
    GET_AxAH(matGRe[1][0], matGIm[1][0], temp0);\
    matGRe[1][1] = _mm512_sub_ps(matBRe[1][1], temp0);\
    __mmask16 k1 = _mm512_cmp_ps_mask(matGRe[1][1],vzero,_CMP_LT_OQ);\
    matD[1] = _mm512_div_ps(vone,_mm512_sqrt_ps(_mm512_mask_blend_ps(k1,matGRe[1][1],_mm512_mul_ps(vneg_one,matGRe[1][1]))));\
    matND[1] = _mm512_sub_ps(vzero, matD[1]);\
}

// get column 1 of matrix G
#define GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, j, temp0, temp1)\
{\
    GET_AxBH(matGRe[j][0], matGIm[j][0], matGRe[1][0], matGIm[1][0], temp0, temp1);\
    matGRe[j][1] = _mm512_sub_ps(matBRe[j][1], temp0);\
    matGIm[j][1] = _mm512_sub_ps(matBIm[j][1], temp1);\
    matGRe[j][1] = _mm512_mul_ps(matGRe[j][1], matD[1]);\
    matGIm[j][1] = _mm512_mul_ps(matGIm[j][1], matD[1]);\
}

// get G22
#define GET_G22(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
    GET_AxAH(matGRe[2][0], matGIm[2][0], temp0);\
    GET_AxAH(matGRe[2][1], matGIm[2][1], temp1);\
    matGRe[2][2] = _mm512_sub_ps(matBRe[2][2], temp0);\
    matGRe[2][2] = _mm512_sub_ps(matGRe[2][2], temp1);\
    matD[2] = _mm512_rsqrt14_ps(matGRe[2][2]);\
    matND[2] = _mm512_sub_ps(vzero, matD[2]);\
}

// get G22
#define GET_G22_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
    GET_AxAH(matGRe[2][0], matGIm[2][0], temp0);\
    GET_AxAH(matGRe[2][1], matGIm[2][1], temp1);\
    matGRe[2][2] = _mm512_sub_ps(matBRe[2][2], temp0);\
    matGRe[2][2] = _mm512_sub_ps(matGRe[2][2], temp1);\
    __mmask16 k2 = _mm512_cmp_ps_mask(matGRe[2][2],vzero,_CMP_LT_OQ);\
    matD[2] = _mm512_rsqrt14_ps(_mm512_mask_blend_ps(k2,matGRe[2][2],_mm512_mul_ps(vneg_one,matGRe[2][2])));\
    matND[2] = _mm512_sub_ps(vzero, matD[2]);\
}


#define GET_G22_SAFE_MIN_APPROX_ERROR(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
    GET_AxAH(matGRe[2][0], matGIm[2][0], temp0);\
    GET_AxAH(matGRe[2][1], matGIm[2][1], temp1);\
    matGRe[2][2] = _mm512_sub_ps(matBRe[2][2], temp0);\
    matGRe[2][2] = _mm512_sub_ps(matGRe[2][2], temp1);\
    __mmask16 k2 = _mm512_cmp_ps_mask(matGRe[2][2],vzero,_CMP_LT_OQ);\
    matD[2] = _mm512_div_ps(vone,_mm512_sqrt_ps(_mm512_mask_blend_ps(k2,matGRe[2][2],_mm512_mul_ps(vneg_one,matGRe[2][2]))));\
    matND[2] = _mm512_sub_ps(vzero, matD[2]);\
}

// get column 2 of matrix G
#define GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, j, temp0, temp1)\
{\
    GET_AxBH(matGRe[j][0], matGIm[j][0], matGRe[2][0], matGIm[2][0], temp0, temp1);\
    matGRe[j][2] = _mm512_sub_ps(matBRe[j][2], temp0);\
    matGIm[j][2] = _mm512_sub_ps(matBIm[j][2], temp1);\
    GET_AxBH(matGRe[j][1], matGIm[j][1], matGRe[2][1], matGIm[2][1], temp0, temp1);\
    matGRe[j][2] = _mm512_sub_ps(matGRe[j][2], temp0);\
    matGIm[j][2] = _mm512_sub_ps(matGIm[j][2], temp1);\
    matGRe[j][2] = _mm512_mul_ps(matGRe[j][2], matD[2]);\
    matGIm[j][2] = _mm512_mul_ps(matGIm[j][2], matD[2]);\
}

// get G33
#define GET_G33(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
    GET_AxAH(matGRe[3][0], matGIm[3][0], temp0);\
    GET_AxAH(matGRe[3][1], matGIm[3][1], temp1);\
    GET_AxAH(matGRe[3][2], matGIm[3][2], temp2);\
    matGRe[3][3] = _mm512_sub_ps(matBRe[3][3], temp0);\
    matGRe[3][3] = _mm512_sub_ps(matGRe[3][3], temp1);\
    matGRe[3][3] = _mm512_sub_ps(matGRe[3][3], temp2);\
    matD[3] = _mm512_rsqrt14_ps(matGRe[3][3]);\
    matND[3] = _mm512_sub_ps(vzero, matD[3]);\
}

// get G33
#define GET_G33_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
    GET_AxAH(matGRe[3][0], matGIm[3][0], temp0);\
    GET_AxAH(matGRe[3][1], matGIm[3][1], temp1);\
    GET_AxAH(matGRe[3][2], matGIm[3][2], temp2);\
    matGRe[3][3] = _mm512_sub_ps(matBRe[3][3], temp0);\
    matGRe[3][3] = _mm512_sub_ps(matGRe[3][3], temp1);\
    matGRe[3][3] = _mm512_sub_ps(matGRe[3][3], temp2);\
    __mmask16 k3 = _mm512_cmp_ps_mask(matGRe[3][3],vzero,_CMP_LT_OQ);\
    matD[3] = _mm512_rsqrt14_ps(_mm512_mask_blend_ps(k3,matGRe[3][3],_mm512_mul_ps(vneg_one,matGRe[3][3])));\
    matND[3] = _mm512_sub_ps(vzero, matD[3]);\
}


#define GET_G33_SAFE_MIN_APPROX_ERROR(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
    GET_AxAH(matGRe[3][0], matGIm[3][0], temp0);\
    GET_AxAH(matGRe[3][1], matGIm[3][1], temp1);\
    GET_AxAH(matGRe[3][2], matGIm[3][2], temp2);\
    matGRe[3][3] = _mm512_sub_ps(matBRe[3][3], temp0);\
    matGRe[3][3] = _mm512_sub_ps(matGRe[3][3], temp1);\
    matGRe[3][3] = _mm512_sub_ps(matGRe[3][3], temp2);\
    __mmask16 k3 = _mm512_cmp_ps_mask(matGRe[3][3],vzero,_CMP_LT_OQ);\
    matD[3] = _mm512_div_ps(vone,_mm512_sqrt_ps(_mm512_mask_blend_ps(k3,matGRe[3][3],_mm512_mul_ps(vneg_one,matGRe[3][3]))));\
    matND[3] = _mm512_sub_ps(vzero, matD[3]);\
}

// get column 3 of matrix G
#define GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, j, temp0, temp1)\
{\
    GET_AxBH(matGRe[j][0], matGIm[j][0], matGRe[3][0], matGIm[3][0], temp0, temp1);\
    matGRe[j][3] = _mm512_sub_ps(matBRe[j][3], temp0);\
    matGIm[j][3] = _mm512_sub_ps(matBIm[j][3], temp1);\
    GET_AxBH(matGRe[j][1], matGIm[j][1], matGRe[3][1], matGIm[3][1], temp0, temp1);\
    matGRe[j][3] = _mm512_sub_ps(matGRe[j][3], temp0);\
    matGIm[j][3] = _mm512_sub_ps(matGIm[j][3], temp1);\
    GET_AxBH(matGRe[j][2], matGIm[j][2], matGRe[3][2], matGIm[3][2], temp0, temp1);\
    matGRe[j][3] = _mm512_sub_ps(matGRe[j][3], temp0);\
    matGIm[j][3] = _mm512_sub_ps(matGIm[j][3], temp1);\
    matGRe[j][3] = _mm512_mul_ps(matGRe[j][3], matD[3]);\
    matGIm[j][3] = _mm512_mul_ps(matGIm[j][3], matD[3]);\
}

#define GET_G44(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
    GET_AxAH(matGRe[4][0], matGIm[4][0], temp0);\
    GET_AxAH(matGRe[4][1], matGIm[4][1], temp1);\
    matGRe[4][4] = _mm512_sub_ps(matBRe[4][4], temp0);\
    matGRe[4][4] = _mm512_sub_ps(matGRe[4][4], temp1);\
    GET_AxAH(matGRe[4][2], matGIm[4][2], temp0);\
    GET_AxAH(matGRe[4][3], matGIm[4][3], temp1);\
    matGRe[4][4] = _mm512_sub_ps(matGRe[4][4], temp0);\
    matGRe[4][4] = _mm512_sub_ps(matGRe[4][4], temp1);\
    matD[4] = _mm512_rsqrt14_ps(matGRe[4][4]);\
    matND[4] = _mm512_sub_ps(vzero, matD[4]);\
}

// get G44
#define GET_G44_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
    GET_AxAH(matGRe[4][0], matGIm[4][0], temp0);\
    GET_AxAH(matGRe[4][1], matGIm[4][1], temp1);\
    matGRe[4][4] = _mm512_sub_ps(matBRe[4][4], temp0);\
    matGRe[4][4] = _mm512_sub_ps(matGRe[4][4], temp1);\
    GET_AxAH(matGRe[4][2], matGIm[4][2], temp0);\
    GET_AxAH(matGRe[4][3], matGIm[4][3], temp1);\
    matGRe[4][4] = _mm512_sub_ps(matGRe[4][4], temp0);\
    matGRe[4][4] = _mm512_sub_ps(matGRe[4][4], temp1);\
    __mmask16 k4 = _mm512_cmp_ps_mask(matGRe[4][4],vzero,_CMP_LT_OQ);\
    matD[4] = _mm512_rsqrt14_ps(_mm512_mask_blend_ps(k4,matGRe[4][4],_mm512_mul_ps(vneg_one,matGRe[4][4])));\
    matND[4] = _mm512_sub_ps(vzero, matD[4]);\
}


#define GET_G44_SAFE_MIN_APPROX_ERROR(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
    GET_AxAH(matGRe[4][0], matGIm[4][0], temp0);\
    GET_AxAH(matGRe[4][1], matGIm[4][1], temp1);\
    matGRe[4][4] = _mm512_sub_ps(matBRe[4][4], temp0);\
    matGRe[4][4] = _mm512_sub_ps(matGRe[4][4], temp1);\
    GET_AxAH(matGRe[4][2], matGIm[4][2], temp0);\
    GET_AxAH(matGRe[4][3], matGIm[4][3], temp1);\
    matGRe[4][4] = _mm512_sub_ps(matGRe[4][4], temp0);\
    matGRe[4][4] = _mm512_sub_ps(matGRe[4][4], temp1);\
    __mmask16 k4 = _mm512_cmp_ps_mask(matGRe[4][4],vzero,_CMP_LT_OQ);\
    matD[4] = _mm512_div_ps(vone,_mm512_sqrt_ps(_mm512_mask_blend_ps(k4,matGRe[4][4],_mm512_mul_ps(vneg_one,matGRe[4][4]))));\
    matND[4] = _mm512_sub_ps(vzero, matD[4]);\
}


#define GET_G55(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
        GET_AxAH(matGRe[5][0], matGIm[5][0], temp0);\
        GET_AxAH(matGRe[5][1], matGIm[5][1], temp1);\
        matGRe[5][5] = _mm512_sub_ps(matBRe[5][5], temp0);\
        matGRe[5][5] = _mm512_sub_ps(matGRe[5][5], temp1);\
        GET_AxAH(matGRe[5][2], matGIm[5][2], temp0);\
        GET_AxAH(matGRe[5][3], matGIm[5][3], temp1);\
        GET_AxAH(matGRe[5][4], matGIm[5][4], temp2);\
        matGRe[5][5] = _mm512_sub_ps(matGRe[5][5], temp0);\
        matGRe[5][5] = _mm512_sub_ps(matGRe[5][5], temp1);\
        matGRe[5][5] = _mm512_sub_ps(matGRe[5][5], temp2);\
        matD[5] = _mm512_rsqrt14_ps(matGRe[5][5]);\
        matND[5] = _mm512_sub_ps(vzero, matD[5]);\
}

// get G55
#define GET_G55_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
        GET_AxAH(matGRe[5][0], matGIm[5][0], temp0);\
        GET_AxAH(matGRe[5][1], matGIm[5][1], temp1);\
        matGRe[5][5] = _mm512_sub_ps(matBRe[5][5], temp0);\
        matGRe[5][5] = _mm512_sub_ps(matGRe[5][5], temp1);\
        GET_AxAH(matGRe[5][2], matGIm[5][2], temp0);\
        GET_AxAH(matGRe[5][3], matGIm[5][3], temp1);\
        GET_AxAH(matGRe[5][4], matGIm[5][4], temp2);\
        matGRe[5][5] = _mm512_sub_ps(matGRe[5][5], temp0);\
        matGRe[5][5] = _mm512_sub_ps(matGRe[5][5], temp1);\
        matGRe[5][5] = _mm512_sub_ps(matGRe[5][5], temp2);\
        __mmask16 k5 = _mm512_cmp_ps_mask(matGRe[5][5],vzero,_CMP_LT_OQ);\
        matD[5] = _mm512_rsqrt14_ps(_mm512_mask_blend_ps(k5,matGRe[5][5],_mm512_mul_ps(vneg_one,matGRe[5][5])));\
        matND[5] = _mm512_sub_ps(vzero, matD[5]);\
}


#define GET_G55_SAFE_MIN_APPROX_ERROR(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
        GET_AxAH(matGRe[5][0], matGIm[5][0], temp0);\
        GET_AxAH(matGRe[5][1], matGIm[5][1], temp1);\
        matGRe[5][5] = _mm512_sub_ps(matBRe[5][5], temp0);\
        matGRe[5][5] = _mm512_sub_ps(matGRe[5][5], temp1);\
        GET_AxAH(matGRe[5][2], matGIm[5][2], temp0);\
        GET_AxAH(matGRe[5][3], matGIm[5][3], temp1);\
        GET_AxAH(matGRe[5][4], matGIm[5][4], temp2);\
        matGRe[5][5] = _mm512_sub_ps(matGRe[5][5], temp0);\
        matGRe[5][5] = _mm512_sub_ps(matGRe[5][5], temp1);\
        matGRe[5][5] = _mm512_sub_ps(matGRe[5][5], temp2);\
        __mmask16 k5 = _mm512_cmp_ps_mask(matGRe[5][5],vzero,_CMP_LT_OQ);\
        matD[5] = _mm512_div_ps(vone,_mm512_sqrt_ps(_mm512_mask_blend_ps(k5,matGRe[5][5],_mm512_mul_ps(vneg_one,matGRe[5][5]))));\
        matND[5] = _mm512_sub_ps(vzero, matD[5]);\
}


#define GET_G66(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
        GET_AxAH(matGRe[6][0], matGIm[6][0], temp0);\
        GET_AxAH(matGRe[6][1], matGIm[6][1], temp1);\
        matGRe[6][6] = _mm512_sub_ps(matBRe[6][6], temp0);\
        matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp1);\
        GET_AxAH(matGRe[6][2], matGIm[6][2], temp0);\
        GET_AxAH(matGRe[6][3], matGIm[6][3], temp1);\
        matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp0);\
        matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp1);\
        GET_AxAH(matGRe[6][4], matGIm[6][4], temp0);\
        GET_AxAH(matGRe[6][5], matGIm[6][5], temp1);\
        matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp0);\
        matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp1);\
        matD[6] = _mm512_rsqrt14_ps(matGRe[6][6]);\
        matND[6] = _mm512_sub_ps(vzero, matD[6]);\
}

// get G66
#define GET_G66_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
        GET_AxAH(matGRe[6][0], matGIm[6][0], temp0);\
        GET_AxAH(matGRe[6][1], matGIm[6][1], temp1);\
        matGRe[6][6] = _mm512_sub_ps(matBRe[6][6], temp0);\
        matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp1);\
        GET_AxAH(matGRe[6][2], matGIm[6][2], temp0);\
        GET_AxAH(matGRe[6][3], matGIm[6][3], temp1);\
        matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp0);\
        matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp1);\
        GET_AxAH(matGRe[6][4], matGIm[6][4], temp0);\
        GET_AxAH(matGRe[6][5], matGIm[6][5], temp1);\
        matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp0);\
        matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp1);\
        __mmask16 k6 = _mm512_cmp_ps_mask(matGRe[6][6],vzero,_CMP_LT_OQ);\
        matD[6] = _mm512_rsqrt14_ps(_mm512_mask_blend_ps(k6,matGRe[6][6],_mm512_mul_ps(vneg_one,matGRe[6][6])));\
        matND[6] = _mm512_sub_ps(vzero, matD[6]);\
}


#define GET_G66_SAFE_MIN_APPROX_ERROR(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
        GET_AxAH(matGRe[6][0], matGIm[6][0], temp0);\
        GET_AxAH(matGRe[6][1], matGIm[6][1], temp1);\
        matGRe[6][6] = _mm512_sub_ps(matBRe[6][6], temp0);\
        matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp1);\
        GET_AxAH(matGRe[6][2], matGIm[6][2], temp0);\
        GET_AxAH(matGRe[6][3], matGIm[6][3], temp1);\
        matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp0);\
        matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp1);\
        GET_AxAH(matGRe[6][4], matGIm[6][4], temp0);\
        GET_AxAH(matGRe[6][5], matGIm[6][5], temp1);\
        matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp0);\
        matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp1);\
        __mmask16 k6 = _mm512_cmp_ps_mask(matGRe[6][6],vzero,_CMP_LT_OQ);\
        matD[6] = _mm512_div_ps(vone,_mm512_sqrt_ps(_mm512_mask_blend_ps(k6,matGRe[6][6],_mm512_mul_ps(vneg_one,matGRe[6][6]))));\
        matND[6] = _mm512_sub_ps(vzero, matD[6]);\
}

#define GET_G77(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
        GET_AxAH(matGRe[7][0], matGIm[7][0], temp0);\
        GET_AxAH(matGRe[7][1], matGIm[7][1], temp1);\
        matGRe[7][7] = _mm512_sub_ps(matBRe[7][7], temp0);\
        matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp1);\
        GET_AxAH(matGRe[7][2], matGIm[7][2], temp0);\
        GET_AxAH(matGRe[7][3], matGIm[7][3], temp1);\
        matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp0);\
        matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp1);\
        GET_AxAH(matGRe[7][4], matGIm[7][4], temp0);\
        GET_AxAH(matGRe[7][5], matGIm[7][5], temp1);\
        GET_AxAH(matGRe[7][6], matGIm[7][6], temp2);\
        matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp0);\
        matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp1);\
        matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp2);\
        matD[7] = _mm512_rsqrt14_ps(matGRe[7][7]);\
        matND[7] = _mm512_sub_ps(vzero, matD[7]);\
}

// get G77
#define GET_G77_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
        GET_AxAH(matGRe[7][0], matGIm[7][0], temp0);\
        GET_AxAH(matGRe[7][1], matGIm[7][1], temp1);\
        matGRe[7][7] = _mm512_sub_ps(matBRe[7][7], temp0);\
        matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp1);\
        GET_AxAH(matGRe[7][2], matGIm[7][2], temp0);\
        GET_AxAH(matGRe[7][3], matGIm[7][3], temp1);\
        matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp0);\
        matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp1);\
        GET_AxAH(matGRe[7][4], matGIm[7][4], temp0);\
        GET_AxAH(matGRe[7][5], matGIm[7][5], temp1);\
        GET_AxAH(matGRe[7][6], matGIm[7][6], temp2);\
        matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp0);\
        matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp1);\
        matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp2);\
        __mmask16 k7 = _mm512_cmp_ps_mask(matGRe[7][7],vzero,_CMP_LT_OQ);\
        matD[7] = _mm512_rsqrt14_ps(_mm512_mask_blend_ps(k7,matGRe[7][7],_mm512_mul_ps(vneg_one,matGRe[7][7])));\
        matND[7] = _mm512_sub_ps(vzero, matD[7]);\
}


#define GET_G77_SAFE_MIN_APPROX_ERROR(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
        GET_AxAH(matGRe[7][0], matGIm[7][0], temp0);\
        GET_AxAH(matGRe[7][1], matGIm[7][1], temp1);\
        matGRe[7][7] = _mm512_sub_ps(matBRe[7][7], temp0);\
        matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp1);\
        GET_AxAH(matGRe[7][2], matGIm[7][2], temp0);\
        GET_AxAH(matGRe[7][3], matGIm[7][3], temp1);\
        matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp0);\
        matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp1);\
        GET_AxAH(matGRe[7][4], matGIm[7][4], temp0);\
        GET_AxAH(matGRe[7][5], matGIm[7][5], temp1);\
        GET_AxAH(matGRe[7][6], matGIm[7][6], temp2);\
        matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp0);\
        matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp1);\
        matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp2);\
        __mmask16 k7 = _mm512_cmp_ps_mask(matGRe[7][7],vzero,_CMP_LT_OQ);\
        matD[7] = _mm512_div_ps(vone,_mm512_sqrt_ps(_mm512_mask_blend_ps(k7,matGRe[7][7],_mm512_mul_ps(vneg_one,matGRe[7][7]))));\
        matND[7] = _mm512_sub_ps(vzero, matD[7]);\
}


#define GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, i)\
{\
        GET_AxAH(matGRe[i][0], matGIm[i][0], temp0);\
        matGRe[i][i] = _mm512_sub_ps(matBRe[i][i], temp0);\
        for (std::int32_t i1 = 1; i1 < i; i1+=2) \
        {\
           GET_AxAH(matGRe[i][i1], matGIm[i][i1], temp0);\
           GET_AxAH(matGRe[i][i1+1], matGIm[i][i1+1], temp1);\
           matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp0);\
           matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp1);\
        }\
        matD[i] = _mm512_rsqrt14_ps(matGRe[i][i]);\
        matND[i] = _mm512_sub_ps(vzero, matD[i]);\
}

// get Gii, odd diagonal element
#define GET_Gii_ODD_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, i)\
{\
        GET_AxAH(matGRe[i][0], matGIm[i][0], temp0);\
        matGRe[i][i] = _mm512_sub_ps(matBRe[i][i], temp0);\
        for (std::int32_t i1 = 1; i1 < i; i1+=2) \
        {\
           GET_AxAH(matGRe[i][i1], matGIm[i][i1], temp0);\
           GET_AxAH(matGRe[i][i1+1], matGIm[i][i1+1], temp1);\
           matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp0);\
           matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp1);\
        }\
        __mmask16 kii_odd = _mm512_cmp_ps_mask(matGRe[i][i],vzero,_CMP_LT_OQ);\
        matD[i] = _mm512_rsqrt14_ps(_mm512_mask_blend_ps(kii_odd,matGRe[i][i],_mm512_mul_ps(vneg_one,matGRe[i][i])));\
        matND[i] = _mm512_sub_ps(vzero, matD[i]);\
}


#define GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, i)\
{\
        GET_AxAH(matGRe[i][0], matGIm[i][0], temp0);\
        matGRe[i][i] = _mm512_sub_ps(matBRe[i][i], temp0);\
        for (std::int32_t i1 = 1; i1 < i; i1+=2) \
        {\
           GET_AxAH(matGRe[i][i1], matGIm[i][i1], temp0);\
           GET_AxAH(matGRe[i][i1+1], matGIm[i][i1+1], temp1);\
           matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp0);\
           matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp1);\
        }\
        matD[i] = _mm512_rsqrt14_ps(matGRe[i][i]);\
        matND[i] = _mm512_sub_ps(vzero, matD[i]);\
}

#define GET_Gii_ODD_SAFE_MIN_APPROX_ERROR(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, i)\
{\
        GET_AxAH(matGRe[i][0], matGIm[i][0], temp0);\
        matGRe[i][i] = _mm512_sub_ps(matBRe[i][i], temp0);\
        for (std::int32_t i1 = 1; i1 < i; i1+=2) \
        {\
           GET_AxAH(matGRe[i][i1], matGIm[i][i1], temp0);\
           GET_AxAH(matGRe[i][i1+1], matGIm[i][i1+1], temp1);\
           matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp0);\
           matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp1);\
        }\
        __mmask16 kii_odd = _mm512_cmp_ps_mask(matGRe[i][i],vzero,_CMP_LT_OQ);\
        matD[i] = _mm512_div_ps(vone,_mm512_sqrt_ps(_mm512_mask_blend_ps(kii_odd,matGRe[i][i],_mm512_mul_ps(vneg_one,matGRe[i][i]))));\
        matND[i] = _mm512_sub_ps(vzero, matD[i]);\
}

#define GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, i)\
{\
        GET_AxAH(matGRe[i][0], matGIm[i][0], temp0);\
        GET_AxAH(matGRe[i][1], matGIm[i][1], temp1);\
        matGRe[i][i] = _mm512_sub_ps(matBRe[i][i], temp0);\
        matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp1);\
        for (std::int32_t i1 = 2; i1 < i; i1+=2) \
        {\
           GET_AxAH(matGRe[i][i1], matGIm[i][i1], temp0);\
           GET_AxAH(matGRe[i][i1+1], matGIm[i][i1+1], temp1);\
           matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp0);\
           matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp1);\
        }\
        matD[i] = _mm512_rsqrt14_ps(matGRe[i][i]);\
        matND[i] = _mm512_sub_ps(vzero, matD[i]);\
}

// get Gii, even diagonal element
#define GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, i)\
{\
        GET_AxAH(matGRe[i][0], matGIm[i][0], temp0);\
        GET_AxAH(matGRe[i][1], matGIm[i][1], temp1);\
        matGRe[i][i] = _mm512_sub_ps(matBRe[i][i], temp0);\
        matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp1);\
        for (std::int32_t i1 = 2; i1 < i; i1+=2) \
        {\
           GET_AxAH(matGRe[i][i1], matGIm[i][i1], temp0);\
           GET_AxAH(matGRe[i][i1+1], matGIm[i][i1+1], temp1);\
           matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp0);\
           matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp1);\
        }\
        __mmask16 kii_even = _mm512_cmp_ps_mask(matGRe[i][i],vzero,_CMP_LT_OQ);\
        matD[i] = _mm512_rsqrt14_ps(_mm512_mask_blend_ps(kii_even,matGRe[i][i],_mm512_mul_ps(vneg_one,matGRe[i][i])));\
        matND[i] = _mm512_sub_ps(vzero, matD[i]);\
}


#define GET_Gii_EVEN_SAFE_MIN_APPROX_ERROR(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, i)\
{\
        GET_AxAH(matGRe[i][0], matGIm[i][0], temp0);\
        GET_AxAH(matGRe[i][1], matGIm[i][1], temp1);\
        matGRe[i][i] = _mm512_sub_ps(matBRe[i][i], temp0);\
        matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp1);\
        for (std::int32_t i1 = 2; i1 < i; i1+=2) \
        {\
           GET_AxAH(matGRe[i][i1], matGIm[i][i1], temp0);\
           GET_AxAH(matGRe[i][i1+1], matGIm[i][i1+1], temp1);\
           matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp0);\
           matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp1);\
        }\
        __mmask16 kii_even = _mm512_cmp_ps_mask(matGRe[i][i],vzero,_CMP_LT_OQ);\
        matD[i] = _mm512_div_ps(vone,_mm512_sqrt_ps(_mm512_mask_blend_ps(kii_even,matGRe[i][i],_mm512_mul_ps(vneg_one,matGRe[i][i]))));\
        matND[i] = _mm512_sub_ps(vzero, matD[i]);\
}


//get odd column n of matrix G, j is row index, i is col index
#define GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, j, i, temp0, temp1)\
{\
        GET_AxBH(matGRe[j][0], matGIm[j][0], matGRe[i][0], matGIm[i][0], temp0, temp1);\
        matGRe[j][i] = _mm512_sub_ps(matBRe[j][i], temp0);\
        matGIm[j][i] = _mm512_sub_ps(matBIm[j][i], temp1);\
        for (std::int32_t i1 = 1; i1 < i; i1+=2)\
        {\
           GET_AxBH(matGRe[j][i1], matGIm[j][i1], matGRe[i][i1], matGIm[i][i1], temp0, temp1);\
           matGRe[j][i] = _mm512_sub_ps(matGRe[j][i], temp0);\
           matGIm[j][i] = _mm512_sub_ps(matGIm[j][i], temp1);\
           GET_AxBH(matGRe[j][i1+1], matGIm[j][i1+1], matGRe[i][i1+1], matGIm[i][i1+1], temp0, temp1);\
           matGRe[j][i] = _mm512_sub_ps(matGRe[j][i], temp0);\
           matGIm[j][i] = _mm512_sub_ps(matGIm[j][i], temp1);\
        }\
        matGRe[j][i] = _mm512_mul_ps(matGRe[j][i], matD[i]);\
        matGIm[j][i] = _mm512_mul_ps(matGIm[j][i], matD[i]);\
}

// get even column n of matrix G, j is row index, i is col index
#define GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, j, i, temp0, temp1)\
{\
        GET_AxBH(matGRe[j][0], matGIm[j][0], matGRe[i][0], matGIm[i][0], temp0, temp1);\
        matGRe[j][i] = _mm512_sub_ps(matBRe[j][i], temp0);\
        matGIm[j][i] = _mm512_sub_ps(matBIm[j][i], temp1);\
        GET_AxBH(matGRe[j][1], matGIm[j][1], matGRe[i][1], matGIm[i][1], temp0, temp1);\
        matGRe[j][i] = _mm512_sub_ps(matGRe[j][i], temp0);\
        matGIm[j][i] = _mm512_sub_ps(matGIm[j][i], temp1);\
        for (std::int32_t i1 = 2; i1 < i; i1+=2)\
        {\
            GET_AxBH(matGRe[j][i1], matGIm[j][i1], matGRe[i][i1], matGIm[i][i1], temp0, temp1);\
            matGRe[j][i] = _mm512_sub_ps(matGRe[j][i], temp0);\
            matGIm[j][i] = _mm512_sub_ps(matGIm[j][i], temp1);\
            GET_AxBH(matGRe[j][i1+1], matGIm[j][i1+1], matGRe[i][i1+1], matGIm[i][i1+1], temp0, temp1);\
            matGRe[j][i] = _mm512_sub_ps(matGRe[j][i], temp0);\
            matGIm[j][i] = _mm512_sub_ps(matGIm[j][i], temp1);\
        }\
        matGRe[j][i] = _mm512_mul_ps(matGRe[j][i], matD[i]);\
        matGIm[j][i] = _mm512_mul_ps(matGIm[j][i], matD[i]);\
}

// set value for Lii, diagonal element
#define SET_Lii(matLRe, matLIm, matD, i)\
{\
    matLRe[i][i] = matD[i];\
    matLIm[i][i] = vzero;\
}

// get element L(i+1, i)
#define GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, j, i)\
{\
    GET_AxRealB(matGRe[j][i], matGIm[j][i], matLRe[i][i],  matLRe[j][i], matLIm[j][i]);\
    matLRe[j][i] = _mm512_mul_ps(matLRe[j][i], matND[j]);\
    matLIm[j][i] = _mm512_mul_ps(matLIm[j][i], matND[j]);\
}

#define GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, j, i, temp0, temp1)\
{\
    GET_AxRealB(matGRe[j][i], matGIm[j][i], matLRe[i][i],  matLRe[j][i], matLIm[j][i]);\
        for (std::int32_t i1 = i+1; i1 < j; i1++)\
        {\
            GET_AxB(matGRe[j][i1], matGIm[j][i1], matLRe[i1][i], matLIm[i1][i], temp0, temp1);\
            matLRe[j][i] = _mm512_add_ps(matLRe[j][i], temp0);\
            matLIm[j][i] = _mm512_add_ps(matLIm[j][i], temp1);\
        }\
    matLRe[j][i] = _mm512_mul_ps(matLRe[j][i], matND[j]);\
    matLIm[j][i] = _mm512_mul_ps(matLIm[j][i], matND[j]);\
}

template<bool use_prefetching,bool mitigate_nan>
void 
gms::math
::mat_inv_cholesky_2x2_16xf32(const __m512 matBRe[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                              const __m512 matBIm[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                              __m512       matInvBRe[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                              __m512       matInvBIm[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2])
{
    __ATTR_ALIGN__(64) 
    __m512 matGRe[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2];
    __ATTR_ALIGN__(64)
    __m512 matGIm[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2];
    __ATTR_ALIGN__(64)
    __m512 matLRe[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2];
    __ATTR_ALIGN__(64)
    __m512 matLIm[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2];
    __ATTR_ALIGN__(64)
    __m512 matD[MAT_SQR_SIZE_2];
    __ATTR_ALIGN__(64)
    __m512 matND[MAT_SQR_SIZE_2];
    __m512 temp0, temp1;
    const __m512 vzero{_mm512_setzero_ps()};
    const __m512 vneg_one{_mm512_set1_ps(-1.0f)};
    //////////////////////////////////////////////////////////////
    if constexpr(static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
    {
        const float * __restrict ptr_matBRe{reinterpret_cast<const float * __restrict>(&matBRe)};
        const float * __restrict ptr_matBIm{reinterpret_cast<const float * __restrict>(&matBIm)};
        _mm_prefetch((const char*)&ptr_matBRe[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBIm[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBRe[16],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBIm[16],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBRe[32],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBIm[32],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBRe[48],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBIm[48],_MM_HINT_T0);
    }

   
    GET_G00(matGRe,matBRe,matD,matND);
    GET_G_COL0(matGRe,matGIm,matBRe,matBIm,matD,1);
    if constexpr(static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true))
    {
        GET_G11_SAFE(matGRe,matGIm,matBRe,matD,matND,temp0);
    }
    else 
    {
        GET_G11(matGRe,matGIm,matBRe,matD,matND,temp0);
    }
    
    SET_Lii(matLRe,matLIm,matD,0);
    GET_L_i1i(matLRe,matLIm,matGRe,matGIm,matND,1,0);
    SET_Lii(matLRe,matLIm,matD,1);

    /////////////////////////////////// get invB = L'*L
    __m512 matLRe00 = matLRe[0][0];
    __m512 matLRe10 = matLRe[1][0];
    __m512 matLRe11 = matLRe[1][1];
    __m512 matLIm10 = matLIm[1][0];
    matInvBRe[0][0] = _mm512_mul_ps(matLRe00,matLRe00);
    temp1           = _mm512_fmadd_ps(matLRe10,matLRe10,_mm512_mul_ps(matLIm10,matLIm10));
    matInvBRe[0][0] = _mm512_add_ps(matInvBRe[0][0], temp1);
    matInvBIm[0][0] = vzero;
    matInvBRe[1][1] = _mm512_mul_ps(matLRe11,matLRe11);
    matInvBIm[1][1] = vzero;
    matInvBRe[0][1] = _mm512_mul_ps(matLRe10,matLRe11);
    matInvBIm[0][1] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm10,matLRe11));
    matInvBRe[1][0] = matInvBRe[0][1];
    matInvBIm[1][0] = _mm512_sub_ps(vzero,matInvBIm[0][1]);
}

template void 
gms::math 
::mat_inv_cholesky_2x2_16xf32<true,true>(const __m512 matBRe[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                                         const __m512 matBIm[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                                         __m512       matInvBRe[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                                         __m512       matInvBIm[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2]);

template void 
gms::math 
::mat_inv_cholesky_2x2_16xf32<false,true>(const __m512 matBRe[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                                         const __m512 matBIm[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                                         __m512       matInvBRe[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                                         __m512       matInvBIm[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2]);

template void 
gms::math 
::mat_inv_cholesky_2x2_16xf32<true,false>(const __m512 matBRe[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                                         const __m512 matBIm[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                                         __m512       matInvBRe[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                                         __m512       matInvBIm[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2]);

template void 
gms::math 
::mat_inv_cholesky_2x2_16xf32<false,false>(const __m512 matBRe[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                                           const __m512 matBIm[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                                           __m512       matInvBRe[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                                           __m512       matInvBIm[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2]);

template<bool use_prefetching,bool mitigate_nan>
void 
gms::math
::mat_inv_cholesky_3x3_16xf32(const __m512 matBRe[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                              const __m512 matBIm[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                              __m512       matInvBRe[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                              __m512       matInvBIm[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3])
{
     __ATTR_ALIGN__(64) 
    __m512 matGRe[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3];
    __ATTR_ALIGN__(64)
    __m512 matGIm[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3];
    __ATTR_ALIGN__(64)
    __m512 matLRe[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3];
    __ATTR_ALIGN__(64)
    __m512 matLIm[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3];
    __ATTR_ALIGN__(64)
    __m512 matD[MAT_SQR_SIZE_3];
    __ATTR_ALIGN__(64)
    __m512 matND[MAT_SQR_SIZE_3];
    __m512 temp0, temp1, temp2;
    __m512 matLReii;
    __m512 matLReki;
    __m512 matLImki;
    std::int32_t i,j,k;
    const __m512 vzero{_mm512_setzero_ps()};
    const __m512 vneg_one{_mm512_set1_ps(-1.0f)};
    //////////////////////////////////////////////////////////////
    if constexpr (static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
    {
        const float * __restrict ptr_matBRe{reinterpret_cast<const float * __restrict>(&matBRe)};
        const float * __restrict ptr_matBIm{reinterpret_cast<const float * __restrict>(&matBIm)};
        _mm_prefetch((const char*)&ptr_matBRe[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBIm[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBRe[16],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBIm[16],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBRe[32],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBIm[32],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBRe[48],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBIm[48],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBRe[64],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBIm[64],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBRe[80],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBIm[80],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBRe[96],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBIm[96],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBRe[112],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBIm[112],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBRe[128],_MM_HINT_T0);
        _mm_prefetch((const char*)&ptr_matBIm[128],_MM_HINT_T0);
    }
    
   
    GET_G00(matGRe,matBRe,matD,matND);
    GET_G_COL0(matGRe,matGIm,matBRe,matBIm,matD,1);
    GET_G_COL0(matGRe,matGIm,matBRe,matBIm,matD,2);
    if constexpr(static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true))
    {
       GET_G11_SAFE(matGRe,matGIm,matBRe,matD,matND,temp0);
    }
    else 
    {
       GET_G11(matGRe,matGIm,matBRe,matD,matND,temp0);
    }
    GET_G_COL1(matGRe,matGIm,matBRe,matBIm,matD,2,temp0,temp1);
    if constexpr(static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true))
    {
       GET_G22_SAFE(matGRe,matGIm,matBRe,matD,matND,temp0,temp1);
    }
    else 
    {
       GET_G22(matGRe,matGIm,matBRe,matD,matND,temp0,temp1);
    }
    SET_Lii(matLRe,matLIm,matD,0);
    GET_L_i1i(matLRe,matLIm,matGRe,matGIm,matND,1,0);
    GET_L_ji(matLRe,matLIm,matGRe,matGIm,matND,2,0,temp0,temp1);
    // Column 1
    SET_Lii(matLRe,matLIm,matD,1);
    GET_L_i1i(matLRe,matLIm,matGRe,matGIm,matND,2,1);
    // Column 2
    SET_Lii(matLRe,matLIm,matD,2);
    
    for(i = 0; i < MAT_SQR_SIZE_3; ++i) 
    {
        matLReii        = matLRe[i][i];
        matInvBRe[i][i] = _mm512_mul_ps(matLReii,matLReii);
        for(k = (i+1); k < MAT_SQR_SIZE_3; ++k) 
        {
            matLReki        = matLRe[k][i];
            matLImki        = matLIm[k][i];
            temp1           = _mm512_fmadd_ps(matLReki,matLReki,_mm512_mul_ps(matLImki,matLImki));
            matInvBRe[i][i] = _mm512_add_ps(matInvBRe[i][i],temp1);
        }
        matInvBIm[i][i] = vzero;
    }

    for(i = 0; i < MAT_SQR_SIZE_3; ++i) 
    {
        for(j = (i+1); j < MAT_SQR_SIZE_3; ++j)  
        {
            matInvBRe[i][j] = _mm512_mul_ps(matLRe[j][i],matLRe[j][j]);
            matInvBIm[i][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][i],matLRe[j][j]));

            for(k = (j+1); k < MAT_SQR_SIZE_3; ++k) 
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][i],matLIm[k][i],temp1,temp2);
                matInvBRe[i][j] = _mm512_add_ps(matInvBRe[i][j], temp1);
                matInvBIm[i][j] = _mm512_add_ps(matInvBIm[i][j], temp2);
            }
            matInvBRe[j][i] = matInvBRe[i][j];
            matInvBIm[j][i] = _mm512_sub_ps(vzero,matInvBIm[i][j]);
        }
    }
}

template void 
gms::math 
::mat_inv_cholesky_3x3_16xf32<true,true>(const __m512 matBRe[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                                         const __m512 matBIm[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                                         __m512       matInvBRe[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                                         __m512       matInvBIm[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3]);

template void 
gms::math 
::mat_inv_cholesky_3x3_16xf32<false,true>(const __m512 matBRe[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                                         const __m512 matBIm[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                                         __m512       matInvBRe[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                                         __m512       matInvBIm[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3]);

template void 
gms::math 
::mat_inv_cholesky_3x3_16xf32<true,false>(const __m512 matBRe[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                                         const __m512 matBIm[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                                         __m512       matInvBRe[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                                         __m512       matInvBIm[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3]);

template void 
gms::math 
::mat_inv_cholesky_3x3_16xf32<false,false>(const __m512 matBRe[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                                           const __m512 matBIm[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                                           __m512       matInvBRe[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                                           __m512       matInvBIm[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3]);


template<bool use_prefetching,bool mitigate_nan>
void 
gms::math
::mat_inv_cholesky_4x4_16xf32(const __m512 matBRe[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                              const __m512 matBIm[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                              __m512       matInvBRe[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                              __m512       matInvBIm[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4])
{
    __ATTR_ALIGN__(64) 
    __m512 matGRe[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4];
    __ATTR_ALIGN__(64)
    __m512 matGIm[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4];
    __ATTR_ALIGN__(64)
    __m512 matLRe[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4];
    __ATTR_ALIGN__(64)
    __m512 matLIm[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4];
    __ATTR_ALIGN__(64)
    __m512 matD[MAT_SQR_SIZE_4];
    __ATTR_ALIGN__(64)
    __m512 matND[MAT_SQR_SIZE_4];
    __m512 temp0, temp1, temp2;
    __m512 matLReii;
    __m512 matLReki;
    __m512 matLImki;
    __m512 matLRejj;
    std::int32_t i,j,k;
    const __m512 vzero{_mm512_setzero_ps()};
    const __m512 vneg_one{_mm512_set1_ps(-1.0f)};
    //////////////////////////////////////////////////////////////
    if constexpr (static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
    {
        constexpr std::int32_t VEC_PS_LEN = 16;
        const float * __restrict ptr_matBRe{reinterpret_cast<const float * __restrict>(&matBRe)};
        const float * __restrict ptr_matBIm{reinterpret_cast<const float * __restrict>(&matBIm)};
        for(i = 0;i != MAT_SQR_SIZE_4; ++i) 
        {
            register std::int32_t outer_idx = i*MAT_SQR_SIZE_4;
            for(j = 0;j != MAT_SQR_SIZE_4; ++j)   
            {
                register std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                _mm_prefetch((const char*)&ptr_matBRe[inner_idx],_MM_HINT_T0);
                _mm_prefetch((const char*)&ptr_matBIm[inner_idx],_MM_HINT_T0);
            }
        }
    }
    
    if constexpr(static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true))
    {
       GET_G00_SAFE(matGRe, matBRe, matD, matND);
    }
    else 
    {
       GET_G00(matGRe, matBRe, matD, matND);
    }
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 1);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 2);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 3);

    // Column 1
    if constexpr(static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true))
    {
       GET_G11_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    else 
    {
       GET_G11(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 2, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);

    // Column 2
    if constexpr(static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true))
    {
       GET_G22_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
       GET_G22(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);

    // Column 3
    if constexpr(static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true))
    {
       GET_G33_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
       GET_G33(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }

    SET_Lii(matLRe, matLIm, matD, 0);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 1, 0);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 2, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 0, temp0, temp1);

    // Column 1
    SET_Lii(matLRe, matLIm, matD, 1);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 2, 1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 1, temp0, temp1);

    // Column 2
    SET_Lii(matLRe, matLIm, matD, 2);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 3, 2);

    // Column 3
    SET_Lii(matLRe, matLIm, matD, 3);

    for(i = 0; i < MAT_SQR_SIZE_4; ++i)
    {
        matLReii        = matLRe[i][i];
        matInvBRe[i][i] = _mm512_mul_ps(matLReii,matLReii);
        for(k = (i+1); k < MAT_SQR_SIZE_4; ++k) 
        {
            matLReki = matLRe[k][i];
            matLImki = matLIm[k][i];
            temp1    = _mm512_fmadd_ps(matLReki,matLReki,_mm512_mul_ps(matLImki,matLImki));
            matInvBRe[i][i] = _mm512_add_ps(matInvBRe[i][i],temp1);
        } 
        matInvBIm[i][i] = vzero;
    }

    // Unrolled out loop (the i-loop)
#if 1
    for(j = 1;j < MAT_SQR_SIZE_4; ++j)   
    {
        matLRejj        = matLRe[j][j];
        matInvBRe[0][j] = _mm512_mul_ps(matLRe[j][0],matLRejj);
        matInvBIm[0][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][0],matLRejj));
        for(k = (j+1); k < MAT_SQR_SIZE_4; k++)
        {
            GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][0],matLIm[k][0],temp1,temp2);
            matInvBRe[0][j] = _mm512_add_ps(matInvBRe[0][j], temp1);
            matInvBIm[0][j] = _mm512_add_ps(matInvBIm[0][j], temp2);
        }
        matInvBRe[j][0] = matInvBRe[0][j];
        matInvBIm[j][0] = _mm512_sub_ps(vzero, matInvBIm[0][j]);
    }

    for(j = 2;j < MAT_SQR_SIZE_4; ++j)   
    {
        matLRejj        = matLRe[j][j];
        matInvBRe[1][j] = _mm512_mul_ps(matLRe[j][1],matLRejj);
        matInvBIm[1][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][1],matLRejj));
        for(k = (j+1); k < MAT_SQR_SIZE_4; k++)
        {
            GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][1],matLIm[k][1],temp1,temp2);
            matInvBRe[1][j] = _mm512_add_ps(matInvBRe[1][j], temp1);
            matInvBIm[1][j] = _mm512_add_ps(matInvBIm[1][j], temp2);
        }
        matInvBRe[j][1] = matInvBRe[1][j];
        matInvBIm[j][1] = _mm512_sub_ps(vzero, matInvBIm[1][j]);
    }

    for(j = 3;j < MAT_SQR_SIZE_4; ++j)   
    {
        matLRejj        = matLRe[j][j];
        matInvBRe[2][j] = _mm512_mul_ps(matLRe[j][2],matLRejj);
        matInvBIm[2][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][2],matLRejj));
        for(k = (j+1); k < MAT_SQR_SIZE_4; k++)
        {
            GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][2],matLIm[k][2],temp1,temp2);
            matInvBRe[2][j] = _mm512_add_ps(matInvBRe[2][j], temp1);
            matInvBIm[2][j] = _mm512_add_ps(matInvBIm[2][j], temp2);
        }
        matInvBRe[j][2] = matInvBRe[0][j];
        matInvBIm[j][2] = _mm512_sub_ps(vzero, matInvBIm[2][j]);
    }

    for(j = 4;j < MAT_SQR_SIZE_4; ++j)   
    {
        matLRejj        = matLRe[j][j];
        matInvBRe[3][j] = _mm512_mul_ps(matLRe[j][3],matLRejj);
        matInvBIm[3][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][3],matLRejj));
        for(k = (j+1); k < MAT_SQR_SIZE_4; k++)
        {
            GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][3],matLIm[k][3],temp1,temp2);
            matInvBRe[3][j] = _mm512_add_ps(matInvBRe[3][j], temp1);
            matInvBIm[3][j] = _mm512_add_ps(matInvBIm[3][j], temp2);
        }
        matInvBRe[j][3] = matInvBRe[3][j];
        matInvBIm[j][3] = _mm512_sub_ps(vzero, matInvBIm[3][j]);
    }
#endif 

#if 0
    for(i = 0; i < MAT_SQR_SIZE_4; ++i) 
    {
        for(j = (i+1);j < MAT_SQR_SIZE_4; ++j)   
        {
            matLRejj        = matLRe[j][j];
            matInvBRe[i][j] = _mm512_mul_ps(matLRe[j][i],matLRejj);
            matInvBIm[i][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][i],matLRejj));
            for(k = (j+1); k < MAT_SQR_SIZE_4; k++)
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][i],matLIm[k][i],temp1,temp2);
                matInvBRe[i][j] = _mm512_add_ps(matInvBRe[i][j], temp1);
                matInvBIm[i][j] = _mm512_add_ps(matInvBIm[i][j], temp2);
            }
            matInvBRe[j][i] = matInvBRe[i][j];
            matInvBIm[j][i] = _mm512_sub_ps(vzero, matInvBIm[i][j]);
        }
    }
#endif 
}

template void 
gms::math 
::mat_inv_cholesky_4x4_16xf32<true,true>(const __m512 matBRe[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                                         const __m512 matBIm[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                                         __m512       matInvBRe[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                                         __m512       matInvBIm[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4]);

template void 
gms::math 
::mat_inv_cholesky_4x4_16xf32<false,true>(const __m512 matBRe[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                                         const __m512 matBIm[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                                         __m512       matInvBRe[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                                         __m512       matInvBIm[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4]);

template void 
gms::math 
::mat_inv_cholesky_4x4_16xf32<true,false>(const __m512 matBRe[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                                         const __m512 matBIm[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                                         __m512       matInvBRe[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                                         __m512       matInvBIm[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4]);

template void 
gms::math 
::mat_inv_cholesky_4x4_16xf32<false,false>(const __m512 matBRe[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                                           const __m512 matBIm[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                                           __m512       matInvBRe[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                                           __m512       matInvBIm[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4]);

template<bool use_prefetching,bool mitigate_nan>
void 
gms::math
::mat_inv_cholesky_5x5_16xf32(const __m512 matBRe[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                              const __m512 matBIm[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                              __m512       matInvBRe[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                              __m512       matInvBIm[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5])
{
    __ATTR_ALIGN__(64) 
    __m512 matGRe[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5];
    __ATTR_ALIGN__(64)
    __m512 matGIm[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5];
    __ATTR_ALIGN__(64)
    __m512 matLRe[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5];
    __ATTR_ALIGN__(64)
    __m512 matLIm[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5];
    __ATTR_ALIGN__(64)
    __m512 matD[MAT_SQR_SIZE_5];
    __ATTR_ALIGN__(64)
    __m512 matND[MAT_SQR_SIZE_5];
    __m512 temp0, temp1, temp2;
    __m512 matLReii;
    __m512 matLReki;
    __m512 matLImki;
    __m512 matLRejj;
    std::int32_t i,j,k;
    const __m512 vzero{_mm512_setzero_ps()};
    const __m512 vneg_one{_mm512_set1_ps(-1.0f)};
    //////////////////////////////////////////////////////////////
    if constexpr (static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
    {
        constexpr std::int32_t VEC_PS_LEN = 16;
        const float * __restrict ptr_matBRe{reinterpret_cast<const float * __restrict>(&matBRe)};
        const float * __restrict ptr_matBIm{reinterpret_cast<const float * __restrict>(&matBIm)};
        for(i = 0;i != MAT_SQR_SIZE_5; ++i) 
        {
            register std::int32_t outer_idx = i*MAT_SQR_SIZE_5;
            for(j = 0;j != MAT_SQR_SIZE_5; ++j)   
            {
                register std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                _mm_prefetch((const char*)&ptr_matBRe[inner_idx],_MM_HINT_T0);
                _mm_prefetch((const char*)&ptr_matBIm[inner_idx],_MM_HINT_T0);
            }
        }
    }
     /////////////////////////////////// get G, B = G*G', G is a lower triangular matrix
    // Column 0
    if constexpr(static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true))
    {
        GET_G00_SAFE(matGRe, matBRe, matD, matND);
    }
    else 
    {
        GET_G00(matGRe, matBRe, matD, matND);
    }
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 1);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 2);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 3);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 4);

    // Column 1
    if constexpr(static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true))
    {
        GET_G11_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    else 
    {
        GET_G11(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 2, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);

    // Column 2
    if constexpr(static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true))
    {
        GET_G22_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G22(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);

    // Column 3
    if constexpr(static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true))
    {
        GET_G33_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G33(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);

    // Column 4
    if constexpr(static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true))
    {
        GET_G44_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G44(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }

    // Column 0
    SET_Lii(matLRe, matLIm, matD, 0);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 1, 0);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 2, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 0, temp0, temp1);

    // Column 1
    SET_Lii(matLRe, matLIm, matD, 1);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 2, 1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 1, temp0, temp1);

    // Column 2
    SET_Lii(matLRe, matLIm, matD, 2);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 3, 2);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 2, temp0, temp1);

    // Column 3
    SET_Lii(matLRe, matLIm, matD, 3);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 4, 3);

    // Column 4
    SET_Lii(matLRe, matLIm, matD, 4);

    /////////////////////////////////// get invB = L'*L
    for(i = 0; i < MAT_SQR_SIZE_5; ++i)
    {
        matLReii        = matLRe[i][i];
        matInvBRe[i][i] = _mm512_mul_ps(matLReii,matLReii);
        for(k = (i+1); k < MAT_SQR_SIZE_5; ++k) 
        {
            matLReki = matLRe[k][i];
            matLImki = matLIm[k][i];
            temp1    = _mm512_fmadd_ps(matLReki,matLReki,_mm512_mul_ps(matLImki,matLImki));
            matInvBRe[i][i] = _mm512_add_ps(matInvBRe[i][i],temp1);
        } 
        matInvBIm[i][i] = vzero;
    }

#if 1
    for(j = 1;j < MAT_SQR_SIZE_5; ++j)   
    {
            matLRejj        = matLRe[j][j];
            matInvBRe[0][j] = _mm512_mul_ps(matLRe[j][0],matLRejj);
            matInvBIm[0][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][0],matLRejj));
            for(k = (j+1); k < MAT_SQR_SIZE_5; k++)
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][0],matLIm[k][0],temp1,temp2);
                matInvBRe[0][j] = _mm512_add_ps(matInvBRe[0][j], temp1);
                matInvBIm[0][j] = _mm512_add_ps(matInvBIm[0][j], temp2);
            }
            matInvBRe[j][0] = matInvBRe[0][j];
            matInvBIm[j][0] = _mm512_sub_ps(vzero, matInvBIm[0][j]);
    }
    
    for(j = 2;j < MAT_SQR_SIZE_5; ++j)   
    {
            matLRejj        = matLRe[j][j];
            matInvBRe[1][j] = _mm512_mul_ps(matLRe[j][1],matLRejj);
            matInvBIm[1][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][1],matLRejj));
            for(k = (j+1); k < MAT_SQR_SIZE_5; k++)
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][1],matLIm[k][1],temp1,temp2);
                matInvBRe[1][j] = _mm512_add_ps(matInvBRe[1][j], temp1);
                matInvBIm[1][j] = _mm512_add_ps(matInvBIm[1][j], temp2);
            }
            matInvBRe[j][1] = matInvBRe[1][j];
            matInvBIm[j][1] = _mm512_sub_ps(vzero, matInvBIm[1][j]);
    }

    for(j = 3;j < MAT_SQR_SIZE_5; ++j)   
    {
            matLRejj        = matLRe[j][j];
            matInvBRe[2][j] = _mm512_mul_ps(matLRe[j][2],matLRejj);
            matInvBIm[2][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][2],matLRejj));
            for(k = (j+1); k < MAT_SQR_SIZE_5; k++)
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][2],matLIm[k][2],temp1,temp2);
                matInvBRe[2][j] = _mm512_add_ps(matInvBRe[2][j], temp1);
                matInvBIm[2][j] = _mm512_add_ps(matInvBIm[2][j], temp2);
            }
            matInvBRe[j][2] = matInvBRe[2][j];
            matInvBIm[j][2] = _mm512_sub_ps(vzero, matInvBIm[2][j]);
    }

    for(j = 4;j < MAT_SQR_SIZE_5; ++j)   
    {
            matLRejj        = matLRe[j][j];
            matInvBRe[3][j] = _mm512_mul_ps(matLRe[j][3],matLRejj);
            matInvBIm[3][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][3],matLRejj));
            for(k = (j+1); k < MAT_SQR_SIZE_5; k++)
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][3],matLIm[k][3],temp1,temp2);
                matInvBRe[3][j] = _mm512_add_ps(matInvBRe[3][j], temp1);
                matInvBIm[3][j] = _mm512_add_ps(matInvBIm[3][j], temp2);
            }
            matInvBRe[j][3] = matInvBRe[3][j];
            matInvBIm[j][3] = _mm512_sub_ps(vzero, matInvBIm[3][j]);
    }

    for(j = 5;j < MAT_SQR_SIZE_5; ++j)   
    {
            matLRejj        = matLRe[j][j];
            matInvBRe[4][j] = _mm512_mul_ps(matLRe[j][4],matLRejj);
            matInvBIm[4][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][4],matLRejj));
            for(k = (j+1); k < MAT_SQR_SIZE_5; k++)
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][4],matLIm[k][4],temp1,temp2);
                matInvBRe[4][j] = _mm512_add_ps(matInvBRe[4][j], temp1);
                matInvBIm[4][j] = _mm512_add_ps(matInvBIm[4][j], temp2);
            }
            matInvBRe[j][4] = matInvBRe[4][j];
            matInvBIm[j][4] = _mm512_sub_ps(vzero, matInvBIm[4][j]);
    }
    
#endif 
#if 0
    for(i = 0; i < MAT_SQR_SIZE_5; ++i) 
    {
        for(j = (i+1);j < MAT_SQR_SIZE_5; ++j)   
        {
            matLRejj        = matLRe[j][j];
            matInvBRe[i][j] = _mm512_mul_ps(matLRe[j][i],matLRejj);
            matInvBIm[i][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][i],matLRejj));
            for(k = (j+1); k < MAT_SQR_SIZE_5; k++)
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][i],matLIm[k][i],temp1,temp2);
                matInvBRe[i][j] = _mm512_add_ps(matInvBRe[i][j], temp1);
                matInvBIm[i][j] = _mm512_add_ps(matInvBIm[i][j], temp2);
            }
            matInvBRe[j][i] = matInvBRe[i][j];
            matInvBIm[j][i] = _mm512_sub_ps(vzero, matInvBIm[i][j]);
        }
    }
#endif 

}

template void 
gms::math 
::mat_inv_cholesky_5x5_16xf32<true,true>(const __m512 matBRe[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                                         const __m512 matBIm[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                                         __m512       matInvBRe[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                                         __m512       matInvBIm[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5]);

template void 
gms::math 
::mat_inv_cholesky_5x5_16xf32<false,true>(const __m512 matBRe[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                                         const __m512 matBIm[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                                         __m512       matInvBRe[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                                         __m512       matInvBIm[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5]);

template void 
gms::math 
::mat_inv_cholesky_5x5_16xf32<true,false>(const __m512 matBRe[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                                         const __m512 matBIm[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                                         __m512       matInvBRe[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                                         __m512       matInvBIm[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5]);

template void 
gms::math 
::mat_inv_cholesky_5x5_16xf32<false,false>(const __m512 matBRe[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                                           const __m512 matBIm[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                                           __m512       matInvBRe[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                                           __m512       matInvBIm[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5]);


template<bool use_prefetching,bool mitigate_nan>
void 
gms::math
::mat_inv_cholesky_6x6_16xf32(const __m512 matBRe[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                              const __m512 matBIm[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                              __m512       matInvBRe[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                              __m512       matInvBIm[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6])
{
    __ATTR_ALIGN__(64) 
    __m512 matGRe[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6];
    __ATTR_ALIGN__(64)
    __m512 matGIm[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6];
    __ATTR_ALIGN__(64)
    __m512 matLRe[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6];
    __ATTR_ALIGN__(64)
    __m512 matLIm[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6];
    __ATTR_ALIGN__(64)
    __m512 matD[MAT_SQR_SIZE_6];
    __ATTR_ALIGN__(64)
    __m512 matND[MAT_SQR_SIZE_6];
    __m512 temp0, temp1, temp2;
    __m512 matLReii;
    __m512 matLReki;
    __m512 matLImki;
    __m512 matLRejj;
    std::int32_t i,j,k;
    const __m512 vzero{_mm512_setzero_ps()};
    const __m512 vneg_one{_mm512_set1_ps(-1.0f)};
    //////////////////////////////////////////////////////////////
    if constexpr (static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
    {
        constexpr std::int32_t VEC_PS_LEN = 16;
        const float * __restrict ptr_matBRe{reinterpret_cast<const float * __restrict>(&matBRe)};
        const float * __restrict ptr_matBIm{reinterpret_cast<const float * __restrict>(&matBIm)};
        for(i = 0;i != MAT_SQR_SIZE_6; ++i) 
        {
            register std::int32_t outer_idx = i*MAT_SQR_SIZE_6;
            for(j = 0;j != MAT_SQR_SIZE_6; ++j)   
            {
                register std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                _mm_prefetch((const char*)&ptr_matBRe[inner_idx],_MM_HINT_T0);
                _mm_prefetch((const char*)&ptr_matBIm[inner_idx],_MM_HINT_T0);
            }
        }
    }
    const bool cond_mitigate_nan = static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true);
     // Column 0
    if constexpr(cond_mitigate_nan) 
    {
        GET_G00_SAFE(matGRe, matBRe, matD, matND);
    }
    else 
    {
        GET_G00(matGRe, matBRe, matD, matND);
    }
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 1);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 2);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 3);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 4);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 5);

    // Column 1
    if constexpr(cond_mitigate_nan)
    {
        GET_G11_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    else 
    {
        GET_G11(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 2, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);

    // Column 2
    if constexpr(cond_mitigate_nan)
    {
        GET_G22_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G22(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);

    // Column 3
    if constexpr(cond_mitigate_nan)
    {
        GET_G33_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G33(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)
    }
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);

    // Column 4
    if constexpr(cond_mitigate_nan)
    {
        GET_G44_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G44(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 5, 4, temp0, temp1);

    // Column 5
    if constexpr(cond_mitigate_nan) 
    {
        GET_G55_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G55(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }

     /////////////////////////////////// get L, L=invG
    // Column 0
    SET_Lii(matLRe, matLIm, matD, 0);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 1, 0);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 2, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 0, temp0, temp1);

    // Column 1
    SET_Lii(matLRe, matLIm, matD, 1);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 2, 1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 1, temp0, temp1);

    // Column 2
    SET_Lii(matLRe, matLIm, matD, 2);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 3, 2);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 2, temp0, temp1);

    // Column 3
    SET_Lii(matLRe, matLIm, matD, 3);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 4, 3);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 3, temp0, temp1);

    // Column 4
    SET_Lii(matLRe, matLIm, matD, 4);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 5, 4);

    // Column 5
    SET_Lii(matLRe, matLIm, matD, 5);

     /////////////////////////////////// get invB = L'*L
    for(i = 0; i < MAT_SQR_SIZE_6; ++i)
    {
        matLReii        = matLRe[i][i];
        matInvBRe[i][i] = _mm512_mul_ps(matLReii,matLReii);
        for(k = (i+1); k < MAT_SQR_SIZE_6; ++k) 
        {
            matLReki = matLRe[k][i];
            matLImki = matLIm[k][i];
            temp1    = _mm512_fmadd_ps(matLReki,matLReki,_mm512_mul_ps(matLImki,matLImki));
            matInvBRe[i][i] = _mm512_add_ps(matInvBRe[i][i],temp1);
        } 
        matInvBIm[i][i] = vzero;
    }

    for(i = 0; i < MAT_SQR_SIZE_6; ++i) 
    {
        for(j = (i+1);j < MAT_SQR_SIZE_6; ++j)   
        {
            matLRejj        = matLRe[j][j];
            matInvBRe[i][j] = _mm512_mul_ps(matLRe[j][i],matLRejj);
            matInvBIm[i][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][i],matLRejj));
            for(k = (j+1); k < MAT_SQR_SIZE_6; k++)
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][i],matLIm[k][i],temp1,temp2);
                matInvBRe[i][j] = _mm512_add_ps(matInvBRe[i][j], temp1);
                matInvBIm[i][j] = _mm512_add_ps(matInvBIm[i][j], temp2);
            }
            matInvBRe[j][i] = matInvBRe[i][j];
            matInvBIm[j][i] = _mm512_sub_ps(vzero, matInvBIm[i][j]);
        }
    }
}

template void 
gms::math 
::mat_inv_cholesky_6x6_16xf32<true,true>(const __m512 matBRe[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                                         const __m512 matBIm[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                                         __m512       matInvBRe[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                                         __m512       matInvBIm[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6]);

template void 
gms::math 
::mat_inv_cholesky_6x6_16xf32<false,true>(const __m512 matBRe[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                                         const __m512 matBIm[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                                         __m512       matInvBRe[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                                         __m512       matInvBIm[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6]);

template void 
gms::math 
::mat_inv_cholesky_6x6_16xf32<true,false>(const __m512 matBRe[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                                         const __m512 matBIm[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                                         __m512       matInvBRe[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                                         __m512       matInvBIm[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6]);

template void 
gms::math 
::mat_inv_cholesky_6x6_16xf32<false,false>(const __m512 matBRe[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                                           const __m512 matBIm[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                                           __m512       matInvBRe[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                                           __m512       matInvBIm[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6]);

template<bool use_prefetching,bool mitigate_nan>
void 
gms::math
::mat_inv_cholesky_7x7_16xf32(const __m512 matBRe[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                              const __m512 matBIm[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                              __m512       matInvBRe[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                              __m512       matInvBIm[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7])
{
    __ATTR_ALIGN__(64) 
    __m512 matGRe[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7];
    __ATTR_ALIGN__(64)
    __m512 matGIm[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7];
    __ATTR_ALIGN__(64)
    __m512 matLRe[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7];
    __ATTR_ALIGN__(64)
    __m512 matLIm[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7];
    __ATTR_ALIGN__(64)
    __m512 matD[MAT_SQR_SIZE_7];
    __ATTR_ALIGN__(64)
    __m512 matND[MAT_SQR_SIZE_7];
    __m512 temp0, temp1, temp2;
    __m512 matLReii;
    __m512 matLReki;
    __m512 matLImki;
    __m512 matLRejj;
    std::int32_t i,j,k;
    const __m512 vzero{_mm512_setzero_ps()};
    const __m512 vneg_one{_mm512_set1_ps(-1.0f)};
    //////////////////////////////////////////////////////////////
    if constexpr (static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
    {
        constexpr std::int32_t VEC_PS_LEN = 16;
        const float * __restrict ptr_matBRe{reinterpret_cast<const float * __restrict>(&matBRe)};
        const float * __restrict ptr_matBIm{reinterpret_cast<const float * __restrict>(&matBIm)};
        for(i = 0;i != MAT_SQR_SIZE_7; ++i) 
        {
            register std::int32_t outer_idx = i*MAT_SQR_SIZE_7;
            for(j = 0;j != MAT_SQR_SIZE_7; ++j)   
            {
                register std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                _mm_prefetch((const char*)&ptr_matBRe[inner_idx],_MM_HINT_T0);
                _mm_prefetch((const char*)&ptr_matBIm[inner_idx],_MM_HINT_T0);
            }
        }
    }
    const bool cond_mitigate_nan = static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true);
    /////////////////////////////////// get G, B = G*G', G is a lower triangular matrix
    // Column 0
    GET_G00(matGRe, matBRe, matD, matND);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 1);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 2);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 3);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 4);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 5);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 6);

    // Column 1
    if constexpr(cond_mitigate_nan)
    {
        GET_G11_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    else 
    {
        GET_G11(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 2, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);

    // Column 2
    if constexpr(cond_mitigate_nan)
    {
        GET_G22_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G22(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);

    // Column 3
    if constexpr(cond_mitigate_nan)
    {
        GET_G33_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G33(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);

    // Column 4
    if constexpr(cond_mitigate_nan)
    {
        GET_G44_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G44(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 5, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 6, 4, temp0, temp1);

    // Column 5
    if constexpr(cond_mitigate_nan)
    {
        GET_G55_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G55(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 6, 5, temp0, temp1);

    // Column 6
    if constexpr(cond_mitigate_nan)
    {
        GET_G66_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G66(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }

     // Column 0
    SET_Lii(matLRe, matLIm, matD, 0);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 1, 0);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 2, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 0, temp0, temp1);

    // Column 1
    SET_Lii(matLRe, matLIm, matD, 1);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 2, 1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 1, temp0, temp1);

    // Column 2
    SET_Lii(matLRe, matLIm, matD, 2);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 3, 2);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 2, temp0, temp1);

    // Column 3
    SET_Lii(matLRe, matLIm, matD, 3);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 4, 3);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 3, temp0, temp1);

    // Column 4
    SET_Lii(matLRe, matLIm, matD, 4);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 5, 4);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 4, temp0, temp1);

    // Column 5
    SET_Lii(matLRe, matLIm, matD, 5);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 6, 5);

    // Column 6
    SET_Lii(matLRe, matLIm, matD, 6);
    
     /////////////////////////////////// get invB = L'*L
    for(i = 0; i < MAT_SQR_SIZE_7; ++i)
    {
        matLReii        = matLRe[i][i];
        matInvBRe[i][i] = _mm512_mul_ps(matLReii,matLReii);
        for(k = (i+1); k < MAT_SQR_SIZE_7; ++k) 
        {
            matLReki = matLRe[k][i];
            matLImki = matLIm[k][i];
            temp1    = _mm512_fmadd_ps(matLReki,matLReki,_mm512_mul_ps(matLImki,matLImki));
            matInvBRe[i][i] = _mm512_add_ps(matInvBRe[i][i],temp1);
        } 
        matInvBIm[i][i] = vzero;
    }

    for(i = 0; i < MAT_SQR_SIZE_7; ++i) 
    {
        for(j = (i+1);j < MAT_SQR_SIZE_7; ++j)   
        {
            matLRejj        = matLRe[j][j];
            matInvBRe[i][j] = _mm512_mul_ps(matLRe[j][i],matLRejj);
            matInvBIm[i][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][i],matLRejj));
            for(k = (j+1); k < MAT_SQR_SIZE_7; k++)
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][i],matLIm[k][i],temp1,temp2);
                matInvBRe[i][j] = _mm512_add_ps(matInvBRe[i][j], temp1);
                matInvBIm[i][j] = _mm512_add_ps(matInvBIm[i][j], temp2);
            }
            matInvBRe[j][i] = matInvBRe[i][j];
            matInvBIm[j][i] = _mm512_sub_ps(vzero, matInvBIm[i][j]);
        }
    }

}

template void 
gms::math 
::mat_inv_cholesky_7x7_16xf32<true,true>(const __m512 matBRe[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                                         const __m512 matBIm[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                                         __m512       matInvBRe[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                                         __m512       matInvBIm[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7]);

template void 
gms::math 
::mat_inv_cholesky_7x7_16xf32<false,true>(const __m512 matBRe[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                                         const __m512 matBIm[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                                         __m512       matInvBRe[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                                         __m512       matInvBIm[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7]);

template void 
gms::math 
::mat_inv_cholesky_7x7_16xf32<true,false>(const __m512 matBRe[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                                         const __m512 matBIm[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                                         __m512       matInvBRe[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                                         __m512       matInvBIm[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7]);

template void 
gms::math 
::mat_inv_cholesky_7x7_16xf32<false,false>(const __m512 matBRe[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                                           const __m512 matBIm[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                                           __m512       matInvBRe[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                                           __m512       matInvBIm[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7]);

template<bool use_prefetching,bool mitigate_nan>
void 
gms::math
::mat_inv_cholesky_8x8_16xf32(const __m512 matBRe[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                              const __m512 matBIm[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                              __m512       matInvBRe[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                              __m512       matInvBIm[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8])
{
    __ATTR_ALIGN__(64) 
    __m512 matGRe[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8];
    __ATTR_ALIGN__(64)
    __m512 matGIm[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8];
    __ATTR_ALIGN__(64)
    __m512 matLRe[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8];
    __ATTR_ALIGN__(64)
    __m512 matLIm[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8];
    __ATTR_ALIGN__(64)
    __m512 matD[MAT_SQR_SIZE_8];
    __ATTR_ALIGN__(64)
    __m512 matND[MAT_SQR_SIZE_8];
    __m512 temp0, temp1, temp2;
    __m512 matLReii;
    __m512 matLReki;
    __m512 matLImki;
    __m512 matLRejj;
    std::int32_t i,j,k;
    const __m512 vzero{_mm512_setzero_ps()};
    const __m512 vneg_one{_mm512_set1_ps(-1.0f)};
    /////////////////////////////////////////////////////////////////////////////////////////////
    if constexpr (static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
    {
        constexpr std::int32_t VEC_PS_LEN = 16;
        const float * __restrict ptr_matBRe{reinterpret_cast<const float * __restrict>(&matBRe)};
        const float * __restrict ptr_matBIm{reinterpret_cast<const float * __restrict>(&matBIm)};
        for(i = 0;i != MAT_SQR_SIZE_8; ++i) 
        {
            register std::int32_t outer_idx = i*MAT_SQR_SIZE_8;
            for(j = 0;j != MAT_SQR_SIZE_8; ++j)   
            {
                register std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                _mm_prefetch((const char*)&ptr_matBRe[inner_idx],_MM_HINT_T0);
                _mm_prefetch((const char*)&ptr_matBIm[inner_idx],_MM_HINT_T0);
            }
        }
    }
    /////////////////////////////////////////////////////////////////////////////
    const bool cond_nan_mitigation = static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true);
     // Column 0
    GET_G00_SAFE(matGRe, matBRe, matD, matND);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 1);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 2);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 3);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 4);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 5);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 6);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 7);

    // Column 1
    if constexpr(cond_nan_mitigation)
    {
        GET_G11_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    else 
    {
        GET_G11(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 2, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);

    // Column 2
    if constexpr(cond_nan_mitigation)
    {
        GET_G22_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G22(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);

    // Column 3
    if constexpr(cond_nan_mitigation)
    {
        GET_G33_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G33(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);

    // Column 4
    if constexpr(cond_nan_mitigation)
    {
        GET_G44_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G44(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }   
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 5, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 6, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 7, 4, temp0, temp1);

    // Column 5
    if constexpr(cond_nan_mitigation)
    {
        GET_G55_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G55(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }   
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 6, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 7, 5, temp0, temp1);

    // Column 6
    if constexpr(cond_nan_mitigation)
    {
        GET_G66_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else
    {
        GET_G66(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 7, 6, temp0, temp1);

    // Column 7
    if constexpr(cond_nan_mitigation)
    {
        GET_G77_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G77(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }

    /////////////////////////////////// get L, L=invG
    // Column 0
    SET_Lii(matLRe, matLIm, matD, 0);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 1, 0);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 2, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 0, temp0, temp1);

    // Column 1
    SET_Lii(matLRe, matLIm, matD, 1);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 2, 1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 1, temp0, temp1);

    // Column 2
    SET_Lii(matLRe, matLIm, matD, 2);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 3, 2);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 2, temp0, temp1);

    // Column 3
    SET_Lii(matLRe, matLIm, matD, 3);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 4, 3);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 3, temp0, temp1);

    // Column 4
    SET_Lii(matLRe, matLIm, matD, 4);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 5, 4);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 4, temp0, temp1);

    // Column 5
    SET_Lii(matLRe, matLIm, matD, 5);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 6, 5);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 5, temp0, temp1);

    // Column 6
    SET_Lii(matLRe, matLIm, matD, 6);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 7, 6);

    // Column 7
    SET_Lii(matLRe, matLIm, matD, 7);

     /////////////////////////////////// get invB = L'*L
    for(i = 0; i < MAT_SQR_SIZE_8; ++i)
    {
        matLReii        = matLRe[i][i];
        matInvBRe[i][i] = _mm512_mul_ps(matLReii,matLReii);
        for(k = (i+1); k < MAT_SQR_SIZE_8; ++k) 
        {
            matLReki = matLRe[k][i];
            matLImki = matLIm[k][i];
            temp1    = _mm512_fmadd_ps(matLReki,matLReki,_mm512_mul_ps(matLImki,matLImki));
            matInvBRe[i][i] = _mm512_add_ps(matInvBRe[i][i],temp1);
        } 
        matInvBIm[i][i] = vzero;
    }

    for(i = 0; i < MAT_SQR_SIZE_8; ++i) 
    {
        for(j = (i+1);j < MAT_SQR_SIZE_8; ++j)   
        {
            matLRejj        = matLRe[j][j];
            matInvBRe[i][j] = _mm512_mul_ps(matLRe[j][i],matLRejj);
            matInvBIm[i][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][i],matLRejj));
            for(k = (j+1); k < MAT_SQR_SIZE_8; k++)
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][i],matLIm[k][i],temp1,temp2);
                matInvBRe[i][j] = _mm512_add_ps(matInvBRe[i][j], temp1);
                matInvBIm[i][j] = _mm512_add_ps(matInvBIm[i][j], temp2);
            }
            matInvBRe[j][i] = matInvBRe[i][j];
            matInvBIm[j][i] = _mm512_sub_ps(vzero, matInvBIm[i][j]);
        }
    }
}

template void 
gms::math 
::mat_inv_cholesky_8x8_16xf32<true,true>(const __m512 matBRe[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                                         const __m512 matBIm[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                                         __m512       matInvBRe[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                                         __m512       matInvBIm[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8]);

template void 
gms::math 
::mat_inv_cholesky_8x8_16xf32<false,true>(const __m512 matBRe[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                                         const __m512 matBIm[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                                         __m512       matInvBRe[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                                         __m512       matInvBIm[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8]);

template void 
gms::math 
::mat_inv_cholesky_8x8_16xf32<true,false>(const __m512 matBRe[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                                         const __m512 matBIm[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                                         __m512       matInvBRe[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                                         __m512       matInvBIm[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8]);

template void 
gms::math 
::mat_inv_cholesky_8x8_16xf32<false,false>(const __m512 matBRe[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                                           const __m512 matBIm[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                                           __m512       matInvBRe[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                                           __m512       matInvBIm[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8]);

template<bool use_prefetching,bool mitigate_nan>
void 
gms::math
::mat_inv_cholesky_9x9_16xf32(const __m512 matBRe[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                              const __m512 matBIm[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                              __m512       matInvBRe[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                              __m512       matInvBIm[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9])
{
    __ATTR_ALIGN__(64) 
    __m512 matGRe[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9];
    __ATTR_ALIGN__(64)
    __m512 matGIm[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9];
    __ATTR_ALIGN__(64)
    __m512 matLRe[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9];
    __ATTR_ALIGN__(64)
    __m512 matLIm[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9];
    __ATTR_ALIGN__(64)
    __m512 matD[MAT_SQR_SIZE_9];
    __ATTR_ALIGN__(64)
    __m512 matND[MAT_SQR_SIZE_9];
    __m512 temp0, temp1, temp2;
    __m512 matLReii;
    __m512 matLReki;
    __m512 matLImki;
    __m512 matLRejj;
    std::int32_t i,j,k;
    const __m512 vzero{_mm512_setzero_ps()};
    const __m512 vneg_one{_mm512_set1_ps(-1.0f)};
    /////////////////////////////////////////////////////////////////////////////////////////////
    if constexpr (static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
    {
        constexpr std::int32_t VEC_PS_LEN = 16;
        const float * __restrict ptr_matBRe{reinterpret_cast<const float * __restrict>(&matBRe)};
        const float * __restrict ptr_matBIm{reinterpret_cast<const float * __restrict>(&matBIm)};
        for(i = 0;i != MAT_SQR_SIZE_9; ++i) 
        {
            register std::int32_t outer_idx = i*MAT_SQR_SIZE_9;
            for(j = 0;j != MAT_SQR_SIZE_9; ++j)   
            {
                register std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                _mm_prefetch((const char*)&ptr_matBRe[inner_idx],_MM_HINT_T0);
                _mm_prefetch((const char*)&ptr_matBIm[inner_idx],_MM_HINT_T0);
            }
        }
    }
    const bool cond_nan_mitigation = static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true);
    /////////////////////////////////// get G, B = G*G', G is a lower triangular matrix
    // Column 0
    GET_G00(matGRe, matBRe, matD, matND);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 1);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 2);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 3);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 4);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 5);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 6);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 7);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 8);

    // Column 1
    if constexpr(cond_nan_mitigation)
    {
        GET_G11_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    else 
    {
        GET_G11(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 2, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);

    // Column 2
    if constexpr(cond_nan_mitigation)
    {
        GET_G22_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G22(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }    
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);

    // Column 3
    if constexpr(cond_nan_mitigation)
    {
        GET_G33_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else
    {
        GET_G33(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);

    // Column 4
    if constexpr(cond_nan_mitigation)
    {
        GET_G44_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G44(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 5, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 6, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 7, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 8, 4, temp0, temp1);

    // Column 5
    if constexpr(cond_nan_mitigation)
    {
         GET_G55_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
         GET_G55(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 6, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 7, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 8, 5, temp0, temp1);

    // Column 6
    if constexpr(cond_nan_mitigation)
    {
        GET_G66_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G66(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 7, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 8, 6, temp0, temp1);

    // Column 7
    if constexpr(cond_nan_mitigation)
    {
        GET_G77_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else
    {
        GET_G77(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 8, 7, temp0, temp1);

    // Column 8
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 8);
    }
    else 
    {
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 8);
    }

    /////////////////////////////////// get L, L=invG
    // Column 0
    SET_Lii(matLRe, matLIm, matD, 0);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 1, 0);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 2, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 0, temp0, temp1);

    // Column 1
    SET_Lii(matLRe, matLIm, matD, 1);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 2, 1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 1, temp0, temp1);
	GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 1, temp0, temp1);

    // Column 2
    SET_Lii(matLRe, matLIm, matD, 2);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 3, 2);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 2, temp0, temp1);
	GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 2, temp0, temp1);

    // Column 3
    SET_Lii(matLRe, matLIm, matD, 3);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 4, 3);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 3, temp0, temp1);
	GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 3, temp0, temp1);

    // Column 4
    SET_Lii(matLRe, matLIm, matD, 4);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 5, 4);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 4, temp0, temp1);
	GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 4, temp0, temp1);

    // Column 5
    SET_Lii(matLRe, matLIm, matD, 5);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 6, 5);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 5, temp0, temp1);
	GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 5, temp0, temp1);

    // Column 6
    SET_Lii(matLRe, matLIm, matD, 6);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 7, 6);
	GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 6, temp0, temp1);

    // Column 7
    SET_Lii(matLRe, matLIm, matD, 7);
	GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 8, 7);

	//Column 8
	SET_Lii(matLRe, matLIm, matD, 8);

      /////////////////////////////////// get invB = L'*L
    for(i = 0; i < MAT_SQR_SIZE_9; ++i)
    {
        matLReii        = matLRe[i][i];
        matInvBRe[i][i] = _mm512_mul_ps(matLReii,matLReii);
        for(k = (i+1); k < MAT_SQR_SIZE_9; ++k) 
        {
            matLReki = matLRe[k][i];
            matLImki = matLIm[k][i];
            temp1    = _mm512_fmadd_ps(matLReki,matLReki,_mm512_mul_ps(matLImki,matLImki));
            matInvBRe[i][i] = _mm512_add_ps(matInvBRe[i][i],temp1);
        } 
        matInvBIm[i][i] = vzero;
    }

    for(i = 0; i < MAT_SQR_SIZE_9; ++i) 
    {
        for(j = (i+1);j < MAT_SQR_SIZE_9; ++j)   
        {
            matLRejj        = matLRe[j][j];
            matInvBRe[i][j] = _mm512_mul_ps(matLRe[j][i],matLRejj);
            matInvBIm[i][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][i],matLRejj));
            for(k = (j+1); k < MAT_SQR_SIZE_9; k++)
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][i],matLIm[k][i],temp1,temp2);
                matInvBRe[i][j] = _mm512_add_ps(matInvBRe[i][j], temp1);
                matInvBIm[i][j] = _mm512_add_ps(matInvBIm[i][j], temp2);
            }
            matInvBRe[j][i] = matInvBRe[i][j];
            matInvBIm[j][i] = _mm512_sub_ps(vzero, matInvBIm[i][j]);
        }
    }
}

template void 
gms::math 
::mat_inv_cholesky_9x9_16xf32<true,true>(const __m512 matBRe[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                                         const __m512 matBIm[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                                         __m512       matInvBRe[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                                         __m512       matInvBIm[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9]);

template void 
gms::math 
::mat_inv_cholesky_9x9_16xf32<false,true>(const __m512 matBRe[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                                         const __m512 matBIm[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                                         __m512       matInvBRe[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                                         __m512       matInvBIm[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9]);

template void 
gms::math 
::mat_inv_cholesky_9x9_16xf32<true,false>(const __m512 matBRe[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                                         const __m512 matBIm[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                                         __m512       matInvBRe[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                                         __m512       matInvBIm[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9]);

template void 
gms::math 
::mat_inv_cholesky_9x9_16xf32<false,false>(const __m512 matBRe[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                                           const __m512 matBIm[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                                           __m512       matInvBRe[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                                           __m512       matInvBIm[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9]);

template<bool use_prefetching,bool mitigate_nan>
void 
gms::math
::mat_inv_cholesky_10x10_16xf32(const __m512 matBRe[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                                const __m512 matBIm[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                                __m512       matInvBRe[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                                __m512       matInvBIm[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10])
{
    __ATTR_ALIGN__(64) 
    __m512 matGRe[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10];
    __ATTR_ALIGN__(64)
    __m512 matGIm[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10];
    __ATTR_ALIGN__(64)
    __m512 matLRe[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10];
    __ATTR_ALIGN__(64)
    __m512 matLIm[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10];
    __ATTR_ALIGN__(64)
    __m512 matD[MAT_SQR_SIZE_10];
    __ATTR_ALIGN__(64)
    __m512 matND[MAT_SQR_SIZE_10];
    __m512 temp0, temp1, temp2;
    __m512 matLReii;
    __m512 matLReki;
    __m512 matLImki;
    __m512 matLRejj;
    std::int32_t i,j,k;
    const __m512 vzero{_mm512_setzero_ps()};
    const __m512 vneg_one{_mm512_set1_ps(-1.0f)};
    /////////////////////////////////////////////////////////////////////////////////////////////
    if constexpr (static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
    {
        constexpr std::int32_t VEC_PS_LEN = 16;
        const float * __restrict ptr_matBRe{reinterpret_cast<const float * __restrict>(&matBRe)};
        const float * __restrict ptr_matBIm{reinterpret_cast<const float * __restrict>(&matBIm)};
        for(i = 0;i != MAT_SQR_SIZE_10; ++i) 
        {
            register std::int32_t outer_idx = i*MAT_SQR_SIZE_10;
            for(j = 0;j != MAT_SQR_SIZE_10; ++j)   
            {
                register std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                _mm_prefetch((const char*)&ptr_matBRe[inner_idx],_MM_HINT_T0);
                _mm_prefetch((const char*)&ptr_matBIm[inner_idx],_MM_HINT_T0);
            }
        }
    }
    const bool cond_nan_mitigation = static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true);
    GET_G00(matGRe, matBRe, matD, matND);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 1);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 2);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 3);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 4);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 5);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 6);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 7);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 8);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 9);

    // Column 1
    if constexpr(cond_nan_mitigation)
    {
        GET_G11_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    else 
    {
        GET_G11(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 2, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);

    // Column 2
    if constexpr(cond_nan_mitigation)
    {
        GET_G22_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G22(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);

    // Column 3
    if constexpr(cond_nan_mitigation)
    {
        GET_G33_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G33(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);

    // Column 4
    if constexpr(cond_nan_mitigation)
    {
        GET_G44_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G44(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 5, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 6, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 7, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 8, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 4, temp0, temp1);

    // Column 5
    if constexpr(cond_nan_mitigation)
    {
        GET_G55_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G55(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 6, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 7, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 8, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 9, 5, temp0, temp1);

    // Column 6
    if constexpr(cond_nan_mitigation)
    {
        GET_G66_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G66(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 7, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 8, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 6, temp0, temp1);

    // Column 7
    if constexpr(cond_nan_mitigation)
    {
        GET_G77_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G77(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 8, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 9, 7, temp0, temp1);

    // Column 8
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 8);
    }
    else 
    {
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 8);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 8, temp0, temp1);

    //Column 9
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_ODD_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 9);
    }
    else 
    {
        GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 9);
    }

    /////////////////////////////////// get L, L=invG
    // Column 0
    SET_Lii(matLRe, matLIm, matD, 0);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 1, 0);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 2, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 0, temp0, temp1);

    // Column 1
    SET_Lii(matLRe, matLIm, matD, 1);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 2, 1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 1, temp0, temp1);

    // Column 2
    SET_Lii(matLRe, matLIm, matD, 2);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 3, 2);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 2, temp0, temp1);
	GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 2, temp0, temp1);

    // Column 3
    SET_Lii(matLRe, matLIm, matD, 3);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 4, 3);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 3, temp0, temp1);
	GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 3, temp0, temp1);

    // Column 4
    SET_Lii(matLRe, matLIm, matD, 4);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 5, 4);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 4, temp0, temp1);
	GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 4, temp0, temp1);

    // Column 5
    SET_Lii(matLRe, matLIm, matD, 5);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 6, 5);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 5, temp0, temp1);
	GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 5, temp0, temp1);

    // Column 6
    SET_Lii(matLRe, matLIm, matD, 6);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 7, 6);
	GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 6, temp0, temp1);

    // Column 7
    SET_Lii(matLRe, matLIm, matD, 7);
	GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 8, 7);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 7, temp0, temp1);

	//Column 8
	SET_Lii(matLRe, matLIm, matD, 8);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 9, 8);

    //Column 9
    SET_Lii(matLRe, matLIm, matD, 9);

      /////////////////////////////////// get invB = L'*L
    for(i = 0; i < MAT_SQR_SIZE_10; ++i)
    {
        matLReii        = matLRe[i][i];
        matInvBRe[i][i] = _mm512_mul_ps(matLReii,matLReii);
        for(k = (i+1); k < MAT_SQR_SIZE_10; ++k) 
        {
            matLReki = matLRe[k][i];
            matLImki = matLIm[k][i];
            temp1    = _mm512_fmadd_ps(matLReki,matLReki,_mm512_mul_ps(matLImki,matLImki));
            matInvBRe[i][i] = _mm512_add_ps(matInvBRe[i][i],temp1);
        } 
        matInvBIm[i][i] = vzero;
    }

    for(i = 0; i < MAT_SQR_SIZE_10; ++i) 
    {
        for(j = (i+1);j < MAT_SQR_SIZE_10; ++j)   
        {
            matLRejj        = matLRe[j][j];
            matInvBRe[i][j] = _mm512_mul_ps(matLRe[j][i],matLRejj);
            matInvBIm[i][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][i],matLRejj));
            for(k = (j+1); k < MAT_SQR_SIZE_10; k++)
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][i],matLIm[k][i],temp1,temp2);
                matInvBRe[i][j] = _mm512_add_ps(matInvBRe[i][j], temp1);
                matInvBIm[i][j] = _mm512_add_ps(matInvBIm[i][j], temp2);
            }
            matInvBRe[j][i] = matInvBRe[i][j];
            matInvBIm[j][i] = _mm512_sub_ps(vzero, matInvBIm[i][j]);
        }
    }
}

template void 
gms::math 
::mat_inv_cholesky_10x10_16xf32<true,true>(const __m512 matBRe[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                                           const __m512 matBIm[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                                           __m512       matInvBRe[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                                           __m512       matInvBIm[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10]);

template void 
gms::math 
::mat_inv_cholesky_10x10_16xf32<false,true>(const __m512 matBRe[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                                            const __m512 matBIm[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                                            __m512       matInvBRe[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                                            __m512       matInvBIm[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10]);

template void 
gms::math 
::mat_inv_cholesky_10x10_16xf32<true,false>(const __m512 matBRe[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                                            const __m512 matBIm[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                                            __m512       matInvBRe[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                                            __m512       matInvBIm[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10]);

template void 
gms::math 
::mat_inv_cholesky_10x10_16xf32<false,false>(const __m512 matBRe[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                                             const __m512 matBIm[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                                             __m512       matInvBRe[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                                             __m512       matInvBIm[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10]);

template<bool use_prefetching,bool mitigate_nan>
void 
gms::math
::mat_inv_cholesky_11x11_16xf32(const __m512 matBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                                const __m512 matBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                                __m512       matInvBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                                __m512       matInvBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11])
{
    __ATTR_ALIGN__(64) 
    __m512 matGRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
    __ATTR_ALIGN__(64)
    __m512 matGIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
    __ATTR_ALIGN__(64)
    __m512 matLRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
    __ATTR_ALIGN__(64)
    __m512 matLIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
    __ATTR_ALIGN__(64)
    __m512 matD[MAT_SQR_SIZE_11];
    __ATTR_ALIGN__(64)
    __m512 matND[MAT_SQR_SIZE_11];
    __m512 temp0, temp1, temp2;
    __m512 matLReii;
    __m512 matLReki;
    __m512 matLImki;
    __m512 matLRejj;
    std::int32_t i,j,k;
    const __m512 vzero{_mm512_setzero_ps()};
    const __m512 vneg_one{_mm512_set1_ps(-1.0f)};
    /////////////////////////////////////////////////////////////////////////////////////////////
    if constexpr (static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
    {
        constexpr std::int32_t VEC_PS_LEN = 16;
        const float * __restrict ptr_matBRe{reinterpret_cast<const float * __restrict>(&matBRe)};
        const float * __restrict ptr_matBIm{reinterpret_cast<const float * __restrict>(&matBIm)};
        for(i = 0;i != MAT_SQR_SIZE_11; ++i) 
        {
            register std::int32_t outer_idx = i*MAT_SQR_SIZE_11;
            for(j = 0;j != MAT_SQR_SIZE_11; ++j)   
            {
                register std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                _mm_prefetch((const char*)&ptr_matBRe[inner_idx],_MM_HINT_T0);
                _mm_prefetch((const char*)&ptr_matBIm[inner_idx],_MM_HINT_T0);
            }
        }
    }
    const bool cond_nan_mitigation = static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true);
     // Column 0
    GET_G00(matGRe, matBRe, matD, matND);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 1);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 2);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 3);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 4);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 5);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 6);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 7);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 8);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 9);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 10);

    // Column 1
    if constexpr(cond_nan_mitigation)
    {
        GET_G11_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    else 
    {
        GET_G11(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 2, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);

    // Column 2
    if constexpr(cond_nan_mitigation)
    {
        GET_G22_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G22(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);

    // Column 3
    if constexpr(cond_nan_mitigation)
    {
        GET_G33_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G33(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);

    // Column 4
    if constexpr(cond_nan_mitigation)
    {
        GET_G44_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G44(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 5, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 6, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 7, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 8, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 10, 4, temp0, temp1);

    // Column 5
    if constexpr(cond_nan_mitigation)
    {
        GET_G55_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G55(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 6, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 7, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 8, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 9, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 10, 5, temp0, temp1);

    // Column 6
    if constexpr(cond_nan_mitigation)
    {
        GET_G66_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G66(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 7, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 8, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 10, 6, temp0, temp1);

    // Column 7
    if constexpr(cond_nan_mitigation)
    {
        GET_G77_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G77(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 8, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 9, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 10, 7, temp0, temp1);

    // Column 8
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 8);
    }
    else 
    {
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 8);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 10, 8, temp0, temp1);

    // Column 9
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_ODD_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 9);
    }
    else 
    {
        GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 9);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 10, 9, temp0, temp1);

    // Column 10
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 10);
    }
    else 
    {
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 10);
    }

    /////////////////////////////////// get L, L=invG
    // Column 0
    SET_Lii(matLRe, matLIm, matD, 0);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 1, 0);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 2, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 0, temp0, temp1);

    // Column 1
    SET_Lii(matLRe, matLIm, matD, 1);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 2, 1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 1, temp0, temp1);

    // Column 2
    SET_Lii(matLRe, matLIm, matD, 2);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 3, 2);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 2, temp0, temp1);

    // Column 3
    SET_Lii(matLRe, matLIm, matD, 3);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 4, 3);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 3, temp0, temp1);

    // Column 4
    SET_Lii(matLRe, matLIm, matD, 4);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 5, 4);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 4, temp0, temp1);

    // Column 5
    SET_Lii(matLRe, matLIm, matD, 5);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 6, 5);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 5, temp0, temp1);

    // Column 6
    SET_Lii(matLRe, matLIm, matD, 6);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 7, 6);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 6, temp0, temp1);

    // Column 7
    SET_Lii(matLRe, matLIm, matD, 7);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 8, 7);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 7, temp0, temp1);

    // Column 8
    SET_Lii(matLRe, matLIm, matD, 8);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 9, 8);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 8, temp0, temp1);

    // Column 9
    SET_Lii(matLRe, matLIm, matD, 9);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 10, 9);

    // Column 10
    SET_Lii(matLRe, matLIm, matD, 10);

      /////////////////////////////////// get invB = L'*L
    for(i = 0; i < MAT_SQR_SIZE_11; ++i)
    {
        matLReii        = matLRe[i][i];
        matInvBRe[i][i] = _mm512_mul_ps(matLReii,matLReii);
        for(k = (i+1); k < MAT_SQR_SIZE_11; ++k) 
        {
            matLReki = matLRe[k][i];
            matLImki = matLIm[k][i];
            temp1    = _mm512_fmadd_ps(matLReki,matLReki,_mm512_mul_ps(matLImki,matLImki));
            matInvBRe[i][i] = _mm512_add_ps(matInvBRe[i][i],temp1);
        } 
        matInvBIm[i][i] = vzero;
    }

    for(i = 0; i < MAT_SQR_SIZE_11; ++i) 
    {
        for(j = (i+1);j < MAT_SQR_SIZE_11; ++j)   
        {
            matLRejj        = matLRe[j][j];
            matInvBRe[i][j] = _mm512_mul_ps(matLRe[j][i],matLRejj);
            matInvBIm[i][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][i],matLRejj));
            for(k = (j+1); k < MAT_SQR_SIZE_11; k++)
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][i],matLIm[k][i],temp1,temp2);
                matInvBRe[i][j] = _mm512_add_ps(matInvBRe[i][j], temp1);
                matInvBIm[i][j] = _mm512_add_ps(matInvBIm[i][j], temp2);
            }
            matInvBRe[j][i] = matInvBRe[i][j];
            matInvBIm[j][i] = _mm512_sub_ps(vzero, matInvBIm[i][j]);
        }
    }
}

template void 
gms::math 
::mat_inv_cholesky_11x11_16xf32<true,true>(const __m512 matBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                                           const __m512 matBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                                           __m512       matInvBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                                           __m512       matInvBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11]);

template void 
gms::math 
::mat_inv_cholesky_11x11_16xf32<false,true>(const __m512 matBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                                            const __m512 matBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                                            __m512       matInvBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                                            __m512       matInvBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11]);

template void 
gms::math 
::mat_inv_cholesky_11x11_16xf32<true,false>(const __m512 matBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                                            const __m512 matBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                                            __m512       matInvBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                                            __m512       matInvBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11]);

template void 
gms::math 
::mat_inv_cholesky_11x11_16xf32<false,false>(const __m512 matBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                                             const __m512 matBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                                             __m512       matInvBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                                             __m512       matInvBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11]);

template<bool use_prefetching,bool mitigate_nan>
void 
gms::math
::mat_inv_cholesky_12x12_16xf32(const __m512 matBRe[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                                const __m512 matBIm[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                                __m512       matInvBRe[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                                __m512       matInvBIm[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12])
{
    __ATTR_ALIGN__(64) 
    __m512 matGRe[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12];
    __ATTR_ALIGN__(64)
    __m512 matGIm[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12];
    __ATTR_ALIGN__(64)
    __m512 matLRe[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12];
    __ATTR_ALIGN__(64)
    __m512 matLIm[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12];
    __ATTR_ALIGN__(64)
    __m512 matD[MAT_SQR_SIZE_12];
    __ATTR_ALIGN__(64)
    __m512 matND[MAT_SQR_SIZE_12];
    __m512 temp0, temp1, temp2;
    __m512 matLReii;
    __m512 matLReki;
    __m512 matLImki;
    __m512 matLRejj;
    std::int32_t i,j,k;
    const __m512 vzero{_mm512_setzero_ps()};
    const __m512 vneg_one{_mm512_set1_ps(-1.0f)};
    /////////////////////////////////////////////////////////////////////////////////////////////
    if constexpr (static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
    {
        constexpr std::int32_t VEC_PS_LEN = 16;
        const float * __restrict ptr_matBRe{reinterpret_cast<const float * __restrict>(&matBRe)};
        const float * __restrict ptr_matBIm{reinterpret_cast<const float * __restrict>(&matBIm)};
        for(i = 0;i != MAT_SQR_SIZE_12; ++i) 
        {
            register std::int32_t outer_idx = i*MAT_SQR_SIZE_12;
            for(j = 0;j != MAT_SQR_SIZE_12; ++j)   
            {
                register std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                _mm_prefetch((const char*)&ptr_matBRe[inner_idx],_MM_HINT_T0);
                _mm_prefetch((const char*)&ptr_matBIm[inner_idx],_MM_HINT_T0);
            }
        }
    }
    const bool cond_nan_mitigation = static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Column 0
    GET_G00(matGRe, matBRe, matD, matND);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 1);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 2);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 3);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 4);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 5);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 6);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 7);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 8);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 9);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 10);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 11);

    // Column 1
    if constexpr(cond_nan_mitigation)
    {
        GET_G11_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    else 
    {
        GET_G11(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 2, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 11, temp0, temp1);

    // Column 2
    if constexpr(cond_nan_mitigation)
    {
        GET_G22_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G22(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 11, temp0, temp1);

    // Column 3
    if constexpr(cond_nan_mitigation)
    {
        GET_G33_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G33(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 11, temp0, temp1);

    // Column 4
    if constexpr(cond_nan_mitigation)
    {
        GET_G44_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G44(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 5, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 6, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 7, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 8, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 10, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 4, temp0, temp1);

    // Column 5
    if constexpr(cond_nan_mitigation)
    {
        GET_G55_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G55(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 6, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 7, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 8, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 9, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 10, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 11, 5, temp0, temp1);

    // Column 6
    if constexpr(cond_nan_mitigation)
    {
        GET_G66_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G66(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 7, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 8, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 10, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 6, temp0, temp1);

    // Column 7
    if constexpr(cond_nan_mitigation)
    {
        GET_G77_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G77(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 8, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 9, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 10, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 11, 7, temp0, temp1);

    // Column 8
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 8);
    }
    else 
    {
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 8);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 10, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 8, temp0, temp1);

    // Column 9
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_ODD_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 9);
    }
    else 
    {
        GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 9);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 10, 9, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 11, 9, temp0, temp1);

    // Column 10
    if constexpr(cond_nan_mitigation) 
    {
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 10);
    }
    else 
    {
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 10);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 10, temp0, temp1);

    //Column 11
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_ODD_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 11);
    }
    else 
    {
        GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 11);
    }

    /////////////////////////////////// get L, L=invG
    // Column 0
    SET_Lii(matLRe, matLIm, matD, 0);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 1, 0);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 2, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 0, temp0, temp1);

    // Column 1
    SET_Lii(matLRe, matLIm, matD, 1);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 2, 1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 1, temp0, temp1);

    // Column 2
    SET_Lii(matLRe, matLIm, matD, 2);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 3, 2);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 2, temp0, temp1);

    // Column 3
    SET_Lii(matLRe, matLIm, matD, 3);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 4, 3);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 3, temp0, temp1);

    // Column 4
    SET_Lii(matLRe, matLIm, matD, 4);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 5, 4);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 4, temp0, temp1);

    // Column 5
    SET_Lii(matLRe, matLIm, matD, 5);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 6, 5);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 5, temp0, temp1);

    // Column 6
    SET_Lii(matLRe, matLIm, matD, 6);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 7, 6);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 6, temp0, temp1);

    // Column 7
    SET_Lii(matLRe, matLIm, matD, 7);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 8, 7);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 7, temp0, temp1);

    // Column 8
    SET_Lii(matLRe, matLIm, matD, 8);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 9, 8);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 8, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 8, temp0, temp1);

    // Column 9
    SET_Lii(matLRe, matLIm, matD, 9);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 10, 9);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 9, temp0, temp1);

    // Column 10
    SET_Lii(matLRe, matLIm, matD, 10);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 11, 10);

    //Column 11 
    SET_Lii(matLRe, matLIm, matD, 11);

      /////////////////////////////////// get invB = L'*L
    for(i = 0; i < MAT_SQR_SIZE_12; ++i)
    {
        matLReii        = matLRe[i][i];
        matInvBRe[i][i] = _mm512_mul_ps(matLReii,matLReii);
        for(k = (i+1); k < MAT_SQR_SIZE_12; ++k) 
        {
            matLReki = matLRe[k][i];
            matLImki = matLIm[k][i];
            temp1    = _mm512_fmadd_ps(matLReki,matLReki,_mm512_mul_ps(matLImki,matLImki));
            matInvBRe[i][i] = _mm512_add_ps(matInvBRe[i][i],temp1);
        } 
        matInvBIm[i][i] = vzero;
    }

    for(i = 0; i < MAT_SQR_SIZE_12; ++i) 
    {
        for(j = (i+1);j < MAT_SQR_SIZE_12; ++j)   
        {
            matLRejj        = matLRe[j][j];
            matInvBRe[i][j] = _mm512_mul_ps(matLRe[j][i],matLRejj);
            matInvBIm[i][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][i],matLRejj));
            for(k = (j+1); k < MAT_SQR_SIZE_12; k++)
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][i],matLIm[k][i],temp1,temp2);
                matInvBRe[i][j] = _mm512_add_ps(matInvBRe[i][j], temp1);
                matInvBIm[i][j] = _mm512_add_ps(matInvBIm[i][j], temp2);
            }
            matInvBRe[j][i] = matInvBRe[i][j];
            matInvBIm[j][i] = _mm512_sub_ps(vzero, matInvBIm[i][j]);
        }
    }
}

template void 
gms::math 
::mat_inv_cholesky_12x12_16xf32<true,true>(const __m512 matBRe[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                                           const __m512 matBIm[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                                           __m512       matInvBRe[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                                           __m512       matInvBIm[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12]);

template void 
gms::math 
::mat_inv_cholesky_12x12_16xf32<false,true>(const __m512 matBRe[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                                            const __m512 matBIm[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                                            __m512       matInvBRe[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                                            __m512       matInvBIm[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12]);

template void 
gms::math 
::mat_inv_cholesky_12x12_16xf32<true,false>(const __m512 matBRe[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                                            const __m512 matBIm[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                                            __m512       matInvBRe[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                                            __m512       matInvBIm[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12]);

template void 
gms::math 
::mat_inv_cholesky_12x12_16xf32<false,false>(const __m512 matBRe[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                                             const __m512 matBIm[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                                             __m512       matInvBRe[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                                             __m512       matInvBIm[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12]);

template<bool use_prefetching,bool mitigate_nan>
void 
gms::math
::mat_inv_cholesky_13x13_16xf32(const __m512 matBRe[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                                const __m512 matBIm[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                                __m512       matInvBRe[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                                __m512       matInvBIm[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13])
{
    __ATTR_ALIGN__(64) 
    __m512 matGRe[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13];
    __ATTR_ALIGN__(64)
    __m512 matGIm[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13];
    __ATTR_ALIGN__(64)
    __m512 matLRe[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13];
    __ATTR_ALIGN__(64)
    __m512 matLIm[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13];
    __ATTR_ALIGN__(64)
    __m512 matD[MAT_SQR_SIZE_13];
    __ATTR_ALIGN__(64)
    __m512 matND[MAT_SQR_SIZE_13];
    __m512 temp0, temp1, temp2;
    __m512 matLReii;
    __m512 matLReki;
    __m512 matLImki;
    __m512 matLRejj;
    std::int32_t i,j,k;
    const __m512 vzero{_mm512_setzero_ps()};
    const __m512 vneg_one{_mm512_set1_ps(-1.0f)};
    /////////////////////////////////////////////////////////////////////////////////////////////
    if constexpr (static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
    {
        constexpr std::int32_t VEC_PS_LEN = 16;
        const float * __restrict ptr_matBRe{reinterpret_cast<const float * __restrict>(&matBRe)};
        const float * __restrict ptr_matBIm{reinterpret_cast<const float * __restrict>(&matBIm)};
        for(i = 0;i != MAT_SQR_SIZE_13; ++i) 
        {
            register std::int32_t outer_idx = i*MAT_SQR_SIZE_13;
            for(j = 0;j != MAT_SQR_SIZE_13; ++j)   
            {
                register std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                _mm_prefetch((const char*)&ptr_matBRe[inner_idx],_MM_HINT_T0);
                _mm_prefetch((const char*)&ptr_matBIm[inner_idx],_MM_HINT_T0);
            }
        }
    }
    const bool cond_nan_mitigation = static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true);
     // Column 0
    GET_G00(matGRe, matBRe, matD, matND);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 1);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 2);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 3);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 4);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 5);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 6);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 7);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 8);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 9);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 10);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 11);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 12);

    // Column 1
    if constexpr(cond_nan_mitigation)
    {
        GET_G11_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    else 
    {
        GET_G11(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 2, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 11, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 12, temp0, temp1);

    // Column 2
    if constexpr(cond_nan_mitigation)
    {
        GET_G22_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G22(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 11, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 12, temp0, temp1);

    // Column 3
    if constexpr(cond_nan_mitigation)
    {
        GET_G33_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G33(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 11, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 12, temp0, temp1);

    // Column 4
    if constexpr(cond_nan_mitigation)
    {
        GET_G44_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G44(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 5, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 6, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 7, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 8, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 10, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 12, 4, temp0, temp1);

    // Column 5
    if constexpr(cond_nan_mitigation)
    {
        GET_G55_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G55(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 6, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 7, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 8, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 9, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 10, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 11, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 12, 5, temp0, temp1);

    // Column 6
    if constexpr(cond_nan_mitigation)
    {
        GET_G66_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G66(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 7, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 8, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 10, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 12, 6, temp0, temp1);

    // Column 7
    if constexpr(cond_nan_mitigation)
    {
        GET_G77_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G77(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 8, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 9, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 10, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 11, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 12, 7, temp0, temp1);

    // Column 8
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 8);
    }
    else 
    {
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 8);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 10, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 12, 8, temp0, temp1);

    // Column 9
    if constexpr(cond_nan_mitigation)
    {
         GET_Gii_ODD_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 9);
    }
    else 
    {
         GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 9);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 10, 9, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 11, 9, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 12, 9, temp0, temp1);

    // Column 10
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 10);
    }
    else 
    {
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 10);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 10, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 12, 10, temp0, temp1);

    //Column 11
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_ODD_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 11);
    }
    else 
    {
        GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 11);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 12, 11, temp0, temp1);

    //Column 12
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 12);
    }
    else 
    {
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 12);
    }

    /////////////////////////////////// get L, L=invG
    // Column 0
    SET_Lii(matLRe, matLIm, matD, 0);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 1, 0);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 2, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 0, temp0, temp1);

    // Column 1
    SET_Lii(matLRe, matLIm, matD, 1);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 2, 1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 1, temp0, temp1);

    // Column 2
    SET_Lii(matLRe, matLIm, matD, 2);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 3, 2);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 2, temp0, temp1);

    // Column 3
    SET_Lii(matLRe, matLIm, matD, 3);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 4, 3);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 3, temp0, temp1);

    // Column 4
    SET_Lii(matLRe, matLIm, matD, 4);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 5, 4);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 4, temp0, temp1);

    // Column 5
    SET_Lii(matLRe, matLIm, matD, 5);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 6, 5);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 5, temp0, temp1);

    // Column 6
    SET_Lii(matLRe, matLIm, matD, 6);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 7, 6);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 6, temp0, temp1);

    // Column 7
    SET_Lii(matLRe, matLIm, matD, 7);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 8, 7);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 7, temp0, temp1);

    // Column 8
    SET_Lii(matLRe, matLIm, matD, 8);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 9, 8);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 8, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 8, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 8, temp0, temp1);

    // Column 9
    SET_Lii(matLRe, matLIm, matD, 9);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 10, 9);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 9, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 9, temp0, temp1);

    // Column 10
    SET_Lii(matLRe, matLIm, matD, 10);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 11, 10);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 10, temp0, temp1);

    //Column 11 
    SET_Lii(matLRe, matLIm, matD, 11);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 12, 11);

    //Column 12
    SET_Lii(matLRe, matLIm, matD, 12);

    for(i = 0; i < MAT_SQR_SIZE_13; ++i)
    {
        matLReii        = matLRe[i][i];
        matInvBRe[i][i] = _mm512_mul_ps(matLReii,matLReii);
        for(k = (i+1); k < MAT_SQR_SIZE_13; ++k) 
        {
            matLReki = matLRe[k][i];
            matLImki = matLIm[k][i];
            temp1    = _mm512_fmadd_ps(matLReki,matLReki,_mm512_mul_ps(matLImki,matLImki));
            matInvBRe[i][i] = _mm512_add_ps(matInvBRe[i][i],temp1);
        } 
        matInvBIm[i][i] = vzero;
    }

    for(i = 0; i < MAT_SQR_SIZE_13; ++i) 
    {
        for(j = (i+1);j < MAT_SQR_SIZE_13; ++j)   
        {
            matLRejj        = matLRe[j][j];
            matInvBRe[i][j] = _mm512_mul_ps(matLRe[j][i],matLRejj);
            matInvBIm[i][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][i],matLRejj));
            for(k = (j+1); k < MAT_SQR_SIZE_13; k++)
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][i],matLIm[k][i],temp1,temp2);
                matInvBRe[i][j] = _mm512_add_ps(matInvBRe[i][j], temp1);
                matInvBIm[i][j] = _mm512_add_ps(matInvBIm[i][j], temp2);
            }
            matInvBRe[j][i] = matInvBRe[i][j];
            matInvBIm[j][i] = _mm512_sub_ps(vzero, matInvBIm[i][j]);
        }
    }
}

template void 
gms::math 
::mat_inv_cholesky_13x13_16xf32<true,true>(const __m512 matBRe[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                                           const __m512 matBIm[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                                           __m512       matInvBRe[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                                           __m512       matInvBIm[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13]);

template void 
gms::math 
::mat_inv_cholesky_13x13_16xf32<false,true>(const __m512 matBRe[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                                            const __m512 matBIm[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                                            __m512       matInvBRe[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                                            __m512       matInvBIm[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13]);

template void 
gms::math 
::mat_inv_cholesky_13x13_16xf32<true,false>(const __m512 matBRe[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                                            const __m512 matBIm[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                                            __m512       matInvBRe[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                                            __m512       matInvBIm[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13]);

template void 
gms::math 
::mat_inv_cholesky_13x13_16xf32<false,false>(const __m512 matBRe[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                                             const __m512 matBIm[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                                             __m512       matInvBRe[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                                             __m512       matInvBIm[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13]);

template<bool use_prefetching,bool mitigate_nan,
         std::int32_t N_VECS_13>
constexpr 
void 
gms::math
::mat_inv_cholesky_13x13_16xf32_v2(const __m512 matBRe[N_VECS_13][N_VECS_13],
                                   const __m512 matBIm[N_VECS_13][N_VECS_13],
                                   __m512       matInvBRe[N_VECS_13][N_VECS_13],
                                   __m512       matInvBIm[N_VECS_13][N_VECS_13])
{
    static_assert(13==N_VECS_13,"***ERROR*** -- N_VECS_13 != 13!!");
    __ATTR_ALIGN__(64) 
    __m512 matGRe[N_VECS_13][N_VECS_13];
    __ATTR_ALIGN__(64)
    __m512 matGIm[N_VECS_13][N_VECS_13];
    __ATTR_ALIGN__(64)
    __m512 matLRe[N_VECS_13][N_VECS_13];
    __ATTR_ALIGN__(64)
    __m512 matLIm[N_VECS_13][N_VECS_13];
    __ATTR_ALIGN__(64)
    __m512 matD[N_VECS_13];
    __ATTR_ALIGN__(64)
    __m512 matND[N_VECS_13];
    __m512 temp0, temp1, temp2;
    __m512 matLReii;
    __m512 matLReki;
    __m512 matLImki;
    __m512 matLRejj;
    std::int32_t i,j,k;
    const __m512 vzero{_mm512_setzero_ps()};
    const __m512 vneg_one{_mm512_set1_ps(-1.0f)};
    /////////////////////////////////////////////////////////////////////////////////////////////
    if constexpr (static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
    {
        constexpr std::int32_t VEC_PS_LEN = 16;
        const float * __restrict ptr_matBRe{reinterpret_cast<const float * __restrict>(&matBRe)};
        const float * __restrict ptr_matBIm{reinterpret_cast<const float * __restrict>(&matBIm)};
        for(i = 0;i != N_VECS_13; ++i) 
        {
            register std::int32_t outer_idx = i*MAT_SQR_SIZE_13;
            for(j = 0;j != N_VECS_13; ++j)   
            {
                register std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                _mm_prefetch((const char*)&ptr_matBRe[inner_idx],_MM_HINT_T0);
                _mm_prefetch((const char*)&ptr_matBIm[inner_idx],_MM_HINT_T0);
            }
        }
    }
    const bool cond_nan_mitigation = static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true);
     // Column 0
    GET_G00(matGRe, matBRe, matD, matND);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 1);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 2);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 3);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 4);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 5);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 6);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 7);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 8);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 9);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 10);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 11);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 12);

    // Column 1
    if constexpr(cond_nan_mitigation)
    {
        GET_G11_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    else 
    {
        GET_G11(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 2, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 11, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 12, temp0, temp1);

    // Column 2
    if constexpr(cond_nan_mitigation)
    {
        GET_G22_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G22(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 11, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 12, temp0, temp1);

    // Column 3
    if constexpr(cond_nan_mitigation)
    {
        GET_G33_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G33(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 11, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 12, temp0, temp1);

    // Column 4
    if constexpr(cond_nan_mitigation)
    {
        GET_G44_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G44(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 5, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 6, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 7, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 8, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 10, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 12, 4, temp0, temp1);

    // Column 5
    if constexpr(cond_nan_mitigation)
    {
        GET_G55_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G55(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 6, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 7, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 8, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 9, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 10, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 11, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 12, 5, temp0, temp1);

    // Column 6
    if constexpr(cond_nan_mitigation)
    {
        GET_G66_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G66(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 7, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 8, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 10, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 12, 6, temp0, temp1);

    // Column 7
    if constexpr(cond_nan_mitigation)
    {
        GET_G77_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G77(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 8, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 9, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 10, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 11, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 12, 7, temp0, temp1);

    // Column 8
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 8);
    }
    else 
    {
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 8);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 10, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 12, 8, temp0, temp1);

    // Column 9
    if constexpr(cond_nan_mitigation)
    {
         GET_Gii_ODD_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 9);
    }
    else 
    {
         GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 9);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 10, 9, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 11, 9, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 12, 9, temp0, temp1);

    // Column 10
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 10);
    }
    else 
    {
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 10);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 10, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 12, 10, temp0, temp1);

    //Column 11
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_ODD_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 11);
    }
    else 
    {
        GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 11);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 12, 11, temp0, temp1);

    //Column 12
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 12);
    }
    else 
    {
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 12);
    }

    /////////////////////////////////// get L, L=invG
    // Column 0
    SET_Lii(matLRe, matLIm, matD, 0);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 1, 0);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 2, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 0, temp0, temp1);

    // Column 1
    SET_Lii(matLRe, matLIm, matD, 1);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 2, 1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 1, temp0, temp1);

    // Column 2
    SET_Lii(matLRe, matLIm, matD, 2);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 3, 2);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 2, temp0, temp1);

    // Column 3
    SET_Lii(matLRe, matLIm, matD, 3);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 4, 3);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 3, temp0, temp1);

    // Column 4
    SET_Lii(matLRe, matLIm, matD, 4);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 5, 4);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 4, temp0, temp1);

    // Column 5
    SET_Lii(matLRe, matLIm, matD, 5);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 6, 5);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 5, temp0, temp1);

    // Column 6
    SET_Lii(matLRe, matLIm, matD, 6);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 7, 6);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 6, temp0, temp1);

    // Column 7
    SET_Lii(matLRe, matLIm, matD, 7);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 8, 7);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 7, temp0, temp1);

    // Column 8
    SET_Lii(matLRe, matLIm, matD, 8);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 9, 8);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 8, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 8, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 8, temp0, temp1);

    // Column 9
    SET_Lii(matLRe, matLIm, matD, 9);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 10, 9);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 9, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 9, temp0, temp1);

    // Column 10
    SET_Lii(matLRe, matLIm, matD, 10);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 11, 10);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 10, temp0, temp1);

    //Column 11 
    SET_Lii(matLRe, matLIm, matD, 11);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 12, 11);

    //Column 12
    SET_Lii(matLRe, matLIm, matD, 12);

    for(i = 0; i < N_VECS_13; ++i)
    {
        matLReii        = matLRe[i][i];
        matInvBRe[i][i] = _mm512_mul_ps(matLReii,matLReii);
        for(k = (i+1); k < N_VECS_13; ++k) 
        {
            matLReki = matLRe[k][i];
            matLImki = matLIm[k][i];
            temp1    = _mm512_fmadd_ps(matLReki,matLReki,_mm512_mul_ps(matLImki,matLImki));
            matInvBRe[i][i] = _mm512_add_ps(matInvBRe[i][i],temp1);
        } 
        matInvBIm[i][i] = vzero;
    }

    for(i = 0; i < N_VECS_13; ++i) 
    {
        for(j = (i+1);j < N_VECS_13; ++j)   
        {
            matLRejj        = matLRe[j][j];
            matInvBRe[i][j] = _mm512_mul_ps(matLRe[j][i],matLRejj);
            matInvBIm[i][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][i],matLRejj));
            for(k = (j+1); k < N_VECS_13; k++)
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][i],matLIm[k][i],temp1,temp2);
                matInvBRe[i][j] = _mm512_add_ps(matInvBRe[i][j], temp1);
                matInvBIm[i][j] = _mm512_add_ps(matInvBIm[i][j], temp2);
            }
            matInvBRe[j][i] = matInvBRe[i][j];
            matInvBIm[j][i] = _mm512_sub_ps(vzero, matInvBIm[i][j]);
        }
    }
}

template void 
gms::math 
::mat_inv_cholesky_13x13_16xf32_v2<true,true,13>(const __m512 matBRe[13][13],
                                                const __m512 matBIm[13][13],
                                                __m512       matInvBRe[13][13],
                                                __m512       matInvBIm[13][13]);

template void 
gms::math 
::mat_inv_cholesky_13x13_16xf32_v2<false,true,13>(const __m512 matBRe[13][13],
                                                  const __m512 matBIm[13][13],
                                                  __m512       matInvBRe[13][13],
                                                  __m512       matInvBIm[13][13]);

template void 
gms::math 
::mat_inv_cholesky_13x13_16xf32_v2<true,false,13>(const __m512 matBRe[13][13],
                                                  const __m512 matBIm[13][13],
                                                  __m512       matInvBRe[13][13],
                                                  __m512       matInvBIm[13][13]);

template void 
gms::math 
::mat_inv_cholesky_13x13_16xf32_v2<false,false,13>(const __m512 matBRe[13][13],
                                                   const __m512 matBIm[13][13],
                                                   __m512       matInvBRe[13][13],
                                                   __m512       matInvBIm[13][13]);

template<bool use_prefetching,bool mitigate_nan>
void 
gms::math
::mat_inv_cholesky_16x16_16xf32(const __m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                const __m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                __m512       matInvBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                __m512       matInvBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16])
{
    __ATTR_ALIGN__(64) 
    __m512 matGRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];
    __ATTR_ALIGN__(64)
    __m512 matGIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];
    __ATTR_ALIGN__(64)
    __m512 matLRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];
    __ATTR_ALIGN__(64)
    __m512 matLIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];
    __ATTR_ALIGN__(64)
    __m512 matD[MAT_SQR_SIZE_16];
    __ATTR_ALIGN__(64)
    __m512 matND[MAT_SQR_SIZE_16];
    __m512 temp0, temp1, temp2;
    __m512 matLReii;
    __m512 matLReki;
    __m512 matLImki;
    __m512 matLRejj;
    std::int32_t i,j,k;
    const __m512 vzero{_mm512_setzero_ps()};
    const __m512 vneg_one{_mm512_set1_ps(-1.0f)};
    /////////////////////////////////////////////////////////////////////////////////////////////
    if constexpr (static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
    {
        constexpr std::int32_t VEC_PS_LEN = 16;
        const float * __restrict ptr_matBRe{reinterpret_cast<const float * __restrict>(&matBRe)};
        const float * __restrict ptr_matBIm{reinterpret_cast<const float * __restrict>(&matBIm)};
        for(i = 0;i != MAT_SQR_SIZE_16; ++i) 
        {
            register std::int32_t outer_idx = i*MAT_SQR_SIZE_16;
            for(j = 0;j != MAT_SQR_SIZE_16; ++j)   
            {
                register std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                _mm_prefetch((const char*)&ptr_matBRe[inner_idx],_MM_HINT_T0);
                _mm_prefetch((const char*)&ptr_matBIm[inner_idx],_MM_HINT_T0);
            }
        }
    }
    const bool cond_nan_mitigation = static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
     // Column 0
    GET_G00(matGRe, matBRe, matD, matND);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 1);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 2);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 3);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 4);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 5);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 6);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 7);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 8);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 9);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 10);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 11);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 12);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 13);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 14);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 15);

    // Column 1
    if constexpr(cond_nan_mitigation)
    {
        GET_G11_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    else 
    {
        GET_G11(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 2, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 11, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 12, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 13, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 14, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 15, temp0, temp1);

    // Column 2
    if constexpr(cond_nan_mitigation)
    {
        GET_G22_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G22(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 11, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 12, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 13, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 14, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 15, temp0, temp1);

    // Column 3
    if constexpr(cond_nan_mitigation)
    {
        GET_G33_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G33(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 11, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 12, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 13, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 14, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 15, temp0, temp1);

    // Column 4
    if constexpr(cond_nan_mitigation)
    {
        GET_G44_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G44(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 5, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 6, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 7, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 8, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 10, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 12, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 13, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 14, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 15, 4, temp0, temp1);

    // Column 5
    if constexpr(cond_nan_mitigation)
    {
        GET_G55_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G55(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 6, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 7, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 8, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 9, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 10, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 11, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 12, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 13, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 14, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 15, 5, temp0, temp1);

    // Column 6
    if constexpr(cond_nan_mitigation)
    {
        GET_G66_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        GET_G66(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 7, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 8, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 10, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 12, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 13, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 14, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 15, 6, temp0, temp1);

    // Column 7
    if constexpr(cond_nan_mitigation)
    {
        GET_G77_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        GET_G77(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 8, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 9, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 10, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 11, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 12, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 13, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 14, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 15, 7, temp0, temp1);

    // Column 8
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 8);
    }
    else 
    {
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 8);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 10, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 12, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 13, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 14, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 15, 8, temp0, temp1);

    // Column 9
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_ODD_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 9);
    }
    else 
    {
        GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 9);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 10, 9, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 11, 9, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 12, 9, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 13, 9, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 14, 9, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 15, 9, temp0, temp1);

    // Column 10
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 10);
    }
    else 
    {
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 10);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 10, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 12, 10, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 13, 10, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 14, 10, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 15, 10, temp0, temp1);

    // Column 11
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_ODD_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 11);
    }
    else 
    {
        GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 11);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 12, 11, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 13, 11, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 14, 11, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 15, 11, temp0, temp1);

    // Column 12
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 12);
    }
    else 
    {
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 12);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 13, 12, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 14, 12, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 15, 12, temp0, temp1);

    // Column 13
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_ODD_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 13);
    }
    else 
    {
        GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 13);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 14, 13, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 15, 13, temp0, temp1);

    // Column 14
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 14);
    }
    else 
    {
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 14);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 15, 14, temp0, temp1);

    // Column 15
    if constexpr(cond_nan_mitigation)
    {
        GET_Gii_ODD_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 15);
    }
    else 
    {
        GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 15);
    }

    /////////////////////////////////// get L, L=invG
    // Column 0
    SET_Lii(matLRe, matLIm, matD, 0);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 1, 0);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 2, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 0, temp0, temp1);

    // Column 1
    SET_Lii(matLRe, matLIm, matD, 1);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 2, 1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 1, temp0, temp1);

    // Column 2
    SET_Lii(matLRe, matLIm, matD, 2);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 3, 2);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 2, temp0, temp1);

    // Column 3
    SET_Lii(matLRe, matLIm, matD, 3);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 4, 3);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 3, temp0, temp1);

    // Column 4
    SET_Lii(matLRe, matLIm, matD, 4);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 5, 4);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 4, temp0, temp1);

    // Column 5
    SET_Lii(matLRe, matLIm, matD, 5);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 6, 5);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 5, temp0, temp1);

    // Column 6
    SET_Lii(matLRe, matLIm, matD, 6);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 7, 6);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 6, temp0, temp1);

    // Column 7
    SET_Lii(matLRe, matLIm, matD, 7);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 8, 7);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 7, temp0, temp1);

    // Column 8
    SET_Lii(matLRe, matLIm, matD, 8);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 9, 8);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 8, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 8, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 8, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 8, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 8, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 8, temp0, temp1);

    // Column 9
    SET_Lii(matLRe, matLIm, matD, 9);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 10, 9);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 9, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 9, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 9, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 9, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 9, temp0, temp1);

    // Column 10
    SET_Lii(matLRe, matLIm, matD, 10);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 11, 10);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 10, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 10, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 10, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 10, temp0, temp1);

    // Column 11
    SET_Lii(matLRe, matLIm, matD, 11);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 12, 11);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 11, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 11, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 11, temp0, temp1);

    // Column 12
    SET_Lii(matLRe, matLIm, matD, 12);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 13, 12);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 12, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 12, temp0, temp1);

    // Column 13
    SET_Lii(matLRe, matLIm, matD, 13);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 14, 13);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 13, temp0, temp1);

    // Column 14
    SET_Lii(matLRe, matLIm, matD, 14);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 15, 14);

    // Column 15
    SET_Lii(matLRe, matLIm, matD, 15);

    for(i = 0; i < MAT_SQR_SIZE_16; ++i)
    {
        matLReii        = matLRe[i][i];
        matInvBRe[i][i] = _mm512_mul_ps(matLReii,matLReii);
        for(k = (i+1); k < MAT_SQR_SIZE_16; ++k) 
        {
            matLReki = matLRe[k][i];
            matLImki = matLIm[k][i];
            temp1    = _mm512_fmadd_ps(matLReki,matLReki,_mm512_mul_ps(matLImki,matLImki));
            matInvBRe[i][i] = _mm512_add_ps(matInvBRe[i][i],temp1);
        } 
        matInvBIm[i][i] = vzero;
    }

    for(i = 0; i < MAT_SQR_SIZE_16; ++i) 
    {
        for(j = (i+1);j < MAT_SQR_SIZE_16; ++j)   
        {
            matLRejj        = matLRe[j][j];
            matInvBRe[i][j] = _mm512_mul_ps(matLRe[j][i],matLRejj);
            matInvBIm[i][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][i],matLRejj));
            for(k = (j+1); k < MAT_SQR_SIZE_16; k++)
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][i],matLIm[k][i],temp1,temp2);
                matInvBRe[i][j] = _mm512_add_ps(matInvBRe[i][j], temp1);
                matInvBIm[i][j] = _mm512_add_ps(matInvBIm[i][j], temp2);
            }
            matInvBRe[j][i] = matInvBRe[i][j];
            matInvBIm[j][i] = _mm512_sub_ps(vzero, matInvBIm[i][j]);
        }
    }
}

template void 
gms::math 
::mat_inv_cholesky_16x16_16xf32<true,true>(const __m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                           const __m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                           __m512       matInvBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                           __m512       matInvBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16]);

template void 
gms::math 
::mat_inv_cholesky_16x16_16xf32<false,true>(const __m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                            const __m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                            __m512       matInvBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                            __m512       matInvBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16]);

template void 
gms::math 
::mat_inv_cholesky_16x16_16xf32<true,false>(const __m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                            const __m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                            __m512       matInvBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                            __m512       matInvBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16]);

template void 
gms::math 
::mat_inv_cholesky_16x16_16xf32<false,false>(const __m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                             const __m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                             __m512       matInvBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                             __m512       matInvBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16]);

/* Dedicated version for the TSC instrumentation profiling*/

#include "GMS_machine_utils.h"

template<bool use_prefetching,bool mitigate_nan>
void 
gms::math
::mat_inv_cholesky_16x16_16xf32_tsc_instr(const __m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                          const __m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                          __m512       matInvBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                          __m512       matInvBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                          mat_inv_chol_tsc_instr::tsc_instrumentation_block_t<MAT_SQR_SIZE_16> &tsc_instr_block)
{
    using namespace gms::common;
    __ATTR_ALIGN__(64) 
    __m512 matGRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];
    __ATTR_ALIGN__(64)
    __m512 matGIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];
    __ATTR_ALIGN__(64)
    __m512 matLRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];
    __ATTR_ALIGN__(64)
    __m512 matLIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];
    __ATTR_ALIGN__(64)
    __m512 matD[MAT_SQR_SIZE_16];
    __ATTR_ALIGN__(64)
    __m512 matND[MAT_SQR_SIZE_16];
    __m512 temp0, temp1, temp2;
    __m512 matLReii;
    __m512 matLReki;
    __m512 matLImki;
    __m512 matLRejj;
    [[maybe_unused]] volatile std::uint64_t rdtscp_fenced_warmup;
    std::int32_t i,j,k;
    const __m512 vzero{_mm512_setzero_ps()};
    const __m512 vneg_one{_mm512_set1_ps(-1.0f)};
    /////////////////////////////////////////////////////////////////////////////////////////////
    rdtscp_fenced_warmup = rdtscp_fenced();
    if constexpr (static_cast<std::int32_t>(use_prefetching)==static_cast<std::int32_t>(true))
    {
        constexpr std::int32_t VEC_PS_LEN = 16;
        const float * __restrict ptr_matBRe{reinterpret_cast<const float * __restrict>(&matBRe)};
        const float * __restrict ptr_matBIm{reinterpret_cast<const float * __restrict>(&matBIm)};
        tsc_instr_block.m_region_prefetch_s = rdtscp_fenced();
        for(i = 0;i != MAT_SQR_SIZE_16; ++i) 
        {
            register std::int32_t outer_idx = i*MAT_SQR_SIZE_16;
            for(j = 0;j != MAT_SQR_SIZE_16; ++j)   
            {
                register std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                _mm_prefetch((const char*)&ptr_matBRe[inner_idx],_MM_HINT_T0);
                _mm_prefetch((const char*)&ptr_matBIm[inner_idx],_MM_HINT_T0);
            }
        }
        tsc_instr_block.m_region_prefetch_e = rdtscp_fenced();
    }
    const bool cond_nan_mitigation = static_cast<std::int32_t>(mitigate_nan)==static_cast<std::int32_t>(true);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
     // Column 0
    tsc_instr_block.m_region_g00_s = rdtscp_fenced();
    GET_G00(matGRe, matBRe, matD, matND);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 1);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 2);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 3);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 4);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 5);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 6);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 7);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 8);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 9);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 10);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 11);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 12);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 13);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 14);
    GET_G_COL0(matGRe, matGIm, matBRe, matBIm, matD, 15);
    tsc_instr_block.m_region_g00_e = rdtscp_fenced();
    // Column 1
    if constexpr(cond_nan_mitigation)
    {   
        tsc_instr_block.m_region_g11_s = rdtscp_fenced();
        GET_G11_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    else 
    {   
        tsc_instr_block.m_region_g11_s = rdtscp_fenced();
        GET_G11(matGRe, matGIm, matBRe, matD, matND, temp0);
    }
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 2, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 11, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 12, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 13, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 14, temp0, temp1);
    GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, 15, temp0, temp1);
    tsc_instr_block.m_region_g11_e = rdtscp_fenced();
    // Column 2
    if constexpr(cond_nan_mitigation)
    {
        tsc_instr_block.m_region_g22_s = rdtscp_fenced();
        GET_G22_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        tsc_instr_block.m_region_g22_s = rdtscp_fenced();
        GET_G22(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 11, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 12, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 13, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 14, temp0, temp1);
    GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, 15, temp0, temp1);
    tsc_instr_block.m_region_g22_e = rdtscp_fenced();
    // Column 3
    if constexpr(cond_nan_mitigation)
    {
        tsc_instr_block.m_region_g33_s = rdtscp_fenced();
        GET_G33_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        tsc_instr_block.m_region_g33_s = rdtscp_fenced();
        GET_G33(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 11, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 12, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 13, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 14, temp0, temp1);
    GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, 15, temp0, temp1);
    tsc_instr_block.m_region_g33_e = rdtscp_fenced();
    // Column 4
    if constexpr(cond_nan_mitigation)
    {
        tsc_instr_block.m_region_g44_s = rdtscp_fenced();
        GET_G44_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        tsc_instr_block.m_region_g44_s = rdtscp_fenced();
        GET_G44(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 5, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 6, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 7, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 8, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 10, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 12, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 13, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 14, 4, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 15, 4, temp0, temp1);
    tsc_instr_block.m_region_g44_e = rdtscp_fenced();
    // Column 5
    if constexpr(cond_nan_mitigation)
    {
        tsc_instr_block.m_region_g55_s = rdtscp_fenced();
        GET_G55_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        tsc_instr_block.m_region_g55_s = rdtscp_fenced();
        GET_G55(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 6, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 7, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 8, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 9, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 10, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 11, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 12, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 13, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 14, 5, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 15, 5, temp0, temp1);
    tsc_instr_block.m_region_g55_e = rdtscp_fenced();
    // Column 6
    if constexpr(cond_nan_mitigation)
    {
        tsc_instr_block.m_region_g66_s = rdtscp_fenced();
        GET_G66_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    else 
    {
        tsc_instr_block.m_region_g66_s = rdtscp_fenced();
        GET_G66(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 7, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 8, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 10, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 12, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 13, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 14, 6, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 15, 6, temp0, temp1);
    tsc_instr_block.m_region_g66_e = rdtscp_fenced();
    // Column 7
    if constexpr(cond_nan_mitigation)
    {
        tsc_instr_block.m_region_g77_s = rdtscp_fenced();
        GET_G77_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    else 
    {
        tsc_instr_block.m_region_g77_s = rdtscp_fenced();
        GET_G77(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 8, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 9, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 10, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 11, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 12, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 13, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 14, 7, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 15, 7, temp0, temp1);
    tsc_instr_block.m_region_g77_e = rdtscp_fenced();
    // Column 8
    if constexpr(cond_nan_mitigation)
    {
        tsc_instr_block.m_region_gii_even_col8_s = rdtscp_fenced();
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 8);
    }
    else 
    {
        tsc_instr_block.m_region_gii_even_col8_s = rdtscp_fenced();
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 8);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 9, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 10, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 12, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 13, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 14, 8, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 15, 8, temp0, temp1);
    tsc_instr_block.m_region_gii_even_col8_e = rdtscp_fenced();
    // Column 9
    if constexpr(cond_nan_mitigation)
    {
        tsc_instr_block.m_region_gii_odd_col9_s = rdtscp_fenced();
        GET_Gii_ODD_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 9);
    }
    else 
    {
        tsc_instr_block.m_region_gii_odd_col9_s = rdtscp_fenced();
        GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 9);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 10, 9, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 11, 9, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 12, 9, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 13, 9, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 14, 9, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 15, 9, temp0, temp1);
    tsc_instr_block.m_region_gii_odd_col9_e = rdtscp_fenced();
    // Column 10
    if constexpr(cond_nan_mitigation)
    {
        tsc_instr_block.m_region_gii_even_col10_s = rdtscp_fenced();
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 10);
    }
    else 
    {
        tsc_instr_block.m_region_gii_even_col10_s = rdtscp_fenced();
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 10);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 11, 10, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 12, 10, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 13, 10, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 14, 10, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 15, 10, temp0, temp1);
    tsc_instr_block.m_region_gii_even_col10_e = rdtscp_fenced();
    // Column 11
    if constexpr(cond_nan_mitigation)
    {
        tsc_instr_block.m_region_gii_odd_col11_s = rdtscp_fenced();
        GET_Gii_ODD_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 11);
    }
    else 
    {
        tsc_instr_block.m_region_gii_odd_col11_s = rdtscp_fenced();
        GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 11);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 12, 11, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 13, 11, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 14, 11, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 15, 11, temp0, temp1);
    tsc_instr_block.m_region_gii_odd_col11_e = rdtscp_fenced();
    // Column 12
    if constexpr(cond_nan_mitigation)
    {
        tsc_instr_block.m_region_gii_even_col12_s = rdtscp_fenced();
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 12);
    }
    else 
    {
        tsc_instr_block.m_region_gii_even_col12_s = rdtscp_fenced();
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 12);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 13, 12, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 14, 12, temp0, temp1);
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 15, 12, temp0, temp1);
    tsc_instr_block.m_region_gii_even_col12_e = rdtscp_fenced();
    // Column 13
    if constexpr(cond_nan_mitigation)
    {
        tsc_instr_block.m_region_gii_odd_col13_s = rdtscp_fenced();
        GET_Gii_ODD_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 13);
    }
    else 
    {   
        tsc_instr_block.m_region_gii_odd_col13_s = rdtscp_fenced();
        GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 13);
    }
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 14, 13, temp0, temp1);
    GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, 15, 13, temp0, temp1);
    tsc_instr_block.m_region_gii_odd_col13_e = rdtscp_fenced();
    // Column 14
    if constexpr(cond_nan_mitigation)
    {
        tsc_instr_block.m_region_gii_even_col14_s = rdtscp_fenced();
        GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 14);
    }
    else 
    {
        tsc_instr_block.m_region_gii_even_col14_s = rdtscp_fenced();
        GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 14);
    }
    GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, 15, 14, temp0, temp1);
    tsc_instr_block.m_region_gii_even_col14_e = rdtscp_fenced();
    // Column 15
    if constexpr(cond_nan_mitigation)
    {
        tsc_instr_block.m_region_gii_odd_col15_s = rdtscp_fenced();
        GET_Gii_ODD_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 15);
    }
    else 
    {
        tsc_instr_block.m_region_gii_odd_col15_s = rdtscp_fenced();
        GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 15);
    }
    tsc_instr_block.m_region_gii_odd_col15_e = rdtscp_fenced();
    /////////////////////////////////// get L, L=invG
    // Column 0
    tsc_instr_block.m_region_lii_col0_s = rdtscp_fenced();
    SET_Lii(matLRe, matLIm, matD, 0);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 1, 0);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 2, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 0, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 0, temp0, temp1);
    tsc_instr_block.m_region_lii_col0_e = rdtscp_fenced();
    // Column 1
    tsc_instr_block.m_region_lii_col1_s = rdtscp_fenced();
    SET_Lii(matLRe, matLIm, matD, 1);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 2, 1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 3, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 1, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 1, temp0, temp1);
    tsc_instr_block.m_region_lii_col1_e = rdtscp_fenced();
    // Column 2
    tsc_instr_block.m_region_lii_col2_s = rdtscp_fenced();
    SET_Lii(matLRe, matLIm, matD, 2);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 3, 2);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 4, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 2, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 2, temp0, temp1);
    tsc_instr_block.m_region_lii_col2_e = rdtscp_fenced();
    // Column 3
    tsc_instr_block.m_region_lii_col3_s = rdtscp_fenced();
    SET_Lii(matLRe, matLIm, matD, 3);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 4, 3);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 5, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 3, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 3, temp0, temp1);
    tsc_instr_block.m_region_lii_col3_e = rdtscp_fenced();
    // Column 4
    tsc_instr_block.m_region_lii_col4_s = rdtscp_fenced();
    SET_Lii(matLRe, matLIm, matD, 4);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 5, 4);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 6, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 4, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 4, temp0, temp1);
    tsc_instr_block.m_region_lii_col4_e = rdtscp_fenced();
    // Column 5
    tsc_instr_block.m_region_lii_col5_s = rdtscp_fenced();
    SET_Lii(matLRe, matLIm, matD, 5);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 6, 5);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 7, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 5, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 5, temp0, temp1);
    tsc_instr_block.m_region_lii_col5_e = rdtscp_fenced();
    // Column 6
    tsc_instr_block.m_region_lii_col6_s = rdtscp_fenced();
    SET_Lii(matLRe, matLIm, matD, 6);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 7, 6);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 8, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 6, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 6, temp0, temp1);
    tsc_instr_block.m_region_lii_col6_e = rdtscp_fenced();
    // Column 7
    tsc_instr_block.m_region_lii_col7_s = rdtscp_fenced();
    SET_Lii(matLRe, matLIm, matD, 7);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 8, 7);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 9, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 7, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 7, temp0, temp1);
    tsc_instr_block.m_region_lii_col7_e = rdtscp_fenced();
    // Column 8
    tsc_instr_block.m_region_lii_col8_s = rdtscp_fenced();
    SET_Lii(matLRe, matLIm, matD, 8);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 9, 8);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 10, 8, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 8, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 8, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 8, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 8, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 8, temp0, temp1);
    tsc_instr_block.m_region_lii_col8_e = rdtscp_fenced();
    // Column 9
    tsc_instr_block.m_region_lii_col9_s = rdtscp_fenced();
    SET_Lii(matLRe, matLIm, matD, 9);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 10, 9);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 11, 9, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 9, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 9, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 9, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 9, temp0, temp1);
    tsc_instr_block.m_region_lii_col9_e = rdtscp_fenced();
    // Column 10
    tsc_instr_block.m_region_lii_col10_s = rdtscp_fenced();
    SET_Lii(matLRe, matLIm, matD, 10);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 11, 10);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 12, 10, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 10, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 10, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 10, temp0, temp1);
    tsc_instr_block.m_region_lii_col10_e = rdtscp_fenced();
    // Column 11
    tsc_instr_block.m_region_lii_col11_s = rdtscp_fenced();
    SET_Lii(matLRe, matLIm, matD, 11);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 12, 11);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 13, 11, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 11, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 11, temp0, temp1);
    tsc_instr_block.m_region_lii_col11_e = rdtscp_fenced();
    // Column 12
    tsc_instr_block.m_region_lii_col12_s = rdtscp_fenced();
    SET_Lii(matLRe, matLIm, matD, 12);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 13, 12);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 14, 12, temp0, temp1);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 12, temp0, temp1);
    tsc_instr_block.m_region_lii_col12_e = rdtscp_fenced();
    // Column 13
    tsc_instr_block.m_region_lii_col13_s = rdtscp_fenced();
    SET_Lii(matLRe, matLIm, matD, 13);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 14, 13);
    GET_L_ji(matLRe, matLIm, matGRe, matGIm, matND, 15, 13, temp0, temp1);
    tsc_instr_block.m_region_lii_col13_e = rdtscp_fenced();
    // Column 14
    tsc_instr_block.m_region_lii_col14_s = rdtscp_fenced();
    SET_Lii(matLRe, matLIm, matD, 14);
    GET_L_i1i(matLRe, matLIm, matGRe, matGIm, matND, 15, 14);
    tsc_instr_block.m_region_lii_col14_e = rdtscp_fenced();
    // Column 15
    tsc_instr_block.m_region_lii_col15_s = rdtscp_fenced();
    SET_Lii(matLRe, matLIm, matD, 15);
    tsc_instr_block.m_region_lii_col15_e = rdtscp_fenced();

    for(i = 0; i < MAT_SQR_SIZE_16; ++i)
    {
        tsc_instr_block.m_region_loop2D_s[i] = rdtscp_fenced();
        matLReii        = matLRe[i][i];
        matInvBRe[i][i] = _mm512_mul_ps(matLReii,matLReii);
        for(k = (i+1); k < MAT_SQR_SIZE_16; ++k) 
        {
            matLReki = matLRe[k][i];
            matLImki = matLIm[k][i];
            temp1    = _mm512_fmadd_ps(matLReki,matLReki,_mm512_mul_ps(matLImki,matLImki));
            matInvBRe[i][i] = _mm512_add_ps(matInvBRe[i][i],temp1);
        } 
        matInvBIm[i][i] = vzero;
        tsc_instr_block.m_region_loop2D_e[i] = rdtscp_fenced();
    }

    for(i = 0; i < MAT_SQR_SIZE_16; ++i) 
    {
        tsc_instr_block.m_region_loop3D_s[i] = rdtscp_fenced();
        for(j = (i+1);j < MAT_SQR_SIZE_16; ++j)   
        {
            matLRejj        = matLRe[j][j];
            matInvBRe[i][j] = _mm512_mul_ps(matLRe[j][i],matLRejj);
            matInvBIm[i][j] = _mm512_sub_ps(vzero,_mm512_mul_ps(matLIm[j][i],matLRejj));
            for(k = (j+1); k < MAT_SQR_SIZE_16; k++)
            {
                GET_AxBH(matLRe[k][j],matLIm[k][j],matLRe[k][i],matLIm[k][i],temp1,temp2);
                matInvBRe[i][j] = _mm512_add_ps(matInvBRe[i][j], temp1);
                matInvBIm[i][j] = _mm512_add_ps(matInvBIm[i][j], temp2);
            }
            matInvBRe[j][i] = matInvBRe[i][j];
            matInvBIm[j][i] = _mm512_sub_ps(vzero, matInvBIm[i][j]);
        }
        tsc_instr_block.m_region_loop3D_e[i] = rdtscp_fenced();
    }
}

template void 
gms::math 
::mat_inv_cholesky_16x16_16xf32_tsc_instr<true,true>(const __m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                           const __m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                           __m512       matInvBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                           __m512       matInvBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                           mat_inv_chol_tsc_instr::tsc_instrumentation_block_t<MAT_SQR_SIZE_16> &tsc_instr_block);

template void 
gms::math 
::mat_inv_cholesky_16x16_16xf32_tsc_instr<false,true>(const __m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                            const __m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                            __m512       matInvBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                            __m512       matInvBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                            mat_inv_chol_tsc_instr::tsc_instrumentation_block_t<MAT_SQR_SIZE_16> &tsc_instr_block);

template void 
gms::math 
::mat_inv_cholesky_16x16_16xf32_tsc_instr<true,false>(const __m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                            const __m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                            __m512       matInvBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                            __m512       matInvBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                            mat_inv_chol_tsc_instr::tsc_instrumentation_block_t<MAT_SQR_SIZE_16> &tsc_instr_block);

template void 
gms::math 
::mat_inv_cholesky_16x16_16xf32_tsc_instr<false,false>(const __m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                             const __m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                             __m512       matInvBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                             __m512       matInvBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                             mat_inv_chol_tsc_instr::tsc_instrumentation_block_t<MAT_SQR_SIZE_16> &tsc_instr_block);