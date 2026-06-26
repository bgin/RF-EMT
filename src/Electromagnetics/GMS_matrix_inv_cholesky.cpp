
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

#define GET_G00_RSQR28(matGRe, natBRe, matD, matND)\
{\
   matD[0]   = _mm512_rsqrt28_ps(matBRe[0][0]);\
   matND[0]  = _mm512_sub_ps(_mm512_setzero_ps(),matD[0]);\
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
    do{\
       GET_AxAH(matGRe[1][0], matGIm[1][0], temp0);\
       matGRe[1][1] = _mm512_sub_ps(matBRe[1][1], temp0);\
       matD[1] = _mm512_rsqrt14_ps(matGRe[1][1]);\
       matND[1] = _mm512_sub_ps(vzero, matD[1]);\
    }while(0);\
}

// get G11
#define GET_G11_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0)\
{\
    do{\
       GET_AxAH(matGRe[1][0], matGIm[1][0], temp0);\
       matGRe[1][1] = _mm512_sub_ps(matBRe[1][1], temp0);\
       __mmask16 k1 = _mm512_cmp_ps_mask(matGRe[1][1],vzero,_CMP_LT_OQ);\
       matD[1] = _mm512_rsqrt14_ps(_mm512_mask_blend_ps(k1,matGRe[1][1],_mm512_mul_ps(vneg_one,matGRe[1][1])));\
       matND[1] = _mm512_sub_ps(vzero, matD[1]);\
    }while(0);\
}

// This version albeit ~4.66 times slower than VRSQRT14 version achieves 
// minimize an approximation error of 2^-14 propagated by the usage of 
// VRSQRT14 machine code instruction.
#define GET_G11_MIN_APPROX_ERROR(matGRe, matGIm, matBRe, matD, matND, temp0)\
{\
    do{\
       GET_AxAH(matGRe[1][0], matGIm[1][0], temp0);\
       matGRe[1][1] = _mm512_sub_ps(matBRe[1][1], temp0);\
       __mmask16 k1 = _mm512_cmp_ps_mask(matGRe[1][1],vzero,_CMP_LT_OQ);\
       matD[1] = _mm512_div_ps(vone,_mm512_sqrt_ps(_mm512_mask_blend_ps(k1,matGRe[1][1],_mm512_mul_ps(vneg_one,matGRe[1][1]))));\
       matND[1] = _mm512_sub_ps(vzero, matD[1]);\
    }while(0);\
}

// get column 1 of matrix G
#define GET_G_COL1(matGRe, matGIm, matBRe, matBIm, matD, j, temp0, temp1)\
{\
    do{\
        GET_AxBH(matGRe[j][0], matGIm[j][0], matGRe[1][0], matGIm[1][0], temp0, temp1);\
        matGRe[j][1] = _mm512_sub_ps(matBRe[j][1], temp0);\
        matGIm[j][1] = _mm512_sub_ps(matBIm[j][1], temp1);\
        matGRe[j][1] = _mm512_mul_ps(matGRe[j][1], matD[1]);\
        matGIm[j][1] = _mm512_mul_ps(matGIm[j][1], matD[1]);\
    }while(0);\
}

// get G22
#define GET_G22(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
    do{\
        GET_AxAH(matGRe[2][0], matGIm[2][0], temp0);\
        GET_AxAH(matGRe[2][1], matGIm[2][1], temp1);\
        matGRe[2][2] = _mm512_sub_ps(matBRe[2][2], temp0);\
        matGRe[2][2] = _mm512_sub_ps(matGRe[2][2], temp1);\
        matD[2] = _mm512_rsqrt14_ps(matGRe[2][2]);\
        matND[2] = _mm512_sub_ps(vzero, matD[2]);\
    }while(0);\
}

// get G22
#define GET_G22_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
    do{\
        GET_AxAH(matGRe[2][0], matGIm[2][0], temp0);\
        GET_AxAH(matGRe[2][1], matGIm[2][1], temp1);\
        matGRe[2][2] = _mm512_sub_ps(matBRe[2][2], temp0);\
        matGRe[2][2] = _mm512_sub_ps(matGRe[2][2], temp1);\
        __mmask16 k2 = _mm512_cmp_ps_mask(matGRe[2][2],vzero,_CMP_LT_OQ);\
        matD[2] = _mm512_rsqrt14_ps(_mm512_mask_blend_ps(k2,matGRe[2][2],_mm512_mul_ps(vneg_one,matGRe[2][2])));\
        matND[2] = _mm512_sub_ps(vzero, matD[2]);\
    }while(0);\
}

#define GET_G22_SAFE_MIN_APPROX_ERROR(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
    do{\
        GET_AxAH(matGRe[2][0], matGIm[2][0], temp0);\
        GET_AxAH(matGRe[2][1], matGIm[2][1], temp1);\
        matGRe[2][2] = _mm512_sub_ps(matBRe[2][2], temp0);\
        matGRe[2][2] = _mm512_sub_ps(matGRe[2][2], temp1);\
        __mmask16 k2 = _mm512_cmp_ps_mask(matGRe[2][2],vzero,_CMP_LT_OQ);\
        matD[2] = _mm512_div_ps(vone,_mm512_sqrt_ps(_mm512_mask_blend_ps(k2,matGRe[2][2],_mm512_mul_ps(vneg_one,matGRe[2][2]))));\
        matND[2] = _mm512_sub_ps(vzero, matD[2]);\
    }while(0);\
}

// get column 2 of matrix G
#define GET_G_COL2(matGRe, matGIm, matBRe, matBIm, matD, j, temp0, temp1)\
{\
    do{\
         GET_AxBH(matGRe[j][0], matGIm[j][0], matGRe[2][0], matGIm[2][0], temp0, temp1);\
         matGRe[j][2] = _mm512_sub_ps(matBRe[j][2], temp0);\
         matGIm[j][2] = _mm512_sub_ps(matBIm[j][2], temp1);\
         GET_AxBH(matGRe[j][1], matGIm[j][1], matGRe[2][1], matGIm[2][1], temp0, temp1);\
         matGRe[j][2] = _mm512_sub_ps(matGRe[j][2], temp0);\
         matGIm[j][2] = _mm512_sub_ps(matGIm[j][2], temp1);\
         matGRe[j][2] = _mm512_mul_ps(matGRe[j][2], matD[2]);\
         matGIm[j][2] = _mm512_mul_ps(matGIm[j][2], matD[2]);\
    }while(0);\
}

// get G33
#define GET_G33(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
    do{\
        GET_AxAH(matGRe[3][0], matGIm[3][0], temp0);\
        GET_AxAH(matGRe[3][1], matGIm[3][1], temp1);\
        GET_AxAH(matGRe[3][2], matGIm[3][2], temp2);\
        matGRe[3][3] = _mm512_sub_ps(matBRe[3][3], temp0);\
        matGRe[3][3] = _mm512_sub_ps(matGRe[3][3], temp1);\
        matGRe[3][3] = _mm512_sub_ps(matGRe[3][3], temp2);\
        matD[3] = _mm512_rsqrt14_ps(matGRe[3][3]);\
        matND[3] = _mm512_sub_ps(vzero, matD[3]);\
    }while(0);\
}

// get G33
#define GET_G33_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
    do{\
        GET_AxAH(matGRe[3][0], matGIm[3][0], temp0);\
        GET_AxAH(matGRe[3][1], matGIm[3][1], temp1);\
        GET_AxAH(matGRe[3][2], matGIm[3][2], temp2);\
        matGRe[3][3] = _mm512_sub_ps(matBRe[3][3], temp0);\
        matGRe[3][3] = _mm512_sub_ps(matGRe[3][3], temp1);\
        matGRe[3][3] = _mm512_sub_ps(matGRe[3][3], temp2);\
        __mmask16 k3 = _mm512_cmp_ps_mask(matGRe[3][3],vzero,_CMP_LT_OQ);\
        matD[3] = _mm512_rsqrt14_ps(_mm512_mask_blend_ps(k3,matGRe[3][3],_mm512_mul_ps(vneg_one,matGRe[3][3])));\
        matND[3] = _mm512_sub_ps(vzero, matD[3]);\
    }while(0);\
}

#define GET_G33_SAFE_MIN_APPROX_ERROR(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
    do{\
        GET_AxAH(matGRe[3][0], matGIm[3][0], temp0);\
        GET_AxAH(matGRe[3][1], matGIm[3][1], temp1);\
        GET_AxAH(matGRe[3][2], matGIm[3][2], temp2);\
        matGRe[3][3] = _mm512_sub_ps(matBRe[3][3], temp0);\
        matGRe[3][3] = _mm512_sub_ps(matGRe[3][3], temp1);\
        matGRe[3][3] = _mm512_sub_ps(matGRe[3][3], temp2);\
        __mmask16 k3 = _mm512_cmp_ps_mask(matGRe[3][3],vzero,_CMP_LT_OQ);\
        matD[3] = _mm512_div_ps(vone,_mm512_sqrt_ps(_mm512_mask_blend_ps(k3,matGRe[3][3],_mm512_mul_ps(vneg_one,matGRe[3][3]))));\
        matND[3] = _mm512_sub_ps(vzero, matD[3]);\
    }while(0);\
}

// get column 3 of matrix G
#define GET_G_COL3(matGRe, matGIm, matBRe, matBIm, matD, j, temp0, temp1)\
{\
    do{\
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
    }while(0);\
}

#define GET_G44(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
    do{\
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
    }while(0);\
}

// get G44
#define GET_G44_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
    do{\
        GET_AxAH(matGRe[4][0], matGIm[4][0], temp0);\
        GET_AxAH(matGRe[4][1], matGIm[4][1], temp1);\
        matGRe[4][4] = _mm512_sub_ps(matBRe[4][4], temp0);\
        matGRe[4][4] = _mm512_sub_ps(matGRe[4][4], temp1);\
        GET_AxAH(matGRe[4][2], matGIm[4][2], temp0);\
        GET_AxAH(matGRe[4][3], matGIm[4][3], temp1);\
        matGRe[4][4] = _mm512_sub_ps(matGRe[4][4], temp0);\
        matGRe[4][4] = _mm512_sub_ps(matGRe[4][4], temp1);\
        __mmask16 k4 = _mm512_cmp_ps_mask(matGRe[4][4],vzero,_CMP_LT_OQ);\
        matD[4] = _mm512_rsqrt14_ps(_mm512_mask_blend_ps(k3,matGRe[4][4],_mm512_mul_ps(vneg_one,matGRe[4][4])));\
        matND[4] = _mm512_sub_ps(vzero, matD[4]);\
    }while(0);\
}

#define GET_G44_SAFE_MIN_APPROX_ERROR(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
    do{\
        GET_AxAH(matGRe[4][0], matGIm[4][0], temp0);\
        GET_AxAH(matGRe[4][1], matGIm[4][1], temp1);\
        matGRe[4][4] = _mm512_sub_ps(matBRe[4][4], temp0);\
        matGRe[4][4] = _mm512_sub_ps(matGRe[4][4], temp1);\
        GET_AxAH(matGRe[4][2], matGIm[4][2], temp0);\
        GET_AxAH(matGRe[4][3], matGIm[4][3], temp1);\
        matGRe[4][4] = _mm512_sub_ps(matGRe[4][4], temp0);\
        matGRe[4][4] = _mm512_sub_ps(matGRe[4][4], temp1);\
        __mmask16 k4 = _mm512_cmp_ps_mask(matGRe[4][4],vzero,_CMP_LT_OQ);\
        matD[4] = _mm512_div_ps(vone,_mm512_sqrt_ps(_mm512_mask_blend_ps(k3,matGRe[4][4],_mm512_mul_ps(vneg_one,matGRe[4][4]))));\
        matND[4] = _mm512_sub_ps(vzero, matD[4]);\
    }while(0);\
}


#define GET_G55(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
    do{\
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
    }while(0);\
}

// get G55
#define GET_G55_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
    do{\
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
    }while(0);\
}

#define GET_G55_SAFE_MIN_APPROX_ERROR(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
    do{\
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
    }while(0);\
}


#define GET_G66(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
    do{\
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
    }while(0);\
}

// get G66
#define GET_G66_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
    do{\
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
    }while(0);\
}

#define GET_G66_SAFE_MIN_APPROX_ERROR(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
    do{\
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
    }while(0);\
}

#define GET_G77(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
    do{\
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
    }while(0);\
}

// get G77
#define GET_G77_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
    do{\
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
    }while(0);\
}

#define GET_G77_SAFE_MIN_APPROX_ERROR(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
    do{\
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
    }while(0);\
}


#define GET_Gii_ODD(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, i)\
{\
    do{\
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
    }while(0);\
}

// get Gii, odd diagonal element
#define GET_Gii_ODD_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, i)\
{\
    do{\
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
    }while(0);\
}

#define GET_Gii_ODD_SAFE_MIN_APPROX_ERROR(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, i)\
{\
    do{\
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
    }while(0);\
}

#define GET_Gii_EVEN(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, i)\
{\
    do{\
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
    }while(0);\
}

// get Gii, even diagonal element
#define GET_Gii_EVEN_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, i)\
{\
    do{\
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
    }while(0);\
}

#define GET_Gii_EVEN_SAFE_MIN_APPROX_ERROR(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, i)\
{\
    do{\
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
    }while(0);\
}


//get odd column n of matrix G, j is row index, i is col index
#define GET_G_COL_ODD(matGRe, matGIm, matBRe, matBIm, matD, j, i, temp0, temp1)\
{\
    do{\
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
    }while(0);\
}

// get even column n of matrix G, j is row index, i is col index
#define GET_G_COL_EVEN(matGRe, matGIm, matBRe, matBIm, matD, j, i, temp0, temp1)\
{\
    do{\
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
    }while(0);\
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
    do{\
        GET_AxRealB(matGRe[j][i], matGIm[j][i], matLRe[i][i],  matLRe[j][i], matLIm[j][i]);\
        for (std::int32_t i1 = i+1; i1 < j; i1++)\
        {\
            GET_AxB(matGRe[j][i1], matGIm[j][i1], matLRe[i1][i], matLIm[i1][i], temp0, temp1);\
            matLRe[j][i] = _mm512_add_ps(matLRe[j][i], temp0);\
            matLIm[j][i] = _mm512_add_ps(matLIm[j][i], temp1);\
        }\
        matLRe[j][i] = _mm512_mul_ps(matLRe[j][i], matND[j]);\
        matLIm[j][i] = _mm512_mul_ps(matLIm[j][i], matND[j]);\
    }while(0);\
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
