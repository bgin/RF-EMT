

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

#ifndef __UNIT_TEST_MAT_INV_CHOLESKY_DEFINES_H__
#define __UNIT_TEST_MAT_INV_CHOLESKY_DEFINES_H__ 260620260818


//static __m512 constZero = _mm512_setzero_ps();

// complex multiplication: A * A'
#define GET_AxAH_REF(re, im, out)\
{\
    out = _mm512_add_ps(_mm512_mul_ps(re, re), _mm512_mul_ps(im, im));\
}

// complex multiplication: A * B'
#define GET_AxBH_REF(are, aim, bre, bim, outre, outim)\
{\
    outre = _mm512_add_ps(_mm512_mul_ps(are, bre), _mm512_mul_ps(aim, bim));\
    outim = _mm512_sub_ps(_mm512_mul_ps(bre, aim), _mm512_mul_ps(are, bim));\
}

// complex multiplication: A * B
#define GET_AxB_REF(are, aim, bre, bim, outre, outim)\
{\
    outre = _mm512_sub_ps(_mm512_mul_ps(are, bre), _mm512_mul_ps(aim, bim));\
    outim = _mm512_add_ps(_mm512_mul_ps(bre, aim), _mm512_mul_ps(are, bim));\
}

// complex multiplication: A * real(B)
#define GET_AxRealB_REF(are, aim, bre, outre, outim)\
{\
    outre = _mm512_mul_ps(are, bre);\
    outim = _mm512_mul_ps(aim, bre);\
}

#define GET_G00_REF(matGRe, matBRe, matD, matND)\
{\
    matD[0] = _mm512_rsqrt14_ps(matBRe[0][0]);\
    matND[0] = _mm512_sub_ps(constZero, matD[0]);\
}

// get column 0 of matrix G
#define GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, j)\
{\
    matGRe[j][0] = _mm512_mul_ps(matBRe[j][0], matD[0]);\
    matGIm[j][0] = _mm512_mul_ps(matBIm[j][0], matD[0]);\
}

// get G11
#define GET_G11_REF(matGRe, matGIm, matBRe, matD, matND, temp0)\
{\
    GET_AxAH_REF(matGRe[1][0], matGIm[1][0], temp0);\
    matGRe[1][1] = _mm512_sub_ps(matBRe[1][1], temp0);\
    matD[1] = _mm512_rsqrt14_ps(matGRe[1][1]);\
    matND[1] = _mm512_sub_ps(constZero, matD[1]);\
}

// get column 1 of matrix G
#define GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, j, temp0, temp1)\
{\
    GET_AxBH_REF(matGRe[j][0], matGIm[j][0], matGRe[1][0], matGIm[1][0], temp0, temp1);\
    matGRe[j][1] = _mm512_sub_ps(matBRe[j][1], temp0);\
    matGIm[j][1] = _mm512_sub_ps(matBIm[j][1], temp1);\
    matGRe[j][1] = _mm512_mul_ps(matGRe[j][1], matD[1]);\
    matGIm[j][1] = _mm512_mul_ps(matGIm[j][1], matD[1]);\
}

// get G22
#define GET_G22_REF(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
    GET_AxAH_REF(matGRe[2][0], matGIm[2][0], temp0);\
    GET_AxAH_REF(matGRe[2][1], matGIm[2][1], temp1);\
    matGRe[2][2] = _mm512_sub_ps(matBRe[2][2], temp0);\
    matGRe[2][2] = _mm512_sub_ps(matGRe[2][2], temp1);\
    matD[2] = _mm512_rsqrt14_ps(matGRe[2][2]);\
    matND[2] = _mm512_sub_ps(constZero, matD[2]);\
}

// get column 2 of matrix G
#define GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, j, temp0, temp1)\
{\
    GET_AxBH_REF(matGRe[j][0], matGIm[j][0], matGRe[2][0], matGIm[2][0], temp0, temp1);\
    matGRe[j][2] = _mm512_sub_ps(matBRe[j][2], temp0);\
    matGIm[j][2] = _mm512_sub_ps(matBIm[j][2], temp1);\
    GET_AxBH_REF(matGRe[j][1], matGIm[j][1], matGRe[2][1], matGIm[2][1], temp0, temp1);\
    matGRe[j][2] = _mm512_sub_ps(matGRe[j][2], temp0);\
    matGIm[j][2] = _mm512_sub_ps(matGIm[j][2], temp1);\
    matGRe[j][2] = _mm512_mul_ps(matGRe[j][2], matD[2]);\
    matGIm[j][2] = _mm512_mul_ps(matGIm[j][2], matD[2]);\
}

// get G33
#define GET_G33_REF(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
    GET_AxAH_REF(matGRe[3][0], matGIm[3][0], temp0);\
    GET_AxAH_REF(matGRe[3][1], matGIm[3][1], temp1);\
    GET_AxAH_REF(matGRe[3][2], matGIm[3][2], temp2);\
    matGRe[3][3] = _mm512_sub_ps(matBRe[3][3], temp0);\
    matGRe[3][3] = _mm512_sub_ps(matGRe[3][3], temp1);\
    matGRe[3][3] = _mm512_sub_ps(matGRe[3][3], temp2);\
    matD[3] = _mm512_rsqrt14_ps(matGRe[3][3]);\
    matND[3] = _mm512_sub_ps(constZero, matD[3]);\
}

// get column 3 of matrix G
#define GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, j, temp0, temp1)\
{\
    GET_AxBH_REF(matGRe[j][0], matGIm[j][0], matGRe[3][0], matGIm[3][0], temp0, temp1);\
    matGRe[j][3] = _mm512_sub_ps(matBRe[j][3], temp0);\
    matGIm[j][3] = _mm512_sub_ps(matBIm[j][3], temp1);\
    GET_AxBH_REF(matGRe[j][1], matGIm[j][1], matGRe[3][1], matGIm[3][1], temp0, temp1);\
    matGRe[j][3] = _mm512_sub_ps(matGRe[j][3], temp0);\
    matGIm[j][3] = _mm512_sub_ps(matGIm[j][3], temp1);\
    GET_AxBH_REF(matGRe[j][2], matGIm[j][2], matGRe[3][2], matGIm[3][2], temp0, temp1);\
    matGRe[j][3] = _mm512_sub_ps(matGRe[j][3], temp0);\
    matGIm[j][3] = _mm512_sub_ps(matGIm[j][3], temp1);\
    matGRe[j][3] = _mm512_mul_ps(matGRe[j][3], matD[3]);\
    matGIm[j][3] = _mm512_mul_ps(matGIm[j][3], matD[3]);\
}

// get G44
#define GET_G44_REF(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
    GET_AxAH_REF(matGRe[4][0], matGIm[4][0], temp0);\
    GET_AxAH_REF(matGRe[4][1], matGIm[4][1], temp1);\
    matGRe[4][4] = _mm512_sub_ps(matBRe[4][4], temp0);\
    matGRe[4][4] = _mm512_sub_ps(matGRe[4][4], temp1);\
    GET_AxAH_REF(matGRe[4][2], matGIm[4][2], temp0);\
    GET_AxAH_REF(matGRe[4][3], matGIm[4][3], temp1);\
    matGRe[4][4] = _mm512_sub_ps(matGRe[4][4], temp0);\
    matGRe[4][4] = _mm512_sub_ps(matGRe[4][4], temp1);\
    matD[4] = _mm512_rsqrt14_ps(matGRe[4][4]);\
    matND[4] = _mm512_sub_ps(constZero, matD[4]);\
}

// get G55
#define GET_G55_REF(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
    GET_AxAH_REF(matGRe[5][0], matGIm[5][0], temp0);\
    GET_AxAH_REF(matGRe[5][1], matGIm[5][1], temp1);\
    matGRe[5][5] = _mm512_sub_ps(matBRe[5][5], temp0);\
    matGRe[5][5] = _mm512_sub_ps(matGRe[5][5], temp1);\
    GET_AxAH_REF(matGRe[5][2], matGIm[5][2], temp0);\
    GET_AxAH_REF(matGRe[5][3], matGIm[5][3], temp1);\
    GET_AxAH_REF(matGRe[5][4], matGIm[5][4], temp2);\
    matGRe[5][5] = _mm512_sub_ps(matGRe[5][5], temp0);\
    matGRe[5][5] = _mm512_sub_ps(matGRe[5][5], temp1);\
    matGRe[5][5] = _mm512_sub_ps(matGRe[5][5], temp2);\
    matD[5] = _mm512_rsqrt14_ps(matGRe[5][5]);\
    matND[5] = _mm512_sub_ps(constZero, matD[5]);\
}

// get G66
#define GET_G66_REF(matGRe, matGIm, matBRe, matD, matND, temp0, temp1)\
{\
    GET_AxAH_REF(matGRe[6][0], matGIm[6][0], temp0);\
    GET_AxAH_REF(matGRe[6][1], matGIm[6][1], temp1);\
    matGRe[6][6] = _mm512_sub_ps(matBRe[6][6], temp0);\
    matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp1);\
    GET_AxAH_REF(matGRe[6][2], matGIm[6][2], temp0);\
    GET_AxAH_REF(matGRe[6][3], matGIm[6][3], temp1);\
    matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp0);\
    matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp1);\
    GET_AxAH_REF(matGRe[6][4], matGIm[6][4], temp0);\
    GET_AxAH_REF(matGRe[6][5], matGIm[6][5], temp1);\
    matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp0);\
    matGRe[6][6] = _mm512_sub_ps(matGRe[6][6], temp1);\
    matD[6] = _mm512_rsqrt14_ps(matGRe[6][6]);\
    matND[6] = _mm512_sub_ps(constZero, matD[6]);\
}

// get G77
#define GET_G77_REF(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2)\
{\
    GET_AxAH_REF(matGRe[7][0], matGIm[7][0], temp0);\
    GET_AxAH_REF(matGRe[7][1], matGIm[7][1], temp1);\
    matGRe[7][7] = _mm512_sub_ps(matBRe[7][7], temp0);\
    matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp1);\
    GET_AxAH_REF(matGRe[7][2], matGIm[7][2], temp0);\
    GET_AxAH_REF(matGRe[7][3], matGIm[7][3], temp1);\
    matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp0);\
    matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp1);\
    GET_AxAH_REF(matGRe[7][4], matGIm[7][4], temp0);\
    GET_AxAH_REF(matGRe[7][5], matGIm[7][5], temp1);\
    GET_AxAH_REF(matGRe[7][6], matGIm[7][6], temp2);\
    matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp0);\
    matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp1);\
    matGRe[7][7] = _mm512_sub_ps(matGRe[7][7], temp2);\
    matD[7] = _mm512_rsqrt14_ps(matGRe[7][7]);\
    matND[7] = _mm512_sub_ps(constZero, matD[7]);\
}

// get Gii, odd diagonal element
#define GET_Gii_ODD_REF(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, i)\
{\
    GET_AxAH_REF(matGRe[i][0], matGIm[i][0], temp0);\
    matGRe[i][i] = _mm512_sub_ps(matBRe[i][i], temp0);\
    for (int32_t i1 = 1; i1 < i; i1+=2) \
    {\
        GET_AxAH_REF(matGRe[i][i1], matGIm[i][i1], temp0);\
        GET_AxAH_REF(matGRe[i][i1+1], matGIm[i][i1+1], temp1);\
        matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp0);\
        matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp1);\
    }\
    matD[i] = _mm512_rsqrt14_ps(matGRe[i][i]);\
    matND[i] = _mm512_sub_ps(constZero, matD[i]);\
}

// get Gii, even diagonal element
#define GET_Gii_EVEN_REF(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, i)\
{\
    GET_AxAH_REF(matGRe[i][0], matGIm[i][0], temp0);\
    GET_AxAH_REF(matGRe[i][1], matGIm[i][1], temp1);\
    matGRe[i][i] = _mm512_sub_ps(matBRe[i][i], temp0);\
    matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp1);\
    for (int32_t i1 = 2; i1 < i; i1+=2) \
    {\
        GET_AxAH_REF(matGRe[i][i1], matGIm[i][i1], temp0);\
        GET_AxAH_REF(matGRe[i][i1+1], matGIm[i][i1+1], temp1);\
        matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp0);\
        matGRe[i][i] = _mm512_sub_ps(matGRe[i][i], temp1);\
    }\
    matD[i] = _mm512_rsqrt14_ps(matGRe[i][i]);\
    matND[i] = _mm512_sub_ps(constZero, matD[i]);\
}

// get odd column n of matrix G, j is row index, i is col index
#define GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, j, i, temp0, temp1)\
{\
    GET_AxBH_REF(matGRe[j][0], matGIm[j][0], matGRe[i][0], matGIm[i][0], temp0, temp1);\
    matGRe[j][i] = _mm512_sub_ps(matBRe[j][i], temp0);\
    matGIm[j][i] = _mm512_sub_ps(matBIm[j][i], temp1);\
    for (int32_t i1 = 1; i1 < i; i1+=2)\
    {\
        GET_AxBH_REF(matGRe[j][i1], matGIm[j][i1], matGRe[i][i1], matGIm[i][i1], temp0, temp1);\
        matGRe[j][i] = _mm512_sub_ps(matGRe[j][i], temp0);\
        matGIm[j][i] = _mm512_sub_ps(matGIm[j][i], temp1);\
        GET_AxBH_REF(matGRe[j][i1+1], matGIm[j][i1+1], matGRe[i][i1+1], matGIm[i][i1+1], temp0, temp1);\
        matGRe[j][i] = _mm512_sub_ps(matGRe[j][i], temp0);\
        matGIm[j][i] = _mm512_sub_ps(matGIm[j][i], temp1);\
    }\
    matGRe[j][i] = _mm512_mul_ps(matGRe[j][i], matD[i]);\
    matGIm[j][i] = _mm512_mul_ps(matGIm[j][i], matD[i]);\
}

// get even column n of matrix G, j is row index, i is col index
#define GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, j, i, temp0, temp1)\
{\
    GET_AxBH_REF(matGRe[j][0], matGIm[j][0], matGRe[i][0], matGIm[i][0], temp0, temp1);\
    matGRe[j][i] = _mm512_sub_ps(matBRe[j][i], temp0);\
    matGIm[j][i] = _mm512_sub_ps(matBIm[j][i], temp1);\
    GET_AxBH_REF(matGRe[j][1], matGIm[j][1], matGRe[i][1], matGIm[i][1], temp0, temp1);\
    matGRe[j][i] = _mm512_sub_ps(matGRe[j][i], temp0);\
    matGIm[j][i] = _mm512_sub_ps(matGIm[j][i], temp1);\
    for (int32_t i1 = 2; i1 < i; i1+=2)\
    {\
        GET_AxBH_REF(matGRe[j][i1], matGIm[j][i1], matGRe[i][i1], matGIm[i][i1], temp0, temp1);\
        matGRe[j][i] = _mm512_sub_ps(matGRe[j][i], temp0);\
        matGIm[j][i] = _mm512_sub_ps(matGIm[j][i], temp1);\
        GET_AxBH_REF(matGRe[j][i1+1], matGIm[j][i1+1], matGRe[i][i1+1], matGIm[i][i1+1], temp0, temp1);\
        matGRe[j][i] = _mm512_sub_ps(matGRe[j][i], temp0);\
        matGIm[j][i] = _mm512_sub_ps(matGIm[j][i], temp1);\
    }\
    matGRe[j][i] = _mm512_mul_ps(matGRe[j][i], matD[i]);\
    matGIm[j][i] = _mm512_mul_ps(matGIm[j][i], matD[i]);\
}

// set value for Lii, diagonal element
#define SET_Lii_REF(matLRe, matLIm, matD, i)\
{\
    matLRe[i][i] = matD[i];\
    matLIm[i][i] = constZero;\
}

// get element L(i+1, i)
#define GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, j, i)\
{\
    GET_AxRealB_REF(matGRe[j][i], matGIm[j][i], matLRe[i][i],  matLRe[j][i], matLIm[j][i]);\
    matLRe[j][i] = _mm512_mul_ps(matLRe[j][i], matND[j]);\
    matLIm[j][i] = _mm512_mul_ps(matLIm[j][i], matND[j]);\
}

#define GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, j, i, temp0, temp1)\
{\
    GET_AxRealB_REF(matGRe[j][i], matGIm[j][i], matLRe[i][i],  matLRe[j][i], matLIm[j][i]);\
    for (int32_t i1 = i+1; i1 < j; i1++)\
    {\
        GET_AxB_REF(matGRe[j][i1], matGIm[j][i1], matLRe[i1][i], matLIm[i1][i], temp0, temp1);\
        matLRe[j][i] = _mm512_add_ps(matLRe[j][i], temp0);\
        matLIm[j][i] = _mm512_add_ps(matLIm[j][i], temp1);\
    }\
    matLRe[j][i] = _mm512_mul_ps(matLRe[j][i], matND[j]);\
    matLIm[j][i] = _mm512_mul_ps(matLIm[j][i], matND[j]);\
}



#define N_2 2 
#define N_3 3
#define N_4 4
#define N_5 5
#define N_6 6
#define N_7 7
#define N_8 8
#define N_9 9
#define N_10 10
#define N_11 11
#define N_12 12
#define N_13 13
#define N_14 14
#define N_15 15
#define N_16 16









#endif /*__UNIT_TEST_MAT_INV_CHOLESKY_DEFINES_H__*/