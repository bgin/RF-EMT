
/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1989 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

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

#ifndef __GMS_VEC_CEPHES_COS_PS_H__
#define __GMS_VEC_CEPHES_COS_PS_H__ 061220250611


namespace file_info 
{

     static const unsigned int GMS_VEC_CEPHES_COS_PS_MAJOR = 1;
     static const unsigned int GMS_VEC_CEPHES_COS_PS_MINOR = 1;
     static const unsigned int GMS_VEC_CEPHES_COS_PS_MICRO = 0;
     static const unsigned int GMS_VEC_CEPHES_COS_PS_FULLVER =
       1000U*GMS_VEC_CEPHES_COS_PS_MAJOR+100U*GMS_VEC_CEPHES_COS_PS_MINOR+
       10U*GMS_VEC_CEPHES_COS_PS_MICRO;
     static const char GMS_VEC_CEPHES_COS_PS_CREATION_DATE[] = "06-12-2025 06:12AM +00200 (SAT 06 DEC 2025 GMT+2)";
     static const char GMS_VEC_CEPHES_COS_PS_BUILD_DATE[]    = __DATE__; 
     static const char GMS_VEC_CEPHES_COS_PS_BUILD_TIME[]    = __TIME__;
     static const char GMS_VEC_CEPHES_COS_PS_SYNOPSIS[]      = "Manual vectorization of CEPHES cosine(float) function.";

}

#include <immintrin.h>
#include <limits>
#include "GMS_config.h"

namespace gms 
{

namespace math 
{

namespace {

__ATTR_ALWAYS_INLINE__		    
static inline 
__m128
_mm_negate_ps(const __m128 v) 
{
    const __m128 nvzero{_mm_set1_ps(-0.0f)};
	return (_mm_xor_ps(v,nvzero));
}

__ATTR_ALWAYS_INLINE__		    
static inline 
__m256
_mm256_negate_ps(const __m256 v) 
{
    const __m256 nvzero{_mm256_set1_ps(-0.0f)};
	return (_mm256_xor_ps(v,nvzero));
}

__ATTR_ALWAYS_INLINE__		    
static inline 
__m512
_mm512_negate_ps(const __m512 v) 
{
    const __m512 nvzero{_mm512_set1_ps(-0.0f)};
	return (_mm512_xor_ps(v,nvzero));
}

}

__ATTR_ALWAYS_INLINE__
static inline 
__m128 _mm_ceph_cosf_ps(const __m128 xx) 
{
	
	    /* These are for a 24-bit significand: 
           constexpr float PIO4F =  0.7853981633974483096f;
           constexpr float FOPI  = 1.27323954473516f;
           constexpr float DP1 = 0.78515625f;
           constexpr float DP2 = 2.4187564849853515625e-4f;
           constexpr float DP3 = 3.77489497744594108e-8f;
           constexpr float lossth = 8192.f;
           constexpr float T24M1 = 16777215.f;
	     */
	const __m128 PIO4F{ _mm_set1_ps(0.7853981633974483096f)};
	const __m128 FOPI{  _mm_set1_ps(1.27323954473516f)};
	const __m128 DP1{   _mm_set1_ps(0.78515625f)};
	const __m128 DP2{   _mm_set1_ps(2.4187564849853515625e-4f)};
	const __m128 DP3{   _mm_set1_ps(3.77489497744594108e-8f)};
	const __m128 lossth{_mm_set1_ps(8192.0f)};
	const __m128 T24M1{ _mm_set1_ps(16777215.0f)};
	const __m128 C19515295891E4{_mm_set1_ps(-1.9515295891E-4)};
	const __m128 C83321608736E3{_mm_set1_ps(8.3321608736E-3)};
	const __m128 C16666654611E1{_mm_set1_ps(1.6666654611E-1)};
	const __m128 C2443315711809948E005{_mm_set1_ps(2.443315711809948E-005)};
	const __m128 C1388731625493765E003{_mm_set1_ps(1.388731625493765E-003)};
	const __m128 C4166664568298827E002{_mm_set1_ps(4.166664568298827E-002)};
	const __m128 C05{                  _mm_set1_ps(0.5f)};
	const __m128 zero{  _mm_setzero_ps()};
	const __m128 infinity{_mm_set1_ps(std::numeric_limits<float>::infinity())};
	const __m128i i_one{_mm_set1_epi32(1)};
    const __m128  f_one{_mm_set1_ps(1.0f)};
	const __m128i  i_7{  _mm_set1_epi32(7)};
	const __m128i  i_3{  _mm_set1_epi32(3)};
	const __m128i  i_4{  _mm_set1_epi32(4)};
	const __m128i  i_0{  _mm_set1_epi32(0)};
	const __m128i  i_2{  _mm_set1_epi32(2)};
    __m128  x;
	__m128  y;
	__m128  z;
	__m128  x_true;
	__m128  x_false;
	__m128  y_true;
	__m128  y_false;
	__m128i  neg_sign;
	__m128i j;
	__m128i j_and_1;
    __m128i sign;
	__mmask8 x_lt_0{};
	__mmask8 x_gt_T24M1{};
	__mmask8 mask_j_and_1{};
	__mmask8 j_gt_3{};
	__mmask8 j_gt_1{};
	__mmask8 x_gt_lossth{};
	__mmask8 j_eq_1{};
	__mmask8 j_eq_2{};
	__mmask8 sign_lt_0;
	sign       = _mm_set1_epi32(1);
	x          = xx;
	neg_sign   = _mm_sub_epi32(i_0,sign);
    x_lt_0     = _mm_cmp_ps_mask(x,zero,_CMP_LT_OQ);
	x          = _mm_mask_blend_ps(x_lt_0,x,_mm_negate_ps(x));
    x_gt_T24M1 = _mm_cmp_ps_mask(x,T24M1,_CMP_GT_OQ);
	if(__builtin_expect(x_gt_T24M1==0xFF,0)) {return infinity;}
    j          = _mm_cvtps_epi32(_mm_mul_ps(FOPI,x));
    y          = _mm_cvtepi32_ps(j);
    j_and_1    = _mm_and_si128(j,i_one);
	mask_j_and_1= _mm_cmp_epi32_mask(j_and_1,i_0,_MM_CMPINT_NE);
    j          = _mm_mask_blend_epi32(mask_j_and_1,j,_mm_add_epi32(j,i_one));
	y          = _mm_mask_blend_ps(mask_j_and_1,y,_mm_add_ps(y,f_one));
    j          = _mm_and_si128(j,i_7);
    j_gt_3     = _mm_cmp_epi32_mask(j,i_3,_MM_CMPINT_GT);
	j          = _mm_mask_blend_epi32(j_gt_3,j,_mm_sub_epi32(j,i_4));
	sign       = _mm_mask_blend_epi32(j_gt_3,sign,neg_sign);
    j_gt_1     = _mm_cmp_epi32_mask(j,i_one,_MM_CMPINT_GT);
    sign       = _mm_mask_blend_epi32(j_gt_1,sign,neg_sign);
	x_gt_lossth= _mm_cmp_ps_mask(x,lossth,_CMP_GT_OQ);
	x_true     = _mm_sub_ps(x,_mm_mul_ps(y,PIO4F));
	//x_false  = ((x-y*DP1)-y*DP2)-y*DP3;
	x_false    = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(x,_mm_mul_ps(y,DP1)),_mm_mul_ps(y,DP2)),_mm_mul_ps(y,DP3));
	x          = _mm_mask_blend_ps(x_gt_lossth,x_false,x_true);
    z          = _mm_mul_ps(x,x);
    j_eq_1     = _mm_cmp_epi32_mask(j,i_one,_MM_CMPINT_EQ);
	j_eq_2     = _mm_cmp_epi32_mask(j,i_2,_MM_CMPINT_EQ);
	const __mmask8 j_eq1_or_j_eq_2{_kor_mask8(j_eq_1,j_eq_2)};
    y_true     = _mm_fmadd_ps(_mm_mul_ps(_mm_fmsub_ps(_mm_fmadd_ps(C19515295891E4,z,C83321608736E3),z,C16666654611E1),z),x,x);
	y_false    = _mm_mul_ps(_mm_mul_ps(_mm_fmadd_ps(_mm_fmsub_ps(C2443315711809948E005,z,C1388731625493765E003),z,C4166664568298827E002),z),z);
	y_false    = _mm_sub_ps(y_false,_mm_mul_ps(C05,z));
	y_false    = _mm_add_ps(y_false,f_one);
	y          = _mm_mask_blend_ps(j_eq1_or_j_eq_2,y_false,y_true);
    sign_lt_0  = _mm_cmp_epi32_mask(sign,i_0,_MM_CMPINT_LT);
	y          = _mm_mask_blend_ps(sign_lt_0,y,_mm_negate_ps(y));
	return (y);
}

__ATTR_ALWAYS_INLINE__
static inline 
__m256 _mm256_ceph_cosf_ps(const __m256 xx) 
{
	
	    /* These are for a 24-bit significand: 
           constexpr float PIO4F =  0.7853981633974483096f;
           constexpr float FOPI  = 1.27323954473516f;
           constexpr float DP1 = 0.78515625f;
           constexpr float DP2 = 2.4187564849853515625e-4f;
           constexpr float DP3 = 3.77489497744594108e-8f;
           constexpr float lossth = 8192.f;
           constexpr float T24M1 = 16777215.f;
	     */
	const __m256 PIO4F{ _mm256_set1_ps(0.7853981633974483096f)};
	const __m256 FOPI{  _mm256_set1_ps(1.27323954473516f)};
	const __m256 DP1{   _mm256_set1_ps(0.78515625f)};
	const __m256 DP2{   _mm256_set1_ps(2.4187564849853515625e-4f)};
	const __m256 DP3{   _mm256_set1_ps(3.77489497744594108e-8f)};
	const __m256 lossth{_mm256_set1_ps(8192.0f)};
	const __m256 T24M1{ _mm256_set1_ps(16777215.0f)};
	const __m256 C19515295891E4{_mm256_set1_ps(-1.9515295891E-4)};
	const __m256 C83321608736E3{_mm256_set1_ps(8.3321608736E-3)};
	const __m256 C16666654611E1{_mm256_set1_ps(1.6666654611E-1)};
	const __m256 C2443315711809948E005{_mm256_set1_ps(2.443315711809948E-005)};
	const __m256 C1388731625493765E003{_mm256_set1_ps(1.388731625493765E-003)};
	const __m256 C4166664568298827E002{_mm256_set1_ps(4.166664568298827E-002)};
	const __m256 C05{                  _mm256_set1_ps(0.5f)};
	const __m256 zero{    _mm256_setzero_ps()};
	const __m256 infinity{_mm256_set1_ps(std::numeric_limits<float>::infinity())};
	const __m256i i_one{  _mm256_set1_epi32(1)};
    const __m256  f_one{  _mm256_set1_ps(1.0f)};
	const __m256i  i_7{   _mm256_set1_epi32(7)};
	const __m256i  i_3{   _mm256_set1_epi32(3)};
	const __m256i  i_4{   _mm256_set1_epi32(4)};
	const __m256i  i_0{   _mm256_set1_epi32(0)};
	const __m256i  i_2{   _mm256_set1_epi32(2)};
    __m256  x;
	__m256  y;
	__m256  z;
	__m256  x_true;
	__m256  x_false;
	__m256  y_true;
	__m256  y_false;
	__m256i  neg_sign;
	__m256i j;
	__m256i j_and_1;
    __m256i sign;
	__mmask8 x_lt_0{};
	__mmask8 x_gt_T24M1{};
	__mmask8 mask_j_and_1{};
	__mmask8 j_gt_3{};
	__mmask8 j_gt_1{};
	__mmask8 x_gt_lossth{};
	__mmask8 j_eq_1{};
	__mmask8 j_eq_2{};
	__mmask8 sign_lt_0;
	sign       = _mm256_set1_epi32(1);
	x          = xx;
	neg_sign   = _mm256_sub_epi32(i_0,sign);
    x_lt_0     = _mm256_cmp_ps_mask(x,zero,_CMP_LT_OQ);
	x          = _mm256_mask_blend_ps(x_lt_0,x,_mm256_negate_ps(x));
    x_gt_T24M1 = _mm256_cmp_ps_mask(x,T24M1,_CMP_GT_OQ);
	if(__builtin_expect(x_gt_T24M1==0xFF,0)) {return infinity;}
    j          = _mm256_cvtps_epi32(_mm256_mul_ps(FOPI,x));
    y          = _mm256_cvtepi32_ps(j);
    j_and_1    = _mm256_and_si256(j,i_one);
	mask_j_and_1= _mm256_cmp_epi32_mask(j_and_1,i_0,_MM_CMPINT_NE);
    j          = _mm256_mask_blend_epi32(mask_j_and_1,j,_mm256_add_epi32(j,i_one));
	y          = _mm256_mask_blend_ps(mask_j_and_1,y,_mm256_add_ps(y,f_one));
    j          = _mm256_and_si256(j,i_7);
    j_gt_3     = _mm256_cmp_epi32_mask(j,i_3,_MM_CMPINT_GT);
	j          = _mm256_mask_blend_epi32(j_gt_3,j,_mm256_sub_epi32(j,i_4));
	sign       = _mm256_mask_blend_epi32(j_gt_3,sign,neg_sign);
    j_gt_1     = _mm256_cmp_epi32_mask(j,i_one,_MM_CMPINT_GT);
    sign       = _mm256_mask_blend_epi32(j_gt_1,sign,neg_sign);
	x_gt_lossth= _mm256_cmp_ps_mask(x,lossth,_CMP_GT_OQ);
	x_true     = _mm256_sub_ps(x,_mm256_mul_ps(y,PIO4F));
	//x_false  = ((x-y*DP1)-y*DP2)-y*DP3;
	x_false    = _mm256_sub_ps(_mm256_sub_ps(_mm256_sub_ps(x,_mm256_mul_ps(y,DP1)),_mm256_mul_ps(y,DP2)),_mm256_mul_ps(y,DP3));
	x          = _mm256_mask_blend_ps(x_gt_lossth,x_false,x_true);
    z          = _mm256_mul_ps(x,x);
    j_eq_1     = _mm256_cmp_epi32_mask(j,i_one,_MM_CMPINT_EQ);
	j_eq_2     = _mm256_cmp_epi32_mask(j,i_2,_MM_CMPINT_EQ);
	const __mmask8 j_eq1_or_j_eq_2{_kor_mask8(j_eq_1,j_eq_2)};
    y_true     = _mm256_fmadd_ps(_mm256_mul_ps(_mm256_fmsub_ps(_mm256_fmadd_ps(C19515295891E4,z,C83321608736E3),z,C16666654611E1),z),x,x);
	y_false    = _mm256_mul_ps(_mm256_mul_ps(_mm256_fmadd_ps(_mm256_fmsub_ps(C2443315711809948E005,z,C1388731625493765E003),z,C4166664568298827E002),z),z);
	y_false    = _mm256_sub_ps(y_false,_mm256_mul_ps(C05,z));
	y_false    = _mm256_add_ps(y_false,f_one);
	y          = _mm256_mask_blend_ps(j_eq1_or_j_eq_2,y_false,y_true);
    sign_lt_0  = _mm256_cmp_epi32_mask(sign,i_0,_MM_CMPINT_LT);
	y          = _mm256_mask_blend_ps(sign_lt_0,y,_mm256_negate_ps(y));
	return (y);
}

__ATTR_ALWAYS_INLINE__
static inline 
__m512 _mm512_ceph_cosf_ps(const __m512 xx) 
{
	
	    /* These are for a 24-bit significand: 
           constexpr float PIO4F =  0.7853981633974483096f;
           constexpr float FOPI  = 1.27323954473516f;
           constexpr float DP1 = 0.78515625f;
           constexpr float DP2 = 2.4187564849853515625e-4f;
           constexpr float DP3 = 3.77489497744594108e-8f;
           constexpr float lossth = 8192.f;
           constexpr float T24M1 = 16777215.f;
	     */
	const __m512 PIO4F{ _mm512_set1_ps(0.7853981633974483096f)};
	const __m512 FOPI{  _mm512_set1_ps(1.27323954473516f)};
	const __m512 DP1{   _mm512_set1_ps(0.78515625f)};
	const __m512 DP2{   _mm512_set1_ps(2.4187564849853515625e-4f)};
	const __m512 DP3{   _mm512_set1_ps(3.77489497744594108e-8f)};
	const __m512 lossth{_mm512_set1_ps(8192.0f)};
	const __m512 T24M1{ _mm512_set1_ps(16777215.0f)};
	const __m512 C19515295891E4{_mm512_set1_ps(-1.9515295891E-4)};
	const __m512 C83321608736E3{_mm512_set1_ps(8.3321608736E-3)};
	const __m512 C16666654611E1{_mm512_set1_ps(1.6666654611E-1)};
	const __m512 C2443315711809948E005{_mm512_set1_ps(2.443315711809948E-005)};
	const __m512 C1388731625493765E003{_mm512_set1_ps(1.388731625493765E-003)};
	const __m512 C4166664568298827E002{_mm512_set1_ps(4.166664568298827E-002)};
	const __m512 C05{                  _mm512_set1_ps(0.5f)};
	const __m512 zero{  _mm512_setzero_ps()};
	const __m512 infinity{_mm512_set1_ps(std::numeric_limits<float>::infinity())};
	const __m512i i_one{_mm512_set1_epi32(1)};
    const __m512  f_one{_mm512_set1_ps(1.0f)};
	const __m512i  i_7{  _mm512_set1_epi32(7)};
	const __m512i  i_3{  _mm512_set1_epi32(3)};
	const __m512i  i_4{  _mm512_set1_epi32(4)};
	const __m512i  i_0{  _mm512_set1_epi32(0)};
	const __m512i  i_2{  _mm512_set1_epi32(2)};
    const __m512i  neg_sign{_mm512_set1_epi32(-1)};
    __m512  x;
	__m512  y;
	__m512  z;
	__m512  x_true;
	__m512  x_false;
	__m512  y_true;
	__m512  y_false;
	//__m512i  neg_sign;
	__m512i j;
	__m512i j_and_1;
    __m512i sign;
	__mmask16 x_lt_0{};
	__mmask16 x_gt_T24M1{};
	__mmask16 mask_j_and_1{};
	__mmask16 j_gt_3{};
	__mmask16 j_gt_1{};
	__mmask16 x_gt_lossth{};
	__mmask16 j_eq_1{};
	__mmask16 j_eq_2{};
	__mmask16 sign_lt_0;
	sign       = _mm512_set1_epi32(1);
	x          = xx;
	x_lt_0     = _mm512_cmp_ps_mask(x,zero,_CMP_LT_OQ);
	x          = _mm512_mask_blend_ps(x_lt_0,x,_mm512_negate_ps(x));
    x_gt_T24M1 = _mm512_cmp_ps_mask(x,T24M1,_CMP_GT_OQ);
	if(__builtin_expect(x_gt_T24M1==0xFFFF,0)) {return infinity;}
    j          = _mm512_cvtps_epi32(_mm512_mul_ps(FOPI,x));
    y          = _mm512_cvtepi32_ps(j);
    j_and_1    = _mm512_and_si512(j,i_one);
	mask_j_and_1= _mm512_cmp_epi32_mask(j_and_1,i_one,_MM_CMPINT_EQ); // was (j_and_1,i_0,_MM_CMP_INT_NE)
    j          = _mm512_mask_blend_epi32(mask_j_and_1,j,_mm512_add_epi32(j,i_one));
	y          = _mm512_mask_blend_ps(mask_j_and_1,y,_mm512_add_ps(y,f_one));
    j          = _mm512_and_si512(j,i_7);
    j_gt_3     = _mm512_cmp_epi32_mask(j,i_3,_MM_CMPINT_GT);
	j          = _mm512_mask_blend_epi32(j_gt_3,j,_mm512_sub_epi32(j,i_4));
	sign       = _mm512_mask_blend_epi32(j_gt_3,sign,neg_sign);
    j_gt_1     = _mm512_cmp_epi32_mask(j,i_one,_MM_CMPINT_GT);
    sign       = _mm512_mask_blend_epi32(j_gt_1,sign,neg_sign);
	x_gt_lossth= _mm512_cmp_ps_mask(x,lossth,_CMP_GT_OQ);
	x_true     = _mm512_sub_ps(x,_mm512_mul_ps(y,PIO4F));
	//x_false  = ((x-y*DP1)-y*DP2)-y*DP3;
	x_false    = _mm512_sub_ps(_mm512_sub_ps(_mm512_sub_ps(x,_mm512_mul_ps(y,DP1)),_mm512_mul_ps(y,DP2)),_mm512_mul_ps(y,DP3));
	x          = _mm512_mask_blend_ps(x_gt_lossth,x_false,x_true);
    z          = _mm512_mul_ps(x,x);
    j_eq_1     = _mm512_cmp_epi32_mask(j,i_one,_MM_CMPINT_EQ);
	j_eq_2     = _mm512_cmp_epi32_mask(j,i_2,_MM_CMPINT_EQ);
	const __mmask16 j_eq_1_or_j_eq_2{_kor_mask16(j_eq_1,j_eq_2)};
    y_true     = _mm512_fmadd_ps(_mm512_mul_ps(_mm512_fmsub_ps(_mm512_fmadd_ps(C19515295891E4,z,C83321608736E3),z,C16666654611E1),z),x,x);
	y_false    = _mm512_mul_ps(_mm512_mul_ps(_mm512_fmadd_ps(_mm512_fmsub_ps(C2443315711809948E005,z,C1388731625493765E003),z,C4166664568298827E002),z),z);
	y_false    = _mm512_sub_ps(y_false,_mm512_mul_ps(C05,z));
	y_false    = _mm512_add_ps(y_false,f_one);
	y          = _mm512_mask_blend_ps(j_eq_1_or_j_eq_2,y_false,y_true);
    sign_lt_0  = _mm512_cmp_epi32_mask(sign,i_0,_MM_CMPINT_LT);
	y          = _mm512_mask_blend_ps(sign_lt_0,y,_mm512_negate_ps(y));
	return (y);
}


} //math

} //gms


#endif /*__GMS_VEC_CEPHES_COS_PS_H__*/