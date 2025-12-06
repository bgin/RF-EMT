
#ifndef __GMS_CEPHES_SIN_COS_H__
#define __GMS_CEPHES_SIN_COS_H__ 061220250824

#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_simd_utils.h"

namespace gms
{
namespace math 
{

#if !defined (DENORMAL)
    #define DENORMAL 1
#endif

__ATTR_ALWAYS_INLINE__
static inline
float ceph_cosf( const float xx ) {
/* Note, these constants are for a 32-bit significand: */
/*
static float DP1 =  0.7853851318359375;
static float DP2 =  1.30315311253070831298828125e-5;
static float DP3 =  3.03855025325309630e-11;
static float lossth = 65536.;
*/

/* These are for a 24-bit significand: */
constexpr float PIO4F =  0.7853981633974483096f;
constexpr float FOPI  = 1.27323954473516f;
constexpr float DP1 = 0.78515625f;
constexpr float DP2 = 2.4187564849853515625e-4f;
constexpr float DP3 = 3.77489497744594108e-8f;
constexpr float lossth = 8192.f;
constexpr float T24M1 = 16777215.0f;
float x, y, z;
int j, sign;

/* make argument positive */
sign = 1;
x = xx;
if( x < 0 )
	x = -x;

if( x > T24M1 )
	{
	//mtherr( "cosf", TLOSS );
	return(std::numeric_limits<float>::infinity());
	}

j = FOPI * x; /* integer part of x/PIO4 */
y = j;
/* integer and fractional part modulo one octant */
if( j & 1 )	/* map zeros to origin */
	{
	j += 1;
	y += 1.0;
	}
j &= 7;
if( j > 3)
	{
	j -=4;
	sign = -sign;
	}

if( j > 1 )
	sign = -sign;

if( x > lossth )
	{
	//mtherr( "cosf", PLOSS );
	x = x - y * PIO4F;
	}
else
/* Extended precision modular arithmetic */
	x = ((x - y * DP1) - y * DP2) - y * DP3;

z = x * x;

if( (j==1) || (j==2) )
	{
	y = (((-1.9515295891E-4 * z
	     + 8.3321608736E-3) * z
	     - 1.6666654611E-1) * z * x)
	     + x;
	}
else
	{
	y = ((  2.443315711809948E-005 * z
	  - 1.388731625493765E-003) * z
	  + 4.166664568298827E-002) * z * z;
	y -= 0.5 * z;
	y += 1.0;
	}
if(sign < 0)
	y = -y;
return( y );
}

__ATTR_ALWAYS_INLINE__
static inline
float ceph_sinf( const float xx ) {
constexpr float FOPI  = 1.27323954473516;
constexpr float PIO4F =  0.7853981633974483096;
/* Note, these constants are for a 32-bit significand: */
/*
static float DP1 =  0.7853851318359375;
static float DP2 =  1.30315311253070831298828125e-5;
static float DP3 =  3.03855025325309630e-11;
static float lossth = 65536.;
*/

/* These are for a 24-bit significand: */
constexpr float DP1 = 0.78515625;
constexpr float DP2 = 2.4187564849853515625e-4;
constexpr float DP3 = 3.77489497744594108e-8;
constexpr float lossth = 8192.;
constexpr float T24M1 = 16777215.;

float sincof[] = {
-1.9515295891E-4,
 8.3321608736E-3,
-1.6666654611E-1
};
float coscof[] = {
 2.443315711809948E-005,
-1.388731625493765E-003,
 4.166664568298827E-002
};
float *p;
float x, y, z;
 unsigned long j;
 int sign;
sign = 1;
x = xx;
if( xx < 0 )
	{
	sign = -1;
	x = -xx;
	}
if( x > T24M1 )
	{
	//mtherr( "sinf", TLOSS );
	return(std::numeric_limits<float>::infinity());
	}
j = FOPI * x; /* integer part of x/(PI/4) */
y = j;
/* map zeros to origin */
if( j & 1 )
	{
	j += 1;
	y += 1.0;
	}
j &= 7; /* octant modulo 360 degrees */
/* reflect in x axis */
if( j > 3)
	{
	sign = -sign;
	j -= 4;
	}

if( x > lossth )
	{
	//mtherr( "sinf", PLOSS );
	x = x - y * PIO4F;
	}
else
	{
/* Extended precision modular arithmetic */
	x = ((x - y * DP1) - y * DP2) - y * DP3;
	}
/*einits();*/
z = x * x;
if( (j==1) || (j==2) )
	{
/* measured relative error in +/- pi/4 is 7.8e-8 */
/*
	y = ((  2.443315711809948E-005 * z
	  - 1.388731625493765E-003) * z
	  + 4.166664568298827E-002) * z * z;
*/
	p = coscof;
	y = *p++;
	y = y * z + *p++;
	y = y * z + *p++;
	y *= z * z;
	y -= 0.5 * z;
	y += 1.0;
	}
else
	{
/* Theoretical relative error = 3.8e-9 in [-pi/4, +pi/4] */
/*
	y = ((-1.9515295891E-4 * z
	     + 8.3321608736E-3) * z
	     - 1.6666654611E-1) * z * x;
	y += x;
*/
	p = sincof;
	y = *p++;
	y = y * z + *p++;
	y = y * z + *p++;
	y *= z * x;
	y += x;
	}
/*einitd();*/
if(sign < 0)
	y = -y;
return( y);
}

///////////// SIMD Kernels /////////////////

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
    __m128i sign;
	__mmask8 x_lt_0{};
	__mmask8 x_gt_T24M1{};
	__mmask8 j_and_1{};
	__mmask8 j_gt_3{};
	__mmask8 j_gt_1{};
	__mmask8 x_gt_lossth{};
	__mmask8 j_eq_1{};
	__mmask8 j_eq_2{};
	__mmask8 sign_lt_0;
	sign       = _mm_set1_epi32(1);
	x          = xx;
	neg_sign   = _mm_sub_ps(i_0,sign);
    x_lt_0     = _mm_cmp_ps_mask(x,zero,_CMP_LT_OQ);
	x          = _mm_mask_blend_ps(x_lt_0,x,negate_xmm4r4(x));
    x_gt_T24M1 = _mm_cmp_ps_mask(x,T24M1,_CMP_GT_OQ);
	if(__builtin_expect(x_gt_T24M1==0xFF,0)) {return infinity;}
    j          = _mm_castpd_si128(_mm_mul_ps(FOPI,x));
    y          = _mm_castsi128_ps(j);
    j_and_1    = _mm_and_si128(j,i_one);
    j          = _mm_mask_blend_epi32(j_and_1,j,_mm_add_epi32(j,i_one));
	y          = _mm_mask_blend_ps(j_and_1,y,_mm_add_ps(y,f_one));
    j          = _mm_and_epi32(j,i_7);
    j_gt_3     = _mm_cmp_epi32_mask(j,i_3,_CMP_GT_OQ);
	j          = _mm_mask_blend_epi32(j_gt_3,j,_mm_sub_epi32(j,i_4));
	sign       = _mm_mask_blend_epi32(j_gt_3,sign,neg_sign);
    j_gt_1     = _mm_cmp_epi32_mask(j,i_1,_CMP_GT_OQ);
    sign       = _mm_mask_blend_epi32(j_gt_1,sign,neg_sign);
	x_gt_lossth= _mm_cmp_ps_mask(x,lossth,_CMP_GT_OQ);
	x_true     = _mm_sub_ps(x,_mm_mul_ps(y,PIO4F));
	//x_false  = ((x-y*DP1)-y*DP2)-y*DP3;
	x_false    = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(x,_mm_mul_ps(y,DP1)),_mm_mul_ps(y,DP2)),_mm_mul_ps(y,DP3));
	x          = _mm_mask_blend_ps(x_gt_lossth,x_false,x_true);
    z          = _mm_mul_ps(x,x);
    j_eq_1     = _mm_cmp_epi32_mask(j,i_one,_CMP_EQ_OQ);
	j_eq_2     = _mm_cmp_epi32_mask(j,i_2,_CMP_EQ_OQ);
	const __mmask8 j_eq1_or_j_eq_2{_kor_mask8(j_eq_1,j_eq_2)};
    y_true     = _mm_fmadd_ps(_mm_mul_ps(_mm_fmsub_ps(_mm_fmadd_ps(C19515295891E4,z,C83321608736E3),z,C16666654611E1),z),x,x);
	y_false    = _mm_mul_ps(_mm_mul_ps(_mm_fmadd_ps(_mm_fmsub_ps(C2443315711809948E005,z,C1388731625493765E003),z,C4166664568298827E002),z),z);
	y_false    = _mm_sub_ps(y_false,_mm_mul_ps(C05,z));
	y_false    = _mm_add_ps(y_false,f_one);
	y          = _mm_mask_blend_ps(j_eq1_or_j_eq_2,y_false,y_true);
    sign_lt_0  = _mm_cmp_epi32_mask(sign,i_0,_CMP_LT_OQ);
	y          = _mm_mask_blend_ps(sign_lt_0,y,xmm4r4_abs(y));
	return (y);
}

__ATTR_ALWAYS_INLINE__
static inline
float ceph_frexpf( float x, int *pw2 ) {
union
  {
    float y;
    unsigned short i[2];
  } u;
int i, k;
unsigned short *q;
u.y = x;
q = &u.i[1];
/* find the exponent (power of 2) */

i  = ( *q >> 7) & 0xff;
if( i == 0 ) {
    if( u.y == 0.0 ) {
	*pw2 = 0;
	 return(0.0);
    }
/* Number is denormal or zero */
#if DENORMAL
/* Handle denormal number. */
	do{
		
	   u.y *= 2.0;
	   i -= 1;
	   k  = ( *q >> 7) & 0xff;
	}
	while( k == 0 );
	i = i + k;
#else
	*pw2 = 0;
	return( 0.0 );
#endif /* DENORMAL */
	}
i -= 0x7e;
*pw2 = i;
*q &= 0x807f;	/* strip all exponent bits */
*q |= 0x3f00;	/* mantissa between 0.5 and 1 */
return( u.y );
}


__ATTR_ALWAYS_INLINE__
static inline
float ceph_ldexpf( float x, int pw2 ) {
#define MEXP 255
union
  {
    float y;
    unsigned short i[2];
  } u;
unsigned short *q;
int e;
u.y = x;
q = &u.i[1];
while( (e = ( *q >> 7) & 0xff) == 0 ) {
	if( u.y == (float )0.0 ) {
		return( 0.0 );
	}
/* Input is denormal. */
	if( pw2 > 0 ) {
	    u.y *= 2.0;
	    pw2 -= 1;
	}
	if( pw2 < 0 ) {
	    if( pw2 < -24 )
		return( 0.0 );
	    u.y *= 0.5;
	    pw2 += 1;
	}
	if( pw2 == 0 )
		return(u.y);
	}
e += pw2;

/* Handle overflow */
if( e > MEXP )
	{
	return( 3.4028235e+38 );
	}

*q &= 0x807f;

/* Handle denormalized results */
if( e < 1 )
	{
#if DENORMAL
	if( e < -24 )
		return( 0.0 );
	*q |= 0x80; /* Set LSB of exponent. */
	/* For denormals, significant bits may be lost even
	   when dividing by 2.  Construct 2^-(1-e) so the result
	   is obtained with only one multiplication.  */
	u.y *= ceph_ldexpf(1.0f, e - 1);
	return(u.y);
#else
	return( 0.0 );
#endif
	}
*q |= (e & 0xff) << 7;
return(u.y);
}


__ATTR_ALWAYS_INLINE__
static inline
float ceph_sqrtf( const float xx ) {
float f, x, y;
int e;
f = xx;
if( f <= 0.0 )
	{
	if( f < 0.0 )
		//mtherr( "sqrtf", DOMAIN );
	return( 0.0 );
	}

x = ceph_frexpf( f, &e );	/* f = x * 2**e,   0.5 <= x < 1.0 */
/* If power of 2 is odd, double x and decrement the power of 2. */
if( e & 1 )
	{
	x = x + x;
	e -= 1;
	}

e >>= 1;	/* The power of 2 of the square root. */

if( x > 1.41421356237 )
	{
/* x is between sqrt(2) and 2. */
	x = x - 2.0;
	y =
	((((( -9.8843065718E-4 * x
	  + 7.9479950957E-4) * x
	  - 3.5890535377E-3) * x
	  + 1.1028809744E-2) * x
	  - 4.4195203560E-2) * x
	  + 3.5355338194E-1) * x
	  + 1.41421356237E0;
	goto sqdon;
	}

if( x > 0.707106781187 )
	{
/* x is between sqrt(2)/2 and sqrt(2). */
	x = x - 1.0;
	y =
	((((( 1.35199291026E-2 * x
	  - 2.26657767832E-2) * x
	  + 2.78720776889E-2) * x
	  - 3.89582788321E-2) * x
	  + 6.24811144548E-2) * x
	  - 1.25001503933E-1) * x * x
	  + 0.5 * x
	  + 1.0;
	goto sqdon;
	}

/* x is between 0.5 and sqrt(2)/2. */
x = x - 0.5;
y =
((((( -3.9495006054E-1 * x
  + 5.1743034569E-1) * x
  - 4.3214437330E-1) * x
  + 3.5310730460E-1) * x
  - 3.5354581892E-1) * x
  + 7.0710676017E-1) * x
  + 7.07106781187E-1;

sqdon:
y = ceph_ldexpf( y, e );  /* y = y * 2**e */
return( y);
}


__ATTR_ALWAYS_INLINE__
static inline 
float ceph_asinf(const float xx) {
float a, x, z;
int sign, flag;
x = xx;
if( x > 0 ) {
  sign = 1;
  a = x;
}
else
    {
  sign = -1;
  a = -x;
}

if( a > 1.0 ) {
    return( 0.0 );
}

if( a < 1.0e-4 ) {
    z = a;
    goto done;
}

if( a > 0.5 ) {
    z = 0.5 * (1.0 - a);
    x = ceph_sqrtf( z );
    flag = 1;
}
else
    {
    x = a;
    z = x * x;
    flag = 0;
}
z =
(((( 4.2163199048E-2 * z
  + 2.4181311049E-2) * z
  + 4.5470025998E-2) * z
  + 7.4953002686E-2) * z
  + 1.6666752422E-1) * z * x
  + x;

if( flag != 0 ) {
    z = z + z;
    z = 1.5707963267948966192 - z;
}
done:
if( sign < 0 )
	z = -z;
return( z );
}

__ATTR_ALWAYS_INLINE__
static inline 
float ceph_acosf(const float x ) {
if( x < -1.0f )
	goto domerr;

if( x < -0.5f) 
	return( 3.14159265358979323846f
                - 2.0f * ceph_asinf(ceph_sqrtf(0.5f*(1.0f+x)) ) );
if( x > 1.0f ) {
domerr:
	return( 0.0f );
}
if( x > 0.5f )
	return( 2.0f * ceph_asinf(  ceph_sqrtf(0.5f*(1.0f-x) ) ) );
	
return(1.5707963267948966192f - ceph_asinf(x) );
}


}
}








#endif /*__GMS_CEPHES_SIN_COS_H__*/