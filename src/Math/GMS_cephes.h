
#ifndef __GMS_CEPHES_H__
#define __GMS_CEPHES_H__

/*
    Slightly adapted version of CEPHES library.
    The rationale: removed the overhead of GLIBC calls for math.h declared
    mathematical functions (elementary and transcendental).
*/

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1989 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

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

#include "GSM_config.h"

#if !defined (DENORMAL)
    #define DENORMAL 1
#endif

#define EXPMSK 0x807f
#define MEXP 255
#define NBITS 24
#define BIG     1.44115188075855872E+17
#define MAXNUMF 3.4028234663852885981170418348451692544e38
#define MAXLOGF 88.72283905206835
#define MINLOGF -103.278929903431851103 /* log(2^-149) */
#define LOG2EF  1.44269504088896341
#define LOGE2F  0.693147180559945309
#define SQRTHF  0.707106781186547524
#define PIF     3.141592653589793238
#define PIO2F   1.5707963267948966192
#define PIO4F   0.7853981633974483096
#define MACHEPF 5.9604644775390625E-8

static int sgngamf;


__ATTR_ALWAYS_INLINE__
static inline
float ceph_ceilf( const float x ){
     float y = 0.0f;
     y = floorf( (float )x );
     if( y < x )
	 y += 1.0;
     return(y);
}

/* Bit clearing masks: */

static unsigned short bmask[] = {
       0xffff,
       0xfffe,
       0xfffc,
       0xfff8,
       0xfff0,
       0xffe0,
       0xffc0,
       0xff80,
       0xff00,
       0xfe00,
       0xfc00,
       0xf800,
       0xf000,
       0xe000,
       0xc000,
       0x8000,
       0x0000,
    };

__ATTR_ALWAYS_INLINE__
static inline
float ceph_floorf( float x ) {
unsigned short *p;
union
  {
    float y;
    unsigned short i[2];
  } u;
int e;
u.y = x;
/* find the exponent (power of 2) */
p = &u.i[1];
e = (( *p >> 7) & 0xff) - 0x7f;
p -= 1;
if( e < 0 ){
    if( u.y < 0 )
	 return( -1.0 );
    else
	 return( 0.0 );
}
e = (NBITS -1) - e;
/* clean out 16 bits at a time */
while( e >= 16 ) {
	*p++ = 0;
	 e -= 16;
}

/* clear the remaining bits */
if( e > 0 )
	*p &= bmask[e];

if( (x < 0) && (u.y != x) )
	u.y -= 1.0;

   return(u.y);
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
short *q;
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
union
  {
    float y;
    unsigned short i[2];
  } u;
short *q;
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
	u.y *= ldexpf(1.0f, e - 1);
	return(u.y);
#else
	return( 0.0 );
#endif
	}
*q |= (e & 0xff) << 7;
return(u.y);
}

/* Return 1 if the sign bit of x is 1, else 0.  */


__ATTR_ALWAYS_INLINE__
static inline
int ceph_signbitf(const float x) {
union
	{
	float f;
	short s[4];
	int i;
	} u;

u.f = x;

if( sizeof(int) == 4 ){
    return( u.i < 0 );
}
else{
     return( u.s[1] < 0 );
   }
	
}

/* Return 1 if x is a number that is Not a Number, else return 0.  */


__ATTR_ALWAYS_INLINE__
static inline
int ceph_isnanf(const float x) {
union
	{
	float f;
	unsigned short s[2];
	unsigned int i;
	} u;

u.f = x;
if( sizeof(int) == 4 ) {
	
     if( ((u.i & 0x7f800000) == 0x7f800000)
	    && ((u.i & 0x007fffff) != 0) )
		return 1;

	return(0);
}
else
    { /* size int not 4 */

	if( (u.s[1] & 0x7f80) == 0x7f80) {
	     if( ((u.s[1] & 0x007f) | u.s[0]) != 0 )
			return(1);
	}

	return(0);
     } /* size int not 4 */

}

/* Return 1 if x is not infinite and is not a NaN.  */

__ATTR_ALWAYS_INLINE__
static inline
int ceph_isfinitef(const float x) {
union{
	
	float f;
	unsigned short s[2];
	unsigned int i;
} u;
u.f = x;
if( sizeof(int) == 4 ) {
	if( (u.i & 0x7f800000) != 0x7f800000)
	   return 1;

	return(0);
}
else
    {

	if( (u.s[1] & 0x7f80) != 0x7f80)
	     return 1;

	return(0);
     }

}

/* Single precision inverse hyperbolic cosine
 * test interval: [1.0, 1.5]
 * trials: 10000
 * peak relative error: 1.7e-7
 * rms relative error: 5.0e-8
 *
 * Copyright (C) 1989 by Stephen L. Moshier.  All rights reserved.
 */


__ATTR_ALWAYS_INLINE__
static inline 
float ceph_acoshf(const float xx) {
float x, z;
x = xx;
if( x < 1.0 ) {
    return(0.0);
}
if( x > 1500.0 )
   return( ceph_logf(x) + 1.44269504088896341  );

z = x - 1.0;
if( z < 0.5 ) {
    z =
	(((( 1.7596881071E-3 * z
	  - 7.5272886713E-3) * z
	  + 2.6454905019E-2) * z
	  - 1.1784741703E-1) * z
	  + 1.4142135263E0) * ceph_sqrtf( z );
}
else
    {
	z = ceph_sqrtf( z*(x+1.0) );
	z = ceph_logf(x + z);
}
return( z );
}

/* Single precision circular arcsine
 * test interval: [-0.5, +0.5]
 * trials: 10000
 * peak relative error: 6.7e-8
 * rms relative error: 2.5e-8
 */

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
    x = sqrtf( z );
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
if( x < -1.0 )
	goto domerr;

if( x < -0.5) 
	return( 3.14159265358979323846
                - 2.0 * ceph_asinf(ceph_sqrtf(0.5*(1.0+x)) ) );
if( x > 1.0 ) {
domerr:
	return( 0.0 );
}
if( x > 0.5 )
	return( 2.0 * ceph_asinf(  ceph_sqrtf(0.5*(1.0-x) ) ) );
	
return(1.5707963267948966192 - ceph_asinf(x) );
}

/* Single precision circular arcsine
 * test interval: [-tan(pi/8), +tan(pi/8)]
 * trials: 10000
 * peak relative error: 7.7e-8
 * rms relative error: 2.9e-8
 */
 

__ATTR_ALWAYS_INLINE__
static inline 
float ceph_atanf(const float xx ) {
float x, y, z;
int sign;
x = xx;

/* make argument positive and save the sign */
if( xx < 0.0 ) {	
    sign = -1;
    x = -xx;
}
else
    {
     sign = 1;
     x = xx;
}
/* range reduction */
if( x > 2.414213562373095 ) {  /* tan 3pi/8 */
    y = 1.5707963267948966192;
    x = -( 1.0/x );
}
else if( x > 0.4142135623730950 ) { /* tan pi/8 */
	y = 0.7853981633974483096;
	x = (x-1.0)/(x+1.0);
}
else
    y = 0.0;

z = x * x;
y +=
((( 8.05374449538e-2 * z
  - 1.38776856032E-1) * z
  + 1.99777106478E-1) * z
  - 3.33329491539E-1) * z * x
  + x;
if( sign < 0 )
	y = -y;
return( y );
}


__ATTR_ALWAYS_INLINE__
static inline 
float ceph_atan2f(const float y,
             const float x ) {
float z, w;
int code;
code = 0;
if( x < 0.0 )
	code = 2;
if( y < 0.0 )
	code |= 1;
if( x == 0.0 ) {
      if( code & 1 ){
	  return( -1.5707963267948966192 );
      }
      if( y == 0.0 )
	     return( 0.0 );
      return(1.5707963267948966192 );
}

if( y == 0.0 ) {
	 if( code & 2 )
		return(3.14159265358979323846);
	 return( 0.0 );
}

switch( code ) {
	
	default:

	case 0:
	case 1: w = 0.0; break;
	case 2: w = 3.14159265358979323846; break;
	case 3: w = -3.14159265358979323846; break;
}
z = ceph_atanf( y/x );
return( w + z );
}


__ATTR_ALWAYS_INLINE__


static inline 
float ceph_atanhf(const float xx) {
float x, z;
x = xx;
if( x < 0 )
	z = -x;
else
	z = x;
if( z >= 1.0 ) {
	if( x == 1.0 )
		return(3.4028234663852885981170418348451692544e38);
	if( x == -1.0 )
		return( -3.4028234663852885981170418348451692544e38);
   	return(3.4028234663852885981170418348451692544e38);
}

if( z < 1.0e-4 )
	return(x);

if( z < 0.5 ) {
	z = x * x;
	z =
	(((( 1.81740078349E-1 * z
	  + 8.24370301058E-2) * z
	  + 1.46691431730E-1) * z
	  + 1.99782164500E-1) * z
	  + 3.33337300303E-1) * z * x
	  + x;
}
else
    {
	z = 0.5 * ceph_logf( (1.0+x)/(1.0-x) );
}
return( z );
}




__ATTR_ALWAYS_INLINE__


static inline 
float ceph_cbrtf(const float xx) {
int e, rem, sign;
float x, z;
x = xx;
if( x == 0 )
	return( 0.0 );
if( x > 0 )
	sign = 1;
else{
    sign = -1;
    x = -x;
}

z = x;
/* extract power of 2, leaving
 * mantissa between 0.5 and 1
 */
x = ceph_frexpf( x, &e );

/* Approximate cube root of number between .5 and 1,
 * peak relative error = 9.2e-6
 */
x = (((-0.13466110473359520655053  * x
      + 0.54664601366395524503440 ) * x
      - 0.95438224771509446525043 ) * x
      + 1.1399983354717293273738  ) * x
      + 0.40238979564544752126924;

/* exponent divided by 3 */
if( e >= 0 ) {
    rem = e;
    e /= 3;
    rem -= 3*e;
    if( rem == 1 )
	x *= 1.25992104989487316477;
    else if( rem == 2 )
		x *= 1.58740105196819947475;
}
/* argument less than 1 */
else {
	
   e = -e;
   rem = e;
   e /= 3;
   rem -= 3*e;
   if( rem == 1 )
       x /= 1.25992104989487316477;
   else if( rem == 2 )
       x /= 1.58740105196819947475;
       e = -e;
}
/* multiply by power of 2 */
x = ceph_ldexpf( x, e );
/* Newton iteration */
x -= ( x - (z/(x*x)) ) * 0.333333333333;
if( sign < 0 )
	x = -x;
return(x);
}


__ATTR_ALWAYS_INLINE__


static inline 
float ceph_coshf(const float xx) {
float x, y;
x = xx;
if( x < 0 )
	x = -x;
if( x > 88.72283905206835)
    return( 3.4028234663852885981170418348451692544e38);
y = ceph_expf(x);
y = y + 1.0/y;
return( 0.5*y );
}



__ATTR_ALWAYS_INLINE__


static inline
float ceph_expf(const float xx) {
float x, z;
int n;
x = xx;
if( x > 88.72283905206835) {
    return(3.4028234663852885981170418348451692544e38);
}
if( x < -103.278929903431851103) {
    return(0.0);
}

/* Express e**x = e**g 2**n
 *   = e**g e**( n loge(2) )
 *   = e**( g + n loge(2) )
 */
z = ceph_floorf(1.44269504088896341 * x + 0.5 ); /* floor() truncates toward -infinity. */
x -= z * 0.693359375;
x -= z * -2.12194440e-4;
n = z;

z = x * x;
/* Theoretical peak relative error in [-0.5, +0.5] is 4.2e-9. */
z =
((((( 1.9875691500E-4  * x
   + 1.3981999507E-3) * x
   + 8.3334519073E-3) * x
   + 4.1665795894E-2) * x
   + 1.6666665459E-1) * x
   + 5.0000001201E-1) * z
   + x
   + 1.0;

/* multiply by power of 2 */
x = ceph_ldexpf( z, n );
return( x );
}


__ATTR_ALWAYS_INLINE__


static inline
float ceph_expx2f(const float x) {
  float u, u1, m;
  if (x < 0.0)
    x = -x;

  /* Represent x as an exact multiple of 1/32 plus a residual.  */
  m = .03125f * ceph_floorf(32.0f * x + 0.5f);
  x -= m;
  /* x**2 = m**2 + 2mf + f**2 */
  u = m * m;
  u1 = 2 * m * x  +  x * x;

  if ((u+u1) > 88.72283905206835 )
    return (3.4028234663852885981170418348451692544e38);

  /* u is exact, u1 is small.  */
  u = ceph_expf(u) * ceph_expf(u1);
  return(u);
}



__ATTR_ALWAYS_INLINE__


static inline
float logf( const float xx ) {
register float y;
float x, z, fe;
int e;
x = xx;
fe = 0.0;
/* Test for domain */
if( x <= 0.0 ) {
    if( x == 0.0 )
	   ;
    else
       return(-103.278929903431851103);
}
x = ceph_frexpf( x, &e );
if( x < 0.707106781186547524) {
    e -= 1;
    x = x + x - 1.0; /*  2x - 1  */
}	
else{
    x = x - 1.0;
}
z = x * x;
/* 3.4e-9 */
/*
p = logfcof;
y = *p++ * x;
for( i=0; i<8; i++ )
	{
	y += *p++;
	y *= x;
	}
y *= z;
*/

y =
(((((((( 7.0376836292E-2 * x
- 1.1514610310E-1) * x
+ 1.1676998740E-1) * x
- 1.2420140846E-1) * x
+ 1.4249322787E-1) * x
- 1.6668057665E-1) * x
+ 2.0000714765E-1) * x
- 2.4999993993E-1) * x
+ 3.3333331174E-1) * x * z;

if( e ){
   fe = e;
   y += -2.12194440e-4 * fe;
}

y +=  -0.5 * z;  /* y - 0.5 x^2 */
z = x + y;   /* ... + x  */
if( e )
	z += 0.693359375 * fe;
return( z );
}

#define ceph_fabsf(x) ( (x) < 0 ? -(x) : (x) )


__ATTR_ALWAYS_INLINE__


static inline
float betaf(const  float aa, const float bb ) {
float a, b, y;
int sign;
sign = 1;
a = aa;
b = bb;
if( a <= 0.0 ) {
     if( a == ceph_floorf(a) )
		goto over;
}
if( b <= 0.0 ) {
     if( b == ceph_floorf(b) )
		goto over;
}

y = a + b;
if( ceph_fabsf(y) > 34.84425627277176174) {
	y = ceph_lgamf(y);
	sign *= sgngamf; /* keep track of the sign */
	y = ceph_lgamf(b) - y;
	sign *= sgngamf;
	y = ceph_lgamf(a) + y;
	sign *= ceph_sgngamf;
	if( y > 88.72283905206835) {
		
over:
	    return( sign * 3.4028234663852885981170418348451692544e38 );
	 }
      return( sign * ceph_expf(y) );
}

y = ceph_gammaf(y);
if( y == 0.0 )
	goto over;

if( a > b ) {
    y =  ceph_gammaf(a)/y;
    y *= ceph_gammaf(b);
}
else{
    y = ceph_gammaf(b)/y;
    y *= ceph_gammaf(a);
}
return(y);
}



__ATTR_ALWAYS_INLINE__


static inline
float ceph_chbevlf( float x, float * __restrict array, int n ) {
float b0, b1, b2, *p;
int i;
p = array;
b0 = *p++;
b1 = 0.0;
i = n - 1;
do{
   	b2 = b1;
	b1 = b0;
	b0 = x * b1  -  b2  + *p++;
}
while( --i );
return( 0.5*(b0-b2) );
}

/*
   DESCRIPTION:
 *
 * Finds the Chi-square argument x such that the integral
 * from x to infinity of the Chi-square density is equal
 * to the given cumulative probability y.
 *
 * This is accomplished using the inverse gamma integral
 * function and the relation
 *
 *    x/2 = igami( df/2, y );
*/




__ATTR_ALWAYS_INLINE__


static inline
float ceph_chdtrcf(const float dff, const float xx){
float df, x;
df = dff;
x = xx;
if( (x < 0.0) || (df < 1.0) ) {
    return(0.0);
}
return( ceph_igamcf( 0.5*df, 0.5*x ) );
}



__ATTR_ALWAYS_INLINE__


static inline
float chdtrf(const float dff, const float xx) {
float df, x;
df = dff;
x = xx;
if( (x < 0.0) || (df < 1.0) ) {
	return(0.0);
}
return( ceph_igamf( 0.5*df, 0.5*x ) );
}



__ATTR_ALWAYS_INLINE__


static inline
float chdtrif( const float dff, const float yy ) {
float y, df, x;
y = yy;
df = dff;
if( (y < 0.0) || (y > 1.0) || (df < 1.0) ) {
	return(0.0);
}
x = ceph_igamif( 0.5 * df, y );
return( 2.0 * x );
}

/*
DESCRIPTION:
 *
 * Approximates the integral
 *
 *                             x
 *                             -
 *                      2     | |        2
 *  dawsn(x)  =  exp( -x  )   |    exp( t  ) dt
 *                          | |
 *                           -
 *                           0
 *
 * Three different rational approximations are employed, for
 * the intervals 0 to 3.25; 3.25 to 6.25; and 6.25 up.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,10        50000       4.4e-7      6.3e-8
*/



//extern float PIF, MACHEPF;




__ATTR_ALWAYS_INLINE__


static inline
float dawsnf( const float xxx ) {
/* Dawson's integral, interval 0 to 3.25 */
const static float AN[10] = {
 1.13681498971755972054E-11,
 8.49262267667473811108E-10,
 1.94434204175553054283E-8,
 9.53151741254484363489E-7,
 3.07828309874913200438E-6,
 3.52513368520288738649E-4,
-8.50149846724410912031E-4,
 4.22618223005546594270E-2,
-9.17480371773452345351E-2,
 9.99999999999999994612E-1,
};
const static float AD[11] = {
 2.40372073066762605484E-11,
 1.48864681368493396752E-9,
 5.21265281010541664570E-8,
 1.27258478273186970203E-6,
 2.32490249820789513991E-5,
 3.25524741826057911661E-4,
 3.48805814657162590916E-3,
 2.79448531198828973716E-2,
 1.58874241960120565368E-1,
 5.74918629489320327824E-1,
 1.00000000000000000539E0,
};

/* interval 3.25 to 6.25 */
const static float BN[11] = {
 5.08955156417900903354E-1,
-2.44754418142697847934E-1,
 9.41512335303534411857E-2,
-2.18711255142039025206E-2,
 3.66207612329569181322E-3,
-4.23209114460388756528E-4,
 3.59641304793896631888E-5,
-2.14640351719968974225E-6,
 9.10010780076391431042E-8,
-2.40274520828250956942E-9,
 3.59233385440928410398E-11,
};
const static float BD[10] = {
/*  1.00000000000000000000E0,*/
-6.31839869873368190192E-1,
 2.36706788228248691528E-1,
-5.31806367003223277662E-2,
 8.48041718586295374409E-3,
-9.47996768486665330168E-4,
 7.81025592944552338085E-5,
-4.55875153252442634831E-6,
 1.89100358111421846170E-7,
-4.91324691331920606875E-9,
 7.18466403235734541950E-11,
};

/* 6.25 to infinity */
const static float CN[5] = {
-5.90592860534773254987E-1,
 6.29235242724368800674E-1,
-1.72858975380388136411E-1,
 1.64837047825189632310E-2,
-4.86827613020462700845E-4,
};
const static float CD[5] = {
/* 1.00000000000000000000E0,*/
-2.69820057197544900361E0,
 1.73270799045947845857E0,
-3.93708582281939493482E-1,
 3.44278924041233391079E-2,
-9.73655226040941223894E-4,
};

float xx, x, y;
int sign;
xx = xxx;
sign = 1;
if( xx < 0.0 ) {
	
      sign = -1;
      xx = -xx;
}

if( xx < 3.25 ) {
      x = xx*xx;
      y = xx * polevlf( x, AN, 9 )/polevlf( x, AD, 10 );
      return( sign * y );
}
x = 1.0/(xx*xx);
if( xx < 6.25 ){
       y = 1.0/xx + x * polevlf( x, BN, 10) / (p1evlf( x, BD, 10) * xx);
       return( sign * 0.5 * y );
}

if( xx > 1.0e9 )
	return( (sign * 0.5)/xx );

/* 6.25 to infinity */
y = 1.0/xx + x * polevlf( x, CN, 4) / (p1evlf( x, CD, 5) * xx);
return( sign * 0.5 * y );
}

/*DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *                phi
 *                 -
 *                | |
 *                |                   2
 * E(phi\m)  =    |    sqrt( 1 - m sin t ) dt
 *                |
 *              | |    
 *               -
 *                0
 *
 * of amplitude phi and modulus m, using the arithmetic -
 * geometric mean algorithm.
*/


__ATTR_ALWAYS_INLINE__


static inline
float ellief( const float phia, const float ma ) {
float phi, m, a, b, c, e, temp;
float lphi, t;
int d, mod;

phi = phia;
m = ma;
if( m == 0.0 )
	return( phi );
if( m == 1.0 )
	return( ceph_sinf(phi) );
lphi = phi;
if( lphi < 0.0 )
	lphi = -lphi;
a = 1.0;
b = 1.0 - m;
b = ceph_sqrtf(b);
c = ceph_sqrtf(m);
d = 1;
e = 0.0;
t = ceph_tanf( lphi );
mod = (lphi +1.5707963267948966192 )/3.141592653589793238;

while( fabsf(c/a) > 5.9604644775390625E-8)
	{
	temp = b/a;
	lphi = lphi + ceph_atanf(t*temp) + mod * 3.141592653589793238;
	mod = (lphi + 1.5707963267948966192)/3.141592653589793238;
	t = t * ( 1.0 + temp )/( 1.0 - temp * t * t );
	c = 0.5 * ( a - b );
	temp = ceph_sqrtf( a * b );
	a = 0.5 * ( a + b );
	b = temp;
	d += d;
	e += c * ceph_sinf(lphi);
	}

b = 1.0 - m;
temp = ellpef(b)/ellpkf(b);
temp *= (ceph_atanf(t) + mod * 3.141592653589793238)/(d * a);
temp += e;
if( phi < 0.0 )
	temp = -temp;
return( temp );
}

/*
    SYNOPSIS:
 *
 * float phi, m, y, ellikf();
 *
 * y = ellikf( phi, m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *                phi
 *                 -
 *                | |
 *                |           dt
 * F(phi\m)  =    |    ------------------
 *                |                   2
 *              | |    sqrt( 1 - m sin t )
 *               -
 *                0
 *
 * of amplitude phi and modulus m, using the arithmetic -
 * geometric mean algorithm.
 *
*/


__ATTR_ALWAYS_INLINE__


static inline
float ellikf(const  float phia, const float ma ){
float phi, m, a, b, c, temp;
float t;
int d, mod, sign;
phi = phia;
m = ma;
if( m == 0.0 )
	return( phi );
if( phi < 0.0 ){
	
	phi = -phi;
	sign = -1;
}
else{
	sign = 0;
}
a = 1.0;
b = 1.0 - m;
if( b == 0.0 )
	return(ceph_logf(ceph_tanf(0.5*(1.5707963267948966192 + phi))));
b = ceph_sqrtf(b);
c = ceph_sqrtf(m);
d = 1;
t = ceph_tanf( phi );
mod = (phi + 1.5707963267948966192)/3.141592653589793238;
while(fabsf(c/a) > 5.9604644775390625E-8) {
	temp = b/a;
	phi = phi + ceph_atanf(t*temp) + mod * PIF;
	mod = (phi + 1.5707963267948966192)/3.141592653589793238;
	t = t * ( 1.0 + temp )/( 1.0 - temp * t * t );
	c = ( a - b )/2.0;
	temp = ceph_sqrtf( a * b );
	a = ( a + b )/2.0;
	b = temp;
	d += d;
	}

temp = (ceph_atanf(t) + mod * 3.141592653589793238)/(d * a);
if( sign < 0 )
	temp = -temp;
return( temp );
}

/*
 SYNOPSIS:
 *
 * float m1, y, ellpef();
 *
 * y = ellpef( m1 );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *            pi/2
 *             -
 *            | |                 2
 * E(m)  =    |    sqrt( 1 - m sin t ) dt
 *          | |    
 *           -
 *            0
 *
 * Where m = 1 - m1, using the approximation
 *
 *      P(x)  -  x log x Q(x).
 *
 * Though there are no singularities, the argument m1 is used
 * rather than m for compatibility with ellpk().
 *
 * E(1) = 1; E(0) = pi/2
*/



__ATTR_ALWAYS_INLINE__


static inline
float ellpef( const float xx){
const static float P[] = {
  1.53552577301013293365E-4,
  2.50888492163602060990E-3,
  8.68786816565889628429E-3,
  1.07350949056076193403E-2,
  7.77395492516787092951E-3,
  7.58395289413514708519E-3,
  1.15688436810574127319E-2,
  2.18317996015557253103E-2,
  5.68051945617860553470E-2,
  4.43147180560990850618E-1,
  1.00000000000000000299E0
};
const static float Q[] = {
  3.27954898576485872656E-5,
  1.00962792679356715133E-3,
  6.50609489976927491433E-3,
  1.68862163993311317300E-2,
  2.61769742454493659583E-2,
  3.34833904888224918614E-2,
  4.27180926518931511717E-2,
  5.85936634471101055642E-2,
  9.37499997197644278445E-2,
  2.49999999999888314361E-1
};

float x;
x = xx;
if( (x <= 0.0) || (x > 1.0) ){
	  if( x == 0.0 )
		return( 1.0 );
          return( 0.0 );
}
return( polevlf(x,P,10) - ceph_logf(x) * (x * polevlf(x,Q,9)) );
}

/*
 	Jacobian Elliptic Functions
 *
 *
 *
 * SYNOPSIS:
 *
 * float u, m, sn, cn, dn, phi;
 * int ellpj();
 *
 * ellpj( u, m, _&sn, _&cn, _&dn, _&phi );
 *
 *
*/


__ATTR_ALWAYS_INLINE__


static inline
int ellpjf( float uu, float mm,
            float * __restrict sn,
	    float * __restrict cn,
	    float * __restrict dn,
	    float * __restrict ph ) {
  float a[10];
  float c[10];
float u, m, ai, b, phi, t, twon;
int i;
u = uu;
m = mm;
/* Check for special cases */
if( m < 0.0 || m > 1.0 ){
     return(-1);
}
if( m < 1.0e-5 ) {
	t = ceph_sinf(u);
	b = ceph_cosf(u);
	ai = 0.25 * m * (u - t*b);
	*sn = t - ai*b;
	*cn = b + ai*t;
	*ph = u - ai;
	*dn = 1.0 - 0.5*m*t*t;
	return(0);
}
if( m >= 0.99999 ) {
	ai = 0.25 * (1.0-m);
	b = ceph_coshf(u);
	t = ceph_tanhf(u);
	phi = 1.0/b;
	twon = b * ceph_sinhf(u);
	*sn = t + ai * (twon - u)/(b*b);
	*ph = 2.0*ceph_atanf(ceph_expf(u)) - 1.5707963267948966192 + ai*(twon - u)/b;
	ai *= t * phi;
	*cn = phi - ai * (twon - u);
	*dn = phi + ai * (twon + u);
	return(0);
}

/*	A. G. M. scale		*/
a[0] = 1.0;
b = ceph_sqrtf(1.0 - m);
c[0] = ceph_sqrtf(m);
twon = 1.0;
i = 0;
while(ceph_fabsf((c[i]/a[i])) > 5.9604644775390625E-8) {
	
	if( i > 8 ) 
	     break;
	ai = a[i];
	++i;
	c[i] = 0.5 * ( ai - b );
	t = ceph_sqrtf( ai * b );
	a[i] = 0.5 * ( ai + b );
	b = t;
	twon += twon;
}


/* backward recurrence */
phi = twon * a[i] * u;
do{
	
	t = c[i] * ceph_sinf(phi) / a[i];
	b = phi;
	phi = 0.5 * (ceph_asinf(t) + phi);
}
while( --i );
t = ceph_sinf(phi);
*sn = t;
*cn = ceph_cosf(phi);
/* Thanks to Hartmut Henkel for reporting a bug here:  */
*dn = ceph_sqrtf(1.0 - m * t * t);
*ph = phi;
return(0);
}

/*
    SYNOPSIS:
 *
 * float m1, y, ellpkf();
 *
 * y = ellpkf( m1 );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *            pi/2
 *             -
 *            | |
 *            |           dt
 * K(m)  =    |    ------------------
 *            |                   2
 *          | |    sqrt( 1 - m sin t )
 *           -
 *            0
 *
 * where m = 1 - m1, using the approximation
 *
 *     P(x)  -  log x Q(x).
*/


__ATTR_ALWAYS_INLINE__


static inline
float ellpkf(const float xx){
  const float P[] =
{
 1.37982864606273237150E-4,
 2.28025724005875567385E-3,
 7.97404013220415179367E-3,
 9.85821379021226008714E-3,
 6.87489687449949877925E-3,
 6.18901033637687613229E-3,
 8.79078273952743772254E-3,
 1.49380448916805252718E-2,
 3.08851465246711995998E-2,
 9.65735902811690126535E-2,
 1.38629436111989062502E0
};

  const float Q[] =
{
 2.94078955048598507511E-5,
 9.14184723865917226571E-4,
 5.94058303753167793257E-3,
 1.54850516649762399335E-2,
 2.39089602715924892727E-2,
 3.01204715227604046988E-2,
 3.73774314173823228969E-2,
 4.88280347570998239232E-2,
 7.03124996963957469739E-2,
 1.24999999999870820058E-1,
 4.99999999999999999821E-1
};
float x;
x = xx;
if( (x < 0.0) || (x > 1.0) ) {
	
	return( 0.0 );
}

if( x > 5.9604644775390625E-8){
    return( polevlf(x,P,10) - ceph_logf(x) * polevlf(x,Q,10) );
}
else
    {
	if( x == 0.0 ) {
	   return(3.4028234663852885981170418348451692544e38);
       }
	else{
	     return(1.3862943611198906188E0 - 0.5 * ceph_logf(x));
	}
    }
}

/*
   SYNOPSIS:
 *
 * float x, y, exp2f();
 *
 * y = exp2f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns 2 raised to the x power.
 *
 * Range reduction is accomplished by separating the argument
 * into an integer k and fraction f such that
 *     x    k  f
 *    2  = 2  2.
 *
 * A polynomial approximates 2**x in the basic range [-0.5, 0.5].
*/


__ATTR_ALWAYS_INLINE__


static inline
float ceph_exp2f( const float xx ){
 const float P[] = {
 1.535336188319500E-004,
 1.339887440266574E-003,
 9.618437357674640E-003,
 5.550332471162809E-002,
 2.402264791363012E-001,
 6.931472028550421E-001
};
float x, px;
int i0;
x = xx;
if( x > 127.0){
	return(3.4028234663852885981170418348451692544e38);
}
if( x < -127.0) {
	return(0.0);
}
/* The following is necessary because range reduction blows up: */
if( x == 0 )
	return(1.0);

/* separate into integer and fractional parts */
px = ceph_floorf(x);
i0 = px;
x = x - px;
if( x > 0.5 ) {
     i0 += 1;
     x -= 1.0;
}

/* rational approximation
 * exp2(x) = 1.0 +  xP(x)
 */
px = 1.0 + x * polevlf( x, P, 5 );
/* scale by power of 2 */
px = ceph_ldexpf( px, i0 );
return(px);
}

/*
   SYNOPSIS:
 *
 * float x, y, exp10f();
 *
 * y = exp10f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns 10 raised to the x power.
 *
 * Range reduction is accomplished by expressing the argument
 * as 10**x = 2**n 10**f, with |f| < 0.5 log10(2).
 * A polynomial approximates 10**f.
 *
*/

__ATTR_ALWAYS_INLINE__


static inline
float ceph_exp10f(const float xx) {
  const float P[] = {
 2.063216740311022E-001,
 5.420251702225484E-001,
 1.171292686296281E+000,
 2.034649854009453E+000,
 2.650948748208892E+000,
 2.302585167056758E+000
};
float x, px, qx;
short n;
x = xx;
if( x > 38.230809449325611792) {
	return(3.4028234663852885981170418348451692544e38);
}
if( x < -38.230809449325611792){ /* Would like to use MINLOG but can't */
	return(0.0);
}
/* The following is necessary because range reduction blows up: */
if( x == 0 )
	return(1.0);

/* Express 10**x = 10**g 2**n
 *   = 10**g 10**( n log10(2) )
 *   = 10**( g + n log10(2) )
 */
px = x * 3.32192809488736234787e0;
qx = ceph_floorf( px + 0.5 );
n = qx;
x -= qx * 3.00781250000000000000E-1;
x -= qx * 2.48745663981195213739E-4;

/* rational approximation for exponential
 * of the fractional part:
 * 10**x - 1  =  2x P(x**2)/( Q(x**2) - P(x**2) )
 */
px = 1.0 + x * polevlf( x, P, 5 );
/* multiply by power of 2 */
x = ceph_ldexpf( px, n );
return(x);
}

/*
 SYNOPSIS:
 *
 * int n;
 * float x, y, expnf();
 *
 * y = expnf( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the exponential integral
 *
 *                 inf.
 *                   -
 *                  | |   -xt
 *                  |    e
 *      E (x)  =    |    ----  dt.
 *       n          |      n
 *                | |     t
 *                 -
 *                  1
 *
 *
 * Both n and x must be nonnegative.
 *
 * The routine employs either a power series, a continued
 * fraction, or an asymptotic formula depending on the
 * relative values of n and x.
*/


__ATTR_ALWAYS_INLINE__


static inline
float expnf( const int n, const float xx ) {
float x, ans, r, t, yk, xk;
float pk, pkm1, pkm2, qk, qkm1, qkm2;
float psi, z;
int i, k;
constexpr float big = 16777216.0;
x = xx;
if( n < 0 )
	goto domerr;
if( x < 0 ) {
	
domerr:
	return(3.4028234663852885981170418348451692544e38);
}
if( x > 88.72283905206835)
	return( 0.0 );

if( x == 0.0 ) {
      if( n < 2 ) {
	  return(3.4028234663852885981170418348451692544e38);
	}
	else
		return( 1.0/(n-1.0) );
	}

if( n == 0 )
	return( ceph_expf(-x)/x );

/*							expn.c	*/
/*		Expansion for large n		*/

if( n > 5000 ) {
	
	xk = x + n;
	yk = 1.0 / (xk * xk);
	t = n;
	ans = yk * t * (6.0 * x * x  -  8.0 * t * x  +  t * t);
	ans = yk * (ans + t * (t  -  2.0 * x));
	ans = yk * (ans + t);
	ans = (ans + 1.0) * ceph_expf( -x ) / xk;
	goto done;
}

if( x > 1.0 )
	goto cfrac;
 
/*							expn.c	*/

/*		Power series expansion		*/

psi = -0.57721566490153286060 - ceph_logf(x);
for( i=1; i<n; i++ )
	psi = psi + 1.0/i;

z = -x;
xk = 0.0;
yk = 1.0;
pk = 1.0 - n;
if( n == 1 )
	ans = 0.0;
else
	ans = 1.0/pk;
do{
	xk += 1.0;
	yk *= z/xk;
	pk += 1.0;
	if( pk != 0.0 ) {
	      ans += yk/pk;
	}
	if( ans != 0.0 )
		t = ceph_fabsf(yk/ans);
	else
		t = 1.0;
}
while(t > 5.9604644775390625E-8);
k = xk;
t = n;
r = n - 1;
ans = (ceph_powf(z, r) * psi / ceph_gammaf(t)) - ans;
goto done;

/*							expn.c	*/
/*		continued fraction		*/
cfrac:
k = 1;
pkm2 = 1.0;
qkm2 = x;
pkm1 = 1.0;
qkm1 = x + n;
ans = pkm1/qkm1;

do{
	
    k += 1;
    if( k & 1 ) {
	yk = 1.0;
	xk = n + (k-1)/2;
    }
    else
	{
	yk = x;
	xk = k/2;
    }
     pk = pkm1 * yk  +  pkm2 * xk;
     qk = qkm1 * yk  +  qkm2 * xk;
     if( qk != 0 ) {
	 r = pk/qk;
	 t = ceph_fabsf( (ans - r)/r );
	 ans = r;
     }
     else
	t = 1.0;
	pkm2 = pkm1;
	pkm1 = pk;
	qkm2 = qkm1;
	qkm1 = qk;
if( ceph_fabsf(pk) > big ) {
      pkm2 *= 5.9604644775390625E-8;
      pkm1 *= 5.9604644775390625E-8;
      qkm2 *= 5.9604644775390625E-8;
      qkm1 *= 5.9604644775390625E-8;
    }
}
while( t > MACHEPF );
ans *= expf( -x );
done:
return( ans );
}

/*
   *
 * SYNOPSIS:
 *
 * double x, y, expx2f();
 *
 * y = expx2f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes y = exp(x*x) while suppressing error amplification
 * that would ordinarily arise from the inexactness of the argument x*x.
 * 
*/


__ATTR_ALWAYS_INLINE__


static inline
float expx2f (const float x) {
  float u, u1, m;
  if (x < 0)
    x = -x;
  /* Represent x as an exact multiple of 1/32 plus a residual.  */
  m = .03125f * ceph_floorf(32.0f * x + 0.5f);
  x -= m;
  /* x**2 = m**2 + 2mf + f**2 */
  u = m * m;
  u1 = 2 * m * x  +  x * x;
  if ((u+u1) > 88.72283905206835)
    return (3.4028234663852885981170418348451692544e38);

  /* u is exact, u1 is small.  */
  u = ceph_expf(u) * ceph_expf(u1);
  return(u);
}

/*
  SYNOPSIS:
 *
 * float x, S, C;
 * void fresnlf();
 *
 * fresnlf( x, _&S, _&C );
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the Fresnel integrals
 *
 *           x
 *           -
 *          | |
 * C(x) =   |   cos(pi/2 t**2) dt,
 *        | |
 *         -
 *          0
 *
 *           x
 *           -
 *          | |
 * S(x) =   |   sin(pi/2 t**2) dt.
 *        | |
 *         -
 *          0
 *
*/


__ATTR_ALWAYS_INLINE__


static inline
void fresnlf(const  float xxa,
             float * __restrict ssa,
	     float *__restrict cca ){
/* S(x) for small x */
 const float sn[7] = {
 1.647629463788700E-009,
-1.522754752581096E-007,
 8.424748808502400E-006,
-3.120693124703272E-004,
 7.244727626597022E-003,
-9.228055941124598E-002,
 5.235987735681432E-001
};

/* C(x) for small x */
 const float cn[7] = {
 1.416802502367354E-008,
-1.157231412229871E-006,
 5.387223446683264E-005,
-1.604381798862293E-003,
 2.818489036795073E-002,
-2.467398198317899E-001,
 9.999999760004487E-001
};


/* Auxiliary function f(x) */
 const float fn[8] = {
-1.903009855649792E+012,
 1.355942388050252E+011,
-4.158143148511033E+009,
 7.343848463587323E+007,
-8.732356681548485E+005,
 8.560515466275470E+003,
-1.032877601091159E+002,
 2.999401847870011E+000
};

/* Auxiliary function g(x) */
 const float gn[8] = {
-1.860843997624650E+011,
 1.278350673393208E+010,
-3.779387713202229E+008,
 6.492611570598858E+006,
-7.787789623358162E+004,
 8.602931494734327E+002,
-1.493439396592284E+001,
 9.999841934744914E-001
};
float f, g, cc, ss, c, s, t, u, x, x2;
/*debug double t1;*/
x = xxa;
x = ceph_fabsf(x);
x2 = x * x;
if( x2 < 2.5625 ) {
	t = x2 * x2;
	ss = x * x2 * polevlf( t, sn, 6);
	cc = x * polevlf( t, cn, 6);
	goto done;
}

if( x > 36974.0 ) {
	
	cc = 0.5;
	ss = 0.5;
	goto done;
}


/*		Asymptotic power series auxiliary functions
 *		for large argument
 */
	x2 = x * x;
	t = 3.141592653589793238 * x2;
	u = 1.0/(t * t);
	t = 1.0/t;
	f = 1.0 - u * polevlf( u, fn, 7);
	g = t * polevlf( u, gn, 7);

	t = 1.5707963267948966192 * x2;
	c = ceph_cosf(t);
	s = ceph_sinf(t);
	t = 3.141592653589793238 * x;
	cc = 0.5  +  (f * s  -  g * c)/t;
	ss = 0.5  -  (f * c  +  g * s)/t;

done:
if( xxa < 0.0 ) {
	cc = -cc;
	ss = -ss;
}

*cca = cc;
*ssa = ss;
return 0;
}

/*
    SYNOPSIS:
 *
 * float x, y, gammaf();
 * extern int sgngamf;
 *
 * y = gammaf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns gamma function of the argument.  The result is
 * correctly signed, and the sign (+1 or -1) is also
 * returned in a global (extern) variable named sgngamf.
 * This same variable is also filled in by the logarithmic
 * gamma function lgam().
 *
 * Arguments between 0 and 10 are reduced by recurrence and the
 * function is approximated by a polynomial function covering
 * the interval (2,3).  Large arguments are handled by Stirling's
 * formula. Negative arguments are made positive using
 * a reflection formula.  
 *
*/

/* Gamma function computed by Stirling's formula,
 * sqrt(2 pi) x^(x-.5) exp(-x) (1 + 1/x P(1/x))
 * The polynomial STIR is valid for 33 <= x <= 172.
 */

__ATTR_ALWAYS_INLINE__


static inline
float stirf(const float xx) {
const float STIR[] = {
-2.705194986674176E-003,
 3.473255786154910E-003,
 8.333331788340907E-002
};
float x, y, w, v;
x = xx;
w = 1.0/x;
w = 1.0 + w * polevlf( w, STIR, 2 );
y = ceph_expf( -x );
if( x > 26.77) {
	 /* Avoid overflow in pow() */
	v = ceph_powf( x, 0.5 * x - 0.25 );
	y *= v;
	y *= v;
}
else{
    y = ceph_powf( x, x - 0.5 ) * y;
}
y = 2.50662827463100050242 * y * w;
return( y );
}



__ATTR_ALWAYS_INLINE__


static inline
float ceph_gammaf(const float xx ) {
/* gamma(x+2), 0 < x < 1 */
 const float P[] = {
 1.536830450601906E-003,
 5.397581592950993E-003,
 4.130370201859976E-003,
 7.232307985516519E-002,
 8.203960091619193E-002,
 4.117857447645796E-001,
 4.227867745131584E-001,
 9.999999822945073E-001
};
float p, q, x, z, nz;
int i, direction, negative;
x = xx;
sgngamf = 1;
negative = 0;
nz = 0.0;
if( x < 0.0 ) {
	negative = 1;
	q = -x;
	p = ceph_floorf(q);
	if( p == q )
		goto goverf;
	i = p;
	if((i & 1) == 0 )
		sgngamf = -1;
	nz = q - p;
	if( nz > 0.5 ) {
	        p += 1.0;
		nz = q - p;
	}
	nz = q * ceph_sinf( 3.141592653589793238 * nz );
	if( nz == 0.0 ) {
		
goverf:
	     return(sgngamf * 3.4028234663852885981170418348451692544e38);
        }
	if( nz < 0 )
		nz = -nz;
	x = q;
	}
if( x >= 10.0 ) {
    z = stirf(x);
}
if( x < 2.0 )
	direction = 1;
else
	direction = 0;
z = 1.0;
while( x >= 3.0 ) {
        x -= 1.0;
	z *= x;
}
/*
while( x < 0.0 )
	{
	if( x > -1.E-4 )
		goto small;
	z *=x;
	x += 1.0;
	}
*/
while( x < 2.0 ) {
	
	if( x < 1.e-4 )
		goto small;
	z *=x;
	x += 1.0;
}

if( direction )
	z = 1.0/z;

if( x == 2.0 )
	return(z);

x -= 2.0;
p = z * polevlf( x, P, 7 );

gdone:

if( negative ) {
    p = sgngamf * 3.141592653589793238/(nz * p );
}
return(p);

small:
if( x == 0.0 ){
	return(3.4028234663852885981170418348451692544e38);
}
else{
	p = z / ((1.0 + 0.5772156649015329 * x) * x);
	goto gdone;
    }
}


__ATTR_ALWAYS_INLINE__


static inline
float lgamf(const float xx ){
/* log gamma(x+2), -.5 < x < .5 */
 const float B[] = {
 6.055172732649237E-004,
-1.311620815545743E-003,
 2.863437556468661E-003,
-7.366775108654962E-003,
 2.058355474821512E-002,
-6.735323259371034E-002,
 3.224669577325661E-001,
 4.227843421859038E-001
};

/* log gamma(x+1), -.25 < x < .25 */
 const float C[] = {
 1.369488127325832E-001,
-1.590086327657347E-001,
 1.692415923504637E-001,
-2.067882815621965E-001,
 2.705806208275915E-001,
-4.006931650563372E-001,
 8.224670749082976E-001,
-5.772156501719101E-001
};
float p, q, w, z, x;
float nx, tx;
int i, direction;
sgngamf = 1;
x = xx;
if( x < 0.0 ) {
	q = -x;
	w = ceph_lgamf(q); /* note this modifies sgngam! */
	p = ceph_floorf(q);
	if( p == q )
		goto loverf;
	i = p;
	if( (i & 1) == 0 )
		sgngamf = -1;
	else
		sgngamf = 1;
	z = q - p;
	if( z > 0.5 ) {
		p += 1.0;
		z = p - q;
	}
	z = q * ceph_sinf(3.141592653589793238 * z);
	if(z == 0.0)
		goto loverf;
	z = -ceph_logf( 0.318309886183790671538*z ) - w;
	return( z );
	}

if( x < 6.5 ) {
	direction = 0;
	z = 1.0;
	tx = x;
	nx = 0.0;
	if(x >= 1.5) {
	       while( tx > 2.5 ){
			nx -= 1.0;
			tx = x + nx;
			z *=tx;
		}
		x += nx - 2.0;
iv1r5:
		p = x * polevlf( x, B, 7 );
		goto cont;
	}
	if( x >= 1.25 ){
		z *= x;
		x -= 1.0; /* x + 1 - 2 */
		direction = 1;
		goto iv1r5;
	}
	if( x >= 0.75 ) {
	        x -= 1.0;
		p = x * polevlf( x, C, 7 );
		q = 0.0;
		goto contz;
	}
	while( tx < 1.5 ){
	       if( tx == 0.0 )
			goto loverf;
		z *=tx;
		nx += 1.0;
		tx = x + nx;
	}
	direction = 1;
	x += nx - 2.0;
	p = x * polevlf( x, B, 7 );

cont:
	if( z < 0.0 ){
		sgngamf = -1;
		z = -z;
	}
	else
	    {
		sgngamf = 1;
	}
	q = ceph_logf(z);
	if( direction )
		q = -q;
contz:
	return( p + q );
	}
if( x > 2.035093e36){
	
loverf:
	return(sgngamf * 3.4028234663852885981170418348451692544e38);
}

/* Note, though an asymptotic formula could be used for x >= 3,
 * there is cancellation error in the following if x < 6.5.  */
q =  0.91893853320467274178 - x;
q += ( x - 0.5 ) * ceph_logf(x);
if( x <= 1.0e4 ){
	z = 1.0/x;
	p = z * z;
	q += ((    6.789774945028216E-004 * p
		 - 2.769887652139868E-003 ) * p
		+  8.333316229807355E-002 ) * z;
}
return( q );
}

/*
	Gauss hypergeometric function   F
 *	                               2 1
 *
 *
 * SYNOPSIS:
 *
 * float a, b, c, x, y, hyp2f1f();
 *
 * y = hyp2f1f( a, b, c, x );
 *
 *
 * DESCRIPTION:
 *
 *
 *  hyp2f1( a, b, c, x )  =   F ( a, b; c; x )
 *                           2 1
 *
 *           inf.
 *            -   a(a+1)...(a+k) b(b+1)...(b+k)   k+1
 *   =  1 +   >   -----------------------------  x   .
 *            -         c(c+1)...(c+k) (k+1)!
 *          k = 0
 *
 *  Cases addressed are
 *	Tests and escapes for negative integer a, b, or c
 *	Linear transformation if c - a or c - b negative integer
 *	Special case c = a or c = b
 *	Linear transformation for  x near +1
 *	Transformation for x < -0.5
 *	Psi function expansion if x > 0.5 and c - a - b integer
 *      Conditionally, a recurrence on c to make c-a-b > 0
 *
*/

#if !defined(ceph_roundf)
    #define ceph_roundf(x) (ceph_floorf((x)+(float )0.5))
#endif


__ATTR_ALWAYS_INLINE__


static inline
float hyp2f1f(const float aa,
              const float bb,
	      const float cc,
	      const float xx ){
float a, b, c, x;
float d, d1, d2, e;
float p, q, r, s, y, ax;
float ia, ib, ic, id, err;
int flag, i, aid;
a = aa;
b = bb;
c = cc;
x = xx;
err = 0.0;
ax = ceph_fabsf(x);
s = 1.0 - x;
flag = 0;
ia = ceph_roundf(a); /* nearest integer to a */
ib = ceph_roundf(b);
if( a <= 0 )
	{
	if( ceph_fabsf(a-ia) < 1.0e-5)		/* a is a negative integer */
		flag |= 1;
	}

if( b <= 0 )
	{
	if( ceph_fabsf(b-ib) < 1.0e-5)		/* b is a negative integer */
		flag |= 2;
	}

if( ax < 1.0 )
	{
	if( ceph_fabsf(b-c) < 1.0e-5)		/* b = c */
		{
		y = ceph_powf( s, -a );	/* s to the -a power */
		goto hypdon;
		}
	if( ceph_fabsf(a-c) < 1.0e-5)		/* a = c */
		{
		y = ceph_powf( s, -b );	/* s to the -b power */
		goto hypdon;
		}
	}

if( c <= 0.0 )
	{
	ic = ceph_roundf(c); 	/* nearest integer to c */
	if( ceph_fabsf(c-ic) < 1.0e-5)		/* c is a negative integer */
		{
		/* check if termination before explosion */
		if( (flag & 1) && (ia > ic) )
			goto hypok;
		if( (flag & 2) && (ib > ic) )
			goto hypok;
		goto hypdiv;
		}
	}

if( flag )			/* function is a polynomial */
	goto hypok;

if( ax > 1.0 )			/* series diverges	*/
	goto hypdiv;

p = c - a;
ia = ceph_roundf(p);
if( (ia <= 0.0) && (ceph_fabsf(p-ia) < 1.0e-5) )	/* negative int c - a */
	flag |= 4;

r = c - b;
ib = ceph_roundf(r); /* nearest integer to r */
if( (ib <= 0.0) && (ceph_fabsf(r-ib) < 1.0e-5) )	/* negative int c - b */
	flag |= 8;

d = c - a - b;
id = ceph_roundf(d); /* nearest integer to d */
q = ceph_fabsf(d-id);

if( ceph_fabsf(ax-1.0) < 1.0e-5)			/* |x| == 1.0	*/
	{
	if( x > 0.0 )
		{
		if( flag & 12 ) /* negative int c-a or c-b */
			{
			if( d >= 0.0 )
				goto hypf;
			else
				goto hypdiv;
			}
		if( d <= 0.0 )
			goto hypdiv;
		y = ceph_gammaf(c)*ceph_gammaf(d)/(ceph_gammaf(p)*ceph_gammaf(r));
		goto hypdon;
		}

	if( d <= -1.0 )
		goto hypdiv;
	}

/* Conditionally make d > 0 by recurrence on c
 * AMS55 #15.2.27
 */
if( d < 0.0 )
	{
/* Try the power series first */
	y = hyt2f1f( a, b, c, x, &err );
	if( err < 1.0e-5)
		goto hypdon;
/* Apply the recurrence if power series fails */
	err = 0.0;
	aid = 2 - id;
	e = c + aid;
	d2 = hyp2f1f(a,b,e,x);
	d1 = hyp2f1f(a,b,e+1.0,x);
	q = a + b + 1.0;
	for( i=0; i<aid; i++ )
		{
		r = e - 1.0;
		y = (e*(r-(2.0*e-q)*x)*d2 + (e-a)*(e-b)*x*d1)/(e*r*s);
		e = r;
		d1 = d2;
		d2 = y;
		}
	goto hypdon;
	}


if( flag & 12 )
	goto hypf; /* negative integer c-a or c-b */

hypok:
y = hyt2f1f( a, b, c, x, &err );

hypdon:
if( err > 1.0e-5)
	{
//	mtherr( "hyp2f1", PLOSS );
/*	printf( "Estimated err = %.2e\n", err );*/
	}
return(y);

/* The transformation for c-a or c-b negative integer
 * AMS55 #15.3.3
 */
hypf:
y = ceph_powf( s, d ) * hys2f1f( c-a, c-b, c, x, &err );
goto hypdon;

/* The alarm exit */
hypdiv:
return(3.4028234663852885981170418348451692544e38);
}

/* Apply transformations for |x| near 1
 * then call the power series
 */


__ATTR_ALWAYS_INLINE__


static inline
float hyt2f1f( const float aa,
               const float bb,
	       const float cc,
	       const float xx,
	       float *loss ) {
float a, b, c, x;
float p, q, r, s, t, y, d, err, err1;
float ax, id, d1, d2, e, y1;
int i, aid;

a = aa;
b = bb;
c = cc;
x = xx;
err = 0.0;
s = 1.0 - x;
if( x < -0.5 )
	{
	if( b > a )
		y = ceph_powf( s, -a ) * hys2f1f( a, c-b, c, -x/s, &err );

	else
		y = ceph_powf( s, -b ) * hys2f1f( c-a, b, c, -x/s, &err );

	goto done;
	}



d = c - a - b;
id = ceph_roundf(d);	/* nearest integer to d */

if( x > 0.8 )
{

if( ceph_fabsf(d-id) > 1.0e-5) /* test for integer c-a-b */
	{
/* Try the power series first */
	y = hys2f1f( a, b, c, x, &err );
	if( err < 1.0e-5)
		goto done;
/* If power series fails, then apply AMS55 #15.3.6 */
	q = hys2f1f( a, b, 1.0-d, s, &err );	
	q *= ceph_gammaf(d) /(ceph_gammaf(c-a) * ceph_gammaf(c-b));
	r = ceph_powf(s,d) * hys2f1f( c-a, c-b, d+1.0, s, &err1 );
	r *= ceph_gammaf(-d)/(ceph_gammaf(a) * ceph_gammaf(b));
	y = q + r;

	q = ceph_fabsf(q); /* estimate cancellation error */
	r = ceph_fabsf(r);
	if( q > r )
		r = q;
	err += err1 + (5.9604644775390625E-8*r)/y;

	y *= ceph_gammaf(c);
	goto done;
	}	
else
	{
/* Psi function expansion, AMS55 #15.3.10, #15.3.11, #15.3.12 */
	if( id >= 0.0 )
		{
		e = d;
		d1 = d;
		d2 = 0.0;
		aid = id;
		}
	else
		{
		e = -d;
		d1 = 0.0;
		d2 = d;
		aid = -id;
		}

	ax = ceph_logf(s);

	/* sum for t = 0 */
	y = psif(1.0) + psif(1.0+e) - psif(a+d1) - psif(b+d1) - ax;
	y /= ceph_gammaf(e+1.0);

	p = (a+d1) * (b+d1) * s / ceph_gammaf(e+2.0);	/* Poch for t=1 */
	t = 1.0;
	do
		{
		r = psif(1.0+t) + psif(1.0+t+e) - psif(a+t+d1)
			- psif(b+t+d1) - ax;
		q = p * r;
		y += q;
		p *= s * (a+t+d1) / (t+1.0);
		p *= (b+t+d1) / (t+1.0+e);
		t += 1.0;
		}
	while( fabsf(q/y) > EPS );


	if( id == 0.0 )
		{
		y *= ceph_gammaf(c)/(ceph_gammaf(a)*ceph_gammaf(b));
		goto psidon;
		}

	y1 = 1.0;

	if( aid == 1 )
		goto nosum;

	t = 0.0;
	p = 1.0;
	for( i=1; i<aid; i++ )
		{
		r = 1.0-e+t;
		p *= s * (a+t+d2) * (b+t+d2) / r;
		t += 1.0;
		p /= t;
		y1 += p;
		}


nosum:
	p = ceph_gammaf(c);
	y1 *= ceph_gammaf(e) * p / (ceph_gammaf(a+d1) * ceph_gammaf(b+d1));
	y *= p / (ceph_gammaf(a+d2) * ceph_gammaf(b+d2));
	if( (aid & 1) != 0 )
		y = -y;

	q = powf( s, id );	/* s to the id power */
	if( id > 0.0 )
		y *= q;
	else
		y1 *= q;

	y += y1;
psidon:
	goto done;
	}
}


/* Use defining power series if no special cases */
y = hys2f1f( a, b, c, x, &err );

done:
*loss = err;
return(y);
}

/* Defining power series expansion of Gauss hypergeometric function */

#if 1

__ATTR_ALWAYS_INLINE__


static inline
float hys2f1f( const float aa,
               const float bb,
	       const float cc,
	       const float xx,
	       float *loss ) {
int i;
float a, b, c, x;
float f, g, h, k, m, s, u, umax;


a = aa;
b = bb;
c = cc;
x = xx;
i = 0;
umax = 0.0;
f = a;
g = b;
h = c;
k = 0.0;
s = 1.0;
u = 1.0;

do
	{
	if( ceph_fabsf(h) < 1.0e-5)
		return(3.4028234663852885981170418348451692544e38);
	m = k + 1.0;
	u = u * ((f+k) * (g+k) * x / ((h+k) * m));
	s += u;
	k = ceph_fabsf(u);  /* remember largest term summed */
	if( k > umax )
		umax = k;
	k = m;
	if( ++i > 10000 ) /* should never happen */
		{
		*loss = 1.0;
		return(s);
		}
	}
while( ceph_fabsf(u/s) > 5.9604644775390625E-8);

/* return estimated relative error */
*loss = (5.9604644775390625E-8*umax)/ceph_fabsf(s) + (5.9604644775390625E-8*i);

return(s);
}


#else /* 0 */




__ATTR_ALWAYS_INLINE__


static inline
float hys2f1f( const float aa,
                      const float bb,
		      const float cc,
		      const float xx,
		      float *loss ) {
int i;
double a, b, c, x;
double f, g, h, k, m, s, u, umax;

a = aa;
b = bb;
c = cc;
x = xx;
i = 0;
umax = 0.0;
f = a;
g = b;
h = c;
k = 0.0;
s = 1.0;
u = 1.0;

do
	{
	if( ceph_fabsf(h) < 1.0e-5)
		{
		*loss = 1.0;
		return(3.4028234663852885981170418348451692544e38);
		}
	m = k + 1.0;
	u = u * ((f+k) * (g+k) * x / ((h+k) * m));
	s += u;
	k = ceph_fabsf(u);  /* remember largest term summed */
	if( k > umax )
		umax = k;
	k = m;
	if( ++i > 10000 ) /* should never happen */
		{
		*loss = 1.0;
		return(s);
		}
	}
while( ceph_fabsf(u/s) > 5.9604644775390625E-8);

/* return estimated relative error */
*loss = (5.9604644775390625E-8*umax)/ceph_fabsf(s) + (5.9604644775390625E-8*i);

return(s);
}
#endif

/*
   	Confluent hypergeometric function
 *
 *
 *
 * SYNOPSIS:
 *
 * float a, b, x, y, hypergf();
 *
 * y = hypergf( a, b, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes the confluent hypergeometric function
 *
 *                          1           2
 *                       a x    a(a+1) x
 *   F ( a,b;x )  =  1 + ---- + --------- + ...
 *  1 1                  b 1!   b(b+1) 2!
 *
 * Many higher transcendental functions are special cases of
 * this power series.
 *
 * As is evident from the formula, b must not be a negative
 * integer or zero unless a is an integer with 0 >= a > b.
 *
 * The routine attempts both a direct summation of the series
 * and an asymptotic expansion.  In each case error due to
 * roundoff, cancellation, and nonconvergence is estimated.
 * The result with smaller estimated error is returned.
*/

 
__ATTR_ALWAYS_INLINE__


static inline
float hypergf( float aa, float bb, float xx ){
float a, b, x, asum, psum, acanc, pcanc, temp;
a = aa;
b = bb;
x = xx;
/* See if a Kummer transformation will help */
temp = b - a;
if( ceph_fabsf(temp) < 0.001 * ceph_fabsf(a) )
	return( ceph_expf(x) * hypergf( temp, b, -x )  );

psum = hy1f1pf( a, b, x, &pcanc );
if( pcanc < 1.0e-6 )
	goto done;


/* try asymptotic series */

asum = hy1f1af( a, b, x, &acanc );


/* Pick the result with less estimated error */

if( acanc < pcanc )
	{
	pcanc = acanc;
	psum = asum;
	}

done:
if( pcanc > 1.0e-3 )
    ;

return( psum );
}

/* Power series summation for confluent hypergeometric function		*/



__ATTR_ALWAYS_INLINE__


static inline
float hy1f1pf( const float aa,
               const float bb,
	       const float xx,
	       float *err ) {
float a, b, x, n, a0, sum, t, u, temp;
float an, bn, maxt, pcanc;
a = aa;
b = bb;
x = xx;
/* set up for power series summation */
an = a;
bn = b;
a0 = 1.0;
sum = 1.0;
n = 1.0;
t = 1.0;
maxt = 0.0;

while( t > 5.9604644775390625E-8)
	{
	if( bn == 0 )			/* check bn first since if both	*/
		{
	
		return(3.4028234663852885981170418348451692544e38 );	/* an and bn are zero it is	*/
		}
	if( an == 0 )			/* a singularity		*/
		return( sum );
	if( n > 200 )
		goto pdone;
	u = x * ( an / (bn * n) );

	/* check for blowup */
	temp = ceph_fabsf(u);
	if( (temp > 1.0 ) && (maxt > (3.4028234663852885981170418348451692544e38/temp)) )
		{
		pcanc = 1.0;	/* estimate 100% error */
		goto blowup;
		}

	a0 *= u;
	sum += a0;
	t = ceph_fabsf(a0);
	if( t > maxt )
		maxt = t;
/*
	if( (maxt/fabsf(sum)) > 1.0e17 )
		{
		pcanc = 1.0;
		goto blowup;
		}
*/
	an += 1.0;
	bn += 1.0;
	n += 1.0;
	}

pdone:

/* estimate error due to roundoff and cancellation */
if( sum != 0.0 )
	maxt /= ceph_fabsf(sum);
maxt *= 5.9604644775390625E-8; 	/* this way avoids multiply overflow */
pcanc = ceph_fabsf( 5.9604644775390625E-8* n  +  maxt );

blowup:

*err = pcanc;

return( sum );
}

/*							hy1f1a()	*/
/* asymptotic formula for hypergeometric function:
 *
 *        (    -a                         
 *  --    ( |z|                           
 * |  (b) ( -------- 2f0( a, 1+a-b, -1/x )
 *        (  --                           
 *        ( |  (b-a)                      
 *
 *
 *                                x    a-b                     )
 *                               e  |x|                        )
 *                             + -------- 2f0( b-a, 1-a, 1/x ) )
 *                                --                           )
 *                               |  (a)                        )
 */


__ATTR_ALWAYS_INLINE__


static inline
float hy1f1af( const
                      const float aa,
		      const float bb,
		      const float xx,
		      float *err ) {
float a, b, x, h1, h2, t, u, temp, acanc, asum, err1, err2;

a = aa;
b = bb;
x = xx;
if( x == 0 )
	{
	acanc = 1.0;
	asum = 3.4028234663852885981170418348451692544e38;
	goto adone;
	}
temp = ceph_logf( ceph_fabsf(x) );
t = x + temp * (a-b);
u = -temp * a;

if( b > 0 )
	{
	temp = ceph_lgamf(b);
	t += temp;
	u += temp;
	}

h1 = hyp2f0f( a, a-b+1, -1.0/x, 1, &err1 );

temp = ceph_expf(u) / ceph_gammaf(b-a);
h1 *= temp;
err1 *= temp;

h2 = hyp2f0f( b-a, 1.0-a, 1.0/x, 2, &err2 );

if( a < 0 )
	temp = ceph_expf(t) / ceph_gammaf(a);
else
	temp = ceph_expf( t - ceph_lgamf(a) );

h2 *= temp;
err2 *= temp;

if( x < 0.0 )
	asum = h1;
else
	asum = h2;

acanc = ceph_fabsf(err1) + ceph_fabsf(err2);


if( b < 0 )
	{
	temp = ceph_gammaf(b);
	asum *= temp;
	acanc *= ceph_fabsf(temp);
	}


if( asum != 0.0 )
	acanc /= ceph_fabsf(asum);

acanc *= 30.0;	/* fudge factor, since error of asymptotic formula
		 * often seems this much larger than advertised */

adone:


*err = acanc;
return( asum );
}

__ATTR_ALWAYS_INLINE__


static inline
float hyp2f0f(const float aa,
              const float bb,
	      const float xx,
	      const int type,
	      const float *err){
float a, b, x, a0, alast, t, tlast, maxt;
float n, an, bn, u, sum, temp;

a = aa;
b = bb;
x = xx;
an = a;
bn = b;
a0 = 1.0;
alast = 1.0;
sum = 0.0;
n = 1.0;
t = 1.0;
tlast = 1.0e9;
maxt = 0.0;

do
	{
	if( an == 0 )
		goto pdone;
	if( bn == 0 )
		goto pdone;

	u = an * (bn * x / n);

	/* check for blowup */
	temp = ceph_fabsf(u);
	if( (temp > 1.0 ) && (maxt > (3.4028234663852885981170418348451692544e38/temp)) )
		goto error;

	a0 *= u;
	t = ceph_fabsf(a0);

	/* terminating condition for asymptotic series */
	if( t > tlast )
		goto ndone;

	tlast = t;
	sum += alast;	/* the sum is one term behind */
	alast = a0;

	if( n > 200 )
		goto ndone;

	an += 1.0;
	bn += 1.0;
	n += 1.0;
	if( t > maxt )
		maxt = t;
	}
while( t > 5.9604644775390625E-8);


pdone:	/* series converged! */

/* estimate error due to roundoff and cancellation */
*err = ceph_fabsf(  5.9604644775390625E-8 * (n + maxt)  );

alast = a0;
goto done;

ndone:	/* series did not converge */

/* The following "Converging factors" are supposed to improve accuracy,
 * but do not actually seem to accomplish very much. */

n -= 1.0;
x = 1.0/x;

switch( type )	/* "type" given as subroutine argument */
{
case 1:
	alast *= ( 0.5 + (0.125 + 0.25*b - 0.5*a + 0.25*x - 0.25*n)/x );
	break;

case 2:
	alast *= 2.0/3.0 - b + 2.0*a + x - n;
	break;

default:
	;
}

/* estimate error due to roundoff, cancellation, and nonconvergence */
*err = 5.9604644775390625E-8* (n + maxt)  +  ceph_fabsf( a0 );


done:
sum += alast;
return( sum );

/* series blew up: */
error:
*err = 3.4028234663852885981170418348451692544e38;
return( sum );
}

/*
    	Modified Bessel function of order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, i0();
 *
 * y = i0f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of order zero of the
 * argument.
 *
 * The function is defined as i0(x) = j0( ix ).
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
*/


__ATTR_ALWAYS_INLINE__


static inline
float i0f( float x ) {
__ATTR_ALIGN__(64) const float A[] =
{
-1.30002500998624804212E-8f,
 6.04699502254191894932E-8f,
-2.67079385394061173391E-7f,
 1.11738753912010371815E-6f,
-4.41673835845875056359E-6f,
 1.64484480707288970893E-5f,
-5.75419501008210370398E-5f,
 1.88502885095841655729E-4f,
-5.76375574538582365885E-4f,
 1.63947561694133579842E-3f,
-4.32430999505057594430E-3f,
 1.05464603945949983183E-2f,
-2.37374148058994688156E-2f,
 4.93052842396707084878E-2f,
-9.49010970480476444210E-2f,
 1.71620901522208775349E-1f,
-3.04682672343198398683E-1f,
 6.76795274409476084995E-1f
};


/* Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
 */

 const float B[] =
{
 3.39623202570838634515E-9f,
 2.26666899049817806459E-8f,
 2.04891858946906374183E-7f,
 2.89137052083475648297E-6f,
 6.88975834691682398426E-5f,
 3.36911647825569408990E-3f,
 8.04490411014108831608E-1f
};
float y;

if( x < 0 )
	x = -x;
if( x <= 8.0f )
	{
	y = 0.5f*x - 2.0f;
	return( ceph_expf(x) * chbevlf( y, A, 18 ) );
	}

return(  ceph_expf(x) * chbevlf( 32.0f/x - 2.0f, B, 7 ) / ceph_sqrtf(x) );
}


__ATTR_ALWAYS_INLINE__


static inline
float i0ef( float x ) {
__ATTR_ALIGN__(64) const float A[] =
{
-1.30002500998624804212E-8f,
 6.04699502254191894932E-8f,
-2.67079385394061173391E-7f,
 1.11738753912010371815E-6f,
-4.41673835845875056359E-6f,
 1.64484480707288970893E-5f,
-5.75419501008210370398E-5f,
 1.88502885095841655729E-4f,
-5.76375574538582365885E-4f,
 1.63947561694133579842E-3f,
-4.32430999505057594430E-3f,
 1.05464603945949983183E-2f,
-2.37374148058994688156E-2f,
 4.93052842396707084878E-2f,
-9.49010970480476444210E-2f,
 1.71620901522208775349E-1f,
-3.04682672343198398683E-1f,
 6.76795274409476084995E-1f
};


/* Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
 */

 const float B[] =
{
 3.39623202570838634515E-9f,
 2.26666899049817806459E-8f,
 2.04891858946906374183E-7f,
 2.89137052083475648297E-6f,
 6.88975834691682398426E-5f,
 3.36911647825569408990E-3f,
 8.04490411014108831608E-1f
};
float y;
if( x < 0 )
	x = -x;
if( x <= 8.0f )
	{
	y = 0.5f*x - 2.0f;
	return( chbevlf( y, A, 18 ) );
	}

return(  chbevlf( 32.0f/x - 2.0f, B, 7 ) / ceph_sqrtf(x) );
}

/*
   	Modified Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, i1f();
 *
 * y = i1f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of order one of the
 * argument.
 *
 * The function is defined as i1(x) = -i j1( ix ).
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
*/


__ATTR_ALWAYS_INLINE__


static inline
float i1f(const float xx) {
/* Chebyshev coefficients for exp(-x) I1(x) / x
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I1(x) / x } = 1/2.
 */

__ATTR_ALIGN__(64) const float A[] =
{
 9.38153738649577178388E-9f,
-4.44505912879632808065E-8f,
 2.00329475355213526229E-7f,
-8.56872026469545474066E-7f,
 3.47025130813767847674E-6f,
-1.32731636560394358279E-5f,
 4.78156510755005422638E-5f,
-1.61760815825896745588E-4f,
 5.12285956168575772895E-4f,
-1.51357245063125314899E-3f,
 4.15642294431288815669E-3f,
-1.05640848946261981558E-2f,
 2.47264490306265168283E-2f,
-5.29459812080949914269E-2f,
 1.02643658689847095384E-1f,
-1.76416518357834055153E-1f,
 2.52587186443633654823E-1f
};


/* Chebyshev coefficients for exp(-x) sqrt(x) I1(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I1(x) } = 1/sqrt(2pi).
 */

  const float B[] =
{
-3.83538038596423702205E-9f,
-2.63146884688951950684E-8f,
-2.51223623787020892529E-7f,
-3.88256480887769039346E-6f,
-1.10588938762623716291E-4f,
-9.76109749136146840777E-3f,
 7.78576235018280120474E-1f,
 0.0F
};
float x, y, z;
x = xx;
z = ceph_fabsf(x);
if( z <= 8.0f )
	{
	y = 0.5f*z - 2.0f;
	z = chbevlf( y, A, 17 ) * z * ceph_expf(z);
	}
else
	{
	z = ceph_expf(z) * chbevlf( 32.0f/z - 2.0f, B, 7 ) / ceph_sqrtf(z);
	}
if( x < 0.0f )
	z = -z;
return( z );
}


__ATTR_ALWAYS_INLINE__


static inline
float i1ef( const float xx ) {
float x, y, z;
__ATTR_ALIGN__(64) const float A[] =
{
 9.38153738649577178388E-9f,
-4.44505912879632808065E-8f,
 2.00329475355213526229E-7f,
-8.56872026469545474066E-7f,
 3.47025130813767847674E-6f,
-1.32731636560394358279E-5f,
 4.78156510755005422638E-5f,
-1.61760815825896745588E-4f,
 5.12285956168575772895E-4f,
-1.51357245063125314899E-3f,
 4.15642294431288815669E-3f,
-1.05640848946261981558E-2f,
 2.47264490306265168283E-2f,
-5.29459812080949914269E-2f,
 1.02643658689847095384E-1f,
-1.76416518357834055153E-1f,
 2.52587186443633654823E-1f
};


/* Chebyshev coefficients for exp(-x) sqrt(x) I1(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I1(x) } = 1/sqrt(2pi).
 */

  const float B[] =
{
-3.83538038596423702205E-9f,
-2.63146884688951950684E-8f,
-2.51223623787020892529E-7f,
-3.88256480887769039346E-6f,
-1.10588938762623716291E-4f,
-9.76109749136146840777E-3f,
 7.78576235018280120474E-1f,
 0.0F
};
x = xx;
z = ceph_fabsf(x);
if( z <= 8.0f )
	{
	y = 0.5f*z - 2.0f;
	z = chbevlf( y, A, 17 ) * z;
	}
else
	{
	z = chbevlf( 32.0f/z - 2.0f, B, 7 ) / sqrtf(z);
	}
if( x < 0.0f )
	z = -z;
return( z );
}

/*
   	Incomplete gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * float a, x, y, igamf();
 *
 * y = igamf( a, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * The function is defined by
 *
 *                           x
 *                            -
 *                   1       | |  -t  a-1
 *  igam(a,x)  =   -----     |   e   t   dt.
 *                  -      | |
 *                 | (a)    -
 *                           0
 *
 *
 * In this implementation both arguments must be positive.
 * The integral is evaluated by either a power series or
 * continued fraction expansion, depending on the relative
 * values of a and x.
*/


__ATTR_ALWAYS_INLINE__


static inline
float igamcf( const float aa, const float xx ) {
float a, x, ans, c, yc, ax, y, z;
float pk, pkm1, pkm2, qk, qkm1, qkm2;
float r, t;
const float big =  16777216.0f;
a = aa;
x = xx;
if( (x <= 0) || ( a <= 0) )
	return( 1.0 );

if( (x < 1.0) || (x < a) )
	return( 1.0 - ceph_igamf(a,x) );

ax = a * ceph_logf(x) - x - ceph_lgamf(a);
if( ax < -88.72283905206835 )
	{

	return( 0.0 );
	}
ax = ceph_expf(ax);

/* continued fraction */
y = 1.0 - a;
z = x + y + 1.0;
c = 0.0;
pkm2 = 1.0;
qkm2 = x;
pkm1 = x + 1.0;
qkm1 = z * x;
ans = pkm1/qkm1;

do
	{
	c += 1.0;
	y += 1.0;
	z += 2.0;
	yc = y * c;
	pk = pkm1 * z  -  pkm2 * yc;
	qk = qkm1 * z  -  qkm2 * yc;
	if( qk != 0 )
		{
		r = pk/qk;
		t = ceph_fabsf( (ans - r)/r );
		ans = r;
		}
	else
		t = 1.0;
	pkm2 = pkm1;
	pkm1 = pk;
	qkm2 = qkm1;
	qkm1 = qk;
	if( eph_fabsf(pk) > big )
		{
		pkm2 *= 5.9604644775390625E-8;
		pkm1 *= 5.9604644775390625E-8;
		qkm2 *= 5.9604644775390625E-8;
		qkm1 *= 5.9604644775390625E-8;
		}
	}
while( t > 5.9604644775390625E-8);

return( ans * ax );
}

/* left tail of incomplete gamma function:
 *
 *          inf.      k
 *   a  -x   -       x
 *  x  e     >   ----------
 *           -     -
 *          k=0   | (a+k+1)
 *
 */


__ATTR_ALWAYS_INLINE__


static inline
float igamf( const float aa, const float xx ){
float a, x, ans, ax, c, r;
a = aa;
x = xx;
if( (x <= 0) || ( a <= 0) )
	return( 0.0 );

if( (x > 1.0) && (x > a ) )
	return( 1.0 - ceph_igamcf(a,x) );

/* Compute  x**a * exp(-x) / gamma(a)  */
ax = a * ceph_logf(x) - x - ceph_lgamf(a);
if( ax < -88.72283905206835;)
	{

	return( 0.0 );
	}
ax = ceph_expf(ax);

/* power series */
r = a;
c = 1.0;
ans = 1.0;

do
	{
	r += 1.0;
	c *= x/r;
	ans += c;
	}
while( c/ans > 5.9604644775390625E-8 );

return( ans * ax/a );
}


__ATTR_ALWAYS_INLINE__


static inline
float igamif( const float aa, const float yy0 ){
float a, y0, d, y, x0, lgm;
int i;

//if( yy0 > 0.5)
//	mtherr( "igamif", PLOSS );
a = aa;
y0 = yy0;
/* approximation to inverse function */
d = 1.0/(9.0*a);
y = ( 1.0 - d - ndtrif(y0) * ceph_sqrtf(d) );
x0 = a * y * y * y;
lgm = ceph_lgamf(a);
for( i=0; i<10; i++ )
	{
	if( x0 <= 0.0 )
		{
	          return(0.0);
		}
	y = igamcf(a,x0);
/* compute the derivative of the function at this point */
	d = (a - 1.0) * ceph_logf(x0) - x0 - lgm;
	if( d < -88.72283905206835 )
		{
	           goto done;
		}
	d = -ceph_expf(d);
/* compute the step to the next approximation of x */
	if( d == 0.0 )
		goto done;
	d = (y - y0)/d;
	x0 = x0 - d;
	if( i < 3 )
		continue;
	if( ceph_fabsf(d/x0) < (2.0 * 5.9604644775390625E-8 ) )
		goto done;
	}

done:
return( x0 );
}

/*
   	Modified Bessel function of noninteger order
 *
 *
 *
 * SYNOPSIS:
 *
 * float v, x, y, ivf();
 *
 * y = ivf( v, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of order v of the
 * argument.  If x is negative, v must be integer valued.
 *
 * The function is defined as Iv(x) = Jv( ix ).  It is
 * here computed in terms of the confluent hypergeometric
 * function, according to the formula
 *
 *              v  -x
 * Iv(x) = (x/2)  e   hyperg( v+0.5, 2v+1, 2x ) / gamma(v+1)
 *
 * If v is a negative integer, then v is replaced by -v.
*/

__ATTR_ALWAYS_INLINE__


static inline
float ivf( const float v, float x ){
int sign;
float t, ax;
/* If v is a negative integer, invoke symmetry */
t = ceph_floorf(v);
if( v < 0.0 )
	{
	if( t == v )
		{
		v = -v;	/* symmetry */
		t = -t;
		}
	}
/* If x is negative, require v to be an integer */
sign = 1;
if( x < 0.0 )
	{
	if( t != v )
		{
	 //	mtherr( "ivf", DOMAIN );
		return( 0.0 );
		}
	if( v != 2.0 * ceph_floorf(v/2.0) )
		sign = -1;
	}

/* Avoid logarithm singularity */
if( x == 0.0 )
	{
	if( v == 0.0 )
		return( 1.0 );
	if( v < 0.0 )
		{
		//mtherr( "ivf", OVERFLOW );
		return( MAXNUMF );
		}
	else
		return( 0.0 );
	}

ax = ceph_fabsf(x);
t = v * ceph_logf( 0.5 * ax )  -  x;
t = sign * ceph_expf(t) / ceph_gammaf( v + 1.0 );
ax = v + 0.5;
return( t * hypergf( ax,  2.0 * ax,  2.0 * x ) );
}

/*
   	Bessel function of order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, j0f();
 *
 * y = j0f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order zero of the argument.
 *
 * The domain is divided into the intervals [0, 2] and
 * (2, infinity). In the first interval the following polynomial
 * approximation is used:
 *
 *
 *        2         2         2
 * (w - r  ) (w - r  ) (w - r  ) P(w)
 *       1         2         3   
 *
*/

__ATTR_ALWAYS_INLINE__


static inline
float j0f( const float xx ){
 const float MO[8] = {
-6.838999669318810E-002f,
 1.864949361379502E-001f,
-2.145007480346739E-001f,
 1.197549369473540E-001f,
-3.560281861530129E-003f,
-4.969382655296620E-002f,
-3.355424622293709E-006f,
 7.978845717621440E-001f
};
 const float PH[8] = {
 3.242077816988247E+001f,
-3.630592630518434E+001f,
 1.756221482109099E+001f,
-4.974978466280903E+000f,
 1.001973420681837E+000f,
-1.939906941791308E-001f,
 6.490598792654666E-002f,
-1.249992184872738E-001f
};
 const float YP[8] = {
 9.454583683980369E-008f,
-9.413212653797057E-006f,
 5.344486707214273E-004f,
-1.584289289821316E-002f,
 1.707584643733568E-001f,
 0.0e+00, // padded to 32-bytes
 0.0e+00,
 0.0e+00
};
_attr_align__(32) const float JP[8] = {
-6.068350350393235E-008f,
 6.388945720783375E-006f,
-3.969646342510940E-004f,
 1.332913422519003E-002f,
-1.729150680240724E-001f,
 0.0e+00, // padded to 32-bytes
 0.0e+00,
 0.0e+00
};
constexpr float YZ1 =  0.43221455686510834878f;
constexpr float YZ2 = 22.401876406482861405f;
constexpr float YZ3 = 64.130620282338755553f;
constexpr float DR1 =  5.78318596294678452118f;
float x, w, z, p, q, xn;
if( xx < 0 )
	x = -xx;
else
	x = xx;

if( x <= 2.0f )
	{
	z = x * x;
	if( x < 1.0e-3f )
		return( 1.0f - 0.25f*z );

	p = (z-DR1) * polevlf( z, JP, 4);
	return( p );
	}

q = 1.0f/x;
w = ceph_sqrtf(q);
p = w * polevlf( q, MO, 7);
w = q*q;
xn = q * polevlf( w, PH, 7) - PIO4F;
p = p * ceph_cosf(xn + x);
return(p);
}

/*							y0() 2	*/
/* Bessel function of second kind, order zero	*/

/* Rational approximation coefficients YP[] are used for x < 6.5.
 * The function computed is  y0(x)  -  2 ln(x) j0(x) / pi,
 * whose value at x = 0 is  2 * ( log(0.5) + EUL ) / pi
 * = 0.073804295108687225 , EUL is Euler's constant.
 */

__ATTR_ALWAYS_INLINE__


static inline
float y0f( const float xx ) {
 const float MO[8] = {
-6.838999669318810E-002f,
 1.864949361379502E-001f,
-2.145007480346739E-001f,
 1.197549369473540E-001f,
-3.560281861530129E-003f,
-4.969382655296620E-002f,
-3.355424622293709E-006f,
 7.978845717621440E-001f
};
 const float PH[8] = {
 3.242077816988247E+001f,
-3.630592630518434E+001f,
 1.756221482109099E+001f,
-4.974978466280903E+000f,
 1.001973420681837E+000f,
-1.939906941791308E-001f,
 6.490598792654666E-002f,
-1.249992184872738E-001f
};
 const float YP[8] = {
 9.454583683980369E-008f,
-9.413212653797057E-006f,
 5.344486707214273E-004f,
-1.584289289821316E-002f,
 1.707584643733568E-001f,
 0.0e+00, // padded to 32-bytes
 0.0e+00,
 0.0e+00
};
constexpr float TWOOPI =  0.636619772367581343075535f;
float x, w, z, p, q, xn;
x = xx;
if( x <= 2.0f )
	{
	if( x <= 0.0f )
		{
		//mtherr( "y0f", DOMAIN );
		return( -MAXNUMF );
		}
	z = x * x;
/*	w = (z-YZ1)*(z-YZ2)*(z-YZ3) * polevlf( z, YP, 4);*/
	w = (z-YZ1) * polevlf( z, YP, 4);
	w += TWOOPI * ceph_logf(x) * j0f(x);
	return( w );
	}

q = 1.0f/x;
w = ceph_sqrtf(q);
p = w * polevlf( q, MO, 7);
w = q*q;
xn = q * polevlf( w, PH, 7) - PIO4F;
p = p * ceph_sinf(xn + x);
return( p );
}

/*
   	Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, j1f();
 *
 * y = j1f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order one of the argument.
 *
 * The domain is divided into the intervals [0, 2] and
 * (2, infinity). In the first interval a polynomial approximation
 *        2 
 * (w - r  ) x P(w)
 *       1  
 *                     2 
 * is used, where w = x  and r is the first zero of the function.
 *
 * In the second interval, the modulus and phase are approximated
 * by polynomials of the form Modulus(x) = sqrt(1/x) Q(1/x)
 * and Phase(x) = x + 1/x R(1/x^2) - 3pi/4.  The function is
 *
 *   j0(x) = Modulus(x) cos( Phase(x) ).
*/


__ATTR_ALWAYS_INLINE__


static inline
float j1f( const float xx ) {
 const float JP[8] = {
-4.878788132172128E-009f,
 6.009061827883699E-007f,
-4.541343896997497E-005f,
 1.937383947804541E-003f,
-3.405537384615824E-002f,
 0.0f,
 0.0f,
 0.0f
};
_ATTR_ALIGN__(32) const float MO1[8] = {
 6.913942741265801E-002f,
-2.284801500053359E-001f,
 3.138238455499697E-001f,
-2.102302420403875E-001f,
 5.435364690523026E-003f,
 1.493389585089498E-001f,
 4.976029650847191E-006f,
 7.978845453073848E-001f
};
 const float PH1[8] = {
-4.497014141919556E+001f,
 5.073465654089319E+001f,
-2.485774108720340E+001f,
 7.222973196770240E+000f,
-1.544842782180211E+000f,
 3.503787691653334E-001f,
-1.637986776941202E-001f,
 3.749989509080821E-001f
};
constexpr float Z1 = 1.46819706421238932572E1f;
constexpr float THPIO4F =  2.35619449019234492885f;
float x, w, z, p, q, xn;
x = xx;
if( x < 0 )
	x = -xx;

if( x <= 2.0f )
	{
	z = x * x;	
	p = (z-Z1) * x * polevlf( z, JP, 4 );
	return( p );
	}

q = 1.0f/x;
w = ceph_sqrtf(q);

p = w * polevlf( q, MO1, 7);
w = q*q;
xn = q * polevlf( w, PH1, 7) - THPIO4F;
p = p * ceph_cosf(xn + x);
return(p);
}


__ATTR_ALWAYS_INLINE__


static inline
float y1f( const float xx ){
 const float YP[8] = {
 8.061978323326852E-009f,
-9.496460629917016E-007f,
 6.719543806674249E-005f,
-2.641785726447862E-003f,
 4.202369946500099E-002f,
 0.0f,
 0.0f,
 0.0f
};
 const float MO1[8] = {
 6.913942741265801E-002f,
-2.284801500053359E-001f,
 3.138238455499697E-001f,
-2.102302420403875E-001f,
 5.435364690523026E-003f,
 1.493389585089498E-001f,
 4.976029650847191E-006f,
 7.978845453073848E-001f
};
 const float PH1[8] = {
-4.497014141919556E+001f,
 5.073465654089319E+001f,
-2.485774108720340E+001f,
 7.222973196770240E+000f,
-1.544842782180211E+000f,
 3.503787691653334E-001f,
-1.637986776941202E-001f,
 3.749989509080821E-001f
};
constexpr  float YO1     =  4.66539330185668857532f;
constexpr  float TWOOPI  =  0.636619772367581343075535f; /* 2/pi */
constexpr  float THPIO4F =  2.35619449019234492885f;    /* 3*pi/4 */
float x, w, z, p, q, xn;

x = xx;
if( x <= 2.0f )
	{
	if( x <= 0.0f )
		{
		//mtherr( "y1f", DOMAIN );
		return( -MAXNUMF );
		}
	z = x * x;
	w = (z - YO1) * x * polevlf( z, YP, 4 );
	w += TWOOPI * ( j1f(x) * ceph_logf(x)  -  1.0f/x );
	return( w );
	}

q = 1.0f/x;
w = ceph_sqrtf(q);

p = w * polevlf( q, MO1, 7);
w = q*q;
xn = q * polevlf( w, PH1, 7) - THPIO4F;
p = p * ceph_sinf(xn + x);
return(p);
}

/*
    	Bessel function of integer order
 *
 *
 *
 * SYNOPSIS:
 *
 * int n;
 * float x, y, jnf();
 *
 * y = jnf( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order n, where n is a
 * (possibly negative) integer.
 *
 * The ratio of jn(x) to j0(x) is computed by backward
 * recurrence.  First the ratio jn/jn-1 is found by a
 * continued fraction expansion.  Then the recurrence
 * relating successive orders is applied until j0 or j1 is
 * reached.
 *
 * If n = 0 or 1 the routine for j0 or j1 is called
 * directly.
 *
*/


__ATTR_ALWAYS_INLINE__


static inline
float jnf( const int n, const float xx ) {
float x, pkm2, pkm1, pk, xk, r, ans, xinv, sign;
int k;
x = xx;
sign = 1.0;
if( n < 0 )
	{
	n = -n;
	if( (n & 1) != 0 )	/* -1**n */
		sign = -1.0;
	}

if( n == 0 )
	return( sign * j0f(x) );
if( n == 1 )
	return( sign * j1f(x) );
if( n == 2 )
	return( sign * (2.0 * j1f(x) / x  -  j0f(x)) );

/*
if( x < MACHEPF )
	return( 0.0 );
*/

/* continued fraction */
k = 24;
pk = 2 * (n + k);
ans = pk;
xk = x * x;

do
	{
	pk -= 2.0;
	ans = pk - (xk/ans);
	}
while( --k > 0 );
/*ans = x/ans;*/

/* backward recurrence */

pk = 1.0;
/*pkm1 = 1.0/ans;*/
xinv = 1.0/x;
pkm1 = ans * xinv;
k = n-1;
r = (float )(2 * k);

do
	{
	pkm2 = (pkm1 * r  -  pk * x) * xinv;
	pk = pkm1;
	pkm1 = pkm2;
	r -= 2.0;
	}
while( --k > 0 );

r = pk;
if( r < 0 )
	r = -r;
ans = pkm1;
if( ans < 0 )
	ans = -ans;

if( r > ans )  /* if( fabs(pk) > fabs(pkm1) ) */
	ans = sign * j1f(x)/pk;
else
	ans = sign * j0f(x)/pkm1;
return( ans );
}

/*
   	Modified Bessel function, third kind, order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, k0f();
 *
 * y = k0f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of the third kind
 * of order zero of the argument.
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
*/


__ATTR_ALWAYS_INLINE__


static inline
float k0f( const float xx )
/* Chebyshev coefficients for K0(x) + log(x/2) I0(x)
 * in the interval [0,2].  The odd order coefficients are all
 * zero; only the even order coefficients are listed.
 * 
 * lim(x->0){ K0(x) + log(x/2) I0(x) } = -EUL.
 */

 const float A[] =
{
 1.90451637722020886025E-9f,
 2.53479107902614945675E-7f,
 2.28621210311945178607E-5f,
 1.26461541144692592338E-3f,
 3.59799365153615016266E-2f,
 3.44289899924628486886E-1f,
-5.35327393233902768720E-1f,
 0.0f
};

/* Chebyshev coefficients for exp(x) sqrt(x) K0(x)
 * in the inverted interval [2,infinity].
 * 
 * lim(x->inf){ exp(x) sqrt(x) K0(x) } = sqrt(pi/2).
 */

__ATTR_ALIGN__(64) const float B[] = {
-1.69753450938905987466E-9f,
 8.57403401741422608519E-9f,
-4.66048989768794782956E-8f,
 2.76681363944501510342E-7f,
-1.83175552271911948767E-6f,
 1.39498137188764993662E-5f,
-1.28495495816278026384E-4f,
 1.56988388573005337491E-3f,
-3.14481013119645005427E-2f,
 2.44030308206595545468E0f,
 0.0f, // padding 64-bytes
 0.0f,
 0.0f,
 0.0f,
 0.0f,
 0.0f
 };
float x, y, z;
x = xx;
if( x <= 0.0f )
	{
	//mtherr( "k0f", DOMAIN );
	return( MAXNUMF );
	}

if( x <= 2.0f )
	{
	y = x * x - 2.0f;
	y = chbevlf( y, A, 7 ) - ceph_logf( 0.5f * x ) * i0f(x);
	return( y );
	}
z = 8.0f/x - 2.0f;
y = ceph_expf(-x) * chbevlf( z, B, 10 ) / ceph_sqrtf(x);
return(y);
}


__ATTR_ALWAYS_INLINE__


static inline
float k0ef( const float xx ){
/* Chebyshev coefficients for K0(x) + log(x/2) I0(x)
 * in the interval [0,2].  The odd order coefficients are all
 * zero; only the even order coefficients are listed.
 * 
 * lim(x->0){ K0(x) + log(x/2) I0(x) } = -EUL.
 */

 const float A[] =
{
 1.90451637722020886025E-9f,
 2.53479107902614945675E-7f,
 2.28621210311945178607E-5f,
 1.26461541144692592338E-3f,
 3.59799365153615016266E-2f,
 3.44289899924628486886E-1f,
-5.35327393233902768720E-1f,
 0.0f
};

/* Chebyshev coefficients for exp(x) sqrt(x) K0(x)
 * in the inverted interval [2,infinity].
 * 
 * lim(x->inf){ exp(x) sqrt(x) K0(x) } = sqrt(pi/2).
 */

__ATTR_ALIGN__(64) const float B[] = {
-1.69753450938905987466E-9f,
 8.57403401741422608519E-9f,
-4.66048989768794782956E-8f,
 2.76681363944501510342E-7f,
-1.83175552271911948767E-6f,
 1.39498137188764993662E-5f,
-1.28495495816278026384E-4f,
 1.56988388573005337491E-3f,
-3.14481013119645005427E-2f,
 2.44030308206595545468E0f,
 0.0f, // padding to 64-bytes
 0.0f,
 0.0f,
 0.0f,
 0.0f,
 0.0f
 };
float x, y;
x = xx;
if( x <= 0.0f )
	{
	//mtherr( "k0ef", DOMAIN );
	return( MAXNUMF );
	}

if( x <= 2.0f )
	{
	y = x * x - 2.0f;
	y = chbevlf( y, A, 7 ) - ceph_logf( 0.5f * x ) * i0f(x);
	return( y * ceph_expf(x) );
	}

y = chbevlf( 8.0f/x - 2.0f, B, 10 ) / ceph_sqrtf(x);
return(y);
}

/*
   	Modified Bessel function, third kind, order one
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, k1f();
 *
 * y = k1f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes the modified Bessel function of the third kind
 * of order one of the argument.
 *
 * The range is partitioned into the two intervals [0,2] and
 * (2, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
*/

#define MINNUMF 6.0e-39


__ATTR_ALWAYS_INLINE__


static inline
float k1f(const float xx) {
 const float A[] =
{
-2.21338763073472585583E-8f,
-2.43340614156596823496E-6f,
-1.73028895751305206302E-4f,
-6.97572385963986435018E-3f,
-1.22611180822657148235E-1f,
-3.53155960776544875667E-1f,
 1.52530022733894777053E0f,
 0.0
};

/* Chebyshev coefficients for exp(x) sqrt(x) K1(x)
 * in the interval [2,infinity].
 *
 * lim(x->inf){ exp(x) sqrt(x) K1(x) } = sqrt(pi/2).
 */

 const float B[] =
{
 2.01504975519703286596E-9f,
-1.03457624656780970260E-8f,
 5.74108412545004946722E-8f,
-3.50196060308781257119E-7f,
 2.40648494783721712015E-6f,
-1.93619797416608296024E-5f,
 1.95215518471351631108E-4f,
-2.85781685962277938680E-3f,
 1.03923736576817238437E-1f,
 2.72062619048444266945E0f,
 0.0f,
 0.0f,
 0.0f,
 0.0f,
 0.0f,
 0.0f
};
float x, y;
x = xx;
if( x <= MINNUMF )
	{
	 //mtherr( "k1f", DOMAIN );
	return( MAXNUMF );
	}

if( x <= 2.0f )
	{
	y = x * x - 2.0f;
	y =  ceph_logf( 0.5f * x ) * i1f(x)  +  chbevlf( y, A, 7 ) / x;
	return( y );
	}

return(  ceph_expf(-x) * chbevlf( 8.0f/x - 2.0f, B, 10 ) / ceph_sqrtf(x) );

}


__ATTR_ALWAYS_INLINE__


static inline
float k1ef( const float xx ) {
 const float A[] =
{
-2.21338763073472585583E-8f,
-2.43340614156596823496E-6f,
-1.73028895751305206302E-4f,
-6.97572385963986435018E-3f,
-1.22611180822657148235E-1f,
-3.53155960776544875667E-1f,
 1.52530022733894777053E0f,
 0.0
};

/* Chebyshev coefficients for exp(x) sqrt(x) K1(x)
 * in the interval [2,infinity].
 *
 * lim(x->inf){ exp(x) sqrt(x) K1(x) } = sqrt(pi/2).
 */

 const float B[] =
{
 2.01504975519703286596E-9f,
-1.03457624656780970260E-8f,
 5.74108412545004946722E-8f,
-3.50196060308781257119E-7f,
 2.40648494783721712015E-6f,
-1.93619797416608296024E-5f,
 1.95215518471351631108E-4f,
-2.85781685962277938680E-3f,
 1.03923736576817238437E-1f,
 2.72062619048444266945E0f,
 0.0f,
 0.0f,
 0.0f,
 0.0f,
 0.0f,
 0.0f
};
float x, y;
x = xx;
if( x <= 0.0f )
	{
	//mtherr( "k1ef", DOMAIN );
	return( MAXNUMF );
	}

if( x <= 2.0f )
	{
	y = x * x - 2.0f;
	y =  ceph_logf( 0.5f * x ) * i1f(x)  +  chbevlf( y, A, 7 ) / x;
	return( y * ceph_expf(x) );
	}

return(  chbevlf( 8.0f/x - 2.0f, B, 10 ) / ceph_sqrtf(x) );

}

/*
   	Modified Bessel function, third kind, integer order
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, knf();
 * int n;
 *
 * y = knf( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of the third kind
 * of order n of the argument.
 *
 * The range is partitioned into the two intervals [0,9.55] and
 * (9.55, infinity).  An ascending power series is used in the
 * low range, and an asymptotic expansion in the high range.
 *
*/


__ATTR_ALWAYS_INLINE__


static inline
float knf( const int nnn, const float xx ) {
float x, k, kf, nk1f, nkf, zn, t, s, z0, z;
float ans, fn, pn, pk, zmn, tlg, tox;
int i, n, nn;
#define EUL 5.772156649015328606065e-1
#define MAXFAC 31
nn = nnn;
x = xx;
if( nn < 0 )
	n = -nn;
else
	n = nn;

if( n > MAXFAC )
	{
overf:
	//mtherr( "knf", OVERFLOW );
	return( MAXNUMF );
	}

if( x <= 0.0 )
	{
	if( x < 0.0 )
	   ;	//mtherr( "knf", DOMAIN );
	else
		//mtherr( "knf", SING );
	return( MAXNUMF );
	}


if( x > 9.55 )
	goto asymp;

ans = 0.0;
z0 = 0.25 * x * x;
fn = 1.0;
pn = 0.0;
zmn = 1.0;
tox = 2.0/x;

if( n > 0 )
	{
	/* compute factorial of n and psi(n) */
	pn = -EUL;
	k = 1.0;
	for( i=1; i<n; i++ )
		{
		pn += 1.0/k;
		k += 1.0;
		fn *= k;
		}

	zmn = tox;

	if( n == 1 )
		{
		ans = 1.0/x;
		}
	else
		{
		nk1f = fn/n;
		kf = 1.0;
		s = nk1f;
		z = -z0;
		zn = 1.0;
		for( i=1; i<n; i++ )
			{
			nk1f = nk1f/(n-i);
			kf = kf * i;
			zn *= z;
			t = nk1f * zn / kf;
			s += t;   
			if( (MAXNUMF - ceph_fabsf(t)) < ceph_fabsf(s) )
				goto overf;
			if( (tox > 1.0) && ((MAXNUMF/tox) < zmn) )
				goto overf;
			zmn *= tox;
			}
		s *= 0.5;
		t = ceph_fabsf(s);
		if( (zmn > 1.0) && ((MAXNUMF/zmn) < t) )
			goto overf;
		if( (t > 1.0) && ((MAXNUMF/t) < zmn) )
			goto overf;
		ans = s * zmn;
		}
	}


tlg = 2.0 * ceph_logf( 0.5 * x );
pk = -EUL;
if( n == 0 )
	{
	pn = pk;
	t = 1.0;
	}
else
	{
	pn = pn + 1.0/n;
	t = 1.0/fn;
	}
s = (pk+pn-tlg)*t;
k = 1.0;
do
	{
	t *= z0 / (k * (k+n));
	pk += 1.0/k;
	pn += 1.0/(k+n);
	s += (pk+pn-tlg)*t;
	k += 1.0;
	}
while( ceph_fabsf(t/s) > MACHEPF );

s = 0.5 * s / zmn;
if( n & 1 )
	s = -s;
ans += s;

return(ans);



/* Asymptotic expansion for Kn(x) */
/* Converges to 1.4e-17 for x > 18.4 */

asymp:

if( x > MAXLOGF )
	{
	//mtherr( "knf", UNDERFLOW );
	return(0.0);
	}
k = n;
pn = 4.0 * k * k;
pk = 1.0;
z0 = 8.0 * x;
fn = 1.0;
t = 1.0;
s = t;
nkf = MAXNUMF;
i = 0;
do
	{
	z = pn - pk * pk;
	t = t * z /(fn * z0);
	nk1f = ceph_fabsf(t);
	if( (i >= n) && (nk1f > nkf) )
		{
		goto adone;
		}
	nkf = nk1f;
	s += t;
	fn += 1.0;
	pk += 2.0;
	i += 1;
	}
while( ceph_fabsf(t/s) > MACHEPF );

adone:
ans = expf(-x) * ceph_sqrtf( PIF/(2.0*x) ) * s;
return(ans);
}

/*
    	Base 2 logarithm
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, log2f();
 *
 * y = log2f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the base 2 logarithm of x.
 *
 * The argument is separated into its exponent and fractional
 * parts.  If the exponent is between -1 and +1, the base e
 * logarithm of the fraction is approximated by
 *
 *     log(1+x) = x - 0.5 x**2 + x**3 P(x)/Q(x).
 *
 * Otherwise, setting  z = 2(x-1)/x+1),
 * 
 *     log(x) = z + z**3 P(z)/Q(z).
*/


__ATTR_ALWAYS_INLINE__


static inline
float ceph_log2f(const xx) {
const float P[] = {
 7.0376836292E-2,
-1.1514610310E-1,
 1.1676998740E-1,
-1.2420140846E-1,
 1.4249322787E-1,
-1.6668057665E-1,
 2.0000714765E-1,
-2.4999993993E-1,
 3.3333331174E-1
};

#define LOG2EA 0.44269504088896340735992
#define SQRTH 0.70710678118654752440
float x, y, z;
int e;
x = xx;
/* Test for domain */
if( x <= 0.0 )
	{
	if( x == 0.0 )
	    ;	//mtherr( fname, SING );
	else
		//mtherr( fname, DOMAIN );
	return( MINLOGF/LOGE2F );
	}

/* separate mantissa from exponent */
x = ceph_frexpf( x, &e );


/* logarithm using log(1+x) = x - .5x**2 + x**3 P(x)/Q(x) */

if( x < SQRTH )
	{
	e -= 1;
	x = 2.0*x - 1.0;
	}	
else
	{
	x = x - 1.0;
	}

z = x*x;
y = x * ( z * polevlf( x, P, 8 ) );
y = y - 0.5 * z;   /*  y - 0.5 * x**2  */


/* Multiply log of fraction by log2(e)
 * and base 2 exponent by 1
 *
 * ***CAUTION***
 *
 * This sequence of operations is critical and it may
 * be horribly defeated by some compiler optimizers.
 */
z = y * LOG2EA;
z += x * LOG2EA;
z += y;
z += x;
z += (float )e;
return( z );
}

/*
   	Common logarithm
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, log10f();
 *
 * y = log10f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns logarithm to the base 10 of x.
 *
 * The argument is separated into its exponent and fractional
 * parts.  The logarithm of the fraction is approximated by
 *
 *     log(1+x) = x - 0.5 x**2 + x**3 P(x).
*/


__ATTR_ALWAYS_INLINE__


static inline
float ceph_log10f(const float xx) {
const float P[] = {
 7.0376836292E-2,
-1.1514610310E-1,
 1.1676998740E-1,
-1.2420140846E-1,
 1.4249322787E-1,
-1.6668057665E-1,
 2.0000714765E-1,
-2.4999993993E-1,
 3.3333331174E-1
};

#define SQRTH 0.70710678118654752440
#define L102A 3.0078125E-1
#define L102B 2.48745663981195213739E-4
#define L10EA 4.3359375E-1
#define L10EB 7.00731903251827651129E-4

constexpr float MAXL10 = 38.230809449325611792;
float x, y, z;
int e;

x = xx;
/* Test for domain */
if( x <= 0.0 )
	{
	if( x == 0.0 )
	     ;	//mtherr( fname, SING );
	else
		//mtherr( fname, DOMAIN );
	return( -MAXL10 );
	}

/* separate mantissa from exponent */

x = ceph_frexpf( x, &e );

/* logarithm using log(1+x) = x - .5x**2 + x**3 P(x) */

if( x < SQRTH )
	{
	e -= 1;
	x = 2.0*x - 1.0;
	}	
else
	{
	x = x - 1.0;
	}


/* rational form */
z = x*x;
y = x * ( z * polevlf( x, P, 8 ) );
y = y - 0.5 * z;   /*  y - 0.5 * x**2  */

/* multiply log of fraction by log10(e)
 * and base 2 exponent by log10(2)
 */
z = (x + y) * L10EB;  /* accumulate terms in order of size */
z += y * L10EA;
z += x * L10EA;
x = e;
z += x * L102B;
z += x * L102A;
return( z );
}

/*SYNOPSIS:
 *
 * float x, y, logf();
 *
 * y = logf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the base e (2.718...) logarithm of x.
 *
 * The argument is separated into its exponent and fractional
 * parts.  If the exponent is between -1 and +1, the logarithm
 * of the fraction is approximated by
 *
 *     log(1+x) = x - 0.5 x**2 + x**3 P(x)
*/


__ATTR_ALWAYS_INLINE__


static inline
float ceph_logf( const float xx ) {
register float y;
float x, z, fe;
int e;

x = xx;
fe = 0.0;
/* Test for domain */
if( x <= 0.0 )
	{
	if( x == 0.0 )
		;//mtherr( "logf", SING );
	else
		//mtherr( "logf", DOMAIN );
	return( MINLOGF );
	}

x = frexpf( x, &e );
if( x < SQRTHF )
	{
	e -= 1;
	x = x + x - 1.0; /*  2x - 1  */
	}	
else
	{
	x = x - 1.0;
	}
z = x * x;
/* 3.4e-9 */
/*
p = logfcof;
y = *p++ * x;
for( i=0; i<8; i++ )
	{
	y += *p++;
	y *= x;
	}
y *= z;
*/

y =
(((((((( 7.0376836292E-2 * x
- 1.1514610310E-1) * x
+ 1.1676998740E-1) * x
- 1.2420140846E-1) * x
+ 1.4249322787E-1) * x
- 1.6668057665E-1) * x
+ 2.0000714765E-1) * x
- 2.4999993993E-1) * x
+ 3.3333331174E-1) * x * z;

if( e )
	{
	fe = e;
	y += -2.12194440e-4 * fe;
	}

y +=  -0.5 * z;  /* y - 0.5 x^2 */
z = x + y;   /* ... + x  */

if( e )
	z += 0.693359375 * fe;

return( z );
}

/*
   Negative binomial distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k, n;
 * float p, y, nbdtrf();
 *
 * y = nbdtrf( k, n, p );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the sum of the terms 0 through k of the negative
 * binomial distribution:
 *
 *   k
 *   --  ( n+j-1 )   n      j
 *   >   (       )  p  (1-p)
 *   --  (   j   )
 *  j=0
 *
*/


__ATTR_ALWAYS_INLINE__


static inline
float nbdtrcf( const int k, const int n, const float pp ) {
float dk, dn, p;
p = pp;
if( (p < 0.0) || (p > 1.0) )
	goto domerr;
if( k < 0 )
	{
domerr:
	//mtherr( "nbdtrf", DOMAIN );
	return( 0.0 );
	}

dk = k+1;
dn = n;
return( incbetf( dk, dn, 1.0 - p ) );
}



__ATTR_ALWAYS_INLINE__


static inline
float nbdtrf( const int k, const int n, const float pp ) {
float dk, dn, p;
p = pp;
if( (p < 0.0) || (p > 1.0) )
	goto domerr;
if( k < 0 )
	{
domerr:
	//mtherr( "nbdtrf", DOMAIN );
	return( 0.0 );
	}
dk = k+1;
dn = n;
return( incbetf( dn, dk, p ) );
}

/*
   	Power function
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, z, powf();
 *
 * z = powf( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes x raised to the yth power.  Analytically,
 *
 *      x**y  =  exp( y log(x) ).
 *
 * Following Cody and Waite, this program uses a lookup table
 * of 2**-i/16 and pseudo extended precision arithmetic to
 * obtain an extra three bits of accuracy in both the logarithm
 * and the exponential.
*/

__ATTR_ALWAYS_INLINE__


static inline
float powf( const float x, const float y ) {
/* 2^(-i/16)
 * The decimal values are rounded to 24-bit precision
 */
 const float A[] = {
  1.00000000000000000000E0,
 9.57603275775909423828125E-1,
 9.17004048824310302734375E-1,
 8.78126084804534912109375E-1,
 8.40896427631378173828125E-1,
 8.05245161056518554687500E-1,
 7.71105408668518066406250E-1,
 7.38413095474243164062500E-1,
 7.07106769084930419921875E-1,
 6.77127778530120849609375E-1,
 6.48419797420501708984375E-1,
 6.20928883552551269531250E-1,
 5.94603538513183593750000E-1,
 5.69394290447235107421875E-1,
 5.45253872871398925781250E-1,
 5.22136867046356201171875E-1,
  5.00000000000000000000E-1
};
/* continuation, for even i only
 * 2^(i/16)  =  A[i] + B[i/2]
 */
const float B[] = {
 0.00000000000000000000E0,
-5.61963907099083340520586E-9,
-1.23776636307969995237668E-8,
 4.03545234539989593104537E-9,
 1.21016171044789693621048E-8,
-2.00949968760174979411038E-8,
 1.89881769396087499852802E-8,
-6.53877009617774467211965E-9,
 0.00000000000000000000E0
};

/* 1 / A[i]
 * The decimal values are full precision
 */
const float Ainv[] = {
 1.00000000000000000000000E0,
 1.04427378242741384032197E0,
 1.09050773266525765920701E0,
 1.13878863475669165370383E0,
 1.18920711500272106671750E0,
 1.24185781207348404859368E0,
 1.29683955465100966593375E0,
 1.35425554693689272829801E0,
 1.41421356237309504880169E0,
 1.47682614593949931138691E0,
 1.54221082540794082361229E0,
 1.61049033194925430817952E0,
 1.68179283050742908606225E0,
 1.75625216037329948311216E0,
 1.83400808640934246348708E0,
 1.91520656139714729387261E0,
 2.00000000000000000000000E0
};
#define F W
#define Fa Wa
#define Fb Wb
#define G W
#define Ga Wa
#define Gb u
#define H W
#define Ha Wb
#define Hb Wb
float u, w, z, W, Wa, Wb, ya, yb;
/* float F, Fa, Fb, G, Ga, Gb, H, Ha, Hb */
int e, i, nflg;
/* Find a multiple of 1/16 that is within 1/16 of x. */
#define reduc(x)  0.0625 * ceph_floorf( 16 * (x) )
#define MEXP 2048.0
#define MNEXP -2400.0
/* log2(e) - 1 */
#define LOG2EA 0.44269504088896340736F
nflg = 0;	/* flag = 1 if x<0 raised to integer power */
w = ceph_floorf(y);
if( w < 0 )
	z = -w;
else
	z = w;
if( (w == y) && (z < 32768.0) )
	{
	i = w;
	w = ceph_powif( x, i );
	return( w );
	}


if( x <= 0.0F )
	{
	if( x == 0.0 )
		{
		if( y == 0.0 )
			return( 1.0 );  /*   0**0   */
		else  
			return( 0.0 );  /*   0**y   */
		}
	else
		{
		if( w != y )
			{ /* noninteger power of negative number */
			//mtherr( fname, DOMAIN );
			return(0.0);
			}
		nflg = 1;
		if( x < 0 )
			x = -x;
		}
	}

/* separate significand from exponent */
x = ceph_frexpf( x, &e );

/* find significand in antilog table A[] */
i = 1;
if( x <= A[9] )
	i = 9;
if( x <= A[i+4] )
	i += 4;
if( x <= A[i+2] )
	i += 2;
if( x >= A[1] )
	i = -1;
i += 1;


/* Find (x - A[i])/A[i]
 * in order to compute log(x/A[i]):
 *
 * log(x) = log( a x/a ) = log(a) + log(x/a)
 *
 * log(x/a) = log(1+v),  v = x/a - 1 = (x-a)/a
 */
x -= A[i];
x -= B[ i >> 1 ];
x *= Ainv[i];


/* rational approximation for log(1+v):
 *
 * log(1+v)  =  v  -  0.5 v^2  +  v^3 P(v)
 * Theoretical relative error of the approximation is 3.5e-11
 * on the interval 2^(1/16) - 1  > v > 2^(-1/16) - 1
 */
z = x*x;
w = (((-0.1663883081054895  * x
      + 0.2003770364206271) * x
      - 0.2500006373383951) * x
      + 0.3333331095506474) * x * z;
w -= 0.5 * z;

/* Convert to base 2 logarithm:
 * multiply by log2(e)
 */
w = w + LOG2EA * w;
/* Note x was not yet added in
 * to above rational approximation,
 * so do it now, while multiplying
 * by log2(e).
 */
z = w + LOG2EA * x;
z = z + x;

/* Compute exponent term of the base 2 logarithm. */
w = -i;
w *= 0.0625;  /* divide by 16 */
w += e;
/* Now base 2 log of x is w + z. */

/* Multiply base 2 log by y, in extended precision. */

/* separate y into large part ya
 * and small part yb less than 1/16
 */
ya = reduc(y);
yb = y - ya;


F = z * y  +  w * yb;
Fa = reduc(F);
Fb = F - Fa;

G = Fa + w * ya;
Ga = reduc(G);
Gb = G - Ga;

H = Fb + Gb;
Ha = reduc(H);
w = 16 * (Ga + Ha);

/* Test the power of 2 for overflow */
if( w > MEXP )
	{
	//mtherr( fname, OVERFLOW );
	return( MAXNUMF );
	}

if( w < MNEXP )
	{
	//mtherr( fname, UNDERFLOW );
	return( 0.0 );
	}

e = w;
Hb = H - Ha;

if( Hb > 0.0 )
	{
	e += 1;
	Hb -= 0.0625;
	}

/* Now the product y * log2(x)  =  Hb + e/16.0.
 *
 * Compute base 2 exponential of Hb,
 * where -0.0625 <= Hb <= 0.
 * Theoretical relative error of the approximation is 2.8e-12.
 */
/*  z  =  2**Hb - 1    */
z = ((( 9.416993633606397E-003 * Hb
      + 5.549356188719141E-002) * Hb
      + 2.402262883964191E-001) * Hb
      + 6.931471791490764E-001) * Hb;

/* Express e/16 as an integer plus a negative number of 16ths.
 * Find lookup table entry for the fractional power of 2.
 */
if( e < 0 )
	i = -( -e >> 4 );
else
	i = (e >> 4) + 1;
e = (i << 4) - e;
w = A[e];
z = w + w * z;      /*    2**-e * ( 1 + (2**Hb-1) )    */
z = ceph_ldexpf( z, i );  /* multiply by integer power of 2 */

if( nflg )
	{
/* For negative x,
 * find out if the integer exponent
 * is odd or even.
 */
	w = 2 * ceph_floorf( (float) 0.5 * w );
	if( w != y )
		z = -z; /* odd exponent */
	}

return( z );
}

/*
   	Real raised to integer power
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, powif();
 * int n;
 *
 * y = powif( x, n );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns argument x raised to the nth power.
 * The routine efficiently decomposes n as a sum of powers of
 * two. The desired power is a product of two-to-the-kth
 * powers of x.  Thus to compute the 32767 power of x requires
 * 28 multiplications instead of 32767 multiplications.
*/

__ATTR_ALWAYS_INLINE__


static inline
float powif( const float x, const int nn ) {
int n, e, sign, asign, lx;
float w, y, s;

if( x == 0.0 )
	{
	if( nn == 0 )
		return( 1.0 );
	else if( nn < 0 )
		return( MAXNUMF );
	else
		return( 0.0 );
	}

if( nn == 0 )
	return( 1.0 );


if( x < 0.0 )
	{
	asign = -1;
	x = -x;
	}
else
	asign = 0;


if( nn < 0 )
	{
	sign = -1;
	n = -nn;
/*
	x = 1.0/x;
*/
	}
else
	{
	sign = 0;
	n = nn;
	}

/* Overflow detection */

/* Calculate approximate logarithm of answer */
s = ceph_frexpf( x, &lx );
e = (lx - 1)*n;
if( (e == 0) || (e > 64) || (e < -64) )
	{
	s = (s - 7.0710678118654752e-1) / (s +  7.0710678118654752e-1);
	s = (2.9142135623730950 * s - 0.5 + lx) * nn * LOGE2F;
	}
else
	{
	s = LOGE2F * e;
	}

if( s > MAXLOGF )
	{
	//mtherr( "powi", OVERFLOW );
	y = MAXNUMF;
	goto done;
	}

if( s < MINLOGF )
	return(0.0);

/* Handle tiny denormal answer, but with less accuracy
 * since roundoff error in 1.0/x will be amplified.
 * The precise demarcation should be the gradual underflow threshold.
 */
if( s < (-MAXLOGF+2.0) )
	{
	x = 1.0/x;
	sign = 0;
	}

/* First bit of the power */
if( n & 1 )
	y = x;
		
else
	{
	y = 1.0;
	asign = 0;
	}

w = x;
n >>= 1;
while( n )
	{
	w = w * w;	/* arg to the 2-to-the-kth power */
	if( n & 1 )	/* if that bit is set, then include in product */
		y *= w;
	n >>= 1;
	}


done:

if( asign )
	y = -y; /* odd power of negative number */
if( sign )
	y = 1.0/y;
return(y);
}

/*
   	Hyperbolic sine and cosine integrals
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, Chi, Shi;
 *
 * shichi( x, &Chi, &Shi );
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integrals
 *
 *                            x
 *                            -
 *                           | |   cosh t - 1
 *   Chi(x) = eul + ln x +   |    -----------  dt,
 *                         | |          t
 *                          -
 *                          0
 *
 *               x
 *               -
 *              | |  sinh t
 *   Shi(x) =   |    ------  dt
 *            | |       t
 *             -
 *             0
 *
 * where eul = 0.57721566490153286061 is Euler's constant.
 * The integrals are evaluated by power series for x < 8
 * and by Chebyshev expansions for x between 8 and 88.
 * For large x, both functions approach exp(x)/2x.
 * Arguments greater than 88 in magnitude return MAXNUM.
 *
*/

__ATTR_ALWAYS_INLINE__


static inline
int shichif( const float xx, float *si, float *ci ) {
/* x exp(-x) shi(x), inverted interval 8 to 18 */
const float S1[] = {
-3.56699611114982536845E-8,
 1.44818877384267342057E-7,
 7.82018215184051295296E-7,
-5.39919118403805073710E-6,
-3.12458202168959833422E-5,
 8.90136741950727517826E-5,
 2.02558474743846862168E-3,
 2.96064440855633256972E-2,
 1.11847751047257036625E0
};

/* x exp(-x) shi(x), inverted interval 18 to 88 */
const float S2[] = {
 1.69050228879421288846E-8,
 1.25391771228487041649E-7,
 1.16229947068677338732E-6,
 1.61038260117376323993E-5,
 3.49810375601053973070E-4,
 1.28478065259647610779E-2,
 1.03665722588798326712E0
};


/* x exp(-x) chin(x), inverted interval 8 to 18 */
const float C1[] = {
 1.31458150989474594064E-8,
-4.75513930924765465590E-8,
-2.21775018801848880741E-7,
 1.94635531373272490962E-6,
 4.33505889257316408893E-6,
-6.13387001076494349496E-5,
-3.13085477492997465138E-4,
 4.97164789823116062801E-4,
 2.64347496031374526641E-2,
 1.11446150876699213025E0
};

/* x exp(-x) chin(x), inverted interval 18 to 88 */
const float C2[] = {
-3.00095178028681682282E-9,
 7.79387474390914922337E-8,
 1.06942765566401507066E-6,
 1.59503164802313196374E-5,
 3.49592575153777996871E-4,
 1.28475387530065247392E-2,
 1.03665693917934275131E0
};
/* Sine and cosine integrals */
#define EUL 0.57721566490153286061
float x, k, z, c, s, a;
short sign;

x = xx;
if( x < 0.0 )
	{
	sign = -1;
	x = -x;
	}
else
	sign = 0;


if( x == 0.0 )
	{
	*si = 0.0;
	*ci = -MAXNUMF;
	return( 0 );
	}

if( x >= 8.0 )
	goto chb;

z = x * x;

/*	Direct power series expansion	*/

a = 1.0;
s = 1.0;
c = 0.0;
k = 2.0;

do
	{
	a *= z/k;
	c += a/k;
	k += 1.0;
	a /= k;
	s += a/k;
	k += 1.0;
	}
while( ceph_fabsf(a/s) > MACHEPF );

s *= x;
goto done;


chb:

if( x < 18.0 )
	{
	a = (576.0/x - 52.0)/10.0;
	k = ceph_expf(x) / x;
	s = k * chbevlf( a, S1, 9 );
	c = k * chbevlf( a, C1, 10 );
	goto done;
	}

if( x <= 88.0 )
	{
	a = (6336.0/x - 212.0)/70.0;
	k = ceph_expf(x) / x;
	s = k * chbevlf( a, S2, 7 );
	c = k * chbevlf( a, C2, 7 );
	goto done;
	}
else
	{
	if( sign )
		*si = -MAXNUMF;
	else
		*si = MAXNUMF;
	*ci = MAXNUMF;
	return(0);
	}
done:
if( sign )
	s = -s;

*si = s;

*ci = EUL + ceph_logf(x) + c;
return(0);
}

/*
   	Sine and cosine integrals
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, Ci, Si;
 *
 * sicif( x, &Si, &Ci );
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the integrals
 *
 *                          x
 *                          -
 *                         |  cos t - 1
 *   Ci(x) = eul + ln x +  |  --------- dt,
 *                         |      t
 *                        -
 *                         0
 *             x
 *             -
 *            |  sin t
 *   Si(x) =  |  ----- dt
 *            |    t
 *           -
 *            0
 *
 * where eul = 0.57721566490153286061 is Euler's constant.
 * The integrals are approximated by rational functions.
 * For x > 8 auxiliary functions f(x) and g(x) are employed
 * such that
 *
 * Ci(x) = f(x) sin(x) - g(x) cos(x)
 * Si(x) = pi/2 - f(x) cos(x) - g(x) sin(x)
*/

__ATTR_ALWAYS_INLINE__


static inline
int sicif( const float xx, float *si, float *ci ) {
const float SN[] = {
-8.39167827910303881427E-11,
 4.62591714427012837309E-8,
-9.75759303843632795789E-6,
 9.76945438170435310816E-4,
-4.13470316229406538752E-2,
 1.00000000000000000302E0,
};
const  float SD[] = {
  2.03269266195951942049E-12,
  1.27997891179943299903E-9,
  4.41827842801218905784E-7,
  9.96412122043875552487E-5,
  1.42085239326149893930E-2,
  9.99999999999999996984E-1,
};

const float CN[] = {
 2.02524002389102268789E-11,
-1.35249504915790756375E-8,
 3.59325051419993077021E-6,
-4.74007206873407909465E-4,
 2.89159652607555242092E-2,
-1.00000000000000000080E0,
};
const float CD[] = {
  4.07746040061880559506E-12,
  3.06780997581887812692E-9,
  1.23210355685883423679E-6,
  3.17442024775032769882E-4,
  5.10028056236446052392E-2,
  4.00000000000000000080E0,
};


const float FN4[] = {
  4.23612862892216586994E0,
  5.45937717161812843388E0,
  1.62083287701538329132E0,
  1.67006611831323023771E-1,
  6.81020132472518137426E-3,
  1.08936580650328664411E-4,
  5.48900223421373614008E-7,
};
const float FD4[] = {
/*  1.00000000000000000000E0,*/
  8.16496634205391016773E0,
  7.30828822505564552187E0,
  1.86792257950184183883E0,
  1.78792052963149907262E-1,
  7.01710668322789753610E-3,
  1.10034357153915731354E-4,
  5.48900252756255700982E-7,
};

const float FN8[] = {
  4.55880873470465315206E-1,
  7.13715274100146711374E-1,
  1.60300158222319456320E-1,
  1.16064229408124407915E-2,
  3.49556442447859055605E-4,
  4.86215430826454749482E-6,
  3.20092790091004902806E-8,
  9.41779576128512936592E-11,
  9.70507110881952024631E-14,
};
const float FD8[] = {
/*  1.00000000000000000000E0,*/
  9.17463611873684053703E-1,
  1.78685545332074536321E-1,
  1.22253594771971293032E-2,
  3.58696481881851580297E-4,
  4.92435064317881464393E-6,
  3.21956939101046018377E-8,
  9.43720590350276732376E-11,
  9.70507110881952025725E-14,
};

const float GN4[] = {
  8.71001698973114191777E-2,
  6.11379109952219284151E-1,
  3.97180296392337498885E-1,
  7.48527737628469092119E-2,
  5.38868681462177273157E-3,
  1.61999794598934024525E-4,
  1.97963874140963632189E-6,
  7.82579040744090311069E-9,
};
const float GD4[] = {
/*  1.00000000000000000000E0,*/
  1.64402202413355338886E0,
  6.66296701268987968381E-1,
  9.88771761277688796203E-2,
  6.22396345441768420760E-3,
  1.73221081474177119497E-4,
  2.02659182086343991969E-6,
  7.82579218933534490868E-9,
};

const float GN8[] = {
  6.97359953443276214934E-1,
  3.30410979305632063225E-1,
  3.84878767649974295920E-2,
  1.71718239052347903558E-3,
  3.48941165502279436777E-5,
  3.47131167084116673800E-7,
  1.70404452782044526189E-9,
  3.85945925430276600453E-12,
  3.14040098946363334640E-15,
};
const float GD8[] = {
/*  1.00000000000000000000E0,*/
  1.68548898811011640017E0,
  4.87852258695304967486E-1,
  4.67913194259625806320E-2,
  1.90284426674399523638E-3,
  3.68475504442561108162E-5,
  3.57043223443740838771E-7,
  1.72693748966316146736E-9,
  3.87830166023954706752E-12,
  3.14040098946363335242E-15,
};

#define EUL 0.57721566490153286061
float x, z, c, s, f, g;
int sign;

x = xx;
if( x < 0.0 )
	{
	sign = -1;
	x = -x;
	}
else
	sign = 0;


if( x == 0.0 )
	{
	*si = 0.0;
	*ci = -MAXNUMF;
	return( 0 );
	}


if( x > 1.0e9 )
	{
	*si = PIO2F - ceph_cosf(x)/x;
	*ci = ceph_sinf(x)/x;
	return( 0 );
	}



if( x > 4.0 )
	goto asympt;

z = x * x;
s = x * polevlf( z, SN, 5 ) / polevlf( z, SD, 5 );
c = z * polevlf( z, CN, 5 ) / polevlf( z, CD, 5 );

if( sign )
	s = -s;
*si = s;
*ci = EUL + ceph_logf(x) + c;	/* real part if x < 0 */
return(0);



/* The auxiliary functions are:
 *
 *
 * *si = *si - PIO2;
 * c = cos(x);
 * s = sin(x);
 *
 * t = *ci * s - *si * c;
 * a = *ci * c + *si * s;
 *
 * *si = t;
 * *ci = -a;
 */


asympt:

s = ceph_sinf(x);
c = ceph_cosf(x);
z = 1.0/(x*x);
if( x < 8.0 )
	{
	f = polevlf( z, FN4, 6 ) / (x * p1evlf( z, FD4, 7 ));
	g = z * polevlf( z, GN4, 7 ) / p1evlf( z, GD4, 7 );
	}
else
	{
	f = polevlf( z, FN8, 8 ) / (x * p1evlf( z, FD8, 8 ));
	g = z * polevlf( z, GN8, 8 ) / p1evlf( z, GD8, 9 );
	}
*si = PIO2F - f * c - g * s;
if( sign )
	*si = -( *si );
*ci = f * s - g * c;

return(0);
}

/*
   	Circular sine of angle in degrees
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, sindgf();
 *
 * y = sindgf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Range reduction is into intervals of 45 degrees.
 *
 * Two polynomial approximating functions are employed.
 * Between 0 and pi/4 the sine is approximated by
 *      x  +  x**3 P(x**2).
 * Between pi/4 and pi/2 the cosine is represented as
 *      1  -  x**2 Q(x**2).
*/


__ATTR_ALWAYS_INLINE__


static inline
float ceph_sindgf( const float xx ) {
float x, y, z;
constexpr float T24M1 = 16777215.;
constexpr float PI180 = 0.0174532925199432957692; /* pi/180 */
long j;
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
	 //mtherr( "sindgf", TLOSS );
	return(0.0);
	}
j = 0.022222222222222222222 * x; /* integer part of x/45 */
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

x = x - y * 45.0;
x *= PI180;	/* multiply by pi/180 to convert to radians */

z = x * x;
if( (j==1) || (j==2) )
	{
/*
	y = ((( 2.4462803166E-5 * z
	  - 1.3887580023E-3) * z
	  + 4.1666650433E-2) * z
	  - 4.9999999968E-1) * z
	  + 1.0;
*/

/* measured relative error in +/- pi/4 is 7.8e-8 */
	y = ((  2.443315711809948E-005 * z
	  - 1.388731625493765E-003) * z
	  + 4.166664568298827E-002) * z * z;
	y -= 0.5 * z;
	y += 1.0;
	}
else
	{
/* Theoretical relative error = 3.8e-9 in [-pi/4, +pi/4] */
	y = ((-1.9515295891E-4 * z
	     + 8.3321608736E-3) * z
	     - 1.6666654611E-1) * z * x;
	y += x;
	}

if(sign < 0)
	y = -y;
return( y);
}


/* Single precision circular cosine
 * test interval: [-pi/4, +pi/4]
 * trials: 10000
 * peak relative error: 8.3e-8
 * rms relative error: 2.2e-8
 */


__ATTR_ALWAYS_INLINE__


static inline
float ceph_cosdgf( const float xx ) {
register float x, y, z;
int j, sign;
/* These are for a 24-bit significand: */
constexpr float T24M1 = 16777215.;
constexpr float PI180 = 0.0174532925199432957692; /* pi/180 */
/* make argument positive */
sign = 1;
x = xx;
if( x < 0 )
	x = -x;

if( x > T24M1 )
	{
	 //mtherr( "cosdgf", TLOSS );
	return(0.0);
	}

j = 0.02222222222222222222222 * x; /* integer part of x/PIO4 */
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

x = x - y * 45.0; /* x mod 45 degrees */
x *= PI180;	/* multiply by pi/180 to convert to radians */

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

/*
    	Circular sine
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, sinf();
 *
 * y = sinf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Range reduction is into intervals of pi/4.  The reduction
 * error is nearly eliminated by contriving an extended precision
 * modular arithmetic.
 *
 * Two polynomial approximating functions are employed.
 * Between 0 and pi/4 the sine is approximated by
 *      x  +  x**3 P(x**2).
 * Between pi/4 and pi/2 the cosine is represented as
 *      1  -  x**2 Q(x**2).
 *
*/


__ATTR_ALWAYS_INLINE__


static inline
float ceph_sinf( const float xx ) {
constexpr float FOPI = 1.27323954473516;
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

const float sincof[] = {
-1.9515295891E-4,
 8.3321608736E-3,
-1.6666654611E-1
};
const float coscof[] = {
 2.443315711809948E-005,
-1.388731625493765E-003,
 4.166664568298827E-002
};
float *p;
float x, y, z;
register unsigned long j;
register int sign;
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
	return(0.0);
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


/* Single precision circular cosine
 * test interval: [-pi/4, +pi/4]
 * trials: 10000
 * peak relative error: 8.3e-8
 * rms relative error: 2.2e-8
 */


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
constexpr float FOPI = 1.27323954473516;
constexpr float DP1 = 0.78515625;
constexpr float DP2 = 2.4187564849853515625e-4;
constexpr float DP3 = 3.77489497744594108e-8;
constexpr float lossth = 8192.;
constexpr float T24M1 = 16777215.;
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
	return(0.0);
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

/*
    	Hyperbolic sine
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, sinhf();
 *
 * y = sinhf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns hyperbolic sine of argument in the range MINLOGF to
 * MAXLOGF.
 *
 * The range is partitioned into two segments.  If |x| <= 1, a
 * polynomial approximation is used.
 * Otherwise the calculation is sinh(x) = ( exp(x) - exp(-x) )/2.
 *
*/

__ATTR_ALWAYS_INLINE__


static inline
float ceph_sinhf( const float xx ) {
register float z;
float x;

x = xx;
if( xx < 0 )
	z = -x;
else
	z = x;

if( z > MAXLOGF )
	{
	//mtherr( "sinhf", DOMAIN );
	;
	if( x > 0 )
		return( MAXNUMF );
	else
		return( -MAXNUMF );
	}
if( z > 1.0 )
	{
	z = ceph_expf(z);
	z = 0.5*z - (0.5/z);
	if( x < 0 )
		z = -z;
	}
else
	{
	z = x * x;
	z =
	(( 2.03721912945E-4 * z
	  + 8.33028376239E-3) * z
	  + 1.66667160211E-1) * z * x
	  + x;
	}
return( z );
}


/*
   	Square root
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, sqrtf();
 *
 * y = sqrtf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the square root of x.
 *
 * Range reduction involves isolating the power of two of the
 * argument and using a polynomial approximation to obtain
 * a rough value for the square root.  Then Heron's iteration
 * is used three times to converge to an accurate value.
 *
*/


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

/*
    	Circular tangent of angle in degrees
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, tandgf();
 *
 * y = tandgf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the circular tangent of the radian argument x.
 *
 * Range reduction is into intervals of 45 degrees.
*/


__ATTR_ALWAYS_INLINE__


static inline
float ceph_tancotf( float xx, int cotflg ) {
float x, y, z, zz;
long j;
int sign;
constexpr  float T24M1 = 16777215.;
constexpr  float PI180 = 0.0174532925199432957692; /* pi/180 */

/* make argument positive but save the sign */
if( xx < 0.0 )
	{
	x = -xx;
	sign = -1;
	}
else
	{
	x = xx;
	sign = 1;
	}

if( x > T24M1 )
	{
	if( cotflg )
		//mtherr( "cotdgf", TLOSS );
	else
		//mtherr( "tandgf", TLOSS );
	return(0.0);
	}

/* compute x mod PIO4 */
j = 0.022222222222222222222 * x; /* integer part of x/45 */
y = j;

/* map zeros and singularities to origin */
if( j & 1 )
	{
	j += 1;
	y += 1.0;
	}

z = x - y * 45.0;
z *= PI180;	/* multiply by pi/180 to convert to radians */

zz = z * z;

if( x > 1.0e-4 )
	{
/* 1.7e-8 relative error in [-pi/4, +pi/4] */
	y =
	((((( 9.38540185543E-3 * zz
	+ 3.11992232697E-3) * zz
	+ 2.44301354525E-2) * zz
	+ 5.34112807005E-2) * zz
	+ 1.33387994085E-1) * zz
	+ 3.33331568548E-1) * zz * z
	+ z;
	}
else
	{
	y = z;
	}

if( j & 2 )
	{
	if( cotflg )
		y = -y;
	else
		{
		if( y != 0.0 )
			{
			y = -1.0/y;
			}
		else
			{
			//mtherr( "tandgf", SING );
			y = MAXNUMF;
			}
		}
	}
else
	{
	if( cotflg )
		{
		if( y != 0.0 )
			y = 1.0/y;
		else
			{
			//mtherr( "cotdgf", SING );
			y = MAXNUMF;
			}
		}
	}

if( sign < 0 )
	y = -y;

return( y );
}



__ATTR_ALWAYS_INLINE__


static inline
float ceph_tandgf( const float x ) {

return( ceph_tancotf(x,0) );
}


__ATTR_ALWAYS_INLINE__


static inline
float cotdgf( const float x ) {
if( x == 0.0 )
	{
	//mtherr( "cotdgf", SING );
	return( MAXNUMF );
	}
return( ceph_tancotf(x,1) );
}

/*
    	Circular tangent
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, tanf();
 *
 * y = tanf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the circular tangent of the radian argument x.
 *
 * Range reduction is modulo pi/4.  A polynomial approximation
 * is employed in the basic interval [0, pi/4].
*/

__ATTR_ALWAYS_INLINE__


static inline
float ceph_tancotf( const float xx, const int cotflg ) {
const float DP1 = 0.78515625;
const float DP2 = 2.4187564849853515625e-4;
const float DP3 = 3.77489497744594108e-8;
const float FOPI = 1.27323954473516;  /* 4/pi */
const float lossth = 8192.;
/*static float T24M1 = 16777215.;*/
float x, y, z, zz;
long j;
int sign;


/* make argument positive but save the sign */
if( xx < 0.0 )
	{
	x = -xx;
	sign = -1;
	}
else
	{
	x = xx;
	sign = 1;
	}

if( x > lossth )
	{
	if( cotflg )
		//mtherr( "cotf", TLOSS );
	else
		//mtherr( "tanf", TLOSS );
	return(0.0);
	}

/* compute x mod PIO4 */
j = FOPI * x; /* integer part of x/(PI/4) */
y = j;

/* map zeros and singularities to origin */
if( j & 1 )
	{
	j += 1;
	y += 1.0;
	}

z = ((x - y * DP1) - y * DP2) - y * DP3;

zz = z * z;

if( x > 1.0e-4 )
	{
/* 1.7e-8 relative error in [-pi/4, +pi/4] */
	y =
	((((( 9.38540185543E-3 * zz
	+ 3.11992232697E-3) * zz
	+ 2.44301354525E-2) * zz
	+ 5.34112807005E-2) * zz
	+ 1.33387994085E-1) * zz
	+ 3.33331568548E-1) * zz * z
	+ z;
	}
else
	{
	y = z;
	}

if( j & 2 )
	{
	if( cotflg )
		y = -y;
	else
		y = -1.0/y;
	}
else
	{
	if( cotflg )
		y = 1.0/y;
	}

if( sign < 0 )
	y = -y;

return( y );
}



__ATTR_ALWAYS_INLINE__


static inline
float ceph_tanf( const float x ) {

return( ceph_tancotf(x,0) );
}


__ATTR_ALWAYS_INLINE__


static inline
float ceph_cotf( const float x ) {
if( x == 0.0 )
	{
	//mtherr( "cotf", SING );
	return( MAXNUMF );
	}
return( ceph_tancotf(x,1) );
}

/*
   	Hyperbolic tangent
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, tanhf();
 *
 * y = tanhf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns hyperbolic tangent of argument in the range MINLOG to
 * MAXLOG.
 *
 * A polynomial approximation is used for |x| < 0.625.
 * Otherwise,
 *
 *    tanh(x) = sinh(x)/cosh(x) = 1  -  2/(exp(2x) + 1).
*/


__ATTR_ALWAYS_INLINE__


static inline
float tanhf( const float xx ) {
float x, z;

if( xx < 0 )
	x = -xx;
else
	x = xx;

if( x > 0.5 * MAXLOGF )
	{
	if( xx > 0 )
		return( 1.0 );
	else
		return( -1.0 );
	}
if( x >= 0.625 )
	{
	x = ceph_expf(x+x);
	z =  1.0  - 2.0/(x + 1.0);
	if( xx < 0 )
		z = -z;
	}
else
	{
	z = x * x;
	z =
	(((( -5.70498872745E-3 * z
	  + 2.06390887954E-2) * z
	  - 5.37397155531E-2) * z
	  + 1.33314422036E-1) * z
	  - 3.33332819422E-1) * z * xx
	  + xx;
	}
return( z );
}




/** SYNOPSIS:
 *
 * int N;
 * float x, y, coef[N+1], polevlf[];
 *
 * y = polevlf( x, coef, N );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates polynomial of degree N:
 *
 *                     2          N
 * y  =  C  + C x + C x  +...+ C x
 *        0    1     2          N
 *
 * Coefficients are stored in reverse order:
 *
 * coef[0] = C  , ..., coef[N] = C  .
 *            N                   0
 *
 *  The function p1evl() assumes that coef[N] = 1.0 and is
 * omitted from the array.  Its calling arguments are
 * otherwise the same as polevl().
*/


__ATTR_ALWAYS_INLINE__


static inline
float polevlf( float xx, float * __restrict coef, int N ) {
float ans, x;
float * __restrict p;
int i;
x = xx;
p = coef;
ans = *p++;
/*
for( i=0; i<N; i++ )
	ans = ans * x  +  *p++;
*/

i = N;
do{
	ans = ans * x  +  *p++;
}
while( --i );
return( ans );
}

/*							p1evl()	*/
/*                                          N
 * Evaluate polynomial when coefficient of x  is 1.0.
 * Otherwise same as polevl.
 */


__ATTR_ALWAYS_INLINE__


static inline
float p1evlf( float xx,
              float * __restrict coef,
	      int N ) {
float ans, x;
float *p;
int i;
x = xx;
p = coef;
ans = x + *p++;
i = N-1;
do{
	ans = ans * x  + *p++;
}
while( --i );

return( ans );
}



#endif /*__GMS_CEPHES_H__*/
