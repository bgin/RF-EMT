
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

/*
    Slightly adapted version of CEPHES library.
    The rationale: removed the overhead of GLIBC calls for math.h declared
    mathematical functions (elementary and transcendental).
*/

#ifndef __GMS_CEPHES_DOUBLE_H__
#define __GMS_CEPHES_DOUBLE_H__ 290720260849

#include <limits>
#include "GMS_config.h"

/* Define if the `long double' type works.  */
#define HAVE_LONG_DOUBLE 1

/* Define as the return type of signal handlers (int or void).  */
#define RETSIGTYPE void

/* Define if you have the ANSI C header files.  */
#define STDC_HEADERS 1

/* Define if your processor stores words with the most significant
   byte first (like Motorola and SPARC, unlike Intel and VAX).  */
/* #undef WORDS_BIGENDIAN */

/* Define if floating point words are bigendian.  */
/* #undef FLOAT_WORDS_BIGENDIAN */

/* The number of bytes in a int.  */
#define SIZEOF_INT 4

/* Define if you have the <string.h> header file.  */
#define HAVE_STRING_H 1

/* Name of package */
//#define PACKAGE "cephes"

/* Version number of package */
//#define VERSION "2.7"

/* Constant definitions for math error conditions
 */

#define DOMAIN		1	/* argument domain error */
#define SING		2	/* argument singularity */
#define OVERFLOW	3	/* overflow range error */
#define UNDERFLOW	4	/* underflow range error */
#define TLOSS		5	/* total loss of precision */
#define PLOSS		6	/* partial loss of precision */

#define EDOM		33
#define ERANGE		34

#define UNK 1

/* If you define UNK, then be sure to set BIGENDIAN properly. */
#ifdef FLOAT_WORDS_BIGENDIAN
#define BIGENDIAN 1
#else
#define BIGENDIAN 0
#endif
/* Define this `volatile' if your compiler thinks
 * that floating point arithmetic obeys the associative
 * and distributive laws.  It will defeat some optimizations
 * (but probably not enough of them).
 *
 * #define VOLATILE volatile
 */
#define VOLATILE

/* For 12-byte long doubles on an i386, pad a 16-bit short 0
 * to the end of real constants initialized by integer arrays.
 *
 * #define XPD 0,
 *
 * Otherwise, the type is 10 bytes long and XPD should be
 * defined blank (e.g., Microsoft C).
 *
 * #define XPD
 */
#define XPD 0,

/* Define to support tiny denormal numbers, else undefine. */
#define DENORMAL 1

/* Define to ask for infinity support, else undefine. */
#define INFINITIES 1

/* Define to ask for support of numbers that are Not-a-Number,
   else undefine.  This may automatically define INFINITIES in some files. */
#define NANS 1

/* Define to distinguish between -0.0 and +0.0.  */
#define MINUSZERO 1

/* Define 1 for ANSI C atan2() function
   See atan.c and clog.c. */
#define ANSIC 1

#define EXPMSK 0x800f
#define MEXP 0x7ff
#define NBITS 53


namespace gms 
{

namespace math 
{

namespace
{

constexpr double MACHEP =  1.38777878078144567553E-17;   /* 2**-56 */
constexpr double UFLOWTHRESH =  2.22507385850720138309E-308; /* 2**-1022 */
constexpr double MAXNUM =  1.79769313486231570815E308;    /* 2**1024*(1-MACHEP) */
constexpr double PI     =  3.14159265358979323846;       /* pi */
constexpr double PIO2   =  1.57079632679489661923;       /* pi/2 */
constexpr double PIO4   =  7.85398163397448309616E-1;    /* pi/4 */
constexpr double SQRT2  =  1.41421356237309504880;       /* sqrt(2) */
constexpr double SQRTH  =  7.07106781186547524401E-1;    /* sqrt(2)/2 */
constexpr double LOG2E  =  1.4426950408889634073599;     /* 1/log(2) */
constexpr double SQ2OPI =  7.9788456080286535587989E-1;  /* sqrt( 2/pi ) */
constexpr double LOGE2  =  6.93147180559945309417E-1;    /* log(2) */
constexpr double LOGSQ2 =  3.46573590279972654709E-1;    /* log(2)/2 */
constexpr double THPIO4 =  2.35619449019234492885;       /* 3*pi/4 */
constexpr double TWOOPI =  6.36619772367581343075535E-1; /* 2/pi */
#ifdef DENORMAL
constexpr double MAXLOG =  7.09782712893383996732E2;     /* log(MAXNUM) */
/* double MINLOG = -7.44440071921381262314E2; */     /* log(2**-1074) */
constexpr double MINLOG = -7.451332191019412076235E2;     /* log(2**-1075) */
#else
double MAXLOG =  7.08396418532264106224E2;     /* log 2**1022 */
double MINLOG = -7.08396418532264106224E2;     /* log 2**-1022 */
#endif 
// Not used
#ifdef INFINITIES
double INFINITY = 1.0/0.0;  /* 99e999; */
#else
double INFINITY =  1.79769313486231570815E308;    /* 2**1024*(1-MACHEP) */
#endif
#ifdef NANS
double NAN = 1.0/0.0 - 1.0/0.0;
#else
double NAN = 0.0;
#endif
#ifdef MINUSZERO
double NEGZERO = -0.0;
#else
double NEGZERO = 0.0;
#endif
/*IBMPC*/
#ifdef IBMPC
			/* 2**-53 =  1.11022302462515654042E-16 */
unsigned short MACHEP[4] = {0x0000,0x0000,0x0000,0x3ca0};
unsigned short UFLOWTHRESH[4] = {0x0000,0x0000,0x0000,0x0010};
#ifdef DENORMAL
			/* log(MAXNUM) =  7.09782712893383996732224E2 */
unsigned short MAXLOG[4] = {0x39ef,0xfefa,0x2e42,0x4086};
			/* log(2**-1074) = - -7.44440071921381262314E2 */
/*unsigned short MINLOG[4] = {0x71c3,0x446d,0x4385,0xc087};*/
unsigned short MINLOG[4] = {0x3052,0xd52d,0x4910,0xc087};
#else
			/* log(2**1022) =   7.08396418532264106224E2 */
unsigned short MAXLOG[4] = {0xbcd2,0xdd7a,0x232b,0x4086};
			/* log(2**-1022) = - 7.08396418532264106224E2 */
unsigned short MINLOG[4] = {0xbcd2,0xdd7a,0x232b,0xc086};
#endif
			/* 2**1024*(1-MACHEP) =  1.7976931348623158E308 */
unsigned short MAXNUM[4] = {0xffff,0xffff,0xffff,0x7fef};
unsigned short PI[4]     = {0x2d18,0x5444,0x21fb,0x4009};
unsigned short PIO2[4]   = {0x2d18,0x5444,0x21fb,0x3ff9};
unsigned short PIO4[4]   = {0x2d18,0x5444,0x21fb,0x3fe9};
unsigned short SQRT2[4]  = {0x3bcd,0x667f,0xa09e,0x3ff6};
unsigned short SQRTH[4]  = {0x3bcd,0x667f,0xa09e,0x3fe6};
unsigned short LOG2E[4]  = {0x82fe,0x652b,0x1547,0x3ff7};
unsigned short SQ2OPI[4] = {0x3651,0x33d4,0x8845,0x3fe9};
unsigned short LOGE2[4]  = {0x39ef,0xfefa,0x2e42,0x3fe6};
unsigned short LOGSQ2[4] = {0x39ef,0xfefa,0x2e42,0x3fd6};
unsigned short THPIO4[4] = {0x21d2,0x7f33,0xd97c,0x4002};
unsigned short TWOOPI[4] = {0xc883,0x6dc9,0x5f30,0x3fe4};
#ifdef INFINITIES
unsigned short INFINITY[4] = {0x0000,0x0000,0x0000,0x7ff0};
#else
unsigned short INFINITY[4] = {0xffff,0xffff,0xffff,0x7fef};
#endif
#ifdef NANS
unsigned short NAN[4] = {0x0000,0x0000,0x0000,0x7ffc};
#else
unsigned short NAN[4] = {0x0000,0x0000,0x0000,0x0000};
#endif
#ifdef MINUSZERO
unsigned short NEGZERO[4] = {0x0000,0x0000,0x0000,0x8000};
#else
unsigned short NEGZERO[4] = {0x0000,0x0000,0x0000,0x0000};
#endif
#endif
}


/* Complex numeral.  */
typedef struct
	{
	double r;
	double i;
	} cmplx;

#ifdef HAVE_LONG_DOUBLE
/* Long double complex numeral.  */
typedef struct
	{
	long double r;
	long double i;
	} cmplxl;
#endif

///////////////////////////////////////////////////////////////////
__ATTR_ALWAYS_INLINE__
static inline 
double polevl( x, coef, N )
double x;
double coef[];
int N;
{
double ans;
int i;
double *p;

p = coef;
ans = *p++;
i = N;

do
	ans = ans * x  +  *p++;
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
double p1evl( double x, double coef[], std::int32_t N )
{
double ans;
double *p;
int i;

p = coef;
ans = x + *p++;
i = N-1;

do
	ans = ans * x  + *p++;
while( --i );

return( ans );
}
///////////////////////////////////////////////////////////////////////////

__ATTR_ALWAYS_INLINE__
static inline
int signbit(double x)
{
union
	{
	double d;
	short s[4];
	int i[2];
	} u;

u.d = x;
if( sizeof(int) == 4 )
	{

return( u.i[1] < 0 );
	}
else
	{

	return( u.s[3] < 0 );
	}
}

__ATTR_ALWAYS_INLINE__
static inline
int isnan(double x)
{
#ifdef NANS
union
	{
	double d;
	unsigned short s[4];
	unsigned int i[2];
	} u;

u.d = x;

if( sizeof(int) == 4 )
	{
	if( ((u.i[1] & 0x7ff00000) == 0x7ff00000)
	    && (((u.i[1] & 0x000fffff) != 0) || (u.i[0] != 0)))
		return 1;
	return(0);
	}
else
	{ /* size int not 4 */
	if( (u.s[3] & 0x7ff0) == 0x7ff0)
		{
		if( ((u.s[3] & 0x000f) | u.s[2] | u.s[1] | u.s[0]) != 0 )
			return(1);
		}
	return(0);
	} /* size int not 4 */
#else
/* No NANS.  */
return(0);
#endif
}

__ATTR_ALWAYS_INLINE__
static inline
int isfinite(double x)
{
#ifdef INFINITIES
union
	{
	double d;
	unsigned short s[4];
	unsigned int i[2];
	} u;

u.d = x;

if( sizeof(int) == 4 )
	{
	if( (u.i[1] & 0x7ff00000) != 0x7ff00000)
		return 1;
	return(0);
	}
else
	{
	if( (u.s[3] & 0x7ff0) != 0x7ff0)
		return 1;
	return(0);
	}
#else
/* No INFINITY.  */
return(1);
#endif
}

__ATTR_ALWAYS_INLINE__
static inline
double fabs(double x)
{
union
  {
    double d;
    short i[4];
  } u;

u.d = x;
u.i[3] &= 0x7fff;
return( u.d );
}

double ceil(double x)
{
double y;
#ifdef NANS
if( isnan(x) )
	return( x );
#endif
#ifdef INFINITIES
if(!isfinite(x))
	return(x);
#endif
y = floor(x);
if( y < x )
	y += 1.0;
#ifdef MINUSZERO
if( y == 0.0 && x < 0.0 )
	return( NEGZERO );
#endif
return(y);
}

__ATTR_ALWAYS_INLINE__
static inline
double floor(double x)
{
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

union
	{
	double y;
	unsigned short sh[4];
	} u;
unsigned short *p;
int e;
#ifdef NANS
if( isnan(x) )
	return( x );
#endif
#ifdef INFINITIES
if(!isfinite(x))
	return(x);
#endif
#ifdef MINUSZERO
if(x == 0.0L)
	return(x);
#endif
u.y = x;
/* find the exponent (power of 2) */
p = (unsigned short *)&u.sh[3];
e = (( *p >> 4) & 0x7ff) - 0x3ff;
p -= 3;

if( e < 0 )
	{
	if( u.y < 0.0 )
		return( -1.0 );
	else
		return( 0.0 );
	}

e = (NBITS -1) - e;
/* clean out 16 bits at a time */
while( e >= 16 )
	{
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
double frexp( double x, int * pw2 )
{
union
	{
	double y;
	unsigned short sh[4];
	} u;
int i;
#ifdef DENORMAL
int k;
#endif
short *q;

u.y = x;
q = (short *)&u.sh[3];
i  = ( *q >> 4) & 0x7ff;
if( i != 0 )
	goto ieeedon;

/* Number is denormal or zero */
#ifdef DENORMAL
if( u.y == 0.0 )
	{
	*pw2 = 0;
	return( 0.0 );
	}
/* Handle denormal number. */
do
	{
	u.y *= 2.0;
	i -= 1;
	k  = ( *q >> 4) & 0x7ff;
	}
while( k == 0 );
i = i + k;
#endif /* DENORMAL */

ieeedon:

i -= 0x3fe;
*pw2 = i;
*q &= 0x800f;
*q |= 0x3fe0;
return( u.y );
}

__ATTR_ALWAYS_INLINE__
static inline
double ldexp( double x, std::int32_t pw2 )
{
union
	{
	double y;
	unsigned short sh[4];
	} u;
short *q;
int e;
u.y = x;
q = (short *)&u.sh[3];
while( (e = (*q & 0x7ff0) >> 4) == 0 )
	{
	if( u.y == 0.0 )
		{
		return( 0.0 );
		}
/* Input is denormal. */
	if( pw2 > 0 )
		{
		u.y *= 2.0;
		pw2 -= 1;
		}
	if( pw2 < 0 )
		{
		if( pw2 < -53 )
			return(0.0);
		u.y /= 2.0;
		pw2 += 1;
		}
	if( pw2 == 0 )
		return(u.y);
	}
#endif /* not DEC */

e += pw2;

/* Handle overflow */
if( e >= MEXP )
	return( 2.0*MAXNUM );


/* Handle denormalized results */
if( e < 1 )
	{
#ifdef DENORMAL
	if( e < -53 )
		return(0.0);
	*q &= 0x800f;
	*q |= 0x10;
	/* For denormals, significant bits may be lost even
	   when dividing by 2.  Construct 2^-(1-e) so the result
	   is obtained with only one multiplication.  */
	u.y *= ldexp(1.0, e-1);
	return(u.y);
#else
	return(0.0);
#endif
	}
else
	{

	*q &= 0x800f;
	*q |= (e & 0x7ff) << 4;
	return(u.y);
	}
}

__ATTR_ALWAYS_INLINE__
static inline
double atan(double x)
{
static short P[20] = {
0x2594,0xa1f7,0x007f,0xbfec,
0x807a,0x5b6b,0x2854,0xc030,
0x0273,0x3688,0xc08c,0xc052,
0xba25,0x2d05,0xb8bf,0xc05e,
0xec8e,0xfd28,0x3669,0xc050,
};
static short Q[20] = {
/* 0x0000,0x0000,0x0000,0x3ff0, */
0x603c,0x5b14,0xdbc4,0x4038,
0xfa25,0x43b8,0xa0dd,0x4064,
0xbe3b,0xd2e2,0x0e18,0x407b,
0x49ea,0x13b0,0x563f,0x407e,
0x62ec,0xfbbd,0x519e,0x4068,
};

/* tan( 3*pi/8 ) = 2.41421356237309504880 */
static unsigned short T3P8A[] = {0x9de6,0x333f,0x504f,0x4003};
#define T3P8 *(double *)T3P8A

#define MOREBITS 6.123233995736765886130E-17

double y, z;
short sign, flag;

#ifdef MINUSZERO
if( x == 0.0 )
	return(x);
#endif
#ifdef INFINITIES
if(x == std::numeric_limits<double>::infinity())
	return(PIO2);
if(x == -std::numeric_limits<double>::infinity())
	return(-PIO2);
#endif
/* make argument positive and save the sign */
sign = 1;
if( x < 0.0 )
	{
	sign = -1;
	x = -x;
	}
/* range reduction */
flag = 0;
if( x > T3P8 )
	{
	y = PIO2;
	flag = 1;
	x = -( 1.0/x );
	}
else if( x <= 0.66 )
	{
	y = 0.0;
	}
else
	{
	y = PIO4;
	flag = 2;
	x = (x-1.0)/(x+1.0);
	}
z = x * x;
z = z * polevl( z, P, 4 ) / p1evl( z, Q, 5 );
z = x * z + x;
if( flag == 2 )
	z += 0.5 * MOREBITS;
else if( flag == 1 )
	z += MOREBITS;
y = y + z;
if( sign < 0 )
	y = -y;
return(y);
}
/////////////////////////////////////////////////////////////////////

__ATTR_ALWAYS_INLINE__
static inline
double exp(double x)
{
static unsigned short P[] = {
0x4be8,0xd5e4,0x89cd,0x3f20,
0x2c7e,0x0cca,0x06d1,0x3f9f,
0x0000,0x0000,0x0000,0x3ff0,
};
static unsigned short Q[] = {
0x5fa0,0xbc36,0x2eb6,0x3ec9,
0xb6c0,0xb508,0xae39,0x3f64,
0xe074,0x9887,0x1709,0x3fcd,
0x0000,0x0000,0x0000,0x4000,
};
static unsigned short sc1[] = {0x0000,0x0000,0x2e40,0x3fe6};
#define C1 (*(double *)sc1)
static unsigned short sc2[] = {0xabca,0xcf79,0xf7d1,0x3eb7};
#define C2 (*(double *)sc2)
double px, xx;
int n;

#ifdef NANS
if( isnan(x) )
	return(x);
#endif
if( x > MAXLOG)
	{
#ifdef INFINITIES
	return( std::numeric_limits<double>::infinity() );
#else
	return( MAXNUM );
#endif
	}

if( x < MINLOG )
	{
	return(0.0);
	}

/* Express e**x = e**g 2**n
 *   = e**g e**( n loge(2) )
 *   = e**( g + n loge(2) )
 */
px = floor( LOG2E * x + 0.5 ); /* floor() truncates toward -infinity. */
n = px;
x -= px * C1;
x -= px * C2;

/* rational approximation for exponential
 * of the fractional part:
 * e**x = 1 + 2x P(x**2)/( Q(x**2) - P(x**2) )
 */
xx = x * x;
px = x * polevl( xx, P, 2 );
x =  px/( polevl( xx, Q, 3 ) - px );
x = 1.0 + 2.0 * x;

/* multiply by power of 2 */
x = ldexp( x, n );
return(x);
}

__ATTR_ALWAYS_INLINE__
static inline
int sincos(double x, double *s, double *c, std::int32_t flg)
{

/* Define one of the following to be 1:
 */
#define ACC5 0
#define ACC11 0
#define ACC17 1

/* Option for linear interpolation when flg = 1
 */
#define LINTERP 1

/* Option for absolute error criterion
 */
#define ABSERR 1

/* Option to include modulo 360 function:
 */
#define MOD360 1

static double sintbl[92] = {
  0.00000000000000000000E0,
  1.74524064372835128194E-2,
  3.48994967025009716460E-2,
  5.23359562429438327221E-2,
  6.97564737441253007760E-2,
  8.71557427476581735581E-2,
  1.04528463267653471400E-1,
  1.21869343405147481113E-1,
  1.39173100960065444112E-1,
  1.56434465040230869010E-1,
  1.73648177666930348852E-1,
  1.90808995376544812405E-1,
  2.07911690817759337102E-1,
  2.24951054343864998051E-1,
  2.41921895599667722560E-1,
  2.58819045102520762349E-1,
  2.75637355816999185650E-1,
  2.92371704722736728097E-1,
  3.09016994374947424102E-1,
  3.25568154457156668714E-1,
  3.42020143325668733044E-1,
  3.58367949545300273484E-1,
  3.74606593415912035415E-1,
  3.90731128489273755062E-1,
  4.06736643075800207754E-1,
  4.22618261740699436187E-1,
  4.38371146789077417453E-1,
  4.53990499739546791560E-1,
  4.69471562785890775959E-1,
  4.84809620246337029075E-1,
  5.00000000000000000000E-1,
  5.15038074910054210082E-1,
  5.29919264233204954047E-1,
  5.44639035015027082224E-1,
  5.59192903470746830160E-1,
  5.73576436351046096108E-1,
  5.87785252292473129169E-1,
  6.01815023152048279918E-1,
  6.15661475325658279669E-1,
  6.29320391049837452706E-1,
  6.42787609686539326323E-1,
  6.56059028990507284782E-1,
  6.69130606358858213826E-1,
  6.81998360062498500442E-1,
  6.94658370458997286656E-1,
  7.07106781186547524401E-1,
  7.19339800338651139356E-1,
  7.31353701619170483288E-1,
  7.43144825477394235015E-1,
  7.54709580222771997943E-1,
  7.66044443118978035202E-1,
  7.77145961456970879980E-1,
  7.88010753606721956694E-1,
  7.98635510047292846284E-1,
  8.09016994374947424102E-1,
  8.19152044288991789684E-1,
  8.29037572555041692006E-1,
  8.38670567945424029638E-1,
  8.48048096156425970386E-1,
  8.57167300702112287465E-1,
  8.66025403784438646764E-1,
  8.74619707139395800285E-1,
  8.82947592858926942032E-1,
  8.91006524188367862360E-1,
  8.98794046299166992782E-1,
  9.06307787036649963243E-1,
  9.13545457642600895502E-1,
  9.20504853452440327397E-1,
  9.27183854566787400806E-1,
  9.33580426497201748990E-1,
  9.39692620785908384054E-1,
  9.45518575599316810348E-1,
  9.51056516295153572116E-1,
  9.56304755963035481339E-1,
  9.61261695938318861916E-1,
  9.65925826289068286750E-1,
  9.70295726275996472306E-1,
  9.74370064785235228540E-1,
  9.78147600733805637929E-1,
  9.81627183447663953497E-1,
  9.84807753012208059367E-1,
  9.87688340595137726190E-1,
  9.90268068741570315084E-1,
  9.92546151641322034980E-1,
  9.94521895368273336923E-1,
  9.96194698091745532295E-1,
  9.97564050259824247613E-1,
  9.98629534754573873784E-1,
  9.99390827019095730006E-1,
  9.99847695156391239157E-1,
  1.00000000000000000000E0,
  9.99847695156391239157E-1,
};

int ix, ssign, csign, xsign;
double y, z, sx, sz, cx, cz;

/* Make argument nonnegative.
 */
xsign = 1;
if( x < 0.0 )
	{
	xsign = -1;
	x = -x;
	}


#if MOD360
x = x  -  360.0 * floor( x/360.0 );
#endif

/* Find nearest integer to x.
 * Note there should be a domain error test here,
 * but this is omitted to gain speed.
 */
ix = x + 0.5;
z = x - ix;		/* the residual */

/* Look up the sine and cosine of the integer.
 */
if( ix <= 180 )
	{
	ssign = 1;
	csign = 1;
	}
else
	{
	ssign = -1;
	csign = -1;
	ix -= 180;
	}

if( ix > 90 )
	{
	csign = -csign;
	ix = 180 - ix;
	}

sx = sintbl[ix];
if( ssign < 0 )
	sx = -sx;
cx = sintbl[ 90-ix ];
if( csign < 0 )
	cx = -cx;

/* If the flag argument is set, then just return
 * the tabulated values for arg to the nearest whole degree.
 */
if( flg )
	{
#if LINTERP
	y = sx + 1.74531263774940077459e-2 * z * cx;
	cx -= 1.74531263774940077459e-2 * z * sx;
	sx = y;
#endif
	if( xsign < 0 )
		sx = -sx;
	*s = sx;	/* sine */
	*c = cx;	/* cosine */
	return 0;
	}

/* Find sine and cosine
 * of the residual angle between -0.5 and +0.5 degree.
 */
#if ACC5
#if ABSERR
/* absolute error = 2.769e-8: */
sz = 1.74531263774940077459e-2 * z;
/* absolute error = 4.146e-11: */
cz = 1.0 - 1.52307909153324666207e-4 * z * z;
#else
/* relative error = 6.346e-6: */
sz = 1.74531817576426662296e-2 * z;
/* relative error = 3.173e-6: */
cz = 1.0 - 1.52308226602566149927e-4 * z * z;
#endif
#else
y = z * z;
#endif


#if ACC11
sz = ( -8.86092781698004819918e-7 * y
      + 1.74532925198378577601e-2     ) * z;

cz = 1.0 - ( -3.86631403698859047896e-9 * y
            + 1.52308709893047593702e-4     ) * y;
#endif


#if ACC17
sz = ((  1.34959795251974073996e-11 * y
       - 8.86096155697856783296e-7     ) * y
       + 1.74532925199432957214e-2          ) * z;

cz = 1.0 - ((  3.92582397764340914444e-14 * y
             - 3.86632385155548605680e-9     ) * y
             + 1.52308709893354299569e-4          ) * y;
#endif


/* Combine the tabulated part and the calculated part
 * by trigonometry.
 */
y = sx * cz  +  cx * sz;
if( xsign < 0 )
	y = - y;
*s = y; /* sine */

*c = cx * cz  -  sx * sz; /* cosine */
return 0;
}



} // math

} // gms

#endif /*__GMS_CEPHES_DOUBLE_H__*/