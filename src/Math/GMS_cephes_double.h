
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
#include <cstdint>
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

namespace cephes_d
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
//const double INFINITY = std::numeric_limits<double>::infinity(); /* 99e999; */
#else
double INFINITY =  1.79769313486231570815E308;    /* 2**1024*(1-MACHEP) */
#endif
#ifdef NANS
//const double NAN = std::numeric_limits<double>::quiet_NaN();
#else
double NAN = 0.0;
#endif
#ifdef MINUSZERO
constexpr double NEGZERO = -0.0;
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

///////////////////////////////////////////////////////////////////
__ATTR_ALWAYS_INLINE__
static inline 
double polevl(double x,double coef[],std::int32_t N)
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

/*
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
*/
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
e += pw2;

/* Handle overflow */
if( e >= MEXP )
	return(MAXNUM );


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
z = z * polevl( z, (double*)&P, 4 ) / p1evl( z, (double*)&Q, 5 );
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
px = x * polevl( xx, (double*)&P, 2 );
x =  px/( polevl( xx, (double*)&Q, 3 ) - px );
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

__ATTR_ALWAYS_INLINE__
static inline
double chbevl(double x,double array[],std::int32_t n)
{
double b0, b1, b2, *p;
int i;
p = array;
b0 = *p++;
b1 = 0.0;
i = n - 1;
do
	{
	b2 = b1;
	b1 = b0;
	b0 = x * b1  -  b2  + *p++;
	}
while( --i );

return( 0.5*(b0-b2) );
}


__ATTR_ALWAYS_INLINE__
static inline
double sqrt(double x)
{
int e;
short *q;
double z, w;

if( x <= 0.0 )
	{
	if( x < 0.0 )
	   return( 0.0 );
	}
w = x;
/* separate exponent and significand */


/* Note, frexp and ldexp are used in order to
 * handle denormal numbers properly.
 */

z = frexp( x, &e );
q = (short *)&x;
q += 3;
/*
e = ((*q >> 4) & 0x0fff) - 0x3fe;
*q &= 0x000f;
*q |= 0x3fe0;
z = x;
*/

/* approximate square root of number between 0.5 and 1
 * relative error of approximation = 7.47e-3
 */
x = 4.173075996388649989089E-1 + 5.9016206709064458299663E-1 * z;

/* adjust for odd powers of 2 */
if( (e & 1) != 0 )
	x *= SQRT2;

/* re-insert exponent */

x = ldexp( x, (e >> 1) );
/*
*q += ((e >>1) & 0x7ff) << 4;
*q &= 077777;
*/

/* Newton iterations: */

/* Note, assume the square root cannot be denormal,
 * so it is safe to use integer exponent operations here.
 */

x += w/x;
*q -= 0x10;
x += w/x;
*q -= 0x10;
x += w/x;
*q -= 0x10;
return(x);
}

__ATTR_ALWAYS_INLINE__
static inline
double i0(double x)
{
static unsigned short A[] = {
0xd0ef,0x2134,0x5cb7,0xbc54,
0xa589,0x977d,0x3362,0x3c83,
0xbbb4,0x721e,0x84eb,0xbcb1,
0x5eba,0x93f6,0xe6d8,0x3cde,
0xfbeb,0xc297,0x5022,0xbd0a,
0x2627,0x4b26,0x9b46,0x3d35,
0x1af0,0x62ee,0x164c,0xbd61,
0xd324,0xe19b,0xfe2f,0x3d89,
0x6abc,0x7a94,0xfc95,0xbdb2,
0x3c10,0xcc74,0x98be,0x3dda,
0x9556,0x13ae,0xd4fe,0xbe01,
0xcb34,0xa454,0xd903,0x3e26,
0x30ab,0x8c0b,0xeaf6,0xbe4b,
0x6435,0x9d4d,0x3b76,0x3e70,
0x7f8d,0x8f22,0xec63,0xbe91,
0xf4ac,0x978c,0xbf24,0x3eb2,
0x6427,0xcba5,0x866f,0xbed2,
0x2859,0xbe9a,0x3f58,0x3ef1,
0x1d5a,0x59c4,0x2b26,0xbf0e,
0x7cab,0x7410,0xb51b,0x3f28,
0xeb52,0x1f15,0xe2fd,0xbf42,
0x100e,0x8a12,0xdc75,0x3f5a,
0xa849,0x201a,0xb65e,0xbf71,
0xe3dd,0xf3dd,0x9961,0x3f85,
0xb6f0,0xf121,0x4e9e,0xbf98,
0xa32d,0xcea8,0x3e8a,0x3fa9,
0x06ea,0x342d,0x4b70,0xbfb8,
0x88c0,0x77ac,0xf7ac,0x3fc5,
0xcd8d,0xc057,0x7feb,0xbfd3,
0xa22a,0x9035,0xa84e,0x3fe5,
};
static unsigned short B[] = {
0x8b19,0x54ca,0xadb7,0xbc60,
0x9130,0x6611,0x46da,0xbc56,
0x8421,0x12d9,0xbe18,0x3c89,
0x41cd,0x0760,0xf3dd,0x3c83,
0x1fe4,0xabd2,0x600b,0xbcb4,
0xde38,0xd908,0xaee7,0xbcb8,
0xfb1f,0xa3ea,0xee7d,0x3cdf,
0xe6d7,0x9094,0x2a91,0x3cf1,
0x629a,0x7e65,0x83fe,0xbd05,
0xbb32,0xcf68,0x5d99,0xbd27,
0xc545,0x0d5f,0x56ff,0x3d11,
0xc073,0x6b83,0x1c8c,0x3d5b,
0x8cec,0xfa26,0x4347,0x3d69,
0x8d66,0x0317,0x9043,0xbd7f,
0x7bf2,0x357e,0x0fd7,0xbdad,
0x7425,0x0839,0x511d,0xbdc1,
0x004f,0xabe8,0x24fe,0x3daa,
0x6f75,0xc0f4,0xf9cc,0x3e00,
0x5b87,0xa922,0x2c64,0x3e2d,
0xd56d,0x80d6,0x5692,0x3e58,
0x616e,0xd9cd,0x8007,0x3e8b,
0xc586,0xc101,0x412b,0x3ec8,
0x9e52,0x7899,0x0fa3,0x3f12,
0x9049,0xa2e5,0x998c,0x3f6b,
0x09cb,0xaca8,0xbe62,0x3fe9
};
double y;

if( x < 0 )
	x = -x;
if( x <= 8.0 )
	{
	y = (x/2.0) - 2.0;
	return( exp(x) * chbevl( y, (double*)&A, 30 ) );
	}

return(  exp(x) * chbevl( 32.0/x - 2.0, (double*)&B, 25 ) / sqrt(x) );

}


__ATTR_ALWAYS_INLINE__
static inline
double i0e( double x )
{
static unsigned short A[] = {
0xd0ef,0x2134,0x5cb7,0xbc54,
0xa589,0x977d,0x3362,0x3c83,
0xbbb4,0x721e,0x84eb,0xbcb1,
0x5eba,0x93f6,0xe6d8,0x3cde,
0xfbeb,0xc297,0x5022,0xbd0a,
0x2627,0x4b26,0x9b46,0x3d35,
0x1af0,0x62ee,0x164c,0xbd61,
0xd324,0xe19b,0xfe2f,0x3d89,
0x6abc,0x7a94,0xfc95,0xbdb2,
0x3c10,0xcc74,0x98be,0x3dda,
0x9556,0x13ae,0xd4fe,0xbe01,
0xcb34,0xa454,0xd903,0x3e26,
0x30ab,0x8c0b,0xeaf6,0xbe4b,
0x6435,0x9d4d,0x3b76,0x3e70,
0x7f8d,0x8f22,0xec63,0xbe91,
0xf4ac,0x978c,0xbf24,0x3eb2,
0x6427,0xcba5,0x866f,0xbed2,
0x2859,0xbe9a,0x3f58,0x3ef1,
0x1d5a,0x59c4,0x2b26,0xbf0e,
0x7cab,0x7410,0xb51b,0x3f28,
0xeb52,0x1f15,0xe2fd,0xbf42,
0x100e,0x8a12,0xdc75,0x3f5a,
0xa849,0x201a,0xb65e,0xbf71,
0xe3dd,0xf3dd,0x9961,0x3f85,
0xb6f0,0xf121,0x4e9e,0xbf98,
0xa32d,0xcea8,0x3e8a,0x3fa9,
0x06ea,0x342d,0x4b70,0xbfb8,
0x88c0,0x77ac,0xf7ac,0x3fc5,
0xcd8d,0xc057,0x7feb,0xbfd3,
0xa22a,0x9035,0xa84e,0x3fe5,
};
static unsigned short B[] = {
0x8b19,0x54ca,0xadb7,0xbc60,
0x9130,0x6611,0x46da,0xbc56,
0x8421,0x12d9,0xbe18,0x3c89,
0x41cd,0x0760,0xf3dd,0x3c83,
0x1fe4,0xabd2,0x600b,0xbcb4,
0xde38,0xd908,0xaee7,0xbcb8,
0xfb1f,0xa3ea,0xee7d,0x3cdf,
0xe6d7,0x9094,0x2a91,0x3cf1,
0x629a,0x7e65,0x83fe,0xbd05,
0xbb32,0xcf68,0x5d99,0xbd27,
0xc545,0x0d5f,0x56ff,0x3d11,
0xc073,0x6b83,0x1c8c,0x3d5b,
0x8cec,0xfa26,0x4347,0x3d69,
0x8d66,0x0317,0x9043,0xbd7f,
0x7bf2,0x357e,0x0fd7,0xbdad,
0x7425,0x0839,0x511d,0xbdc1,
0x004f,0xabe8,0x24fe,0x3daa,
0x6f75,0xc0f4,0xf9cc,0x3e00,
0x5b87,0xa922,0x2c64,0x3e2d,
0xd56d,0x80d6,0x5692,0x3e58,
0x616e,0xd9cd,0x8007,0x3e8b,
0xc586,0xc101,0x412b,0x3ec8,
0x9e52,0x7899,0x0fa3,0x3f12,
0x9049,0xa2e5,0x998c,0x3f6b,
0x09cb,0xaca8,0xbe62,0x3fe9
};

double y;

if( x < 0 )
	x = -x;
if( x <= 8.0 )
	{
	y = (x/2.0) - 2.0;
	return( chbevl( y, (double*)&A, 30 ) );
	}

return(  chbevl( 32.0/x - 2.0, (double*)&B, 25 ) / sqrt(x) );
}

__ATTR_ALWAYS_INLINE__
static inline
double powi(double x, std::int32_t nn )
{
int n, e, sign, asign, lx;
double w, y, s;

/* See pow.c for these tests.  */
if( x == 0.0 )
	{
	if( nn == 0 )
		return( 1.0 );
	else if( nn < 0 )
	    return( std::numeric_limits<double>::infinity() );
	else
	  {
	    if( nn & 1 )
	      return( x );
	    else
	      return( 0.0 );
	  }
	}

if( nn == 0 )
	return( 1.0 );

if( nn == -1 )
	return( 1.0/x );

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
	}
else
	{
	sign = 1;
	n = nn;
	}

/* Even power will be positive. */
if( (n & 1) == 0 )
	asign = 0;

/* Overflow detection */

/* Calculate approximate logarithm of answer */
s = frexp( x, &lx );
e = (lx - 1)*n;
if( (e == 0) || (e > 64) || (e < -64) )
	{
	s = (s - 7.0710678118654752e-1) / (s +  7.0710678118654752e-1);
	s = (2.9142135623730950 * s - 0.5 + lx) * nn * LOGE2;
	}
else
	{
	s = LOGE2 * e;
	}

if( s > MAXLOG )
	{
	y = std::numeric_limits<double>::infinity();
	goto done;
	}

#if DENORMAL
if( s < MINLOG )
	{
	y = 0.0;
	goto done;
	}

/* Handle tiny denormal answer, but with less accuracy
 * since roundoff error in 1.0/x will be amplified.
 * The precise demarcation should be the gradual underflow threshold.
 */
if( (s < (-MAXLOG+2.0)) && (sign < 0) )
	{
	x = 1.0/x;
	sign = -sign;
	}
#else
/* do not produce denormal answer */
if( s < -MAXLOG )
	return(0.0);
#endif


/* First bit of the power */
if( n & 1 )
	y = x;
		
else
	y = 1.0;

w = x;
n >>= 1;
while( n )
	{
	w = w * w;	/* arg to the 2-to-the-kth power */
	if( n & 1 )	/* if that bit is set, then include in product */
		y *= w;
	n >>= 1;
	}

if( sign < 0 )
	y = 1.0/y;

done:

if( asign )
	{
	/* odd power of negative number */
	if( y == 0.0 )
		y = NEGZERO;
	else
		y = -y;
	}
return(y);
}

__ATTR_ALWAYS_INLINE__
static inline 
double reduc(double x)
{
double t;

t = ldexp( x, 4 );
t = floor( t );
t = ldexp( t, -4 );
return(t);
}

__ATTR_ALWAYS_INLINE__
static inline
double pow(double x,double y)
{

static unsigned short P[] = {
0x5cf0,0x7f5b,0xdb99,0x3fdf,
0xdf15,0xea9e,0xddef,0x400d,
0xeb6f,0x7f78,0xccbd,0x401e,
0x9b74,0xb65c,0xaa83,0x4012,
};
static unsigned short Q[] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0x914e,0x9b20,0xaab4,0x4022,
0xc9f5,0x41c1,0xffff,0x403b,
0x6402,0x1b17,0xccbc,0x4040,
0xe92e,0x918a,0xffc5,0x402b,
};
static unsigned short A[] = {
0x0000,0x0000,0x0000,0x3ff0,
0x90da,0xa2a4,0xa4af,0x3fee,
0xa487,0xdcfb,0x5818,0x3fed,
0x529c,0xdd85,0x199b,0x3fec,
0xd3ad,0x995a,0xe89f,0x3fea,
0xf090,0x82a3,0xc491,0x3fe9,
0xa0db,0x422a,0xace5,0x3fe8,
0x0187,0x73eb,0xa114,0x3fe7,
0x3bcd,0x667f,0xa09e,0x3fe6,
0x5429,0xdd48,0xab07,0x3fe5,
0x2a27,0xd536,0xbfda,0x3fe4,
0x3422,0x4c12,0xdea6,0x3fe3,
0xb715,0x0a31,0x06fe,0x3fe3,
0x6238,0x6e75,0x387a,0x3fe2,
0x517b,0x3c7d,0x72b8,0x3fe1,
0x890f,0x6cf9,0xb558,0x3fe0,
0x0000,0x0000,0x0000,0x3fe0
};
static unsigned short B[] = {
0x0000,0x0000,0x0000,0x0000,
0x3707,0xd75b,0xed02,0x3c72,
0xcc81,0x345d,0xa1cd,0x3c87,
0x4b27,0x5686,0xe9f1,0x3c86,
0x6456,0x13b2,0xdd34,0xbc8b,
0x42e2,0xafec,0x4397,0x3c6d,
0x82e4,0xd231,0xf46a,0x3c76,
0x8a76,0xb9d7,0x9041,0xbc71,
0x0000,0x0000,0x0000,0x0000
};
static unsigned short R[] = {
0x937f,0xd7f2,0x6307,0x3eef,
0x9259,0x60fc,0x2fbe,0x3f24,
0xef1d,0xc84a,0xd87e,0x3f55,
0x33b7,0x6ef1,0xb2ab,0x3f83,
0x1a92,0xd704,0x6b08,0x3fac,
0xc56d,0xff82,0xbfbd,0x3fce,
0x39ef,0xfefa,0x2e42,0x3fe6
};
#define SQRTH 0.70710678118654752440
#define douba(k) (*(double *)&A[(k)<<2])
#define doubb(k) (*(double *)&B[(k)<<2])
//#define MEXP 16383.0
#ifdef DENORMAL
#define MNEXP -17183.0
#else
#define MNEXP -16383.0
#endif

double w, z, W, Wa, Wb, ya, yb, u;
double F, Fa, Fb, G, Ga, Gb, H, Ha, Hb;
double aw, ay, wy;
int e, i, nflg, iyflg, yoddint;

if( y == 0.0 )
	return( 1.0 );
#ifdef NANS
if( isnan(x) )
	return( x );
if( isnan(y) )
	return( y );
#endif
if( y == 1.0 )
	return( x );


#ifdef INFINITIES
if( !isfinite(y) && (x == 1.0 || x == -1.0) )
	{
#ifdef NANS
	return( std::numeric_limits<double>::quiet_NaN());
#else
	return( std::numeric_limits<double>::infinity() );
#endif
	}
#endif

if( x == 1.0 )
	return( 1.0 );

if( y >= MAXNUM )
	{
#ifdef INFINITIES
	if( x > 1.0 )
		return( std::numeric_limits<double>::infinity() );
#else
	if( x > 1.0 )
		return( MAXNUM );
#endif
	if( x > 0.0 && x < 1.0 )
		return( 0.0);
	if( x < -1.0 )
		{
#ifdef INFINITIES
		return( std::numeric_limits<double>::infinity() );
#else
		return( MAXNUM );
#endif
		}
	if( x > -1.0 && x < 0.0 )
		return( 0.0 );
	}
if( y <= -MAXNUM )
	{
	if( x > 1.0 )
		return( 0.0 );
#ifdef INFINITIES
	if( x > 0.0 && x < 1.0 )
		return( std::numeric_limits<double>::infinity() );
#else
	if( x > 0.0 && x < 1.0 )
		return( MAXNUM );
#endif
	if( x < -1.0 )
		return( 0.0 );
#ifdef INFINITIES
	if( x > -1.0 && x < 0.0 )
		return( std::numeric_limits<double>::infinity() );
#else
	if( x > -1.0 && x < 0.0 )
		return( MAXNUM );
#endif
	}
if( x >= MAXNUM )
	{
#if INFINITIES
	if( y > 0.0 )
		return( std::numeric_limits<double>::infinity() );
#else
	if( y > 0.0 )
		return( MAXNUM );
#endif
	return(0.0);
	}
/* Set iyflg to 1 if y is an integer.  */
iyflg = 0;
w = floor(y);
if( w == y )
	iyflg = 1;

/* Test for odd integer y.  */
yoddint = 0;
if( iyflg )
	{
	ya = fabs(y);
	ya = floor(0.5 * ya);
	yb = 0.5 * fabs(w);
	if( ya != yb )
		yoddint = 1;
	}

if( x <= -MAXNUM )
	{
	if( y > 0.0 )
		{
#ifdef INFINITIES
		if( yoddint )
			return( -std::numeric_limits<double>::infinity() );
		return( std::numeric_limits<double>::infinity() );
#else
		if( yoddint )
			return( -MAXNUM );
		return( MAXNUM );
#endif
		}
	if( y < 0.0 )
		{
#ifdef MINUSZERO
		if( yoddint )
			return( NEGZERO );
#endif
		return( 0.0 );
		}
 	}

nflg = 0;	/* flag = 1 if x<0 raised to integer power */
if( x <= 0.0 )
	{
	if( x == 0.0 )
		{
		if( y < 0.0 )
			{
#ifdef MINUSZERO
			if( signbit(x) && yoddint )
				return( -std::numeric_limits<double>::infinity() );
#endif
#ifdef INFINITIES
			return( std::numeric_limits<double>::infinity() );
#else
			return( MAXNUM );
#endif
			}
		if( y > 0.0 )
			{
#ifdef MINUSZERO
			if( signbit(x) && yoddint )
				return( NEGZERO );
#endif
			return( 0.0 );
			}
		return( 1.0 );
		}
	else
		{
		if( iyflg == 0 )
			{ /* noninteger power of negative number */
#ifdef NANS
			return(std::numeric_limits<double>::quiet_NaN());
#else
			return(0.0L);
#endif
			}
		nflg = 1;
		}
	}

/* Integer power of an integer.  */

if( iyflg )
	{
	i = w;
	w = floor(x);
	if( (w == x) && (fabs(y) < 32768.0) )
		{
		w = powi( x, (int) y );
		return( w );
		}
	}

if( nflg )
	x = fabs(x);

/* For results close to 1, use a series expansion.  */
w = x - 1.0;
aw = fabs(w);
ay = fabs(y);
wy = w * y;
ya = fabs(wy);
if((aw <= 1.0e-3 && ay <= 1.0)
   || (ya <= 1.0e-3 && ay >= 1.0))
	{
	z = (((((w*(y-5.)/720. + 1./120.)*w*(y-4.) + 1./24.)*w*(y-3.)
		+ 1./6.)*w*(y-2.) + 0.5)*w*(y-1.) )*wy + wy + 1.;
	goto done;
	}
/* These are probably too much trouble.  */
#if 0
w = y * log(x);
if (aw > 1.0e-3 && fabs(w) < 1.0e-3)
  {
    z = ((((((
    w/7. + 1.)*w/6. + 1.)*w/5. + 1.)*w/4. + 1.)*w/3. + 1.)*w/2. + 1.)*w + 1.;
    goto done;
  }

if(ya <= 1.0e-3 && aw <= 1.0e-4)
  {
    z = (((((
	     wy*1./720.
	     + (-w*1./48. + 1./120.) )*wy
	    + ((w*17./144. - 1./12.)*w + 1./24.) )*wy
	   + (((-w*5./16. + 7./24.)*w - 1./4.)*w + 1./6.) )*wy
	  + ((((w*137./360. - 5./12.)*w + 11./24.)*w - 1./2.)*w + 1./2.) )*wy
	 + (((((-w*1./6. + 1./5.)*w - 1./4)*w + 1./3.)*w -1./2.)*w ) )*wy
	   + wy + 1.0;
    goto done;
  }
#endif

/* separate significand from exponent */
x = frexp( x, &e );

#if 0
/* For debugging, check for gross overflow. */
if( (e * y)  > (16383.0 + 1024) )
	goto overflow;
#endif

/* Find significand of x in antilog table A[]. */
i = 1;
if( x <= douba(9) )
	i = 9;
if( x <= douba(i+4) )
	i += 4;
if( x <= douba(i+2) )
	i += 2;
if( x >= douba(1) )
	i = -1;
i += 1;


/* Find (x - A[i])/A[i]
 * in order to compute log(x/A[i]):
 *
 * log(x) = log( a x/a ) = log(a) + log(x/a)
 *
 * log(x/a) = log(1+v),  v = x/a - 1 = (x-a)/a
 */
x -= douba(i);
x -= doubb(i/2);
x /= douba(i);


/* rational approximation for log(1+v):
 *
 * log(1+v)  =  v  -  v**2/2  +  v**3 P(v) / Q(v)
 */
z = x*x;
w = x * ( z * polevl( x, (double*)&P, 3 ) / p1evl( x, (double*)&Q, 4 ) );
w = w - ldexp( z, -1 );   /*  w - 0.5 * z  */
#define LOG2EA 0.44269504088896340736
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
w = ldexp( w, -4 );	/* divide by 16 */
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
w = ldexp( Ga+Ha, 4 );

/* Test the power of 2 for overflow */
if( w > 16383.0 )
	{
#ifdef INFINITIES
	if( nflg && yoddint )
	  return( -std::numeric_limits<double>::infinity() );
	return( std::numeric_limits<double>::infinity() );
#else
	if( nflg && yoddint )
	  return( -MAXNUM );
	return( MAXNUM );
#endif
	}

if( w < (MNEXP - 1) )
	{
#ifdef MINUSZERO
	if( nflg && yoddint )
	  return( NEGZERO );
#endif
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
 */
z = Hb * polevl( Hb, (double*)&R, 6 );  /*    z  =  2**Hb - 1    */

/* Express e/16 as an integer plus a negative number of 16ths.
 * Find lookup table entry for the fractional power of 2.
 */
if( e < 0 )
	i = 0;
else
	i = 1;
i = e/16 + i;
e = 16*i - e;
w = douba( e );
z = w + w * z;      /*    2**-e * ( 1 + (2**Hb-1) )    */
z = ldexp( z, i );  /* multiply by integer power of 2 */

done:

/* Negate if odd integer power of negative number */
if( nflg && yoddint )
	{
#ifdef MINUSZERO
	if( z == 0.0 )
		z = NEGZERO;
	else
#endif
		z = -z;
	}
return( z );
}

__ATTR_ALWAYS_INLINE__
static inline
double sin(double x)
{
constexpr unsigned short sincof[] = {
0x9ccd,0x1fd1,0xd8fd,0x3de5,
0x1f5d,0xa929,0xe5e5,0xbe5a,
0x48a1,0x567d,0x1de3,0x3ec7,
0xdf03,0x19bf,0x01a0,0xbf2a,
0xf7d0,0x1110,0x1111,0x3f81,
0x5548,0x5555,0x5555,0xbfc5,
};
constexpr unsigned short coscof[24] = {
0x1a9b,0xa086,0xfa49,0xbda8,
0x3f05,0x7b4e,0xee9d,0x3e21,
0x4bc6,0x7eac,0x7e4f,0xbe92,
0x44f5,0x19c8,0x01a0,0x3efa,
0x4f91,0x16c1,0xc16c,0xbf56,
0x554b,0x5555,0x5555,0x3fa5,
};
constexpr unsigned short P1[] = {0x0000,0x4000,0x21fb,0x3fe9};
constexpr unsigned short P2[] = {0x0000,0x0000,0x442d,0x3e64};
constexpr unsigned short P3[] = {0x5170,0x98cc,0x4698,0x3ce8};
#define DP1 *(double *)P1
#define DP2 *(double *)P2
#define DP3 *(double *)P3
constexpr double lossth = 1.073741824e9;
double y, z, zz;
int j, sign;

#ifdef MINUSZERO
if( x == 0.0 )
	return(x);
#endif
#ifdef NANS
if( isnan(x) )
	return(x);
if( !isfinite(x) )
	{
	return(std::numeric_limits<double>::quiet_NaN());
	}
#endif
/* make argument positive but save the sign */
sign = 1;
if( x < 0 )
	{
	x = -x;
	sign = -1;
	}

if( x > lossth )
	{
	return(0.0);
	}

y = floor( x/PIO4 ); /* integer part of x/PIO4 */

/* strip high bits of integer part to prevent integer overflow */
z = ldexp( y, -4 );
z = floor(z);           /* integer part of y/8 */
z = y - ldexp( z, 4 );  /* y - 16 * (y/16) */

j = z; /* convert to integer for tests on the phase angle */
/* map zeros to origin */
if( j & 1 )
	{
	j += 1;
	y += 1.0;
	}
j = j & 07; /* octant modulo 360 degrees */
/* reflect in x axis */
if( j > 3)
	{
	sign = -sign;
	j -= 4;
	}

/* Extended precision modular arithmetic */
z = ((x - y * DP1) - y * DP2) - y * DP3;

zz = z * z;

if( (j==1) || (j==2) )
	{
	y = 1.0 - ldexp(zz,-1) + zz * zz * polevl( zz, (double*)&coscof, 5 );
	}
else
	{
/*	y = z  +  z * (zz * polevl( zz, sincof, 5 ));*/
	y = z  +  z * z * z * polevl( zz, (double*)&sincof, 5 );
	}

if(sign < 0)
	y = -y;

return(y);
}

__ATTR_ALWAYS_INLINE__
static inline
double cos(double x)
{
constexpr unsigned short sincof[] = {
0x9ccd,0x1fd1,0xd8fd,0x3de5,
0x1f5d,0xa929,0xe5e5,0xbe5a,
0x48a1,0x567d,0x1de3,0x3ec7,
0xdf03,0x19bf,0x01a0,0xbf2a,
0xf7d0,0x1110,0x1111,0x3f81,
0x5548,0x5555,0x5555,0xbfc5,
};
constexpr unsigned short coscof[24] = {
0x1a9b,0xa086,0xfa49,0xbda8,
0x3f05,0x7b4e,0xee9d,0x3e21,
0x4bc6,0x7eac,0x7e4f,0xbe92,
0x44f5,0x19c8,0x01a0,0x3efa,
0x4f91,0x16c1,0xc16c,0xbf56,
0x554b,0x5555,0x5555,0x3fa5,
};
constexpr unsigned short P1[] = {0x0000,0x4000,0x21fb,0x3fe9};
constexpr unsigned short P2[] = {0x0000,0x0000,0x442d,0x3e64};
constexpr unsigned short P3[] = {0x5170,0x98cc,0x4698,0x3ce8};
#define DP1 *(double *)P1
#define DP2 *(double *)P2
#define DP3 *(double *)P3
constexpr double lossth = 1.073741824e9;
double y, z, zz;
long i;
int j, sign;

#ifdef NANS
if( isnan(x) )
	return(x);
if( !isfinite(x) )
	{
	return(std::numeric_limits<double>::quiet_NaN());
	}
#endif

/* make argument positive */
sign = 1;
if( x < 0 )
	x = -x;

if( x > lossth )
	{
	return(0.0);
	}

y = floor( x/PIO4 );
z = ldexp( y, -4 );
z = floor(z);		/* integer part of y/8 */
z = y - ldexp( z, 4 );  /* y - 16 * (y/16) */

/* integer and fractional part modulo one octant */
i = z;
if( i & 1 )	/* map zeros to origin */
	{
	i += 1;
	y += 1.0;
	}
j = i & 07;
if( j > 3)
	{
	j -=4;
	sign = -sign;
	}

if( j > 1 )
	sign = -sign;

/* Extended precision modular arithmetic */
z = ((x - y * DP1) - y * DP2) - y * DP3;

zz = z * z;

if( (j==1) || (j==2) )
	{
/*	y = z  +  z * (zz * polevl( zz, sincof, 5 ));*/
	y = z  +  z * z * z * polevl( zz, (double*)&sincof, 5 );
	}
else
	{
	y = 1.0 - ldexp(zz,-1) + zz * zz * polevl( zz, (double*)&coscof, 5 );
	}

if(sign < 0)
	y = -y;

return(y);
}





/* Degrees, minutes, seconds to radians: */

/* 1 arc second, in radians = 4.8481368110953599358991410e-5 */
__ATTR_ALWAYS_INLINE__
static inline
double radian(double d,double m,double s)
{
constexpr double P64800 = 4.8481368110953599358991410e-5;
return( ((d*60.0 + m)*60.0 + s)*P64800 );
}

/* Find a multiple of 1/16 that is within 1/16 of x. */



} // cephes_d

} // math

} // gms

#endif /*__GMS_CEPHES_DOUBLE_H__*/