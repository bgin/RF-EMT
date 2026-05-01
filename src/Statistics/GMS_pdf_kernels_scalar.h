
/* Copyright (C) J. Racine, 1995-2001 
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

/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef __GMS_PDF_KERNELS_SCALAR_H__
#define __GMS_PDF_KERNELS_SCALAR_H__ 010520260831 

namespace file_info 
{

     const unsigned int GMS_PDF_KERNELS_SCALAR_MAJOR = 1;
     const unsigned int GMS_PDF_KERNELS_SCALAR_MINOR = 0;
     const unsigned int GMS_PDF_KERNELS_SCALAR_MICRO = 0;
     const unsigned int GMS_PDF_KERNELS_SCALAR_FULLVER =
       1000U*GMS_PDF_KERNELS_SCALAR_MAJOR+100U*GMS_PDF_KERNELS_SCALAR_MINOR+
       10U*GMS_PDF_KERNELS_SCALAR_MICRO;
     const char * const GMS_PDF_KERNELS_SCALAR_CREATION_DATE = "01-05-2026 08:31 +00200 (FRI 01 MAY 2026 GMT+2)";
     const char * const GMS_PDF_KERNELS_SCALAR_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_PDF_KERNELS_SCALAR_SYNOPSIS      = "The Density estimation kernel functions.Original implementation by J. Raccine modified by \
                                                                Bernard Gingold.";
}

#include <cstdint>
#include <cmath>
#include "GMS_config.h"

namespace math 
{

namespace np_standalone_funcs
{
    // Anonymous namespace for exponential function implementation 
namespace 
{

__ATTR_ALWAYS_INLINE__
static inline
float ceph_floorf( float x ) {
constexpr std::int32_t NBITS 24;
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
float ceph_ldexpf( float x, int pw2 ) {
#define DENORMAL 1
constexpr std::int32_t MEXP 255;
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
float ceph_expf(const float xx) {
float x, z;
int n;
x = xx;
if( x > 88.72283905206835f) {
    return(3.4028234663852885981170418348451692544e38f);
}
if( x < -103.278929903431851103f) {
    return(0.0);
}

/* Express e**x = e**g 2**n
 *   = e**g e**( n loge(2) )
 *   = e**( g + n loge(2) )
 */
z = ceph_floorf(1.44269504088896341f * x + 0.5 ); /* floor() truncates toward -infinity. */
x -= z * 0.693359375f;
x -= z * -2.12194440e-4f;
n = z;

z = x * x;
/* Theoretical peak relative error in [-0.5, +0.5] is 4.2e-9. */
z =
((((( 1.9875691500E-4f  * x
   + 1.3981999507E-3f) * x
   + 8.3334519073E-3f) * x
   + 4.1665795894E-2f) * x
   + 1.6666665459E-1f) * x
   + 5.0000001201E-1f) * z
   + x
   + 1.0;

/* multiply by power of 2 */
x = ceph_ldexpf( z, n );
return( x );
}

} // anonymous namespace

__ATTR_ALWAYS_INLINE__
static inline 
float  np_gauss2_r4(const float z)
{   
    constexpr float C0398942280401432677939946059934{0.398942280401432677939946059934f};
    return C0398942280401432677939946059934*ceph_expf(-0.5f*z*z);
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_gauss4_r4(const float z)
{
     constexpr float C0398942280401432677939946059934{0.398942280401432677939946059934f};
     const float exp_arg{0.5f*z*z};
     return C0398942280401432677939946059934*(1.5f-exp_arg)*ceph_expf(-exp_arg);
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_gauss6_r4(const float z)
{
     constexpr float C0398942280401432677939946059934{0.398942280401432677939946059934f};
     const float z2{z*z};
     const float term1{z2*0.125f-1.25f};
     return C0398942280401432677939946059934*std::fma(z2,term1,1.875)*ceph_expf(-0.5f*z2);
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_gauss8_r4(const float z)
{
    constexpr float C0398942280401432677939946059934{0.398942280401432677939946059934f};
    const float z2{z*z};
    const float term1{0.4375f-z2*0.02083333333f};
    return C0398942280401432677939946059934*std::fma(z2,std::fma(z2,term1,-2.1875f),2.1875f)*ceph_expf(-0.5f*z2);
}

} //np_standalone_funcs

} // math
















#endif /*__GMS_PDF_KERNELS_SCALAR_H__*/