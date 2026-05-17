
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
#if (USE_OPENMP) == 1
#include <omp.h>
#endif
#include "GMS_config.h"
#include "GMS_compare_fp_safe.h"

namespace gms
{

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

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float  np_gauss2_r4(const float z)
{   
    constexpr float C0398942280401432677939946059934{0.398942280401432677939946059934f};
    return C0398942280401432677939946059934*ceph_expf(-0.5f*z*z);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_gauss4_r4(const float z)
{
     constexpr float C0398942280401432677939946059934{0.398942280401432677939946059934f};
     const float exp_arg{0.5f*z*z};
     return C0398942280401432677939946059934*(1.5f-exp_arg)*ceph_expf(-exp_arg);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_gauss6_r4(const float z)
{
     constexpr float C0398942280401432677939946059934{0.398942280401432677939946059934f};
     const float z2{z*z};
     const float term1{z2*0.125f-1.25f};
     return C0398942280401432677939946059934*std::fma(z2,term1,1.875)*ceph_expf(-0.5f*z2);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_gauss8_r4(const float z)
{
    constexpr float C0398942280401432677939946059934{0.398942280401432677939946059934f};
    const float z2{z*z};
    const float term1{0.4375f-z2*0.02083333333f};
    return C0398942280401432677939946059934*std::fma(z2,std::fma(z2,term1,-2.1875f),2.1875f)*ceph_expf(-0.5f*z2);
}

/*
   For Epanechnikov kernels the input check i.e. z^2 < 5.0 was removed 
   in order to remove the branch-misprediction.
*/

/*
    The input range is (z*z < 5.0)?result:0.0
*/
#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_epanechnikov2_r4(const float z)
{
    const float zz{z*z};
    return (static_cast<float>(0.33541019662496845446f-0.067082039324993690892f*zz));
   
}

/*
    The input range is (z*z < 5.0)?result:0.0
*/
#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_epanechnikov4_r4(const float z)
{
    const float zz{z*z};
    return (static_cast<float>(0.008385254916f*std::fma(7.0f,zz,-15.0f)*(-5.0+zz)));
    
}

/*
    The input range is (z*z < 5.0)?result:0.0
*/
#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_epanechnikov6_r4(const float z)
{
    const float zz{z*z};
    const float term1{static_cast<float>(0.33541019662496845446f*std::fma(std::fma(0.721875f,zz,-3.28125f),zz,2.734375f)*
                                         (1.0f-0.2f*zz))}; 
    return (term1);
   
                  
}

/*
    The input range is (z*z < 5.0)?result:0.0
*/
#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_epanechnikov8_r4(const float z)
{
    const float zz{z*z};
    const float fma_arg{4.1056640625f-0.5865234375f*zz};
    const float fma_term{std::fma(std::fma(fma_arg,zz,-7.8955078125f),zz,3.5888671875f)};
    return (static_cast<float>(0.33541019662496845446f*fma_term*(1.0f-0.2f*zz)));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_rect_r4(const float z)
{
    return (z*z < 1.0f)?0.5f:0.0f;
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_uaa_r4(const bool same_cat,const float lambda, const std::int32_t c)
{
    return (same_cat)?(1.0f-lambda):lambda/((static_cast<float>(c)-1.0f));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_score_uaa_r4(const bool same_cat,const float lambda, const std::int32_t c)
{
    return (same_cat)?-1.0f:(1.0f/(static_cast<float>(c)-1.0f)); 
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_uli_racine_r4(const bool same_cat,const float lambda)
{
    return (same_cat)?1.0f:lambda;
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_unli_racine_r4(const bool same_cat,const float lambda, const std::int32_t c)
{
     return ((same_cat)?1.0f:lambda)/((static_cast<float>(c)-1.0f)*lambda + 1.0f);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_score_uli_racine_r4(const bool same_cat,const float lambda)
{
     return (same_cat)?0.0f:1.0f;
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_score_unli_racine_r4(const bool same_cat,const float lambda, const std::int32_t c)
{
     const float inorm{1.0f/(std::fma((static_cast<float>(c)-1.0f),lambda,1.0f))};
     const float term{1.0f-static_cast<float>(c)};
     return (same_cat)?term*inorm*inorm:inorm*std::fma(lambda*term,inorm,1.0f);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_owang_van_ryzin_r4(const float x, const float y, 
                            const float lambda,
                            const float epsilon)
{
     const float term{1.0f-lambda};
     const float ipow_arg{static_cast<float>(static_cast<std::int32_t>(std::fabs(x-y)))};
     return (gms::approximatelyEqual(x,y,epsilon))?term:std::pow(lambda,ipow_arg)*term*0.5f;
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_score_owang_van_ryzin_r4(const float x, const float y, 
                                  const float lambda,
                                  const float epsilon)
{
    const float fabs_term{std::fabs(x-y)};
    const float ipow_arg{static_cast<float>(static_cast<std::int32_t>(std::fabs(x-y)))};
    return (gms::approximatelyEqual(x,y,epsilon))?-1.0f:(0.5*std::pow(lambda,ipow_arg)*(fabs_term/lambda-2.0f));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_oli_racine_r4(const float x, const float y, 
                       const float lambda)           
{
      return (std::pow(lambda,static_cast<float>(static_cast<std::int32_t>(std::fabs(x-y)))));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_score_oli_racine_r4(const float x, const float y, 
                             const float lambda)
{
      const float fabs_term{std::fabs(x-y)};
      return (fabs_term*std::pow(lambda,static_cast<float>(static_cast<std::int32_t>(fabs_term))));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_onli_racine_r4(const float x, const float y, 
                        const float lambda)
{
      const float ratio{(1.0f-lambda)/(1.0f+lambda)};
      const float fabs_term{std::fabs(x-y)};
      return std::pow(lambda,static_cast<float>(static_cast<std::int32_t>(fabs_term)))*ratio;
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_score_onli_racine_r4(const float x, const float y, 
                              const float lambda)
{
      const std::int32_t cxy{static_cast<float>(static_cast<std::int32_t>(std::fabs(x-y)))};
      float cxy_f{static_cast<float>(cxy)};
      return ((cxy!=0) || (lambda!=0.0f))?std::pow(lambda,cxy_f-1.0f)*(cxy_f*(1.0f-lambda*lambda)-
                                                   2.0f*lambda):-2.0f;
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_econvol_rect_r4(const float z)
{
      return ((std::fabs(z)<2.0f)?0.25f:0.0f);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_econvol_gauss2_r4(const float z)
{
      return(0.28209479177387814348f*ceph_expf(-0.25f*z*z));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_econvol_gauss4_r4(const float z)
{
      const float zz{z*z};
      return(0.0044077311214668459918f*ceph_expf(-0.25f*zz)*std::fma(zz,zz-28.0f,108.0f));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_econvol_gauss6_r4(const float z)
{
      const float zz{z*z};
      return(0.00001721769969f*ceph_expf(-0.25f*zz)*std::fma(std::fma(std::fma((-88.0f+zz),zz,2312.0f),zz,-19360.0f),zz,36240.0));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_econvol_gauss8_r4(const float z)
{
      const float zz{z*z};
      return(0.2989183974E-7f*ceph_exp(-0.25f*zz)*std::fma(std::fma(std::fma(std::fma(std::fma((-180.0+zz),zz,11604.0f),
                                                           zz,-331680.0f),zz,4202352.0f),zz,-20462400.0f),zz,25018560.0f));
}

// Derivative kernels
#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_deriv_gauss2_r4(const float z)
{
      constexpr float C0398942280401432677939946059934{0.398942280401432677939946059934f};
      return (-z*C0398942280401432677939946059934*ceph_expf(-0.5f*z*z)); 
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_deriv_gauss4_r4(const float z)
{
     const float zz{z*z};
     constexpr float C0398942280401432677939946059934{0.398942280401432677939946059934f};
     return (-C0398942280401432677939946059934*z*(2.5f-0.5f*zz)*ceph_expf(-0.5f*zz));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_deriv_gauss6_r4(const float z)
{
       const float zz{z*z};
       return (-0.049867785050179084743f*z*ceph_expf(-0.5f*zz)*std::fma(-14.0+zz,zz,35.0f));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_deriv_gauss8_r4(const float z)
{
      const float zz{z*z};
      constexpr float C0398942280401432677939946059934{0.398942280401432677939946059934f};
      return (-C0398942280401432677939946059934*z*std::fma(std::fma((0.5625f-0.02083333333f*zz),zz,-3.9375f),zz,6.5625)*
               ceph_expf(-0.5f*zz));
}

/* **WARNING** The test (z^2 < 5.0f) was removed*/
#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_deriv_epan2_r4(const float z)
{
      return (-0.13416407864998738178f*z);
}

/* **WARNING** The test (z^2 < 5.0f) was removed*/
#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_deriv_epan4_r4(const float z)
{
    const float zz{z*z};
    return z*(2.347871374742824e-1f*zz-8.385254921942804e-1f);
}

/* **WARNING** The test (z^2 < 5.0f) was removed*/
#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_deriv_epan6_r4(const float z)
{
    const float zz{z*z};
    return z*std::fma((1.848948710641142f-2.905490831007508e-1f*zz),zz,-2.567984320334919f);
}

/* **WARNING** The test (z^2 < 5.0f) was removed*/
#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_deriv_epan8_r4(const float z)
{
    const float zz{z*z};
    return z*std::fma(std::fma(std::fma(zz,3.147615066924801e-1f,-2.83285356023232f),zz,7.626913431394709),zz,-5.777964720753567);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_deriv_rect_r4() { return 0.0f;}

// CDF kernels
#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_gauss2_r4(const float z) 
{
    return (std::fma(0.5f,std::erf(0.7071067810f*z),0.5f));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_gauss4_r4(const float z) 
{
     const float expf_term{ceph_expf(-0.5f*z*z)};
     return (std::fma(0.5f,std::erf(0.7071067810*z),std::fma(0.1994711401f*z,expf_term,0.5f)));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_gauss6_r4(const float z)
{
    const float zz{z*z};
    const float expf_term{z*ceph_expf(-0.5f*zz)};
    return (std::fma(0.5f,std::erf(0.7071067810f*z),
                          std::fma(expf_term,(0.3490744952-0.04986778504*zz),0.5f)));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_gauss8_r4(const float z)
{
    const float zz{z*z};
    const float expf_term{z*ceph_expf(-0.5f*zz)};
    return (std::fma(0.5f,std::erf(0.7071067810f*z),
                          std::fma(expf_term,std::fma(std::fma(0.008311297511f*zz,-0.1329807601f),zz,0.4737439578f),0.5f)));
}
#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_epan2_r4(const float z) 
{
    //SQRT(5)
    constexpr float C2236067977499789696409173668731{2.236067977499789696409173668731f};
    return (z < -C2236067977499789696409173668731) ? 0.0f : (z > C2236067977499789696409173668731) ? 1.0f : 
    (std::fma(z,(0.3354101967f-0.02236067978f*z*z),0.5));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_epan4_r4(const float z) 
{
    const float zz{z*z};
    constexpr float C2236067977499789696409173668731{2.236067977499789696409173668731f};
    return (z < -C2236067977499789696409173668731) ? 0.0f : (z > C2236067977499789696409173668731) ? 1.0f : 
    (std::fma(std::fma(std::fma(0.01173935688,zz,-0.1397542486f),zz,0.6288941188f),z,0.5f));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_epan6_r4(const float z) 
{
    const float zz{z*z};
    constexpr float C2236067977499789696409173668731{2.236067977499789696409173668731f};
    return (z < -C2236067977499789696409173668731) ? 0.0f : (z > C2236067977499789696409173668731) ? 1.0f : 
    (std::fma(std::fma(std::fma((0.09244743547f-0.006917835307f*zz),zz,-0.4279973864f),zz,0.9171372566f),z,0.5f));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_epan8_r4(const float z) 
{
    const float zz{z*z};
    constexpr float C2236067977499789696409173668731{2.236067977499789696409173668731f};
    return (z < -C2236067977499789696409173668731) ? 0.0f : (z > C2236067977499789696409173668731) ? 1.0f : 
    (std::fma(std::fma(std::fma(std::fma(std::fma(zz,0.004371687590f,-0.06744889424f),zz,0.3813456714f),zz,-0.9629941194f),zz,1.203742649),z,0.5f));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_rect_r4(const float z)
{
    return (z < -1.0f) ? 0.0f : (z > 1.0f) ? 1.0f : (std::fma(0.5f,z,0.5f));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_owang_van_ryzin_r4(const float y, const float x, 
                                const float lambda, const float epsilon)
{
    if(__builtin_expect(gms::approximatelyEqual(x,y,epsilon),0)) {return (1.0f-0.5f*lambda);}
    const float cxy{ceph_floorf(std::fabs(x-y))};
    const float gee{std::pow(lambda,cxy)};
    return (gms::definitelyLessThan(x,y,epsilon)) ? 0.5f*gee : (1.0f-gee);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_oli_racine_r4(const float y, const float x, 
                           const float lambda, const float ch,
                           const float cl, const float epsilon)
{
    const float xh = gms::definitelyGreaterThan(x,ch,epsilon) ? ch : x;
    const float cxy{ceph_floorf(std::fabs(xh-y))};
    const float one_m_lambda{1.0f-lambda};
    const float gee{std::pow(lambda,cxy)/one_m_lambda};
    if(gms::definitelyLessThan(x,y,epsilon))
    {
        const float true_branch{0.0f};
        const float false_branch{gee*(1.0f-std::pow(lambda,x-cl+1.0f))};
        return (gms::definitelyLessThan(x,cl,epsilon))?true_branch:false_branch;
    }
    else 
    {
        const float one_p_lambda{1.0f+lambda};
        const float lambda_mul_gee{lambda*gee};
        return (one_p_lambda-std::pow(lambda,y-cl+1.0f))/one_m_lambda-lambda_mul_gee;
    }
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_onli_rancine_r4(const float y, const float x, 
                             const float lambda, const float epsilon)
{
    const float cxy{ceph_floorf(std::fabs(x-y))};
    const float gee{std::pow(lambda,cxy)/(1.0f+lambda)};
    return (gms::definitelyLessThan(x,y,epsilon))?gee:1.0f-lambda*gee;
}

// Adaptive convolution kernels.
#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_adaptconvol_gauss2_r4(const float x, const float y,
                               const float hx,const float hy)
{
    const float h2{std::fma(hx,hy,hy*hy)};
    const float z2{(x*y)*(x-y)/h2};
    const float sqrt_h2{std::sqrt(h2)};
    const float ret_term{0.3989422803f*hx*hy};
    return (ret_term*ceph_expf(-0.5f*z2)/sqrt_h2);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float np_adaptconvol_epan2_total_r4(const float x, const float y,
                                    const float hx,const float hy)
{
    constexpr float C6708203932499369089227521006194{6.708203932499369089227521006194f};
    const float ayy{-C6708203932499369089227521006194*y*y};
    const float hl{std::max(hx,hy)};
    const float hs{std::min(hx,hy)};
    const float axy{2.0f*C6708203932499369089227521006194*x*y};
    const float hlhl{hl*hl};
    const float axx{C6708203932499369089227521006194*x*x};   
    const float ahlhl{5.0f*C6708203932499369089227521006194*hlhl};
    const float ahshs{C6708203932499369089227521006194*hs*hs};
    const float left_term{ayy+axy-axx+ahlhl-ahshs};
    return (left_term*hs/(100.0f*hlhl));
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float kernel_ordered_van_ryzin_wang(const float x,const float y,
                                    const float lambda,const float epsilon)
{
    float ret_val{};
    const float one_m_lambda{1.0f-lambda};
    const bool x_eq_y = gms::approximatelyEqual(x,y,epsilon);
    ret_val = (x_eq_y) ? one_m_lambda : std::pow(lambda,ceph_floorf(std::fabs(x-y)))*one_m_lambda*0.5f;
    return (ret_val);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float kernel_ordered_qi(const float x,const float y,
                        const float lambda,const float epsilon)
{
    float ret_val{};
    const bool x_eq_y = gms::approximatelyEqual(x,y,epsilon);
    ret_val = (x_eq_y) ? 1.0f : std::pow(lambda,ceph_floorf(std::fabs(x-y)));
    return (ret_val);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float kernel_unordered_aitchenson_aitken_0(const float x,const float y,const float lambda,
                                           const float c,const float epsilon)
{
    float ret_val{};
    const bool x_eq_y = gms::approximatelyEqual(x,y,epsilon);
    ret_val = (x_eq_y) ? 1.0f-lambda : lambda/(c-1.0f);
    return (ret_val);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif
__ATTR_ALWAYS_INLINE__
static inline 
float kernel_unordered_aitchenson_aitken_1(const float x,const float y,
                                           const float lambda,const float epsilon)
{
    float ret_val{};
    const bool x_eq_y = gms::approximatelyEqual(x,y,epsilon);
    ret_val = (x_eq_y) ? 1.0f : lambda;
    return (ret_val);
}


enum class kernel_functions_types : int32_t 
{
     SECOND_ORDER_GAUSSIAN,
     FOURTH_ORDER_GAUSSIAN,
     SIXTH_ORDER_GAUSSIAN,
     EIGHT_ORDER_GAUSSIAN,
     SECOND_ORDER_EPANECHNIKOV,
     FOURTH_ORDER_EPANECHNIKOV,
     SIXTH_ORDER_EPANECHNIKOV,
     EIGHT_ORDER_EPANECHNIKOV,
     RECTANGULAR_KERNEL
};

enum class kernel_functions_cdf : int32_t 
{
     SECOND_ORDER_GAUSSIAN,
     FOURTH_ORDER_GAUSSIAN,
     SIXTH_ORDER_GAUSSIAN,
     EIGHT_ORDER_GAUSSIAN,
     SECOND_ORDER_EPANECHNIKOV,
     FOURTH_ORDER_EPANECHNIKOV,
     SIXTH_ORDER_EPANECHNIKOV,
     EIGHT_ORDER_EPANECHNIKOV,
     RECTANGULAR_KERNEL
};

enum class kernel_functions_derivative: int32_t 
{
     SECOND_ORDER_GAUSSIAN,
     FOURTH_ORDER_GAUSSIAN,
     SIXTH_ORDER_GAUSSIAN,
     EIGHT_ORDER_GAUSSIAN,
     SECOND_ORDER_EPANECHNIKOV,
     FOURTH_ORDER_EPANECHNIKOV,
     SIXTH_ORDER_EPANECHNIKOV,
     EIGHT_ORDER_EPANECHNIKOV,
     RECTANGULAR_KERNEL
};

enum class kernel_convolution_functions : int32_t 
{
     SECOND_ORDER_GAUSSIAN,
     FOURTH_ORDER_GAUSSIAN,
     SIXTH_ORDER_GAUSSIAN,
     EIGHT_ORDER_GAUSSIAN,
     SECOND_ORDER_EPANECHNIKOV,
     FOURTH_ORDER_EPANECHNIKOV,
     SIXTH_ORDER_EPANECHNIKOV,
     EIGHT_ORDER_EPANECHNIKOV,
     RECTANGULAR_KERNEL
};

enum class kernel_regression_asymptotic_constants : int32_t 
{
     SECOND_ORDER_GAUSSIAN,
     FOURTH_ORDER_GAUSSIAN,
     SIXTH_ORDER_GAUSSIAN,
     EIGHT_ORDER_GAUSSIAN,
     SECOND_ORDER_EPANECHNIKOV,
     FOURTH_ORDER_EPANECHNIKOV,
     SIXTH_ORDER_EPANECHNIKOV,
     EIGHT_ORDER_EPANECHNIKOV,
     RECTANGULAR_KERNEL
};

enum class kernel_density_asymptotic_constants : int32_t 
{
     SECOND_ORDER_GAUSSIAN,
     FOURTH_ORDER_GAUSSIAN,
     SIXTH_ORDER_GAUSSIAN,
     EIGHT_ORDER_GAUSSIAN,
     SECOND_ORDER_EPANECHNIKOV,
     FOURTH_ORDER_EPANECHNIKOV,
     SIXTH_ORDER_EPANECHNIKOV,
     EIGHT_ORDER_EPANECHNIKOV,
     RECTANGULAR_KERNEL
};

enum class kernel_ordered_funcs : int32_t 
{
     WANG_VAN_RYZIN_KERNEL,
     QI_KERNEL
};

enum class kernel_unordered_funcs : int32_t 
{
     AITCHENSON_AITKEN_KERNEL_0,
     AITCHENSON_AITKEN_KERNEL_1
};

template<kernel_ordered_funcs kernel_ordered>
float kernel_ordered_convolution(float * __restrict c_vals,const std::int32_t c_len,
                                 const float x, const float y, const float lambda,
                                 const float epsilon)
{
    
#if defined(__INTEL_COMPILER) || defined(__ICC)
    float * __restrict ptr_c_vals = &c_vals[0];
    assume_aligned(ptr_c_vals,16);   
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
    float * __restrict ptr_c_vals = (float*)__builtin_assume_aligned(c_vals,16);    
#endif  
    float kernel_sum{0.0f}; 
    if constexpr(kernel_ordered==kernel_ordered_funcs::WANG_VAN_RYZIN_KERNEL)
    {
#if (USE_OPENMP) == 1
#pragma simd simdlen(4) linear(i:1) private(curr_val,left_kern_val,right_kern_val) reduction(+:kernel_sum)
#endif 
        for(std::int32_t i = 0;i != c_len; ++i) 
        {
            const float curr_val{ptr_c_vals[i]};
            const float left_kern_val{kernel_ordered_van_ryzin_wang(x,curr_val,lambda,epsilon)};
            const float right_kern_val{kernel_ordered_van_ryzin_wang(y,curr_val,lambda,epsilon)};
            kernel_sum += (left_kern_val*right_kern_val);
        }
        return (kernel_sum);
    }
    else if constexpr(kernel_ordered==kernel_ordered_funcs::QI_KERNEL)
    {
#if (USE_OPENMP) == 1
#pragma simd simdlen(4) linear(i:1) private(curr_val,left_kern_val,right_kern_val) reduction(+:kernel_sum)
#endif 
        for(std::int32_t i = 0;i != c_len; ++i) 
        {
            const float curr_val{ptr_c_vals[i]};
            const float left_kern_val{kernel_ordered_qi(x,curr_val,lambda,epsilon)};
            const float right_kern_val{kernel_ordered_qi(y,curr_val,lambda,epsilon)};
            kernel_sum += (left_kern_val*right_kern_val);
        }
        return (kernel_sum);
    }

}

template<kernel_unordered_funcs kernel_unordered>
float cdf_kernel_unordered(float * __restrict categorical_vals,const std::int32_t vals_len,const float c,
                           const float x,const float y,const float lambda,const float epsilon)
{
#if defined(__INTEL_COMPILER) || defined(__ICC)
    float * __restrict ptr_categorical_vals = &categorical_vals[0];
    assume_aligned(ptr_categorical_vals,16);   
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
    float * __restrict ptr_categorical_vals = (float*)__builtin_assume_aligned(categorical_vals,16);    
#endif  
    float kernel_sum{0.0f}; 
    if constexpr(kernel_unordered==kernel_unordered_funcs::AITCHENSON_AITKEN_KERNEL_0)
    {
#if (USE_OPENMP) == 1
#pragma simd simdlen(4) linear(i:1) private(current_value,current_value_le_x) reduction(+:kernel_sum)
#endif        
        for(std::int32_t i = 0;i != vals_len; ++i) 
        {
            const float current_value{ptr_categorical_vals[i]};
            const bool  current_value_le_x = gms::definitelyLessThan(current_value,x,epsilon);
            if(current_value_le_x)
            {
                kernel_sum += kernel_unordered_aitchenson_aitken_0(current_value,y,lambda,c,epsilon);
            }
        }
        return (kernel_sum);
    }
    else if constexpr(kernel_unordered==kernel_unordered_funcs::AITCHENSON_AITKEN_KERNEL_1)
    {
#if (USE_OPENMP) == 1
#pragma simd simdlen(4) linear(i:1) private(current_value,current_value_le_x) reduction(+:kernel_sum)
#endif        
        for(std::int32_t i = 0;i != vals_len; ++i) 
        {
            const float current_value{ptr_categorical_vals[i]};
            const bool  current_value_le_x = gms::definitelyLessThan(current_value,x,epsilon);
            if(current_value_le_x)
            {
                kernel_sum += kernel_unordered_aitchenson_aitken_1(current_value,y,lambda,epsilon);
            }
        }
        return (kernel_sum);
    }
}

template<kernel_unordered_funcs kernel_unordered>
float kernel_unordered_convolution(float * __restrict c_vals,const std::int32_t c_vals_len,
                                   const float x,const float y,const float lambda,const float c,
                                   const float epsilon)
{
#if defined(__INTEL_COMPILER) || defined(__ICC)
    float * __restrict ptr_c_vals = &c_vals[0];
    assume_aligned(ptr_c_vals,16);   
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
    float * __restrict ptr_c_vals = (float*)__builtin_assume_aligned(c_vals,16);    
#endif  
    float kernel_sum{0.0f}; 
    if constexpr(kernel_unordered==kernel_unordered_funcs::AITCHENSON_AITKEN_KERNEL_0)
    {
#if (USE_OPENMP) == 1
#pragma simd simdlen(4) linear(i:1) private(current_value,left_kern_val,right_kern_val) reduction(+:kernel_sum)
#endif          
        for(std::int32_t i = 0;i != c_vals_len; ++i) 
        {
            const float current_value{ptr_c_vals[i]};
            const float left_kern_val{kernel_unordered_aitchenson_aitken_0(x,current_value,lambda,c,epsilon)};
            const float right_kern_val{kernel_unordered_aitchenson_aitken_0(y,current_value,lambda,c,epsilon)};
            kernel_sum += (left_kern_val*right_kern_val);
        }
        return (kernel_sum);
    }
    else if constexpr(kernel_unordered==kernel_unordered_funcs::AITCHENSON_AITKEN_KERNEL_1)
    {
#if (USE_OPENMP) == 1
#pragma simd simdlen(4) linear(i:1) private(current_value,left_kern_val,right_kern_val) reduction(+:kernel_sum)
#endif          
        for(std::int32_t i = 0;i != c_vals_len; ++i) 
        {
            const float current_value{ptr_c_vals[i]};
            const float left_kern_val{kernel_unordered_aitchenson_aitken_1(x,current_value,lambda,epsilon)};
            const float right_kern_val{kernel_unordered_aitchenson_aitken_1(y,current_value,lambda,epsilon)};
            kernel_sum += (left_kern_val*right_kern_val);
        }
        return (kernel_sum);
    }
}


template<kernel_regression_asymptotic_constants regress_constants>
void init_kernel_regression_asymptotic_constants(const float num_reg_continous,
                                                 float & INT_KERNEL_P,
                                                 float & K_INT_KERNEL_P,
                                                 float & INT_KERNEL_PM_HALF,
                                                 float & DIFF_KER_PPM)
{
    if constexpr(regress_constants==kernel_regression_asymptotic_constants::SECOND_ORDER_GAUSSIAN)
    {
        INT_KERNEL_P = 0.28209479177387814348f;
		K_INT_KERNEL_P = std::pow(INT_KERNEL_P, num_reg_continuous);
		INT_KERNEL_PM_HALF = 0.21969564473386119853f;
		DIFF_KER_PPM = (2.0f*(K_INT_KERNEL_P/INT_KERNEL_P)*0.062399147040017f);
        return;
    }
    else if constexpr(regress_constants==kernel_regression_asymptotic_constants::FOURTH_ORDER_GAUSSIAN)
    {
        INT_KERNEL_P = 0.47603496111841936711f;
		K_INT_KERNEL_P = std::pow(INT_KERNEL_P,num_reg_continuous);
		INT_KERNEL_PM_HALF = 0.27805230036629307938f;
		DIFF_KER_PPM = (2.0f*(K_INT_KERNEL_P/INT_KERNEL_P)*0.19798266075212628773f);
        return;
    }
    else if constexpr(regress_constants==kernel_regression_asymptotic_constants::SIXTH_ORDER_GAUSSIAN)
    {
        INT_KERNEL_P = 0.62396943688265038571f;
		K_INT_KERNEL_P = std::pow(INT_KERNEL_P, num_reg_continuous);
		INT_KERNEL_PM_HALF = 0.25618196366213489976f;
		DIFF_KER_PPM = (2.0f*(K_INT_KERNEL_P/INT_KERNEL_P)*0.36778747322051548595f);
        return;
    }
    else if constexpr(regress_constants==kernel_regression_asymptotic_constants::EIGHT_ORDER_GAUSSIAN)
    {
        INT_KERNEL_P = 0.74785078617543927990f;
		K_INT_KERNEL_P = std::pow(INT_KERNEL_P, num_reg_continuous);
		INT_KERNEL_PM_HALF = 0.19644083574560137818f;
		DIFF_KER_PPM = (2.0f*(K_INT_KERNEL_P/INT_KERNEL_P)*0.55140995042983790172f);
        return;
    }
    else if constexpr(regress_constants==kernel_regression_asymptotic_constants::SECOND_ORDER_EPANECHNIKOV)
    {
        INT_KERNEL_P = 0.26832815729997476357f;
		K_INT_KERNEL_P = std::pow(INT_KERNEL_P, num_reg_continuous);
		INT_KERNEL_PM_HALF = 0.20250390621232470438f;
		DIFF_KER_PPM = (2.0f*(K_INT_KERNEL_P/INT_KERNEL_P)*0.0658242510876501f);
        return;
    }
    else if constexpr(regress_constants==kernel_regression_asymptotic_constants::FOURTH_ORDER_EPANECHNIKOV)
    {
        INT_KERNEL_P = 0.55901699437494742410f;
		K_INT_KERNEL_P = std::pow(INT_KERNEL_P, num_reg_continuous);
		INT_KERNEL_PM_HALF = 0.25635637709255874475f;
		DIFF_KER_PPM = (2.0f*(K_INT_KERNEL_P/INT_KERNEL_P)*0.30266061728238867935f);
        return;
    }
    else if constexpr(regress_constants==kernel_regression_asymptotic_constants::SIXTH_ORDER_EPANECHNIKOV)
    {
        INT_KERNEL_P = 0.84658823667359826246f;
		K_INT_KERNEL_P = std::pow(INT_KERNEL_P, num_reg_continuous);
		INT_KERNEL_PM_HALF = 0.27428761935713012265f;
		DIFF_KER_PPM = (2.0f*(K_INT_KERNEL_P/INT_KERNEL_P)*0.57230061731646813981f);
        return;
    }
    else if constexpr(regress_constants==kernel_regression_asymptotic_constants::EIGHT_ORDER_EPANECHNIKOV)
    {
        INT_KERNEL_P = 1.1329342579014329689f;
		K_INT_KERNEL_P = std::pow(INT_KERNEL_P, num_reg_continuous);
		INT_KERNEL_PM_HALF = 0.15585854498586945817f;
		DIFF_KER_PPM = (2.0f*(K_INT_KERNEL_P/INT_KERNEL_P)*0.97707571291556351073f);
        return;
    }
    else if constexpr(regress_constants==kernel_regression_asymptotic_constants::RECTANGULAR_KERNEL)
    {

			INT_KERNEL_P = 0.5f;
			K_INT_KERNEL_P = std::pow(INT_KERNEL_P, num_reg_continuous);
			INT_KERNEL_PM_HALF = 0.25f;
			DIFF_KER_PPM = (2.0f*(K_INT_KERNEL_P/INT_KERNEL_P)*0.25f);
            return;
    }
}

template<kernel_density_asymptotic_constants density_constants>
void init_kernel_density_asymptotic_constants(const float num_var_continous,
                                              float & INT_KERNEL_P,
                                              float & K_INT_KERNEL_P)
{
    /* Initialize constants for various kernels required for asymptotic standard errors */
    if constexpr(density_constants==kernel_density_asymptotic_constants::SECOND_ORDER_GAUSSIAN)
    {
         INT_KERNEL_P = 0.28209479177387814348f;
		 K_INT_KERNEL_P = std::pow(INT_KERNEL_P, num_var_continuous);
         return;
    }
    else if constexpr(density_constants==kernel_density_asymptotic_constants::FOURTH_ORDER_GAUSSIAN)
    {
         INT_KERNEL_P = 0.47603496111841936711f;
		 K_INT_KERNEL_P = std::pow(INT_KERNEL_P, num_var_continuous);
         return; 
    }
    else if constexpr(density_constants==kernel_density_asymptotic_constants::SIXTH_ORDER_GAUSSIAN)
    {
         INT_KERNEL_P = 0.62396943688265038571f;
		 K_INT_KERNEL_P = std::pow(INT_KERNEL_P, num_var_continuous);
         return;
    }
    else if constexpr(density_constants==kernel_density_asymptotic_constants::EIGHT_ORDER_GAUSSIAN)
    {
         INT_KERNEL_P = 0.74785078617543927990f;
		 K_INT_KERNEL_P = std::pow(INT_KERNEL_P, num_var_continuous);
         return;
    }
    else if constexpr(density_constants==kernel_density_asymptotic_constants::SECOND_ORDER_EPANECHNIKOV)
    {
         INT_KERNEL_P = 0.26832815729997476357f;
		 K_INT_KERNEL_P = std::pow(INT_KERNEL_P, num_var_continuous);
         return;
    }
    else if constexpr(density_constants==kernel_density_asymptotic_constants::FOURTH_ORDER_EPANECHNIKOV)
    {
         INT_KERNEL_P = 0.55901699437494742410f;
		 K_INT_KERNEL_P = std::pow(INT_KERNEL_P, num_var_continuous);
         return;
    }
    else if constexpr(density_constants==kernel_density_asymptotic_constants::SIXTH_ORDER_EPANECHNIKOV)
    {
         INT_KERNEL_P = 0.84658823667359826246f;
		 K_INT_KERNEL_P = std::pow(INT_KERNEL_P, num_var_continuous);
         return;
    }
    else if constexpr(density_constants==kernel_density_asymptotic_constants::EIGHT_ORDER_EPANECHNIKOV)
    {
         INT_KERNEL_P = 1.1329342579014329689f;
		 K_INT_KERNEL_P = std::pow(INT_KERNEL_P, num_var_continuous);
         return;
    }
    else if constexpr(density_constants==kernel_density_asymptotic_constants::RECTANGULAR_KERNEL)
    {
         INT_KERNEL_P = 0.5f;
		 K_INT_KERNEL_P = std::pow(INT_KERNEL_P, num_var_continuous);
         return;
    }
}




} //np_standalone_funcs

} // math

} // gms

#endif /*__GMS_PDF_KERNELS_SCALAR_H__*/