
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

/*
   For Epanechnikov kernels the input check i.e. z^2 < 5.0 was removed 
   in order to remove the branch-misprediction.
*/

/*
    The input range is (z*z < 5.0)?result:0.0
*/
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
__ATTR_ALWAYS_INLINE__
static inline 
float np_epanechnikov8_r4(const float z)
{
    const float zz{z*z};
    const float fma_arg{4.1056640625f-0.5865234375f*zz};
    const float fma_term{std::fma(std::fma(fma_arg,zz,-7.8955078125f),zz,3.5888671875f)};
    return (static_cast<float>(0.33541019662496845446f*fma_term*(1.0f-0.2f*zz)));
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_rect_r4(const float z)
{
    return (z*z < 1.0f)?0.5f:0.0f;
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_uaa_r4(const bool same_cat,const float lambda, const std::int32_t c)
{
    return (same_cat)?(1.0f-lambda):lambda/((static_cast<float>(c)-1.0f));
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_score_uaa_r4(const bool same_cat,const float lambda, const std::int32_t c)
{
    return (same_cat)?-1.0f:(1.0f/(static_cast<float>(c)-1.0f)); 
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_uli_racine_r4(const bool same_cat,const float lambda)
{
    return (same_cat)?1.0f:lambda;
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_unli_racine_r4(const bool same_cat,const float lambda, const std::int32_t c)
{
     return ((same_cat)?1.0f:lambda)/((static_cast<float>(c)-1.0f)*lambda + 1.0f);
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_score_uli_racine_r4(const bool same_cat,const float lambda)
{
     return (same_cat)?0.0f:1.0f;
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_score_unli_racine_r4(const bool same_cat,const float lambda, const std::int32_t c)
{
     const float inorm{1.0f/(std::fma((static_cast<float>(c)-1.0f),lambda,1.0f))};
     const float term{1.0f-static_cast<float>(c)};
     return (same_cat)?term*inorm*inorm:inorm*std::fma(lambda*term,inorm,1.0f);
}

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

__ATTR_ALWAYS_INLINE__
static inline 
float np_oli_racine_r4(const float x, const float y, 
                       const float lambda)           
{
      return (std::pow(lambda,static_cast<float>(static_cast<std::int32_t>(std::fabs(x-y)))));
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_score_oli_racine_r4(const float x, const float y, 
                             const float lambda)
{
      const float fabs_term{std::fabs(x-y)};
      return (fabs_term*std::pow(lambda,static_cast<float>(static_cast<std::int32_t>(fabs_term))));
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_onli_racine_r4(const float x, const float y, 
                        const float lambda)
{
      const float ratio{(1.0f-lambda)/(1.0f+lambda)};
      const float fabs_term{std::fabs(x-y)};
      return std::pow(lambda,static_cast<float>(static_cast<std::int32_t>(fabs_term)))*ratio;
}

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

__ATTR_ALWAYS_INLINE__
static inline 
float np_econvol_rect_r4(const float z)
{
      return ((std::fabs(z)<2.0f)?0.25f:0.0f);
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_econvol_gauss2_r4(const float z)
{
      return(0.28209479177387814348f*ceph_expf(-0.25f*z*z));
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_econvol_gauss4_r4(const float z)
{
      const float zz{z*z};
      return(0.0044077311214668459918f*ceph_expf(-0.25f*zz)*std::fma(zz,zz-28.0f,108.0f));
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_econvol_gauss6_r4(const float z)
{
      const float zz{z*z};
      return(0.00001721769969f*ceph_expf(-0.25f*zz)*std::fma(std::fma(std::fma((-88.0f+zz),zz,2312.0f),zz,-19360.0f),zz,36240.0));
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_econvol_gauss8_r4(const float z)
{
      const float zz{z*z};
      return(0.2989183974E-7f*ceph_exp(-0.25f*zz)*std::fma(std::fma(std::fma(std::fma(std::fma((-180.0+zz),zz,11604.0f),
                                                           zz,-331680.0f),zz,4202352.0f),zz,-20462400.0f),zz,25018560.0f));
}

// Derivative kernels
__ATTR_ALWAYS_INLINE__
static inline 
float np_deriv_gauss2_r4(const float z)
{
      constexpr float C0398942280401432677939946059934{0.398942280401432677939946059934f};
      return (-z*C0398942280401432677939946059934*ceph_expf(-0.5f*z*z)); 
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_deriv_gauss4_r4(const float z)
{
     const float zz{z*z};
     constexpr float C0398942280401432677939946059934{0.398942280401432677939946059934f};
     return (-C0398942280401432677939946059934*z*(2.5f-0.5f*zz)*ceph_expf(-0.5f*zz));
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_deriv_gauss6_r4(const float z)
{
       const float zz{z*z};
       return (-0.049867785050179084743f*z*ceph_expf(-0.5f*zz)*std::fma(-14.0+zz,zz,35.0f));
}

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
__ATTR_ALWAYS_INLINE__
static inline 
float np_deriv_epan2_r4(const float z)
{
      return (-0.13416407864998738178f*z);
}

/* **WARNING** The test (z^2 < 5.0f) was removed*/
__ATTR_ALWAYS_INLINE__
static inline 
float np_deriv_epan4_r4(const float z)
{
    const float zz{z*z};
    return z*(2.347871374742824e-1f*zz-8.385254921942804e-1f);
}

/* **WARNING** The test (z^2 < 5.0f) was removed*/
__ATTR_ALWAYS_INLINE__
static inline 
float np_deriv_epan6_r4(const float z)
{
    const float zz{z*z};
    return z*std::fma((1.848948710641142f-2.905490831007508e-1f*zz),zz,-2.567984320334919f);
}

/* **WARNING** The test (z^2 < 5.0f) was removed*/
__ATTR_ALWAYS_INLINE__
static inline 
float np_deriv_epan8_r4(const float z)
{
    const float zz{z*z};
    return z*std::fma(std::fma(std::fma(zz,3.147615066924801e-1f,-2.83285356023232f),zz,7.626913431394709),zz,-5.777964720753567);
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_deriv_rect_r4() { return 0.0f;}

// CDF kernels

__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_gauss2_r4(const float z) 
{
    return (std::fma(0.5f,std::erf(0.7071067810f*z),0.5f));
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_gauss4_r4(const float z) 
{
     const float expf_term{ceph_expf(-0.5f*z*z)};
     return (std::fma(0.5f,std::erf(0.7071067810*z),std::fma(0.1994711401f*z,expf_term,0.5f)));
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_gauss6_r4(const float z)
{
    const float zz{z*z};
    const float expf_term{z*ceph_expf(-0.5f*zz)};
    return (std::fma(0.5f,std::erf(0.7071067810f*z),
                          std::fma(expf_term,(0.3490744952-0.04986778504*zz),0.5f)));
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_gauss8_r4(const float z)
{
    const float zz{z*z};
    const float expf_term{z*ceph_expf(-0.5f*zz)};
    return (std::fma(0.5f,std::erf(0.7071067810f*z),
                          std::fma(expf_term,std::fma(std::fma(0.008311297511f*zz,-0.1329807601f),zz,0.4737439578f),0.5f)));
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_epan2_r4(const float z) 
{
    //SQRT(5)
    constexpr float C2236067977499789696409173668731{2.236067977499789696409173668731f};
    return (z < -C2236067977499789696409173668731) ? 0.0f : (z > C2236067977499789696409173668731) ? 1.0f : 
    (std::fma(z,(0.3354101967f-0.02236067978f*z*z),0.5));
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_epan4_r4(const float z) 
{
    const float zz{z*z};
    constexpr float C2236067977499789696409173668731{2.236067977499789696409173668731f};
    return (z < -C2236067977499789696409173668731) ? 0.0f : (z > C2236067977499789696409173668731) ? 1.0f : 
    (std::fma(std::fma(std::fma(0.01173935688,zz,-0.1397542486f),zz,0.6288941188f),z,0.5f));
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_epan6_r4(const float z) 
{
    const float zz{z*z};
    constexpr float C2236067977499789696409173668731{2.236067977499789696409173668731f};
    return (z < -C2236067977499789696409173668731) ? 0.0f : (z > C2236067977499789696409173668731) ? 1.0f : 
    (std::fma(std::fma(std::fma((0.09244743547f-0.006917835307f*zz),zz,-0.4279973864f),zz,0.9171372566f),z,0.5f));
}

__ATTR_ALWAYS_INLINE__
static inline 
float np_cdf_epan8_r4(const float z) 
{
    const float zz{z*z};
    constexpr float C2236067977499789696409173668731{2.236067977499789696409173668731f};
    return (z < -C2236067977499789696409173668731) ? 0.0f : (z > C2236067977499789696409173668731) ? 1.0f : 
    (std::fma(std::fma(std::fma(std::fma(std::fma(zz,0.004371687590f,-0.06744889424f),zz,0.3813456714f),zz,-0.9629941194f),zz,1.203742649),z,0.5f));
}

} //np_standalone_funcs

} // math

} // gms















#endif /*__GMS_PDF_KERNELS_SCALAR_H__*/