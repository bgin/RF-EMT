
/*
   Copyright 2021 Bruce Dawson

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
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

#ifndef __GMS_COMPARE_FP_SAFE_H__
#define __GMS_COMPARE_FP_SAFE_H__ 050920220913





namespace file_version 
{

    const unsigned int GMS_COMPARE_FP_SAFE_MAJOR = 1U;
    const unsigned int GMS_COMPARE_FP_SAFE_MINOR = 0U;
    const unsigned int GMS_COMPARE_FP_SAFE_MICRO = 0U;
    const unsigned int GMS_COMPARE_FP_SAFE_FULLVER =
      1000U*GMS_COMPARE_FP_SAFE_MAJOR+
      100U*GMS_COMPARE_FP_SAFE_MINOR+
      10U*GMS_COMPARE_FP_SAFE_MICRO;
    const char * const GMS_COMPARE_FP_SAFE_CREATION_DATE = "05-09-2022 09:13 AM +00200 (MON 05 SEP 2022 GMT+2)";
    const char * const GMS_COMPARE_FP_SAFE_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_COMPARE_FP_SAFE_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_COMPARE_FP_SAFE_DESCRIPTION   = "Based on book: Mathematical Theory of Electro-Optical Sensors (rus).";

}


#include <cmath>
#include <assert.h>
#include <float.h>
#include <cstdint>
#include "GMS_config.h"
#if (USE_OPENMP) == 1
#include <omp.h>
#endif

namespace gms
{

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif 
__ATTR_ALWAYS_INLINE__
static inline
bool approximatelyEqual(const float a,
		                const float b,
					    const float epsilon) 
{
	const float fabsa = std::fabs(a);
	const float fabsb = std::fabs(b);
    return std::fabs(a - b) <=
        ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
}


#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(8)
#endif 
__ATTR_ALWAYS_INLINE__
static inline
bool approximatelyEqual(const double a,
		                const double b,
					    const double epsilon) 
{
	const double fabsa = std::fabs(a);
	const double fabsb = std::fabs(b);
    return fabs(a - b) <=
        ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif 
__ATTR_ALWAYS_INLINE__
static inline
bool essentiallyEqual(const float a,
		              const float b,
					  const float epsilon) 
{
    const float fabsa = std::fabs(a);
	const float fabsb = std::fabs(b);
    return fabs(a - b) <=
		((fabsa > fabsb ? fabsb : fabsa) * epsilon);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(8)
#endif 
__ATTR_ALWAYS_INLINE__
static inline
bool essentiallyEqual(const double a,
		              const double b,
					  const double epsilon) 
{
    const double fabsa = std::fabs(a);
	const double fabsb = std::fabs(b);
    return fabs(a - b) <=
		((fabsa > fabsb ? fabsb : fabsa) * epsilon);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif 
__ATTR_ALWAYS_INLINE__
static inline		   
bool definitelyGreaterThan(const float a,
		                   const float b,
					       const float epsilon) 
{
    const float fabsa = std::fabs(a);
	const float fabsb = std::fabs(b);
    return (a - b) >
		((fabsa < fabsb ? fabsb : fabsa) * epsilon);
}


#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(8)
#endif 
__ATTR_ALWAYS_INLINE__
static inline
bool definitelyGreaterThan(const double a,
		                   const double b,
					       const double epsilon) 
{
    const double fabsa = std::fabs(a);
	const double fabsb = std::fabs(b);
    return (a - b) >
		((fabsa < fabsb ? fabsb : fabsa) * epsilon);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif 
__ATTR_ALWAYS_INLINE__
static inline
bool definitelyLessThan(const float a,
		                const float b,
					    const float epsilon) 
{
    const float fabsa = std::fabs(a);
	const float fabsb = std::fabs(b);
    return (b - a) >
		((fabsa < fabsb ? fabsb : fabsa) * epsilon);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(8)
#endif 
__ATTR_ALWAYS_INLINE__
static inline
bool definitelyLessThan(const double a,
		                const double b,
					    const double epsilon) 
{
    const double fabsa = std::fabs(a);
	const double fabsb = std::fabs(b);
    return (b - a) >
		((fabsa < fabsb ? fabsb : fabsa) * epsilon);
}

/*
    Adapted from "CompareAsInt.cpp" source file.
    Generally shall be used for any floating-point comparison!!
*/
/*
   Disassembly of this function, GCC-16.1 -O3 

 almostEqual2sComplement(float, float, int):
 mov    ecx,0x80000000
 movd   eax,xmm0
 mov    edx,ecx
 sub    edx,eax
 test   eax,eax
 cmovs  eax,edx
 movd   edx,xmm1
 sub    ecx,edx
 test   edx,edx
 cmovs  edx,ecx
 sub    eax,edx
 mov    edx,eax
 neg    edx
 cmovns eax,edx
 cmp    eax,edi
 setle  al
 ret
*/
// Usable AlmostEqual function
#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(8)
#endif 
__ATTR_ALWAYS_INLINE__
static inline
bool 
almostEqual2sComplement(float A,float B, std::int32_t maxUlps) 
{
  // Make sure maxUlps is non-negative and small enough that the
  // default NAN won't compare as equal to anything.
  // This check disabled for now so that my tests run properly.
  // This assert should be enabled for normal use.
  // assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
  int aInt =  *reinterpret_cast<std::int32_t*>(&A);   //*(int *)&A;
  // Make aInt lexicographically ordered as a twos-complement int
  if (aInt < 0)
    aInt = 0x80000000 - aInt;
  // Make bInt lexicographically ordered as a twos-complement int
  int bInt =   *reinterpret_cast<std::int32_t*>(&B);  //*(int *)&B;
  if (bInt < 0)
    bInt = 0x80000000 - bInt;
  int intDiff = std::abs(aInt - bInt);
  if (intDiff <= maxUlps)
    return true;
  return false;
}

// More comprehence version for the final usage in e.g. production code

#if !defined(COMPARE_FP_SAFE_INFINITYCHECK)
#define COMPARE_FP_SAFE_INFINITYCHECK
#endif 
#if !defined(COMPARE_FP_SAFE_NANCHECK)
#define COMPARE_FP_SAFE_NANCHECK
#endif 
#if !defined(COMPARE_FP_SAFE_SIGNCHECK)
#define COMPARE_FP_SAFE_SIGNCHECK
#endif 

__ATTR_ALWAYS_INLINE__
static inline
bool isInfinite(float A) 
{
  constexpr std::int32_t kInfAsInt = 0x7F800000;

  // An infinity has an exponent of 255 (shift left 23 positions) and
  // a zero mantissa. There are two infinities - positive and negative.
  //if ((*(int *)&A & 0x7FFFFFFF) == kInfAsInt)
    if((*reinterpret_cast<std::int32_t*>(&A) & 0x7FFFFFFF) == kInfAsInt)
    return true;
  return false;
}

__ATTR_ALWAYS_INLINE__
static inline
bool isNan(float A) 
{
  // A NAN has an exponent of 255 (shifted left 23 positions) and
  // a non-zero mantissa.
  //int exp = *(int *)&A & 0x7F800000;
  //int mantissa = *(int *)&A & 0x007FFFFF;
  std::int32_t exp      = *reinterpret_cast<std::int32_t*>(&A) & 0x7F800000;
  std::int32_t mantissa = *reinterpret_cast<std::int32_t*>(&A) & 0x007FFFFF; 
  if (exp == 0x7F800000 && mantissa != 0)
    return true;
  return false;
}

__ATTR_ALWAYS_INLINE__
static inline
std::int32_t sign_bit(float A) {
  // The sign bit of a number is the high bit.
  //return (*(int *)&A) & 0x80000000;
  return (*reinterpret_cast<std::int32_t*>(&A)) & 0x80000000;
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(8)
#endif 
__ATTR_ALWAYS_INLINE__
static inline
bool 
almostEqualUlpsFinal(float A, float B, 
                     const std::int32_t maxUlps) 
{
  
#ifdef COMPARE_FP_SAFE_INFINITYCHECK
  // If A or B are infinity (positive or negative) then
  // only return true if they are exactly equal to each other -
  // that is, if they are both infinities of the same sign.
  // This check is only needed if you will be generating
  // infinities and you don't want them 'close' to numbers
  // near FLT_MAX.
  if (isInfinite(A) || isInfinite(B))
    return A == B;
#endif

#ifdef COMPARE_FP_SAFE_NANCHECK
  // If A or B are a NAN, return false. NANs are equal to nothing,
  // not even themselves.
  // This check is only needed if you will be generating NANs
  // and you use a maxUlps greater than 4 million or you want to
  // ensure that a NAN does not equal itself.
  if (isNan(A) || isNan(B))
    return false;
#endif

#ifdef COMPARE_FP_SAFE_SIGNCHECK
  // After adjusting floats so their representations are lexicographically
  // ordered as twos-complement integers a very small positive number
  // will compare as 'close' to a very small negative number. If this is
  // not desireable, and if you are on a platform that supports
  // subnormals (which is the only place the problem can show up) then
  // you need this check.
  // The check for A == B is because zero and negative zero have different
  // signs but are equal to each other.
  if (sign_bit(A) != sign_bit(B))
    return A == B;
#endif

  //int aInt = *(int *)&A;
  std::int32_t aInt = *reinterpret_cast<std::int32_t*>(&A);
  // Make aInt lexicographically ordered as a twos-complement int
  if (aInt < 0)
    aInt = 0x80000000 - aInt;
  // Make bInt lexicographically ordered as a twos-complement int
  //int bInt = *(int *)&B;
  std::int32_t bInt = *reinterpret_cast<std::int32_t*>(&B);
  if (bInt < 0)
    bInt = 0x80000000 - bInt;

  // Now we can compare aInt and bInt to find out how far apart A and B
  // are.
  int intDiff = std::abs(aInt - bInt);
  if (intDiff <= maxUlps)
    return true;
  return false;
}

/*
    The 'offset' argument may be set to 0.
*/
__ATTR_ALWAYS_INLINE__
static inline 
void 
printNumber(float f, const std::int32_t offset) 
{
  (*((int *)&f)) += offset;
  std::printf("%+1.11g,0x%08X,%d\n", f, *(int *)&f, *(int *)&f);
}

__ATTR_ALWAYS_INLINE__
static inline 
void 
test_1_almostEqual2sComplement()
{
    constexpr const float x{0.3989423f};
    constexpr const float y{0.3989429f};
    bool is_equal{almostEqual2sComplement(x,y,10)};
    std::cout << std::boolalpha << is_equal <<"\n";
    printNumber(x,0);
    printNumber(y,0);
}

__ATTR_ALWAYS_INLINE__
static inline 
void 
test_2_almostEqual2sComplement(float x, float y, 
                                    const std::int32_t maxUlps)
{
    bool is_equal{almostEqual2sComplement(x,y,10)};
    std::cout << std::boolalpha << is_equal <<"\n";
    printNumber(x,0);
    printNumber(y,0);
}



} // gms


#endif /*__GMS_COMPARE_FP_SAFE_H__*/