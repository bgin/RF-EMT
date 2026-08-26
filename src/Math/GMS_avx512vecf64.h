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
#ifndef __GMS_AVX512VECF64_H__
#define __GMS_AVX512VECF64_H__ 131020191146

namespace file_info {


  	const unsigned int gGMS_AVX512VECF64_MAJOR = 1;

	const unsigned int gGMS_AVX512VECF64_MINOR = 0;

	const unsigned int gGMS_AVX512VECF64_MICRO = 1;

	const unsigned int gGMS_AVX512VECF64_FULLVER = 
	    1000U*gGMS_AVX512VECF64_MAJOR + 100U*gGMS_AVX512VECF64_MINOR + 10U*gGMS_AVX512VECF64_MICRO;

	const  char * const  pgGMS_AVX512VEC64_CREATE_DATE = "13-10-2019 11:46 +00200 (SUN 13 OCT 2019 GMT+2)";

	const  char * const  pgGMS_AVX512VECG64_BUILD_DATE = __DATE__":"__TIME__;

	const  char * const  pgGMS_AVX512VECF64_AUTHOR     = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const  char * const  pgGMS_AVX512VECF64_DESCRIPT   =  "Wrapper class around __m512d data type";
}



#include <iostream>
#include <cstdint>
#include <zmmintrin.h>

namespace gms {
	namespace math {


		__attribute__((align(64))) struct AVX512Vec8 {

		

		   //
		   // Class Constructors and Destructor
		   //
		   AVX512Vec8();

		   AVX4512Vec8(const double[8]); // Unaligned load

		   AVX512Vec8(const double);

		   AVX512Vec8(const double,
				const double,
				const double,
				const double,
				const double,
				const double,
				const double,
				const double);

		   AVX512Vec8(const __m512d);

		   AVX512Vec8(const __m512i);

		   AVX512Vec8(const AVX512Vec8 &);

		   AVX512Vec8(const __m256d,
			        const __m256d);

		   ~AVX512Vec8() = default;

		   //
		   // getters
		   //
		   const __m512d get_v8() const;

		   __m512d get_v8();

		   __m256d lo_part() const;

		   __m256d hi_part() const;

		   //
		   //	Load/store functions
		   //

		   // Address argument should be aligned on 64-byte boundary
		   AVX512Vec8 & load_a(const double* __restrict);

		   AVX512Vec8 & load_u(const double* __restrict);

		   // Address argument should be aligned on 64-byte boundary
		   void store_a(double* __restrict) const;

		   void store_u(double* __restrict) const;

		   void stream_store(double* __restrict) const;

		   double extract_scalar(const uint32_t) const;

		   //
		   // Few member and friend operators
		   //

		   AVX512Vec8 & operator=(const AVX512Vec8 &);

		   // Type-cast operator
		   operator __m512d () const;

		   double operator[](const uint32_t) const;

		   friend std::ostream & operator<<(std::ostream &,
						    const AVX512Vec8 &);


		  

	                 __m512d m_v8;

		};

		// Global (namespace) static functions

		// Extract __m256d part only
		static inline __m256d extract(const AVX512Vec8 ,
					      const int32_t);

		// Branchless conditional selection
		static inline AVX512Vec8 select(const AVX512Vec8 ,
						  const AVX512Vec8 ,
						  const __mmask8);
		//
		//	Arithmetic and mathematical operations
		//

		// SIMD max
		static inline AVX512Vec8 simd_max(const AVX512Vec8 ,
						    const AVX512Vec8 );

		// SIMD min
		static inline AVX512Vec8 simd_min(const AVX512Vec8 ,
						    const AVX512Vec8 );

	    // SIMD abs
		static inline AVX512Vec8 abs(const AVX512Vec8 );

	    // SIMD sqrt
		static inline AVX512Vec8 sqrt(const AVX512Vec8 ); 

	    // SIMD rsqrt
		static inline AVX512Vec8 rsqrt(const AVX512Vec8);

	    // SIMD cbrt
		static inline AVX512Vec8 cbrt(const AVX512Vec8); 

	    // SIMD squared
		static inline AVX512Vec8 sqr(const AVX512Vec8);

		// Horizontal reduction by addition
		static inline double reduce_add(const AVX512Vec8);
		
		// Horizontal reduction by multiplication
		static inline double reduce_mul(const AVX512Vec8);

		// Horizontal reduction by maximum
		static inline double reduce_max(const AVX512Vec8);

		// Horizontal reduction by minimmum
		static inline double reduce_min(const AVX512Vec8);

	    // SIMD ceil
		static inline AVX512Vec8 ceil(const AVX512Vec8);

	    // SIMD floor
		static inline AVX512Vec8 floor(const AVX512Vec8 );

	    // SIMD round
		static inline AVX512Vec8 round(const AVX512Vec8 ,
						 const int32_t);

	    // SVML sin
		static inline AVX512Vec8 sin(const AVX512Vec8 );

		// SVML sind (degree)
		static inline AVX512Vec8 sind(const AVX512Vec8 );

	    // SVML cos
		static inline AVX512Vec8 cos(const AVX512Vec8 );

		// SVML cosd (degree)
		static inline AVX512Vec8 cosd(const AVX512Vec8 );

	    // SVML sinh
		static inline AVX512Vec8 sinh(const AVX512Vec8 );

	    //  SVML cosh
		static inline AVX512Vec8 cosh(const AVX512Vec8 );

        // SVML tan
		static inline AVX512Vec8 tan(const AVX512Vec8 );

	    // SVML tanh
		static inline AVX512Vec8 tanh(const AVX512Vec8 );

		// SVML asin
		static inline AVX512Vec8 asin(const AVX512Vec8 );

		// SVML asinh
		static inline AVX512Vec8 asinh(const AVX512Vec8 );

		// SVML acos
		static inline AVX512Vec8 acos(const AVX512Vec8 );

		// SVML acosh
		static inline AVX512Vec8 acosh(const AVX512Vec8 );

		// SVML atan
		static inline AVX512Vec8 atan(const AVX512Vec8 );

		// SVML atanh
		static inline AVX512Vec8 atanh(const AVX512Vec8 );

		// SVML log
		static inline AVX512Vec8 log(const AVX512Vec8 );

		// SVML exp
		static inline AVX512Vec8 exp(const AVX512Vec8 );

		// SVML atan2
		static inline AVX512Vec8 atan2(const AVX512Vec8 ,
						 const AVX512Vec8 );

		// SVML hypot
		static inline AVX512Vec8 hypot(const AVX512Vec8 ,
						 const AVX512Vec8  );

		// FMA functions

		// fmadd
		static inline AVX512Vec8 fmadd(const AVX512Vec8 ,
						 const AVX512Vec8 ,
						 const AVX512Vec8 );

		// fmadsubb
		static inline AVX512Vec8 fmadsubb(const AVX512Vec8 ,
						    const AVX512Vec8 ,
						    const AVX512Vec8 );

		// fmsub
		static inline AVX512Vec8 fmsub(const AVX512Vec8 ,
						 const AVX512Vec8 ,
						 const AVX512Vec8 );

		// fmsubadd
		static inline AVX512Vec8 fmsubadd(const AVX512Vec8 ,
						    const AVX512Vec8 ,
						    const AVX512Vec8 );

		// fnmadd
		static inline AVX512Vec8 fnmadd(const AVX512Vec8 ,
						  const AVX512Vec8 ,
						  const AVX512Vec8 );

		// fnmsub
		static inline AVX512Vec8 fnmsub(const AVX512Vec8 ,
						  const AVX512Vec8 ,
						  const AVX512Vec8 );

		//
		// Global (namespace) operators
		//

		// C = A+B, vector + vector
		static inline AVX512Vec8 operator+(const AVX512Vec8 ,
						     const AVX512Vec8 );

		// C = A+B, vector + scalar
		static inline AVX512Vec8 operator+(const AVX512Vec8 ,
						     const double);

		// C = A+B, scalar + vector
		static inline AVX512Vec8 operator+(const double,
						     const AVX512Vec8 );

		// A = A+B (in-place)
		static inline AVX512Vec8 operator+=(AVX512Vec8 ,
						      const AVX512Vec8 );

		// C = A-B, vector-vector
		static inline AVX512Vec8 operator-(const AVX512Vec8 ,
						     const AVX512Vec8 );

		// C = A-B, vector - scalar
		static inline AVX512Vec8 operator-(const AVX512Vec8 ,
						     const double);

		// C = A-B, scalar - vector
		static inline AVX512Vec8 operator-(const double,
						     const AVX512Vec8 );

		// A = A-B (in-place)
		static inline AVX512Vec8 operator-=(AVX512Vec8 ,
						      const AVX512Vec8 );

		// C = A*B, vector * vector
		static inline AVX512Vec8 operator*(const AVX512Vec8 ,
						      const AVX512Vec8 );

		// C = A*B, vector * scalar
		static inline AVX512Vec8 operator*(const AVX512Vec8 ,
						     const double);

		// C = A*B, scalar * vector
		static inline AVX512Vec8 operator*(const double,
						     const AVX512Vec8 );

		// A = A*B (in-place)
		static inline AVX512Vec8 operator*=(AVX512Vec8 ,
						      const AVX512Vec8 );

		// C = A/B, vector / vector
		static inline AVX512Vec8 operator/(const AVX512Vec8 ,
						     const AVX512Vec8 );

		// C = A/B, vector / scalar
		static inline AVX512Vec8 operator/(const AVX512Vec8 ,
						     const double);

		// C = A/B, scalar / vector
		static inline AVX512Vec8 operator/(const double,
						     const AVX512Vec8 );

		// A = A/B (in-place)
		static inline AVX512Vec8 operator/=(AVX512Vec8 ,
						      const AVX512Vec8 ;

		// C = A==B, C is of type __mmask8
		static inline __mmask8 operator==(const AVX512Vec8 ,
						  const AVX512Vec8 );

		// C = A!=B, C is of type __mmask8
		static inline __mmask8 operator!=(const AVX512Vec8 ,
						  const AVX512Vec8 );

		// C = A>B, C is of type __mmask8
		static inline __mmask8 operator>(const AVX512Vec8 ,
						 const AVX512Vec8 );

		// C = A<B, C is of type __mmask8
		static inline __mmask8 operator<(const AVX512Vec8 ,
						 const AVX512Vec8 );

		// C = A>=B, C is of type __mmask8
		static inline __mmask8 operator>=(const AVX512Vec8 ,
						 const AVX512Vec8 );

		// C = A<=B, C is of type __mmask8
		static inline __mmask8 operator<=(const AVX512Vec8 ,
						  const AVX512Vec8 );

		// C = A&B
		static inline AVX512Vec8 operator&(const AVX512Vec8 ,
						     const AVX512Vec8 );

		// A = A&B
		static inline AVX512Vec8 operator&=(AVX512Vec8 ,
						      const AVX512Vec8 );

		// C = A|B
		static inline AVX512Vec8 operator|(const AVX512Vec8 ,
						     const AVX512Vec8 );

		// A = A|B
		static inline AVX512Vec8 operator|=(AVX512Vec8 ,
						       const AVX512Vec8 );

		// C = A^B
		static inline AVX512Vec8 operator^(const AVX512Vec8 ,
						    const AVX512Vec8 );

		// A = A^B
		static inline AVX512Vec8 operator^=(AVX512Vec8 ,
						      const AVX512Vec8 );

		// A = A + 1.0
		static inline AVX512Vec8 operator++(AVX512Vec8 );

		// A = A - 1.0
		static inline AVX512Vec8 operator--(AVX512Vec8);


#include "GMS_avx512vecf64.inl"
	}
}



#endif /*__GMS_AVX512VECF64_H__*/
