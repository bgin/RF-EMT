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
#ifndef __GMS_AVX512VCPROCEDURES_H__
#define __GMS_AVX512VCPROCEDURES_H__ 121020191812

namespace file_info {


	const unsigned int gGMS_AVX512PROCEDURES_MAJOR = 1;

	const unsigned int gGMS_AVX512PROCEDURES_MINOR = 0;

	const unsigned int gGMS_AVX512PROCEDURES_MICRO = 1;

	const unsigned int gGMS_AVX512PROCEDURES_FULLVER = 
		1000U*gGMS_AVX512COMPLEX_MAJOR+100U*gGMS_AVX512COMPLEX_MINOR+10U*gGMS_AVX512COMPLEX_MICRO;

	const char * const pgGMS_AVX512PROCEDURES_CREATE_DATE = "12-10-2019 18:12 +00200 (SAT  12 OCT  2019 GMT+2)";

	const char * const pgGMS_AVX512PROCEDURES_BUILD_DATE = __DATE__":"__TIME__;

	const char * const pgGMS_AVX512PROCEDURES_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_AVX512PROCEDURES_SYNOPSIS = "Global (namespace) operator-like procedures.";
}


namespace gms {
	namespace math {

class AVX512VComplex1D;
			
		

		
		void
		avx512vcomplex_add(AVX512VComplex1D    &,
				   const AVX512VComplex1D &,
				   const AVX512VComplex1D &,
				   const bool ) noexcept(true);

		
		void
		avx512vcomplex_add(AVX512VComplex1D &,
				   const AVX512VComplex1D &,
				   const double * __restrict,
				   const bool) noexcept(true);

		
		void
		avx512vcomplex_sub(AVX512VComplex1D &,
				   const AVX512VComplex1D &,
				   const AVX512VComplex1D &,
				   const bool) noexcept(true);

	
		void
		avx512vcomplex_sub(AVX512VComplex1D &,
				   const AVX512VComplex1D &,
				   const double * __restrict,
				   const bool) noexcept(true);

	
		void
		avx512vcomplex_mul(AVX512VComplex1D &,
				   const AVX512VComplex1D &,
				   const AVX512VComplex1D &,
				   const bool) noexcept(true);

	
		void
		avx512vcomplex_mul(AVX512VComplex1D &,
				   const AVX512VComplex1D &,
				   const double * __restrict,
				   const bool) noexcept(true);

		void
		avx512vcomplex_div(AVX512VComplex1D	 &, 
				   const AVX512VComplex1D &,
				   const AVX512VComplex1D &,
				   const bool);
						

		
		void
		avx512vcomplex_div(AVX512VComplex1D    &, 
				   const AVX512VComplex1D &,
				   const double * __restrict,
				   const bool) noexcept(true);

	 
	}
}



#endif /*__GMS_AVX512VCPROCEDURES_H__*/
