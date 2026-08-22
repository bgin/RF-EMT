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
#ifndef __GMS_FPEXCEPTIONS_H__
#define __GMS_FPEXCEPTIONS_H__



namespace file_info {

  const unsigned int gGMS_FPEXCEPTIONS_MAJOR = 1U;
  const unsigned int gGMS_FPEXCEPTIONS_MINOR = 1U;
  const unsigned int gGMS_FPEXCEPTIONS_MICRO = 0U;
  const unsigned int gGMS_FPEXCEPTIONS_FULLVER =
    1000U*gGMS_FPEXCEPTIONS_MAJOR + 100U*gGMS_FPEXCEPTIONS_MINOR + 10U*gGMS_FPEXCEPTIONS_MICRO;
  const char * const pgGMS_FPEXCEPTIONS_CREATE_DATE = "02-10-2019 19:04 +00200 (WED 02 OCT 2019 GMT+2)";
  const char * const pgGMS_FPEXCEPTIONS_BUILD_DATE  = __DATE__ ":" __TIME__
  const char * const pgGMS_FPEXCEPTIONS_AUTHOR      = "Programmer: Bernard Gingold contact: beniekg@gmail.com";
  const char * const pgGMS_FPEXCEPTIONS_SYNOPSIS    = "Access floating-point environment state."
}

#include <cstdint>

namespace gms {
	namespace math {

		/*
		   @Description: Checks for existance of denormal values in
		                 input range i.e. for some |x| < FLT_MIN
		                 Upon detecting error condition 'domain range'
		                 third argument accumulates number of denormal
						 occurrences in array.
		                 Fast version without input checking
		@Params:  1st source array 1D holding floats (32-bit)
		@Params:  none
		@Params:  length of  source array.
		@Params:  none
		@Params:  pointer to variable holding 'Domain Range'error value.
		@Returns: Nothing
		@Throws:  Nothing
		@Calls:   fpclassify
		*/
		void is_denormalf32_present(const float * __restrict,
					    const int32_t, 
					    uint32_t *,
					    const bool);

		/*
		    @Description: Checks for existance of denormal values in
		                  input range i.e. for some |x| < FLT_MIN
		                  Upon detecting error condition 'domain range'
		                  third argument accumulates number of denormal
						  occurrences in array.
		                  Fast version without input checking
		@Params:  1st source array 2D holding floats (32-bit)
		@Params:  none
		@Params:  length of  source array.
		@Params:  none
		@Params:  pointer to variable holding 'Domain Range'error value.
		@Returns: Nothing
		@Throws:  Nothing
		@Calls:   fpclassify
		*/
		void is_denormalf32_present(const float * __restrict,
					    const int32_t,
					    const int32_t,
					    uint32_t *,
					    const bool );

		/*
			@Description: Checks for existance of denormal values in
						  input range i.e. for some |x| < FLT_MIN
						  Upon detecting error condition 'domain range'
						  fifth argument accumulates number of denormal
						  occurrences in array.
		                  Fast version without input checking
		@Params:  1st source array 3D holding floats (32-bit)
		@Params:  none
		@Params:  length of  source array.
		@Params:  none
		@Params:  pointer to variable holding 'Domain Range'error value.
		@Returns: Nothing
		@Throws:  Nothing
		@Calls:   fpclassify
		*/
		void is_denormalf32_present(const float * __restrict,
					    const int32_t,
					    const int32_t,
					    const int32_t,
					    uint32_t *,
					    const bool );

		/*
		    @Description: Checks for existance of denormal values in
		                  input range i.e. for some |x| < FLT_MIN
		                  Upon detecting error condition 'domain range'
		                  sixth argument accumulates number of denormal
		                  occurrences in array.
		Fast version without input checking
		@Params:  1st source array 4D holding floats (32-bit)
		@Params:  none
		@Params:  length of  source array.
		@Params:  none
		@Params:  pointer to variable holding 'Domain Range'error value.
		@Returns: Nothing
		@Throws:  Nothing
		@Calls:   fpclassify
		*/
		void is_denormalf32_present(const float * __restrict,
					    const int32_t,
					    const int32_t,
					    const int32_t,
					    const int32_t,
					    uint32_t *,
					    const bool );

		/*
			@Description: Checks for existance of denormal values in
						  input range i.e. for some |x| < DBL_MIN
						  Upon detecting error condition 'domain range'
						  third argument accumulates number of denormal
						  occurrences in array.
						  Fast version without input checking
		@Params:  1st source array 1D holding doubles (64-bit)
		@Params:  none
		@Params:  length of  source array.
		@Params:  none
		@Params:  pointer to variable holding 'Domain Range'error value.
		@Returns: Nothing
		@Throws:  Nothing
		@Calls:   fpclassify
		*/
		void is_denormalf64_present(const double * __restrict,
					    const int32_t,
					    uint32_t *,
					    const bool );


		/*
			@Description: Checks for existance of denormal values in
						  input range i.e. for some |x| < DBL_MIN
						  Upon detecting error condition 'domain range'
						  fourth argument accumulates number of denormal
						  occurrences in array.
		Fast version without input checking
		@Params:  1st source array 2D holding doubles (64-bit)
		@Params:  none
		@Params:  length of  source array.
		@Params:  none
		@Params:  pointer to variable holding 'Domain Range'error value.
		@Returns: Nothing
		@Throws:  Nothing
		@Calls:   fpclassify
		*/
		void is_denormalf64_present(const double * __restrict,
					    const int32_t,
					    const int32_t,
					    uint32_t *,
					    const bool );

		/*
			@Description: Checks for existance of denormal values in
						  input range i.e. for some |x| < DBL_MIN
						  Upon detecting error condition 'domain range'
						  fivth argument accumulates number of denormal
						  occurrences in array.
						  Fast version without input checking
		@Params:  1st source array 3D holding doubles (64-bit)
		@Params:  none
		@Params:  length of  source array.
		@Params:  none
		@Params:  pointer to variable holding 'Domain Range'error value.
		@Returns: Nothing
		@Throws:  Nothing
		@Calls:   fpclassify
		*/
		void is_denormalf64_present(const double * __restrict,
					    const int32_t ,
					    const int32_t,
					    const int32_t,
					    uint32_t *,
					    const bool );

		/*
			@Description: Checks for existance of denormal values in
						  input range i.e. for some |x| < DBL_MIN
						  Upon detecting error condition 'domain range'
						  sixth argument accumulates number of denormal
						  occurrences in array.
						  Fast version without input checking
		@Params:  1st source array 3D holding doubles (64-bit)
		@Params:  none
		@Params:  length of  source array.
		@Params:  none
		@Params:  pointer to variable holding 'Domain Range'error value.
		@Returns: Nothing
		@Throws:  Nothing
		@Calls:   fpclassify
		*/
		void is_denormalf64_present(const double * __restrict,
					    const int32_t,
					    const int32_t,
					    const int32_t,
					    const int32_t,
					    uint32_t * ,
					    const bool);

		/*
			@Description: Checks for existance of invalid values
						  in domain range of dimension 1D
						  Upon detecting error condition 'domain range'
						  third argument accumulates number of denormal
						  occurrences in array.
						  Checks for:
						  INF
						  DENORMAL
						  NAN
						  Fast version without input checking
		@Params:  1st source array 1D holding floats (32-bit)
		@Params:  none
		@Params:  length of  source array.
		@Params:  none
		@Params:  pointer to variable holding 'Domain Range'error value.
		@Returns: Nothing
		@Throws:  Nothing
		@Calls:   fpclassify
		*/
		void is_abnormalf32(   const float * __restrict,
				       const int32_t,
				       uint32_t * ,
				       const bool,
				       const uint32_t);

		/*
			@Description:	Checks for existance of invalid values
							in domain range of dimension 2D
							Upon detecting error condition 'domain range'
							fourth argument accumulates number of denormal
							occurrences in array.
		Checks for:
		INF
		DENORMAL
		NAN
		Fast version without input checking
		@Params:  1st source array 2D holding floats (32-bit)
		@Params:  none
		@Params:  length of  source array.
		@Params:  none
		@Params:  pointer to variable holding 'Domain Range'error value.
		@Returns: Nothing
		@Throws:  Nothing
		@Calls:   fpclassify
		*/
		void is_abnormalf32(const float * __restrict,
				    const int32_t,
				    const int32_t ,
				    uint32_t * ,
				    const bool,
				    const uint32_t);

		/*
		@Description:	Checks for existance of invalid values
					    in domain range of dimension 3D
						Upon detecting error condition 'domain range'
						fifth argument accumulates number of denormal
						occurrences in array.
		Checks for:
		INF
		DENORMAL
		NAN
		Fast version without input checking
		@Params:  1st source array 3D holding floats (32-bit)
		@Params:  none
		@Params:  length of  source array.
		@Params:  none
		@Params:  pointer to variable holding 'Domain Range'error value.
		@Returns: Nothing
		@Throws:  Nothing
		@Calls:   fpclassify
		*/
		void is_abnormalf32(const float * __restrict,
				    const int32_t,
				    const int32_t,
				    const int32_t,
				    uint32_t * ,
				    const bool,
				    const uint32_t );

		/*
		@Description:	Checks for existance of invalid values
						in domain range of dimension 4D
						Upon detecting error condition 'domain range'
						sixth argument accumulates number of denormal
						occurrences in array.
		Checks for:
		INF
		DENORMAL
		NAN
		Fast version without input checking
		@Params:  1st source array 4D holding floats (32-bit)
		@Params:  none
		@Params:  length of  source array.
		@Params:  none
		@Params:  pointer to variable holding 'Domain Range'error value.
		@Returns: Nothing
		@Throws:  Nothing
		@Calls:   fpclassify
		*/
		void is_abnormalf32(const float * __restrict,
				    const int32_t,
				    const int32_t,
				    const int32_t,
				    const int32_t,
				    uint32_t * ,
				    const bool,
				    const uint32_t );

		/*
		@Description: Checks for existance of invalid values
				      in domain range of dimension 1D
					  Upon detecting error condition 'domain range'
					  third argument accumulates number of denormal
		`			  occurrences in array.
		Checks for:
		INF
		DENORMAL
		NAN
		Fast version without input checking
		@Params:  1st source array 1D holding doubles (64-bit)
		@Params:  none
		@Params:  length of  source array.
		@Params:  none
		@Params:  pointer to variable holding 'Domain Range'error value.
		@Returns: Nothing
		@Throws:  Nothing
		@Calls:   fpclassify
		*/
		void is_abnormalf64(const double * __restrict,
				    const int32_t ,
				    uint32_t *,
				    const bool,
				    const uint32_t );

		/*
		@Description: Checks for existance of invalid values
					  in domain range of dimension 2D
					  Upon detecting error condition 'domain range'
					  fourth argument accumulates number of denormal
					  occurrences in array.
		Checks for:
		INF
		DENORMAL
		NAN
		Fast version without input checking
		@Params:  1st source array 2D holding doubles (64-bit)
		@Params:  none
		@Params:  length of  source array.
		@Params:  none
		@Params:  pointer to variable holding 'Domain Range'error value.
		@Returns: Nothing
		@Throws:  Nothing
		@Calls:   fpclassify
		*/
		void is_abnormalf64(const double * __restrict,
				    const int32_t,
				    const int32_t,
				    uint32_t * ,
				    const bool,
				    const uint32_t );

		/*
		@Description: Checks for existance of invalid values
					  in domain range of dimension 3D
					  Upon detecting error condition 'domain range'
					  fifth argument accumulates number of denormal
					  occurrences in array.
		Checks for:
		INF
		DENORMAL
		NAN
		Fast version without input checking
		@Params:  1st source array 3D holding doubles (64-bit)
		@Params:  none
		@Params:  length of  source array.
		@Params:  none
		@Params:  pointer to variable holding 'Domain Range'error value.
		@Returns: Nothing
		@Throws:  Nothing
		@Calls:   fpclassify
		*/
		void is_abnormalf64(const double * __restrict,
				    const int32_t ,
				    const int32_t ,
				    const int32_t,
				    uint32_t * ,
				    const bool,
				    const uint32_t );

		/*
		@Description: Checks for existance of invalid values
					  in domain range of dimension 4D
					  Upon detecting error condition 'domain range'
					  sixth argument accumulates number of denormal
					  occurrences in array.
		Checks for:
		INF
		DENORMAL
		NAN
		Fast version without input checking
		@Params:  1st source array 4D holding doubles (64-bit)
		@Params:  none
		@Params:  length of  source array.
		@Params:  none
		@Params:  pointer to variable holding 'Domain Range'error value.
		@Returns: Nothing
		@Throws:  Nothing
		@Calls:   fpclassify
		*/
		void is_abnormalf64(const double * __restrict,
				    const int32_t ,
				    const int32_t,
				    const int32_t,
				    const int32_t,
				    uint32_t * ,
				    const bool,
				    const uint32_t );

		/*
		@Description: Clears all floating-point state exceptions
		Exception cleared:
		FE_DENORMAL, FE_INVALID, FE_INEXACT, FE_UNDERFLOW
		FE_OVERFLOW.
		Scalar version
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Returns: integer 'err' which indicates success or error as a return
		value from library 'feclearexcept' function
		Non-zero value means error.
		@Throws:  Nothing
		@Calls:   'feclearexcept'
		*/
		int32_t clear_fpexcepts(void);

		/*
		@Description: Clears only FE_DENORMAL exception
		Exception cleared:
		FE_DENORMAL
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Returns: integer 'err' which indicates success or error as a return
		value from library 'feclearexcept' function
		Non-zero value means error.
		@Throws:  Nothing
		@Calls:   'feclearexcept'
		*/
		int32_t clear_fedenormal(void);

		/*
		@Description: Clears only FE_INEXACT exception
		Exception cleared:
		FE_INEXACT
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Returns: integer 'err' which indicates success or error as a return
		value from library 'feclearexcept' function
		Non-zero value means error.
		@Throws:  Nothing
		@Calls:   'feclearexcept'
		*/
		int32_t clear_feinexact(void);

		/*
		@Description: Clears only FE_INVALID exception
		Exception cleared:
		FE_INVALID
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Returns: integer 'err' which indicates success or error as a return
		value from library 'feclearexcept' function
		Non-zero value means error.
		@Throws:  Nothing
		@Calls:   'feclearexcept'
		*/
		int32_t clear_feinvalid(void);

		/*
		@Description: Clears only FE_DIVBYZERO exception
		Exception cleared:
		FE_DIVBYZERO
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Returns: integer 'err' which indicates success or error as a return
		value from library 'feclearexcept' function
		Non-zero value means error.
		@Throws:  Nothing
		@Calls:   'feclearexcept'
		*/
		int32_t clear_fedivbyzero(void);

		/*
		@Description: Clears only FE_OVERFLOW exception
		Exception cleared:
		FE_OVERFLOW
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Returns: integer 'err' which indicates success or error as a return
		value from library 'feclearexcept' function
		Non-zero value means error.
		@Throws:  Nothing
		@Calls:   'feclearexcept'
		*/
		int32_t clear_feoverflow(void);

		/*
		@Description: Clears only FE_UNDERFLOW exception
		Exception cleared:
		FE_UNDERFLOW
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Returns: integer 'err' which indicates success or error as a return
		value from library 'feclearexcept' function
		Non-zero value means error.
		@Throws:  Nothing
		@Calls:   'feclearexcept'
		*/
		int32_t clear_feunderflow(void);

		/*
		@Description: Tests all floating-point exceptions.

		@Params:  All 7 floating-point exception types (exception values must be or'ed).
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Returns: integer 'val' which indicates success or error as a return
		value from library 'fetestexcept' function

		@Throws:  Nothing
		@Calls:   'fetestexcept'
		*/
		int32_t test_feexcepts(const int32_t);

		/*
		@Description: Tests for existance of FE_INVALID exception.

		@Params:  argument FE_INVALID macro.
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Returns: integer 'val' which indicates success or error as a return
		value from library 'fetestexcept' function

		@Throws:  Nothing
		@Calls:   'fetestexcept'
		*/
		int32_t test_feinvalid(const int32_t);

		/*
		@Description: Tests for existance of FE_INEXACT exception.

		@Params:  argument FE_INEXACT macro.
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Returns: integer 'val' which indicates success or error as a return
		value from library 'fetestexcept' function

		@Throws:  Nothing
		@Calls:   'fetestexcept'
		*/
		int32_t test_feinexact(_In_ const int32_t);

		/*
		@Description: Tests for existance of FE_DIVBYZERO exception.

		@Params:  argument FE_DIVBYZERO macro.
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Returns: integer 'val' which indicates success or error as a return
		value from library 'fetestexcept' function

		@Throws:  Nothing
		@Calls:   'fetestexcept'
		*/
		int32_t test_fedivbyzero(const int32_t);

		/*
		@Description: Tests for existance of FE_UNNORMAL exception.

		@Params:  argument FE_UNNORMAL macro.
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Returns: integer 'val' which indicates success or error as a return
		value from library 'fetestexcept' function

		@Throws:  Nothing
		@Calls:   'fetestexcept'
		*/
		int32_t test_feunormal(const int32_t);

		/*
		@Description: Tests for existance of FE_OVERFLOW exception.

		@Params:  argument FE_OVERFLOW macro.
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Returns: integer 'val' which indicates success or error as a return
		value from library 'fetestexcept' function

		@Throws:  Nothing
		@Calls:   'fetestexcept'
		*/
		int32_t test_feoverflow(const int32_t);

		/*
		@Description: Tests for existance of FE_UNDERFLOW exception.

		@Params:  argument FE_UNDERFLOW macro.
		@Params:  none
		@Params:  none
		@Params:  none
		@Params:  none
		@Returns: integer 'val' which indicates success or error as a return
		value from library 'fetestexcept' function

		@Throws:  Nothing
		@Calls:   'fetestexcept'
		*/
		int32_t test_feunderflow(const int32_t);


		void rise_fedenormal(const bool, 
				     const int32_t,
				     int32_t &);

		void rise_feinvalid(const bool,
				    const int32_t,
				    int32_t &);

		void rise_feinexact(const bool,
				    const int32_t,
				    int32_t &);

		void rise_fedivbyzero(const bool,
				      const int32_t,
				      int32_t &);

		void rise_feoverflow(const bool,
				     const int32_t,
				     int32_t &);

		void rise_feundeflow(const bool,
				     const int32_t,
				     int32_t & );
							



	} // math
} 

#endif /*__GMS_FPEXCEPTIONS_H__*/
