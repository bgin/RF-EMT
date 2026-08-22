
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

#ifndef __GMS_NRLMSISE00_SPLINE_INTEGRATORS_H_21_03_16__
#define __GMS_NRLMSISE00_SPLINE_INTEGRATORS_H_21_03_16__


#include "TypeTraits.hpp"

namespace ad = atmosphere::detail;
namespace  atmosphere {

	/* ------------------------------------------------------------------- */
	/* ------------------------------- SPLINI ---------------------------- */
	/* ------------------------------------------------------------------- */
	template<typename T> struct SPLINE_INTEGRATOR {
               
		static std::enable_if<ad::is_single_precision<T>::value ||
		                               ad::is_double_precision<T>::value, void>::type
	 	                                splinint(_In_ const T*__restrict, _In_ const T* __restrict, _In_ const T* __restrict, 
	 	                                                        _In_ const int, _In_ const T, _Out_ T*__restrict)noexcept(false);

	};

	/* ------------------------------------------------------------------- */
	/* ------------------------------- SPLINT ---------------------------- */
	/* ------------------------------------------------------------------- */
	template<typename T> struct SPLINE_INTERPOLATOR {

		static std::enable_if<ad::is_single_precision<T>::value || 
		                    ad::is_double_precision<T>::value, void>::type
	                          	splinterp(_In_ const T* __restrict, _In_ const T* __restrict, 
		                            _In_ const T* __restrict, _In_ const int, _In_ const T, _Out_  T*__restrict)noexcept(false);
	};

	/* ------------------------------------------------------------------- */
	/* ------------------------------- SPLINE ---------------------------- */
	/* ------------------------------------------------------------------- */
	template<typename T> struct SPLINE_CUBIC_DERIV {

		static std::enable_if<ad::is_single_precision<T>::value || 
		                      ad::is_double_precision<T>::value, void>::type
	                          	splinderiv(_In_ const T*__restrict, _In_ const T*__restrict, _In_ const int, 
		                                        _In_ const T, _In_ const T, _Out_ T*__restrict)noexcept(false);
	};
#include "GMS_NRLMSISE00_SPLINE_INTEGRATORS.inl"
}
#endif /*__GMS_NRLMSISE00_SPLINE_INTEGRATORS_H_21_03_16__*/
