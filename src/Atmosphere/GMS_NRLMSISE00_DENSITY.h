#ifndef __GMS_NRLMSISE00_DENSITY_H_22_03_16__
#define __GMS_NRLMSISE00_DENSITY_H_22_03_16__

/* -------------------------------------------------------------------- */
/* ---------  N R L M S I S E - 0 0    M O D E L    2 0 0 1  ---------- */
/* -------------------------------------------------------------------- */

/*  This file is part of C++ port of NRLMSISE-00 implemented in C.
* @File  NRLMSISE00_DENSITY.h
* The NRLMSISE-00 model was developed by Mike Picone, Alan Hedin, and
* Doug Drob. They also wrote a NRLMSISE-00 distribution package in
* FORTRAN which is available at
* http://uap-www.nrl.navy.mil/models_web/msis/msis_home.htm
*
* Dominik Brodowski implemented and maintains this C version. You can
* reach him at mail@brodo.de. See the file "DOCUMENTATION" for details,
* and check http://www.brodo.de/english/pub/nrlmsise/index.html for
* updated releases of this package.
*
* Adapted from the work of Dominik Brodowski by Bernard Gingold
*/


#include "GMS_NRLMSISE00_SHARED.h"
#include "GMS_NRLMSISE00_SPLINE_INTEGRATORS.h"

namespace ad = atmosphere::detail;

namespace   atmosphere {

	/* ------------------------------------------------------------------- */
	/* ------------------------------- DENSM ----------------------------- */
	/* ------------------------------------------------------------------- */
	/* ------------------------------------------------------------------- */
	/* ------------------------------- DENSU ----------------------------- */
	/* ------------------------------------------------------------------- */

	template<typename T, class SPLINE = atmosphere::SPLINE_CUBIC_DERIV<T>,
		class SPLINT = atmosphere::SPLINE_INTERPOLATOR<T>, class SPLINI = 
			atmosphere::SPLINE_INTEGRATOR < T >> struct DENSM {

				/* inline helper function*/
				static inline  auto  zeta(_In_ const T, _In_ const T)->T;

				/* Calculate Temperature and Density Profiles for 
				   lower Atmosphere */
				std::enable_if<ad::is_single_precision<T>::value || ad::is_double_precision<T>::value, T>::type
					densprof(_In_ const T, _In_ const T, _In_ const T, _In_ T*, _In_ const int,
					_In_ const T*, _In_ const T*, _In_ const T*, _In_ const int, _In_ const T*,
					_In_ const T*, _In_ const T*);

				/*  Calculate Temperature and Density Profiles for MSIS model
				   New lower thermo polynomial
				*/
				std::enable_if<ad::is_single_precision<T>::value || ad::is_double_precision<T>::value, T>::type
					densuprof(_In_ const T, _In_ const T, _In_ const T, _In_ const T, _In_ const T,
					_In_ const T, _In_ const T*, _In_ const T, _In_ const T, _In_ const int, _In_ const T*,
					_In_ const T*, _In_ const T*);
			};






#include "GMS_NRLMSISE00_DENSITY.inl"
}
#endif /*__GMS_NRLMSISE00_DENSITY_H_22_03_16__*/
