#ifndef __GMS_NRLMSISE00_CCOR_H_18_03_16__
#define __GMS_NRLMSISE00_CCOR_H_18_03_16__
/* -------------------------------------------------------------------- */
/* ---------  N R L M S I S E - 0 0    M O D E L    2 0 0 1  ---------- */
/* -------------------------------------------------------------------- */

/*  This file is part of C++ port of NRLMSISE-00 implemented in C.
* @File  NRLMSISE00_CCOR.h
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

#include "TypeTraits.hpp"

namespace  ad = atmosphere::detail;

namespace  atmosphere {

	/* ------------------------------------------------------------------- */
	/* ------------------------------ CCOR ------------------------------- */
	/* ------------------------------------------------------------------- */

	template<typename T>  struct CCOR {

		std::enable_if<ad::is_single_precision<T>::value || ad::is_double_precision<T>::value, T>::type
		operator()(_In_ const T, _In_ const T, _In_ const T, _In_ const T);
	};

	/* ------------------------------------------------------------------- */
	/* ------------------------------ CCOR2 ------------------------------- */
	/* ------------------------------------------------------------------- */

	template<typename T>  struct CCOR2 {

		std::enable_if<ad::is_single_precision<T>::value || ad::is_double_precision<T>::value, T>::type
		operator()(_In_ const T, _In_ const T, _In_ const T, _In_ const T, _In_ const T);
	};
#include "GMS_NRLMSISE00_CCOR.inl"
}
#endif  /*__GMS_NRMLSISE00_CCOR_H_18_03_16__*/
