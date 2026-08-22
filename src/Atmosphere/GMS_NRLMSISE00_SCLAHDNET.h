
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

#ifndef __GMS_NRLMSISE00_SCLAH_DNET_H_18_03_16__
#define __GMS_NRLMSISE00_SCLAH_DNET_H_18_03_16__

/* -------------------------------------------------------------------- */
/* ---------  N R L M S I S E - 0 0    M O D E L    2 0 0 1  ---------- */
/* -------------------------------------------------------------------- */

/*  This file is part of C++ port of NRLMSISE-00 implemented in C.
* @File  NRLMSISE00_SCLAH_DNET.h
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
namespace ad = atmosphere::detail;
namespace ang = atmosphere::nrlmsise00_globals;
namespace   atmosphere {

	/* ------------------------------------------------------------------- */
	/* ------------------------------- SCALH ----------------------------- */
	/* ------------------------------------------------------------------- */
	template<typename T>   struct  SCALH {

		std::enable_if<ad::is_single_precision<T>::value || ad::is_double_precision<T>::value, T>::type
		operator()(_In_ const T, _In_ const T, _In_ const T);
	};

	/* ------------------------------------------------------------------- */
	/* -------------------------------- DNET ----------------------------- */
	/* ------------------------------------------------------------------- */
	template<typename T>  struct DNET {
		/*       TURBOPAUSE CORRECTION FOR MSIS MODELS
		*        Root mean density
		*         DD - diffusive density
		*         DM - full mixed density
		*         ZHM - transition scale length
		*         XMM - full mixed molecular weight
		*         XM  - species molecular weight
		*         DNET - combined density
		*/
		std::enable_if<ad::is_single_precision<T>::value || ad::is_double_precision<T>::value, T>::type
		operator()(_In_  T, _In_ const T, _In_ const T, _In_ const T, _In_ const T);
	};

#include "GMS_NRLMSISE00_SCLAHDNET.inl"
}
#endif  /*__GMS_NRLMSISE00_SCLAH_DNET_H_18_03_16__*/
