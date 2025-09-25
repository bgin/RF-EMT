#ifndef __GMS_NRLMSISE00_INOUT_H_24_03_16__
#define __GMS_NRLMSISE00_INOUT_H_24_03_16__

/* -------------------------------------------------------------------- */
/* ---------  N R L M S I S E - 0 0    M O D E L    2 0 0 1  ---------- */
/* -------------------------------------------------------------------- */

/*  This file is part of C++ port of NRLMSISE-00 implemented in C.
* @File  NRLMSISE00_INOUT.h
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

namespace   atmosphere {

	/* ------------------------------------------------------------------- */
	/* ------------------------------- INPUT ----------------------------- */
	/* ------------------------------------------------------------------- */

	template<typename T> struct NRLMSISE_FLAGS {

		int switches[24];
		T sw[24];
		T swc[24];
	};
	/*
	*   Switches: to turn on and off particular variations use these switches.
	*   0 is off, 1 is on, and 2 is main effects off but cross terms on.
	*
	*   Standard values are 0 for switch 0 and 1 for switches 1 to 23. The
	*   array "switches" needs to be set accordingly by the calling program.
	*   The arrays sw and swc are set internally.
	*
	*   switches[i]:
	*    i - explanation
	*   -----------------
	*    0 - output in meters and kilograms instead of centimeters and grams
	*    1 - F10.7 effect on mean
	*    2 - time independent
	*    3 - symmetrical annual
	*    4 - symmetrical semiannual
	*    5 - asymmetrical annual
	*    6 - asymmetrical semiannual
	*    7 - diurnal
	*    8 - semidiurnal
	*    9 - daily ap [when this is set to -1 (!) the pointer
	*                  ap_a in struct nrlmsise_input must
	*                  point to a struct ap_array]
	*   10 - all UT/long effects
	*   11 - longitudinal
	*   12 - UT and mixed UT/long
	*   13 - mixed AP/UT/LONG
	*   14 - terdiurnal
	*   15 - departures from diffusive equilibrium
	*   16 - all TINF var
	*   17 - all TLB var
	*   18 - all TN1 var
	*   19 - all S var
	*   20 - all TN2 var
	*   21 - all NLB var
	*   22 - all TN3 var
	*   23 - turbo scale height var
	*/

	template<typename T> struct AP_VALUES {

		T a[7];
	};
	/* Array containing the following magnetic values:
	*   0 : daily AP
	*   1 : 3 hr AP index for current time
	*   2 : 3 hr AP index for 3 hrs before current time
	*   3 : 3 hr AP index for 6 hrs before current time
	*   4 : 3 hr AP index for 9 hrs before current time
	*   5 : Average of eight 3 hr AP indicies from 12 to 33 hrs
	*           prior to current time
	*   6 : Average of eight 3 hr AP indicies from 36 to 57 hrs
	*           prior to current time
	*/

	template<typename T> struct NRLMSISE_INPUT {
		int year; // year, currently ignored.//
		int doy;  // day of year //
		T sec;    // seconds in day //
		T alt;    // altitude in kilometers //
		T g_lat;  // geodetic latitude //
		T g_long; // geodetic longtitude //
		T lst;    // local apparent solar time //
		T F107A;  // 81 day average of F10.7 flux (centered on doy) //
		T F107;   // daily F10.7 flux for previous day //
		T ap;     // magnetic index(daily)
		AP_VALUES<T> ap_values;
	};

	/*
	*   NOTES ON INPUT VARIABLES:
	*      UT, Local Time, and Longitude are used independently in the
	*      model and are not of equal importance for every situation.
	*      For the most physically realistic calculation these three
	*      variables should be consistent (lst=sec/3600 + g_long/15).
	*      The Equation of Time departures from the above formula
	*      for apparent local time can be included if available but
	*      are of minor importance.
	*
	*      f107 and f107A values used to generate the model correspond
	*      to the 10.7 cm radio flux at the actual distance of the Earth
	*      from the Sun rather than the radio flux at 1 AU. The following
	*      site provides both classes of values:
	*      ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
	*
	*      f107, f107A, and ap effects are neither large nor well
	*      established below 80 km and these parameters should be set to
	*      150., 150., and 4. respectively.
	*/

	/* ------------------------------------------------------------------- */
	/* ------------------------------ OUTPUT ----------------------------- */
	/* ------------------------------------------------------------------- */

	template<typename T> struct NRLMSISE_OUTPUT {
		T d[9]; //** densities **//
		T t[2]; //** temperatures **// 
	};
	/*
	*   OUTPUT VARIABLES:
	*      d[0] - HE NUMBER DENSITY(CM-3)
	*      d[1] - O NUMBER DENSITY(CM-3)
	*      d[2] - N2 NUMBER DENSITY(CM-3)
	*      d[3] - O2 NUMBER DENSITY(CM-3)
	*      d[4] - AR NUMBER DENSITY(CM-3)
	*      d[5] - TOTAL MASS DENSITY(GM/CM3) [includes d[8] in td7d]
	*      d[6] - H NUMBER DENSITY(CM-3)
	*      d[7] - N NUMBER DENSITY(CM-3)
	*      d[8] - Anomalous oxygen NUMBER DENSITY(CM-3)
	*      t[0] - EXOSPHERIC TEMPERATURE
	*      t[1] - TEMPERATURE AT ALT
	*
	*
	*      O, H, and N are set to zero below 72.5 km
	*
	*      t[0], Exospheric temperature, is set to global average for
	*      altitudes below 120 km. The 120 km gradient is left at global
	*      average value for altitudes below 72 km.
	*
	*      d[5], TOTAL MASS DENSITY, is NOT the same for subroutines GTD7
	*      and GTD7D
	*/
}

#endif /*__GMS_NRLMSISE00_INOUT_H_24_03_16__*/
