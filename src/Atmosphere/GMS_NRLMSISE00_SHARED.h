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

#ifndef __GMS_NRLMSISE00_SHARED_H_21_03_16__
#define __GMS_NRLMSISE00_SHARED_H_21_03_16__


/*
           Shared Variables
*/

namespace atmosphere {

	namespace nrlmsise00_globals {

		// As evil as it is:)

		/****************************************************
		            PARMB  Globals
					double precision by default.
		*****************************************************/
		  double gsurf;
		  double re;

		/****************************************************
		          GTS3C Global variables defaulted to double
				  precision.
		*****************************************************/
		 double dd;

		/***************************************************
		          DMIX Global variables.
				  double precision by default.
		****************************************************/
		 double dm04;
		 double dm16;
		 double dm28;
		 double dm32;
		 double dm40;
		 double dm01;
		 double dm14;

		/****************************************************
		             MESO7 Global static arrays
					 double precision by default
		*****************************************************/
		 double meso_tn1[5];
		 double meso_tn2[4];
		 double meso_tn3[5];
		 double meso_tgn1[2];
		 double meso_tgn2[2];
		 double meso_tgn3[2];

		/****************************************************
		            LPOLY Global static arrays and scalar 
					variables, double precision by default.
		*****************************************************/
		 double dfa;
		 double plg[4][9];
		 double ctloc;
		 double stloc;
		 double c2tloc;
		 double s2tloc;
		 double c3tloc;
		 double s3tloc;
		 double apdf;
		 double apt[4];
	}
}
#endif /*__GMS_NRLMSISE00_SHARED_H_21_03_16__*/
