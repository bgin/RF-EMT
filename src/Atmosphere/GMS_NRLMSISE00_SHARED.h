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
