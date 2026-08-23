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

#ifndef __GMS_RCS_HPP__
#define __GMS_RCS_HPP__ 230120221416

namespace file_info {

     const unsigned int GMS_RCS_MAJOR = 1;
     const unsigned int GMS_RCS_MINOR = 1;
     const unsigned int GMS_RCS_MICRO = 0;
     const unsigned int GMS_RCS_FULLVER =
       1000U*GMS_RCS_MAJOR+100U*GMS_RCS_MINOR+
       10U*GMS_RCS_MICRO;
     const char * const GMS_RCS_CREATION_DATE = "23-01-2022 14:16 +00200 (SUN 23 JAN 2022 GMT+2)";
     const char * const GMS_RCS_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_RCS_SYNOPSIS      = "Calculation of various aspects of the RCS."

}

#include <cstdint>
#include <cmath> // for double precision version
#include "GMS_config.h"
#include "GMS_cephes.h" // for single precision version

namespace gms {


             namespace radiolocation {


	           namespace {
                        // Free-space target RCS [m^2]
                        const float   sig0_r4 = 1.0f;
                        const double  sig0_r8 = 1.0;
			const float   PIr4   = 3.14159265358979323846264338328f;
			const double  PIr8   = 3.14159265358979323846264338328;
			const float   PI2r4  = 6.283185307179586476925286766559f;
			const double  PI2r8  = 6.283185307179586476925286766559;
			const float   PI4r4  = 12.566370614359172953850573533118f;
			const double  PI4r8  = 12.566370614359172953850573533118;
                        const float   zr4    = 0.017453292519943295769236907685f;
			const double  zr8    = 0.017453292519943295769236907685; //coeff deg-to-rad conversion
		      
		   }

		   /*
                       This function computes the empirical K parameter when unknown apriori.
                       Forward-scatterer path when beta=180.
                    */

		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float empirical_K_r4_1(const float sig0,  // m^2, RCS monostatic target
		                            const float A,     // m^2, Area of target projected to the normal radar beam
					    const float gamma) { // m,   wavelength
			   constexpr float PI24 = 0.74159265358979323846264338328f;
                           float K              = 0.0f;
			   const float A2       = A*A;
			   const float gamm2    = gamma*gamma*sig0;
			   const float t0       = PI4r4*(A2/gamm2);
			   K = cepehes_logf(t0)/PI24;
			   return (K);
		    }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double empirical_K_r4_1(const double sig0,  // m^2, RCS monostatic target
		                             const double A,     // m^2, Area of target projected to the normal radar beam
					     const double gamma) { // m,   wavelength
			   constexpr double PI24 = 0.74159265358979323846264338328;
                           double K              = 0.0;
			   const double A2       = A*A;
			   const double gamm2    = gamma*gamma*sig0;
			   const double t0       = PI4r8*(A2/gamm2);
			   K = std::log(t0)/PI24;
			   return (K);
		    }


		    // Helper function wave-number calculation
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float wavenumber_r4_1(const float gamma) { //m,wavelength

		            return (PI2r4/gamma);
		    }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double wavenumber_r8_1(const double gamma) { //m,wavelength

		            return (PI2r8/gamma);
		    }


		    /*
                        
                      */

	             __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void effective_rcs_r4_1(const float gamma, // m, wavelength
		                             const float R,     // nm, target range
					     const float h_a,   // ft, antenna height
					     const float h_t,   // ft, target height
					     float * __restrict sig_eff,    // m, RCS effective
					     float * __restrict sig_eff_db) { // db, RCS effective

                         const float t0 = PI2r4*h_a*h_t;
			 const float t1 = R*gamma;
			 const float F  = 2.0f*cephes_sin(t0/t1);
			 const float F4 = F*F*F*F;
			 *sig_eff       = sig0_r4*F4;
			 *sig_eff_db    = 10.0f*cephes_log10f(*sig_eff);
		   }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void effective_rcs_r8_1(const double gamma, // m, wavelength
		                             const double R,     // nm, target range
					     const double h_a,   // ft, antenna height
					     const double h_t,   // ft, target height
					     double * __restrict sig_eff,    // m, RCS effective
					     double * __restrict sig_eff_db) { // db, RCS effective

                         const double t0 = PI2r8*h_a*h_t;
			 const double t1 = R*gamma;
			 const double F  = 2.0*std::sin(t0/t1);
			 const double F4 = F*F*F*F;
			 *sig_eff        = sig0_r8*F4;
			 *sig_eff_db     = 10.0*std::log10(*sig_eff);
		   }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
                    void
		    bistatic_target_rcs_r4_1(const float sig0, // m^2, RCS monostatic target
		                             const float K,    // empirical constant defined by target configuration and complexity  // use function emprirical_K for computation
					     const float Beta, // deg, bistatic angle (for input=1)
					     const float R1,   // m, transmitter - target range (input=2)
					     const float R2,   // m, receiver    - target range (input=2)
					     const float B,    // m, baseline
					     const int32_t type, // input switch (1,2)
					     float * __restrict sigma, // RCS                  
                                             float * __restrict sigma_db) {   // dBsm of Target
                         
                         if(type==1) {
                            const float alpha = Beta*zr4;
			    const float t0    = K*std::abs(alpha)-2.4f*K-1.0f;
			    const float t1    = 1.0f+cephes_expf(t0);
			    *sigma            = sig0*t1;
			    *sigma_db         = 10.0*cephes_log10f(*sigma+0.00000000001f);
			 }
			 else if(type==2) {
                            const float t0    = 1.0f/(2.0*R1*R2);
			    const float R12   = R1*R1;
			    const float R22   = R2*R2;
			    const float B2    = B*B;
			    const float t1    = R12+R22+B2;
			    const float gam   = cephes_acosf(t0*t1);
			    const float alpha = gam;
			    const float t2    = K*std::abs(alpha)-2.4f*K-1.0f;
			    const float t3    = 1.0f+cephes_expf(t0);
			    *sigma            = sig0*t1;
			    *sigma_db         = 10.0*cephes_log10f(*sigma+0.00000000001f);
			 }
		    }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
                    void
		    bistatic_target_rcs_r8_1(const double sig0, // m^2, RCS monostatic target
		                             const double K,    // empirical constant defined by target configuration and complexity  // use function emprirical_K for computation
					     const double Beta, // deg, bistatic angle (for input=1)
					     const double R1,   // m, transmitter - target range (input=2)
					     const double R2,   // m, receiver    - target range (input=2)
					     const double B,    // m, baseline
					     const int32_t type, // input switch (1,2)
					     double * __restrict sigma, // RCS                  
                                             double * __restrict sigma_db) {   // dBsm of Target
                         
                         if(type==1) {
                            const double alpha = Beta*zr8;
			    const double t0    = K*std::abs(alpha)-2.4*K-1.0;
			    const double t1    = 1.0+std::exp(t0);
			    *sigma            = sig0*t1;
			    *sigma_db         = 10.0*std::log10(*sigma+0.00000000001);
			 }
			 else if(type==2) {
                            const double t0    = 1.0/(2.0*R1*R2);
			    const double R12   = R1*R1;
			    const double R22   = R2*R2;
			    const double B2    = B*B;
			    const double t1    = R12+R22+B2;
			    const double gam   = std::acos(t0*t1);
			    const double alpha = gam;
			    const double t2    = K*std::abs(alpha)-2.4*K-1.0;
			    const double t3    = 1.0+std::exp(t0);
			    *sigma             = sig0*t1;
			    *sigma_db          = 10.0*std::log10(*sigma+0.00000000001);
			 }
		    }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
                    void antenna_rcs_r4_1(const float Ae,   // m^2, antenna efective aperture
		                          const float gam,  // m, wavelength
					  const float G,    // refelcetion coefficient
					  float * __restrict sigma, // m^2, radar cross section
					  float * __restrict sigma_db) { //dBsm, rcs with repsect to 1m^2

			  const float Ae2  = Ae*Ae;
			  const float gam2 = gam*gam;
			  const float t0   = Ae2/gam2;
			  *sigma           = PI4r4*t0*G;
			  *sigma_db        = 10.0f*cephes_log10f(*sigma);
		   }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
                    void antenna_rcs_r8_1(const double Ae,   // m^2, antenna efective aperture
		                          const double gam,  // m, wavelength
					  const double G,    // reflection coefficient
					  double * __restrict sigma, // m^2, radar cross section
					  double * __restrict sigma_db) { //dBsm, rcs with repsect to 1m^2

			  const double Ae2  = Ae*Ae;
			  const double gam2 = gam*gam;
			  const double t0   = Ae2/gam2;
			  *sigma            = PI4r8*t0*G;
			  *sigma_db         = 10.0*std::log10(*sigma);
		   }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
		    float bird_insect_rcs_r4_1(const float W) { //gram, weight of bird or insect

                           return (-46.0f+5.8f*cephes_log10f(W));
		    }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
		    double bird_insect_rcs_r8_1(const double W) { //gram, weight of bird or insect

                           return (-46.0+5.8*std::log10(W));
		    }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
		    float cone_ogive_rcs_r4_1(const float gamma, // m, wavelength
		                        const float cha) {   // deg, cone half-angle, incidence nose-on
                         // This model applies to infintely large cone, or ogive larger than incident wavelength
			 float sig1           = 0.0f; 
                         const float gamm2    = gamma*gamma;
			 constexpr float PI16 = 50.2654824574366918154023f;
			 const float ratio    = gamm2/PI16;
			 //const float targ     = cephes_tanf(zr4*cha);
			 const float targ     = zr4*cha;
			 const float targ4    = targ*targ*targ*targ;
			 const float tan4     = cephes_tanf(targ4);
			 sig1                 = ratio*tan4;
			 return (sig1);
		    }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
		    double cone_ogive_rcs_r8_1(const double gamma, // m, wavelength
		                               const double cha) {   // deg, cone half-angle, incidence nose-on
                         // This model applies to infintely large cone, or ogive larger than incident wavelength
			 double sig1           = 0.0; //rcs
                         const double gamm2    = gamma*gamma;
			 constexpr double PI16 = 50.2654824574366918154023;
			 const double ratio    = gamm2/PI16;
			 const double targ     = zr4*cha;
			 const double targ4    = targ*targ*targ*targ;
			 const double tan4     = std::tan(targ4);
			 sig1                  = ratio*tan4;
			 return (sig1);
		    }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
		    float cylinder_rcs_r4_1(const float Rcyl,  //m, radius of cylinder
		                            const float Lcyl,  //m, length of cylinder
					    const float gamm,  //m, wavelength
                                            const float theta) { //deg, angle between cylinder axis and radar line-of-sight

			   float sig        = 0.0f; //m^2, rcs
			   const float t0   = PI4r4*Rcyl;
			   const float sinz = cephes_sinf(zr4*theta);
		           const float term = (t0*sinz)/gamm;
			   if(term<1.0f) {
                              return (-1.0f);
			   }
			   const float cosz  = cephes_cosf(zr4*theta);
			   const float term1 = (Rcyl*gamm*sinz)/PI2r4;
			   const float sarg1 = (PI2r4*Lcyl)/gamm;
			   const float sarg2 = sarg1*cosz;
			   const float term2 = cephes_sinf(sarg2/cosz);
			   const float term22= term2*term2;
			   sig               = term1*term22;
			   return (sig);
                   }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
		    double cylinder_rcs_r8_1(const double Rcyl,  //m, radius of cylinder
		                             const double Lcyl,  //m, length of cylinder
					     const double gamm,  //m, wavelength
                                             const double theta) { //deg, angle between cylinder axis and radar line-of-sight

			   double sig        = 0.0; //m^2, rcs
			   const double t0   = PI4r8*Rcyl;
			   const double sinz = std::sin(zr8*theta);
		           const double term = (t0*sinz)/gamm;
			   if(term<1.0) {
                              return (-1.0);
			   }
			   const double cosz  = std::cos(zr4*theta);
			   const double term1 = (Rcyl*gamm*sinz)/PI2r8;
			   const double sarg1 = (PI2r4*Lcyl)/gamm;
			   const double sarg2 = sarg1*cosz;
			   const double term2 = std::sin(sarg2/cosz);
			   const double term22= term2*term2;
			   sig                = term1*term22;
			   return (sig);
                   }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
		    float disk_rcs_r4_1(const float Rd,    //m^2, radius of disk
		                        const float theta, //deg, angle relative to disk normal
					const float gamm,  //m, wavelength
					const int32_t type) { // disk type: 1=electrically large disk, 2=electrically small circular disk

			   float sig         = 0.0f;
			   if(type==1) {
			      const float k     = wavenumber_r4_1(gamm);
			      const float gamm2 = gamm*gamm;
			      const float carg  = zr4*theta;
			      const float coszt = cephes_cosf(carg*carg);
                              const float t0    = PIr4*Rd*Rd;
			      const float t1    = PI4r4*t0*t0;
			      const float tleft = t1/gamm2;
			      const float sinzt = cephes_sinf(zr4*theta);
			      const float jarg  = 2.0f*k*Rd*sinzt;
			      const float bj1   = cephes_j1f(jarg);
			      const float term  = 2.0f*bj1/jarg;
			      const float tright= term*term;
			      sig               = tleft*tright*coszt;
			      return (sig);
			   }
			   else if(type==2) {
                              const float Rd2   = Rd*Rd;
			      const float tleft = 3534.2917352885173932704738f*Rd2;
			      const float ratio = Rd/gamm;
			      const float ratio4= ratio*ratio*ratio*ratio;
			      sig               = tleft*ratio4;
			      return (sig);
                               
			   }
			   
                   }

#include <math.h> // bessel function j1
		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
		    double disk_rcs_r8_1(const double Rd,    //m^2, radius of disk
		                         const double theta, //deg, angle relative to disk normal
					 const double gamm,  //m, wavelength
					 const int32_t type) { // disk type: 1=electrically large disk, 2=electrically small circular disk

			   double sig         = 0.0;
			   if(type==1) {
			      const double k     = wavenumber_r8_1(gamm);
			      const double gamm2 = gamm*gamm;
			      const double carg  = zr8*theta;
			      const double coszt = std::cos(carg*carg);
                              const double t0    = PIr8*Rd*Rd;
			      const double t1    = PI4r8*t0*t0;
			      const double tleft = t1/gamm2;
			      const double sinzt = std::sin(zr4*theta);
			      const double jarg  = 2.0*k*Rd*sinzt;
			      const double bj1   = j1(jarg);
			      const double term  = 2.0*bj1/jarg;
			      const double tright= term*term;
			      sig               = tleft*tright*coszt;
			      return (sig);
			   }
			   else if(type==2) {
                              const double Rd2   = Rd*Rd;
			      const double tleft = 3534.2917352885173932704738*Rd2;
			      const double ratio = Rd/gamm;
			      const double ratio4= ratio*ratio*ratio*ratio;
			      sig                = tleft*ratio4;
			      return (sig);
                               
			   }
			   
                   }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
		    float curved_edge_rcs_r4_1(const float Redg, //m, radius of edge contour, incidence edge-perpendicular
		                               const float gamm) { //m, wavelength

                           float sig   = 0.0f; //m^2, rcs
			   sig         = 0.5f*(Redg*gamm);
			   return (sig);
		   }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
		    float curved_edge_rcs_r8_1(const double Redg, //m, radius of edge contour, incidence edge-perpendicular
		                               const double gamm) { //m, wavelength

                           double sig   = 0.0; //m^2, rcs
			   sig         = 0.5*(Redg*gamm);
			   return (sig);
		   }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
		    float straight_edge_rcs_r4_1(const float Ledg) { //m, length of straight edge

		            float sig   = 0.0f;
			    sig         = (Ledg*Ledg)/PIr4;
			    return (sig);
                   }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
		    double straight_edge_rcs_r8_1(const double Ledg) { //m, length of straight edge

		            double sig   = 0.0;
			    sig          = (Ledg*Ledg)/PIr8;
			    return (sig);
                   }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
		    float ellipsoid_rcs_r4_1(const float a, //m, semimajor axis
		                             const float b) { //m, semiminor axis

                            float sig      = 0.0f;
			    const float b4 = b*b*b*b;
			    const float a2 = a*a;
			    sig            = (PIr4*b4)/a2;
			    return (sig);
		   }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
		    double ellipsoid_rcs_r4_1(const double a, //m, semimajor axis, incidence along the major axis
		                              const double b) { //m, semiminor axis

                            double sig      = 0.0;
			    const double b4 = b*b*b*b;
			    const double a2 = a*a;
			    sig             = (PIr8*b4)/a2;
			    return (sig);
		   }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
		    float plate_rcs_r4_1(const float x, //m, plate length
		                         const float y, //m, plate width
					 const float theta, //deg, angle relative to plate normal
					 const float phi,    //deg, angle between the plate normal and radar los
					 const float gamm) { //m, wavelength

			    float sig         = 0.0f;
			    const float gamm2 = gamm*gamm;
			    const float sinzt = cephes_sinf(zr4*theta);
			    const float xy2   = PI4r4*(x*y*x*y);
			    const float coszt = cephes_cosf(zr4*phi);
			    const float zth2  = zr4*theta*zr4*theta;
			    const float u     = wavenumber_r4_1(gamm)*x*sinzt*coszt;
			    const float sinu  = cephes_sinf(u)/u;
			    const float v     = wavenumber_r4_1(gamm)*y*sinzt*coszt;
			    const float sinv  = cephes_sinf(v)/v;
			    const float tleft = xy2/gamm2;
			    const float tmid  = sinu*sinv*sinu*sinv;
			    sig               = tleft*tmid*coszt;
			    return (sig);
		   }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
		    double plate_rcs_r8_1(const double x, //m, plate length
		                          const double y, //m, plate width
					  const double theta, //deg, angle relative to plate normal
					  const double phi,    //deg, angle between the plate normal and radar los
					  const double gamm) { //m, wavelength

			    double sig         = 0.0;
			    const double gamm2 = gamm*gamm;
			    const double sinzt = std::sin(zr8*theta);
			    const double xy2   = PI4r8*(x*y*x*y);
			    const double coszt = std::cos(zr8*phi);
			    const double zth2  = zr8*theta*zr8*theta;
			    const double u     = wavenumber_r8_1(gamm)*x*sinzt*coszt;
			    const double sinu  = std::sin(u)/u;
			    const double v     = wavenumber_r8_1(gamm)*y*sinzt*coszt;
			    const double sinv  = std::sin(v)/v;
			    const double tleft = xy2/gamm2;
			    const double tmid  = sinu*sinv*sinu*sinv;
			    sig                = tleft*tmid*coszt;
			    return (sig);
		   }


     }//radiolocation

} //gms











#endif /*__GMS_RCS_HPP__*/
