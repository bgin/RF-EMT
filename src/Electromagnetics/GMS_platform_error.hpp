
#ifndef __GMS_PLATFORM_ERROR_HPP__
#define __GMS_PLATFORM_ERROR_HPP__ 270120221613
 

namespace file_info {

     const unsigned int GMS_PLATFORM_ERROR_MAJOR = 1;
     const unsigned int GMS_PLATFORM_ERROR_MINOR = 1;
     const unsigned int GMS_PLATFORM_ERROR_MICRO = 0;
     const unsigned int GMS_PLATFORM_ERROR_FULLVER =
       1000U*GMS_PLATFORM_ERROR_MAJOR+100U*GMS_PLATFORM_ERROR_MINOR+
       10U*GMS_PLATFORM_ERROR_MICRO;
     const char * const GMS_PLATFORM_ERROR_CREATION_DATE = "27-01-2022 16:13 +00200 (THR 27 JAN 2022 GMT+2)";
     const char * const GMS_PLATFORM_ERROR_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_PLATFORM_ERROR_SYNOPSIS      = "Radar Jamming Equations."

}


#include <cstdint>
#include <cmath> //for double precision
#include "GMS_cephes.h" // single precision
#include "GMS_config.h"
#include "GMS_radar_types.h"

namespace gms {

          namespace radiolocation {


	              // L1D cache data prefetch
		     __ATTR_ALWAYS_INLINE__
		     static
		     inline
		     bool data_prefetch_r4_1(const PlatformErrAoS_R4_1 &pe) {

		            _mm_prefetch((const char*)&pe,0);
			    return (true);
		     }

		     __ATTR_ALWAYS_INLINE__
		     static
		     inline
		     bool data_prefetch_r8_1(const PlatformErrAoS_R8_1 &pe) {

		            _mm_prefetch((const char*)&pe,0);
			    return (true);
		     }


	            // Initialize PlatformErrorAoS_R4_1 data type
		     __ATTR_ALWAYS_INLINE__
		     static
		     inline
		     void initPlatformErrAoS_R4_1(const float xR,
		                                  const float xpsi,
						  const float xtheta,
						  const float xda1,
						  const float xda2,
						  const float xda3,
						  const float xdx1,
						  const float xdx2,
						  const float xdx3,
						  PlatformErrAoS_R4_1 &pe) {

                pe.R    = xR;
			    pe.psi  = xpsi;
			    pe.theta= xtheta;
			    pe.da1  = xda1;
			    pe.da2  = xda2;
			    pe.da3  = xda3;
			    pe.dx1  = xdx1;
			    pe.dx2  = xdx2;
			    pe.dx3  = xdx3;
		    }


		     // Initialize PlatformErrorAoS_R8_1 data type
		     __ATTR_ALWAYS_INLINE__
		     static
		     inline
		     void initPlatformErrAoS_R8_1(const double xR,
		                                  const double xpsi,
						  const double xtheta,
						  const double xda1,
						  const double xda2,
						  const double xda3,
						  const double xdx1,
						  const double xdx2,
						  const double xdx3,
						  PlatformErrAoS_R8_1 &pe) {

                pe.R    = xR;
			    pe.psi  = xpsi;
			    pe.theta= xtheta;
			    pe.da1  = xda1;
			    pe.da2  = xda2;
			    pe.da3  = xda3;
			    pe.dx1  = xdx1;
			    pe.dx2  = xdx2;
			    pe.dx3  = xdx3;
		    }



		    

		     
	             // Compute Radar platform azimut and elevation
		     // measurement errors
	         __ATTR_ALWAYS_INLINE__
		     static
		     inline
		     void platform_orient_errors_r4_1(const PlatformErrAoS_R4_1 &pe,
		                                      float &dpsi1,      //ang,min, azimuth measurement error
					              float &dth1) {  //ang,min, elevation measurement error

			   constexpr float z  = 0.0174532925199432957692f;
			   constexpr float PI = 3.1415926535897932384626f;
			   constexpr float PI2= 6.2831853071795864769253F;
			   const float xpsi  = pe.psi;
			   const float xth   = pe.theta;
			   const float czpsi = cephes_cosf(z*xpsi);
			   const float szpsi = cephes_sinf(z*xpsi);
			   const float thc   = 90.0f-xth;
                           const float xda1  = pe.da1;
			   const float czth  = cephes_cosf(z*xth);
			   const float szth  = cephes_sinf(z*xth);
			   const float arg1  = z*xda1/60.0f;
			   const float c1    = cephes_cosf(arg1);
			   const float s1    = cephes_sinf(arg1);
			   const float xda2  = pe.da2;
			   const float arg2  = z*xda2/60.0f;
			   const float c2    = cephes_cosf(arg2);
			   const float s2    = cephes_sinf(arg2);
			   const float xda3  = pe.da3;
			   const float arg3  = z*xda3/60.0f;
			   const float c3    = cephes_cosf(arg3);
			   const float s3    = cephes_sinf(arg3);
			   const float A0    = czpsi*szth*c1*c2;
			   const float B0    = -czpsi*szth*(s1*c3+c1*s2*s3);
			   const float C0    = czpsi*szth*(s1*s3-c1*s2*c3);
			   const float A1    = szpsi*szth*s1*c2;
			   const float B1    = -szpsi*szth*(s1*s2*s3-c1*c3);
			   const float C1    = -szpsi*szth*(c1*s3+s1*s2*c3);
			   const float A2    = czth*s2;
			   const float B2    = czth*c2*s3;
			   const float C2    = czth*c2*c3;
			   const float A012  = A0+A1+A2;
			   const float Asqr  = A012*A012;
			   const float B012  = B0+B1+B2;
			   const float Bsqr  = B012*B012;
			   const float C012  = C0+C1+C2;
			   const float Csqr  = C012*C012;
			   const float U     = A012/(Asqr+Bsqr);
			   if(zpsi>=0.0f && z<=PI) {
			      dpsi1 = (cephes_acosf(U)-zpsi)/z*60.0f;
			   }
			   else {
                              dpsi1 = (PI2-cephes_acosf(U)-zpsi)/z*60.0f; 
			   }
			   const float V     = C012/(Asqr+Bsqr+Csqr);
			   dthi = -((cephes_acosf(V)-z*xth)/z*60.0f);
			   
		   }


		     __ATTR_ALWAYS_INLINE__
		     static
		     inline
		     void platform_orient_errors_r8_1(const PlatformErrAoS_R8_1 &pe,
		                                      double &dpsi1,      //ang,min, azimuth measurement error
					              double &dth1) {  //ang,min, elevation measurement error

			   constexpr double z  = 0.0174532925199432957692;
			   constexpr double PI = 3.1415926535897932384626;
			   constexpr double PI2= 6.2831853071795864769253;
			   const double xpsi  = pe.psi;
			   const double xth   = pe.theta;
			   const double czpsi = std::cos(z*xpsi);
			   const double szpsi = std::sin(z*xpsi);
			   const double thc   = 90.0-xth;
               const double xda1  = pe.da1;
			   const double czth  = std::cos(z*xth);
			   const double szth  = std::sin(z*xth);
			   const double arg1  = z*xda1/60.0;
			   const double c1    = std::cos(arg1);
			   const double s1    = std::sin(arg1);
			   const double xda2  = pe.da2;
			   const double arg2  = z*xda2/60.0;
			   const double c2    = std::cos(arg2);
			   const double s2    = std::sin(arg2);
			   const double xda3  = pe.da3;
			   const double arg3  = z*xda3/60.0;
			   const double c3    = std::cos(arg3);
			   const double s3    = std::sin(arg3);
			   const double A0    = czpsi*szth*c1*c2;
			   const double B0    = -czpsi*szth*(s1*c3+c1*s2*s3);
			   const double C0    = czpsi*szth*(s1*s3-c1*s2*c3);
			   const double A1    = szpsi*szth*s1*c2;
			   const double B1    = -szpsi*szth*(s1*s2*s3-c1*c3);
			   const double C1    = -szpsi*szth*(c1*s3+s1*s2*c3);
			   const double A2    = czth*s2;
			   const double B2    = czth*c2*s3;
			   const double C2    = czth*c2*c3;
			   const double A012  = A0+A1+A2;
			   const double Asqr  = A012*A012;
			   const double B012  = B0+B1+B2;
			   const double Bsqr  = B012*B012;
			   const double C012  = C0+C1+C2;
			   const double Csqr  = C012*C012;
			   const double U     = A012/(Asqr+Bsqr);
			   if(zpsi>=0.0f && z<=PI) {
			      dpsi1 = (std::acos(U)-zpsi)/z*60.0;
			   }
			   else {
                    dpsi1 = (PI2-std::acos(U)-zpsi)/z*60.0; 
			   }
			   const float V     = C012/(Asqr+Bsqr+Csqr);
			   dthi = -((std::acos(V)-z*xth)/z*60.0);
			   
		   }


		   // Platform position errors caused by the angular and range
		   // measurements.
		     __ATTR_ALWAYS_INLINE__
		     static
		     inline
		     void platform_pos_errors_r4_1(PlatformErrAoS_R4_1 &pe,
		                                   float &dpsi2, //min, azimuth error
                                                   float &dth2,  //min, elevation error
						   float &dR) {  //m,   range measurement error

			  constexpr float tb  = 0.0000000001f;
			  constexpr float z   = 0.0174532925199432957692f;
			  constexpr float PI = 3.1415926535897932384626f;
			  constexpr float PI2= 6.2831853071795864769253f;
			  const float xR    = pe.R;
			  const float xpsi  = pe.psi;
			  const float zpsi  = z*xpsi;
			  const float xth   = pe.theta;
			  const float szth  = cephes_sinf(z*xth);
			  const float czth  = cephes_cosf(z*xth);
                          const float xdx1  = pe.dx1;
			  const float xdx12 = xdx1*xdx1;
			  const float xdx2  = pe.dx2;
			  const float xdx22 = xdx2*xdx2;
			  const float xdx3  = pe.dx3;
			  const float xdx32 = xdx3*xdx3;
			  const float eta1  = cephes_sqrtf(xdx12+xdx22+xdx32);
			  const float eta1R = eta1/xR;
			  const float eta1R2= eta1R*eta1R;
			  const float eta2  = cephes_atanf(xdx2/(xdx1+tb));
			  const float num   = cephes_sqrtf(xdx12+xdx22);
			  const float eta3  = cephes_atanf(num/xdx3+tb);
			  const float seta3 = cephes_sinf(eta3);
			  const float ceta3 = cephes_cosf(eta3);
			  const float lterm = cephes_cosf(zpsi)*cephes_sinf(zpsi);
			  const float rterm = eta1R*cephes_cosf(eta2)*seta3;
			  const float D0    = lterm-rterm;
			  const float lterm2= cephes_sinf(zpsi)*szth;
			  const float rterm2= eta1R*cephes_sinf(eta2)*seta3;
			  const float D1    = lterm2-rterm2;
			  const float D2    = czth-eta1R*ceta3;
			  const float Dk01  = D0*D0+D1*D1;
			  const float U     = D0/cephes_sqrtf(Dk01);
			  if(zpsi>=0.0f && zpsi<=PI) {
                             dpsi2 = (cephes_acosf(U)-zpsi)/z*60.0f;
			  }
			  else {
                             dpsi2 = (cephes_acosf(U)-zpsi)/z*60.0f;
			  }
			  const float Dk012 = DK01+D2*D2;
			  const float V     = D2/cephes_sqrtf(Dk012);
			  dth2 = -((cephes_acosf(V)-z*xth)/z*60.0f);
			  const float lterm3= szth*cephes_cosf(zpsi-eta2)*seta3;
			  const float rterm3= czth*ceta3;
			  const float aR    = lterm3+rterm3;
			  const float lterm4= 1.0f-2.0f*eta1R*aR+eta1R2;
			  const float delR  = cephes_sqrtf(lterm4);
			  dR                = xR*(delR-1.0f);
			  
		    }


		     __ATTR_ALWAYS_INLINE__
		     static
		     inline
		     void platform_pos_errors_r8_1(PlatformErrAoS_R8_1 &pe,
		                                   double &dpsi2, //min, azimuth error
                                                   double &dth2,  //min, elevation error
						   double &dR) {  //m,   range measurement error

			  constexpr double tb  = 0.0000000001;
			  constexpr double z   = 0.0174532925199432957692;
			  constexpr double PI = 3.1415926535897932384626;
			  constexpr double PI2= 6.2831853071795864769253;
			  const double xR    = pe.R;
			  const double xpsi  = pe.psi;
			  const double zpsi  = z*xpsi;
			  const double xth   = pe.theta;
			  const double szth  = std::sin(z*xth);
			  const double czth  = std::cos(z*xth);
                          const double xdx1  = pe.dx1;
			  const double xdx12 = xdx1*xdx1;
			  const double xdx2  = pe.dx2;
			  const double xdx22 = xdx2*xdx2;
			  const double xdx3  = pe.dx3;
			  const double xdx32 = xdx3*xdx3;
			  const double eta1  = std::sqrt(xdx12+xdx22+xdx32);
			  const double eta1R = eta1/xR;
			  const double eta1R2= eta1R*eta1R;
			  const double eta2  = std::atan(xdx2/(xdx1+tb));
			  const double num   = std::sqrt(xdx12+xdx22);
			  const doubl eta3   = std::atan(num/xdx3+tb);
			  const double seta3 = std::sin(eta3);
			  const double ceta3 = std::cos(eta3);
			  const double lterm = std::cos(zpsi)*std::sin(zpsi);
			  const double rterm = eta1R*std::cos(eta2)*seta3;
			  const double D0    = lterm-rterm;
			  const double lterm2= cephes_sinf(zpsi)*szth;
			  const double rterm2= eta1R*std::sin(eta2)*seta3;
			  const double D1    = lterm2-rterm2;
			  const double D2    = czth-eta1R*ceta3;
			  const double Dk01  = D0*D0+D1*D1;
			  const double U     = D0/std::sqrt(Dk01);
			  if(zpsi>=0.0 && zpsi<=PI) {
                             dpsi2 = (std::acos(U)-zpsi)/z*60.0;
			  }
			  else {
                             dpsi2 = (std::acos(U)-zpsi)/z*60.0;
			  }
			  const double Dk012 = DK01+D2*D2;
			  const double V     = D2/std::sqrt(Dk012);
			  dth2 = -((std::acos(V)-z*xth)/z*60.0);
			  const double lterm3= szth*std::cos(zpsi-eta2)*seta3;
			  const double rterm3= czth*ceta3;
			  const double aR    = lterm3+rterm3;
			  const double lterm4= 1.0f-2.0f*eta1R*aR+eta1R2;
			  const double delR  = std::sqrt(lterm4);
			  dR                 = xR*(delR-1.0f);
			  
		    }


		    // Monopulse antenna pattern as a voltage sum pattern
		    // of the angle from the (beamwidth) normalized axis pattern
		     __ATTR_ALWAYS_INLINE__
		     static
		     inline
		     float v_sum_pattern_r4_1(const float u) {

		             const float uu = u*u;
                             return (cephes_expf(-1.3866f*uu));
			     
		     }

		     __ATTR_ALWAYS_INLINE__
		     static
		     inline
		     double v_sum_pattern_r8_1(const double u) {

		             const double uu = u*u;
                             return (std::exp(-1.3866*uu));
			     
		     }


		    // Monopulse antenna pattern as a voltage differential pattern
		    // of the angle from the (beamwidth) normalized axis pattern
		     __ATTR_ALWAYS_INLINE__		    
		     static
		     inline
		     float v_diff_pattern_r4_1(const float u) {

                            return (2.354f*u*v_sum_pattern_r4_1(u));
		    }


		     __ATTR_ALWAYS_INLINE__		   
		     static
		     inline
		     double v_diff_pattern_r8_1(const double u) {

                            return (2.354*u*v_sum_pattern_r8_1(u));
		    }

                    // Calculation of the elevation multipath error.
		     __ATTR_ALWAYS_INLINE__
		     static
		     inline
		     float elev_multipath_err_r4_1(const PropagationErrAoS_R4_1 &pe,
		                                   float &Sigem) { // deg, elevation multipath error

			    constexpr float ae   = 8493333.0f;
			    constexpr float PI   = 3.1415926535897932384626f;
			    constexpr float sqt2 = 1.4142135623730950488017f;
			    const float   rad1 = 57.2957795130823208767982f;
			    const float   ae2  = ae+ae;
                            const float   xR   = pe.R;
			    const float   xha  = pe.ha;
			    const float   xht  = pe.ht;
			    const float   xthm = pe.thmax;
			    const float   xth3 = pe.th3;
			    const float   xps0 = pe.psi0;
			    const float   xpss = pe.psis;
			    const float   xpsv = pe.psiv;
			    const float   xkme = pe.kme;
			    const float   RR   = xR*xR;
			    const float   dp   = (4.0f*ae*(xha+xht)+RR)*0.3333333333333333333f;
			    const float   p    = cephes_sqrtf(dp);
			    const float   p3   = p*p*p;
			    const float   dph  = (ae2*(xht-xha)*xR)/p3;
			    const float   ph   = cephes_acosf(dph);
			    const float   cph  = cephes_cosf((ph+PI)*0.3333333333333333333f);
			    const float   d1   = xR*0.5f-p*cph;
			    const float   H1   = xha-(d1*d1)/ae2;
			    const float   H2   = ((xR-d1)*H1)/d1;
			    const float   tht  = rad1*cephes_asinf((H2-H1)/xR));
			    const float   ps0  = rad1*cephes_asinf((H2+H1)/xR));
			    const float   thd  = -xthm+tht;
			    const float   fs   = v_sum_pattern_r4_1(thd/xth3);
			    const float   fs2  = fs*fs;
			    const float   thr  = -xthm-ps0;
			    const float   fd   = v_diff_pattern_r4_1(thr/xth3);
			    const float   Fdel = fd*xps0*xpss*xpsv;
			    const float   Fdel2= Fdel*Fdel;
			    const float   Strm = th3/(sqt2*xkme);
			    const float   S_es = Strm*cephes_sqrtf(Fdel2/fs2);
			    const float   S_esm= th3*0.5f;
			    Sigem              = (S_es<S_esm)?S_es:S_esm;
			    
		    }


		     __ATTR_ALWAYS_INLINE__
		     static
		     inline
		     double elev_multipath_err_r8_1(const PropagationErrAoS_R8_1 &pe,
		                                    double &Sigem) { // deg, elevation multipath error

			    constexpr double ae   = 8493333.0;
			    constexpr double PI   = 3.1415926535897932384626;
			    constexpr double sqt2 = 1.4142135623730950488017;
			    const double   rad1 = 57.2957795130823208767982;
			    const double   ae2  = ae+ae;
                            const double   xR   = pe.R;
			    const double   xha  = pe.ha;
			    const double   xht  = pe.ht;
			    const double   xthm = pe.thmax;
			    const double   xth3 = pe.th3;
			    const double   xps0 = pe.psi0;
			    const double   xpss = pe.psis;
			    const double   xpsv = pe.psiv;
			    const double   xkme = pe.kme;
			    const double   RR   = xR*xR;
			    const double   dp   = (4.0*ae*(xha+xht)+RR)*0.3333333333333333333;
			    const double   p    = std::sqrt(dp);
			    const double   p3   = p*p*p;
			    const double   dph  = (ae2*(xht-xha)*xR)/p3;
			    const double   ph   = std::acos(dph);
			    const double   cph  = std::cos((ph+PI)*0.3333333333333333333);
			    const double   d1   = xR*0.5f-p*cph;
			    const double   H1   = xha-(d1*d1)/ae2;
			    const double   H2   = ((xR-d1)*H1)/d1;
			    const double   tht  = rad1*std::asin((H2-H1)/xR));
			    const double   ps0  = rad1*std::asin((H2+H1)/xR));
			    const double   thd  = -xthm+tht;
			    const double   fs   = v_sum_pattern_r8_1(thd/xth3);
			    const double   fs2  = fs*fs;
			    const double   thr  = -xthm-ps0;
			    const double   fd   = v_diff_pattern_r8_1(thr/xth3);
			    const double   Fdel = fd*xps0*xpss*xpsv;
			    const double   Fdel2= Fdel*Fdel;
			    const double   Strm = th3/(sqt2*xkme);
			    const double   S_es = Strm*std::sqrt(Fdel2/fs2);
			    const double   S_esm= th3*0.5;
			    Sigem              = (S_es<S_esm)?S_es:S_esm;
			    
		    }


		    // Elevation error due to atmospheric refraction
		     __ATTR_ALWAYS_INLINE__
		     static
		     inline
		     float elev_refract_error_r4_1(const PropagationErrAoS_R4_1 &pe,
		                                   float &Siger) {
						   //deg, rms error of elevation measurement
                            float C1             = 0.0f;
			    float a              = 0.0f;
			    float b              = 0.0f;
			    float dth            = 0.0f;
			    constexpr float ae   = 8493333.0f;
			    constexpr float PI   = 3.1415926535897932384626f;
			    constexpr float sqt2 = 1.4142135623730950488017f;
			    const float   rad1 = 57.2957795130823208767982f;
			    const float   ae2  = ae+ae;
                            const float   xR   = pe.R;
			    const float   xha  = pe.ha;
			    const float   xht  = pe.ht;
			    const float   xthm = pe.thmax;
			    const float   xth3 = pe.th3;
			    const float   xps0 = pe.psi0;
			    const float   xpss = pe.psis;
			    const float   xpsv = pe.psiv;
			    const float   xkme = pe.kme;
			    const float   xNs  = pe.Ns;
			    const float   RR   = xR*xR;
			    const float   dp   = (4.0f*ae*(xha+xht)+RR)*0.3333333333333333333f;
			    switch(pe.flucts) {
                                 case 1: C1=0.000001f;
				 break;
				 case 2: C1=2.0f*0.000001f;
				 break;
				 case 3: C1=4.0f*0.000001f;
				 break;
			    } 
			    const float   p    = cephes_sqrtf(dp);
			    const float   p3   = p*p*p;
			    const float   dph  = (ae2*(xht-xha)*xR)/p3;
			    const float   ph   = cephes_acosf(dph);
			    const float   cph  = cephes_cosf((ph+PI)*0.3333333333333333333f);
			    const float   d1   = xR*0.5f-p*cph;
			    const float   H1   = xha-(d1*d1)/ae2;
			    const float   H2   = ((xR-d1)*H1)/d1;
			    const float   tht  = rad1*cephes_asinf((H2-H1)/xR));
			    const float   ztht = rad1*tht;
			    const float   thef = ztht+(0,00025f/(ztht+0.028f));
			    const float   sthf = cephes_sinf(thef);
			    const float   sig0 = C1/rad1*cephes_sqrtf(3/sthf);
			    const float   cthf = cephes_cotf(tht);
			    const float   lgth = cephes_logf(tht);
			    if(tht<0.1f) {
                               b = 0.0000001f*cthf*(14.5+4.5f*lgth);
			    }
			    else {
                               b = 0.000001g*cthf;
			    }
			    if(tht<0.02f) {
                               a = −0.00015f*cthf*(1+0.3*lgth);
			    }
			    else if(tht>=0.02f && tht<=0.2f) {
                               a = 0.00005f*cthf*(1.0f+1.4f*lgth);
			    }
			    else if(tht>0.2f) {
                               a = 0.0f;
			    }
			    const float rtrm = 1.0f+(0.001f/((0.01f+tht))*(0.01f+tht));
			    if(tht>=0.0f && tht<=0.001f) {
                               dth = 0.019f;
			    }
			    else {
                               b*xNs+a*(xht/(11000.0+xht);
			    }
			    const float sig1 = dth/0.3490658503988659153847f*rtrm;
			    const float sig12 = sig1*sig1;
			    Siger = sig12+sig12;
		    }


		      // Elevation error due to atmospheric refraction
		     __ATTR_ALWAYS_INLINE__
		     static
		     inline
		     double elev_refract_error_r8_1(const PropagationErrAoS_R8_1 &pe,
		                                    double &Siger) {
						   //deg, rms error of elevation measurement
                            double C1             = 0.0;
			    double a              = 0.0;
			    double b              = 0.0;
			    double dth            = 0.0;
			    constexpr double ae   = 8493333.0;
			    constexpr double PI   = 3.1415926535897932384626;
			    constexpr double sqt2 = 1.4142135623730950488017;
			    const double   rad1 = 57.2957795130823208767982;
			    const double   ae2  = ae+ae;
                            const double   xR   = pe.R;
			    const double   xha  = pe.ha;
			    const double   xht  = pe.ht;
			    const double   xthm = pe.thmax;
			    const double   xth3 = pe.th3;
			    const double   xps0 = pe.psi0;
			    const double   xpss = pe.psis;
			    const double   xpsv = pe.psiv;
			    const double   xkme = pe.kme;
			    const double   xNs  = pe.Ns;
			    const double   RR   = xR*xR;
			    const double   dp   = (4.0*ae*(xha+xht)+RR)*0.3333333333333333333;
			    switch(pe.flucts) {
                                 case 1: C1=0.000001;
				 break;
				 case 2: C1=2.0*0.000001;
				 break;
				 case 3: C1=4.0*0.000001;
				 break;
			    } 
			    const double   p    = std::sqrt(dp);
			    const double   p3   = p*p*p;
			    const double   dph  = (ae2*(xht-xha)*xR)/p3;
			    const double   ph   = std::acos(dph);
			    const double   cph  = std::cos((ph+PI)*0.3333333333333333333);
			    const double   d1   = xR*0.5-p*cph;
			    const double   H1   = xha-(d1*d1)/ae2;
			    const double   H2   = ((xR-d1)*H1)/d1;
			    const double   tht  = rad1*std::asin((H2-H1)/xR));
			    const double   ztht = rad1*tht;
			    const double   thef = ztht+(0.00025/(ztht+0.028));
			    const double   sthf = std::sin(thef);
			    const double   sig0 = C1/rad1*std::sqrt(3/sthf);
			    const double   cthf = std::cot(tht);
			    const double   lgth = std::log(tht);
			    if(tht<0.1) {
                               b = 0.0000001*cthf*(14.5+4.5*lgth);
			    }
			    else {
                               b = 0.000001g*cthf;
			    }
			    if(tht<0.02f) {
                               a = −0.00015*cthf*(1+0.3*lgth);
			    }
			    else if(tht>=0.02 && tht<=0.2) {
                               a = 0.00005*cthf*(1.0+1.4*lgth);
			    }
			    else if(tht>0.2) {
                               a = 0.0;
			    }
			    const double rtrm = 1.0+(0.001/((0.01+tht))*(0.01+tht));
			    if(tht>=0.0 && tht<=0.001) {
                               dth = 0.019;
			    }
			    else {
                               b*xNs+a*(xht/(11000.0+xht);
			    }
			    const double sig1 = dth/0.3490658503988659153847*rtrm;
			    const double sig12 = sig1*sig1;
			    Siger = sig12+sig12;
		    }
		    

    }

}

#endif /*__GMS_PLATFORM_ERROR_HPP__*/
