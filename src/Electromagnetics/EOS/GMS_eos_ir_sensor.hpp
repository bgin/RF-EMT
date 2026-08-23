
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
#ifndef __GMS_EOS_IR_SENSOR_HPP__
#define __GMS_EOS_IR_SENSOR_HPP__ 050920220913

namespace file_version {

    const unsigned int GMS_EOS_IR_SENSOR_MAJOR = 1U;
    const unsigned int GMS_EOS_IR_SENSOR_MINOR = 0U;
    const unsigned int GMS_EOS_IR_SENSOR_MICRO = 0U;
    const unsigned int GMS_EOS_IR_SENSOR_FULLVER =
      1000U*GMS_EOS_IR_SENSOR_MAJOR+
      100U*GMS_EOS_IR_SENSOR_MINOR+
      10U*GMS_EOS_IR_SENSOR_MICRO;
    const char * const GMS_EOS_IR_SENSOR_CREATION_DATE = "05-09-2022 09:13 AM +00200 (MON 05 SEP 2022 GMT+2)";
    const char * const GMS_EOS_IR_SENSOR_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_EOS_IR_SENSOR_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_EOS_IR_SENSOR_DESCRIPTION   = "Based on book: Mathematical Theory of Electro-Optical Sensors (rus)."

}


/*
   Various characteristics of Electro-Optical Sensors   
 ! Based mainly on Based mainly on Miroshenko M.M book (rus):          
 ! "Mathematical Theory of Electro-Optical Sensors".
*/

#include <cstdint>
#include <omp.h>
#include <cmath>
#include <limits>
#include "GMS_config.h"

namespace gms {

         namespace eos {



	      namespace {
                  // Taken from StackOverflow article.
		  // https://stackoverflow.com/questions/17333/what-is-the-most-effective-way-for-float-and-double-comparison
#pragma omp declare simd simdlen(16)
                  bool approximatelyEqual(const float a,
		                          const float b,
					  const float epsilon) {
			   const float fabsa = std::fabs(a);
			   const float fabsb = std::fabs(b);
                           return std::fabs(a - b) <=
                           ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
                   }
#pragma omp declare simd simdlen(8)
		   bool approximatelyEqual(const double a,
		                           const double b,
					   const double epsilon) {
			   const double fabsa = std::fabs(a);
			   const double fabsb = std::fabs(b);
                           return fabs(a - b) <=
                           ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
                   }
#pragma omp declare simd simdlen(16)
                  bool essentiallyEqual(const float a,
		                        const float b,
					const float epsilon) {
                           const float fabsa = std::fabs(a);
			   const float fabsb = std::fabs(b);
                           return fabs(a - b) <=
			   ((fabsa > fabsb ? fabsb : fabsa) * epsilon);
                   }
#pragma omp declare simd simdlen(8)
                   bool essentiallyEqual(const double a,
		                         const double b,
					 const double epsilon) {
                           const double fabsa = std::fabs(a);
			   const double fabsb = std::fabs(b);
                           return fabs(a - b) <=
			   ((fabsa > fabsb ? fabsb : fabsa) * epsilon);
                   }
#pragma omp declare simd simdlen(16)		   
                  bool definitelyGreaterThan(const float a,
		                             const float b,
					     const float epsilon) {
                           const float fabsa = std::fabs(a);
			   const float fabsb = std::fabs(b);
                           return (a - b) >
			   ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
		  }
#pragma omp declare simd simdlen(8)
		  bool definitelyGreaterThan(const double a,
		                             const double b,
					     const double epsilon) {
                           const double fabsa = std::fabs(a);
			   const double fabsb = std::fabs(b);
                           return (a - b) >
			   ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
		  }
#pragma omp declare simd simdlen(16)
                  bool definitelyLessThan(const float a,
		                          const float b,
					  const float epsilon) {
                           const float fabsa = std::fabs(a);
			   const float fabsb = std::fabs(b);
                           return (b - a) >
			   ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
                 }
#pragma omp declare simd simdlen(8)
		 bool definitelyLessThan( const double a,
		                          const double b,
					  const double epsilon) {
                           const double fabsa = std::fabs(a);
			   const double fabsb = std::fabs(b);
                           return (b - a) >
			   ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
                 }


	    }


	      __ATTR_ALWAYS_INLINE__
	     
	     
	      static inline
#pragma omp declare simd simdlen(16)
	      float param_gamma(const float phi) {

	            float gamma = 0.0f;
		    gamma       = 0.5f*phi*0.5f;
		    return (gamma);
	   }


              __ATTR_ALWAYS_INLINE__
	     
	     
	      static inline
#pragma omp declare simd simdlen(8)
	      double param_gamma(const double phi) {

	            double gamma = 0.0;
		    gamma        = 0.5*phi*0.5;
		    return (gamma);
	   }


	   //! Formula 1, p.54
           //!Тогда длина перпендикуляра SN, опущенного из 
           //!светящейся точки на плоскость зеркала
	      __ATTR_ALWAYS_INLINE__
	     
	     
	      static inline
#pragma omp declare simd simdlen(16)
	      float compute_SN(const float R,
	                       const float phi,
			       const float gamma,
			       const bool present) {

		     float SN = 0.0f;
		     if(present) {  //if true compute phi, otherwise compute gamma
                        SN = R*std::sin(phi);
		     }
		     else {
                        SN = R*std::cos(gamma);
		     }
		     return (SN);
	    }
	    

	      __ATTR_ALWAYS_INLINE__
	     
	     
	      static inline
#pragma omp declare simd simdlen(8)
	      double compute_SN(const double R,
	                        const double phi,
			        const double gamma,
			        const bool   present) {

		     double SN = 0.0;
		     if(present) {  //if true compute phi, otherwise compute gamma
                        SN = R*std::sin(phi);
		     }
		     else {
                        SN = R*std::cos(gamma);
		     }
		     return (SN);
	    }


	    //! Formula 2, p. 54
            //! расстояние SM от светящейся точки до ее изображения
	      __ATTR_ALWAYS_INLINE__
	     
	     
	      static inline
#pragma omp declare simd simdlen(16)
              float compute_SM(const float R,
	                       const float phi,
			       const float gamma,
			       const bool present) {

                  float SM = 0.0f;
		  SM       = 2.0f*compute_SN(R,phi,gamma,present);
		  return (SM);
	    }


	      __ATTR_ALWAYS_INLINE__
	     
	     
	      static inline
#pragma omp declare simd simdlen(8)
              double compute_SM(const double R,
	                        const double phi,
			        const double gamma,
			        const bool present) {

                  double SM = 0.0;
		  SM        = 2.0*compute_SN(R,phi,gamma,present);
		  return (SM);
	    }


	   //!Сканирующее зеркало для обеспечения осмотра всего поля
           //!обзора ф необходимо повернуть на угол, обеспечивающий 
           //!совмещение края изображения источника излучения с отверстием 
           //!диафрагмы а, находящимся в центре поля. Для этого необходимо 
           //!повернуть изображение светящейся точки S на угол ф/2
           //! Formula 1, p. 56
	      __ATTR_ALWAYS_INLINE__
	     
	     
	      static inline
#pragma omp declare simd simdlen(16)
              float ratio_FH(const float psi,
	                     const float phi) {

                    float FH = 0.0f;
		    float hpsi,hphi;
		    hpsi     = std::tan(0.5f*psi)
		    hphi     = std::tan(0.5f*phi);
		    FH       = hpsi/hphi;
		    return (FH);
	    }


	      __ATTR_ALWAYS_INLINE__
	     
	     
	      static inline
#pragma omp declare simd simdlen(8)
              double ratio_FH(const double psi,
	                      const double phi) {
			      
                    double FH = 0.0;
		    double hpsi,hphi;
		    hpsi     = std::tan(0.5*psi)
		    hphi     = std::tan(0.5*phi);
		    FH       = hpsi/hphi;
		    return (FH);
	    }


	    //! следовательно, угол установки сканирующего зеркала
            //! Formula 4, p. 56
	     __ATTR_ALWAYS_INLINE__
	    
	    
	     static inline
#pragma omp declare simd simdlen(16)
	     float scan_mirror_ang(const float gam0,
	                           const float psi,
				   const float phi,
				   const int32_t dir) { // value '1' -- positive, value '2' negative direction

                float gamma,t0,t1;
		t1     = 0.5f*phi*0.5f;
		gamma  = 0.0f;
		if(dir==1) {
                   t0    = gam0+t1;
		   gamma = t0*ratio_FH(psi,phi); 
		}
		else if(dir==2) {
                   t0    = gam0-t1;
		   gamma = t0*ratio_FH(psi,phi);
		}
		return (gamma);
	    }


	     __ATTR_ALWAYS_INLINE__
	    
	    
	     static inline
#pragma omp declare simd simdlen(8)
	     double scan_mirror_ang(const double gam0,
	                            const double psi,
				    const double phi,
				    const int32_t dir) { // value '1' -- positive, value '2' negative direction

                double gamma,t0,t1;
		t1     = 0.5*phi*0.5;
		gamma  = 0.0;
		if(dir==1) {
                   t0    = gam0+t1;
		   gamma = t0*ratio_FH(psi,phi); 
		}
		else if(dir==2) {
                   t0    = gam0-t1;
		   gamma = t0*ratio_FH(psi,phi);
		}
		return (gamma);
	    }


	    //! Maximum size of (verical) diameter of scanning mirror.
            //! Formula 2, page. 56, part: 1.3
            //! Anax = [h tg (6/2) + do6/2] [2 cos y' + sin yf {tg (у' + 6/2) +
            //!+ tg(Y'-6/2)}].
	     __ATTR_ALWAYS_INLINE__
	    
	    
	     static inline
#pragma omp declare simd simdlen(16)
             float compute_Dmax(const float h,
	                        const float delta,
				const float d_ob,
				const float gamma) {
             
                  float Dmax = 0.0f;
		  float delta2,d_ob2,cosg,sing,t0,t1,t2,tant0,tant1;
		  float t3,t4,t5;
		  delta2 = 0.5f*delta;
                  d_ob2  = 0.5f*d_ob;
                  cosg   = std::cos(gamma);
                  if(delta<=gamma) {
                     t0  = h*delta+d_ob;
                     Dmax= t0/cosg;
                 }
                  t0     = gamma+delta2;
                  t1     = gamma-delta2;
                  sing   = std::sin(gamma);
                  tant1  = std::tan(t0);
                  tant2  = std::tan(t1);
                  t3     = h*std::tan(delta2)+d_ob2;
                  t4     = 2.0f*cosg+sing;
                  t5     = tant1+tant2;
                  Dmax   = t3*t4*t5;
		  return (Dmax);
	   }


	   
             __ATTR_ALWAYS_INLINE__
	    
	    
	     static inline
#pragma omp declare simd simdlen(8)
             double compute_Dmax(const double h,
	                         const double delta,
				 const double d_ob,
				 const double gamma) {
             
                  double Dmax = 0.0;
		  double delta2,d_ob2,cosg,sing,t0,t1,t2,tant0,tant1;
		  double t3,t4,t5;
		  delta2 = 0.5*delta;
                  d_ob2  = 0.5*d_ob;
                  cosg   = std::cos(gamma);
                  if(delta<=gamma) {
                     t0  = h*delta+d_ob;
                     Dmax= t0/cosg;
                 }
                  t0     = gamma+delta2;
                  t1     = gamma-delta2;
                  sing   = std::sin(gamma);
                  tant1  = std::tan(t0);
                  tant2  = std::tan(t1);
                  t3     = h*std::tan(delta2)+d_ob2;
                  t4     = 2.0*cosg+sing;
                  t5     = tant1+tant2;
                  Dmax   = t3*t4*t5;
		  return (Dmax);
	   }


	   //! Размер зеркала в направлении, перпендикулярном плоскости
           //! чертежа, приблизительно равен
           //! Formula 2, p. 58
             __ATTR_ALWAYS_INLINE__
	    
	    
	     static inline
#pragma omp declare simd simdlen(16)
             float compute_Dmin(const float h,
	                        const float delta,
				const float d_ob) {

                    float Dmin = 0.0f;
		    Dmin       = h*delta+d_ob;
		    return (Dmin);
	  }


	     __ATTR_ALWAYS_INLINE__
	    
	    
	     static inline
#pragma omp declare simd simdlen(8)
             double compute_Dmin(const double h,
	                         const double delta,
				 const double d_ob) {

                    double Dmin = 0.0;
		    Dmin       = h*delta+d_ob;
		    return (Dmin);
	  }



          //!Если зеркало осуществляет сканирование в пространстве
          //!изображений его размеры
          //! Formula 7, p. 58
             __ATTR_ALWAYS_INLINE__
	    
	    
	     static inline
#pragma omp declare simd simdlen(16)
             float Dmax_imag_scan(const float H,
	                          const float F,
				  const float B,
				  const float d_ob,
				  const float gamma,
				  const float psi,
				  const float phi,
				  const float d) {

                 float Dmax = 0.0f;
		 float t0,t1,t2,t3;
		 float cosg,sing,tanp1,tanp2,psi2,phi2;
		 psi2  = 0.5f*psi;
                 if(psi2<=gamma && B<=d) {
                     phi2 = 0.5f*phi;
                     t0   = (F+F)*std::tan(phi2);
                     t1   = (H/F)*d_ob;
                     t2   = std::sin(gamma);
                     Dmax = (t0+t1)*t2;
                 }
                 t0    = (H/F)*(d_ob-B)+B;
                 cosg  = std::cos(gamma);
                 tanp1 = gamma+psi2;
                 tanp2 = gamma-psi2;
                 sing  = std::sin(gamma);
                 t1    = 2.0f*cosg+sing;
                 t2    = std::tan(tanp1)+std::tan(tanp2);
                 t3    = 0.5f*t1*t2;
                 Dmax  = t0*t3;
		 return (Dmax);
	  }


	     __ATTR_ALWAYS_INLINE__
	    
	    
	     static inline
#pragma omp declare simd simdlen(8)
             double Dmax_imag_scan(const double H,
	                           const double F,
				   const double B,
				   const double d_ob,
				   const double gamma,
				   const double psi,
				   const double phi,
				   const double d) {

                 double Dmax = 0.0;
		 double t0,t1,t2,t3;
		 double cosg,sing,tanp1,tanp2,psi2,phi2;
		 psi2  = 0.5*psi;
                 if(psi2<=gamma && B<=d) {
                     phi2 = 0.5*phi;
                     t0   = (F+F)*std::tan(phi2);
                     t1   = (H/F)*d_ob;
                     t2   = std::sin(gamma);
                     Dmax = (t0+t1)*t2;
                 }
                 t0    = (H/F)*(d_ob-B)+B;
                 cosg  = std::cos(gamma);
                 tanp1 = gamma+psi2;
                 tanp2 = gamma-psi2;
                 sing  = std::sin(gamma);
                 t1    = 2.0*cosg+sing;
                 t2    = std::tan(tanp1)+std::tan(tanp2);
                 t3    = 0.5*t1*t2;
                 Dmax  = t0*t3;
		 return (Dmax);
	  }


	     __ATTR_ALWAYS_INLINE__
	    
	    
	     static inline
#pragma omp declare simd simdlen(16)
             float Dmin_imag_scan(const float H,
	                          const float F,
				  const float d_ob,
				  const float B) {

                   float Dmin = 0.0f;
		   float t0,t1;
		   t0   = H/F;
                   t1   = (d_ob-B)+B;
                   Dmin = t0*t1;
		   return (Dmin);
	   }


	     __ATTR_ALWAYS_INLINE__
	    
	    
	     static inline
#pragma omp declare simd simdlen(16)
             double Dmin_imag_scan(const double H,
	                           const double F,
				   const double d_ob,
				   const double B) {

                   double Dmin = 0.0;
		   double t0,t1;
		   t0   = H/F;
                   t1   = (d_ob-B)+B;
                   Dmin = t0*t1;
		   return (Dmin);
	   }


	   //!величина расфокусировки
           //!Formula 1, p. 59
	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            float defocus_cof(const float l2,
	                      const float alpha,
			      const float O,
			      const int32_t inf) {

                 float df = 0.0f;
		 float cos2a,icos;
		 cos2a = std::cos(alpha+alpha);
                 icos  = 1.0f/cos2a;
                 if(inf) 
                    df    = l2/(icos-1.0f)*O;
                 else
                    df    = l2/(icos-1.0f);
                 return (df);
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            double defocus_cof(const double l2,
	                       const double alpha,
			       const double O,
			       const int32_t inf) {

                 double df = 0.0;
		 double cos2a,icos;
		 cos2a = std::cos(alpha+alpha);
                 icos  = 1.0/cos2a;
                 if(inf) 
                    df    = l2/(icos-1.0)*O;
                 else
                    df    = l2/(icos-1.0);
                 return (df);
	  }


	  //! Диаметр кружка рассеяния р
          //! Formula 3, p.59
	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            float circle_dispersion(const float d,
	                            const float l1,
				    const float l2,
				    const float alpha,
				    const float O,
				    const int32_t inf) {

                  float rho = 0.0f;
		  float t0,t1;
		  t0  = d/(l1+l2);
                  t1  = defocus_cof(l2,alpha,O,inf);
                  rho = t0*t1;
		  return (rho);
	   }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            double circle_dispersion(const double d,
	                             const double l1,
				     const double l2,
				     const double alpha,
				     const double O,
				     const int32_t inf) {

                  double rho = 0.0;
		  double t0,t1;
		  t0  = d/(l1+l2);
                  t1  = defocus_cof(l2,alpha,O,inf);
                  rho = t0*t1;
		  return (rho);
	   }


	   //!Formula 2, p. 59
	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            float circ_dispers_diam(const float l1,
	                            const float l2,
				    const float alpha,
				    const float O,
				    const int32_t inf) {

                 float ratio = 0.0f;
		 float t0,t1;
		 t0    = l1+l2;
                 t1    = defocus_cos(l2,alpha,O,inf);
                 ratio = t1/t0;
		 return (ratio);
	   }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            double circ_dispers_diam(const double l1,
	                             const double l2,
				     const double alpha,
				     const double O,
				     const int32_t inf) {

                 double ratio = 0.0;
		 double t0,t1;
		 t0    = l1+l2;
                 t1    = defocus_cos(l2,alpha,O,inf);
                 ratio = t1/t0;
		 return (ratio);
	   }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            float defocus_small_ang(const float O,
	                            const float l2,
				    const float alpha) {

                float rho = 0.0f;
		float t0,t1,t2,alpha2;
		const float eps = std::numeric_limits<float>::epsilon();
		alpha2 = alpha+alpha;
		t0     = std::cos(alpha2);
		t1     = 1.0f-(alpha2*alpha2*0.5f);
		if(approximatelyEqual(t0,t1,eps)) {
                   t2  = l2*0.5f;
		   rho = O*t2*alpha2*alpha2;
		   return (rho);
		}
		else {
                   rho = std::numeric_limits<float>::quiet_NaN();
		   return (rho);
		}
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            double defocus_small_ang(const double O,
	                             const double l2,
				     const double alpha) {

                double rho = 0.0;
		double t0,t1,t2,alpha2;
		const double eps = std::numeric_limits<double>::epsilon();
		alpha2 = alpha+alpha;
		t0     = std::cos(alpha2);
		t1     = 1.0-(alpha2*alpha2*0.5);
		if(approximatelyEqual(t0,t1,eps)) {
                   t2  = l2*0.5;
		   rho = O*t2*alpha2*alpha2;
		   return (rho);
		}
		else {
                   rho = std::numeric_limits<double>::quiet_NaN();
		   return (rho);
		}
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            float traj_scan_dxdt(const float dx[2],
	                         const float dt[2]) {
 
                  float dxdt = 0.0f;
		  dxdt = (dx[1]-dx[0])/(dt[1]-dt[0]);
		  return (dxdt);
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            double traj_scan_dxdt(const double dx[2],
	                         const double dt[2]) {
 
                  double dxdt = 0.0;
		  dxdt = (dx[1]-dx[0])/(dt[1]-dt[0]);
		  return (dxdt);
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            float traj_scan_dydt(const float dy[2],
	                         const float dt[2]) {
 
                  float dydt = 0.0f;
		  dxdt = (dy[1]-dy[0])/(dt[1]-dt[0]);
		  return (dydt);
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            double traj_scan_dydt(const double dy[2],
	                          const double dt[2]) {
 
                  double dydt = 0.0;
		  dxdt = (dy[1]-dy[0])/(dt[1]-dt[0]);
		  return (dydt);
	  }


	  //! СКАНИРОВАНИЕ ЗЕРКАЛОМ, ВРАЩАЮЩИМСЯ
          //! ВОКРУГ ОСИ, НЕПЕРПЕНДИКУЛЯРНОЙ К НЕМУ
          //! Formula 1, p. 100
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
	    float fov_x_axis(const float H,
	                     const float delta,
			     const float gamma) {

               float ax = 0.0f;
	       float gamm2,tdel;
	       gamm2    = 0.5f*gamma;
	       tdel     = std::tan(delta);
	       ax       = H*tdel*std::cos(gamm2);
	       return (ax);
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
	    double fov_x_axis(const double H,
	                      const double delta,
			      const double gamma) {

               double ax = 0.0;
	       double gamm2,tdel;
	       gamm2    = 0.5*gamma;
	       tdel     = std::tan(delta);
	       ax       = H*tdel*std::cos(gamm2);
	       return (ax);
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
	    float fov_y_axis(const float H,
	                     const float delta,
			     const float gamma) {

               float ay = 0.0f;
	       float ax,t0;
	       t0       = 0.5f*gamma;
	       ax       = fov_x_axis(H,delta,gamma);
	       ay       = t0*ax;
	       return (ay);
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
	    double fov_y_axis(const double H,
	                      const double delta,
			      const double gamma) {

               double ay = 0.0;
	       double ax,t0;
	       t0       = 0.5*gamma;
	       ax       = fov_x_axis(H,delta,gamma);
	       ay       = t0*ax;
	       return (ay);
	  }


	 //!Если рабочая зона сканирования ограничена углом G, то
         //!ширина захвата
         //!Formula 3, p. 100
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            float scan_width(const float H,
	                     const float gamma,
			     const float theta) {

                float B = 0.0f;
		float gam2,th2,t0,t1;
		gam2  = 0.5f*gamma;
                th2   = 0.5f*theta;
                t0    = std::tan(gam2);
                t1    = std::sin(th2);
                B     = (H+H)*t0*t1;
		return (B);
	   }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            double scan_width(const double H,
	                      const double gamma,
			      const double theta) {

                double B = 0.0;
		double gam2,th2,t0,t1;
		gam2  = 0.5*gamma;
                th2   = 0.5*theta;
                t0    = std::tan(gam2);
                t1    = std::sin(th2);
                B     = (H+H)*t0*t1;
		return (B);
	   }


	   //!Плоскопараллельная пластинка, установленная за 
           //!объективом, изменяет ход лучей таким образом, что изображение
           //! светящейся точки отодвигается и его положение зависит от угла у
           //!между оптической осью и нормалью N к поверхности пластинки
           //! Formula 7,8 p. 106
	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            float refract_shift(const float i1,
	                        const float delta,
				const float alfa,
				const float gamma,
				const float n) {

                 float l = 0.0f;
		 float ag,num,den,sin2,sag,t0,t1;
		 const float eps = std::numeric_limits<float>::epsilon();
		 if(approximatelyEqual(i1,ag,eps)) {
                     sag  = std::sin(ag);
                     t0   = delta*sag;
                     sin2 = sag*sag;
                     num  = 1.0f-sag;
                     den  = n*n-sag;
                     t1   = 1.0f-std::sqrt(num/den);
                     l    = t0*t2;
		 }
		 else if(alfa==0.0f) {
                     sag  = std::sin(gamma);
                     t0   = -delta*sag;
                     sin2 = sag*sag;
                     num  = 1.0f-sin2;
                     den  = n*n-sin2;
                     t1   = 1.0f-sqrt(num/den);
                     l    = t0*t1;
		 }
		 else {
                     sag  = std::sin(i1);
                     t0   = delta*sag;
                     sin2 = sag*sag;
                     num  = 1.0f-sin2;
                     den  = n*n-sin2;
                     t1   = 1.0f-sqrt(num/den);
                     l    = t0*t1;
		 }
		 return (l);
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            double refract_shift(const double i1,
	                         const double delta,
				 const double alfa,
				 const double gamma,
				 const double n) {

                 double l = 0.0;
		 double ag,num,den,sin2,sag,t0,t1;
		 const double eps = std::numeric_limits<double>::epsilon();
		 if(approximatelyEqual(i1,ag,eps)) {
                     sag  = std::sin(ag);
                     t0   = delta*sag;
                     sin2 = sag*sag;
                     num  = 1.0-sag;
                     den  = n*n-sag;
                     t1   = 1.0-std::sqrt(num/den);
                     l    = t0*t2;
		 }
		 else if(alfa==0.0) {
                     sag  = std::sin(gamma);
                     t0   = -delta*sag;
                     sin2 = sag*sag;
                     num  = 1.0-sin2;
                     den  = n*n-sin2;
                     t1   = 1.0-sqrt(num/den);
                     l    = t0*t1;
		 }
		 else {
                     sag  = std::sin(i1);
                     t0   = delta*sag;
                     sin2 = sag*sag;
                     num  = 1.0-sin2;
                     den  = n*n-sin2;
                     t1   = 1.0-sqrt(num/den);
                     l    = t0*t1;
		 }
		 return (l);
	  }


	  //!Formula 1, p. 108
	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            void project_xy_axis( const float l,
	                          const float alpha,
				  float &xl,
				  float &yl) {

                float absl = 0.0f;
		absl       = std::abs(l);
		xl         = std::cos(alpha);
		yl         = std::sin(alpha);
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            void project_xy_axis( const double l,
	                          const double alpha,
				  double &xl,
				  double &yl) {

                double absl = 0.0;
		absl       = std::abs(l);
		xl         = std::cos(alpha);
		yl         = std::sin(alpha);
	  }


	 //!Величину смещения луча s вдоль перпендикуляра к 
         //!поверхности пластинки
         //! Formula 2, p. 108
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            float s_shift(const float l,
	                  const float alpha,
		          const float gamma) {

              float s = 0.0f;
	      float ag,sag;
	      ag = alpha-gamma;
              sag= std::sin(ag);
              s  = l/sag;
	      return (s);
	 }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            double s_shift(const double l,
	                   const double alpha,
		           const double gamma) {

              double s = 0.0;
	      double ag,sag;
	      ag = alpha-gamma;
              sag= std::sin(ag);
              s  = l/sag;
	      return (s);
	 }


	 //! Проекции s на оси координат равны
         //! Formula 4, p. 108
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            void  project_s_xy(const float s,
	                       const float gamma,
			       float &xs,
			       float &ys) {

                 xs = s*std::cos(gamma);
                 ys = s*std::sin(gamma);
	 }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            void  project_s_xy(const double s,
	                       const double gamma,
			       double &xs,
			       double &ys) {

                 xs = s*std::cos(gamma);
                 ys = s*std::sin(gamma);
	 }


      //! что расстояния от начала координат О до точек
      //! пересечения лучей, образующих с горизонталью угла ±а, с 
      //! перпендикуляром к пластинке
      //! Formula 1, p. 110
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            float ray_intercept_pa(const float delta,
	                           const float alpha,
				   const float gamma,
				   const float n) {

               float sp = 0.0f;
	       float ag,num,den,n2,sag,sin2;
	       ag  = abs(alpha)-gamma;
               n2  = n*n;
               num = std::cos(ag);
               sag = std::sin(ag);
               sin2= std::sag*sag;
               den = std::sqrt(n2-sin2);
               sp  = std::delta*1.0f-(num/den);
	       return (sp);
	   }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            double ray_intercept_pa(const double delta,
	                            const double alpha,
				    const double gamma,
				    const double n) {

               double sp = 0.0;
	       double ag,num,den,n2,sag,sin2;
	       ag  = abs(alpha)-gamma;
               n2  = n*n;
               num = std::cos(ag);
               sag = std::sin(ag);
               sin2= std::sag*sag;
               den = std::sqrt(n2-sin2);
               sp  = std::delta*1.0-(num/den);
	       return (sp);
	   }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            float ray_intercept_na( const float delta,
	                            const float alpha,
				    const float gamma,
				    const float n) {

               float sn = 0.0f;
	       float ag,num,den,n2,sag,sin2;
	       ag  = abs(alpha)+gamma;
               n2  = n*n;
               num = std::cos(ag);
               sag = std::sin(ag);
               sin2= sag*sag;
               den = sqrt(n2-sin2);
               sn  = delta*1.0f-(num/den);
	       return (sn);
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            double ray_intercept_na( const double delta,
	                             const double alpha,
				     const double gamma,
				     const double n) {

               double sn = 0.0;
	       double ag,num,den,n2,sag,sin2;
	       ag  = abs(alpha)+gamma;
               n2  = n*n;
               num = std::cos(ag);
               sag = std::sin(ag);
               sin2= sag*sag;
               den = sqrt(n2-sin2);
               sn  = delta*1.0-(num/den);
	       return (sn);
	  }


	  //! Formula 3, p. 110
	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            float ray_diff(const float delta,
	                   const float alpha,
			   const float gamma,
			   const float n,
			   const float u) {

              float ds = 0.0f;
	      float t0,t1,u2,u2g,su2,sg,t2,n2,t3,t4,t5;
	      n   = n*n;
              u2  = u*0.5f;
              u2g = u2-gamma;
              t2  = std::sin(u2g);
              su2 = t2*t2;
              if(n2>=su2){
                 t3 = (-2.0f*delta)/n;
                 t4 = std::sin(u2);
                 t5 = std::sin(gamma);
                 ds = t3*t4*t5;
	      }
              else {
                 t0  = ray_intercept_pa(delta,alpha,gamma,n)
                 t1  = ray_intercept_na(delta,alpha,gamma,n)
                 ds  = t0-t1;
             }
	     return (ds);
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            double ray_diff(const double delta,
	                    const double alpha,
			    const double gamma,
			    const double n,
			    const double u) {

              double ds = 0.0;
	      double t0,t1,u2,u2g,su2,sg,t2,n2,t3,t4,t5;
	      n   = n*n;
              u2  = u*0.5;
              u2g = u2-gamma;
              t2  = std::sin(u2g);
              su2 = t2*t2;
              if(n2>=su2){
                 t3 = (-2.0*delta)/n;
                 t4 = std::sin(u2);
                 t5 = std::sin(gamma);
                 ds = t3*t4*t5;
	      }
              else {
                 t0  = ray_intercept_pa(delta,alpha,gamma,n)
                 t1  = ray_intercept_na(delta,alpha,gamma,n)
                 ds  = t0-t1;
             }
	     return (ds);
	  }


	  //!Поле точек пересечения лучей, преломленных пластинкой,
          //!относительно оси Ох (рис. 87) имеет симметрию, поэтому 
          //!упростим обозначения и выполним расчет соответствующих 
          //!координат на основании
          //! Formula 6,7, p. 111
	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            void compute_dxdy(const float alpha,
	                      const float beta,
			      const float delta,
			      const float gamma,
			      const float n,
			      const float u,
			      float &dx,
			      float &dy) {

                 float ag,ds,t0,t1,t2;
		 ag  = alpha+gamma;
                 ds  = ray_diff(delta,alfa,gamma,n,u);
                 t0  = sin(ag);
                 t1  = 2.0f*std::sin(alpha);
                 t2  = 2.0f*std::cos(alpha);
                 dx  = t0/t1*ds;
                 dy  = t0/t2*ds;
	   }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            void compute_dxdy(const double alpha,
	                      const double beta,
			      const double delta,
			      const double gamma,
			      const double n,
			      const double u,
			      double &dx,
			      double &dy) {

                 double ag,ds,t0,t1,t2;
		 ag  = alpha+gamma;
                 ds  = ray_diff(delta,alfa,gamma,n,u);
                 t0  = sin(ag);
                 t1  = 2.0*std::sin(alpha);
                 t2  = 2.0*std::cos(alpha);
                 dx  = t0/t1*ds;
                 dy  = t0/t2*ds;
	   }


	   //! Formula 7,8  p. 111
	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            void compute_xy(  const float alpha,
	                      const float beta,
			      const float delta,
			      const float gamma,
			      const float n,
			      const float u,
			      float &x,
			      float &y) {

                float sag,cag,pa,dx,dy,xs,ys;
		sag  = std::sin(gamma);
                cag  = std::cos(gamma);
                pa   = ray_intercept_pa(delta,alpha,gamma,n);
                xs   = pa*sag;
                ys   = pa*cag;
                compute_dxdy(alpha,beta,delta,gamma,n,u,dx,dy);
                x    = xs+dx;
                y    = ys+dx;
	   }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            void compute_xy(  const double alpha,
	                      const double beta,
			      const double delta,
			      const double gamma,
			      const double n,
			      const double u,
			      double &x,
			      double &y) {

                double sag,cag,pa,dx,dy,xs,ys;
		sag  = std::sin(gamma);
                cag  = std::cos(gamma);
                pa   = ray_intercept_pa(delta,alpha,gamma,n);
                xs   = pa*sag;
                ys   = pa*cag;
                compute_dxdy(alpha,beta,delta,gamma,n,u,dx,dy);
                x    = xs+dx;
                y    = ys+dx;
	   }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
	    void compute_xdyd(const float gamma,
	                      const float u,
			      const float n,
			      float &xd,
			      float &yd) {

                float cosg,sing,sin2s,sin2d,u2,u2gs,u2gd,n2,t0,t1,t2,t3,t4;
		cosg = std::cos(gamma);
                sing = std::sin(gamma);
                u2   = u*0.5f;
                n2   = n*n;
                u2gs = u2+gamma;
                ungd = u2-gamma;
                t0   = sin(u2gs);
                sin2s= t0*t0;
                t1   = std::sin(u2gd);
                sin2d= t1*t1;
                t2   = 1.0f/(4.0f*std::sin(u2));
                t3   = std::sqrt(n2-sin2s);
                t0   = sin2s/t3;
                t4   = std::sqrt(n2-sin2d);
                t1   = sin2d/t4;
                dx   = cosg-t2*(t0+t1);
                t2   = 1.0f/(4.0f*std::cos(u2));
                dy   = sing-t2*(t0-t1);
	   }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
	    void compute_xdyd(const double gamma,
	                      const double u,
			      const double n,
			      double &xd,
			      double &yd) {

                double cosg,sing,sin2s,sin2d,u2,u2gs,u2gd,n2,t0,t1,t2,t3,t4;
		cosg = std::cos(gamma);
                sing = std::sin(gamma);
                u2   = u*0.5;
                n2   = n*n;
                u2gs = u2+gamma;
                ungd = u2-gamma;
                t0   = sin(u2gs);
                sin2s= t0*t0;
                t1   = std::sin(u2gd);
                sin2d= t1*t1;
                t2   = 1.0f/(4.0*std::sin(u2));
                t3   = std::sqrt(n2-sin2s);
                t0   = sin2s/t3;
                t4   = std::sqrt(n2-sin2d);
                t1   = sin2d/t4;
                dx   = cosg-t2*(t0+t1);
                t2   = 1.0/(4.0f*std::cos(u2));
                dy   = sing-t2*(t0-t1);
	   }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            void paraxial_xdyd(const float gamma,
	                       const float alpha,
			       const float n,
			       float &xd,
			       float &yd) {

               float n2,cosg,sing,sin4g,sin2g,num,den,cos2g,t0,t1,n2ss;
	       n2    = n*n;
               cosg  = std::cos(gamma);
	       n2ss  = n2-sing*sing;
               cos2g = std::cos(gamma+gamma);
               sing  = std::sin(gamma);
               sin4g = sing*sing*sing*sing;
               num   = n2*cos2g+sin4g;
               den   = std::pow(n2ss,1.5f);
               xd    = cosg-num/den;
               t0    = sqrt(n2ss);
               t1    = 1.0f-cosg/t0;
               yd    = sing*t1;
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            void paraxial_xdyd(const double gamma,
	                       const double alpha,
			       const double n,
			       double &xd,
			       double &yd) {

               double n2,cosg,sing,sin4g,sin2g,num,den,cos2g,t0,t1,n2ss;
	       n2    = n*n;
               cosg  = std::cos(gamma);
	       n2ss  = n2-sing*sing;
               cos2g = std::cos(gamma+gamma);
               sing  = std::sin(gamma);
               sin4g = sing*sing*sing*sing;
               num   = n2*cos2g+sin4g;
               den   = std::pow(n2ss,1.5);
               xd    = cosg-num/den;
               t0    = sqrt(n2ss);
               t1    = 1.0-cosg/t0;
               yd    = sing*t1;
	  }


	  //!СКАНИРОВАНИЕ ВРАЩАЮЩИМИСЯ ОБЪЕКТИВАМИ
          //!Formula 1, p. 121
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            void fov_axay(const float H,
	                  const float delx,
			  const float dely,
			  const float phi,
			  float &ax,
			  float &ay) {

               float sec2,phi2,t0,t1,sec;
	       phi2  = 0.5f*phi;
               sec   = 1.0f/std::cos(phi2);
               sec2  = sec*sec;
               ax    = H*delx*sec2;
               ay    = H*dely*sec;
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            void fov_axay(const double H,
	                  const doube delx,
			  const double dely,
			  const double phi,
			  double &ax,
			  double &ay) {

               double sec2,phi2,t0,t1,sec;
	       phi2  = 0.5*phi;
               sec   = 1.0/std::cos(phi2);
               sec2  = sec*sec;
               ax    = H*delx*sec2;
               ay    = H*dely*sec;
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            void fov_dxdy(const float x,
	                  const float y,
			  const float F,
			  const float phi,
			  float &dx,
			  float &dy) {

              float d0x,d0y,phi2;
	      d0y   = y/F;
              phi2  = 0.5f*phi;
              d0x   = x/F;
              dy    = d0y;
              dx    = d0x*std::cos(phi2);
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            void fov_dxdy(const double x,
	                  const double y,
			  const double F,
			  const double phi,
			  double &dx,
			  double &dy) {

              double d0x,d0y,phi2;
	      d0y   = y/F;
              phi2  = 0.5*phi;
              d0x   = x/F;
              dy    = d0y;
              dx    = d0x*std::cos(phi2);
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            void volt_impulse_uxuy(const float u,
	                           const float om1,
				   const float om2,
				   const float t,
				   float &ux,
				   float &uy) {

                float om1t,om2t,t0,t1;
		om1t = om1*t;
                om2t = om2*t;
                t0   = std::sin(om1t)+std::sin(om2t);
                t1   = std::cos(om1t)+std::cos(om2t);
                ux   = u*t0;
                uy   = u*t1;
	   }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            void volt_impulse_uxuy(const double u,
	                           const double om1,
				   const double om2,
				   const double t,
				   double &ux,
				   double &uy) {

                double om1t,om2t,t0,t1;
		om1t = om1*t;
                om2t = om2*t;
                t0   = std::sin(om1t)+std::sin(om2t);
                t1   = std::cos(om1t)+std::cos(om2t);
                ux   = u*t0;
                uy   = u*t1;
	   }


	//! Phase Modulation
        //! Formula 1, p. 143
        //! растрового анализатора со 
        //!скрещивающимися осями, выполненного в виде надетой на вращающийся 
        //!барабан тонкой пленки, прозрачность которой изменяется по 
        //!синусоидальному закону
	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(16)
            float raster_transparency(const float rho_avg,
	                              const float rho_max,
				      const float rho_min,
				      const float l,
				      const float L,
				      const float N) {

               float rho = 0.0f;
	       constexpr float twopi = 6.283185307179586476925286766559f;
	       float t0,t1,t2;
	       t0  = 0.5f*(rho_max-rho_min);
               t1  = L/N;
               t2  = std::sin(twopi*l*t1);
               rho = rho_avg+t0*t2;
	       return (rho);
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
#pragma omp declare simd simdlen(8)
            double raster_transparency(const double rho_avg,
	                               const double rho_max,
				       const double rho_min,
				       const double l,
				       const double L,
				       const double N) {

               double rho = 0.0;
	       constexpr double twopi = 6.283185307179586476925286766559;
	       double t0,t1,t2;
	       t0  = 0.5*(rho_max-rho_min);
               t1  = L/N;
               t2  = std::sin(twopi*l*t1);
               rho = rho_avg+t0*t2;
	       return (rho);
	  }


#include "GMS_avint.hpp"
 

             //!СТРУКТУРА И СПЕКТР МОДУЛИРОВАННОГО ПОТОКА
             //!ИЗЛУЧЕНИЯ
             //!Formula 1, p. 178
             //! Ф(*) = Int rp(z,t)E(z,t) dsig
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void raster_flux_integral_omp(const float * __restrict rhoE,
	                                  const float * __restrict absc,
					  const int32_t n,
					  const int32_t t,
					  const float * __restrict xlo,
					  const float * __restrict xup,
					  float * __restrict Phit,
					  float * __restrict ier) {

                 float ans_x   = 0.0f;
		 int32_t err_x = 0;
#pragma omp parallel for default(none) schedule(runtime) \
                 private(i,ans,err) shared(t,rhoE,absc,n,xlo,xup,ier)
		 for(int32_t i = 0, i < t; ++i) {
                     avint(&rhoE[i*n],&absc[0],n,xlo,xup,ans,err);
		     Phit[i] = ans;
		     ier[i]  = err;
		 }
	   }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void raster_flux_integral_omp(const double * __restrict rhoE, // points to memory of size: (0:n-1,t)
	                                  const double * __restrict absc,
					  const int32_t n,
					  const int32_t t,
					  const double * __restrict xlo,
					  const double * __restrict xup,
					  double * __restrict Phit, // results data size: 0:t-1
					  double * __restrict ier) {

                 double ans_x   = 0.0f;
		 int32_t err_x = 0;
#pragma omp parallel for default(none) schedule(runtime) \
                 private(i,ans,err) shared(t,rhoE,absc,n,xlo,xup,ier)
		 for(int32_t i = 0, i < t; ++i) {
                     avint(&rhoE[i*n],&absc[0],n,xlo,xup,ans,err);
		     Phit[i] = ans;
		     ier[i]  = err;
		 }
	   }


	     __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void raster_flux_integral(    const float * __restrict rhoE,
	                                  const float * __restrict absc,
					  const int32_t n,
					  const int32_t t,
					  const float * __restrict xlo,
					  const float * __restrict xup,
					  float * __restrict Phit,
					  float * __restrict ier) {

                 float ans_x   = 0.0f;
		 int32_t err_x = 0;
		 for(int32_t i = 0, i < t; ++i) {
                     avint(&rhoE[i*n],&absc[0],n,xlo,xup,ans,err);
		     Phit[i] = ans;
		     ier[i]  = err;
		 }
	   }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void raster_flux_integral(    const double * __restrict rhoE, // points to memory of size: (0:n-1,t)
	                                  const double * __restrict absc,
					  const int32_t n,
					  const int32_t t,
					  const double * __restrict xlo,
					  const double * __restrict xup,
					  double * __restrict Phit, // results data size: 0:t-1
					  double * __restrict ier) {

                 double ans_x   = 0.0f;
		 int32_t err_x = 0;
		 for(int32_t i = 0, i < t; ++i) {
                     avint(&rhoE[i*n],&absc[0],n,xlo,xup,ans,err);
		     Phit[i] = ans;
		     ier[i]  = err;
		 }
	   }


	   //!! Formula 3, p. 180
	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void raster_opacity_integral_omp(const float invs,
	                                     const float * __restrict rhophi,
					     const float * __restrict absc,
					     const int32_t n,
					     const int32_t t,
					     const float xlo,
					     const float xup,
					     float * __restrict rho,
					     int32_t * __restrict ier) {

                 float ans   = 0.0f;
		 int32_t err = 0;
#pragma omp parallel for default(none) schedule(runtime) \
                 private(i,ans,err) \
		 shared(t,rhophi,n,xlo,xup,invs,absc,ier)
		 for(int32_t i = 0; i < t; ++i) {
                     avint(&rhophi[i*n],&absc[0],n,xlo,xup,ans,err);
		     rho[i] = invs*ans;
		     ier[i] = err;
		 }
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void raster_opacity_integral(    const float invs,
	                                     const float * __restrict rhophi,
					     const float * __restrict absc,
					     const int32_t n,
					     const int32_t t,
					     const float xlo,
					     const float xup,
					     float * __restrict rho,
					     int32_t * __restrict ier) {

                 float ans   = 0.0f;
		 int32_t err = 0;
		 for(int32_t i = 0; i < t; ++i) {
                     avint(&rhophi[i*n],&absc[0],n,xlo,xup,ans,err);
		     rho[i] = invs*ans;
		     ier[i] = err;
		 }
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void raster_opacity_integral_omp(const double invs,
	                                     const double * __restrict rhophi,
					     const double * __restrict absc,
					     const int32_t n,
					     const int32_t t,
					     const double xlo,
					     const double xup,
					     double * __restrict rho,
					     int32_t * __restrict ier) {

                 double ans   = 0.0;
		 int32_t err = 0;
#pragma omp parallel for default(none) schedule(runtime) \
                 private(i,ans,err) \
		 shared(t,rhophi,n,xlo,xup,invs,absc,ier)
		 for(int32_t i = 0; i < t; ++i) {
                     avint(&rhophi[i*n],&absc[0],n,xlo,xup,ans,err);
		     rho[i] = invs*ans;
		     ier[i] = err;
		 }
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void raster_opacity_integral(    const double invs,
	                                     const double * __restrict rhophi,
					     const float * __restrict absc,
					     const int32_t n,
					     const int32_t t,
					     const double xlo,
					     const double xup,
					     double * __restrict rho,
					     int32_t * __restrict ier) {

                 double ans   = 0.0;
		 int32_t err = 0;
		 for(int32_t i = 0; i < t; ++i) {
                     avint(&rhophi[i*n],&absc[0],n,xlo,xup,ans,err);
		     rho[i] = invs*ans;
		     ier[i] = err;
		 }
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void cos_series_unroll_16x(const float om0,
	                               const int32_t n,
				       float * __restrict __ATTR_ALIGN__(64) coss,
				       const float k) {

                float arg0,arg1,arg2,arg3,arg4,arg,arg6,arg7;
		float arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15;
		float t0,t1,t2,t3,t4,t5,t6,t7;
		float t8,t9,t10,t11,t12,t13,t14,t15;
		float kom0;
		int32_t i,m,m1;
		kom0 = k*om0;
		m    = n%16;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0      = (float)i;
		       arg0    = kom0*t0;
		       coss[i] = std::cos(arg0);
		   }
		   if(n<16) return;
		}
		m1 = m+1;
	        __assume_aligned(coss,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(4)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	      for(i = m1; i != n; i += 16) {
                   t0        = (float)i;
                   arg0      = kom0*t0;
                   coss[i]   = std::cos(arg0);
                   t1        = (float)i+1;
                   arg1      = kom0*t1;
                   coss[i+1] = std::cos(arg1);
                   t2        = (float)i+2;
                   arg2      = kom0*t2;
                   coss[i+2] = std::cos(arg2);
                   t3        = (float)i+3;
                   arg3      = kom0*t3;
                   coss[i+3] = std::cos(arg3);
                   t4        = (float)i+4;
                   arg4      = kom0*t4;
                   coss[i+4] = std::cos(arg4);
                   t5        = (float)i+5;
                   arg5      = kom0*t5;
                   coss[i+5] = std::cos(arg5);
                   t6        = (float)i+6;
                   arg6      = kom0*t6;
                   coss[i+6] = std::cos(arg6);
                   t7        = (float)i+7;
                   arg7      = kom0*t7;
                   coss[i+7] = std::cos(arg7);
                   t8        = (float)i+8;
                   arg8      = kom0*t8;
                   coss[i+8] = std::cos(arg8);
                   t9        = (float)i+9;
                   arg9      = kom0*t9;
                   coss[i+9] = std::cos(arg9);
                   t10       = (float)i+10;
                   arg10     = kom0*t10;
                   coss[i+10]= std::cos(arg10);
                   t11       = (float)i+11;
                   arg11     = kom0*t11;
                   coss[i+11]= std::cos(arg11);
                   t12       = (float)i+12;
                   arg12     = kom0*t12;
                   coss[i+12]= std::cos(arg12);
                   t13       = (float)i+13;
                   arg13     = kom0*t13;
                   coss[i+13]= std::cos(arg13);
                   t14       = (float)i+14;
                   arg14     = kom0*t14;
                   coss[i+14]= std::cos(arg14);
                   t15       = (float)i+15;
                   arg15     = kom0*t15;
                   coss[i+15]= std::cos(arg15);
	      }
	  }


	     __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void cos_series_unroll_16x(const double om0,
	                               const int32_t n,
				       double * __restrict __ATTR_ALIGN__(64) coss,
				       const float k) {

                double arg0,arg1,arg2,arg3,arg4,arg,arg6,arg7;
		double arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15;
		double t0,t1,t2,t3,t4,t5,t6,t7;
		double t8,t9,t10,t11,t12,t13,t14,t15;
		double kom0;
		int32_t i,m,m1;
		kom0 = k*om0;
		m    = n%16;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0      = (double)i;
		       arg0    = kom0*t0;
		       coss[i] = std::cos(arg0);
		   }
		   if(n<16) return;
		}
		m1 = m+1;
	        __assume_aligned(coss,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(8)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	      for(i = m1; i != n; i += 16) {
                   t0        = (double)i;
                   arg0      = kom0*t0;
                   coss[i]   = std::cos(arg0);
                   t1        = (double)i+1;
                   arg1      = kom0*t1;
                   coss[i+1] = std::cos(arg1);
                   t2        = (double)i+2;
                   arg2      = kom0*t2;
                   coss[i+2] = std::cos(arg2);
                   t3        = (double)i+3;
                   arg3      = kom0*t3;
                   coss[i+3] = std::cos(arg3);
                   t4        = (double)i+4;
                   arg4      = kom0*t4;
                   coss[i+4] = std::cos(arg4);
                   t5        = (double)i+5;
                   arg5      = kom0*t5;
                   coss[i+5] = std::cos(arg5);
                   t6        = (double)i+6;
                   arg6      = kom0*t6;
                   coss[i+6] = std::cos(arg6);
                   t7        = (double)i+7;
                   arg7      = kom0*t7;
                   coss[i+7] = std::cos(arg7);
                   t8        = (double)i+8;
                   arg8      = kom0*t8;
                   coss[i+8] = std::cos(arg8);
                   t9        = (double)i+9);
                   arg9      = kom0*t9;
                   coss[i+9] = std::cos(arg9);
                   t10       = (double)i+10;
                   arg10     = kom0*t10;
                   coss[i+10]= std::cos(arg10);
                   t11       = (double)i+11;
                   arg11     = kom0*t11;
                   coss[i+11]= std::cos(arg11);
                   t12       = (double)i+12;
                   arg12     = kom0*t12;
                   coss[i+12]= std::cos(arg12);
                   t13       = (double)i+13;
                   arg13     = kom0*t13;
                   coss[i+13]= std::cos(arg13);
                   t14       = (double)i+14;
                   arg14     = kom0*t14;
                   coss[i+14]= std::cos(arg14);
                   t15       = (double)i+15;
                   arg15     = kom0*t15;
                   coss[i+15]= std::cos(arg15);
	      }
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void cos_series_unroll_8x(const float om0,
	                               const int32_t n,
				       float * __restrict __ATTR_ALIGN__(64) coss,
				       const float k) {

                float arg0,arg1,arg2,arg3,arg4,arg,arg6,arg7;
		float t0,t1,t2,t3,t4,t5,t6,t7;
		float kom0;
		int32_t i,m,m1;
		kom0 = k*om0;
		m    = n%8;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0      = (float)i;
		       arg0    = kom0*t0;
		       coss[i] = std::cos(arg0);
		   }
		   if(n<8) return;
		}
		m1 = m+1;
	        __assume_aligned(coss,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(4)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	      for(i = m1; i != n; i += 8) {
                   t0        = (float)i;
                   arg0      = kom0*t0;
                   coss[i]   = std::cos(arg0);
                   t1        = (float)i+1;
                   arg1      = kom0*t1;
                   coss[i+1] = std::cos(arg1);
                   t2        = (float)i+2;
                   arg2      = kom0*t2;
                   coss[i+2] = std::cos(arg2);
                   t3        = (float)i+3;
                   arg3      = kom0*t3;
                   coss[i+3] = std::cos(arg3);
                   t4        = (float)i+4;
                   arg4      = kom0*t4;
                   coss[i+4] = std::cos(arg4);
                   t5        = (float)i+5;
                   arg5      = kom0*t5;
                   coss[i+5] = std::cos(arg5);
                   t6        = (float)i+6;
                   arg6      = kom0*t6;
                   coss[i+6] = std::cos(arg6);
                   t7        = (float)i+7;
                   arg7      = kom0*t7;
                   coss[i+7] = std::cos(arg7);
              }
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void cos_series_unroll_8x(const double om0,
	                               const int32_t n,
				       double * __restrict __ATTR_ALIGN__(64) coss,
				       const float k) {

                double arg0,arg1,arg2,arg3,arg4,arg,arg6,arg7;
		double t0,t1,t2,t3,t4,t5,t6,t7;
		double kom0;
		int32_t i,m,m1;
		kom0 = k*om0;
		m    = n%8;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0      = (double)i;
		       arg0    = kom0*t0;
		       coss[i] = std::cos(arg0);
		   }
		   if(n<8) return;
		}
		m1 = m+1;
	        __assume_aligned(coss,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(8)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	      for(i = m1; i != n; i += 8) {
                   t0        = (double)i;
                   arg0      = kom0*t0;
                   coss[i]   = std::cos(arg0);
                   t1        = (double)i+1;
                   arg1      = kom0*t1;
                   coss[i+1] = std::cos(arg1);
                   t2        = (double)i+2;
                   arg2      = kom0*t2;
                   coss[i+2] = std::cos(arg2);
                   t3        = (double)i+3;
                   arg3      = kom0*t3;
                   coss[i+3] = std::cos(arg3);
                   t4        = (double)i+4;
                   arg4      = kom0*t4;
                   coss[i+4] = std::cos(arg4);
                   t5        = (double)i+5;
                   arg5      = kom0*t5;
                   coss[i+5] = std::cos(arg5);
                   t6        = (double)i+6;
                   arg6      = kom0*t6;
                   coss[i+6] = std::cos(arg6);
                   t7        = (double)i+7;
                   arg7      = kom0*t7;
                   coss[i+7] = std::cos(arg7);
              }
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void cos_series_unroll_4x( const float om0,
	                               const int32_t n,
				       float * __restrict __ATTR_ALIGN__(64) coss,
				       const float k) {

                float arg0,arg1,arg2,arg3;
		float t0,t1,t2,t3;
		float kom0;
		int32_t i,m,m1;
		kom0 = k*om0;
		m    = n%4;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0      = (float)i;
		       arg0    = kom0*t0;
		       coss[i] = std::cos(arg0);
		   }
		   if(n<4) return;
		}
		m1 = m+1;
	        __assume_aligned(coss,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(4)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	      for(i = m1; i != n; i += 4) {
                   t0        = (float)i;
                   arg0      = kom0*t0;
                   coss[i]   = std::cos(arg0);
                   t1        = (float)i+1;
                   arg1      = kom0*t1;
                   coss[i+1] = std::cos(arg1);
                   t2        = (float)i+2;
                   arg2      = kom0*t2;
                   coss[i+2] = std::cos(arg2);
                   t3        = (float)i+3;
                   arg3      = kom0*t3;
                   coss[i+3] = std::cos(arg3);
              }
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void cos_series_unroll_4x(const double om0,
	                               const int32_t n,
				       double * __restrict __ATTR_ALIGN__(64) coss,
				       const float k) {

                double arg0,arg1,arg2,arg3;
		double t0,t1,t2,t3;
		double kom0;
		int32_t i,m,m1;
		kom0 = k*om0;
		m    = n%4;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0      = (double)i;
		       arg0    = kom0*t0;
		       coss[i] = std::cos(arg0);
		   }
		   if(n<4) return;
		}
		m1 = m+1;
	        __assume_aligned(coss,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(8)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	      for(i = m1; i != n; i += 4) {
                   t0        = (double)i;
                   arg0      = kom0*t0;
                   coss[i]   = std::cos(arg0);
                   t1        = (double)i+1;
                   arg1      = kom0*t1;
                   coss[i+1] = std::cos(arg1);
                   t2        = (double)i+2;
                   arg2      = kom0*t2;
                   coss[i+2] = std::cos(arg2);
                   t3        = (double)i+3;
                   arg3      = kom0*t3;
                   coss[i+3] = std::cos(arg3);
              }
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void cos_series_unroll_2x( const float om0,
	                               const int32_t n,
				       float * __restrict __ATTR_ALIGN__(64) coss,
				       const float k) {

                float arg0,arg1;
		float t0,t1;
		float kom0;
		int32_t i,m,m1;
		kom0 = k*om0;
		m    = n%2;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0      = (float)i;
		       arg0    = kom0*t0;
		       coss[i] = std::cos(arg0);
		   }
		   if(n<2) return;
		}
		m1 = m+1;
	        __assume_aligned(coss,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(4)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	      for(i = m1; i != n; i += 2) {
                   t0        = (float)i;
                   arg0      = kom0*t0;
                   coss[i]   = std::cos(arg0);
                   t1        = (float)i+1;
                   arg1      = kom0*t1;
                   coss[i+1] = std::cos(arg1);
              }
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void cos_series_unroll_2x(const double om0,
	                               const int32_t n,
				       double * __restrict __ATTR_ALIGN__(64) coss,
				       const float k) {

                double arg0,arg1;
		double t0,t1;
		double kom0;
		int32_t i,m,m1;
		kom0 = k*om0;
		m    = n%2;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0      = (double)i;
		       arg0    = kom0*t0;
		       coss[i] = std::cos(arg0);
		   }
		   if(n<2) return;
		}
		m1 = m+1;
	        __assume_aligned(coss,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(8)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	      for(i = m1; i != n; i += 2) {
                   t0        = (double)i;
                   arg0      = kom0*t0;
                   coss[i]   = std::cos(arg0);
                   t1        = (double)i+1;
                   arg1      = kom0*t1;
                   coss[i+1] = std::cos(arg1);
               }
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void sin_series_unroll_16x(const float om0,
	                               const int32_t n,
				       float * __restrict __ATTR_ALIGN__(64) sins,
				       const float k) {

                float arg0,arg1,arg2,arg3,arg4,arg,arg6,arg7;
		float arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15;
		float t0,t1,t2,t3,t4,t5,t6,t7;
		float t8,t9,t10,t11,t12,t13,t14,t15;
		float kom0;
		int32_t i,m,m1;
		kom0 = k*om0;
		m    = n%16;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0      = (float)i;
		       arg0    = kom0*t0;
		       sins[i] = std::sin(arg0);
		   }
		   if(n<16) return;
		}
		m1 = m+1;
	        __assume_aligned(sins,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(4)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	      for(i = m1; i != n; i += 16) {
                   t0        = (float)i;
                   arg0      = kom0*t0;
                   sins[i]   = std::sin(arg0);
                   t1        = (float)i+1;
                   arg1      = kom0*t1;
                   sins[i+1] = std::sin(arg1);
                   t2        = (float)i+2;
                   arg2      = kom0*t2;
                   sins[i+2] = std::sin(arg2);
                   t3        = (float)i+3;
                   arg3      = kom0*t3;
                   sins[i+3] = std::sin(arg3);
                   t4        = (float)i+4;
                   arg4      = kom0*t4;
                   sins[i+4] = std::sins(arg4);
                   t5        = (float)i+5;
                   arg5      = kom0*t5;
                   sins[i+5] = std::sin(arg5);
                   t6        = (float)i+6;
                   arg6      = kom0*t6;
                   sins[i+6] = std::sin(arg6);
                   t7        = (float)i+7;
                   arg7      = kom0*t7;
                   sins[i+7] = std::sin(arg7);
                   t8        = (float)i+8;
                   arg8      = kom0*t8;
                   sins[i+8] = std::sin(arg8);
                   t9        = (float)i+9;
                   arg9      = kom0*t9;
                   sins[i+9] = std::sin(arg9);
                   t10       = (float)i+10;
                   arg10     = kom0*t10;
                   sins[i+10]= std::sin(arg10);
                   t11       = (float)i+11;
                   arg11     = kom0*t11;
                   sins[i+11]= std::sin(arg11);
                   t12       = (float)i+12;
                   arg12     = kom0*t12;
                   sins[i+12]= std::sin(arg12);
                   t13       = (float)i+13;
                   arg13     = kom0*t13;
                   sins[i+13]= std::sin(arg13);
                   t14       = (float)i+14;
                   arg14     = kom0*t14;
                   sins[i+14]= std::sin(arg14);
                   t15       = (float)i+15;
                   arg15     = kom0*t15;
                   sins[i+15]= std::sin(arg15);
	      }
	  }



	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void sin_series_unroll_16x(const double om0,
	                               const int32_t n,
				       double * __restrict __ATTR_ALIGN__(64) sins,
				       const double k) {

                double arg0,arg1,arg2,arg3,arg4,arg,arg6,arg7;
		double arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15;
		double t0,t1,t2,t3,t4,t5,t6,t7;
		double t8,t9,t10,t11,t12,t13,t14,t15;
		double kom0;
		int32_t i,m,m1;
		kom0 = k*om0;
		m    = n%16;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0      = (double)i;
		       arg0    = kom0*t0;
		       sins[i] = std::sin(arg0);
		   }
		   if(n<16) return;
		}
		m1 = m+1;
	        __assume_aligned(sins,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(8)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	      for(i = m1; i != n; i += 16) {
                   t0        = (double)i;
                   arg0      = kom0*t0;
                   sins[i]   = std::sin(arg0);
                   t1        = (double)i+1;
                   arg1      = kom0*t1;
                   sins[i+1] = std::sin(arg1);
                   t2        = (double)i+2;
                   arg2      = kom0*t2;
                   sins[i+2] = std::sin(arg2);
                   t3        = (double)i+3;
                   arg3      = kom0*t3;
                   sins[i+3] = std::sin(arg3);
                   t4        = (double)i+4;
                   arg4      = kom0*t4;
                   sins[i+4] = std::sins(arg4);
                   t5        = (double)i+5;
                   arg5      = kom0*t5;
                   sins[i+5] = std::sin(arg5);
                   t6        = (double)i+6;
                   arg6      = kom0*t6;
                   sins[i+6] = std::sin(arg6);
                   t7        = (double)i+7;
                   arg7      = kom0*t7;
                   sins[i+7] = std::sin(arg7);
                   t8        = (double)i+8;
                   arg8      = kom0*t8;
                   sins[i+8] = std::sin(arg8);
                   t9        = (double)i+9;
                   arg9      = kom0*t9;
                   sins[i+9] = std::sin(arg9);
                   t10       = (double)i+10;
                   arg10     = kom0*t10;
                   sins[i+10]= std::sin(arg10);
                   t11       = (double)i+11;
                   arg11     = kom0*t11;
                   sins[i+11]= std::sin(arg11);
                   t12       = (double)i+12;
                   arg12     = kom0*t12;
                   sins[i+12]= std::sin(arg12);
                   t13       = (double)i+13;
                   arg13     = kom0*t13;
                   sins[i+13]= std::sin(arg13);
                   t14       = (double)i+14;
                   arg14     = kom0*t14;
                   sins[i+14]= std::sin(arg14);
                   t15       = (double)i+15;
                   arg15     = kom0*t15;
                   sins[i+15]= std::sin(arg15);
	      }
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void sin_series_unroll_8x(const float om0,
	                               const int32_t n,
				       float * __restrict __ATTR_ALIGN__(64) sins,
				       const float k) {

                float arg0,arg1,arg2,arg3,arg4,arg,arg6,arg7;
		float t0,t1,t2,t3,t4,t5,t6,t7;
		float kom0;
		int32_t i,m,m1;
		kom0 = k*om0;
		m    = n%8;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0      = (float)i;
		       arg0    = kom0*t0;
		       sins[i] = std::sin(arg0);
		   }
		   if(n<8) return;
		}
		m1 = m+1;
	        __assume_aligned(sins,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(4)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	      for(i = m1; i != n; i += 8) {
                   t0        = (float)i;
                   arg0      = kom0*t0;
                   sins[i]   = std::sin(arg0);
                   t1        = (float)i+1;
                   arg1      = kom0*t1;
                   sins[i+1] = std::sin(arg1);
                   t2        = (float)i+2;
                   arg2      = kom0*t2;
                   sins[i+2] = std::sin(arg2);
                   t3        = (float)i+3;
                   arg3      = kom0*t3;
                   sins[i+3] = std::sin(arg3);
                   t4        = (float)i+4;
                   arg4      = kom0*t4;
                   sins[i+4] = std::sins(arg4);
                   t5        = (float)i+5;
                   arg5      = kom0*t5;
                   sins[i+5] = std::sin(arg5);
                   t6        = (float)i+6;
                   arg6      = kom0*t6;
                   sins[i+6] = std::sin(arg6);
                   t7        = (float)i+7;
                   arg7      = kom0*t7;
                   sins[i+7] = std::sin(arg7);
              }
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void sin_series_unroll_8x( const double om0,
	                               const int32_t n,
				       double * __restrict __ATTR_ALIGN__(64) sins,
				       const double k) {

                double arg0,arg1,arg2,arg3,arg4,arg,arg6,arg7;
		double t0,t1,t2,t3,t4,t5,t6,t7;
		double kom0;
		int32_t i,m,m1;
		kom0 = k*om0;
		m    = n%8;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0      = (double)i;
		       arg0    = kom0*t0;
		       sins[i] = std::sin(arg0);
		   }
		   if(n<8) return;
		}
		m1 = m+1;
	        __assume_aligned(sins,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(8)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	      for(i = m1; i != n; i += 8) {
                   t0        = (double)i;
                   arg0      = kom0*t0;
                   sins[i]   = std::sin(arg0);
                   t1        = (double)i+1;
                   arg1      = kom0*t1;
                   sins[i+1] = std::sin(arg1);
                   t2        = (double)i+2;
                   arg2      = kom0*t2;
                   sins[i+2] = std::sin(arg2);
                   t3        = (double)i+3;
                   arg3      = kom0*t3;
                   sins[i+3] = std::sin(arg3);
                   t4        = (double)i+4;
                   arg4      = kom0*t4;
                   sins[i+4] = std::sins(arg4);
                   t5        = (double)i+5;
                   arg5      = kom0*t5;
                   sins[i+5] = std::sin(arg5);
                   t6        = (double)i+6;
                   arg6      = kom0*t6;
                   sins[i+6] = std::sin(arg6);
                   t7        = (double)i+7;
                   arg7      = kom0*t7;
                   sins[i+7] = std::sin(arg7);
               }
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void sin_series_unroll_4x(const float om0,
	                               const int32_t n,
				       float * __restrict __ATTR_ALIGN__(64) sins,
				       const float k) {

                float arg0,arg1,arg2,arg3;
		float t0,t1,t2,t3;
		float kom0;
		int32_t i,m,m1;
		kom0 = k*om0;
		m    = n%4;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0      = (float)i;
		       arg0    = kom0*t0;
		       sins[i] = std::sin(arg0);
		   }
		   if(n<4) return;
		}
		m1 = m+1;
	        __assume_aligned(sins,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(4)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	      for(i = m1; i != n; i += 4) {
                   t0        = (float)i;
                   arg0      = kom0*t0;
                   sins[i]   = std::sin(arg0);
                   t1        = (float)i+1;
                   arg1      = kom0*t1;
                   sins[i+1] = std::sin(arg1);
                   t2        = (float)i+2;
                   arg2      = kom0*t2;
                   sins[i+2] = std::sin(arg2);
                   t3        = (float)i+3;
                   arg3      = kom0*t3;
                   sins[i+3] = std::sin(arg3);
               }
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void sin_series_unroll_4x( const double om0,
	                               const int32_t n,
				       double * __restrict __ATTR_ALIGN__(64) sins,
				       const double k) {

                double arg0,arg1,arg2,arg3;
		double t0,t1,t2,t3;
		double kom0;
		int32_t i,m,m1;
		kom0 = k*om0;
		m    = n%4;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0      = (double)i;
		       arg0    = kom0*t0;
		       sins[i] = std::sin(arg0);
		   }
		   if(n<4) return;
		}
		m1 = m+1;
	        __assume_aligned(sins,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(8)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	      for(i = m1; i != n; i += 4) {
                   t0        = (double)i;
                   arg0      = kom0*t0;
                   sins[i]   = std::sin(arg0);
                   t1        = (double)i+1;
                   arg1      = kom0*t1;
                   sins[i+1] = std::sin(arg1);
                   t2        = (double)i+2;
                   arg2      = kom0*t2;
                   sins[i+2] = std::sin(arg2);
                   t3        = (double)i+3;
                   arg3      = kom0*t3;
                   sins[i+3] = std::sin(arg3);
               }
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void sin_series_unroll_2x(const float om0,
	                               const int32_t n,
				       float * __restrict __ATTR_ALIGN__(64) sins,
				       const float k) {

                float arg0,arg1;
		float t0,t1;
		float kom0;
		int32_t i,m,m1;
		kom0 = k*om0;
		m    = n%2;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0      = (float)i;
		       arg0    = kom0*t0;
		       sins[i] = std::sin(arg0);
		   }
		   if(n<2) return;
		}
		m1 = m+1;
	        __assume_aligned(sins,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(4)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	      for(i = m1; i != n; i += 2) {
                   t0        = (float)i;
                   arg0      = kom0*t0;
                   sins[i]   = std::sin(arg0);
                   t1        = (float)i+1;
                   arg1      = kom0*t1;
                   sins[i+1] = std::sin(arg1);
              }
	  }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void sin_series_unroll_2x( const double om0,
	                               const int32_t n,
				       double * __restrict __ATTR_ALIGN__(64) sins,
				       const double k) {

                double arg0,arg1;
		double t0,t1;
		double kom0;
		int32_t i,m,m1;
		kom0 = k*om0;
		m    = n%2;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0      = (double)i;
		       arg0    = kom0*t0;
		       sins[i] = std::sin(arg0);
		   }
		   if(n<4) return;
		}
		m1 = m+1;
	        __assume_aligned(sins,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(8)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	      for(i = m1; i != n; i += 2) {
                   t0        = (double)i;
                   arg0      = kom0*t0;
                   sins[i]   = std::sin(arg0);
                   t1        = (double)i+1;
                   arg1      = kom0*t1;
                   sins[i+1] = std::sin(arg1);
              }
	  }


	   //! форма импульса потока излучения описывается,
           //! например, косинус-квадратной зависимостью
           //! Formula 3, p. 184
	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void squared_cos_flux_unroll_16x(float * __restrict __ATTR_ALIGN__(64) Phi0t,
	                                     const float Phi0,
					     const int32_t n,
					     const float tin) {

                constexpr float hpi = 1.57079632679489661923132169164f;
		float tin2,t0,t1,t2,t3,t4,t5,t6,t7;
                float t8,t9,t10,t11,t12,t13,t14,t15;
                float arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7;
                float arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15;
                float c0,c1,c2,c3,c4,c5,c6,c7;
                float c8,c9,c10,c11,c12,c13,c14,c15;
                tin2 = 0.5f*tin;
		int32_t i,m,m1;
		m = n%16;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0       = (float)i;
		       arg0     = hpi*t0/tin2;
		       c0       = std::cos(arg0);
		       Phi0t[i] = c0*c0;
		   }
		   if(n<16) return;
		}
		m1 = m+1;
		__assume_aligned(Phi0t,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(4)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	        for(i = m1; i != n; i += 16) {
                      t0         = (float)i;
		      arg0       = hpi*t0/tin2;
                      c0         = std::cos(arg0);
                      Phi0t[i]   = c0*c0;
                      t1         = (float)i+1;
                      arg1       = hpi*t1/tin2;
                      c1         = std::cos(arg1);
                      Phi0t[i+1] = c1*c1;
                      t2         = (float)i+2;
                      arg2       = hpi*t2/tin2;
                      c2         = std::cos(arg2);
                      Phi0t[i+2] = c2*c2;
                      t3         = (float)i+3;
                      arg3       = hpi*t3/tin2;
                      c3         = std::cos(arg3);
                      Phi0t[i+3] = c3*c3;
                      t4         = (float)i+4;
                      arg4       = hpi*t4/tin2;
                      c4         = std::cos(arg4);
                      Phi0t[i+4] = c4*c4;
                      t5         = (float)i+5;
                      arg5       = hpi*t5/tin2;
                      c5         = std::cos(arg5);
                      Phi0t[i+5] = c5*c5;
                      t6         = (float)i+6;
                      arg6       = hpi*t6/tin2;
                      c6         = std::cos(arg6);
                      Phi0t[i+6] = c6*c6;
                      t7         = (float)i+7;
		      arg7       = hpi*t7/tin2;
                      c7         = std::cos(arg7);
                      Phi0t[i+7] = c7*c7;
                      t8         = (float)i+8;
                      arg8       = hpi*t8/tin2;
                      c8         = std::cos(arg8);
                      Phi0t[i+8] = c8*c8;
                      t9         = (float)i+9;
                      arg9       = hpi*t9/tin2;
                      c9         = std::cos(arg9);
                      Phi0t[i+9] = c9*c9;
                      t10        = (float)i+10;
                      arg10      = hpi*t10/tin2;
                      c10        = std::cos(arg10);
                      Phi0t[i+10]= c10*c10;
                      t11        = (float)i+11;
                      arg11      = hpi*t11/tin2;
                      c11        = std::cos(arg11);
                      Phi0t[i+11]= c11*c11;
                      t12        = (float)i+12;
                      arg12      = hpi*t12/tin2;
                      c12        = std::cos(arg12);
                      Phi0t[i+12]= c12*c12;
                      t13        = (float)i+13;
                      arg13      = hpi*t13/tin2;
                      c13        = std::cos(arg13);
                      Phi0t[i+13]= c13*c13;
                      t14        = (float)i+14;
                      arg14      = hpi*t14/tin2;
                      c14        = std::cos(arg14);
                      Phi0t[i+14]= c14*c14;
                      t15        = (float)i+15;
                      arg15      = hpi*t15/tin2;
                      c15        = std::cos(arg15);
                      Phi0t[i+15]= c15*c15;
       

		}
	 }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void squared_cos_flux_unroll_16x(double * __restrict __ATTR_ALIGN__(64) Phi0t,
	                                     const double Phi0,
					     const int32_t n,
					     const double tin) {

                constexpr double hpi = 1.57079632679489661923132169164;
		double tin2,t0,t1,t2,t3,t4,t5,t6,t7;
                double t8,t9,t10,t11,t12,t13,t14,t15;
                double arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7;
                double arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15;
                double c0,c1,c2,c3,c4,c5,c6,c7;
                double c8,c9,c10,c11,c12,c13,c14,c15;
                tin2 = 0.5f*tin;
		int32_t i,m,m1;
		m = n%16;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0       = (double)i;
		       arg0     = hpi*t0/tin2;
		       c0       = std::cos(arg0);
		       Phi0t[i] = c0*c0;
		   }
		   if(n<16) return;
		}
		m1 = m+1;
		__assume_aligned(Phi0t,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(8)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	        for(i = m1; i != n; i += 16) {
                      t0         = (double)i;
		      arg0       = hpi*t0/tin2;
                      c0         = std::cos(arg0);
                      Phi0t[i]   = c0*c0;
                      t1         = (double)i+1;
                      arg1       = hpi*t1/tin2;
                      c1         = std::cos(arg1);
                      Phi0t[i+1] = c1*c1;
                      t2         = (double)i+2;
                      arg2       = hpi*t2/tin2;
                      c2         = std::cos(arg2);
                      Phi0t[i+2] = c2*c2;
                      t3         = (double)i+3;
                      arg3       = hpi*t3/tin2;
                      c3         = std::cos(arg3);
                      Phi0t[i+3] = c3*c3;
                      t4         = (double)i+4;
                      arg4       = hpi*t4/tin2;
                      c4         = std::cos(arg4);
                      Phi0t[i+4] = c4*c4;
                      t5         = (double)i+5;
                      arg5       = hpi*t5/tin2;
                      c5         = std::cos(arg5);
                      Phi0t[i+5] = c5*c5;
                      t6         = (double)i+6;
                      arg6       = hpi*t6/tin2;
                      c6         = std::cos(arg6);
                      Phi0t[i+6] = c6*c6;
                      t7         = (double)i+7;
		      arg7       = hpi*t7/tin2;
                      c7         = std::cos(arg7);
                      Phi0t[i+7] = c7*c7;
                      t8         = (double)i+8;
                      arg8       = hpi*t8/tin2;
                      c8         = std::cos(arg8);
                      Phi0t[i+8] = c8*c8;
                      t9         = (double)i+9;
                      arg9       = hpi*t9/tin2;
                      c9         = std::cos(arg9);
                      Phi0t[i+9] = c9*c9;
                      t10        = (double)i+10;
                      arg10      = hpi*t10/tin2;
                      c10        = std::cos(arg10);
                      Phi0t[i+10]= c10*c10;
                      t11        = (double)i+11;
                      arg11      = hpi*t11/tin2;
                      c11        = std::cos(arg11);
                      Phi0t[i+11]= c11*c11;
                      t12        = (double)i+12;
                      arg12      = hpi*t12/tin2;
                      c12        = std::cos(arg12);
                      Phi0t[i+12]= c12*c12;
                      t13        = (double)i+13;
                      arg13      = hpi*t13/tin2;
                      c13        = std::cos(arg13);
                      Phi0t[i+13]= c13*c13;
                      t14        = (double)i+14;
                      arg14      = hpi*t14/tin2;
                      c14        = std::cos(arg14);
                      Phi0t[i+14]= c14*c14;
                      t15        = (double)i+15;
                      arg15      = hpi*t15/tin2;
                      c15        = std::cos(arg15);
                      Phi0t[i+15]= c15*c15;
       

		}
	 }


            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void squared_cos_flux_unroll_8x(float * __restrict __ATTR_ALIGN__(64) Phi0t,
	                                     const float Phi0,
					     const int32_t n,
					     const float tin) {

                constexpr float hpi = 1.57079632679489661923132169164f;
		float tin2,t0,t1,t2,t3,t4,t5,t6,t7;
                float arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7;
                float c0,c1,c2,c3,c4,c5,c6,c7;
                tin2 = 0.5f*tin;
		int32_t i,m,m1;
		m = n%8;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0       = (float)i;
		       arg0     = hpi*t0/tin2;
		       c0       = std::cos(arg0);
		       Phi0t[i] = c0*c0;
		   }
		   if(n<8) return;
		}
		m1 = m+1;
		__assume_aligned(Phi0t,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(4)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	        for(i = m1; i != n; i += 8) {
                      t0         = (float)i;
		      arg0       = hpi*t0/tin2;
                      c0         = std::cos(arg0);
                      Phi0t[i]   = c0*c0;
                      t1         = (float)i+1;
                      arg1       = hpi*t1/tin2;
                      c1         = std::cos(arg1);
                      Phi0t[i+1] = c1*c1;
                      t2         = (float)i+2;
                      arg2       = hpi*t2/tin2;
                      c2         = std::cos(arg2);
                      Phi0t[i+2] = c2*c2;
                      t3         = (float)i+3;
                      arg3       = hpi*t3/tin2;
                      c3         = std::cos(arg3);
                      Phi0t[i+3] = c3*c3;
                      t4         = (float)i+4;
                      arg4       = hpi*t4/tin2;
                      c4         = std::cos(arg4);
                      Phi0t[i+4] = c4*c4;
                      t5         = (float)i+5;
                      arg5       = hpi*t5/tin2;
                      c5         = std::cos(arg5);
                      Phi0t[i+5] = c5*c5;
                      t6         = (float)i+6;
                      arg6       = hpi*t6/tin2;
                      c6         = std::cos(arg6);
                      Phi0t[i+6] = c6*c6;
                      t7         = (float)i+7;
		      arg7       = hpi*t7/tin2;
                      c7         = std::cos(arg7);
                      Phi0t[i+7] = c7*c7;
                
       
		}
	 }


	    __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void squared_cos_flux_unroll_8x(double * __restrict __ATTR_ALIGN__(64) Phi0t,
	                                     const double Phi0,
					     const int32_t n,
					     const double tin) {

                constexpr double hpi = 1.57079632679489661923132169164;
		double tin2,t0,t1,t2,t3,t4,t5,t6,t7;
                double arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7;
                double c0,c1,c2,c3,c4,c5,c6,c7;
                tin2 = 0.5f*tin;
		int32_t i,m,m1;
		m = n%8;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0       = (double)i;
		       arg0     = hpi*t0/tin2;
		       c0       = std::cos(arg0);
		       Phi0t[i] = c0*c0;
		   }
		   if(n<8) return;
		}
		m1 = m+1;
		__assume_aligned(Phi0t,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(8)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	        for(i = m1; i != n; i += 8) {
                      t0         = (double)i;
		      arg0       = hpi*t0/tin2;
                      c0         = std::cos(arg0);
                      Phi0t[i]   = c0*c0;
                      t1         = (double)i+1;
                      arg1       = hpi*t1/tin2;
                      c1         = std::cos(arg1);
                      Phi0t[i+1] = c1*c1;
                      t2         = (double)i+2;
                      arg2       = hpi*t2/tin2;
                      c2         = std::cos(arg2);
                      Phi0t[i+2] = c2*c2;
                      t3         = (double)i+3;
                      arg3       = hpi*t3/tin2;
                      c3         = std::cos(arg3);
                      Phi0t[i+3] = c3*c3;
                      t4         = (double)i+4;
                      arg4       = hpi*t4/tin2;
                      c4         = std::cos(arg4);
                      Phi0t[i+4] = c4*c4;
                      t5         = (double)i+5;
                      arg5       = hpi*t5/tin2;
                      c5         = std::cos(arg5);
                      Phi0t[i+5] = c5*c5;
                      t6         = (double)i+6;
                      arg6       = hpi*t6/tin2;
                      c6         = std::cos(arg6);
                      Phi0t[i+6] = c6*c6;
                      t7         = (double)i+7;
		      arg7       = hpi*t7/tin2;
                      c7         = std::cos(arg7);
                      Phi0t[i+7] = c7*c7;
                  
		}
	 }


            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void squared_cos_flux_unroll_4x( float * __restrict __ATTR_ALIGN__(64) Phi0t,
	                                     const float Phi0,
					     const int32_t n,
					     const float tin) {

                constexpr float hpi = 1.57079632679489661923132169164f;
		float tin2,t0,t1,t2,t3;
                float arg0,arg1,arg2,arg3;
                float c0,c1,c2,c3;
                tin2 = 0.5f*tin;
		int32_t i,m,m1;
		m = n%4;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0       = (float)i;
		       arg0     = hpi*t0/tin2;
		       c0       = std::cos(arg0);
		       Phi0t[i] = c0*c0;
		   }
		   if(n<4) return;
		}
		m1 = m+1;
		__assume_aligned(Phi0t,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(4)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	        for(i = m1; i != n; i += 4) {
                      t0         = (float)i;
		      arg0       = hpi*t0/tin2;
                      c0         = std::cos(arg0);
                      Phi0t[i]   = c0*c0;
                      t1         = (float)i+1;
                      arg1       = hpi*t1/tin2;
                      c1         = std::cos(arg1);
                      Phi0t[i+1] = c1*c1;
                      t2         = (float)i+2;
                      arg2       = hpi*t2/tin2;
                      c2         = std::cos(arg2);
                      Phi0t[i+2] = c2*c2;
                      t3         = (float)i+3;
                      arg3       = hpi*t3/tin2;
                      c3         = std::cos(arg3);
                      Phi0t[i+3] = c3*c3;
                                   
       
		}
	 }


         
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void squared_cos_flux_unroll_4x(double * __restrict __ATTR_ALIGN__(64) Phi0t,
	                                     const double Phi0,
					     const int32_t n,
					     const double tin) {

                constexpr double hpi = 1.57079632679489661923132169164;
		double tin2,t0,t1,t2,t3;
                double arg0,arg1,arg2,arg3;
                double c0,c1,c2,c3;
                tin2 = 0.5f*tin;
		int32_t i,m,m1;
		m = n%4;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0       = (double)i;
		       arg0     = hpi*t0/tin2;
		       c0       = std::cos(arg0);
		       Phi0t[i] = c0*c0;
		   }
		   if(n<4) return;
		}
		m1 = m+1;
		__assume_aligned(Phi0t,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(8)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	        for(i = m1; i != n; i += 4) {
                      t0         = (double)i;
		      arg0       = hpi*t0/tin2;
                      c0         = std::cos(arg0);
                      Phi0t[i]   = c0*c0;
                      t1         = (double)i+1;
                      arg1       = hpi*t1/tin2;
                      c1         = std::cos(arg1);
                      Phi0t[i+1] = c1*c1;
                      t2         = (double)i+2;
                      arg2       = hpi*t2/tin2;
                      c2         = std::cos(arg2);
                      Phi0t[i+2] = c2*c2;
                      t3         = (double)i+3;
                      arg3       = hpi*t3/tin2;
                      c3         = std::cos(arg3);
                      Phi0t[i+3] = c3*c3;
              }
	 }


            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void squared_cos_flux_unroll_2x( float * __restrict __ATTR_ALIGN__(64) Phi0t,
	                                     const float Phi0,
					     const int32_t n,
					     const float tin) {

                constexpr float hpi = 1.57079632679489661923132169164f;
		float tin2,t0,t1;
                float arg0,arg1;
                float c0,c1;
                tin2 = 0.5f*tin;
		int32_t i,m,m1;
		m = n%2;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0       = (float)i;
		       arg0     = hpi*t0/tin2;
		       c0       = std::cos(arg0);
		       Phi0t[i] = c0*c0;
		   }
		   if(n<2) return;
		}
		m1 = m+1;
		__assume_aligned(Phi0t,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(4)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	        for(i = m1; i != n; i += 2) {
                      t0         = (float)i;
		      arg0       = hpi*t0/tin2;
                      c0         = std::cos(arg0);
                      Phi0t[i]   = c0*c0;
                      t1         = (float)i+1;
                      arg1       = hpi*t1/tin2;
                      c1         = std::cos(arg1);
                      Phi0t[i+1] = c1*c1;
                                                    
       
		}
	 }


            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void squared_cos_flux_unroll_2x(double * __restrict __ATTR_ALIGN__(64) Phi0t,
	                                     const double Phi0,
					     const int32_t n,
					     const double tin) {

                constexpr double hpi = 1.57079632679489661923132169164;
		double tin2,t0,t1;
                double arg0,arg1;
                double c0,c1;
                tin2 = 0.5f*tin;
		int32_t i,m,m1;
		m = n%2;
		if(m!=0) {
                   for(i = 0; i != m; ++i) {
                       t0       = (double)i;
		       arg0     = hpi*t0/tin2;
		       c0       = std::cos(arg0);
		       Phi0t[i] = c0*c0;
		   }
		   if(n<2) return;
		}
		m1 = m+1;
		__assume_aligned(Phi0t,64);
              #pragma vector aligned
	      #pragma ivdep
	      #pragma vector vectorlength(8)
	      #pragma vector multiple_gather_scatter_by_shuffles
	      #pragma vector always
	        for(i = m1; i != n; i += 2) {
                      t0         = (double)i;
		      arg0       = hpi*t0/tin2;
                      c0         = std::cos(arg0);
                      Phi0t[i]   = c0*c0;
                      t1         = (double)i+1;
                      arg1       = hpi*t1/tin2;
                      c1         = std::cos(arg1);
                      Phi0t[i+1] = c1*c1;
                    
              }
	 }


            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void const_flux_spectr_unroll_16x(float * __restrict __ATTR_ALIGN__(64) Phi0f,
                                              const float Phi0,
                                              const float * __restrict __ATTR_ALIGN__(64) freq,
                                              const int32_t n,
                                              const float T) {

                   constexpr float twopi = 6.283185307179586476925286766559f;
                   float arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7;
                   float arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15;
                   float f0,f1,f2,f3,f4,f5,f6,f7;
                   float f8,f9,f10,f11,f12,f13,f14,f15;
                   float c0,c1,c2,c3,c4,c5,c6,c7;
                   float c8,c9,c10,c11,c12,c13,c14,c15;
                   float sa0,sa1,sa2,sa3,sa4,sa,sa6,sa7;
                   float sa8,sa9,sa10,sa11,sa12,sa13,sa14,sa15;
                   float hT,Phi0T;
                   int32_t i,m,m1;
                   hT    = 0.5f*T;
                   Phi0T = Phi0*hT;
                   m     = n%16;
                   if(m!=0) {
                      for(i = 0; i != m; ++i) {
                           f0         = freq[i];
                           c0         = twopi*f0*hT;
                           sa0        = std::sin(c0)/c0;
                           arg0       = 1.0f-(f0*hT*f0*hT);
                           Phi0f[i]   = Phi0T*(sa0/arg0);
                      }
                      if(m<16) return;
                   }
                   m1 = m+1;
                   __assume_aligned(Phi0f,64);
                   __assume_aligned(freq,64);
                  #pragma vector aligned
	          #pragma ivdep
	          #pragma vector vectorlength(4)
	          #pragma vector multiple_gather_scatter_by_shuffles
	          #pragma vector always
                  for(i = m1; i != n; i += 16) {
                       f0         = freq[i];
                       c0         = twopi*f0*hT;
                       sa0        = std::sin(c0)/c0;
                       arg0       = 1.0f-(f0*hT*f0*hT);
                       Phi0f[i]   = Phi0T*(sa0/arg0);
                       f1         = freq[i+1];
                       c1         = twopi*f1*hT;
                       sa1        = std::sin(c1)/c1;
                       arg1       = 1.0f-(f0*hT*f1*hT);
                       Phi0f[i+1] = Phi0T*(sa1/arg1);
                       f2         = freq[i+2];
                       c2         = twopi*f2*hT;
                       sa2        = std::sin(c2)/c2;
                       arg2       = 1.0f-(f2*hT*f2*hT);
                       Phi0f[i+2] = Phi0T*(sa2/arg2);
                       f3         = freq[i+3];
                       c3         = twopi*f3*hT;
                       sa3        = std::sin(c3)/c3;
                       arg3       = 1.0f-(f3*hT*f3*hT);
                       Phi0f[i+3] = Phi0T*(sa3/arg3);
                       f4         = freq[i+4];
                       c4         = twopi*f4*hT;
                       sa4        = std::sin(c4)/c4;
                       arg4       = 1.0f-(f4*hT*f4*hT);
                       Phi0f[i+4] = Phi0T*(sa4/arg4);
                       f5         = freq[i+5];
                       c5         = twopi*f5*hT;
                       sa5        = std::sin(c5)/c5;
                       arg5       = 1.0f-(f5*hT*f5*hT);
                       Phi0f[i+5] = Phi0T*(sa5/arg5);
                       f6         = freq(i+6);
                       c6         = twopi*f6*hT;
                       sa6        = sin(c6)/c6;
                       arg6       = 1.0f-(f6*hT*f6*hT);
                       Phi0f[i+6] = Phi0T*(sa6/arg6);
                       f7         = freq[i+7];
                       c7         = twopi*f7*hT;
                       sa7        = std::sin(c7)/c7;
                       arg7       = 1.0f-(f7*hT*f7*hT);
                       Phi0f[i+7] = Phi0T*(sa7/arg7);
                       f8         = freq[i+8];
                       c8         = twopi*f8*hT;
                       sa8        = std::sin(c8)/c8;
                       arg8       = 1.0f-(f8*hT*f8*hT);
                       Phi0f[i+8] = Phi0T*(sa8/arg8);
                       f9         = freq[i+9];
                       c9         = twopi*f9*hT;
                       sa9        = sstd::in(c9)/c9;
                       arg9       = 1.0f-(f9*hT*f9*hT);
                       Phi0f[i+9] = Phi0T*(sa9/arg9);
                       f10        = freq[i+10];
                       c10        = twopi*f10*hT;
                       sa10       = std::sin(c10)/c10;
                       arg10      = 1.0f-(f10*hT*f10*hT);
                       Phi0f[i+10]= Phi0T*(sa10/arg10);
                       f11        = freq[i+11];
                       c11        = twopi*f11*hT;
                       sa11       = std::sin(c11)/c11;
                       arg11      = 1.0f-(f11*hT*f11*hT);
                       Phi0f[i+11]= Phi0T*(sa11/arg11);
                       f12        = freq[i+12];
                       c12        = twopi*f12*hT;
                       sa12       = std::sin(c12)/c12;
                       arg12      = 1.0f-(f12*hT*f12*hT);
                       Phi0f[i+12]= Phi0T*(sa12/arg12);
                       f13        = freq(i+13);
                       c13        = twopi*f13*hT;
                       sa13       = std::sin(c13)/c13;
                       arg13      = 1.0f-(f13*hT*f13*hT);
                       Phi0f[i+13]= Phi0T*(sa13/arg13);
                       f14        = freq[i+14];
                       c14        = twopi*f14*hT;
                       sa14       = std::sin(c14)/c14;
                       arg14      = 1.0f-(f14*hT*f14*hT);
                       Phi0f[i+14]= Phi0T*(sa14/arg14);
                       f15        = freq[i+15];
                       c15        = twopi*f15*h;
                       sa15       = std::sin(c15)/c15;
                       arg15      = 1.0f-(f15*hT*f15*hT);
                       Phi0f[i+15]= Phi0T*(sa15/arg15);
                  }
           }


            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void const_flux_spectr_unroll_16x(double * __restrict __ATTR_ALIGN__(64) Phi0f,
                                              const double Phi0,
                                              const double * __restrict __ATTR_ALIGN__(64) freq,
                                              const int32_t n,
                                              const double T) {

                   constexpr double twopi = 6.283185307179586476925286766559;
                   double arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7;
                   double arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15;
                   double f0,f1,f2,f3,f4,f5,f6,f7;
                   double f8,f9,f10,f11,f12,f13,f14,f15;
                   double c0,c1,c2,c3,c4,c5,c6,c7;
                   double c8,c9,c10,c11,c12,c13,c14,c15;
                   double sa0,sa1,sa2,sa3,sa4,sa,sa6,sa7;
                   double sa8,sa9,sa10,sa11,sa12,sa13,sa14,sa15;
                   double hT,Phi0T;
                   int32_t i,m,m1;
                   hT    = 0.5f*T;
                   Phi0T = Phi0*hT;
                   m     = n%16;
                   if(m!=0) {
                      for(i = 0; i != m; ++i) {
                           f0         = freq[i];
                           c0         = twopi*f0*hT;
                           sa0        = std::sin(c0)/c0;
                           arg0       = 1.0-(f0*hT*f0*hT);
                           Phi0f[i]   = Phi0T*(sa0/arg0);
                      }
                      if(m<16) return;
                   }
                   m1 = m+1;
                   __assume_aligned(Phi0f,64);
                   __assume_aligned(freq,64);
                  #pragma vector aligned
	          #pragma ivdep
	          #pragma vector vectorlength(8)
	          #pragma vector multiple_gather_scatter_by_shuffles
	          #pragma vector always
                  for(i = m1; i != n; i += 16) {
                       f0         = freq[i];
                       c0         = twopi*f0*hT;
                       sa0        = std::sin(c0)/c0;
                       arg0       = 1.0-(f0*hT*f0*hT);
                       Phi0f[i]   = Phi0T*(sa0/arg0);
                       f1         = freq[i+1];
                       c1         = twopi*f1*hT;
                       sa1        = std::sin(c1)/c1;
                       arg1       = 1.0-(f0*hT*f1*hT);
                       Phi0f[i+1] = Phi0T*(sa1/arg1);
                       f2         = freq[i+2];
                       c2         = twopi*f2*hT;
                       sa2        = std::sin(c2)/c2;
                       arg2       = 1.0-(f2*hT*f2*hT);
                       Phi0f[i+2] = Phi0T*(sa2/arg2);
                       f3         = freq[i+3];
                       c3         = twopi*f3*hT;
                       sa3        = std::sin(c3)/c3;
                       arg3       = 1.0-(f3*hT*f3*hT);
                       Phi0f[i+3] = Phi0T*(sa3/arg3);
                       f4         = freq[i+4];
                       c4         = twopi*f4*hT;
                       sa4        = std::sin(c4)/c4;
                       arg4       = 1.0-(f4*hT*f4*hT);
                       Phi0f[i+4] = Phi0T*(sa4/arg4);
                       f5         = freq[i+5];
                       c5         = twopi*f5*hT;
                       sa5        = std::sin(c5)/c5;
                       arg5       = 1.0-(f5*hT*f5*hT);
                       Phi0f[i+5] = Phi0T*(sa5/arg5);
                       f6         = freq(i+6);
                       c6         = twopi*f6*hT;
                       sa6        = sin(c6)/c6;
                       arg6       = 1.0-(f6*hT*f6*hT);
                       Phi0f[i+6] = Phi0T*(sa6/arg6);
                       f7         = freq[i+7];
                       c7         = twopi*f7*hT;
                       sa7        = std::sin(c7)/c7;
                       arg7       = 1.0-(f7*hT*f7*hT);
                       Phi0f[i+7] = Phi0T*(sa7/arg7);
                       f8         = freq[i+8];
                       c8         = twopi*f8*hT;
                       sa8        = std::sin(c8)/c8;
                       arg8       = 1.0-(f8*hT*f8*hT);
                       Phi0f[i+8] = Phi0T*(sa8/arg8);
                       f9         = freq[i+9];
                       c9         = twopi*f9*hT;
                       sa9        = sstd::in(c9)/c9;
                       arg9       = 1.0-(f9*hT*f9*hT);
                       Phi0f[i+9] = Phi0T*(sa9/arg9);
                       f10        = freq[i+10];
                       c10        = twopi*f10*hT;
                       sa10       = std::sin(c10)/c10;
                       arg10      = 1.0-(f10*hT*f10*hT);
                       Phi0f[i+10]= Phi0T*(sa10/arg10);
                       f11        = freq[i+11];
                       c11        = twopi*f11*hT;
                       sa11       = std::sin(c11)/c11;
                       arg11      = 1.0-(f11*hT*f11*hT);
                       Phi0f[i+11]= Phi0T*(sa11/arg11);
                       f12        = freq[i+12];
                       c12        = twopi*f12*hT;
                       sa12       = std::sin(c12)/c12;
                       arg12      = 1.0-(f12*hT*f12*hT);
                       Phi0f[i+12]= Phi0T*(sa12/arg12);
                       f13        = freq(i+13);
                       c13        = twopi*f13*hT;
                       sa13       = std::sin(c13)/c13;
                       arg13      = 1.0-(f13*hT*f13*hT);
                       Phi0f[i+13]= Phi0T*(sa13/arg13);
                       f14        = freq[i+14];
                       c14        = twopi*f14*hT;
                       sa14       = std::sin(c14)/c14;
                       arg14      = 1.0-(f14*hT*f14*hT);
                       Phi0f[i+14]= Phi0T*(sa14/arg14);
                       f15        = freq[i+15];
                       c15        = twopi*f15*h;
                       sa15       = std::sin(c15)/c15;
                       arg15      = 1.0-(f15*hT*f15*hT);
                       Phi0f[i+15]= Phi0T*(sa15/arg15);
                  }
           }


           
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void const_flux_spectr_unroll_8x(float * __restrict __ATTR_ALIGN__(64) Phi0f,
                                              const float Phi0,
                                              const float * __restrict __ATTR_ALIGN__(64) freq,
                                              const int32_t n,
                                              const float T) {

                   constexpr float twopi = 6.283185307179586476925286766559f;
                   float arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7;
                   float f0,f1,f2,f3,f4,f5,f6,f7;
                   float c0,c1,c2,c3,c4,c5,c6,c7;
                   float sa0,sa1,sa2,sa3,sa4,sa,sa6,sa7;
                   float hT,Phi0T;
                   int32_t i,m,m1;
                   hT    = 0.5f*T;
                   Phi0T = Phi0*hT;
                   m     = n%8;
                   if(m!=0) {
                      for(i = 0; i != m; ++i) {
                           f0         = freq[i];
                           c0         = twopi*f0*hT;
                           sa0        = std::sin(c0)/c0;
                           arg0       = 1.0f-(f0*hT*f0*hT);
                           Phi0f[i]   = Phi0T*(sa0/arg0);
                      }
                      if(m<8) return;
                   }
                   m1 = m+1;
                   __assume_aligned(Phi0f,64);
                   __assume_aligned(freq,64);
                  #pragma vector aligned
	          #pragma ivdep
	          #pragma vector vectorlength(4)
	          #pragma vector multiple_gather_scatter_by_shuffles
	          #pragma vector always
                  for(i = m1; i != n; i += 8) {
                       f0         = freq[i];
                       c0         = twopi*f0*hT;
                       sa0        = std::sin(c0)/c0;
                       arg0       = 1.0f-(f0*hT*f0*hT);
                       Phi0f[i]   = Phi0T*(sa0/arg0);
                       f1         = freq[i+1];
                       c1         = twopi*f1*hT;
                       sa1        = std::sin(c1)/c1;
                       arg1       = 1.0f-(f0*hT*f1*hT);
                       Phi0f[i+1] = Phi0T*(sa1/arg1);
                       f2         = freq[i+2];
                       c2         = twopi*f2*hT;
                       sa2        = std::sin(c2)/c2;
                       arg2       = 1.0f-(f2*hT*f2*hT);
                       Phi0f[i+2] = Phi0T*(sa2/arg2);
                       f3         = freq[i+3];
                       c3         = twopi*f3*hT;
                       sa3        = std::sin(c3)/c3;
                       arg3       = 1.0f-(f3*hT*f3*hT);
                       Phi0f[i+3] = Phi0T*(sa3/arg3);
                       f4         = freq[i+4];
                       c4         = twopi*f4*hT;
                       sa4        = std::sin(c4)/c4;
                       arg4       = 1.0f-(f4*hT*f4*hT);
                       Phi0f[i+4] = Phi0T*(sa4/arg4);
                       f5         = freq[i+5];
                       c5         = twopi*f5*hT;
                       sa5        = std::sin(c5)/c5;
                       arg5       = 1.0f-(f5*hT*f5*hT);
                       Phi0f[i+5] = Phi0T*(sa5/arg5);
                       f6         = freq(i+6);
                       c6         = twopi*f6*hT;
                       sa6        = sin(c6)/c6;
                       arg6       = 1.0f-(f6*hT*f6*hT);
                       Phi0f[i+6] = Phi0T*(sa6/arg6);
                       f7         = freq[i+7];
                       c7         = twopi*f7*hT;
                       sa7        = std::sin(c7)/c7;
                       arg7       = 1.0f-(f7*hT*f7*hT);
                       Phi0f[i+7] = Phi0T*(sa7/arg7);
                     
                  }
           }


           
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void const_flux_spectr_unroll_8x( double * __restrict __ATTR_ALIGN__(64) Phi0f,
                                              const double Phi0,
                                              const double * __restrict __ATTR_ALIGN__(64) freq,
                                              const int32_t n,
                                              const double T) {

                   constexpr double twopi = 6.283185307179586476925286766559;
                   double arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7;
                   double f0,f1,f2,f3,f4,f5,f6,f7;
                   double c0,c1,c2,c3,c4,c5,c6,c7;
                   double sa0,sa1,sa2,sa3,sa4,sa,sa6,sa7;
                   double hT,Phi0T;
                   int32_t i,m,m1;
                   hT    = 0.5*T;
                   Phi0T = Phi0*hT;
                   m     = n%8;
                   if(m!=0) {
                      for(i = 0; i != m; ++i) {
                           f0         = freq[i];
                           c0         = twopi*f0*hT;
                           sa0        = std::sin(c0)/c0;
                           arg0       = 1.0-(f0*hT*f0*hT);
                           Phi0f[i]   = Phi0T*(sa0/arg0);
                      }
                      if(m<8) return;
                   }
                   m1 = m+1;
                   __assume_aligned(Phi0f,64);
                   __assume_aligned(freq,64);
                  #pragma vector aligned
	          #pragma ivdep
	          #pragma vector vectorlength(8)
	          #pragma vector multiple_gather_scatter_by_shuffles
	          #pragma vector always
                  for(i = m1; i != n; i += 8) {
                       f0         = freq[i];
                       c0         = twopi*f0*hT;
                       sa0        = std::sin(c0)/c0;
                       arg0       = 1.0-(f0*hT*f0*hT);
                       Phi0f[i]   = Phi0T*(sa0/arg0);
                       f1         = freq[i+1];
                       c1         = twopi*f1*hT;
                       sa1        = std::sin(c1)/c1;
                       arg1       = 1.0-(f0*hT*f1*hT);
                       Phi0f[i+1] = Phi0T*(sa1/arg1);
                       f2         = freq[i+2];
                       c2         = twopi*f2*hT;
                       sa2        = std::sin(c2)/c2;
                       arg2       = 1.0-(f2*hT*f2*hT);
                       Phi0f[i+2] = Phi0T*(sa2/arg2);
                       f3         = freq[i+3];
                       c3         = twopi*f3*hT;
                       sa3        = std::sin(c3)/c3;
                       arg3       = 1.0-(f3*hT*f3*hT);
                       Phi0f[i+3] = Phi0T*(sa3/arg3);
                       f4         = freq[i+4];
                       c4         = twopi*f4*hT;
                       sa4        = std::sin(c4)/c4;
                       arg4       = 1.0-(f4*hT*f4*hT);
                       Phi0f[i+4] = Phi0T*(sa4/arg4);
                       f5         = freq[i+5];
                       c5         = twopi*f5*hT;
                       sa5        = std::sin(c5)/c5;
                       arg5       = 1.0-(f5*hT*f5*hT);
                       Phi0f[i+5] = Phi0T*(sa5/arg5);
                       f6         = freq(i+6);
                       c6         = twopi*f6*hT;
                       sa6        = sin(c6)/c6;
                       arg6       = 1.0-(f6*hT*f6*hT);
                       Phi0f[i+6] = Phi0T*(sa6/arg6);
                       f7         = freq[i+7];
                       c7         = twopi*f7*hT;
                       sa7        = std::sin(c7)/c7;
                       arg7       = 1.0-(f7*hT*f7*hT);
                       Phi0f[i+7] = Phi0T*(sa7/arg7);
                      
                  }
           }



            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void const_flux_spectr_unroll_4x( float * __restrict __ATTR_ALIGN__(64) Phi0f,
                                              const float Phi0,
                                              const float * __restrict __ATTR_ALIGN__(64) freq,
                                              const int32_t n,
                                              const float T) {

                   constexpr float twopi = 6.283185307179586476925286766559f;
                   float arg0,arg1,arg2,arg3;
                   float f0,f1,f2,f3;
                   float c0,c1,c2,c3,c4,c5,c6,c7;
                   float sa0,sa1,sa2,sa3;
                   float hT,Phi0T;
                   int32_t i,m,m1;
                   hT    = 0.5f*T;
                   Phi0T = Phi0*hT;
                   m     = n%4;
                   if(m!=0) {
                      for(i = 0; i != m; ++i) {
                           f0         = freq[i];
                           c0         = twopi*f0*hT;
                           sa0        = std::sin(c0)/c0;
                           arg0       = 1.0f-(f0*hT*f0*hT);
                           Phi0f[i]   = Phi0T*(sa0/arg0);
                      }
                      if(m<4) return;
                   }
                   m1 = m+1;
                   __assume_aligned(Phi0f,64);
                   __assume_aligned(freq,64);
                  #pragma vector aligned
	          #pragma ivdep
	          #pragma vector vectorlength(4)
	          #pragma vector multiple_gather_scatter_by_shuffles
	          #pragma vector always
                  for(i = m1; i != n; i += 4) {
                       f0         = freq[i];
                       c0         = twopi*f0*hT;
                       sa0        = std::sin(c0)/c0;
                       arg0       = 1.0f-(f0*hT*f0*hT);
                       Phi0f[i]   = Phi0T*(sa0/arg0);
                       f1         = freq[i+1];
                       c1         = twopi*f1*hT;
                       sa1        = std::sin(c1)/c1;
                       arg1       = 1.0f-(f0*hT*f1*hT);
                       Phi0f[i+1] = Phi0T*(sa1/arg1);
                       f2         = freq[i+2];
                       c2         = twopi*f2*hT;
                       sa2        = std::sin(c2)/c2;
                       arg2       = 1.0f-(f2*hT*f2*hT);
                       Phi0f[i+2] = Phi0T*(sa2/arg2);
                       f3         = freq[i+3];
                       c3         = twopi*f3*hT;
                       sa3        = std::sin(c3)/c3;
                       arg3       = 1.0f-(f3*hT*f3*hT);
                       Phi0f[i+3] = Phi0T*(sa3/arg3);
                                           
                  }
           }


            
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void const_flux_spectr_unroll_4x( double * __restrict __ATTR_ALIGN__(64) Phi0f,
                                              const double Phi0,
                                              const double * __restrict __ATTR_ALIGN__(64) freq,
                                              const int32_t n,
                                              const double T) {

                   constexpr double twopi = 6.283185307179586476925286766559;
                   double arg0,arg1,arg2,arg3;
                   double f0,f1,f2,f3;
                   double c0,c1,c2,c3;
                   double sa0,sa1,sa2,sa3;
                   double hT,Phi0T;
                   int32_t i,m,m1;
                   hT    = 0.5*T;
                   Phi0T = Phi0*hT;
                   m     = n%4;
                   if(m!=0) {
                      for(i = 0; i != m; ++i) {
                           f0         = freq[i];
                           c0         = twopi*f0*hT;
                           sa0        = std::sin(c0)/c0;
                           arg0       = 1.0-(f0*hT*f0*hT);
                           Phi0f[i]   = Phi0T*(sa0/arg0);
                      }
                      if(m<4) return;
                   }
                   m1 = m+1;
                   __assume_aligned(Phi0f,64);
                   __assume_aligned(freq,64);
                  #pragma vector aligned
	          #pragma ivdep
	          #pragma vector vectorlength(8)
	          #pragma vector multiple_gather_scatter_by_shuffles
	          #pragma vector always
                  for(i = m1; i != n; i += 4) {
                       f0         = freq[i];
                       c0         = twopi*f0*hT;
                       sa0        = std::sin(c0)/c0;
                       arg0       = 1.0-(f0*hT*f0*hT);
                       Phi0f[i]   = Phi0T*(sa0/arg0);
                       f1         = freq[i+1];
                       c1         = twopi*f1*hT;
                       sa1        = std::sin(c1)/c1;
                       arg1       = 1.0-(f0*hT*f1*hT);
                       Phi0f[i+1] = Phi0T*(sa1/arg1);
                       f2         = freq[i+2];
                       c2         = twopi*f2*hT;
                       sa2        = std::sin(c2)/c2;
                       arg2       = 1.0-(f2*hT*f2*hT);
                       Phi0f[i+2] = Phi0T*(sa2/arg2);
                       f3         = freq[i+3];
                       c3         = twopi*f3*hT;
                       sa3        = std::sin(c3)/c3;
                       arg3       = 1.0-(f3*hT*f3*hT);
                       Phi0f[i+3] = Phi0T*(sa3/arg3);
                 }
           }


            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void const_flux_spectr_unroll_2x( float * __restrict __ATTR_ALIGN__(64) Phi0f,
                                              const float Phi0,
                                              const float * __restrict __ATTR_ALIGN__(64) freq,
                                              const int32_t n,
                                              const float T) {

                   constexpr float twopi = 6.283185307179586476925286766559f;
                   float arg0,arg1;
                   float f0,f1;
                   float c0,c1;
                   float sa0,sa1;
                   float hT,Phi0T;
                   int32_t i,m,m1;
                   hT    = 0.5f*T;
                   Phi0T = Phi0*hT;
                   m     = n%2;
                   if(m!=0) {
                      for(i = 0; i != m; ++i) {
                           f0         = freq[i];
                           c0         = twopi*f0*hT;
                           sa0        = std::sin(c0)/c0;
                           arg0       = 1.0f-(f0*hT*f0*hT);
                           Phi0f[i]   = Phi0T*(sa0/arg0);
                      }
                      if(m<2) return;
                   }
                   m1 = m+1;
                   __assume_aligned(Phi0f,64);
                   __assume_aligned(freq,64);
                  #pragma vector aligned
	          #pragma ivdep
	          #pragma vector vectorlength(4)
	          #pragma vector multiple_gather_scatter_by_shuffles
	          #pragma vector always
                  for(i = m1; i != n; i += 2) {
                       f0         = freq[i];
                       c0         = twopi*f0*hT;
                       sa0        = std::sin(c0)/c0;
                       arg0       = 1.0f-(f0*hT*f0*hT);
                       Phi0f[i]   = Phi0T*(sa0/arg0);
                       f1         = freq[i+1];
                       c1         = twopi*f1*hT;
                       sa1        = std::sin(c1)/c1;
                       arg1       = 1.0f-(f0*hT*f1*hT);
                       Phi0f[i+1] = Phi0T*(sa1/arg1);
                                                                 
                  }
           }


            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void const_flux_spectr_unroll_2x( double * __restrict __ATTR_ALIGN__(64) Phi0f,
                                              const double Phi0,
                                              const double * __restrict __ATTR_ALIGN__(64) freq,
                                              const int32_t n,
                                              const double T) {

                   constexpr double twopi = 6.283185307179586476925286766559;
                   double arg0,arg1;
                   double f0,f1;
                   double c0,c1;
                   double sa0,sa1;
                   double hT,Phi0T;
                   int32_t i,m,m1;
                   hT    = 0.5*T;
                   Phi0T = Phi0*hT;
                   m     = n%2;
                   if(m!=0) {
                      for(i = 0; i != m; ++i) {
                           f0         = freq[i];
                           c0         = twopi*f0*hT;
                           sa0        = std::sin(c0)/c0;
                           arg0       = 1.0-(f0*hT*f0*hT);
                           Phi0f[i]   = Phi0T*(sa0/arg0);
                      }
                      if(m<2) return;
                   }
                   m1 = m+1;
                   __assume_aligned(Phi0f,64);
                   __assume_aligned(freq,64);
                  #pragma vector aligned
	          #pragma ivdep
	          #pragma vector vectorlength(8)
	          #pragma vector multiple_gather_scatter_by_shuffles
	          #pragma vector always
                  for(i = m1; i != n; i += 2) {
                       f0         = freq[i];
                       c0         = twopi*f0*hT;
                       sa0        = std::sin(c0)/c0;
                       arg0       = 1.0-(f0*hT*f0*hT);
                       Phi0f[i]   = Phi0T*(sa0/arg0);
                       f1         = freq[i+1];
                       c1         = twopi*f1*hT;
                       sa1        = std::sin(c1)/c1;
                       arg1       = 1.0-(f0*hT*f1*hT);
                       Phi0f[i+1] = Phi0T*(sa1/arg1);
                 }
           }


        //  !Идеальный гармонический модулятор
        //!Formula 1,2 p. 187
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void ideal_modulator_unroll_16x(float * __restrict __ATTR_ALIGN__(64) rhot_s,
                                            float * __restrict __ATTR_ALIGN__(64) rhot_c,
                                            const int32_t n,
                                            const float f0,
                                            const float phi0,
                                            const float rho0,
                                            const float rho1) {

                 constexpr float twopi = 6.283185307179586476925286766559f;
                 float t0,t1,t2,t3,t4,t5,t6,t7;
                 float t8,t9,t10,t11,t12,t13,t14,t15;
                 float s0,s1,s2,s3,s4,s5,s6,s7;
                 float s8,s9,s10,s11,s12,s13,s14,s15;
                 float c0,c1,c2,c3,c4,c5,c6,c7;
                 float c8,c9,c10,c11,c12,c13,c14,c15;
                 float psi0,psi1,psi2,psi3,psi4,psi5,psi6,psi7;
                 float psi8,psi9,psi10,psi11,psi12,psi13,psi14,psi15;
                 float om0;
                 int32_t i,m,m1;
                 om0 = twopi*f0;
                 m   = n%16;
                 if(m!=0) {
                    for(i = 0; i != m; ++i) {
                        t0        = (float)i;
                        psi0      = om0*t0+phi0;
                        s0        = rho0+rho1*std::sin(psi0);
                        rhot_s[i] = s0
                        c0        = rho0+rho1*std::cos(psi0);
                        rhot_c[i] = c0;
                    }
                    if(n<16) return;
                 }
                 m1 = m+1;
                 __assume_aligned(rhot_s,64);
                 __assume_aligned(rhot_c,64);
                 #pragma vector aligned
	         #pragma ivdep
	         #pragma vector vectorlength(4)
	         #pragma vector multiple_gather_scatter_by_shuffles
	         #pragma vector always
                 for(i = m1; i != n; i += 16) {
                      t0          = (float)i;
                      psi0        = om0*t0+phi0;
                      s0          = rho0+rho1*std::sin(psi0);
                      rhot_s[i]   = s0;
                      c0          = rho0+rho1*std::cos(psi0);
                      rhot_c[i]   = c0;
                      t1          = (float)i+1;
                      psi1        = om0*t1+phi0;
                      s1          = rho0+rho1*std::sin(psi1);
                      rhot_s[i+1] = s1;
                      c0          = rho0+rho1*std::cos(psi1);
                      rhot_c[i+1] = c1;
                      t2          = (float)i+2;
                      psi2        = om0*t2+phi0;
                      s2          = rho0+rho1*std::sin(psi2);
                      rhot_s[i+2] = s2;
                      c2          = rho0+rho1*std::cos(psi2);
                      rhot_c[i+2] = c2;
                      t3          = (float)i+3;
                      psi3        = om0*t3+phi0;
                      s3          = rho0+rho1*std::sin(psi3);
                      rhot_s[i+3] = s3;
                      c3          = rho0+rho1*std::cos(psi3);
                      rhot_c[i+3] = c3;
                      t4          = (float)i+4;
                      psi4        = om0*t4+phi0;
                      s4          = rho0+rho1*std::sin(psi4);
                      rhot_s[i+4] = s4;
                      c4          = rho0+rho1*std::cos(psi4);
                      rhot_c[i+4] = c4;
                      t5          = (float)i+5;
                      psi5        = om0*t5+phi0;
                      s5          = rho0+rho1*std::sin(psi5);
                      rhot_s[i+5] = s5;
                      c5          = rho0+rho1*std::cos(psi5)
                      rhot_c[i+5] = c5
                      t6          = (float)i+6;
                      psi6        = om0*t6+phi0
                      s6          = rho0+rho1*std::sin(psi6);
                      rhot_s[i+6] = s6;
                      c6          = rho0+rho1*cos(psi6);
                      rhot_c[i+6] = c6;
                      t7          = (float)i+7;
                      psi7        = om0*t7+phi0;
                      s7          = rho0+rho1*std::sin(psi7);
                      rhot_s[i+7] = s7;
                      c7          = rho0+rho1*std::cos(psi7);
                      rhot_c[i+7] = c7;
                      t8          = (float)i+8;
                      psi8        = om0*t8+phi0;
                      s8          = rho0+rho1*std::sin(psi8);
                      rhot_s[i+8] = s8;
                      c8          = rho0+rho1*std::cos(psi8);
                      rhot_c[i+8] = c8;
                      t9          = (float)i+9;
                      psi9        = om0*t9+phi0;
                      s9          = rho0+rho1*std::sin(psi9);
                      rhot_s[i+9] = s9;
                      c9          = rho0+rho1*std::cos(psi9);
                      rhot_c[i+9] = c9;
                      t10         = (float)i+10;
                      psi10       = om0*t10+phi0;
                      s10         = rho0+rho1*std::sin(psi10);
                      rhot_s[i+10]= s10;
                      c10         = rho0+rho1*std::cos(psi10);
                      rhot_c[i+10]= c10;
                      t11         = (float)i+11;
                      psi11       = om0*t11+phi0;
                      s11         = rho0+rho1*std::sin(psi11);
                      rhot_s[i+11]= s11;
                      c11         = rho0+rho1*std::cos(psi11);
                      rhot_c[i+11]= c11;
                      t12         = (float)i+12;
                      psi12       = om0*t12+phi0;
                      s12         = rho0+rho1*std::sin(psi12);
                      rhot_s[i+12]= s12;
                      c12         = rho0+rho1*std::cos(psi12);
                      rhot_c[i+12]= c12;
                      t13         = (float)i+13;
                      psi13       = om0*t13+phi0;
                      s13         = rho0+rho1*std::sin(psi13);
                      rhot_s[i+13]= s13;
                      c13         = rho0+rho1*std::cos(psi13);
                      rhot_c[i+13]= c13;
                      t14         = (float)i+14;
                      psi14       = om0*t14+phi0;
                      s14         = rho0+rho1*std::sin(psi14);
                      rhot_s[i+14]= s14;
                      c14         = rho0+rho1*std::cos(psi14);
                      rhot_c[i+14]= c14;
                      t15         = (float)i+15;
                      psi15       = om0*t15+phi0;
                      s15         = rho0+rho1*std::sin(psi15)
                      rhot_s[i+15]= s15;
                      c15         = rho0+rho1*std::cos(psi15)
                      rhot_c[i+15]= c15;
                 }
        }


   
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void ideal_modulator_unroll_16x(double * __restrict __ATTR_ALIGN__(64) rhot_s,
                                            double * __restrict __ATTR_ALIGN__(64) rhot_c,
                                            const int32_t n,
                                            const double f0,
                                            const double phi0,
                                            const double rho0,
                                            const double rho1) {

                 constexpr double twopi = 6.283185307179586476925286766559;
                 double t0,t1,t2,t3,t4,t5,t6,t7;
                 double t8,t9,t10,t11,t12,t13,t14,t15;
                 double s0,s1,s2,s3,s4,s5,s6,s7;
                 double s8,s9,s10,s11,s12,s13,s14,s15;
                 double c0,c1,c2,c3,c4,c5,c6,c7;
                 double c8,c9,c10,c11,c12,c13,c14,c15;
                 double psi0,psi1,psi2,psi3,psi4,psi5,psi6,psi7;
                 double psi8,psi9,psi10,psi11,psi12,psi13,psi14,psi15;
                 double om0;
                 int32_t i,m,m1;
                 om0 = twopi*f0;
                 m   = n%16;
                 if(m!=0) {
                    for(i = 0; i != m; ++i) {
                        t0        = (double)i;
                        psi0      = om0*t0+phi0;
                        s0        = rho0+rho1*std::sin(psi0);
                        rhot_s[i] = s0
                        c0        = rho0+rho1*std::cos(psi0);
                        rhot_c[i] = c0;
                    }
                    if(n<16) return;
                 }
                 m1 = m+1;
                 __assume_aligned(rhot_s,64);
                 __assume_aligned(rhot_c,64);
                 #pragma vector aligned
	         #pragma ivdep
	         #pragma vector vectorlength(8)
	         #pragma vector multiple_gather_scatter_by_shuffles
	         #pragma vector always
                 for(i = m1; i != n; i += 16) {
                      t0          = (double)i;
                      psi0        = om0*t0+phi0;
                      s0          = rho0+rho1*std::sin(psi0);
                      rhot_s[i]   = s0;
                      c0          = rho0+rho1*std::cos(psi0);
                      rhot_c[i]   = c0;
                      t1          = (double)i+1;
                      psi1        = om0*t1+phi0;
                      s1          = rho0+rho1*std::sin(psi1);
                      rhot_s[i+1] = s1;
                      c0          = rho0+rho1*std::cos(psi1);
                      rhot_c[i+1] = c1;
                      t2          = (double)i+2;
                      psi2        = om0*t2+phi0;
                      s2          = rho0+rho1*std::sin(psi2);
                      rhot_s[i+2] = s2;
                      c2          = rho0+rho1*std::cos(psi2);
                      rhot_c[i+2] = c2;
                      t3          = (double)i+3;
                      psi3        = om0*t3+phi0;
                      s3          = rho0+rho1*std::sin(psi3);
                      rhot_s[i+3] = s3;
                      c3          = rho0+rho1*std::cos(psi3);
                      rhot_c[i+3] = c3;
                      t4          = (double)i+4;
                      psi4        = om0*t4+phi0;
                      s4          = rho0+rho1*std::sin(psi4);
                      rhot_s[i+4] = s4;
                      c4          = rho0+rho1*std::cos(psi4);
                      rhot_c[i+4] = c4;
                      t5          = (double)i+5;
                      psi5        = om0*t5+phi0;
                      s5          = rho0+rho1*std::sin(psi5);
                      rhot_s[i+5] = s5;
                      c5          = rho0+rho1*std::cos(psi5)
                      rhot_c[i+5] = c5
                      t6          = (double)i+6;
                      psi6        = om0*t6+phi0
                      s6          = rho0+rho1*std::sin(psi6);
                      rhot_s[i+6] = s6;
                      c6          = rho0+rho1*cos(psi6);
                      rhot_c[i+6] = c6;
                      t7          = (double)i+7;
                      psi7        = om0*t7+phi0;
                      s7          = rho0+rho1*std::sin(psi7);
                      rhot_s[i+7] = s7;
                      c7          = rho0+rho1*std::cos(psi7);
                      rhot_c[i+7] = c7;
                      t8          = (double)i+8;
                      psi8        = om0*t8+phi0;
                      s8          = rho0+rho1*std::sin(psi8);
                      rhot_s[i+8] = s8;
                      c8          = rho0+rho1*std::cos(psi8);
                      rhot_c[i+8] = c8;
                      t9          = (double)i+9;
                      psi9        = om0*t9+phi0;
                      s9          = rho0+rho1*std::sin(psi9);
                      rhot_s[i+9] = s9;
                      c9          = rho0+rho1*std::cos(psi9);
                      rhot_c[i+9] = c9;
                      t10         = (double)i+10;
                      psi10       = om0*t10+phi0;
                      s10         = rho0+rho1*std::sin(psi10);
                      rhot_s[i+10]= s10;
                      c10         = rho0+rho1*std::cos(psi10);
                      rhot_c[i+10]= c10;
                      t11         = (double)i+11;
                      psi11       = om0*t11+phi0;
                      s11         = rho0+rho1*std::sin(psi11);
                      rhot_s[i+11]= s11;
                      c11         = rho0+rho1*std::cos(psi11);
                      rhot_c[i+11]= c11;
                      t12         = (double)i+12;
                      psi12       = om0*t12+phi0;
                      s12         = rho0+rho1*std::sin(psi12);
                      rhot_s[i+12]= s12;
                      c12         = rho0+rho1*std::cos(psi12);
                      rhot_c[i+12]= c12;
                      t13         = (double)i+13;
                      psi13       = om0*t13+phi0;
                      s13         = rho0+rho1*std::sin(psi13);
                      rhot_s[i+13]= s13;
                      c13         = rho0+rho1*std::cos(psi13);
                      rhot_c[i+13]= c13;
                      t14         = (double)i+14;
                      psi14       = om0*t14+phi0;
                      s14         = rho0+rho1*std::sin(psi14);
                      rhot_s[i+14]= s14;
                      c14         = rho0+rho1*std::cos(psi14);
                      rhot_c[i+14]= c14;
                      t15         = (double)i+15;
                      psi15       = om0*t15+phi0;
                      s15         = rho0+rho1*std::sin(psi15)
                      rhot_s[i+15]= s15;
                      c15         = rho0+rho1*std::cos(psi15)
                      rhot_c[i+15]= c15;
                 }
        }



            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void ideal_modulator_unroll_8x(float * __restrict __ATTR_ALIGN__(64) rhot_s,
                                            float * __restrict __ATTR_ALIGN__(64) rhot_c,
                                            const int32_t n,
                                            const float f0,
                                            const float phi0,
                                            const float rho0,
                                            const float rho1) {

                 constexpr float twopi = 6.283185307179586476925286766559f;
                 float t0,t1,t2,t3,t4,t5,t6,t7;
                 float s0,s1,s2,s3,s4,s5,s6,s7;
                 float c0,c1,c2,c3,c4,c5,c6,c7;
                 float psi0,psi1,psi2,psi3,psi4,psi5,psi6,psi7;
                 float om0;
                 int32_t i,m,m1;
                 om0 = twopi*f0;
                 m   = n%8;
                 if(m!=0) {
                    for(i = 0; i != m; ++i) {
                        t0        = (float)i;
                        psi0      = om0*t0+phi0;
                        s0        = rho0+rho1*std::sin(psi0);
                        rhot_s[i] = s0
                        c0        = rho0+rho1*std::cos(psi0);
                        rhot_c[i] = c0;
                    }
                    if(n<8) return;
                 }
                 m1 = m+1;
                 __assume_aligned(rhot_s,64);
                 __assume_aligned(rhot_c,64);
                 #pragma vector aligned
	         #pragma ivdep
	         #pragma vector vectorlength(4)
	         #pragma vector multiple_gather_scatter_by_shuffles
	         #pragma vector always
                 for(i = m1; i != n; i += 8) {
                      t0          = (float)i;
                      psi0        = om0*t0+phi0;
                      s0          = rho0+rho1*std::sin(psi0);
                      rhot_s[i]   = s0;
                      c0          = rho0+rho1*std::cos(psi0);
                      rhot_c[i]   = c0;
                      t1          = (float)i+1;
                      psi1        = om0*t1+phi0;
                      s1          = rho0+rho1*std::sin(psi1);
                      rhot_s[i+1] = s1;
                      c0          = rho0+rho1*std::cos(psi1);
                      rhot_c[i+1] = c1;
                      t2          = (float)i+2;
                      psi2        = om0*t2+phi0;
                      s2          = rho0+rho1*std::sin(psi2);
                      rhot_s[i+2] = s2;
                      c2          = rho0+rho1*std::cos(psi2);
                      rhot_c[i+2] = c2;
                      t3          = (float)i+3;
                      psi3        = om0*t3+phi0;
                      s3          = rho0+rho1*std::sin(psi3);
                      rhot_s[i+3] = s3;
                      c3          = rho0+rho1*std::cos(psi3);
                      rhot_c[i+3] = c3;
                      t4          = (float)i+4;
                      psi4        = om0*t4+phi0;
                      s4          = rho0+rho1*std::sin(psi4);
                      rhot_s[i+4] = s4;
                      c4          = rho0+rho1*std::cos(psi4);
                      rhot_c[i+4] = c4;
                      t5          = (float)i+5;
                      psi5        = om0*t5+phi0;
                      s5          = rho0+rho1*std::sin(psi5);
                      rhot_s[i+5] = s5;
                      c5          = rho0+rho1*std::cos(psi5)
                      rhot_c[i+5] = c5
                      t6          = (float)i+6;
                      psi6        = om0*t6+phi0
                      s6          = rho0+rho1*std::sin(psi6);
                      rhot_s[i+6] = s6;
                      c6          = rho0+rho1*cos(psi6);
                      rhot_c[i+6] = c6;
                      t7          = (float)i+7;
                      psi7        = om0*t7+phi0;
                      s7          = rho0+rho1*std::sin(psi7);
                      rhot_s[i+7] = s7;
                      c7          = rho0+rho1*std::cos(psi7);
                      rhot_c[i+7] = c7;
                 }
        }



          


            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void ideal_modulator_unroll_8x( double * __restrict __ATTR_ALIGN__(64) rhot_s,
                                            double * __restrict __ATTR_ALIGN__(64) rhot_c,
                                            const int32_t n,
                                            const double f0,
                                            const double phi0,
                                            const double rho0,
                                            const double rho1) {

                 constexpr double twopi = 6.283185307179586476925286766559;
                 double t0,t1,t2,t3,t4,t5,t6,t7;
                 double s0,s1,s2,s3,s4,s5,s6,s7;
                 double c0,c1,c2,c3,c4,c5,c6,c7;
                 double psi0,psi1,psi2,psi3,psi4,psi5,psi6,psi7;
                 double om0;
                 int32_t i,m,m1;
                 om0 = twopi*f0;
                 m   = n%8;
                 if(m!=0) {
                    for(i = 0; i != m; ++i) {
                        t0        = (double)i;
                        psi0      = om0*t0+phi0;
                        s0        = rho0+rho1*std::sin(psi0);
                        rhot_s[i] = s0
                        c0        = rho0+rho1*std::cos(psi0);
                        rhot_c[i] = c0;
                    }
                    if(n<8) return;
                 }
                 m1 = m+1;
                 __assume_aligned(rhot_s,64);
                 __assume_aligned(rhot_c,64);
                 #pragma vector aligned
	         #pragma ivdep
	         #pragma vector vectorlength(8)
	         #pragma vector multiple_gather_scatter_by_shuffles
	         #pragma vector always
                 for(i = m1; i != n; i += 8) {
                      t0          = (double)i;
                      psi0        = om0*t0+phi0;
                      s0          = rho0+rho1*std::sin(psi0);
                      rhot_s[i]   = s0;
                      c0          = rho0+rho1*std::cos(psi0);
                      rhot_c[i]   = c0;
                      t1          = (double)i+1;
                      psi1        = om0*t1+phi0;
                      s1          = rho0+rho1*std::sin(psi1);
                      rhot_s[i+1] = s1;
                      c0          = rho0+rho1*std::cos(psi1);
                      rhot_c[i+1] = c1;
                      t2          = (double)i+2;
                      psi2        = om0*t2+phi0;
                      s2          = rho0+rho1*std::sin(psi2);
                      rhot_s[i+2] = s2;
                      c2          = rho0+rho1*std::cos(psi2);
                      rhot_c[i+2] = c2;
                      t3          = (double)i+3;
                      psi3        = om0*t3+phi0;
                      s3          = rho0+rho1*std::sin(psi3);
                      rhot_s[i+3] = s3;
                      c3          = rho0+rho1*std::cos(psi3);
                      rhot_c[i+3] = c3;
                      t4          = (double)i+4;
                      psi4        = om0*t4+phi0;
                      s4          = rho0+rho1*std::sin(psi4);
                      rhot_s[i+4] = s4;
                      c4          = rho0+rho1*std::cos(psi4);
                      rhot_c[i+4] = c4;
                      t5          = (double)i+5;
                      psi5        = om0*t5+phi0;
                      s5          = rho0+rho1*std::sin(psi5);
                      rhot_s[i+5] = s5;
                      c5          = rho0+rho1*std::cos(psi5)
                      rhot_c[i+5] = c5
                      t6          = (double)i+6;
                      psi6        = om0*t6+phi0
                      s6          = rho0+rho1*std::sin(psi6);
                      rhot_s[i+6] = s6;
                      c6          = rho0+rho1*cos(psi6);
                      rhot_c[i+6] = c6;
                      t7          = (double)i+7;
                      psi7        = om0*t7+phi0;
                      s7          = rho0+rho1*std::sin(psi7);
                      rhot_s[i+7] = s7;
                      c7          = rho0+rho1*std::cos(psi7);
                      rhot_c[i+7] = c7;
                   
                 }
        }


        
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void ideal_modulator_unroll_4x( float * __restrict __ATTR_ALIGN__(64) rhot_s,
                                            float * __restrict __ATTR_ALIGN__(64) rhot_c,
                                            const int32_t n,
                                            const float f0,
                                            const float phi0,
                                            const float rho0,
                                            const float rho1) {

                 constexpr float twopi = 6.283185307179586476925286766559f;
                 float t0,t1,t2,t3;
                 float s0,s1,s2,s3;
                 float c0,c1,c2,c3;
                 float psi0,psi1,psi2,psi3;
                 float om0;
                 int32_t i,m,m1;
                 om0 = twopi*f0;
                 m   = n%4;
                 if(m!=0) {
                    for(i = 0; i != m; ++i) {
                        t0        = (float)i;
                        psi0      = om0*t0+phi0;
                        s0        = rho0+rho1*std::sin(psi0);
                        rhot_s[i] = s0
                        c0        = rho0+rho1*std::cos(psi0);
                        rhot_c[i] = c0;
                    }
                    if(n<4) return;
                 }
                 m1 = m+1;
                 __assume_aligned(rhot_s,64);
                 __assume_aligned(rhot_c,64);
                 #pragma vector aligned
	         #pragma ivdep
	         #pragma vector vectorlength(4)
	         #pragma vector multiple_gather_scatter_by_shuffles
	         #pragma vector always
                 for(i = m1; i != n; i += 4) {
                      t0          = (float)i;
                      psi0        = om0*t0+phi0;
                      s0          = rho0+rho1*std::sin(psi0);
                      rhot_s[i]   = s0;
                      c0          = rho0+rho1*std::cos(psi0);
                      rhot_c[i]   = c0;
                      t1          = (float)i+1;
                      psi1        = om0*t1+phi0;
                      s1          = rho0+rho1*std::sin(psi1);
                      rhot_s[i+1] = s1;
                      c0          = rho0+rho1*std::cos(psi1);
                      rhot_c[i+1] = c1;
                      t2          = (float)i+2;
                      psi2        = om0*t2+phi0;
                      s2          = rho0+rho1*std::sin(psi2);
                      rhot_s[i+2] = s2;
                      c2          = rho0+rho1*std::cos(psi2);
                      rhot_c[i+2] = c2;
                      t3          = (float)i+3;
                      psi3        = om0*t3+phi0;
                      s3          = rho0+rho1*std::sin(psi3);
                      rhot_s[i+3] = s3;
                      c3          = rho0+rho1*std::cos(psi3);
                      rhot_c[i+3] = c3;
                     
                 }
        }



            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void ideal_modulator_unroll_4x( double * __restrict __ATTR_ALIGN__(64) rhot_s,
                                            double * __restrict __ATTR_ALIGN__(64) rhot_c,
                                            const int32_t n,
                                            const double f0,
                                            const double phi0,
                                            const double rho0,
                                            const double rho1) {

                 constexpr double twopi = 6.283185307179586476925286766559;
                 double t0,t1,t2,t3;
                 double s0,s1,s2,s3;
                 double c0,c1,c2,c3;
                 double psi0,psi1,psi2,psi3;
                 double om0;
                 int32_t i,m,m1;
                 om0 = twopi*f0;
                 m   = n%4;
                 if(m!=0) {
                    for(i = 0; i != m; ++i) {
                        t0        = (double)i;
                        psi0      = om0*t0+phi0;
                        s0        = rho0+rho1*std::sin(psi0);
                        rhot_s[i] = s0
                        c0        = rho0+rho1*std::cos(psi0);
                        rhot_c[i] = c0;
                    }
                    if(n<4) return;
                 }
                 m1 = m+1;
                 __assume_aligned(rhot_s,64);
                 __assume_aligned(rhot_c,64);
                 #pragma vector aligned
	         #pragma ivdep
	         #pragma vector vectorlength(8)
	         #pragma vector multiple_gather_scatter_by_shuffles
	         #pragma vector always
                 for(i = m1; i != n; i += 4) {
                      t0          = (double)i;
                      psi0        = om0*t0+phi0;
                      s0          = rho0+rho1*std::sin(psi0);
                      rhot_s[i]   = s0;
                      c0          = rho0+rho1*std::cos(psi0);
                      rhot_c[i]   = c0;
                      t1          = (double)i+1;
                      psi1        = om0*t1+phi0;
                      s1          = rho0+rho1*std::sin(psi1);
                      rhot_s[i+1] = s1;
                      c0          = rho0+rho1*std::cos(psi1);
                      rhot_c[i+1] = c1;
                      t2          = (double)i+2;
                      psi2        = om0*t2+phi0;
                      s2          = rho0+rho1*std::sin(psi2);
                      rhot_s[i+2] = s2;
                      c2          = rho0+rho1*std::cos(psi2);
                      rhot_c[i+2] = c2;
                      t3          = (double)i+3;
                      psi3        = om0*t3+phi0;
                      s3          = rho0+rho1*std::sin(psi3);
                      rhot_s[i+3] = s3;
                      c3          = rho0+rho1*std::cos(psi3);
                      rhot_c[i+3] = c3;
                                       
                 }
        }



            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void ideal_modulator_unroll_2x( float * __restrict __ATTR_ALIGN__(64) rhot_s,
                                            float * __restrict __ATTR_ALIGN__(64) rhot_c,
                                            const int32_t n,
                                            const float f0,
                                            const float phi0,
                                            const float rho0,
                                            const float rho1) {

                 constexpr float twopi = 6.283185307179586476925286766559f;
                 float t0,t1;
                 float s0,s1;
                 float c0,c1;
                 float psi0,psi1;
                 float om0;
                 int32_t i,m,m1;
                 om0 = twopi*f0;
                 m   = n%2;
                 if(m!=0) {
                    for(i = 0; i != m; ++i) {
                        t0        = (float)i;
                        psi0      = om0*t0+phi0;
                        s0        = rho0+rho1*std::sin(psi0);
                        rhot_s[i] = s0
                        c0        = rho0+rho1*std::cos(psi0);
                        rhot_c[i] = c0;
                    }
                    if(n<2) return;
                 }
                 m1 = m+1;
                 __assume_aligned(rhot_s,64);
                 __assume_aligned(rhot_c,64);
                 #pragma vector aligned
	         #pragma ivdep
	         #pragma vector vectorlength(4)
	         #pragma vector multiple_gather_scatter_by_shuffles
	         #pragma vector always
                 for(i = m1; i != n; i += 2) {
                      t0          = (float)i;
                      psi0        = om0*t0+phi0;
                      s0          = rho0+rho1*std::sin(psi0);
                      rhot_s[i]   = s0;
                      c0          = rho0+rho1*std::cos(psi0);
                      rhot_c[i]   = c0;
                      t1          = (float)i+1;
                      psi1        = om0*t1+phi0;
                      s1          = rho0+rho1*std::sin(psi1);
                      rhot_s[i+1] = s1;
                      c0          = rho0+rho1*std::cos(psi1);
                      rhot_c[i+1] = c1;
                                          
                 }
        }

  

            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void ideal_modulator_unroll_2x( double * __restrict __ATTR_ALIGN__(64) rhot_s,
                                            double * __restrict __ATTR_ALIGN__(64) rhot_c,
                                            const int32_t n,
                                            const double f0,
                                            const double phi0,
                                            const double rho0,
                                            const double rho1) {

                 constexpr double twopi = 6.283185307179586476925286766559;
                 double t0,t1;
                 double s0,s1;
                 double c0,c1;
                 double psi0,psi1;
                 double om0;
                 int32_t i,m,m1;
                 om0 = twopi*f0;
                 m   = n%2;
                 if(m!=0) {
                    for(i = 0; i != m; ++i) {
                        t0        = (double)i;
                        psi0      = om0*t0+phi0;
                        s0        = rho0+rho1*std::sin(psi0);
                        rhot_s[i] = s0
                        c0        = rho0+rho1*std::cos(psi0);
                        rhot_c[i] = c0;
                    }
                    if(n<2) return;
                 }
                 m1 = m+1;
                 __assume_aligned(rhot_s,64);
                 __assume_aligned(rhot_c,64);
                 #pragma vector aligned
	         #pragma ivdep
	         #pragma vector vectorlength(8)
	         #pragma vector multiple_gather_scatter_by_shuffles
	         #pragma vector always
                 for(i = m1; i != n; i += 2) {
                      t0          = (double)i;
                      psi0        = om0*t0+phi0;
                      s0          = rho0+rho1*std::sin(psi0);
                      rhot_s[i]   = s0;
                      c0          = rho0+rho1*std::cos(psi0);
                      rhot_c[i]   = c0;
                      t1          = (double)i+1;
                      psi1        = om0*t1+phi0;
                      s1          = rho0+rho1*std::sin(psi1);
                      rhot_s[i+1] = s1;
                      c0          = rho0+rho1*std::cos(psi1);
                      rhot_c[i+1] = c1;
                                                           
                 }
        }



        // Ошибни изготовления растра —
        //!модулятора излучения
        //!Formula 2, p. 189
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void rect_pulse_flux_unroll_16x(float * __restrict __ATTR_ALIGN__(64) Phik,
                                            const float * __restrict __ATTR_ALIGN__(64) fk,
                                            const float Phi0,
                                            const int32_t n,
                                            const float Tin) {

                  constexpr float twopi = 6.283185307179586476925286766559f;
                  float fk0,fk1,fk2,fk3,fk4,fk5,fk6,fk7;
                  float fk8,fk9,fk10,fk11,fk12,fk13,fk14,fk15;
                  float sinc0,sinc1,sinc2,sinc3,sinc4,sinc5,sinc6,sinc7;
                  float sinc8,sinc9,sinc10,sinc11,sinc12,sinc13,sinc14,sinc15;
                  float arg0,arg1,arg2,arg3,arg4,arg,arg6,arg7;
                  float arg1,arg9,arg10,arg11,arg12,arg13,arg14,arg15;
                  float hTin,Phi0fk;
                  int32_t i,m,m1;
                  hTin  = 0.5f*Tin;
                  Phi0fk= Phi0*Tin;
                  m     = n%16;
                  if(m!=0) {
                     for(i = 0; i != m; ++i) {
                          fk0       = fk[i];
                          arg0      = twopi*fk0*hTin;
                          sinc0     = std::sin(arg0)/arg0;
                          Phik[i]   = Phi0fk*sinc0;
                     }
                     if(n<16) return;
                 }
                  m1 = m+1;
                  __assume_aligned(Phik,64);
                  __assume_aligned(fk,64);
                 #pragma vector aligned
	         #pragma ivdep
	         #pragma vector vectorlength(4)
	         #pragma vector multiple_gather_scatter_by_shuffles
	         #pragma vector always
                 for(i = m1; i != n; i += 16) {
                     fk0       = fk[i];
                     arg0      = twopi*fk0*hTin;
                     sinc0     = std::sin(arg0)/arg0;
                     Phik[i]   = Phi0fk*sinc0;
                     fk1       = fk[i+1];
                     arg1      = twopi*fk1*hTin;
                     sinc1     = std::sin(arg1)/arg1;
                     Phik[i+1] = Phi0fk*sinc1;
                     fk2       = fk[i+2];
                     arg2      = twopi*fk2*hTin;
                     sinc2     = std::sin(arg2)/arg2;
                     Phik[i+2] = Phi0fk*sinc2;
                     fk3       = fk[i+3];
                     arg3      = twopi*fk3*hTin;
                     sinc3     = std::sin(arg3)/arg3;
                     Phik[i+3] = Phi0fk*sinc3;
                     fk4       = fk[i+4];
                     arg4      = twopi*fk4*hTin;
                     sinc4     = std::sin(arg4)/arg4;
                     Phik[i+4] = Phi0fk*sinc4;
                     fk5       = fk[i+5];
                     arg5      = twopi*fk5*hTin;
                     sinc5     = std::sin(arg5)/arg5;
                     Phik[i+5] = Phi0fk*sinc5;
                     fk6       = fk[i+6];
                     arg6      = twopi*fk6*hTin;
                     sinc6     = std::sin(arg6)/arg6;
                     Phik(i+6) = Phi0fk*sinc6;
                     fk7       = fk[i+7];
                     arg7      = twopi*fk7*hTin;
                     sinc7     = std::sin(arg7)/arg7;
                     Phik[i+7] = Phi0fk*sinc7;
                     fk8       = fk[i+8];
                     arg8      = twopi*fk8*hTin;
                     sinc8     = std::sin(arg8)/arg8;
                     Phik[i+8] = Phi0fk*sinc8;
                     fk9       = fk[i+9];
                     arg9      = twopi*fk9*hTin;
                     sinc9     = std::sin(arg9)/arg9;
                     Phik[i+9] = Phi0fk*sinc1;
                     fk10      = fk[i+10];
                     arg10     = twopi*fk10*hTin;
                     sinc10    = std::sin(arg10)/arg10;
                     Phik[i+10]= Phi0fk*sinc10;
                     fk11      = fk[i+11];
                     arg11     = twopi*fk11*hTin;
                     sinc11    = std::sin(arg11)/arg11;
                     Phik[i+11]= Phi0fk*sinc3;
                     fk12      = fk[i+12];
                     arg12     = twopi*fk12*hTin;
                     sinc12    = std::sin(arg12)/arg12;
                     Phik[i+12]= Phi0fk*sinc12;
                     fk13      = fk[i+13];
                     arg13     = twopi*fk13*hTin;
                     sinc13    = std::sin(arg13)/arg13;
                     Phik[i+13]= Phi0fk*sinc13;
                     fk14      = fk[i+14];
                     arg14     = twopi*fk14*hTin;
                     sinc14    = std::sin(arg14)/arg14;
                     Phik[i+14]= Phi0fk*sinc14;
                     fk15      = fk[i+15];
                     arg15     = twopi*fk15*hTin;
                     sinc15    = std::sin(arg15)/arg15;
                     Phik[i+15]= Phi0fk*sinc15 ;
                 }
         }



            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void rect_pulse_flux_unroll_16x(double * __restrict __ATTR_ALIGN__(64) Phik,
                                            const double * __restrict __ATTR_ALIGN__(64) fk,
                                            const double Phi0,
                                            const int32_t n,
                                            const double Tin) {

                  constexpr double twopi = 6.283185307179586476925286766559;
                  double fk0,fk1,fk2,fk3,fk4,fk5,fk6,fk7;
                  double fk8,fk9,fk10,fk11,fk12,fk13,fk14,fk15;
                  double sinc0,sinc1,sinc2,sinc3,sinc4,sinc5,sinc6,sinc7;
                  double sinc8,sinc9,sinc10,sinc11,sinc12,sinc13,sinc14,sinc15;
                  double arg0,arg1,arg2,arg3,arg4,arg,arg6,arg7;
                  double arg1,arg9,arg10,arg11,arg12,arg13,arg14,arg15;
                  double hTin,Phi0fk;
                  int32_t i,m,m1;
                  hTin  = 0.5*Tin;
                  Phi0fk= Phi0*Tin;
                  m     = n%16;
                  if(m!=0) {
                     for(i = 0; i != m; ++i) {
                          fk0       = fk[i];
                          arg0      = twopi*fk0*hTin;
                          sinc0     = std::sin(arg0)/arg0;
                          Phik[i]   = Phi0fk*sinc0;
                     }
                     if(n<16) return;
                 }
                  m1 = m+1;
                  __assume_aligned(Phik,64);
                  __assume_aligned(fk,64);
                 #pragma vector aligned
	         #pragma ivdep
	         #pragma vector vectorlength(8)
	         #pragma vector multiple_gather_scatter_by_shuffles
	         #pragma vector always
                 for(i = m1; i != n; i += 16) {
                     fk0       = fk[i];
                     arg0      = twopi*fk0*hTin;
                     sinc0     = std::sin(arg0)/arg0;
                     Phik[i]   = Phi0fk*sinc0;
                     fk1       = fk[i+1];
                     arg1      = twopi*fk1*hTin;
                     sinc1     = std::sin(arg1)/arg1;
                     Phik[i+1] = Phi0fk*sinc1;
                     fk2       = fk[i+2];
                     arg2      = twopi*fk2*hTin;
                     sinc2     = std::sin(arg2)/arg2;
                     Phik[i+2] = Phi0fk*sinc2;
                     fk3       = fk[i+3];
                     arg3      = twopi*fk3*hTin;
                     sinc3     = std::sin(arg3)/arg3;
                     Phik[i+3] = Phi0fk*sinc3;
                     fk4       = fk[i+4];
                     arg4      = twopi*fk4*hTin;
                     sinc4     = std::sin(arg4)/arg4;
                     Phik[i+4] = Phi0fk*sinc4;
                     fk5       = fk[i+5];
                     arg5      = twopi*fk5*hTin;
                     sinc5     = std::sin(arg5)/arg5;
                     Phik[i+5] = Phi0fk*sinc5;
                     fk6       = fk[i+6];
                     arg6      = twopi*fk6*hTin;
                     sinc6     = std::sin(arg6)/arg6;
                     Phik(i+6) = Phi0fk*sinc6;
                     fk7       = fk[i+7];
                     arg7      = twopi*fk7*hTin;
                     sinc7     = std::sin(arg7)/arg7;
                     Phik[i+7] = Phi0fk*sinc7;
                     fk8       = fk[i+8];
                     arg8      = twopi*fk8*hTin;
                     sinc8     = std::sin(arg8)/arg8;
                     Phik[i+8] = Phi0fk*sinc8;
                     fk9       = fk[i+9];
                     arg9      = twopi*fk9*hTin;
                     sinc9     = std::sin(arg9)/arg9;
                     Phik[i+9] = Phi0fk*sinc1;
                     fk10      = fk[i+10];
                     arg10     = twopi*fk10*hTin;
                     sinc10    = std::sin(arg10)/arg10;
                     Phik[i+10]= Phi0fk*sinc10;
                     fk11      = fk[i+11];
                     arg11     = twopi*fk11*hTin;
                     sinc11    = std::sin(arg11)/arg11;
                     Phik[i+11]= Phi0fk*sinc3;
                     fk12      = fk[i+12];
                     arg12     = twopi*fk12*hTin;
                     sinc12    = std::sin(arg12)/arg12;
                     Phik[i+12]= Phi0fk*sinc12;
                     fk13      = fk[i+13];
                     arg13     = twopi*fk13*hTin;
                     sinc13    = std::sin(arg13)/arg13;
                     Phik[i+13]= Phi0fk*sinc13;
                     fk14      = fk[i+14];
                     arg14     = twopi*fk14*hTin;
                     sinc14    = std::sin(arg14)/arg14;
                     Phik[i+14]= Phi0fk*sinc14;
                     fk15      = fk[i+15];
                     arg15     = twopi*fk15*hTin;
                     sinc15    = std::sin(arg15)/arg15;
                     Phik[i+15]= Phi0fk*sinc15 ;
                 }
         }


            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void rect_pulse_flux_unroll_8x(float * __restrict __ATTR_ALIGN__(64) Phik,
                                            const float * __restrict __ATTR_ALIGN__(64) fk,
                                            const float Phi0,
                                            const int32_t n,
                                            const float Tin) {

                  constexpr float twopi = 6.283185307179586476925286766559f;
                  float fk0,fk1,fk2,fk3,fk4,fk5,fk6,fk7;
                  float sinc0,sinc1,sinc2,sinc3,sinc4,sinc5,sinc6,sinc7;
                  float arg0,arg1,arg2,arg3,arg4,arg,arg6,arg7;
                  float hTin,Phi0fk;
                  int32_t i,m,m1;
                  hTin  = 0.5f*Tin;
                  Phi0fk= Phi0*Tin;
                  m     = n%8;
                  if(m!=0) {
                     for(i = 0; i != m; ++i) {
                          fk0       = fk[i];
                          arg0      = twopi*fk0*hTin;
                          sinc0     = std::sin(arg0)/arg0;
                          Phik[i]   = Phi0fk*sinc0;
                     }
                     if(n<8) return;
                 }
                  m1 = m+1;
                  __assume_aligned(Phik,64);
                  __assume_aligned(fk,64);
                 #pragma vector aligned
	         #pragma ivdep
	         #pragma vector vectorlength(4)
	         #pragma vector multiple_gather_scatter_by_shuffles
	         #pragma vector always
                 for(i = m1; i != n; i += 8) {
                     fk0       = fk[i];
                     arg0      = twopi*fk0*hTin;
                     sinc0     = std::sin(arg0)/arg0;
                     Phik[i]   = Phi0fk*sinc0;
                     fk1       = fk[i+1];
                     arg1      = twopi*fk1*hTin;
                     sinc1     = std::sin(arg1)/arg1;
                     Phik[i+1] = Phi0fk*sinc1;
                     fk2       = fk[i+2];
                     arg2      = twopi*fk2*hTin;
                     sinc2     = std::sin(arg2)/arg2;
                     Phik[i+2] = Phi0fk*sinc2;
                     fk3       = fk[i+3];
                     arg3      = twopi*fk3*hTin;
                     sinc3     = std::sin(arg3)/arg3;
                     Phik[i+3] = Phi0fk*sinc3;
                     fk4       = fk[i+4];
                     arg4      = twopi*fk4*hTin;
                     sinc4     = std::sin(arg4)/arg4;
                     Phik[i+4] = Phi0fk*sinc4;
                     fk5       = fk[i+5];
                     arg5      = twopi*fk5*hTin;
                     sinc5     = std::sin(arg5)/arg5;
                     Phik[i+5] = Phi0fk*sinc5;
                     fk6       = fk[i+6];
                     arg6      = twopi*fk6*hTin;
                     sinc6     = std::sin(arg6)/arg6;
                     Phik(i+6) = Phi0fk*sinc6;
                     fk7       = fk[i+7];
                     arg7      = twopi*fk7*hTin;
                     sinc7     = std::sin(arg7)/arg7;
                     Phik[i+7] = Phi0fk*sinc7;
                    
                 }
         }


 
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void rect_pulse_flux_unroll_8x(double * __restrict __ATTR_ALIGN__(64) Phik,
                                            const double * __restrict __ATTR_ALIGN__(64) fk,
                                            const double Phi0,
                                            const int32_t n,
                                            const double Tin) {

                  constexpr double twopi = 6.283185307179586476925286766559;
                  double fk0,fk1,fk2,fk3,fk4,fk5,fk6,fk7;
                  double sinc0,sinc1,sinc2,sinc3,sinc4,sinc5,sinc6,sinc7;
                  double arg0,arg1,arg2,arg3,arg4,arg,arg6,arg7;
                  double hTin,Phi0fk;
                  int32_t i,m,m1;
                  hTin  = 0.5*Tin;
                  Phi0fk= Phi0*Tin;
                  m     = n%8;
                  if(m!=0) {
                     for(i = 0; i != m; ++i) {
                          fk0       = fk[i];
                          arg0      = twopi*fk0*hTin;
                          sinc0     = std::sin(arg0)/arg0;
                          Phik[i]   = Phi0fk*sinc0;
                     }
                     if(n<8) return;
                 }
                  m1 = m+1;
                  __assume_aligned(Phik,64);
                  __assume_aligned(fk,64);
                 #pragma vector aligned
	         #pragma ivdep
	         #pragma vector vectorlength(8)
	         #pragma vector multiple_gather_scatter_by_shuffles
	         #pragma vector always
                 for(i = m1; i != n; i += 8) {
                     fk0       = fk[i];
                     arg0      = twopi*fk0*hTin;
                     sinc0     = std::sin(arg0)/arg0;
                     Phik[i]   = Phi0fk*sinc0;
                     fk1       = fk[i+1];
                     arg1      = twopi*fk1*hTin;
                     sinc1     = std::sin(arg1)/arg1;
                     Phik[i+1] = Phi0fk*sinc1;
                     fk2       = fk[i+2];
                     arg2      = twopi*fk2*hTin;
                     sinc2     = std::sin(arg2)/arg2;
                     Phik[i+2] = Phi0fk*sinc2;
                     fk3       = fk[i+3];
                     arg3      = twopi*fk3*hTin;
                     sinc3     = std::sin(arg3)/arg3;
                     Phik[i+3] = Phi0fk*sinc3;
                     fk4       = fk[i+4];
                     arg4      = twopi*fk4*hTin;
                     sinc4     = std::sin(arg4)/arg4;
                     Phik[i+4] = Phi0fk*sinc4;
                     fk5       = fk[i+5];
                     arg5      = twopi*fk5*hTin;
                     sinc5     = std::sin(arg5)/arg5;
                     Phik[i+5] = Phi0fk*sinc5;
                     fk6       = fk[i+6];
                     arg6      = twopi*fk6*hTin;
                     sinc6     = std::sin(arg6)/arg6;
                     Phik(i+6) = Phi0fk*sinc6;
                     fk7       = fk[i+7];
                     arg7      = twopi*fk7*hTin;
                     sinc7     = std::sin(arg7)/arg7;
                     Phik[i+7] = Phi0fk*sinc7;
                    
                 }
         }




            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void rect_pulse_flux_unroll_4x(float * __restrict __ATTR_ALIGN__(64) Phik,
                                            const float * __restrict __ATTR_ALIGN__(64) fk,
                                            const float Phi0,
                                            const int32_t n,
                                            const float Tin) {

                  constexpr float twopi = 6.283185307179586476925286766559f;
                  float fk0,fk1,fk2,fk3;
                  float sinc0,sinc1,sinc2,sinc3;
                  float arg0,arg1,arg2,arg3;
                  float hTin,Phi0fk;
                  int32_t i,m,m1;
                  hTin  = 0.5f*Tin;
                  Phi0fk= Phi0*Tin;
                  m     = n%4;
                  if(m!=0) {
                     for(i = 0; i != m; ++i) {
                          fk0       = fk[i];
                          arg0      = twopi*fk0*hTin;
                          sinc0     = std::sin(arg0)/arg0;
                          Phik[i]   = Phi0fk*sinc0;
                     }
                     if(n<4) return;
                 }
                  m1 = m+1;
                  __assume_aligned(Phik,64);
                  __assume_aligned(fk,64);
                 #pragma vector aligned
	         #pragma ivdep
	         #pragma vector vectorlength(4)
	         #pragma vector multiple_gather_scatter_by_shuffles
	         #pragma vector always
                 for(i = m1; i != n; i += 4) {
                     fk0       = fk[i];
                     arg0      = twopi*fk0*hTin;
                     sinc0     = std::sin(arg0)/arg0;
                     Phik[i]   = Phi0fk*sinc0;
                     fk1       = fk[i+1];
                     arg1      = twopi*fk1*hTin;
                     sinc1     = std::sin(arg1)/arg1;
                     Phik[i+1] = Phi0fk*sinc1;
                     fk2       = fk[i+2];
                     arg2      = twopi*fk2*hTin;
                     sinc2     = std::sin(arg2)/arg2;
                     Phik[i+2] = Phi0fk*sinc2;
                     fk3       = fk[i+3];
                     arg3      = twopi*fk3*hTin;
                     sinc3     = std::sin(arg3)/arg3;
                     Phik[i+3] = Phi0fk*sinc3;
                                        
                 }
         }



            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void rect_pulse_flux_unroll_4x(double * __restrict __ATTR_ALIGN__(64) Phik,
                                            const double * __restrict __ATTR_ALIGN__(64) fk,
                                            const double Phi0,
                                            const int32_t n,
                                            const double Tin) {

                  constexpr double twopi = 6.283185307179586476925286766559;
                  double fk0,fk1,fk2,fk3;
                  double sinc0,sinc1,sinc2,sinc3;
                  double arg0,arg1,arg2,arg3;
                  double hTin,Phi0fk;
                  int32_t i,m,m1;
                  hTin  = 0.5*Tin;
                  Phi0fk= Phi0*Tin;
                  m     = n%4;
                  if(m!=0) {
                     for(i = 0; i != m; ++i) {
                          fk0       = fk[i];
                          arg0      = twopi*fk0*hTin;
                          sinc0     = std::sin(arg0)/arg0;
                          Phik[i]   = Phi0fk*sinc0;
                     }
                     if(n<4) return;
                 }
                  m1 = m+1;
                  __assume_aligned(Phik,64);
                  __assume_aligned(fk,64);
                 #pragma vector aligned
	         #pragma ivdep
	         #pragma vector vectorlength(8)
	         #pragma vector multiple_gather_scatter_by_shuffles
	         #pragma vector always
                 for(i = m1; i != n; i += 4) {
                     fk0       = fk[i];
                     arg0      = twopi*fk0*hTin;
                     sinc0     = std::sin(arg0)/arg0;
                     Phik[i]   = Phi0fk*sinc0;
                     fk1       = fk[i+1];
                     arg1      = twopi*fk1*hTin;
                     sinc1     = std::sin(arg1)/arg1;
                     Phik[i+1] = Phi0fk*sinc1;
                     fk2       = fk[i+2];
                     arg2      = twopi*fk2*hTin;
                     sinc2     = std::sin(arg2)/arg2;
                     Phik[i+2] = Phi0fk*sinc2;
                     fk3       = fk[i+3];
                     arg3      = twopi*fk3*hTin;
                     sinc3     = std::sin(arg3)/arg3;
                     Phik[i+3] = Phi0fk*sinc3;
                                       
                 }
         }



            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void rect_pulse_flux_unroll_2x(float * __restrict __ATTR_ALIGN__(64) Phik,
                                            const float * __restrict __ATTR_ALIGN__(64) fk,
                                            const float Phi0,
                                            const int32_t n,
                                            const float Tin) {

                  constexpr float twopi = 6.283185307179586476925286766559f;
                  float fk0,fk1,fk2,fk3;
                  float sinc0,sinc1;
                  float arg0,arg1;
                  float hTin,Phi0fk;
                  int32_t i,m,m1;
                  hTin  = 0.5f*Tin;
                  Phi0fk= Phi0*Tin;
                  m     = n%2;
                  if(m!=0) {
                     for(i = 0; i != m; ++i) {
                          fk0       = fk[i];
                          arg0      = twopi*fk0*hTin;
                          sinc0     = std::sin(arg0)/arg0;
                          Phik[i]   = Phi0fk*sinc0;
                     }
                     if(n<2) return;
                 }
                  m1 = m+1;
                  __assume_aligned(Phik,64);
                  __assume_aligned(fk,64);
                 #pragma vector aligned
	         #pragma ivdep
	         #pragma vector vectorlength(4)
	         #pragma vector multiple_gather_scatter_by_shuffles
	         #pragma vector always
                 for(i = m1; i != n; i += 2) {
                     fk0       = fk[i];
                     arg0      = twopi*fk0*hTin;
                     sinc0     = std::sin(arg0)/arg0;
                     Phik[i]   = Phi0fk*sinc0;
                     fk1       = fk[i+1];
                     arg1      = twopi*fk1*hTin;
                     sinc1     = std::sin(arg1)/arg1;
                     Phik[i+1] = Phi0fk*sinc1;
                                                            
                 }
         }



            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void rect_pulse_flux_unroll_2x(double * __restrict __ATTR_ALIGN__(64) Phik,
                                            const double * __restrict __ATTR_ALIGN__(64) fk,
                                            const double Phi0,
                                            const int32_t n,
                                            const double Tin) {

                  constexpr double twopi = 6.283185307179586476925286766559;
                  double fk0,fk1;
                  double sinc0,sinc1;
                  double arg0,arg1;
                  double hTin,Phi0fk;
                  int32_t i,m,m1;
                  hTin  = 0.5*Tin;
                  Phi0fk= Phi0*Tin;
                  m     = n%2;
                  if(m!=0) {
                     for(i = 0; i != m; ++i) {
                          fk0       = fk[i];
                          arg0      = twopi*fk0*hTin;
                          sinc0     = std::sin(arg0)/arg0;
                          Phik[i]   = Phi0fk*sinc0;
                     }
                     if(n<2) return;
                 }
                  m1 = m+1;
                  __assume_aligned(Phik,64);
                  __assume_aligned(fk,64);
                 #pragma vector aligned
	         #pragma ivdep
	         #pragma vector vectorlength(8)
	         #pragma vector multiple_gather_scatter_by_shuffles
	         #pragma vector always
                 for(i = m1; i != n; i += 2) {
                     fk0       = fk[i];
                     arg0      = twopi*fk0*hTin;
                     sinc0     = std::sin(arg0)/arg0;
                     Phik[i]   = Phi0fk*sinc0;
                     fk1       = fk[i+1];
                     arg1      = twopi*fk1*hTin;
                     sinc1     = std::sin(arg1)/arg1;
                     Phik[i+1] = Phi0fk*sinc1;
                   
                 }
         }


         
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void rect_pulse_amp_unroll_16x(float * __restrict __ATTR_ALIGN__(64) Ak,
                                           const float * __restrict __ATTR_ALIGN__(64) Phik,
                                           const float Phi0,
                                           const int32_t n,
                                           const float T,
                                           const float * __restrict __ATTR_ALIGN__(64) k,
                                           const float tin) {

                 constexpr float pi2 = 1.5707963267948966192313216916398f;
                 float phik0,phik1,phik2,phik3,phik4,phik5,phik6,phik7;
                 float phik8,phik9,phik10,phik11,phik12,phik13,phik14,phik15;
                 float arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7;
                 float arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15;
                 float k0,k1,k2,k3,k4,k5,k6,k7;
                 float k8,k9,k10,k11,k12,k13,k14,k15;
                 float twoT,kpi2;
                 int32_t i,m,m1;
                 twoT = 2.0f/T;
                 if(approximatelyEqual(tin,twoT,
                             std::numeric_limits<float>::epsilon())) {
                      m = n%16;
                      if(m!=0) {
                         for(i = 0; i != m; ++i) {
                              k0    = k[i];
                              arg0  = k0*pi2;
                              Ak[i] = Phi0*(std::sin(arg0)/arg0); 
                          }
                          if(n<16) return;
                      }
                      m1 = m+1;
                      __assume_aligned(Ak,64);
                      __assume_aligned(k,64);
                      #pragma vector aligned
	              #pragma ivdep
	              #pragma vector vectorlength(4)
	              #pragma vector multiple_gather_scatter_by_shuffles
	              #pragma vector always
                      for(i = m1; i != n; i += 16) {
                          k0      = k[i];
                          arg0    = k0*pi2;
                          Ak[i] = Phi0*(std::sin(arg0)/arg0);
                          k1      = k(i+1);
                          arg1    = k1*pi2;
                          Ak[i+1] = Phi0*(std::sin(arg1)/arg1); 
                          k2      = k[i+2];
                          arg2    = k2*pi2;
                          Ak[i+2] = Phi0*(std::sin(arg2)/arg2); 
                          k3      = k[i+3];
                          arg3    = k3*pi2;
                          Ak[i+3] = Phi0*(std::sin(arg3)/arg3); 
                          k4      = k[i+4];
                          arg4    = k4*pi2;
                          Ak[i+4] = Phi0*(std::sin(arg4)/arg4); 
                          k5      = k[i+5];
                          arg5    = k5*pi2;
                          Ak[i+5] = Phi0*(std::sin(arg5)/arg5); 
                          k6      = k[i+6];
                          arg6    = k6*pi2;
                          Ak[i+6] = Phi0*(std::sin(arg6)/arg6);
                          k7      = k[i+7];
                          arg7    = k7*pi2;
                          Ak[i+7] = Phi0*(std::sin(arg7)/arg7);   
                          k8      = k[i+8];
                          arg8    = k8*pi2;
                          Ak[i+8] = Phi0*(std::sin(arg8)/arg8); 
                          k9      = k(i+9);
                          arg9    = k9*pi2;
                          Ak[i+9] = Phi0*(std::sin(arg9)/arg9); 
                          k10     = k[i+10];
                          arg10   = k10*pi2;
                          Ak[i+10]= Phi0*(std::sin(arg10)/arg10); 
                          k11     = k[i+11];
                          arg11   = k11*pi2;
                          Ak[i+11]= Phi0*(std::sin(arg11)/arg11); 
                          k12     = k[i+12];
                          arg12   = k12*pi2;
                          Ak[i+12]= Phi0*(std::sin(arg12)/arg12); 
                          k13     = k[i+13];
                          arg13   = k13*pi2;
                          Ak[i+13]= Phi0*(std::sin(arg13)/arg13); 
                          k14     = k(i+14);
                          arg14   = k14*pi2;
                          Ak[i+14]= Phi0*(std::sin(arg14)/arg14);
                          k15     = k[i+15];
                          arg15   = k15*pi2;
                          Ak[i+15]= Phi0*(std::sin(arg15)/arg15);  
                      }
                      
                     }
                     else {
                          m = n%16;
                          if(m!=0) {
                             for(i = 0; i != m; ++i) {
                                 phik0 = Phik[i];
                                 Ak[i] = twoT*phik0;
                             }
                             if(n<16) return;
                          }
                          m1 = m+1;
                          __assume_aligned(Ak,64);
                          __assume_aligned(k,64);
                          #pragma vector aligned
	                  #pragma ivdep
	                  #pragma vector vectorlength(4)
	                  #pragma vector multiple_gather_scatter_by_shuffles
	                  #pragma vector always
                          for(i = m1; i != n; i += 16) {
                              phik0   = Phik[i];
                              Ak[i]   = twoT*phik0;
                              phik1   = Phik[i+1];
                              Ak[i+1] = twoT*phik1;
                              phik2   = Phik[i+2];
                              Ak[i+2] = twoT*phik2;
                              phik3   = Phik[i+3];
                              Ak[i+3] = twoT*phik3;
                              phik4   = Phik[i+4];
                              Ak[i+4] = twoT*phik4;
                              phik5   = Phik[i+5];
                              Ak[i+5] = twoT*phik5;
                              phik6   = Phik[i+6];
                              Ak[i+6] = twoT*phik6;
                              phik7   = Phik[i+7];
                              Ak[i+7] = twoT*phik7;
                              phik8   = Phik[i+8];
                              Ak[i+8] = twoT*phik8;
                              phik9   = Phik[i+9];
                              Ak[i+9] = twoT*phik9;
                              phik10  = Phik[i+10];
                              Ak[i+10]= twoT*phik10;
                              phik11  = Phik[i+11];
                              Ak[i+11]= twoT*phik11;
                              phik12  = Phik[i+12];
                              Ak[i+12]= twoT*phik12;
                              phik13  = Phik[i+13];
                              Ak[i+13]= twoT*phik13;
                              phik14  = Phik[i+14];
                              Ak[i+14]= twoT*phik14;
                              phik15  = Phik[i+15];
                              Ak[i+15] = twoT*phik15;
                          }
                 }
          }


          
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void rect_pulse_amp_unroll_16x(double * __restrict __ATTR_ALIGN__(64) Ak,
                                           const double * __restrict __ATTR_ALIGN__(64) Phik,
                                           const double Phi0,
                                           const int32_t n,
                                           const double T,
                                           const double * __restrict __ATTR_ALIGN__(64) k,
                                           const double tin) {

                 constexpr double pi2 = 1.5707963267948966192313216916398;
                 double phik0,phik1,phik2,phik3,phik4,phik5,phik6,phik7;
                 double phik8,phik9,phik10,phik11,phik12,phik13,phik14,phik15;
                 double arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7;
                 double arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15;
                 double k0,k1,k2,k3,k4,k5,k6,k7;
                 double k8,k9,k10,k11,k12,k13,k14,k15;
                 double twoT,kpi2;
                 int32_t i,m,m1;
                 twoT = 2.0/T;
                 if(approximatelyEqual(tin,twoT,
                             std::numeric_limits<float>::epsilon())) {
                      m = n%16;
                      if(m!=0) {
                         for(i = 0; i != m; ++i) {
                              k0    = k[i];
                              arg0  = k0*pi2;
                              Ak[i] = Phi0*(std::sin(arg0)/arg0); 
                          }
                          if(n<16) return;
                      }
                      m1 = m+1;
                      __assume_aligned(Ak,64);
                      __assume_aligned(k,64);
                      #pragma vector aligned
	              #pragma ivdep
	              #pragma vector vectorlength(8)
	              #pragma vector multiple_gather_scatter_by_shuffles
	              #pragma vector always
                      for(i = m1; i != n; i += 16) {
                          k0      = k[i];
                          arg0    = k0*pi2;
                          Ak[i] = Phi0*(std::sin(arg0)/arg0);
                          k1      = k(i+1);
                          arg1    = k1*pi2;
                          Ak[i+1] = Phi0*(std::sin(arg1)/arg1); 
                          k2      = k[i+2];
                          arg2    = k2*pi2;
                          Ak[i+2] = Phi0*(std::sin(arg2)/arg2); 
                          k3      = k[i+3];
                          arg3    = k3*pi2;
                          Ak[i+3] = Phi0*(std::sin(arg3)/arg3); 
                          k4      = k[i+4];
                          arg4    = k4*pi2;
                          Ak[i+4] = Phi0*(std::sin(arg4)/arg4); 
                          k5      = k[i+5];
                          arg5    = k5*pi2;
                          Ak[i+5] = Phi0*(std::sin(arg5)/arg5); 
                          k6      = k[i+6];
                          arg6    = k6*pi2;
                          Ak[i+6] = Phi0*(std::sin(arg6)/arg6);
                          k7      = k[i+7];
                          arg7    = k7*pi2;
                          Ak[i+7] = Phi0*(std::sin(arg7)/arg7);   
                          k8      = k[i+8];
                          arg8    = k8*pi2;
                          Ak[i+8] = Phi0*(std::sin(arg8)/arg8); 
                          k9      = k(i+9);
                          arg9    = k9*pi2;
                          Ak[i+9] = Phi0*(std::sin(arg9)/arg9); 
                          k10     = k[i+10];
                          arg10   = k10*pi2;
                          Ak[i+10]= Phi0*(std::sin(arg10)/arg10); 
                          k11     = k[i+11];
                          arg11   = k11*pi2;
                          Ak[i+11]= Phi0*(std::sin(arg11)/arg11); 
                          k12     = k[i+12];
                          arg12   = k12*pi2;
                          Ak[i+12]= Phi0*(std::sin(arg12)/arg12); 
                          k13     = k[i+13];
                          arg13   = k13*pi2;
                          Ak[i+13]= Phi0*(std::sin(arg13)/arg13); 
                          k14     = k(i+14);
                          arg14   = k14*pi2;
                          Ak[i+14]= Phi0*(std::sin(arg14)/arg14);
                          k15     = k[i+15];
                          arg15   = k15*pi2;
                          Ak[i+15]= Phi0*(std::sin(arg15)/arg15);  
                      }
                      
                     }
                     else {
                          m = n%16;
                          if(m!=0) {
                             for(i = 0; i != m; ++i) {
                                 phik0 = Phik[i];
                                 Ak[i] = twoT*phik0;
                             }
                             if(n<16) return;
                          }
                          m1 = m+1;
                          __assume_aligned(Ak,64);
                          __assume_aligned(k,64);
                          #pragma vector aligned
	                  #pragma ivdep
	                  #pragma vector vectorlength(8)
	                  #pragma vector multiple_gather_scatter_by_shuffles
	                  #pragma vector always
                          for(i = m1; i != n; i += 16) {
                              phik0   = Phik[i];
                              Ak[i]   = twoT*phik0;
                              phik1   = Phik[i+1];
                              Ak[i+1] = twoT*phik1;
                              phik2   = Phik[i+2];
                              Ak[i+2] = twoT*phik2;
                              phik3   = Phik[i+3];
                              Ak[i+3] = twoT*phik3;
                              phik4   = Phik[i+4];
                              Ak[i+4] = twoT*phik4;
                              phik5   = Phik[i+5];
                              Ak[i+5] = twoT*phik5;
                              phik6   = Phik[i+6];
                              Ak[i+6] = twoT*phik6;
                              phik7   = Phik[i+7];
                              Ak[i+7] = twoT*phik7;
                              phik8   = Phik[i+8];
                              Ak[i+8] = twoT*phik8;
                              phik9   = Phik[i+9];
                              Ak[i+9] = twoT*phik9;
                              phik10  = Phik[i+10];
                              Ak[i+10]= twoT*phik10;
                              phik11  = Phik[i+11];
                              Ak[i+11]= twoT*phik11;
                              phik12  = Phik[i+12];
                              Ak[i+12]= twoT*phik12;
                              phik13  = Phik[i+13];
                              Ak[i+13]= twoT*phik13;
                              phik14  = Phik[i+14];
                              Ak[i+14]= twoT*phik14;
                              phik15  = Phik[i+15];
                              Ak[i+15] = twoT*phik15;
                          }
                 }
          }



            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void rect_pulse_amp_unroll_8x(float * __restrict __ATTR_ALIGN__(64) Ak,
                                           const float * __restrict __ATTR_ALIGN__(64) Phik,
                                           const float Phi0,
                                           const int32_t n,
                                           const float T,
                                           const float * __restrict __ATTR_ALIGN__(64) k,
                                           const float tin) {

                 constexpr float pi2 = 1.5707963267948966192313216916398f;
                 float phik0,phik1,phik2,phik3,phik4,phik5,phik6,phik7;
                 float arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7;
                 float k0,k1,k2,k3,k4,k5,k6,k7;
                 float twoT,kpi2;
                 int32_t i,m,m1;
                 twoT = 2.0f/T;
                 if(approximatelyEqual(tin,twoT,
                             std::numeric_limits<float>::epsilon())) {
                      m = n%8;
                      if(m!=0) {
                         for(i = 0; i != m; ++i) {
                              k0    = k[i];
                              arg0  = k0*pi2;
                              Ak[i] = Phi0*(std::sin(arg0)/arg0); 
                          }
                          if(n<8) return;
                      }
                      m1 = m+1;
                      __assume_aligned(Ak,64);
                      __assume_aligned(k,64);
                      #pragma vector aligned
	              #pragma ivdep
	              #pragma vector vectorlength(4)
	              #pragma vector multiple_gather_scatter_by_shuffles
	              #pragma vector always
                      for(i = m1; i != n; i += 8) {
                          k0      = k[i];
                          arg0    = k0*pi2;
                          Ak[i]   = Phi0*(std::sin(arg0)/arg0);
                          k1      = k(i+1);
                          arg1    = k1*pi2;
                          Ak[i+1] = Phi0*(std::sin(arg1)/arg1); 
                          k2      = k[i+2];
                          arg2    = k2*pi2;
                          Ak[i+2] = Phi0*(std::sin(arg2)/arg2); 
                          k3      = k[i+3];
                          arg3    = k3*pi2;
                          Ak[i+3] = Phi0*(std::sin(arg3)/arg3); 
                          k4      = k[i+4];
                          arg4    = k4*pi2;
                          Ak[i+4] = Phi0*(std::sin(arg4)/arg4); 
                          k5      = k[i+5];
                          arg5    = k5*pi2;
                          Ak[i+5] = Phi0*(std::sin(arg5)/arg5); 
                          k6      = k[i+6];
                          arg6    = k6*pi2;
                          Ak[i+6] = Phi0*(std::sin(arg6)/arg6);
                          k7      = k[i+7];
                          arg7    = k7*pi2;
                          Ak[i+7] = Phi0*(std::sin(arg7)/arg7);   
                        
                      }
                      
                     }
                     else {
                          m = n%8;
                          if(m!=0) {
                             for(i = 0; i != m; ++i) {
                                 phik0 = Phik[i];
                                 Ak[i] = twoT*phik0;
                             }
                             if(n<8) return;
                          }
                          m1 = m+1;
                          __assume_aligned(Ak,64);
                          __assume_aligned(k,64);
                          #pragma vector aligned
	                  #pragma ivdep
	                  #pragma vector vectorlength(4)
	                  #pragma vector multiple_gather_scatter_by_shuffles
	                  #pragma vector always
                          for(i = m1; i != n; i += 8) {
                              phik0   = Phik[i];
                              Ak[i]   = twoT*phik0;
                              phik1   = Phik[i+1];
                              Ak[i+1] = twoT*phik1;
                              phik2   = Phik[i+2];
                              Ak[i+2] = twoT*phik2;
                              phik3   = Phik[i+3];
                              Ak[i+3] = twoT*phik3;
                              phik4   = Phik[i+4];
                              Ak[i+4] = twoT*phik4;
                              phik5   = Phik[i+5];
                              Ak[i+5] = twoT*phik5;
                              phik6   = Phik[i+6];
                              Ak[i+6] = twoT*phik6;
                              phik7   = Phik[i+7];
                              Ak[i+7] = twoT*phik7;
                              
                          }
                 }
          }



            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void rect_pulse_amp_unroll_8x(double * __restrict __ATTR_ALIGN__(64) Ak,
                                           const double * __restrict __ATTR_ALIGN__(64) Phik,
                                           const double Phi0,
                                           const int32_t n,
                                           const double T,
                                           const double * __restrict __ATTR_ALIGN__(64) k,
                                           const double tin) {

                 constexpr double pi2 = 1.5707963267948966192313216916398;
                 double phik0,phik1,phik2,phik3,phik4,phik5,phik6,phik7;
                 double arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7;
                 double k0,k1,k2,k3,k4,k5,k6,k7;
                 double twoT,kpi2;
                 int32_t i,m,m1;
                 twoT = 2.0/T;
                 if(approximatelyEqual(tin,twoT,
                             std::numeric_limits<float>::epsilon())) {
                      m = n%8;
                      if(m!=0) {
                         for(i = 0; i != m; ++i) {
                              k0    = k[i];
                              arg0  = k0*pi2;
                              Ak[i] = Phi0*(std::sin(arg0)/arg0); 
                          }
                          if(n<8) return;
                      }
                      m1 = m+1;
                      __assume_aligned(Ak,64);
                      __assume_aligned(k,64);
                      #pragma vector aligned
	              #pragma ivdep
	              #pragma vector vectorlength(8)
	              #pragma vector multiple_gather_scatter_by_shuffles
	              #pragma vector always
                      for(i = m1; i != n; i += 8) {
                          k0      = k[i];
                          arg0    = k0*pi2;
                          Ak[i]   = Phi0*(std::sin(arg0)/arg0);
                          k1      = k(i+1);
                          arg1    = k1*pi2;
                          Ak[i+1] = Phi0*(std::sin(arg1)/arg1); 
                          k2      = k[i+2];
                          arg2    = k2*pi2;
                          Ak[i+2] = Phi0*(std::sin(arg2)/arg2); 
                          k3      = k[i+3];
                          arg3    = k3*pi2;
                          Ak[i+3] = Phi0*(std::sin(arg3)/arg3); 
                          k4      = k[i+4];
                          arg4    = k4*pi2;
                          Ak[i+4] = Phi0*(std::sin(arg4)/arg4); 
                          k5      = k[i+5];
                          arg5    = k5*pi2;
                          Ak[i+5] = Phi0*(std::sin(arg5)/arg5); 
                          k6      = k[i+6];
                          arg6    = k6*pi2;
                          Ak[i+6] = Phi0*(std::sin(arg6)/arg6);
                          k7      = k[i+7];
                          arg7    = k7*pi2;
                          Ak[i+7] = Phi0*(std::sin(arg7)/arg7);   
                         
                      }
                      
                     }
                     else {
                          m = n%8;
                          if(m!=0) {
                             for(i = 0; i != m; ++i) {
                                 phik0 = Phik[i];
                                 Ak[i] = twoT*phik0;
                             }
                             if(n<8) return;
                          }
                          m1 = m+1;
                          __assume_aligned(Ak,64);
                          __assume_aligned(k,64);
                          #pragma vector aligned
	                  #pragma ivdep
	                  #pragma vector vectorlength(8)
	                  #pragma vector multiple_gather_scatter_by_shuffles
	                  #pragma vector always
                          for(i = m1; i != n; i += 8) {
                              phik0   = Phik[i];
                              Ak[i]   = twoT*phik0;
                              phik1   = Phik[i+1];
                              Ak[i+1] = twoT*phik1;
                              phik2   = Phik[i+2];
                              Ak[i+2] = twoT*phik2;
                              phik3   = Phik[i+3];
                              Ak[i+3] = twoT*phik3;
                              phik4   = Phik[i+4];
                              Ak[i+4] = twoT*phik4;
                              phik5   = Phik[i+5];
                              Ak[i+5] = twoT*phik5;
                              phik6   = Phik[i+6];
                              Ak[i+6] = twoT*phik6;
                              phik7   = Phik[i+7];
                              Ak[i+7] = twoT*phik7;
                             
                          }
                 }
          }



            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void rect_pulse_amp_unroll_4x(float * __restrict __ATTR_ALIGN__(64) Ak,
                                           const float * __restrict __ATTR_ALIGN__(64) Phik,
                                           const float Phi0,
                                           const int32_t n,
                                           const float T,
                                           const float * __restrict __ATTR_ALIGN__(64) k,
                                           const float tin) {

                 constexpr float pi2 = 1.5707963267948966192313216916398f;
                 float phik0,phik1,phik2,phik3;
                 float arg0,arg1,arg2,arg3;
                 float k0,k1,k2,k3;
                 float twoT,kpi2;
                 int32_t i,m,m1;
                 twoT = 2.0f/T;
                 if(approximatelyEqual(tin,twoT,
                             std::numeric_limits<float>::epsilon())) {
                      m = n%4;
                      if(m!=0) {
                         for(i = 0; i != m; ++i) {
                              k0    = k[i];
                              arg0  = k0*pi2;
                              Ak[i] = Phi0*(std::sin(arg0)/arg0); 
                          }
                          if(n<4) return;
                      }
                      m1 = m+1;
                      __assume_aligned(Ak,64);
                      __assume_aligned(k,64);
                      #pragma vector aligned
	              #pragma ivdep
	              #pragma vector vectorlength(4)
	              #pragma vector multiple_gather_scatter_by_shuffles
	              #pragma vector always
                      for(i = m1; i != n; i += 4) {
                          k0      = k[i];
                          arg0    = k0*pi2;
                          Ak[i]   = Phi0*(std::sin(arg0)/arg0);
                          k1      = k(i+1);
                          arg1    = k1*pi2;
                          Ak[i+1] = Phi0*(std::sin(arg1)/arg1); 
                          k2      = k[i+2];
                          arg2    = k2*pi2;
                          Ak[i+2] = Phi0*(std::sin(arg2)/arg2); 
                          k3      = k[i+3];
                          arg3    = k3*pi2;
                          Ak[i+3] = Phi0*(std::sin(arg3)/arg3); 
                                               
                      }
                      
                     }
                     else {
                          m = n%4;
                          if(m!=0) {
                             for(i = 0; i != m; ++i) {
                                 phik0 = Phik[i];
                                 Ak[i] = twoT*phik0;
                             }
                             if(n<4) return;
                          }
                          m1 = m+1;
                          __assume_aligned(Ak,64);
                          __assume_aligned(k,64);
                          #pragma vector aligned
	                  #pragma ivdep
	                  #pragma vector vectorlength(4)
	                  #pragma vector multiple_gather_scatter_by_shuffles
	                  #pragma vector always
                          for(i = m1; i != n; i += 4) {
                              phik0   = Phik[i];
                              Ak[i]   = twoT*phik0;
                              phik1   = Phik[i+1];
                              Ak[i+1] = twoT*phik1;
                              phik2   = Phik[i+2];
                              Ak[i+2] = twoT*phik2;
                              phik3   = Phik[i+3];
                              Ak[i+3] = twoT*phik3;
                                                           
                          }
                 }
          }


           
            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void rect_pulse_amp_unroll_4x(double * __restrict __ATTR_ALIGN__(64) Ak,
                                           const double * __restrict __ATTR_ALIGN__(64) Phik,
                                           const double Phi0,
                                           const int32_t n,
                                           const double T,
                                           const double * __restrict __ATTR_ALIGN__(64) k,
                                           const double tin) {

                 constexpr double pi2 = 1.5707963267948966192313216916398;
                 double phik0,phik1,phik2,phik3;
                 double arg0,arg1,arg2,arg3;
                 double k0,k1,k2,k3;
                 double twoT,kpi2;
                 int32_t i,m,m1;
                 twoT = 2.0/T;
                 if(approximatelyEqual(tin,twoT,
                             std::numeric_limits<float>::epsilon())) {
                      m = n%4;
                      if(m!=0) {
                         for(i = 0; i != m; ++i) {
                              k0    = k[i];
                              arg0  = k0*pi2;
                              Ak[i] = Phi0*(std::sin(arg0)/arg0); 
                          }
                          if(n<4) return;
                      }
                      m1 = m+1;
                      __assume_aligned(Ak,64);
                      __assume_aligned(k,64);
                      #pragma vector aligned
	              #pragma ivdep
	              #pragma vector vectorlength(8)
	              #pragma vector multiple_gather_scatter_by_shuffles
	              #pragma vector always
                      for(i = m1; i != n; i += 4) {
                          k0      = k[i];
                          arg0    = k0*pi2;
                          Ak[i]   = Phi0*(std::sin(arg0)/arg0);
                          k1      = k(i+1);
                          arg1    = k1*pi2;
                          Ak[i+1] = Phi0*(std::sin(arg1)/arg1); 
                          k2      = k[i+2];
                          arg2    = k2*pi2;
                          Ak[i+2] = Phi0*(std::sin(arg2)/arg2); 
                          k3      = k[i+3];
                          arg3    = k3*pi2;
                          Ak[i+3] = Phi0*(std::sin(arg3)/arg3); 
                                                  
                      }
                      
                     }
                     else {
                          m = n%4;
                          if(m!=0) {
                             for(i = 0; i != m; ++i) {
                                 phik0 = Phik[i];
                                 Ak[i] = twoT*phik0;
                             }
                             if(n<4) return;
                          }
                          m1 = m+1;
                          __assume_aligned(Ak,64);
                          __assume_aligned(k,64);
                          #pragma vector aligned
	                  #pragma ivdep
	                  #pragma vector vectorlength(8)
	                  #pragma vector multiple_gather_scatter_by_shuffles
	                  #pragma vector always
                          for(i = m1; i != n; i += 4) {
                              phik0   = Phik[i];
                              Ak[i]   = twoT*phik0;
                              phik1   = Phik[i+1];
                              Ak[i+1] = twoT*phik1;
                              phik2   = Phik[i+2];
                              Ak[i+2] = twoT*phik2;
                              phik3   = Phik[i+3];
                              Ak[i+3] = twoT*phik3;
                                                          
                          }
                 }
          }



            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void rect_pulse_amp_unroll_2x(float * __restrict __ATTR_ALIGN__(64) Ak,
                                           const float * __restrict __ATTR_ALIGN__(64) Phik,
                                           const float Phi0,
                                           const int32_t n,
                                           const float T,
                                           const float * __restrict __ATTR_ALIGN__(64) k,
                                           const float tin) {

                 constexpr float pi2 = 1.5707963267948966192313216916398f;
                 float phik0,phik1;
                 float arg0,arg1;
                 float k0,k1;
                 float twoT,kpi2;
                 int32_t i,m,m1;
                 twoT = 2.0f/T;
                 if(approximatelyEqual(tin,twoT,
                             std::numeric_limits<float>::epsilon())) {
                      m = n%2;
                      if(m!=0) {
                         for(i = 0; i != m; ++i) {
                              k0    = k[i];
                              arg0  = k0*pi2;
                              Ak[i] = Phi0*(std::sin(arg0)/arg0); 
                          }
                          if(n<2) return;
                      }
                      m1 = m+1;
                      __assume_aligned(Ak,64);
                      __assume_aligned(k,64);
                      #pragma vector aligned
	              #pragma ivdep
	              #pragma vector vectorlength(4)
	              #pragma vector multiple_gather_scatter_by_shuffles
	              #pragma vector always
                      for(i = m1; i != n; i += 2) {
                          k0      = k[i];
                          arg0    = k0*pi2;
                          Ak[i]   = Phi0*(std::sin(arg0)/arg0);
                          k1      = k(i+1);
                          arg1    = k1*pi2;
                          Ak[i+1] = Phi0*(std::sin(arg1)/arg1); 
                                                                        
                      }
                      
                     }
                     else {
                          m = n%2;
                          if(m!=0) {
                             for(i = 0; i != m; ++i) {
                                 phik0 = Phik[i];
                                 Ak[i] = twoT*phik0;
                             }
                             if(n<2) return;
                          }
                          m1 = m+1;
                          __assume_aligned(Ak,64);
                          __assume_aligned(k,64);
                          #pragma vector aligned
	                  #pragma ivdep
	                  #pragma vector vectorlength(4)
	                  #pragma vector multiple_gather_scatter_by_shuffles
	                  #pragma vector always
                          for(i = m1; i != n; i += 2) {
                              phik0   = Phik[i];
                              Ak[i]   = twoT*phik0;
                              phik1   = Phik[i+1];
                              Ak[i+1] = twoT*phik1;
                             
                                                           
                          }
                 }
          }



            __ATTR_ALWAYS_INLINE__
	   
	   
	    static inline
	    void rect_pulse_amp_unroll_2x(double * __restrict __ATTR_ALIGN__(64) Ak,
                                           const double * __restrict __ATTR_ALIGN__(64) Phik,
                                           const double Phi0,
                                           const int32_t n,
                                           const double T,
                                           const double * __restrict __ATTR_ALIGN__(64) k,
                                           const double tin) {

                 constexpr double pi2 = 1.5707963267948966192313216916398;
                 double phik0,phik1;
                 double arg0,arg1;
                 double k0,k1;
                 double twoT,kpi2;
                 int32_t i,m,m1;
                 twoT = 2.0/T;
                 if(approximatelyEqual(tin,twoT,
                             std::numeric_limits<float>::epsilon())) {
                      m = n%2;
                      if(m!=0) {
                         for(i = 0; i != m; ++i) {
                              k0    = k[i];
                              arg0  = k0*pi2;
                              Ak[i] = Phi0*(std::sin(arg0)/arg0); 
                          }
                          if(n<2) return;
                      }
                      m1 = m+1;
                      __assume_aligned(Ak,64);
                      __assume_aligned(k,64);
                      #pragma vector aligned
	              #pragma ivdep
	              #pragma vector vectorlength(8)
	              #pragma vector multiple_gather_scatter_by_shuffles
	              #pragma vector always
                      for(i = m1; i != n; i += 2) {
                          k0      = k[i];
                          arg0    = k0*pi2;
                          Ak[i]   = Phi0*(std::sin(arg0)/arg0);
                          k1      = k(i+1);
                          arg1    = k1*pi2;
                          Ak[i+1] = Phi0*(std::sin(arg1)/arg1); 
                        
                                                  
                      }
                      
                     }
                     else {
                          m = n%2;
                          if(m!=0) {
                             for(i = 0; i != m; ++i) {
                                 phik0 = Phik[i];
                                 Ak[i] = twoT*phik0;
                             }
                             if(n<2) return;
                          }
                          m1 = m+1;
                          __assume_aligned(Ak,64);
                          __assume_aligned(k,64);
                          #pragma vector aligned
	                  #pragma ivdep
	                  #pragma vector vectorlength(8)
	                  #pragma vector multiple_gather_scatter_by_shuffles
	                  #pragma vector always
                          for(i = m1; i != n; i += 2) {
                              phik0   = Phik[i];
                              Ak[i]   = twoT*phik0;
                              phik1   = Phik[i+1];
                              Ak[i+1] = twoT*phik1;
                            
                                                          
                          }
                 }
          }
















         








         






         






   
         

           


           










           






            






         










	  






  
	  







	  











    } //eos

}// gms


#endif /*__GMS_EOS_IR_SENSOR_HPP__*/
