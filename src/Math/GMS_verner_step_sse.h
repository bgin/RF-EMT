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

#ifndef __GMS_VERNER_STEP_SSE_H__
#define __GMS_VERNER_STEP_SSE_H__ 050620221002


/*
   Adapted from: http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_verner.html
   Manually vectorized by Bernard Gingold, beniekg@gmail.com
*/

namespace file_version {

    const unsigned int GMS_VERNER_STEP_SSE_MAJOR = 1U;
    const unsigned int GMS_VERNER_STEP_SSE_MINOR = 0U;
    const unsigned int GMS_VERNER_STEP_SSE_MICRO = 0U;
    const unsigned int GMS_VERNER_STEP_SSE_FULLVER =
      1000U*GMS_VERNER_STEP_SSE_MAJOR+
      100U*GMS_VERNER_STEP_SSE_MINOR+
      10U*GMS_VERNER_STEP_SSE_MICRO;
    const char * const GMS_VERNER_STEP_SSE_CREATION_DATE = "05-06-2022 10:002 AM +00200 (SUN 05 JUN 2022 GMT+2)";
    const char * const GMS_VERNER_STEP_SSE_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_VERNER_STEP_SSE_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_VERNER_STEP_SSE_DESCRIPTION   = "Vectorized (SSE) Runge-Kutta-Verner order 8 step."

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace gms {

       namespace math {
////////////////////////////////////////////////////////////////////////////////
//  const  __m128d Runge_Kutta_Verner( const  __m128d (*f)(const  __m128d, const  __m128d), const  __m128d y0,        //
//                               const  __m128d x0, const  __m128d h, int number_of_steps );  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the Runge-Kutta_Verner method described above to     //
//     approximate the solution at x = x0 + h * number_of_steps of the initial//
//     value problem y'=f(x,y), y(x0) = y0.                                   //
//                                                                            //
//  Arguments:                                                                //
//     const  __m128d *f                                                              //
//            Pointer to the function which returns the slope at (x,y) of the //
//            integral curve of the differential equation y' = f(x,y) which   //
//            passes through the point (x0,y0).                               //
//     const  __m128d y0                                                              //
//            The initial value of y at x = x0.                               //
//     const  __m128d x0                                                              //
//            The initial value of x.                                         //
//     const  __m128d h                                                               //
//            The step size.                                                  //
//     int    number_of_steps                                                 //
//            The number of steps. Must be a nonnegative integer.             //
//                                                                            //
//  Return Values:                                                            //
//     The solution of the initial value problem y' = f(x,y), y(x0) = y0 at   //
//     x = x0 + number_of_steps * h.                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


                    
			__m128d
			verner_step_xmm2r8(__m128d(*f)(__m128d,
			                             __m128d),
					    __m128d y0,
					    __m128d x0,
					    const __m128d h,
					    const int32_t n); 
					    
			__m128
			verner_step_xmm4r4(__m128(*f)(__m128,
			                             __m128),
					    __m128 y0,
					    __m128 x0,
					    const __m128 h,
					    const int32_t n); 
					    
					    
////////////////////////////////////////////////////////////////////////////////
//  double Runge_Kutta_Verner_Richardson( double (*f)(double, double),        //
//                      double y0, double x0, double h, int number_of_steps,  //
//                                                   int richarson_columns);  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the Runge-Kutta-Verner method described above        //
//     together with Richardson extrapolation to approximate the solution at  //
//     x = x0 + h * number_of_steps of the initial value problem y'=f(x,y),   //
//     y(x0) = y0.                                                            //
//                                                                            //
//  Arguments:                                                                //
//     double *f                                                              //
//            Pointer to the function which returns the slope at (x,y) of the //
//            integral curve of the differential equation y' = f(x,y) which   //
//            passes through the point (x0,y0).                               //
//     double y0                                                              //
//            The initial value of y at x = x0.                               //
//     double x0                                                              //
//            The initial value of x.                                         //
//     double h                                                               //
//            The step size.                                                  //
//     int    number_of_steps                                                 //
//            The number of steps. Must be nonnegative.                       //
//     int    richardson_columns                                              //
//            The maximum number of columns to use in the Richardson          //
//            extrapolation to the limit.                                     //
//                                                                            //
//  Return Values:                                                            //
//     The solution of the initial value problem y' = f(x,y), y(x0) = y0 at   //
//     x = x0 + number_of_steps * h.                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////	 



                     
			__m128d
			verner_richardson_xmm2r8(__m128d(*f)(__m128d,
			                                     __m128d),
						 __m128d y0,
						 __m128d x0,
						 const __m128d h,
						 const int32_t n,
						 int32_t n_cols); 


		       
			__m128
			verner_richardson_xmm4r4(__m128(*f)(__m128,
			                                     __m2556),
						 __m128 y0,
						 __m128 x0,
						 const __m128 h,
						 const int32_t n,
						 int32_t n_cols); 


	 
    }


}







