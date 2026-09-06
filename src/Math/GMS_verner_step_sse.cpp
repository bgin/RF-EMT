




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



/*
   Adapted from: http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_verner.html
   Manually vectorized by Bernard Gingold, beniekg@gmail.com
*/


#include "GMS_verner_step_sse.h"



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
			gms::math::verner_step_xmm2r8(__m128d(*f)(__m128d,
			                             __m128d),
					    __m128d y0,
					    __m128d x0,
					    const __m128d h,
					    const int32_t n) {
#define sqrt21 4.58257569495584000680
                               const  __m128d c1 = _mm_set1_pd(1.0 / 2.0);    
                               const  __m128d c2 = _mm_set1_pd((7.0 + sqrt21 ) / 14.0);
                               const  __m128d c3 = _mm_set1_pd((7.0 - sqrt21 ) / 14.0);

                               const  __m128d a21 =  _mm_set1_pd(1.0 / 2.0);
                               const  __m128d a31 =  _mm_set1_pd(1.0 / 4.0);
                               const  __m128d a32 =  _mm_set1_pd(1.0 / 4.0);
                               const  __m128d a41 =  _mm_set1_pd(1.0 / 7.0);
                               const  __m128d a42 =  _mm_set1_pd(-(7.0 + 3.0 * sqrt21) / 98.0);
                               const  __m128d a43 =  _mm_set1_pd((21.0 + 5.0 * sqrt21) / 49.0);
                               const  __m128d a51 =  _mm_set1_pd((11.0 + sqrt21) / 84.0);
                               const  __m128d a53 =  _mm_set1_pd((18.0 + 4.0 * sqrt21) / 63.0);
                               const  __m128d a54 =  _mm_set1_pd((21.0 - sqrt21) / 252.0);
                               const  __m128d a61 =  _mm_set1_pd((5.0 + sqrt21) / 48.0);
                               const  __m128d a63 =  _mm_set1_pd((9.0 + sqrt21) / 36.0);
                               const  __m128d a64 =  _mm_set1_pd((-231.0 + 14.0 * sqrt21) / 360.0);
                               const  __m128d a65 =  _mm_set1_pd((63.0 - 7.0 * sqrt21) / 80.0);
                               const  __m128d a71 =  _mm_set1_pd((10.0 - sqrt21) / 42.0);
                               const  __m128d a73 =  _mm_set1_pd((-432.0 + 92.0 * sqrt21) / 315.0);
                               const  __m128d a74 =  _mm_set1_pd((633.0 - 145.0 * sqrt21) / 90.0);
                               const  __m128d a75 =  _mm_set1_pd((-504.0 + 115.0 * sqrt21) / 70.0);
                               const  __m128d a76 =  _mm_set1_pd((63.0 - 13.0 * sqrt21) / 35.0);
                               const  __m128d a81 =  _mm_set1_pd(1.0 / 14.0);
                               const  __m128d a85 =  _mm_set1_pd((14.0 - 3.0 * sqrt21) / 126.0);
                               const  __m128d a86 =  _mm_set1_pd((13.0 - 3.0 * sqrt21) / 63.0);
                               const  __m128d a87 =  _mm_set1_pd(1.0 / 9.0);
                               const  __m128d a91 =  _mm_set1_pd(1.0 / 32.0);
                               const  __m128d a95 =  _mm_set1_pd((91.0 - 21.0 * sqrt21) / 576.0);
                               const  __m128d a96 =  _mm_set1_pd(11.0 / 72.0);
                               const  __m128d a97 =  _mm_set1_pd(-(385.0 + 75.0 * sqrt21) / 1152.0);
                               const  __m128d a98 =  _mm_set1_pd((63.0 + 13.0 * sqrt21) / 128.0);
                               const  __m128d a10_1 =  _mm_set1_pd(1.0 / 14.0);
                               const  __m128d a10_5 =  _mm_set1_pd(1.0 / 9.0);
                               const  __m128d a10_6 =  _mm_set1_pd(-(733.0 + 147.0 * sqrt21) / 2205.0);
                               const  __m128d a10_7 =  _mm_set1_pd((515.0 + 111.0 * sqrt21) / 504.0);
                               const  __m128d a10_8 =  _mm_set1_pd(-(51.0 + 11.0 * sqrt21) / 56.0);
                               const  __m128d a10_9 =  _mm_set1_pd((132.0 + 28.0 * sqrt21) / 245.0);
                               const  __m128d a11_5 =  _mm_set1_pd((-42.0 + 7.0 * sqrt21) / 18.0);
                               const  __m128d a11_6 =  _mm_set1_pd((-18.0 + 28.0 * sqrt21) / 45.0);
                               const  __m128d a11_7 =  _mm_set1_pd(-(273.0 + 53.0 * sqrt21) / 72.0);
                               const  __m128d a11_8 =  _mm_set1_pd((301.0 + 53.0 * sqrt21) / 72.0);
                               const  __m128d a11_9 =  _mm_set1_pd((28.0 - 28.0 * sqrt21) / 45.0);
                               const  __m128d a11_10 = _mm_set1_pd((49.0 - 7.0 * sqrt21) / 18.0);

                               const  __m128d  b1  = _mm_set1_pd(9.0 / 180.0);
                               const  __m128d  b8  = _mm_set1_pd(49.0 / 180.0);
                               const  __m128d  b9  = _mm_set1_pd(64.0 / 180.0);

                               const  __m128d c1h = _mm_mul_pd(c1,h);
			       const  __m128d c2h = _mm_mul_pd(c2,h);
			       const  __m128d c3h = _mm_mul_pd(c3,h);
			       __m128d k1,k2,k3,k4,k,k6,k7,k8,k9,k10;
			      
			       while(--n >= 0) {
                                     k1 = _mm_mul_pd(h,f(x0,y0));
				     k2 = _mm_mul_pd(h,f(_mm_add_pd(x0,c1h),
				            _mm_fmadd_pd(a21,k1,y0)));
				     k3 = _mm_mul_pd(h,f(_mm_add_pd(x0,c1h),
				            _mm_add_pd(y0,
					              _mm_fmadd_pd(a31,k1,
						                  _mm_mul_pd(a32,k2)))));
				     k4 = _mm_mul_pd(h,f(_mm_add_pd(x0,c2h),
				            		  _mm_add_pd(y0,
						                 _mm_fmadd_pd(a41,k1,
								             _mm_fmadd_pd(a42,k2,
									                 _mm_mul_pd(a43,k3))))));
				     k5 = _mm_mul_pd(h,f(_mm_add_pd(x0,c2h),
				             		   _mm_add_pd(y0,
						                 _mm_fmadd_pd(a51,k1,
								             _mm_fmadd_pd(a53,k3,
									                 _mm_mul_pd(a54,k4))))));
				     k6 = _mm_mul_pd(h,f(_mm_add_pd(x0,c1h),
				                          _mm_add_pd(y0,
							         _mm_fmadd_pd(a61,k1,
								           _mm_fmadd_pd(a63,k3,
									           _mm_fmadd_pd(a64,k4,
										              _mm_mul_pd(a65,k5)))))));
				     k7 = _mm_mul_pd(h,f(_mm_add_pd(x0,c3h),
				                          _mm_add_pd(y0,
							         _mm_fmadd_pd(a71,k1,
								       _mm_fmadd_pd(a73,k3,
								              _mm_fmadd_pd(a74,k4,
									             _mm_fmadd_pd(a75,k5,
										                _mm_mul_pd(a76,k6))))))));
				     k8 = _mm_mul_pd(h,f(_mm_add_pd(x0,c3h),
				                          _mm_add_pd(y0,
							          _mm_fmadd_pd(a81,k1,
								           _mm_fmadd_pd(a85,k5,
									           _mm_fmadd_pd(a86,k6,
										              _mm_mul_pd(a87,k7)))))));
				     k9 = _mm_mul_pd(h,f(_mm_add_pd(x0,c1h),
				                          _mm_add_pd(y0,
							         _mm_fmadd_pd(a91,k1,
								       _mm_fmadd_pd(a95,k5,
								              _mm_fmadd_pd(a96,k6,
									             _mm_fmadd_pd(a97,k7,
										                _mm_mul_pd(a98,k8))))))));
				     k10 = _mm_mul_pd(h,f(_mm_add_pd(x0,c2h),
				                           _mm_add_pd(y0,
							           _mm_fmadd_pd(a10_1,k1,
								         _mm_fmadd_pd(a10_5,k5,
									       _mm_fmadd_pd(a10_6,k6,
									             _mm_fmadd_pd(a10_7,k7,
										          _mm_fmadd_pd(a10_8,k8,
											        _mm_mul_pd(a10_9,k9)))))))));
				     x0 = _mm_add_pd(x0,h);
				     k10 = _mm_mul_pd(h,f(x0,
				                           _mm_add_pd(y0,
							           _mm_fmadd_pd(a11_5,k5,
								         _mm_fmadd_pd(a11_6,k6,
									       _mm_fmadd_pd(a11_7,k7,
									             _mm_fmadd_pd(a11_8,k8,
										          _mm_fmadd_pd(a11_9,k9,
											        _mm_mul_pd(a11_10,k10)))))))));
				     y0 = _mm_add_pd(y0,
				                    _mm_fmadd_pd(b1,k1,
						          _mm_fmadd_pd(b8,k8,
							        _mm_fmadd_pd(b9,k9,
								      _mm_fmadd_pd(b8,k10,
								            _mm_mul_pd(b1,k11))))));
			       }
			       return (y0);
		      }


		    
			__m128
			gms::math::verner_step_xmm4r4(__m128(*f)(__m128,
			                             __m128),
					    __m128 y0,
					    __m128 x0,
					    const __m128 h,
					    const int32_t n) {
#define sqrt21 4.58257569495584000680f
                               const  __m128 c1 = _mm_set1_ps(1.0f / 2.0f);    
                               const  __m128 c2 = _mm_set1_ps((7.0f + sqrt21 ) / 14.0f);
                               const  __m128 c3 = _mm_set1_ps((7.0f - sqrt21 ) / 14.0f);

                               const  __m128 a21 =  _mm_set1_ps(1.0f / 2.0f);
                               const  __m128 a31 =  _mm_set1_ps(1.0f / 4.0f);
                               const  __m128 a32 =  _mm_set1_ps(1.0f / 4.0f);
                               const  __m128 a41 =  _mm_set1_ps(1.0f / 7.0f);
                               const  __m128 a42 =  _mm_set1_ps(-(7.0f + 3.0f * sqrt21) / 98.0f);
                               const  __m128 a43 =  _mm_set1_ps((21.0f + 5.0f * sqrt21) / 49.0f);
                               const  __m128 a51 =  _mm_set1_ps((11.0f + sqrt21) / 84.0f);
                               const  __m128 a53 =  _mm_set1_ps((18.0f + 4.0f * sqrt21) / 63.0f);
                               const  __m128 a54 =  _mm_set1_ps((21.0f - sqrt21) / 252.0f);
                               const  __m128 a61 =  _mm_set1_ps((5.0f + sqrt21) / 48.0f);
                               const  __m128 a63 =  _mm_set1_ps((9.0f + sqrt21) / 36.0f);
                               const  __m128 a64 =  _mm_set1_ps((-231.0f + 14.0f * sqrt21) / 360.0f);
                               const  __m128 a65 =  _mm_set1_ps((63.0f - 7.0f * sqrt21) / 80.0f);
                               const  __m128 a71 =  _mm_set1_ps((10.0f - sqrt21) / 42.0f);
                               const  __m128 a73 =  _mm_set1_ps((-432.0f + 92.0f * sqrt21) / 315.0f);
                               const  __m128 a74 =  _mm_set1_ps((633.0f - 145.0f * sqrt21) / 90.0f);
                               const  __m128 a75 =  _mm_set1_ps((-504.0f + 115.0f * sqrt21) / 70.0f);
                               const  __m128 a76 =  _mm_set1_ps((63.0f - 13.0f * sqrt21) / 35.0f);
                               const  __m128 a81 =  _mm_set1_ps(1.0f / 14.0f);
                               const  __m128 a85 =  _mm_set1_ps((14.0f - 3.0f * sqrt21) / 126.0f);
                               const  __m128 a86 =  _mm_set1_ps((13.0f - 3.0f * sqrt21) / 63.0f);
                               const  __m128 a87 =  _mm_set1_ps(1.0f / 9.0f);
                               const  __m128 a91 =  _mm_set1_ps(1.0f / 32.0f);
                               const  __m128 a95 =  _mm_set1_ps((91.0f - 21.0f * sqrt21) / 576.0f);
                               const  __m128 a96 =  _mm_set1_ps(11.0f / 72.0f);
                               const  __m128 a97 =  _mm_set1_ps(-(385.0f + 75.0f * sqrt21) / 1152.0f);
                               const  __m128 a98 =  _mm_set1_ps((63.0f + 13.0f * sqrt21) / 128.0f);
                               const  __m128 a10_1 =  _mm_set1_ps(1.0f / 14.0f);
                               const  __m128 a10_5 =  _mm_set1_ps(1.0f / 9.0f);
                               const  __m128 a10_6 =  _mm_set1_ps(-(733.0f + 147.0f * sqrt21) / 2205.0f);
                               const  __m128 a10_7 =  _mm_set1_ps((515.0f + 111.0f * sqrt21) / 504.0f);
                               const  __m128 a10_8 =  _mm_set1_ps(-(51.0f + 11.0f * sqrt21) / 56.0f);
                               const  __m128 a10_9 =  _mm_set1_ps((132.0f + 28.0f * sqrt21) / 245.0f);
                               const  __m128 a11_5 =  _mm_set1_ps((-42.0f + 7.0f * sqrt21) / 18.0f);
                               const  __m128 a11_6 =  _mm_set1_ps((-18.0f + 28.0f * sqrt21) / 45.0f);
                               const  __m128 a11_7 =  _mm_set1_ps(-(273.0f + 53.0f * sqrt21) / 72.0f);
                               const  __m128 a11_8 =  _mm_set1_ps((301.0f + 53.0f * sqrt21) / 72.0f);
                               const  __m128 a11_9 =  _mm_set1_ps((28.0f - 28.0f * sqrt21) / 45.0f);
                               const  __m128 a11_10 = _mm_set1_ps((49.0f - 7.0f * sqrt21) / 18.0f);

                               const  __m128  b1  = _mm_set1_ps(9.0f / 180.0f);
                               const  __m128  b8  = _mm_set1_ps(49.0f / 180.0f);
                               const  __m128  b9  = _mm_set1_ps(64.0f / 180.0f);

                               const  __m128 c1h = _mm_mul_ps(c1,h);
			       const  __m128 c2h = _mm_mul_ps(c2,h);
			       const  __m128 c3h = _mm_mul_ps(c3,h);
			       __m128 k1,k2,k3,k4,k,k6,k7,k8,k9,k10;
			      
			       while(--n >= 0) {
                                     k1 = _mm_mul_ps(h,f(x0,y0));
				     k2 = _mm_mul_ps(h,f(_mm_add_ps(x0,c1h),
				            _mm_fmadd_ps(a21,k1,y0)));
				     k3 = _mm_mul_ps(h,f(_mm_add_ps(x0,c1h),
				            _mm_add_ps(y0,
					              _mm_fmadd_ps(a31,k1,
						                  _mm_mul_ps(a32,k2)))));
				     k4 = _mm_mul_ps(h,f(_mm_add_ps(x0,c2h),
				            		  _mm_add_ps(y0,
						                 _mm_fmadd_ps(a41,k1,
								             _mm_fmadd_ps(a42,k2,
									                 _mm_mul_ps(a43,k3))))));
				     k5 = _mm_mul_ps(h,f(_mm_add_ps(x0,c2h),
				             		   _mm_add_ps(y0,
						                 _mm_fmadd_ps(a51,k1,
								             _mm_fmadd_ps(a53,k3,
									                 _mm_mul_ps(a54,k4))))));
				     k6 = _mm_mul_ps(h,f(_mm_add_ps(x0,c1h),
				                          _mm_add_ps(y0,
							         _mm_fmadd_ps(a61,k1,
								           _mm_fmadd_ps(a63,k3,
									           _mm_fmadd_ps(a64,k4,
										              _mm_mul_ps(a65,k5)))))));
				     k7 = _mm_mul_ps(h,f(_mm_add_ps(x0,c3h),
				                          _mm_add_ps(y0,
							         _mm_fmadd_ps(a71,k1,
								       _mm_fmadd_ps(a73,k3,
								              _mm_fmadd_ps(a74,k4,
									             _mm_fmadd_ps(a75,k5,
										                _mm_mul_ps(a76,k6))))))));
				     k8 = _mm_mul_ps(h,f(_mm_add_ps(x0,c3h),
				                          _mm_add_ps(y0,
							          _mm_fmadd_ps(a81,k1,
								           _mm_fmadd_ps(a85,k5,
									           _mm_fmadd_ps(a86,k6,
										              _mm_mul_ps(a87,k7)))))));
				     k9 = _mm_mul_ps(h,f(_mm_add_ps(x0,c1h),
				                          _mm_add_ps(y0,
							         _mm_fmadd_ps(a91,k1,
								       _mm_fmadd_ps(a95,k5,
								              _mm_fmadd_ps(a96,k6,
									             _mm_fmadd_ps(a97,k7,
										                _mm_mul_ps(a98,k8))))))));
				     k10 = _mm_mul_ps(h,f(_mm_add_ps(x0,c2h),
				                           _mm_add_ps(y0,
							           _mm_fmadd_ps(a10_1,k1,
								         _mm_fmadd_ps(a10_5,k5,
									       _mm_fmadd_ps(a10_6,k6,
									             _mm_fmadd_ps(a10_7,k7,
										          _mm_fmadd_ps(a10_8,k8,
											        _mm_mul_ps(a10_9,k9)))))))));
				     x0 = _mm_add_ps(x0,h);
				     k10 = _mm_mul_ps(h,f(x0,
				                           _mm_add_ps(y0,
							           _mm_fmadd_ps(a11_5,k5,
								         _mm_fmadd_ps(a11_6,k6,
									       _mm_fmadd_ps(a11_7,k7,
									             _mm_fmadd_ps(a11_8,k8,
										          _mm_fmadd_ps(a11_9,k9,
											        _mm_mul_ps(a11_10,k10)))))))));
				     y0 = _mm_add_ps(y0,
				                    _mm_fmadd_ps(b1,k1,
						          _mm_fmadd_ps(b8,k8,
							        _mm_fmadd_ps(b9,k9,
								      _mm_fmadd_ps(b8,k10,
								            _mm_mul_ps(b1,k11))))));
			       }
			       return (y0);
		      }

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

#include <algorithm> 

                     
			__m128d
			gms::math::verner_richardson_xmm2r8(__m128d(*f)(__m128d,
			                                     __m128d),
						 __m128d y0,
						 __m128d x0,
						 const __m128d h,
						 const int32_t n,
						 int32_t n_cols) {

			   __ATTR_ALIGN__(32) const __m128d richardson[] = {_mm_set1_pd(1.0/255.0),
			                                                    _mm_set1_pd(1.0/511.0),
								            _mm_set1_pd(1.0/1023.0),
								            _mm_set1_pd(1.0/2047.0),
								            _mm_set1_pd(1.0/4095.0)};
			   constexpr int32_t MAX_COLS = 6;
			   __ATTR_ALIGN__(16) __m128d dt[MAX_COLS];
			   const __m128d _1_2 = _mm_set1_pd(0.5);
			   __m128d integral,delta,h_used;
			   int32_t n_subints;
			   n_cols = std::max(1,std::min(MAX_COLS,n_cols));
			   while(--n >= 0) {
                                 h_used = h;
				 n_subints = 1;
				 #pragma loop_count min(1),max(6),avg(3)
				 for(int32_t j = 0; j < n_cols; ++j) {
                                     integral = verner_step_xmm2r8(f,y0,x0,h_used,n_subints);
				     for(int32_t k = 0; k < j; ++k) {
                                         delta = _mm_sub_pd(integral,dt[k]);
					 dt[k] = integral;
					 integral = _mm_fmadd_pd(richardson[k],delta,integral);
				     }
				     dt[j] = integral;
				     h_used = _mm_mul_pd(h_used,_1_2);
				     n_subints += n_subints;
				 }
				 y0 = integral;
				 x0 = _mm_add_pd(x0,h);
			   }
			   return (y0);
		    }


		       
			__m128
			gms::math::verner_richardson_xmm4r4(__m128(*f)(__m128,
			                                     __m2556),
						 __m128 y0,
						 __m128 x0,
						 const __m128 h,
						 const int32_t n,
						 int32_t n_cols) {

			   __ATTR_ALIGN__(16) const __m128 richardson[] = { _mm_set1_ps(1.0f/255.0f),
			                                                    _mm_set1_ps(1.0f/511.0f),
								            _mm_set1_ps(1.0f/1023.0f),
								            _mm_set1_ps(1.0f/2047.0f),
								            _mm_set1_ps(1.0f/4095.0f)};
			   constexpr int32_t MAX_COLS = 6;
			   __ATTR_ALIGN__(16) __m128 dt[MAX_COLS];
			   const __m128 _1_2 = _mm_set1_ps(0.5f);
			   __m128 integral,delta,h_used;
			   int32_t n_subints;
			   n_cols = std::max(1,std::min(MAX_COLS,n_cols));
			   while(--n >= 0) {
                                 h_used = h;
				 n_subints = 1;
				 #pragma loop_count min(1),max(6),avg(3)
				 for(int32_t j = 0; j < n_cols; ++j) {
                                     integral = verner_step_xmm4r4(f,y0,x0,h_used,n_subints);
				     for(int32_t k = 0; k < j; ++k) {
                                         delta = _mm_sub_ps(integral,dt[k]);
					 dt[k] = integral;
					 integral = _mm_fmadd_ps(richardson[k],delta,integral);
				     }
				     dt[j] = integral;
				     h_used = _mm_mul_ps(h_used,_1_2);
				     n_subints += n_subints;
				 }
				 y0 = integral;
				 x0 = _mm_add_ps(x0,h);
			   }
			   return (y0);
		    }


	 
  
