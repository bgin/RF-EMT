
#ifndef __GMS_NDIM_CUBATURE_HPP__
#define __GMS_NDIM_CUBATURE_HPP__ 111220221413

/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

namespace file_version {

    const unsigned int GMS_NDIM_CUBATURE_MAJOR = 1U;
    const unsigned int GMS_NDIM_CUBATURE_MINOR = 0U;
    const unsigned int GMS_NDIM_CUBATURE2_MICRO = 0U;
    const unsigned int GMS_NDIM_CUBATURE_FULLVER =
      1000U*GMS_NDIM_CUBATURE_MAJOR+
      100U*GMS_NDIM_CUBATURE_MINOR+
      10U*GMS_NDIM_CUBATURE_MICRO;
    const char * const GMS_NDIM_CUBATURE_CREATION_DATE = "11-12-2022 14:13 AM +00200 (SUN 11 DEC 2022 GMT+2)";
    const char * const GMS_NDIM_CUBATURE_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_NDIM_CUBATURE_AUTHOR        = "Programmer: Adapted (from the original work) by Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_NDIM_CUBATURE_DESCRIPTION   = "N-dimensional cubatures by Arthur Stroud."

}


#include <valarray>
#include <memory>
#include <utility>
#include <math.h>

#ifndef FUSE_FOR_LOOP
#define FUSE_FOR_LOOP 0x1
#endif

namespace gms {
	namespace math {

		

		// Global apply functions.

		template<class Fx, typename Args>
		auto apply(_In_ Fx fx, _In_ Args&& args)->std::result_of<Fx(Args)> { return (fx(args)); }

		template<template<class,class> BinOp,
		         class Fx1, class Fx2,
					 typename Args1, typename Args2> auto apply(_In_ BinOp<class,class> &op,
					                                            _In_ Fx1 fx1,
					                                            _In_ Args1 &&args1,
					                                            _In_ Fx2 fx2,
					                                            _In_ Args2 &&args2)->std::result_of<BinOp<Fx1(Args1), Fx2(Args2)>> {
						 return (op.operator()(fx1(args1), fx2(args2)));
					 }

		template<template<class,class,class> TernOp,
		                 class Fx1, class Fx2, class Fx3,
				           typename Args1, typename Args2, typename Args3>
				                 auto apply(_In_ TernOp<class, class, class> &op,
				                            _In_ Fx1 fx1,
				                            _In_ Args1 &&args1,
				                            _In_ Fx2 fx2,
				                            _In_ Args2 &&args2,
				                            _In_ Fx3 fx3,
				                            _In_ Args3 &&args3)->std::result_of<TernOp<Fx1(Args1), Fx2(Args2), Fx3(Args3)>> {

					 return (op.operator()(fx1(args1), fx2(args2), fx3(args3)));
				 }


		

		// Cubature constants.

		template<typename Real_t,
			typename = std::enable_if<std::is_floating_point<Real_t>
			                             ::value>::type> struct CubConsts {

				  static constexpr Real_t _0_0_{ static_cast<Real_t>(0.0L) };

				  static constexpr Real_t _1_0_{ static_cast<Real_t>(1.0L) };

				  static constexpr Real_t _2_0_{ static_cast<Real_t>(2.0L) };

				  static constexpr Real_t _4_0_0{ static_cast<Real_t>(4.0L) };

				  static constexpr Real_t _8_0_{ static_cast<Real_t>(8.0L) };

				  static constexpr Real_t _pi_{ static_cast<Real_t>(3.1415926535897932384626433832795L) };

				  static constexpr Real_t _3_div_16_{ static_cast<Real_t>(0.1875L) };

				  static constexpr Real_t _1_div_3360_{ static_cast<Real_t>(2.9761904761904761904761904761905e-4L) };

				  static constexpr Real_t _3_div_32_{ static_cast<Real_t>(0.09375L) };

				  static constexpr Real_t _pi_div_8_{ static_cast<Real_t>(0.39269908169872415480783042290994L) };

				  static constexpr Real_t _4_div_3_{ static_cast<Real_t>(1.3333333333333333333333333333333L) };

				  static constexpr Real_t _2_div_3_{ static_cast<Real_t>(0.66666666666666666666666666666667L) };

				  static constexpr Real_t _1_div_8_{ static_cast<Real_t>(0.125L) };

				  static constexpr Real_t _1_div_4_{ static_cast<Real_t>(0.25L) };

				  static constexpr Real_t _sqrt2_div_2_{ static_cast<Real_t>(0.7071067811865475L) };

				  static constexpr Real_t _1_div_6_{ static_cast<Real_t>(0.16666666666666666666666666666667L) };

				  static constexpr Real_t _4_div_6_{ static_cast<Real_t>(0.66666666666666666666666666666667L) };

				  static constexpr Real_t _1_div_52_{ static_cast<Real_t>(0.01923076923076923076923076923077L) };

				  static constexpr Real_t _1_div_1566_{ static_cast<Real_t>(6.3856960408684546615581098339719e-4L) };
                  };

		//****************************************************************************
		//
		//  Purpose:
		//
		//    BALL_F1_ND approximates an integral inside a ball in ND.
		//
		//  Integration region:
		//
		//    sum ( X(1:N) - CENTER(1:N) )^2 <= R * R.
		//
		//  Discussion:
		//
		//    An (N+1)*2^N point 5-th degree formula is used, Stroud number SN:5-6.
		//
		//  
		//
		//  
		//
		//  Author:
		//
		//    Bernard Gingold who adapted this implementation from the work of original author
		//    John Burkardt.
		//  
		//   
		//
		//  Reference:
		//
		//    Arthur Stroud,
		//    Approximate Calculation of Multiple Integrals,
		//    Prentice Hall, 1971,
		//    ISBN: 0130438936,
		//    LC: QA311.S85.
		///////////////////////////////////////////////////////////     

		template<class I4_POWER, class BALL_VOLUME_ND, typename Real_t,
			     typename = std::enable_if<std::is_floating_point<Real_t>
				                             ::value>::type> struct Ball_f1_nd {
                  

				 
			    
				 static  Real_t  ball_f1_nd(_In_ Real_t(*Fx)(const int n, std::unique_ptr<Real_t,void(*)(Real_t*)> &x),
										    _In_ const int n, 
											_In_ std::unique_ptr<Real_t,void(*)(Real_t*)> &center,
										    _In_ const Real_t r) {


					 int  ihi{},  k{}, khi{}, ktemp{}, itemp{};
					 Real_t quad{}, result{}, t{};
					 Real_t temp{}, u{}, u2{}, v{}, volume{}, w{}, y{};
					
					 
					 Real_t y{};

					 if (r == CubConsts<Real_t, typename
						      std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_) {
						 result = CubConsts<Real_t, typename
							     std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_; 
						 return (result);
					 }
					 std::unique_ptr<Real_t, void(*)(Real_t*)> x{ new Real_t[n], [](Real_t *p)->void {delete[]p;} };
					 u2 = (CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ - CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_2_0_ * std::sqrt(CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ / (static_cast<Real_t>(n + 4) ))) / 
						 static_cast<Real_t>(n + 2);
					 u = std::sqrt(u2);
					 for (int i{ 0 }; i != n; ++i) {
						 (&*x)[i] = (&*center)[i] - r * u;
					 }
					 w = CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ / (static_cast<Real_t>((n + 1) * 
						 I4_POWER::i4_power(2,n)));
					 ihi = I4_POWER::i4_power(2, n);
					 for (int i{ 0 }; i != ihi; ++i) {
						 itemp = 1;
						 for (int j{ 0 }; j != n; ++j) {
							 u = ((&*center)[j] - (&*x)[j]) / r;
							 if ((itemp % 2) == 1) {
								 (&*x)[j] = (&*center)[j] - std::fabs((&*x)[j] - (&*center)[j]);
							 }
							 else {
								 (&*x)[j] = (&*center)[j] + std::fabs((&*x)[j] - (&*center)[j]);
							 }
							 itemp /= 2;
						 }
						 quad += w * Fx(n,x);
					 }
					 temp = std::sqrt(static_cast<Real_t>(n + 4));
					 t = std::sqrt(CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_2_0_ * 
						 static_cast<Real_t>(n + 1) / static_cast<Real_t>(n + 2)) / 
						 static_cast<Real_t>(n) * temp;
					 y = (CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ + 
						 CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_2_0_ / (
						 static_cast<Real_t>(n)* temp)) / static_cast<Real_t>(n + 2);
					 v = std::sqrt(y - t);
					 u = std::sqrt(y + static_cast<Real_t>(n - 1) * t);

					 khi = I4_POWER::i4_power(2, n);

					 for (int i{ 0 }; i != n; ++i) {
						 for (int j{ 0 }; j != n; ++j) {
							 (&*x)[j] = (&*center)[j] - r * v;
						 }
						 (&*x)[i] = (&*center)[i] - r * v;

						 for (int k{ 0 }; k != khi; ++k) {
							 ktemp = k;
							 for (int j{ 0 }; j != n; ++j) {
								 if ((ktemp % 2) == 1) {
									 (&*x)[j] = (&*center)[j] - std::fabs((&*x)[j] - (&*center)[j]);
								 }
								 else{
									 (&*x)[j] = (&*center)[j] + std::fabs((&*x)[j] - (&*center)[j]);
								 }
								 ktemp /= 2;
							 }
							 quad += w * Fx(n,x);
						 }
						 (&*x)[i] = (&*center)[i] - r * v;
					 }

					 volume = BALL_VOLUME_ND::ball_volume_nd(n,r);
					 result = quad * volume;

					 return (result);
			}

			// overloaded valarray version
				 
				 static Real_t  ball_f1_nd(_In_ Real_t(*Fx)(std::size_t n, std::valarray<Real_t>&),
					                       _In_ const std::size_t n,
					                       _In_ const std::valarray<Real_t> &center,
					                       _In_ const Real_t r) {

					 std::size_t  
					       ihi{}, 
						   k{}, 
						   khi{}, 
						   ktemp{};
					 Real_t quad{}, 
					        result{}, 
							t{},
					        temp{}, 
							u{}, 
							u2{},
							v{},
							volume{},
						    w{}, 
							y{};

					 Real_t y{};

					 if (r == CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_) {
						 result = CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;
						 return (result);
					 }
					 std::valarray<Real_t> x;
					 x.resize(n);
					 u2 = (CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ - CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_2_0_ * std::sqrt(CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ / (static_cast<Real_t>(n + 4)))) /
						 static_cast<Real_t>(n + 2);
					 u = std::sqrt(u2);

					 for (std::size_t i{ 0ULL }; i != n; ++i) {
						 x.operator[](i) = center.operator[](i) - r * v;
					 }
					 w = CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ / (static_cast<Real_t>((n + 1) *
						 I4_POWER::i4_power(2, n)));
					 ihi = I4_POWER::i4_power(2, n);

					 for (std::size_t i{ 0ULL }; i != n; ++i) {
						 itemp = i;
						 for (std::size_t j{ 0ULL }; j != n; ++j) {
							 u = (center.operator[](j) - x.operator[](j)) / r;
							 if ((itemp % 2) == 1) {
								 x.operator[](j) = center.operator[](j) - std::fabs(x[j] - center[j]);
							 }
							 else {
								 x.operator[](j) = center.operator[](j) + std::fabs(x[j] - center[j]);
							 }
							 itemp /= 2;
						 }
						 quad += w  * Fx(n,x);
					 }

					 temp = std::sqrt(static_cast<Real_t>(n + 4));
					 t = std::sqrt(CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_2_0_ *
						 static_cast<Real_t>(n + 1) / static_cast<Real_t>(n + 2)) /
						 static_cast<Real_t>(n)* temp;
					 y = (CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ +
						 CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_2_0_ / (
						 static_cast<Real_t>(n)* temp)) / static_cast<Real_t>(n + 2);
					 v = std::sqrt(y - t);
					 u = std::sqrt(y + static_cast<Real_t>(n - 1) * t);

					 khi = I4_POWER::i4_power(2, n);

					 for (std::size_t i{ 0ULL }; i != n; ++i) {
						 for (std::size_t j{ 0ULL }; j != n; ++j) {
							 x[j] = center[j] - r * v;
						 }
						  x[i] =  center[i] - r * v;

						 for (std::size_t k{ 0ULL }; k != khi; ++k) {
							 ktemp = k;
							 for (std::size_t j{ 0ULL }; j != n; ++j) {
								 if ((ktemp % 2) == 1) {
									 x[j] = center[j] - std::fabs(x[j] - center[j]);
								 }
								 else{
									 x[j] = center[j] + std::fabs(x[j] - center[j]);
								 }
								 ktemp /= 2;
							 }
							 quad += w * Fx(n, x);
						 }
						x[i] = center[i] - r * v;
					 }

					 volume = BALL_VOLUME_ND::ball_volume_nd(n, r);
					 result = quad * volume;

					 return (result);
				 }
};

//****************************************************************************
//
//  Purpose:
//
//    BALL_F3_ND approximates an integral inside a ball in ND.
//
//  Integration region:
//
//    sum ( X(1:N) - CENTER(1:N) )^2 <= R * R.
//
//  Discussion:
//
//    A 2**(N+1)-1 point 5-th degree formula is used, Stroud number SN:5-4.
//
//  
//  Author:
//
//    Adapted by Bernard Gingold from the John Burkardt implementation.
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the user supplied
//    function which evaluates F at the N-vector X.
//
//    Input, int N, the dimension of the space.
//
//    Input, double CENTER[N], the center of the ball.
//
//    Input, double R, the radius of the ball.
//
//    Output, double BALL_F3_ND, the approximate integral of the function.
/////////////////////////////////////////////////////////////////////////////////

		template<class I4_POWER,class BALL_VOLUME_ND,
		                              typename Real_t,
				    typename = std::enable_if<std::is_floating_point<Real_t>
				                           ::value>::type> struct Ball_f3_nd {

					 //std::unique_ptr version.
					 static Real_t  ball_f3_nd(_In_ Real_t(*Fx)(std::size_t n, std::unique_ptr<Real_t, void(*)(Real_t*)> &x),
											   _In_ std::size_t n,
											   _In_ std::unique_ptr<Real_t, void(*)(Real_t*)> &center,
											   _In_ const Real_t r) {

						 int   jhi{}, 
						       jtemp{}, 
							   k{};
						 Real_t quad{},
							    result{},
							    ri{},
							    s{},
							    volume{},
							    weight{};

						 if (r == CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_) {
							 result = CubConsts<Real_t, typename
								 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;
							 return (result);
						 }

						 std::unique_ptr<Real_t, void(*)(Real_t*)> x{ new Real_t[n], [](Real_t *p)->void {delete[]p;} };
						 quad = CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;

						 //  The first point is the center of the ball.
						 for (int i{ 0 }; i != n; ++i) {
							 (&*x)[i] = (&*center)[i];
						 }
						 weight = CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_4_0_ / 
							 static_cast<Real_t>(I4_POWER::i4_power(n+2,2));
						 quad += weight * Fx(n,x);

						 s = CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ / 
							 std::sqrt(static_cast<Real_t>(n + 4));

						 for (int i{ 0 }; i != n; ++i) {
							 ri = std::sqrt(static_cast<Real_t>(i + 3) / static_cast<Real_t>(n + 4));
							 //  Set up the first point, with (I) zeroes, RI, and then N-I-1 S's.
							 for (int j{ 0 }; j != n; ++j) {
								 if (j < i) {
									 (&*x)[j] = (&*center)[j];
								 }
								 else if (j == i) {
									 (&*x)[j] = (&*center)[j] + r * v;
								 }
								 else {
									 (&*x)[j] = (&*center)[j] + r * s;
								 }
							 }
							 weight = std::pow(CubConsts<Real_t, typename
								 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_2_0_, i + 1 - n) * 
								 static_cast<Real_t>(n + 4) / static_cast<Real_t>((i + 2) * (i + 3) * (n + 2));
							 //  Now go through all sign permutations of the basic point.
							 jhi = I4_POWER::i4_power(2, n - i);
							 for (int j{ 0 }; j != jhi; ++j) {
								 jtemp = j;
								 for (int k{ 0 }; k != n; ++k) {
									 
									 (jtemp % 2) == 1 ? (&*x)[k] = (&*center)[k] - std::fabs((&*x)[k] - (&*center)[k]) : 
										                (&*x)[k] = (&*center)[k] - std::fabs((&*x)[k] - (&*center)[k]);
									 jtemp /= 2;
								 }
								 quad += weight * Fx(n, x);
							 }
							 
						 }
						 volume = BALL_VOLUME_ND::ball_volume_nd(n,r);
						 result = quad * volume;
						 return (result);
							
	   }


					 // valarray based version.
					 static Real_t 	ball_f3_nd(_In_ Real_t(*Fx)(std::size_t n, std::valarray<Real_t> &x),
						                       _In_ const std::size_t n,
						                       _In_ const std::valarray<Real_t> &center,
											   _In_ const Real_t r) {

						 std::size_t jhi{},
							         jtemp{},
							         k{};

						 Real_t quad{},
							 result{},
							 ri{},
							 s{},
							 volume{},
							 weight{};

						 if (r == CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_) {
							 result = CubConsts<Real_t, typename
								 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;
							 return (result);
						 }
						 std::valarray<Real_t> x;
						 x.resize(n);
						 quad = CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;

						 //  The first point is the center of the ball.
						 for (std::size_t i{ 0ULL }; i != n; ++i) {
							 x[i] = center[i];
						 }
						 weight = CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_4_0_ /
							 static_cast<Real_t>(I4_POWER::i4_power(n + 2, 2));
						 quad += weight * Fx(n, x);

						 s = CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ /
							 std::sqrt(static_cast<Real_t>(n + 4));

						 for (std::size_t i{ 0ULL }; i != n; ++i) {
							 ri = std::sqrt(static_cast<Real_t>(i + 3) / static_cast<Real_t>(n + 4));
							 //  Set up the first point, with (I) zeroes, RI, and then N-I-1 S's.
							 for (std::size_t j{ 0 }; j != n; ++j) {
								 
								 (j < i) ? x[j] = center[j] : (j == i ? x[j] = center[j] + r * v : 
									                           x[j] = center[j] + r * s);
							 }
							 weight = std::pow(CubConsts<Real_t, typename
								 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_2_0_, i + 1 - n) *
								 static_cast<Real_t>(n + 4) / static_cast<Real_t>((i + 2) * (i + 3) * (n + 2));
							 //  Now go through all sign permutations of the basic point.
							 jhi = I4_POWER::i4_power(2, n - i);
							 for (std::size_t j{ 0ULL }; j != jhi; ++j) {
								 jtemp = j;
								 for (std::size_t k{ 0 }; k != n; ++k) {

									 (jtemp % 2) == 1 ? x[k] = center[k] - std::fabs( x[k] - center[k]) :
										 *x[k] = center[k] - std::fabs( x[k] - center[k]);
									 jtemp /= 2;
								 }
								 quad += weight * Fx(n, x);
							 }

						 }
						 volume = BALL_VOLUME_ND::ball_volume_nd(n,r);
						 result = quad * volume;
						 return (result);

					 }
};

///////////////////////////////////////////////////////////////////////////////
//
//  Purpose:
//
//    BALL_MONOMIAL_ND integrates a monomial on a ball in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) <= R * R
//
//  
//
//   
//
//  License:
//         
//         This code is distributed under the GNU LGPL license. 
//
//  Author:
//
//    Adapted by Bernard Gingold from the John Burkardt implementation.
//
//  Reference:
//
//    Gerald Folland,
//    How to Integrate a Polynomial Over a Sphere,
//    American Mathematical Monthly,
//    Volume 108, Number 5, May 2001, pages 446-448.
//
//  Parameters:
//
//    Input, int N, the dimension of the space.
//
//    Input, int P[N], the exponents of X(1) through X(N) in the monomial.
//    The exponents must be nonnegative.
//
//    Input, double R, the radius of the ball.
///////////////////////////////////////////////////////////////////////////////				                                              
		        

		template<class SPHERE_UNIT_MONOMIAL_ND,
			typename Real_t, typename = std::enable_if<
			   std::is_floating_point<Real_t>::value>::type> struct Ball_monomial_nd {
								
                        
					//std::unique_ptr version.
				   static Real_t ball_monomial_nd(_In_ const int n,
				                                  _In_ std::unique_ptr<int, void(*)(int*)> &p,
												  _In_ const Real_t r) {
                       
					   
					   Real_t power{},
						      value{};

					   power = static_cast<Real_t>(n);
					   for (int i{ 0 }; i != n; ++i) {
						   power += static_cast<Real_t>((&*p)[i]);
					   }
					   value = SPHERE_UNIT_MONOMIAL_ND::sphere_unit_monomial_nd(n,p) * 
						       std::pow(r,power) / power;

					   return (value);
		 }

		           //valarray version.
				   static Real_t ball_monomial_nd(_In_ const std::size_t n,
					                              _In_ const std::valarray<int> &p,
												  _In_ const Real_t r) {

					   Real_t power{},
						      value{};
					   power = static_cast<Real_t>(n);
					   for (std::size_t i{ 0ULL }; i != n; ++i){
						   power += static_cast<Real_t>(p[i]);
					   }
					   value = SPHERE_UNIT_MONOMIAL_ND::sphere_unit_monomial_nd(n,p) * 
						   std::pow(r,power) / power;

					   return (value);
		}

};

///////////////////////////////////////////////////////////////////////////////
//
//  Purpose:
//
//    BALL_UNIT_07_3D approximates an integral inside the unit ball in 3D.
//
//  Integration region:
//
//    X*X + Y*Y + Z*Z <= 1.
//
//  Discussion:
//
//    A 64 point 7-th degree formula is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  
//
//    
//
//  Author:
//
//    Adapted by Bernard Gingold from the John Burkardt implementation.
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of 
//    the user supplied function which evaluates F(X,Y,Z).
//
//    Output, double BALL_UNIT_07_3D, the approximate integral of the function.
///////////////////////////////////////////////////////////////////////////////

		template<class LEGENDRE_SET,
		         class BALL_UNIT_VOLUME_3D,
			     typename Real_t,
			     typename = std::enable_if<std::is_floating_point<Real_t>
				                              ::value>::type> struct Ball_unit_07_3d {

				
				// Only single function, no overloded instances.
				   static Real_t ball_unit_07_3d(_In_ Real_t(*Fx)(_In_ Real_t x, _In_ Real_t y, _In_ Real_t z)) {

					   constexpr int order{ 4 };
					   constexpr int size{4};
					   Real_t angle{},
						   quad{},
						   result{},
						   volume{},
						   w{},
						   x{},
						   y{},
						   x{};
					  
					   constexpr Real_t weight1[size] = { Real_t(0.19455533421780251826L),
						                                  Real_t(0.13877799911553081506L),
													      Real_t(0.13877799911553081506L),
													      Real_t(0.19455533421780251826L) };

					   Real_t weight2[size];
					   Real_t weight3[size];
					   constexpr Real_t xtab1[size] = { Real_t(-0.906179845938663992797626878299L),
						                                Real_t(-0.538469310105683091036314420700L),
														Real_t(0.538469310105683091036314420700L),
														Real_t(0.906179845938663992797626878299L) };
					   Real_t xtab2[size];
					   Real_t xtab3[size];

					   // Set xtab2 and weight2.
					   for (int j{ 0 }; j != order; ++j) {
						   angle = CubConsts<Real_t, typename
							   std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_ * 
							   static_cast<Real_t>(2 * j - 1) / static_cast<Real_t>(2 * order);
					   }

					   for (int j{ 0 }; j != order; ++j) {
						   weight2[j] = CubConsts<Real_t, typename
							   std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_;
					   }

					   // Set xtab3 and weight3 for the interval [-1,1].
					   LEGENDRE_SET::legendre_set(order,xtab3,weight3);
					   w = CubConsts<Real_t, typename
						   std::enable_if<std::is_floating_point<Real_t>::value>::type>::_3_div_16_;
					   quad = CubConsts<Real_t, typename
						   std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;

					   for (int i{ 0 }; i != order; ++i) {
						   for (int j{ 0 }; j != order; ++j) {
							   for (int k{ 0 }; k != order; ++k) {
								   x = xtab1[i] - std::sqrt(CubConsts<Real_t, typename
									   std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ - xtab2[j] * xtab2[j]) * 
									   std::sqrt(CubConsts<Real_t, typename
									   std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ - xtab3[k] * xtab3[k]);
								   y = xtab1[i] * xtab2[j] * std::sqrt(CubConsts<Real_t, typename
									   std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ - xtab3[k] * xtab3[k]);
								   z = xtab1[i] * xtab3[k];

								   quad += w * weight1[i] * weight2[j] * weight3[k] * Fx(x,y,z);
							   
							   }
						   }
					   }
					   volume = BALL_UNIT_VOLUME_3D::ball_unit_volume_3d();
					   result = quad * volume;
					   return (result);
					   
	    }
};

///////////////////////////////////////////////////////////////////////////////
//
//  Purpose:
//
//    BALL_UNIT_14_3D approximates an integral inside the unit ball in 3D.
//
//  Integration region:
//
//    X*X + Y*Y + Z*Z <= 1.
//
//  Discussion:
//
//    A 288 point 14-th degree formula is used, Stroud number S3:14-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  
//
//   
//
//  Author:
//
//    Adapted by Bernard Gingold from the John Burkardt implementation.
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of 
//    the user supplied function which evaluates F(X,Y,Z).
//
//    Output, double RESULT, the approximate integral of the function.
///////////////////////////////////////////////////////////////////////////////

		template<class BALL_UNIT_VOLUME_3D,
		         typename Real_t,
				 typename = std::enable_if<std::is_floating_point<Real_t>
				                              ::value>::type> struct Ball_unit_14_3d {

				 
				 // Only single instance.
				 static Real_t ball_unit_14_3d(_In_ Real_t(*Fx)(_In_ Real_t x, _In_ Real_t y, _In_ Real_t z)) {

					 constexpr int size{ 4 };
					 Real_t quad{},
						    result{},
						    temp{},
						    volume{},
						    w1{},
						    w2{},
						    x{},
						    y{},
							z{};
					 constexpr Real_t xval{ Real_t(0.525731112) };
					 constexpr Real_t yval{ Real_t(0.850650808) };
					 constexpr Real_t w1val{ Real_t(125.0) };
					 constexpr Real_t w2val{ Real_t(143.0) };
					 constexpr Real_t r[size] =       {  Real_t(0.968160240),
						                                 Real_t(0.836031107),
												         Real_t(0.613371433),
												         Real_t(0.324253423) };
					 constexpr Real_t weight[size] = {   Real_t(0.076181268),
						                                 Real_t(0.126263673),
													     Real_t(0.098048133),
													     Real_t(0.032840260) };
					 constexpr Real_t xtab[size + 1] = { Real_t(-0.151108275),
						                                 Real_t(0.315838353),
														 Real_t(0.346307112),
														 Real_t(-0.101808787),
														 Real_t(-0.409228403) };
					 constexpr Real_t ytab[size + 1] = { Real_t(0.155240600),
						                                 Real_t(0.257049387),
														 Real_t(0.666277790),
														 Real_t(0.817386065),
														 Real_t(0.501547712) };
					 constexpr Real_t ztab[size + 1] = { Real_t(0.976251323),
						                                 Real_t(0.913330032),
														 Real_t(0.660412970),
														 Real_t(0.567022920),
														 Real_t(0.762221757) };

					 quad = CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;

					 for (int m{ 0 }; m != size; ++m) {
						 w1 = w1val * weight[m] * CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_div_3360_;
						 x = xval * r[m];
						 y = zval * r[m];
						 z = CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;
						 for (int j{ 0 }; j != 2; ++j){
							 x = -x;
							 for (int k{ 0 }; k != 2; ++k) {
								 y = -y;
								 for (int l{ 0 }; l != 3; ++l) {
									 temp = z;
									 z = y;
									 y = x;
									 x = temp;
									 quad += w1 * Fx(x,y,z);
								 }
							 }
						 }
						 w2 = w2val * weight[m] * CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_div_3360_;
						 for (int n{ 0 }; n != 5; ++n) {
							 x = xtab[n] * r[m];
							 y = ytab[n] * r[m];
							 z = ztab[n] * r[m];
							 for (int i{ 0 }; i != 3; ++i) {
								 temp = x;
								 x = z;
								 z = -y;
								 y = -temp;
								 for (int j{ 0 }; j != 3; ++j) {
									 temp = z;
									 z = y;
									 y = x;
									 x = temp;

									 quad += w2 * Fx(x,y,z);
								 }
								 y = -y;
								 z = -z;
								 quad += w2 * Fx(x,y,z);
							 }
						 }
					 }

					 volume = BALL_UNIT_VOLUME_3D::ball_unit_volume_3d();
					 result = quad * volume;
					 return (result);
		 }
};

///////////////////////////////////////////////////////////////////////////////
//
//  Purpose:
//
//    BALL_UNIT_15_3D approximates an integral inside the unit ball in 3D.
//
//  Integration region:
//
//    X * X + Y * Y + Z * Z <= 1.
//
//  Discussion:
//
//    A 512 point 15-th degree formula is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  
//
//   
//
//  Author:
//
//    Adapted by Bernard Gingold from the John Burkardt implementation.
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function which evaluates F(X,Y,Z).
//
//    Output, double BALL_UNIT_15_3D, the approximate integral of the function.
///////////////////////////////////////////////////////////////////////////////

		template<class LEGENDRE_SET,
		         class BALL_UNIT_VOLUME_3D,
				 typename Real_t,
				 typename = std::enable_if<std::is_floating_point<Real_t>
				                        ::value>::type> struct Ball_unit_15_3d{


				 // Single instance only.
				 static Real_t ball_unit_15_3d(_In_ Real_t(*Fx)(_In_ const Real_t x, 
				                                                _In_ const Real_t y, 
			                                                        _In_ const Real_t z)) {

					 constexpr int order1{ 4 };
					 constexpr int order2{ 8 };
					 constexpr int _16_iter_{16};
					 Real_t cj{},
						 ck{},
						 quad{},
						 result{},
						 sj{},
						 sk{},
						 volume{},
						 w{},
						 x{},
						 y{},
						 z{};
					 
					 constexpr Real_t weight[order1] = { Real_t(0.0328402599),
						                                 Real_t(0.0980481327),
														 Real_t(0.1262636728),
														 Real_t(0.0761812678) };
					 constexpr Real_t xtab1[order1] = {  Real_t(0.3242534234),
						                                 Real_t(0.6133714327),
														 Real_t(0.8360311073),
														 Real_t(0.9681602395) };
					 Real_t weight2[order2] = { CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_ };
					 Real_t   xtab2[order2] = { CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_ };

					 LEGENDRE_SET::legendre_set(order2,xtab2,weight2);

					 w = CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_3_div_32_;
					 quad = CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;

					 for (int i{ 0 }; i != order1; ++i) {
						 for (int j{ 0 }; j != order2; ++j) {
							 sj = xtab2[j];
							 cj = std::sqrt(CubConsts<Real_t, typename
								 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ - 
								 sj * sj);
							 for (int k{ 1 }; k != _16_iter_; ++k) {
								 sk = std::sin(static_cast<Real_t>(k)*CubConsts<Real_t, typename
									 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_div_8_);
								 ck = std::cos(static_cast<Real_t>(k)* CubConsts<Real_t, typename
									 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_div_8_);
								 x = xtab1[i] * cj * ck;
								 y = xtab1[i] * cj * sk;
								 z = xtab1[i] * sj;
								 quad += w * weight1[i] * weight2[j] * Fx(x,y,z);
							 }
						 }
					 }

					 volume = BALL_UNIT_VOLUME_3D::ball_unit_volume_3d();
					 result = quad * volume;
					 return (result);

		}
};

///////////////////////////////////////////////////////////////////////////////
//
//  Purpose:
//
//    BALL_UNIT_F1_ND approximates an integral inside the unit ball in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) <= 1.
//
//  Discussion:
//
//    An (N+1)*2^N point 5-th degree formula is used, Stroud number SN:5-6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  
//
//   
//
//  Author:
//
//    Adapted by Bernard Gingold from John Burkardt implementation.
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the 
//    user supplied function which evaluates F at the N-vector X.
//
//    Input, int N, the dimension of the space.
//
//    Output, double BALL_UNIT_F1_ND, the approximate integral of the function.
///////////////////////////////////////////////////////////////////////////////

		template<class I4_POWER,
		         class BALL_UNIT_VOLUME_ND,
		         typename Real_t,
				 typename = std::enable_if<std::is_floating_point<Real_t>
				                              ::value>::type> struct Ball_unit_f1_nd {

                  template<typename Real_t> using V = std::valarray<Real_t>;
				 // std::unique_ptr version.
				 static Real_t ball_unit_f1_nd(_In_ Real_t(*Fx)(_In_ const int n,
											   _In_ std::unique_ptr<Real_t, void(*)(Real_t*)> &x),
											   _In_ const int n) {

					 int ihi{},
						 itemp{},
						 khi{},
						 ktemp{};
					 Real_t quad{},
						 result{},
						 t{},
						 temp{},
						 u{},
						 u2{},
						 v{},
						 volume{},
						 w{},
						 y{};
					 
					 std::unique_ptr<Real_t, void(*)(Real_t*)> x{ new Real_t[n], [](Real_t *p)->void {delete[]p;} };
					 u2 = (CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ - 
						 CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_2_0_ * std::sqrt(CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ / static_cast<Real_t>(n+1))) / 
						 static_cast<Real_t>(n + 2);
					 u = std::sqrt(u2);
					 for (int i{ 0 }; i != n; ++i) {
						 (&*x)[i] = -u;
					 }
					 w = CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ / 
						 static_cast<Real_t>(n + 1) * I4_POWER::i4_power(2,n);
					 quad = CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;
					 ihi = I4_POWER::i4_power(2, n);

					 for (int i{ 0 }; i != ihi; ++i) {
						 itemp = i;
						 for (int j{ 0 }; j != n; ++j) {
							 (itemp % 2) == 1 ? (&*x)[j] = -std::fabs((&*x)[j]) : 
								                (&*x)[j] = std::fabs((&*x)[j]);
							 itemp /= 2;

						 }
						 quad += w * Fx(n,x);
					 }
					 temp = std::sqrt(static_cast<Real_t>(n + 4));
					 t = std::sqrt(CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_2_0_ * 
						 static_cast<Real_t>(n + 1) / static_cast<Real_t>(n + 2)) / 
						 static_cast<Real_t>(n) * temp;
					 y = (CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ + 
						 CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_2_0_ / 
						 static_cast<Real_t>(n)* temp) / static_cast<Real_t>(n + 2);
					 v = std::sqrt(y - t);
					 u = std::sqrt(y + static_cast<Real_t>(n - 1) * t);
					 khi = I4_POWER::i4_power(2, n);

					 for (int i{ 0 }; i != n; ++i) {
						 for (int j{ 0 }; j != n; ++j) {
							 (&*x)[j] = -v;
						 }
						 (&*x)[i] = -u;

						 for (int k{ 0 }; k != khi; ++k) {
							 ktemp = k;
							 for (int j{ 0 }; j != n; ++j) {
								 (ktemp % 2) == 2 ? (&*x)[j] = -std::fabs((&*x)[j]) : 
									                (&*x)[j] = std::fabs((&*x)[j]);
								 ktemp /= 2;
							 }
							 quad += w * Fx(n,x);
						 }
						 (&*x)[i] = -v;
					 }

					 volume = BALL_UNIT_VOLUME_ND::ball_unit_volume_nd();
					 result = quad * volume;
					 return (result);
		 }

		          //std::valarray version.
				 static Real_t ball_unit_f1_nd(_In_ Real_t(*Fx)(_In_ const std::size_t n,
					                           _In_ const V &x),
					                           _In_ const std::size_t n) {

			 std::size_t ihi{},
						 itemp{},
						 khi{},
						 ktemp{};
					 Real_t quad{},
						 result{},
						 t{},
						 temp{},
						 u{},
						 u2{},
						 v{},
						 volume{},
						 w{},
						 y{};

					 V<Real_t> x;
					 x.resize(n);
					 u2 = (CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ -
						 CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_2_0_ * std::sqrt(CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ / static_cast<Real_t>(n + 1))) /
						 static_cast<Real_t>(n + 2);
					 u = std::sqrt(u2);
					 for (std::size_t i{ 0 }; i != n; ++i) {
						 x[i] = -u;
					 }
					 w = CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ /
						 static_cast<Real_t>(n + 1) * I4_POWER::i4_power(2, n);
					 quad = CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;
					 ihi = I4_POWER::i4_power(2, n);
					 for (std::size_t i{ 0 }; i != ihi; ++i) {
						 itemp = i;
						 for (std::size_t j{ 0 }; j != n; ++j) {
							 (itemp % 2) == 1 ? x[j] = -std::fabs(x[j]) :
								 x[j] = std::fabs(x[j]);
							 itemp /= 2;

						 }
						 quad += w * Fx(n, x);
					 }
					 temp = std::sqrt(static_cast<Real_t>(n + 4));
					 t = std::sqrt(CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_2_0_ *
						 static_cast<Real_t>(n + 1) / static_cast<Real_t>(n + 2)) /
						 static_cast<Real_t>(n) * temp;
					 y = (CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ +
						 CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_2_0_ /
						 static_cast<Real_t>(n)* temp) / static_cast<Real_t>(n + 2);
					 v = std::sqrt(y - t);
					 u = std::sqrt(y + static_cast<Real_t>(n - 1) * t);
					 khi = I4_POWER::i4_power(2, n);

					 for (std::size_t i{ 0 }; i != n; ++i) {
						 for (std::size_t j{ 0 }; j != n; ++j) {
							 x[j] = -v;
						 }
						  x[i] = -u;

						 for (std::size_t k{ 0 }; k != khi; ++k) {
							 ktemp = k;
							 for (std::size_t j{ 0 }; j != n; ++j) {
								 (ktemp % 2) == 2 ? x[j] = -std::fabs(x[j]) :
									                x[j] = std::fabs(x[j]);
								 ktemp /= 2;
							 }
							 quad += w * Fx(n, x);
						 }
						  x[i] = -v;
					 }

					 volume = BALL_UNIT_VOLUME_ND::ball_unit_volume_nd();
					 result = quad * volume;
					 return (result);

				 }
};

///////////////////////////////////////////////////////////////////////////////
//
//  Purpose:
//
//    BALL_UNIT_F3_ND approximates an integral inside the unit ball in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) <= 1.
//
//  Discussion:
//
//    A 2^(N+1)-1 point 5-th degree formula is used, Stroud number SN:5-4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  
//
//    
//
//  Author:
//
//    Adapted by Bernard Gingold from John Burkardt original implementation.
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the 
//    user supplied function which evaluates F at the N-vector X.
//
//    Input, int N, the dimension of the space.
//
//    Output, double BALL_UNIT_F3_ND, the approximate integral of the function.
///////////////////////////////////////////////////////////////////////////////

		template<class I4_POWER,
		         class BALL_UNIT_VOLUME_ND,
				 typename Real_t,
				 typename = std::enable_if<std::is_floating_point<Real_t>
				                              ::value>::type> struct Ball_unit_f3_nd {

					 
					 template<typename Real_t> using V = std::valarray<Real_t>;
					 
					 //std::unique_ptr implementation.
					 static Real_t ball_unit_f3_nd(_In_ Real_t(*Fx)(_In_ const int n,
						                           _In_ std::unique_ptr<Real_t, void(*)(Real_t*)> &x),
												   _In_ const int n) {


						 int jtemp{};
						 Real_t quad{},
							 result{},
							 ri{},
							 s{},
							 volume{},
							 weight{};

						 std::unique_ptr<Real_t, void(*)(Real_t*)> x{ new Real_t[n], [](Real_t *p)->void {delete[]p;} };
						 quad = CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;
						 for (int i{ 0 }; i != n; ++i) {
							 (&*x)[i] = CubConsts<Real_t, typename
								 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;
						 }
						 weight = CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_4_0_ / 
							 static_cast<Real_t>((n + 2) * (n + 2));
						 quad += weight * Fx(n,x);
						 s = CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ / 
							 std::sqrt(static_cast<Real_t>(n + 4));

						 for (int i{ 0 }; i != n; ++i) {
							 ri = std::sqrt(static_cast<Real_t>(i + 3) / static_cast<Real_t>(n + 4));
							 for (int j{ 0 }; j != n; ++j) {
								 j < i ? (&*x)[j] = CubConsts<Real_t, typename
									 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_ : 
									 (j == i ? (&*x)[j] = ri : (&*x)[j] = s );
									 
							 }
							 weight = std::pow(CubConsts<Real_t, typename
								 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_2_0_,
								 i + 1 - n) * static_cast<Real_t>(n + 4) / static_cast<Real_t>((i + 2) * 
								 (i + 3) * (n + 2));

							 for (int j{ 0 }; j != I4_POWER::i4_power(2, n - 1); ++j) {
								 jtemp = j;
								 for (int k{ 0 }; k != n; ++k) {
									 (jtemp % 2) == 1 ? (&*x)[k] = -std::fabs(&*x)[k] : 
										                (&*x)[k] = std::fabs((&*x)[k]);
									  jtemp /= 2;
								 }
								 quad += weight * Fx(n,x);
							 }
						 }

						 volume = BALL_UNIT_VOLUME_ND::ball_unit_volume_nd(n);
						 result = quad * volume;
						 return (result);
					 }

					 // valarray implementation.
					 static Real_t ball_unit_f3_nd(_In_ Real_t(*Fx)(_In_ const std::size_t n,
						                           _In_ const V &x),
						                           _In_ const std::size_t n) {

							
							std::size_t jtemp;
							Real_t quad{},
								result{},
								ri{},
								s{},
								volume{},
								weight{};

							V<Real_t> x;
							x.resize(n);
							quad = CubConsts<Real_t, typename
								std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;
							for (std::size_t i{ 0 }; i != n; ++i) {
								x[i] = CubConsts<Real_t, typename
									std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;
							}
							weight = CubConsts<Real_t, typename
								std::enable_if<std::is_floating_point<Real_t>::value>::type>::_4_0_ /
								static_cast<Real_t>((n + 2) * (n + 2));
							quad += weight * Fx(n, x);
							s = CubConsts<Real_t, typename
								std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ /
								std::sqrt(static_cast<Real_t>(n + 4));

							for (std::size_t i{ 0 }; i != n; ++i) {
								ri = std::sqrt(static_cast<Real_t>(i + 3) / static_cast<Real_t>(n + 4));
								for (std::size_t j{ 0 }; j != n; ++j) {
									j < i ? x[j] = CubConsts<Real_t, typename
										std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_ :
										(j == i ? x[j] = ri : x[j] = s);

								}
								weight = std::pow(CubConsts<Real_t, typename
									std::enable_if<std::is_floating_point<Real_t>::value>::type>::_2_0_,
									i + 1 - n) * static_cast<Real_t>(n + 4) / static_cast<Real_t>((i + 2) *
									(i + 3) * (n + 2));

								for (std::size_t j{ 0 }; j != I4_POWER::i4_power(2, n - 1); ++j) {
									jtemp = j;
									for (int k{ 0 }; k != n; ++k) {
										(jtemp % 2ULL) == 1ULL ? x[k] = -std::fabs(x[k]) :
											x[k] = std::fabs(x[k]);
										jtemp /= 2;
									}
									quad += weight * Fx(n, x);
								}
							}

							volume = BALL_UNIT_VOLUME_ND::ball_unit_volume_nd(n);
							result = quad * volume;
							return (result);
					 }
};

///////////////////////////////////////////////////////////////////////////////
//
//  Purpose:
//
//    BALL_UNIT_VOLUME_3D computes the volume of the unit ball in 3D.
//
//  Integration region:
//
//    X * X + Y * Y + Z * Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
// 
//
//   
//
//  Author:
//
//    Adapted by Bernard Gingold from John Burkardt's original implementation.
//
//  Parameters:
//
//    Output, double BALL_UNIT_VOLUME_3D, the volume of the ball.
///////////////////////////////////////////////////////////////////////////////


		template<typename Real_t,
			     typename = std::enable_if<std::is_floating_point<Real_t>
			                        ::value>::type> struct Ball_unit_volume_3d {

				
				static Real_t ball_unit_volume_3d() {

					return (CubConsts<Real_t, typename
						std::enable_if<std::is_floating_point<Real_t>::value>::type>::_4_div_3_ * 
						CubConsts<Real_t, typename
						std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_);

	    }
		       

};

///////////////////////////////////////////////////////////////////////////////
//
//  Purpose:
//
//    BALL_UNIT_VOLUME_ND computes the volume of the unit ball in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) <= 1.
//
//  Discussion:
//
//    N  Volume
//
//    2             PI
//    3  (4/3)    * PI
//    4  (1/2)    * PI^2
//    5  (8/15)   * PI^2
//    6  (1/6)    * PI^3
//    7  (16/105) * PI^3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  
//
//    
//
//  Author:
//
//    Adapted by Bernard Gingold from John Burkardt's original implementation.
//
//  Parameters:
//
//    Input, int N, the dimension of the space.
//
//    Output, double BALL_UNIT_VOLUME_ND, the volume of the ball.
///////////////////////////////////////////////////////////////////////////////


		template<class I4_POWER,
		         typename Real_t,
			     typename = std::enable_if<std::is_floating_point<Real_t>
			                             ::value>::type> struct Ball_unit_volume_nd {

				
				 static Real_t ball_unit_volume_nd(_In_ const int n) {

					 Real_t volume{};
					 int m{};

					 if ((n % 2) == 0) {
						 m = n / 2;
						 volume = std::pow(CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_, m);
						 for (int i{ 1 }; i <= m; ++i) volume /= static_cast<Real_t>(i);
					 }
					 else {
						 m = (n - 1) / 2;
						 volume = std::pow(CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_,m) * 
							 I4_POWER::i4_power(2, n);
						 for (int i{ m + 1 }; i <= 2 * m + 1; ++i) volume /= static_cast<Real_t>(i);
					 }
					 return (volume);
		}

				 static Real_t ball_unit_volume_nd(_In_ const std::size_t n) {
				      
					 Real_t volume{};
					 std::size_t m{};

					 if ((n % 2) == 0) {
						 m = n / 2;
						 volume = std::pow(CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_, m);
						 for (std:;size_t i{ 1 }; i <= m; ++i) volume /= static_cast<Real_t>(i);
					 }
					 else {
						 m = (n - 1) / 2;
						 volume = std::pow(CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_, m) *
							 I4_POWER::i4_power(2, n);
						 for (std::size_t i{ m + 1 }; i <= 2 * m + 1; ++i) volume /= static_cast<Real_t>(i);
					 }
					 return (volume);
				 }
};

///////////////////////////////////////////////////////////////////////////////
//
//  Purpose:
//
//    BALL_VOLUME_3D computes the volume of a ball in 3D.
//
//  Integration region:
//
//    X*X + Y*Y + Z*Z <= R * R
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  
//
//   
//
//  Author:
//
//   Adapted by Bernard Gingold from  John Burkardt's original implementation.
//
//  Parameters:
//
//    Input, double R, the radius of the ball.
//
//    Output, double BALL_VOLUME_3D, the volume of the ball.
///////////////////////////////////////////////////////////////////////////////

		template<typename Real_t,
			     typename = std::enable_if<std::is_floating_point<Real_t>
			                                  ::value>::type> struct Ball_volume_3d {

				 
				 static Real_t ball_volume_3d(_In_ const Real_t r) {
				        
					 return (CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_4_div_3_ * 
					         CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_ * r * r * r);
		 }
};

///////////////////////////////////////////////////////////////////////////////
//
//  Purpose:
//
//    BALL_VOLUME_ND computes the volume of a ball in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) <= R * R
//
//  Discussion:
//
//    N  Volume
//
//    2             PI   * R^2
//    3  (4/3)    * PI   * R^3
//    4  (1/2)    * PI^2 * R^4
//    5  (8/15)   * PI^2 * R^5
//    6  (1/6)    * PI^3 * R^6
//    7  (16/105) * PI^3 * R^7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  
//
//   
//
//  Author:
//
//   Adapted by Bernard Gingold from  John Burkardt's original implementation.
//
//  Parameters:
//
//    Input, int N, the dimension of the space.
//
//    Input, double R, the radius of the ball.
//
//    Output, double BALL_VOLUME_ND, the volume of the ball.
///////////////////////////////////////////////////////////////////////////////

		template<class BALL_UNIT_VOLUME_ND,
		         template Real_t,
				 template = std::enable_if<std::is_floating_point<Real_t>
				                              ::value>::type> struct Ball_volume_nd {

				
				  static Real_t ball_volume_nd(_In_ const int n,
											   _In_ const Real_t r) {
						 
					  return (BALL_UNIT_VOLUME_ND::ball_unit_volume_nd(n) * std::pow(r,n));
	 }

				  static Real_t ball_volume_nd(_In_ const std::size_t n,
					                           _In_ const Real_t r) {

					  return (BALL_UNIT_VOLUME_ND::ball_unit_volume_nd(n) * std::pow(r,n));
				  }
};

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_ANNULUS approximates an integral in an annulus.
//
//  Integration region:
//
//    RADIUS1^2 <= ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= RADIUS2^2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  
//
//    
//
//  Author:
//
//        Adapted by Bernard Gingold from  John Burkardt's original implementation.
//
//  Reference:
//
//    William Peirce,
//    Numerical Integration Over the Planar Annulus,
//    Journal of the Society for Industrial and Applied Mathematics,
//    Volume 5, Number 2, June 1957, pages 66-73.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y ), the name of the user supplied 
//    function of two variables which is to be integrated.
//
//    Input, double CENTER[2], the center of the circle.
//
//    Input, double RADIUS1, RADIUS2, the radii of the circles.
//
//    Input, int NR, the order of the rule.  This quantity specifies
//    the number of distinct radii to use.  The number of angles used will
//    be 4*NR, for a total of 4*NR*NR points.
//
//    Output, double CIRCLE_ANNULUS, the approximation to the integral.
///////////////////////////////////////////////////////////////////////////////

		template<class LEGENDRE_SET,
		         class RULE_ADJUST,
				 class CIRCLE_ANNULUS_AREA_2D,
				 typename Real_t,
				 typename = std::enable_if<std::is_floating_point<Real_t>
				                              ::value>::type> struct Circle_annulus {

					
				 static Real_t circle_annulus(_In_ Real_t(*Fx)(_In_ const Real_t x,
										      _In_ const Real_t y),
											  _In_ Real_t center[2],
											  _In_ Real_t const radius1,
											  _In_ Real_t const radius2,
											  _In_ int nr) {


					 Real_t a{},
						 area{},
						 b{},
						 c{},
						 d{},
						 quad{},
						 result{},
						 t{},
						 tw{},
						 x{},
						 y{};
					 int nt{};

					 std::unique_ptr<Real_t, void(*)(Real_t*)> ra{ new Real_t[nr], [](Real_t *p)->void {delete[]p; } };
					 std::unique_ptr<Real_t, void(*)(Real_t*)> rw{ new Real_t[nr], [](Real_t *p)->void {delete[]p;} };

					 LEGENDRE_SET::legendre_set(nr,ra,rw);
					 a = -CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_;
					 b = +CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_;
					 c = radius1 * radius1;
					 d = radius2 * radius2;
					 RULE_ADJUST::rule_adjust(a, b, c, d, nr, ra, rw);
					 // loop fusion here!
					 for (int i{ 0 }; i != nr; ++i) {
						 (&*ra)[i] = std::sqrt((&*ra)[i]);
						 (&*rw)[i] = (&*rw)[i] / (radius2 - radius1) / (radius2 + radius1);
					 }
					 //  Set angular abscissas and weights.
					 nt = 4 * nr;
					 tw = CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ / 
						 static_cast<Real_t>(nt);
					 
					 // Approximate the Integral.
					 quad = CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;

					 for (int i{ 0 }; i != nt; ++i) {
						 t = CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_2_0_ * 
							 CubConsts<Real_t, typename
							 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_ * 
							 static_cast<Real_t>(i - 1) / static_cast<Real_t>(nt);
						 for (int j{ 0 }; j != nr; ++j) {
							 x = center[0] + (&*ra)[j] * std::cos(t);
							 y = center[1] + (&*ra)[j] * std::sin(t);
							 quad += tw * (&*rw)[j] * Fx(x,y);
						 }
					 }
					 
					 area = CIRCLE_ANNULUS_AREA_2D::circle_annulus_area_2d(radius1,radius2);
					 result = quad * area;
					 return (result);

	  }
				


};

///////////////////////////////////////////////////////////////////////////////
//
//  Purpose:
//
//    CIRCLE_ANNULUS_AREA_2D returns the area of a circular annulus in 2D.
//
//  Integration region:
//
//    RADIUS1^2 <= ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= RADIUS2^2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  
//
//    
//
//  Author:
//
//   Adapted by Bernard Gingold from  John Burkardt's original implementation.
//
//  Parameters:
//
//    Input, double RADIUS1, RADIUS2, the radii of the circles.
//
//    Output, double CIRCLE_ANNULUS_AREA_2D, the area of the annulus.
///////////////////////////////////////////////////////////////////////////////

		template<typename Real_t,
			typename = std::enable_if<std::is_floating_point<Real_t>
			              ::value>::type> struct Circle_annulus_area_2d {

			
			 static Real_t circle_annulus_area_2d(_In_ const Real_t radius1,
			                                      _In_ const Real_t radius2) {

				 return (CubConsts<Real_t, typename
					 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_ * 
					 (radius1 + radius2) * (radius2 - radius1));
		  }
};

///////////////////////////////////////////////////////////////////////////////
//
//  Purpose:
//
//    CIRCLE_ANNULUS_SECTOR approximates an integral in a circular annulus sector.
//
//  Discussion:
//
//    A circular annulus sector comprises the area between two concentric
//    circles and two concentric rays.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  
//
//   
//
//  Author:
//
//   Adapted by Bernard Gingold from  John Burkardt's original implementation.
//
//  Reference:
//
//    William Peirce,
//    Numerical Integration Over the Planar Annulus,
//    Journal of the Society for Industrial and Applied Mathematics,
//    Volume 5, Number 2, June 1957, pages 66-73.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y ), the name of the 
//    user supplied function to be integrated.
//
//    Input, double CENTER[2], the center of the circle.
//
//    Input, double RADIUS1, RADIUS2, the radii of the circles.
//
//    Input, double THETA1, THETA2, the angles defining the sector.
//    The sector is measured from THETA1 to THETA2.
//
//    Input, int NR, the order of the rule.  This quantity specifies
//    the number of distinct radii to use.  The number of angles used will
//    be 4*NR, for a total of 4*NR*NR points.
//
//    Output, double CIRCLE_ANNULUS_SECTOR, the approximation to the integral.
///////////////////////////////////////////////////////////////////////////////


		template<class LEGENDRE_SET,
		         class RULE_ADJUST,
				 class TVEC_EVEN_BRACKET3,
				 class CIRCLE_ANNULUS_SECTOR_AREA_2D,
				 typename Real_t,
				 typename = std::enable_if<std::is_floating_point<Real_t>
				               ::value>::type> struct Circle_annulus_sector {


				static Real_t circle_annulus_sector(_In_ Real_t(*Fx)(_In_ const Real_t x,
								                                     _In_ const Real_t y),
								                    _In_ Real_t center[2],
								                    _In_ const Real_t radius1,
								                    _In_ const Real_t radius2,
								                    _In_ const Real_t theta1,
								                    _In_ const Real_t theta2,
													_In_ const int nr) {

					int nt{};
					Real_t a{},
						area{},
						b{},
						c{},
						d{},
						quad{},
						result{},
						x{},
						y{};

					std::unique_ptr<Real_t, void(*)(Real_t*)> ra{ new Real_t[nr], [](Real_t *p)->void {delete[]p; } },
						                                      rw{ new Real_t[nr], [](Real_t *p)->void {delete[]p; } },
						                                      ta{ new Real_t[nr], [](Real_t *p)->void {delete[]p; } },
						                                      tw{ new Real_t[nr], [](Real_t *p)->void {delete[]p;} };

					LEGENDRE_SET::legendre_set(nr,ra,rw);
					a = -CubConsts<Real_t, typename
						std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_;
					b = +CubConsts<Real_t, typename
						std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_;
					c = radius1 * radius1;
					d = radius2 * radius2;
					RULE_ADJUST::rule_adjust(a, b, c, d, nr, ra, rw);
					// Fusing for loops.
					// TODO:
					// Add #else branch for undefined FUSE_FOR_LOOP
#ifdef FUSE_FOR_LOOP
					for (int i{ 0 }; i != nr; ++i) {
						(&*ra)[i] = std::sqrt((&*ra)[i]);
						(&*rw)[i] = (&*rw)[i] / (radius2 - radius1) / (radius2 + radius1);
					}
#endif

					nt = 4 * nr;
					ta = TVEC_EVEN_BRACKET3::tvec_even_bracket3(nt,theta1,theta2);
					quad = CubConsts<Real_t, typename
						std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;

					// Fusing second loop
					// TODO:
					// Add #else branch for undefined FUSE_FOR_LOOP
#ifdef FUSE_FOR_LOOP
					for (int i{ 0 }; i != nt; ++i) {
						(&*tw)[i] =  CubConsts<Real_t, typename
							std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ / 
							static_cast<Real_t>(nt);
						for (int j{ 0 }; j != nr; ++j) {
							x = center[0] + (&*ra)[j] * std::cos((&*ta)[i]);
							y = center[1] + (&*ra)[j] * std::sin((&*ta)[i]);
							quad += (&*tw)[i] + (&*rw)[j] + Fx(x,y);
						}
						
					}
#endif			
					area = CIRCLE_ANNULUS_SECTOR_AREA_2D::circle_annulus_sector_area_2d(radius1,
					                                                                    radius2,
																						theta1,
																						theta2);
					result = quad * area;
					return (result;)

	  }


};

///////////////////////////////////////////////////////////////////////////////
//
//  Purpose:
//
//    CIRCLE_ANNULUS_SECTOR_AREA_2D returns the area of a circular annulus sector in 2D.
//
//  Discussion:
//
//    A circular annulus sector comprises the area between two concentric
//    circles and two concentric rays.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  
//
//  
//
//  Author:
//
//   Adapted by Bernard Gingold from  John Burkardt's original implementation.
//
//  Parameters:
//
//    Input, double RADIUS1, RADIUS2, the radii of the circles.
//
//    Input, double THETA1, THETA2, the angles of the rays.
//    Ordinarily, (THETA2-THETA1) is between 0 and 2*PI.
//
//    Output, double CIRCLE_ANNULUS_SECTOR_AREA_2D, the area of the
//    circular annulus sector.
///////////////////////////////////////////////////////////////////////////////


		template<typename Real_t,
			     typename = std::enable_if<std::is_floating_point<Real_t>
			                             ::value>::type> struct Circle_annulus_sector_area_2d {

				 
				 static Real_t circle_annulus_sector_area_2d(_In_ const Real_t radius1,
											                 _In_ const Real_t radius2,
											                 _In_ const Real_t theta1,
											                 _In_ const Real_t theta2) {

						
					 return (CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_div_2_ * 
						 (radius1 + radius2) * (radius2 - radius1) * (theta2 - theta1));
	 }
							
};

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_AREA_2D returns the area of a circle in 2D.
//
//  Integration region:
//
//    ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= R * R
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  
//
//    
//
//  Author:
//
//   Adapted by Bernard Gingold from  John Burkardt's original implementation.
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//
//    Output, double CIRCLE_AREA_2D, the area of the circle.
//


		template<typename Real_t,
			     typename = std::enable_if<std::is_floating_point<Real_t>
				                      ::value>::type> struct Circle_area_2d {
				
				 static Real_t circle_area_2d(_In_ const Real_t r) {
						
					 return (CubConsts<Real_t, typename
						 std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ * r * r);
	 }
};

///////////////////////////////////////////////////////////////////////////////
//
//  Purpose:
//
//    CIRCLE_CUM approximates an integral on the circumference of a circle in 2D.
//
//  Integration region:
//
//    ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= R * R
//
//  Discussion:
//
//    An ORDER point, (ORDER-1)-th degree formula is used, 
//    Stroud number U2:M-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  
//
//    
//
//  Author:
//
//   Adapted by Bernard Gingold from John Burkardt's original implementation.
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y ), the name of the 
//    user supplied function to be integrated.
//
//    Input, double CENTER[2], the coordinates of the center of 
//    the circle.
//
//    Input, double RADIUS, the radius of the circle.
//
//    Input, int ORDER, the number of points to use.
//
//    Output, double CIRCLE_CUM, the approximate integral of the function.
///////////////////////////////////////////////////////////////////////////////

		template<typename Real_t,
			     typename = std::enable_if<std::is_floating_point<Real_t>
				                      ::value>::type> struct Circle_circumference{


				static Real_t circle_circumference(_In_ Real_t(*Fx)(_In_ const Real_t x,
										           _In_ const Real_t y),
										           _In_ Real_t center[2],
										           _In_ Real_t const radius,
												   _In_ const int order) {



					Real_t angle{},
						quad{},
						result{},
						volume{},
						x{},
						y{},
						quad = CubConsts<Real_t, typename
						std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;

					for (int i{ 0 }; i != order; ++i) {
						angle = static_cast<Real_t>(2 * i) * CubConsts<Real_t, typename
							std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_ / 
							static_cast<Real_t>(order);
						x = center[0] + radius * std::cos(angle);
						y = center[1] + radius * std::sin(angle);
						quad += Fx(x,y);
					}
					quad /= static_cast<Real_t>(order);
					volume = CubConsts<Real_t, typename
						std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_ * 
						radius * radius;
					result = quad * volume;
					return (result);

			  }

};


///////////////////////////////////////////////////////////////////////////////
//
//  Purpose:
//
//    CIRCLE_RT_SET sets an R, THETA product quadrature rule in the unit circle.
//
//  Discussion:
//
//    For a given value of RULE, here are the number of points used at the
//    center (NC), the number of points along the radial direction (NR) and
//    the number of points along the theta direction (NT).  The total number
//    of points in the rule will be 
//
//      Total = NC + NR * NT.
//
//    The user, when choosing RULE, must allocate enough space in the arrays
//    RA, RW, TA and TW for the resulting values of NR and NT.
//
//    RULE  NC  NR  NT  Total
//    ----  --  --  --  -----
//       1   1   0   0      1
//       2   0   1   4      4
//       3   1   1   4      5
//       4   1   1   6      7
//       5   1   2   4      9
//       6   0   3   4     12
//       7   1   2  10     21
//       8   0   4  16     64
//       9   0   5  20    120
//
//    The integral of F(X,Y) over the unit circle is approximated by
//
//      Integral ( X*X + Y*Y <= 1 ) F(X,Y) dx dy 
//      = Integral ( 0 <= R <= 1, 0 <= T <= 2PI ) F(R*cos(T),R*sin(T)) r dr dt
//      = approximately
//        CW * F(0,0) 
//        + sum ( 1 <= I <= NR ) Sum ( 1 <= J <= NT )
//        RW(I) * TW(J) * F ( R(I) * cos ( TA(J) ), R(I) * sin ( TA(J) ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  
//
//    
//
//  Author:
//
//   Adapted by Bernard Gingold from  John Burkardt's original implementation.
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int RULE, the rule desired.
//
//    Input, int NR, the number of R abscissas.
//
//    Input, int NT, the number of Theta abscissas.
//
//    Input, int NC, the number of center abscissas (0 or 1 ).
//
//    Output, double RA[NR], RW[NR], the R abscissas and weights.
//
//    Output, double TA[NT], TW[NT], the THETA abscissas and weights.
//
//    Output, double *ZW, the weight to use for the center.
///////////////////////////////////////////////////////////////////////////////

		template<class LEGENDRE_SET,
		         class RULE_ADJUST,
				 typename Real_t,
				 typename = std::enable_if<std::is_floating_point<Real_t>
				                      ::value>::type> struct Circle_rt_set{


			    static void circle_rt_set(_In_ const int rule,
										  _In_ const int nr,
										  _In_ const int nt,
										  _In_ const int nc,
										  _Inout_ std::unique_ptr<Real_t, void(*)(Real_t*)> &ra,
										  _Inout_ std::unique_ptr<Real_t, void(*)(Real_t*)> &rw,
										  _In_    std::unique_ptr<Real_t, void(*)(Real_t*)> &ta,
										  _In_    std::unique_ptr<Real_t, void(*)(Real_t*)> &tw,
										  _Inout_ Real_t &cw) {


					Real_t a{},
						b{},
						c{},
						d{},
						u{},
						v{},
						w{};

					switch (rule) {

					case 1:
					{
						cw = CubConsts<Real_t, typename
							std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_;
					}
						break;
					case 2:
					{
						 (&*ra)[0] = CubConsts<Real_t, typename
							std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_div_2;
						 (&*rw)[0] = CubConsts<Real_t, typename
							std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_;
#ifdef FUSE_FOR_LOOP
						for (int i{ 0 }; i != nt; ++i) {
							(&*ta)[i] = static_cast<Real_t>(2 * i) * CubConsts<Real_t, typename
								std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_ / 
								static_cast<Real_t>(nt);
							(&*tw)[i] = CubConsts<Real_t, typename
								std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ / 
								static_cast<Real_t>(nt);
						}
#endif
						cw = CubConsts<Real_t, typename
							std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;
                    }
					    break;
					case 3:
					{
							  (&*ra)[0] = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_;
							  (&*rw)[0] = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_;
#ifdef FUSE_FOR_LOOP
							  for (int i{ 0 }; i != nt; ++i) {
								  (&*ta)[i] = static_cast<Real_t>(2 * i) * 
									  CubConsts<Real_t, typename
									  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_ / 
									  static_cast<Real_t>(nt);
								  (&*tw)[i] = CubConsts<Real_t, typename
									  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_div_8_;
							  }
#endif
							  cw = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_div_2_;
					}
					    break;
					case 4:
					{
							  (&*ra)[0] = std::sqrt(CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_2_div_3_);
							  (&*rw)[0] = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_;
#ifdef FUSE_FOR_LOOP
							  for (int i{ 0 }; i != nt; ++i) {
								  (&*ta)[i] = static_cast<Real_t>(2 * i) * CubConsts<Real_t, typename
									  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_ / 
									  static_cast<Real_t>(nt);
								  (&*tw)[0] = CubConsts<Real_t, typename
									  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_div_8_;
							  }
#endif						
							  cw = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_div_4_
					}
					   break;
					case 5:
					{
							  a = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_;
							  b = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_sqrt2_div_2_;
							  u = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_div_6_;
							  v = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_4_div_6_;

							  (&*ra)[0] = a;
							  (&*ra)[1] = b;
							  (&*rw)[0] = u;
							  (&*rw)[1] = v;
#ifdef FUSE_FOR_LOOP
							  for (int i{ 0 }; i != nt; ++i) {
								  (&*ta)[i] = static_cast<Real_t>(2 * i) * 
									  CubConsts<Real_t, typename
									  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_ / 
									  static_cast<Real_t>(nt);
								  (&*tw)[i] = CubConsts<Real_t, typename
									  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ / 
									  static_cast<Real_t>(nt);
							  }
#endif                         
							  cw = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_div_6_;
					}
					   break;
					case 6: 
				    {
							  a = std::sqrt(static_cast<Real_t>(3.0L)) * static_cast<Real_t>(0.5L);
							  b = std::sqrt((static_cast<Real_t>(27.0L - 3.0L * 5.385164807134504L)) * CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_div_52_);
							  c = std::sqrt((static_cast<Real_t>(27.0L + 3.0L * 5.385164807134504L)) * CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_div_52_);
							  v = static_cast<Real_t>(0.2962962962962962962962962962963L);
							  u = (static_cast<Real_t>(551.0L + 41.0L * 5.385164807134504)) * CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_div_1566_;
							  w = (static_cast<Real_t>(551.0L - 41.0L * 5.385164807134504)) * CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_div_1566_;
							  (&*ra)[0] = a;
							  (&*ra)[1] = b;
							  (&*ra)[2] = c;
							  (&*rw)[0] = u;
							  (&*rw)[1] = v;
							  (&*rw)[2] = w;
							  for (int i{ 0 }; i != nt; ++i) {
#ifdef FUSE_FOR_LOOP
								  (&*ta)[i] = static_cast<Real_t>(2 * i) *
									  CubConsts<Real_t, typename
									  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_ /
									  static_cast<Real_t>(nt);
								  (&*tw)[i] = CubConsts<Real_t, typename
									  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ /
									  static_cast<Real_t>(nt);
							  }
#endif
							  cw = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_;
				    }
					   break;
					case 7:
					{
							  a = static_cast<Real_t>(0.5958615826865181L);
							  b = static_cast<Real_t>(0.9192110607898046L);
							  u = static_cast<Real_t>(0.5124858261884216L);
							  v = static_cast<Real_t>(0.37640306270046727777777777777778L);

							  (&*ra)[0] = a;
							  (&*ra)[1] = b;
							  (&*rw)[0] = u;
							  (&*rw)[1] = v;
#ifdef FUSE_FOR_LOOP
							  for (int i{ 0 }; i != nt; ++i) {
								  (&*ta)[i] = static_cast<Real_t>(2 * i) *
									  CubConsts<Real_t, typename
									  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_ /
									  static_cast<Real_t>(nt);
								  (&*tw)[i] = CubConsts<Real_t, typename
									  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ /
									  static_cast<Real_t>(nt);
							  }
#endif
							  cw = static_cast<Real_t>(0.11111111111111111111111111111111L);
					}
					   break;
					case 8:
					{
							  LEGENDRE_SET::legendre_set(nr,ra,rw);
							  a = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_;
							  b = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_;
							  c = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;
							  d = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_;
							  RULE_ADJUST::rule_adjust(a,b,c,d,nr,na,rw);

							  for (int i{ 0 }; i != nr; ++i) {
								  (&*ra)[i] = std::sqrt((&*ra)[i]);
							  }
#ifdef FUSE_FOR_LOOP

							  for (int i{ 0 }; i != nt; ++i) {
								  (&*ta)[i] = static_cast<Real_t>(2 * i) *
									  CubConsts<Real_t, typename
									  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_ /
									  static_cast<Real_t>(nt);
								  (&*tw)[i] = CubConsts<Real_t, typename
									  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ /
									  static_cast<Real_t>(nt);
							  }

							  cw = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;
								 
#endif
					}
					  break;
					case 9:
					{
							  LEGENDRE_SET::legendre_set(nr,ra,rw);
							  a = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_;
							  b = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_;
							  c = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;
							  d = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_;
							  RULE_ADJUST::rule_adjust(a,b,c,d,nr,ra,rw);
							  for (int i{ 0 }; i != nr; ++i) {
								  (&*ra)[i] = std::sqrt((&*ra)[i]);
							  }
#ifdef FUSE_FOR_LOOP

							  for (int i{ 0 }; i != nt; ++i) {
								  (&*ta)[i] = static_cast<Real_t>(2 * i) *
									  CubConsts<Real_t, typename
									  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_pi_ /
									  static_cast<Real_t>(nt);
								  (&*tw)[i] = CubConsts<Real_t, typename
									  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_1_0_ /
									  static_cast<Real_t>(nt);
							  }

							  cw = CubConsts<Real_t, typename
								  std::enable_if<std::is_floating_point<Real_t>::value>::type>::_0_0_;

#endif
					}
					   break:
					default:
					{
							   std::cerr << "\n";
							   std::cerr << "Circle_rt_set::circle_rt_set - Fatal error!\n";
							   std::cerr << " No rule of index: " << rule << "\n";
							   std::exit(1);
					}

			}
	 }

};

	}
}


#endif /*__GMS_NDIM_CUBATURE_HPP__*/
