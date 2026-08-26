
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

#ifndef __GMS_AVX512C8F64_H__
#define __GMS_AVX512C8F64_H__

namespace file_info {


        const unsigned int gGMS_AVX512C8F64_MAJOR = 1;

	const unsigned int gGMS_AVX512C8F64_MINOR = 1;

	const unsigned int gGMS_AVX512C8F64_MICRO = 0;

	const unsigned int gGMS_AVX512C8F64_FULLVER = 
		1000U*gGMS_AVX512C8F64_MAJOR+100U*gGMS_AVX3C8F64_MINOR+10U*gGMS_AVX3C8F64_MICRO;

	const char * const pgGMS_AVX512C8F64_CREATE_DATE = "29-09-2018 11:13 +00200 (SAT 29 SEP 2018 GMT+2)";

	const char * const pgGMS_AVX512C8F64_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const pgGMS_AVX512C8F64_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_AVX512C8F64_SYNOPSIS = "AVX512 complex number class decomposed into real and imaginary parts stored as 8-tuple.";
}

#include <cstdint>
#include <iostream>
#include <immintrin.h>
#include <complex> // for comparison operators


#include "GMS_config.h"


namespace gms {
	namespace math {



	        struct 	  __attribute__((aligned(64))) ZMM8c8Payload {

			double re0,re1,re2,re3,re4,re5,re6,re7;
			double im0,im1,im2,im3,im4,im5,im6,im7;
		};

	


		//static AVX3C8f64 CZERO{ _mm512_setzero_pd(),
		//						_mm512_setzero_pd() };


	        struct 	  __attribute__((aligned(64))) ZMM8c8 {

				__m512d m_re;
				
				__m512d m_im;

				static const ZMM8c8 CZERO;
				
				// Default Ctor (sets components to 0.0)
				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
		                 __ATTR_VECTORCALL__
				ZMM8c8();
					
				 __ATTR_ALWAYS_INLINE__
                                 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
		                 __ATTR_VECTORCALL__
				 ZMM8c8(const ZMM8c8Payload);
					
				
				// Arrays: Re and Im must be aligned on 64-bytes boundary, otherwise GP will be signalled.
				  __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
		                 __ATTR_VECTORCALL__
				 ZMM8c8(const double * __restrict ,
					 const double * __restrict);
					
				
				// From single complex number
				  __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 ZMM8c8(const double ,
					 const double ); 
					
				
				// From scalar
				  __ATTR_ALWAYS_INLINE__
				  __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				explicit  ZMM8c8(const double ); 
					
				

				// Real parts only
				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 __ATTR_VECTORCALL__  
				 ZMM8c8(   const double ,
				            const double ,
				            const double ,
					    const double,
				            const double ,
					    const double ,
					    const double ,
					    const double );
					
				
                                 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 __ATTR_VECTORCALL__ 
				 ZMM8c8(const __m512d,
					 const __m512d);

		               
				
                                 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 __ATTR_VECTORCALL__ 
				 ZMM8c8(const  ZMM8c8 &);
					
				
                                 
				~ZMM8c8() noexcept(true) = default;

				
				// Load-store functions
				// Load aligned
				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 __ATTR_VECTORCALL__ 
				 ZMM8c8 & load_a(const double * __restrict ,
						     const double * __restrict );
					
				
				// Load unaligned
				  __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 __ATTR_VECTORCALL__ 
				 ZMM8c8 & load_u(const double * __restrict ,
						     const double * __restrict );
					
				

				// Store aligned
				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 __ATTR_VECTORCALL__ 
				void store_a(double * __restrict ,
					     double * __restrict ) const;
					
				

				// Store unaligned
				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 __ATTR_VECTORCALL__ 
				void store_u(double * __restrict ,
				             double * __restrict ) const; 
					
				
                                 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 __ATTR_VECTORCALL__ 
				void stream(double * __restrict ,
					    double * __restrict ) const;
					
                                  __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				double extract_f64(const int32_t) const;
					
                                 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				std::pair<double, double> 
				extract_2f64(const int32_t ,
					     const int32_t );
					
				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 ZMM8c8 & insert_1f64(const int32_t ,
						          const double );
					
                                  __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 ZMM8c8 & insert_2f64(const int32_t ,
						          const int32_t ,
						          const double ,
						          const double );
					

				// Length of 16 doubles
				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				  __ATTR_VECTORCALL__ 
				void concatenate_a(double * __restrict ) const;
					

				// Length of 16 doubles
				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				  __ATTR_VECTORCALL__ 
				void concatenate_u(double * __restrict ) const;

				   __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				  __ATTR_VECTORCALL__ 
				 ZMM8c8 & partial_loadu(const double const * __restrict,
							    const int32_t,
							    const double const * __restrict,
							    const int32_t );

				  __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				  __ATTR_VECTORCALL__ 			    
				 ZMM8c8 & partial_loada(const double const * __restrict,
							    const int32_t,
							    const double const * __restrict,
							    const int32_t);

				  __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				  __ATTR_VECTORCALL__ 			    
				void partial_storeu(double * __restrict,
					            const int32_t,
						    double * __restrict,
					            const int32_t);

				  __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				  __ATTR_VECTORCALL__ 		    
				void partial_storea(double * __restrict,
						    const int32_t,
						    double * __restrict,
					            const int32_t);
					
                                  __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				  __ATTR_VECTORCALL__ 
				 ZMM8c8 & expand(const  ZMM8c8,
						     const __mmask8); 

				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				  __ATTR_VECTORCALL__ 		     
				 ZMM8c8 & expand_load(const  ZMM8c8,
						          const __mmask8 ,
						          const double * __restrict ,
						          const double * __restrict );


                                 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 ZMM8c8 & permute(const __mmask8 ,
						      const int32_t);
					
                                 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				__m256d re_low2() const; 

				
                                 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				__m256d re_hi2() const; 

				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				__m256d im_low2() const;

				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				__m256d im_hi2() const;

				

				 ZMM8c8 & operator=(const  ZMM8c8);

				
		};

	


		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
		static inline  ZMM8c8
		        conj(const  ZMM8c8);

	         __ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 		
		static inline  ZMM8c8
		        polar(const __m512d,
			      const __m512d);

	        __ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 		      
		static inline __m512d
		        carg(const  ZMM8c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	
		static inline __m512d 
		        carg(const double,
			     const double);

		 __ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
		static inline  ZMM8c8 
			csin(const  ZMM8c8);

		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 		
		static inline  ZMM8c8
		        csin(const double,
			     const double);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__     
		static inline  ZMM8c8
		        csinh(const  ZMM8c8);

		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
		static inline  ZMM8c8
		        csinh(const double,
			      const double);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	      
		static inline  ZMM8c8
			ccos(const  ZMM8c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	
		static inline  ZMM8c8
		        ccos(const double,
			     const double);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	     
		static inline  ZMM8c8
		        ccosh(const  ZMM8c8);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	
		static inline  ZMM8c8
		        ccosh(const double,
			      const double);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	      
		static inline  ZMM8c8
			cexp(const  ZMM8c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		
		static inline  ZMM8c8
		        cexp(const double,
			     const double);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		     
		static inline __m512d
			cabs(const  ZMM8c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		
		static inline __m512d
		        cabs(const double,
			     const double);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		     
		static inline  ZMM8c8
		        cpowi(const ZMM8c8,
			      const double );

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		      
		static inline ZMM8c8
		        clog(const  ZMM8c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		
		static inline  ZMM8c8
		        clog(const double,
			     const double);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		     
		static inline  ZMM8c8
		        csqrt(const  ZMM8c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		
		static inline  ZMM8c8
		        csqrt(const double,
			      const double);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		      
		static inline  ZMM8c8
		        ctan(const  ZMM8c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		
		static inline  ZMM8c8
		        ctan(const double,
			     const double);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		     
		static inline  ZMM8c8
		        ctanh(const  ZMM8c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		
		static inline  ZMM8c8
		        ctanh(const double,
			      const double);
			      
			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	
		static inline  ZMM8c8
		select(const  ZMM8c8 ,
		       const  ZMM8c8,
		       const __mmask8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	       
		static inline  ZMM8c8
		cdiv_smith(const ZMM8c8 ,
			   const  ZMM8c8 );

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		   
		static inline  ZMM8c8
		operator+(const  ZMM8c8,
			  const  ZMM8c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		  
		static inline  ZMM8c8
		operator+(const  ZMM8c8,
			  const __m512d);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator+(const __m512d,
			  const  ZMM8c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator+(const  ZMM8c8,
		          const double );
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator+(const double ,
		          const  ZMM8c8); 
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator+=( ZMM8c8,
		           const  ZMM8c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator+=( ZMM8c8,
			   const __m512d);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator+=(const __m512d,
			    ZMM8c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator+=( ZMM8c8 ,
		          const double);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator+=(const double,
			    ZMM8c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator-(const  ZMM8c8,
			  const  ZMM8c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator-(const  ZMM8c8,
			  const __m512d);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline ZMM8c8
		operator-(const __m512d,
			  const  ZMM8c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator-(const  ZMM8c8,
			  const double );
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator-(const double ,
			  const ZMM8c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator-( ZMM8c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator-=( ZMM8c8,
			   const  ZMM8c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator-=( ZMM8c8,
			   const __m512d);
		  
		static inline  ZMM8c8
		operator-=(const __m512d,
			    ZMM8c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator-=( ZMM8c8,
			   const double );
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator-=(const double,
		            ZMM8c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline ZMM8c8
		operator*(const  ZMM8c8,
			  const  ZMM8c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator*(const  ZMM8c8,
			  const __m512d);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator*(const __m512d,
			  const  ZMM8c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator*(const  ZMM8c8,
		          const double);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		  
		static inline  ZMM8c8
		operator*(const double,
		          const  ZMM8c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		  
		static inline  ZMM8c8
		operator*=( ZMM8c8,
			   const  ZMM8c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		   
		static inline  ZMM8c8
		operator*=( ZMM8c8,
			   const __m512d);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		   
		static inline  ZMM8c8
		operator*=(const __m512d,
			    ZMM8c8);

		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__			    
	        static inline  ZMM8c8
		operator*=( ZMM8c8,
			   const double);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator*=(const double,
		           ZMM8c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator/(const  ZMM8c8,
		          const  ZMM8c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator/(const ZMM8c8,
			  const __m512d);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator/(const __m512d,
			  const  ZMM8c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM8c8
		operator/(const  ZMM8c8,
			  const double);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline ZMM8c8
		operator/(const double,
			  const  ZMM8c8);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	  
		static inline  ZMM8c8
		operator/=( ZMM8c8,
			   const  ZMM8c8);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	   
		static inline  ZMM8c8
		operator/=( ZMM8c8,
			   const __m512d);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	   
		static inline  ZMM8c8
		operator/=(const __m512d,
			    ZMM8c8);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	    
	        static inline  ZMM8c8
		operator/=( ZMM8c8,
			   const double);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	   
		static inline  ZMM8c8
		operator/=(const double,
		            ZMM8c8);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	    
		static inline  ZMM8c8
		operator~( ZMM8c8);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline
			std::pair<__mmask8, __mmask8>
		operator==(const  ZMM8c8,
			   const  ZMM8c8);

		

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	   
		static inline
		        std::pair<__mmask8,__mmask8>
		operator==(const  ZMM8c8,
		           std::complex<double>);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	   
		static inline
		        std::pair<__mmask8,__mmask8>
		operator==(const std::complex<double>,
		           const  ZMM8c8)
		           

			   	__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
		static inline 
			std::pair<__mmask8, __mmask8>
		operator!=(const  ZMM8c8,
		           const  ZMM8c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	   
		static inline
		        std::pair<__mmask8,__mmask8>
		operator!=(const  ZMM8c8,
		           const std::complex<double>);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	   
		static inline
		        std::pair<__mmask8,__mmask8>
		operator!=(const std::complex<double>,
		           const  ZMM8c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	   
		static inline 
			std::pair<__mmask8, __mmask8>
		operator>(const  ZMM8c8,
			  const  ZMM8c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	  
		static inline
                        std::pair<__mmask8,__mmask8>
		operator>(const  ZMM8c8,
		          const std::complex<double>);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	  
		static inline	  
			std::pair<__mmask8,__mmask8>
		operator>(const std::complex<double>,
		          const  ZMM8c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	  
		static inline 
			std::pair<__mmask8, __mmask8>
		operator<(const  ZMM8c8,
			  const  ZMM8c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	  
		static inline	  
                        std::pair<__mmask8,__mmask8>
		operator<(const  ZMM8c8,
		          std::complex<double>);

			  	__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
                static inline
		        std::pair<__mmask8,__mmask8>
		operator<(const std::complex<double>,
		          const  ZMM8c8);

			  	__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
		static inline 
			std::pair<__mmask8, __mmask8>
		operator>=(const  ZMM8c8,
		           const  ZMM8c8);

			   	__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
		static inline
		        std::pair<__mmask8,__mmask8>
		operator>=(const  ZMM8c8,
		           const std::complex<double>);

			   	__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
		static inline
		        std::pair<__mmask8,__mmask8>
		operator>=(const std::complex<double>,
		           const  ZMM8c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	   
		static inline 
			std::pair<__mmask8, __mmask8>
	        operator<=(const  ZMM8c8,
			   const  ZMM8c8);

			   	__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
		static inline
		        std::pair<__mmask8,__mmask8>
		operator<=(const  ZMM8c8,
		           const std::complex<double>);
			   
                __ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
		static inline
		        std::pair<__mmask8,__mmask8>
		operator<=(const std::complex<double>,
		           const  ZMM8c8);

#include "GMS_avx512c8f64.inl"	

	}
}



#endif /*__GMS_AVX512C8F64_H__*/
