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
#ifndef __GMS_INDEXING_OPERATORS_HPP__
#define __GMS_INDEXING_OPERATORS_HPP__

//Code: Function definitions Col 8(Tab*2)
//Code: Function variables Col 32 (Tab*8)

// File version granularity.

namespace file_version {

    const unsigned int INDEXING_OPERATORS_MAJOR = 1U;
    const unsigned int INDEXING_OPERATORS_MINOR = 0U;
    const unsigned int INDEXING_OPERATORS_MICRO = 0U;
    const unsigned int INDEXING_OPERATORS_FULLVER =
      1000U*INDEXING_OPERATORS_MAJOR+
      100U*INDEXING_OPERATORS_MINOR+
      10U*INDEXING_OPERATORS_MICRO;
    const char * const INDEXING_OPERATORS_CREATION_DATE = "03-11-2023 09:06 AM +00200 (FRI 03 NOV 2023 GMT+2)";
    const char * const INDEXING_OPERATORS_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const INDEXING_OPERATORS_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const INDEXING_OPERATORS_DESCRIPTION   = "Multidimensional array subscriptor (indexing) operators."

}

#include <type_traits>

namespace gms {
	namespace math {



		enum class ArrayDimension : unsigned int{

		                DIM_1D = 1U,

				DIM_2D = 2U,

				DIM_3D = 3U,

				DIM_4D = 4U,

				DIM_5D = 5U,

				DIM_6D = 6U,

				DIM_7D = 7U,

				DIM_8D = 8U,

				DIM_9D = 9U,

				DIM_10D = 10U
		};


		
	    
	   template<typename Index_t,
		         typename = std::enable_if<std::is_integral<Index_t>::
				            value && !std::is_floating_point<Index_t>::value>::type> struct Subscriptor1D {

								Subscriptor1D() = default;

								__forceinline auto operator()(_In_ Index_t i)->Index_t { return (i); }

								inline constexpr Index_t dim(){ return static_cast<Index_t>(ArrayDimension::DIM_1D); }

				 };

		template<typename Index_t,
		          typename = std::enable_if<std::is_integral<Index_t>::
				             value && !std::is_floating_point<Index_t>::value>::type> struct Subscriptor2D {

								 Subscriptor2D() = default;

								__forceinline auto operator()(_In_ Index_t i,
									                          _In_ const Index_t isize,
									                          _In_ Index_t j)->Index_t {

									 return (i + isize * j);
								 }

							    inline constexpr Index_t dim() { return (static_cast<Index_t>(ArrayDimension::DIM_2D)); }
				  };

		template<typename Index_t,
		          typename = std::enable_if<std::is_integral<Index_t>::
				             value && !std::is_floating_point<Index_t>::value>::type> struct Subscriptor3D {

								 Subscriptor3D() = default;

								__forceinline auto operator()(_In_ Index_t i,
									             _In_ const Index_t isize,
									             _In_ Index_t j,
									             _In_ const Index_t jsize,
									             _In_ Index_t k)->Index_t {

									 return (i + isize * (j + jsize * k));
								 }

								inline constexpr Index_t dim() { return (static_cast<Index_t>(ArrayDimension::DIM_3D)); }
				  };

		template<typename Index_t,
		          typename = std::enable_if<std::is_integral<Index_t>::
				             value && !std::is_floating_point<Index_t>::value>::type> struct Subscriptor4D {

								 Subscriptor4D() = default;

								__forceinline auto operator()(_In_ Index_t i,
									                          _In_ const Index_t isize,
									                          _In_ Index_t j,
									                          _In_ const Index_t jsize,
									                          _In_ Index_t k,
									                          _In_ const Index_t ksize,
									                          _In_ Index_t l)->Index_t {

									 return (i + isize * (j + jsize * (k + ksize * l)));
								 }

								inline constexpr Index_t dim() { return (static_cast<Index_t>(ArrayDimension::DIM_4D)) };
				  };

		template<typename Index_t,
		          typename = std::enable_if<std::is_integral<Index_t>::
				             value && !std::is_floating_point<Index_t>::value>::type> struct Subscriptor5D {

								 Subscriptor5D() = default;

								 __forceinline auto operator()(_In_ Index_t i,
									                           _In_ const Index_t isize,
									                           _In_ Index_t j,
									                           _In_ const Index_t jsize,
									                           _In_ Index_t k,
									                           _In_ const Index_t ksize,
									                           _In_ Index_t l,
									                           _In_ const Index_t lsize,
									                           _In_ Index_t m)->Index_t {


									 return (i + isize * (j + jsize * (k + ksize * (l + lsize * m))));
								 }

								 inline constexpr Index_t dim() { return (static_cast<Index_t>(ArrayDimension::DIM_5D)) };
				  };

		template<typename Index_t,
		          typename = std::enable_if<std::is_integral<Index_t>::
				             value && !std::is_floating_point<Index_t>::value>::type> struct Subscriptor6D {

								 Subscriptor6D() = default;

								 __forceinline auto operator()(_In_ Index_t i,
									                           _In_ const Index_t isize,
									                           _In_ Index_t j,
									                           _In_ const Index_t jsize,
									                           _In_ Index_t k,
									                           _In_ const Index_t ksize,
									                           _In_ Index_t l,
									                           _In_ const Index_t lsize,
									                           _In_ Index_t m,
									                           _In_ const Index_t msize,
									                           _In_ Index_t n)->Index_t {

									 return (i + isize * (j + jsize * (k + ksize * (l + lsize * (m + msize * n)))));
								 }

								 inline constexpr Index_t dim() { return (static_cast<Index_t>(ArrayDimension::DIM_6D)) };
				  };

		template<typename Index_t,
		          typename = std::enable_if<std::is_integral<Index_t>::
				             value && !std::is_floating_point<Index_t>::value>::type> struct Subscriptor7D {

								 Subscriptor7D() = default;

								 __forceinline auto operator()(_In_ Index_t i,
									                           _In_ const Index_t isize,
									                           _In_ Index_t j,
									                           _In_ const Index_t jsize,
									                           _In_ Index_t k,
									                           _In_ const Index_t ksize,
									                           _In_ Index_t l,
									                           _In_ const Index_t lsize,
									                           _In_ Index_t m,
									                           _In_ const Index_t msize,
									                           _In_ Index_t n,
									                           _In_ const Index_t nsize,
									                           _In_ Index_t p)->Index_t {

									 return (i + isize * (j + jsize * (k + ksize * (l + lsize * 
									                                  (m + msize * (n + nsize * p))))));

                                   }
								   
								 inline constexpr Index_t dim() { return (static_cast<Index_t>(ArrayDimension::DIM_7D)); }
				  };

		template<typename Index_t,
		          typename = std::enable_if<std::is_integral<Index_t>::
				              value && !std::is_floating_point<Index_t>::value>::type> struct Subscriptor8D {

								  Subscriptor8D() = default;

								  __forceinline auto operator()(_In_ Index_t i,
									                            _In_ const Index_t isize,
									                            _In_ Index_t j,
									                            _In_ const Index_t jsize,
									                            _In_ Index_t k,
									                            _In_ const Index_t ksize,
									                            _In_ Index_t l,
									                            _In_ const Index_t lsize,
									                            _In_ Index_t m,
									                            _In_ const Index_t msize,
									                            _In_ Index_t n,
									                            _In_ const Index_t nsize,
									                            _In_ Index_t p,
									                            _In_ const Index_t psize,
									                            _In_ Index_t r)->Index_t {

									  return (i + isize * (j + jsize * (k + ksize * (l + lsize * 
									                      (m + msize * (n + nsize * (p + psize * r)))))));
								  }

								  inline constexpr Index_t dim() { return (static_cast<Index_t>(ArrayDimension::DIM_8D)); }
				  };

		template<typename Index_t,
		          typename = std::enable_if<std::is_integral<Index_t>::
				               value && !std::is_floating_point<Index_t>::value>::type> struct Subscriptor9D {

								   Subscriptor9D() = default;

								   __forceinline auto operator()(_In_ Index_t i,
									                             _In_ const Index_t isize,
									                             _In_ Index_t j,
									                             _In_ const Index_t jsize,
									                             _In_ Index_t k,
									                             _In_ const Index_t ksize,
									                             _In_ Index_t l,
									                             _In_ const Index_t lsize,
									                             _In_ Index_t m,
									                             _In_ const Index_t msize,
									                             _In_ Index_t n,
									                             _In_ const Index_t nsize,
									                             _In_ Index_t p,
									                             _In_ const psize,
									                             _In_ Index_t r,
									                             _In_ const Index_t rsize,
									                             _In_ Index_t q)->Index_t {

									   return (i + isize * (j + jsize * (k + ksize * (l + lsize * (m + msize * 
									                                     (n + nsize * (p + psize * (r + rsize * q))))))));
								   }

								   inline constexpr Index_t dim() { return (static_cast<Index_t>(ArrayDimension::DIM_9D)); }
				 };

		template<typename Index_t,
		         typename = std::enable_if<std::is_integral<Index_t>::
				             value && !std::is_floating_point<Index_t>::value>::type> struct Subscriptor10D{

								 Subscriptor10D() = default;

								 __forceinline auto operator()(_In_ Index_t i,
									                           _In_ const Index_t isize,
									                           _In_ Index_t j,
									                           _In_ const Index_t jsize,
									                           _In_ Index_t k,
									                           _In_ const Index_t ksize,
									                           _In_ Index_t l,
									                           _In_ const Index_t lsize,
									                           _In_ Index_t m,
									                           _In_ const Index_t msize,
									                           _In_ Index_t n,
									                           _In_ const Index_t nsize,
									                           _In_ Index_t p,
									                           _In_ const Index_t psize,
									                           _In_ Index_t r,
									                           _In_ const Index_t rsize,
									                           _In_ Index_t q,
									                           _In_ const Index_t qsize,
									                           _In_ Index_t s)->Index_t {

									 return (i + isize * (j + jsize * (k + ksize * (l + lsize * (m + msize *
									                     (n + nsize * (p + psize * (r + rsize * (q + qsize + s))))))));
								 }

								 inline constexpr Index_t dim() { return (static_cast<Index_t>(ArrayDimension::DIM_10D)); }
				 };

		

	}
}



#endif /*__GMS_INDEXING_OPERATORS_HPP__*/
