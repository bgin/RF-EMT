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
#ifndef __GMS_COMMON_H__
#define __GMS_COMMON_H__

namespace file_info {

	
        const unsigned int GMS_COMMON_MAJOR = 1;

	const unsigned int GMS_COMMON_MINOR = 1;

	const unsigned int GMS_COMMON_MICRO = 0;

	const unsigned int GMS_COMMON_FULLVER = 
	 1000U*GMS_COMMON_MAJOR+100U*GMS_COMMON_MINOR+10U*GMS_COMMON_MICRO;

	const char * const GMS_COMMON_CREATE_DATE = "27-09-2019 20:05 +00200 (FRI 27 SEP 2019 GMT+2)";

	const char * const GMS_COMMON_BUILD_DATE = __DATE__ ":" __TIME__;

	const char * const GMS_COMMON_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const GMS_COMMON_SYNOPSIS = "Common inlined functions.";
}





#include <cstdint>
#include "GMS_config.h"



namespace gms{
	

		

		// D. Knuth safe floating-point comparison

                static
                inline
		bool approximately_equalf64( const double a,
					     const double b,
					     const double eps) {
		     return std::fabs(a - b) <= ( (std::fabs(a) < std::fabs(b) ? 
		                        std::fabs(b) : std::fabs(a)) * eps);		 
		}
 
                static
                inline
		bool essentialy_equalf64( const double a,
					  const double b,
					  const double eps) {
		     return std::fabs(a - b) <= ((std::fabs(a) > std::fabs(b) ? 
		                        std::fabs(b) : std::fabs(a)) * eps);			  
		}

                static
                inline
		bool definitely_greaterf64( const double a,
					    const double b,
					    const double eps) {
		    return std::fabs(a - b) > ((std::fabs(a) < std::fabs(b) ? 
		                        std::fabs(b) : std::fabs(a)) * eps);			    
		}

                static
                inline
		bool definitely_lessf64( const double a,
					 const double b,
					 const double eps) {
			return std::fabs(b - a) > ((std::fabs(a) < std::fabs(b) ? 
			                std::fabs(b) : std::fabs(a)) * eps);		 
		}

                static
                inline
		bool approximately_equalf32(const float a,
					    const float b,
					    const float eps) {
		      	return std::fabsf(a - b) <= ((std::fabsf(a) > std::fabsf(b) ? 
		      	                  std::fabsf(b) : std::fabsf(a)) * eps);		    
                }

                static
                inline
		bool essentialy_equalf32(const float a,
					 const float b,
					 const float eps) {
			return std::fabsf(a - b) <= ((std::fabsf(a) > std::fabsf(b) ? 
			                   std::fabsf(b) : std::fabsf(a)) * eps);		 
                }

                static
                inline
		bool definitely_greaterf32(const float a,
					   const float b,
					   const float eps) {
		     return std::fabsf(a - b) > ((std::fabsf(a) < std::fabsf(b) ?
		                      std::fabsf(b) : std::fabsf(a)) * eps);		   
		}

                static
                inline
		bool definitely_lessf32(const float a,
					const float b,
					const float eps) {
		      return std::fabsf(b - a) > ((std::fabsf(a) < std::fabsf(b) ? 
			                std::fabsf(b) : std::fabsf(a)) * eps);			
	       }

		

		// Pointer alignment
		static
                inline
		bool
                Is_ptr_aligned32(const double * __restrict x) {
	               if ((reinterpret_cast<uintptr_t>(x)& 0x1F) == 0ULL) {
		            return (true);
	               }
	               else { return (false); }
                }

               static
               inline
               bool
               Is_ptr_aligned32(const int64_t * __restrict x) {
	              if ((reinterpret_cast<uintptr_t>(x) & 0x1F) == 0ULL) {
		          return (true);
	              }
	              else { return (false); }
               }

              static
              inline
              bool
              Is_ptr_aligned64(const double * __restrict x) {
	              if ((reinterpret_cast<uintptr_t>(x)& 0x3F) == 0ULL) {
		           return (true);
	              }
	              else { return (false); }
              }

             static
             inline
             bool
             Is_ptr_aligned64(const int64_t * __restrict x) {
	             if ((reinterpret_cast<uintptr_t>(x)& 0x3F) == 0ULL) {
		          return (true);
	             }
	             else { return (false); }
              }

		template<typename PTR, uint32_t Alignment,
			     typename = std::enable_if<(std::is_pointer<PTR>::value &&
							std::is_floating_point<PTR>::value) ||
							(std::is_pointer<PTR>::value &&
							std::is_integral<PTR>::value),bool>::type>
							Is_ptr_aligned(PTR * ptr) {
			       if ((reinterpret_cast<uintptr_t>(ptr) & Alignment) == 0ULL){
							  return (true);
					}
						 else {
							 return (false);
						 }
			}

}

#endif /*__GMS_COMMON_H__*/
