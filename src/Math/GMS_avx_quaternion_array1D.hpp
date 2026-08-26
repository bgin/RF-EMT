
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
#ifndef __GMS_AVX_QUATERNION_ARRAY1D_HPP__
#define __GMS_AVX_QUATERNION_ARRAY1D_HPP__


namespace file_info {

const unsigned int gGMS_AVX_QUATERNION_ARRAY1D_MAJOR = 1U;
const unsigned int gGMS_AVX_QUATERNION_ARRAY1D_MINOR = 1U;
const unsigned int gGMS_AVX_QUATERNION_ARRAY1D_MICRO = 0U;
const unsigned int gGMS_AVX_QUATERNION_ARRAY1D_FULLVER =
   1000U*gGMS_AVX_QUATERNION_ARRAY1D_MAJOR+
   100U*gGMS_AVX_QUATERNION_ARRAY1D_MINOR+
   10U*gGMS_AVX_QUATERNION_ARRAY1D_MICRO;

const char * const pgGMS_AVX_QUATERNION_ARRAY1D_CREATION_DATE = "01-08-2021 09:58 +00200 (SUN 01 AUG 2021 GMT+2)";
const char * const pgGMS_AVX_QUATERNION_ARRAY1D_BUILD_DATE = __DATE__ ":" __TIME__;
const char * const pgGMS_AVX_QUATERNION_ARRAY1D_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
const char * const pgGMS_AVX_QUATERNION_ARRAY1D_DESCRIPTION = "Single-dimensional array of quaternionic vectors."
}


#include <cstdint>
#include <omp.h>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_malloc.h"
#include "GMS_common.h"
#include "GMS_constants.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined(USE_AVXQUATERNION_ARRAY1D_NT_STORES)
    #define USE_AVXQUATERNION_ARRAY1D_NT_STORES 0
#endif

// Enable/disable software prefetching.
#if !defined(AVXQUATERNION_ARRAY1D_PREFETCH_FROM_OBJ)
#define AVXQUATERNION_ARRAY1D_PREFETCH_FROM_OBJ(obj,idx,off,hint)  \
   _mm_prefetch(reinterpret_cast<const char*>(&(obj).m_x[idx+off],hint);
   _mm_prefetch(reinterpret_cast<const char*>(&(obj).m_y[idx+off],hint);
   _mm_prefetch(reinterpret_cast<const char*>(&(obj).m_z[idx+off],hint);
   _mm_prefetch(reinterpret_cast<const char*>(&(obj).m_w[idx+off],hint);
#endif

#if !defined(AVXQUATERNION_ARRAY1D_PREFETCH_FROM_PTR)
#define AVXQUATERNION_ARRAY1D_PREFETCH_FROM_PTR(ptr,idx,off,hint)  \
   _mm_prefetch(reinterpret_cast<const char*>(&(ptr)[idx+off],hint);
#endif


namespace gms {

          namespace math {


              struct AVXQuatArray1D __ATTR_ALIGN__(32) {

                     __ATTR_ALIGN__(8) double * __restrict m_x; // x - part
		     __ATTR_ALIGN__(8) double * __restrict m_y; // y - part
		     __ATTR_ALIGN__(8) double * __restrict m_z; // z - part
		     __ATTR_ALIGN__(8) double * __restrict m_w; // w - part
		     int32_t m_size; // size of quaternion array 1D.
#if(USE_STRUCT_PADDING) == 1
                     PAD_TO_ALIGNED(8,0,28)
#endif
                     /*
                          Default Ctor
                      */
		     __ATTR_ALWAYS_INLINE__
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
                     AVXQuatArray1D() {

                         m_x = nullptr;
			 m_y = nullptr;
			 m_z = nullptr;
			 m_w = nullptr;
			 m_size = 0;
		     }

		     __ATTR_ALWAYS_INLINE__
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     AVXQuatArray1D(const int32_t size,
		                    const bool lockm) {
                        using namespace gms::common;
			m_size = size;
			const std::size_t ullsize = static_cast<m_size>;
#if (USE_MMAP_2MiB) == 1
#pragma omp parallel sections
{
                #pragma omp section
		{
		  m_x = (double*)gms_mmap_2MiB(ullsize,
			                      PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE,
					       -1,0);
		}
		#pragma omp section
		{
		  m_y = (double*)gms_mmap_2MiB(ullsize,
			                      PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE,
					       -1,0);
		}
		#pragma omp section
		{
		  m_z = (double*)gms_mmap_2MiB(ullsize,
			                      PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE,
					       -1,0);
		}
		#pragma omp section
		{
		  m_w = (double*)gms_mmap_2MiB(ullsize,
			                      PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE,
					       -1,0);
		}
}
#else
#pragma omp parallel sections
{
                #pragma omp section
		{
		  m_x = (double*)gms_mm_malloc(ullsize,align32B);
		}
		#pragma omp section
		{
		  m_y = (double*)gms_mm_malloc(ullsize,align32B);
		}
		#pragma omp section
		m_z = (double*)gms_mm_malloc(ullsize,align32B);
		}
		#pragma omp section
		{
		  m_w = (double*)gms_mm_malloc(ullsize,align32B);
		}
}
#endif
                        // Zero initialization.
#if (ZERO_INIT_ARRAYS) == 1

             	      avx256_init_unroll8x_pd(&m_x[0],m_size,0.0);
	      	      avx256_init_unroll8x_pd(&m_y[0],m_size,0.0);
	      	      avx256_init_unroll8x_pd(&m_z[0],m_size,0.0);
	      	      avx256_init_unroll8x_pd(&m_w[0],m_size,0.0);

#endif
		 }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     AVXQuatArray1D(const double * __restrict __ATTR_ALIGN__(32) x,
		                    const double * __restrict __ATTR_ALIGN__(32) y,
				    const double * __restrict __ATTR_ALIGN__(32) z,
				    const double * __restrict __ATTR_ALIGN__(32) w,
				    const int32_t size,
				    const bool lockm) {
                          using namespace gms::common;
			  m_size = size;
                          const std::size_t ullsize = static_cast<m_size>;
#if (USE_MMAP_2MiB) == 1
#pragma omp parallel sections
{
                    #pragma omp section
		    {
		      m_x = (double*)gms_mmap_2MiB(ullsize,
			                    PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE,
					    -1,0);
		    }
		    #pragma omp section
		    {
		      m_y = (double*)gms_mmap_2MiB(ullsize,
			                    PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE,
					    -1,0);
		    }
		    #pragma omp section
		    {
		      m_z = (double*)gms_mmap_2MiB(ullsize,
			                    PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE,
					    -1,0);
		    }
		    #pragma omp section
		    {
		      m_w = (double*)gms_mmap_2MiB(ullsize,
			                    PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE,
					    -1,0);
		    }
}
#else
#pragma omp parallel sections
{
                  #pragma omp section
		  {
		    m_x = (double*)gms_mm_malloc(ullsize,align32B);
		  }
		  #pragma omp section
		  {
		    m_y = (double*)gms_mm_malloc(ullsize,align32B);
		  }
		  #pragma omp section
		  {
		    m_z = (double*)gms_mm_malloc(ullsize,align32B);
		  }
		  #pragma omp section
		  {
		    m_w = (double*)gms_mm_malloc(ullsize,align32B);
		  }
}
#endif
#if (USE_AVXQUATERNION_ARRAY1D_NT_STORES) == 1
                        avx256_uncached_memmove(&m_x[0],&x[0],m_size);
			avx256_uncached_memmove(&m_y[0],&y[0],m_size);
			avx256_uncached_memmove(&m_z[0],&z[0],m_size);
			avx256_uncached_memmove(&m_w[0],&w[0],m_size);
#else
                        avx256_cached_memmove(&m_x[0],&x[0],m_size);
			avx256_cached_memmove(&m_y[0],&y[0],m_size);
			avx256_cached_memmove(&m_z[0],&z[0],m_size);
			avx256_cached_memmove(&m_w[0],&w[0],m_size);
#endif
		    }


		    AVXQuatArray1D(const AVXQuatArray1D &x) = delete;

		    __ATTR_ALWAYS_INLINE__
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    AVXQuatArray1D(AVXQuatArray1D &&rhs) {

                        m_x = &rhs.m_x[0];
			m_y = &rhs.m_y[0];
			m_z = &rhs.m_z[0];
			m_w = &rhs.m_w[0];
			m_size = rhs.m_size;
			rhs.m_x = NULL;
			rhs.m_y = NULL;
			rhs.m_z = NULL;
			rhs.m_w = NULL;
			rhs.m_size = 0;
		    }


		    __ATTR_ALWAYS_INLINE__
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    ~AVXQuatArray1D() {
                       if(NULL != m_x) _mm_free(m_x); m_x = NULL;
		       if(NULL != m_y) _mm_free(m_y); m_y = NULL;
		       if(NULL != m_z) _mm_free(m_z); m_z = NULL;
		       if(NULL != m_w) _mm_free(m_w); m_w = NULL;
                       
		    }

		    AVXQuatArray1D & operator=(const AVXQuatArray1D &rhs) = delete;

		    __ATTR_ALWAYS_INLINE__
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    AVXQuatArray1D &
		    operator=(AVXQuatArray1D &&rhs) {
                        if(this == &rhs) return (*this);
			_mm_free(m_x);
			_mm_free(m_y);
			_mm_free(m_z);
			_mm_free(m_w);
			m_x = &rhs.m_x[0];
			rhs.m_x = NULL;
			m_y = &rhs.m_y[0];
			rhs.m_y = NULL;
			m_z = &rhs.m_z[0];
			rhs.m_z = NULL;
			m_w = &rhs.m_w[0];
			rhs.m_w = NULL;
			m_size = rhs.m_size;
			rhs.m_size = 0;
			return (*this);
		    }

		     
		      
	      }; // AVXQuatArray1D


	       __ATTR_ALWAYS_INLINE__
	       __ATTR_HOT__
	       __ATTR_ALIGN__(32)
	       static inline
	       void AVXQuat_Add_AVXQuat(AVXQuatArray1D &c,
	                                const AVXQuatArray1D &b,
					const AVXQuatArray1D &a) {
					
		      if(__builtin_expect(b.m_size != a.m_size,0) ||
		         __builtin_expect(c.m_size != a.m_size,0)) { return;}
                      int32_t i,last_i;
		      last_i = 0;
		      const int32_t size = c.m_size;
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma code_align(32)
#endif
#pragma omp parallel for schedule(static,32) default(none) private(i) \
                      lastprivate(last_i) shared(c,b,a,size) if(size>=5000)
		      for(i = 0; i != ROUND_TO_FOUR(size,4); i += 16) {

		          last_i = i;
			  _mm_prefetch((const char*)&b.m_x[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_x[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&b.m_y[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_y[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&b.m_z[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_z[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&b.m_w[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_w[i+16],_MM_HINT_T0);
			  _mm256_store_pd(&c.m_x[i+0],
			                  _mm256_add_pd(_mm256_load_pd(&b.m_x[i+0]),
					                  _mm256_load_pd(&a.m_x[i+0])));
			  _mm256_store_pd(&c.m_y[i+0],
			                  _mm256_add_pd(_mm256_load_pd(&b.m_y[i+0]),
					                 _mm256_load_pd(&a.m_y[i+0])));
			  _mm256_store_pd(&c.m_z[i+0],
			                  _mm256_add_pd(_mm256_load_pd(&b.m_z[i+0]),
					                 _mm256_load_pd(&a.m_z[i+0]) ));
			  _mm256_store_pd(&c.m_w[i+0],
			                  _mm256_add_pd(_mm256_load_pd(&b.m_w[i+0]),
					                _mm256_load_pd(&a.m_w[i+0])));
			  _mm256_store_pd(&c.m_x[i+4],
			                  _mm256_add_pd(_mm256_load_pd(&b.m_x[i+4]),
					                 _mm256_load_pd(&a.m_x[i+4])));
			  _mm256_store_pd(&c.m_y[i+4],
			                  _mm256_add_pd(_mm256_load_pd(&b.m_y[i+4]),
					                  _mm256_load_pd(&a.m_y[i+4])));
			  _mm256_store_pd(&c.m_z[i+4],
			                  _mm256_add_pd(_mm256_load_pd(&b.m_z[i+4]),
					                 _mm256_load_pd(&a.m_z[i+4])));
			  _mm256_store_pd(&c.m_w[i+4],
			                  _mm256_add_pd(_mm256_load_pd(&b.m_w[i+4]),
					                 _mm256_load_pd(&a.m_w[i+4])));
			  _mm256_store_pd(&c.m_x[i+8],
			                  _mm256_add_pd(_mm256_load_pd(&b.m_x[i+8]),
					                _mm256_load_pd(&a.m_x[i+8]) ));
			  _mm256_store_pd(&c.m_y[i+8],
			                  _mm256_add_pd(_mm256_load_pd(&b.m_y[i+8]),
					                _mm256_load_pd(&a.m_y[i+8])));
			  _mm256_store_pd(&c.m_z[i+8],
			                  _mm256_add_pd(_mm256_load_pd(&b.m_z[i+8]),
					                _mm256_load_pd(&a.m_z[i+8])));
			  _mm256_store_pd(&c.m_w[i+8],
			                  _mm256_add_pd(_mm256_load_pd(&b.m_w[i+8]),
					                 _mm256_load_pd(&a.m_w[i+8])));
			  _mm256_store_pd(&c.m_x[i+12],
			                  _mm256_add_pd(_mm256_load_pd(&b.m_x[i+12]),
					                _mm256_load_pd(&a.m_x[i+12])));
			  _mm256_store_pd(&c.m_y[i+12],
			                  _mm256_add_pd(_mm256_load_pd(&b.m_y[i+12]),
					                _mm256_load_pd(&a.m_y[i+12])));
			  _mm256_store_pd(&c.m_z[i+12],
			                  _mm256_add_pd(_mm256_load_pd(&b.m_z[i+12]),
					                 _mm256_load_pd(&a.m_z[i+12])));
			  _mm256_store_pd(&c.m_w[i+12],
			                  _mm256_add_pd(_mm256_load_pd(&b.m_w[i+12]),
					                _mm256_load_pd(&a.m_w[i+12])));
			  
		      }
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma loop_count min(1),avg(2),max(3)
#endif
                      for(; last_i != size; ++last_i) {
                           c.m_x[last_i] = b.m_x[last_i]+a.m_x[last_i];
			   c.m_y[last_i] = b.m_y[last_i]+a.m_y[last_i];
			   c.m_z[last_i] = b.m_z[last_i]+a.m_z[last_i];
			   c.m_w[last_i] = b.m_w[last_i]+a.m_w[last_i];
		      }
	       }


               __ATTR_ALWAYS_INLINE__
	       __ATTR_HOT__
	       __ATTR_ALIGN__(32)
	       static inline
	       void AVXQuat_Add_AVXQuat(AVXQuatArray1D &c,
	                                AVXQuatArray1D &b,
	                                const double * __restrict __ATTR_ALIGN__(32) re,
					const double * __restrict __ATTR_ALIGN__(32) im,
					const int32_t n) {

		     if(__builtin_expect(c.m_size != n,0) ||
		        __builtin_expect(c.m_size != b.m_size,0)) { return;}
                        int32_t i,last_i;
		        last_i = 0;
			const int32_t size = c.m_size;
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma code_align(32)
#endif
#pragma omp parallel for schedule(static,32) default(none) private(i) \
                      lastprivate(last_i) shared(c,b,re,im,size) if(size>=5000)
                     for(i = 0; i != ROUND_TO_FOUR(size,4); i += 16) {
                         last_i = i;
			 _mm_prefetch((const char*)&b.m_x[i+16],_MM_HINT_T0);
			 _mm_prefetch((const char*)&b.m_y[i+16],_MM_HINT_T0);
			 _mm_prefetch((const char*)&b.m_z[i+16],_MM_HINT_T0);
			 _mm_prefetch((const char*)&b.m_w[i+16],_MM_HINT_T0);
			 _mm_prefetch((const char*)&re[i+16],_MM_HINT_T0);
			 _mm_prefetch((const char*)&im[i+16],_MM_HINT_T0);
			 _mm256_store_pd(&c.m_x[i+0],
			                     _mm256_add_pd(_mm256_load_pd(&b.m_x[i+0]),
					                   _mm256_load_pd(&re[i+0])));
			 _mm256_store_pd(&c.m_y[i+0],
			                     _mm256_add_pd(_mm256_load_pd(&b.m_y[i+0]),
					                   _mm256_load_pd(&im[i+0])));
			 _mm256_store_pd(&c.m_z[i+0],
			                     _mm256_add_pd(_mm256_load_pd(&b.m_z[i+0]),
					                   _mm256_load_pd(&re[i+0])));
			 _mm256_store_pd(&c.m_w[i+0],
			                     _mm256_add_pd(_mm256_load_pd(&b.m_w[i+0]),
					                   _mm256_load_pd(&im[i+0])));
			 _mm256_store_pd(&c.m_x[i+4],
			                     _mm256_add_pd(_mm256_load_pd(&b.m_x[i+4]),
					                   _mm256_load_pd(&re[i+4])));
			 _mm256_store_pd(&c.m_y[i+4],
			                     _mm256_add_pd(_mm256_load_pd(&b.m_y[i+4]),
					                   _mm256_load_pd(&im[i+4])));
			 _mm256_store_pd(&c.m_z[i+4],
			                     _mm256_add_pd(_mm256_load_pd(&b.m_z[i+4]),
					                   _mm256_load_pd(&re[i+4])));
			 _mm256_store_pd(&c.m_w[i+4],
			                     _mm256_add_pd(_mm256_load_pd(&b.m_w[i+4]),
					                   _mm256_load_pd(&im[i+4])));
			 _mm256_store_pd(&c.m_x[i+8],
			                     _mm256_add_pd(_mm256_load_pd(&b.m_x[i+8]),
					                   _mm256_load_pd(&re[i+8])));
			 _mm256_store_pd(&c.m_y[i+8],
			                     _mm256_add_pd(_mm256_load_pd(&b.m_y[i+8]),
					                   _mm256_load_pd(&im[i+8])));
			 _mm256_store_pd(&c.m_z[i+8],
			                     _mm256_add_pd(_mm256_load_pd(&b.m_z[i+8]),
					                   _mm256_load_pd(&re[i+8])));
			 _mm256_store_pd(&c.m_w[i+8],
			                     _mm256_add_pd(_mm256_load_pd(&b.m_w[i+8]),
					                   _mm256_load_pd(&im[i+8])));
			 _mm256_store_pd(&c.m_x[i+12],
			                     _mm256_add_pd(_mm256_load_pd(&b.m_x[i+12]),
					                   _mm256_load_pd(&re[i+12])));
			 _mm256_store_pd(&c.m_y[i+12],
			                     _mm256_add_pd(_mm256_load_pd(&b.m_y[i+12]),
					                   _mm256_load_pd(&im[i+12])));
			 _mm256_store_pd(&c.m_z[i+12],
			                     _mm256_add_pd(_mm256_load_pd(&b.m_z[i+12]),
					                   _mm256_load_pd(&re[i+12])));
			 _mm256_store_pd(&c.m_w[i+12],
			                     _mm256_add_pd(_mm256_load_pd(&b.m_w[i+12]),
					                   _mm256_load_pd(&im[i+12])));
			 
		     }
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma loop_count min(1),avg(2),max(3)
#endif
                      for(; last_i != size; ++last_i) {
                          c.m_x[last_i] = b.m_x[last_i]+re[last_i];
			  c.m_y[last_i] = b.m_y[last_i]+im[last_i];
			  c.m_z[last_i] = b.m_z[last_i]+re[last_i];
			  c.m_w[last_i] = b.m_w[last_i]+im[last_i];
		      }
		     
	      }


	       __ATTR_ALWAYS_INLINE__
	       __ATTR_HOT__
	       __ATTR_ALIGN__(32)
	       static inline
	       void AVXQuat_Add_AVXQuat(AVXQuatArray1D &b,
	                                const AVXQuatArray1D &a,
	                                const double * __restrict __ATTR_ALIGN__(32) re,
					const int32_t n) {

                     if(__builtin_expect(b.m_size != a.m_size,0) ||
		        __builtin_expect(b.m_size != n,0)) { return;}
		     __m256d re0,re4,re8,re12;
		     re0 = _mm256_setzero_pd();
		     re4 = re0;
		     re8 = re0;
		     re12= re0;
		     int32_t i, last_i;
		     last_i = 0;
		     const int32_t size = a.m_size;
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma code_align(32)
#endif
#pragma omp parallel for schedule(static,32) default(none) private(i,re0,re4,re8,re12) \
                      lastprivate(last_i) shared(a,b,re,size) if(size>=5000)
		      for(i = 0; i != ROUND_TO_FOUR(size,4); i += 16) {
                          last_i = i;
			  _mm_prefetch((const char*)&a.m_x[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_y[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_z[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_w[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&re[i+16],_MM_HINT_T0);
			  re0 = _mm256_load_pd(&re[i+0]);
			  _mm256_store_pd(&b.m_x[i+0],
			                     _mm256_add_pd(_mm256_load_pd(&a.m_x[i+0]),re0));
					                  
			  _mm256_store_pd(&b.m_y[i+0],
			                     _mm256_add_pd(_mm256_load_pd(&a.m_y[i+0]),re0));
					                  
			  _mm256_store_pd(&b.m_z[i+0],
			                     _mm256_add_pd(_mm256_load_pd(&a.m_z[i+0]),re0));
					                   
			  _mm256_store_pd(&b.m_w[i+0],
			                     _mm256_add_pd(_mm256_load_pd(&a.m_w[i+0]),re0));
			  re4 = _mm256_load_pd(&re[i+4]);		                  
		          _mm256_store_pd(&b.m_x[i+4],
			                     _mm256_add_pd(_mm256_load_pd(&a.m_x[i+4]),re4));
					                  
			  _mm256_store_pd(&b.m_y[i+4],
			                     _mm256_add_pd(_mm256_load_pd(&a.m_y[i+4]),re4));
					                  
			  _mm256_store_pd(&b.m_z[i+4],
			                     _mm256_add_pd(_mm256_load_pd(&a.m_z[i+4]),re4));
					                   
			  _mm256_store_pd(&b.m_w[i+4],
			                     _mm256_add_pd(_mm256_load_pd(&a.m_w[i+4]),re4));
			  re8 = _mm256_load_pd(&re[i+8]);
			  _mm256_store_pd(&b.m_x[i+8],
			                     _mm256_add_pd(_mm256_load_pd(&a.m_x[i+8]),re8));
					                  
			  _mm256_store_pd(&b.m_y[i+8],
			                     _mm256_add_pd(_mm256_load_pd(&a.m_y[i+8]),re8));
					                  
			  _mm256_store_pd(&b.m_z[i+8],
			                     _mm256_add_pd(_mm256_load_pd(&a.m_z[i+8]),re8));
					                   
			  _mm256_store_pd(&b.m_w[i+8],
			                     _mm256_add_pd(_mm256_load_pd(&a.m_w[i+8]),re8));
			  re12 = _mm256_load_pd(&re[i+12]);
			  _mm256_store_pd(&b.m_x[i+12],
			                     _mm256_add_pd(_mm256_load_pd(&a.m_x[i+12]),re12));
					                  
			  _mm256_store_pd(&b.m_y[i+12],
			                     _mm256_add_pd(_mm256_load_pd(&a.m_y[i+12]),re12));
					                  
			  _mm256_store_pd(&b.m_z[i+12],
			                     _mm256_add_pd(_mm256_load_pd(&a.m_z[i+12]),re12));
					                   
			  _mm256_store_pd(&b.m_w[i+12],
			                     _mm256_add_pd(_mm256_load_pd(&a.m_w[i+12]),re12));

		      }
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma loop_count min(1),avg(2),max(3)
#endif
                      for(; last_i != size; ++last_i) {
                          b.m_x[last_i] = a.m_x[last_i]+re[last_i];
			  b.m_y[last_i] = a.m_y[last_i]+re[last_i];
			  b.m_z[last_i] = a.m_z[last_i]+re[last_i];
			  b.m_w[last_i] = a.m_w[last_i]+re[last_i];
		      }
		      
	      }

              //
	      // In-place quaternionic addition
	      //

	       __ATTR_ALWAYS_INLINE__
	       __ATTR_HOT__
	       __ATTR_ALIGN__(32)
	       static inline
	       void AVXQuat_Add_AVXQuat(AVXQuatArray1D &b,
	                                const AVXQuatArray1D &a) {

                     if(__builtin_expect(b.m_size != a.m_size,0)) {return;}
		     int32_t i,last_i;
		     last_i = 0;
		     const int32_t size = b.m_size;
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma code_align(32)
#endif
#pragma omp parallel for schedule(static,32) default(none) private(i) \
                      lastprivate(last_i) shared(b,a,size) if(size>=5000)
                      for(i = 0; i != ROUND_TO_FOUR(size,4); i += 16) {
                          last_i = i;
			  _mm_prefetch((const char*)&b.m_x[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_x[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&b.m_y[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_y[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&b.m_z[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_z[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&b.m_w[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_w[i+16],_MM_HINT_T0);
			  _mm256_store_pd(&b.m_x[i+0],
			                      _mm256_add_pd(_mm256_load_pd(&b.m_x[i+0]),
					                    _mm256_load_pd(&a.m_x[i+0])));
			  _mm256_store_pd(&b.m_y[i+0],
			                      _mm256_add_pd(_mm256_load_pd(&b.m_y[i+0]),
					                    _mm256_load_pd(&a.m_y[i+0])));
			  _mm256_store_pd(&b.m_z[i+0],
			                      _mm256_add_pd(_mm256_load_pd(&b.m_z[i+0]),
					                    _mm256_load_pd(&a.m_z[i+0])));
			  _mm256_store_pd(&b.m_w[i+0],
			                      _mm256_add_pd(_mm256_load_pd(&b.m_w[i+0]),
					                    _mm256_load_pd(&a.m_w[i+0])));
			  _mm256_store_pd(&b.m_x[i+4],
			                      _mm256_add_pd(_mm256_load_pd(&b.m_x[i+4]),
					                    _mm256_load_pd(&a.m_x[i+4])));
			  _mm256_store_pd(&b.m_y[i+4],
			                      _mm256_add_pd(_mm256_load_pd(&b.m_y[i+4]),
					                    _mm256_load_pd(&a.m_y[i+4])));
			  _mm256_store_pd(&b.m_z[i+4],
			                      _mm256_add_pd(_mm256_load_pd(&b.m_z[i+4]),
					                    _mm256_load_pd(&a.m_z[i+4])));
			  _mm256_store_pd(&b.m_w[i+4],
			                      _mm256_add_pd(_mm256_load_pd(&b.m_w[i+4]),
					                    _mm256_load_pd(&a.m_w[i+4])));
			  _mm256_store_pd(&b.m_x[i+8],
			                      _mm256_add_pd(_mm256_load_pd(&b.m_x[i+8]),
					                    _mm256_load_pd(&a.m_x[i+8])));
			  _mm256_store_pd(&b.m_y[i+8],
			                      _mm256_add_pd(_mm256_load_pd(&b.m_y[i+8]),
					                    _mm256_load_pd(&a.m_y[i+8])));
			  _mm256_store_pd(&b.m_z[i+8],
			                      _mm256_add_pd(_mm256_load_pd(&b.m_z[i+8]),
					                    _mm256_load_pd(&a.m_z[i+8])));
			  _mm256_store_pd(&b.m_w[i+8],
			                      _mm256_add_pd(_mm256_load_pd(&b.m_w[i+8]),
					                    _mm256_load_pd(&a.m_w[i+8])));
			  _mm256_store_pd(&b.m_x[i+12],
			                      _mm256_add_pd(_mm256_load_pd(&b.m_x[i+12]),
					                    _mm256_load_pd(&a.m_x[i+12])));
			  _mm256_store_pd(&b.m_y[i+12],
			                      _mm256_add_pd(_mm256_load_pd(&b.m_y[i+12]),
					                    _mm256_load_pd(&a.m_y[i+12])));
			  _mm256_store_pd(&b.m_z[i+12],
			                      _mm256_add_pd(_mm256_load_pd(&b.m_z[i+12]),
					                    _mm256_load_pd(&a.m_z[i+12])));
			  _mm256_store_pd(&b.m_w[i+12],
			                      _mm256_add_pd(_mm256_load_pd(&b.m_w[i+12]),
					                    _mm256_load_pd(&a.m_w[i+12])));
			  
		      }
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma loop_count min(1),avg(2),max(3)
#endif
                      for(; last_i != size; ++last_i) {
                          b.m_x[last_i] += a.m_x[last_i];
			  b.m_y[last_i] += a.m_y[last_i];
			  b.m_z[last_i] += a.m_z[last_i];
			  b.m_w[last_i] += a.m_w[last_i];
		      }
                }
  
	       
	       __ATTR_ALWAYS_INLINE__
	       __ATTR_HOT__
	       __ATTR_ALIGN__(32)
	       static inline
	       void AVXQuat_Add_AVXQuat(AVXQuatArray1D &a,
	                                const double * __restrict __ATTR_ALIGN__(32) re,
					const double * __restrict __ATTR_ALIGN__(32) im,
					const int32_t n) {

                     if(__builtin_expect(a.m_size != n,0)) { return;}
		     __m256d re0,re4,re8,re12;
		     __m256d im0,im4,im8,im12;
		     int32_t i,last_i;
		     last_i = 0;
		     const int32_t size = a.m_size;
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma code_align(32)
#endif
#pragma omp parallel for schedule(static,32) default(none) private(i,re0,re4,re8,re12,im0,im4,im8,im12) \
                      lastprivate(last_i) shared(a,re,im,size) if(size>=5000)
		     for(i = 0; i != ROUND_TO_FOUR(size,4); i += 16) {
                         last_i = i;
			 _mm_prefetch((const char*)&a.m_x[i+16],_MM_HINT_T0);
			 _mm_prefetch((const char*)&re[i+16],   _MM_HINT_T0);
			 _mm_prefetch((const char*)&im[i+16],   _MM_HINT_T0);
			 _mm_prefetch((const char*)&a.m_y[i+16],_MM_HINT_T0);
			 _mm_prefetch((const char*)&a.m_z[i+16],_MM_HINT_T0);
			 _mm_prefetch((const char*)&a.m_w[i+16],_MM_HINT_T0);
			 re0 = _mm256_load_pd(&re[i+0]);
			 im0 = _mm256_load_pd(&im[i+0]);
                         _mm256_store_pd(&a.m_x[i+0],
			                      _mm256_add_pd(_mm256_load_pd(&a.m_x[i+0]),re0));
					                    
			 _mm256_store_pd(&a.m_y[i+0],
			                      _mm256_add_pd(_mm256_load_pd(&a.m_y[i+0]),im0));
					                   
			 _mm256_store_pd(&a.m_z[i+0],
			                      _mm256_add_pd(_mm256_load_pd(&a.m_z[i+0]),re0));
					                   
			 _mm256_store_pd(&a.m_w[i+0],
			                      _mm256_add_pd(_mm256_load_pd(&a.m_w[i+0]),im0));
			 re4 = _mm256_load_pd(&re[i+4]);
			 im4 = _mm256_load_pd(&im[i+4]);
			 _mm256_store_pd(&a.m_x[i+4],
			                      _mm256_add_pd(_mm256_load_pd(&a.m_x[i+4]),re4));
					                    
			 _mm256_store_pd(&a.m_y[i+4],
			                      _mm256_add_pd(_mm256_load_pd(&a.m_y[i+4]),im4));
					                   
			 _mm256_store_pd(&a.m_z[i+4],
			                      _mm256_add_pd(_mm256_load_pd(&a.m_z[i+4]),re4));
					                   
			 _mm256_store_pd(&a.m_w[i+4],
			                      _mm256_add_pd(_mm256_load_pd(&a.m_w[i+4]),im4));
			 re8 = _mm256_load_pd(&re[i+8]);
			 im8 = _mm256_load_pd(&im[i+8]);
			 _mm256_store_pd(&a.m_x[i+8],
			                      _mm256_add_pd(_mm256_load_pd(&a.m_x[i+8]),re8));
					                    
			 _mm256_store_pd(&a.m_y[i+8],
			                      _mm256_add_pd(_mm256_load_pd(&a.m_y[i+8]),im8));
					                   
			 _mm256_store_pd(&a.m_z[i+8],
			                      _mm256_add_pd(_mm256_load_pd(&a.m_z[i+8]),re8));
					                   
			 _mm256_store_pd(&a.m_w[i+8],
			                      _mm256_add_pd(_mm256_load_pd(&a.m_w[i+8]),im8));
			 re12 = _mm256_load_pd(&re[i+12]);
			 im12 = _mm256_load_pd(&im[i+12]);
			 _mm256_store_pd(&a.m_x[i+12],
			                      _mm256_add_pd(_mm256_load_pd(&a.m_x[i+12]),re12));
					                    
			 _mm256_store_pd(&a.m_y[i+12],
			                      _mm256_add_pd(_mm256_load_pd(&a.m_y[i+12]),im12));
					                   
			 _mm256_store_pd(&a.m_z[i+12],
			                      _mm256_add_pd(_mm256_load_pd(&a.m_z[i+12]),re12));
					                   
			 _mm256_store_pd(&a.m_w[i+12],
			                      _mm256_add_pd(_mm256_load_pd(&a.m_w[i+12]),im12));
					                   
		     }
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma loop_count min(1),avg(2),max(3)
#endif
                      for(; last_i != size; ++last_i) {
                           a.m_x[last_i] += re[last_i];
			   a.m_y[last_i] += im[last_i];
			   a.m_z[last_i] += re[last_i];
			   a.m_w[last_i] += im[last_i];
		      }
		     
	      }


	      //
	      // Quaternion array subtraction
	      //

	       __ATTR_ALWAYS_INLINE__
	       __ATTR_HOT__
	       __ATTR_ALIGN__(32)
	       static inline
	       void AVXQuat_Sub_AVXQuat(AVXQuatArray1D &c,
	                                const AVXQuatArray1D &b,
					const AVXQuatArray1D &a) {
					
		      if(__builtin_expect(b.m_size != a.m_size,0) ||
		         __builtin_expect(c.m_size != a.m_size,0)) { return;}
                      int32_t i,last_i;
		      last_i = 0;
		      const int32_t size = c.m_size;
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma code_align(32)
#endif
#pragma omp parallel for schedule(static,32) default(none) private(i) \
                      lastprivate(last_i) shared(c,b,a,size) if(size>=5000)
		      for(i = 0; i != ROUND_TO_FOUR(size,4); i += 16) {

		          last_i = i;
			  _mm_prefetch((const char*)&b.m_x[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_x[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&b.m_y[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_y[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&b.m_z[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_z[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&b.m_w[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_w[i+16],_MM_HINT_T0);
			  _mm256_store_pd(&c.m_x[i+0],
			                  _mm256_sub_pd(_mm256_load_pd(&b.m_x[i+0]),
					                  _mm256_load_pd(&a.m_x[i+0])));
			  _mm256_store_pd(&c.m_y[i+0],
			                  _mm256_sub_pd(_mm256_load_pd(&b.m_y[i+0]),
					                 _mm256_load_pd(&a.m_y[i+0])));
			  _mm256_store_pd(&c.m_z[i+0],
			                  _mm256_sub_pd(_mm256_load_pd(&b.m_z[i+0]),
					                 _mm256_load_pd(&a.m_z[i+0]) ));
			  _mm256_store_pd(&c.m_w[i+0],
			                  _mm256_sub_pd(_mm256_load_pd(&b.m_w[i+0]),
					                _mm256_load_pd(&a.m_w[i+0])));
			  _mm256_store_pd(&c.m_x[i+4],
			                  _mm256_sub_pd(_mm256_load_pd(&b.m_x[i+4]),
					                 _mm256_load_pd(&a.m_x[i+4])));
			  _mm256_store_pd(&c.m_y[i+4],
			                  _mm256_sub_pd(_mm256_load_pd(&b.m_y[i+4]),
					                  _mm256_load_pd(&a.m_y[i+4])));
			  _mm256_store_pd(&c.m_z[i+4],
			                  _mm256_sub_pd(_mm256_load_pd(&b.m_z[i+4]),
					                 _mm256_load_pd(&a.m_z[i+4])));
			  _mm256_store_pd(&c.m_w[i+4],
			                  _mm256_sub_pd(_mm256_load_pd(&b.m_w[i+4]),
					                 _mm256_load_pd(&a.m_w[i+4])));
			  _mm256_store_pd(&c.m_x[i+8],
			                  _mm256_sub_pd(_mm256_load_pd(&b.m_x[i+8]),
					                _mm256_load_pd(&a.m_x[i+8]) ));
			  _mm256_store_pd(&c.m_y[i+8],
			                  _mm256_sub_pd(_mm256_load_pd(&b.m_y[i+8]),
					                _mm256_load_pd(&a.m_y[i+8])));
			  _mm256_store_pd(&c.m_z[i+8],
			                  _mm256_sub_pd(_mm256_load_pd(&b.m_z[i+8]),
					                _mm256_load_pd(&a.m_z[i+8])));
			  _mm256_store_pd(&c.m_w[i+8],
			                  _mm256_sub_pd(_mm256_load_pd(&b.m_w[i+8]),
					                 _mm256_load_pd(&a.m_w[i+8])));
			  _mm256_store_pd(&c.m_x[i+12],
			                  _mm256_sub_pd(_mm256_load_pd(&b.m_x[i+12]),
					                _mm256_load_pd(&a.m_x[i+12])));
			  _mm256_store_pd(&c.m_y[i+12],
			                  _mm256_sub_pd(_mm256_load_pd(&b.m_y[i+12]),
					                _mm256_load_pd(&a.m_y[i+12])));
			  _mm256_store_pd(&c.m_z[i+12],
			                  _mm256_sub_pd(_mm256_load_pd(&b.m_z[i+12]),
					                 _mm256_load_pd(&a.m_z[i+12])));
			  _mm256_store_pd(&c.m_w[i+12],
			                  _mm256_sub_pd(_mm256_load_pd(&b.m_w[i+12]),
					                _mm256_load_pd(&a.m_w[i+12])));
			  
		      }
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma loop_count min(1),avg(2),max(3)
#endif
                      for(; last_i != size; ++last_i) {
                           c.m_x[last_i] = b.m_x[last_i]-a.m_x[last_i];
			   c.m_y[last_i] = b.m_y[last_i]-a.m_y[last_i];
			   c.m_z[last_i] = b.m_z[last_i]-a.m_z[last_i];
			   c.m_w[last_i] = b.m_w[last_i]-a.m_w[last_i];
		      }
	       }


               __ATTR_ALWAYS_INLINE__
	       __ATTR_HOT__
	       __ATTR_ALIGN__(32)
	       static inline
	       void AVXQuat_Sub_AVXQuat(AVXQuatArray1D &c,
	                                AVXQuatArray1D &b,
	                                const double * __restrict __ATTR_ALIGN__(32) re,
					const double * __restrict __ATTR_ALIGN__(32) im,
					const int32_t n) {

		     if(__builtin_expect(c.m_size != n,0) ||
		        __builtin_expect(c.m_size != b.m_size,0)) { return;}
                        int32_t i,last_i;
		        last_i = 0;
			const int32_t size = c.m_size;
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma code_align(32)
#endif
#pragma omp parallel for schedule(static,32) default(none) private(i) \
                      lastprivate(last_i) shared(c,b,re,im,size) if(size>=5000)
                     for(i = 0; i != ROUND_TO_FOUR(size,4); i += 16) {
                         last_i = i;
			 _mm_prefetch((const char*)&b.m_x[i+16],_MM_HINT_T0);
			 _mm_prefetch((const char*)&b.m_y[i+16],_MM_HINT_T0);
			 _mm_prefetch((const char*)&b.m_z[i+16],_MM_HINT_T0);
			 _mm_prefetch((const char*)&b.m_w[i+16],_MM_HINT_T0);
			 _mm_prefetch((const char*)&re[i+16],_MM_HINT_T0);
			 _mm_prefetch((const char*)&im[i+16],_MM_HINT_T0);
			 _mm256_store_pd(&c.m_x[i+0],
			                     _mm256_sub_pd(_mm256_load_pd(&b.m_x[i+0]),
					                   _mm256_load_pd(&re[i+0])));
			 _mm256_store_pd(&c.m_y[i+0],
			                     _mm256_sub_pd(_mm256_load_pd(&b.m_y[i+0]),
					                   _mm256_load_pd(&im[i+0])));
			 _mm256_store_pd(&c.m_z[i+0],
			                     _mm256_sub_pd(_mm256_load_pd(&b.m_z[i+0]),
					                   _mm256_load_pd(&re[i+0])));
			 _mm256_store_pd(&c.m_w[i+0],
			                     _mm256_sub_pd(_mm256_load_pd(&b.m_w[i+0]),
					                   _mm256_load_pd(&im[i+0])));
			 _mm256_store_pd(&c.m_x[i+4],
			                     _mm256_sub_pd(_mm256_load_pd(&b.m_x[i+4]),
					                   _mm256_load_pd(&re[i+4])));
			 _mm256_store_pd(&c.m_y[i+4],
			                     _mm256_sub_pd(_mm256_load_pd(&b.m_y[i+4]),
					                   _mm256_load_pd(&im[i+4])));
			 _mm256_store_pd(&c.m_z[i+4],
			                     _mm256_sub_pd(_mm256_load_pd(&b.m_z[i+4]),
					                   _mm256_load_pd(&re[i+4])));
			 _mm256_store_pd(&c.m_w[i+4],
			                     _mm256_sub_pd(_mm256_load_pd(&b.m_w[i+4]),
					                   _mm256_load_pd(&im[i+4])));
			 _mm256_store_pd(&c.m_x[i+8],
			                     _mm256_sub_pd(_mm256_load_pd(&b.m_x[i+8]),
					                   _mm256_load_pd(&re[i+8])));
			 _mm256_store_pd(&c.m_y[i+8],
			                     _mm256_sub_pd(_mm256_load_pd(&b.m_y[i+8]),
					                   _mm256_load_pd(&im[i+8])));
			 _mm256_store_pd(&c.m_z[i+8],
			                     _mm256_sub_pd(_mm256_load_pd(&b.m_z[i+8]),
					                   _mm256_load_pd(&re[i+8])));
			 _mm256_store_pd(&c.m_w[i+8],
			                     _mm256_sub_pd(_mm256_load_pd(&b.m_w[i+8]),
					                   _mm256_load_pd(&im[i+8])));
			 _mm256_store_pd(&c.m_x[i+12],
			                     _mm256_sub_pd(_mm256_load_pd(&b.m_x[i+12]),
					                   _mm256_load_pd(&re[i+12])));
			 _mm256_store_pd(&c.m_y[i+12],
			                     _mm256_sub_pd(_mm256_load_pd(&b.m_y[i+12]),
					                   _mm256_load_pd(&im[i+12])));
			 _mm256_store_pd(&c.m_z[i+12],
			                     _mm256_sub_pd(_mm256_load_pd(&b.m_z[i+12]),
					                   _mm256_load_pd(&re[i+12])));
			 _mm256_store_pd(&c.m_w[i+12],
			                     _mm256_sub_pd(_mm256_load_pd(&b.m_w[i+12]),
					                   _mm256_load_pd(&im[i+12])));
			 
		     }
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma loop_count min(1),avg(2),max(3)
#endif
                      for(; last_i != size; ++last_i) {
                          c.m_x[last_i] = b.m_x[last_i]-re[last_i];
			  c.m_y[last_i] = b.m_y[last_i]-im[last_i];
			  c.m_z[last_i] = b.m_z[last_i]-re[last_i];
			  c.m_w[last_i] = b.m_w[last_i]-im[last_i];
		      }
		     
	      }


	       __ATTR_ALWAYS_INLINE__
	       __ATTR_HOT__
	       __ATTR_ALIGN__(32)
	       static inline
	       void AVXQuat_Sub_AVXQuat(AVXQuatArray1D &b,
	                                const AVXQuatArray1D &a,
	                                const double * __restrict __ATTR_ALIGN__(32) re,
					const int32_t n) {

                     if(__builtin_expect(b.m_size != a.m_size,0) ||
		        __builtin_expect(b.m_size != n,0)) { return;}
		     __m256d re0,re4,re8,re12;
		     re0 = _mm256_setzero_pd();
		     re4 = re0;
		     re8 = re0;
		     re12= re0;
		     int32_t i, last_i;
		     last_i = 0;
		     const int32_t size = a.m_size;
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma code_align(32)
#endif
#pragma omp parallel for schedule(static,32) default(none) private(i,re0,re4,re8,re12) \
                      lastprivate(last_i) shared(a,b,re,size) if(size>=5000)
		      for(i = 0; i != ROUND_TO_FOUR(size,4); i += 16) {
                          last_i = i;
			  _mm_prefetch((const char*)&a.m_x[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_y[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_z[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_w[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&re[i+16],_MM_HINT_T0);
			  re0 = _mm256_load_pd(&re[i+0]);
			  _mm256_store_pd(&b.m_x[i+0],
			                     _mm256_sub_pd(_mm256_load_pd(&a.m_x[i+0]),re0));
					                  
			  _mm256_store_pd(&b.m_y[i+0],
			                     _mm256_sub_pd(_mm256_load_pd(&a.m_y[i+0]),re0));
					                  
			  _mm256_store_pd(&b.m_z[i+0],
			                     _mm256_sub_pd(_mm256_load_pd(&a.m_z[i+0]),re0));
					                   
			  _mm256_store_pd(&b.m_w[i+0],
			                     _mm256_sub_pd(_mm256_load_pd(&a.m_w[i+0]),re0));
			  re4 = _mm256_load_pd(&re[i+4]);		                  
		          _mm256_store_pd(&b.m_x[i+4],
			                     _mm256_sub_pd(_mm256_load_pd(&a.m_x[i+4]),re4));
					                  
			  _mm256_store_pd(&b.m_y[i+4],
			                     _mm256_sub_pd(_mm256_load_pd(&a.m_y[i+4]),re4));
					                  
			  _mm256_store_pd(&b.m_z[i+4],
			                     _mm256_sub_pd(_mm256_load_pd(&a.m_z[i+4]),re4));
					                   
			  _mm256_store_pd(&b.m_w[i+4],
			                     _mm256_sub_pd(_mm256_load_pd(&a.m_w[i+4]),re4));
			  re8 = _mm256_load_pd(&re[i+8]);
			  _mm256_store_pd(&b.m_x[i+8],
			                     _mm256_sub_pd(_mm256_load_pd(&a.m_x[i+8]),re8));
					                  
			  _mm256_store_pd(&b.m_y[i+8],
			                     _mm256_sub_pd(_mm256_load_pd(&a.m_y[i+8]),re8));
					                  
			  _mm256_store_pd(&b.m_z[i+8],
			                     _mm256_sub_pd(_mm256_load_pd(&a.m_z[i+8]),re8));
					                   
			  _mm256_store_pd(&b.m_w[i+8],
			                     _mm256_sub_pd(_mm256_load_pd(&a.m_w[i+8]),re8));
			  re12 = _mm256_load_pd(&re[i+12]);
			  _mm256_store_pd(&b.m_x[i+12],
			                     _mm256_sub_pd(_mm256_load_pd(&a.m_x[i+12]),re12));
					                  
			  _mm256_store_pd(&b.m_y[i+12],
			                     _mm256_sub_pd(_mm256_load_pd(&a.m_y[i+12]),re12));
					                  
			  _mm256_store_pd(&b.m_z[i+12],
			                     _mm256_sub_pd(_mm256_load_pd(&a.m_z[i+12]),re12));
					                   
			  _mm256_store_pd(&b.m_w[i+12],
			                     _mm256_sub_pd(_mm256_load_pd(&a.m_w[i+12]),re12));

		      }
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma loop_count min(1),avg(2),max(3)
#endif
                      for(; last_i != size; ++last_i) {
                          b.m_x[last_i] = a.m_x[last_i]-re[last_i];
			  b.m_y[last_i] = a.m_y[last_i]-re[last_i];
			  b.m_z[last_i] = a.m_z[last_i]-re[last_i];
			  b.m_w[last_i] = a.m_w[last_i]-re[last_i];
		      }
		      
	      }

              //
	      // In-place quaternionic addition
	      //

	       __ATTR_ALWAYS_INLINE__
	       __ATTR_HOT__
	       __ATTR_ALIGN__(32)
	       static inline
	       void AVXQuat_Sub_AVXQuat(AVXQuatArray1D &b,
	                                const AVXQuatArray1D &a) {

                     if(__builtin_expect(b.m_size != a.m_size,0)) {return;}
		     int32_t i,last_i;
		     last_i = 0;
		     const int32_t size = b.m_size;
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma code_align(32)
#endif
#pragma omp parallel for schedule(static,32) default(none) private(i) \
                      lastprivate(last_i) shared(b,a,size) if(size>=5000)
                      for(i = 0; i != ROUND_TO_FOUR(size,4); i += 16) {
                          last_i = i;
			  _mm_prefetch((const char*)&b.m_x[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_x[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&b.m_y[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_y[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&b.m_z[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_z[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&b.m_w[i+16],_MM_HINT_T0);
			  _mm_prefetch((const char*)&a.m_w[i+16],_MM_HINT_T0);
			  _mm256_store_pd(&b.m_x[i+0],
			                      _mm256_sub_pd(_mm256_load_pd(&b.m_x[i+0]),
					                    _mm256_load_pd(&a.m_x[i+0])));
			  _mm256_store_pd(&b.m_y[i+0],
			                      _mm256_sub_pd(_mm256_load_pd(&b.m_y[i+0]),
					                    _mm256_load_pd(&a.m_y[i+0])));
			  _mm256_store_pd(&b.m_z[i+0],
			                      _mm256_sub_pd(_mm256_load_pd(&b.m_z[i+0]),
					                    _mm256_load_pd(&a.m_z[i+0])));
			  _mm256_store_pd(&b.m_w[i+0],
			                      _mm256_sub_pd(_mm256_load_pd(&b.m_w[i+0]),
					                    _mm256_load_pd(&a.m_w[i+0])));
			  _mm256_store_pd(&b.m_x[i+4],
			                      _mm256_sub_pd(_mm256_load_pd(&b.m_x[i+4]),
					                    _mm256_load_pd(&a.m_x[i+4])));
			  _mm256_store_pd(&b.m_y[i+4],
			                      _mm256_sub_pd(_mm256_load_pd(&b.m_y[i+4]),
					                    _mm256_load_pd(&a.m_y[i+4])));
			  _mm256_store_pd(&b.m_z[i+4],
			                      _mm256_sub_pd(_mm256_load_pd(&b.m_z[i+4]),
					                    _mm256_load_pd(&a.m_z[i+4])));
			  _mm256_store_pd(&b.m_w[i+4],
			                      _mm256_sub_pd(_mm256_load_pd(&b.m_w[i+4]),
					                    _mm256_load_pd(&a.m_w[i+4])));
			  _mm256_store_pd(&b.m_x[i+8],
			                      _mm256_sub_pd(_mm256_load_pd(&b.m_x[i+8]),
					                    _mm256_load_pd(&a.m_x[i+8])));
			  _mm256_store_pd(&b.m_y[i+8],
			                      _mm256_sub_pd(_mm256_load_pd(&b.m_y[i+8]),
					                    _mm256_load_pd(&a.m_y[i+8])));
			  _mm256_store_pd(&b.m_z[i+8],
			                      _mm256_sub_pd(_mm256_load_pd(&b.m_z[i+8]),
					                    _mm256_load_pd(&a.m_z[i+8])));
			  _mm256_store_pd(&b.m_w[i+8],
			                      _mm256_sub_pd(_mm256_load_pd(&b.m_w[i+8]),
					                    _mm256_load_pd(&a.m_w[i+8])));
			  _mm256_store_pd(&b.m_x[i+12],
			                      _mm256_sub_pd(_mm256_load_pd(&b.m_x[i+12]),
					                    _mm256_load_pd(&a.m_x[i+12])));
			  _mm256_store_pd(&b.m_y[i+12],
			                      _mm256_sub_pd(_mm256_load_pd(&b.m_y[i+12]),
					                    _mm256_load_pd(&a.m_y[i+12])));
			  _mm256_store_pd(&b.m_z[i+12],
			                      _mm256_sub_pd(_mm256_load_pd(&b.m_z[i+12]),
					                    _mm256_load_pd(&a.m_z[i+12])));
			  _mm256_store_pd(&b.m_w[i+12],
			                      _mm256_sub_pd(_mm256_load_pd(&b.m_w[i+12]),
					                    _mm256_load_pd(&a.m_w[i+12])));
			  
		      }
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma loop_count min(1),avg(2),max(3)
#endif
                      for(; last_i != size; ++last_i) {
                          b.m_x[last_i] -= a.m_x[last_i];
			  b.m_y[last_i] -= a.m_y[last_i];
			  b.m_z[last_i] -= a.m_z[last_i];
			  b.m_w[last_i] -= a.m_w[last_i];
		      }
                }
  
	       
	       __ATTR_ALWAYS_INLINE__
	       __ATTR_HOT__
	       __ATTR_ALIGN__(32)
	       static inline
	       void AVXQuat_Sub_AVXQuat(AVXQuatArray1D &a,
	                                const double * __restrict __ATTR_ALIGN__(32) re,
					const double * __restrict __ATTR_ALIGN__(32) im,
					const int32_t n) {

                     if(__builtin_expect(a.m_size != n,0)) { return;}
		     __m256d re0,re4,re8,re12;
		     __m256d im0,im4,im8,im12;
		     int32_t i,last_i;
		     last_i = 0;
		     const int32_t size = a.m_size;
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma code_align(32)
#endif
#pragma omp parallel for schedule(static,32) default(none) private(i,re0,re4,re8,re12,im0,im4,im8,im12) \
                      lastprivate(last_i) shared(a,re,im,size) if(size>=5000)
		     for(i = 0; i != ROUND_TO_FOUR(size,4); i += 16) {
                         last_i = i;
			 _mm_prefetch((const char*)&a.m_x[i+16],_MM_HINT_T0);
			 _mm_prefetch((const char*)&re[i+16],   _MM_HINT_T0);
			 _mm_prefetch((const char*)&im[i+16],   _MM_HINT_T0);
			 _mm_prefetch((const char*)&a.m_y[i+16],_MM_HINT_T0);
			 _mm_prefetch((const char*)&a.m_z[i+16],_MM_HINT_T0);
			 _mm_prefetch((const char*)&a.m_w[i+16],_MM_HINT_T0);
			 re0 = _mm256_load_pd(&re[i+0]);
			 im0 = _mm256_load_pd(&im[i+0]);
                         _mm256_store_pd(&a.m_x[i+0],
			                      _mm256_sub_pd(_mm256_load_pd(&a.m_x[i+0]),re0));
					                    
			 _mm256_store_pd(&a.m_y[i+0],
			                      _mm256_sub_pd(_mm256_load_pd(&a.m_y[i+0]),im0));
					                   
			 _mm256_store_pd(&a.m_z[i+0],
			                      _mm256_sub_pd(_mm256_load_pd(&a.m_z[i+0]),re0));
					                   
			 _mm256_store_pd(&a.m_w[i+0],
			                      _mm256_sub_pd(_mm256_load_pd(&a.m_w[i+0]),im0));
			 re4 = _mm256_load_pd(&re[i+4]);
			 im4 = _mm256_load_pd(&im[i+4]);
			 _mm256_store_pd(&a.m_x[i+4],
			                      _mm256_sub_pd(_mm256_load_pd(&a.m_x[i+4]),re4));
					                    
			 _mm256_store_pd(&a.m_y[i+4],
			                      _mm256_sub_pd(_mm256_load_pd(&a.m_y[i+4]),im4));
					                   
			 _mm256_store_pd(&a.m_z[i+4],
			                      _mm256_sub_pd(_mm256_load_pd(&a.m_z[i+4]),re4));
					                   
			 _mm256_store_pd(&a.m_w[i+4],
			                      _mm256_sub_pd(_mm256_load_pd(&a.m_w[i+4]),im4));
			 re8 = _mm256_load_pd(&re[i+8]);
			 im8 = _mm256_load_pd(&im[i+8]);
			 _mm256_store_pd(&a.m_x[i+8],
			                      _mm256_sub_pd(_mm256_load_pd(&a.m_x[i+8]),re8));
					                    
			 _mm256_store_pd(&a.m_y[i+8],
			                      _mm256_sub_pd(_mm256_load_pd(&a.m_y[i+8]),im8));
					                   
			 _mm256_store_pd(&a.m_z[i+8],
			                      _mm256_sub_pd(_mm256_load_pd(&a.m_z[i+8]),re8));
					                   
			 _mm256_store_pd(&a.m_w[i+8],
			                      _mm256_sub_pd(_mm256_load_pd(&a.m_w[i+8]),im8));
			 re12 = _mm256_load_pd(&re[i+12]);
			 im12 = _mm256_load_pd(&im[i+12]);
			 _mm256_store_pd(&a.m_x[i+12],
			                      _mm256_sub_pd(_mm256_load_pd(&a.m_x[i+12]),re12));
					                    
			 _mm256_store_pd(&a.m_y[i+12],
			                      _mm256_sub_pd(_mm256_load_pd(&a.m_y[i+12]),im12));
					                   
			 _mm256_store_pd(&a.m_z[i+12],
			                      _mm256_sub_pd(_mm256_load_pd(&a.m_z[i+12]),re12));
					                   
			 _mm256_store_pd(&a.m_w[i+12],
			                      _mm256_sub_pd(_mm256_load_pd(&a.m_w[i+12]),im12));
					                   
		     }
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma loop_count min(1),avg(2),max(3)
#endif
                      for(; last_i != size; ++last_i) {
                           a.m_x[last_i] -= re[last_i];
			   a.m_y[last_i] -= im[last_i];
			   a.m_z[last_i] -= re[last_i];
			   a.m_w[last_i] -= im[last_i];
		      }
		     
	      }



	      
	      
	  
     } // math


}// gms






#endif /*__GMS_AVX_QUATERNION_ARRAY1D__HPP__*/
