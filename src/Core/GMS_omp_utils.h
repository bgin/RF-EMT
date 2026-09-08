
/*
   The LLVM Project is under the Apache License v2.0 with LLVM Exceptions:
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

/*
 * Copyright (C) Bernard Gingold, 2020-2026 
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

#ifndef __GMS_OMP_UTILS_H__
#define __GMS_OMP_UTILS_H__ 070920260928

#include <cstdint>
#include <sched.h>
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_OMP_UTILS_MAJOR = 1;
     static const unsigned int GMS_OMP_UTILS_MINOR = 1;
     static const unsigned int GMS_OMP_UTILS_MICRO = 0;
     static const unsigned int GMS_OMP_UTILS_FULLVER =
       1000U*GMS_OMP_UTILS_MAJOR+100U*GMS_OMP_UTILS_MINOR+
       10U*GMS_OMP_UTILS_MICRO;
     static const char GMS_OMP_UTILS_CREATION_DATE[] = "07-09-2026 09:31AM +00200 (MON 07 SEP 2026 GMT+2)";
     static const char GMS_OMP_UTILS_BUILD_DATE[]    = __DATE__; 
     static const char GMS_OMP_UTILS_BUILD_TIME[]    = __TIME__;
     static const char GMS_OMP_UTILS_SYNOPSIS[]      = "Collection of the OpenMP utility routines, adapted and modified by Bernard Gingold. \
                                                        Original work is due to LLVM OMP-Lib authors.";

}

namespace gms 
{

namespace common
{

typedef struct __attribute__((aligned(16))) affinity_mask_t 
{
  std::size_t setsize;
  cpu_set_t * set;
}affinity_mask_t;

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t 
set_thread_affinity_and_priority(const std::int32_t,const std::int32_t);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
print_omp_thread_affinity();

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
print_affinity_environment();

/*
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
verify_thread_limit(const char * ,const std::int32_t,const std::int32_t);
*/


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
check_nesting_level(const std::int32_t,const std::int32_t);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void
check_thread_num(const std::int32_t,const std::int32_t);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void
check_num_threads(const std::int32_t,const std::int32_t);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void
check_host_helper();

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void
check_thread_ancestor(const std::int32_t,const std::int32_t);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void
check_process_id(const std::int32_t,const std::int32_t);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void
check_host(const std::int32_t);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
affinity_mask_zero(affinity_mask_t *);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
affinity_mask_t * affinity_mask_allocate();

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
affinity_mask_free(affinity_mask_t *);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
affinity_mask_copy(affinity_mask_t * __restrict__,
                   const affinity_mask_t * __restrict__);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
affinity_mask_set(affinity_mask_t * ,std::int32_t);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
affinity_mask_clear(affinity_mask_t * ,std::int32_t);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t 
affinity_mask_isset(affinity_mask_t,std::int32_t);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t 
affinity_mask_count(affinity_mask_t * ,std::int32_t);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t 
affinity_mask_equal(const affinity_mask_t * __restrict__,
                    const affinity_mask_t * __restrict__);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
get_thread_affinity(affinity_mask_t *);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
set_thread_affinity(const affinity_mask_t * __restrict__);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
affinity_update_snprintf_values(char ** __restrict__,std::size_t * __restrict__,
                                std::size_t,std::size_t * __restrict__);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::size_t 
affinity_mask_snprintf(char * __restrict__,std::size_t,
                       affinity_mask_t * __restrict__);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
display_affinity_environment();



} // common

} // gms

#endif /*__GMS_OMP_UTILS_H__*/