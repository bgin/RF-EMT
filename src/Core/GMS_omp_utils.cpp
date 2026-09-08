
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

#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <sys/types.h>
#include <sys/resource.h> //getpid 
#include <signal.h>
#include <unistd.h>
#include <string.h>
#include <immintrin.h>
#include <omp.h> 
#include "GMS_omp_utils.h"

std::int32_t 
gms::common
::set_thread_affinity_and_priority(const std::int32_t cpu, const std::int32_t priority)
{
    cpu_set_t cpu_set;
    sched_param sp;
    int32_t status{-1};
    CPU_ZERO(&cpu_set);
    CPU_SET(cpu,&cpu_set);
    if(sched_setaffinity(0,sizeof(cpu_set), &cpu_set) < 0) 
    {
         status = 1;
         return status;
    }
    [[maybe_unused]] std::int32_t printf_ret = std::printf("Affinity set to cpu: %d\n",cpu);
    __builtin_memset(&sp,0,sizeof(sp));
    sp.sched_priority = priority; //99
    if((sched_setscheduler(0,SCHED_FIFO,&sp)) == -1) 
    {
        status = 2;
        return status;
    }
    status = 0;
    return status;
}

void 
gms::common
::print_omp_thread_affinity()
{
#define BUFFER_STORE_SIZE 80
#define FORMAT_STORE_SIZE 80
    char default_format[FORMAT_STORE_SIZE];
    char format_specifier[] = "host=%20H tid=%0.4n binds_to=%A";
    char buffer[BUFFER_STORE_SIZE];
    std::size_t nchars{};
    std::size_t diff;
    nchars = omp_get_affinity_format(default_format,(std::size_t)FORMAT_STORE_SIZE);
    diff   = nchars-(std::size_t)FORMAT_STORE_SIZE;
    if(diff>0ull)
        nchars += diff;
    omp_set_affinity_format(format_specifier);
    nchars = omp_capture_affinity(&buffer[0],(std::size_t)BUFFER_STORE_SIZE,NULL);
    [[maybe_unused]] std::int32_t printf_ret = printf("tid=%d affinity:%s\n",omp_get_thread_num(),buffer);
    
}

/*
   Usage: 
   int cl_thread_limit;

  if (argc != 2) {
    usage();
    return 1;
  }
  cl_thread_limit = atoi(argv[1]);

  omp_set_dynamic(0);
  if (omp_get_thread_limit() != cl_thread_limit) {
    fprintf(stderr, "omp_get_thread_limit failed with %d, should be%d\n",
            omp_get_thread_limit(), cl_thread_limit);
    return 1;
  }
  else if (omp_get_max_threads() > cl_thread_limit) {
#if _OPENMP
    int team_size = cl_thread_limit;
#else
    int team_size = 1;
#endif
    omp_set_num_threads(19);
    verify(__FILE__, __LINE__, 1);
#pragma omp parallel
    {
      verify(__FILE__, __LINE__, team_size);
      verify(__FILE__, __LINE__, team_size);
    }
    verify(__FILE__, __LINE__, 1);

    omp_set_nested(1);
#pragma omp parallel num_threads(3)
    {
      verify(__FILE__, __LINE__, 3);
#pragma omp master
#pragma omp parallel num_threads(21)
      {
        verify(__FILE__, __LINE__, team_size-2);
        verify(__FILE__, __LINE__, team_size-2);
      }
    }
    verify(__FILE__, __LINE__, 1);

    return failed;
  } else {
    fprintf(stderr, "This test is not applicable for max num_threads='%d'\n",
            omp_get_max_threads());
    return 0;

*/

/*
void 
gms::common
::verify_thread_limit(const char * file_name,
                      const std::int32_t line_number,
                      const std::int32_t team_size)
{
  std::int32_t num_threads = omp_get_num_threads();
  if (team_size != num_threads) 
  {
#pragma omp critical(A)
    {
      char label[256];
      snprintf(label, sizeof(label), "%s:%d", file_name, line_number);
      failed = 1;
      std::printf("failed: %s: team_size(%d) != omp_get_num_threads(%d)\n",
                 label, team_size, num_threads);
    }
  }
}
*/

#define XSTR(x) #x
#define STR(x) XSTR(x)

#define streqls(s1, s2) (!strcmp(s1, s2))

#define check(condition)\
  if (!(condition)) {\
    fprintf(stderr, "error: %s: %d: " STR(condition) "\n", __FILE__,\
            __LINE__);\
    exit(1);\
  }

__ATTR_ALWAYS_INLINE__ 
static inline 
std::int32_t get_integer()
{
   std::int32_t n,retval;
   constexpr std::int32_t BUFFER_SIZE = 256;
   char buf[BUFFER_SIZE];
   size_t needed = omp_capture_affinity(buf, BUFFER_SIZE, NULL);
   check(needed < BUFFER_SIZE);
   n = sscanf(buf, "%d", &retval);
   check(n == 1);
   return retval;
}

__ATTR_ALWAYS_INLINE__
static inline 
char * get_string()
{
    constexpr std::int32_t BUFFER_SIZE = 256;
    char buf[BUFFER_SIZE];
    size_t needed = omp_capture_affinity(buf, BUFFER_SIZE, NULL);
    check(needed < BUFFER_SIZE);
    return strdup(buf);
}

__ATTR_ALWAYS_INLINE__
static inline
void check_integer(const char * formats[2], std::int32_t(*func)(),
                   const std::int32_t n_threads_outer,
                   const std::int32_t n_threads_inner)
{
    std::int32_t i;
    for (i = 0; i < 2; ++i) 
    {
    omp_set_affinity_format(formats[i]);
    #pragma omp parallel num_threads(n_threads_outer) // was 8
    {
      check(get_integer() == func());
      #pragma omp parallel num_threads(n_threads_inner) // was 3
      {
        check(get_integer() == func());
      }
      check(get_integer() == func());
    }
  }
}

void 
gms::common
::check_nesting_level(const std::int32_t n_threads_outer,
                      const std::int32_t n_threads_inner)
{
    const char* formats[2] = {"%{nesting_level}", "%L"};
    check_integer(formats, omp_get_level,n_threads_outer,n_threads_inner);
                   
}

void 
gms::common 
::check_thread_num(const std::int32_t n_threads_outer,
                   const std::int32_t n_threads_inner)
{
    const char* formats[2] = {"%{thread_num}", "%n"};
    check_integer(formats, omp_get_thread_num,n_threads_outer,n_threads_inner);
}

void 
gms::common
::check_num_threads(const std::int32_t n_threads_outer,
                    const std::int32_t n_threads_inner)
{
    const char* formats[2] = {"%{num_threads}", "%N"};
    check_integer(formats, omp_get_num_threads,n_threads_outer,n_threads_inner);
}

__ATTR_ALWAYS_INLINE__
static inline 
std::int32_t 
ancestor_helper() {return omp_get_ancestor_thread_num(omp_get_level() - 1); }

void 
gms::common 
::check_thread_ancestor(const std::int32_t n_threads_outer,
                        const std::int32_t n_threads_inner)
{
    const char* formats[2] = {"%{ancestor_tnum}", "%a"};
    check_integer(formats, ancestor_helper,n_threads_outer,n_threads_inner);
}

__ATTR_ALWAYS_INLINE__
static inline 
std::int32_t get_pid_helper() { return (static_cast<std::int32_t>(getpid()));}

void 
gms::common
::check_process_id(const std::int32_t n_threads_outer,
                   const std::int32_t n_threads_inner)
{
    const char* formats[2] = {"%{process_id}", "%P"};
    check_integer(formats, get_pid_helper,n_threads_outer,n_threads_inner);
}

void 
gms::common 
::check_host(const std::int32_t n_threads)
{
    typedef std::int32_t buffer_size_t;
    std::int32_t i;
    buffer_size_t buffer_size = 256;
    const char* formats[2] = {"%{host}", "%H"};
    char hostname[256];
    gethostname(hostname, buffer_size);
    for (i = 0; i < 2; ++i)
    {
        omp_set_affinity_format(formats[i]);
        #pragma omp parallel num_threads(n_threads)
        {
           char* host = get_string();
           check(streqls(host, hostname));
           std::free(host);
        }
    }
}

#ifdef __ve__ 
#define AFFINITY_MAX_CPUS 512
#else
#define AFFINITY_MAX_CPUS (32 * 64)
#endif

void
gms::common
::affinity_mask_zero(affinity_mask_t * mask)
{
  CPU_ZERO_S(mask->setsize, mask->set);
}

gms::common::affinity_mask_t * 
gms::common
::affinity_mask_allocate()
{
  size_t setsize = CPU_ALLOC_SIZE(AFFINITY_MAX_CPUS);
  cpu_set_t *set = CPU_ALLOC(AFFINITY_MAX_CPUS);
  affinity_mask_t *retval = (affinity_mask_t *)malloc(sizeof(affinity_mask_t));
  if(retval==NULL) {return (NULL);}
  retval->setsize = setsize;
  retval->set = set;
  affinity_mask_zero(retval);
  return retval;
}

void 
gms::common
::affinity_mask_free(affinity_mask_t * mask) 
{
   CPU_FREE(mask->set);
}

void 
gms::common
::affinity_mask_copy(affinity_mask_t * __restrict__ dest,
                     const affinity_mask_t * __restrict__ src)
{
   memcpy(dest->set, src->set, dest->setsize);
}

void 
gms::common
::affinity_mask_set(affinity_mask_t * mask,std::int32_t cpu)
{
  CPU_SET_S(cpu, mask->setsize, mask->set);
}

void 
gms::common 
::affinity_mask_clear(affinity_mask_t * mask,std::int32_t cpu)
{
  CPU_CLR_S(cpu, mask->setsize, mask->set);
}

std::int32_t 
gms::common
::affinity_mask_isset(const affinity_mask_t * mask,std::int32_t cpu)
{
   return CPU_ISSET_S(cpu, mask->setsize, mask->set);
}

std::int32_t 
gms::common
::affinity_mask_count(const affinity_mask_t * mask)
{
  return CPU_COUNT_S(mask->setsize, mask->set);
}

std::int32_t 
gms::common
::affinity_mask_equal(const affinity_mask_t * mask1,
                      const affinity_mask_t * mask2)
{
  return CPU_EQUAL_S(mask1->setsize, mask1->set, mask2->set);
}

void 
gms::common 
::get_thread_affinity(affinity_mask_t * mask)
{
  if (sched_getaffinity(0, mask->setsize, mask->set) != 0) 
  {
    perror("sched_getaffinity()");
    exit(EXIT_FAILURE);
  }
}

void 
gms::common 
::set_thread_affinity(const affinity_mask_t * mask)
{
  if (sched_setaffinity(0, mask->setsize, mask->set) != 0) 
  {
    perror("sched_setaffinity()");
    exit(EXIT_FAILURE);
  }
}

void 
gms::common
::affinity_update_snprintf_values(char ** __restrict__ ptr,size_t * __restrict__ remaining,
                                  std::size_t n,std::size_t * __restrict__ retval)
{
  if (n > *remaining && *remaining > 0) 
  {
    *ptr += *remaining;
    *remaining = 0;
  } else {
    *ptr += n;
    *remaining -= n;
  }
  *retval += n;
}

std::size_t 
gms::common 
::affinity_mask_snprintf(char * __restrict__ buf,std::size_t bufsize,
                         const affinity_mask_t * __restrict__ mask)
{
  std::int32_t cpu, need_comma, begin, end;
  std::size_t n;
  char *ptr = buf;
  std::size_t remaining = bufsize;
  std::size_t retval = 0;

  n = snprintf(ptr, remaining, "%c", '{');
  affinity_update_snprintf_values(&ptr, &remaining, n, &retval);

  need_comma = 0;
  for (cpu = 0; cpu < AFFINITY_MAX_CPUS; cpu++) {
    if (!affinity_mask_isset(mask, cpu))
      continue;
    if (need_comma) {
      n = snprintf(ptr, remaining, "%c", ',');
      affinity_update_snprintf_values(&ptr, &remaining, n, &retval);
    }
    begin = cpu;
    for (end = begin + 1; end < AFFINITY_MAX_CPUS; ++end) {
      if (!affinity_mask_isset(mask, end))
        break;
    }
    end--;

    if (end - begin >= 2) {
      n = snprintf(ptr, remaining, "%d-%d", begin, end);
      affinity_update_snprintf_values(&ptr, &remaining, n, &retval);
    } else if (end - begin == 1) {
      n = snprintf(ptr, remaining, "%d,%d", begin, end);
      affinity_update_snprintf_values(&ptr, &remaining, n, &retval);
    } else if (end - begin == 0) {
      n = snprintf(ptr, remaining, "%d", begin);
      affinity_update_snprintf_values(&ptr, &remaining, n, &retval);
    }
    need_comma = 1;
    cpu = end;
  }
  n = snprintf(ptr, remaining, "%c", '}');
  affinity_update_snprintf_values(&ptr, &remaining, n, &retval);
  return retval;
}

void 
gms::common
::display_affinity_environment()
{ 
  __attribute__((unused)) std::int32_t printf_ret;
  printf_ret = printf("  Affinity Environment:\n");
  printf_ret = printf("  OMP_PROC_BIND=%s\n", getenv("OMP_PROC_BIND"));
  printf_ret = printf("  OMP_PLACES=%s\n", getenv("OMP_PLACES"));
  printf_ret = printf("  KMP_AFFINITY=%s\n", getenv("KMP_AFFINITY"));
}


