
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

#ifndef __GMS_MACHINE_UTILS_H__
#define __GMS_MACHINE_UTILS_H__ 020720260221


#include <immintrin.h>
#include <cstdint>
#include <sys/sysinfo.h>
#include <sys/utsname.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <cstdio>
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_MACHINE_UTILS_MAJOR = 1;
     static const unsigned int GMS_MACHINE_UTILS_MINOR = 1;
     static const unsigned int GMS_MACHINE_UTILS_MICRO = 0;
     static const unsigned int GMS_MACHINE_UTILS_FULLVER =
       1000U*GMS_MACHINE_UTILS_MAJOR+100U*GMS_MACHINE_UTILS_MINOR+
       10U*GMS_MACHINE_UTILS_MICRO;
     static const char GMS_MACHINE_UTILS_CREATION_DATE[] = "02-07-2026 02:21PM +00200 (THR 02 JUL 2026 GMT+2)";
     static const char GMS_MACHINE_UTILS_BUILD_DATE[]    = __DATE__; 
     static const char GMS_MACHINE_UTILS_BUILD_TIME[]    = __TIME__;
     static const char GMS_MACHINE_UTILS_SYNOPSIS[]      = "Collection of low-level machine and system utilities..";

}

namespace gms 
{

namespace common 
{

__ATTR_ALWAYS_INLINE__
inline static 
std::size_t 
get_physical_memory_size()
{
    struct sysinfo sinfo;
    sysinfo(&sinfo);
    std::size_t phy_mem_size{static_cast<std::size_t>(sinfo.totalram)*
                             static_cast<std::size_t>(sinfo.mem_unit)};
    return (phy_mem_size);
}

__ATTR_ALWAYS_INLINE__
inline static 
std::uint64_t 
rdtscp_fenced()
{
  uint64_t a = 0ull, d = 0ull;
  asm volatile("mfence");
  asm volatile("rdtscp" : "=a"(a), "=d"(d) :: "rcx");
  a = (d << 32) | a;
  asm volatile("mfence");
  return a;
}

__ATTR_ALWAYS_INLINE__
inline static 
std::uint64_t 
rdtsc_fenced_start()
{
   asm volatile("lfence");
   return (__rdtsc());
}

__ATTR_ALWAYS_INLINE__
inline static 
std::uint64_t 
rdtsc_fenced_stop()
{
   std::uint64_t tsc{};
   tsc = __rdtsc();
   asm volatile("lfence");
   return (tsc);
}

__ATTR_ALWAYS_INLINE__
inline static 
std::uint64_t
get_phys_mem_addr(std::size_t vaddr)
{
  int fd = open("/proc/self/pagemap", O_RDONLY);
  uint64_t virtual_addr = (uint64_t)vaddr;
  size_t value = 0;
  off_t offset = (virtual_addr / 4096ull) * sizeof(value);
  int got = pread(fd, &value, sizeof(value), offset);
  if(got != sizeof(value)) 
  {
        return 0ull;
  }
  close(fd);
  return (value << 12) | ((size_t)vaddr & 0xFFFULL);
}

__ATTR_ALWAYS_INLINE__
inline static 
std::size_t 
get_direct_physical_map() 
{
  struct utsname buf;
  uname(&buf);
  std::int32_t major = atoi(strtok(buf.release, "."));
  std::int32_t minor = atoi(strtok(NULL, ".")); 
  if((major == 4 && minor < 19) || major < 4) 
  {
      return 0xffff880000000000ull;
  } 
  else 
  {
        return 0xffff888000000000ull;
  }
}

/* Do not optimize implementation*/

template<class T>
static inline 
void do_not_optimize(T const& val)
{
     asm volatile("" : : "m,r"(val) : "memory");
}

template<class T>
static inline 
void do_not_optimize(T &val)
{
    asm volatile("" : : "m,r"(val) : "memory");
}

__ATTR_ALWAYS_INLINE__
inline static 
void clobber()
{
     asm volatile("" : : : "memory");
}

#define RDTSC_START_v1(cycles)\
    std::uint32_t cycs_high, cycs_low;\
    __asm volatile("cpuid\n"\
                   "rdtsc\n"\
                   "mov %%edx, %0\n"\
                   "mov %%eax, %1"\
                   : "=r"(cycs_high), "=r"(cycs_low)\
                   :\
                   :                              /* no read only */\
                   "%rax", "%rbx", "%rcx", "%rdx" /* clobbers */\
    );\
    (cycles) = ((std::uint64_t)cycs_high << 32) | cycs_low;                             


#define RDTSC_STOP_v1(cycles)\
    std::uint32_t cyce_high, cyce_low;\
    __asm volatile("rdtscp\n"\
                   "mov %%edx, %0\n"\
                   "mov %%eax, %1\n"\
                   "cpuid"\
                   : "=r"(cyce_high), "=r"(cyce_low)\
                   : /* no read only registers */\
                   : "%rax", "%rbx", "%rcx", "%rdx" /* clobbers */\
    );\
    (cycles) = ((std::uint64_t)cyce_high << 32) | cyce_low;  

#define RDTSC_START_v2(cycles)\
  do{\
    std::uint32_t cyc_high, cyc_low;\
    __asm volatile("cpuid\n"\
                   "rdtsc\n"\
                   "mov %%edx, %0\n"\
                   "mov %%eax, %1"\
                   : "=r"(cyc_high), "=r"(cyc_low)\
                   :\
                   :                              /* no read only */\
                   "%rax", "%rbx", "%rcx", "%rdx" /* clobbers */\
    );\
    (cycles) = ((std::uint64_t)cyc_high << 32) | cyc_low;\
  }while(0);
                           


#define RDTSC_STOP_v2(cycles)\
  do{\
    std::uint32_t cyc_high, cyc_low;\
    __asm volatile("rdtscp\n"\
                   "mov %%edx, %0\n"\
                   "mov %%eax, %1\n"\
                   "cpuid"\
                   : "=r"(cyc_high), "=r"(cyc_low)\
                   : /* no read only registers */\
                   : "%rax", "%rbx", "%rcx", "%rdx" /* clobbers */\
    );\
    (cycles) = ((std::uint64_t)cyc_high << 32) | cyc_low;\
  }while(0);

#define RDTSC_START_v3(cycles)\
    __asm volatile("cpuid\n"\
                   "rdtsc\n"\
                   "mov %%edx, %0\n"\
                   "mov %%eax, %1"\
                   : "=r"(cycs_high), "=r"(cycs_low)\
                   :\
                   :                              /* no read only */\
                   "%rax", "%rbx", "%rcx", "%rdx" /* clobbers */\
    );\
    (cycles) = ((std::uint64_t)cycs_high << 32) | cycs_low;                             

#define RDTSC_STOP_v3(cycles)\
    __asm volatile("rdtscp\n"\
                   "mov %%edx, %0\n"\
                   "mov %%eax, %1\n"\
                   "cpuid"\
                   : "=r"(cyce_high), "=r"(cyce_low)\
                   : /* no read only registers */\
                   : "%rax", "%rbx", "%rcx", "%rdx" /* clobbers */\
    );\
    (cycles) = ((std::uint64_t)cyce_high << 32) | cyce_low;  

__ATTR_ALWAYS_INLINE__
inline static 
std::uint64_t 
rdtsc_serialized_start()
{
  std::uint64_t tsc_readout{};
  std::uint32_t cyc_high, cyc_low;
  __asm volatile(  "cpuid\n"\
                   "rdtsc\n"\
                   "mov %%edx, %0\n"\
                   "mov %%eax, %1"\
                   : "=r"(cyc_high), "=r"(cyc_low)\
                   :\
                   :                              /* no read only */\
                   "%rax", "%rbx", "%rcx", "%rdx" /* clobbers */\
    );
    tsc_readout = ((std::uint64_t)cyc_high << 32) | cyc_low;
    return (tsc_readout);
}

__ATTR_ALWAYS_INLINE__
inline static 
std::uint64_t 
rdtsc_serialized_stop()
{
  std::uint64_t tsc_readout{};
  std::uint32_t cyce_high, cyce_low;
  __asm volatile(  "rdtscp\n"\
                   "mov %%edx, %0\n"\
                   "mov %%eax, %1\n"\
                   "cpuid"\
                   : "=r"(cyce_high), "=r"(cyce_low)\
                   : /* no read only registers */\
                   : "%rax", "%rbx", "%rcx", "%rdx" /* clobbers */\
    );
    tsc_readout = ((std::uint64_t)cyce_high << 32) | cyce_low;
    return (tsc_readout);
}

__ATTR_ALWAYS_INLINE__
inline static 
void 
check_prefetch_data_size()
{
  std::int32_t input = 0x2, eax, ebx, ecx, edx;
  asm volatile ("movl %0, %%eax;"::"r"(input));
  asm volatile ("cpuid;");
  asm volatile ("movl %%eax, %0;":"=r" (eax));
  asm volatile ("movl %%ebx, %0;":"=r" (ebx));
  asm volatile ("movl %%ecx, %0;":"=r" (ecx));
  asm volatile ("movl %%edx, %0;":"=r" (edx));
  /* EBX shall contain bits 23-16 = F0 = Prefetch : 64-Byte prefetching*/
  [[maybe_unused]] std::int32_t printf_ret = 
  std::printf("eax = 0x%08.8X, \nebx = 0x%08.8X, \necx = 0x%08.8X, \nedc = 0x%08.8X \n", eax, ebx, ecx, edx);
}

} // common

} // gms


#endif /*__GMS_MACHINE_UTILS_H__*/