
/*MIT License
!Copyright (c) 2020 Bernard Gingold
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
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

#define RDTSC_START(cycles)\
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


#define RDTSC_STOP(cycles)\
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


/*
push    rbx
cpuid
rdtsc
mov %edx, esi
mov %eax, edi
        rdtscp
mov %edx, esi
mov %eax, edi
cpuid
*/

} // common

} // gms


#endif /*__GMS_MACHINE_UTILS_H__*/
