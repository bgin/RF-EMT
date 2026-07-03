
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



} // common

} // gms


#endif /*__GMS_MACHINE_UTILS_H__*/