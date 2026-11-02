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

#ifndef __GMS_SIMD_INTRIN_PERF_HPP__
#define __GMS_SIMD_INTRIN_PERF_HPP__ 011120260730

#include <map>
#include <memory>
#include <string>
#include <cstdint>
#include "GMS_config.h"

namespace file_info 
{

     const unsigned int GMS_SIMD_INTRIN_PERF_MAJOR = 1;
     const unsigned int GMS_SIMD_INTRIN_PERF_MINOR = 1;
     const unsigned int GMS_SIMD_INTRIN_PERF_MICRO = 0;
     const unsigned int GMS_SIMD_INTRIN_PERF_FULLVER =
       1000U*GMS_SIMD_INTRIN_PERF_MAJOR+100U*GMS_SIMD_INTRIN_PERF_MINOR+
       10U*GMS_SIMD_INTRIN_PERF_MICRO;
     const char * const GMS_SIMD_INTRIN_PERF_CREATION_DATE = "01-10-2026 07:30 AM  +00200 (SUN 01 NOV 2026 GMT+2)";
     const char * const GMS_SIMD_INTRIN_PERF_BUILD_DATE    = __DATE__ " " __TIME__;
     static const char  GMS_SIMD_INTRIN_PERF_AUTHOR[]      = "Programmer: Bernard Gingold, beniekg@gmail.com"; 
     const char * const GMS_SIMD_INTRIN_PERF_SYNOPSIS      = "Configuration ADTs for the SIMD intrinsics performance measurements.";
}

namespace gms 
{

namespace hw_profiling
{

typedef std::map<const std::string,std::pair<const std::string,const std::string>> TypeMap;
			 
typedef std::map<const std::string, std::pair<const std::string,double>>           UopUnfusedMap;

typedef std::map<const std::string, double>                                        UopFusedMap;			

typedef std::unique_ptr<UopFusedMap> UptrFusedMap;

typedef std::unique_ptr<UopUnfusedMap> UptrUnfusedMap;

template <typename T> struct intrinsics_config_t 
{

  T CpuCycles;

  T ReciprocalThroughput;

  T RelativeError;

  /**********************************************
      CPU Latency and Throughput Pipeline Cycles.
  ***********************************************/
  static constexpr T Lat5PipeCycles = static_cast<T>(5.0);

  static constexpr T Lat6PipeCycles = static_cast<T>(6.0);

  static constexpr T Lat2PipeCycles = static_cast<T>(2.0);

  static constexpr T Lat4PipeCycles = static_cast<T>(4.0);

  static constexpr T Lat7PipeCycles = static_cast<T>(7.0);

  static constexpr T Thr1PipeCycle = static_cast<T>(1.0);

  static constexpr T Lat3PipeCycles = static_cast<T>(3.0);

  static constexpr T Thr3PipeCycles = static_cast<T>(3.0);

  static constexpr T Thr0_5PipeCycle = static_cast<T>(0.5);

  static constexpr T Lat35PipeCycles = static_cast<T>(35.0);

  static constexpr T Lat43PipeCycles = static_cast<T>(43.0);

  static constexpr T Thr25PipeCycles = static_cast<T>(25.0);

  static constexpr T Thr28PipeCycles = static_cast<T>(28.0);

  static constexpr T Thr44PipeCycles = static_cast<T>(44.0);

  static constexpr T Lat29PipeCycles = static_cast<T>(29.0);

  static constexpr T Lat21PipeCycles = static_cast<T>(21.0);

  static constexpr T Thr24PipeCycles = static_cast<T>(24.0);

  static constexpr T Thr14PipeCycles = static_cast<T>(14.0);

  static constexpr T Thr13PipeCycles = static_cast<T>(13.0);

  /* Used to describe zero - throughput instruction.*/
  static constexpr T Thr0PipeCycle = static_cast<T>(0.0);

  static constexpr T RCPRelError = static_cast<T>(0.0003662109375);

  static constexpr T RSQRTRelError = static_cast<T>(0.0003662109375);

  static constexpr T RelErrorUnk = static_cast<T>(std::numeric_limits<T>::quiet_Nan());
  
  static constexpr bool LoadSource{false};

  static auto is_memory_load() -> std::int32_t { return ((LoadSource) ? 1 : 0); }
};

struct alignas(64) port_decomp_helper_sse_t 
{

  static auto set_port_value(std::string &val1, std::int32_t val2)
      -> std::pair<std::string, std::int32_t> {
    return (std::make_pair(val1, val2));
  }

  static const constexpr auto uops_min_pd_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm_min_pd per Port unfused decomposition:"});
  }

  static const constexpr auto uops_min_ps_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm_min_ps per Port unfused decomposition:"});
  }

  static const constexpr auto uops_max_pd_descr() -> std::string {
    return (std::string{
        "Micro-ops  of _mm_max_pd per Port unfused decomposition: "});
  }

  static const constexpr auto uops_max_ps_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm_max_ps per Port unfused decomposition: "});
  }

  static const constexpr auto uops_ceil_ps_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm_ceil_ps per Port unfused decomposition: "});
  }

  static const constexpr auto uops_ceil_pd_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm_ceil_pd per Port unfused decomposition: "});
  }

  static const constexpr auto uops_floor_ps_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm_floor_ps per Port unfused decomposition: "});
  }

  static const constexpr auto uops_floor_pd_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm_floor_pd per Port unfused decomposition: "});
  }

  static const constexpr auto uops_zero_ps_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm_setzero_ps per Port unfused decomposition: "});
  }

  static const constexpr auto uops_zero_pd_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm_setzero_pd per Port unfused decomposition: "});
  }


};	


struct alignas(64) port_decomp_helper_avx_t 
{

  static auto set_port_value(std::string &val1, std::int32_t val2)
      -> std::pair<std::string, std::int32_t> {
    return (std::make_pair(val1, val2));
  }

  static const constexpr auto uops_min_pd_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm256_min_pd per Port unfused decomposition:"});
  }

  static const constexpr auto uops_min_ps_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm256_min_ps per Port unfused decomposition:"});
  }

  static const constexpr auto uops_max_pd_descr() -> std::string {
    return (std::string{
        "Micro-ops  of _mm256_max_pd per Port unfused decomposition: "});
  }

  static const constexpr auto uops_max_ps_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm256_max_ps per Port unfused decomposition: "});
  }

  static const constexpr auto uops_ceil_ps_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm256_ceil_ps per Port unfused decomposition: "});
  }

  static const constexpr auto uops_ceil_pd_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm256_ceil_pd per Port unfused decomposition: "});
  }

  static const constexpr auto uops_floor_ps_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm256_floor_ps per Port unfused decomposition: "});
  }

  static const constexpr auto uops_floor_pd_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm256_floor_pd per Port unfused decomposition: "});
  }

  static const constexpr auto uops_zero_ps_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm256_setzero_ps per Port unfused decomposition: "});
  }

  static const constexpr auto uops_zero_pd_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm256_setzero_pd per Port unfused decomposition: "});
  }

};	


struct alignas(64) port_decomp_helper_avx512_t 
{

  static auto set_port_value(std::string &val1, std::int32_t val2)
      -> std::pair<std::string, std::int32_t> {
    return (std::make_pair(val1, val2));
  }

  static const constexpr auto uops_min_pd_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm512_min_pd per Port unfused decomposition:"});
  }

  static const constexpr auto uops_min_ps_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm512_min_ps per Port unfused decomposition:"});
  }

  static const constexpr auto uops_max_pd_descr() -> std::string {
    return (std::string{
        "Micro-ops  of _mm512_max_pd per Port unfused decomposition: "});
  }

  static const constexpr auto uops_max_ps_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm512_max_ps per Port unfused decomposition: "});
  }

  static const constexpr auto uops_ceil_ps_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm512_ceil_ps per Port unfused decomposition: "});
  }

  static const constexpr auto uops_ceil_pd_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm512_ceil_pd per Port unfused decomposition: "});
  }

  static const constexpr auto uops_floor_ps_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm512_floor_ps per Port unfused decomposition: "});
  }

  static const constexpr auto uops_floor_pd_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm512_floor_pd per Port unfused decomposition: "});
  }

  static const constexpr auto uops_zero_ps_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm512_setzero_ps per Port unfused decomposition: "});
  }

  static const constexpr auto uops_zero_pd_descr() -> std::string {
    return (std::string{
        "Micro-ops of _mm512_setzero_pd per Port unfused decomposition: "});
  }

};	

struct RegisterType {

  enum class Types {

    XMM = 0,

    YMM = 1,

    ZMM = 2
  };
};

template <class Type = RegisterType>
struct Register_XMM final : public RegisterType {

  enum class reg_type { XMM = Type::Types::XMM };

  static constexpr std::uint32_t m_BitWidth = 1 << 7;

  static constexpr std::uint32_t m_NumScalarsSP = 1 << 2;

  static constexpr std::uint32_t m_NumScalarsDP = 1 << 1;

  static auto getBitWidth() -> constexpr std::uint32_t { return (m_BitWidth); }

  static auto getNumScalarsSP() -> constexpr std::uint32_t {
    return (m_NumScalarsSP);
  }

  static auto getNumScalarsDP() -> constexpr std::uint32_t {
    return (m_NumScalarsDP);
  }
};

template <class Type = RegisterType>
struct Register_YMM final : public RegisterType {

  enum class reg_type { YMM = Type::Types::YMM };

  static constexpr std::uint32_t m_BitWidth = 1 << 8;

  static constexpr std::uint32_t m_NumScalar32bit = 1 << 3;

  static constexpr std::uint32_t m_NumScalar64bit = 1 << 2;

  static auto getBitWidth() -> constexpr std::uint32_t { return (m_BitWidth); }

  static auto getNumScalar32bit() -> constexpr std::uint32_t {
    return (m_NumScalar32bit);
  }

  static auto getNumScalar64bit() -> constexpr std::uint32_t {
    return (m_NumScalar64bit);
  }
};

template <class Type = RegisterType>
struct Register_ZMM final : public RegisterType {

  enum class reg_type { ZMM = Type::Types::ZMM };

  static constexpr std::uint32_t m_BitWidth = 1 << 9;

  static constexpr std::uint32_t m_NumScalar32bit = 1 << 4;

  static constexpr std::uint32_t m_NumScalar64bit = 1 << 3;

  static auto getBitWidth() -> constexpr std::uint32_t { return (m_BitWidth); }

  static auto getNumScalar32bit() -> constexpr std::uint32_t {
    return (m_NumScalar32bit);
  }

  static auto getNumScalar64bit() -> constexpr std::uint32_t {
    return (m_NumScalar64bit);
  }
};

			

} // hw_profiling

} // gms

#endif /*__GMS_SIMD_INTRIN_PERF_HPP__*/