

/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef __GMS_EULER_MATRIX_R8_HPP__
#define __GMS_EULER_MATRIX_R8_HPP__ 240420250657

namespace file_info {

     const unsigned int GMS_EULER_MATRIX_R4_MAJOR = 1;
     const unsigned int GMS_EULER_MATRIX_R4_MINOR = 1;
     const unsigned int GMS_EULER_MATRIX_R4_MICRO = 0;
     const unsigned int GMS_EULER_MATRIX_R4_FULLVER =
       1000U*GMS_EULER_MATRIX_R4_MAJOR+100U*GMS_EULER_MATRIX_R4_MINOR+
       10U*GMS_EULER_MATRIX_R4_MICRO;
     const char * const GMS_EULER_MATRIX_R4_CREATION_DATE = "24-04-2025 06:57 AM  +00200 (THR 24 APR 2025 GMT+2)";
     const char * const GMS_EULER_MATRIX_R4_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_EULER_MATRIX_R4_SYNOPSIS      = "Dynamically allocated, 64-byte aligned Euler Matrix (generic representation) double-precision type.";

}

#include <cstdint>
#include <array>
#include <cassert>
#include <vector>
#include "GMS_config.h"
#include "GMS_dyn_array.h"
#include "GMS_sse_memset.h"
#include "GMS_func_plotter.h
#include "GMS_compare_fp_safe.h"

#if !defined(EULER_MATRIX_R8_USE_PMC_INSTRUMENTATION)
#define EULER_MATRIX_R8_USE_PMC_INSTRUMENTATION 0
#endif 

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(EULER_MATRIX_R8_INIT_STORAGE)
#define EULER_MATRIX_R8_INIT_STORAGE 1
#endif 
#endif

#if (EULER_MATRIX_R8_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_EULER_MATRIX_R8_NT_STORES)
#define USE_EULER_MATRIX_R8_NT_STORES 0
#endif

#if (EULER_MATRIX_R8_USE_PMC_INSTRUMENTATION) == 1
#include "GMS_hw_perf_macros.h"

#define PMC_VARS\
      uint64_t prog_counters_start[4] = {};\
      uint64_t prog_counters_stop[4]  = {};\
      uint64_t tsc_start,tsc_stop;\
      uint64_t act_cyc_start,act_cyc_stop;\
      uint64_t ref_cyc_start,ref_cyc_stop;\
      [[maybe_unused]] uint64_t dummy1;\
      [[maybe_unused]] uint64_t dummy2;\
      [[maybe_unused]] uint64_t dummy3;\
      int32_t core_counter_width;\
      double utilization,nom_ghz,avg_ghz;
#endif 

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_EULER_MATRIX_R8_NT_STORES)
#define USE_EULER_MATRIX_R8_NT_SORES 0
#endif

namespace gms {

        namespace radio_navigation {

                       

                      struct alignas(64) em_r8_comparison_holder_t 
                      {
                             std::vector<bool> comp_r1;
                             std::vector<bool> comp_r2;
                             std::vector<bool> comp_r3;
                             std::vector<bool> comp_r4;
                             std::vector<bool> comp_r5;
                             std::vector<bool> comp_r6;
                             std::vector<bool> comp_r7;
                             std::vector<bool> comp_r8;
                             std::vector<bool> comp_r9;
                             double comp_epsilon{};
                      };
                        /*
                              Generic representation of the Euler Matrices (the total number of 12 matrices). 
                              During this class instantiation the proper type will be named and
                              the standalone vectorized procedures will compute the results.
                        */
                       struct alignas(64) euler_matrix_r8_t
                       {
                              std::size_t        m_n;
                              darray_r8_t        m_r1;
                              darray_r8_t        m_r2;
                              darray_r8_t        m_r3;
                              darray_r8_t        m_r4;
                              darray_r8_t        m_r5;
                              darray_r8_t        m_r6;
                              darray_r8_t        m_r7;
                              darray_r8_t        m_r8;
                              darray_r8_t        m_r9;
                              std::array<char,6> m_seq; // Axis rotation sequence, eg. "1,2,3"
              
                              __ATTR_ALWAYS_INLINE__ 
                              inline 
                              euler_matrix_r8_t() noexcept(true)
                              :
                              m_n{},
                              m_r1(),
                              m_r2(),
                              m_r3(),
                              m_r4(),
                              m_r5(),
                              m_r6(),
                              m_r7(),
                              m_r8(),
                              m_r9(),
                              m_seq{}
                              {

                              }

                              __ATTR_ALWAYS_INLINE__ 
                              inline
                              euler_matrix_r8_t(const std::size_t n,
                                                const std::array<char,6> seq) noexcept(false)
                              :
                              m_n(n),
                              m_r1(darray_r8_t(m_n)),
                              m_r2(darray_r8_t(m_n)),
                              m_r3(darray_r8_t(m_n)),
                              m_r4(darray_r8_t(m_n)),
                              m_r5(darray_r8_t(m_n)),
                              m_r6(darray_r8_t(m_n)),
                              m_r7(darray_r8_t(m_n)),
                              m_r8(darray_r8_t(m_n)),
                              m_r9(darray_r8_t(m_n)),
                              m_seq(seq)
                              {

                              }

                              __ATTR_ALWAYS_INLINE__
                              inline 
                              euler_matrix_r8_t(const euler_matrix_r8_t &other) noexcept(false) 
                              :
                              m_n(other.m_n),
                              m_r1(other.m_r1),
                              m_r2(other.m_r2),
                              m_r3(other.m_r3),
                              m_r4(other.m_r4),
                              m_r5(other.m_r5),
                              m_r6(other.m_r6),
                              m_r7(other.m_r7),
                              m_r8(other.m_r8),
                              m_r9(other.m_r9),
                              m_seq(other.m_seq)
                              {

                              }

                              __ATTR_ALWAYS_INLINE__
                              inline 
                              euler_matrix_r8_t(euler_matrix_r8_t &&other) noexcept(true) 
                              :
                              m_n(other.m_n),
                              m_r1(std::move(other.m_r1)),
                              m_r2(std::move(other.m_r2)),
                              m_r3(std::move(other.m_r3)),
                              m_r4(std::move(other.m_r4)),
                              m_r5(std::move(other.m_r5)),
                              m_r6(std::move(other.m_r6)),
                              m_r7(std::move(other.m_r7)),
                              m_r8(std::move(other.m_r8)),
                              m_r9(std::move(other.m_r9)),
                              m_seq(std::move(other.m_seq))
                              {

                              }

                              __ATTR_ALWAYS_INLINE__ 
                              inline 
                              ~euler_matrix_r8_t() noexcept(true) 
                              {

                              }

                              __ATTR_ALWAYS_INLINE__ 
                              inline
                              euler_matrix_r8_t & 
                              operator=(const euler_matrix_r8_t &other) noexcept(false) 
                              {
                                   if(__builtin_expect(this==&other,0)) { return (*this);}
                                   this->m_n  = other.m_n;
                                   this->m_r1 = other.m_r1;
                                   this->m_r2 = other.m_r2;
                                   this->m_r3 = other.m_r3;
                                   this->m_r4 = other.m_r4;
                                   this->m_r5 = other.m_r5;
                                   this->m_r6 = other.m_r6;
                                   this->m_r7 = other.m_r7;
                                   this->m_r8 = other.m_r8;
                                   this->m_r9 = other.m_r9;
                                   this->m_seq= other.m_seq; 
                                   return (*this);
                              }

                              __ATTR_ALWAYS_INLINE__ 
                              inline
                              euler_matrix_r8_t & 
                              operator=(euler_matrix_r8_t &&other) noexcept(true) 
                              {
                                   if(__builtin_expect(this==&other,0)) { return (*this);}
                                   this->m_n  = other.m_n;
                                   this->m_r1 = std::move(other.m_r1);
                                   this->m_r2 = std::move(other.m_r2);
                                   this->m_r3 = std::move(other.m_r3);
                                   this->m_r4 = std::move(other.m_r4);
                                   this->m_r5 = std::move(other.m_r5);
                                   this->m_r6 = std::move(other.m_r6);
                                   this->m_r7 = std::move(other.m_r7);
                                   this->m_r8 = std::move(other.m_r8);
                                   this->m_r9 = std::move(other.m_r9);
                                   this->m_seq= std::move(other.m_seq);
                                   return (*this); 
                              }

                              __ATTR_ALWAYS_INLINE__ 
                              inline
                              void init_r1(const double fill) noexcept(true)
                              {
#if (INIT_BY_STD_FILL) == 0
                                 using namespace gms::common;
	                             sse_memset_unroll8x_pd(this->m_r1.m_data,filler,this->m_n);
#else 
                                 std::fill(this->m_r1.m_data,this->m_r1.m_data+this->m_n,filler);
#endif
                              }

                               __ATTR_ALWAYS_INLINE__ 
                              inline
                              void init_r2(const double fill) noexcept(true)
                              {
#if (INIT_BY_STD_FILL) == 0
                                 using namespace gms::common;
	                             sse_memset_unroll8x_pd(this->m_r2.m_data,filler,this->m_n);
#else 
                                 std::fill(this->m_r2.m_data,this->m_r2.m_data+this->m_n,filler);
#endif
                              }

                               __ATTR_ALWAYS_INLINE__ 
                              inline
                              void init_r3(const double fill) noexcept(true)
                              {
#if (INIT_BY_STD_FILL) == 0
                                 using namespace gms::common;
	                             sse_memset_unroll8x_pd(this->m_r3.m_data,filler,this->m_n);
#else 
                                 std::fill(this->m_r3.m_data,this->m_r3.m_data+this->m_n,filler);
#endif
                              }

                               __ATTR_ALWAYS_INLINE__ 
                              inline
                              void init_r4(const double fill) noexcept(true)
                              {
#if (INIT_BY_STD_FILL) == 0
                                 using namespace gms::common;
	                             sse_memset_unroll8x_pd(this->m_r4.m_data,filler,this->m_n);
#else 
                                 std::fill(this->m_r4.m_data,this->m_r4.m_data+this->m_n,filler);
#endif
                              }

                               __ATTR_ALWAYS_INLINE__ 
                              inline
                              void init_r5(const double fill) noexcept(true)
                              {
#if (INIT_BY_STD_FILL) == 0
                                 using namespace gms::common;
	                             sse_memset_unroll8x_pd(this->m_r5.m_data,filler,this->m_n);
#else 
                                 std::fill(this->m_r5.m_data,this->m_r5.m_data+this->m_n,filler);
#endif
                              }

                               __ATTR_ALWAYS_INLINE__ 
                              inline
                              void init_r6(const double fill) noexcept(true)
                              {
#if (INIT_BY_STD_FILL) == 0
                                 using namespace gms::common;
	                             sse_memset_unroll8x_pd(this->m_r6.m_data,filler,this->m_n);
#else 
                                 std::fill(this->m_r6.m_data,this->m_r6.m_data+this->m_n,filler);
#endif
                              }

                               __ATTR_ALWAYS_INLINE__ 
                              inline
                              void init_r7(const double fill) noexcept(true)
                              {
#if (INIT_BY_STD_FILL) == 0
                                 using namespace gms::common;
	                             sse_memset_unroll8x_pd(this->m_r7.m_data,filler,this->m_n);
#else 
                                 std::fill(this->m_r7.m_data,this->m_r7.m_data+this->m_n,filler);
#endif
                              }

                               __ATTR_ALWAYS_INLINE__ 
                              inline
                              void init_r8(const double fill) noexcept(true)
                              {
#if (INIT_BY_STD_FILL) == 0
                                 using namespace gms::common;
	                             sse_memset_unroll8x_pd(this->m_r8.m_data,filler,this->m_n);
#else 
                                 std::fill(this->m_r8.m_data,this->m_r8.m_data+this->m_n,filler);
#endif
                              }

                               __ATTR_ALWAYS_INLINE__ 
                              inline
                              void init_r9(const double fill) noexcept(true)
                              {
#if (INIT_BY_STD_FILL) == 0
                                 using namespace gms::common;
	                             sse_memset_unroll8x_pd(this->m_r9.m_data,filler,this->m_n);
#else 
                                 std::fill(this->m_r9.m_data,this->m_r9.m_data+this->m_n,filler);
#endif
                              }

                              __ATTR_ALWAYS_INLINE__ 
                              inline 
                              std::int32_t plot_euler_matrix() const 
                              {
                                   std::int32_t ret_status = -9999;
                                   const std::string struct_header(typeid(*this).name());
                                   const std::string r1_axis("Eulerian r1-rotation axis.");
                                   const std::string r2_axis("Eulerian r2-rotation axis.");
                                   const std::string r3_axis("Eulerian r3-rotation axis.");
                                   const std::string r4_axis("Eulerian r4-rotation axis.");
                                   const std::string r5_axis("Eulerian r5-rotation axis.");
                                   const std::string r6_axis("Eulerian r6-rotation axis.");
                                   const std::string r7_axis("Eulerian r7-rotation axis.");
                                   const std::string r8_axis("Eulerian r8-rotation axis.");
                                   const std::string r9_axis("Eulerian r9-rotation axis.");
                                   gms::common::create_function_xy_plot(this->m_n,nullptr,this->m_r1.m_data,struct_header,r1_axis,false);
                                   gms::common::create_function_xy_plot(this->m_n,nullptr,this->m_r2.m_data,struct_header,r2_axis,false);
                                   gms::common::create_function_xy_plot(this->m_n,nullptr,this->m_r3.m_data,struct_header,r3_axis,false);
                                   gms::common::create_function_xy_plot(this->m_n,nullptr,this->m_r4.m_data,struct_header,r4_axis,false);
                                   gms::common::create_function_xy_plot(this->m_n,nullptr,this->m_r5.m_data,struct_header,r5_axis,false);
                                   gms::common::create_function_xy_plot(this->m_n,nullptr,this->m_r6.m_data,struct_header,r6_axis,false);
                                   gms::common::create_function_xy_plot(this->m_n,nullptr,this->m_r7.m_data,struct_header,r7_axis,false);
                                   gms::common::create_function_xy_plot(this->m_n,nullptr,this->m_r8.m_data,struct_header,r8_axis,false);
                                   gms::common::create_function_xy_plot(this->m_n,nullptr,this->m_r9.m_data,struct_header,r9_axis,false);
                                   ret_status = 0;
                                   return (ret_status);
                              }

                              
                                 
                       };


        }
}

#endif /*__GMS_EULER_MATRIX_R8_HPP__*/