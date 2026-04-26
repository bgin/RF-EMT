

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

#ifndef __GMS_EARTH_MOTION_HPP__
#define __GMS_EARTH_MOTION_HPP__ 200420251241


namespace file_info {

     const unsigned int GMS_EARTH_MOTION_MAJOR = 2;
     const unsigned int GMS_EARTH_MOTION_MINOR = 0;
     const unsigned int GMS_EARTH_MOTION_MICRO = 0;
     const unsigned int GMS_EARTH_MOTION_FULLVER =
       1000U*GMS_EARTH_MOTION_MAJOR+100U*GMS_EARTH_MOTION_MINOR+
       10U*GMS_EARTH_MOTION_MICRO;
     const char * const GMS_EARTH_MOTION_CREATION_DATE = "20-04-2025 12:41 +00200 (SUN 20 APR 2025 GMT+2)";
     const char * const GMS_EARTH_MOTION_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_EARTH_MOTION_SYNOPSIS      = "Dynamically allocated, 64-byte aligned Earth Motion (single-precision) type.";

}

#include <cstdint>
#include <cassert>
#include <fstream>
#include <iomanip>
#include "GMS_config.h"
#include "GMS_dyn_array.h"
#include "GMS_sse_memset.h"

#if !defined(EARTH_MOTION_USE_PMC_INSTRUMENTATION)
#define EARTH_MOTION_USE_PMC_INSTRUMENTATION 0
#endif 

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(EARTH_MOTION_INIT_STORAGE)
#define EARTH_MOTION_INIT_STORAGE 1
#endif 
#endif

#if (EARTH_MOTION_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_GMS_EARTH_MOTION_NT_STORES)
#define USE_GMS_EARTH_MOTION_NT_STORES 0
#endif

#if (EARTH_MOTION_USE_PMC_INSTRUMENTATION) == 1
#include "GMS_hw_perf_macros.h"

#define PMC_VARS                      \
uint64_t prog_counters_start[4] = {}; \
uint64_t prog_counters_stop[4]  = {}; \
uint64_t tsc_start,tsc_stop;          \
uint64_t act_cyc_start,act_cyc_stop;  \
uint64_t ref_cyc_start,ref_cyc_stop;  \
[[maybe_unused]] uint64_t dummy1;     \
[[maybe_unused]] uint64_t dummy2;     \
[[maybe_unused]] uint64_t dummy3;     \
int32_t core_counter_width;           \
double utilization,nom_ghz,avg_ghz;
#endif 

namespace gms {

    namespace fdm {

               /*
                    Earth Motion angular velocity
                    Bernard Etkin, "Dynamics of Atmospheric Flight", page: 124.
               */
              struct alignas(64) earth_motion_t final 
              {
                     // Only non-zero vector-elements are defined
                     std::size_t m_nsec;   // number of seconds of Earth's rotation
                     darray_r4_t m_OE; //rate of rotation 7.27*10^-5 rad/sec
                     darray_r4_t m_cosLat; // cosine of Latitude
                     darray_r4_t m_sinLat; // sine of Latitude
                     darray_r4_t m_cosLon; // cosine of Longitude
                     darray_r4_t m_sinLon; // sine of Longtitude
                     
                     inline earth_motion_t() noexcept(true)
                     : 
                     m_nsec{},
                     m_OE{darray_r4_t()},
                     m_cosLat{darray_r4_t()},
                     m_sinLat{darray_r4_t()},
                     m_cosLon{darray_r4_t()},
                     m_sinLon{darray_r4_t()}
                     {

                     }

                     inline  earth_motion_t(const std::size_t nsec) noexcept(false)
                     :
                     m_nsec(nsec<=1ull?4ull:nsec),
                     m_OE(darray_r4_t(this->m_nsec)),
                     m_cosLat(darray_r4_t(this->m_nsec)),    
                     m_sinLat(darray_r4_t(this->m_nsec)),
                     m_cosLon(darray_r4_t(this->m_nsec)),
                     m_sinLon(darray_r4_t(this->m_nsec))
                     {

                     }

                     inline earth_motion_t(const earth_motion_t &other) noexcept(false)
                     :
                     m_nsec(other.m_nsec),
                     m_OE(other.m_OE),
                     m_cosLat(other.m_cosLat),
                     m_sinLat(other.m_sinLat),
                     m_cosLon(other.m_cosLon),
                     m_sinLon(other.m_sinLon)
                     {

                     }

                     inline earth_motion_t(earth_motion_t &&other) noexcept(true)
                     :
                     m_nsec(other.m_nsec),
                     m_OE(std::move(other.m_OE)),
                     m_cosLat(std::move(other.m_cosLat)),
                     m_sinLat(std::move(other.m_sinLat)),
                     m_cosLon(std::move(other.m_cosLon)),
                     m_sinLon(std::move(other.m_sinLon))
                     {

                     }

                     inline ~earth_motion_t() noexcept(true)
                     {

                     }

                     inline earth_motion_t &
                     operator=(const earth_motion_t &other) noexcept(false)
                     {
                         if(__builtin_expect(this==&other,0)) { return (*this);}
                         this->m_nsec = other.m_nsec;
                         this->m_OE.operator=(other.m_OE);
                         this->m_cosLat.operator=(other.m_cosLat);
                         this->m_sinLat.operator=(other.m_sinLat);
                         this->m_cosLon.operator=(other.m_cosLon);
                         this->m_sinLon.operator=(other.m_sinLon);
                         return (*this);
                     }

                     inline earth_motion_t &
                     operator=(earth_motion_t &&other) noexcept(true)
                     {
                         if(__builtin_expect(this==&other,0)) { return (*this);}
                         this->m_nsec = other.m_nsec;
                         this->m_OE.operator=(std::move(other.m_OE));
                         this->m_cosLat.operator=(std::move(other.m_cosLat));
                         this->m_sinLat.operator=(std::move(other.m_sinLat));
                         this->m_cosLon.operator=(std::move(other.m_cosLon));
                         this->m_sinLon.operator=(std::move(other.m_sinLon));
                         return (*this);
                     }

                    __ATTR_ALWAYS_INLINE__ 
                    inline void init_OE(const float filler) noexcept(true)
                    {
#if (INIT_BY_STD_FILL) == 0
                      using namespace gms::common;
	                 sse_memset_unroll8x_ps(this->m_OE.m_data,filler,this->m_nsec);
#else 
                      std::fill(this->m_OE.m_data,this->m_OE.m_data+this->m_nsec,filler);
#endif   
                    }

                    __ATTR_ALWAYS_INLINE__ 
                    inline void init_cosLat(const float filler) noexcept(true)
                    {
#if (INIT_BY_STD_FILL) == 0
                      using namespace gms::common;
	                 sse_memset_unroll8x_ps(this->m_cosLat.m_data,filler,this->m_nsec);
#else 
                      std::fill(this->m_cosLat.m_data,this->m_cosLat.m_data+this->m_nsec,filler);
#endif   
                     }

                     __ATTR_ALWAYS_INLINE__ 
                    inline void init_sinLat(const float filler) noexcept(true)
                    {
#if (INIT_BY_STD_FILL) == 0
                      using namespace gms::common;
	                 sse_memset_unroll8x_ps(this->m_sinLat.m_data,filler,this->m_nsec);
#else 
                      std::fill(this->m_sinLat.m_data,this->m_sinLat.m_data+this->m_nsec,filler);
#endif   
                    }

                    __ATTR_ALWAYS_INLINE__ 
                    inline void init_cosLon(const float filler) noexcept(true)
                    {
#if (INIT_BY_STD_FILL) == 0
                      using namespace gms::common;
	                 sse_memset_unroll8x_ps(this->m_cosLon.m_data,filler,this->m_nsec);
#else 
                      std::fill(this->m_cosLon.m_data,this->m_cosLon.m_data+this->m_nsec,filler);
#endif   
                    }

                    __ATTR_ALWAYS_INLINE__ 
                    inline void init_sinLon(const float filler) noexcept(true)
                    {
#if (INIT_BY_STD_FILL) == 0
                      using namespace gms::common;
	                 sse_memset_unroll8x_ps(this->m_sinLon.m_data,filler,this->m_nsec);
#else 
                      std::fill(this->m_sinLon.m_data,this->m_sinLon.m_data+this->m_nsec,filler);
#endif   
                    }

                    static void 
                    create_sequence_plot(const std::uint32_t n_samp,
                                         const float * __restrict sig_arg,
                                         const float * __restrict sig_val,
                                         const std::string &header,
                                         const std::string &title,
                                         const bool is_sig_arg_present)
                    {
                         std::string plot_fname;
                         std::string sig_fname;
                         std::ofstream plot_unit;
                         std::ofstream sig_unit;
                         sig_fname = header+"_plot.txt";
                         sig_unit.open(sig_fname.c_str());
                         if(is_sig_arg_present==true)
                         {
                              for(std::size_t __i{0ull}; __i != n_samp; ++__i)
                              {
                                     sig_unit << " " << sig_arg[__i] << " "
                                                     << sig_val[__i] << "\n";
                              }
                         }
                         else
                         {
                              for(std::size_t __i{0ull}; __i != n_samp; ++__i)
                              {
                                      sig_unit << " " << sig_arg[__i] << "\n";
         
                              }
                         }
                         sig_unit.close();
                         std::cout << "Created signal data file \"" << sig_fname << "\".\n";
                         plot_fname = header+"plot_commands.txt";
                         plot_unit.open(plot_fname.c_str());
                         plot_unit << "#" << plot_fname << "\n";
                         plot_unit << "#\n";
                         plot_unit << "# Usage:\n";
                         plot_unit << "# gnuplot < " << plot_fname << "\n";
                         plot_unit << "#\n";
                         plot_unit << "set term png\n";
                         plot_unit << "set output \"" << header << ".png\"\n";
                         plot_unit << "set xlabel 't'\n";
                         plot_unit << "set ylabel 'y(t)'\n";
                         plot_unit << "set title '" << title << "'\n";
                         plot_unit << "set grid\n";
                         plot_unit << "set style data lines\n";
                         if(is_sig_arg_present==true)
                         {
                                   plot_unit << "plot \"" << sig_fname << "\" using 1:2 lw 1 linecolor rgb \"red\"\n";
                         }
                         else
                         {
                                   plot_unit << "plot \"" << sig_fname << "\" lw 1 linecolor rgb \"red\"\n";
                         }
                         plot_unit << "quit\n";
                         plot_unit.close();
                         std::cout << " Created signal data file \"" << plot_fname << "\"\n";
                    }
                     
                   
                     
                     
              };

               auto 
               operator<<(std::ostream &os,const earth_motion_t &rhs)->std::ostream &
               {
                    std::cout << typeid(rhs).name() << "Begin: object state dump." << std::endl;
                    std::cout << "m_nsec: " << std::cout << rhs.m_nsec << "\n";
                    std::cout << "Dumping the components of: m_OE" << "\n";
                    for(std::size_t i{0ull}; i!=rhs.m_nsec; ++i) 
                    {
                         os << std::fixed << std::setprecision(7) << rhs.m_OE.m_data[i] << "\n";
                    }  
                    std::cout << "Dumping the components of: m_cosLat" << "\n";
                    for(std::size_t i{0ull}; i!=rhs.m_nsec; ++i) 
                    {
                         os << std::fixed << std::setprecision(7) << rhs.m_cosLat.m_data[i] << "\n";
                    }      
                    std::cout << "Dumping the components of: m_sinLat" << "\n";
                    for(std::size_t i{0ull}; i!=rhs.m_nsec; ++i) 
                    {
                         os << std::fixed << std::setprecision(7) << rhs.m_sinLat.m_data[i] << "\n";
                    }  
                    std::cout << "Dumping the components of: m_cosLon" << "\n";
                    for(std::size_t i{0ull}; i!=rhs.m_nsec; ++i) 
                    {
                         os << std::fixed << std::setprecision(7) << rhs.m_cosLon.m_data[i] << "\n";
                    }  
                    std::cout << "Dumping the components of: m_sinLon" << "\n";
                    for(std::size_t i{0ull}; i!=rhs.m_nsec; ++i) 
                    {
                         os << std::fixed << std::setprecision(7) << rhs.m_sinLon.m_data[i] << "\n";
                    } 
                    std::cout << typeid(rhs).name() << "End: object state dump." << std::endl;
                    return (os);
               }
    } // fdm

} // gms











#endif /*__GMS_EARTH_MOTION_R8_HPP__*/