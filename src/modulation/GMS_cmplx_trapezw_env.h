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


#ifndef __GMS_CMPLX_TRAPEZW_ENV_H__
#define __GMS_CPPLX_TRAPEZW_ENV_H__ 220920250641

namespace file_info 
{

     static const unsigned int GMS_CMPLX_TRAPEZW_ENV_MAJOR = 1;
     static const unsigned int GMS_CMPLX_TRAPEZW_ENV_MINOR = 1;
     static const unsigned int GMS_CMPLX_TRAPEZW_ENV_MICRO = 0;
     static const unsigned int GMS_CMPLX_TRAPEZW_ENV_FULLVER =
       1000U*GMS_CMPLX_TRAPEZW_ENV_MAJOR+100U*GMS_CMPLX_TRAPEZW_ENV_MINOR+
       10U*GMS_CMPLX_TRAPEZW_ENV_MICRO;
     static const char GMS_CMPLX_TRAPEZW_ENV_CREATION_DATE[] = "22-09-2025 06:41 +00200 (MON 22 SEP 2025 GMT+2)";
     static const char GMS_CMPLX_TRAPEZW_ENV_BUILD_DATE[]    = __DATE__; 
     static const char GMS_CMPLX_TRAPEZW_ENV_BUILD_TIME[]    = __TIME__;
     static const char GMS_CMPLX_TRAPEZW_ENV_SYNOPSIS[]      = "Complex Trapezoid Envelope Signal.";

}

#include <cstdint>
#include <string>
#include <iostream>
#include "GMS_config.h"
#include "GMS_dyn_array.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (CMPLX_TRAPEZW_ENV_USE_NT_STORES)
#define CMPLX_TRAPEZW_ENV_USE_NT_STORES 0
#endif

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(CMPLX_TRAPEZW_ENV_INIT_STORAGE)
#define CMPLX_TRAPEZW_ENV_INIT_STORAGE 1
#endif 
#endif

#if (CMPLX_TRAPEZW_ENV_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(CMPLX_TRAPEZW_ENV_USE_CEPHES)
#define CMPLX_TRAPEZW_ENV_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (CMPLX_TRAPEZW_ENV_USE_PMC_INSTRUMENTATION)
#define CMPLX_TRAPEZW_ENV_USE_PMC_INSTRUMENTATION 1
#endif 

#if (CMPLX_TRAPEZW_ENV_USE_PMC_INSTRUMENTATION) == 1
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

#if (CMPLX_TRAPEZW_ENV_USE_CEPHES) == 0
#include <cmath>
#else 
#include "GMS_cephes_sin_cos.h"
#endif 


namespace gms 
{

namespace radiolocation 
{
             enum class cmplx_trapezw_env_rand_gens : int32_t 
             {
                     rg_minstd_rand0,
                     rg_minstd_rand,
                     rg_mt19937,
                     rg_mt19937_64,
                     rg_ranlux24_base,
                     rg_ranlux48_base,
                     rg_ranlux24,
                     rg_ranlux48,
                     rg_knuth_b
             };

             struct alignas(64) cmplx_trapezw_env_pdf_params_t 
             {
                   float uni_real_a_r;
                   float uni_real_a_i;
                   float uni_real_b_r;
                   float uni_real_b_i;
                   float norm_mu_r;
                   float norm_mu_i;
                   float norm_sigma_r;
                   float norm_sigma_i;
                   float cauchy_a_r;
                   float cauchy_a_i;
                   float cauchy_b_r;
                   float cauchy_b_i;
                   float log_norm_m_r;
                   float log_norm_m_i;
                   float log_norm_s_r;
                   float log_norm_s_i;
                   float expo_gamma_r;
                   float expo_gamma_i;
                   float weibull_a_r;
                   float weibull_a_i;
                   float weibull_b_r;
                   float weibull_b_i;
                   float gamma_alph_r;
                   float gamma_alph_i;
                   float gamma_bet_r;
                   float gamma_bet_i;
                   int32_t poisson_mu;

            };

            /*
                 Complex trapezoidal waveform envelope (baseband signal)
            */
            struct alignas(64) cmplx_trapezw_env_t final 
            {
                   std::size_t     __I_n_samples__; 
                   std::size_t     __Q_n_samples__; 
                   std::uint32_t   __I_n_K__; // number of inner loop wave components
                   std::uint32_t   __Q_n_K__; // number of inner loop wave components 
                   // envelope parameters
                   float           __I_a__;
                   float           __I_l__;
                   float           __I_c__;
                   float           __I_m__;
                   float           __Q_a__;
                   float           __Q_l__;
                   float           __Q_c__;
                   float           __Q_m__;
                   darray_r4_t     __I_chan__; //In-phase channel
                   darray_r4_t     __Q_chan__; //Quadrature channel

                   cmplx_trapezw_env_t() = delete;

                   cmplx_trapezw_env_t(const std::size_t,
                                       const std::size_t,
                                       const std::uint32_t,
                                       const std::uint32_t,
                                       const float,
                                       const float,
                                       const float,
                                       const float,
                                       const float,
                                       const float,
                                       const float,
                                       const float) noexcept(false);

                    cmplx_trapezw_env_t(cmplx_trapezw_env_t &&);

                    cmplx_trapezw_env_t(const cmplx_trapezw_env_t &) = delete;

                    ~cmplx_trapezw_env_t() noexcept(false);

                    cmplx_trapezw_env_t & operator=(const cmplx_trapezw_env_t &) = delete;

                    cmplx_trapezw_env_t & operator=(cmplx_trapezw_env_t &&);

                    void init_storage(const float,
                                      const float);

                    static void create_signal_plot(const std::uint32_t,
                                                   const float * __restrict,
                                                   const float * __restrict,
                                                   const std::string &,
                                                   const std::string &,
                                                   const bool );

                    template <typename T> std::int32_t signum(T val)
                    {
                             return (T(0)<val)-(val<T(0));
                    }

                    __ATTR_ALWAYS_INLINE__
                    inline float sin_squarew_I_sample(const float t,
                                                     const float invT) // the value of either __I_m__ or __Q_m__
                    {
                          using namespace gms::math;
                          constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                          const float arg{this->__I_m__*C6283185307179586476925286766559*t*invT};
                          float sin_val{};
#if (CMPLX_TRAPEZW_ENV_USE_CEPHES) == 1
                          sin_val = ceph_sinf(arg);
                          return static_cast<float>(signum<float>(sin_val));
#else                     
                          sin_val = std::sin(arg);
                          return static_cast<float>(signum<float>(sin_val));
#endif 
                    }

                     __ATTR_ALWAYS_INLINE__
                    inline float sin_squarew_Q_sample(const float t,
                                                    const float invT) // the value of either __I_m__ or __Q_m__
                    {
                          using namespace gms::math;
                          constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                          const float arg{this->__Q_m__*C6283185307179586476925286766559*t*invT};
                          float sin_val{};
#if (CMPLX_TRAPEZW_ENV_USE_CEPHES) == 1
                          sin_val = ceph_sinf(arg);
                          return static_cast<float>(signum<float>(sin_val));
#else                     
                          sin_val = std::sin(arg);
                          return static_cast<float>(signum<float>(sin_val));
#endif 
                    }



                     __ATTR_ALWAYS_INLINE__
                    inline float cos_squarew_I_sample(const float t,
                                                    const float invT) // the value of either __I_m__ or __Q_m__
                    {
                          using namespace gms::math;
                          constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                          const float arg{this->__I_m__*C6283185307179586476925286766559*t*invT};
                          float cos_val{};
#if (CMPLX_TRAPEZW_ENV_USE_CEPHES) == 1
                          cos_val = ceph_cosf(arg);
                          return static_cast<float>(signum<float>(cos_val));
#else 
                          cos_val = std::cos(arg);
                          return static_cast<float>(signum<float>(cos_val));
#endif 
                    }

                    __ATTR_ALWAYS_INLINE__
                    inline float cos_squarew_Q_sample(const float t,
                                                    const float invT) // the value of either __I_m__ or __Q_m__
                    {
                          using namespace gms::math;
                          constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                          const float arg{this->__Q_m__*C6283185307179586476925286766559*t*invT};
                          float cos_val{};
#if (CMPLX_TRAPEZW_ENV_USE_CEPHES) == 1
                          cos_val = ceph_cosf(arg);
                          return static_cast<float>(signum<float>(cos_val));
#else 
                          cos_val = std::cos(arg);
                          return static_cast<float>(signum<float>(cos_val));
#endif 
                    }

                    __ATTR_ALWAYS_INLINE__ 
                    inline float I_sample(const float t)
                    {
                          using namespace gms::math;
                          constexpr float invPI{0.318309886183790671537767526745f};
                          constexpr float PI{3.14159265358979323846264338328f};
	                     float     a_over_PI{this->__I_a__*invPI};
	                     float     PI_over_m{PI/this->__I_m__};
	                     const float arg{PI_over_m*t+this->__I_l__};
#if (CMPLX_TRAPEZW_ENV_USE_CEPHES) == 1
                           const float t_as{ceph_asinf(ceph_sinf(arg))};
		                 const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                           const float t_as{std::asin(std::sin(arg))};
		                 const float t_ac{std::acos(std::cos(arg))};
#endif 
                           const float t_0{a_over_PI*(t_as+t_ac)};
		                 const float t_1{t_0-5.0f+this->__I_c__};
		                 return (t_1);
                    }
                    
                   __ATTR_ALWAYS_INLINE__ 
                    inline float Q_sample(const float t)
                    {
                          using namespace gms::math;
                          constexpr float invPI{0.318309886183790671537767526745f};
                          constexpr float PI{3.14159265358979323846264338328f};
	                     float     a_over_PI{this->__Q_a__*invPI};
	                     float     PI_over_m{PI/this->__Q_m__};
	                     const float arg{PI_over_m*t+this->__Q_l__};
#if (CMPLX_TRAPEZW_ENV_USE_CEPHES) == 1
                           const float t_as{ceph_asinf(ceph_sinf(arg))};
		               const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                           const float t_as{std::asin(std::sin(arg))};
		               const float t_ac{std::acos(std::cos(arg))};
#endif 
                           const float t_0{a_over_PI*(t_as+t_ac)};
		               const float t_1{t_0-5.0f+this->__Q_c__};
		               return (t_1);
                    }

                    /*
                         Square-wave modulated (no data symbol)
                    */
                    std::int32_t chan_I_squarew_modulated(const std::uint32_t);

                    std::int32_t chan_Q_squarew_modulated(const std::uint32_t);
                    
                    /*
                          Outer loop unrolled four times
                    */
                    std::int32_t chan_I_squarew_modulated_u4x(const std::uint32_t);
                    /*
                          Outer loop unrolled four times
                    */
                    std::int32_t chan_Q_squarew_modulated_u4x(const std::uint32_t);

                    
                    /*  Data symbol-transmitted*/
                    std::int32_t chan_I_data_symbol(const float * __restrict__,
                                                     const std::uint32_t,
                                                     const std::uint32_t);
                    
                     /*  Data symbol-transmitted*/
                    std::int32_t chan_Q_data_symbol(const float * __restrict__,
                                                     const std::uint32_t,
                                                     const std::uint32_t);

                     /*  Data symbol-transmitted, outer loop unrolled four times.*/
                    std::int32_t chan_I_data_symbol_u4x(const float * __restrict__,
                                                        const std::uint32_t,
                                                        const std::uint32_t);
                    
                     /*  Data symbol-transmitted, outer loop unrolled four times.*/
                    std::int32_t chan_Q_data_symbol_u4x(const float * __restrict__,
                                                        const std::uint32_t,
                                                        const std::uint32_t);

                    /*  Data symbol-transmitted random*/
                    
                    template<class Functor> std::int32_t 
                    chan_I_rand_data_symbol(Functor & f)
                    {
                         const float T{static_cast<float>(this->__I_n_samples__)};
                         float sum;
                         for(std::uint32_t __t{0ull}; __t != this->__I_n_samples__; ++__t) 
                         {
                              const float t{static_cast<float>(__t)};
                              sum = 0.0f;
                              for(std::uint32_t __k{0}; __k != this->__I_n_K__; ++__k) 
                              {
                                  const float k{static_cast<float>(__k)};
                                  const float arg{t-k*T};
                                  sum += I_sample(arg)*f();
                              }
                             this->__I_chan__.m_data[__t] = sum;
                         }
                        return (0);
                    }
                    
                     /*  Data symbol-transmitted*/
                    template<class Functor>
                    std::int32_t chan_Q_rand_data_symbol(Functor & f)
                    {
                         const float T{static_cast<float>(this->__Q_n_samples__)};
                         float sum;
                         for(std::uint32_t __t{0ull}; __t != this->__Q_n_samples__; ++__t) 
                         {
                            const float t{static_cast<float>(__t)};
                            sum = 0.0f;
                            for(std::uint32_t __k{0}; __k != this->__Q_n_K__; ++__k) 
                            {
                               const float k{static_cast<float>(__k)};
                               const float arg{t-k*T};
                               sum += Q_sample(arg)*f();
                            }
                            this->__Q_chan__.m_data[__t] = sum;
                         }
                         return (0);
                    }                                 
 
            };

            auto operator<<(std::ostream &,
                            cmplx_trapezw_env_t &)->std::ostream &;


}

}



















#endif /*__GMS_CMPLX_TRAPEZW_ENV_H__*/