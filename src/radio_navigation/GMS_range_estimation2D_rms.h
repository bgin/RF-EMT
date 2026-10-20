
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

#ifndef __GMS_RANGE_ESTIMATION_2D_RMS_H__
#define __GMS_RANGE_ESTIMATION_2D_RMS_H__

#include <cstdint>
#include <iostream>
#include <vector>
#include <optional>
#include "GMS_config.h"
#include "GMS_dyn_array.h"

namespace file_info 
{

     static const unsigned int GMS_RANGE_ESTIMATION_2D_RMS_MAJOR = 1;
     static const unsigned int GMS_RANGE_ESTIMATION_2D_RMS_MINOR = 1;
     static const unsigned int GMS_RANGE_ESTIMATION_2D_RMS_MICRO = 0;
     static const unsigned int GMS_RANGE_ESTIMATION_2D_RMS_FULLVER =
       1000U*GMS_RANGE_ESTIMATION_2D_RMS_MAJOR+100U*GMS_RANGE_ESTIMATION_2D_RMS_MINOR+
       10U*GMS_RANGE_ESTIMATION_2D_RMS_MICRO;
     static const char GMS_RANGE_ESTIMATION_2D_RMS_CREATION_DATE[] = "06-10-2026 05:55PM +00200 (TUE 06 OCT 2026 GMT+2)";
     static const char GMS_RANGE_ESTIMATION_2D_RMS_BUILD_DATE[]    = __DATE__; 
     static const char GMS_RANGE_ESTIMATION_2D_RMS_BUILD_TIME[]    = __TIME__;
     static const char GMS_RANGE_ESTIMATION_2D_RMS_AUTHOR[]        = "Programmer: Bernard Gingold, beniekg@gmail.com";
     static const char GMS_RANGE_ESTIMATION_2D_RMS_SYNOPSIS[]      = "Based upon the chapter 2 of the book titled (Rus): Статистическая теория радионавигации Ярлыков, Михаил Семенович,\
                                                                      (Eng):  The Statistical Theory of Radio-Navigation, authorship of Mikhail Semyonovich Yarlykov";
}

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(RANGE_ESTIMATION_2D_RMS_USE_CEPHES)
#define RANGE_ESTIMATION_2D_RMS_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (RANGE_ESTIMATION_2D_RMS_USE_PMC_INSTRUMENTATION)
#define RANGE_ESTIMATION_2D_RMS_USE_PMC_INSTRUMENTATION 0
#endif 

#if (RANGE_ESTIMATION_2D_RMS_USE_PMC_INSTRUMENTATION) == 1
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

#if (RANGE_ESTIMATION_2D_RMS_USE_CEPHES) == 0
#include <cmath>
#else 
#include "GMS_cephes_sin_cos.h"
#endif 

namespace gms 
{

namespace radio_navigation
{

struct alignas(64) range_estimation_2D_rms_t 
{
    std::size_t  m_nsteps; // in
    darray_r4_t  m_radio_range_buoy_a; // 'x' in
    darray_r4_t  m_radio_range_buoy_b; //  'y' in
    darray_r4_t  m_radio_buoys_target_angle; // in
    darray_r4_t  m_rmse_signal_delay; // in
    darray_r4_t  m_hyperbolic_pose_shift; // out
    darray_r4_t  m_hyperbolic_pose_rms; // out: the result of calculation
    float        m_radio_buoys_ab_base; // in

    range_estimation_2D_rms_t() = default;

    range_estimation_2D_rms_t(const std::size_t,
                              const darray_r4_t &,
                              const darray_r4_t &,
                              const darray_r4_t &,
                              const darray_r4_t &,
                              const float) noexcept(false);

    range_estimation_2D_rms_t(const std::size_t,
                              darray_r4_t &&,
                              darray_r4_t &&,
                              darray_r4_t &&,
                              darray_r4_t &&,
                              const float) noexcept(true);

    range_estimation_2D_rms_t(range_estimation_2D_rms_t &&) noexcept(true);

    range_estimation_2D_rms_t(const range_estimation_2D_rms_t &) noexcept(false);

    ~range_estimation_2D_rms_t() noexcept(true);

    range_estimation_2D_rms_t & operator=(const range_estimation_2D_rms_t &) noexcept(false);

    range_estimation_2D_rms_t & operator=(range_estimation_2D_rms_t &&) noexcept(true);

    bool analytic_grad_D_param(const darray_r4_t &,
                               const darray_r4_t &,
                               darray_r4_t &);

    bool compute_hyberbolic_pose_rms( const darray_r4_t &,
                                      const darray_r4_t &,
                                      const darray_r4_t &,
                                      std::optional<float>);    
    
    __ATTR_ALWAYS_INLINE__
    inline 
    bool reset_hyperbolic_pose_shift(const float fill) noexcept
    {
        bool b_state = false;
        std::fill(this->m_hyperbolic_pose_shift.m_data,this->m_hyperbolic_pose_shift+this->m_nsteps,fill);
        b_state = true;
        return (b_state);
    }

    __ATTR_ALWAYS_INLINE__
    inline
    bool reset_hyperbolic_pose_rms(const float fill) noexcept 
    {
        bool b_state = false;
        std::fill(this->m_hyperbolic_pose_rms.m_data,this->m_hyperbolic_pose_rms.m_data+this->m_nsteps,fill);
        b_state = true;
        return (b_state);
    }

};

    auto operator<<(std::ostream &,
                    range_estimation_2D_rms_t &)->std::ostream &;

    auto operator==(const range_estimation_2D_rms_t &,
                    const range_estimation_2D_rms_t &)->std::vector<bool>;

    auto operator!=(const range_estimation_2D_rms_t &,
                    const range_estimation_2D_rms_t &)->std::vector<bool>;

}

}

#endif /*__GMS_RANGE_ESTIMATION_RMS_H__*/