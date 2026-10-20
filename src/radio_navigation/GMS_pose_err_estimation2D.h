
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

#ifndef __GMS_POSE_ERR_ESTIMATION_2D_H__
#define __GMS_POSE_ERR_ESTIMATION_2D_H__ 201020260822

#include <cstdint>
#include <iostream>
#include <vector>
#include <optional>
#include "GMS_config.h"
#include "GMS_dyn_array.h"

namespace file_info 
{

     static const unsigned int GMS_POSE_ERR_ESTIMATION_2D_MAJOR = 1;
     static const unsigned int GMS_POSE_ERR_ESTIMATION_2D_MINOR = 1;
     static const unsigned int GMS_POSE_ERR_ESTIMATION_2D_MICRO = 0;
     static const unsigned int GMS_POSE_ERR_ESTIMATION_2D_FULLVER =
       1000U*GMS_POSE_ERR_ESTIMATION_2D_MAJOR+100U*GMS_POSE_ERR_ESTIMATION_2D_MINOR+
       10U*GMS_POSE_ERR_ESTIMATION_2D_MICRO;
     static const char GMS_POSE_ERR_ESTIMATION_2D_CREATION_DATE[]  = "20-10-2026 08:22AM +00200 (TUE 20 OCT 2026 GMT+2)";
     static const char GMS_POSE_ERR_ESTIMATION_2D_BUILD_DATE[]    = __DATE__; 
     static const char GMS_POSE_ERR_ESTIMATION_2D_BUILD_TIME[]    = __TIME__;
     static const char GMS_POSE_ERR_ESTIMATION_2D_AUTHOR[]        = "Programmer: Bernard Gingold, beniekg@gmail.com";
     static const char GMS_POSE_ERR_ESTIMATION_2D_SYNOPSIS[]      = "Based upon the chapter 2.2.2 of the book titled (Rus): Статистическая теория радионавигации Ярлыков, Михаил Семенович,\
                                                                      (Eng):  The Statistical Theory of Radio-Navigation, authorship of Mikhail Semyonovich Yarlykov";
}

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(POSE_ERR_ESTIMATION_2D_USE_CEPHES)
#define POSE_ERR_ESTIMATION_2D_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (POSE_ERR_ESTIMATION_2D_USE_PMC_INSTRUMENTATION)
#define POSE_ERR_ESTIMATION_2D_USE_PMC_INSTRUMENTATION 0
#endif 

#if (POSE_ERR_ESTIMATION_2D_USE_PMC_INSTRUMENTATION) == 1
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

#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 0
#include <cmath>
#else 
#include "GMS_cephes_sin_cos.h"
#endif 

namespace gms
{

namespace radio_navigation
{

alignas(64) struct pose_err_estimation_2D_compare_t
{
        std::size_t m_nvals;
        std::vector<bool> m_r1_state_comp;
        std::vector<bool> m_r2_state_comp;
        std::vector<bool> m_l1_est_pdf_state_comp;
        std::vector<bool> m_l2_est_pdf_state_comp;
        std::vector<bool> m_l12_est_joint_pdf_state_comp;
        std::vector<bool> m_pose_prob_state_comp;

        pose_err_estimation_2D_compare_t() noexcept(true) = default;

        pose_err_estimation_2D_compare_t(const std::size_t nvals) 
        :
        m_nvals(nvals),
        m_r1_state_comp(m_nvals),
        m_r2_state_comp(m_nvals),
        m_l1_est_pdf_state_comp(m_nvals),
        m_l2_est_pdf_state_comp(m_nvals),
        m_l12_est_joint_pdf_state_comp(m_nvals),
        m_pose_prob_state_comp(m_nvals)
        {

        }

        pose_err_estimation_2D_compare_t(const pose_err_estimation_2D_compare_t &other) noexcept(false) 
        :
        m_nvals(other.m_nvals),
        m_r1_state_comp(other.m_r1_state_comp),
        m_r2_state_comp(other.m_r2_state_comp),
        m_l1_est_pdf_state_comp(other.m_l1_est_pdf_state_comp),
        m_l2_est_pdf_state_comp(other.m_l2_est_pdf_state_comp),
        m_l12_est_joint_pdf_state_comp(other.m_l12_est_joint_pdf_state_comp),
        m_pose_prob_state_comp(other.m_pose_prob_state_comp)
        {

        }

        pose_err_estimation_2D_compare_t(pose_err_estimation_2D_compare_t &&other) noexcept(true)
        :
        m_nvals(other.m_nvals),
        m_r1_state_comp(std::move(other.m_r1_state_comp)),
        m_r2_state_comp(std::move(other.m_r2_state_comp)),
        m_l1_est_pdf_state_comp(std::move(other.m_l1_est_pdf_state_comp)),
        m_l2_est_pdf_state_comp(std::move(other.m_l2_est_pdf_state_comp)),
        m_l12_est_joint_pdf_state_comp(std::move(other.m_l12_est_joint_pdf_state_comp)),
        m_pose_prob_state_comp(std::move(other.m_pose_prob_state_comp))
        {

        }

        ~pose_err_estimation_2D_compare_t() = default;

        pose_err_estimation_2D_compare_t & 
        operator=(const pose_err_estimation_2D_compare_t &rhs) noexcept(false)
        {
            if(__builtin_expect(this==&rhs,0)) { return (*this);}
            
            this->m_nvals = other.m_nvals;
            this->m_r1_state_comp.operator=(rhs.m_r1_state_comp);
            this->m_r2_state_comp.operator=(rhs.m_r2_state_comp);
            this->m_l1_est_pdf_state_comp.operator=(rhs.m_l1_est_pdf_state_comp);
            this->m_l2_est_pdf_state_comp.operator=(rhs.m_l2_est_pdf_state_comp);
            this->m_l12_est_joint_pdf_state_comp.operator=(rhs.m_l12_est_joint_pdf_state_comp);
            this->m_pose_prob_state_comp.operator=(rhs.m_pose_prob_state_comp);

            return (*this);
        }

        pose_err_estimation_2D_compare_t & 
        operator=(pose_err_estimation_2D_compare_t &&rhs) noexcept(true)
        {
            if(__builtin_expect(this==&rhs,0)) { return (*this);}
            
            this->m_nvals = other.m_nvals;
            this->m_r1_state_comp.operator=(rhs.m_r1_state_comp);
            this->m_r2_state_comp.operator=(rhs.m_r2_state_comp);
            this->m_l1_est_pdf_state_comp.operator=(rhs.m_l1_est_pdf_state_comp);
            this->m_l2_est_pdf_state_comp.operator=(rhs.m_l2_est_pdf_state_comp);
            this->m_l12_est_joint_pdf_state_comp.operator=(rhs.m_l12_est_joint_pdf_state_comp);
            this->m_pose_prob_state_comp.operator=(rhs.m_pose_prob_state_comp);

            return (*this);
        }

};

/* Valid only for the Gaussian PDF prior assumption!!*/
alignas(64) struct pose_err_estimation_2D_t
{
    std::size_t    m_nsteps;
    darray_r4_t    m_r1; // true pose radius       (output-state)
    darray_r4_t    m_r2; // errornous pose radius  (output-state)
    darray_r4_t    m_l1; // AB-line estimation error (input)
    darray_r4_t    m_l2; // CD-line estimation error (input)
    darray_r4_t    m_alpham; // a variable angle between the line CD and line AB (input)
    darray_r4_t    m_sigma_l1; //RMS error estimation of AB (input)
    darray_r4_t    m_sigma_l2; //RMS error estimation of CD (input)
    darray_r4_t    m_rho; // cross-correl coefficient 'rho' of mutually measured lines AB and CD (input) 
    darray_r4_t    m_l1_est_pdf; // Normal PDF of l1 estimation (output-state)
    darray_r4_t    m_l2_est_pdf; // Normal PDF of l2 estimation (output-state)
    darray_r4_t    m_l12_joint_est_pdf; // Joint Normal PDF of l1,l2 error estimation (output-state)
    darray_r4_t    m_pose_error_prob; // Gaussian PDF INtegrated (CDF) of position error ellipsoid (output-state)

   
    
    pose_err_estimation_2D_t() noexcept(true) = default;

    pose_err_estimation_2D_t(const std::size_t,
                             const darray_r4_t &,
                             const darray_r4_t &,
                             const darray_r4_t &,
                             const darray_r4_t &,
                             const darray_r4_t &,
                             const darray_r4_t &) noexcept(false);

    pose_err_estimation_2D_t(const std::size_t,
                             darray_r4_t &&,
                             darray_r4_t &&,
                             darray_r4_t &&,
                             darray_r4_t &&,
                             darray_r4_t &&,
                             darray_r4_t &&) noexcept(true);

    pose_err_estimation_2D_t(const pose_err_estimation_2D_t &) noexcept(false);

    pose_err_estimation_2D_t(pose_err_estimation_2D_t &&) noexcept(true);

    ~pose_err_estimation_2D_t();

    pose_err_estimation_2D_t & operator=(const pose_err_estimation_2D_t &) noexcept(false);

    pose_err_estimation_2D_t & operator=(pose_err_estimation_2D_t &&) noexcept(true);

    __ATTR_ALWAYS_INLINE__
    inline 
    bool reset_r1(const float fill) noexcept(true)
    {
        bool ret_state = false;
        std::fill(this->m_r1.m_data,this->m_r1.m_data+this->m_nsteps,fill);
        ret_state = true;
        return (ret_state);
    }

    __ATTR_ALWAYS_INLINE__
    inline 
    bool reset_r2(const float fill) noexcept(true)
    {
        bool ret_state = false;
        std::fill(this->m_r2.m_data,this->m_r2.m_data+this->m_nsteps,fill);
        ret_state = true;
        return (ret_state);
    }

    __ATTR_ALWAYS_INLINE__
    inline 
    bool reset_l1_est_pdf(const float fill) noexcept(true)
    {
        bool ret_state = false;
        std::fill(this->m_l1_est_pdf.m_data,this->m_l1_est_pdf.m_data+this->m_nsteps,fill);
        ret_state = true;
        return (ret_state);
    }

    __ATTR_ALWAYS_INLINE__
    inline 
    bool reset_l2_est_pdf(const float fill) noexcept(true)
    {
        bool ret_state = false;
        std::fill(this->m_l2_est_pdf.m_data,this->m_l2_est_pdf.m_data+this->m_nsteps,fill);
        ret_state = true;
        return (ret_state);
    }
    
    __ATTR_ALWAYS_INLINE__
    inline 
    bool reset_l12_joint_est_pdf(const float fill) noexcept(true)
    {
        bool ret_state = false;
        std::fill(this->m_l12_joint_est_pdf.m_data,this->m_l12_joint_est_pdf.m_data+this->m_nsteps,fill);
        ret_state = true;
        return (ret_state);
    }

    __ATTR_ALWAYS_INLINE__
    inline 
    bool reset_pose_err_prob(const float fill) noexcept(true)
    {
        bool ret_state = false;
        std::fill(this->m_pose_error_prob.m_data,this->m_pose_error_prob.m_data+this->m_nsteps,fill);
        ret_state = true;
        return (ret_state);
    }

    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    std::int32_t compute_l1_pdf();

    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    std::int32_t compute_l2_pdf();

    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    std::int32_t compute_l12_joint_pdf();

    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    std::int32_t compue_pose_error_prob();

};

auto operator<<(std::ostream &,const pose_err_estimation_2D_t &)->std::ostream &;

auto operator==(const pose_err_estimation_2D_t &,
                const pose_err_estimation_2D_t &)->pose_err_estimation_2D_compare_t;

auto operator!=(const pose_err_estimation_2D_t &,
                const pose_err_estimation_2D_t &)->pose_err_estimation_2D_compare_t;

}

}

#endif /*__GMS_POSE_ERR_ESTIMATION_2D_H__*/