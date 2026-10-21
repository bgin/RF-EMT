
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

#include "GMS_pose_err_estimation2D.h"
#include "GMS_compare_fp_safe.h"


gms::radio_navigation
::pose_err_estimation_2D_t
::pose_err_estimation_2D_t(const std::size_t nsteps,
                           const darray_r4_t &l1,
                           const darray_r4_t &l2,
                           const darray_r4_t &alpham,
                           const darray_r4_t &sigma_l1,
                           const darray_r4_t &sigma_l2,
                           const darray_r4_t &rho)
:
m_nsteps(nsteps),
m_r1(m_nsteps),
m_r2(m_nsteps),
m_l1(l1),
m_l2(l2),
m_alpham(alpham),
m_sigma_l1(sigma_l1),
m_sigma_l2(sigma_l2),
m_rho(rho),
m_l1_est_pdf(m_nsteps),
m_l2_est_pdf(m_nsteps),
m_l12_joint_est_pdf(m_nsteps),
m_pose_error_prob(m_nsteps)
{

}

gms::radio_navigation
::pose_err_estimation_2D_t
::pose_err_estimation_2D_t(const std::size_t nsteps,
                           darray_r4_t &&l1,
                           darray_r4_t &&l2,
                           darray_r4_t &&alpham,
                           darray_r4_t &&sigma_l1,
                           darray_r4_t &&sigma_l2,
                           darray_r4_t &&rho)
:
m_nsteps(nsteps),
m_r1(m_nsteps),
m_r2(m_nsteps),
m_l1(std::move(l1)),
m_l2(std::move(l2)),
m_alpham(std::move(alpham)),
m_sigma_l1(std::move(sigma_l1)),
m_sigma_l2(std::move(sigma_l2)),
m_rho(std::move(rho)),
m_l1_est_pdf(m_nsteps),
m_l2_est_pdf(m_nsteps),
m_l12_joint_est_pdf(m_nsteps),
m_pose_error_prob(m_nsteps)
{

}

gms::radio_navigation
::pose_err_estimation_2D_t
::pose_err_estimation_2D_t(const pose_err_estimation_2D_t &other) noexcept(false)
:
m_nsteps(other.m_nsteps),
m_r1(other.m_r1),
m_r2(other.m_r2),
m_l1(other.m_l1),
m_l2(other.m_l2),
m_alpham(other.m_alpham),
m_sigma_l1(other.m_sigma_l1),
m_sigma_l2(other.m_sigma_l2),
m_rho(other.m_rho),
m_l1_est_pdf(other.m_l1_est_pdf),
m_l2_est_pdf(other.m_l2_est_pdf),
m_l12_joint_est_pdf(other.m_l12_joint_est_pdf),
m_pose_error_prob(other.m_pose_error_prob)
{

}

gms::radio_navigation
::pose_err_estimation_2D_t
::pose_err_estimation_2D_t(pose_err_estimation_2D_t &&other) noexcept(true)
:
m_nsteps(other.m_nsteps),
m_r1(std::move(other.m_r1)),
m_r2(std::move(other.m_r2)),
m_l1(std::move(other.m_l1)),
m_l2(std::move(other.m_l2)),
m_alpham(std::move(other.m_alpham)),
m_sigma_l1(std::move(other.m_sigma_l1)),
m_sigma_l2(std::move(other.m_sigma_l2)),
m_rho(std::move(other.m_rho)),
m_l1_est_pdf(std::move(other.m_l1_est_pdf)),
m_l2_est_pdf(std::move(other.m_l2_est_pdf)),
m_l12_joint_est_pdf(std::move(other.m_l12_joint_est_pdf)),
m_pose_error_prob(std::move(other.m_pose_error_prob))
{

}

gms::radio_navigation
::pose_err_estimation_2D_t
::~pose_err_estimation_2D_t()
{

}

gms::radio_navigation
::pose_err_estimation_2D_t &
gms::radio_navigation
::pose_err_estimation_2D_t
::operator=(const pose_err_estimation_2D_t &other) noexcept(false)
{
    if(__builtin_expect(this==&other,0)) { return (*this);}

    this->m_nsteps = other.m_nsteps;
    this->m_r1.operator=(other.m_r1);
    this->m_r2.operator=(other.m_r2);
    this->m_l1.operator=(other.m_l1);
    this->m_l2.operator=(other.m_l2);
    this->m_alpham.operator=(other.m_alpham);
    this->m_sigma_l1.operator=(other.m_sigma_l1);
    this->m_sigma_l2.operator=(other.m_sigma_l2);
    this->m_rho.operator=(other.m_rho);
    this->m_l1_est_pdf.operator=(other.m_l1_est_pdf);
    this->m_l2_est_pdf.operator=(other.m_l2_est_pdf);
    this->m_l12_joint_est_pdf.operator=(other.m_l12_joint_est_pdf);
    this->m_pose_error_prob.operator=(other.m_pose_error_prob);

    return (*this);
}

gms::radio_navigation
::pose_err_estimation_2D_t &
gms::radio_navigation
::pose_err_estimation_2D_t
::operator=(pose_err_estimation_2D_t &&other) noexcept(false)
{
    if(__builtin_expect(this==&other,0)) { return (*this);}

    this->m_nsteps = other.m_nsteps;
    this->m_r1.operator=(std::move(other.m_r1));
    this->m_r2.operator=(std::move(other.m_r2));
    this->m_l1.operator=(std::move(other.m_l1));
    this->m_l2.operator=(std::move(other.m_l2));
    this->m_alpham.operator=(std::move(other.m_alpham));
    this->m_sigma_l1.operator=(std::move(other.m_sigma_l1));
    this->m_sigma_l2.operator=(std::move(other.m_sigma_l2));
    this->m_rho.operator=(std::move(other.m_rho));
    this->m_l1_est_pdf.operator=(std::move(other.m_l1_est_pdf));
    this->m_l2_est_pdf.operator=(std::move(other.m_l2_est_pdf));
    this->m_l12_joint_est_pdf.operator=(std::move(other.m_l12_joint_est_pdf));
    this->m_pose_error_prob.operator=(std::move(other.m_pose_error_prob));

    return (*this);
}

std::int32_t 
gms::radio_navigation
::pose_err_estimation_2D_t
::compute_l1_pdf()
{
    constexpr float sqrt2Pi = 2.506628274631000502415765284811;
    std::size_t i,j;
    float val_l1;
    float val_l1s;
    float val_sigl1;
    float val_sigl1s;
    float fac1;
    float exp_arg;
    float exp_val;
    float pdf_res;
    std::int32_t ret_state = -9999;
    if(__builtin_expect(this->m_nsteps>=1ull,0) && 
       __builtin_expect(this->m_nsteps<=16ull,0)) 
    {
       for(i = 0ull; i<this->m_nsteps; ++i) 
       {
           val_l1    = this->m_l1.m_data[i];
           val_l1s   = val_l1*val_l1;
           val_sigl1 = this->m_sigma_l1.m_data[i];
           val_sigl1s= val_sigl1*val_sigl1;
           fac1      = 1.0f/(val_sigl1*sqrt2Pi);
           exp_arg   = -val_l1s/(val_sigl1s+val_sigl1s);
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
           exp_val   = gms::math::ceph_expf(exp_arg);
#else 
           exp_val   = std::exp(exp_arg);
#endif 
           pdf_res   = fac1*exp_val;
           this->m_l1_est_pdf.m_data[i];
       }
       ret_state = 1;
       return (ret_state);
    }
    else if(__builtin_expect(this->m_nsteps>16ull,1)) 
    {
        for(i = 0; i!=ROUND_DOWN(this->m_nsteps,3ull); i += 4ull) 
        {
            const std::size_t idx1 = i+1ull;
            const std::size_t idx2 = i+2ull;
            const std::size_t idx3 = i+3ull;
            val_l1    = this->m_l1.m_data[i];
            val_l1s   = val_l1*val_l1;
            val_sigl1 = this->m_sigma_l1.m_data[i];
            val_sigl1s= val_sigl1*val_sigl1;
            fac1      = 1.0f/(val_sigl1*sqrt2Pi);
            exp_arg   = -val_l1s/(val_sigl1s+val_sigl1s);
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val   = gms::math::ceph_expf(exp_arg);
#else 
            exp_val   = std::exp(exp_arg);
#endif 
            pdf_res   = fac1*exp_val;
            this->m_l1_est_pdf.m_data[i];

            val_l1    = this->m_l1.m_data[idx1];
            val_l1s   = val_l1*val_l1;
            val_sigl1 = this->m_sigma_l1.m_data[idx1];
            val_sigl1s= val_sigl1*val_sigl1;
            fac1      = 1.0f/(val_sigl1*sqrt2Pi);
            exp_arg   = -val_l1s/(val_sigl1s+val_sigl1s);
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val   = gms::math::ceph_expf(exp_arg);
#else 
            exp_val   = std::exp(exp_arg);
#endif 
            pdf_res   = fac1*exp_val;
            this->m_l1_est_pdf.m_data[idx1];

            val_l1    = this->m_l1.m_data[idx2];
            val_l1s   = val_l1*val_l1;
            val_sigl1 = this->m_sigma_l1.m_data[idx2];
            val_sigl1s= val_sigl1*val_sigl1;
            fac1      = 1.0f/(val_sigl1*sqrt2Pi);
            exp_arg   = -val_l1s/(val_sigl1s+val_sigl1s);
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val   = gms::math::ceph_expf(exp_arg);
#else 
            exp_val   = std::exp(exp_arg);
#endif 
            pdf_res   = fac1*exp_val;
            this->m_l1_est_pdf.m_data[idx2];

            val_l1    = this->m_l1.m_data[idx3];
            val_l1s   = val_l1*val_l1;
            val_sigl1 = this->m_sigma_l1.m_data[idx3];
            val_sigl1s= val_sigl1*val_sigl1;
            fac1      = 1.0f/(val_sigl1*sqrt2Pi);
            exp_arg   = -val_l1s/(val_sigl1s+val_sigl1s);
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val   = gms::math::ceph_expf(exp_arg);
#else 
            exp_val   = std::exp(exp_arg);
#endif 
            pdf_res   = fac1*exp_val;
            this->m_l1_est_pdf.m_data[idx3];
        }

        for(j = i; j<this->m_nsteps; ++j)  
        {
            val_l1    = this->m_l1.m_data[j];
            val_l1s   = val_l1*val_l1;
            val_sigl1 = this->m_sigma_l1.m_data[j];
            val_sigl1s= val_sigl1*val_sigl1;
            fac1      = 1.0f/(val_sigl1*sqrt2Pi);
            exp_arg   = -val_l1s/(val_sigl1s+val_sigl1s);
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val   = gms::math::ceph_expf(exp_arg);
#else 
            exp_val   = std::exp(exp_arg);
#endif 
            pdf_res   = fac1*exp_val;
            this->m_l1_est_pdf.m_data[j];
        }
        ret_state = 2;
        return (ret_state);
    }
}

std::int32_t 
gms::radio_navigation
::pose_err_estimation_2D_t
::compute_l2_pdf()
{
    constexpr float sqrt2Pi = 2.506628274631000502415765284811;
    std::size_t i,j;
    float val_l2;
    float val_l2s;
    float val_sigl2;
    float val_sigl2s;
    float fac2;
    float exp_arg;
    float exp_val;
    float pdf_res;
    std::int32_t ret_state = -9999;
    if(__builtin_expect(this->m_nsteps>=1ull,0) && 
       __builtin_expect(this->m_nsteps<=16ull,0)) 
    {
       for(i = 0ull; i<this->m_nsteps; ++i) 
       {
           val_l2    = this->m_l2.m_data[i];
           val_l2s   = val_l2*val_l2;
           val_sigl2 = this->m_sigma_l2.m_data[i];
           val_sigl2s= val_sigl2*val_sigl2;
           fac2      = 1.0f/(val_sigl2*sqrt2Pi);
           exp_arg   = -val_l2s/(val_sigl2s+val_sigl2s);
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
           exp_val   = gms::math::ceph_expf(exp_arg);
#else 
           exp_val   = std::exp(exp_arg);
#endif 
           pdf_res   = fac2*exp_val;
           this->m_l2_est_pdf.m_data[i];
       }
       ret_state = 1;
       return (ret_state);
    }
    else if(__builtin_expect(this->m_nsteps>16ull,1)) 
    {
        for(i = 0; i!=ROUND_DOWN(this->m_nsteps,3ull); i += 4ull) 
        {
            const std::size_t idx1 = i+1ull;
            const std::size_t idx2 = i+2ull;
            const std::size_t idx3 = i+3ull;
            val_l2    = this->m_l2.m_data[i];
            val_l2s   = val_l2*val_l2;
            val_sigl2 = this->m_sigma_l2.m_data[i];
            val_sigl2s= val_sigl2*val_sigl2;
            fac2      = 1.0f/(val_sigl2*sqrt2Pi);
            exp_arg   = -val_l2s/(val_sigl2s+val_sigl2s);
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val   = gms::math::ceph_expf(exp_arg);
#else 
            exp_val   = std::exp(exp_arg);
#endif 
            pdf_res   = fac2*exp_val;
            this->m_l2_est_pdf.m_data[i];

            val_l2    = this->m_l2.m_data[idx1];
            val_l2s   = val_l2*val_l2;
            val_sigl2 = this->m_sigma_l2.m_data[idx1];
            val_sigl2s= val_sigl2*val_sigl2;
            fac2      = 1.0f/(val_sigl2*sqrt2Pi);
            exp_arg   = -val_l2s/(val_sigl2s+val_sigl2s);
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val   = gms::math::ceph_expf(exp_arg);
#else 
            exp_val   = std::exp(exp_arg);
#endif 
            pdf_res   = fac2*exp_val;
            this->m_l2_est_pdf.m_data[idx1];

            val_l2    = this->m_l2.m_data[idx2];
            val_l2s   = val_l2*val_l2;
            val_sigl2 = this->m_sigma_l2.m_data[idx2];
            val_sigl2s= val_sigl2*val_sigl2;
            fac2      = 1.0f/(val_sigl2*sqrt2Pi);
            exp_arg   = -val_l2s/(val_sigl2s+val_sigl2s);
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val   = gms::math::ceph_expf(exp_arg);
#else 
            exp_val   = std::exp(exp_arg);
#endif 
            pdf_res   = fac2*exp_val;
            this->m_l2_est_pdf.m_data[idx2];

            val_l2    = this->m_l2.m_data[idx3];
            val_l2s   = val_l2*val_l2;
            val_sigl2 = this->m_sigma_l2.m_data[idx3];
            val_sigl2s= val_sigl2*val_sigl2;
            fac2      = 1.0f/(val_sigl2*sqrt2Pi);
            exp_arg   = -val_l2s/(val_sigl2s+val_sigl2s);
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val   = gms::math::ceph_expf(exp_arg);
#else 
            exp_val   = std::exp(exp_arg);
#endif 
            pdf_res   = fac2*exp_val;
            this->m_l2_est_pdf.m_data[idx3];
        }

        for(j = i; j<this->m_nsteps; ++j)  
        {
            val_l2    = this->m_l2.m_data[j];
            val_l2s   = val_l2*val_l2;
            val_sigl2 = this->m_sigma_l2.m_data[j];
            val_sigl2s= val_sigl2*val_sigl2;
            fac2      = 1.0f/(val_sigl2*sqrt2Pi);
            exp_arg   = -val_l2s/(val_sigl2s+val_sigl2s);
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val   = gms::math::ceph_expf(exp_arg);
#else 
            exp_val   = std::exp(exp_arg);
#endif 
            pdf_res   = fac2*exp_val;
            this->m_l2_est_pdf.m_data[j];
        }
        ret_state = 2;
        return (ret_state);
    }
}

std::int32_t 
gms::radio_navigation
::pose_err_estimation_2D_t
::compute_l12_joint_pdf()
{
    constexpr float twoPi = 6.283185307179586476925286766559f;
    std::size_t i,j;
    float val_l1;
    float val_l2;
    float val_l1s;
    float val_l2s;
    float val_sigl1;
    float val_sigl2;
    float val_sigl1s;
    float val_sigl2s;
    float fac;
    float left_rat;
    float right_rat;
    float exp_arg;
    float exp_val;
    float joint_pdf;
    std::int32_t ret_state = -9999;
    if(__builtin_expect(this->m_nsteps>=1ull,0) && 
       __builtin_expect(this->m_nsteps<=16ull,0))
    {
        for(i = 0ull; i<this->m_nsteps; ++i) 
        {
            val_l1     = this->m_l1.m_data[i];
            val_l1s    = val_l1*val_l1;
            val_l2     = this->m_l2.m_data[i];
            val_l2s    = val_l2*val_l2;
            val_sigl1  = this->m_sigma_l1.m_data[i];
            val_sigl1s = val_sigl1*val_sigl1;
            val_sigl2  = this->m_sigma_l2.m_data[i];
            val_sigl2s = val_sigl2*val_sigl2;
            fac        = 1.0f/(twoPi*val_sigl1*val_sigl2);
            left_rat   = val_l1s/val_sigl1s;
            right_rat  = val_l2s/val_sigl2s;
            exp_arg    = -0.5f*(left_rat+right_rat);
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val    = gms::math::ceph_expf(exp_arg);
#else 
            exp_val    = std::exp(exp_arg);
#endif
            joint_pdf  = fac*exp_val;
            this->m_l12_joint_est_pdf.m_data[i] = joint_pdf;
        }
        ret_state = 1;
        return (ret_state);
    }
    else if(__builtin_expect(this->m_nsteps>16ull,1))
    {
        for(i = 0ull; i!=ROUND_DOWN(this->m_nsteps,3ull); i += 4ull) 
        {
            const std::size_t idx1 = i+1ull;
            const std::size_t idx2 = i+2ull;
            const std::size_t idx3 = i+3ull;
            val_l1     = this->m_l1.m_data[i];
            val_l1s    = val_l1*val_l1;
            val_l2     = this->m_l2.m_data[i];
            val_l2s    = val_l2*val_l2;
            val_sigl1  = this->m_sigma_l1.m_data[i];
            val_sigl1s = val_sigl1*val_sigl1;
            val_sigl2  = this->m_sigma_l2.m_data[i];
            val_sigl2s = val_sigl2*val_sigl2;
            fac        = 1.0f/(twoPi*val_sigl1*val_sigl2);
            left_rat   = val_l1s/val_sigl1s;
            right_rat  = val_l2s/val_sigl2s;
            exp_arg    = -0.5f*(left_rat+right_rat);
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val    = gms::math::ceph_expf(exp_arg);
#else 
            exp_val    = std::exp(exp_arg);
#endif
            joint_pdf  = fac*exp_val;
            this->m_l12_joint_est_pdf.m_data[i] = joint_pdf;

            val_l1     = this->m_l1.m_data[idx1];
            val_l1s    = val_l1*val_l1;
            val_l2     = this->m_l2.m_data[idx1];
            val_l2s    = val_l2*val_l2;
            val_sigl1  = this->m_sigma_l1.m_data[idx1];
            val_sigl1s = val_sigl1*val_sigl1;
            val_sigl2  = this->m_sigma_l2.m_data[idx1];
            val_sigl2s = val_sigl2*val_sigl2;
            fac        = 1.0f/(twoPi*val_sigl1*val_sigl2);
            left_rat   = val_l1s/val_sigl1s;
            right_rat  = val_l2s/val_sigl2s;
            exp_arg    = -0.5f*(left_rat+right_rat);
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val    = gms::math::ceph_expf(exp_arg);
#else 
            exp_val    = std::exp(exp_arg);
#endif
            joint_pdf  = fac*exp_val;
            this->m_l12_joint_est_pdf.m_data[idx1] = joint_pdf;

            val_l1     = this->m_l1.m_data[idx2];
            val_l1s    = val_l1*val_l1;
            val_l2     = this->m_l2.m_data[idx2];
            val_l2s    = val_l2*val_l2;
            val_sigl1  = this->m_sigma_l1.m_data[idx2];
            val_sigl1s = val_sigl1*val_sigl1;
            val_sigl2  = this->m_sigma_l2.m_data[idx2];
            val_sigl2s = val_sigl2*val_sigl2;
            fac        = 1.0f/(twoPi*val_sigl1*val_sigl2);
            left_rat   = val_l1s/val_sigl1s;
            right_rat  = val_l2s/val_sigl2s;
            exp_arg    = -0.5f*(left_rat+right_rat);
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val    = gms::math::ceph_expf(exp_arg);
#else 
            exp_val    = std::exp(exp_arg);
#endif
            joint_pdf  = fac*exp_val;
            this->m_l12_joint_est_pdf.m_data[idx2] = joint_pdf;

            val_l1     = this->m_l1.m_data[idx3];
            val_l1s    = val_l1*val_l1;
            val_l2     = this->m_l2.m_data[idx3];
            val_l2s    = val_l2*val_l2;
            val_sigl1  = this->m_sigma_l1.m_data[idx3];
            val_sigl1s = val_sigl1*val_sigl1;
            val_sigl2  = this->m_sigma_l2.m_data[idx3];
            val_sigl2s = val_sigl2*val_sigl2;
            fac        = 1.0f/(twoPi*val_sigl1*val_sigl2);
            left_rat   = val_l1s/val_sigl1s;
            right_rat  = val_l2s/val_sigl2s;
            exp_arg    = -0.5f*(left_rat+right_rat);
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val    = gms::math::ceph_expf(exp_arg);
#else 
            exp_val    = std::exp(exp_arg);
#endif
            joint_pdf  = fac*exp_val;
            this->m_l12_joint_est_pdf.m_data[idx3] = joint_pdf;
        }

        for(j = i; j<this->m_nsteps; ++j)  
        {
            val_l1     = this->m_l1.m_data[j];
            val_l1s    = val_l1*val_l1;
            val_l2     = this->m_l2.m_data[j];
            val_l2s    = val_l2*val_l2;
            val_sigl1  = this->m_sigma_l1.m_data[j];
            val_sigl1s = val_sigl1*val_sigl1;
            val_sigl2  = this->m_sigma_l2.m_data[j];
            val_sigl2s = val_sigl2*val_sigl2;
            fac        = 1.0f/(twoPi*val_sigl1*val_sigl2);
            left_rat   = val_l1s/val_sigl1s;
            right_rat  = val_l2s/val_sigl2s;
            exp_arg    = -0.5f*(left_rat+right_rat);
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val    = gms::math::ceph_expf(exp_arg);
#else 
            exp_val    = std::exp(exp_arg);
#endif
            joint_pdf  = fac*exp_val;
            this->m_l12_joint_est_pdf.m_data[j] = joint_pdf;
        }
        ret_state = 2;
        return (ret_state);
    }
}

std::int32_t 
gms::radio_navigation
::pose_err_estimation_2D_t
::compute_pose_error_prob()
{
    std::size_t i,j;
    float val_l1;
    float val_l2;
    float val_l1s;
    float val_l2s;
    float val_sigl1;
    float val_sigl2;
    float val_sigl1s;
    float val_sigl2s;
    float left_rat;
    float right_rat;
    float c;
    float exp_val;
    float cdf;
    std::int32_t ret_state = -9999;
    if(__builtin_expect(this->m_nsteps>=1ull,0) &&
       __builtin_expect(this->m_nsteps<=16ull,0)) 
    {
#if (POSE_ERR_ESTIMATION_2D_VERBOSE) == 1
     [[nodiscard]] std::int32_t printf_ret = std::printf("[***WARNING***]: %s,loc=%d, msg=Integration of the Gaussian PDF!!\n",__func__,__LINE__);
#endif
        for(i = 0ull; i<this->m_nsteps; ++i) 
        {
            val_l1     = this->m_l1.m_data[i];
            val_l1s    = val_l1*val_l1;
            val_l2     = this->m_l2.m_data[i];
            val_l2s    = val_l2*val_l2;
            val_sigl1  = this->m_sigma_l1.m_data[i];
            val_sigl1s = val_sigl1*val_sigl1;
            val_sigl2  = this->m_sigma_l2.m_data[i];
            val_sigl2s = val_sigl2*val_sigl2;
            left_rat   = val_l1s/val_sigl1s;
            right_rat  = val_l2s/val_sigl2s;
            c          = left_rat+right_rat;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val    = gms::math::ceph_expf(-0.5f*c);
#else 
            exp_val    = std::exp(-0.5f*c);
#endif
            cdf        = 1.0f-exp_val;
            this->m_pose_error_prob.m_data[i] = cdf;
        }
        ret_state = 1;
        return (ret_state);
    }
    else if(__builtin_expect(this->m_nsteps>16ull,1))
    {
#if (POSE_ERR_ESTIMATION_2D_VERBOSE) == 1
     [[nodiscard]] std::int32_t printf_ret = std::printf("[***WARNING***]: %s,loc=%d, msg=Integration of the Gaussian PDF!!\n",__func__,__LINE__);
#endif
        for(i = 0ull; i!=ROUND_DOWN(this->m_nsteps,3ull); i += 4ull) 
        {
            const std::size_t idx1 = i+1ull;
            const std::size_t idx2 = i+2ull;
            const std::size_t idx3 = i+3ull;
            val_l1     = this->m_l1.m_data[i];
            val_l1s    = val_l1*val_l1;
            val_l2     = this->m_l2.m_data[i];
            val_l2s    = val_l2*val_l2;
            val_sigl1  = this->m_sigma_l1.m_data[i];
            val_sigl1s = val_sigl1*val_sigl1;
            val_sigl2  = this->m_sigma_l2.m_data[i];
            val_sigl2s = val_sigl2*val_sigl2;
            left_rat   = val_l1s/val_sigl1s;
            right_rat  = val_l2s/val_sigl2s;
            c          = left_rat+right_rat;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val    = gms::math::ceph_expf(-0.5f*c);
#else 
            exp_val    = std::exp(-0.5f*c);
#endif
            cdf        = 1.0f-exp_val;
            this->m_pose_error_prob.m_data[i] = cdf;

            val_l1     = this->m_l1.m_data[idx1];
            val_l1s    = val_l1*val_l1;
            val_l2     = this->m_l2.m_data[idx1];
            val_l2s    = val_l2*val_l2;
            val_sigl1  = this->m_sigma_l1.m_data[idx1];
            val_sigl1s = val_sigl1*val_sigl1;
            val_sigl2  = this->m_sigma_l2.m_data[idx1];
            val_sigl2s = val_sigl2*val_sigl2;
            left_rat   = val_l1s/val_sigl1s;
            right_rat  = val_l2s/val_sigl2s;
            c          = left_rat+right_rat;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val    = gms::math::ceph_expf(-0.5f*c);
#else 
            exp_val    = std::exp(-0.5f*c);
#endif
            cdf        = 1.0f-exp_val;
            this->m_pose_error_prob.m_data[idx1] = cdf;

            val_l1     = this->m_l1.m_data[idx2];
            val_l1s    = val_l1*val_l1;
            val_l2     = this->m_l2.m_data[idx2];
            val_l2s    = val_l2*val_l2;
            val_sigl1  = this->m_sigma_l1.m_data[idx2];
            val_sigl1s = val_sigl1*val_sigl1;
            val_sigl2  = this->m_sigma_l2.m_data[idx2];
            val_sigl2s = val_sigl2*val_sigl2;
            left_rat   = val_l1s/val_sigl1s;
            right_rat  = val_l2s/val_sigl2s;
            c          = left_rat+right_rat;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val    = gms::math::ceph_expf(-0.5f*c);
#else 
            exp_val    = std::exp(-0.5f*c);
#endif
            cdf        = 1.0f-exp_val;
            this->m_pose_error_prob.m_data[idx2] = cdf;

            val_l1     = this->m_l1.m_data[idx3];
            val_l1s    = val_l1*val_l1;
            val_l2     = this->m_l2.m_data[idx3];
            val_l2s    = val_l2*val_l2;
            val_sigl1  = this->m_sigma_l1.m_data[idx3];
            val_sigl1s = val_sigl1*val_sigl1;
            val_sigl2  = this->m_sigma_l2.m_data[idx3];
            val_sigl2s = val_sigl2*val_sigl2;
            left_rat   = val_l1s/val_sigl1s;
            right_rat  = val_l2s/val_sigl2s;
            c          = left_rat+right_rat;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
            exp_val    = gms::math::ceph_expf(-0.5f*c);
#else 
            exp_val    = std::exp(-0.5f*c);
#endif
            cdf        = 1.0f-exp_val;
            this->m_pose_error_prob.m_data[idx3] = cdf;
        }
        ret_state = 2;
        return (ret_state);
    }
}

std::int32_t 
gms::radio_navigation
::pose_err_estimation_2D_t
::compute_radial_pose_err_rms(const bool is_rho_zero)
{
    std::size_t i,j;
    float val_l1;
    float val_l2;
    float val_l1s;
    float val_l2s;
    float val_sigl1;
    float val_sigl2;
    float val_sigl1s;
    float val_sigl2s;
    float val_alpham;
    float cosam;
    float sinam;
    float val_rho;
    float r_rtherm;
    float sig_rtherm;
    float sqrt_r;
    float sqrt_sig;
    std::int32_t ret_state = -9999;
    if(__builtin_expect(this->m_nsteps>=1ull,0) && 
       __builtin_expect(this->m_nsteps<=16ull,0))
    {
        if(is_rho_zero==true)
        {
            for(int i = 0ull; i<this->m_nsteps; ++i)
            {
                val_alpham = this->m_alpham.m_data[i];
                val_l1     = this->m_l1.m_data[i];
                val_l1s    = val_l1*val_l1;
                val_l2     = this->m_l2.m_data[i];
                val_l2s    = val_l2*val_l2;
                val_sigl1  = this->m_sigma_l1.m_data[i];
                val_sigl1s = val_sigl1*val_sigl1;
                val_sigl2  = this->m_sigma_l2.m_data[i];
                val_sigl2s = val_sigl2*val_sigl2;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                cosam      = gms::math::ceph_cosf(val_alpham);
                sinam      = gms::math::ceph_sinf(val_alpham);
#else 
                cosam      = std::cos(val_alpham);
                sinam      = std::sin(val_alpham);
#endif 
                r_rtherm   = 2.0f*val_l1*val_l2*cosam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_sig   = gms::math::ceph_sqrtf(val_sigl1s+val_sigl2s);
#else 
                sqrt_sig   = std::sqrt(val_sigl1s+val_sigl2s);
#endif 
                this->m_r2.m_data[i] = sqrt_sig/sinam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_r     = gms::math::ceph_sqrtf(val_l1s+val_l2s+r_rtherm);
#else 
                sqrt_r     = std::sqrt(val_l1s+val_l2s+r_rtherm);
#endif 
                this->m_r1.m_data[i] = sqrt_r/sinam;
            }
        }
        else 
        {
            for(int i = 0ull; i<this->m_nsteps; ++i)
            {
                val_alpham = this->m_alpham.m_data[i];
                val_l1     = this->m_l1.m_data[i];
                val_rho    = this->m_rho.m_data[i];
                val_l1s    = val_l1*val_l1;
                val_l2     = this->m_l2.m_data[i];
                val_l2s    = val_l2*val_l2;
                val_sigl1  = this->m_sigma_l1.m_data[i];
                val_sigl1s = val_sigl1*val_sigl1;
                val_sigl2  = this->m_sigma_l2.m_data[i];
                val_sigl2s = val_sigl2*val_sigl2;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                cosam      = gms::math::ceph_cosf(val_alpham);
                sinam      = gms::math::ceph_sinf(val_alpham);
#else 
                cosam      = std::cos(val_alpham);
                sinam      = std::sin(val_alpham);
#endif 
                r_rtherm   = 2.0f*val_l1*val_l2*cosam;
                sig_rtherm = 2.0f*val_sigl1*val_sigl2*val_rho*cosam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_sig   = gms::math::ceph_sqrtf(val_sigl1s+val_sigl2s+sig_rtherm);
#else 
                sqrt_sig   = std::sqrt(val_sigl1s+val_sigl2s+sig_rtherm);
#endif 
                this->m_r2.m_data[i] = sqrt_sig/sinam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_r     = gms::math::ceph_sqrtf(val_l1s+val_l2s+r_rtherm);
#else 
                sqrt_r     = std::sqrt(val_l1s+val_l2s+r_rtherm);
#endif 
                this->m_r1.m_data[i] = sqrt_r/sinam;
            }
        }
        ret_state = 1;
        return (ret_state);
    }
    else if(__builtin_expect(this->m_nsteps>16ull,1))
    {
        if(is_rho_zero==true)
        {
            for(i = 0ull; i!=ROUND_DOWN(this->m_nsteps,3ull); i += 4ull)
            {
                const std::size_t idx1 = i+1ull;
                const std::size_t idx2 = i+2ull;
                const std::size_t idx3 = i+3ull;
                val_alpham = this->m_alpham.m_data[i];
                val_l1     = this->m_l1.m_data[i];
                val_l1s    = val_l1*val_l1;
                val_l2     = this->m_l2.m_data[i];
                val_l2s    = val_l2*val_l2;
                val_sigl1  = this->m_sigma_l1.m_data[i];
                val_sigl1s = val_sigl1*val_sigl1;
                val_sigl2  = this->m_sigma_l2.m_data[i];
                val_sigl2s = val_sigl2*val_sigl2;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                cosam      = gms::math::ceph_cosf(val_alpham);
                sinam      = gms::math::ceph_sinf(val_alpham);
#else 
                cosam      = std::cos(val_alpham);
                sinam      = std::sin(val_alpham);
#endif 
                r_rtherm   = 2.0f*val_l1*val_l2*cosam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_sig   = gms::math::ceph_sqrtf(val_sigl1s+val_sigl2s);
#else 
                sqrt_sig   = std::sqrt(val_sigl1s+val_sigl2s);
#endif 
                this->m_r2.m_data[i] = sqrt_sig/sinam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_r     = gms::math::ceph_sqrtf(val_l1s+val_l2s+r_rtherm);
#else 
                sqrt_r     = std::sqrt(val_l1s+val_l2s+r_rtherm);
#endif 
                this->m_r1.m_data[i] = sqrt_r/sinam;

                val_alpham = this->m_alpham.m_data[idx1];
                val_l1     = this->m_l1.m_data[idx1];
                val_l1s    = val_l1*val_l1;
                val_l2     = this->m_l2.m_data[idx1];
                val_l2s    = val_l2*val_l2;
                val_sigl1  = this->m_sigma_l1.m_data[idx1];
                val_sigl1s = val_sigl1*val_sigl1;
                val_sigl2  = this->m_sigma_l2.m_data[idx1];
                val_sigl2s = val_sigl2*val_sigl2;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                cosam      = gms::math::ceph_cosf(val_alpham);
                sinam      = gms::math::ceph_sinf(val_alpham);
#else 
                cosam      = std::cos(val_alpham);
                sinam      = std::sin(val_alpham);
#endif 
                r_rtherm   = 2.0f*val_l1*val_l2*cosam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_sig   = gms::math::ceph_sqrtf(val_sigl1s+val_sigl2s);
#else 
                sqrt_sig   = std::sqrt(val_sigl1s+val_sigl2s);
#endif 
                this->m_r2.m_data[idx1] = sqrt_sig/sinam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_r     = gms::math::ceph_sqrtf(val_l1s+val_l2s+r_rtherm);
#else 
                sqrt_r     = std::sqrt(val_l1s+val_l2s+r_rtherm);
#endif 
                this->m_r1.m_data[idx1] = sqrt_r/sinam;

                val_alpham = this->m_alpham.m_data[idx2];
                val_l1     = this->m_l1.m_data[idx2];
                val_l1s    = val_l1*val_l1;
                val_l2     = this->m_l2.m_data[idx2];
                val_l2s    = val_l2*val_l2;
                val_sigl1  = this->m_sigma_l1.m_data[idx2];
                val_sigl1s = val_sigl1*val_sigl1;
                val_sigl2  = this->m_sigma_l2.m_data[idx2];
                val_sigl2s = val_sigl2*val_sigl2;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                cosam      = gms::math::ceph_cosf(val_alpham);
                sinam      = gms::math::ceph_sinf(val_alpham);
#else 
                cosam      = std::cos(val_alpham);
                sinam      = std::sin(val_alpham);
#endif 
                r_rtherm   = 2.0f*val_l1*val_l2*cosam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_sig   = gms::math::ceph_sqrtf(val_sigl1s+val_sigl2s);
#else 
                sqrt_sig   = std::sqrt(val_sigl1s+val_sigl2s);
#endif 
                this->m_r2.m_data[idx2] = sqrt_sig/sinam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_r     = gms::math::ceph_sqrtf(val_l1s+val_l2s+r_rtherm);
#else 
                sqrt_r     = std::sqrt(val_l1s+val_l2s+r_rtherm);
#endif 
                this->m_r1.m_data[idx2] = sqrt_r/sinam;

                val_alpham = this->m_alpham.m_data[idx3];
                val_l1     = this->m_l1.m_data[idx3];
                val_l1s    = val_l1*val_l1;
                val_l2     = this->m_l2.m_data[idx3];
                val_l2s    = val_l2*val_l2;
                val_sigl1  = this->m_sigma_l1.m_data[idx3];
                val_sigl1s = val_sigl1*val_sigl1;
                val_sigl2  = this->m_sigma_l2.m_data[idx3];
                val_sigl2s = val_sigl2*val_sigl2;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                cosam      = gms::math::ceph_cosf(val_alpham);
                sinam      = gms::math::ceph_sinf(val_alpham);
#else 
                cosam      = std::cos(val_alpham);
                sinam      = std::sin(val_alpham);
#endif 
                r_rtherm   = 2.0f*val_l1*val_l2*cosam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_sig   = gms::math::ceph_sqrtf(val_sigl1s+val_sigl2s);
#else 
                sqrt_sig   = std::sqrt(val_sigl1s+val_sigl2s);
#endif 
                this->m_r2.m_data[idx3] = sqrt_sig/sinam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_r     = gms::math::ceph_sqrtf(val_l1s+val_l2s+r_rtherm);
#else 
                sqrt_r     = std::sqrt(val_l1s+val_l2s+r_rtherm);
#endif 
                this->m_r1.m_data[idx3] = sqrt_r/sinam;
            }

            for(j = i; j<this->m_nsteps; ++j)  
            {
                val_alpham = this->m_alpham.m_data[j];
                val_l1     = this->m_l1.m_data[j];
                val_l1s    = val_l1*val_l1;
                val_l2     = this->m_l2.m_data[j];
                val_l2s    = val_l2*val_l2;
                val_sigl1  = this->m_sigma_l1.m_data[j];
                val_sigl1s = val_sigl1*val_sigl1;
                val_sigl2  = this->m_sigma_l2.m_data[j];
                val_sigl2s = val_sigl2*val_sigl2;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                cosam      = gms::math::ceph_cosf(val_alpham);
                sinam      = gms::math::ceph_sinf(val_alpham);
#else 
                cosam      = std::cos(val_alpham);
                sinam      = std::sin(val_alpham);
#endif 
                r_rtherm   = 2.0f*val_l1*val_l2*cosam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_sig   = gms::math::ceph_sqrtf(val_sigl1s+val_sigl2s);
#else 
                sqrt_sig   = std::sqrt(val_sigl1s+val_sigl2s);
#endif 
                this->m_r2.m_data[j] = sqrt_sig/sinam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_r     = gms::math::ceph_sqrtf(val_l1s+val_l2s+r_rtherm);
#else 
                sqrt_r     = std::sqrt(val_l1s+val_l2s+r_rtherm);
#endif 
                this->m_r1.m_data[j] = sqrt_r/sinam;
            }
        }
        else 
        {
            for(i = 0ull; i!=ROUND_DOWN(this->m_nsteps,3ull); i += 4ull)
            {
                const std::size_t idx1 = i+1ull;
                const std::size_t idx2 = i+2ull;
                const std::size_t idx3 = i+3ull;
                val_alpham = this->m_alpham.m_data[i];
                val_l1     = this->m_l1.m_data[i];
                val_rho    = this->m_rho.m_data[i];
                val_l1s    = val_l1*val_l1;
                val_l2     = this->m_l2.m_data[i];
                val_l2s    = val_l2*val_l2;
                val_sigl1  = this->m_sigma_l1.m_data[i];
                val_sigl1s = val_sigl1*val_sigl1;
                val_sigl2  = this->m_sigma_l2.m_data[i];
                val_sigl2s = val_sigl2*val_sigl2;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                cosam      = gms::math::ceph_cosf(val_alpham);
                sinam      = gms::math::ceph_sinf(val_alpham);
#else 
                cosam      = std::cos(val_alpham);
                sinam      = std::sin(val_alpham);
#endif 
                r_rtherm   = 2.0f*val_l1*val_l2*cosam;
                sig_rtherm = 2.0f*val_sigl1*val_sigl2*val_rho*cosam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_sig   = gms::math::ceph_sqrtf(val_sigl1s+val_sigl2s+sig_rtherm);
#else 
                sqrt_sig   = std::sqrt(val_sigl1s+val_sigl2s+sig_rtherm);
#endif 
                this->m_r2.m_data[i] = sqrt_sig/sinam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_r     = gms::math::ceph_sqrtf(val_l1s+val_l2s+r_rtherm);
#else 
                sqrt_r     = std::sqrt(val_l1s+val_l2s+r_rtherm);
#endif 
                this->m_r1.m_data[i] = sqrt_r/sinam;

                val_alpham = this->m_alpham.m_data[idx1];
                val_l1     = this->m_l1.m_data[idx1];
                val_rho    = this->m_rho.m_data[idx1];
                val_l1s    = val_l1*val_l1;
                val_l2     = this->m_l2.m_data[idx1];
                val_l2s    = val_l2*val_l2;
                val_sigl1  = this->m_sigma_l1.m_data[idx1];
                val_sigl1s = val_sigl1*val_sigl1;
                val_sigl2  = this->m_sigma_l2.m_data[idx1];
                val_sigl2s = val_sigl2*val_sigl2;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                cosam      = gms::math::ceph_cosf(val_alpham);
                sinam      = gms::math::ceph_sinf(val_alpham);
#else 
                cosam      = std::cos(val_alpham);
                sinam      = std::sin(val_alpham);
#endif 
                r_rtherm   = 2.0f*val_l1*val_l2*cosam;
                sig_rtherm = 2.0f*val_sigl1*val_sigl2*val_rho*cosam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_sig   = gms::math::ceph_sqrtf(val_sigl1s+val_sigl2s+sig_rtherm);
#else 
                sqrt_sig   = std::sqrt(val_sigl1s+val_sigl2s+sig_rtherm);
#endif 
                this->m_r2.m_data[idx1] = sqrt_sig/sinam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_r     = gms::math::ceph_sqrtf(val_l1s+val_l2s+r_rtherm);
#else 
                sqrt_r     = std::sqrt(val_l1s+val_l2s+r_rtherm);
#endif 
                this->m_r1.m_data[idx1] = sqrt_r/sinam;

                val_alpham = this->m_alpham.m_data[idx2];
                val_l1     = this->m_l1.m_data[idx2];
                val_rho    = this->m_rho.m_data[idx2];
                val_l1s    = val_l1*val_l1;
                val_l2     = this->m_l2.m_data[idx2];
                val_l2s    = val_l2*val_l2;
                val_sigl1  = this->m_sigma_l1.m_data[idx2];
                val_sigl1s = val_sigl1*val_sigl1;
                val_sigl2  = this->m_sigma_l2.m_data[idx2];
                val_sigl2s = val_sigl2*val_sigl2;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                cosam      = gms::math::ceph_cosf(val_alpham);
                sinam      = gms::math::ceph_sinf(val_alpham);
#else 
                cosam      = std::cos(val_alpham);
                sinam      = std::sin(val_alpham);
#endif 
                r_rtherm   = 2.0f*val_l1*val_l2*cosam;
                sig_rtherm = 2.0f*val_sigl1*val_sigl2*val_rho*cosam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_sig   = gms::math::ceph_sqrtf(val_sigl1s+val_sigl2s+sig_rtherm);
#else 
                sqrt_sig   = std::sqrt(val_sigl1s+val_sigl2s+sig_rtherm);
#endif 
                this->m_r2.m_data[idx2] = sqrt_sig/sinam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_r     = gms::math::ceph_sqrtf(val_l1s+val_l2s+r_rtherm);
#else 
                sqrt_r     = std::sqrt(val_l1s+val_l2s+r_rtherm);
#endif 
                this->m_r1.m_data[idx2] = sqrt_r/sinam;

                val_alpham = this->m_alpham.m_data[idx3];
                val_l1     = this->m_l1.m_data[idx3];
                val_rho    = this->m_rho.m_data[idx3];
                val_l1s    = val_l1*val_l1;
                val_l2     = this->m_l2.m_data[idx3];
                val_l2s    = val_l2*val_l2;
                val_sigl1  = this->m_sigma_l1.m_data[idx3];
                val_sigl1s = val_sigl1*val_sigl1;
                val_sigl2  = this->m_sigma_l2.m_data[idx3];
                val_sigl2s = val_sigl2*val_sigl2;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                cosam      = gms::math::ceph_cosf(val_alpham);
                sinam      = gms::math::ceph_sinf(val_alpham);
#else 
                cosam      = std::cos(val_alpham);
                sinam      = std::sin(val_alpham);
#endif 
                r_rtherm   = 2.0f*val_l1*val_l2*cosam;
                sig_rtherm = 2.0f*val_sigl1*val_sigl2*val_rho*cosam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_sig   = gms::math::ceph_sqrtf(val_sigl1s+val_sigl2s+sig_rtherm);
#else 
                sqrt_sig   = std::sqrt(val_sigl1s+val_sigl2s+sig_rtherm);
#endif 
                this->m_r2.m_data[idx3] = sqrt_sig/sinam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_r     = gms::math::ceph_sqrtf(val_l1s+val_l2s+r_rtherm);
#else 
                sqrt_r     = std::sqrt(val_l1s+val_l2s+r_rtherm);
#endif 
                this->m_r1.m_data[idx3] = sqrt_r/sinam;
            }

            for(j = i; j<this->m_nsteps; ++j)   
            {
                val_alpham = this->m_alpham.m_data[j];
                val_l1     = this->m_l1.m_data[j];
                val_rho    = this->m_rho.m_data[j];
                val_l1s    = val_l1*val_l1;
                val_l2     = this->m_l2.m_data[j];
                val_l2s    = val_l2*val_l2;
                val_sigl1  = this->m_sigma_l1.m_data[j];
                val_sigl1s = val_sigl1*val_sigl1;
                val_sigl2  = this->m_sigma_l2.m_data[j];
                val_sigl2s = val_sigl2*val_sigl2;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                cosam      = gms::math::ceph_cosf(val_alpham);
                sinam      = gms::math::ceph_sinf(val_alpham);
#else 
                cosam      = std::cos(val_alpham);
                sinam      = std::sin(val_alpham);
#endif 
                r_rtherm   = 2.0f*val_l1*val_l2*cosam;
                sig_rtherm = 2.0f*val_sigl1*val_sigl2*val_rho*cosam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_sig   = gms::math::ceph_sqrtf(val_sigl1s+val_sigl2s+sig_rtherm);
#else 
                sqrt_sig   = std::sqrt(val_sigl1s+val_sigl2s+sig_rtherm);
#endif 
                this->m_r2.m_data[j] = sqrt_sig/sinam;
#if (POSE_ERR_ESTIMATION_2D_USE_CEPHES) == 1
                sqrt_r     = gms::math::ceph_sqrtf(val_l1s+val_l2s+r_rtherm);
#else 
                sqrt_r     = std::sqrt(val_l1s+val_l2s+r_rtherm);
#endif 
                this->m_r1.m_data[j] = sqrt_r/sinam;
            }
        }
        ret_state = 2;
        return (ret_state);
    }
}

auto 
gms::radio_navigation
::operator<<(std::ostream &os, const pose_err_estimation_2D_t &rhs)->std::ostream &
{
    std::cout << typeid(rhs).name()     << "Begin: object state dump." << std::endl;
    std::cout << "m_nsteps:              " << rhs.m_nsteps << std::endl;
    std::cout << std::setw(4) << "r1" << std::setw(7) << "r2" << std::setw(17) << "l1_est_pdf" << std::setw(27) << "l2_est_pdf" 
              << std::setw(44) << "l12_joint_est_pdf" << std::setw(59) << "pose_error_prob" << std::endl;
    for(std::size_t i = 0ull; i<rhs.m_nsteps; ++i) 
    {
        os << std::fixed 
           << std::setprecision(7) 
           << std::setw(4) 
           << std::setw(7) 
           << std::setw(17)  
           << std::setw(27)  
           << std::setw(44) 
           << std::setw(59) 
           << rhs.m_r1.m_data[i]
           << rhs.m_r2.m_data[i]
           << rhs.m_l1_est_pdf.m_data[i]
           << rhs.m_l2_est_pdf.m_data[i]
           << rhs.m_l12_joint_est_pdf.m_data[i]
           << rhs.m_pose_error_prob.m_data[i]
           << std::endl;
    } 
    std::cout << typeid(rhs).name() << "End: object state dump." << std::endl;
    return (os);
}

auto 
gms::radio_navigation
::operator==(const pose_err_estimation_2D_t &lhs,
             const pose_err_estimation_2D_t &rhs)->pose_err_estimation_2D_compare_t
{
    if(__builtin_expect(lhs.m_nsteps!=rhs.m_nsteps,0)) { return (pose_err_estimation_2D_compare_t());}
    pose_err_estimation_2D_compare_t pose_err_comp(lhs.m_nsteps);
    float lhs_val;
    float rhs_val;
    bool  are_equal;
    for(std::size_t i = 0ull; i<lhs.m_nsteps; ++i)
    {
       lhs_val   = lhs.m_r1.m_data[i];
       rhs_val   = rhs.m_r1.m_data[i];
       are_equal = almostEqual2sComplement(lhs_val,rhs_val,7); 
       pose_err_comp.m_r1_state_comp.operator[](i) = are_equal;
       lhs_val   = lhs.m_r2.m_data[i];
       rhs_val   = rhs.m_r2.m_data[i];
       are_equal = almostEqual2sComplement(lhs_val,rhs_val,7); 
       pose_err_comp.m_r2_state_comp.operator[](i) = are_equal;
       lhs_val   = lhs.m_l1_est_pdf.m_data[i];
       rhs_val   = rhs.m_l1_est_pdf.m_data[i];
       are_equal = almostEqual2sComplement(lhs_val,rhs_val,7); 
       pose_err_comp.m_l1_est_pdf_state_comp.operator[](i) = are_equal;
       lhs_val   = lhs.m_l2_est_pdf.m_data[i];
       rhs_val   = rhs.m_l2_est_pdf.m_data[i];
       are_equal = almostEqual2sComplement(lhs_val,rhs_val,7); 
       pose_err_comp.m_l2_est_pdf_state_comp.operator[](i) = are_equal;
       lhs_val   = lhs.m_l12_joint_est_pdf.m_data[i];
       rhs_val   = rhs.m_l12_joint_est_pdf.m_data[i];
       are_equal = almostEqual2sComplement(lhs_val,rhs_val,7); 
       pose_err_comp.m_l12_est_joint_pdf_state_comp.operator[](i) = are_equal;
       lhs_val   = lhs.m_pose_error_prob.m_data[i];
       rhs_val   = rhs.m_pose_error_prob.m_data[i];
       are_equal = almostEqual2sComplement(lhs_val,rhs_val,7); 
       pose_err_comp.m_pose_prob_state_comp.operator[](i) = are_equal;
    }
    return (pose_err_comp);
}