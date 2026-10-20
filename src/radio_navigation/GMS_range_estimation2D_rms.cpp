
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
#include "GMS_range_estimation2D_rms.h"
#include "GMS_compare_fp_safe.h"


gms::radio_navigation
::range_estimation_2D_rms_t
::range_estimation_2D_rms_t(const std::size_t nsteps,
                            const darray_r4_t & radio_range_buoy_a,
                            const darray_r4_t & radio_range_buoy_b,
                            const darray_r4_t & radio_buoys_target_angle,
                            const darray_r4_t & rmse_signal_delay,
                            const float radio_buoys_ab_base) noexcept(false)
:
m_nsteps(nsteps),
m_radio_range_buoy_a(radio_range_buoy_a),
m_radio_range_buoy_b(radio_range_buoy_b),
m_radio_buoys_target_angle(radio_buoys_target_angle),
m_rmse_signal_delay(rmse_signal_delay),
m_hyperbolic_pose_shift(m_nsteps),
m_hyperbolic_pose_rms(m_nsteps),
m_radio_buoys_ab_base(radio_buoys_ab_base)
{

}

gms::radio_navigation
::range_estimation_2D_rms_t
::range_estimation_2D_rms_t(const std::size_t nsteps,
                            darray_r4_t && radio_range_buoy_a,
                            darray_r4_t && radio_range_buoy_b,
                            darray_r4_t && radio_buoys_target_angle,
                            darray_r4_t && rmse_signal_delay,
                            const float    radio_buoys_ab_base) noexcept(true)
:
m_nsteps(nsteps),
m_radio_range_buoy_a(std::move(radio_range_buoy_a)),
m_radio_range_buoy_b(std::move(radio_range_buoy_b)),
m_radio_buoys_target_angle(std::move(radio_buoys_target_angle)),
m_hyperbolic_pose_shift(std::move(hyperbolic_pose_shift)),
m_rmse_signal_delay(std::move(rmse_signal_delay)),
m_hyperbolic_pose_shift(m_nsteps),
m_hyperbolic_pose_rms(m_nsteps),
m_radio_buoys_ab_base(radio_buoys_ab_base) 
{

}

gms::radio_navigation
::range_estimation_2D_rms_t
::range_estimation_2D_rms_t(range_estimation_2D_rms_t && other) noexcept(true) 
:
m_nsteps(other.m_nsteps),
m_radio_range_buoy_a(std::move(other.m_radio_range_buoy_a)),
m_radio_range_buoy_b(std::move(other.m_radio_range_buoy_b)),
m_radio_buoys_target_angle(std::move(other.m_radio_buoys_target_angle)),
m_rmse_signal_delay(std::move(other.m_rmse_signal_delay)),
m_hyperbolic_pose_shift(std::move(other.m_hyperbolic_pose_shift)),
m_hyperbolic_pose_rms(std::move(other.m_hyperbolic_pose_rms)),
m_radio_buoys_ab_base(other.m_radio_buoys_ab_base)
{

}

gms::radio_navigation
::range_estimation_2D_rms_t
::range_estimation_2D_rms_t(const range_estimation_2D_rms_t & other) noexcept(false) 
:
m_nsteps(other.m_nsteps),
m_radio_range_buoy_a(other.m_radio_range_buoy_a),
m_radio_range_buoy_b(other.m_radio_range_buoy_b),
m_radio_buoys_target_angle(other.m_radio_buoys_target_angle),
m_rmse_signal_delay(other.m_rmse_signal_delay),
m_hyperbolic_pose_shift(other.m_hyperbolic_pose_shift),
m_hyperbolic_pose_rms(other.m_hyperbolic_pose_rms),
m_radio_buoys_ab_base(other.m_radio_buoys_ab_base)
{

}

gms::radio_navigation
::range_estimation_2D_rms_t
::~range_estimation_2D_rms_t()
{

}

gms::radio_navigation
::range_estimation_2D_rms_t &
gms::radio_navigation
::range_estimation_2D_rms_t
::operator=(const range_estimation_2D_rms_t & other) noexcept(false)
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_nsteps = other.m_nsteps;
    this->m_radio_range_buoy_a.operator=(other.m_radio_range_buoy_a);
    this->m_radio_range_buoy_b.operator=(other.m_radio_range_buoy_b);
    this->m_radio_buoys_target_angle.operator=(other.m_radio_buoys_target_angle);
    this->m_rmse_signal_delay.operator=(other.m_rmse_signal_delay);
    this->m_hyperbolic_pose_shift.operator=(other.m_hyperbolic_pose_shift);
    this->m_hyperbolic_pose_rms.operator=(other.m_hyperbolic_pose_rms);
    this->m_radio_buoys_ab_base = other.m_radio_buoys_ab_base;
    return (*this);
}

gms::radio_navigation
::range_estimation_2D_rms_t &
gms::radio_navigation
::range_estimation_2D_rms_t
::operator=(range_estimation_2D_rms_t && other) noexcept(true)
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_nsteps = other.m_nsteps;
    this->m_radio_range_buoy_a.operator=(std::move(other.m_radio_range_buoy_a));
    this->m_radio_range_buoy_b.operator=(std::move(other.m_radio_range_buoy_b));
    this->m_radio_buoys_target_angle.operator=(std::move(other.m_radio_buoys_target_angle));
    this->m_rmse_signal_delay.operator=(std::move(other.m_rmse_signal_delay));
    this->m_hyperbolic_pose_shift.operator=(std::move(other.m_hyperbolic_pose_shift));
    this->m_hyperbolic_pose_rms.operator=(std::move(other.m_hyperbolic_pose_rms));
    this->m_radio_buoys_ab_base = other.m_radio_buoys_ab_base;
    return (*this);
}

bool 
gms::radio_navigation
::range_estimation_2D_rms_t
::analytic_grad_D_param(const darray_r4_t & coord_x,
                        const darray_r4_t & coord_y,
                        darray_r4_t & grad_D) 
{
    bool are_ranges_equal = ((this->m_nsteps==coord_x.mnx) && 
                             (this->m_nsteps==coord_y.mnx) && 
                             (this->m_nsteps==grad_D.mnx));
    if(__builtin_expect(are_ranges_equal!=true,0)) { return (false);}
    const float * __restrict__ p_grad_D  = &grad_D.m_data[0];
    const float * __restrict__ p_coord_x = &coord_x.m_data[0];
    const float * __restrict__ p_coord_y = &coord_y.m_data[0]; 
    std::size_t i,j;
    const float d_sqr = 0.25f*(this->m_radio_buoys_ab_base*this->m_radio_buoys_ab_base);
    float val_grad_D;
    float val_coord_x;
    float val_coord_y;
    float val_Da;
    float val_Db;
    float xsqr;
    float ysqr;
    float DaDb;
    float num;
    float ratio;
    float sqrt_val;
    if(this->m_nsteps>=1ull && this->m_nsteps<=16ull)
    {
        for(i = 0ull; i<this->m_nsteps; ++i) 
        {
            val_grad_D       = p_grad_D[i];
            val_coord_x      = p_coord_x[i];
            val_coord_y      = p_coord_y[i];
            val_Da           = this->m_radio_range_buoy_a[i];
            val_Db           = this->m_radio_range_buoy_b[i];
            xsqr             = val_coord_x*val_coord_x;
            ysqr             = val_coord_y*val_coord_y;
            DaDb             = val_Da*val_Db;
            num              = xsqr+ysqr-d_sqr;
            ratio            = 2.0f*(1.0f*(num/DaDb));
            sqrt_val         = std::sqrt(ratio);
            p_grad_D[i]      = sqrt_val;
        }
        return (true);
    }
    else if(this->m_nsteps>16ull)
    {
        for(i = 0ull; i!=ROUND_DOWN(this->m_nsteps,3ull); i += 4ull) 
        {
            const std::size_t idx1 = i+1ull;
            const std::size_t idx2 = i+2ull;
            const std::size_t idx3 = i+3ull;
            val_grad_D       = p_grad_D[i];
            val_coord_x      = p_coord_x[i];
            val_coord_y      = p_coord_y[i];
            val_Da           = this->m_radio_range_buoy_a[i];
            val_Db           = this->m_radio_range_buoy_b[i];
            xsqr             = val_coord_x*val_coord_x;
            ysqr             = val_coord_y*val_coord_y;
            DaDb             = val_Da*val_Db;
            num              = xsqr+ysqr-d_sqr;
            ratio            = 2.0f*(1.0f*(num/DaDb));
            sqrt_val         = std::sqrt(ratio);
            p_grad_D[i] = sqrt_val;
            val_grad_D       = p_grad_D[idx1];
            val_coord_x      = p_coord_x[idx1];
            val_coord_y      = p_coord_y[idx1];
            val_Da           = this->m_radio_range_buoy_a[idx1];
            val_Db           = this->m_radio_range_buoy_b[idx1];
            xsqr             = val_coord_x*val_coord_x;
            ysqr             = val_coord_y*val_coord_y;
            DaDb             = val_Da*val_Db;
            num              = xsqr+ysqr-d_sqr;
            ratio            = 2.0f*(1.0f*(num/DaDb));
            sqrt_val         = std::sqrt(ratio);
            p_grad_D[idx1] = sqrt_val;
            val_grad_D       = p_grad_D[idx2];
            val_coord_x      = p_coord_x[idx2];
            val_coord_y      = p_coord_y[idx2];
            val_Da           = this->m_radio_range_buoy_a[idx2];
            val_Db           = this->m_radio_range_buoy_b[idx2];
            xsqr             = val_coord_x*val_coord_x;
            ysqr             = val_coord_y*val_coord_y;
            DaDb             = val_Da*val_Db;
            num              = xsqr+ysqr-d_sqr;
            ratio            = 2.0f*(1.0f*(num/DaDb));
            sqrt_val         = std::sqrt(ratio);
            p_grad_D[idx2] = sqrt_val;
            val_grad_D       = p_grad_D[idx3];
            val_coord_x      = p_coord_x[idx3];
            val_coord_y      = p_coord_y[idx3];
            val_Da           = this->m_radio_range_buoy_a[idx3];
            val_Db           = this->m_radio_range_buoy_b[idx3];
            xsqr             = val_coord_x*val_coord_x;
            ysqr             = val_coord_y*val_coord_y;
            DaDb             = val_Da*val_Db;
            num              = xsqr+ysqr-d_sqr;
            ratio            = 2.0f*(1.0f*(num/DaDb));
            sqrt_val         = std::sqrt(ratio);
            p_grad_D[idx3] = sqrt_val;
        }
        
        for(j = i; j<this->m_nsteps; ++j)   
        {
            val_grad_D       = p_grad_D[j];
            val_coord_x      = p_coord_x[j];
            val_coord_y      = p_coord_y[j];
            val_Da           = this->m_radio_range_buoy_a[j];
            val_Db           = this->m_radio_range_buoy_b[j];
            xsqr             = val_coord_x*val_coord_x;
            ysqr             = val_coord_y*val_coord_y;
            DaDb             = val_Da*val_Db;
            num              = xsqr+ysqr-d_sqr;
            ratio            = 2.0f*(1.0f*(num/DaDb));
            sqrt_val         = std::sqrt(ratio);
            p_grad_D[j] = sqrt_val;
        }
        return (true);
    }
}


bool 
gms::radio_navigation
::range_estimation_2D_rms_t
::compute_hyberbolic_pose_rms( const darray_r4_t & coord_x,
                               const darray_r4_t & coord_y,
                               const darray_r4_t & dist_err_meas,
                               std::optional<float> c_speed = std::nullopt)
{
    bool are_ranges_equal = ((this->m_nsteps==coord_x.mnx) && 
                             (this->m_nsteps==coord_y.mnx) && 
                             (this->m_nsteps==dist_err_meas.mnx));
    if(__builtin_expect(are_ranges_equal!=true,0)) { return (false);}
    const float * __restrict__ p_coord_x       = &coord_x.m_data[0];
    const float * __restrict__ p_coord_y       = &coord_y.m_data[0]; 
    const float * __restrict__ p_dist_err_meas = &dist_err_meas.m_data[0];
    std::size_t i,j;
    constexpr float half = 0.5f;
    float val_x;
    float x_term;
    float left_sqr;
    float right_sqr;
    float y_term;
    float y_sqr;
    float psi;
    float val_y;
    float meas_err;
    float Dp;
    float psi;
    float sin_val;
    float sin_term;
    float sigma_l;
    float sigma_d;
    float temp;
    float c;
    if(c_speed) {c = *c_speed;}
    else {c = 299792458.0f;}
    const float half_d = 0.5f*this->m_radio_buoys_ab_base;
    if(this->m_nsteps>=1 && this->m_nsteps<=16)
    {
        for(i = 0ull; i<this->m_nsteps; ++i) 
        {
           val_x       = p_coord_x[i];
           val_y       = p_coord_y[i];
           psi         = this->m_radio_buoys_target_angle[i];
           meas_err    = p_dist_err_meas[i];
           temp        = val_x+half_d;
           y_sqr       = val_y*val_y;
           sin_term    = 2.0f*gms::math::ceph_sinf(0.5f*psi);
           x_term      = val_x+half_d;
           y_term      = val_y-half_d;
           left_sqr    = gms::math::ceph_sqrtf(__builtin_fmaf(x_term,x_term,y_sqr));
           right_sqr   = gms::math::ceph_sqrtf(__builtin_fmaf(y_term,y_term,y_sqr));
           Dp          = left_sqr-right_sqr;
           sigma_d     = (meas_err+Dp)/sin_term;
           this->m_hyperbolic_pose_shift[i] = sigma_d;
           sigma_l     = (c*sigma_l)/sin_term;
           this->m_hyperbolic_pose_rms[i] = sigma_l;
        }
        return (true);
    }
    else if(this->m_nsteps>16ull)
    {
        for(i = 0ull; i!=ROUND_DOWN(this->m_nsteps,3ull); i += 4ull)
        {
            const std::size_t idx1 = i+1ull;
            const std::size_t idx2 = i+2ull;
            const std::size_t idx3 = i+3ull;
            val_x       = p_coord_x[i];
            val_y       = p_coord_y[i];
            psi         = this->m_radio_buoys_target_angle[i];
            meas_err    = p_dist_err_meas[i];
            temp        = val_x+half_d;
            y_sqr       = val_y*val_y;
            sin_term    = 2.0f*gms::math::ceph_sinf(0.5f*psi);
            x_term      = val_x+half_d;
            y_term      = val_y-half_d;
            left_sqr    = gms::math::ceph_sqrtf(__builtin_fmaf(x_term,x_term,y_sqr));
            right_sqr   = gms::math::ceph_sqrtf(__builtin_fmaf(y_term,y_term,y_sqr));
            Dp          = left_sqr-right_sqr;
            sigma_d     = (meas_err+Dp)/sin_term;
            this->m_hyperbolic_pose_shift[i] = sigma_d;
            sigma_l     = (c*sigma_l)/sin_term;
            this->m_hyperbolic_pose_rms[i] = sigma_l;
            val_x       = p_coord_x[idx1];
            val_y       = p_coord_y[idx1];
            psi         = this->m_radio_buoys_target_angle[idx1];
            meas_err    = p_dist_err_meas[idx1];
            temp        = val_x+half_d;
            y_sqr       = val_y*val_y;
            sin_term    = 2.0f*gms::math::ceph_sinf(0.5f*psi);
            x_term      = val_x+half_d;
            y_term      = val_y-half_d;
            left_sqr    = gms::math::ceph_sqrtf(__builtin_fmaf(x_term,x_term,y_sqr));
            right_sqr   = gms::math::ceph_sqrtf(__builtin_fmaf(y_term,y_term,y_sqr));
            Dp          = left_sqr-right_sqr;
            sigma_d     = (meas_err+Dp)/sin_term;
            this->m_hyperbolic_pose_shift[idx1] = sigma_d;
            sigma_l     = (c*sigma_l)/sin_term;
            this->m_hyperbolic_pose_rms[idx1] = sigma_l;
            val_x       = p_coord_x[idx2];
            val_y       = p_coord_y[idx2];
            psi         = this->m_radio_buoys_target_angle[idx2];
            meas_err    = p_dist_err_meas[idx2];
            temp        = val_x+half_d;
            y_sqr       = val_y*val_y;
            sin_term    = 2.0f*gms::math::ceph_sinf(0.5f*psi);
            x_term      = val_x+half_d;
            y_term      = val_y-half_d;
            left_sqr    = gms::math::ceph_sqrtf(__builtin_fmaf(x_term,x_term,y_sqr));
            right_sqr   = gms::math::ceph_sqrtf(__builtin_fmaf(y_term,y_term,y_sqr));
            Dp          = left_sqr-right_sqr;
            sigma_d     = (meas_err+Dp)/sin_term;
            this->m_hyperbolic_pose_shift[idx2] = sigma_d;
            sigma_l     = (c*sigma_l)/sin_term;
            this->m_hyperbolic_pose_rms[idx2] = sigma_l;
            val_x       = p_coord_x[idx3];
            val_y       = p_coord_y[idx3];
            psi         = this->m_radio_buoys_target_angle[idx3];
            meas_err    = p_dist_err_meas[idx3];
            temp        = val_x+half_d;
            y_sqr       = val_y*val_y;
            sin_term    = 2.0f*gms::math::ceph_sinf(0.5f*psi);
            x_term      = val_x+half_d;
            y_term      = val_y-half_d;
            left_sqr    = gms::math::ceph_sqrtf(__builtin_fmaf(x_term,x_term,y_sqr));
            right_sqr   = gms::math::ceph_sqrtf(__builtin_fmaf(y_term,y_term,y_sqr));
            Dp          = left_sqr-right_sqr;
            sigma_d     = (meas_err+Dp)/sin_term;
            this->m_hyperbolic_pose_shift[idx3] = sigma_d;
            sigma_l     = (c*sigma_l)/sin_term;
            this->m_hyperbolic_pose_rms[idx3] = sigma_l;
        }

        for(j = i; j<this->m_nsteps; ++j)   
        {
            val_x       = p_coord_x[j];
            val_y       = p_coord_y[j];
            psi         = this->m_radio_buoys_target_angle[j];
            meas_err    = p_dist_err_meas[j];
            temp        = val_x+half_d;
            y_sqr       = val_y*val_y;
            sin_term    = 2.0f*gms::math::ceph_sinf(0.5f*psi);
            x_term      = val_x+half_d;
            y_term      = val_y-half_d;
            left_sqr    = gms::math::ceph_sqrtf(__builtin_fmaf(x_term,x_term,y_sqr));
            right_sqr   = gms::math::ceph_sqrtf(__builtin_fmaf(y_term,y_term,y_sqr));
            Dp          = left_sqr-right_sqr;
            sigma_d     = (meas_err+Dp)/sin_term;
            this->m_hyperbolic_pose_shift[j] = sigma_d;
            sigma_l     = (c*sigma_l)/sin_term;
            this->m_hyperbolic_pose_rms[j] = sigma_l;
        }
        return (true);
    }
}

auto 
gms::radio_navigation
::operator<<(std::ostream &os,
             gms::radio_navigation::range_estimation_2D_rms_t &rhs)->std::ostream &
{
    std::cout << typeid(rhs).name()     << "Begin: object state dump." << std::endl;
    std::cout << "m_nsteps:              " << rhs.m_nsteps << std::endl;
    std::cout << "m_radio_buoys_base_ab: " << std::fixed << std::setprecision(7) << rhs.m_radio_buoys_ab_base << std::endl;
    std::cout << std::setw(8) << "Hyperbolic Pose Shift" << std::setw(22) << "Hyperbolic Pose RMS" << std::endl;
    for(std::size_t i = 0ull; i<rhs.m_nsteps; ++i) 
    {
        os << std::fixed << std::setprecision(7) << std::setw(8) << std::setw(22) << 
        rhs.m_hyperbolic_pose_shift.m_data[i] << 
        rhs.m_hyperbolic_pose_rms.m_data[i] << std::endl;
    }
    std::cout << typeid(rhs).name() << "End: object state dump." << std::endl;
    return (os);
}

auto 
gms::radio_navigation
::operator==(const gms::radio_navigation::range_estimation_2D_rms_t &lhs,
             const gms::radio_navigation::range_estimation_2D_rms_t &rhs)->std::vector<bool>
{
    const bool are_ranges_equal = lhs.m_nsteps==rhs.m_nsteps;
    if(__builtin_expect(are_ranges_equal!=,0)) {return (std::vector<bool>{});}
    const float * __restrict__ p_lhs = &lhs.m_hyperbolic_pose_rms.m_data[0];
    const float * __restrict__ p_rhs = &rhs.m_hyperbolic_pose_rms.m_data[0];
    std::vector<bool> vbool(lhs.m_nsteps);
    for(std::size_t i = 0ull; i<lhs.m_nsteps; ++i) 
    {
        float lhs_val = p_lhs[i];
        float rhs_val = p_rhs[i];
        const bool are_equal = almostEqual2sComplement(lhs_val,rhs_val,7);
        vbool.operator[](i)  = are_equal;
    }
    return (vbool);
}

auto 
gms::radio_navigation
::operator!=(const gms::radio_navigation::range_estimation_2D_rms_t &lhs,
             const gms::radio_navigation::range_estimation_2D_rms_t &rhs)->std::vector<bool>
{
    const bool are_ranges_equal = lhs.m_nsteps==rhs.m_nsteps;
    if(__builtin_expect(are_ranges_equal!=,0)) {return (std::vector<bool>{});}
    const float * __restrict__ p_lhs = &lhs.m_hyperbolic_pose_rms.m_data[0];
    const float * __restrict__ p_rhs = &rhs.m_hyperbolic_pose_rms.m_data[0];
    std::vector<bool> vbool(lhs.m_nsteps);
    for(std::size_t i = 0ull; i<lhs.m_nsteps; ++i) 
    {
        float lhs_val = p_lhs[i];
        float rhs_val = p_rhs[i];
        const bool are_not_equal = !almostEqual2sComplement(lhs_val,rhs_val,7);
        vbool.operator[](i)  = are_not_equal;
    }
    return (vbool);
}


