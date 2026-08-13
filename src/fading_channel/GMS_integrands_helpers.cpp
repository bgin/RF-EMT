/*
 * Copyright (C) Bernard Gingold, 2020-2026 
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

#include <immintrin.h>
#include <random>
#include <functional>
#include <cmath>
#include <algorithm>
#include "GMS_integrands_helpers.h"

namespace 
{
static inline 
std::int32_t 
print_double(const char * val_name,double f, const std::int32_t offset) 
{
  (*((std::uint64_t *)&f)) += offset;
  std::int32_t printf_ret = std::printf("value: %s,%1.17g,%017llx,%llu\n",val_name,f, *(std::uint64_t *)&f, *(std::uint64_t *)&f);
  return (printf_ret);
}
}


std::int32_t 
gms::fading_channel
::integrand_4_1_gauss_Q_arg_checker(double * __restrict__ arg_range,
                                    const std::int32_t n_func_args,
                                    const std::int32_t n_gauss_Q_vals,
                                    const double lo1,
                                    const double hi1,
                                    const double hi2)
                                        
{
    thread_local std::uniform_real_distribution<double> rv_func_arg;
    thread_local std::mt19937 rv_func_arg_gen;
    thread_local std::uint64_t seed_func_arg{};
    thread_local std::uniform_real_distribution<double> rv_gauss_q;
    thread_local std::mt19937 rv_gauss_q_gen;
    thread_local std::uint64_t seed_gauss_q{};
    rv_func_arg = std::uniform_real_distribution<double>(lo1,hi1);
    seed_func_arg = __rdtsc();
    rv_func_arg_gen = std::mt19937(seed_func_arg);
    const bool  is_arg_range_null = (arg_range==nullptr);
    if(!is_arg_range_null)
    {
        for(std::int32_t i{0}; i<n_func_args; ++i)
        {
            const std::int32_t outer_idx = i*n_gauss_Q_vals;
            const double x = rv_func_arg.operator()(rv_func_arg_gen);
            rv_gauss_q = std::uniform_real_distribution<double>(x,hi2);
            seed_gauss_q = __rdtsc();
            rv_gauss_q_gen = std::mt19937(seed_gauss_q);
            for(std::int32_t j{0}; j<n_gauss_Q_vals; ++j)   
            {
                const std::int32_t inner_idx = outer_idx+j;
                const double y = rv_gauss_q.operator()(rv_gauss_q_gen);
                const double exp_arg = 0.5*(y*y);
                arg_range[inner_idx] = exp_arg;
            }
        }
    }
    else 
    {
        for(std::int32_t i{0}; i<n_func_args; ++i)
        {
            const double x = rv_func_arg.operator()(rv_func_arg_gen);
            rv_gauss_q = std::uniform_real_distribution<double>(x,hi2);
            seed_gauss_q = __rdtsc();
            rv_gauss_q_gen = std::mt19937(seed_gauss_q);
            for(std::int32_t j{0}; j<n_gauss_Q_vals; ++j)   
            {
                const double y = rv_gauss_q.operator()(rv_gauss_q_gen);
                const double exp_arg = 0.5*(y*y);
                [[maybe_unused]] std::int32_t printf_ret = print_double("exp_arg",exp_arg,0);
            }
        }
    }
    return (0);
}

std::int32_t 
gms::fading_channel
::integrand_4_2_gauss_Q_arg_checker(double * __restrict__ arg_range,
                                    const std::int32_t n_func_args,
                                    const std::int32_t n_gauss_Q_vals,
                                    const double lo1,
                                    const double hi1,
                                    const double lo2,
                                    const double hi2)
{
    thread_local std::uniform_real_distribution<double> rv_func_arg;
    thread_local std::mt19937 rv_func_arg_gen;
    thread_local std::uint64_t seed_func_arg{};
    thread_local std::uniform_real_distribution<double> rv_gauss_q;
    thread_local std::mt19937 rv_gauss_q_gen;
    thread_local std::uint64_t seed_gauss_q{};
    rv_func_arg = std::uniform_real_distribution<double>(lo1,hi1);
    seed_func_arg = __rdtsc();
    rv_func_arg_gen = std::mt19937(seed_func_arg);
    rv_gauss_q = std::uniform_real_distribution<double>(lo2,hi2);
    seed_gauss_q = __rdtsc();
    rv_gauss_q_gen = std::mt19937(seed_gauss_q);
    const bool  is_arg_range_null = (arg_range==nullptr);
    if(!is_arg_range_null)
    {
       for(std::int32_t i{0}; i<n_func_args; ++i) 
       {
           const std::int32_t outer_idx = i*n_gauss_Q_vals;
           const double x = rv_func_arg.operator()(rv_func_arg_gen);
           const double xsqr{x*x};
           for(std::int32_t j{0}; j<n_gauss_Q_vals; ++j)   
           {
               const std::int32_t inner_idx = outer_idx+j;
               const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
               const double sin_val{std::sin(theta)};
               const double sin_denom{2.0*(sin_val*sin_val)};
               const double exp_arg{xsqr/sin_denom};
               arg_range[inner_idx] = exp_arg;
           }
       }
    }
    else 
    {
       for(std::int32_t i{0}; i<n_func_args; ++i) 
       {
           const double x = rv_func_arg.operator()(rv_func_arg_gen);
           const double xsqr{x*x};
           for(std::int32_t j{0}; j<n_gauss_Q_vals; ++j)   
           {
               const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
               const double sin_val{std::sin(theta)};
               const double sin_denom{2.0*(sin_val*sin_val)};
               const double exp_arg{xsqr/sin_denom};
               [[maybe_unused]] std::int32_t printf_ret = print_double("exp_arg",exp_arg,0);
           }
       }
    }
    return (0);
}

std::int32_t 
gms::fading_channel
::integrand_4_6_gauss_Q_sin_arg_checker(double * __restrict__ arg_range,
                                       const std::int32_t n_func_args,
                                       const std::int32_t n_gauss_Q_vals,
                                       const double lo1,
                                       const double hi1,
                                       const double lo2,
                                       const double hi2,
                                       const double lo3,
                                       const double hi3,
                                       const double lo4,
                                       const double hi4)
{
    thread_local std::uniform_real_distribution<double> rv_func_arg_x1;
    thread_local std::mt19937 rv_func_arg_x1_gen;
    thread_local std::uint64_t seed_func_x1_arg{};
    thread_local std::uniform_real_distribution<double> rv_func_arg_y1;
    thread_local std::mt19937 rv_func_arg_y1_gen;
    thread_local std::uint64_t seed_func_y1_arg{};
    thread_local std::uniform_real_distribution<double> rv_func_arg_rho;
    thread_local std::mt19937 rv_func_arg_rho_gen;
    thread_local std::uint64_t seed_func_rho_arg{};
    thread_local std::uniform_real_distribution<double> rv_gauss_q;
    thread_local std::mt19937 rv_gauss_q_gen;
    thread_local std::uint64_t seed_gauss_q{};
    rv_func_arg_x1 = std::uniform_real_distribution<double>(lo1,hi1);
    seed_func_x1_arg = __rdtsc();
    rv_func_arg_x1_gen = std::mt19937(seed_func_x1_arg);
    rv_func_arg_y1 = std::uniform_real_distribution<double>(lo2,hi2);
    seed_func_y1_arg = __rdtsc();
    rv_func_arg_y1_gen = std::mt19937(seed_func_y1_arg);
    rv_func_arg_rho = std::uniform_real_distribution<double>(lo3,hi3);
    seed_func_rho_arg = __rdtsc();
    rv_func_arg_rho_gen = std::mt19937(seed_func_rho_arg);
    const bool  is_arg_range_null = (arg_range==nullptr);
    if(!is_arg_range_null)
    {
        for(std::int32_t i{0}; i<n_func_args; ++i) 
        {
            const std::int32_t outer_idx = i*n_gauss_Q_vals;
            const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
            const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
            const double rho= rv_func_arg_rho.operator()(rv_func_arg_rho_gen);
            const double one_m_rho{1.0f-(rho*rho)};
            const double sqr_1_m_rho{std::sqrt(one_m_rho)};
            const double phi_arg = y1/x1;
            const double phi_s = std::atan(phi_arg);
            const double tmp{std::sin(phi_s)};
            const double sin_or_cos_sqr = tmp*tmp;
            rv_gauss_q = std::uniform_real_distribution<double>(lo4,hi4-phi_s);
            seed_gauss_q = __rdtsc();
            rv_gauss_q_gen = std::mt19937(seed_gauss_q);
            const double squared_sum{(x1*x1)+(y1*y1)}; 
            const double S_hat{std::sqrt(squared_sum)};
            for(std::int32_t j{0}; j<n_gauss_Q_vals; ++j)   
            {
                const std::int32_t inner_idx = outer_idx+j;
                const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
                double tsin{std::sin(theta)};
                const double sin_sqr_tht{tsin*tsin};
                const double sin2theta{std::sin(theta+theta)};
                const double one_m_sin2theta{1.0-(rho*sin2theta)};
                const double S_hat_half{0.5*(S_hat*S_hat)};
                const double ratio1{sqr_1_m_rho/one_m_sin2theta};
                const double sin_sin_or_cos_rat = sin_or_cos_sqr/sin_sqr_tht;
                const double ratio2{one_m_sin2theta/one_m_rho};
                const double exp_arg{S_hat_half*ratio2*sin_sin_or_cos_rat};
                arg_range[inner_idx] = exp_arg;
            }
        }
        
    }
    else 
    {
        for(std::int32_t i{0}; i<n_func_args; ++i) 
        {
            const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
            const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
            const double rho= rv_func_arg_rho.operator()(rv_func_arg_rho_gen);
            const double one_m_rho{1.0f-(rho*rho)};
            const double sqr_1_m_rho{std::sqrt(one_m_rho)};
            const double phi_arg = y1/x1;
            const double phi_s = std::atan(phi_arg);
            const double tmp{std::sin(phi_s)};
            const double sin_or_cos_sqr = tmp*tmp;
            rv_gauss_q = std::uniform_real_distribution<double>(lo4,hi4-phi_s);
            seed_gauss_q = __rdtsc();
            rv_gauss_q_gen = std::mt19937(seed_gauss_q);
            const double squared_sum{(x1*x1)+(y1*y1)}; 
            const double S_hat{std::sqrt(squared_sum)};
            for(std::int32_t j{0}; j<n_gauss_Q_vals; ++j)   
            {
                const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
                double tsin{std::sin(theta)};
                const double sin_sqr_tht{tsin*tsin};
                const double sin2theta{std::sin(theta+theta)};
                const double one_m_sin2theta{1.0-(rho*sin2theta)};
                const double S_hat_half{0.5*(S_hat*S_hat)};
                const double ratio1{sqr_1_m_rho/one_m_sin2theta};
                const double sin_sin_or_cos_rat = sin_or_cos_sqr/sin_sqr_tht;
                const double ratio2{one_m_sin2theta/one_m_rho};
                const double exp_arg{S_hat_half*ratio2*sin_sin_or_cos_rat};
                [[maybe_unused]] std::int32_t printf_ret = print_double("exp_arg",exp_arg,0);
            }
        }
    }
    return (0);
}

std::int32_t 
gms::fading_channel
::integrand_4_6_gauss_Q_cos_arg_checker(double * __restrict__ arg_range,
                                       const std::int32_t n_func_args,
                                       const std::int32_t n_gauss_Q_vals,
                                       const double lo1,
                                       const double hi1,
                                       const double lo2,
                                       const double hi2,
                                       const double lo3,
                                       const double hi3,
                                       const double lo4,
                                       const double hi4)
{
    thread_local std::uniform_real_distribution<double> rv_func_arg_x1;
    thread_local std::mt19937 rv_func_arg_x1_gen;
    thread_local std::uint64_t seed_func_x1_arg{};
    thread_local std::uniform_real_distribution<double> rv_func_arg_y1;
    thread_local std::mt19937 rv_func_arg_y1_gen;
    thread_local std::uint64_t seed_func_y1_arg{};
    thread_local std::uniform_real_distribution<double> rv_func_arg_rho;
    thread_local std::mt19937 rv_func_arg_rho_gen;
    thread_local std::uint64_t seed_func_rho_arg{};
    thread_local std::uniform_real_distribution<double> rv_gauss_q;
    thread_local std::mt19937 rv_gauss_q_gen;
    thread_local std::uint64_t seed_gauss_q{};
    rv_func_arg_x1 = std::uniform_real_distribution<double>(lo1,hi1);
    seed_func_x1_arg = __rdtsc();
    rv_func_arg_x1_gen = std::mt19937(seed_func_x1_arg);
    rv_func_arg_y1 = std::uniform_real_distribution<double>(lo2,hi2);
    seed_func_y1_arg = __rdtsc();
    rv_func_arg_y1_gen = std::mt19937(seed_func_y1_arg);
    rv_func_arg_rho = std::uniform_real_distribution<double>(lo3,hi3);
    seed_func_rho_arg = __rdtsc();
    rv_func_arg_rho_gen = std::mt19937(seed_func_rho_arg);
    const bool  is_arg_range_null = (arg_range==nullptr);
    if(!is_arg_range_null)
    {
        for(std::int32_t i{0}; i<n_func_args; ++i) 
        {
            const std::int32_t outer_idx = i*n_gauss_Q_vals;
            const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
            const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
            const double rho= rv_func_arg_rho.operator()(rv_func_arg_rho_gen);
            const double one_m_rho{1.0f-(rho*rho)};
            const double sqr_1_m_rho{std::sqrt(one_m_rho)};
            const double phi_arg = y1/x1;
            const double phi_s = std::atan(phi_arg);
            const double tmp{std::cos(phi_s)};
            const double sin_or_cos_sqr = tmp*tmp;
            rv_gauss_q = std::uniform_real_distribution<double>(lo4,hi4-phi_s);
            seed_gauss_q = __rdtsc();
            rv_gauss_q_gen = std::mt19937(seed_gauss_q);
            const double squared_sum{(x1*x1)+(y1*y1)}; 
            const double S_hat{std::sqrt(squared_sum)};
            for(std::int32_t j{0}; j<n_gauss_Q_vals; ++j)   
            {
                const std::int32_t inner_idx = outer_idx+j;
                const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
                double tsin{std::sin(theta)};
                const double sin_sqr_tht{tsin*tsin};
                const double sin2theta{std::sin(theta+theta)};
                const double one_m_sin2theta{1.0-(rho*sin2theta)};
                const double S_hat_half{0.5*(S_hat*S_hat)};
                const double ratio1{sqr_1_m_rho/one_m_sin2theta};
                const double sin_sin_or_cos_rat = sin_or_cos_sqr/sin_sqr_tht;
                const double ratio2{one_m_sin2theta/one_m_rho};
                const double exp_arg{S_hat_half*ratio2*sin_sin_or_cos_rat};
                arg_range[inner_idx] = exp_arg;
            }
        }
        
    }
    else 
    {
        for(std::int32_t i{0}; i<n_func_args; ++i) 
        {
            const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
            const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
            const double rho= rv_func_arg_rho.operator()(rv_func_arg_rho_gen);
            const double one_m_rho{1.0f-(rho*rho)};
            const double sqr_1_m_rho{std::sqrt(one_m_rho)};
            const double phi_arg = y1/x1;
            const double phi_s = std::atan(phi_arg);
            const double tmp{std::sin(phi_s)};
            const double sin_or_cos_sqr = tmp*tmp;
            rv_gauss_q = std::uniform_real_distribution<double>(lo4,hi4-phi_s);
            seed_gauss_q = __rdtsc();
            rv_gauss_q_gen = std::mt19937(seed_gauss_q);
            const double squared_sum{(x1*x1)+(y1*y1)}; 
            const double S_hat{std::sqrt(squared_sum)};
            for(std::int32_t j{0}; j<n_gauss_Q_vals; ++j)   
            {
                const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
                double tsin{std::sin(theta)};
                const double sin_sqr_tht{tsin*tsin};
                const double sin2theta{std::sin(theta+theta)};
                const double one_m_sin2theta{1.0-(rho*sin2theta)};
                const double S_hat_half{0.5*(S_hat*S_hat)};
                const double ratio1{sqr_1_m_rho/one_m_sin2theta};
                const double sin_sin_or_cos_rat = sin_or_cos_sqr/sin_sqr_tht;
                const double ratio2{one_m_sin2theta/one_m_rho};
                const double exp_arg{S_hat_half*ratio2*sin_sin_or_cos_rat};
                [[maybe_unused]] std::int32_t printf_ret = print_double("exp_arg",exp_arg,0);
            }
        }
    }
    return (0);
}

std::int32_t 
gms::fading_channel
::integrand_4_10_marcum_Q_arg_checker( double * __restrict__ arg_range,
                                       const std::int32_t n_func_args,
                                       const std::int32_t n_gauss_Q_vals,
                                       const double lo1,
                                       const double hi1,
                                       const double lo2,
                                       const double hi2)
{
    thread_local std::uniform_real_distribution<double> rv_s_param;
    thread_local std::mt19937 rv_s_param_gen;
    thread_local std::uint64_t seed_s_param{};
    thread_local std::uniform_real_distribution<double> rv_marcum_q;
    thread_local std::mt19937 rv_marcum_q_gen;
    thread_local std::uint64_t seed_marcum_q{};
    rv_s_param = std::uniform_real_distribution<double>(lo1,hi1);
    seed_s_param = __rdtsc();
    rv_s_param_gen = std::mt19937(seed_s_param);
    rv_marcum_q = std::uniform_real_distribution<double>(lo2,hi2);
    seed_marcum_q = __rdtsc();
    rv_marcum_q_gen = std::mt19937(seed_marcum_q);
    const bool  is_arg_range_null = (arg_range==nullptr);
    if(!is_arg_range_null)
    {
        for(std::int32_t i{0}; i<n_func_args; ++i) 
        {
            const std::int32_t outer_idx = i*n_gauss_Q_vals;
            const double s = rv_s_param.operator()(rv_s_param_gen);
            const double ss = s*s;
            for(std::int32_t j{0}; j<n_gauss_Q_vals; ++j)   
            {
                const std::int32_t inner_idx = outer_idx+j;
                const double x = rv_marcum_q.operator()(rv_marcum_q_gen);
                const double sx{s*x};
                const double xspow2{std::fma(x,x,ss)};
                const double I0_val{std::cyl_bessel_i(0,sx)};
                const double exp_arg{0.5*xspow2};
                arg_range[inner_idx] = exp_arg;
            }
        }
    }
    else 
    {
        for(std::int32_t i{0}; i<n_func_args; ++i) 
        {
            const double s = rv_s_param.operator()(rv_s_param_gen);
            const double ss = s*s;
            for(std::int32_t j{0}; j<n_gauss_Q_vals; ++j)   
            {
                const double x = rv_marcum_q.operator()(rv_marcum_q_gen);
                const double sx{s*x};
                const double xspow2{std::fma(x,x,ss)};
                const double I0_val{std::cyl_bessel_i(0,sx)};
                const double exp_arg{0.5*xspow2};
                [[maybe_unused]] std::int32_t printf_ret = print_double("exp_arg",exp_arg,0);
            }
        }
    }
    return (0);
}