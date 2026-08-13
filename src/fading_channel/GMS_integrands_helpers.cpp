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