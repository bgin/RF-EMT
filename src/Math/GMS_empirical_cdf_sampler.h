
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

#ifndef __GMS_EMPIRICAL_CDF_SAMPLER_H__
#define __GMS_EMPIRICAL_CDF_SAMPLER_H__ 061120260841

#include <iostream>
#include <valarray>
#include <algorithm>
#include <random>
#include <numeric>
#include "GMS_config.h"

namespace file_info 
{

     const unsigned int GMS_EMPIRICAL_CDF_SAMPLER_MAJOR = 1;
     const unsigned int GMS_EMPIRICAL_CDF_SAMPLER_MINOR = 0;
     const unsigned int GMS_EMPIRICAL_CDF_SAMPLER_MICRO = 0;
     const unsigned int GMS_EMPIRICAL_CDF_SAMPLER_FULLVER =
       1000U*GMS_EMPIRICAL_CDF_SAMPLER_MAJOR+100U*GMS_EMPIRICAL_CDF_SAMPLER_MINOR+
       10U*GMS_EMPIRICAL_CDF_SAMPLER_MICRO;
     const char * const GMS_EMPIRICAL_CDF_SAMPLER_CREATION_DATE = "06-11-2026 08:43 AM +00200 (FRI 05 NOV 2026 GMT+2)";
     const char * const GMS_EMPIRICAL_CDF_SAMPLER_BUILD_DATE    = __DATE__ " " __TIME__;
     static const char  GMS_EMPIRICAL_CDF_SAMPLER_AUTHOR[]      = "Programmer: Bernard Gingold, beniekg@gmail.com"; 
     const char * const GMS_EMPIRICAL_CDF_SAMPLER_SYNOPSIS      = "Empirical CDF Inverse Sampler.";

}

namespace gms 
{

namespace math 
{

template<typename T>
struct alignas(64) ecdf_inverse_sampler_t
{
    std::size_t      m_num_samples;
    std::valarray<T> m_samples;
    std::valarray<T> m_cdf_weights;

    ecdf_inverse_sampler_t() noexcept(true)
    :
    m_num_samples{},
    m_samples(),
    m_cdf_weights()
    {

    }

    ecdf_inverse_sampler_t(const std::size_t num_samples) noexcept(false)
    :
    m_num_samples((num_samples<=1ull)?2ull:num_samples),
    m_samples(m_num_samples),
    m_cdf_weights(m_num_samples)
    {

    }

    ecdf_inverse_sampler_t(const std::valarray<T> &data) noexcept(false)
    :
    m_num_samples(data.size()),
    m_samples(data),
    m_cdf_weights(m_samples.size())
    {

    }

    ecdf_inverse_sampler_t(std::valarray<T> &&data) noexcept(true)
    :
    m_num_samples(data.size()),
    m_samples(data),
    m_cdf_weights(m_samples.size())
    {

    }

    ~ecdf_inverse_sampler_t() noexcept(true)
    {

    }

    ecdf_inverse_sampler_t & 
    operator=(const ecdf_inverse_sampler_t &other) noexcept(false) 
    {
         if(__builtin_expect(this==&other)) { return (*this);}
         this->m_num_samples = other.m_num_samples;
         this->m_samples.operator=(other.m_samples);
         this->m_cdf_weights.operator=(other.m_cdf_weights);
         return (*this);
    }

    ecdf_inverse_sampler_t & 
    operator=(ecdf_inverse_sampler_t &&other) noexcept(true) 
    {
         if(__builtin_expect(this==&other)) { return (*this);}
         this->m_num_samples = other.m_num_samples;
         this->m_samples.operator=(other.m_samples);
         this->m_cdf_weights.operator=(other.m_cdf_weights);
         return (*this);
    }

    bool compute_weights()
    {
        bool ret_status = false;
        std::sort(this->m_samples.begin(),this->m_samples.end());
        const T inv_nsamples = 1.0/static_cast<T>(this->m_num_samples);
        for(std::size_t i = 1ull; i<=this->m_num_samples; ++i)
        {
            this->m_cdf_weights.operator[](i) = static_cast<T>(i)*inv_nsamples;
        }
        ret_status = true;
        return (ret_status);
    }

    template<typename URNG>
    T sample(URNG &g)
    {
        std::uniform_real_distribution<T> dist(static_cast<T>(0.0), static_cast<T>(1.0));
        T u = dist(g);
        auto it = std::lower_bound(this->m_cdf_weights.begin(), this->m_cdf_weights.end(),u);
        size_t index = std::distance(this->m_cdf_weights.begin(), it);
        if (index >= this->m_samples.size()) {
            index = this->m_samples.size() - 1ull;
        }
        return this->m_samples[index];
    }
};

}

}

#endif /*__GMS_EMPIRICAL_CDF_SAMPLER__*/