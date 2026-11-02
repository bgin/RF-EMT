
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

#ifndef __GMS_AMPLITUDE_MODULATOR_H__
#define __GMS_AMPLITUDE_MODULATOR_H__ 281020260635

#include <random>
#include "GMS_config.h"

namespace file_info 
{

     const unsigned int GMS_AMPLITUDE_MODULATOR_MAJOR = 1;
     const unsigned int GMS_AMPLITUDE_MODULATOR_MINOR = 0;
     const unsigned int GMS_AMPLITUDE_MODULATOR_MICRO = 0;
     const unsigned int GMS_AMPLITUDE_MODULATOR_FULLVER =
       1000U*GMS_AMPLITUDE_MODULATOR_MAJOR+100U*GMS_AMPLITUDE_MODULATOR_MINOR+
       10U*GMS_AMPLITUDE_MODULATOR_MICRO;
     const char * const GMS_AMPLITUDE_MODULATOR_CREATION_DATE = "28-10-2026 06:35 AM +00200 (WED 28 OCT 2026 GMT+2)";
     const char * const GMS_AMPLITUDE_MODULATOR_BUILD_DATE    = __DATE__ " " __TIME__;
     static const char  GMS_AMPLITUDE_MODULATOR_AUTHOR[]      = "Programmer: Bernard Gingold, beniekg@gmail.com"; 
     const char * const GMS_AMPLITUDE_MODULATOR_SYNOPSIS      = "Basic AM modulator classes (float and double).";

}

namespace gms 
{

namespace radiolocation 
{


struct alignas(64) amplitude_modulator_f32_t 
{
    
    float m_carrier_freq;
    float m_sampling_rate;
    float m_time;
    
    thread_local std::mt19937 m_rng;
    thread_local std::uniform_real_distribution<float> m_dist;

    amplitude_modulator_f32_t() = default;

    amplitude_modulator_f32_t(const float,const float);

    amplitude_modulator_f32_t(const amplitude_modulator_f32_t &);

    amplitude_modulator_f32_t( amplitude_modulator_f32_t &&) noexcept(true);
    
    ~amplitude_modulator_f32_t();

    amplitude_modulator_f32_t &
             operator=(const amplitude_modulator_f32_t &);

    amplitude_modulator_f32_t &
             operator=(amplitude_modulator_f32_t &&) noexcept(true);

    float generate_next_sample(const float);

    float generate_random_message_value();

   
};

struct alignas(64) amplitude_modulator_f64_t 
{
    
    double m_carrier_freq;
    double m_sampling_rate;
    double m_time;
    
    thread_local std::mt19937 m_rng;
    thread_local std::uniform_real_distribution<double> m_dist;

    amplitude_modulator_f64_t() = default;

    amplitude_modulator_f64_t(const double,const double);

    amplitude_modulator_f64_t(const amplitude_modulator_f64_t &);

    amplitude_modulator_f64_t(amplitude_modulator_f64_t &&) noexcept(true);

    ~amplitude_modulator_f64_t();

    amplitude_modulator_f64_t &
             operator=(amplitude_modulator_f64_t &);

    amplitude_modulator_f64_t &
             operator=(amplitude_modulator_f64_t &&) noexcept(true);

    double generate_next_sample(const double);

    double generate_random_message_value();

   
};

}

}

#endif // __GMS_AMPLITUDE_MODULATOR_H__
