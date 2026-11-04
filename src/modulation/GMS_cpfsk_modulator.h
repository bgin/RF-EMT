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

#ifndef __GMS_CPFSK_MODULATOR_H__
#define __GMS_CPFSK_MODULATOR_H__ 231020260156

#include <valarray>
#include <vector>
#include <cmath>
#include <numeric>
#include <cassert>
#include <cstdint>
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_CPFSK_MODULATOR_MAJOR = 1;
     static const unsigned int GMS_CPFSK_MODULATOR_MINOR = 0;
     static const unsigned int GMS_CPFSK_MODULATOR_MICRO = 0;
     static const unsigned int GMS_CPFSK_MODULATOR_FULLVER =
       1000U*GMS_CPFSK_MODULATOR_MAJOR+100U*GMS_CPFSK_MODULATOR_MINOR+
       10U*GMS_CPFSK_MODULATOR_MICRO;
     static const char GMS_CPFSK_MODULATOR_CREATION_DATE[] = "23-10-2026 01:56PM +00200 (FRI 23 OCT 2026 GMT+2)";
     static const char GMS_CPFSK_MODULATOR_BUILD_DATE[]    = __DATE__; 
     static const char GMS_CPFSK_MODULATOR_BUILD_TIME[]    = __TIME__;
     static const char GMS_CPFSK_MODULATOR_SYNOPSIS[]      = "CPFSK Modulator std::valarray based.";

}


namespace gms 
{

namespace radiolocation 
{

template<typename T>
struct alignas(64) cpfsk_modulator_t 
{
    std::size_t      m_nsymbols;
    T                m_sample_rate;
    T                m_carrier_freq;
    T                m_symbol_rate;
    std::size_t      m_samples_per_symbol;
    T                m_modulation_index;
    T                m_amplitude;
    T                m_tstep;
    T                m_phase_scale_factor;
    std::valarray<T> m_time_domain_signal;
    
    __ATTR_ALWAYS_INLINE__
    inline 
    cpfsk_modulator_t() noexcept(true)
    :
    m_nsymbols{},
    m_sample_rate{},
    m_carrier_freq{},
    m_symbol_rate{},
    m_samples_per_symbol{},
    m_modulation_index{},
    m_amplitude{},
    m_tstep{},
    m_phase_scale_factor{},
    m_time_domain_signal()
    {

    }

    __ATTR_ALWAYS_INLINE__
    inline 
    cpfsk_modulator_t(const std::size_t nsymbols,
                      const T sample_rate,
                      const T carrier_freq,
                      const T symbol_rate,
                      const T modulation_index,
                      const T amplitude) noexcept(false)
    :
    m_nsymbols(nsymbols),
    m_sample_rate(sample_rate),
    m_carrier_freq(carrier_freq),
    m_symbol_rate(symbol_rate),
    m_samples_per_symbol(static_cast<std::size_t>(std::round(m_sample_rate/m_symbol_rate))),
    m_modulation_index(modulation_index),
    m_amplitude(amplitude),
    m_tstep(static_cast<T>(1.0)/m_sample_rate),
    m_phase_scale_factor(static_cast<T>(3.1415926535897932384626433832795)*m_modulation_index*m_symbol_rate),
    m_time_domain_signal(m_nsymbols*m_samples_per_symbol)
    {

    }

    __ATTR_ALWAYS_INLINE__
    inline 
    cpfsk_modulator_t(const cpfsk_modulator_t &other) noexcept(false) 
    :
    m_nsymbols(other.m_nsymbols),
    m_sample_rate(other.m_sample_rate),
    m_carrier_freq(other.m_carrier_freq),
    m_symbol_rate(other.m_symbol_rate),
    m_samples_per_symbol(other.m_samples_per_symbol),
    m_modulation_index(other.m_modulation_index),
    m_amplitude(other.m_amplitude),
    m_tstep(other.m_tstep),
    m_phase_scale_factor(other.m_phase_scale_factor),
    m_time_domain_signal(other.m_time_domain_signal)
    {

    }

     __ATTR_ALWAYS_INLINE__
    inline 
    cpfsk_modulator_t(cpfsk_modulator_t &&other) noexcept(true) 
    :
    m_nsymbols(other.m_nsymbols),
    m_sample_rate(other.m_sample_rate),
    m_carrier_freq(other.m_carrier_freq),
    m_symbol_rate(other.m_symbol_rate),
    m_samples_per_symbol(other.m_samples_per_symbol),
    m_modulation_index(other.m_modulation_index),
    m_amplitude(other.m_amplitude),
    m_tstep(other.m_tstep),
    m_phase_scale_factor(other.m_phase_scale_factor),
    m_time_domain_signal(std::move(other.m_time_domain_signal))
    {

    }

     __ATTR_ALWAYS_INLINE__
    inline
    ~cpfsk_modulator_t() noexcept(true)
    {

    }

     __ATTR_ALWAYS_INLINE__
    inline
    cpfsk_modulator_t & 
    operator=(const cpfsk_modulator_t &other) noexcept(false) 
    {
        if(__builtin_expect(this==&other,0)) { return (*this);}
        this->m_nsymbols           = other.m_nsymbols;
        this->m_sample_rate        = other.m_sample_rate;
        this->m_carrier_freq       = other.m_carrier_freq;
        this->m_symbol_rate        = other.m_symbol_rate;
        this->m_samples_per_symbol = other.m_samples_per_symbol;
        this->m_modulation_index   = other.m_modulation_index;
        this->m_amplitude          = other.m_amplitude;
        this->m_tstep              = other.m_tstep;
        this->m_phase_scale_factor = other.m_phase_scale_factor;
        this->m_time_domain_signal.operator=(other.m_time_domain_signal);
        return (*this);
    }

     __ATTR_ALWAYS_INLINE__
    inline
    cpfsk_modulator_t & 
    operator=(cpfsk_modulator_t &&other) noexcept(true) 
    {
        if(__builtin_expect(this==&other,0)) { return (*this);}
        this->m_nsymbols           = other.m_nsymbols;
        this->m_sample_rate        = other.m_sample_rate;
        this->m_carrier_freq       = other.m_carrier_freq;
        this->m_symbol_rate        = other.m_symbol_rate;
        this->m_samples_per_symbol = other.m_samples_per_symbol;
        this->m_modulation_index   = other.m_modulation_index;
        this->m_amplitude          = other.m_amplitude;
        this->m_tstep              = other.m_tstep;
        this->m_phase_scale_factor = other.m_phase_scale_factor;
        this->m_time_domain_signal.operator=(std::move(other.m_time_domain_signal));
        return (*this);
    }

     __ATTR_ALWAYS_INLINE__
    inline
    std::int32_t modulate(const std::valarray<T> &symbols) noexcept(false) 
    {
        std::int32_t ret_status = -9999;
        if(__builtin_expect(this->m_nsymbols!=symbols.size(),0))
        {
            ret_status = -1;
            return (ret_status);
        }
        
        std::valarray<T> slice_t(this->m_samples_per_symbol);
        T phase_memory{static_cast<T>(0.0)};
        for(std::size_t i = 0ull; i<this->m_samples_per_symbol; ++i) 
        {
            const T i_val     = static_cast<T>(i);
            const T tstep_val = i_val*this->m_tstep;
            slice_t[i]        = tstep_val;
        }
        constexpr T TwoPI = static_cast<T>(6.283185307179586476925286766559);
        const std::valarray<T> carrier_phase_slice = slice_t*(TwoPI*this->m_carrier_freq);
        std::valarray<T> current_phase_trajectory(this->m_samples_per_symbol);
        std::valarray<T> total_phase(this->m_samples_per_symbol);
        std::valarray<T> symbol_waveform(this->m_samples_per_symbol);
        constexpr T PI = static_cast<T>(3.1415926535897932384626433832795);
        
        for(std::size_t k = 0ull; k<this->m_nsymbols; ++k) 
        {
            const T current_symbol   = symbols[k];
            current_phase_trajectory = slice_t*(this->m_phase_scale_factor*current_symbol);
            total_phase              = carrier_phase_slice+current_phase_trajectory+phase_memory;
            symbol_waveform          = this->m_amplitude * std::cos(total_phase);
            this->m_time_domain_signal[std::slice(k*m_samples_per_symbol,m_samples_per_symbol,1)] = symbol_waveform;
            T phase_accumulated      = PI*this->m_modulation_index*current_symbol;
            T carrier_advance        = TwoPI*this->m_carrier_freq*(static_cast<T>(this->m_samples_per_symbol)*this->m_tstep);
            phase_memory             += (carrier_advance+phase_accumulated);
            phase_memory             = std::remainder(phase_memory,TwoPI);
        }
        ret_status = 0;
        return (ret_status);
    }


    
};

}

}

#endif /*__GMS_CPFSK_MODULATOR_H__*/