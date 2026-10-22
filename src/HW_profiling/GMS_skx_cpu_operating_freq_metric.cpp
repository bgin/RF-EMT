
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

#include <immintrin.h>
#include "GMS_skx_cpu_operating_freq_metric.h"
#include "GMS_skx_clk_hw_events_metrics.hpp"


__ATTR_ALWAYS_INLINE__
static inline 
std::int32_t 
convert_uint64_to_fp64_u8x(const std::valarray<std::uint64_t> &input,
                           std::valarray<double> &output) noexcept(true) 
{
    const std::uint64_t * __restrict__ p_in = reinterpret_cast<const std::uint64_t* __restrict__>(&input.operator[](0ull));
    double              * __restrict__ p_out= reinterpret_cast<double* __restrict__>(&output.operator[](0ull));
    std::size_t i,j;
    std::uint64_t inval0;
    std::uint64_t inval1;
    std::uint64_t inval2;
    std::uint64_t inval3;
    double        outval0;
    double        outval1;
    double        outval2;
    double        outval3;
    std::int32_t ret_status = -9999;
    const std::size_t in_size = input.size();
    if(__builtin_expect(in_size>=1ull,0) && 
       __builtin_expect(in_size<=8ull,0))
    {
        for(i = 0ull; i<in_size; ++i)
        {
            const std::uint64_t inval = p_in[i];
            const double outval       = static_cast<double>(inval);
            p_out[i]                  = outval;
        }
        ret_status = 1;
        return (ret_status);
    }
    else if(__builtin_expect(in_size>8ull,1))
    {
        for(i = 0ull; i!=ROUND_DOWN(in_size,7ull); i += 8ull)
        {
            const std::size_t idx1 = i+1ull;
            const std::size_t idx2 = i+2ull;
            const std::size_t idx3 = i+3ull;
            const std::size_t idx4 = i+4ull;
            const std::size_t idx5 = i+5ull;
            const std::size_t idx6 = i+6ull;
            const std::size_t idx7 = i+7ull;

            inval0  = p_in[i];
            outval0 = static_cast<double>(inval0);
            p_out[i]= outval0;
            inval1  = p_in[idx1];
            outval1 = static_cast<double>(inval1);
            p_out[idx1]= outval1;
            inval2  = p_in[idx2];
            outval2 = static_cast<double>(inval2);
            p_out[idx2]= outval2;
            inval3  = p_in[idx3];
            outval3 = static_cast<double>(inval3);
            p_out[idx3]= outval3;
            inval0  = p_in[idx4];
            outval0 = static_cast<double>(inval0);
            p_out[idx4]= outval0;
            inval1  = p_in[idx5];
            outval1 = static_cast<double>(inval1);
            p_out[idx5]= outval1;
            inval2  = p_in[idx6];
            outval2 = static_cast<double>(inval2);
            p_out[idx6]= outval2;
            inval3  = p_in[idx7];
            outval3 = static_cast<double>(inval3);
            p_out[idx7]= outval3;
        }

        for(j = i; j<in_size; ++j)   
        {
            inval0  = p_in[j];
            outval0 = static_cast<double>(inval0);
            p_out[j]= outval0;
        }
        ret_status = 2;
        return (ret_status);
    }
    /*Maybe unreachable!!*/
    return (ret_status);
}

__ATTR_ALWAYS_INLINE__
static inline
void 
vectorized_uint64_to_double_avx512(const std::uint64_t* __restrict src, double* __restrict dst, std::size_t size) {
    std::size_t i = 0ull;

    std::size_t dst_addr = reinterpret_cast<std::size_t>(dst);
    if (dst_addr % 64ull != 0 && size > 0) {
        std::size_t alignment_offset = (64ull - (dst_addr % 64ull)) / sizeof(double);
        std::size_t peel_limit = (alignment_offset < size) ? alignment_offset : size;
        
        for (; i < peel_limit; ++i) {
            dst[i] = static_cast<double>(src[i]);
        }
    }

    constexpr std::size_t UNROLL_FACTOR = 8ull;
    constexpr std::size_t ELEMENTS_PER_REG = 8ull;
    constexpr std::size_t LOOP_STRIDE = UNROLL_FACTOR * ELEMENTS_PER_REG; // 64 elements

    if (size >= i + LOOP_STRIDE) {
        std::size_t main_limit = size - LOOP_STRIDE;
        for (; i <= main_limit; i += LOOP_STRIDE) {
             std::size_t idx0  = i+0ull;
             std::size_t idx8  = i+8ull;
             std::size_t idx16 = i+16ull;
             std::size_t idx24 = i+24ull;
             std::size_t idx32 = i+32ull;
             std::size_t idx40 = i+40ull;
             std::size_t idx48 = i+48ull;
             std::size_t idx56 = i+56ull;
            // Load source chunks (using unaligned loads as source may not match destination alignment)
            __m512i src0 = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + idx0));
            __m512i src1 = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + idx8));
            __m512i src2 = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + idx16));
            __m512i src3 = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + idx24));
            __m512i src4 = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + idx32));
            __m512i src5 = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + idx40));
            __m512i src6 = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + idx48));
            __m512i src7 = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + idx56));

            __m512d dst0 = _mm512_cvtepu64_pd(src0);
            __m512d dst1 = _mm512_cvtepu64_pd(src1);
            __m512d dst2 = _mm512_cvtepu64_pd(src2);
            __m512d dst3 = _mm512_cvtepu64_pd(src3);
            __m512d dst4 = _mm512_cvtepu64_pd(src4);
            __m512d dst5 = _mm512_cvtepu64_pd(src5);
            __m512d dst6 = _mm512_cvtepu64_pd(src6);
            __m512d dst7 = _mm512_cvtepu64_pd(src7);

            _mm512_store_pd(dst + idx0, dst0);
            _mm512_store_pd(dst + idx8, dst1);
            _mm512_store_pd(dst + idx16, dst2);
            _mm512_store_pd(dst + idx24, dst3);
            _mm512_store_pd(dst + idx32, dst4);
            _mm512_store_pd(dst + idx40, dst5);
            _mm512_store_pd(dst + idx48, dst6);
            _mm512_store_pd(dst + idx56, dst7);
        }
    }

    for (; i + ELEMENTS_PER_REG <= size; i += ELEMENTS_PER_REG) {
        __m512i src_rem = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + i));
        __m512d dst_rem = _mm512_cvtepu64_pd(src_rem);
        _mm512_storeu_pd(dst + i, dst_rem); // Use unaligned store here as loop stride is smaller than alignment constraints
    }

    if (i < size) {
        std::size_t remaining = size - i;
        // Generate a tail mask based on the remaining elements
        __mmask8 mask = (__mmask8)((1U << remaining) - 1U);

        __m512i src_tail = _mm512_maskz_loadu_epi64(mask, src + i);
        __m512d dst_tail = _mm512_cvtepu64_pd(src_tail);
        _mm512_mask_storeu_pd(dst + i, mask, dst_tail);
    }
}



/*
 std::size_t           m_nruns;
    std::size_t           m_nsamples;
    std::size_t           m_totsamples;
    std::valarray<double> m_a; // a - event 
    std::valarray<double> m_b; // b - event
    std::valarray<double> m_c; // c - event
    std::valarray<double> m_d; // d - metric calculation result
    std::string           m_cpu_arch;
    std::string           m_date; // date of processing
    std::string           m_time; // time of processing
    double                m_a_robust_mean; // Mathematical-Expectation
    double                m_b_robust_mean; // Mathematical-Expectation
    double                m_c_robust_mean; // Mathematical-Expectation
    double                m_d_robust_mean; // Mathematical-Expectation
    double                m_a_normality; // Shapiro-Wilk normality test
    double                m_b_normality; // Shapiro-Wilk normality test
    double                m_c_normality; // Shapiro-Wilk normality test
    double                m_d_normality; // Shapiro-Wilk normality test
*/


gms::hw_profiling
::SKX_cpu_operating_frequency_t
::SKX_cpu_operating_frequency_t() noexcept(true) 
:
m_nruns{},
m_nsamples{},
m_totsamples{},
m_a(),
m_b(),
m_c(),
m_d(),
m_cpu_arch(""),
m_date(""),
m_time(""),
m_a_robust_mean{},
m_b_robust_mean{},
m_c_robust_mean{},
m_d_robust_mean{},
m_a_normality{},
m_b_normality{},
m_c_normality{},
m_d_normality{}
{

}

gms::hw_profiling
::SKX_cpu_operating_frequency_t
::SKX_cpu_operating_frequency_t(const std::size_t nruns,
                                const std::size_t nsamples,
                                const std::string &cpu_arch,
                                const std::string date,
                                const std::string time) noexcept(false)
:
m_nruns((nruns<=1ull)?8ull:nruns),
m_nsamples((nsamples<=1ull)?8ull:nsamples),
m_totsamples(m_nruns*m_nsamples),
m_a(m_totsamples),
m_b(m_totsamples),
m_c(m_totsamples),
m_d(m_totsamples),
m_cpu_arch(cpu_arch),
m_date(date),
m_time(time),
m_a_robust_mean{},
m_b_robust_mean{},
m_c_robust_mean{},
m_d_robust_mean{},
m_a_normality{},
m_b_normality{},
m_c_normality{},
m_d_normality{}
{

}

gms::hw_profiling
::SKX_cpu_operating_frequency_t
::SKX_cpu_operating_frequency_t(const SKX_cpu_operating_frequency_t &other) noexcept(false)
:
m_nruns(other.m_nruns),
m_nsamples(other.m_nsamples),
m_totsamples(other.m_totsamples),
m_a(other.m_a),
m_b(other.m_b),
m_c(other.m_c),
m_d(other.m_d),
m_cpu_arch(other.m_cpu_arch),
m_date(other.m_date),
m_time(other.m_time),
m_a_robust_mean(other.m_a_robust_mean),
m_b_robust_mean(other.m_b_robust_mean),
m_c_robust_mean(other.m_c_robust_mean),
m_d_robust_mean(other.m_d_robust_mean),
m_a_normality(other.m_a_normality),
m_b_normality(other.m_b_normality),
m_c_normality(other.m_c_normality),
m_d_normality(other.m_d_normality)
{

}

gms::hw_profiling
::SKX_cpu_operating_frequency_t
::SKX_cpu_operating_frequency_t(SKX_cpu_operating_frequency_t &&other) noexcept(true)
:
m_nruns(other.m_nruns),
m_nsamples(other.m_nsamples),
m_totsamples(other.m_totsamples),
m_a(std::move(other.m_a)),
m_b(std::move(other.m_b)),
m_c(std::move(other.m_c)),
m_d(std::move(other.m_d)),
m_cpu_arch(std::move(other.m_cpu_arch)),
m_date(std::move(other.m_date)),
m_time(std::move(other.m_time)),
m_a_robust_mean(other.m_a_robust_mean),
m_b_robust_mean(other.m_b_robust_mean),
m_c_robust_mean(other.m_c_robust_mean),
m_d_robust_mean(other.m_d_robust_mean),
m_a_normality(other.m_a_normality),
m_b_normality(other.m_b_normality),
m_c_normality(other.m_c_normality),
m_d_normality(other.m_d_normality)
{

}

gms::hw_profiling
::SKX_cpu_operating_frequency_t
::~SKX_cpu_operating_frequency_t()
{

}

gms::hw_profiling
::SKX_cpu_operating_frequency_t &
gms::hw_profiling
::SKX_cpu_operating_frequency_t
::operator=(const SKX_cpu_operating_frequency_t &other) noexcept(false) 
{
    if(__builtin_expect(this==&other,0)) { return (*this);}

    this->m_nruns      = other.m_nruns;
    this->m_nsamples   = other.m_nsamples;
    this->m_totsamples = other.m_totsamples; 
    this->m_a.operator=(other.m_a);
    this->m_b.operator=(other.m_b);
    this->m_c.operator=(other.m_c);
    this->m_d.operator=(other.m_d);
    this->m_cpu_arch.operator=(other.m_cpu_arch);
    this->m_date.operator=(other.m_date);
    this->m_time.operator=(other.m_time);
    this->m_a_robust_mean=other.m_a_robust_mean;
    this->m_b_robust_mean=other.m_b_robust_mean;
    this->m_c_robust_mean=other.m_c_robust_mean;
    this->m_d_robust_mean=other.m_d_robust_mean;
    this->m_a_normality  = other.m_a_normality;
    this->m_b_normality  = other.m_b_normality;
    this->m_c_normality  = other.m_c_normality;
    this->m_d_normality  = other.m_d_normality;

    return (*this);
}

gms::hw_profiling
::SKX_cpu_operating_frequency_t &
gms::hw_profiling
::SKX_cpu_operating_frequency_t
::operator=(SKX_cpu_operating_frequency_t &&other) noexcept(true) 
{
    if(__builtin_expect(this==&other,0)) { return (*this);}

    this->m_nruns      = other.m_nruns;
    this->m_nsamples   = other.m_nsamples;
    this->m_totsamples = other.m_totsamples; 
    this->m_a.operator=(std::move(other.m_a));
    this->m_b.operator=(std::move(other.m_b));
    this->m_c.operator=(std::move(other.m_c));
    this->m_d.operator=(std::move(other.m_d));
    this->m_cpu_arch.operator=(std::move(other.m_cpu_arch));
    this->m_date.operator=(std::move(other.m_date));
    this->m_time.operator=(std::move(other.m_time));
    this->m_a_robust_mean=other.m_a_robust_mean;
    this->m_b_robust_mean=other.m_b_robust_mean;
    this->m_c_robust_mean=other.m_c_robust_mean;
    this->m_d_robust_mean=other.m_d_robust_mean;
    this->m_a_normality  = other.m_a_normality;
    this->m_b_normality  = other.m_b_normality;
    this->m_c_normality  = other.m_c_normality;
    this->m_d_normality  = other.m_d_normality;

    return (*this);
}

auto 
gms::hw_profiling
::SKX_cpu_operating_frequency_t
::operator==(const SKX_cpu_operating_frequency_t &rhs)->std::valarray<bool>
{
    if(__builtin_expect(this->m_totsamples!=rhs.m_totsamples,0)) { return (std::valarray<bool>());}
    std::valarray<bool> vequal = this->m_d == rhs.m_d;
    return (vequal);
}

auto 
gms::hw_profiling
::SKX_cpu_operating_frequency_t
::operator!=(const SKX_cpu_operating_frequency_t &rhs)->std::valarray<bool> 
{
    if(__builtin_expect(this->m_totsamples!=rhs.m_totsamples,0)) { return (std::valarray<bool>());}
    std::valarray<bool> vnequal = this->m_d != rhs.m_d;
    return (vnequal);
}

std::int32_t 
gms::hw_profiling
::SKX_cpu_operating_frequency_t
::convert_m_a_to_fp64_u8x(const std::valarray<std::uint64_t> &ui_a) noexcept(true)
{   
    std::int32_t ret_status = -1;
    if(__builtin_expect(this->m_a.size()!=ui_a.size(),0))
    {
        ret_status = 0;
        return (ret_status);
    }
    ret_status = convert_uint64_to_fp64_u8x(ui_a,this->m_a);
    return (ret_status);
}

std::int32_t 
gms::hw_profiling
::SKX_cpu_operating_frequency_t
::convert_m_b_to_fp64_u8x(const std::valarray<std::uint64_t> &ui_b) noexcept(true)
{   
    std::int32_t ret_status = -1;
    if(__builtin_expect(this->m_b.size()!=ui_b.size(),0))
    {
        ret_status = 0;
        return (ret_status);
    }
    ret_status = convert_uint64_to_fp64_u8x(ui_b,this->m_b);
    return (ret_status);
}

std::int32_t 
gms::hw_profiling
::SKX_cpu_operating_frequency_t
::convert_m_c_to_fp64_u8x(const std::valarray<std::uint64_t> &ui_c) noexcept(true)
{   
    std::int32_t ret_status = -1;
    if(__builtin_expect(this->m_c.size()!=ui_c.size(),0))
    {
        ret_status = 0;
        return (ret_status);
    }
    ret_status = convert_uint64_to_fp64_u8x(ui_c,this->m_c);
    return (ret_status);
}

std::int32_t 
gms::hw_profiling
::SKX_cpu_operating_frequency_t
::convert_m_a_to_fp64_u8x_avx512(const std::valarray<std::uint64_t> &ui_a) noexcept(true)
{
    std::int32_t ret_status = -1;
    if(__builtin_expect(this->m_a.size()!=ui_a.size(),0))
    {
        ret_status = 0;
        return (ret_status);
    }
    const std::uint64_t * __restrict__ p_ui_a = reinterpret_cast<const std::uint64_t* __restrict__>(&ui_a.operator[](0ull));
    double              * __restrict__ p_a    = reinterpret_cast<double * __restrict__>(&this->m_a.operator[](0ull));
    vectorized_uint64_to_double_avx512(&p_ui_a[0],&p_a[0],this->m_a.size());
    ret_status = 1;
    return (ret_status);
}

std::int32_t 
gms::hw_profiling
::SKX_cpu_operating_frequency_t
::convert_m_b_to_fp64_u8x_avx512(const std::valarray<std::uint64_t> &ui_b) noexcept(true)
{
    std::int32_t ret_status = -1;
    if(__builtin_expect(this->m_b.size()!=ui_b.size(),0))
    {
        ret_status = 0;
        return (ret_status);
    }
    const std::uint64_t * __restrict__ p_ui_b = reinterpret_cast<const std::uint64_t* __restrict__>(&ui_b.operator[](0ull));
    double              * __restrict__ p_b    = reinterpret_cast<double * __restrict__>(&this->m_b.operator[](0ull));
    vectorized_uint64_to_double_avx512(&p_ui_b[0],&p_b[0],this->m_b.size());
    ret_status = 1;
    return (ret_status);
}

std::int32_t 
gms::hw_profiling
::SKX_cpu_operating_frequency_t
::convert_m_c_to_fp64_u8x_avx512(const std::valarray<std::uint64_t> &ui_c) noexcept(true)
{
    std::int32_t ret_status = -1;
    if(__builtin_expect(this->m_c.size()!=ui_c.size(),0))
    {
        ret_status = 0;
        return (ret_status);
    }
    const std::uint64_t * __restrict__ p_ui_c = reinterpret_cast<const std::uint64_t* __restrict__>(&ui_c.operator[](0ull));
    double              * __restrict__ p_c    = reinterpret_cast<double * __restrict__>(&this->m_c.operator[](0ull));
    vectorized_uint64_to_double_avx512(&p_ui_c[0],&p_c[0],this->m_c.size());
    ret_status = 1;
    return (ret_status);
}

std::int32_t 
gms::hw_profiling
::SKX_cpu_operating_frequency_t
::compute_metric_u4x() noexcept(true) 
{
    std::size_t i,j;
    double a0;
    double a1;
    double a2;
    double a3;
    double b0;
    double b1;
    double b2;
    double b3;
    double c0;
    double c1;
    double c2;
    double c3;
    std::int32_t ret_status = -9999;
    if(__builtin_expect(this->m_totsamples>=1ull,0) && 
       __builtin_expect(this->m_totsamples<=4ull,0))
    {
        for(i = 0ull; i<this->m_totsamples; ++i) 
        {
            a0 = this->m_a.operator[](i);
            b0 = this->m_b.operator[](i);
            c0 = this->m_c.operator[](i);
            this->m_d.operator[](i) = skx_cpu_operating_freq(a0,b0,c0);
        }
        ret_status = 1;
        return (ret_status);
    }
    else if(__builtin_expect(this->m_totsamples>4ull,1))
    {
        for(i = 0ull; i!=ROUND_DOWN(this->m_totsamples,3ull); i += 4ull) 
        {
            std::size_t idx1 = i+1ull;
            std::size_t idx2 = i+2ull;
            std::size_t idx3 = i+3ull;
            a0 = this->m_a.operator[](i);
            b0 = this->m_b.operator[](i);
            c0 = this->m_c.operator[](i);
            this->m_d.operator[](i) = skx_cpu_operating_freq(a0,b0,c0);
            a1 = this->m_a.operator[](idx1);
            b1 = this->m_b.operator[](idx1);
            c1 = this->m_c.operator[](idx1);
            this->m_d.operator[](idx1) = skx_cpu_operating_freq(a1,b1,c1);
            a2 = this->m_a.operator[](idx2);
            b2 = this->m_b.operator[](idx2);
            c2 = this->m_c.operator[](idx2);
            this->m_d.operator[](idx2) = skx_cpu_operating_freq(a2,b2,c2);
            a3 = this->m_a.operator[](idx3);
            b3 = this->m_b.operator[](idx3);
            c3 = this->m_c.operator[](idx3);
            this->m_d.operator[](idx3) = skx_cpu_operating_freq(a3,b3,c3);
        }

        for(j = i; j<this->m_totsamples; ++j)  
        {
            a0 = this->m_a.operator[](j);
            b0 = this->m_b.operator[](j);
            c0 = this->m_c.operator[](j);
            this->m_d.operator[](j) = skx_cpu_operating_freq(a0,b0,c0);
        }
        ret_status = 2;
        return (ret_status);
    }
    return (ret_status);
}









