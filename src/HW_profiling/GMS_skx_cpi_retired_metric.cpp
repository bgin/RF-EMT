
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
#include <cpuid.h>
#if (SKX_CPI_RETIRED_METRIC_VERBOSE) == 1
#include <iostream>
#include <iomanip>
#endif
#include "GMS_skx_cpi_retired_metric.h"
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

    // 1. PEEL LOOP: Align the destination pointer to 64 bytes to guarantee optimal streaming stores
    // Aligning 'dst' is generally preferred since write-combining and store ports are highly sensitive to cache-line splits.
    std::size_t dst_addr = reinterpret_cast<std::size_t>(dst);
    if (dst_addr % 64ull != 0 && size > 0) {
        std::size_t alignment_offset = (64ull - (dst_addr % 64ull)) / sizeof(double);
        std::size_t peel_limit = (alignment_offset < size) ? alignment_offset : size;
        
        for (; i < peel_limit; ++i) {
            dst[i] = static_cast<double>(src[i]);
        }
    }

    // 2. MAIN UNROLLED KERNEL: Process 64 elements per loop iteration
    // Each ZMM register holds 8 uint64_t/double elements. Unrolling 8x uses ZMM0-ZMM7.
    // This depth hides the ~4-5 cycle latency of _mm512_cvtepu64_pd.
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

            // Native hardware translation: uint64_t -> double via vcvtuqq2pd
            __m512d dst0 = _mm512_cvtepu64_pd(src0);
            __m512d dst1 = _mm512_cvtepu64_pd(src1);
            __m512d dst2 = _mm512_cvtepu64_pd(src2);
            __m512d dst3 = _mm512_cvtepu64_pd(src3);
            __m512d dst4 = _mm512_cvtepu64_pd(src4);
            __m512d dst5 = _mm512_cvtepu64_pd(src5);
            __m512d dst6 = _mm512_cvtepu64_pd(src6);
            __m512d dst7 = _mm512_cvtepu64_pd(src7);

            // Aligned streaming or standard stores (since dst is forced to 64-byte alignment via peel)
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

    // 3. VECTOR REMAINDER LOOP: Convert chunks of 8 elements at a time
    for (; i + ELEMENTS_PER_REG <= size; i += ELEMENTS_PER_REG) {
        __m512i src_rem = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + i));
        __m512d dst_rem = _mm512_cvtepu64_pd(src_rem);
        _mm512_storeu_pd(dst + i, dst_rem); // Use unaligned store here as loop stride is smaller than alignment constraints
    }

    // 4. SCALAR MASK / TAIL CLEANUP: Handle remaining < 8 elements
    if (i < size) {
        std::size_t remaining = size - i;
        // Generate a tail mask based on the remaining elements
        __mmask8 mask = (__mmask8)((1U << remaining) - 1U);

        // Masked load, conversion, and store to completely avoid faults and scalar fallback branches
        __m512i src_tail = _mm512_maskz_loadu_epi64(mask, src + i);
        __m512d dst_tail = _mm512_cvtepu64_pd(src_tail);
        _mm512_mask_storeu_pd(dst + i, mask, dst_tail);
    }
}


gms::hw_profiling
::SKX_cpu_cpi_retired_t
::SKX_cpu_cpi_retired_t() noexcept(true)
:
m_nruns{},
m_nsamples{},
m_totsamples{},
m_a(),
m_b(),
m_c(),
m_date(""),
m_time(""),
m_a_robust_mean{},
m_b_robust_mean{},
m_c_robust_mean{},
m_a_normality{},
m_b_normality{},
m_c_normality{}
{
    
}

gms::hw_profiling
::SKX_cpu_cpi_retired_t
::SKX_cpu_cpi_retired_t(const std::size_t nruns,
                             const std::size_t nsamples,
                             const std::string date,
                             const std::string time,
                             const std::uint32_t * __restrict__ p_cpuid_basic) noexcept(false)
:
m_nruns((nruns<=1ull)?8ull:nruns),
m_nsamples((nsamples<=1ull)?8ull:nsamples),
m_totsamples(m_nruns*m_nsamples),
m_a(m_totsamples),
m_b(m_totsamples),
m_c(m_totsamples),
m_date(date),
m_time(time),
m_a_robust_mean{},
m_b_robust_mean{},
m_c_robust_mean{},
m_a_normality{},
m_b_normality{},
m_c_normality{}
{
   if(nullptr!=p_cpuid_basic)
   {
      _mm256_storeu_si256((__m256i*)&this->m_cpuid_basic[0],
                           _mm256_load_si256((__m256i*)&p_cpuid_basic[0]));
   }
}

gms::hw_profiling
::SKX_cpu_cpi_retired_t
::SKX_cpu_cpi_retired_t(const SKX_cpu_cpi_retired_t &other) noexcept(false)
:
m_nruns(other.m_nruns),
m_nsamples(other.m_nsamples),
m_totsamples(other.m_totsamples),
m_a(other.m_a),
m_b(other.m_b),
m_c(other.m_c),
m_date(other.m_date),
m_time(other.m_time),
m_a_robust_mean(other.m_a_robust_mean),
m_b_robust_mean(other.m_b_robust_mean),
m_c_robust_mean(other.m_c_robust_mean),
m_a_normality(other.m_a_normality),
m_b_normality(other.m_b_normality),
m_c_normality(other.m_c_normality)
{
    _mm256_storeu_si256((__m256i*)&this->m_cpuid_basic[0],
                           _mm256_loadu_si256((__m256i*)&other.m_cpuid_basic[0]));
}

gms::hw_profiling
::SKX_cpu_cpi_retired_t
::SKX_cpu_cpi_retired_t(SKX_cpu_cpi_retired_t &&other) noexcept(true)
:
m_nruns(other.m_nruns),
m_nsamples(other.m_nsamples),
m_totsamples(other.m_totsamples),
m_a(std::move(other.m_a)),
m_b(std::move(other.m_b)),
m_c(std::move(other.m_c)),
m_date(std::move(other.m_date)),
m_time(std::move(other.m_time)),
m_a_robust_mean(other.m_a_robust_mean),
m_b_robust_mean(other.m_b_robust_mean),
m_c_robust_mean(other.m_c_robust_mean),
m_a_normality(other.m_a_normality),
m_b_normality(other.m_b_normality),
m_c_normality(other.m_c_normality)
{
    _mm256_storeu_si256((__m256i*)&this->m_cpuid_basic[0],
                           _mm256_loadu_si256((__m256i*)&other.m_cpuid_basic[0]));
}

gms::hw_profiling
::SKX_cpu_cpi_retired_t
::~SKX_cpu_cpi_retired_t()
{

}

gms::hw_profiling
::SKX_cpu_cpi_retired_t &
gms::hw_profiling
::SKX_cpu_cpi_retired_t
::operator=(const SKX_cpu_cpi_retired_t &other) noexcept(false)
{
    if(__builtin_expect(this==&other,0)) { return (*this);}

    this->m_nruns      = other.m_nruns;
    this->m_nsamples   = other.m_nsamples;
    this->m_totsamples = other.m_totsamples; 
    this->m_a.operator=(other.m_a);
    this->m_b.operator=(other.m_b);
    this->m_c.operator=(other.m_c);
    this->m_date.operator=(other.m_date);
    this->m_time.operator=(other.m_time);
    this->m_a_robust_mean=other.m_a_robust_mean;
    this->m_b_robust_mean=other.m_b_robust_mean;
    this->m_c_robust_mean=other.m_c_robust_mean;
    this->m_a_normality  = other.m_a_normality;
    this->m_b_normality  = other.m_b_normality;
    this->m_c_normality  = other.m_c_normality;
    _mm256_storeu_si256((__m256i*)&this->m_cpuid_basic[0],
                           _mm256_loadu_si256((__m256i*)&other.m_cpuid_basic[0]));
    return (*this);
}

gms::hw_profiling
::SKX_cpu_cpi_retired_t &
gms::hw_profiling
::SKX_cpu_cpi_retired_t
::operator=(SKX_cpu_cpi_retired_t &&other) noexcept(true)
{
    if(__builtin_expect(this==&other,0)) { return (*this);}

    this->m_nruns      = other.m_nruns;
    this->m_nsamples   = other.m_nsamples;
    this->m_totsamples = other.m_totsamples; 
    this->m_a.operator=(std::move(other.m_a));
    this->m_b.operator=(std::move(other.m_b));
    this->m_c.operator=(std::move(other.m_c));
    this->m_date.operator=(std::move(other.m_date));
    this->m_time.operator=(std::move(other.m_time));
    this->m_a_robust_mean=other.m_a_robust_mean;
    this->m_b_robust_mean=other.m_b_robust_mean;
    this->m_c_robust_mean=other.m_c_robust_mean;
    this->m_a_normality  = other.m_a_normality;
    this->m_b_normality  = other.m_b_normality;
    this->m_c_normality  = other.m_c_normality;
    _mm256_storeu_si256((__m256i*)&this->m_cpuid_basic[0],
                           _mm256_loadu_si256((__m256i*)&other.m_cpuid_basic[0]));
    return (*this);
}

auto 
gms::hw_profiling
::SKX_cpu_cpi_retired_t
::operator==(const SKX_cpu_cpi_retired_t &rhs)->std::valarray<bool>
{
    if(__builtin_expect(this->m_totsamples!=rhs.m_totsamples,0)) { return (std::valarray<bool>());}
    std::valarray<bool> vequal = this->m_c == rhs.m_c;
    return (vequal);
}

auto 
gms::hw_profiling
::SKX_cpu_cpi_retired_t
::operator!=(const SKX_cpu_cpi_retired_t &rhs)->std::valarray<bool> 
{
    if(__builtin_expect(this->m_totsamples!=rhs.m_totsamples,0)) { return (std::valarray<bool>());}
    std::valarray<bool> vnequal = this->m_c != rhs.m_c;
    return (vnequal);
}

std::int32_t 
gms::hw_profiling
::SKX_cpu_cpi_retired_t
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
::SKX_cpu_cpi_retired_t
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
::SKX_cpu_cpi_retired_t
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
::SKX_cpu_cpi_retired_t
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
::SKX_cpu_cpi_retired_t
::parse_cpuid_basic_info(const bool is_this_cpu_profiled)
{
    std::int32_t ret_status = -9999;
    if(is_this_cpu_profiled == true)
    {
        std::uint32_t registers[4] = {0};
        ret_status = __get_cpuid(1, &registers[0], &registers[1], &registers[2], &registers[3]);
        if(0==ret_status) {return (ret_status);}
        uint32_t eax = registers[0];

        uint32_t base_family     = (eax >> 8) & 0xF;
        uint32_t base_model      = (eax >> 4) & 0xF;
        uint32_t stepping        = eax & 0xF;
        uint32_t extended_model  = (eax >> 16) & 0xF;
        uint32_t extended_family = (eax >> 20) & 0xFF;
    // If base_family is 0xF, Display Family = base_family + extended_family
        uint32_t display_family = (base_family == 0xF) ? (base_family + extended_family) : base_family;
    // For Intel processors, if display_family is 6 or 15, 
    // Display Model = (extended_model << 4) | base_model
        uint32_t display_model = base_model;
        if (display_family == 6 || display_family == 15) 
        {
             display_model = (extended_model << 4) | base_model;
        }

#if (SKX_CPI_RETIRED_METRIC_VERBOSE) == 1
       std::cout << "====================================================\n";
       std::cout << "      INTEL CPUID DECODER (MICROARCHITECTURE)       \n";
       std::cout << "====================================================\n";
       std::cout << "Raw EAX Signature: 0x" << std::hex << std::uppercase << std::setw(8) << std::setfill('0') << eax << std::dec << "\n\n";
    
       std::cout << "Parsed Raw Fields:\n";
       std::cout << "  - Base Family:     " << base_family << "\n";
       std::cout << "  - Base Model:      " << base_model << " (0x" << std::hex << base_model << std::dec << ")\n";
       std::cout << "  - Extended Model:  " << extended_model << " (0x" << std::hex << extended_model << std::dec << ")\n";
       std::cout << "  - Extended Family: " << extended_family << "\n";
       std::cout << "  - Stepping ID:     " << stepping << "\n\n";

       std::cout << "Calculated Values for Signature Matching:\n";
       std::cout << "  -> Display Family: " << display_family << "\n";
       std::cout << "  -> Display Model:  " << display_model << " (0x" << std::hex << display_model << std::dec << ")\n\n";

    // Explicit validation check against Skylake-X Server attributes
       std::cout << "Microarchitecture Identification Target:\n";
       if (display_family == 6 && display_model == 85) { // 85 == 0x55
            std::cout << "  [MATCH] Detected: Intel Skylake-X Server / Skylake-SP Microarchitecture.\n";
       } else {
           std::cout << "  [MISMATCH] Host processor is not a Skylake-X Server chip.\n";
       }
       std::cout << "====================================================\n";
#endif 
    this->m_cpuid_basic[0] = base_family;
    this->m_cpuid_basic[1] = base_model;
    this->m_cpuid_basic[2] = stepping;
    this->m_cpuid_basic[3] = extended_model;
    this->m_cpuid_basic[4] = extended_family;
    this->m_cpuid_basic[5] = display_family;
    this->m_cpuid_basic[6] = display_model;
    return (ret_status);
    }
    else 
    {
       ret_status = -1;
       return (ret_status);
    }
}

std::int32_t 
gms::hw_profiling
::SKX_cpu_cpi_retired_t
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
    std::int32_t ret_status = -9999;
    if(__builtin_expect(this->m_totsamples>=1ull,0) && 
       __builtin_expect(this->m_totsamples<=4ull,0))
    {
        for(i = 0ull; i<this->m_totsamples; ++i) 
        {
               a0 = this->m_a.operator[](i);
               b0 = this->m_b.operator[](i);
               this->m_c.operator[](i) = skx_cycles_per_instr(a0,b0);
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
                this->m_c.operator[](i) = skx_cycles_per_instr(a0,b0);
                a1 = this->m_a.operator[](idx1);
                b1 = this->m_b.operator[](idx1);
                this->m_c.operator[](idx1) = skx_cycles_per_instr(a1,b1);
                a2 = this->m_a.operator[](idx2);
                b2 = this->m_b.operator[](idx2);
                this->m_c.operator[](idx2) = skx_cycles_per_instr(a2,b2);
                a3 = this->m_a.operator[](idx3);
                b3 = this->m_b.operator[](idx3);
                this->m_c.operator[](idx3) = skx_cycles_per_instr(a3,b3);
        }

        for(j = i; j<this->m_totsamples; ++j)  
        {
                a0 = this->m_a.operator[](j);
                b0 = this->m_b.operator[](j);
                this->m_c.operator[](j) = skx_cycles_per_instr(a0,b0);
        }
        ret_status = 2;
        return (ret_status);
    }
    return (ret_status);
}

