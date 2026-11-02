
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

#ifndef __GMS_OFDM_SIGNAL_MODULATOR_H__
#define __GMS_OFDM_SIGNAL_MODULATOR_H__ 021120260635

#include <algorithm>
#include <complex>
#include <cstdint>
#include <immintrin.h>
#include <random>
#include "GMS_config.h"
#include "GMS_dyn_array.h"
#include "GMS_sse_memset.h"

namespace file_info 
{

     const unsigned int GMS_OFDM_SIGNAL_MODULATOR_MAJOR = 2;
     const unsigned int GMS_OFDM_SIGNAL_MODULATOR_MINOR = 0;
     const unsigned int GMS_OFDM_SIGNAL_MODULATOR_MICRO = 0;
     const unsigned int GMS_OFDM_SIGNAL_MODULATOR_FULLVER =
       1000U*GMS_OFDM_SIGNAL_MODULATOR_MAJOR+100U*GMS_OFDM_SIGNAL_MODULATOR_MINOR+
       10U*GMS_OFDM_SIGNAL_MODULATOR_MICRO;
     const char * const GMS_OFDM_SIGNAL_MODULATOR_CREATION_DATE = "02-11-2026 06:35 AM +00200 (MON 02 NOV 2026 GMT+2)";
     const char * const GMS_OFDM_SIGNAL_MODULATOR_BUILD_DATE    = __DATE__ " " __TIME__;
     static const char  GMS_OFDM_SIGNAL_MODULATOR_AUTHOR[]      = "Programmer: Bernard Gingold, beniekg@gmail.com"; 
     const char * const GMS_OFDM_SIGNAL_MODULATOR_SYNOPSIS      = "Basic OFDM signal modulator classes (float and double).";

}

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(OFDM_SIGNAL_MODULATOR_USE_CEPHES)
#define OFDM_SIGNAL_MODULATOR_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (OFDM_SIGNAL_MODULATOR_USE_PMC_INSTRUMENTATION)
#define OFDM_SIGNAL_MODULATOR_USE_PMC_INSTRUMENTATION 0
#endif 

#if (OFDM_SIGNAL_MODULATOR_USE_PMC_INSTRUMENTATION) == 1
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

#if (OFDM_SIGNAL_MODULATOR_USE_CEPHES) == 0
#include <cmath>
#else 
#include "GMS_cephes_sin_cos.h"
#endif 

namespace gms 
{

namespace radiolocation
{


enum class PrefetchMode : std::int32_t { Disabled, Enabled };


enum class QamScheme : std::int32_t 
{
  QAM16 = 16,
  QAM32 = 32,
  QAM64 = 64,
  QAM128 = 128,
  QAM256 = 256
};

template<PrefetchMode PMode = PrefetchMode::Disabled>
struct alignas(64) ofdm_signal_modulator_f32_t
{
    // Total active users/subcarriers (Up to 1024)
    // FFT Size (typically 1024, 2048, or 4096)
    // CP length in samples
    std::size_t m_num_channels;
    std::size_t m_fft_size;
    std::size_t m_cyc_prefix_len;

    darray_r4_t m_subcarriers_re;
    darray_r4_t m_subcarriers_im;
    darray_r4_t m_time_signal_re;
    darray_r4_t m_time_signal_im;

    __ATTR_ALWAYS_INLINE__
    inline 
    ofdm_signal_modulator_f32_t() noexcept(true)
    :
    m_num_channels{},
    m_fft_size{},
    m_cyc_prefix_len{},
    m_subcarriers_re(),
    m_subcarriers_im(),
    m_time_signal_re(),
    m_time_signal_im()
    {

    }
    

    __ATTR_ALWAYS_INLINE__
    inline
    ofdm_signal_modulator_f32_t(const std::size_t num_channels,
                                const std::size_t fft_size,
                                const std::size_t cyc_prefix_len) noexcept(false)
    :
    m_num_channels((num_channels>1024ull)?1024ull:num_channels),
    m_fft_size(fft_size),
    m_cyc_prefix_len(cyc_prefix_len),
    m_subcarriers_re(m_fft_size),
    m_subcarriers_im(m_fft_size),
    m_time_signal_re(m_fft_size+m_cyc_prefix_len),
    m_time_signal_im(m_fft_size+m_cyc_prefix_len)
    {

    }

    __ATTR_ALWAYS_INLINE__
    inline
    ofdm_signal_modulator_f32_t(const ofdm_signal_modulator_f32_t &other) noexcept(false) 
    :
    m_num_channels(other.m_num_channels),
    m_fft_size(other.m_fft_size),
    m_cyc_prefix_len(other.m_cyc_prefix_len),
    m_subcarriers_re(other.m_subcarriers_re),
    m_subcarriers_im(other.m_subcarriers_im),
    m_time_signal_re(other.m_time_signal_re),
    m_time_signal_im(other.m_time_signal_im)
    {

    }

    __ATTR_ALWAYS_INLINE__
    inline
    ofdm_signal_modulator_f32_t(ofdm_signal_modulator_f32_t &&other) noexcept(true) 
    :
    m_num_channels(other.m_num_channels),
    m_fft_size(other.m_fft_size),
    m_cyc_prefix_len(other.m_cyc_prefix_len),
    m_subcarriers_re(std::move(other.m_subcarriers_re)),
    m_subcarriers_im(std::move(other.m_subcarriers_im)),
    m_time_signal_re(std::move(other.m_time_signal_re)),
    m_time_signal_im(std::move(other.m_time_signal_im))
    {

    }

     __ATTR_ALWAYS_INLINE__
    inline
    ~ofdm_signal_modulator_f32_t() noexcept(true) 
    {

    }
    
    __ATTR_ALWAYS_INLINE__
    inline
    ofdm_signal_modulator_f32_t & 
        operator=(const ofdm_signal_modulator_f32_t &other) noexcept(false) 
    {
        if(__builtin_expect(this==&other,0)) { return (*this);}
        this->m_num_channels   = other.m_num_channels;
        this->m_fft_size       = other.m_fft_size;
        this->m_cyc_prefix_len = other.m_cyc_prefix_len;
        this->m_subcarriers_re.operator=(other.m_subcarriers_re);
        this->m_subcarriers_im.operator=(other.m_subcarriers_im);
        this->m_time_signal_re.operator=(other.m_time_signal_re);
        this->m_time_signal_im.operator=(other.m_time_signal_im);
        return (*this);
    }

    __ATTR_ALWAYS_INLINE__
    inline
    ofdm_signal_modulator_f32_t & 
        operator=(ofdm_signal_modulator_f32_t &&other) noexcept(true) 
    {
        if(__builtin_expect(this==&other,0)) { return (*this);}
        this->m_num_channels   = other.m_num_channels;
        this->m_fft_size       = other.m_fft_size;
        this->m_cyc_prefix_len = other.m_cyc_prefix_len;
        this->m_subcarriers_re.operator=(std::move(other.m_subcarriers_re));
        this->m_subcarriers_im.operator=(std::move(other.m_subcarriers_im));
        this->m_time_signal_re.operator=(std::move(other.m_time_signal_re));
        this->m_time_signal_im.operator=(std::move(other.m_time_signal_im));
        return (*this);
    }

    __ATTR_ALWAYS_INLINE__
    inline
    std::size_t get_total_output_size() const 
    {
        return (this->m_fft_size+this->m_cyc_prefix_len);
    }

    __ATTR_ALWAYS_INLINE__
    inline
    bool clear_subcarriers_re(const float fill) noexcept(true) 
    {
        using namespace gms::common;
        bool ret_status = false;
        sse_memset_unroll8x_ps(this->m_subcarriers_re.m_data,filler,this->m_fft_size);
        ret_status = true;
        return (ret_status);
    }

    __ATTR_ALWAYS_INLINE__
    inline
    bool clear_subcarriers_im(const float fill) noexcept(true) 
    {
        using namespace gms::common;
        bool ret_status = false;
        sse_memset_unroll8x_ps(this->m_subcarriers_im.m_data,filler,this->m_fft_size);
        ret_status = true;
        return (ret_status);
    }

    __ATTR_ALWAYS_INLINE__
    inline
    bool clear_time_signal_re(const float fill) noexcept(true) 
    {
        using namespace gms::common;
        bool ret_status = false;
        sse_memset_unroll8x_ps(this->m_time_signal_re.m_data,fill,get_total_output_size());
        ret_status = true;
        return (ret_status);
    }

    __ATTR_ALWAYS_INLINE__
    inline
    bool clear_time_signal_im(const float fill) noexcept(true) 
    {
        using namespace gms::common;
        bool ret_status = false;
        sse_memset_unroll8x_ps(this->m_time_signal_im.m_data,fill,get_total_output_size());
        ret_status = true;
        return (ret_status);
    }

    /*
       * Generates random binary streams, maps to specified QAM, and loads
       * into subcarriers.
    */
    __ATTR_ALWAYS_INLINE__
    inline
    bool 
    generate_channels_randomly(QamScheme scheme) noexcept(true)
    {
        bool ret_status = false;
        std::int32_t side_levels = std::sqrt(static_cast<int>(scheme));
#if (OFDM_SIGNAL_MODULATOR_USE_CEPHES) == 1
        float normalization_factor =
        gms::math::ceph_sqrtf(2.0f * (static_cast<int>(scheme) - 1) / 3.0f);
#else 
        float normalization_factor =
        std::sqrt(2.0f * (static_cast<int>(scheme) - 1) / 3.0f);
#endif 
        thread_local std::uniform_int_distribution<std::int32_t> uni_dist;
        thread_local std::uint64_t seed_val{};
        thread_local std::mt19937_64 rv_uni_gen;
        std::size_t i,j;
        seed_val   = __rdtsc();
        uni_dist   = std::uniform_int_distribution<std::int32_t>(0,side_levels-1);
        rv_uni_gen = std::mt19937_64(seed_val);
        const float inv_normalization_factor = 1.0f/normalization_factor;
        const float f_side_levels            = static_cast<float>(side_levels-1);
        float rv_re;
        float rv_im;
        float rv_qam_re;
        float rv_qam_im;
        for(i = 0ull; i!=ROUND_DOWN(this->m_num_channels,3ull); i += 4ull) 
        {
            const std::size_t idx1 = i+1ull;
            const std::size_t idx2 = i+2ull;
            const std::size_t idx3 = i+3ull;

            rv_re = uni_dist.operator()(rv_uni_gen);
            rv_im = uni_dist.operator()(rv_uni_gen);
            rv_qam_re = (2.0f*rv_re-f_side_levels)*inv_normalization_factor;
            this->m_subcarriers_re.m_data[i] = rv_qam_re;
            rv_qam_im = (2.0f*rv_im-f_side_levels)*inv_normalization_factor;
            this->m_subcarriers_im.m_data[i] = rv_qam_im;

            rv_re = uni_dist.operator()(rv_uni_gen);
            rv_im = uni_dist.operator()(rv_uni_gen);
            rv_qam_re = (2.0f*rv_re-f_side_levels)*inv_normalization_factor;
            this->m_subcarriers_re.m_data[idx1] = rv_qam_re;
            rv_qam_im = (2.0f*rv_im-f_side_levels)*inv_normalization_factor;
            this->m_subcarriers_im.m_data[idx1] = rv_qam_im;

            rv_re = uni_dist.operator()(rv_uni_gen);
            rv_im = uni_dist.operator()(rv_uni_gen);
            rv_qam_re = (2.0f*rv_re-f_side_levels)*inv_normalization_factor;
            this->m_subcarriers_re.m_data[idx2] = rv_qam_re;
            rv_qam_im = (2.0f*rv_im-f_side_levels)*inv_normalization_factor;
            this->m_subcarriers_im.m_data[idx2] = rv_qam_im;

            rv_re = uni_dist.operator()(rv_uni_gen);
            rv_im = uni_dist.operator()(rv_uni_gen);
            rv_qam_re = (2.0f*rv_re-f_side_levels)*inv_normalization_factor;
            this->m_subcarriers_re.m_data[idx3] = rv_qam_re;
            rv_qam_im = (2.0f*rv_im-f_side_levels)*inv_normalization_factor;
            this->m_subcarriers_im.m_data[idx3] = rv_qam_im;
        } 

        for(j = i; j<this->m_num_channels; ++j)  
        {
            rv_re = uni_dist.operator()(rv_uni_gen);
            rv_im = uni_dist.operator()(rv_uni_gen);
            rv_qam_re = (2.0f*rv_re-f_side_levels)*inv_normalization_factor;
            this->m_subcarriers_re.m_data[j] = rv_qam_re;
            rv_qam_im = (2.0f*rv_im-f_side_levels)*inv_normalization_factor;
            this->m_subcarriers_im.m_data[j] = rv_qam_im;
        }
        ret_status = true;
        return (ret_status);
    }

    __ATTR_ALWAYS_INLINE__
    inline
    bool
    compute_idft_and_apply_cycle_prefix()
    {
        bool ret_status = false;
        const std::size_t total_samples = this->get_total_output_size();
        const float f_fft_size = static_cast<float>(this->m_fft_size);
#if (OFDM_SIGNAL_MODULATOR_USE_CEPHES) == 1
        const float inv_fft = 1.0f / gms::math::ceph_sqrtf(f_fft_size);
#else 
        const float inv_fft = 1.0f / std::sqrt(f_fft_size);
#endif 
        __m512 inv_fft_vec = _mm512_set1_ps(inv_fft);
        float twiddle_base = 2.0f * M_PI / f_fft_size;
       for (std::size_t n = 0ull; n < this->m_fft_size; ++n) 
       {
           __m512 sum_r = _mm512_setzero_ps();
           __m512 sum_i = _mm512_setzero_ps();
           std::size_t k = 0;
           std::size_t vector_limit = (this->m_num_channels & ~15);

          for (; k < vector_limit; k += 16ull) 
          {
              if constexpr (PMode == PrefetchMode::Enabled) 
              {
                 _mm_prefetch(
                 reinterpret_cast<const char *>(&this->m_subcarriers_re.m_data[k + 32ull]),
                 _MM_HINT_T0);
                 _mm_prefetch(
                 reinterpret_cast<const char *>(&this->m_subcarriers_im.m_data[k + 32ull]),
                 _MM_HINT_T0);
              }
              __m512 sub_r = _mm512_load_ps(&this->m_subcarriers_re.m_data[k]);
              __m512 sub_i = _mm512_load_ps(&this->m_subcarriers_im.m_data[k]);

              alignas(64) float angle_array[16];
              for (int m = 0; m < 16; ++m) 
              {
                   angle_array[m] = twiddle_base * n * (k + m);
             }
             __m512 alpha = _mm512_load_ps(angle_array);
              alignas(64) float cos_arr[16], sin_arr[16];
              for (int m = 0; m < 16; ++m) 
              {
#if (OFDM_SIGNAL_MODULATOR_USE_CEPHES) == 1
                  cos_arr[m] = gms::math::ceph_cosf(angle_array[m]);
                  sin_arr[m] = gms::math::ceph_sinf(angle_array[m]);
#else 
                  cos_arr[m] = std::cos(angle_array[m]);
                  sin_arr[m] = std::sin(angle_array[m]);
#endif
              }
              __m512 W_r = _mm512_load_ps(cos_arr);
              __m512 W_i = _mm512_load_ps(sin_arr);
        // Real part = sub_r * W_r - sub_i * W_i
              __m512 r_prod = _mm512_mul_ps(sub_r, W_r);
              sum_r = _mm512_fmsubadd_ps(sub_i, W_i,r_prod);
        // Imaginary part = sub_r * W_i + sub_i * W_r
              __m512 i_prod = _mm512_mul_ps(sub_r, W_i);
              sum_i = _mm512_fmadd_ps(sub_i, W_r, i_prod);
          }

          std::size_t rem = this->m_num_channels - k;
          if (rem > 0ull) 
          {
              __mmask16 mask = (1U << rem) - 1ull;
              __m512 sub_r = _mm512_maskz_load_ps(mask, &this->m_subcarriers_re.m_data[k]);
              __m512 sub_i = _mm512_maskz_load_ps(mask, &this->m_subcarriers_im.m_data[k]);
              alignas(64) float angle_array[16] = {0};
              for (uint32_t m = 0; m < rem; ++m) 
              {
                  angle_array[m] = twiddle_base * n * (k + m);
              }
              alignas(64) float cos_arr[16] = {0}, sin_arr[16] = {0};
              for (uint32_t m = 0; m < rem; ++m) 
              {
#if (OFDM_SIGNAL_MODULATOR_USE_CEPHES) == 1
                  cos_arr[m] = gms::math::ceph_cosf(angle_array[m]);
                  sin_arr[m] = gms::math::ceph_sinf(angle_array[m]);
#else 
                  cos_arr[m] = std::cos(angle_array[m]);
                  sin_arr[m] = std::sin(angle_array[m]);
#endif 
              }
               __m512 W_r = _mm512_load_ps(cos_arr);
               __m512 W_i = _mm512_load_ps(sin_arr);
               __m512 r_prod = _mm512_mul_ps(sub_r, W_r);
               sum_r = _mm512_mask_fmsubadd_ps(sum_r, mask, sub_i, r_prod);
               __m512 i_prod = _mm512_mul_ps(sub_r, W_i);
               sum_i = _mm512_mask_fmadd_ps(sum_i, mask, sub_i, i_prod);
           }
            // Horizontal reduction of sums across the 16 registers lanes
               float final_r = _mm512_reduce_add_ps(sum_r) * inv_fft;
               float final_i = _mm512_reduce_add_ps(sum_i) * inv_fft;
     
               this->m_time_signal_re.m_data[this->m_cyc_prefix_len + n] = final_r;
               this->m_time_signal_im.m_data[this->m_cyc_prefix_len + n] = final_i;
        }
        std::size_t cp_idx = 0ull;
        std::size_t cp_vector_limit = (this->m_cyc_prefix_len & ~15);
        float * __restrict__ src_cp_r = &this->m_time_signal_re.m_data[this->m_fft_size];
        float * __restrict__ src_cp_i = &this->m_time_signal_im.m_data[this->m_fft_size];
        for (; cp_idx < cp_vector_limit; cp_idx += 16ull) 
        {
             __m512 val_r = _mm512_load_ps(&src_cp_r[cp_idx]);
             __m512 val_i = _mm512_load_ps(&src_cp_i[cp_idx]);
             _mm512_store_ps(&this->m_subcarriers_re.m_data[cp_idx], val_r);
             _mm512_store_ps(&this->m_time_signal_im.m_data[cp_idx], val_i);
        }
         std::size_t cp_rem = this->m_cyc_prefix_len - cp_idx;
         if (cp_rem > 0) 
         {
             __mmask16 mask = (1U << cp_rem) - 1U;
             __m512 val_r = _mm512_maskz_load_ps(mask, &src_cp_r[cp_idx]);
             __m512 val_i = _mm512_maskz_load_ps(mask, &src_cp_i[cp_idx]);
             _mm512_mask_store_ps(&this->m_time_signal_re.m_data[cp_idx], mask, val_r);
             _mm512_mask_store_ps(&this->m_time_signal_im.m_data[cp_idx], mask, val_i);
        }
        ret_status = true;
        return (ret_status);
    }
    
};
    
    

} // radiolocation 

} // gms

#endif /*__GMS_OFDM_SIGNAL_MODULATOR_H__*/