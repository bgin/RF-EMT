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

#ifndef __GMS_AM_BB_CMPLX_COS_SEQUENCE_H__
#define __GMS_AM_BB_CMPLX_COS_SEQUENCE_H__ 181120250648


namespace file_info 
{

     static const unsigned int GMS_AM_BB_CMPLX_COS_SEQUENCE_MAJOR = 1;
     static const unsigned int GMS_AM_BB_CMPLX_COS_SEQUENCE_MINOR = 1;
     static const unsigned int GMS_AM_BB_CMPLX_COS_SEQUENCE_MICRO = 0;
     static const unsigned int GMS_AM_BB_CMPLX_COS_SEQUENCE_FULLVER =
       1000U*GMS_AM_BB_CMPLX_COS_SEQUENCE_MAJOR+100U*GMS_AM_BB_CMPLX_COS_SEQUENCE_MINOR+
       10U*GMS_AM_BB_CMPLX_COS_SEQUENCE_MICRO;
     static const char GMS_AM_BB_CMPLX_COS_SEQUENCE_CREATION_DATE[] = "18-11-2025 06:48 +00200 (TUE 18 NOV 2025 GMT+2)";
     static const char GMS_AM_BB_CMPLX_COS_SEQUENCE_BUILD_DATE[]    = __DATE__; 
     static const char GMS_AM_BB_CMPLX_COS_SEQUENCE_BUILD_TIME[]    = __TIME__;
     static const char GMS_AM_BB_CMPLX_COS_SEQUENCE_SYNOPSIS[]      = "AM baseband complex cosine signal sequence (vertically summed).";

}

#include <cstdint>
#include <string>
#include <vector>
#include <valarray>
#include "GMS_config.h"
#include "GMS_am_bb_cmplx_cos_signal.h"

#if !defined(AM_BB_CMPLX_COS_SEQUENCE_SOFT_PREFETCH)
#define AM_BB_CMPLX_COS_SEQUENCE_SOFT_PREFETCH 1 
#endif

namespace gms 
{

namespace radiolocation 
{

            struct alignas(64) am_bb_cmplx_cos_signal_sequence_t final
            {
                   std::size_t                           m_nsamples; 
                   std::size_t                           m_nsignals;
                   std::vector<std::uint32_t>            m_nK_values;
                   std::vector<float>                    m_n_values;
                   std::vector<float>                    m_A_values;
                   std::vector<float>                    m_P_values;
                   std::vector<am_bb_cmplx_cos_signal_t> m_bb_cmplx_cos_signals;
                   darray_r4_t                           m_re_sequence;
                   darray_r4_t                           m_im_sequence;

                   am_bb_cmplx_cos_signal_sequence_t() = default;

                   am_bb_cmplx_cos_signal_sequence_t(const std::size_t,
                                                     const std::size_t,
                                                     const std::vector<std::uint32_t> &,
                                                     const std::vector<float> &,
                                                     const std::vector<float> &,
                                                     const std::vector<float> &) noexcept(false);

                   am_bb_cmplx_cos_signal_sequence_t(am_bb_cmplx_cos_signal_sequence_t &&) noexcept(true);

                   am_bb_cmplx_cos_signal_sequence_t(const am_bb_cmplx_cos_signal_sequence_t &) noexcept(false);

                   ~am_bb_cmplx_cos_signal_sequence_t() noexcept(false);

                   am_bb_cmplx_cos_signal_sequence_t & 
                                           operator=(const am_bb_cmplx_cos_signal_sequence_t &) noexcept(false);

                   am_bb_cmplx_cos_signal_sequence_t &
                                           operator=(am_bb_cmplx_cos_signal_sequence_t &&) noexcept(true);

                   static void create_sequence_plot(  const std::uint32_t,
                                                   const float * __restrict,
                                                   const float * __restrict,
                                                   const std::string &,
                                                   const std::string &,
                                                   const bool );
                /*
                    Signal sequence generators -- SSE,AVX,AVX512 vectorized.
                    ***WARNING***
                                   Expected large memory-bound due to
                                   low arithemtic-to-load-store operations
                                   ratio.
                */
                  std::int32_t 
                signal_sequence_sse42_u16x(const float         * __restrict__,
                                           am_bb_cmplx_cos_signal_pdf_params_t &,
                                           am_bb_cmplx_cos_signal_rand_distr,
                                           const float          * __restrict__,
                                           std::int32_t * __restrict__,
                                           const std::uint32_t);                                                   

                std::int32_t 
                signal_sequence_sse42_u10x(const float         * __restrict__,
                                           am_bb_cmplx_cos_signal_pdf_params_t &,
                                           am_bb_cmplx_cos_signal_rand_distr,
                                           const float          * __restrict__,
                                           std::int32_t * __restrict__,
                                           const std::uint32_t);                                                   

                std::int32_t 
                signal_sequence_sse42_u6x(const float         * __restrict__,
                                          am_bb_cmplx_cos_signal_pdf_params_t &,
                                          am_bb_cmplx_cos_signal_rand_distr,
                                          const float          * __restrict__,
                                          std::int32_t * __restrict__,
                                          const std::uint32_t);

                std::int32_t 
                signal_sequence_avx_u16x(const float         * __restrict__,
                                         am_bb_cmplx_cos_signal_pdf_params_t &,
                                         am_bb_cmplx_cos_signal_rand_distr,
                                         const float          * __restrict__,
                                         std::int32_t * __restrict__,
                                         const std::uint32_t);                                                   

                std::int32_t 
                signal_sequence_avx_u10x(const float         * __restrict__,
                                         am_bb_cmplx_cos_signal_pdf_params_t &,
                                         am_bb_cmplx_cos_signal_rand_distr,
                                         const float          * __restrict__,
                                         std::int32_t * __restrict__,
                                         const std::uint32_t);                                                   

                std::int32_t 
                signal_sequence_avx_u6x(const float         * __restrict__,
                                        am_bb_cmplx_cos_signal_pdf_params_t &,
                                        am_bb_cmplx_cos_signal_rand_distr,
                                        const float          * __restrict__,
                                        std::int32_t * __restrict__,
                                        const std::uint32_t);

                std::int32_t 
                signal_sequence_avx512_u16x(const float         * __restrict__,
                                            am_bb_cmplx_cos_signal_pdf_params_t &,
                                            am_bb_cmplx_cos_signal_rand_distr,
                                            const float          * __restrict__,
                                            std::int32_t * __restrict__,
                                            const std::uint32_t);                                                   

                std::int32_t 
                signal_sequence_avx512_u10x(const float         * __restrict__,
                                            am_bb_cmplx_cos_signal_pdf_params_t &,
                                            am_bb_cmplx_cos_signal_rand_distr,
                                            const float          * __restrict__,
                                            std::int32_t * __restrict__,
                                            const std::uint32_t);                                                   

                std::int32_t 
                signal_sequence_avx512_u6x(const float         * __restrict__,
                                           am_bb_cmplx_cos_signal_pdf_params_t &,
                                           am_bb_cmplx_cos_signal_rand_distr,
                                           const float          * __restrict__,
                                           std::int32_t * __restrict__,
                                           const std::uint32_t);
                                             
                   
            };

            auto 
             operator<<(std::ostream &,
                        const am_bb_cmplx_cos_signal_sequence_t &)->std::ostream &;

}

}


















#endif /*__GMS_AM_BB_CMPLX_COS_SEQUENCE_H__*/