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

#ifndef __GMS_PLATFORM_ERROR_KERNELS_CUH__
#define __GMS_PLATFORM_ERROR_KERNELS_CUH__

#include <cstdint>

extern "C"
void platform_orient_err_cuda(const float * __restrict__ ,
                              const float * __restrict__ ,
                              const float * __restrict__ , 
                              const float * __restrict__ ,
                              const float * __restrict__ ,
                              float * __restrict__ ,//ang,min, azimuth measurement error
                              float * __restrict__ , //ang,min, elevation measurement error
                              const uint32_t,
                              const uint32_t,
                              const uint32_t);


extern "C"
void platform_pos_err_cuda(const float * __restrict__,
                           const float * __restrict__,
                           const float * __restrict__,
                           const float * __restrict__,
                           const float * __restrict__,
                           const float * __restrict__,
                           float * __restrict__ , //min, azimuth error
                           float * __restrict__ ,  //min, elevation error
                           float * __restrict__,    //m,   range error
                           const uint32_t,
                           const uint32_t,
                           const uint32_t);


extern "C"
void elev_mpath_err_cuda(   const float * __restrict__,
                            const float * __restrict__,
                            const float * __restrict__,
                            const float,
                            const float,
                            const float,
                            const float * __restrict__,
                            const float * __restrict__,
                            const float * __restrict__,
                            float * __restrict__, // deg, elevation multipath error
                            const uint32_t,
                            const uint32_t,
                            const uint32_t);


extern "C"
void elev_refract_err_cuda(     const float * __restrict__,
                                const float * __restrict__,
                                const float * __restrict__,
                                const float,
                                const float,
                                const float,
                                const float * __restrict__,
                                const float * __restrict__,
                                const float * __restrict__, 
                                const float,
                                const int32_t,
                                float * __restrict__, //deg, rms error of elevation measurement
                                const uint32_t,
                                const uint32_t,
                                const uint32_t);














#endif /*__GMS_PLATFORM_ERROR_KERNELS_CUH__*/
