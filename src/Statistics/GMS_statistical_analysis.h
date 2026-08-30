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
#ifndef __GMS_STATISTICAL_ANALYSIS_H__
#define __GMS_STATISTICAL_ANALYSIS_H__

#include <stddef.h>
            

void run_stat_analysis( const float   * __restrict__  , //ftab
                         const int32_t * __restrict__ , // abscissas
                         float         * __restrict__ ,  // denset
                         float         * __restrict__ , // denset
                         const float     , // denset
                         const float     , // denset
                         const float     ,  // denset
                         const float     , // cubint
                         const float     , // cubint
                         const int32_t   , // number of samples must be a power of 2.
                         const int32_t   , // number of FFT sample points;
                         const char     * __restrict__ );















#endif /*__GMS_STATISTICAL_ANALYSIS_H__*/
