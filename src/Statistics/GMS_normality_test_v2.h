
/*
    Adapted from: 
    [freely-available stand-alone C++
    implementation](https://www.johndcook.com/blog/cpp_phi_inverse/).
*/

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

#ifndef __GMS_NORMALITY_TEST_V2_H__
#define __GMS_NORMALITY_TEST_V2_H__ 141020251201

namespace file_info 
{

     static const unsigned int GMS_NORMALITY_TEST_V2_MAJOR = 1;
     static const unsigned int GMS_NORMALITY_TEST_V2_MINOR = 1;
     static const unsigned int GMS_NORMALITY_TEST_V2_MICRO = 0;
     static const unsigned int GMS_NORMALITY_TEST_V2_FULLVER =
       1000U*GMS_NORMALITY_TEST_V2_MAJOR+100U*GMS_NORMALITY_TEST_V2_MINOR+
       10U*GMS_NORMALITY_TEST_V2_MICRO;
     static const char GMS_NORMALITY_TEST_V2_CREATION_DATE[] =   "14-10-2025 12:02PM +00200 (TUE 14 OCT 2025 GMT+2)";
     static const char GMS_NORMALITY_TEST_V2_BUILD_DATE[]    = __DATE__; 
     static const char GMS_NORMALITY_TEST_V2_BUILD_TIME[]    = __TIME__;
     static const char GMS_NORMALITY_TEST_V2_SYNOPSIS[]      = "Adapted from: https://www.johndcook.com/blog/cpp_phi_inverse/";

}

#include <cstdint>
#include "GMS_config.h"
#include "GMS_dyn_array.h"

namespace gms 
{

namespace math 
{

__ATTR_HOT__
__ATTR_OPTIMIZE_03__
double 
calculate_W_statistics(const darray_r8_t &);

}

}

#endif /*__GMS_NORMALITY_TEST_V2_H__*/
