

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

#ifndef __GMS_BOOST_MULTIARRAY_TYPES_H__
#define __GMS_BOOST_MULTIARRAY_TYPES_H__ 

#include <boost/multi_array.hpp>
#include <boost/align/aligned_allocator.hpp>
#include <complex>
#include "GMS_half.h"

namespace gms 
{

namespace common 
{

//Create a 2D array of floats with forced 64-byte memory block alignment
//typedef boost::multi_array<float, 2, boost::alignment::aligned_allocator<float, 64>> AlignedArray2D;

typedef boost::multi_array<half_float::half,1, boost::alignment::aligned_allocator<half_float::half,16>> Aligned16Vec1DF16;

typedef boost::multi_array<half_float::half,1, boost::alignment::aligned_allocator<half_float::half,32>> Aligned32Vec1DF16;

typedef boost::multi_array<half_float::half,1, boost::alignment::aligned_allocator<half_float::half,64>> Aligned64Vec1DF16;

typedef boost::multi_array<float,1, boost::alignment::aligned_allocator<float,16>> Aligned16Vec1DF32;

typedef boost::multi_array<float,1, boost::alignment::aligned_allocator<float,32>> Aligned32Vec1DF32;

typedef boost::multi_array<float,1, boost::alignment::aligned_allocator<float,64>> Aligned64Vec1DF32;

typedef boost::multi_array<double,1, boost::alignment::aligned_allocator<float,16>> Aligned16Vec1DF64;

typedef boost::multi_array<double,1, boost::alignment::aligned_allocator<float,32>> Aligned32Vec1DF64;

typedef boost::multi_array<double,1, boost::alignment::aligned_allocator<float,64>> Aligned64Vec1DF64;

typedef boost::multi_array<std::complex<half_float::half>,1, boost::alignment::aligned_allocator<std::complex<half_float::half>,16>> Aligned16CVec1DF16;

typedef boost::multi_array<std::complex<half_float::half>,1, boost::alignment::aligned_allocator<std::complex<half_float::half>,32>> Aligned32CVec1DF16;

typedef boost::multi_array<std::complex<half_float::half>,1, boost::alignment::aligned_allocator<std::complex<half_float::half>,64>> Aligned64CVec1DF16;

typedef boost::multi_array<std::complex<float>,1, boost::alignment::aligned_allocator<std::complex<float>,16>> Aligned16CVec1DF32;

typedef boost::multi_array<std::complex<float>,1, boost::alignment::aligned_allocator<std::complex<float>,32>> Aligned32CVec1DF32;

typedef boost::multi_array<std::complex<float>,1, boost::alignment::aligned_allocator<std::complex<float>,64>> Aligned64CVec1DF32;

typedef boost::multi_array<std::complex<double>,1, boost::alignment::aligned_allocator<std::complex<double>,16>> Aligned16CVec1DF64;

typedef boost::multi_array<std::complex<double>,1, boost::alignment::aligned_allocator<std::complex<double>,32>> Aligned32CVec1DF64;

typedef boost::multi_array<std::complex<double>,1, boost::alignment::aligned_allocator<std::complex<double>,64>> Aligned64CVec1DF64;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef boost::multi_array<half_float::half,2, boost::alignment::aligned_allocator<half_float::half,16>> Aligned16Mat2DF16;

typedef boost::multi_array<half_float::half,2, boost::alignment::aligned_allocator<half_float::half,32>> Aligned32Mat21DF16;

typedef boost::multi_array<half_float::half,2, boost::alignment::aligned_allocator<half_float::half,64>> Aligned64Mat21DF16;

typedef boost::multi_array<float,2, boost::alignment::aligned_allocator<float,16>> Aligned16Mat2DDF32;

typedef boost::multi_array<float,2, boost::alignment::aligned_allocator<float,32>> Aligned32Mat2DF32;

typedef boost::multi_array<float,2, boost::alignment::aligned_allocator<float,64>> Aligned64Mat2DF32;

typedef boost::multi_array<double,2, boost::alignment::aligned_allocator<float,16>> Aligned16Mat2DF64;

typedef boost::multi_array<double,2, boost::alignment::aligned_allocator<float,32>> Aligned32Mat2DF64;

typedef boost::multi_array<double,2, boost::alignment::aligned_allocator<float,64>> Aligned64Mat2DF64;

typedef boost::multi_array<std::complex<half_float::half>,2, boost::alignment::aligned_allocator<std::complex<half_float::half>,16>> Aligned16CMat2DF16;

typedef boost::multi_array<std::complex<half_float::half>,2, boost::alignment::aligned_allocator<std::complex<half_float::half>,32>> Aligned32CMat2DF16;

typedef boost::multi_array<std::complex<half_float::half>,2, boost::alignment::aligned_allocator<std::complex<half_float::half>,64>> Aligned64CMat2DF16;

typedef boost::multi_array<std::complex<float>,2, boost::alignment::aligned_allocator<std::complex<float>,16>> Aligned16CMat2DF32;

typedef boost::multi_array<std::complex<float>,2, boost::alignment::aligned_allocator<std::complex<float>,32>> Aligned32CMat2DF32;

typedef boost::multi_array<std::complex<float>,2, boost::alignment::aligned_allocator<std::complex<float>,64>> Aligned64CMat2DF32;

typedef boost::multi_array<std::complex<double>,2, boost::alignment::aligned_allocator<std::complex<double>,16>> Aligned16CMat2DF64;

typedef boost::multi_array<std::complex<double>,2, boost::alignment::aligned_allocator<std::complex<double>,32>> Aligned32CMat2DF64;

typedef boost::multi_array<std::complex<double>,2, boost::alignment::aligned_allocator<std::complex<double>,64>> Aligned64CMat2DF64;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}

}

#endif /*__GMS_BOOST_MULTIARRAY_TYPES_H__*/