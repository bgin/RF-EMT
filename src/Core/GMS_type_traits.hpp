
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

#ifndef __GMS_TYPE_TRAITS_HPP__
#define __GMS_TYPE_TRAITS_HPP__ 311020261246

#include <boost/multi_array.hpp>
#include <vector>
#include <valarray>
#include <type_traits>

namespace gms
{

namespace common
{

// mulit_array alignment.
/*
#include <boost/multi_array.hpp>
#include <boost/align/aligned_allocator.hpp>

// Create a 2D array of floats with forced 64-byte memory block alignment
typedef boost::multi_array<float, 2, boost::alignment::aligned_allocator<float, 64>> AlignedArray2D;

AlignedArray2D A(boost::extents[10][10]); 
// A.origin() is now guaranteed to start on a 64-byte boundary

*/

template <typename T, class Container_t>
struct Is_Container_std_vector : public std::integral_constant<bool, false> {};

template <typename T>
struct Is_Container_std_vector<T, std::vector<T>>
    : public std::integral_constant<bool, true> {};

template <typename T, class Container_t>
struct Is_Container_std_valarray 
    : public std::integral_constant<bool, false> {};

template <typename T>
struct Is_Container_std_valarray<T, std::valarray<T>>
    : public std::integral_constant<bool, true> {};

template <typename T, class Container_t>
struct Is_Container_boost_multi_array1d
    : public std::integral_constant<bool, false> {};

template <typename T>
struct Is_Container_boost_multi_array1d<T, boost::multi_array<T, 1>>
    : public std::integral_constant<bool, true> {};

class primitive_type_tag {};

class float32_t_tag final : public primitive_type_tag {};

class float64_t_tag final : public primitive_type_tag {};

/**********************************************
                 Array types tags
                                 Dynamic arrays
***********************************************/

class dyn_container_type_tag {};

class dyn_container_boost_t_tag final : public dyn_container_type_tag {};

class dyn_container_valarray_t_tag final : public dyn_container_type_tag {};

class dyn_container_vector_t_tag final : public dyn_container_type_tag {};

class dyn_container_blaze_t_tag final : public dyn_container_type_tag {};

class dyn_container_armadillo_t_tag final : public dyn_container_type_tag {};

class dyn_container_raw_ptr_t_tag final : public dyn_container_type_tag {};

class dyn_container_smart_ptr_t_tag final : public dyn_container_type_tag {};

class dyn_container_own_impl_t_tag final : public dyn_container_type_tag {};

/*****************************************************
                  Array types tags
                                  Static arrays.
******************************************************/

class stat_container_type_tag {};

class stat_container_std_array_t_tag final : public stat_container_type_tag {};

class stat_container_builtin_array_t_tag final : public stat_container_type_tag {};
    
class stat_container_rf_emt_half_impl_t_tag final : public stat_container_type_tag {};

class stat_container_rf_emt_float_impl_t_tag final : public stat_container_type_tag {};

class stat_container_rf_emt_double_impl_t_tag final : public stat_container_type_tag {};

class stat_container_rf_emt_complex_half_impl_t_tag final : public stat_container_type_tag {};

class stat_container_rf_emt_complex_float_impl_t_tag final : public stat_container_type_tag {};

class stat_container_rf_emt_complex_double_impl_t_tag final : public stat_container_type_tag {};

}

}
#endif /*__GMS_TYPE_TRAITS_HPP__ */