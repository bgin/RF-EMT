
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

/*******************************************************
              CPUID tags types, denotes
              different CPU architectures
*******************************************************/

class cpu_arch_t_tag {};

class cpu_arch_0F_03_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_0F_02_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_2A_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_2D_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_3A_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_3E_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_3C_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_45_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_46_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_25_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_2C_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_1A_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_1E_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_1F_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_2E_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_2F_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_17_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_1D_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_0F_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_0E_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_06_0D_t_tag final : public cpu_arch_t_tag {};

/* Represnts following versions of CPUs
  06_2A/2D/3A/3E/3C/45/46
*/
class cpu_arch_common_models_1_t_tag final : public cpu_arch_t_tag {};

/* Represents following versions of CPUs
  06_25/2C/1A/1E/1F/2E/2F
*/
class cpu_arch_common_models_2_t_tag final : public cpu_arch_t_tag {};

/* Combined models: 06_3C/45/46*/
class cpu_arch_06_3C_45_46_t_tag final : public cpu_arch_t_tag {};

/* Combined models: 06_2A/2D/3A/3E */
class cpu_arch_06_2A_2D_3A_3E_t_tag final : public cpu_arch_t_tag {};

/* Combined models: 06_25/2C/1A/1E/1F/2E/2F */
class cpu_arch_06_25_2C_1A_1E_1F_2E_2F_t_tag final : public cpu_arch_t_tag {};

/***********************************
  Condensed to 3 architecture types
*************************************/

class cpu_arch_sandy_bridge_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_ivy_bridge_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_haswell_t_tag final : public cpu_arch_t_tag {};

/*
   These tags are not implemented
   because of lack of appropriate data:
   1) Latency.
   2) Reciprocal Throughput.
   3) Unfused Micro-ops decomposition.
   4) Fused Micro-ops decomposition.
*/
class cpu_arch_broadwell_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_sky_lake_t_tag final : public cpu_arch_t_tag {};

class cpu_arch_mic_t_tag final : public cpu_arch_t_tag {};

/* Combined tag for 3 architectures
   Not used.*/
class cpu_arch_sb_ib_haswell_t_tag final : public cpu_arch_t_tag {};

/***********************************
   Primitive type data size tags
************************************/

class cpu_arch_data_size_t_tag {};

class cpu_arch_float32_t_tag final : public cpu_arch_data_size_t_tag {};

class cpu_arch_float64_t_tag final : public cpu_arch_data_size_t_tag {};

class cpu_arch_int8_t_tag final : public cpu_arch_data_size_t_tag {};

class cpu_arch_int16_t_tag final : public cpu_arch_data_size_t_tag {};

class cpu_arch_int32_t_tag final : public cpu_arch_data_size_t_tag {};

class cpu_arch_int64_t_tag final : public cpu_arch_data_size_t_tag {};

}

}
#endif /*__GMS_TYPE_TRAITS_HPP__ */