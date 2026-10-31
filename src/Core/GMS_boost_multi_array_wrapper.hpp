
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

#ifndef __GMS_BOOST_MULTI_ARRAY_WRAPPER_HPP__
#define __GMS_BOOST_MULTI_ARRAY_WRAPPER_HPP__ 311020261325

#include <boost/multi_array.hpp>
#include <boost/align/aligned_allocator.hpp>
#include <complex>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <algorithm>
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_BOOST_MULTI_ARRAY_WRAPPER_MAJOR = 1;
     static const unsigned int GMS_BOOST_MULTI_ARRAY_WRAPPER_MINOR = 1;
     static const unsigned int GMS_BOOST_MULTI_ARRAY_WRAPPER_MICRO = 0;
     static const unsigned int GMS_BOOST_MULTI_ARRAY_WRAPPER_FULLVER =
       1000U*GMS_BOOST_MULTI_ARRAY_WRAPPER_MAJOR+100U*GMS_BOOST_MULTI_ARRAY_WRAPPER_MINOR+
       10U*GMS_BOOST_MULTI_ARRAY_WRAPPER_MICRO;
     static const char GMS_BOOST_MULTI_ARRAY_WRAPPER_CREATION_DATE[]  = "31-10-2026 13:26PM +00200 (SAT 31 OCT 2026 GMT+2)";
     static const char GMS_BOOST_MULTI_ARRAY_WRAPPER_BUILD_DATE[]    = __DATE__; 
     static const char GMS_BOOST_MULTI_ARRAY_WRAPPER_BUILD_TIME[]    = __TIME__;
     static const char GMS_BOOST_MULTI_ARRAY_WRAPPER_AUTHOR[]        = "Programmer: Bernard Gingold, beniekg@gmail.com";            
     static const char GMS_BOOST_MULTI_ARRAY_WRAPPER_SYNOPSIS[]      = "Boost-MultiArray simple wrapper.";
}

namespace gms
{

namespace common
{

// Traits to validate supported types and alignments
template <typename T>
struct is_supported_type : std::false_type {};

template <> struct is_supported_type<float> : std::true_type {};
template <> struct is_supported_type<double> : std::true_type {};
template <> struct is_supported_type<std::complex<float>> : std::true_type {};
template <> struct is_supported_type<std::complex<double>> : std::true_type {};

template <std::size_t Align>
struct is_supported_alignment : std::false_type {};

template <> struct is_supported_alignment<16> : std::true_type {};
template <> struct is_supported_alignment<32> : std::true_type {};
template <> struct is_supported_alignment<64> : std::true_type {};

template <std::size_t Dims>
struct is_supported_dim : std::false_type {};

template <> struct is_supported_dim<1> : std::true_type {};
template <> struct is_supported_dim<2> : std::true_type {};
template <> struct is_supported_dim<3> : std::true_type {};

template <typename T, std::size_t Dims, std::size_t Alignment>
class AlignedMultiArray {
    static_assert(is_supported_type<T>::value, "Unsupported data type for AlignedMultiArray.");
    static_assert(is_supported_dim<Dims>::value, "Unsupported dimension. Only 1D, 2D, and 3D are supported.");
    static_assert(is_supported_alignment<Alignment>::value, "Unsupported alignment. Only 16, 32, and 64 bytes are supported.");

public:

    using value_type = T;
    using allocator_type = boost::alignment::aligned_allocator<T, Alignment>;
    using internal_array_type = boost::multi_array<T, Dims, allocator_type>;
    using size_type = typename internal_array_type::size_type;
    using index = typename internal_array_type::index;
    using reference = typename internal_array_type::reference;
    using const_reference = typename internal_array_type::const_reference;
    using element_iterator = typename internal_array_type::element_iterator;
    using const_element_iterator = typename internal_array_type::const_element_iterator;
    
    AlignedMultiArray() = default;

    template <typename... Extents, 
              typename = std::enable_if_t<sizeof...(Extents) == Dims>>
    explicit AlignedMultiArray(Extents... extents) {
        m_array.resize(boost::extents[extents...]);
    }

    explicit AlignedMultiArray(const std::array<size_type, Dims>& extents) {
        m_array.resize(extents);
    }

    AlignedMultiArray(const AlignedMultiArray&) = default;
    AlignedMultiArray& operator=(const AlignedMultiArray&) = default;
    AlignedMultiArray(AlignedMultiArray&&) noexcept = default;
    AlignedMultiArray& operator=(AlignedMultiArray&&) noexcept = default;

    template <typename... Extents, 
              typename = std::enable_if_t<sizeof...(Extents) == Dims>>
    void resize(Extents... extents) {
        m_array.resize(boost::extents[extents...]);
    }

    void resize(const std::array<size_type, Dims>& extents) {
        m_array.resize(extents);
    }

    // Note: boost::multi_array operator[] returns a reference to a sub-array (or element if 1D)
    auto operator[](index idx) -> decltype(std::declval<internal_array_type&>()[idx]) {
        return m_array[idx];
    }

    auto operator[](index idx) const -> decltype(std::declval<const internal_array_type&>()[idx]) {
        return m_array[idx];
    }

    // Parentheses operator for direct element access across all dimensions
    template <typename... Indices, typename = std::enable_if_t<sizeof...(Indices) == Dims>>
    reference operator()(Indices... indices) {
        return m_array(indices...);
    }

    template <typename... Indices, typename = std::enable_if_t<sizeof...(Indices) == Dims>>
    const_reference operator()(Indices... indices) const {
        return m_array(indices...);
    }

    T* data() noexcept { return m_array.data(); }
    const T* data() const noexcept { return m_array.data(); }

    T* origin() noexcept { return m_array.origin(); }
    const T* origin() const noexcept { return m_array.origin(); }

    element_iterator num_elements_begin() { return m_array.num_elements_begin(); }
    element_iterator num_elements_end() { return m_array.num_elements_end(); }

    const_element_iterator num_elements_begin() const { return m_array.num_elements_begin(); }
    const_element_iterator num_elements_end() const { return m_array.num_elements_end(); }

    size_type num_elements() const noexcept { return m_array.num_elements(); }
    size_type num_dimensions() const noexcept { return Dims; }
    
    const size_type* shape() const noexcept { return m_array.shape(); }
    const index* strides() const noexcept { return m_array.strides(); }
    const index* index_bases() const noexcept { return m_array.index_bases(); }

    void fill(const T& value)
     {
        std::fill_n(m_array.origin(), m_array.num_elements(), value);
    }

    internal_array_type& get_internal() noexcept { return m_array; }
    const internal_array_type& get_internal() const noexcept { return m_array; }

private:
    internal_array_type m_array;
};

// Aliases for 1D, 2D, 3D and specific alignments for ease of use
template <typename T, std::size_t Alignment = 32>
using AlignedArray1D = AlignedMultiArray<T, 1, Alignment>;

template <typename T, std::size_t Alignment = 32>
using AlignedArray2D = AlignedMultiArray<T, 2, Alignment>;

template <typename T, std::size_t Alignment = 32>
using AlignedArray3D = AlignedMultiArray<T, 3, Alignment>;

} // common

} // gms

#endif // __GMS_BOOST_MULTI_ARRAY_WRAPPER_HPP__
