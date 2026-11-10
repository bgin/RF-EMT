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

#ifndef __GMS_TERRAIN_TOPOGRAPHY_H__
#define __GMS_TERRAIN_TOPOGRAPHY_H__

#include <cstdint>
#include <iostream>
#include <iomainp>
#include <boost/multi_array.hpp>
#include <boost/align/aligned_allocator.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/systems/si.hpp>
#include "GMS_config.h"


namespace file_info 
{

     static const unsigned int GMS_TERRAIN_TOPOGRAPHY_MAJOR = 1;
     static const unsigned int GMS_TERRAIN_TOPOGRAPHY_MINOR = 1;
     static const unsigned int GMS_TERRAIN_TOPOGRAPHY_MICRO = 0;
     static const unsigned int GMS_TERRAIN_TOPOGRAPHY_FULLVER =
       1000U*GMS_TERRAIN_TOPOGRAPHY_MAJOR+100U*GMS_TERRAIN_TOPOGRAPHY_MINOR+
       10U*GMS_TERRAIN_TOPOGRAPHY_MICRO;
     static const char GMS_TERRAIN_TOPOGRAPHY_CREATION_DATE[] = "10-11-2025 07:03AM +00200 (TUE 10 NOV 2026 GMT+2)";
     static const char GMS_TERRAIN_TOPOGRAPHY_BUILD_DATE[]    = __DATE__; 
     static const char GMS_TERRAIN_TOPOGRAPHY_BUILD_TIME[]    = __TIME__;
     static const char GMS_TERRAIN_TOPOGRAPHY_SYNOPSIS[]      = "Simple Terrain Topography implementation.";

}

namespace gms 
{

namespace radiolocation 
{



*/
/******************************************
    Terrain topography as 2D height map.
	Based on boost::units::quantity<length>
	in meters above sea level.
*******************************************/

template<typename T,
        const std::size_t IMS,
        const std::size_t IME,
        const std::size_t JMS,
        const std::size_t JME>
struct alignas(64) terrain_topography_t 
{
    static_assert(IME>IMS,"[IME>IMS] -- ***FAILED***");
    static_assert(JME>JMS,"[JME>JMS] -- ***FAILED***");
     /* Class member variable m_terrain_map,
		denotes 2D height map.*/
    typedef boost::units::si::length height;
    boost::multi_array<boost::units::quantity<height,T>,2,boost::alignment::aligned_allocator<T,16>> m_terrain_map;
    
    terrain_topography_t()
      : 
    m_terrain_map(boost::extents[IME-IMS][JME-JMS])
    {
        constexpr T zero = static_cast<T>(0.0);
        for(std::size_t i = IMS; i<=IME; ++i)
        {
            for(std::size_t j = JMS; j<=JME; ++j)  
            {
                this->m_terrain_map[i][j] = boost::units::quantity<height,T>(zero*boost::units::si::meters);
            }
        }
    }
    
    // No error checking on the user passed argument!!
    terrain_topography_t(const boost::multi_array<T,2,boost::alignment::aligned_allocator<T,16>> &height_map) 
      : 
    m_terrain_map(height_map.shape()[0],height_map.shape()[1])
    {
        for(std::size_t i = IMS; i<=IME; ++i)
        {
            for(std::size_t j = JMS; j<=JME; ++j)
            {
                this->m_terrain_map[i][j] = boost::units::quantity<height, T>
						(height_map[i][j] * boost::units::si::meters);
            }
        }
    }

    terrain_topography_t(const terrain_topography_t &other) 
      : 
    m_terrain_map(other.m_terrain_map) 
    {

    }

    terrain_topography_t(terrain_topography_t &&other) 
      : 
    m_terrain_map(std::move(other.m_terrain_map)) 
    {

    }

    auto operator=(const terrain_topography_t &other)->terrain_topography_t &
    {
        if(__builtin_expect(this==&other,0)) { return (*this);}
        this->m_terrain_map.operator=(other.m_terrain_map);
        return (*this);
    }

    auto operator=(terrain_topography_t &&other)->terrain_topography_t &
    {
        if(__builtin_expect(this==&other,0)) { return (*this);}
        this->m_terrain_map.operator=(std::move(other.m_terrain_map));
        return (*this);
    }
 
    auto operator()(const std::size_t I,
                    const std::size_t J)const-> const T 
    {
        return (this->m_terrain_map[I][J].value());
    }

    auto operator==(const terrain_topography_t &other)const->const bool 
    {
        if(__builtin_expect(this->m_terrain_map.num_elements()!= other.m_terrain_map.num_elements(),0)) {return false;}
        return (this->m_terrain_map.operator==(other.m_terrain_map));
    }

    auto operator!=(const terrain_topography_t &other)const->const bool 
    {
        if(__builtin_expect(this->m_terrain_map.num_elements()!= other.m_terrain_map.num_elements(),0)) {return false;}
        return (this->m_terrain_map.operator!=(other.m_terrain_map));
    }

    auto operator+=(const boost::multi_array<T,2,boost::alignment::aligned_allocator<T,16>> &height_map)->terrain_topography_t &
    {
        if(__builtin_expect(this->m_terrain_map.num_elements()!=height_map.num_elements(),0)) {return terrain_topography_t();}
        for(std::size_t i = IMS; i<=IME; ++i)
        {
            for(std::size_t j = JMS; j<=JME; ++j)   
            {
                this->m_terrain_map[i][j].value() += height_map[i][j];
            }
        }
        return (*this);
    }

    auto get_terrain_map()const->boost::multi_array<boost::units::quantity<height,T>,2,boost::alignment::aligned_allocator<T,16>> const &
    {
        return (this->m_terrain_map);
    }

    auto increase_height(const std::size_t I,
                         const std::size_t J,
                         const T value)->terrain_topography_t &
    {   
        const bool I_dim = (I<IMS) || (I>IME);
        const bool J_dim = (J<JMS) || (J>JME);
        if(__builtin_expect(I_dim,0) || 
           __builtin_expect(J_dim,0)) { return (terrain_topography_t());}
        this->m_terrain_map[I][J].value() += value;
        return (*this);
    }

    auto increase_height(const std::size_t I0,
                         const std::size_t I1,
                         const std::size_t J0,
                         const std::size_t J1,
                         const boost::multi_array<T,2,boost::alignment::aligned_allocator<T,16>> &height_map)
    {
        const bool I_dim = (I0<IMS) || (I1>IME);
        const bool J_dim = (J0<JMS) || (J1>JME);
        if(__builtin_expect(I_dim,0) || 
           __builtin_expect(J_dim,0)) { return (terrain_topography_t());}
        for(std::size_t i = I0; i<I1; ++i)
        {
            for(std::size_t j = J0; j<J1; ++j)   
            {
                this->m_terrain_map[i][j].value() += height_map[i][j];
            }
        }
        return (*this);
    }

    auto operator<<(std::ostream &os)->std::ostream &
    {
        std::cout << typeid(*this).name() << " object state dump: " << std::endl;
        for(std::size_t i = IMS; i<=IME; ++i)
        {
            for(std::size_t j = JMS; j<=JME; ++j)   
            {
                os << std::fixed << std::setprecision(15) <<
						this->operator()(i, j) << std::endl;
            }
            os << "------------------" << std::endl;
        }
        return (os);
    }
};
                     

}

}
#endif /*__GMS_TERRAIN_TOPOGRAPHY_H__*/