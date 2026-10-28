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

#ifndef __GMS_ODE_TO_STATE_SPACE_H__
#define __GMS_ODE_TO_STATE_SPACE_H__ 281020260600

#include <Eigen/Dense>
#include <array>
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_ODE_TO_STATE_SPACE_MAJOR = 1;
     static const unsigned int GMS_ODE_TO_STATE_SPACE_MINOR = 1;
     static const unsigned int GMS_ODE_TO_STATE_SPACE_MICRO = 0;
     static const unsigned int GMS_ODE_TO_STATE_SPACE_FULLVER =
       1000U*GMS_ODE_TO_STATE_SPACE_MAJOR+100U*GMS_ODE_TO_STATE_SPACE_MINOR+
       10U*GMS_ODE_TO_STATE_SPACE_MICRO;
     static const char GMS_ODE_TO_STATE_SPACE_CREATION_DATE[]  = "24-10-2026 07:45PM +00200 (SAT 24 OCT 2026 GMT+2)";
     static const char GMS_ODE_TO_STATE_SPACE_BUILD_DATE[]    = __DATE__; 
     static const char GMS_ODE_TO_STATE_SPACE_BUILD_TIME[]    = __TIME__;
     static const char GMS_ODE_TO_STATE_SPACE_AUTHOR[]        = "Programmer: Bernard Gingold, beniekg@gmail.com";   
     static const char GMS_ODE_TO_STATE_SPACE_SYNOPSIS[]      = "Ordinary Differential Equation to State-Space model converting class.";
}


namespace gms
{

namespace state_space
{


struct alignas(64) ode_to_state_space_t
{
    Eigen::Index m_norder;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_A;
    Eigen::VectorXd                                                        m_B;
    Eigen::RowVectorXd                                                     m_C;
    double                                                                 m_D;

    __ATTR_ALWAYS_INLINE__ 
    inline ode_to_state_space_t() noexcept(true)
    :
    m_norder{},
    m_A(),
    m_B(),
    m_C(),
    m_D{}
    {

    }

    __ATTR_ALWAYS_INLINE__ 
    inline ode_to_state_space_t(const Eigen::Index norder) noexcept(false)
    :
    m_norder(norder),
    m_A(m_norder,m_norder),
    m_B(m_norder),
    m_C(m_norder),
    m_D{}
    {

    }

    __ATTR_ALWAYS_INLINE__ 
    inline ode_to_state_space_t(const ode_to_state_space_t &other) noexcept(false)
    :
    m_norder(other.m_norder),
    m_A(other.m_A),
    m_B(other.m_B),
    m_C(other.m_C),
    m_D(other.m_D)
    {

    }

    __ATTR_ALWAYS_INLINE__ 
    inline ode_to_state_space_t(ode_to_state_space_t &&other) noexcept(true) 
    :
    m_norder(other.m_norder),
    m_A(std::move(other.m_A)),
    m_B(std::move(other.m_B)),
    m_C(std::move(other.m_C)),
    m_D(other.m_D)
    {

    }

    __ATTR_ALWAYS_INLINE__
    inline ~ode_to_state_space_t() noexcept(true) 
    {

    }

    __ATTR_ALWAYS_INLINE__ 
    ode_to_state_space_t & 
    operator=(const ode_to_state_space_t &other) noexcept(false) 
    {
        if(__builtin_expect(this==&other,0)) { return (*this);}
        this->m_norder = other.m_norder;
        this->m_A.operator=(other.m_A);
        this->m_B.operator=(other.m_B);
        this->m_C.operator=(other.m_C);
        this->m_D      = other.m_D;

        return (*this);
    }

    __ATTR_ALWAYS_INLINE__ 
    ode_to_state_space_t & 
    operator=(ode_to_state_space_t &&other) noexcept(true) 
    {
        if(__builtin_expect(this==&other,0)) { return (*this);}
        this->m_norder = other.m_norder;
        this->m_A.operator=(std::move(other.m_A));
        this->m_B.operator=(std::move(other.m_B));
        this->m_C.operator=(std::move(other.m_C));
        this->m_D      = other.m_D;

        return (*this);
    }

    __ATTR_ALWAYS_INLINE__ 
    inline
    const double * __restrict__ get_m_A_raw_ptr() const
    {
        return (this->m_A.data());
    } 

    __ATTR_ALWAYS_INLINE__ 
    inline 
    const double * __restrict__ get_m_B_raw_ptr() const 
    {
        return (this->m_B.data());
    }

    __ATTR_ALWAYS_INLINE__ 
    inline 
    const double * __restrict__ get_m_C_raw_ptr() const 
    {
        return (this->m_C.data());
    }

    // True if not-equal, i.e. the a.size() != b.size()+1
    template<std::size_t a_order,std::size_t b_order>
    __ATTR_ALWAYS_INLINE__ 
    static
    inline bool check_ode_order(const std::array<double,a_order> &a,
                                const std::array<double,b_order> &b) 
    {
        bool ret_status = false;
        constexpr std::size_t an = a.size();
        constexpr std::size_t bn = b.size();
        const is_bn_neq          = (bn != (an+1ull));
        if(__builtin_expect(is_bn_neq==true,0)) { ret_status = true; return (ret_status);}
    } 

    template<Eigen::Index a_order,Eigen::Index b_order>
    __ATTR_ALWAYS_INLINE__ 
    static
    inline bool check_ode_order(const double a[a_order],
                                const double b[b_order])                      
    {
        bool ret_status = false;
        const is_bn_neq          = (b_order != (a_order+1ull));
        if(__builtin_expect(is_bn_neq==true,0)) { ret_status = true; return (ret_status);}
    }

    template<std::size_t a_order,std::size_t b_order>
    __ATTR_ALWAYS_INLINE__ 
    std::int32_t convert_ode_to_ss(const std::array<double,a_order> &a,
                                   const std::array<double,b_order> &b) 
    {
        std::int32_t ret_status = -9999;
        this->m_D               = b.operator[](this->m_norder);
        for(Eigen::Index i = 0; i<this->m_norder-1; ++i) 
        {
            this->m_A.operator[](i,i+1) = 1.0;
        }
        for(Eigen::Index j = 0; j<this->m_norder; ++j)   
        {
            this->m_A.operator[](this->m_norder-1,j) = -a.operator[](j);
        }
        const double bn = b.operator[](this->m_norder);
        for(Eigen::Index i = 0; i<this->m_norder; ++i) 
        {   
            const double b_i = b.operator[](i);
            const double a_i = a.operator[](i);
            this->m_C.operator[](i) = b_i-bn*a_i;
        }
        this->m_B.operator[](this->m_norder-1) = 1.0;
        ret_status = 1;
        return (ret_status);
    }

    template<std::size_t a_order,std::size_t b_order>
    __ATTR_ALWAYS_INLINE__ 
    std::int32_t convert_ode_to_ss(const double a[a_order],
                                   const double b[b_order]) 
    {
        std::int32_t ret_status = -9999;
        this->m_D               = b[this->m_norder];
        for(Eigen::Index i = 0; i<this->m_norder-1; ++i) 
        {
            this->m_A.operator[](i,i+1) = 1.0;
        }
        for(Eigen::Index j = 0; j<this->m_norder; ++j)   
        {
            this->m_A.operator[](this->m_norder-1,j) = -a[j];
        }
        const double bn = b[this->m_norder];
        for(Eigen::Index i = 0; i<this->m_norder; ++i) 
        {   
            const double b_i = b[i];
            const double a_i = a[i];
            this->m_C.operator[](i) = b_i-bn*a_i;
        }
        this->m_B.operator[](this->m_norder-1) = 1.0;
        ret_status = 1;
        return (ret_status);
    }

};

}

}

#endif /*__GMS_ODE_TO_STATE_SPACE_H__*/