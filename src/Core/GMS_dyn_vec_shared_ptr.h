/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef __GMS_DYN_VEC_SHARED_PTR_H__
#define __GMS_DYN_VEC_SHARED_PTR_H__ 290620260600

#include <cstdint>
#include <memory>
#include <initializer_list>
#include <iterator>
#include <algorithm>
#include <utility>
#include <iostream>
#include "GMS_config.h"

namespace file_version 
{

    const unsigned int GMS_DYN_VEC_SHARED_PTR_MAJOR = 1U;
    const unsigned int GMS_DYN_VEC_SHARED_PTR_MINOR = 0U;
    const unsigned int GMS_DYN_VEC_SHARED_PTR_MICRO = 0U;
    const unsigned int GMS_DYN_VEC_SHARED_PTR_FULLVER =
      1000U*GMS_DYN_VEC_SHARED_PTR_MAJOR+
      100U*GMS_DYN_VEC_SHARED_PTR_MINOR+
      10U*GMS_DYN_VEC_SHARED_PTR_MICRO;
    const char * const GMS_DYN_VEC_SHARED_PTR_CREATION_DATE = "29-06-2026 06:00 AM +00200 (MON 29 JUN 2026 GMT+2)";
    const char * const GMS_DYN_VEC_SHARED_PTR_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_DYN_VEC_SHARED_PTR_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_DYN_VEC_SHARED_PTR_DESCRIPTION   = "Dynamic vector based on shared pointer memory management.";

}

namespace gms
{

namespace common 
{

template<typename T>
class gms_vec_sptr_t 
{

public:
using value_type      = T;
using size_type       = std::size_t;
using pointer         = T* __restrict;
using const_pointer   = const T* __restrict;
using reference       = T&;
using const_reference = const T&;

private:
std::unique_ptr<value_type[]> m_data;
size_type m_data_len{};
size_type m_capacity{};

public:
size_type size() const 
{
   return this->m_data_len;
}
size_type capacity() const 
{
    return this->m_capacity;
}
bool empty() const 
{
    return (this->size() == 0ull);
}
private:
bool full() const 
{
    return (this->size == this->capacity());
}

public:

using iterator       = pointer;
using const_iterator = const_pointer;
iterator begin()
{
    return (this->m_data.get());
}
const_iterator begin() const
{
    return (this->m_data.get());
}
const_iterator cbegin() const 
{
    return (this->begin());
}
iterator end()
{
    return (this->begin()+this->size());
}
const_iterator end() const 
{
    return (this->begin()+this->size());
}

void swap_wrapper(const gms_vec_sptr_t &other) noexcept 
{
    std::swap(this->m_data,other.m_data);
    std::swap(this->m_data_len,other.m_data_len);
    std::swap(this->m_capacity,other.m_capacity);
}

gms_vec_sptr_t() = default;

gms_vec_sptr_t(const size_type data_len)
:
m_data{new value_type[data_len]},
m_data_len{data_len},
m_capacity{data_len}
{

}

gms_vec_sptr_t(const size_type data_len,const_reference init_val)
:
m_data{new value_type[data_len]},
m_data_len{data_len},
m_capacity{data_len}
{
    std::fill(this->begin(),this->end(),init_val);
}

gms_vec_sptr_t(const gms_vec_sptr_t &other)
:
m_data{new value_type[other.size()]},
m_data_len{other.m_data_len},
m_capacity{other.m_capacity}
{
    std::copy(other.begin(),other.end(),this->begin());
}

gms_vec_sptr_t(const float * __restrict other_ptr,
               const std::size_t other_data_len)
:
m_data{new value_type[other_data_len]},
m_data_len{other_data_len},
m_capacity{other_data_len}
{
    std::copy(other_ptr,other_ptr+other_data_len,this->begin());
}

gms_vec_sptr_t(gms_vec_sptr_t &&other) = default;

gms_vec_sptr_t(std::initializer_list<T> list)
:
m_data{new value_type[list.size()]},
m_data_len{list.size()},
m_capacity{list.size()}
{
    std::copy(list.begin(),list.end(),this->begin());
}

~gms_vec_sptr_t() = default;

gms_vec_sptr_t & 
operator=(const gms_vec_sptr_t &other)
{
    gms_vec_sptr_t{other.swap(*this)};
    return (*this);
}

gms_vec_sptr_t &
operator=(gms_vec_sptr_t &&other) = default;

reference operator[](size_type nth_data_elem) 
{
    return (this->m_data[nth_data_elem]);
}

const_reference operator[](size_type nth_data_elem) const 
{
    return (this->m_data[nth_data_elem]);
}

reference front()
{
    return ((*this)[0]);
}

const_reference front() const 
{
    return ((*this)[0]);
}

reference back()
{
    return ((*this)[this->size()-1ull]);
}

const_reference back() const 
{
    return ((*this)[this->size()-1ull]);
}

bool 
operator==(const gms_vec_sptr_t &other) const 
{
    bool cond = (this->size()==other.size()) && 
                std::equal(this->begin(),this->end(),other.begin());
    return (cond);
}

bool 
operator!=(const gms_vec_sptr_t &other) const 
{
    return !(*this==other);
}

void push_back(const_reference value) 
{
    if(this->full())
    {
       this->grow();
    }
    this->m_data[this->size()] = value;
    ++this->m_data_len;
}

void push_back(T&& value) 
{
    if(this->full())
    {
       this->grow();
    }
    this->m_data[this->size()] = std::move(value);
    ++this->m_data_len;
}

template<class ... Args>
reference emplace_back(Args &&...args) 
{
    if(this->full())
    {
        this->grow();
    }
    this->m_data[this->size()] = value_type(std::forward<Args>(args)...);
    ++this->m_data_len;
    return (this->back());
}

private:

void grow() 
{
    this->resize(this->capacity()*2ull);
}
public:

void resize(const size_type new_capacity) 
{
    if(new_capacity<=this->capacity()) {return;}
    std::shared_ptr<value_type[]> sptr = std::make_unique<value_type[]>(new_capacity);
    if constexpr(std::is_nothrow_move_assignable_v<T>)
    {
        std::move(this->begin(),this->end(),sptr.get());
    }
    else 
    {
        std::copy(this->begin(),this->end(),sptr.get());
    }
    this->m_data.reset(sptr.release());
    this->m_capacity = new_capacity;
}

template<typename Iter>
iterator insert(const_iterator given_pos,Iter first,Iter last)
{
    iterator give_pos_c = const_cast<iterator>(given_pos);
    const std::size_t remaining = this->capacity()-this->size();
    const std::size_t n_dist    = std::distance(first,last);
    if(remaining<n_dist)
    {
        auto idx = std::distance(this->begin(),given_ps_c);
        this->resize(this->capacity()+n_dist-remaining);
        give_pos_c = std::next(this->begin(),idx);
    }
    std::copy_backward(given_pos_c,this->end(),this->end()+n_dist);
    std::copy(first,last,given_pos_c);
    this->m_data_len += n_dist;
    return (give_pos_c);
}

iterator erase(const_iterator given_pos) 
{
    iterator given_pos_c = const_cast<iterator>(given_pos);
    if (given_pos_c==this->end()) return given_pos_c;
    std::copy(std::next(given_pos_c),this->end(),given_pos_c);
    *std::prev(this->end()) = {};
    --this->m_data_len;
    return given_pos_c;
}

};

template <class T>
std::ostream& operator<<(std::ostream& os, 
                         const gms_vec_sptr_t<T>& vec) 
{
   if (vec.empty()) return os;
   os << vec.front();
   for (auto p = std::next(vec.begin()); p != vec.end(); ++p)
      os << ',' << *p;
   return os;
}


} // common

} // gms

#endif /*__GMS_DYN_VEC_SHARED_PTR_H__*/