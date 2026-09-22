/*
 * Copyright (C) Bernard Gingold, 2020-2026 
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

#ifndef __GMS_MARCUM_Q_APPROX_H__
#define __GMS_MARCUM_Q_APPROX_H__ 220920260715 


#include <cstdint>
#include "GMS_config.h"
#include "GMS_gaussian_Q_approx.h"

namespace file_info 
{

     static const unsigned int GMS_MARCUM_Q_APPROX_MAJOR = 1;
     static const unsigned int GMS_MARCUM_Q_APPROX_MINOR = 1;
     static const unsigned int GMS_MARCUM_Q_APPROX_MICRO = 0;
     static const unsigned int GMS_MARCUM_Q_APPROX_FULLVER =
       1000U*GMS_MARCUM_Q_APPROX_MAJOR+100U*GMS_MARCUM_Q_APPROX_MINOR+
       10U*GMS_MARCUM_Q_APPROX_MICRO;
     static const char GMS_MARCUM_Q_APPROX_CREATION_DATE[] = "16-09-2026 13:45AM +00200 (WED 16 AUG 2026 GMT+2)";
     static const char GMS_MARCUM_Q_APPROX_BUILD_DATE[]    = __DATE__; 
     static const char GMS_MARCUM_Q_APPROX_BUILD_TIME[]    = __TIME__;
     static const char GMS_MARCUM_Q_APPROX_AUTHOR[]        = "Bernard Gingold, beniekg@gmail.com";
     static const char GMS_MARCUM_Q_APPROX_SYNOPSIS[]      = "Marcum-Q function approximations, based on the Wikipedia article: https://en.wikipedia.org/wiki/Marcum_Q-function#Asymptotic_forms ";

}

namespace gms
{

namespace math 
{

/* Using Gaussian-Q approx by Chiani method*/
__ATTR_ALWAYS_INLINE__
static inline
double marcum_Q_approx_chiani(const double mu,const double a, const double b)
{
    const double pow_term = std::pow(b/a,mu-0.5);
    if(b>a)  
       return (pow_term*gaussian_Q_approx_chiani(b-a));
    else if(a>b)
       return (1.0-pow_term*gaussian_Q_approx_chiani(a-b));    
}

/* Using Gaussian-Q approx by Loskot method*/
__ATTR_ALWAYS_INLINE__
static inline
double marcum_Q_approx_loskot_2T(const double mu,const double a, const double b)
{
    const double pow_term = std::pow(b/a,mu-0.5);
    if(b>a)  
       return (pow_term*gaussian_Q_approx_loskot_2T(b-a));
    else if(a>b)
       return (1.0-pow_term*gaussian_Q_approx_loskot_2T(a-b));    
}

}

}

#endif /*__GMS_MARCUM_Q_APPROX_H__*/