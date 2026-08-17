

#ifndef __GMS_TABULATED_QUADRATURE_H__
#define __GMS_TABULATED_QUADRATURE_H__ 170820260619

#include <cstdint>
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_TABULATED_QUADRATURE_MAJOR = 1;
     static const unsigned int GMS_TABULATED_QUADRATURE_MINOR = 1;
     static const unsigned int GMS_TABULATED_QUADRATURE_MICRO = 0;
     static const unsigned int GMS_TABULATED_QUADRATURE_FULLVER =
       1000U*GMS_TABULATED_QUADRATURE_MAJOR+100U*GMS_TABULATED_QUADRATURE_MINOR+
       10U*GMS_TABULATED_QUADRATURE_MICRO;
     static const char GMS_TABULATED_QUADRATURE_CREATION_DATE[] = "17-08-2026 06:07AM +00200 (MON 17 AUG 2026 GMT+2)";
     static const char GMS_TABULATED_QUADRATURE_BUILD_DATE[]    = __DATE__; 
     static const char GMS_TABULATED_QUADRATURE_BUILD_TIME[]    = __TIME__;
     static const char GMS_TABULATED_QUADRATURE_SYNOPSIS[]      = "Tabulated data numerical integrators collection. Ported by Bernard Gingold from the Fortran 77 implementation, i.e. \
                                                                   intlib.f90 file (library, modified by J. Burkardt).";

}

/*
   Integrators expects an equally spaced input data!!
*/

namespace gms
{

namespace math 
{

__ATTR_HOT__
__ATTR_ALIGN__(32)			
void hiordq(const std::int32_t,
            const double ,
            const double * __restrict,
            double * __restrict,
            double &);

__ATTR_HOT__
__ATTR_ALIGN__(32)
void simpn(const std::int32_t ,
           const double ,
           double * __restrict ,
           double &);


__ATTR_HOT__
__ATTR_ALIGN__(32)			
void wedint(const std::int32_t ,
            const double ,
            const double * __restrict ,
            double &); 

__ATTR_HOT__
__ATTR_ALIGN__(32)			
void avint(  const std::int32_t ,
             const double * __restrict ,
             const double * __restrict ,
             const double,
             const double,
             double &);

}


}

#endif /*__GMS_TABULATED_QUADRATURE_H__*/