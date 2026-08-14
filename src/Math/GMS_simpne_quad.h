

#ifndef __GMS_SIMPNE_QUAD_H__
#define __GMS_SIMPNE_QUAD_H__ 081020221553


namespace file_info 
{


        const unsigned int GMS_SIMPNE_QUAD_MAJOR = 1;

	const unsigned int GMS_SIMPNE_QUAD_MINOR = 1;

	const unsigned int GMS_SIMPNE_QUAD_MICRO = 0;

	const unsigned int GMS_SIMPNE_QUAD_FULLVER = 
		1000U*GMS_SIMPNE_QUAD_MAJOR+100U*GMS_SIMPNE_QUAD_MINOR+10U*GMS_SIMPNE_QUAD_MICRO;

	const char * const GMS_SIMPNE_QUAD_CREATE_DATE = "08-10-2022 15:51 +00200 (SAT 08 OCT 2022 GMT+2)";

	const char * const GMS_SIMPNE_QUAD_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const GMS_SIMPNE_QUAD_AUTHOR = "Tabulated data integrator ported from Fortran version";

}


/*
!*****************************************************************************80
!
!! SIMPNE approximates the integral of unevenly spaced data.
!
!  Discussion:
!
!    The routine repeatedly interpolates a 3-point Lagrangian polynomial 
!    to the data and integrates that exactly.
!
!  Modified:
!
!    10 February 2006
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, number of data points.  NTAB must be at least 3.
!
!    Input, real ( kind = 8 ) X(NTAB), contains the X values of the data,
!    in order.
!
!    Input, real ( kind = 8 ) Y(NTAB), contains the Y values of the data.
!
!    Output, real ( kind = 8 ) RESULT.
!    RESULT is the approximate value of the integra

*/

#include <cstdint>
#include "GMS_config.h"


namespace  gms {

            namespace math {
                   
__ATTR_HOT__
__ATTR_ALIGN__(32)
void simpne(const int32_t ntab,
            double * __restrict  x,
            double * __restrict  y,
            double &result); 


         } // math

} // gms


#endif /*__GMS_SIMPNE_QUAD_H__*/
