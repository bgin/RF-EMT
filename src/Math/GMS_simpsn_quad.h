
#ifndef __GMS_SIMPSN_QUAD_H__
#define __GMS_SIMPSN_QUAD_H__ 091020221053


namespace file_info {


        const unsigned int GMS_SIMPN_QUAD_MAJOR = 1;

	const unsigned int GMS_SIMPN_QUAD_MINOR = 1;

	const unsigned int GMS_SIMPN_QUAD_MICRO = 0;

	const unsigned int GMS_SIMPN_QUAD_FULLVER = 
		1000U*GMS_SIMPN_QUAD_MAJOR+100U*GMS_SIMPN_QUAD_MINOR+10U*GMS_SIMPN_QUAD_MICRO;

	const char * const GMS_SIMPN_QUAD_CREATE_DATE = "09-10-2022 10:53 +00200 (SUN 09 OCT 2022 GMT+2)";

	const char * const GMS_SIMPN_QUAD_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const GMS_SIMPN_QUAD_AUTHOR = "Tabulated data integrator ported from Fortran version";

}


/*

!*****************************************************************************80
!
!! SIMPSN approximates the integral of evenly spaced data.
!
!  Discussion:
!
!    Simpson's rule is used.
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
!    Input, integer ( kind = 4 ) NTAB, the number of data points.  NTAB must be at least 3.
!
!    Input, real ( kind = 8 ) H, specifies the increment between the
!    X values.  Note that the actual X values are not needed,
!    just the constant spacing!
!
!    Input, real ( kind = 8 ) Y(NTAB), the data.
!
!    Output, real ( kind = 8 ) RESULT, the value of the integral
!    from the first to the last point.
*/


#include <cstdint>
#include "GMS_config.h"


namespace  gms {
      
           namespace  math {
                 
                      
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
                        void simpn(const int32_t ntab,
                                   const double h,
                                   double * __restrict y,
                                   double &result);
    } // math

} // gms


#endif /*__GMS_SIMPSN_QUAD_H__*/
