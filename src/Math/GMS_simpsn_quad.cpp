



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



#include <limits>
#include "GMS_simpsn_quad.h"



void gms::math::simpn(const int32_t ntab,
                      const double h,
                      double * __restrict  y,
                      double &result) 
{
    if(__builtin_expect(ntab<=2,0)) 
    {
      result = std::numeric_limits<double>::quiet_NaN();
      return;
    }
    double del[3];
    double pi[3];
    double g[3];
    double f,e,sum1,f3,h2;
    int32_t i,n;
    bool b;
    result = 0.0;           
    b = ((ntab%2)==0);
    if(b)
      n = ntab-1;
    else
      n = ntab;
    result = y[0]+y[n]+4.0*y[n-1];
    for(i = 1; i != n-2; i += 2) 
    {
      const double yi = y[i];
      const double yi1= y[i+1];
      result = result+4.0*yi+2.0*yi1;
    }
    result = h*result/3.0;
    const bool b2 = ((ntab%2)==1);
    if(b2) { return;}
    f      = h*h*h;
    del[0] = h;
    del[1] = -2.0*h;
    del[2] = h;
    g[0]   = h;
    g[1]   = 0.0;
    g[2]   = -h;
    pi[0]  = 0.0;
    pi[1]  = -h*h;
    pi[2]  = 0.0;
    n      = n-1;
    sum1 = 0.0;
    for(i = 0; i<3; ++i) 
    {
       sum1 = sum1+y[n-1+i]*del[i]*(f/3.0-g[i]*0.5*h*h+pi[i]*h);
    }
    result = result+0.5*sum1/(h*h*h);
}
