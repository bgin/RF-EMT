



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
!    RESULT is the approximate value of the integral

*/

#include <limits>
#include "GMS_simpne_quad.h"


void 
gms::math 
::simpne(const std::int32_t ntab,
         double * __restrict__ x,
         double * __restrict__ y,
         double &result)
{
   if(__builtin_expect(ntab<=2,0)) 
   {
      result = std::numeric_limits<double>::quiet_NaN();
      return;
   }
   double del[3];
   double g[3];
   double pi[3];
   double e,f,feints,sum1;
   double x1,x2,x3;
   std::int32_t i,n;

   result = 0.0;
   n      = 1;
   do 
   {
      x1 = x[n];
      x2 = x[n+1];
      x3 = x[n+2];
      e  = x3*x3-x1*x1;
      f  = x3*x3*x3-x1*x1*x1;
      feints = x3-x1;
      del[0] = x3-x2;
      del[1] = x1-x3;
      del[2] = x2-x1;
      g[0]   = x2+x3;
      g[1]   = x1+x3;
      g[2]   = x1+x2;
      pi[0]  = x2*x3;
      pi[1]  = x1*x3;
      pi[2]  = x1*x2;
      sum1   = 0.0;
      for(std::int32_t i{0}; i<3; ++i) 
      {
         sum1 = sum1+y[n-1+i]*del[i]*(f/3.0-g[i]*0.5*e+pi[i]*feints);
      }
      result = result-sum1/(del[0]*del[1]*del[2]);
      n = n+2;
   }while(ntab<=(n+1));
   const bool b = ((ntab%2) != 0);
   if(b) return;
   n = ntab-2;
   x3 = x[ntab];
   x2 = x[ntab-1];
   x1 = x[ntab-2];
   e = x3*x3-x2*x2;
   f = x3*x3*x3*-x2*x2*x2;
   feints = x3-x2;
   del[0] = x3-x2;
   del[1] = x1-x3;
   del[2] = x2-x1;
   g[0]   = x2+x3;
   g[1]   = x1+x3;
   g[2]   = x1+x2;
   pi[0]  = x2*x3;
   pi[1]  = x1*x3;
   pi[2]  = x1*x2;
   sum1 = 0.0;
   for(std::int32_t i{0}; i<3; ++i) 
   {
      sum1 = sum1+y[n-1+i]*del[i]*(f/3.0-g[i]*0.5*e+pi[i]*feints);
   }
   result = result-sum1/(del[0]*del[1]*del[2]);
}




 
