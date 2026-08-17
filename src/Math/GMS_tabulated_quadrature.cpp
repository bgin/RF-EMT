
#include <limits>
#include <cstdio>
#include "GMS_tabulated_quadrature.h"


/*
!*****************************************************************************80
!
!! HIORDQ approximates the integral of a function using equally spaced data.
!
!  Discussion:
!
!    The method applies the trapezoidal rule to various subsets of the
!    data, and then applies Richardson extrapolation.
!
!  Modified:
!
!    10 February 2006
!
!  Author:
!
!    Alan Kaylor Cline,
!    Department of Computer Science,
!    University of Texas at Austin.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, number of data points.
!
!    Input, real ( kind = 8 ) DELT, the spacing between the X values of the
!    data.  The actual X values are not needed!
!
!    Input, real ( kind = 8 ) Y(NTAB), the Y values of the data.
!
!    Work array, real ( kind = 8 ) WORK(2*(NTAB-1)).  The actual minimum amount
!    of workspace required is two times the number of integer
!    divisors of NTAB-1.
!
!    Output, real ( kind = 8 ) RESULT, the approximation to the integral.
!
*/

void 
gms::math
::hiordq(const std::int32_t ntab,
         const double delt,
         const double * __restrict  y,
         double * __restrict  work,
         double &result)
{
double fac,sum2,sum1;
std::int32_t i,j,jbak,jj,k;
std::int32_t ntab1,ntab1k;
ntab1 = ntab-1;         
sum1 = 0.5*(y[0]+y[ntab-1]);
j = -1;
for(k = 1; k != ntab1; ++k) 
{
    //ntab1k = ntab1/k;
    if(((ntab1/k)*k)==ntab1) 
    {
        sum2 = -sum1;
        for(i = 0; i < ntab; i += (ntab1/k)) 
        {
            sum2 += y[i];
        }
        j += 2;
        const double t0 = static_cast<double>(ntab1/k);
        work[j]   = delt*sum2*t0;
        work[j+1] = t0*t0;
        if(k != 1) 
        {
            for(jj = 3; jj < j; jj += 2) 
            {
                jbak = j+1-jj;
                fac = work[j+1]/(work[j+1]-work[jbak+1]);
                work[jbak] = work[jbak+2]+fac*(work[jbak]-work[jbak+2]);
            }
        }
    }
}
    result = work[1];
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

void gms::math
::simpn(const int32_t ntab,
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

/*
 *****************************************************************************80
!
!! WEDINT uses Weddle's rule to integrate data at equally spaced points.
!
!  Modified:
!
!    10 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, is the number of data points.  
!    (NTAB-1) must be divisible by 6.
!
!    Input, real ( kind = 8 ) H, is the spacing between the points at which
!    the data was evaluated.
!
!    Input, real ( kind = 8 ) FTAB(NTAB), contains the tabulated data values.
!
!    Output, real ( kind = 8 ) RESULT, is the approximation to the integral.   
*/

 void gms::math
 ::wedint(const int32_t ntab,
          const double h,
          const double * __restrict ftab,
          double &result) 
{

    if(__builtin_expect(ntab<=1,0)) 
    {
        result = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    if((ntab%6)!=0) 
    {
        result = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    std::int32_t i;
#pragma omp simd reduction(+:result)
    for(i = 0; i != ntab-5; i += 6)
    {
        result = result+ftab[i]+5.0*ftab[i+1]+
                 ftab[i+2]+6.0*ftab[i+3]     +
                 ftab[i+4]+5.0*ftab[i+5]     +
                 ftab[i+6];
    }
    result = 3.0*h*result/10.0;
}

/*
!*****************************************************************************80
!
!! AVINT estimates the integral of unevenly spaced data.
!
!  Discussion:
!
!    The data is given as NTAB pairs of values 
!    ( XTAB(1:NTAB), YTAB(1:NTAB) ).
!
!    The quadrature method uses overlapping parabolas and smoothing.
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
!    Paul Hennion,
!    Algorithm 77:
!    Interpolation, Differentiation and Integration,
!    Communications of the ACM,
!    Volume 5, page 96, 1962.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of entries in XTAB and
!    YTAB.  NTAB must be at least 2.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), the abscissas at which the
!    function values are given.  The XTAB's must be distinct
!    and in ascending order.
!
!    Input, real ( kind = 8 ) YTAB(NTAB), the function values,
!    YTAB(I) = F(XTAB(I)).
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.  A should
!    be, but need not be, near one endpoint of the interval
!    (X(1), X(NTAB)).
!
!    Input, real ( kind = 8 ) B, the upper limit of integration.  B should
!    be, but need not be, near one endpoint of the interval
!    (X(1), X(NTAB)).
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
!!
*/

void
gms::math
::avint(const std::int32_t ntab,
        const double * __restrict xtab,
        const double * __restrict ytab,  
        const double a,
        const double b,
        double &result)        
{

    double fa,fb,ca,cb,cc,slope,sum,syl,syl2,syl3,syu, 
           syu2 , syu3 , term1 , term2 , term3 , x1 , 
           x12 , x13 , x2 , x23 , x3, total,ba,bb,bc;
    std::int32_t i , inlft , inrt , istart , istop;
    result = 0.0;
    // Error check and trivial cases
    if(a==b) 
     {
        result =  std::numeric_limits<double>::quiet_NaN(); 
        return;
     }
    if(a>b) 
    {
        result =  std::numeric_limits<double>::quiet_NaN();
        return;
    }
    if(ntab<2) 
    {
        result =  std::numeric_limits<double>::quiet_NaN(); 
        return;
    }
   for(i = 1; i < ntab; ++i) 
   {
       if(xtab[i]<=xtab[i-1]) 
       {
          result =  std::numeric_limits<double>::quiet_NaN();
          return;
       }                    
    }
    // Special case fort NTAB=2     
    if(ntab==2)
    {
       slope = (ytab[1] - ytab[0]) / ( xtab[1]-xtab[0]);
       fa    =  ytab[0] + slope * ( a - xtab[0]);
       fb    =  ytab[1]  + slope * ( b - xtab[1] );
       result = 0.5 * ( fa + fb ) * ( b - a );
       return;
    }  
    if(xtab[ntab-2]<0.0 || b<xtab[3])
    {
        result = std::numeric_limits<double>::quiet_NaN(); 
        return;
    }   

    i = 0;
    for(; i<ntab; ++i) 
    {
       if(a<xtab[i]) 
       {
           result = std::numeric_limits<double>::quiet_NaN();
           return;

        }
    } 
    inlft = i;
    i = ntab;
    for(; i>0; --i) 
    {
       if(xtab[i]>b)
       {
        result = std::numeric_limits<double>::quiet_NaN();
        return; 
        }
    }  
    inrt = i;

    if(inrt-inlft<2) 
    {
        result = std::numeric_limits<double>::quiet_NaN();
       return;
    }
    if(inlft==1)
       istart = 2;
    else 
       istart = inlft;
    if(inrt==ntab)
       istop = ntab-1;
    else 
       istop = inrt;

    total = 0.0;
    syl = a;
    syl2 = syl*syl;
    syl3 = syl2*syl;
    for(i=istart; i<istop; ++i) 
    {
        x1 = xtab[i-1];
        x2 = xtab[i];
        x3 = xtab[i+1];
        x12 = x1-x2;
        x13 = x1-x3;
        x23 = x2-x3;
        term1 =  (ytab[i-1]/(x12*x13));
        term2 = -(ytab[i]/(x12*x23));
        term3 =  (ytab[i+1]/(x13*x23));
        ba = term1+term2+term3;
        bb = -( x2 + x3 ) * term1 - ( x1 + x3 ) * term2 - ( x1 + x2 ) * term3;
        bc =  x2 * x3 * term1 + x1 * x3 * term2 + x1 * x2 * term3;

        if(i=istart)
        {
            ca = ba;
            cb = bb;
            cc = bc;
        }
        else 
        {
            ca = 0.5*(ba+ca);
            cb = 0.5*(bb+cb);
            cc = 0.5*(bc+cc);
        }
        syu   = x2;
        syu2  = syu*syu;
        syu3  = syu2*syu;
        total = total + ca * ( syu3 - syl3 ) / 3.0D+00
                      + cb * ( syu2 - syl2 ) / 2.0D+00 
                      + cc * ( syu  - syl );
        ca = ba;  
        cb = bb;
        cc = bc;
        syl = syu;
        syl2 = syu2;
        syl3 = syu3;
    }
    syu = b;
    syu2 = syu * syu;
    syu3 = syu2 * syu;
    result = total + ca * ( syu3 - syl3 ) / 3.0D+00
                   + cb * ( syu2 - syl2 ) / 2.0D+00 
                   + cc * ( syu  - syl  );
}
