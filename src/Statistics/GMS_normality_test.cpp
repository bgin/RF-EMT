

#include <cmath>
#include <limits>
#include <cstdio>
#include "GMS_normality_test.h"

__ATTR_ALWAYS_INLINE__
static inline
void ppnd7(const float p,
           float &normal_dev,
           std::int32_t &ifault)
{
    constexpr float zero   = 0.0f;
    constexpr float one    = 1.0f;
    constexpr float half   = 0.5f;
    constexpr float split1 = 0.425f;
    constexpr float split2 = 5.0f;
    constexpr float const1 = 0.180625;
    constexpr float const2 = 1.6f;
    float q,r;
    // Coefficients for P close to 0.5
    constexpr float a0 = 3.3871327179E+00f, a1 = 5.0434271938E+01f,
                    a2 = 1.5929113202E+02f, a3 = 5.9109374720E+01f,
                    b1 = 1.7895169469E+01f, b2 = 7.8757757664E+01f,
                    b3 = 6.7187563600E+01f;
    // ! Coefficients for P not close to 0, 0.5 or 1.
    constexpr float c0 = 1.4234372777E+00f, c1 = 2.7568153900E+00f, 
                    c2 = 1.3067284816E+00f, c3 = 1.7023821103E-01f, 
                    d1 = 7.3700164250E-01f, d2 = 1.2021132975E-01f;
    // Coefficients for P near 0 or 1.
    constexpr float e0 = 6.6579051150E+00f, e1 = 3.0812263860E+00f, 
                    e2 = 4.2868294337E-01f, e3 = 1.7337203997E-02f, 
                    f1 = 2.4197894225E-01f, f2 = 1.2258202635E-02f;

    ifault = 0;
    q = p - half;
    if(std::abs(q) <= split1)
    {
       r = const1 - q * q;
       normal_dev = q * (((a3 * r + a2) * r + a1) * r + a0) / 
                    (((b3 * r + b2) * r + b1) * r + one);
       return;
    }
    else 
    {
      if(q < zero)
          r = p;
      else
          r = one - p;
      if(r <= zero)
      {
         ifault = 1;
         normal_dev = zero;
         return;
      }
      r = std::sqrt(-std::log(r));
      if(r <= split2)
      {
        r = r - const2;
        normal_dev = (((c3 * r + c2) * r + c1) * r + c0) / ((d2 * r + d1) * r + one);
      }
      else 
      {
        r = r - split2;
        normal_dev = (((e3 * r + e2) * r + e1) * r + e0) / ((f2 * r + f1) * r + one);
      }
      if(q < zero) normal_dev = -normal_dev;
      return;
    }
}

__ATTR_ALWAYS_INLINE__
static inline
float poly(const float * __restrict c,
           const int32_t nord,
           const float x) {
float p,result;
int32_t n2,j;
result = c[0];
if(nord==1) { return std::numeric_limits<float>::quiet_NaN();}
p = x*c[nord];
if(nord==2) {goto label20;}
n2 = nord-2;
j = n2+1;
for(int32_t i = 0; i != n2; ++i) 
{
     p = (p+c[j])*x;
     j = j-1;
}
label20:
result += p;
return (result);
}

__ATTR_ALWAYS_INLINE__
static inline
double alnorm(const double x,
              const bool upper) 
{
    /*
     Evaluates the tail area of the standardised normal curve
c       from x to infinity if upper is .true. or
c       from minus infinity to x if upper is .false.
    */
constexpr double zero{0.0e+0};
constexpr double one{1.0e+0};
constexpr double half{0.5e+0};
constexpr double ltone{7.0e+0};
constexpr double utzero{18.66e+0};
constexpr double con{1.28e+00};
constexpr double p{3.98942280443999985894976134659E-1};
constexpr double q{3.99903485040000006289773182289E-1};
constexpr double r{3.98942280385000014319984984468E-1};
constexpr double a1{5.75885480458000031944720831234E0};
constexpr double a2{2.62433121678999992099079463514E0};
constexpr double a3{5.92885724437999961367040668847E0};
constexpr double b1{-2.98213557806999993715635355329E1};
constexpr double b2{4.86959930692000000362895661965E1};
constexpr double c1{-3.80520000000000021871322027772E-8};
constexpr double c2{3.98064793999999973986986256946E-4};
constexpr double c3{-1.51679116634999999746469256934E-1};
constexpr double c4{4.83859128080000022720241759089E0};
constexpr double c5{7.4238092402699995542292299433E-1};
constexpr double c6{3.99019417011000010475640920049E0};
constexpr double d1{1.00000615301999995487847172626E0};
constexpr double d2{1.98615381364000009867254448181E0};
constexpr double d3{5.29330324926000006513504558825E0};
constexpr double d4{-1.51508972450999994663334291545E1};
constexpr double d5{3.07899330340000005890033207834E1};
double z,y,result;
bool up;
z = x;
if(z>=zero) {goto label_10;}
up = !up;
z = -z;
label_10:
if(z<=ltone || (up && z <= utzero)) {goto label_20;}
result = 0.0;
goto label_40;
label_20:
y = half*z*z;
if(z>con) {goto label_30;}
result = half-z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))));
goto label_40;
label_30:
result = std::exp(-y)/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))));
label_40:
if(!up) {result = result-one;}
return (result);
}


void        
gms::math::shapiro_wilk(bool init,
                        float * __restrict x,
                        const int32_t n,
                        const int32_t n1,
                        const int32_t n2,
                        float * __restrict a,
                        float &w,
                        float &pw,
                        int32_t &ifault) 
{
    // Automatic constants
    const float C1[6] = {0.0,0.221157E0,-0.147981E0,
                         -0.207119E1,0.4434685E1,-0.2706056E1};
    const float C2[6] = {0.0,0.42981E-1,-0.293762E0,-0.1752461E1,
                        0.5682633E1, -0.3582633E1};
    const float C3[4] = {0.5440E0,-0.39978E0,0.25054E-1,-0.6714E-3};
    const float C4[4] = {0.13822E1,-0.77857E0,0.62767E-1,-0.20322E-2};
    const float C5[4] = {-0.15861E1,-0.31082E0,-0.83751E-1,0.38915E-2};
    const float C6[3] = {-0.4803E0, -0.82676E-1, 0.30302E-2};
    const float C7[2] = {0.164E0, 0.533E0};
    const float C8[2] = {0.1736E0, 0.315E0};
    const float C9[2] = {0.256E0, -0.635E-2};
    const float G[2]  = {-0.2273E1, 0.459E0};
    constexpr float Z90 = 0.12816E1;
    constexpr float Z95 = 0.16449E1;
    constexpr float Z99 = 0.23263E1;
    constexpr float ZM  = 0.17509E1;
    constexpr float ZSS = 0.56268E0;
    constexpr float BF1 = 0.8378E0;
    constexpr float XX90 = 0.556E0;
    constexpr float XX95 = 0.622E0;
    constexpr float ZERO = 0.0E0;
    constexpr float ONE  = 1.0E0;
    constexpr float TWO  = 2.0E0;
    constexpr float THREE = 3.0E0;
    constexpr float SQRTH = 0.70711E0;
    constexpr float QTR = 0.25E0;
    constexpr float TH = 0.375E0;
    constexpr float SMALL = 1.0E-19;
    constexpr float PI6 = 0.1909859E1;
    constexpr float STQR = 0.1047198E1;
    constexpr bool upper = true;
    float summ2,ssumm2,fac,rsn,an,an25,a1,a2,delta,range;
    float sa,sx,ssx,ssa,sax,asa,xsx,ssassx,w1,y,xx,xi;
    float gamma,m,s,ld,bf,z90f,z95f,z99f,zfm,zsd,zbar;
    int32_t ncens,nn2,i1,j;
    pw = ONE;
    if(w>=ZERO) {w = ONE;}
    an = n;
    ifault = 3;
    nn2 = n/2;
    if(n2<nn2) { return;}
    ifault = 1;
    if(n<3) { return;}

    // If INIT is false, calculates coefficients for the test
    if(!init) 
    {
        if(n==3) 
        {
            a[0] = SQRTH;
        }
        else 
        {
            an25 = an+QTR;
            summ2 = ZERO;

            for(int32_t i = 1; i != n2; ++i) 
            {
                ppnd7((i-TH)/an25,a[i],ifault);
                summ2 += a[i] * a[i];
            }
            summ2 *= TWO;
                      
            ssumm2 = std::sqrt(summ2);
            rsn = ONE/std::sqrt(an);
            a1 = poly(C1,6,rsn)-a[0]/ssumm2;
            //         Normalize coefficients
            if(n>5) 
            {
                i1 = 3;
                a2 = -a[1]/ssumm2+poly(C2,6,rsn);
                const float t0 = a[0]*a[0];
                const float t1 = a[1]*a[1];
                float arg1 = (summ2-TWO*t0-TWO*t1)/
                              (ONE-TWO*t0-TWO*t1);
                if(arg1<0.0f) arg1 *= -1.0f;
                fac = std::sqrt(arg1);
                //fac = std::sqrt(-1.0f*((summ2-TWO*t0-TWO*t1)/
                //         (ONE-TWO*t0-TWO*t1)) );
                a[0] = a1;
                a[1] = a2;
            }
            else 
            {
                i1 = 2;
                float arg2 = (summ2-TWO*a[0]*a[0])/
                             (ONE-TWO*a[1]*a[1]);
                if(arg2<0.0f) arg2 *= -1.0f;
                fac = std::sqrt(arg2);
                //fac = std::sqrt(-1.0f*((summ2-TWO*a[0]*a[0])/
                //                (ONE-TWO*a[1]*a[1])));
                a[0] = a1;
            }
            
            for(int32_t i = i1; i != nn2; ++i) 
            {
                const float t = -a[i];
                a[i] = t/fac;
            }
        }
        init = true;
    }
    if(n1<3) {return;}
    ncens = n-n1;
    ifault = 4;
    const bool b0 = ncens>0 && n<20;
    if(ncens<0 || b0) { return;}
    ifault = 5;
    delta = ((float)ncens)/an;
    if(delta>0.8f) {return;}
    //  If W input as negative, calculate significance level of -W
    if(w<ZERO) {
        w1 = ONE+w;
        ifault = 0;
        goto label_70;
    }
    //  Check for zero range
    ifault = 6;
    range = x[n1]-x[0];
    if(range<SMALL) {return;}
    // Check for correct sort order on range - scaled X
    ifault = 7;
    xx = x[0]/range;
    sx = xx;
    sa = -a[0];
    //std::printf("sa=%.7f\n",sa);
    j = n-1;
    
    for(int32_t i = 1; i != n2; ++i) 
    {
        xi = x[i]/range;
        if(xx-xi>SMALL) {return;}
        sx += xi;
        if(i!=j)
        { 
            //std::printf("copysignf=%.7f\n",copysignf((float)1,(float)i-j));
            //std::printf("a[std::min(i,j)]=%.7f\n",a[std::min(i,j)]);
            sa = sa+copysignf((float)1,(float)i-j)
                             *a[std::min(i,j)];
            //std::printf("sa=%.7f\n",sa);
        }
        xx = xi;
        j = j-1;
    }
    ifault = 0;
    if(n>5000) { ifault=2;}
    //
    // Calculate W statistic as squared correlation
    //        between data and coefficients
    //
    
    sa = sa/n1;
    //std::printf("sa=%.7f\n",sa);
    sx = sx/n1;
    ssa = ZERO;
    ssx = ZERO;
    sax = ZERO;
    j = n;
//#pragma loop_count min(2),avg(1000),max(5000)  

    for(int32_t i = 1; i != n1; ++i) {
        if(i!=j) 
        {
            asa = std::copysignf((float)1,(float)i-j)*
                                 a[std::min(i,j)]-sa;
        }
        else 
        {
            asa = -sa;
        }
        xsx = x[i]/range-sx;
        ssa = ssa+asa*asa;
        //std::printf("ssa=%.7f\n",ssa);
        ssx = ssx+xsx*xsx;
        sax = sax+asa*xsx;
        j = j-1;
    }
    //
    //  W1 equals (1-W) claculated to avoid excessive rounding error
    //       for W very near 1 (a potential problem in very large samples)
    // 
    //std::printf("%.17f\n",ssa*ssx);
    ssassx = std::sqrt(ssa*ssx);
    w1 = (ssassx-sax)*(ssassx+sax)/(ssa*ssx);
label_70:
    w = ONE-w1;
    //
    // Calculate significance level for W (exact for N=3)
    //
    if(n==3) 
    {
        pw = PI6*(std::asin(std::sqrt(w))-STQR);
        return;
    }
    y = std::log(w1);
    xx = std::log(an);
    m = ZERO;
    s = ONE;
    if(n<=11) 
    {
        gamma = poly(G,2,an);
        if(y>=gamma) {
            pw=SMALL;
            return;
        }
        y = -std::log(gamma-y);
        m = poly(C3,4,an);
        s = std::exp(poly(C4,4,an));
    }
    else 
    {
        m = poly(C5,4,xx);
        s = std::exp(poly(C6,3,xx));
    }
    if(ncens>0) 
    {
        // Censoring by proportion NCENS/N.  Calculate mean and sd
        //       of normal equivalent deviate of W.
        ld = -log(delta);
        bf = ONE+xx*BF1;
        z90f = Z90+bf*std::pow(poly(C7,2,std::pow(XX90,xx)),ld);
        z95f = Z95+bf*std::pow(poly(C8,2,std::pow(XX95,xx)),ld);
        z99f = Z99+bf*std::pow(poly(C9,2,xx),ld);
        //  Regress Z90F,...,Z99F on normal deviates Z90,...,Z99 to get
        // pseudo-mean and pseudo-sd of z as the slope and intercept
        zfm = (z90f+z95f+z99f)/THREE;
        zsd = (Z90*(z90f-zfm)+Z95*(z95f-zfm)+Z99*(z99f-zfm))/ZSS;
        zbar = zfm-zsd*ZM;
        m = m+zbar*s;
        s = s*zsd;
    }
    pw = alnorm((y-m)/s,upper);
} 
