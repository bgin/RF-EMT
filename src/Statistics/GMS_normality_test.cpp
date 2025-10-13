

#include <cmath>
#include <limits>
#include <cstdio>
#include "GMS_normality_test.h"

__ATTR_ALWAYS_INLINE__
static inline
double ppnd(double p)
{
constexpr double split{4.19999999999999984456877655248E-1};
constexpr double A0{2.50662823883999985596915394126E0};
constexpr double A1{-1.86150006252900013237194798421E1};
constexpr double A2{4.13911977353400004631112096831E1};
constexpr double A3{-2.54410604963699995550996391103E1};
constexpr double B1{-8.47351093089999984897531248862E0};
constexpr double B2{2.308336743742999885853350861E1};
constexpr double B3{-2.1062241018259999947304095258E1};
constexpr double B4{3.13082909833000000432434717368E0};
constexpr double C0{-2.78718931138000014513522728521E0};
constexpr double C1{-2.2979647913400000902583997231E0};
constexpr double C2{4.85014127135000006063592081773E0};
constexpr double C3{2.32121276858000014087224371906E0};
constexpr double D1{3.54388924762000012691487427219E0};
constexpr double D2{1.63706781897000008818565675028E0};
constexpr double ZERO{0.0};
constexpr double ONE{1.0};
constexpr double HALF{0.5};
double q,r,result;
q = p-HALF;
if(std::fabs(q)>split) { goto label10;}
    // 0.08 < P < 0.92
    r = q*q;
    result = q*(((A3*r+A2)*r+A1)*r+A0)/((((B4*r+B3)*r+B2)*r+B1)*r+ONE);
    return (result);
    // P < 0.08 OR P > 0.92, SET R = MIN(P,1-P)
label10:
    r = p;
    if(q>ZERO) {r = ONE-p;}
    if(r<=ZERO) {goto label20;}
    r = std::sqrt(-std::log(r));
    result = (((C3*r+C2)*r+C1)*r+C0)/((D2*r+D1)*r+ONE);
    if(q<ZERO) {result = -result;}
    return (result);
label20:
    result = 0.0;
    return (result);
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
                a[i] = ppnd((i-TH)/an25);
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
                fac = std::sqrt((summ2-TWO*t0-TWO*t1)/
                         (ONE-TWO*t0-TWO*t1) );
                a[0] = a1;
                a[1] = a2;
            }
            else 
            {
                i1 = 2;
                fac = std::sqrt((summ2-TWO*a[0]*a[0])/
                                (ONE-TWO*a[1]*a[1]));
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
    delta = (float)ncens/an;
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
    j = n-1;

    for(int32_t i = 1; i != n2; ++i) 
    {
        xi = x[i]/range;
        if(xx-xi>SMALL) {return;}
        sx += xi;
        if(i!=j)
        { 
            sa = sa+copysignf((float)1,(float)i-j)
                              *a[std::min(i,j)];
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
        xsx = x[i]/range;
        ssa = ssa+asa*asa;
        ssx = ssx+xsx*xsx;
        sax = sax+asa*xsx;
        j = j-1;
    }
    //
    //  W1 equals (1-W) claculated to avoid excessive rounding error
    //       for W very near 1 (a potential problem in very large samples)
    // 
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
