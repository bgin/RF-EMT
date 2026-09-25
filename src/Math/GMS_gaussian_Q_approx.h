
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

#ifndef __GMS_GAUSSIAN_Q_APPROX_H__
#define __GMS_GAUSSIAN_Q_APPROX_H__ 160920261342

#include <cstdint>
#include "GMS_config.h"
#include "GMS_fast_simd_funcs_approx.h"

namespace file_info 
{

     static const unsigned int GMS_GAUSSIAN_Q_APPROX_MAJOR = 1;
     static const unsigned int GMS_GAUSSIAN_Q_APPROX_MINOR = 1;
     static const unsigned int GMS_GAUSSIAN_Q_APPROX_MICRO = 0;
     static const unsigned int GMS_GAUSSIAN_Q_APPROX_FULLVER =
       1000U*GMS_GAUSSIAN_Q_APPROX_MAJOR+100U*GMS_GAUSSIAN_Q_APPROX_MINOR+
       10U*GMS_GAUSSIAN_Q_APPROX_MICRO;
     static const char GMS_GAUSSIAN_Q_APPROX_CREATION_DATE[] = "16-09-2026 13:45AM +00200 (WED 16 AUG 2026 GMT+2)";
     static const char GMS_GAUSSIAN_Q_APPROX_BUILD_DATE[]    = __DATE__; 
     static const char GMS_GAUSSIAN_Q_APPROX_BUILD_TIME[]    = __TIME__;
     static const char GMS_GAUSSIAN_Q_APPROX_SYNOPSIS[]      = "Gaussian Q function approximations, based on the paper: Gaussian Q Function Approximation in Wireless \
                                                                Communication System Design: A \
                                                                Gradient-Based Optimization Approach ";

}



namespace gms 
{

namespace math 
{

__ATTR_ALWAYS_INLINE__
static inline
float gaussian_Q_approx_chiani(const float x)
{
    constexpr float C0083333333333333333333333333 = 0.083333333333333333333333333;
    constexpr float C05                           = 0.5f;
    constexpr float C03                           = 0.333333333333333333333333333;
    constexpr float C025                          = 0.25f;
    const float xx            = x*x;
    const float left_exp_val  = C0083333333333333333333333333*expapprox(-C05*xx);
    const float right_exp_val = C025*expapprox(-C03*(2.0*xx));
    return (left_exp_val+right_exp_val);
}

__ATTR_ALWAYS_INLINE__
static inline
double gaussian_Q_approx_chiani(const double x)
{
    constexpr double C0083333333333333333333333333 = 0.083333333333333333333333333;
    constexpr double C05                           = 0.5;
    constexpr double C03                           = 0.333333333333333333333333333;
    constexpr double C025                          = 0.25;
    const double xx            = x*x;
    const double left_exp_val  = C0083333333333333333333333333*expapprox_d(-C05*xx);
    const double right_exp_val = C025*expapprox_d(-C03*(2.0*xx));
    return (left_exp_val+right_exp_val);
}

__ATTR_ALWAYS_INLINE__
static inline 
float gaussian_Q_approx_loskot_2T(const float x)
{
    constexpr float C0208  = 0.208f;
    constexpr float C0971  = 0.971f;
    constexpr float C0147  = 0.147f;
    constexpr float C0525  = 0.525f;
    const float xx = x*x;
    return (__builtin_fmaf(C0208,expapprox(-C0971*xx),C0147*expapprox(-C0525*xx)));

}

__ATTR_ALWAYS_INLINE__
static inline 
double gaussian_Q_approx_loskot_2T(const double x)
{
    constexpr double C0208  = 0.208;
    constexpr double C0971  = 0.971;
    constexpr double C0147  = 0.147;
    constexpr double C0525  = 0.525;
    const double xx = x*x;
    return (__builtin_fma(C0208,expapprox_d(-C0971*xx),C0147*expapprox_d(-C0525*xx)));
}

__ATTR_ALWAYS_INLINE__
static inline 
float gaussian_Q_approx_loskot_3T(const float x)
{
    constexpr float C0168 = 0.168f;
    constexpr float C0876 = 0.876f;
    constexpr float C0144 = 0.144f;
    constexpr float C0525 = 0.525f;
    constexpr float C0002 = 0.002f;
    constexpr float C0603 = 0.603f;
    const float xx = x*x;

    return (__builtin_fmaf(C0168,expapprox(-C0876*xx),
                           __builtin_fmaf(C0144,expapprox(-C0525*xx),C0002*expapprox(-C0603*xx)))); 
}

__ATTR_ALWAYS_INLINE__
static inline 
double gaussian_Q_approx_loskot_3T(const double x)
{
    constexpr double C0168 = 0.168;
    constexpr double C0876 = 0.876;
    constexpr double C0144 = 0.144;
    constexpr double C0525 = 0.525;
    constexpr double C0002 = 0.002;
    constexpr double C0603 = 0.603;
    const double xx = x*x;
    return (__builtin_fma(C0168,expapprox_d(-C0876*xx),
                          __builtin_fma(C0144,expapprox_d(-C0525*xx),C0002*expapprox_d(-C0603*xx))));

}

__ATTR_ALWAYS_INLINE__
static inline 
float gaussian_Q_approx_sadhwani_1T(const float x)
{
    constexpr float C025 = 0.25f;
    constexpr float C05  = 0.5f;
    const float xx = x*x;
    return (C025*expapprox(-C05*xx));
}

__ATTR_ALWAYS_INLINE__
static inline 
double gaussian_Q_approx_sadhwani_1T(const double x)
{
    constexpr double C025 = 0.25;
    constexpr double C05  = 0.5;
    const double xx = x*x;
    return (C025*expapprox_d(-C05*xx));
}

__ATTR_ALWAYS_INLINE__
static inline 
float gaussian_Q_approx_sadhwani_2T(const float x)
{
    constexpr float C0125 = 0.125f;
    constexpr float C025  = 0.25f;
    constexpr float C05   = 0.5f;
    const float xx  = x*x;
    return (__builtin_fmaf(C0125,expapprox(-xx),C025*expapprox(-C05*xx)));
}

__ATTR_ALWAYS_INLINE__
static inline 
double gaussian_Q_approx_sadhwani_2T(const double x)
{
    constexpr double C0125 = 0.125;
    constexpr double C025  = 0.25;
    constexpr double C05   = 0.5;
    const double xx = x*x;
    return (__builtin_fma(C0125,expapprox_d(-xx),C025*expapprox_d(-C05*xx)));
}

/*
gaussian_Q_approx_loskot_3T(double):
 vmulsd xmm10,xmm0,xmm0
 sub    rsp,0x48
 vmovsd xmm3,QWORD PTR [rip+0x1408]        # 403088 <_IO_stdin_used+0x88>
 vmulsd xmm0,xmm10,QWORD PTR [rip+0x1488]        # 403110 <_IO_stdin_used+0x110>
 vmulsd xmm0,xmm0,QWORD PTR [rip+0x13e8]        # 403078 <_IO_stdin_used+0x78>
 vaddsd xmm0,xmm0,QWORD PTR [rip+0x13e8]        # 403080 <_IO_stdin_used+0x80>
 vcomisd xmm3,xmm0
 jbe    401ef0 <gaussian_Q_approx_loskot_3T(double)+0x280>
 vxorpd xmm1,xmm1,xmm1
 vmovsd xmm3,QWORD PTR [rip+0x13e2]        # 403090 <_IO_stdin_used+0x90>
 movabs rdx,0xfffffffffffff
 movabs rcx,0x3ff0000000000000
 vmaxsd xmm0,xmm0,xmm1
 vmovsd xmm9,QWORD PTR [rip+0x13ca]        # 403098 <_IO_stdin_used+0x98>
 vmovsd xmm8,QWORD PTR [rip+0x13ca]        # 4030a0 <_IO_stdin_used+0xa0>
 vmovsd xmm7,QWORD PTR [rip+0x13ca]        # 4030a8 <_IO_stdin_used+0xa8>
 vmovsd xmm6,QWORD PTR [rip+0x13ca]        # 4030b0 <_IO_stdin_used+0xb0>
 vmovsd xmm5,QWORD PTR [rip+0x13ca]        # 4030b8 <_IO_stdin_used+0xb8>
 vmovsd xmm4,QWORD PTR [rip+0x13ca]        # 4030c0 <_IO_stdin_used+0xc0>
 vcvttsd2si rax,xmm0
 and    rdx,rax
 or     rdx,rcx
 vmovq  xmm0,rdx
 movabs rdx,0x7ff0000000000000
 vmulsd xmm2,xmm0,xmm3
 and    rax,rdx
 vmovq  xmm1,rax
 vsubsd xmm2,xmm2,xmm9
 vmulsd xmm2,xmm2,xmm0
 vaddsd xmm2,xmm2,xmm8
 vmulsd xmm2,xmm2,xmm0
 vaddsd xmm2,xmm2,xmm7
 vmulsd xmm2,xmm2,xmm0
 vaddsd xmm2,xmm2,xmm6
 vmulsd xmm2,xmm2,xmm0
 vaddsd xmm2,xmm2,xmm5
 vmulsd xmm2,xmm2,xmm0
 vmulsd xmm0,xmm10,QWORD PTR [rip+0x13a2]        # 4030f0 <_IO_stdin_used+0xf0>
 vmulsd xmm0,xmm0,QWORD PTR [rip+0x1322]        # 403078 <_IO_stdin_used+0x78>
 vaddsd xmm0,xmm0,QWORD PTR [rip+0x1322]        # 403080 <_IO_stdin_used+0x80>
 vaddsd xmm2,xmm2,xmm4
 vmulsd xmm2,xmm2,xmm1
 vmulsd xmm2,xmm2,QWORD PTR [rip+0x13aa]        # 403118 <_IO_stdin_used+0x118>
 vxorpd xmm1,xmm1,xmm1
 vmovsd QWORD PTR [rsp],xmm4
 movabs rdx,0xfffffffffffff
 movabs rcx,0x3ff0000000000000
 vmaxsd xmm0,xmm0,xmm1
 vmovsd QWORD PTR [rsp+0x38],xmm10
 vmovsd QWORD PTR [rsp+0x30],xmm3
 vmovsd QWORD PTR [rsp+0x28],xmm9
 vcvttsd2si rax,xmm0
 vmovsd QWORD PTR [rsp+0x20],xmm8
 vmovsd QWORD PTR [rsp+0x18],xmm7
 vmovsd QWORD PTR [rsp+0x10],xmm6
 vmovsd QWORD PTR [rsp+0x8],xmm5
 and    rdx,rax
 or     rdx,rcx
 vmovq  xmm1,rdx
 movabs rdx,0x7ff0000000000000
 vmulsd xmm0,xmm1,xmm3
 and    rax,rdx
 vsubsd xmm0,xmm0,xmm9
 vmulsd xmm0,xmm0,xmm1
 vaddsd xmm0,xmm0,xmm8
 vmulsd xmm0,xmm0,xmm1
 vaddsd xmm0,xmm0,xmm7
 vmulsd xmm0,xmm0,xmm1
 vaddsd xmm0,xmm0,xmm6
 vmulsd xmm0,xmm0,xmm1
 vaddsd xmm0,xmm0,xmm5
 vmulsd xmm0,xmm0,xmm1
 vmovsd xmm1,QWORD PTR [rip+0x1314]        # 403120 <_IO_stdin_used+0x120>
 vaddsd xmm0,xmm0,xmm4
 vmovq  xmm4,rax
 vmulsd xmm0,xmm0,xmm4
 call   401030 <fma@plt>
 vmovsd xmm10,QWORD PTR [rsp+0x38]
 vmovsd xmm4,QWORD PTR [rsp]
 vmovsd xmm5,QWORD PTR [rsp+0x8]
 vmovsd xmm6,QWORD PTR [rsp+0x10]
 vmovapd xmm2,xmm0
 vmulsd xmm10,xmm10,QWORD PTR [rip+0x12e7]        # 403128 <_IO_stdin_used+0x128>
 vmovsd xmm7,QWORD PTR [rsp+0x18]
 vmulsd xmm10,xmm10,QWORD PTR [rip+0x1229]        # 403078 <_IO_stdin_used+0x78>
 vmovsd xmm8,QWORD PTR [rsp+0x20]
 vaddsd xmm1,xmm10,QWORD PTR [rip+0x1223]        # 403080 <_IO_stdin_used+0x80>
 vmovsd xmm9,QWORD PTR [rsp+0x28]
 vmovsd xmm3,QWORD PTR [rsp+0x30]
 movabs rdx,0xfffffffffffff
 vxorpd xmm0,xmm0,xmm0
 movabs rcx,0x3ff0000000000000
 vmaxsd xmm1,xmm1,xmm0
 vcvttsd2si rax,xmm1
 and    rdx,rax
 or     rdx,rcx
 vmovq  xmm1,rdx
 movabs rdx,0x7ff0000000000000
 vmulsd xmm3,xmm1,xmm3
 and    rax,rdx
 vsubsd xmm3,xmm3,xmm9
 vmulsd xmm3,xmm3,xmm1
 vaddsd xmm3,xmm3,xmm8
 vmulsd xmm3,xmm3,xmm1
 vaddsd xmm3,xmm3,xmm7
 vmulsd xmm3,xmm3,xmm1
 vaddsd xmm3,xmm3,xmm6
 vmulsd xmm3,xmm3,xmm1
 vaddsd xmm3,xmm3,xmm5
 vmovq  xmm5,rax
 vmulsd xmm3,xmm3,xmm1
 vaddsd xmm3,xmm3,xmm4
 vmulsd xmm0,xmm3,xmm5
 vmovsd xmm1,QWORD PTR [rip+0x124b]        # 403130 <_IO_stdin_used+0x130>
 add    rsp,0x48
 jmp    401030 <fma@plt>
 xchg   ax,ax
 vmulsd xmm0,xmm10,QWORD PTR [rip+0x11f8]        # 4030f0 <_IO_stdin_used+0xf0>
 vmulsd xmm0,xmm0,QWORD PTR [rip+0x1178]        # 403078 <_IO_stdin_used+0x78>
 vaddsd xmm0,xmm0,QWORD PTR [rip+0x1178]        # 403080 <_IO_stdin_used+0x80>
 vcomisd xmm3,xmm0
 ja     401f68 <gaussian_Q_approx_loskot_3T(double)+0x2f8>
 vmovsd xmm2,QWORD PTR [rip+0x115a]        # 403070 <_IO_stdin_used+0x70>
 vmovsd xmm1,QWORD PTR [rip+0x1202]        # 403120 <_IO_stdin_used+0x120>
 vmovsd QWORD PTR [rsp],xmm10
 vmovapd xmm0,xmm2
 call   401030 <fma@plt> <---------------------------------------------------------------- Inserted trampoline jump into the @.plt section!!
 vmovsd xmm10,QWORD PTR [rsp]
 vmulsd xmm10,xmm10,QWORD PTR [rip+0x11ef]        # 403128 <_IO_stdin_used+0x128>
 vmulsd xmm10,xmm10,QWORD PTR [rip+0x1137]        # 403078 <_IO_stdin_used+0x78>
 vmovsd xmm3,QWORD PTR [rip+0x113f]        # 403088 <_IO_stdin_used+0x88>
 vmovapd xmm2,xmm0
 vaddsd xmm1,xmm10,QWORD PTR [rip+0x112b]        # 403080 <_IO_stdin_used+0x80>
 vcomisd xmm3,xmm1
 ja     401fad <gaussian_Q_approx_loskot_3T(double)+0x33d>
 vmovsd xmm0,QWORD PTR [rip+0x110d]        # 403070 <_IO_stdin_used+0x70>
 jmp    401edd <gaussian_Q_approx_loskot_3T(double)+0x26d>
 vmovsd xmm2,QWORD PTR [rip+0x1100]        # 403070 <_IO_stdin_used+0x70>
 vmovsd xmm3,QWORD PTR [rip+0x1118]        # 403090 <_IO_stdin_used+0x90>
 vmovsd xmm9,QWORD PTR [rip+0x1118]        # 403098 <_IO_stdin_used+0x98>
 vmovsd xmm8,QWORD PTR [rip+0x1118]        # 4030a0 <_IO_stdin_used+0xa0>
 vmovsd xmm7,QWORD PTR [rip+0x1118]        # 4030a8 <_IO_stdin_used+0xa8>
 vmovsd xmm6,QWORD PTR [rip+0x1118]        # 4030b0 <_IO_stdin_used+0xb0>
 vmovsd xmm5,QWORD PTR [rip+0x1118]        # 4030b8 <_IO_stdin_used+0xb8>
 vmovsd xmm4,QWORD PTR [rip+0x1118]        # 4030c0 <_IO_stdin_used+0xc0>
 jmp    401d6e <gaussian_Q_approx_loskot_3T(double)+0xfe>
 vmovsd xmm3,QWORD PTR [rip+0x10db]        # 403090 <_IO_stdin_used+0x90>
 vmovsd xmm9,QWORD PTR [rip+0x10db]        # 403098 <_IO_stdin_used+0x98>
 vmovsd xmm8,QWORD PTR [rip+0x10db]        # 4030a0 <_IO_stdin_used+0xa0>
 vmovsd xmm7,QWORD PTR [rip+0x10db]        # 4030a8 <_IO_stdin_used+0xa8>
 vmovsd xmm6,QWORD PTR [rip+0x10db]        # 4030b0 <_IO_stdin_used+0xb0>
 vmovsd xmm5,QWORD PTR [rip+0x10db]        # 4030b8 <_IO_stdin_used+0xb8>
 vmovsd xmm4,QWORD PTR [rip+0x10db]        # 4030c0 <_IO_stdin_used+0xc0>
 jmp    401e69 <gaussian_Q_approx_loskot_3T(double)+0x1f9>
 nop    WORD PTR [rax+rax*1+0x0]
*/

__ATTR_ALWAYS_INLINE__
static inline 
float gaussian_Q_approx_sadhwani_4T(const float x)
{
    constexpr float C1 = 0.0625f;
    constexpr float C2 = 0.125f;
    constexpr float C3 = 0.5f;
    constexpr float C4 = 0.33333333333333333333333f;
    constexpr float C5 = 0.058823529411764705882352941f;
    const float xx = x*x;
    const float tenxx = 10.0f*xx;
    return (__builtin_fmaf(C1,expapprox(-C3*xx),
                __builtin_fmaf(C2,expapprox(-xx),
                    __builtin_fmaf(C2,expapprox(-C4*tenxx),C2*expapprox(-C5*tenxx)))));
}

__ATTR_ALWAYS_INLINE__
static inline 
double gaussian_Q_approx_sadhwani_4T(const double x)
{
    constexpr double C1 = 0.0625;
    constexpr double C2 = 0.125;
    constexpr double C3 = 0.5;
    constexpr double C4 = 0.33333333333333333333333;
    constexpr double C5 = 0.058823529411764705882352941;
    const double xx = x*x;
    const double tenxx = 10.0*xx;
    return (__builtin_fma(C1,expapprox_d(-C3*xx),
                __builtin_fma(C2,expapprox_d(-xx),
                    __builtin_fma(C2,expapprox_d(-C4*tenxx),C2*expapprox_d(-C5*tenxx)))));
}

/*
    Bound-based/exponential polynomial based approximations
*/

__ATTR_ALWAYS_INLINE__
static inline
float gaussian_Q_approx_cooper(const float x)
{
    constexpr float C1 = 0.5f;
    constexpr float C2 = 6.283185307179586476925286767f;
    const float xx     = x*x;
    const float sqrt_2pix = std::sqrt(C2*x);
    const float right_term= 1.0f-(1.0f/(xx+xx));
    const float exp_term  = expapprox(-C1*xx);
    const float ratio     = exp_term/sqrt_2pix;
    return (ratio*right_term);
}

__ATTR_ALWAYS_INLINE__
static inline
double gaussian_Q_approx_cooper(const double x)
{
    constexpr double C1 = 0.5;
    constexpr double C2 = 6.283185307179586476925286767;
    const double xx     = x*x;
    const double sqrt_2pix = std::sqrt(C2*x);
    const double right_term= 1.0-(1.0/(xx+xx));
    const double exp_term  = expapprox_d(-C1*xx);
    const double ratio     = exp_term/sqrt_2pix;
    return (ratio*right_term);
}

__ATTR_ALWAYS_INLINE__
static inline 
float gaussian_Q_approx_borjesson(const float x)
{
    constexpr float C1 = 0.5f;
    constexpr float C2 = 2.506628274631000502415765285f;
    const float xx     = x*x;
    const float exp_val= C2*expapprox(-C1*xx);
    const float sqrt_term  = 0.339f*std::sqrt(xx+5.510f);
    const float right_term = 1.0f/(__builtin_fmaf(0.661f,x,sqrt_term));
    return (exp_val*right_term);
}

__ATTR_ALWAYS_INLINE__
static inline 
double gaussian_Q_approx_borjesson(const double x)
{
    constexpr double C1 = 0.5;
    constexpr double C2 = 2.506628274631000502415765285;
    const double xx     = x*x;
    const double exp_val= C2*expapprox_d(-C1*xx);
    const double sqrt_term  = 0.339*std::sqrt(xx+5.510);
    const double right_term = 1.0/(__builtin_fma(0.661,x,sqrt_term));
    return (exp_val*right_term);
}

__ATTR_ALWAYS_INLINE__
static inline
float gaussian_Q_approx_sadhwani_summed(const float x,const std::int32_t n) 
{
if(__builtin_expect(n>128,0)) { return (-1.0f);}
constexpr const float theta_lut[128] = 
{
26561.07370058031301596,
2951.52672978394593883,
1062.76301068144402961,
542.38938216885594557,
328.24391161267163852,
219.84388820771997075,
157.49790211713116150,
118.38163379128440056,
92.23973597198599350,
73.90969946516314337,
60.56287935539454281,
50.54400438065793111,
42.83209252697486136,
36.76965720093080847,
31.91778586717441257,
27.97440695989148551,
24.72611246534243534,
22.01866773467314076,
19.73836125017757581,
17.79986866518916599,
16.13815476591664932,
14.70294783067476274,
13.45489318988675542,
12.36282713537543643,
11.40181318331400817,
10.55170643021994792,
9.79608972017488178,
9.12147551075540264,
8.51670021771579755,
7.97245976035782089,
7.48094990150482975,
7.03558520545476096,
6.63077757000577961,
6.26176032547270278,
5.92444749260821357,
5.61532039144371531,
5.33133569098266236,
5.06985038851206493,
4.82856024770480818,
4.60544900516174227,
4.39874624527744107,
4.20689229309963242,
4.02850882009587341,
3.86237412456220319,
3.70740225596646766,
3.56262531497762280,
3.42717838884132009,
3.30028668303544404,
3.18125449075313638,
3.06945570625878084,
2.96432564001850762,
2.86535393539436001,
2.77207842067833932,
2.68407975793799203,
2.60097677280667572,
2.52242236796785813,
2.44809993843393503,
2.37772021942528200,
2.31101850820766996,
2.24775221004035153,
2.18769866574001037,
2.13065322453304118,
2.07642753105719757,
2.02484799975171637,
1.97575445357971491,
1.92899890717015055,
1.88444447714107866,
1.84196440464714883,
1.80144117714498853,
1.76276573804219461,
1.72583677433227445,
1.69056007355504945,
1.65684794248971023,
1.62461868091121531,
1.59379610454104137,
1.56430911201842116,
1.53609129132299893,
1.50908056160717741,
1.48321884685699712,
1.45845177820343341,
1.43472842205929152,
1.41200103156710211,
1.39022481911630291,
1.36935774792838338,
1.34936034092082524,
1.33019550524813823,
1.31182837108428330,
1.29422614335786368,
1.27735796528209100,
1.26119479263764811,
1.24570927786994146,
1.23087566315441310,
1.21666968166584666,
1.20306846636115816,
1.19005046565099137,
1.17759536539444132,
1.16568401670416888,
1.15429836909675965,
1.14342140856596108,
1.13303710019499548,
1.12313033495892256,
1.11368688039941000,
1.10469333488267329,
1.09613708517704822,
1.08800626710994641,
1.08028972908511101,
1.07297699826029014,
1.06605824920296688,
1.05952427485770095,
1.05336645967322906,
1.04757675475074774,
1.04214765488702521,
1.03707217739714941,
1.03234384261201928,
1.02795665595514785,
1.02390509151210085,
1.02018407701397673,
1.01678898016385832,
1.01371559624214091,
1.01096013693317843,
1.00851922032179320,
1.00638986201394287,
1.00456946734126218,
1.00305582461434328,
1.00184709939451322,
1.00094182975856194,
1.00033892253539203,
1.00003765049793447
};

   const float inv2n = 1.0/static_cast<float>(n+n);
   const float xx    = x*x;
   const float C1    = 0.5f;
   float gauss_Q_sum = 0.0f;
   for(std::int32_t j = 1;j <= n; ++j)  
   {
      const float theta_lut_val = theta_lut[j];
      const float exp_arg       = -theta_lut_val*xx*C1;
      const float exp_val       = expapprox(exp_arg);
      gauss_Q_sum               += exp_val;
   }
   return (inv2n*gauss_Q_sum);
}

__ATTR_ALWAYS_INLINE__
static inline
double gaussian_Q_approx_sadhwani_summed(const double x,const std::int32_t n) 
{
if(__builtin_expect(n>128,0)) { return (-1.0);}
constexpr const double theta_lut[128] = 
{
26561.07370058031301596,
2951.52672978394593883,
1062.76301068144402961,
542.38938216885594557,
328.24391161267163852,
219.84388820771997075,
157.49790211713116150,
118.38163379128440056,
92.23973597198599350,
73.90969946516314337,
60.56287935539454281,
50.54400438065793111,
42.83209252697486136,
36.76965720093080847,
31.91778586717441257,
27.97440695989148551,
24.72611246534243534,
22.01866773467314076,
19.73836125017757581,
17.79986866518916599,
16.13815476591664932,
14.70294783067476274,
13.45489318988675542,
12.36282713537543643,
11.40181318331400817,
10.55170643021994792,
9.79608972017488178,
9.12147551075540264,
8.51670021771579755,
7.97245976035782089,
7.48094990150482975,
7.03558520545476096,
6.63077757000577961,
6.26176032547270278,
5.92444749260821357,
5.61532039144371531,
5.33133569098266236,
5.06985038851206493,
4.82856024770480818,
4.60544900516174227,
4.39874624527744107,
4.20689229309963242,
4.02850882009587341,
3.86237412456220319,
3.70740225596646766,
3.56262531497762280,
3.42717838884132009,
3.30028668303544404,
3.18125449075313638,
3.06945570625878084,
2.96432564001850762,
2.86535393539436001,
2.77207842067833932,
2.68407975793799203,
2.60097677280667572,
2.52242236796785813,
2.44809993843393503,
2.37772021942528200,
2.31101850820766996,
2.24775221004035153,
2.18769866574001037,
2.13065322453304118,
2.07642753105719757,
2.02484799975171637,
1.97575445357971491,
1.92899890717015055,
1.88444447714107866,
1.84196440464714883,
1.80144117714498853,
1.76276573804219461,
1.72583677433227445,
1.69056007355504945,
1.65684794248971023,
1.62461868091121531,
1.59379610454104137,
1.56430911201842116,
1.53609129132299893,
1.50908056160717741,
1.48321884685699712,
1.45845177820343341,
1.43472842205929152,
1.41200103156710211,
1.39022481911630291,
1.36935774792838338,
1.34936034092082524,
1.33019550524813823,
1.31182837108428330,
1.29422614335786368,
1.27735796528209100,
1.26119479263764811,
1.24570927786994146,
1.23087566315441310,
1.21666968166584666,
1.20306846636115816,
1.19005046565099137,
1.17759536539444132,
1.16568401670416888,
1.15429836909675965,
1.14342140856596108,
1.13303710019499548,
1.12313033495892256,
1.11368688039941000,
1.10469333488267329,
1.09613708517704822,
1.08800626710994641,
1.08028972908511101,
1.07297699826029014,
1.06605824920296688,
1.05952427485770095,
1.05336645967322906,
1.04757675475074774,
1.04214765488702521,
1.03707217739714941,
1.03234384261201928,
1.02795665595514785,
1.02390509151210085,
1.02018407701397673,
1.01678898016385832,
1.01371559624214091,
1.01096013693317843,
1.00851922032179320,
1.00638986201394287,
1.00456946734126218,
1.00305582461434328,
1.00184709939451322,
1.00094182975856194,
1.00033892253539203,
1.00003765049793447
};

   const double inv2n = 1.0/static_cast<double>(n+n);
   const double xx    = x*x;
   const double C1    = 0.5;
   double gauss_Q_sum = 0.0;
   for(std::int32_t j = 1;j <= n; ++j)  
   {
      const double theta_lut_val = theta_lut[j];
      const double exp_arg       = -theta_lut_val*xx*C1;
      const double exp_val       = expapprox_d(exp_arg);
      gauss_Q_sum                += exp_val;
   }
   return (inv2n*gauss_Q_sum);
}

} //math

} //gms

#endif /*__GMS_GAUSSIAN_Q_APPROX_H__*/