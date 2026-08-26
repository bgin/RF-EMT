
#include <cmath> // bessel J0
#include <cstdio>
#include "GMS_cephes_double.h"
#include "GMS_integrands_func_ch5.h"


namespace
{
    /*1/SQRT(2*PI)*/
    static constexpr const double C03989422804014326779399460599344 = 0.3989422804014326779399460599344;
    /* 1/PI*/
    static constexpr const double C0318309886183790671537767526745  = 0.318309886183790671537767526745;
    /*0.5*PI*/
    static constexpr const double C15707963267948966192313216916398 = 1.5707963267948966192313216916398;

}

/*
   Two L1D and L1I issues will plague the execution of this function
   on its first invocation.
   The simpliest workaround is to issue a dummy call before starting
   the numerical integration.
*/

double 
gms
::fading_channel
::detail
::warmup_sin(const double x)
{
    return (std::sin(x));
}

double 
gms
::fading_channel
::detail
::warmup_cos(const double x)
{
    return (std::cos(x));
}

double 
gms
::fading_channel
::detail
::warmup_tan(const double x)
{
    return (std::tan(x));
}

double 
gms 
::fading_channel
::detail
::warmup_exp(double x)
{
    return (std::exp(x));
}

double 
gms 
::fading_channel
::detail
::warmup_cyl_bess_j(const std::int32_t n,const double x)
{
    return (std::cyl_bessel_j(n,x));
}

/* a -- modulation/detection threshold constant,
   gamma -- instantenous average SNR per bit
   theta -- integrand argument [0,0.5*PI]
*/
double 
gms::fading_channel
::integrand_Rayleigh_chan_5_6(const double theta,const double a,const double gamma)
{
#if (INTEGRANDS_FUNC_CH5_USE_CEPHES_DOUBLE) == 0
    const double asqr{a*a};
    const double sintht{std::sin(theta)};
    const double num{asqr*gamma};
    const double tmp{sintht*sintht};
    const double denom{tmp+tmp};
    const double ratio{1.0+(num/denom)};
    const double result{1.0/ratio};
    return (result);
#else 
    const double asqr{a*a};
    const double sintht{gms::math::cephes_d::sin(theta)};
    const double num{asqr*gamma};
    const double tmp{sintht*sintht};
    const double denom{tmp+tmp};
    const double ratio{1.0+(num/denom)};
    const double result{1.0/ratio};
    return (result);
#endif
}

double 
gms::fading_channel
::integrand_Rayleigh_chan_5_6_iface(const double theta,void * __restrict__ user_data)
{
    func_args_ch5_payload_t * __restrict__ p_payload{reinterpret_cast<func_args_ch5_payload_t* __restrict__>(user_data)};
    const double xa     = p_payload->arg1d;
    const double xgamma = p_payload->arg2d;
    return (integrand_Rayleigh_chan_5_6(theta,xa,xgamma));
}

/* a -- modulation/detection threshold constant,
   gamma -- instantenous average SNR per bit
   theta -- integrand argument [0,0.5*PI]
   q -- Nakagami-q (Hoyt) PDF parameter [0<q<1]
*/
double 
gms::fading_channel
::integrand_Hoyt_chan_5_9(const double theta,const double a,
                          const double gamma,const double q)
{
#if (INTEGRANDS_FUNC_CH5_USE_CEPHES_DOUBLE) == 0
    const double asqr{a*a};
    const double qsqr{q*q};
    const double gammasqr{gamma*gamma};
    const double sintht{std::sin(theta)};
    const double sinthtp2{sintht*sintht};
    const double tmp1{1.0+qsqr};
    const double den_fac1{tmp1*tmp1};
    const double left_ratio{1.0+(asqr/sinthtp2)*gamma};
    const double num{qsqr*asqr*gammasqr};
    const double den_fac2{den_fac1*(sinthtp2*sinthtp2)};
    const double right_ratio{num/den_fac2};
    const double sqrt_arg{left_ratio+right_ratio};
    const double result{1.0/std::sqrt(sqrt_arg)};
    return (result);
#else 
    const double asqr{a*a};
    const double qsqr{q*q};
    const double gammasqr{gamma*gamma};
    const double sintht{gms::math::cephes_d::sin(theta)};
    const double sinthtp2{sintht*sintht};
    const double tmp1{1.0+qsqr};
    const double den_fac1{tmp1*tmp1};
    const double left_ratio{1.0+(asqr/sinthtp2)*gamma};
    const double num{qsqr*asqr*gammasqr};
    const double den_fac2{den_fac1*(sinthtp2*sinthtp2)};
    const double right_ratio{num/den_fac2};
    const double sqrt_arg{left_ratio+right_ratio};
    const double result{1.0/gms::math::cephes_d::sqrt(sqrt_arg)};
    return (result);
#endif
}

double 
gms::fading_channel
::integrand_Hoyt_chan_5_9_iface(const double theta,void * __restrict__ user_data)
{
    func_args_ch5_payload_t * __restrict__ p_payload{reinterpret_cast<func_args_ch5_payload_t* __restrict__>(user_data)};
    const double xa = p_payload->arg1d;
    const double xgamma = p_payload->arg2d;
    const double xq = p_payload->arg3d;
    return (integrand_Hoyt_chan_5_9(theta,xa,xgamma,xq));
}

/* a -- modulation/detection threshold constant,
   gamma -- instantenous average SNR per bit
   theta -- integrand argument [0,0.5*PI]
   n -- Nakagami-n (Rice) PDF parameter [0<n<1]
*/
double 
gms::fading_channel
::integrand_Rice_chan_5_12(const double theta,const double a,
                           const double n,const double gamma)
{
#if (INTEGRANDS_FUNC_CH5_USE_CEPHES_DOUBLE) == 0
    const double nsqr{n*n};
    const double asqr{a*a};
    const double halfgamma{0.5*gamma};
    const double one_p_nsqr{1.0+nsqr};
    const double asqrgamm{asqr*halfgamma};
    const double sintht{std::sin(theta)};
    const double sinthtp2{sintht*sintht};
    const double num_den_fac{one_p_nsqr*sinthtp2};
    const double den_fac{num_den_fac+asqrgamm};
    const double ratio_left(one_p_nsqr/den_fac);
    const double ratio_right{((nsqr*asqrgamm)/den_fac)};
    const double exp_val{std::exp(-ratio_right)};
    const double result{ratio_left*exp_val};
    return (result);
#else 
    const double nsqr{n*n};
    const double asqr{a*a};
    const double halfgamma{0.5*gamma};
    const double one_p_nsqr{1.0+nsqr};
    const double asqrgamm{asqr*halfgamma};
    const double sintht{gms::math::cephes_d::sin(theta)};
    const double sinthtp2{sintht*sintht};
    const double num_den_fac{one_p_nsqr*sinthtp2};
    const double den_fac{num_den_fac+asqrgamm};
    const double ratio_left(one_p_nsqr/den_fac);
    const double ratio_right{((nsqr*asqrgamm)/den_fac)};
    const double exp_val{gms::math::cephes_d::exp(-ratio_right)};
    const double result{ratio_left*exp_val};
    return (result);
#endif 
}

double 
gms::fading_channel
::integrand_Rice_chan_5_12_iface(const double theta,void * __restrict__ user_data)
{
    func_args_ch5_payload_t * __restrict__ p_payload{reinterpret_cast<func_args_ch5_payload_t* __restrict__>(user_data)};
    const double xa = p_payload->arg1d;
    const double xn = p_payload->arg2d;
    const double xgamma = p_payload->arg3d;
    return (integrand_Rice_chan_5_12(theta,xa,xn,xgamma));
}

/* a -- modulation/detection threshold constant,
   gamma -- instantenous average SNR per bit
   theta -- integrand argument [0,0.5*PI]
   m -- Nakagami-m PDF parameter [m>=0.5]
*/
double 
gms::fading_channel
::integrand_Nakagami_m_chan_5_16(const double theta,const double a,
                                 const double m,const double gamma)
{
#if (INTEGRANDS_FUNC_CH5_USE_CEPHES_DOUBLE) == 0
    const double asqrgamm{a*a*gamma};
    const double sintht{std::sin(theta)};
    const double sinthtp2{sintht*sintht};
    const double denom{(m+m)*sinthtp2};
    const double inv_arg{1.0+(asqrgamm/denom)};
    const double invpow_val{1.0/std::pow(inv_arg)};
    return (invpow_val);
#else 
    const double asqrgamm{a*a*gamma};
    const double sintht{gms::math::cephes_d::sin(theta)};
    const double sinthtp2{sintht*sintht};
    const double denom{(m+m)*sinthtp2};
    const double inv_arg{1.0+(asqrgamm/denom)};
    const double invpow_val{1.0/gms::math::cephes_d::pow(inv_arg)};
    return (invpow_val);
#endif 
}

double 
gms::fading_channel
::integrand_Nakagami_m_chan_5_16_iface(const double theta,void * __restrict__ user_data)
{
    func_args_ch5_payload_t * __restrict__ p_payload{reinterpret_cast<func_args_ch5_payload_t* __restrict__>(user_data)};
    const double xa   = p_payload->arg1d;
    const double xm   = p_payload->arg2d;
    const double xgamma = p_payload->arg3d;
    return (integrand_Nakagami_m_chan_5_16(theta,xa,xm,xgamma));
}

/* a -- modulation/detection threshold constant,
   mu -- 10*log10(gamma) [db] i.e. instantenous average SNR per bit
   sigma -- logarithmic standev of shadowing [db]
   theta -- outer integral integration variable [0,0.5*PI]
*/
double 
gms::fading_channel
::integrand_Log_Norm_Shadow_chan_5_20(const double x,const double theta,
                                      const double a,const double mu,
                                      const double sigma)
{
#if (INTEGRANDS_FUNC_CH5_USE_CEPHES_DOUBLE) == 0
    const double asqr{a*a};
    const xsqr{x*x};
    const double sintht{std::sin(theta)};
    const double sinthtp2{sintht*sintht};
    const double exp_val2{std::exp(-xsqr)};
    const double exp_arg1{asqr/(sinthtp2+sinthtp2)};
    const double tmp{x*std::sqrt(sigma+sigma)+mu};
    const double ten_arg{tmp*0.1};
    const double ten_to_pow{std::pow(10.0,ten_arg)};
    const double exp_val1{std::exp(-exp_arg1*ten_to_pow)};
    const double result(exp_val1*exp_val2);
    return (result);
#else 
    const double asqr{a*a};
    const xsqr{x*x};
    const double sintht{gms::math::cephes_d::sin(theta)};
    const double sinthtp2{sintht*sintht};
    const double exp_val2{gms::math::cephes_d::exp(-xsqr)};
    const double exp_arg1{asqr/(sinthtp2+sinthtp2)};
    const double tmp{x*gms::math::cephes_d::sqrt(sigma+sigma)+mu};
    const double ten_arg{tmp*0.1};
    const double ten_to_pow{gms::math::cephes_d::pow(10.0,ten_arg)};
    const double exp_val1{gms::math::cephes_d::exp(-exp_arg1*ten_to_pow)};
    const double result(exp_val1*exp_val2);
    return (result);
#endif 
}

double 
gms::fading_channel
::integrand_Log_Norm_Shadow_chan_5_20_iface(const double x,void * __restrict__ user_data)
{
    func_args_ch5_payload_t * __restrict__ p_payload{reinterpret_cast<func_args_ch5_payload_t* __restrict__>(user_data)};
    const double xtheta = p_payload->arg1d;
    const double xa     = p_payload->arg2d;
    const double xmu    = p_payload->arg3d;
    const double xsigma = p_payload->arg4d;
    return (integrand_Log_Norm_Shadow_chan_5_20(x,xtheta,xa,xmu,xsigma));
}



