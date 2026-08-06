
#include <cmath> // bessel J0
#include <cstdio>
#include "GMS_cephes_double.h"
#include "GMS_integrands_func_ch4.h"


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
::warmup_sin(const double x)
{
    return (std::sin(x));
}

double 
gms
::fading_channel
::warmup_cos(const double x)
{
    return (std::cos(x));
}

double 
gms
::fading_channel
::warmup_tan(const double x)
{
    return (std::tan(x));
}

double 
gms 
::fading_channel
::warmup_exp(double x)
{
    return (std::exp(x));
}

double 
gms 
::fading_channel
::warmup_cyl_bess_j(const std::int32_t n,const double x)
{
    return (std::cyl_bessel_j(n,x));
}

double
gms
::fading_channel
::integrand_4_1_gauss_Q_func(const double y)
{
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 0
       const double exp_arg{0.5*(y*y)};
       const double exp_val{std::exp(-exp_arg)};
       return (C03989422804014326779399460599344*exp_val);
#else
       const double exp_arg{0.5*(y*y)};
       const double exp_val{gms::math::cephes_d::exp(-exp_arg)};
       return (C03989422804014326779399460599344*exp_val);
#endif
}



double 
gms
::fading_channel
::integrand_4_2_gauss_Q_func(const double x,const double theta)
{
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 0
       const double xsqr{x*x};
       const double sin_val{std::sin(theta)};
       const double sin_denom{2.0*(sin_val*sin_val)};
       const double exp_arg{xsqr/sin_denom};
       const double exp_val{std::exp(-exp_arg)};
       return (exp_val);
#else
       const double xsqr{x*x};
       const double sin_val{gms::math::cephes_d::sin(theta)};
       const double sin_denom{2.0*(sin_val*sin_val)};
       const double exp_arg{xsqr/sin_denom};
       const double exp_val{gms::math::cephes_d::exp(-exp_arg)};
       return (exp_val);
#endif 
}

double 
gms
::fading_channel
::integrand_4_6_gauss_Q_func(const double x1,const double y1,
                             const double rho,const double theta,
                             const std::int32_t choose_sin_or_cos)
{

#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 0
        double sin_or_cos_sqr{};
        double sin_sin_or_cos_rat{};
        const double squared_sum{(x1*x1)+(y1*y1)};
        double tsin{std::sin(theta)};
        const double sin_sqr_tht{tsin*tsin};
        const double S_hat{std::sqrt(squared_sum)};
        const double phi_s{std::atan(y1/x1)};
        const double one_m_rho{1.0f-(rho*rho)};
        const double sin2theta{std::sin(theta+theta)};
        const double one_m_sin2theta{1.0-(rho*sin2theta)};
        const double sqr_1_m_rho{std::sqrt(one_m_rho)};
        const double S_hat_half{0.5*(S_hat*S_hat)};
        const double ratio1{sqr_1_m_rho/one_m_sin2theta};
        if (choose_sin_or_cos == 1)
        {
           const double tmp{std::cos(phi_s)};
           sin_or_cos_sqr = tmp*tmp;
           sin_sin_or_cos_rat = sin_or_cos_sqr/sin_sqr_tht;
        }
        else if(choose_sin_or_cos == 2)
        {
           const double tmp{std::sin(phi_s)};
           sin_or_cos_sqr = tmp*tmp;
           sin_sin_or_cos_rat = sin_or_cos_sqr/sin_sqr_tht;
        }
        const double ratio2{one_m_sin2theta/one_m_rho};
        const double exp_arg{S_hat_half*ratio2*sin_sin_or_cos_rat};
        const double exp_val{std::exp(-exp_arg)};
        return(ratio1*exp_val);
#else 
        double sin_or_cos_sqr{};
        double sin_sin_or_cos_rat{};
        const double squared_sum{(x1*x1)+(y1*y1)}; 
        double tsin{gms::math::cephes_d::sin(theta)};
        const double sin_sqr_tht{tsin*tsin};
        const double S_hat{gms::math::cephes_d::sqrt(squared_sum)};
        const double phi_s{gms::math::cephes_d::atan(y1/x1)};
        const double one_m_rho{1.0f-(rho*rho)};
        const double sin2theta{gms::math::cephes_d::sin(theta+theta)};
        const double one_m_sin2theta{1.0-(rho*sin2theta)};
        const double sqr_1_m_rho{gms::math::cephes_d::sqrt(one_m_rho)};
        const double S_hat_half{0.5*(S_hat*S_hat)};
        const double ratio1{sqr_1_m_rho/one_m_sin2theta};
        if (choose_sin_or_cos == 1)
        {
           const double cos_tmp{gms::math::cephes_d::cos(phi_s)};
           sin_or_cos_sqr = cos_tmp*cos_tmp;
           sin_sin_or_cos_rat = sin_or_cos_sqr/sin_sqr_tht;
        }
        else if(choose_sin_or_cos == 2)
        {
           const double tmp{gms::math::cephes_d::sin(phi_s)};
           sin_or_cos_sqr = tmp*tmp;
           sin_sin_or_cos_rat = sin_or_cos_sqr/sin_sqr_tht;
        }
        const double ratio2{one_m_sin2theta/one_m_rho};
        const double exp_arg{S_hat_half*ratio2*sin_sin_or_cos_rat};
        const double exp_val{gms::math::cephes_d::exp(-exp_arg)};
        return(ratio1*exp_val);
#endif 
}

#if (INTEGRANDS_FUNC_CH4_SPLIT_MATH_IMPL_PERF_TEST) == 1

double
gms
::fading_channel
::integrand_4_6_gauss_Q_func_cephes(const double x1,const double y1,
                             const double rho,const double theta,
                             const std::int32_t choose_sin_or_cos)
{
     double sin_or_cos_sqr{};
     double sin_sin_or_cos_rat{};
     const double squared_sum{(x1*x1)+(y1*y1)}; 
     double tsin{gms::math::cephes_d::sin(theta)};
     const double sin_sqr_tht{tsin*tsin};
     const double S_hat{gms::math::cephes_d::sqrt(squared_sum)};
     const double phi_s{gms::math::cephes_d::atan(y1/x1)};
     const double one_m_rho{1.0f-(rho*rho)};
     const double sin2theta{gms::math::cephes_d::sin(theta+theta)};
     const double one_m_sin2theta{1.0-(rho*sin2theta)};
     const double sqr_1_m_rho{gms::math::cephes_d::sqrt(one_m_rho)};
     const double S_hat_half{0.5*(S_hat*S_hat)};
     const double ratio1{sqr_1_m_rho/one_m_sin2theta};
     if (choose_sin_or_cos == 1)
     {
           const double cos_tmp{gms::math::cephes_d::cos(phi_s)};
           sin_or_cos_sqr = cos_tmp*cos_tmp;
           sin_sin_or_cos_rat = sin_or_cos_sqr/sin_sqr_tht;
     }
     else if(choose_sin_or_cos == 2)
     {
           const double tmp{gms::math::cephes_d::sin(phi_s)};
           sin_or_cos_sqr = tmp*tmp;
           sin_sin_or_cos_rat = sin_or_cos_sqr/sin_sqr_tht;
     }
     const double ratio2{one_m_sin2theta/one_m_rho};
     const double exp_arg{S_hat_half*ratio2*sin_sin_or_cos_rat};
     const double exp_val{gms::math::cephes_d::exp(-exp_arg)};
     return(ratio1*exp_val);    
}

double
gms
::fading_channel
::integrand_4_6_gauss_Q_func_stdlib(const double x1,const double y1,
                                    const double rho,const double theta,
                                    const std::int32_t choose_sin_or_cos)
{
        double sin_or_cos_sqr{};
        double sin_sin_or_cos_rat{};
        const double squared_sum{(x1*x1)+(y1*y1)};
        double tsin{std::sin(theta)};
        const double sin_sqr_tht{tsin*tsin};
        const double S_hat{std::sqrt(squared_sum)};
        const double phi_s{std::atan(y1/x1)};
        const double one_m_rho{1.0f-(rho*rho)};
        const double sin2theta{std::sin(theta+theta)};
        const double one_m_sin2theta{1.0-(rho*sin2theta)};
        const double sqr_1_m_rho{std::sqrt(one_m_rho)};
        const double S_hat_half{0.5*(S_hat*S_hat)};
        const double ratio1{sqr_1_m_rho/one_m_sin2theta};
        if (choose_sin_or_cos == 1)
        {
           const double tmp{std::cos(phi_s)};
           sin_or_cos_sqr = tmp*tmp;
           sin_sin_or_cos_rat = sin_or_cos_sqr/sin_sqr_tht;
        }
        else if(choose_sin_or_cos == 2)
        {
           const double tmp{std::sin(phi_s)};
           sin_or_cos_sqr = tmp*tmp;
           sin_sin_or_cos_rat = sin_or_cos_sqr/sin_sqr_tht;
        }
        const double ratio2{one_m_sin2theta/one_m_rho};
        const double exp_arg{S_hat_half*ratio2*sin_sin_or_cos_rat};
        const double exp_val{std::exp(-exp_arg)};
        return(ratio1*exp_val);
}

#endif 

double 
gms
::fading_channel
::integrand_4_7_gauss_Q_func(const double x1,const double y1,
                             const double rho,const double theta,
                             const std::int32_t choose_x1_or_y1)
{
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 0
   
       double x1_or_y1{};
       const double rho_m_sin2tht{rho*std::sin(theta+theta)};
       const double one_m_rhosqr{1.0-(rho*rho)};
       const double tmps{std::sin(theta)};
       const double sqr_one_m_rhosqr{std::sqrt(one_m_rhosqr)};
       const double one_m_rhosintht{1.0-rho_m_sin2tht};
       const double sintht_pow2{tmps*tmps};
       const double lead_factor{sqr_one_m_rhosqr/one_m_rhosintht};
       if (choose_x1_or_y1==1)
       {
          x1_or_y1 = 0.5*(x1*x1);
       }
       else if (choose_x1_or_y1==2)
       {
          x1_or_y1 = 0.5*(y1*y1);
       }
       const double ratio{one_m_rhosintht/(one_m_rhosqr*sintht_pow2)};
       const double exp_arg{x1_or_y1*ratio};
       const double exp_val{std::exp(-exp_arg)};
       return (lead_factor*exp_val);
#else 
        double x1_or_y1{};
        const double rho_m_sin2tht{rho*gms::math::cephes_d::sin(theta+theta)};
        const double one_m_rhosqr{1.0-(rho*rho)};
        const double tmps{gms::math::cephes_d::sin(theta)};
        const double sqr_one_m_rhosqr{gms::math::cephes_d::sqrt(one_m_rhosqr)};
        const double one_m_rhosintht{1.0-rho_m_sin2tht};
        const double sintht_pow2{tmps*tmps};
        const double lead_factor{sqr_one_m_rhosqr/one_m_rhosintht};
        if (choose_x1_or_y1==1)
        {
          x1_or_y1 = 0.5*(x1*x1);
        }
        else if (choose_x1_or_y1==2)
        {
          x1_or_y1 = 0.5*(y1*y1);
        }
        const double ratio{one_m_rhosintht/(one_m_rhosqr*sintht_pow2)};
        const double exp_arg{x1_or_y1*ratio};
        const double exp_val{gms::math::cephes_d::exp(-exp_arg)};
        return (lead_factor*exp_val);
#endif 
}

double 
gms
::fading_channel
::integrand_4_8_gauss_Q_func(const double x1,const double y1,
                             const double theta,const std::int32_t choose_x1_or_y1)
{
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 0
        double x1_or_y1{};
        if(choose_x1_or_y1==1) {x1_or_y1 = x1*x1;}
        else if(choose_x1_or_y1==2) {x1_or_y1 = y1*y1;}
        const double tms{std::sin(theta)};
        const double two_sinthtsqr{tms+tms};
        const double exp_arg{-(x1_or_y1/two_sinthtsqr)};
        const double exp_val{std::exp(exp_arg)};
        return (exp_val);
#else 
        double x1_or_y1{};
        if(choose_x1_or_y1==1) {x1_or_y1 = x1*x1;}
        else if(choose_x1_or_y1==2) {x1_or_y1 = y1*y1;}
        const double tms{gms::math::cephes_d::sin(theta)};
        const double two_sinthtsqr{tms+tms};
        const double exp_arg{-(x1_or_y1/two_sinthtsqr)};
        const double exp_val{gms::math::cephes_d::exp(exp_arg)};
        return (exp_val);
#endif
}

double 
gms 
::fading_channel
::integrand_4_10_marcum_Q_func(const double x,const double s)
{
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 0
   
        const double sx{s*x};
        const double xspow2{(x*x)+(s*s)};
        const double I0_val{std::cyl_bessel_i(0,sx)};
        const double exp_arg{0.5*xspow2};
        const double exp_val{std::exp(-exp_arg)};
        return (x*exp_val*I0_val);
#else
        const double sx{s*x};
        const double xspow2{(x*x)+(s*s)};
        const double I0_val{gms::math::cephes_d::i0(sx)};
        const double exp_arg{0.5*xspow2};
        const double exp_val{gms::math::cephes_d::exp(-exp_arg)};
        return (x*exp_val*I0_val);
#endif
}

double 
gms 
::fading_channel
::integrand_4_16_marcum_Q_func(const double beta,
                               const double psi,const double theta)
{
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 0
        const double sintht{psi*std::sin(theta)};
        const double psip2{psi*psi};
        const double lead_fac_num{1.0+sintht};
        const double half_beta{0.5*(beta*beta)};
        const double lead_fac_den{1.0+(2.0*sintht)+psip2};
        const double lead_fac_rat{lead_fac_num/lead_fac_den};
        const double exp_arg{half_beta*lead_fac_den};
        const double exp_val{std::exp(-exp_arg)};
        return (lead_fac_rat*exp_val);
#else 
        const double sintht{psi*gms::math::cephes_d::sin(theta)};
        const double psip2{psi*psi};
        const double lead_fac_num{1.0+sintht};
        const double half_beta{0.5*(beta*beta)};
        const double lead_fac_den{1.0+(2.0*sintht)+psip2};
        const double lead_fac_rat{lead_fac_num/lead_fac_den};
        const double exp_arg{half_beta*lead_fac_den};
        const double exp_val{gms::math::cephes_d::exp(-exp_arg)};
        return (lead_fac_rat*exp_val);
#endif 
}

double 
gms 
::fading_channel
::integrand_4_20_marcum_Q_func_lo(const double beta,
                                  const double psi,const double theta)
{
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 0
        const double costht{psi*std::cos(theta)};
        const double psip2{psi*psi};
        const double lead_fac_num{1.0-costht};
        const double half_beta{0.5*(beta*beta)};
        const double lead_fac_den{1.0-(2.0*costht)+psip2};
        const double lead_fac_rat{lead_fac_num/lead_fac_den};
        const double exp_arg{half_beta*lead_fac_den};
        const double exp_val{std::exp(-exp_arg)};
        return (lead_fac_rat*exp_val);
#else
        double costht{psi*gms::math::cephes_d::cos(theta)};
        const double psip2{psi*psi};
        const double lead_fac_num{1.0-costht};
        const double half_beta{0.5*(beta*beta)};
        const double lead_fac_den{1.0-(2.0*costht)+psip2};
        const double lead_fac_rat{lead_fac_num/lead_fac_den};
        const double exp_arg{half_beta*lead_fac_den};
        const double exp_val{gms::math::cephes_d::exp(-exp_arg)};
        return (lead_fac_rat*exp_val);
#endif
}

double 
gms 
::fading_channel
::integrand_4_20_marcum_Q_func_hi(const double beta,
                                  const double psi,const double theta)
{
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 0
        const double costht{psi*std::cos(theta)};
        const double psip2{psi*psi};
        const double lead_fac_num{1.0+costht};
        const double half_beta{0.5*(beta*beta)};
        const double lead_fac_den{1.0+(2.0*costht)+psip2};
        const double lead_fac_rat{lead_fac_num/lead_fac_den};
        const double exp_arg{half_beta*lead_fac_den};
        const double exp_val{std::exp(-exp_arg)};
        return (lead_fac_rat*exp_val);
#else
        const double costht{psi*gms::math::cephes_d::cos(theta)};
        const double psip2{psi*psi};
        const double lead_fac_num{1.0+costht};
        const double half_beta{0.5*(beta*beta)};
        const double lead_fac_den{1.0+(2.0*costht)+psip2};
        const double lead_fac_rat{lead_fac_num/lead_fac_den};
        const double exp_arg{half_beta*lead_fac_den};
        const double exp_val{gms::math::cephes_d::exp(-exp_arg)};
        return (lead_fac_rat*exp_val);
#endif 
}

double 
gms 
::fading_channel
::integrand_4_26_marcum_Q_func(const double beta,const double psi,
                               const double theta)
{
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 0
        const double psip2{psi*psi};
        const double halfbeta{0.5*(beta*beta)};
        const double psi_m_sintht{psi*std::sin(theta)};
        const double tmp{1.0-psip2};
        const double num{tmp*tmp};
        const double denom{1.0+(2.0*psi_m_sintht)+psip2};
        const double exp_arg_1{-(halfbeta*denom)};
        const double ratio{num/denom};
        const double left_exp_val{std::exp(-exp_arg_1)};
        const double exp_arg_2{halfbeta*ratio};
        const double right_exp_val{std::exp(-exp_arg_2)};
        return (left_exp_val+right_exp_val);
#else
        const double psip2{psi*psi};
        const double halfbeta{0.5*(beta*beta)};
        const double psi_m_sintht{psi*gms::math::cephes_d::sin(theta)};
        const double tmp{1.0-psip2};
        const double num{tmp*tmp};
        const double denom{1.0+(2.0*psi_m_sintht)+psip2};
        const double exp_arg_1{halfbeta*denom};
        const double ratio{num/denom};
        const double left_exp_val{gms::math::cephes_d::exp(-exp_arg_1)};
        const double exp_arg_2{halfbeta*ratio};
        const double right_exp_val{gms::math::cephes_d::exp(-exp_arg_2)};
        return (left_exp_val+right_exp_val);
#endif 
}

double 
gms 
::fading_channel
::integrand_4_27_marcum_Q_func(const double alpha,const double psi,
                               const double theta)
{
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 0
        const double psip2{psi*psi};
        const double halfalpha{0.5*(alpha*alpha)};
        const double psi_m_sintht{psi*std::sin(theta)};
        const double tmp{1.0-psip2};
        const double num{tmp*tmp};
        const double denom{1.0+(2.0*psi_m_sintht)+psip2};
        const double exp_arg_1{halfalpha*denom};
        const double ratio{num/denom};
        const double left_exp_val{std::exp(-exp_arg_1)};
        const double exp_arg_2{halfalpha*ratio};
        const double right_exp_val{std::exp(-exp_arg_2)};
        return (left_exp_val+right_exp_val);
#else
        const double psip2{psi*psi};
        const double halfalpha{0.5*(alpha*alpha)};
        const double psi_m_sintht{psi*gms::math::cephes_d::sin(theta)};
        const double tmp{1.0-psip2};
        const double num{tmp*tmp};
        const double denom{1.0+(2.0*psi_m_sintht)+psip2};
        const double exp_arg_1{halfalpha*denom};
        const double ratio{num/denom};
        const double left_exp_val{gms::math::cephes_d::exp(-exp_arg_1)};
        const double exp_arg_2{halfalpha*ratio};
        const double right_exp_val{gms::math::cephes_d::exp(-exp_arg_2)};
        return (left_exp_val+right_exp_val);
#endif 
}

double 
gms 
::fading_channel
::integrand_4_32_marcum_Q_m_func(const double x,const double s,
                                 const double d_m,const std::int32_t i_m)
{
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 0
        const double xx{x*x};
        const double ss{s*s};
        const double sx{s*x};
        const double x_to_m{std::pow(x,d_m)};
        const double cyl_bes_val{std::cyl_bessel_i(i_m,sx)};
        const double exp_arg{0.5*(xx+ss)};
        const double exp_val{std::exp(exp_arg)};
        return (x_to_m*exp_val*cyl_bes_val);
#else 
        const double xx{x*x};
        const double ss{s*s};
        const double sx{s*x};
        const double x_to_m{gms::math::cephes_d::pow(x,d_m)};
        const double cyl_bes_val{std::cyl_bessel_i(i_m,sx)};
        const double exp_arg{0.5*(xx+ss)};
        const double exp_val{gms::math::cephes_d::exp(-exp_arg)};
        return (x_to_m*exp_val*cyl_bes_val);
#endif 
}

double
gms
::fading_channel
::integrand_4_42_marcum_Q_m_func(const double beta,const double psi,
                                 const double m,const double theta)
{
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 0
       const double psip2{psi*psi};
       const double psinpm{1.0/std::pow(psi,m-1.0)};
       const double halfbeta{0.5*(beta*beta)};
       const double tht_p_pi2{theta+C15707963267948966192313216916398};
       const double cos_arg_left{m-1.0*tht_p_pi2};
       const double cos_val_left{std::cos(cos_arg_left)};
       const double cos_arg_right{m*tht_p_pi2};
       const double cos_val_right{std::cos(cos_arg_right)};
       const double psi_m_sintht{psi*std::sin(theta)};
       const double sin_factor{1.0+(2.0*psi_m_sintht)+psip2};
       const double exp_arg{halfbeta*sin_factor};
       const double cos_num_fac{psinpm*(cos_val_left-psi*cos_val_right)};
       const double exp_val{std::exp(-exp_arg)};
       const double ratio_factor{cos_num_fac/sin_factor};
       return (ratio_factor*exp_val);
#else 
      const double psip2{psi*psi};
       const double psinpm{1.0/gms::math::cephes_d::pow(psi,m-1.0)};
       const double halfbeta{0.5*(beta*beta)};
       const double tht_p_pi2{theta+C15707963267948966192313216916398};
       const double cos_arg_left{m-1.0*tht_p_pi2};
       const double cos_val_left{gms::math::cephes_d::cos(cos_arg_left)};
       const double cos_arg_right{m*tht_p_pi2};
       const double cos_val_right{gms::math::cephes_d::cos(cos_arg_right)};
       const double psi_m_sintht{psi*gms::math::cephes_d::sin(theta)};
       const double sin_factor{1.0+(2.0*psi_m_sintht)+psip2};
       const double exp_arg{halfbeta*sin_factor};
       const double cos_num_fac{psinpm*(cos_val_left-psi*cos_val_right)};
       const double exp_val{gms::math::cephes_d::exp(-exp_arg)};
       const double ratio_factor{cos_num_fac/sin_factor};
       return (ratio_factor*exp_val);
#endif 
}

double
gms
::fading_channel
::integrand_4_45_marcum_Q_m_func(const double theta,const double beta,const double m)
{
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 0
        const double pow_arg{1.0+2.0*m};
        const double betap2{beta*beta};
        const double tmps{std::sin(theta)};
        const double sinp2{tmps*tmps};
        const double costht{std::cos(theta)};
        const double denom{2.0*sinp2};
        const double sinpm{std::pow(tmps,pow_arg)};
        const double exp_arg{betap2/denom};
        const double cot_ratio{costht/sinpm};
        const double exp_val{std::exp(-exp_arg)};
        return (cot_ratio*exp_val);
#else 
        const double pow_arg{1.0+2.0*m};// FMA(1,2,m), 1+m+m
        const double betap2{beta*beta};
        const double tmps{gms::math::cephes_d::sin(theta)};
        const double sinp2{tmps*tmps};
        const double costht{gms::math::cephes_d::cos(theta)};
        const double denom{2.0*sinp2};
        const double sinpm{gms::math::cephes_d::pow(tmps,pow_arg)};
        const double exp_arg{betap2/denom};
        const double cot_ratio{costht/sinpm};
        const double exp_val{gms::math::cephes_d::exp(-exp_arg)};
        return (cot_ratio*exp_val);
#endif 
}

double
gms
::fading_channel
::integrand_4_66_pawula_func(const double psi,const double delphi,const double t,
                             const double A,const double sigmasqr) 
{
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 0
        const double phase_diff{delphi-psi};
        const double A_factor{(A*A)/(sigmasqr+sigmasqr)};
        const double cost{std::cos(t)};
        const double cosphdif{std::cos(phase_diff)};
        const double cosfdiff_fac{1.0-cosphdif*cost};
        const double exp_arg{A_factor*cosfdiff_fac};
        const double inv_cosdiff{1.0/cosfdiff_fac};
        const double exp_val{-exp_arg};
        return (inv_cosdiff*exp_val);
#else
        const double phase_diff{delphi-psi};
        const double A_factor{(A*A)/(sigmasqr+sigmasqr)};
        const double cost{gms::math::cephes_d::cos(t)};
        const double cosphdif{gms::math::cephes_d::cos(phase_diff)};
        const double cosfdiff_fac{1.0-cosphdif*cost};
        const double exp_arg{A_factor*cosfdiff_fac};
        const double inv_cosdiff{1.0/cosfdiff_fac};
        const double exp_val{-exp_arg};
        return (inv_cosdiff*exp_val);
#endif 
}

double
gms
::fading_channel
::integrand_4_67_pawula_func(const double psi,const double t,
                             const double A,const double sigmasqr)
{
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 0
        const double A_factor{(A*A)/(sigmasqr+sigmasqr)};
        const double cost{std::cos(t)};
        const double cospsi{std::cos(psi)};
        const double cos_term{1.0-cospsi*cost};
        const double exp_arg{A_factor*cos_term};
        const double lead_factor{1.0/cos_term};
        const double exp_value{std::exp(-exp_arg)};
        return (lead_factor*exp_value);
#else
        const double A_factor{(A*A)/(sigmasqr+sigmasqr)};
        const double cost{gms::math::cephes_d::cos(t)};
        const double cospsi{gms::math::cephes_d::cos(psi)};
        const double cos_term{1.0-cospsi*cost};
        const double exp_arg{A_factor*cos_term};
        const double lead_factor{1.0/cos_term};
        const double exp_value{gms::math::cephes_d::exp(-exp_arg)};
        return (lead_factor*exp_value);
#endif 
}

