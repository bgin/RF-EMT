
#include <cmath>
#include <algorithm>
#include <numeric>
#include "GMS_normality_test_v2.h"

__ATTR_ALWAYS_INLINE__
static inline double 
rational_approximation(const double t)
{
    const double c[] = {2.515517, 0.802853, 0.010328};
    const double d[] = {1.432788, 0.189269, 0.001308};
    return t - ((c[2]*t + c[1])*t + c[0]) /
               (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}

__ATTR_ALWAYS_INLINE__
static inline double
normal_cdf_inverse(const double p)
{
    if(__builtin_expect(p<=0.0 || p>=1.0,0)) { return (p);}
    if (p < 0.5) 
    {
        return -rational_approximation(std::sqrt(-2.0*std::log(p)));
    } 
    else 
    {
        return rational_approximation(std::sqrt(-2.0*std::log(1.0-p)));
    }
}

double 
gms::math
::calculate_W_statistics(const darray_r8_t &x)
{
    const std::size_t n{x.mnx};
    const double d_n{1.0/(static_cast<double>(n)+0.25)};
    darray_r8_t m_vec(std::move(x));
    double m_val{0.0};
    for(std::size_t __i{0ull}; __i != n; ++__i) 
    {
        const double t__i{static_cast<double>(__i+1ull)};
        m_vec.m_data[__i] = normal_cdf_inverse((t__i-0.375)*d_n);
        const double t__d{m_vec.m_data[__i]};
        m_val += t__d*t__d;
    }

    const double pow_m_val{1.0/std::sqrt(m_val)};
    const double u{1.0/static_cast<double>(n)};
    const double u2{u*u};
    const double u3{u2*u};
    const double u4{u3*u};
    const double u5{u4*u};
    darray_r8_t coeff(n);

    coeff.m_data[n-1ull] = -2.706056*u5+4.434685*u4-2.2071190*u3-
                            0.147981*u2+0.221157*u+m_vec.m_data[n-1ull]*pow_m_val;   /*pow (m_val, -0.5);*/
    coeff.m_data[0ull]   = -1.0*coeff.m_data[n-1ull];
    coeff.m_data[n-2ull] = -3.582633*u5+ 5.682633*u4-1.752461*u3-
                            0.293762*u2+0.042981*u+m_vec.m_data[n-2ull]*pow_m_val;
    coeff.m_data[1ull]   = -1.0*coeff.m_data[n-2ull];

    const double s_mv1{m_vec.m_data[n-1ull]*m_vec.m_data[n-1ull]};
    const double s_mv2{m_vec.m_data[n-2ull]*m_vec.m_data[n-2ull]};
    const double s_dc1{coeff.m_data[n-1ull]*coeff.m_data[n-1ull]};
    const double s_dc2{coeff.m_data[n-2ull]*coeff.m_data[n-2ull]};
    const double mv_eps{m_val-(2.0*s_mv1)-(2.0*s_mv2)};
    const double dc_eps{1.0-(2.0*s_dc1)-(2.0*s_dc2)};
    const double eps{mv_eps/dc_eps};
    const double sqrt_eps{1.0/std::sqrt(eps)};
    for(std::size_t __i{1ull}; __i != n-2ull; ++__i) { coeff.m_data[__i] = m_vec.m_data[__i]*sqrt_eps; }

    double denom{0.0};
    double num{0.0};
    const double accum{std::accumulate(x.begin(),x.end(),0.0)};
    const double mean{accum/static_cast<double>(n)};
    darray_r8_t cpy_x(std::move(x));   
    std::sort(cpy_x.begin(),cpy_x.end());
    for(std::size_t __i{0ull}; __i != n; ++__i) 
    {
        const double t__cpy_x{cpy_x.m_data[__i]};
        num   += coeff.m_data[__i]*t__cpy_x;
        denom += (t__cpy_x*t__cpy_x)-(mean*mean); 
    }
    
    double W{(num*num)/denom};
    return (W);
    
}