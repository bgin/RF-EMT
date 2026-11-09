
#include <iomanip>
#include "GMS_rand_pm_gaussian_signal_f64.h"
#include "GMS_compare_fp_safe.h"

gms::radiolocation
::rand_pm_gaussian_signal_f64_t
::rand_pm_gaussian_signal_f64_t() noexcept(true) 
:
m_nsamples{},
m_num_gaussians{},
m_tstart{},
m_tend{},
m_min_sigma{},
m_max_sigma{},
m_max_amplitude{},
m_signal_samples()
{

}

gms::radiolocation
::rand_pm_gaussian_signal_f64_t
::rand_pm_gaussian_signal_f64_t(const std::size_t nsamples,
                                const std::size_t num_gaussians,
                                const double tstart,
                                const double tend,
                                const double min_sigma,
                                const double max_sigma,
                                const double max_amplitude)
:
m_nsamples(nsamples),
m_num_gaussians((num_gaussians>1000ull)?1000ull:num_gaussians),
m_tstart((tstart<=0.0f)?0.0:tstart),
m_tend((tend<=m_tstart)?m_tstart+2.0f:tend),
m_min_sigma(min_sigma),
m_max_sigma(max_sigma),
m_max_amplitude(max_amplitude),
m_signal_samples(m_nsamples)
{

}

gms::radiolocation
::rand_pm_gaussian_signal_f64_t
::rand_pm_gaussian_signal_f64_t(const rand_pm_gaussian_signal_f64_t &other) noexcept(false) 
:
m_nsamples(other.m_nsamples),
m_num_gaussians(other.m_num_gaussians),
m_tstart(other.m_tstart),
m_tend(other.m_tend),
m_min_sigma(other.m_min_sigma),
m_max_sigma(other.m_max_sigma),
m_max_amplitude(other.m_max_amplitude),
m_signal_samples(other.m_signal_samples)
{

}

gms::radiolocation
::rand_pm_gaussian_signal_f64_t
::rand_pm_gaussian_signal_f64_t(rand_pm_gaussian_signal_f64_t &&other) noexcept(true) 
:
m_nsamples(other.m_nsamples),
m_num_gaussians(other.m_num_gaussians),
m_tstart(other.m_tstart),
m_tend(other.m_tend),
m_min_sigma(other.m_min_sigma),
m_max_sigma(other.m_max_sigma),
m_max_amplitude(other.m_max_amplitude),
m_signal_samples(std::move(other.m_signal_samples))
{

}

gms::radiolocation
::rand_pm_gaussian_signal_f64_t
::~rand_pm_gaussian_signal_f64_t() noexcept(true) 
{

}

gms::radiolocation
::rand_pm_gaussian_signal_f64_t &
gms::radiolocation
::rand_pm_gaussian_signal_f64_t
::operator=(const rand_pm_gaussian_signal_f64_t &other) noexcept(false) 
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_nsamples      = other.m_nsamples;
    this->m_num_gaussians = other.m_num_gaussians;
    this->m_tstart        = other.m_tstart;
    this->m_tend          = other.m_tend;
    this->m_min_sigma     = other.m_min_sigma;
    this->m_max_sigma     = other.m_max_sigma;
    this->m_max_amplitude = other.m_max_amplitude;
    this->m_signal_samples.operator=(other.m_signal_samples);
    return (*this);
}

gms::radiolocation
::rand_pm_gaussian_signal_f64_t &
gms::radiolocation
::rand_pm_gaussian_signal_f64_t
::operator=(rand_pm_gaussian_signal_f64_t &&other) noexcept(true) 
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_nsamples      = other.m_nsamples;
    this->m_num_gaussians = other.m_num_gaussians;
    this->m_tstart        = other.m_tstart;
    this->m_tend          = other.m_tend;
    this->m_min_sigma     = other.m_min_sigma;
    this->m_max_sigma     = other.m_max_sigma;
    this->m_max_amplitude = other.m_max_amplitude;
    this->m_signal_samples.operator=(std::move(other.m_signal_samples));
    return (*this);
}

auto 
gms::radiolocation 
::rand_pm_gaussian_signal_f64_t
::operator==(const rand_pm_gaussian_signal_f64_t &other)->std::vector<bool> 
{
    if(__builtin_expect(this->m_nsamples!=other.m_nsamples,0)) {return (std::vector<bool>());}
    std::vector<bool> eq_comp_holder;
    double * __restrict__ p_other = &other.m_signal_samples.m_data[0];
    constexpr double eps = std::numeric_limits<double>::epsilon();
    double x;
    double y;
    for(std::size_t i = 0ull; i<this->m_nsamples; ++i) 
    {
        x = this->m_signal_samples.m_data[i];
        y = p_other[i];
        const bool xy_eq = gms::approximatelyEqual(x,y,eps);
        eq_comp_holder.operator[](i) = xy_eq;
    }
    return (eq_comp_holder);
}

auto 
gms::radiolocation
::rand_pm_gaussian_signal_f64_t
::operator!=(const rand_pm_gaussian_signal_f64_t &other)->std::vector<bool>                                         
{   // Just call the operator== and inspect the boolean values.
    return (this->operator==(other));
}

auto 
gms::radiolocation
::rand_pm_gaussian_signal_f64_t
::operator>=(const rand_pm_gaussian_signal_f64_t &other)->std::vector<bool>                                               
{
    if(__builtin_expect(this->m_nsamples!=other.m_nsamples,0)) {return (std::vector<bool>());}
    std::vector<bool> gt_comp_holder;
    double * __restrict__ p_other = &other.m_signal_samples.m_data[0];
    constexpr double eps = std::numeric_limits<double>::epsilon();
    double x;
    double y;
    for(std::size_t i = 0ull; i<this->m_nsamples; ++i) 
    {
        x = this->m_signal_samples.m_data[i];
        y = p_other[i];
        const bool xy_gt = gms::definitelyGreaterThan(x,y,eps);
        gt_comp_holder.operator[](i) = xy_gt;
    }
    return (gt_comp_holder);
}

auto 
gms::radiolocation
::rand_pm_gaussian_signal_f64_t
::operator<=(const rand_pm_gaussian_signal_f64_t &other)->std::vector<bool>
{
    if(__builtin_expect(this->m_nsamples!=other.m_nsamples,0)) {return (std::vector<bool>());}
    std::vector<bool> lt_comp_holder;
    double * __restrict__ p_other = &other.m_signal_samples.m_data[0];
    constexpr double eps = std::numeric_limits<double>::epsilon();
    double x;
    double y;
    for(std::size_t i = 0ull; i<this->m_nsamples; ++i) 
    {
        x = this->m_signal_samples.m_data[i];
        y = p_other[i];
        const bool xy_lt = gms::definitelyLessThan(x,y,eps);
        lt_comp_holder.operator[](i) = xy_lt;
    }
    return (lt_comp_holder);
}

auto 
gms::radiolocation
::operator<<(std::ostream &os,
             const rand_pm_gaussian_signal_f64_t &rhs)->std::ostream &
{
    std::cout << typeid(rhs).name() << "Begin: object state dump." << std::endl;
    std::cout << "m_nsamples     : "      << rhs.m_nsamples << std::endl;
    std::cout << "m_num_gaussians: "      << rhs.m_num_gaussians       << std::endl;
    std::cout << "m_tstart       : "      << std::fixed << std::setprecision(17) << rhs.m_tstart << std::endl;
    std::cout << "m_tend         : "      << std::fixed << std::setprecision(17) << rhs.m_tend << std::endl;
    std::cout << "m_min_sigma    : "      << std::fixed << std::setprecision(17) << rhs.m_min_sigma << std::endl;
    std::cout << "m_max_sigma    : "      << std::fixed << std::setprecision(17) << rhs.m_max_sigma << std::endl;
    std::cout << "m_max_amplitude: "      << std::fixed << std::setprecision(17) << rhs.m_max_amplitude << std::endl;
    std::cout << "Randomly PM Gaussian Signal Samples:" << std::endl;
    for(std::size_t i{0ull}; i != rhs.m_nsamples; ++i)
    {
         os << std::fixed << std::setprecision(17) << 
                             rhs.m_signal_samples.m_data[i] << std::endl;
    }
    std::cout << typeid(rhs).name() << "End: object state dump." << std::endl;
    return (os);
}
