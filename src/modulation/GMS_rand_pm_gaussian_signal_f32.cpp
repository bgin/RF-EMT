
#include "GMS_rand_pm_gaussian_signal_f32.h"
#include "GMS_sse_memset.h"
#include "GMS_compare_fp_safe.h"

gms::radiolocation
::rand_pm_gaussian_signal_f32_t
::rand_pm_gaussian_signal_f32_t() noexcept(true) 
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
::rand_pm_gaussian_signal_f32_t
::rand_pm_gaussian_signal_f32_t(const std::size_t nsamples,
                                const std::size_t num_gaussians,
                                const float tstart,
                                const float tend,
                                const float min_sigma,
                                const float max_sigma,
                                const float max_amplitude)
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
::rand_pm_gaussian_signal_f32_t
::rand_pm_gaussian_signal_f32_t(const rand_pm_gaussian_signal_f32_t &other) noexcept(false) 
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
::rand_pm_gaussian_signal_f32_t
::rand_pm_gaussian_signal_f32_t(rand_pm_gaussian_signal_f32_t &&other) noexcept(true) 
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
::rand_pm_gaussian_signal_f32_t
::~rand_pm_gaussian_signal_f32_t() noexcept(true) 
{

}

gms::radiolocation
::rand_pm_gaussian_signal_f32_t &
gms::radiolocation
::rand_pm_gaussian_signal_f32_t
::operator=(const rand_pm_gaussian_signal_f32_t &other) noexcept(false) 
{
    if(__builtin_except(this==&other)) { return (*this);}
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
::rand_pm_gaussian_signal_f32_t &
gms::radiolocation
::rand_pm_gaussian_signal_f32_t
::operator=(rand_pm_gaussian_signal_f32_t &&other) noexcept(true) 
{
    if(__builtin_except(this==&other)) { return (*this);}
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
::rand_pm_gaussian_signal_f32_t
::operator==(const rand_pm_gaussian_signal_f32_t &other)->std::vector<bool> 
{
    if(__builtin_expect(this->m_nsamples!=other.m_signal_samples,0)) {return (std::vector<bool>());}
    std::vector<bool> eq_comp_holder;
    float * __restrict__ p_other = &other.m_signal_samples.m_data[0];
    constexpr float eps = std::numeric_limits<float>::epsilon();
    float x;
    float y;
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
::rand_pm_gaussian_signal_f32_t
::operator!=(const rand_pm_gaussian_signal_f32_t &other)->std::vector<bool>                                         
{   // Just call the operator== and inspect the boolean values.
    return (this->operator==(other));
}

auto 
gms::radiolocation
::rand_pm_gaussian_signal_f32_t
::operator>=(const rand_pm_gaussian_signal_f32_t &other)->std::vector<bool>                                               
{
    std::vector<bool> gt_comp_holder;
    double * __restrict__ p_other = &other.m_sig_samples.m_data[0];
    constexpr double eps = std::numeric_limits<double>::epsilon();
    double x;
    double y;
    for(std::size_t i = 0ull; i<this->m_nsamples; ++i) 
    {
        x = this->m_sigNAL_samples.m_data[i];
        y = p_other[i];
        const bool xy_gt = gms::definitelyGreaterThan(x,y,eps);
        gt_comp_holder.operator[](i) = xy_gt;
    }
    return (gt_comp_holder);
}

auto 
gms::radiolocation
::rand_pm_gaussian_signal_f32_t
::operator<=(const rand_pm_gaussian_signal_f32_t &other)->std::vector<bool>                                                             
{
    std::vector<bool> lt_comp_holder;
    double * __restrict__ p_other = &other.m_sig_samples.m_data[0];
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

template<std::int32_t N>
std::int32_t 
gms::radiolocation
::rand_pm_gaussian_signal_f32_t
::generate()
{
    std::int32_t ret_status = -9999;
    if(__builtin_expect(this->m_num_gaussians>N,0)) { ret_status = -1; return (ret_status);}
    struct alignas(16) gaussian_components_t 
    {
        float amplitude;
        float mu;
        float two_sigma_sq;
    };
    alignas(64) gaussian_components_t gc_store[N];
    thread_local std::uniform_real_distribution<float> dist_mu;
    thread_local std::uint64_t seed_mu{};
    thread_local std::mt19937  dist_mu_gen;
    thread_local std::uniform_real_distribution<float> dist_sigma;
    thread_local std::uint64_t seed_sigma{};
    thread_local std::mt19937 dist_sigma_gen;
    thread_local std::uniform_real_distribution<float> dist_amp;
    thread_local std::uint64_t seed_amp{};
    thread_local std::mt19937 dist_amp_gen;
    const float  padding = 3.0f*this->m_max_sigma;
    seed_mu = __rdtsc();
    dist_mu = std::uniform_real_distribution<float>(this->m_tstart-padding,this->m_tend+padding);
    dist_mu_gen = std::mt19937(seed_mu);
    seed_sigma = __rdtsc();
    dist_sigma = std::uniform_real_distribution<float>(this->m_min_sigma,this->m_max_sigma);
    dist_sigma_gen = std::mt19937(seed_sigma);
    seed_amp   = __rdtsc();
    dist_amp   = std::uniform_real_distribution<float>(-this->m_max_amplitude,this->m_max_amplitude);
    dist_amp_gen   = std::mt19937(seed_amp);
    for(std::size_t i = 0; i<this->m_num_gaussians; ++i) 
    {
        gc_store[i].amplitude    = dist_amp.operator()(dist_amp_gen);
        gc_store[i].mu           = dist_mu.operator()(dist_mu_gen);
        float sigma_tmp          = dist_sigma.operator()(dist_sigma_gen);
        gc_store[i].two_sigma_sq = 2.0f*sigma_tmp*sigma_tmp; 

    }
    const float dt = (this->m_tend-this->m_tstart)/static_cast<float>(this->m_nsamples-1ull);
    for(std::size_t j = 0ull; j<this->m_nsamples; ++j)   
    {
        const float t_j = static_cast<float>(j);
        const float t   = this->m_tstart+t_j*dt;
        float phase_accumulator = 0.0f;
        for(std::size_t k = 0; k<this->m_num_gaussians; ++k) 
        {   
            const float tmu        = gc_store[k].mu;
            float diff             = t-tmu;
            const float tamplitude = gc_store[k].amplitude;
            const float tsigma     = gc_store[k].two_sigma_sq;
#if (RAND_PM_GAUSSIAN_SIGNAL_F32_USE_CEPHES) == 1
            const float exp_val    = gms::math::ceph_expf(-(diff*diff)/tsigma);
#else 
            const float exp_val    = std::exp(-(diff*diff)/tsigma);
#endif 
            phase_accumulator += (tamplitude*exp_val);
        }
        this->m_signal_samples.m_data[j] = phase_accumulator;
    }
    ret_status = 0;
    return (ret_status);
}