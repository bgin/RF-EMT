#include <fstream>
#include <iomanip>
#include <random>
#include "GMS_pm_bb_cmplx_sine_signal.h"
#include "GMS_sse_memset.h"
#include "GMS_indices.h"

gms::radiolocation
::pm_bb_cmplx_sin_signal_t
::pm_bb_cmplx_sin_signal_t(const std::size_t nsamples,
                           const std::uint32_t nK,
                           const float w0,
                           const float ph,
                           const float n,
                           const float A,
                           const float P)
:
m_nsamples{nsamples},
m_nK{nK},
m_w0{w0},
m_ph{ph},
m_n{n},
m_A{A},
m_P{P},
m_sig_samples{darray_c4_t(m_nsamples)}
{

}

gms::radiolocation
::pm_bb_cmplx_sin_signal_t
::pm_bb_cmplx_sin_signal_t(pm_bb_cmplx_sin_signal_t &&other) noexcept(true)
:
m_nsamples{std::move(other.m_nsamples)},
m_nK{std::move(other.m_nK)},
m_n{std::move(other.m_n)},
m_w0{std::move(other.m_w0)},
m_ph{std::move(other.m_ph)},
m_A{std::move(other.m_A)},
m_P{std::move(other.m_P)},
m_sig_samples{std::move(other.m_sig_samples)}
{
    
}   

gms::radiolocation 
::pm_bb_cmplx_sin_signal_t
::~pm_bb_cmplx_sin_signal_t()
{

}

gms::radiolocation
::pm_bb_cmplx_sin_signal_t &
gms::radiolocation
::pm_bb_cmplx_sin_signal_t
::operator=(pm_bb_cmplx_sin_signal_t &&other) noexcept(true)
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_nsamples  = std::move(other.m_nsamples);
    this->m_nK        = std::move(other.m_nK);
    this->m_w0        = std::move(other.m_w0);
    this->m_ph        = std::move(other.m_ph);
    this->m_n         = std::move(other.m_n);
    this->m_A         = std::move(other.m_A);
    this->m_P         = std::move(other.m_P);
    this->m_sig_samples.operator=(std::move(other.m_sig_samples));
    return (*this);
}

void 
gms::radiolocation
::pm_bb_cmplx_sin_signal_t
::init_storage(const float filler)
{
#if (INIT_BY_STD_FILL) == 0
     using namespace gms::common;
	 sse_memset_unroll8x_ps(reinterpret_cast<float*>(this->m_sig_samples.m_data),filler,2ull*this->m_nsamples);
#else 
     std::fill(this->m_sig_samples.m_data,this->m_sig_samples.m_data+this->m_nsamples,filler);
#endif
}

void 
gms::radiolocation
::pm_bb_cmplx_sin_signal_t
::create_signal_plot(const std::uint32_t n_samp,
                     const float * __restrict sig_arg,
                     const float * __restrict sig_val,
                     const std::string &header,
                     const std::string &title,
                     const bool is_sig_arg_present)
{
    std::string plot_fname;
    std::string sig_fname;
    std::ofstream plot_unit;
    std::ofstream sig_unit;
    sig_fname = header+"_plot.txt";
    sig_unit.open(sig_fname.c_str());
    if(is_sig_arg_present==true)
    {
        for(std::size_t __i{0ull}; __i != n_samp; ++__i)
        {
            sig_unit << " " << sig_arg[__i] << " "
                            << sig_val[__i] << "\n";
        }
    }
    else
    {
        for(std::size_t __i{0ull}; __i != n_samp; ++__i)
        {
            sig_unit << " " << sig_arg[__i] << "\n";
         
        }
    }
    sig_unit.close();
    std::cout << "Created signal data file \"" << sig_fname << "\".\n";
    plot_fname = header+"plot_commands.txt";
    plot_unit.open(plot_fname.c_str());
    plot_unit << "#" << plot_fname << "\n";
    plot_unit << "#\n";
    plot_unit << "# Usage:\n";
    plot_unit << "# gnuplot < " << plot_fname << "\n";
    plot_unit << "#\n";
    plot_unit << "set term png\n";
    plot_unit << "set output \"" << header << ".png\"\n";
    plot_unit << "set xlabel 't'\n";
    plot_unit << "set ylabel 'y(t)'\n";
    plot_unit << "set title '" << title << "'\n";
    plot_unit << "set grid\n";
    plot_unit << "set style data lines\n";
    if(is_sig_arg_present==true)
    {
            plot_unit << "plot \"" << sig_fname << "\" using 1:2 lw 1 linecolor rgb \"red\"\n";
    }
    else
    {
            plot_unit << "plot \"" << sig_fname << "\" lw 1 linecolor rgb \"red\"\n";
    }
    plot_unit << "quit\n";
    plot_unit.close();
    std::cout << " Created signal data file \"" << plot_fname << "\"\n";
}

std::int32_t 
gms::radiolocation
::pm_bb_cmplx_sin_signal_t
::create_signal_user_data(const float * __restrict__ sym_in, // size of m_nsamples*m_nK
                          const std::uint32_t n_T,
                          const std::uint32_t n_K)
{
      using namespace gms::math;
      if(__builtin_expect(static_cast<std::uint32_t>(this->m_nsamples)!=n_T,0) || 
         __builtin_expect(this->m_nK!=n_K,0)) { return (-1);}

       const float T{static_cast<float>(this->m_nsamples)};
       const float invTc{1.0f/T};
       constexpr float C141421356237309504880168872421{1.41421356237309504880168872421f};
       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
       const float invT{this->m_P/T};
       std::complex<float> sum;
       for(std::uint32_t __t{0}; __t != n_T; ++__t) 
       {
            const float t{static_cast<float>(__t)};
            const float carrier_arg{this->m_ph+(C6283185307179586476925286766559*this->m_w0*t*invTc)};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1
            const float cos_carrier{ceph_cosf(carrier_arg)};
            const float sin_carrier{ceph_sinf(carrier_arg)};
#else 
            const float cos_carrier{std::cos(carrier_arg)};
            const float sin_carrier{std::sin(carrier_arg)};
#endif
            sum = (0.0f,0.0f);
            for(std::uint32_t __k{0}; __k != n_K; ++__k) 
            {
                 const float k{static_cast<float>(__k)};
                 const float arg{t-k*T};
                 const float sym{sym_in[Ix2D(__t,n_K,__k)]};
                 const float pm{sin_sample(arg,invT)*sym};
                 //const float re_im{1.0f-(2.0f*sym)};
                 //const std::complex<float> csym{re_im,re_im};
                 //sum += sin_sample(arg,invT)*(C141421356237309504880168872421*csym);
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1                
                const float  I{ceph_cosf(pm)*cos_carrier};
                const float  Q{ceph_sinf(pm)*sin_carrier};
#else 
                const float  I{std::cos(pm)*cos_carrier};
                const float  Q{std::sin(pm)*sin_carrier};
#endif
                const std::complex<float> sample{I,-Q};
                sum += sample;
            }
            this->m_sig_samples.m_data[__t] = sum;
        }
        return (0);
}

std::int32_t 
gms::radiolocation
::pm_bb_cmplx_sin_signal_t
::create_signal_user_data_u4x(const float * __restrict__ sym_in, // size of m_nsamples*m_nK
                              const std::uint32_t n_T,
                              const std::uint32_t n_K)
{
      using namespace gms::math;
      if(__builtin_expect(static_cast<std::uint32_t>(this->m_nsamples)!=n_T,0) || 
         __builtin_expect(this->m_nK!=n_K,0)) { return (-1);}

        const float T{static_cast<float>(this->m_nsamples)};
        const float invTc{1.0f/T};
        constexpr float C141421356237309504880168872421{1.41421356237309504880168872421f};
        constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
        const float invT{this->m_P/T};
        std::complex<float> sum0;
        std::complex<float> sum1;
        std::complex<float> sum2;
        std::complex<float> sum3;
        float cos_carrier_0;
        float sin_carrier_0;
        float cos_carrier_1;
        float sin_carrier_1;
        float cos_carrier_2;
        float sin_carrier_2;
        float cos_carrier_3;
        float sin_carrier_3;
        float carrier_arg_0;
        float carrier_arg_1;
        float carrier_arg_2;
        float carrier_arg_3;
        float t__i_0{-4.0f};
        float t__i_1{1.0f};
        float t__i_2{2.0f};
        float t__i_3{3.0f};
        std::uint32_t __i;
        std::uint32_t __j;
        for(__i = 0; __i != ROUND_TO_FOUR(n_T,4); __i += 4u)
        {
            t__i_0 += 4.0f;
            const float carrier_arg_0{this->m_ph+(C6283185307179586476925286766559*this->m_w0*t__i_0*invTc)};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1
            cos_carrier_0 = ceph_cosf(carrier_arg_0);
            sin_carrier_0 = ceph_sinf(carrier_arg_0);
#else 
           cos_carrier_0 = std::cos(carrier_arg_0);
           sin_carrier_0 = std::sin(carrier_arg_0);
#endif
            sum0    = (0.0f,0.0f);
            for(std::uint32_t __k{0}; __k != n_K; ++__k) 
            {
                const float k{static_cast<float>(__k)};
                const float arg{t__i_0-k*T};
                const float sym{sym_in[Ix2D(__i+0,n_K,__k)]};
                //const float re_im{1.0f-(2.0f*sym)};
                //const std::complex<float> csym{re_im,re_im};
                //sum0 += sin_sample(arg,invT)*(C141421356237309504880168872421*csym);
                const float pm{sin_sample(arg,invT)*sym};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1                
                const float  I{ceph_cosf(pm)*cos_carrier_0};
                const float  Q{ceph_sinf(pm)*sin_carrier_0};
#else 
                const float  I{std::cos(pm)*cos_carrier_0};
                const float  Q{std::sin(pm)*sin_carrier_0};
#endif
                const std::complex<float> sample{I,-Q};
                sum0 += sample;
            }
            this->m_sig_samples.m_data[__i+0] = sum0;
            t__i_1 += 4.0f;
            const float carrier_arg_1{this->m_ph+(C6283185307179586476925286766559*this->m_w0*t__i_1*invTc)};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1
            cos_carrier_1 = ceph_cosf(carrier_arg_1);
            sin_carrier_1 = ceph_sinf(carrier_arg_1);
#else 
           cos_carrier_1 = std::cos(carrier_arg_1);
           sin_carrier_1 = std::sin(carrier_arg_1);
#endif
            sum1    = (0.0f,0.0f);
            for(std::uint32_t __k{0}; __k != n_K; ++__k) 
            {
                const float k{static_cast<float>(__k)};
                const float arg{t__i_1-k*T};
                const float sym{sym_in[Ix2D(__i+1,n_K,__k)]};
                //const float re_im{1.0f-(2.0f*sym)};
               // const std::complex<float> csym{re_im,re_im};
                //sum1 += sin_sample(arg,invT)*(C141421356237309504880168872421*csym);
                const float pm{sin_sample(arg,invT)*sym};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1                
                const float  I{ceph_cosf(pm)*cos_carrier_1};
                const float  Q{ceph_sinf(pm)*sin_carrier_1};
#else 
                const float  I{std::cos(pm)*cos_carrier_1};
                const float  Q{std::sin(pm)*sin_carrier_1};
#endif
                const std::complex<float> sample{I,-Q};
                sum1 += sample;
            }
            this->m_sig_samples.m_data[__i+1] = sum1;
            t__i_2  += 4.0f;
            const float carrier_arg_2{this->m_ph+(C6283185307179586476925286766559*this->m_w0*t__i_2*invTc)};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1
            cos_carrier_2 = ceph_cosf(carrier_arg_2);
            sin_carrier_2 = ceph_sinf(carrier_arg_2);
#else 
           cos_carrier_2 = std::cos(carrier_arg_2);
           sin_carrier_2 = std::sin(carrier_arg_2);
#endif
            sum2     = (0.0f,0.0f);
            for(std::uint32_t __k{0}; __k != n_K; ++__k) 
            {
                const float k{static_cast<float>(__k)};
                const float arg{t__i_2-k*T};
                const float sym{sym_in[Ix2D(__i+2,n_K,__k)]};
                //const float re_im{1.0f-(2.0f*sym)};
                //const std::complex<float> csym{re_im,re_im};
                //sum2 += sin_sample(arg,invT)*(C141421356237309504880168872421*csym);
                 const float pm{sin_sample(arg,invT)*sym};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1                
                 const float  I{ceph_cosf(pm)*cos_carrier_2};
                 const float  Q{ceph_sinf(pm)*sin_carrier_2};
#else 
                 const float  I{std::cos(pm)*cos_carrier_2};
                 const float  Q{std::sin(pm)*sin_carrier_2};
#endif
                 const std::complex<float> sample{I,-Q};
                 sum2 += sample;
            }
            this->m_sig_samples.m_data[__i+2] = sum2;
            t__i_3  += 4.0f;
            const float carrier_arg_3{this->m_ph+(C6283185307179586476925286766559*this->m_w0*t__i_3*invTc)};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1
            cos_carrier_3 = ceph_cosf(carrier_arg_3);
            sin_carrier_3 = ceph_sinf(carrier_arg_3);
#else 
           cos_carrier_3 = std::cos(carrier_arg_3);
           sin_carrier_3 = std::sin(carrier_arg_3);
#endif
            sum3     = (0.0f,0.0f);
            for(std::uint32_t __k{0}; __k != n_K; ++__k) 
            {
                const float k{static_cast<float>(__k)};
                const float arg{t__i_3-k*T};
                const float sym{sym_in[Ix2D(__i+3,n_K,__k)]};
                //const float re_im{1.0f-(2.0f*sym)};
                //const std::complex<float> csym{re_im,re_im};
                //sum3 += sin_sample(arg,invT)*(C141421356237309504880168872421*csym);
                 const float pm{sin_sample(arg,invT)*sym};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1                
                 const float  I{ceph_cosf(pm)*cos_carrier_3};
                 const float  Q{ceph_sinf(pm)*sin_carrier_3};
#else 
                 const float  I{std::cos(pm)*cos_carrier_3};
                 const float  Q{std::sin(pm)*sin_carrier_3};
#endif
                 const std::complex<float> sample{I,-Q};
                 sum3 += sample;
            }
            this->m_sig_samples.m_data[__i+3] = sum3;
            }
            for(__j = __i; __j != n_T; ++__j) 
            {
                const float t{static_cast<float>(__j)};
                const float carrier_arg_0{this->m_ph+(C6283185307179586476925286766559*this->m_w0*t*invTc)};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1
                cos_carrier_0 = ceph_cosf(carrier_arg_0);
                sin_carrier_0 = ceph_sinf(carrier_arg_0);
#else 
                cos_carrier_0 = std::cos(carrier_arg_0);
                sin_carrier_0 = std::sin(carrier_arg_0);
#endif                
                sum0 = (0.0f,0.0f);
                for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                {
                    const float k{static_cast<float>(__k)};
                    const float arg{t-k*T};
                    const float sym{sym_in[Ix2D(__j,n_K,__k)]};
                    const float pm{sin_sample(arg,invT)*sym};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1                
                     const float  I{ceph_cosf(pm)*cos_carrier_0};
                     const float  Q{ceph_sinf(pm)*sin_carrier_0};
#else 
                     const float  I{std::cos(pm)*cos_carrier_0};
                     const float  Q{std::sin(pm)*sin_carrier_0};
#endif
                     const std::complex<float> sample{I,-Q};
                     sum0 += sample;
                }
                this->m_sig_samples.m_data[__j] = sum0;
            }
            return (0);
}

/*
  Defaulted to uniform random noise sample.
*/
std::int32_t 
gms::radiolocation
::pm_bb_cmplx_sin_signal_t
::create_noisy_signal_user_data(pm_bb_cmplx_sin_signal_pdf_params_t & pdf_params,
                                const float scale,
                                const float * __restrict__ sym_in, // size of m_nsamples*m_nK values [0,1]
                                const std::uint32_t n_T,
                                const std::uint32_t n_K)
{
      using namespace gms::math;
      if(__builtin_expect(static_cast<std::uint32_t>(this->m_nsamples)!=n_T,0) || 
         __builtin_expect(this->m_nK!=n_K,0)) { return (-1);}

      const float T{static_cast<float>(this->m_nsamples)};
      const float invT{this->m_P/T};
      const float invTc{1.0f/T};
      constexpr float C141421356237309504880168872421{1.41421356237309504880168872421f};
      constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
      thread_local std::uint64_t tsc{};
      tsc = __rdtsc();
      auto uni_rand{std::bind(std::uniform_real_distribution<float>(pdf_params.uni_real_a_r,pdf_params.uni_real_b_r),
                                                                    std::mt19937(tsc))};
      std::complex<float> sum;
      for(std::uint32_t __t{0}; __t != n_T; ++__t) 
      {
            const float t{static_cast<float>(__t)};
            const float carrier_arg{this->m_ph+(C6283185307179586476925286766559*this->m_w0*t*invTc)};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1
            const float cos_carrier{ceph_cosf(carrier_arg)};
            const float sin_carrier{ceph_sinf(carrier_arg)};
#else 
            const float cos_carrier{std::cos(carrier_arg)};
            const float sin_carrier{std::sin(carrier_arg)};
#endif
            sum = (0.0f,0.0f);
            for(std::uint32_t __k{0}; __k != n_K; ++__k) 
            {
                 const float k{static_cast<float>(__k)};
                 const float arg{t-k*T};
                 const float sym{sym_in[Ix2D(__t,n_K,__k)]};
                 const float pm{sin_sample(arg,invT)*sym};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1                
                 const float  I{ceph_cosf(pm)*cos_carrier};
                 const float  Q{ceph_sinf(pm)*sin_carrier};
#else 
                 const float  I{std::cos(pm)*cos_carrier};
                 const float  Q{std::sin(pm)*sin_carrier};
#endif
                 const std::complex<float> sample{I,-Q};
                 sum += sample;
            }
            this->m_sig_samples.m_data[__t] = sum;
        }
        return (0);
}

std::int32_t 
gms::radiolocation
::pm_bb_cmplx_sin_signal_t
::create_noisy_signal_user_data(pm_bb_cmplx_sin_signal_pdf_params_t & pdf_params,
                                pm_bb_cmplx_sin_signal_rand_distr rd_enum,
                                const float scale,
                                const float * __restrict__ sym_in, // size of m_nsamples*m_nK values [0,1]
                                const std::uint32_t n_T,
                                const std::uint32_t n_K)
{
     using namespace gms::math;
     if(__builtin_expect(static_cast<std::uint32_t>(this->m_nsamples)!=n_T,0) || 
         __builtin_expect(this->m_nK!=n_K,0)) { return (-1);}

      const float T{static_cast<float>(this->m_nsamples)};
      const float invT{this->m_P/T};
      const float invTc{1.0f/T};
      constexpr float C141421356237309504880168872421{1.41421356237309504880168872421f};
      constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
      thread_local std::uint64_t tsc{};
      std::complex<float> sum;
      switch (rd_enum)
      {
        case pm_bb_cmplx_sin_signal_rand_distr::uniform : 
              tsc = __rdtsc();
              auto uni_rand{std::bind(std::uniform_real_distribution<float>(pdf_params.uni_real_a_r,pdf_params.uni_real_b_r),std::mt19937(tsc))};
              for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
              {
                     const float t{static_cast<float>(__t)};
                     const float carrier_arg{this->m_ph+(C6283185307179586476925286766559*this->m_w0*t*invTc)};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1
                     const float cos_carrier{ceph_cosf(carrier_arg)};
                     const float sin_carrier{ceph_sinf(carrier_arg)};
#else 
                     const float cos_carrier{std::cos(carrier_arg)};
                     const float sin_carrier{std::sin(carrier_arg)};
#endif
                     sum =  (0.0f,0.0f);
                     for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                     {
                          const float k{static_cast<float>(__k)};
                          const float arg{t-k*T};
                          const float sym{sym_in[Ix2D(__t,n_K,__k)]};
                          const float pm{sin_sample_noisy(arg,invT,scale,uni_rand)*sym};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1                
                          const float  I{ceph_cosf(pm)*cos_carrier};
                          const float  Q{ceph_sinf(pm)*sin_carrier};
#else 
                          const float  I{std::cos(pm)*cos_carrier};
                          const float  Q{std::sin(pm)*sin_carrier};
#endif
                          const std::complex<float> sample{I,-Q};
                          sum += sample;
                     }
                     this->m_sig_samples.m_data[__t] = sum;
               }
        break;
        case pm_bb_cmplx_sin_signal_rand_distr::normal : 
              tsc = __rdtsc();
              auto normal_rand{std::bind(std::normal_distribution<float>(pdf_params.norm_mu_r,pdf_params.norm_sigma_r),
                                                                         std::mt19937(tsc))};
              for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
              {
                     const float t{static_cast<float>(__t)};
                     const float carrier_arg{this->m_ph+(C6283185307179586476925286766559*this->m_w0*t*invTc)};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1
                     const float cos_carrier{ceph_cosf(carrier_arg)};
                     const float sin_carrier{ceph_sinf(carrier_arg)};
#else 
                     const float cos_carrier{std::cos(carrier_arg)};
                     const float sin_carrier{std::sin(carrier_arg)};
#endif
                     sum =  (0.0f,0.0f);
                     for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                     {
                          const float k{static_cast<float>(__k)};
                          const float arg{t-k*T};
                          const float sym{sym_in[Ix2D(__t,n_K,__k)]};
                          const float pm{sin_sample_noisy(arg,invT,scale,normal_rand)*sym};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1                
                          const float  I{ceph_cosf(pm)*cos_carrier};
                          const float  Q{ceph_sinf(pm)*sin_carrier};
#else 
                          const float  I{std::cos(pm)*cos_carrier};
                          const float  Q{std::sin(pm)*sin_carrier};
#endif
                          const std::complex<float> sample{I,-Q};
                          sum += sample;
                     }
                     this->m_sig_samples.m_data[__t] = sum;
               }
         break;
         case pm_bb_cmplx_sin_signal_rand_distr::cauchy : 
                tsc = __rdtsc();
                auto cauchy_rand{std::bind(std::cauchy_distribution<float>(pdf_params.cauchy_a_r,pdf_params.cauchy_b_r),
                                                                           std::mt19937(tsc))};
                for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
                {
                     const float t{static_cast<float>(__t)};
                     const float carrier_arg{this->m_ph+(C6283185307179586476925286766559*this->m_w0*t*invTc)};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1
                     const float cos_carrier{ceph_cosf(carrier_arg)};
                     const float sin_carrier{ceph_sinf(carrier_arg)};
#else 
                     const float cos_carrier{std::cos(carrier_arg)};
                     const float sin_carrier{std::sin(carrier_arg)};
#endif
                     sum =  (0.0f,0.0f);
                     for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                     {
                          const float k{static_cast<float>(__k)};
                          const float arg{t-k*T};
                          const float sym{sym_in[Ix2D(__t,n_K,__k)]};
                          const float pm{sin_sample_noisy(arg,invT,scale,cauchy_rand)*sym};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1                
                          const float  I{ceph_cosf(pm)*cos_carrier};
                          const float  Q{ceph_sinf(pm)*sin_carrier};
#else 
                          const float  I{std::cos(pm)*cos_carrier};
                          const float  Q{std::sin(pm)*sin_carrier};
#endif
                          const std::complex<float> sample{I,-Q};
                          sum += sample;
                     }
                     this->m_sig_samples.m_data[__t] = sum;
               }
         break;
         case pm_bb_cmplx_sin_signal_rand_distr::log_norm : 
               tsc = __rdtsc();
               auto lognorm_rand{std::bind(std::lognormal_distribution<float>(pdf_params.log_norm_m_r,pdf_params.log_norm_s_r),
                                                                               std::mt19937(tsc))};
                for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
                {
                     const float t{static_cast<float>(__t)};
                     const float carrier_arg{this->m_ph+(C6283185307179586476925286766559*this->m_w0*t*invTc)};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1
                     const float cos_carrier{ceph_cosf(carrier_arg)};
                     const float sin_carrier{ceph_sinf(carrier_arg)};
#else 
                     const float cos_carrier{std::cos(carrier_arg)};
                     const float sin_carrier{std::sin(carrier_arg)};
#endif
                     sum =  (0.0f,0.0f);
                     for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                     {
                          const float k{static_cast<float>(__k)};
                          const float arg{t-k*T};
                          const float sym{sym_in[Ix2D(__t,n_K,__k)]};
                          const float pm{sin_sample_noisy(arg,invT,scale,lognorm_rand)*sym};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1                
                          const float  I{ceph_cosf(pm)*cos_carrier};
                          const float  Q{ceph_sinf(pm)*sin_carrier};
#else 
                          const float  I{std::cos(pm)*cos_carrier};
                          const float  Q{std::sin(pm)*sin_carrier};
#endif
                          const std::complex<float> sample{I,-Q};
                          sum += sample;
                     }
                     this->m_sig_samples.m_data[__t] = sum;
               }
         break;
         case pm_bb_cmplx_sin_signal_rand_distr::expo_gamma : 
               tsc = __rdtsc();
               auto expo_rand{std::bind(std::exponential_distribution<float>(pdf_params.expo_gamma_r),
                                                                             std::mt19937(tsc))};
                for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
                {
                     const float t{static_cast<float>(__t)};
                     const float carrier_arg{this->m_ph+(C6283185307179586476925286766559*this->m_w0*t*invTc)};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1
                     const float cos_carrier{ceph_cosf(carrier_arg)};
                     const float sin_carrier{ceph_sinf(carrier_arg)};
#else 
                     const float cos_carrier{std::cos(carrier_arg)};
                     const float sin_carrier{std::sin(carrier_arg)};
#endif
                     sum =  (0.0f,0.0f);
                     for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                     {
                          const float k{static_cast<float>(__k)};
                          const float arg{t-k*T};
                          const float sym{sym_in[Ix2D(__t,n_K,__k)]};
                          const float pm{sin_sample_noisy(arg,invT,scale,expo_rand)*sym};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1                
                          const float  I{ceph_cosf(pm)*cos_carrier};
                          const float  Q{ceph_sinf(pm)*sin_carrier};
#else 
                          const float  I{std::cos(pm)*cos_carrier};
                          const float  Q{std::sin(pm)*sin_carrier};
#endif
                          const std::complex<float> sample{I,-Q};
                          sum += sample;
                     }
                     this->m_sig_samples.m_data[__t] = sum;
               }
          break;
          case pm_bb_cmplx_sin_signal_rand_distr::weibull : 
                 tsc = __rdtsc();
                 auto weibull_rand{std::bind(std::weibull_distribution<float>(pdf_params.weibull_a_r,pdf_params.weibull_b_r),
                                                                              std::mt19937(tsc))};
                 for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
                 {
                     const float t{static_cast<float>(__t)};
                     const float carrier_arg{this->m_ph+(C6283185307179586476925286766559*this->m_w0*t*invTc)};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1
                     const float cos_carrier{ceph_cosf(carrier_arg)};
                     const float sin_carrier{ceph_sinf(carrier_arg)};
#else 
                     const float cos_carrier{std::cos(carrier_arg)};
                     const float sin_carrier{std::sin(carrier_arg)};
#endif                     
                     sum =  (0.0f,0.0f);
                     for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                     {
                          const float k{static_cast<float>(__k)};
                          const float arg{t-k*T};
                          const float sym{sym_in[Ix2D(__t,n_K,__k)]};
                          const float pm{sin_sample_noisy(arg,invT,scale,weibull_rand)*sym};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1                
                          const float  I{ceph_cosf(pm)*cos_carrier};
                          const float  Q{ceph_sinf(pm)*sin_carrier};
#else 
                          const float  I{std::cos(pm)*cos_carrier};
                          const float  Q{std::sin(pm)*sin_carrier};
#endif
                          const std::complex<float> sample{I,-Q};
                          sum += sample;
                     }
                     this->m_sig_samples.m_data[__t] = sum;
               }
          break;
          case pm_bb_cmplx_sin_signal_rand_distr::gamma : 
                 tsc = __rdtsc();
                 auto gamma_rand{std::bind(std::gamma_distribution<float>(pdf_params.gamma_alph_r,pdf_params.gamma_bet_r),
                                                                          std::mt19937(tsc))};
                 for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
                 {
                     const float t{static_cast<float>(__t)};
                     const float carrier_arg{this->m_ph+(C6283185307179586476925286766559*this->m_w0*t*invTc)};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1
                     const float cos_carrier{ceph_cosf(carrier_arg)};
                     const float sin_carrier{ceph_sinf(carrier_arg)};
#else 
                     const float cos_carrier{std::cos(carrier_arg)};
                     const float sin_carrier{std::sin(carrier_arg)};
#endif                  
                     sum =  (0.0f,0.0f);
                     for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                     {
                          const float k{static_cast<float>(__k)};
                          const float arg{t-k*T};
                          const float sym{sym_in[Ix2D(__t,n_K,__k)]};
                          const float pm{sin_sample_noisy(arg,invT,scale,gamma_rand)*sym};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1                
                          const float  I{ceph_cosf(pm)*cos_carrier};
                          const float  Q{ceph_sinf(pm)*sin_carrier};
#else 
                          const float  I{std::cos(pm)*cos_carrier};
                          const float  Q{std::sin(pm)*sin_carrier};
#endif
                          const std::complex<float> sample{I,-Q};
                          sum += sample;
                     }
                     this->m_sig_samples.m_data[__t] = sum;
               }
           break;
              default : 
                  return (-2);
      }
      return (0);
}

auto
gms::radiolocation
::operator<<(std::ostream &os,
             gms::radiolocation::pm_bb_cmplx_sin_signal_t &rhs)->std::ostream &
{
    std::cout << typeid(rhs).name() << "Begin: object state dump." << std::endl;
    std::cout << "m_nsamples : "      << rhs.m_nsamples << std::endl;
    std::cout << "m_nK       : "      << rhs.m_nK       << std::endl;
    std::cout << "m_w0       : "      << std::fixed << std::setprecision(7) << rhs.m_w0 << std::endl;
    std::cout << "m_ph       : "      << std::fixed << std::setprecision(7) << rhs.m_ph << std::endl;
    std::cout << "m_n        : "      << std::fixed << std::setprecision(7) << rhs.m_n << std::endl;
    std::cout << "m_A        : "      << std::fixed << std::setprecision(7) << rhs.m_A << std::endl;
    std::cout << "m_P        : "      << std::fixed << std::setprecision(7) << rhs.m_P << std::endl;
    std::cout << "Signal-samples:" << std::endl;
    for(std::size_t __i{0ull}; __i != rhs.m_nsamples; ++__i)
    {
        os << std::fixed << std::setprecision(7) << "("
                          << rhs.m_sig_samples.m_data[__i].real() << ",+j" << rhs.m_sig_samples.m_data[__i].imag() << ")" << std::endl;
    }
    std::cout << typeid(rhs).name() << "End: object state dump." << std::endl;
    return (os);
}

