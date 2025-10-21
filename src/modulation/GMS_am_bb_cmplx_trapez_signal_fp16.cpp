#include <fstream>
#include <iomanip>
#include <random>
#include "GMS_am_bb_cmplx_trapez_signal_fp16.h"
#include "GMS_sse_memset.h"
#include "GMS_indices.h"


gms::radiolocation 
::am_bb_cmplx_trapez_signal_fp16_t
::am_bb_cmplx_trapez_signal_fp16_t(const std::size_t nsamples,
                              const std::uint32_t nK,
                              const float a,
                              const float l,
                              const float c,
                              const float m)
:
m_nsamples{nsamples},
m_nK{nK},
m_a{a},
m_l{l},
m_c{c},
m_m{m},
m_sig_samples{darray_c2_t(m_nsamples)}
{

}

gms::radiolocation 
::am_bb_cmplx_trapez_signal_fp16_t
::am_bb_cmplx_trapez_signal_fp16_t(am_bb_cmplx_trapez_signal_fp16_t &&other) noexcept(true)
:
m_nsamples{std::move(other.m_nsamples)},
m_nK{std::move(other.m_nK)},
m_a{std::move(other.m_a)},
m_l{std::move(other.m_l)},
m_c{std::move(other.m_c)},
m_m{std::move(other.m_m)},
m_sig_samples{std::move(other.m_sig_samples)}
{

}

gms::radiolocation 
::am_bb_cmplx_trapez_signal_fp16_t
::am_bb_cmplx_trapez_signal_fp16_t(const am_bb_cmplx_trapez_signal_fp16_t &other) noexcept(false)
:
m_nsamples{other.m_nsamples},
m_nK{other.m_nK},
m_a{other.m_a},
m_l{other.m_l},
m_c{other.m_c},
m_m{other.m_m},
m_sig_samples{other.m_sig_samples}
{

}

gms::radiolocation 
::am_bb_cmplx_trapez_signal_fp16_t
::~am_bb_cmplx_trapez_signal_fp16_t()
{

}

gms::radiolocation 
::am_bb_cmplx_trapez_signal_fp16_t & 
gms::radiolocation 
::am_bb_cmplx_trapez_signal_fp16_t
::operator=(am_bb_cmplx_trapez_signal_fp16_t &&other) noexcept(true)
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_nsamples            = std::move(other.m_nsamples);
    this->m_nK                  = std::move(other.m_nK);
    this->m_a                   = std::move(other.m_a);
    this->m_l                   = std::move(other.m_l);
    this->m_c                   = std::move(other.m_c);
    this->m_m                   = std::move(other.m_m);
    this->m_sig_samples.operator=(std::move(other.m_sig_samples));
    return (*this);
}

void 
gms::radiolocation
::am_bb_cmplx_trapez_signal_fp16_t
::init_storage(const half_float::half filler)
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
::am_bb_cmplx_trapez_signal_fp16_t
::create_signal_plot(const std::uint32_t n_samp,
                    const half_float::half * __restrict sig_arg,
                    const half_float::half * __restrict sig_val,
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
::am_bb_cmplx_trapez_signal_fp16_t
::create_signal_user_data(const float * __restrict__ sym_in, // size of m_nsamples*m_nK values [0,1]
                          const std::uint32_t n_T,
                          const std::uint32_t n_K)
{
      using namespace half_float;
      if(__builtin_expect(static_cast<std::uint32_t>(this->m_nsamples)!=n_T,0) || 
         __builtin_expect(this->m_nK!=n_K,0)) { return (-1);}

      const float T{static_cast<float>(this->m_nsamples)};
      constexpr float C141421356237309504880168872421{1.41421356237309504880168872421f};
      std::complex<float> sum;
      for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
      {
            const float t{static_cast<float>(__t)};
            sum =  (0.0f,0.0f);
            for(std::uint32_t __k{0}; __k != n_K; ++__k) 
            {
                const float k{static_cast<float>(__k)};
                const float arg{t-k*T};
                const float sym{sym_in[Ix2D(__t,n_K,__k)]};
                const float re_im{1.0f-(2.0f*sym)};
                const std::complex<float> c_sym{re_im,re_im};
                sum += trapezoid_sample(arg)*(C141421356237309504880168872421*c_sym);
                
            }
            std::complex<half> sum_fp16(half_cast<half>(sum.real()),half_cast<half>(sum.imag()));
            this->m_sig_samples.m_data[__t] = sum_fp16;
      }
      return (0);
}

std::int32_t 
gms::radiolocation
::am_bb_cmplx_trapez_signal_fp16_t
::create_signal_user_data_u4x(const float * __restrict__ sym_in, // size of __I_n_samples__ * __I_n_K__
                              const std::uint32_t n_T,
                              const std::uint32_t n_K)
{
      using namespace half_float;
      if(__builtin_expect(static_cast<std::uint32_t>(this->m_nsamples)!=n_T,0) || 
         __builtin_expect(this->m_nK!=n_K,0)) { return (-2);}

      const float T{static_cast<float>(this->m_nsamples)};
      constexpr float C141421356237309504880168872421{1.41421356237309504880168872421f};
      std::complex<float> sum0;
      std::complex<float> sum1;
      std::complex<float> sum2;
      std::complex<float> sum3;
      float t__i_0{-4.0f};
      float t__i_1{1.0f};
      float t__i_2{2.0f};
      float t__i_3{3.0f};
      std::uint32_t __i;
      std::uint32_t __j;
      for(__i = 0; __i != ROUND_TO_FOUR(n_T,4); __i += 4) 
      {
           t__i_0 += 4.0f;
           sum0   = (0.0f,0.0f);
           for(std::uint32_t __k{0}; __k != n_K; ++__k) 
           {
                const float k{static_cast<float>(__k)};
                const float arg{t__i_0-k*T};
                const float sym{sym_in[Ix2D(__i+0,n_K,__k)]};
                const float im_re{1.0f-(2.0f*sym)};
                const std::complex<float> c_sym{im_re,im_re};
                sum0 += trapezoid_sample(arg)*(C141421356237309504880168872421*c_sym);
               
           }
           std::complex<half> tmp_sum0(half_cast<half>(sum0.real()),half_cast<half>(sum0.imag()));
           this->m_sig_samples.m_data[__i+0] = tmp_sum0;
           t__i_1 += 4.0f;
           sum1   = (0.0f,0.0f);
           for(std::uint32_t __k{0}; __k != n_K; ++__k) 
           {
                const float k{static_cast<float>(__k)};
                const float arg{t__i_1-k*T};
                const float sym{sym_in[Ix2D(__i+1,n_K,__k)]};
                const float im_re{1.0f-(2.0f*sym)};
                const std::complex<float> c_sym{im_re,im_re};
                sum1 += trapezoid_sample(arg)*(C141421356237309504880168872421*c_sym);
                            
           }
           std::complex<half> tmp_sum1(half_cast<half>(sum1.real()),half_cast<half>(sum1.imag()));
           this->m_sig_samples.m_data[__i+1] = tmp_sum1;
           t__i_2 += 4.0f;
           sum2    = (0.0f,0.0f);
           for(std::uint32_t __k{0}; __k != n_K; ++__k) 
           {
                const float k{static_cast<float>(__k)};
                const float arg{t__i_2-k*T};
                const float sym{sym_in[Ix2D(__i+2,n_K,__k)]};
                const float im_re{1.0f-(2.0f*sym)};
                const std::complex<float> c_sym{im_re,im_re};
                sum2 += trapezoid_sample(arg)*(C141421356237309504880168872421*c_sym);
                
           }
           std::complex<half> tmp_sum2(half_cast<half>(sum2.real()),half_cast<half>(sum2.imag()));
           this->m_sig_samples.m_data[__i+2] = tmp_sum2;
           t__i_3 += 4.0f;
           sum3   =  (0.0f,0.0f);
           for(std::uint32_t __k{0}; __k != n_K; ++__k) 
           {
                const float k{static_cast<float>(__k)};
                const float arg{t__i_3-k*T};
                const float sym{sym_in[Ix2D(__i+3,n_K,__k)]};
                const float im_re{1.0f-(2.0f*sym)};
                const std::complex<float> c_sym{im_re,im_re};
                sum3 += trapezoid_sample(arg)*(C141421356237309504880168872421*c_sym);
                     
           }
           std::complex<half> tmp_sum3(half_cast<half>(sum3.real()),half_cast<half>(sum3.imag()));
           this->m_sig_samples.m_data[__i+3] = tmp_sum3;
       }
       for(__j = __i; __j != n_T; ++__j) 
       {
            const float t{static_cast<float>(__j)};
            sum0 = (0.0f,0.0f);
            for(std::uint32_t __k{0}; __k != n_K; ++__k) 
            {
                const float k{static_cast<float>(__k)};
                const float arg{t-k*T};
                const float sym{sym_in[Ix2D(__j,n_K,__k)]};
                const float im_re{1.0f-(2.0f*sym)};
                const std::complex<float> c_sym{im_re,im_re};
                sum0 += trapezoid_sample(arg)*(C141421356237309504880168872421*c_sym);
               
            }
            std::complex<half> sum0_rem(half_cast<half>(sum0.real()),half_cast<half>(sum0.imag()));
            this->m_sig_samples.m_data[__j] = sum0_rem;
        }
        return (0);
}

std::int32_t 
gms::radiolocation
::am_bb_cmplx_trapez_signal_fp16_t
::create_noisy_signal_user_data(am_bb_cmplx_trapez_signal_fp16_pdf_params_t & pdf_params,
                                const float scale,
                                const float * __restrict__ sym_in, // size of m_nsamples*m_nK values [0,1]
                                const std::uint32_t n_T,
                                const std::uint32_t n_K)
{
      using namespace half_float;
      if(__builtin_expect(static_cast<std::uint32_t>(this->m_nsamples)!=n_T,0) || 
         __builtin_expect(this->m_nK!=n_K,0)) { return (-1);}

      const float T{static_cast<float>(this->m_nsamples)};
      constexpr float C141421356237309504880168872421{1.41421356237309504880168872421f};
      thread_local std::uint64_t tsc{};
      tsc = __rdtsc();
      auto uni_rand{std::bind(std::uniform_real_distribution<float>(pdf_params.uni_real_a_r,pdf_params.uni_real_b_r),std::mt19937(tsc))};
      std::complex<float> sum;
      
      for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
      {
            const float t{static_cast<float>(__t)};
            sum =  (0.0f,0.0f);
            for(std::uint32_t __k{0}; __k != n_K; ++__k) 
            {
                const float k{static_cast<float>(__k)};
                const float arg{t-k*T};
                const float sym{sym_in[Ix2D(__t,n_K,__k)]};
                const float re_im{1.0f-(2.0f*sym)};
                const std::complex<float> c_sym{re_im,re_im};
                sum += trapezoid_sample_noisy(arg,scale,uni_rand)*(C141421356237309504880168872421*c_sym);
            }
            std::complex<half> tmp_sum(half_cast<half>(sum.real()),half_cast<half>(sum.imag()));
            this->m_sig_samples.m_data[__t] = tmp_sum;
      }
      return (0);
}

std::int32_t 
gms::radiolocation
::am_bb_cmplx_trapez_signal_fp16_t
::create_noisy_signal_user_data(am_bb_cmplx_trapez_signal_fp16_pdf_params_t & pdf_params,
                                am_bb_cmplx_trapez_signal_fp16_rand_distr rd_enum,
                                const float scale,
                                const float * __restrict__ sym_in, // size of m_nsamples*m_nK values [0,1]
                                const std::uint32_t n_T,
                                const std::uint32_t n_K)
{
      using namespace half_float;
      if(__builtin_expect(static_cast<std::uint32_t>(this->m_nsamples)!=n_T,0) || 
         __builtin_expect(this->m_nK!=n_K,0)) { return (-1);}

      const float T{static_cast<float>(this->m_nsamples)};
      constexpr float C141421356237309504880168872421{1.41421356237309504880168872421f};
      thread_local std::uint64_t tsc{};
      
      std::complex<float> sum;
      switch (rd_enum)
      {
        case am_bb_cmplx_trapez_signal_fp16_rand_distr::uniform : 
        {
              tsc = __rdtsc();
              auto uni_rand{std::bind(std::uniform_real_distribution<float>(pdf_params.uni_real_a_r,pdf_params.uni_real_b_r),std::mt19937(tsc))};
              for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
              {
                     const float t{static_cast<float>(__t)};
                     sum =  (0.0f,0.0f);
                     for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                     {
                          const float k{static_cast<float>(__k)};
                          const float arg{t-k*T};
                          const float sym{sym_in[Ix2D(__t,n_K,__k)]};
                          const float re_im{1.0f-(2.0f*sym)};
                          const std::complex<float> c_sym{re_im,re_im};
                          sum += trapezoid_sample_noisy(arg,scale,uni_rand)*(C141421356237309504880168872421*c_sym);
                     }
                     std::complex<half> tmp_sum(half_cast<half>(sum.real()),half_cast<half>(sum.imag()));
                     this->m_sig_samples.m_data[__t] = tmp_sum;
               }
        }
        break;
        case am_bb_cmplx_trapez_signal_fp16_rand_distr::normal : 
        {
              tsc = __rdtsc();
              auto normal_rand{std::bind(std::normal_distribution<float>(pdf_params.norm_mu_r,pdf_params.norm_sigma_r),
                                                                         std::mt19937(tsc))};
              for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
              {
                     const float t{static_cast<float>(__t)};
                     sum =  (0.0f,0.0f);
                     for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                     {
                          const float k{static_cast<float>(__k)};
                          const float arg{t-k*T};
                          const float sym{sym_in[Ix2D(__t,n_K,__k)]};
                          const float re_im{1.0f-(2.0f*sym)};
                          const std::complex<float> c_sym{re_im,re_im};
                          sum += trapezoid_sample_noisy(arg,scale,normal_rand)*(C141421356237309504880168872421*c_sym);
                     }
                     std::complex<half> tmp_sum(half_cast<half>(sum.real()),half_cast<half>(sum.imag()));
                     this->m_sig_samples.m_data[__t] = tmp_sum;
               }
        }
         break;
         case am_bb_cmplx_trapez_signal_fp16_rand_distr::cauchy : 
         {
                tsc = __rdtsc();
                auto cauchy_rand{std::bind(std::cauchy_distribution<float>(pdf_params.cauchy_a_r,pdf_params.cauchy_b_r),
                                                                           std::mt19937(tsc))};
                for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
                {
                     const float t{static_cast<float>(__t)};
                     sum =  (0.0f,0.0f);
                     for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                     {
                          const float k{static_cast<float>(__k)};
                          const float arg{t-k*T};
                          const float sym{sym_in[Ix2D(__t,n_K,__k)]};
                          const float re_im{1.0f-(2.0f*sym)};
                          const std::complex<float> c_sym{re_im,re_im};
                          sum += trapezoid_sample_noisy(arg,scale,cauchy_rand)*(C141421356237309504880168872421*c_sym);
                     }
                     std::complex<half> tmp_sum(half_cast<half>(sum.real()),half_cast<half>(sum.imag()));
                     this->m_sig_samples.m_data[__t] = sum;
               }
         }
         break;
         case am_bb_cmplx_trapez_signal_fp16_rand_distr::log_norm : 
         {
               tsc = __rdtsc();
               auto lognorm_rand{std::bind(std::lognormal_distribution<float>(pdf_params.log_norm_m_r,pdf_params.log_norm_s_r),
                                                                               std::mt19937(tsc))};
                for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
                {
                     const float t{static_cast<float>(__t)};
                     sum =  (0.0f,0.0f);
                     for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                     {
                          const float k{static_cast<float>(__k)};
                          const float arg{t-k*T};
                          const float sym{sym_in[Ix2D(__t,n_K,__k)]};
                          const float re_im{1.0f-(2.0f*sym)};
                          const std::complex<float> c_sym{re_im,re_im};
                          sum += trapezoid_sample_noisy(arg,scale,lognorm_rand)*(C141421356237309504880168872421*c_sym);
                     }
                     std::complex<half> tmp_sum(half_cast<half>(sum.real()),half_cast<half>(sum.imag()));
                     this->m_sig_samples.m_data[__t] = tmp_sum;
               }
         }
         break;
         case am_bb_cmplx_trapez_signal_fp16_rand_distr::expo_gamma : 
         {
               tsc = __rdtsc();
               auto expo_rand{std::bind(std::exponential_distribution<float>(pdf_params.expo_gamma_r),
                                                                             std::mt19937(tsc))};
                for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
                {
                     const float t{static_cast<float>(__t)};
                     sum =  (0.0f,0.0f);
                     for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                     {
                          const float k{static_cast<float>(__k)};
                          const float arg{t-k*T};
                          const float sym{sym_in[Ix2D(__t,n_K,__k)]};
                          const float re_im{1.0f-(2.0f*sym)};
                          const std::complex<float> c_sym{re_im,re_im};
                          sum += trapezoid_sample_noisy(arg,scale,expo_rand)*(C141421356237309504880168872421*c_sym);
                     }
                     std::complex<half> tmp_sum(half_cast<half>(sum.real()),half_cast<half>(sum.imag()));
                     this->m_sig_samples.m_data[__t] = tmp_sum;
               }
         }
          break;
          case am_bb_cmplx_trapez_signal_fp16_rand_distr::weibull : 
          {
                 tsc = __rdtsc();
                 auto weibull_rand{std::bind(std::weibull_distribution<float>(pdf_params.weibull_a_r,pdf_params.weibull_b_r),
                                                                              std::mt19937(tsc))};
                 for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
                 {
                     const float t{static_cast<float>(__t)};
                     sum =  (0.0f,0.0f);
                     for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                     {
                          const float k{static_cast<float>(__k)};
                          const float arg{t-k*T};
                          const float sym{sym_in[Ix2D(__t,n_K,__k)]};
                          const float re_im{1.0f-(2.0f*sym)};
                          const std::complex<float> c_sym{re_im,re_im};
                          sum += trapezoid_sample_noisy(arg,scale,weibull_rand)*(C141421356237309504880168872421*c_sym);
                     }
                     std::complex<half> tmp_sum(half_cast<half>(sum.real()),half_cast<half>(sum.imag()));
                     this->m_sig_samples.m_data[__t] = tmp_sum;
               }
          }
          break;
          case am_bb_cmplx_trapez_signal_fp16_rand_distr::gamma :
          { 
                 tsc = __rdtsc();
                 auto gamma_rand{std::bind(std::gamma_distribution<float>(pdf_params.gamma_alph_r,pdf_params.gamma_bet_r),
                                                                          std::mt19937(tsc))};
                 for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
                 {
                     const float t{static_cast<float>(__t)};
                     sum =  (0.0f,0.0f);
                     for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                     {
                          const float k{static_cast<float>(__k)};
                          const float arg{t-k*T};
                          const float sym{sym_in[Ix2D(__t,n_K,__k)]};
                          const float re_im{1.0f-(2.0f*sym)};
                          const std::complex<float> c_sym{re_im,re_im};
                          sum += trapezoid_sample_noisy(arg,scale,gamma_rand)*(C141421356237309504880168872421*c_sym);
                     }
                     std::complex<half> tmp_sum(half_cast<half>(sum.real()),half_cast<half>(sum.imag()));
                     this->m_sig_samples.m_data[__t] = tmp_sum;
               }
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
             gms::radiolocation::am_bb_cmplx_trapez_signal_fp16_t &rhs)->std::ostream &
{
    std::cout << typeid(rhs).name() << "Begin: object state dump." << std::endl;
    std::cout << "m_nsamples : "      << rhs.m_nsamples << std::endl;
    std::cout << "m_nK       : "      << rhs.m_nK       << std::endl;
    std::cout << "m_a        : "      << std::fixed << std::setprecision(7) << rhs.m_a << std::endl;
    std::cout << "m_m        : "      << std::fixed << std::setprecision(7) << rhs.m_m << std::endl;
    std::cout << "m_l        : "      << std::fixed << std::setprecision(7) << rhs.m_l << std::endl;
    std::cout << "m_c        : "      << std::fixed << std::setprecision(7) << rhs.m_c << std::endl;
    std::cout << "Signal-samples:" << std::endl;
    for(std::size_t __i{0ull}; __i != rhs.m_nsamples; ++__i)
    {
         os << std::fixed << std::setprecision(4) << "("
                          << rhs.m_sig_samples.m_data[__i].real() << ",+j" << rhs.m_sig_samples.m_data[__i].imag() << ")" << std::endl;
    }
    std::cout << typeid(rhs).name() << "End: object state dump." << std::endl;
    return (os);
}
