#include <fstream>
#include <iomanip>
#include "GMS_am_bb_square_signal_fp16.h"
#include "GMS_sse_memset.h"
#include "GMS_indices.h"

gms::radiolocation
::am_bb_square_signal_fp16_t
::am_bb_square_signal_fp16_t(const std::size_t nsamples,
                             const std::uint32_t nK,
                             const float n,
                             const float A,
                             const float P)
:
m_nsamples{nsamples},
m_nK{nK},
m_n{n},
m_A{A},
m_P{P},
m_sig_samples{darray_r2_t(m_nsamples)}
{

}

gms::radiolocation
::am_bb_square_signal_fp16_t
::am_bb_square_signal_fp16_t(am_bb_square_signal_fp16_t &&other) noexcept(true)
:
m_nsamples{std::move(other.m_nsamples)},
m_nK{std::move(other.m_nK)},
m_n{std::move(other.m_n)},
m_A{std::move(other.m_A)},
m_P{std::move(other.m_P)},
m_sig_samples{std::move(other.m_sig_samples)}
{

}

gms::radiolocation
::am_bb_square_signal_fp16_t
::am_bb_square_signal_fp16_t(const am_bb_square_signal_fp16_t &other) noexcept(false)
:
m_nsamples{other.m_nsamples},
m_nK{other.m_nK},
m_n{other.m_n},
m_A{other.m_A},
m_P{other.m_P},
m_sig_samples{other.m_sig_samples}
{

}

gms::radiolocation
::am_bb_square_signal_fp16_t
::~am_bb_square_signal_fp16_t()
{

}

gms::radiolocation
::am_bb_square_signal_fp16_t &
gms::radiolocation
::am_bb_square_signal_fp16_t
::operator=(am_bb_square_signal_fp16_t &&other) noexcept(true)
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_nsamples  = std::move(other.m_nsamples);
    this->m_nK        = std::move(other.m_nK);
    this->m_n         = std::move(other.m_n);
    this->m_A         = std::move(other.m_A);
    this->m_P         = std::move(other.m_P);
    this->m_sig_samples.operator=(std::move(other.m_sig_samples));
    return (*this);
}

void 
gms::radiolocation
::am_bb_square_signal_fp16_t
::init_storage(const half_float::half filler)
{         
#if (INIT_BY_STD_FILL) == 0
     using namespace gms::common;
	 sse_memset_unroll8x_ps(reinterpret_cast<float*>(this->m_sig_samples.m_data),filler,this->m_nsamples);
#else 
     std::fill(this->m_sig_samples.m_data,this->m_sig_samples.m_data+this->m_nsamples,filler);
#endif 
}

void 
gms::radiolocation 
::am_bb_square_signal_fp16_t
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
::am_bb_square_signal_fp16_t
::create_signal_user_data(const float * __restrict__ sym_in, // size of m_nsamples*m_nK
                          const std::uint32_t sin_or_cos,
                          const std::uint32_t n_T,
                          const std::uint32_t n_K)
{
      using namespace half_float;
      if(__builtin_expect(static_cast<std::uint32_t>(this->m_nsamples)!=n_T,0) || 
         __builtin_expect(this->m_nK!=n_K,0)) { return (-1);}
      
       const float T{static_cast<float>(this->m_nsamples)};
       const float invT{this->m_P/T};
       float sum;
       switch (sin_or_cos)
       {
          case 0 : 
             
              for(std::uint32_t __t{0}; __t != n_T; ++__t) 
              {
                   const float t{static_cast<float>(__t)};
                   sum = 0.0f;
                   for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                   {
                        const float k{static_cast<float>(__k)};
                        const float arg{t-k*T};
                        sum += sin_squarew_sample(arg,invT)*sym_in[Ix2D(__t,n_K,__k)];
                   }
                   this->m_sig_samples.m_data[__t] = half_cast<half>(sum);
               }
          break;
          case 1 :
              
              for(std::uint32_t __t{0}; __t != n_T; ++__t) 
              {
                   const float t{static_cast<float>(__t)};
                   sum = 0.0f;
                   for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                   {
                        const float k{static_cast<float>(__k)};
                        const float arg{t-k*T};
                        sum += cos_squarew_sample(arg,invT)*sym_in[Ix2D(__t,n_K,__k)];
                   }
                   this->m_sig_samples.m_data[__t] = half_cast<half>(sum);
               }
            break;
            default : 
                return (-2);
      }
      return (0);
}

std::int32_t 
gms::radiolocation
::am_bb_square_signal_fp16_t
::create_signal_user_data_u4x(const float * __restrict__ sym_in, // size of m_nsamples*m_nK
                              const std::uint32_t sin_or_cos,
                              const std::uint32_t n_T,
                              const std::uint32_t n_K)
{
       using namespace half_float;
       if(__builtin_expect(static_cast<std::uint32_t>(this->m_nsamples)!=n_T,0) || 
          __builtin_expect(this->m_nK!=n_K,0)) { return (-1);}
      
        const float T{static_cast<float>(this->m_nsamples)};
        const float invT{this->m_P/T};
        float sum0;
        float sum1;
        float sum2;
        float sum3;
        float t__i_0{-4.0f};
        float t__i_1{1.0f};
        float t__i_2{2.0f};
        float t__i_3{3.0f};
        std::uint32_t __i;
        std::uint32_t __j;

        switch (sin_or_cos)
        {
              case 0 : 
                  
                  for(__i = 0; __i != ROUND_TO_FOUR(n_T,4); __i += 4u)
                  {
                        t__i_0 += 4.0f;
                        sum0    = 0.0f;
                        for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                        {
                            const float k{static_cast<float>(__k)};
                            const float arg{t__i_0-k*T};
                            sum0 += sin_squarew_sample(arg,invT)*sym_in[Ix2D(__i+0,n_K,__k)];
                        }
                        this->m_sig_samples.m_data[__i+0] = half_cast<half>(sum0);
                        t__i_1 += 4.0f;
                        sum1    = 0.0f;
                        for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                        {
                            const float k{static_cast<float>(__k)};
                            const float arg{t__i_1-k*T};
                            sum1 += sin_squarew_sample(arg,invT)*sym_in[Ix2D(__i+1,n_K,__k)];
                        }
                        this->m_sig_samples.m_data[__i+1] = half_cast<half>(sum1);
                        t__i_2  += 4.0f;
                        sum2     = 0.0f;
                        for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                        {
                            const float k{static_cast<float>(__k)};
                            const float arg{t__i_2-k*T};
                            sum2 += sin_squarew_sample(arg,invT)*sym_in[Ix2D(__i+2,n_K,__k)];
                        }
                        this->m_sig_samples.m_data[__i+2] = half_cast<half>(sum2);
                        t__i_3  += 4.0f;
                        sum3     = 0.0f;
                        for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                        {
                            const float k{static_cast<float>(__k)};
                            const float arg{t__i_3-k*T};
                            sum3 += sin_squarew_sample(arg,invT)*sym_in[Ix2D(__i+3,n_K,__k)];
                        }
                        this->m_sig_samples.m_data[__i+3] = half_cast<half>(sum3);
                  }
                  for(__j = __i; __j != n_T; ++__j) 
                  {
                      const float t{static_cast<float>(__j)};
                      sum0 = 0.0f;
                      for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                      {
                           const float k{static_cast<float>(__k)};
                           const float arg{t-k*T};
                           sum0 += sin_squarew_sample(arg,invT)*sym_in[Ix2D(__j,n_K,__k)];
                      }
                       this->m_sig_samples.m_data[__j] = half_cast<half>(sum0);
                  }
             break;
             case 1 : 

                   for(__i = 0; __i != ROUND_TO_FOUR(n_T,4); __i += 4u)
                   {
                        t__i_0 += 4.0f;
                        sum0    = 0.0f;
                        for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                        {
                            const float k{static_cast<float>(__k)};
                            const float arg{t__i_0-k*T};
                            sum0 += cos_squarew_sample(arg,invT)*sym_in[Ix2D(__i+0,n_K,__k)];
                        }
                        this->m_sig_samples.m_data[__i+0] = half_cast<half>(sum0);
                        t__i_1 += 4.0f;
                        sum1    = 0.0f;
                        for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                        {
                            const float k{static_cast<float>(__k)};
                            const float arg{t__i_1-k*T};
                            sum1 += cos_squarew_sample(arg,invT)*sym_in[Ix2D(__i+1,n_K,__k)];
                        }
                        this->m_sig_samples.m_data[__i+1] = half_cast<half>(sum1);
                        t__i_2  += 4.0f;
                        sum2     = 0.0f;
                        for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                        {
                            const float k{static_cast<float>(__k)};
                            const float arg{t__i_2-k*T};
                            sum2 += cos_squarew_sample(arg,invT)*sym_in[Ix2D(__i+2,n_K,__k)];
                        }
                        this->m_sig_samples.m_data[__i+2] = half_cast<half>(sum2);
                        t__i_3  += 4.0f;
                        sum3     = 0.0f;
                        for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                        {
                            const float k{static_cast<float>(__k)};
                            const float arg{t__i_3-k*T};
                            sum3 += cos_squarew_sample(arg,invT)*sym_in[Ix2D(__i+3,n_K,__k)];
                        }
                        this->m_sig_samples.m_data[__i+3] = half_cast<half>(sum3);
                  }
                  for(__j = __i; __j != n_T; ++__j) 
                  {
                      const float t{static_cast<float>(__j)};
                      sum0 = 0.0f;
                      for(std::uint32_t __k{0}; __k != n_K; ++__k) 
                      {
                           const float k{static_cast<float>(__k)};
                           const float arg{t-k*T};
                           sum0 += cos_squarew_sample(arg,invT)*sym_in[Ix2D(__j,n_K,__k)];
                      }
                       this->m_sig_samples.m_data[__j] = half_cast<half>(sum0);
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
             gms::radiolocation::am_bb_square_signal_fp16_t &rhs)->std::ostream &
{
    std::cout << typeid(rhs).name() << "Begin: object state dump." << std::endl;
    std::cout << "m_nsamples : "      << rhs.m_nsamples << std::endl;
    std::cout << "m_nK       : "      << rhs.m_nK       << std::endl;
    std::cout << "m_n        : "      << std::fixed << std::setprecision(7) << rhs.m_n << std::endl;
    std::cout << "m_A        : "      << std::fixed << std::setprecision(7) << rhs.m_A << std::endl;
    std::cout << "m_P        : "      << std::fixed << std::setprecision(7) << rhs.m_P << std::endl;
    std::cout << "Signal-samples:" << std::endl;
    for(std::size_t __i{0ull}; __i != rhs.m_nsamples; ++__i)
    {
         os << std::fixed << std::setprecision(4) << 
                             rhs.m_sig_samples.m_data[__i] << std::endl;
    }
    std::cout << typeid(rhs).name() << "End: object state dump." << std::endl;
    return (os);
}

