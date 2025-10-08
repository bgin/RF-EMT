
#include <fstream>
#include <iomanip>
#include "GMS_cmplx_trapezw_env.h"
#include "GMS_sse_memset.h"
#include "GMS_indices.h"


gms::radiolocation
::cmplx_trapezw_env_t
::cmplx_trapezw_env_t(const std::size_t   I_n_samples,
                      const std::size_t   Q_n_samples,
                      const std::uint32_t I_n_K,
                      const std::uint32_t Q_n_K,
                      const float         I_a,
                      const float         I_l,
                      const float         I_c,
                      const float         I_m,
                      const float         Q_a,
                      const float         Q_l,
                      const float         Q_c,
                      const float         Q_m)
:
__I_n_samples__{I_n_samples},
__Q_n_samples__{Q_n_samples},
__I_n_K__{I_n_K},
__Q_n_K__{Q_n_K},
__I_a__{I_a},
__I_l__{I_l},
__I_c__{I_c},
__I_m__{I_m},
__Q_a__{Q_a},
__Q_l__{Q_l},
__Q_c__{Q_c},
__Q_m__{Q_m},
__I_chan__{darray_r4_t(__I_n_samples__)},
__Q_chan__{darray_r4_t(__Q_n_samples__)}
{

}

gms::radiolocation
::cmplx_trapezw_env_t
::cmplx_trapezw_env_t(cmplx_trapezw_env_t &&other)
:
__I_n_samples__{std::move(other.__I_n_samples__)},
__Q_n_samples__{std::move(other.__Q_n_samples__)},
__I_n_K__{      std::move(other.__I_n_K__)},
__Q_n_K__{      std::move(other.__Q_n_K__)},
__I_a__{        std::move(other.__I_a__)},
__I_l__{        std::move(other.__I_l__)},
__I_c__{        std::move(other.__I_c__)},
__I_m__{        std::move(other.__I_m__)},
__Q_a__{        std::move(other.__Q_a__)},
__Q_l__{        std::move(other.__Q_l__)},
__Q_c__{        std::move(other.__Q_c__)},
__Q_m__{        std::move(other.__Q_m__)},
__I_chan__{     std::move(other.__I_chan__)},
__Q_chan__{     std::move(other.__Q_chan__)}
{

}

gms::radiolocation
::cmplx_trapezw_env_t
::~cmplx_trapezw_env_t()
{

}

gms::radiolocation
::cmplx_trapezw_env_t &
gms::radiolocation
::cmplx_trapezw_env_t
::operator=(cmplx_trapezw_env_t &&other)
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->__I_n_samples__ = std::move(other.__I_n_samples__);
    this->__Q_n_samples__ = std::move(other.__Q_n_samples__);
    this->__I_n_K__       = std::move(other.__I_n_K__);
    this->__Q_n_K__       = std::move(other.__Q_n_K__);
    this->__I_a__         = std::move(other.__I_a__);
    this->__I_l__         = std::move(other.__I_l__);
    this->__I_c__         = std::move(other.__I_c__);
    this->__I_m__         = std::move(other.__I_m__);
    this->__Q_a__         = std::move(other.__Q_a__);
    this->__Q_l__         = std::move(other.__Q_l__);
    this->__Q_c__         = std::move(other.__Q_c__);
    this->__Q_m__         = std::move(other.__Q_m__);
    this->__I_chan__.operator=(std::move(other.__I_chan__));
    this->__Q_chan__.operator=(std::move(other.__Q_chan__));
    return (*this);
}

void 
gms::radiolocation
::cmplx_trapezw_env_t
::init_storage(const float I_filler,
               const float Q_filler)
{
#if (INIT_BY_STD_FILL) == 0
     using namespace gms::common;
	 sse_memset_unroll8x_ps(&this->__I_chan__.m_data[0],I_filler,this->__I_n_samples__);
     sse_memset_unroll8x_ps(&this->__Q_chan__.m_data[0],Q_filler,this->__Q_n_samples__);
#else 
     std::fill(this->__I_chan__.m_data,this->__I_chan__.m_data+this->__I_n_samples__,I_filler);
     std::fill(this->__Q_chan__.m_data,this->__Q_chan__.m_data+this->__Q_n_samples__,Q_filler);
#endif 
}

void 
gms::radiolocation 
::cmplx_trapezw_env_t
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
::cmplx_trapezw_env_t
::chan_I_squarew_modulated(const std::uint32_t  which_squarew)
{
         const float T{static_cast<float>(this->__I_n_samples__)};
         const float invT{1.0f/this->__I_m__};
         float sum;
         switch (which_squarew)
         {
            case 0 :
               for(std::size_t __t{0ull}; __t != this->__I_n_samples__; ++__t) 
               {
                    const float t{static_cast<float>(__t)};
                    sum = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__I_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t-k*T};
                        sum += I_sample(arg)*sin_squarew_I_sample(arg,invT);
                    }
                    this->__I_chan__.m_data[__t] = sum;
               }
            break;
            case 1 : 
               for(std::size_t __t{0ull}; __t != this->__I_n_samples__; ++__t) 
               {
                    const float t{static_cast<float>(__t)};
                    sum = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__I_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t-k*T};
                        sum += I_sample(arg)*cos_squarew_I_sample(arg,invT);
                    }
                    this->__I_chan__.m_data[__t] = sum;
               }
            break;
            default : 
                return (-1);
        }
        return (0);
}

std::int32_t  
gms::radiolocation
::cmplx_trapezw_env_t
::chan_Q_squarew_modulated(const std::uint32_t  which_squarew)
{
        const float T{static_cast<float>(this->__Q_n_samples__)};
        const float invT{1.0f/this->__Q_m__};
        float sum;
        switch (which_squarew)
        {
            case 0 :
               for(std::size_t __t{0ull}; __t != this->__Q_n_samples__; ++__t) 
               {
                    const float t{static_cast<float>(__t)};
                    sum = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__Q_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t-k*T};
                        sum += Q_sample(arg)*sin_squarew_Q_sample(arg,invT);
                    }
                    this->__Q_chan__.m_data[__t] = sum;
               }
            break;
            case 1 : 
               for(std::size_t __t{0ull}; __t != this->__Q_n_samples__; ++__t) 
               {
                    const float t{static_cast<float>(__t)};
                    sum = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__I_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t-k*T};
                        sum += Q_sample(arg)*cos_squarew_Q_sample(arg,invT);
                    }
                    this->__Q_chan__.m_data[__t] = sum;
               }
            break;
            default : 
                return (-1);
        }
        return (0);
}

std::int32_t 
gms::radiolocation
::cmplx_trapezw_env_t
::chan_I_squarew_modulated_u4x(const std::uint32_t which_squarew) 
{
        const float T{static_cast<float>(this->__I_n_samples__)};
        const std::uint32_t n_samples{static_cast<std::uint32_t>(this->__I_n_samples__)};
        const float invT{1.0f/this->__I_m__};
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
        switch (which_squarew)
        {
            case 0 :
               for(__i = 0; __i != ROUND_TO_FOUR(n_samples,4); __i += 4) 
               {
                    t__i_0 += 4.0f;
                    sum0   = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__I_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t__i_0-k*T};
                        sum0 += I_sample(arg)*sin_squarew_I_sample(arg,invT);
                    }
                    this->__I_chan__.m_data[__i+0] = sum0;
                    t__i_1 += 4.0f;
                    sum1   = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__I_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t__i_1-k*T};
                        sum1 += I_sample(arg)*sin_squarew_I_sample(arg,invT);
                    }
                    this->__I_chan__.m_data[__i+1] = sum1;
                    t__i_2 += 4.0f;
                    sum2    = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__I_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t__i_2-k*T};
                        sum2 += I_sample(arg)*sin_squarew_I_sample(arg,invT);
                    }
                    this->__I_chan__.m_data[__i+2] = sum2;
                    t__i_3 += 4.0f;
                    sum3   =  0.0f;
                    for(std::uint32_t __k{0}; __k != this->__I_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t__i_3-k*T};
                        sum3 += I_sample(arg)*sin_squarew_I_sample(arg,invT);
                    }
                    this->__I_chan__.m_data[__i+3] = sum3;
               }
               for(__j = __i; __j != n_samples; ++__j) 
               {
                    const float t{static_cast<float>(__j)};
                    sum0 = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__I_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t-k*T};
                        sum0 += I_sample(arg)*sin_squarew_I_sample(arg,invT);
                    }
                    this->__I_chan__.m_data[__j] = sum0;
               }
            break;
            case 1 : 
               for(__i = 0; __i != ROUND_TO_FOUR(n_samples,4); __i += 4) 
               {
                    t__i_0 += 4.0f;
                    sum0   = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__I_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t__i_0-k*T};
                        sum0 += I_sample(arg)*cos_squarew_I_sample(arg,invT);
                    }
                    this->__I_chan__.m_data[__i+0] = sum0;
                    t__i_1 += 4.0f;
                    sum1   = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__I_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t__i_1-k*T};
                        sum1 += I_sample(arg)*cos_squarew_I_sample(arg,invT);
                    }
                    this->__I_chan__.m_data[__i+1] = sum1;
                    t__i_2 += 4.0f;
                    sum2    = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__I_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t__i_2-k*T};
                        sum2 += I_sample(arg)*cos_squarew_I_sample(arg,invT);
                    }
                    this->__I_chan__.m_data[__i+2] = sum2;
                    t__i_3 += 4.0f;
                    sum3   =  0.0f;
                    for(std::uint32_t __k{0}; __k != this->__I_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t__i_3-k*T};
                        sum3 += I_sample(arg)*cos_squarew_I_sample(arg,invT);
                    }
                    this->__I_chan__.m_data[__i+3] = sum3;
               }
               for(__j = __i; __j != n_samples; ++__j) 
               {
                    const float t{static_cast<float>(__j)};
                    sum0 = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__I_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t-k*T};
                        sum0 += I_sample(arg)*cos_squarew_I_sample(arg,invT);
                    }
                    this->__I_chan__.m_data[__j] = sum0;
               }
            break;
            default : 
                return (-1);
        }
        return (0);
}

std::int32_t 
gms::radiolocation
::cmplx_trapezw_env_t
::chan_Q_squarew_modulated_u4x(const std::uint32_t which_squarew)
{
        const float T{static_cast<float>(this->__Q_n_samples__)};
        const std::uint32_t n_samples{static_cast<std::uint32_t>(this->__Q_n_samples__)};
        const float invT{1.0f/this->__Q_m__};
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
        switch (which_squarew)
        {
            case 0 :
               for(__i = 0; __i != ROUND_TO_FOUR(n_samples,4); __i += 4) 
               {
                    t__i_0 += 4.0f;
                    sum0   = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__Q_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t__i_0-k*T};
                        sum0 += Q_sample(arg)*sin_squarew_Q_sample(arg,invT);
                    }
                    this->__Q_chan__.m_data[__i+0] = sum0;
                    t__i_1 += 4.0f;
                    sum1   = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__Q_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t__i_1-k*T};
                        sum1 += Q_sample(arg)*sin_squarew_Q_sample(arg,invT);
                    }
                    this->__Q_chan__.m_data[__i+1] = sum1;
                    t__i_2 += 4.0f;
                    sum2    = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__Q_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t__i_2-k*T};
                        sum2 += Q_sample(arg)*sin_squarew_Q_sample(arg,invT);
                    }
                    this->__Q_chan__.m_data[__i+2] = sum2;
                    t__i_3 += 4.0f;
                    sum3   =  0.0f;
                    for(std::uint32_t __k{0}; __k != this->__Q_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t__i_3-k*T};
                        sum3 += Q_sample(arg)*sin_squarew_Q_sample(arg,invT);
                    }
                    this->__Q_chan__.m_data[__i+3] = sum3;
               }
               for(__j = __i; __j != n_samples; ++__j) 
               {
                    const float t{static_cast<float>(__j)};
                    sum0 = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__Q_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t-k*T};
                        sum0 += Q_sample(arg)*sin_squarew_Q_sample(arg,invT);
                    }
                    this->__Q_chan__.m_data[__j] = sum0;
               }
            break;
            case 1 : 
               for(__i = 0; __i != ROUND_TO_FOUR(n_samples,4); __i += 4) 
               {
                    t__i_0 += 4.0f;
                    sum0   = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__Q_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t__i_0-k*T};
                        sum0 += Q_sample(arg)*cos_squarew_Q_sample(arg,invT);
                    }
                    this->__Q_chan__.m_data[__i+0] = sum0;
                    t__i_1 += 4.0f;
                    sum1   = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__Q_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t__i_1-k*T};
                        sum1 += Q_sample(arg)*cos_squarew_Q_sample(arg,invT);
                    }
                    this->__Q_chan__.m_data[__i+1] = sum1;
                    t__i_2 += 4.0f;
                    sum2    = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__Q_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t__i_2-k*T};
                        sum2 += Q_sample(arg)*cos_squarew_Q_sample(arg,invT);
                    }
                    this->__Q_chan__.m_data[__i+2] = sum2;
                    t__i_3 += 4.0f;
                    sum3   =  0.0f;
                    for(std::uint32_t __k{0}; __k != this->__Q_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t__i_3-k*T};
                        sum3 += Q_sample(arg)*cos_squarew_Q_sample(arg,invT);
                    }
                    this->__Q_chan__.m_data[__i+3] = sum3;
               }
               for(__j = __i; __j != n_samples; ++__j) 
               {
                    const float t{static_cast<float>(__j)};
                    sum0 = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__Q_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        const float arg{t-k*T};
                        sum0 += Q_sample(arg)*cos_squarew_Q_sample(arg,invT);
                    }
                    this->__Q_chan__.m_data[__j] = sum0;
               }
            break;
            default : 
                return (-1);
        }
        return (0);
}

std::int32_t 
gms::radiolocation
::cmplx_trapezw_env_t
::chan_I_data_symbol(const float * __restrict__ I_sym_in, // size of __I_n_samples__ * __I_n_K__
                     const std::uint32_t n_T,
                     const std::uint32_t n_K)
{
      if(__builtin_expect(static_cast<std::uint32_t>(this->__I_n_samples__)!=n_T,0) || 
         __builtin_expect(this->__I_n_K__!=n_K,0)) { return (-2);}

      const float T{static_cast<float>(this->__I_n_samples__)};
      float sum;
      for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
      {
            const float t{static_cast<float>(__t)};
            sum = 0.0f;
            for(std::uint32_t __k{0}; __k != n_K; ++__k) 
            {
                const float k{static_cast<float>(__k)};
                const float arg{t-k*T};
                sum += I_sample(arg)*I_sym_in[Ix2D(__t,n_K,__k)];
            }
            this->__I_chan__.m_data[__t] = sum;
      }
      return (0);
}

std::int32_t 
gms::radiolocation 
::cmplx_trapezw_env_t
::chan_Q_data_symbol(const float * __restrict__ Q_sym_in, // size of __Q_n_samples__ * __Q_n_K__
                     const std::uint32_t n_T,
                     const std::uint32_t n_K)
{
      if(__builtin_expect(static_cast<std::uint32_t>(this->__Q_n_samples__)!=n_T,0) || 
         __builtin_expect(this->__Q_n_K__!=n_K,0)) { return (-2);}

      const float T{static_cast<float>(this->__Q_n_samples__)};
      float sum;
      for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
      {
            const float t{static_cast<float>(__t)};
            sum = 0.0f;
            for(std::uint32_t __k{0}; __k != n_K; ++__k) 
            {
                const float k{static_cast<float>(__k)};
                const float arg{t-k*T};
                sum += Q_sample(arg)*Q_sym_in[Ix2D(__t,n_K,__k)];
            }
            this->__Q_chan__.m_data[__t] = sum;
      }
      return (0);
}

std::int32_t 
gms::radiolocation
::cmplx_trapezw_env_t
::chan_I_data_symbol_u4x(const float * __restrict__ I_sym_in, // size of __I_n_samples__ * __I_n_K__
                         const std::uint32_t n_T,
                         const std::uint32_t n_K)
{
      if(__builtin_expect(static_cast<std::uint32_t>(this->__I_n_samples__)!=n_T,0) || 
         __builtin_expect(this->__I_n_K__!=n_K,0)) { return (-2);}

      const float T{static_cast<float>(this->__I_n_samples__)};
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
      for(__i = 0; __i != ROUND_TO_FOUR(n_T,4); __i += 4) 
      {
           t__i_0 += 4.0f;
           sum0   = 0.0f;
           for(std::uint32_t __k{0}; __k != n_K; ++__k) 
           {
                const float k{static_cast<float>(__k)};
                const float arg{t__i_0-k*T};
                sum0 += I_sample(arg)*I_sym_in[Ix2D(__i+0,n_K,__k)];
           }
           this->__I_chan__.m_data[__i+0] = sum0;
           t__i_1 += 4.0f;
           sum1   = 0.0f;
           for(std::uint32_t __k{0}; __k != n_K; ++__k) 
           {
                const float k{static_cast<float>(__k)};
                const float arg{t__i_1-k*T};
                sum1 += I_sample(arg)*I_sym_in[Ix2D(__i+1,n_K,__k)];
           }
           this->__I_chan__.m_data[__i+1] = sum1;
           t__i_2 += 4.0f;
           sum2    = 0.0f;
           for(std::uint32_t __k{0}; __k != n_K; ++__k) 
           {
                const float k{static_cast<float>(__k)};
                const float arg{t__i_2-k*T};
                sum2 += I_sample(arg)*I_sym_in[Ix2D(__i+2,n_K,__k)];
           }
           this->__I_chan__.m_data[__i+2] = sum2;
           t__i_3 += 4.0f;
           sum3   =  0.0f;
           for(std::uint32_t __k{0}; __k != n_K; ++__k) 
           {
                const float k{static_cast<float>(__k)};
                const float arg{t__i_3-k*T};
                sum3 += I_sample(arg)*I_sym_in[Ix2D(__i+3,n_K,__k)];
           }
           this->__I_chan__.m_data[__i+3] = sum3;
       }
       for(__j = __i; __j != n_T; ++__j) 
       {
            const float t{static_cast<float>(__j)};
            sum0 = 0.0f;
            for(std::uint32_t __k{0}; __k != n_K; ++__k) 
            {
                const float k{static_cast<float>(__k)};
                const float arg{t-k*T};
                sum0 += I_sample(arg)*I_sym_in[Ix2D(__j,n_K,__k)];
            }
            this->__I_chan__.m_data[__j] = sum0;
        }
        return (0);
}

std::int32_t 
gms::radiolocation
::cmplx_trapezw_env_t
::chan_Q_data_symbol_u4x(const float * __restrict__ Q_sym_in, // size of __Q_n_samples__ * __Q_n_K__
                         const std::uint32_t n_T,
                         const std::uint32_t n_K)
{
      if(__builtin_expect(static_cast<std::uint32_t>(this->__Q_n_samples__)!=n_T,0) || 
         __builtin_expect(this->__Q_n_K__!=n_K,0)) { return (-2);}

      const float T{static_cast<float>(this->__Q_n_samples__)};
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
      for(__i = 0; __i != ROUND_TO_FOUR(n_T,4); __i += 4) 
      {
           t__i_0 += 4.0f;
           sum0   = 0.0f;
           for(std::uint32_t __k{0}; __k != n_K; ++__k) 
           {
                const float k{static_cast<float>(__k)};
                const float arg{t__i_0-k*T};
                sum0 += Q_sample(arg)*Q_sym_in[Ix2D(__i+0,n_K,__k)];
           }
           this->__Q_chan__.m_data[__i+0] = sum0;
           t__i_1 += 4.0f;
           sum1   = 0.0f;
           for(std::uint32_t __k{0}; __k != n_K; ++__k) 
           {
                const float k{static_cast<float>(__k)};
                const float arg{t__i_1-k*T};
                sum1 += Q_sample(arg)*Q_sym_in[Ix2D(__i+1,n_K,__k)];
           }
           this->__Q_chan__.m_data[__i+1] = sum1;
           t__i_2 += 4.0f;
           sum2    = 0.0f;
           for(std::uint32_t __k{0}; __k != n_K; ++__k) 
           {
                const float k{static_cast<float>(__k)};
                const float arg{t__i_2-k*T};
                sum2 += Q_sample(arg)*Q_sym_in[Ix2D(__i+2,n_K,__k)];
           }
           this->__Q_chan__.m_data[__i+2] = sum2;
           t__i_3 += 4.0f;
           sum3   =  0.0f;
           for(std::uint32_t __k{0}; __k != n_K; ++__k) 
           {
                const float k{static_cast<float>(__k)};
                const float arg{t__i_3-k*T};
                sum3 += Q_sample(arg)*Q_sym_in[Ix2D(__i+3,n_K,__k)];
           }
           this->__Q_chan__.m_data[__i+3] = sum3;
       }
       for(__j = __i; __j != n_T; ++__j) 
       {
            const float t{static_cast<float>(__j)};
            sum0 = 0.0f;
            for(std::uint32_t __k{0}; __k != n_K; ++__k) 
            {
                const float k{static_cast<float>(__k)};
                const float arg{t-k*T};
                sum0 += Q_sample(arg)*Q_sym_in[Ix2D(__j,n_K,__k)];
            }
            this->__Q_chan__.m_data[__j] = sum0;
        }
        return (0);
}



auto
gms::radiolocation
::operator<<(std::ostream &os,
             gms::radiolocation::cmplx_trapezw_env_t &rhs)->std::ostream &
{
    std::cout << typeid(rhs).name() << "Begin: object state dump." << std::endl;
    std::cout << "__I_n_samples__: " << rhs.__I_n_samples__ << std::endl;
    std::cout << "__Q_n_samples__: " << rhs.__Q_n_samples__ << std::endl;
    std::cout << "__I_n_K__      : " << rhs.__I_n_K__       << std::endl;
    std::cout << "__Q_n_K__      :"  << rhs.__Q_n_K__       << std::endl;
    std::cout << "__I_a__        :"  << std::fixed << std::setprecision(7) << rhs.__I_a__ << std::endl;
    std::cout << "__Q_a__        :"  << std::fixed << std::setprecision(7) << rhs.__Q_a__ << std::endl;
    std::cout << "__I_m__        :"  << std::fixed << std::setprecision(7) << rhs.__I_m__ << std::endl;
    std::cout << "__Q_m__        :"  << std::fixed << std::setprecision(7) << rhs.__Q_m__ << std::endl;
    std::cout << "__I_l__        :"  << std::fixed << std::setprecision(7) << rhs.__I_l__ << std::endl;
    std::cout << "__Q_l__        :"  << std::fixed << std::setprecision(7) << rhs.__Q_l__ << std::endl;
    std::cout << "__I_c__        :"  << std::fixed << std::setprecision(7) << rhs.__I_c__ << std::endl;
    std::cout << "__Q_c__        :"  << std::fixed << std::setprecision(7) << rhs.__Q_c__ << std::endl;
    std::cout << "Channel: In-Phase(I)" << std::endl;
    for(std::size_t __i{0ull}; __i != rhs.__I_n_samples__; ++__i)
    {
         os << std::fixed << std::setprecision(7) << 
                             rhs.__I_chan__.m_data[__i] << std::endl;
    }
    std::cout <<"Channel: Quadrature(Q)" << std::endl;
    for(std::size_t __i{0ull}; __i != rhs.__Q_n_samples__; ++__i) 
    {
          os << std::fixed << std::setprecision(7) << 
                             rhs.__Q_chan__.m_data[__i] << std::endl;
    }
    std::cout << typeid(rhs).name() << "End: object state dump." << std::endl;
    return (os);
}
    






