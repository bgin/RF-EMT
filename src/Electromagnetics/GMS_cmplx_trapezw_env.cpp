
#include <fstream>
#include "GMS_cmplx_trapezw_env.h"
#include "GMS_sse_memset.h"



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
::cmplx_env_I_chan(const float * __restrict I_seq_1_0,
                   const float * __restrict I_seq_1_1,
                   const std::uint32_t  n_seq,
                   const std::uint32_t  which_seq)
{
        if(__builtin_expect(this->__I_n_K__!=n_seq,0)) {return (-2);}
        const float T{static_cast<float>(this->__I_n_samples__)};
        float sum;
         switch (which_seq)
         {
            case 0 :
               for(std::size_t __t{0ull}; __t != this->__I_n_samples__; ++__t) 
               {
                    const float t{static_cast<float>(__t)};
                    sum = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__I_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        sum += I_sample(t-k*T)*I_seq_1_0[__k];
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
                        sum += I_sample(t-k*T)*I_seq_1_1[__k];
                    }
                    this->__I_chan__.m_data[__t] = sum;
               }
            break;
            default : 
                return (-2);
        }
        return (0);
}

std::int32_t  
gms::radiolocation
::cmplx_trapezw_env_t
::cmplx_env_Q_chan(const float * __restrict Q_seq_1_0,
                   const float * __restrict Q_seq_1_1,
                   const std::uint32_t  n_seq,
                   const std::uint32_t  which_seq)
{
        if(__builtin_expect(this->__Q_n_K__!=n_seq,0)) {return (-2);}
        const float T{static_cast<float>(this->__Q_n_samples__)};
        float sum;
        std::int32_t status{9999};
        switch (which_seq)
        {
            case 0 :
               for(std::size_t __t{0ull}; __t != this->__Q_n_samples__; ++__t) 
               {
                    const float t{static_cast<float>(__t)};
                    sum = 0.0f;
                    for(std::uint32_t __k{0}; __k != this->__Q_n_K__; ++__k) 
                    {
                        const float k{static_cast<float>(__k)};
                        sum += Q_sample(t-k*T)*Q_seq_1_0[__k];
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
                        sum += Q_sample(t-k*T)*Q_seq_1_1[__k];
                    }
                    this->__Q_chan__.m_data[__t] = sum;
               }
            break;
            default : 
                return (-2);
        }
        return (0);
}

void 
gms::radiolocation
::cmplx_trapezw_env_t
::cmplx_env_I_chan_uncoded()
{
       const float T{static_cast<float>(this->__I_n_samples__)};
       float sum;
       for(std::size_t __t{0ull}; __t != this->__I_n_samples__; ++__t) 
       {
            const float t{static_cast<float>(__t)};
            sum = 0.0f;
            for(std::uint32_t __k{0}; __k != this->__I_n_K__; ++__k) 
            {
                const float k{static_cast<float>(__k)};
                sum += I_sample(t-k*T);
            }
            this->__I_chan__.m_data[__t] = sum;
        }
}

void 
gms::radiolocation
::cmplx_trapezw_env_t
::cmplx_env_Q_chan_uncoded()
{
       const float T{static_cast<float>(this->__Q_n_samples__)};
       float sum;
       for(std::size_t __t{0ull}; __t != this->__Q_n_samples__; ++__t) 
       {
            const float t{static_cast<float>(__t)};
            sum = 0.0f;
            for(std::uint32_t __k{0}; __k != this->__Q_n_K__; ++__k) 
            {
                const float k{static_cast<float>(__k)};
                sum += Q_sample(t-k*T);
            }
            this->__Q_chan__.m_data[__t] = sum;
        }
}




