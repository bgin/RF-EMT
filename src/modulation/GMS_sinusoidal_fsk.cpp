
#include <fstream>
#include <iomanip>
#include "GMS_sinusoidal_fsk.h"

gms::radiolocation
::sinusoidal_fsk_t
::sinusoidal_fsk_t() noexcept(true)
:
m_sfsk_nsamples{0ull},
m_I_ch_nsamples{0ull},
m_Q_ch_nsamples{0ull},
m_psfunc_nsamples{0ull},
m_k{0.0f},
m_T{0.0f},
m_Ac{0.0f},
m_fc{0.0f},
m_ph0{0.0f},
m_I_fc{0.0f},
m_Q_fc{0.0f},
m_I_ph0{0.0f},
m_Q_ph0{0.0f},
m_I_ch_bitstream{darray_r4_t()},
m_Q_ch_bitstream{darray_r4_t()},
m_psfunc{darray_r4_t()},
m_sfsk_signal{darray_r4_t()}
{

}

gms::radiolocation
::sinusoidal_fsk_t
::sinusoidal_fsk_t(const std::size_t sfsk_nsamples,
                   const std::size_t I_ch_nsamples,
                   const std::size_t Q_ch_nsamples,
                   const std::size_t psfunc_nsamples,
                   const float       k,
                   const float       T,
                   const float       Ac,
                   const float       fc,
                   const float       ph0,
                   const float       I_fc,
                   const float       Q_fc,
                   const float       I_ph0,
                   const float       Q_ph0) noexcept(false) 
:
m_sfsk_nsamples{sfsk_nsamples},
m_I_ch_nsamples{I_ch_nsamples},
m_Q_ch_nsamples{Q_ch_nsamples},
m_psfunc_nsamples{psfunc_nsamples},
m_k{k},
m_T{T},
m_Ac{Ac},
m_fc{fc},
m_ph0{ph0},
m_I_fc{I_fc},
m_Q_fc{Q_fc},
m_I_ph0{I_ph0},
m_Q_ph0{Q_ph0},
m_I_ch_bitstream{darray_r4_t(this->m_I_ch_nsamples)},
m_Q_ch_bitstream{darray_r4_t(this->m_Q_ch_nsamples)},
m_psfunc{darray_r4_t(this->m_psfunc_nsamples)},
m_sfsk_signal{darray_r4_t(this->m_sfsk_nsamples)}
{

}

gms::radiolocation 
::sinusoidal_fsk_t
::sinusoidal_fsk_t(const sinusoidal_fsk_t &other) noexcept(false)
:
m_sfsk_nsamples{other.m_sfsk_nsamples},
m_I_ch_nsamples{other.m_I_ch_nsamples},
m_Q_ch_nsamples{other.m_Q_ch_nsamples},
m_psfunc_nsamples{other.m_psfunc_nsamples},
m_k{other.m_k},
m_T{other.m_T},
m_Ac{other.m_Ac},
m_fc{other.m_fc},
m_ph0{other.m_ph0},
m_I_fc{other.m_I_fc},
m_Q_fc{other.m_Q_fc},
m_I_ph0{other.m_I_ph0},
m_Q_ph0{other.m_Q_ph0},
m_I_ch_bitstream{other.m_I_ch_bitstream},
m_Q_ch_bitstream{other.m_Q_ch_bitstream},
m_psfunc{other.m_psfunc},
m_sfsk_signal{other.m_sfsk_signal}
{

}

gms::radiolocation 
::sinusoidal_fsk_t
::sinusoidal_fsk_t(sinusoidal_fsk_t &&other) noexcept(true)
:
m_sfsk_nsamples{other.m_sfsk_nsamples},
m_I_ch_nsamples{other.m_I_ch_nsamples},
m_Q_ch_nsamples{other.m_Q_ch_nsamples},
m_psfunc_nsamples{other.m_psfunc_nsamples},
m_k{other.m_k},
m_T{other.m_T},
m_Ac{other.m_Ac},
m_fc{other.m_fc},
m_ph0{other.m_ph0},
m_I_fc{other.m_I_fc},
m_Q_fc{other.m_Q_fc},
m_I_ph0{other.m_I_ph0},
m_Q_ph0{other.m_Q_ph0},
m_I_ch_bitstream{std::move(other.m_I_ch_bitstream)},
m_Q_ch_bitstream{std::move(other.m_Q_ch_bitstream)},
m_psfunc{std::move(other.m_psfunc)},
m_sfsk_signal{std::move(other.m_sfsk_signal)}
{

}

gms::radiolocation 
::sinusoidal_fsk_t
::~sinusoidal_fsk_t() noexcept(true)
{

}

gms::radiolocation
::sinusoidal_fsk_t &
gms::radiolocation 
::sinusoidal_fsk_t
::operator=(const sinusoidal_fsk_t &other) noexcept(false) 
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_sfsk_nsamples   = other.m_sfsk_nsamples;
    this->m_I_ch_nsamples   = other.m_I_ch_nsamples;
    this->m_Q_ch_nsamples   = other.m_Q_ch_nsamples;
    this->m_psfunc_nsamples = other.m_psfunc_nsamples;
    this->m_k               = other.m_k;
    this->m_T               = other.m_T;
    this->m_Ac              = other.m_Ac;
    this->m_fc              = other.m_fc;
    this->m_ph0             = other.m_ph0;
    this->m_I_fc            = other.m_I_fc;
    this->m_Q_fc            = other.m_Q_fc;
    this->m_I_ph0           = other.m_I_ph0;
    this->m_Q_ph0           = other.m_Q_ph0;
    this->m_I_ch_bitstream.operator=(other.m_I_ch_bitstream);
    this->m_Q_ch_bitstream.operator=(other.m_Q_ch_bitstream);
    this->m_psfunc.operator=(other.m_psfunc);
    this->m_sfsk_signal.operator=(other.m_sfsk_signal);
    return (*this);
}

gms::radiolocation
::sinusoidal_fsk_t &
gms::radiolocation 
::sinusoidal_fsk_t
::operator=(sinusoidal_fsk_t &&other) noexcept(true) 
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_sfsk_nsamples   = other.m_sfsk_nsamples;
    this->m_I_ch_nsamples   = other.m_I_ch_nsamples;
    this->m_Q_ch_nsamples   = other.m_Q_ch_nsamples;
    this->m_psfunc_nsamples = other.m_psfunc_nsamples;
    this->m_k               = other.m_k;
    this->m_T               = other.m_T;
    this->m_Ac              = other.m_Ac;
    this->m_fc              = other.m_fc;
    this->m_ph0             = other.m_ph0;
    this->m_I_fc            = other.m_I_fc;
    this->m_Q_fc            = other.m_Q_fc;
    this->m_I_ph0           = other.m_I_ph0;
    this->m_Q_ph0           = other.m_Q_ph0;
    this->m_I_ch_bitstream.operator=(std::move( other.m_I_ch_bitstream));
    this->m_Q_ch_bitstream.operator=(std::move(other.m_Q_ch_bitstream));
    this->m_psfunc.operator=(std::move(other.m_psfunc));
    this->m_sfsk_signal.operator=(std::move(other.m_sfsk_signal));
    return (*this);
}

void 
gms::radiolocation 
::sinusoidal_fsk_t
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



