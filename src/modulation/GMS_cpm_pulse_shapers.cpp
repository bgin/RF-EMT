
#include <fstream>
#include <iomanip>
#include <immintrin.h>
#include "GMS_cpm_pulse_shapers.h"

gms::radiolocation
::cpm_pulse_shapers_t
::cpm_pulse_shapers_t() noexcept(true)
:
m_nTLsamples{0ull},
m_nTsamples{0ull},
m_T{0.0f},
m_L{0.0f},
m_beta{0.0f},
m_BbT{0.0f},
m_lrec_pulse{darray_r4_t()},
m_lrc_pulse{darray_r4_t()},
m_lsrc_pulse{darray_r4_t()},
m_tfm_pulse{darray_r4_t()},
m_gmsk_pulse{darray_r4_t()}
{

}

gms::radiolocation
::cpm_pulse_shapers_t
::cpm_pulse_shapers_t(const std::size_t nTLsamples,
                      const std::size_t nTsamples,
                      const float T,
                      const float L,
                      const float beta,
                      const float BbT) noexcept(false)
:
m_nTLsamples{nTLsamples},
m_nTsamples{nTsamples},
m_T{T},
m_L{L},
m_beta{beta},
m_BbT{BbT},
m_lrec_pulse{darray_r4_t(this->m_nTLsamples)},
m_lrc_pulse{darray_r4_t(this->m_nTLsamples)},
m_lsrc_pulse{darray_r4_t(this->m_nTLsamples)},
m_tfm_pulse{darray_r4_t(this->m_nTsamples)},
m_gmsk_pulse{darray_r4_t(this->m_nTsamples)}
{

}

gms::radiolocation 
::cpm_pulse_shapers_t
::cpm_pulse_shapers_t(const cpm_pulse_shapers_t &other) noexcept(false)
:
m_nTLsamples{other.m_nTLsamples},
m_nTsamples{other.m_nTsamples},
m_T{other.m_T},
m_L{other.m_L},
m_beta{other.m_beta},
m_BbT{other.m_BbT},
m_lrec_pulse{other.m_lrec_pulse},
m_lrc_pulse{other.m_lrc_pulse},
m_lsrc_pulse{other.m_lsrc_pulse},
m_tfm_pulse{other.m_tfm_pulse},
m_gmsk_pulse{other.m_gmsk_pulse} 
{

}

gms::radiolocation 
::cpm_pulse_shapers_t
::cpm_pulse_shapers_t(cpm_pulse_shapers_t &&other) noexcept(true)
:
m_nTLsamples{other.m_nTLsamples},
m_nTsamples{other.m_nTsamples},
m_T{other.m_T},
m_L{other.m_L},
m_beta{other.m_beta},
m_BbT{other.m_BbT},
m_lrec_pulse{std::move(other.m_lrec_pulse)},
m_lrc_pulse{std::move(other.m_lrc_pulse)},
m_lsrc_pulse{std::moive(other.m_lsrc_pulse)},
m_tfm_pulse{std::move(other.m_tfm_pulse)},
m_gmsk_pulse{std::move(other.m_gmsk_pulse)} 
{

}

gms::radiolocation
::cpm_pulse_shapers_t
::~cpm_pulse_shapers_t() noexcept(true) 
{

}

gms::radiolocation
::cpm_pulse_shapers_t &
gms::radiolocation
::cpm_pulse_shapers_t
::operator=(const cpm_pulse_shapers_t &other) noexcept(false) 
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_nTLsamples = other.m_nTLsamples;
    this->m_nTsamples  = other.m_nTsamples;
    this->m_T          = other.m_T;
    this->m_L          = other.m_L;
    this->m_beta       = other.m_beta;
    this->m_BbT        = other.m_BbT; 
    this->m_lrec_pulse.operator=(other.m_lrec_pulse);
    this->m_lrc_pulse.operator=(other.m_lrc_pulse);
    this->m_lsrc_pulse.operator=(other.m_lsrc_pulse);
    this->m_tfm_pulse.operator=(other.m_tfm_pulse);
    this->m_gmsk_pulse.operator=(other.m_gmsk_pulse);
    return (*this);
}

gms::radiolocation
::cpm_pulse_shapers_t &
gms::radiolocation
::cpm_pulse_shapers_t
::operator=(cpm_pulse_shapers_t &&other) noexcept(true) 
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_nTLsamples = other.m_nTLsamples;
    this->m_nTsamples  = other.m_nTsamples;
    this->m_T          = other.m_T;
    this->m_L          = other.m_L;
    this->m_beta       = other.m_beta;
    this->m_BbT        = other.m_BbT; 
    this->m_lrec_pulse.operator=(std::move(other.m_lrec_pulse));
    this->m_lrc_pulse.operator=(std::move(other.m_lrc_pulse));
    this->m_lsrc_pulse.operator=(std::move(other.m_lsrc_pulse));
    this->m_tfm_pulse.operator=(std::move(other.m_tfm_pulse));
    this->m_gmsk_pulse.operator=(std::move(other.m_gmsk_pulse));
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



