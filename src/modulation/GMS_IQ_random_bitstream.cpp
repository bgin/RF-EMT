
#include <fstream>
#include <iomanip>
#include "GMS_IQ_random_bitstream.h"


gms::radiolocation 
::iq_random_bitstream_t
::iq_random_bitstream_t(const float I_duration,
                        const float I_sample_rate,
                        const float I_w0,
                        const float I_ph0,
                        const float Q_duration,
                        const float Q_sample_rate,
                        const float Q_w0,
                        const float Q_ph0) noexcept(false)
:
m_I_duration{(I_duration<1.0f)?1.0f:I_duration}, // sanity check!!
m_I_sample_rate{I_sample_rate},
m_I_w0{I_w0},
m_I_ph0{I_ph0},
m_I_nsamples{static_cast<std::size_t>(this->m_I_duration*this->m_I_sample_rate)},
m_Q_duration{(Q_duration<1.0f)?1.0f:Q_duration}, // sanity check!!
m_Q_sample_rate{Q_sample_rate},
m_Q_w0{Q_w0},
m_Q_ph0{Q_ph0},
m_Q_nsamples{static_cast<std::size_t>(this->m_Q_duration*this->m_Q_sample_rate)},
m_I_bitstream{darray_r4_t(this->m_I_nsamples)},
m_Q_bitstream{darray_r4_t(this->m_Q_nsamples)}
{

}

gms::radiolocation 
::iq_random_bitstream_t
::iq_random_bitstream_t(const iq_random_bitstream_t &other) noexcept(false) 
:
m_I_duration{other.m_I_duration},
m_I_sample_rate{other.m_I_sample_rate},
m_I_w0{other.m_I_w0},
m_I_ph0{other.m_I_ph0},
m_I_nsamples{other.m_I_nsamples},
m_Q_duration{other.m_Q_duration},
m_Q_sample_rate{other.m_Q_sample_rate},
m_Q_w0{other.m_Q_w0},
m_Q_ph0{other.m_Q_ph0},
m_Q_nsamples{other.m_Q_nsamples},
m_I_bitstream{other.m_I_bitstream},
m_Q_bitstream{other.m_Q_bitstream} 
{

}

gms::radiolocation 
::iq_random_bitstream_t
::iq_random_bitstream_t(iq_random_bitstream_t &&other) noexcept(true) 
:
m_I_duration{std::move(other.m_I_duration)},
m_I_sample_rate{std::move(other.m_I_sample_rate)},
m_I_w0{std::move(other.m_I_w0)},
m_I_ph0{std::move(other.m_I_ph0)},
m_I_nsamples{std::move(other.m_I_nsamples)},
m_Q_duration{std::move(other.m_Q_duration)},
m_Q_sample_rate{std::move(other.m_Q_sample_rate)},
m_Q_w0{std::move(other.m_Q_w0)},
m_Q_ph0{std::move(other.m_Q_ph0)},
m_Q_nsamples{std::move(other.m_Q_nsamples)},
m_I_bitstream{std::move(other.m_I_bitstream)},
m_Q_bitstream{std::move(other.m_Q_bitstream)} 
{

}

gms::radiolocation
::iq_random_bitstream_t
::~iq_random_bitstream_t() noexcept(false) 
{

}

gms::radiolocation
::iq_random_bitstream_t &
gms::radiolocation
::iq_random_bitstream_t
::operator=(const iq_random_bitstream_t &other) noexcept(false) 
{
    if(__builtin_expect(this==&other,0)) {return(*this);}
    this->m_I_duration     = other.m_I_duration;
    this->m_I_sample_rate  = other.m_I_sample_rate;
    this->m_I_w0           = other.m_I_w0;
    this->m_I_ph0          = other.m_I_ph0;
    this->m_I_nsamples     = other.m_I_nsamples;
    this->m_Q_duration     = other.m_Q_duration;
    this->m_Q_sample_rate  = other.m_Q_sample_rate;
    this->m_Q_w0           = other.m_Q_w0;
    this->m_Q_ph0          = other.m_Q_ph0;
    this->m_Q_nsamples     = other.m_Q_nsamples;
    this->m_I_bitstream.operator=(other.m_I_bitstream);
    this->m_Q_bitstream.operator=(other.m_Q_bitstream);
    return (*this);
}

gms::radiolocation
::iq_random_bitstream_t &
gms::radiolocation
::iq_random_bitstream_t
::operator=(iq_random_bitstream_t &&other) noexcept(true) 
{
    if(__builtin_expect(this==&other,0)) {return(*this);}
    this->m_I_duration     = std::move(other.m_I_duration);
    this->m_I_sample_rate  = std::move(other.m_I_sample_rate);
    this->m_I_w0           = std::move(other.m_I_w0);
    this->m_I_ph0          = std::move(other.m_I_ph0);
    this->m_I_nsamples     = std::move(other.m_I_nsamples);
    this->m_Q_duration     = std::move(other.m_Q_duration);
    this->m_Q_sample_rate  = std::move(other.m_Q_sample_rate);
    this->m_Q_w0           = std::move(other.m_Q_w0);
    this->m_Q_ph0          = std::move(other.m_Q_ph0);
    this->m_Q_nsamples     = std::move(other.m_Q_nsamples);
    this->m_I_bitstream.operator=(std::move(other.m_I_bitstream));
    this->m_Q_bitstream.operator=(std::move(other.m_Q_bitstream));
    return (*this);
}

void 
gms::radiolocation 
::iq_random_bitstream_t
::create_signal_plot(const std::uint32_t n_samp,
                     const float * __restrict sig_arg,
                     const float * __restrict sig_val,
                     const std::string &header,
                     const std::string &title,
                     const bool is_sig_arg_present) noexcept(false)
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

auto 
gms::radiolocation 
::operator<<(std::ostream &os, 
             const iq_random_bitstream_t &rhs)->std::ostream & 
{
     std::cout << typeid(rhs).name() << "Begin: object state dump." << std::endl;
     std::cout << "m_I_duration    :"   << std::fixed << std::setprecision(7) << rhs.m_I_duration    << "\n";
     std::cout << "m_I_sample_rate :"   << std::fixed << std::setprecision(7) << rhs.m_I_sample_rate << "\n";
     std::cout << "m_I_w0          :"   << std::fixed << std::setprecision(7) << rhs.m_I_w0          << "\n";
     std::cout << "m_I_ph0         :"   << std::fixed << std::setprecision(7) << rhs.m_I_ph0         << "\n";
     std::cout << "m_I_nsamples    :"                                         << rhs.m_I_nsamples    << "\n";
     std::cout << "m_Q_duration    :"   << std::fixed << std::setprecision(7) << rhs.m_Q_duration    << "\n";
     std::cout << "m_Q_sample_rate :"   << std::fixed << std::setprecision(7) << rhs.m_Q_sample_rate << "\n";
     std::cout << "m_Q_w0          :"   << std::fixed << std::setprecision(7) << rhs.m_Q_w0          << "\n";
     std::cout << "m_Q_ph0         :"   << std::fixed << std::setprecision(7) << rhs.m_Q_ph0         << "\n";
     std::cout << "m_Q_nsamples    :"                                         << rhs.m_Q_nsamples    << "\n";  
     std::cout << "m_I_bitstream size and alignment info: " << "\n\n";
     rhs.m_I_bitstream.info_size_alignment();
     std::cout << "m_Q_bitstream size and alignment info: " << "\n\n";
     rhs.m_Q_bitstream.info_size_alignment(); 
     std::cout << "I-Bitstream -- uniformely distributed (PWM) samples:" << "\n";
     for(std::size_t i{0ull}; i != rhs.m_I_nsamples; ++i) 
     {
       os << std::fixed << std::setprecision(7) << rhs.m_I_bitstream.m_data[i] << "\n";
     }
     std::cout << "Q-Bitstream -- uniformely distributed (PWM) samples:" << "\n";
     for(std::size_t i{0ull}; i != rhs.m_Q_nsamples; ++i) 
     {
       os << std::fixed << std::setprecision(7) << rhs.m_Q_bitstream.m_data[i] << "\n";
     }     
     std::cout << typeid(rhs).name() << "End: object state dump." << std::endl;
     return (os);                           
}


