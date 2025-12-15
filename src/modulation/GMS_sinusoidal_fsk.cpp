
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

std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_pulse_shaping_function_u8x()
{
    using namespace gms::math;
    constexpr float C314159265358979323846264338328{3.14159265358979323846264338328f};
    constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
    const float f_psfunc_nsamples{static_cast<float>(this->m_psfunc_nsamples)};
    const float invT{1.0f/f_psfunc_nsamples};
    const float inv2T{1.0f/(2.0f*f_psfunc_nsamples)};
    constexpr std::size_t LUT_loop_idx_threshold{2257ull};
    std::size_t i,j;
    
    if(this->m_psfunc_nsamples>LUT_loop_idx_threshold)
    {
          for(i = 0ull; i != ROUND_TO_EIGHT(this->m_psfunc_nsamples,8ull); i += 8ull)
          {
             const float t_i_0{static_cast<float>(i)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_0{0.25f*ceph_sinf(C6283185307179586476925286766559*t_i_0*invT)};
             const float cos_arg_0{C314159265358979323846264338328*t_i_0*inv2T};
             const float psf_val_0{ceph_cosf(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i+0ull] = psf_val_0;
#else 
             const float sin_arg_0{0.25f*std::sin(C6283185307179586476925286766559*t_i_0*invT)};
             const float cos_arg_0{C314159265358979323846264338328*t_i_0*inv2T};
             const float psf_val_0{std::cos(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i+0ull] = psf_val_0;
#endif 

              const float t_i_1{static_cast<float>(i+1ull)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_1{0.25f*ceph_sinf(C6283185307179586476925286766559*t_i_1*invT)};
             const float cos_arg_1{C314159265358979323846264338328*t_i_1*inv2T};
             const float psf_val_1{ceph_cosf(cos_arg_1-sin_arg_1)};
             this->m_psfunc.m_data[i+1ull] = psf_val_1;
#else 
             const float sin_arg_1{0.25f*std::sin(C6283185307179586476925286766559*t_i_1*invT)};
             const float cos_arg_1{C314159265358979323846264338328*t_i_1*inv2T};
             const float psf_val_1{std::cos(cos_arg_1-sin_arg_1)};
             this->m_psfunc.m_data[i+1ull] = psf_val_1;
#endif 

              const float t_i_2{static_cast<float>(i+2ull)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_2{0.25f*ceph_sinf(C6283185307179586476925286766559*t_i_2*invT)};
             const float cos_arg_2{C314159265358979323846264338328*t_i_2*inv2T};
             const float psf_val_2{ceph_cosf(cos_arg_2-sin_arg_2)};
             this->m_psfunc.m_data[i+2ull] = psf_val_2;
#else 
             const float sin_arg_2{0.25f*std::sin(C6283185307179586476925286766559*t_i_2*invT)};
             const float cos_arg_2{C314159265358979323846264338328*t_i_2*inv2T};
             const float psf_val_2{std::cos(cos_arg_2-sin_arg_2)};
             this->m_psfunc.m_data[i+2ull] = psf_val_2;
#endif 

              const float t_i_3{static_cast<float>(i+3ull)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_3{0.25f*ceph_sinf(C6283185307179586476925286766559*t_i_3*invT)};
             const float cos_arg_3{C314159265358979323846264338328*t_i_3*inv2T};
             const float psf_val_3{ceph_cosf(cos_arg_3-sin_arg_3)};
             this->m_psfunc.m_data[i+3ull] = psf_val_3;
#else 
             const float sin_arg_3{0.25f*std::sin(C6283185307179586476925286766559*t_i_3*invT)};
             const float cos_arg_3{C314159265358979323846264338328*t_i_3*inv2T};
             const float psf_val_3{std::cos(cos_arg_3-sin_arg_3)};
             this->m_psfunc.m_data[i+3ull] = psf_val_3;
#endif 

              const float t_i_4{static_cast<float>(i+4ull)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_4{0.25f*ceph_sinf(C6283185307179586476925286766559*t_i_4*invT)};
             const float cos_arg_4{C314159265358979323846264338328*t_i_4*inv2T};
             const float psf_val_4{ceph_cosf(cos_arg_4-sin_arg_4)};
             this->m_psfunc.m_data[i+4ull] = psf_val_4;
#else 
             const float sin_arg_4{0.25f*std::sin(C6283185307179586476925286766559*t_i_4*invT)};
             const float cos_arg_4{C314159265358979323846264338328*t_i_4*inv2T};
             const float psf_val_4{std::cos(cos_arg_4-sin_arg_4)};
             this->m_psfunc.m_data[i+1ull] = psf_val_4;
#endif 

              const float t_i_5{static_cast<float>(i+5ull)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_5{0.25f*ceph_sinf(C6283185307179586476925286766559*t_i_5*invT)};
             const float cos_arg_5{C314159265358979323846264338328*t_i_5*inv2T};
             const float psf_val_5{ceph_cosf(cos_arg_5-sin_arg_5)};
             this->m_psfunc.m_data[i+5ull] = psf_val_5;
#else 
             const float sin_arg_5{0.25f*std::sin(C6283185307179586476925286766559*t_i_5*invT)};
             const float cos_arg_5{C314159265358979323846264338328*t_i_1*inv2T};
             const float psf_val_5{std::cos(cos_arg_5-sin_arg_5)};
             this->m_psfunc.m_data[i+5ull] = psf_val_5;
#endif 

             const float t_i_6{static_cast<float>(i+6ull)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_6{0.25f*ceph_sinf(C6283185307179586476925286766559*t_i_6*invT)};
             const float cos_arg_6{C314159265358979323846264338328*t_i_6*inv2T};
             const float psf_val_6{ceph_cosf(cos_arg_6-sin_arg_6)};
             this->m_psfunc.m_data[i+6ull] = psf_val_6;
#else 
             const float sin_arg_6{0.25f*std::sin(C6283185307179586476925286766559*t_i_6*invT)};
             const float cos_arg_6{C314159265358979323846264338328*t_i_6*inv2T};
             const float psf_val_6{std::cos(cos_arg_6-sin_arg_6)};
             this->m_psfunc.m_data[i+6ull] = psf_val_6;
#endif 

              const float t_i_7{static_cast<float>(i+7ull)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_7{0.25f*ceph_sinf(C6283185307179586476925286766559*t_i_7*invT)};
             const float cos_arg_7{C314159265358979323846264338328*t_i_7*inv2T};
             const float psf_val_7{ceph_cosf(cos_arg_7-sin_arg_7)};
             this->m_psfunc.m_data[i+7ull] = psf_val_7;
#else 
             const float sin_arg_7{0.25f*std::sin(C6283185307179586476925286766559*t_i_7*invT)};
             const float cos_arg_7{C314159265358979323846264338328*t_i_7*inv2T};
             const float psf_val_7{std::cos(cos_arg_7-sin_arg_7)};
             this->m_psfunc.m_data[i+7ull] = psf_val_7;
#endif 
          }

          for(j = i; j != this->m_psfunc_nsamples; ++j)  
          {
             const float t_j{static_cast<float>(j)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_0{0.25f*ceph_sinf(C6283185307179586476925286766559*t_j*invT)};
             const float cos_arg_0{C314159265358979323846264338328*t_j*inv2T};
             const float psf_val_0{ceph_cosf(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[j] = psf_val_0;
#else 
             const float sin_arg_0{0.25f*std::sin(C6283185307179586476925286766559*t_j*invT)};
             const float cos_arg_0{C314159265358979323846264338328*t_j*inv2T};
             const float psf_val_0{std::cos(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[j] = psf_val_0;
#endif 
          }
    }
    else 
    {
          for(i = 0ull; i != ROUND_TO_EIGHT(this->m_psfunc_nsamples,8ull); i += 8ull) 
          {
             const float t_i_0{gms::math::LUT_loop_indices_2257_align16[i+0ull]};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_0{0.25f*ceph_sinf(C6283185307179586476925286766559*t_i_0*invT)};
             const float cos_arg_0{C314159265358979323846264338328*t_i_0*inv2T};
             const float psf_val_0{ceph_cosf(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i+0ull] = psf_val_0;
#else 
             const float sin_arg_0{0.25f*std::sin(C6283185307179586476925286766559*t_i_0*invT)};
             const float cos_arg_0{C314159265358979323846264338328*t_i_0*inv2T};
             const float psf_val_0{std::cos(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i+0ull] = psf_val_0;
#endif 

              const float t_i_1{gms::math::LUT_loop_indices_2257_align16[i+1ull]};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_1{0.25f*ceph_sinf(C6283185307179586476925286766559*t_i_1*invT)};
             const float cos_arg_1{C314159265358979323846264338328*t_i_1*inv2T};
             const float psf_val_1{ceph_cosf(cos_arg_1-sin_arg_1)};
             this->m_psfunc.m_data[i+1ull] = psf_val_1;
#else 
             const float sin_arg_1{0.25f*std::sin(C6283185307179586476925286766559*t_i_1*invT)};
             const float cos_arg_1{C314159265358979323846264338328*t_i_1*inv2T};
             const float psf_val_1{std::cos(cos_arg_1-sin_arg_1)};
             this->m_psfunc.m_data[i+1ull] = psf_val_1;
#endif 

              const float t_i_2{gms::math::LUT_loop_indices_2257_align16[i+2ull]};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_2{0.25f*ceph_sinf(C6283185307179586476925286766559*t_i_2*invT)};
             const float cos_arg_2{C314159265358979323846264338328*t_i_2*inv2T};
             const float psf_val_2{ceph_cosf(cos_arg_2-sin_arg_2)};
             this->m_psfunc.m_data[i+2ull] = psf_val_2;
#else 
             const float sin_arg_2{0.25f*std::sin(C6283185307179586476925286766559*t_i_2*invT)};
             const float cos_arg_2{C314159265358979323846264338328*t_i_2*inv2T};
             const float psf_val_2{std::cos(cos_arg_2-sin_arg_2)};
             this->m_psfunc.m_data[i+2ull] = psf_val_2;
#endif 

              const float t_i_3{gms::math::LUT_loop_indices_2257_align16[i+3ull]};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_3{0.25f*ceph_sinf(C6283185307179586476925286766559*t_i_3*invT)};
             const float cos_arg_3{C314159265358979323846264338328*t_i_3*inv2T};
             const float psf_val_3{ceph_cosf(cos_arg_3-sin_arg_3)};
             this->m_psfunc.m_data[i+3ull] = psf_val_3;
#else 
             const float sin_arg_3{0.25f*std::sin(C6283185307179586476925286766559*t_i_3*invT)};
             const float cos_arg_3{C314159265358979323846264338328*t_i_3*inv2T};
             const float psf_val_3{std::cos(cos_arg_3-sin_arg_3)};
             this->m_psfunc.m_data[i+3ull] = psf_val_3;
#endif 

              const float t_i_4{gms::math::LUT_loop_indices_2257_align16[i+4ull]};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_4{0.25f*ceph_sinf(C6283185307179586476925286766559*t_i_4*invT)};
             const float cos_arg_4{C314159265358979323846264338328*t_i_4*inv2T};
             const float psf_val_4{ceph_cosf(cos_arg_4-sin_arg_4)};
             this->m_psfunc.m_data[i+4ull] = psf_val_4;
#else 
             const float sin_arg_4{0.25f*std::sin(C6283185307179586476925286766559*t_i_4*invT)};
             const float cos_arg_4{C314159265358979323846264338328*t_i_4*inv2T};
             const float psf_val_4{std::cos(cos_arg_4-sin_arg_4)};
             this->m_psfunc.m_data[i+1ull] = psf_val_4;
#endif 

              const float t_i_5{gms::math::LUT_loop_indices_2257_align16[i+5ull]};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_5{0.25f*ceph_sinf(C6283185307179586476925286766559*t_i_5*invT)};
             const float cos_arg_5{C314159265358979323846264338328*t_i_5*inv2T};
             const float psf_val_5{ceph_cosf(cos_arg_5-sin_arg_5)};
             this->m_psfunc.m_data[i+5ull] = psf_val_5;
#else 
             const float sin_arg_5{0.25f*std::sin(C6283185307179586476925286766559*t_i_5*invT)};
             const float cos_arg_5{C314159265358979323846264338328*t_i_1*inv2T};
             const float psf_val_5{std::cos(cos_arg_5-sin_arg_5)};
             this->m_psfunc.m_data[i+5ull] = psf_val_5;
#endif 

             const float t_i_6{gms::math::LUT_loop_indices_2257_align16[i+6ull]};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_6{0.25f*ceph_sinf(C6283185307179586476925286766559*t_i_6*invT)};
             const float cos_arg_6{C314159265358979323846264338328*t_i_6*inv2T};
             const float psf_val_6{ceph_cosf(cos_arg_6-sin_arg_6)};
             this->m_psfunc.m_data[i+6ull] = psf_val_6;
#else 
             const float sin_arg_6{0.25f*std::sin(C6283185307179586476925286766559*t_i_6*invT)};
             const float cos_arg_6{C314159265358979323846264338328*t_i_6*inv2T};
             const float psf_val_6{std::cos(cos_arg_6-sin_arg_6)};
             this->m_psfunc.m_data[i+6ull] = psf_val_6;
#endif 

              const float t_i_7{gms::math::LUT_loop_indices_2257_align16[i+7ull]};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_7{0.25f*ceph_sinf(C6283185307179586476925286766559*t_i_7*invT)};
             const float cos_arg_7{C314159265358979323846264338328*t_i_7*inv2T};
             const float psf_val_7{ceph_cosf(cos_arg_7-sin_arg_7)};
             this->m_psfunc.m_data[i+7ull] = psf_val_7;
#else 
             const float sin_arg_7{0.25f*std::sin(C6283185307179586476925286766559*t_i_7*invT)};
             const float cos_arg_7{C314159265358979323846264338328*t_i_7*inv2T};
             const float psf_val_7{std::cos(cos_arg_7-sin_arg_7)};
             this->m_psfunc.m_data[i+7ull] = psf_val_7;
#endif 
          }

          for(j = i; j != this->m_psfunc_nsamples; ++j)  
          {
             const float t_j{gms::math::LUT_loop_indices_2257_align16[j]};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_0{0.25f*ceph_sinf(C6283185307179586476925286766559*t_j*invT)};
             const float cos_arg_0{C314159265358979323846264338328*t_j*inv2T};
             const float psf_val_0{ceph_cosf(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[j] = psf_val_0;
#else 
             const float sin_arg_0{0.25f*std::sin(C6283185307179586476925286766559*t_j*invT)};
             const float cos_arg_0{C314159265358979323846264338328*t_j*inv2T};
             const float psf_val_0{std::cos(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[j] = psf_val_0;
#endif 
          }
    }

    return (0);
}

std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_pulse_shaping_function_sse_u4x(const bool do_prefetch_constants)
{
     using namespace gms::math;
     const float f_psfunc_nsamples{static_cast<float>(this->m_psfunc_nsamples)};
     const float invT{1.0f/f_psfunc_nsamples};
     const float inv2T{1.0f/(2.0f*f_psfunc_nsamples)};
     constexpr std::size_t LUT_loop_idx_threshold{2257ull};
     const __m128 vinvT{_mm_set1_ps(invT)};
     const __m128 vinv2T{_mm_set1_ps(inv2T)};
     __m128 c0_3,c4_7;
     __m128 c8_11,c12_15;
     __m128 vpi,v2pi;
     __m128 v025;
     std::size_t i;
     float j;
     if(!do_prefetch_constants)
     {
          c0_3  = _mm_setr_ps(0.0f,1.0f,2.0f,3.0f);
          c4_7  = _mm_setr_ps(4.0f,5.0f,6.0f,7.0f);
          c8_11 = _mm_setr_ps(8.0f,9.0f,10.0f,11.0f);
          c12_15= _mm_setr_ps(12.0f,13.0f,14.0f,15.0f);
          vpi   = _mm_set1_ps(3.14159265358979323846264338328f);
          v2pi  = _mm_set1_ps(6.283185307179586476925286766559f); 
          v025  = _mm_set1_ps(0.25f);
     }
     else 
     {
          __ATTR_ALIGN__(16) const float prefetched_constants[32] = {
                                        0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,
                                        8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f,
                                        3.14159265358979323846264338328f,
                                        3.14159265358979323846264338328f,
                                        3.14159265358979323846264338328f,
                                        3.14159265358979323846264338328f,
                                        6.283185307179586476925286766559f,
                                        6.283185307179586476925286766559f,
                                        6.283185307179586476925286766559f,
                                        6.283185307179586476925286766559f,
                                        0.25f,0.25f,0.25f,0.25f,0.0f,0.0f,0.0f,0.0f};
          _mm_prefetch((const char*)&prefetched_constants[0]);
          _mm_prefetch((const char*)&prefetched_constants[16]);
          c0_3  = _mm_load_ps(&prefetched_constants[0]);
          c4_7  = _mm_load_ps(&prefetched_constants[4]);
          c8_11 = _mm_load_ps(&prefetched_constants[8]);
          c12_15= _mm_load_ps(&prefetched_constants[12]);
          vpi   = _mm_load_ps(&prefetched_constants[16]);
          v2pi  = _mm_load_ps(&prefetched_constants[20]);
          v025  = _mm_load_ps(&prefetched_constants[24]);
     }

     if(this->m_psfunc_nsamples>LUT_loop_idx_threshold)
     {
          for(i = 0ull,j = 0.0f; (i+15ull) < this->m_psfunc_nsamples; i += 16ull;j += 16.0f) 
          {
               const __m128 vt_i_0{_mm_add_ps(_mm_set1_ps(j),c0_3)};
               const __m128 vsin_arg_0{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_0,vinvT))))};
               const __m128 vcos_arg_0{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_0,vinv2T))};
               const __m128 vpsf_value_0{_mm_cos_ps(_mm_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m128 vt_i_1{_mm_add_ps(_mm_set1_ps(j),c4_7)};
               const __m128 vsin_arg_1{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_1,vinvT))))};
               const __m128 vcos_arg_1{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_1,vinv2T))};
               const __m128 vpsf_value_1{_mm_cos_ps(_mm_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm_store_ps(&this->m_psfunc.m_data[i+4ull],vpsf_value_1);
               const __m128 vt_i_2{_mm_add_ps(_mm_set1_ps(j),c8_11)};
               const __m128 vsin_arg_2{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_2,vinvT))))};
               const __m128 vcos_arg_2{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_2,vinv2T))};
               const __m128 vpsf_value_2{_mm_cos_ps(_mm_sub_ps(vcos_arg_2,vsin_arg_2))};
               _mm_store_ps(&this->m_psfunc.m_data[i+8ull],vpsf_value_2);
               const __m128 vt_i_3{_mm_add_ps(_mm_set1_ps(j),c12_15)};
               const __m128 vsin_arg_3{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_3,vinvT))))};
               const __m128 vcos_arg_3{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_3,vinv2T))};
               const __m128 vpsf_value_3{_mm_cos_ps(_mm_sub_ps(vcos_arg_3,vsin_arg_3))};
               _mm_store_ps(&this->m_psfunc.m_data[i+12ull],vpsf_value_3);
          }

          for(; (i+11ull) < this->m_psfunc_nsamples; i += 12ull,j += 12.0f) 
          {
               const __m128 vt_i_0{_mm_add_ps(_mm_set1_ps(j),c0_3)};
               const __m128 vsin_arg_0{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_0,vinvT))))};
               const __m128 vcos_arg_0{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_0,vinv2T))};
               const __m128 vpsf_value_0{_mm_cos_ps(_mm_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m128 vt_i_1{_mm_add_ps(_mm_set1_ps(j),c4_7)};
               const __m128 vsin_arg_1{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_1,vinvT))))};
               const __m128 vcos_arg_1{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_1,vinv2T))};
               const __m128 vpsf_value_1{_mm_cos_ps(_mm_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm_store_ps(&this->m_psfunc.m_data[i+4ull],vpsf_value_1);
               const __m128 vt_i_2{_mm_add_ps(_mm_set1_ps(j),c8_11)};
               const __m128 vsin_arg_2{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_2,vinvT))))};
               const __m128 vcos_arg_2{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_2,vinv2T))};
               const __m128 vpsf_value_2{_mm_cos_ps(_mm_sub_ps(vcos_arg_2,vsin_arg_2))};
               _mm_store_ps(&this->m_psfunc.m_data[i+8ull],vpsf_value_2);
          }

          for(; (i+7ull) < this->m_psfunc_nsamples; i += 8ull,j += 8.0f) 
          {
               const __m128 vt_i_0{_mm_add_ps(_mm_set1_ps(j),c0_3)};
               const __m128 vsin_arg_0{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_0,vinvT))))};
               const __m128 vcos_arg_0{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_0,vinv2T))};
               const __m128 vpsf_value_0{_mm_cos_ps(_mm_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m128 vt_i_1{_mm_add_ps(_mm_set1_ps(j),c4_7)};
               const __m128 vsin_arg_1{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_1,vinvT))))};
               const __m128 vcos_arg_1{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_1,vinv2T))};
               const __m128 vpsf_value_1{_mm_cos_ps(_mm_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm_store_ps(&this->m_psfunc.m_data[i+4ull],vpsf_value_1);
          }

          for(; (i+3ull) < this->m_psfunc_nsamples; i += 4ull,j += 4.0f) 
          {
               const __m128 vt_i_0{_mm_add_ps(_mm_set1_ps(j),c0_3)};
               const __m128 vsin_arg_0{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_0,vinvT))))};
               const __m128 vcos_arg_0{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_0,vinv2T))};
               const __m128 vpsf_value_0{_mm_cos_ps(_mm_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0)
          }

          for(; (i+0ull) < this->m_psfunc_nsamples; i += 1ull,j += 1.0f)
          {
                //const float t_j{j};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_0{0.25f*ceph_sinf(6.283185307179586476925286766559f*j*invT)};
             const float cos_arg_0{3.14159265358979323846264338328f*j*inv2T};
             const float psf_val_0{ceph_cosf(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i] = psf_val_0;
#else 
             const float sin_arg_0{0.25f*std::sin(6.283185307179586476925286766559f*j*invT)};
             const float cos_arg_0{C.314159265358979323846264338328f*j*inv2T};
             const float psf_val_0{std::cos(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i] = psf_val_0;
#endif 
          }
     }
     else 
     {
           for(i = 0ull; (i+15ull) < this->m_psfunc_nsamples; i += 16ull) 
          {
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i],_MM_HINT_T0);
               const __m128 vt_i_0{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull])};
               const __m128 vsin_arg_0{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_0,vinvT))))};
               const __m128 vcos_arg_0{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_0,vinv2T))};
               const __m128 vpsf_value_0{_mm_cos_ps(_mm_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m128 vt_i_1{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull])};
               const __m128 vsin_arg_1{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_1,vinvT))))};
               const __m128 vcos_arg_1{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_1,vinv2T))};
               const __m128 vpsf_value_1{_mm_cos_ps(_mm_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm_store_ps(&this->m_psfunc.m_data[i+4ull],vpsf_value_1);
               const __m128 vt_i_2{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull])};
               const __m128 vsin_arg_2{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_2,vinvT))))};
               const __m128 vcos_arg_2{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_2,vinv2T))};
               const __m128 vpsf_value_2{_mm_cos_ps(_mm_sub_ps(vcos_arg_2,vsin_arg_2))};
               _mm_store_ps(&this->m_psfunc.m_data[i+8ull],vpsf_value_2);
               const __m128 vt_i_3{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+12ull])};
               const __m128 vsin_arg_3{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_3,vinvT))))};
               const __m128 vcos_arg_3{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_3,vinv2T))};
               const __m128 vpsf_value_3{_mm_cos_ps(_mm_sub_ps(vcos_arg_3,vsin_arg_3))};
               _mm_store_ps(&this->m_psfunc.m_data[i+12ull],vpsf_value_3);
          }

          for(; (i+11ull) < this->m_psfunc_nsamples; i += 12ull) 
          {
               const __m128 vt_i_0{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull])};
               const __m128 vsin_arg_0{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_0,vinvT))))};
               const __m128 vcos_arg_0{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_0,vinv2T))};
               const __m128 vpsf_value_0{_mm_cos_ps(_mm_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m128 vt_i_1{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull])};
               const __m128 vsin_arg_1{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_1,vinvT))))};
               const __m128 vcos_arg_1{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_1,vinv2T))};
               const __m128 vpsf_value_1{_mm_cos_ps(_mm_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm_store_ps(&this->m_psfunc.m_data[i+4ull],vpsf_value_1);
               const __m128 vt_i_2{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull])};
               const __m128 vsin_arg_2{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_2,vinvT))))};
               const __m128 vcos_arg_2{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_2,vinv2T))};
               const __m128 vpsf_value_2{_mm_cos_ps(_mm_sub_ps(vcos_arg_2,vsin_arg_2))};
               _mm_store_ps(&this->m_psfunc.m_data[i+8ull],vpsf_value_2);
          }

          for(; (i+7ull) < this->m_psfunc_nsamples; i += 8ull) 
          {
               const __m128 vt_i_0{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull])};
               const __m128 vsin_arg_0{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_0,vinvT))))};
               const __m128 vcos_arg_0{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_0,vinv2T))};
               const __m128 vpsf_value_0{_mm_cos_ps(_mm_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m128 vt_i_1{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull])};
               const __m128 vsin_arg_1{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_1,vinvT))))};
               const __m128 vcos_arg_1{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_1,vinv2T))};
               const __m128 vpsf_value_1{_mm_cos_ps(_mm_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm_store_ps(&this->m_psfunc.m_data[i+4ull],vpsf_value_1);
          }

          for(; (i+3ull) < this->m_psfunc_nsamples; i += 4ull) 
          {
               const __m128 vt_i_0{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull])};
               const __m128 vsin_arg_0{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(v2pi,_mm_mul_ps(vt_i_0,vinvT))))};
               const __m128 vcos_arg_0{_mm_mul_ps(vpi,_mm_mul_ps(vt_i_0,vinv2T))};
               const __m128 vpsf_value_0{_mm_cos_ps(_mm_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0)
          }

          for(; (i+0ull) < this->m_psfunc_nsamples; i += 1ull)
          {
             const float t_i{static_cast<float>(gms::math::LUT_loop_indices_2257_align16[i])};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_0{0.25f*ceph_sinf(6.283185307179586476925286766559f*t_i*invT)};
             const float cos_arg_0{3.14159265358979323846264338328f*t_i*inv2T};
             const float psf_val_0{ceph_cosf(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i] = psf_val_0;
#else 
             const float sin_arg_0{0.25f*std::sin(6.283185307179586476925286766559f*t_i*invT)};
             const float cos_arg_0{C.314159265358979323846264338328f*t_i*inv2T};
             const float psf_val_0{std::cos(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i] = psf_val_0;
#endif 
          }
     }

     return (0);

}

std::int32_t
gms::radiolocation
::sinusoidal_fsk_t
::generate_pulse_shaping_function_avx_u4x(const bool do_prefetch_constants)
{
     using namespace gms::math;
     const float f_psfunc_nsamples{static_cast<float>(this->m_psfunc_nsamples)};
     const float invT{1.0f/f_psfunc_nsamples};
     const float inv2T{1.0f/(2.0f*f_psfunc_nsamples)};
     constexpr std::size_t LUT_loop_idx_threshold{2257ull};
     const __m256 vinvT{_mm256_set1_ps(invT)};
     const __m256 vinv2T{_mm256_set1_ps(inv2T)};
     __m256 c0_7,c8_15;
     __m256 c16_23,c24_31;
     __m256 vpi,v2pi;
     __m256 v025;
     std::size_t i;
     float j;
     if(!do_prefetch_constants)
     {
          c0_7   = _mm256_setr_ps(0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f);
          c8_15  = _mm256_setr_ps(8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f);
          c16_23 = _mm256_setr_ps(16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f);
          c24_31 = _mm256_setr_ps(24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f);
          vpi    = _mm256_set1_ps(3.14159265358979323846264338328f);
          v2pi   = _mm256_set1_ps(6.283185307179586476925286766559f); 
          v025   = _mm256_set1_ps(0.25f);
     }
     else 
     {
          __ATTR_ALIGN__(32) const float prefetched_constants[64] = {
                                    0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,          
                                    8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f,
                                    16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f,
                                    24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f,
                                    3.14159265358979323846264338328f,3.14159265358979323846264338328f,
                                    3.14159265358979323846264338328f,3.14159265358979323846264338328f,
                                    3.14159265358979323846264338328f,3.14159265358979323846264338328f,
                                    3.14159265358979323846264338328f,3.14159265358979323846264338328f,
                                    6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                    6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                    6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                    6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                    0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,
                                    0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
          _mm_prefetch((const char*)&prefetched_constants[0]);
          _mm_prefetch((const char*)&prefetched_constants[16]);
          _mm_prefetch((const char*)&prefetched_constants[32]);
          _mm_prefetch((const char*)&prefetched_constants[48]);
          c0_7   = _mm256_load_ps(&prefetched_constants[0]);
          c8_15  = _mm256_load_ps(&prefetched_constants[8]);
          c16_23 = _mm256_load_ps(&prefetched_constants[16]);
          c24_31 = _mm256_load_ps(&prefetched_constants[24]);
          vpi    = _mm256_load_ps(&prefetched_constants[32]);
          v2pi   = _mm256_load_ps(&prefetched_constants[40]);
          v025   = _mm256_load_ps(&prefetched_constants[48]);
     }

     if(this->m_psfunc_nsamples>LUT_loop_idx_threshold)
     {
          for(i = 0ull,j = 0.0f; (i+31ull) < this->m_psfunc_nsamples; i += 32ull,j += 32.0f) 
          {
               const __m256 vt_i_0{_mm256_add_ps(_mm256_set1_ps(j),c0_7)};
               const __m256 vsin_arg_0{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(v2pi,_mm256_mul_ps(vt_i_0,vinvT))))};
               const __m256 vcos_arg_0{_mm256_mul_ps(vpi,_mm256_mul_ps(vt_i_0,vinv2T))};
               const __m256 vpsf_value_0{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m256 vt_i_1{_mm256_add_ps(_mm256_set1_ps(j),c8_15)};
               const __m256 vsin_arg_1{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(v2pi,_mm256_mul_ps(vt_i_1,vinvT))))};
               const __m256 vcos_arg_1{_mm256_mul_ps(vpi,_mm256_mul_ps(vt_i_1,vinv2T))};
               const __m256 vpsf_value_1{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+8ull],vpsf_value_1);
               const __m256 vt_i_2{_mm256_add_ps(_mm256_set1_ps(j),c16_23)};
               const __m256 vsin_arg_2{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(v2pi,_mm256_mul_ps(vt_i_2,vinvT))))};
               const __m256 vcos_arg_2{_mm256_mul_ps(vpi,_mm256_mul_ps(vt_i_2,vinv2T))};
               const __m256 vpsf_value_2{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_2,vsin_arg_2))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+16ull],vpsf_value_2);
               const __m256 vt_i_3{_mm256_add_ps(_mm256_set1_ps(j),c16_23)};
               const __m256 vsin_arg_3{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(v2pi,_mm256_mul_ps(vt_i_3,vinvT))))};
               const __m256 vcos_arg_3{_mm256_mul_ps(vpi,_mm256_mul_ps(vt_i_3,vinv2T))};
               const __m256 vpsf_value_3{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_3,vsin_arg_3))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+24ull],vpsf_value_3);
          }

          for(; (i+15ull) < this->m_psfunc_nsamples; i += 16ull,j += 16.0f) 
          {
                const __m256 vt_i_0{_mm256_add_ps(_mm256_set1_ps(j),c0_7)};
               const __m256 vsin_arg_0{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(v2pi,_mm256_mul_ps(vt_i_0,vinvT))))};
               const __m256 vcos_arg_0{_mm256_mul_ps(vpi,_mm256_mul_ps(vt_i_0,vinv2T))};
               const __m256 vpsf_value_0{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m256 vt_i_1{_mm256_add_ps(_mm256_set1_ps(j),c8_15)};
               const __m256 vsin_arg_1{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(v2pi,_mm256_mul_ps(vt_i_1,vinvT))))};
               const __m256 vcos_arg_1{_mm256_mul_ps(vpi,_mm256_mul_ps(vt_i_1,vinv2T))};
               const __m256 vpsf_value_1{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+8ull],vpsf_value_1);
          }

          for(; (i+7ull) < this->m_psfunc_nsamples; i += 8ull,j += 8.0f) 
          {
               const __m256 vt_i_0{_mm256_add_ps(_mm256_set1_ps(j),c0_7)};
               const __m256 vsin_arg_0{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(v2pi,_mm256_mul_ps(vt_i_0,vinvT))))};
               const __m256 vcos_arg_0{_mm256_mul_ps(vpi,_mm256_mul_ps(vt_i_0,vinv2T))};
               const __m256 vpsf_value_0{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0)
          }

          for(; (i+0ull) < this->m_psfunc_nsamples; i += 1ull,j += 1.0f)
          {
                //const float t_j{j};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_0{0.25f*ceph_sinf(6.283185307179586476925286766559f*j*invT)};
             const float cos_arg_0{3.14159265358979323846264338328f*j*inv2T};
             const float psf_val_0{ceph_cosf(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i] = psf_val_0;
#else 
             const float sin_arg_0{0.25f*std::sin(6.283185307179586476925286766559f*j*invT)};
             const float cos_arg_0{C.314159265358979323846264338328f*j*inv2T};
             const float psf_val_0{std::cos(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i] = psf_val_0;
#endif 
          }
     }
     else 
     {
           for(i = 0ull; (i+31ull) < this->m_psfunc_nsamples; i += 32ull) 
           {
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+0ull],_MM_HINT_T0);
               const __m256 vt_i_0{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull])};
               const __m256 vsin_arg_0{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(v2pi,_mm256_mul_ps(vt_i_0,vinvT))))};
               const __m256 vcos_arg_0{_mm256_mul_ps(vpi,_mm256_mul_ps(vt_i_0,vinv2T))};
               const __m256 vpsf_value_0{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m256 vt_i_1{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull])};
               const __m256 vsin_arg_1{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(v2pi,_mm256_mul_ps(vt_i_1,vinvT))))};
               const __m256 vcos_arg_1{_mm256_mul_ps(vpi,_mm256_mul_ps(vt_i_1,vinv2T))};
               const __m128 vpsf_value_1{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm_store_ps(&this->m_psfunc.m_data[i+8ull],vpsf_value_1);
              _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+16ull],_MM_HINT_T0);
               const __m256 vt_i_2{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+16ull])};
               const __m256 vsin_arg_2{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(v2pi,_mm256_mul_ps(vt_i_2,vinvT))))};
               const __m256 vcos_arg_2{_mm256_mul_ps(vpi,_mm256_mul_ps(vt_i_2,vinv2T))};
               const __m256 vpsf_value_2{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_2,vsin_arg_2))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+16ull],vpsf_value_2);
               const __m256 vt_i_3{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+24ull])};
               const __m256 vsin_arg_3{_mm_256mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(v2pi,_mm256_mul_ps(vt_i_3,vinvT))))};
               const __m256 vcos_arg_3{_mm256_mul_ps(vpi,_mm256_mul_ps(vt_i_3,vinv2T))};
               const __m256 vpsf_value_3{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_3,vsin_arg_3))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+24ull],vpsf_value_3);
          }

          for(; (i+15ull) < this->m_psfunc_nsamples; i += 16ull) 
          {
                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+0ull],_MM_HINT_T0);
               const __m256 vt_i_0{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull])};
               const __m256 vsin_arg_0{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(v2pi,_mm256_mul_ps(vt_i_0,vinvT))))};
               const __m256 vcos_arg_0{_mm256_mul_ps(vpi,_mm256_mul_ps(vt_i_0,vinv2T))};
               const __m256 vpsf_value_0{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m256 vt_i_1{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull])};
               const __m256 vsin_arg_1{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(v2pi,_mm256_mul_ps(vt_i_1,vinvT))))};
               const __m256 vcos_arg_1{_mm256_mul_ps(vpi,_mm256_mul_ps(vt_i_1,vinv2T))};
               const __m128 vpsf_value_1{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm_store_ps(&this->m_psfunc.m_data[i+8ull],vpsf_value_1);
          }

          for(; (i+7ull) < this->m_psfunc_nsamples; i += 8ull) 
          {
               const __m256 vt_i_0{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull])};
               const __m256 vsin_arg_0{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(v2pi,_mm256_mul_ps(vt_i_0,vinvT))))};
               const __m256 vcos_arg_0{_mm256_mul_ps(vpi,_mm256_mul_ps(vt_i_0,vinv2T))};
               const __m256 vpsf_value_0{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
          }

          for(; (i+0ull) < this->m_psfunc_nsamples; i += 1ull)
          {
             const float t_i{static_cast<float>(gms::math::LUT_loop_indices_2257_align16[i])};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_0{0.25f*ceph_sinf(6.283185307179586476925286766559f*t_i*invT)};
             const float cos_arg_0{3.14159265358979323846264338328f*t_i*inv2T};
             const float psf_val_0{ceph_cosf(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i] = psf_val_0;
#else 
             const float sin_arg_0{0.25f*std::sin(6.283185307179586476925286766559f*t_i*invT)};
             const float cos_arg_0{C.314159265358979323846264338328f*t_i*inv2T};
             const float psf_val_0{std::cos(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i] = psf_val_0;
#endif 
          }
     }

     return (0);
}





