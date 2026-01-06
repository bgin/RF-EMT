
#include <fstream>
#include <iomanip>
#include "GMS_sinusoidal_fsk.h"
#include "GMS_sse_memcpy.h"
#include "GMS_avx_memcpy.h"
#include "GMS_avx512_memcpy.h"

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
             const float sin_arg_0{0.25f*ceph_sinf(this->m_k*C6283185307179586476925286766559*t_i_0*invT)};
             const float cos_arg_0{this->m_k*C314159265358979323846264338328*t_i_0*inv2T};
             const float psf_val_0{ceph_cosf(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i+0ull] = psf_val_0;
#else 
             const float sin_arg_0{0.25f*std::sin(this->m_k*C6283185307179586476925286766559*t_i_0*invT)};
             const float cos_arg_0{this->m_k*C314159265358979323846264338328*t_i_0*inv2T};
             const float psf_val_0{std::cos(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i+0ull] = psf_val_0;
#endif 

              const float t_i_1{static_cast<float>(i+1ull)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_1{0.25f*ceph_sinf(this->m_k*C6283185307179586476925286766559*t_i_1*invT)};
             const float cos_arg_1{this->m_k*C314159265358979323846264338328*t_i_1*inv2T};
             const float psf_val_1{ceph_cosf(cos_arg_1-sin_arg_1)};
             this->m_psfunc.m_data[i+1ull] = psf_val_1;
#else 
             const float sin_arg_1{0.25f*std::sin(this->m_k*C6283185307179586476925286766559*t_i_1*invT)};
             const float cos_arg_1{this->m_k*C314159265358979323846264338328*t_i_1*inv2T};
             const float psf_val_1{std::cos(cos_arg_1-sin_arg_1)};
             this->m_psfunc.m_data[i+1ull] = psf_val_1;
#endif 

              const float t_i_2{static_cast<float>(i+2ull)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_2{0.25f*ceph_sinf(this->m_k*C6283185307179586476925286766559*t_i_2*invT)};
             const float cos_arg_2{this->m_k*C314159265358979323846264338328*t_i_2*inv2T};
             const float psf_val_2{ceph_cosf(cos_arg_2-sin_arg_2)};
             this->m_psfunc.m_data[i+2ull] = psf_val_2;
#else 
             const float sin_arg_2{0.25f*std::sin(this->m_k*C6283185307179586476925286766559*t_i_2*invT)};
             const float cos_arg_2{this->m_k*C314159265358979323846264338328*t_i_2*inv2T};
             const float psf_val_2{std::cos(cos_arg_2-sin_arg_2)};
             this->m_psfunc.m_data[i+2ull] = psf_val_2;
#endif 

              const float t_i_3{static_cast<float>(i+3ull)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_3{0.25f*ceph_sinf(this->m_k*C6283185307179586476925286766559*t_i_3*invT)};
             const float cos_arg_3{this->m_k*C314159265358979323846264338328*t_i_3*inv2T};
             const float psf_val_3{ceph_cosf(cos_arg_3-sin_arg_3)};
             this->m_psfunc.m_data[i+3ull] = psf_val_3;
#else 
             const float sin_arg_3{0.25f*std::sin(this->m_k*C6283185307179586476925286766559*t_i_3*invT)};
             const float cos_arg_3{this->m_k*C314159265358979323846264338328*t_i_3*inv2T};
             const float psf_val_3{std::cos(cos_arg_3-sin_arg_3)};
             this->m_psfunc.m_data[i+3ull] = psf_val_3;
#endif 

              const float t_i_4{static_cast<float>(i+4ull)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_4{0.25f*ceph_sinf(this->m_k*C6283185307179586476925286766559*t_i_4*invT)};
             const float cos_arg_4{this->m_k*C314159265358979323846264338328*t_i_4*inv2T};
             const float psf_val_4{ceph_cosf(cos_arg_4-sin_arg_4)};
             this->m_psfunc.m_data[i+4ull] = psf_val_4;
#else 
             const float sin_arg_4{0.25f*std::sin(this->m_k*C6283185307179586476925286766559*t_i_4*invT)};
             const float cos_arg_4{this->m_k*C314159265358979323846264338328*t_i_4*inv2T};
             const float psf_val_4{std::cos(cos_arg_4-sin_arg_4)};
             this->m_psfunc.m_data[i+1ull] = psf_val_4;
#endif 

              const float t_i_5{static_cast<float>(i+5ull)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_5{0.25f*ceph_sinf(this->m_k*C6283185307179586476925286766559*t_i_5*invT)};
             const float cos_arg_5{this->m_k*C314159265358979323846264338328*t_i_5*inv2T};
             const float psf_val_5{ceph_cosf(cos_arg_5-sin_arg_5)};
             this->m_psfunc.m_data[i+5ull] = psf_val_5;
#else 
             const float sin_arg_5{0.25f*std::sin(this->m_k*C6283185307179586476925286766559*t_i_5*invT)};
             const float cos_arg_5{this->m_k*C314159265358979323846264338328*t_i_1*inv2T};
             const float psf_val_5{std::cos(cos_arg_5-sin_arg_5)};
             this->m_psfunc.m_data[i+5ull] = psf_val_5;
#endif 

             const float t_i_6{static_cast<float>(i+6ull)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_6{0.25f*ceph_sinf(this->m_k*C6283185307179586476925286766559*t_i_6*invT)};
             const float cos_arg_6{this->m_k*C314159265358979323846264338328*t_i_6*inv2T};
             const float psf_val_6{ceph_cosf(cos_arg_6-sin_arg_6)};
             this->m_psfunc.m_data[i+6ull] = psf_val_6;
#else 
             const float sin_arg_6{0.25f*std::sin(this->m_k*C6283185307179586476925286766559*t_i_6*invT)};
             const float cos_arg_6{this->m_k*C314159265358979323846264338328*t_i_6*inv2T};
             const float psf_val_6{std::cos(cos_arg_6-sin_arg_6)};
             this->m_psfunc.m_data[i+6ull] = psf_val_6;
#endif 

              const float t_i_7{static_cast<float>(i+7ull)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_7{0.25f*ceph_sinf(this->m_k*C6283185307179586476925286766559*t_i_7*invT)};
             const float cos_arg_7{this->m_k*C314159265358979323846264338328*t_i_7*inv2T};
             const float psf_val_7{ceph_cosf(cos_arg_7-sin_arg_7)};
             this->m_psfunc.m_data[i+7ull] = psf_val_7;
#else 
             const float sin_arg_7{0.25f*std::sin(this->m_k*C6283185307179586476925286766559*t_i_7*invT)};
             const float cos_arg_7{this->m_k*C314159265358979323846264338328*t_i_7*inv2T};
             const float psf_val_7{std::cos(cos_arg_7-sin_arg_7)};
             this->m_psfunc.m_data[i+7ull] = psf_val_7;
#endif 
          }

          for(j = i; j != this->m_psfunc_nsamples; ++j)  
          {
             const float t_j{static_cast<float>(j)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_0{0.25f*ceph_sinf(this->m_k*C6283185307179586476925286766559*t_j*invT)};
             const float cos_arg_0{this->m_k*C314159265358979323846264338328*t_j*inv2T};
             const float psf_val_0{ceph_cosf(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[j] = psf_val_0;
#else 
             const float sin_arg_0{0.25f*std::sin(this->m_k*C6283185307179586476925286766559*t_j*invT)};
             const float cos_arg_0{this->m_k*C314159265358979323846264338328*t_j*inv2T};
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
             const float sin_arg_0{0.25f*ceph_sinf(this->m_k*C6283185307179586476925286766559*t_i_0*invT)};
             const float cos_arg_0{this->m_k*C314159265358979323846264338328*t_i_0*inv2T};
             const float psf_val_0{ceph_cosf(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i+0ull] = psf_val_0;
#else 
             const float sin_arg_0{0.25f*std::sin(this->m_k*C6283185307179586476925286766559*t_i_0*invT)};
             const float cos_arg_0{this->m_k*C314159265358979323846264338328*t_i_0*inv2T};
             const float psf_val_0{std::cos(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i+0ull] = psf_val_0;
#endif 

              const float t_i_1{gms::math::LUT_loop_indices_2257_align16[i+1ull]};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_1{0.25f*ceph_sinf(this->m_k*C6283185307179586476925286766559*t_i_1*invT)};
             const float cos_arg_1{this->m_k*C314159265358979323846264338328*t_i_1*inv2T};
             const float psf_val_1{ceph_cosf(cos_arg_1-sin_arg_1)};
             this->m_psfunc.m_data[i+1ull] = psf_val_1;
#else 
             const float sin_arg_1{0.25f*std::sin(this->m_k*C6283185307179586476925286766559*t_i_1*invT)};
             const float cos_arg_1{this->m_k*C314159265358979323846264338328*t_i_1*inv2T};
             const float psf_val_1{std::cos(cos_arg_1-sin_arg_1)};
             this->m_psfunc.m_data[i+1ull] = psf_val_1;
#endif 

              const float t_i_2{gms::math::LUT_loop_indices_2257_align16[i+2ull]};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_2{0.25f*ceph_sinf(this->m_k*C6283185307179586476925286766559*t_i_2*invT)};
             const float cos_arg_2{this->m_k*C314159265358979323846264338328*t_i_2*inv2T};
             const float psf_val_2{ceph_cosf(cos_arg_2-sin_arg_2)};
             this->m_psfunc.m_data[i+2ull] = psf_val_2;
#else 
             const float sin_arg_2{0.25f*std::sin(this->m_k*C6283185307179586476925286766559*t_i_2*invT)};
             const float cos_arg_2{this->m_k*C314159265358979323846264338328*t_i_2*inv2T};
             const float psf_val_2{std::cos(cos_arg_2-sin_arg_2)};
             this->m_psfunc.m_data[i+2ull] = psf_val_2;
#endif 

              const float t_i_3{gms::math::LUT_loop_indices_2257_align16[i+3ull]};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_3{0.25f*ceph_sinf(this->m_k*C6283185307179586476925286766559*t_i_3*invT)};
             const float cos_arg_3{this->m_k*C314159265358979323846264338328*t_i_3*inv2T};
             const float psf_val_3{ceph_cosf(cos_arg_3-sin_arg_3)};
             this->m_psfunc.m_data[i+3ull] = psf_val_3;
#else 
             const float sin_arg_3{0.25f*std::sin(this->m_k*C6283185307179586476925286766559*t_i_3*invT)};
             const float cos_arg_3{this->m_k*C314159265358979323846264338328*t_i_3*inv2T};
             const float psf_val_3{std::cos(cos_arg_3-sin_arg_3)};
             this->m_psfunc.m_data[i+3ull] = psf_val_3;
#endif 

              const float t_i_4{gms::math::LUT_loop_indices_2257_align16[i+4ull]};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_4{0.25f*ceph_sinf(this->m_k*C6283185307179586476925286766559*t_i_4*invT)};
             const float cos_arg_4{this->m_k*C314159265358979323846264338328*t_i_4*inv2T};
             const float psf_val_4{ceph_cosf(cos_arg_4-sin_arg_4)};
             this->m_psfunc.m_data[i+4ull] = psf_val_4;
#else 
             const float sin_arg_4{0.25f*std::sin(this->m_k*C6283185307179586476925286766559*t_i_4*invT)};
             const float cos_arg_4{this->m_k*C314159265358979323846264338328*t_i_4*inv2T};
             const float psf_val_4{std::cos(cos_arg_4-sin_arg_4)};
             this->m_psfunc.m_data[i+1ull] = psf_val_4;
#endif 

              const float t_i_5{gms::math::LUT_loop_indices_2257_align16[i+5ull]};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_5{0.25f*ceph_sinf(this->m_k*C6283185307179586476925286766559*t_i_5*invT)};
             const float cos_arg_5{this->m_k*C314159265358979323846264338328*t_i_5*inv2T};
             const float psf_val_5{ceph_cosf(cos_arg_5-sin_arg_5)};
             this->m_psfunc.m_data[i+5ull] = psf_val_5;
#else 
             const float sin_arg_5{0.25f*std::sin(this->m_k*C6283185307179586476925286766559*t_i_5*invT)};
             const float cos_arg_5{this->m_k*C314159265358979323846264338328*t_i_1*inv2T};
             const float psf_val_5{std::cos(cos_arg_5-sin_arg_5)};
             this->m_psfunc.m_data[i+5ull] = psf_val_5;
#endif 

             const float t_i_6{gms::math::LUT_loop_indices_2257_align16[i+6ull]};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_6{0.25f*ceph_sinf(this->m_k*C6283185307179586476925286766559*t_i_6*invT)};
             const float cos_arg_6{this->m_k*C314159265358979323846264338328*t_i_6*inv2T};
             const float psf_val_6{ceph_cosf(cos_arg_6-sin_arg_6)};
             this->m_psfunc.m_data[i+6ull] = psf_val_6;
#else 
             const float sin_arg_6{0.25f*std::sin(this->m_k*C6283185307179586476925286766559*t_i_6*invT)};
             const float cos_arg_6{this->m_k*C314159265358979323846264338328*t_i_6*inv2T};
             const float psf_val_6{std::cos(cos_arg_6-sin_arg_6)};
             this->m_psfunc.m_data[i+6ull] = psf_val_6;
#endif 

              const float t_i_7{gms::math::LUT_loop_indices_2257_align16[i+7ull]};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_7{0.25f*ceph_sinf(this->m_k*C6283185307179586476925286766559*t_i_7*invT)};
             const float cos_arg_7{this->m_k*C314159265358979323846264338328*t_i_7*inv2T};
             const float psf_val_7{ceph_cosf(cos_arg_7-sin_arg_7)};
             this->m_psfunc.m_data[i+7ull] = psf_val_7;
#else 
             const float sin_arg_7{0.25f*std::sin(this->m_k*C6283185307179586476925286766559*t_i_7*invT)};
             const float cos_arg_7{this->m_k*C314159265358979323846264338328*t_i_7*inv2T};
             const float psf_val_7{std::cos(cos_arg_7-sin_arg_7)};
             this->m_psfunc.m_data[i+7ull] = psf_val_7;
#endif 
          }

          for(j = i; j != this->m_psfunc_nsamples; ++j)  
          {
             const float t_j{gms::math::LUT_loop_indices_2257_align16[j]};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_0{0.25f*ceph_sinf(this->m_k*C6283185307179586476925286766559*t_j*invT)};
             const float cos_arg_0{this->m_k*C314159265358979323846264338328*t_j*inv2T};
             const float psf_val_0{ceph_cosf(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[j] = psf_val_0;
#else 
             const float sin_arg_0{0.25f*std::sin(this->m_k*C6283185307179586476925286766559*t_j*invT)};
             const float cos_arg_0{this->m_k*C314159265358979323846264338328*t_j*inv2T};
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
     const __m128 vk{_mm_set1_ps(this->m_k)};
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
          _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
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
          for(i = 0ull,j = 0.0f; (i+15ull) < this->m_psfunc_nsamples; i += 16ull,j += 16.0f) 
          {
               const __m128 vt_i_0{_mm_add_ps(_mm_set1_ps(j),c0_3)};
               const __m128 vsin_arg_0{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_0,vinvT))))};
               const __m128 vcos_arg_0{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_0,vinv2T))};
               const __m128 vpsf_value_0{_mm_cos_ps(_mm_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m128 vt_i_1{_mm_add_ps(_mm_set1_ps(j),c4_7)};
               const __m128 vsin_arg_1{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_1,vinvT))))};
               const __m128 vcos_arg_1{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_1,vinv2T))};
               const __m128 vpsf_value_1{_mm_cos_ps(_mm_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm_store_ps(&this->m_psfunc.m_data[i+4ull],vpsf_value_1);
               const __m128 vt_i_2{_mm_add_ps(_mm_set1_ps(j),c8_11)};
               const __m128 vsin_arg_2{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_2,vinvT))))};
               const __m128 vcos_arg_2{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_2,vinv2T))};
               const __m128 vpsf_value_2{_mm_cos_ps(_mm_sub_ps(vcos_arg_2,vsin_arg_2))};
               _mm_store_ps(&this->m_psfunc.m_data[i+8ull],vpsf_value_2);
               const __m128 vt_i_3{_mm_add_ps(_mm_set1_ps(j),c12_15)};
               const __m128 vsin_arg_3{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_3,vinvT))))};
               const __m128 vcos_arg_3{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_3,vinv2T))};
               const __m128 vpsf_value_3{_mm_cos_ps(_mm_sub_ps(vcos_arg_3,vsin_arg_3))};
               _mm_store_ps(&this->m_psfunc.m_data[i+12ull],vpsf_value_3);
          }

          for(; (i+11ull) < this->m_psfunc_nsamples; i += 12ull,j += 12.0f) 
          {
               const __m128 vt_i_0{_mm_add_ps(_mm_set1_ps(j),c0_3)};
               const __m128 vsin_arg_0{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_0,vinvT))))};
               const __m128 vcos_arg_0{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_0,vinv2T))};
               const __m128 vpsf_value_0{_mm_cos_ps(_mm_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m128 vt_i_1{_mm_add_ps(_mm_set1_ps(j),c4_7)};
               const __m128 vsin_arg_1{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_1,vinvT))))};
               const __m128 vcos_arg_1{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_1,vinv2T))};
               const __m128 vpsf_value_1{_mm_cos_ps(_mm_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm_store_ps(&this->m_psfunc.m_data[i+4ull],vpsf_value_1);
               const __m128 vt_i_2{_mm_add_ps(_mm_set1_ps(j),c8_11)};
               const __m128 vsin_arg_2{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_2,vinvT))))};
               const __m128 vcos_arg_2{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_2,vinv2T))};
               const __m128 vpsf_value_2{_mm_cos_ps(_mm_sub_ps(vcos_arg_2,vsin_arg_2))};
               _mm_store_ps(&this->m_psfunc.m_data[i+8ull],vpsf_value_2);
          }

          for(; (i+7ull) < this->m_psfunc_nsamples; i += 8ull,j += 8.0f) 
          {
               const __m128 vt_i_0{_mm_add_ps(_mm_set1_ps(j),c0_3)};
               const __m128 vsin_arg_0{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_0,vinvT))))};
               const __m128 vcos_arg_0{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_0,vinv2T))};
               const __m128 vpsf_value_0{_mm_cos_ps(_mm_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m128 vt_i_1{_mm_add_ps(_mm_set1_ps(j),c4_7)};
               const __m128 vsin_arg_1{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_1,vinvT))))};
               const __m128 vcos_arg_1{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_1,vinv2T))};
               const __m128 vpsf_value_1{_mm_cos_ps(_mm_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm_store_ps(&this->m_psfunc.m_data[i+4ull],vpsf_value_1);
          }

          for(; (i+3ull) < this->m_psfunc_nsamples; i += 4ull,j += 4.0f) 
          {
               const __m128 vt_i_0{_mm_add_ps(_mm_set1_ps(j),c0_3)};
               const __m128 vsin_arg_0{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_0,vinvT))))};
               const __m128 vcos_arg_0{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_0,vinv2T))};
               const __m128 vpsf_value_0{_mm_cos_ps(_mm_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
          }

          for(; (i+0ull) < this->m_psfunc_nsamples; i += 1ull,j += 1.0f)
          {
                //const float t_j{j};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_0{0.25f*ceph_sinf(this->m_k*6.283185307179586476925286766559f*j*invT)};
             const float cos_arg_0{this->m_k*3.14159265358979323846264338328f*j*inv2T};
             const float psf_val_0{ceph_cosf(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i] = psf_val_0;
#else 
             const float sin_arg_0{0.25f*std::sin(this->m_k*6.283185307179586476925286766559f*j*invT)};
             const float cos_arg_0{this->m_k*3.14159265358979323846264338328f*j*inv2T};
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
               const __m128 vsin_arg_0{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_0,vinvT))))};
               const __m128 vcos_arg_0{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_0,vinv2T))};
               const __m128 vpsf_value_0{_mm_cos_ps(_mm_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m128 vt_i_1{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull])};
               const __m128 vsin_arg_1{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_1,vinvT))))};
               const __m128 vcos_arg_1{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_1,vinv2T))};
               const __m128 vpsf_value_1{_mm_cos_ps(_mm_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm_store_ps(&this->m_psfunc.m_data[i+4ull],vpsf_value_1);
               const __m128 vt_i_2{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull])};
               const __m128 vsin_arg_2{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_2,vinvT))))};
               const __m128 vcos_arg_2{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_2,vinv2T))};
               const __m128 vpsf_value_2{_mm_cos_ps(_mm_sub_ps(vcos_arg_2,vsin_arg_2))};
               _mm_store_ps(&this->m_psfunc.m_data[i+8ull],vpsf_value_2);
               const __m128 vt_i_3{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+12ull])};
               const __m128 vsin_arg_3{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_3,vinvT))))};
               const __m128 vcos_arg_3{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_3,vinv2T))};
               const __m128 vpsf_value_3{_mm_cos_ps(_mm_sub_ps(vcos_arg_3,vsin_arg_3))};
               _mm_store_ps(&this->m_psfunc.m_data[i+12ull],vpsf_value_3);
          }

          for(; (i+11ull) < this->m_psfunc_nsamples; i += 12ull) 
          {
               const __m128 vt_i_0{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull])};
               const __m128 vsin_arg_0{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_0,vinvT))))};
               const __m128 vcos_arg_0{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_0,vinv2T))};
               const __m128 vpsf_value_0{_mm_cos_ps(_mm_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m128 vt_i_1{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull])};
               const __m128 vsin_arg_1{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_1,vinvT))))};
               const __m128 vcos_arg_1{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_1,vinv2T))};
               const __m128 vpsf_value_1{_mm_cos_ps(_mm_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm_store_ps(&this->m_psfunc.m_data[i+4ull],vpsf_value_1);
               const __m128 vt_i_2{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull])};
               const __m128 vsin_arg_2{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_2,vinvT))))};
               const __m128 vcos_arg_2{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_2,vinv2T))};
               const __m128 vpsf_value_2{_mm_cos_ps(_mm_sub_ps(vcos_arg_2,vsin_arg_2))};
               _mm_store_ps(&this->m_psfunc.m_data[i+8ull],vpsf_value_2);
          }

          for(; (i+7ull) < this->m_psfunc_nsamples; i += 8ull) 
          {
               const __m128 vt_i_0{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull])};
               const __m128 vsin_arg_0{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_0,vinvT))))};
               const __m128 vcos_arg_0{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_0,vinv2T))};
               const __m128 vpsf_value_0{_mm_cos_ps(_mm_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m128 vt_i_1{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull])};
               const __m128 vsin_arg_1{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_1,vinvT))))};
               const __m128 vcos_arg_1{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_1,vinv2T))};
               const __m128 vpsf_value_1{_mm_cos_ps(_mm_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm_store_ps(&this->m_psfunc.m_data[i+4ull],vpsf_value_1);
          }

          for(; (i+3ull) < this->m_psfunc_nsamples; i += 4ull) 
          {
               const __m128 vt_i_0{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull])};
               const __m128 vsin_arg_0{_mm_mul_ps(v025,_mm_sin_ps(_mm_mul_ps(_mm_mul_ps(vk,v2pi),_mm_mul_ps(vt_i_0,vinvT))))};
               const __m128 vcos_arg_0{_mm_mul_ps(_mm_mul_ps(vk,vpi),_mm_mul_ps(vt_i_0,vinv2T))};
               const __m128 vpsf_value_0{_mm_cos_ps(_mm_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
          }

          for(; (i+0ull) < this->m_psfunc_nsamples; i += 1ull)
          {
             const float t_i{static_cast<float>(gms::math::LUT_loop_indices_2257_align16[i])};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_0{0.25f*ceph_sinf(this->m_k*6.283185307179586476925286766559f*t_i*invT)};
             const float cos_arg_0{this->m_k*3.14159265358979323846264338328f*t_i*inv2T};
             const float psf_val_0{ceph_cosf(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i] = psf_val_0;
#else 
             const float sin_arg_0{0.25f*std::sin(this->m_k*6.283185307179586476925286766559f*t_i*invT)};
             const float cos_arg_0{this->m_k*3.14159265358979323846264338328f*t_i*inv2T};
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
     const __m256 vk{_mm256_set1_ps(this->m_k)};
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
          _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[32],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[48],_MM_HINT_T0);
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
               const __m256 vsin_arg_0{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(_mm256_mul_ps(vk,v2pi),_mm256_mul_ps(vt_i_0,vinvT))))};
               const __m256 vcos_arg_0{_mm256_mul_ps(_mm256_mul_ps(vk,vpi),_mm256_mul_ps(vt_i_0,vinv2T))};
               const __m256 vpsf_value_0{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m256 vt_i_1{_mm256_add_ps(_mm256_set1_ps(j),c8_15)};
               const __m256 vsin_arg_1{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(_mm256_mul_ps(vk,v2pi),_mm256_mul_ps(vt_i_1,vinvT))))};
               const __m256 vcos_arg_1{_mm256_mul_ps(_mm256_mul_ps(vk,vpi),_mm256_mul_ps(vt_i_1,vinv2T))};
               const __m256 vpsf_value_1{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+8ull],vpsf_value_1);
               const __m256 vt_i_2{_mm256_add_ps(_mm256_set1_ps(j),c16_23)};
               const __m256 vsin_arg_2{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(_mm256_mul_ps(vk,v2pi),_mm256_mul_ps(vt_i_2,vinvT))))};
               const __m256 vcos_arg_2{_mm256_mul_ps(_mm256_mul_ps(vk,vpi),_mm256_mul_ps(vt_i_2,vinv2T))};
               const __m256 vpsf_value_2{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_2,vsin_arg_2))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+16ull],vpsf_value_2);
               const __m256 vt_i_3{_mm256_add_ps(_mm256_set1_ps(j),c16_23)};
               const __m256 vsin_arg_3{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(_mm256_mul_ps(vk,v2pi),_mm256_mul_ps(vt_i_3,vinvT))))};
               const __m256 vcos_arg_3{_mm256_mul_ps(_mm256_mul_ps(vk,vpi),_mm256_mul_ps(vt_i_3,vinv2T))};
               const __m256 vpsf_value_3{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_3,vsin_arg_3))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+24ull],vpsf_value_3);
          }

          for(; (i+15ull) < this->m_psfunc_nsamples; i += 16ull,j += 16.0f) 
          {
                const __m256 vt_i_0{_mm256_add_ps(_mm256_set1_ps(j),c0_7)};
               const __m256 vsin_arg_0{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(_mm256_mul_ps(vk,v2pi),_mm256_mul_ps(vt_i_0,vinvT))))};
               const __m256 vcos_arg_0{_mm256_mul_ps(_mm256_mul_ps(vk,vpi),_mm256_mul_ps(vt_i_0,vinv2T))};
               const __m256 vpsf_value_0{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m256 vt_i_1{_mm256_add_ps(_mm256_set1_ps(j),c8_15)};
               const __m256 vsin_arg_1{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(_mm256_mul_ps(vk,v2pi),_mm256_mul_ps(vt_i_1,vinvT))))};
               const __m256 vcos_arg_1{_mm256_mul_ps(_mm256_mul_ps(vk,vpi),_mm256_mul_ps(vt_i_1,vinv2T))};
               const __m256 vpsf_value_1{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+8ull],vpsf_value_1);
          }

          for(; (i+7ull) < this->m_psfunc_nsamples; i += 8ull,j += 8.0f) 
          {
               const __m256 vt_i_0{_mm256_add_ps(_mm256_set1_ps(j),c0_7)};
               const __m256 vsin_arg_0{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(_mm256_mul_ps(vk,v2pi),_mm256_mul_ps(vt_i_0,vinvT))))};
               const __m256 vcos_arg_0{_mm256_mul_ps(_mm256_mul_ps(vk,vpi),_mm256_mul_ps(vt_i_0,vinv2T))};
               const __m256 vpsf_value_0{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
          }

          for(; (i+0ull) < this->m_psfunc_nsamples; i += 1ull,j += 1.0f)
          {
                //const float t_j{j};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_0{0.25f*ceph_sinf(this->m_k*6.283185307179586476925286766559f*j*invT)};
             const float cos_arg_0{this->m_k*3.14159265358979323846264338328f*j*inv2T};
             const float psf_val_0{ceph_cosf(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i] = psf_val_0;
#else 
             const float sin_arg_0{0.25f*std::sin(this->m_k*6.283185307179586476925286766559f*j*invT)};
             const float cos_arg_0{this->m_k*3.14159265358979323846264338328f*j*inv2T};
             const float psf_val_0{std::cos(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i] = psf_val_0;
#endif 
          }
     }
     else 
     {
           for(i = 0ull; (i+31ull) < this->m_psfunc_nsamples; i += 32ull) 
           {
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align32[i+0ull],_MM_HINT_T0);
               const __m256 vt_i_0{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+0ull])};
               const __m256 vsin_arg_0{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(_mm256_mul_ps(vk,v2pi),_mm256_mul_ps(vt_i_0,vinvT))))};
               const __m256 vcos_arg_0{_mm256_mul_ps(_mm256_mul_ps(vk,vpi),_mm256_mul_ps(vt_i_0,vinv2T))};
               const __m256 vpsf_value_0{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m256 vt_i_1{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+8ull])};
               const __m256 vsin_arg_1{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(_mm256_mul_ps(vk,v2pi),_mm256_mul_ps(vt_i_1,vinvT))))};
               const __m256 vcos_arg_1{_mm256_mul_ps(_mm256_mul_ps(vk,vpi),_mm256_mul_ps(vt_i_1,vinv2T))};
               const __m256 vpsf_value_1{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+8ull],vpsf_value_1);
              _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+16ull],_MM_HINT_T0);
               const __m256 vt_i_2{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+16ull])};
               const __m256 vsin_arg_2{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(_mm256_mul_ps(vk,v2pi),_mm256_mul_ps(vt_i_2,vinvT))))};
               const __m256 vcos_arg_2{_mm256_mul_ps(_mm256_mul_ps(vk,vpi),_mm256_mul_ps(vt_i_2,vinv2T))};
               const __m256 vpsf_value_2{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_2,vsin_arg_2))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+16ull],vpsf_value_2);
               const __m256 vt_i_3{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+24ull])};
               const __m256 vsin_arg_3{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(_mm256_mul_ps(vk,v2pi),_mm256_mul_ps(vt_i_3,vinvT))))};
               const __m256 vcos_arg_3{_mm256_mul_ps(_mm256_mul_ps(vk,vpi),_mm256_mul_ps(vt_i_3,vinv2T))};
               const __m256 vpsf_value_3{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_3,vsin_arg_3))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+24ull],vpsf_value_3);
          }

          for(; (i+15ull) < this->m_psfunc_nsamples; i += 16ull) 
          {
                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+0ull],_MM_HINT_T0);
               const __m256 vt_i_0{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+0ull])};
               const __m256 vsin_arg_0{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(_mm256_mul_ps(vk,v2pi),_mm256_mul_ps(vt_i_0,vinvT))))};
               const __m256 vcos_arg_0{_mm256_mul_ps(_mm256_mul_ps(vk,vpi),_mm256_mul_ps(vt_i_0,vinv2T))};
               const __m256 vpsf_value_0{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m256 vt_i_1{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+8ull])};
               const __m256 vsin_arg_1{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(_mm256_mul_ps(vk,v2pi),_mm256_mul_ps(vt_i_1,vinvT))))};
               const __m256 vcos_arg_1{_mm256_mul_ps(_mm256_mul_ps(vk,vpi),_mm256_mul_ps(vt_i_1,vinv2T))};
               const __m256 vpsf_value_1{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+8ull],vpsf_value_1);
          }

          for(; (i+7ull) < this->m_psfunc_nsamples; i += 8ull) 
          {
               const __m256 vt_i_0{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+0ull])};
               const __m256 vsin_arg_0{_mm256_mul_ps(v025,_mm256_sin_ps(_mm256_mul_ps(_mm256_mul_ps(vk,v2pi),_mm256_mul_ps(vt_i_0,vinvT))))};
               const __m256 vcos_arg_0{_mm256_mul_ps(_mm256_mul_ps(vk,vpi),_mm256_mul_ps(vt_i_0,vinv2T))};
               const __m256 vpsf_value_0{_mm256_cos_ps(_mm256_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm256_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
          }

          for(; (i+0ull) < this->m_psfunc_nsamples; i += 1ull)
          {
             const float t_i{static_cast<float>(gms::math::LUT_loop_indices_2257_align16[i])};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_0{0.25f*ceph_sinf(this->m_k*6.283185307179586476925286766559f*t_i*invT)};
             const float cos_arg_0{this->m_k*3.14159265358979323846264338328f*t_i*inv2T};
             const float psf_val_0{ceph_cosf(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i] = psf_val_0;
#else 
             const float sin_arg_0{0.25f*std::sin(this->m_k*6.283185307179586476925286766559f*t_i*invT)};
             const float cos_arg_0{this->m_k*3.14159265358979323846264338328f*t_i*inv2T};
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
::generate_pulse_shaping_function_avx512_u4x(const bool do_prefetch_constants)
{
     using namespace gms::math;
     const float f_psfunc_nsamples{static_cast<float>(this->m_psfunc_nsamples)};
     const float invT{1.0f/f_psfunc_nsamples};
     const float inv2T{1.0f/(2.0f*f_psfunc_nsamples)};
     constexpr std::size_t LUT_loop_idx_threshold{2257ull};
     const __m512 vinvT{_mm512_set1_ps(invT)};
     const __m512 vinv2T{_mm512_set1_ps(inv2T)};
     const __m512 vk{_mm512_set1_ps(this->m_k)};
     __m512 c0_15,c16_31;
     __m512 c32_47,c48_63;
     __m512 vpi,v2pi;
     __m512 v025;
     std::size_t i;
     float j;
     if(!do_prefetch_constants)
     {
         c0_15  = _mm512_setr_ps(0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f);
         c16_31 = _mm512_setr_ps(16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f,24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f);
         c32_47 = _mm512_setr_ps(32.0f,33.0f,34.0f,35.0f,36.0f,37.0f,38.0f,39.0f,40.0f,41.0f,42.0f,43.0f,44.0f,45.0f,46.0f,47.0f);
         c48_63 = _mm512_setr_ps(48.0f,49.0f,50.0f,51.0f,52.0f,53.0f,54.0f,55.0f,56.0f,57.0f,58.0f,59.0f,60.0f,61.0f,62.0f,63.0f);
         vpi    = _mm512_set1_ps(3.14159265358979323846264338328f);
         v2pi   = _mm512_set1_ps(6.283185307179586476925286766559f); 
         v025   = _mm512_set1_ps(0.25f);
     }
     else 
     {
          __ATTR_ALIGN__(64) const float prefetched_constants[112] = {
                                0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f,
                                16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f,24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f,
                                32.0f,33.0f,34.0f,35.0f,36.0f,37.0f,38.0f,39.0f,40.0f,41.0f,42.0f,43.0f,44.0f,45.0f,46.0f,47.0f,
                                48.0f,49.0f,50.0f,51.0f,52.0f,53.0f,54.0f,55.0f,56.0f,57.0f,58.0f,59.0f,60.0f,61.0f,62.0f,63.0f,
                                3.14159265358979323846264338328f,3.14159265358979323846264338328f,3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f,3.14159265358979323846264338328f,3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f,3.14159265358979323846264338328f,3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f,3.14159265358979323846264338328f,3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f,3.14159265358979323846264338328f,3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f,6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,
                                0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f};

          _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[32],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[48],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[64],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[80],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[96],_MM_HINT_T0);
          c0_15  = _mm512_load_ps(&prefetched_constants[0]);
          c16_31 = _mm512_load_ps(&prefetched_constants[16]); 
          c32_47 = _mm512_load_ps(&prefetched_constants[32]);
          c48_63 = _mm512_load_ps(&prefetched_constants[48]);
          vpi    = _mm512_load_ps(&prefetched_constants[64]);
          v2pi   = _mm512_load_ps(&prefetched_constants[80]);
          v025   = _mm512_load_ps(&prefetched_constants[96]);
     }

     if(this->m_psfunc_nsamples>LUT_loop_idx_threshold)
     {
          for(i = 0ull,j = 0.0f; (i+63ull) < this->m_psfunc_nsamples; i += 64ull,j += 64.0f) 
          {
               const __m512 vt_i_0{_mm512_add_ps(_mm512_set1_ps(j),c0_15)};
               const __m512 vsin_arg_0{_mm512_mul_ps(v025,_mm512_sin_ps(_mm512_mul_ps(_mm512_mul_ps(vk,v2pi),_mm512_mul_ps(vt_i_0,vinvT))))};
               const __m512 vcos_arg_0{_mm512_mul_ps(_mm512_mul_ps(vk,vpi),_mm512_mul_ps(vt_i_0,vinv2T))};
               const __m512 vpsf_value_0{_mm512_cos_ps(_mm512_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm512_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m512 vt_i_1{_mm512_add_ps(_mm512_set1_ps(j),c16_31)};
               const __m512 vsin_arg_1{_mm512_mul_ps(v025,_mm512_sin_ps(_mm512_mul_ps(_mm512_mul_ps(vk,v2pi),_mm512_mul_ps(vt_i_1,vinvT))))};
               const __m512 vcos_arg_1{_mm512_mul_ps(_mm512_mul_ps(vk,vpi),_mm512_mul_ps(vt_i_1,vinv2T))};
               const __m512 vpsf_value_1{_mm512_cos_ps(_mm512_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm512_store_ps(&this->m_psfunc.m_data[i+16ull],vpsf_value_1);
               const __m512 vt_i_2{_mm512_add_ps(_mm512_set1_ps(j),c32_47)};
               const __m512 vsin_arg_2{_mm512_mul_ps(v025,_mm512_sin_ps(_mm512_mul_ps(_mm512_mul_ps(vk,v2pi),_mm512_mul_ps(vt_i_2,vinvT))))};
               const __m512 vcos_arg_2{_mm512_mul_ps(_mm512_mul_ps(vk,vpi),_mm512_mul_ps(vt_i_2,vinv2T))};
               const __m512 vpsf_value_2{_mm512_cos_ps(_mm512_sub_ps(vcos_arg_2,vsin_arg_2))};
               _mm512_store_ps(&this->m_psfunc.m_data[i+32ull],vpsf_value_2);
               const __m512 vt_i_3{_mm512_add_ps(_mm512_set1_ps(j),c48_63)};
               const __m512 vsin_arg_3{_mm512_mul_ps(v025,_mm512_sin_ps(_mm512_mul_ps(_mm512_mul_ps(vk,v2pi),_mm512_mul_ps(vt_i_3,vinvT))))};
               const __m512 vcos_arg_3{_mm512_mul_ps(_mm512_mul_ps(vk,vpi),_mm512_mul_ps(vt_i_3,vinv2T))};
               const __m512 vpsf_value_3{_mm512_cos_ps(_mm512_sub_ps(vcos_arg_3,vsin_arg_3))};
               _mm512_store_ps(&this->m_psfunc.m_data[i+48ull],vpsf_value_3);
          }

          for(; (i+31ull) < this->m_psfunc_nsamples; i += 32ull,j += 32.0f) 
          {
               const __m512 vt_i_0{_mm512_add_ps(_mm512_set1_ps(j),c0_15)};
               const __m512 vsin_arg_0{_mm512_mul_ps(v025,_mm512_sin_ps(_mm512_mul_ps(_mm512_mul_ps(vk,v2pi),_mm512_mul_ps(vt_i_0,vinvT))))};
               const __m512 vcos_arg_0{_mm512_mul_ps(_mm512_mul_ps(vk,vpi),_mm512_mul_ps(vt_i_0,vinv2T))};
               const __m512 vpsf_value_0{_mm512_cos_ps(_mm512_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm512_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m512 vt_i_1{_mm512_add_ps(_mm512_set1_ps(j),c16_31)};
               const __m512 vsin_arg_1{_mm512_mul_ps(v025,_mm512_sin_ps(_mm512_mul_ps(_mm512_mul_ps(vk,v2pi),_mm512_mul_ps(vt_i_1,vinvT))))};
               const __m512 vcos_arg_1{_mm512_mul_ps(_mm512_mul_ps(vk,vpi),_mm512_mul_ps(vt_i_1,vinv2T))};
               const __m512 vpsf_value_1{_mm512_cos_ps(_mm512_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm512_store_ps(&this->m_psfunc.m_data[i+16ull],vpsf_value_1);
          }

          for(; (i+15ull) < this->m_psfunc_nsamples; i += 16ull,j += 16.0f) 
          {
               const __m512 vt_i_0{_mm512_add_ps(_mm512_set1_ps(j),c0_15)};
               const __m512 vsin_arg_0{_mm512_mul_ps(v025,_mm512_sin_ps(_mm512_mul_ps(_mm512_mul_ps(vk,v2pi),_mm512_mul_ps(vt_i_0,vinvT))))};
               const __m512 vcos_arg_0{_mm512_mul_ps(_mm512_mul_ps(vk,vpi),_mm512_mul_ps(vt_i_0,vinv2T))};
               const __m512 vpsf_value_0{_mm512_cos_ps(_mm512_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm512_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
          }

          for(; (i+0ull) < this->m_psfunc_nsamples; i += 1ull,j += 1.0f)
          {
                //const float t_j{j};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_0{0.25f*ceph_sinf(this->m_k*6.283185307179586476925286766559f*j*invT)};
             const float cos_arg_0{this->m_k*3.14159265358979323846264338328f*j*inv2T};
             const float psf_val_0{ceph_cosf(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i] = psf_val_0;
#else 
             const float sin_arg_0{0.25f*std::sin(this->m_k*6.283185307179586476925286766559f*j*invT)};
             const float cos_arg_0{this->m_k*3.14159265358979323846264338328f*j*inv2T};
             const float psf_val_0{std::cos(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i] = psf_val_0;
#endif 
          }
     }
     else 
     {
           for(i = 0ull; (i+63ull) < this->m_psfunc_nsamples; i += 64ull) 
           {
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+0ull],_MM_HINT_T0);
               const __m512 vt_i_0{_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+0ull])};
               const __m512 vsin_arg_0{_mm512_mul_ps(v025,_mm512_sin_ps(_mm512_mul_ps(_mm512_mul_ps(vk,v2pi),_mm512_mul_ps(vt_i_0,vinvT))))};
               const __m512 vcos_arg_0{_mm512_mul_ps(_mm512_mul_ps(vk,vpi),_mm512_mul_ps(vt_i_0,vinv2T))};
               const __m512 vpsf_value_0{_mm512_cos_ps(_mm512_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+16ull],_MM_HINT_T0);
               _mm512_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m512 vt_i_1{_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+16ull])};
               const __m512 vsin_arg_1{_mm512_mul_ps(v025,_mm512_sin_ps(_mm512_mul_ps(_mm512_mul_ps(vk,v2pi),_mm512_mul_ps(vt_i_1,vinvT))))};
               const __m512 vcos_arg_1{_mm512_mul_ps(_mm512_mul_ps(vk,vpi),_mm512_mul_ps(vt_i_1,vinv2T))};
               const __m512 vpsf_value_1{_mm512_cos_ps(_mm512_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+32ull],_MM_HINT_T0);
               _mm512_store_ps(&this->m_psfunc.m_data[i+16ull],vpsf_value_1);
               const __m512 vt_i_2{_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+32ull])};
               const __m512 vsin_arg_2{_mm512_mul_ps(v025,_mm512_sin_ps(_mm512_mul_ps(_mm512_mul_ps(vk,v2pi),_mm512_mul_ps(vt_i_2,vinvT))))};
               const __m512 vcos_arg_2{_mm512_mul_ps(_mm512_mul_ps(vk,vpi),_mm512_mul_ps(vt_i_2,vinv2T))};
               const __m512 vpsf_value_2{_mm512_cos_ps(_mm512_sub_ps(vcos_arg_2,vsin_arg_2))};
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+48ull],_MM_HINT_T0);
               _mm512_store_ps(&this->m_psfunc.m_data[i+32ull],vpsf_value_2);
               const __m512 vt_i_3{_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+48ull])};
               const __m512 vsin_arg_3{_mm512_mul_ps(v025,_mm512_sin_ps(_mm512_mul_ps(_mm512_mul_ps(vk,v2pi),_mm512_mul_ps(vt_i_3,vinvT))))};
               const __m512 vcos_arg_3{_mm512_mul_ps(_mm512_mul_ps(vk,vpi),_mm512_mul_ps(vt_i_3,vinv2T))};
               const __m512 vpsf_value_3{_mm512_cos_ps(_mm512_sub_ps(vcos_arg_3,vsin_arg_3))};
               _mm512_store_ps(&this->m_psfunc.m_data[i+48ull],vpsf_value_3);
          }

          for(; (i+31ull) < this->m_psfunc_nsamples; i += 32ull) 
          {
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+0ull],_MM_HINT_T0);
               const __m512 vt_i_0{_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+0ull])};
               const __m512 vsin_arg_0{_mm512_mul_ps(v025,_mm512_sin_ps(_mm512_mul_ps(_mm512_mul_ps(vk,v2pi),_mm512_mul_ps(vt_i_0,vinvT))))};
               const __m512 vcos_arg_0{_mm512_mul_ps(_mm512_mul_ps(vk,vpi),_mm512_mul_ps(vt_i_0,vinv2T))};
               const __m512 vpsf_value_0{_mm512_cos_ps(_mm512_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+16ull],_MM_HINT_T0);
               _mm512_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
               const __m512 vt_i_1{_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+16ull])};
               const __m512 vsin_arg_1{_mm512_mul_ps(v025,_mm512_sin_ps(_mm512_mul_ps(_mm512_mul_ps(vk,v2pi),_mm512_mul_ps(vt_i_1,vinvT))))};
               const __m512 vcos_arg_1{_mm512_mul_ps(_mm512_mul_ps(vk,vpi),_mm512_mul_ps(vt_i_1,vinv2T))};
               const __m512 vpsf_value_1{_mm512_cos_ps(_mm512_sub_ps(vcos_arg_1,vsin_arg_1))};
               _mm512_store_ps(&this->m_psfunc.m_data[i+16ull],vpsf_value_1);
          }

          for(; (i+15ull) < this->m_psfunc_nsamples; i += 16ull) 
          {
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+0ull],_MM_HINT_T0);
               const __m512 vt_i_0{_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+0ull])};
               const __m512 vsin_arg_0{_mm512_mul_ps(v025,_mm512_sin_ps(_mm512_mul_ps(_mm512_mul_ps(vk,v2pi),_mm512_mul_ps(vt_i_0,vinvT))))};
               const __m512 vcos_arg_0{_mm512_mul_ps(_mm512_mul_ps(vk,vpi),_mm512_mul_ps(vt_i_0,vinv2T))};
               const __m512 vpsf_value_0{_mm512_cos_ps(_mm512_sub_ps(vcos_arg_0,vsin_arg_0))};
               _mm512_store_ps(&this->m_psfunc.m_data[i+0ull],vpsf_value_0);
          }

          for(; (i+0ull) < this->m_psfunc_nsamples; i += 1ull)
          {
             const float t_i{static_cast<float>(gms::math::LUT_loop_indices_2257_align16[i])};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
             const float sin_arg_0{0.25f*ceph_sinf(this->m_k*6.283185307179586476925286766559f*t_i*invT)};
             const float cos_arg_0{this->m_k*3.14159265358979323846264338328f*t_i*inv2T};
             const float psf_val_0{ceph_cosf(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i] = psf_val_0;
#else 
             const float sin_arg_0{0.25f*std::sin(this->m_k*6.283185307179586476925286766559f*t_i*invT)};
             const float cos_arg_0{this->m_k*3.14159265358979323846264338328f*t_i*inv2T};
             const float psf_val_0{std::cos(cos_arg_0-sin_arg_0)};
             this->m_psfunc.m_data[i] = psf_val_0;
#endif 
          }
     }

     return (0);
}

template<gms::radiolocation
          ::sinusoidal_fsk_t::I_channel_bitstream_optimized_path I_ch_path,
          gms::radiolocation
          ::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path Q_ch_path,
          gms::radiolocation
          ::sinusoidal_fsk_t::pulse_shaping_function_optimized_path psf_path>
std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_u8x(const iq_rectw_bitstream_vsequence_t &iq_bitstream_vsequence,
                          const std::int32_t I_ch_samples_len,
                          const std::int32_t Q_ch_samples_len,
                          const bool I_ch_do_const_prefetch,
                          const bool Q_ch_do_const_prefetch,
                          const bool psfunc_do_const_prefetch,
                          const bool use_iq_bitstream_vsequence)
{
     using namespace gms::math;
     if(__builtin_expect(iq_bitstream_vsequence.m_I_nsamples!=this->m_I_ch_nsamples,0)) { return (-19);}
     if(__builtin_expect(iq_bitstream_vsequence.m_Q_nsamples!=this->m_Q_ch_nsamples,0)) { return (-20);};

     if(__builtin_expect(use_iq_bitstream_vsequence==true,1))
     {
          std::memcpy(&this->m_I_ch_bitstream.m_data[0],&iq_bitstream_vsequence.m_I_vsequence.m_data[0],
                      this->m_I_ch_nsamples*sizeof(float));
          std::memcpy(&this->m_Q_ch_bitstream.m_data[0],&iq_bitstream_vsequence.m_Q_vsequence.m_data[0],
                      this->m_Q_ch_nsamples*sizeof(float));
     }
     else 
     {

     
          if constexpr (I_ch_path == I_channel_bitstream_optimized_path::default_scalar_path)
          {
              std::int32_t ret_stat_scalar = this->generate_I_channel_bitstream(I_ch_samples_len);
              if(__builtin_expect(ret_stat_scalar<0,0)) { return (-1);}
          }
          else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_sse_path)
          {
              std::int32_t ret_stat_sse    = this->generate_I_channel_bitstream_sse();
              if(__builtin_expect(ret_stat_sse<0,0))    { return (-2);}
          }
          else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_sse_u4x_path)
          {
              std::int32_t ret_stat_sse_u4x = this->generate_I_channel_bitstream_sse_u4x();
              if(__builtin_expect(ret_stat_sse_u4x<0,0)) {return (-3);}
          }
          else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_avx_path)
          {
              std::int32_t ret_stat_avx     = this->generate_I_channel_bitstream_avx(I_ch_do_const_prefetch);
              if(__builtin_expect(ret_stat_avx<0,0))     {return (-4);}
          }
          else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_avx_u4x_path)
          {
              std::int32_t ret_stat_avx_u4x = this->generate_I_channel_bitstream_avx_u4x(I_ch_do_const_prefetch);
              if(__builtin_expect(ret_stat_avx_u4x<0,0)) { return (-5);}
          }
          else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_avx512_path)
          {
              std::int32_t ret_stat_avx512  = this->generate_I_channel_bitstream_avx512(I_ch_do_const_prefetch);
              if(__builtin_expect(ret_stat_avx512<0,0))  { return (-6);}
          }
          else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_avx512_u4x_path)
          {
              std::int32_t ret_stat_avx512_u4x= this->generate_I_channel_bitstream_avx512_u4x(I_ch_do_const_prefetch);
              if(__builtin_expect(ret_stat_avx512_u4x<0,0)){ return (-7);}
          }
          else
          {
               return (-9999);
          }

          if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::default_scalar_path)
          {
              std::int32_t ret_stat_scalar = this->generate_Q_channel_bitstream(Q_ch_samples_len);
              if(__builtin_expect(ret_stat_scalar<0,0)) { return (-8);}
          }
          else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_sse_path)
          {
              std::int32_t ret_stat_sse    = this->generate_Q_channel_bitstream_sse();
              if(__builtin_expect(ret_stat_sse<0,0))    { return (-9);}
          }
          else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_sse_u4x_path)
          {
              std::int32_t ret_stat_sse_u4x = this->generate_Q_channel_bitstream_sse_u4x();
              if(__builtin_expect(ret_stat_sse_u4x<0,0)) {return (-10);}
          }
          else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_avx_path)
          {
              std::int32_t ret_stat_avx     = this->generate_Q_channel_bitstream_avx(I_ch_do_const_prefetch);
              if(__builtin_expect(ret_stat_avx<0,0))     {return (-11);}
          }
          else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_avx_u4x_path)
          {
              std::int32_t ret_stat_avx_u4x = this->generate_Q_channel_bitstream_avx_u4x(I_ch_do_const_prefetch);
              if(__builtin_expect(ret_stat_avx_u4x<0,0)) { return (-12);}
          }
          else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_avx512_path)
          {
              std::int32_t ret_stat_avx512  = this->generate_Q_channel_bitstream_avx512(I_ch_do_const_prefetch);
              if(__builtin_expect(ret_stat_avx512<0,0))  { return (-13);}
          }
          else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_avx512_u4x_path)
          {
              std::int32_t ret_stat_avx512_u4x= this->generate_Q_channel_bitstream_avx512_u4x(I_ch_do_const_prefetch);
              if(__builtin_expect(ret_stat_avx512_u4x<0,0)){ return (-14);}
          }
          else 
          {
             return (-9999);
          }
     }

          if constexpr (psf_path == pulse_shaping_function_optimized_path::default_scalar_path)
          {
              std::int32_t ret_stat_scalar = this->generate_pulse_shaping_function_u8x();
              if(__builtin_expect(ret_stat_scalar<0,0)) { return (-15);}
          }
          else if constexpr (psf_path == pulse_shaping_function_optimized_path::vector_sse_path)
          {
              std::int32_t ret_stat_sse    = this->generate_pulse_shaping_function_sse_u4x(psfunc_do_const_prefetch);
              if(__builtin_expect(ret_stat_sse<0,0))    { return (-16);}
          }
          else if constexpr (psf_path == pulse_shaping_function_optimized_path::vector_avx_path)
          {
              std::int32_t ret_stat_avx    = this->generate_pulse_shaping_function_avx_u4x(psfunc_do_const_prefetch);
              if(__builtin_expect(ret_stat_avx<0,0))    { return (-17);}
          }
          else if constexpr (psf_path == pulse_shaping_function_optimized_path::vector_avx512_path)
          {
              std::int32_t ret_stat_avx512 = this->generate_pulse_shaping_function_avx512_u4x(psfunc_do_const_prefetch);
              if(__builtin_expect(ret_stat_avx512<0,0))  { return (-18);}
          }
          else 
          {
               return (-9999);
          }
     constexpr std::size_t LUT_loop_idx_threshold{2257ull};
     constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
     constexpr float C314159265358979323846264338328{3.14159265358979323846264338328f};
     std::size_t i,j;
     if(this->m_sfsk_nsamples>LUT_loop_idx_threshold)
     {
          for(i = 0ull; i != ROUND_TO_EIGHT(this->m_sfsk_nsamples,8ull); i += 8ull) 
          {
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i],_MM_HINT_T0);
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i],_MM_HINT_T0);
               const float ik_0{this->m_I_ch_bitstream.m_data[i+0ull]};
               const float qk_0{this->m_Q_ch_bitstream.m_data[i+0ull]};
               const float dk_0{-ik_0*qk_0};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i],_MM_HINT_T0);
               const float t_i_0{static_cast<float>(i)};
               const float psfunc_0{this->m_psfunc.m_data[i+0ull]};
               const float Phik_0{ik_0==1.0f?0.0f:C314159265358979323846264338328};
               const float carrier_term1_0{C6283185307179586476925286766559*t_i_0};
               const float carrier_term2_0{std::fma(dk_0,psfunc_0,Phik_0)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_0{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_0+carrier_term2_0)};
#else 
               const float sfsk_sample_0{this->m_Ac*std::cos(this->m_ph0+carrier_term1_0+carrier_term2_0)};
#endif 
               this->m_sfsk_signal.m_data[i+0ull] = sfsk_sample_0;

               const float ik_1{this->m_I_ch_bitstream.m_data[i+1ull]};
               const float qk_1{this->m_Q_ch_bitstream.m_data[i+1ull]};
               const float dk_1{-ik_1*qk_1};
               const float t_i_1{static_cast<float>(i+1ull)};
               const float psfunc_1{this->m_psfunc.m_data[i+1ull]};
               const float Phik_1{ik_1==1.0f?0.0f:C314159265358979323846264338328};
               const float carrier_term1_1{C6283185307179586476925286766559*t_i_1};
               const float carrier_term2_1{std::fma(dk_1,psfunc_1,Phik_1)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_1{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#else 
               const float sfsk_sample_1{this->m_Ac*std::cos(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#endif 
               this->m_sfsk_signal.m_data[i+1ull] = sfsk_sample_1; 

               const float ik_2{this->m_I_ch_bitstream.m_data[i+2ull]};
               const float qk_2{this->m_Q_ch_bitstream.m_data[i+2ull]};
               const float dk_2{-ik_2*qk_2};
               const float t_i_2{static_cast<float>(i+2ull)};
               const float psfunc_2{this->m_psfunc.m_data[i+2ull]};
               const float Phik_2{ik_2==1.0f?0.0f:C314159265358979323846264338328};
               const float carrier_term1_2{C6283185307179586476925286766559*t_i_2};
               const float carrier_term2_2{std::fma(dk_2,psfunc_2,Phik_2)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_2{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_2+carrier_term2_2)};
#else 
               const float sfsk_sample_2{this->m_Ac*std::cos(this->m_ph0+carrier_term1_2+carrier_term2_2)};
#endif 
               this->m_sfsk_signal.m_data[i+2ull] = sfsk_sample_2;

               const float ik_3{this->m_I_ch_bitstream.m_data[i+3ull]};
               const float qk_3{this->m_Q_ch_bitstream.m_data[i+3ull]};
               const float dk_3{-ik_3*qk_3};
               const float t_i_3{static_cast<float>(i+3ull)};
               const float psfunc_3{this->m_psfunc.m_data[i+3ull]};
               const float Phik_3{ik_3==1.0f?0.0f:C314159265358979323846264338328};
               const float carrier_term1_3{C6283185307179586476925286766559*t_i_3};
               const float carrier_term2_3{std::fma(dk_3,psfunc_3,Phik_3)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_3{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_3+carrier_term2_3)};
#else 
               const float sfsk_sample_3{this->m_Ac*std::cos(this->m_ph0+carrier_term1_3+carrier_term2_3)};
#endif 
               this->m_sfsk_signal.m_data[i+3ull] = sfsk_sample_3;

               const float ik_4{this->m_I_ch_bitstream.m_data[i+4ull]};
               const float qk_4{this->m_Q_ch_bitstream.m_data[i+4ull]};
               const float dk_4{-ik_4*qk_4};
               const float t_i_4{static_cast<float>(i+4ull)};
               const float psfunc_4{this->m_psfunc.m_data[i+4ull]};
               const float Phik_4{ik_4==1.0f?0.0f:C314159265358979323846264338328};
               const float carrier_term1_4{C6283185307179586476925286766559*t_i_4};
               const float carrier_term2_4{std::fma(dk_4,psfunc_4,Phik_4)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_4{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_4+carrier_term2_4)};
#else 
               const float sfsk_sample_4{this->m_Ac*std::cos(this->m_ph0+carrier_term1_4+carrier_term2_4)};
#endif 
               this->m_sfsk_signal.m_data[i+4ull] = sfsk_sample_4;

               const float ik_5{this->m_I_ch_bitstream.m_data[i+5ull]};
               const float qk_5{this->m_Q_ch_bitstream.m_data[i+5ull]};
               const float dk_5{-ik_5*qk_5};
               const float t_i_5{static_cast<float>(i+5ull)};
               const float psfunc_5{this->m_psfunc.m_data[i+5ull]};
               const float Phik_5{ik_5==1.0f?0.0f:C314159265358979323846264338328};
               const float carrier_term1_5{C6283185307179586476925286766559*t_i_5};
               const float carrier_term2_5{std::fma(dk_5,psfunc_5,Phik_5)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_5{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_5+carrier_term2_5)};
#else 
               const float sfsk_sample_5{this->m_Ac*std::cos(this->m_ph0+carrier_term1_5+carrier_term2_5)};
#endif 
               this->m_sfsk_signal.m_data[i+5ull] = sfsk_sample_5;

               const float ik_6{this->m_I_ch_bitstream.m_data[i+6ull]};
               const float qk_6{this->m_Q_ch_bitstream.m_data[i+6ull]};
               const float dk_6{-ik_6*qk_6};
               const float t_i_6{static_cast<float>(i+6ull)};
               const float psfunc_6{this->m_psfunc.m_data[i+6ull]};
               const float Phik_6{ik_6==1.0f?0.0f:C314159265358979323846264338328};
               const float carrier_term1_6{C6283185307179586476925286766559*t_i_6};
               const float carrier_term2_6{std::fma(dk_6,psfunc_6,Phik_6)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_6{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_6+carrier_term2_6)};
#else 
               const float sfsk_sample_1{this->m_Ac*std::cos(this->m_ph0+carrier_term1_6+carrier_term2_6)};
#endif 
               this->m_sfsk_signal.m_data[i+6ull] = sfsk_sample_6;

               const float ik_7{this->m_I_ch_bitstream.m_data[i+7ull]};
               const float qk_7{this->m_Q_ch_bitstream.m_data[i+7ull]};
               const float dk_7{-ik_7*qk_7};
               const float t_i_7{static_cast<float>(i+7ull)};
               const float psfunc_7{this->m_psfunc.m_data[i+7ull]};
               const float Phik_7{ik_7==1.0f?0.0f:C314159265358979323846264338328};
               const float carrier_term1_7{C6283185307179586476925286766559*t_i_7};
               const float carrier_term2_7{std::fma(dk_7,psfunc_7,Phik_7)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_7{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_7+carrier_term2_7)};
#else 
               const float sfsk_sample_7{this->m_Ac*std::cos(this->m_ph0+carrier_term1_7+carrier_term2_7)};
#endif 
               this->m_sfsk_signal.m_data[i+7ull] = sfsk_sample_7;
          }

          for(j = i; j != this->m_sfsk_nsamples; ++j)  
          {
               const float ik_1{this->m_I_ch_bitstream.m_data[j]};
               const float qk_1{this->m_Q_ch_bitstream.m_data[j]};
               const float dk_1{-ik_1*qk_1};
               const float t_i_1{static_cast<float>(j)};
               const float psfunc_1{this->m_psfunc.m_data[i+1ull]};
               const float Phik_1{ik_1==1.0f?0.0f:C314159265358979323846264338328};
               const float carrier_term1_1{C6283185307179586476925286766559*t_i_1};
               const float carrier_term2_1{std::fma(dk_1,psfunc_1,Phik_1)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_1{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#else 
               const float sfsk_sample_1{this->m_Ac*std::cos(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#endif 
               this->m_sfsk_signal.m_data[j] = sfsk_sample_1; 

          }
     }
     else 
     {
          for(i = 0ull; i != ROUND_TO_EIGHT(this->m_sfsk_nsamples,8ull); i += 8ull) 
          {
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i],_MM_HINT_T0);
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i],_MM_HINT_T0);
               const float ik_0{this->m_I_ch_bitstream.m_data[i+0ull]};
               const float qk_0{this->m_Q_ch_bitstream.m_data[i+0ull]};
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i],_MM_HINT_T0);
               const float dk_0{-ik_0*qk_0};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i],_MM_HINT_T0);
               const float t_i_0{gms::math::LUT_loop_indices_2257_align16[i+0ull]};
               const float psfunc_0{this->m_psfunc.m_data[i+0ull]};
               const float Phik_0{ik_0==1.0f?0.0f:C314159265358979323846264338328};
               const float carrier_term1_0{C6283185307179586476925286766559*t_i_0};
               const float carrier_term2_0{std::fma(dk_0,psfunc_0,Phik_0)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_0{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_0+carrier_term2_0)};
#else 
               const float sfsk_sample_0{this->m_Ac*std::cos(this->m_ph0+carrier_term1_0+carrier_term2_0)};
#endif 
               this->m_sfsk_signal.m_data[i+0ull] = sfsk_sample_0;

               const float ik_1{this->m_I_ch_bitstream.m_data[i+1ull]};
               const float qk_1{this->m_Q_ch_bitstream.m_data[i+1ull]};
               const float dk_1{-ik_1*qk_1};
               const float t_i_1{gms::math::LUT_loop_indices_2257_align16[i+1ull]};
               const float psfunc_1{this->m_psfunc.m_data[i+1ull]};
               const float Phik_1{ik_1==1.0f?0.0f:C314159265358979323846264338328};
               const float carrier_term1_1{C6283185307179586476925286766559*t_i_1};
               const float carrier_term2_1{std::fma(dk_1,psfunc_1,Phik_1)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_1{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#else 
               const float sfsk_sample_1{this->m_Ac*std::cos(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#endif 
               this->m_sfsk_signal.m_data[i+1ull] = sfsk_sample_1; 

               const float ik_2{this->m_I_ch_bitstream.m_data[i+2ull]};
               const float qk_2{this->m_Q_ch_bitstream.m_data[i+2ull]};
               const float dk_2{-ik_2*qk_2};
               const float t_i_2{gms::math::LUT_loop_indices_2257_align16[i+2ull]};
               const float psfunc_2{this->m_psfunc.m_data[i+2ull]};
               const float Phik_2{ik_2==1.0f?0.0f:C314159265358979323846264338328};
               const float carrier_term1_2{C6283185307179586476925286766559*t_i_2};
               const float carrier_term2_2{std::fma(dk_2,psfunc_2,Phik_2)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_2{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_2+carrier_term2_2)};
#else 
               const float sfsk_sample_2{this->m_Ac*std::cos(this->m_ph0+carrier_term1_2+carrier_term2_2)};
#endif 
               this->m_sfsk_signal.m_data[i+2ull] = sfsk_sample_2;

               const float ik_3{this->m_I_ch_bitstream.m_data[i+3ull]};
               const float qk_3{this->m_Q_ch_bitstream.m_data[i+3ull]};
               const float dk_3{-ik_3*qk_3};
               const float t_i_3{gms::math::LUT_loop_indices_2257_align16[i+3ull]};
               const float psfunc_3{this->m_psfunc.m_data[i+3ull]};
               const float Phik_3{ik_3==1.0f?0.0f:C314159265358979323846264338328};
               const float carrier_term1_3{C6283185307179586476925286766559*t_i_3};
               const float carrier_term2_3{std::fma(dk_3,psfunc_3,Phik_3)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_3{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_3+carrier_term2_3)};
#else 
               const float sfsk_sample_3{this->m_Ac*std::cos(this->m_ph0+carrier_term1_3+carrier_term2_3)};
#endif 
               this->m_sfsk_signal.m_data[i+3ull] = sfsk_sample_3;

               const float ik_4{this->m_I_ch_bitstream.m_data[i+4ull]};
               const float qk_4{this->m_Q_ch_bitstream.m_data[i+4ull]};
               const float dk_4{-ik_4*qk_4};
               const float t_i_4{gms::math::LUT_loop_indices_2257_align16[i+4ull]};
               const float psfunc_4{this->m_psfunc.m_data[i+4ull]};
               const float Phik_4{ik_4==1.0f?0.0f:C314159265358979323846264338328};
               const float carrier_term1_4{C6283185307179586476925286766559*t_i_4};
               const float carrier_term2_4{std::fma(dk_4,psfunc_4,Phik_4)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_4{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_4+carrier_term2_4)};
#else 
               const float sfsk_sample_4{this->m_Ac*std::cos(this->m_ph0+carrier_term1_4+carrier_term2_4)};
#endif 
               this->m_sfsk_signal.m_data[i+4ull] = sfsk_sample_4;

               const float ik_5{this->m_I_ch_bitstream.m_data[i+5ull]};
               const float qk_5{this->m_Q_ch_bitstream.m_data[i+5ull]};
               const float dk_5{-ik_5*qk_5};
               const float t_i_5{gms::math::LUT_loop_indices_2257_align16[i+5ull]};
               const float psfunc_5{this->m_psfunc.m_data[i+5ull]};
               const float Phik_5{ik_5==1.0f?0.0f:C314159265358979323846264338328};
               const float carrier_term1_5{C6283185307179586476925286766559*t_i_5};
               const float carrier_term2_5{std::fma(dk_5,psfunc_5,Phik_5)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_5{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_5+carrier_term2_5)};
#else 
               const float sfsk_sample_5{this->m_Ac*std::cos(this->m_ph0+carrier_term1_5+carrier_term2_5)};
#endif 
               this->m_sfsk_signal.m_data[i+5ull] = sfsk_sample_5;

               const float ik_6{this->m_I_ch_bitstream.m_data[i+6ull]};
               const float qk_6{this->m_Q_ch_bitstream.m_data[i+6ull]};
               const float dk_6{-ik_6*qk_6};
               const float t_i_6{gms::math::LUT_loop_indices_2257_align16[i+6ull]};
               const float psfunc_6{this->m_psfunc.m_data[i+6ull]};
               const float Phik_6{ik_6==1.0f?0.0f:C314159265358979323846264338328};
               const float carrier_term1_6{C6283185307179586476925286766559*t_i_6};
               const float carrier_term2_6{std::fma(dk_6,psfunc_6,Phik_6)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_6{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_6+carrier_term2_6)};
#else 
               const float sfsk_sample_1{this->m_Ac*std::cos(this->m_ph0+carrier_term1_6+carrier_term2_6)};
#endif 
               this->m_sfsk_signal.m_data[i+6ull] = sfsk_sample_6;

               const float ik_7{this->m_I_ch_bitstream.m_data[i+7ull]};
               const float qk_7{this->m_Q_ch_bitstream.m_data[i+7ull]};
               const float dk_7{-ik_7*qk_7};
               const float t_i_7{gms::math::LUT_loop_indices_2257_align16[i+7ull]};
               const float psfunc_7{this->m_psfunc.m_data[i+7ull]};
               const float Phik_7{ik_7==1.0f?0.0f:C314159265358979323846264338328};
               const float carrier_term1_7{C6283185307179586476925286766559*t_i_7};
               const float carrier_term2_7{std::fma(dk_7,psfunc_7,Phik_7)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_7{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_7+carrier_term2_7)};
#else 
               const float sfsk_sample_7{this->m_Ac*std::cos(this->m_ph0+carrier_term1_7+carrier_term2_7)};
#endif 
               this->m_sfsk_signal.m_data[i+7ull] = sfsk_sample_7;
          }

          for(j = i; j != this->m_sfsk_nsamples; ++j)  
          {
               const float ik_1{this->m_I_ch_bitstream.m_data[j]};
               const float qk_1{this->m_Q_ch_bitstream.m_data[j]};
               const float dk_1{-ik_1*qk_1};
               const float t_i_1{gms::math::LUT_loop_indices_2257_align16[j]};
               const float psfunc_1{this->m_psfunc.m_data[i+1ull]};
               const float Phik_1{ik_1==1.0f?0.0f:C314159265358979323846264338328};
               const float carrier_term1_1{C6283185307179586476925286766559*t_i_1};
               const float carrier_term2_1{std::fma(dk_1,psfunc_1,Phik_1)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_1{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#else 
               const float sfsk_sample_1{this->m_Ac*std::cos(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#endif 
               this->m_sfsk_signal.m_data[j] = sfsk_sample_1; 

          }
     }
     
     return (0);

}

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_u8x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::default_scalar_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::default_scalar_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::default_scalar_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                                                   const std::int32_t,
                                                                                                                                   const std::int32_t,
                                                                                                                                   const bool,
                                                                                                                                   const bool,
                                                                                                                                   const bool,
                                                                                                                                   const bool);
                                                                                                              

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_u8x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_sse_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_sse_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_sse_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                                               const std::int32_t,
                                                                                                                               const std::int32_t,
                                                                                                                               const bool,
                                                                                                                               const bool,
                                                                                                                               const bool,
                                                                                                                               const bool);
                                                                                                             

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_u8x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_sse_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_sse_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_sse_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                                               const std::int32_t,
                                                                                                                               const std::int32_t,
                                                                                                                               const bool,
                                                                                                                               const bool,
                                                                                                                               const bool,
                                                                                                                               const bool);
                                                                                                              

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_u8x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                                               const std::int32_t,
                                                                                                                               const std::int32_t,
                                                                                                                               const bool,
                                                                                                                               const bool,
                                                                                                                               const bool,
                                                                                                                               const bool);
                                                                                                              

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_u8x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                                               const std::int32_t,
                                                                                                                               const std::int32_t,
                                                                                                                               const bool,
                                                                                                                               const bool,
                                                                                                                               const bool,
                                                                                                                               const bool);
                                                                                                              

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_u8x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx512_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx512_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx512_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                                                  const std::int32_t,
                                                                                                                                  const std::int32_t,
                                                                                                                                  const bool,
                                                                                                                                  const bool,
                                                                                                                                  const bool,
                                                                                                                                  const bool);
                                                                                                             

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_u8x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx512_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx512_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx512_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                                                  const std::int32_t,
                                                                                                                                  const std::int32_t,
                                                                                                                                  const bool,
                                                                                                                                  const bool,
                                                                                                                                  const bool,
                                                                                                                                  const bool);
                                                                                                              


namespace 
{

__ATTR_ALWAYS_INLINE__		    
static inline 
__m128
_mm_negate_ps(const __m128 v) 
{
    const __m128 nvzero{_mm_set1_ps(-0.0f)};
	return (_mm_xor_ps(v,nvzero));
}

__ATTR_ALWAYS_INLINE__		    
static inline 
__m256
_mm256_negate_ps(const __m256 v) 
{
    const __m256 nvzero{_mm256_set1_ps(-0.0f)};
	return (_mm256_xor_ps(v,nvzero));
}

__ATTR_ALWAYS_INLINE__		    
static inline 
__m512
_mm512_negate_ps(const __m512 v) 
{
    const __m512 nvzero{_mm512_set1_ps(-0.0f)};
	return (_mm512_xor_ps(v,nvzero));
}

}

template<gms::radiolocation
          ::sinusoidal_fsk_t::I_channel_bitstream_optimized_path I_ch_path,
          gms::radiolocation
          ::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path Q_ch_path,
          gms::radiolocation
          ::sinusoidal_fsk_t::pulse_shaping_function_optimized_path psf_path>
std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_sse_u4x(const iq_rectw_bitstream_vsequence_t &iq_bitstream_vsequence,
                              const std::int32_t I_ch_samples_len,
                              const std::int32_t Q_ch_samples_len,
                              const bool I_ch_do_const_prefetch,
                              const bool Q_ch_do_const_prefetch,
                              const bool psfunc_do_const_prefetch,
                              const bool do_const_prefetch,
                              const bool use_iq_bitstream_vsequence)
{
     using namespace gms::math;
     using namespace gms::common; 
     if(__builtin_expect(iq_bitstream_vsequence.m_I_nsamples!=this->m_I_ch_nsamples,0)) { return (-19);}
     if(__builtin_expect(iq_bitstream_vsequence.m_Q_nsamples!=this->m_Q_ch_nsamples,0)) { return (-20);};

     if(__builtin_expect(use_iq_bitstream_vsequence==true,1)) 
     {          
          sse_memcpy_unroll8x_ps(&this->m_I_ch_bitstream.m_data[0],
                                 &iq_bitstream_vsequence.m_I_vsequence.m_data[0],
                                 this->m_I_ch_nsamples);
          sse_memcpy_unroll8x_ps(&this->m_Q_ch_bitstream.m_data[0],
                                 &iq_bitstream_vsequence.m_Q_vsequence.m_data[0],
                                 this->m_Q_ch_nsamples);
     }
     else 
     {
     
          if constexpr (I_ch_path == I_channel_bitstream_optimized_path::default_scalar_path)
          {
                std::int32_t ret_stat_scalar = this->generate_I_channel_bitstream(I_ch_samples_len);
                if(__builtin_expect(ret_stat_scalar<0,0)) { return (-1);}
          }
          else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_sse_path)
          {
                std::int32_t ret_stat_sse    = this->generate_I_channel_bitstream_sse();
                if(__builtin_expect(ret_stat_sse<0,0))    { return (-2);}
          }
          else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_sse_u4x_path)
          {
                std::int32_t ret_stat_sse_u4x = this->generate_I_channel_bitstream_sse_u4x();
                if(__builtin_expect(ret_stat_sse_u4x<0,0)) {return (-3);}
          }
          else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_avx_path)
          {
                std::int32_t ret_stat_avx     = this->generate_I_channel_bitstream_avx(I_ch_do_const_prefetch);
                if(__builtin_expect(ret_stat_avx<0,0))     {return (-4);}
          }
          else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_avx_u4x_path)
          {
                std::int32_t ret_stat_avx_u4x = this->generate_I_channel_bitstream_avx_u4x(I_ch_do_const_prefetch);
                if(__builtin_expect(ret_stat_avx_u4x<0,0)) { return (-5);}
          }
          else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_avx512_path)
          {
                std::int32_t ret_stat_avx512  = this->generate_I_channel_bitstream_avx512(I_ch_do_const_prefetch);
                if(__builtin_expect(ret_stat_avx512<0,0))  { return (-6);}
          }
          else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_avx512_u4x_path)
          {
                std::int32_t ret_stat_avx512_u4x= this->generate_I_channel_bitstream_avx512_u4x(I_ch_do_const_prefetch);
                if(__builtin_expect(ret_stat_avx512_u4x<0,0)){ return (-7);}
          }
          else
          {
               return (-9999);
          }

          if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::default_scalar_path)
          {
               std::int32_t ret_stat_scalar = this->generate_Q_channel_bitstream(Q_ch_samples_len);
               if(__builtin_expect(ret_stat_scalar<0,0)) { return (-8);}
          }
          else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_sse_path)
          {
               std::int32_t ret_stat_sse    = this->generate_Q_channel_bitstream_sse();
               if(__builtin_expect(ret_stat_sse<0,0))    { return (-9);}
          }
          else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_sse_u4x_path)
          {
               std::int32_t ret_stat_sse_u4x = this->generate_Q_channel_bitstream_sse_u4x();
               if(__builtin_expect(ret_stat_sse_u4x<0,0)) {return (-10);}
          }
          else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_avx_path)
          {
               std::int32_t ret_stat_avx     = this->generate_Q_channel_bitstream_avx(I_ch_do_const_prefetch);
               if(__builtin_expect(ret_stat_avx<0,0))     {return (-11);}
          }
          else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_avx_u4x_path)
          {
               std::int32_t ret_stat_avx_u4x = this->generate_Q_channel_bitstream_avx_u4x(I_ch_do_const_prefetch);
               if(__builtin_expect(ret_stat_avx_u4x<0,0)) { return (-12);}
          }
          else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_avx512_path)
          {
               std::int32_t ret_stat_avx512  = this->generate_Q_channel_bitstream_avx512(I_ch_do_const_prefetch);
               if(__builtin_expect(ret_stat_avx512<0,0))  { return (-13);}
          }
          else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_avx512_u4x_path)
          {
               std::int32_t ret_stat_avx512_u4x= this->generate_Q_channel_bitstream_avx512_u4x(I_ch_do_const_prefetch);
               if(__builtin_expect(ret_stat_avx512_u4x<0,0)){ return (-14);}
          }
          else 
          {
              return (-9999);
          }

     }

     if constexpr (psf_path == pulse_shaping_function_optimized_path::default_scalar_path)
     {
          std::int32_t ret_stat_scalar = this->generate_pulse_shaping_function_u8x();
          if(__builtin_expect(ret_stat_scalar<0,0)) { return (-15);}
     }
     else if constexpr (psf_path == pulse_shaping_function_optimized_path::vector_sse_path)
     {
          std::int32_t ret_stat_sse    = this->generate_pulse_shaping_function_sse_u4x(psfunc_do_const_prefetch);
          if(__builtin_expect(ret_stat_sse<0,0))    { return (-16);}
     }
     else if constexpr (psf_path == pulse_shaping_function_optimized_path::vector_avx_path)
     {
          std::int32_t ret_stat_avx    = this->generate_pulse_shaping_function_avx_u4x(psfunc_do_const_prefetch);
          if(__builtin_expect(ret_stat_avx<0,0))    { return (-17);}
     }
     else if constexpr (psf_path == pulse_shaping_function_optimized_path::vector_avx512_path)
     {
          std::int32_t ret_stat_avx512 = this->generate_pulse_shaping_function_avx512_u4x(psfunc_do_const_prefetch);
          if(__builtin_expect(ret_stat_avx512<0,0))  { return (-18);}
     }
     else 
     {
          return (-9999);
     }

     constexpr std::size_t LUT_loop_idx_threshold{2257ull};
     const __m128 vk{_mm_set1_ps(this->m_k)};
     const __m128 vAc{_mm_set1_ps(this->m_Ac)};
     const __m128 vph0{_mm_set1_ps(this->m_ph0)};
     __m128 c0_3,c4_7;
     __m128 c8_11,c12_15;
     __m128 vpi,v2pi;
     __m128 vzero,vone;  
     std::size_t i;
     float j;
     if(!do_const_prefetch)
     {
          c0_3  = _mm_setr_ps(0.0f,1.0f,2.0f,3.0f);
          c4_7  = _mm_setr_ps(4.0f,5.0f,6.0f,7.0f);
          c8_11 = _mm_setr_ps(8.0f,9.0f,10.0f,11.0f);
          c12_15= _mm_setr_ps(12.0f,13.0f,14.0f,15.0f);
          vpi   = _mm_set1_ps(3.14159265358979323846264338328f);
          v2pi  = _mm_set1_ps(6.283185307179586476925286766559f);
          vzero = _mm_setzero_ps();
          vone  = _mm_set1_ps(1.0f);
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
                                        0.0f,0.0f,0.0f,0.0f,1.0f,1.0f,1.0f,1.0f};        
          _mm_prefetch((const char*)&prefetched_constants[0], _MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
          c0_3  = _mm_load_ps(&prefetched_constants[0]);
          c4_7  = _mm_load_ps(&prefetched_constants[4]);
          c8_11 = _mm_load_ps(&prefetched_constants[8]);
          c12_15= _mm_load_ps(&prefetched_constants[12]);
          vpi   = _mm_load_ps(&prefetched_constants[16]);
          v2pi  = _mm_load_ps(&prefetched_constants[20]);
          vzero = _mm_load_ps(&prefetched_constants[24]);
          vone  = _mm_load_ps(&prefetched_constants[28]);
     }
     if(this->m_sfsk_nsamples>LUT_loop_idx_threshold)
     {
          for(i = 0ull,j = 0.0f; (i+15ull) < this->m_sfsk_nsamples; i += 16ull,j += 16.0f) 
          {
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i],_MM_HINT_T0);
               const __m128 vt_i_0{_mm_add_ps(_mm_set1_ps(j),c0_3)};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i],_MM_HINT_T0);
               const __m128 vik_0{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m128 vqk_0{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i],_MM_HINT_T0);
               const __m128 vdk_0{_mm_mul_ps(_mm_negate_ps(vik_0),vqk_0)};
               const __m128 vpsfunc_0{_mm_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask8 vik_eq_1_0{_mm_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_0{_mm_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m128 vcarrier_term1_0{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_0))};
               const __m128 vcarrier_term2_0{_mm_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m128 vsfsk_sample_0{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0);
               //////////////////////////////////////////////////////////////////
               const __m128 vt_i_1{_mm_add_ps(_mm_set1_ps(j),c4_7)};
               const __m128 vik_1{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+4ull])};
               const __m128 vqk_1{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+4ull])};
               const __m128 vdk_1{_mm_mul_ps(_mm_negate_ps(vik_1),vqk_1)};
               const __m128 vpsfunc_1{_mm_load_ps(&this->m_psfunc.m_data[i+4ull])};
               const __mmask8 vik_eq_1_1{_mm_cmp_ps_mask(vik_1,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_1{_mm_mask_blend_ps(vik_eq_1_1,vzero,vpi)};
               const __m128 vcarrier_term1_1{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_1))};
               const __m128 vcarrier_term2_1{_mm_fmadd_ps(vdk_1,vpsfunc_1,vPhik_1)};
               const __m128 vsfsk_sample_1{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_1,vcarrier_term2_1))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+4ull],vsfsk_sample_1);
               /////////////////////////////////////////////////////////////////////
               const __m128 vt_i_2{_mm_add_ps(_mm_set1_ps(j),c8_11)};
               const __m128 vik_2{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+8ull])};
               const __m128 vqk_2{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+8ull])};
               const __m128 vdk_2{_mm_mul_ps(_mm_negate_ps(vik_2),vqk_2)};
               const __m128 vpsfunc_2{_mm_load_ps(&this->m_psfunc.m_data[i+8ull])};
               const __mmask8 vik_eq_1_2{_mm_cmp_ps_mask(vik_2,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_2{_mm_mask_blend_ps(vik_eq_1_2,vzero,vpi)};
               const __m128 vcarrier_term1_2{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_2))};
               const __m128 vcarrier_term2_2{_mm_fmadd_ps(vdk_2,vpsfunc_2,vPhik_2)};
               const __m128 vsfsk_sample_2{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_2,vcarrier_term2_2))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+8ull],vsfsk_sample_2);
               //////////////////////////////////////////////////////////////////////
               const __m128 vt_i_3{_mm_add_ps(_mm_set1_ps(j),c12_15)};
               const __m128 vik_3{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+12ull])};
               const __m128 vqk_3{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+12ull])};
               const __m128 vdk_3{_mm_mul_ps(_mm_negate_ps(vik_3),vqk_3)};
               const __m128 vpsfunc_3{_mm_load_ps(&this->m_psfunc.m_data[i+12ull])};
               const __mmask8 vik_eq_1_3{_mm_cmp_ps_mask(vik_3,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_3{_mm_mask_blend_ps(vik_eq_1_3,vzero,vpi)};
               const __m128 vcarrier_term1_3{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_3))};
               const __m128 vcarrier_term2_3{_mm_fmadd_ps(vdk_3,vpsfunc_3,vPhik_3)};
               const __m128 vsfsk_sample_3{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_3,vcarrier_term2_3))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+12ull],vsfsk_sample_3);
          }

          for(; (i+11ull) < this->m_sfsk_nsamples; i += 12ull,j += 12.0f) 
          {
               
               const __m128 vt_i_0{_mm_add_ps(_mm_set1_ps(j),c0_3)};
               const __m128 vik_0{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m128 vqk_0{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               const __m128 vdk_0{_mm_mul_ps(_mm_negate_ps(vik_0),vqk_0)};
               const __m128 vpsfunc_0{_mm_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask8 vik_eq_1_0{_mm_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_0{_mm_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m128 vcarrier_term1_0{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_0))};
               const __m128 vcarrier_term2_0{_mm_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m128 vsfsk_sample_0{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0);
               //////////////////////////////////////////////////////////////////
               const __m128 vt_i_1{_mm_add_ps(_mm_set1_ps(j),c4_7)};
               const __m128 vik_1{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+4ull])};
               const __m128 vqk_1{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+4ull])};
               const __m128 vdk_1{_mm_mul_ps(_mm_negate_ps(vik_1),vqk_1)};
               const __m128 vpsfunc_1{_mm_load_ps(&this->m_psfunc.m_data[i+4ull])};
               const __mmask8 vik_eq_1_1{_mm_cmp_ps_mask(vik_1,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_1{_mm_mask_blend_ps(vik_eq_1_1,vzero,vpi)};
               const __m128 vcarrier_term1_1{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_1))};
               const __m128 vcarrier_term2_1{_mm_fmadd_ps(vdk_1,vpsfunc_1,vPhik_1)};
               const __m128 vsfsk_sample_1{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_1,vcarrier_term2_1))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+4ull],vsfsk_sample_1);
               /////////////////////////////////////////////////////////////////////
               const __m128 vt_i_2{_mm_add_ps(_mm_set1_ps(j),c8_11)};
               const __m128 vik_2{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+8ull])};
               const __m128 vqk_2{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+8ull])};
               const __m128 vdk_2{_mm_mul_ps(_mm_negate_ps(vik_2),vqk_2)};
               const __m128 vpsfunc_2{_mm_load_ps(&this->m_psfunc.m_data[i+8ull])};
               const __mmask8 vik_eq_1_2{_mm_cmp_ps_mask(vik_2,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_2{_mm_mask_blend_ps(vik_eq_1_2,vzero,vpi)};
               const __m128 vcarrier_term1_2{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_2))};
               const __m128 vcarrier_term2_2{_mm_fmadd_ps(vdk_2,vpsfunc_2,vPhik_2)};
               const __m128 vsfsk_sample_2{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_2,vcarrier_term2_2))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+8ull],vsfsk_sample_2);
          }

          for(; (i+7ull) < this->m_sfsk_nsamples; i += 8ull,j += 8.0f) 
          {
                
               const __m128 vt_i_0{_mm_add_ps(_mm_set1_ps(j),c0_3)};
               const __m128 vik_0{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m128 vqk_0{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               const __m128 vdk_0{_mm_mul_ps(_mm_negate_ps(vik_0),vqk_0)};
               const __m128 vpsfunc_0{_mm_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask8 vik_eq_1_0{_mm_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_0{_mm_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m128 vcarrier_term1_0{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_0))};
               const __m128 vcarrier_term2_0{_mm_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m128 vsfsk_sample_0{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0);
               //////////////////////////////////////////////////////////////////
               const __m128 vt_i_1{_mm_add_ps(_mm_set1_ps(j),c4_7)};
               const __m128 vik_1{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+4ull])};
               const __m128 vqk_1{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+4ull])};
               const __m128 vdk_1{_mm_mul_ps(_mm_negate_ps(vik_1),vqk_1)};
               const __m128 vpsfunc_1{_mm_load_ps(&this->m_psfunc.m_data[i+4ull])};
               const __mmask8 vik_eq_1_1{_mm_cmp_ps_mask(vik_1,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_1{_mm_mask_blend_ps(vik_eq_1_1,vzero,vpi)};
               const __m128 vcarrier_term1_1{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_1))};
               const __m128 vcarrier_term2_1{_mm_fmadd_ps(vdk_1,vpsfunc_1,vPhik_1)};
               const __m128 vsfsk_sample_1{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_1,vcarrier_term2_1))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+4ull],vsfsk_sample_1);
          }

          for(; (i+3ull) < this->m_sfsk_nsamples; i += 4ull,j += 4.0f) 
          {
               const __m128 vt_i_0{_mm_add_ps(_mm_set1_ps(j),c0_3)};
               const __m128 vik_0{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m128 vqk_0{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               const __m128 vdk_0{_mm_mul_ps(_mm_negate_ps(vik_0),vqk_0)};
               const __m128 vpsfunc_0{_mm_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask8 vik_eq_1_0{_mm_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_0{_mm_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m128 vcarrier_term1_0{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_0))};
               const __m128 vcarrier_term2_0{_mm_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m128 vsfsk_sample_0{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0); 
          }

          for(; (i+0ull) < this->m_sfsk_nsamples; i += 1ull,j += 1.0f) 
          {
               const float ik_1{this->m_I_ch_bitstream.m_data[i]};
               const float qk_1{this->m_Q_ch_bitstream.m_data[i]};
               const float dk_1{-ik_1*qk_1};
               const float psfunc_1{this->m_psfunc.m_data[i+1ull]};
               const float Phik_1{ik_1==1.0f?0.0f:3.14159265358979323846264338328};
               const float carrier_term1_1{this->m_k*6.283185307179586476925286766559*j};
               const float carrier_term2_1{std::fma(dk_1,psfunc_1,Phik_1)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_1{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#else 
               const float sfsk_sample_1{this->m_Ac*std::cos(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#endif 
               this->m_sfsk_signal.m_data[i] = sfsk_sample_1; 
          }
     }
     else 
     {
          for(i = 0ull; (i+15ull) < this->m_sfsk_nsamples; i += 16ull) 
          {
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i],_MM_HINT_T0);
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i],_MM_HINT_T0);
               const __m128 vt_i_0{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull])};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i],_MM_HINT_T0);
               const __m128 vik_0{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m128 vqk_0{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i],_MM_HINT_T0);
               const __m128 vdk_0{_mm_mul_ps(_mm_negate_ps(vik_0),vqk_0)};
               const __m128 vpsfunc_0{_mm_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask8 vik_eq_1_0{_mm_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_0{_mm_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m128 vcarrier_term1_0{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_0))};
               const __m128 vcarrier_term2_0{_mm_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m128 vsfsk_sample_0{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0);
               //////////////////////////////////////////////////////////////////
               const __m128 vt_i_1{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull])};
               const __m128 vik_1{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+4ull])};
               const __m128 vqk_1{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+4ull])};
               const __m128 vdk_1{_mm_mul_ps(_mm_negate_ps(vik_1),vqk_1)};
               const __m128 vpsfunc_1{_mm_load_ps(&this->m_psfunc.m_data[i+4ull])};
               const __mmask8 vik_eq_1_1{_mm_cmp_ps_mask(vik_1,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_1{_mm_mask_blend_ps(vik_eq_1_1,vzero,vpi)};
               const __m128 vcarrier_term1_1{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_1))};
               const __m128 vcarrier_term2_1{_mm_fmadd_ps(vdk_1,vpsfunc_1,vPhik_1)};
               const __m128 vsfsk_sample_1{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_1,vcarrier_term2_1))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+4ull],vsfsk_sample_1);
               /////////////////////////////////////////////////////////////////////
               const __m128 vt_i_2{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull])};
               const __m128 vik_2{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+8ull])};
               const __m128 vqk_2{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+8ull])};
               const __m128 vdk_2{_mm_mul_ps(_mm_negate_ps(vik_2),vqk_2)};
               const __m128 vpsfunc_2{_mm_load_ps(&this->m_psfunc.m_data[i+8ull])};
               const __mmask8 vik_eq_1_2{_mm_cmp_ps_mask(vik_2,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_2{_mm_mask_blend_ps(vik_eq_1_2,vzero,vpi)};
               const __m128 vcarrier_term1_2{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_2))};
               const __m128 vcarrier_term2_2{_mm_fmadd_ps(vdk_2,vpsfunc_2,vPhik_2)};
               const __m128 vsfsk_sample_2{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_2,vcarrier_term2_2))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+8ull],vsfsk_sample_2);
               //////////////////////////////////////////////////////////////////////
               const __m128 vt_i_3{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+12ull])};
               const __m128 vik_3{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+12ull])};
               const __m128 vqk_3{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+12ull])};
               const __m128 vdk_3{_mm_mul_ps(_mm_negate_ps(vik_3),vqk_3)};
               const __m128 vpsfunc_3{_mm_load_ps(&this->m_psfunc.m_data[i+12ull])};
               const __mmask8 vik_eq_1_3{_mm_cmp_ps_mask(vik_3,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_3{_mm_mask_blend_ps(vik_eq_1_3,vzero,vpi)};
               const __m128 vcarrier_term1_3{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_3))};
               const __m128 vcarrier_term2_3{_mm_fmadd_ps(vdk_3,vpsfunc_3,vPhik_3)};
               const __m128 vsfsk_sample_3{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_3,vcarrier_term2_3))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+12ull],vsfsk_sample_3);
          }

          for(; (i+11ull) < this->m_sfsk_nsamples; i += 12ull) 
          {
               
               const __m128 vt_i_0{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull])};
               const __m128 vik_0{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m128 vqk_0{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               const __m128 vdk_0{_mm_mul_ps(_mm_negate_ps(vik_0),vqk_0)};
               const __m128 vpsfunc_0{_mm_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask8 vik_eq_1_0{_mm_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_0{_mm_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m128 vcarrier_term1_0{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_0))};
               const __m128 vcarrier_term2_0{_mm_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m128 vsfsk_sample_0{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0);
               //////////////////////////////////////////////////////////////////
               const __m128 vt_i_1{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull])};
               const __m128 vik_1{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+4ull])};
               const __m128 vqk_1{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+4ull])};
               const __m128 vdk_1{_mm_mul_ps(_mm_negate_ps(vik_1),vqk_1)};
               const __m128 vpsfunc_1{_mm_load_ps(&this->m_psfunc.m_data[i+4ull])};
               const __mmask8 vik_eq_1_1{_mm_cmp_ps_mask(vik_1,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_1{_mm_mask_blend_ps(vik_eq_1_1,vzero,vpi)};
               const __m128 vcarrier_term1_1{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_1))};
               const __m128 vcarrier_term2_1{_mm_fmadd_ps(vdk_1,vpsfunc_1,vPhik_1)};
               const __m128 vsfsk_sample_1{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_1,vcarrier_term2_1))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+4ull],vsfsk_sample_1);
               /////////////////////////////////////////////////////////////////////
               const __m128 vt_i_2{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull])};
               const __m128 vik_2{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+8ull])};
               const __m128 vqk_2{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+8ull])};
               const __m128 vdk_2{_mm_mul_ps(_mm_negate_ps(vik_2),vqk_2)};
               const __m128 vpsfunc_2{_mm_load_ps(&this->m_psfunc.m_data[i+8ull])};
               const __mmask8 vik_eq_1_2{_mm_cmp_ps_mask(vik_2,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_2{_mm_mask_blend_ps(vik_eq_1_2,vzero,vpi)};
               const __m128 vcarrier_term1_2{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_2))};
               const __m128 vcarrier_term2_2{_mm_fmadd_ps(vdk_2,vpsfunc_2,vPhik_2)};
               const __m128 vsfsk_sample_2{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_2,vcarrier_term2_2))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+8ull],vsfsk_sample_2);
          }

          for(; (i+7ull) < this->m_sfsk_nsamples; i += 8ull) 
          {
                
               const __m128 vt_i_0{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull])};
               const __m128 vik_0{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m128 vqk_0{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               const __m128 vdk_0{_mm_mul_ps(_mm_negate_ps(vik_0),vqk_0)};
               const __m128 vpsfunc_0{_mm_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask8 vik_eq_1_0{_mm_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_0{_mm_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m128 vcarrier_term1_0{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_0))};
               const __m128 vcarrier_term2_0{_mm_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m128 vsfsk_sample_0{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0);
               //////////////////////////////////////////////////////////////////
               const __m128 vt_i_1{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull])};
               const __m128 vik_1{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+4ull])};
               const __m128 vqk_1{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+4ull])};
               const __m128 vdk_1{_mm_mul_ps(_mm_negate_ps(vik_1),vqk_1)};
               const __m128 vpsfunc_1{_mm_load_ps(&this->m_psfunc.m_data[i+4ull])};
               const __mmask8 vik_eq_1_1{_mm_cmp_ps_mask(vik_1,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_1{_mm_mask_blend_ps(vik_eq_1_1,vzero,vpi)};
               const __m128 vcarrier_term1_1{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_1))};
               const __m128 vcarrier_term2_1{_mm_fmadd_ps(vdk_1,vpsfunc_1,vPhik_1)};
               const __m128 vsfsk_sample_1{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_1,vcarrier_term2_1))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+4ull],vsfsk_sample_1);
          }

          for(; (i+3ull) < this->m_sfsk_nsamples; i += 4ull) 
          {
               const __m128 vt_i_0{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull])};
               const __m128 vik_0{_mm_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m128 vqk_0{_mm_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               const __m128 vdk_0{_mm_mul_ps(_mm_negate_ps(vik_0),vqk_0)};
               const __m128 vpsfunc_0{_mm_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask8 vik_eq_1_0{_mm_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m128 vPhik_0{_mm_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m128 vcarrier_term1_0{_mm_mul_ps(vk,_mm_mul_ps(v2pi,vt_i_0))};
               const __m128 vcarrier_term2_0{_mm_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m128 vsfsk_sample_0{_mm_mul_ps(vAc,_mm_sin_ps(_mm_add_ps(vph0,_mm_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0); 
          }

          for(; (i+0ull) < this->m_sfsk_nsamples; i += 1ull) 
          {
               const float t_i{gms::math::LUT_loop_indices_2257_align16[i]};
               const float ik_1{this->m_I_ch_bitstream.m_data[i]};
               const float qk_1{this->m_Q_ch_bitstream.m_data[i]};
               const float dk_1{-ik_1*qk_1};
               const float psfunc_1{this->m_psfunc.m_data[i+1ull]};
               const float Phik_1{ik_1==1.0f?0.0f:3.14159265358979323846264338328};
               const float carrier_term1_1{this->m_k*6.283185307179586476925286766559*t_i};
               const float carrier_term2_1{std::fma(dk_1,psfunc_1,Phik_1)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_1{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#else 
               const float sfsk_sample_1{this->m_Ac*std::cos(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#endif 
               this->m_sfsk_signal.m_data[i] = sfsk_sample_1; 
          }
     }

     return (0);
}

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_sse_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::default_scalar_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::default_scalar_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::default_scalar_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_sse_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_sse_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_sse_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_sse_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_sse_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_sse_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_sse_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_sse_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_sse_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_sse_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_sse_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx512_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx512_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx512_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                             const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_sse_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx512_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx512_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx512_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);
template<gms::radiolocation
          ::sinusoidal_fsk_t::I_channel_bitstream_optimized_path I_ch_path,
          gms::radiolocation
          ::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path Q_ch_path,
          gms::radiolocation
          ::sinusoidal_fsk_t::pulse_shaping_function_optimized_path psf_path>
std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_avx_u4x(const iq_rectw_bitstream_vsequence_t &iq_bitstream_vsequence,
                              const std::int32_t I_ch_samples_len,
                              const std::int32_t Q_ch_samples_len,
                              const bool I_ch_do_const_prefetch,
                              const bool Q_ch_do_const_prefetch,
                              const bool psfunc_do_const_prefetch,
                              const bool do_const_prefetch,
                              const bool use_iq_bitstream_vsequence)
{
     using namespace gms::math;
     using namespace gms::common;
     if(__builtin_expect(iq_bitstream_vsequence.m_I_nsamples!=this->m_I_ch_nsamples,0)) { return (-19);}
     if(__builtin_expect(iq_bitstream_vsequence.m_Q_nsamples!=this->m_Q_ch_nsamples,0)) { return (-20);};

     if(__builtin_expect(use_iq_bitstream_vsequence==true,1)) 
     {          
          avx_memcpy_unroll8x_ps(&this->m_I_ch_bitstream.m_data[0],
                                 &iq_bitstream_vsequence.m_I_vsequence.m_data[0],
                                 this->m_I_ch_nsamples);
          avx_memcpy_unroll8x_ps(&this->m_Q_ch_bitstream.m_data[0],
                                 &iq_bitstream_vsequence.m_Q_vsequence.m_data[0],
                                 this->m_Q_ch_nsamples);
     }
     else 
     {
     
           if constexpr (I_ch_path == I_channel_bitstream_optimized_path::default_scalar_path)
           {
                 std::int32_t ret_stat_scalar = this->generate_I_channel_bitstream(I_ch_samples_len);
                 if(__builtin_expect(ret_stat_scalar<0,0)) { return (-1);}
           }
           else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_sse_path)
           {
                 std::int32_t ret_stat_sse    = this->generate_I_channel_bitstream_sse();
                 if(__builtin_expect(ret_stat_sse<0,0))    { return (-2);}
           }
           else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_sse_u4x_path)
           {
                 std::int32_t ret_stat_sse_u4x = this->generate_I_channel_bitstream_sse_u4x();
                 if(__builtin_expect(ret_stat_sse_u4x<0,0)) {return (-3);}
           }
           else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_avx_path)
           {
                 std::int32_t ret_stat_avx     = this->generate_I_channel_bitstream_avx(I_ch_do_const_prefetch);
                 if(__builtin_expect(ret_stat_avx<0,0))     {return (-4);}
           }
           else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_avx_u4x_path)
           {
                 std::int32_t ret_stat_avx_u4x = this->generate_I_channel_bitstream_avx_u4x(I_ch_do_const_prefetch);
                 if(__builtin_expect(ret_stat_avx_u4x<0,0)) { return (-5);}
           }
           else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_avx512_path)
           {
                 std::int32_t ret_stat_avx512  = this->generate_I_channel_bitstream_avx512(I_ch_do_const_prefetch);
                 if(__builtin_expect(ret_stat_avx512<0,0))  { return (-6);}
           }
           else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_avx512_u4x_path)
           {
                 std::int32_t ret_stat_avx512_u4x= this->generate_I_channel_bitstream_avx512_u4x(I_ch_do_const_prefetch);
                 if(__builtin_expect(ret_stat_avx512_u4x<0,0)){ return (-7);}
           }
           else
           {
                 return (-9999);
           }

           if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::default_scalar_path)
           {
                  std::int32_t ret_stat_scalar = this->generate_Q_channel_bitstream(Q_ch_samples_len);
                  if(__builtin_expect(ret_stat_scalar<0,0)) { return (-8);}
           }
           else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_sse_path)
           {
                  std::int32_t ret_stat_sse    = this->generate_Q_channel_bitstream_sse();
                  if(__builtin_expect(ret_stat_sse<0,0))    { return (-9);}
           }
           else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_sse_u4x_path)
           {
                  std::int32_t ret_stat_sse_u4x = this->generate_Q_channel_bitstream_sse_u4x();
                  if(__builtin_expect(ret_stat_sse_u4x<0,0)) {return (-10);}
           }
           else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_avx_path)
           {
                 std::int32_t ret_stat_avx     = this->generate_Q_channel_bitstream_avx(I_ch_do_const_prefetch);
                 if(__builtin_expect(ret_stat_avx<0,0))     {return (-11);}
           }
           else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_avx_u4x_path)
           {
                 std::int32_t ret_stat_avx_u4x = this->generate_Q_channel_bitstream_avx_u4x(I_ch_do_const_prefetch);
                 if(__builtin_expect(ret_stat_avx_u4x<0,0)) { return (-12);}
           }
           else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_avx512_path)
           {
                 std::int32_t ret_stat_avx512  = this->generate_Q_channel_bitstream_avx512(I_ch_do_const_prefetch);
                 if(__builtin_expect(ret_stat_avx512<0,0))  { return (-13);}
           }
           else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_avx512_u4x_path)
           {
                 std::int32_t ret_stat_avx512_u4x= this->generate_Q_channel_bitstream_avx512_u4x(I_ch_do_const_prefetch);
                 if(__builtin_expect(ret_stat_avx512_u4x<0,0)){ return (-14);}
           }
           else 
           {
                return (-9999);
           }
     
     }

     if constexpr (psf_path == pulse_shaping_function_optimized_path::default_scalar_path)
     {
          std::int32_t ret_stat_scalar = this->generate_pulse_shaping_function_u8x();
          if(__builtin_expect(ret_stat_scalar<0,0)) { return (-15);}
     }
     else if constexpr (psf_path == pulse_shaping_function_optimized_path::vector_sse_path)
     {
          std::int32_t ret_stat_sse    = this->generate_pulse_shaping_function_sse_u4x(psfunc_do_const_prefetch);
          if(__builtin_expect(ret_stat_sse<0,0))    { return (-16);}
     }
     else if constexpr (psf_path == pulse_shaping_function_optimized_path::vector_avx_path)
     {
          std::int32_t ret_stat_avx    = this->generate_pulse_shaping_function_avx_u4x(psfunc_do_const_prefetch);
          if(__builtin_expect(ret_stat_avx<0,0))    { return (-17);}
     }
     else if constexpr (psf_path == pulse_shaping_function_optimized_path::vector_avx512_path)
     {
          std::int32_t ret_stat_avx512 = this->generate_pulse_shaping_function_avx512_u4x(psfunc_do_const_prefetch);
          if(__builtin_expect(ret_stat_avx512<0,0))  { return (-18);}
     }
     else 
     {
          return (-9999);
     }
     constexpr std::size_t LUT_loop_idx_threshold{2257ull};
     const __m256 vk{_mm256_set1_ps(this->m_k)};
     const __m256 vAc{_mm256_set1_ps(this->m_Ac)};
     const __m256 vph0{_mm256_set1_ps(this->m_ph0)};
     __m256 c0_7,c8_15;
     __m256 c16_23,c24_31;
     __m256 vpi,v2pi;
     __m256 vzero,vone;  
     std::size_t i;
     float j;
     if(!do_const_prefetch)
     {
          
          c0_7   = _mm256_setr_ps(0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f);
          c8_15  = _mm256_setr_ps(8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f);
          c16_23 = _mm256_setr_ps(16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f);
          c24_31 = _mm256_setr_ps(24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f);
          vpi    = _mm256_set1_ps(3.14159265358979323846264338328f);
          v2pi   = _mm256_set1_ps(6.283185307179586476925286766559f); 
          vzero  = _mm256_setzero_ps();
          vone   = _mm256_set1_ps(1.0f);
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
                                    0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
                                    1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f};
          _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[32],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[48],_MM_HINT_T0);       
          _mm_prefetch((const char*)&prefetched_constants[0], _MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
          c0_7   = _mm256_load_ps(&prefetched_constants[0]);
          c8_15  = _mm256_load_ps(&prefetched_constants[8]);
          c16_23 = _mm256_load_ps(&prefetched_constants[16]);
          c24_31 = _mm256_load_ps(&prefetched_constants[24]);
          vpi    = _mm256_load_ps(&prefetched_constants[32]);
          v2pi   = _mm256_load_ps(&prefetched_constants[40]);
          vzero  = _mm256_load_ps(&prefetched_constants[48]);
          vone   = _mm256_load_ps(&prefetched_constants[56]);
     }
     if(this->m_sfsk_nsamples>LUT_loop_idx_threshold)
     {
          for(i = 0ull,j = 0.0f; (i+31ull) < this->m_sfsk_nsamples; i += 32ull,j += 32.0f) 
          {
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m256 vt_i_0{_mm256_add_ps(_mm256_set1_ps(j),c0_7)};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m256 vik_0{_mm256_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m256 vqk_0{_mm256_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i+0ull],_MM_HINT_T0);
               const __m256 vdk_0{_mm256_mul_ps(_mm256_negate_ps(vik_0),vqk_0)};
               const __m256 vpsfunc_0{_mm256_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask8 vik_eq_1_0{_mm256_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m256 vPhik_0{_mm256_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m256 vcarrier_term1_0{_mm256_mul_ps(vk,_mm256_mul_ps(v2pi,vt_i_0))};
               const __m256 vcarrier_term2_0{_mm256_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m256 vsfsk_sample_0{_mm256_mul_ps(vAc,_mm256_sin_ps(_mm256_add_ps(vph0,_mm256_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm256_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0);
               //////////////////////////////////////////////////////////////////
               const __m256 vt_i_1{_mm256_add_ps(_mm256_set1_ps(j),c8_15)};
               const __m256 vik_1{_mm256_load_ps(&this->m_I_ch_bitstream.m_data[i+8ull])};
               const __m256 vqk_1{_mm256_load_ps(&this->m_Q_ch_bitstream.m_data[i+8ull])};
               const __m256 vdk_1{_mm256_mul_ps(_mm256_negate_ps(vik_1),vqk_1)};
               const __m256 vpsfunc_1{_mm256_load_ps(&this->m_psfunc.m_data[i+8ull])};
               const __mmask8 vik_eq_1_1{_mm256_cmp_ps_mask(vik_1,vone,_CMP_EQ_OQ)};
               const __m256 vPhik_1{_mm256_mask_blend_ps(vik_eq_1_1,vzero,vpi)};
               const __m256 vcarrier_term1_1{_mm256_mul_ps(vk,_mm256_mul_ps(v2pi,vt_i_1))};
               const __m256 vcarrier_term2_1{_mm256_fmadd_ps(vdk_1,vpsfunc_1,vPhik_1)};
               const __m256 vsfsk_sample_1{_mm256_mul_ps(vAc,_mm256_sin_ps(_mm256_add_ps(vph0,_mm256_add_ps(vcarrier_term1_1,vcarrier_term2_1))))};
               _mm256_store_ps(&this->m_sfsk_signal.m_data[i+8ull],vsfsk_sample_1);
               /////////////////////////////////////////////////////////////////////
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T0);
               const __m256 vt_i_2{_mm256_add_ps(_mm256_set1_ps(j),c16_23)};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T0);
               const __m256 vik_2{_mm256_load_ps(&this->m_I_ch_bitstream.m_data[i+16ull])};
               const __m256 vqk_2{_mm256_load_ps(&this->m_Q_ch_bitstream.m_data[i+16ull])};
                _mm_prefetch((const char*)&this->m_psfunc.m_data[i+16ull],_MM_HINT_T0);
               const __m256 vdk_2{_mm256_mul_ps(_mm256_negate_ps(vik_2),vqk_2)};
               const __m256 vpsfunc_2{_mm256_load_ps(&this->m_psfunc.m_data[i+16ull])};
               const __mmask8 vik_eq_1_2{_mm256_cmp_ps_mask(vik_2,vone,_CMP_EQ_OQ)};
               const __m256 vPhik_2{_mm256_mask_blend_ps(vik_eq_1_2,vzero,vpi)};
               const __m256 vcarrier_term1_2{_mm256_mul_ps(vk,_mm256_mul_ps(v2pi,vt_i_2))};
               const __m256 vcarrier_term2_2{_mm256_fmadd_ps(vdk_2,vpsfunc_2,vPhik_2)};
               const __m256 vsfsk_sample_2{_mm256_mul_ps(vAc,_mm256_sin_ps(_mm256_add_ps(vph0,_mm256_add_ps(vcarrier_term1_2,vcarrier_term2_2))))};
               _mm256_store_ps(&this->m_sfsk_signal.m_data[i+16ull],vsfsk_sample_2);
               //////////////////////////////////////////////////////////////////////
               const __m256 vt_i_3{_mm256_add_ps(_mm256_set1_ps(j),c24_31)};
               const __m256 vik_3{_mm256_load_ps(&this->m_I_ch_bitstream.m_data[i+24ull])};
               const __m256 vqk_3{_mm256_load_ps(&this->m_Q_ch_bitstream.m_data[i+24ull])};
               const __m256 vdk_3{_mm256_mul_ps(_mm256_negate_ps(vik_3),vqk_3)};
               const __m256 vpsfunc_3{_mm256_load_ps(&this->m_psfunc.m_data[i+24ull])};
               const __mmask8 vik_eq_1_3{_mm256_cmp_ps_mask(vik_3,vone,_CMP_EQ_OQ)};
               const __m256 vPhik_3{_mm256_mask_blend_ps(vik_eq_1_3,vzero,vpi)};
               const __m256 vcarrier_term1_3{_mm256_mul_ps(vk,_mm256_mul_ps(v2pi,vt_i_3))};
               const __m256 vcarrier_term2_3{_mm256_fmadd_ps(vdk_3,vpsfunc_3,vPhik_3)};
               const __m256 vsfsk_sample_3{_mm256_mul_ps(vAc,_mm256_sin_ps(_mm256_add_ps(vph0,_mm256_add_ps(vcarrier_term1_3,vcarrier_term2_3))))};
               _mm256_store_ps(&this->m_sfsk_signal.m_data[i+24ull],vsfsk_sample_3);
          }

          for(; (i+15ull) < this->m_sfsk_nsamples; i += 16ull,j += 16.0f) 
          {
               
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m256 vt_i_0{_mm256_add_ps(_mm256_set1_ps(j),c0_7)};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m256 vik_0{_mm256_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m256 vqk_0{_mm256_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i+0ull],_MM_HINT_T0);
               const __m256 vdk_0{_mm256_mul_ps(_mm256_negate_ps(vik_0),vqk_0)};
               const __m256 vpsfunc_0{_mm256_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask8 vik_eq_1_0{_mm256_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m256 vPhik_0{_mm256_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m256 vcarrier_term1_0{_mm256_mul_ps(vk,_mm256_mul_ps(v2pi,vt_i_0))};
               const __m256 vcarrier_term2_0{_mm256_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m256 vsfsk_sample_0{_mm256_mul_ps(vAc,_mm256_sin_ps(_mm256_add_ps(vph0,_mm256_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm256_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0);
               //////////////////////////////////////////////////////////////////
               const __m256 vt_i_1{_mm256_add_ps(_mm256_set1_ps(j),c8_15)};
               const __m256 vik_1{_mm256_load_ps(&this->m_I_ch_bitstream.m_data[i+8ull])};
               const __m256 vqk_1{_mm256_load_ps(&this->m_Q_ch_bitstream.m_data[i+8ull])};
               const __m256 vdk_1{_mm256_mul_ps(_mm256_negate_ps(vik_1),vqk_1)};
               const __m256 vpsfunc_1{_mm256_load_ps(&this->m_psfunc.m_data[i+8ull])};
               const __mmask8 vik_eq_1_1{_mm256_cmp_ps_mask(vik_1,vone,_CMP_EQ_OQ)};
               const __m256 vPhik_1{_mm256_mask_blend_ps(vik_eq_1_1,vzero,vpi)};
               const __m256 vcarrier_term1_1{_mm256_mul_ps(vk,_mm256_mul_ps(v2pi,vt_i_1))};
               const __m256 vcarrier_term2_1{_mm256_fmadd_ps(vdk_1,vpsfunc_1,vPhik_1)};
               const __m256 vsfsk_sample_1{_mm256_mul_ps(vAc,_mm256_sin_ps(_mm256_add_ps(vph0,_mm256_add_ps(vcarrier_term1_1,vcarrier_term2_1))))};
               _mm256_store_ps(&this->m_sfsk_signal.m_data[i+8ull],vsfsk_sample_1);
          }

          for(; (i+7ull) < this->m_sfsk_nsamples; i += 8ull,j += 8.0f) 
          {
                
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m256 vt_i_0{_mm256_add_ps(_mm256_set1_ps(j),c0_7)};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m256 vik_0{_mm256_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m256 vqk_0{_mm256_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i+0ull],_MM_HINT_T0);
               const __m256 vdk_0{_mm256_mul_ps(_mm256_negate_ps(vik_0),vqk_0)};
               const __m256 vpsfunc_0{_mm256_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask8 vik_eq_1_0{_mm256_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m256 vPhik_0{_mm256_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m256 vcarrier_term1_0{_mm256_mul_ps(vk,_mm256_mul_ps(v2pi,vt_i_0))};
               const __m256 vcarrier_term2_0{_mm256_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m256 vsfsk_sample_0{_mm256_mul_ps(vAc,_mm256_sin_ps(_mm256_add_ps(vph0,_mm256_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm256_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0);
          }

          for(; (i+0ull) < this->m_sfsk_nsamples; i += 1ull,j += 1.0f) 
          {
               const float ik_1{this->m_I_ch_bitstream.m_data[i]};
               const float qk_1{this->m_Q_ch_bitstream.m_data[i]};
               const float dk_1{-ik_1*qk_1};
               const float psfunc_1{this->m_psfunc.m_data[i+1ull]};
               const float Phik_1{ik_1==1.0f?0.0f:3.14159265358979323846264338328};
               const float carrier_term1_1{this->m_k*6.283185307179586476925286766559*j};
               const float carrier_term2_1{std::fma(dk_1,psfunc_1,Phik_1)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_1{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#else 
               const float sfsk_sample_1{this->m_Ac*std::cos(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#endif 
               this->m_sfsk_signal.m_data[i] = sfsk_sample_1; 
          }
     }
     else 
     {
          for(i = 0ull; (i+31ull) < this->m_sfsk_nsamples; i += 32ull) 
          {
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align32[i+0ull],_MM_HINT_T0);
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m256 vt_i_0{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+0ull])};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m256 vik_0{_mm256_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m256 vqk_0{_mm256_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i+0ull],_MM_HINT_T0);
               const __m256 vdk_0{_mm256_mul_ps(_mm256_negate_ps(vik_0),vqk_0)};
               const __m256 vpsfunc_0{_mm256_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask8 vik_eq_1_0{_mm256_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m256 vPhik_0{_mm256_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m256 vcarrier_term1_0{_mm256_mul_ps(vk,_mm256_mul_ps(v2pi,vt_i_0))};
               const __m256 vcarrier_term2_0{_mm256_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m256 vsfsk_sample_0{_mm256_mul_ps(vAc,_mm256_sin_ps(_mm256_add_ps(vph0,_mm256_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm256_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0);
               //////////////////////////////////////////////////////////////////
               const __m256 vt_i_1{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+8ull])};
               const __m256 vik_1{_mm256_load_ps(&this->m_I_ch_bitstream.m_data[i+8ull])};
               const __m256 vqk_1{_mm256_load_ps(&this->m_Q_ch_bitstream.m_data[i+8ull])};
               const __m256 vdk_1{_mm256_mul_ps(_mm256_negate_ps(vik_1),vqk_1)};
               const __m256 vpsfunc_1{_mm256_load_ps(&this->m_psfunc.m_data[i+8ull])};
               const __mmask8 vik_eq_1_1{_mm256_cmp_ps_mask(vik_1,vone,_CMP_EQ_OQ)};
               const __m256 vPhik_1{_mm256_mask_blend_ps(vik_eq_1_1,vzero,vpi)};
               const __m256 vcarrier_term1_1{_mm256_mul_ps(vk,_mm256_mul_ps(v2pi,vt_i_1))};
               const __m256 vcarrier_term2_1{_mm256_fmadd_ps(vdk_1,vpsfunc_1,vPhik_1)};
               const __m256 vsfsk_sample_1{_mm256_mul_ps(vAc,_mm256_sin_ps(_mm256_add_ps(vph0,_mm256_add_ps(vcarrier_term1_1,vcarrier_term2_1))))};
               _mm256_store_ps(&this->m_sfsk_signal.m_data[i+8ull],vsfsk_sample_1);
               /////////////////////////////////////////////////////////////////////
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align32[i+16ull],_MM_HINT_T0);
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T0);
               const __m256 vt_i_2{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+16ull])};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T0);
               const __m256 vik_2{_mm256_load_ps(&this->m_I_ch_bitstream.m_data[i+16ull])};
               const __m256 vqk_2{_mm256_load_ps(&this->m_Q_ch_bitstream.m_data[i+16ull])};
                _mm_prefetch((const char*)&this->m_psfunc.m_data[i+16ull],_MM_HINT_T0);
               const __m256 vdk_2{_mm256_mul_ps(_mm256_negate_ps(vik_2),vqk_2)};
               const __m256 vpsfunc_2{_mm256_load_ps(&this->m_psfunc.m_data[i+16ull])};
               const __mmask8 vik_eq_1_2{_mm256_cmp_ps_mask(vik_2,vone,_CMP_EQ_OQ)};
               const __m256 vPhik_2{_mm256_mask_blend_ps(vik_eq_1_2,vzero,vpi)};
               const __m256 vcarrier_term1_2{_mm256_mul_ps(vk,_mm256_mul_ps(v2pi,vt_i_2))};
               const __m256 vcarrier_term2_2{_mm256_fmadd_ps(vdk_2,vpsfunc_2,vPhik_2)};
               const __m256 vsfsk_sample_2{_mm256_mul_ps(vAc,_mm256_sin_ps(_mm256_add_ps(vph0,_mm256_add_ps(vcarrier_term1_2,vcarrier_term2_2))))};
               _mm256_store_ps(&this->m_sfsk_signal.m_data[i+16ull],vsfsk_sample_2);
               //////////////////////////////////////////////////////////////////////
               const __m256 vt_i_3{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+24ull])};
               const __m256 vik_3{_mm256_load_ps(&this->m_I_ch_bitstream.m_data[i+24ull])};
               const __m256 vqk_3{_mm256_load_ps(&this->m_Q_ch_bitstream.m_data[i+24ull])};
               const __m256 vdk_3{_mm256_mul_ps(_mm256_negate_ps(vik_3),vqk_3)};
               const __m256 vpsfunc_3{_mm256_load_ps(&this->m_psfunc.m_data[i+24ull])};
               const __mmask8 vik_eq_1_3{_mm256_cmp_ps_mask(vik_3,vone,_CMP_EQ_OQ)};
               const __m256 vPhik_3{_mm256_mask_blend_ps(vik_eq_1_3,vzero,vpi)};
               const __m256 vcarrier_term1_3{_mm256_mul_ps(vk,_mm256_mul_ps(v2pi,vt_i_3))};
               const __m256 vcarrier_term2_3{_mm256_fmadd_ps(vdk_3,vpsfunc_3,vPhik_3)};
               const __m256 vsfsk_sample_3{_mm256_mul_ps(vAc,_mm256_sin_ps(_mm256_add_ps(vph0,_mm256_add_ps(vcarrier_term1_3,vcarrier_term2_3))))};
               _mm256_store_ps(&this->m_sfsk_signal.m_data[i+24ull],vsfsk_sample_3);
          }

          for(; (i+15ull) < this->m_sfsk_nsamples; i += 16ull) 
          {
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align32[i+0ull],_MM_HINT_T0);
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m256 vt_i_0{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+0ull])};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m256 vik_0{_mm256_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m256 vqk_0{_mm256_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i+0ull],_MM_HINT_T0);
               const __m256 vdk_0{_mm256_mul_ps(_mm256_negate_ps(vik_0),vqk_0)};
               const __m256 vpsfunc_0{_mm256_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask8 vik_eq_1_0{_mm256_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m256 vPhik_0{_mm256_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m256 vcarrier_term1_0{_mm256_mul_ps(vk,_mm256_mul_ps(v2pi,vt_i_0))};
               const __m256 vcarrier_term2_0{_mm256_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m256 vsfsk_sample_0{_mm256_mul_ps(vAc,_mm256_sin_ps(_mm256_add_ps(vph0,_mm256_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm256_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0);
               //////////////////////////////////////////////////////////////////
               const __m256 vt_i_1{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+8ull])};
               const __m256 vik_1{_mm256_load_ps(&this->m_I_ch_bitstream.m_data[i+8ull])};
               const __m256 vqk_1{_mm256_load_ps(&this->m_Q_ch_bitstream.m_data[i+8ull])};
               const __m256 vdk_1{_mm256_mul_ps(_mm256_negate_ps(vik_1),vqk_1)};
               const __m256 vpsfunc_1{_mm256_load_ps(&this->m_psfunc.m_data[i+8ull])};
               const __mmask8 vik_eq_1_1{_mm256_cmp_ps_mask(vik_1,vone,_CMP_EQ_OQ)};
               const __m256 vPhik_1{_mm256_mask_blend_ps(vik_eq_1_1,vzero,vpi)};
               const __m256 vcarrier_term1_1{_mm256_mul_ps(vk,_mm256_mul_ps(v2pi,vt_i_1))};
               const __m256 vcarrier_term2_1{_mm256_fmadd_ps(vdk_1,vpsfunc_1,vPhik_1)};
               const __m256 vsfsk_sample_1{_mm256_mul_ps(vAc,_mm256_sin_ps(_mm256_add_ps(vph0,_mm256_add_ps(vcarrier_term1_1,vcarrier_term2_1))))};
               _mm256_store_ps(&this->m_sfsk_signal.m_data[i+8ull],vsfsk_sample_1);
              
          }

          for(; (i+7ull) < this->m_sfsk_nsamples; i += 8ull) 
          {
                
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align32[i+0ull],_MM_HINT_T0);
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m256 vt_i_0{_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+0ull])};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m256 vik_0{_mm256_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m256 vqk_0{_mm256_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i+0ull],_MM_HINT_T0);
               const __m256 vdk_0{_mm256_mul_ps(_mm256_negate_ps(vik_0),vqk_0)};
               const __m256 vpsfunc_0{_mm256_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask8 vik_eq_1_0{_mm256_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m256 vPhik_0{_mm256_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m256 vcarrier_term1_0{_mm256_mul_ps(vk,_mm256_mul_ps(v2pi,vt_i_0))};
               const __m256 vcarrier_term2_0{_mm256_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m256 vsfsk_sample_0{_mm256_mul_ps(vAc,_mm256_sin_ps(_mm256_add_ps(vph0,_mm256_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm256_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0);
          }

          for(; (i+0ull) < this->m_sfsk_nsamples; i += 1ull) 
          {
               const float t_i{gms::math::LUT_loop_indices_2257_align32[i]};
               const float ik_1{this->m_I_ch_bitstream.m_data[i]};
               const float qk_1{this->m_Q_ch_bitstream.m_data[i]};
               const float dk_1{-ik_1*qk_1};
               const float psfunc_1{this->m_psfunc.m_data[i+1ull]};
               const float Phik_1{ik_1==1.0f?0.0f:3.14159265358979323846264338328};
               const float carrier_term1_1{this->m_k*6.283185307179586476925286766559*t_i};
               const float carrier_term2_1{std::fma(dk_1,psfunc_1,Phik_1)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_1{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#else 
               const float sfsk_sample_1{this->m_Ac*std::cos(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#endif 
               this->m_sfsk_signal.m_data[i] = sfsk_sample_1; 
          }
     }

     return (0);
}

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_avx512_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::default_scalar_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::default_scalar_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::default_scalar_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_avx512_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_sse_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_sse_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_sse_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_avx512_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_sse_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_sse_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_sse_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_avx_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_avx_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_avx_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx512_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx512_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx512_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_avx_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx512_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx512_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx512_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);


 template<gms::radiolocation
          ::sinusoidal_fsk_t::I_channel_bitstream_optimized_path I_ch_path,
          gms::radiolocation
          ::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path Q_ch_path,
          gms::radiolocation
          ::sinusoidal_fsk_t::pulse_shaping_function_optimized_path psf_path>
std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_avx512_u4x(const iq_rectw_bitstream_vsequence_t &iq_bitstream_vsequence,
                                 const std::int32_t I_ch_samples_len,
                                 const std::int32_t Q_ch_samples_len,
                                 const bool I_ch_do_const_prefetch,
                                 const bool Q_ch_do_const_prefetch,
                                 const bool psfunc_do_const_prefetch,
                                 const bool do_const_prefetch,
                                 const bool use_iq_bitstream_vsequence)
{
     using namespace gms::math;
     using namespace gms::common;
     if(__builtin_expect(iq_bitstream_vsequence.m_I_nsamples!=this->m_I_ch_nsamples,0)) { return (-19);}
     if(__builtin_expect(iq_bitstream_vsequence.m_Q_nsamples!=this->m_Q_ch_nsamples,0)) { return (-20);};

     if(__builtin_expect(use_iq_bitstream_vsequence==true,1)) 
     {          
          avx512_memcpy_unroll8x_ps(&this->m_I_ch_bitstream.m_data[0],
                                 &iq_bitstream_vsequence.m_I_vsequence.m_data[0],
                                 this->m_I_ch_nsamples);
          avx512_memcpy_unroll8x_ps(&this->m_Q_ch_bitstream.m_data[0],
                                 &iq_bitstream_vsequence.m_Q_vsequence.m_data[0],
                                 this->m_Q_ch_nsamples);
     }
     else 
     {
     
          if constexpr (I_ch_path == I_channel_bitstream_optimized_path::default_scalar_path)
          {
               std::int32_t ret_stat_scalar = this->generate_I_channel_bitstream(I_ch_samples_len);
               if(__builtin_expect(ret_stat_scalar<0,0)) { return (-1);}
          }
          else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_sse_path)
          {
               std::int32_t ret_stat_sse    = this->generate_I_channel_bitstream_sse();
               if(__builtin_expect(ret_stat_sse<0,0))    { return (-2);}
          }
          else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_sse_u4x_path)
          {
               std::int32_t ret_stat_sse_u4x = this->generate_I_channel_bitstream_sse_u4x();
               if(__builtin_expect(ret_stat_sse_u4x<0,0)) {return (-3);}
          }
          else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_avx_path)
          {
               std::int32_t ret_stat_avx     = this->generate_I_channel_bitstream_avx(I_ch_do_const_prefetch);
               if(__builtin_expect(ret_stat_avx<0,0))     {return (-4);}
          }
          else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_avx_u4x_path)
          {
               std::int32_t ret_stat_avx_u4x = this->generate_I_channel_bitstream_avx_u4x(I_ch_do_const_prefetch);
               if(__builtin_expect(ret_stat_avx_u4x<0,0)) { return (-5);}
          }
          else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_avx512_path)
          {
               std::int32_t ret_stat_avx512  = this->generate_I_channel_bitstream_avx512(I_ch_do_const_prefetch);
               if(__builtin_expect(ret_stat_avx512<0,0))  { return (-6);}
          }
          else if constexpr (I_ch_path == I_channel_bitstream_optimized_path::vector_avx512_u4x_path)
          {
               std::int32_t ret_stat_avx512_u4x= this->generate_I_channel_bitstream_avx512_u4x(I_ch_do_const_prefetch);
               if(__builtin_expect(ret_stat_avx512_u4x<0,0)){ return (-7);}
          }
          else
          {
               return (-9999);
          }

          if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::default_scalar_path)
          {
               std::int32_t ret_stat_scalar = this->generate_Q_channel_bitstream(Q_ch_samples_len);
               if(__builtin_expect(ret_stat_scalar<0,0)) { return (-8);}
          }
          else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_sse_path)
          {
               std::int32_t ret_stat_sse    = this->generate_Q_channel_bitstream_sse();
               if(__builtin_expect(ret_stat_sse<0,0))    { return (-9);}
          }
          else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_sse_u4x_path)
          {
               std::int32_t ret_stat_sse_u4x = this->generate_Q_channel_bitstream_sse_u4x();
               if(__builtin_expect(ret_stat_sse_u4x<0,0)) {return (-10);}
          }
          else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_avx_path)
          {
               std::int32_t ret_stat_avx     = this->generate_Q_channel_bitstream_avx(I_ch_do_const_prefetch);
               if(__builtin_expect(ret_stat_avx<0,0))     {return (-11);}
          }
          else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_avx_u4x_path)
          {
               std::int32_t ret_stat_avx_u4x = this->generate_Q_channel_bitstream_avx_u4x(I_ch_do_const_prefetch);
               if(__builtin_expect(ret_stat_avx_u4x<0,0)) { return (-12);}
          }
          else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_avx512_path)
          {
               std::int32_t ret_stat_avx512  = this->generate_Q_channel_bitstream_avx512(I_ch_do_const_prefetch);
               if(__builtin_expect(ret_stat_avx512<0,0))  { return (-13);}
          }
          else if constexpr (Q_ch_path == Q_channel_bitstream_optimized_path::vector_avx512_u4x_path)
          {
               std::int32_t ret_stat_avx512_u4x= this->generate_Q_channel_bitstream_avx512_u4x(I_ch_do_const_prefetch);
               if(__builtin_expect(ret_stat_avx512_u4x<0,0)){ return (-14);}
          }
          else 
          {
             return (-9999);
          }

     }

     if constexpr (psf_path == pulse_shaping_function_optimized_path::default_scalar_path)
     {
          std::int32_t ret_stat_scalar = this->generate_pulse_shaping_function_u8x();
          if(__builtin_expect(ret_stat_scalar<0,0)) { return (-15);}
     }
     else if constexpr (psf_path == pulse_shaping_function_optimized_path::vector_sse_path)
     {
          std::int32_t ret_stat_sse    = this->generate_pulse_shaping_function_sse_u4x(psfunc_do_const_prefetch);
          if(__builtin_expect(ret_stat_sse<0,0))    { return (-16);}
     }
     else if constexpr (psf_path == pulse_shaping_function_optimized_path::vector_avx_path)
     {
          std::int32_t ret_stat_avx    = this->generate_pulse_shaping_function_avx_u4x(psfunc_do_const_prefetch);
          if(__builtin_expect(ret_stat_avx<0,0))    { return (-17);}
     }
     else if constexpr (psf_path == pulse_shaping_function_optimized_path::vector_avx512_path)
     {
          std::int32_t ret_stat_avx512 = this->generate_pulse_shaping_function_avx512_u4x(psfunc_do_const_prefetch);
          if(__builtin_expect(ret_stat_avx512<0,0))  { return (-18);}
     }
     else 
     {
          return (-9999);
     }

     constexpr std::size_t LUT_loop_idx_threshold{2257ull};
     const __m512 vk{_mm512_set1_ps(this->m_k)};
     const __m512 vAc{_mm512_set1_ps(this->m_Ac)};
     const __m512 vph0{_mm512_set1_ps(this->m_ph0)};
     __m512 c0_15,c16_31;
     __m512 c32_47,c48_63;
     __m512 vpi,v2pi;
     __m512 vzero,vone;  
     std::size_t i;
     float j;
     if(!do_const_prefetch)
     {
         c0_15  = _mm512_setr_ps(0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f);
         c16_31 = _mm512_setr_ps(16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f,24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f);
         c32_47 = _mm512_setr_ps(32.0f,33.0f,34.0f,35.0f,36.0f,37.0f,38.0f,39.0f,40.0f,41.0f,42.0f,43.0f,44.0f,45.0f,46.0f,47.0f);
         c48_63 = _mm512_setr_ps(48.0f,49.0f,50.0f,51.0f,52.0f,53.0f,54.0f,55.0f,56.0f,57.0f,58.0f,59.0f,60.0f,61.0f,62.0f,63.0f);
         vpi    = _mm512_set1_ps(3.14159265358979323846264338328f);
         v2pi   = _mm512_set1_ps(6.283185307179586476925286766559f); 
         vzero  = _mm512_setzero_ps();
         vone   = _mm512_set1_ps(1.0f);
     }
     else 
     {
          __ATTR_ALIGN__(64) const float prefetched_constants[128] = {
                                0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f,
                                16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f,24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f,
                                32.0f,33.0f,34.0f,35.0f,36.0f,37.0f,38.0f,39.0f,40.0f,41.0f,42.0f,43.0f,44.0f,45.0f,46.0f,47.0f,
                                48.0f,49.0f,50.0f,51.0f,52.0f,53.0f,54.0f,55.0f,56.0f,57.0f,58.0f,59.0f,60.0f,61.0f,62.0f,63.0f,
                                3.14159265358979323846264338328f,3.14159265358979323846264338328f,3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f,3.14159265358979323846264338328f,3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f,3.14159265358979323846264338328f,3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f,3.14159265358979323846264338328f,3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f,3.14159265358979323846264338328f,3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f,6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
                                0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,
                                1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f};

          _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[32],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[48],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[64],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[80],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[96],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[112],_MM_HINT_T0);
          c0_15  = _mm512_load_ps(&prefetched_constants[0]);
          c16_31 = _mm512_load_ps(&prefetched_constants[16]); 
          c32_47 = _mm512_load_ps(&prefetched_constants[32]);
          c48_63 = _mm512_load_ps(&prefetched_constants[48]);
          vpi    = _mm512_load_ps(&prefetched_constants[64]);
          v2pi   = _mm512_load_ps(&prefetched_constants[80]);
          vzero  = _mm512_load_ps(&prefetched_constants[96]);
          vone   = _mm512_load_ps(&prefetched_constants[112]);
     }
     
     if(this->m_sfsk_nsamples>LUT_loop_idx_threshold)
     {
          for(i = 0ull,j = 0.0f; (i+63ull) < this->m_sfsk_nsamples; i += 64ull,j += 64.0f) 
          {
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m512 vt_i_0{_mm512_add_ps(_mm512_set1_ps(j),c0_15)};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m512 vik_0{_mm512_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m512 vqk_0{_mm512_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i+0ull],_MM_HINT_T0);
               const __m512 vdk_0{_mm512_mul_ps(_mm512_negate_ps(vik_0),vqk_0)};
               const __m512 vpsfunc_0{_mm512_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask16 vik_eq_1_0{_mm512_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m512 vPhik_0{_mm512_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m512 vcarrier_term1_0{_mm512_mul_ps(vk,_mm512_mul_ps(v2pi,vt_i_0))};
               const __m512 vcarrier_term2_0{_mm512_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m512 vsfsk_sample_0{_mm512_mul_ps(vAc,_mm512_sin_ps(_mm512_add_ps(vph0,_mm512_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm512_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0);
               //////////////////////////////////////////////////////////////////
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T0);
               const __m512 vt_i_1{_mm512_add_ps(_mm512_set1_ps(j),c16_31)};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T0);
               const __m512 vik_1{_mm512_load_ps(&this->m_I_ch_bitstream.m_data[i+16ull])};
               const __m512 vqk_1{_mm512_load_ps(&this->m_Q_ch_bitstream.m_data[i+16ull])};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i+16ull],_MM_HINT_T0);
               const __m512 vdk_1{_mm512_mul_ps(_mm512_negate_ps(vik_1),vqk_1)};
               const __m512 vpsfunc_1{_mm512_load_ps(&this->m_psfunc.m_data[i+16ull])};
               const __mmask16 vik_eq_1_1{_mm512_cmp_ps_mask(vik_1,vone,_CMP_EQ_OQ)};
               const __m512 vPhik_1{_mm512_mask_blend_ps(vik_eq_1_1,vzero,vpi)};
               const __m512 vcarrier_term1_1{_mm512_mul_ps(vk,_mm512_mul_ps(v2pi,vt_i_1))};
               const __m512 vcarrier_term2_1{_mm512_fmadd_ps(vdk_1,vpsfunc_1,vPhik_1)};
               const __m512 vsfsk_sample_1{_mm512_mul_ps(vAc,_mm512_sin_ps(_mm512_add_ps(vph0,_mm512_add_ps(vcarrier_term1_1,vcarrier_term2_1))))};
               _mm512_store_ps(&this->m_sfsk_signal.m_data[i+16ull],vsfsk_sample_1);
               /////////////////////////////////////////////////////////////////////
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+32ull],_MM_HINT_T0);
               const __m512 vt_i_2{_mm512_add_ps(_mm512_set1_ps(j),c32_47)};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+32ull],_MM_HINT_T0);
               const __m512 vik_2{_mm512_load_ps(&this->m_I_ch_bitstream.m_data[i+32ull])};
               const __m512 vqk_2{_mm512_load_ps(&this->m_Q_ch_bitstream.m_data[i+32ull])};
                _mm_prefetch((const char*)&this->m_psfunc.m_data[i+32ull],_MM_HINT_T0);
               const __m512 vdk_2{_mm512_mul_ps(_mm512_negate_ps(vik_2),vqk_2)};
               const __m512 vpsfunc_2{_mm512_load_ps(&this->m_psfunc.m_data[i+32ull])};
               const __mmask16 vik_eq_1_2{_mm512_cmp_ps_mask(vik_2,vone,_CMP_EQ_OQ)};
               const __m512 vPhik_2{_mm512_mask_blend_ps(vik_eq_1_2,vzero,vpi)};
               const __m512 vcarrier_term1_2{_mm512_mul_ps(vk,_mm512_mul_ps(v2pi,vt_i_2))};
               const __m512 vcarrier_term2_2{_mm512_fmadd_ps(vdk_2,vpsfunc_2,vPhik_2)};
               const __m512 vsfsk_sample_2{_mm512_mul_ps(vAc,_mm512_sin_ps(_mm512_add_ps(vph0,_mm512_add_ps(vcarrier_term1_2,vcarrier_term2_2))))};
               _mm512_store_ps(&this->m_sfsk_signal.m_data[i+32ull],vsfsk_sample_2);
               //////////////////////////////////////////////////////////////////////
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+48ull],_MM_HINT_T0);
               const __m512 vt_i_3{_mm512_add_ps(_mm512_set1_ps(j),c48_63)};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+48ull],_MM_HINT_T0);
               const __m512 vik_3{_mm512_load_ps(&this->m_I_ch_bitstream.m_data[i+48ull])};
               const __m512 vqk_3{_mm512_load_ps(&this->m_Q_ch_bitstream.m_data[i+48ull])};
               const __m512 vdk_3{_mm512_mul_ps(_mm512_negate_ps(vik_3),vqk_3)};
               const __m512 vpsfunc_3{_mm512_load_ps(&this->m_psfunc.m_data[i+48ull])};
               const __mmask16 vik_eq_1_3{_mm512_cmp_ps_mask(vik_3,vone,_CMP_EQ_OQ)};
               const __m512 vPhik_3{_mm512_mask_blend_ps(vik_eq_1_3,vzero,vpi)};
               const __m512 vcarrier_term1_3{_mm512_mul_ps(vk,_mm512_mul_ps(v2pi,vt_i_3))};
               const __m512 vcarrier_term2_3{_mm512_fmadd_ps(vdk_3,vpsfunc_3,vPhik_3)};
               const __m512 vsfsk_sample_3{_mm512_mul_ps(vAc,_mm512_sin_ps(_mm512_add_ps(vph0,_mm512_add_ps(vcarrier_term1_3,vcarrier_term2_3))))};
               _mm512_store_ps(&this->m_sfsk_signal.m_data[i+48ull],vsfsk_sample_3);
          }

          for(; (i+31ull) < this->m_sfsk_nsamples; i += 32ull,j += 32.0f) 
          {
                _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m512 vt_i_0{_mm512_add_ps(_mm512_set1_ps(j),c0_15)};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m512 vik_0{_mm512_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m512 vqk_0{_mm512_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i+0ull],_MM_HINT_T0);
               const __m512 vdk_0{_mm512_mul_ps(_mm512_negate_ps(vik_0),vqk_0)};
               const __m512 vpsfunc_0{_mm512_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask16 vik_eq_1_0{_mm512_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m512 vPhik_0{_mm512_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m512 vcarrier_term1_0{_mm512_mul_ps(vk,_mm512_mul_ps(v2pi,vt_i_0))};
               const __m512 vcarrier_term2_0{_mm512_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m512 vsfsk_sample_0{_mm512_mul_ps(vAc,_mm512_sin_ps(_mm512_add_ps(vph0,_mm512_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm512_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0);
               //////////////////////////////////////////////////////////////////
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T0);
               const __m512 vt_i_1{_mm512_add_ps(_mm512_set1_ps(j),c16_31)};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T0);
               const __m512 vik_1{_mm512_load_ps(&this->m_I_ch_bitstream.m_data[i+16ull])};
               const __m512 vqk_1{_mm512_load_ps(&this->m_Q_ch_bitstream.m_data[i+16ull])};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i+16ull],_MM_HINT_T0);
               const __m512 vdk_1{_mm512_mul_ps(_mm512_negate_ps(vik_1),vqk_1)};
               const __m512 vpsfunc_1{_mm512_load_ps(&this->m_psfunc.m_data[i+16ull])};
               const __mmask16 vik_eq_1_1{_mm512_cmp_ps_mask(vik_1,vone,_CMP_EQ_OQ)};
               const __m512 vPhik_1{_mm512_mask_blend_ps(vik_eq_1_1,vzero,vpi)};
               const __m512 vcarrier_term1_1{_mm512_mul_ps(vk,_mm512_mul_ps(v2pi,vt_i_1))};
               const __m512 vcarrier_term2_1{_mm512_fmadd_ps(vdk_1,vpsfunc_1,vPhik_1)};
               const __m512 vsfsk_sample_1{_mm512_mul_ps(vAc,_mm512_sin_ps(_mm512_add_ps(vph0,_mm512_add_ps(vcarrier_term1_1,vcarrier_term2_1))))};
               _mm512_store_ps(&this->m_sfsk_signal.m_data[i+16ull],vsfsk_sample_1);
              
          }

          for(; (i+15ull) < this->m_sfsk_nsamples; i += 16ull,j += 16.0f) 
          {
                
              _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m512 vt_i_0{_mm512_add_ps(_mm512_set1_ps(j),c0_15)};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m512 vik_0{_mm512_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m512 vqk_0{_mm512_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i+0ull],_MM_HINT_T0);
               const __m512 vdk_0{_mm512_mul_ps(_mm512_negate_ps(vik_0),vqk_0)};
               const __m512 vpsfunc_0{_mm512_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask16 vik_eq_1_0{_mm512_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m512 vPhik_0{_mm512_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m512 vcarrier_term1_0{_mm512_mul_ps(vk,_mm512_mul_ps(v2pi,vt_i_0))};
               const __m512 vcarrier_term2_0{_mm512_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m512 vsfsk_sample_0{_mm512_mul_ps(vAc,_mm512_sin_ps(_mm512_add_ps(vph0,_mm512_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm512_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0);
          }

          for(; (i+0ull) < this->m_sfsk_nsamples; i += 1ull,j += 1.0f) 
          {
               const float ik_1{this->m_I_ch_bitstream.m_data[i]};
               const float qk_1{this->m_Q_ch_bitstream.m_data[i]};
               const float dk_1{-ik_1*qk_1};
               const float psfunc_1{this->m_psfunc.m_data[i+1ull]};
               const float Phik_1{ik_1==1.0f?0.0f:3.14159265358979323846264338328};
               const float carrier_term1_1{this->m_k*6.283185307179586476925286766559*j};
               const float carrier_term2_1{std::fma(dk_1,psfunc_1,Phik_1)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_1{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#else 
               const float sfsk_sample_1{this->m_Ac*std::cos(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#endif 
               this->m_sfsk_signal.m_data[i] = sfsk_sample_1; 
          }
     }
     else 
     {
          for(i = 0ull; (i+63ull) < this->m_sfsk_nsamples; i += 64ull) 
          {
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+0ull],_MM_HINT_T0); 
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m512 vt_i_0{_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+0ull])};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m512 vik_0{_mm512_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m512 vqk_0{_mm512_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i+0ull],_MM_HINT_T0);
               const __m512 vdk_0{_mm512_mul_ps(_mm512_negate_ps(vik_0),vqk_0)};
               const __m512 vpsfunc_0{_mm512_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask16 vik_eq_1_0{_mm512_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m512 vPhik_0{_mm512_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m512 vcarrier_term1_0{_mm512_mul_ps(vk,_mm512_mul_ps(v2pi,vt_i_0))};
               const __m512 vcarrier_term2_0{_mm512_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m512 vsfsk_sample_0{_mm512_mul_ps(vAc,_mm512_sin_ps(_mm512_add_ps(vph0,_mm512_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm512_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0);
               //////////////////////////////////////////////////////////////////
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+16ull],_MM_HINT_T0);
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T0);
               const __m512 vt_i_1{_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+16ull])};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T0);
               const __m512 vik_1{_mm512_load_ps(&this->m_I_ch_bitstream.m_data[i+16ull])};
               const __m512 vqk_1{_mm512_load_ps(&this->m_Q_ch_bitstream.m_data[i+16ull])};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i+16ull],_MM_HINT_T0);
               const __m512 vdk_1{_mm512_mul_ps(_mm512_negate_ps(vik_1),vqk_1)};
               const __m512 vpsfunc_1{_mm512_load_ps(&this->m_psfunc.m_data[i+16ull])};
               const __mmask16 vik_eq_1_1{_mm512_cmp_ps_mask(vik_1,vone,_CMP_EQ_OQ)};
               const __m512 vPhik_1{_mm512_mask_blend_ps(vik_eq_1_1,vzero,vpi)};
               const __m512 vcarrier_term1_1{_mm512_mul_ps(vk,_mm512_mul_ps(v2pi,vt_i_1))};
               const __m512 vcarrier_term2_1{_mm512_fmadd_ps(vdk_1,vpsfunc_1,vPhik_1)};
               const __m512 vsfsk_sample_1{_mm512_mul_ps(vAc,_mm512_sin_ps(_mm512_add_ps(vph0,_mm512_add_ps(vcarrier_term1_1,vcarrier_term2_1))))};
               _mm512_store_ps(&this->m_sfsk_signal.m_data[i+16ull],vsfsk_sample_1);
               /////////////////////////////////////////////////////////////////////
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+32ull],_MM_HINT_T0);
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+32ull],_MM_HINT_T0);
               const __m512 vt_i_2{_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+32ull])};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+32ull],_MM_HINT_T0);
               const __m512 vik_2{_mm512_load_ps(&this->m_I_ch_bitstream.m_data[i+32ull])};
               const __m512 vqk_2{_mm512_load_ps(&this->m_Q_ch_bitstream.m_data[i+32ull])};
                _mm_prefetch((const char*)&this->m_psfunc.m_data[i+32ull],_MM_HINT_T0);
               const __m512 vdk_2{_mm512_mul_ps(_mm512_negate_ps(vik_2),vqk_2)};
               const __m512 vpsfunc_2{_mm512_load_ps(&this->m_psfunc.m_data[i+32ull])};
               const __mmask16 vik_eq_1_2{_mm512_cmp_ps_mask(vik_2,vone,_CMP_EQ_OQ)};
               const __m512 vPhik_2{_mm512_mask_blend_ps(vik_eq_1_2,vzero,vpi)};
               const __m512 vcarrier_term1_2{_mm512_mul_ps(vk,_mm512_mul_ps(v2pi,vt_i_2))};
               const __m512 vcarrier_term2_2{_mm512_fmadd_ps(vdk_2,vpsfunc_2,vPhik_2)};
               const __m512 vsfsk_sample_2{_mm512_mul_ps(vAc,_mm512_sin_ps(_mm512_add_ps(vph0,_mm512_add_ps(vcarrier_term1_2,vcarrier_term2_2))))};
               _mm512_store_ps(&this->m_sfsk_signal.m_data[i+32ull],vsfsk_sample_2);
               //////////////////////////////////////////////////////////////////////
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+48ull],_MM_HINT_T0);
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+48ull],_MM_HINT_T0);
               const __m512 vt_i_3{_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+48ull])};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+48ull],_MM_HINT_T0);
               const __m512 vik_3{_mm512_load_ps(&this->m_I_ch_bitstream.m_data[i+48ull])};
               const __m512 vqk_3{_mm512_load_ps(&this->m_Q_ch_bitstream.m_data[i+48ull])};
               const __m512 vdk_3{_mm512_mul_ps(_mm512_negate_ps(vik_3),vqk_3)};
               const __m512 vpsfunc_3{_mm512_load_ps(&this->m_psfunc.m_data[i+48ull])};
               const __mmask16 vik_eq_1_3{_mm512_cmp_ps_mask(vik_3,vone,_CMP_EQ_OQ)};
               const __m512 vPhik_3{_mm512_mask_blend_ps(vik_eq_1_3,vzero,vpi)};
               const __m512 vcarrier_term1_3{_mm512_mul_ps(vk,_mm512_mul_ps(v2pi,vt_i_3))};
               const __m512 vcarrier_term2_3{_mm512_fmadd_ps(vdk_3,vpsfunc_3,vPhik_3)};
               const __m512 vsfsk_sample_3{_mm512_mul_ps(vAc,_mm512_sin_ps(_mm512_add_ps(vph0,_mm512_add_ps(vcarrier_term1_3,vcarrier_term2_3))))};
               _mm512_store_ps(&this->m_sfsk_signal.m_data[i+48ull],vsfsk_sample_3);
          }

          for(; (i+31ull) < this->m_sfsk_nsamples; i += 32ull) 
          {
              
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+0ull],_MM_HINT_T0); 
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m512 vt_i_0{_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+0ull])};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m512 vik_0{_mm512_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m512 vqk_0{_mm512_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i+0ull],_MM_HINT_T0);
               const __m512 vdk_0{_mm512_mul_ps(_mm512_negate_ps(vik_0),vqk_0)};
               const __m512 vpsfunc_0{_mm512_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask16 vik_eq_1_0{_mm512_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m512 vPhik_0{_mm512_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m512 vcarrier_term1_0{_mm512_mul_ps(vk,_mm512_mul_ps(v2pi,vt_i_0))};
               const __m512 vcarrier_term2_0{_mm512_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m512 vsfsk_sample_0{_mm512_mul_ps(vAc,_mm512_sin_ps(_mm512_add_ps(vph0,_mm512_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm512_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0);
               //////////////////////////////////////////////////////////////////
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+16ull],_MM_HINT_T0);
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T0);
               const __m512 vt_i_1{_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+16ull])};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T0);
               const __m512 vik_1{_mm512_load_ps(&this->m_I_ch_bitstream.m_data[i+16ull])};
               const __m512 vqk_1{_mm512_load_ps(&this->m_Q_ch_bitstream.m_data[i+16ull])};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i+16ull],_MM_HINT_T0);
               const __m512 vdk_1{_mm512_mul_ps(_mm512_negate_ps(vik_1),vqk_1)};
               const __m512 vpsfunc_1{_mm512_load_ps(&this->m_psfunc.m_data[i+16ull])};
               const __mmask16 vik_eq_1_1{_mm512_cmp_ps_mask(vik_1,vone,_CMP_EQ_OQ)};
               const __m512 vPhik_1{_mm512_mask_blend_ps(vik_eq_1_1,vzero,vpi)};
               const __m512 vcarrier_term1_1{_mm512_mul_ps(vk,_mm512_mul_ps(v2pi,vt_i_1))};
               const __m512 vcarrier_term2_1{_mm512_fmadd_ps(vdk_1,vpsfunc_1,vPhik_1)};
               const __m512 vsfsk_sample_1{_mm512_mul_ps(vAc,_mm512_sin_ps(_mm512_add_ps(vph0,_mm512_add_ps(vcarrier_term1_1,vcarrier_term2_1))))};
               _mm512_store_ps(&this->m_sfsk_signal.m_data[i+16ull],vsfsk_sample_1);
          }

          for(; (i+15ull) < this->m_sfsk_nsamples; i += 16ull) 
          {
                
              _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+0ull],_MM_HINT_T0); 
               _mm_prefetch((const char*)&this->m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m512 vt_i_0{_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+0ull])};
               _mm_prefetch((const char*)&this->m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T0);
               const __m512 vik_0{_mm512_load_ps(&this->m_I_ch_bitstream.m_data[i+0ull])};
               const __m512 vqk_0{_mm512_load_ps(&this->m_Q_ch_bitstream.m_data[i+0ull])};
               _mm_prefetch((const char*)&this->m_psfunc.m_data[i+0ull],_MM_HINT_T0);
               const __m512 vdk_0{_mm512_mul_ps(_mm512_negate_ps(vik_0),vqk_0)};
               const __m512 vpsfunc_0{_mm512_load_ps(&this->m_psfunc.m_data[i+0ull])};
               const __mmask16 vik_eq_1_0{_mm512_cmp_ps_mask(vik_0,vone,_CMP_EQ_OQ)};
               const __m512 vPhik_0{_mm512_mask_blend_ps(vik_eq_1_0,vzero,vpi)};
               const __m512 vcarrier_term1_0{_mm512_mul_ps(vk,_mm512_mul_ps(v2pi,vt_i_0))};
               const __m512 vcarrier_term2_0{_mm512_fmadd_ps(vdk_0,vpsfunc_0,vPhik_0)};
               const __m512 vsfsk_sample_0{_mm512_mul_ps(vAc,_mm512_sin_ps(_mm512_add_ps(vph0,_mm512_add_ps(vcarrier_term1_0,vcarrier_term2_0))))};
               _mm512_store_ps(&this->m_sfsk_signal.m_data[i+0ull],vsfsk_sample_0);
          }

          for(; (i+0ull) < this->m_sfsk_nsamples; i += 1ull) 
          {
               const float t_i{gms::math::LUT_loop_indices_2257_align64[i]};
               const float ik_1{this->m_I_ch_bitstream.m_data[i]};
               const float qk_1{this->m_Q_ch_bitstream.m_data[i]};
               const float dk_1{-ik_1*qk_1};
               const float psfunc_1{this->m_psfunc.m_data[i+1ull]};
               const float Phik_1{ik_1==1.0f?0.0f:3.14159265358979323846264338328};
               const float carrier_term1_1{this->m_k*6.283185307179586476925286766559*t_i};
               const float carrier_term2_1{std::fma(dk_1,psfunc_1,Phik_1)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
               const float sfsk_sample_1{this->m_Ac*ceph_cosf(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#else 
               const float sfsk_sample_1{this->m_Ac*std::cos(this->m_ph0+carrier_term1_1+carrier_term2_1)};
#endif 
               this->m_sfsk_signal.m_data[i] = sfsk_sample_1; 
          }
     }

     return (0);

}

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_avx512_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::default_scalar_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::default_scalar_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::default_scalar_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_avx512_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_sse_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_sse_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_sse_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_avx512_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_sse_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_sse_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_sse_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_avx512_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_avx512_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_avx512_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx512_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx512_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx512_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

template std::int32_t 
gms::radiolocation
::sinusoidal_fsk_t
::generate_fsk_signal_avx512_u4x<gms::radiolocation::sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx512_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx512_u4x_path,
                                 gms::radiolocation::sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx512_path>(const iq_rectw_bitstream_vsequence_t &,
                                                                                                              const std::int32_t,
                                                                                                              const std::int32_t,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool,
                                                                                                              const bool);

auto 
gms::radiolocation 
::operator<<(std::ostream &os, 
             const sinusoidal_fsk_t &rhs)->std::ostream & 
{
    std::cout << typeid(rhs).name() << "Begin: object state dump." << std::endl;
    std::cout << "m_sfsk_nsamples    : "      << rhs.m_sfsk_nsamples  << std::endl;
    std::cout << "m_I_ch_nsamples    : "      << rhs.m_I_ch_nsamples  << std::endl;
    std::cout << "m_Q_ch_nsamples    : "      << rhs.m_Q_ch_nsamples  << std::endl;
    std::cout << "m_psfunc_nsamples  : "      << rhs.m_psfunc_nsamples << std::endl;
    std::cout << "m_k                : "      << std::fixed << std::setprecision(7) << rhs.m_k << std::endl;
    std::cout << "m_T                : "      << std::fixed << std::setprecision(7) << rhs.m_T << std::endl;
    std::cout << "m_Ac               : "      << std::fixed << std::setprecision(7) << rhs.m_Ac << std::endl;
    std::cout << "m_fc               : "      << std::fixed << std::setprecision(7) << rhs.m_fc << std::endl;
    std::cout << "m_ph0              : "      << std::fixed << std::setprecision(7) << rhs.m_ph0 << std::endl;
    std::cout << "m_I_fc             : "      << std::fixed << std::setprecision(7) << rhs.m_I_fc << std::endl;
    std::cout << "m_Q_fc             : "      << std::fixed << std::setprecision(7) << rhs.m_Q_fc << std::endl;
    std::cout << "m_I_ph0            : "      << std::fixed << std::setprecision(7) << rhs.m_I_ph0 << std::endl;
    std::cout << "m_Q_ph0            : "      << std::fixed << std::setprecision(7) << rhs.m_Q_ph0 << std::endl;
   
    std::cout << "Sinusoidal FSK -- Samples:" << "\n\n";
    for(std::size_t i{0ull}; i != rhs.m_sfsk_nsamples; ++i)
    {
        os << std::fixed << std::setprecision(7) << rhs.m_sfsk_signal.m_data[i] << std::endl;
    
    }
    std::cout << "I-Bitstream -- Samples:" << "\n\n";
    for(std::size_t i{0ull}; i != rhs.m_I_ch_nsamples; ++i) 
    {
       os << std::fixed << std::setprecision(7) << rhs.m_I_ch_bitstream.m_data[i] << std::endl;
    }
    std::cout << "Q-Bitstream -- Samples:" << "\n\n";
    for(std::size_t i{0ull}; i != rhs.m_Q_ch_nsamples; ++i) 
    {
       os << std::fixed << std::setprecision(7) << rhs.m_Q_ch_bitstream.m_data[i] << std::endl;
    }
    std::cout << "Pulse Shaping Function -- Samples:" << "\n\n";
    for(std::size_t i{0ull}; i != rhs.m_psfunc_nsamples; ++i) 
    {
       os << std::fixed << std::setprecision(7) << rhs.m_psfunc.m_data[i] << std::endl;
    }
    std::cout << typeid(rhs).name() << "End: object state dump." << std::endl;
    return (os);
}





