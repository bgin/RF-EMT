
#include <fstream>
#include <iomanip>
#include <random>
#include "GMS_sinusoidal_weighted_oqpsk.h"
#include "GMS_sse_memset.h"


gms::radiolocation 
::sinusoidal_weighted_oqpsk_t
::sinusoidal_weighted_oqpsk_t(const std::size_t nsamples,
                              const std::size_t I_ch_nsamples,
                              const std::size_t Q_ch_nsamples,
                              const float A_I,
                              const float A_Q,
                              const float T,
                              const float cw0,
                              const float sw0,
                              const float cph0,
                              const float sph0)
:
m_nsamples{nsamples},
m_I_ch_nsamples{I_ch_nsamples},
m_Q_ch_nsamples{Q_ch_nsamples},
m_A_I{A_I},
m_A_Q{A_Q},
m_T{T},
m_cw0{cw0},
m_sw0{sw0},
m_cph0{cph0},
m_sph0{sph0},
m_have_msk_samples{false},
m_I_bitstream{darray_r4_t(this->m_I_ch_nsamples)},
m_Q_bitstream{darray_r4_t(this->m_Q_ch_nsamples)},
m_I_channel{darray_r4_t(this->m_nsamples)},
m_Q_channel{darray_r4_t(this->m_nsamples)},
m_msk_signal{darray_r4_t(this->m_nsamples)} 
{

}

gms::radiolocation 
::sinusoidal_weighted_oqpsk_t
::sinusoidal_weighted_oqpsk_t(const sinusoidal_weighted_oqpsk_t &other) noexcept(false) 
: 
m_nsamples{other.m_nsamplex},
m_I_ch_nsamples{other.m_I_ch_nsamples},
m_Q_ch_nsamples{other.m_Q_ch_nsamples},
m_A_I{other.m_A_I},
m_A_Q{other.m_A_Q},
m_T{other.m_T},
m_cw0{other.m_cw0},
m_sw0{other.m_sw0},
m_cph0{other.m_cph0},
m_sph0{other.m_sph0},
m_have_msk_samples{other.m_have_msk_samples},
m_I_bitstream{other.m_I_bitstream},
m_Q_bitstream{other.m_Q_bitstream},
m_I_channel{other.m_I_channel},
m_Q_channel{other.m_Q_channel},
m_msk_signal{other.m_msk_signal} 
{

}

gms::radiolocation 
::sinusoidal_weighted_oqpsk_t
::sinusoidal_weighted_oqpsk_t(sinusoidal_weighted_oqpsk_t &&other) noexcept(true) 
: 
m_nsamples{std::move(other.m_nsamplex)},
m_I_ch_nsamples{std::move(other.m_I_ch_nsamples)},
m_Q_ch_nsamples{std::move(other.m_Q_ch_nsamples)},
m_A_I{std::move(other.m_A_I)},
m_A_Q{std::move(other.m_A_Q)},
m_T{std::move(other.m_T)},
m_cw0{std::move(other.m_cw0)},
m_sw0{std::move(other.m_sw0)},
m_cph0{std::move(other.m_cph0)},
m_sph0{std::move(other.m_sph0)},
m_have_msk_samples{std::move(other.m_have_msk_samples)},
m_I_bitstream{std::move(other.m_I_bitstream)},
m_Q_bitstream{std::move(other.m_Q_bitstream)},
m_I_channel{std::move(other.m_I_channel)},
m_Q_channel{std::move(other.m_Q_channel)},
m_msk_signal{std::move(other.m_msk_signal)} 
{

}

gms::radiolocation 
::sinusoidal_weighted_oqpsk_t
::~sinusoidal_weighted_oqpsk_t() noexcept(true) 
{

}

gms::radiolocation 
::sinusoidal_weighted_oqpsk_t &
gms::radiolocation 
::sinusoidal_weighted_oqpsk_t
::operator=(const sinusoidal_weighted_oqpsk_t &other) noexcept(false) 
{
      if(__builtin_expect(this==&other,0)) {return (*this)};
      this->m_nsamples         = other.m_nsamples;
      this->m_I_ch_nsamples    = other.m_I_ch_nsamples;
      this->m_Q_ch_nsamples    = other.m_Q_ch_nsamples;
      this->m_A_I              = other.m_A_I;
      this->m_A_Q              = other.m_A_Q;
      this->m_T                = other.m_T;
      this->m_cw0              = other.m_cw0;
      this->m_sw0              = other.m_sw0;
      this->m_cph0             = other.m_cph0;
      this->m_sph0             = other.m_sph0;
      this->m_have_msk_samples = other.m_have_msk_samples;
      this->m_I_bitstream.operator=(other.m_I_bitstream);
      this->m_Q_bitstream.operator=(other.m_Q_bitstream);
      this->m_I_channel.operator=(other.m_I_channel);
      this->m_Q_channel.operator=(other.m_Q_channel);
      this->m_msk_signal.operator=(other.m_msk_signal);
      return (*this);
}

gms::radiolocation 
::sinusoidal_weighted_oqpsk_t &
gms::radiolocation 
::sinusoidal_weighted_oqpsk_t
::operator=(sinusoidal_weighted_oqpsk_t &&other) noexcept(true) 
{
      if(__builtin_expect(this==&other,0)) {return (*this)};
      this->m_nsamples      = std::move(other.m_nsamples);
      this->m_I_ch_nsamples = std::move(other.m_I_ch_nsamples);
      this->m_Q_ch_nsamples = std::move(other.m_Q_ch_nsamples);
      this->m_A_I           = std::move(other.m_A_I);
      this->m_A_Q           = std::move(other.m_A_Q);
      this->m_T             = std::move(other.m_T);
      this->m_cw0           = std::move(other.m_cw0);
      this->m_sw0           = std::move(other.m_sw0);
      this->m_cph0          = std::move(other.m_cph0);
      this->m_sph0          = std::move(other.m_sph0);
      this->m_have_msk_samples = std::move(other.m_have_msk_samples);
      this->m_I_bitstream.operator=(std::move(other.m_I_bitstream));
      this->m_Q_bitstream.operator=(std::move(other.m_Q_bitstream));
      this->m_I_channel.operator=(std::move(other.m_I_channel));
      this->m_Q_channel.operator=(std::move(other.m_Q_channel));
      this->m_msk_signal.operator=(std::move(other.m_msk_signal));
      return (*this);
}

void 
gms::radiolocation 
::sinusoidal_weighted_oqpsk_t
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
::sinusoidal_weighted_oqpsk_t
::generate_oqpsk_signal( const float I_duration, // user passed
                         const float I_w0,       // user passed
                         const float I_ph0,      // user passed
                         const float I_sample_rate,
                         const float Q_duration, // user passed
                         const float Q_w0,       // user passed
                         const float Q_ph0,      // user passed
                         const float Q_sample_rate)
{
      std::int32_t I_ch_bitstream_stat{1};
      std::int32_t Q_ch_bitstream_stat{1};
      I_ch_bitstream_stat = this->generate_I_channel_bitstream(I_duration,I_w0,I_ph0,I_sample_rate);
      if(I_ch_bitstream_stat < 0) { return (-1);}
      Q_ch_bitstream_stat = this->generate_Q_channel_bitstream(Q_duration,Q_w0,Q_ph0,Q_sample_rate); 
      if(Q_ch_bitstream_stat < 0) {return (-2);}
      
      constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
      constexpr float C314159265358979323846264338328{3.14159265358979323846264338328f};
      std::size_t i,j;
      float cos_carrier0;
      float cos_carrier1;
      float cos_carrier2;
      float cos_carrier3;
      float sin_carrier0;
      float sin_carrier1;
      float sin_carrier2;
      float sin_carrier3;
      float cos_weight0;
      float cos_weight1;
      float cos_weight2;
      float cos_weight3;
      float sin_weight0;
      float sin_weight1;
      float sin_weight2;
      float sin_weight3;
      const float invTc{1.0f/this->m_T};
      const float inv2T{1.0f/(this->m_T+this->m_T)};

      for(i = 0ull; i != ROUND_TO_FOUR(this->m_nsamples,4ull); i += 4ull) 
      {
            const float t_i_0{static_cast<std::size_t>(i)};
            const float I_t0{this->m_I_bitstream.m_data[i+0ull]};
            const float I_data0{this->m_A_I*I_t0};
#if (SINUSOIDAL_WEIGHTED_OQPSK_USE_CEPHES) == 1 
            cos_carrier0 = ceph_cosf(this->m_cph0+(C6283185307179586476925286766559*this->m_cw0*t_i_0*invTc));
            cos_weight0  = ceph_cosf(C314159265358979323846264338328*t_i_0*inv2T);
#else 
            cos_carrier0 = std::cos(this->m_cph0+(C6283185307179586476925286766559*this->m_cw0*t_i_0*invTc));
            cos_weight0  = std::cos(C314159265358979323846264338328*t_i_0*inv2T);
#endif 
            const float cos_term0 = I_data0*cos_carrier0*cos_weight0;
            this->m_I_channel.m_data[i+0ull] = cos_term0;

            const float Q_t0{this->m_Q_bitstream.m_data[i+0ull]};
            const float Q_data0{this->m_A_Q*Q_t0};
#if (SINUSOIDAL_WEIGHTED_OQPSK_USE_CEPHES) == 1 
            sin_carrier0 = ceph_sinf(this->m_sph0+(C6283185307179586476925286766559*this->m_sw0*t_i_0*invTc));
            sin_weight0  = ceph_sinf(C314159265358979323846264338328*t_i_0*inv2T);
#else 
            sin_carrier0 = std::sin(this->m_sph0+(C6283185307179586476925286766559*this->m_sw0*t_i_0*invTc));
            sin_weight0  = std::sin(C314159265358979323846264338328*t_i_0*inv2T);
#endif 
            const float sin_term0 = Q_data0*sin_carrier0*sin_weight0;
            this->m_Q_channel.m_data[i+0ull] = sin_term0;
            this->m_msk_signal.m_data[i+0ull]= cos_term0+sin_term0;
            /////////////////////////////////////////////////////////////////////////////////////////////////
            const float t_i_1{static_cast<std::size_t>(i+1ull)};
            const float I_t1{this->m_I_bitstream.m_data[i+1ull]};
            const float I_data1{this->m_A_I*I_t1};
#if (SINUSOIDAL_WEIGHTED_OQPSK_USE_CEPHES) == 1 
            cos_carrier1 = ceph_cosf(this->m_cph0+(C6283185307179586476925286766559*this->m_cw0*t_i_1*invTc));
            cos_weight1  = ceph_cosf(C314159265358979323846264338328*t_i_1*inv2T);
#else 
            cos_carrier1 = std::cos(this->m_cph0+(C6283185307179586476925286766559*this->m_cw0*t_i_1*invTc));
            cos_weight1  = std::cos(C314159265358979323846264338328*t_i_1*inv2T);
#endif 
            const float cos_term1 = I_data1*cos_carrier1*cos_weight1;
            this->m_I_channel.m_data[i+1ull] = cos_term1;

            const float Q_t1{this->m_Q_bitstream.m_data[i+1ull]};
            const float Q_data1{this->m_A_Q*Q_t1};
#if (SINUSOIDAL_WEIGHTED_OQPSK_USE_CEPHES) == 1 
            sin_carrier1 = ceph_sinf(this->m_sph0+(C6283185307179586476925286766559*this->m_sw0*t_i_1*invTc));
            sin_weight1  = ceph_sinf(C314159265358979323846264338328*t_i_1*inv2T);
#else 
            sin_carrier1 = std::sin(this->m_sph0+(C6283185307179586476925286766559*this->m_sw0*t_i_1*invTc));
            sin_weight1  = std::sin(C314159265358979323846264338328*t_i_1*inv2T);
#endif 
            const float sin_term1 = Q_data1*sin_carrier1*sin_weight1;
            this->m_Q_channel.m_data[i+1ull] = sin_term1;
            this->m_msk_signal.m_data[i+1ull]= cos_term1+sin_term1;
            ////////////////////////////////////////////////////////////////////////////////////////////////
            const float t_i_2{static_cast<std::size_t>(i+2ull)};
            const float I_t2{this->m_I_bitstream.m_data[i+2ull]};
            const float I_data2{this->m_A_I*I_t2};
#if (SINUSOIDAL_WEIGHTED_OQPSK_USE_CEPHES) == 1 
            cos_carrier2 = ceph_cosf(this->m_cph0+(C6283185307179586476925286766559*this->m_cw0*t_i_2*invTc));
            cos_weight2  = ceph_cosf(C314159265358979323846264338328*t_i_2*inv2T);
#else 
            cos_carrier2 = std::cos(this->m_cph0+(C6283185307179586476925286766559*this->m_cw0*t_i_2*invTc));
            cos_weight2  = std::cos(C314159265358979323846264338328*t_i_2*inv2T);
#endif 
            const float cos_term2 = I_data2*cos_carrier2*cos_weight2;
            this->m_I_channel.m_data[i+2ull] = cos_term2;

            const float Q_t2{this->m_Q_bitstream.m_data[i+2ull]};
            const float Q_data2{this->m_A_Q*Q_t2};
#if (SINUSOIDAL_WEIGHTED_OQPSK_USE_CEPHES) == 1 
            sin_carrier2 = ceph_sinf(this->m_sph0+(C6283185307179586476925286766559*this->m_sw0*t_i_2*invTc));
            sin_weight2  = ceph_sinf(C314159265358979323846264338328*t_i_2*inv2T);
#else 
            sin_carrier2 = std::sin(this->m_sph0+(C6283185307179586476925286766559*this->m_sw0*t_i_2*invTc));
            sin_weight2  = std::sin(C314159265358979323846264338328*t_i_2*inv2T);
#endif 
            const float sin_term2 = Q_data2*sin_carrier2*sin_weight2;
            this->m_Q_channel.m_data[i+2ull] = sin_term2;
            this->m_msk_signal.m_data[i+2ull]= cos_term2+sin_term2;
            /////////////////////////////////////////////////////////////////////////////////////////////////
            const float t_i_3{static_cast<std::size_t>(i+3ull)};
            const float I_t3{this->m_I_bitstream.m_data[i+3ull]};
            const float I_data3{this->m_A_I*I_t3};
#if (SINUSOIDAL_WEIGHTED_OQPSK_USE_CEPHES) == 1 
            cos_carrier3 = ceph_cosf(this->m_cph0+(C6283185307179586476925286766559*this->m_cw0*t_i_3*invTc));
            cos_weight3  = ceph_cosf(C314159265358979323846264338328*t_i_3*inv2T);
#else 
            cos_carrier3 = std::cos(this->m_cph0+(C6283185307179586476925286766559*this->m_cw0*t_i_3*invTc));
            cos_weight3  = std::cos(C314159265358979323846264338328*t_i_3*inv2T);
#endif 
            const float cos_term3 = I_data3*cos_carrier3*cos_weight3;
            this->m_I_channel.m_data[i+3ull] = cos_term3;

            const float Q_t3{this->m_Q_bitstream.m_data[i+3ull]};
            const float Q_data3{this->m_A_Q*Q_t3};
#if (SINUSOIDAL_WEIGHTED_OQPSK_USE_CEPHES) == 1 
            sin_carrier3 = ceph_sinf(this->m_sph0+(C6283185307179586476925286766559*this->m_sw0*t_i_3*invTc));
            sin_weight3  = ceph_sinf(C314159265358979323846264338328*t_i_3*inv2T);
#else 
            sin_carrier3 = std::sin(this->m_sph0+(C6283185307179586476925286766559*this->m_sw0*t_i_3*invTc));
            sin_weight3  = std::sin(C314159265358979323846264338328*t_i_3*inv2T);
#endif 
            const float sin_term3 = Q_data3*sin_carrier3*sin_weight3;
            this->m_Q_channel.m_data[i+3ull] = sin_term3;
            this->m_msk_signal.m_data[i+3ull]= cos_term3+sin_term3;
      }

      for(j = i; j != this->m_nsamples; ++j)
      {
            const float t_j{static_cast<std::size_t>(j)};
            const float I_t0{this->m_I_bitstream.m_data[j]};
            const float I_data0{this->m_A_I*I_t0};
#if (SINUSOIDAL_WEIGHTED_OQPSK_USE_CEPHES) == 1 
            cos_carrier0 = ceph_cosf(this->m_cph0+(C6283185307179586476925286766559*this->m_cw0*t_j*invTc));
            cos_weight0  = ceph_cosf(C314159265358979323846264338328*t_j*inv2T);
#else 
            cos_carrier0 = std::cos(this->m_cph0+(C6283185307179586476925286766559*this->m_cw0*t_j*invTc));
            cos_weight0  = std::cos(C314159265358979323846264338328*t_j*inv2T);
#endif 
            const float cos_term0 = I_data0*cos_carrier0*cos_weight0;
            this->m_I_channel.m_data[j] = cos_term0;

            const float Q_t0{this->m_Q_bitstream.m_data[j]};
            const float Q_data0{this->m_A_Q*Q_t0};
#if (SINUSOIDAL_WEIGHTED_OQPSK_USE_CEPHES) == 1 
            sin_carrier0 = ceph_sinf(this->m_sph0+(C6283185307179586476925286766559*this->m_sw0*t_j*invTc));
            sin_weight0  = ceph_sinf(C314159265358979323846264338328*t_j*inv2T);
#else 
            sin_carrier0 = std::sin(this->m_sph0+(C6283185307179586476925286766559*this->m_sw0*t_j*invTc));
            sin_weight0  = std::sin(C314159265358979323846264338328*t_j*inv2T);
#endif 
            const float sin_term0 = Q_data0*sin_carrier0*sin_weight0;
            this->m_Q_channel.m_data[j] = sin_term0;
            this->m_msk_signal.m_data[j]= cos_term0+sin_term0;
      }
      // sanity check.
      if(__builtin_expect(this->m_have_msk_samples==false,0)) { this->m_have_msk_samples = true;}
      return (0);
}


                              