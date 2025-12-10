
#include <fstream>
#include <iomanip>
#include "GMS_sinusoidal_weighted_oqpsk.h"



gms::radiolocation 
::sinusoidal_weighted_oqpsk_t
::sinusoidal_weighted_oqpsk_t(const std::size_t nsamples,
                              const std::size_t I_ch_nsamples,
                              const std::size_t Q_ch_nsamples,
                              const std::size_t nfrequencies,
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
m_nfrequencies{nfrequencies},
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
m_Psi_f{darray_r4_t(this->m_nfrequencies)},
m_msk_signal{darray_r4_t(this->m_nsamples)} 
{

}

gms::radiolocation 
::sinusoidal_weighted_oqpsk_t
::sinusoidal_weighted_oqpsk_t(const sinusoidal_weighted_oqpsk_t &other) noexcept(false) 
: 
m_nsamples{other.m_nsamples},
m_I_ch_nsamples{other.m_I_ch_nsamples},
m_Q_ch_nsamples{other.m_Q_ch_nsamples},
m_nfrequencies{other.m_nfrequencies},
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
m_Psi_f{other.m_Psi_f},
m_msk_signal{other.m_msk_signal} 
{

}

gms::radiolocation 
::sinusoidal_weighted_oqpsk_t
::sinusoidal_weighted_oqpsk_t(sinusoidal_weighted_oqpsk_t &&other) noexcept(true) 
: 
m_nsamples{std::move(other.m_nsamples)},
m_I_ch_nsamples{std::move(other.m_I_ch_nsamples)},
m_Q_ch_nsamples{std::move(other.m_Q_ch_nsamples)},
m_nfrequencies{std::move(other.m_nfrequencies)},
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
m_Psi_f{std::move(other.m_Psi_f)},
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
      if(__builtin_expect(this==&other,0)) {return (*this);}
      this->m_nsamples         = other.m_nsamples;
      this->m_I_ch_nsamples    = other.m_I_ch_nsamples;
      this->m_Q_ch_nsamples    = other.m_Q_ch_nsamples;
      this->m_nfrequencies     = other.m_nfrequencies;
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
      this->m_Psi_f.operator=(other.m_Psi_f);
      this->m_msk_signal.operator=(other.m_msk_signal);
      return (*this);
}

gms::radiolocation 
::sinusoidal_weighted_oqpsk_t &
gms::radiolocation 
::sinusoidal_weighted_oqpsk_t
::operator=(sinusoidal_weighted_oqpsk_t &&other) noexcept(true) 
{
      if(__builtin_expect(this==&other,0)) {return (*this);}
      this->m_nsamples      = std::move(other.m_nsamples);
      this->m_I_ch_nsamples = std::move(other.m_I_ch_nsamples);
      this->m_Q_ch_nsamples = std::move(other.m_Q_ch_nsamples);
      this->m_nfrequencies  = std::move(other.m_nfrequencies);
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
      this->m_Psi_f.operator=(std::move(other.m_Psi_f));
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
      using namespace gms::math;
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

std::int32_t 
gms::radiolocation 
::sinusoidal_weighted_oqpsk_t
::generate_oqpsk_signal_additive_noise(const float I_duration, // user passed
                                       const float I_w0,       // user passed
                                       const float I_ph0,      // user passed
                                       const float I_sample_rate,
                                       const float Q_duration, // user passed
                                       const float Q_w0,       // user passed
                                       const float Q_ph0,      // user passed
                                       const float Q_sample_rate,
                                       sinusoidal_weighted_oqpsk_pdf_params_t & pdf_params,
                                       sinusoidal_weighted_oqpsk_rand_distr rd_enum)
{
      using namespace gms::math;
      std::int32_t I_ch_bitstream_stat{1};
      std::int32_t Q_ch_bitstream_stat{1};
      I_ch_bitstream_stat = this->generate_I_channel_bitstream(I_duration,I_w0,I_ph0,I_sample_rate);
      if(I_ch_bitstream_stat < 0) { return (-1);}
      Q_ch_bitstream_stat = this->generate_Q_channel_bitstream(Q_duration,Q_w0,Q_ph0,Q_sample_rate); 
      if(Q_ch_bitstream_stat < 0) {return (-2);}    
      constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
      constexpr float C314159265358979323846264338328{3.14159265358979323846264338328f};
      const float invTc{1.0f/this->m_T};
      const float inv2T{1.0f/(this->m_T+this->m_T)};

      switch (rd_enum) 
      {
            case sinusoidal_weighted_oqpsk_rand_distr::uniform : 
            {
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
                         const float cos_term0 = channel_I_sample_noise_n1_0_1(+0.0f,+1.0f)+(I_data0*cos_carrier0*cos_weight0);
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
                         const float sin_term0 = channel_Q_sample_noise_n1_0_1(+0.0f,+1.0f)+(Q_data0*sin_carrier0*sin_weight0);
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
                         const float cos_term1 = channel_I_sample_noise_n1_0_1(+0.0f,+1.0f)+(I_data1*cos_carrier1*cos_weight1);
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
                         const float sin_term1 = channel_Q_sample_noise_n1_0_1(+0.0f,+1.0f)+(Q_data1*sin_carrier1*sin_weight1);
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
                         const float cos_term2 = channel_I_sample_noise_n1_0_1(+0.0f,+1.0f)+(I_data2*cos_carrier2*cos_weight2);
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
                         const float sin_term2 = channel_Q_sample_noise_n1_0_1(+0.0f,+1.0f)+(Q_data2*sin_carrier2*sin_weight2);
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
                         const float cos_term3 = channel_I_sample_noise_n1_0_1(+0.0f,+1.0f)+(I_data3*cos_carrier3*cos_weight3);
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
                         const float sin_term3 = channel_Q_sample_noise_n1_0_1(+0.0f,+1.0f)+(Q_data3*sin_carrier3*sin_weight3);
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
                        const float cos_term0 = channel_I_sample_noise_n1_0_1(+0.0f,+1.0f)+(I_data0*cos_carrier0*cos_weight0);
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
                        const float sin_term0 = channel_Q_sample_noise_n1_0_1(+0.0f,+1.0f)+(Q_data0*sin_carrier0*sin_weight0);
                        this->m_Q_channel.m_data[j] = sin_term0;
                        this->m_msk_signal.m_data[j]= cos_term0+sin_term0;
                  }
                  // sanity check.
                  if(__builtin_expect(this->m_have_msk_samples==false,0)) { this->m_have_msk_samples = true;}
            }
            break;
            case sinusoidal_weighted_oqpsk_rand_distr::normal : 
            {
                  thread_local unsigned long long I_seed{0ull};
                  thread_local unsigned long long Q_seed{0ull};
                  thread_local std::normal_distribution<float> I_rand;
                  thread_local std::normal_distribution<float> Q_rand;
                  I_seed = __rdtsc();
                  I_rand = std::normal_distribution<float>(pdf_params.norm_mu_r,pdf_params.norm_sigma_r);
                  auto I_rengine{std::mt19937(I_seed)};
                  Q_seed = __rdtsc();
                  Q_rand = std::normal_distribution<float>(pdf_params.norm_mu_r,pdf_params.norm_sigma_r);
                  auto Q_rengine{std::mt19937(Q_seed)};
                  
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
                         const float cos_term0 = I_rand.operator()(I_rengine)+(I_data0*cos_carrier0*cos_weight0);
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
                         const float sin_term0 = Q_rand.operator()(Q_rengine)+(Q_data0*sin_carrier0*sin_weight0);
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
                         const float cos_term1 = I_rand.operator()(I_rengine)+(I_data1*cos_carrier1*cos_weight1);
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
                         const float sin_term1 = Q_rand.operator()(Q_rengine)+(Q_data1*sin_carrier1*sin_weight1);
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
                         const float cos_term2 = I_rand.operator()(I_rengine)+(I_data2*cos_carrier2*cos_weight2);
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
                         const float sin_term2 = Q_rand.operator()(Q_rengine)+(Q_data2*sin_carrier2*sin_weight2);
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
                         const float cos_term3 = I_rand.operator()(I_rengine)+(I_data3*cos_carrier3*cos_weight3);
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
                         const float sin_term3 = Q_rand.operator()(Q_rengine)+(Q_data3*sin_carrier3*sin_weight3);
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
                        const float cos_term0 = I_rand.operator()(I_rengine)+(I_data0*cos_carrier0*cos_weight0);
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
                        const float sin_term0 = Q_rand.operator()(Q_rengine)+(Q_data0*sin_carrier0*sin_weight0);
                        this->m_Q_channel.m_data[j] = sin_term0;
                        this->m_msk_signal.m_data[j]= cos_term0+sin_term0;
                  }
                  // sanity check.
                  if(__builtin_expect(this->m_have_msk_samples==false,0)) { this->m_have_msk_samples = true;}
            }
            break;
            case sinusoidal_weighted_oqpsk_rand_distr::cauchy : 
            {
                  thread_local unsigned long long I_seed{0ull};
                  thread_local unsigned long long Q_seed{0ull};
                  thread_local std::cauchy_distribution<float> I_rand;
                  thread_local std::cauchy_distribution<float> Q_rand;
                  I_seed = __rdtsc();
                  I_rand = std::cauchy_distribution<float>(pdf_params.cauchy_a_r,pdf_params.cauchy_b_r);
                  auto I_rengine{std::mt19937(I_seed)};
                  Q_seed = __rdtsc();
                  Q_rand = std::cauchy_distribution<float>(pdf_params.cauchy_a_r,pdf_params.cauchy_b_r);
                  auto Q_rengine{std::mt19937(Q_seed)};

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
                         const float cos_term0 = I_rand.operator()(I_rengine)+(I_data0*cos_carrier0*cos_weight0);
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
                         const float sin_term0 = Q_rand.operator()(Q_rengine)+(Q_data0*sin_carrier0*sin_weight0);
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
                         const float cos_term1 = I_rand.operator()(I_rengine)+(I_data1*cos_carrier1*cos_weight1);
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
                         const float sin_term1 = Q_rand.operator()(Q_rengine)+(Q_data1*sin_carrier1*sin_weight1);
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
                         const float cos_term2 = I_rand.operator()(I_rengine)+(I_data2*cos_carrier2*cos_weight2);
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
                         const float sin_term2 = Q_rand.operator()(Q_rengine)+(Q_data2*sin_carrier2*sin_weight2);
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
                         const float cos_term3 = I_rand.operator()(I_rengine)+(I_data3*cos_carrier3*cos_weight3);
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
                         const float sin_term3 = Q_rand.operator()(Q_rengine)+(Q_data3*sin_carrier3*sin_weight3);
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
                        const float cos_term0 = I_rand.operator()(I_rengine)+(I_data0*cos_carrier0*cos_weight0);
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
                        const float sin_term0 = Q_rand.operator()(Q_rengine)+(Q_data0*sin_carrier0*sin_weight0);
                        this->m_Q_channel.m_data[j] = sin_term0;
                        this->m_msk_signal.m_data[j]= cos_term0+sin_term0;
                  }
                  // sanity check.
                  if(__builtin_expect(this->m_have_msk_samples==false,0)) { this->m_have_msk_samples = true;}
            }
            break;
            case sinusoidal_weighted_oqpsk_rand_distr::log_norm : 
            {
                  thread_local unsigned long long I_seed{0ull};
                  thread_local unsigned long long Q_seed{0ull};
                  thread_local std::lognormal_distribution<float> I_rand;
                  thread_local std::lognormal_distribution<float> Q_rand;
                  I_seed = __rdtsc();
                  I_rand = std::lognormal_distribution<float>(pdf_params.log_norm_m_r,pdf_params.log_norm_s_r);
                  auto I_rengine{std::mt19937(I_seed)};
                  Q_seed = __rdtsc();
                  Q_rand = std::lognormal_distribution<float>(pdf_params.log_norm_m_r,pdf_params.log_norm_s_r);
                  auto Q_rengine{std::mt19937(Q_seed)};

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
                         const float cos_term0 = I_rand.operator()(I_rengine)+(I_data0*cos_carrier0*cos_weight0);
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
                         const float sin_term0 = Q_rand.operator()(Q_rengine)+(Q_data0*sin_carrier0*sin_weight0);
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
                         const float cos_term1 = I_rand.operator()(I_rengine)+(I_data1*cos_carrier1*cos_weight1);
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
                         const float sin_term1 = Q_rand.operator()(Q_rengine)+(Q_data1*sin_carrier1*sin_weight1);
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
                         const float cos_term2 = I_rand.operator()(I_rengine)+(I_data2*cos_carrier2*cos_weight2);
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
                         const float sin_term2 = Q_rand.operator()(Q_rengine)+(Q_data2*sin_carrier2*sin_weight2);
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
                         const float cos_term3 = I_rand.operator()(I_rengine)+(I_data3*cos_carrier3*cos_weight3);
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
                         const float sin_term3 = Q_rand.operator()(Q_rengine)+(Q_data3*sin_carrier3*sin_weight3);
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
                        const float cos_term0 = I_rand.operator()(I_rengine)+(I_data0*cos_carrier0*cos_weight0);
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
                        const float sin_term0 = Q_rand.operator()(Q_rengine)+(Q_data0*sin_carrier0*sin_weight0);
                        this->m_Q_channel.m_data[j] = sin_term0;
                        this->m_msk_signal.m_data[j]= cos_term0+sin_term0;
                  }
                  // sanity check.
                  if(__builtin_expect(this->m_have_msk_samples==false,0)) { this->m_have_msk_samples = true;}
            }
            break;
            case sinusoidal_weighted_oqpsk_rand_distr::expo_gamma : 
            {
                  thread_local unsigned long long I_seed{0ull};
                  thread_local unsigned long long Q_seed{0ull};
                  thread_local std::exponential_distribution<float> I_rand;
                  thread_local std::exponential_distribution<float> Q_rand;
                  I_seed = __rdtsc();
                  I_rand = std::exponential_distribution<float>(pdf_params.expo_gamma_r);
                  auto I_rengine{std::mt19937(I_seed)};
                  Q_seed = __rdtsc();
                  Q_rand = std::exponential_distribution<float>(pdf_params.expo_gamma_r);
                  auto Q_rengine{std::mt19937(Q_seed)};

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
                         const float cos_term0 = I_rand.operator()(I_rengine)+(I_data0*cos_carrier0*cos_weight0);
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
                         const float sin_term0 = Q_rand.operator()(Q_rengine)+(Q_data0*sin_carrier0*sin_weight0);
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
                         const float cos_term1 = I_rand.operator()(I_rengine)+(I_data1*cos_carrier1*cos_weight1);
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
                         const float sin_term1 = Q_rand.operator()(Q_rengine)+(Q_data1*sin_carrier1*sin_weight1);
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
                         const float cos_term2 = I_rand.operator()(I_rengine)+(I_data2*cos_carrier2*cos_weight2);
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
                         const float sin_term2 = Q_rand.operator()(Q_rengine)+(Q_data2*sin_carrier2*sin_weight2);
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
                         const float cos_term3 = I_rand.operator()(I_rengine)+(I_data3*cos_carrier3*cos_weight3);
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
                         const float sin_term3 = Q_rand.operator()(Q_rengine)+(Q_data3*sin_carrier3*sin_weight3);
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
                        const float cos_term0 = I_rand.operator()(I_rengine)+(I_data0*cos_carrier0*cos_weight0);
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
                        const float sin_term0 = Q_rand.operator()(Q_rengine)+(Q_data0*sin_carrier0*sin_weight0);
                        this->m_Q_channel.m_data[j] = sin_term0;
                        this->m_msk_signal.m_data[j]= cos_term0+sin_term0;
                  }
                  // sanity check.
                  if(__builtin_expect(this->m_have_msk_samples==false,0)) { this->m_have_msk_samples = true;}
            }
            break;
            case sinusoidal_weighted_oqpsk_rand_distr::weibull : 
            {
                  thread_local unsigned long long I_seed{0ull};
                  thread_local unsigned long long Q_seed{0ull};
                  thread_local std::weibull_distribution<float> I_rand;
                  thread_local std::weibull_distribution<float> Q_rand;
                  I_seed = __rdtsc();
                  I_rand = std::weibull_distribution<float>(pdf_params.weibull_a_r,pdf_params.weibull_b_r);
                  auto I_rengine{std::mt19937(I_seed)};
                  Q_seed = __rdtsc();
                  Q_rand = std::weibull_distribution<float>(pdf_params.weibull_a_r,pdf_params.weibull_b_r);
                  auto Q_rengine{std::mt19937(Q_seed)};

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
                         const float cos_term0 = I_rand.operator()(I_rengine)+(I_data0*cos_carrier0*cos_weight0);
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
                         const float sin_term0 = Q_rand.operator()(Q_rengine)+(Q_data0*sin_carrier0*sin_weight0);
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
                         const float cos_term1 = I_rand.operator()(I_rengine)+(I_data1*cos_carrier1*cos_weight1);
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
                         const float sin_term1 = Q_rand.operator()(Q_rengine)+(Q_data1*sin_carrier1*sin_weight1);
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
                         const float cos_term2 = I_rand.operator()(I_rengine)+(I_data2*cos_carrier2*cos_weight2);
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
                         const float sin_term2 = Q_rand.operator()(Q_rengine)+(Q_data2*sin_carrier2*sin_weight2);
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
                         const float cos_term3 = I_rand.operator()(I_rengine)+(I_data3*cos_carrier3*cos_weight3);
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
                         const float sin_term3 = Q_rand.operator()(Q_rengine)+(Q_data3*sin_carrier3*sin_weight3);
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
                        const float cos_term0 = I_rand.operator()(I_rengine)+(I_data0*cos_carrier0*cos_weight0);
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
                        const float sin_term0 = Q_rand.operator()(Q_rengine)+(Q_data0*sin_carrier0*sin_weight0);
                        this->m_Q_channel.m_data[j] = sin_term0;
                        this->m_msk_signal.m_data[j]= cos_term0+sin_term0;
                  }
                  // sanity check.
                  if(__builtin_expect(this->m_have_msk_samples==false,0)) { this->m_have_msk_samples = true;}
            }
            break;
            case sinusoidal_weighted_oqpsk_rand_distr::gamma : 
            {
                  thread_local unsigned long long I_seed{0ull};
                  thread_local unsigned long long Q_seed{0ull};
                  thread_local std::gamma_distribution<float> I_rand;
                  thread_local std::gamma_distribution<float> Q_rand;
                  I_seed = __rdtsc();
                  I_rand = std::gamma_distribution<float>(pdf_params.gamma_alph_r,pdf_params.gamma_bet_r);
                  auto I_rengine{std::mt19937(I_seed)};
                  Q_seed = __rdtsc();
                  Q_rand = std::gamma_distribution<float>(pdf_params.gamma_alph_r,pdf_params.gamma_bet_r);
                  auto Q_rengine{std::mt19937(Q_seed)};

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
                         const float cos_term0 = I_rand.operator()(I_rengine)+(I_data0*cos_carrier0*cos_weight0);
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
                         const float sin_term0 = Q_rand.operator()(Q_rengine)+(Q_data0*sin_carrier0*sin_weight0);
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
                         const float cos_term1 = I_rand.operator()(I_rengine)+(I_data1*cos_carrier1*cos_weight1);
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
                         const float sin_term1 = Q_rand.operator()(Q_rengine)+(Q_data1*sin_carrier1*sin_weight1);
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
                         const float cos_term2 = I_rand.operator()(I_rengine)+(I_data2*cos_carrier2*cos_weight2);
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
                         const float sin_term2 = Q_rand.operator()(Q_rengine)+(Q_data2*sin_carrier2*sin_weight2);
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
                         const float cos_term3 = I_rand.operator()(I_rengine)+(I_data3*cos_carrier3*cos_weight3);
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
                         const float sin_term3 = Q_rand.operator()(Q_rengine)+(Q_data3*sin_carrier3*sin_weight3);
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
                        const float cos_term0 = I_rand.operator()(I_rengine)+(I_data0*cos_carrier0*cos_weight0);
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
                        const float sin_term0 = Q_rand.operator()(Q_rengine)+(Q_data0*sin_carrier0*sin_weight0);
                        this->m_Q_channel.m_data[j] = sin_term0;
                        this->m_msk_signal.m_data[j]= cos_term0+sin_term0;
                  }
                  // sanity check.
                  if(__builtin_expect(this->m_have_msk_samples==false,0)) { this->m_have_msk_samples = true;}
            }
            break;
            default : 
                return (-3);
      }

      return (0);
}

std::int32_t 
gms::radiolocation
::sinusoidal_weighted_oqpsk_t
::fourier_transform_iq_shape_pulse_u4x(const darray_r4_t &frequencies)                            
{
      using namespace gms::math;
      if(__builtin_expect(frequencies.mnx!=this->m_nfrequencies,0)) { return (-1);}
      constexpr float C010132118364233777144387946321{0.10132118364233777144387946321f};
      constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
      const     float fourT{4.0f*this->m_T};
      const     float sqrA{this->m_A_I*this->m_A_I}; //same for both channels
      const     float lead_term{16.0f*sqrA*this->m_T*C010132118364233777144387946321};
      const     std::size_t n_freqs{frequencies.mnx};
      std::size_t i,j;
      for(i = 0ull; i != ROUND_TO_FOUR(n_freqs,4ull); i += 4ull) 
      {
          const float freq_sample0{frequencies.m_data[i+0ull]};
#if (SINUSOIDAL_WEIGHTED_OQPSK_USE_CEPHES) == 1 
          const float cos_sample0{ceph_cosf(C6283185307179586476925286766559*this->m_T*freq_sample0)};
#else 
          const float cos_sample0{std::cos(C6283185307179586476925286766559*this->m_T*freq_sample0)};
#endif 
          const float t0{fourT*freq_sample0};
          const float denom0{1.0f-(t0*t0)};
          const float cos_term0{cos_sample0/denom0};
          const float FT_sample0{lead_term*(cos_term0*cos_term0)};
          this->m_Psi_f.m_data[i+0ull] = FT_sample0;
          const float freq_sample1{frequencies.m_data[i+1ull]};
#if (SINUSOIDAL_WEIGHTED_OQPSK_USE_CEPHES) == 1 
          const float cos_sample1{ceph_cosf(C6283185307179586476925286766559*this->m_T*freq_sample1)};
#else 
          const float cos_sample1{std::cos(C6283185307179586476925286766559*this->m_T*freq_sample1)};
#endif 
          const float t1{fourT*freq_sample1};
          const float denom1{1.0f-(t1*t1)};
          const float cos_term1{cos_sample1/denom1};
          const float FT_sample1{lead_term*(cos_term1*cos_term1)};
          this->m_Psi_f.m_data[i+1ull] = FT_sample1;
          const float freq_sample2{frequencies.m_data[i+2ull]};
#if (SINUSOIDAL_WEIGHTED_OQPSK_USE_CEPHES) == 1 
          const float cos_sample2{ceph_cosf(C6283185307179586476925286766559*this->m_T*freq_sample2)};
#else 
          const float cos_sample2{std::cos(C6283185307179586476925286766559*this->m_T*freq_sample2)};
#endif 
          const float t2{fourT*freq_sample2};
          const float denom2{1.0f-(t0*t0)};
          const float cos_term2{cos_sample2/denom2};
          const float FT_sample2{lead_term*(cos_term2*cos_term2)};
          this->m_Psi_f.m_data[i+2ull] = FT_sample2;
          const float freq_sample3{frequencies.m_data[i+3ull]};
#if (SINUSOIDAL_WEIGHTED_OQPSK_USE_CEPHES) == 1 
          const float cos_sample3{ceph_cosf(C6283185307179586476925286766559*this->m_T*freq_sample3)};
#else 
          const float cos_sample3{std::cos(C6283185307179586476925286766559*this->m_T*freq_sample3)};
#endif 
          const float t3{fourT*freq_sample3};
          const float denom3{1.0f-(t3*t3)};
          const float cos_term3{cos_sample3/denom3};
          const float FT_sample3{lead_term*(cos_term3*cos_term3)};
          this->m_Psi_f.m_data[i+3ull] = FT_sample3;
      }

      for(j = i; j != n_freqs; ++j)  
      {
          const float freq_sample0{frequencies.m_data[j]};
#if (SINUSOIDAL_WEIGHTED_OQPSK_USE_CEPHES) == 1 
          const float cos_sample0{ceph_cosf(C6283185307179586476925286766559*this->m_T*freq_sample0)};
#else 
          const float cos_sample0{std::cos(C6283185307179586476925286766559*this->m_T*freq_sample0)};
#endif 
          const float t0{fourT*freq_sample0};
          const float denom0{1.0f-(t0*t0)};
          const float cos_term0{cos_sample0/denom0};
          const float FT_sample0{lead_term*(cos_term0*cos_term0)};
          this->m_Psi_f.m_data[j] = FT_sample0; 
      }

      return (0);
} 

auto
gms::radiolocation
::operator<<(std::ostream &os,
             gms::radiolocation::sinusoidal_weighted_oqpsk_t &rhs)->std::ostream &
{
    std::cout << typeid(rhs).name() << "Begin: object state dump." << std::endl;
    std::cout << "m_nsamples         : "      << rhs.m_nsamples       << std::endl;
    std::cout << "m_I_ch_nsamples    : "      << rhs.m_I_ch_nsamples  << std::endl;
    std::cout << "m_Q_ch_nsamples    : "      << rhs.m_Q_ch_nsamples  << std::endl;
    std::cout << "m_nfrequencies     : "      << rhs.m_nfrequencies   << std::endl;
    std::cout << "m_A_I              : "      << std::fixed << std::setprecision(7) << rhs.m_A_I << std::endl;
    std::cout << "m_A_Q              : "      << std::fixed << std::setprecision(7) << rhs.m_A_Q << std::endl;
    std::cout << "m_cw0              : "      << std::fixed << std::setprecision(7) << rhs.m_cw0 << std::endl;
    std::cout << "m_sw0              : "      << std::fixed << std::setprecision(7) << rhs.m_sw0 << std::endl;
    std::cout << "m_cph0             : "      << std::fixed << std::setprecision(7) << rhs.m_cph0 << std::endl;
    std::cout << "m_sph0             : "      << std::fixed << std::setprecision(7) << rhs.m_sph0 << std::endl;
    std::cout << std::boolalpha << std::endl;
    std::cout << "m_have_msk_samples : "      << rhs.m_have_msk_samples << std::endl;    
    std::cout << "OQPSK[MSK] -- Samples:" << "\n\n";
    for(std::size_t i{0ull}; i != rhs.m_nsamples; ++i)
    {
        os << std::fixed << std::setprecision(7) << rhs.m_msk_signal.m_data[i] << std::endl;
    
    }
    std::cout << "I-Bitstream -- Samples:" << "\n\n";
    for(std::size_t i{0ull}; i != rhs.m_I_ch_nsamples; ++i) 
    {
       os << std::fixed << std::setprecision(7) << rhs.m_I_bitstream.m_data[i] << std::endl;
    }
    std::cout << "Q-Bitstream -- Samples:" << "\n\n";
    for(std::size_t i{0ull}; i != rhs.m_Q_ch_nsamples; ++i) 
    {
       os << std::fixed << std::setprecision(7) << rhs.m_Q_bitstream.m_data[i] << std::endl;
    }
    std::cout << "I-Channel -- Samples:" << "\n\n";
    for(std::size_t i{0ull}; i != rhs.m_nsamples; ++i) 
    {
       os << std::fixed << std::setprecision(7) << rhs.m_I_channel.m_data[i] << std::endl;
    }
    std::cout << "Q-Channel -- Samples:" << "\n\n";
    for(std::size_t i{0ull}; i != rhs.m_nsamples; ++i) 
    {
       os << std::fixed << std::setprecision(7) << rhs.m_Q_channel.m_data[i] << std::endl;
    }
    std::cout << "Fourier-Transform -- I/Q shaping pulses:" << "\n\n";
    for(std::size_t i{0ull}; i != rhs.m_nsamples; ++i) 
    {
       os << std::fixed << std::setprecision(7) << rhs.m_Psi_f.m_data[i] << std::endl;
    }
    std::cout << typeid(rhs).name() << "End: object state dump." << std::endl;
    return (os);
}
                              