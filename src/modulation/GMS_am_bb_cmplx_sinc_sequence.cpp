
#include <immintrin.h>
#include <fstream>
#include <iomanip>
#include "GMS_am_bb_cmplx_sinc_sequence.h"

gms::radiolocation
::am_bb_cmplx_sinc_signal_sequence_t
::am_bb_cmplx_sinc_signal_sequence_t(  const std::size_t nsamples,
                                       const std::size_t nsignals,
                                       const std::vector<std::uint32_t> &nK_values,
                                       const std::vector<float> &n_values,
                                       const std::vector<float> &A_values,
                                       const std::vector<float> &P_values)
:
m_nsamples{nsamples},
m_nsignals{nsignals},
m_nK_values(std::move(nK_values)),
m_n_values(std::move(n_values)),
m_A_values(std::move(A_values)),
m_P_values(std::move(P_values)),
m_bb_cmplx_sinc_signals(std::vector<am_bb_cmplx_sinc_signal_t>(this->m_nsignals)),
m_re_sequence(darray_r4_t(this->m_nsamples)),
m_im_sequence(darray_r4_t(this->m_nsamples))
{
    for(std::size_t i{0ull}; i != this->m_nsignals; ++i) 
    {
         const std::uint32_t nK_vals{this->m_nK_values.operator[](i)};
         const float n_vals{this->m_n_values.operator[](i)};
         const float A_vals{this->m_A_values.operator[](i)};
         const float P_vals{this->m_P_values.operator[](i)};
         this->m_bb_cmplx_sinc_signals.operator[](i) = am_bb_cmplx_sinc_signal_t(this->m_nsamples,
                                                                                nK_vals,
                                                                                n_vals,A_vals,P_vals);
    }
}

gms::radiolocation
::am_bb_cmplx_sinc_signal_sequence_t
::am_bb_cmplx_sinc_signal_sequence_t(am_bb_cmplx_sinc_signal_sequence_t &&other) noexcept(true) 
:
m_nsamples{std::move(other.m_nsamples)},
m_nsignals{std::move(other.m_nsignals)},
m_nK_values{std::move(other.m_nK_values)},
m_n_values{std::move(other.m_n_values)},
m_A_values{std::move(other.m_A_values)},
m_P_values{std::move(other.m_P_values)},
m_bb_cmplx_sinc_signals{std::move(other.m_bb_cmplx_sinc_signals)},
m_re_sequence{std::move(other.m_re_sequence)},
m_im_sequence{std::move(other.m_im_sequence)}
{

}

gms::radiolocation
::am_bb_cmplx_sinc_signal_sequence_t
::am_bb_cmplx_sinc_signal_sequence_t(const am_bb_cmplx_sinc_signal_sequence_t &other) noexcept(false)
:
m_nsamples{other.m_nsamples},
m_nsignals{other.m_nsignals},
m_nK_values{other.m_nK_values},
m_n_values{other.m_n_values},
m_A_values{other.m_A_values},
m_P_values{other.m_P_values},
m_bb_cmplx_sinc_signals{other.m_bb_cmplx_sinc_signals},
m_re_sequence{other.m_re_sequence},
m_im_sequence{other.m_im_sequence}
{

}

gms::radiolocation
::am_bb_cmplx_sinc_signal_sequence_t
::~am_bb_cmplx_sinc_signal_sequence_t() noexcept(false) 
{

}

gms::radiolocation
::am_bb_cmplx_sinc_signal_sequence_t &
gms::radiolocation
::am_bb_cmplx_sinc_signal_sequence_t
::operator=(const am_bb_cmplx_sinc_signal_sequence_t &other) noexcept(false) 
{
      if(__builtin_expect(this==&other,0)) { return (*this);}
      this->m_nsamples         = other.m_nsamples;
      this->m_nsignals         = other.m_nsignals;
      this->m_nK_values.operator=(other.m_nK_values);
      this->m_n_values.operator=(other.m_n_values);
      this->m_A_values.operator=(other.m_A_values);
      this->m_P_values.operator=(other.m_P_values);
      this->m_bb_cmplx_sinc_signals.operator=(other.m_bb_cmplx_sinc_signals);
      this->m_re_sequence.operator=(other.m_re_sequence);
      this->m_im_sequence.operator=(other.m_im_sequence);
      return (*this);
}

gms::radiolocation
::am_bb_cmplx_sinc_signal_sequence_t &
gms::radiolocation 
::am_bb_cmplx_sinc_signal_sequence_t 
::operator=(am_bb_cmplx_sinc_signal_sequence_t &&other) noexcept(true) 
{
      if(__builtin_expect(this==&other,0)) { return (*this);}
      this->m_nsamples         = std::move(other.m_nsamples);
      this->m_nsignals         = std::move(other.m_nsignals);
      this->m_nK_values.operator=(std::move(other.m_nK_values));
      this->m_n_values.operator=(std::move(other.m_n_values));
      this->m_A_values.operator=(std::move(other.m_A_values));
      this->m_P_values.operator=(std::move(other.m_P_values));
      this->m_bb_cmplx_sinc_signals.operator=(std::move(other.m_bb_cmplx_sinc_signals));
      this->m_re_sequence.operator=(std::move(other.m_re_sequence));
      this->m_im_sequence.operator=(std::move(other.m_im_sequence));
      return (*this);
}

void 
gms::radiolocation 
::am_bb_cmplx_sinc_signal_sequence_t 
::create_sequence_plot(const std::uint32_t n_samp,
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
::am_bb_cmplx_sinc_signal_sequence_t
::signal_sequence_sse42_u16x(const float * __restrict__           sym_in ,
                             am_bb_cmplx_sinc_signal_pdf_params_t &pdf_params,
                             am_bb_cmplx_sinc_signal_rand_distr    rd_enum,
                             const float          * __restrict__  noise_scales,
                             std::int32_t * __restrict__          ret_code,
                             const std::uint32_t                  which_generator)
{
        if(__builtin_expect(nullptr==sym_in,0)) {return (-1);}
        if(__builtin_expect(this->m_nsignals>256ull,0)) {return (-2);}
        if(__builtin_expect((this->m_nsamples%2ull)!=0ull,0)) {return (-3);}

#if (AM_BB_CMPLX_SINC_SEQUENCE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif
        constexpr std::size_t max_buff_size{256ull};
        __ATTR_ALIGN__(16)
        __m128 xmm0_r[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm1_r[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm2_r[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm3_r[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm4_r[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm5_r[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm6_r[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm7_r[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm8_r[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm9_r[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm10_r[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm11_r[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm12_r[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm13_r[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm14_r[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm15_r[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm0_i[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm1_i[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm2_i[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm3_i[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm4_i[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm5_i[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm6_i[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm7_i[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm8_i[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm9_i[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm10_i[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm11_i[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm12_i[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm13_i[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm14_i[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm15_i[max_buff_size];
         float scal_rem_r[max_buff_size];
         float scal_rem_i[max_buff_size];
        __m128 vzero{_mm_setzero_ps()};
        __m128 sum0_r{vzero};
        __m128 sum1_r{sum0_r};
        __m128 sum2_r{sum0_r};
        __m128 sum3_r{sum0_r};
        __m128 sum4_r{sum0_r};
        __m128 sum5_r{sum0_r};
        __m128 sum6_r{sum0_r};
        __m128 sum7_r{sum0_r};
        __m128 sum8_r{sum0_r};
        __m128 sum9_r{sum0_r};
        __m128 sum10_r{sum0_r};
        __m128 sum11_r{sum0_r};
        __m128 sum12_r{sum0_r};
        __m128 sum13_r{sum0_r};
        __m128 sum14_r{sum0_r};
        __m128 sum15_r{sum0_r};
        __m128 sum0_i{sum0_r};
        __m128 sum1_i{sum0_r};
        __m128 sum2_i{sum0_r};
        __m128 sum3_i{sum0_r};
        __m128 sum4_i{sum0_r};
        __m128 sum5_i{sum0_r};
        __m128 sum6_i{sum0_r};
        __m128 sum7_i{sum0_r};
        __m128 sum8_i{sum0_r};
        __m128 sum9_i{sum0_r};
        __m128 sum10_i{sum0_r};
        __m128 sum11_i{sum0_r};
        __m128 sum12_i{sum0_r};
        __m128 sum13_i{sum0_r};
        __m128 sum14_i{sum0_r};
        __m128 sum15_i{sum0_r}; 
        std::size_t i;
        float sum_rem_r{0.0f};
        float sum_rem_i{0.0f};
        constexpr std::size_t exceeds_L1D_size{32768u};
        const     std::size_t n_samples_halved{this->m_nsamples/2ull};

        switch(which_generator) 
        {
             case 0 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_signal_user_data(sym_in,
                                                                              static_cast<std::uint32_t>(this->m_nsamples),
                                                                              this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 1 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_signal_user_data_u4x(sym_in,
                                                                              static_cast<std::uint32_t>(this->m_nsamples),
                                                                              this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 2 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {    
                        const float scale{noise_scales[k]};
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_noisy_signal_user_data(pdf_params,
                                                                                      scale,sym_in,
                                                                                      static_cast<std::uint32_t>(this->m_nsamples),
                                                                                      this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 3 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {    
                        const float scale{noise_scales[k]};
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_noisy_signal_user_data(pdf_params,
                                                                                      rd_enum,scale,sym_in,
                                                                                      static_cast<std::uint32_t>(this->m_nsamples),
                                                                                      this->m_nK_values.operator[](k));
                   } 
             }
             break;
             default : 
                   return (-3);
        }

        if(this->m_nsamples > exceeds_L1D_size)
        {
             for(i = 0ull; (i+63ull) < n_samples_halved; i += 64ull) 
             {
                 sum0_r = vzero;
                 sum0_i = vzero;
                 sum1_r = vzero;
                 sum1_i = vzero;
                 sum2_r = vzero;
                 sum2_i = vzero;
                 sum3_r = vzero;
                 sum3_i = vzero;
                 sum4_r = vzero;
                 sum4_i = vzero;
                 sum5_r = vzero;
                 sum5_i = vzero;
                 sum6_r = vzero;
                 sum6_i = vzero;
                 sum7_r = vzero;
                 sum7_i = vzero;
                 sum8_r = vzero;
                 sum8_i = vzero;
                 sum9_r = vzero;
                 sum9_i = vzero;
                 sum10_r = vzero;
                 sum10_i = vzero;
                 sum11_r = vzero;
                 sum11_i = vzero;
                 sum12_r = vzero;
                 sum12_i = vzero;
                 sum13_r = vzero;
                 sum13_i = vzero;
                 sum14_r = vzero;
                 sum14_i = vzero;
                 sum15_r = vzero;
                 sum15_i = vzero;
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
#if (AM_BB_CMPLX_SINC_SEQUENCE_SOFT_PREFETCH) == 1    
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull],_MM_HINT_T1);
                     
#endif 
                      xmm0_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm_add_ps(xmm0_r[j],sum0_r);
                      xmm0_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                      xmm1_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+4ull]);
                      sum1_r    = _mm_add_ps(xmm1_r[j],sum1_r);
                      xmm1_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+4ull]);
                      sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                      xmm2_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum2_r    = _mm_add_ps(xmm2_r[j],sum2_r);
                      xmm2_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                      xmm3_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+12ull]);
                      sum3_r    = _mm_add_ps(xmm3_r[j],sum3_r);
                      xmm3_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+12ull]);
                      sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                      xmm4_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum4_r    = _mm_add_ps(xmm4_r[j],sum4_r);
                      xmm4_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);
                      xmm5_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+20ull]);
                      sum5_r    = _mm_add_ps(xmm5_r[j],sum5_r);
                      xmm5_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+20ull]);
                      sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);
                      xmm6_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                      sum6_r    = _mm_add_ps(xmm6_r[j],sum6_r);
                      xmm6_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                      sum6_i    = _mm_add_ps(xmm6_i[j],sum6_i);
                      xmm7_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+28ull]);
                      sum7_r    = _mm_add_ps(xmm7_r[j],sum7_r);
                      xmm7_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+28ull]);
                      sum7_i    = _mm_add_ps(xmm7_i[j],sum7_i);
                      xmm8_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum8_r    = _mm_add_ps(xmm8_r[j],sum8_r);
                      xmm8_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum8_i    = _mm_add_ps(xmm8_i[j],sum8_i);
                      xmm9_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+36ull]);
                      sum9_r    = _mm_add_ps(xmm9_r[j],sum9_r);
                      xmm9_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+36ull]);
                      sum9_i    = _mm_add_ps(xmm9_i[j],sum9_i);
                      xmm10_r[j]= _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+40ull]);
                      sum10_r   = _mm_add_ps(xmm10_r[j],sum10_r);
                      xmm10_i[j]= _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+40ull]);
                      sum10_i   = _mm_add_ps(xmm10_i[j],sum10_i);
                      xmm11_r[j]= _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+44ull]);
                      sum11_r    = _mm_add_ps(xmm11_r[j],sum11_r);
                      xmm11_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+44ull]);
                      sum11_i    = _mm_add_ps(xmm11_i[j],sum11_i);
                      xmm12_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum12_r    = _mm_add_ps(xmm12_r[j],sum12_r);
                      xmm12_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum12_i    = _mm_add_ps(xmm12_i[j],sum12_i);
                      xmm13_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+52ull]);
                      sum13_r    = _mm_add_ps(xmm13_r[j],sum13_r);
                      xmm13_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+52ull]);
                      sum13_i    = _mm_add_ps(xmm13_i[j],sum13_i);
                      xmm14_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+56ull]);
                      sum14_r    = _mm_add_ps(xmm14_r[j],sum14_r);
                      xmm14_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+56ull]);
                      sum14_i    = _mm_add_ps(xmm14_i[j],sum14_i);
                      xmm15_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+60ull]);
                      sum15_r    = _mm_add_ps(xmm15_r[j],sum15_r);
                      xmm15_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+60ull]);
                      sum15_i    = _mm_add_ps(xmm15_i[j],sum15_i);

                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+4ull],sum1_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+4ull],sum1_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum2_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum2_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+12ull],sum3_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+12ull],sum3_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum4_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum4_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+20ull],sum5_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+20ull],sum5_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum6_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum6_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+28ull],sum7_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+28ull],sum7_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum8_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum8_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+36ull],sum9_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+36ull],sum9_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+40ull],sum10_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+40ull],sum10_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+44ull],sum11_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+44ull],sum11_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum12_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum12_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+52ull],sum13_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+52ull],sum13_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+56ull],sum14_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+56ull],sum14_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+60ull],sum15_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+60ull],sum15_i);
                 }

             }
        }
        else 
        {
                for(i = 0ull; (i+63ull) < n_samples_halved; i += 64ull) 
                {
                   sum0_r = vzero;
                   sum0_i = vzero;
                   sum1_r = vzero;
                   sum1_i = vzero;
                   sum2_r = vzero;
                   sum2_i = vzero;
                   sum3_r = vzero;
                   sum3_i = vzero;
                   sum4_r = vzero;
                   sum4_i = vzero;
                   sum5_r = vzero;
                   sum5_i = vzero;
                   sum6_r = vzero;
                   sum6_i = vzero;
                   sum7_r = vzero;
                   sum7_i = vzero;
                   sum8_r = vzero;
                   sum8_i = vzero;
                   sum9_r = vzero;
                   sum9_i = vzero;
                   sum10_r = vzero;
                   sum10_i = vzero;
                   sum11_r = vzero;
                   sum11_i = vzero;
                   sum12_r = vzero;
                   sum12_i = vzero;
                   sum13_r = vzero;
                   sum13_i = vzero;
                   sum14_r = vzero;
                   sum14_i = vzero;
                   sum15_r = vzero;
                   sum15_i = vzero;
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                      xmm0_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm_add_ps(xmm0_r[j],sum0_r);
                      xmm0_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                      xmm1_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+4ull]);
                      sum1_r    = _mm_add_ps(xmm1_r[j],sum1_r);
                      xmm1_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+4ull]);
                      sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                      xmm2_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum2_r    = _mm_add_ps(xmm2_r[j],sum2_r);
                      xmm2_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                      xmm3_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+12ull]);
                      sum3_r    = _mm_add_ps(xmm3_r[j],sum3_r);
                      xmm3_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+12ull]);
                      sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                      xmm4_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum4_r    = _mm_add_ps(xmm4_r[j],sum4_r);
                      xmm4_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);
                      xmm5_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+20ull]);
                      sum5_r    = _mm_add_ps(xmm5_r[j],sum5_r);
                      xmm5_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+20ull]);
                      sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);
                      xmm6_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                      sum6_r    = _mm_add_ps(xmm6_r[j],sum6_r);
                      xmm6_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                      sum6_i    = _mm_add_ps(xmm6_i[j],sum6_i);
                      xmm7_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+28ull]);
                      sum7_r    = _mm_add_ps(xmm7_r[j],sum7_r);
                      xmm7_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+28ull]);
                      sum7_i    = _mm_add_ps(xmm7_i[j],sum7_i);
                      xmm8_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum8_r    = _mm_add_ps(xmm8_r[j],sum8_r);
                      xmm8_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum8_i    = _mm_add_ps(xmm8_i[j],sum8_i);
                      xmm9_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+36ull]);
                      sum9_r    = _mm_add_ps(xmm9_r[j],sum9_r);
                      xmm9_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+36ull]);
                      sum9_i    = _mm_add_ps(xmm9_i[j],sum9_i);
                      xmm10_r[j]= _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+40ull]);
                      sum10_r   = _mm_add_ps(xmm10_r[j],sum10_r);
                      xmm10_i[j]= _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+40ull]);
                      sum10_i   = _mm_add_ps(xmm10_i[j],sum10_i);
                      xmm11_r[j]= _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+44ull]);
                      sum11_r    = _mm_add_ps(xmm11_r[j],sum11_r);
                      xmm11_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+44ull]);
                      sum11_i    = _mm_add_ps(xmm11_i[j],sum11_i);
                      xmm12_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum12_r    = _mm_add_ps(xmm12_r[j],sum12_r);
                      xmm12_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum12_i    = _mm_add_ps(xmm12_i[j],sum12_i);
                      xmm13_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+52ull]);
                      sum13_r    = _mm_add_ps(xmm13_r[j],sum13_r);
                      xmm13_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+52ull]);
                      sum13_i    = _mm_add_ps(xmm13_i[j],sum13_i);
                      xmm14_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+56ull]);
                      sum14_r    = _mm_add_ps(xmm14_r[j],sum14_r);
                      xmm14_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+56ull]);
                      sum14_i    = _mm_add_ps(xmm14_i[j],sum14_i);
                      xmm15_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+60ull]);
                      sum15_r    = _mm_add_ps(xmm15_r[j],sum15_r);
                      xmm15_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+60ull]);
                      sum15_i    = _mm_add_ps(xmm15_i[j],sum15_i);

                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+4ull],sum1_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+4ull],sum1_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum2_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum2_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+12ull],sum3_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+12ull],sum3_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum4_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum4_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+20ull],sum5_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+20ull],sum5_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum6_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum6_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+28ull],sum7_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+28ull],sum7_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum8_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum8_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+36ull],sum9_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+36ull],sum9_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+40ull],sum10_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+40ull],sum10_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+44ull],sum11_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+44ull],sum11_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum12_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum12_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+52ull],sum13_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+52ull],sum13_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+56ull],sum14_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+56ull],sum14_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+60ull],sum15_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+60ull],sum15_i);
                 }

             }
        }
        

              for(; (i+39ull) < n_samples_halved; i += 40ull) 
              {
                       sum0_r = vzero;
                       sum0_i = vzero;
                       sum1_r = vzero;
                       sum1_i = vzero;
                       sum2_r = vzero;
                       sum2_i = vzero;
                       sum3_r = vzero;
                       sum3_i = vzero;
                       sum4_r = vzero;
                       sum4_i = vzero;
                       sum5_r = vzero;
                       sum5_i = vzero;
                       sum6_r = vzero;
                       sum6_i = vzero;
                       sum7_r = vzero;
                       sum7_i = vzero;
                       sum8_r = vzero;
                       sum8_i = vzero;
                       sum9_r = vzero;
                       sum9_i = vzero;
                       for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                       {
                            xmm0_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                            sum0_r    = _mm_add_ps(xmm0_r[j],sum0_r);
                            xmm0_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                            sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                            xmm1_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+4ull]);
                            sum1_r    = _mm_add_ps(xmm1_r[j],sum1_r);
                            xmm1_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+4ull]);
                            sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                            xmm2_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                            sum2_r    = _mm_add_ps(xmm2_r[j],sum2_r);
                            xmm2_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                            sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                            xmm3_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+12ull]);
                            sum3_r    = _mm_add_ps(xmm3_r[j],sum3_r);
                            xmm3_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+12ull]);
                            sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                            xmm4_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                            sum4_r    = _mm_add_ps(xmm4_r[j],sum4_r);
                            xmm4_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                            sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);
                            xmm5_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+20ull]);
                            sum5_r    = _mm_add_ps(xmm5_r[j],sum5_r);
                            xmm5_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+20ull]);
                            sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);
                            xmm6_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                            sum6_r    = _mm_add_ps(xmm6_r[j],sum6_r);
                            xmm6_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                            sum6_i    = _mm_add_ps(xmm6_i[j],sum6_i);
                            xmm7_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+28ull]);
                            sum7_r    = _mm_add_ps(xmm7_r[j],sum7_r);
                            xmm7_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+28ull]);
                            sum7_i    = _mm_add_ps(xmm7_i[j],sum7_i);
                            xmm8_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                            sum8_r    = _mm_add_ps(xmm8_r[j],sum8_r);
                            xmm8_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                            sum8_i    = _mm_add_ps(xmm8_i[j],sum8_i);
                            xmm9_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+36ull]);
                            sum9_r    = _mm_add_ps(xmm9_r[j],sum9_r);
                            xmm9_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+36ull]);
                            sum9_i    = _mm_add_ps(xmm9_i[j],sum9_i); 

                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+4ull],sum1_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+4ull],sum1_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum2_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum2_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+12ull],sum3_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+12ull],sum3_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum4_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum4_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+20ull],sum5_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+20ull],sum5_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum6_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum6_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+28ull],sum7_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+28ull],sum7_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum8_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum8_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+36ull],sum9_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+36ull],sum9_i);
                       }
               }

               for(; (i+31ull) < n_samples_halved; i += 32ull) 
               {
                       sum0_r = vzero;
                       sum0_i = vzero;
                       sum1_r = vzero;
                       sum1_i = vzero;
                       sum2_r = vzero;
                       sum2_i = vzero;
                       sum3_r = vzero;
                       sum3_i = vzero;
                       sum4_r = vzero;
                       sum4_i = vzero;
                       sum5_r = vzero;
                       sum5_i = vzero;
                       sum6_r = vzero;
                       sum6_i = vzero;
                       sum7_r = vzero;
                       sum7_i = vzero;
                       for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                       {
                            xmm0_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                            sum0_r    = _mm_add_ps(xmm0_r[j],sum0_r);
                            xmm0_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                            sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                            xmm1_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+4ull]);
                            sum1_r    = _mm_add_ps(xmm1_r[j],sum1_r);
                            xmm1_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+4ull]);
                            sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                            xmm2_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                            sum2_r    = _mm_add_ps(xmm2_r[j],sum2_r);
                            xmm2_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                            sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                            xmm3_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+12ull]);
                            sum3_r    = _mm_add_ps(xmm3_r[j],sum3_r);
                            xmm3_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+12ull]);
                            sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                            xmm4_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                            sum4_r    = _mm_add_ps(xmm4_r[j],sum4_r);
                            xmm4_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                            sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);
                            xmm5_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+20ull]);
                            sum5_r    = _mm_add_ps(xmm5_r[j],sum5_r);
                            xmm5_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+20ull]);
                            sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);
                            xmm6_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                            sum6_r    = _mm_add_ps(xmm6_r[j],sum6_r);
                            xmm6_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                            sum6_i    = _mm_add_ps(xmm6_i[j],sum6_i);
                            xmm7_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+28ull]);
                            sum7_r    = _mm_add_ps(xmm7_r[j],sum7_r);
                            xmm7_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+28ull]);
                            sum7_i    = _mm_add_ps(xmm7_i[j],sum7_i);

                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+4ull],sum1_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+4ull],sum1_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum2_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum2_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+12ull],sum3_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+12ull],sum3_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum4_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum4_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+20ull],sum5_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+20ull],sum5_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum6_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum6_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+28ull],sum7_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+28ull],sum7_i);
                      }
              
               }

               for(; (i+15ull) < n_samples_halved; i += 16ull) 
               {
                       sum0_r = vzero;
                       sum0_i = vzero;
                       sum1_r = vzero;
                       sum1_i = vzero;
                       sum2_r = vzero;
                       sum2_i = vzero;
                       sum3_r = vzero;
                       sum3_i = vzero;
                       for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                       {
                            xmm0_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                            sum0_r    = _mm_add_ps(xmm0_r[j],sum0_r);
                            xmm0_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                            sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                            xmm1_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+4ull]);
                            sum1_r    = _mm_add_ps(xmm1_r[j],sum1_r);
                            xmm1_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+4ull]);
                            sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                            xmm2_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                            sum2_r    = _mm_add_ps(xmm2_r[j],sum2_r);
                            xmm2_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                            sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                            xmm3_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+12ull]);
                            sum3_r    = _mm_add_ps(xmm3_r[j],sum3_r);
                            xmm3_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+12ull]);
                            sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                            
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+4ull],sum1_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+4ull],sum1_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum2_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum2_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+12ull],sum3_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+12ull],sum3_i);
                           
                      }
              
               }

               for(; (i+7ull) < n_samples_halved; i += 8ull) 
               {
                       sum0_r = vzero;
                       sum0_i = vzero;
                       sum1_r = vzero;
                       sum1_i = vzero;
                       for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                       {
                            xmm0_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                            sum0_r    = _mm_add_ps(xmm0_r[j],sum0_r);
                            xmm0_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                            sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                            xmm1_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+4ull]);
                            sum1_r    = _mm_add_ps(xmm1_r[j],sum1_r);
                            xmm1_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+4ull]);
                            sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                                                        
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+4ull],sum1_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+4ull],sum1_i);
                                                      
                      }
              
               }

              for(; (i+3ull) < n_samples_halved; i += 4ull) 
              {
                       sum0_r = vzero;
                       sum0_i = vzero;
                       for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                       {
                            xmm0_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                            sum0_r    = _mm_add_ps(xmm0_r[j],sum0_r);
                            xmm0_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                            sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                                                                                   
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                                                                                  
                      }
              }

              for(; (i+0ull) < n_samples_halved; i += 1ull) 
              {
                   sum_rem_r = 0.0f;
                   sum_rem_i = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_nsignals/2; ++j)  
                   {
                       scal_rem_r[j] = this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i];
                       sum_rem_r += scal_rem_r[j];
                       scal_rem_i[j] = this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i];
                       sum_rem_i += scal_rem_i[j];

                       this->m_re_sequence.m_data[i] = sum_rem_r;
                       this->m_im_sequence.m_data[i] = sum_rem_i;
                   }
              }
#if (AM_BB_CMPLX_SINC_SEQUENCE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif
              return (0);

}

std::int32_t 
gms::radiolocation
::am_bb_cmplx_sinc_signal_sequence_t
::signal_sequence_sse42_u10x(const float * __restrict__           sym_in ,
                             am_bb_cmplx_sinc_signal_pdf_params_t &pdf_params,
                             am_bb_cmplx_sinc_signal_rand_distr    rd_enum,
                             const float          * __restrict__  noise_scales,
                             std::int32_t * __restrict__          ret_code,
                             const std::uint32_t                  which_generator)
{
        if(__builtin_expect(nullptr==sym_in,0)) {return (-1);}
        if(__builtin_expect(this->m_nsignals>256ull,0)) {return (-2);}
        if(__builtin_expect((this->m_nsamples%2ull)!=0ull,0)) { return (-3);}     

#if (AM_BB_CMPLX_SINC_SEQUENCE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif
        constexpr std::size_t max_buff_size{256ull};
        __ATTR_ALIGN__(16)
        __m128 xmm0_r[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm1_r[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm2_r[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm3_r[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm4_r[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm5_r[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm6_r[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm7_r[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm8_r[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm9_r[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm0_i[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm1_i[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm2_i[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm3_i[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm4_i[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm5_i[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm6_i[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm7_i[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm8_i[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm9_i[max_buff_size];
         float scal_rem_r[max_buff_size];
         float scal_rem_i[max_buff_size];
        __m128 vzero{_mm_setzero_ps()};
        __m128 sum0_r{vzero};
        __m128 sum1_r{sum0_r};
        __m128 sum2_r{sum0_r};
        __m128 sum3_r{sum0_r};
        __m128 sum4_r{sum0_r};
        __m128 sum5_r{sum0_r};
        __m128 sum6_r{sum0_r};
        __m128 sum7_r{sum0_r};
        __m128 sum8_r{sum0_r};
        __m128 sum9_r{sum0_r};
        __m128 sum0_i{sum0_r};
        __m128 sum1_i{sum0_r};
        __m128 sum2_i{sum0_r};
        __m128 sum3_i{sum0_r};
        __m128 sum4_i{sum0_r};
        __m128 sum5_i{sum0_r};
        __m128 sum6_i{sum0_r};
        __m128 sum7_i{sum0_r};
        __m128 sum8_i{sum0_r};
        __m128 sum9_i{sum0_r};
        std::size_t i;
        float sum_rem_r{0.0f};
        float sum_rem_i{0.0f};
        constexpr std::size_t exceeds_L1D_size{32768u};
        const     std::size_t n_samples_halved{this->m_nsamples/2ull};
        switch(which_generator) 
        {
             case 0 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_signal_user_data(sym_in,
                                                                              static_cast<std::uint32_t>(this->m_nsamples),
                                                                              this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 1 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_signal_user_data_u4x(sym_in,
                                                                              static_cast<std::uint32_t>(this->m_nsamples),
                                                                              this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 2 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {    
                        const float scale{noise_scales[k]};
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_noisy_signal_user_data(pdf_params,
                                                                                      scale,sym_in,
                                                                                      static_cast<std::uint32_t>(this->m_nsamples),
                                                                                      this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 3 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {    
                        const float scale{noise_scales[k]};
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_noisy_signal_user_data(pdf_params,
                                                                                      rd_enum,scale,sym_in,
                                                                                      static_cast<std::uint32_t>(this->m_nsamples),
                                                                                      this->m_nK_values.operator[](k));
                   } 
             }
             break;
             default : 
                   return (-3);
        }

        if(this->m_nsamples > exceeds_L1D_size)
        {
            for(i = 0ull; (i+39ull) < n_samples_halved; i += 40ull) 
            {
              sum0_r = vzero;
              sum0_i = vzero;
              sum1_r = vzero;
              sum1_i = vzero;
              sum2_r = vzero;
              sum2_i = vzero;
              sum3_r = vzero;
              sum3_i = vzero;
              sum4_r = vzero;
              sum4_i = vzero;
              sum5_r = vzero;
              sum5_i = vzero;
              sum6_r = vzero;
              sum6_i = vzero;
              sum7_r = vzero;
              sum7_i = vzero;
              sum8_r = vzero;
              sum8_i = vzero;
              sum9_r = vzero;
              sum9_i = vzero;
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
#if (AM_BB_CMPLX_SINC_SEQUENCE_SOFT_PREFETCH) == 1    
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull],_MM_HINT_T1);
                      
                     
#endif 
                      xmm0_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm_add_ps(xmm0_r[j],sum0_r);
                      xmm0_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                      xmm1_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+4ull]);
                      sum1_r    = _mm_add_ps(xmm1_r[j],sum1_r);
                      xmm1_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+4ull]);
                      sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                      xmm2_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum2_r    = _mm_add_ps(xmm2_r[j],sum2_r);
                      xmm2_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                      xmm3_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+12ull]);
                      sum3_r    = _mm_add_ps(xmm3_r[j],sum3_r);
                      xmm3_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+12ull]);
                      sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                      xmm4_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum4_r    = _mm_add_ps(xmm4_r[j],sum4_r);
                      xmm4_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);
                      xmm5_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+20ull]);
                      sum5_r    = _mm_add_ps(xmm5_r[j],sum5_r);
                      xmm5_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+20ull]);
                      sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);
                      xmm6_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                      sum6_r    = _mm_add_ps(xmm6_r[j],sum6_r);
                      xmm6_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                      sum6_i    = _mm_add_ps(xmm6_i[j],sum6_i);
                      xmm7_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+28ull]);
                      sum7_r    = _mm_add_ps(xmm7_r[j],sum7_r);
                      xmm7_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+28ull]);
                      sum7_i    = _mm_add_ps(xmm7_i[j],sum7_i);
                      xmm8_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum8_r    = _mm_add_ps(xmm8_r[j],sum8_r);
                      xmm8_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum8_i    = _mm_add_ps(xmm8_i[j],sum8_i);
                      xmm9_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+36ull]);
                      sum9_r    = _mm_add_ps(xmm9_r[j],sum9_r);
                      xmm9_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+36ull]);
                      sum9_i    = _mm_add_ps(xmm9_i[j],sum9_i);
                     
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+4ull],sum1_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+4ull],sum1_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum2_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum2_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+12ull],sum3_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+12ull],sum3_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum4_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum4_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+20ull],sum5_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+20ull],sum5_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum6_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum6_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+28ull],sum7_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+28ull],sum7_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum8_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum8_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+36ull],sum9_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+36ull],sum9_i);
                  }

             }

        }
        else 
        {
             for(i = 0ull; (i+39ull) < n_samples_halved; i += 40ull) 
        {
              sum0_r = vzero;
              sum0_i = vzero;
              sum1_r = vzero;
              sum1_i = vzero;
              sum2_r = vzero;
              sum2_i = vzero;
              sum3_r = vzero;
              sum3_i = vzero;
              sum4_r = vzero;
              sum4_i = vzero;
              sum5_r = vzero;
              sum5_i = vzero;
              sum6_r = vzero;
              sum6_i = vzero;
              sum7_r = vzero;
              sum7_i = vzero;
              sum8_r = vzero;
              sum8_i = vzero;
              sum9_r = vzero;
              sum9_i = vzero;
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                      xmm0_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm_add_ps(xmm0_r[j],sum0_r);
                      xmm0_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                      xmm1_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+4ull]);
                      sum1_r    = _mm_add_ps(xmm1_r[j],sum1_r);
                      xmm1_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+4ull]);
                      sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                      xmm2_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum2_r    = _mm_add_ps(xmm2_r[j],sum2_r);
                      xmm2_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                      xmm3_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+12ull]);
                      sum3_r    = _mm_add_ps(xmm3_r[j],sum3_r);
                      xmm3_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+12ull]);
                      sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                      xmm4_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum4_r    = _mm_add_ps(xmm4_r[j],sum4_r);
                      xmm4_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);
                      xmm5_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+20ull]);
                      sum5_r    = _mm_add_ps(xmm5_r[j],sum5_r);
                      xmm5_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+20ull]);
                      sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);
                      xmm6_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                      sum6_r    = _mm_add_ps(xmm6_r[j],sum6_r);
                      xmm6_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                      sum6_i    = _mm_add_ps(xmm6_i[j],sum6_i);
                      xmm7_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+28ull]);
                      sum7_r    = _mm_add_ps(xmm7_r[j],sum7_r);
                      xmm7_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+28ull]);
                      sum7_i    = _mm_add_ps(xmm7_i[j],sum7_i);
                      xmm8_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum8_r    = _mm_add_ps(xmm8_r[j],sum8_r);
                      xmm8_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum8_i    = _mm_add_ps(xmm8_i[j],sum8_i);
                      xmm9_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+36ull]);
                      sum9_r    = _mm_add_ps(xmm9_r[j],sum9_r);
                      xmm9_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+36ull]);
                      sum9_i    = _mm_add_ps(xmm9_i[j],sum9_i);
                     
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+4ull],sum1_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+4ull],sum1_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum2_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum2_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+12ull],sum3_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+12ull],sum3_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum4_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum4_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+20ull],sum5_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+20ull],sum5_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum6_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum6_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+28ull],sum7_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+28ull],sum7_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum8_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum8_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+36ull],sum9_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+36ull],sum9_i);
                  }

             }

        }
                      
               for(; (i+31ull) < n_samples_halved; i += 32ull) 
               {
                       sum0_r = vzero;
                       sum0_i = vzero;
                       sum1_r = vzero;
                       sum1_i = vzero;
                       sum2_r = vzero;
                       sum2_i = vzero;
                       sum3_r = vzero;
                       sum3_i = vzero;
                       sum4_r = vzero;
                       sum4_i = vzero;
                       sum5_r = vzero;
                       sum5_i = vzero;
                       sum6_r = vzero;
                       sum6_i = vzero;
                       sum7_r = vzero;
                       sum7_i = vzero;
                       for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                       {
                            xmm0_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                            sum0_r    = _mm_add_ps(xmm0_r[j],sum0_r);
                            xmm0_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                            sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                            xmm1_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+4ull]);
                            sum1_r    = _mm_add_ps(xmm1_r[j],sum1_r);
                            xmm1_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+4ull]);
                            sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                            xmm2_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                            sum2_r    = _mm_add_ps(xmm2_r[j],sum2_r);
                            xmm2_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                            sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                            xmm3_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+12ull]);
                            sum3_r    = _mm_add_ps(xmm3_r[j],sum3_r);
                            xmm3_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+12ull]);
                            sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                            xmm4_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                            sum4_r    = _mm_add_ps(xmm4_r[j],sum4_r);
                            xmm4_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                            sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);
                            xmm5_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+20ull]);
                            sum5_r    = _mm_add_ps(xmm5_r[j],sum5_r);
                            xmm5_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+20ull]);
                            sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);
                            xmm6_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                            sum6_r    = _mm_add_ps(xmm6_r[j],sum6_r);
                            xmm6_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                            sum6_i    = _mm_add_ps(xmm6_i[j],sum6_i);
                            xmm7_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+28ull]);
                            sum7_r    = _mm_add_ps(xmm7_r[j],sum7_r);
                            xmm7_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+28ull]);
                            sum7_i    = _mm_add_ps(xmm7_i[j],sum7_i);

                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+4ull],sum1_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+4ull],sum1_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum2_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum2_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+12ull],sum3_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+12ull],sum3_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum4_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum4_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+20ull],sum5_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+20ull],sum5_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum6_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum6_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+28ull],sum7_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+28ull],sum7_i);
                      }
              
               }

               for(; (i+15ull) < n_samples_halved; i += 16ull) 
               {
                       sum0_r = vzero;
                       sum0_i = vzero;
                       sum1_r = vzero;
                       sum1_i = vzero;
                       sum2_r = vzero;
                       sum2_i = vzero;
                       sum3_r = vzero;
                       sum3_i = vzero;
                       for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                       {
                            xmm0_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                            sum0_r    = _mm_add_ps(xmm0_r[j],sum0_r);
                            xmm0_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                            sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                            xmm1_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+4ull]);
                            sum1_r    = _mm_add_ps(xmm1_r[j],sum1_r);
                            xmm1_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+4ull]);
                            sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                            xmm2_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                            sum2_r    = _mm_add_ps(xmm2_r[j],sum2_r);
                            xmm2_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                            sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                            xmm3_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+12ull]);
                            sum3_r    = _mm_add_ps(xmm3_r[j],sum3_r);
                            xmm3_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+12ull]);
                            sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                            
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+4ull],sum1_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+4ull],sum1_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum2_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum2_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+12ull],sum3_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+12ull],sum3_i);
                           
                      }
              
               }

               for(; (i+7ull) < n_samples_halved; i += 8ull) 
               {
                       sum0_r = vzero;
                       sum0_i = vzero;
                       sum1_r = vzero;
                       sum1_i = vzero;
                       for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                       {
                            xmm0_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                            sum0_r    = _mm_add_ps(xmm0_r[j],sum0_r);
                            xmm0_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                            sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                            xmm1_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+4ull]);
                            sum1_r    = _mm_add_ps(xmm1_r[j],sum1_r);
                            xmm1_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+4ull]);
                            sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                                                        
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+4ull],sum1_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+4ull],sum1_i);
                                                      
                      }
              
               }

              for(; (i+3ull) < n_samples_halved; i += 4ull) 
              {
                       sum0_r = vzero;
                       sum0_i = vzero;
                       for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                       {
                            xmm0_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                            sum0_r    = _mm_add_ps(xmm0_r[j],sum0_r);
                            xmm0_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                            sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                                                                                   
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                                                                                  
                      }
              }

              for(; (i+0ull) < n_samples_halved; i += 1ull) 
              {
                   sum_rem_r = 0.0f;
                   sum_rem_i = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_nsignals/2; ++j)  
                   {
                       scal_rem_r[j] = this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i];
                       sum_rem_r += scal_rem_r[j];
                       scal_rem_i[j] = this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i];
                       sum_rem_i += scal_rem_i[j];

                       this->m_re_sequence.m_data[i] = sum_rem_r;
                       this->m_im_sequence.m_data[i] = sum_rem_i;
                   }
              }

#if (AM_BB_CMPLX_SINC_SEQUENCE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif
              return (0);

}

std::int32_t 
gms::radiolocation
::am_bb_cmplx_sinc_signal_sequence_t
::signal_sequence_sse42_u6x(const float * __restrict__           sym_in ,
                             am_bb_cmplx_sinc_signal_pdf_params_t &pdf_params,
                             am_bb_cmplx_sinc_signal_rand_distr    rd_enum,
                             const float          * __restrict__  noise_scales,
                             std::int32_t * __restrict__          ret_code,
                             const std::uint32_t                  which_generator)
{
        if(__builtin_expect(nullptr==sym_in,0)) {return (-1);}
        if(__builtin_expect(this->m_nsignals>256ull,0)) {return (-2);}
        if(__builtin_expect((this->m_nsamples%2ull)!=0,0)) {return (-3);}

#if (AM_BB_CMPLX_SINC_SEQUENCE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif
        constexpr std::size_t max_buff_size{256ull};
        __ATTR_ALIGN__(16)
        __m128 xmm0_r[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm1_r[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm2_r[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm3_r[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm4_r[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm5_r[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm0_i[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm1_i[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm2_i[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm3_i[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm4_i[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm5_i[max_buff_size];
         float scal_rem_r[max_buff_size];
         float scal_rem_i[max_buff_size];
        __m128 vzero{_mm_setzero_ps()};
        __m128 sum0_r{vzero};
        __m128 sum1_r{sum0_r};
        __m128 sum2_r{sum0_r};
        __m128 sum3_r{sum0_r};
        __m128 sum4_r{sum0_r};
        __m128 sum5_r{sum0_r};
        __m128 sum0_i{sum0_r};
        __m128 sum1_i{sum0_r};
        __m128 sum2_i{sum0_r};
        __m128 sum3_i{sum0_r};
        __m128 sum4_i{sum0_r};
        __m128 sum5_i{sum0_r};
        std::size_t i;
        float sum_rem_r{0.0f};
        float sum_rem_i{0.0f};
        constexpr std::size_t exceeds_L1D_size{32768ull};
        const     std::size_t n_samples_halved{this->m_nsamples/2ull};

        switch(which_generator) 
        {
             case 0 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_signal_user_data(sym_in,
                                                                              static_cast<std::uint32_t>(this->m_nsamples),
                                                                              this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 1 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_signal_user_data_u4x(sym_in,
                                                                              static_cast<std::uint32_t>(this->m_nsamples),
                                                                              this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 2 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {    
                        const float scale{noise_scales[k]};
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_noisy_signal_user_data(pdf_params,
                                                                                      scale,sym_in,
                                                                                      static_cast<std::uint32_t>(this->m_nsamples),
                                                                                      this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 3 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {    
                        const float scale{noise_scales[k]};
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_noisy_signal_user_data(pdf_params,
                                                                                      rd_enum,scale,sym_in,
                                                                                      static_cast<std::uint32_t>(this->m_nsamples),
                                                                                      this->m_nK_values.operator[](k));
                   } 
             }
             break;
             default : 
                   return (-3);
        }

        if(this->m_nsamples > exceeds_L1D_size) 
        {
             for(i = 0ull; (i+23ull) < n_samples_halved; i += 24ull) 
            {
              sum0_r = vzero;
              sum0_i = vzero;
              sum1_r = vzero;
              sum1_i = vzero;
              sum2_r = vzero;
              sum2_i = vzero;
              sum3_r = vzero;
              sum3_i = vzero;
              sum4_r = vzero;
              sum4_i = vzero;
              sum5_r = vzero;
              sum5_i = vzero;
             
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
#if (AM_BB_CMPLX_SINC_SEQUENCE_SOFT_PREFETCH) == 1    
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull],_MM_HINT_T1);
                               
                     
#endif 
                      xmm0_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm_add_ps(xmm0_r[j],sum0_r);
                      xmm0_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                      xmm1_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+4ull]);
                      sum1_r    = _mm_add_ps(xmm1_r[j],sum1_r);
                      xmm1_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+4ull]);
                      sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                      xmm2_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum2_r    = _mm_add_ps(xmm2_r[j],sum2_r);
                      xmm2_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                      xmm3_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+12ull]);
                      sum3_r    = _mm_add_ps(xmm3_r[j],sum3_r);
                      xmm3_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+12ull]);
                      sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                      xmm4_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum4_r    = _mm_add_ps(xmm4_r[j],sum4_r);
                      xmm4_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);
                      xmm5_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+20ull]);
                      sum5_r    = _mm_add_ps(xmm5_r[j],sum5_r);
                      xmm5_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+20ull]);
                      sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);
                     
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+4ull],sum1_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+4ull],sum1_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum2_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum2_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+12ull],sum3_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+12ull],sum3_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum4_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum4_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+20ull],sum5_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+20ull],sum5_i);
                      
                  }

             }
        }
        else 
        {
            for(i = 0ull; (i+23ull) < n_samples_halved; i += 24ull) 
        {
              sum0_r = vzero;
              sum0_i = vzero;
              sum1_r = vzero;
              sum1_i = vzero;
              sum2_r = vzero;
              sum2_i = vzero;
              sum3_r = vzero;
              sum3_i = vzero;
              sum4_r = vzero;
              sum4_i = vzero;
              sum5_r = vzero;
              sum5_i = vzero;
             
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                      xmm0_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm_add_ps(xmm0_r[j],sum0_r);
                      xmm0_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                      xmm1_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+4ull]);
                      sum1_r    = _mm_add_ps(xmm1_r[j],sum1_r);
                      xmm1_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+4ull]);
                      sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                      xmm2_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum2_r    = _mm_add_ps(xmm2_r[j],sum2_r);
                      xmm2_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                      xmm3_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+12ull]);
                      sum3_r    = _mm_add_ps(xmm3_r[j],sum3_r);
                      xmm3_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+12ull]);
                      sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                      xmm4_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum4_r    = _mm_add_ps(xmm4_r[j],sum4_r);
                      xmm4_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);
                      xmm5_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+20ull]);
                      sum5_r    = _mm_add_ps(xmm5_r[j],sum5_r);
                      xmm5_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+20ull]);
                      sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);
                     
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+4ull],sum1_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+4ull],sum1_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum2_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum2_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+12ull],sum3_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+12ull],sum3_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum4_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum4_i);
                      _mm_storeu_ps(&this->m_re_sequence.m_data[i+20ull],sum5_r);
                      _mm_storeu_ps(&this->m_im_sequence.m_data[i+20ull],sum5_i);
                      
                  }

             }
        }
                     
            for(; (i+15ull) < n_samples_halved; i += 16ull) 
            {
                       sum0_r = vzero;
                       sum0_i = vzero;
                       sum1_r = vzero;
                       sum1_i = vzero;
                       sum2_r = vzero;
                       sum2_i = vzero;
                       sum3_r = vzero;
                       sum3_i = vzero;
                       for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                       {
                            xmm0_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                            sum0_r    = _mm_add_ps(xmm0_r[j],sum0_r);
                            xmm0_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                            sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                            xmm1_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+4ull]);
                            sum1_r    = _mm_add_ps(xmm1_r[j],sum1_r);
                            xmm1_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+4ull]);
                            sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                            xmm2_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                            sum2_r    = _mm_add_ps(xmm2_r[j],sum2_r);
                            xmm2_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                            sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                            xmm3_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+12ull]);
                            sum3_r    = _mm_add_ps(xmm3_r[j],sum3_r);
                            xmm3_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+12ull]);
                            sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                            
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+4ull],sum1_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+4ull],sum1_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum2_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum2_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+12ull],sum3_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+12ull],sum3_i);
                           
                      }
              
               }

               for(; (i+7ull) < n_samples_halved; i += 8ull) 
               {
                       sum0_r = vzero;
                       sum0_i = vzero;
                       sum1_r = vzero;
                       sum1_i = vzero;
                       for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                       {
                            xmm0_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                            sum0_r    = _mm_add_ps(xmm0_r[j],sum0_r);
                            xmm0_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                            sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                            xmm1_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+4ull]);
                            sum1_r    = _mm_add_ps(xmm1_r[j],sum1_r);
                            xmm1_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+4ull]);
                            sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                                                        
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+4ull],sum1_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+4ull],sum1_i);
                                                      
                      }
              
               }

              for(; (i+3ull) < n_samples_halved; i += 4ull) 
              {
                       sum0_r = vzero;
                       sum0_i = vzero;
                       for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                       {
                            xmm0_r[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                            sum0_r    = _mm_add_ps(xmm0_r[j],sum0_r);
                            xmm0_i[j] = _mm_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                            sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                                                                                   
                            _mm_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                            _mm_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                                                                                  
                      }
              }

              for(; (i+0ull) < n_samples_halved; i += 1ull) 
              {
                   sum_rem_r = 0.0f;
                   sum_rem_i = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_nsignals/2; ++j)  
                   {
                       scal_rem_r[j] = this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i];
                       sum_rem_r += scal_rem_r[j];
                       scal_rem_i[j] = this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i];
                       sum_rem_i += scal_rem_i[j];

                       this->m_re_sequence.m_data[i] = sum_rem_r;
                       this->m_im_sequence.m_data[i] = sum_rem_i;
                   }
              }

#if (AM_BB_CMPLX_SINC_SEQUENCE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif
              return (0);

}

std::int32_t 
gms::radiolocation
::am_bb_cmplx_sinc_signal_sequence_t
::signal_sequence_avx_u16x(  const float * __restrict__           sym_in ,
                             am_bb_cmplx_sinc_signal_pdf_params_t &pdf_params,
                             am_bb_cmplx_sinc_signal_rand_distr    rd_enum,
                             const float          * __restrict__  noise_scales,
                             std::int32_t * __restrict__          ret_code,
                             const std::uint32_t                  which_generator)
{
        if(__builtin_expect(nullptr==sym_in,0)) {return (-1);}
        if(__builtin_expect(this->m_nsignals>256ull,0)) {return (-2);}
        if(__builtin_expect((this->m_nsamples%2ull)!=0ull,0)) { return (-3);}

#if (AM_BB_CMPLX_SINC_SEQUENCE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif
        constexpr std::size_t max_buff_size{256ull};
        __ATTR_ALIGN__(32)
        __m256 ymm0_r[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm1_r[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm2_r[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm3_r[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm4_r[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm5_r[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm6_r[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm7_r[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm8_r[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm9_r[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm10_r[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm11_r[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm12_r[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm13_r[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm14_r[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm15_r[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm0_i[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm1_i[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm2_i[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm3_i[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm4_i[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm5_i[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm6_i[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm7_i[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm8_i[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm9_i[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm10_i[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm11_i[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm12_i[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm13_i[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm14_i[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm15_i[max_buff_size];
         float scal_rem_r[max_buff_size];
         float scal_rem_i[max_buff_size];
        __m256 vzero{_mm256_setzero_ps()};
        __m256 sum0_r{vzero};
        __m256 sum1_r{sum0_r};
        __m256 sum2_r{sum0_r};
        __m256 sum3_r{sum0_r};
        __m256 sum4_r{sum0_r};
        __m256 sum5_r{sum0_r};
        __m256 sum6_r{sum0_r};
        __m256 sum7_r{sum0_r};
        __m256 sum8_r{sum0_r};
        __m256 sum9_r{sum0_r};
        __m256 sum10_r{sum0_r};
        __m256 sum11_r{sum0_r};
        __m256 sum12_r{sum0_r};
        __m256 sum13_r{sum0_r};
        __m256 sum14_r{sum0_r};
        __m256 sum15_r{sum0_r};
        __m256 sum0_i{sum0_r};
        __m256 sum1_i{sum0_r};
        __m256 sum2_i{sum0_r};
        __m256 sum3_i{sum0_r};
        __m256 sum4_i{sum0_r};
        __m256 sum5_i{sum0_r};
        __m256 sum6_i{sum0_r};
        __m256 sum7_i{sum0_r};
        __m256 sum8_i{sum0_r};
        __m256 sum9_i{sum0_r};
        __m256 sum10_i{sum0_r};
        __m256 sum11_i{sum0_r};
        __m256 sum12_i{sum0_r};
        __m256 sum13_i{sum0_r};
        __m256 sum14_i{sum0_r};
        __m256 sum15_i{sum0_r}; 
        std::size_t i;
        float sum_rem_r{0.0f};
        float sum_rem_i{0.0f};
        constexpr std::size_t exceeds_L1D_size{32768u};
        const     std::size_t n_samples_halved{this->m_nsamples/2ull};
        switch(which_generator) 
        {
             case 0 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_signal_user_data(sym_in,
                                                                              static_cast<std::uint32_t>(this->m_nsamples),
                                                                              this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 1 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_signal_user_data_u4x(sym_in,
                                                                              static_cast<std::uint32_t>(this->m_nsamples),
                                                                              this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 2 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {    
                        const float scale{noise_scales[k]};
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_noisy_signal_user_data(pdf_params,
                                                                                      scale,sym_in,
                                                                                      static_cast<std::uint32_t>(this->m_nsamples),
                                                                                      this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 3 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {    
                        const float scale{noise_scales[k]};
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_noisy_signal_user_data(pdf_params,
                                                                                      rd_enum,scale,sym_in,
                                                                                      static_cast<std::uint32_t>(this->m_nsamples),
                                                                                      this->m_nK_values.operator[](k));
                   } 
             }
             break;
             default : 
                   return (-3);
        }

        if(this->m_nsamples > exceeds_L1D_size)
        {

            for(i = 0ull; (i+127ull) < n_samples_halved; i += 128ull) 
            {
              sum0_r = vzero;
              sum0_i = vzero;
              sum1_r = vzero;
              sum1_i = vzero;
              sum2_r = vzero;
              sum2_i = vzero;
              sum3_r = vzero;
              sum3_i = vzero;
              sum4_r = vzero;
              sum4_i = vzero;
              sum5_r = vzero;
              sum5_i = vzero;
              sum6_r = vzero;
              sum6_i = vzero;
              sum7_r = vzero;
              sum7_i = vzero;
              sum8_r = vzero;
              sum8_i = vzero;
              sum9_r = vzero;
              sum9_i = vzero;
              sum10_r = vzero;
              sum10_i = vzero;
              sum11_r = vzero;
              sum11_i = vzero;
              sum12_r = vzero;
              sum12_i = vzero;
              sum13_r = vzero;
              sum13_i = vzero;
              sum14_r = vzero;
              sum14_i = vzero;
              sum15_r = vzero;
              sum15_i = vzero;
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
#if (AM_BB_CMPLX_SINC_SEQUENCE_SOFT_PREFETCH) == 1    
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+80ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+80ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+96ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+96ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+112ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+112ull],_MM_HINT_T1);
#endif                       
                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                      ymm1_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum1_r    = _mm256_add_ps(ymm1_r[j],sum1_r);
                      ymm1_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum1_i    = _mm256_add_ps(ymm1_i[j],sum1_i);
                      ymm2_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum2_r    = _mm256_add_ps(ymm2_r[j],sum2_r);
                      ymm2_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum2_i    = _mm256_add_ps(ymm2_i[j],sum2_i);
                      ymm3_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                      sum3_r    = _mm256_add_ps(ymm3_r[j],sum3_r);
                      ymm3_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                      sum3_i    = _mm256_add_ps(ymm3_i[j],sum3_i);
                      ymm4_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum4_r    = _mm256_add_ps(ymm4_r[j],sum4_r);
                      ymm4_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum4_i    = _mm256_add_ps(ymm4_i[j],sum4_i);
                      ymm5_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+40ull]);
                      sum5_r    = _mm256_add_ps(ymm5_r[j],sum5_r);
                      ymm5_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+40ull]);
                      sum5_i    = _mm256_add_ps(ymm5_i[j],sum5_i);
                      ymm6_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum6_r    = _mm256_add_ps(ymm6_r[j],sum6_r);
                      ymm6_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum6_i    = _mm256_add_ps(ymm6_i[j],sum6_i);
                      ymm7_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+56ull]);
                      sum7_r    = _mm256_add_ps(ymm7_r[j],sum7_r);
                      ymm7_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+56ull]);
                      sum7_i    = _mm256_add_ps(ymm7_i[j],sum7_i);
                      ymm8_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull]);
                      sum8_r    = _mm256_add_ps(ymm8_r[j],sum8_r);
                      ymm8_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull]);
                      sum8_i    = _mm256_add_ps(ymm8_i[j],sum8_i);
                      ymm9_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+72ull]);
                      sum9_r    = _mm256_add_ps(ymm9_r[j],sum9_r);
                      ymm9_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+72ull]);
                      sum9_i    = _mm256_add_ps(ymm9_i[j],sum9_i);
                      ymm10_r[j]= _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+80ull]);
                      sum10_r   = _mm256_add_ps(ymm10_r[j],sum10_r);
                      ymm10_i[j]= _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+80ull]);
                      sum10_i   = _mm256_add_ps(ymm10_i[j],sum10_i);
                      ymm11_r[j]= _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+88ull]);
                      sum11_r    = _mm256_add_ps(ymm11_r[j],sum11_r);
                      ymm11_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+88ull]);
                      sum11_i    = _mm256_add_ps(ymm11_i[j],sum11_i);
                      ymm12_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+96ull]);
                      sum12_r    = _mm256_add_ps(ymm12_r[j],sum12_r);
                      ymm12_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+96ull]);
                      sum12_i    = _mm256_add_ps(ymm12_i[j],sum12_i);
                      ymm13_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+104ull]);
                      sum13_r    = _mm256_add_ps(ymm13_r[j],sum13_r);
                      ymm13_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+104ull]);
                      sum13_i    = _mm256_add_ps(ymm13_i[j],sum13_i);
                      ymm14_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+112ull]);
                      sum14_r    = _mm256_add_ps(ymm14_r[j],sum14_r);
                      ymm14_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+112ull]);
                      sum14_i    = _mm256_add_ps(ymm14_i[j],sum14_i);
                      ymm15_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+120ull]);
                      sum15_r    = _mm256_add_ps(ymm15_r[j],sum15_r);
                      ymm15_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+120ull]);
                      sum15_i    = _mm256_add_ps(ymm15_i[j],sum15_i);

                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum1_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum1_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum2_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum2_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum3_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum3_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum4_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum4_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+40ull],sum5_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+40ull],sum5_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum6_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum6_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+56ull],sum7_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+56ull],sum7_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+64ull],sum8_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+64ull],sum8_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+72ull],sum9_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+72ull],sum9_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+80ull],sum10_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+80ull],sum10_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+88ull],sum11_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+88ull],sum11_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+96ull],sum12_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+96ull],sum12_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+104ull],sum13_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+104ull],sum13_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+112ull],sum14_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+112ull],sum14_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+120ull],sum15_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+120ull],sum15_i);
                 }

             }
        }
        else 
        {
            for(i = 0ull; (i+127ull) < n_samples_halved; i += 128ull) 
        {
              sum0_r = vzero;
              sum0_i = vzero;
              sum1_r = vzero;
              sum1_i = vzero;
              sum2_r = vzero;
              sum2_i = vzero;
              sum3_r = vzero;
              sum3_i = vzero;
              sum4_r = vzero;
              sum4_i = vzero;
              sum5_r = vzero;
              sum5_i = vzero;
              sum6_r = vzero;
              sum6_i = vzero;
              sum7_r = vzero;
              sum7_i = vzero;
              sum8_r = vzero;
              sum8_i = vzero;
              sum9_r = vzero;
              sum9_i = vzero;
              sum10_r = vzero;
              sum10_i = vzero;
              sum11_r = vzero;
              sum11_i = vzero;
              sum12_r = vzero;
              sum12_i = vzero;
              sum13_r = vzero;
              sum13_i = vzero;
              sum14_r = vzero;
              sum14_i = vzero;
              sum15_r = vzero;
              sum15_i = vzero;
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                      ymm1_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum1_r    = _mm256_add_ps(ymm1_r[j],sum1_r);
                      ymm1_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum1_i    = _mm256_add_ps(ymm1_i[j],sum1_i);
                      ymm2_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum2_r    = _mm256_add_ps(ymm2_r[j],sum2_r);
                      ymm2_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum2_i    = _mm256_add_ps(ymm2_i[j],sum2_i);
                      ymm3_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                      sum3_r    = _mm256_add_ps(ymm3_r[j],sum3_r);
                      ymm3_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                      sum3_i    = _mm256_add_ps(ymm3_i[j],sum3_i);
                      ymm4_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum4_r    = _mm256_add_ps(ymm4_r[j],sum4_r);
                      ymm4_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum4_i    = _mm256_add_ps(ymm4_i[j],sum4_i);
                      ymm5_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+40ull]);
                      sum5_r    = _mm256_add_ps(ymm5_r[j],sum5_r);
                      ymm5_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+40ull]);
                      sum5_i    = _mm256_add_ps(ymm5_i[j],sum5_i);
                      ymm6_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum6_r    = _mm256_add_ps(ymm6_r[j],sum6_r);
                      ymm6_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum6_i    = _mm256_add_ps(ymm6_i[j],sum6_i);
                      ymm7_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+56ull]);
                      sum7_r    = _mm256_add_ps(ymm7_r[j],sum7_r);
                      ymm7_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+56ull]);
                      sum7_i    = _mm256_add_ps(ymm7_i[j],sum7_i);
                      ymm8_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull]);
                      sum8_r    = _mm256_add_ps(ymm8_r[j],sum8_r);
                      ymm8_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull]);
                      sum8_i    = _mm256_add_ps(ymm8_i[j],sum8_i);
                      ymm9_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+72ull]);
                      sum9_r    = _mm256_add_ps(ymm9_r[j],sum9_r);
                      ymm9_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+72ull]);
                      sum9_i    = _mm256_add_ps(ymm9_i[j],sum9_i);
                      ymm10_r[j]= _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+80ull]);
                      sum10_r   = _mm256_add_ps(ymm10_r[j],sum10_r);
                      ymm10_i[j]= _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+80ull]);
                      sum10_i   = _mm256_add_ps(ymm10_i[j],sum10_i);
                      ymm11_r[j]= _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+88ull]);
                      sum11_r    = _mm256_add_ps(ymm11_r[j],sum11_r);
                      ymm11_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+88ull]);
                      sum11_i    = _mm256_add_ps(ymm11_i[j],sum11_i);
                      ymm12_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+96ull]);
                      sum12_r    = _mm256_add_ps(ymm12_r[j],sum12_r);
                      ymm12_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+96ull]);
                      sum12_i    = _mm256_add_ps(ymm12_i[j],sum12_i);
                      ymm13_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+104ull]);
                      sum13_r    = _mm256_add_ps(ymm13_r[j],sum13_r);
                      ymm13_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+104ull]);
                      sum13_i    = _mm256_add_ps(ymm13_i[j],sum13_i);
                      ymm14_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+112ull]);
                      sum14_r    = _mm256_add_ps(ymm14_r[j],sum14_r);
                      ymm14_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+112ull]);
                      sum14_i    = _mm256_add_ps(ymm14_i[j],sum14_i);
                      ymm15_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+120ull]);
                      sum15_r    = _mm256_add_ps(ymm15_r[j],sum15_r);
                      ymm15_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+120ull]);
                      sum15_i    = _mm256_add_ps(ymm15_i[j],sum15_i);

                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum1_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum1_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum2_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum2_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum3_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum3_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum4_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum4_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+40ull],sum5_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+40ull],sum5_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum6_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum6_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+56ull],sum7_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+56ull],sum7_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+64ull],sum8_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+64ull],sum8_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+72ull],sum9_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+72ull],sum9_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+80ull],sum10_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+80ull],sum10_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+88ull],sum11_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+88ull],sum11_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+96ull],sum12_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+96ull],sum12_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+104ull],sum13_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+104ull],sum13_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+112ull],sum14_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+112ull],sum14_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+120ull],sum15_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+120ull],sum15_i);
                 }

             }
        }
        
              for(; (i+95ull) < n_samples_halved; i += 96ull) 
              {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                        sum2_r = vzero;
                        sum2_i = vzero;
                        sum3_r = vzero;
                        sum3_i = vzero;
                        sum4_r = vzero;
                        sum4_i = vzero;
                        sum5_r = vzero;
                        sum5_i = vzero;
                        sum6_r = vzero;
                        sum6_i = vzero;
                        sum7_r = vzero;
                        sum7_i = vzero;
                        sum8_r = vzero;
                        sum8_i = vzero;
                        sum9_r = vzero;
                        sum9_i = vzero;
                        sum10_r = vzero;
                        sum10_i = vzero;
                        sum11_r = vzero;
                        sum11_i = vzero;
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                      ymm1_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum1_r    = _mm256_add_ps(ymm1_r[j],sum1_r);
                      ymm1_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum1_i    = _mm256_add_ps(ymm1_i[j],sum1_i);
                      ymm2_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum2_r    = _mm256_add_ps(ymm2_r[j],sum2_r);
                      ymm2_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum2_i    = _mm256_add_ps(ymm2_i[j],sum2_i);
                      ymm3_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                      sum3_r    = _mm256_add_ps(ymm3_r[j],sum3_r);
                      ymm3_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                      sum3_i    = _mm256_add_ps(ymm3_i[j],sum3_i);
                      ymm4_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum4_r    = _mm256_add_ps(ymm4_r[j],sum4_r);
                      ymm4_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum4_i    = _mm256_add_ps(ymm4_i[j],sum4_i);
                      ymm5_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+40ull]);
                      sum5_r    = _mm256_add_ps(ymm5_r[j],sum5_r);
                      ymm5_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+40ull]);
                      sum5_i    = _mm256_add_ps(ymm5_i[j],sum5_i);
                      ymm6_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum6_r    = _mm256_add_ps(ymm6_r[j],sum6_r);
                      ymm6_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum6_i    = _mm256_add_ps(ymm6_i[j],sum6_i);
                      ymm7_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+56ull]);
                      sum7_r    = _mm256_add_ps(ymm7_r[j],sum7_r);
                      ymm7_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+56ull]);
                      sum7_i    = _mm256_add_ps(ymm7_i[j],sum7_i);
                      ymm8_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull]);
                      sum8_r    = _mm256_add_ps(ymm8_r[j],sum8_r);
                      ymm8_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull]);
                      sum8_i    = _mm256_add_ps(ymm8_i[j],sum8_i);
                      ymm9_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+72ull]);
                      sum9_r    = _mm256_add_ps(ymm9_r[j],sum9_r);
                      ymm9_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+72ull]);
                      sum9_i    = _mm256_add_ps(ymm9_i[j],sum9_i);
                      ymm10_r[j]= _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+80ull]);
                      sum10_r   = _mm256_add_ps(ymm10_r[j],sum10_r);
                      ymm10_i[j]= _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+80ull]);
                      sum10_i   = _mm256_add_ps(ymm10_i[j],sum10_i);
                      ymm11_r[j]= _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+88ull]);
                      sum11_r    = _mm256_add_ps(ymm11_r[j],sum11_r);
                      ymm11_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+88ull]);
                      sum11_i    = _mm256_add_ps(ymm11_i[j],sum11_i);

                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum1_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum1_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum2_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum2_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum3_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum3_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum4_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum4_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+40ull],sum5_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+40ull],sum5_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum6_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum6_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+56ull],sum7_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+56ull],sum7_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+64ull],sum8_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+64ull],sum8_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+72ull],sum9_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+72ull],sum9_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+80ull],sum10_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+80ull],sum10_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+88ull],sum11_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+88ull],sum11_i);
                  }
            }

            for(; (i+79ull) < n_samples_halved; i += 80ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                        sum2_r = vzero;
                        sum2_i = vzero;
                        sum3_r = vzero;
                        sum3_i = vzero;
                        sum4_r = vzero;
                        sum4_i = vzero;
                        sum5_r = vzero;
                        sum5_i = vzero;
                        sum6_r = vzero;
                        sum6_i = vzero;
                        sum7_r = vzero;
                        sum7_i = vzero;
                        sum8_r = vzero;
                        sum8_i = vzero;
                        sum9_r = vzero;
                        sum9_i = vzero;
                       
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                      ymm1_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum1_r    = _mm256_add_ps(ymm1_r[j],sum1_r);
                      ymm1_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum1_i    = _mm256_add_ps(ymm1_i[j],sum1_i);
                      ymm2_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum2_r    = _mm256_add_ps(ymm2_r[j],sum2_r);
                      ymm2_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum2_i    = _mm256_add_ps(ymm2_i[j],sum2_i);
                      ymm3_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                      sum3_r    = _mm256_add_ps(ymm3_r[j],sum3_r);
                      ymm3_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                      sum3_i    = _mm256_add_ps(ymm3_i[j],sum3_i);
                      ymm4_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum4_r    = _mm256_add_ps(ymm4_r[j],sum4_r);
                      ymm4_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum4_i    = _mm256_add_ps(ymm4_i[j],sum4_i);
                      ymm5_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+40ull]);
                      sum5_r    = _mm256_add_ps(ymm5_r[j],sum5_r);
                      ymm5_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+40ull]);
                      sum5_i    = _mm256_add_ps(ymm5_i[j],sum5_i);
                      ymm6_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum6_r    = _mm256_add_ps(ymm6_r[j],sum6_r);
                      ymm6_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum6_i    = _mm256_add_ps(ymm6_i[j],sum6_i);
                      ymm7_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+56ull]);
                      sum7_r    = _mm256_add_ps(ymm7_r[j],sum7_r);
                      ymm7_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+56ull]);
                      sum7_i    = _mm256_add_ps(ymm7_i[j],sum7_i);
                      ymm8_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull]);
                      sum8_r    = _mm256_add_ps(ymm8_r[j],sum8_r);
                      ymm8_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull]);
                      sum8_i    = _mm256_add_ps(ymm8_i[j],sum8_i);
                      ymm9_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+72ull]);
                      sum9_r    = _mm256_add_ps(ymm9_r[j],sum9_r);
                      ymm9_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+72ull]);
                      sum9_i    = _mm256_add_ps(ymm9_i[j],sum9_i);
                     

                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum1_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum1_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum2_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum2_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum3_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum3_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum4_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum4_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+40ull],sum5_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+40ull],sum5_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum6_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum6_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+56ull],sum7_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+56ull],sum7_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+64ull],sum8_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+64ull],sum8_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+72ull],sum9_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+72ull],sum9_i);
                      
                  }
            }

            
            for(; (i+63ull) < n_samples_halved; i += 64ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                        sum2_r = vzero;
                        sum2_i = vzero;
                        sum3_r = vzero;
                        sum3_i = vzero;
                        sum4_r = vzero;
                        sum4_i = vzero;
                        sum5_r = vzero;
                        sum5_i = vzero;
                        sum6_r = vzero;
                        sum6_i = vzero;
                        sum7_r = vzero;
                        sum7_i = vzero;
                        
                       
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                      ymm1_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum1_r    = _mm256_add_ps(ymm1_r[j],sum1_r);
                      ymm1_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum1_i    = _mm256_add_ps(ymm1_i[j],sum1_i);
                      ymm2_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum2_r    = _mm256_add_ps(ymm2_r[j],sum2_r);
                      ymm2_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum2_i    = _mm256_add_ps(ymm2_i[j],sum2_i);
                      ymm3_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                      sum3_r    = _mm256_add_ps(ymm3_r[j],sum3_r);
                      ymm3_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                      sum3_i    = _mm256_add_ps(ymm3_i[j],sum3_i);
                      ymm4_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum4_r    = _mm256_add_ps(ymm4_r[j],sum4_r);
                      ymm4_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum4_i    = _mm256_add_ps(ymm4_i[j],sum4_i);
                      ymm5_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+40ull]);
                      sum5_r    = _mm256_add_ps(ymm5_r[j],sum5_r);
                      ymm5_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+40ull]);
                      sum5_i    = _mm256_add_ps(ymm5_i[j],sum5_i);
                      ymm6_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum6_r    = _mm256_add_ps(ymm6_r[j],sum6_r);
                      ymm6_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum6_i    = _mm256_add_ps(ymm6_i[j],sum6_i);
                      ymm7_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+56ull]);
                      sum7_r    = _mm256_add_ps(ymm7_r[j],sum7_r);
                      ymm7_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+56ull]);
                      sum7_i    = _mm256_add_ps(ymm7_i[j],sum7_i);
                                           

                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum1_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum1_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum2_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum2_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum3_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum3_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum4_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum4_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+40ull],sum5_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+40ull],sum5_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum6_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum6_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+56ull],sum7_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+56ull],sum7_i);
                                            
                  }
            }

            
            for(; (i+31ull) < n_samples_halved; i += 32ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                        sum2_r = vzero;
                        sum2_i = vzero;
                        sum3_r = vzero;
                        sum3_i = vzero;
                                            
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                      ymm1_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum1_r    = _mm256_add_ps(ymm1_r[j],sum1_r);
                      ymm1_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum1_i    = _mm256_add_ps(ymm1_i[j],sum1_i);
                      ymm2_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum2_r    = _mm256_add_ps(ymm2_r[j],sum2_r);
                      ymm2_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum2_i    = _mm256_add_ps(ymm2_i[j],sum2_i);
                      ymm3_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                      sum3_r    = _mm256_add_ps(ymm3_r[j],sum3_r);
                      ymm3_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                      sum3_i    = _mm256_add_ps(ymm3_i[j],sum3_i);
                                                              
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum1_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum1_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum2_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum2_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum3_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum3_i);
                                                                 
                  }
            }

               
            for(; (i+15ull) < n_samples_halved; i += 16ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                                                                    
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                      ymm1_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum1_r    = _mm256_add_ps(ymm1_r[j],sum1_r);
                      ymm1_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum1_i    = _mm256_add_ps(ymm1_i[j],sum1_i);
                                                                                    
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum1_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum1_i);
                                                                                       
                  }
            }


            for(; (i+7ull) < n_samples_halved; i += 8ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                                                                                          
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                                                                                                         
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                                                                                                            
                  }
            } 

              for(; (i+0ull) < n_samples_halved; i += 1ull) 
              {
                   sum_rem_r = 0.0f;
                   sum_rem_i = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                   {
                       scal_rem_r[j] = this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i];
                       sum_rem_r += scal_rem_r[j];
                       scal_rem_i[j] = this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i];
                       sum_rem_i += scal_rem_i[j];

                       this->m_re_sequence.m_data[i] = sum_rem_r;
                       this->m_im_sequence.m_data[i] = sum_rem_i;
                   }
              }

#if (AM_BB_CMPLX_SINC_SEQUENCE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif
              return (0);

}

std::int32_t 
gms::radiolocation
::am_bb_cmplx_sinc_signal_sequence_t
::signal_sequence_avx_u10x(  const float * __restrict__           sym_in ,
                             am_bb_cmplx_sinc_signal_pdf_params_t &pdf_params,
                             am_bb_cmplx_sinc_signal_rand_distr    rd_enum,
                             const float          * __restrict__  noise_scales,
                             std::int32_t * __restrict__          ret_code,
                             const std::uint32_t                  which_generator)
{
        if(__builtin_expect(nullptr==sym_in,0)) {return (-1);}
        if(__builtin_expect(this->m_nsignals>256ull,0)) {return (-2);}
        if(__builtin_expect((this->m_nsignals%2ull)!=0ull,0)) {return (-3);}

#if (AM_BB_CMPLX_SINC_SEQUENCE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif
        constexpr std::size_t max_buff_size{256ull};
        __ATTR_ALIGN__(32)
        __m256 ymm0_r[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm1_r[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm2_r[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm3_r[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm4_r[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm5_r[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm6_r[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm7_r[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm8_r[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm9_r[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm0_i[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm1_i[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm2_i[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm3_i[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm4_i[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm5_i[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm6_i[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm7_i[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm8_i[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm9_i[max_buff_size];
        float scal_rem_r[max_buff_size];
         float scal_rem_i[max_buff_size];
        __m256 vzero{_mm256_setzero_ps()};
        __m256 sum0_r{vzero};
        __m256 sum1_r{sum0_r};
        __m256 sum2_r{sum0_r};
        __m256 sum3_r{sum0_r};
        __m256 sum4_r{sum0_r};
        __m256 sum5_r{sum0_r};
        __m256 sum6_r{sum0_r};
        __m256 sum7_r{sum0_r};
        __m256 sum8_r{sum0_r};
        __m256 sum9_r{sum0_r};
        __m256 sum0_i{sum0_r};
        __m256 sum1_i{sum0_r};
        __m256 sum2_i{sum0_r};
        __m256 sum3_i{sum0_r};
        __m256 sum4_i{sum0_r};
        __m256 sum5_i{sum0_r};
        __m256 sum6_i{sum0_r};
        __m256 sum7_i{sum0_r};
        __m256 sum8_i{sum0_r};
        __m256 sum9_i{sum0_r};
        std::size_t i;
        float sum_rem_r{0.0f};
        float sum_rem_i{0.0f};
        constexpr std::size_t exceeds_L1D_size{32768u};
        const     std::size_t n_samples_halved{this->m_nsamples/2ull};
        switch(which_generator) 
        {
             case 0 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_signal_user_data(sym_in,
                                                                              static_cast<std::uint32_t>(this->m_nsamples),
                                                                              this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 1 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_signal_user_data_u4x(sym_in,
                                                                              static_cast<std::uint32_t>(this->m_nsamples),
                                                                              this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 2 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {    
                        const float scale{noise_scales[k]};
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_noisy_signal_user_data(pdf_params,
                                                                                      scale,sym_in,
                                                                                      static_cast<std::uint32_t>(this->m_nsamples),
                                                                                      this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 3 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {    
                        const float scale{noise_scales[k]};
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_noisy_signal_user_data(pdf_params,
                                                                                      rd_enum,scale,sym_in,
                                                                                      static_cast<std::uint32_t>(this->m_nsamples),
                                                                                      this->m_nK_values.operator[](k));
                   } 
             }
             break;
             default : 
                   return (-3);
        }

        if(this->m_nsamples > exceeds_L1D_size) 
        {
            for(i = 0ull; (i+79ull) < n_samples_halved; i += 80ull) 
            {
              sum0_r = vzero;
              sum0_i = vzero;
              sum1_r = vzero;
              sum1_i = vzero;
              sum2_r = vzero;
              sum2_i = vzero;
              sum3_r = vzero;
              sum3_i = vzero;
              sum4_r = vzero;
              sum4_i = vzero;
              sum5_r = vzero;
              sum5_i = vzero;
              sum6_r = vzero;
              sum6_i = vzero;
              sum7_r = vzero;
              sum7_i = vzero;
              sum8_r = vzero;
              sum8_i = vzero;
              sum9_r = vzero;
              sum9_i = vzero;
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
#if (AM_BB_CMPLX_SINC_SEQUENCE_SOFT_PREFETCH) == 1    
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull],_MM_HINT_T1);
                     
                    
#endif 
                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                      ymm1_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum1_r    = _mm256_add_ps(ymm1_r[j],sum1_r);
                      ymm1_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum1_i    = _mm256_add_ps(ymm1_i[j],sum1_i);
                      ymm2_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum2_r    = _mm256_add_ps(ymm2_r[j],sum2_r);
                      ymm2_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum2_i    = _mm256_add_ps(ymm2_i[j],sum2_i);
                      ymm3_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                      sum3_r    = _mm256_add_ps(ymm3_r[j],sum3_r);
                      ymm3_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                      sum3_i    = _mm256_add_ps(ymm3_i[j],sum3_i);
                      ymm4_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum4_r    = _mm256_add_ps(ymm4_r[j],sum4_r);
                      ymm4_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum4_i    = _mm256_add_ps(ymm4_i[j],sum4_i);
                      ymm5_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+40ull]);
                      sum5_r    = _mm256_add_ps(ymm5_r[j],sum5_r);
                      ymm5_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+40ull]);
                      sum5_i    = _mm256_add_ps(ymm5_i[j],sum5_i);
                      ymm6_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum6_r    = _mm256_add_ps(ymm6_r[j],sum6_r);
                      ymm6_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum6_i    = _mm256_add_ps(ymm6_i[j],sum6_i);
                      ymm7_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+56ull]);
                      sum7_r    = _mm256_add_ps(ymm7_r[j],sum7_r);
                      ymm7_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+56ull]);
                      sum7_i    = _mm256_add_ps(ymm7_i[j],sum7_i);
                      ymm8_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull]);
                      sum8_r    = _mm256_add_ps(ymm8_r[j],sum8_r);
                      ymm8_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull]);
                      sum8_i    = _mm256_add_ps(ymm8_i[j],sum8_i);
                      ymm9_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+72ull]);
                      sum9_r    = _mm256_add_ps(ymm9_r[j],sum9_r);
                      ymm9_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+72ull]);
                      sum9_i    = _mm256_add_ps(ymm9_i[j],sum9_i);
                      

                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum1_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum1_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum2_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum2_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum3_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum3_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum4_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum4_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+40ull],sum5_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+40ull],sum5_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum6_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum6_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+56ull],sum7_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+56ull],sum7_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+64ull],sum8_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+64ull],sum8_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+72ull],sum9_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+72ull],sum9_i);
                     
                 }

             }
        }
        else 
        {
             for(i = 0ull; (i+79ull) < n_samples_halved; i += 80ull) 
        {
              sum0_r = vzero;
              sum0_i = vzero;
              sum1_r = vzero;
              sum1_i = vzero;
              sum2_r = vzero;
              sum2_i = vzero;
              sum3_r = vzero;
              sum3_i = vzero;
              sum4_r = vzero;
              sum4_i = vzero;
              sum5_r = vzero;
              sum5_i = vzero;
              sum6_r = vzero;
              sum6_i = vzero;
              sum7_r = vzero;
              sum7_i = vzero;
              sum8_r = vzero;
              sum8_i = vzero;
              sum9_r = vzero;
              sum9_i = vzero;
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                      ymm1_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum1_r    = _mm256_add_ps(ymm1_r[j],sum1_r);
                      ymm1_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum1_i    = _mm256_add_ps(ymm1_i[j],sum1_i);
                      ymm2_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum2_r    = _mm256_add_ps(ymm2_r[j],sum2_r);
                      ymm2_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum2_i    = _mm256_add_ps(ymm2_i[j],sum2_i);
                      ymm3_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                      sum3_r    = _mm256_add_ps(ymm3_r[j],sum3_r);
                      ymm3_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                      sum3_i    = _mm256_add_ps(ymm3_i[j],sum3_i);
                      ymm4_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum4_r    = _mm256_add_ps(ymm4_r[j],sum4_r);
                      ymm4_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum4_i    = _mm256_add_ps(ymm4_i[j],sum4_i);
                      ymm5_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+40ull]);
                      sum5_r    = _mm256_add_ps(ymm5_r[j],sum5_r);
                      ymm5_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+40ull]);
                      sum5_i    = _mm256_add_ps(ymm5_i[j],sum5_i);
                      ymm6_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum6_r    = _mm256_add_ps(ymm6_r[j],sum6_r);
                      ymm6_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum6_i    = _mm256_add_ps(ymm6_i[j],sum6_i);
                      ymm7_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+56ull]);
                      sum7_r    = _mm256_add_ps(ymm7_r[j],sum7_r);
                      ymm7_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+56ull]);
                      sum7_i    = _mm256_add_ps(ymm7_i[j],sum7_i);
                      ymm8_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull]);
                      sum8_r    = _mm256_add_ps(ymm8_r[j],sum8_r);
                      ymm8_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull]);
                      sum8_i    = _mm256_add_ps(ymm8_i[j],sum8_i);
                      ymm9_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+72ull]);
                      sum9_r    = _mm256_add_ps(ymm9_r[j],sum9_r);
                      ymm9_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+72ull]);
                      sum9_i    = _mm256_add_ps(ymm9_i[j],sum9_i);
                      

                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum1_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum1_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum2_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum2_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum3_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum3_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum4_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum4_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+40ull],sum5_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+40ull],sum5_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum6_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum6_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+56ull],sum7_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+56ull],sum7_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+64ull],sum8_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+64ull],sum8_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+72ull],sum9_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+72ull],sum9_i);
                     
                 }

             }
        }
                 
            for(; (i+63ull) < n_samples_halved; i += 64ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                        sum2_r = vzero;
                        sum2_i = vzero;
                        sum3_r = vzero;
                        sum3_i = vzero;
                        sum4_r = vzero;
                        sum4_i = vzero;
                        sum5_r = vzero;
                        sum5_i = vzero;
                        sum6_r = vzero;
                        sum6_i = vzero;
                        sum7_r = vzero;
                        sum7_i = vzero;
                        
                       
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                      ymm1_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum1_r    = _mm256_add_ps(ymm1_r[j],sum1_r);
                      ymm1_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum1_i    = _mm256_add_ps(ymm1_i[j],sum1_i);
                      ymm2_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum2_r    = _mm256_add_ps(ymm2_r[j],sum2_r);
                      ymm2_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum2_i    = _mm256_add_ps(ymm2_i[j],sum2_i);
                      ymm3_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                      sum3_r    = _mm256_add_ps(ymm3_r[j],sum3_r);
                      ymm3_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                      sum3_i    = _mm256_add_ps(ymm3_i[j],sum3_i);
                      ymm4_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum4_r    = _mm256_add_ps(ymm4_r[j],sum4_r);
                      ymm4_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum4_i    = _mm256_add_ps(ymm4_i[j],sum4_i);
                      ymm5_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+40ull]);
                      sum5_r    = _mm256_add_ps(ymm5_r[j],sum5_r);
                      ymm5_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+40ull]);
                      sum5_i    = _mm256_add_ps(ymm5_i[j],sum5_i);
                      ymm6_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum6_r    = _mm256_add_ps(ymm6_r[j],sum6_r);
                      ymm6_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum6_i    = _mm256_add_ps(ymm6_i[j],sum6_i);
                      ymm7_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+56ull]);
                      sum7_r    = _mm256_add_ps(ymm7_r[j],sum7_r);
                      ymm7_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+56ull]);
                      sum7_i    = _mm256_add_ps(ymm7_i[j],sum7_i);
                                           

                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum1_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum1_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum2_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum2_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum3_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum3_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum4_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum4_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+40ull],sum5_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+40ull],sum5_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum6_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum6_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+56ull],sum7_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+56ull],sum7_i);
                                            
                  }
            }

            
            for(; (i+31ull) < n_samples_halved; i += 32ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                        sum2_r = vzero;
                        sum2_i = vzero;
                        sum3_r = vzero;
                        sum3_i = vzero;
                                            
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                      ymm1_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum1_r    = _mm256_add_ps(ymm1_r[j],sum1_r);
                      ymm1_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum1_i    = _mm256_add_ps(ymm1_i[j],sum1_i);
                      ymm2_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum2_r    = _mm256_add_ps(ymm2_r[j],sum2_r);
                      ymm2_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum2_i    = _mm256_add_ps(ymm2_i[j],sum2_i);
                      ymm3_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                      sum3_r    = _mm256_add_ps(ymm3_r[j],sum3_r);
                      ymm3_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                      sum3_i    = _mm256_add_ps(ymm3_i[j],sum3_i);
                                                              
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum1_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum1_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum2_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum2_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum3_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum3_i);
                                                                 
                  }
            }

               
            for(; (i+15ull) < n_samples_halved; i += 16ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                                                                    
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                      ymm1_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum1_r    = _mm256_add_ps(ymm1_r[j],sum1_r);
                      ymm1_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum1_i    = _mm256_add_ps(ymm1_i[j],sum1_i);
                                                                                    
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum1_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum1_i);
                                                                                       
                  }
            }


            for(; (i+7ull) < n_samples_halved; i += 8ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                                                                                          
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                                                                                                         
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                                                                                                            
                  }
            } 

              for(; (i+0ull) < n_samples_halved; i += 1ull) 
              {
                   sum_rem_r = 0.0f;
                   sum_rem_i = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                   {
                       scal_rem_r[j] = this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i];
                       sum_rem_r += scal_rem_r[j];
                       scal_rem_i[j] = this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i];
                       sum_rem_i += scal_rem_i[j];

                       this->m_re_sequence.m_data[i] = sum_rem_r;
                       this->m_im_sequence.m_data[i] = sum_rem_i;
                   }
              }

#if (AM_BB_CMPLX_SINC_SEQUENCE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif
              return (0);

}

std::int32_t 
gms::radiolocation
::am_bb_cmplx_sinc_signal_sequence_t
::signal_sequence_avx_u6x(  const float * __restrict__           sym_in ,
                             am_bb_cmplx_sinc_signal_pdf_params_t &pdf_params,
                             am_bb_cmplx_sinc_signal_rand_distr    rd_enum,
                             const float          * __restrict__  noise_scales,
                             std::int32_t * __restrict__          ret_code,
                             const std::uint32_t                  which_generator)
{
        if(__builtin_expect(nullptr==sym_in,0)) {return (-1);}
        if(__builtin_expect(this->m_nsignals>256ull,0)) {return (-2);}
        if(__builtin_expect(this->m_nsignals%2ull!=0ull,0)) { return (-3);}

#if (AM_BB_CMPLX_SINC_SEQUENCE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif
        constexpr std::size_t max_buff_size{256ull};
        __ATTR_ALIGN__(32)
        __m256 ymm0_r[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm1_r[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm2_r[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm3_r[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm4_r[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm5_r[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm0_i[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm1_i[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm2_i[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm3_i[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm4_i[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm5_i[max_buff_size];
         float scal_rem_r[max_buff_size];
         float scal_rem_i[max_buff_size];
        __m256 vzero{_mm256_setzero_ps()};
        __m256 sum0_r{vzero};
        __m256 sum1_r{sum0_r};
        __m256 sum2_r{sum0_r};
        __m256 sum3_r{sum0_r};
        __m256 sum4_r{sum0_r};
        __m256 sum5_r{sum0_r};
        __m256 sum0_i{sum0_r};
        __m256 sum1_i{sum0_r};
        __m256 sum2_i{sum0_r};
        __m256 sum3_i{sum0_r};
        __m256 sum4_i{sum0_r};
        __m256 sum5_i{sum0_r};
        std::size_t i;
        float sum_rem_r{0.0f};
        float sum_rem_i{0.0f};
        constexpr std::size_t exceeds_L1D_size{32768u};
        const     std::size_t n_samples_halved{this->m_nsamples/2ull};

        switch(which_generator) 
        {
             case 0 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_signal_user_data(sym_in,
                                                                              static_cast<std::uint32_t>(this->m_nsamples),
                                                                              this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 1 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_signal_user_data_u4x(sym_in,
                                                                              static_cast<std::uint32_t>(this->m_nsamples),
                                                                              this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 2 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {    
                        const float scale{noise_scales[k]};
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_noisy_signal_user_data(pdf_params,
                                                                                      scale,sym_in,
                                                                                      static_cast<std::uint32_t>(this->m_nsamples),
                                                                                      this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 3 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {    
                        const float scale{noise_scales[k]};
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_noisy_signal_user_data(pdf_params,
                                                                                      rd_enum,scale,sym_in,
                                                                                      static_cast<std::uint32_t>(this->m_nsamples),
                                                                                      this->m_nK_values.operator[](k));
                   } 
             }
             break;
             default : 
                   return (-3);
        }

        if(this->m_nsignals > exceeds_L1D_size) 
        {
             for(i = 0ull; (i+47ull) < n_samples_halved; i += 48ull) 
             {
              sum0_r = vzero;
              sum0_i = vzero;
              sum1_r = vzero;
              sum1_i = vzero;
              sum2_r = vzero;
              sum2_i = vzero;
              sum3_r = vzero;
              sum3_i = vzero;
              sum4_r = vzero;
              sum4_i = vzero;
              sum5_r = vzero;
              sum5_i = vzero;
             
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
#if (AM_BB_CMPLX_SINC_SEQUENCE_SOFT_PREFETCH) == 1    
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull],_MM_HINT_T1);
                     
                                          
#endif 
                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                      ymm1_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum1_r    = _mm256_add_ps(ymm1_r[j],sum1_r);
                      ymm1_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum1_i    = _mm256_add_ps(ymm1_i[j],sum1_i);
                      ymm2_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum2_r    = _mm256_add_ps(ymm2_r[j],sum2_r);
                      ymm2_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum2_i    = _mm256_add_ps(ymm2_i[j],sum2_i);
                      ymm3_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                      sum3_r    = _mm256_add_ps(ymm3_r[j],sum3_r);
                      ymm3_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                      sum3_i    = _mm256_add_ps(ymm3_i[j],sum3_i);
                      ymm4_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum4_r    = _mm256_add_ps(ymm4_r[j],sum4_r);
                      ymm4_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum4_i    = _mm256_add_ps(ymm4_i[j],sum4_i);
                      ymm5_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+40ull]);
                      sum5_r    = _mm256_add_ps(ymm5_r[j],sum5_r);
                      ymm5_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+40ull]);
                      sum5_i    = _mm256_add_ps(ymm5_i[j],sum5_i);
                                           

                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum1_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum1_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum2_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum2_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum3_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum3_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum4_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum4_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+40ull],sum5_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+40ull],sum5_i);
                                          
                 }

             }
        }
        else 
        {
             for(i = 0ull; (i+47ull) < n_samples_halved; i += 48ull) 
             {
              sum0_r = vzero;
              sum0_i = vzero;
              sum1_r = vzero;
              sum1_i = vzero;
              sum2_r = vzero;
              sum2_i = vzero;
              sum3_r = vzero;
              sum3_i = vzero;
              sum4_r = vzero;
              sum4_i = vzero;
              sum5_r = vzero;
              sum5_i = vzero;
             
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                      ymm1_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum1_r    = _mm256_add_ps(ymm1_r[j],sum1_r);
                      ymm1_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum1_i    = _mm256_add_ps(ymm1_i[j],sum1_i);
                      ymm2_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum2_r    = _mm256_add_ps(ymm2_r[j],sum2_r);
                      ymm2_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum2_i    = _mm256_add_ps(ymm2_i[j],sum2_i);
                      ymm3_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                      sum3_r    = _mm256_add_ps(ymm3_r[j],sum3_r);
                      ymm3_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                      sum3_i    = _mm256_add_ps(ymm3_i[j],sum3_i);
                      ymm4_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum4_r    = _mm256_add_ps(ymm4_r[j],sum4_r);
                      ymm4_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum4_i    = _mm256_add_ps(ymm4_i[j],sum4_i);
                      ymm5_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+40ull]);
                      sum5_r    = _mm256_add_ps(ymm5_r[j],sum5_r);
                      ymm5_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+40ull]);
                      sum5_i    = _mm256_add_ps(ymm5_i[j],sum5_i);
                                           

                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum1_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum1_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum2_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum2_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum3_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum3_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum4_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum4_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+40ull],sum5_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+40ull],sum5_i);
                                          
                 }

             }
        }
            
            for(; (i+31ull) < n_samples_halved; i += 32ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                        sum2_r = vzero;
                        sum2_i = vzero;
                        sum3_r = vzero;
                        sum3_i = vzero;
                                            
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                      ymm1_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum1_r    = _mm256_add_ps(ymm1_r[j],sum1_r);
                      ymm1_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum1_i    = _mm256_add_ps(ymm1_i[j],sum1_i);
                      ymm2_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum2_r    = _mm256_add_ps(ymm2_r[j],sum2_r);
                      ymm2_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum2_i    = _mm256_add_ps(ymm2_i[j],sum2_i);
                      ymm3_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+24ull]);
                      sum3_r    = _mm256_add_ps(ymm3_r[j],sum3_r);
                      ymm3_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+24ull]);
                      sum3_i    = _mm256_add_ps(ymm3_i[j],sum3_i);
                                                              
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum1_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum1_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum2_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum2_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+24ull],sum3_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+24ull],sum3_i);
                                                                 
                  }
            }

               
            for(; (i+15ull) < n_samples_halved; i += 16ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                                                                    
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                      ymm1_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+8ull]);
                      sum1_r    = _mm256_add_ps(ymm1_r[j],sum1_r);
                      ymm1_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+8ull]);
                      sum1_i    = _mm256_add_ps(ymm1_i[j],sum1_i);
                                                                                    
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+8ull],sum1_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+8ull],sum1_i);
                                                                                       
                  }
            }


            for(; (i+7ull) < n_samples_halved; i += 8ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                                                                                          
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      ymm0_r[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm256_add_ps(ymm0_r[j],sum0_r);
                      ymm0_i[j] = _mm256_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm256_add_ps(ymm0_i[j],sum0_i);
                                                                                                         
                      _mm256_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm256_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                                                                                                            
                  }
            } 

              for(; (i+0ull) < n_samples_halved; i += 1ull) 
              {
                   sum_rem_r = 0.0f;
                   sum_rem_i = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                   {
                       scal_rem_r[j] = this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i];
                       sum_rem_r += scal_rem_r[j];
                       scal_rem_i[j] = this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i];
                       sum_rem_i += scal_rem_i[j];

                       this->m_re_sequence.m_data[i] = sum_rem_r;
                       this->m_im_sequence.m_data[i] = sum_rem_i;
                   }
              }

#if (AM_BB_CMPLX_SINC_SEQUENCE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif
              return (0);

}

std::int32_t 
gms::radiolocation
::am_bb_cmplx_sinc_signal_sequence_t
::signal_sequence_avx512_u16x(  const float * __restrict__           sym_in ,
                                am_bb_cmplx_sinc_signal_pdf_params_t &pdf_params,
                                am_bb_cmplx_sinc_signal_rand_distr    rd_enum,
                                const float          * __restrict__  noise_scales,
                                std::int32_t * __restrict__          ret_code,
                                const std::uint32_t                  which_generator)
{
        if(__builtin_expect(nullptr==sym_in,0)) {return (-1);}
        if(__builtin_expect(this->m_nsignals>256ull,0)) {return (-2);}
        if(__builtin_expect((this->m_nsamples%2ull)!=0ull,0)) {return (-3);}

#if (AM_BB_CMPLX_SINC_SEQUENCE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif
        constexpr std::size_t max_buff_size{256ull};
        __ATTR_ALIGN__(64)
        __m512 zmm0_r[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm1_r[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm2_r[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm3_r[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm4_r[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm5_r[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm6_r[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm7_r[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm8_r[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm9_r[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm10_r[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm11_r[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm12_r[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm13_r[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm14_r[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm15_r[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm0_i[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm1_i[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm2_i[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm3_i[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm4_i[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm5_i[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm6_i[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm7_i[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm8_i[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm9_i[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm10_i[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm11_i[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm12_i[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm13_i[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm14_i[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm15_i[max_buff_size];
         float scal_rem_r[max_buff_size];
         float scal_rem_i[max_buff_size];
        __m512 vzero{_mm512_setzero_ps()};
        __m512 sum0_r{vzero};
        __m512 sum1_r{sum0_r};
        __m512 sum2_r{sum0_r};
        __m512 sum3_r{sum0_r};
        __m512 sum4_r{sum0_r};
        __m512 sum5_r{sum0_r};
        __m512 sum6_r{sum0_r};
        __m512 sum7_r{sum0_r};
        __m512 sum8_r{sum0_r};
        __m512 sum9_r{sum0_r};
        __m512 sum10_r{sum0_r};
        __m512 sum11_r{sum0_r};
        __m512 sum12_r{sum0_r};
        __m512 sum13_r{sum0_r};
        __m512 sum14_r{sum0_r};
        __m512 sum15_r{sum0_r};
        __m512 sum0_i{sum0_r};
        __m512 sum1_i{sum0_r};
        __m512 sum2_i{sum0_r};
        __m512 sum3_i{sum0_r};
        __m512 sum4_i{sum0_r};
        __m512 sum5_i{sum0_r};
        __m512 sum6_i{sum0_r};
        __m512 sum7_i{sum0_r};
        __m512 sum8_i{sum0_r};
        __m512 sum9_i{sum0_r};
        __m512 sum10_i{sum0_r};
        __m512 sum11_i{sum0_r};
        __m512 sum12_i{sum0_r};
        __m512 sum13_i{sum0_r};
        __m512 sum14_i{sum0_r};
        __m512 sum15_i{sum0_r}; 
        std::size_t i;
        float sum_rem_r{0.0f};
        float sum_rem_i{0.0f};
        constexpr std::size_t exceeds_L1D_size{32768u};
        const     std::size_t n_samples_halved{this->m_nsamples/2ull};

        switch(which_generator) 
        {
             case 0 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_signal_user_data(sym_in,
                                                                              static_cast<std::uint32_t>(this->m_nsamples),
                                                                              this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 1 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_signal_user_data_u4x(sym_in,
                                                                              static_cast<std::uint32_t>(this->m_nsamples),
                                                                              this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 2 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {    
                        const float scale{noise_scales[k]};
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_noisy_signal_user_data(pdf_params,
                                                                                      scale,sym_in,
                                                                                      static_cast<std::uint32_t>(this->m_nsamples),
                                                                                      this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 3 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {    
                        const float scale{noise_scales[k]};
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_noisy_signal_user_data(pdf_params,
                                                                                      rd_enum,scale,sym_in,
                                                                                      static_cast<std::uint32_t>(this->m_nsamples),
                                                                                      this->m_nK_values.operator[](k));
                   } 
             }
             break;
             default : 
                   return (-3);
        }

        if(this->m_nsamples > exceeds_L1D_size) 
        {
           for(i = 0ull; (i+255ull) < n_samples_halved; i += 256ull) 
           {
              sum0_r = vzero;
              sum0_i = vzero;
              sum1_r = vzero;
              sum1_i = vzero;
              sum2_r = vzero;
              sum2_i = vzero;
              sum3_r = vzero;
              sum3_i = vzero;
              sum4_r = vzero;
              sum4_i = vzero;
              sum5_r = vzero;
              sum5_i = vzero;
              sum6_r = vzero;
              sum6_i = vzero;
              sum7_r = vzero;
              sum7_i = vzero;
              sum8_r = vzero;
              sum8_i = vzero;
              sum9_r = vzero;
              sum9_i = vzero;
              sum10_r = vzero;
              sum10_i = vzero;
              sum11_r = vzero;
              sum11_i = vzero;
              sum12_r = vzero;
              sum12_i = vzero;
              sum13_r = vzero;
              sum13_i = vzero;
              sum14_r = vzero;
              sum14_i = vzero;
              sum15_r = vzero;
              sum15_i = vzero;
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
#if (AM_BB_CMPLX_SINC_SEQUENCE_SOFT_PREFETCH) == 1    
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+80ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+80ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+96ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+96ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+112ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+112ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+128ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+128ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+144ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+144ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+160ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+160ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+176ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+176ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+192ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+192ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+208ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+208ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+224ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+224ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+240ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+240ull],_MM_HINT_T1);
#endif 
                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);
                      zmm1_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum1_r    = _mm512_add_ps(zmm1_r[j],sum1_r);
                      zmm1_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum1_i    = _mm512_add_ps(zmm1_i[j],sum1_i);
                      zmm2_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum2_r    = _mm512_add_ps(zmm2_r[j],sum2_r);
                      zmm2_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum2_i    = _mm512_add_ps(zmm2_i[j],sum2_i);
                      zmm3_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum3_r    = _mm512_add_ps(zmm3_r[j],sum3_r);
                      zmm3_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum3_i    = _mm512_add_ps(zmm3_i[j],sum3_i);
                      zmm4_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull]);
                      sum4_r    = _mm512_add_ps(zmm4_r[j],sum4_r);
                      zmm4_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull]);
                      sum4_i    = _mm512_add_ps(zmm4_i[j],sum4_i);
                      zmm5_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+80ull]);
                      sum5_r    = _mm512_add_ps(zmm5_r[j],sum5_r);
                      zmm5_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+80ull]);
                      sum5_i    = _mm512_add_ps(zmm5_i[j],sum5_i);
                      zmm6_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+96ull]);
                      sum6_r    = _mm512_add_ps(zmm6_r[j],sum6_r);
                      zmm6_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+96ull]);
                      sum6_i    = _mm512_add_ps(zmm6_i[j],sum6_i);
                      zmm7_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+112ull]);
                      sum7_r    = _mm512_add_ps(zmm7_r[j],sum7_r);
                      zmm7_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+112ull]);
                      sum7_i    = _mm512_add_ps(zmm7_i[j],sum7_i);
                      zmm8_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+128ull]);
                      sum8_r    = _mm512_add_ps(zmm8_r[j],sum8_r);
                      zmm8_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+128ull]);
                      sum8_i    = _mm512_add_ps(zmm8_i[j],sum8_i);
                      zmm9_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+144ull]);
                      sum9_r    = _mm512_add_ps(zmm9_r[j],sum9_r);
                      zmm9_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+144ull]);
                      sum9_i    = _mm512_add_ps(zmm9_i[j],sum9_i);
                      zmm10_r[j]= _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+160ull]);
                      sum10_r   = _mm512_add_ps(zmm10_r[j],sum10_r);
                      zmm10_i[j]= _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+160ull]);
                      sum10_i   = _mm512_add_ps(zmm10_i[j],sum10_i);
                      zmm11_r[j]= _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+176ull]);
                      sum11_r    = _mm512_add_ps(zmm11_r[j],sum11_r);
                      zmm11_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+176ull]);
                      sum11_i    = _mm512_add_ps(zmm11_i[j],sum11_i);
                      zmm12_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+192ull]);
                      sum12_r    = _mm512_add_ps(zmm12_r[j],sum12_r);
                      zmm12_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+192ull]);
                      sum12_i    = _mm512_add_ps(zmm12_i[j],sum12_i);
                      zmm13_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+208ull]);
                      sum13_r    = _mm512_add_ps(zmm13_r[j],sum13_r);
                      zmm13_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+208ull]);
                      sum13_i    = _mm512_add_ps(zmm13_i[j],sum13_i);
                      zmm14_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+224ull]);
                      sum14_r    = _mm512_add_ps(zmm14_r[j],sum14_r);
                      zmm14_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+224ull]);
                      sum14_i    = _mm512_add_ps(zmm14_i[j],sum14_i);
                      zmm15_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+240ull]);
                      sum15_r    = _mm512_add_ps(zmm15_r[j],sum15_r);
                      zmm15_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+240ull]);
                      sum15_i    = _mm512_add_ps(zmm15_i[j],sum15_i);

                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum1_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum1_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum2_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum2_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum3_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum3_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+64ull],sum4_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+64ull],sum4_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+80ull],sum5_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+80ull],sum5_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+96ull],sum6_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+96ull],sum6_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+112ull],sum7_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+112ull],sum7_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+128ull],sum8_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+128ull],sum8_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+144ull],sum9_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+144ull],sum9_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+160ull],sum10_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+160ull],sum10_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+176ull],sum11_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+176ull],sum11_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+192ull],sum12_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+192ull],sum12_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+208ull],sum13_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+208ull],sum13_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+224ull],sum14_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+224ull],sum14_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+240ull],sum15_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+240ull],sum15_i);
                 }

             }
        }
        else 
        {
            for(i = 0ull; (i+255ull) < n_samples_halved; i += 256ull) 
            {
              sum0_r = vzero;
              sum0_i = vzero;
              sum1_r = vzero;
              sum1_i = vzero;
              sum2_r = vzero;
              sum2_i = vzero;
              sum3_r = vzero;
              sum3_i = vzero;
              sum4_r = vzero;
              sum4_i = vzero;
              sum5_r = vzero;
              sum5_i = vzero;
              sum6_r = vzero;
              sum6_i = vzero;
              sum7_r = vzero;
              sum7_i = vzero;
              sum8_r = vzero;
              sum8_i = vzero;
              sum9_r = vzero;
              sum9_i = vzero;
              sum10_r = vzero;
              sum10_i = vzero;
              sum11_r = vzero;
              sum11_i = vzero;
              sum12_r = vzero;
              sum12_i = vzero;
              sum13_r = vzero;
              sum13_i = vzero;
              sum14_r = vzero;
              sum14_i = vzero;
              sum15_r = vzero;
              sum15_i = vzero;
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);
                      zmm1_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum1_r    = _mm512_add_ps(zmm1_r[j],sum1_r);
                      zmm1_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum1_i    = _mm512_add_ps(zmm1_i[j],sum1_i);
                      zmm2_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum2_r    = _mm512_add_ps(zmm2_r[j],sum2_r);
                      zmm2_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum2_i    = _mm512_add_ps(zmm2_i[j],sum2_i);
                      zmm3_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum3_r    = _mm512_add_ps(zmm3_r[j],sum3_r);
                      zmm3_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum3_i    = _mm512_add_ps(zmm3_i[j],sum3_i);
                      zmm4_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull]);
                      sum4_r    = _mm512_add_ps(zmm4_r[j],sum4_r);
                      zmm4_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull]);
                      sum4_i    = _mm512_add_ps(zmm4_i[j],sum4_i);
                      zmm5_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+80ull]);
                      sum5_r    = _mm512_add_ps(zmm5_r[j],sum5_r);
                      zmm5_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+80ull]);
                      sum5_i    = _mm512_add_ps(zmm5_i[j],sum5_i);
                      zmm6_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+96ull]);
                      sum6_r    = _mm512_add_ps(zmm6_r[j],sum6_r);
                      zmm6_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+96ull]);
                      sum6_i    = _mm512_add_ps(zmm6_i[j],sum6_i);
                      zmm7_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+112ull]);
                      sum7_r    = _mm512_add_ps(zmm7_r[j],sum7_r);
                      zmm7_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+112ull]);
                      sum7_i    = _mm512_add_ps(zmm7_i[j],sum7_i);
                      zmm8_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+128ull]);
                      sum8_r    = _mm512_add_ps(zmm8_r[j],sum8_r);
                      zmm8_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+128ull]);
                      sum8_i    = _mm512_add_ps(zmm8_i[j],sum8_i);
                      zmm9_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+144ull]);
                      sum9_r    = _mm512_add_ps(zmm9_r[j],sum9_r);
                      zmm9_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+144ull]);
                      sum9_i    = _mm512_add_ps(zmm9_i[j],sum9_i);
                      zmm10_r[j]= _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+160ull]);
                      sum10_r   = _mm512_add_ps(zmm10_r[j],sum10_r);
                      zmm10_i[j]= _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+160ull]);
                      sum10_i   = _mm512_add_ps(zmm10_i[j],sum10_i);
                      zmm11_r[j]= _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+176ull]);
                      sum11_r    = _mm512_add_ps(zmm11_r[j],sum11_r);
                      zmm11_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+176ull]);
                      sum11_i    = _mm512_add_ps(zmm11_i[j],sum11_i);
                      zmm12_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+192ull]);
                      sum12_r    = _mm512_add_ps(zmm12_r[j],sum12_r);
                      zmm12_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+192ull]);
                      sum12_i    = _mm512_add_ps(zmm12_i[j],sum12_i);
                      zmm13_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+208ull]);
                      sum13_r    = _mm512_add_ps(zmm13_r[j],sum13_r);
                      zmm13_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+208ull]);
                      sum13_i    = _mm512_add_ps(zmm13_i[j],sum13_i);
                      zmm14_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+224ull]);
                      sum14_r    = _mm512_add_ps(zmm14_r[j],sum14_r);
                      zmm14_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+224ull]);
                      sum14_i    = _mm512_add_ps(zmm14_i[j],sum14_i);
                      zmm15_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+240ull]);
                      sum15_r    = _mm512_add_ps(zmm15_r[j],sum15_r);
                      zmm15_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+240ull]);
                      sum15_i    = _mm512_add_ps(zmm15_i[j],sum15_i);

                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum1_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum1_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum2_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum2_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum3_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum3_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+64ull],sum4_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+64ull],sum4_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+80ull],sum5_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+80ull],sum5_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+96ull],sum6_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+96ull],sum6_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+112ull],sum7_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+112ull],sum7_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+128ull],sum8_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+128ull],sum8_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+144ull],sum9_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+144ull],sum9_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+160ull],sum10_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+160ull],sum10_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+176ull],sum11_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+176ull],sum11_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+192ull],sum12_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+192ull],sum12_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+208ull],sum13_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+208ull],sum13_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+224ull],sum14_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+224ull],sum14_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+240ull],sum15_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+240ull],sum15_i);
                 }

             }
        }

              for(; (i+191ull) < n_samples_halved; i += 192ull) 
              {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                        sum2_r = vzero;
                        sum2_i = vzero;
                        sum3_r = vzero;
                        sum3_i = vzero;
                        sum4_r = vzero;
                        sum4_i = vzero;
                        sum5_r = vzero;
                        sum5_i = vzero;
                        sum6_r = vzero;
                        sum6_i = vzero;
                        sum7_r = vzero;
                        sum7_i = vzero;
                        sum8_r = vzero;
                        sum8_i = vzero;
                        sum9_r = vzero;
                        sum9_i = vzero;
                        sum10_r = vzero;
                        sum10_i = vzero;
                        sum11_r = vzero;
                        sum11_i = vzero;
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);
                      zmm1_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum1_r    = _mm512_add_ps(zmm1_r[j],sum1_r);
                      zmm1_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum1_i    = _mm512_add_ps(zmm1_i[j],sum1_i);
                      zmm2_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum2_r    = _mm512_add_ps(zmm2_r[j],sum2_r);
                      zmm2_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum2_i    = _mm512_add_ps(zmm2_i[j],sum2_i);
                      zmm3_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum3_r    = _mm512_add_ps(zmm3_r[j],sum3_r);
                      zmm3_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum3_i    = _mm512_add_ps(zmm3_i[j],sum3_i);
                      zmm4_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull]);
                      sum4_r    = _mm512_add_ps(zmm4_r[j],sum4_r);
                      zmm4_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull]);
                      sum4_i    = _mm512_add_ps(zmm4_i[j],sum4_i);
                      zmm5_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+80ull]);
                      sum5_r    = _mm512_add_ps(zmm5_r[j],sum5_r);
                      zmm5_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+80ull]);
                      sum5_i    = _mm512_add_ps(zmm5_i[j],sum5_i);
                      zmm6_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+96ull]);
                      sum6_r    = _mm512_add_ps(zmm6_r[j],sum6_r);
                      zmm6_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+96ull]);
                      sum6_i    = _mm512_add_ps(zmm6_i[j],sum6_i);
                      zmm7_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+112ull]);
                      sum7_r    = _mm512_add_ps(zmm7_r[j],sum7_r);
                      zmm7_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+112ull]);
                      sum7_i    = _mm512_add_ps(zmm7_i[j],sum7_i);
                      zmm8_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+128ull]);
                      sum8_r    = _mm512_add_ps(zmm8_r[j],sum8_r);
                      zmm8_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+128ull]);
                      sum8_i    = _mm512_add_ps(zmm8_i[j],sum8_i);
                      zmm9_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+144ull]);
                      sum9_r    = _mm512_add_ps(zmm9_r[j],sum9_r);
                      zmm9_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+144ull]);
                      sum9_i    = _mm512_add_ps(zmm9_i[j],sum9_i);
                      zmm10_r[j]= _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+160ull]);
                      sum10_r   = _mm512_add_ps(zmm10_r[j],sum10_r);
                      zmm10_i[j]= _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+160ull]);
                      sum10_i   = _mm512_add_ps(zmm10_i[j],sum10_i);
                      zmm11_r[j]= _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+176ull]);
                      sum11_r    = _mm512_add_ps(zmm11_r[j],sum11_r);
                      zmm11_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+176ull]);
                      sum11_i    = _mm512_add_ps(zmm11_i[j],sum11_i);

                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum1_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum1_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum2_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum2_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum3_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum3_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+64ull],sum4_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+64ull],sum4_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+80ull],sum5_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+80ull],sum5_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+96ull],sum6_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+96ull],sum6_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+112ull],sum7_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+112ull],sum7_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+128ull],sum8_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+128ull],sum8_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+144ull],sum9_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+144ull],sum9_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+160ull],sum10_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+160ull],sum10_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+176ull],sum11_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+176ull],sum11_i);
                  }
            }

            for(; (i+159ull) < n_samples_halved; i += 160ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                        sum2_r = vzero;
                        sum2_i = vzero;
                        sum3_r = vzero;
                        sum3_i = vzero;
                        sum4_r = vzero;
                        sum4_i = vzero;
                        sum5_r = vzero;
                        sum5_i = vzero;
                        sum6_r = vzero;
                        sum6_i = vzero;
                        sum7_r = vzero;
                        sum7_i = vzero;
                        sum8_r = vzero;
                        sum8_i = vzero;
                        sum9_r = vzero;
                        sum9_i = vzero;
                       
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      
                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);
                      zmm1_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum1_r    = _mm512_add_ps(zmm1_r[j],sum1_r);
                      zmm1_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum1_i    = _mm512_add_ps(zmm1_i[j],sum1_i);
                      zmm2_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum2_r    = _mm512_add_ps(zmm2_r[j],sum2_r);
                      zmm2_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum2_i    = _mm512_add_ps(zmm2_i[j],sum2_i);
                      zmm3_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum3_r    = _mm512_add_ps(zmm3_r[j],sum3_r);
                      zmm3_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum3_i    = _mm512_add_ps(zmm3_i[j],sum3_i);
                      zmm4_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull]);
                      sum4_r    = _mm512_add_ps(zmm4_r[j],sum4_r);
                      zmm4_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull]);
                      sum4_i    = _mm512_add_ps(zmm4_i[j],sum4_i);
                      zmm5_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+80ull]);
                      sum5_r    = _mm512_add_ps(zmm5_r[j],sum5_r);
                      zmm5_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+80ull]);
                      sum5_i    = _mm512_add_ps(zmm5_i[j],sum5_i);
                      zmm6_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+96ull]);
                      sum6_r    = _mm512_add_ps(zmm6_r[j],sum6_r);
                      zmm6_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+96ull]);
                      sum6_i    = _mm512_add_ps(zmm6_i[j],sum6_i);
                      zmm7_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+112ull]);
                      sum7_r    = _mm512_add_ps(zmm7_r[j],sum7_r);
                      zmm7_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+112ull]);
                      sum7_i    = _mm512_add_ps(zmm7_i[j],sum7_i);
                      zmm8_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+128ull]);
                      sum8_r    = _mm512_add_ps(zmm8_r[j],sum8_r);
                      zmm8_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+128ull]);
                      sum8_i    = _mm512_add_ps(zmm8_i[j],sum8_i);
                      zmm9_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+144ull]);
                      sum9_r    = _mm512_add_ps(zmm9_r[j],sum9_r);
                      zmm9_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+144ull]);
                      sum9_i    = _mm512_add_ps(zmm9_i[j],sum9_i);

                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum1_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum1_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum2_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum2_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum3_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum3_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+64ull],sum4_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+64ull],sum4_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+80ull],sum5_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+80ull],sum5_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+96ull],sum6_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+96ull],sum6_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+112ull],sum7_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+112ull],sum7_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+128ull],sum8_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+128ull],sum8_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+144ull],sum9_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+144ull],sum9_i);
                      
                  }
            }

            
            for(; (i+127ull) < n_samples_halved; i += 128ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                        sum2_r = vzero;
                        sum2_i = vzero;
                        sum3_r = vzero;
                        sum3_i = vzero;
                        sum4_r = vzero;
                        sum4_i = vzero;
                        sum5_r = vzero;
                        sum5_i = vzero;
                        sum6_r = vzero;
                        sum6_i = vzero;
                        sum7_r = vzero;
                        sum7_i = vzero;
                        
                       
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      
                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);
                      zmm1_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum1_r    = _mm512_add_ps(zmm1_r[j],sum1_r);
                      zmm1_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum1_i    = _mm512_add_ps(zmm1_i[j],sum1_i);
                      zmm2_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum2_r    = _mm512_add_ps(zmm2_r[j],sum2_r);
                      zmm2_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum2_i    = _mm512_add_ps(zmm2_i[j],sum2_i);
                      zmm3_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum3_r    = _mm512_add_ps(zmm3_r[j],sum3_r);
                      zmm3_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum3_i    = _mm512_add_ps(zmm3_i[j],sum3_i);
                      zmm4_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull]);
                      sum4_r    = _mm512_add_ps(zmm4_r[j],sum4_r);
                      zmm4_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull]);
                      sum4_i    = _mm512_add_ps(zmm4_i[j],sum4_i);
                      zmm5_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+80ull]);
                      sum5_r    = _mm512_add_ps(zmm5_r[j],sum5_r);
                      zmm5_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+80ull]);
                      sum5_i    = _mm512_add_ps(zmm5_i[j],sum5_i);
                      zmm6_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+96ull]);
                      sum6_r    = _mm512_add_ps(zmm6_r[j],sum6_r);
                      zmm6_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+96ull]);
                      sum6_i    = _mm512_add_ps(zmm6_i[j],sum6_i);
                      zmm7_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+112ull]);
                      sum7_r    = _mm512_add_ps(zmm7_r[j],sum7_r);
                      zmm7_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+112ull]);
                      sum7_i    = _mm512_add_ps(zmm7_i[j],sum7_i);                    

                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum1_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum1_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum2_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum2_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum3_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum3_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+64ull],sum4_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+64ull],sum4_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+80ull],sum5_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+80ull],sum5_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+96ull],sum6_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+96ull],sum6_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+112ull],sum7_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+112ull],sum7_i);
                                            
                  }
            }

            
            for(; (i+63ull) < n_samples_halved; i += 64ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                        sum2_r = vzero;
                        sum2_i = vzero;
                        sum3_r = vzero;
                        sum3_i = vzero;
                                            
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);
                      zmm1_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum1_r    = _mm512_add_ps(zmm1_r[j],sum1_r);
                      zmm1_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum1_i    = _mm512_add_ps(zmm1_i[j],sum1_i);
                      zmm2_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum2_r    = _mm512_add_ps(zmm2_r[j],sum2_r);
                      zmm2_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum2_i    = _mm512_add_ps(zmm2_i[j],sum2_i);
                      zmm3_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum3_r    = _mm512_add_ps(zmm3_r[j],sum3_r);
                      zmm3_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum3_i    = _mm512_add_ps(zmm3_i[j],sum3_i);
                                                              
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum1_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum1_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum2_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum2_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum3_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum3_i);
                                                                 
                  }
            }

               
            for(; (i+31ull) < n_samples_halved; i += 32ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                                                                    
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      
                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);
                      zmm1_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum1_r    = _mm512_add_ps(zmm1_r[j],sum1_r);
                      zmm1_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum1_i    = _mm512_add_ps(zmm1_i[j],sum1_i);                                                              
                      
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum1_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum1_i);
                                                                                       
                  }
            }


            for(; (i+15ull) < n_samples_halved; i += 16ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                                                                                          
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      
                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);  

                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                                                                                                            
                  }
            } 

              for(; (i+0ull) < n_samples_halved; i += 1ull) 
              {
                   sum_rem_r = 0.0f;
                   sum_rem_i = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                   {
                       scal_rem_r[j] = this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i];
                       sum_rem_r += scal_rem_r[j];
                       scal_rem_i[j] = this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i];
                       sum_rem_i += scal_rem_i[j];

                       this->m_re_sequence.m_data[i] = sum_rem_r;
                       this->m_im_sequence.m_data[i] = sum_rem_i;
                   }
              }

#if (AM_BB_CMPLX_SINC_SEQUENCE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif
              return (0);

}

std::int32_t 
gms::radiolocation
::am_bb_cmplx_sinc_signal_sequence_t
::signal_sequence_avx512_u10x(  const float * __restrict__           sym_in ,
                                am_bb_cmplx_sinc_signal_pdf_params_t &pdf_params,
                                am_bb_cmplx_sinc_signal_rand_distr    rd_enum,
                                const float          * __restrict__  noise_scales,
                                std::int32_t * __restrict__          ret_code,
                                const std::uint32_t                  which_generator)
{
        if(__builtin_expect(nullptr==sym_in,0)) {return (-1);}
        if(__builtin_expect(this->m_nsignals>256ull,0)) {return (-2);}
        if(__builtin_expect(this->m_nsamples%2ull!=0ull,0)) { return (-3);}

#if (AM_BB_CMPLX_SINC_SEQUENCE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif
        constexpr std::size_t max_buff_size{256ull};
        __ATTR_ALIGN__(64)
        __m512 zmm0_r[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm1_r[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm2_r[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm3_r[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm4_r[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm5_r[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm6_r[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm7_r[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm8_r[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm9_r[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm0_i[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm1_i[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm2_i[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm3_i[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm4_i[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm5_i[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm6_i[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm7_i[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm8_i[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm9_i[max_buff_size];
         float scal_rem_r[max_buff_size];
         float scal_rem_i[max_buff_size];
        __m512 vzero{_mm512_setzero_ps()};
        __m512 sum0_r{vzero};
        __m512 sum1_r{sum0_r};
        __m512 sum2_r{sum0_r};
        __m512 sum3_r{sum0_r};
        __m512 sum4_r{sum0_r};
        __m512 sum5_r{sum0_r};
        __m512 sum6_r{sum0_r};
        __m512 sum7_r{sum0_r};
        __m512 sum8_r{sum0_r};
        __m512 sum9_r{sum0_r};
        __m512 sum0_i{sum0_r};
        __m512 sum1_i{sum0_r};
        __m512 sum2_i{sum0_r};
        __m512 sum3_i{sum0_r};
        __m512 sum4_i{sum0_r};
        __m512 sum5_i{sum0_r};
        __m512 sum6_i{sum0_r};
        __m512 sum7_i{sum0_r};
        __m512 sum8_i{sum0_r};
        __m512 sum9_i{sum0_r};
        std::size_t i;
        float sum_rem_r{0.0f};
        float sum_rem_i{0.0f};
        constexpr std::size_t exceeds_L1D_size{32768u};
        const     std::size_t n_samples_halved{this->m_nsamples/2ull};

        switch(which_generator) 
        {
             case 0 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_signal_user_data(sym_in,
                                                                              static_cast<std::uint32_t>(this->m_nsamples),
                                                                              this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 1 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_signal_user_data_u4x(sym_in,
                                                                              static_cast<std::uint32_t>(this->m_nsamples),
                                                                              this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 2 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {    
                        const float scale{noise_scales[k]};
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_noisy_signal_user_data(pdf_params,
                                                                                      scale,sym_in,
                                                                                      static_cast<std::uint32_t>(this->m_nsamples),
                                                                                      this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 3 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {    
                        const float scale{noise_scales[k]};
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_noisy_signal_user_data(pdf_params,
                                                                                      rd_enum,scale,sym_in,
                                                                                      static_cast<std::uint32_t>(this->m_nsamples),
                                                                                      this->m_nK_values.operator[](k));
                   } 
             }
             break;
             default : 
                   return (-3);
        }

        if(this->m_nsamples > exceeds_L1D_size) 
        {
             for(i = 0ull; (i+159ull) < n_samples_halved; i += 160ull) 
            {
              sum0_r = vzero;
              sum0_i = vzero;
              sum1_r = vzero;
              sum1_i = vzero;
              sum2_r = vzero;
              sum2_i = vzero;
              sum3_r = vzero;
              sum3_i = vzero;
              sum4_r = vzero;
              sum4_i = vzero;
              sum5_r = vzero;
              sum5_i = vzero;
              sum6_r = vzero;
              sum6_i = vzero;
              sum7_r = vzero;
              sum7_i = vzero;
              sum8_r = vzero;
              sum8_i = vzero;
              sum9_r = vzero;
              sum9_i = vzero;
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
#if (AM_BB_CMPLX_SINC_SEQUENCE_SOFT_PREFETCH) == 1    
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+80ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+80ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+96ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+96ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+112ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+112ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+128ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+128ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+144ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+144ull],_MM_HINT_T1);
                     
#endif 
                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);
                      zmm1_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum1_r    = _mm512_add_ps(zmm1_r[j],sum1_r);
                      zmm1_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum1_i    = _mm512_add_ps(zmm1_i[j],sum1_i);
                      zmm2_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum2_r    = _mm512_add_ps(zmm2_r[j],sum2_r);
                      zmm2_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum2_i    = _mm512_add_ps(zmm2_i[j],sum2_i);
                      zmm3_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum3_r    = _mm512_add_ps(zmm3_r[j],sum3_r);
                      zmm3_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum3_i    = _mm512_add_ps(zmm3_i[j],sum3_i);
                      zmm4_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull]);
                      sum4_r    = _mm512_add_ps(zmm4_r[j],sum4_r);
                      zmm4_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull]);
                      sum4_i    = _mm512_add_ps(zmm4_i[j],sum4_i);
                      zmm5_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+80ull]);
                      sum5_r    = _mm512_add_ps(zmm5_r[j],sum5_r);
                      zmm5_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+80ull]);
                      sum5_i    = _mm512_add_ps(zmm5_i[j],sum5_i);
                      zmm6_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+96ull]);
                      sum6_r    = _mm512_add_ps(zmm6_r[j],sum6_r);
                      zmm6_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+96ull]);
                      sum6_i    = _mm512_add_ps(zmm6_i[j],sum6_i);
                      zmm7_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+112ull]);
                      sum7_r    = _mm512_add_ps(zmm7_r[j],sum7_r);
                      zmm7_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+112ull]);
                      sum7_i    = _mm512_add_ps(zmm7_i[j],sum7_i);
                      zmm8_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+128ull]);
                      sum8_r    = _mm512_add_ps(zmm8_r[j],sum8_r);
                      zmm8_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+128ull]);
                      sum8_i    = _mm512_add_ps(zmm8_i[j],sum8_i);
                      zmm9_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+144ull]);
                      sum9_r    = _mm512_add_ps(zmm9_r[j],sum9_r);
                      zmm9_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+144ull]);
                      sum9_i    = _mm512_add_ps(zmm9_i[j],sum9_i);
                     

                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum1_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum1_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum2_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum2_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum3_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum3_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+64ull],sum4_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+64ull],sum4_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+80ull],sum5_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+80ull],sum5_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+96ull],sum6_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+96ull],sum6_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+112ull],sum7_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+112ull],sum7_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+128ull],sum8_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+128ull],sum8_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+144ull],sum9_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+144ull],sum9_i);
                    
                 }

             }
        }
        else 
        {
             for(i = 0ull; (i+159ull) < n_samples_halved; i += 160ull) 
        {
              sum0_r = vzero;
              sum0_i = vzero;
              sum1_r = vzero;
              sum1_i = vzero;
              sum2_r = vzero;
              sum2_i = vzero;
              sum3_r = vzero;
              sum3_i = vzero;
              sum4_r = vzero;
              sum4_i = vzero;
              sum5_r = vzero;
              sum5_i = vzero;
              sum6_r = vzero;
              sum6_i = vzero;
              sum7_r = vzero;
              sum7_i = vzero;
              sum8_r = vzero;
              sum8_i = vzero;
              sum9_r = vzero;
              sum9_i = vzero;
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);
                      zmm1_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum1_r    = _mm512_add_ps(zmm1_r[j],sum1_r);
                      zmm1_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum1_i    = _mm512_add_ps(zmm1_i[j],sum1_i);
                      zmm2_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum2_r    = _mm512_add_ps(zmm2_r[j],sum2_r);
                      zmm2_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum2_i    = _mm512_add_ps(zmm2_i[j],sum2_i);
                      zmm3_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum3_r    = _mm512_add_ps(zmm3_r[j],sum3_r);
                      zmm3_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum3_i    = _mm512_add_ps(zmm3_i[j],sum3_i);
                      zmm4_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull]);
                      sum4_r    = _mm512_add_ps(zmm4_r[j],sum4_r);
                      zmm4_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull]);
                      sum4_i    = _mm512_add_ps(zmm4_i[j],sum4_i);
                      zmm5_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+80ull]);
                      sum5_r    = _mm512_add_ps(zmm5_r[j],sum5_r);
                      zmm5_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+80ull]);
                      sum5_i    = _mm512_add_ps(zmm5_i[j],sum5_i);
                      zmm6_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+96ull]);
                      sum6_r    = _mm512_add_ps(zmm6_r[j],sum6_r);
                      zmm6_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+96ull]);
                      sum6_i    = _mm512_add_ps(zmm6_i[j],sum6_i);
                      zmm7_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+112ull]);
                      sum7_r    = _mm512_add_ps(zmm7_r[j],sum7_r);
                      zmm7_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+112ull]);
                      sum7_i    = _mm512_add_ps(zmm7_i[j],sum7_i);
                      zmm8_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+128ull]);
                      sum8_r    = _mm512_add_ps(zmm8_r[j],sum8_r);
                      zmm8_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+128ull]);
                      sum8_i    = _mm512_add_ps(zmm8_i[j],sum8_i);
                      zmm9_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+144ull]);
                      sum9_r    = _mm512_add_ps(zmm9_r[j],sum9_r);
                      zmm9_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+144ull]);
                      sum9_i    = _mm512_add_ps(zmm9_i[j],sum9_i);
                     

                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum1_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum1_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum2_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum2_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum3_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum3_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+64ull],sum4_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+64ull],sum4_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+80ull],sum5_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+80ull],sum5_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+96ull],sum6_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+96ull],sum6_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+112ull],sum7_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+112ull],sum7_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+128ull],sum8_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+128ull],sum8_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+144ull],sum9_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+144ull],sum9_i);
                    
                 }

             }
        }
            

            for(; (i+127ull) < n_samples_halved; i += 128ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                        sum2_r = vzero;
                        sum2_i = vzero;
                        sum3_r = vzero;
                        sum3_i = vzero;
                        sum4_r = vzero;
                        sum4_i = vzero;
                        sum5_r = vzero;
                        sum5_i = vzero;
                        sum6_r = vzero;
                        sum6_i = vzero;
                        sum7_r = vzero;
                        sum7_i = vzero;
                        
                       
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      
                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);
                      zmm1_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum1_r    = _mm512_add_ps(zmm1_r[j],sum1_r);
                      zmm1_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum1_i    = _mm512_add_ps(zmm1_i[j],sum1_i);
                      zmm2_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum2_r    = _mm512_add_ps(zmm2_r[j],sum2_r);
                      zmm2_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum2_i    = _mm512_add_ps(zmm2_i[j],sum2_i);
                      zmm3_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum3_r    = _mm512_add_ps(zmm3_r[j],sum3_r);
                      zmm3_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum3_i    = _mm512_add_ps(zmm3_i[j],sum3_i);
                      zmm4_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull]);
                      sum4_r    = _mm512_add_ps(zmm4_r[j],sum4_r);
                      zmm4_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull]);
                      sum4_i    = _mm512_add_ps(zmm4_i[j],sum4_i);
                      zmm5_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+80ull]);
                      sum5_r    = _mm512_add_ps(zmm5_r[j],sum5_r);
                      zmm5_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+80ull]);
                      sum5_i    = _mm512_add_ps(zmm5_i[j],sum5_i);
                      zmm6_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+96ull]);
                      sum6_r    = _mm512_add_ps(zmm6_r[j],sum6_r);
                      zmm6_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+96ull]);
                      sum6_i    = _mm512_add_ps(zmm6_i[j],sum6_i);
                      zmm7_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+112ull]);
                      sum7_r    = _mm512_add_ps(zmm7_r[j],sum7_r);
                      zmm7_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+112ull]);
                      sum7_i    = _mm512_add_ps(zmm7_i[j],sum7_i);                    

                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum1_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum1_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum2_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum2_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum3_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum3_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+64ull],sum4_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+64ull],sum4_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+80ull],sum5_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+80ull],sum5_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+96ull],sum6_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+96ull],sum6_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+112ull],sum7_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+112ull],sum7_i);
                                            
                  }
            }

            
            for(; (i+63ull) < n_samples_halved; i += 64ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                        sum2_r = vzero;
                        sum2_i = vzero;
                        sum3_r = vzero;
                        sum3_i = vzero;
                                            
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);
                      zmm1_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum1_r    = _mm512_add_ps(zmm1_r[j],sum1_r);
                      zmm1_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum1_i    = _mm512_add_ps(zmm1_i[j],sum1_i);
                      zmm2_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum2_r    = _mm512_add_ps(zmm2_r[j],sum2_r);
                      zmm2_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum2_i    = _mm512_add_ps(zmm2_i[j],sum2_i);
                      zmm3_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum3_r    = _mm512_add_ps(zmm3_r[j],sum3_r);
                      zmm3_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum3_i    = _mm512_add_ps(zmm3_i[j],sum3_i);
                                                              
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum1_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum1_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum2_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum2_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum3_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum3_i);
                                                                 
                  }
            }

               
            for(; (i+31ull) < n_samples_halved; i += 32ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                                                                    
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      
                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);
                      zmm1_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum1_r    = _mm512_add_ps(zmm1_r[j],sum1_r);
                      zmm1_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum1_i    = _mm512_add_ps(zmm1_i[j],sum1_i);                                                              
                      
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum1_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum1_i);
                                                                                       
                  }
            }


            for(; (i+15ull) < n_samples_halved; i += 16ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                                                                                          
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      
                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);  

                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                                                                                                            
                  }
            } 

              for(; (i+0ull) < n_samples_halved; i += 1ull) 
              {
                   sum_rem_r = 0.0f;
                   sum_rem_i = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                   {
                       scal_rem_r[j] = this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i];
                       sum_rem_r += scal_rem_r[j];
                       scal_rem_i[j] = this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i];
                       sum_rem_i += scal_rem_i[j];

                       this->m_re_sequence.m_data[i] = sum_rem_r;
                       this->m_im_sequence.m_data[i] = sum_rem_i;
                   }
              }

#if (AM_BB_CMPLX_SINC_SEQUENCE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif
              return (0);

}

std::int32_t 
gms::radiolocation
::am_bb_cmplx_sinc_signal_sequence_t
::signal_sequence_avx512_u6x(  const float * __restrict__           sym_in ,
                                am_bb_cmplx_sinc_signal_pdf_params_t &pdf_params,
                                am_bb_cmplx_sinc_signal_rand_distr    rd_enum,
                                const float          * __restrict__  noise_scales,
                                std::int32_t * __restrict__          ret_code,
                                const std::uint32_t                  which_generator)
{
        if(__builtin_expect(nullptr==sym_in,0)) {return (-1);}
        if(__builtin_expect(this->m_nsignals>256ull,0)) {return (-2);}
        if(__builtin_expect((this->m_nsignals!=2ull)!=0ull,0)) { return (-3);}

#if (AM_BB_CMPLX_SINC_SEQUENCE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif
        constexpr std::size_t max_buff_size{256ull};
        __ATTR_ALIGN__(64)
        __m512 zmm0_r[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm1_r[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm2_r[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm3_r[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm4_r[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm5_r[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm0_i[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm1_i[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm2_i[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm3_i[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm4_i[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm5_i[max_buff_size];
         float scal_rem_r[max_buff_size];
         float scal_rem_i[max_buff_size];
        __m512 vzero{_mm512_setzero_ps()};
        __m512 sum0_r{vzero};
        __m512 sum1_r{sum0_r};
        __m512 sum2_r{sum0_r};
        __m512 sum3_r{sum0_r};
        __m512 sum4_r{sum0_r};
        __m512 sum5_r{sum0_r};
        __m512 sum0_i{sum0_r};
        __m512 sum1_i{sum0_r};
        __m512 sum2_i{sum0_r};
        __m512 sum3_i{sum0_r};
        __m512 sum4_i{sum0_r};
        __m512 sum5_i{sum0_r};
        std::size_t i;
        float sum_rem_r{0.0f};
        float sum_rem_i{0.0f};
        constexpr std::size_t exceeds_L1D_size{32768u};
        const     std::size_t n_samples_halved{this->m_nsamples/2ull};

        switch(which_generator) 
        {
             case 0 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_signal_user_data(sym_in,
                                                                              static_cast<std::uint32_t>(this->m_nsamples),
                                                                              this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 1 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_signal_user_data_u4x(sym_in,
                                                                              static_cast<std::uint32_t>(this->m_nsamples),
                                                                              this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 2 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {    
                        const float scale{noise_scales[k]};
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_noisy_signal_user_data(pdf_params,
                                                                                      scale,sym_in,
                                                                                      static_cast<std::uint32_t>(this->m_nsamples),
                                                                                      this->m_nK_values.operator[](k));
                   }
             }
             break;
             case 3 : 
             {
                   for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                   {    
                        const float scale{noise_scales[k]};
                        ret_code[k] = this->m_bb_cmplx_sinc_signals.operator[](k).create_noisy_signal_user_data(pdf_params,
                                                                                      rd_enum,scale,sym_in,
                                                                                      static_cast<std::uint32_t>(this->m_nsamples),
                                                                                      this->m_nK_values.operator[](k));
                   } 
             }
             break;
             default : 
                   return (-3);
        }

        if(this->m_nsamples > exceeds_L1D_size) 
        {
             for(i = 0ull; (i+95ull) < n_samples_halved; i += 96ull) 
            {
              sum0_r = vzero;
              sum0_i = vzero;
              sum1_r = vzero;
              sum1_i = vzero;
              sum2_r = vzero;
              sum2_i = vzero;
              sum3_r = vzero;
              sum3_i = vzero;
              sum4_r = vzero;
              sum4_i = vzero;
              sum5_r = vzero;
              sum5_i = vzero;
              
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
#if (AM_BB_CMPLX_SINC_SEQUENCE_SOFT_PREFETCH) == 1    
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+80ull],_MM_HINT_T1);
                      _mm_prefetch((const char*)&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+80ull],_MM_HINT_T1);
                                          
#endif 
                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);
                      zmm1_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum1_r    = _mm512_add_ps(zmm1_r[j],sum1_r);
                      zmm1_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum1_i    = _mm512_add_ps(zmm1_i[j],sum1_i);
                      zmm2_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum2_r    = _mm512_add_ps(zmm2_r[j],sum2_r);
                      zmm2_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum2_i    = _mm512_add_ps(zmm2_i[j],sum2_i);
                      zmm3_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum3_r    = _mm512_add_ps(zmm3_r[j],sum3_r);
                      zmm3_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum3_i    = _mm512_add_ps(zmm3_i[j],sum3_i);
                      zmm4_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull]);
                      sum4_r    = _mm512_add_ps(zmm4_r[j],sum4_r);
                      zmm4_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull]);
                      sum4_i    = _mm512_add_ps(zmm4_i[j],sum4_i);
                      zmm5_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+80ull]);
                      sum5_r    = _mm512_add_ps(zmm5_r[j],sum5_r);
                      zmm5_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+80ull]);
                      sum5_i    = _mm512_add_ps(zmm5_i[j],sum5_i);
                      
                     

                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum1_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum1_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum2_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum2_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum3_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum3_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+64ull],sum4_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+64ull],sum4_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+80ull],sum5_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+80ull],sum5_i);
                                         
                 }

             }
        }
       else 
       {
             for(i = 0ull; (i+95ull) < n_samples_halved; i += 96ull) 
        {
              sum0_r = vzero;
              sum0_i = vzero;
              sum1_r = vzero;
              sum1_i = vzero;
              sum2_r = vzero;
              sum2_i = vzero;
              sum3_r = vzero;
              sum3_i = vzero;
              sum4_r = vzero;
              sum4_i = vzero;
              sum5_r = vzero;
              sum5_i = vzero;
              
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);
                      zmm1_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum1_r    = _mm512_add_ps(zmm1_r[j],sum1_r);
                      zmm1_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum1_i    = _mm512_add_ps(zmm1_i[j],sum1_i);
                      zmm2_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum2_r    = _mm512_add_ps(zmm2_r[j],sum2_r);
                      zmm2_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum2_i    = _mm512_add_ps(zmm2_i[j],sum2_i);
                      zmm3_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum3_r    = _mm512_add_ps(zmm3_r[j],sum3_r);
                      zmm3_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum3_i    = _mm512_add_ps(zmm3_i[j],sum3_i);
                      zmm4_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+64ull]);
                      sum4_r    = _mm512_add_ps(zmm4_r[j],sum4_r);
                      zmm4_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+64ull]);
                      sum4_i    = _mm512_add_ps(zmm4_i[j],sum4_i);
                      zmm5_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+80ull]);
                      sum5_r    = _mm512_add_ps(zmm5_r[j],sum5_r);
                      zmm5_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+80ull]);
                      sum5_i    = _mm512_add_ps(zmm5_i[j],sum5_i);
                      
                     

                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum1_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum1_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum2_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum2_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum3_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum3_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+64ull],sum4_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+64ull],sum4_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+80ull],sum5_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+80ull],sum5_i);
                                         
                 }

             }
       }
                        
            for(; (i+63ull) < n_samples_halved; i += 64ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                        sum2_r = vzero;
                        sum2_i = vzero;
                        sum3_r = vzero;
                        sum3_i = vzero;
                                            
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);
                      zmm1_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum1_r    = _mm512_add_ps(zmm1_r[j],sum1_r);
                      zmm1_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum1_i    = _mm512_add_ps(zmm1_i[j],sum1_i);
                      zmm2_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+32ull]);
                      sum2_r    = _mm512_add_ps(zmm2_r[j],sum2_r);
                      zmm2_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+32ull]);
                      sum2_i    = _mm512_add_ps(zmm2_i[j],sum2_i);
                      zmm3_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+48ull]);
                      sum3_r    = _mm512_add_ps(zmm3_r[j],sum3_r);
                      zmm3_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+48ull]);
                      sum3_i    = _mm512_add_ps(zmm3_i[j],sum3_i);
                                                              
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum1_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum1_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+32ull],sum2_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+32ull],sum2_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+48ull],sum3_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+48ull],sum3_i);
                                                                 
                  }
            }

               
            for(; (i+31ull) < n_samples_halved; i += 32ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                        sum1_r = vzero;
                        sum1_i = vzero;
                                                                    
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      
                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);
                      zmm1_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+16ull]);
                      sum1_r    = _mm512_add_ps(zmm1_r[j],sum1_r);
                      zmm1_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+16ull]);
                      sum1_i    = _mm512_add_ps(zmm1_i[j],sum1_i);                                                              
                      
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+16ull],sum1_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+16ull],sum1_i);
                                                                                       
                  }
            }


            for(; (i+15ull) < n_samples_halved; i += 16ull) 
            {
                        sum0_r = vzero;
                        sum0_i = vzero;
                                                                                          
                  for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                  {
                           
                      
                      zmm0_r[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i+0ull]);
                      sum0_r    = _mm512_add_ps(zmm0_r[j],sum0_r);
                      zmm0_i[j] = _mm512_loadu_ps(&this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i+0ull]);
                      sum0_i    = _mm512_add_ps(zmm0_i[j],sum0_i);  

                      _mm512_storeu_ps(&this->m_re_sequence.m_data[i+0ull],sum0_r);
                      _mm512_storeu_ps(&this->m_im_sequence.m_data[i+0ull],sum0_i);
                                                                                                            
                  }
            } 

              for(; (i+0ull) < n_samples_halved; i += 1ull) 
              {
                   sum_rem_r = 0.0f;
                   sum_rem_i = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
                   {
                       scal_rem_r[j] = this->m_bb_cmplx_sinc_signals.operator[](j).get_I_channel()[i];
                       sum_rem_r += scal_rem_r[j];
                       scal_rem_i[j] = this->m_bb_cmplx_sinc_signals.operator[](j).get_Q_channel()[i];
                       sum_rem_i += scal_rem_i[j];

                       this->m_re_sequence.m_data[i] = sum_rem_r;
                       this->m_im_sequence.m_data[i] = sum_rem_i;
                   }
              }
              
#if (AM_BB_CMPLX_SINC_SEQUENCE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif
              return (0);

}






