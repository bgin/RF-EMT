
#include <immintrin.h>
#include <fstream>
#include <iomanip>
#include "GMS_am_bb_sine_sequence.h"


gms::radiolocation
::am_bb_sine_signal_sequence_t
::am_bb_sine_signal_sequence_t(const std::size_t nsamples,
                              const std::size_t nsignals,
                              const std::vector<std::uint32_t> & nK_values,
                              const std::vector<float> &n_values,
                              const std::vector<float> &A_values,
                              const std::vector<float> &P_values) 
:
m_nsamples{nsamples},
m_nsignals{nsignals}
{
    this->m_nK_values        = std::vector<std::uint32_t>(nK_values);
    this->m_n_values         = std::vector<float>(n_values);
    this->m_A_values         = std::vector<float>(A_values);
    this->m_P_values         = std::vector<float>(P_values);
    this->m_bb_sine_signals  = std::vector<am_bb_sine_signal_t>(this->m_nsignals);
    this->m_vsummed_sequence = std::valarray<float>(this->m_nsamples);
    for(std::size_t i{0ull}; i != this->m_nsignals; ++i) 
    {
        const std::uint32_t nK_vals{this->m_nK_values.operator[](i)};
        const float n_vals{this->m_n_values.operator[](i)};
        const float A_vals{this->m_A_values.operator[](i)};
        const float P_vals{this->m_P_values.operator[](i)};
        this->m_bb_sine_signals.operator[](i) = am_bb_sine_signal_t(this->m_nsamples,
                                                                     nK_vals,n_vals,
                                                                     A_vals,P_vals);
    }
}

gms::radiolocation
::am_bb_sine_signal_sequence_t
::am_bb_sine_signal_sequence_t(am_bb_sine_signal_sequence_t &&other) noexcept(true) 
:
m_nsamples{std::move(other.m_nsamples)},
m_nsignals{std::move(other.m_nsignals)},
m_nK_values{std::move(other.m_nK_values)},
m_n_values{std::move(other.m_n_values)},
m_A_values{std::move(other.m_A_values)},
m_P_values{std::move(other.m_P_values)},
m_bb_sine_signals{std::move(other.m_bb_sine_signals)},
m_vsummed_sequence{std::move(other.m_vsummed_sequence)}
{

}

gms::radiolocation
::am_bb_sine_signal_sequence_t
::am_bb_sine_signal_sequence_t(const am_bb_sine_signal_sequence_t &other) noexcept(false) 
:
m_nsamples{other.m_nsamples},
m_nsignals{other.m_nsignals},
m_nK_values{other.m_nK_values},
m_n_values{other.m_n_values},
m_A_values{other.m_A_values},
m_P_values{other.m_P_values},
m_bb_sine_signals{other.m_bb_sine_signals},
m_vsummed_sequence{other.m_vsummed_sequence}
{

}

gms::radiolocation
::am_bb_sine_signal_sequence_t
::~am_bb_sine_signal_sequence_t() noexcept(false)
{

}

gms::radiolocation
::am_bb_sine_signal_sequence_t &
gms::radiolocation 
::am_bb_sine_signal_sequence_t
::operator=(const am_bb_sine_signal_sequence_t &other) noexcept(false) 
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_nsamples = other.m_nsamples;
    this->m_nsignals = other.m_nsignals;
    this->m_nK_values.operator=(other.m_nK_values);
    this->m_n_values.operator=(other.m_n_values);
    this->m_A_values.operator=(other.m_A_values);
    this->m_P_values.operator=(other.m_P_values);
    this->m_bb_sine_signals.operator=(other.m_bb_sine_signals);
    this->m_vsummed_sequence.operator=(other.m_vsummed_sequence);
    return (*this);
}

gms::radiolocation
::am_bb_sine_signal_sequence_t &
gms::radiolocation 
::am_bb_sine_signal_sequence_t
::operator=(am_bb_sine_signal_sequence_t &&other) noexcept(true) 
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_nsamples = std::move(other.m_nsamples);
    this->m_nsignals = std::move(other.m_nsignals);
    this->m_nK_values.operator=(std::move(other.m_nK_values));
    this->m_n_values.operator=(std::move(other.m_n_values));
    this->m_A_values.operator=(std::move(other.m_A_values));
    this->m_P_values.operator=(std::move(other.m_P_values));
    this->m_bb_sine_signals.operator=(std::move(other.m_bb_sine_signals));
    this->m_vsummed_sequence.operator=(std::move(other.m_vsummed_sequence));
    return (*this);
}

void 
gms::radiolocation
::am_bb_sine_signal_sequence_t
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
::am_bb_sine_signal_sequence_t
::signal_sequence_sse42_u16x(const float * __restrict__ sym_in ,
                               std::int32_t * __restrict__ ret_code,
                               const std::uint32_t which_generator) 
{
        if(__builtin_expect(nullptr==sym_in,0)) {return (-1);}
        if(__builtin_expect(this->m_nsignals>256ull,0)) {return (-2);}

        constexpr std::size_t max_buff_size{256ull};
        __ATTR_ALIGN__(16)
        __m128 xmm0[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm1[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm2[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm3[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm4[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm5[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm6[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm7[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm8[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm9[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm10[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm11[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm12[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm13[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm14[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm15[max_buff_size];
        float scal_rem[max_buff_size];
        __m128 vzero{_mm_setzero_ps()};
        __m128 sum0{vzero};
        __m128 sum1{sum0};
        __m128 sum2{sum0};
        __m128 sum3{sum0};
        __m128 sum4{sum0};
        __m128 sum5{sum0};
        __m128 sum6{sum0};
        __m128 sum7{sum0};
        __m128 sum8{sum0};
        __m128 sum9{sum0};
        __m128 sum10{sum0};
        __m128 sum11{sum0};
        __m128 sum12{sum0};
        __m128 sum13{sum0};
        __m128 sum14{sum0};
        __m128 sum15{sum0};    
        std::size_t i;
        float sum_rem{0.0f};

        switch (which_generator)
        {
             case 0 :
             {
                  for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                  {
            
                            ret_code[k] = this->m_bb_sine_signals.operator[](k).create_signal_user_data(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
                  }      
             }     
             break;
             case 1 :
             {
                  for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                  {
            
                            ret_code[k] = this->m_bb_sine_signals.operator[](k).create_signal_user_data_u4x(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
                  } 
             }
             break;
             default : 
                  return (-3);     
         }

        for(i = 0ull; (i+63ull) < this->m_nsamples; i += 64ull) 
        {
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;
              sum4  = vzero;
              sum5  = vzero;
              sum6  = vzero;
              sum7  = vzero;
              sum8  = vzero;
              sum9  = vzero;
              sum10 = vzero;
              sum11 = vzero;
              sum12 = vzero;
              sum13 = vzero;
              sum14 = vzero;
              sum15 = vzero;  
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
              {
#if (AM_BB_SINE_SEQUENCE_SOFT_PREFETCH) == 1
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+4ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+12ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+20ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+28ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+36ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+40ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+44ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+52ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+56ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+60ull],_MM_HINT_T1);
#endif 
               xmm0[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
               sum0  =  _mm_add_ps(xmm0[j],sum0);
               xmm1[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
               sum1  =  _mm_add_ps(xmm1[j],sum1);
               xmm2[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
               sum2  =  _mm_add_ps(xmm2[j],sum2);
               xmm3[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+12ull]);
               sum3  =  _mm_add_ps(xmm3[j],sum3);
               xmm4[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
               sum4  =  _mm_add_ps(xmm4[j],sum4);
               xmm5[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+20ull]);
               sum5  =  _mm_add_ps(xmm5[j],sum5);
               xmm6[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
               sum6  =  _mm_add_ps(xmm6[j],sum6);
               xmm7[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+28ull]);
               sum7  =  _mm_add_ps(xmm7[j],sum7);
               xmm8[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
               sum8  =  _mm_add_ps(xmm8[j],sum8);
               xmm9[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+36ull]);
               sum9  =  _mm_add_ps(xmm9[j],sum9);
               xmm10[j] =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+40ull]);
               sum10 =  _mm_add_ps(xmm10[j],sum10);
               xmm11[j] =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+44ull]);
               sum11 =  _mm_add_ps(xmm11[j],sum11);
               xmm1[j] =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
               sum12 =  _mm_add_ps(xmm12[j],sum12);
               xmm13[j] =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+52ull]);
               sum13 =  _mm_add_ps(xmm13[j],sum13);
               xmm14[j] =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+56ull]);
               sum14 =  _mm_add_ps(xmm14[j],sum14);
               xmm15[j] =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+60ull]);
               sum15 =  _mm_add_ps(xmm15[j],sum15);
               
               _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum2);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+12ull],sum3);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum4);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+20ull],sum5);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum6);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+28ull],sum7);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum8);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+36ull],sum9);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+40ull],sum10);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+44ull],sum11);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum12);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+52ull],sum13);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+56ull],sum14);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+60ull],sum15);
            }
             
        }

        for(; (i+39ull) < this->m_nsamples; i += 40ull) 
        {
           
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;
              sum4  = vzero;
              sum5  = vzero;
              sum6  = vzero;
              sum7  = vzero;
              sum8  = vzero;
              sum9  = vzero;
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

                 xmm0[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                 sum0  =  _mm_add_ps(xmm0[j],sum0);
                 xmm1[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
                 sum1  =  _mm_add_ps(xmm1[j],sum1);
                 xmm2[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                 sum2  =  _mm_add_ps(xmm2[j],sum2);
                 xmm3[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+12ull]);
                 sum3  =  _mm_add_ps(xmm3[j],sum3);
                 xmm4[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                 sum4  =  _mm_add_ps(xmm4[j],sum4);
                 xmm5[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+20ull]);
                 sum5  =  _mm_add_ps(xmm5[j],sum5);
                 xmm6[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                 sum6  =  _mm_add_ps(xmm6[j],sum6);
                 xmm7[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+28ull]);
                 sum7  =  _mm_add_ps(xmm7[j],sum7);
                 xmm8[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                 sum8  =  _mm_add_ps(xmm8[j],sum8);
                 xmm9[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+36ull]);
                 sum9  =  _mm_add_ps(xmm9[j],sum9);

                 _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                 _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);
                 _mm_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum2);
                 _mm_storeu_ps(&this->m_vsummed_sequence[i+12ull],sum3);
                 _mm_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum4);
                 _mm_storeu_ps(&this->m_vsummed_sequence[i+20ull],sum5);
                 _mm_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum6);
                 _mm_storeu_ps(&this->m_vsummed_sequence[i+28ull],sum7);
                 _mm_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum8);
                 _mm_storeu_ps(&this->m_vsummed_sequence[i+36ull],sum9);    
            }
             
        }

        for(; (i+31ull) < this->m_nsamples; i += 32ull) 
        {
            
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;
              sum4  = vzero;
              sum5  = vzero;
              sum6  = vzero;
              sum7  = vzero;          
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

                 xmm0[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                 sum0  =  _mm_add_ps(xmm0[j],sum0);
                 xmm1[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
                 sum1  =  _mm_add_ps(xmm1[j],sum1);
                 xmm2[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                 sum2  =  _mm_add_ps(xmm2[j],sum2);
                 xmm3[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+12ull]);
                 sum3  =  _mm_add_ps(xmm3[j],sum3);
                 xmm4[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                 sum4  =  _mm_add_ps(xmm4[j],sum4);
                 xmm5[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+20ull]);
                 sum5  =  _mm_add_ps(xmm5[j],sum5);
                 xmm6[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                 sum6  =  _mm_add_ps(xmm6[j],sum6);
                 xmm7[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+28ull]);
                 sum7  =  _mm_add_ps(xmm7[j],sum7);

                  _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum2);
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+12ull],sum3);
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum4);
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+20ull],sum5);
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum6);
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+28ull],sum7);
                                        
            }
            
        }

        for(; (i+15ull) < this->m_nsamples; i += 16ull) 
        {
            
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;                      
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {
                 xmm0[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                 sum0  =  _mm_add_ps(xmm0[j],sum0);
                 xmm1[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
                 sum1  =  _mm_add_ps(xmm1[j],sum1);
                 xmm2[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                 sum2  =  _mm_add_ps(xmm2[j],sum2);
                 xmm3[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+12ull]);
                 sum3  =  _mm_add_ps(xmm3[j],sum3);

                  _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum2);
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+12ull],sum3);                                      
            }
            
        }

        for(; (i+7ull) < this->m_nsamples; i += 8ull) 
        {
           
              sum0  = vzero;
              sum1  = vzero;                                   
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {
                 xmm0[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                 sum0  =  _mm_add_ps(xmm0[j],sum0);
                 xmm1[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
                 sum1  =  _mm_add_ps(xmm1[j],sum1);

                 _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                 _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);                                   
            }
             
        }

        for(; (i+3ull) < this->m_nsamples; i += 4ull) 
        {
             
              sum0  = vzero;                                          
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {
                 xmm0[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                 sum0  =  _mm_add_ps(xmm0[j],sum0);

                  _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);                                              
            }
            
        }

        for(; (i+0ull) < this->m_nsamples; i += 1ull) 
        {
            sum_rem = 0.0f;
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               scal_rem[j]   =  this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i];
               sum_rem  +=  scal_rem[j];

                this->m_vsummed_sequence[i] = sum_rem;                                               
            }
           
        }

        return (0);
}


std::int32_t  
gms::radiolocation
::am_bb_sine_signal_sequence_t
::signal_sequence_sse42_u10x(const float * __restrict__ sym_in ,
                               std::int32_t * __restrict__ ret_code,
                               const std::uint32_t which_generator) 
{
        if(__builtin_expect(nullptr==sym_in,0)) {return (-1);}
        if(__builtin_expect(this->m_nsignals>256ull,0)) {return (-2);}

        constexpr std::size_t max_buff_size{256ull};
        __ATTR_ALIGN__(16)
        __m128 xmm0[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm1[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm2[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm3[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm4[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm5[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm6[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm7[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm8[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm9[max_buff_size];
        float scal_rem[max_buff_size];
        __m128 vzero{_mm_setzero_ps()};
        __m128 sum0{vzero};
        __m128 sum1{sum0};
        __m128 sum2{sum0};
        __m128 sum3{sum0};
        __m128 sum4{sum0};
        __m128 sum5{sum0};
        __m128 sum6{sum0};
        __m128 sum7{sum0};
        __m128 sum8{sum0};
        __m128 sum9{sum0};
        float sum_rem{0.0f};    
        std::size_t i;

        switch (which_generator)
        {
             case 0 :
             {
                  for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                  {
            
                            ret_code[k] = this->m_bb_sine_signals.operator[](k).create_signal_user_data(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
                  }      
             }     
             break;
             case 1 :
             {
                  for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                  {
            
                            ret_code[k] = this->m_bb_sine_signals.operator[](k).create_signal_user_data_u4x(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
                  } 
             }
             break;
             default : 
                  return (-3);     
         }

        for(i = 0ull; (i+39ull) < this->m_nsamples; i += 40ull) 
        {
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;
              sum4  = vzero;
              sum5  = vzero;
              sum6  = vzero;
              sum7  = vzero;
              sum8  = vzero;
              sum9  = vzero;       
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {
#if (AM_BB_SINE_SEQUENCE_SOFT_PREFETCH) == 1
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+4ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+12ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+20ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+28ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+36ull],_MM_HINT_T1);
               
#endif 
               xmm0[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
               sum0  =  _mm_add_ps(xmm0[j],sum0);
               xmm1[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
               sum1  =  _mm_add_ps(xmm1[j],sum1);
               xmm2[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
               sum2  =  _mm_add_ps(xmm2[j],sum2);
               xmm3[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+12ull]);
               sum3  =  _mm_add_ps(xmm3[j],sum3);
               xmm4[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
               sum4  =  _mm_add_ps(xmm4[j],sum4);
               xmm5[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+20ull]);
               sum5  =  _mm_add_ps(xmm5[j],sum5);
               xmm6[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
               sum6  =  _mm_add_ps(xmm6[j],sum6);
               xmm7[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+28ull]);
               sum7  =  _mm_add_ps(xmm7[j],sum7);
               xmm8[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
               sum8  =  _mm_add_ps(xmm8[j],sum8);
               xmm9[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+36ull]);
               sum9  =  _mm_add_ps(xmm9[j],sum9);
               
               _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum2);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+12ull],sum3);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum4);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+20ull],sum5);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum6);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+28ull],sum7);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum8);
               _mm_storeu_ps(&this->m_vsummed_sequence[i+36ull],sum9);
            }
                         
        }
    
        for(; (i+31ull) < this->m_nsamples; i += 32ull) 
        {
            
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;
              sum4  = vzero;
              sum5  = vzero;
              sum6  = vzero;
              sum7  = vzero;         
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

                 xmm0[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                 sum0  =  _mm_add_ps(xmm0[j],sum0);
                 xmm1[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
                 sum1  =  _mm_add_ps(xmm1[j],sum1);
                 xmm2[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                 sum2  =  _mm_add_ps(xmm2[j],sum2);
                 xmm3[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+12ull]);
                 sum3  =  _mm_add_ps(xmm3[j],sum3);
                 xmm4[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                 sum4  =  _mm_add_ps(xmm4[j],sum4);
                 xmm5[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+20ull]);
                 sum5  =  _mm_add_ps(xmm5[j],sum5);
                 xmm6[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                 sum6  =  _mm_add_ps(xmm6[j],sum6);
                 xmm7[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+28ull]);
                 sum7  =  _mm_add_ps(xmm7[j],sum7);

                 _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                 _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);
                 _mm_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum2);
                 _mm_storeu_ps(&this->m_vsummed_sequence[i+12ull],sum3);
                 _mm_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum4);
                 _mm_storeu_ps(&this->m_vsummed_sequence[i+20ull],sum5);
                 _mm_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum6);
                 _mm_storeu_ps(&this->m_vsummed_sequence[i+28ull],sum7);                       
            }
             
        }

        for(; (i+15ull) < this->m_nsamples; i += 16ull) 
        {          
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;                     
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {
                  xmm0[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                  sum0  =  _mm_add_ps(xmm0[j],sum0);
                  xmm1[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
                  sum1  =  _mm_add_ps(xmm1[j],sum1);
                  xmm2[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                  sum2  =  _mm_add_ps(xmm2[j],sum2);
                  xmm3[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+12ull]);
                  sum3  =  _mm_add_ps(xmm3[j],sum3);

                  _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum2);
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+12ull],sum3);
                                                       
            }            
        }

        for(; (i+7ull) < this->m_nsamples; i += 8ull) 
        {            
               sum0  = vzero;
               sum1  = vzero;                                 
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

                  xmm0[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                  sum0  =  _mm_add_ps(xmm0[j],sum0);
                  xmm1[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
                  sum1  =  _mm_add_ps(xmm1[j],sum1);

                  _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);                           
            }             
        }

        for(; (i+3ull) < this->m_nsamples; i += 4ull) 
        {
           
               sum0  = vzero;                                         
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {
                  xmm0[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                  sum0  =  _mm_add_ps(xmm0[j],sum0);

                  _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);                                               
            }            
        }

        for(; (i+0ull) < this->m_nsamples; i += 1ull) 
        {
            sum_rem = 0.0f;
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               scal_rem[j]       =  this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i];
               sum_rem  +=  scal_rem[j];

               this->m_vsummed_sequence[i] = sum_rem;                                             
            }
           
        }

        return (0);
}


std::int32_t  
gms::radiolocation
::am_bb_sine_signal_sequence_t
::signal_sequence_sse42_u6x(const float * __restrict__ sym_in ,
                               std::int32_t * __restrict__ ret_code,
                               const std::uint32_t which_generator) // 0 for rolled, 1 for unrolled.
{
        if(__builtin_expect(nullptr==sym_in,0)) {return (-1);}
        if(__builtin_expect(this->m_nsignals>256ull,0)) {return (-2);}

        constexpr std::size_t max_buff_size{256ull};
        __ATTR_ALIGN__(16)
        __m128 xmm0[max_buff_size];
        __ATTR_ALIGN__(16)
        __m128 xmm1[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm2[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm3[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm4[max_buff_size];
         __ATTR_ALIGN__(16)
        __m128 xmm5[max_buff_size];
        float scal_rem[max_buff_size];
        __m128 vzero{_mm_setzero_ps()};
        __m128 sum0{vzero};
        __m128 sum1{sum0};
        __m128 sum2{sum0};
        __m128 sum3{sum0};
        __m128 sum4{sum0};
        __m128 sum5{sum0};              
        float sum_rem{0.0f};
        std::size_t i;

        switch (which_generator)
        {
             case 0 :
             {
                  for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                  {
            
                            ret_code[k] = this->m_bb_sine_signals.operator[](k).create_signal_user_data(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
                  }      
             }     
             break;
             case 1 :
             {
                  for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                  {
            
                            ret_code[k] = this->m_bb_sine_signals.operator[](k).create_signal_user_data_u4x(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
                  } 
             }
             break;
             default : 
                  return (-3);     
         }

        for(i = 0ull; (i+23ull) < this->m_nsamples; i += 24ull) 
        {
           
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;
              sum4  = vzero;
              sum5  = vzero;        
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {
#if (AM_BB_SINE_SEQUENCE_SOFT_PREFETCH) == 1
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+4ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+12ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+20ull],_MM_HINT_T1);
                             
#endif 
                 xmm0[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                 sum0  =  _mm_add_ps(xmm0[j],sum0);
                 xmm1[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
                 sum1  =  _mm_add_ps(xmm1[j],sum1);
                 xmm2[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                 sum2  =  _mm_add_ps(xmm2[j],sum2);
                 xmm3[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+12ull]);
                 sum3  =  _mm_add_ps(xmm3[j],sum3);
                 xmm4[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                 sum4  =  _mm_add_ps(xmm4[j],sum4);
                 xmm5[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+20ull]);
                 sum5  =  _mm_add_ps(xmm5[j],sum5);
                            
                _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);
                _mm_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum2);
                _mm_storeu_ps(&this->m_vsummed_sequence[i+12ull],sum3);
                _mm_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum4);
                _mm_storeu_ps(&this->m_vsummed_sequence[i+20ull],sum5);
            }            
                     
        }  
       
        for(; (i+15ull) < this->m_nsamples; i += 16ull) 
        {            
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;                     
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {
                 xmm0[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                 sum0  =  _mm_add_ps(xmm0[j],sum0);
                 xmm1[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
                 sum1  =  _mm_add_ps(xmm1[j],sum1);
                 xmm2[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                 sum2  =  _mm_add_ps(xmm2[j],sum2);
                 xmm3[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+12ull]);
                 sum3  =  _mm_add_ps(xmm3[j],sum3);
               
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum2);
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+12ull],sum3);                       
            }           
        }

        for(; (i+7ull) < this->m_nsamples; i += 8ull) 
        {
            
            sum0  = vzero;
            sum1  = vzero;                                
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {
                  xmm0[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                  sum0  =  _mm_add_ps(xmm0[j],sum0);
                  xmm1[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
                  sum1  =  _mm_add_ps(xmm1[j],sum1);

                  _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                  _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);                                    
            }             
        }

        for(; (i+3ull) < this->m_nsamples; i += 4ull) 
        {           
               sum0  = vzero;                                         
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {
                  xmm0[j]  =  _mm_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                  sum0  =  _mm_add_ps(xmm0[j],sum0);

                   _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);                                                 
            }            
        }

        for(; (i+0ull) < this->m_nsamples; i += 1ull) 
        {
            sum_rem = 0.0f;
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               scal_rem[j]       =  this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i];
               sum_rem           +=  scal_rem[j];

               this->m_vsummed_sequence[i] = sum_rem; 
            }            
        }

        return (0);
}

std::int32_t
gms::radiolocation 
::am_bb_sine_signal_sequence_t
::signal_sequence_avx_u16x(const float * __restrict__ sym_in ,
                             std::int32_t * __restrict__ ret_codes,
                             const std::uint32_t which_generator) 
{
        if(__builtin_expect(nullptr==sym_in,0)) {return (-1);}
        if(__builtin_expect(this->m_nsignals>256ull,0)) {return (-2);}

        constexpr std::size_t max_buff_size{256ull};
        __ATTR_ALIGN__(32)
        __m256 ymm0[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm1[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm2[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm3[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm4[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm5[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm6[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm7[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm8[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm9[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm10[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm11[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm12[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm13[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm14[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm15[max_buff_size];
        float scal_rem[max_buff_size];
        __m256 vzero{_mm256_setzero_ps()};
        __m256 sum0{vzero};
        __m256 sum1{sum0};
        __m256 sum2{sum0};
        __m256 sum3{sum0};
        __m256 sum4{sum0};
        __m256 sum5{sum0};
        __m256 sum6{sum0};
        __m256 sum7{sum0};
        __m256 sum8{sum0};
        __m256 sum9{sum0};
        __m256 sum10{sum0};
        __m256 sum11{sum0};
        __m256 sum12{sum0};
        __m256 sum13{sum0};
        __m256 sum14{sum0};
        __m256 sum15{sum0};
        float sum_rem{0.0f};
        std::size_t i;
        // Create a base-band cosine signal components of this sequence.
       
        switch (which_generator)
        {
             case 0 :
             {
                  for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                  {
            
                            ret_codes[k] = this->m_bb_sine_signals.operator[](k).create_signal_user_data(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
                  }      
             }     
             break;
             case 1 :
             {
                  for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                  {
            
                            ret_codes[k] = this->m_bb_sine_signals.operator[](k).create_signal_user_data_u4x(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
                  } 
             }
             break;
             default : 
                  return (-3);     
      
        }

        for(i = 0ull; (i+127ull) < this->m_nsamples; i += 128ull) 
        {
           
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;
              sum4  = vzero;
              sum5  = vzero;
              sum6  = vzero;
              sum7  = vzero;
              sum8  = vzero;
              sum9  = vzero;
              sum10 = vzero;
              sum11 = vzero;
              sum12 = vzero;
              sum13 = vzero;
              sum14 = vzero;
              sum15 = vzero;              
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
#if (AM_BB_SINE_SEQUENCE_SOFT_PREFETCH) == 1
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+40ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+56ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+64ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+72ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+80ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+88ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+96ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+104ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+112ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+120ull],_MM_HINT_T1);
                                          
#endif 
                   
                    ymm0[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0[j],sum0);
                    ymm1[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1[j],sum1);
                    ymm2[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2[j],sum2);
                    ymm3[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3[j],sum3);
                    ymm4[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum4  =  _mm256_add_ps(ymm4[j],sum4);
                    ymm5[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+40ull]);
                    sum5  =  _mm256_add_ps(ymm5[j],sum5);
                    ymm6[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum6  =  _mm256_add_ps(ymm6[j],sum6); 
                    ymm7[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+56ull]);
                    sum7  =  _mm256_add_ps(ymm7[j],sum7);
                    ymm8[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+64ull]);
                    sum8  =  _mm256_add_ps(ymm8[j],sum8);
                    ymm9[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+72ull]);
                    sum9  =  _mm256_add_ps(ymm9[j],sum9);
                    ymm10[j] =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+80ull]);
                    sum10 =  _mm256_add_ps(ymm10[j],sum10);
                    ymm11[j] =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+88ull]);
                    sum11 =  _mm256_add_ps(ymm11[j],sum11);
                    ymm12[j] =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+96ull]);
                    sum12 =  _mm256_add_ps(ymm12[j],sum12);
                    ymm13[j] =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+104ull]);
                    sum13 =  _mm256_add_ps(ymm13[j],sum13);
                    ymm14[j] =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+112ull]);
                    sum14 =  _mm256_add_ps(ymm14[j],sum14);
                    ymm15[j]  = _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+120ull]);
                    sum15  = _mm256_add_ps(ymm15[j],sum15);
                   
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum2);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum3);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum4);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+40ull],sum5);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum6);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+56ull],sum7);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+64ull],sum8);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+72ull],sum9);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+80ull],sum10);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+88ull],sum11);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+96ull],sum12);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+104ull],sum13);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+112ull],sum14);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+120ull],sum15);
              }
                   
        }
       //////////////////////////////////////////////////////
        for(; (i+95ull) < this->m_nsamples; i += 96ull) 
        {
            
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;
              sum4  = vzero;
              sum5  = vzero;
              sum6  = vzero;
              sum7  = vzero;
              sum8  = vzero;
              sum9  = vzero;
              sum10 = vzero;
              sum11 = vzero;
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
                   
                    ymm0[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0[j],sum0);
                    ymm1[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1[j],sum1);
                    ymm2[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2[j],sum2);
                    ymm3[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3[j],sum3);
                    ymm4[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum4  =  _mm256_add_ps(ymm4[j],sum4);
                    ymm5[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+40ull]);
                    sum5  =  _mm256_add_ps(ymm5[j],sum5);
                    ymm6[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum6  =  _mm256_add_ps(ymm6[j],sum6); 
                    ymm7[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+56ull]);
                    sum7  =  _mm256_add_ps(ymm7[j],sum7);
                    ymm8[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+64ull]);
                    sum8  =  _mm256_add_ps(ymm8[j],sum8);
                    ymm9[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+72ull]);
                    sum9  =  _mm256_add_ps(ymm9[j],sum9);
                    ymm10[j] =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+80ull]);
                    sum10 =  _mm256_add_ps(ymm10[j],sum10);
                    ymm11[j] =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+88ull]);
                    sum11 =  _mm256_add_ps(ymm11[j],sum11);
                    
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum2);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum3);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum4);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+40ull],sum5);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum6);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+56ull],sum7);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+64ull],sum8);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+72ull],sum9);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+80ull],sum10);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+88ull],sum11);
              }
                           
        }
        //////////////////////////////////////////////////////////////////
        for(; (i+79ull) < this->m_nsamples; i += 80ull) 
        {
            
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;
              sum4  = vzero;
              sum5  = vzero;
              sum6  = vzero;
              sum7  = vzero;
              sum8  = vzero;
              sum9  = vzero;       
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                    ymm0[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0[j],sum0);
                    ymm1[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1[j],sum1);
                    ymm2[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2[j],sum2);
                    ymm3[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3[j],sum3);
                    ymm4[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum4  =  _mm256_add_ps(ymm4[j],sum4);
                    ymm5[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+40ull]);
                    sum5  =  _mm256_add_ps(ymm5[j],sum5);
                    ymm6[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum6  =  _mm256_add_ps(ymm6[j],sum6); 
                    ymm7[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+56ull]);
                    sum7  =  _mm256_add_ps(ymm7[j],sum7);
                    ymm8[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+64ull]);
                    sum8  =  _mm256_add_ps(ymm8[j],sum8);
                    ymm9[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+72ull]);
                    sum9  =  _mm256_add_ps(ymm9[j],sum9);
                   
                   _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                   _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);
                   _mm256_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum2);
                   _mm256_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum3);
                   _mm256_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum4);
                   _mm256_storeu_ps(&this->m_vsummed_sequence[i+40ull],sum5);
                   _mm256_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum6);
                   _mm256_storeu_ps(&this->m_vsummed_sequence[i+56ull],sum7);
                   _mm256_storeu_ps(&this->m_vsummed_sequence[i+64ull],sum8);
                   _mm256_storeu_ps(&this->m_vsummed_sequence[i+72ull],sum9); 
                            
              }
                                             
        }
        ///////////////////////////////////////////////////////////////
        for(; (i+63ull) < this->m_nsamples; i += 64ull) 
        {
            
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;
              sum4  = vzero;
              sum5  = vzero;
              sum6  = vzero;
              sum7  = vzero;                
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                    ymm0[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0[j],sum0);
                    ymm1[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1[j],sum1);
                    ymm2[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2[j],sum2);
                    ymm3[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3[j],sum3);
                    ymm4[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum4  =  _mm256_add_ps(ymm4[j],sum4);
                    ymm5[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+40ull]);
                    sum5  =  _mm256_add_ps(ymm5[j],sum5);
                    ymm6[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum6  =  _mm256_add_ps(ymm6[j],sum6); 
                    ymm7[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+56ull]);
                    sum7  =  _mm256_add_ps(ymm7[j],sum7);                 

                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum2);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum3);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum4);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+40ull],sum5);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum6);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+56ull],sum7);                            
              }              
                                            
        }
        /////////////////////////////////////////////////////////////////
        for(; (i+31ull) < this->m_nsamples; i += 32ull) 
        {
            
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;       
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
                    ymm0[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0[j],sum0);
                    ymm1[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1[j],sum1);
                    ymm2[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2[j],sum2);
                    ymm3[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3[j],sum3);

                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum2);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum3);                                                
              }
                                                            
        }
        //////////////////////////////////////////////////////////////////
        for(; (i+15ull) < this->m_nsamples; i += 16ull) 
        {            
              sum0  = vzero;
              sum1  = vzero;                  
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
                    ymm0[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0[j],sum0);
                    ymm1[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1[j],sum1);  

                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);                                                                          
              }
                                                            
        }
        /////////////////////////////////////////////////////////////////
        for(; (i+7ull) < this->m_nsamples; i += 8ull) 
        {
               sum0  = vzero;                           
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
                    ymm0[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0[j],sum0);

                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);                                                                    
              }
                                                                           
        }
        ///////////////////////////////////////////////////////////////
        for(; (i+0ull) < this->m_nsamples; i += 1ull) 
        {
            sum_rem = 0.0f;
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {
               scal_rem[j]  =  this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i];
               sum_rem      +=  scal_rem[j];
               this->m_vsummed_sequence[i] = sum_rem;                                             
            }
           
        }
        return (0);
}

std::int32_t 
gms::radiolocation 
::am_bb_sine_signal_sequence_t
::signal_sequence_avx_u10x(const float * __restrict__ sym_in ,
                             std::int32_t * __restrict__ ret_codes,
                             const std::uint32_t which_generator) 
{
        if(__builtin_expect(nullptr==sym_in,0)) {return (-1);}
        if(__builtin_expect(this->m_nsignals>256ull,0)) {return (-2);}

        constexpr std::size_t max_buff_size{256ull};
        __ATTR_ALIGN__(32)
        __m256 ymm0[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm1[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm2[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm3[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm4[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm5[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm6[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm7[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm8[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm9[max_buff_size];
        float scal_rem[max_buff_size];
        __m256 vzero{_mm256_setzero_ps()};
        __m256 sum0{vzero};
        __m256 sum1{sum0};
        __m256 sum2{sum0};
        __m256 sum3{sum0};
        __m256 sum4{sum0};
        __m256 sum5{sum0};
        __m256 sum6{sum0};
        __m256 sum7{sum0};
        __m256 sum8{sum0};
        __m256 sum9{sum0};
        float sum_rem{0.0f};
       
        std::size_t i;
        // Create a base-band cosine signal components of this sequence.
        switch (which_generator)
        {
             case 0 :
             {
                  for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                  {
            
                            ret_codes[k] = this->m_bb_sine_signals.operator[](k).create_signal_user_data(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
                  }      
             }     
             break;
             case 1 :
             {
                  for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                  {
            
                            ret_codes[k] = this->m_bb_sine_signals.operator[](k).create_signal_user_data_u4x(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
                  } 
             }
             break;
             default : 
                  return (-3);     
         }

        for(i = 0ull; (i+79ull) < this->m_nsamples; i += 80ull) 
        {
              sum0 = vzero;
              sum1 = vzero;
              sum2 = vzero;
              sum3 = vzero;
              sum4 = vzero;
              sum5 = vzero;
              sum6 = vzero;
              sum7 = vzero;
              sum8 = vzero;
              sum9 = vzero;
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
#if (AM_BB_SINE_SEQUENCE_SOFT_PREFETCH) == 1
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+40ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+56ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+64ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+72ull],_MM_HINT_T1);
                                                         
#endif 
                   
                    ymm0[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0[j],sum0);
                    ymm1[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1[j],sum1);
                    ymm2[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2[j],sum2);
                    ymm3[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3[j],sum3);
                    ymm4[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum4  =  _mm256_add_ps(ymm4[j],sum4);
                    ymm5[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+40ull]);
                    sum5  =  _mm256_add_ps(ymm5[j],sum5);
                    ymm6[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum6  =  _mm256_add_ps(ymm6[j],sum6); 
                    ymm7[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+56ull]);
                    sum7  =  _mm256_add_ps(ymm7[j],sum7);
                    ymm8[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+64ull]);
                    sum8  =  _mm256_add_ps(ymm8[j],sum8);
                    ymm9[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+72ull]);
                    sum9  =  _mm256_add_ps(ymm9[j],sum9);

                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum2);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum3);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum4);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+40ull],sum5);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum6);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+56ull],sum7);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+64ull],sum8);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+72ull],sum9);
              }

              
        }
           
        ///////////////////////////////////////////////////////////////
        for(; (i+63ull) < this->m_nsamples; i += 64ull) 
        {
              
              sum0 = vzero;
              sum1 = vzero;
              sum2 = vzero;
              sum3 = vzero;
              sum4 = vzero;
              sum5 = vzero;
              sum6 = vzero;
              sum7 = vzero;                   
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                    ymm0[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0[j],sum0);
                    ymm1[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1[j],sum1);
                    ymm2[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2[j],sum2);
                    ymm3[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3[j],sum3);
                    ymm4[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum4  =  _mm256_add_ps(ymm4[j],sum4);
                    ymm5[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+40ull]);
                    sum5  =  _mm256_add_ps(ymm5[j],sum5);
                    ymm6[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum6  =  _mm256_add_ps(ymm6[j],sum6); 
                    ymm7[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+56ull]);
                    sum7  =  _mm256_add_ps(ymm7[j],sum7);
                   
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum2);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum3);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum4);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+40ull],sum5);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum6);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+56ull],sum7);                                   
              }

               
                                            
        }
        /////////////////////////////////////////////////////////////////
        for(; (i+31ull) < this->m_nsamples; i += 32ull) 
        {
           
              sum0 = vzero;
              sum1 = vzero;
              sum2 = vzero;
              sum3 = vzero;       
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                    ymm0[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0[j],sum0);
                    ymm1[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1[j],sum1);
                    ymm2[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2[j],sum2);
                    ymm3[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3[j],sum3);

                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum2);
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum3);                                           
              }              
                                             
        }
        //////////////////////////////////////////////////////////////////
        for(; (i+15ull) < this->m_nsamples; i += 16ull) 
        {
            
              sum0 = vzero;
              sum1 = vzero;                  
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                    ymm0[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0[j],sum0);
                    ymm1[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1[j],sum1);
                   
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                     _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);                                                       
              }
              
                                                            
        }
        /////////////////////////////////////////////////////////////////
        for(; (i+7ull) < this->m_nsamples; i += 8ull) 
        {
              sum0 = vzero;                              
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                    ymm0[j]  =  _mm256_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0[j],sum0);  
                    _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);                                                                                        
              }              
                                                           
        }
        ///////////////////////////////////////////////////////////////
        for(; (i+0ull) < this->m_nsamples; i += 1ull) 
        {
            sum_rem = 0.0f;
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               scal_rem[j] =  this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i];
               sum_rem     += scal_rem[j];
               this->m_vsummed_sequence[i] = sum_rem;
                                                              
            }
           
        }
        return (0);
}


std::int32_t 
gms::radiolocation 
::am_bb_sine_signal_sequence_t
::signal_sequence_avx_u6x(const float * __restrict__ sym_in ,
                             std::int32_t * __restrict__ ret_codes,
                             const std::uint32_t which_generator)  
{
        if(__builtin_expect(nullptr==sym_in,0)) {return (-1);}
        if(__builtin_expect(this->m_nsignals>256ull,0)) {return (-2);}

        constexpr std::size_t max_buff_size{256ull};
        __ATTR_ALIGN__(32)
        __m256 ymm0[max_buff_size];
        __ATTR_ALIGN__(32)
        __m256 ymm1[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm2[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm3[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm4[max_buff_size];
         __ATTR_ALIGN__(32)
        __m256 ymm5[max_buff_size];
        float scal_rem[max_buff_size];
        __m256 vzero{_mm256_setzero_ps()};
        __m256 sum0{vzero};
        __m256 sum1{sum0};
        __m256 sum2{sum0};
        __m256 sum3{sum0};
        __m256 sum4{sum0};
        __m256 sum5{sum0};
        
        float sum_rem{0.0f};

        std::size_t i;
        // Create a base-band cosine signal components of this sequence.
        switch (which_generator)
        {
             case 0 :
             {
                  for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                  {
            
                            ret_codes[k] = this->m_bb_sine_signals.operator[](k).create_signal_user_data(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
                  }      
             }     
             break;
             case 1 :
             {
                  for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                  {
            
                            ret_codes[k] = this->m_bb_sine_signals.operator[](k).create_signal_user_data_u4x(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
                  } 
             }
             break;
             default : 
                  return (-3);     
         }

        for(i = 0ull; (i+47ull) < this->m_nsamples; i += 48ull) 
        {
            
              sum0 = vzero;
              sum1 = vzero;
              sum2 = vzero;
              sum3 = vzero;
              sum4 = vzero; 
              sum5 = vzero;     
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
#if (AM_BB_SINE_SEQUENCE_SOFT_PREFETCH) == 1
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+40ull],_MM_HINT_T1);
                                                                        
#endif 
                    
                    ymm0[j]  =  _mm256_load_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0[j],sum0);
                    ymm1[j]  =  _mm256_load_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1[j],sum1);
                    ymm2[j]  =  _mm256_load_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2[j],sum2);
                    ymm3[j]  =  _mm256_load_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3[j],sum3);
                    ymm4[j]  =  _mm256_load_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum4  =  _mm256_add_ps(ymm4[j],sum4);
                    ymm5[j]  =  _mm256_load_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+40ull]);
                    sum5  =  _mm256_add_ps(ymm5[j],sum5);
                       
               _mm256_store_ps(&this->m_vsummed_sequence[i+0ull], sum0);
               _mm256_store_ps(&this->m_vsummed_sequence[i+8ull], sum1);
               _mm256_store_ps(&this->m_vsummed_sequence[i+16ull],sum2);
               _mm256_store_ps(&this->m_vsummed_sequence[i+24ull],sum3);
               _mm256_store_ps(&this->m_vsummed_sequence[i+32ull],sum4);
               _mm256_store_ps(&this->m_vsummed_sequence[i+40ull],sum5);
                  
              }
                          
        }
           
        /////////////////////////////////////////////////////////////////
        for(; (i+31ull) < this->m_nsamples; i += 32ull) 
        {
                        
              
              sum0 = vzero;
              sum1 = vzero;
              sum2 = vzero;
              sum3 = vzero;      
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
                                     
                    ymm0[j]  =  _mm256_load_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0[j],sum0);
                    ymm1[j]  =  _mm256_load_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1[j],sum1);
                    ymm2[j]  =  _mm256_load_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2[j],sum2);
                    ymm3[j]  =  _mm256_load_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3[j],sum3);

                     _mm256_store_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                     _mm256_store_ps(&this->m_vsummed_sequence[i+8ull], sum1);
                     _mm256_store_ps(&this->m_vsummed_sequence[i+16ull],sum2);
                     _mm256_store_ps(&this->m_vsummed_sequence[i+24ull],sum3);                               
              }

              
                                             
        }
        //////////////////////////////////////////////////////////////////
        for(; (i+15ull) < this->m_nsamples; i += 16ull) 
        {
                         
              sum0 = vzero;
              sum1 = vzero;             
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
                                                                            
                    ymm0[j]  =  _mm256_load_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0[j],sum0);
                    ymm1[j]  =  _mm256_load_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1[j],sum1);

                     _mm256_store_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                     _mm256_store_ps(&this->m_vsummed_sequence[i+8ull], sum1);                                                       
              }
                                                                         
        }
        /////////////////////////////////////////////////////////////////
        for(; (i+7ull) < this->m_nsamples; i += 8ull) 
        {
              sum0 = vzero;                         
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
                   
                    ymm0[j]  =  _mm256_load_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0[j],sum0);

                     _mm256_store_ps(&this->m_vsummed_sequence[i+0ull], sum0);                                                                     
              }              
                                                           
        }
        ///////////////////////////////////////////////////////////////
        for(; (i+0ull) < this->m_nsamples; i += 1ull) 
        {
            sum_rem = 0.0f;
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {
               
               scal_rem[j] = this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i];
               sum_rem  +=  scal_rem[j];
               this->m_vsummed_sequence[i] = sum_rem;                                             
            }
            
        }
        return (0);
}

std::int32_t
gms::radiolocation 
::am_bb_sine_signal_sequence_t
::signal_sequence_avx512_u16x(const float * __restrict__ sym_in ,
                                std::int32_t * __restrict__ ret_codes,
                                const std::uint32_t which_generator)  
{
        if(__builtin_expect(nullptr==sym_in,0)) {return (-1);}
        if(__builtin_expect(this->m_nsignals>256ull,0)) {return (-2);}

        constexpr std::size_t max_buff_size{256ull};
        __ATTR_ALIGN__(64)
        __m512 zmm0[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm1[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm2[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm3[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm4[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm5[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm6[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm7[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm8[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm9[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm10[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm11[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm12[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm13[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm14[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm15[max_buff_size];
        float scal_rem[max_buff_size];
        __m512 vzero{_mm512_setzero_ps()};
        __m512 sum0{vzero};
        __m512 sum1{sum0};
        __m512 sum2{sum0};
        __m512 sum3{sum0};
        __m512 sum4{sum0};
        __m512 sum5{sum0};
        __m512 sum6{sum0};
        __m512 sum7{sum0};
        __m512 sum8{sum0};
        __m512 sum9{sum0};
        __m512 sum10{sum0};
        __m512 sum11{sum0};
        __m512 sum12{sum0};
        __m512 sum13{sum0};
        __m512 sum14{sum0};
        __m512 sum15{sum0};
        float sum_rem{0.0f};

        std::size_t i;
        // Create a base-band cosine signal components of this sequence.
        switch (which_generator)
        {
             case 0 :
             {
                  for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                  {
            
                            ret_codes[k] = this->m_bb_sine_signals.operator[](k).create_signal_user_data(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
                  }      
             }     
             break;
             case 1 :
             {
                  for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                  {
            
                            ret_codes[k] = this->m_bb_sine_signals.operator[](k).create_signal_user_data_u4x(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
                  } 
             }
             break;
             default : 
                  return (-3);     
         }

        for(i = 0ull; (i+255ull) < this->m_nsamples; i += 256ull) 
        {
           
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;
              sum4  = vzero;
              sum5  = vzero;
              sum6  = vzero;
              sum7  = vzero;
              sum8  = vzero;
              sum9  = vzero;
              sum10 = vzero;
              sum11 = vzero;
              sum12 = vzero;
              sum13 = vzero;
              sum14 = vzero;
              sum15 = vzero;              
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
#if (AM_BB_SINE_SEQUENCE_SOFT_PREFETCH) == 1
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+64ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+80ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+96ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+112ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+128ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+144ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+160ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+176ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+192ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+208ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+224ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+240ull],_MM_HINT_T1);
                                          
#endif 
                   
                    zmm0[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm512_add_ps(zmm0[j],sum0);
                    zmm1[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum1  =  _mm512_add_ps(zmm1[j],sum1);
                    zmm2[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum2  =  _mm512_add_ps(zmm2[j],sum2);
                    zmm3[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum3  =  _mm512_add_ps(zmm3[j],sum3);
                    zmm4[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+64ull]);
                    sum4  =  _mm512_add_ps(zmm4[j],sum4);
                    zmm5[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+80ull]);
                    sum5  =  _mm512_add_ps(zmm5[j],sum5);
                    zmm6[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+96ull]);
                    sum6  =  _mm512_add_ps(zmm6[j],sum6); 
                    zmm7[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+112ull]);
                    sum7  =  _mm512_add_ps(zmm7[j],sum7);
                    zmm8[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+128ull]);
                    sum8  =  _mm512_add_ps(zmm8[j],sum8);
                    zmm9[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+144ull]);
                    sum9  =  _mm512_add_ps(zmm9[j],sum9);
                    zmm10[j] =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+160ull]);
                    sum10 =  _mm512_add_ps(zmm10[j],sum10);
                    zmm11[j] =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+176ull]);
                    sum11 =  _mm512_add_ps(zmm11[j],sum11);
                    zmm12[j] =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+192ull]);
                    sum12 =  _mm512_add_ps(zmm12[j],sum12);
                    zmm13[j] =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+208ull]);
                    sum13 =  _mm512_add_ps(zmm13[j],sum13);
                    zmm14[j] =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+224ull]);
                    sum14 =  _mm512_add_ps(zmm14[j],sum14);
                    zmm15[j]  = _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+240ull]);
                    sum15  = _mm512_add_ps(zmm15[j],sum15);
                   
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+16ull], sum1);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum2);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum3);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+64ull],sum4);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+80ull],sum5);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+96ull],sum6);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+112ull],sum7);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+128ull],sum8);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+144ull],sum9);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+160ull],sum10);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+176ull],sum11);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+192ull],sum12);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+208ull],sum13);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+224ull],sum14);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+240ull],sum15);
              }
                   
        }
       //////////////////////////////////////////////////////
        for(; (i+191ull) < this->m_nsamples; i += 192ull) 
        {
            
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;
              sum4  = vzero;
              sum5  = vzero;
              sum6  = vzero;
              sum7  = vzero;
              sum8  = vzero;
              sum9  = vzero;
              sum10 = vzero;
              sum11 = vzero;
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
                    zmm0[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm512_add_ps(zmm0[j],sum0);
                    zmm1[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum1  =  _mm512_add_ps(zmm1[j],sum1);
                    zmm2[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum2  =  _mm512_add_ps(zmm2[j],sum2);
                    zmm3[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum3  =  _mm512_add_ps(zmm3[j],sum3);
                    zmm4[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+64ull]);
                    sum4  =  _mm512_add_ps(zmm4[j],sum4);
                    zmm5[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+80ull]);
                    sum5  =  _mm512_add_ps(zmm5[j],sum5);
                    zmm6[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+96ull]);
                    sum6  =  _mm512_add_ps(zmm6[j],sum6); 
                    zmm7[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+112ull]);
                    sum7  =  _mm512_add_ps(zmm7[j],sum7);
                    zmm8[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+128ull]);
                    sum8  =  _mm512_add_ps(zmm8[j],sum8);
                    zmm9[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+144ull]);
                    sum9  =  _mm512_add_ps(zmm9[j],sum9);
                    zmm10[j] =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+160ull]);
                    sum10 =  _mm512_add_ps(zmm10[j],sum10);
                    zmm11[j] =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+176ull]);
                    sum11 =  _mm512_add_ps(zmm11[j],sum11);

                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+16ull], sum1);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum2);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum3);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+64ull],sum4);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+80ull],sum5);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+96ull],sum6);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+112ull],sum7);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+128ull],sum8);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+144ull],sum9);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+160ull],sum10);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+176ull],sum11);
                    
              }
                           
        }
        //////////////////////////////////////////////////////////////////
        for(; (i+159ull) < this->m_nsamples; i += 160ull) 
        {
            
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;
              sum4  = vzero;
              sum5  = vzero;
              sum6  = vzero;
              sum7  = vzero;
              sum8  = vzero;
              sum9  = vzero;       
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                    zmm0[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm512_add_ps(zmm0[j],sum0);
                    zmm1[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum1  =  _mm512_add_ps(zmm1[j],sum1);
                    zmm2[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum2  =  _mm512_add_ps(zmm2[j],sum2);
                    zmm3[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum3  =  _mm512_add_ps(zmm3[j],sum3);
                    zmm4[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+64ull]);
                    sum4  =  _mm512_add_ps(zmm4[j],sum4);
                    zmm5[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+80ull]);
                    sum5  =  _mm512_add_ps(zmm5[j],sum5);
                    zmm6[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+96ull]);
                    sum6  =  _mm512_add_ps(zmm6[j],sum6); 
                    zmm7[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+112ull]);
                    sum7  =  _mm512_add_ps(zmm7[j],sum7);
                    zmm8[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+128ull]);
                    sum8  =  _mm512_add_ps(zmm8[j],sum8);
                    zmm9[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+144ull]);
                    sum9  =  _mm512_add_ps(zmm9[j],sum9);
                                     
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+16ull], sum1);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum2);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum3);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+64ull],sum4);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+80ull],sum5);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+96ull],sum6);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+112ull],sum7);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+128ull],sum8);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+144ull],sum9);        
              }
                                             
        }
        ///////////////////////////////////////////////////////////////
        for(; (i+127ull) < this->m_nsamples; i += 128ull) 
        {
            
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;
              sum4  = vzero;
              sum5  = vzero;
              sum6  = vzero;
              sum7  = vzero;                
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                    zmm0[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm512_add_ps(zmm0[j],sum0);
                    zmm1[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum1  =  _mm512_add_ps(zmm1[j],sum1);
                    zmm2[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum2  =  _mm512_add_ps(zmm2[j],sum2);
                    zmm3[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum3  =  _mm512_add_ps(zmm3[j],sum3);
                    zmm4[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+64ull]);
                    sum4  =  _mm512_add_ps(zmm4[j],sum4);
                    zmm5[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+80ull]);
                    sum5  =  _mm512_add_ps(zmm5[j],sum5);
                    zmm6[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+96ull]);
                    sum6  =  _mm512_add_ps(zmm6[j],sum6); 
                    zmm7[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+112ull]);
                    sum7  =  _mm512_add_ps(zmm7[j],sum7); 

                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+16ull], sum1);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum2);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum3);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+64ull],sum4);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+80ull],sum5);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+96ull],sum6);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+112ull],sum7);                      
              }              
                                            
        }
        /////////////////////////////////////////////////////////////////
        for(; (i+63ull) < this->m_nsamples; i += 64ull) 
        {
            
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;       
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
                    zmm0[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm512_add_ps(zmm0[j],sum0);
                    zmm1[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum1  =  _mm512_add_ps(zmm1[j],sum1);
                    zmm2[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum2  =  _mm512_add_ps(zmm2[j],sum2);
                    zmm3[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum3  =  _mm512_add_ps(zmm3[j],sum3);   

                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+16ull], sum1);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum2);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum3);                                          
              }
                                                            
        }
        //////////////////////////////////////////////////////////////////
        for(; (i+31ull) < this->m_nsamples; i += 32ull) 
        {            
              sum0  = vzero;
              sum1  = vzero;                  
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
                    zmm0[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm512_add_ps(zmm0[j],sum0);
                    zmm1[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum1  =  _mm512_add_ps(zmm1[j],sum1);  

                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+16ull], sum1);                                                                     
              }
                                                            
        }
        /////////////////////////////////////////////////////////////////
        for(; (i+15ull) < this->m_nsamples; i += 16ull) 
        {
               sum0  = vzero;                           
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
                    zmm0[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm512_add_ps(zmm0[j],sum0);

                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);                                                                  
              }                                                                           
        }
        ///////////////////////////////////////////////////////////////
        for(; (i+0ull) < this->m_nsamples; i += 1ull) 
        {
            sum_rem = 0.0f;
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {
               scal_rem[j]  =  this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i];
               sum_rem      +=  scal_rem[j];
               this->m_vsummed_sequence[i] = sum_rem;                                             
            }
           
        }
        return (0);
}

std::int32_t
gms::radiolocation 
::am_bb_sine_signal_sequence_t
::signal_sequence_avx512_u10x(const float * __restrict__ sym_in ,
                              std::int32_t * __restrict__ ret_codes,
                              const std::uint32_t which_generator) 
{
        if(__builtin_expect(nullptr==sym_in,0)) {return (-1);}
        if(__builtin_expect(this->m_nsignals>256ull,0)) {return (-2);}

        constexpr std::size_t max_buff_size{256ull};
        __ATTR_ALIGN__(64)
        __m512 zmm0[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm1[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm2[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm3[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm4[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm5[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm6[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm7[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm8[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm9[max_buff_size];
        float scal_rem[max_buff_size];
        __m512 vzero{_mm512_setzero_ps()};
        __m512 sum0{vzero};
        __m512 sum1{sum0};
        __m512 sum2{sum0};
        __m512 sum3{sum0};
        __m512 sum4{sum0};
        __m512 sum5{sum0};
        __m512 sum6{sum0};
        __m512 sum7{sum0};
        __m512 sum8{sum0};
        __m512 sum9{sum0};
        float sum_rem{0.0f};

        std::size_t i;
        // Create a base-band cosine signal components of this sequence.
        switch (which_generator)
        {
             case 0 :
             {
                  for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                  {
            
                            ret_codes[k] = this->m_bb_sine_signals.operator[](k).create_signal_user_data(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
                  }      
             }     
             break;
             case 1 :
             {
                  for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                  {
            
                            ret_codes[k] = this->m_bb_sine_signals.operator[](k).create_signal_user_data_u4x(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
                  } 
             }
             break;
             default : 
                  return (-3);     
         }

        for(i = 0ull; (i+159ull) < this->m_nsamples; i += 160ull) 
        {
           
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;
              sum4  = vzero;
              sum5  = vzero;
              sum6  = vzero;
              sum7  = vzero;
              sum8  = vzero;
              sum9  = vzero;
                        
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
#if (AM_BB_SINE_SEQUENCE_SOFT_PREFETCH) == 1
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+64ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+80ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+96ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+112ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+128ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+144ull],_MM_HINT_T1);
                                                     
#endif 
                   
                    zmm0[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm512_add_ps(zmm0[j],sum0);
                    zmm1[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum1  =  _mm512_add_ps(zmm1[j],sum1);
                    zmm2[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum2  =  _mm512_add_ps(zmm2[j],sum2);
                    zmm3[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum3  =  _mm512_add_ps(zmm3[j],sum3);
                    zmm4[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+64ull]);
                    sum4  =  _mm512_add_ps(zmm4[j],sum4);
                    zmm5[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+80ull]);
                    sum5  =  _mm512_add_ps(zmm5[j],sum5);
                    zmm6[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+96ull]);
                    sum6  =  _mm512_add_ps(zmm6[j],sum6); 
                    zmm7[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+112ull]);
                    sum7  =  _mm512_add_ps(zmm7[j],sum7);
                    zmm8[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+128ull]);
                    sum8  =  _mm512_add_ps(zmm8[j],sum8);
                    zmm9[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+144ull]);
                    sum9  =  _mm512_add_ps(zmm9[j],sum9);
                                      
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+16ull], sum1);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum2);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum3);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+64ull],sum4);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+80ull],sum5);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+96ull],sum6);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+112ull],sum7);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+128ull],sum8);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+144ull],sum9);
            }
                   
        }
             
        ///////////////////////////////////////////////////////////////
        for(; (i+127ull) < this->m_nsamples; i += 128ull) 
        {
            
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;
              sum4  = vzero;
              sum5  = vzero;
              sum6  = vzero;
              sum7  = vzero;                
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                    zmm0[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm512_add_ps(zmm0[j],sum0);
                    zmm1[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum1  =  _mm512_add_ps(zmm1[j],sum1);
                    zmm2[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum2  =  _mm512_add_ps(zmm2[j],sum2);
                    zmm3[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum3  =  _mm512_add_ps(zmm3[j],sum3);
                    zmm4[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+64ull]);
                    sum4  =  _mm512_add_ps(zmm4[j],sum4);
                    zmm5[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+80ull]);
                    sum5  =  _mm512_add_ps(zmm5[j],sum5);
                    zmm6[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+96ull]);
                    sum6  =  _mm512_add_ps(zmm6[j],sum6); 
                    zmm7[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+112ull]);
                    sum7  =  _mm512_add_ps(zmm7[j],sum7); 

                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+16ull], sum1);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum2);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum3);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+64ull],sum4);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+80ull],sum5);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+96ull],sum6);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+112ull],sum7);                      
              }              
                                            
        }
        /////////////////////////////////////////////////////////////////
        for(; (i+63ull) < this->m_nsamples; i += 64ull) 
        {
            
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;       
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
                    zmm0[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm512_add_ps(zmm0[j],sum0);
                    zmm1[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum1  =  _mm512_add_ps(zmm1[j],sum1);
                    zmm2[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum2  =  _mm512_add_ps(zmm2[j],sum2);
                    zmm3[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum3  =  _mm512_add_ps(zmm3[j],sum3);   

                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+16ull], sum1);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum2);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum3);                                          
              }
                                                            
        }
        //////////////////////////////////////////////////////////////////
        for(; (i+31ull) < this->m_nsamples; i += 32ull) 
        {            
              sum0  = vzero;
              sum1  = vzero;                  
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
                    zmm0[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm512_add_ps(zmm0[j],sum0);
                    zmm1[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum1  =  _mm512_add_ps(zmm1[j],sum1);  

                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+16ull], sum1);                                                                     
              }
                                                            
        }
        /////////////////////////////////////////////////////////////////
        for(; (i+15ull) < this->m_nsamples; i += 16ull) 
        {
               sum0  = vzero;                           
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
                    zmm0[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm512_add_ps(zmm0[j],sum0);

                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);                                                                  
              }                                                                           
        }
        ///////////////////////////////////////////////////////////////
        for(; (i+0ull) < this->m_nsamples; i += 1ull) 
        {
            sum_rem = 0.0f;
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {
               scal_rem[j]  =  this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i];
               sum_rem      +=  scal_rem[j];
               this->m_vsummed_sequence[i] = sum_rem;                                             
            }
           
        }
        return (0);
}

std::int32_t
gms::radiolocation 
::am_bb_sine_signal_sequence_t
::signal_sequence_avx512_u6x(const float * __restrict__ sym_in ,
                              std::int32_t * __restrict__ ret_codes,
                              const std::uint32_t which_generator) 
{
        if(__builtin_expect(nullptr==sym_in,0)) {return (-1);}
        if(__builtin_expect(this->m_nsignals>256ull,0)) {return (-2);}

        constexpr std::size_t max_buff_size{256ull};
        __ATTR_ALIGN__(64)
        __m512 zmm0[max_buff_size];
        __ATTR_ALIGN__(64)
        __m512 zmm1[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm2[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm3[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm4[max_buff_size];
         __ATTR_ALIGN__(64)
        __m512 zmm5[max_buff_size];
        float scal_rem[max_buff_size];
        __m512 vzero{_mm512_setzero_ps()};
        __m512 sum0{vzero};
        __m512 sum1{sum0};
        __m512 sum2{sum0};
        __m512 sum3{sum0};
        __m512 sum4{sum0};
        __m512 sum5{sum0};
        float sum_rem{0.0f};

        std::size_t i;
        // Create a base-band cosine signal components of this sequence.        
        switch (which_generator)
        {
             case 0 :
             {
                  for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                  {
            
                            ret_codes[k] = this->m_bb_sine_signals.operator[](k).create_signal_user_data(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
                  }      
             }     
             break;
             case 1 :
             {
                  for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
                  {
            
                            ret_codes[k] = this->m_bb_sine_signals.operator[](k).create_signal_user_data_u4x(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
                  } 
             }
             break;
             default : 
                  return (-3);     
         }

        for(i = 0ull; (i+95ull) < this->m_nsamples; i += 96ull) 
        {
           
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;
              sum4  = vzero;
              sum5  = vzero;
                                      
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
#if (AM_BB_SINE_SEQUENCE_SOFT_PREFETCH) == 1
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+64ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+80ull],_MM_HINT_T1);
                                                                   
#endif 
                   
                    zmm0[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm512_add_ps(zmm0[j],sum0);
                    zmm1[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum1  =  _mm512_add_ps(zmm1[j],sum1);
                    zmm2[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum2  =  _mm512_add_ps(zmm2[j],sum2);
                    zmm3[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum3  =  _mm512_add_ps(zmm3[j],sum3);
                    zmm4[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+64ull]);
                    sum4  =  _mm512_add_ps(zmm4[j],sum4);
                    zmm5[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+80ull]);
                    sum5  =  _mm512_add_ps(zmm5[j],sum5);
                                                         
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+16ull], sum1);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum2);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum3);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+64ull],sum4);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+80ull],sum5);
            }
                   
        }
           
        /////////////////////////////////////////////////////////////////
        for(; (i+63ull) < this->m_nsamples; i += 64ull) 
        {
            
              sum0  = vzero;
              sum1  = vzero;
              sum2  = vzero;
              sum3  = vzero;       
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
                    zmm0[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm512_add_ps(zmm0[j],sum0);
                    zmm1[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum1  =  _mm512_add_ps(zmm1[j],sum1);
                    zmm2[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum2  =  _mm512_add_ps(zmm2[j],sum2);
                    zmm3[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum3  =  _mm512_add_ps(zmm3[j],sum3);   

                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+16ull], sum1);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum2);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum3);                                          
              }
                                                            
        }
        //////////////////////////////////////////////////////////////////
        for(; (i+31ull) < this->m_nsamples; i += 32ull) 
        {            
              sum0  = vzero;
              sum1  = vzero;                  
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
                    zmm0[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm512_add_ps(zmm0[j],sum0);
                    zmm1[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum1  =  _mm512_add_ps(zmm1[j],sum1);  

                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+16ull], sum1);                                                                     
              }
                                                            
        }
        /////////////////////////////////////////////////////////////////
        for(; (i+15ull) < this->m_nsamples; i += 16ull) 
        {
               sum0  = vzero;                           
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
                    zmm0[j]  =  _mm512_loadu_ps(&this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm512_add_ps(zmm0[j],sum0);

                    _mm512_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);                                                                  
              }                                                                           
        }
        ///////////////////////////////////////////////////////////////
        for(; (i+0ull) < this->m_nsamples; i += 1ull) 
        {
            sum_rem = 0.0f;
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {
               scal_rem[j]  =  this->m_bb_sine_signals.operator[](j).m_sig_samples.m_data[i];
               sum_rem      +=  scal_rem[j];
               this->m_vsummed_sequence[i] = sum_rem;                                             
            }
           
        }
        return (0);
}





auto
gms::radiolocation
::operator<<(std::ostream &os,
             gms::radiolocation::am_bb_sine_signal_sequence_t &rhs)->std::ostream &
{
    std::cout << typeid(rhs).name() << "Begin: object state dump." << std::endl;
    std::cout << "Sequence of cosine signal-samples:" << std::endl;
    for(std::size_t __i{0ull}; __i != rhs.m_nsamples; ++__i)
    {
         os << std::fixed << std::setprecision(7) << 
                            rhs.m_vsummed_sequence.operator[](__i) << std::endl;
    }
    std::cout << "Number of narrowband cosine signals:" << std::endl;
    for(std::size_t __i{0ull}; __i != rhs.m_nsignals; ++__i) 
    {
         os << rhs.m_nK_values.operator[](__i) << std::endl;
    }
    std::cout << "Number of 2*PI cycles:" << std::endl;
    for(std::size_t __i{0ull}; __i != rhs.m_nsignals; ++__i) 
    {
        os << std::fixed << std::setprecision(7) << 
                           rhs.m_nK_values.operator[](__i) << std::endl;
    }
    std::cout << "Narrowband signals amplitude values:" << std::endl;
    for(std::size_t __i{0ull}; __i != rhs.m_nsignals; ++__i) 
    {
        os << std::fixed << std::setprecision(7) << 
                           rhs.m_A_values.operator[](__i) << std::endl;
    }
    std::cout << "Narrowband signals period values:" << std::endl;
    for(std::size_t __i{0ull}; __i != rhs.m_nsignals; ++__i) 
    {
        os << std::fixed << std::setprecision(7) << 
                           rhs.m_P_values.operator[](__i) << std::endl;
    }
    std::cout << typeid(rhs).name() << "End: object state dump." << std::endl;
    return (os);
}
