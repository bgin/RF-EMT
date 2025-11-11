
#include <immintrin.h>
#include <fstream>
#include <iomanip>
#include "GMS_am_bb_cosine_sequence.h"


gms::radiolocation
::am_bb_cos_signal_sequence_t
::am_bb_cos_signal_sequence_t(const std::size_t nsamples,
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
    this->m_bb_cos_signals   = std::vector<am_bb_cosine_signal_t>(this->m_nsignals);
    this->m_vsummed_sequence = std::valarray<float>(this->m_nsamples);
    for(std::size_t i{0ull}; i != this->m_nsignals; ++i) 
    {
        const std::uint32_t nK_vals{this->m_nK_values.operator[](i)};
        const float n_vals{this->m_n_values.operator[](i)};
        const float A_vals{this->m_A_values.operator[](i)};
        const float P_vals{this->m_P_values.operator[](i)};
        this->m_bb_cos_signals.operator[](i) = am_bb_cosine_signal_t(this->m_nsamples,
                                                                     nK_vals,n_vals,
                                                                     A_vals,P_vals);
    }
}

gms::radiolocation
::am_bb_cos_signal_sequence_t
::am_bb_cos_signal_sequence_t(am_bb_cos_signal_sequence_t &&other) noexcept(true) 
:
m_nsamples{std::move(other.m_nsamples)},
m_nsignals{std::move(other.m_nsignals)},
m_nK_values{std::move(other.m_nK_values)},
m_n_values{std::move(other.m_n_values)},
m_A_values{std::move(other.m_A_values)},
m_P_values{std::move(other.m_P_values)},
m_bb_cos_signals{std::move(other.m_bb_cos_signals)},
m_vsummed_sequence{std::move(other.m_vsummed_sequence)}
{

}

gms::radiolocation
::am_bb_cos_signal_sequence_t
::am_bb_cos_signal_sequence_t(const am_bb_cos_signal_sequence_t &other) noexcept(false) 
:
m_nsamples{other.m_nsamples},
m_nsignals{other.m_nsignals},
m_nK_values{other.m_nK_values},
m_n_values{other.m_n_values},
m_A_values{other.m_A_values},
m_P_values{other.m_P_values},
m_bb_cos_signals{other.m_bb_cos_signals},
m_vsummed_sequence{other.m_vsummed_sequence}
{

}

gms::radiolocation
::am_bb_cos_signal_sequence_t
::~am_bb_cos_signal_sequence_t() noexcept(false)
{

}

gms::radiolocation
::am_bb_cos_signal_sequence_t &
gms::radiolocation 
::am_bb_cos_signal_sequence_t
::operator=(const am_bb_cos_signal_sequence_t &other) noexcept(false) 
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_nsamples = other.m_nsamples;
    this->m_nsignals = other.m_nsignals;
    this->m_nK_values.operator=(other.m_nK_values);
    this->m_n_values.operator=(other.m_n_values);
    this->m_A_values.operator=(other.m_A_values);
    this->m_P_values.operator=(other.m_P_values);
    this->m_bb_cos_signals.operator=(other.m_bb_cos_signals);
    this->m_vsummed_sequence.operator=(other.m_vsummed_sequence);
    return (*this);
}

gms::radiolocation
::am_bb_cos_signal_sequence_t &
gms::radiolocation 
::am_bb_cos_signal_sequence_t
::operator=(am_bb_cos_signal_sequence_t &&other) noexcept(true) 
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_nsamples = std::move(other.m_nsamples);
    this->m_nsignals = std::move(other.m_nsignals);
    this->m_nK_values.operator=(std::move(other.m_nK_values));
    this->m_n_values.operator=(std::move(other.m_n_values));
    this->m_A_values.operator=(std::move(other.m_A_values));
    this->m_P_values.operator=(std::move(other.m_P_values));
    this->m_bb_cos_signals.operator=(std::move(other.m_bb_cos_signals));
    this->m_vsummed_sequence.operator=(std::move(other.m_vsummed_sequence));
    return (*this);
}

void 
gms::radiolocation
::am_bb_cos_signal_sequence_t
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

void 
gms::radiolocation
::am_bb_cos_signal_sequence_t
::signal_sequence_sse42_u16x_u(const float * __restrict__ sym_in ,
                               std::int32_t * __restrict__ ret_code) 
{
        if(__builtin_expect(nullptr==sym_in,0)) {return;}
        
        std::size_t i;
        __m128 sum0{_mm_setzero_ps()};
        __m128 sum1{_mm_setzero_ps()};
        __m128 sum2{_mm_setzero_ps()};
        __m128 sum3{_mm_setzero_ps()};
        __m128 sum4{_mm_setzero_ps()};
        __m128 sum5{_mm_setzero_ps()};
        __m128 sum6{_mm_setzero_ps()};
        __m128 sum7{_mm_setzero_ps()};
        __m128 sum8{_mm_setzero_ps()};
        __m128 sum9{_mm_setzero_ps()};
        __m128 sum10{_mm_setzero_ps()};
        __m128 sum11{_mm_setzero_ps()};
        __m128 sum12{_mm_setzero_ps()};
        __m128 sum13{_mm_setzero_ps()};
        __m128 sum14{_mm_setzero_ps()};
        __m128 sum15{_mm_setzero_ps()};
        float sum_rem{0.0f};

        for(std::size_t i{0ull}; i != this->m_nsignals; ++i) 
        {
            
            ret_code[i] = this->m_bb_cos_signals.operator[](i).create_signal_user_data(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](i));
           
        }

        for(i = 0ull; (i+63ull) < this->m_nsamples; i += 64ull) 
        {
            __m128 xmm0;
            __m128 xmm1;
            __m128 xmm2;
            __m128 xmm3;
            __m128 xmm4;
            __m128 xmm5;
            __m128 xmm6;
            __m128 xmm7;
            __m128 xmm8;
            __m128 xmm9;
            __m128 xmm10;
            __m128 xmm11;
            __m128 xmm12;
            __m128 xmm13;
            __m128 xmm14;
            __m128 xmm15;
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {
#if (AM_BB_COSINE_SEQUENCE_SOFT_PREFETCH) == 1
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+4ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+12ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+20ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+28ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+32ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+36ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+40ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+44ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+48ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+52ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+56ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+60ull],_MM_HINT_T1);
#endif 
               xmm0  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
               sum0  =  _mm_add_ps(xmm0,sum0);
               xmm1  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
               sum1  =  _mm_add_ps(xmm1,sum1);
               xmm2  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
               sum2  =  _mm_add_ps(xmm2,sum2);
               xmm3  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+12ull]);
               sum3  =  _mm_add_ps(xmm3,sum3);
               xmm4  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
               sum4  =  _mm_add_ps(xmm4,sum4);
               xmm5  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+20ull]);
               sum5  =  _mm_add_ps(xmm5,sum5);
               xmm6  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
               sum6  =  _mm_add_ps(xmm6,sum6);
               xmm7  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+28ull]);
               sum7  =  _mm_add_ps(xmm7,sum7);
               xmm8  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
               sum8  =  _mm_add_ps(xmm8,sum8);
               xmm9  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+36ull]);
               sum9  =  _mm_add_ps(xmm9,sum9);
               xmm10 =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+40ull]);
               sum10 =  _mm_add_ps(xmm10,sum10);
               xmm11 =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+44ull]);
               sum11 =  _mm_add_ps(xmm11,sum11);
               xmm12 =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
               sum12 =  _mm_add_ps(xmm12,sum12);
               xmm13 =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+52ull]);
               sum13 =  _mm_add_ps(xmm13,sum13);
               xmm14 =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+56ull]);
               sum14 =  _mm_add_ps(xmm14,sum14);
               xmm15 =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+60ull]);
               sum15 =  _mm_add_ps(xmm15,sum15);
               
            }
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

        for(; (i+39ull) < this->m_nsamples; i += 40ull) 
        {
            __m128 xmm0;
            __m128 xmm1;
            __m128 xmm2;
            __m128 xmm3;
            __m128 xmm4;
            __m128 xmm5;
            __m128 xmm6;
            __m128 xmm7;
            __m128 xmm8;
            __m128 xmm9;
           
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               xmm0  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
               sum0  =  _mm_add_ps(xmm0,sum0);
               xmm1  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
               sum1  =  _mm_add_ps(xmm1,sum1);
               xmm2  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
               sum2  =  _mm_add_ps(xmm2,sum2);
               xmm3  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+12ull]);
               sum3  =  _mm_add_ps(xmm3,sum3);
               xmm4  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
               sum4  =  _mm_add_ps(xmm4,sum4);
               xmm5  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+20ull]);
               sum5  =  _mm_add_ps(xmm5,sum5);
               xmm6  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
               sum6  =  _mm_add_ps(xmm6,sum6);
               xmm7  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+28ull]);
               sum7  =  _mm_add_ps(xmm7,sum7);
               xmm8  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
               sum8  =  _mm_add_ps(xmm8,sum8);
               xmm9  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+36ull]);
               sum9  =  _mm_add_ps(xmm9,sum9);
                          
            }
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

        for(; (i+31ull) < this->m_nsamples; i += 32ull) 
        {
            __m128 xmm0;
            __m128 xmm1;
            __m128 xmm2;
            __m128 xmm3;
            __m128 xmm4;
            __m128 xmm5;
            __m128 xmm6;
            __m128 xmm7;
                       
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               xmm0  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
               sum0  =  _mm_add_ps(xmm0,sum0);
               xmm1  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
               sum1  =  _mm_add_ps(xmm1,sum1);
               xmm2  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
               sum2  =  _mm_add_ps(xmm2,sum2);
               xmm3  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+12ull]);
               sum3  =  _mm_add_ps(xmm3,sum3);
               xmm4  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
               sum4  =  _mm_add_ps(xmm4,sum4);
               xmm5  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+20ull]);
               sum5  =  _mm_add_ps(xmm5,sum5);
               xmm6  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
               sum6  =  _mm_add_ps(xmm6,sum6);
               xmm7  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+28ull]);
               sum7  =  _mm_add_ps(xmm7,sum7);
                                        
            }
             _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum2);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+12ull],sum3);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum4);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+20ull],sum5);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum6);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+28ull],sum7);
        }

        for(; (i+15ull) < this->m_nsamples; i += 16ull) 
        {
            __m128 xmm0;
            __m128 xmm1;
            __m128 xmm2;
            __m128 xmm3;
                                   
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               xmm0  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
               sum0  =  _mm_add_ps(xmm0,sum0);
               xmm1  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
               sum1  =  _mm_add_ps(xmm1,sum1);
               xmm2  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
               sum2  =  _mm_add_ps(xmm2,sum2);
               xmm3  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+12ull]);
               sum3  =  _mm_add_ps(xmm3,sum3);
               
                                        
            }
             _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum2);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+12ull],sum3);
        }

        for(; (i+7ull) < this->m_nsamples; i += 8ull) 
        {
            __m128 xmm0;
            __m128 xmm1;
                                               
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               xmm0  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
               sum0  =  _mm_add_ps(xmm0,sum0);
               xmm1  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
               sum1  =  _mm_add_ps(xmm1,sum1);
                                                     
            }
             _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);
        }

        for(; (i+3ull) < this->m_nsamples; i += 4ull) 
        {
             __m128 xmm0;
                                                        
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               xmm0  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
               sum0  =  _mm_add_ps(xmm0,sum0);
                                                                   
            }
             _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
        }

        for(; (i+0ull) < this->m_nsamples; i += 1ull) 
        {
            float x0;
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               x0       =  this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i];
               sum_rem  =  x0+sum_rem;
                                                              
            }
            this->m_vsummed_sequence[i] = sum_rem;
        }
}


void 
gms::radiolocation
::am_bb_cos_signal_sequence_t
::signal_sequence_sse42_u10x_u(const float * __restrict__ sym_in ,
                               std::int32_t * __restrict__ ret_code) 
{
        if(__builtin_expect(nullptr==sym_in,0)) {return;}
        
        std::size_t i;
        __m128 sum0{_mm_setzero_ps()};
        __m128 sum1{_mm_setzero_ps()};
        __m128 sum2{_mm_setzero_ps()};
        __m128 sum3{_mm_setzero_ps()};
        __m128 sum4{_mm_setzero_ps()};
        __m128 sum5{_mm_setzero_ps()};
        __m128 sum6{_mm_setzero_ps()};
        __m128 sum7{_mm_setzero_ps()};
        __m128 sum8{_mm_setzero_ps()};
        __m128 sum9{_mm_setzero_ps()};
        float sum_rem{0.0f};

        for(std::size_t i{0ull}; i != this->m_nsignals; ++i) 
        {
            
            ret_code[i] = this->m_bb_cos_signals.operator[](i).create_signal_user_data(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](i));
           
        }

        for(i = 0ull; (i+39ull) < this->m_nsamples; i += 40ull) 
        {
            __m128 xmm0;
            __m128 xmm1;
            __m128 xmm2;
            __m128 xmm3;
            __m128 xmm4;
            __m128 xmm5;
            __m128 xmm6;
            __m128 xmm7;
            __m128 xmm8;
            __m128 xmm9;
           
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {
#if (AM_BB_COSINE_SEQUENCE_SOFT_PREFETCH) == 1
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+4ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+12ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+20ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+28ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+32ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+36ull],_MM_HINT_T1);
               
#endif 
               xmm0  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
               sum0  =  _mm_add_ps(xmm0,sum0);
               xmm1  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
               sum1  =  _mm_add_ps(xmm1,sum1);
               xmm2  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
               sum2  =  _mm_add_ps(xmm2,sum2);
               xmm3  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+12ull]);
               sum3  =  _mm_add_ps(xmm3,sum3);
               xmm4  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
               sum4  =  _mm_add_ps(xmm4,sum4);
               xmm5  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+20ull]);
               sum5  =  _mm_add_ps(xmm5,sum5);
               xmm6  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
               sum6  =  _mm_add_ps(xmm6,sum6);
               xmm7  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+28ull]);
               sum7  =  _mm_add_ps(xmm7,sum7);
               xmm8  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
               sum8  =  _mm_add_ps(xmm8,sum8);
               xmm9  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+36ull]);
               sum9  =  _mm_add_ps(xmm9,sum9);
              
            }
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

    
        for(; (i+31ull) < this->m_nsamples; i += 32ull) 
        {
            __m128 xmm0;
            __m128 xmm1;
            __m128 xmm2;
            __m128 xmm3;
            __m128 xmm4;
            __m128 xmm5;
            __m128 xmm6;
            __m128 xmm7;
                       
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               xmm0  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
               sum0  =  _mm_add_ps(xmm0,sum0);
               xmm1  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
               sum1  =  _mm_add_ps(xmm1,sum1);
               xmm2  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
               sum2  =  _mm_add_ps(xmm2,sum2);
               xmm3  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+12ull]);
               sum3  =  _mm_add_ps(xmm3,sum3);
               xmm4  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
               sum4  =  _mm_add_ps(xmm4,sum4);
               xmm5  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+20ull]);
               sum5  =  _mm_add_ps(xmm5,sum5);
               xmm6  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
               sum6  =  _mm_add_ps(xmm6,sum6);
               xmm7  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+28ull]);
               sum7  =  _mm_add_ps(xmm7,sum7);
                                        
            }
             _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum2);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+12ull],sum3);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum4);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+20ull],sum5);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum6);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+28ull],sum7);
        }

        for(; (i+15ull) < this->m_nsamples; i += 16ull) 
        {
            __m128 xmm0;
            __m128 xmm1;
            __m128 xmm2;
            __m128 xmm3;
                                   
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               xmm0  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
               sum0  =  _mm_add_ps(xmm0,sum0);
               xmm1  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
               sum1  =  _mm_add_ps(xmm1,sum1);
               xmm2  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
               sum2  =  _mm_add_ps(xmm2,sum2);
               xmm3  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+12ull]);
               sum3  =  _mm_add_ps(xmm3,sum3);
               
                                        
            }
             _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum2);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+12ull],sum3);
        }

        for(; (i+7ull) < this->m_nsamples; i += 8ull) 
        {
            __m128 xmm0;
            __m128 xmm1;
                                               
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               xmm0  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
               sum0  =  _mm_add_ps(xmm0,sum0);
               xmm1  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
               sum1  =  _mm_add_ps(xmm1,sum1);
                                                     
            }
             _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);
        }

        for(; (i+3ull) < this->m_nsamples; i += 4ull) 
        {
             __m128 xmm0;
                                                        
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               xmm0  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
               sum0  =  _mm_add_ps(xmm0,sum0);
                                                                   
            }
             _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
        }

        for(; (i+0ull) < this->m_nsamples; i += 1ull) 
        {
            float x0;
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               x0       =  this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i];
               sum_rem  =  x0+sum_rem;
                                                              
            }
            this->m_vsummed_sequence[i] = sum_rem;
        }
}


void 
gms::radiolocation
::am_bb_cos_signal_sequence_t
::signal_sequence_sse42_u6x_u(const float * __restrict__ sym_in ,
                               std::int32_t * __restrict__ ret_code) 
{
        if(__builtin_expect(nullptr==sym_in,0)) {return;}
        
        std::size_t i;
        __m128 sum0{_mm_setzero_ps()};
        __m128 sum1{_mm_setzero_ps()};
        __m128 sum2{_mm_setzero_ps()};
        __m128 sum3{_mm_setzero_ps()};
        __m128 sum4{_mm_setzero_ps()};
        __m128 sum5{_mm_setzero_ps()};
        float sum_rem{0.0f};

        for(std::size_t i{0ull}; i != this->m_nsignals; ++i) 
        {
            
            ret_code[i] = this->m_bb_cos_signals.operator[](i).create_signal_user_data(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](i));
           
        }

        for(i = 0ull; (i+23ull) < this->m_nsamples; i += 24ull) 
        {
            __m128 xmm0;
            __m128 xmm1;
            __m128 xmm2;
            __m128 xmm3;
            __m128 xmm4;
            __m128 xmm5;
                       
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {
#if (AM_BB_COSINE_SEQUENCE_SOFT_PREFETCH) == 1
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+4ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+12ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+20ull],_MM_HINT_T1);
                             
#endif 
               xmm0  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
               sum0  =  _mm_add_ps(xmm0,sum0);
               xmm1  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
               sum1  =  _mm_add_ps(xmm1,sum1);
               xmm2  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
               sum2  =  _mm_add_ps(xmm2,sum2);
               xmm3  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+12ull]);
               sum3  =  _mm_add_ps(xmm3,sum3);
               xmm4  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
               sum4  =  _mm_add_ps(xmm4,sum4);
               xmm5  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+20ull]);
               sum5  =  _mm_add_ps(xmm5,sum5);
                            
            }
             _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum2);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+12ull],sum3);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum4);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+20ull],sum5);
                     
        }
   
       
        for(; (i+15ull) < this->m_nsamples; i += 16ull) 
        {
            __m128 xmm0;
            __m128 xmm1;
            __m128 xmm2;
            __m128 xmm3;
                                   
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               xmm0  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
               sum0  =  _mm_add_ps(xmm0,sum0);
               xmm1  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
               sum1  =  _mm_add_ps(xmm1,sum1);
               xmm2  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
               sum2  =  _mm_add_ps(xmm2,sum2);
               xmm3  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+12ull]);
               sum3  =  _mm_add_ps(xmm3,sum3);
               
                                        
            }
             _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum2);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+12ull],sum3);
        }

        for(; (i+7ull) < this->m_nsamples; i += 8ull) 
        {
            __m128 xmm0;
            __m128 xmm1;
                                               
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               xmm0  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
               sum0  =  _mm_add_ps(xmm0,sum0);
               xmm1  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+4ull]);
               sum1  =  _mm_add_ps(xmm1,sum1);
                                                     
            }
             _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
             _mm_storeu_ps(&this->m_vsummed_sequence[i+4ull], sum1);
        }

        for(; (i+3ull) < this->m_nsamples; i += 4ull) 
        {
             __m128 xmm0;
                                                        
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               xmm0  =  _mm_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
               sum0  =  _mm_add_ps(xmm0,sum0);
                                                                   
            }
             _mm_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
        }

        for(; (i+0ull) < this->m_nsamples; i += 1ull) 
        {
            float x0;
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               x0       =  this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i];
               sum_rem  =  x0+sum_rem;
                                                              
            }
            this->m_vsummed_sequence[i] = sum_rem;
        }
}

void 
gms::radiolocation 
::am_bb_cos_signal_sequence_t
::signal_sequence_avx_u16x_u(const float * __restrict__ sym_in ,
                             std::int32_t * __restrict__ ret_codes) 
{
        if(__builtin_expect(nullptr==sym_in,0)) {return;}

        __m256 sum0{_mm256_setzero_ps()};
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
        for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
        {
            //std::printf("%d\n",this->m_nK_values.operator[](i));
            ret_codes[k] = this->m_bb_cos_signals.operator[](k).create_signal_user_data(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
        }

        for(i = 0ull; (i+127ull) < this->m_nsamples; i += 128ull) 
        {
            __m256 ymm0;
            __m256 ymm1;
            __m256 ymm2;
            __m256 ymm3;
            __m256 ymm4;
            __m256 ymm5;
            __m256 ymm6;
            __m256 ymm7;
            __m256 ymm8;
            __m256 ymm9;
            __m256 ymm10;
            __m256 ymm11;
            __m256 ymm12;
            __m256 ymm13;
            __m256 ymm14;
            __m256 ymm15;

              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
#if (AM_BB_COSINE_SEQUENCE_SOFT_PREFETCH) == 1
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+32ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+40ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+48ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+56ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+64ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+72ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+80ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+88ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+96ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+104ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+112ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+120ull],_MM_HINT_T1);
                                          
#endif 
                   
                    ymm0  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0,sum0);
                    ymm1  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1,sum1);
                    ymm2  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2,sum2);
                    ymm3  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3,sum3);
                    ymm4  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum4  =  _mm256_add_ps(ymm4,sum4);
                    ymm5  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+40ull]);
                    sum5  =  _mm256_add_ps(ymm5,sum5);
                    ymm6  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum6  =  _mm256_add_ps(ymm6,sum6); 
                    ymm7  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+56ull]);
                    sum7  =  _mm256_add_ps(ymm7,sum7);
                    ymm8  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+64ull]);
                    sum8  =  _mm256_add_ps(ymm8,sum8);
                    ymm9  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+72ull]);
                    sum9  =  _mm256_add_ps(ymm9,sum9);
                    ymm10 =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+80ull]);
                    sum10 =  _mm256_add_ps(ymm10,sum10);
                    ymm11 =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+88ull]);
                    sum11 =  _mm256_add_ps(ymm11,sum11);
                    ymm12 =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+96ull]);
                    sum12 =  _mm256_add_ps(ymm12,sum12);
                    ymm13 =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+104ull]);
                    sum13 =  _mm256_add_ps(ymm13,sum13);
                    ymm14 =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+112ull]);
                    sum14 =  _mm256_add_ps(ymm14,sum14);
                    ymm15  = _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+120ull]);
                    sum15  = _mm256_add_ps(ymm15,sum15);
              }

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
       //////////////////////////////////////////////////////
        for(; (i+95ull) < this->m_nsamples; i += 96ull) 
        {
            __m256 ymm0;
            __m256 ymm1;
            __m256 ymm2;
            __m256 ymm3;
            __m256 ymm4;
            __m256 ymm5;
            __m256 ymm6;
            __m256 ymm7;
            __m256 ymm8;
            __m256 ymm9;
            __m256 ymm10;
            __m256 ymm11;
           
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                   
                    ymm0  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0,sum0);
                    ymm1  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1,sum1);
                    ymm2  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2,sum2);
                    ymm3  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3,sum3);
                    ymm4  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum4  =  _mm256_add_ps(ymm4,sum4);
                    ymm5  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+40ull]);
                    sum5  =  _mm256_add_ps(ymm5,sum5);
                    ymm6  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum6  =  _mm256_add_ps(ymm6,sum6); 
                    ymm7  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+56ull]);
                    sum7  =  _mm256_add_ps(ymm7,sum7);
                    ymm8  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+64ull]);
                    sum8  =  _mm256_add_ps(ymm8,sum8);
                    ymm9  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+72ull]);
                    sum9  =  _mm256_add_ps(ymm9,sum9);
                    ymm10 =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+80ull]);
                    sum10 =  _mm256_add_ps(ymm10,sum10);
                    ymm11 =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+88ull]);
                    sum11 =  _mm256_add_ps(ymm11,sum11);
                  
              }

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
        //////////////////////////////////////////////////////////////////
        for(; (i+79ull) < this->m_nsamples; i += 80ull) 
        {
            __m256 ymm0;
            __m256 ymm1;
            __m256 ymm2;
            __m256 ymm3;
            __m256 ymm4;
            __m256 ymm5;
            __m256 ymm6;
            __m256 ymm7;
            __m256 ymm8;
            __m256 ymm9;
                      
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                   
                    ymm0  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0,sum0);
                    ymm1  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1,sum1);
                    ymm2  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2,sum2);
                    ymm3  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3,sum3);
                    ymm4  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum4  =  _mm256_add_ps(ymm4,sum4);
                    ymm5  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+40ull]);
                    sum5  =  _mm256_add_ps(ymm5,sum5);
                    ymm6  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum6  =  _mm256_add_ps(ymm6,sum6); 
                    ymm7  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+56ull]);
                    sum7  =  _mm256_add_ps(ymm7,sum7);
                    ymm8  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+64ull]);
                    sum8  =  _mm256_add_ps(ymm8,sum8);
                    ymm9  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+72ull]);
                    sum9  =  _mm256_add_ps(ymm9,sum9);
                                    
              }

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
        ///////////////////////////////////////////////////////////////
        for(; (i+63ull) < this->m_nsamples; i += 64ull) 
        {
            __m256 ymm0;
            __m256 ymm1;
            __m256 ymm2;
            __m256 ymm3;
            __m256 ymm4;
            __m256 ymm5;
            __m256 ymm6;
            __m256 ymm7;
                                  
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                   
                    ymm0  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0,sum0);
                    ymm1  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1,sum1);
                    ymm2  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2,sum2);
                    ymm3  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3,sum3);
                    ymm4  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum4  =  _mm256_add_ps(ymm4,sum4);
                    ymm5  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+40ull]);
                    sum5  =  _mm256_add_ps(ymm5,sum5);
                    ymm6  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum6  =  _mm256_add_ps(ymm6,sum6); 
                    ymm7  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+56ull]);
                    sum7  =  _mm256_add_ps(ymm7,sum7);
                                                       
              }

               _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum2);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum3);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum4);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+40ull],sum5);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum6);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+56ull],sum7);
                                            
        }
        /////////////////////////////////////////////////////////////////
        for(; (i+31ull) < this->m_nsamples; i += 32ull) 
        {
            __m256 ymm0;
            __m256 ymm1;
            __m256 ymm2;
            __m256 ymm3;
            __m256 ymm4;
            __m256 ymm5;
            __m256 ymm6;
            __m256 ymm7;
            __m256 ymm8;
            __m256 ymm9;
                      
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                   
                    ymm0  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0,sum0);
                    ymm1  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1,sum1);
                    ymm2  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2,sum2);
                    ymm3  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3,sum3);
                                                       
              }

               _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum2);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum3);
                                             
        }
        //////////////////////////////////////////////////////////////////
        for(; (i+15ull) < this->m_nsamples; i += 16ull) 
        {
            __m256 ymm0;
            __m256 ymm1;
                                  
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                   
                    ymm0  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0,sum0);
                    ymm1  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1,sum1);
                                                                           
              }

               _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);
                                                            
        }
        /////////////////////////////////////////////////////////////////
        for(; (i+7ull) < this->m_nsamples; i += 8ull) 
        {
            __m256 ymm0;
                                  
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                   
                    ymm0  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0,sum0);
                                                                         
              }

               _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                                                           
        }
        ///////////////////////////////////////////////////////////////
        for(; (i+0ull) < this->m_nsamples; i += 1ull) 
        {
            float x0;
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               x0       =  this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i];
               sum_rem  =  x0+sum_rem;
                                                              
            }
            this->m_vsummed_sequence[i] = sum_rem;
        }
}

void 
gms::radiolocation 
::am_bb_cos_signal_sequence_t
::signal_sequence_avx_u10x_u(const float * __restrict__ sym_in ,
                             std::int32_t * __restrict__ ret_codes) 
{
        if(__builtin_expect(nullptr==sym_in,0)) {return;}

        __m256 sum0{_mm256_setzero_ps()};
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
        for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
        {
            
            ret_codes[k] = this->m_bb_cos_signals.operator[](k).create_signal_user_data(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
        }

        for(i = 0ull; (i+79ull) < this->m_nsamples; i += 80ull) 
        {
            __m256 ymm0;
            __m256 ymm1;
            __m256 ymm2;
            __m256 ymm3;
            __m256 ymm4;
            __m256 ymm5;
            __m256 ymm6;
            __m256 ymm7;
            __m256 ymm8;
            __m256 ymm9;
            
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
#if (AM_BB_COSINE_SEQUENCE_SOFT_PREFETCH) == 1
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+32ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+40ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+48ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+56ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+64ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+72ull],_MM_HINT_T1);
                                                         
#endif 
                   
                    ymm0  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0,sum0);
                    ymm1  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1,sum1);
                    ymm2  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2,sum2);
                    ymm3  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3,sum3);
                    ymm4  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum4  =  _mm256_add_ps(ymm4,sum4);
                    ymm5  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+40ull]);
                    sum5  =  _mm256_add_ps(ymm5,sum5);
                    ymm6  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum6  =  _mm256_add_ps(ymm6,sum6); 
                    ymm7  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+56ull]);
                    sum7  =  _mm256_add_ps(ymm7,sum7);
                    ymm8  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+64ull]);
                    sum8  =  _mm256_add_ps(ymm8,sum8);
                    ymm9  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+72ull]);
                    sum9  =  _mm256_add_ps(ymm9,sum9);
                   
              }

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
           
        ///////////////////////////////////////////////////////////////
        for(; (i+63ull) < this->m_nsamples; i += 64ull) 
        {
            __m256 ymm0;
            __m256 ymm1;
            __m256 ymm2;
            __m256 ymm3;
            __m256 ymm4;
            __m256 ymm5;
            __m256 ymm6;
            __m256 ymm7;
                                  
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                   
                    ymm0  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0,sum0);
                    ymm1  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1,sum1);
                    ymm2  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2,sum2);
                    ymm3  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3,sum3);
                    ymm4  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum4  =  _mm256_add_ps(ymm4,sum4);
                    ymm5  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+40ull]);
                    sum5  =  _mm256_add_ps(ymm5,sum5);
                    ymm6  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+48ull]);
                    sum6  =  _mm256_add_ps(ymm6,sum6); 
                    ymm7  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+56ull]);
                    sum7  =  _mm256_add_ps(ymm7,sum7);
                                                       
              }

               _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum2);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum3);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum4);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+40ull],sum5);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+48ull],sum6);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+56ull],sum7);
                                            
        }
        /////////////////////////////////////////////////////////////////
        for(; (i+31ull) < this->m_nsamples; i += 32ull) 
        {
            __m256 ymm0;
            __m256 ymm1;
            __m256 ymm2;
            __m256 ymm3;
            __m256 ymm4;
            __m256 ymm5;
            __m256 ymm6;
            __m256 ymm7;
            __m256 ymm8;
            __m256 ymm9;
                      
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                   
                    ymm0  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0,sum0);
                    ymm1  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1,sum1);
                    ymm2  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2,sum2);
                    ymm3  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3,sum3);
                                                       
              }

               _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum2);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum3);
                                             
        }
        //////////////////////////////////////////////////////////////////
        for(; (i+15ull) < this->m_nsamples; i += 16ull) 
        {
            __m256 ymm0;
            __m256 ymm1;
                                  
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                   
                    ymm0  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0,sum0);
                    ymm1  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1,sum1);
                                                                           
              }

               _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);
                                                            
        }
        /////////////////////////////////////////////////////////////////
        for(; (i+7ull) < this->m_nsamples; i += 8ull) 
        {
            __m256 ymm0;
                                  
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                   
                    ymm0  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0,sum0);
                                                                         
              }

               _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                                                           
        }
        ///////////////////////////////////////////////////////////////
        for(; (i+0ull) < this->m_nsamples; i += 1ull) 
        {
            float x0;
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               x0       =  this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i];
               sum_rem  =  x0+sum_rem;
                                                              
            }
            this->m_vsummed_sequence[i] = sum_rem;
        }
}

void 
gms::radiolocation 
::am_bb_cos_signal_sequence_t
::signal_sequence_avx_u6x_u(const float * __restrict__ sym_in ,
                             std::int32_t * __restrict__ ret_codes) 
{
        if(__builtin_expect(nullptr==sym_in,0)) {return;}

        __m256 sum0{_mm256_setzero_ps()};
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
        for(std::size_t k{0ull}; k != this->m_nsignals; ++k) 
        {
            ret_codes[k] = this->m_bb_cos_signals.operator[](k).create_signal_user_data(sym_in,
                                                 static_cast<std::uint32_t>(this->m_nsamples),
                                                 this->m_nK_values.operator[](k));
        }

        for(i = 0ull; (i+47ull) < this->m_nsamples; i += 48ull) 
        {
            __m256 ymm0;
            __m256 ymm1;
            __m256 ymm2;
            __m256 ymm3;
            __m256 ymm4;
            __m256 ymm5;
            __m256 ymm6;
            __m256 ymm7;
            __m256 ymm8;
            __m256 ymm9;
            
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {
#if (AM_BB_COSINE_SEQUENCE_SOFT_PREFETCH) == 1
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+32ull],_MM_HINT_T1);
               _mm_prefetch((const char*)&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+40ull],_MM_HINT_T1);
                                                                        
#endif 
                   
                    ymm0  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0,sum0);
                    ymm1  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1,sum1);
                    ymm2  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2,sum2);
                    ymm3  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3,sum3);
                    ymm4  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+32ull]);
                    sum4  =  _mm256_add_ps(ymm4,sum4);
                    ymm5  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+40ull]);
                    sum5  =  _mm256_add_ps(ymm5,sum5);
                                     
              }

               _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum2);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum3);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+32ull],sum4);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+40ull],sum5);
             
        }
           
        /////////////////////////////////////////////////////////////////
        for(; (i+31ull) < this->m_nsamples; i += 32ull) 
        {
            __m256 ymm0;
            __m256 ymm1;
            __m256 ymm2;
            __m256 ymm3;
            __m256 ymm4;
            __m256 ymm5;
            __m256 ymm6;
            __m256 ymm7;
            __m256 ymm8;
            __m256 ymm9;
                      
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                   
                    ymm0  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0,sum0);
                    ymm1  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1,sum1);
                    ymm2  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+16ull]);
                    sum2  =  _mm256_add_ps(ymm2,sum2);
                    ymm3  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+24ull]);
                    sum3  =  _mm256_add_ps(ymm3,sum3);
                                                       
              }

               _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+16ull],sum2);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+24ull],sum3);
                                             
        }
        //////////////////////////////////////////////////////////////////
        for(; (i+15ull) < this->m_nsamples; i += 16ull) 
        {
            __m256 ymm0;
            __m256 ymm1;
                                  
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                   
                    ymm0  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0,sum0);
                    ymm1  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+8ull]);
                    sum1  =  _mm256_add_ps(ymm1,sum1);
                                                                           
              }

               _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
               _mm256_storeu_ps(&this->m_vsummed_sequence[i+8ull], sum1);
                                                            
        }
        /////////////////////////////////////////////////////////////////
        for(; (i+7ull) < this->m_nsamples; i += 8ull) 
        {
            __m256 ymm0;
                                  
              for(std::size_t j{0ull}; j != this->m_nsignals; ++j)  
              {

                   
                    ymm0  =  _mm256_loadu_ps(&this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i+0ull]);
                    sum0  =  _mm256_add_ps(ymm0,sum0);
                                                                         
              }

               _mm256_storeu_ps(&this->m_vsummed_sequence[i+0ull], sum0);
                                                           
        }
        ///////////////////////////////////////////////////////////////
        for(; (i+0ull) < this->m_nsamples; i += 1ull) 
        {
            float x0;
            for(std::size_t j{0ull}; j != this->m_nsignals; ++j) 
            {

               x0       =  this->m_bb_cos_signals.operator[](j).m_sig_samples.m_data[i];
               sum_rem  =  x0+sum_rem;
                                                              
            }
            this->m_vsummed_sequence[i] = sum_rem;
        }
}



auto
gms::radiolocation
::operator<<(std::ostream &os,
             gms::radiolocation::am_bb_cos_signal_sequence_t &rhs)->std::ostream &
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
