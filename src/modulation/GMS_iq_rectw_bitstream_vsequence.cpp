
#include <immintrin.h>
#include <fstream>
#include <iomanip>
#include "GMS_iq_rectw_bitstream_vsequence.h"

gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::iq_rectw_bitstream_vsequence_t() noexcept(true)
:
m_I_nsamples{0ull},
m_IQ_nsignals{0ull},
m_I_k_values{darray_r4_t()},
m_I_T_values{darray_r4_t()},
m_I_Ac_values{darray_r4_t()},
m_I_fc_values{darray_r4_t()},
m_I_ph0_values{darray_r4_t()},
m_Q_nsamples{0ull},
m_Q_k_values{darray_r4_t()},
m_Q_T_values{darray_r4_t()},
m_Q_Ac_values{darray_r4_t()},
m_Q_fc_values{darray_r4_t()},
m_Q_ph0_values{darray_r4_t()},
m_iq_rectw_bitstreams{std::vector<iq_rectwave_bitstream_t>()},
m_I_rw_bitstream_pdf_params{std::vector<I_rectw_bitstream_pdf_params_t>()},
m_Q_rw_bitstream_pdf_params{std::vector<Q_rectw_bitstream_pdf_params_t>()},
m_I_vsequence{darray_r4_t()},
m_Q_vsequence{darray_r4_t()}
{

}

gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::iq_rectw_bitstream_vsequence_t(const std::size_t   I_nsamples,
                                 const std::size_t   IQ_nsignals,
                                 const darray_r4_t & I_k_values,
                                 const darray_r4_t & I_T_values,
                                 const darray_r4_t & I_Ac_values,
                                 const darray_r4_t & I_fc_values,
                                 const darray_r4_t & I_ph0_values,
                                 const std::size_t   Q_nsamples,
                                 const darray_r4_t & Q_k_values,
                                 const darray_r4_t & Q_T_values,
                                 const darray_r4_t & Q_Ac_values,
                                 const darray_r4_t & Q_fc_values,
                                 const darray_r4_t & Q_ph0_values,
                                 const std::vector<I_rectw_bitstream_pdf_params_t> & I_rw_bitstream_pdf_params,
                                 const std::vector<Q_rectw_bitstream_pdf_params_t> & Q_rw_bitstream_pdf_params) noexcept(false)
:
m_I_nsamples{I_nsamples},
m_IQ_nsignals{IQ_nsignals},
m_I_k_values(I_k_values),
m_I_T_values(I_T_values),
m_I_Ac_values(I_Ac_values),
m_I_fc_values(I_fc_values),
m_I_ph0_values(I_ph0_values),
m_Q_nsamples{Q_nsamples},
m_Q_k_values(Q_k_values),
m_Q_T_values(Q_T_values),
m_Q_Ac_values(Q_Ac_values),
m_Q_fc_values(Q_fc_values),
m_Q_ph0_values(Q_ph0_values),
m_iq_rectw_bitstreams(std::vector<iq_rectwave_bitstream_t>(this->m_IQ_nsignals)),
m_I_rw_bitstream_pdf_params(I_rw_bitstream_pdf_params),
m_Q_rw_bitstream_pdf_params(Q_rw_bitstream_pdf_params),
m_I_vsequence{darray_r4_t(this->m_I_nsamples)},
m_Q_vsequence{darray_r4_t(this->m_Q_nsamples)}
{
    for(std::size_t i{0ull}; i != this->m_IQ_nsignals; ++i) 
    {
        const float s_I_k_value{this->m_I_k_values.m_data[i]};
        const float s_I_T_value{this->m_I_T_values.m_data[i]};
        const float s_I_Ac_value{this->m_I_Ac_values.m_data[i]};
        const float s_I_fc_value{this->m_I_fc_values.m_data[i]};
        const float s_I_ph0_value{this->m_I_ph0_values.m_data[i]};
        const float s_Q_k_value{this->m_Q_k_values.m_data[i]};
        const float s_Q_T_value{this->m_Q_T_values.m_data[i]};
        const float s_Q_Ac_value{this->m_Q_Ac_values.m_data[i]};
        const float s_Q_fc_value{this->m_Q_fc_values.m_data[i]};
        const float s_Q_ph0_value{this->m_Q_ph0_values.m_data[i]};

        this->m_iq_rectw_bitstreams.operator[](i) = iq_rectwave_bitstream_t(s_I_k_value,s_I_T_value,
                                                                            s_I_Ac_value,s_I_fc_value,
                                                                            s_I_ph0_value,s_Q_k_value,
                                                                            s_Q_T_value,s_Q_Ac_value,
                                                                            s_Q_fc_value,s_Q_ph0_value);
    }
}

gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::iq_rectw_bitstream_vsequence_t(const iq_rectw_bitstream_vsequence_t &other) noexcept(false)
:
m_I_nsamples{other.m_I_nsamples},
m_IQ_nsignals{other.m_IQ_nsignals},
m_I_k_values{other.m_I_k_values},
m_I_T_values{other.m_I_T_values},
m_I_Ac_values{other.m_I_Ac_values},
m_I_fc_values{other.m_I_fc_values},
m_I_ph0_values{other.m_I_ph0_values},
m_Q_nsamples{other.m_Q_nsamples},
m_Q_k_values{other.m_Q_k_values},
m_Q_T_values{other.m_Q_T_values},
m_Q_Ac_values{other.m_Q_Ac_values},
m_Q_fc_values{other.m_Q_fc_values},
m_Q_ph0_values{other.m_Q_ph0_values},
m_iq_rectw_bitstreams{other.m_iq_rectw_bitstreams},
m_I_rw_bitstream_pdf_params{other.m_I_rw_bitstream_pdf_params},
m_Q_rw_bitstream_pdf_params{other.m_Q_rw_bitstream_pdf_params},
m_I_vsequence{other.m_I_vsequence},
m_Q_vsequence{other.m_Q_vsequence}
{

}

gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::iq_rectw_bitstream_vsequence_t(iq_rectw_bitstream_vsequence_t &&other) noexcept(true)
:
m_I_nsamples{other.m_I_nsamples},
m_IQ_nsignals{other.m_IQ_nsignals},
m_I_k_values{std::move(other.m_I_k_values)},
m_I_T_values{std::move(other.m_I_T_values)},
m_I_Ac_values{std::move(other.m_I_Ac_values)},
m_I_fc_values{std::move(other.m_I_fc_values)},
m_I_ph0_values{std::move(other.m_I_ph0_values)},
m_Q_nsamples{other.m_Q_nsamples},
m_Q_k_values{std::move(other.m_Q_k_values)},
m_Q_T_values{std::move(other.m_Q_T_values)},
m_Q_Ac_values{std::move(other.m_Q_Ac_values)},
m_Q_fc_values{std::move(other.m_Q_fc_values)},
m_Q_ph0_values{std::move(other.m_Q_ph0_values)},
m_iq_rectw_bitstreams{std::move(other.m_iq_rectw_bitstreams)},
m_I_rw_bitstream_pdf_params{std::move(other.m_I_rw_bitstream_pdf_params)},
m_Q_rw_bitstream_pdf_params{std::move(other.m_Q_rw_bitstream_pdf_params)},
m_I_vsequence{std::move(other.m_I_vsequence)},
m_Q_vsequence{std::move(other.m_Q_vsequence)}
{

}

gms::radiolocation 
::iq_rectw_bitstream_vsequence_t
::~iq_rectw_bitstream_vsequence_t()
{

}

gms::radiolocation 
::iq_rectw_bitstream_vsequence_t &
gms::radiolocation 
::iq_rectw_bitstream_vsequence_t
::operator=(const iq_rectw_bitstream_vsequence_t &other) noexcept(false)
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_I_nsamples         = other.m_I_nsamples;
    this->m_IQ_nsignals        = other.m_IQ_nsignals;
    this->m_I_k_values.operator=(other.m_I_k_values);
    this->m_I_T_values.operator=(other.m_I_T_values);
    this->m_I_Ac_values.operator=(other.m_I_Ac_values);
    this->m_I_fc_values.operator=(other.m_I_fc_values);
    this->m_I_ph0_values.operator=(other.m_I_ph0_values);
    this->m_Q_nsamples           = other.m_Q_nsamples;
    this->m_Q_k_values.operator=(other.m_Q_k_values);
    this->m_Q_T_values.operator=(other.m_Q_T_values);
    this->m_Q_Ac_values.operator=(other.m_Q_Ac_values);
    this->m_Q_fc_values.operator=(other.m_Q_fc_values);
    this->m_Q_ph0_values.operator=(other.m_Q_ph0_values);
    this->m_iq_rectw_bitstreams.operator=(other.m_iq_rectw_bitstreams);
    this->m_I_rw_bitstream_pdf_params.operator=(other.m_I_rw_bitstream_pdf_params);
    this->m_Q_rw_bitstream_pdf_params.operator=(other.m_Q_rw_bitstream_pdf_params);
    this->m_I_vsequence.operator=(other.m_I_vsequence);
    this->m_Q_vsequence.operator=(other.m_Q_vsequence);
    return (*this);
}

gms::radiolocation 
::iq_rectw_bitstream_vsequence_t &
gms::radiolocation 
::iq_rectw_bitstream_vsequence_t
::operator=(iq_rectw_bitstream_vsequence_t &&other) noexcept(true)
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_I_nsamples         = other.m_I_nsamples;
    this->m_IQ_nsignals        = other.m_IQ_nsignals;
    this->m_I_k_values.operator=(std::move(other.m_I_k_values));
    this->m_I_T_values.operator=(std::move(other.m_I_T_values));
    this->m_I_Ac_values.operator=(std::move(other.m_I_Ac_values));
    this->m_I_fc_values.operator=(std::move(other.m_I_fc_values));
    this->m_I_ph0_values.operator=(std::move(other.m_I_ph0_values));
    this->m_Q_nsamples           = other.m_Q_nsamples;
    this->m_Q_k_values.operator=(std::move(other.m_Q_k_values));
    this->m_Q_T_values.operator=(std::move(other.m_Q_T_values));
    this->m_Q_Ac_values.operator=(std::move(other.m_Q_Ac_values));
    this->m_Q_fc_values.operator=(std::move(other.m_Q_fc_values));
    this->m_Q_ph0_values.operator=(std::move(other.m_Q_ph0_values));
    this->m_iq_rectw_bitstreams.operator=(std::move(other.m_iq_rectw_bitstreams));
    this->m_I_rw_bitstream_pdf_params.operator=(std::move(other.m_I_rw_bitstream_pdf_params));
    this->m_Q_rw_bitstream_pdf_params.operator=(std::move(other.m_Q_rw_bitstream_pdf_params));
    this->m_I_vsequence.operator=(std::move(other.m_I_vsequence));
    this->m_Q_vsequence.operator=(std::move(other.m_Q_vsequence));
    return (*this);
}

void 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
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

template<gms::radiolocation::iq_rectw_bitstream_vsequence_t
         ::IQ_rectw_bitstreams_optim_path  optim_path,
         gms::radiolocation::iq_rectw_bitstream_vsequence_t
         ::IQ_rectw_bitstream_max_signal_components max_signals>       
std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x(std::int32_t * I_ch_ret_codes,
                                                std::int32_t * Q_ch_ret_codes)
{
       if(__builtin_expect(nullptr==I_ch_ret_codes,0)) { return (-1);}
       if(__builtin_expect(nullptr==Q_ch_ret_codes,0)) { return (-2);}
       //if(__builtin_expect(this->m_IQ_nsignals>256ull,0)) { return (-2);}

#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
           PMC_VARS
           HW_PMC_COLLECTION_PROLOGE_BODY
#endif 
        constexpr std::size_t max_signals_buff_size{static_cast<std::size_t>(max_signals)};
        static_assert(max_signals_buff_size>0ull,"Invalid [IQ_rectw_bitstream_max_signal] Parameter!!");
        __ATTR_ALIGN__(16) 
        __m128 xmm0_i[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm1_i[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm2_i[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm3_i[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm4_i[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm5_i[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm6_i[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm7_i[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm8_i[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm9_i[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm10_i[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm11_i[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm12_i[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm13_i[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm14_i[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm15_i[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm0_q[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm1_q[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm2_q[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm3_q[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm4_q[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm5_q[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm6_q[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm7_q[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm8_q[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm9_q[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm10_q[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm11_q[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm12_q[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm13_q[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm14_q[max_signals_buff_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm15_q[max_signals_buff_size];
        float scal_rem_i[max_signals_buff_size];
        float scal_rem_q[max_signals_buff_size];
        __m128 vzero{_mm_setzero_ps()};
        __m128 sum0_i{vzero};
        __m128 sum1_i{sum0_i};
        __m128 sum2_i{sum0_i};
        __m128 sum3_i{sum0_i};
        __m128 sum4_i{sum0_i};
        __m128 sum5_i{sum0_i};
        __m128 sum6_i{sum0_i};
        __m128 sum7_i{sum0_i};
        __m128 sum8_i{sum0_i};
        __m128 sum9_i{sum0_i};
        __m128 sum10_i{sum0_i};
        __m128 sum11_i{sum0_i};
        __m128 sum12_i{sum0_i};
        __m128 sum13_i{sum0_i};
        __m128 sum14_i{sum0_i};
        __m128 sum15_i{sum0_i};
        __m128 sum0_q{sum0_i};
        __m128 sum1_q{sum0_i};
        __m128 sum2_q{sum0_i};
        __m128 sum3_q{sum0_i};
        __m128 sum4_q{sum0_i};
        __m128 sum5_q{sum0_i};
        __m128 sum6_q{sum0_i};
        __m128 sum7_q{sum0_i};
        __m128 sum8_q{sum0_i};
        __m128 sum9_q{sum0_i};
        __m128 sum10_q{sum0_i};
        __m128 sum11_q{sum0_i};
        __m128 sum12_q{sum0_i};
        __m128 sum13_q{sum0_i};
        __m128 sum14_q{sum0_i};
        __m128 sum15_q{sum0_i}; 
        std::size_t i;
        float sum_rem_i{0.0f};
        float sum_rem_q{0.0f};

        if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SCALAR_PATH)
        {          
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    I_ch_ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_scalar();
                    Q_ch_ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_Q_channel_bitstream_scalar();
                }            
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    I_ch_ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_sse();
                    Q_ch_ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_Q_channel_bitstream_sse();
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    I_ch_ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_sse_u4x();
                    Q_ch_ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_Q_channel_bitstream_sse_u4x();
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    I_ch_ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx(true);
                    Q_ch_ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_Q_channel_bitstream_avx(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    I_ch_ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx_u4x(true);
                    Q_ch_ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_Q_channel_bitstream_avx_u4x(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    I_ch_ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx_u4x(true);
                    Q_ch_ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_Q_channel_bitstream_avx_u4x(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    I_ch_ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx512_u4x(true);
                    Q_ch_ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_Q_channel_bitstream_avx512_u4x(true);
                } 
        }

        if(this->m_I_nsamples==this->m_Q_nsamples)
        {
              for(i = 0ull; (i+63ull) < this->m_I_nsamples; i += 64ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    sum6_i = vzero;
                    sum6_q = vzero;
                    sum7_i = vzero;
                    sum7_q = vzero;
                    sum8_i = vzero;
                    sum8_q = vzero;
                    sum9_i = vzero;
                    sum9_q = vzero;
                    sum10_i = vzero;
                    sum10_q = vzero;
                    sum11_i = vzero;
                    sum11_q = vzero;
                    sum12_i = vzero;
                    sum12_q = vzero;
                    sum13_i = vzero;
                    sum13_q = vzero;
                    sum14_i = vzero;
                    sum14_q = vzero;
                    sum15_i = vzero;
                    sum15_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);
#endif 
                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);
                        xmm2_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum2_q    = _mm_add_ps(xmm2_q[j],sum2_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum2_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum2_q);
                        xmm3_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+12ull]);
                        sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+12ull]);
                        sum3_q    = _mm_add_ps(xmm3_q[j],sum3_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+12ull],sum3_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+12ull],sum3_q);
                        xmm4_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum4_q    = _mm_add_ps(xmm4_q[j],sum4_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum4_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum4_q);
                        xmm5_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+20ull]);
                        sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+20ull]);
                        sum5_q    = _mm_add_ps(xmm5_q[j],sum5_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+20ull],sum5_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+20ull],sum5_q);
                        xmm6_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum6_i    = _mm_add_ps(xmm6_i[j],sum6_i);
                        xmm6_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum6_q    = _mm_add_ps(xmm6_q[j],sum6_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum6_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum6_q);
                        xmm7_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+28ull]);
                        sum7_i    = _mm_add_ps(xmm7_i[j],sum7_i);
                        xmm7_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+28ull]);
                        sum7_q    = _mm_add_ps(xmm7_q[j],sum7_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+28ull],sum7_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+28ull],sum7_q);
                        xmm8_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum8_i    = _mm_add_ps(xmm8_i[j],sum8_i);
                        xmm8_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum8_q    = _mm_add_ps(xmm8_q[j],sum8_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum8_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum8_q);
                        xmm9_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+36ull]);
                        sum9_i    = _mm_add_ps(xmm9_i[j],sum9_i);
                        xmm9_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+36ull]);
                        sum9_q    = _mm_add_ps(xmm9_q[j],sum9_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+36ull],sum9_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+36ull],sum9_q);
                        xmm10_i[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+40ull]);
                        sum10_i   = _mm_add_ps(xmm10_i[j],sum10_i);
                        xmm10_q[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+40ull]);
                        sum10_q   = _mm_add_ps(xmm10_q[j],sum10_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+40ull],sum10_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+40ull],sum10_q);
                        xmm11_i[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+44ull]);
                        sum11_i   = _mm_add_ps(xmm11_i[j],sum11_i);
                        xmm11_q[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+44ull]);
                        sum11_q   = _mm_add_ps(xmm11_q[j],sum11_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+44ull],sum11_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+44ull],sum11_q);
                        xmm12_i[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum12_i   = _mm_add_ps(xmm12_i[j],sum12_i);
                        xmm12_q[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum12_q   = _mm_add_ps(xmm12_q[j],sum12_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum12_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum12_q);
                        xmm13_i[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+52ull]);
                        sum13_i   = _mm_add_ps(xmm13_i[j],sum13_i);
                        xmm13_q[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+52ull]);
                        sum13_q   = _mm_add_ps(xmm13_q[j],sum13_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+52ull],sum13_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+52ull],sum13_q);
                        xmm14_i[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+56ull]);
                        sum14_i   = _mm_add_ps(xmm14_i[j],sum14_i);
                        xmm14_q[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+56ull]);
                        sum14_q   = _mm_add_ps(xmm14_q[j],sum14_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+56ull],sum14_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+56ull],sum14_q);
                        xmm15_i[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+60ull]);
                        sum15_i   = _mm_add_ps(xmm15_i[j],sum15_i);
                        xmm15_q[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+60ull]);
                        sum15_q   = _mm_add_ps(xmm15_q[j],sum15_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+60ull],sum15_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+60ull],sum15_q);
                       
                    }
              }

              for(; (i+39ull) < this->m_I_nsamples; i += 40ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    sum6_i = vzero;
                    sum6_q = vzero;
                    sum7_i = vzero;
                    sum7_q = vzero;
                    sum8_i = vzero;
                    sum8_q = vzero;
                    sum9_i = vzero;
                    sum9_q = vzero;
                   
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);
                        xmm2_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum2_q    = _mm_add_ps(xmm2_q[j],sum2_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum2_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum2_q);
                        xmm3_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+12ull]);
                        sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+12ull]);
                        sum3_q    = _mm_add_ps(xmm3_q[j],sum3_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+12ull],sum3_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+12ull],sum3_q);
                        xmm4_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum4_q    = _mm_add_ps(xmm4_q[j],sum4_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum4_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum4_q);
                        xmm5_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+20ull]);
                        sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+20ull]);
                        sum5_q    = _mm_add_ps(xmm5_q[j],sum5_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+20ull],sum5_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+20ull],sum5_q);
                        xmm6_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum6_i    = _mm_add_ps(xmm6_i[j],sum6_i);
                        xmm6_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum6_q    = _mm_add_ps(xmm6_q[j],sum6_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum6_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum6_q);
                        xmm7_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+28ull]);
                        sum7_i    = _mm_add_ps(xmm7_i[j],sum7_i);
                        xmm7_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+28ull]);
                        sum7_q    = _mm_add_ps(xmm7_q[j],sum7_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+28ull],sum7_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+28ull],sum7_q);
                        xmm8_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum8_i    = _mm_add_ps(xmm8_i[j],sum8_i);
                        xmm8_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum8_q    = _mm_add_ps(xmm8_q[j],sum8_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum8_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum8_q);
                        xmm9_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+36ull]);
                        sum9_i    = _mm_add_ps(xmm9_i[j],sum9_i);
                        xmm9_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+36ull]);
                        sum9_q    = _mm_add_ps(xmm9_q[j],sum9_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+36ull],sum9_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+36ull],sum9_q);
                                               
                    }
              }

               for(; (i+31ull) < this->m_I_nsamples; i += 32ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    sum6_i = vzero;
                    sum6_q = vzero;
                    sum7_i = vzero;
                    sum7_q = vzero;
                                       
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);
                        xmm2_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum2_q    = _mm_add_ps(xmm2_q[j],sum2_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum2_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum2_q);
                        xmm3_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+12ull]);
                        sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+12ull]);
                        sum3_q    = _mm_add_ps(xmm3_q[j],sum3_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+12ull],sum3_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+12ull],sum3_q);
                        xmm4_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum4_q    = _mm_add_ps(xmm4_q[j],sum4_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum4_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum4_q);
                        xmm5_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+20ull]);
                        sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+20ull]);
                        sum5_q    = _mm_add_ps(xmm5_q[j],sum5_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+20ull],sum5_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+20ull],sum5_q);
                        xmm6_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum6_i    = _mm_add_ps(xmm6_i[j],sum6_i);
                        xmm6_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum6_q    = _mm_add_ps(xmm6_q[j],sum6_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum6_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum6_q);
                        xmm7_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+28ull]);
                        sum7_i    = _mm_add_ps(xmm7_i[j],sum7_i);
                        xmm7_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+28ull]);
                        sum7_q    = _mm_add_ps(xmm7_q[j],sum7_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+28ull],sum7_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+28ull],sum7_q);
                                                                      
                    }
              }

               for(; (i+15ull) < this->m_I_nsamples; i += 16ull) 
               {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                                      
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);
                        xmm2_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum2_q    = _mm_add_ps(xmm2_q[j],sum2_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum2_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum2_q);
                        xmm3_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+12ull]);
                        sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+12ull]);
                        sum3_q    = _mm_add_ps(xmm3_q[j],sum3_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+12ull],sum3_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+12ull],sum3_q);
                       
                    }
              }

              for(; (i+7ull) < this->m_I_nsamples; i += 8ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                                      
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);
                                                                     
                    }
              }

              for(; (i+3ull) < this->m_I_nsamples; i += 4ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                                                          
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                                                                                             
                    }
              }

              for(; (i+0ull) < this->m_I_nsamples; i += 1ull) 
              {
                   sum_rem_i = 0.0f;
                   sum_rem_q = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_i[j] = this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i];
                       sum_rem_i     += scal_rem_i[j];
                       this->m_I_vsequence.m_data[i] = sum_rem_i;
                       scal_rem_q[j] = this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i];
                       sum_rem_q     += scal_rem_q[j];
                       this->m_Q_vsequence.m_data[i] = sum_rem_q;
                   }
              }
        }
        else 
        {
             for(i = 0ull; (i+63ull) < this->m_I_nsamples; i += 64ull) 
             {
                    sum0_i = vzero;
                    sum1_i = vzero;
                    sum2_i = vzero; 
                    sum3_i = vzero;                   
                    sum4_i = vzero;                    
                    sum5_i = vzero;                   
                    sum6_i = vzero;                    
                    sum7_i = vzero;                    
                    sum8_i = vzero;                    
                    sum9_i = vzero;                  
                    sum10_i = vzero;                    
                    sum11_i = vzero;                    
                    sum12_i = vzero;                   
                    sum13_i = vzero;                   
                    sum14_i = vzero;                   
                    sum15_i = vzero;
                    
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);                      
#endif 
                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);                       
                        xmm2_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i); 
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum2_i);                        
                        xmm3_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+12ull]);
                        sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+12ull],sum3_i);                        
                        xmm4_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum4_i);                        
                        xmm5_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+20ull]);
                        sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+20ull],sum5_i);                       
                        xmm6_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum6_i    = _mm_add_ps(xmm6_i[j],sum6_i);                      
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum6_i);                       
                        xmm7_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+28ull]);
                        sum7_i    = _mm_add_ps(xmm7_i[j],sum7_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+28ull],sum7_i);                       
                        xmm8_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum8_i    = _mm_add_ps(xmm8_i[j],sum8_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum8_i);                       
                        xmm9_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+36ull]);
                        sum9_i    = _mm_add_ps(xmm9_i[j],sum9_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+36ull],sum9_i);                        
                        xmm10_i[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+40ull]);
                        sum10_i   = _mm_add_ps(xmm10_i[j],sum10_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+40ull],sum10_i);                       
                        xmm11_i[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+44ull]);
                        sum11_i   = _mm_add_ps(xmm11_i[j],sum11_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+44ull],sum11_i);                       
                        xmm12_i[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum12_i   = _mm_add_ps(xmm12_i[j],sum12_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum12_i);                        
                        xmm13_i[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+52ull]);
                        sum13_i   = _mm_add_ps(xmm13_i[j],sum13_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+52ull],sum13_i);                       
                        xmm14_i[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+56ull]);
                        sum14_i   = _mm_add_ps(xmm14_i[j],sum14_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+56ull],sum14_i);                        
                        xmm15_i[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+60ull]);
                        sum15_i   = _mm_add_ps(xmm15_i[j],sum15_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+60ull],sum15_i);
                                               
                    }
              }

               for(; (i+39ull) < this->m_I_nsamples; i += 40ull) 
               {
                    sum0_i = vzero;
                    sum1_i = vzero;
                    sum2_i = vzero; 
                    sum3_i = vzero;                   
                    sum4_i = vzero;                    
                    sum5_i = vzero;                   
                    sum6_i = vzero;                    
                    sum7_i = vzero;                    
                    sum8_i = vzero;                    
                    sum9_i = vzero;                  
                                        
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);                       
                        xmm2_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i); 
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum2_i);                        
                        xmm3_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+12ull]);
                        sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+12ull],sum3_i);                        
                        xmm4_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum4_i);                        
                        xmm5_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+20ull]);
                        sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+20ull],sum5_i);                       
                        xmm6_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum6_i    = _mm_add_ps(xmm6_i[j],sum6_i);                      
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum6_i);                       
                        xmm7_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+28ull]);
                        sum7_i    = _mm_add_ps(xmm7_i[j],sum7_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+28ull],sum7_i);                       
                        xmm8_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum8_i    = _mm_add_ps(xmm8_i[j],sum8_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum8_i);                       
                        xmm9_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+36ull]);
                        sum9_i    = _mm_add_ps(xmm9_i[j],sum9_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+36ull],sum9_i);                        
                                                                       
                    }
              }

               for(; (i+31ull) < this->m_I_nsamples; i += 32ull) 
               {
                    sum0_i = vzero;
                    sum1_i = vzero;
                    sum2_i = vzero; 
                    sum3_i = vzero;                   
                    sum4_i = vzero;                    
                    sum5_i = vzero;                   
                    sum6_i = vzero;                    
                    sum7_i = vzero;                    
                                       
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);                       
                        xmm2_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i); 
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum2_i);                        
                        xmm3_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+12ull]);
                        sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+12ull],sum3_i);                        
                        xmm4_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum4_i);                        
                        xmm5_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+20ull]);
                        sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+20ull],sum5_i);                       
                        xmm6_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum6_i    = _mm_add_ps(xmm6_i[j],sum6_i);                      
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum6_i);                       
                        xmm7_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+28ull]);
                        sum7_i    = _mm_add_ps(xmm7_i[j],sum7_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+28ull],sum7_i);                       
                                                                      
                    }
              }

               for(; (i+15ull) < this->m_I_nsamples; i += 16ull) 
               {
                    sum0_i = vzero;
                    sum1_i = vzero;
                    sum2_i = vzero; 
                    sum3_i = vzero;                   
                                        
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);                       
                        xmm2_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i); 
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum2_i);                        
                        xmm3_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+12ull]);
                        sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+12ull],sum3_i);                        
                                                                      
                    }
              }

               for(; (i+7ull) < this->m_I_nsamples; i += 8ull) 
               {
                    sum0_i = vzero;
                    sum1_i = vzero;                 
                                        
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);                                        
                                                                      
                    }
              }

               for(; (i+3ull) < this->m_I_nsamples; i += 4ull) 
               {
                    sum0_i = vzero;    
                                        
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                                      
                                                                      
                    }
              }
              
              for(; (i+0ull) < this->m_I_nsamples; i += 1ull) 
              {
                   sum_rem_i = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_i[j] = this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i];
                       sum_rem_i     += scal_rem_i[j];
                       this->m_I_vsequence.m_data[i] = sum_rem_i;
                       
                   }
              }

              ///////////////////////////// Q-channel //////////////////////////////////////////////
              
               for(i = 0ull; (i+63ull) < this->m_Q_nsamples; i += 64ull) 
               {
                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                    
                    sum6_q = vzero;                    
                    sum7_q = vzero;                    
                    sum8_q = vzero;                    
                    sum9_q = vzero;                    
                    sum10_q = vzero;                    
                    sum11_q = vzero;                    
                    sum12_q = vzero;                   
                    sum13_q = vzero;                    
                    sum14_q = vzero;                   
                    sum15_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);
#endif 
                       
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);                      
                        xmm2_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum2_q    = _mm_add_ps(xmm2_q[j],sum2_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum2_q);                       
                        xmm3_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+12ull]);
                        sum3_q    = _mm_add_ps(xmm3_q[j],sum3_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+12ull],sum3_q);                      
                        xmm4_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum4_q    = _mm_add_ps(xmm4_q[j],sum4_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum4_q);                        
                        xmm5_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+20ull]);
                        sum5_q    = _mm_add_ps(xmm5_q[j],sum5_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+20ull],sum5_q);                        
                        xmm6_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum6_q    = _mm_add_ps(xmm6_q[j],sum6_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum6_q);                      
                        xmm7_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+28ull]);
                        sum7_q    = _mm_add_ps(xmm7_q[j],sum7_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+28ull],sum7_q);                       
                        xmm8_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum8_q    = _mm_add_ps(xmm8_q[j],sum8_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum8_q);                       
                        xmm9_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+36ull]);
                        sum9_q    = _mm_add_ps(xmm9_q[j],sum9_q);                      
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+36ull],sum9_q);                       
                        xmm10_q[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+40ull]);
                        sum10_q   = _mm_add_ps(xmm10_q[j],sum10_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+40ull],sum10_q);                        
                        xmm11_q[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+44ull]);
                        sum11_q   = _mm_add_ps(xmm11_q[j],sum11_q);                      
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+44ull],sum11_q);                       
                        xmm12_q[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum12_q   = _mm_add_ps(xmm12_q[j],sum12_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum12_q);                        
                        xmm13_q[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+52ull]);
                        sum13_q   = _mm_add_ps(xmm13_q[j],sum13_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+52ull],sum13_q);                       
                        xmm14_q[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+56ull]);
                        sum14_q   = _mm_add_ps(xmm14_q[j],sum14_q);                      
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+56ull],sum14_q);                       
                        xmm15_q[j]= _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+60ull]);
                        sum15_q   = _mm_add_ps(xmm15_q[j],sum15_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+60ull],sum15_q);
                       
                    }
              }
                
              for(; (i+39ull) < this->m_Q_nsamples; i += 40ull) 
              {
                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                    
                    sum6_q = vzero;                    
                    sum7_q = vzero;                    
                    sum8_q = vzero;                    
                    sum9_q = vzero;                                       
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {                       
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);                      
                        xmm2_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum2_q    = _mm_add_ps(xmm2_q[j],sum2_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum2_q);                       
                        xmm3_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+12ull]);
                        sum3_q    = _mm_add_ps(xmm3_q[j],sum3_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+12ull],sum3_q);                      
                        xmm4_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum4_q    = _mm_add_ps(xmm4_q[j],sum4_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum4_q);                        
                        xmm5_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+20ull]);
                        sum5_q    = _mm_add_ps(xmm5_q[j],sum5_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+20ull],sum5_q);                        
                        xmm6_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum6_q    = _mm_add_ps(xmm6_q[j],sum6_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum6_q);                      
                        xmm7_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+28ull]);
                        sum7_q    = _mm_add_ps(xmm7_q[j],sum7_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+28ull],sum7_q);                       
                        xmm8_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum8_q    = _mm_add_ps(xmm8_q[j],sum8_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum8_q);                       
                        xmm9_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+36ull]);
                        sum9_q    = _mm_add_ps(xmm9_q[j],sum9_q);                      
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+36ull],sum9_q);                       
                                              
                    }
              }

              for(; (i+31ull) < this->m_Q_nsamples; i += 32ull) 
              {
                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                    
                    sum6_q = vzero;                    
                    sum7_q = vzero;                                                   
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {                       
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);                      
                        xmm2_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum2_q    = _mm_add_ps(xmm2_q[j],sum2_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum2_q);                       
                        xmm3_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+12ull]);
                        sum3_q    = _mm_add_ps(xmm3_q[j],sum3_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+12ull],sum3_q);                      
                        xmm4_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum4_q    = _mm_add_ps(xmm4_q[j],sum4_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum4_q);                        
                        xmm5_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+20ull]);
                        sum5_q    = _mm_add_ps(xmm5_q[j],sum5_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+20ull],sum5_q);                        
                        xmm6_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum6_q    = _mm_add_ps(xmm6_q[j],sum6_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum6_q);                      
                        xmm7_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+28ull]);
                        sum7_q    = _mm_add_ps(xmm7_q[j],sum7_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+28ull],sum7_q);                       
                                                                                      
                    }
              }

              for(; (i+15ull) < this->m_Q_nsamples; i += 16ull) 
              {
                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                                                
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {                       
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);                      
                        xmm2_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum2_q    = _mm_add_ps(xmm2_q[j],sum2_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum2_q);                       
                        xmm3_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+12ull]);
                        sum3_q    = _mm_add_ps(xmm3_q[j],sum3_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+12ull],sum3_q);                                      
                                              
                    }
              }

              for(; (i+7ull) < this->m_Q_nsamples; i += 8ull) 
              {                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                                                      
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                       
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);                                            
                                              
                    }
              }

              for(; (i+3ull) < this->m_Q_nsamples; i += 4ull) 
              {                    
                    sum0_q = vzero;                                                       
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {                       
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                                             
                                              
                    }
              }

               for(; (i+0ull) < this->m_Q_nsamples; i += 1ull) 
               {                 
                   sum_rem_q = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_q[j] = this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i];
                       sum_rem_q     += scal_rem_q[j];
                       this->m_Q_vsequence.m_data[i] = sum_rem_q;
                   }
              }
        }
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
           HW_PMC_COLLECTION_EPILOGE_BODY
           HW_PMC_SHOW_INFO
#endif 
     return (0);
}

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SCALAR_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_25>(std::int32_t * __restrict__,
                                                                                                                  std::int32_t * __restrict__);
                                              

                                              
template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_25>(std::int32_t * __restrict__,
                                                                                                                  std::int32_t * __restrict__);
                                             


template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_25>(std::int32_t * __restrict__,
                                                                                                                  std::int32_t * __restrict__);
                                             

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_25>(std::int32_t * __restrict__,
                                                                                                                  std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_25>(std::int32_t * __restrict__,
                                                                                                                  std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_25>(std::int32_t * __restrict__,
                                                                                                                  std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_25>(std::int32_t * __restrict__,
                                                                                                                  std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SCALAR_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_50>(std::int32_t * __restrict__,
                                                                                                                  std::int32_t * __restrict__);
                                              

                                              
template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_50>(std::int32_t * __restrict__,
                                                                                                                  std::int32_t * __restrict__);
                                             


template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_50>(std::int32_t * __restrict__,
                                                                                                                  std::int32_t * __restrict__);
                                             

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_50>(std::int32_t * __restrict__,
                                                                                                                  std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_50>(std::int32_t * __restrict__,
                                                                                                                  std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_50>(std::int32_t * __restrict__,
                                                                                                                  std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_50>(std::int32_t * __restrict__,
                                                                                                                  std::int32_t * __restrict__);


template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SCALAR_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_100>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                              

                                              
template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_100>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                             


template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_100>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                             

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_100>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_100>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_100>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_100>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SCALAR_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_150>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                              

                                              
template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_150>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                             


template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_150>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                             

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_150>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_150>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_150>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_150>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SCALAR_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_200>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                              

                                              
template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_200>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                             


template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_200>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                             

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_200>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_200>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_200>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_200>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);  

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SCALAR_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_250>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                              

                                              
template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_250>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                             


template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_250>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                             

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_250>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_250>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_250>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH,
                                              gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstream_max_signal_components::MAX_SIGNAL_COMPONENTS_250>(std::int32_t * __restrict__,
                                                                                                                   std::int32_t * __restrict__);                                            



                                              


template<gms::radiolocation::iq_rectw_bitstream_vsequence_t
         ::IQ_rectw_bitstreams_optim_path  optim_path>        
std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u10x(std::int32_t * ret_codes)
{
       if(__builtin_expect(nullptr==ret_codes,0)) { return (-1);}
       if(__builtin_expect(this->m_IQ_nsignals>256ull,0)) { return (-2);}
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
           PMC_VARS
           HW_PMC_COLLECTION_PROLOGE_BODY
#endif 

        constexpr std::size_t nsignals_size{256ull};
        __ATTR_ALIGN__(16) 
        __m128 xmm0_i[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm1_i[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm2_i[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm3_i[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm4_i[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm5_i[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm6_i[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm7_i[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm8_i[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm9_i[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm0_q[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm1_q[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm2_q[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm3_q[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm4_q[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm5_q[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm6_q[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm7_q[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm8_q[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm9_q[nsignals_size];
        float scal_rem_i[nsignals_size];
        float scal_rem_q[nsignals_size];
        __m128 vzero{_mm_setzero_ps()};
        __m128 sum0_i{vzero};
        __m128 sum1_i{sum0_i};
        __m128 sum2_i{sum0_i};
        __m128 sum3_i{sum0_i};
        __m128 sum4_i{sum0_i};
        __m128 sum5_i{sum0_i};
        __m128 sum6_i{sum0_i};
        __m128 sum7_i{sum0_i};
        __m128 sum8_i{sum0_i};
        __m128 sum9_i{sum0_i};
        __m128 sum0_q{sum0_i};
        __m128 sum1_q{sum0_i};
        __m128 sum2_q{sum0_i};
        __m128 sum3_q{sum0_i};
        __m128 sum4_q{sum0_i};
        __m128 sum5_q{sum0_i};
        __m128 sum6_q{sum0_i};
        __m128 sum7_q{sum0_i};
        __m128 sum8_q{sum0_i};
        __m128 sum9_q{sum0_i};
        std::size_t i;
        float sum_rem_i{0.0f};
        float sum_rem_q{0.0f};

        if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SCALAR_PATH)
        {          
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_scalar();
                }            
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_sse();
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_sse_u4x();
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx_u4x(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx_u4x(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx512_u4x(true);
                } 
        }

        if(this->m_I_nsamples==this->m_Q_nsamples)
        {
              for(i = 0ull; (i+39ull) < this->m_I_nsamples; i += 40ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    sum6_i = vzero;
                    sum6_q = vzero;
                    sum7_i = vzero;
                    sum7_q = vzero;
                    sum8_i = vzero;
                    sum8_q = vzero;
                    sum9_i = vzero;
                    sum9_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);
#endif 
                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);
                        xmm2_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum2_q    = _mm_add_ps(xmm2_q[j],sum2_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum2_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum2_q);
                        xmm3_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+12ull]);
                        sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+12ull]);
                        sum3_q    = _mm_add_ps(xmm3_q[j],sum3_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+12ull],sum3_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+12ull],sum3_q);
                        xmm4_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum4_q    = _mm_add_ps(xmm4_q[j],sum4_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum4_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum4_q);
                        xmm5_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+20ull]);
                        sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+20ull]);
                        sum5_q    = _mm_add_ps(xmm5_q[j],sum5_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+20ull],sum5_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+20ull],sum5_q);
                        xmm6_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum6_i    = _mm_add_ps(xmm6_i[j],sum6_i);
                        xmm6_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum6_q    = _mm_add_ps(xmm6_q[j],sum6_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum6_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum6_q);
                        xmm7_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+28ull]);
                        sum7_i    = _mm_add_ps(xmm7_i[j],sum7_i);
                        xmm7_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+28ull]);
                        sum7_q    = _mm_add_ps(xmm7_q[j],sum7_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+28ull],sum7_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+28ull],sum7_q);
                        xmm8_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum8_i    = _mm_add_ps(xmm8_i[j],sum8_i);
                        xmm8_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum8_q    = _mm_add_ps(xmm8_q[j],sum8_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum8_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum8_q);
                        xmm9_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+36ull]);
                        sum9_i    = _mm_add_ps(xmm9_i[j],sum9_i);
                        xmm9_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+36ull]);
                        sum9_q    = _mm_add_ps(xmm9_q[j],sum9_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+36ull],sum9_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+36ull],sum9_q);
                                              
                    }
              }

              for(; (i+31ull) < this->m_I_nsamples; i += 32ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    sum6_i = vzero;
                    sum6_q = vzero;
                    sum7_i = vzero;
                    sum7_q = vzero;
                                       
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);
                        xmm2_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum2_q    = _mm_add_ps(xmm2_q[j],sum2_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum2_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum2_q);
                        xmm3_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+12ull]);
                        sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+12ull]);
                        sum3_q    = _mm_add_ps(xmm3_q[j],sum3_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+12ull],sum3_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+12ull],sum3_q);
                        xmm4_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum4_q    = _mm_add_ps(xmm4_q[j],sum4_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum4_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum4_q);
                        xmm5_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+20ull]);
                        sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+20ull]);
                        sum5_q    = _mm_add_ps(xmm5_q[j],sum5_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+20ull],sum5_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+20ull],sum5_q);
                        xmm6_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum6_i    = _mm_add_ps(xmm6_i[j],sum6_i);
                        xmm6_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum6_q    = _mm_add_ps(xmm6_q[j],sum6_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum6_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum6_q);
                        xmm7_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+28ull]);
                        sum7_i    = _mm_add_ps(xmm7_i[j],sum7_i);
                        xmm7_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+28ull]);
                        sum7_q    = _mm_add_ps(xmm7_q[j],sum7_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+28ull],sum7_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+28ull],sum7_q);
                                                                      
                    }
              }

               for(; (i+15ull) < this->m_I_nsamples; i += 16ull) 
               {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                                      
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);
                        xmm2_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum2_q    = _mm_add_ps(xmm2_q[j],sum2_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum2_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum2_q);
                        xmm3_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+12ull]);
                        sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+12ull]);
                        sum3_q    = _mm_add_ps(xmm3_q[j],sum3_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+12ull],sum3_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+12ull],sum3_q);
                       
                    }
              }

              for(; (i+7ull) < this->m_I_nsamples; i += 8ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                                      
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);
                                                                     
                    }
              }

              for(; (i+3ull) < this->m_I_nsamples; i += 4ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                                                          
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                                                                                             
                    }
              }

              for(; (i+0ull) < this->m_I_nsamples; i += 1ull) 
              {
                   sum_rem_i = 0.0f;
                   sum_rem_q = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_i[j] = this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i];
                       sum_rem_i     += scal_rem_i[j];
                       this->m_I_vsequence.m_data[i] = sum_rem_i;
                       scal_rem_q[j] = this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i];
                       sum_rem_q     += scal_rem_q[j];
                       this->m_Q_vsequence.m_data[i] = sum_rem_q;
                   }
              }
        }
        else 
        {
             for(i = 0ull; (i+39ull) < this->m_I_nsamples; i += 40ull) 
             {
                    sum0_i = vzero;
                    sum1_i = vzero;
                    sum2_i = vzero; 
                    sum3_i = vzero;                   
                    sum4_i = vzero;                    
                    sum5_i = vzero;                   
                    sum6_i = vzero;                    
                    sum7_i = vzero;                    
                    sum8_i = vzero;                    
                    sum9_i = vzero;                  
                                        
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);                      
#endif 
                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);                       
                        xmm2_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i); 
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum2_i);                        
                        xmm3_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+12ull]);
                        sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+12ull],sum3_i);                        
                        xmm4_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum4_i);                        
                        xmm5_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+20ull]);
                        sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+20ull],sum5_i);                       
                        xmm6_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum6_i    = _mm_add_ps(xmm6_i[j],sum6_i);                      
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum6_i);                       
                        xmm7_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+28ull]);
                        sum7_i    = _mm_add_ps(xmm7_i[j],sum7_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+28ull],sum7_i);                       
                        xmm8_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum8_i    = _mm_add_ps(xmm8_i[j],sum8_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum8_i);                       
                        xmm9_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+36ull]);
                        sum9_i    = _mm_add_ps(xmm9_i[j],sum9_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+36ull],sum9_i);                        
                                                                       
                    }
              }
              
              for(; (i+31ull) < this->m_I_nsamples; i += 32ull) 
              {
                    sum0_i = vzero;
                    sum1_i = vzero;
                    sum2_i = vzero; 
                    sum3_i = vzero;                   
                    sum4_i = vzero;                    
                    sum5_i = vzero;                   
                    sum6_i = vzero;                    
                    sum7_i = vzero;                    
                                       
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);                       
                        xmm2_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i); 
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum2_i);                        
                        xmm3_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+12ull]);
                        sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+12ull],sum3_i);                        
                        xmm4_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum4_i);                        
                        xmm5_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+20ull]);
                        sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+20ull],sum5_i);                       
                        xmm6_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum6_i    = _mm_add_ps(xmm6_i[j],sum6_i);                      
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum6_i);                       
                        xmm7_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+28ull]);
                        sum7_i    = _mm_add_ps(xmm7_i[j],sum7_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+28ull],sum7_i);                       
                                                                      
                    }
              }

               for(; (i+15ull) < this->m_I_nsamples; i += 16ull) 
               {
                    sum0_i = vzero;
                    sum1_i = vzero;
                    sum2_i = vzero; 
                    sum3_i = vzero;                   
                                        
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);                       
                        xmm2_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i); 
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum2_i);                        
                        xmm3_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+12ull]);
                        sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+12ull],sum3_i);                        
                                                                      
                    }
              }

               for(; (i+7ull) < this->m_I_nsamples; i += 8ull) 
               {
                    sum0_i = vzero;
                    sum1_i = vzero;                 
                                        
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);                                        
                                                                      
                    }
              }

               for(; (i+3ull) < this->m_I_nsamples; i += 4ull) 
               {
                    sum0_i = vzero;    
                                        
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                                      
                                                                      
                    }
              }
              
              for(; (i+0ull) < this->m_I_nsamples; i += 1ull) 
              {
                   sum_rem_i = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_i[j] = this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i];
                       sum_rem_i     += scal_rem_i[j];
                       this->m_I_vsequence.m_data[i] = sum_rem_i;
                       
                   }
              }

              ///////////////////////////// Q-channel //////////////////////////////////////////////
              
               for(i = 0ull; (i+39ull) < this->m_Q_nsamples; i += 40ull) 
               {
                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                    
                    sum6_q = vzero;                    
                    sum7_q = vzero;                    
                    sum8_q = vzero;                    
                    sum9_q = vzero;                    
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);
#endif 
                       
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);                      
                        xmm2_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum2_q    = _mm_add_ps(xmm2_q[j],sum2_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum2_q);                       
                        xmm3_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+12ull]);
                        sum3_q    = _mm_add_ps(xmm3_q[j],sum3_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+12ull],sum3_q);                      
                        xmm4_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum4_q    = _mm_add_ps(xmm4_q[j],sum4_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum4_q);                        
                        xmm5_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+20ull]);
                        sum5_q    = _mm_add_ps(xmm5_q[j],sum5_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+20ull],sum5_q);                        
                        xmm6_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum6_q    = _mm_add_ps(xmm6_q[j],sum6_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum6_q);                      
                        xmm7_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+28ull]);
                        sum7_q    = _mm_add_ps(xmm7_q[j],sum7_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+28ull],sum7_q);                       
                        xmm8_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum8_q    = _mm_add_ps(xmm8_q[j],sum8_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum8_q);                       
                        xmm9_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+36ull]);
                        sum9_q    = _mm_add_ps(xmm9_q[j],sum9_q);                      
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+36ull],sum9_q);                       
                                              
                    }
              }
                
              for(; (i+31ull) < this->m_Q_nsamples; i += 32ull) 
              {
                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                    
                    sum6_q = vzero;                    
                    sum7_q = vzero;                                                   
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {                       
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);                      
                        xmm2_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum2_q    = _mm_add_ps(xmm2_q[j],sum2_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum2_q);                       
                        xmm3_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+12ull]);
                        sum3_q    = _mm_add_ps(xmm3_q[j],sum3_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+12ull],sum3_q);                      
                        xmm4_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum4_q    = _mm_add_ps(xmm4_q[j],sum4_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum4_q);                        
                        xmm5_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+20ull]);
                        sum5_q    = _mm_add_ps(xmm5_q[j],sum5_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+20ull],sum5_q);                        
                        xmm6_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum6_q    = _mm_add_ps(xmm6_q[j],sum6_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum6_q);                      
                        xmm7_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+28ull]);
                        sum7_q    = _mm_add_ps(xmm7_q[j],sum7_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+28ull],sum7_q);                       
                                                                                      
                    }
              }

              for(; (i+15ull) < this->m_Q_nsamples; i += 16ull) 
              {
                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                                                
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {                       
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);                      
                        xmm2_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum2_q    = _mm_add_ps(xmm2_q[j],sum2_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum2_q);                       
                        xmm3_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+12ull]);
                        sum3_q    = _mm_add_ps(xmm3_q[j],sum3_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+12ull],sum3_q);                                      
                                              
                    }
              }

              for(; (i+7ull) < this->m_Q_nsamples; i += 8ull) 
              {                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                                                      
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                       
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);                                            
                                              
                    }
              }

              for(; (i+3ull) < this->m_Q_nsamples; i += 4ull) 
              {                    
                    sum0_q = vzero;                                                       
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {                       
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                                             
                                              
                    }
              }

               for(; (i+0ull) < this->m_Q_nsamples; i += 1ull) 
               {                 
                   sum_rem_q = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_q[j] = this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i];
                       sum_rem_q     += scal_rem_q[j];
                       this->m_Q_vsequence.m_data[i] = sum_rem_q;
                   }
              }
        }
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
           HW_PMC_COLLECTION_EPILOGE_BODY
           HW_PMC_SHOW_INFO
#endif 
     return (0);
}

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SCALAR_PATH>(std::int32_t * __restrict__);
                                              
template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH>(std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH>(std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH>(std::int32_t * __restrict__);
                                             

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH>(std::int32_t * __restrict__);
                                              


template<gms::radiolocation::iq_rectw_bitstream_vsequence_t
         ::IQ_rectw_bitstreams_optim_path  optim_path>         
std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u6x(std::int32_t * ret_codes)
{
       if(__builtin_expect(nullptr==ret_codes,0)) { return (-1);}
       if(__builtin_expect(this->m_IQ_nsignals>256ull,0)) { return (-2);}
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
           PMC_VARS
           HW_PMC_COLLECTION_PROLOGE_BODY
#endif 
       
       constexpr std::size_t nsignals_size{256ull};
        __ATTR_ALIGN__(16) 
        __m128 xmm0_i[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm1_i[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm2_i[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm3_i[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm4_i[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm5_i[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm0_q[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm1_q[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm2_q[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm3_q[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm4_q[nsignals_size];
        __ATTR_ALIGN__(16) 
        __m128 xmm5_q[nsignals_size];
        float scal_rem_i[nsignals_size];
        float scal_rem_q[nsignals_size];
        __m128 vzero{_mm_setzero_ps()};
        __m128 sum0_i{vzero};
        __m128 sum1_i{sum0_i};
        __m128 sum2_i{sum0_i};
        __m128 sum3_i{sum0_i};
        __m128 sum4_i{sum0_i};
        __m128 sum5_i{sum0_i};
        __m128 sum0_q{sum0_i};
        __m128 sum1_q{sum0_i};
        __m128 sum2_q{sum0_i};
        __m128 sum3_q{sum0_i};
        __m128 sum4_q{sum0_i};
        __m128 sum5_q{sum0_i};
        std::size_t i;
        float sum_rem_i{0.0f};
        float sum_rem_q{0.0f};

        if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SCALAR_PATH)
        {          
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_scalar();
                }            
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_sse();
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_sse_u4x();
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx_u4x(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx_u4x(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx512_u4x(true);
                } 
        }

        if(this->m_I_nsamples==this->m_Q_nsamples)
        {
              for(i = 0ull; (i+23ull) < this->m_I_nsamples; i += 24ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);
#endif 
                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);
                        xmm2_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum2_q    = _mm_add_ps(xmm2_q[j],sum2_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum2_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum2_q);
                        xmm3_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+12ull]);
                        sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+12ull]);
                        sum3_q    = _mm_add_ps(xmm3_q[j],sum3_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+12ull],sum3_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+12ull],sum3_q);
                        xmm4_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum4_q    = _mm_add_ps(xmm4_q[j],sum4_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum4_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum4_q);
                        xmm5_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+20ull]);
                        sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+20ull]);
                        sum5_q    = _mm_add_ps(xmm5_q[j],sum5_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+20ull],sum5_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+20ull],sum5_q);
                                                                     
                    }
              }

               for(; (i+15ull) < this->m_I_nsamples; i += 16ull) 
               {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                                      
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);
                        xmm2_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum2_q    = _mm_add_ps(xmm2_q[j],sum2_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum2_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum2_q);
                        xmm3_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+12ull]);
                        sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+12ull]);
                        sum3_q    = _mm_add_ps(xmm3_q[j],sum3_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+12ull],sum3_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+12ull],sum3_q);
                       
                    }
              }

              for(; (i+7ull) < this->m_I_nsamples; i += 8ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                                      
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);
                                                                     
                    }
              }

              for(; (i+3ull) < this->m_I_nsamples; i += 4ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                                                          
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                                                                                             
                    }
              }

              for(; (i+0ull) < this->m_I_nsamples; i += 1ull) 
              {
                   sum_rem_i = 0.0f;
                   sum_rem_q = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_i[j] = this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i];
                       sum_rem_i     += scal_rem_i[j];
                       this->m_I_vsequence.m_data[i] = sum_rem_i;
                       scal_rem_q[j] = this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i];
                       sum_rem_q     += scal_rem_q[j];
                       this->m_Q_vsequence.m_data[i] = sum_rem_q;
                   }
              }
        }
        else 
        {
             for(i = 0ull; (i+23ull) < this->m_I_nsamples; i += 24ull) 
             {
                    sum0_i = vzero;
                    sum1_i = vzero;
                    sum2_i = vzero; 
                    sum3_i = vzero;                   
                    sum4_i = vzero;                    
                    sum5_i = vzero;                   
                                                      
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);                       
                                              
#endif 
                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);                       
                        xmm2_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i); 
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum2_i);                        
                        xmm3_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+12ull]);
                        sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+12ull],sum3_i);                        
                        xmm4_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum4_i    = _mm_add_ps(xmm4_i[j],sum4_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum4_i);                        
                        xmm5_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+20ull]);
                        sum5_i    = _mm_add_ps(xmm5_i[j],sum5_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+20ull],sum5_i);                                          
                                                                       
                    }
              }
              
              for(; (i+15ull) < this->m_I_nsamples; i += 16ull) 
              {
                    sum0_i = vzero;
                    sum1_i = vzero;
                    sum2_i = vzero; 
                    sum3_i = vzero;                   
                                        
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);                       
                        xmm2_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum2_i    = _mm_add_ps(xmm2_i[j],sum2_i); 
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum2_i);                        
                        xmm3_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+12ull]);
                        sum3_i    = _mm_add_ps(xmm3_i[j],sum3_i);                        
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+12ull],sum3_i);                        
                                                                      
                    }
              }

               for(; (i+7ull) < this->m_I_nsamples; i += 8ull) 
               {
                    sum0_i = vzero;
                    sum1_i = vzero;                 
                                        
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+4ull]);
                        sum1_i    = _mm_add_ps(xmm1_i[j],sum1_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+4ull],sum1_i);                                        
                                                                      
                    }
              }

               for(; (i+3ull) < this->m_I_nsamples; i += 4ull) 
               {
                    sum0_i = vzero;    
                                        
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm_add_ps(xmm0_i[j],sum0_i);                       
                        _mm_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                                      
                                                                      
                    }
              }
              
              for(; (i+0ull) < this->m_I_nsamples; i += 1ull) 
              {
                   sum_rem_i = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_i[j] = this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i];
                       sum_rem_i     += scal_rem_i[j];
                       this->m_I_vsequence.m_data[i] = sum_rem_i;
                       
                   }
              }

              ///////////////////////////// Q-channel //////////////////////////////////////////////
              
               for(i = 0ull; (i+23ull) < this->m_Q_nsamples; i += 24ull) 
               {
                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                    
                                    
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);                        
                        
#endif 
                       
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);                      
                        xmm2_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum2_q    = _mm_add_ps(xmm2_q[j],sum2_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum2_q);                       
                        xmm3_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+12ull]);
                        sum3_q    = _mm_add_ps(xmm3_q[j],sum3_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+12ull],sum3_q);                      
                        xmm4_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum4_q    = _mm_add_ps(xmm4_q[j],sum4_q);                        
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum4_q);                        
                        xmm5_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+20ull]);
                        sum5_q    = _mm_add_ps(xmm5_q[j],sum5_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+20ull],sum5_q);                        
                                                            
                    }
              }
                
              for(; (i+15ull) < this->m_Q_nsamples; i += 16ull) 
              {
                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                                                
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {                       
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);                      
                        xmm2_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum2_q    = _mm_add_ps(xmm2_q[j],sum2_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum2_q);                       
                        xmm3_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+12ull]);
                        sum3_q    = _mm_add_ps(xmm3_q[j],sum3_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+12ull],sum3_q);                                      
                                              
                    }
              }

              for(; (i+7ull) < this->m_Q_nsamples; i += 8ull) 
              {                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                                                      
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                       
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+4ull]);
                        sum1_q    = _mm_add_ps(xmm1_q[j],sum1_q);                       
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+4ull],sum1_q);                                            
                                              
                    }
              }

              for(; (i+3ull) < this->m_Q_nsamples; i += 4ull) 
              {                    
                    sum0_q = vzero;                                                       
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {                       
                        xmm0_q[j] = _mm_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm_add_ps(xmm0_q[j],sum0_q);
                        _mm_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                                             
                                              
                    }
              }

               for(; (i+0ull) < this->m_Q_nsamples; i += 1ull) 
               {                 
                   sum_rem_q = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_q[j] = this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i];
                       sum_rem_q     += scal_rem_q[j];
                       this->m_Q_vsequence.m_data[i] = sum_rem_q;
                   }
              }
        }
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
           HW_PMC_COLLECTION_EPILOGE_BODY
           HW_PMC_SHOW_INFO
#endif 
     return (0);
}

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SCALAR_PATH>(std::int32_t * __restrict__);
                                             

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH>(std::int32_t * __restrict__);
                                             

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH>(std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH>(std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH>(std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH>(std::int32_t * __restrict__);
                                              

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_sse_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);


template<gms::radiolocation::iq_rectw_bitstream_vsequence_t
         ::IQ_rectw_bitstreams_optim_path  optim_path>         
std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u16x(std::int32_t * ret_codes)
{
       if(__builtin_expect(nullptr==ret_codes,0)) { return (-1);}
       if(__builtin_expect(this->m_IQ_nsignals>256ull,0)) { return (-2);}
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
           PMC_VARS
           HW_PMC_COLLECTION_PROLOGE_BODY
#endif 

        constexpr std::size_t nsignals_size{256ull};
        __ATTR_ALIGN__(32) 
        __m256 xmm0_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm1_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm2_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm3_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm4_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm5_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm6_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm7_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm8_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm9_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm10_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm11_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm12_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm13_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm14_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm15_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm0_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm1_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm2_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm3_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm4_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm5_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm6_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm7_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm8_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm9_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm10_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm11_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm12_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm13_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm14_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm15_q[nsignals_size];
        float scal_rem_i[nsignals_size];
        float scal_rem_q[nsignals_size];
        __m256 vzero{_mm256_setzero_ps()};
        __m256 sum0_i{vzero};
        __m256 sum1_i{sum0_i};
        __m256 sum2_i{sum0_i};
        __m256 sum3_i{sum0_i};
        __m256 sum4_i{sum0_i};
        __m256 sum5_i{sum0_i};
        __m256 sum6_i{sum0_i};
        __m256 sum7_i{sum0_i};
        __m256 sum8_i{sum0_i};
        __m256 sum9_i{sum0_i};
        __m256 sum10_i{sum0_i};
        __m256 sum11_i{sum0_i};
        __m256 sum12_i{sum0_i};
        __m256 sum13_i{sum0_i};
        __m256 sum14_i{sum0_i};
        __m256 sum15_i{sum0_i};
        __m256 sum0_q{sum0_i};
        __m256 sum1_q{sum0_i};
        __m256 sum2_q{sum0_i};
        __m256 sum3_q{sum0_i};
        __m256 sum4_q{sum0_i};
        __m256 sum5_q{sum0_i};
        __m256 sum6_q{sum0_i};
        __m256 sum7_q{sum0_i};
        __m256 sum8_q{sum0_i};
        __m256 sum9_q{sum0_i};
        __m256 sum10_q{sum0_i};
        __m256 sum11_q{sum0_i};
        __m256 sum12_q{sum0_i};
        __m256 sum13_q{sum0_i};
        __m256 sum14_q{sum0_i};
        __m256 sum15_q{sum0_i}; 
        std::size_t i;
        float sum_rem_i{0.0f};
        float sum_rem_q{0.0f};

        if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SCALAR_PATH)
        {          
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_scalar();
                }            
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_sse();
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_sse_u4x();
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx_u4x(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx_u4x(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx512_u4x(true);
                } 
        }

        if(this->m_I_nsamples==this->m_Q_nsamples)
        {
              for(i = 0ull; (i+127ull) < this->m_I_nsamples; i += 128ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    sum6_i = vzero;
                    sum6_q = vzero;
                    sum7_i = vzero;
                    sum7_q = vzero;
                    sum8_i = vzero;
                    sum8_q = vzero;
                    sum9_i = vzero;
                    sum9_q = vzero;
                    sum10_i = vzero;
                    sum10_q = vzero;
                    sum11_i = vzero;
                    sum11_q = vzero;
                    sum12_i = vzero;
                    sum12_q = vzero;
                    sum13_i = vzero;
                    sum13_q = vzero;
                    sum14_i = vzero;
                    sum14_q = vzero;
                    sum15_i = vzero;
                    sum15_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull],_MM_HINT_T1);
#endif 
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);
                        xmm4_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum4_i    = _mm256_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum4_q    = _mm256_add_ps(xmm4_q[j],sum4_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum4_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum4_q);
                        xmm5_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+40ull]);
                        sum5_i    = _mm256_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+40ull]);
                        sum5_q    = _mm256_add_ps(xmm5_q[j],sum5_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+40ull],sum5_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+40ull],sum5_q);
                        xmm6_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum6_i    = _mm256_add_ps(xmm6_i[j],sum6_i);
                        xmm6_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum6_q    = _mm256_add_ps(xmm6_q[j],sum6_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum6_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum6_q);
                        xmm7_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+56ull]);
                        sum7_i    = _mm256_add_ps(xmm7_i[j],sum7_i);
                        xmm7_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+56ull]);
                        sum7_q    = _mm256_add_ps(xmm7_q[j],sum7_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+56ull],sum7_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+56ull],sum7_q);
                        xmm8_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum8_i    = _mm256_add_ps(xmm8_i[j],sum8_i);
                        xmm8_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum8_q    = _mm256_add_ps(xmm8_q[j],sum8_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum8_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum8_q);
                        xmm9_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+72ull]);
                        sum9_i    = _mm256_add_ps(xmm9_i[j],sum9_i);
                        xmm9_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+72ull]);
                        sum9_q    = _mm256_add_ps(xmm9_q[j],sum9_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+72ull],sum9_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+72ull],sum9_q);
                        xmm10_i[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull]);
                        sum10_i   = _mm256_add_ps(xmm10_i[j],sum10_i);
                        xmm10_q[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull]);
                        sum10_q   = _mm256_add_ps(xmm10_q[j],sum10_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+80ull],sum10_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+80ull],sum10_q);
                        xmm11_i[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+88ull]);
                        sum11_i   = _mm256_add_ps(xmm11_i[j],sum11_i);
                        xmm11_q[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+88ull]);
                        sum11_q   = _mm256_add_ps(xmm11_q[j],sum11_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+88ull],sum11_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+88ull],sum11_q);
                        xmm12_i[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull]);
                        sum12_i   = _mm256_add_ps(xmm12_i[j],sum12_i);
                        xmm12_q[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull]);
                        sum12_q   = _mm256_add_ps(xmm12_q[j],sum12_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+96ull],sum12_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+96ull],sum12_q);
                        xmm13_i[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+104ull]);
                        sum13_i   = _mm256_add_ps(xmm13_i[j],sum13_i);
                        xmm13_q[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+104ull]);
                        sum13_q   = _mm256_add_ps(xmm13_q[j],sum13_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+104ull],sum13_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+104ull],sum13_q);
                        xmm14_i[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull]);
                        sum14_i   = _mm256_add_ps(xmm14_i[j],sum14_i);
                        xmm14_q[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull]);
                        sum14_q   = _mm256_add_ps(xmm14_q[j],sum14_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+112ull],sum14_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+112ull],sum14_q);
                        xmm15_i[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+120ull]);
                        sum15_i   = _mm256_add_ps(xmm15_i[j],sum15_i);
                        xmm15_q[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+120ull]);
                        sum15_q   = _mm256_add_ps(xmm15_q[j],sum15_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+120ull],sum15_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+120ull],sum15_q);
                       
                    }
              }
              
              for(; (i+95ull) < this->m_I_nsamples; i += 96ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    sum6_i = vzero;
                    sum6_q = vzero;
                    sum7_i = vzero;
                    sum7_q = vzero;
                    sum8_i = vzero;
                    sum8_q = vzero;
                    sum9_i = vzero;
                    sum9_q = vzero;
                    sum10_i = vzero;
                    sum10_q = vzero;
                    sum11_i = vzero;
                    sum11_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);
                        xmm4_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum4_i    = _mm256_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum4_q    = _mm256_add_ps(xmm4_q[j],sum4_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum4_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum4_q);
                        xmm5_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+40ull]);
                        sum5_i    = _mm256_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+40ull]);
                        sum5_q    = _mm256_add_ps(xmm5_q[j],sum5_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+40ull],sum5_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+40ull],sum5_q);
                        xmm6_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum6_i    = _mm256_add_ps(xmm6_i[j],sum6_i);
                        xmm6_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum6_q    = _mm256_add_ps(xmm6_q[j],sum6_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum6_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum6_q);
                        xmm7_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+56ull]);
                        sum7_i    = _mm256_add_ps(xmm7_i[j],sum7_i);
                        xmm7_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+56ull]);
                        sum7_q    = _mm256_add_ps(xmm7_q[j],sum7_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+56ull],sum7_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+56ull],sum7_q);
                        xmm8_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum8_i    = _mm256_add_ps(xmm8_i[j],sum8_i);
                        xmm8_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum8_q    = _mm256_add_ps(xmm8_q[j],sum8_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum8_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum8_q);
                        xmm9_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+72ull]);
                        sum9_i    = _mm256_add_ps(xmm9_i[j],sum9_i);
                        xmm9_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+72ull]);
                        sum9_q    = _mm256_add_ps(xmm9_q[j],sum9_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+72ull],sum9_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+72ull],sum9_q);
                        xmm10_i[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull]);
                        sum10_i   = _mm256_add_ps(xmm10_i[j],sum10_i);
                        xmm10_q[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull]);
                        sum10_q   = _mm256_add_ps(xmm10_q[j],sum10_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+80ull],sum10_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+80ull],sum10_q);
                        xmm11_i[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+88ull]);
                        sum11_i   = _mm256_add_ps(xmm11_i[j],sum11_i);
                        xmm11_q[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+88ull]);
                        sum11_q   = _mm256_add_ps(xmm11_q[j],sum11_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+88ull],sum11_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+88ull],sum11_q);
                                              
                    }
              }
             
              for(; (i+79ull) < this->m_I_nsamples; i += 80ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    sum6_i = vzero;
                    sum6_q = vzero;
                    sum7_i = vzero;
                    sum7_q = vzero;
                    sum8_i = vzero;
                    sum8_q = vzero;
                    sum9_i = vzero;
                    sum9_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);
                        xmm4_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum4_i    = _mm256_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum4_q    = _mm256_add_ps(xmm4_q[j],sum4_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum4_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum4_q);
                        xmm5_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+40ull]);
                        sum5_i    = _mm256_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+40ull]);
                        sum5_q    = _mm256_add_ps(xmm5_q[j],sum5_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+40ull],sum5_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+40ull],sum5_q);
                        xmm6_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum6_i    = _mm256_add_ps(xmm6_i[j],sum6_i);
                        xmm6_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum6_q    = _mm256_add_ps(xmm6_q[j],sum6_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum6_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum6_q);
                        xmm7_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+56ull]);
                        sum7_i    = _mm256_add_ps(xmm7_i[j],sum7_i);
                        xmm7_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+56ull]);
                        sum7_q    = _mm256_add_ps(xmm7_q[j],sum7_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+56ull],sum7_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+56ull],sum7_q);
                        xmm8_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum8_i    = _mm256_add_ps(xmm8_i[j],sum8_i);
                        xmm8_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum8_q    = _mm256_add_ps(xmm8_q[j],sum8_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum8_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum8_q);
                        xmm9_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+72ull]);
                        sum9_i    = _mm256_add_ps(xmm9_i[j],sum9_i);
                        xmm9_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+72ull]);
                        sum9_q    = _mm256_add_ps(xmm9_q[j],sum9_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+72ull],sum9_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+72ull],sum9_q);
                                              
                    }
              }
              
              for(; (i+63ull) < this->m_I_nsamples; i += 64ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    sum6_i = vzero;
                    sum6_q = vzero;
                    sum7_i = vzero;
                    sum7_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);
                        xmm4_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum4_i    = _mm256_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum4_q    = _mm256_add_ps(xmm4_q[j],sum4_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum4_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum4_q);
                        xmm5_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+40ull]);
                        sum5_i    = _mm256_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+40ull]);
                        sum5_q    = _mm256_add_ps(xmm5_q[j],sum5_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+40ull],sum5_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+40ull],sum5_q);
                        xmm6_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum6_i    = _mm256_add_ps(xmm6_i[j],sum6_i);
                        xmm6_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum6_q    = _mm256_add_ps(xmm6_q[j],sum6_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum6_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum6_q);
                        xmm7_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+56ull]);
                        sum7_i    = _mm256_add_ps(xmm7_i[j],sum7_i);
                        xmm7_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+56ull]);
                        sum7_q    = _mm256_add_ps(xmm7_q[j],sum7_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+56ull],sum7_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+56ull],sum7_q);
                                                                     
                    }
              }
               
             for(; (i+31ull) < this->m_I_nsamples; i += 32ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);
                                                                      
                    }
              }
              
              for(; (i+15ull) < this->m_I_nsamples; i += 16ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);
                                                                     
                    }
              }
                            
             for(; (i+7ull) < this->m_I_nsamples; i += 8ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                                                                                            
                    }
              }

              for(; (i+0ull) < this->m_I_nsamples; i += 1ull) 
              {
                   sum_rem_i = 0.0f;
                   sum_rem_q = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_i[j] = this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i];
                       sum_rem_i     += scal_rem_i[j];
                       this->m_I_vsequence.m_data[i] = sum_rem_i;
                       scal_rem_q[j] = this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i];
                       sum_rem_q     += scal_rem_q[j];
                       this->m_Q_vsequence.m_data[i] = sum_rem_q;
                   }
              }

        }
        else 
        {
              for(i = 0ull; (i+127ull) < this->m_I_nsamples; i += 128ull) 
              {
                    sum0_i = vzero;                
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                   
                    sum4_i = vzero;                   
                    sum5_i = vzero;                    
                    sum6_i = vzero;                  
                    sum7_i = vzero;                   
                    sum8_i = vzero;                    
                    sum9_i = vzero;                   
                    sum10_i = vzero;                    
                    sum11_i = vzero;                    
                    sum12_i = vzero;                  
                    sum13_i = vzero;                    
                    sum14_i = vzero;                   
                    sum15_i = vzero;                  
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull],_MM_HINT_T1);                      
#endif 
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);                       
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);                       
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);                        
                        xmm4_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum4_i    = _mm256_add_ps(xmm4_i[j],sum4_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum4_i);                        
                        xmm5_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+40ull]);
                        sum5_i    = _mm256_add_ps(xmm5_i[j],sum5_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+40ull],sum5_i);                        
                        xmm6_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum6_i    = _mm256_add_ps(xmm6_i[j],sum6_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum6_i);                       
                        xmm7_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+56ull]);
                        sum7_i    = _mm256_add_ps(xmm7_i[j],sum7_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+56ull],sum7_i);                        
                        xmm8_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum8_i    = _mm256_add_ps(xmm8_i[j],sum8_i);                      
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum8_i);                       
                        xmm9_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+72ull]);
                        sum9_i    = _mm256_add_ps(xmm9_i[j],sum9_i);                     
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+72ull],sum9_i);                        
                        xmm10_i[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull]);
                        sum10_i   = _mm256_add_ps(xmm10_i[j],sum10_i);                      
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+80ull],sum10_i);                       
                        xmm11_i[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+88ull]);
                        sum11_i   = _mm256_add_ps(xmm11_i[j],sum11_i);                    
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+88ull],sum11_i);                      
                        xmm12_i[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull]);
                        sum12_i   = _mm256_add_ps(xmm12_i[j],sum12_i);                      
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+96ull],sum12_i);                       
                        xmm13_i[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+104ull]);
                        sum13_i   = _mm256_add_ps(xmm13_i[j],sum13_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+104ull],sum13_i);                       
                        xmm14_i[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull]);
                        sum14_i   = _mm256_add_ps(xmm14_i[j],sum14_i);                     
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+112ull],sum14_i);                       
                        xmm15_i[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+120ull]);
                        sum15_i   = _mm256_add_ps(xmm15_i[j],sum15_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+120ull],sum15_i);                       
                       
                    }
              }

              for(; (i+95ull) < this->m_I_nsamples; i += 96ull) 
              {
                    sum0_i = vzero;                
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                   
                    sum4_i = vzero;                   
                    sum5_i = vzero;                    
                    sum6_i = vzero;                  
                    sum7_i = vzero;                   
                    sum8_i = vzero;                    
                    sum9_i = vzero;                   
                    sum10_i = vzero;                    
                    sum11_i = vzero;                                   
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);                       
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);                       
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);                        
                        xmm4_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum4_i    = _mm256_add_ps(xmm4_i[j],sum4_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum4_i);                        
                        xmm5_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+40ull]);
                        sum5_i    = _mm256_add_ps(xmm5_i[j],sum5_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+40ull],sum5_i);                        
                        xmm6_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum6_i    = _mm256_add_ps(xmm6_i[j],sum6_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum6_i);                       
                        xmm7_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+56ull]);
                        sum7_i    = _mm256_add_ps(xmm7_i[j],sum7_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+56ull],sum7_i);                        
                        xmm8_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum8_i    = _mm256_add_ps(xmm8_i[j],sum8_i);                      
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum8_i);                       
                        xmm9_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+72ull]);
                        sum9_i    = _mm256_add_ps(xmm9_i[j],sum9_i);                     
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+72ull],sum9_i);                        
                        xmm10_i[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull]);
                        sum10_i   = _mm256_add_ps(xmm10_i[j],sum10_i);                      
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+80ull],sum10_i);                       
                        xmm11_i[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+88ull]);
                        sum11_i   = _mm256_add_ps(xmm11_i[j],sum11_i);                    
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+88ull],sum11_i);                                      
                       
                    }
              }

              for(; (i+79ull) < this->m_I_nsamples; i += 80ull) 
              {
                    sum0_i = vzero;                
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                   
                    sum4_i = vzero;                   
                    sum5_i = vzero;                    
                    sum6_i = vzero;                  
                    sum7_i = vzero;                   
                    sum8_i = vzero;                    
                    sum9_i = vzero;                                                
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);                       
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);                       
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);                        
                        xmm4_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum4_i    = _mm256_add_ps(xmm4_i[j],sum4_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum4_i);                        
                        xmm5_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+40ull]);
                        sum5_i    = _mm256_add_ps(xmm5_i[j],sum5_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+40ull],sum5_i);                        
                        xmm6_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum6_i    = _mm256_add_ps(xmm6_i[j],sum6_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum6_i);                       
                        xmm7_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+56ull]);
                        sum7_i    = _mm256_add_ps(xmm7_i[j],sum7_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+56ull],sum7_i);                        
                        xmm8_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum8_i    = _mm256_add_ps(xmm8_i[j],sum8_i);                      
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum8_i);                       
                        xmm9_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+72ull]);
                        sum9_i    = _mm256_add_ps(xmm9_i[j],sum9_i);                     
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+72ull],sum9_i);                                                         
                       
                    }
              }

               for(; (i+63ull) < this->m_I_nsamples; i += 64ull) 
              {
                    sum0_i = vzero;                
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                   
                    sum4_i = vzero;                   
                    sum5_i = vzero;                    
                    sum6_i = vzero;                  
                    sum7_i = vzero;                                                  
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);                       
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);                       
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);                        
                        xmm4_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum4_i    = _mm256_add_ps(xmm4_i[j],sum4_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum4_i);                        
                        xmm5_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+40ull]);
                        sum5_i    = _mm256_add_ps(xmm5_i[j],sum5_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+40ull],sum5_i);                        
                        xmm6_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum6_i    = _mm256_add_ps(xmm6_i[j],sum6_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum6_i);                       
                        xmm7_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+56ull]);
                        sum7_i    = _mm256_add_ps(xmm7_i[j],sum7_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+56ull],sum7_i);                        
                                                      
                    }
              }


              for(; (i+31ull) < this->m_I_nsamples; i += 32ull) 
              {
                    sum0_i = vzero;                
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                                                                
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);                       
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);                       
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);                        
                                                                            
                    }
              }

              for(; (i+15ull) < this->m_I_nsamples; i += 16ull) 
              {
                    sum0_i = vzero;                
                    sum1_i = vzero;                                                                    
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);                       
                                                                           
                    }
              } 
                            
             for(; (i+7ull) < this->m_I_nsamples; i += 8ull) 
             {
                    sum0_i = vzero;                                                                                                     
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                                                                                               
                    }
              } 

              for(; (i+0ull) < this->m_I_nsamples; i += 1ull) 
              {
                   sum_rem_i = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_i[j] = this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i];
                       sum_rem_i     += scal_rem_i[j];
                       this->m_I_vsequence.m_data[i] = sum_rem_i;
                   }
              }

              ///////////////////////////// Q-channel //////////////////////////////////////////////
              for(i = 0ull; (i+127ull) < this->m_Q_nsamples; i += 128ull) 
              {
                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                    
                    sum6_q = vzero;                    
                    sum7_q = vzero;                   
                    sum8_q = vzero;                   
                    sum9_q = vzero;                   
                    sum10_q = vzero;                    
                    sum11_q = vzero;                   
                    sum12_q = vzero;                    
                    sum13_q = vzero;                    
                    sum14_q = vzero;                    
                    sum15_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull],_MM_HINT_T1);                    
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull],_MM_HINT_T1);
#endif 
                        
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                      
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);                       
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);                        
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);                       
                        xmm4_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum4_q    = _mm256_add_ps(xmm4_q[j],sum4_q);                     
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum4_q);                       
                        xmm5_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+40ull]);
                        sum5_q    = _mm256_add_ps(xmm5_q[j],sum5_q);                      
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+40ull],sum5_q);                       
                        xmm6_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum6_q    = _mm256_add_ps(xmm6_q[j],sum6_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum6_q);                       
                        xmm7_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+56ull]);
                        sum7_q    = _mm256_add_ps(xmm7_q[j],sum7_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+56ull],sum7_q);                       
                        xmm8_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum8_q    = _mm256_add_ps(xmm8_q[j],sum8_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum8_q);                        
                        xmm9_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+72ull]);
                        sum9_q    = _mm256_add_ps(xmm9_q[j],sum9_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+72ull],sum9_q);                      
                        xmm10_q[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull]);
                        sum10_q   = _mm256_add_ps(xmm10_q[j],sum10_q);                      
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+80ull],sum10_q);                       
                        xmm11_q[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+88ull]);
                        sum11_q   = _mm256_add_ps(xmm11_q[j],sum11_q);                      
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+88ull],sum11_q);                      
                        xmm12_q[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull]);
                        sum12_q   = _mm256_add_ps(xmm12_q[j],sum12_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+96ull],sum12_q);                       
                        xmm13_q[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+104ull]);
                        sum13_q   = _mm256_add_ps(xmm13_q[j],sum13_q);                      
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+104ull],sum13_q);                     
                        xmm14_q[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull]);
                        sum14_q   = _mm256_add_ps(xmm14_q[j],sum14_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+112ull],sum14_q);                       
                        xmm15_q[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+120ull]);
                        sum15_q   = _mm256_add_ps(xmm15_q[j],sum15_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+120ull],sum15_q);
                       
                    }
              }
              
              for(; (i+95ull) < this->m_Q_nsamples; i += 96ull) 
              {                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                    
                    sum6_q = vzero;                    
                    sum7_q = vzero;                   
                    sum8_q = vzero;                   
                    sum9_q = vzero;                   
                    sum10_q = vzero;                    
                    sum11_q = vzero;                   
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {                        
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                      
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);                       
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);                        
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);                       
                        xmm4_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum4_q    = _mm256_add_ps(xmm4_q[j],sum4_q);                     
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum4_q);                       
                        xmm5_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+40ull]);
                        sum5_q    = _mm256_add_ps(xmm5_q[j],sum5_q);                      
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+40ull],sum5_q);                       
                        xmm6_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum6_q    = _mm256_add_ps(xmm6_q[j],sum6_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum6_q);                       
                        xmm7_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+56ull]);
                        sum7_q    = _mm256_add_ps(xmm7_q[j],sum7_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+56ull],sum7_q);                       
                        xmm8_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum8_q    = _mm256_add_ps(xmm8_q[j],sum8_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum8_q);                        
                        xmm9_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+72ull]);
                        sum9_q    = _mm256_add_ps(xmm9_q[j],sum9_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+72ull],sum9_q);                      
                        xmm10_q[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull]);
                        sum10_q   = _mm256_add_ps(xmm10_q[j],sum10_q);                      
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+80ull],sum10_q);                       
                        xmm11_q[j]= _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+88ull]);
                        sum11_q   = _mm256_add_ps(xmm11_q[j],sum11_q);                      
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+88ull],sum11_q);                      
                       
                    }
              }

              for(; (i+79ull) < this->m_Q_nsamples; i += 80ull) 
              {
                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                    
                    sum6_q = vzero;                    
                    sum7_q = vzero;                   
                    sum8_q = vzero;                   
                    sum9_q = vzero;                                     
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {                        
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                      
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);                       
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);                        
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);                       
                        xmm4_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum4_q    = _mm256_add_ps(xmm4_q[j],sum4_q);                     
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum4_q);                       
                        xmm5_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+40ull]);
                        sum5_q    = _mm256_add_ps(xmm5_q[j],sum5_q);                      
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+40ull],sum5_q);                       
                        xmm6_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum6_q    = _mm256_add_ps(xmm6_q[j],sum6_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum6_q);                       
                        xmm7_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+56ull]);
                        sum7_q    = _mm256_add_ps(xmm7_q[j],sum7_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+56ull],sum7_q);                       
                        xmm8_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum8_q    = _mm256_add_ps(xmm8_q[j],sum8_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum8_q);                        
                        xmm9_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+72ull]);
                        sum9_q    = _mm256_add_ps(xmm9_q[j],sum9_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+72ull],sum9_q);                                           
                    }
              }

              for(; (i+63ull) < this->m_Q_nsamples; i += 64ull) 
              {
                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                    
                    sum6_q = vzero;                    
                    sum7_q = vzero;                                    
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                      
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);                       
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);                        
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);                       
                        xmm4_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum4_q    = _mm256_add_ps(xmm4_q[j],sum4_q);                     
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum4_q);                       
                        xmm5_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+40ull]);
                        sum5_q    = _mm256_add_ps(xmm5_q[j],sum5_q);                      
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+40ull],sum5_q);                       
                        xmm6_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum6_q    = _mm256_add_ps(xmm6_q[j],sum6_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum6_q);                       
                        xmm7_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+56ull]);
                        sum7_q    = _mm256_add_ps(xmm7_q[j],sum7_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+56ull],sum7_q);                       
                    }
              }

              for(; (i+31ull) < this->m_Q_nsamples; i += 32ull) 
              {
                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                                     
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {                        
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                      
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);                       
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);                        
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);                       
                    }
              }

              for(; (i+15ull) < this->m_Q_nsamples; i += 16ull) 
              {
                    sum0_q = vzero;                    
                    sum1_q = vzero;                                    
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                      
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);                       
                    }
              }

              for(; (i+7ull) < this->m_Q_nsamples; i += 8ull) 
              {                    
                    sum0_q = vzero;                                                     
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                                            
                    }
              }
               
               for(; (i+0ull) < this->m_Q_nsamples; i += 1ull) 
               {                 
                   sum_rem_q = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_q[j] = this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i];
                       sum_rem_q     += scal_rem_q[j];
                       this->m_Q_vsequence.m_data[i] = sum_rem_q;
                   }
              }
        }
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
           HW_PMC_COLLECTION_EPILOGE_BODY
           HW_PMC_SHOW_INFO
#endif 
     return (0);
}

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SCALAR_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);


template<gms::radiolocation::iq_rectw_bitstream_vsequence_t
         ::IQ_rectw_bitstreams_optim_path  optim_path>        
std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u10x(std::int32_t * ret_codes)
{
       if(__builtin_expect(nullptr==ret_codes,0)) { return (-1);}
       if(__builtin_expect(this->m_IQ_nsignals>256ull,0)) { return (-2);}
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
           PMC_VARS
           HW_PMC_COLLECTION_PROLOGE_BODY
#endif 

        constexpr std::size_t nsignals_size{256ull};
        __ATTR_ALIGN__(32) 
        __m256 xmm0_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm1_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm2_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm3_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm4_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm5_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm6_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm7_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm8_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm9_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm0_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm1_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm2_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm3_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm4_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm5_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm6_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm7_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm8_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm9_q[nsignals_size];
        float scal_rem_i[nsignals_size];
        float scal_rem_q[nsignals_size];
        __m256 vzero{_mm256_setzero_ps()};
        __m256 sum0_i{vzero};
        __m256 sum1_i{sum0_i};
        __m256 sum2_i{sum0_i};
        __m256 sum3_i{sum0_i};
        __m256 sum4_i{sum0_i};
        __m256 sum5_i{sum0_i};
        __m256 sum6_i{sum0_i};
        __m256 sum7_i{sum0_i};
        __m256 sum8_i{sum0_i};
        __m256 sum9_i{sum0_i};
        __m256 sum0_q{sum0_i};
        __m256 sum1_q{sum0_i};
        __m256 sum2_q{sum0_i};
        __m256 sum3_q{sum0_i};
        __m256 sum4_q{sum0_i};
        __m256 sum5_q{sum0_i};
        __m256 sum6_q{sum0_i};
        __m256 sum7_q{sum0_i};
        __m256 sum8_q{sum0_i};
        __m256 sum9_q{sum0_i}; 
        std::size_t i;
        float sum_rem_i{0.0f};
        float sum_rem_q{0.0f};

        if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SCALAR_PATH)
        {          
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_scalar();
                }            
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_sse();
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_sse_u4x();
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx_u4x(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx_u4x(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx512_u4x(true);
                } 
        }

        if(this->m_I_nsamples==this->m_Q_nsamples)
        {
              for(i = 0ull; (i+79ull) < this->m_I_nsamples; i += 80ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    sum6_i = vzero;
                    sum6_q = vzero;
                    sum7_i = vzero;
                    sum7_q = vzero;
                    sum8_i = vzero;
                    sum8_q = vzero;
                    sum9_i = vzero;
                    sum9_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);
#endif 
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);
                        xmm4_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum4_i    = _mm256_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum4_q    = _mm256_add_ps(xmm4_q[j],sum4_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum4_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum4_q);
                        xmm5_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+40ull]);
                        sum5_i    = _mm256_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+40ull]);
                        sum5_q    = _mm256_add_ps(xmm5_q[j],sum5_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+40ull],sum5_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+40ull],sum5_q);
                        xmm6_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum6_i    = _mm256_add_ps(xmm6_i[j],sum6_i);
                        xmm6_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum6_q    = _mm256_add_ps(xmm6_q[j],sum6_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum6_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum6_q);
                        xmm7_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+56ull]);
                        sum7_i    = _mm256_add_ps(xmm7_i[j],sum7_i);
                        xmm7_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+56ull]);
                        sum7_q    = _mm256_add_ps(xmm7_q[j],sum7_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+56ull],sum7_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+56ull],sum7_q);
                        xmm8_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum8_i    = _mm256_add_ps(xmm8_i[j],sum8_i);
                        xmm8_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum8_q    = _mm256_add_ps(xmm8_q[j],sum8_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum8_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum8_q);
                        xmm9_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+72ull]);
                        sum9_i    = _mm256_add_ps(xmm9_i[j],sum9_i);
                        xmm9_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+72ull]);
                        sum9_q    = _mm256_add_ps(xmm9_q[j],sum9_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+72ull],sum9_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+72ull],sum9_q);
                    }
              }
              
              for(; (i+63ull) < this->m_I_nsamples; i += 64ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    sum6_i = vzero;
                    sum6_q = vzero;
                    sum7_i = vzero;
                    sum7_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);
                        xmm4_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum4_i    = _mm256_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum4_q    = _mm256_add_ps(xmm4_q[j],sum4_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum4_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum4_q);
                        xmm5_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+40ull]);
                        sum5_i    = _mm256_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+40ull]);
                        sum5_q    = _mm256_add_ps(xmm5_q[j],sum5_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+40ull],sum5_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+40ull],sum5_q);
                        xmm6_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum6_i    = _mm256_add_ps(xmm6_i[j],sum6_i);
                        xmm6_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum6_q    = _mm256_add_ps(xmm6_q[j],sum6_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum6_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum6_q);
                        xmm7_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+56ull]);
                        sum7_i    = _mm256_add_ps(xmm7_i[j],sum7_i);
                        xmm7_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+56ull]);
                        sum7_q    = _mm256_add_ps(xmm7_q[j],sum7_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+56ull],sum7_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+56ull],sum7_q);
                                                                     
                    }
              }
               
             for(; (i+31ull) < this->m_I_nsamples; i += 32ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);
                                                                      
                    }
              }
              
              for(; (i+15ull) < this->m_I_nsamples; i += 16ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);
                                                                     
                    }
              }
                            
             for(; (i+7ull) < this->m_I_nsamples; i += 8ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                                                                                            
                    }
              }

              for(; (i+0ull) < this->m_I_nsamples; i += 1ull) 
              {
                   sum_rem_i = 0.0f;
                   sum_rem_q = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_i[j] = this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i];
                       sum_rem_i     += scal_rem_i[j];
                       this->m_I_vsequence.m_data[i] = sum_rem_i;
                       scal_rem_q[j] = this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i];
                       sum_rem_q     += scal_rem_q[j];
                       this->m_Q_vsequence.m_data[i] = sum_rem_q;
                   }
              }

        }
        else 
        {
              for(i = 0ull; (i+79ull) < this->m_I_nsamples; i += 80ull) 
              {
                    sum0_i = vzero;                
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                   
                    sum4_i = vzero;                   
                    sum5_i = vzero;                    
                    sum6_i = vzero;                  
                    sum7_i = vzero;                   
                    sum8_i = vzero;                    
                    sum9_i = vzero;                                    
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);                                          
#endif 
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);                       
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);                       
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);                        
                        xmm4_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum4_i    = _mm256_add_ps(xmm4_i[j],sum4_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum4_i);                        
                        xmm5_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+40ull]);
                        sum5_i    = _mm256_add_ps(xmm5_i[j],sum5_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+40ull],sum5_i);                        
                        xmm6_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum6_i    = _mm256_add_ps(xmm6_i[j],sum6_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum6_i);                       
                        xmm7_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+56ull]);
                        sum7_i    = _mm256_add_ps(xmm7_i[j],sum7_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+56ull],sum7_i);                        
                        xmm8_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum8_i    = _mm256_add_ps(xmm8_i[j],sum8_i);                      
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum8_i);                       
                        xmm9_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+72ull]);
                        sum9_i    = _mm256_add_ps(xmm9_i[j],sum9_i);                     
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+72ull],sum9_i);                               
                    }
              }

               for(; (i+63ull) < this->m_I_nsamples; i += 64ull) 
              {
                    sum0_i = vzero;                
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                   
                    sum4_i = vzero;                   
                    sum5_i = vzero;                    
                    sum6_i = vzero;                  
                    sum7_i = vzero;                                                  
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);                       
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);                       
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);                        
                        xmm4_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum4_i    = _mm256_add_ps(xmm4_i[j],sum4_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum4_i);                        
                        xmm5_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+40ull]);
                        sum5_i    = _mm256_add_ps(xmm5_i[j],sum5_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+40ull],sum5_i);                        
                        xmm6_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum6_i    = _mm256_add_ps(xmm6_i[j],sum6_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum6_i);                       
                        xmm7_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+56ull]);
                        sum7_i    = _mm256_add_ps(xmm7_i[j],sum7_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+56ull],sum7_i);                        
                                                      
                    }
              }

              for(; (i+31ull) < this->m_I_nsamples; i += 32ull) 
              {
                    sum0_i = vzero;                
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                                                                
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);                       
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);                       
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);                        
                                                                            
                    }
              }

              for(; (i+15ull) < this->m_I_nsamples; i += 16ull) 
              {
                    sum0_i = vzero;                
                    sum1_i = vzero;                                                                    
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);                       
                                                                           
                    }
              } 
                            
             for(; (i+7ull) < this->m_I_nsamples; i += 8ull) 
             {
                    sum0_i = vzero;                                                                                                     
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                                                                                               
                    }
              } 

              for(; (i+0ull) < this->m_I_nsamples; i += 1ull) 
              {
                   sum_rem_i = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_i[j] = this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i];
                       sum_rem_i     += scal_rem_i[j];
                       this->m_I_vsequence.m_data[i] = sum_rem_i;
                   }
              }

              ///////////////////////////// Q-channel //////////////////////////////////////////////
              for(i = 0ull; (i+79ull) < this->m_Q_nsamples; i += 80ull) 
              {
                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                    
                    sum6_q = vzero;                    
                    sum7_q = vzero;                   
                    sum8_q = vzero;                   
                    sum9_q = vzero;                   
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);                        
#endif 
                        
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                      
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);                       
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);                        
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);                       
                        xmm4_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum4_q    = _mm256_add_ps(xmm4_q[j],sum4_q);                     
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum4_q);                       
                        xmm5_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+40ull]);
                        sum5_q    = _mm256_add_ps(xmm5_q[j],sum5_q);                      
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+40ull],sum5_q);                       
                        xmm6_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum6_q    = _mm256_add_ps(xmm6_q[j],sum6_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum6_q);                       
                        xmm7_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+56ull]);
                        sum7_q    = _mm256_add_ps(xmm7_q[j],sum7_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+56ull],sum7_q);                       
                        xmm8_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum8_q    = _mm256_add_ps(xmm8_q[j],sum8_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum8_q);                        
                        xmm9_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+72ull]);
                        sum9_q    = _mm256_add_ps(xmm9_q[j],sum9_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+72ull],sum9_q);                      
                    }
              }
              
              for(; (i+63ull) < this->m_Q_nsamples; i += 64ull) 
              {
                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                    
                    sum6_q = vzero;                    
                    sum7_q = vzero;                                    
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                      
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);                       
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);                        
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);                       
                        xmm4_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum4_q    = _mm256_add_ps(xmm4_q[j],sum4_q);                     
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum4_q);                       
                        xmm5_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+40ull]);
                        sum5_q    = _mm256_add_ps(xmm5_q[j],sum5_q);                      
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+40ull],sum5_q);                       
                        xmm6_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum6_q    = _mm256_add_ps(xmm6_q[j],sum6_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum6_q);                       
                        xmm7_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+56ull]);
                        sum7_q    = _mm256_add_ps(xmm7_q[j],sum7_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+56ull],sum7_q);                       
                    }
              }

              for(; (i+31ull) < this->m_Q_nsamples; i += 32ull) 
              {
                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                                     
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {                        
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                      
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);                       
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);                        
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);                       
                    }
              }

              for(; (i+15ull) < this->m_Q_nsamples; i += 16ull) 
              {
                    sum0_q = vzero;                    
                    sum1_q = vzero;                                    
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                      
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);                       
                    }
              }

              for(; (i+7ull) < this->m_Q_nsamples; i += 8ull) 
              {                    
                    sum0_q = vzero;                                                     
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                                            
                    }
              }
               
               for(; (i+0ull) < this->m_Q_nsamples; i += 1ull) 
               {                 
                   sum_rem_q = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_q[j] = this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i];
                       sum_rem_q     += scal_rem_q[j];
                       this->m_Q_vsequence.m_data[i] = sum_rem_q;
                   }
              }
        }
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
           HW_PMC_COLLECTION_EPILOGE_BODY
           HW_PMC_SHOW_INFO
#endif 
     return (0);
}

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SCALAR_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);


template<gms::radiolocation::iq_rectw_bitstream_vsequence_t
         ::IQ_rectw_bitstreams_optim_path  optim_path>         
std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u6x(std::int32_t * ret_codes)
{
       if(__builtin_expect(nullptr==ret_codes,0)) { return (-1);}
       if(__builtin_expect(this->m_IQ_nsignals>256ull,0)) { return (-2);}
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
           PMC_VARS
           HW_PMC_COLLECTION_PROLOGE_BODY
#endif 

        constexpr std::size_t nsignals_size{256ull};
        __ATTR_ALIGN__(32) 
        __m256 xmm0_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm1_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm2_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm3_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm4_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm5_i[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm0_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm1_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm2_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm3_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm4_q[nsignals_size];
        __ATTR_ALIGN__(32) 
        __m256 xmm5_q[nsignals_size];
        float scal_rem_i[nsignals_size];
        float scal_rem_q[nsignals_size];
        __m256 vzero{_mm256_setzero_ps()};
        __m256 sum0_i{vzero};
        __m256 sum1_i{sum0_i};
        __m256 sum2_i{sum0_i};
        __m256 sum3_i{sum0_i};
        __m256 sum4_i{sum0_i};
        __m256 sum5_i{sum0_i};
        __m256 sum0_q{sum0_i};
        __m256 sum1_q{sum0_i};
        __m256 sum2_q{sum0_i};
        __m256 sum3_q{sum0_i};
        __m256 sum4_q{sum0_i};
        __m256 sum5_q{sum0_i};
        std::size_t i;
        float sum_rem_i{0.0f};
        float sum_rem_q{0.0f};

        if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SCALAR_PATH)
        {          
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_scalar();
                }            
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_sse();
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_sse_u4x();
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx_u4x(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx_u4x(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx512_u4x(true);
                } 
        }

        if(this->m_I_nsamples==this->m_Q_nsamples)
        {
              for(i = 0ull; (i+47ull) < this->m_I_nsamples; i += 48ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);
#endif 
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);
                        xmm4_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum4_i    = _mm256_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum4_q    = _mm256_add_ps(xmm4_q[j],sum4_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum4_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum4_q);
                        xmm5_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+40ull]);
                        sum5_i    = _mm256_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+40ull]);
                        sum5_q    = _mm256_add_ps(xmm5_q[j],sum5_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+40ull],sum5_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+40ull],sum5_q);
                   }
              }
              
               
             for(; (i+31ull) < this->m_I_nsamples; i += 32ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);
                                                                      
                    }
              }
              
              for(; (i+15ull) < this->m_I_nsamples; i += 16ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);
                                                                     
                    }
              }
                            
             for(; (i+7ull) < this->m_I_nsamples; i += 8ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                                                                                            
                    }
              }

              for(; (i+0ull) < this->m_I_nsamples; i += 1ull) 
              {
                   sum_rem_i = 0.0f;
                   sum_rem_q = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_i[j] = this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i];
                       sum_rem_i     += scal_rem_i[j];
                       this->m_I_vsequence.m_data[i] = sum_rem_i;
                       scal_rem_q[j] = this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i];
                       sum_rem_q     += scal_rem_q[j];
                       this->m_Q_vsequence.m_data[i] = sum_rem_q;
                   }
              }

        }
        else 
        {
              for(i = 0ull; (i+47ull) < this->m_I_nsamples; i += 48ull) 
              {
                    sum0_i = vzero;                
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                   
                    sum4_i = vzero;                   
                    sum5_i = vzero;                                                     
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);                        
                                                                
#endif 
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);                       
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);                       
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);                        
                        xmm4_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum4_i    = _mm256_add_ps(xmm4_i[j],sum4_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum4_i);                        
                        xmm5_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+40ull]);
                        sum5_i    = _mm256_add_ps(xmm5_i[j],sum5_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+40ull],sum5_i);                                                    
                    }
              }

             
              for(; (i+31ull) < this->m_I_nsamples; i += 32ull) 
              {
                    sum0_i = vzero;                
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                                                                
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);                       
                        xmm2_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum2_i    = _mm256_add_ps(xmm2_i[j],sum2_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum2_i);                       
                        xmm3_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+24ull]);
                        sum3_i    = _mm256_add_ps(xmm3_i[j],sum3_i);                       
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+24ull],sum3_i);                        
                                                                            
                    }
              }

              for(; (i+15ull) < this->m_I_nsamples; i += 16ull) 
              {
                    sum0_i = vzero;                
                    sum1_i = vzero;                                                                    
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+8ull]);
                        sum1_i    = _mm256_add_ps(xmm1_i[j],sum1_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+8ull],sum1_i);                       
                                                                           
                    }
              } 
                            
             for(; (i+7ull) < this->m_I_nsamples; i += 8ull) 
             {
                    sum0_i = vzero;                                                                                                     
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm256_add_ps(xmm0_i[j],sum0_i);                        
                        _mm256_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                                                                                               
                    }
              } 

              for(; (i+0ull) < this->m_I_nsamples; i += 1ull) 
              {
                   sum_rem_i = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_i[j] = this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i];
                       sum_rem_i     += scal_rem_i[j];
                       this->m_I_vsequence.m_data[i] = sum_rem_i;
                   }
              }

              ///////////////////////////// Q-channel //////////////////////////////////////////////
              for(i = 0ull; (i+47ull) < this->m_Q_nsamples; i += 48ull) 
              {
                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                                    
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);                                            
#endif 
                        
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                      
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);                       
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);                        
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);                       
                        xmm4_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum4_q    = _mm256_add_ps(xmm4_q[j],sum4_q);                     
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum4_q);                       
                        xmm5_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+40ull]);
                        sum5_q    = _mm256_add_ps(xmm5_q[j],sum5_q);                      
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+40ull],sum5_q);                       
                                            
                    }
              }
              
            
              for(; (i+31ull) < this->m_Q_nsamples; i += 32ull) 
              {
                    
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                                     
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {                        
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                      
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);                       
                        xmm2_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum2_q    = _mm256_add_ps(xmm2_q[j],sum2_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum2_q);                        
                        xmm3_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+24ull]);
                        sum3_q    = _mm256_add_ps(xmm3_q[j],sum3_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+24ull],sum3_q);                       
                    }
              }

              for(; (i+15ull) < this->m_Q_nsamples; i += 16ull) 
              {
                    sum0_q = vzero;                    
                    sum1_q = vzero;                                    
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                      
                        xmm1_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+8ull]);
                        sum1_q    = _mm256_add_ps(xmm1_q[j],sum1_q);                       
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+8ull],sum1_q);                       
                    }
              }

              for(; (i+7ull) < this->m_Q_nsamples; i += 8ull) 
              {                    
                    sum0_q = vzero;                                                     
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        
                        xmm0_q[j] = _mm256_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm256_add_ps(xmm0_q[j],sum0_q);                        
                        _mm256_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                                            
                    }
              }
               
               for(; (i+0ull) < this->m_Q_nsamples; i += 1ull) 
               {                 
                   sum_rem_q = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_q[j] = this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i];
                       sum_rem_q     += scal_rem_q[j];
                       this->m_Q_vsequence.m_data[i] = sum_rem_q;
                   }
              }
        }
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
           HW_PMC_COLLECTION_EPILOGE_BODY
           HW_PMC_SHOW_INFO
#endif 
     return (0);
}

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SCALAR_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);


template<gms::radiolocation::iq_rectw_bitstream_vsequence_t
         ::IQ_rectw_bitstreams_optim_path  optim_path>
std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u16x(std::int32_t * ret_codes)
{
       if(__builtin_expect(nullptr==ret_codes,0)) { return (-1);}
       if(__builtin_expect(this->m_IQ_nsignals>256ull,0)) { return (-2);}
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
           PMC_VARS
           HW_PMC_COLLECTION_PROLOGE_BODY
#endif 
 
        constexpr std::size_t nsignals_size{256ull};
        __ATTR_ALIGN__(64) 
        __m512 xmm0_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm1_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm2_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm3_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm4_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm5_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm6_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm7_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm8_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm9_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm10_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm11_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm12_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm13_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm14_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm15_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm0_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm1_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm2_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm3_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm4_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm5_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm6_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm7_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm8_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm9_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm10_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm11_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm12_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm13_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm14_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm15_q[nsignals_size];
        float scal_rem_i[nsignals_size];
        float scal_rem_q[nsignals_size];
        __m512 vzero{_mm512_setzero_ps()};
        __m512 sum0_i{vzero};
        __m512 sum1_i{sum0_i};
        __m512 sum2_i{sum0_i};
        __m512 sum3_i{sum0_i};
        __m512 sum4_i{sum0_i};
        __m512 sum5_i{sum0_i};
        __m512 sum6_i{sum0_i};
        __m512 sum7_i{sum0_i};
        __m512 sum8_i{sum0_i};
        __m512 sum9_i{sum0_i};
        __m512 sum10_i{sum0_i};
        __m512 sum11_i{sum0_i};
        __m512 sum12_i{sum0_i};
        __m512 sum13_i{sum0_i};
        __m512 sum14_i{sum0_i};
        __m512 sum15_i{sum0_i};
        __m512 sum0_q{sum0_i};
        __m512 sum1_q{sum0_i};
        __m512 sum2_q{sum0_i};
        __m512 sum3_q{sum0_i};
        __m512 sum4_q{sum0_i};
        __m512 sum5_q{sum0_i};
        __m512 sum6_q{sum0_i};
        __m512 sum7_q{sum0_i};
        __m512 sum8_q{sum0_i};
        __m512 sum9_q{sum0_i};
        __m512 sum10_q{sum0_i};
        __m512 sum11_q{sum0_i};
        __m512 sum12_q{sum0_i};
        __m512 sum13_q{sum0_i};
        __m512 sum14_q{sum0_i};
        __m512 sum15_q{sum0_i}; 
        std::size_t i;
        float sum_rem_i{0.0f};
        float sum_rem_q{0.0f};

        if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SCALAR_PATH)
        {          
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_scalar();
                }            
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_sse();
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_sse_u4x();
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx_u4x(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx_u4x(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx512_u4x(true);
                } 
        }

        if(this->m_I_nsamples==this->m_Q_nsamples)
        {
              for(i = 0ull; (i+255ull) < this->m_I_nsamples; i += 256ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    sum6_i = vzero;
                    sum6_q = vzero;
                    sum7_i = vzero;
                    sum7_q = vzero;
                    sum8_i = vzero;
                    sum8_q = vzero;
                    sum9_i = vzero;
                    sum9_q = vzero;
                    sum10_i = vzero;
                    sum10_q = vzero;
                    sum11_i = vzero;
                    sum11_q = vzero;
                    sum12_i = vzero;
                    sum12_q = vzero;
                    sum13_i = vzero;
                    sum13_q = vzero;
                    sum14_i = vzero;
                    sum14_q = vzero;
                    sum15_i = vzero;
                    sum15_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+128ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+128ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+144ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+144ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+160ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+160ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+176ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+176ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+192ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+192ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+208ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+208ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+224ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+224ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+240ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+240ull],_MM_HINT_T1);
#endif 
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);
                        xmm4_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum4_i    = _mm512_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum4_q    = _mm512_add_ps(xmm4_q[j],sum4_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum4_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum4_q);
                        xmm5_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull]);
                        sum5_i    = _mm512_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull]);
                        sum5_q    = _mm512_add_ps(xmm5_q[j],sum5_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+80ull],sum5_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+80ull],sum5_q);
                        xmm6_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull]);
                        sum6_i    = _mm512_add_ps(xmm6_i[j],sum6_i);
                        xmm6_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull]);
                        sum6_q    = _mm512_add_ps(xmm6_q[j],sum6_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+96ull],sum6_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+96ull],sum6_q);
                        xmm7_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull]);
                        sum7_i    = _mm512_add_ps(xmm7_i[j],sum7_i);
                        xmm7_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull]);
                        sum7_q    = _mm512_add_ps(xmm7_q[j],sum7_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+112ull],sum7_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+112ull],sum7_q);
                        xmm8_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+128ull]);
                        sum8_i    = _mm512_add_ps(xmm8_i[j],sum8_i);
                        xmm8_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+128ull]);
                        sum8_q    = _mm512_add_ps(xmm8_q[j],sum8_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+128ull],sum8_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+128ull],sum8_q);
                        xmm9_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+144ull]);
                        sum9_i    = _mm512_add_ps(xmm9_i[j],sum9_i);
                        xmm9_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+144ull]);
                        sum9_q    = _mm512_add_ps(xmm9_q[j],sum9_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+144ull],sum9_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+144ull],sum9_q);
                        xmm10_i[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+160ull]);
                        sum10_i   = _mm512_add_ps(xmm10_i[j],sum10_i);
                        xmm10_q[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+160ull]);
                        sum10_q   = _mm512_add_ps(xmm10_q[j],sum10_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+160ull],sum10_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+160ull],sum10_q);
                        xmm11_i[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+176ull]);
                        sum11_i   = _mm512_add_ps(xmm11_i[j],sum11_i);
                        xmm11_q[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+176ull]);
                        sum11_q   = _mm512_add_ps(xmm11_q[j],sum11_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+176ull],sum11_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+176ull],sum11_q);
                        xmm12_i[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+192ull]);
                        sum12_i   = _mm512_add_ps(xmm12_i[j],sum12_i);
                        xmm12_q[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+192ull]);
                        sum12_q   = _mm512_add_ps(xmm12_q[j],sum12_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+192ull],sum12_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+192ull],sum12_q);
                        xmm13_i[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+208ull]);
                        sum13_i   = _mm512_add_ps(xmm13_i[j],sum13_i);
                        xmm13_q[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+208ull]);
                        sum13_q   = _mm512_add_ps(xmm13_q[j],sum13_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+208ull],sum13_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+208ull],sum13_q);
                        xmm14_i[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+224ull]);
                        sum14_i   = _mm512_add_ps(xmm14_i[j],sum14_i);
                        xmm14_q[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+224ull]);
                        sum14_q   = _mm512_add_ps(xmm14_q[j],sum14_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+224ull],sum14_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+224ull],sum14_q);
                        xmm15_i[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+240ull]);
                        sum15_i   = _mm512_add_ps(xmm15_i[j],sum15_i);
                        xmm15_q[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+240ull]);
                        sum15_q   = _mm512_add_ps(xmm15_q[j],sum15_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+240ull],sum15_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+240ull],sum15_q);
                       
                    }
              }
              
             for(; (i+191ull) < this->m_I_nsamples; i += 192ull) 
             {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    sum6_i = vzero;
                    sum6_q = vzero;
                    sum7_i = vzero;
                    sum7_q = vzero;
                    sum8_i = vzero;
                    sum8_q = vzero;
                    sum9_i = vzero;
                    sum9_q = vzero;
                    sum10_i = vzero;
                    sum10_q = vzero;
                    sum11_i = vzero;
                    sum11_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);
                        xmm4_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum4_i    = _mm512_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum4_q    = _mm512_add_ps(xmm4_q[j],sum4_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum4_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum4_q);
                        xmm5_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull]);
                        sum5_i    = _mm512_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull]);
                        sum5_q    = _mm512_add_ps(xmm5_q[j],sum5_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+80ull],sum5_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+80ull],sum5_q);
                        xmm6_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull]);
                        sum6_i    = _mm512_add_ps(xmm6_i[j],sum6_i);
                        xmm6_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull]);
                        sum6_q    = _mm512_add_ps(xmm6_q[j],sum6_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+96ull],sum6_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+96ull],sum6_q);
                        xmm7_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull]);
                        sum7_i    = _mm512_add_ps(xmm7_i[j],sum7_i);
                        xmm7_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull]);
                        sum7_q    = _mm512_add_ps(xmm7_q[j],sum7_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+112ull],sum7_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+112ull],sum7_q);
                        xmm8_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+128ull]);
                        sum8_i    = _mm512_add_ps(xmm8_i[j],sum8_i);
                        xmm8_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+128ull]);
                        sum8_q    = _mm512_add_ps(xmm8_q[j],sum8_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+128ull],sum8_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+128ull],sum8_q);
                        xmm9_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+144ull]);
                        sum9_i    = _mm512_add_ps(xmm9_i[j],sum9_i);
                        xmm9_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+144ull]);
                        sum9_q    = _mm512_add_ps(xmm9_q[j],sum9_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+144ull],sum9_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+144ull],sum9_q);
                        xmm10_i[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+160ull]);
                        sum10_i   = _mm512_add_ps(xmm10_i[j],sum10_i);
                        xmm10_q[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+160ull]);
                        sum10_q   = _mm512_add_ps(xmm10_q[j],sum10_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+160ull],sum10_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+160ull],sum10_q);
                        xmm11_i[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+176ull]);
                        sum11_i   = _mm512_add_ps(xmm11_i[j],sum11_i);
                        xmm11_q[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+176ull]);
                        sum11_q   = _mm512_add_ps(xmm11_q[j],sum11_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+176ull],sum11_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+176ull],sum11_q);
                       
                    }
              }
             
             for(; (i+159ull) < this->m_I_nsamples; i += 160ull) 
             {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    sum6_i = vzero;
                    sum6_q = vzero;
                    sum7_i = vzero;
                    sum7_q = vzero;
                    sum8_i = vzero;
                    sum8_q = vzero;
                    sum9_i = vzero;
                    sum9_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);
                        xmm4_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum4_i    = _mm512_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum4_q    = _mm512_add_ps(xmm4_q[j],sum4_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum4_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum4_q);
                        xmm5_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull]);
                        sum5_i    = _mm512_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull]);
                        sum5_q    = _mm512_add_ps(xmm5_q[j],sum5_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+80ull],sum5_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+80ull],sum5_q);
                        xmm6_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull]);
                        sum6_i    = _mm512_add_ps(xmm6_i[j],sum6_i);
                        xmm6_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull]);
                        sum6_q    = _mm512_add_ps(xmm6_q[j],sum6_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+96ull],sum6_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+96ull],sum6_q);
                        xmm7_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull]);
                        sum7_i    = _mm512_add_ps(xmm7_i[j],sum7_i);
                        xmm7_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull]);
                        sum7_q    = _mm512_add_ps(xmm7_q[j],sum7_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+112ull],sum7_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+112ull],sum7_q);
                        xmm8_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+128ull]);
                        sum8_i    = _mm512_add_ps(xmm8_i[j],sum8_i);
                        xmm8_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+128ull]);
                        sum8_q    = _mm512_add_ps(xmm8_q[j],sum8_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+128ull],sum8_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+128ull],sum8_q);
                        xmm9_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+144ull]);
                        sum9_i    = _mm512_add_ps(xmm9_i[j],sum9_i);
                        xmm9_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+144ull]);
                        sum9_q    = _mm512_add_ps(xmm9_q[j],sum9_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+144ull],sum9_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+144ull],sum9_q);
                        
                    }
              }
              
              for(; (i+127ull) < this->m_I_nsamples; i += 128ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    sum6_i = vzero;
                    sum6_q = vzero;
                    sum7_i = vzero;
                    sum7_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);
                        xmm4_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum4_i    = _mm512_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum4_q    = _mm512_add_ps(xmm4_q[j],sum4_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum4_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum4_q);
                        xmm5_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull]);
                        sum5_i    = _mm512_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull]);
                        sum5_q    = _mm512_add_ps(xmm5_q[j],sum5_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+80ull],sum5_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+80ull],sum5_q);
                        xmm6_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull]);
                        sum6_i    = _mm512_add_ps(xmm6_i[j],sum6_i);
                        xmm6_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull]);
                        sum6_q    = _mm512_add_ps(xmm6_q[j],sum6_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+96ull],sum6_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+96ull],sum6_q);
                        xmm7_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull]);
                        sum7_i    = _mm512_add_ps(xmm7_i[j],sum7_i);
                        xmm7_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull]);
                        sum7_q    = _mm512_add_ps(xmm7_q[j],sum7_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+112ull],sum7_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+112ull],sum7_q);
                        
                    }
              }
               
              for(; (i+63ull) < this->m_I_nsamples; i += 64ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);
                        
                    }
              }

              for(; (i+31ull) < this->m_I_nsamples; i += 32ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);
                        
                    }
              }

               for(; (i+15ull) < this->m_I_nsamples; i += 16ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                     
                    }
              }

                         
              for(; (i+0ull) < this->m_I_nsamples; i += 1ull) 
              {
                   sum_rem_i = 0.0f;
                   sum_rem_q = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_i[j] = this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i];
                       sum_rem_i     += scal_rem_i[j];
                       this->m_I_vsequence.m_data[i] = sum_rem_i;
                       scal_rem_q[j] = this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i];
                       sum_rem_q     += scal_rem_q[j];
                       this->m_Q_vsequence.m_data[i] = sum_rem_q;
                   }
              }

        }
        else 
        {
             for(i = 0ull; (i+255ull) < this->m_I_nsamples; i += 256ull) 
             {
                    sum0_i = vzero;                    
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                    
                    sum4_i = vzero;                    
                    sum5_i = vzero;                   
                    sum6_i = vzero;                    
                    sum7_i = vzero;                    
                    sum8_i = vzero;
                    sum9_i = vzero;                   
                    sum10_i = vzero;                    
                    sum11_i = vzero;                    
                    sum12_i = vzero;                   
                    sum13_i = vzero;                    
                    sum14_i = vzero;                    
                    sum15_i = vzero;                    
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+128ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+144ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+160ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+176ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+192ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+208ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+224ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+240ull],_MM_HINT_T1);                      
#endif 
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);                        
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);                      
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);                       
                        xmm4_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum4_i    = _mm512_add_ps(xmm4_i[j],sum4_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum4_i);                       
                        xmm5_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull]);
                        sum5_i    = _mm512_add_ps(xmm5_i[j],sum5_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+80ull],sum5_i);                        
                        xmm6_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull]);
                        sum6_i    = _mm512_add_ps(xmm6_i[j],sum6_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+96ull],sum6_i);                        
                        xmm7_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull]);
                        sum7_i    = _mm512_add_ps(xmm7_i[j],sum7_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+112ull],sum7_i);                      
                        xmm8_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+128ull]);
                        sum8_i    = _mm512_add_ps(xmm8_i[j],sum8_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+128ull],sum8_i);                        
                        xmm9_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+144ull]);
                        sum9_i    = _mm512_add_ps(xmm9_i[j],sum9_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+144ull],sum9_i);                       
                        xmm10_i[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+160ull]);
                        sum10_i   = _mm512_add_ps(xmm10_i[j],sum10_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+160ull],sum10_i);                        
                        xmm11_i[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+176ull]);
                        sum11_i   = _mm512_add_ps(xmm11_i[j],sum11_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+176ull],sum11_i);                        
                        xmm12_i[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+192ull]);
                        sum12_i   = _mm512_add_ps(xmm12_i[j],sum12_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+192ull],sum12_i);                        
                        xmm13_i[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+208ull]);
                        sum13_i   = _mm512_add_ps(xmm13_i[j],sum13_i);                      
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+208ull],sum13_i);                       
                        xmm14_i[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+224ull]);
                        sum14_i   = _mm512_add_ps(xmm14_i[j],sum14_i);                      
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+224ull],sum14_i);                       
                        xmm15_i[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+240ull]);
                        sum15_i   = _mm512_add_ps(xmm15_i[j],sum15_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+240ull],sum15_i);
                                              
                    }
              }
              
              for(; (i+191ull) < this->m_I_nsamples; i += 192ull) 
              {
                    sum0_i = vzero;                    
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                    
                    sum4_i = vzero;                    
                    sum5_i = vzero;                   
                    sum6_i = vzero;                    
                    sum7_i = vzero;                    
                    sum8_i = vzero;
                    sum9_i = vzero;                   
                    sum10_i = vzero;                    
                    sum11_i = vzero;                                     
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);                        
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);                      
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);                       
                        xmm4_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum4_i    = _mm512_add_ps(xmm4_i[j],sum4_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum4_i);                       
                        xmm5_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull]);
                        sum5_i    = _mm512_add_ps(xmm5_i[j],sum5_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+80ull],sum5_i);                        
                        xmm6_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull]);
                        sum6_i    = _mm512_add_ps(xmm6_i[j],sum6_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+96ull],sum6_i);                        
                        xmm7_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull]);
                        sum7_i    = _mm512_add_ps(xmm7_i[j],sum7_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+112ull],sum7_i);                      
                        xmm8_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+128ull]);
                        sum8_i    = _mm512_add_ps(xmm8_i[j],sum8_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+128ull],sum8_i);                        
                        xmm9_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+144ull]);
                        sum9_i    = _mm512_add_ps(xmm9_i[j],sum9_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+144ull],sum9_i);                       
                        xmm10_i[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+160ull]);
                        sum10_i   = _mm512_add_ps(xmm10_i[j],sum10_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+160ull],sum10_i);                        
                        xmm11_i[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+176ull]);
                        sum11_i   = _mm512_add_ps(xmm11_i[j],sum11_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+176ull],sum11_i);                        
                                              
                    }
              }

               for(; (i+159ull) < this->m_I_nsamples; i += 160ull) 
               {
                    sum0_i = vzero;                    
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                    
                    sum4_i = vzero;                    
                    sum5_i = vzero;                   
                    sum6_i = vzero;                    
                    sum7_i = vzero;                    
                    sum8_i = vzero;
                    sum9_i = vzero;                                                   
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);                        
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);                      
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);                       
                        xmm4_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum4_i    = _mm512_add_ps(xmm4_i[j],sum4_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum4_i);                       
                        xmm5_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull]);
                        sum5_i    = _mm512_add_ps(xmm5_i[j],sum5_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+80ull],sum5_i);                        
                        xmm6_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull]);
                        sum6_i    = _mm512_add_ps(xmm6_i[j],sum6_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+96ull],sum6_i);                        
                        xmm7_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull]);
                        sum7_i    = _mm512_add_ps(xmm7_i[j],sum7_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+112ull],sum7_i);                      
                        xmm8_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+128ull]);
                        sum8_i    = _mm512_add_ps(xmm8_i[j],sum8_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+128ull],sum8_i);                        
                        xmm9_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+144ull]);
                        sum9_i    = _mm512_add_ps(xmm9_i[j],sum9_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+144ull],sum9_i);                       
                                                                 
                    }
              }

               for(; (i+127ull) < this->m_I_nsamples; i += 128ull) 
               {
                    sum0_i = vzero;                    
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                    
                    sum4_i = vzero;                    
                    sum5_i = vzero;                   
                    sum6_i = vzero;                    
                    sum7_i = vzero;                                                                    
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);                        
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);                      
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);                       
                        xmm4_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum4_i    = _mm512_add_ps(xmm4_i[j],sum4_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum4_i);                       
                        xmm5_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull]);
                        sum5_i    = _mm512_add_ps(xmm5_i[j],sum5_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+80ull],sum5_i);                        
                        xmm6_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull]);
                        sum6_i    = _mm512_add_ps(xmm6_i[j],sum6_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+96ull],sum6_i);                        
                        xmm7_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull]);
                        sum7_i    = _mm512_add_ps(xmm7_i[j],sum7_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+112ull],sum7_i);                      
                                                                               
                    }
              }

              for(; (i+63ull) < this->m_I_nsamples; i += 64ull) 
              {
                    sum0_i = vzero;                    
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                                                                                     
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);                        
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);                      
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);                       
                                                                                               
                    }
              }

              for(; (i+31ull) < this->m_I_nsamples; i += 32ull) 
              {
                    sum0_i = vzero;                    
                    sum1_i = vzero;                                                                                                       
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);                        
                                                                                                                
                    }
              }

              for(; (i+15ull) < this->m_I_nsamples; i += 16ull) 
              {
                    sum0_i = vzero;                                                                                                                        
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                                                                                                                                  
                    }
              
              }

              for(; (i+0ull) < this->m_I_nsamples; i += 1ull) 
              {
                   sum_rem_i = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_i[j] = this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i];
                       sum_rem_i     += scal_rem_i[j];
                       this->m_I_vsequence.m_data[i] = sum_rem_i;
                   }
              }

              ///////////////////////////// Q-channel //////////////////////////////////////////////
              
              for(i = 0ull; (i+255ull) < this->m_Q_nsamples; i += 256ull) 
              {
                   
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                    
                    sum6_q = vzero;                    
                    sum7_q = vzero;                    
                    sum8_q = vzero;                   
                    sum9_q = vzero;                    
                    sum10_q = vzero;                    
                    sum11_q = vzero;                    
                    sum12_q = vzero;                    
                    sum13_q = vzero;                   
                    sum14_q = vzero;                    
                    sum15_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+128ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+144ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+160ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+176ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+192ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+208ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+224ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+240ull],_MM_HINT_T1);
#endif 
                       
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);                       
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);                        
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);                        
                        xmm4_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum4_q    = _mm512_add_ps(xmm4_q[j],sum4_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum4_q);                        
                        xmm5_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull]);
                        sum5_q    = _mm512_add_ps(xmm5_q[j],sum5_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+80ull],sum5_q);                        
                        xmm6_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull]);
                        sum6_q    = _mm512_add_ps(xmm6_q[j],sum6_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+96ull],sum6_q);                        
                        xmm7_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull]);
                        sum7_q    = _mm512_add_ps(xmm7_q[j],sum7_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+112ull],sum7_q);                        
                        xmm8_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+128ull]);
                        sum8_q    = _mm512_add_ps(xmm8_q[j],sum8_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+128ull],sum8_q);                        
                        xmm9_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+144ull]);
                        sum9_q    = _mm512_add_ps(xmm9_q[j],sum9_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+144ull],sum9_q);                        
                        xmm10_q[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+160ull]);
                        sum10_q   = _mm512_add_ps(xmm10_q[j],sum10_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+160ull],sum10_q);                       
                        xmm11_q[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+176ull]);
                        sum11_q   = _mm512_add_ps(xmm11_q[j],sum11_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+176ull],sum11_q);                       
                        xmm12_q[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+192ull]);
                        sum12_q   = _mm512_add_ps(xmm12_q[j],sum12_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+192ull],sum12_q);                       
                        xmm13_q[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+208ull]);
                        sum13_q   = _mm512_add_ps(xmm13_q[j],sum13_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+208ull],sum13_q);                      
                        xmm14_q[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+224ull]);
                        sum14_q   = _mm512_add_ps(xmm14_q[j],sum14_q);                     
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+224ull],sum14_q);                       
                        xmm15_q[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+240ull]);
                        sum15_q   = _mm512_add_ps(xmm15_q[j],sum15_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+240ull],sum15_q);
                       
                    }
              }

              for(; (i+191ull) < this->m_Q_nsamples; i += 192ull) 
              {
                   
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                    
                    sum6_q = vzero;                    
                    sum7_q = vzero;                    
                    sum8_q = vzero;                   
                    sum9_q = vzero;                    
                    sum10_q = vzero;                    
                    sum11_q = vzero;                    
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                       
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);                       
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);                        
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);                        
                        xmm4_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum4_q    = _mm512_add_ps(xmm4_q[j],sum4_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum4_q);                        
                        xmm5_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull]);
                        sum5_q    = _mm512_add_ps(xmm5_q[j],sum5_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+80ull],sum5_q);                        
                        xmm6_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull]);
                        sum6_q    = _mm512_add_ps(xmm6_q[j],sum6_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+96ull],sum6_q);                        
                        xmm7_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull]);
                        sum7_q    = _mm512_add_ps(xmm7_q[j],sum7_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+112ull],sum7_q);                        
                        xmm8_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+128ull]);
                        sum8_q    = _mm512_add_ps(xmm8_q[j],sum8_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+128ull],sum8_q);                        
                        xmm9_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+144ull]);
                        sum9_q    = _mm512_add_ps(xmm9_q[j],sum9_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+144ull],sum9_q);                        
                        xmm10_q[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+160ull]);
                        sum10_q   = _mm512_add_ps(xmm10_q[j],sum10_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+160ull],sum10_q);                       
                        xmm11_q[j]= _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+176ull]);
                        sum11_q   = _mm512_add_ps(xmm11_q[j],sum11_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+176ull],sum11_q);                       
                       
                    }
              }

              for(; (i+159ull) < this->m_Q_nsamples; i += 160ull) 
              {
                   
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                    
                    sum6_q = vzero;                    
                    sum7_q = vzero;                    
                    sum8_q = vzero;                   
                    sum9_q = vzero;                                 
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                       
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);                       
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);                        
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);                        
                        xmm4_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum4_q    = _mm512_add_ps(xmm4_q[j],sum4_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum4_q);                        
                        xmm5_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull]);
                        sum5_q    = _mm512_add_ps(xmm5_q[j],sum5_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+80ull],sum5_q);                        
                        xmm6_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull]);
                        sum6_q    = _mm512_add_ps(xmm6_q[j],sum6_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+96ull],sum6_q);                        
                        xmm7_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull]);
                        sum7_q    = _mm512_add_ps(xmm7_q[j],sum7_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+112ull],sum7_q);                        
                        xmm8_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+128ull]);
                        sum8_q    = _mm512_add_ps(xmm8_q[j],sum8_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+128ull],sum8_q);                        
                        xmm9_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+144ull]);
                        sum9_q    = _mm512_add_ps(xmm9_q[j],sum9_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+144ull],sum9_q);                        
                                          
                    }
              }

              for(; (i+127ull) < this->m_Q_nsamples; i += 128ull) 
              {
                   
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                    
                    sum6_q = vzero;                    
                    sum7_q = vzero;                                                   
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                       
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);                       
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);                        
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);                        
                        xmm4_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum4_q    = _mm512_add_ps(xmm4_q[j],sum4_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum4_q);                        
                        xmm5_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull]);
                        sum5_q    = _mm512_add_ps(xmm5_q[j],sum5_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+80ull],sum5_q);                        
                        xmm6_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull]);
                        sum6_q    = _mm512_add_ps(xmm6_q[j],sum6_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+96ull],sum6_q);                        
                        xmm7_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull]);
                        sum7_q    = _mm512_add_ps(xmm7_q[j],sum7_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+112ull],sum7_q);                        
                                                                 
                    }
              }

              for(; (i+63ull) < this->m_Q_nsamples; i += 64ull) 
              {
                   
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                                                                     
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                       
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);                       
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);                        
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);                        
                                                                                      
                    }
              }

              for(; (i+31ull) < this->m_Q_nsamples; i += 32ull) 
              {
                   
                    sum0_q = vzero;                    
                    sum1_q = vzero;                                                                                       
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                       
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);                       
                                                                                                            
                    }
              }

              for(; (i+15ull) < this->m_Q_nsamples; i += 16ull) 
              {                  
                    sum0_q = vzero;                                                                                                                           
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                       
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                                                                                                                                  
                    }
              }
               
               for(; (i+0ull) < this->m_Q_nsamples; i += 1ull) 
               {                 
                   sum_rem_q = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_q[j] = this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i];
                       sum_rem_q     += scal_rem_q[j];
                       this->m_Q_vsequence.m_data[i] = sum_rem_q;
                   }
              }
        }
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
           HW_PMC_COLLECTION_EPILOGE_BODY
           HW_PMC_SHOW_INFO
#endif 
     return (0);
}

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SCALAR_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u16x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);

template<gms::radiolocation::iq_rectw_bitstream_vsequence_t
         ::IQ_rectw_bitstreams_optim_path  optim_path>         
std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u10x(std::int32_t * ret_codes)
{
       if(__builtin_expect(nullptr==ret_codes,0)) { return (-1);}
       if(__builtin_expect(this->m_IQ_nsignals>256ull,0)) { return (-2);}
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
           PMC_VARS
           HW_PMC_COLLECTION_PROLOGE_BODY
#endif 

        constexpr std::size_t nsignals_size{256ull};
        __ATTR_ALIGN__(64) 
        __m512 xmm0_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm1_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm2_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm3_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm4_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm5_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm6_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm7_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm8_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm9_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm0_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm1_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm2_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm3_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm4_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm5_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm6_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm7_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm8_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm9_q[nsignals_size];
        float scal_rem_i[nsignals_size];
        float scal_rem_q[nsignals_size];
        __m512 vzero{_mm512_setzero_ps()};
        __m512 sum0_i{vzero};
        __m512 sum1_i{sum0_i};
        __m512 sum2_i{sum0_i};
        __m512 sum3_i{sum0_i};
        __m512 sum4_i{sum0_i};
        __m512 sum5_i{sum0_i};
        __m512 sum6_i{sum0_i};
        __m512 sum7_i{sum0_i};
        __m512 sum8_i{sum0_i};
        __m512 sum9_i{sum0_i};
        __m512 sum0_q{sum0_i};
        __m512 sum1_q{sum0_i};
        __m512 sum2_q{sum0_i};
        __m512 sum3_q{sum0_i};
        __m512 sum4_q{sum0_i};
        __m512 sum5_q{sum0_i};
        __m512 sum6_q{sum0_i};
        __m512 sum7_q{sum0_i};
        __m512 sum8_q{sum0_i};
        __m512 sum9_q{sum0_i}; 
        std::size_t i;
        float sum_rem_i{0.0f};
        float sum_rem_q{0.0f};

        if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SCALAR_PATH)
        {          
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_scalar();
                }            
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_sse();
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_sse_u4x();
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx_u4x(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx_u4x(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx512_u4x(true);
                } 
        }

        if(this->m_I_nsamples==this->m_Q_nsamples)
        {
              for(i = 0ull; (i+159ull) < this->m_I_nsamples; i += 160ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    sum6_i = vzero;
                    sum6_q = vzero;
                    sum7_i = vzero;
                    sum7_q = vzero;
                    sum8_i = vzero;
                    sum8_q = vzero;
                    sum9_i = vzero;
                    sum9_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+128ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+128ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+144ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+144ull],_MM_HINT_T1);
#endif 
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);
                        xmm4_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum4_i    = _mm512_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum4_q    = _mm512_add_ps(xmm4_q[j],sum4_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum4_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum4_q);
                        xmm5_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull]);
                        sum5_i    = _mm512_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull]);
                        sum5_q    = _mm512_add_ps(xmm5_q[j],sum5_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+80ull],sum5_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+80ull],sum5_q);
                        xmm6_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull]);
                        sum6_i    = _mm512_add_ps(xmm6_i[j],sum6_i);
                        xmm6_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull]);
                        sum6_q    = _mm512_add_ps(xmm6_q[j],sum6_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+96ull],sum6_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+96ull],sum6_q);
                        xmm7_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull]);
                        sum7_i    = _mm512_add_ps(xmm7_i[j],sum7_i);
                        xmm7_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull]);
                        sum7_q    = _mm512_add_ps(xmm7_q[j],sum7_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+112ull],sum7_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+112ull],sum7_q);
                        xmm8_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+128ull]);
                        sum8_i    = _mm512_add_ps(xmm8_i[j],sum8_i);
                        xmm8_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+128ull]);
                        sum8_q    = _mm512_add_ps(xmm8_q[j],sum8_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+128ull],sum8_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+128ull],sum8_q);
                        xmm9_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+144ull]);
                        sum9_i    = _mm512_add_ps(xmm9_i[j],sum9_i);
                        xmm9_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+144ull]);
                        sum9_q    = _mm512_add_ps(xmm9_q[j],sum9_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+144ull],sum9_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+144ull],sum9_q);
                       
                    }
              }
              
              
              for(; (i+127ull) < this->m_I_nsamples; i += 128ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    sum6_i = vzero;
                    sum6_q = vzero;
                    sum7_i = vzero;
                    sum7_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);
                        xmm4_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum4_i    = _mm512_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum4_q    = _mm512_add_ps(xmm4_q[j],sum4_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum4_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum4_q);
                        xmm5_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull]);
                        sum5_i    = _mm512_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull]);
                        sum5_q    = _mm512_add_ps(xmm5_q[j],sum5_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+80ull],sum5_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+80ull],sum5_q);
                        xmm6_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull]);
                        sum6_i    = _mm512_add_ps(xmm6_i[j],sum6_i);
                        xmm6_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull]);
                        sum6_q    = _mm512_add_ps(xmm6_q[j],sum6_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+96ull],sum6_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+96ull],sum6_q);
                        xmm7_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull]);
                        sum7_i    = _mm512_add_ps(xmm7_i[j],sum7_i);
                        xmm7_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull]);
                        sum7_q    = _mm512_add_ps(xmm7_q[j],sum7_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+112ull],sum7_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+112ull],sum7_q);
                        
                    }
              }
               
              for(; (i+63ull) < this->m_I_nsamples; i += 64ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);
                        
                    }
              }

              for(; (i+31ull) < this->m_I_nsamples; i += 32ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);
                        
                    }
              }

               for(; (i+15ull) < this->m_I_nsamples; i += 16ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                     
                    }
              }

                         
              for(; (i+0ull) < this->m_I_nsamples; i += 1ull) 
              {
                   sum_rem_i = 0.0f;
                   sum_rem_q = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_i[j] = this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i];
                       sum_rem_i     += scal_rem_i[j];
                       this->m_I_vsequence.m_data[i] = sum_rem_i;
                       scal_rem_q[j] = this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i];
                       sum_rem_q     += scal_rem_q[j];
                       this->m_Q_vsequence.m_data[i] = sum_rem_q;
                   }
              }

        }
        else 
        {
             for(i = 0ull; (i+159ull) < this->m_I_nsamples; i += 160ull) 
             {
                    sum0_i = vzero;                    
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                    
                    sum4_i = vzero;                    
                    sum5_i = vzero;                   
                    sum6_i = vzero;                    
                    sum7_i = vzero;                    
                    sum8_i = vzero;
                    sum9_i = vzero;                                      
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+128ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+144ull],_MM_HINT_T1);                                          
#endif 
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);                        
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);                      
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);                       
                        xmm4_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum4_i    = _mm512_add_ps(xmm4_i[j],sum4_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum4_i);                       
                        xmm5_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull]);
                        sum5_i    = _mm512_add_ps(xmm5_i[j],sum5_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+80ull],sum5_i);                        
                        xmm6_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull]);
                        sum6_i    = _mm512_add_ps(xmm6_i[j],sum6_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+96ull],sum6_i);                        
                        xmm7_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull]);
                        sum7_i    = _mm512_add_ps(xmm7_i[j],sum7_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+112ull],sum7_i);                      
                        xmm8_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+128ull]);
                        sum8_i    = _mm512_add_ps(xmm8_i[j],sum8_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+128ull],sum8_i);                        
                        xmm9_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+144ull]);
                        sum9_i    = _mm512_add_ps(xmm9_i[j],sum9_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+144ull],sum9_i);                       
                                              
                    }
              }
              
               for(; (i+127ull) < this->m_I_nsamples; i += 128ull) 
               {
                    sum0_i = vzero;                    
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                    
                    sum4_i = vzero;                    
                    sum5_i = vzero;                   
                    sum6_i = vzero;                    
                    sum7_i = vzero;                                                                    
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);                        
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);                      
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);                       
                        xmm4_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum4_i    = _mm512_add_ps(xmm4_i[j],sum4_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum4_i);                       
                        xmm5_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull]);
                        sum5_i    = _mm512_add_ps(xmm5_i[j],sum5_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+80ull],sum5_i);                        
                        xmm6_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+96ull]);
                        sum6_i    = _mm512_add_ps(xmm6_i[j],sum6_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+96ull],sum6_i);                        
                        xmm7_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+112ull]);
                        sum7_i    = _mm512_add_ps(xmm7_i[j],sum7_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+112ull],sum7_i);                      
                                                                               
                    }
              }

              for(; (i+63ull) < this->m_I_nsamples; i += 64ull) 
              {
                    sum0_i = vzero;                    
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                                                                                     
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);                        
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);                      
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);                       
                                                                                               
                    }
              }

              for(; (i+31ull) < this->m_I_nsamples; i += 32ull) 
              {
                    sum0_i = vzero;                    
                    sum1_i = vzero;                                                                                                       
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);                        
                                                                                                                
                    }
              }

              for(; (i+15ull) < this->m_I_nsamples; i += 16ull) 
              {
                    sum0_i = vzero;                                                                                                                        
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                                                                                                                                  
                    }
              
              }

              for(; (i+0ull) < this->m_I_nsamples; i += 1ull) 
              {
                   sum_rem_i = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_i[j] = this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i];
                       sum_rem_i     += scal_rem_i[j];
                       this->m_I_vsequence.m_data[i] = sum_rem_i;
                   }
              }

              ///////////////////////////// Q-channel //////////////////////////////////////////////
              
              for(i = 0ull; (i+159ull) < this->m_Q_nsamples; i += 160ull) 
              {
                   
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                    
                    sum6_q = vzero;                    
                    sum7_q = vzero;                    
                    sum8_q = vzero;                   
                    sum9_q = vzero;                    
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+128ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+144ull],_MM_HINT_T1);                       
#endif 
                       
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);                       
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);                        
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);                        
                        xmm4_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum4_q    = _mm512_add_ps(xmm4_q[j],sum4_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum4_q);                        
                        xmm5_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull]);
                        sum5_q    = _mm512_add_ps(xmm5_q[j],sum5_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+80ull],sum5_q);                        
                        xmm6_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull]);
                        sum6_q    = _mm512_add_ps(xmm6_q[j],sum6_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+96ull],sum6_q);                        
                        xmm7_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull]);
                        sum7_q    = _mm512_add_ps(xmm7_q[j],sum7_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+112ull],sum7_q);                        
                        xmm8_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+128ull]);
                        sum8_q    = _mm512_add_ps(xmm8_q[j],sum8_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+128ull],sum8_q);                        
                        xmm9_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+144ull]);
                        sum9_q    = _mm512_add_ps(xmm9_q[j],sum9_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+144ull],sum9_q);                        
                       
                    }
              }

              for(; (i+127ull) < this->m_Q_nsamples; i += 128ull) 
              {
                   
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                    
                    sum6_q = vzero;                    
                    sum7_q = vzero;                                                   
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                       
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);                       
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);                        
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);                        
                        xmm4_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum4_q    = _mm512_add_ps(xmm4_q[j],sum4_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum4_q);                        
                        xmm5_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull]);
                        sum5_q    = _mm512_add_ps(xmm5_q[j],sum5_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+80ull],sum5_q);                        
                        xmm6_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+96ull]);
                        sum6_q    = _mm512_add_ps(xmm6_q[j],sum6_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+96ull],sum6_q);                        
                        xmm7_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+112ull]);
                        sum7_q    = _mm512_add_ps(xmm7_q[j],sum7_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+112ull],sum7_q);                        
                                                                 
                    }
              }

              for(; (i+63ull) < this->m_Q_nsamples; i += 64ull) 
              {
                   
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                                                                     
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                       
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);                       
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);                        
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);                        
                                                                                      
                    }
              }

              for(; (i+31ull) < this->m_Q_nsamples; i += 32ull) 
              {
                   
                    sum0_q = vzero;                    
                    sum1_q = vzero;                                                                                       
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                       
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);                       
                                                                                                            
                    }
              }

              for(; (i+15ull) < this->m_Q_nsamples; i += 16ull) 
              {                  
                    sum0_q = vzero;                                                                                                                           
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                       
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                                                                                                                                  
                    }
              }
               
               for(; (i+0ull) < this->m_Q_nsamples; i += 1ull) 
               {                 
                   sum_rem_q = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_q[j] = this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i];
                       sum_rem_q     += scal_rem_q[j];
                       this->m_Q_vsequence.m_data[i] = sum_rem_q;
                   }
              }
        }
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
           HW_PMC_COLLECTION_EPILOGE_BODY
           HW_PMC_SHOW_INFO
#endif 
     return (0);
}

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SCALAR_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u10x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);

template<gms::radiolocation::iq_rectw_bitstream_vsequence_t
         ::IQ_rectw_bitstreams_optim_path  optim_path>       
std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u6x(std::int32_t * ret_codes)
{
       if(__builtin_expect(nullptr==ret_codes,0)) { return (-1);}
       if(__builtin_expect(this->m_IQ_nsignals>256ull,0)) { return (-2);}
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
           PMC_VARS
           HW_PMC_COLLECTION_PROLOGE_BODY
#endif 

        constexpr std::size_t nsignals_size{256ull};
        __ATTR_ALIGN__(64) 
        __m512 xmm0_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm1_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm2_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm3_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm4_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm5_i[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm0_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm1_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm2_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm3_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm4_q[nsignals_size];
        __ATTR_ALIGN__(64) 
        __m512 xmm5_q[nsignals_size];
        float scal_rem_i[nsignals_size];
        float scal_rem_q[nsignals_size];
        __m512 vzero{_mm512_setzero_ps()};
        __m512 sum0_i{vzero};
        __m512 sum1_i{sum0_i};
        __m512 sum2_i{sum0_i};
        __m512 sum3_i{sum0_i};
        __m512 sum4_i{sum0_i};
        __m512 sum5_i{sum0_i};
        __m512 sum0_q{sum0_i};
        __m512 sum1_q{sum0_i};
        __m512 sum2_q{sum0_i};
        __m512 sum3_q{sum0_i};
        __m512 sum4_q{sum0_i};
        __m512 sum5_q{sum0_i};
        std::size_t i;
        float sum_rem_i{0.0f};
        float sum_rem_q{0.0f};

        if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SCALAR_PATH)
        {          
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_scalar();
                }            
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_sse();
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_sse_u4x();
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx_u4x(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx_u4x(true);
                } 
        }
        else if constexpr(optim_path == IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH)
        {
                for(std::size_t k{0ull}; k != this->m_IQ_nsignals; ++k) 
                {
                    ret_codes[k] = this->m_iq_rectw_bitstreams.operator[](k).generate_I_channel_bitstream_avx512_u4x(true);
                } 
        }

        if(this->m_I_nsamples==this->m_Q_nsamples)
        {
              for(i = 0ull; (i+95ull) < this->m_I_nsamples; i += 96ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    sum4_i = vzero;
                    sum4_q = vzero;
                    sum5_i = vzero;
                    sum5_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull],_MM_HINT_T1);
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull],_MM_HINT_T1);
#endif 
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);
                        xmm4_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum4_i    = _mm512_add_ps(xmm4_i[j],sum4_i);
                        xmm4_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum4_q    = _mm512_add_ps(xmm4_q[j],sum4_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum4_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum4_q);
                        xmm5_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull]);
                        sum5_i    = _mm512_add_ps(xmm5_i[j],sum5_i);
                        xmm5_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull]);
                        sum5_q    = _mm512_add_ps(xmm5_q[j],sum5_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+80ull],sum5_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+80ull],sum5_q);
                    }
              }
              
                
              for(; (i+63ull) < this->m_I_nsamples; i += 64ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    sum2_i = vzero;
                    sum2_q = vzero;
                    sum3_i = vzero;
                    sum3_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);
                        
                    }
              }

              for(; (i+31ull) < this->m_I_nsamples; i += 32ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    sum1_i = vzero;
                    sum1_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);
                        
                    }
              }

               for(; (i+15ull) < this->m_I_nsamples; i += 16ull) 
              {
                    sum0_i = vzero;
                    sum0_q = vzero;
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);
                     
                    }
              }

                         
              for(; (i+0ull) < this->m_I_nsamples; i += 1ull) 
              {
                   sum_rem_i = 0.0f;
                   sum_rem_q = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_i[j] = this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i];
                       sum_rem_i     += scal_rem_i[j];
                       this->m_I_vsequence.m_data[i] = sum_rem_i;
                       scal_rem_q[j] = this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i];
                       sum_rem_q     += scal_rem_q[j];
                       this->m_Q_vsequence.m_data[i] = sum_rem_q;
                   }
              }

        }
        else 
        {
             for(i = 0ull; (i+95ull) < this->m_I_nsamples; i += 96ull) 
             {
                    sum0_i = vzero;                    
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                    
                    sum4_i = vzero;                    
                    sum5_i = vzero;                                                       
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull],_MM_HINT_T1);                                                               
#endif 
                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);                        
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);                      
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);                       
                        xmm4_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+64ull]);
                        sum4_i    = _mm512_add_ps(xmm4_i[j],sum4_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+64ull],sum4_i);                       
                        xmm5_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+80ull]);
                        sum5_i    = _mm512_add_ps(xmm5_i[j],sum5_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+80ull],sum5_i);                        
                                                                    
                    }
              }
              
              
              for(; (i+63ull) < this->m_I_nsamples; i += 64ull) 
              {
                    sum0_i = vzero;                    
                    sum1_i = vzero;                    
                    sum2_i = vzero;                    
                    sum3_i = vzero;                                                                                     
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);                        
                        xmm2_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+32ull]);
                        sum2_i    = _mm512_add_ps(xmm2_i[j],sum2_i);                        
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+32ull],sum2_i);                      
                        xmm3_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+48ull]);
                        sum3_i    = _mm512_add_ps(xmm3_i[j],sum3_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+48ull],sum3_i);                       
                                                                                               
                    }
              }

              for(; (i+31ull) < this->m_I_nsamples; i += 32ull) 
              {
                    sum0_i = vzero;                    
                    sum1_i = vzero;                                                                                                       
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                        xmm1_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+16ull]);
                        sum1_i    = _mm512_add_ps(xmm1_i[j],sum1_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+16ull],sum1_i);                        
                                                                                                                
                    }
              }

              for(; (i+15ull) < this->m_I_nsamples; i += 16ull) 
              {
                    sum0_i = vzero;                                                                                                                        
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {

                        xmm0_i[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i+0ull]);
                        sum0_i    = _mm512_add_ps(xmm0_i[j],sum0_i);                       
                        _mm512_store_ps(&this->m_I_vsequence.m_data[i+0ull],sum0_i);                        
                                                                                                                                  
                    }
              
              }

              for(; (i+0ull) < this->m_I_nsamples; i += 1ull) 
              {
                   sum_rem_i = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_i[j] = this->m_iq_rectw_bitstreams.operator[](j).m_I_ch_bitstream.m_data[i];
                       sum_rem_i     += scal_rem_i[j];
                       this->m_I_vsequence.m_data[i] = sum_rem_i;
                   }
              }

              ///////////////////////////// Q-channel //////////////////////////////////////////////
              
              for(i = 0ull; (i+95ull) < this->m_Q_nsamples; i += 96ull) 
              {                  
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                   
                    sum4_q = vzero;                   
                    sum5_q = vzero;                                                        
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH) == 1
                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull],_MM_HINT_T1);                        
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull],_MM_HINT_T1);                      
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull],_MM_HINT_T1);                       
                        _mm_prefetch((const char*)&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull],_MM_HINT_T1);                                          
#endif 
                       
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);                       
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);                        
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);                        
                        xmm4_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+64ull]);
                        sum4_q    = _mm512_add_ps(xmm4_q[j],sum4_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+64ull],sum4_q);                        
                        xmm5_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+80ull]);
                        sum5_q    = _mm512_add_ps(xmm5_q[j],sum5_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+80ull],sum5_q);                        
                                            
                    }
              }

              for(; (i+63ull) < this->m_Q_nsamples; i += 64ull) 
              {                   
                    sum0_q = vzero;                    
                    sum1_q = vzero;                   
                    sum2_q = vzero;                    
                    sum3_q = vzero;                                                                     
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                       
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);                       
                        xmm2_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+32ull]);
                        sum2_q    = _mm512_add_ps(xmm2_q[j],sum2_q);                        
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+32ull],sum2_q);                        
                        xmm3_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+48ull]);
                        sum3_q    = _mm512_add_ps(xmm3_q[j],sum3_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+48ull],sum3_q);                        
                                                                                      
                    }
              }

              for(; (i+31ull) < this->m_Q_nsamples; i += 32ull) 
              {                   
                    sum0_q = vzero;                    
                    sum1_q = vzero;                                                                                       
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                       
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                        xmm1_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+16ull]);
                        sum1_q    = _mm512_add_ps(xmm1_q[j],sum1_q);                       
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+16ull],sum1_q);                       
                                                                                                            
                    }
              }

              for(; (i+15ull) < this->m_Q_nsamples; i += 16ull) 
              {                  
                    sum0_q = vzero;                                                                                                                           
                    for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)  
                    {
                       
                        xmm0_q[j] = _mm512_load_ps(&this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i+0ull]);
                        sum0_q    = _mm512_add_ps(xmm0_q[j],sum0_q);                      
                        _mm512_store_ps(&this->m_Q_vsequence.m_data[i+0ull],sum0_q);                       
                                                                                                                                  
                    }
              }
               
               for(; (i+0ull) < this->m_Q_nsamples; i += 1ull) 
               {                 
                   sum_rem_q = 0.0f;
                   for(std::size_t j{0ull}; j != this->m_IQ_nsignals; ++j)   
                   {
                       scal_rem_q[j] = this->m_iq_rectw_bitstreams.operator[](j).m_Q_ch_bitstream.m_data[i];
                       sum_rem_q     += scal_rem_q[j];
                       this->m_Q_vsequence.m_data[i] = sum_rem_q;
                   }
              }
        }
#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
           HW_PMC_COLLECTION_EPILOGE_BODY
           HW_PMC_SHOW_INFO
#endif 
     return (0);
}

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SCALAR_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::SSE_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_VEC_PATH
                                              >(std::int32_t * __restrict__);

template std::int32_t 
gms::radiolocation
::iq_rectw_bitstream_vsequence_t
::generate_iq_rectw_bitstream_sequence_avx512_u6x<gms::radiolocation::
                                              iq_rectw_bitstream_vsequence_t::
                                              IQ_rectw_bitstreams_optim_path::AVX512_U4X_VEC_PATH
                                              >(std::int32_t * __restrict__);
