
#include <fstream>
#include <iomanip>
#include <immintrin.h>
#include "GMS_iq_rectwave_bitstream.h"

gms::radiolocation
::iq_rectwave_bitstream_t
::iq_rectwave_bitstream_t() noexcept(true) 
:
m_I_k{0.0f},
m_I_T{0.0f},
m_I_nsamples{0ull},
m_I_Ac{0.0f},
m_I_fc{0.0f},
m_I_ph0{0.0f},
m_Q_k{0.0f},
m_Q_T{0.0f},
m_Q_nsamples{0ull},
m_Q_Ac{0.0f},
m_Q_fc{0.0f},
m_Q_ph0{0.0f},
m_I_ch_bitstream{darray_r4_t()},
m_Q_ch_bitstream{darray_r4_t()}
{

}

gms::radiolocation
::iq_rectwave_bitstream_t
::iq_rectwave_bitstream_t(const float I_k,
                          const float I_T,
                          const float I_Ac,
                          const float I_fc,
                          const float I_ph0,
                          const float Q_k,
                          const float Q_T,
                          const float Q_Ac,
                          const float Q_fc,
                          const float Q_ph0) noexcept(true) 
:
m_I_k{I_k},
m_I_T{I_T},
m_I_nsamples{static_cast<std::size_t>(this->m_I_k*this->m_I_T)},
m_I_Ac{I_Ac},
m_I_fc{I_fc},
m_I_ph0{I_ph0},
m_Q_k{Q_k},
m_Q_T{Q_T},
m_Q_nsamples{static_cast<std::size_t>(this->m_Q_k*this->m_Q_T)},
m_Q_Ac{Q_Ac},
m_Q_fc{Q_fc},
m_Q_ph0{Q_ph0},
m_I_ch_bitstream{darray_r4_t(this->m_I_nsamples)},
m_Q_ch_bitstream{darray_r4_t(this->m_Q_nsamples)}
{

}

gms::radiolocation 
::iq_rectwave_bitstream_t
::iq_rectwave_bitstream_t(const iq_rectwave_bitstream_t &other) noexcept(true) 
:
m_I_k{other.m_I_k},
m_I_T{other.m_I_T},
m_I_nsamples{other.m_I_nsamples},
m_I_Ac{other.m_I_Ac},
m_I_fc{other.m_I_fc},
m_I_ph0{other.m_I_ph0},
m_Q_k{other.m_Q_k},
m_Q_T{other.m_Q_T},
m_Q_nsamples{other.m_Q_nsamples},
m_Q_Ac{other.m_Q_Ac},
m_Q_fc{other.m_Q_fc},
m_Q_ph0{other.m_Q_ph0},
m_I_ch_bitstream(other.m_I_ch_bitstream),
m_Q_ch_bitstream(other.m_Q_ch_bitstream)
{

}

gms::radiolocation 
::iq_rectwave_bitstream_t
::iq_rectwave_bitstream_t(iq_rectwave_bitstream_t &&other) noexcept(true) 
:
m_I_k{other.m_I_k},
m_I_T{other.m_I_T},
m_I_nsamples{other.m_I_nsamples},
m_I_Ac{other.m_I_Ac},
m_I_fc{other.m_I_fc},
m_I_ph0{other.m_I_ph0},
m_Q_k{other.m_Q_k},
m_Q_T{other.m_Q_T},
m_Q_nsamples{other.m_Q_nsamples},
m_Q_Ac{other.m_Q_Ac},
m_Q_fc{other.m_Q_fc},
m_Q_ph0{other.m_Q_ph0},
m_I_ch_bitstream(std::move(other.m_I_ch_bitstream)),
m_Q_ch_bitstream(std::move(other.m_Q_ch_bitstream))
{
    
}

gms::radiolocation 
::iq_rectwave_bitstream_t
::~iq_rectwave_bitstream_t() noexcept(true) 
{

}

gms::radiolocation 
::iq_rectwave_bitstream_t &
gms::radiolocation 
::iq_rectwave_bitstream_t
::operator=(const iq_rectwave_bitstream_t &other) noexcept(true) 
{
    if(__builtin_expect(this==&other)) { return (*this);}
    this->m_I_k        = other.m_I_k;
    this->m_I_T        = other.m_I_T;
    this->m_I_nsamples = other.m_I_nsamples;
    this->m_I_Ac       = other.m_I_Ac;
    this->m_I_fc       = other.m_I_fc;
    this->m_I_ph0      = other.m_I_ph0;
    this->m_Q_k        = other.m_Q_k;
    this->m_Q_T        = other.m_Q_T;
    this->m_Q_nsamples = other.m_Q_nsamples;
    this->m_Q_Ac       = other.m_Q_Ac;
    this->m_Q_fc       = other.m_Q_fc;
    this->m_Q_ph0      = other.m_Q_ph0;
    this->m_I_ch_bitstream.operator=(other.m_I_ch_bitstream);
    this->m_Q_ch_bitstream.operator=(other.m_Q_ch_bitstream);
    return (*this);
}

gms::radiolocation 
::iq_rectwave_bitstream_t &
gms::radiolocation 
::iq_rectwave_bitstream_t
::operator=(iq_rectwave_bitstream_t &&other) noexcept(true) 
{
    if(__builtin_expect(this==&other)) { return (*this);}
    this->m_I_k        = other.m_I_k;
    this->m_I_T        = other.m_I_T;
    this->m_I_nsamples = other.m_I_nsamples;
    this->m_I_Ac       = other.m_I_Ac;
    this->m_I_fc       = other.m_I_fc;
    this->m_I_ph0      = other.m_I_ph0;
    this->m_Q_k        = other.m_Q_k;
    this->m_Q_T        = other.m_Q_T;
    this->m_Q_nsamples = other.m_Q_nsamples;
    this->m_Q_Ac       = other.m_Q_Ac;
    this->m_Q_fc       = other.m_Q_fc;
    this->m_Q_ph0      = other.m_Q_ph0;
    this->m_I_ch_bitstream.operator=(std::move(other.m_I_ch_bitstream));
    this->m_Q_ch_bitstream.operator=(std::move(other.m_Q_ch_bitstream));
    return (*this);
}

void
gms::radiolocation
::iq_rectwave_bitstream_t
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
::iq_rectwave_bitstream_t
::generate_I_channel_bitstream_scalar()
{
     using namespace gms::math;
     constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
     const float invT{1.0f/this->m_I_nsamples};
     for(std::size_t i{0ull}; i != this->m_I_nsamples; ++i) 
     {
         const float t_i{static_cast<float>(i*invT)};
#if (IQ_RECTWAVE_BITSTREAM_USE_CEPHES) == 1
         const float cos_val{this->m_I_Ac*ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_i))};
         this->m_I_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
#else 
         const float cos_val{this->m_I_Ac*std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_i))};
         this->m_I_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif 
     }
     return (0);
}

std::int32_t 
gms::radiolocation 
::iq_rectwave_bitstream_t
::generate_Q_channel_bitstream_scalar()
{
     using namespace gms::math;
     constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
     const float invT{1.0f/this->m_Q_nsamples};
     for(std::size_t i{0ull}; i != this->m_Q_nsamples; ++i) 
     {
         const float t_i{static_cast<float>(i*invT)};
#if (IQ_RECTWAVE_BITSTREAM_USE_CEPHES) == 1
         const float sin_val{this->m_Q_Ac*ceph_sinf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_i))};
         this->m_Q_ch_bitstream.m_data[i] = (sin_val>=0.0f)?+1.0f:-1.0f;
#else 
         const float sin_val{this->m_Q_Ac*std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_i))};
         this->m_Q_ch_bitstream.m_data[i] = (sin_val>=0.0f)?+1.0f:-1.0f;
#endif 
     }
     return (0);
}

std::int32_t 
gms::radiolocation
::iq_rectwave_bitstream_t
::generate_I_channel_bitstream_sse()
{
     using namespace gms::math;
     std::size_t i,j;
     const float inv_sr{1.0f/this->m_I_nsamples};
     const __m128 vinv_sr{_mm_set1_ps(inv_sr)};
     const __m128 vfc{_mm_set1_ps(this->m_I_fc)};
     const __m128 vph0{_mm_set1_ps(this->m_I_ph0)};
     const __m128 vAc{_mm_set1_ps(this->m_I_Ac)};
     const __m128 v2pi{_mm_set1_ps(6.283185307179586476925286766559f)};
     const __m128 vzero{_mm_setzero_ps()};
     const __m128 vpone{_mm_set1_ps(+1.0f)};
     const __m128 vnone{_mm_set1_ps(-1.0f)};
     float jj;
     constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
     constexpr std::size_t LUT_loop_idx_threshold{2257ull};
     if(this->m_I_nsamples>LUT_loop_idx_threshold)
     {
        for(i = 0ull,jj = 0.0f; i != ROUND_TO_FOUR(this->m_I_nsamples,4ull); i += 4ull,jj += 4.0f) 
        {
            const __m128 vt_i{_mm_setr_ps(jj*inv_sr,(jj+1.0f)*inv_sr,(jj+2.0f)*inv_sr,(jj+3.0f)*inv_sr)};
            const __m128 vcos_val{_mm_mul_ps(vAc,_mm_cos_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i),vph0)))};
            const __mmask8 vcos_ge_0{_mm_cmp_ps_mask(vcos_val,vzero,_CMP_GE_OQ)};
            _mm_store_ps(&this->m_I_ch_bitstream.m_data[i], _mm_mask_blend_ps(vcos_ge_0,vnone,vpone));   
        }
                       
        for(j = i; j != this->m_I_nsamples; ++j) 
        {
            const float t_j{static_cast<float>(j*inv_sr)};
#if (IQ_RECTWAVE_BITSTREAM_USE_CEPHES) == 1 
            const float cos_val = this->m_I_Ac*ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_j));
            this->m_I_ch_bitstream.m_data[j] = (cos_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
            const float cos_val = this->m_I_Ac*std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_j));
            this->m_I_ch_bitstream.m_data[j] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif
        }
    }
    else 
    {
        for(i = 0ull; i != ROUND_TO_FOUR(this->m_I_nsamples,4ull); i += 4ull) 
        {
            _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i],_MM_HINT_T0);
            const __m128 vt_i{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]),vinv_sr)};
            const __m128 vcos_val{_mm_mul_ps(this->m_I_Ac,_mm_cos_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i),vph0)))};
            const __mmask8 vcos_ge_0{_mm_cmp_ps_mask(vcos_val,vzero,_CMP_GE_OQ)};
            _mm_store_ps(&this->m_I_ch_bitstream.m_data[i], _mm_mask_blend_ps(vcos_ge_0,vnone,vpone));   
        }
                       
        for(j = i; j != this->m_I_nsamples; ++j) 
        {
            const float t_j{gms::math::LUT_loop_indices_2257_align16[j]*inv_sr};
#if (IQ_RECTWAVE_BITSTREAM_USE_CEPHES) == 1 
            const float cos_val = this->m_I_Ac*ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_j));
            this->m_I_ch_bitstream.m_data[j] = (cos_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
            const float cos_val = this->m_I_Ac*std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_j));
            this->m_I_ch_bitstream.m_data[j] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif
        }
    }

    return (0);
}

std::int32_t 
gms::radiolocation
::iq_rectwave_bitstream_t
::generate_Q_channel_bitstream_sse()
{
     using namespace gms::math;
     std::size_t i,j;
     const float inv_sr{1.0f/this->m_Q_nsamples};
     const __m128 vinv_sr{_mm_set1_ps(inv_sr)};
     const __m128 vfc{_mm_set1_ps(this->m_Q_fc)};
     const __m128 vph0{_mm_set1_ps(this->m_Q_ph0)};
     const __m128 vAc{_mm_set1_ps(this->m_Q_Ac)};
     const __m128 v2pi{_mm_set1_ps(6.283185307179586476925286766559f)};
     const __m128 vzero{_mm_setzero_ps()};
     const __m128 vpone{_mm_set1_ps(+1.0f)};
     const __m128 vnone{_mm_set1_ps(-1.0f)};
     float jj;
     constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
     constexpr std::size_t LUT_loop_idx_threshold{2257ull};
     if(this->m_I_nsamples>LUT_loop_idx_threshold)
     {
        for(i = 0ull,jj = 0.0f; i != ROUND_TO_FOUR(this->m_Q_nsamples,4ull); i += 4ull,jj += 4.0f) 
        {
            const __m128 vt_i{_mm_setr_ps(jj*inv_sr,(jj+1.0f)*inv_sr,(jj+2.0f)*inv_sr,(jj+3.0f)*inv_sr)};
            const __m128 vsin_val{_mm_mul_ps(vAc,_mm_sin_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i),vph0)))};
            const __mmask8 vsin_ge_0{_mm_cmp_ps_mask(vsin_val,vzero,_CMP_GE_OQ)};
            _mm_store_ps(&this->m_Q_ch_bitstream.m_data[i], _mm_mask_blend_ps(vsin_ge_0,vnone,vpone));   
        }
                       
        for(j = i; j != this->m_Q_nsamples; ++j) 
        {
            const float t_j{static_cast<float>(j*inv_sr)};
#if (IQ_RECTWAVE_BITSTREAM_USE_CEPHES) == 1 
            const float sin_val = this->m_Q_Ac*ceph_sinf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_j));
            this->m_Q_ch_bitstream.m_data[j] = (sin_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
            const float sin_val = this->m_Q_Ac*std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_j));
            this->m_Q_ch_bitstream.m_data[j] = (sin_val>=0.0f)?+1.0f:-1.0f;
#endif
        }
    }
    else 
    {
        for(i = 0ull; i != ROUND_TO_FOUR(this->m_Q_nsamples,4ull); i += 4ull) 
        {
            _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i],_MM_HINT_T0);
            const __m128 vt_i{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]),vinv_sr)};
            const __m128 vsin_val{_mm_mul_ps(this->m_Q_Ac,_mm_sin_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i),vph0)))};
            const __mmask8 vsin_ge_0{_mm_cmp_ps_mask(vsin_val,vzero,_CMP_GE_OQ)};
            _mm_store_ps(&this->m_Q_ch_bitstream.m_data[i], _mm_mask_blend_ps(vsin_ge_0,vnone,vpone));   
        }
                       
        for(j = i; j != this->m_Q_nsamples; ++j) 
        {
            const float t_j{gms::math::LUT_loop_indices_2257_align16[j]*inv_sr};
#if (IQ_RECTWAVE_BITSTREAM_USE_CEPHES) == 1 
            const float sin_val = this->m_Q_Ac*ceph_sinf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_j));
            this->m_Q_ch_bitstream.m_data[j] = (sin_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
            const float sin_val = this->m_Q_Ac*std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_j));
            this->m_Q_ch_bitstream.m_data[j] = (sin_val>=0.0f)?+1.0f:-1.0f;
#endif
        }
    }

    return (0);
}





