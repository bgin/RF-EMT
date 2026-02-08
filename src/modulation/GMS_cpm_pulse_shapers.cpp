
#include <fstream>
#include <iomanip>
#include <immintrin.h>
#include "GMS_cpm_pulse_shapers.h"

gms::radiolocation
::cpm_pulse_shapers_t
::cpm_pulse_shapers_t() noexcept(true)
:
m_nTLsamples{0ull},
m_nTsamples{0ull},
m_T{0.0f},
m_L{0.0f},
m_beta{0.0f},
m_BbT{0.0f},
m_lrec_pulse{darray_r4_t()},
m_lrc_pulse{darray_r4_t()},
m_lsrc_pulse{darray_r4_t()},
m_gmsk_pulse{darray_r4_t()}
{

}

gms::radiolocation
::cpm_pulse_shapers_t
::cpm_pulse_shapers_t(const std::size_t nTLsamples,
                      const std::size_t nTsamples,
                      const float T,
                      const float L,
                      const float beta,
                      const float BbT) noexcept(false)
:
m_nTLsamples{nTLsamples},
m_nTsamples{nTsamples},
m_T{T},
m_L{L},
m_beta{beta},
m_BbT{BbT},
m_lrec_pulse{darray_r4_t(this->m_nTLsamples)},
m_lrc_pulse{darray_r4_t(this->m_nTLsamples)},
m_lsrc_pulse{darray_r4_t(this->m_nTLsamples)},
m_gmsk_pulse{darray_r4_t(this->m_nTsamples)}
{

}

gms::radiolocation 
::cpm_pulse_shapers_t
::cpm_pulse_shapers_t(const cpm_pulse_shapers_t &other) noexcept(false)
:
m_nTLsamples{other.m_nTLsamples},
m_nTsamples{other.m_nTsamples},
m_T{other.m_T},
m_L{other.m_L},
m_beta{other.m_beta},
m_BbT{other.m_BbT},
m_lrec_pulse{other.m_lrec_pulse},
m_lrc_pulse{other.m_lrc_pulse},
m_lsrc_pulse{other.m_lsrc_pulse},
m_gmsk_pulse{other.m_gmsk_pulse} 
{

}

gms::radiolocation 
::cpm_pulse_shapers_t
::cpm_pulse_shapers_t(cpm_pulse_shapers_t &&other) noexcept(true)
:
m_nTLsamples{other.m_nTLsamples},
m_nTsamples{other.m_nTsamples},
m_T{other.m_T},
m_L{other.m_L},
m_beta{other.m_beta},
m_BbT{other.m_BbT},
m_lrec_pulse{std::move(other.m_lrec_pulse)},
m_lrc_pulse{std::move(other.m_lrc_pulse)},
m_lsrc_pulse{std::move(other.m_lsrc_pulse)},
m_gmsk_pulse{std::move(other.m_gmsk_pulse)} 
{

}

gms::radiolocation
::cpm_pulse_shapers_t
::~cpm_pulse_shapers_t() noexcept(true) 
{

}

gms::radiolocation
::cpm_pulse_shapers_t &
gms::radiolocation
::cpm_pulse_shapers_t
::operator=(const cpm_pulse_shapers_t &other) noexcept(false) 
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_nTLsamples = other.m_nTLsamples;
    this->m_nTsamples  = other.m_nTsamples;
    this->m_T          = other.m_T;
    this->m_L          = other.m_L;
    this->m_beta       = other.m_beta;
    this->m_BbT        = other.m_BbT; 
    this->m_lrec_pulse.operator=(other.m_lrec_pulse);
    this->m_lrc_pulse.operator=(other.m_lrc_pulse);
    this->m_lsrc_pulse.operator=(other.m_lsrc_pulse);
    this->m_gmsk_pulse.operator=(other.m_gmsk_pulse);
    return (*this);
}

gms::radiolocation
::cpm_pulse_shapers_t &
gms::radiolocation
::cpm_pulse_shapers_t
::operator=(cpm_pulse_shapers_t &&other) noexcept(true) 
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_nTLsamples = other.m_nTLsamples;
    this->m_nTsamples  = other.m_nTsamples;
    this->m_T          = other.m_T;
    this->m_L          = other.m_L;
    this->m_beta       = other.m_beta;
    this->m_BbT        = other.m_BbT; 
    this->m_lrec_pulse.operator=(std::move(other.m_lrec_pulse));
    this->m_lrc_pulse.operator=(std::move(other.m_lrc_pulse));
    this->m_lsrc_pulse.operator=(std::move(other.m_lsrc_pulse));
    this->m_gmsk_pulse.operator=(std::move(other.m_gmsk_pulse));
    return (*this);
}

void 
gms::radiolocation 
::cpm_pulse_shapers_t
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
::cpm_pulse_shapers_t
::generate_lrc_pulse_sse_u8x() noexcept(true)
{
    using namespace gms::math;
    constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
    constexpr std::size_t LUT_loop_idx_threshold{2257ull};
    const float twoLT{2.0f*this->m_L*this->m_T};
    const float inv2LT{1.0f/twoLT};
    const float invLT{1.0f/(this->m_L*this->m_T)};
    const __m128 vC6283185307179586476925286766559{_mm_set1_ps(C6283185307179586476925286766559)};
    const __m128 vinv2LT{_mm_set1_ps(inv2LT)};
    const __m128 vinvLT{_mm_set1_ps(invLT)};
    const __m128 vone{_mm_set1_ps(1.0f)};
    std::size_t i;
    float j;
    if(__builtin_expect(this->m_nTLsamples>LUT_loop_idx_threshold,0)) 
    {
         for(i = 0ull,j = 0.0f; (i+31ull) < this->m_nTLsamples; i += 32ull,j += 32.0f) 
         {
                 const __m128 vt_i_0{_mm_setr_ps(j*invLT,(j+1.0f)*invLT,(j+2.0f)*invLT,(j+3.0f)*invLT)};
                 const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0))};
                 const __m128 vlrc_value_0{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_0))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value_0);
                 const __m128 vt_i_1{_mm_setr_ps((j+4.0f)*invLT,(j+5.0f)*invLT,(j+6.0f)*invLT,(j+7.0f)*invLT)};
                 const __m128 vcos_term_1{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_1))};
                 const __m128 vlrc_value_1{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_1))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+4ull],vlrc_value_1);
                 const __m128 vt_i_2{_mm_setr_ps((j+8.0f)*invLT,(j+9.0f)*invLT,(j+10.0f)*invLT,(j+11.0f)*invLT)};
                 const __m128 vcos_term_2{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_2))};
                 const __m128 vlrc_value_2{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_2))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+8ull],vlrc_value_2);
                 const __m128 vt_i_3{_mm_setr_ps((j+12.0f)*invLT,(j+13.0f)*invLT,(j+14.0f)*invLT,(j+15.0f)*invLT)};
                 const __m128 vcos_term_3{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_3))};
                 const __m128 vlrc_value_3{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_3))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+12ull],vlrc_value_3);
                 const __m128 vt_i_4{_mm_setr_ps((j+16.0f)*invLT,(j+17.0f)*invLT,(j+18.0f)*invLT,(j+19.0f)*invLT)};
                 const __m128 vcos_term_4{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_4))};
                 const __m128 vlrc_value_4{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_4))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+16ull],vlrc_value_4);
                 const __m128 vt_i_5{_mm_setr_ps((j+20.0f)*invLT,(j+21.0f)*invLT,(j+22.0f)*invLT,(j+23.0f)*invLT)};
                 const __m128 vcos_term_5{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_5))};
                 const __m128 vlrc_value_5{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_5))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+20ull],vlrc_value_5);
                 const __m128 vt_i_6{_mm_setr_ps((j+24.0f)*invLT,(j+25.0f)*invLT,(j+26.0f)*invLT,(j+27.0f)*invLT)};
                 const __m128 vcos_term_6{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_6))};
                 const __m128 vlrc_value_6{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_6))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+24ull],vlrc_value_6);
                 const __m128 vt_i_7{_mm_setr_ps((j+28.0f)*invLT,(j+29.0f)*invLT,(j+30.0f)*invLT,(j+31.0f)*invLT)};
                 const __m128 vcos_term_7{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_7))};
                 const __m128 vlrc_value_7{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_7))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+28ull],vlrc_value_7);
         }

         for(; (i+23ull) < this->m_nTLsamples; i += 24ull,j += 24.0f) 
         {
                 const __m128 vt_i_0{_mm_setr_ps(j*invLT,(j+1.0f)*invLT,(j+2.0f)*invLT,(j+3.0f)*invLT)};
                 const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0))};
                 const __m128 vlrc_value_0{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_0))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value_0);
                 const __m128 vt_i_1{_mm_setr_ps((j+4.0f)*invLT,(j+5.0f)*invLT,(j+6.0f)*invLT,(j+7.0f)*invLT)};
                 const __m128 vcos_term_1{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_1))};
                 const __m128 vlrc_value_1{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_1))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+4ull],vlrc_value_1);
                 const __m128 vt_i_2{_mm_setr_ps((j+8.0f)*invLT,(j+9.0f)*invLT,(j+10.0f)*invLT,(j+11.0f)*invLT)};
                 const __m128 vcos_term_2{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_2))};
                 const __m128 vlrc_value_2{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_2))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+8ull],vlrc_value_2);
                 const __m128 vt_i_3{_mm_setr_ps((j+12.0f)*invLT,(j+13.0f)*invLT,(j+14.0f)*invLT,(j+15.0f)*invLT)};
                 const __m128 vcos_term_3{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_3))};
                 const __m128 vlrc_value_3{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_3))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+12ull],vlrc_value_3);
                 const __m128 vt_i_4{_mm_setr_ps((j+16.0f)*invLT,(j+17.0f)*invLT,(j+18.0f)*invLT,(j+19.0f)*invLT)};
                 const __m128 vcos_term_4{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_4))};
                 const __m128 vlrc_value_4{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_4))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+16ull],vlrc_value_4);
                 const __m128 vt_i_5{_mm_setr_ps((j+20.0f)*invLT,(j+21.0f)*invLT,(j+22.0f)*invLT,(j+23.0f)*invLT)};
                 const __m128 vcos_term_5{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_5))};
                 const __m128 vlrc_value_5{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_5))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+20ull],vlrc_value_5);
         }

         for(; (i+15ull) < this->m_nTLsamples; i += 16ull,j += 16.0f) 
         {
                 const __m128 vt_i_0{_mm_setr_ps(j*invLT,(j+1.0f)*invLT,(j+2.0f)*invLT,(j+3.0f)*invLT)};
                 const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0))};
                 const __m128 vlrc_value_0{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_0))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value_0);
                 const __m128 vt_i_1{_mm_setr_ps((j+4.0f)*invLT,(j+5.0f)*invLT,(j+6.0f)*invLT,(j+7.0f)*invLT)};
                 const __m128 vcos_term_1{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_1))};
                 const __m128 vlrc_value_1{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_1))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+4ull],vlrc_value_1);
                 const __m128 vt_i_2{_mm_setr_ps((j+8.0f)*invLT,(j+9.0f)*invLT,(j+10.0f)*invLT,(j+11.0f)*invLT)};
                 const __m128 vcos_term_2{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_2))};
                 const __m128 vlrc_value_2{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_2))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+8ull],vlrc_value_2);
                 const __m128 vt_i_3{_mm_setr_ps((j+12.0f)*invLT,(j+13.0f)*invLT,(j+14.0f)*invLT,(j+15.0f)*invLT)};
                 const __m128 vcos_term_3{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_3))};
                 const __m128 vlrc_value_3{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_3))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+12ull],vlrc_value_3); 
         }

         for(; (i+7ull) < this->m_nTLsamples; i += 8ull,j += 8.0f) 
         {
                 const __m128 vt_i_0{_mm_setr_ps(j*invLT,(j+1.0f)*invLT,(j+2.0f)*invLT,(j+3.0f)*invLT)};
                 const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0))};
                 const __m128 vlrc_value_0{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_0))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value_0);
                 const __m128 vt_i_1{_mm_setr_ps((j+4.0f)*invLT,(j+5.0f)*invLT,(j+6.0f)*invLT,(j+7.0f)*invLT)};
                 const __m128 vcos_term_1{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_1))};
                 const __m128 vlrc_value_1{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_1))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+4ull],vlrc_value_1); 
         }

         for(; (i+3ull) < this->m_nTLsamples; i += 4ull,j += 4.0f) 
         {
                 const __m128 vt_i_0{_mm_setr_ps(j*invLT,(j+1.0f)*invLT,(j+2.0f)*invLT,(j+3.0f)*invLT)};
                 const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0))};
                 const __m128 vlrc_value_0{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_0))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value_0);
         }

         for(; (i+0ull) < this->m_nTLsamples; i += 1ull,j += 1.0f) 
         {
                const float t_i{j*invLT};
#if (CPM_PULSE_SHAPERS_USE_CEPHES)
                const float cos_term_0{ceph_cosf(C6283185307179586476925286766559*t_i*invLT)};
#else
                const float cos_term_0{std::cos(C6283185307179586476925286766559*t_i*invLT)};
#endif
                const float lrc_val_0{inv2LT*(1.0f-cos_term_0)};
                this->m_lrc_pulse.m_data[i] = lrc_val_0;
         }
    }
    else 
    {
          for(i = 0ull; (i+31ull) < this->m_nTLsamples; i += 32ull) 
          {
                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
                 const __m128 vt_i_0{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinvLT)};
                 const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0))};
                 const __m128 vlrc_value_0{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_0))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value_0);
                 const __m128 vt_i_1{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull]),vinvLT)};
                 const __m128 vcos_term_1{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_1))};
                 const __m128 vlrc_value_1{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_1))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+4ull],vlrc_value_1);
                 const __m128 vt_i_2{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]),vinvLT)};
                 const __m128 vcos_term_2{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_2))};
                 const __m128 vlrc_value_2{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_2))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+8ull],vlrc_value_2);
                 const __m128 vt_i_3{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+12ull]),vinvLT)};
                 const __m128 vcos_term_3{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_3))};
                 const __m128 vlrc_value_3{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_3))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+12ull],vlrc_value_3);
                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+16ull],_MM_HINT_T0);
                 const __m128 vt_i_4{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+16ull]),vinvLT)};
                 const __m128 vcos_term_4{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_4))};
                 const __m128 vlrc_value_4{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_4))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+16ull],vlrc_value_4);
                 const __m128 vt_i_5{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+20ull]),vinvLT)};
                 const __m128 vcos_term_5{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_5))};
                 const __m128 vlrc_value_5{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_5))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+20ull],vlrc_value_5);
                 const __m128 vt_i_6{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+24ull]),vinvLT)};
                 const __m128 vcos_term_6{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_6))};
                 const __m128 vlrc_value_6{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_6))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+24ull],vlrc_value_6);
                 const __m128 vt_i_7{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+28ull]),vinvLT)};
                 const __m128 vcos_term_7{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_7))};
                 const __m128 vlrc_value_7{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_7))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+28ull],vlrc_value_7);
          }

          for(; (i+23ull) < this->m_nTLsamples; i += 24ull) 
          {
                  _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
                 const __m128 vt_i_0{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinvLT)};
                 const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0))};
                 const __m128 vlrc_value_0{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_0))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value_0);
                 const __m128 vt_i_1{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull]),vinvLT)};
                 const __m128 vcos_term_1{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_1))};
                 const __m128 vlrc_value_1{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_1))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+4ull],vlrc_value_1);
                 const __m128 vt_i_2{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]),vinvLT)};
                 const __m128 vcos_term_2{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_2))};
                 const __m128 vlrc_value_2{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_2))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+8ull],vlrc_value_2);
                 const __m128 vt_i_3{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+12ull]),vinvLT)};
                 const __m128 vcos_term_3{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_3))};
                 const __m128 vlrc_value_3{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_3))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+12ull],vlrc_value_3);
                 const __m128 vt_i_4{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+16ull]),vinvLT)};
                 const __m128 vcos_term_4{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_4))};
                 const __m128 vlrc_value_4{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_4))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+16ull],vlrc_value_4);
                 const __m128 vt_i_5{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+20ull]),vinvLT)};
                 const __m128 vcos_term_5{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_5))};
                 const __m128 vlrc_value_5{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_5))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+20ull],vlrc_value_5);
          }

          for(; (i+15ull) < this->m_nTLsamples; i += 16ull) 
          {
                  _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
                 const __m128 vt_i_0{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinvLT)};
                 const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0))};
                 const __m128 vlrc_value_0{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_0))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value_0);
                 const __m128 vt_i_1{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull]),vinvLT)};
                 const __m128 vcos_term_1{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_1))};
                 const __m128 vlrc_value_1{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_1))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+4ull],vlrc_value_1);
                 const __m128 vt_i_2{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]),vinvLT)};
                 const __m128 vcos_term_2{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_2))};
                 const __m128 vlrc_value_2{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_2))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+8ull],vlrc_value_2);
                 const __m128 vt_i_3{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+12ull]),vinvLT)};
                 const __m128 vcos_term_3{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_3))};
                 const __m128 vlrc_value_3{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_3))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+12ull],vlrc_value_3);
          }

          for(; (i+7ull) < this->m_nTLsamples; i += 8ull) 
          {
                 const __m128 vt_i_0{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinvLT)};
                 const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0))};
                 const __m128 vlrc_value_0{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_0))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value_0);
                 const __m128 vt_i_1{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull]),vinvLT)};
                 const __m128 vcos_term_1{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_1))};
                 const __m128 vlrc_value_1{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_1))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+4ull],vlrc_value_1); 
          }

          for(; (i+3ull) < this->m_nTLsamples; i += 4ull) 
          {
                 const __m128 vt_i_0{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinvLT)};
                 const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0))};
                 const __m128 vlrc_value_0{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term_0))};
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value_0); 
          }

          for(; (i+0ull) < this->m_nTLsamples; ++i) 
          {
                const float t_i{gms::math::LUT_loop_indices_2257_align16[i]};
#if (CPM_PULSE_SHAPERS_USE_CEPHES)
                const float cos_term_0{ceph_cosf(C6283185307179586476925286766559*t_i*invLT)};
#else
                const float cos_term_0{std::cos(C6283185307179586476925286766559*t_i*invLT)};
#endif
                const float lrc_val_0{inv2LT*(1.0f-cos_term_0)};
                this->m_lrc_pulse.m_data[i] = lrc_val_0;
          }
    }
    return (0);
}

std::int32_t 
gms::radiolocation 
::cpm_pulse_shapers_t
::generate_lsrc_pulse_sse_u8x() noexcept(true) 
{
        using namespace gms::math; 
        constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
        constexpr std::size_t LUT_loop_idx_threshold{2257ull};
        std::size_t i;
        float       jj;
        const float invLT{1.0f/(this->m_L*this->m_T)};
        const __m128 vinvLT{_mm_set1_ps(invLT)};
        const __m128 vC6283185307179586476925286766559{_mm_set1_ps(C6283185307179586476925286766559)};
        const __m128 vone{_mm_set1_ps(1.0f)};
        const __m128 vfour{_mm_set1_ps(4.0f)};
        const __m128 vbeta{_mm_set1_ps(this->m_beta)};
        if(__builtin_expect(this->m_nTLsamples>LUT_loop_idx_threshold,0)) 
        {
                for(i = 0ull,jj = 0.0f; (i+31ull) < this->m_nTLsamples; i += 32ull,jj = 32.0f) 
                {
                        const __m128 vt_i_0{_mm_setr_ps(jj*invLT,(jj+1.0f)*invLT,(jj+2.0f)*invLT,(jj+3.0f)*invLT)};
                        const __m128 varg_term_0{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0)};
                        const __m128 vcosinc_denom_0{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_0))};
                        const __m128 vsin_term_0{_mm_sin_ps(varg_term_0)};
                        const __m128 vsinc_term_0{_mm_div_ps(vsin_term_0,varg_term_0)};
                        const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_0))};
                        const __m128 vcosinc_term_0{_mm_div_ps(vcos_term_0,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_0,vcosinc_denom_0)))};
                        const __m128 vlsrc_sample_0{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_0,vcosinc_term_0))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+0ull], vlsrc_sample_0);
                        const __m128 vt_i_1{_mm_setr_ps((jj+4.0f)*invLT,(jj+5.0f)*invLT,(jj+6.0f)*invLT,(jj+7.0f)*invLT)};
                        const __m128 varg_term_1{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_1)};
                        const __m128 vcosinc_denom_1{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_1))};
                        const __m128 vsin_term_1{_mm_sin_ps(varg_term_1)};
                        const __m128 vsinc_term_1{_mm_div_ps(vsin_term_1,varg_term_1)};
                        const __m128 vcos_term_1{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_1))};
                        const __m128 vcosinc_term_1{_mm_div_ps(vcos_term_1,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_1,vcosinc_denom_1)))};
                        const __m128 vlsrc_sample_1{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_1,vcosinc_term_1))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+4ull], vlsrc_sample_1);
                        const __m128 vt_i_2{_mm_setr_ps((jj+8.0f)*invLT,(jj+9.0f)*invLT,(jj+10.0f)*invLT,(jj+11.0f)*invLT)};
                        const __m128 varg_term_2{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_2)};
                        const __m128 vcosinc_denom_2{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_2))};
                        const __m128 vsin_term_2{_mm_sin_ps(varg_term_2)};
                        const __m128 vsinc_term_2{_mm_div_ps(vsin_term_2,varg_term_2)};
                        const __m128 vcos_term_2{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_2))};
                        const __m128 vcosinc_term_2{_mm_div_ps(vcos_term_2,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_2,vcosinc_denom_2)))};
                        const __m128 vlsrc_sample_2{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_2,vcosinc_term_2))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+8ull], vlsrc_sample_2);
                        const __m128 vt_i_3{_mm_setr_ps((jj+12.0f)*invLT,(jj+13.0f)*invLT,(jj+14.0f)*invLT,(jj+15.0f)*invLT)};
                        const __m128 varg_term_3{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_3)};
                        const __m128 vcosinc_denom_3{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_3))};
                        const __m128 vsin_term_3{_mm_sin_ps(varg_term_3)};
                        const __m128 vsinc_term_3{_mm_div_ps(vsin_term_3,varg_term_3)};
                        const __m128 vcos_term_3{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_3))};
                        const __m128 vcosinc_term_3{_mm_div_ps(vcos_term_3,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_3,vcosinc_denom_3)))};
                        const __m128 vlsrc_sample_3{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_3,vcosinc_term_3))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+12ull], vlsrc_sample_3);
                        const __m128 vt_i_4{_mm_setr_ps((jj+16.0f)*invLT,(jj+17.0f)*invLT,(jj+18.0f)*invLT,(jj+19.0f)*invLT)};
                        const __m128 varg_term_4{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_4)};
                        const __m128 vcosinc_denom_4{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_4))};
                        const __m128 vsin_term_4{_mm_sin_ps(varg_term_4)};
                        const __m128 vsinc_term_4{_mm_div_ps(vsin_term_4,varg_term_4)};
                        const __m128 vcos_term_4{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_4))};
                        const __m128 vcosinc_term_4{_mm_div_ps(vcos_term_4,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_4,vcosinc_denom_4)))};
                        const __m128 vlsrc_sample_4{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_4,vcosinc_term_4))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+16ull], vlsrc_sample_4);
                        const __m128 vt_i_5{_mm_setr_ps((jj+20.0f)*invLT,(jj+21.0f)*invLT,(jj+22.0f)*invLT,(jj+23.0f)*invLT)};
                        const __m128 varg_term_5{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_5)};
                        const __m128 vcosinc_denom_5{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_5))};
                        const __m128 vsin_term_5{_mm_sin_ps(varg_term_5)};
                        const __m128 vsinc_term_5{_mm_div_ps(vsin_term_5,varg_term_5)};
                        const __m128 vcos_term_5{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_5))};
                        const __m128 vcosinc_term_5{_mm_div_ps(vcos_term_5,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_5,vcosinc_denom_5)))};
                        const __m128 vlsrc_sample_5{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_5,vcosinc_term_5))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+20ull], vlsrc_sample_5);
                        const __m128 vt_i_6{_mm_setr_ps((jj+24.0f)*invLT,(jj+25.0f)*invLT,(jj+26.0f)*invLT,(jj+27.0f)*invLT)};
                        const __m128 varg_term_6{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_6)};
                        const __m128 vcosinc_denom_6{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_6))};
                        const __m128 vsin_term_6{_mm_sin_ps(varg_term_6)};
                        const __m128 vsinc_term_6{_mm_div_ps(vsin_term_6,varg_term_6)};
                        const __m128 vcos_term_6{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_6))};
                        const __m128 vcosinc_term_6{_mm_div_ps(vcos_term_6,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_6,vcosinc_denom_6)))};
                        const __m128 vlsrc_sample_6{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_6,vcosinc_term_6))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+24ull], vlsrc_sample_6);
                        const __m128 vt_i_7{_mm_setr_ps((jj+28.0f)*invLT,(jj+29.0f)*invLT,(jj+30.0f)*invLT,(jj+31.0f)*invLT)};
                        const __m128 varg_term_7{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_7)};
                        const __m128 vcosinc_denom_7{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_7))};
                        const __m128 vsin_term_7{_mm_sin_ps(varg_term_7)};
                        const __m128 vsinc_term_7{_mm_div_ps(vsin_term_7,varg_term_7)};
                        const __m128 vcos_term_7{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_7))};
                        const __m128 vcosinc_term_7{_mm_div_ps(vcos_term_7,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_7,vcosinc_denom_7)))};
                        const __m128 vlsrc_sample_7{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_7,vcosinc_term_7))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+28ull], vlsrc_sample_7);
                }

                for(; (i+23ull) < this->m_nTLsamples; i += 24ull,jj = 24.0f) 
                {
                        const __m128 vt_i_0{_mm_setr_ps(jj*invLT,(jj+1.0f)*invLT,(jj+2.0f)*invLT,(jj+3.0f)*invLT)};
                        const __m128 varg_term_0{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0)};
                        const __m128 vcosinc_denom_0{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_0))};
                        const __m128 vsin_term_0{_mm_sin_ps(varg_term_0)};
                        const __m128 vsinc_term_0{_mm_div_ps(vsin_term_0,varg_term_0)};
                        const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_0))};
                        const __m128 vcosinc_term_0{_mm_div_ps(vcos_term_0,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_0,vcosinc_denom_0)))};
                        const __m128 vlsrc_sample_0{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_0,vcosinc_term_0))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+0ull], vlsrc_sample_0);
                        const __m128 vt_i_1{_mm_setr_ps((jj+4.0f)*invLT,(jj+5.0f)*invLT,(jj+6.0f)*invLT,(jj+7.0f)*invLT)};
                        const __m128 varg_term_1{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_1)};
                        const __m128 vcosinc_denom_1{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_1))};
                        const __m128 vsin_term_1{_mm_sin_ps(varg_term_1)};
                        const __m128 vsinc_term_1{_mm_div_ps(vsin_term_1,varg_term_1)};
                        const __m128 vcos_term_1{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_1))};
                        const __m128 vcosinc_term_1{_mm_div_ps(vcos_term_1,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_1,vcosinc_denom_1)))};
                        const __m128 vlsrc_sample_1{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_1,vcosinc_term_1))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+4ull], vlsrc_sample_1);
                        const __m128 vt_i_2{_mm_setr_ps((jj+8.0f)*invLT,(jj+9.0f)*invLT,(jj+10.0f)*invLT,(jj+11.0f)*invLT)};
                        const __m128 varg_term_2{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_2)};
                        const __m128 vcosinc_denom_2{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_2))};
                        const __m128 vsin_term_2{_mm_sin_ps(varg_term_2)};
                        const __m128 vsinc_term_2{_mm_div_ps(vsin_term_2,varg_term_2)};
                        const __m128 vcos_term_2{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_2))};
                        const __m128 vcosinc_term_2{_mm_div_ps(vcos_term_2,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_2,vcosinc_denom_2)))};
                        const __m128 vlsrc_sample_2{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_2,vcosinc_term_2))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+8ull], vlsrc_sample_2);
                        const __m128 vt_i_3{_mm_setr_ps((jj+12.0f)*invLT,(jj+13.0f)*invLT,(jj+14.0f)*invLT,(jj+15.0f)*invLT)};
                        const __m128 varg_term_3{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_3)};
                        const __m128 vcosinc_denom_3{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_3))};
                        const __m128 vsin_term_3{_mm_sin_ps(varg_term_3)};
                        const __m128 vsinc_term_3{_mm_div_ps(vsin_term_3,varg_term_3)};
                        const __m128 vcos_term_3{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_3))};
                        const __m128 vcosinc_term_3{_mm_div_ps(vcos_term_3,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_3,vcosinc_denom_3)))};
                        const __m128 vlsrc_sample_3{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_3,vcosinc_term_3))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+12ull], vlsrc_sample_3);
                        const __m128 vt_i_4{_mm_setr_ps((jj+16.0f)*invLT,(jj+17.0f)*invLT,(jj+18.0f)*invLT,(jj+19.0f)*invLT)};
                        const __m128 varg_term_4{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_4)};
                        const __m128 vcosinc_denom_4{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_4))};
                        const __m128 vsin_term_4{_mm_sin_ps(varg_term_4)};
                        const __m128 vsinc_term_4{_mm_div_ps(vsin_term_4,varg_term_4)};
                        const __m128 vcos_term_4{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_4))};
                        const __m128 vcosinc_term_4{_mm_div_ps(vcos_term_4,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_4,vcosinc_denom_4)))};
                        const __m128 vlsrc_sample_4{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_4,vcosinc_term_4))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+16ull], vlsrc_sample_4);
                        const __m128 vt_i_5{_mm_setr_ps((jj+20.0f)*invLT,(jj+21.0f)*invLT,(jj+22.0f)*invLT,(jj+23.0f)*invLT)};
                        const __m128 varg_term_5{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_5)};
                        const __m128 vcosinc_denom_5{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_5))};
                        const __m128 vsin_term_5{_mm_sin_ps(varg_term_5)};
                        const __m128 vsinc_term_5{_mm_div_ps(vsin_term_5,varg_term_5)};
                        const __m128 vcos_term_5{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_5))};
                        const __m128 vcosinc_term_5{_mm_div_ps(vcos_term_5,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_5,vcosinc_denom_5)))};
                        const __m128 vlsrc_sample_5{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_5,vcosinc_term_5))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+20ull], vlsrc_sample_5);
                }

                for(; (i+15ull) < this->m_nTLsamples; i += 16ull,jj = 16.0f) 
                {
                        const __m128 vt_i_0{_mm_setr_ps(jj*invLT,(jj+1.0f)*invLT,(jj+2.0f)*invLT,(jj+3.0f)*invLT)};
                        const __m128 varg_term_0{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0)};
                        const __m128 vcosinc_denom_0{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_0))};
                        const __m128 vsin_term_0{_mm_sin_ps(varg_term_0)};
                        const __m128 vsinc_term_0{_mm_div_ps(vsin_term_0,varg_term_0)};
                        const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_0))};
                        const __m128 vcosinc_term_0{_mm_div_ps(vcos_term_0,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_0,vcosinc_denom_0)))};
                        const __m128 vlsrc_sample_0{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_0,vcosinc_term_0))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+0ull], vlsrc_sample_0);
                        const __m128 vt_i_1{_mm_setr_ps((jj+4.0f)*invLT,(jj+5.0f)*invLT,(jj+6.0f)*invLT,(jj+7.0f)*invLT)};
                        const __m128 varg_term_1{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_1)};
                        const __m128 vcosinc_denom_1{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_1))};
                        const __m128 vsin_term_1{_mm_sin_ps(varg_term_1)};
                        const __m128 vsinc_term_1{_mm_div_ps(vsin_term_1,varg_term_1)};
                        const __m128 vcos_term_1{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_1))};
                        const __m128 vcosinc_term_1{_mm_div_ps(vcos_term_1,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_1,vcosinc_denom_1)))};
                        const __m128 vlsrc_sample_1{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_1,vcosinc_term_1))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+4ull], vlsrc_sample_1);
                        const __m128 vt_i_2{_mm_setr_ps((jj+8.0f)*invLT,(jj+9.0f)*invLT,(jj+10.0f)*invLT,(jj+11.0f)*invLT)};
                        const __m128 varg_term_2{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_2)};
                        const __m128 vcosinc_denom_2{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_2))};
                        const __m128 vsin_term_2{_mm_sin_ps(varg_term_2)};
                        const __m128 vsinc_term_2{_mm_div_ps(vsin_term_2,varg_term_2)};
                        const __m128 vcos_term_2{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_2))};
                        const __m128 vcosinc_term_2{_mm_div_ps(vcos_term_2,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_2,vcosinc_denom_2)))};
                        const __m128 vlsrc_sample_2{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_2,vcosinc_term_2))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+8ull], vlsrc_sample_2);
                        const __m128 vt_i_3{_mm_setr_ps((jj+12.0f)*invLT,(jj+13.0f)*invLT,(jj+14.0f)*invLT,(jj+15.0f)*invLT)};
                        const __m128 varg_term_3{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_3)};
                        const __m128 vcosinc_denom_3{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_3))};
                        const __m128 vsin_term_3{_mm_sin_ps(varg_term_3)};
                        const __m128 vsinc_term_3{_mm_div_ps(vsin_term_3,varg_term_3)};
                        const __m128 vcos_term_3{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_3))};
                        const __m128 vcosinc_term_3{_mm_div_ps(vcos_term_3,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_3,vcosinc_denom_3)))};
                        const __m128 vlsrc_sample_3{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_3,vcosinc_term_3))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+12ull], vlsrc_sample_3);
                }

                for(; (i+7ull) < this->m_nTLsamples; i += 8ull,jj = 8.0f) 
                {
                        const __m128 vt_i_0{_mm_setr_ps(jj*invLT,(jj+1.0f)*invLT,(jj+2.0f)*invLT,(jj+3.0f)*invLT)};
                        const __m128 varg_term_0{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0)};
                        const __m128 vcosinc_denom_0{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_0))};
                        const __m128 vsin_term_0{_mm_sin_ps(varg_term_0)};
                        const __m128 vsinc_term_0{_mm_div_ps(vsin_term_0,varg_term_0)};
                        const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_0))};
                        const __m128 vcosinc_term_0{_mm_div_ps(vcos_term_0,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_0,vcosinc_denom_0)))};
                        const __m128 vlsrc_sample_0{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_0,vcosinc_term_0))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+0ull], vlsrc_sample_0);
                        const __m128 vt_i_1{_mm_setr_ps((jj+4.0f)*invLT,(jj+5.0f)*invLT,(jj+6.0f)*invLT,(jj+7.0f)*invLT)};
                        const __m128 varg_term_1{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_1)};
                        const __m128 vcosinc_denom_1{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_1))};
                        const __m128 vsin_term_1{_mm_sin_ps(varg_term_1)};
                        const __m128 vsinc_term_1{_mm_div_ps(vsin_term_1,varg_term_1)};
                        const __m128 vcos_term_1{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_1))};
                        const __m128 vcosinc_term_1{_mm_div_ps(vcos_term_1,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_1,vcosinc_denom_1)))};
                        const __m128 vlsrc_sample_1{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_1,vcosinc_term_1))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+4ull], vlsrc_sample_1);
                }

                for(; (i+3ull) < this->m_nTLsamples; i += 4ull,jj = 4.0f) 
                {
                        const __m128 vt_i_0{_mm_setr_ps(jj*invLT,(jj+1.0f)*invLT,(jj+2.0f)*invLT,(jj+3.0f)*invLT)};
                        const __m128 varg_term_0{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0)};
                        const __m128 vcosinc_denom_0{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_0))};
                        const __m128 vsin_term_0{_mm_sin_ps(varg_term_0)};
                        const __m128 vsinc_term_0{_mm_div_ps(vsin_term_0,varg_term_0)};
                        const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_0))};
                        const __m128 vcosinc_term_0{_mm_div_ps(vcos_term_0,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_0,vcosinc_denom_0)))};
                        const __m128 vlsrc_sample_0{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_0,vcosinc_term_0))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+0ull], vlsrc_sample_0);
                }

                 for(; (i+0ull) < this->m_nTLsamples; i += 1ull,jj += 1.0f)  
                 {
                                const float t_i{jj*invLT};
                                const float arg_term_0{C6283185307179586476925286766559*t_i};
                                const float cosinc_denom_0{4.0f*this->m_beta*t_i};
#if (CPM_PULSE_SHAPERS_USE_CEPHES)
                                const float sin_term_0{ceph_sinf(arg_term_0)};
                                const float cos_term_0{ceph_cosf(this->m_beta*arg_term_0)};
#else 
                                const float sin_term_0{std::sin(arg_term_0)};
                                const float cos_term_0{std::cos(this->m_beta*arg_term_0)};
#endif
                                const float sinc_term_0{sin_term_0/arg_term_0};
                                const float cosinc_term_0{cos_term_0/(1.0f-(cosinc_denom_0*cosinc_denom_0))};
                                const float lsrc_sample_0{invLT*sinc_term_0*cosinc_term_0};
                                this->m_lsrc_pulse.m_data[i] = lsrc_sample_0;
                }

        }
        else 
        {
                for(i = 0ull;(i+31ull) < this->m_nTLsamples;i += 32ull) 
                {
                        _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
                        const __m128 vt_i_0{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinvLT)};
                        const __m128 varg_term_0{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0)};
                        const __m128 vcosinc_denom_0{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_0))};
                        const __m128 vsin_term_0{_mm_sin_ps(varg_term_0)};
                        const __m128 vsinc_term_0{_mm_div_ps(vsin_term_0,varg_term_0)};
                        const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_0))};
                        const __m128 vcosinc_term_0{_mm_div_ps(vcos_term_0,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_0,vcosinc_denom_0)))};
                        const __m128 vlsrc_sample_0{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_0,vcosinc_term_0))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+0ull], vlsrc_sample_0);
                        const __m128 vt_i_1{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull]),vinvLT)};
                        const __m128 varg_term_1{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_1)};
                        const __m128 vcosinc_denom_1{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_1))};
                        const __m128 vsin_term_1{_mm_sin_ps(varg_term_1)};
                        const __m128 vsinc_term_1{_mm_div_ps(vsin_term_1,varg_term_1)};
                        const __m128 vcos_term_1{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_1))};
                        const __m128 vcosinc_term_1{_mm_div_ps(vcos_term_1,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_1,vcosinc_denom_1)))};
                        const __m128 vlsrc_sample_1{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_1,vcosinc_term_1))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+4ull], vlsrc_sample_1);
                        const __m128 vt_i_2{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]),vinvLT)};
                        const __m128 varg_term_2{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_2)};
                        const __m128 vcosinc_denom_2{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_2))};
                        const __m128 vsin_term_2{_mm_sin_ps(varg_term_2)};
                        const __m128 vsinc_term_2{_mm_div_ps(vsin_term_2,varg_term_2)};
                        const __m128 vcos_term_2{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_2))};
                        const __m128 vcosinc_term_2{_mm_div_ps(vcos_term_2,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_2,vcosinc_denom_2)))};
                        const __m128 vlsrc_sample_2{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_2,vcosinc_term_2))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+8ull], vlsrc_sample_2);
                        const __m128 vt_i_3{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+12ull]),vinvLT)};
                        const __m128 varg_term_3{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_3)};
                        const __m128 vcosinc_denom_3{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_3))};
                        const __m128 vsin_term_3{_mm_sin_ps(varg_term_3)};
                        const __m128 vsinc_term_3{_mm_div_ps(vsin_term_3,varg_term_3)};
                        const __m128 vcos_term_3{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_3))};
                        const __m128 vcosinc_term_3{_mm_div_ps(vcos_term_3,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_3,vcosinc_denom_3)))};
                        const __m128 vlsrc_sample_3{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_3,vcosinc_term_3))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+12ull], vlsrc_sample_3);
                        _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+16ull],_MM_HINT_T0);
                        const __m128 vt_i_4{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+16ull]),vinvLT)};
                        const __m128 varg_term_4{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_4)};
                        const __m128 vcosinc_denom_4{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_4))};
                        const __m128 vsin_term_4{_mm_sin_ps(varg_term_4)};
                        const __m128 vsinc_term_4{_mm_div_ps(vsin_term_4,varg_term_4)};
                        const __m128 vcos_term_4{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_4))};
                        const __m128 vcosinc_term_4{_mm_div_ps(vcos_term_4,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_4,vcosinc_denom_4)))};
                        const __m128 vlsrc_sample_4{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_4,vcosinc_term_4))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+16ull], vlsrc_sample_4);
                        const __m128 vt_i_5{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+20ull]),vinvLT)};
                        const __m128 varg_term_5{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_5)};
                        const __m128 vcosinc_denom_5{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_5))};
                        const __m128 vsin_term_5{_mm_sin_ps(varg_term_5)};
                        const __m128 vsinc_term_5{_mm_div_ps(vsin_term_5,varg_term_5)};
                        const __m128 vcos_term_5{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_5))};
                        const __m128 vcosinc_term_5{_mm_div_ps(vcos_term_5,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_5,vcosinc_denom_5)))};
                        const __m128 vlsrc_sample_5{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_5,vcosinc_term_5))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+20ull], vlsrc_sample_5);
                        const __m128 vt_i_6{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+24ull]),vinvLT)};
                        const __m128 varg_term_6{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_6)};
                        const __m128 vcosinc_denom_6{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_6))};
                        const __m128 vsin_term_6{_mm_sin_ps(varg_term_6)};
                        const __m128 vsinc_term_6{_mm_div_ps(vsin_term_6,varg_term_6)};
                        const __m128 vcos_term_6{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_6))};
                        const __m128 vcosinc_term_6{_mm_div_ps(vcos_term_6,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_6,vcosinc_denom_6)))};
                        const __m128 vlsrc_sample_6{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_6,vcosinc_term_6))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+24ull], vlsrc_sample_6);
                        const __m128 vt_i_7{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+28ull]),vinvLT)};
                        const __m128 varg_term_7{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_7)};
                        const __m128 vcosinc_denom_7{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_7))};
                        const __m128 vsin_term_7{_mm_sin_ps(varg_term_7)};
                        const __m128 vsinc_term_7{_mm_div_ps(vsin_term_7,varg_term_7)};
                        const __m128 vcos_term_7{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_7))};
                        const __m128 vcosinc_term_7{_mm_div_ps(vcos_term_7,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_7,vcosinc_denom_7)))};
                        const __m128 vlsrc_sample_7{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_7,vcosinc_term_7))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+28ull], vlsrc_sample_7);
                }

                for(; (i+23ull) < this->m_nTLsamples; i += 24ull) 
                {
                        _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
                        const __m128 vt_i_0{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinvLT)};
                        const __m128 varg_term_0{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0)};
                        const __m128 vcosinc_denom_0{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_0))};
                        const __m128 vsin_term_0{_mm_sin_ps(varg_term_0)};
                        const __m128 vsinc_term_0{_mm_div_ps(vsin_term_0,varg_term_0)};
                        const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_0))};
                        const __m128 vcosinc_term_0{_mm_div_ps(vcos_term_0,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_0,vcosinc_denom_0)))};
                        const __m128 vlsrc_sample_0{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_0,vcosinc_term_0))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+0ull], vlsrc_sample_0);
                        const __m128 vt_i_1{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull]),vinvLT)};
                        const __m128 varg_term_1{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_1)};
                        const __m128 vcosinc_denom_1{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_1))};
                        const __m128 vsin_term_1{_mm_sin_ps(varg_term_1)};
                        const __m128 vsinc_term_1{_mm_div_ps(vsin_term_1,varg_term_1)};
                        const __m128 vcos_term_1{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_1))};
                        const __m128 vcosinc_term_1{_mm_div_ps(vcos_term_1,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_1,vcosinc_denom_1)))};
                        const __m128 vlsrc_sample_1{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_1,vcosinc_term_1))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+4ull], vlsrc_sample_1);
                        const __m128 vt_i_2{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]),vinvLT)};
                        const __m128 varg_term_2{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_2)};
                        const __m128 vcosinc_denom_2{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_2))};
                        const __m128 vsin_term_2{_mm_sin_ps(varg_term_2)};
                        const __m128 vsinc_term_2{_mm_div_ps(vsin_term_2,varg_term_2)};
                        const __m128 vcos_term_2{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_2))};
                        const __m128 vcosinc_term_2{_mm_div_ps(vcos_term_2,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_2,vcosinc_denom_2)))};
                        const __m128 vlsrc_sample_2{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_2,vcosinc_term_2))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+8ull], vlsrc_sample_2);
                        const __m128 vt_i_3{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+12ull]),vinvLT)};
                        const __m128 varg_term_3{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_3)};
                        const __m128 vcosinc_denom_3{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_3))};
                        const __m128 vsin_term_3{_mm_sin_ps(varg_term_3)};
                        const __m128 vsinc_term_3{_mm_div_ps(vsin_term_3,varg_term_3)};
                        const __m128 vcos_term_3{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_3))};
                        const __m128 vcosinc_term_3{_mm_div_ps(vcos_term_3,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_3,vcosinc_denom_3)))};
                        const __m128 vlsrc_sample_3{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_3,vcosinc_term_3))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+12ull], vlsrc_sample_3);
                        _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+16ull],_MM_HINT_T0);
                        const __m128 vt_i_4{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+16ull]),vinvLT)};
                        const __m128 varg_term_4{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_4)};
                        const __m128 vcosinc_denom_4{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_4))};
                        const __m128 vsin_term_4{_mm_sin_ps(varg_term_4)};
                        const __m128 vsinc_term_4{_mm_div_ps(vsin_term_4,varg_term_4)};
                        const __m128 vcos_term_4{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_4))};
                        const __m128 vcosinc_term_4{_mm_div_ps(vcos_term_4,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_4,vcosinc_denom_4)))};
                        const __m128 vlsrc_sample_4{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_4,vcosinc_term_4))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+16ull], vlsrc_sample_4);
                        const __m128 vt_i_5{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+20ull]),vinvLT)};
                        const __m128 varg_term_5{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_5)};
                        const __m128 vcosinc_denom_5{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_5))};
                        const __m128 vsin_term_5{_mm_sin_ps(varg_term_5)};
                        const __m128 vsinc_term_5{_mm_div_ps(vsin_term_5,varg_term_5)};
                        const __m128 vcos_term_5{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_5))};
                        const __m128 vcosinc_term_5{_mm_div_ps(vcos_term_5,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_5,vcosinc_denom_5)))};
                        const __m128 vlsrc_sample_5{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_5,vcosinc_term_5))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+20ull], vlsrc_sample_5);
                
                }

                for(; (i+15ull) < this->m_nTLsamples; i += 16ull) 
                {
                        _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
                        const __m128 vt_i_0{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinvLT)};
                        const __m128 varg_term_0{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0)};
                        const __m128 vcosinc_denom_0{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_0))};
                        const __m128 vsin_term_0{_mm_sin_ps(varg_term_0)};
                        const __m128 vsinc_term_0{_mm_div_ps(vsin_term_0,varg_term_0)};
                        const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_0))};
                        const __m128 vcosinc_term_0{_mm_div_ps(vcos_term_0,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_0,vcosinc_denom_0)))};
                        const __m128 vlsrc_sample_0{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_0,vcosinc_term_0))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+0ull], vlsrc_sample_0);
                        const __m128 vt_i_1{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull]),vinvLT)};
                        const __m128 varg_term_1{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_1)};
                        const __m128 vcosinc_denom_1{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_1))};
                        const __m128 vsin_term_1{_mm_sin_ps(varg_term_1)};
                        const __m128 vsinc_term_1{_mm_div_ps(vsin_term_1,varg_term_1)};
                        const __m128 vcos_term_1{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_1))};
                        const __m128 vcosinc_term_1{_mm_div_ps(vcos_term_1,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_1,vcosinc_denom_1)))};
                        const __m128 vlsrc_sample_1{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_1,vcosinc_term_1))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+4ull], vlsrc_sample_1);
                        const __m128 vt_i_2{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]),vinvLT)};
                        const __m128 varg_term_2{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_2)};
                        const __m128 vcosinc_denom_2{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_2))};
                        const __m128 vsin_term_2{_mm_sin_ps(varg_term_2)};
                        const __m128 vsinc_term_2{_mm_div_ps(vsin_term_2,varg_term_2)};
                        const __m128 vcos_term_2{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_2))};
                        const __m128 vcosinc_term_2{_mm_div_ps(vcos_term_2,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_2,vcosinc_denom_2)))};
                        const __m128 vlsrc_sample_2{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_2,vcosinc_term_2))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+8ull], vlsrc_sample_2);
                        const __m128 vt_i_3{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+12ull]),vinvLT)};
                        const __m128 varg_term_3{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_3)};
                        const __m128 vcosinc_denom_3{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_3))};
                        const __m128 vsin_term_3{_mm_sin_ps(varg_term_3)};
                        const __m128 vsinc_term_3{_mm_div_ps(vsin_term_3,varg_term_3)};
                        const __m128 vcos_term_3{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_3))};
                        const __m128 vcosinc_term_3{_mm_div_ps(vcos_term_3,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_3,vcosinc_denom_3)))};
                        const __m128 vlsrc_sample_3{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_3,vcosinc_term_3))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+12ull], vlsrc_sample_3);
                                       
                }

                for(; (i+7ull) < this->m_nTLsamples; i += 8ull) 
                {
                        const __m128 vt_i_0{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinvLT)};
                        const __m128 varg_term_0{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0)};
                        const __m128 vcosinc_denom_0{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_0))};
                        const __m128 vsin_term_0{_mm_sin_ps(varg_term_0)};
                        const __m128 vsinc_term_0{_mm_div_ps(vsin_term_0,varg_term_0)};
                        const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_0))};
                        const __m128 vcosinc_term_0{_mm_div_ps(vcos_term_0,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_0,vcosinc_denom_0)))};
                        const __m128 vlsrc_sample_0{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_0,vcosinc_term_0))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+0ull], vlsrc_sample_0);
                        const __m128 vt_i_1{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull]),vinvLT)};
                        const __m128 varg_term_1{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_1)};
                        const __m128 vcosinc_denom_1{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_1))};
                        const __m128 vsin_term_1{_mm_sin_ps(varg_term_1)};
                        const __m128 vsinc_term_1{_mm_div_ps(vsin_term_1,varg_term_1)};
                        const __m128 vcos_term_1{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_1))};
                        const __m128 vcosinc_term_1{_mm_div_ps(vcos_term_1,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_1,vcosinc_denom_1)))};
                        const __m128 vlsrc_sample_1{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_1,vcosinc_term_1))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+4ull], vlsrc_sample_1);
                }

                for(; (i+3ull) < this->m_nTLsamples; i += 4ull) 
                {
                        const __m128 vt_i_0{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinvLT)};
                        const __m128 varg_term_0{_mm_mul_ps(vC6283185307179586476925286766559,vt_i_0)};
                        const __m128 vcosinc_denom_0{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i_0))};
                        const __m128 vsin_term_0{_mm_sin_ps(varg_term_0)};
                        const __m128 vsinc_term_0{_mm_div_ps(vsin_term_0,varg_term_0)};
                        const __m128 vcos_term_0{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term_0))};
                        const __m128 vcosinc_term_0{_mm_div_ps(vcos_term_0,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom_0,vcosinc_denom_0)))};
                        const __m128 vlsrc_sample_0{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term_0,vcosinc_term_0))};
                        _mm_store_ps(&this->m_lsrc_pulse.m_data[i+0ull], vlsrc_sample_0);
                }

                 for(; (i+0ull) < this->m_nTLsamples; i += 1ull)  
                 {
                                const float t_i{gms::math::LUT_loop_indices_2257_align16[i]*invLT};
                                const float arg_term_0{C6283185307179586476925286766559*t_i};
                                const float cosinc_denom_0{4.0f*this->m_beta*t_i};
#if (CPM_PULSE_SHAPERS_USE_CEPHES)
                                const float sin_term_0{ceph_sinf(arg_term_0)};
                                const float cos_term_0{ceph_cosf(this->m_beta*arg_term_0)};
#else 
                                const float sin_term_0{std::sin(arg_term_0)};
                                const float cos_term_0{std::cos(this->m_beta*arg_term_0)};
#endif
                                const float sinc_term_0{sin_term_0/arg_term_0};
                                const float cosinc_term_0{cos_term_0/(1.0f-(cosinc_denom_0*cosinc_denom_0))};
                                const float lsrc_sample_0{invLT*sinc_term_0*cosinc_term_0};
                                this->m_lsrc_pulse.m_data[i] = lsrc_sample_0;
                }

        }
        return (0);
}

std::int32_t 
gms::radiolocation
::cpm_pulse_shapers_t
::generate_lrc_pulse_avx_u8x(const bool do_consts_prefetch)
{
        using namespace gms::math;
        constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
        const float twoLT{2.0f*this->m_L*this->m_T};
        const float inv2LT{1.0f/twoLT};
        const float invLT{1.0f/(this->m_L*this->m_T)};
        const __m256 vinv2LT{_mm256_set1_ps(inv2LT)};
        const __m256 vinvLT{_mm256_set1_ps(invLT)};
        __m256 c0_7,c8_15,c16_23,c24_31;
        __m256 c32_39,c40_47,c48_55,c56_63;
        __m256 vC6283185307179586476925286766559;
        __m256 vone;
        __m256 vt_i;
        __m256 vcos_term;
        __m256 vlrc_value;
        if(!do_consts_prefetch)
        {
                c0_7   = _mm256_setr_ps(0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f);
                c8_15  = _mm256_setr_ps(8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f);
                c16_23 = _mm256_setr_ps(16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f);
                c24_31 = _mm256_setr_ps(24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f);
                c32_39 = _mm256_setr_ps(32.0f,33.0f,34.0f,35.0f,36.0f,37.0f,38.0f,39.0f);
                c40_47 = _mm256_setr_ps(40.0f,41.0f,42.0f,43.0f,44.0f,45.0f,46.0f,47.0f);
                c48_55 = _mm256_setr_ps(48.0f,49.0f,50.0f,51.0f,52.0f,53.0f,54.0f,55.0f);
                c56_63 = _mm256_setr_ps(56.0f,57.0f,58.0f,59.0f,60.0f,61.0f,62.0f,63.0f);
                vC6283185307179586476925286766559 = _mm256_set1_ps(C6283185307179586476925286766559);
                vone                              = _mm256_set1_ps(1.0f);
        }
        else 
        {
                __ATTR_ALIGN__(32) const float prefetched_constants[80] = 
                {
                        0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,
                        8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f,
                        16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f,
                        24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f,
                        32.0f,33.0f,34.0f,35.0f,36.0f,37.0f,38.0f,39.0f,
                        40.0f,41.0f,42.0f,43.0f,44.0f,45.0f,46.0f,47.0f,
                        48.0f,49.0f,50.0f,51.0f,52.0f,53.0f,54.0f,55.0f,
                        56.0f,57.0f,58.0f,59.0f,60.0f,61.0f,62.0f,63.0f,
                        6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                        6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                        6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                        6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                        +1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f                       
                };
                _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
                _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
                _mm_prefetch((const char*)&prefetched_constants[32],_MM_HINT_T0);
                _mm_prefetch((const char*)&prefetched_constants[48],_MM_HINT_T0);
                _mm_prefetch((const char*)&prefetched_constants[64],_MM_HINT_T0);
                c0_7   = _mm256_load_ps(&prefetched_constants[0]);
                c8_15  = _mm256_load_ps(&prefetched_constants[8]);
                c16_23 = _mm256_load_ps(&prefetched_constants[16]);
                c24_31 = _mm256_load_ps(&prefetched_constants[24]);
                c32_39 = _mm256_load_ps(&prefetched_constants[32]);
                c40_47 = _mm256_load_ps(&prefetched_constants[40]);
                c48_55 = _mm256_load_ps(&prefetched_constants[48]);
                c56_63 = _mm256_load_ps(&prefetched_constants[56]);
                vC6283185307179586476925286766559 = _mm256_load_ps(&prefetched_constants[64]);
                vone  = _mm256_load_ps(&prefetched_constants[72]);
        }
        std::size_t i;
        float j;
        constexpr std::size_t LUT_loop_idx_threshold{2257ull};
        if(__builtin_expect(this->m_nTLsamples>LUT_loop_idx_threshold,0))
        {
             for(i = 0ull,j = 0.0f;(i+63ull) < this->m_nTLsamples;i += 64ull,j += 64.0f) 
             {
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c0_7),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c8_15),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+8ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c16_23),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+16ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c24_31),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+24ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c32_39),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+32ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c40_47),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+40ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c48_55),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+48ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c56_63),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+56ull],vlrc_value);
             }

             for(; (i+47ull) < this->m_nTLsamples;i += 48ull,j += 48.0f) 
             {
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c0_7),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c8_15),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+8ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c16_23),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+16ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c24_31),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+24ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c32_39),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+32ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c40_47),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+40ull],vlrc_value);
             }

             for(; (i+31ull) < this->m_nTLsamples;i += 32ull,j += 32.0f) 
             {
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c0_7),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c8_15),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+8ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c16_23),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+16ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c24_31),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+24ull],vlrc_value);
              }

              for(; (i+15ull) < this->m_nTLsamples;i += 16ull,j += 16.0f) 
              {
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c0_7),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c8_15),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+8ull],vlrc_value);
              }

              for(; (i+7ull) < this->m_nTLsamples;i += 8ull,j += 8.0f) 
              {
                 vt_i       = _mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c0_7),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value);
              }

              for(; (i+0ull) < this->m_nTLsamples; i += 1ull)   
              {
                   const float t_i{static_cast<float>(i)};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                   const float cos_term_0{ceph_cosf(C6283185307179586476925286766559*t_i*invLT)};
#else
                   const float cos_term_0{std::cos(C6283185307179586476925286766559*t_i*invLT)};
#endif
                   const float lrc_val_0{inv2LT*(1.0f-cos_term_0)};
                   this->m_lrc_pulse.m_data[i] = lrc_val_0;
              }
        }
        else 
        {
             for(i = 0ull;(i+63ull) < this->m_nTLsamples;i += 64ull) 
             {
                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align32[i+0ull],_MM_HINT_T0);
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+8ull],vlrc_value);
                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align32[i+16ull],_MM_HINT_T0);
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+16ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+16ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+24ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+24ull],vlrc_value);
                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align32[i+32ull],_MM_HINT_T0);
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+32ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+32ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+40ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+40ull],vlrc_value);
                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align32[i+48ull],_MM_HINT_T0);
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+48ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+48ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+56ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+56ull],vlrc_value);
             }

             for(; (i+47ull) < this->m_nTLsamples;i += 48ull) 
             {
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+8ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+16ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+16ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+24ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+24ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+32ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+32ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+40ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+40ull],vlrc_value);
             }

             for(; (i+31ull) < this->m_nTLsamples;i += 32ull) 
             {
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+8ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+16ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+16ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+24ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+24ull],vlrc_value);
             }

              for(; (i+15ull) < this->m_nTLsamples;i += 16ull) 
              {
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value);
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+8ull],vlrc_value);
              }

              for(; (i+7ull) < this->m_nTLsamples;i += 8ull) 
              {
                 vt_i       = _mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinvLT);
                 vcos_term  = _mm256_cos_ps(_mm256_mul_ps(vC6283185307179586476925286766559,vt_i));
                 vlrc_value = _mm256_mul_ps(vinv2LT,_mm256_sub_ps(vone,vcos_term));
                 _mm256_store_ps(&this->m_lrc_pulse.m_data[i+0ull],vlrc_value);
              }

              for(; (i+0ull) < this->m_nTLsamples; i += 1ull)   
              {
                   const float t_i{gms::math::LUT_loop_indices_2257_align16[i]};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                   const float cos_term_0{ceph_cosf(C6283185307179586476925286766559*t_i*invLT)};
#else
                   const float cos_term_0{std::cos(C6283185307179586476925286766559*t_i*invLT)};
#endif
                   const float lrc_val_0{inv2LT*(1.0f-cos_term_0)};
                   this->m_lrc_pulse.m_data[i] = lrc_val_0;
              }
        }
        return (0);
}

auto 
gms::radiolocation 
::operator<<(std::ostream &os,
             cpm_pulse_shapers_t &rhs)->std::ostream &
{
        std::cout << typeid(rhs).name() << "Begin: object state dump." << std::endl;
        std::cout << "m_nTLsamples    : "      << rhs.m_nTLsamples  << std::endl;
        std::cout << "m_nTsamples     : "      << rhs.m_nTsamples   << std::endl;   
        std::cout << "m_T             : "      << std::fixed << std::setprecision(7) << rhs.m_T << std::endl;
        std::cout << "m_L             : "      << std::fixed << std::setprecision(7) << rhs.m_L << std::endl;
        std::cout << "m_beta          : "      << std::fixed << std::setprecision(7) << rhs.m_beta << std::endl;
        std::cout << "m_BbT           : "      << std::fixed << std::setprecision(7) << rhs.m_BbT << std::endl;
        std::cout << "L-Rectangular Pulse -- Samples:" << "\n";
        for(std::size_t i{0ull}; i != rhs.m_nTLsamples; ++i)
        {
           os << std::fixed << std::setprecision(7) << rhs.m_lrec_pulse.m_data[i] << std::endl;
        }
        std::cout << "\n";
        std::cout << "L-Raised Cosine Pulse -- Samples:" << "\n";
        for(std::size_t i{0ull}; i != rhs.m_nTLsamples; ++i)
        {
           os << std::fixed << std::setprecision(7) << rhs.m_lrc_pulse.m_data[i] << std::endl;
        }
        std::cout << "\n";
        std::cout << "L-Spectrally Raised Cosine Pulse -- Samples:" << "\n";
        for(std::size_t i{0ull}; i != rhs.m_nTLsamples; ++i)
        {
           os << std::fixed << std::setprecision(7) << rhs.m_lsrc_pulse.m_data[i] << std::endl;
        }
        std::cout << "\n";
        std::cout << "Gaussian(MSK) Pulse -- Samples:" << "\n";
        for(std::size_t i{0ull}; i != rhs.m_nTsamples; ++i)
        {
           os << std::fixed << std::setprecision(7) << rhs.m_gmsk_pulse.m_data[i] << std::endl;
        }
        std::cout << "\n";
        std::cout << typeid(rhs).name() << "End: object state dump." << std::endl;
        return (os);
}



