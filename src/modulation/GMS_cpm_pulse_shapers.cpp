
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
m_tfm_pulse{darray_r4_t()},
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
m_tfm_pulse{darray_r4_t(this->m_nTsamples)},
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
m_tfm_pulse{other.m_tfm_pulse},
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
m_tfm_pulse{std::move(other.m_tfm_pulse)},
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
    this->m_tfm_pulse.operator=(other.m_tfm_pulse);
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
    this->m_tfm_pulse.operator=(std::move(other.m_tfm_pulse));
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
    if(__builtin_expect(this->m_nTLsamples>LUT_loop_indices_2257_align16,0)) 
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
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+20ull],vlrc_value_5) 
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
                 _mm_store_ps(&this->m_lrc_pulse.m_data[i+12ull],vlrc_value_3)
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
                this->m_lrc_pulse.m_data[j] = lrc_val_0;
          }
    }
    return (0);
}



