
#include <fstream>
#include <iomanip>
#include <immintrin.h>
#include "GMS_cpfsk_signal.h"
#include "GMS_loop_indices_LUT.h"

gms::radiolocation 
::cpfsk_signal_t
::cpfsk_signal_t() noexcept(true)
:
m_T{},
m_nfreqs{},
m_k{},
m_A{},
m_h{},
m_fc{},
m_ph0{},
m_Phik{},
m_cpfsk_signal{darray_r4_t()},
m_cpfsk_psd{darray_r4_t()}
{

}

gms::radiolocation
::cpfsk_signal_t
::cpfsk_signal_t(const std::size_t T,
                 const std::size_t nfreqs,
                 const std::uint32_t k,
                 const float A,
                 const float h,
                 const float fc,
                 const float ph0) noexcept(false)
:
m_T(T),
m_nfreqs{nfreqs},
m_k(k),
m_A(A),
m_h(h),
m_fc(fc),
m_ph0(ph0),
m_Phik{},
m_cpfsk_signal(darray_r4_t(this->m_T)),
m_cpfsk_psd(darray_r4_t(this->m_nfreqs))
{

}

gms::radiolocation 
::cpfsk_signal_t
::cpfsk_signal_t(const cpfsk_signal_t &rhs) noexcept(false)
:
m_T(rhs.m_T),
m_nfreqs(rhs.m_nfreqs)
m_k(rhs.m_k),
m_A(rhs.m_A),
m_h(rhs.m_h),
m_fc(rhs.m_fc),
m_ph0(rhs.m_ph0),
m_Phik(rhs.m_Phik),
m_cpfsk_signal(rhs.m_cpfsk_signal),
m_cpfsk_psd(rhs.m_cpfsk_psd)
{

}

gms::radiolocation 
::cpfsk_signal_t
::cpfsk_signal_t(cpfsk_signal_t &&rhs) noexcept(true) 
:
m_T(rhs.m_T),
m_nfreqs(rhs.m_nfreqs),
m_k(rhs.m_k),
m_A(rhs.m_A),
m_h(rhs.m_h),
m_fc(rhs.m_fc),
m_ph0(rhs.m_ph0),
m_Phik(rhs.m_Phik),
m_cpfsk_signal(std::move(rhs.m_cpfsk_signal)) 
m_cpfsk_psd(std::move(rhs.m_cpfsk_psd))
{

}

gms::radiolocation 
::cpfsk_signal_t
::~cpfsk_signal_t() noexcept(true)
{

}

gms::radiolocation
::cpfsk_signal_t &
gms::radiolocation
::cpfsk_signal_t
::operator=(const cpfsk_signal_t &other) noexcept(false) 
{
     if(__builtin_expect(this==&other,0)) { return (*this);}
     this->m_T      = other.m_T;
     this->m_nfreqs = other.m_nfreqs
     this->m_k      = other.m_k;
     this->m_A      = other.m_A;
     this->m_h      = other.m_h;
     this->m_fc     = other.m_fc;
     this->m_ph0    = other.m_ph0;
     this->m_Phik   =other.m_Phik;
     this->m_cpfsk_signal.operator=(other.m_cpfsk_signal);
     this->m_cpfsk_psd.operator=(other.m_cpfsk_psd);
     return (*this);
}

gms::radiolocation
::cpfsk_signal_t &
gms::radiolocation
::cpfsk_signal_t
::operator=(cpfsk_signal_t &&other) noexcept(true) 
{
     if(__builtin_expect(this==&other,0)) { return (*this);}
     this->m_T  = other.m_T;
     this->m_nfreqs = other.m_nfreqs;
     this->m_k  = other.m_k;
     this->m_A  = other.m_A;
     this->m_h  = other.m_h;
     this->m_fc = other.m_fc;
     this->m_ph0= other.m_ph0;
     this->m_Phik=other.m_Phik;
     this->m_cpfsk_signal.operator=(std::move(other.m_cpfsk_signal));
     this->m_cpfsk_psd.operator=(std::move(other.m_cpfsk_psd));
     return (*this);
}

void 
gms::radiolocation 
::cpfsk_signal_t
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
::cpfsk_signal_t
::compute_cpfsk_psd_analytic(const darray_r4_t &freqs,
                             const std::int32_t M) noexcept(true) 
{
     if(__builtin_expect(this->m_cpfsk_psd.mnx!=freqs.mnx,0)) { return (-1);}
     if(__builtin_expect(M<=1,0)) {return (-2);}
     constexpr float C628318530717958647692528676656{6.28318530717958647692528676656f};
     constexpr float C314159265358979323846264338328{3.14159265358979323846264338328f};
     const half_h{0.5f*this->m_h};
     const float sqr_A{this->m_A*this->m_A};
     const float f_M{static_cast<float>(M)};
     const float f_T{static_cast<float>(this->m_T)};
     const float invf_M{1.0f/f_M};
     const float lead_term{sqr_A*f_T*invf_M};
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
     const float sin_num{ceph_sin(f_M*C314159265358979323846264338328*this->m_h)};
     const float sin_den{ceph_sin(C314159265358979323846264338328*this->m_h)};
     const float Ca{f_M*(sin_num/sin_den)};
#else 
     const float sin_num{std::sin(f_M*C314159265358979323846264338328*this->m_h)};
     const float sin_den{std::sin(C314159265358979323846264338328*this->m_h)};
     const float Ca{f_M*(sin_num/sin_den)};
#endif
     float j_loop_term{0.0f};
     float i_loop_term{0.0f};
     float i_sum{}, gamma_i{},gamma_j{},cos2PIft{};
     float twoPIft{};
     for(std::size_t k{0ull}; k!=this->m_nfreqs; ++k) 
     {    
          const float freq = freqs.m_data[k];
          twoPIft          = C628318530717958647692528676656*freq*f_T;
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
          cos2PIft = 2.0f*Ca*ceph_cosf(twoPIft);
#else 
          cos2PIft = 2.0f*Ca*std::cos(twoPIft);
#endif
          i_loop_term = 0.0f;
          for(std::int32_t i{1}; i!=M; ++i) 
          {
               const float t_i{2.0f*static_cast<float>(i)};
               gamma_i = 314159265358979323846264338328*
                                   (freq*f_T-(t_i-f_M-1.0f)*half_h);
               const float sqr_gamm_i = gamma_i*gamma_i;
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
               const float sin_term{ceph_sin(gamma_i)};
               const float sinc_term{0.5f*((sin_term*sin_term)/sqr_gamm_i)};
               
#else 
               const float sin_term{std::sin(gamma_i)};
               const float sinc_term{0.5f*((sin_term*sin_term)/sqr_gamm_i)};
#endif 
               j_loop_term = 0.0f;
               for(std::int32_t j{0}; j!=M; ++j)  
               {
                    const float t_j{2.0f*static_cast<float>(j)};
                    gamma_j = 314159265358979323846264338328*
                                   (freq*f_T-(t_j-f_M-1.0f)*half_h);
                    const float gamma_ij_sum = gamma_i+gamma_j;
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
                    const float cos1_term = ceph_cosf(gamma_ij_sum);
                    const float cos2_term = Ca*ceph_cosf(gamma_ij_sum-twoPIft);
                    const float sinc1_term= sin_term/gamma_i;
                    const float sinc2_term= ceph_sinf(gamma_j)/gamma_j;
#else 
                    const float cos1_term = std::cos(gamma_ij_sum);
                    const float cos2_term = std::cos(gamma_ij_sum-twoPIft);
                    const float sinc1_term= sin_term/gamma_i;
                    const float sinc2_term= std::sin(gamma_j)/gamma_j;
#endif 
                    const float den =  1.0f+(Ca*Ca)-cos2PIft;
                    const float num =  cos1_term-cos2_term;
                    const float Aij =  num/den;
                    j_loop_term     += Aij*sinc1_term*sinc2_term;
               }
               i_loop_term += std::fma(invf_M,j_loop_term,sinc_term);
          }
          this->m_cpfsk_psd.m_data[k] = lead_term*i_loop_term;
     }
     return (0);
}

std::int32_t 
gms::radiolocation
::cpfsk_signal_t 
::cpfsk_signal_scalar_u8x(const float ak) noexcept(true)
{
     using namespace gms::math;
     constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
     constexpr float C314159265358979323846264338328{3.14159265358979323846264338328f};
     std::size_t     i,j;
     const     float invT{1.0f/(2.0f*this->m_T)};
     float     t_i{};
     float     freq_term{};
     float     phase_term{};
     float     cpfsk_sample{};
     for(i = 0ull;i != ROUND_TO_EIGHT(this->m_T,8ull);i += 8ull) 
     {
          t_i        = static_cast<float>(i+0ull);
          freq_term  = C6283185307179586476925286766559*this->m_fc*t_i;
          phase_term = std::fma(ak,C314159265358979323846264338328*t_i*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
          cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
          cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
          this->m_cpfsk_signal.m_data[i+0ull] = cpfsk_sample;
          t_i        = static_cast<float>(i+1ull);
          freq_term  = C6283185307179586476925286766559*this->m_fc*t_i;
          phase_term = std::fma(ak,C314159265358979323846264338328*t_i*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
          cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
          cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
          this->m_cpfsk_signal.m_data[i+1ull] = cpfsk_sample;
          t_i        = static_cast<float>(i+2ull);
          freq_term  = C6283185307179586476925286766559*this->m_fc*t_i;
          phase_term = std::fma(ak,C314159265358979323846264338328*t_i*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
          cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
          cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
          this->m_cpfsk_signal.m_data[i+2ull] = cpfsk_sample;
          t_i        = static_cast<float>(i+3ull);
          freq_term  = C6283185307179586476925286766559*this->m_fc*t_i;
          phase_term = std::fma(ak,C314159265358979323846264338328*t_i*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
          cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
          cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
          this->m_cpfsk_signal.m_data[i+3ull] = cpfsk_sample;
          t_i        = static_cast<float>(i+4ull);
          freq_term  = C6283185307179586476925286766559*this->m_fc*t_i;
          phase_term = std::fma(ak,C314159265358979323846264338328*t_i*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
          cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
          cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
          this->m_cpfsk_signal.m_data[i+4ull] = cpfsk_sample;
          t_i        = static_cast<float>(i+5ull);
          freq_term  = C6283185307179586476925286766559*this->m_fc*t_i;
          phase_term = std::fma(ak,C314159265358979323846264338328*t_i*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
          cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
          cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
          this->m_cpfsk_signal.m_data[i+5ull] = cpfsk_sample;
          t_i        = static_cast<float>(i+6ull);
          freq_term  = C6283185307179586476925286766559*this->m_fc*t_i;
          phase_term = std::fma(ak,C314159265358979323846264338328*t_i*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
          cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
          cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
          this->m_cpfsk_signal.m_data[i+6ull] = cpfsk_sample;
          t_i        = static_cast<float>(i+7ull);
          freq_term  = C6283185307179586476925286766559*this->m_fc*t_i;
          phase_term = std::fma(ak,C314159265358979323846264338328*t_i*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
          cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
          cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
          this->m_cpfsk_signal.m_data[i+7ull] = cpfsk_sample;
     }

     for(j = i;j != this->m_T; ++j)  
     {
          t_i        = static_cast<float>(j);
          freq_term  = C6283185307179586476925286766559*this->m_fc*t_i;
          phase_term = std::fma(ak,C314159265358979323846264338328*t_i*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
          cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
          cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
          this->m_cpfsk_signal.m_data[j] = cpfsk_sample;
     }
     return (0);

}

std::int32_t 
gms::radiolocation 
::cpfsk_signal_t
::generate_cpfsk_signal_sse(const float ak) noexcept(true) 
{
     using namespace gms::math;
     constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
     constexpr float C314159265358979323846264338328{3.14159265358979323846264338328f};
     constexpr std::size_t LUT_loop_idx_threshold{2257ull};
     std::size_t     i,j;
     const float     invT{1.0f/(2.0f*this->m_T)};
     float           jj;
     const __m128    vinvT{_mm_set1_ps(invT)};
     const __m128    vC6283185307179586476925286766559{_mm_set1_ps(C6283185307179586476925286766559)};
     const __m128    vC314159265358979323846264338328{_mm_set1_ps(C314159265358979323846264338328)};
     const __m128    vfc{_mm_set1_ps(this->m_fc)};
     const __m128    vh{_mm_set1_ps(this->m_h)};
     const __m128    vA{_mm_set1_ps(this->m_A)};
     const __m128    vak{_mm_set1_ps(ak)};
     const __m128    vPhik{_mm_set1_ps(this->m_Phik)};
     if(__builtin_expect(this->m_T>LUT_loop_idx_threshold,0))
     {
          for(i=0ull,jj=0.0f;i != ROUND_TO_FOUR(this->m_T,4ull); i+=4ull,jj+=4.0f) 
          {
              const __m128 vt_i{_mm_setr_ps(jj,jj+1.0f,jj+2.0f,jj+3.0f)};
              const __m128 freq_term{_mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i))};
              const __m128 ph_term1{_mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT))};
              const __m128 ph_term2{_mm_fmadd_ps(vak,ph_term1,vPhik)};
              const __m128 cpfsk_sample{_mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)))};
              _mm_store_ps(&this->m_cpfsk_signal.m_data[i],cpfsk_sample);
          }

          for(j = i;j != this->m_T; ++j)  
          {
              const float t_j        = static_cast<float>(j);
              const float freq_term  = C6283185307179586476925286766559*this->m_fc*t_j;
              const float phase_term = std::fma(ak,C314159265358979323846264338328*t_j*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
              const float cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
              cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
              this->m_cpfsk_signal.m_data[j] = cpfsk_sample;
          }
     }
     else 
     {
          for(i=0ull;i != ROUND_TO_FOUR(this->m_T,4ull); i+=4ull) 
          {
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i],_MM_HINT_T0);
               const __m128 vt_i{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i])};
               const __m128 freq_term{_mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i))};
               const __m128 ph_term1{_mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT))};
               const __m128 ph_term2{_mm_fmadd_ps(vak,ph_term1,vPhik)};
               const __m128 cpfsk_sample{_mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)))};
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i],cpfsk_sample);
          }

          for(j = i;j != this->m_T; ++j)  
          {
              const float t_j        = gms::math::LUT_loop_indices_2257_align16[j];
              const float freq_term  = C6283185307179586476925286766559*this->m_fc*t_j;
              const float phase_term = std::fma(ak,C314159265358979323846264338328*t_j*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
              const float cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
              cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
              this->m_cpfsk_signal.m_data[j] = cpfsk_sample;
          }
     }
     return (0);
}

std::int32_t 
gms::radiolocation
::cpfsk_signal_t
::generate_cpfsk_signal_sse_u4x(const float ak) noexcept(true)                            
{
     using namespace gms::math;
     constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
     constexpr float C314159265358979323846264338328{3.14159265358979323846264338328f};
     constexpr std::size_t LUT_loop_idx_threshold{2257ull};
     std::size_t     i;
     const float     invT{1.0f/(2.0f*this->m_T)};
     float           jj;
     const __m128    vinvT{_mm_set1_ps(invT)};
     const __m128    vC6283185307179586476925286766559{_mm_set1_ps(C6283185307179586476925286766559)};
     const __m128    vC314159265358979323846264338328{_mm_set1_ps(C314159265358979323846264338328)};
     const __m128    vfc{_mm_set1_ps(this->m_fc)};
     const __m128    vh{_mm_set1_ps(this->m_h)};
     const __m128    vA{_mm_set1_ps(this->m_A)};
     const __m128    vak{_mm_set1_ps(ak)};
     const __m128    vPhik{_mm_set1_ps(this->m_Phik)};
     __m128          vt_i,freq_term;
     __m128          ph_term1,ph_term2;
     __m128          cpfsk_sample;
     if(__builtin_expect(this->m_T>LUT_loop_idx_threshold,0))  
     {
           for(i=0ull,jj=0.0f;(i+31ull) < this->m_T;i+=32ull,jj+=32.0f) 
           {
                vt_i         = _mm_setr_ps(jj,jj+1.0f,jj+2.0f,jj+3.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
                vt_i         = _mm_setr_ps(jj+4.0f,jj+5.0f,jj+6.0f,jj+7.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+4ull],cpfsk_sample);
                vt_i         = _mm_setr_ps(jj+8.0f,jj+9.0f,jj+10.0f,jj+11.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+8ull],cpfsk_sample);
                vt_i         = _mm_setr_ps(jj+12.0f,jj+13.0f,jj+14.0f,jj+15.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+12ull],cpfsk_sample);
                vt_i         = _mm_setr_ps(jj+16.0f,jj+17.0f,jj+18.0f,jj+19.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+16ull],cpfsk_sample);
                vt_i         = _mm_setr_ps(jj+20.0f,jj+21.0f,jj+22.0f,jj+23.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+20ull],cpfsk_sample);
                vt_i         = _mm_setr_ps(jj+24.0f,jj+25.0f,jj+26.0f,jj+27.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+24ull],cpfsk_sample);
                vt_i         = _mm_setr_ps(jj+28.0f,jj+29.0f,jj+30.0f,jj+31.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+28ull],cpfsk_sample);
           }

           for(; (i+23ull) < this->m_T; i+=24ull, jj+=24.0f) 
           {
                vt_i         = _mm_setr_ps(jj,jj+1.0f,jj+2.0f,jj+3.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
                vt_i         = _mm_setr_ps(jj+4.0f,jj+5.0f,jj+6.0f,jj+7.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+4ull],cpfsk_sample);
                vt_i         = _mm_setr_ps(jj+8.0f,jj+9.0f,jj+10.0f,jj+11.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+8ull],cpfsk_sample);
                vt_i         = _mm_setr_ps(jj+12.0f,jj+13.0f,jj+14.0f,jj+15.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+12ull],cpfsk_sample);
                vt_i         = _mm_setr_ps(jj+16.0f,jj+17.0f,jj+18.0f,jj+19.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+16ull],cpfsk_sample);
                vt_i         = _mm_setr_ps(jj+20.0f,jj+21.0f,jj+22.0f,jj+23.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+20ull],cpfsk_sample);
           }

           for(; (i+15ull) < this->m_T; i+=16ull, jj+=16.0f) 
           {
                vt_i         = _mm_setr_ps(jj,jj+1.0f,jj+2.0f,jj+3.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
                vt_i         = _mm_setr_ps(jj+4.0f,jj+5.0f,jj+6.0f,jj+7.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+4ull],cpfsk_sample);
                vt_i         = _mm_setr_ps(jj+8.0f,jj+9.0f,jj+10.0f,jj+11.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+8ull],cpfsk_sample);
                vt_i         = _mm_setr_ps(jj+12.0f,jj+13.0f,jj+14.0f,jj+15.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+12ull],cpfsk_sample);
           }

           for(; (i+7ull) < this->m_T; i+=8ull, jj+=8.0f) 
           {
                vt_i         = _mm_setr_ps(jj,jj+1.0f,jj+2.0f,jj+3.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
                vt_i         = _mm_setr_ps(jj+4.0f,jj+5.0f,jj+6.0f,jj+7.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+4ull],cpfsk_sample);
           }

           for(; (i+3ull) < this->m_T; i+=4ull, jj+=4.0f) 
           {
                vt_i         = _mm_setr_ps(jj,jj+1.0f,jj+2.0f,jj+3.0f);
                freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
                ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
                ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
                _mm_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
           }

           for(; (i+0ull) < this->m_T; i+=1ull)  
           {
              const float t_i        = static_cast<float>(i);
              const float freq_term  = C6283185307179586476925286766559*this->m_fc*t_i;
              const float phase_term = std::fma(ak,C314159265358979323846264338328*t_i*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
              const float cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
              cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
              this->m_cpfsk_signal.m_data[i] = cpfsk_sample;
          }
     }
     else 
     {
          for(i = 0ull;(i+31ull) < this->m_T; i+=32ull) 
          {
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+4ull],cpfsk_sample);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+8ull],cpfsk_sample);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+12ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+12ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+16ull],_MM_HINT_T0);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+16ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+16ull],cpfsk_sample);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+20ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+20ull],cpfsk_sample);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+24ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+24ull],cpfsk_sample);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+28ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+28ull],cpfsk_sample);
          }

          for(; (i+23ull) < this->m_T; i+=24ull) 
          {
                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+4ull],cpfsk_sample);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+8ull],cpfsk_sample);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+12ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+12ull],cpfsk_sample);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+16ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+16ull],cpfsk_sample);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+20ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+20ull],cpfsk_sample);
          }

          for(; (i+15ull) < this->m_T; i+=16ull) 
          {
                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+4ull],cpfsk_sample);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+8ull],cpfsk_sample);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+12ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+12ull],cpfsk_sample);
          }

          for(; (i+7ull) < this->m_T; i+=8ull) 
          {
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+4ull],cpfsk_sample);
          }

          for(; (i+3ull) < this->m_T; i+=4ull) 
          {
               vt_i         = _mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]);
               freq_term    = _mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vfc,vt_i));
               ph_term1     = _mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt_i,vinvT));
               ph_term2     = _mm_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm_mul_ps(vA,_mm_cos_ps(_mm_add_ps(freq_term,ph_term2)));
               _mm_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
          }

          for(; (i+0ull) < this->m_T; i+=1ull)  
          {
              const float t_i        = gms::math::LUT_loop_indices_2257_align16[i];
              const float freq_term  = C6283185307179586476925286766559*this->m_fc*t_i;
              const float phase_term = std::fma(ak,C314159265358979323846264338328*t_i*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
              const float cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
              cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
              this->m_cpfsk_signal.m_data[i] = cpfsk_sample;
          }
     }

     return (0);
}

std::int32_t 
gms::radiolocation
::cpfsk_signal_t
::generate_cpfsk_signal_avx(const float ak,
                            const bool do_const_prefetch) noexcept(true) 
{
     using namespace gms::math;         
     const float     invT{1.0f/(2.0f*this->m_T)};    
     const __m256    vinvT{_mm256_set1_ps(invT)};   
     const __m256    vfc{_mm256_set1_ps(this->m_fc)};
     const __m256    vh{_mm256_set1_ps(this->m_h)};
     const __m256    vA{_mm256_set1_ps(this->m_A)};
     const __m256    vak{_mm256_set1_ps(ak)};
     const __m256    vPhik{_mm256_set1_ps(this->m_Phik)};
     __m256          vC6283185307179586476925286766559;
     __m256          vC314159265358979323846264338328;
     __m256          c0_8;
     __m256          vt_i,freq_term;
     __m256          ph_term1,ph_term2;
     __m256          cpfsk_sample;
     if(!do_const_prefetch) 
     {
         vC6283185307179586476925286766559 = _mm256_set1_ps(6.283185307179586476925286766559f);
         vC314159265358979323846264338328  = _mm256_set1_ps(3.14159265358979323846264338328f);
         c0_8                              = _mm256_setr_ps(0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f);
     }
     else 
     {
         __ATTR_ALIGN__(32) const float prefetched_constants[32] = {
                                0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                                0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f}; // the last eight floats serve as a padding
          _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
          c0_8                              = _mm256_load_ps(&prefetched_constants[0]);
          vC6283185307179586476925286766559 = _mm256_load_ps(&prefetched_constants[8]);
          vC314159265358979323846264338328  = _mm256_load_ps(&prefetched_constants[16]);
     }
     std::size_t     i,j;
     float           jj;
     constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
     constexpr float C314159265358979323846264338328{3.14159265358979323846264338328f};
     constexpr std::size_t LUT_loop_idx_threshold{2257ull}; 
     if(this->m_T>LUT_loop_idx_threshold) 
     {
          for(i=0ull,jj=0.0f; i != ROUND_TO_EIGHT(this->m_T,8ull); i+=8ull,jj+=8.0f) 
          {
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c0_8);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
              _mm256_store_ps(&this->m_cpfsk_signal.m_data[i],cpfsk_sample);
          }

          for(j = i;j != this->m_T; ++j)  
          {
              const float t_j        = static_cast<float>(j);
              const float freq_term  = C6283185307179586476925286766559*this->m_fc*t_j;
              const float phase_term = std::fma(ak,C314159265358979323846264338328*t_j*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
              const float cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
              cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
              this->m_cpfsk_signal.m_data[j] = cpfsk_sample;
          }
     }
     else
     {
          for(i=0ull; i != ROUND_TO_EIGHT(this->m_T,8ull); i+=8ull) 
          {
                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i],_MM_HINT_T0);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i],cpfsk_sample);
          }

          for(j = i;j != this->m_T; ++j)  
          {
              const float t_j        = gms::math::LUT_loop_indices_2257_align16[j];
              const float freq_term  = C6283185307179586476925286766559*this->m_fc*t_j;
              const float phase_term = std::fma(ak,C314159265358979323846264338328*t_j*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
              const float cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
              cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
              this->m_cpfsk_signal.m_data[j] = cpfsk_sample;
          }
     }

     return (0);
}

std::int32_t 
gms::radiolocation
::cpfsk_signal_t
::generate_cpfsk_signal_avx_u8x(const float ak,
                                const bool do_const_prefetch) noexcept(true) 
{
     using namespace gms::math;         
     const float     invT{1.0f/(2.0f*this->m_T)};    
     const __m256    vinvT{_mm256_set1_ps(invT)};   
     const __m256    vfc{_mm256_set1_ps(this->m_fc)};
     const __m256    vh{_mm256_set1_ps(this->m_h)};
     const __m256    vA{_mm256_set1_ps(this->m_A)};
     const __m256    vak{_mm256_set1_ps(ak)};
     const __m256    vPhik{_mm256_set1_ps(this->m_Phik)};
     __m256          vC6283185307179586476925286766559;
     __m256          vC314159265358979323846264338328;
     __m256          c0_7,c8_15,c16_23,c24_31;
     __m256          c32_39,c40_47,c48_55,c56_63;
     __m256          vt_i,freq_term;
     __m256          ph_term1,ph_term2;
     __m256          cpfsk_sample;
     if(do_const_prefetch) 
     {
          c0_7   = _mm256_setr_ps(0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f);
          c8_15  = _mm256_setr_ps(8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f);
          c16_23 = _mm256_setr_ps(16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f);
          c24_31 = _mm256_setr_ps(24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f);
          c32_39 = _mm256_setr_ps(32.0f,33.0f,34.0f,35.0f,36.0f,37.0f,38.0f,39.0f);
          c40_47 = _mm256_setr_ps(40.0f,41.0f,42.0f,43.0f,44.0f,45.0f,46.0f,47.0f);
          c48_55 = _mm256_setr_ps(48.0f,49.0f,50.0f,51.0f,52.0f,53.0f,54.0f,55.0f);
          c56_63 = _mm256_setr_ps(56.0f,57.0f,58.0f,59.0f,60.0f,61.0f,62.0f,63.0f);
          vC6283185307179586476925286766559 = _mm256_set1_ps(6.283185307179586476925286766559f);
          vC314159265358979323846264338328  = _mm256_set1_ps(3.14159265358979323846264338328f);
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
                        3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                        3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                        3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                        3.14159265358979323846264338328f, 3.14159265358979323846264338328f
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
          vC314159265358979323846264338328  = _mm256_load_ps(&prefetched_constants[72]);
     }
     std::size_t     i,j;
     float           jj;
     constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
     constexpr float C314159265358979323846264338328{3.14159265358979323846264338328f};
     constexpr std::size_t LUT_loop_idx_threshold{2257ull}; 
     if(__builtin_expect(this->m_T>LUT_loop_idx_threshold,0))
     {
          for(i=0ull,jj=0.0f; (i+63ull) < this->m_T; i+=64ull,jj+=64.0f) 
          {
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c0_7);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c8_15);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+8ull],cpfsk_sample);
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c16_23);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+16ull],cpfsk_sample);
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c24_31);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+24ull],cpfsk_sample);
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c32_39);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+32ull],cpfsk_sample);
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c40_47);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+40ull],cpfsk_sample);
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c48_55);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+48ull],cpfsk_sample);
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c56_63);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+56ull],cpfsk_sample);
          }

          for(; (i+47ull) < this->m_T; i+=48ull,jj+48.0f) 
          {
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c0_7);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c8_15);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+8ull],cpfsk_sample);
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c16_23);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+16ull],cpfsk_sample);
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c24_31);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+24ull],cpfsk_sample);
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c32_39);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+32ull],cpfsk_sample);
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c40_47);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+40ull],cpfsk_sample);
          }

          for(; (i+31ull) < this->m_T; i+=32ull,jj+32.0f) 
          {
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c0_7);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c8_15);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+8ull],cpfsk_sample);
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c16_23);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+16ull],cpfsk_sample);
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c24_31);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+24ull],cpfsk_sample);
          }

          for(; (i+15ull) < this->m_T; i+=16ull,jj+=16.0f) 
          {
                vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c0_7);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
               vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c8_15);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+8ull],cpfsk_sample);
          }

          for(; (i+7ull) < this->m_T; i+=8ull,jj+=8.0f) 
          {
                vt_i         = _mm256_add_ps(_mm256_set1_ps(jj),c0_7);
               freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
               ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
               ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
          }

          for(; (i+0ull) < this->m_T; i+=1ull) 
          {               
              const float t_i        = static_cast<float>(i);
              const float freq_term  = C6283185307179586476925286766559*this->m_fc*t_i;
              const float phase_term = std::fma(ak,C314159265358979323846264338328*t_i*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
              const float cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
              cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
              this->m_cpfsk_signal.m_data[i] = cpfsk_sample;
          
          }
     }
     else  
     {
          for(i=0ull; (i+63ull) < this->m_T; i+=64ull) 
          {
                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+8ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+16ull],_MM_HINT_T0);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+16ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+16ull],cpfsk_sample);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+24ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+24ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+32ull],_MM_HINT_T0);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+32ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+32ull],cpfsk_sample);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+40ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+40ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+48ull],_MM_HINT_T0);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+48ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+48ull],cpfsk_sample);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+56ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+56ull],cpfsk_sample);
          }

          for(;(i+47ull) < this->m_T; i+=48ull) 
          {
                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+8ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+16ull],_MM_HINT_T0);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+16ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+16ull],cpfsk_sample);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+24ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+24ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+32ull],_MM_HINT_T0);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+32ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+32ull],cpfsk_sample);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+40ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+40ull],cpfsk_sample);
          }

          for(;(i+31ull) < this->m_T; i+=32ull) 
          {
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+8ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+16ull],_MM_HINT_T0);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+16ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+16ull],cpfsk_sample);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+24ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+24ull],cpfsk_sample);
          }

          for(;(i+15ull) < this->m_T; i+=16ull) 
          {
                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+8ull],cpfsk_sample);
          }

          for(;(i+7ull) < this->m_T; i+=8ull) 
          {
                vt_i         = _mm256_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]);
                freq_term    = _mm256_mul_ps(vC6283185307179586476925286766559,_mm256_mul_ps(vfc,vt_i));
                ph_term1     = _mm256_mul_ps(vC314159265358979323846264338328,_mm256_mul_ps(vt_i,vinvT));
                ph_term2     = _mm256_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm256_mul_ps(vA,_mm256_cos_ps(_mm256_add_ps(freq_term,ph_term2)));
               _mm256_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
          }

          for(;(i+0ull) < this->m_T; i+=1ull) 
          {
              const float t_i        = gms::math::LUT_loop_indices_2257_align16[i];
              const float freq_term  = C6283185307179586476925286766559*this->m_fc*t_i;
              const float phase_term = std::fma(ak,C314159265358979323846264338328*t_i*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
              const float cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
              cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
              this->m_cpfsk_signal.m_data[i] = cpfsk_sample;
          }
     }

     return (0);
}

std::int32_t 
gms::radiolocation
::cpfsk_signal_t
::generate_cpfsk_signal_avx512(const float ak,
                               const bool do_const_prefetch) noexcept(true) 
{
     using namespace gms::math;         
     const float     invT{1.0f/(2.0f*this->m_T)};    
     const __m512    vinvT{_mm512_set1_ps(invT)};   
     const __m512    vfc{_mm512_set1_ps(this->m_fc)};
     const __m512    vh{_mm512_set1_ps(this->m_h)};
     const __m512    vA{_mm512_set1_ps(this->m_A)};
     const __m512    vak{_mm512_set1_ps(ak)};
     const __m512    vPhik{_mm512_set1_ps(this->m_Phik)};
     __m512          vC6283185307179586476925286766559;
     __m512          vC314159265358979323846264338328;
     __m512          c0_15;
     __m512          vt_i,freq_term;
     __m512          ph_term1,ph_term2;
     __m512          cpfsk_sample;
     if(!do_const_prefetch) 
     {
         c0_15                             = _mm512_setr_ps(0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f);
         vC6283185307179586476925286766559 = _mm512_set1_ps(6.283185307179586476925286766559f);
         vC314159265358979323846264338328  = _mm512_set1_ps(3.14159265358979323846264338328f);
     } 
     else 
     {
         __ATTR_ALIGN__(64) const float prefetched_constants[48] = 
         {
                                0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,  
                                3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                                3.14159265358979323846264338328f, 3.14159265358979323846264338328f
         };
          _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
          _mm_prefetch((const char*)&prefetched_constants[32],_MM_HINT_T0);
          c0_15                             = _mm512_load_ps(&prefetched_constants[0]);
          vC6283185307179586476925286766559 = _mm512_load_ps(&prefetched_constants[16]);
          vC314159265358979323846264338328  = _mm512_load_ps(&prefetched_constants[32]);
     }
     std::size_t     i,j;
     float           jj;
     constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
     constexpr float C314159265358979323846264338328{3.14159265358979323846264338328f};
     constexpr std::size_t LUT_loop_idx_threshold{2257ull};
     if(__builtin_expect(this->m_T>LUT_loop_idx_threshold,0)) 
     {
          for(i=0ull,jj=0.0f; i != ROUND_TO_SIXTEEN(this->m_T,16ull); i+=16ull,jj+=16.0f) 
          {
               vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c0_15);
               freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
               ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
               ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
               cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
              _mm512_store_ps(&this->m_cpfsk_signal.m_data[i],cpfsk_sample);
          }

          for(j = i;j != this->m_T; ++j)  
          {
              const float t_j        = static_cast<float>(j);
              const float freq_term  = C6283185307179586476925286766559*this->m_fc*t_j;
              const float phase_term = std::fma(ak,C314159265358979323846264338328*t_j*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
              const float cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
              cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
              this->m_cpfsk_signal.m_data[j] = cpfsk_sample;
          }
     }
     else 
     {
          for(i=0ull;i != ROUND_TO_SIXTEEN(this->m_T,16ull);i+=16ull) 
          {
                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i],cpfsk_sample);
          }

          for(j = i;j != this->m_T; ++j)  
          {
              const float t_j        = gms::math::LUT_loop_indices_2257_align16[j];
              const float freq_term  = C6283185307179586476925286766559*this->m_fc*t_j;
              const float phase_term = std::fma(ak,C314159265358979323846264338328*t_j*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
              const float cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
              cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
              this->m_cpfsk_signal.m_data[j] = cpfsk_sample;
          }
     }

     return (0);
}

std::int32_t 
gms::radiolocation
::cpfsk_signal_t
::generate_cpfsk_signal_avx512_u8x(const float ak,
                                   const bool do_const_prefetch) noexcept(true) 
{
     using namespace gms::math;         
     const float     invT{1.0f/(2.0f*this->m_T)};    
     const __m512    vinvT{_mm512_set1_ps(invT)};   
     const __m512    vfc{_mm512_set1_ps(this->m_fc)};
     const __m512    vh{_mm512_set1_ps(this->m_h)};
     const __m512    vA{_mm512_set1_ps(this->m_A)};
     const __m512    vak{_mm512_set1_ps(ak)};
     const __m512    vPhik{_mm512_set1_ps(this->m_Phik)};
     __m512          vC6283185307179586476925286766559;
     __m512          vC314159265358979323846264338328;
     __m512 c0_15,c16_31;
     __m512 c32_47,c48_63;
     __m512 c64_79,c80_95;
     __m512 c96_111,c112_127;
     __m512 vt_i,freq_term;
     __m512 ph_term1,ph_term2;
     __m512 cpfsk_sample;
     if(!do_const_prefetch) 
     {
         c0_15   = _mm512_setr_ps(0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f);
         c16_31  = _mm512_setr_ps(16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f,24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f);
         c32_47  = _mm512_setr_ps(32.0f,33.0f,34.0f,35.0f,36.0f,37.0f,38.0f,39.0f,40.0f,41.0f,42.0f,43.0f,44.0f,45.0f,46.0f,47.0f);
         c48_63  = _mm512_setr_ps(48.0f,49.0f,50.0f,51.0f,52.0f,53.0f,54.0f,55.0f,56.0f,57.0f,58.0f,59.0f,60.0f,61.0f,62.0f,63.0f);
         c64_79  = _mm512_setr_ps(64.0f,65.0f,66.0f,67.0f,68.0f,69.0f,70.0f,71.0f,72.0f,73.0f,74.0f,75.0f,76.0f,77.0f,78.0f,79.0f);
         c80_95  = _mm512_setr_ps(80.0f,81.0f,82.0f,83.0f,84.0f,85.0f,86.0f,87.0f,88.0f,89.0f,90.0f,91.0f,92.0f,93.0f,94.0f,95.0f);
         c96_111 = _mm512_setr_ps(96.0f,97.0f,98.0f,99.0f,100.0f,101.0f,102.0f,103.0f,104.0f,105.0f,106.0f,107.0f,108.0f,109.0f,110.0f,111.0f);
         c112_127= _mm512_setr_ps(112.0f,113.0f,114.0f,115.0f,116.0f,117.0f,118.0f,119.0f,120.0f,121.0f,122.0f,123.0f,124.0f,125.0f,126.0f,127.0f);
         vC6283185307179586476925286766559 = _mm512_set1_ps(6.283185307179586476925286766559f);
         vC314159265358979323846264338328  = _mm512_set1_ps(3.14159265358979323846264338328f);
     }
     else 
     {
         __ATTR_ALIGN__(64) const float prefetched_constants[160] = 
         {
                            0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f,
                            16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f,24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f,
                            32.0f,33.0f,34.0f,35.0f,36.0f,37.0f,38.0f,39.0f,40.0f,41.0f,42.0f,43.0f,44.0f,45.0f,46.0f,47.0f,
                            48.0f,49.0f,50.0f,51.0f,52.0f,53.0f,54.0f,55.0f,56.0f,57.0f,58.0f,59.0f,60.0f,61.0f,62.0f,63.0f,
                            64.0f,65.0f,66.0f,67.0f,68.0f,69.0f,70.0f,71.0f,72.0f,73.0f,74.0f,75.0f,76.0f,77.0f,78.0f,79.0f,
                            80.0f,81.0f,82.0f,83.0f,84.0f,85.0f,86.0f,87.0f,88.0f,89.0f,90.0f,91.0f,92.0f,93.0f,94.0f,95.0f,
                            96.0f,97.0f,98.0f,99.0f,100.0f,101.0f,102.0f,103.0f,104.0f,105.0f,106.0f,107.0f,108.0f,109.0f,110.0f,111.0f,
                            112.0f,113.0f,114.0f,115.0f,116.0f,117.0f,118.0f,119.0f,120.0f,121.0f,122.0f,123.0f,124.0f,125.0f,126.0f,127.0f,
                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,  
                            3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                            3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                            3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                            3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                            3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                            3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                            3.14159265358979323846264338328f, 3.14159265358979323846264338328f,
                            3.14159265358979323846264338328f, 3.14159265358979323846264338328f
         };
        _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[32],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[48],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[64],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[80],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[96],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[112],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[128],_MM_HINT_T0);
        _mm_prefetch((const char*)&prefetched_constants[144],_MM_HINT_T0);
        c0_15                             =   _mm512_load_ps(&prefetched_constants[0]);
        c16_31                            =   _mm512_load_ps(&prefetched_constants[16]);
        c32_47                            =   _mm512_load_ps(&prefetched_constants[32]);
        c48_63                            =   _mm512_load_ps(&prefetched_constants[48]);
        c64_79                            =   _mm512_load_ps(&prefetched_constants[64]);
        c80_95                            =   _mm512_load_ps(&prefetched_constants[80]);
        c96_111                           =   _mm512_load_ps(&prefetched_constants[96]);
        c112_127                          =   _mm512_load_ps(&prefetched_constants[112]);
        vC6283185307179586476925286766559 =   _mm512_load_ps(&prefetched_constants[128]);
        vC314159265358979323846264338328  =   _mm512_load_ps(&prefetched_constants[144]);
     }
     std::size_t     i,j;
     float           jj;
     constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
     constexpr float C314159265358979323846264338328{3.14159265358979323846264338328f};
     constexpr std::size_t LUT_loop_idx_threshold{2257ull};
      if(__builtin_expect(this->m_T>LUT_loop_idx_threshold,0))
      {
            for(i=0ull,jj=0.0f;(i+127ull) < this->m_T; i+=128ull,jj+=128.0f) 
            {
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c0_15);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c16_31);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+16ull],cpfsk_sample);
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c32_47);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+32ull],cpfsk_sample);
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c48_63);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+48ull],cpfsk_sample);
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c64_79);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+64ull],cpfsk_sample);
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c80_95);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+80ull],cpfsk_sample);
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c96_111);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+96ull],cpfsk_sample);
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c112_127);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+112ull],cpfsk_sample);
            }

            for(; (i+95ull) < this->m_T; i+=96ull,jj+=96.0f) 
            {
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c0_15);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c16_31);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+16ull],cpfsk_sample);
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c32_47);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+32ull],cpfsk_sample);
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c48_63);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+48ull],cpfsk_sample);
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c64_79);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+64ull],cpfsk_sample);
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c80_95);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+80ull],cpfsk_sample);
            }

            for(; (i+63ull) < this->m_T; i+=64ull,jj+=64.0f) 
            {
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c0_15);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c16_31);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+16ull],cpfsk_sample);
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c32_47);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+32ull],cpfsk_sample);
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c48_63);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+48ull],cpfsk_sample);
            }

            for(; (i+31ull) < this->m_T; i+=32ull,jj+=32.0f) 
            {
                 vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c0_15);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c16_31);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+16ull],cpfsk_sample);
            }

            for(; (i+15ull) < this->m_T; i+=16ull,jj+=16.0f) 
            {
                vt_i         = _mm512_add_ps(_mm512_set1_ps(jj),c0_15);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
                _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
            }

            for(; (i+0ull) < this->m_T; i+=1ull) 
            {               
                const float t_i        = static_cast<float>(i);
                const float freq_term  = C6283185307179586476925286766559*this->m_fc*t_i;
                const float phase_term = std::fma(ak,C314159265358979323846264338328*t_i*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
                const float cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
                cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
                this->m_cpfsk_signal.m_data[i] = cpfsk_sample;
          
            }
      }
      else 
      {
            for(i=0ull; (i+127ull) < this->m_T; i+=128ull) 
            {
                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+16ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+16ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+32ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+32ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+48ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+48ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+64ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+64ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+80ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+80ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+96ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+96ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+112ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+112ull],cpfsk_sample);
            }

            for(; (i+95ull) < this->m_T; i+=96ull) 
            {
                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+16ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+16ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+32ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+32ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+48ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+48ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+64ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+64ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+80ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+80ull],cpfsk_sample);
            }

            for(; (i+63ull) < this->m_T; i+=64ull) 
            {
                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+16ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+16ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+32ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+32ull],cpfsk_sample);
               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+48ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+48ull],cpfsk_sample);
            }

            for(; (i+31ull) < this->m_T; i+=32ull) 
            {
                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+16ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+16ull],cpfsk_sample);
            }

            for(; (i+15ull) < this->m_T; i+=16ull) 
            {
                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
                vt_i         = _mm512_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]);
                freq_term    = _mm512_mul_ps(vC6283185307179586476925286766559,_mm512_mul_ps(vfc,vt_i));
                ph_term1     = _mm512_mul_ps(vC314159265358979323846264338328,_mm512_mul_ps(vt_i,vinvT));
                ph_term2     = _mm512_fmadd_ps(vak,ph_term1,vPhik);
                cpfsk_sample = _mm512_mul_ps(vA,_mm512_cos_ps(_mm512_add_ps(freq_term,ph_term2)));
               _mm512_store_ps(&this->m_cpfsk_signal.m_data[i+0ull],cpfsk_sample);
            }

            for(;(i+0ull) < this->m_T; i+=1ull) 
            {
              const float t_i        = gms::math::LUT_loop_indices_2257_align16[i];
              const float freq_term  = C6283185307179586476925286766559*this->m_fc*t_i;
              const float phase_term = std::fma(ak,C314159265358979323846264338328*t_i*invT,this->m_Phik);
#if (CPFSK_SIGNAL_USE_CEPHES) == 1
              const float cpfsk_sample = this->m_A*ceph_cosf(freq_term+phase_term);
#else 
              cpfsk_sample = this->m_A*std::cos(freq_term+phase_term);
#endif 
              this->m_cpfsk_signal.m_data[i] = cpfsk_sample;
           }
      }

      return (0);
}

auto 
gms::radiolocation
::operator<<(std::ostream &os,
             const cpfsk_signal_t &rhs)->std::ostream &
{
     std::cout << typeid(rhs).name() << "Begin: object state dump." << std::endl;
     std::cout << "m_T          :"   << rhs.m_T                     << std::endl;
     std::cout << "m_k          :"   << rhs.m_k                     << std::endl;
     std::cout << "m_A          : "  << std::fixed << std::setprecision(7) << rhs.m_A << std::endl;
     std::cout << "m_h          : "  << std::fixed << std::setprecision(7) << rhs.m_h << std::endl;
     std::cout << "m_fc         : "  << std::fixed << std::setprecision(7) << rhs.m_I_fc << std::endl;
     std::cout << "m_ph0        : "  << std::fixed << std::setprecision(7) << rhs.m_ph0 << std::endl;
     std::cout << "m_Phik       : "  << std::fixed << std::setprecision(7) << rhs.m_Phik << std::endl;
     std::cout << "CPFKS signal -- samples:" <<"\n\n";
     for(std::size_t i{0ull}; i!=this->m_T; ++i) 
     {
          os << std::fixed << std::setprecision(7) << rhs.m_cpfsk_signal.m_data[i] << std::endl;
     }
     std::cout << typeid(rhs).name() << "End: object state dump." << std::endl;
     return (os);
}

