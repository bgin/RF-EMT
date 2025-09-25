

#include <fstream>
#include <iomanip>
#include "GMS_sawtooth_waveform.h"
#include "GMS_sse_memset.h"
#if (SAWTOOTH_WAVEFORM_USE_CEPHES) == 0
#include <cmath>
#endif

#if (SAWTOOTH_WAVEFORM_USE_CEPHES) == 1

__ATTR_ALWAYS_INLINE__
static inline
float ceph_sinf( const float xx ) {
constexpr float FOPI  = 1.27323954473516;
constexpr float PIO4F =  0.7853981633974483096;
/* Note, these constants are for a 32-bit significand: */
/*
static float DP1 =  0.7853851318359375;
static float DP2 =  1.30315311253070831298828125e-5;
static float DP3 =  3.03855025325309630e-11;
static float lossth = 65536.;
*/

/* These are for a 24-bit significand: */
constexpr float DP1 = 0.78515625;
constexpr float DP2 = 2.4187564849853515625e-4;
constexpr float DP3 = 3.77489497744594108e-8;
constexpr float lossth = 8192.;
constexpr float T24M1 = 16777215.;

float sincof[] = {
-1.9515295891E-4,
 8.3321608736E-3,
-1.6666654611E-1
};
float coscof[] = {
 2.443315711809948E-005,
-1.388731625493765E-003,
 4.166664568298827E-002
};
float *p;
float x, y, z;
unsigned long j;
int sign;
sign = 1;
x = xx;
if( xx < 0 )
	{
	sign = -1;
	x = -xx;
	}
if( x > T24M1 )
	{
	//mtherr( "sinf", TLOSS );
	return(0.0);
	}
j = FOPI * x; /* integer part of x/(PI/4) */
y = j;
/* map zeros to origin */
if( j & 1 )
	{
	j += 1;
	y += 1.0;
	}
j &= 7; /* octant modulo 360 degrees */
/* reflect in x axis */
if( j > 3)
	{
	sign = -sign;
	j -= 4;
	}

if( x > lossth )
	{
	//mtherr( "sinf", PLOSS );
	x = x - y * PIO4F;
	}
else
	{
/* Extended precision modular arithmetic */
	x = ((x - y * DP1) - y * DP2) - y * DP3;
	}
/*einits();*/
z = x * x;
if( (j==1) || (j==2) )
	{
/* measured relative error in +/- pi/4 is 7.8e-8 */
/*
	y = ((  2.443315711809948E-005 * z
	  - 1.388731625493765E-003) * z
	  + 4.166664568298827E-002) * z * z;
*/
	p = coscof;
	y = *p++;
	y = y * z + *p++;
	y = y * z + *p++;
	y *= z * z;
	y -= 0.5 * z;
	y += 1.0;
	}
else
	{
/* Theoretical relative error = 3.8e-9 in [-pi/4, +pi/4] */
/*
	y = ((-1.9515295891E-4 * z
	     + 8.3321608736E-3) * z
	     - 1.6666654611E-1) * z * x;
	y += x;
*/
	p = sincof;
	y = *p++;
	y = y * z + *p++;
	y = y * z + *p++;
	y *= z * x;
	y += x;
	}
/*einitd();*/
if(sign < 0)
	y = -y;
return( y);
}

#endif 

gms::radiolocation
::sawtooth_waveform_t
::sawtooth_waveform_t(const std::size_t n_samples,
                      const std::uint32_t n_waves,
                      const float A,
                      const float T)
:
__n_samples__{n_samples},
__n_waves__{n_waves},
__A__{A},
__T__{T},
__sw_samples__{darray_r4_t(__n_samples__)}
{

}

gms::radiolocation
::sawtooth_waveform_t
::sawtooth_waveform_t(sawtooth_waveform_t &&other) noexcept(true)
:
__n_samples__{std::move(other.__n_samples__)},
__n_waves__{std::move(other.__n_waves__)},
__A__{std::move(other.__A__)},
__T__{std::move(other.__T__)},
__sw_samples__{std::move(other.__sw_samples__)}
{

}

gms::radiolocation
::sawtooth_waveform_t
::~sawtooth_waveform_t() 
{

}

gms::radiolocation
::sawtooth_waveform_t & 
gms::radiolocation
::sawtooth_waveform_t
::operator=(sawtooth_waveform_t &&other)
{
    if(__builtin_expect(this==&other,0)) {return (*this);}
    this->__n_samples__ = std::move(other.__n_samples__);
    this->__n_waves__   = std::move(other.__n_waves__);
    this->__A__         = std::move(other.__A__);
    this->__T__         = std::move(other.__T__);
    this->__sw_samples__.operator=(std::move(other.__sw_samples__));
    return (*this);
}

void 
gms::radiolocation 
::sawtooth_waveform_t
::init_storage(const float filler)
{
#if (INIT_BY_STD_FILL) == 0
     using namespace gms::common;
	 sse_memset_unroll8x_ps(&this->__sw_samples__.m_data[0],filler,this->__n_samples__);
#else 
     std::fill(this->__sw_samples__.m_data,this->__sw_samples__.m_data+this->__n_samples__,filler);
#endif
}

void 
gms::radiolocation 
::sawtooth_waveform_t
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
::sawtooth_waveform_t
::fourier_series_expansion(const std::uint32_t which_dir)
{
       if(__builtin_expect(this->__n_waves__>128,0)) { return (-2);}

static constexpr float k_coff_buf[128] = {
                -1.0000000f,
                 0.5000000f,
                -0.3333333f,
                 0.2500000f,
                -0.2000000f,
                 0.1666667f,
                -0.1428571f,
                 0.1250000f,
                -0.1111111f,
                 0.1000000f,
                -0.0909091f,
                 0.0833333f,
                -0.0769231f,
                 0.0714286f,
                -0.0666667f,
                 0.0625000f,
                -0.0588235f,
                 0.0555556f,
                -0.0526316f,
                 0.0500000f,
                -0.0476190f,
                 0.0454545f,
                -0.0434783f,
                 0.0416667f,
                -0.0400000f,
                 0.0384615f,
                -0.0370370f,
                 0.0357143f,
                -0.0344828f,
                 0.0333333f,
                -0.0322581f,
                 0.0312500f,
                -0.0303030f,
                 0.0294118f,
                -0.0285714f,
                 0.0277778f,
                -0.0270270f,
                 0.0263158f,
                -0.0256410f,
                 0.0250000f,
                -0.0243902f,
                 0.0238095f,
                -0.0232558f,
                 0.0227273f,
                -0.0222222f,
                 0.0217391f,
                -0.0212766f,
                 0.0208333f,
                -0.0204082f,
                 0.0200000f,
                -0.0196078f,
                 0.0192308f,
                -0.0188679f,
                 0.0185185f,
                -0.0181818f,
                 0.0178571f,
                -0.0175439f,
                 0.0172414f,
                -0.0169492f,
                 0.0166667f,
                -0.0163934f,
                 0.0161290f,
                -0.0158730f,
                 0.0156250f,
                -0.0153846f,
                 0.0151515f,
                -0.0149254f,
                 0.0147059f,
                -0.0144928f,
                 0.0142857f,
                -0.0140845f,
                 0.0138889f,
                -0.0136986f,
                 0.0135135f,
                -0.0133333f,
                 0.0131579f,
                -0.0129870f,
                 0.0128205f,
                -0.0126582f,
                 0.0125000f,
                -0.0123457f,
                 0.0121951f,
                -0.0120482f,
                 0.0119048f,
                -0.0117647f,
                 0.0116279f,
                -0.0114943f,
                 0.0113636f,
                -0.0112360f,
                 0.0111111f,
                -0.0109890f,
                 0.0108696f,
                -0.0107527f,
                 0.0106383f,
                -0.0105263f,
                 0.0104167f,
                -0.0103093f,
                 0.0102041f,
                -0.0101010f,
                 0.0100000f,
                -0.0099010f,
                 0.0098039f,
                -0.0097087f,
                 0.0096154f,
                -0.0095238f,
                 0.0094340f,
                -0.0093458f,
                 0.0092593f,
                -0.0091743f,
                 0.0090909f,
                -0.0090090f,
                 0.0089286f,
                -0.0088496f,
                 0.0087719f,
                -0.0086957f,
                 0.0086207f,
                -0.0085470f,
                 0.0084746f,
                -0.0084034f,
                 0.0083333f,
                -0.0082645f,
                 0.0081967f,
                -0.0081301f,
                 0.0080645f,
                -0.0080000f,
                 0.0079365f,
                -0.0078740f,
                 0.0078125f};
        
        constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
        const float     lead_coff{2.0f*0.318309886183790671537767526745f};
        const std::uint32_t n_samples{static_cast<std::uint32_t>(this->__n_samples__)};
        const float invT{this->__T__/static_cast<float>(n_samples)};
        float sum;
        int32_t idx;
        switch (which_dir)
        {
            case 0 : // direct sawtooth 
                
                          
                     for(std::uint32_t __t{0}; __t != n_samples; ++__t) 
                     {
                         const float t__t{static_cast<float>(__t)};
                         sum = 0.0f;
                         idx = 0;
                         for(std::uint32_t __k{1}; __k != this->__n_waves__; ++__k) 
                         {
                             const float t__k{static_cast<float>(__k)};
                             const float arg{C6283185307179586476925286766559*t__k*t__t*invT};
#if (SAWTOOTH_WAVEFORM_USE_CEPHES) == 1 
                             const float t_sin{ceph_sinf(arg)};
#else 
                             const float t_sin{std::sin(arg)};
#endif 
                             const float sample{k_coff_buf[idx++]*this->__A__*t_sin};
                             sum += sample;
                         }
                         this->__sw_samples__.m_data[__t] = -lead_coff*sum;
                     }

                
                
           break;
           case 1 : // reverse sawtooth
               
                            
                     for(std::uint32_t __t{0}; __t != n_samples; ++__t) 
                     {
                         const float t__t{static_cast<float>(__t)};
                         sum = 0.0f;
                         idx = 0;
                         for(std::uint32_t __k{1}; __k != this->__n_waves__; ++__k) 
                         {
                             const float t__k{static_cast<float>(__k)};
                             const float arg{C6283185307179586476925286766559*t__k*t__t*invT};
#if (SAWTOOTH_WAVEFORM_USE_CEPHES) == 1 
                             const float t_sin{ceph_sinf(arg)};
#else 
                             const float t_sin{std::sin(arg)};
#endif 
                             const float sample{k_coff_buf[idx++]*this->__A__*t_sin};
                             sum += sample;
                         }
                          this->__sw_samples__.m_data[__t] = lead_coff*sum;
                     }
           break;     
           default : 
                 return (-1);
     }
     return (0);
             
}

std::int32_t 
gms::radiolocation
::sawtooth_waveform_t
::fourier_series_expansion_u4x(const std::uint32_t which_dir)
{
       if(__builtin_expect(this->__n_waves__>128,0)) { return (-2);}

static constexpr float k_coff_buf[128] = {
                -1.0000000f,
                 0.5000000f,
                -0.3333333f,
                 0.2500000f,
                -0.2000000f,
                 0.1666667f,
                -0.1428571f,
                 0.1250000f,
                -0.1111111f,
                 0.1000000f,
                -0.0909091f,
                 0.0833333f,
                -0.0769231f,
                 0.0714286f,
                -0.0666667f,
                 0.0625000f,
                -0.0588235f,
                 0.0555556f,
                -0.0526316f,
                 0.0500000f,
                -0.0476190f,
                 0.0454545f,
                -0.0434783f,
                 0.0416667f,
                -0.0400000f,
                 0.0384615f,
                -0.0370370f,
                 0.0357143f,
                -0.0344828f,
                 0.0333333f,
                -0.0322581f,
                 0.0312500f,
                -0.0303030f,
                 0.0294118f,
                -0.0285714f,
                 0.0277778f,
                -0.0270270f,
                 0.0263158f,
                -0.0256410f,
                 0.0250000f,
                -0.0243902f,
                 0.0238095f,
                -0.0232558f,
                 0.0227273f,
                -0.0222222f,
                 0.0217391f,
                -0.0212766f,
                 0.0208333f,
                -0.0204082f,
                 0.0200000f,
                -0.0196078f,
                 0.0192308f,
                -0.0188679f,
                 0.0185185f,
                -0.0181818f,
                 0.0178571f,
                -0.0175439f,
                 0.0172414f,
                -0.0169492f,
                 0.0166667f,
                -0.0163934f,
                 0.0161290f,
                -0.0158730f,
                 0.0156250f,
                -0.0153846f,
                 0.0151515f,
                -0.0149254f,
                 0.0147059f,
                -0.0144928f,
                 0.0142857f,
                -0.0140845f,
                 0.0138889f,
                -0.0136986f,
                 0.0135135f,
                -0.0133333f,
                 0.0131579f,
                -0.0129870f,
                 0.0128205f,
                -0.0126582f,
                 0.0125000f,
                -0.0123457f,
                 0.0121951f,
                -0.0120482f,
                 0.0119048f,
                -0.0117647f,
                 0.0116279f,
                -0.0114943f,
                 0.0113636f,
                -0.0112360f,
                 0.0111111f,
                -0.0109890f,
                 0.0108696f,
                -0.0107527f,
                 0.0106383f,
                -0.0105263f,
                 0.0104167f,
                -0.0103093f,
                 0.0102041f,
                -0.0101010f,
                 0.0100000f,
                -0.0099010f,
                 0.0098039f,
                -0.0097087f,
                 0.0096154f,
                -0.0095238f,
                 0.0094340f,
                -0.0093458f,
                 0.0092593f,
                -0.0091743f,
                 0.0090909f,
                -0.0090090f,
                 0.0089286f,
                -0.0088496f,
                 0.0087719f,
                -0.0086957f,
                 0.0086207f,
                -0.0085470f,
                 0.0084746f,
                -0.0084034f,
                 0.0083333f,
                -0.0082645f,
                 0.0081967f,
                -0.0081301f,
                 0.0080645f,
                -0.0080000f,
                 0.0079365f,
                -0.0078740f,
                 0.0078125f};
        
        constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
        const float     lead_coff{2.0f*0.318309886183790671537767526745f};
        const std::uint32_t n_samples{static_cast<std::uint32_t>(this->__n_samples__)};
        const float invT{this->__T__/static_cast<float>(n_samples)};
        float sum0;
        float sum1;
        float sum2;
        float sum3;
        float t__i_0{-4.0f};
        float t__i_1{1.0f};
        float t__i_2{2.0f};
        float t__i_3{3.0f};
        std::uint32_t __i;
        std::uint32_t __j;
        int32_t idx;
        switch (which_dir)
        {
            case 0 : // direct sawtooth 
                
                          
                     for( __i = 0; __i != ROUND_TO_FOUR(n_samples,4); __i += 4) 
                     {
                         //const float t__t{static_cast<float>(__t)};
                         idx = 0;
                         t__i_0 += 4.0f;
                         sum0 = 0.0f;
                         for(std::uint32_t __k{1}; __k != this->__n_waves__; ++__k) 
                         {
                             const float t__k{static_cast<float>(__k)};
                             const float arg{C6283185307179586476925286766559*t__k*t__i_0*invT};
#if (SAWTOOTH_WAVEFORM_USE_CEPHES) == 1 
                             const float t_sin{ceph_sinf(arg)};
#else 
                             const float t_sin{std::sin(arg)};
#endif 
                             const float sample{k_coff_buf[idx++]*this->__A__*t_sin};
                             sum0 += sample;
                         }
                         this->__sw_samples__.m_data[__i+0] = -lead_coff*sum0;
                         idx = 0;
                         t__i_1 += 4.0f;
                         sum1   = 0.0f;
                         for(std::uint32_t __k{1}; __k != this->__n_waves__; ++__k) 
                         {
                             const float t__k{static_cast<float>(__k)};
                             const float arg{C6283185307179586476925286766559*t__k*t__i_1*invT};
#if (SAWTOOTH_WAVEFORM_USE_CEPHES) == 1 
                             const float t_sin{ceph_sinf(arg)};
#else 
                             const float t_sin{std::sin(arg)};
#endif 
                             const float sample{k_coff_buf[idx++]*this->__A__*t_sin};
                             sum1 += sample;
                         }
                         this->__sw_samples__.m_data[__i+1] = -lead_coff*sum1;
                         idx = 0;
                         t__i_2 += 4.0f;
                         sum2   = 0.0f;
                         for(std::uint32_t __k{1}; __k != this->__n_waves__; ++__k) 
                         {
                             const float t__k{static_cast<float>(__k)};
                             const float arg{C6283185307179586476925286766559*t__k*t__i_2*invT};
#if (SAWTOOTH_WAVEFORM_USE_CEPHES) == 1 
                             const float t_sin{ceph_sinf(arg)};
#else 
                             const float t_sin{std::sin(arg)};
#endif 
                             const float sample{k_coff_buf[idx++]*this->__A__*t_sin};
                             sum2 += sample;
                         }
                         this->__sw_samples__.m_data[__i+2] = -lead_coff*sum2;
                         idx = 0;
                         t__i_3 += 4.0f;
                         sum3   =  0.0f;
                         for(std::uint32_t __k{1}; __k != this->__n_waves__; ++__k) 
                         {
                             const float t__k{static_cast<float>(__k)};
                             const float arg{C6283185307179586476925286766559*t__k*t__i_3*invT};
#if (SAWTOOTH_WAVEFORM_USE_CEPHES) == 1 
                             const float t_sin{ceph_sinf(arg)};
#else 
                             const float t_sin{std::sin(arg)};
#endif 
                             const float sample{k_coff_buf[idx++]*this->__A__*t_sin};
                             sum3 += sample;
                         }
                         this->__sw_samples__.m_data[__i+3] = -lead_coff*sum3;
                     }

                     for(__j = __i; __j != n_samples; ++__j) 
                     {
                         const float t__t{static_cast<float>(__j)};
                         sum0 = 0.0f;
                         idx = 0;
                         for(std::uint32_t __k{1}; __k != this->__n_waves__; ++__k) 
                         {
                             const float t__k{static_cast<float>(__k)};
                             const float arg{C6283185307179586476925286766559*t__k*t__t*invT};
#if (SAWTOOTH_WAVEFORM_USE_CEPHES) == 1 
                             const float t_sin{ceph_sinf(arg)};
#else 
                             const float t_sin{std::sin(arg)};
#endif 
                             const float sample{k_coff_buf[idx++]*this->__A__*t_sin};
                             sum0 += sample;
                         }
                         this->__sw_samples__.m_data[__j] = -lead_coff*sum0;
                     }

                
                
           break;
           case 1 : // reverse sawtooth
               
                     for( __i = 0; __i != ROUND_TO_FOUR(n_samples,4); __i += 4) 
                     {
                         //const float t__t{static_cast<float>(__t)};
                         idx = 0;
                         t__i_0 += 4.0f;
                         sum0 = 0.0f;
                         for(std::uint32_t __k{1}; __k != this->__n_waves__; ++__k) 
                         {
                             const float t__k{static_cast<float>(__k)};
                             const float arg{C6283185307179586476925286766559*t__k*t__i_0*invT};
#if (SAWTOOTH_WAVEFORM_USE_CEPHES) == 1 
                             const float t_sin{ceph_sinf(arg)};
#else 
                             const float t_sin{std::sin(arg)};
#endif 
                             const float sample{k_coff_buf[idx++]*this->__A__*t_sin};
                             sum0 += sample;
                         }
                         this->__sw_samples__.m_data[__i+0] = lead_coff*sum0;
                         idx = 0;
                         t__i_1 += 4.0f;
                         sum1   = 0.0f;
                         for(std::uint32_t __k{1}; __k != this->__n_waves__; ++__k) 
                         {
                             const float t__k{static_cast<float>(__k)};
                             const float arg{C6283185307179586476925286766559*t__k*t__i_1*invT};
#if (SAWTOOTH_WAVEFORM_USE_CEPHES) == 1 
                             const float t_sin{ceph_sinf(arg)};
#else 
                             const float t_sin{std::sin(arg)};
#endif 
                             const float sample{k_coff_buf[idx++]*this->__A__*t_sin};
                             sum1 += sample;
                         }
                         this->__sw_samples__.m_data[__i+1] = lead_coff*sum1;
                         idx = 0;
                         t__i_2 += 4.0f;
                         sum2   = 0.0f;
                         for(std::uint32_t __k{1}; __k != this->__n_waves__; ++__k) 
                         {
                             const float t__k{static_cast<float>(__k)};
                             const float arg{C6283185307179586476925286766559*t__k*t__i_2*invT};
#if (SAWTOOTH_WAVEFORM_USE_CEPHES) == 1 
                             const float t_sin{ceph_sinf(arg)};
#else 
                             const float t_sin{std::sin(arg)};
#endif 
                             const float sample{k_coff_buf[idx++]*this->__A__*t_sin};
                             sum2 += sample;
                         }
                         this->__sw_samples__.m_data[__i+2] = lead_coff*sum2;
                         idx = 0;
                         t__i_3 += 4.0f;
                         sum3   =  0.0f;
                         for(std::uint32_t __k{1}; __k != this->__n_waves__; ++__k) 
                         {
                             const float t__k{static_cast<float>(__k)};
                             const float arg{C6283185307179586476925286766559*t__k*t__i_3*invT};
#if (SAWTOOTH_WAVEFORM_USE_CEPHES) == 1 
                             const float t_sin{ceph_sinf(arg)};
#else 
                             const float t_sin{std::sin(arg)};
#endif 
                             const float sample{k_coff_buf[idx++]*this->__A__*t_sin};
                             sum3 += sample;
                         }
                         this->__sw_samples__.m_data[__i+3] = lead_coff*sum3;
                     }

                     for(__j = __i; __j != n_samples; ++__j) 
                     {
                         const float t__t{static_cast<float>(__j)};
                         sum0 = 0.0f;
                         idx = 0;
                         for(std::uint32_t __k{1}; __k != this->__n_waves__; ++__k) 
                         {
                             const float t__k{static_cast<float>(__k)};
                             const float arg{C6283185307179586476925286766559*t__k*t__t*invT};
#if (SAWTOOTH_WAVEFORM_USE_CEPHES) == 1 
                             const float t_sin{ceph_sinf(arg)};
#else 
                             const float t_sin{std::sin(arg)};
#endif 
                             const float sample{k_coff_buf[idx++]*this->__A__*t_sin};
                             sum0 += sample;
                         }
                         this->__sw_samples__.m_data[__j] = lead_coff*sum0;
                     }
     
                    
           break;     
           default : 
                 return (-1);
     }
     return (0);
             
}

auto 
gms::radiolocation
::operator<<(std::ostream &os, 
             gms::radiolocation::sawtooth_waveform_t &rhs)->std::ostream &
{
    std::cout << typeid(rhs).name() << "Begin: object state dump." << std::endl;
    std::cout << "__n_samples__ :"  << rhs.__n_samples__ << std::endl;
    std::cout << "__n_waves__   :"  << rhs.__n_waves__   << std::endl;
    std::cout << "__A__         :"  << rhs.__A__         << std::endl;
    std::cout << "__T__         :"  << rhs.__T__         << std::endl;
    std::cout << "Sawtooth-Wave samples (Fourier-expansion) of " << rhs.__n_waves__ << " :waves" << std::endl;
     for(std::size_t __i{0ull}; __i != rhs.__n_samples__; ++__i)
    {
         os << std::fixed << std::setprecision(7) << 
                       rhs.__sw_samples__.m_data[__i] << std::endl;
    }
    std::cout << typeid(rhs).name() << "End: object state dump." << std::endl;
    return (os);
}

