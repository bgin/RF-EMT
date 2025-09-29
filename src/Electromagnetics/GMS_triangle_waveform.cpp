
#include <fstream>
#include <iomanip>
#include "GMS_triangle_waveform.h"
#include "GMS_sse_memset.h"
#if (TRIANGLE_WAVEFORM_USE_CEPHES) == 0
#include <cmath>
#endif

#if (TRIANGLE_WAVEFORM_USE_CEPHES) == 1

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
register unsigned long j;
register int sign;
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
::triangle_waveform_t
::triangle_waveform_t(const std::size_t n_samples,
                    const std::uint32_t n_waves,
                    const float A,
                    const float T) 
:
__n_samples__{n_samples},
__n_waves__{n_waves},
__A__{A},
__T__{T},
__tw_samples__{darray_r4_t(__n_samples__)}
{

}

gms::radiolocation
::triangle_waveform_t
::triangle_waveform_t(triangle_waveform_t &&other) noexcept(true)
:
__n_samples__{std::move(other.__n_samples__)},
__n_waves__{std::move(other.__n_waves__)},
__A__{std::move(other.__A__)},
__T__{std::move(other.__T__)},
__tw_samples__{std::move(other.__tw_samples__)}
{

}

gms::radiolocation
::triangle_waveform_t
::~triangle_waveform_t()
{

}

gms::radiolocation
::triangle_waveform_t & 
gms::radiolocation
::triangle_waveform_t
::operator=(triangle_waveform_t &&other)
{
    if(__builtin_expect(this==&other,0)) {return (*this);}
    this->__n_samples__ = std::move(other.__n_samples__);
    this->__n_waves__   = std::move(other.__n_waves__);
    this->__A__         = std::move(other.__A__);
    this->__T__         = std::move(other.__T__);
    this->__tw_samples__.operator=(std::move(other.__tw_samples__));
    return (*this);
}

void 
gms::radiolocation 
::triangle_waveform_t
::init_storage(const float filler)
{
#if (INIT_BY_STD_FILL) == 0
     using namespace gms::common;
	 sse_memset_unroll8x_ps(&this->__tw_samples__.m_data[0],filler,this->__n_samples__);
#else 
     std::fill(this->__tw_samples__.m_data,this->__sw_samples__.m_data+this->__n_samples__,filler);
#endif
}

void 
gms::radiolocation 
::triangle_waveform_t
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

void 
gms::radiolocation
::triangle_waveform_t
::fourier_series_expansion()
{
    static constexpr float k_coff_buf[64] = {
                   -0.012345679104328155518f,
                    0.008264462463557720184f,
                   -0.005917159840464591980f,
                    0.004444444552063941956f,
                   -0.003460207721218466759f,
                    0.002770083025097846985f,
                   -0.002267573727294802666f,
                    0.001890359213575720787f,
                   -0.001599999959580600262f,
                    0.001371742109768092632f,
                   -0.001189060625620186329f,
                    0.001040582777932286263f,
                   -0.000918273639399558306f,
                    0.000816326530184596777f,
                   -0.000730460218619555235f,
                    0.000657462223898619413f,
                   -0.000594883982557803392f,
                    0.000540832872502505779f,
                   -0.000493827159516513348f,
                    0.000452693522674962878f,
                   -0.000416493130614981055f,
                    0.000384467508411034942f,
                   -0.000355998578015714884f,
                    0.000330578506691381335f,
                   -0.000307787006022408605f,
                    0.000287273782305419445f,
                   -0.000268744974164292216f,
                    0.000251952646067366004f,
                   -0.000236686391872353852f,
                    0.000222766757360659540f,
                   -0.000210039914236404002f,
                    0.000198373338207602501f,
                   -0.000187652462045662105f,
                    0.000177777779754251242f,
                   -0.000168662503710947931f,
                    0.000160230731125921011f,
                   -0.000152415785123594105f,
                    0.000145158948726020753f,
                   -0.000138408300699666142f,
                    0.000132117842440493405f,
                   -0.000126246683066710830f,
                    0.000120758362754713744f,
                   -0.000115620299766305834f,
                    0.000110803324787411839f,
                   -0.000106281222542747855f,
                    0.000102030404377728701f,
                   -0.000098029602668248117f,
                    0.000094259594334289432f,
                   -0.000090702946181409061f,
                    0.000087343869381584227f,
                   -0.000084168001194484532f,
                    0.000081162244896404445f,
                   -0.000078314667916856706f,
                    0.000075614363595377654f,
                   -0.000073051356594078243f,
                    0.000070616479206364602f,
                   -0.000068301342253107578f,
                    0.000066098218667320907f,
                   -0.000063999999838415533f,
                    0.000062000122852623463f};

     constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
     constexpr float lead_coff{-8.0f/(C6283185307179586476925286766559*
                                C6283185307179586476925286766559)};
     const std::uint32_t n_samples{static_cast<std::uint32_t>(this->__n_samples__)};
     //const float period{8.0f/static_cast<float>(n_samples)};
     const float invT{this->__T__/static_cast<float>(n_samples)};
     float sum;
     if(__builtin_expect(this->__n_waves__<=64,1))
     {
          for(std::uint32_t __i{0}; __i != n_samples; ++__i) 
          {
                const float t__t{static_cast<float>(__i)};
                sum = 0.0f;
                for(std::uint32_t __k{1}; __k != this->__n_waves__; ++__k) 
                {
                     const float t__k{static_cast<float>(__k)};
                     const float two_k_1{2.0f*t__k-1.0f};
                     const float arg{C6283185307179586476925286766559*two_k_1*t__t*invT};
#if (TRIANGLE_WAVEFORM_USE_CEPHES) == 1 
                  const float t_sin{ceph_sinf(arg)};
#else 
                  const float t_sin{std::sin(arg)};
#endif 
                  const float sample{k_coff_buf[__k]*this->__A__*t_sin};
                  sum += sample;
                }
                this->__tw_samples__.m_data[__i] = lead_coff*sum;
          }
     }
     else 
     {
           for(std::uint32_t __i{0}; __i != n_samples; ++__i) 
           {
                const float t__t{static_cast<float>(__i)};
                sum = 0.0f;
                for(std::uint32_t __k{1}; __k != this->__n_waves__; ++__k) 
                {
                     const float t__k{static_cast<float>(__k)};
                     const float two_k_1{2.0f*t__k-1.0f};
                     const float tmp{two_k_1*two_k_1};
                     const float one_over_tk1{std::pow(-1.0f,t__k)/tmp};
                     const float arg{C6283185307179586476925286766559*two_k_1*t__t*invT};
#if (TRIANGLE_WAVEFORM_USE_CEPHES) == 1 
                     const float t_sin{ceph_sinf(arg)};
#else 
                     const float t_sin{std::sin(arg)};
#endif 
                     const float sample{one_over_tk1*this->__A__*t_sin};
                     sum += sample;
                }
                this->__tw_samples__.m_data[__i] = lead_coff*sum;
          }
     }

    
}

auto 
gms::radiolocation
::operator<<(std::ostream &os, 
             gms::radiolocation::triangle_waveform_t &rhs)->std::ostream &
{
    std::cout << typeid(rhs).name() << "Begin: object state dump." << std::endl;
    std::cout << "__n_samples__ :"  << rhs.__n_samples__ << std::endl;
    std::cout << "__n_waves__   :"  << rhs.__n_waves__   << std::endl;
    std::cout << "__A__         :"  << rhs.__A__         << std::endl;
    std::cout << "__T__         :"  << rhs.__T__         << std::endl;
    std::cout << "Triangle-Wave samples (Fourier-expansion) of " << rhs.__n_waves__ << "waves" << std::endl;
     for(std::size_t __i{0ull}; __i != rhs.__n_samples__; ++__i)
    {
         os << std::fixed << std::setprecision(7) << 
                       rhs.__tw_samples__.m_data[__i] << std::endl;
    }
    std::cout << typeid(rhs).name() << "End: object state dump." << std::endl;
    return (os);
}
