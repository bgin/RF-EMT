#include <fstream>
#include "GMS_rectangular_waveform.h"
#include "GMS_sse_memset.h"
#if (RECTANGULAR_WAVEFORM_USE_CEPHES) == 0
#include <cmath>
#endif 

#if (RECTANGULAR_WAVEFORM_USE_CEPHES) == 1

     
__ATTR_ALWAYS_INLINE__
static inline
float ceph_cosf( const float xx ) {
/* Note, these constants are for a 32-bit significand: */
/*
static float DP1 =  0.7853851318359375;
static float DP2 =  1.30315311253070831298828125e-5;
static float DP3 =  3.03855025325309630e-11;
static float lossth = 65536.;
*/

/* These are for a 24-bit significand: */
constexpr float PIO4F =  0.7853981633974483096f;
constexpr float FOPI  = 1.27323954473516f;
constexpr float DP1 = 0.78515625f;
constexpr float DP2 = 2.4187564849853515625e-4f;
constexpr float DP3 = 3.77489497744594108e-8f;
constexpr float lossth = 8192.f;
constexpr float T24M1 = 16777215.f;
float x, y, z;
int j, sign;

/* make argument positive */
sign = 1;
x = xx;
if( x < 0 )
	x = -x;

if( x > T24M1 )
	{
	//mtherr( "cosf", TLOSS );
	return(0.0);
	}

j = FOPI * x; /* integer part of x/PIO4 */
y = j;
/* integer and fractional part modulo one octant */
if( j & 1 )	/* map zeros to origin */
	{
	j += 1;
	y += 1.0;
	}
j &= 7;
if( j > 3)
	{
	j -=4;
	sign = -sign;
	}

if( j > 1 )
	sign = -sign;

if( x > lossth )
	{
	//mtherr( "cosf", PLOSS );
	x = x - y * PIO4F;
	}
else
/* Extended precision modular arithmetic */
	x = ((x - y * DP1) - y * DP2) - y * DP3;

z = x * x;

if( (j==1) || (j==2) )
	{
	y = (((-1.9515295891E-4 * z
	     + 8.3321608736E-3) * z
	     - 1.6666654611E-1) * z * x)
	     + x;
	}
else
	{
	y = ((  2.443315711809948E-005 * z
	  - 1.388731625493765E-003) * z
	  + 4.166664568298827E-002) * z * z;
	y -= 0.5 * z;
	y += 1.0;
	}
if(sign < 0)
	y = -y;
return( y );
}

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
::rectangular_waveform_t
::rectangular_waveform_t(const std::size_t n_samples,
                         const std::uint32_t n_waves,
                         const float T,
                         const float rho,
                         const float A)
:
__n_samples__{n_samples},
__n_waves__{n_waves},
__T__{T},
__rho__{rho},
__A__{A},
__rw_samples__{__n_samples__}
{
   
}


gms::radiolocation
::rectangular_waveform_t
::rectangular_waveform_t(rectangular_waveform_t &&other)
:
__n_samples__{std::move(other.__n_samples__)},
__n_waves__{std::move(other.__n_waves__)},
__T__{std::move(other.__T__)},
__rho__{std::move(other.__rho__)},
__A__{std::move(other.__A__)},
__rw_samples__{std::move(other.__rw_samples__)}
{
     
}


gms::radiolocation
::rectangular_waveform_t
::~rectangular_waveform_t()
{

}

gms::radiolocation
::rectangular_waveform_t &
gms::radiolocation
::rectangular_waveform_t
::operator=(rectangular_waveform_t &&other)
{
    if(this==&other) { return (*this);}
    this->__n_samples__          = std::move(other.__n_samples__);
    this->__n_waves__            = std::move(other.__n_waves__);
    this->__T__                  = std::move(other.__T__);
    this->__rho__                = std::move(other.__rho__);
    this->__A__                  = std::move(other.__A__);
    this->__rw_samples__.operator=(std::move(other.__rw_samples__));
    return (*this);
}

void
gms::radiolocation 
::rectangular_waveform_t
::init_storage(const float filler)
{
#if (INIT_BY_STD_FILL) == 0
     using namespace gms::common;
	 sse_memset_unroll8x_ps(&this->__rw_samples__.m_data[0],filler,this->__n_samples__);
#else 
     std::fill(this->__rw_samples__.m_data,this->__rw_samples__+this->__n_samples__,filler);
#endif 
}

void
gms::radiolocation
::rectangular_waveform_t
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
::rectangular_waveform_t 
::fourier_series_expansion(const std::uint32_t which_index)
{
    constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
    constexpr float C314159265358979323846264338328{3.14159265358979323846264338328f};
    const float rho_over_T{this->__rho__/this->__T__};
    const float f{1.0f/this->__T__};
    const float twoA_over_PI{(this->__A__+this->__A__)*0.318309886183790671537767526745f};
    const std::uint32_t n_samples{static_cast<std::uint32_t>(this->__n_samples__)};
    float sum;
    switch (which_index)
    {
         case 0 : 
             for(std::uint32_t __i{0}; __i != n_samples; ++__i) 
             {
                   const float t__i{static_cast<float>(__i)};
                   sum = 0.0f;
                   for(std::uint32_t __n{1}; __n <= this->__n_waves__; ++__n) 
                   {
                        const float t__n{static_cast<float>(__n)};
                        const float i_t__n{1.0f/t__n};
                        const float sin_arg{C314159265358979323846264338328*t__n*rho_over_T};
                        const float cos_arg{C6283185307179586476925286766559*t__n*f*t__i};
#if (RECTANGULAR_WAVEFORM_USE_CEPHES) == 1
                        const float t_sin{i_t__n*ceph_sinf(sin_arg)};
                        const float t_cos{ceph_cosf(cos_arg)};
#else 
                        const float t_sin{std::sin(sin_arg)};
                        const float t_cos{std::cos(cos_arg)};
#endif 
                        const float mul_part{t_sin*t_cos};
                        sum += mul_part;
                   }
                  this->__rw_samples__.m_data[__i] = (this->__A__*rho_over_T)+twoA_over_PI*sum;
             }
         break;
         case 1 : 
             const float f_n_samples{1.0f/static_cast<float>(n_samples)};
             for(std::uint32_t __i{0}; __i != n_samples; ++__i) 
             {
                   const float t__i{static_cast<float>(__i)};
                   const float t__t{t__i*this->__T__*f_n_samples};
                   sum = 0.0f;
                   for(std::uint32_t __n{1}; __n <= this->__n_waves__; ++__n) 
                   {
                        const float t__n{static_cast<float>(__n)};
                        const float i_t__n{1.0f/t__n};
                        const float sin_arg{C314159265358979323846264338328*t__n*rho_over_T};
                        const float cos_arg{C6283185307179586476925286766559*t__n*f*t__t};
#if (RECTANGULAR_WAVEFORM_USE_CEPHES) == 1
                        const float t_sin{i_t__n*ceph_sinf(sin_arg)};
                        const float t_cos{ceph_cosf(cos_arg)};
#else 
                        const float t_sin{std::sin(sin_arg)};
                        const float t_cos{std::cos(cos_arg)};
#endif 
                        const float mul_part{t_sin*t_cos};
                        sum += mul_part;
                   }
                  this->__rw_samples__.m_data[__i] = (this->__A__*rho_over_T)+twoA_over_PI*sum;
             }
         break;
         default : 
            return;
    }
    
}

void 
gms::radiolocation
::rectangular_waveform_t 
::fourier_series_expansion_optim(const std::uint32_t which_index)
{
    if(__builtin_expect(this->__n_waves__>256u,0)) { return;}

    constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
    constexpr float C314159265358979323846264338328{3.14159265358979323846264338328f};
    __ATTR_ALIGN__(64) float sin_buf[256];
    const float rho_over_T{this->__rho__/this->__T__};
    const float f{1.0f/this->__T__};
    const float twoA_over_PI{(this->__A__+this->__A__)*0.318309886183790671537767526745f};
    const std::uint32_t n_samples{static_cast<std::uint32_t>(this->__n_samples__)};
    std::int32_t idx{-1};

    for(std::uint32_t __n{1}; __n <= this->__n_waves__; ++__n) 
    {
         //idx += 1;
         const float t__n{static_cast<float>(__n)};
         const float i_t__n{1.0f/t__n};
         const float sin_arg{C314159265358979323846264338328*t__n*rho_over_T};
#if (RECTANGULAR_WAVEFORM_USE_CEPHES) == 1
         const float t_sin{i_t__n*ceph_sinf(sin_arg)};
#else 
         const float t_sin{std::sin(sin_arg)};
#endif 
         sin_buf[idx++] = t_sin; 
    }

    float sum;
    switch (which_index)
    {
         case 0 : 
             for(std::uint32_t __i{0}; __i != n_samples; ++__i) 
             {
                   const float t__i{static_cast<float>(__i)};
                   sum = 0.0f;
                   idx = -1;
                   for(std::uint32_t __n{1}; __n <= this->__n_waves__; ++__n) 
                   {
                        const float t__n{static_cast<float>(__n)};
                        const float cos_arg{C6283185307179586476925286766559*t__n*f*t__i};
#if (RECTANGULAR_WAVEFORM_USE_CEPHES) == 1
                        
                        const float t_cos{ceph_cosf(cos_arg)};
#else 
                       
                        const float t_cos{std::cos(cos_arg)};
#endif 
                        const float mul_part{sin_buf[idx++]*t_cos};
                        sum += mul_part;
                   }
                  this->__rw_samples__.m_data[__i] = (this->__A__*rho_over_T)+twoA_over_PI*sum;
             }
         break;
         case 1 : 
             const float f_n_samples{1.0f/static_cast<float>(n_samples)};
             for(std::uint32_t __i{0}; __i != n_samples; ++__i) 
             {
                   const float t__i{static_cast<float>(__i)};
                   const float t__t{t__i*this->__T__*f_n_samples};
                   sum = 0.0f;
                   idx = -1;
                   for(std::uint32_t __n{1}; __n <= this->__n_waves__; ++__n) 
                   {
                        const float t__n{static_cast<float>(__n)};
                        const float cos_arg{C6283185307179586476925286766559*t__n*f*t__t};
#if (RECTANGULAR_WAVEFORM_USE_CEPHES) == 1
                       
                        const float t_cos{ceph_cosf(cos_arg)};
#else 
                       
                        const float t_cos{std::cos(cos_arg)};
#endif 
                        const float mul_part{sin_buf[idx++]*t_cos};
                        sum += mul_part;
                   }
                  this->__rw_samples__.m_data[__i] = (this->__A__*rho_over_T)+twoA_over_PI*sum;
             }
         break;
         default : 
            return;
    }
    
}
