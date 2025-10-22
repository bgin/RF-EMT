

#include <fstream>
#include "GMS_trapezoid_waveform_fp16.h"
#include "GMS_sse_memset.h"
#if (TRAPEZOID_WAVEFORM_FP16_USE_CEPHES) == 0
#include <cmath>
#endif 

#if (TRAPEZOID_WAVEFORM_FP16_USE_CEPHES) == 1

#if !defined (DENORMAL)
    #define DENORMAL 1
#endif


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

__ATTR_ALWAYS_INLINE__
static inline
float ceph_frexpf( float x, int *pw2 ) {
union
  {
    float y;
    unsigned short i[2];
  } u;
int i, k;
unsigned short *q;
u.y = x;
q = &u.i[1];
/* find the exponent (power of 2) */

i  = ( *q >> 7) & 0xff;
if( i == 0 ) {
    if( u.y == 0.0 ) {
	*pw2 = 0;
	 return(0.0);
    }
/* Number is denormal or zero */
#if DENORMAL
/* Handle denormal number. */
	do{
		
	   u.y *= 2.0;
	   i -= 1;
	   k  = ( *q >> 7) & 0xff;
	}
	while( k == 0 );
	i = i + k;
#else
	*pw2 = 0;
	return( 0.0 );
#endif /* DENORMAL */
	}
i -= 0x7e;
*pw2 = i;
*q &= 0x807f;	/* strip all exponent bits */
*q |= 0x3f00;	/* mantissa between 0.5 and 1 */
return( u.y );
}


__ATTR_ALWAYS_INLINE__
static inline
float ceph_ldexpf( float x, int pw2 ) {
#define MEXP 255
union
  {
    float y;
    unsigned short i[2];
  } u;
unsigned short *q;
int e;
u.y = x;
q = &u.i[1];
while( (e = ( *q >> 7) & 0xff) == 0 ) {
	if( u.y == (float )0.0 ) {
		return( 0.0 );
	}
/* Input is denormal. */
	if( pw2 > 0 ) {
	    u.y *= 2.0;
	    pw2 -= 1;
	}
	if( pw2 < 0 ) {
	    if( pw2 < -24 )
		return( 0.0 );
	    u.y *= 0.5;
	    pw2 += 1;
	}
	if( pw2 == 0 )
		return(u.y);
	}
e += pw2;

/* Handle overflow */
if( e > MEXP )
	{
	return( 3.4028235e+38 );
	}

*q &= 0x807f;

/* Handle denormalized results */
if( e < 1 )
	{
#if DENORMAL
	if( e < -24 )
		return( 0.0 );
	*q |= 0x80; /* Set LSB of exponent. */
	/* For denormals, significant bits may be lost even
	   when dividing by 2.  Construct 2^-(1-e) so the result
	   is obtained with only one multiplication.  */
	u.y *= ceph_ldexpf(1.0f, e - 1);
	return(u.y);
#else
	return( 0.0 );
#endif
	}
*q |= (e & 0xff) << 7;
return(u.y);
}


__ATTR_ALWAYS_INLINE__
static inline
float ceph_sqrtf( const float xx ) {
float f, x, y;
int e;
f = xx;
if( f <= 0.0 )
	{
	if( f < 0.0 )
		//mtherr( "sqrtf", DOMAIN );
	return( 0.0 );
	}

x = ceph_frexpf( f, &e );	/* f = x * 2**e,   0.5 <= x < 1.0 */
/* If power of 2 is odd, double x and decrement the power of 2. */
if( e & 1 )
	{
	x = x + x;
	e -= 1;
	}

e >>= 1;	/* The power of 2 of the square root. */

if( x > 1.41421356237 )
	{
/* x is between sqrt(2) and 2. */
	x = x - 2.0;
	y =
	((((( -9.8843065718E-4 * x
	  + 7.9479950957E-4) * x
	  - 3.5890535377E-3) * x
	  + 1.1028809744E-2) * x
	  - 4.4195203560E-2) * x
	  + 3.5355338194E-1) * x
	  + 1.41421356237E0;
	goto sqdon;
	}

if( x > 0.707106781187 )
	{
/* x is between sqrt(2)/2 and sqrt(2). */
	x = x - 1.0;
	y =
	((((( 1.35199291026E-2 * x
	  - 2.26657767832E-2) * x
	  + 2.78720776889E-2) * x
	  - 3.89582788321E-2) * x
	  + 6.24811144548E-2) * x
	  - 1.25001503933E-1) * x * x
	  + 0.5 * x
	  + 1.0;
	goto sqdon;
	}

/* x is between 0.5 and sqrt(2)/2. */
x = x - 0.5;
y =
((((( -3.9495006054E-1 * x
  + 5.1743034569E-1) * x
  - 4.3214437330E-1) * x
  + 3.5310730460E-1) * x
  - 3.5354581892E-1) * x
  + 7.0710676017E-1) * x
  + 7.07106781187E-1;

sqdon:
y = ceph_ldexpf( y, e );  /* y = y * 2**e */
return( y);
}


__ATTR_ALWAYS_INLINE__
static inline 
float ceph_asinf(const float xx) {
float a, x, z;
int sign, flag;
x = xx;
if( x > 0 ) {
  sign = 1;
  a = x;
}
else
    {
  sign = -1;
  a = -x;
}

if( a > 1.0 ) {
    return( 0.0 );
}

if( a < 1.0e-4 ) {
    z = a;
    goto done;
}

if( a > 0.5 ) {
    z = 0.5 * (1.0 - a);
    x = ceph_sqrtf( z );
    flag = 1;
}
else
    {
    x = a;
    z = x * x;
    flag = 0;
}
z =
(((( 4.2163199048E-2 * z
  + 2.4181311049E-2) * z
  + 4.5470025998E-2) * z
  + 7.4953002686E-2) * z
  + 1.6666752422E-1) * z * x
  + x;

if( flag != 0 ) {
    z = z + z;
    z = 1.5707963267948966192 - z;
}
done:
if( sign < 0 )
	z = -z;
return( z );
}

__ATTR_ALWAYS_INLINE__
static inline 
float ceph_acosf(const float x ) {
if( x < -1.0f )
	goto domerr;

if( x < -0.5f) 
	return( 3.14159265358979323846f
                - 2.0f * ceph_asinf(ceph_sqrtf(0.5f*(1.0f+x)) ) );
if( x > 1.0f ) {
domerr:
	return( 0.0f );
}
if( x > 0.5f )
	return( 2.0f * ceph_asinf(  ceph_sqrtf(0.5f*(1.0f-x) ) ) );
	
return(1.5707963267948966192f - ceph_asinf(x) );
}

#endif 


gms::radiolocation
::trapezoid_waveform_fp16_t
::trapezoid_waveform_fp16_t(const std::size_t n_samples,
                       const std::uint32_t n_waves,
                       const std::uint32_t n_param_a,
                       const std::uint32_t n_param_l,
                       const std::uint32_t n_param_c,
                       const std::uint32_t n_param_m)
:
__n_samples__{n_samples},
__n_waves__{n_waves},
__n_param_a__{n_param_a},
__n_param_l__{n_param_l},
__n_param_c__{n_param_c},
__n_param_m__{n_param_m},
__trapezw_samples__{__n_samples__}
{
   
}


gms::radiolocation
::trapezoid_waveform_fp16_t
::trapezoid_waveform_fp16_t(trapezoid_waveform_fp16_t &&other) noexcept(true)
:
__n_samples__{std::move(other.__n_samples__)},
__n_waves__{std::move(other.__n_waves__)},
__n_param_a__{std::move(other.__n_param_a__)},
__n_param_l__{std::move(other.__n_param_l__)},
__n_param_c__{std::move(other.__n_param_c__)},
__n_param_m__{std::move(other.__n_param_m__)},
__trapezw_samples__{std::move(other.__trapezw_samples__)}
{
     
}

gms::radiolocation
::trapezoid_waveform_fp16_t
::trapezoid_waveform_fp16_t(const trapezoid_waveform_fp16_t &other) noexcept(false)
:
__n_samples__{other.__n_samples__},
__n_waves__{other.__n_waves__},
__n_param_a__{other.__n_param_a__},
__n_param_l__{other.__n_param_l__},
__n_param_c__{other.__n_param_c__},
__n_param_m__{other.__n_param_m__},
__trapezw_samples__{other.__trapezw_samples__}
{
     
}


gms::radiolocation
::trapezoid_waveform_fp16_t
::~trapezoid_waveform_fp16_t()
{

}

gms::radiolocation
::trapezoid_waveform_fp16_t &
gms::radiolocation
::trapezoid_waveform_fp16_t
::operator=(trapezoid_waveform_fp16_t &&other) noexcept(true)
{
    if(this==&other) { return (*this);}
    this->__n_samples__               = std::move(other.__n_samples__);
    this->__n_waves__                 = std::move(other.__n_waves__);
    this->__n_param_a__               = std::move(other.__n_param_a__);
    this->__n_param_c__               = std::move(other.__n_param_c__);
    this->__n_param_l__               = std::move(other.__n_param_l__);
    this->__n_param_m__               = std::move(other.__n_param_m__);
    this->__trapezw_samples__.operator=(std::move(other.__trapezw_samples__));
    return (*this);
}

void
gms::radiolocation 
::trapezoid_waveform_fp16_t
::init_storage(const half_float::half filler)
{
#if (INIT_BY_STD_FILL) == 0
     using namespace gms::common;
	 sse_memset_unroll8x_ps(reinterpret_cast<float*>(this->__trapezw_samples__.m_data),filler,this->__n_samples__);
#else 
     std::fill(this->__trapezw_samples__,this->__trapezw_samples__+this->__n_samples__,filler);
#endif 
}

void
gms::radiolocation
::trapezoid_waveform_fp16_t
::create_signal_plot(const std::uint32_t n_samp,
                    const half_float::half * __restrict sig_arg,
                    const half_float::half * __restrict sig_val,
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

float
gms::radiolocation
::trapezoid_waveform_fp16_t
::sample_of_trapezoid_wave(const float t,
	                       const float a,
                           const float m,
						   const float l,
					       const float c)
{
	 constexpr float PI{3.14159265358979323846264338328f};
	 float     a_over_PI{a/PI};
	 float     PI_over_m{PI/m};
	 const float arg{PI_over_m*t+l};
#if (TRAPEZOID_WAVEFORM_USE_CEPHES) == 1
         const float t_as{ceph_asinf(ceph_sinf(arg))};
		 const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
         const float t_as{std::asin(std::sin(arg))};
		 const float t_ac{std::acos(std::cos(arg))};
#endif 
         const float t_0{a_over_PI*(t_as+t_ac)};
		 const float t_1{t_0-5.0f+c};
		 return (t_1);
}

void 
gms::radiolocation
::trapezoid_waveform_fp16_t
::single_trapezoid_wave(const float a,
                        const float m,
						const float l,
					    const float c)
{
     using namespace half_float;
	 constexpr float PI{3.14159265358979323846264338328f};
	 float     a_over_PI{a/PI};
	 float     PI_over_m{PI/m};
	 const  std::uint32_t n_samples{static_cast<std::uint32_t>(this->__n_samples__)};
     for(std::uint32_t __i{0}; __i != n_samples; ++__i)
	 {
		 const float t__i{static_cast<float>(__i)};
		 const float arg{PI_over_m*t__i+l};
#if (TRAPEZOID_WAVEFORM_USE_CEPHES) == 1
         const float t_as{ceph_asinf(ceph_sinf(arg))};
		 const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
         const float t_as{std::asin(std::sin(arg))};
		 const float t_ac{std::acos(std::cos(arg))};
#endif 
         const float t_0{a_over_PI*(t_as+t_ac)};
		 const float t_1{t_0-5.0f+c};
		 this->__trapezw_samples__.m_data[__i] = half_cast<half>(t_1);
	 }
}


void 
gms::radiolocation
::trapezoid_waveform_fp16_t
::series_of_trapezoid_waves(const float a,
                            const float m,
						    const float l,
						    const float c,
							const std::uint32_t shaping)
{
     using namespace half_float;
     constexpr float PI{3.14159265358979323846264338328f};
	 float     a_over_PI{a/PI};
	 float     PI_over_m{PI/m};
	 float     sum;
	 const  std::uint32_t n_samples{static_cast<std::uint32_t>(this->__n_samples__)};
	 switch (shaping)
	 {
		   case 0 : 
		      for(std::uint32_t __i{0}; __i != n_samples; ++__i)
			  {
				   const float t__i{static_cast<float>(__i)};
				   sum = 0.0f;
				   for(std::uint32_t __j{0}; __j != this->__n_waves__; ++__j) 
				   {
					    const float arg{PI_over_m*t__i+l};
#if (TRAPEZOID_WAVEFORM_FP16_USE_CEPHES) == 1
                        const float t_as{ceph_asinf(ceph_sinf(arg))};
		                const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                        const float t_as{std::asin(std::sin(arg))};
		                const float t_ac{std::acos(std::cos(arg))};
#endif 
                        sum = sum+a_over_PI*(t_as+t_ac)-5.0f+c;
				   }
				   this->__trapezw_samples__.m_data[__i] = half_cast<half>(sum);
			  }
		  break;
		  case 1 : 
		      for(std::uint32_t __i{0}; __i != n_samples; ++__i)
			  {
				   const float t__i{static_cast<float>(__i)};
				   sum = 0.0f;
				   for(std::uint32_t __j{0}; __j != this->__n_waves__; ++__j) 
				   {
					    const float t__j = static_cast<float>(__j);
					    const float arg{PI_over_m*(t__i+t__j)+l};
#if (TRAPEZOID_WAVEFORM_USE_CEPHES) == 1
                        const float t_as{ceph_asinf(ceph_sinf(arg))};
		                const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                        const float t_as{std::asin(std::sin(arg))};
		                const float t_ac{std::acos(std::cos(arg))};
#endif 
                        sum = sum+a_over_PI*(t_as+t_ac)-5.0f+c;
				   }
				   this->__trapezw_samples__.m_data[__i] = half_cast<half>(sum);
			  }
		  break;
		  case 2 : 
		      for(std::uint32_t __i{0}; __i != n_samples; ++__i)
			  {
				   const float t__i{static_cast<float>(__i)};
				   sum = 0.0f;
				   for(std::uint32_t __j{0}; __j != this->__n_waves__; ++__j) 
				   {
					    const float t__j = static_cast<float>(__j);
					    const float arg{PI_over_m*t__i*t__j+l};
#if (TRAPEZOID_WAVEFORM_USE_CEPHES) == 1
                        const float t_as{ceph_asinf(ceph_sinf(arg))};
		                const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                        const float t_as{std::asin(std::sin(arg))};
		                const float t_ac{std::acos(std::cos(arg))};
#endif 
                        sum = sum+a_over_PI*(t_as+t_ac)-5.0f+c;
				   }
				   this->__trapezw_samples__.m_data[__i] = half_cast<half>(sum);
			break;
			default : 
		         return;
		}
		   
	 
	 }
}

void 
gms::radiolocation
::trapezoid_waveform_fp16_t
::series_of_trapezoid_waves_u2x(const float a,
                                const float m,
						        const float l,
						        const float c)
{
     using namespace half_float;
	 constexpr float PI{3.14159265358979323846264338328f};
	 float     a_over_PI{a/PI};
	 float     PI_over_m{PI/m};
	 float     sum0;
	 float     sum1;
	 float     t__i_0{-2.0f};
	 float     t__i_1{-1.0f};
	 
	 const  std::uint32_t n_samples{static_cast<std::uint32_t>(this->__n_samples__)};
	 std::uint32_t __i;
	 std::uint32_t __j;
	 for(__i = 0; __i != ROUND_TO_FOUR(n_samples,2); __i += 2)
	 {
		 t__i_0  += 2.0f;
		 sum0     = 0.0f;
		 for(std::uint32_t __j{0}; __j != this->__n_waves__; ++__j) 
		 {
					    const float t__j = static_cast<float>(__j);
					    const float arg{PI_over_m*t__i_0*t__j+l};
#if (TRAPEZOID_WAVEFORM_FP16_USE_CEPHES) == 1
                        const float t_as{ceph_asinf(ceph_sinf(arg))};
		                const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                        const float t_as{std::asin(std::sin(arg))};
		                const float t_ac{std::acos(std::cos(arg))};
#endif 
                        sum0 = sum0+a_over_PI*(t_as+t_ac)-5.0f+c;
		 }
		 this->__trapezw_samples__.m_data[__i+0] = half_cast<half>(sum0);
		 t__i_1  += 2.0f;
		 sum1    =  0.0f;
         for(std::uint32_t __j{0}; __j != this->__n_waves__; ++__j) 
		 {
					    const float t__j = static_cast<float>(__j);
					    const float arg{PI_over_m*t__i_1*t__j+l};
#if (TRAPEZOID_WAVEFORM_FP16_USE_CEPHES) == 1
                        const float t_as{ceph_asinf(ceph_sinf(arg))};
		                const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                        const float t_as{std::asin(std::sin(arg))};
		                const float t_ac{std::acos(std::cos(arg))};
#endif 
                        sum1 = sum1+a_over_PI*(t_as+t_ac)-5.0f+c;
		 }
		 this->__trapezw_samples__.m_data[__i+1] = half_cast<half>(sum1);
		
	 }

	 for(__j = __i; __j != n_samples; ++__j)
	 {
		  const float t__j{static_cast<float>(__j)};
		  sum0 = 0.0f;
		  for(std::uint32_t __k{0}; __k != this->__n_waves__; ++__k) 
		  {
					    const float t__k = static_cast<float>(__k);
					    const float arg{PI_over_m*t__j*t__k+l};
#if (TRAPEZOID_WAVEFORM_FP16_USE_CEPHES) == 1
                        const float t_as{ceph_asinf(ceph_sinf(arg))};
		                const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                        const float t_as{std::asin(std::sin(arg))};
		                const float t_ac{std::acos(std::cos(arg))};
#endif 
                        sum0 = sum0+a_over_PI*(t_as+t_ac)-5.0f+c;
		  }
		 this->__trapezw_samples__.m_data[__j] = half_cast<half>(sum0);
	 }
}

void 
gms::radiolocation
::trapezoid_waveform_fp16_t
::series_of_trapezoid_waves_u4x(const float a,
                                const float m,
						        const float l,
						        const float c)
{
     using namespace half_float;
	 constexpr float PI{3.14159265358979323846264338328f};
	 float     a_over_PI{a/PI};
	 float     PI_over_m{PI/m};
	 float     sum0;
	 float     sum1;
	 float     sum2;
	 float     sum3;
	 float     t__i_0{-4.0f};
	 float     t__i_1{1.0f};
	 float     t__i_2{2.0f};
	 float     t__i_3{3.0f};
	 const  std::uint32_t n_samples{static_cast<std::uint32_t>(this->__n_samples__)};
	 std::uint32_t __i;
	 std::uint32_t __j;
	 for(__i = 0; __i != ROUND_TO_FOUR(n_samples,4); __i += 4)
	 {
		 t__i_0  += 4.0f;
		 sum0     = 0.0f;
		 for(std::uint32_t __j{0}; __j != this->__n_waves__; ++__j) 
		 {
					    const float t__j = static_cast<float>(__j);
					    const float arg{PI_over_m*t__i_0*t__j+l};
#if (TRAPEZOID_WAVEFORM_FP16_USE_CEPHES) == 1
                        const float t_as{ceph_asinf(ceph_sinf(arg))};
		                const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                        const float t_as{std::asin(std::sin(arg))};
		                const float t_ac{std::acos(std::cos(arg))};
#endif 
                        sum0 = sum0+a_over_PI*(t_as+t_ac)-5.0f+c;
		 }
		 this->__trapezw_samples__.m_data[__i+0] = half_cast<half>(sum0);
		 t__i_1  += 4.0f;
		 sum1    =  0.0f;
         for(std::uint32_t __j{0}; __j != this->__n_waves__; ++__j) 
		 {
					    const float t__j = static_cast<float>(__j);
					    const float arg{PI_over_m*t__i_1*t__j+l};
#if (TRAPEZOID_WAVEFORM_FP16_USE_CEPHES) == 1
                        const float t_as{ceph_asinf(ceph_sinf(arg))};
		                const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                        const float t_as{std::asin(std::sin(arg))};
		                const float t_ac{std::acos(std::cos(arg))};
#endif 
                        sum1 = sum1+a_over_PI*(t_as+t_ac)-5.0f+c;
		 }
		 this->__trapezw_samples__.m_data[__i+1] = half_cast<half>(sum1);
		 t__i_2  += 4.0f;
		 sum2    =  0.0f;
		 for(std::uint32_t __j{0}; __j != this->__n_waves__; ++__j) 
		 {
					    const float t__j = static_cast<float>(__j);
					    const float arg{PI_over_m*t__i_2*t__j+l};
#if (TRAPEZOID_WAVEFORM_FP16_USE_CEPHES) == 1
                        const float t_as{ceph_asinf(ceph_sinf(arg))};
		                const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                        const float t_as{std::asin(std::sin(arg))};
		                const float t_ac{std::acos(std::cos(arg))};
#endif 
                        sum2 = sum2+a_over_PI*(t_as+t_ac)-5.0f+c;
		 }
		  this->__trapezw_samples__.m_data[__i+2] = half_cast<half>(sum2);
		  t__i_3  += 4.0f;
		  sum3    =  0.0f;
		  for(std::uint32_t __j{0}; __j != this->__n_waves__; ++__j) 
		  {
					    const float t__j = static_cast<float>(__j);
					    const float arg{PI_over_m*t__i_3*t__j+l};
#if (TRAPEZOID_WAVEFORM_FP16_USE_CEPHES) == 1
                        const float t_as{ceph_asinf(ceph_sinf(arg))};
		                const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                        const float t_as{std::asin(std::sin(arg))};
		                const float t_ac{std::acos(std::cos(arg))};
#endif 
                        sum3 = sum3+a_over_PI*(t_as+t_ac)-5.0f+c;
		 }
		  this->__trapezw_samples__.m_data[__i+3] = half_cast<half>(sum3);
	 }

	 for(__j = __i; __j != n_samples; ++__j)
	 {
		  const float t__j{static_cast<float>(__j)};
		  sum0 = 0.0f;
		  for(std::uint32_t __k{0}; __k != this->__n_waves__; ++__k) 
		  {
					    const float t__k = static_cast<float>(__k);
					    const float arg{PI_over_m*t__j*t__k+l};
#if (TRAPEZOID_WAVEFORM_FP16_USE_CEPHES) == 1
                        const float t_as{ceph_asinf(ceph_sinf(arg))};
		                const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                        const float t_as{std::asin(std::sin(arg))};
		                const float t_ac{std::acos(std::cos(arg))};
#endif 
                        sum0 = sum0+a_over_PI*(t_as+t_ac)-5.0f+c;
		  }
		 this->__trapezw_samples__.m_data[__j] = half_cast<half>(sum0);
	 }
}


void 
gms::radiolocation
::trapezoid_waveform_fp16_t
::single_trapezoid_wave_v2(const float a,
                           const float m,
						   const float l,
						   const float c,
						   const std::uint32_t shaping)
{
     using namespace half_float;
	 constexpr float PI{3.14159265358979323846264338328f};
	 //float     a_over_PI{a*0.318309886183790671537767526745f};
	 //float     PI_over_m{PI/m};
	 float __a{a};
	 float __m{m};
	 float __l{l};
	 float __c{c};
	 const  std::uint32_t n_samples{static_cast<std::uint32_t>(this->__n_samples__)};
	 float  sum{0.0f};
	 this->init_storage(half_cast<half>(0.0f));
     switch (shaping)
	 {
		case 0 : 
		      for(std::uint32_t __i{1}; __i != this->__n_waves__; ++__i) 
			  {
				  const float t__i{static_cast<float>(__i)};
				  __a += t__i;
				  __m += t__i;
				  __l += t__i;
				  __c += t__i;
                  float a_over_PI{__a*0.318309886183790671537767526745f};
				  float PI_over_m{PI/__m};
				  for(std::uint32_t __j{0}; __j != n_samples; ++__j) 
				  {
                        const float arg{PI_over_m*t__i+__l};
#if (TRAPEZOID_WAVEFORM_FP16_USE_CEPHES) == 1
                        const float t_as{ceph_asinf(ceph_sinf(arg))};
		                const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                        const float t_as{std::asin(std::sin(arg))};
		                const float t_ac{std::acos(std::cos(arg))};
#endif 
                        sum = a_over_PI*(t_as+t_ac)-5.0f+__c;
						this->__trapezw_samples__.m_data[__j] += half_cast<half>(sum);
				  }
			  }
		break;
		case 1 : 
		      for(std::uint32_t __i{1}; __i != this->__n_waves__; ++__i) 
			  {
				  const float t__i{static_cast<float>(__i)};
				  __a += t__i;
				  __m += t__i;
				  __l += t__i;
				  __c += t__i;
                  float a_over_PI{__a*0.318309886183790671537767526745f};
				  float PI_over_m{PI/__m};
				  for(std::uint32_t __j{0}; __j != n_samples; ++__j) 
				  {
					    const float t__j{static_cast<float>(__j)};
                        const float arg{PI_over_m*(t__i+t__j)+__l};
#if (TRAPEZOID_WAVEFORM_FP16_USE_CEPHES) == 1
                        const float t_as{ceph_asinf(ceph_sinf(arg))};
		                const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                        const float t_as{std::asin(std::sin(arg))};
		                const float t_ac{std::acos(std::cos(arg))};
#endif 
                        sum = a_over_PI*(t_as+t_ac)-5.0f+__c;
						this->__trapezw_samples__.m_data[__j] += half_cast<half>(sum);
				  }
			  }
		break;
		case 2 : 
		      for(std::uint32_t __i{1}; __i != this->__n_waves__; ++__i) 
			  {
				  const float t__i{static_cast<float>(__i)};
				  __a += t__i;
				  __m += t__i;
				  __l += t__i;
				  __c += t__i;
                  float a_over_PI{__a*0.318309886183790671537767526745f};
				  float PI_over_m{PI/__m};
				  for(std::uint32_t __j{0}; __j != n_samples; ++__j) 
				  {
					    const float t__j{static_cast<float>(__j)};
                        const float arg{PI_over_m*t__i*t__j+__l};
#if (TRAPEZOID_WAVEFORM_FP16_USE_CEPHES) == 1
                        const float t_as{ceph_asinf(ceph_sinf(arg))};
		                const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                        const float t_as{std::asin(std::sin(arg))};
		                const float t_ac{std::acos(std::cos(arg))};
#endif 
                        sum = a_over_PI*(t_as+t_ac)-5.0f+__c;
						this->__trapezw_samples__.m_data[__j] += half_cast<half>(sum);
				  }
			  }
		break;
		default : 
		      return;
	 }
}

void 
gms::radiolocation 
::trapezoid_waveform_fp16_t
::single_trapezoid_wave_v2(const float * __restrict p_a,
                           const float * __restrict p_m,
						   const float * __restrict p_l,
						   const float * __restrict p_c,
						   const std::uint32_t shaping)
{
    using namespace half_float;
	constexpr float PI{3.14159265358979323846264338328f};
	 //float     a_over_PI{a*0.318309886183790671537767526745f};
	 //float     PI_over_m{PI/m};
	 const float * __restrict __p_a{p_a};
	 const float * __restrict __p_m{p_m};
	 const float * __restrict __p_l{p_l};
	 const float * __restrict __p_c{p_c};
	 const  std::uint32_t n_samples{static_cast<std::uint32_t>(this->__n_samples__)};
	 float  sum{0.0f};
	 this->init_storage(half_cast<half>(0.0f));
     switch (shaping)
	 {
		case 0 : 
		      for(std::uint32_t __i{0}; __i != this->__n_waves__; ++__i) 
			  {
				  const float t__i{static_cast<float>(__i)};
				  const float __a{__p_a[__i]}; 
				  const float __m{__p_m[__i]};
				  const float __l{__p_l[__i]};
				  const float __c{__p_c[__i]};
                  float a_over_PI{__a*0.318309886183790671537767526745f};
				  float PI_over_m{PI/__m};
				  for(std::uint32_t __j{0}; __j != n_samples; ++__j) 
				  {
                        const float arg{PI_over_m*t__i+__l};
#if (TRAPEZOID_WAVEFORM_FP16_USE_CEPHES) == 1
                        const float t_as{ceph_asinf(ceph_sinf(arg))};
		                const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                        const float t_as{std::asin(std::sin(arg))};
		                const float t_ac{std::acos(std::cos(arg))};
#endif 
                        sum = a_over_PI*(t_as+t_ac)-5.0f+__c;
						this->__trapezw_samples__.m_data[__j] += half_cast<half>(sum);
				  }
			  }
		break;
		case 1 : 
		      for(std::uint32_t __i{0}; __i != this->__n_waves__; ++__i) 
			  {
				  const float t__i{static_cast<float>(__i)};
				  const float __a{__p_a[__i]}; 
				  const float __m{__p_m[__i]};
				  const float __l{__p_l[__i]};
				  const float __c{__p_c[__i]};
                  float a_over_PI{__a*0.318309886183790671537767526745f};
				  float PI_over_m{PI/__m};
				  for(std::uint32_t __j{0}; __j != n_samples; ++__j) 
				  {
					    const float t__j{static_cast<float>(__j)};
                        const float arg{PI_over_m*(t__i+t__j)+__l};
#if (TRAPEZOID_WAVEFORM_FP16_USE_CEPHES) == 1
                        const float t_as{ceph_asinf(ceph_sinf(arg))};
		                const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                        const float t_as{std::asin(std::sin(arg))};
		                const float t_ac{std::acos(std::cos(arg))};
#endif 
                        sum = a_over_PI*(t_as+t_ac)-5.0f+__c;
						this->__trapezw_samples__.m_data[__j] += half_cast<half>(sum);
				  }
			  }
		break;
		case 2 : 
		      for(std::uint32_t __i{0}; __i != this->__n_waves__; ++__i) 
			  {
				  const float t__i{static_cast<float>(__i)};
				  const float __a{__p_a[__i]}; 
				  const float __m{__p_m[__i]};
				  const float __l{__p_l[__i]};
				  const float __c{__p_c[__i]};
                  float a_over_PI{__a*0.318309886183790671537767526745f};
				  float PI_over_m{PI/__m};
				  for(std::uint32_t __j{0}; __j != n_samples; ++__j) 
				  {
					    const float t__j{static_cast<float>(__j)};
                        const float arg{PI_over_m*t__i*t__j+__l};
#if (TRAPEZOID_WAVEFORM_FP16_USE_CEPHES) == 1
                        const float t_as{ceph_asinf(ceph_sinf(arg))};
		                const float t_ac{ceph_acosf(ceph_cosf(arg))};
#else 
                        const float t_as{std::asin(std::sin(arg))};
		                const float t_ac{std::acos(std::cos(arg))};
#endif 
                        sum = a_over_PI*(t_as+t_ac)-5.0f+__c;
						this->__trapezw_samples__.m_data[__j] += half_cast<half>(sum);
				  }
			  }
		break;
		default : 
		      return;
	 }
}
