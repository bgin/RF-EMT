#include <fstream>
#include <iomanip>
#include "GMS_am_bb_trapez_signal_f64.h"
#include "GMS_sse_memset.h"
#include "GMS_indices.h"
#include "GMS_compare_fp_safe.h"

gms::radiolocation 
::am_bb_trapez_signal_f64_t
::am_bb_trapez_signal_f64_t(const std::size_t nsamples,
                        const std::uint32_t nK,
                        const double a,
                        const double l,
                        const double c,
                        const double m)
:
m_nsamples{nsamples},
m_nK{nK},
m_a{a},
m_l{l},
m_c{c},
m_m{m},
m_sig_samples{darray_r8_t(m_nsamples)}
{

}

gms::radiolocation 
::am_bb_trapez_signal_f64_t
::am_bb_trapez_signal_f64_t(am_bb_trapez_signal_f64_t &&other)
:
m_nsamples{std::move(other.m_nsamples)},
m_nK{std::move(other.m_nK)},
m_a{std::move(other.m_a)},
m_l{std::move(other.m_l)},
m_c{std::move(other.m_c)},
m_m{std::move(other.m_m)},
m_sig_samples{std::move(other.m_sig_samples)}
{

}

gms::radiolocation 
::am_bb_trapez_signal_f64_t
::~am_bb_trapez_signal_f64_t()
{

}

gms::radiolocation 
::am_bb_trapez_signal_f64_t & 
gms::radiolocation 
::am_bb_trapez_signal_f64_t
::operator=(am_bb_trapez_signal_f64_t &&other)
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_nsamples            = std::move(other.m_nsamples);
    this->m_nK                  = std::move(other.m_nK);
    this->m_a                   = std::move(other.m_a);
    this->m_l                   = std::move(other.m_l);
    this->m_c                   = std::move(other.m_c);
    this->m_m                   = std::move(other.m_m);
    this->m_sig_samples.operator=(std::move(other.m_sig_samples));
    return (*this);
}

void 
gms::radiolocation
::am_bb_trapez_signal_f64_t
::init_storage(const double filler)
{         
#if (INIT_BY_STD_FILL) == 0
     using namespace gms::common;
	 sse_memset_unroll8x_pd(&this->m_sig_samples.m_data[0],filler,this->m_nsamples);
#else 
     std::fill(this->m_sig_samples.m_data,this->m_sig_samples.m_data+this->m_nsamples,filler);
#endif 
}

void 
gms::radiolocation 
::am_bb_trapez_signal_f64_t
::create_signal_plot(const std::uint32_t n_samp,
                    const double * __restrict sig_arg,
                    const double * __restrict sig_val,
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
::am_bb_trapez_signal_f64_t
::create_signal_user_data(const double * __restrict__ sym_in, // size of m_nsamples*m_nK
                          const std::uint32_t n_T,
                          const std::uint32_t n_K)
{
      if(__builtin_expect(static_cast<std::uint32_t>(this->m_nsamples)!=n_T,0) || 
         __builtin_expect(this->m_nK!=n_K,0)) { return (-1);}

      const double T{static_cast<double>(this->m_nsamples)};
      double sum;
      for(std::uint32_t __t{0ull}; __t != n_T; ++__t) 
      {
            const double t{static_cast<double>(__t)};
            sum = 0.0;
            for(std::uint32_t __k{0}; __k != n_K; ++__k) 
            {
                const double k{static_cast<double>(__k)};
                const double arg{t-k*T};
                sum += trapezoid_sample(arg)*sym_in[Ix2D(__t,n_K,__k)];
            }
            this->m_sig_samples.m_data[__t] = sum;
      }
      return (0);
}

std::int32_t 
gms::radiolocation
::am_bb_trapez_signal_f64_t
::create_signal_user_data_u4x(const double * __restrict__ sym_in, // size of __I_n_samples__ * __I_n_K__
                              const std::uint32_t n_T,
                              const std::uint32_t n_K)
{
      if(__builtin_expect(static_cast<std::uint32_t>(this->m_nsamples)!=n_T,0) || 
         __builtin_expect(this->m_nK!=n_K,0)) { return (-2);}

      const double T{static_cast<double>(this->m_nsamples)};
      double sum0;
      double sum1;
      double sum2;
      double sum3;
      double t__i_0{-4.0};
      double t__i_1{1.0};
      double t__i_2{2.0};
      double t__i_3{3.0};
      std::uint32_t __i;
      std::uint32_t __j;
      for(__i = 0; __i != ROUND_TO_FOUR(n_T,3); __i += 4) 
      {
           t__i_0 += 4.0;
           sum0   = 0.0;
           for(std::uint32_t __k{0}; __k != n_K; ++__k) 
           {
                const double k{static_cast<double>(__k)};
                const double arg{t__i_0-k*T};
                sum0 += trapezoid_sample(arg)*sym_in[Ix2D(__i+0,n_K,__k)];
           }
           this->m_sig_samples.m_data[__i+0] = sum0;
           t__i_1 += 4.0;
           sum1   = 0.0;
           for(std::uint32_t __k{0}; __k != n_K; ++__k) 
           {
                const double k{static_cast<double>(__k)};
                const double arg{t__i_1-k*T};
                sum1 += trapezoid_sample(arg)*sym_in[Ix2D(__i+1,n_K,__k)];
           }
           this->m_sig_samples.m_data[__i+1] = sum1;
           t__i_2 += 4.0;
           sum2    = 0.0;
           for(std::uint32_t __k{0}; __k != n_K; ++__k) 
           {
                const double k{static_cast<double>(__k)};
                const double arg{t__i_2-k*T};
                sum2 += trapezoid_sample(arg)*sym_in[Ix2D(__i+2,n_K,__k)];
           }
           this->m_sig_samples.m_data[__i+2] = sum2;
           t__i_3 += 4.0;
           sum3   =  0.0;
           for(std::uint32_t __k{0}; __k != n_K; ++__k) 
           {
                const double k{static_cast<double>(__k)};
                const double arg{t__i_3-k*T};
                sum3 += trapezoid_sample(arg)*sym_in[Ix2D(__i+3,n_K,__k)];
           }
           this->m_sig_samples.m_data[__i+3] = sum3;
       }
       for(__j = __i; __j != n_T; ++__j) 
       {
            const double t{static_cast<double>(__j)};
            sum0 = 0.0;
            for(std::uint32_t __k{0}; __k != n_K; ++__k) 
            {
                const double k{static_cast<double>(__k)};
                const double arg{t-k*T};
                sum0 += trapezoid_sample(arg)*sym_in[Ix2D(__j,n_K,__k)];
            }
            this->m_sig_samples.m_data[__j] = sum0;
        }
        return (0);
}

auto 
gms::radiolocation
::am_bb_trapez_signal_f64_t
::operator==(const am_bb_trapez_signal_f64_t &other)->std::vector<bool>                                                       
{
    if(__builtin_expect(this==&other,0)|| 
       __builtin_expect(this->m_nsamples!=other.m_nsamples,0)) { return std::vector<bool>();}
    
    std::vector<bool> eq_comp_holder;
    double * __restrict__ p_other = &other.m_sig_samples.m_data[0];
    constexpr double eps = std::numeric_limits<double>::epsilon();
    double x;
    double y;
    for(std::size_t i = 0ull; i<this->m_nsamples; ++i) 
    {
        x = this->m_sig_samples.m_data[i];
        y = p_other[i];
        const bool xy_eq = gms::approximatelyEqual(x,y,eps);
        eq_comp_holder.operator[](i) = xy_eq;
    }
    return (eq_comp_holder);
}

auto 
gms::radiolocation
::am_bb_trapez_signal_f64_t
::operator!=(const am_bb_trapez_signal_f64_t &other)->std::vector<bool>                                         
{   // Just call the operator== and inspect the boolean values.
    return (this->operator==(other));
}

auto 
gms::radiolocation
::am_bb_trapez_signal_f64_t
::operator>=(const am_bb_trapez_signal_f64_t &other)->std::vector<bool>                                               
{
    std::vector<bool> gt_comp_holder;
    double * __restrict__ p_other = &other.m_sig_samples.m_data[0];
    constexpr double eps = std::numeric_limits<double>::epsilon();
    double x;
    double y;
    for(std::size_t i = 0ull; i<this->m_nsamples; ++i) 
    {
        x = this->m_sig_samples.m_data[i];
        y = p_other[i];
        const bool xy_gt = gms::definitelyGreaterThan(x,y,eps);
        gt_comp_holder.operator[](i) = xy_gt;
    }
    return (gt_comp_holder);
}

auto 
gms::radiolocation
::am_bb_trapez_signal_f64_t
::operator<=(const am_bb_trapez_signal_f64_t &other)->std::vector<bool>                                                             
{
    std::vector<bool> lt_comp_holder;
    double * __restrict__ p_other = &other.m_sig_samples.m_data[0];
    constexpr double eps = std::numeric_limits<double>::epsilon();
    double x;
    double y;
    for(std::size_t i = 0ull; i<this->m_nsamples; ++i) 
    {
        x = this->m_sig_samples.m_data[i];
        y = p_other[i];
        const bool xy_lt = gms::definitelyLessThan(x,y,eps);
        lt_comp_holder.operator[](i) = xy_lt;
    }
    return (lt_comp_holder);
}

auto
gms::radiolocation
::operator<<(std::ostream &os,
             gms::radiolocation::am_bb_trapez_signal_f64_t &rhs)->std::ostream &
{
    std::cout << typeid(rhs).name() << "Begin: object state dump." << std::endl;
    std::cout << "m_nsamples : "      << rhs.m_nsamples << std::endl;
    std::cout << "m_nK       : "      << rhs.m_nK       << std::endl;
    std::cout << "m_a        : "      << std::fixed << std::setprecision(17) << rhs.m_a << std::endl;
    std::cout << "m_m        : "      << std::fixed << std::setprecision(17) << rhs.m_m << std::endl;
    std::cout << "m_l        : "      << std::fixed << std::setprecision(17) << rhs.m_l << std::endl;
    std::cout << "m_c        : "      << std::fixed << std::setprecision(17) << rhs.m_c << std::endl;
    std::cout << "Signal-samples:" << std::endl;
    for(std::size_t __i{0ull}; __i != rhs.m_nsamples; ++__i)
    {
         os << std::fixed << std::setprecision(17) << 
                             rhs.m_sig_samples.m_data[__i] << std::endl;
    }
    std::cout << typeid(rhs).name() << "End: object state dump." << std::endl;
    return (os);
}
