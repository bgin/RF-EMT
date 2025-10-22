/*MIT License
!Copyright (c) 2020 Bernard Gingold
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
*/

#ifndef __GMS_ATMOS_REFRACTION_H__
#define __GMS_ATMOS_REFRACTION_H__ 


namespace file_info 
{

     static const unsigned int GMS_ATMOS_REFRACTION_MAJOR = 1;
     static const unsigned int GMS_ATMOS_REFRACTION_MINOR = 0;
     static const unsigned int GMS_ATMOS_REFRACTION_MICRO = 0;
     static const unsigned int GMS_ATMOS_REFRACTION_FULLVER =
       1000U*GMS_ATMOS_REFRACTION_MAJOR+100U*GMS_ATMOS_REFRACTION_MINOR+
       10U*GMS_ATMOS_REFRACTION_MICRO;
     static const char GMS_ATMOS_REFRACTION_CREATION_DATE[] = "22-10-2025 01:11PM +00200 (WED 22 OCT 2025 GMT+2)";
     static const char GMS_ATMOS_REFRACTION_BUILD_DATE[]    = __DATE__; 
     static const char GMS_ATMOS_REFRACTION_BUILD_TIME[]    = __TIME__;
     static const char GMS_ATMOS_REFRACTION_SYNOPSIS[]      = "Calculation of EM Wave atmospheric refraction.
                                                               Various characteristics and formulae of atmospheric refraction 
                                                               (radio waves and visible light/IR wavelengths)  
                                                               Based mainly on Колосов М.А., Шабельников А.В. - 
                                                               Рефракция электромагнитных волн в атмосферах Земли, 
                                                               Венеры и Марса-Советское Радио (1976)";

}

#include <cstdint>
#include <string>
#include <iostream>
#include "GMS_config.h"
#include "GMS_dyn_array.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (ATMOS_REFRACTION_USE_NT_STORES)
#define ATMOS_REFRACTION_USE_NT_STORES 0
#endif

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(ATMOS_REFRACTION_INIT_STORAGE)
#define ATMOS_REFRACTION_INIT_STORAGE 1
#endif 
#endif

#if (ATMOS_REFRACTION_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(ATMOS_REFRACTION_USE_CEPHES)
#define ATMOS_REFRACTION_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (ATMOS_REFRACTION_USE_PMC_INSTRUMENTATION)
#define ATMOS_REFRACTION_USE_PMC_INSTRUMENTATION 1
#endif 

#if (ATMOS_REFRACTION_USE_PMC_INSTRUMENTATION) == 1
#include "GMS_hw_perf_macros.h"

#define PMC_VARS                      \
uint64_t prog_counters_start[4] = {}; \
uint64_t prog_counters_stop[4]  = {}; \
uint64_t tsc_start,tsc_stop;          \
uint64_t act_cyc_start,act_cyc_stop;  \
uint64_t ref_cyc_start,ref_cyc_stop;  \
[[maybe_unused]] uint64_t dummy1;     \
[[maybe_unused]] uint64_t dummy2;     \
[[maybe_unused]] uint64_t dummy3;     \
int32_t core_counter_width;           \
double utilization,nom_ghz,avg_ghz;
#endif 

#if (ATMOS_REFRACTION_USE_CEPHES) == 0
#include <cmath>
#else 
#include "GMS_cephes_sin_cos.h"
#endif 

namespace gms
{

namespace atmos
{
          
          struct alignas(64) atmos_refraction_size_params_t 
          {
                  std::size_t  n_nidx;     //!size of refractive index source
                  std::size_t  n_n0idx;    //!size of refractive index dest
                  std::size_t  n_z;        //!size of 'z' angle (dest)
                  std::size_t  n_z0;       //!size of 'z0' angle (source)
                  std::size_t  n_r;        //!upper limit of integration (radius)
                  std::size_t  n_R0;       //!lower limit of integration (radius)
                  std::size_t  n_phi;      //!'phi' angle at dest
                  std::size_t  n_phi0;     //!'phi' angle at source 
                  std::size_t  n_ntht;     //!n-index over theta (result)
                  std::size_t  n_nphi;     //!n-index over phi (result)
                  std::size_t  n_dndr;     //!n-index derivative over r
                  std::size_t  n_rho;      //!result of rad_ray_curvature_f251_r4
                  std::size_t  n_dn0;      //!coefficient of refreaction near the Earth surface i.e. dn0 = (240*10e-6->380*10e-6)
                  std::size_t  n_beta;     //!coefficient describing the diminishing of 'n' as function of height, i.e. 0.10->0.14 1/km
                  std::size_t  n_navgh;    //!result of n_avg_h_f145_r4
                  std::size_t  n_H;        //! height
                  std::size_t  n_f345;     //! result of refraction_angle_f345_r4
                  std::size_t  n_f351;     //! result of refraction_angle_n90_f351_r4
                  std::size_t  n_f352;     //! result of refraction_angle_at90_f352_r4
                  std::size_t  n_nh;       //! n-index over height 
                  std::size_t  n_f41;      //! result of refraction_angle_for_gl5cm_f41_r4
                  std::size_t  n_f;        //! central signal frequency 
                  std::size_t  n_d;        //! height at maximum of layer F2 
                  std::size_t  n_Nmf;      //! electron density in layer F2 
                  std::size_t  n_f412;     //! result of refractive_idx_lo_ionosphere_f412_r4
                  std::size_t  n_f413;     //! result of refractive_idx_hi_ionosphere_f413_r4
                  std::size_t  n_D1; 
                  std::size_t  n_f415;     //!result of refraction_angle_in_ionosphere_f415_r4
                  std::size_t  n_f425;     //!result of refraction_angle_ionosphere_z0le60_f425_r4
                  std::size_t  n_f428;     //!result of refraction_angle_ionosphere_z0eq90_f428_r4
                  std::size_t  n_f429;    //!result of n_avg_0_h_H1_f429_r4
                  std::size_t  n_H1;       //! H1<=h<=H2
                  std::size_t  n_H2;      //!  H1<=h<=H2
                  std::size_t  n_f430;     //! n_avg_H1_h_H2_f430_r4
                  std::size_t  n_H3;       //! H2<=h<=H3 
                  std::size_t  n_deln0; 
                  std::size_t  n_f431;     //! result of n_avg_H2_h_H3_f431_r4 
                  std::size_t  n_f438;     //! result of refraction_angle_atmos_L1_lo_f438_r4
                  std::size_t  n_f442;     //! result of refraction_angle_atmos_L2_med_f442_r4 
                  std::size_t  n_f445;     //! result of refraction_angle_atmos_L3_upper_f445_r4
                  std::size_t  n_g;        //! formula: 4.50, page: 80  
                  std::size_t  n_f450;     //! result of refraction_angle_z0le60_med_atmos_f450_r4
                  std::size_t  n_f451;     //! result of refraction_angle_z0eq90_med_atmos_f451_r4
                  std::size_t  n_Hc0;      //! formula 5.4, page: 93 
                  std::size_t  n_delnA;    //! formula 5.4, page: 93 
                  std::size_t  n_na;       //! formula 5.4, page: 93 
                  std::size_t  n_nc;       //! formula 5.4, page: 93 
                  std::size_t  n_f53;      //! result of refraction_angle_tropo_wvle5cm_f53_r4
                  std::size_t  n_f517;     //! result of refraction_angle_tropo_wvle5cm_f517_r4
                  std::size_t  n_tht;      //! theta angles 
                  std::size_t  n_f531;     //! result of ray_traj_inhomogenous_atmos_f531_r4
                  std::size_t  n_thtc;     //! Formula: 5.34, page: 100
                  std::size_t  n_f534;     //! result of analytic_sol_L_atmos_wvle5cm_f534_r4
                  std::size_t  n_Rc;       //! !Formula: 5.35, page: 101
                  std::size_t  n_f535;     //! result of refraction_angle_atmos_wvle5cm_f535_r4
                  std::size_t  n_f538;     //! result of refraction_angle_astronomical_wvle5cm_f538_r4
                  std::size_t  n_f539;     //! result of refraction_angle_astronomical_wvle5cm_z0eq90_f539_r4
                  std::size_t  n_f541;     //! result of refractive_idx_lo_ionosphere_wv5cm3m_f541_r4
                  std::size_t  n_H0;       //! described in: analytic_sol_L2_whole_atmosphere_f552_r4
                  std::size_t  n_Hc;       //! ! Formula: 5.49, page: 105
                  std::size_t  n_f543;     //! result of refraction_angle_whole_atmos_vw5cm3m_f543_r4
                  std::size_t  n_H10;      //! described in ! Formula: 5.65, page: 108
                  std::size_t  n_f554;     //! result of  refraction_angle_whole_atmos_vw5cm3m_f554_r4
                  std::size_t  n_H20;     //! described in ! Formula: 5.72, page: 110
                  std::size_t  n_f572;     //! result of refractive_idx_hi_ionosphere_approx_f572_r4
                  std::size_t  n_f576;     //! result of analytic_sol_L31_whole_atmos_wv5cm3m_f576_r4
                  std::size_t  n_f577;     //! result of analytic_sol_L32_whole_atmos_wv5cm3m_f577_r4
                  std::size_t  n_f578;     //! result of analytic_sol_L33_whole_atmos_wv5cm3m_f578_r4
                  std::size_t  n_f579;     //! result of analytic_sol_L34_whole_atmos_wv5cm3m_f579_r4
                  std::size_t  n_f590;     //! result of refraction_angle_lo_ionospere_wv5cm3m_f590_r4
                  std::size_t  n_R2;       //! described in  !Formula: 5.91, page: 115
                  std::size_t  n_f591;     //! result of refraction_angle_lo_ionospere_wv5cm3m_f591_r4
                  std::size_t  n_f593;     //! result of refraction_angle_lo_ionospere_wv5cm3m_f593_r4
                  std::size_t  n_f595;     //! refraction_angle_lo_ionospere_wv5cm3m_f595_r4
                  std::size_t  n_HB;       //! described in  !Formula: 6.3, page: 119
                  std::size_t  n_f62;      //! result of refraction_angle_B_whole_atmos_vwl5cm_f62_r4
                  std::size_t  n_HC2;      // ! described in   !Formula: 6.7, page: 120
                  std::size_t  n_f66;     // ! result of refraction_angle_C_whole_atmos_vwl5cm_f66_r4
                  std::size_t  n_f61;     // ! result of refraction_angle_whole_atmos_vwl5cm_f61_r4
                  std::size_t  n_f61b;    // ! result of refraction_angle_whole_atmos_vwl5cm_f61b_r4
                  std::size_t  n_Bf61b;   // ! result of refraction_angle_whole_atmos_case_B_vwl5cm_f61b_r4
                  std::size_t  n_f619;    // ! result of deriv_alpha_over_R0_f619_r4
                  std::size_t  n_Lc;      // ! described in !Formula: 6.18, page: 122
                  std::size_t  n_Lb;      // ! described in !Formula: 6.18, page: 122
                  std::size_t  n_gamma;   // ! described in  !Formula: 6.18, page: 122
                  std::size_t  n_f618;    // ! result of refracted_signal_weakening_Vp_f618_r4
                  std::size_t  n_f620;    // ! result of refracted_signal_weakening_case_1_Vp_f620_r4
                  std::size_t  n_f621;     //! result of refracted_signal_weakening_case_2_Vp_f621_r4
                  std::size_t  n_Lh;      // ! described in !Formula: 6.22, page: 123
                  std::size_t  n_f622;    // ! result of refracted_signal_weakening_case_3_Vp_f622_r4
                  std::size_t  n_f623;    // ! result of refraction_angle_C_earth_atmos_case_1_wv5cm3m_f623_r4
                  std::size_t  n_f625;    // ! result of refraction_angle_C_earth_atmos_case_2_wv5cm3m_f625_r4
                  std::size_t  n_f627;    // ! result of refraction_angle_C_earth_atmos_case_3_wv5cm3m_f627_r4
                  std::size_t  n_Hh;      // ! described in !Formula: 7.2, page: 132
                  std::size_t  n_f72;     // ! result of refraction_angle_C_earth_atmos_stratified_f72_r4
                  std::size_t  n_delnb;   // ! described in  !Formula: 7.4, page: 132
                  std::size_t  n_delnc;   // ! described in  !Formula: 7.4, page: 132
                  std::size_t  n_delnh;   // ! described in  !Formula: 7.4, page: 132
                  std::size_t  n_f74;     // ! result of refraction_angle_C_earth_atmos_stratified_case_1_f74_r4
                  std::size_t  n_f75;     // ! result of refraction_angle_C_earth_atmos_stratified_case_2_f75_r4
                  std::size_t  n_f714;    // ! result of refraction_angle_delta_ionosphere_strata_f714_r4
                  std::size_t  n_1f714;   // ! result of refraction_angle_delta_ionosphere_strata_case_1_f714_r4
                  std::size_t  n_f739;    // ! result of a_gamma_coeff_f739_r4
                  std::size_t  n_f741;    // ! result of refraction_angle_atmos_2D_stratified_f741_r4
                  std::size_t  n_f743;    // ! result of refraction_angle_atmos_2D_stratified_f743_r4
                  std::size_t  n_f744;     //! result of refraction_angle_atmos_2D_stratified_f744_r4
                  std::size_t  n_f747;    // ! result of refraction_angle_atmos_2D_stratified_f747_r8
                  std::size_t  n_f96;     // ! result of analytic_sol_phase_to_geo_path_f96_r4
                  std::size_t  n_f910;    // ! result of analytic_sol_phase_to_geo_path_case_1_f910_r4
                  std::size_t  n_2f96;    // ! result of analytic_sol_phase_to_geo_path_case_2_f96_r4
                  std::size_t  n_f914;    // ! result of  analytic_sol_phase_shift_ionosphere_to_earth_f914_r4
                  std::size_t  n_L;       // ! described in  !Formula: 9.24, page: 157
                  std::size_t  n_f924;     //! result of emitter_height_delta_atmos_refraction_f924_r4 
          };


} // atmos

} //gms
























#endif /*__GMS_ATMOS_REFRACTION_H__*/