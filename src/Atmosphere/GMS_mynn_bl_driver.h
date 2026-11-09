/* Copyright (C) Bernard Gingold, 2020-2026 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/


#ifndef __GMS_MYNN_BL_DRIVER_H__
#define __GMS_MYNN_BL_DRIVER_H__ 071120260733

#include <cstdint>
#include <optional>
#include "GMS_config.h"
#include "GMS_dyn_array.h"


namespace file_info 
{

     const unsigned int GMS_MYNN_BL_DRIVER_MAJOR = 1;
     const unsigned int GMS_MYNN_BL_DRIVER_MINOR = 0;
     const unsigned int GMS_MYNN_BL_DRIVER_MICRO = 0;
     const unsigned int GMS_MYNN_BL_DRIVER_FULLVER =
       1000U*GMS_MYNN_BL_DRIVER_MAJOR+100U*GMS_MYNN_BL_DRIVER_MINOR+
       10U*GMS_MYNN_BL_DRIVER_MICRO;
     const char * const GMS_MYNN_BL_DRIVER_CREATION_DATE = "07-11-2026 10:03 AM +00200 (SAT 07 NOV 2026 GMT+2)";
     const char * const GMS_MYNN_BL_DRIVER_BUILD_DATE    = __DATE__ " " __TIME__;
     static const char  GMS_MYNN_BL_DRIVER_AUTHOR[]      = "Programmer: Bernard Gingold, beniekg@gmail.com"; 
     const char * const GMS_MYNN_BL_DRIVER_SYNOPSIS      = "C++ wrapper for the Fortran 90 mynn_bl_driver subroutine.";

}

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (MYNN_BL_DRIVER_USE_NT_STORES)
#define MYNN_BL_DRIVER_USE_NT_STORES 0
#endif

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(MYNN_BL_DRIVER_INIT_STORAGE)
#define MYNN_BL_DRIVER_INIT_STORAGE 1
#endif 
#endif

#define WRF_CHEM 0

#if (MYNN_BL_DRIVER_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (MYNN_BL_DRIVER_USE_PMC_INSTRUMENTATION)
#define MYNN_BL_DRIVER_USE_PMC_INSTRUMENTATION 1
#endif 

#if (MYNN_BL_DRIVER_USE_PMC_INSTRUMENTATION) == 1
#include "GMS_hw_perf_macros.h"

#define PMC_VARS\
    uint64_t prog_counters_start[4] = {};\
    uint64_t prog_counters_stop[4]  = {};\
    uint64_t tsc_start,tsc_stop;\
    uint64_t act_cyc_start,act_cyc_stop;\
    uint64_t ref_cyc_start,ref_cyc_stop;\
    [[maybe_unused]] uint64_t dummy1;\
    [[maybe_unused]] uint64_t dummy2;\
    [[maybe_unused]] uint64_t dummy3;\
    int32_t core_counter_width;\
    double utilization,nom_ghz,avg_ghz;
#endif 

namespace gms 
{

namespace atmosphere
{

/*
An example of using std::optional for holding the "optional" member.
#include <valarray>
#include <optional>

#include <optional>
#include <valarray>

class DataProcessor {
public:
    std::optional<std::valarray<double>> weights;
    
    DataProcessor() noexcept(true)
      : 
    weights() {}
    // Passes std::in_place, then the size (5) and default value (1.0) for the valarray
    DataProcessor(size_t size) 
        : weights(std::in_place, size) {}
};

int main()
{
// Usage: creates a valarray of 5 elements.
DataProcessor processor(5); 
double * __restrict__ ptr = &processor.weights->operator[](0);
ptr[0] = 3.14;
return 0;
}
*/

struct alignas(64) mynn_bl_driver_args_holder_t
{
    std::int32_t m_initflag;
    std::int32_t m_grav_settling;
		std::int32_t m_bl_mynn_tkebudget;
		std::int32_t m_bl_mynn_cloudpdf;
		std::int32_t m_bl_mixlength;
		std::int32_t m_bl_mynn_edmf;
		std::int32_t m_bl_mynn_tkeadvect;
		std::int32_t m_bl_mynn_edmf_mom;
		std::int32_t m_bl_mynn_edmf_tke;
		std::int32_t m_bl_mynn_edmf_part;
		std::int32_t m_bl_mynn_cloudmix;
		std::int32_t m_bl_mynn_mixqt;
		std::int32_t m_icloud_bl;
		std::int32_t m_FLAG_QI;
		std::int32_t m_FLAG_QNI;
		std::int32_t m_FLAG_QC;
		std::int32_t m_FLAG_QNC;
		std::int32_t m_IDS;
		std::int32_t m_IDE;
		std::int32_t m_JDS;
		std::int32_t m_JDE;
		std::int32_t m_KDS;
		std::int32_t m_KDE;
		std::int32_t m_IMS;
		std::int32_t m_IME;
		std::int32_t m_JMS;
		std::int32_t m_JME;
		std::int32_t m_KMS;
		std::int32_t m_KME;
		std::int32_t m_ITS;
		std::int32_t m_ITE;
		std::int32_t m_JTS;
		std::int32_t m_JTE;
		std::int32_t m_KTS;
		std::int32_t m_KTE;
#if WRF_CHEM == 1
		std::int32_t m_nchem;
		std::int32_t m_kdvel;
		std::int32_t m_ndvel;
		std::int32_t m_num_vert_mix;
#endif
		// Scalar real variables.
		float m_delt;
		float m_dx;

    std::size_t size3d_col_major() const 
    {
        std::int32_t total_size = (m_JME-m_JMS)*(m_KME-m_KMS)*(m_IME-m_IMS);
        return (static_cast<std::size_t>(total_size));
    }

    std::size_t size3d_row_major() const 
    {
        std::int32_t total_size = (m_IME-m_IMS)*(m_KME-m_KMS)*(m_JME-m_JMS);
        return (static_cast<std::size_t>(total_size));
    }

     std::size_t size2d_col_major() const 
    {
        std::int32_t total_size = (m_JME-m_JMS)*(m_IME-m_IMS);
        return (static_cast<std::size_t>(total_size));
    }

    std::size_t size2d_row_major() const 
    {
        std::int32_t total_size = (m_IME-m_IMS)*(m_JME-m_JMS);
        return (static_cast<std::size_t>(total_size));
    }

#if WRF_CHEM == 1
    std::size_t size4D_row_major_type1() const 
    {
        std::int32_t total_size = ((m_IME-m_IMS)*(m_KME-m_KMS)*(m_JME-m_JMS)*m_nchem);
        return (static_cast<std::size_t>(total_size));
    }

    std::size_t size4D_col_major_type1() const 
    {
        std::int32_t total_size = ((m_JME-m_JMS)*(m_KME-m_KMS)*(m_IME-m_IMS)*m_nchem);
        return (static_cast<std::size_t>(total_size));
    }

    std::size_t size4D_row_major_type2() const 
    {  
        const std::int32_t dim_I = m_IME-m_IMS;
        const std::int32_t dim_J = m_JME-m_JMS;
        std::int32_t total_size = (dim_I*m_kdvel*dim_J*m_ndvel);
        return (static_cast<std::size_t>(total_size));
    }

    std::size_t size4D_col_major_type2() const 
    {  
        const std::int32_t dim_I = m_IME-m_IMS;
        const std::int32_t dim_J = m_JME-m_JMS;
        std::int32_t total_size = (dim_J*m_ndvel*dim_I*m_kdvel);
        return (static_cast<std::size_t>(total_size));
    }
#endif 
};

struct alignas(64) mynn_bl_driver_wrapper_t
{
    std::int32_t m_initflag;
    std::int32_t m_grav_settling;
		std::int32_t m_bl_mynn_tkebudget;
		std::int32_t m_bl_mynn_cloudpdf;
		std::int32_t m_bl_mixlength;
		std::int32_t m_bl_mynn_edmf;
		std::int32_t m_bl_mynn_tkeadvect;
		std::int32_t m_bl_mynn_edmf_mom;
		std::int32_t m_bl_mynn_edmf_tke;
		std::int32_t m_bl_mynn_edmf_part;
		std::int32_t m_bl_mynn_cloudmix;
		std::int32_t m_bl_mynn_mixqt;
		std::int32_t m_icloud_bl;
		std::int32_t m_FLAG_QI;
		std::int32_t m_FLAG_QNI;
		std::int32_t m_FLAG_QC;
		std::int32_t m_FLAG_QNC;
		std::int32_t m_IDS;
		std::int32_t m_IDE;
		std::int32_t m_JDS;
		std::int32_t m_JDE;
		std::int32_t m_KDS;
		std::int32_t m_KDE;
		std::int32_t m_IMS;
		std::int32_t m_IME;
		std::int32_t m_JMS;
		std::int32_t m_JME;
		std::int32_t m_KMS;
		std::int32_t m_KME;
		std::int32_t m_ITS;
		std::int32_t m_ITE;
		std::int32_t m_JTS;
		std::int32_t m_JTE;
		std::int32_t m_KTS;
		std::int32_t m_KTE;
#if WRF_CHEM == 1
		std::int32_t m_nchem;
		std::int32_t m_kdvel;
		std::int32_t m_ndvel;
		std::int32_t m_num_vert_mix;
#endif
		// Scalar real variables.
		float m_delt;
		float m_dx;
		// Array real variables.
		// Arrays 3D - Input.
    darray_r4_t  m_dz;
    darray_r4_t  m_u;
    darray_r4_t  m_v;
    darray_r4_t  m_w;
    darray_r4_t  m_th;
    darray_r4_t  m_qv;
    darray_r4_t  m_p;
    darray_r4_t  m_exner;
    darray_r4_t  m_rho;
    darray_r4_t  m_T3D;
    darray_r4_t  m_qc;
    darray_r4_t  m_qi;
    darray_r4_t  m_qni;
    darray_r4_t  m_qnc;
    // Arrays 2D input
    darray_r4_t  m_xland;
    darray_r4_t  m_ust;
    darray_r4_t  m_ch;
    darray_r4_t  m_rmol;
    darray_r4_t  m_ts;
    darray_r4_t  m_qsfc;
    darray_r4_t  m_qcg;
    darray_r4_t  m_ps;
    darray_r4_t  m_hfx;
    darray_r4_t  m_qfx;
    darray_r4_t  m_wspd;
    darray_r4_t  m_uoce;
    darray_r4_t  m_voce;
    darray_r4_t  m_vdfg;
    darray_r4_t  m_znt;
// Arrays 3D - Input/Output.
    darray_r4_t  m_QKE;
    darray_r4_t  m_Tsq;
    darray_r4_t  m_Qsq;
    darray_r4_t  m_Cov;
    darray_r4_t  m_tke_pbl;
    darray_r4_t  m_qke_adv;
    darray_r4_t  m_RUBLTEN;
    darray_r4_t  m_RVBLTEN;
    darray_r4_t  m_RTHBLTEN;
    darray_r4_t  m_RQVBLTEN;
    darray_r4_t  m_RQCBLTEN;
    darray_r4_t  m_RQIBLTEN;
    darray_r4_t  m_RQNBLTEN;
    // Arrays 3D - Output.
    darray_r4_t m_exch_h;
    darray_r4_t m_exch_m;
   // Arrays 3D - Input/Output , [optional]
    std::optional<darray_r4_t> m_edmf_a;
    std::optional<darray_r4_t> m_edmf_w;
    std::optional<darray_r4_t> m_edmf_qt;
    std::optional<darray_r4_t> m_edmf_thl;
    std::optional<darray_r4_t> m_edmf_ent;
    std::optional<darray_r4_t> m_edmf_qc;
		// Arrays 2D - Input/Output.
    darray_r4_t m_Pblh;
    darray_r4_t m_wstar;
    darray_r4_t m_delta;
    darray_r4_t m_KPBL;
  // Arrays 3D - Input/Output and Output.
    darray_r4_t m_el_pbl;
    darray_r4_t m_qWT;
    darray_r4_t m_qSHEAR;
    darray_r4_t m_qBUOY;
    darray_r4_t m_qDISS;
    darray_r4_t m_dqke;
	// Arrays  3D.
    darray_r4_t m_Sh3D;
    darray_r4_t m_qc_bl;
    darray_r4_t m_cldfra_bl;
	  // defined if WRF_CHEM == 1 (F90 source)
#if WRF_CHEM == 1
    darray_r4_t m_chem3D;
    darray_r4_t m_vd3d;
#endif
		static constexpr const std::int32_t m_total_arrays   = 63;
		static constexpr const std::int32_t m_total_arrays2D = 19;
		static constexpr const std::int32_t m_total_arrays3D = 44;

    mynn_bl_driver_wrapper_t() noexcept(true)
        : 
      m_initflag{},
      m_grav_settling{},
		  m_bl_mynn_tkebudget{},
		  m_bl_mynn_cloudpdf{},
		  m_bl_mixlength{},
		  m_bl_mynn_edmf{},
		  m_bl_mynn_tkeadvect{},
		  m_bl_mynn_edmf_mom{},
		  m_bl_mynn_edmf_tke{},
		  m_bl_mynn_edmf_part{},
		  m_bl_mynn_cloudmix{},
		  m_bl_mynn_mixqt{},
		  m_icloud_bl{},
		  m_FLAG_QI{},
		  m_FLAG_QNI{},
		  m_FLAG_QC{},
		  m_FLAG_QNC{},
		  m_IDS{},
		  m_IDE{},
		  m_JDS{},
		  m_JDE{},
		  m_KDS{},
		  m_KDE{},
		  m_IMS{},
		  m_IME{},
		  m_JMS{},
		  m_JME{},
		  m_KMS{},
		  m_KME{},
		  m_ITS{},
		  m_ITE{},
		  m_JTS{},
		  m_JTE{},
		  m_KTS{},
		  m_KTE{},
#if WRF_CHEM == 1
		  m_nchem{},
		  m_kdvel{},
		  m_ndvel{},
		  m_num_vert_mix{},
#endif
		// Scalar real variables.
		  m_delt{},
		  m_dx{},
      m_dz(),
      m_u(),
      m_v(),
      m_w(),
      m_th(),
      m_qv(),
      m_p(),
      m_exner(),
      m_rho(),
      m_T3D(),
      m_qc(),
      m_qi(),
      m_qni(),
      m_qnc(),
      m_xland(),
      m_ust(),
      m_ch(),
      m_rmol(),
      m_ts(),
      m_qsfc(),
      m_qcg(),
      m_ps(),
      m_hfx(),
      m_qfx(),
      m_wspd(),
      m_uoce(),
      m_voce(),
      m_vdfg(),
      m_znt(),
      m_QKE(),
      m_Tsq(),
      m_Qsq(),
      m_Cov(),
      m_tke_pbl(),
      m_qke_adv(),
      m_RUBLTEN(),
      m_RVBLTEN(),
      m_RTHBLTEN(),
      m_RQVBLTEN(),
      m_RQCBLTEN(),
      m_RQIBLTEN(),
      m_RQNBLTEN(),
      m_exch_h(),
      m_exch_m(),
      m_edmf_a(),
      m_edmf_w(),
      m_edmf_qt(),
      m_edmf_thl(),
      m_edmf_ent(),
      m_edmf_qc(),
      m_Pblh(),
      m_wstar(),
      m_delta(),
      m_KPBL(),
      m_el_pbl(),
      m_qWT(),
      m_qSHEAR(),
      m_qBUOY(),
      m_qDISS(),
      m_dqke(),
      m_Sh3D(),
      m_qc_bl(),
      m_cldfra_bl(),
#if WRF_CHEM == 1
      m_chem3D(),
      m_vd3d()
#endif
      {

      }

      mynn_bl_driver_wrapper_t(const mynn_bl_driver_args_holder_t &holder) noexcept(false) 
         : 
      m_initflag(holder.m_initflag),
      m_grav_settling(holder.m_grav_settling),
		  m_bl_mynn_tkebudget(holder.m_bl_mynn_tkebudget),
		  m_bl_mynn_cloudpdf(holder.m_bl_mynn_cloudpdf),
		  m_bl_mixlength(holder.m_bl_mixlength),
		  m_bl_mynn_edmf(holder.m_bl_mynn_edmf),
		  m_bl_mynn_tkeadvect(holder.m_bl_mynn_tkeadvect),
		  m_bl_mynn_edmf_mom(holder.m_bl_mynn_edmf_mom),
		  m_bl_mynn_edmf_tke(holder.m_bl_mynn_edmf_tke),
		  m_bl_mynn_edmf_part(holder.m_bl_mynn_edmf_part),
		  m_bl_mynn_cloudmix(holder.m_bl_mynn_cloudmix),
		  m_bl_mynn_mixqt(holder.m_bl_mynn_mixqt),
		  m_icloud_bl(holder.m_icloud_bl),
		  m_FLAG_QI(holder.m_FLAG_QI),
		  m_FLAG_QNI(holder.m_FLAG_QNI),
		  m_FLAG_QC(holder.m_FLAG_QC),
		  m_FLAG_QNC(holder.m_FLAG_QNC),
		  m_IDS(holder.m_IDS),
		  m_IDE(holder.m_IDE),
		  m_JDS(holder.m_JDS),
		  m_JDE(holder.m_JDE),
		  m_KDS(holder.m_KDS),
		  m_KDE(holder.m_KDE),
		  m_IMS(holder.m_IMS),
		  m_IME(holder.m_IME),
		  m_JMS(holder.m_JMS),
		  m_JME(holder.m_JME),
		  m_KMS(holder.m_KMS),
		  m_KME(holder.m_KME),
		  m_ITS(holder.m_ITS),
		  m_ITE(holder.m_ITE),
		  m_JTS(holder.m_JTS),
		  m_JTE(holder.m_JTE),
		  m_KTS(holder.m_KTS),
		  m_KTE(holder.m_KTE),
#if WRF_CHEM == 1
		  m_nchem(holder.m_nchem),
		  m_kdvel(holder.m_ldvel),
		  m_ndvel(holder.m_ndvel),
		  m_num_vert_mix(holder.m_num_vert_mix),
#endif
      m_delt(holder.m_delt),
		  m_dx(holder.m_dx),
      m_dz(holder.size3d_row_major()),
      m_u(holder.size3d_row_major()),
      m_v(holder.size3d_row_major()),
      m_w(holder.size3d_row_major()),
      m_th(holder.size3d_row_major()),
      m_qv(holder.size3d_row_major()),
      m_p(holder.size3d_row_major()),
      m_exner(holder.size3d_row_major()),
      m_rho(holder.size3d_row_major()),
      m_T3D(holder.size3d_row_major()),
      m_qc(holder.size3d_row_major()),
      m_qi(holder.size3d_row_major()),
      m_qni(holder.size3d_row_major()),
      m_qnc(holder.size3d_row_major()),
      m_xland(holder.size2d_row_major()),
      m_ust(holder.size2d_row_major()),
      m_ch(holder.size2d_row_major()),
      m_rmol(holder.size2d_row_major()),
      m_ts(holder.size2d_row_major()),
      m_qsfc(holder.size2d_row_major()),
      m_qcg(holder.size2d_row_major()),
      m_ps(holder.size2d_row_major()),
      m_hfx(holder.size2d_row_major()),
      m_qfx(holder.size2d_row_major()),
      m_wspd(holder.size2d_row_major()),
      m_uoce(holder.size2d_row_major()),
      m_voce(holder.size2d_row_major()),
      m_vdfg(holder.size2d_row_major()),
      m_znt(holder.size2d_row_major()),
      m_QKE(holder.size3d_row_major()),
      m_Tsq(holder.size3d_row_major()),
      m_Qsq(holder.size3d_row_major()),
      m_Cov(holder.size3d_row_major()),
      m_tke_pbl(holder.size3d_row_major()),
      m_qke_adv(holder.size3d_row_major()),
      m_RUBLTEN(holder.size3d_row_major()),
      m_RVBLTEN(holder.size3d_row_major()),
      m_RTHBLTEN(holder.size3d_row_major()),
      m_RQVBLTEN(holder.size3d_row_major()),
      m_RQCBLTEN(holder.size3d_row_major()),
      m_RQIBLTEN(holder.size3d_row_major()),
      m_RQNBLTEN(holder.size3d_row_major()),
      m_exch_h(holder.size3d_row_major()),
      m_exch_m(holder.size3d_row_major()),
      m_edmf_a(std::in_place,holder.size3d_row_major()),
      m_edmf_w(std::in_place,holder.size3d_row_major()),
      m_edmf_qt(std::in_place,holder.size3d_row_major()),
      m_edmf_thl(std::in_place,holder.size3d_row_major()),
      m_edmf_ent(std::in_place,holder.size3d_row_major()),
      m_edmf_qc(std::in_place,holder.size3d_row_major()),
      m_Pblh(holder.size2d_row_major()),
      m_wstar(holder.size2d_row_major()),
      m_delta(holder.size2d_row_major()),
      m_KPBL(holder.size2d_row_major()),
      m_el_pbl(holder.size3d_row_major()),
      m_qWT(holder.size3d_row_major()),
      m_qSHEAR(holder.size3d_row_major()),
      m_qBUOY(holder.size3d_row_major()),
      m_qDISS(holder.size3d_row_major()),
      m_dqke(holder.size3d_row_major()),
      m_Sh3D(holder.size3d_row_major()),
      m_qc_bl(holder.size3d_row_major()),
      m_cldfra_bl(holder.size3d_row_major()),
#if WRF_CHEM == 1
      m_chem3D(holder.size4D_row_major_type1()),
      m_vd3d(holder.size4D_row_major_type2())
#endif
      {

      }

      ~mynn_bl_driver_wrapper_t() 
      {

      }

      mynn_bl_driver_wrapper_t(const mynn_bl_driver_wrapper_t &) = delete;

      mynn_bl_driver_wrapper_t(mynn_bl_driver_wrapper_t &&) = delete;

      mynn_bl_driver_wrapper_t &
                    operator=(const mynn_bl_driver_wrapper_t &) = delete;

      mynn_bl_driver_wrapper_t &
                    operator=(mynn_bl_driver_wrapper_t &&) = delete;


      std::int32_t call_mynn_bl_driver() 
      {
           std::int32_t ret_status = -9999;
           const bool is_present_edmf_a   = this->m_edmf_a.has_value();
           const bool is_present_edmf_w   = this->m_edmf_w.has_value();
           const bool is_present_edmf_qt  = this->m_edmf_qt.has_value();
           const bool is_present_edmf_thl = this->m_edmf_thl.has_value();
           const bool is_present_edmf_ent = this->m_edmf_ent.has_value();
           const bool is_present_edmf_qc  = this->m_edmf_qc.has_value();
           if(is_present_edmf_a   && 
              is_present_edmf_w   && 
              is_present_edmf_qt  && 
              is_present_edmf_thl && 
              is_present_edmf_qc  && 
              is_present_edmf_ent) 
            {  
 
//#if defined(__INTEL_COMPILER) || defined(__ICC)

#if (MYNN_BL_DRIVER_USE_PMC_INSTRUMENTATION) == 1
           PMC_VARS
           HW_PMC_COLLECTION_PROLOGE_BODY
#endif
           MODULE_BL_MYNN_mp_MYNN_BL_DRIVER(&this->m_bl_mixlength, &this->m_grav_settling,
							 &this->m_delt, &this->m_dz.m_data[0], &this->m_dx, &this->m_znt.m_data[0],
							 &this->m_u.m_data[0], &this->m_v.m_data[0], &this->m_w.m_data[0], &this->m_th.m_data[0], &this->m_qv.m_data[0], &this->m_qc.m_data[0],
							 &this->m_qi.m_data[0], &this->m_qni.m_data[0], &this->m_qnc.m_data[0], &this->m_p.m_data[0], &this->m_exner.m_data[0], &this->m_rho.m_data[0],
							 &this->m_T3D.m_data[0], &this->m_xland.m_data[0], &this->m_ts.m_data[0], &this->m_qsfc.m_data[0], &this->m_qcg.m_data[0], &this->m_ps.m_data[0],
							 &this->m_ust.m_data[0], &this->m_ch.m_data[0], &this->m_hfx.m_data[0], &this->m_qfx.m_data[0], &this->m_rmol.m_data[0], &this->m_wspd.m_data[0],
							 &this->m_uoce.m_data[0], &this->m_voce.m_data[0], &this->m_vdfg.m_data[0], &this->m_Qke.m_data[0], &this->m_tke_pbl.m_data[0],
							 &this->m_qke_adv.m_data[0], &this->m_bl_mynn_tkeadvect,
#if WRF_CHEM == 1
							 &this->m_chem3D.m_data[0], &this->m_vd3d.m_data[0], &this->m_nchem,
							 &this->m_kdvel, &this->m_ndvel, &this->m_num_vert_mix,
#endif							
							 &this->m_Tsq.m_data[0], &this->m_Qsq.m_data[0], &this->m_Cov.m_data[0],
							 &this->m_RUBLTEN.m_data[0], &this->m_RVBLTEN.m_data[0], &this->m_RTHBLTEN.m_data[0],
							 &this->m_RQVBLTEN.m_data[0], &this->m_RQCBLTEN.m_data[0], &this->m_RQIBLTEN.m_data[0],
							 &this->m_RQNBLTEN.m_data[0], &this->m_exch_h.m_data[0], &this->m_exch_m.m_data[0],
							 &this->m_Pblh.m_data[0], &this->m_KPBL.m_data[0],
							 &this->m_el_pbl.m_data[0], 
							 &this->m_dqke.m_data[0], &this->m_qWT.m_data[0], &this->m_qSHEAR.m_data[0], &this->m_qBUOY.m_data[0], &this->m_qDISS.m_data[0],
							 &this->m_wstar.m_data[0], &this->m_delta.m_data[0],
							 &this->m_bl_mynn_tkebudget, &this->m_bl_mynn_cloudpdf, &this->m_Sh3D.m_data[0],
							 &this->m_bl_mixlength, &this->m_icloud_bl, &this->m_qc_bl.m_data[0], &this->m_cldfra_bl.m_data[0],
							 &this->m_bl_mynn_edmf, &this->m_bl_mynn_edmf_mom, &this->m_bl_mynn_edmf_tke,
							 &this->m_bl_mynn_edmf_part, &this->m_bl_mynn_cloudmix, &this->m_bl_mynn_mixqt,
							 &this->m_edmf_a.operator->()->m_data[0], &this->m_edmf_w.operator->()->m_data[0], &this->m_edmf_qt.operator->()->m_data[0],
							 &this->m_edmf_thl.operator->()->m_data[0], &this->m_edmf_ent.operator->()->m_data[0], &this->m_edmf_qc.operator->()->m_data[0],
							 &this->m_FLAG_QI, &this->m_FLAG_QNI, &this->m_FLAG_QC, &this->m_FLAG_QNC,
							 &this->m_IDS, &this->m_IDE, &this->m_JDS, &this->m_JDE, &this->m_KDS, &this->m_KDE,
							 &this->m_IMS, &this->m_IME, &this->m_JMS, &this->m_JME, &this->m_KMS, &this->m_KME,
							 &this->m_ITS, &this->m_ITE, &this->m_JTS, &this->m_JTE, &this->m_KTS, &this->m_KTE);
#if (MYNN_BL_DRIVER_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif 
               ret_status = 0;
               return (ret_status);
          }
          else
          {
              ret_status = -1;
              return (ret_status);
          }
      }
};


}

}


#endif /*__GMS_MYNN_BL_DRIVER_H__*/