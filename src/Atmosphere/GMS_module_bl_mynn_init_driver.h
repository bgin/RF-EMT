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

#ifndef __GMS_MODULE_BL_MYNN_INIT_DRIVER_H__
#define __GMS_MODULE_BL_MYNN_INIT_DRIVER_H__ 

#include <cstdint>
#include "GMS_config.h"
#include "GMS_module_bl_mynn_c_iface.h"
#include "GMS_dyn_array.h"

namespace file_info 
{

     const unsigned int GMS_MODULE_BL_MYNN_INIT_DRIVER_MAJOR = 1;
     const unsigned int GMS_MODULE_BL_MYNN_INIT_DRIVER_MINOR = 0;
     const unsigned int GMS_MODULE_BL_MYNN_INIT_DRIVER_MICRO = 0;
     const unsigned int GMS_MODULE_BL_MYNN_INIT_DRIVER_FULLVER =
       1000U*GMS_MODULE_BL_MYNN_INIT_DRIVER_MAJOR+100U*GMS_MODULE_BL_MYNN_INIT_DRIVER_MINOR+
       10U*GMS_MODULE_BL_MYNN_INIT_DRIVER_MICRO;
     const char * const GMS_MODULE_BL_MYNN_INIT_DRIVER_CREATION_DATE = "05-11-2026 10:03 AM +00200 (THR 11 NOV 2026 GMT+2)";
     const char * const GMS_MODULE_BL_MYNN_INIT_DRIVER_BUILD_DATE    = __DATE__ " " __TIME__;
     static const char  GMS_MODULE_BL_MYNN_INIT_DRIVER_AUTHOR[]      = "Programmer: Bernard Gingold, beniekg@gmail.com"; 
     const char * const GMS_MODULE_BL_MYNN_INIT_DRIVER_SYNOPSIS      = "Fortran 90 Module module_bl_mynn driver, C++ initialization wrapper.";

}


// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (MODULE_BL_MYNN_INIT_DRIVER_USE_NT_STORES)
#define MODULE_BL_MYNN_INIT_DRIVER_USE_NT_STORES 1
#endif

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(MODULE_BL_MYNN_INIT_DRIVER_INIT_STORAGE)
#define MODULE_BL_MYNN_INIT_DRIVER_INIT_STORAGE 1
#endif 
#endif

#if (MODULE_BL_MYNN_INIT_DRIVER_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 



// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (MODULE_BL_MYNN_INIT_DRIVER_USE_PMC_INSTRUMENTATION)
#define MODULE_BL_MYNN_INIT_DRIVER_USE_PMC_INSTRUMENTATION 1
#endif 

#if (MODULE_BL_MYNN_INIT_DRIVER_USE_PMC_INSTRUMENTATION) == 1
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

struct alignas(64) module_bl_mynn_init_driver_t 
{
    std::int32_t  m_ids;
    std::int32_t  m_ide;
    std::int32_t  m_jds;
    std::int32_t  m_jde;
    std::int32_t  m_kds;
    std::int32_t  m_kde;
    std::int32_t  m_ims;
    std::int32_t  m_ime;
    std::int32_t  m_jms;
    std::int32_t  m_jme;
    std::int32_t  m_kms;
    std::int32_t  m_kme;
    std::int32_t  m_its;
    std::int32_t  m_ite;
    std::int32_t  m_jts;
    std::int32_t  m_jte;
    std::int32_t  m_kts;
    std::int32_t  m_kte;
    std::int32_t  m_allowed_to_read;
    std::int32_t  m_restart;
    std::int32_t  m_level;
/*
   Size of every array is: m_IME * m_KME * m_JME
*/
    darray_r4_t   m_RUBLTEN;
    darray_r4_t   m_RVBLTEN;
    darray_r4_t   m_RTHBLTEN;
    darray_r4_t   m_RQVBLTEN;
    darray_r4_t   m_RQCBLTEN;
    darray_r4_t   m_RQIBLTEN;
    darray_r4_t   m_QKE;
    darray_r4_t   m_TKE_PBL;
    darray_r4_t   m_EXCH_H;

    module_bl_mynn_init_driver_t() noexcept(true) 
      : 
    m_ids{},
    m_ide{},
    m_jds{},
    m_jde{},
    m_kds{},
    m_kde{},
    m_ims{},
    m_ime{},
    m_jms{},
    m_jme{},
    m_kms{},
    m_kme{},
    m_its{},
    m_ite{},
    m_jts{},
    m_jte{},
    m_kts{},
    m_kte{},
    m_allowed_to_read{},
    m_restart{},
    m_level{},
    m_RUBLTEN(),
    m_RVBLTEN(),
    m_RTHBLTEN(),
    m_RQVBLTEN(),
    m_RQCBLTEN(),
    m_RQIBLTEN(),
    m_QKE(),
    m_TKE_PBL(),
    m_EXCH_H()
    {
      
    }
};

}

}

#endif /*__GMS_MODULE_BL_MYNN_INIT_DRIVER_H__*/