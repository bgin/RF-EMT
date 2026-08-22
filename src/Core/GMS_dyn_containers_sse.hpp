
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

#ifndef __GMS_DYN_CONTAINERS_SSE_HPP__
#define __GMS_DYN_CONTAINERS_SSE_HPP__ 091220230819


namespace file_info {

     const unsigned int GMS_DYN_CONTAINERS_SSE_MAJOR = 1;
     const unsigned int GMS_DYN_CONTAINERS_SSE_MINOR = 1;
     const unsigned int GMS_DYN_CONTAINERS_SSE_MICRO = 0;
     const unsigned int GMS_DYN_CONTAINERS_SSE_FULLVER =
       1000U*GMS_DYN_CONTAINERS_SSE_MAJOR+100U*GMS_DYN_CONTAINERS_SSE_MINOR+
       10U*GMS_DYN_CONTAINERS_SSE_MICRO;
     const char * const GMS_DYN_CONTAINERS_SSE_CREATION_DATE = "09-12-2023 08:19 AM +00200 (SAT 09 DEC 2023 GMT+2)";
     const char * const GMS_DYN_CONTAINERS_SSE_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_DYN_CONTAINERS_SSE_SYNOPSIS      = "Dynamically allocated, 64-byte aligned dense SSE containers."

}





#include <cstdint>
#include "GMS_config.h"
#include "GMS_complex_xmm2r8.hpp"
#include "GMS_complex_xmm4r4.hpp"
#include "GMS_malloc.h"
#include "GMS_simd_memops.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_GMS_DYN_CONTAINERS_SSE_NT_STORES)
#define USE_GMS_DYN_CONTAINERS_SSE_NT_STORES 0
#endif

 
namespace gms {


           struct __ATTR_ALIGN__(64) DC3D_xmm2c8_t {
                      // Complmx electric field.
                      xmm2c8_t * __restrict mx;
                      xmm2c8_t * __restrict my;
                      xmm2c8_t * __restrict mz;
                      std::size_t                      mnx;
                      std::size_t                      mny;
                      std::size_t                      mnz;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,15)
#endif
                     inline DC3D_xmm2c8_t()  {
                          
                          this->mnx  = 0ULL;
                          this->mny  = 0ULL;
                          this->mnz  = 0ULL;
                          this->mx  = NULL;
                          this->my  = NULL;
                          this->mz  = NULL;
                      } 
                          
                     inline DC3D_xmm2c8_t(const std::size_t _mnx,
                                   const std::size_t _mny,
                                   const std::size_t _mnz) {
                          using namespace gms::common;
                          this->mnx = _mnx;
                          this->mny = _mny;
                          this->mnz = _mnz;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline DC3D_xmm2c8_t(const std::size_t _mnx,
                                   const std::size_t _mny,
                                   const std::size_t _mnz,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->mnx = _mnx;
                             this->mny = _mny;
                             this->mnz = _mnz;
                             switch (fsize) {
                                 case 0:
                                      this->mx = (xmm2c8_t*)
                                                 gms_mmap_4KiB<xmm2c8_t>(this->mnx,prot,flags,fd,offset);
                                      this->my = (xmm2c8_t*)
                                                 gms_mmap_4KiB<xmm2c8_t>(this->mny,prot,flags,fd,offset);
                                      this->mz = (xmm2c8_t*)
                                                 gms_mmap_4KiB<xmm2c8_t>(this->mnz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 
                                 break;
                                 case 1:
                                      this->mx = (xmm2c8_t*)
                                                 gms_mmap_2MiB<xmm2c8_t>(this->mnx,prot,flags,fd,offset);
                                      this->my = (xmm2c8_t*)
                                                 gms_mmap_2MiB<xmm2c8_t>(this->mny,prot,flags,fd,offset);
                                      this->mz = (xmm2c8_t*)
                                                 gms_mmap_2MiB<xmm2c8_t>(this->mnz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 2:
                                      this->mx = (xmm2c8_t*)
                                                 gms_mmap_1GiB<xmm2c8_t>(this->mnx,prot,flags,fd,offset);
                                      this->my = (xmm2c8_t*)
                                                 gms_mmap_1GiB<xmm2c8_t>(this->mny,prot,flags,fd,offset);
                                      this->mz = (xmm2c8_t*)
                                                 gms_mmap_1GiB<xmm2c8_t>(this->mnz,prot,flags,fd,offset); 
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                                                 
                            
                      
                     inline DC3D_xmm2c8_t(const std::size_t _mnx,
                                          const std::size_t _mny,
                                          const std::size_t _mnz,
                                          const  xmm2c8_t * __restrict m_x,
                                          const  xmm2c8_t * __restrict m_y,
                                          const  xmm2c8_t * __restrict m_z) {
                          using namespace gms::common;
                          this->mnx = _mnx;
                          this->mny = _mny;
                          this->mnz = _mnz;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_DYN_CONTAINERS_SSE_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_uncached_memmove(&this->my[0],&m_y[0],this->mny);
	                  avx512_uncached_memmove(&this->mz[0],&m_z[0],this->mnz);
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_cached_memmove(&this->my[0],&m_y[0],this->mny);
	                  avx512_cached_memmove(&this->mz[0],&m_z[0],this->mnz);
#endif
                   }  
                   
                    inline  DC3D_xmm2c8_t(DC3D_xmm2c8_t && rhs) {
                          
                          this->mnx = rhs.mnx;
                          this->mny = rhs.mny;
                          this->mnz = rhs.mnz;
                          this->mx   = &rhs.mx[0];
                          this->my   = &rhs.my[0];
                          this->mz   = &rhs.mz[0];
                          rhs.mnx    = 0ULL;
                          rhs.mny    = 0ULL;
                          rhs.mnz    = 0ULL;
                          rhs.mx     = NULL;
                          rhs.my     = NULL;
                          rhs.mz     = NULL;
                      }
                                 
                     DC3D_xmm2c8_t(const DC3D_xmm2c8_t &)             = delete;
                      
                     inline ~DC3D_xmm2c8_t() {
                      
                          using namespace gms::common;
                          if(this->issmap) {
                             gms_unmap<xmm2c8_t>(this->mx,this->mnx);
                             gms_unmap<xmm2c8_t>(this->my,this->mny);
                             gms_unmap<xmm2c8_t>(this->mz,this->mnz);
                          }
                          else {
                              gms_mm_free(this->mx);
                              gms_mm_free(this->my);
                              gms_mm_free(this->mz);
                          }
                      }
                      
                     DC3D_xmm2c8_t & operator=(const DC3D_xmm2c8_t &) = delete;
                      
                     inline DC3D_xmm2c8_t & operator=(DC3D_xmm2c8_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           if(this->ismmap) {
                             gms_unmap<xmm2c8_t>(this->mx,this->mnx);
                             gms_unmap<xmm2c8_t>(this->my,this->mny);
                             gms_unmap<xmm2c8_t>(this->mz,this->mnz);
                           }
                           else {
                             gms_mm_free(this->mx);
                             gms_mm_free(this->my);
                             gms_mm_free(this->mz);
                           }
                           this->mnx   = rhs.mnx;
                           this->mny   = rhs.mny;
                           this->mnz   = rhs.mnz;
                           this->mx    = &rhs.mx[0];
                           this->my    = &rhs.my[0];
                           this->mz    = &rhs.mz[0];
                           rhs.mnx     = 0ULL;
                           rhs.mny     = 0ULL;
                           rhs.mnz     = 0ULL;
                           rhs.mx      = NULL;
                           rhs.my      = NULL;
                           rhs.mz      = NULL;
                           return (*this);
                      }
                      
                      inline void allocate() {
                          using namespace gms::common;
                          this->mx  = (xmm2c8_t*)
                                         gms_mm_malloc( sizeof(xmm2c8_t)*this->mnx,64ULL);
                          this->my  = (xmm2c8_t*)
                                         gms_mm_malloc( sizeof(xmm2c8_t)*this->mny,64ULL);
                          this->mz  = (xmm2c8_t*)
                                         gms_mm_malloc( sizeof(xmm2c8_t)*this->mny,64ULL);
                      }
                      
                     
                    
              };
              
              
              struct __ATTR_ALIGN__(64) DC2D_xmm2c8_t {
                     
                      xmm2c8_t * __restrict mx
                      xmm2c8_t * __restrict my;
                      std::size_t                      mnx;
                      std::size_t                      mny;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,23)
#endif
                     inline DC2D_xmm2c8_t()  {
                          
                          this->mnx  = 0ULL;
                          this->mny  = 0ULL;
                          this->mx  = NULL;
                          this->my  = NULL;
                    } 
                          
                     inline DC2D_xmm2c8_t(const std::size_t _mnx,
                                   const std::size_t _mny) {
                                   
                          using namespace gms::common;
                          this->mnx = _mnx;
                          this->mny = _mny;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline DC2D_xmm2c8_t(const std::size_t _mnx,
                                   const std::size_t _mny,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->mnx = _mnx;
                             this->mny = _mny;
                              switch (fsize) {
                                 case:0
                                      this->mx = (xmm2c8_t*)
                                                 gms_mmap_4KiB<xmm2c8_t>(this->mnx,prot,flags,fd,offset);
                                      this->my = (xmm2c8_t*)
                                                 gms_mmap_4KiB<xmm2c8_t>(this->mny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mx = (xmm2c8_t*)
                                                 gms_mmap_2MiB<xmm2c8_t>(this->mnx,prot,flags,fd,offset);
                                      this->my = (xmm2c8_t*)
                                                 gms_mmap_2MiB<xmm2c8_t>(this->mny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mx = (xmm2c8_t*)
                                                 gms_mmap_1GiB<xmm2c8_t>(this->mnx,prot,flags,fd,offset);
                                      this->my = (xmm2c8_t*)
                                                 gms_mmap_1GiB<xmm2c8_t>(this->mny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                                                 
                            
                      
                     inline DC2D_xmm2c8_t(const std::size_t _mnx,
                                          const std::size_t _mny,
                                          const  xmm2c8_t * __restrict m_x,
                                          const  xmm2c8_t * __restrict m_y) {
                                          
                          using namespace gms::common;
                          this->mnx = _mnx;
                          this->mny = _mny;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_DYN_CONTAINERS_SSE_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_uncached_memmove(&this->my[0],&m_y[0],this->mny);
	                  
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_cached_memmove(&this->my[0],&m_y[0],this->mny);
	               
#endif
                   }  
                   
                    inline  DC2D_xmm2c8_t(DC2D_xmm2c8_t && rhs) {
                          
                          this->mnx = rhs.mnx;
                          this->mny = rhs.mny;
                          this->mx   = &rhs.mx[0];
                          this->my   = &rhs.my[0];
                          rhs.mnx    = 0ULL;
                          rhs.mny    = 0ULL;
                          rhs.mx     = NULL;
                          rhs.my     = NULL;
                        
                      }
                                 
                     DC2D_xmm2c8_t(const DC2D_xmm2c8_t &)             = delete;
                      
                     inline ~DC2D_xmm2c8_t() {
                      
                          using namespace gms::common;
                          if(this->issmap) {
                             gms_unmap<xmm2c8_t>(this->mx,this->mnx);
                             gms_unmap<xmm2c8_t>(this->my,this->mny);
                             
                          }
                          else {
                              gms_mm_free(this->mx);
                              gms_mm_free(this->my);
                              
                          }
                      }
                      
                     DC2D_xmm2c8_t & operator=(const DC2D_xmm2c8_t &) = delete;
                      
                     inline DC2D_xmm2c8_t & operator=(DC2D_xmm2c8_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           if(this->ismmap) {
                             gms_unmap<xmm2c8_t>(this->mx,this->mnx);
                             gms_unmap<xmm2c8_t>(this->my,this->mny);
                            
                           }
                           else {
                             gms_mm_free(this->mx);
                             gms_mm_free(this->my);
                           
                           }
                           this->mnx   = rhs.mnx;
                           this->mny   = rhs.mny;
                           this->mx    = &rhs.mx[0];
                           this->my    = &rhs.my[0];
                           rhs.mnx     = 0ULL;
                           rhs.mny     = 0ULL;
                           rhs.mx      = NULL;
                           rhs.my      = NULL;
                           return (*this);
                      }
                      
                      inline void allocate() {
                          using namespace gms::common;
                          this->mx  = (xmm2c8_t*)
                                         gms_mm_malloc( sizeof(xmm2c8_t)*this->mnx,64ULL);
                          this->my  = (xmm2c8_t*)
                                         gms_mm_malloc( sizeof(xmm2c8_t)*this->mny,64ULL);
                    }
                      
                     
                    
              };
              
              
               struct __ATTR_ALIGN__(64) DC1D_xmm2c8_t {
                     
                      xmm2c8_t * __restrict mx
                      std::size_t                      mnx;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,31)
#endif
                     inline DC1D_xmm2c8_t()  {
                          
                          this->mnx  = 0ULL;
                          this->mx  = NULL;
                         
                    } 
                          
                     inline explicit DC1D_xmm2c8_t(const std::size_t _mnx) {
                                 
                          using namespace gms::common;
                          this->mnx = _mnx;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline DC1D_xmm2c8_t(const std::size_t _mnx,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->mnx = _mnx;
                             switch (fsize) {
                                 case:0
                                      this->mx = (xmm2c8_t*)
                                                 gms_mmap_4KiB<xmm2c8_t>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mx = (xmm2c8_t*)
                                                 gms_mmap_2MiB<xmm2c8_t>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mx = (xmm2c8_t*)
                                                 gms_mmap_1GiB<xmm2c8_t>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                                                 
                            
                      
                     inline DC1D_xmm2c8_t(const std::size_t _mnx,
                                          const  xmm2c8_t * __restrict m_x) {
                                     
                          using namespace gms::common;
                          this->mnx = _mnx;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_DYN_CONTAINERS_SSE_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                 
	                  
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                 
	               
#endif
                   }  
                   
                    inline  DC1D_xmm2c8_t(DC1D_xmm2c8_t && rhs) {
                          
                          this->mnx = rhs.mnx;
                          this->mx   = &rhs.mx[0];
                          rhs.mnx    = 0ULL;
                          rhs.mx     = NULL;
                                                 
                      }
                                 
                     DC1D_xmm2c8_t(const DC1D_xmm2c8_t &)             = delete;
                      
                     inline ~DC1D_xmm2c8_t() {
                      
                          using namespace gms::common;
                          if(this->issmap) 
                             gms_unmap<xmm2c8_t>(this->mx,this->mnx);
                          else 
                              gms_mm_free(this->mx);   
                     }        
                          
                                                 
                                  
                     DC1D_xmm2c8_t & operator=(const DC1D_xmm2c8_t &) = delete;
                      
                     inline DC1D_xmm2c8_t & operator=(DC1D_xmm2c8_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           if(this->ismmap) 
                              gms_unmap<xmm2c8_t>(this->mx,this->mnx);
                           else 
                             gms_mm_free(this->mx); 
                                                                                                                         
                           this->mnx   = rhs.mnx;
                           this->mx    = &rhs.mx[0];
                           rhs.mnx     = 0ULL;
                           rhs.mx      = NULL;
                           return (*this);
                      }
                      
                      inline void allocate() {
                          using namespace gms::common;
                          this->mx  = (xmm2c8_t*)
                                         gms_mm_malloc( sizeof(xmm2c8_t)*this->mnx,64ULL);
                         
                    }
                      
                     
                    
              };
              
              
          //////////////////////////////////////////////////////////////////////////////////////////////
          
           struct __ATTR_ALIGN__(64) DC3D_xmm4c4_t {
                      
                      xmm4c4_t * __restrict mx
                      xmm4c4_t * __restrict my;
                      xmm4c4_t * __restrict mz;
                      std::size_t                      mnx;
                      std::size_t                      mny;
                      std::size_t                      mnz;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,15)
#endif
                     inline DC3D_xmm4c4_t()  {
                          
                          this->mnx  = 0ULL;
                          this->mny  = 0ULL;
                          this->mnz  = 0ULL;
                          this->mx  = NULL;
                          this->my  = NULL;
                          this->mz  = NULL;
                      } 
                          
                     inline DC3D_xmm4c4_t(const std::size_t _mnx,
                                   const std::size_t _mny,
                                   const std::size_t _mnz) {
                          using namespace gms::common;
                          this->mnx = _mnx;
                          this->mny = _mny;
                          this->mnz = _mnz;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline DC3D_xmm4c4_t(const std::size_t _mnx,
                                   const std::size_t _mny,
                                   const std::size_t _mnz,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->mnx = _mnx;
                             this->mny = _mny;
                             this->mnz = _mnz;
                             switch (fsize) {
                                 case:0
                                      this->mx = (xmm4c4_t*)
                                                 gms_mmap_4KiB<xmm4c4_t>(this->mnx,prot,flags,fd,offset);
                                      this->my = (xmm4c4_t*)
                                                 gms_mmap_4KiB<xmm4c4_t>(this->mny,prot,flags,fd,offset);
                                      this->mz = (xmm4c4_t*)
                                                 gms_mmap_4KiB<xmm4c4_t>(this->mnz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mx = (xmm4c4_t*)
                                                 gms_mmap_2MiB<xmm4c4_t>(this->mnx,prot,flags,fd,offset);
                                      this->my = (xmm4c4_t*)
                                                 gms_mmap_2MiB<xmm4c4_t>(this->mny,prot,flags,fd,offset);
                                      this->mz = (xmm4c4_t*)
                                                 gms_mmap_2MiB<xmm4c4_t>(this->mnz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mx = (xmm4c4_t*)
                                                 gms_mmap_1GiB<xmm4c4_t>(this->mnx,prot,flags,fd,offset);
                                      this->my = (xmm4c4_t*)
                                                 gms_mmap_1GiB<xmm4c4_t>(this->mny,prot,flags,fd,offset);
                                      this->mz = (xmm4c4_t*)
                                                 gms_mmap_1GiB<xmm4c4_t>(this->mnz,prot,flags,fd,offset); 
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                                                 
                            
                      
                     inline DC3D_xmm4c4_t(const std::size_t _mnx,
                                          const std::size_t _mny,
                                          const std::size_t _mnz,
                                          const  xmm4c4_t * __restrict m_x,
                                          const  xmm4c4_t * __restrict m_y,
                                          const  xmm4c4_t * __restrict m_z) {
                          using namespace gms::common;
                          this->mnx = _mnx;
                          this->mny = _mny;
                          this->mnz = _mnz;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_DYN_CONTAINERS_SSE_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_uncached_memmove(&this->my[0],&m_y[0],this->mny);
	                  avx512_uncached_memmove(&this->mz[0],&m_z[0],this->mnz);
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_cached_memmove(&this->my[0],&m_y[0],this->mny);
	                  avx512_cached_memmove(&this->mz[0],&m_z[0],this->mnz);
#endif
                   }  
                   
                    inline  DC3D_xmm4c4_t(DC3D_xmm4c4_t && rhs) {
                          
                          this->mnx = rhs.mnx;
                          this->mny = rhs.mny;
                          this->mnz = rhs.mnz;
                          this->mx   = &rhs.mx[0];
                          this->my   = &rhs.my[0];
                          this->mz   = &rhs.mz[0];
                          rhs.mnx    = 0ULL;
                          rhs.mny    = 0ULL;
                          rhs.mnz    = 0ULL;
                          rhs.mx     = NULL;
                          rhs.my     = NULL;
                          rhs.mz     = NULL;
                      }
                                 
                     DC3D_xmm4c4_t(const DC3D_xmm4c4_t &)             = delete;
                      
                     inline ~DC3D_xmm4c4_t() {
                      
                          using namespace gms::common;
                          if(this->issmap) {
                             gms_unmap<xmm4c4_t>(this->mx,this->mnx);
                             gms_unmap<xmm4c4_t>(this->my,this->mny);
                             gms_unmap<xmm4c4_t>(this->mz,this->mnz);
                          }
                          else {
                              gms_mm_free(this->mx);
                              gms_mm_free(this->my);
                              gms_mm_free(this->mz);
                          }
                      }
                      
                     DC3D_xmm4c4_t & operator=(const DC3D_xmm4c4_t &) = delete;
                      
                     inline DC3D_xmm4c4_t & operator=(DC3D_xmm4c4_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           if(this->ismmap) {
                             gms_unmap<xmm4c4_t>(this->mx,this->mnx);
                             gms_unmap<xmm4c4_t>(this->my,this->mny);
                             gms_unmap<xmm4c4_t>(this->mz,this->mnz);
                           }
                           else {
                             gms_mm_free(this->mx);
                             gms_mm_free(this->my);
                             gms_mm_free(this->mz);
                           }
                           this->mnx   = rhs.mnx;
                           this->mny   = rhs.mny;
                           this->mnz   = rhs.mnz;
                           this->mx    = &rhs.mx[0];
                           this->my    = &rhs.my[0];
                           this->mz    = &rhs.mz[0];
                           rhs.mnx     = 0ULL;
                           rhs.mny     = 0ULL;
                           rhs.mnz     = 0ULL;
                           rhs.mx      = NULL;
                           rhs.my      = NULL;
                           rhs.mz      = NULL;
                           return (*this);
                      }
                      
                      inline void allocate() {
                          using namespace gms::common;
                          this->mx  = (xmm4c4_t*)
                                         gms_mm_malloc( sizeof(xmm4c4_t)*this->mnx,64ULL);
                          this->my  = (xmm4c4_t*)
                                         gms_mm_malloc( sizeof(xmm4c4_t)*this->mny,64ULL);
                          this->mz  = (xmm4c4_t*)
                                         gms_mm_malloc( sizeof(xmm4c4_t)*this->mny,64ULL);
                      }
                      
                     
                    
              };
              
              
              
             struct __ATTR_ALIGN__(64) DC2D_xmm4c4_t {
                     
                      xmm4c4_t * __restrict mx
                      xmm4c4_t * __restrict my;
                      std::size_t                      mnx;
                      std::size_t                      mny;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,23)
#endif
                     inline DC2D_xmm4c4_t()  {
                          
                          this->mnx  = 0ULL;
                          this->mny  = 0ULL;
                          this->mx  = NULL;
                          this->my  = NULL;
                    } 
                          
                     inline DC2D_xmm4c4_t(const std::size_t _mnx,
                                          const std::size_t _mny) {
                                   
                          using namespace gms::common;
                          this->mnx = _mnx;
                          this->mny = _mny;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline DC2D_xmm4c4_t(const std::size_t _mnx,
                                   const std::size_t _mny,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->mnx = _mnx;
                             this->mny = _mny;
                              switch (fsize) {
                                 case:0
                                      this->mx = (xmm4c4_t*)
                                                 gms_mmap_4KiB<xmm4c4_t>(this->mnx,prot,flags,fd,offset);
                                      this->my = (xmm4c4_t*)
                                                 gms_mmap_4KiB<xmm4c4_t>(this->mny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mx = (xmm4c4_t*)
                                                 gms_mmap_2MiB<xmm4c4_t>(this->mnx,prot,flags,fd,offset);
                                      this->my = (xmm4c4_t*)
                                                 gms_mmap_2MiB<xmm4c4_t>(this->mny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mx = (xmm4c4_t*)
                                                 gms_mmap_1GiB<xmm4c4_t>(this->mnx,prot,flags,fd,offset);
                                      this->my = (xmm2c8_t*)
                                                 gms_mmap_1GiB<xmm4c4_t>(this->mny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                                                 
                            
                      
                     inline DC2D_xmm4c4_t(const std::size_t _mnx,
                                          const std::size_t _mny,
                                          const  xmm4c4_t * __restrict m_x,
                                          const  xmm4c4_t * __restrict m_y) {
                                          
                          using namespace gms::common;
                          this->mnx = _mnx;
                          this->mny = _mny;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_DYN_CONTAINERS_SSE_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_uncached_memmove(&this->my[0],&m_y[0],this->mny);
	                  
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_cached_memmove(&this->my[0],&m_y[0],this->mny);
	               
#endif
                   }  
                   
                    inline  DC2D_xmm4c4_t(DC2D_xmm4c4_t && rhs) {
                          
                          this->mnx = rhs.mnx;
                          this->mny = rhs.mny;
                          this->mx   = &rhs.mx[0];
                          this->my   = &rhs.my[0];
                          rhs.mnx    = 0ULL;
                          rhs.mny    = 0ULL;
                          rhs.mx     = NULL;
                          rhs.my     = NULL;
                        
                      }
                                 
                     DC2D_xmm4c4_t(const DC2D_xmm4c4_t &)             = delete;
                      
                     inline ~DC2D_xmm4c4_t() {
                      
                          using namespace gms::common;
                          if(this->issmap) {
                             gms_unmap<xmm4c4_t>(this->mx,this->mnx);
                             gms_unmap<xmm4c4_t>(this->my,this->mny);
                             
                          }
                          else {
                              gms_mm_free(this->mx);
                              gms_mm_free(this->my);
                              
                          }
                      }
                      
                     DC2D_xmm4c4_t & operator=(const DC2D_xmm4c4_t &) = delete;
                      
                     inline DC2D_xmm4c4_t & operator=(DC2D_xmm4c4_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           if(this->ismmap) {
                             gms_unmap<xmm4c4_t>(this->mx,this->mnx);
                             gms_unmap<xmm4c4_t>(this->my,this->mny);
                            
                           }
                           else {
                             gms_mm_free(this->mx);
                             gms_mm_free(this->my);
                           
                           }
                           this->mnx   = rhs.mnx;
                           this->mny   = rhs.mny;
                           this->mx    = &rhs.mx[0];
                           this->my    = &rhs.my[0];
                           rhs.mnx     = 0ULL;
                           rhs.mny     = 0ULL;
                           rhs.mx      = NULL;
                           rhs.my      = NULL;
                           return (*this);
                      }
                      
                      inline void allocate() {
                          using namespace gms::common;
                          this->mx  = (xmm4c4_t*)
                                         gms_mm_malloc( sizeof(xmm4c4_t)*this->mnx,64ULL);
                          this->my  = (xmm4c4_t*)
                                         gms_mm_malloc( sizeof(xmm4c4_t)*this->mny,64ULL);
                    }
                      
                     
                    
              }; 
              
              
              struct __ATTR_ALIGN__(64) DC1D_xmm4c4_t {
                     
                      xmm4c4_t * __restrict mx
                      std::size_t                      mnx;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,31)
#endif
                     inline DC1D_xmm4c4_t()  {
                          
                          this->mnx  = 0ULL;
                          this->mx  = NULL;
                         
                    } 
                          
                     inline explicit DC1D_xmm4c4_t(const std::size_t _mnx) {
                                 
                          using namespace gms::common;
                          this->mnx = _mnx;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline DC1D_xmm4c4_t(const std::size_t _mnx,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->mnx = _mnx;
                             switch (fsize) {
                                 case:0
                                      this->mx = (xmm4c4_t*)
                                                 gms_mmap_4KiB<xmm4c4_t>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mx = (xmm4c4_t*)
                                                 gms_mmap_2MiB<xmm4c4_t>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mx = (xmm4c4_t*)
                                                 gms_mmap_1GiB<xmm4c4_t>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                                                 
                            
                      
                     inline DC1D_xmm4c4_t(const std::size_t _mnx,
                                          const  xmm4c4_t * __restrict m_x) {
                                     
                          using namespace gms::common;
                          this->mnx = _mnx;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_DYN_CONTAINERS_SSE_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                 
	                  
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                 
	               
#endif
                   }  
                   
                    inline  DC1D_xmm4c4_t(DC1D_xmm4c4_t && rhs) {
                          
                          this->mnx = rhs.mnx;
                          this->mx   = &rhs.mx[0];
                          rhs.mnx    = 0ULL;
                          rhs.mx     = NULL;
                                                 
                      }
                                 
                     DC1D_xmm4c4_t(const DC1D_xmm4c4_t &)             = delete;
                      
                     inline ~DC1D_xmm4c4_t() {
                      
                          using namespace gms::common;
                          if(this->issmap) 
                             gms_unmap<xmm4c4_t>(this->mx,this->mnx);
                          else 
                              gms_mm_free(this->mx);   
                     }        
                          
                                                 
                                  
                     DC1D_xmm4c4_t & operator=(const DC1D_xmm4c4_t &) = delete;
                      
                     inline DC1D_xmm4c4_t & operator=(DC1D_xmm4c4_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           if(this->ismmap) 
                              gms_unmap<xmm4c4_t>(this->mx,this->mnx);
                           else 
                             gms_mm_free(this->mx); 
                                                                                                                         
                           this->mnx   = rhs.mnx;
                           this->mx    = &rhs.mx[0];
                           rhs.mnx     = 0ULL;
                           rhs.mx      = NULL;
                           return (*this);
                      }
                      
                      inline void allocate() {
                          using namespace gms::common;
                          this->mx  = (xmm4c4_t*)
                                         gms_mm_malloc( sizeof(xmm4c4_t)*this->mnx,64ULL);
                         
                    }
                      
                     
                    
              };
                   
              
                struct __ATTR_ALIGN__(64) DC1D_m128_t {
                     
                      __m128 * __restrict              mx
                      std::size_t                      mnx;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,31)
#endif
                     inline DC1D_m128_t()  {
                          
                          this->mnx  = 0ULL;
                          this->mx  = NULL;
                         
                    } 
                          
                     inline explicit DC1D_m128_t(const std::size_t _mnx) {
                                 
                          using namespace gms::common;
                          this->mnx = _mnx;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline DC1D_m128_t(const std::size_t _mnx,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->mnx = _mnx;
                             switch (fsize) {
                                 case:0
                                      this->mx = (__m128*)
                                                 gms_mmap_4KiB<__m128>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mx = (__m128*)
                                                 gms_mmap_2MiB<__m128>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mx = (__m128*)
                                                 gms_mmap_1GiB<__m128>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                                                 
                            
                      
                     inline DC1D_m128_t(const std::size_t _mnx,
                                          const  __m128 * __restrict m_x) {
                                     
                          using namespace gms::common;
                          this->mnx = _mnx;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_DYN_CONTAINERS_SSE_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                 
	                  
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                 
	               
#endif
                   }  
                   
                    inline  DC1D_m128_t(DC1D_m128_t && rhs) {
                          
                          this->mnx = rhs.mnx;
                          this->mx   = &rhs.mx[0];
                          rhs.mnx    = 0ULL;
                          rhs.mx     = NULL;
                                                 
                      }
                                 
                     DC1D_m128_t(const DC1D_m128_t &)             = delete;
                      
                     inline ~DC1D_m128_t() {
                      
                          using namespace gms::common;
                          if(this->issmap) 
                             gms_unmap<__m128>(this->mx,this->mnx);
                          else 
                              gms_mm_free(this->mx);   
                     }        
                          
                                                 
                                  
                     DC1D_m128_t & operator=(const DC1D_m128_t &) = delete;
                      
                     inline DC1D_m128_t & operator=(DC1D_m128_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           if(this->ismmap) 
                              gms_unmap<__m128>(this->mx,this->mnx);
                           else 
                             gms_mm_free(this->mx); 
                                                                                                                         
                           this->mnx   = rhs.mnx;
                           this->mx    = &rhs.mx[0];
                           rhs.mnx     = 0ULL;
                           rhs.mx      = NULL;
                           return (*this);
                      }
                      
                      inline void allocate() {
                          using namespace gms::common;
                          this->mx  = (__m128*)
                                         gms_mm_malloc( sizeof(__m128)*this->mnx,64ULL);
                         
                    }
                      
                     
                    
              };
              
              
              
                 struct __ATTR_ALIGN__(64) DC1D_m128d_t {
                     
                      __m128d * __restrict              mx
                      std::size_t                      mnx;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,31)
#endif
                     inline DC1D_m128d_t()  {
                          
                          this->mnx  = 0ULL;
                          this->mx  = NULL;
                         
                    } 
                          
                     inline explicit DC1D_m128d_t(const std::size_t _mnx) {
                                 
                          using namespace gms::common;
                          this->mnx = _mnx;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline DC1D_m128d_t(const std::size_t _mnx,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->mnx = _mnx;
                             switch (fsize) {
                                 case:0
                                      this->mx = (__m128d*)
                                                 gms_mmap_4KiB<__m128d>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mx = (__m128d*)
                                                 gms_mmap_2MiB<__m128d>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mx = (__m128d*)
                                                 gms_mmap_1GiB<__m128d>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                                                 
                            
                      
                     inline DC1D_m128d_t(const std::size_t _mnx,
                                          const  __m128d * __restrict m_x) {
                                     
                          using namespace gms::common;
                          this->mnx = _mnx;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_DYN_CONTAINERS_SSE_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                 
	                  
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                 
	               
#endif
                   }  
                   
                    inline  DC1D_m128d_t(DC1D_m128d_t && rhs) {
                          
                          this->mnx = rhs.mnx;
                          this->mx   = &rhs.mx[0];
                          rhs.mnx    = 0ULL;
                          rhs.mx     = NULL;
                                                 
                      }
                                 
                     DC1D_m128d_t(const DC1D_m128d_t &)             = delete;
                      
                     inline ~DC1D_m128d_t() {
                      
                          using namespace gms::common;
                          if(this->issmap) 
                             gms_unmap<__m128d>(this->mx,this->mnx);
                          else 
                              gms_mm_free(this->mx);   
                     }        
                          
                                                 
                                  
                     DC1D_m128d_t & operator=(const DC1D_m128d_t &) = delete;
                      
                     inline DC1D_m128d_t & operator=(DC1D_m128d_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           if(this->ismmap) 
                              gms_unmap<__m128d>(this->mx,this->mnx);
                           else 
                             gms_mm_free(this->mx); 
                                                                                                                         
                           this->mnx   = rhs.mnx;
                           this->mx    = &rhs.mx[0];
                           rhs.mnx     = 0ULL;
                           rhs.mx      = NULL;
                           return (*this);
                      }
                      
                      inline void allocate() {
                          using namespace gms::common;
                          this->mx  = (__m128d*)
                                         gms_mm_malloc( sizeof(__m128d)*this->mnx,64ULL);
                         
                    }
                      
                     
                    
              };
              
               
            
              
              
           
           
           
           
              
              

              
              




} //gms

















#endif /*__GMS_DYN_CONTAINERS_SSE_HPP__*/
