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

#ifndef __GMS_DYN_CONTAINERS_HPP__
#define __GMS_DYN_CONTAINERS_HPP__ 051220230457


namespace file_info {

     static const unsigned int GMS_DYN_CONTAINERS_MAJOR = 1;
     static const unsigned int GMS_DYN_CONTAINERS_MINOR = 1;
     static const unsigned int GMS_DYN_CONTAINERS_MICRO = 0;
     static const unsigned int GMS_DYN_CONTAINERS_FULLVER =
       1000U*GMS_DYN_CONTAINERS_MAJOR+100U*GMS_DYN_CONTAINERS_MINOR+
       10U*GMS_DYN_CONTAINERS_MICRO;
     static const char GMS_DYN_CONTAINERS_CREATION_DATE[] = "05-12-2023 04:57 +00200 (TUE 05 DEC 2023 GMT+2)";
     static const char GMS_DYN_CONTAINERS_BUILD_DATE[]    = __DATE__; 
     static const char GMS_DYN_CONTAINERS_BUILD_TIME[]  = __TIME__;
     static const char GMS_DYN_CONTAINERS_SYNOPSIS[]    = "Dynamically allocated, 64-byte aligned dense containers.";

}





#include <cstdint>
#include <complex>
#include <vector>
#include <exception> //std::terminate
#include <valarray>
#include <cstring> // std::memcpy
#include <iostream>
#include "GMS_config.h"
#include "GMS_malloc.h"


// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (USE_GMS_DYN_CONTAINERS_NT_STORES)
#define USE_GMS_DYN_CONTAINERS_NT_STORES 0
#endif


namespace gms {


           struct alignas(64) DC3D_c4_t final 
           {
                      // Complex electric field.
                      std::complex<float> * __restrict m_Ex;
                      std::complex<float> * __restrict m_Ey;
                      std::complex<float> * __restrict m_Ez;
                      std::size_t                      mnx;
                      std::size_t                      mny;
                      std::size_t                      mnz;
                      bool                             ismmap;

                     inline DC3D_c4_t() noexcept(true) 
                     {
                          
                          this->mnx  = 0ULL;
                          this->mny  = 0ULL;
                          this->mnz  = 0ULL;
                          this->m_Ex  = NULL;
                          this->m_Ey  = NULL;
                          this->m_Ez  = NULL;
                      } 
                          
                     inline  DC3D_c4_t(const std::size_t nx,
                                       const std::size_t ny,
                                       const std::size_t nz) noexcept(false)
                     {
                          using namespace gms::common;
                          this->mnx = nx;
                          this->mny = ny;
                          this->mnz = nz;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline  DC3D_c4_t(const std::size_t nx,
                                       const std::size_t ny,
                                       const std::size_t nz,
                                       const int32_t prot,
                                       const int32_t flags,
                                       const int32_t fd,
                                       const long offset,
                                       const int32_t fsize) noexcept(false)
                     {
                             using namespace gms::common;
                             this->mnx = nx;
                             this->mny = ny;
                             this->mnz = nz;
                             switch (fsize) {
                                 case 0:
                                      this->m_Ex = (std::complex<float>*)
                                                 gms_mmap_4KiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->m_Ey = (std::complex<float>*)
                                                 gms_mmap_4KiB<std::complex<float>>(this->mny,prot,flags,fd,offset);
                                      this->m_Ez = (std::complex<float>*)
                                                 gms_mmap_4KiB<std::complex<float>>(this->mnz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 1:
                                      this->m_Ex = (std::complex<float>*)
                                                 gms_mmap_2MiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->m_Ey = (std::complex<float>*)
                                                 gms_mmap_2MiB<std::complex<float>>(this->mny,prot,flags,fd,offset);
                                      this->m_Ez = (std::complex<float>*)
                                                 gms_mmap_2MiB<std::complex<float>>(this->mnz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 2:
                                      this->m_Ex = (std::complex<float>*)
                                                 gms_mmap_1GiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->m_Ey = (std::complex<float>*)
                                                 gms_mmap_1GiB<std::complex<float>>(this->mny,prot,flags,fd,offset);
                                      this->m_Ez = (std::complex<float>*)
                                                 gms_mmap_1GiB<std::complex<float>>(this->mnz,prot,flags,fd,offset); 
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                     inline DC3D_c4_t(const std::vector<std::complex<float>> &Ex,    //shall be of the same size (no error checking implemented)
                                      const std::vector<std::complex<float>> &Ey,    //shall be of the same size (no error checking implemented)
                                      const std::vector<std::complex<float>> &Ez) noexcept(false) {  //shall be of the same size (no error checking implemented)
                          using namespace gms::common;
                          
                          this->mnx = Ex.size();
                          this->mny = Ey.size();
                          this->mnz = Ez.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          const std::size_t leny{bytes_mny()};
                          const std::size_t lenz{bytes_mnz()};
                          std::memcpy(this->m_Ex,&Ex[0],lenx);
                          std::memcpy(this->m_Ey,&Ey[0],leny);
                          std::memcpy(this->m_Ez,&Ez[0],lenz);       
                     }
                     
                     inline DC3D_c4_t(const std::valarray<std::complex<float>> &Ex,
                                      const std::valarray<std::complex<float>> &Ey,
                                      const std::valarray<std::complex<float>> &Ez) noexcept(false) {
                          using namespace gms::common;
                          
                          this->mnx = Ex.size();
                          this->mny = Ey.size();
                          this->mnz = Ez.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          const std::size_t leny{bytes_mny()};
                          const std::size_t lenz{bytes_mnz()};
                          std::memcpy(this->m_Ex,&Ex[0],lenx);
                          std::memcpy(this->m_Ey,&Ey[0],leny);
                          std::memcpy(this->m_Ez,&Ez[0],lenz);           
                     }
                             
                                           
                     inline explicit DC3D_c4_t(const std::size_t nx,
                                               const std::size_t ny,
                                               const std::size_t nz,
                                               const std::complex<float> * __restrict Ex,
                                               const std::complex<float> * __restrict Ey,
                                               const std::complex<float> * __restrict Ez) noexcept(false) {
                          using namespace gms::common;
                          this->mnx = nx;
                          this->mny = ny;
                          this->mnz = nz;
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          const std::size_t leny{bytes_mny()};
                          const std::size_t lenz{bytes_mnz()};
                          std::memcpy(this->m_Ex,&Ex[0],lenx);
                          std::memcpy(this->m_Ey,&Ey[0],leny);
                          std::memcpy(this->m_Ez,&Ez[0],lenz);  

                   }  
                   
                    inline  DC3D_c4_t(DC3D_c4_t && rhs) noexcept(true) {
                          
                          this->mnx    = rhs.mnx;
                          this->mny    = rhs.mny;
                          this->mnz    = rhs.mnz;
                          this->ismmap = rhs.ismmap;
                          this->m_Ex   = &rhs.m_Ex[0];
                          this->m_Ey   = &rhs.m_Ey[0];
                          this->m_Ez   = &rhs.m_Ez[0];
                          rhs.mnx      = 0ULL;
                          rhs.mny      = 0ULL;
                          rhs.mnz      = 0ULL;
                          rhs.m_Ex     = NULL;
                          rhs.m_Ey     = NULL;
                          rhs.m_Ez     = NULL;
                      }
                                 
                     DC3D_c4_t(const DC3D_c4_t &)             = delete;
                      
                     inline ~DC3D_c4_t() noexcept(false) {
                      
                          using namespace gms::common;
                          if(this->ismmap) {
                              int32_t err1{}, err2{}, err3{};
                             err1 = gms_ummap<std::complex<float>>(this->m_Ex,this->mnx); this->m_Ex = NULL;
                             err2 = gms_ummap<std::complex<float>>(this->m_Ey,this->mny); this->m_Ey = NULL;
                             err3 = gms_ummap<std::complex<float>>(this->m_Ez,this->mnz); this->m_Ez = NULL;
                             if(__builtin_expect(err1==-1,0) || 
                                __builtin_expect(err2==-1,0) || 
                                __builtin_expect(err3==-1,0))
                             {
#if (FAST_TERMINATE_ON_CRITICAL_ERROR) == 1
                                __builtin_trap();
#else
                                std::terminate();
#endif
                             }
                          }
                          else {
                              gms_mm_free(this->m_Ex); this->m_Ex = NULL;
                              gms_mm_free(this->m_Ey); this->m_Ey = NULL;
                              gms_mm_free(this->m_Ez); this->m_Ez = NULL;
                          }
                      }
                      
                     DC3D_c4_t & operator=(const DC3D_c4_t &) = delete;
                      
                     inline DC3D_c4_t & operator=(DC3D_c4_t &&rhs) noexcept(true)
                     {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           
                           gms_swap(this->mnx,rhs.mnx);
                           gms_swap(this->mny,rhs.mny);
                           gms_swap(this->mnz,rhs.mnz);
                           gms_swap(this->ismmap,rhs.ismmap);
                           gms_swap(this->m_Ex,rhs.m_Ex);
                           gms_swap(this->m_Ey,rhs.m_Ey);
                           gms_swap(this->m_Ez,rhs.m_Ez);
                           rhs.mnx       = 0ULL;
                           rhs.mny       = 0ULL;
                           rhs.mnz       = 0ULL;
                           rhs.m_Ex      = NULL;
                           rhs.m_Ey      = NULL;
                           rhs.m_Ez      = NULL;
                           return (*this);
                      }

                    inline std::size_t bytes_mnx() const noexcept(true)
                    {
                         return (sizeof(std::complex<float>)*this->mnx);
                    }

                    inline std::size_t bytes_mny() const noexcept(true)
                    {
                         return (sizeof(std::complex<float>)*this->mny);
                    }

                    inline std::size_t bytes_mnz() const noexcept(true)
                    {
                         return (sizeof(std::complex<float>)*this->mnz);
                    }

                      private:
                      inline void allocate() noexcept(true)
                      {
                          using namespace gms::common;
                          this->m_Ex  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->mnx,64ULL);
                          this->m_Ey  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->mny,64ULL);
                          this->m_Ez  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->mny,64ULL);
                      }
                     public:
                     inline void info_size_alignment() const 
                     {
                         std::cout << "alignof(struct DC3D_c4_t) = " << alignof(DC3D_c4_t) << '\n'
                                  << "sizeof(struct  DC3D_c4_t) = " << sizeof(DC3D_c4_t)  << '\n'
                                  << std::hex << std::showbase      << '\n'
                                  << "&this->m_Ex               =" << (void*)&this->m_Ex  << "\n"
                                     "&this->m_Ey               =" << (void*)&this->m_Ey  << "\n"
                                     "&this->m_Ez               =" << (void*)&this->m_Ez  << "\n"
                                     "&this->mnx                =" << &this->mnx   << "\n"
                                     "&this->mny                =" << &this->mny   << "\n"
                                     "&this->mnz                =" << &this->mnz   << "\n"
                                     "&this->ismmap             =" << &this->ismmap << "\n";

                     }
                    
              };
              
              
               struct alignas(64) DC2D_c4_t final 
               {
                      // Complex electric field.
                      std::complex<float> * __restrict m_Ex;
                      std::complex<float> * __restrict m_Ey;
                      std::size_t                      mnx;
                      std::size_t                      mny;
                      bool                             ismmap;

                     inline DC2D_c4_t() noexcept(true) 
                     {
                          
                          this->mnx   = 0ULL;
                          this->mny   = 0ULL;
                          this->m_Ex  = NULL;
                          this->m_Ey  = NULL;
                         
                      } 
                          
                     inline  DC2D_c4_t(const std::size_t nx,
                                       const std::size_t ny) noexcept(false)
                    {
                          using namespace gms::common;
                          this->mnx = nx;
                          this->mny = ny;
                          allocate();
                          this->ismmap = false;
                    }  
                      
                     inline  DC2D_c4_t(const std::size_t nx,
                                       const std::size_t ny,
                                       const int32_t prot,
                                       const int32_t flags,
                                       const int32_t fd,
                                       const long offset,
                                       const int32_t fsize) noexcept(false)
                    {
                             using namespace gms::common;
                             this->mnx = nx;
                             this->mny = ny;
                             switch (fsize) {
                                 case 0:
                                      this->m_Ex = (std::complex<float>*)
                                                 gms_mmap_4KiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->m_Ey = (std::complex<float>*)
                                                 gms_mmap_4KiB<std::complex<float>>(this->mny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 1:
                                      this->m_Ex = (std::complex<float>*)
                                                 gms_mmap_2MiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->m_Ey = (std::complex<float>*)
                                                 gms_mmap_2MiB<std::complex<float>>(this->mny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 2:
                                      this->m_Ex = (std::complex<float>*)
                                                 gms_mmap_1GiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->m_Ey = (std::complex<float>*)
                                                 gms_mmap_1GiB<std::complex<float>>(this->mny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                     inline DC2D_c4_t(const std::vector<std::complex<float>> &Ex,    //shall be of the same size (no error checking implemented)
                                      const std::vector<std::complex<float>> &Ey) noexcept(false)
                    {    //shall be of the same size (no error checking implemented)
                                  
                          using namespace gms::common;
                          this->mnx = Ex.size();
                          this->mny = Ey.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          const std::size_t leny{bytes_mny()};
                          std::memcpy(this->m_Ex,&Ex[0],lenx);
                          std::memcpy(this->m_Ey,&Ey[0],leny);
                        
                     }
                     
                     inline DC2D_c4_t(const std::valarray<std::complex<float>> &Ex,
                                      const std::valarray<std::complex<float>> &Ey) noexcept(false)
                     {
                                 
                          using namespace gms::common;
                          this->mnx = Ex.size();
                          this->mny = Ey.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          const std::size_t leny{bytes_mny()};
                          std::memcpy(this->m_Ex,&Ex[0],lenx);
                          std::memcpy(this->m_Ey,&Ey[0],leny);
                             
                     }
                             
                             
                      
                     inline explicit DC2D_c4_t(const std::size_t nx,
                                               const std::size_t ny,
                                               const std::complex<float> * __restrict Ex,
                                               const std::complex<float> * __restrict Ey) noexcept(false)
                     {
                                
                          using namespace gms::common;
                          this->mnx = nx;
                          this->mny = ny;
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          const std::size_t leny{bytes_mny()};
                          std::memcpy(this->m_Ex,&Ex[0],lenx);
                          std::memcpy(this->m_Ey,&Ey[0],leny);
                   }  
                   
                    inline  DC2D_c4_t(DC2D_c4_t && rhs) noexcept(true)
                    {
                          
                          this->mnx  = rhs.mnx;
                          this->mny  = rhs.mny;
                          this->m_Ex = &rhs.m_Ex[0];
                          this->m_Ey = &rhs.m_Ey[0];
                          rhs.mnx    = 0ULL;
                          rhs.mny    = 0ULL;
                          rhs.m_Ex     = NULL;
                          rhs.m_Ey     = NULL;
                         
                      }
                                 
                     DC2D_c4_t(const DC2D_c4_t &)             = delete;
                      
                     inline ~DC2D_c4_t() noexcept(false)
                     {
                      
                          using namespace gms::common;
                          if(this->ismmap) 
                          {
                             int32_t err1{}, err2{};
                             err1 = gms_ummap<std::complex<float>>(this->m_Ex,this->mnx); this->m_Ex = NULL;
                             err2 = gms_ummap<std::complex<float>>(this->m_Ey,this->mny); this->m_Ey = NULL;
                             if(__builtin_expect(err1==-1,0) || 
                                __builtin_expect(err2==-1,0))
                             {
#if (FAST_TERMINATE_ON_CRITICAL_ERROR) == 1
                                __builtin_trap();
#else
                                std::terminate();
#endif                                 
                             }
                          }
                          else 
                          {
                              gms_mm_free(this->m_Ex); this->m_Ex = NULL;
                              gms_mm_free(this->m_Ey); this->m_Ey = NULL; 
                             
                          }
                      }
                      
                     DC2D_c4_t & operator=(const DC2D_c4_t &) = delete;
                      
                     inline DC2D_c4_t & operator=(DC2D_c4_t &&rhs) noexcept(true)
                      {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           
                           gms_swap(this->mnx,rhs.mnx);
                           gms_swap(this->mny,rhs.mny);
                           gms_swap(this->m_Ex,rhs.m_Ex);
                           gms_swap(this->m_Ey,rhs.m_Ey);
                           rhs.mnx       = 0ULL;
                           rhs.mny       = 0ULL;
                           rhs.m_Ex      = NULL;
                           rhs.m_Ey      = NULL;
                           return (*this);
                      }

                    inline std::size_t bytes_mnx() const noexcept(true)
                    {
                         return (sizeof(std::complex<float>)*this->mnx);
                    }

                    inline std::size_t bytes_mny() const noexcept(true)
                    {
                         return (sizeof(std::complex<float>)*this->mny);
                    }
                    
                    private:
                      inline void allocate() 
                      {
                          using namespace gms::common;
                          this->m_Ex  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->mnx,64ULL);
                          this->m_Ey  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->mny,64ULL);
                     }

                     public:
                     inline void info_size_alignment() const 
                     {
                          std::cout << "alignof(struct DC2D_c4_t) = " << alignof(DC2D_c4_t) << '\n'
                                   << "sizeof(struct  DC2D_c4_t) = " << sizeof(DC2D_c4_t)  << '\n'
                                   << std::hex << std::showbase      << '\n'
                                   << "&this->m_Ex               ="  << (void*)&this->m_Ex  << "\n"
                                      "&this->m_Ey               ="  << (void*)&this->m_Ey  << "\n"
                                      "&this->mnx                ="  << &this->mnx   << "\n"
                                      "&this->mny                ="  << &this->mny   << "\n"
                                      "&this->ismmap             ="  << &this->ismmap << "\n";
                     }
                      
               };
               
               
               struct alignas(64) DC1D_c4_t final 
               {
                      // Complex electric field.
                      std::complex<float> * __restrict m_Ex;
                      std::size_t                      mnx;
                      bool                             ismmap;

                     inline DC1D_c4_t() noexcept(true)
                     {
                          
                          this->mnx   = 0ULL;
                          this->m_Ex  = NULL;
                                      
                      } 
                          
                     inline explicit DC1D_c4_t(const std::size_t nx) noexcept(false)
                     {
                          using namespace gms::common;
                          this->mnx = nx;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline  DC1D_c4_t(const std::size_t nx,
                                       const int32_t prot,
                                       const int32_t flags,
                                       const int32_t fd,
                                       const long offset,
                                       const int32_t fsize) noexcept(false)
                     {
                             using namespace gms::common;
                             this->mnx = nx;
                             switch (fsize) {
                                 case 0:
                                      this->m_Ex = (std::complex<float>*)
                                                 gms_mmap_4KiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 1:
                                      this->m_Ex = (std::complex<float>*)
                                                 gms_mmap_2MiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 2:
                                      this->m_Ex = (std::complex<float>*)
                                                 gms_mmap_1GiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                     inline DC1D_c4_t(const std::vector<std::complex<float>> &Ex) noexcept(false)
                     {    //shall be of the same size (no error checking implemented)
                                  
                          using namespace gms::common;
                          this->mnx = Ex.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_Ex,&Ex[0],lenx);
                                                
                     }
                     
                     inline DC1D_c4_t(const std::valarray<std::complex<float>> &Ex) noexcept(false)
                    {
                                   
                          using namespace gms::common;
                          this->mnx = Ex.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_Ex,&Ex[0],lenx);
                                                    
                     }
                             
                             
                      
                     inline  DC1D_c4_t(const std::size_t nx,
                                       const std::complex<float> * __restrict Ex) noexcept(false)
                     {
                                 
                          using namespace gms::common;
                          this->mnx = nx;
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<float>)*this->mnx;
                          std::memcpy(this->m_Ex,&Ex[0],lenx);
                   }  
                   
                    inline  DC1D_c4_t(DC1D_c4_t && rhs) noexcept(true)
                    {
                          
                          this->mnx  = rhs.mnx;
                          this->m_Ex = &rhs.m_Ex[0];
                          rhs.mnx    = 0ULL;
                          rhs.m_Ex   = NULL;
                                                
                      }
                                 
                     DC1D_c4_t(const DC1D_c4_t &)             = delete;
                      
                     inline ~DC1D_c4_t() noexcept(false)
                     {
                          using namespace gms::common;
                          if(this->ismmap) 
                          { 
                             int32_t err1{};
                             err1 = gms_ummap<std::complex<float>>(this->m_Ex,this->mnx); this->m_Ex = NULL;
                             if(__builtin_expect(err1==-1,0))
                             {
#if (FAST_TERMINATE_ON_CRITICAL_ERROR) == 1
                                __builtin_trap();
#else
                                std::terminate();
#endif                               
                             }
                          }
                          else 
                          {
                             gms_mm_free(this->m_Ex); this->m_Ex = NULL;  
                          }  
                     }    
                         
                                                                 
                     DC1D_c4_t & operator=(const DC1D_c4_t &) = delete;
                      
                     inline DC1D_c4_t & operator=(DC1D_c4_t &&rhs) noexcept(true)
                     {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                                                                                                                
                           gms_swap(this->mnx,rhs.mnx);
                           gms_swap(this->m_Ex,rhs.m_Ex);
                           rhs.mnx     = 0ULL;
                           rhs.m_Ex    = NULL;
                           return (*this);
                      }

                    inline std::size_t bytes_mnx() const noexcept(true)
                    {
                         return (sizeof(std::complex<float>)*this->mnx);
                    }

                      private:
                      inline void allocate() noexcept(false)
                      {
                          using namespace gms::common;
                          this->m_Ex  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->mnx,64ULL);
                       
                     }
                     public: 
                     inline void info_size_alignment() const 
                     {
                          std::cout << "alignof(struct DC1D_c4_t) = " << alignof(DC1D_c4_t) << '\n'
                                   << "sizeof(struct  DC1D_c4_t) = " << sizeof(DC1D_c4_t)  << '\n'
                                   << std::hex << std::showbase      << '\n'
                                   << "&this->m_Ex               ="  << (void*)&this->m_Ex  << "\n"
                                      "&this->mnx                ="  << &this->mnx   << "\n"
                                      "&this->ismmap             ="  << &this->ismmap << "\n";
                     }
                      
               };
               
               
               struct alignas(64) DC3D_r4_t final 
               {
                     
                      float * __restrict m_Exr;
                      float * __restrict m_Exi;
                      float * __restrict m_Eyr;
                      float * __restrict m_Eyi;
                      float * __restrict m_Ezr;
                      float * __restrict m_Ezi;
                      std::size_t        mnx;
                      std::size_t        mny;
                      std::size_t        mnz;
                      bool               ismmap;

                      inline DC3D_r4_t() noexcept(true)
                      {
                         this->mnx    = 0ULL;
                         this->mny    = 0ULL;
                         this->mnz    = 0ULL;
                         this->m_Exr  = NULL;
                         this->m_Exi  = NULL;
                         this->m_Eyr  = NULL;
                         this->m_Eyi  = NULL;
                         this->m_Ezr  = NULL;
                         this->m_Ezi  = NULL;
                      }                    
                      
                      inline  DC3D_r4_t(const std::size_t nx,
                                        const std::size_t ny,
                                        const std::size_t nz) noexcept(false)
                     {
                             this->mnx = nx;
                             this->mny = ny;
                             this->mnz = nz;
                             allocate();
                             this->ismmap = false;
                      }  
                      
                      inline  DC3D_r4_t(const std::size_t nx,
                                        const std::size_t ny,
                                        const std::size_t nz,
                                        const int32_t prot,
                                        const int32_t flags,
                                        const int32_t fd,
                                        const long offset,
                                        const int32_t fsize) noexcept(false)
                      {
                             using namespace gms::common;
                             this->mnx = nx;
                             this->mny = ny;
                             this->mnz = nz;
                             switch (fsize) {
                                 case 0:
                                      this->m_Exr = (float*)
                                                  gms_mmap_4KiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->m_Exi = (float*)
                                                  gms_mmap_4KiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->m_Eyr = (float*)
                                                  gms_mmap_4KiB<float>(this->mny,prot,flags,fd,offset);
                                      this->m_Eyi = (float*)
                                                  gms_mmap_4KiB<float>(this->mny,prot,flags,fd,offset);
                                      this->m_Ezr = (float*)
                                                  gms_mmap_4KiB<float>(this->mnz,prot,flags,fd,offset);
                                      this->m_Ezi = (float*)
                                                  gms_mmap_4KiB<float>(this->mnz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 1:
                                      this->m_Exr = (float*)
                                                  gms_mmap_2MiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->m_Exi = (float*)
                                                  gms_mmap_2MiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->m_Eyr = (float*)
                                                  gms_mmap_2MiB<float>(this->mny,prot,flags,fd,offset);
                                      this->m_Eyi = (float*)
                                                  gms_mmap_2MiB<float>(this->mny,prot,flags,fd,offset);
                                      this->m_Ezr = (float*)
                                                  gms_mmap_2MiB<float>(this->mnz,prot,flags,fd,offset);
                                      this->m_Ezi = (float*)
                                                  gms_mmap_2MiB<float>(this->mnz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 2:
                                      this->m_Exr = (float*)
                                                  gms_mmap_1GiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->m_Exi = (float*)
                                                  gms_mmap_1GiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->m_Eyr = (float*)
                                                  gms_mmap_1GiB<float>(this->mny,prot,flags,fd,offset);
                                      this->m_Eyi = (float*)
                                                  gms_mmap_1GiB<float>(this->mny,prot,flags,fd,offset);
                                      this->m_Ezr = (float*)
                                                  gms_mmap_1GiB<float>(this->mnz,prot,flags,fd,offset);
                                      this->m_Ezi = (float*)
                                                  gms_mmap_1GiB<float>(this->mnz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     } 
                       
                      //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline DC3D_r4_t(const std::vector<float> &Exr,
                                       const std::vector<float> &Exi,
                                       const std::vector<float> &Eyr,
                                       const std::vector<float> &Eyi,
                                       const std::vector<float> &Ezr,
                                       const std::vector<float> &Ezi) noexcept(false)
                      {
                               
                               this->mnx = Exr.size(); 
                               this->mny = Eyr.size();
                               this->mnz = Ezr.size();
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx{bytes_mnx()};
                               const std::size_t leny{bytes_mny()};
                               const std::size_t lenz{bytes_mnz()};
                               std::memcpy(this->m_Exr,&Exr[0],lenx);
                               std::memcpy(this->m_Exi,&Exi[0],lenx);
                               std::memcpy(this->m_Eyr,&Eyr[0],leny);
                               std::memcpy(this->m_Eyi,&Eyi[0],leny);
                               std::memcpy(this->m_Ezr,&Ezr[0],lenz);
                               std::memcpy(this->m_Ezi,&Ezi[0],lenz);     
                      }
                      
                      
                       //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline DC3D_r4_t(const std::valarray<float> &Exr,
                                       const std::valarray<float> &Exi,
                                       const std::valarray<float> &Eyr,
                                       const std::valarray<float> &Eyi,
                                       const std::valarray<float> &Ezr,
                                       const std::valarray<float> &Ezi) noexcept(false) 
                      {
                               
                               this->mnx = Exr.size(); 
                               this->mny = Eyr.size();
                               this->mnz = Ezr.size();
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(float)*this->mnx;
                               const std::size_t leny = sizeof(float)*this->mny;
                               const std::size_t lenz = sizeof(float)*this->mnz;
                               std::memcpy(this->m_Exr,&Exr[0],lenx);
                               std::memcpy(this->m_Exi,&Exi[0],lenx);
                               std::memcpy(this->m_Eyr,&Eyr[0],leny);
                               std::memcpy(this->m_Eyi,&Eyi[0],leny);
                               std::memcpy(this->m_Ezr,&Ezr[0],lenz);
                               std::memcpy(this->m_Ezi,&Ezi[0],lenz);     
                      }
                      
                      inline explicit DC3D_r4_t(const std::size_t nx,
                                                const std::size_t ny,
                                                const std::size_t nz,
                                                const float * __restrict Exr,   
                                                const float * __restrict Exi,
                                                const float * __restrict Eyr,
                                                const float * __restrict Eyi,
                                                const float * __restrict Ezr,
                                                const float * __restrict Ezi) noexcept(false)
                     {
                          using namespace gms::common;
                          this->mnx = nx;
                          this->mny = ny;
                          this->mnz = nz;
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          const std::size_t leny{bytes_mny()};
                          const std::size_t lenz{bytes_mnz()};
                          std::memcpy(this->m_Exr,&Exr[0],lenx);
                          std::memcpy(this->m_Exi,&Exi[0],lenx);
                          std::memcpy(this->m_Eyr,&Eyr[0],leny);
                          std::memcpy(this->m_Eyi,&Eyi[0],leny);
                          std::memcpy(this->m_Ezr,&Ezr[0],lenz);
                          std::memcpy(this->m_Ezi,&Ezi[0],lenz);  
                      }
                      
                      inline DC3D_r4_t(DC3D_r4_t &&rhs) noexcept(true)
                       {
                           
                          this->mnx    = rhs.mnx;
                          this->mny    = rhs.mny;
                          this->mnz    = rhs.mnz;
                          this->m_Exr  = &rhs.m_Exr[0];
                          this->m_Exi  = &rhs.m_Exi[0];
                          this->m_Eyr  = &rhs.m_Eyr[0];
                          this->m_Eyi  = &rhs.m_Eyi[0];
                          this->m_Ezr  = &rhs.m_Ezr[0];
                          this->m_Ezi  = &rhs.m_Ezi[0];
                          rhs.mnx      = 0ULL;
                          rhs.mny      = 0ULL;
                          rhs.mnz      = 0ULL;
                          rhs.m_Exr    = NULL;
                          rhs.m_Exi    = NULL;
                          rhs.m_Eyr    = NULL;
                          rhs.m_Eyi    = NULL;
                          rhs.m_Ezr    = NULL;
                          rhs.m_Ezi    = NULL;
                      }   
                        
                      DC3D_r4_t(const DC3D_r4_t &)             = delete;
                      
                      inline ~DC3D_r4_t() noexcept(false)
                      {
                           using namespace gms::common;
                           if(this->ismmap) {
                               int32_t err1{}, err2{}, err3{},
                                               err4{}, err5{}, err6{};
                              err1 = gms_ummap<float>(this->m_Exr,this->mnx);this->m_Exr = NULL;
                              err2 = gms_ummap<float>(this->m_Exi,this->mnx);this->m_Exi = NULL;
                              err3 = gms_ummap<float>(this->m_Eyr,this->mny);this->m_Eyr = NULL;
                              err4 = gms_ummap<float>(this->m_Eyi,this->mny);this->m_Eyi = NULL;
                              err4 = gms_ummap<float>(this->m_Ezr,this->mnz);this->m_Ezr = NULL;
                              err5 = gms_ummap<float>(this->m_Ezi,this->mnz);this->m_Ezi = NULL;
                              if(__builtin_expect(err1==-1,0) || 
                                 __builtin_expect(err2==-1,0) || 
                                 __builtin_expect(err3==-1,0) || 
                                 __builtin_expect(err4==-1,0) || 
                                 __builtin_expect(err5==-1,0) || 
                                 __builtin_expect(err6==-1,0))
                              {
#if (FAST_TERMINATE_ON_CRITICAL_ERROR) == 1
                                __builtin_trap();
#else
                                std::terminate();
#endif  
                              }
                           }
                           else {
                               gms_mm_free(this->m_Exr);this->m_Exr = NULL;
                               gms_mm_free(this->m_Exi);this->m_Exi = NULL;
                               gms_mm_free(this->m_Eyr);this->m_Eyr = NULL;
                               gms_mm_free(this->m_Eyi);this->m_Eyi = NULL;
                               gms_mm_free(this->m_Ezr);this->m_Ezr = NULL;
                               gms_mm_free(this->m_Ezi);this->m_Ezi = NULL;
                           }
                      }
                      
                      DC3D_r4_t & operator=(const DC3D_r4_t &) = delete;
                      
                      inline DC3D_r4_t & operator=(DC3D_r4_t &&rhs) noexcept(true)
                      {
                            using namespace gms::common;
                            if(this==&rhs) return (*this);
                           
                            gms_swap(this->mnx,rhs.mnx);
                            gms_swap(this->mny,rhs.mny);
                            gms_swap(this->mnz,rhs.mnz);
                            gms_swap(this->m_Exr,rhs.m_Exr);
                            gms_swap(this->m_Exi,rhs.m_Exi);
                            gms_swap(this->m_Eyr,rhs.m_Eyr);
                            gms_swap(this->m_Eyi,rhs.m_Eyi);
                            gms_swap(this->m_Ezr,rhs.m_Ezr);
                            gms_swap(this->m_Ezi,rhs.m_Ezi);
                            rhs.mnx       = 0ULL;
                            rhs.mny       = 0ULL;
                            rhs.mnz       = 0ULL;
                            rhs.m_Exr     = NULL;
                            rhs.m_Exi     = NULL;
                            rhs.m_Eyr     = NULL;
                            rhs.m_Eyi     = NULL;
                            rhs.m_Ezr     = NULL;
                            rhs.m_Ezi     = NULL;
                            return (*this);
                      }

                   inline  std::size_t bytes_mnx() const noexcept(true)
                    {
                         return (sizeof(float)*this->mnx);
                    }

                   inline  std::size_t bytes_mny() const noexcept(true)
                    {
                         return (sizeof(float)*this->mny);
                    }

                   inline std::size_t bytes_mnz() const noexcept(true)
                    {
                         return (sizeof(float)*this->mnz);
                    }
                   
                   private:  
                   inline void allocate() noexcept(false)
                   {
                        using namespace gms::common;
                        this->m_Exr  = (float*)gms_mm_malloc(sizeof(float)*this->mnx,64ULL);
                        this->m_Exi  = (float*)gms_mm_malloc(sizeof(float)*this->mnx,64ULL);
                        this->m_Eyr  = (float*)gms_mm_malloc(sizeof(float)*this->mny,64ULL);
                        this->m_Eyi  = (float*)gms_mm_malloc(sizeof(float)*this->mny,64ULL);
                        this->m_Ezr  = (float*)gms_mm_malloc(sizeof(float)*this->mnz,64ULL);
                        this->m_Ezi  = (float*)gms_mm_malloc(sizeof(float)*this->mnz,64ULL);
                   }    
                   public:
                   inline void info_size_alignment() const
                   {
                        std::cout  << "alignof(struct DC3D_r4_t) = "  << alignof(DC3D_r4_t) << '\n'
                                  << "sizeof(struct  DC3D_r4_t) = "  << sizeof(DC3D_r4_t)  << '\n'
                                  << std::hex << std::showbase       << '\n'
                                  << "&this->m_Exr              ="  << (void*)&this->m_Exr << "\n"
                                     "&this->m_Exi              ="  << (void*)&this->m_Exi << "\n"
                                     "&this->m_Eyr              ="  << (void*)&this->m_Eyr << "\n"
                                     "&this->m_Eyi              ="  << (void*)&this->m_Eyi << "\n"
                                     "&this->m_Ezr              ="  << (void*)&this->m_Ezr << "\n"
                                     "&this->m_Ezi              ="  << (void*)&this->m_Ezi << "\n"
                                     "&this->mnx                ="  << &this->mnx   <<  "\n"
                                     "&this->mny                ="  << &this->mny   <<  "\n"
                                     "&this->mnz                ="  << &this->mnz   <<  "\n"
                                     "&this->ismmap             ="  << &this->ismmap << "\n";
                   }
                   
              };
              
              
             struct alignas(64) DC2D_r4_t final 
             {
                     
                      float * __restrict m_Exr;
                      float * __restrict m_Exi;
                      float * __restrict m_Eyr;
                      float * __restrict m_Eyi;
                      std::size_t        mnx;
                      std::size_t        mny;
                      bool              ismmap;

                      inline DC2D_r4_t() noexcept(true)
                      {
                      
                         this->mnx    = 0ULL;
                         this->mny    = 0ULL;
                         this->m_Exr  = NULL;
                         this->m_Exi  = NULL;
                         this->m_Eyr  = NULL;
                         this->m_Eyi  = NULL;
                    }                    
                      
                      inline  DC2D_r4_t(const std::size_t nx,
                                        const std::size_t ny) noexcept(false)
                     {
                                     
                             this->mnx = nx;
                             this->mny = ny;
                             allocate();
                             this->ismmap = false;
                      }  
                      
                      inline  DC2D_r4_t(const std::size_t nx,
                                        const std::size_t ny,
                                        const int32_t prot,
                                        const int32_t flags,
                                        const int32_t fd,
                                        const long offset,
                                        const int32_t fsize) noexcept(false)
                      {
                             using namespace gms::common;
                             this->mnx = nx;
                             this->mny = ny;
                             switch (fsize) {
                                 case 0:
                                      this->m_Exr = (float*)
                                                  gms_mmap_4KiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->m_Exi = (float*)
                                                  gms_mmap_4KiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->m_Eyr = (float*)
                                                  gms_mmap_4KiB<float>(this->mny,prot,flags,fd,offset);
                                      this->m_Eyi = (float*)
                                                  gms_mmap_4KiB<float>(this->mny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 1:
                                      this->m_Exr = (float*)
                                                  gms_mmap_2MiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->m_Exi = (float*)
                                                  gms_mmap_2MiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->m_Eyr = (float*)
                                                  gms_mmap_2MiB<float>(this->mny,prot,flags,fd,offset);
                                      this->m_Eyi = (float*)
                                                  gms_mmap_2MiB<float>(this->mny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 2:
                                      this->m_Exr = (float*)
                                                  gms_mmap_1GiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->m_Exi = (float*)
                                                  gms_mmap_1GiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->m_Eyr = (float*)
                                                  gms_mmap_1GiB<float>(this->mny,prot,flags,fd,offset);
                                      this->m_Eyi = (float*)
                                                  gms_mmap_1GiB<float>(this->mny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     } 
                       
                      //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline DC2D_r4_t(const std::vector<float> &Exr,
                                       const std::vector<float> &Exi,
                                       const std::vector<float> &Eyr,
                                       const std::vector<float> &Eyi) noexcept(false)
                      {
                                    
                               this->mnx = Exr.size(); 
                               this->mny = Eyr.size();
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx{bytes_mnx()};
                               const std::size_t leny{bytes_mny()};
                               std::memcpy(this->m_Exr,&Exr[0],lenx);
                               std::memcpy(this->m_Exi,&Exi[0],lenx);
                               std::memcpy(this->m_Eyr,&Eyr[0],leny);
                               std::memcpy(this->m_Eyi,&Eyi[0],leny);
                                
                      }
                      
                      
                       //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline DC2D_r4_t(const std::valarray<float> &Exr,
                                       const std::valarray<float> &Exi,
                                       const std::valarray<float> &Eyr,
                                       const std::valarray<float> &Eyi) noexcept(false)
                     {
                                   
                               this->mnx = Exr.size(); 
                               this->mny = Eyr.size();
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx{bytes_mnx()};
                               const std::size_t leny{bytes_mny()};
                               std::memcpy(this->m_Exr,&Exr[0],lenx);
                               std::memcpy(this->m_Exi,&Exi[0],lenx);
                               std::memcpy(this->m_Eyr,&Eyr[0],leny);
                               std::memcpy(this->m_Eyi,&Eyi[0],leny);
                                 
                      }
                      
                      inline  DC2D_r4_t(const std::size_t nx,
                                        const std::size_t ny,
                                        const float * __restrict Exr,   
                                        const float * __restrict Exi,
                                        const float * __restrict Eyr,
                                        const float * __restrict Eyi) noexcept(false)
                      {                                   
                          using namespace gms::common;
                          this->mnx = nx;
                          this->mny = ny;
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          const std::size_t leny{bytes_mny()};
                          std::memcpy(this->m_Exr,&Exr[0],lenx);
                          std::memcpy(this->m_Exi,&Exi[0],lenx);
                          std::memcpy(this->m_Eyr,&Eyr[0],leny);
                          std::memcpy(this->m_Eyi,&Eyi[0],leny);
                      }
                      
                      inline DC2D_r4_t(DC2D_r4_t &&rhs) noexcept(true)
                      {
                          this->mnx    = rhs.mnx;
                          this->mny    = rhs.mny;
                          this->m_Exr  = &rhs.m_Exr[0];
                          this->m_Exi  = &rhs.m_Exi[0];
                          this->m_Eyr  = &rhs.m_Eyr[0];
                          this->m_Eyi  = &rhs.m_Eyi[0];
                          rhs.mnx      = 0ULL;
                          rhs.mny      = 0ULL;
                          rhs.m_Exr    = NULL;
                          rhs.m_Exi    = NULL;
                          rhs.m_Eyr    = NULL;
                          rhs.m_Eyi    = NULL;
                          
                      }   
                        
                      DC2D_r4_t(const DC2D_r4_t &)             = delete;
                      
                      inline ~DC2D_r4_t() noexcept(false)
                      {
                           using namespace gms::common;
                           if(this->ismmap) {
                              int32_t  err1{}, err2{}, err3{}, err4{};
                              err1 = gms_ummap<float>(this->m_Exr,this->mnx);this->m_Exr = NULL;
                              err2 = gms_ummap<float>(this->m_Exi,this->mnx);this->m_Exi = NULL;
                              err3 = gms_ummap<float>(this->m_Eyr,this->mny);this->m_Eyr = NULL;
                              err4 = gms_ummap<float>(this->m_Eyi,this->mny);this->m_Eyi = NULL;
                              if(__builtin_expect(err1==-1,0) || 
                                 __builtin_expect(err2==-1,0) || 
                                 __builtin_expect(err3==-1,0) || 
                                 __builtin_expect(err4==-1,0))
                              {
 #if (FAST_TERMINATE_ON_CRITICAL_ERROR) == 1
                                __builtin_trap();
#else
                                std::terminate();
#endif                                   
                              }
                           }
                           else {
                               gms_mm_free(this->m_Exr);this->m_Exr = NULL;
                               gms_mm_free(this->m_Exi);this->m_Exi = NULL;
                               gms_mm_free(this->m_Eyr);this->m_Eyr = NULL;
                               gms_mm_free(this->m_Eyi);this->m_Eyi = NULL;
                               
                           }
                      }
                      
                      DC2D_r4_t & operator=(const DC2D_r4_t &) = delete;
                      
                      inline DC2D_r4_t & operator=(DC2D_r4_t &&rhs) noexcept(true)
                       {
                            using namespace gms::common;
                            if(this==&rhs) return (*this);
                            
                            gms_swap(this->mnx,rhs.mnx);
                            gms_swap(this->mny,rhs.mny);
                            gms_swap(this->m_Exr,rhs.m_Exr);
                            gms_swap(this->m_Exi,rhs.m_Exi);
                            gms_swap(this->m_Eyr,rhs.m_Eyr);
                            gms_swap(this->m_Eyi,rhs.m_Eyi);
                            rhs.mnx       = 0ULL;
                            rhs.mny       = 0ULL;
                            rhs.m_Exr     = NULL;
                            rhs.m_Exi     = NULL;
                            rhs.m_Eyr     = NULL;
                            rhs.m_Eyi     = NULL;
                            return (*this);
                      }

                   inline std::size_t bytes_mnx() const noexcept(true)
                    {
                         return (sizeof(float)*this->mnx);
                    }

                   inline std::size_t bytes_mny() const noexcept(true)
                    {
                         return (sizeof(float)*this->mny);
                    }

                   private: 
                   inline void allocate() noexcept(false)
                   {
                        using namespace gms::common;
                        this->m_Exr  = (float*)gms_mm_malloc(sizeof(float)*this->mnx,64ULL);
                        this->m_Exi  = (float*)gms_mm_malloc(sizeof(float)*this->mnx,64ULL);
                        this->m_Eyr  = (float*)gms_mm_malloc(sizeof(float)*this->mny,64ULL);
                        this->m_Eyi  = (float*)gms_mm_malloc(sizeof(float)*this->mny,64ULL);
                       
                   }  
                   public: 
                   inline void info_size_alignment() const 
                   {
                        std::cout  << "alignof(struct DC2D_r4_t) = "  << alignof(DC2D_r4_t) << '\n'
                                  << "sizeof(struct  DC2D_r4_t) = "  << sizeof(DC2D_r4_t)  << '\n'
                                  << std::hex << std::showbase       << '\n'
                                  << "&this->m_Exr              =" << (void*)&this->m_Exr << "\n"
                                     "&this->m_Exi              =" << (void*)&this->m_Exi << "\n"
                                     "&this->m_Eyr              =" << (void*)&this->m_Eyr << "\n"
                                     "&this->m_Eyi              =" << (void*)&this->m_Eyi << "\n"
                                     "&this->mnx                =" << &this->mnx   <<  "\n"
                                     "&this->mny                =" << &this->mny   <<  "\n"
                                     "&this->ismmap             =" << &this->ismmap << "\n";
                   }  
                   
              };
              
              
              
               struct alignas(64) DC1D_r4_t final
               {
                     
                      float * __restrict m_Exr;
                      float * __restrict m_Exi;
                      std::size_t        mnx;
                      bool              ismmap;

                      inline DC1D_r4_t() noexcept(true)
                      {
                      
                         this->mnx   = 0ULL;
                         this->m_Exr  = NULL;
                         this->m_Exi  = NULL;
                        
                      }                    
                      
                      inline explicit DC1D_r4_t(const std::size_t nx) noexcept(false)
                       {
                                       
                             this->mnx = nx;
                             allocate();
                             this->ismmap = false;
                      }  
                      
                      inline  DC1D_r4_t(const std::size_t nx,
                                        const int32_t prot,
                                        const int32_t flags,
                                        const int32_t fd,
                                        const long offset,
                                        const int32_t fsize) noexcept(false)
                     {
                             using namespace gms::common;
                             this->mnx = nx;
                             switch (fsize) {
                                 case 0:
                                      this->m_Exr = (float*)
                                                  gms_mmap_4KiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->m_Exi = (float*)
                                                  gms_mmap_4KiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 1:
                                      this->m_Exr = (float*)
                                                  gms_mmap_2MiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->m_Exi = (float*)
                                                  gms_mmap_2MiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 2:
                                      this->m_Exr = (float*)
                                                  gms_mmap_1GiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->m_Exi = (float*)
                                                  gms_mmap_1GiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     } 
                       
                      //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline DC1D_r4_t(const std::vector<float> &Exr,
                                       const std::vector<float> &Exi) noexcept(false)
                      {
                                    
                               this->mnx = Exr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(float)*this->mnx;
                               std::memcpy(this->m_Exr,&Exr[0],lenx);
                               std::memcpy(this->m_Exi,&Exi[0],lenx);
                                                              
                      }
                      
                      
                       //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline DC1D_r4_t(const std::valarray<float> &Exr,
                                       const std::valarray<float> &Exi) noexcept(false)
                      {
                               
                               this->mnx = Exr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(float)*this->mnx;
                               std::memcpy(this->m_Exr,&Exr[0],lenx);
                               std::memcpy(this->m_Exi,&Exi[0],lenx);
                                                               
                      }
                      
                      inline  DC1D_r4_t(const std::size_t nx,
                                        const float * __restrict Exr,   
                                        const float * __restrict Exi) noexcept(false)
                      {
                          using namespace gms::common;                                 
                          this->mnx = nx;
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_Exr,&Exr[0],lenx);
                          std::memcpy(this->m_Exi,&Exi[0],lenx);
                      }
                      
                      inline DC1D_r4_t(DC1D_r4_t &&rhs) noexcept(true)
                     {
                           
                          this->mnx    = rhs.mnx;
                          this->m_Exr  = &rhs.m_Exr[0];
                          this->m_Exi  = &rhs.m_Exi[0];
                          rhs.mnx      = 0ULL;
                          rhs.m_Exr    = NULL;
                          rhs.m_Exi    = NULL;
                                                    
                      }   
                        
                      DC1D_r4_t(const DC1D_r4_t &)             = delete;
                      
                      inline ~DC1D_r4_t() noexcept(false)
                      {
                           using namespace gms::common;
                           if(this->ismmap) {
                               int32_t err1{}, err2{};
                              err1 = gms_ummap<float>(this->m_Exr,this->mnx); this->m_Exr = NULL;
                              err2 = gms_ummap<float>(this->m_Exi,this->mnx); this->m_Exi = NULL;
                              if(__builtin_expect(err1==-1,0) || 
                                 __builtin_expect(err2==-1,0))
                              {
#if (FAST_TERMINATE_ON_CRITICAL_ERROR) == 1
                                __builtin_trap();
#else
                                std::terminate();
#endif                                       
                              }                             
                           }
                           else {
                               gms_mm_free(this->m_Exr);this->m_Exr = NULL;
                               gms_mm_free(this->m_Exi);this->m_Exi = NULL;
                                                           
                           }
                      }
                      
                      DC1D_r4_t & operator=(const DC1D_r4_t &) = delete;
                      
                      inline DC1D_r4_t & operator=(DC1D_r4_t &&rhs) noexcept(true)
                      {
                            using namespace gms::common;
                            if(this==&rhs) return (*this);
                            
                            gms_swap(this->mnx,rhs.mnx);
                            gms_swap(this->m_Exr,rhs.m_Exr);
                            gms_swap(this->m_Exi,rhs.m_Exi);
                            rhs.mnx         = 0ULL;
                            rhs.m_Exr       = NULL;
                            rhs.m_Exi       = NULL;
                            return (*this);
                      }

                   inline std::size_t bytes_mnx() const noexcept(true)
                    {
                         return (sizeof(float)*this->mnx);
                    }

                   private: 
                   inline void allocate() noexcept(false)
                   {
                        using namespace gms::common;
                        this->m_Exr  = (float*)gms_mm_malloc(sizeof(float)*this->mnx,64ULL);
                        this->m_Exi  = (float*)gms_mm_malloc(sizeof(float)*this->mnx,64ULL);
                                             
                   }    
                   public: 
                   inline void info_size_alignment() const 
                   {
                         std::cout  << "alignof(struct DC1D_r4_t) = "   << alignof(DC1D_r4_t) << '\n'
                                   << "sizeof(struct  DC1D_r4_t) = "   << sizeof(DC1D_r4_t)  << '\n'
                                   << std::hex << std::showbase        << '\n'
                                   << "&this->m_Exr              ="    << (void*)&this->m_Exr  << "\n"
                                      "&this->m_Exi              ="    << (void*)&this->m_Exi  << "\n"
                                      "&this->mnx                ="    << &this->mnx    <<  "\n"
                                      "&this->ismmap             ="    << &this->ismmap << "\n";
                   }
                   
              };
              
              
            template<typename T, typename Enable = void>
            struct DC3D_t {};  
              
            template<typename T>
            struct alignas(64) DC3D_t<T,
            typename std::enable_if<std::is_floating_point<T>::value>::type> {
                     
                      T * __restrict m_Ex;
                      T * __restrict m_Ey;
                      T * __restrict m_Ez;
                      std::size_t        m_nx; 
                      std::size_t        m_ny; 
                      std::size_t        m_nz;
                      bool               ismmap;

                     inline DC3D_t() noexcept(true) 
                     {
                          
                          this->m_nx      = 0ULL;
                          this->m_ny      = 0ULL;
                          this->m_nz      = 0ULL
                          this->m_Ex      = NULL;
                          this->m_Ey      = NULL;
                          this->m_Ez      = NULL;
                      } 
                          
                     inline  DC3D_t(const std::size_t nx,
                                    const std::size_t ny,
                                    const std::size_t nz) noexcept(false)
                    {
                                 
                          this->m_nx = _nx;
                          this->m_ny = _ny;
                          this->m_nz = _nz;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline  DC3D_t(  const std::size_t nx,
                                      const std::size_t ny,
                                      const std::size_t nz,
                                      const int32_t prot,
                                      const int32_t flags,
                                      const int32_t fd,
                                      const long offset,
                                      const int32_t fsize) noexcept(false)
                     {
                             using namespace gms::common;
                             this->m_nx = nx;
                             this->m_ny = ny;
                             this->m_nz = nz;
                             switch (fsize) {
                                 case 0:
                                      this->m_Ex = (T*)
                                                 gms_mmap_4KiB<T>(this->m_nx,prot,flags,fd,offset);
                                      this->m_Ey = (T*)
                                                 gms_mmap_4KiB<T>(this->m_ny,prot,flags,fd,offset);
                                      this->m_Ez = (T*)
                                                 gms_mmap_4KiB<T>(this->m_nz,prot,flags,fd,offset);           
                                      this->ismmap = true;
                                 break;
                                 case 1:
                                      this->m_Ex = (T*)
                                                 gms_mmap_2MiB<T>(this->m_nx,prot,flags,fd,offset);
                                      this->m_Ey = (T*)
                                                 gms_mmap_2MiB<T>(this->m_ny,prot,flags,fd,offset);
                                      this->m_Ez = (T*)
                                                 gms_mmap_2MiB<T>(this->m_nz,prot,flags,fd,offset);    
                                      this->ismmap = true;
                                 break;
                                 case 2:
                                      this->m_Ex = (T*)
                                                 gms_mmap_1GiB<T>(this->m_nx,prot,flags,fd,offset);
                                      this->m_Ey = (T*)
                                                 gms_mmap_1GiB<T>(this->m_ny,prot,flags,fd,offset);
                                      this->m_Ez = (T*)
                                                 gms_mmap_1GiB<T>(this->m_nz,prot,flags,fd,offset);  
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                      
                    inline   DC3D_t(   const std::vector<T> &Ex,
                                       const std::vector<T> &Ey,
                                       const std::vector<T> &Ez) noexcept(false)
                    {   
                                      
                          this->m_nx = Ex.size();
                          this->m_ny = Ey.size();
                          this->m_nz = Ez.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(T)*this->m_nx;
                          const std::size_t leny = sizeof(T)*this->m_ny;
                          const std::size_t lenz = sizeof(T)*this->m_nz;
                          std::memcpy(this->m_Ex,&Ex[0],lenx);
                          std::memcpy(this->m_Ey,&Ey[0],leny);  
                          std::memcpy(this->m_Ez,&Ez[0],lenz);                           
                     }
                     
                    inline   DC3D_t(   const std::valarray<T> &Ex,
                                       const std::valarray<T> &Ey,
                                       const std::valarray<T> &Ez) noexcept(false)
                    {   
                                      
                          this->m_nx = Ex.size();
                          this->m_ny = Ey.size();
                          this->m_nz = Ez.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(T)*this->m_nx;
                          const std::size_t leny = sizeof(T)*this->m_ny;
                          const std::size_t lenz = sizeof(T)*this->m_nz;
                          std::memcpy(this->m_Ex,&Ex[0],lenx);
                          std::memcpy(this->m_Ey,&Ey[0],leny);  
                          std::memcpy(this->m_Ez,&Ez[0],lenz);       
                                      
                  }    
                                           
                                                  
                   inline    DC3D_t(  const std::size_t nx,
                                      const std::size_t ny,
                                      const std::size_t nz,
                                      const T * __restrict Ex,
                                      const T * __restrict Ey,
                                      const T * __restrict Ez ) {
                                    
                          this->m_nx  = nx;
                          this->m_ny  = ny;
                          this->m_nz  = nz;
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(T)*this->m_nx;
                          const std::size_t leny = sizeof(T)*this->m_ny;
                          const std::size_t lenz = sizeof(T)*this->m_nz;
                          std::memcpy(this->m_Ex,&Ex[0],lenx);
                          std::memcpy(this->m_Ey,&Ey[0],leny);  
                          std::memcpy(this->m_Ez,&Ez[0],lenz);
                   }  
                   
                  inline  DC3D_t(DC3D_t && rhs) noexcept(true)
                  {
                          
                          this->m_nx     = rhs.m_nx;
                          this->m_ny     = rhs.m_ny;
                          this->m_nz     = rhs.m_nz
                          this->m_Ex     = &rhs.m_Ex[0];
                          this->m_Ey     = &rhs.m_Ey[0];
                          this->m_Ez     = &rhs.m_Ez[0];
                          rhs.m_nx       = 0ULL;
                          rhs.m_ny       = 0ULL;
                          rhs.m_nz       = 0ULL;
                          rhs.m_Ex       = NULL;
                          rhs.m_Ey       = NULL;
                          rhs.m_Ez       = NULL;     
                 }
                                 
                   DC3D_t(const DC3D_t &)     = delete;
                      
                   inline   ~DC3D_t() noexcept(false)
                   {
                      
                          using namespace gms::common;
                          if(this->ismmap) 
                          {
                             int32_t err1{}, err2{}, err3{};
                             err1 = gms_ummap<T>(this->m_Ex,this->m_nx); this->m_Ex = NULL;
                             err2 = gms_ummap<T>(this->m_Ey,this->m_ny); this->m_Ey = NULL;
                             err3 = gms_ummap<T>(this->m_Ez,this->m_nz); this->m_Ez = NULL; 
                             if(__builtin_expect(err1==-1,0) || 
                                __builtin_expect(err2==-1,0) ||
                                __builtin_expect(err3==-1,0))
                              {
#if (FAST_TERMINATE_ON_CRITICAL_ERROR) == 1
                                __builtin_trap();
#else
                                std::terminate();
#endif
                              }
                          }
                          else 
                          {
                              gms_mm_free(this->m_nx); this->m_Ex = NULL;
                              gms_mm_free(this->m_ny); this->m_Ey = NULL;
                              gms_mm_free(this->m_nz); this->m_Ez = NULL;
                          }
                      }
                      
                    DC3D_t & operator=(const DC3D_t &) = delete;
                      
                    inline  DC3D_t & operator=(DC3D_t &&rhs) noexcept(true)
                    {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           
                           gms_swap(this->m_nx,rhs.m_nx);
                           gms_swap(this->m_ny,rhs.m_ny);
                           gms_swap(this->m_nz,rhs.m_nz);
                           gms_swap(this->m_Ex,rhs.m_Ex);
                           gms_swap(this->m_Ey,rhs.m_Ey);
                           gms_swap(this->m_Ez,rhs.m_Ez);
                           rhs.m_nx        = 0ULL;
                           rhs.m_ny        = 0ULL;
                           rhs.m_nz        = 0ULL;
                           rhs.m_Ex        = NULL;
                           rhs.m_Ey        = NULL;
                           rhs.m_Ez        = NULL;
                           return (*this);
                      }
                      
                    inline void allocate() noexcept(true)
                    {
                        using namespace gms::common;
                        this->m_Ex =   (T*)
                                         gms_mm_malloc(sizeof(T)*this->m_nx,64ULL);
                        this->m_Ey =   (T*)
                                         gms_mm_malloc(sizeof(T)*this->m_ny,64ULL);
                        this->m_Ez =   (T*)
                                         gms_mm_malloc(sizeof(T)*this->m_nz,64ULL);
                       
                    }

                    
                    
           };
           
           
           template<typename T, typename Enable = void>
           struct DC2D_t {};
           
            template<typename T>
            struct alignas(64) DC2D_t<T,
            typename std::enable_if<std::is_floating_point<T>::value>::type> 
            {
                                  
                      T * __restrict m_Ex;
                      T * __restrict m_Ey;
                      std::size_t    nx; 
                      std::size_t    ny; 
                      bool           ismmap;

                     inline DC2D_t() noexcept(true) 
                     {
                          
                          this->nx    = 0ULL;
                          this->ny    = 0ULL;
                          this->m_Ex    = NULL;
                          this->m_Ey    = NULL;
                    } 
                          
                     inline  DC2D_t(const std::size_t _nx,
                                    const std::size_t _ny) noexcept(false)
                    {                                  
                                 
                          this->nx = _nx;
                          this->ny = _ny;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline  DC2D_t(  const std::size_t _nx,
                                      const std::size_t _ny,
                                      const int32_t prot,
                                      const int32_t flags,
                                      const int32_t fd,
                                      const long offset,
                                      const int32_t fsize) noexcept(false)
                    {
                             using namespace gms::common;
                             this->nx = _nx;
                             this->ny = _ny;
                             switch (fsize) {
                                 case 0:
                                      this->m_Ex = (T*)
                                                 gms_mmap_4KiB<T>(this->nx,prot,flags,fd,offset);
                                      this->m_Ey = (T*)
                                                 gms_mmap_4KiB<T>(this->ny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 1:
                                      this->m_Ex = (T*)
                                                 gms_mmap_2MiB<T>(this->nx,prot,flags,fd,offset);
                                      this->m_Ey = (T*)
                                                 gms_mmap_2MiB<T>(this->ny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 2:
                                      this->m_Ex = (T*)
                                                 gms_mmap_1GiB<T>(this->nx,prot,flags,fd,offset);
                                      this->m_Ey = (T*)
                                                 gms_mmap_1GiB<T>(this->ny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                      
                    inline   DC2D_t(   const std::vector<T> &Ex,
                                       const std::vector<T> &Ey) noexcept(true)
                    {
                                       
                          this->nx = Ex.size();
                          this->ny = Ey.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(T)*this->nx;
                          const std::size_t leny = sizeof(T)*this->ny;
                          std::memcpy(this->m_Ex,&Ex[0],lenx);
                          std::memcpy(this->m_Ey,&Ey[0],leny);  
                                           
                     }
                     
                    inline   DC2D_t(   const std::valarray<T> &Ex,
                                       const std::valarray<T> &Ey) noexcept(true)
                    {
                                      
                          this->nx = Ex.size();
                          this->ny = Ey.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(T)*this->nx;
                          const std::size_t leny = sizeof(T)*this->ny;
                          std::memcpy(this->m_Ex,&Ex[0],lenx);
                          std::memcpy(this->m_Ey,&Ey[0],leny);  
                                                   
                  }    
                                           
                                                  
                   inline   DC2D_t(   const std::size_t _nx,
                                      const std::size_t _ny,
                                      const T * __restrict Ex,
                                      const T * __restrict Ey) noexcept(true)
                    {
                                    
                          this->nx  = _nx;
                          this->ny  = _ny;
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(T)*this->nx;
                          const std::size_t leny = sizeof(T)*this->ny;
                          const std::size_t lenx = sizeof(T)*this->nx;
                          const std::size_t leny = sizeof(T)*this->ny;
                          std::memcpy(this->m_Ex,&Ex[0],lenx);
                          std::memcpy(this->m_Ey,&Ey[0],leny); 
                   }  
                   
                  inline  DC2D_t(DC2D_t && rhs) noexcept(false)
                  {
                          
                          this->nx       = rhs.nx;
                          this->ny       = rhs.ny;
                          this->m_Ex     = &rhs.m_Ex[0];
                          this->m_Ey     = &rhs.m_Ey[0];
                          rhs.nx         = 0ULL;
                          rhs.ny         = 0ULL;
                          rhs.m_Ex       = NULL;
                          rhs.m_Ey       = NULL;
                 }
                                 
                   DC2D_t(const DC2D_t &)     = delete;
                      
                   inline   ~DC2D_t() noexcept(false) 
                   {
                          using namespace gms::common;
                          if(this->ismmap)
                          {
                             int32_t err1{}, err2{};
                             err1 = gms_ummap<T>(this->m_Ex,this->nx); this->m_Ex = NULL;
                             err2 = gms_ummap<T>(this->m_Ey,this->ny); this->m_Ey = NULL;
                             if(__builtin_expect(err1==-1,0) || 
                                __builtin_expect(err2==-1,0))
                              {
#if (FAST_TERMINATE_ON_CRITICAL_ERROR) == 1
                                __builtin_trap();
#else
                                std::terminate();
#endif                                   
                              }
                          }
                          else 
                          {
                              gms_mm_free(this->m_Ex); this->m_Ex = NULL;
                              gms_mm_free(this->m_Ey); this->m_Ey = NULL; 
                           
                          }
                      }
                      
                    DC2D_t & operator=(const DC2D_t &) = delete;
                      
                    inline  DC2D_t & operator=(DC2D_t &&rhs) noexcept(true)
                    {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           
                           gms_swap(this->nx,rhs.nx);
                           gms_swap(this->ny,rhs.ny);
                           gms_swap(this->m_Ex,rhs.m_Ex);
                           gms_swap(this->m_Ey,rhs.m_Ey);
                           rhs.nx        = 0ULL;
                           rhs.ny        = 0ULL;
                           rhs.m_Ex      = NULL;
                           rhs.m_Ey      = NULL;
                           return (*this);
                      }
                      
                    inline void allocate() noexcept(true)
                    {
                        using namespace gms::common;
                        this->m_Ex =   (T*)
                                         gms_mm_malloc(sizeof(T)*this->nx,64ULL);
                        this->m_Ey =   (T*)
                                         gms_mm_malloc(sizeof(T)*this->ny,64ULL);
                                              
                    }
           };
           
            template<typename T,
            typename  Enable = void>
            struct  DC1D_t {};
           
           template<typename T>
            struct alignas(64) DC1D_t<T,
            typename std::enable_if<std::is_floating_point<T>::value>::type> {
                     
                      T * __restrict     m_Ex;
                      std::size_t        nx; 
                      bool               ismmap;

                     inline DC1D_t() noexcept(true) 
                     {
                          
                          this->nx      = 0ULL;
                          this->m_Ex    = NULL;
                     } 
                          
                     inline explicit DC1D_t(const std::size_t _nx) noexcept(false)
                     {
                                                             
                          this->nx = _nx;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline DC1D_t(   const std::size_t _nx,
                                      const int32_t prot,
                                      const int32_t flags,
                                      const int32_t fd,
                                      const long offset,
                                      const int32_t fsize) noexcept(false)
                     {
                             using namespace gms::common;
                             this->nx = _nx;
                             switch (fsize) {
                                 case 0:
                                      this->m_Ex = (T*)
                                                 gms_mmap_4KiB<T>(this->nx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 1:
                                      this->m_Ex = (T*)
                                                 gms_mmap_2MiB<T>(this->nx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 2:
                                      this->m_Ex = (T*)
                                                 gms_mmap_1GiB<T>(this->nx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                      
                    inline   DC1D_t(   const std::vector<T> &Ex) noexcept(false)
                    {
                          this->nx = Ex.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(T)*this->nx;
                          std::memcpy(this->m_Ex,&Ex[0],lenx);
                                                             
                     }
                     
                    inline   DC1D_t(   const std::valarray<T> &Ex) noexcept(false)
                    {
                          this->nx = Ex.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(T)*this->nx;
                          std::memcpy(this->m_Ex,&Ex[0],lenx);
                                                                    
                  }    
                                           
                                                  
                   inline   DC1D_t(   const std::size_t _nx,
                                      const T * __restrict Ex) noexcept(false)
                  {
                                      
                          this->nx  = _nx;
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(T)*this->nx;
                          std::memcpy(this->m_Ex,&Ex[0],lenx);
                   }  
                   
                  inline  DC1D_t(DC1D_t && rhs) noexcept(true)
                  {
                          
                          this->nx       = rhs.nx;
                          this->m_Ex     = &rhs.m_Ex[0];
                          rhs.nx         = 0ULL;
                          rhs.m_Ex       = NULL;
                  }
                                 
                   DC1D_t(const DC1D_t &)     = delete;
                      
                   inline   ~DC1D_t() noexcept(false)
                   {
                          using namespace gms::common;
                          if(this->ismmap) 
                          {  
                             int32_t err{};
                             err = gms_ummap<T>(this->m_Ex,this->nx); this->m_Ex = NULL;
                             if(__builtin_expect(err==-1,0))
                             {
#if (FAST_TERMINATE_ON_CRITICAL_ERROR) == 1
                                __builtin_trap();
#else
                                std::terminate();
#endif                                  
                             }
                          }
                          else 
                          {
                              gms_mm_free(this->m_Ex); this->m_Ex = NULL;
                          }                            
                                                                                        
                   }
                      
                    DC1D_t & operator=(const DC1D_t &) = delete;
                      
                    inline  DC1D_t & operator=(DC1D_t &&rhs) noexcept(true)
                    {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                                                     
                           gms_swap(this->nx,rhs.nx);
                           gms_swap(this->m_Ex,rhs.m_Ex);
                           rhs.nx        = 0ULL;
                           rhs.m_Ex      = NULL;
                           return (*this);
                      }
                      
                    inline void allocate() noexcept(false)
                    {
                        using namespace gms::common;
                        this->m_Ex =   (T*)
                                         gms_mm_malloc(sizeof(T)*this->nx,64ULL);
                                                                    
                    }
           };
           
           
         
           
              
              

              
              




} //gms

















#endif /*__GMS_DYN_CONTAINERS_HPP__*/
