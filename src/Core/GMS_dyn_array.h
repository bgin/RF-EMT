


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
!*/

#ifndef __GMS_DYN_ARRAY_H__
#define __GMS_DYN_ARRAY_H__ 180820251017



namespace file_info 
{

     static const unsigned int GMS_DYN_ARRAY_MAJOR = 1;
     static const unsigned int GMS_DYN_ARRAY_MINOR = 1;
     static const unsigned int GMS_DYN_ARRAY_MICRO = 0;
     static const unsigned int GMS_DYN_ARRAY_FULLVER =
       1000U*GMS_DYN_ARRAY_MAJOR+100U*GMS_DYN_ARRAY_MINOR+
       10U*GMS_DYN_ARRAY_MICRO;
     static const char GMS_DYN_ARRAY_CREATION_DATE[] = "18-08-2025 10:17 +00200 (MON  18 AUG 2025 GMT+2)";
     static const char GMS_DYN_ARRAY_BUILD_DATE[]    = __DATE__; 
     static const char GMS_DYN_ARRAY_BUILD_TIME[]  = __TIME__;
     static const char GMS_DYN_ARRAY_SYNOPSIS[]    = "Dynamically allocated, 64-byte aligned array.";

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
#include "GMS_half.h"


// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (USE_GMS_DYN_ARRAY_NT_STORES)
#define USE_GMS_DYN_ARRAY_NT_STORES 0
#endif

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (DYN_ARRAY_USE_PMC_INSTRUMENTATION)
#define DYN_ARRAY_USE_PMC_INSTRUMENTATION 0
#endif 

#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
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

namespace gms 
{

               // Simple dynamically allocated arrays without vector-like functionality (growing,resizing, ...etc)
               // Used mainly for the baseband/narrowband signal representation and henceforth processing.
               //using namespace half_float;

               struct alignas(32) darray_c2_t final 
               {
                      
                      std::complex<half_float::half> * __restrict m_data;
                      std::size_t                      mnx;
                     
                     inline darray_c2_t() noexcept(true)
                     {
                          
                          this->mnx     = 0ULL;
                          this->m_data  = NULL;
                                  
                     } 
                          
                     inline explicit darray_c2_t(const std::size_t nx) noexcept(false)
                     {
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif
                          this->mnx = nx;
                          allocate();
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif
                      }  
                      
                                         
                     inline darray_c2_t(const std::vector<std::complex<half_float::half>> &data) noexcept(false)
                     {    //shall be of the same size (no error checking implemented)
                                  
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                          
                          this->mnx = data.size();
                          allocate();
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                                                  
                     }
                     
                     inline darray_c2_t(const std::valarray<std::complex<half_float::half>> &data) noexcept(false)
                    {
                                   
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                 
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                           
                          this->mnx = data.size();
                          allocate();
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
                   HW_PMC_COLLECTION_EPILOGE_BODY

                   HW_PMC_SHOW_INFO
#endif                                                    
                     }
                             
                             
                      
                     inline  darray_c2_t(const std::size_t nx,
                                         const std::complex<half_float::half> * __restrict data) noexcept(false)
                     {
                                 
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                           
                          this->mnx = nx;
                          allocate();
                          const std::size_t lenx = sizeof(std::complex<half_float::half>)*this->mnx;
                          std::memcpy(this->m_data,&data[0],lenx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
                HW_PMC_COLLECTION_EPILOGE_BODY

                HW_PMC_SHOW_INFO
#endif
                     }  
                   
                     inline  darray_c2_t(darray_c2_t && rhs) noexcept(true)
                     {
                          
                          this->mnx    =  rhs.mnx;
                          this->m_data =  &rhs.m_data[0];
                          rhs.mnx      =  0ULL;
                          rhs.m_data   =  NULL;
                                             
                      }
                                 
                     inline darray_c2_t(const darray_c2_t &rhs)  noexcept(false)
                     {
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                  
                           this->mnx = rhs.mnx;
                           this->allocate();
                           std::memcpy(this->m_data,&rhs.m_data[0],this->mnx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                                
                     }
                      
                     inline ~darray_c2_t() noexcept(false)
                     {
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                          
                          if(__builtin_expect(this->m_data!=NULL,1))
                          {     
                             gms_mm_free(this->m_data); this->m_data = NULL;  
                          } 
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                          
                     }    
                         
                                                                 
                     darray_c2_t & operator=(const darray_c2_t &rhs) noexcept(false) 
                     {
                           using namespace gms::common;
                           if(__builtin_expect(this==&rhs,0)) { return (*this);}

                           gms_mm_free(this->m_data);
                           this->mnx = rhs.mnx;
                           this->allocate();
                           std::memcpy(this->m_data,&rhs.m_data[0],this->mnx);
                           return (*this);
                     } 
                      
                     inline darray_c2_t & operator=(darray_c2_t &&rhs) noexcept(true)
                     {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                                                                                                                
                           gms_swap(this->mnx,rhs.mnx);
                           gms_swap(this->m_data,rhs.m_data);
                           rhs.mnx       = 0ULL;
                           rhs.m_data    = NULL;
                           return (*this);
                      }

                    inline std::size_t bytes_mnx() const noexcept(true)
                    {
                         return (sizeof(std::complex<half_float::half>)*this->mnx);
                    }

                    inline std::complex<half_float::half> * begin() const noexcept
                    {
                         return this->m_data;
                    }

                    inline std::complex<half_float::half> * end() const noexcept
                    {
                         return this->m_data+this->mnx;
                    }

                    inline constexpr std::size_t size() const noexcept
                    {
                           return this->mnx;
                    }

                    inline bool empty() const noexcept
                    {
                         return (this->size() == 0ull);
                    }

                    private:
                    inline void allocate() noexcept(false)
                    {
                          using namespace gms::common;
                          this->m_data  = (std::complex<half_float::half>*)
                                         gms_mm_malloc( sizeof(std::complex<half_float::half>)*this->mnx,64ULL);
                       
                     }
                     public: 
                     inline void info_size_alignment() const 
                     {
                          std::cout << "alignof(struct darray_c2_t) = " << alignof(darray_c2_t) << '\n'
                                    << "sizeof(struct  darray_c2_t) = " << sizeof(darray_c2_t)  << '\n'
                                    << std::hex << std::showbase        << '\n'
                                    << "&this->m_data               ="  << (void*)&this->m_data  << "\n"
                                       "&this->mnx                  ="  << &this->mnx            << "\n";
                     
                     }
                      
               };

               struct alignas(32) darray_c4_t final 
               {
                      
                      std::complex<float> * __restrict m_data;
                      std::size_t                      mnx;
                      
                     inline darray_c4_t() noexcept(true)
                     {
                          
                          this->mnx     = 0ULL;
                          this->m_data  = NULL;
                                   
                     } 
                          
                     inline explicit darray_c4_t(const std::size_t nx) noexcept(false)
                     {
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif
                          this->mnx = nx;
                          allocate();
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif
                      }  
                      
                                          
                     inline darray_c4_t(const std::vector<std::complex<float>> &data) noexcept(false)
                     {    //shall be of the same size (no error checking implemented)
                                  
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                          
                          this->mnx = data.size();
                          allocate();
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                                                  
                     }
                     
                     inline darray_c4_t(const std::valarray<std::complex<float>> &data) noexcept(false)
                    {
                                   
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                 
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                           
                          this->mnx = data.size();
                          allocate();
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
                   HW_PMC_COLLECTION_EPILOGE_BODY

                   HW_PMC_SHOW_INFO
#endif                                                    
                     }
                             
                             
                      
                     inline  darray_c4_t(const std::size_t nx,
                                         const std::complex<float> * __restrict data) noexcept(false)
                     {
                                 
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                           
                          this->mnx = nx;
                          allocate();
                          const std::size_t lenx = sizeof(std::complex<float>)*this->mnx;
                          std::memcpy(this->m_data,&data[0],lenx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
                HW_PMC_COLLECTION_EPILOGE_BODY

                HW_PMC_SHOW_INFO
#endif
                     }  
                   
                     inline  darray_c4_t(darray_c4_t && rhs) noexcept(true)
                     {
                          
                          this->mnx    =  rhs.mnx;
                          this->m_data =  &rhs.m_data[0];
                          rhs.mnx      =  0ULL;
                          rhs.m_data   =  NULL;
                                             
                      }
                                 
                     darray_c4_t(const darray_c4_t &rhs)  noexcept(false)
                     {
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                  
                           this->mnx = rhs.mnx;
                           this->allocate();
                           std::memcpy(this->m_data,&rhs.m_data[0],this->mnx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                                
                     }
                      
                     inline ~darray_c4_t() noexcept(false)
                     {
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                          
                          if(__builtin_expect(this->m_data!=NULL,1))
                          {
                             gms_mm_free(this->m_data); this->m_data = NULL;  
                          }  
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                          
                     }    
                         
                                                                 
                     darray_c4_t & operator=(const darray_c4_t &rhs) noexcept(false)
                     {
                          using namespace gms::common;
                          if(__builtin_expect(this==&rhs,0)) {return (*this);} 
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif  
                          gms_mm_free(this->m_data);
                          this->mnx = rhs.mnx;
                          this->allocate();
                          std::memcpy(this->m_data,&rhs.m_data[0],this->mnx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                            
                          return (*this); 
                     }
                      
                     inline darray_c4_t & operator=(darray_c4_t &&rhs) noexcept(true)
                     {
                           using namespace gms::common;
                           if(__builtin_expect(this==&rhs,0)) {return (*this);}
                                                                                                                
                           gms_swap(this->mnx,rhs.mnx);
                           gms_swap(this->m_data,rhs.m_data);
                           rhs.mnx       = 0ULL;
                           rhs.m_data    = NULL;
                           return (*this);
                      }

                    inline std::size_t bytes_mnx() const noexcept(true)
                    {
                         return (sizeof(std::complex<float>)*this->mnx);
                    }

                    inline std::complex<float> * begin() const noexcept
                    {
                         return this->m_data;
                    }

                    inline std::complex<float> * end() const noexcept
                    {
                         return this->m_data+this->mnx;
                    }

                    inline constexpr std::size_t size() const noexcept
                    {
                           return this->mnx;
                    }

                    inline bool empty() const noexcept
                    {
                         return (this->size() == 0ull);
                    }

                    private:
                    inline void allocate() noexcept(false)
                    {
                          using namespace gms::common;
                          this->m_data  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->mnx,64ULL);
                       
                     }
                     public: 
                     inline void info_size_alignment() const 
                     {
                          std::cout << "alignof(struct darray_c4_t) = " << alignof(darray_c4_t) << '\n'
                                    << "sizeof(struct  darray_c4_t) = " << sizeof(darray_c4_t)  << '\n'
                                    << std::hex << std::showbase        << '\n'
                                    << "&this->m_data               ="  << (void*)&this->m_data  << "\n"
                                       "&this->mnx                  ="  << &this->mnx            << "\n";
                                      
                     }
                      
               };


               struct alignas(32) darray_c8_t final 
               {
                      
                      std::complex<double> * __restrict m_data;
                      std::size_t                      mnx;
                     
                     inline darray_c8_t() noexcept(true)
                     {
                          
                          this->mnx     = 0ULL;
                          this->m_data  = NULL;
                                
                     } 
                          
                     inline explicit darray_c8_t(const std::size_t nx) noexcept(false)
                     {
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif          
                          this->mnx = nx;
                          allocate();
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                          
                      }  
                      
                                          
                     inline darray_c8_t(const std::vector<std::complex<double>> &data) noexcept(false)
                     {    //shall be of the same size (no error checking implemented)
                                  
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                 
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif  
                          this->mnx = data.size();
                          allocate();
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                                                
                     }
                     
                     inline darray_c8_t(const std::valarray<std::complex<double>> &data) noexcept(false)
                    {
                                   
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif  
                          this->mnx = data.size();
                          allocate();
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                                                    
                     }
                             
                             
                      
                     inline  darray_c8_t(const std::size_t nx,
                                         const std::complex<double> * __restrict data) noexcept(false)
                     {
                                 
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                   
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                          
                          this->mnx = nx;
                          allocate();
                          const std::size_t lenx = sizeof(std::complex<double>)*this->mnx;
                          std::memcpy(this->m_data,&data[0],lenx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif
                     }  
                   
                     inline  darray_c8_t(darray_c8_t && rhs) noexcept(true)
                     {
                          
                          this->mnx    =  rhs.mnx;
                          this->m_data =  &rhs.m_data[0];
                          rhs.mnx      =  0ULL;
                          rhs.m_data   =  NULL;
                                             
                      }
                                 
                     darray_c8_t(const darray_c8_t &rhs)  noexcept(false)
                     {
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                  
                           this->mnx = rhs.mnx;
                           this->allocate();
                           std::memcpy(this->m_data,&rhs.m_data[0],this->mnx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                                
                     }
                      
                     inline ~darray_c8_t() noexcept(false)
                     {
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif
                          if(__builtin_expect(this->m_data!=NULL,1))
                          {
                             gms_mm_free(this->m_data); this->m_data = NULL;  
                          }  
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif
                     }    
                         
                                                                 
                     darray_c8_t & operator=(const darray_c8_t &rhs) noexcept(false) 
                     {
                           using namespace gms::common;
                           if(__builtin_expect(this==&rhs,0)) {return (*this);}
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif
                            gms_mm_free(this->m_data);
                            this->mnx = rhs.mnx;
                            this->allocate();
                            std::memcpy(this->m_data,&rhs.m_data[0],this->mnx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif
                             return (*this);                            
                     }
                      
                     inline darray_c8_t & operator=(darray_c8_t &&rhs) noexcept(true)
                     {
                           using namespace gms::common;
                           if(__builtin_expect(this==&rhs,0)) {return (*this);}
                                                                                                          
                           gms_swap(this->mnx,rhs.mnx);
                           gms_swap(this->m_data,rhs.m_data);
                           rhs.mnx       = 0ULL;
                           rhs.m_data    = NULL;
                           return (*this);
                      }

                    inline std::size_t bytes_mnx() const noexcept(true)
                    {
                         return (sizeof(std::complex<double>)*this->mnx);
                    }

                    inline std::complex<double> * begin() const noexcept
                    {
                         return this->m_data;
                    }

                    inline std::complex<double> * end() const noexcept
                    {
                         return this->m_data+this->mnx;
                    }

                    inline constexpr std::size_t size() const noexcept
                    {
                           return this->mnx;
                    }

                    inline bool empty() const noexcept
                    {
                         return (this->size() == 0ull);
                    }

                    private:
                    inline void allocate() noexcept(false)
                    {
                          using namespace gms::common;
                          this->m_data  = (std::complex<double>*)
                                         gms_mm_malloc( sizeof(std::complex<double>)*this->mnx,64ULL);
                       
                    }
                    public: 
                    inline void info_size_alignment() const 
                    {
                          std::cout << "alignof(struct darray_c8_t) = " << alignof(darray_c8_t) << '\n'
                                    << "sizeof(struct  darray_c8_t) = " << sizeof(darray_c8_t)  << '\n'
                                    << std::hex << std::showbase        << '\n'
                                    << "&this->m_data               ="  << (void*)&this->m_data  << "\n"
                                       "&this->mnx                  ="  << &this->mnx            << "\n";
                                      
                    }
                      
               };


                struct alignas(32) darray_r2_t final 
                {
                      
                      half_float::half * __restrict    m_data;
                      std::size_t                      mnx;
                     

                     inline darray_r2_t() noexcept(true)
                     {
                          
                          this->mnx     = 0ULL;
                          this->m_data  = NULL;
                                
                     } 
                          
                     inline explicit darray_r2_t(const std::size_t nx) noexcept(false)
                     {
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                          
                          this->mnx = nx;
                          allocate();
                          
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif
                      }  
                      
                    
                      
                     inline darray_r2_t(const std::vector<half_float::half> &data) noexcept(false)
                     {    //shall be of the same size (no error checking implemented)
                                  
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                 
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                          
                          this->mnx = data.size();
                          allocate();
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                                                
                     }
                     
                    inline darray_r2_t(const std::valarray<half_float::half> &data) noexcept(false)
                    {
                                   
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                 
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                         
                          this->mnx = data.size();
                          allocate();
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                                                  
                     }
                             
                             
                      
                     inline  darray_r2_t(const std::size_t nx,
                                         const half_float::half * __restrict data) noexcept(false)
                     {
                                 
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                   
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                         
                          this->mnx = nx;
                          allocate();
                          const std::size_t lenx = sizeof(half_float::half)*this->mnx;
                          std::memcpy(this->m_data,&data[0],lenx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                         
                     }  
                   
                     inline  darray_r2_t(darray_r2_t && rhs) noexcept(true)
                     {
                          
                          this->mnx    =  rhs.mnx;
                          this->m_data =  &rhs.m_data[0];
                          rhs.mnx      =  0ULL;
                          rhs.m_data   =  NULL;
                                               
                      }
                                 
                     darray_r2_t(const darray_r2_t &rhs)   noexcept(false)
                     {
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                  
                           this->mnx = rhs.mnx;
                           this->allocate();
                           std::memcpy(this->m_data,&rhs.m_data[0],this->mnx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                                
                     }
                      
                     inline ~darray_r2_t() noexcept(false)
                     {
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                   
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif

                          if(__builtin_expect(this->m_data!=NULL,1))
                          {
                             gms_mm_free(this->m_data); this->m_data = NULL;  
                          }  
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif
                     }    
                         
                                                                 
                     darray_r2_t & operator=(const darray_r2_t &rhs) noexcept(false) 
                     {
                           using namespace gms::common;
                           if(__builtin_expect(this==&rhs,0))  {return (*this);}
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                   
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif
                           gms_mm_free(this->m_data);
                           this->mnx = rhs.mnx;
                           this->allocate();
                           std::memcpy(this->m_data,&rhs.m_data[0],this->mnx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif
                            return (*this);
                     }
                      
                     inline darray_r2_t & operator=(darray_r2_t &&rhs) noexcept(true)
                     {
                           using namespace gms::common;
                           if(__builtin_expect(this==&rhs,0))  {return (*this);}
                                                                           
                           gms_swap(this->mnx,rhs.mnx);
                           gms_swap(this->m_data,rhs.m_data);
                           rhs.mnx       = 0ULL;
                           rhs.m_data    = NULL;
                           return (*this);
                      }

                    inline std::size_t bytes_mnx() const noexcept(true)
                    {
                         return (sizeof(float)*this->mnx);
                    }

                    inline half_float::half * begin() const noexcept
                    {
                         return this->m_data;
                    }

                    inline half_float::half * end() const noexcept
                    {
                         return this->m_data+this->mnx;
                    }

                    inline constexpr std::size_t size() const noexcept
                    {
                           return this->mnx;
                    }

                    inline bool empty() const noexcept
                    {
                         return (this->size() == 0ull);
                    }

                    private:
                    inline void allocate() noexcept(false)
                    {
                          using namespace gms::common;
                          this->m_data  = (half_float::half*)
                                         gms_mm_malloc( sizeof(half_float::half)*this->mnx,64ULL);
                       
                     }
                     public: 
                     inline void info_size_alignment() const 
                     {
                          std::cout << "alignof(struct darray_r2_t) = " << alignof(darray_r2_t) << '\n'
                                    << "sizeof(struct  darray_r2_t) = " << sizeof(darray_r2_t)  << '\n'
                                    << std::hex << std::showbase        << '\n'
                                    << "&this->m_data               ="  << (void*)&this->m_data  << "\n"
                                       "&this->mnx                  ="  << &this->mnx            << "\n";
                                      
                     }
                      
               };

               struct alignas(32) darray_r4_t final 
               {
                      
                      float * __restrict               m_data;
                      std::size_t                      mnx;
                                            
                     inline darray_r4_t() noexcept(true)
                     {
                          
                          this->mnx     = 0ULL;
                          this->m_data  = NULL;
                                
                     } 
                          
                     inline explicit darray_r4_t(const std::size_t nx) noexcept(false)
                     {
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                          
                          this->mnx = nx;
                          allocate();
                         
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif
                      }  
                      
                    
                      
                     inline darray_r4_t(const std::vector<float> &data) noexcept(false)
                     {    //shall be of the same size (no error checking implemented)
                                  
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                 
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                          
                          this->mnx = data.size();
                          allocate();
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                                                
                     }
                     
                    inline darray_r4_t(const std::valarray<float> &data) noexcept(false)
                    {
                                   
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                 
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                         
                          this->mnx = data.size();
                          allocate();
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                                                  
                     }
                             
                             
                      
                     inline  darray_r4_t(const std::size_t nx,
                                         const float * __restrict data) noexcept(false)
                     {
                                 
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                   
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                         
                          this->mnx = nx;
                          allocate();
                          const std::size_t lenx = sizeof(float)*this->mnx;
                          std::memcpy(this->m_data,&data[0],lenx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                         
                     }  
                   
                     inline  darray_r4_t(darray_r4_t && rhs) noexcept(true)
                     {
                          
                          this->mnx     =  rhs.mnx;
                          this->m_data  =  &rhs.m_data[0];
                          rhs.mnx       =  0ULL;
                          rhs.m_data    =  NULL;
                                           
                      }
                                 
                     darray_r4_t(const darray_r4_t &rhs)   noexcept(false)
                     {
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                  
                           this->mnx = rhs.mnx;
                           this->allocate();
                           std::memcpy(this->m_data,&rhs.m_data[0],this->mnx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                                
                     }
                      
                     inline ~darray_r4_t() noexcept(false)
                     {
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                   
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif

                          if(__builtin_expect(this->m_data!=NULL,1))
                          {
                             gms_mm_free(this->m_data); this->m_data = NULL;  
                          }  
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif
                     }    
                         
                     
                                                                                    
                     darray_r4_t & operator=(const darray_r4_t &rhs) noexcept(false) 
                     {
                           using namespace gms::common;
                           if(__builtin_expect(this==&rhs,0)) { return (*this);}
 #if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                   
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                          
                           gms_mm_free(this->m_data);
                           this->mnx = rhs.mnx;
                           this->allocate();
                           std::memcpy(this->m_data,&rhs.m_data[0],this->mnx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                           
                           return (*this);
                     }


                     inline darray_r4_t & operator=(darray_r4_t &&rhs) noexcept(true)
                     {
                           using namespace gms::common;
                           if(__builtin_expect(this==&rhs,0)) {return (*this);}
                                                                           
                           gms_swap(this->mnx,rhs.mnx);
                           gms_swap(this->m_data,rhs.m_data);
                           rhs.mnx       = 0ULL;
                           rhs.m_data    = NULL;
                         
                           return (*this);
                      }

                    inline std::size_t bytes_mnx() const noexcept(true)
                    {
                         return (sizeof(float)*this->mnx);
                    }

                    inline float * begin() const noexcept
                    {
                         return this->m_data;
                    }

                    inline float * end() const noexcept
                    {
                         return this->m_data+this->mnx;
                    }

                    inline constexpr std::size_t size() const noexcept
                    {
                           return this->mnx;
                    }

                    inline bool empty() const noexcept
                    {
                         return (this->size() == 0ull);
                    }

                    private:
                    inline void allocate() noexcept(false)
                    {
                          using namespace gms::common;
                          this->m_data  = (float*)
                                         gms_mm_malloc( sizeof(float)*this->mnx,64ULL);
                       
                     }
                     public: 
                     inline void info_size_alignment() const 
                     {
                          std::cout << "alignof(struct darray_r4_t) = " << alignof(darray_r4_t) << '\n'
                                    << "sizeof(struct  darray_r4_t) = " << sizeof(darray_r4_t)  << '\n'
                                    << std::hex << std::showbase        << '\n'
                                    << "&this->m_data               ="  << (void*)&this->m_data  << "\n"
                                       "&this->mnx                  ="  << &this->mnx            << "\n";
                                     
                     }
                      
               };


               struct alignas(32) darray_r8_t final 
               {
                      
                      double * __restrict              m_data;
                      std::size_t                      mnx;
                      

                     inline darray_r8_t() noexcept(true)
                     {
                          
                          this->mnx     = 0ULL;
                          this->m_data  = NULL;
                                   
                     } 
                          
                     inline explicit darray_r8_t(const std::size_t nx) noexcept(false)
                     {
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                         
                          this->mnx = nx;
                          allocate();
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif
                      }  
                      
                   
                      
                     inline darray_r8_t(const std::vector<double> &data) noexcept(false)
                     {    //shall be of the same size (no error checking implemented)
                                  
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                         
                          this->mnx = data.size();
                          allocate();
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                                                
                     }
                     
                    inline darray_r8_t(const std::valarray<double> &data) noexcept(false)
                    {
                                   
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif
                          this->mnx = data.size();
                          allocate();
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                                                  
                     }
                             
                             
                      
                     inline  darray_r8_t(const std::size_t nx,
                                         const float * __restrict data) noexcept(false)
                     {
                                 
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                          
                          this->mnx = nx;
                          allocate();
                          const std::size_t lenx = sizeof(float)*this->mnx;
                          std::memcpy(this->m_data,&data[0],lenx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif
                     }  
                   
                     inline  darray_r8_t(darray_r8_t && rhs) noexcept(true)
                     {
                          
                          this->mnx    =  rhs.mnx;
                          this->m_data =  &rhs.m_data[0];
                          rhs.mnx      =  0ULL;
                          rhs.m_data   =  NULL;
                                              
                      }
                                 
                     darray_r8_t(const darray_r8_t &rhs) noexcept(false)
                     {
                           using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif                  
                           this->mnx = rhs.mnx;
                           this->allocate();
                           std::memcpy(this->m_data,&rhs.m_data[0],this->mnx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                                
                     }
                      
                     inline ~darray_r8_t() noexcept(false)
                     {
                          using namespace gms::common;
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif
                          
                          if(__builtin_expect(this->m_data!=NULL,1)) 
                          {
                             gms_mm_free(this->m_data); this->m_data = NULL;  
                          }  
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif                          
                     }    
                         
                                                                 
                     darray_r8_t & operator=(const darray_r8_t &rhs) noexcept(false) 
                     {
                           using namespace gms::common;
                           if(__builtin_expect(this==&rhs,0))  {return (*this);}
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               PMC_VARS                  
               HW_PMC_COLLECTION_PROLOGE_BODY
#endif
                           gms_mm_free(this->m_data);
                           this->mnx = rhs.mnx;
                           this->allocate();
                           std::memcpy(this->m_data,&rhs.m_data[0],this->mnx);
#if (DYN_ARRAY_USE_PMC_INSTRUMENTATION) == 1
               HW_PMC_COLLECTION_EPILOGE_BODY

               HW_PMC_SHOW_INFO
#endif
                           return (*this);                           
                     }
                      
                     inline darray_r8_t & operator=(darray_r8_t &&rhs) noexcept(true)
                     {
                           using namespace gms::common;
                           if(__builtin_expect(this==&rhs,0))  {return (*this);}
                                                                           
                           gms_swap(this->mnx,rhs.mnx);
                           gms_swap(this->m_data,rhs.m_data);
                        
                           rhs.mnx       = 0ULL;
                           rhs.m_data    = NULL;
                           
                           return (*this);
                      }

                    inline std::size_t bytes_mnx() const noexcept(true)
                    {
                         return (sizeof(double)*this->mnx);
                    }

                    inline double * begin() const noexcept
                    {
                         return this->m_data;
                    }

                    inline double * end() const noexcept
                    {
                         return this->m_data+this->mnx;
                    }

                    inline constexpr std::size_t size() const noexcept
                    {
                           return this->mnx;
                    }

                    inline bool empty() const noexcept
                    {
                         return (this->size() == 0ull);
                    }

                    private:
                    inline void allocate() noexcept(false)
                    {
                          using namespace gms::common;
                          this->m_data  = (double*)
                                         gms_mm_malloc( sizeof(double)*this->mnx,64ULL);
                       
                     }
                     public: 
                     inline void info_size_alignment() const 
                     {
                          std::cout << "alignof(struct darray_r8_t) = " << alignof(darray_r8_t) << '\n'
                                    << "sizeof(struct  darray_r8_t) = " << sizeof(darray_r8_t)  << '\n'
                                    << std::hex << std::showbase        << '\n'
                                    << "&this->m_data               ="  << (void*)&this->m_data  << "\n"
                                       "&this->mnx                  ="  << &this->mnx            << "\n";
                                      
                     }
                      
               };


} //gms 


#endif /*__GMS_DYN_ARRAY_H__*/