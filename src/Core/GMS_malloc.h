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

#ifndef __GMS_MALLOC_H__
#define __GMS_MALLOC_H__ 011020191914



namespace  file_info {

    const unsigned int GMS_MALLOC_MAJOR = 2U;
    const unsigned int GMS_MALLOC_MINOR = 1U;
    const unsigned int GMS_MALLOC_MICRO = 0U;
    const unsigned int GMS_MALLOC_FULLVER =
      1000U*GMS_MALLOC_MAJOR+100U*GMS_MALLOC_MINOR+10U*GMS_MALLOC_MICRO;
    const char * const GMS_MALLOC_CREATE_DATE = "01-10-2019 19:14 +00200 (TUE 01 OCT 2019 GMT+2)";
    const char * const GMS_MALLOC_BUILD_DATE  = __DATE__ ":" __TIME__;
    const char * const GMS_MALLOC_AUTHOR      =  "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";
    const char * const GMS_MALLOC_DESCRIPT    =  "Malloc wrappers.";
}

#include <cstdint>
#include <malloc.h>
#include <sys/mman.h>
#include <linux/mman.h>
#include <cassert>
#include <omp.h> // OMP allocators
#include <alloca.h>
#include <iostream>
#include <cerrno>
#include <cstring>
#include <cxxabi.h> // for demangle function 
#include "GMS_config.h"
#include "tbb/scalable_allocator.h"



//Handy macros for page alignment (mmap)
#if !defined(NORMAL_PAGE_4KiB)
    #define NORMAL_PAGE_4KiB (4* 1024)
#endif

#if !defined(ALIGN_TO_PAGE_4KiB)
     #define ALIGN_TO_PAGE_4KiB(x) \
        (((x) + NORMAL_PAGE_4KiB - 1) / NORMAL_PAGE_4KiB * NORMAL_PAGE_4KiB)
#endif

#if !defined(HUGE_PAGE_2MiB)
    #define HUGE_PAGE_2MiB (2 * 1024 * 1024)
#endif

#if !defined(ALIGN_TO_PAGE_2MiB)
     #define ALIGN_TO_PAGE_2MiB(x) \
        (((x) + HUGE_PAGE_2MiB - 1) / HUGE_PAGE_2MiB * HUGE_PAGE_2MiB)
#endif

#if !defined(HUGE_PAGE_1GiB)
    #define HUGE_PAGE_1GiB (1024 * 1024 * 1024)
#endif

#if !defined(ALIGN_TO_PAGE_1GiB)
    #define ALIGN_TO_PAGE_1GiB(x) \
       (((x) + HUGE_PAGE_1GiB - 1) / HUGE_PAGE_1GiB * HUGE_PAGE_1GiB)
#endif




namespace gms {
	namespace common {


	  /*
               Removing Windows support!!
               Bumping version:  major to 2 (complete rewrite).
               Adding mmap allocations aligned to various page size
               boundaries.
               Modified: 11-10-2020 11:14AM +00200
               Modified: 30-08-2021 09:44AM +00200
            */

#if defined(__INTEL_COMPILER) || defined(__ICC)	
	    __forceinline 
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
       __attribute__((always_inline))
#endif       
	   	std::size_t
	    add_cacheline_pad_float(const std::size_t len) {
        std::size_t misalign = 0ULL;
		std::size_t pad = 0ULL;
		misalign = len % 16;
		pad = len + (misalign == 0 ? 0 : 16 - misalign);
		return (pad);
	    }


#if defined(__INTEL_COMPILER) || defined(__ICC)	
	    __forceinline 
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
       __attribute__((always_inline))
#endif 	    
	    std::size_t
	    add_cacheline_pad_double(const std::size_t len) {
        std::size_t misalign = 0ULL;
		std::size_t pad = 0ULL;
		misalign = len % 8;
		pad = len + (misalign == 0 ? 0 : 8 - misalign);
		return (pad);
   	    }

		
		
#if defined(__INTEL_COMPILER) || defined(__ICC)	
	    __forceinline 
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
       __attribute__((always_inline))
#endif 
	    void * gms_mm_malloc(const std::size_t len,
	                           const std::size_t alignment) {

		     void * __restrict ptr = NULL;
	         ptr = _mm_malloc(len, alignment);
             if (NULL == ptr && len != 0ULL) 
			 {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	                 std::cerr << " Display of the call-stack not implemented yet!!" << "\n";
#endif
	                 std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << "\n";
					 std::cerr <<  "[errno]: --"  << std::strerror(errno) << "\n";
	                 std::cerr << "at: " __FILE__ ":"    << __LINE__ << " in: " <<  __PRETTY_FUNCTION__ << "\n";
	                 std::exit(EXIT_FAILURE);
	         }
                     return (ptr);
	      }



#if defined(__INTEL_COMPILER) || defined(__ICC)	
	    __forceinline 
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
       __attribute__((always_inline))
#endif 	       
	      void gms_mm_free(void * __restrict ptr) {
                   _mm_free(ptr);
	      }

		  inline char * 
		  demangle(const char * mangled_name,
		           int32_t & status)
		  {
			   std::size_t buf_size{512};
			   char * out_buf{reinterpret_cast<char*>(gms_mm_malloc(buf_size,64ULL))};
			   char  * result{NULL};
               int32_t stat{1};
               result = abi::__cxa_demangle(mangled_name,out_buf,&buf_size,&stat);
			   status = stat;
			   return result;
		  }

#if !defined(__GNUC__)
		  /*TBB-based allocators*/
#if defined(__INTEL_COMPILER) || defined(__ICC)	
	    __forceinline 
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
       __attribute__((always_inline))
#endif           
	      void * __restrict gms_tbb_malloc(const std::size_t nbytes,
		                        const std::size_t alignment)
		  {
			     assert(nbytes>0);
				 void * __restrict ptr = NULL;
                 ptr                   = scalable_aligned_malloc(nbytes,alignment);
				 if (NULL == ptr && nbytes != 0ULL) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	                  std::cerr << " Not implemented yet!!" << "\n";
#endif
	                   std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << std::strerror(errno) << "\n";
	                   std::cerr << "at: " __FILE__ ":" <<  __LINE__ << " in: " << __PRETTY_FUNCTION__ << "\n";
	                   std::exit(EXIT_FAILURE);
	            }
				return ptr;
		  }

#if defined(__INTEL_COMPILER) || defined(__ICC)	
	    __forceinline 
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
       __attribute__((always_inline))
#endif 		   
		   void gms_tbb_free(void * __restrict ptr)
		   {
			    if(ptr == NULL) return;
				scalable_aligned_free(ptr);
		   }

#endif 

#define ALIGNED_ALLOCA(p, s, a) { \
  char *c_ptr = alloca((s) + (a)); \
  int64 off1 = (a) - 1; \
  int64 i_ptr64 = ((int64)c_ptr) + off1; \
  (p) = (void*)(i_ptr64 - (i_ptr64 & off1)); \
}

             


	        //
		// Using MMAP allocations for 4KiB page, 2MiB page and 1GiB page.
		// Subroutines include an error handling
		// See: https://man7.org/linux/man-pages/man2/mmap.2.html
		//
		// void *mmap(void *addr, size_t length, int prot, int flags, std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << "MEMORY ALLOCATION FAILURE!!" << "\n";
	                 
                //  int fd, off_t offset);



		 
		 template<typename T>
		 void * gms_mmap_4KiB(const std::size_t length,
		                      const int32_t prot,
				      const int32_t flags,
				      const int32_t fd,
				      const off_t offset) {

                      void * ptr = nullptr; 
                      std::size_t totmem = sizeof(T)*length;   
                      std::size_t nlargep= totmem/(4096ULL);
                      if(totmem != nlargep*4096ULL) nlargep++;
                      totmem = nlargep*4096ULL;           
                      ptr = mmap(NULL,totmem,prot,flags,fd,offset); 
                      if((ptr == (void*)(-1))) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	                        std::cerr << "Requested stack-backtrace -- not implemented yet!!" << "\n";
#endif
                            std::cerr << "value of ptr="  << std::hex << ptr << "\n";
		                    std::cerr << "[" << __DATE__ << " : " << __TIME__ << " ] " << "[ERROR]:" << std::strerror(errno) << "\n";
	                        std::cerr << "at: " __FILE__ ":" << __LINE__ << " in: " << __PRETTY_FUNCTION__ << "\n";
	                        std::exit(EXIT_FAILURE);
                        }
			          return (ptr);
		 }
				          

       
		template<typename T>
		void *  gms_mmap_2MiB(const std::size_t length,
		                      const int32_t prot,
				      const int32_t flags,
				      const int32_t fd,
				      const off_t offset) {
				     
				      
                      void * ptr = nullptr; 
                      std::size_t totmem = sizeof(T)*length;   
                      std::size_t nlargep= totmem/(2097152ULL);
                      if(totmem != nlargep*2097152ULL) nlargep++;
                      totmem = nlargep*2097152ULL;           
                      ptr = mmap(NULL,totmem,prot,flags,fd,offset); 
		              if((ptr == (void*)(-1))) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	                           std::cerr << "Requested stack-backtrace -- not implemented yet!!" << "\n";
#endif
                               std::cerr << "value of ptr= " <<  std::hex << ptr << "\n";
		                       std::cerr << "[ " << __DATE__ << " : " << __TIME__ << "] " << std::strerror(errno) << "\n";
	                           std::cerr << "at: " __FILE__ ":" << __LINE__ << " in: " << __PRETTY_FUNCTION__ << "\n";
	                           std::exit(EXIT_FAILURE);
                       }
                       return (ptr);
		}


		
		template<typename T>
                void * gms_mmap_1GiB(const std::size_t length,
		                     const int32_t prot,
				     const int32_t flags,
				     const int32_t fd,
				     const off_t offset) {

		              void * ptr = nullptr; 
                      std::size_t totmem = sizeof(T)*length;   
                      std::size_t nlargep= totmem/(1073741824ULL);
                      if(totmem != nlargep*1073741824ULL) nlargep++;
                      totmem = nlargep*1073741824ULL;           
                      ptr = mmap(NULL,totmem,prot,flags,fd,offset); 
		              if((ptr == (void*)(-1))) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	                         std::cerr << "Requested stack-backtrace -- not implemented yet!!" << "\n";
#endif
                             std::cerr << "value of ptr= " <<  std::hex << ptr << "\n";
		                     std::cerr << "[ " << __DATE__ << " : " << __TIME__ << "] " << std::strerror(errno) << "\n";
	                         std::cerr << "at: " __FILE__ ":" << __LINE__ << " in: " << __PRETTY_FUNCTION__ << "\n";
	                         std::exit(EXIT_FAILURE);
                       }
		              return (ptr);
		}

		
		
		template<typename T>		       
                int gms_ummap(void * __restrict ptr,
		              const std::size_t len) {
                     return (munmap(ptr,sizeof(T)*len));
		}



	    
		
		template<typename T> 
		void gms_swap(T &a, T&b) 
		{
			 T tmp = std::move(a);
			 a     = std::move(b);
			 b     = std::move(tmp);
		}

		/* TBB Memory Pool class.*/


             
              
#if __cplusplus >= 202002L	

#ifndef TBB_PREVIEW_MEMORY_POOL
#define TBB_PREVIEW_MEMORY_POOL 1
#endif 
#include <type_traits>
#include <concepts>
#include <memory>
#include "oneapi/tbb/scalable_allocator.h"
#include "oneapi/tbb/memory_pool.h"

			     /*Slightly modified from the tit/core/par/memory_pool.h*/
			     template<typename T>
				 requires std::is_object_v<T> && std::is_trivially_destructible_v<T>
				 struct TbbMemoryPool final 
				 {
                       
                       public:

                       template<typename... Args>
					   requires std::constructible_from<T, Args&&...>
                       auto create_mem_pool(Args&&... args)->T* __restrict noexcept(false)
					   {
						   assert(m_pool != nullptr);
						   auto * __restrict const ptr = 
						             static_cast<T* __restrict>(m_pool->malloc(sizeof(T)));
						   if (nullptr == ptr) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	                           std::cerr << " Not implemented yet!!";
#endif
	                           ABORT_ON_ERROR(__PRETTY_FUNCTION__, MALLOC_FAILED)
	                       }

						   return std::construct_at(ptr,std::forward<Args>(args)...);
					   }
					   private:

					   std::unique_ptr<tbb::memory_pool<tbb::scalable_allocator<T>>> m_pool = 
					   std::make_unique<typename decltype(m_pool)::element_type();
				 };
#endif
				 /*example of usage*/
				 /*
				       TEST_CASE("par::MemoryPool") {
                       struct Struct {
                         int data_1;
                         int data_2;
                       };
                      par::MemoryPool<Struct> pool{};
                      auto* const root = pool.create(10, 20);
                     CHECK(root->data_1 == 10);
                     CHECK(root->data_2 == 20);
}

				 */

				/*Helper enum classes for tbb and mm_malloc functions calls*/
                enum class MEM_ALLOC_FUNC_TYPE : int32_t
				{
                           CALL_MM_MALLOC,
						   CALL_TBB_MALLOC
				};

				enum class MEM_FREE_FUNC_TYPE : int32_t 
				{
					       CALL_MM_FREE,
						   CALL_TBB_FREE
				};


	} // common
} // gms


#endif /*__GMS_MALLOC_H__*/
