
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

#ifndef __GMS_DEVICE_VECTOR_CUH__
#define __GMS_DEVICE_VECTOR_CUH__

#include <cstdint>
#include "GMS_gpu_config.cuh"


__global__
void fill_r4_kernel(float * __restrict__ data,
                    const std::size_t len,
                    const float value) {
   
    const uint32_t tid    = blockDim.x*blockIdx.x+threadIdx.x;
    const std::size_t idx = (std::size_t)tid;
    if(idx < len) {
        data[idx] = value;
    }
}

__global__
void fill_c4_kernel(cuComplex * __restrict__ data,
                    const std::size_t len,
                    const cuComplex value) {
  
     const uint32_t    tid  = blockDim.x*blockIdx.x+threadIdx.x;
     const std::size_t idx  = (std::size_t)tid;
     if(idx < len) {
         data[idx] = value;
     }
}

__global__
void fill_i4_kernel(int32_t * __restrict__ data,
                    const std::size_t len,
                    const int32_t value) {

     const uint32_t    tid = blockDim.x*blockIdx.x+threadIdx.x;
     const std::size_t idx = (std::size_t)tid;
     if(idx < len) {
         data[idx] = value;
     }
}

 


enum class Device_Memory {
      global  = 0,
      managed = 1
};


struct CudaVectorR4 {

    public:
    std::size_t m_len;
    float * __restrict__ m_data;
    std::size_t m_memsize; // size in bytes.
    Device_Memory m_memtype;
   
    
    public:
    
    CudaVectorR4() {
       m_len     = 0ULL;
       m_data    = nullptr;
       m_memsize = 0ULL;
       m_memtype = Device_Memory::global;
      
    }

   
    
    CudaVectorR4(const std::size_t len,
                 const Device_Memory memtype = Device_Memory::global) {
            
           allocate(len,memtype);
    }

    CudaVectorR4(const std::size_t len,
                 const Device_Memory memtype = Device_Memory::global,
                 const float value) {
          
            allocate(len,value,memtype);
    }

    CudaVectorR4(CudaVectorR4 &&other) {
        m_len           = other.m_len;
        m_data          = other.m_data;
        m_memsize       = other.m_memsize;
        m_memtype       = other.m_memtype;
        other.m_len     = 0ULL;
        other.m_data    = nullptr;
        other.m_memsize = 0ULL;
    }

    ~CudaVectorR4() {
        if(m_data) {
          cudaFree(m_data);
          m_data = nullptr;
        }
     }


    void allocate(const std::size_t len,
                  const Device_Memory memtype = Device_Memory::global) {
        
         m_len     = len;
         m_memsize = m_len*sizeof(float);
         m_memtype = memtype;
         cudaError_t status;
         if(m_memtype == Device_Memory::global) {
            GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&m_data,m_memsize));
         }
         else {
            GMS_CUDA_DEBUG_CHECK(cudaMallocManaged((void**)&m_data,m_memsize));
         }
          return;
        Error:
              if(__builtin_expect(status==cudaErrorMemoryAllocation,0)){
	         fatal_gpu_error(__PRETTY_FUNCTION__, status);
	       }             
        } 

     void allocate(const std::size_t len,
                   const float value,
                   const Device_Memory memtype = Device_Memory::global) {

         m_len     = len;
         m_memsize = m_len*sizeof(float);
         m_memtype = memtype;
         cudaError_t status;
         if(m_memtype == Device_Memory::global) {
            GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&m_data,m_memsize));
         }
         else {
            GMS_CUDA_DEBUG_CHECK(cudaMallocManaged((void**)&m_data,m_memsize));
         }
         fill(m_data,m_len,value);
         return;
         Error:
               if(__builtin_expect(status==cudaErrorMemoryAllocation,0)){
	          fatal_gpu_error(__PRETTY_FUNCTION__, status);
	       }  
     }
    
     void fill(const float value) {
          if(m_memtype==Device_Memory::global) {
             const uint32_t blockSize = 32;
             const uint32_t gridSize  = (m_size+blockSize-1)/blockSize;
             fill_r4_kernel<<<gridSize,blockSize>>>(m_data,m_len,value);
          }
          else {
              for(std::size_t i = 0ULL; i < m_len; ++i) {
                  m_data[i] = value;
              }
          }
     }

     void copy_cpu_gpu(const float * __restrict__ h_ptr,
                       cudaError_t &err) {
            cudaError_t status;
            GMS_CUDA_DEBUG_CHECK(cudaMemcpy(m_data,h_ptr,m_memsize,cudaMemcpyHostToDevice));
            return;
            Error:
                   err = status;
                   return;
     }

     void copy_gpu_cpu(float * __restrict__ h_ptr,
                       cudaError_t &err) {
            cudaError_t status;
            GMS_CUDA_DEBUG_CHECK(cudaMemcpy(h_ptr,m_data,m_memsize,cudaMemcpyDeviceToHost));
            return;
            Error:
                   err = status;
                   return;
     }

     void copy_gpu_gpu(float * __restrict__ d_ptr,
                       cudaError_t &err) {
             cudaError_t status;
             GMS_CUDA_DEBUG_CHECK(cudaMemcpy(d_ptr,m_data,m_memsize,cudaMemcpyDeviceToDevice));
             return;
             Error:
                   err = status;
                   return;
     }

     CudaVectorR4 & 
     operator=(CudaVectorR4 &&other) {
          if(this == &other) return (*this);
          cudaFree(m_data);
          m_len           = other.m_len;
          m_data          = other.m_data;
          m_memsize       = other.m_memsize; 
          m_memtype       = other.m_memtype;
          other.m_len     = 0ULL;
          other.m_data    = nullptr;
          other.m_memsize = 0ULL;
          return (*this);
     }

     float & operator[](const std::size_t idx) {
              return (m_data[idx]);
     }
   
};


struct CudaVectorC4 {

    public:
    std::size_t m_len;
    cuComplex * __restrict__ m_data;
    std::size_t m_memsize; // size in bytes.
    Device_Memory m_memtype;
   
    
    public:
    
    CudaVectorC4() {
       m_len     = 0ULL;
       m_data    = nullptr;
       m_memsize = 0ULL;
       m_memtype = Device_Memory::global;
      
    }

   
    
    CudaVectorC4(const std::size_t len,
                 const Device_Memory memtype = Device_Memory::global) {
            
           allocate(len,memtype);
    }

    CudaVectorC4(const std::size_t len,
                 const Device_Memory memtype = Device_Memory::global,
                 const float value) {
          
            allocate(len,value,memtype);
    }

    CudaVectorC4(CudaVectorR4 &&other) {
        m_len           = other.m_len;
        m_data          = other.m_data;
        m_memsize       = other.m_memsize;
        m_memtype       = other.m_memtype;
        other.m_len     = 0ULL;
        other.m_data    = nullptr;
        other.m_memsize = 0ULL;
    }

    ~CudaVectorC4() {
        if(m_data) {
          cudaFree(m_data);
          m_data = nullptr;
        }
     }


    void allocate(const std::size_t len,
                  const Device_Memory memtype = Device_Memory::global) {
        
         m_len     = len;
         m_memsize = m_len*sizeof(cuComplex);
         m_memtype = memtype;
         cudaError_t status;
         if(m_memtype == Device_Memory::global) {
            GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&m_data,m_memsize));
         }
         else {
            GMS_CUDA_DEBUG_CHECK(cudaMallocManaged((void**)&m_data,m_memsize));
         }
          return;
        Error:
              if(__builtin_expect(status==cudaErrorMemoryAllocation,0)){
	         fatal_gpu_error(__PRETTY_FUNCTION__, status);
	       }             
        } 

     void allocate(const std::size_t len,
                   const cuComplex value,
                   const Device_Memory memtype = Device_Memory::global) {

         m_len     = len;
         m_memsize = m_len*sizeof(cuComplex);
         m_memtype = memtype;
         cudaError_t status;
         if(m_memtype == Device_Memory::global) {
            GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&m_data,m_memsize));
         }
         else {
            GMS_CUDA_DEBUG_CHECK(cudaMallocManaged((void**)&m_data,m_memsize));
         }
         fill(m_data,m_len,value);
         return;
         Error:
               if(__builtin_expect(status==cudaErrorMemoryAllocation,0)){
	          fatal_gpu_error(__PRETTY_FUNCTION__, status);
	       }  
     }
    
     void fill(const cuComplex value) {
          if(m_memtype==Device_Memory::global) {
             const uint32_t blockSize = 32;
             const uint32_t gridSize  = (m_size+blockSize-1)/blockSize;
             fill_c4_kernel<<<gridSize,blockSize>>>(m_data,m_len,value);
          }
          else {
              for(std::size_t i = 0ULL; i < m_len; ++i) {
                  m_data[i] = value;
              }
          }
     }

     void copy_cpu_gpu(const cuComplex * __restrict__ h_ptr,
                       cudaError_t &err) {
            cudaError_t status;
            GMS_CUDA_DEBUG_CHECK(cudaMemcpy(m_data,h_ptr,m_memsize,cudaMemcpyHostToDevice));
            return;
            Error:
                   err = status;
                   return;
     }

     void copy_gpu_cpu(cuComplex * __restrict__ h_ptr,
                       cudaError_t &err) {
            cudaError_t status;
            GMS_CUDA_DEBUG_CHECK(cudaMemcpy(h_ptr,m_data,m_memsize,cudaMemcpyDeviceToHost));
            return;
            Error:
                   err = status;
                   return;
     }

     void copy_gpu_gpu(cuComplex * __restrict__ d_ptr,
                       cudaError_t &err) {
             cudaError_t status;
             GMS_CUDA_DEBUG_CHECK(cudaMemcpy(d_ptr,m_data,m_memsize,cudaMemcpyDeviceToDevice));
             return;
             Error:
                   err = status;
                   return;
     }

     CudaVectorC4 & 
     operator=(CudaVectorC4 &&other) {
          if(this == &other) return (*this);
          cudaFree(m_data);
          m_len           = other.m_len;
          m_data          = other.m_data;
          m_memsize       = other.m_memsize; 
          m_memtype       = other.m_memtype;
          other.m_len     = 0ULL;
          other.m_data    = nullptr;
          other.m_memsize = 0ULL;
          return (*this);
     }

     cuComplex & operator[](const std::size_t idx) {
              return (m_data[idx]);
     }
   
};


struct CudaVectorI4 {

    public:
    std::size_t m_len;
    int32_t * __restrict__ m_data;
    std::size_t m_memsize; // size in bytes.
    Device_Memory m_memtype;
   
    
    public:
    
    CudaVectorI4() {
       m_len     = 0ULL;
       m_data    = nullptr;
       m_memsize = 0ULL;
       m_memtype = Device_Memory::global;
      
    }

   
    
    CudaVectorI4(const std::size_t len,
                 const Device_Memory memtype = Device_Memory::global) {
            
           allocate(len,memtype);
    }

    CudaVectorI4(const std::size_t len,
                 const Device_Memory memtype = Device_Memory::global,
                 const float value) {
          
            allocate(len,value,memtype);
    }

    CudaVectorI4(CudaVectorI4 &&other) {
        m_len           = other.m_len;
        m_data          = other.m_data;
        m_memsize       = other.m_memsize;
        m_memtype       = other.m_memtype;
        other.m_len     = 0ULL;
        other.m_data    = nullptr;
        other.m_memsize = 0ULL;
    }

    ~CudaVectorI4() {
        if(m_data) {
          cudaFree(m_data);
          m_data = nullptr;
        }
     }


    void allocate(const std::size_t len,
                  const Device_Memory memtype = Device_Memory::global) {
        
         m_len     = len;
         m_memsize = m_len*sizeof(int32_t);
         m_memtype = memtype;
         cudaError_t status;
         if(m_memtype == Device_Memory::global) {
            GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&m_data,m_memsize));
         }
         else {
            GMS_CUDA_DEBUG_CHECK(cudaMallocManaged((void**)&m_data,m_memsize));
         }
          return;
        Error:
              if(__builtin_expect(status==cudaErrorMemoryAllocation,0)){
	         fatal_gpu_error(__PRETTY_FUNCTION__, status);
	       }             
        } 

     void allocate(const std::size_t len,
                   const int32_t value,
                   const Device_Memory memtype = Device_Memory::global) {

         m_len     = len;
         m_memsize = m_len*sizeof(int32_t);
         m_memtype = memtype;
         cudaError_t status;
         if(m_memtype == Device_Memory::global) {
            GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&m_data,m_memsize));
         }
         else {
            GMS_CUDA_DEBUG_CHECK(cudaMallocManaged((void**)&m_data,m_memsize));
         }
         fill(m_data,m_len,value);
         return;
         Error:
               if(__builtin_expect(status==cudaErrorMemoryAllocation,0)){
	          fatal_gpu_error(__PRETTY_FUNCTION__, status);
	       }  
     }
    
     void fill(const int32_t value) {
          if(m_memtype==Device_Memory::global) {
             const uint32_t blockSize = 32;
             const uint32_t gridSize  = (m_size+blockSize-1)/blockSize;
             fill_i4_kernel<<<gridSize,blockSize>>>(m_data,m_len,value);
          }
          else {
              for(std::size_t i = 0ULL; i < m_len; ++i) {
                  m_data[i] = value;
              }
          }
     }

     void copy_cpu_gpu(const int32_t * __restrict__ h_ptr,
                       cudaError_t &err) {
            cudaError_t status;
            GMS_CUDA_DEBUG_CHECK(cudaMemcpy(m_data,h_ptr,m_memsize,cudaMemcpyHostToDevice));
            return;
            Error:
                   err = status;
                   return;
     }

     void copy_gpu_cpu(int32_t * __restrict__ h_ptr,
                       cudaError_t &err) {
            cudaError_t status;
            GMS_CUDA_DEBUG_CHECK(cudaMemcpy(h_ptr,m_data,m_memsize,cudaMemcpyDeviceToHost));
            return;
            Error:
                   err = status;
                   return;
     }

     void copy_gpu_gpu(int32_t * __restrict__ d_ptr,
                       cudaError_t &err) {
             cudaError_t status;
             GMS_CUDA_DEBUG_CHECK(cudaMemcpy(d_ptr,m_data,m_memsize,cudaMemcpyDeviceToDevice));
             return;
             Error:
                   err = status;
                   return;
     }

     CudaVectorI4 & 
     operator=(CudaVectorI4 &&other) {
          if(this == &other) return (*this);
          cudaFree(m_data);
          m_len           = other.m_len;
          m_data          = other.m_data;
          m_memsize       = other.m_memsize; 
          m_memtype       = other.m_memtype;
          other.m_len     = 0ULL;
          other.m_data    = nullptr;
          other.m_memsize = 0ULL;
          return (*this);
     }

     int32_t & operator[](const std::size_t idx) {
              return (m_data[idx]);
     }
   
};












#endif /*__GMS_DEVICE_VECTOR_CUH__*/
