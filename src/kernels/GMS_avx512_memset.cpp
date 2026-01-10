



#include <immintrin.h>
#include "GMS_avx512_memset.h"

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx512")
#endif
void 
gms::common::
avx512_memset_unroll8x_ps(float * __restrict__ dst,const float filler,std::size_t size)
{
   
    const __m512 vfiller{_mm512_set1_ps(filler)};
    float * __restrict__ p_dst{dst};
    const float cfiller{filler};

    while(((uintptr_t)p_dst & 63) && size)
    {
        *p_dst++ = cfiller;
         --size;
    }

    while(size >= 128ULL)
    {
         
        _mm512_storeu_ps(p_dst+0, vfiller);
        _mm512_storeu_ps(p_dst+16, vfiller);
        _mm512_storeu_ps(p_dst+32, vfiller);
        _mm512_storeu_ps(p_dst+48,vfiller);
        _mm512_storeu_ps(p_dst+64,vfiller);
        _mm512_storeu_ps(p_dst+80,vfiller);
        _mm512_storeu_ps(p_dst+96,vfiller);
        _mm512_storeu_ps(p_dst+112,vfiller);
         
        size  -= 128ULL;
        p_dst += 128ULL;
    }

    while(size >= 96ULL)
    {
        _mm512_storeu_ps(p_dst+0, vfiller);
        _mm512_storeu_ps(p_dst+16, vfiller);
        _mm512_storeu_ps(p_dst+32, vfiller);
        _mm512_storeu_ps(p_dst+48,vfiller);
        _mm512_storeu_ps(p_dst+64,vfiller);
        _mm512_storeu_ps(p_dst+80,vfiller);

        size  -= 96ULL;
        p_dst += 96ULL;
    }

    while(size >= 64ULL)
    {
        _mm512_storeu_ps(p_dst+0, vfiller);
        _mm512_storeu_ps(p_dst+16, vfiller);
        _mm512_storeu_ps(p_dst+32, vfiller);
        _mm512_storeu_ps(p_dst+48,vfiller);

        size  -= 64ULL;
        p_dst += 64ULL;
    }

    while(size >= 32ULL)
    {
        _mm512_storeu_ps(p_dst+0, vfiller);
        _mm512_storeu_ps(p_dst+16, vfiller);

        size  -= 32ULL;
        p_dst += 32ULL;
    }

    while(size >= 16ULL)
    {
        _mm512_storeu_ps(p_dst+0, vfiller);

        size  -= 16ULL;
        p_dst += 16ULL;
    }

    while(size)
    {
        *p_dst = cfiller;
        size--;
        p_dst++;
    }
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx512")
#endif
void
gms::common::
avx512_memset_unroll16x_ps(float * __restrict__ dst,const float filler,std::size_t size)
{
    const __m512 vfiller{_mm512_set1_ps(filler)};
    float * __restrict__ p_dst = dst;
    std::size_t __i;
    const float cfiller{filler};

   
    for(__i = 0ull; (__i+255ull) < size; __i += 256ull)
    {
          _mm512_storeu_ps(&p_dst[__i+0ull],  vfiller);
          _mm512_storeu_ps(&p_dst[__i+16ull],  vfiller);
          _mm512_storeu_ps(&p_dst[__i+32ull],  vfiller);
          _mm512_storeu_ps(&p_dst[__i+48ull], vfiller);
          _mm512_storeu_ps(&p_dst[__i+64ull], vfiller);
          _mm512_storeu_ps(&p_dst[__i+80ULL], vfiller);
          _mm512_storeu_ps(&p_dst[__i+96ULL], vfiller);
          _mm512_storeu_ps(&p_dst[__i+112ULL], vfiller);
          _mm512_storeu_ps(&p_dst[__i+128ULL], vfiller);
          _mm512_storeu_ps(&p_dst[__i+144ULL], vfiller);
          _mm512_storeu_ps(&p_dst[__i+160ULL], vfiller);
          _mm512_storeu_ps(&p_dst[__i+176ULL], vfiller);
          _mm512_storeu_ps(&p_dst[__i+192ULL], vfiller);
          _mm512_storeu_ps(&p_dst[__i+208ULL], vfiller);
          _mm512_storeu_ps(&p_dst[__i+224ULL], vfiller);
          _mm512_storeu_ps(&p_dst[__i+240ULL], vfiller);
    }

    for(; (__i+191ull) < size;  __i += 192ull)
    {
          _mm512_storeu_ps(&p_dst[__i+0ull],  vfiller);
          _mm512_storeu_ps(&p_dst[__i+16ull],  vfiller);
          _mm512_storeu_ps(&p_dst[__i+32ull],  vfiller);
          _mm512_storeu_ps(&p_dst[__i+48ull], vfiller);
          _mm512_storeu_ps(&p_dst[__i+64ull], vfiller);
          _mm512_storeu_ps(&p_dst[__i+80ULL], vfiller);
          _mm512_storeu_ps(&p_dst[__i+96ULL], vfiller);
          _mm512_storeu_ps(&p_dst[__i+112ULL], vfiller);
          _mm512_storeu_ps(&p_dst[__i+128ULL], vfiller);
          _mm512_storeu_ps(&p_dst[__i+144ULL], vfiller);
          _mm512_storeu_ps(&p_dst[__i+160ULL], vfiller);
          _mm512_storeu_ps(&p_dst[__i+176ULL], vfiller);
    }

    for(; (__i+127ull) < size; __i += 128ull)
    {
          _mm512_storeu_ps(&p_dst[__i+0ull],  vfiller);
          _mm512_storeu_ps(&p_dst[__i+16ull],  vfiller);
          _mm512_storeu_ps(&p_dst[__i+32ull],  vfiller);
          _mm512_storeu_ps(&p_dst[__i+48ull], vfiller);
          _mm512_storeu_ps(&p_dst[__i+64ull], vfiller);
          _mm512_storeu_ps(&p_dst[__i+80ULL], vfiller);
          _mm512_storeu_ps(&p_dst[__i+96ULL], vfiller);
          _mm512_storeu_ps(&p_dst[__i+112ULL], vfiller);
    }

    for(; (__i+95ull) < size; __i += 96ull)
    {
          _mm512_storeu_ps(&p_dst[__i+0ull],  vfiller);
          _mm512_storeu_ps(&p_dst[__i+16ull],  vfiller);
          _mm512_storeu_ps(&p_dst[__i+32ull],  vfiller);
          _mm512_storeu_ps(&p_dst[__i+48ull], vfiller);
          _mm512_storeu_ps(&p_dst[__i+64ull], vfiller);
          _mm512_storeu_ps(&p_dst[__i+80ULL], vfiller);
    }

    for(; (__i+63ull) < size; __i += 64ull)
    {
          _mm512_storeu_ps(&p_dst[__i+0ull],  vfiller);
          _mm512_storeu_ps(&p_dst[__i+16ull],  vfiller);
          _mm512_storeu_ps(&p_dst[__i+32ull],  vfiller);
          _mm512_storeu_ps(&p_dst[__i+48ull], vfiller);
    }

    for(; (__i+31ull) < size; __i += 32ull)
    {
           _mm512_storeu_ps(&p_dst[__i+0ull],  vfiller);
           _mm512_storeu_ps(&p_dst[__i+16ull],  vfiller);
    }

    for(; (__i+15ull) < size; __i += 16ull)
    {
           _mm512_storeu_ps(&p_dst[__i+0ull],  vfiller);
    }

    for(; (__i+0ull) < size; __i += 1ull)
    {
          p_dst[__i] = cfiller;
    }

    
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx512")
#endif
void 
gms::common::
avx512_memset_unroll8x_pd(double * __restrict__ dst,const double filler,std::size_t size)
{
    const __m512d vfiller{_mm512_set1_pd(filler)};
    double * __restrict__ p_dst{dst};
    const double cfiller{filler};

    while(((uintptr_t)p_dst & 63) && size)
    {
        *p_dst++ = cfiller;
         --size;
    }

    while(size >= 64ULL)
    {
         
        _mm512_store_pd(p_dst+0ULL, vfiller);
        _mm512_store_pd(p_dst+8ULL, vfiller);
        _mm512_store_pd(p_dst+16ULL, vfiller);
        _mm512_store_pd(p_dst+24ULL, vfiller);
        _mm512_store_pd(p_dst+32ULL, vfiller);
        _mm512_store_pd(p_dst+40ULL,vfiller);
        _mm512_store_pd(p_dst+48ULL,vfiller);
        _mm512_store_pd(p_dst+56ULL,vfiller);
         
        p_dst += 64ULL;
        size  -= 64ULL;
        
    }

    while(size >= 48ULL)
    {
        _mm512_store_pd(p_dst+0ULL, vfiller);
        _mm512_store_pd(p_dst+8ULL, vfiller);
        _mm512_store_pd(p_dst+16ULL, vfiller);
        _mm512_store_pd(p_dst+24ULL, vfiller);
        _mm512_store_pd(p_dst+32ULL, vfiller);
        _mm512_store_pd(p_dst+40ULL,vfiller);
        
        p_dst += 48ULL;
        size  -= 48ULL;
        
    }

    while(size >= 32ULL)
    {
        _mm512_store_pd(p_dst+0ULL, vfiller);
        _mm512_store_pd(p_dst+8ULL, vfiller);
        _mm512_store_pd(p_dst+16ULL, vfiller);
        _mm512_store_pd(p_dst+24ULL, vfiller);
        
        p_dst += 32ULL;
        size  -= 32ULL;
       
    }

    while(size >= 16ULL)
    {
        _mm512_store_pd(p_dst+0ULL, vfiller);
        _mm512_store_pd(p_dst+8ULL, vfiller);
        
        p_dst += 16ULL;
        size  -= 16ULL;
        
    }

    while(size >= 8ULL)
    {
        _mm512_store_pd(p_dst+0ULL, vfiller);
        
        p_dst += 8ULL;
        size  -= 8ULL;
        
    }

    while(size)
    {
        *p_dst = cfiller;
         p_dst++;
         size--;
        
    }
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx512")
#endif
void
gms::common::
avx512_memset_unroll16x_pd(double * __restrict__ dst,const double filler,std::size_t size)
{
     const __m512d vfiller{_mm512_set1_pd(filler)};
     double * __restrict__ p_dst{dst};
     double cfiller{filler};

     while(((uintptr_t)p_dst & 63) && size)
     {
          *p_dst++ = cfiller;
          --size;
     }

     while(size >= 128ULL)
     {
          _mm512_store_pd(p_dst+0ULL, vfiller);
          _mm512_store_pd(p_dst+8ULL, vfiller);
          _mm512_store_pd(p_dst+16ULL, vfiller);
          _mm512_store_pd(p_dst+24ULL, vfiller);
          _mm512_store_pd(p_dst+32ULL, vfiller);
          _mm512_store_pd(p_dst+40ULL,vfiller);
          _mm512_store_pd(p_dst+48ULL,vfiller);
          _mm512_store_pd(p_dst+56ULL,vfiller);
          _mm512_store_pd(p_dst+64ULL,vfiller);
          _mm512_store_pd(p_dst+72ULL,vfiller);
          _mm512_store_pd(p_dst+80ULL,vfiller);
          _mm512_store_pd(p_dst+88ULL,vfiller);
          _mm512_store_pd(p_dst+96ULL,vfiller);
          _mm512_store_pd(p_dst+104ULL,vfiller);
          _mm512_store_pd(p_dst+112ULL,vfiller);
          _mm512_store_pd(p_dst+120ULL,vfiller);

          size  -= 128ULL;
          p_dst += 128ULL;
      }

      while(size >= 96ULL)
      {
          
          _mm512_store_pd(p_dst+0ULL, vfiller);
          _mm512_store_pd(p_dst+8ULL, vfiller);
          _mm512_store_pd(p_dst+16ULL, vfiller);
          _mm512_store_pd(p_dst+24ULL, vfiller);
          _mm512_store_pd(p_dst+32ULL, vfiller);
          _mm512_store_pd(p_dst+40ULL,vfiller);
          _mm512_store_pd(p_dst+48ULL,vfiller);
          _mm512_store_pd(p_dst+56ULL,vfiller);
          _mm512_store_pd(p_dst+64ULL,vfiller);
          _mm512_store_pd(p_dst+72ULL,vfiller);
          _mm512_store_pd(p_dst+80ULL,vfiller);
          _mm512_store_pd(p_dst+88ULL,vfiller);
          size  -= 96ULL;
          p_dst += 96ULL;
      }

      while(size >= 64ULL)
      {
          
          _mm512_store_pd(p_dst+0ULL, vfiller);
          _mm512_store_pd(p_dst+8ULL, vfiller);
          _mm512_store_pd(p_dst+16ULL, vfiller);
          _mm512_store_pd(p_dst+24ULL, vfiller);
          _mm512_store_pd(p_dst+32ULL, vfiller);
          _mm512_store_pd(p_dst+40ULL,vfiller);
          _mm512_store_pd(p_dst+48ULL,vfiller);
          _mm512_store_pd(p_dst+56ULL,vfiller);
          size  -= 64ULL;
          p_dst += 64ULL;
      }

      while(size >= 32ULL)
      {
          _mm512_store_pd(p_dst+0ULL, vfiller);
          _mm512_store_pd(p_dst+8ULL, vfiller);
          _mm512_store_pd(p_dst+16ULL, vfiller);
          _mm512_store_pd(p_dst+24ULL, vfiller);

          size  -= 32ULL;
          p_dst += 32ULL;
      }

      while(size >= 16ULL)
      {
          _mm512_store_pd(p_dst+0ULL, vfiller);
          _mm512_store_pd(p_dst+8ULL, vfiller);

          size  -= 16ULL;
          p_dst += 16ULL;
      }

      while(size >= 8ULL)
      {
          _mm512_store_pd(p_dst+0ULL, vfiller);

          size  -= 8ULL;
          p_dst += 8ULL;
      }

      while(size)
      {
          *p_dst = cfiller;
          size--;
          p_dst++;
      }
}