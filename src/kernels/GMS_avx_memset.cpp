


#include <immintrin.h>
#include "GMS_avx_memset.h"

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
void 
gms::common::
avx_memset_unroll8x_ps(float * __restrict__ dst,const float filler,std::size_t size)
{
   
    const __m256 vfiller{_mm256_set1_ps(filler)};
    float * __restrict__ p_dst{dst};
    const float cfiller{filler};

    while(((uintptr_t)p_dst & 31) && size)
    {
        *p_dst++ = cfiller;
         --size;
    }

    while(size >= 64ULL)
    {
         
        _mm256_storeu_ps(p_dst+0, vfiller);
        _mm256_storeu_ps(p_dst+8, vfiller);
        _mm256_storeu_ps(p_dst+16, vfiller);
        _mm256_storeu_ps(p_dst+24,vfiller);
        _mm256_storeu_ps(p_dst+32,vfiller);
        _mm256_storeu_ps(p_dst+40,vfiller);
        _mm256_storeu_ps(p_dst+48,vfiller);
        _mm256_storeu_ps(p_dst+56,vfiller);
         
        size  -= 64ULL;
        p_dst += 64ULL;
    }

    while(size >= 48ULL)
    {
        _mm256_storeu_ps(p_dst+0, vfiller);
        _mm256_storeu_ps(p_dst+8, vfiller);
        _mm256_storeu_ps(p_dst+16, vfiller);
        _mm256_storeu_ps(p_dst+24,vfiller);
        _mm256_storeu_ps(p_dst+32,vfiller);
        _mm256_storeu_ps(p_dst+40,vfiller);

        size  -= 48ULL;
        p_dst += 48ULL;
    }

    while(size >= 32ULL)
    {
        _mm256_storeu_ps(p_dst+0, vfiller);
        _mm256_storeu_ps(p_dst+8, vfiller);
        _mm256_storeu_ps(p_dst+16, vfiller);
        _mm256_storeu_ps(p_dst+24,vfiller);

        size  -= 32ULL;
        p_dst += 32ULL;
    }

    while(size >= 16ULL)
    {
        _mm256_storeu_ps(p_dst+0, vfiller);
        _mm256_storeu_ps(p_dst+8, vfiller);

        size  -= 16ULL;
        p_dst += 16ULL;
    }

    while(size >= 8ULL)
    {
        _mm256_storeu_ps(p_dst+0, vfiller);

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

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
void
gms::common::
avx_memset_unroll16x_ps(float * __restrict__ dst,const float filler,std::size_t size)
{
    const __m256 vfiller{_mm256_set1_ps(filler)};
    float * __restrict__ p_dst = dst;
    std::size_t __i;
    const float cfiller{filler};

   
    for(__i = 0ull; (__i+127ull) < size; __i += 128ull)
    {
          _mm256_storeu_ps(&p_dst[__i+0ull],  vfiller);
          _mm256_storeu_ps(&p_dst[__i+8ull],  vfiller);
          _mm256_storeu_ps(&p_dst[__i+16ull],  vfiller);
          _mm256_storeu_ps(&p_dst[__i+24ull], vfiller);
          _mm256_storeu_ps(&p_dst[__i+32ull], vfiller);
          _mm256_storeu_ps(&p_dst[__i+40ULL], vfiller);
          _mm256_storeu_ps(&p_dst[__i+48ULL], vfiller);
          _mm256_storeu_ps(&p_dst[__i+56ULL], vfiller);
          _mm256_storeu_ps(&p_dst[__i+64ULL], vfiller);
          _mm256_storeu_ps(&p_dst[__i+72ULL], vfiller);
          _mm256_storeu_ps(&p_dst[__i+80ULL], vfiller);
          _mm256_storeu_ps(&p_dst[__i+88ULL], vfiller);
          _mm256_storeu_ps(&p_dst[__i+96ULL], vfiller);
          _mm256_storeu_ps(&p_dst[__i+104ULL], vfiller);
          _mm256_storeu_ps(&p_dst[__i+112ULL], vfiller);
          _mm256_storeu_ps(&p_dst[__i+120ULL], vfiller);
    }

    for(; (__i+95ull) < size;  __i += 96ull)
    {
          _mm256_storeu_ps(&p_dst[__i+0ull],  vfiller);
          _mm256_storeu_ps(&p_dst[__i+8ull],  vfiller);
          _mm256_storeu_ps(&p_dst[__i+16ull],  vfiller);
          _mm256_storeu_ps(&p_dst[__i+24ull], vfiller);
          _mm256_storeu_ps(&p_dst[__i+32ull], vfiller);
          _mm256_storeu_ps(&p_dst[__i+40ULL], vfiller);
          _mm256_storeu_ps(&p_dst[__i+48ULL], vfiller);
          _mm256_storeu_ps(&p_dst[__i+56ULL], vfiller);
          _mm256_storeu_ps(&p_dst[__i+64ULL], vfiller);
          _mm256_storeu_ps(&p_dst[__i+72ULL], vfiller);
          _mm256_storeu_ps(&p_dst[__i+80ULL], vfiller);
          _mm256_storeu_ps(&p_dst[__i+88ULL], vfiller);
    }

    for(; (__i+63ull) < size; __i += 64ull)
    {
          _mm256_storeu_ps(&p_dst[__i+0ull],  vfiller);
          _mm256_storeu_ps(&p_dst[__i+8ull],  vfiller);
          _mm256_storeu_ps(&p_dst[__i+16ull],  vfiller);
          _mm256_storeu_ps(&p_dst[__i+24ull], vfiller);
          _mm256_storeu_ps(&p_dst[__i+32ull], vfiller);
          _mm256_storeu_ps(&p_dst[__i+40ULL], vfiller);
          _mm256_storeu_ps(&p_dst[__i+48ULL], vfiller);
          _mm256_storeu_ps(&p_dst[__i+56ULL], vfiller);
    }

    for(; (__i+47ull) < size; __i += 48ull)
    {
          _mm256_storeu_ps(&p_dst[__i+0ull],  vfiller);
          _mm256_storeu_ps(&p_dst[__i+8ull],  vfiller);
          _mm256_storeu_ps(&p_dst[__i+16ull],  vfiller);
          _mm256_storeu_ps(&p_dst[__i+24ull], vfiller);
          _mm256_storeu_ps(&p_dst[__i+32ull], vfiller);
          _mm256_storeu_ps(&p_dst[__i+40ULL], vfiller);
    }

    for(; (__i+31ull) < size; __i += 32ull)
    {
          _mm256_storeu_ps(&p_dst[__i+0ull],  vfiller);
          _mm256_storeu_ps(&p_dst[__i+8ull],  vfiller);
          _mm256_storeu_ps(&p_dst[__i+16ull],  vfiller);
          _mm256_storeu_ps(&p_dst[__i+24ull], vfiller);
    }

    for(; (__i+15ull) < size; __i += 16ull)
    {
           _mm256_storeu_ps(&p_dst[__i+0ull],  vfiller);
           _mm256_storeu_ps(&p_dst[__i+8ull],  vfiller);
    }

    for(; (__i+7ull) < size; __i += 8ull)
    {
           _mm256_storeu_ps(&p_dst[__i+0ull],  vfiller);
    }

    for(; (__i+0ull) < size; __i += 1ull)
    {
          p_dst[__i] = cfiller;
    }

    
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
void 
gms::common::
avx_memset_unroll8x_pd(double * __restrict__ dst,const double filler,std::size_t size)
{
    const __m256d vfiller{_mm256_set1_pd(filler)};
    double * __restrict__ p_dst{dst};
    const double cfiller{filler};

    while(((uintptr_t)p_dst & 31) && size)
    {
        *p_dst++ = cfiller;
         --size;
    }

    while(size >= 32ULL)
    {
         
        _mm256_store_pd(p_dst+0ULL, vfiller);
        _mm256_store_pd(p_dst+4ULL, vfiller);
        _mm256_store_pd(p_dst+8ULL, vfiller);
        _mm256_store_pd(p_dst+12ULL, vfiller);
        _mm256_store_pd(p_dst+16ULL, vfiller);
        _mm256_store_pd(p_dst+20ULL,vfiller);
        _mm256_store_pd(p_dst+24ULL,vfiller);
        _mm256_store_pd(p_dst+28ULL,vfiller);
         
        p_dst += 32ULL;
        size  -= 32ULL;
        
    }

    while(size >= 24ULL)
    {
        _mm256_store_pd(p_dst+0ULL, vfiller);
        _mm256_store_pd(p_dst+4ULL, vfiller);
        _mm256_store_pd(p_dst+8ULL, vfiller);
        _mm256_store_pd(p_dst+12ULL, vfiller);
        _mm256_store_pd(p_dst+16ULL, vfiller);
        _mm256_store_pd(p_dst+20ULL,vfiller);
        
        p_dst += 24ULL;
        size  -= 24ULL;
        
    }

    while(size >= 16ULL)
    {
        _mm256_store_pd(p_dst+0ULL, vfiller);
        _mm256_store_pd(p_dst+4ULL, vfiller);
        _mm256_store_pd(p_dst+8ULL, vfiller);
        _mm256_store_pd(p_dst+12ULL, vfiller);
        
        p_dst += 16ULL;
        size  -= 16ULL;
       
    }

    while(size >= 8ULL)
    {
        _mm256_store_pd(p_dst+0ULL, vfiller);
        _mm256_store_pd(p_dst+4ULL, vfiller);
        
        p_dst += 8ULL;
        size  -= 8ULL;
        
    }

    while(size >= 4ULL)
    {
        _mm256_store_pd(p_dst+0ULL, vfiller);
        
        p_dst += 4ULL;
        size  -= 4ULL;
        
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
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
void
gms::common::
avx_memset_unroll16x_pd(double * __restrict__ dst,const double filler,std::size_t size)
{
     const __m256d vfiller{_mm256_set1_pd(filler)};
     double * __restrict__ p_dst{dst};
     double cfiller{filler};

     while(((uintptr_t)p_dst & 31) && size)
     {
          *p_dst++ = cfiller;
          --size;
     }

     while(size >= 64ULL)
     {
          _mm256_store_pd(p_dst+0ULL, vfiller);
          _mm256_store_pd(p_dst+4ULL, vfiller);
          _mm256_store_pd(p_dst+8ULL, vfiller);
          _mm256_store_pd(p_dst+12ULL, vfiller);
          _mm256_store_pd(p_dst+16ULL, vfiller);
          _mm256_store_pd(p_dst+20ULL,vfiller);
          _mm256_store_pd(p_dst+24ULL,vfiller);
          _mm256_store_pd(p_dst+28ULL,vfiller);
          _mm256_store_pd(p_dst+32ULL,vfiller);
          _mm256_store_pd(p_dst+36ULL,vfiller);
          _mm256_store_pd(p_dst+40ULL,vfiller);
          _mm256_store_pd(p_dst+44ULL,vfiller);
          _mm256_store_pd(p_dst+48ULL,vfiller);
          _mm256_store_pd(p_dst+52ULL,vfiller);
          _mm256_store_pd(p_dst+56ULL,vfiller);
          _mm256_store_pd(p_dst+60ULL,vfiller);

          size  -= 64ULL;
          p_dst += 64ULL;
      }

      while(size >= 48ULL)
      {
          _mm256_store_pd(p_dst+0ULL, vfiller);
          _mm256_store_pd(p_dst+4ULL, vfiller);
          _mm256_store_pd(p_dst+8ULL, vfiller);
          _mm256_store_pd(p_dst+12ULL, vfiller);
          _mm256_store_pd(p_dst+16ULL, vfiller);
          _mm256_store_pd(p_dst+20ULL,vfiller);
          _mm256_store_pd(p_dst+24ULL,vfiller);
          _mm256_store_pd(p_dst+28ULL,vfiller);
          _mm256_store_pd(p_dst+32ULL,vfiller);
          _mm256_store_pd(p_dst+36ULL,vfiller);
          _mm256_store_pd(p_dst+40ULL,vfiller);
          _mm256_store_pd(p_dst+44ULL,vfiller);

          size  -= 48ULL;
          p_dst += 48ULL;
      }

      while(size >= 32ULL)
      {
          _mm256_store_pd(p_dst+0ULL, vfiller);
          _mm256_store_pd(p_dst+4ULL, vfiller);
          _mm256_store_pd(p_dst+8ULL, vfiller);
          _mm256_store_pd(p_dst+12ULL, vfiller);
          _mm256_store_pd(p_dst+16ULL, vfiller);
          _mm256_store_pd(p_dst+20ULL,vfiller);
          _mm256_store_pd(p_dst+24ULL,vfiller);
          _mm256_store_pd(p_dst+28ULL,vfiller);

          size  -= 32ULL;
          p_dst += 32ULL;
      }

      while(size >= 16ULL)
      {
          _mm256_store_pd(p_dst+0ULL, vfiller);
          _mm256_store_pd(p_dst+4ULL, vfiller);
          _mm256_store_pd(p_dst+8ULL, vfiller);
          _mm256_store_pd(p_dst+12ULL, vfiller);

          size  -= 16ULL;
          p_dst += 16ULL;
      }

      while(size >= 8ULL)
      {
          _mm256_store_pd(p_dst+0ULL, vfiller);
          _mm256_store_pd(p_dst+4ULL, vfiller);

          size  -= 8ULL;
          p_dst += 8ULL;
      }

      while(size >= 4ULL)
      {
          _mm256_store_pd(p_dst+0ULL, vfiller);

          size  -= 4ULL;
          p_dst += 4ULL;
      }

      while(size)
      {
          *p_dst = cfiller;
          size--;
          p_dst++;
      }
}