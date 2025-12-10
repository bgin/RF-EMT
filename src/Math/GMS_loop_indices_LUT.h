



#ifndef __GMS_LOOP_INDICES_LUT_H__
#define __GMS_LOOP_INDICES_LUT_H__ 031220250814 

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include "GMS_config.h"

namespace gms 
{

namespace math 
{

/*
     Loop indices LUT ranging from 0,...2256
*/
extern const float LUT_loop_indices_2257_align16[2257];

extern const float LUT_loop_indices_2257_align32[2257];

extern const float LUT_loop_indices_2257_align64[2257];

extern const float LUT_loop_indices_1000000_align64[1000000];

__ATTR_ALWAYS_INLINE__
static inline 
std::int32_t 
create_LUT_loop_indices(const std::size_t niters,
                        const char * fname)
{
     if(__builtin_expect(niters<1ull,0)) { return (-1);}
     if(__builtin_expect(fname==NULL,0)) { return (-2);}
     FILE * fp{NULL};
     fp = fopen(fname,"w+");
     if(!fp)
     {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           return (-3);
     }
     for(std::size_t i{0ull}; i != niters; ++i)
     {
           fprintf(fp,"%.7ff,\n", static_cast<float>(i));
    }
    if(fp) fclose(fp);
    return (0);
}



}

}


#endif /*__GMS_LOOP_INDICES_LUT_H__*/