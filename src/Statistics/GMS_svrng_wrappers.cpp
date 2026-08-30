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
#include <svrng.h>
#include "GMS_svrng_wrappers.h"

#if !defined(GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CASE_NO_DISTR)
     GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CASE_NO_DISTR          \
                 svrng_engine_t engine;                          \
		 int32_t err = -9999;                            \
                 uint32_t seed = 0U;                             \
		 int32_t result = -9999;                         \
		 result = _rdrand32_step(&seed);                 \
		 if(!result) seed = 154625984U;                  \
		 engine = svrng_new_mt19937_engine(seed);        \
		 err = svrng_get_status();                       \
		 if(err != SVRNG_STATUS_OK) {                    \
                    status = err;                                \
		    return;                                      \
		 }
#endif

#if !defined(GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CMPLX_CASE_NO_DISTR)
    GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CMPLX_CASE_NO_DISTR        \
                 svrng_engine_t re_eng,im_eng;                      \
                 int32_t err = -9999;                               \
		 uint32_t seedre = 0U;                              \
		 uint32_t seedim = 0U;                              \
		 int32_t result = -9999;                            \
		 result = _rdrand32_step(&seedre);                  \
		 if(!result) seedre = 125654897U;                   \
		 re_eng = svrng_new_mt19937_engine(seedre);         \
		 result = _rdrand32_step(&seedim);                  \
		 if(!result) seedim = 256987415U;                   \
		 im_eng = svrng_new_mt19937_engine(seedim);         \
		 err = svrng_get_status();                          \
		 if(err != SVRNG_STATUS_OK) {                       \
                    status = err;                                   \
		    return;                                         \
		 }
#endif

#if !defined(GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CASE_NORMAL)
    GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CASE_NORMAL              \
                 svrng_engine_t engine;                           \
                 svrng_distribution_t normal;                     \
		 int32_t err = -9999;                             \
		 uint32_t seed = 0U;                              \
		 int32_t result = -9999;                          \
		 result = _rdrand32_step(&seed);                  \
		 if(!result) seed = 256984512U;                   \
		 engine  = svrng_new_mt19937_engine(seed);        \
		 normal = svrng_new_normal_distribution(lo,hi);   \
		 err = svrng_get_status();                        \
		 if(err != SVRNG_STATUS_OK) {                     \
                    status = err;                                 \
		    return;                                       \
		 }
#endif

#if !defined(GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CMPLX_CASE_NORMAL)
    GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CMPLX_CASE_NORMAL              \
                 svrng_engine_t re_eng,im_eng;                          \
		 svrng_distribution_t re_norm,im_norm;                  \ 
		 int32_t err = -9999;                                   \
		 uint32_t seedre = 0U;                                  \ 
		 uint32_t seedim = 0U;                                  \
		 int32_t result = -9999;                                \
		 result = _rdrand32_step(&seedre);                      \
                 if(!result) seedre = 235698756U;                       \
		 re_eng = svrng_new_mt19937_engine(seedre);             \
		 re_norm = svrng_new_normal_distribution(relo,rehi);    \
		 result = _rdrand32_step(&seedim);                      \
		 if(!result) seedim = 112565498U;                       \
		 im_eng = svrng_new_mt19937_engine(seedim);             \
		 im_norm = svrng_new_normal_distribution(imlo,imhi);    \
		 err = svrng_get_status();                              \
		 if(err != SVRNG_STATUS_OK) {                           \
                    status = err;                                       \
		    return;                                             \
		 }
#endif

#if !defined(GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CASE_UNIFORM)
    GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CASE_UNIFORM             \
                 svrng_engine_t engine;                          \
                 svrng_distribution_t uniform;                    \
		 int32_t err = -9999;                             \
		 uint32_t seed = 0U;                              \
		 int32_t result = -9999;                          \
		 result = _rdrand32_step(&seed);                  \
		 if(!result) seed = 612984845U;                   \
		 engine  = svrng_new_mt19937_engine(seed);        \
		 uniform = svrng_new_uniform_distribution(lo,hi); \
		 err = svrng_get_status();                        \
		 if(err != SVRNG_STATUS_OK) {                     \
                    status = err;                                 \
		    return;                                       \
		 }
#endif

#if !defined(GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CMPLX_CASE_UNIFORM)
    GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CMPLX_CASE_UNIFORM             \
                 svrng_engine_t re_eng,im_eng;                          \
		 svrng_distribution_t re_uni,im_uni;                    \ 
		 int32_t err = -9999;                                   \
		 uint32_t seedre = 0U;                                  \ 
		 uint32_t seedim = 0U;                                  \
		 int32_t result = -9999;                                \
		 result = _rdrand32_step(&seedre);                      \
                 if(!result) seedre = 235698756U;                       \
		 re_eng = svrng_new_mt19937_engine(seedre);             \
		 re_uni = svrng_new_uniform_distribution(relo,rehi);    \
		 result = _rdrand32_step(&seedim);                      \
		 if(!result) seedim = 112565498U;                       \
		 im_eng = svrng_new_mt19937_engine(seedim);             \
		 im_uni = svrng_new_uniform_distribution(imlo,imhi);    \
		 err = svrng_get_status();                              \
		 if(err != SVRNG_STATUS_OK) {                           \
                    status = err;                                       \
		    return;                                             \
		 }


void
gms::math::stat::
svrng_wrapper_mt19937_init_float8(float * __restrict data,
				  const in32_t length, // must have a length Mod(len,8) == 0
				  const float lo,
				  const float hi,
				  const int32_t type,
				  int32_t & status) {
     switch(type) {

         case 0: { // no distribution
	         GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CASE_NO_DISTR
		 for(int32_t i = 0; i != length; i += 8) {
                     *((svrng_float8_t*)&data[i]) =
		            svrng_generate8_float(engine,NULL);
		 }
		 err = -9999;
		 err = svrng_get_status();
		 if(err != SVRNG_STATUS_OK) {
		    svrng_delete_engine(engine);
                    status = err;
		    return;
		 }
		 svrng_delete_engine(engine);
		 break;
	   }
	 case 1: {  // normal distribution
                 GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CASE_NORMAL
		 for(int32_t i = 0; i != length; i += 8) {
                     *((svrng_float8_t)&data[i]) =
		            svrng_generate8_float(engine,normal);
		 }
		 err = -9999;
		 err = svrng_get_status();
		 if(err != SVRNG_STATUS_OK) {
                    svrng_delete_distribution(normal);
		    svrng_delete_engine(engine);
		    status = err;
		    return;
		 }
		  svrng_delete_distribution(normal);
		  svrng_delete_engine(engine);
		  break;
	   }
	 case 2: { // uniform distribution
                 GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CASE_UNIFORM
		 for(int32_t i = 0; i != length; i += 8) {
                     *((svrng_float8_t)&data[i]) =
		            svrng_generate8_float(engine,normal);
		 }
		 err = -9999;
		 err = svrng_get_status();
		 if(err != SVRNG_STATUS_OK) {
                    svrng_delete_distribution(uniform);
		    svrng_delete_engine(engine);
		    status = err;
		    return;
		 }
		  svrng_delete_distribution(uniform);
		  svrng_delete_engine(engine);
		  break; 

	   }
	  default : {
                      status = -99; // invalid switch argument
		      return;
	  }
     }
}

void
gms::math::stat::
svrng_wrapper_mt19937_init_double4(double * __restrict data,
                                   const int32_t length,
				   const double lo,
				   const double hi,
				   const int32_t type,
				   int32_t & status) {
     switch(type) {

         case 0: { // no distribution
              GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CASE_NO_DISTR
	      for(int32_t i = 0; i != length; i += 4) {
                     *((svrng_double4_t*)&data[i]) =
		            svrng_generate4_double(engine,NULL);
		 }
	       err = -9999;
	       err = svrng_get_status();
	       if(err != SVRNG_STATUS_OK) {
		    svrng_delete_engine(engine);
                    status = err;
		    return;
	       }
	       svrng_delete_engine(engine);
	       break;
	  }
	 case 1: { // normal distribution
              GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CASE_NORMAL
	      for(int32_t i = 0; i != length; i += 4) {
                     *((svrng_double4_t)&data[i]) =
		            svrng_generate4_double(engine,normal);
	       }
	       err = -9999;
	       err = svrng_get_status();
	       if(err != SVRNG_STATUS_OK) {
                    svrng_delete_distribution(normal);
		    svrng_delete_engine(engine);
		    status = err;
		    return;
		 }
		svrng_delete_distribution(normal);
		svrng_delete_engine(engine);
		break;
	  }
	case 2: { // uniform distribution
             GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CASE_UNIFORM
	     for(int32_t i = 0; i != length; i += 4) {
                     *((svrng_double4_t)&data[i]) =
		            svrng_generate4_double(engine,normal);
	      }
	      err = -9999;
	      err = svrng_get_status();
	      if(err != SVRNG_STATUS_OK) {
                    svrng_delete_distribution(uniform);
		    svrng_delete_engine(engine);
		    status = err;
		    return;
	       }
	       svrng_delete_distribution(uniform);
	       svrng_delete_engine(engine);
	       break; 
	 }
       default : {
                    status = -99;
		    return;
         }
    }
}

void
gms::math::stat::
svrng_wrapper_mt19937_init_avxvec8(AVXVec8 * __restrict data,
                                   const int64_t length,
				   const float lo,
				   const float hi,
				   const int32_t type,
				   int32_t & status) {
     switch(type) {

         case 0: {
             GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CASE_NO_DISTR
	     svrng_float8_t rand;
	     for(int_t i = 0; i != length; ++i) {
                 rand = svrng_generate8_float(engine,NULL);
		 data[i] = *(AVXVec8*)&rand;
	     }
             err = -9999;
	     err = svrng_get_status();
	     if(err != SVRNG_STATUS_OK) {
		    svrng_delete_engine(engine);
                    status = err;
		    return;
	       }
	      svrng_delete_engine(engine);
	      break;
          }
	case 1: {
            GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CASE_NORMAL
	    svrng_float8_t rand;
	    for(int32_t i = 0; i != length; ++i) {
                rand = svrng_generate8_float(engine,normal);
		data[i] = *(AVXVec8*)&rand;
	    }
	    err = -9999;
	    err = svrng_get_status();
	    if(err != SVRNG_STATUS_OK) {
                    svrng_delete_distribution(normal);
		    svrng_delete_engine(engine);
		    status = err;
		    return;
	     }
	     svrng_delete_distribution(normal);
	     svrng_delete_engine(engine);
	     break;
	  }
	case 2: {
            GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CASE_UNIFORM
	    svrng_float8_t rand;
	    for(int32_t i = 0; i != length; ++i) {
                rand = svrng_generate8_float(engine,uniform);
		data[i] = *(AVXVec8*)&rand;
	    }
	    err = -9999;
	    err = svrng_get_status();
	    if(err != SVRNG_STATUS_OK) {
                    svrng_delete_distribution(uniform);
		    svrng_delete_engine(engine);
		    status = err;
		    return;
	     }
	     svrng_delete_distribution(uniform);
	     svrng_delete_engine(engine);
	     break; 
	}
       default : {
                   statuse = -99;
		   return;
          }
     }  
}

void
gms::math::stat::
svnrg_wrapper_mt19937_init_avx512c4f32(AVX512c4f32 * __restrict data,
                                       const int32_t length,
				       const float relo,
				       const float rehi,
				       const float imlo,
				       const float imhi,
				       const int32_t type,
				       int32_t & status) {

     switch(type) {

         case 0: {
              GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CMPLX_CASE_NO_DISTR
	      svrng_float16_t re_rand;
	      svrng_float16_t im_rand;
	      for(int32_t i = 0; i != length; ++i) {
                  re_rand = svrng_generate16_float(re_eng,NULL);
		  data[i].m_re = *(__m512*)&re_rand;
		  im_rand = svrng_generate16_float(im_eng,NULL);
		  data[i].m_im = *(__m512*)&im_rand;
	      }
	      err = -9999;
	      err = svrng_get_status();
	      if(err != SVRNG_STATUS_OK) {
                 svrng_delete_engine(re_eng);
		 svrng_delete_engine(im_eng);
		 status = err;
		 return;
	      }
	      svrng_delete_engine(re_eng);
	      svrng_delete_engine(im_eng);
	      break;
	  }
        case 1: {
             GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CMPLX_CASE_NORMAL
	     svrng_float16_t re_rand;
	     svrng_float16_t im_rand;
	     for(int32_t i = 0; i != length; ++i) {
                 re_rand = svrng_generate16_float(re_eng,re_norm);
		 data[i].m_re = *(__512*)&re_rand;
		 im_rand = svrng_generate16_float(im_eng,im_norm);
		 data[i].m_im = *(__m512*)&im_rand;
	     }
	     err = -9999;
	     err = svrng_get_status();
	     if(err != SVRNG_STATUS_OK) {
                svrng_delete_engine(re_eng);
		svrng_delete_distribution(re_norm);
		svrng_delete_engine(im_eng);
		svrng_delete_distribution(im_norm);
		status = err;
		return;
	     }
	     svrng_delete_engine(re_eng);
	     svrng_delete_distribution(re_norm);
	     svrng_delete_engine(im_eng);
	     svrn_delete_distribution(im_norm);
	     break;
	  }
	case 2: {
             GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CMPLX_CASE_UNIFORM
	     svrng_float16_t re_rand;
	     svrng_float16_t im_rand;
	     for(int32_t i = 0; i != length; ++i) {
                 re_rand = svrng_generate16_float(re_eng,re_uni);
		 data[i].m_re = *(__m512*)&re_rand;
		 im_rand = svrng_generate16_float(im_eng,im_uni);
		 data[i].m_im = *(__m512*)&im_rand;
	     }
	     err = -9999;
	     err = svrng_get_status();
	     if(err != SVRNG_STATUS_OK) {
                svrng_delete_engine(re_eng);
		svrng_delete_distribution(re_uni);
		svrng_delete_engine(im_eng);
		svrng_delete_distribution(im_uni);
		status = err;
		return;
	     }
	     svrng_delete_engine(re_eng);
	     svrng_delete_distribution(re_uni);
	     svrng_delete_engine(im_eng);
	     svrng_delete_distribution(im_uni);
	     break;
	}
      default : {
                   status = -99;
		   return;
        }
    }
}

void
gms::math::stat::
svrng_wrapper_mt19937_init_avxc8f32(AVXc8f32 * __restrict data,
				    const int32_t length,
				    const float   relo,
				    const float   rehi,
				    const float   imlo,
				    const float   imhi,
				    const int32_t type,
				    int32_t &     status)  {

      switch(type) {

         case 0: {
              GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CMPLX_CASE_NO_DISTR
              svrng_float8_t re_rand;
	      svrng_float8_t im_rand;
	      for(int32_t i = 0; i != length; ++i) {
                  re_rand = svrng_generate8_float(re_eng,NULL);
		  data[i].m_re = *(__m256*)&re_rand;
		  im_rand = svrng_generate8_float(im_eng,NULL);
		  data[i].m_im = *(__m256*)&im_rand;
	      }
	      err = -9999;
	      err = svrng_get_status();
	      if(err != SVRNG_STATUS_OK) {
                 svrng_delete_engine(re_eng);
		 svrng_delete_engine(im_eng);
		 status = err;
		 return;
	      }
	      svrng_delete_engine(re_eng);
	      svrng_delete_engine(im_eng);
	      break;
	  }
        case 1: {
             GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CMPLX_CASE_NORMAL
	     svrng_float8 re_rand;
	     svrng_float8 im_rand;
             for(int32_t i = 0; i != length; ++i) {
                 re_rand = svrng_generate8_float(re_eng,re_norm);
                 data[i].m_re = *(__m256*)&re_rand;
		 im_rand = svrng_generate8_float(im_eng,im_norm);
		 data[i].m_im = *(__m256*)&im_rand;
	     }
             err = -9999;
	     err = svrng_get_status();
	     if(err != SVRNG_STATUS_OK) {
                svrng_delete_engine(re_eng);
		svrng_delete_distribution(re_norm);
		svrng_delete_engine(im_eng);
		svrng_delete_distribution(im_norm);
		status = err;
		return;
	     }
	     svrng_delete_engine(re_eng);
	     svrng_delete_distribution(re_norm);
	     svrng_delete_engine(im_eng);
	     svrn_delete_distribution(im_norm);
	     break;
	  }
       case 2: {
            GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CMPLX_CASE_UNIFORM
	    svrng_float8 re_rand;
	    svrng_float8 im_rand;
            for(int32_t i = 0; i != length; ++i) {
                re_rand = svrng_generate8_float(re_eng,re_uni);
                data[i].m_re = *(__m256*)&re_rand;
		im_rand = svrng_generate8_float(im_eng,im_uni);
		data[i].m_im = *(__m256*)&im_rand;
	    }
	     err = -9999;
	     err = svrng_get_status();
	     if(err != SVRNG_STATUS_OK) {
                svrng_delete_engine(re_eng);
		svrng_delete_distribution(re_uni);
		svrng_delete_engine(im_eng);
		svrng_delete_distribution(im_uni);
		status = err;
		return;
	     }
	     svrng_delete_engine(re_eng);
	     svrng_delete_distribution(re_uni);
	     svrng_delete_engine(im_eng);
	     svrng_delete_distribution(im_uni);
	     break;
         }
       default : {
                     status = -99;
		     return;
           } 
    }

}

void
gms::math::stat::
svrng_wrapper_mt19937_init_avxc4f64(AVXc4f64 * __restrict data,
                                    const int32_t length,
                                    const double relo,
				    const double rehi,
				    const double imlo,
				    const double imhi,
				    const int32_t type,
				    int32_t & status) {

        switch(type) {

             case 0: {
                  GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CMPLX_CASE_NO_DISTR
                  svrng_double4_t re_rand;
                  svrng_double4_t im_rand;
                  __assume_aligned(data,64);
		  for(int_t i = 0; i != length; ++i) {
                      re_rand = svrng_generate4_double(re_eng,NULL);
                      data[i].m_re = *(__m256d*)&re_rand;
		      im_rand = svrng_generate4_double(im_eng,NULL);
		      data[i].m_im = *(__m256*)&im_rand;
		  }
		  err = -9999;
	          err = svrng_get_status();
	          if(err != SVRNG_STATUS_OK) {
                       svrng_delete_engine(re_eng);
		       svrng_delete_engine(im_eng);
		       status = err;
		       return;
	          }
	          svrng_delete_engine(re_eng);
	          svrng_delete_engine(im_eng);
	          break;
	     }
          case 1: {
               GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CMPLX_CASE_NORMAL
               svrng_double4_t re_rand;
	       svrng_double4_t im_rand;
	       __assume_aligned(data,64);
	       for(int_t i = 0; i != length; ++i) {
                   re_rand = svrng_generate4_double(re_eng,re_norm);
		   data[i].m_re = *(__m256d*)&re_rand;
		   im_rand = svrng_generate4_double(im_eng,im_norm);
		   data[i].m_im = *(__m256d*)&im_rand;
                }
	        err = -9999;
	        err = svrng_get_status();
	        if(err != SVRNG_STATUS_OK) {
                    svrng_delete_engine(re_eng);
		    svrng_delete_distribution(re_norm);
		    svrng_delete_engine(im_eng);
		    svrng_delete_distribution(im_norm);
		    status = err;
		    return;
	         }
	         svrng_delete_engine(re_eng);
	         svrng_delete_distribution(re_norm);
	         svrng_delete_engine(im_eng);
	         svrn_delete_distribution(im_norm);
	         break;
	  }
	case 2: {
             GMS_SVRNG_WRAPPERS_MT19937_FUNC_BODY_CMPLX_CASE_UNIFORM
             svrng_double_t re_rand;
	     svrng_double_t im_rand;
	     __assume_aligned(data,64);
	     for(int32_t i = 0; i != length; ++i) {
                 re_rand = svrng_generate4_double(re_eng,re_uni);
                 data[i].m_re = *(__m256d*)&re_rand;
		 im_rand = svrng_generate4_double(im_eng,im_uni);
		 data[i].m_im = *(__m256d*)&im_rand;
	     }
	     err = -9999;
	     err = svrng_get_status();
	     if(err != SVRNG_STATUS_OK) {
                svrng_delete_engine(re_eng);
		svrng_delete_distribution(re_uni);
		svrng_delete_engine(im_eng);
		svrng_delete_distribution(im_uni);
		status = err;
		return;
	     }
	     svrng_delete_engine(re_eng);
	     svrng_delete_distribution(re_uni);
	     svrng_delete_engine(im_eng);
	     svrng_delete_distribution(im_uni);
	     break;
         }
       default : {
                     status = -99;
		     return;
           } 
     }

}

