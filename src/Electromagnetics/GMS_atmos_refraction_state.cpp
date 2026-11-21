
#include <omp.h>
#include <iostream>
#include "GMS_atmos_refraction_state.h"

gms::atmos
::atmos_refraction_state_t
::atmos_refraction_state_t(const atmos_refraction_size_params_t &sz_params,
                           const bool do_omp_alloc)
:
m_n_nidx{sz_params.n_nidx},
m_n_n0idx{sz_params.n_n0idx},
m_n_z{sz_params.n_z},
m_n_z0{sz_params.n_z0},
m_n_r{sz_params.n_r},
m_n_R0{sz_params.n_R0},
m_n_phi{sz_params.n_phi},
m_n_phi0{sz_params.n_phi0},
m_n_ntht{sz_params.n_ntht},
m_n_nphi{sz_params.n_nphi},
m_n_dndr{sz_params.n_dndr},
m_n_rho{sz_params.n_rho},
m_n_dn0{sz_params.n_dn0},
m_n_beta{sz_params.n_beta},
m_n_navgh{sz_params.n_navgh},
m_n_H{sz_params.n_H},
m_n_f345{sz_params.n_f345},
m_n_f351{sz_params.n_f351},
m_n_f352{sz_params.n_f352},
m_n_nh{sz_params.n_nh},
m_n_f41{sz_params.n_f41},
m_n_f{sz_params.n_f},
m_n_d{sz_params.n_d},
m_n_Nmf{sz_params.n_Nmf},
m_n_f412{sz_params.n_f412},
m_n_f413{sz_params.n_f413},
m_n_D1{sz_params.n_D1},
m_n_f415{sz_params.n_f415},
m_n_f425{sz_params.n_f425},
m_n_f428{sz_params.n_f428},
m_n_f429{sz_params.n_f429},
m_n_H1{sz_params.n_H1},
m_n_H2{sz_params.n_H2},
m_n_f430{sz_params.n_f430},
m_n_H3{sz_params.n_H3},
m_n_deln0{sz_params.n_deln0},
m_n_f431{sz_params.n_f431},
m_n_f438{sz_params.n_f438},
m_n_f442{sz_params.n_f442},
m_n_f445{sz_params.n_f445},
m_n_g{sz_params.n_g},
m_n_f450{sz_params.n_f450},
m_n_f451{sz_params.n_f451},
m_n_Hc0{sz_params.n_Hc0},
m_n_delnA{sz_params.n_delnA},
m_n_na{sz_params.n_na},
m_n_nc{sz_params.n_nc},
m_n_f53{sz_params.n_f53},
m_n_f517{sz_params.n_f517},
m_n_tht{sz_params.n_tht},
m_n_f531{sz_params.n_f531},
m_n_thtc{sz_params.n_thtc},
m_n_f534{sz_params.n_f534},
m_n_Rc{sz_params.n_Rc},
m_n_f535{sz_params.n_f535},
m_n_f538{sz_params.n_f538},
m_n_f539{sz_params.n_f539},
m_n_f541{sz_params.n_f541},
m_n_H0{sz_params.n_H0},
m_n_Hc{sz_params.n_Hc},
m_n_f543{sz_params.n_f543},
m_n_H10{sz_params.n_H10},
m_n_f554{sz_params.n_f554},
m_n_H20{sz_params.n_H20},
m_n_f572{sz_params.n_f572},
m_n_f576{sz_params.n_f576},
m_n_f577{sz_params.n_f577},
m_n_f578{sz_params.n_f578},
m_n_f579{sz_params.n_f579},
m_n_f590{sz_params.n_f590},
m_n_R2{sz_params.n_R2},
m_n_f591{sz_params.n_f591},
m_n_f593{sz_params.n_f593},
m_n_f595{sz_params.n_f595},
m_n_HB{sz_params.n_HB},
m_n_f62{sz_params.n_f62},
m_n_HC2{sz_params.n_HC2},
m_n_f66{sz_params.n_f66},
m_n_f61{sz_params.n_f61},
m_n_f61b{sz_params.n_f61b},
m_n_Bf61b{sz_params.n_Bf61b},
m_n_f619{sz_params.n_f619},
m_n_Lc{sz_params.n_Lc},
m_n_Lb{sz_params.n_Lb},
m_n_gamma{sz_params.n_gamma},
m_n_f618{sz_params.n_f618},
m_n_f620{sz_params.n_f620},
m_n_f621{sz_params.n_f621},
m_n_Lh{sz_params.n_Lh},
m_n_f622{sz_params.n_f622},
m_n_f623{sz_params.n_f623},
m_n_f625{sz_params.n_f625},
m_n_f627{sz_params.n_f627},
m_n_Hh{sz_params.n_Hh},
m_n_f72{sz_params.n_f72},
m_n_delnb{sz_params.n_delnb},
m_n_delnc{sz_params.n_delnc},
m_n_delnh{sz_params.n_delnh},
m_n_f74{sz_params.n_f74},
m_n_f75{sz_params.n_f75},
m_n_f714{sz_params.n_f714},
m_n_1f714{sz_params.n_1f714},   
m_n_f739{sz_params.n_f739}, 
m_n_f741{sz_params.n_f741},     
m_n_f743{sz_params.n_f743},     
m_n_f744{sz_params.n_f744},     
m_n_f747{sz_params.n_f747},    
m_n_f96{sz_params.n_f96},    
m_n_f910{sz_params.n_f910},     
m_n_2f96{sz_params.n_2f96},  
m_n_f914{sz_params.n_f914},    
m_n_L{sz_params.n_L},      
m_n_f924{sz_params.n_f924}
{
    if(do_omp_alloc)
    {
        std::cout << "Num-threads= " << omp_get_num_threads() << std::endl;
        std::cout << "Max-threads= " << omp_get_max_threads() << std::endl;

#pragma omp parallel sections
{
        #pragma omp section 
        {
            std::cout << "Executing-thread= " << omp_get_thread_num() << std::endl;
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif
            this->m_data_nidx   = darray_r8_t(this->m_n_nidx);
            this->m_data_n0idx  = darray_r8_t(this->m_n_n0idx);
            this->m_data_z      = darray_r8_t(this->m_n_z);
            this->m_data_z0     = darray_r8_t(this->m_n_z0);
            this->m_data_r      = darray_r8_t(this->m_n_r);
            this->m_data_R0     = darray_r8_t(this->m_n_R0);
            this->m_data_phi    = darray_r8_t(this->m_n_phi);
            this->m_data_phi0   = darray_r8_t(this->m_n_phi0);
            this->m_data_ntht   = darray_r8_t(this->m_n_ntht);
            this->m_data_nphi   = darray_r8_t(this->m_n_nphi);
            this->m_data_dndr   = darray_r8_t(this->m_n_dndr);
            this->m_data_rho    = darray_r8_t(this->m_n_rho);
            this->m_data_dn0    = darray_r8_t(this->m_n_dn0);
            this->m_data_beta   = darray_r8_t(this->m_n_beta);
            this->m_data_navgh  = darray_r8_t(this->m_n_navgh);
            this->m_data_H      = darray_r8_t(this->m_n_H);
            this->m_data_f345   = darray_r8_t(this->m_n_f345);
            this->m_data_f351   = darray_r8_t(this->m_n_f351);
            this->m_data_f352   = darray_r8_t(this->m_n_f352);
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif            
        }

        #pragma omp section 
        {
            std::cout << "Executing-thread= " << omp_get_thread_num() << std::endl;
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif            
            this->m_data_nh    = darray_r8_t(this->m_n_nh);
            this->m_data_f41   = darray_r8_t(this->m_n_f41);
            this->m_data_f     = darray_r8_t(this->m_n_f);
            this->m_data_d     = darray_r8_t(this->m_n_d);
            this->m_data_Nmf   = darray_r8_t(this->m_n_Nmf);
            this->m_data_f412  = darray_r8_t(this->m_n_f412);
            this->m_data_f413  = darray_r8_t(this->m_n_f413);
            this->m_data_D1    = darray_r8_t(this->m_n_D1);
            this->m_data_f415  = darray_r8_t(this->m_n_f415);
            this->m_data_f425  = darray_r8_t(this->m_n_f425);
            this->m_data_f428  = darray_r8_t(this->m_n_f428);
            this->m_data_f429  = darray_r8_t(this->m_n_f429);
            this->m_data_H1    = darray_r8_t(this->m_n_H1);
            this->m_data_H2    = darray_r8_t(this->m_n_H2);
            this->m_data_f430  = darray_r8_t(this->m_n_f430);
            this->m_data_H3    = darray_r8_t(this->m_n_H3);
            this->m_data_deln0 = darray_r8_t(this->m_n_deln0);
            this->m_data_f431  = darray_r8_t(this->m_n_f431);
            this->m_data_f438  = darray_r8_t(this->m_n_f438);
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif 
        }

        #pragma omp section 
        {
            std::cout << "Executing-thread= " << omp_get_thread_num() << std::endl;
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif  
            this->m_data_f442  = darray_r8_t(this->m_n_f442);
            this->m_data_f445  = darray_r8_t(this->m_n_f445);
            this->m_data_g     = darray_r8_t(this->m_n_g);
            this->m_data_f450  = darray_r8_t(this->m_n_f450);
            this->m_data_f451  = darray_r8_t(this->m_n_f451);
            this->m_data_Hc0   = darray_r8_t(this->m_n_Hc0);
            this->m_data_delnA = darray_r8_t(this->m_n_delnA);
            this->m_data_na    = darray_r8_t(this->m_n_na);
            this->m_data_nc    = darray_r8_t(this->m_n_nc);
            this->m_data_f53   = darray_r8_t(this->m_n_f53);
            this->m_data_f517  = darray_r8_t(this->m_n_f517);
            this->m_data_tht   = darray_r8_t(this->m_n_tht);
            this->m_data_f531  = darray_r8_t(this->m_n_f531);
            this->m_data_thtc  = darray_r8_t(this->m_n_thtc);
            this->m_data_f534  = darray_r8_t(this->m_n_f534);
            this->m_data_Rc    = darray_r8_t(this->m_n_Rc);
            this->m_data_f535  = darray_r8_t(this->m_n_f535);
            this->m_data_f538  = darray_r8_t(this->m_n_f538);
            this->m_data_f539  = darray_r8_t(this->m_n_f539);
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif             
        }

        #pragma omp section 
        {
            std::cout << "Executing-thread= " << omp_get_thread_num() << std::endl;
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif  
            this->m_data_f541  = darray_r8_t(this->m_n_f541);
            this->m_data_H0    = darray_r8_t(this->m_n_H0);
            this->m_data_Hc    = darray_r8_t(this->m_n_Hc);
            this->m_data_f543  = darray_r8_t(this->m_n_f543);
            this->m_data_H10   = darray_r8_t(this->m_n_H10);
            this->m_data_f554  = darray_r8_t(this->m_n_f554);
            this->m_data_H20   = darray_r8_t(this->m_n_H20);
            this->m_data_f572  = darray_r8_t(this->m_n_f572);
            this->m_data_f576  = darray_r8_t(this->m_n_f576);
            this->m_data_f577  = darray_r8_t(this->m_n_f577);
            this->m_data_f578  = darray_r8_t(this->m_n_f578);
            this->m_data_f579  = darray_r8_t(this->m_n_f579);
            this->m_data_f590  = darray_r8_t(this->m_n_f590);
            this->m_data_R2    = darray_r8_t(this->m_n_R2);
            this->m_data_f591  = darray_r8_t(this->m_n_f591);
            this->m_data_f593  = darray_r8_t(this->m_n_f593);
            this->m_data_f595  = darray_r8_t(this->m_n_f595);
            this->m_data_HB    = darray_r8_t(this->m_n_HB);
            this->m_data_f62   = darray_r8_t(this->m_n_f62);
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif              
        }

        #pragma omp section 
        {
            std::cout << "Executing-thread= " << omp_get_thread_num() << std::endl;
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif              
            this->m_data_HC2   = darray_r8_t(this->m_n_HC2);
            this->m_data_f66   = darray_r8_t(this->m_n_f66);
            this->m_data_f61   = darray_r8_t(this->m_n_f61);
            this->m_data_f61b  = darray_r8_t(this->m_n_f61b);
            this->m_data_Bf61b = darray_r8_t(this->m_n_Bf61b);
            this->m_data_f619  = darray_r8_t(this->m_n_f619);
            this->m_data_Lc    = darray_r8_t(this->m_n_Lc);
            this->m_data_Lb    = darray_r8_t(this->m_n_Lb);
            this->m_data_gamma = darray_r8_t(this->m_n_gamma);
            this->m_data_f618  = darray_r8_t(this->m_n_f618);
            this->m_data_f620  = darray_r8_t(this->m_n_f620);
            this->m_data_f621  = darray_r8_t(this->m_n_f621);
            this->m_data_Lh    = darray_r8_t(this->m_n_Lh);
            this->m_data_f622  = darray_r8_t(this->m_n_f622);
            this->m_data_f623  = darray_r8_t(this->m_n_f623);
            this->m_data_f625  = darray_r8_t(this->m_n_f625);
            this->m_data_f627  = darray_r8_t(this->m_n_f627);
            this->m_data_Hh    = darray_r8_t(this->m_n_Hh);
            this->m_data_f72   = darray_r8_t(this->m_n_f72);
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif

       }

       #pragma omp section 
       {
              std::cout << "Executing-thread= " << omp_get_thread_num() << std::endl;
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif         
            this->m_data_delnb   = darray_r8_t(this->m_n_delnb);
            this->m_data_delnc   = darray_r8_t(this->m_n_delnc);
            this->m_data_delnh   = darray_r8_t(this->m_n_delnh);
            this->m_data_f74     = darray_r8_t(this->m_n_f74);
            this->m_data_f75     = darray_r8_t(this->m_n_f75);
            this->m_data_f714    = darray_r8_t(this->m_n_f714);
            this->m_data_1f714   = darray_r8_t(this->m_n_1f714);
            this->m_data_f739    = darray_r8_t(this->m_n_f739);
            this->m_data_f741    = darray_r8_t(this->m_n_f741);
            this->m_data_f743    = darray_r8_t(this->m_n_f743);
            this->m_data_f744    = darray_r8_t(this->m_n_f744);
            this->m_data_f747    = darray_r8_t(this->m_n_f747);
            this->m_data_f96     = darray_r8_t(this->m_n_f96);
            this->m_data_f910    = darray_r8_t(this->m_n_f910);
            this->m_data_f2f96   = darray_r8_t(this->m_n_2f96);
            this->m_data_f914    = darray_r8_t(this->m_n_f914);
            this->m_data_L       = darray_r8_t(this->m_n_L);
            this->m_data_f924    = darray_r8_t(this->m_n_f924);
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif            
       }
}

    }
    else 
    {
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif
            this->m_data_nidx   = darray_r8_t(this->m_n_nidx);
            this->m_data_n0idx  = darray_r8_t(this->m_n_n0idx);
            this->m_data_z      = darray_r8_t(this->m_n_z);
            this->m_data_z0     = darray_r8_t(this->m_n_z0);
            this->m_data_r      = darray_r8_t(this->m_n_r);
            this->m_data_R0     = darray_r8_t(this->m_n_R0);
            this->m_data_phi    = darray_r8_t(this->m_n_phi);
            this->m_data_phi0   = darray_r8_t(this->m_n_phi0);
            this->m_data_ntht   = darray_r8_t(this->m_n_ntht);
            this->m_data_nphi   = darray_r8_t(this->m_n_nphi);
            this->m_data_dndr   = darray_r8_t(this->m_n_dndr);
            this->m_data_rho    = darray_r8_t(this->m_n_rho);
            this->m_data_dn0    = darray_r8_t(this->m_n_dn0);
            this->m_data_beta   = darray_r8_t(this->m_n_beta);
            this->m_data_navgh  = darray_r8_t(this->m_n_navgh);
            this->m_data_H      = darray_r8_t(this->m_n_H);
            this->m_data_f345   = darray_r8_t(this->m_n_f345);
            this->m_data_f351   = darray_r8_t(this->m_n_f351);
            this->m_data_f352   = darray_r8_t(this->m_n_f352);

            this->m_data_nh    = darray_r8_t(this->m_n_nh);
            this->m_data_f41   = darray_r8_t(this->m_n_f41);
            this->m_data_f     = darray_r8_t(this->m_n_f);
            this->m_data_d     = darray_r8_t(this->m_n_d);
            this->m_data_Nmf   = darray_r8_t(this->m_n_Nmf);
            this->m_data_f412  = darray_r8_t(this->m_n_f412);
            this->m_data_f413  = darray_r8_t(this->m_n_f413);
            this->m_data_D1    = darray_r8_t(this->m_n_D1);
            this->m_data_f415  = darray_r8_t(this->m_n_f415);
            this->m_data_f425  = darray_r8_t(this->m_n_f425);
            this->m_data_f428  = darray_r8_t(this->m_n_f428);
            this->m_data_f429  = darray_r8_t(this->m_n_f429);
            this->m_data_H1    = darray_r8_t(this->m_n_H1);
            this->m_data_H2    = darray_r8_t(this->m_n_H2);
            this->m_data_f430  = darray_r8_t(this->m_n_f430);
            this->m_data_H3    = darray_r8_t(this->m_n_H3);
            this->m_data_deln0 = darray_r8_t(this->m_n_deln0);
            this->m_data_f431  = darray_r8_t(this->m_n_f431);
            this->m_data_f438  = darray_r8_t(this->m_n_f438);

            this->m_data_f442  = darray_r8_t(this->m_n_f442);
            this->m_data_f445  = darray_r8_t(this->m_n_f445);
            this->m_data_g     = darray_r8_t(this->m_n_g);
            this->m_data_f450  = darray_r8_t(this->m_n_f450);
            this->m_data_f451  = darray_r8_t(this->m_n_f451);
            this->m_data_Hc0   = darray_r8_t(this->m_n_Hc0);
            this->m_data_delnA = darray_r8_t(this->m_n_delnA);
            this->m_data_na    = darray_r8_t(this->m_n_na);
            this->m_data_nc    = darray_r8_t(this->m_n_nc);
            this->m_data_f53   = darray_r8_t(this->m_n_f53);
            this->m_data_f517  = darray_r8_t(this->m_n_f517);
            this->m_data_tht   = darray_r8_t(this->m_n_tht);
            this->m_data_f531  = darray_r8_t(this->m_n_f531);
            this->m_data_thtc  = darray_r8_t(this->m_n_thtc);
            this->m_data_f534  = darray_r8_t(this->m_n_f534);
            this->m_data_Rc    = darray_r8_t(this->m_n_Rc);
            this->m_data_f535  = darray_r8_t(this->m_n_f535);
            this->m_data_f538  = darray_r8_t(this->m_n_f538);
            this->m_data_f539  = darray_r8_t(this->m_n_f539);

            this->m_data_f541  = darray_r8_t(this->m_n_f541);
            this->m_data_H0    = darray_r8_t(this->m_n_H0);
            this->m_data_Hc    = darray_r8_t(this->m_n_Hc);
            this->m_data_f543  = darray_r8_t(this->m_n_f543);
            this->m_data_H10   = darray_r8_t(this->m_n_H10);
            this->m_data_f554  = darray_r8_t(this->m_n_f554);
            this->m_data_H20   = darray_r8_t(this->m_n_H20);
            this->m_data_f572  = darray_r8_t(this->m_n_f572);
            this->m_data_f576  = darray_r8_t(this->m_n_f576);
            this->m_data_f577  = darray_r8_t(this->m_n_f577);
            this->m_data_f578  = darray_r8_t(this->m_n_f578);
            this->m_data_f579  = darray_r8_t(this->m_n_f579);
            this->m_data_f590  = darray_r8_t(this->m_n_f590);
            this->m_data_R2    = darray_r8_t(this->m_n_R2);
            this->m_data_f591  = darray_r8_t(this->m_n_f591);
            this->m_data_f593  = darray_r8_t(this->m_n_f593);
            this->m_data_f595  = darray_r8_t(this->m_n_f595);
            this->m_data_HB    = darray_r8_t(this->m_n_HB);
            this->m_data_f62   = darray_r8_t(this->m_n_f62);

            this->m_data_HC2   = darray_r8_t(this->m_n_HC2);
            this->m_data_f66   = darray_r8_t(this->m_n_f66);
            this->m_data_f61   = darray_r8_t(this->m_n_f61);
            this->m_data_f61b  = darray_r8_t(this->m_n_f61b);
            this->m_data_Bf61b = darray_r8_t(this->m_n_Bf61b);
            this->m_data_f619  = darray_r8_t(this->m_n_f619);
            this->m_data_Lc    = darray_r8_t(this->m_n_Lc);
            this->m_data_Lb    = darray_r8_t(this->m_n_Lb);
            this->m_data_gamma = darray_r8_t(this->m_n_gamma);
            this->m_data_f618  = darray_r8_t(this->m_n_f618);
            this->m_data_f620  = darray_r8_t(this->m_n_f620);
            this->m_data_f621  = darray_r8_t(this->m_n_f621);
            this->m_data_Lh    = darray_r8_t(this->m_n_Lh);
            this->m_data_f622  = darray_r8_t(this->m_n_f622);
            this->m_data_f623  = darray_r8_t(this->m_n_f623);
            this->m_data_f625  = darray_r8_t(this->m_n_f625);
            this->m_data_f627  = darray_r8_t(this->m_n_f627);
            this->m_data_Hh    = darray_r8_t(this->m_n_Hh);
            this->m_data_f72   = darray_r8_t(this->m_n_f72);

            this->m_data_delnb   = darray_r8_t(this->m_n_delnb);
            this->m_data_delnc   = darray_r8_t(this->m_n_delnc);
            this->m_data_delnh   = darray_r8_t(this->m_n_delnh);
            this->m_data_f74     = darray_r8_t(this->m_n_f74);
            this->m_data_f75     = darray_r8_t(this->m_n_f75);
            this->m_data_f714    = darray_r8_t(this->m_n_f714);
            this->m_data_1f714   = darray_r8_t(this->m_n_1f714);
            this->m_data_f739    = darray_r8_t(this->m_n_f739);
            this->m_data_f741    = darray_r8_t(this->m_n_f741);
            this->m_data_f743    = darray_r8_t(this->m_n_f743);
            this->m_data_f744    = darray_r8_t(this->m_n_f744);
            this->m_data_f747    = darray_r8_t(this->m_n_f747);
            this->m_data_f96     = darray_r8_t(this->m_n_f96);
            this->m_data_f910    = darray_r8_t(this->m_n_f910);
            this->m_data_f2f96    = darray_r8_t(this->m_n_2f96);
            this->m_data_f914    = darray_r8_t(this->m_n_f914);
            this->m_data_L       = darray_r8_t(this->m_n_L);
            this->m_data_f924    = darray_r8_t(this->m_n_f924);
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif  
    }

}  

gms::atmos
::atmos_refraction_state_t
::atmos_refraction_state_t(const atmos_refraction_size_params_t &sz_params,
                           const int32_t prot,
                           const int32_t flags,
                           const int32_t fd,
                           const long offset,
                           const int32_t fsize,
                           const bool do_omp_alloc)
:
m_n_nidx{sz_params.n_nidx},
m_n_n0idx{sz_params.n_n0idx},
m_n_z{sz_params.n_z},
m_n_z0{sz_params.n_z0},
m_n_r{sz_params.n_r},
m_n_R0{sz_params.n_R0},
m_n_phi{sz_params.n_phi},
m_n_phi0{sz_params.n_phi0},
m_n_ntht{sz_params.n_ntht},
m_n_nphi{sz_params.n_nphi},
m_n_dndr{sz_params.n_dndr},
m_n_rho{sz_params.n_rho},
m_n_dn0{sz_params.n_dn0},
m_n_beta{sz_params.n_beta},
m_n_navgh{sz_params.n_navgh},
m_n_H{sz_params.n_H},
m_n_f345{sz_params.n_f345},
m_n_f351{sz_params.n_f351},
m_n_f352{sz_params.n_f352},
m_n_nh{sz_params.n_nh},
m_n_f41{sz_params.n_f41},
m_n_f{sz_params.n_f},
m_n_d{sz_params.n_d},
m_n_Nmf{sz_params.n_Nmf},
m_n_f412{sz_params.n_f412},
m_n_f413{sz_params.n_f413},
m_n_D1{sz_params.n_D1},
m_n_f415{sz_params.n_f415},
m_n_f425{sz_params.n_f425},
m_n_f428{sz_params.n_f428},
m_n_f429{sz_params.n_f429},
m_n_H1{sz_params.n_H1},
m_n_H2{sz_params.n_H2},
m_n_f430{sz_params.n_f430},
m_n_H3{sz_params.n_H3},
m_n_deln0{sz_params.n_deln0},
m_n_f431{sz_params.n_f431},
m_n_f438{sz_params.n_f438},
m_n_f442{sz_params.n_f442},
m_n_f445{sz_params.n_f445},
m_n_g{sz_params.n_g},
m_n_f450{sz_params.n_f450},
m_n_f451{sz_params.n_f451},
m_n_Hc0{sz_params.n_Hc0},
m_n_delnA{sz_params.n_delnA},
m_n_na{sz_params.n_na},
m_n_nc{sz_params.n_nc},
m_n_f53{sz_params.n_f53},
m_n_f517{sz_params.n_f517},
m_n_tht{sz_params.n_tht},
m_n_f531{sz_params.n_f531},
m_n_thtc{sz_params.n_thtc},
m_n_f534{sz_params.n_f534},
m_n_Rc{sz_params.n_Rc},
m_n_f535{sz_params.n_f535},
m_n_f538{sz_params.n_f538},
m_n_f539{sz_params.n_f539},
m_n_f541{sz_params.n_f541},
m_n_H0{sz_params.n_H0},
m_n_Hc{sz_params.n_Hc},
m_n_f543{sz_params.n_f543},
m_n_H10{sz_params.n_H10},
m_n_f554{sz_params.n_f554},
m_n_H20{sz_params.n_H20},
m_n_f572{sz_params.n_f572},
m_n_f576{sz_params.n_f576},
m_n_f577{sz_params.n_f577},
m_n_f578{sz_params.n_f578},
m_n_f579{sz_params.n_f579},
m_n_f590{sz_params.n_f590},
m_n_R2{sz_params.n_R2},
m_n_f591{sz_params.n_f591},
m_n_f593{sz_params.n_f593},
m_n_f595{sz_params.n_f595},
m_n_HB{sz_params.n_HB},
m_n_f62{sz_params.n_f62},
m_n_HC2{sz_params.n_HC2},
m_n_f66{sz_params.n_f66},
m_n_f61{sz_params.n_f61},
m_n_f61b{sz_params.n_f61b},
m_n_Bf61b{sz_params.n_Bf61b},
m_n_f619{sz_params.n_f619},
m_n_Lc{sz_params.n_Lc},
m_n_Lb{sz_params.n_Lb},
m_n_gamma{sz_params.n_gamma},
m_n_f618{sz_params.n_f618},
m_n_f620{sz_params.n_f620},
m_n_f621{sz_params.n_f621},
m_n_Lh{sz_params.n_Lh},
m_n_f622{sz_params.n_f622},
m_n_f623{sz_params.n_f623},
m_n_f625{sz_params.n_f625},
m_n_f627{sz_params.n_f627},
m_n_Hh{sz_params.n_Hh},
m_n_f72{sz_params.n_f72},
m_n_delnb{sz_params.n_delnb},
m_n_delnc{sz_params.n_delnc},
m_n_delnh{sz_params.n_delnh},
m_n_f74{sz_params.n_f74},
m_n_f75{sz_params.n_f75},
m_n_f714{sz_params.n_f714},
m_n_1f714{sz_params.n_1f714},   
m_n_f739{sz_params.n_f739}, 
m_n_f741{sz_params.n_f741},     
m_n_f743{sz_params.n_f743},     
m_n_f744{sz_params.n_f744},     
m_n_f747{sz_params.n_f747},    
m_n_f96{sz_params.n_f96},    
m_n_f910{sz_params.n_f910},     
m_n_2f96{sz_params.n_2f96},  
m_n_f914{sz_params.n_f914},    
m_n_L{sz_params.n_L},      
m_n_f924{sz_params.n_f924}
{
    if(do_omp_alloc)
    {
        std::cout << "Num-threads= " << omp_get_num_threads() << std::endl;
        std::cout << "Max-threads= " << omp_get_max_threads() << std::endl;

#pragma omp parallel sections
{
        #pragma omp section 
        {
            std::cout << "Executing-thread= " << omp_get_thread_num() << std::endl;
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif
            this->m_data_nidx   = darray_r8_t(this->m_n_nidx,
                                              prot,flags,fd,offset,fsize);
            this->m_data_n0idx  = darray_r8_t(this->m_n_n0idx,
                                              prot,flags,fd,offset,fsize);
            this->m_data_z      = darray_r8_t(this->m_n_z,
                                              prot,flags,fd,offset,fsize);
            this->m_data_z0     = darray_r8_t(this->m_n_z0,
                                              prot,flags,fd,offset,fsize);
            this->m_data_r      = darray_r8_t(this->m_n_r,
                                              prot,flags,fd,offset,fsize);
            this->m_data_R0     = darray_r8_t(this->m_n_R0,
                                              prot,flags,fd,offset,fsize);
            this->m_data_phi    = darray_r8_t(this->m_n_phi,
                                              prot,flags,fd,offset,fsize);
            this->m_data_phi0   = darray_r8_t(this->m_n_phi0,
                                              prot,flags,fd,offset,fsize);
            this->m_data_ntht   = darray_r8_t(this->m_n_ntht,
                                              prot,flags,fd,offset,fsize);
            this->m_data_nphi   = darray_r8_t(this->m_n_nphi,
                                              prot,flags,fd,offset,fsize);
            this->m_data_dndr   = darray_r8_t(this->m_n_dndr,
                                              prot,flags,fd,offset,fsize);
            this->m_data_rho    = darray_r8_t(this->m_n_rho,
                                              prot,flags,fd,offset,fsize);
            this->m_data_dn0    = darray_r8_t(this->m_n_dn0,
                                              prot,flags,fd,offset,fsize);
            this->m_data_beta   = darray_r8_t(this->m_n_beta,
                                              prot,flags,fd,offset,fsize);
            this->m_data_navgh  = darray_r8_t(this->m_n_navgh,
                                              prot,flags,fd,offset,fsize);
            this->m_data_H      = darray_r8_t(this->m_n_H,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f345   = darray_r8_t(this->m_n_f345,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f351   = darray_r8_t(this->m_n_f351,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f352   = darray_r8_t(this->m_n_f352,
                                              prot,flags,fd,offset,fsize);
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif            
        }

        #pragma omp section 
        {
            std::cout << "Executing-thread= " << omp_get_thread_num() << std::endl;
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif            
            this->m_data_nh    = darray_r8_t(this->m_n_nh,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f41   = darray_r8_t(this->m_n_f41,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f     = darray_r8_t(this->m_n_f,
                                              prot,flags,fd,offset,fsize);
            this->m_data_d     = darray_r8_t(this->m_n_d,
                                              prot,flags,fd,offset,fsize);
            this->m_data_Nmf   = darray_r8_t(this->m_n_Nmf,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f412  = darray_r8_t(this->m_n_f412,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f413  = darray_r8_t(this->m_n_f413,
                                              prot,flags,fd,offset,fsize);
            this->m_data_D1    = darray_r8_t(this->m_n_D1,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f415  = darray_r8_t(this->m_n_f415,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f425  = darray_r8_t(this->m_n_f425,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f428  = darray_r8_t(this->m_n_f428,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f429  = darray_r8_t(this->m_n_f429,
                                              prot,flags,fd,offset,fsize);
            this->m_data_H1    = darray_r8_t(this->m_n_H1,
                                              prot,flags,fd,offset,fsize);
            this->m_data_H2    = darray_r8_t(this->m_n_H2,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f430  = darray_r8_t(this->m_n_f430,
                                              prot,flags,fd,offset,fsize);
            this->m_data_H3    = darray_r8_t(this->m_n_H3,
                                              prot,flags,fd,offset,fsize);
            this->m_data_deln0 = darray_r8_t(this->m_n_deln0,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f431  = darray_r8_t(this->m_n_f431,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f438  = darray_r8_t(this->m_n_f438,
                                              prot,flags,fd,offset,fsize);
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif 
        }

        #pragma omp section 
        {
            std::cout << "Executing-thread= " << omp_get_thread_num() << std::endl;
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif  
            this->m_data_f442  = darray_r8_t(this->m_n_f442,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f445  = darray_r8_t(this->m_n_f445,
                                              prot,flags,fd,offset,fsize);
            this->m_data_g     = darray_r8_t(this->m_n_g,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f450  = darray_r8_t(this->m_n_f450,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f451  = darray_r8_t(this->m_n_f451,
                                              prot,flags,fd,offset,fsize);
            this->m_data_Hc0   = darray_r8_t(this->m_n_Hc0,
                                              prot,flags,fd,offset,fsize);
            this->m_data_delnA = darray_r8_t(this->m_n_delnA,
                                              prot,flags,fd,offset,fsize);
            this->m_data_na    = darray_r8_t(this->m_n_na,
                                              prot,flags,fd,offset,fsize);
            this->m_data_nc    = darray_r8_t(this->m_n_nc,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f53   = darray_r8_t(this->m_n_f53,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f517  = darray_r8_t(this->m_n_f517,
                                              prot,flags,fd,offset,fsize);
            this->m_data_tht   = darray_r8_t(this->m_n_tht,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f531  = darray_r8_t(this->m_n_f531,
                                              prot,flags,fd,offset,fsize);
            this->m_data_thtc  = darray_r8_t(this->m_n_thtc,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f534  = darray_r8_t(this->m_n_f534,
                                              prot,flags,fd,offset,fsize);
            this->m_data_Rc    = darray_r8_t(this->m_n_Rc,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f535  = darray_r8_t(this->m_n_f535,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f538  = darray_r8_t(this->m_n_f538,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f539  = darray_r8_t(this->m_n_f539,
                                              prot,flags,fd,offset,fsize);
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif             
        }

        #pragma omp section 
        {
            std::cout << "Executing-thread= " << omp_get_thread_num() << std::endl;
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif  
            this->m_data_f541  = darray_r8_t(this->m_n_f541,
                                              prot,flags,fd,offset,fsize);
            this->m_data_H0    = darray_r8_t(this->m_n_H0,
                                              prot,flags,fd,offset,fsize);
            this->m_data_Hc    = darray_r8_t(this->m_n_Hc,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f543  = darray_r8_t(this->m_n_f543,
                                              prot,flags,fd,offset,fsize);
            this->m_data_H10   = darray_r8_t(this->m_n_H10,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f554  = darray_r8_t(this->m_n_f554,
                                              prot,flags,fd,offset,fsize);
            this->m_data_H20   = darray_r8_t(this->m_n_H20,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f572  = darray_r8_t(this->m_n_f572,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f576  = darray_r8_t(this->m_n_f576,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f577  = darray_r8_t(this->m_n_f577,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f578  = darray_r8_t(this->m_n_f578,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f579  = darray_r8_t(this->m_n_f579,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f590  = darray_r8_t(this->m_n_f590,
                                              prot,flags,fd,offset,fsize);
            this->m_data_R2    = darray_r8_t(this->m_n_R2,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f591  = darray_r8_t(this->m_n_f591,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f593  = darray_r8_t(this->m_n_f593,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f595  = darray_r8_t(this->m_n_f595,
                                              prot,flags,fd,offset,fsize);
            this->m_data_HB    = darray_r8_t(this->m_n_HB,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f62   = darray_r8_t(this->m_n_f62,
                                              prot,flags,fd,offset,fsize);
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif              
        }

        #pragma omp section 
        {
            std::cout << "Executing-thread= " << omp_get_thread_num() << std::endl;
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif              
            this->m_data_HC2   = darray_r8_t(this->m_n_HC2,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f66   = darray_r8_t(this->m_n_f66,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f61   = darray_r8_t(this->m_n_f61,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f61b  = darray_r8_t(this->m_n_f61b,
                                              prot,flags,fd,offset,fsize);
            this->m_data_Bf61b = darray_r8_t(this->m_n_Bf61b,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f619  = darray_r8_t(this->m_n_f619,
                                              prot,flags,fd,offset,fsize);
            this->m_data_Lc    = darray_r8_t(this->m_n_Lc,
                                              prot,flags,fd,offset,fsize);
            this->m_data_Lb    = darray_r8_t(this->m_n_Lb,
                                              prot,flags,fd,offset,fsize);
            this->m_data_gamma = darray_r8_t(this->m_n_gamma,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f618  = darray_r8_t(this->m_n_f618,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f620  = darray_r8_t(this->m_n_f620,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f621  = darray_r8_t(this->m_n_f621,
                                              prot,flags,fd,offset,fsize);
            this->m_data_Lh    = darray_r8_t(this->m_n_Lh,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f622  = darray_r8_t(this->m_n_f622,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f623  = darray_r8_t(this->m_n_f623,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f625  = darray_r8_t(this->m_n_f625,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f627  = darray_r8_t(this->m_n_f627,
                                              prot,flags,fd,offset,fsize);
            this->m_data_Hh    = darray_r8_t(this->m_n_Hh,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f72   = darray_r8_t(this->m_n_f72,
                                              prot,flags,fd,offset,fsize);
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif

       }

       #pragma omp section 
       {
              std::cout << "Executing-thread= " << omp_get_thread_num() << std::endl;
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif         
            this->m_data_delnb   = darray_r8_t(this->m_n_delnb,
                                              prot,flags,fd,offset,fsize);
            this->m_data_delnc   = darray_r8_t(this->m_n_delnc,
                                              prot,flags,fd,offset,fsize);
            this->m_data_delnh   = darray_r8_t(this->m_n_delnh,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f74     = darray_r8_t(this->m_n_f74,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f75     = darray_r8_t(this->m_n_f75,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f714    = darray_r8_t(this->m_n_f714,
                                              prot,flags,fd,offset,fsize);
            this->m_data_1f714   = darray_r8_t(this->m_n_1f714,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f739    = darray_r8_t(this->m_n_f739,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f741    = darray_r8_t(this->m_n_f741,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f743    = darray_r8_t(this->m_n_f743,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f744    = darray_r8_t(this->m_n_f744,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f747    = darray_r8_t(this->m_n_f747,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f96     = darray_r8_t(this->m_n_f96,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f910    = darray_r8_t(this->m_n_f910,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f2f96   = darray_r8_t(this->m_n_2f96,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f914    = darray_r8_t(this->m_n_f914,
                                              prot,flags,fd,offset,fsize);
            this->m_data_L       = darray_r8_t(this->m_n_L,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f924    = darray_r8_t(this->m_n_f924,
                                              prot,flags,fd,offset,fsize);
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif            
       }
}

    }
    else 
    {
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              PMC_VARS
              HW_PMC_COLLECTION_PROLOGE_BODY
#endif
            this->m_data_nidx   = darray_r8_t(this->m_n_nidx,
                                              prot,flags,fd,offset,fsize);
            this->m_data_n0idx  = darray_r8_t(this->m_n_n0idx,
                                              prot,flags,fd,offset,fsize);
            this->m_data_z      = darray_r8_t(this->m_n_z,
                                              prot,flags,fd,offset,fsize);
            this->m_data_z0     = darray_r8_t(this->m_n_z0,
                                              prot,flags,fd,offset,fsize);
            this->m_data_r      = darray_r8_t(this->m_n_r,
                                              prot,flags,fd,offset,fsize);
            this->m_data_R0     = darray_r8_t(this->m_n_R0,
                                              prot,flags,fd,offset,fsize);
            this->m_data_phi    = darray_r8_t(this->m_n_phi,
                                              prot,flags,fd,offset,fsize);
            this->m_data_phi0   = darray_r8_t(this->m_n_phi0,
                                              prot,flags,fd,offset,fsize);
            this->m_data_ntht   = darray_r8_t(this->m_n_ntht,
                                              prot,flags,fd,offset,fsize);
            this->m_data_nphi   = darray_r8_t(this->m_n_nphi,
                                              prot,flags,fd,offset,fsize);
            this->m_data_dndr   = darray_r8_t(this->m_n_dndr,
                                              prot,flags,fd,offset,fsize);
            this->m_data_rho    = darray_r8_t(this->m_n_rho,
                                              prot,flags,fd,offset,fsize);
            this->m_data_dn0    = darray_r8_t(this->m_n_dn0,
                                              prot,flags,fd,offset,fsize);
            this->m_data_beta   = darray_r8_t(this->m_n_beta,
                                              prot,flags,fd,offset,fsize);
            this->m_data_navgh  = darray_r8_t(this->m_n_navgh,
                                              prot,flags,fd,offset,fsize);
            this->m_data_H      = darray_r8_t(this->m_n_H,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f345   = darray_r8_t(this->m_n_f345,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f351   = darray_r8_t(this->m_n_f351,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f352   = darray_r8_t(this->m_n_f352,
                                              prot,flags,fd,offset,fsize);

            this->m_data_nh    = darray_r8_t(this->m_n_nh,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f41   = darray_r8_t(this->m_n_f41,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f     = darray_r8_t(this->m_n_f,
                                              prot,flags,fd,offset,fsize);
            this->m_data_d     = darray_r8_t(this->m_n_d,
                                              prot,flags,fd,offset,fsize);
            this->m_data_Nmf   = darray_r8_t(this->m_n_Nmf,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f412  = darray_r8_t(this->m_n_f412,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f413  = darray_r8_t(this->m_n_f413,
                                              prot,flags,fd,offset,fsize);
            this->m_data_D1    = darray_r8_t(this->m_n_D1,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f415  = darray_r8_t(this->m_n_f415,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f425  = darray_r8_t(this->m_n_f425,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f428  = darray_r8_t(this->m_n_f428,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f429  = darray_r8_t(this->m_n_f429,
                                              prot,flags,fd,offset,fsize);
            this->m_data_H1    = darray_r8_t(this->m_n_H1,
                                              prot,flags,fd,offset,fsize);
            this->m_data_H2    = darray_r8_t(this->m_n_H2,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f430  = darray_r8_t(this->m_n_f430,
                                              prot,flags,fd,offset,fsize);
            this->m_data_H3    = darray_r8_t(this->m_n_H3,
                                              prot,flags,fd,offset,fsize);
            this->m_data_deln0 = darray_r8_t(this->m_n_deln0,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f431  = darray_r8_t(this->m_n_f431,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f438  = darray_r8_t(this->m_n_f438,
                                              prot,flags,fd,offset,fsize);

            this->m_data_f442  = darray_r8_t(this->m_n_f442,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f445  = darray_r8_t(this->m_n_f445,
                                              prot,flags,fd,offset,fsize);
            this->m_data_g     = darray_r8_t(this->m_n_g,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f450  = darray_r8_t(this->m_n_f450,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f451  = darray_r8_t(this->m_n_f451,
                                              prot,flags,fd,offset,fsize);
            this->m_data_Hc0   = darray_r8_t(this->m_n_Hc0,
                                              prot,flags,fd,offset,fsize);
            this->m_data_delnA = darray_r8_t(this->m_n_delnA,
                                              prot,flags,fd,offset,fsize);
            this->m_data_na    = darray_r8_t(this->m_n_na,
                                              prot,flags,fd,offset,fsize);
            this->m_data_nc    = darray_r8_t(this->m_n_nc,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f53   = darray_r8_t(this->m_n_f53,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f517  = darray_r8_t(this->m_n_f517,
                                              prot,flags,fd,offset,fsize);
            this->m_data_tht   = darray_r8_t(this->m_n_tht,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f531  = darray_r8_t(this->m_n_f531,
                                              prot,flags,fd,offset,fsize);
            this->m_data_thtc  = darray_r8_t(this->m_n_thtc,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f534  = darray_r8_t(this->m_n_f534,
                                              prot,flags,fd,offset,fsize);
            this->m_data_Rc    = darray_r8_t(this->m_n_Rc,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f535  = darray_r8_t(this->m_n_f535,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f538  = darray_r8_t(this->m_n_f538,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f539  = darray_r8_t(this->m_n_f539,
                                              prot,flags,fd,offset,fsize);

            this->m_data_f541  = darray_r8_t(this->m_n_f541,
                                              prot,flags,fd,offset,fsize);
            this->m_data_H0    = darray_r8_t(this->m_n_H0,
                                              prot,flags,fd,offset,fsize);
            this->m_data_Hc    = darray_r8_t(this->m_n_Hc,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f543  = darray_r8_t(this->m_n_f543,
                                              prot,flags,fd,offset,fsize);
            this->m_data_H10   = darray_r8_t(this->m_n_H10,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f554  = darray_r8_t(this->m_n_f554,
                                              prot,flags,fd,offset,fsize);
            this->m_data_H20   = darray_r8_t(this->m_n_H20,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f572  = darray_r8_t(this->m_n_f572,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f576  = darray_r8_t(this->m_n_f576,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f577  = darray_r8_t(this->m_n_f577,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f578  = darray_r8_t(this->m_n_f578,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f579  = darray_r8_t(this->m_n_f579,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f590  = darray_r8_t(this->m_n_f590,
                                              prot,flags,fd,offset,fsize);
            this->m_data_R2    = darray_r8_t(this->m_n_R2,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f591  = darray_r8_t(this->m_n_f591,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f593  = darray_r8_t(this->m_n_f593,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f595  = darray_r8_t(this->m_n_f595,
                                              prot,flags,fd,offset,fsize);
            this->m_data_HB    = darray_r8_t(this->m_n_HB,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f62   = darray_r8_t(this->m_n_f62,
                                              prot,flags,fd,offset,fsize);

            this->m_data_HC2   = darray_r8_t(this->m_n_HC2,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f66   = darray_r8_t(this->m_n_f66,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f61   = darray_r8_t(this->m_n_f61,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f61b  = darray_r8_t(this->m_n_f61b,
                                              prot,flags,fd,offset,fsize);
            this->m_data_Bf61b = darray_r8_t(this->m_n_Bf61b,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f619  = darray_r8_t(this->m_n_f619,
                                              prot,flags,fd,offset,fsize);
            this->m_data_Lc    = darray_r8_t(this->m_n_Lc,
                                              prot,flags,fd,offset,fsize);
            this->m_data_Lb    = darray_r8_t(this->m_n_Lb,
                                              prot,flags,fd,offset,fsize);
            this->m_data_gamma = darray_r8_t(this->m_n_gamma,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f618  = darray_r8_t(this->m_n_f618,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f620  = darray_r8_t(this->m_n_f620,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f621  = darray_r8_t(this->m_n_f621,
                                              prot,flags,fd,offset,fsize);
            this->m_data_Lh    = darray_r8_t(this->m_n_Lh,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f622  = darray_r8_t(this->m_n_f622,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f623  = darray_r8_t(this->m_n_f623,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f625  = darray_r8_t(this->m_n_f625,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f627  = darray_r8_t(this->m_n_f627,
                                              prot,flags,fd,offset,fsize);
            this->m_data_Hh    = darray_r8_t(this->m_n_Hh,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f72   = darray_r8_t(this->m_n_f72,
                                              prot,flags,fd,offset,fsize);

            this->m_data_delnb   = darray_r8_t(this->m_n_delnb,
                                              prot,flags,fd,offset,fsize);
            this->m_data_delnc   = darray_r8_t(this->m_n_delnc,
                                              prot,flags,fd,offset,fsize);
            this->m_data_delnh   = darray_r8_t(this->m_n_delnh,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f74     = darray_r8_t(this->m_n_f74,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f75     = darray_r8_t(this->m_n_f75,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f714    = darray_r8_t(this->m_n_f714,
                                              prot,flags,fd,offset,fsize);
            this->m_data_1f714   = darray_r8_t(this->m_n_1f714,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f739    = darray_r8_t(this->m_n_f739,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f741    = darray_r8_t(this->m_n_f741,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f743    = darray_r8_t(this->m_n_f743,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f744    = darray_r8_t(this->m_n_f744,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f747    = darray_r8_t(this->m_n_f747,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f96     = darray_r8_t(this->m_n_f96,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f910    = darray_r8_t(this->m_n_f910,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f2f96   = darray_r8_t(this->m_n_2f96,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f914    = darray_r8_t(this->m_n_f914,
                                              prot,flags,fd,offset,fsize);
            this->m_data_L       = darray_r8_t(this->m_n_L,
                                              prot,flags,fd,offset,fsize);
            this->m_data_f924    = darray_r8_t(this->m_n_f924,
                                              prot,flags,fd,offset,fsize);
#if (ATMOS_REFRACTION_STATE_USE_PMC_INSTRUMENTATION) == 1
              HW_PMC_COLLECTION_EPILOGE_BODY

              HW_PMC_SHOW_INFO
	            
#endif  
    }

}  

gms::atmos
::atmos_refraction_state_t
::~atmos_refraction_state_t() noexcept(false)
{

}



