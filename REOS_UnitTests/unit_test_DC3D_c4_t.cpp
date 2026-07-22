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
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_containers.h"

/*
   icpc -o unit_test_DC3D_c4_t -fp-model fast=2 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_dyn_containers.h unit_test_DC3D_c4_t.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -march=skylake-avx512 -mavx512f -falign-functions=32   GMS_config.h GMS_malloc.h GMS_dyn_containers.h unit_test_DC3D_c4_t.cpp

*/



void unit_test_DC3D_c4_t_2nd_Ctor();

void unit_test_DC3D_c4_t_2nd_Ctor()
{
       using namespace gms;
       constexpr std::size_t nx{4096};
       constexpr std::size_t ny{1024};
       constexpr std::size_t nz{512};
       int32_t status;
       const char * fname{"UNIT_TEST_Output_Random_Data_DC3D_c4_t_2nd_Ctor.csv"};
       std::clock_t seeds[6] = {0ULL};
       std::uniform_real_distribution<float> distros[6]{};
       FILE * fp{NULL};
       DC3D_c4_t test_2nd_Ctor = DC3D_c4_t(nx,ny,nz);
       char * ctor_name{gms::common::demangle(typeid(test_2nd_Ctor).name(),status)};
       if(status==0 && ctor_name != NULL)
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
       }
       else
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(test_2nd_Ctor).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
       }
       test_2nd_Ctor.info_size_alignment();
       fp = fopen(fname,"w+");
       if(!fp)
       {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           std::exit(EXIT_FAILURE);
       }
       printf("[UNIT-TEST]: -- Start of random data generation and initialization.\n");
       for(int __i{0}; __i != 6; ++__i)
       {
           seeds[__i] = std::clock();
           distros[__i]  = std::uniform_real_distribution<float>(0.0f,1.0f);
       }
       auto rdev1{std::mt19937(seeds[0])};
       auto rdev2{std::mt19937(seeds[1])};
       for(std::size_t __i{0ULL}; __i != test_2nd_Ctor.mnx; ++__i)
       {
           const std::complex<float> rc{distros[0].operator()(rdev1),distros[1].operator()(rdev2)};
           test_2nd_Ctor.m_Ex[__i] = rc;
       }
       auto rdev3{std::mt19937(seeds[2])};
       auto rdev4{std::mt19937(seeds[3])};
       for(std::size_t __i{0ULL}; __i != test_2nd_Ctor.mny; ++__i)
       {
           const std::complex<float> rc{distros[2].operator()(rdev3),distros[3].operator()(rdev4)};
           test_2nd_Ctor.m_Ey[__i] = rc;
       }
       auto rdev5{std::mt19937(seeds[4])};
       auto rdev6{std::mt19937(seeds[5])};
       for(std::size_t __i{0ULL}; __i != test_2nd_Ctor.mnz; ++__i)
       {
           const std::complex<float> rc{distros[4].operator()(rdev5),distros[5].operator()(rdev6)};
           test_2nd_Ctor.m_Ez[__i] = rc;
       }
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       printf("[UNIT-TEST]: Dumping: m_Ex \n");
       for(std::size_t __i{0}; __i != test_2nd_Ctor.mnx; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", test_2nd_Ctor.m_Ex[__i].real(),
                                        test_2nd_Ctor.m_Ex[__i].imag());
       }
      printf("[UNIT-TEST]: Dumping: m_Ey \n");
      for(std::size_t __i{0}; __i != test_2nd_Ctor.mny; ++__i)
      {
           fprintf(fp,"[%.7f,j%.7f]\n", test_2nd_Ctor.m_Ey[__i].real(),
                                        test_2nd_Ctor.m_Ey[__i].imag());
      }
      printf("[UNIT-TEST]: Dumping: m_Ez \n");
      for(std::size_t __i{0}; __i != test_2nd_Ctor.mnz; ++__i)
      {
           fprintf(fp,"[%.7f,j%.7f]\n", test_2nd_Ctor.m_Ez[__i].real(),
                                        test_2nd_Ctor.m_Ez[__i].imag());
      }
      fclose(fp);
}

void unit_test_DC3D_c4_t_3rd_Ctor();

void unit_test_DC3D_c4_t_3rd_Ctor()
{
       using namespace gms;
       constexpr std::size_t nx{4*4096};
       constexpr std::size_t ny{3*4096};
       constexpr std::size_t nz{2*4096};
       constexpr int32_t prot{PROT_READ | PROT_WRITE};
       constexpr int32_t flags{MAP_ANONYMOUS | MAP_PRIVATE};
       constexpr int32_t fd{-1};
       constexpr off_t offset{0};
       constexpr int32_t fsize{0};
       const char * fname{"UNIT_TEST_Output_Random_Data_DC3D_c4_t_3rd_Ctor.csv"};
       std::clock_t seeds[6] = {0ULL};
       std::uniform_real_distribution<float> distros[6]{};
       FILE * fp{NULL};
       int32_t status;
       DC3D_c4_t test_3rd_Ctor = DC3D_c4_t(nx,ny,nz,prot,flags,fd,offset,fsize);
       char * ctor_name{gms::common::demangle(typeid(test_3rd_Ctor).name(),status)};
       if(status==0 && ctor_name != NULL)
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
       }
       else
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(test_3rd_Ctor).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
       }
       test_3rd_Ctor.info_size_alignment();
       fp = fopen(fname,"w+");
       if(!fp)
       {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           std::exit(EXIT_FAILURE);
       }
       printf("[UNIT-TEST]: -- Start of random data generation and initialization.\n");                                
       for(int32_t __i{0}; __i != 6; ++__i)
       {
           seeds[__i]   = std::clock();
           distros[__i] = std::uniform_real_distribution<float>(0.0f,1.0f); 
       }  
       auto rdev1{std::mt19937(seeds[0])};
       auto rdev2{std::mt19937(seeds[1])};
       for(std::size_t __i{0ULL}; __i != test_3rd_Ctor.mnx; ++__i)
       {
           test_3rd_Ctor.m_Ex[__i] = std::complex<float>{distros[0].operator()(rdev1),
                                                         distros[1].operator()(rdev2)};
       }
       auto rdev3{std::mt19937(seeds[2])};
       auto rdev4{std::mt19937(seeds[3])};
       for(std::size_t __i{0ULL}; __i != test_3rd_Ctor.mny; ++__i)
       {
           test_3rd_Ctor.m_Ey[__i] = std::complex<float>{distros[2].operator()(rdev3),
                                                         distros[3].operator()(rdev4)};
       }
       auto rdev5{std::mt19937(seeds[4])};
       auto rdev6{std::mt19937(seeds[5])};
       for(std::size_t __i{0ULL}; __i != test_3rd_Ctor.mnz; ++__i)
       {
           test_3rd_Ctor.m_Ez[__i] = std::complex<float>{distros[4].operator()(rdev5),
                                                         distros[5].operator()(rdev6)};
       }
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       printf("[UNIT-TEST]: Dumping: m_Ex \n");
       for(std::size_t __i{0}; __i != test_3rd_Ctor.mnx; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", test_3rd_Ctor.m_Ex[__i].real(),
                                        test_3rd_Ctor.m_Ex[__i].imag());
       }
       printf("[UNIT-TEST]: Dumping: m_Ey \n");
       for(std::size_t __i{0}; __i != test_3rd_Ctor.mny; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", test_3rd_Ctor.m_Ey[__i].real(),
                                        test_3rd_Ctor.m_Ey[__i].imag());
       }
       printf("[UNIT-TEST]: Dumping: m_Ez \n");
       for(std::size_t __i{0}; __i != test_3rd_Ctor.mnz; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", test_3rd_Ctor.m_Ez[__i].real(),
                                        test_3rd_Ctor.m_Ez[__i].imag());
       }
       fclose(fp);
}


void unit_test_DC3D_c4_t_4th_Ctor();

void unit_test_DC3D_c4_t_4th_Ctor()
{
       using namespace gms;
       constexpr std::size_t nx{4096};
       constexpr std::size_t ny{1024};
       constexpr std::size_t nz{512};
       int32_t status;
       const char * fname{"UNIT_TEST_Output_Random_Data_DC3D_c4_t_4th_Ctor.csv"};
       std::clock_t seeds[6] = {0ULL};
       std::uniform_real_distribution<float> distros[6]{};
       FILE * fp{NULL};
       std::vector<std::complex<float>> Ex = std::vector<std::complex<float>>(nx);
       std::vector<std::complex<float>> Ey = std::vector<std::complex<float>>(ny);
       std::vector<std::complex<float>> Ez = std::vector<std::complex<float>>(nz);
       printf("[UNIT-TEST]: -- Start of random data generation and initialization.\n");                                
       for(int32_t __i{0}; __i != 6; ++__i)
       {
           seeds[__i]   = std::clock();
           distros[__i] = std::uniform_real_distribution<float>(0.0f,1.0f); 
       }  
       auto rdev1{std::mt19937(seeds[0])};
       auto rdev2{std::mt19937(seeds[1])};
       for(std::size_t __i{0ULL}; __i != nx; ++__i)
       {
           Ex[__i] = std::complex<float>{distros[0].operator()(rdev1),
                                         distros[1].operator()(rdev2)};
       }
       auto rdev3{std::mt19937(seeds[2])};
       auto rdev4{std::mt19937(seeds[3])};
       for(std::size_t __i{0ULL}; __i != ny; ++__i)
       {
           Ey[__i] = std::complex<float>{distros[2].operator()(rdev3),
                                         distros[3].operator()(rdev4)};
       }
       auto rdev5{std::mt19937(seeds[4])};
       auto rdev6{std::mt19937(seeds[5])};
       for(std::size_t __i{0ULL}; __i != nz; ++__i)
       {
           Ez[__i] = std::complex<float>{distros[4].operator()(rdev5),
                                         distros[5].operator()(rdev6)};
       }
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       DC3D_c4_t test_4th_Ctor = DC3D_c4_t(Ex,Ey,Ez);
       char * ctor_name{gms::common::demangle(typeid(test_4th_Ctor).name(),status)};
       if(status==0 && ctor_name != NULL)
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
       }
       else
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(test_4th_Ctor).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
       }
       test_4th_Ctor.info_size_alignment();
       fp = fopen(fname,"w+");
       if(!fp)
       {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           std::exit(EXIT_FAILURE);
       }
       printf("[UNIT-TEST]: Dumping: m_Ex \n");
       for(std::size_t __i{0}; __i != test_4th_Ctor.mnx; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", test_4th_Ctor.m_Ex[__i].real(),
                                        test_4th_Ctor.m_Ex[__i].imag());
       }
       printf("[UNIT-TEST]: Dumping: m_Ey \n");
       for(std::size_t __i{0}; __i != test_4th_Ctor.mny; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", test_4th_Ctor.m_Ey[__i].real(),
                                        test_4th_Ctor.m_Ey[__i].imag());
       }
       printf("[UNIT-TEST]: Dumping: m_Ez \n");
       for(std::size_t __i{0}; __i != test_4th_Ctor.mnz; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", test_4th_Ctor.m_Ez[__i].real(),
                                        test_4th_Ctor.m_Ez[__i].imag());
       }
       fclose(fp);

}

void unit_test_DC3D_c4_t_5th_Ctor();

void unit_test_DC3D_c4_t_5th_Ctor()
{
       using namespace gms;
       constexpr std::size_t nx{4096};
       constexpr std::size_t ny{1024};
       constexpr std::size_t nz{512};
       int32_t status;
       const char * fname{"UNIT_TEST_Output_Random_Data_DC3D_c4_t_5th_Ctor.csv"};
       std::clock_t seeds[6] = {0ULL};
       std::uniform_real_distribution<float> distros[6]{};
       FILE * fp{NULL};
       std::valarray<std::complex<float>> Ex = std::valarray<std::complex<float>>(nx);
       std::valarray<std::complex<float>> Ey = std::valarray<std::complex<float>>(ny);
       std::valarray<std::complex<float>> Ez = std::valarray<std::complex<float>>(nz);
       printf("[UNIT-TEST]: -- Start of random data generation and initialization.\n");                                
       for(int32_t __i{0}; __i != 6; ++__i)
       {
           seeds[__i]   = std::clock();
           distros[__i] = std::uniform_real_distribution<float>(0.0f,1.0f); 
       }  
       auto rdev1{std::mt19937(seeds[0])};
       auto rdev2{std::mt19937(seeds[1])};
       for(std::size_t __i{0ULL}; __i != nx; ++__i)
       {
           Ex[__i] = std::complex<float>{distros[0].operator()(rdev1),
                                         distros[1].operator()(rdev2)};
       }
       auto rdev3{std::mt19937(seeds[2])};
       auto rdev4{std::mt19937(seeds[3])};
       for(std::size_t __i{0ULL}; __i != ny; ++__i)
       {
           Ey[__i] = std::complex<float>{distros[2].operator()(rdev3),
                                         distros[3].operator()(rdev4)};
       }
       auto rdev5{std::mt19937(seeds[4])};
       auto rdev6{std::mt19937(seeds[5])};
       for(std::size_t __i{0ULL}; __i != nz; ++__i)
       {
           Ez[__i] = std::complex<float>{distros[4].operator()(rdev5),
                                         distros[5].operator()(rdev6)};
       }
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       DC3D_c4_t test_5th_Ctor = DC3D_c4_t(Ex,Ey,Ez);
       char * ctor_name{gms::common::demangle(typeid(test_5th_Ctor).name(),status)};
       if(status==0 && ctor_name != NULL)
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
       }
       else
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(test_5th_Ctor).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
       }
       test_5th_Ctor.info_size_alignment();
       fp = fopen(fname,"w+");
       if(!fp)
       {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           std::exit(EXIT_FAILURE);
       }
       printf("[UNIT-TEST]: Dumping: m_Ex \n");
       for(std::size_t __i{0}; __i != test_5th_Ctor.mnx; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", test_5th_Ctor.m_Ex[__i].real(),
                                        test_5th_Ctor.m_Ex[__i].imag());
       }
       printf("[UNIT-TEST]: Dumping: m_Ey \n");
       for(std::size_t __i{0}; __i != test_5th_Ctor.mny; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", test_5th_Ctor.m_Ey[__i].real(),
                                        test_5th_Ctor.m_Ey[__i].imag());
       }
       printf("[UNIT-TEST]: Dumping: m_Ez \n");
       for(std::size_t __i{0}; __i != test_5th_Ctor.mnz; ++__i)
       {
           fprintf(fp,"[%.7f,j%.7f]\n", test_5th_Ctor.m_Ez[__i].real(),
                                        test_5th_Ctor.m_Ez[__i].imag());
       }
       fclose(fp);

}





int main()
{   
    unit_test_DC3D_c4_t_2nd_Ctor();
    unit_test_DC3D_c4_t_3rd_Ctor();
    unit_test_DC3D_c4_t_4th_Ctor();
    unit_test_DC3D_c4_t_5th_Ctor();
    return 0;
}
