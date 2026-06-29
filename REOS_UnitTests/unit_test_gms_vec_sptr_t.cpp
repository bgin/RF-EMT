#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_malloc.h"
#include "GMS_dyn_vec_unique_ptr.h"

/*
   icpc -o unit_test_gms_vec_sptr_t -fp-model fast=2 -ftz -ggdb -ipo  --std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_dyn_vec_unique_ptr.h unit_test_gms_vec_sptr_t.cpp 
   ASM: 
   icpc -S -fverbose-asm -masm=intel  --std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32  GMS_config.h GMS_malloc.h GMS_dyn_vec_unique_ptr.h unit_test_gms_vec_sptr_t.cpp

*/


void unit_test_1_obj_construct();

void unit_test_1_obj_construct()
{
       using namespace gms::common;
       constexpr std::size_t nx{18ull*4096ull};
       int32_t status;
       const char * fname{"UNIT_TEST_Output_Random_Data_test_1_obj_construct.csv"};
       std::clock_t seed{0ULL};
       std::uniform_real_distribution<float> distro{};
       FILE * fp{NULL};
       //darray_r4_t test_2nd_Ctor = darray_r4_t(nx);
       printf("[UNIT-TEST #1]: test_1_obj_construct -- ***START***!!\n");
       gms_vec_sptr_t<float> ctor_1_arg(nx);
       char * ctor_name{gms::common::demangle(typeid(ctor_1_arg).name(),status)};
       if(status==0 && ctor_name != NULL)
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
       }
       else
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(ctor_1_arg).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
       }
       fp = fopen(fname,"w+");
       if(!fp)
       {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           std::exit(EXIT_FAILURE);
       }
       printf("[UNIT-TEST]: -- Start of random data generation and initialization.\n");
       seed = std::clock();
       distro  = std::uniform_real_distribution<float>(0.0f,1.0f);
       auto rdev1{std::mt19937(seed)};
       for(std::size_t __i{0ULL}; __i != ctor_1_arg.size(); ++__i)
       {
           const float rc{distro.operator()(rdev1)};
           ctor_1_arg.operator[](__i) = rc;
       }
       
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       printf("[UNIT-TEST]: Dumping: m_data to file=%s\n",fname);
       for(std::size_t __i{0ull}; __i != ctor_1_arg.size(); ++__i)
       {
           fprintf(fp,"[%.7f]\n", ctor_1_arg.operator[](__i));
                   
       }
#if 0
       std::cout << ctor_1_arg << "\n";
#endif     
       fclose(fp);
       printf("[UNIT-TEST #1]: test_1_obj_construct -- ***END***!!\n");
}

void unit_test_2_obj_init_construct();

void unit_test_2_obj_init_construct()
{
       using namespace gms::common;
       constexpr std::size_t nx{4096ull};
       int32_t status;
       const char * fname{"UNIT_TEST_Output_Random_Data_test_2_obj_init_construct.csv"};
       std::clock_t seed{0ULL};
       std::uniform_real_distribution<float> distro{};
       FILE * fp{NULL};
       //darray_r4_t test_2nd_Ctor = darray_r4_t(nx);
       gms_vec_sptr_t<float> ctor_2_arg;
       printf("[UNIT-TEST #2]: test_2_obj_init_construct -- ***START***!!\n");
       char * ctor_name{gms::common::demangle(typeid(ctor_2_arg).name(),status)};
       if(status==0 && ctor_name != NULL)
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
       }
       else
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(ctor_2_arg).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
       }
       fp = fopen(fname,"w+");
       if(!fp)
       {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           std::exit(EXIT_FAILURE);
       }
       printf("[UNIT-TEST]: -- Start of random data generation and initialization.\n");
       seed    = std::clock();
       distro  = std::uniform_real_distribution<float>(0.0f,1.0f);
       auto rdev1{std::mt19937(seed)};
       ctor_2_arg = gms_vec_sptr_t<float>(nx,distro.operator()(rdev1));       
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       printf("[UNIT-TEST]: Dumping: m_data to file=%s\n",fname);
       for(std::size_t __i{0ull}; __i != ctor_2_arg.size(); ++__i)
       {
           fprintf(fp,"[%.7f]\n", ctor_2_arg.operator[](__i));
                   
       }
#if 0
       std::cout << ctor_2_arg << "\n";
#endif     
       fclose(fp);
       printf("[UNIT-TEST #2]: test_2_obj_init_construct -- ***END***!!\n");
}

void unit_test_3_obj_construct();

void unit_test_3_obj_construct()
{
       using namespace gms::common;
       constexpr std::size_t nx{4096ull};
       int32_t status;
       const char * fname{"UNIT_TEST_Output_Random_Data_test_3_obj_construct.csv"};
       std::clock_t seed{0ULL};
       std::uniform_real_distribution<float> distro{};
       float * __restrict ptr_rand_buf = NULL;
       FILE * fp{NULL};
       //darray_r4_t test_2nd_Ctor = darray_r4_t(nx);
       printf("[UNIT-TEST #1]: test_3_obj_construct -- ***START***!!\n");
       ptr_rand_buf = reinterpret_cast<float * __restrict>(gms_mm_malloc(nx*sizeof(float),64));
       gms_vec_sptr_t<float> ctor_2_arg;
       char * ctor_name{gms::common::demangle(typeid(ctor_2_arg).name(),status)};
       if(status==0 && ctor_name != NULL)
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
       }
       else
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(ctor_2_arg).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
       }
       fp = fopen(fname,"w+");
       if(!fp)
       {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           std::exit(EXIT_FAILURE);
       }
       printf("[UNIT-TEST]: -- Start of random data generation and initialization.\n");
       seed = std::clock();
       distro  = std::uniform_real_distribution<float>(0.0f,1.0f);
       auto rdev1{std::mt19937(seed)};
       for(std::size_t __i{0ULL}; __i != nx; ++__i)
       {
           const float rc{distro.operator()(rdev1)};
           ptr_rand_buf[__i] = rc;
       }
       ctor_2_arg = gms_vec_sptr_t<float>(ptr_rand_buf,nx);
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       printf("[UNIT-TEST]: Dumping: m_data to file=%s\n",fname);
       for(std::size_t __i{0ull}; __i != ctor_2_arg.size(); ++__i)
       {
           fprintf(fp,"[%.7f]\n", ctor_2_arg.operator[](__i));
                   
       }
#if 0
       std::cout << ctor_1_arg << "\n";
#endif     
       fclose(fp);
       if(ptr_rand_buf && nx!=0ull) { gms_mm_free(ptr_rand_buf);}
       printf("[UNIT-TEST #1]: test_3_obj_construct -- ***END***!!\n");
}

void unit_test_4_obj_construct_resize();

void unit_test_4_obj_construct_resize()
{
       using namespace gms::common;
       constexpr std::size_t init_nx{4096ull};
       constexpr std::size_t resize_nx{1024ull};
       int32_t status;
       const char * fname{"UNIT_TEST_Output_Random_Data_test_4_obj_construct_resize.csv"};
       std::clock_t seed{0ULL};
       std::uniform_real_distribution<float> distro{};
       float * __restrict ptr_rand_buf  = NULL;
       float * __restrict ptr_rand_buf2 = NULL;
       FILE * fp{NULL};
       //darray_r4_t test_2nd_Ctor = darray_r4_t(nx);
       printf("[UNIT-TEST #1]: test_4_obj_construct_resize -- ***START***!!\n");
       ptr_rand_buf  = reinterpret_cast<float * __restrict>(gms_mm_malloc(init_nx*sizeof(float),64ull));
       ptr_rand_buf2 = reinterpret_cast<float * __restrict>(gms_mm_malloc(resize_nx*sizeof(float),64ull)); 
       gms_vec_sptr_t<float> ctor_2_arg;
       char * ctor_name{gms::common::demangle(typeid(ctor_2_arg).name(),status)};
       if(status==0 && ctor_name != NULL)
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
       }
       else
       {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(ctor_2_arg).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
       }
       fp = fopen(fname,"w+");
       if(!fp)
       {
           std::perror("fopen failed to open a file -- TERMINATING!!");
           std::exit(EXIT_FAILURE);
       }
       printf("[UNIT-TEST]: -- Start of random data generation and initialization.\n");
       seed = std::clock();
       distro  = std::uniform_real_distribution<float>(0.0f,1.0f);
       auto rdev1{std::mt19937(seed)};
       for(std::size_t __i{0ULL}; __i != init_nx; ++__i)
       {
           const float rc{distro.operator()(rdev1)};
           ptr_rand_buf[__i] = rc;
       }
       for(std::size_t __i{0ULL}; __i != resize_nx; ++__i)
       {
           const float rc{distro.operator()(rdev1)};
           ptr_rand_buf2[__i] = rc;
       }
       ctor_2_arg = gms_vec_sptr_t<float>(ptr_rand_buf,init_nx);
       for(std::size_t __i{0ULL}; __i != resize_nx; ++__i)
       {   const float rc = ptr_rand_buf2[__i];
           ctor_2_arg.push_back(rc);
       }
       std::cout << "size=" << ctor_2_arg.size() << ",capacity=" << ctor_2_arg.capacity() << "\n";
       printf("[UNIT-TEST]: End of random data generation and initialization.\n");
       printf("[UNIT-TEST]: Dumping: m_data to file=%s\n",fname);
       for(std::size_t __i{0ull}; __i != ctor_2_arg.size(); ++__i)
       {
           fprintf(fp,"[%.7f]\n", ctor_2_arg.operator[](__i));
                   
       }
       
#if 0
       std::cout << ctor_1_arg << "\n";
#endif     
       fclose(fp);
       if(ptr_rand_buf  && init_nx!=0ull)   { gms_mm_free(ptr_rand_buf);}
       if(ptr_rand_buf2 && resize_nx!=0ull) { gms_mm_free(ptr_rand_buf2);}
       printf("[UNIT-TEST #1]: test_4_obj_construct_resize -- ***END***!!\n");
}


int main()
{
   unit_test_1_obj_construct();
   unit_test_2_obj_init_construct();
   unit_test_3_obj_construct();
   unit_test_4_obj_construct_resize();
   return 0;
}
