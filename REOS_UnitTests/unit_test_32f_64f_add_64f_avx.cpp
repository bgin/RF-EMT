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
#include "GMS_32f_64f_add_64f_avx.h" 


static const char *funcs[] = {"dsp_32f_64f_add_64f_u_avx", "dsp_32f_64f_add_64f_a_avx"};

void unit_test_dsp_32f_64f_add_64f_u_avx()
{
     
     constexpr int32_t npts{1000};
     double buf_c[npts];
     float  buf_a[npts];
     double buf_b[npts];
     std::clock_t seed_a{0L};
     std::clock_t seed_b{0L};
     const char * __restrict fname = "UNIT_TEST_OutputData_dsp_32f_64f_add_64f_u_avx.csv";
     bool  fopen_failed = false;
     FILE * fp{fopen(fname,"w+")};
     if(!fp)
     {
         std::perror("[ERROR]: Failed to open a file!!");
         fopen_failed = true;
     }
     printf("Unit-Test]: -- Start of data initialization\n");
     seed_a = std::clock();
     auto srand_a{std::bind(std::uniform_real_distribution<float>{0.0f,1.0f},
                           std::mt19937(seed_a))};
     seed_b = std::clock();
     auto srand_b{std::bind(std::uniform_real_distribution<double>{0.0,1.0},
                           std::mt19937(seed_b))};
     for(int32_t i = 0; i != npts; ++i)
     {
         const float tra{srand_a()};
         buf_a[i] = tra;
         const double trb{srand_b()};
         buf_b[i] = trb;
     }
     printf("[Unit-Test]: -- End of data initialization\n");
     
     printf("[Unit-Test]: %s() -- START\n",funcs[0]);
     dsp_32f_64f_add_64f_u_avx(&buf_c[0],&buf_a[0],&buf_b[0],npts);
     if(fopen_failed)
     {
        printf("Dumping results to screen\n");
        for(int32_t i = 0; i != npts; ++i)
        {
            printf("%.16f\n",buf_c[i]);
        }
     }
     for(int32_t i = 0; i != npts; ++i)
     {
         fprintf(fp,"%.16f\n",buf_c[i]);
     }
     printf("[Unit-Test]: %s() -- END\n",funcs[0]);
     fclose(fp);

}

void unit_test_dsp_32f_64f_add_64f_a_avx()
{
     
     constexpr int32_t npts{1000};
     __ATTR_ALIGN__(32) double buf_c[npts];
     __ATTR_ALIGN__(32) float  buf_a[npts];
     __ATTR_ALIGN__(32) double buf_b[npts];
     std::clock_t seed_a{0L};
     std::clock_t seed_b{0L};
     const char * __restrict fname = "UNIT_TEST_OutputData_dsp_32f_64f_add_64f_a_avx.csv";
     bool  fopen_failed = false;
     FILE * fp{fopen(fname,"w+")};
     if(!fp)
     {
         std::perror("[ERROR]: Failed to open a file!!");
         fopen_failed = true;
     }
     printf("Unit-Test]: -- Start of data initialization\n");
     seed_a = std::clock();
     auto srand_a{std::bind(std::uniform_real_distribution<float>{0.0f,1.0f},
                           std::mt19937(seed_a))};
     seed_b = std::clock();
     auto srand_b{std::bind(std::uniform_real_distribution<double>{0.0,1.0},
                           std::mt19937(seed_b))};
     for(int32_t i = 0; i != npts; ++i)
     {
         const float tra{srand_a()};
         buf_a[i] = tra;
         const double trb{srand_b()};
         buf_b[i] = trb;
     }
     printf("[Unit-Test]: -- End of data initialization\n");
     
     printf("[Unit-Test]: %s() -- START\n",funcs[1]);
     dsp_32f_64f_add_64f_u_avx(&buf_c[0],&buf_a[0],&buf_b[0],npts);
     if(fopen_failed)
     {
        printf("Dumping results to screen\n");
        for(int32_t i = 0; i != npts; ++i)
        {
            printf("%.16f\n",buf_c[i]);
        }
     }
     for(int32_t i = 0; i != npts; ++i)
     {
         fprintf(fp,"%.16f\n",buf_c[i]);
     }
     printf("[Unit-Test]: %s() -- END\n",funcs[1]);
     fclose(fp);

}


int main()
{

    unit_test_dsp_32f_64f_add_64f_u_avx();
    unit_test_dsp_32f_64f_add_64f_a_avx();
    return 0;
}


