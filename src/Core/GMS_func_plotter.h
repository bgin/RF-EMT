
/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef __GMS_FUNC_PLOTTER_H__
#define __GMS_FUNC_PLOTTER_H__ 301020260646 

#include <cstdint>
#include <cassert>
#include <fstream>
#include <iomanip>

namespace file_info {

     const unsigned int GMS_FUNC_PLOTTER_MAJOR = 2;
     const unsigned int GMS_FUNC_PLOTTER_MINOR = 0;
     const unsigned int GMS_FUNC_PLOTTER_MICRO = 0;
     const unsigned int GMS_FUNC_PLOTTER_FULLVER =
       1000U*GMS_FUNC_PLOTTER_MAJOR+100U*GMS_FUNC_PLOTTER_MINOR+
       10U*GMS_FUNC_PLOTTER_MICRO;
     const char * const GMS_FUNC_PLOTTER_CREATION_DATE = "30-10-2026 06:48 +00200 (FRI 30 OCT 2026 GMT+2)";
     const char * const GMS_FUNC_PLOTTER_BUILD_DATE    = __DATE__ " " __TIME__;
     static const char GMS_FUNC_PLOTTER_AUTHOR[]        = "Programmer: Bernard Gingold, beniekg@gmail.com"; 
     const char * const GMS_FUNC_PLOTTER_SYNOPSIS      = "Simplistic function (abscissae,ordinatae) plotter.";

}


namespace gms 
{

namespace common
{
                    inline 
                    static void 
                    create_function_xy_plot(const std::uint32_t n_samp,
                                         const float * __restrict func_arg,
                                         const float * __restrict func_val,
                                         const std::string &header,
                                         const std::string &title,
                                         const bool is_func_arg_present)
                    {
                         std::string plot_fname;
                         std::string sig_fname;
                         std::ofstream plot_unit;
                         std::ofstream sig_unit;
                         sig_fname = header+"_plot.txt";
                         sig_unit.open(sig_fname.c_str());
                         if(is_func_arg_present==true)
                         {
                              for(std::size_t __i{0ull}; __i != n_samp; ++__i)
                              {
                                     sig_unit << " " << func_arg[__i] << " "
                                                     << func_val[__i] << "\n";
                              }
                         }
                         else
                         {
                              for(std::size_t __i{0ull}; __i != n_samp; ++__i)
                              {
                                      sig_unit << " " << func_val[__i] << "\n";
         
                              }
                         }
                         sig_unit.close();
                         std::cout << "Created function data file \"" << sig_fname << "\".\n";
                         plot_fname = header+"plot_commands.txt";
                         plot_unit.open(plot_fname.c_str());
                         plot_unit << "#" << plot_fname << "\n";
                         plot_unit << "#\n";
                         plot_unit << "# Usage:\n";
                         plot_unit << "# gnuplot < " << plot_fname << "\n";
                         plot_unit << "#\n";
                         plot_unit << "set term png\n";
                         plot_unit << "set output \"" << header << ".png\"\n";
                         plot_unit << "set xlabel 't'\n";
                         plot_unit << "set ylabel 'y(t)'\n";
                         plot_unit << "set title '" << title << "'\n";
                         plot_unit << "set grid\n";
                         plot_unit << "set style data lines\n";
                         if(is_func_arg_present==true)
                         {
                                   plot_unit << "plot \"" << sig_fname << "\" using 1:2 lw 1 linecolor rgb \"red\"\n";
                         }
                         else
                         {
                                   plot_unit << "plot \"" << sig_fname << "\" lw 1 linecolor rgb \"red\"\n";
                         }
                         plot_unit << "quit\n";
                         plot_unit.close();
                         std::cout << " Created function data file \"" << plot_fname << "\"\n";
                    }

                    inline 
                    static void 
                    create_function_xy_plot(const std::uint32_t n_samp,
                                         const double * __restrict func_arg,
                                         const double * __restrict func_val,
                                         const std::string &header,
                                         const std::string &title,
                                         const bool is_func_arg_present)
                    {
                         std::string plot_fname;
                         std::string sig_fname;
                         std::ofstream plot_unit;
                         std::ofstream sig_unit;
                         sig_fname = header+"_plot.txt";
                         sig_unit.open(sig_fname.c_str());
                         if(is_func_arg_present==true)
                         {
                              for(std::size_t __i{0ull}; __i != n_samp; ++__i)
                              {
                                     sig_unit << " " << func_arg[__i] << " "
                                                     << func_val[__i] << "\n";
                              }
                         }
                         else
                         {
                              for(std::size_t __i{0ull}; __i != n_samp; ++__i)
                              {
                                      sig_unit << " " << func_val[__i] << "\n";
         
                              }
                         }
                         sig_unit.close();
                         std::cout << "Created function data file \"" << sig_fname << "\".\n";
                         plot_fname = header+"plot_commands.txt";
                         plot_unit.open(plot_fname.c_str());
                         plot_unit << "#" << plot_fname << "\n";
                         plot_unit << "#\n";
                         plot_unit << "# Usage:\n";
                         plot_unit << "# gnuplot < " << plot_fname << "\n";
                         plot_unit << "#\n";
                         plot_unit << "set term png\n";
                         plot_unit << "set output \"" << header << ".png\"\n";
                         plot_unit << "set xlabel 't'\n";
                         plot_unit << "set ylabel 'y(t)'\n";
                         plot_unit << "set title '" << title << "'\n";
                         plot_unit << "set grid\n";
                         plot_unit << "set style data lines\n";
                         if(is_func_arg_present==true)
                         {
                                   plot_unit << "plot \"" << sig_fname << "\" using 1:2 lw 1 linecolor rgb \"red\"\n";
                         }
                         else
                         {
                                   plot_unit << "plot \"" << sig_fname << "\" lw 1 linecolor rgb \"red\"\n";
                         }
                         plot_unit << "quit\n";
                         plot_unit.close();
                         std::cout << " Created function data file \"" << plot_fname << "\"\n";
                    }
                     

}

}

#endif /*__GMS_FUNC_PLOT_H__*/