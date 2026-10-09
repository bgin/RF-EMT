/*
 * Copyright (C) Bernard Gingold, 2020-2026 
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

#ifndef __GMS_AKIMA_QUADRATURE_C_H__
#define __GMS_AKIMA_QUADRATURE_C_H__ 081020260144

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define AKIMA_QUADRATURE_C_MAX_STAT_BUF 1
#if (AKIMA_QUADRATURE_C_MAX_STAT_BUF) == 1
#define MAX_BUF_SIZE 4093
#else 
#error "AKIMA_QUADRATURE_C_MAX_STAT_BUF -- Is Not Equal to 1!!"
#endif 

// Struct representing the tabular grid array

typedef struct alignas(32) 
{
    double * __restrict__ x;
    double * __restrict__ y;
    int size;
} tabular_data_t;

// Highly stable piecewise integration using Akima local slope boundaries
double integrate_tabular_akima(const tabular_data_t * __restrict__ data) {
    if (data == NULL || data->size < 5) return __builtin_nan("***[FATAL]*** -- data==NULL || size<5!!"); // Akima requires at least 5 points

    int n = data->size;
    double * __restrict__ x = data->x;
    double * __restrict__ y = data->y;

    // 1. Allocate arrays for consecutive finite differences (slopes)
    double * __restrict__ m = (double * __restrict__)malloc((n + 3) * sizeof(double));
    if(NULL==m && (n+3)>0) { return __builtin_nan("***[FATAL]*** -- NULL==m!!"); }
    
    double * __restrict__ t = (double * __restrict__)malloc(n * sizeof(double));
    if(NULL==t && n>0)     { return __builtin_nan("***[FATAL]*** -- NULL==t!!");}
    
    // Shift pointer so index matches mathematical convention m to m[n-2]
    double * __restrict__ m_shifted = m + 2;

    // Compute standard differences for internal steps
    for (int i = 0; i < n - 1; i++) {
        m_shifted[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]);
    }

    // Extrapolate endpoints boundary conditions stably
    m_shifted[-1] = 2.0 * m_shifted[0] - m_shifted[1];
    m_shifted[-2] = 2.0 * m_shifted[-1] - m_shifted[0];
    m_shifted[n-1] = 2.0 * m_shifted[n-2] - m_shifted[n-3];
    m_shifted[n]   = 2.0 * m_shifted[n-1] - m_shifted[n-2];

    // 2. Calculate stable Akima weights to prevent oscillation
    for (int i = 0; i < n; i++) {
        double ne0 = fabs(m_shifted[i+1] - m_shifted[i]);
        double ne1 = fabs(m_shifted[i-1] - m_shifted[i-2]);
        if ((ne0 + ne1) == 0.0) {
            t[i] = 0.5 * (m_shifted[i-1] + m_shifted[i]);
        } else {
            t[i] = (ne0 * m_shifted[i-1] + ne1 * m_shifted[i]) / (ne0 + ne1);
        }
    }

    // 3. Perform analytical piecewise cubic polynomial integration over each subinterval
    double total_integral = 0.0;
    for (int i = 0; i < n - 1; i++) {
        double h = x[i+1] - x[i];
        double p0 = y[i];
        double p1 = t[i];
        double p2 = (3.0 * m_shifted[i] - 2.0 * t[i] - t[i+1]) / h;
        double p3 = (t[i] + t[i+1] - 2.0 * m_shifted[i]) / (h * h);

        // Term-by-term analytical integration across interval 'h'
        /*
           double interval_area = p0 * h + 
                               p1 * (h * h) / 2.0 + 
                               p2 * (h * h * h) / 3.0 + 
                               p3 * (h * h * h * h) / 4.0;   
        */
        double interval_area = p0 * h + 
                               p1 * (h * h) * 0.5 + 
                               p2 * (h * h * h) * 0.33333333333333 + 
                               p3 * (h * h * h * h) * 0.25;
        total_integral += interval_area;
    }

    // Free local structural allocations
    free(m);
    free(t);

    return total_integral;
}


double 
integrate_tabular_akima_v2(const tabular_data_t * __restrict__ data)
{

    if (__builtin_expect(data==NULL,0)   || 
        __builtin_expect(data->size<5,0) ||
        __builtin_expect(data->size>MAX_BUF_SIZE,0)) {return __builtin_nan("***[FATAL]*** -- data==NULL || size<5 || size>MAX_BUF_SIZE"); } // Akima requires at least 5 points
    
    alignas(8) double t[MAX_BUF_SIZE+3];
    alignas(8) double m[MAX_BUF_SIZE];
    int n = data->size;
    double * __restrict__ x         = data->x;
    double * __restrict__ y         = data->y;
    double * __restrict__ m_shifted = &m[0]+2;
    double * __restrict__ p_t       = &t[0];
     // Compute standard differences for internal steps
    for (int i = 0; i < n - 1; i++) {
        m_shifted[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]);
    }

    // Extrapolate endpoints boundary conditions stably
    m_shifted[-1] = 2.0 * m_shifted[0] - m_shifted[1];
    m_shifted[-2] = 2.0 * m_shifted[-1] - m_shifted[0];
    m_shifted[n-1] = 2.0 * m_shifted[n-2] - m_shifted[n-3];
    m_shifted[n]   = 2.0 * m_shifted[n-1] - m_shifted[n-2];

    // 2. Calculate stable Akima weights to prevent oscillation
    for (int i = 0; i < n; i++) 
    {
        double ne0 = fabs(m_shifted[i+1] - m_shifted[i]);
        double ne1 = fabs(m_shifted[i-1] - m_shifted[i-2]);
        if ((ne0 + ne1) == 0.0) {
            p_t[i] = 0.5 * (m_shifted[i-1] + m_shifted[i]);
        } else {
            p_t[i] = (ne0 * m_shifted[i-1] + ne1 * m_shifted[i]) / (ne0 + ne1);
        }
    }

    // 3. Perform analytical piecewise cubic polynomial integration over each subinterval
    double total_integral = 0.0;
    for (int i = 0; i < n - 1; i++) 
    {
        double h = x[i+1] - x[i];
        double p0 = y[i];
        double p1 = t[i];
        double p2 = (3.0 * m_shifted[i] - 2.0 * p_t[i] - p_t[i+1]) / h;
        double p3 = (p_t[i] + p_t[i+1] - 2.0 * m_shifted[i]) / (h * h);

        // Term-by-term analytical integration across interval 'h'
        /*
           double interval_area = p0 * h + 
                               p1 * (h * h) / 2.0 + 
                               p2 * (h * h * h) / 3.0 + 
                               p3 * (h * h * h * h) / 4.0;   
        */
        double interval_area = p0 * h + 
                               p1 * (h * h) * 0.5 + 
                               p2 * (h * h * h) * 0.33333333333333 + 
                               p3 * (h * h * h * h) * 0.25;
        total_integral += interval_area;
    }
    return (total_integral);
} 

double integrate_tabular_akima_iface(const tabular_data_t * __restrict__ data) 
{
    const int n = data->size;
    if(MAX_BUF_SIZE<=n)
       return (integrate_tabular_akima_v2(data));
    else 
       return (integrate_tabular_akima(data));
}

/*
int akima_quadarture_example 
{
    // Example: Integrating y = x^2 from x = 0 to x = 2
    double x_vals[] = {0.0, 0.5, 1.0, 1.5, 2.0};
    double y_vals[] = {0.0, 0.25, 1.0, 2.25, 4.0};

    tabular_data_t table = {x_vals, y_vals, 5};
    double area = integrate_tabular_akima(&table);

    printf("Calculated Integral Area: %f\n", area); // Theoretical value = 2.666667
    return 0;
}
*/

#endif /*__GMS_AKIMA_QUADRATURE_C_H__*/
