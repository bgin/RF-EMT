
/*
       NUMERICAL METHODS: FORTRAN Programs, (c) John H. Mathews 1994
!C     To accompany the text:
!C     NUMERICAL METHODS for Mathematics, Science and Engineering, 2nd Ed, 1992
!C     Prentice Hall, Englewood Cliffs, New Jersey, 07632, U.S.A.
!C     This free software is complements of the author.
*/
// Algorithm 6.1 (Differentiation Using Limits)
// Section 6.1, Approximating the Derivative, Page 326

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

#ifndef __GMS_NUMERICAL_METHODS_H__
#define __GMS_NUMERICAL_METHODS_H__ 161020261157

#include <cmath>
#include <functional>
#include "GMS_config.h"

namespace gms
{

namespace math
{

__ATTR_HOT__
__ATTR_ALIGN__(32)
void DLIMIT(const std::function<double(double)>& F, double D[], double DX[], double E[], 
            double X, double& H, int& M, double Tol);

__ATTR_HOT__
__ATTR_ALIGN__(32)
void XLIMIT(const std::function<double(double)>& F, double D[], double DX[], double E[], 
            double X, double& H, int& M, double Tol);

__ATTR_HOT__
__ATTR_ALIGN__(32)
void EXTRAP(const std::function<double(double)>& F, double D[][16], double DX[], 
            double X, double& H, double& Error, int& N, double Tol, double Delta);

__ATTR_HOT__
__ATTR_ALIGN__(32)
void XEXTRAP(const std::function<double(double)>& F, double D[][16], double DX[], 
             double X, double& H, double& Error, int& N, double Tol, double Delta);

__ATTR_HOT__
__ATTR_ALIGN__(32)
void DIFPOLY(double A[], const double X[], const double Y[], int N, double& Df);

__ATTR_HOT__
__ATTR_ALIGN__(32)
void XDIFPOLY(double A[], const double X[], const double Y[], int N, double& Df);

}

}

#endif // __GMS_NUMERICAL_METHODS_H__
