

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

#include "GMS_numerical_differentiation.h"


void gms::math::DLIMIT(const std::function<double(double)>& F, double D[], double DX[], double E[], 
            double X, double& H, int& M, double Tol) {
    const int Max = 100;
    double R[Max + 1];
    double Small = 1E-9;
    
    H = 1.0;
    DX[0] = H;
    D[0] = 0.5 * (F(X + H) - F(X - H)) / H;
    
    for (int N = 1; N <= 2; ++N) {
        H = H / 2.0;
        DX[N] = H;
        D[N] = 0.5 * (F(X + H) - F(X - H)) / H;
        E[N] = std::abs(D[N] - D[N - 1]);
        R[N] = 2.0 * E[N] / (std::abs(D[N]) + std::abs(D[N - 1]) + Small);
    }
    
    int N = 1;
    while ((E[N] > E[N + 1] || R[N] > Tol) && (N < Max)) {
        H = H / 2.0;
        DX[N + 2] = H;
        D[N + 2] = 0.5 * (F(X + H) - F(X - H)) / H;
        E[N + 2] = std::abs(D[N + 2] - D[N + 1]);
        R[N + 2] = 2.0 * E[N + 2] / (std::abs(D[N + 2]) + std::abs(D[N + 1]) + Small);
        N = N + 1;
    }
    M = N;
}

void gms::math::XLIMIT(const std::function<double(double)>& F, double D[], double DX[], double E[], 
            double X, double& H, int& M, double Tol) {
    const int Max = 100;
    double R[Max + 1];
    double Small = 1E-9;
    
    H = 1.0;
    DX[0] = H;
    D[0] = 0.5 * (F(X + H) - F(X - H)) / H;
    
    for (int N = 1; N <= 2; ++N) {
        H = H / 2.0;
        DX[N] = H;
        D[N] = 0.5 * (F(X + H) - F(X - H)) / H;
        E[N] = std::abs(D[N] - D[N - 1]);
        R[N] = 2.0 * E[N] / (std::abs(D[N]) + std::abs(D[N - 1]) + Small);
    }
    
    int N = 1;
    while ((E[N] > E[N + 1] || R[N] > Tol) && (N < Max)) {
        H = H / 2.0;
        DX[N + 2] = H;
        D[N + 2] = 0.5 * (F(X + H) - F(X - H)) / H;
        E[N + 2] = std::abs(D[N + 2] - D[N + 1]);
        R[N + 2] = 2.0 * E[N + 2] / (std::abs(D[N + 2]) + std::abs(D[N + 1]) + Small);
        N = N + 1;
    }
    M = N;
}

void gms::math::EXTRAP(const std::function<double(double)>& F, double D[][16], double DX[], 
            double X, double& H, double& Error, int& N, double Tol, double Delta) {
    double RelErr;
    double Small = 1E-7;
    H = 1.0;
    DX[0] = H;
    int J = 1;
    Error = 1.0;
    RelErr = 1.0;
    D[0][0] = 0.5 * (F(X + H) - F(X - H)) / H;
    
    while (RelErr > Tol && Error > Delta && J < 16) {
        H = H / 2.0;
        DX[J] = H;
        D[J][0] = 0.5 * (F(X + H) - F(X - H)) / H;
        for (int K = 1; K <= J; ++K) {
            D[J][K] = D[J][K - 1] + (D[J][K - 1] - D[J - 1][K - 1]) / (std::pow(4, K) - 1.0);
        }
        Error = std::abs(D[J][J] - D[J - 1][J - 1]);
        RelErr = 2.0 * Error / (std::abs(D[J][J]) + std::abs(D[J - 1][J - 1]) + Small);
        N = J;
        J = J + 1;
    }
}

void gms::math::XEXTRAP(const std::function<double(double)>& F, double D[][16], double DX[], 
             double X, double& H, double& Error, int& N, double Tol, double Delta) {
    double RelErr;
    double Small = 1E-7;
    H = 1.0;
    DX[0] = H;
    int J = 1;
    Error = 1.0;
    RelErr = 1.0;
    D[0][0] = 0.5 * (F(X + H) - F(X - H)) / H;
    
    while (RelErr > Tol && Error > Delta && J < 16) {
        H = H / 2.0;
        DX[J] = H;
        D[J][0] = 0.5 * (F(X + H) - F(X - H)) / H;
        for (int K = 1; K <= J; ++K) {
            D[J][K] = D[J][K - 1] + (D[J][K - 1] - D[J - 1][K - 1]) / (std::pow(4, K) - 1.0);
        }
        Error = std::abs(D[J][J] - D[J - 1][J - 1]);
        RelErr = 2.0 * Error / (std::abs(D[J][J]) + std::abs(D[J - 1][J - 1]) + Small);
        N = J;
        J = J + 1;
    }
}

void gms::math::DIFPOLY(double A[], const double X[], const double Y[], int N, double& Df) {
    for (int K = 0; K <= N; ++K) {
        A[K] = Y[K];
    }
    for (int J = 1; J <= N; ++J) {
        for (int K = N; K >= J; --K) {
            A[K] = (A[K] - A[K - 1]) / (X[K] - X[K - J]);
        }
    }
    double T = X[0];
    Df = A[1];
    double Prod = 1.0;
    for (int K = 2; K <= N; ++K) {
        Prod = Prod * (T - X[K - 1]);
        Df = Df + Prod * A[K];
    }
}

void gms::math::XDIFPOLY(double A[], const double X[], const double Y[], int N, double& Df) {
    for (int K = 0; K <= N; ++K) {
        A[K] = Y[K];
    }
    for (int J = 1; J <= N; ++J) {
        for (int K = N; K >= J; --K) {
            A[K] = (A[K] - A[K - 1]) / (X[K] - X[K - J]);
        }
    }
    double T = X[0];
    Df = A[1];
    double Prod = 1.0;
    for (int K = 2; K <= N; ++K) {
        Prod = Prod * (T - X[K - 1]);
        Df = Df + Prod * A[K];
    }
}
