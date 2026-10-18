
/*BSD 2-Clause License

Copyright (c) 2020, Matt Landreman
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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


#include <iostream>
#include <cassert>
#include "GMS_valarray_tensors.h"


extern "C" {
  // Matrix-vector multiply:
  void dgemv_(char* TRANS, int* M, int* N, double* ALPHA, double* A, int* LDA,
	      double* X, int* INCX, double* BETA, double* Y, int* INCY);

  void sgemv_(char* TRANS, int* M, int* N, float* ALPHA, float* A, int* LDA,
	      float* X, int* INCX, float* BETA, float* Y, int* INCY);

  // Matrix-matrix multiply:
  // void dgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* ALPHA,
  //	      double* A, int* LDA, double* B, int* LDB, double* BETA, double* C, int* LDC);

  // Solve linear system:
  void dgesv_(int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);

  void sgesv_(int* N, int* NRHS, float* A, int* LDA, int* IPIV, float* B, int* LDB, int* INFO);
}

#ifdef SINGLE

#define gemv_ sgemv_
#define gesv_ sgesv_

#else

#define gemv_ dgemv_
#define gesv_ dgesv_

#endif

using namespace gms;


Matrix::Matrix()
  : std::valarray<float>(1) // Call constructor of base class.
{
  nrows_ = 1;
  ncols_ = 1;
  len_ = 1;
}

Matrix::Matrix(std::size_t nrows_in, std::size_t ncols_in)
  : std::valarray<float>(nrows_in * ncols_in) // Call constructor of base class.
{
  nrows_ = nrows_in;
  ncols_ = ncols_in;
  len_ = nrows_ * ncols_;
}

void Matrix::resize(std::size_t nrows_in, std::size_t ncols_in, float v) {
  nrows_ = nrows_in;
  ncols_ = ncols_in;
  len_ = nrows_ * ncols_;
  std::valarray<float>::resize(nrows_ * ncols_, v);
}
	       
std::ostream& gms::operator<< (std::ostream& os, Vector& v) {
  for (std::size_t j = 0; j < v.size(); j++) {
    if (j > 0) os << " ";
    os << v[j];
  }
  return os;
}

std::ostream& gms::operator<< (std::ostream& os, Matrix& m) {
  for (std::size_t j = 0; j < m.nrows(); j++) {
    for (std::size_t k = 0; k < m.ncols(); k++) {
      if (k > 0) os << " ";
      os << m[j + m.nrows() * k];
    }
    os << std::endl;
  }
  return os;
}

std::ostream& gms::operator<< (std::ostream& os, std::valarray<bool>& v) {
  for (std::size_t j = 0; j < v.size(); j++) {
    if (j > 0) os << " ";
    os << v[j];
  }
  return os;
}

float gms::dot_product(Vector& u, Vector& v) {
  assert(u.size() == v.size());
  float sum = 0;
  for (int j = 0; j < u.size(); j++) {
    sum += u[j] * v[j];
  }
  return sum;
}

Vector gms::operator*(Matrix& m, Vector& v) {
  assert(m.ncols() == v.size());
  Vector result(m.nrows());

  /*
  // Non-BLAS simple-minded method:
  float sum = 0;
  for (int j = 0; j < m.nrows(); j++) {
    sum = 0;
    for (int k = 0; k < v.size(); k++) sum += m(j, k) * v[k];
    result[j] = sum;
  }
  */

  char TRANS = 'N';
  int INC = 1;
  int nrows = m.nrows();
  int ncols = m.ncols();
  float ALPHA = 1.0;
  float BETA = 0.0;
  if (single) {
    std::cout << "calling SINGLE precision BLAS" << std::endl;
  } else {
    std::cout << "calling DOUBLE precision BLAS" << std::endl;
  }
  gemv_(&TRANS, &nrows, &ncols, &ALPHA, &m(0, 0), &nrows,
	      &v[0], &INC, &BETA, &result[0], &INC);
  return result;
}


void gms::matrix_vector_product(Matrix& m, Vector& v, Vector& result) {
  assert(m.ncols() == v.size());
  assert(m.nrows() == result.size());
  /*
  // Non-BLAS simple-minded method:
  float sum = 0;
  for (int j = 0; j < m.nrows(); j++) {
    sum = 0;
    for (int k = 0; k < v.size(); k++) sum += m(j, k) * v[k];
    result[j] = sum;
  }
  */

  char TRANS = 'N';
  int INC = 1;
  int nrows = m.nrows();
  int ncols = m.ncols();
  float ALPHA = 1.0;
  float BETA = 0.0;
  gemv_(&TRANS, &nrows, &ncols, &ALPHA, &m(0, 0), &nrows,
	      &v[0], &INC, &BETA, &result[0], &INC);
}

/** Solve a linear system A x = b for x.
 *
 *  Like LAPACK's *gesv, this subroutine over-writes the matrix with
 *  the LU factors, and over-writes the vector with the solution.
 *
 *  Also like LAPACK, we require the integer work array IPIV to be
 *  provided. This lets the calling function have control over memory
 *  allocation.
 */
void gms::linear_solve(Matrix& m, Vector& v, std::valarray<int>& IPIV) {
  assert(m.ncols() == v.size());
  assert(m.ncols() == IPIV.size());
  assert(m.ncols() == m.nrows());

  int n = m.nrows();
  int INFO = 0;
  int one = 1;
  /*
    if (single) {
      std::cout << "calling SINGLE precision LAPACK" << std::endl;
    } else {
      std::cout << "calling DOUBLE precision LAPACK" << std::endl;
    }
  */
  gesv_(&n, &one, &m(0, 0), &n, &IPIV[0], &v[0], &n, &INFO);
  if (INFO != 0) {
    throw std::runtime_error("LAPACK error in *gesv");
  }
}

////////////////////////////////////////////////////


Rank3Tensor::Rank3Tensor()
  : std::valarray<float>(1) // Call constructor of base class.
{
  d1_ = 1;
  d2_ = 1;
  d3_ = 1;
  len_ = 1;
}

Rank3Tensor::Rank3Tensor(std::size_t d1_in, std::size_t d2_in, std::size_t d3_in)
  : std::valarray<float>(d1_in * d2_in * d3_in) // Call constructor of base class.
{
  d1_ = d1_in;
  d2_ = d2_in;
  d3_ = d3_in;
  len_ = d1_ * d2_ * d3_;
}

void Rank3Tensor::resize(std::size_t d1_in, std::size_t d2_in, std::size_t d3_in, float v) {
  d1_ = d1_in;
  d2_ = d2_in;
  d3_ = d3_in;
  len_ = d1_ * d2_ * d3_;
  std::valarray<float>::resize(len_, v);
}

void Rank3Tensor::set_row(Vector& v, std::size_t j2, std::size_t j3) {
  assert(v.size() == d1_);
  for (int j = 0; j < d1_; j++) {
    (*this)[j + d1_ * (j2 + d2_ * j3)] = v[j];
  }
}

////////////////////////////////////////////////////


Rank4Tensor::Rank4Tensor()
  : std::valarray<float>(1) // Call constructor of base class.
{
  d1_ = 1;
  d2_ = 1;
  d3_ = 1;
  d4_ = 1;
  len_ = 1;
}

Rank4Tensor::Rank4Tensor(std::size_t d1_in, std::size_t d2_in, std::size_t d3_in, std::size_t d4_in)
  : std::valarray<float>(d1_in * d2_in * d3_in * d4_in) // Call constructor of base class.
{
  d1_ = d1_in;
  d2_ = d2_in;
  d3_ = d3_in;
  d4_ = d4_in;
  len_ = d1_ * d2_ * d3_ * d4_;
}

void Rank4Tensor::resize(std::size_t d1_in, std::size_t d2_in, std::size_t d3_in, std::size_t d4_in, float v) {
  d1_ = d1_in;
  d2_ = d2_in;
  d3_ = d3_in;
  d4_ = d4_in;
  len_ = d1_ * d2_ * d3_ * d4_;
  std::valarray<float>::resize(len_, v);
}

void Rank4Tensor::set_row(Vector& v, std::size_t j2, std::size_t j3, std::size_t j4) {
  assert(v.size() == d1_);
  for (int j = 0; j < d1_; j++) {
    (*this)[j + d1_ * (j2 + d2_ * (j3 + d3_ * j4))] = v[j];
  }
}

