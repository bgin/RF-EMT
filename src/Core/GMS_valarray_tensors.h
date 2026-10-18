
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

#ifndef __GMS_VALARRAY_TENSORS_H__
#define __GMS_VALARRAY_TENSORS_H__ 181020260811

#include <valarray>
#include <iostream>
#include <cstdint>


namespace gms
{


  
  typedef std::valarray<float> Vector;


  class Matrix : public std::valarray<float> {

  private:
    std::size_t nrows_, ncols_, len_;

  public:
    Matrix();
    Matrix(std::size_t, std::size_t);
    std::size_t nrows();
    std::size_t ncols();
    void resize(std::size_t, std::size_t, float);

    // https://isocpp.org/wiki/faq/operator-overloading#matrix-subscript-op
    float& operator()(std::size_t, std::size_t);
    float  operator()(std::size_t, std::size_t) const;
    Matrix& operator=(const float);
    Matrix& operator=(const Matrix&);
  };

  void matrix_vector_product(Matrix&, Vector&, Vector&);
  float dot_product(Vector&, Vector&);
  void linear_solve(Matrix&, Vector&, std::valarray<int>&);
  
  std::ostream& operator<<(std::ostream&, Vector&);
  std::ostream& operator<<(std::ostream&, Matrix&);
  std::ostream& operator<<(std::ostream&, std::valarray<bool>&);
  Vector operator*(Matrix&, Vector&);

  
  inline Vector operator*(int j, Vector& v) {
    //return float(j) * v;
    return std::operator*(float(j), v);
  }


  inline Vector operator*(float j, Vector& v) {
    return std::operator*(j, v);
  }
  
  /*
  inline Vector operator*(float j, Vector& v) {
    return std::operator*(float(j), v);
  }
  
  inline Vector operator*(double j, Vector& v) {
    return std::operator*(float(j), v);
  }
  */

 
  inline std::size_t Matrix::nrows() {
    return nrows_;
  }

  inline std::size_t Matrix::ncols() {
    return ncols_;
  }

  inline float& Matrix::operator()(std::size_t m, std::size_t n) {
    return (*this)[m + nrows_ * n];
  }

  inline float Matrix::operator()(std::size_t m, std::size_t n) const {
    return (*this)[m + nrows_ * n];
  }

  inline Matrix& Matrix::operator=(const float v) {
    // Delegate to parent class:
    std::valarray<float>::operator=(v);
    return *this;
  }

  inline Matrix& Matrix::operator=(const Matrix &v) {
    // Delegate to parent class:
    std::valarray<float>::operator=(v);
    return *this;
  }

  //////////////////////////////////////////////////////////
  
  class Rank3Tensor : public std::valarray<float> {

  private:
    std::size_t d1_, d2_, d3_, len_;

  public:
    Rank3Tensor();
    Rank3Tensor(std::size_t, std::size_t, std::size_t);
    void resize(std::size_t, std::size_t, std::size_t, float);
    float& operator()(std::size_t, std::size_t, std::size_t);
    float  operator()(std::size_t, std::size_t, std::size_t) const;
    void set_row(Vector&, std::size_t, std::size_t);
  };

  inline float& Rank3Tensor::operator()(std::size_t j1, std::size_t j2, std::size_t j3) {
    return (*this)[j1 + d1_ * (j2 + d2_ * j3)];
  }

  inline float Rank3Tensor::operator()(std::size_t j1, std::size_t j2, std::size_t j3) const {
    return (*this)[j1 + d1_ * (j2 + d2_ * j3)];
  }

  //////////////////////////////////////////////////////////
  
  class Rank4Tensor : public std::valarray<float> {

  private:
    std::size_t d1_, d2_, d3_, d4_, len_;

  public:
    Rank4Tensor();
    Rank4Tensor(std::size_t, std::size_t, std::size_t, std::size_t);
    void resize(std::size_t, std::size_t, std::size_t, std::size_t, float);
    float& operator()(std::size_t, std::size_t, std::size_t, std::size_t);
    float  operator()(std::size_t, std::size_t, std::size_t, std::size_t) const;
    void set_row(Vector&, std::size_t, std::size_t, std::size_t);
  };

  inline float& Rank4Tensor::operator()(std::size_t j1, std::size_t j2, std::size_t j3, std::size_t j4) {
    return (*this)[j1 + d1_ * (j2 + d2_ * (j3 + d3_ * j4))];
  }

  inline float Rank4Tensor::operator()(std::size_t j1, std::size_t j2, std::size_t j3, std::size_t j4) const {
    return (*this)[j1 + d1_ * (j2 + d2_ * (j3 + d3_ * j4))];
  }


} // namespace gms


#endif /*__GMS_VALARRAY_TENSORS_H__*/
