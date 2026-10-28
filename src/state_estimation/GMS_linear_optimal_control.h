

/*
The author's address is:
A. E. Bryson
761 Mayfield Ave.
Stanford, CA 94305

The author's email address is:
brysonae@leland.stanford.edu

The user is welcome to modify these files provided he or she
does not sell or publish them without acknowledging the book.
Applied Linear Optimal Control: Examples and Algorithms
*/

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

/*
   Translated to the modern C++17 by Bernard Gingold, beniekg@gmail.com
   from the original source code "OPTEST", distributed 
   as a part of the author's book: "Applied Linear Optimal Control: Examples and Algorithms".
*/



#ifndef __GMS_LINEAR_OPTIMAL_CONTROL_H__
#define __GMS_LINEAR_OPTIMAL_CONTROL_H__ 241020260745

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>
/*This compiles*/
#include<unsupported/Eigen/KroneckerProduct>
#include <cmath>
#include <random>
#include <stdexcept>
#include <tuple>
#include <vector>
#include <complex>
#include <immintrin.h> 
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_LINEAR_OPTIMAL_CONTROL_MAJOR = 1;
     static const unsigned int GMS_LINEAR_OPTIMAL_CONTROL_MINOR = 1;
     static const unsigned int GMS_LINEAR_OPTIMAL_CONTROL_MICRO = 0;
     static const unsigned int GMS_LINEAR_OPTIMAL_CONTROL_FULLVER =
       1000U*GMS_LINEAR_OPTIMAL_CONTROL_MAJOR+100U*GMS_LINEAR_OPTIMAL_CONTROL_MINOR+
       10U*GMS_LINEAR_OPTIMAL_CONTROL_MICRO;
     static const char GMS_LINEAR_OPTIMAL_CONTROL_CREATION_DATE[]  = "24-10-2026 07:45PM +00200 (SAT 24 OCT 2026 GMT+2)";
     static const char GMS_LINEAR_OPTIMAL_CONTROL_BUILD_DATE[]    = __DATE__; 
     static const char GMS_LINEAR_OPTIMAL_CONTROL_BUILD_TIME[]    = __TIME__;
     static const char GMS_LINEAR_OPTIMAL_CONTROL_AUTHOR[]        = "Programmer: Bernard Gingold, beniekg@gmail.com, \
                                                                     Adapted from the Matlab original implementation, being a part of the book: \
                                                                     Applied Linear Optimal Control: Examples and Algorithms, by Bryson A.E.";
     static const char GMS_LINEAR_OPTIMAL_CONTROL_SYNOPSIS[]      = "Linear optimal control essential routines.";
}

namespace gms
{

namespace state_space
{

using MatrixXd = Eigen::MatrixXd;
using VectorXd = Eigen::VectorXd;
using MatrixXi = Eigen::MatrixXi;
using VectorXi = Eigen::VectorXi;

// row-major type alias
using MatXdRM  = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using MatXiRM  = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using MatXcdRM = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>; 

// ---------------------------
// Basic MATLAB-like helpers
// ---------------------------


// Header for AVX-512 intrinsics



__ATTR_ALWAYS_INLINE__ 
inline 
void lightweight_gemm_avx512(const MatXdRM& A, const MatXdRM& B, MatXdRM& C) 
{
    const Eigen::Index M = A.rows();
    const Eigen::Index K = A.cols();
    const Eigen::Index N = B.cols();

    C.setZero(M, N);

    const double  * __restrict__ a_ptr = A.data();
    const double  * __restrict__ b_ptr = B.data();
    double        * __restrict__ c_ptr = C.data();

    for (Eigen::Index i = 0; i < M; ++i) {
        for (Eigen::Index k = 0; k < K; ++k) {
            double a_val = a_ptr[i * K + k]; 
            __m512d a_vec = _mm512_set1_pd(a_val);
            Eigen::Index j = 0;
            for (; j <= N - 8; j += 8) {
               
                Eigen::Index b_idx = k * N + j;
                Eigen::Index c_idx = i * N + j;
                __m512d b_vec = _mm512_loadu_pd(&b_ptr[b_idx]);
                __m512d c_vec = _mm512_loadu_pd(&c_ptr[c_idx]);
                c_vec = _mm512_fmadd_pd(a_vec, b_vec, c_vec);
                _mm512_storeu_pd(&c_ptr[c_idx], c_vec);
            }

            for (; j < N; ++j) {
                c_ptr[i * N + j] += a_val * b_ptr[k * N + j];
            }
        }
    }
}

__ATTR_ALWAYS_INLINE__
inline MatrixXd zeros(Eigen::Index r, Eigen::Index c)
{
    return MatrixXd::Zero(r, c);
}

__ATTR_ALWAYS_INLINE__
inline VectorXd zeros_vec(Eigen::Index n)
{
    return VectorXd::Zero(n);
}

__ATTR_ALWAYS_INLINE__
inline MatrixXd eye(Eigen::Index n)
{
    return MatrixXd::Identity(n, n);
}

__ATTR_ALWAYS_INLINE__
inline MatrixXd ones(Eigen::Index r, Eigen::Index c)
{
    return MatrixXd::Ones(r, c);
}

__ATTR_ALWAYS_INLINE__
inline VectorXd ones_vec(Eigen::Index n)
{
    return VectorXd::Ones(n);
}

__ATTR_ALWAYS_INLINE__
inline MatrixXd inv(const MatrixXd& A)
{
    if (A.rows() != A.cols())
        throw std::invalid_argument("inv(): matrix must be square");

    Eigen::FullPivLU<MatrixXd> lu(A);
    if (lu.rank() < A.rows())
        throw std::runtime_error("inv(): matrix is singular or rank deficient");

    return lu.inverse();
}

__ATTR_ALWAYS_INLINE__
inline MatrixXd inv(MatrixXd&& A)
{
    if (A.rows() != A.cols())
        throw std::invalid_argument("inv(): matrix must be square");

    Eigen::FullPivLU<MatrixXd> lu(A);
    if (lu.rank() < A.rows())
        throw std::runtime_error("inv(): matrix is singular or rank deficient");

    return lu.inverse();
}

__ATTR_ALWAYS_INLINE__
inline MatrixXd transpose(const MatrixXd& A)
{
    return A.transpose();
}

__ATTR_ALWAYS_INLINE__
inline MatrixXd transpose(MatrixXd&& A)
{
    return A.transpose();
}

__ATTR_ALWAYS_INLINE__
inline VectorXd diag(const MatrixXd& M)
{
    return M.diagonal();
}

__ATTR_ALWAYS_INLINE__
inline VectorXd diag(MatrixXd&& M)
{
    return M.diagonal();
}

__ATTR_ALWAYS_INLINE__
inline MatrixXd kron(const MatrixXd& A, const MatrixXd& B)
{
    return Eigen::kroneckerProduct(A, B);
}

__ATTR_ALWAYS_INLINE__
inline MatrixXd kron(MatrixXd&& A, MatrixXd&& B)
{
    return Eigen::kroneckerProduct(A, B);
}

__ATTR_ALWAYS_INLINE__
inline MatrixXd sqrtm(const MatrixXd& A)
{
    // For the PSD matrices used in the original M-files, this is appropriate.
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(A);
    if (es.info() != Eigen::Success)
        throw std::runtime_error("sqrtm(): eigen decomposition failed");
    return es.operatorSqrt();
}

__ATTR_ALWAYS_INLINE__
inline MatrixXd randn(Eigen::Index rows, Eigen::Index cols)
{
    static std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> dist(0.0, 1.0);
    MatrixXd M(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            M(i, j) = dist(rng);
    return M;
}

__ATTR_ALWAYS_INLINE__
inline void randn(MatrixXd& M)
{
    static std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> dist(0.0, 1.0);
    const Eigen::Index rows = M.rows();
    const Eigen::Index cols = M.cols();
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            M(i, j) = dist(rng);
}

__ATTR_ALWAYS_INLINE__
inline MatrixXd randn(Eigen::Index rows, Eigen::Index cols,
                      const double mu, const double sigma)
{
    static std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> dist(mu, sigma);
    MatrixXd M(rows, cols);
    for (Eigen::Index i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            M(i, j) = dist(rng);
    return M;
}

__ATTR_ALWAYS_INLINE__
inline VectorXd randn_vec(Eigen::Index n)
{
    static std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> dist(0.0, 1.0);
    VectorXd v(n);
    for (int i = 0; i < n; ++i)
        v(i) = dist(rng);
    return v;
}

__ATTR_ALWAYS_INLINE__
inline 
void randn_vec(VectorXd& v)
{
    static std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> dist(0.0, 1.0);
    const Eigen::Index n = v.rows();
    for (int i = 0; i < n; ++i)
        v(i) = dist(rng);
}

__ATTR_ALWAYS_INLINE__
inline VectorXd randn_vec(Eigen::Index n,
                          const double mu,
                          const double sigma)
{
    static std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> dist(mu, sigma);
    VectorXd v(n);
    for (Eigen::Index i = 0; i < n; ++i)
        v(i) = dist(rng);
    return v;
}

__ATTR_ALWAYS_INLINE__
inline MatrixXd diff(const MatrixXd& M)
{
    // MATLAB diff(A) along columns for row-major sampled data.
    MatrixXd out = M;
    if (M.cols() > 1) {
        out = M.rightCols(M.cols() - 1) - M.leftCols(M.cols() - 1);
    } else {
        out = MatrixXd::Zero(M.rows(), 1);
    }
    return out;
}

__ATTR_ALWAYS_INLINE__
inline MatrixXd diff(MatrixXd&& M) noexcept(true)
{
    // MATLAB diff(A) along columns for row-major sampled data.
    MatrixXd out = M;
    if (M.cols() > 1) {
        out = M.rightCols(M.cols() - 1) - M.leftCols(M.cols() - 1);
    } else {
        out = MatrixXd::Zero(M.rows(), 1);
    }
    return out;
}

// ---------------------------
// Proper state_space_holder_t wrapper
// ---------------------------

/*
• A (State matrix): n × n (square matrix)
• B (Input-to-state matrix): n × m
• C (State-to-output matrix): p × n
• D (Feedthrough / direct-transmission matrix): p × m

Summary Table


Matrix	Description	Size (Rows × Columns)
A	State	n × n (n states)
B	Input	n × m (m inputs)
C	Output	p × n (p outputs)
D	Feedthrough	p × m (p outputs, m inputs)


Associated Vectors

• x(t) or \(x_{k}\) (State vector): n × 1
• u(t) or \(u_{k}\) (Input vector): m × 1
• y(t) or \(y_{k}\) (Output vector): p × 1
*/

/* Matrix based limited state-space ADT*/
struct alignas(64) state_space_holder_t
{
    /* state matrix: size= n*n*/
    MatrixXd A; 
    /* Input-to-state matrix: size= n*m*/
    MatrixXd B; 
    /* State-to-output matrix: size=p*n*/
    MatrixXd C; 
    /* Feedthtough matrix: size=p*m*/
    MatrixXd D;

    state_space_holder_t() = default;



    state_space_holder_t(const MatrixXd& A_, 
                  const MatrixXd& B_, 
                  const MatrixXd& C_, 
                  const MatrixXd& D_) noexcept(false)
        : 
        A(A_), 
        B(B_), 
        C(C_), 
        D(D_)
    {

    }

    state_space_holder_t(MatrixXd &&A_, 
               MatrixXd &&B_, 
               MatrixXd &&C_, 
               MatrixXd &&D_) noexcept(true)
        :
        A(std::move(A_)),
        B(std::move(B_)),
        C(std::move(C_)),
        D(std::move(D_))
    {

    }

    ~state_space_holder_t()
    {

    }
        

    std::tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd> ssdata() const
    {
        return {A, B, C, D};
    }

};

struct alignas(64) ss_model_ptr_holder_t
{
    double * __restrict__ m_A_ptr{nullptr};
    double * __restrict__ m_B_ptr{nullptr};
    double * __restrict__ m_C_ptr{nullptr};
    double * __restrict__ m_D_ptr{nullptr};
    double * __restrict__ m_X_ptr{nullptr};
    double * __restrict__ m_U_ptr{nullptr};
    double * __restrict__ m_Y_ptr{nullptr};
};

struct alignas(64) state_space_model_t
{
    Eigen::Index m_n;
    Eigen::Index m_m;
    Eigen::Index m_p;
    /* state matrix: size= n*n*/
    MatXdRM m_A; 
    /* Input-to-state matrix: size= n*m*/
    MatXdRM m_B; 
    /* State-to-output matrix: size=p*n*/
    MatXdRM m_C; 
    /* Feedthtough matrix: size=p*m*/
    MatXdRM m_D;
    /*x(t) or \(x_{k}\) (State vector): n × 1*/
    VectorXd m_X;
    /*u(t) or \(u_{k}\) (Input vector): m × 1*/
    VectorXd m_U;
    /*y(t) or \(y_{k}\) (Output vector): p × 1*/
    VectorXd m_Y;

    __ATTR_ALWAYS_INLINE__ 
    inline 
    state_space_model_t()
    :
    m_n{},
    m_m{},
    m_p{},
    m_A(),
    m_B(),
    m_C(),
    m_D(),
    m_X(),
    m_U(),
    m_Y()
    {

    }

    template<std::int32_t initialization_kind>
    __ATTR_ALWAYS_INLINE__
    inline 
    state_space_model_t(const Eigen::Index n,
                        const Eigen::Index m,
                        const Eigen::Index p)         
    :
    m_n((n<(std::ptrdiff_t)4)?(std::ptrdiff_t)32:n),
    m_m((m<(std::ptrdiff_t)4)?(std::ptrdiff_t)32:m),
    m_p((m<(std::ptrdiff_t)4)?(std::ptrdiff_t)32:p)
    {
       if constexpr(initialization_kind==1)
       {
           this->m_A = MatXdRM::Zero(m_n,m_n);
           this->m_B = MatXdRM::Zero(m_n,m_m);
           this->m_C = MatXdRM::Zero(m_p,m_n);
           this->m_D = MatXdRM::Zero(m_p,m_m);
           this->m_X = VectorXd::Zero(m_n);
           this->m_U = VectorXd::Zero(m_m);
           this->m_Y = VectorXd::Zero(m_p);
       }
       else if constexpr(initialization_kind==2)
       {
           this->m_A = MatXdRM::Random(m_n,m_n);
           this->m_B = MatXdRM::Random(m_n,m_m);
           this->m_C = MatXdRM::Random(m_p,m_n);
           this->m_D = MatXdRM::Random(m_p,m_m);
           this->m_X = VectorXd::Random(m_n);
           this->m_U = VectorXd::Random(m_m);
           this->m_Y = VectorXd::Random(m_p);
       }
    }

    __ATTR_ALWAYS_INLINE__
    inline 
    state_space_model_t(const Eigen::Index n,
                        const Eigen::Index m,
                        const Eigen::Index p)
    :
    m_n((n<(std::ptrdiff_t)4)?(std::ptrdiff_t)32:n),
    m_m((m<(std::ptrdiff_t)4)?(std::ptrdiff_t)32:m),
    m_p((m<(std::ptrdiff_t)4)?(std::ptrdiff_t)32:p),
    m_A(m_n,m_n),
    m_B(m_n,m_m),
    m_C(m_p,m_n),
    m_D(m_p,m_m),
    m_X(m_n),
    m_U(m_m),
    m_Y(m_p)
    {

    }
    
    __ATTR_ALWAYS_INLINE__ 
    inline state_space_model_t(const Eigen::Index n,
                               const Eigen::Index m,
                               const Eigen::Index p,
                               const ss_model_ptr_holder_t & ptr_holder)
    :
    m_n((n<(std::ptrdiff_t)4)?(std::ptrdiff_t)32:n),
    m_m((m<(std::ptrdiff_t)4)?(std::ptrdiff_t)32:m),
    m_p((m<(std::ptrdiff_t)4)?(std::ptrdiff_t)32:p),
    m_A(Eigen::Map<MatXdRM>(ptr_holder.m_A_ptr,m_n,m_n)),
    m_B(Eigen::Map<MatXdRM>(ptr_holder.m_B_ptr,m_n,m_m)),
    m_C(Eigen::Map<MatXdRM>(ptr_holder.m_C_ptr,m_p,m_n)),
    m_D(Eigen::Map<MatXdRM>(ptr_holder.m_D_ptr,m_p,m_m)),
    m_X(Eigen::Map<VectorXd>(ptr_holder.m_X_ptr,m_n)),
    m_U(Eigen::Map<VectorXd>(ptr_holder.m_U_ptr,m_m)),
    m_Y(Eigen::Map<VectorXd>(ptr_holder.m_Y_ptr,m_p))
    {

    }

    __ATTR_ALWAYS_INLINE__ 
    inline state_space_model_t(const state_space_model_t &other) 
    :
    m_n(other.m_n),
    m_m(other.m_m),
    m_p(other.m_p),
    m_A(other.m_A),
    m_B(other.m_B),
    m_C(other.m_C),
    m_D(other.m_D),
    m_X(other.m_X),
    m_U(other.m_U),
    m_Y(other.m_Y) 
    {

    }

    __ATTR_ALWAYS_INLINE__ 
    inline state_space_model_t(state_space_model_t &&other) noexcept(true)
    :
    m_n(other.m_n),
    m_m(other.m_m),
    m_p(other.m_p),
    m_A(std::move(other.m_A)),
    m_B(std::move(other.m_B)),
    m_C(std::move(other.m_C)),
    m_D(std::move(other.m_D)),
    m_X(std::move(other.m_X)),
    m_U(std::move(other.m_U)),
    m_Y(std::move(other.m_Y)) 
    {

    }

    __ATTR_ALWAYS_INLINE__ 
    inline   
    ~state_space_model_t() noexcept(true)
    {

    }   

    state_space_model_t & operator=(const state_space_model_t &other) 
    {
        if(__builtin_expect(this==&other,0)) { return (*this);}

        this->m_n = other.m_n;
        this->m_m = other.m_m;
        this->m_p = other.m_p;
        this->m_A.operator=(other.m_A);
        this->m_B.operator=(other.m_B);
        this->m_C.operator=(other.m_C);
        this->m_D.operator=(other.m_D);
        this->m_X.operator=(other.m_X);
        this->m_U.operator=(other.m_U);
        this->m_Y.operator=(other.m_Y);

        return (*this);
    } 

    state_space_model_t & operator=(state_space_model_t &&other) 
    {
        if(__builtin_expect(this==&other,0)) { return (*this);}
        
        this->m_n = other.m_n;
        this->m_m = other.m_m;
        this->m_p = other.m_p;
        this->m_A.operator=(std::move(other.m_A));
        this->m_B.operator=(std::move(other.m_B));
        this->m_C.operator=(std::move(other.m_C));
        this->m_D.operator=(std::move(other.m_D));
        this->m_X.operator=(std::move(other.m_X));
        this->m_U.operator=(std::move(other.m_U));
        this->m_Y.operator=(std::move(other.m_Y));

        return (*this);
    } 

    __ATTR_ALWAYS_INLINE__ 
    inline static   
    MatXdRM matd_cpy_from_ptr(double * __restrict__ ptr, 
                              const Eigen::Index rows,
                              const Eigen::Index cols) 
    {
        MatXdRM mat_from_ptr = Eigen::Map<MatXdRM>(ptr,rows,cols);
        return (mat_from_ptr);
    }

    __ATTR_ALWAYS_INLINE__  
    inline static 
    VectorXd vecd_cpy_from_ptr(double * __restrict__ ptr,
                               const std::int32_t rows)
    {
        VectorXd vec_from_ptr = Eigen::Map<VectorXd>(ptr,rows);
        return (vec_from_ptr);
    }

/*
    
2. Lifetime and Dangling Pointers

The raw pointer returned by .data() is only valid as long as the parent Eigen::MatrixXd object remains alive and 
is not resized. If the matrix goes out of scope or you 
call functions like .resize(), the pointer becomes a dangling pointer, leading to undefined behavior.

3. Strict Memory Alignment

Eigen may align its internal memory buffers for SIMD vectorization. 
If you are passing this raw pointer to legacy C functions, external libraries (like LAPACK/BLAS), 
or computing strides yourself, you can query the memory stride using mat.outerStride() and mat.innerStride() t
o ensure your pointer offsets remain mathematically accurate
*/

    __ATTR_ALWAYS_INLINE__ 
    inline
    const double * __restrict__ get_m_A_raw_ptr() const
    {
        return (this->m_A.data());
    } 

    __ATTR_ALWAYS_INLINE__ 
    inline 
    const double * __restrict__ get_m_B_raw_ptr() const 
    {
        return (this->m_B.data());
    }

    __ATTR_ALWAYS_INLINE__ 
    inline 
    const double * __restrict__ get_m_C_raw_ptr() const 
    {
        return (this->m_C.data());
    }

    __ATTR_ALWAYS_INLINE__ 
    inline 
    const double * __restrict__ get_m_D_raw_ptr() const 
    {
        return (this->m_D.data());
    }

    __ATTR_ALWAYS_INLINE__ 
    inline 
    const double * __restrict__ get_m_X_raw_ptr() const 
    {
        return (this->m_X.data());
    }

    __ATTR_ALWAYS_INLINE__ 
    inline 
    const double * __restrict__ get_m_U_raw_ptr() const 
    {
        return (this->m_U.data());
    }

    __ATTR_ALWAYS_INLINE__ 
    inline 
    const double * __restrict__ get_m_Y_raw_ptr() const 
    {
        return (this->m_Y.data());
    }

    __ATTR_ALWAYS_INLINE__ 
    inline static
    void check_matd_alignment(bool results[3],const MatXdRM &mat)
    {
        const double * __restrict__ mat_ptr = mat.data();
        const std::uintptr_t mem_addr       = reinterpret_cast<std::uintptr_t>(mat_ptr);
        constexpr std::uintptr_t val64      = (std::uintptr_t)64;
        constexpr std::uintptr_t val32      = (std::uintptr_t)32;
        constexpr std::uintptr_t val16      = (std::uintptr_t)16;  
        if(mem_addr%val64 == 0)
           results[0] = true;
        else 
           results[0] = false;
        if(mem_addr%val32 == 0)
           results[1] = true;
        else 
           results[1] = false;
        if(mem_addr%16 == 0)
           results[2] = true;
        else 
           results[2] = false;
    }  

    __ATTR_ALWAYS_INLINE__ 
    inline static 
    bool check_matd_rowmajor_contiguity(const MatXdRM &mat)
    {
        const std::ptrdiff_t inner = mat.innerStride();
        const std::ptrdiff_t outer = mat.outerStride();
        bool is_contiguous = false;
        is_contiguous      = (inner == 1 && outer == mat.rows());
        return (is_contiguous);
    }

    __ATTR_ALWAYS_INLINE__ 
    inline static 
    bool check_matd_colmajor_contiguity(const MatXdRM &mat)
    {
        const std::ptrdiff_t inner = mat.innerStride();
        const std::ptrdiff_t outer = mat.outerStride();
        bool is_contiguous = false;
        is_contiguous      = (inner == 1 && outer == mat.cols());
        return (is_contiguous);
    }

};

inline std::tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd> ssdata(const state_space_holder_t& sys)
{
    return sys.ssdata();
}

// discretization.
inline state_space_holder_t c2d(const state_space_holder_t& sys, double Ts, bool zeroOrderHold = true)
{
    const int n = sys.A.rows();
    const int m = sys.B.cols();
    const int p = sys.C.rows();

    if (n != sys.A.cols())
        throw std::invalid_argument("c2d(): A must be square");

    // Exact discretization using augmented matrix exponential:
    // [A  B]
    // [0  0] * Ts
    // expm( [A B; 0 0]*Ts ) = [Phi Gamma; 0 I]
    MatrixXd M = MatrixXd::Zero(n + m, n + m);
    M.block(0, 0, n, n) = sys.A;
    M.block(0, n, n, m) = sys.B;
    M = M * Ts;

    MatrixXd E = M.exp(); // from unsupported Eigen::MatrixFunctions
    MatrixXd Phi = E.block(0, 0, n, n);
    MatrixXd Gamma = E.block(0, n, n, m);

    // For direct feedthrough, D is preserved for discrete time at sample times.
    // Exact ZOH gives y[k] = C * x[k] + D * u[k]
    MatrixXd Dd = sys.D;
    return state_space_holder_t(Phi, Gamma, sys.C, Dd);
}

// ---------------------------
// Discrete-time LQR
// ---------------------------
// Solves K = argmin sum x^T Q x + u^T R u + 2 x^T N u
// assuming discrete-time dynamics x_{k+1} = A x_k + B u_k
inline MatrixXd dlqr(const MatrixXd& A,
                     const MatrixXd& B,
                     const MatrixXd& Q,
                     const MatrixXd& R,
                     const MatrixXd& N = MatrixXd())
{
    const int n = A.rows();
    const int m = B.cols();

    MatrixXd Np = N;
    if (Np.rows() == 0 || Np.cols() == 0) {
        Np = MatrixXd::Zero(n, m);
    }

    // Solve DARE using fixed-point iteration:
    // S = Q + A' S A - (A' S B + N) (R + B' S B)^{-1} (B' S A + N')
    MatrixXd S = Q;
    const int max_iter = 200;
    const double tol = 1e-12;

    for (int iter = 0; iter < max_iter; ++iter) {
        MatrixXd M = R + B.transpose() * S * B;
        MatrixXd T = B.transpose() * S * A + Np.transpose();

        MatrixXd S_next = Q + A.transpose() * S * A
            - (A.transpose() * S * B + Np) * inv(M) * T;

        if ((S_next - S).norm() < tol)
            break;

        S = S_next;
    }

    MatrixXd M = R + B.transpose() * S * B;
    MatrixXd T = B.transpose() * S * A + Np.transpose();

    return inv(M) * T;
}

// MATLAB-style alias (same as dlqr)
inline MatrixXd lqr(const MatrixXd& A,
                    const MatrixXd& B,
                    const MatrixXd& Q,
                    const MatrixXd& R,
                    const MatrixXd& N = MatrixXd())
{
    return dlqr(A, B, Q, R, N);
}

// ---------------------------
// Discrete-time Kalman filter / LQE
// ---------------------------
// For x_{k+1} = A x_k + G w_k
//     y_k = C x_k + v_k
// with covariances W and V
inline MatrixXd dlqe(const MatrixXd& A,
                     const MatrixXd& G,
                     const MatrixXd& C,
                     const MatrixXd& W,
                     const MatrixXd& V)
{
    const int n = A.rows();
    const int p = C.rows();

    // Solve steady-state Riccati:
    // P = A P A' + G W G' - A P C' (C P C' + V)^{-1} C P A'
    MatrixXd P = G * W * G.transpose();
    const int max_iter = 200;
    const double tol = 1e-12;

    for (int iter = 0; iter < max_iter; ++iter) {
        MatrixXd S = C * P * C.transpose() + V;
        MatrixXd K = P * C.transpose() * inv(S);

        MatrixXd P_next = A * P * A.transpose()
            + G * W * G.transpose()
            - (A * P * C.transpose()) * inv(S) * (C * P * A.transpose());

        if ((P_next - P).norm() < tol)
            break;

        P = P_next;
    }

    MatrixXd S = C * P * C.transpose() + V;
    return A * P * C.transpose() * inv(S);
}

// MATLAB-style alias
inline MatrixXd lqe(const MatrixXd& A,
                    const MatrixXd& G,
                    const MatrixXd& C,
                    const MatrixXd& W,
                    const MatrixXd& V)
{
    return dlqe(A, G, C, W, V);
}

// ---------------------------
// Discrete-time state-space simulation
// ---------------------------
// x[k+1] = A x[k] + B u[k]
// y[k]   = C x[k] + D u[k]
// U is m x N, x0 is n x 1
inline std::pair<MatrixXd, MatrixXd> lsim(const state_space_holder_t& sys,
                                         const MatrixXd& U,
                                         const VectorXd& x0)
{
    const int n = sys.A.rows();
    const int p = sys.C.rows();
    const int m = sys.B.cols();
    const int N = U.cols();

    if (U.rows() != m)
        throw std::invalid_argument("lsim(): U row count must equal B.cols()");
    if (x0.rows() != n)
        throw std::invalid_argument("lsim(): x0 length must equal A.rows()");

    MatrixXd X = MatrixXd::Zero(n, N + 1);
    MatrixXd Y = MatrixXd::Zero(p, N + 1);

    X.col(0) = x0;
    Y.col(0) = sys.C * x0;

    for (int k = 0; k < N; ++k) {
        VectorXd uk = U.col(k);
        X.col(k + 1) = sys.A * X.col(k) + sys.B * uk;
        Y.col(k + 1) = sys.C * X.col(k + 1) + sys.D * uk;
    }

    return {X, Y};
}

// MATLAB "initial" wrapper, approximated via discrete simulation over a time grid
inline std::pair<MatrixXd, MatrixXd> initial(const state_space_holder_t& sys,
                                           const VectorXd& x0,
                                           double tf,
                                           int Ns = 100)
{
    const int n = sys.A.rows();
    const int p = sys.C.rows();
    const int m = sys.B.cols();

    if (sys.A.rows() != sys.A.cols())
        throw std::invalid_argument("initial(): A must be square");

    if (x0.rows() != n)
        throw std::invalid_argument("initial(): x0 length mismatch");

    // Use zero input for initial response.
    MatrixXd U = MatrixXd::Zero(m, Ns + 1);
    return lsim(sys, U, x0);
}

// ---------------------------
// MATLAB::ss() construction
// ---------------------------

inline state_space_holder_t ss(const MatrixXd& A,
                    const MatrixXd& B,
                    const MatrixXd& C,
                    const MatrixXd& D)
{
    return state_space_holder_t(A, B, C, D);
}

// ---------------------------
// Utility: matrix exponential helper
// ---------------------------

inline MatrixXd expm(const MatrixXd& A)
{
    return A.exp();
}

// ---------------------------
// MATLAB-like time-varying gain storage (3D arrays)
// ---------------------------
// In Eigen, the natural representation is std::vector<MatrixXd>.
// This is used for the original MATLAB arrays such as Kd(:,:,k), Sb(:,:,k), etc.
using Matrix3D = std::vector<MatrixXd>;

// ---------------------------
// 1) dbat_u.m translation
// ---------------------------

struct DbatResult
{
    MatrixXd xs;
    MatrixXd ws;
    MatrixXd vss;
    MatrixXd Rxs;
    MatrixXd Rw;
    MatrixXd Rv;
};

inline DbatResult dbat_u(const MatrixXd& Ph,
                        const MatrixXd& Gu,
                        const MatrixXd& Gw,
                        const MatrixXd& Wd,
                        const MatrixXd& C0,
                        const MatrixXd& V0,
                        const VectorXd& y0,
                        const MatrixXd& Cs,
                        const MatrixXd& Vd,
                        const MatrixXd& ys,
                        const MatrixXd& Cf,
                        const MatrixXd& Vf,
                        const VectorXd& yf,
                        const MatrixXd& ud)
{
    const int n0 = C0.rows();
    const int ns = Ph.rows();
    const int nm = ys.rows();
    const int nyf = Cf.rows();
    const int nd = Gw.cols();
    const int nc = Gu.cols();
    const int N1 = ys.cols();
    const int Ns = N1 - 1;

    // Build Cyx = [Cs*Ph; Cs*Ph^2; ...; Cs*Ph^Ns]
    MatrixXd Cyx = MatrixXd::Zero(nm * Ns, ns);
    MatrixXd Z = Ph;
    for (int k = 0; k < Ns; ++k) {
        MatrixXd block = Cs * Z;
        Cyx.block(k * nm, 0, nm, ns) = block;
        Z = Z * Ph;
    }

    // Build Cyw = [Cs*Ph^(Ns-1)*Gw, ..., Cs*Gw] but in the Matlab ordering used there.
    MatrixXd Cyw = MatrixXd::Zero(nm * Ns, nd * Ns);
    MatrixXd Z2 = Ph;
    std::vector<MatrixXd> Cpieces;
    Cpieces.push_back(Cs * Gw);
    for (int k = 1; k < Ns; ++k) {
        Z2 = Z2 * Ph;
        Cpieces.push_back(Cs * Z2 * Gw);
    }
    for (int k = 0; k < Ns; ++k) {
        // This follows the same row-block/column-block pattern used in MATLAB
        // by placing the appropriate Cpiece into the right-side blocks.
        int col0 = (Ns - 1 - k) * nd;
        int row0 = k * nm;
        Cyw.block(row0, col0, nm, nd) = Cpieces[k];
    }

    // Cfw = [Cf*Gw, Cf*Ph*Gw, ..., Cf*Ph^(Ns-1)*Gw]
    MatrixXd Cfw = MatrixXd::Zero(nyf, nd * Ns);
    Z = Cf;
    std::vector<MatrixXd> Fpieces;
    Fpieces.push_back(Z * Gw);
    for (int k = 1; k < Ns; ++k) {
        Z = Z * Ph;
        Fpieces.push_back(Z * Gw);
    }
    for (int k = 0; k < Ns; ++k) {
        Cfw.block(0, k * nd, nyf, nd) = Fpieces[k];
    }

    // Build Cfx = Cf * Ph^Ns
    MatrixXd PhPow = MatrixXd::Identity(ns, ns);
    for (int i = 0; i < Ns; ++i)
        PhPow = PhPow * Ph;
    MatrixXd Cfx = Cf * PhPow;

    // Cw = [C0 0; Cs 0; Cyx Cyw; Cfx Cfw; 0 I]
    MatrixXd Cw = MatrixXd::Zero(n0 + (Ns + 1) * nm + nyf + nd * Ns,
                                ns + nd * Ns);
    Cw.block(0, 0, n0, ns) = C0;
    Cw.block(n0, 0, nm, ns) = Cs;
    Cw.block(n0 + nm, 0, nm * Ns, ns) = Cyx;
    Cw.block(n0 + nm + nm * Ns, 0, nyf, ns) = Cfx;
    Cw.block(n0 + nm + nm * Ns + nyf, ns, nd * Ns, nd * Ns)
        = MatrixXd::Identity(nd * Ns, nd * Ns);

    Cw.block(0, ns, n0, nd * Ns) = MatrixXd::Zero(n0, nd * Ns);
    Cw.block(n0, ns, nm, nd * Ns) = MatrixXd::Zero(nm, nd * Ns);
    Cw.block(n0 + nm, ns, nm * Ns, nd * Ns) = Cyw;
    Cw.block(n0 + nm + nm * Ns, ns, nyf, nd * Ns) = Cfw;

    // ya = [y0; ys(:,1)...ys(:,Ns+1); yf; 0...0]
    VectorXd ya = VectorXd::Zero(n0 + (Ns + 1) * nm + nyf + nd * Ns);
    ya.segment(0, n0) = y0;
    for (int k = 0; k <= Ns; ++k)
        ya.segment(n0 + k * nm, nm) = ys.col(k);

    ya.segment(n0 + (Ns + 1) * nm, nyf) = yf;

    // Build Cu = [0; 0; Cyu; Cfu; 0]
    MatrixXd Cyu = MatrixXd::Zero(nm * Ns, nc * Ns);
    // Build Cnu = [Cs*Gu, Cs*Ph*Gu, ..., Cs*Ph^(Ns-1)*Gu]
    MatrixXd Cnu = MatrixXd::Zero(nm * Ns, nc * Ns);
    Z = Ph;
    std::vector<MatrixXd> Ppieces;
    Ppieces.push_back(Cs * Gu);
    for (int k = 1; k < Ns; ++k) {
        Z = Z * Ph;
        Ppieces.push_back(Cs * Z * Gu);
    }
    for (int k = 0; k < Ns; ++k) {
        // Use same block ordering as the MATLAB code
        int row0 = k * nm;
        int col0 = (Ns - 1 - k) * nc;
        Cnu.block(row0, col0, nm, nc) = Ppieces[k];
    }
    for (int k = 0; k < Ns; ++k) {
        int row0 = k * nm;
        int col0 = 0;
        int ncols = (k + 1) * nc;
        Cyu.block(row0, 0, nm, ncols) = Cnu.block(row0, (Ns - k - 1) * nc, nm, ncols);
    }

    // Cfu = [Cf*Gu, Cf*Ph*Gu, ..., Cf*Ph^(Ns-1)*Gu]
    MatrixXd Cfu = MatrixXd::Zero(nyf, nc * Ns);
    Z = Cf;
    std::vector<MatrixXd> FUpieces;
    FUpieces.push_back(Z * Gu);
    for (int k = 1; k < Ns; ++k) {
        Z = Z * Ph;
        FUpieces.push_back(Z * Gu);
    }
    for (int k = 0; k < Ns; ++k)
        Cfu.block(0, k * nc, nyf, nc) = FUpieces[k];

    MatrixXd Cu = MatrixXd::Zero(n0 + (Ns + 1) * nm + nyf + nd * Ns,
                                nc * Ns);
    Cu.block(n0 + nm, 0, nm * Ns, nc * Ns) = Cyu;
    Cu.block(n0 + nm + nm * Ns, 0, nyf, nc * Ns) = Cfu;

    // ua = [ud(:,1); ...; ud(:,Ns)]
    VectorXd ua = VectorXd::Zero(nc * Ns);
    for (int k = 0; k < Ns; ++k)
        ua.segment(k * nc, nc) = ud.col(k);

    // Vai = block-diagonal weight matrix
    MatrixXd Vai = MatrixXd::Zero(n0 + (Ns + 1) * nm + nyf + Ns * nd,
                                  n0 + (Ns + 1) * nm + nyf + Ns * nd);
    Vai.block(0, 0, n0, n0) = inv(V0);
    for (int k = 0; k <= Ns; ++k) {
        int base = n0 + k * nm;
        Vai.block(base, base, nm, nm) = inv(Vd);
    }
    Vai.block(n0 + (Ns + 1) * nm, n0 + (Ns + 1) * nm, nyf, nyf) = inv(Vf);
    for (int k = 0; k < Ns; ++k) {
        int base = n0 + (Ns + 1) * nm + nyf + k * nd;
        Vai.block(base, base, nd, nd) = inv(Wd);
    }

    MatrixXd Wa = inv(Cw.transpose() * Vai * Cw);
    VectorXd wa = Wa * Cw.transpose() * Vai * (ya - Cu * ua);

    MatrixXd xs = MatrixXd::Zero(ns, Ns + 1);
    MatrixXd ws = MatrixXd::Zero(nd, Ns);
    xs.col(0) = wa.segment(0, ns);

    for (int k = 0; k < Ns; ++k) {
        ws.col(k) = wa.segment(ns + k * nd, nd);
        xs.col(k + 1) = Ph * xs.col(k) + Gu * ud.col(k) + Gw * ws.col(k);
    }

    MatrixXd vss = ys - Cs * xs;

    // Covariance approx
    MatrixXd Pha = MatrixXd::Zero((Ns + 1) * ns, ns + nd * Ns);
    MatrixXd Sx = MatrixXd::Identity(ns, ns);
    Pha.block(0, 0, ns, ns) = Sx;

    MatrixXd Pow = MatrixXd::Identity(ns, ns);
    for (int k = 1; k <= Ns; ++k) {
        Pha.block(k * ns, 0, ns, ns) = Pow;
        for (int j = 0; j < k; ++j) {
            MatrixXd Gpow = MatrixXd::Identity(ns, ns);
            for (int t = 0; t < j; ++t)
                Gpow = Gpow * Ph;
            Pha.block(k * ns, ns + j * nd, ns, nd) = Gpow * Gw;
        }
        Pow = Pow * Ph;
    }

    MatrixXd Xcov = Pha * Wa * Pha.transpose();

    DbatResult out;
    out.xs = xs;
    out.ws = ws;
    out.vss = vss;

    out.Rxs = MatrixXd::Zero(ns, Ns + 1);
    out.Rw = MatrixXd::Zero(nd, Ns);
    out.Rv = MatrixXd::Zero(nm, Ns + 1);

    VectorXd Rx1 = Xcov.diagonal().cwiseSqrt();
    for (int k = 0; k <= Ns; ++k)
        out.Rxs.col(k) = Rx1.segment(k * ns, ns);

    VectorXd W1 = Wa.diagonal().cwiseSqrt();
    VectorXd Rw1 = W1.segment(ns, Ns * nd);
    for (int k = 0; k < Ns; ++k)
        out.Rw.col(k) = Rw1.segment(k * nd, nd);

    for (int k = 0; k <= Ns; ++k) {
        MatrixXd blk = Xcov.block(k * ns, k * ns, ns, ns);
        out.Rv.col(k) = (Cs * blk * Cs.transpose()).diagonal().cwiseSqrt();
    }

    return out;
}

// ---------------------------
// 2) dbat.m translation
// ---------------------------

// Same structure as dbat_u but without u input.
inline DbatResult dbat(const MatrixXd& Ph,
                      const MatrixXd& Gw,
                      const MatrixXd& Wd,
                      const MatrixXd& C0,
                      const MatrixXd& V0,
                      const VectorXd& y0,
                      const MatrixXd& Cs,
                      const MatrixXd& Vd,
                      const MatrixXd& ys,
                      const MatrixXd& Cf,
                      const MatrixXd& Vf,
                      const VectorXd& yf)
{
    const int n0 = C0.rows();
    const int ns = Ph.rows();
    const int nm = ys.rows();
    const int nd = Gw.cols();
    const int nyf = Cf.rows();
    const int N1 = ys.cols();
    const int Ns = N1 - 1;

    // Cyx and Cyw
    MatrixXd Cyx = MatrixXd::Zero(nm * Ns, ns);
    MatrixXd Z = Ph;
    for (int k = 0; k < Ns; ++k) {
        Cyx.block(k * nm, 0, nm, ns) = Cs * Z;
        Z = Z * Ph;
    }

    MatrixXd Cyw = MatrixXd::Zero(nm * Ns, nd * Ns);
    std::vector<MatrixXd> Cpieces;
    Z = Ph;
    Cpieces.push_back(Cs * Gw);
    for (int k = 1; k < Ns; ++k) {
        Z = Z * Ph;
        Cpieces.push_back(Cs * Z * Gw);
    }
    for (int k = 0; k < Ns; ++k) {
        int col0 = (Ns - 1 - k) * nd;
        int row0 = k * nm;
        Cyw.block(row0, col0, nm, nd) = Cpieces[k];
    }

    MatrixXd Cfw = MatrixXd::Zero(nyf, nd * Ns);
    Z = Cf;
    std::vector<MatrixXd> Fpieces;
    Fpieces.push_back(Z * Gw);
    for (int k = 1; k < Ns; ++k) {
        Z = Z * Ph;
        Fpieces.push_back(Z * Gw);
    }
    for (int k = 0; k < Ns; ++k)
        Cfw.block(0, k * nd, nyf, nd) = Fpieces[k];

    MatrixXd PhPow = MatrixXd::Identity(ns, ns);
    for (int i = 0; i < Ns; ++i)
        PhPow = PhPow * Ph;
    MatrixXd Cfx = Cf * PhPow;

    MatrixXd Cw = MatrixXd::Zero(n0 + (Ns + 1) * nm + nyf + nd * Ns,
                                ns + nd * Ns);
    Cw.block(0, 0, n0, ns) = C0;
    Cw.block(n0, 0, nm, ns) = Cs;
    Cw.block(n0 + nm, 0, nm * Ns, ns) = Cyx;
    Cw.block(n0 + nm + nm * Ns, 0, nyf, ns) = Cfx;
    Cw.block(n0 + nm + nm * Ns + nyf, ns, nd * Ns, nd * Ns)
        = MatrixXd::Identity(nd * Ns, nd * Ns);

    Cw.block(0, ns, n0, nd * Ns) = MatrixXd::Zero(n0, nd * Ns);
    Cw.block(n0, ns, nm, nd * Ns) = MatrixXd::Zero(nm, nd * Ns);
    Cw.block(n0 + nm, ns, nm * Ns, nd * Ns) = Cyw;
    Cw.block(n0 + nm + nm * Ns, ns, nyf, nd * Ns) = Cfw;

    VectorXd ya = VectorXd::Zero(n0 + (Ns + 1) * nm + nyf + nd * Ns);
    ya.segment(0, n0) = y0;
    for (int k = 0; k <= Ns; ++k)
        ya.segment(n0 + k * nm, nm) = ys.col(k);
    ya.segment(n0 + (Ns + 1) * nm, nyf) = yf;

    MatrixXd Vai = MatrixXd::Zero(n0 + (Ns + 1) * nm + nyf + Ns * nd,
                                  n0 + (Ns + 1) * nm + nyf + Ns * nd);
    Vai.block(0, 0, n0, n0) = inv(V0);
    for (int k = 0; k <= Ns; ++k) {
        int base = n0 + k * nm;
        Vai.block(base, base, nm, nm) = inv(Vd);
    }
    Vai.block(n0 + (Ns + 1) * nm, n0 + (Ns + 1) * nm, nyf, nyf) = inv(Vf);
    for (int k = 0; k < Ns; ++k) {
        int base = n0 + (Ns + 1) * nm + nyf + k * nd;
        Vai.block(base, base, nd, nd) = inv(Wd);
    }

    MatrixXd Wa = inv(Cw.transpose() * Vai * Cw);
    VectorXd wa = Wa * Cw.transpose() * Vai * ya;

    MatrixXd xs = MatrixXd::Zero(ns, Ns + 1);
    MatrixXd ws = MatrixXd::Zero(nd, Ns);
    xs.col(0) = wa.segment(0, ns);

    for (int k = 0; k < Ns; ++k) {
        ws.col(k) = wa.segment(ns + k * nd, nd);
        xs.col(k + 1) = Ph * xs.col(k) + Gw * ws.col(k);
    }

    MatrixXd vss = Cs * xs - ys;

    // Covariance approximation
    MatrixXd Pha = MatrixXd::Zero((Ns + 1) * ns, ns + nd * Ns);
    Pha.block(0, 0, ns, ns) = MatrixXd::Identity(ns, ns);
    MatrixXd Pow = MatrixXd::Identity(ns, ns);
    for (int k = 1; k <= Ns; ++k) {
        Pha.block(k * ns, 0, ns, ns) = Pow;
        for (int j = 0; j < k; ++j) {
            MatrixXd Gpow = MatrixXd::Identity(ns, ns);
            for (int t = 0; t < j; ++t)
                Gpow = Gpow * Ph;
            Pha.block(k * ns, ns + j * nd, ns, nd) = Gpow * Gw;
        }
        Pow = Pow * Ph;
    }

    MatrixXd Xcov = Pha * Wa * Pha.transpose();

    DbatResult out;
    out.xs = xs;
    out.ws = ws;
    out.vss = vss;

    out.Rxs = MatrixXd::Zero(ns, Ns + 1);
    out.Rw = MatrixXd::Zero(nd, Ns);
    out.Rv = MatrixXd::Zero(nm, Ns + 1);

    VectorXd Rx1 = Xcov.diagonal().cwiseSqrt();
    for (int k = 0; k <= Ns; ++k)
        out.Rxs.col(k) = Rx1.segment(k * ns, ns);

    VectorXd W1 = Wa.diagonal().cwiseSqrt();
    VectorXd Rw1 = W1.segment(ns, Ns * nd);
    for (int k = 0; k < Ns; ++k)
        out.Rw.col(k) = Rw1.segment(k * nd, nd);

    for (int k = 0; k <= Ns; ++k) {
        MatrixXd blk = Xcov.block(k * ns, k * ns, ns, ns);
        out.Rv.col(k) = (Cs * blk * Cs.transpose()).diagonal().cwiseSqrt();
    }

    return out;
}

// ---------------------------
// 3) tdfcf_cn.m translation
// ---------------------------

struct TdfcfCnResult
{
    MatrixXd x;
    MatrixXd ys;
    VectorXd yf;
    MatrixXd xb;
    MatrixXd Rxb;
    MatrixXd xh;
    MatrixXd Rxh;
    MatrixXd xs;
    MatrixXd Rxs;
};

inline TdfcfCnResult tdfcf_cn(const MatrixXd& Ph,
                             const MatrixXd& Gw,
                             const MatrixXd& Wd,
                             const MatrixXd& Cs,
                             const MatrixXd& Th,
                             const MatrixXd& V,
                             int Ns,
                             const VectorXd& x0,
                             const MatrixXd& X0,
                             const MatrixXd& Cf,
                             const MatrixXd& Vf,
                             int flg)
{
    const int ns = Ph.rows();
    const int nd = Gw.cols();
    const int nm = Cs.rows();
    const int nf = Cf.rows();

    MatrixXd Psi = V - Th * V * Th.transpose();

    MatrixXd w = sqrtm(Wd) * randn(nd, Ns);
    MatrixXd psi = sqrtm(Psi) * randn(nm, Ns + 1);

    MatrixXd x = MatrixXd::Zero(ns, Ns + 1);
    MatrixXd v = MatrixXd::Zero(nm, Ns + 1);
    MatrixXd ys = MatrixXd::Zero(nm, Ns);

    x.col(0) = x0 + sqrtm(X0) * randn_vec(ns);
    v.col(0) = sqrtm(V) * randn_vec(nm);

    for (int k = 0; k < Ns; ++k) {
        x.col(k + 1) = Ph * x.col(k) + Gw * w.col(k);
        v.col(k + 1) = Th * v.col(k) + psi.col(k);
        ys.col(k) = Cs * x.col(k) + v.col(k);
    }

    VectorXd yf = Cf * x.col(Ns) + sqrtm(Vf) * randn_vec(nf);

    TdfcfCnResult out;
    out.x = x;
    out.ys = ys;
    out.yf = yf;

    if (flg == 1) {
        const int naug = ns + nm;
        MatrixXd z1 = MatrixXd::Zero(ns, nm);
        MatrixXd z2 = MatrixXd::Zero(nd, nm);
        MatrixXd z3 = MatrixXd::Zero(ns, nm);

        MatrixXd Pha = MatrixXd::Zero(naug, naug);
        Pha.block(0, 0, ns, ns) = Ph;
        Pha.block(0, ns, ns, nm) = z1;
        Pha.block(ns, 0, nm, ns) = z1.transpose();
        Pha.block(ns, ns, nm, nm) = Th;

        MatrixXd Gwa = MatrixXd::Zero(naug, nd + nm);
        Gwa.block(0, 0, ns, nd) = Gw;
        Gwa.block(0, nd, ns, nm) = z3;
        Gwa.block(ns, 0, nm, nd) = z2.transpose();
        Gwa.block(ns, nd, nm, nm) = MatrixXd::Identity(nm, nm);

        MatrixXd Ca = MatrixXd::Zero(nm, naug);
        Ca.block(0, 0, nm, ns) = Cs;
        Ca.block(0, ns, nm, nm) = MatrixXd::Identity(nm, nm);

        MatrixXd Wa = MatrixXd::Zero(nd + nm, nd + nm);
        Wa.block(0, 0, nd, nd) = Wd;
        Wa.block(0, nd, nd, nm) = z2;
        Wa.block(nd, 0, nm, nd) = z2.transpose();
        Wa.block(nd, nd, nm, nm) = Psi;

        VectorXd xa0 = VectorXd::Zero(naug);
        xa0.segment(0, ns) = x0;
        xa0.segment(ns, nm).setZero();

        MatrixXd X0a = MatrixXd::Zero(naug, naug);
        X0a.block(0, 0, ns, ns) = X0;
        X0a.block(0, ns, ns, nm) = z1;
        X0a.block(ns, 0, nm, ns) = z1.transpose();
        X0a.block(ns, ns, nm, nm) = V;

        MatrixXd xb = MatrixXd::Zero(naug, Ns + 1);
        MatrixXd xh = MatrixXd::Zero(naug, Ns + 1);
        MatrixXd Rxb = MatrixXd::Zero(naug, Ns + 1);
        MatrixXd Rxh = MatrixXd::Zero(naug, Ns + 1);

        MatrixXd Pba = X0a;
        xb.col(0) = xa0;
        Rxb.col(0) = Pba.diagonal().cwiseSqrt();

        for (int k = 0; k < Ns; ++k) {
            MatrixXd S = V + Ca * Pba * Ca.transpose();
            MatrixXd L = Pba * Ca.transpose() * inv(S);

            xh.col(k) = xb.col(k) + L * (ys.col(k) - Ca * xb.col(k));
            MatrixXd Pht = Pba - L * S * L.transpose();
            Rxh.col(k) = Pht.diagonal().cwiseSqrt();

            xb.col(k + 1) = Pha * xh.col(k);
            Pba = Pha * Pht * Pha.transpose() + Gwa * Wa * Gwa.transpose();
            Rxb.col(k + 1) = Pba.diagonal().cwiseSqrt();
        }

        MatrixXd Sfinal = Vf + Cf * Pba * Cf.transpose();
        MatrixXd Lf = Pba * Cf.transpose() * inv(Sfinal);
        xh.col(Ns) = xb.col(Ns) + Lf * (yf - Cf * xb.col(Ns));
        MatrixXd Pht = Pba - Lf * Sfinal * Lf.transpose();
        Rxh.col(Ns) = Pht.diagonal().cwiseSqrt();

        out.xb = xb;
        out.Rxb = Rxb;
        out.xh = xh;
        out.Rxh = Rxh;
        out.xs = MatrixXd::Zero(ns, 1);
        out.Rxs = MatrixXd::Zero(ns, 1);
    }
    else if (flg == 2) {
        MatrixXd Vb = Psi + Cs * Gw * Wd * Gw.transpose() * Cs.transpose();
        MatrixXd Cb = Cs * Ph - Th * Cs;
        MatrixXd D1 = Gw * Wd * Gw.transpose() * Cs.transpose() * inv(Vb);

        MatrixXd yp = MatrixXd::Zero(nm, Ns - 1);
        for (int k = 0; k < Ns - 1; ++k) {
            yp.col(k) = ys.col(k + 1) - Th * ys.col(k);
        }

        MatrixXd xh = MatrixXd::Zero(ns, Ns);
        MatrixXd xs = MatrixXd::Zero(ns, Ns);
        MatrixXd Rxh = MatrixXd::Zero(ns, Ns);
        MatrixXd Rxs = MatrixXd::Zero(ns, Ns);

        MatrixXd L0 = X0 * Cs.transpose() * inv(V + Cs * X0 * Cs.transpose());
        xh.col(0) = x0 + L0 * (ys.col(0) - Cs * x0);
        MatrixXd Pt = X0 - L0 * Cs * X0;
        Rxh.col(0) = Pt.diagonal().cwiseSqrt();

        for (int k = 0; k < Ns - 1; ++k) {
            MatrixXd L = Pt * Cb.transpose() * inv(Vb + Cb * Pt * Cb.transpose());
            xs.col(k) = xh.col(k) + L * (yp.col(k) - Cb * xh.col(k));
            MatrixXd Psi_k = Pt - L * (Vb + Cb * Pt * Cb.transpose()) * L.transpose();
            Rxs.col(k) = Psi_k.diagonal().cwiseSqrt();

            xh.col(k + 1) = Ph * xs.col(k) + D1 * (yp.col(k) - Cb * xs.col(k));
            Pt = (Ph - D1 * Cb) * Psi_k * (Ph - D1 * Cb).transpose()
                 + Gw * Wd * Gw.transpose()
                 - D1 * Vb * D1.transpose();
            Rxh.col(k + 1) = Pt.diagonal().cwiseSqrt();
        }

        MatrixXd Lf = Pt * Cf.transpose() * inv(Vf + Cf * Pt * Cf.transpose());
        xs.col(Ns - 1) = xh.col(Ns - 1) + Lf * (yf - Cf * xh.col(Ns - 1));
        MatrixXd Psi_f = Pt - Lf * (Vf + Cf * Pt * Cf.transpose()) * Lf.transpose();
        Rxs.col(Ns - 1) = Psi_f.diagonal().cwiseSqrt();

        out.xb = MatrixXd::Zero(ns, 1);
        out.Rxb = MatrixXd::Zero(ns, 1);
        out.xh = xh;
        out.Rxh = Rxh;
        out.xs = xs;
        out.Rxs = Rxs;
    }
    else {
        throw std::invalid_argument("tdfcf_cn: flg must be 1 or 2");
    }

    return out;
}

// ---------------------------
// 4) lqg.m translation
// ---------------------------

struct LqgResult
{
    MatrixXd x;
    MatrixXd u;
    MatrixXd ys;
    MatrixXd xh;
    MatrixXd xm;
    MatrixXd um;
    VectorXd t;
    VectorXd Rx;
    VectorXd Re;
    VectorXd Ry;
    VectorXd Rey;
    VectorXd Ri;
    VectorXd Ru;
    MatrixXd K;
    MatrixXd L;
    VectorXd ev;
};

// ---------------------------
// A small helper for Lyapunov-like steady state solves
// ---------------------------

inline MatrixXd lyapunov_discrete(const MatrixXd& A, const MatrixXd& Q)
{
    // Solve X = A X A^T + Q
    MatrixXd X = Q;
    for (int it = 0; it < 100; ++it) {
        MatrixXd Xnext = A * X * A.transpose() + Q;
        if ((Xnext - X).norm() < 1e-12)
            return Xnext;
        X = Xnext;
    }
    return X;
}

inline MatrixXd GausssianRiccatiApprox(const MatrixXd& A,
                                      const MatrixXd& G,
                                      const MatrixXd& C,
                                      const MatrixXd& W,
                                      const MatrixXd& V)
{
    MatrixXd P = G * W * G.transpose();
    for (int it = 0; it < 100; ++it) {
        MatrixXd S = C * P * C.transpose() + V;
        MatrixXd Pnext = A * P * A.transpose()
            + G * W * G.transpose()
            - A * P * C.transpose() * inv(S) * C * P * A.transpose();
        if ((Pnext - P).norm() < 1e-12)
            return Pnext;
        P = Pnext;
    }
    return P;
}

inline LqgResult lqg(const MatrixXd& A,
                     const MatrixXd& Bu,
                     const MatrixXd& Cc,
                     const MatrixXd& Qc,
                     const MatrixXd& R,
                     const MatrixXd& Bw,
                     const MatrixXd& Cs,
                     const MatrixXd& W,
                     const MatrixXd& V,
                     const VectorXd& x0,
                     double tf,
                     int Ns)
{
    const int ns = A.rows();
    const int nc = Bu.cols();
    const int nm = Cs.rows();

    MatrixXd K = dlqr(A, Bu, Cc.transpose() * Qc * Cc, R);
    MatrixXd L = dlqe(A, Bw, Cs, W, V);

    // Continuous-time closed-loop covariance approximations
    // This mirrors the MATLAB logic in a discrete-time approximation.
    MatrixXd P = GausssianRiccatiApprox(A, Bw, Cs, W, V);
    MatrixXd Xh = lyapunov_discrete(A - Bu * K, L * V * L.transpose());
    MatrixXd X = Xh + P;

    LqgResult out;
    out.K = K;
    out.L = L;
    out.Re = P.diagonal().cwiseSqrt();
    out.Rx = X.diagonal().cwiseSqrt();
    out.Rey = (Cc * P * Cc.transpose()).diagonal().cwiseSqrt();
    out.Ry = (Cc * X * Cc.transpose()).diagonal().cwiseSqrt();
    out.Ri = (V + Cs * P * Cs.transpose()).diagonal().cwiseSqrt();
    out.Ru = (K * Xh * K.transpose()).diagonal().cwiseSqrt();

    out.t = VectorXd::LinSpaced(Ns + 1, 0.0, tf);
    out.x = MatrixXd::Zero(ns, Ns + 1);
    out.u = MatrixXd::Zero(nc, Ns + 1);
    out.ys = MatrixXd::Zero(nm, Ns + 1);
    out.xh = MatrixXd::Zero(ns, Ns + 1);
    out.xm = MatrixXd::Zero(ns, Ns + 1);
    out.um = MatrixXd::Zero(nc, Ns + 1);
    out.ev = VectorXd::Zero(2 * ns);

    return out;
}

} // namespace state_space
} // namespace gms

#endif /*__GMS_LINEAR_OPTIMAL_CONTROL_H__*/
