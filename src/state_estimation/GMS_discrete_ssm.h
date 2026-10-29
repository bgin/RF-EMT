
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

#ifndef __GMS_DISCRETE_SSM_H__
#define __GMS_DISCRETE_SSM_H__ 281020261321

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
#include "GMS_eigen_defines.h"

namespace file_info 
{

     static const unsigned int GMS_DISCRETE_SSM_MAJOR = 1;
     static const unsigned int GMS_DISCRETE_SSM_MINOR = 1;
     static const unsigned int GMS_DISCRETE_SSM_MICRO = 0;
     static const unsigned int GMS_DISCRETE_SSM_FULLVER =
       1000U*GMS_DISCRETE_SSM_MAJOR+100U*GMS_DISCRETE_SSM_MINOR+
       10U*GMS_DISCRETE_SSM_MICRO;
     static const char GMS_DISCRETE_SSM_CREATION_DATE[]  = "28-10-2026 13:26PM +00200 (WED 28 OCT 2026 GMT+2)";
     static const char GMS_DISCRETE_SSM_BUILD_DATE[]    = __DATE__; 
     static const char GMS_DISCRETE_SSM_BUILD_TIME[]    = __TIME__;
     static const char GMS_DISCRETE_SSM_AUTHOR[]        = "Programmer: Bernard Gingold, beniekg@gmail.com";            
     static const char GMS_DISCRETE_SSM_SYNOPSIS[]      = "Discrete State-Space Model class implementation.";
}


namespace gms
{

namespace state_space
{

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
    static 
    inline MatrixXd zeros(Eigen::Index r, Eigen::Index c)
    {
        return MatrixXd::Zero(r, c);
    }

    __ATTR_ALWAYS_INLINE__
    static 
    inline VectorXd zeros_vec(Eigen::Index n)
    {
        return VectorXd::Zero(n);
    }

    __ATTR_ALWAYS_INLINE__
    static 
    inline MatrixXd eye(Eigen::Index n)
    {
        return MatrixXd::Identity(n, n);
    }

    __ATTR_ALWAYS_INLINE__
    static
    inline MatrixXd ones(Eigen::Index r, Eigen::Index c)
    {
        return MatrixXd::Ones(r, c);
    }

    __ATTR_ALWAYS_INLINE__
    static 
    inline VectorXd ones_vec(Eigen::Index n)
    {
        return VectorXd::Ones(n);
    }

    __ATTR_ALWAYS_INLINE__
    inline MatrixXd inv()
    {
        Eigen::FullPivLU<MatrixXd> lu(this->m_A);
        if (lu.rank() < this->m_A.rows()) {return MatrixXd();}
        return lu.inverse();
    }

    __ATTR_ALWAYS_INLINE__ 
    inline MatrixXd transpose()
    {
        return this->m_A.transpose();
    }

    __ATTR_ALWAYS_INLINE__ 
    inline void transpose(MatrixXd &A)
    {
        A.operator=(this->m_A.transpose());
    }

    __ATTR_ALWAYS_INLINE__ 
    inline MatrixXd  diag() 
    {
        return this->m_A.diagonal();
    }

    __ATTR_ALWAYS_INLINE__ 
    inline void diag(MatrixXd &A)  
    {
        A.operator=(this->m_A.diagonal());
    }

    __ATTR_ALWAYS_INLINE__ 
    inline  MatrixXd  kron()
    {
        return Eigen::kroneckerProduct(this->m_A,this->m_B);
    }

    __ATTR_ALWAYS_INLINE__ 
    inline void kron(MatrixXd &C) 
    {
        C.operator=(Eigen::kroneckerProduct(this->m_A,this->m_B));
    }

    __ATTR_ALWAYS_INLINE__
    inline MatrixXd sqrtm()
    {
       // For the PSD matrices used in the original M-files, this is appropriate.
       Eigen::SelfAdjointEigenSolver<MatrixXd> es(this->m_A);
       if (es.info() != Eigen::Success) {return (MatrixXd());}
       return es.operatorSqrt();
    }

    __ATTR_ALWAYS_INLINE__
    inline sqrtm(MatrixXd &A)
    {
       // For the PSD matrices used in the original M-files, this is appropriate.
       Eigen::SelfAdjointEigenSolver<MatrixXd> es(this->m_A);
       if (es.info() != Eigen::Success) {return (MatrixXd());}
       A.operator=(es.operatorSqrt());
    }

    __ATTR_ALWAYS_INLINE__
    static
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

    template<typename Func>
    __ATTR_ALWAYS_INLINE__
    static
    inline MatrixXd rand_parametric_distr(Func &dist,
                                          Eigen::Index rows, Eigen::Index cols)
    {
        MatrixXd M(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                M(i, j) = dist();
        return M;
    }

    __ATTR_ALWAYS_INLINE__
    static
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

    template<typename Func> 
    __ATTR_ALWAYS_INLINE__
    static
    inline void rand_parametric_distr(MatrixXd& M,
                                      Func &dist)
    {
       const Eigen::Index rows = M.rows();
       const Eigen::Index cols = M.cols();
       for (int i = 0; i < rows; ++i)
           for (int j = 0; j < cols; ++j)
               M(i, j) = dist();
    }

    __ATTR_ALWAYS_INLINE__
    static
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
    static
    inline VectorXd randn_vec(Eigen::Index n)
    {
        static std::mt19937 rng(std::random_device{}());
        std::normal_distribution<double> dist(0.0, 1.0);
        VectorXd v(n);
        for (int i = 0; i < n; ++i)
            v(i) = dist(rng);
        return v;
}

};

}

}

#endif /*__GMS_DISCRETE_SSM_H__*/